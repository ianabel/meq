#include "SurfaceFit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

/*
 * THE LINEAR ALGEBRA IS HAND-WRITTEN AND IT IS ABOUT A HUNDRED AND FIFTY LINES,
 * WHICH IS A DECISION AND NOT AN ACCIDENT.
 *
 * The header commits this file to including no MFEM, for the reason Zernike.hpp
 * and Profiles.hpp give: a fit against a point cloud is arithmetic on plain
 * doubles, so it should be testable and buildable where the finite element
 * library is not. That rules out mfem::DenseMatrix and its LAPACK-backed
 * decompositions, and Boost.Math -- the only other dependency meq_core has --
 * carries special functions and no dense linear algebra at all. So the choice
 * is between a hundred and fifty lines here and a new third-party dependency
 * for one least-squares solve.
 *
 * WHAT IS WRITTEN IS THE STABLE ROUTE AND NOT THE SHORT ONE. The short one is
 * the normal equations: form A^T A, factor it by Cholesky, done in thirty
 * lines. It is also wrong for this problem in a way that has no symptom.
 * Zernike is orthogonal on the disc under the measure discRadius d(discRadius),
 * so a sample set drawn that way gives a design matrix at condition number of
 * order one and normal equations would be perfectly safe -- but the sample set
 * is whatever surfaces the caller traced. Equispaced levels, or an annulus with
 * a hole in the middle where nobody traces, push the design matrix to 1e6 and
 * beyond; the normal matrix is then at 1e12 and Cholesky returns a coefficient
 * vector with no correct digits, through a route whose only outward sign is a
 * fit that looks merely poor. So:
 *
 *   * HOUSEHOLDER QR reduces the tall design matrix to a triangular factor
 *     without forming A^T A, so nothing is squared;
 *   * A ONE-SIDED JACOBI SWEEP on that triangular factor gives the singular
 *     values OF THE DESIGN MATRIX, to high relative accuracy even for the tiny
 *     ones -- which is exactly what an eigenvalue routine on A^T A cannot do,
 *     since a singular value of 1e-8 sits at 1e-16 there and is indistinguishable
 *     from round-off;
 *   * and the solve goes through that decomposition, so a direction the samples
 *     do not determine is DROPPED rather than inverted, and the count is
 *     reported.
 *
 * The condition number in SurfaceFitDiagnostics is therefore the design
 * matrix's own and not the normal matrix's. Those differ by a square. Quoting
 * the second where the first is meant halves the number of digits the reader
 * thinks are at risk, which is the wrong direction for a diagnostic to be
 * wrong in.
 *
 * COST. The QR is 2 m n^2 for m samples and n modes, the Jacobi sweep is n^3
 * per sweep on the small factor, and both right-hand sides -- r and z -- share
 * one factorisation. At the sizes IN-3 measures, a few thousand samples against
 * a couple of hundred modes, that is a tenth of a second. It is not on any
 * inner loop: a fit is built once and evaluated many times.
 */

namespace
{

	/// Exact for the small non-negative exponents a fit degree produces, and
	/// with no special case at base = 0 to think about. Lifted deliberately from
	/// Zernike.cpp, which needs the same thing for the same reason.
	double integerPower( double base, int exponent )
	{
		double result = 1.0;

		for ( int i = 0; i < exponent; ++i )
			result *= base;

		return result;
	}

	/**
	 * Householder QR of an m x n matrix stored row-major, m >= n, applying the
	 * same reflectors to two right-hand sides.
	 *
	 * On return the upper n x n triangle of @a a is R and the first n entries of
	 * each right-hand side are the corresponding entries of Q^T b. The reflectors
	 * themselves are not kept: this file never needs Q.
	 */
	void householderQr( std::vector<double> &a, std::size_t rows,
	                    std::size_t columns, std::vector<double> &firstRhs,
	                    std::vector<double> &secondRhs )
	{
		std::vector<double> reflector( rows, 0.0 );

		for ( std::size_t k = 0; k < columns; ++k )
		{
			// The reflector that zeroes a( k+1.., k ).
			double norm = 0.0;
			for ( std::size_t i = k; i < rows; ++i )
			{
				double const value = a[ i*columns + k ];
				reflector[ i ] = value;
				norm += value*value;
			}
			norm = std::sqrt( norm );

			if ( !( norm > 0.0 ) )
				continue;

			// Choose the sign that avoids cancellation in the leading entry.
			double const alpha = ( reflector[ k ] >= 0.0 ) ? -norm : norm;
			reflector[ k ] -= alpha;

			double reflectorNorm = 0.0;
			for ( std::size_t i = k; i < rows; ++i )
				reflectorNorm += reflector[ i ]*reflector[ i ];

			if ( !( reflectorNorm > 0.0 ) )
				continue;

			auto apply = [ & ]( double *column, std::size_t stride )
			{
				double dot = 0.0;
				for ( std::size_t i = k; i < rows; ++i )
					dot += reflector[ i ]*column[ i*stride ];

				double const factor = 2.0*dot/reflectorNorm;
				for ( std::size_t i = k; i < rows; ++i )
					column[ i*stride ] -= factor*reflector[ i ];
			};

			for ( std::size_t j = k; j < columns; ++j )
				apply( a.data() + j, columns );

			apply( firstRhs.data(), 1 );
			apply( secondRhs.data(), 1 );
		}
	}

	/**
	 * One-sided Jacobi on the columns of an n x n matrix held column-major in
	 * @a work, accumulating the right factor in @a rightFactor.
	 *
	 * On return the columns of @a work are U_j sigma_j and @a rightFactor is V,
	 * so the input was U Sigma V^T. The singular values come out UNSORTED, which
	 * costs nothing here: the caller wants their extremes and a floor, not an
	 * order.
	 */
	void oneSidedJacobi( std::vector<double> &work, std::vector<double> &rightFactor,
	                     std::size_t n, std::vector<double> &singularValues )
	{
		rightFactor.assign( n*n, 0.0 );
		for ( std::size_t j = 0; j < n; ++j )
			rightFactor[ j*n + j ] = 1.0;

		double const tolerance = 1.0e-15;
		int const maxSweeps = 60;

		for ( int sweep = 0; sweep < maxSweeps; ++sweep )
		{
			std::size_t rotations = 0;

			for ( std::size_t i = 0; i + 1 < n; ++i )
			{
				for ( std::size_t j = i + 1; j < n; ++j )
				{
					double alpha = 0.0;
					double beta = 0.0;
					double gamma = 0.0;

					for ( std::size_t k = 0; k < n; ++k )
					{
						double const left = work[ i*n + k ];
						double const right = work[ j*n + k ];
						alpha += left*left;
						beta += right*right;
						gamma += left*right;
					}

					if ( !( std::abs( gamma ) > tolerance*std::sqrt( alpha*beta ) ) )
						continue;

					++rotations;

					double const zeta = ( beta - alpha )/( 2.0*gamma );
					double const sign = ( zeta >= 0.0 ) ? 1.0 : -1.0;
					double const t = sign/( std::abs( zeta )
					                        + std::sqrt( 1.0 + zeta*zeta ) );
					double const cosine = 1.0/std::sqrt( 1.0 + t*t );
					double const sine = cosine*t;

					for ( std::size_t k = 0; k < n; ++k )
					{
						double const left = work[ i*n + k ];
						double const right = work[ j*n + k ];
						work[ i*n + k ] = cosine*left - sine*right;
						work[ j*n + k ] = sine*left + cosine*right;

						double const vLeft = rightFactor[ i*n + k ];
						double const vRight = rightFactor[ j*n + k ];
						rightFactor[ i*n + k ] = cosine*vLeft - sine*vRight;
						rightFactor[ j*n + k ] = sine*vLeft + cosine*vRight;
					}
				}
			}

			if ( rotations == 0 )
				break;
		}

		singularValues.assign( n, 0.0 );
		for ( std::size_t j = 0; j < n; ++j )
		{
			double sum = 0.0;
			for ( std::size_t k = 0; k < n; ++k )
				sum += work[ j*n + k ]*work[ j*n + k ];
			singularValues[ j ] = std::sqrt( sum );
		}
	}

}

namespace meq
{

	char const *fitBasisName( FitBasis basis )
	{
		switch ( basis )
		{
			case FitBasis::Zernike:       return "Zernike";
			case FitBasis::TensorProduct: return "tensor product";
		}

		return "unknown";
	}

	char const *fitRadialCoordinateName( FitRadialCoordinate coordinate )
	{
		switch ( coordinate )
		{
			case FitRadialCoordinate::DiscRadius:     return "disc radius";
			case FitRadialCoordinate::NormalisedFlux: return "normalised flux";
		}

		return "unknown";
	}

	namespace
	{

		/// The eigen-decomposition of a symmetric positive or negative definite
		/// 2 x 2 form, returned as the AxisShape the relabelling wants. The sign
		/// is irrelevant -- a maximum and a minimum have the same level-set
		/// ellipse -- but INDEFINITE is not, and is refused: an indefinite form
		/// is a SADDLE, its level sets are hyperbolae, and an ellipse fitted to
		/// them is meaningless rather than merely inaccurate.
		AxisShape shapeOfForm( double dRR, double dRZ, double dZZ, char const *where )
		{
			double const trace = dRR + dZZ;
			double const determinant = dRR*dZZ - dRZ*dRZ;

			if ( !( determinant > 0.0 ) || !std::isfinite( determinant )
			     || !std::isfinite( trace ) )
				throw std::invalid_argument(
					std::string( where ) + ": the quadratic form has determinant "
					+ std::to_string( determinant ) + ", so it is degenerate or"
					" indefinite. At a magnetic axis it is neither; an indefinite"
					" form is a SADDLE, whose level sets are hyperbolae and have no"
					" axis ellipse at all" );

			// The two SIGNED eigenvalues, and then the one of larger magnitude.
			// The sign itself is irrelevant here -- psi's axis is a minimum for
			// one sign of F and a maximum for the other, and the level-set
			// ellipse is the same either way, which CriticalPoints.hpp records at
			// length as a trap in the other direction. What matters is which
			// eigenvalue is the bigger, because a bigger second derivative climbs
			// out of the well faster and so belongs to the SHORT semi-axis.
			// gap = sqrt( half^2 - determinant ) is the textbook form and it is
			// WRONG HERE BY HALF THE DIGITS. On a nearly circular axis the two
			// eigenvalues nearly coincide, so half^2 and determinant nearly
			// cancel and the square root of their difference carries about
			// sqrt( eps ) relative error -- measured on an exactly circular
			// family, the recovered semi-axis ratio came out 1 + 1.3e-08 rather
			// than 1, and the relabelling below was then not quite the identity
			// where it must be exactly so. The equivalent form below has no
			// subtraction of nearly equal quantities in it at all and returns an
			// exact zero for an exact circle.
			double const half = 0.5*trace;
			double const gap = std::hypot( 0.5*( dRR - dZZ ), dRZ );
			double const upper = half + gap;
			double const lower = half - gap;

			double const dominant = ( std::abs( upper ) >= std::abs( lower ) )
				? upper : lower;
			double const other = ( std::abs( upper ) >= std::abs( lower ) )
				? lower : upper;

			// The eigenvector of `dominant`: ( lambda - dZZ, dRZ ), which is
			// exact for a symmetric form and degenerates only when dRZ vanishes,
			// where the axes are already the coordinate ones.
			double directionR = 1.0;
			double directionZ = 0.0;

			if ( std::abs( dRZ ) > 0.0 )
			{
				directionR = dominant - dZZ;
				directionZ = dRZ;
				double const norm = std::sqrt( directionR*directionR
				                               + directionZ*directionZ );
				if ( norm > 0.0 )
				{
					directionR /= norm;
					directionZ /= norm;
				}
			}
			else if ( std::abs( dZZ ) > std::abs( dRR ) )
			{
				directionR = 0.0;
				directionZ = 1.0;
			}

			AxisShape shape;
			shape.tilt = std::atan2( directionZ, directionR );
			shape.semiAxisRatio = std::sqrt( std::abs( other )/std::abs( dominant ) );

			return shape;
		}

	}

	AxisShape axisShapeFromHessian( double dRR, double dRZ, double dZZ )
	{
		return shapeOfForm( dRR, dRZ, dZZ, "axisShapeFromHessian" );
	}

	AxisShape axisShapeFromSamples( std::vector<SurfaceSample> const &samples,
	                                double axisR, double axisZ )
	{
		// The innermost surface, taken as every sample within a whisker of the
		// smallest normalised flux present. A whisker rather than an exact
		// comparison because a caller may have computed its levels by two routes.
		double smallest = std::numeric_limits<double>::infinity();
		for ( SurfaceSample const &sample : samples )
			smallest = std::min( smallest, sample.normalisedFlux );

		if ( !std::isfinite( smallest ) )
			throw std::invalid_argument(
				"axisShapeFromSamples: the sample set is empty" );

		// The form: p a^2 + 2 q a b + s b^2 = 1 through the innermost points, by
		// linear least squares in ( p, q, s ). Three unknowns and a normal matrix
		// of order three, which is the one place in this file where normal
		// equations are safe -- the design is 3 x 3 and the columns are the
		// second moments of a curve enclosing the origin, so it is nowhere near
		// singular unless the surface is a straight line.
		double matrix[ 3 ][ 3 ] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
		                            { 0.0, 0.0, 0.0 } };
		double rhs[ 3 ] = { 0.0, 0.0, 0.0 };
		std::size_t used = 0;

		for ( SurfaceSample const &sample : samples )
		{
			if ( sample.normalisedFlux > smallest*( 1.0 + 1.0e-9 ) + 1.0e-15 )
				continue;

			double const a = sample.r - axisR;
			double const b = sample.z - axisZ;
			double const row[ 3 ] = { a*a, 2.0*a*b, b*b };
			++used;

			for ( int i = 0; i < 3; ++i )
			{
				for ( int j = 0; j < 3; ++j )
					matrix[ i ][ j ] += row[ i ]*row[ j ];
				rhs[ i ] += row[ i ];
			}
		}

		if ( used < 3 )
			throw std::invalid_argument(
				"axisShapeFromSamples: the innermost surface carries "
				+ std::to_string( used ) + " points, which cannot determine a"
				" quadratic form" );

		// Gaussian elimination with partial pivoting, on three rows.
		for ( int k = 0; k < 3; ++k )
		{
			int pivot = k;
			for ( int i = k + 1; i < 3; ++i )
				if ( std::abs( matrix[ i ][ k ] ) > std::abs( matrix[ pivot ][ k ] ) )
					pivot = i;

			if ( pivot != k )
			{
				for ( int j = 0; j < 3; ++j )
					std::swap( matrix[ k ][ j ], matrix[ pivot ][ j ] );
				std::swap( rhs[ k ], rhs[ pivot ] );
			}

			if ( !( std::abs( matrix[ k ][ k ] ) > 0.0 ) )
				throw std::invalid_argument(
					"axisShapeFromSamples: the innermost surface does not determine"
					" a quadratic form -- its points are degenerate about the axis"
					" given" );

			for ( int i = k + 1; i < 3; ++i )
			{
				double const factor = matrix[ i ][ k ]/matrix[ k ][ k ];
				for ( int j = k; j < 3; ++j )
					matrix[ i ][ j ] -= factor*matrix[ k ][ j ];
				rhs[ i ] -= factor*rhs[ k ];
			}
		}

		double solution[ 3 ] = { 0.0, 0.0, 0.0 };
		for ( int k = 2; k >= 0; --k )
		{
			double sum = rhs[ k ];
			for ( int j = k + 1; j < 3; ++j )
				sum -= matrix[ k ][ j ]*solution[ j ];
			solution[ k ] = sum/matrix[ k ][ k ];
		}

		return shapeOfForm( solution[ 0 ], solution[ 1 ], solution[ 2 ],
		                    "axisShapeFromSamples" );
	}

	double shapedPoloidalAngle( AxisShape const &shape, double geometricAngle )
	{
		double const turned = geometricAngle - shape.tilt;
		double const along = std::cos( turned );
		double const across = std::sin( turned );

		// tan( label ) = ( short/long ) tan( angle ) in the ellipse's own frame,
		// which is the relation between an ellipse's parameter and the geometric
		// angle its point subtends. The tilt goes back on so that a circular
		// axis leaves the angle exactly alone.
		return std::atan2( shape.semiAxisRatio*across, along ) + shape.tilt;
	}

	double geometricPoloidalAngle( AxisShape const &shape, double label )
	{
		double const turned = label - shape.tilt;

		return std::atan2( std::sin( turned ),
		                   shape.semiAxisRatio*std::cos( turned ) ) + shape.tilt;
	}

	std::vector<SurfaceSample> relabelByAxisShape(
		std::vector<SurfaceSample> const &samples, AxisShape const &shape )
	{
		std::vector<SurfaceSample> out = samples;

		for ( SurfaceSample &sample : out )
			sample.theta = shapedPoloidalAngle( shape, sample.theta );

		return out;
	}

	std::size_t tensorProductModeCount( int maxDegree )
	{
		if ( maxDegree < 0 || maxDegree > maxZernikeDegree )
			throw std::invalid_argument(
				"tensorProductModeCount: the maximum degree must lie between 0 and "
				+ std::to_string( maxZernikeDegree ) + ", got "
				+ std::to_string( maxDegree ) );

		std::size_t const count = static_cast<std::size_t>( maxDegree ) + 1;
		return count*count;
	}

	SurfaceFit::SurfaceFit( int maxDegree,
	                        std::vector<SurfaceSample> const &samples,
	                        SurfaceFitOptions const &options )
		: degree( maxDegree ), option( options )
	{
		if ( maxDegree < 0 || maxDegree > maxZernikeDegree )
			throw std::invalid_argument(
				"SurfaceFit: the maximum degree must lie between 0 and "
				+ std::to_string( maxZernikeDegree ) + ", got "
				+ std::to_string( maxDegree ) );

		if ( !( options.discEdge > 0.0 ) || !std::isfinite( options.discEdge ) )
			throw std::invalid_argument(
				"SurfaceFit: discEdge is the normalised flux that maps to the edge"
				" of the disc and must be positive and finite" );

		// The modes. Both orders are by degree ascending and then by m ascending,
		// so both are PREFIX STABLE and a truncation is a truncation of the
		// coefficient vector -- which is what a study in mode number needs and is
		// the property zernikeModes() is documented to have.
		if ( options.basis == FitBasis::Zernike )
		{
			std::vector<ZernikeMode> const zernikeList = zernikeModes( maxDegree );
			modeList.reserve( zernikeList.size() );
			for ( ZernikeMode mode : zernikeList )
				modeList.push_back( FitMode{ mode.l, mode.m } );
		}
		else
		{
			modeList.reserve( tensorProductModeCount( maxDegree ) );
			for ( int d = 0; d <= maxDegree; ++d )
				for ( int m = -d; m <= d; ++m )
					modeList.push_back( FitMode{ d - ( m < 0 ? -m : m ), m } );
		}

		std::size_t const modes = modeList.size();
		std::size_t const rows = samples.size();

		if ( rows < modes )
			throw std::invalid_argument(
				"SurfaceFit: " + std::to_string( rows ) + " samples cannot determine "
				+ std::to_string( modes ) + " modes. Trace more surfaces, take more"
				" angles, or drop the degree" );

		diagnostic.samples = rows;
		diagnostic.modes = modes;
		diagnostic.smallestSampledRadius = std::numeric_limits<double>::infinity();
		diagnostic.largestSampledRadius = 0.0;

		// The design matrix, row-major, with the least-squares weights folded in
		// as their square roots so that the QR below minimises the weighted
		// residual and nothing downstream has to remember they are there.
		std::vector<double> matrix( rows*modes, 0.0 );
		std::vector<double> rhsR( rows, 0.0 );
		std::vector<double> rhsZ( rows, 0.0 );

		for ( std::size_t i = 0; i < rows; ++i )
		{
			SurfaceSample const &sample = samples[ i ];

			if ( !std::isfinite( sample.normalisedFlux ) || sample.normalisedFlux < 0.0 )
				throw std::invalid_argument(
					"SurfaceFit: a sample carries a negative or non-finite normalised"
					" flux. Psi_N is a fraction of the flux between the axis and the"
					" boundary and cannot be negative" );

			// A tenth of an ulp of slack, so that a caller that computed its
			// outermost level and its discEdge by two routes is not refused for
			// the last bit.
			if ( sample.normalisedFlux > options.discEdge*( 1.0 + 1.0e-12 ) )
				throw std::invalid_argument(
					"SurfaceFit: a sample at Psi_N = "
					+ std::to_string( sample.normalisedFlux )
					+ " lies outside the disc, whose edge is at Psi_N = "
					+ std::to_string( options.discEdge )
					+ ". Raise discEdge, or drop the surface" );

			if ( !( sample.weight > 0.0 ) || !std::isfinite( sample.weight ) )
				throw std::invalid_argument(
					"SurfaceFit: a sample carries a weight that is not positive and"
					" finite" );

			double const argument = argumentOf( sample.normalisedFlux );
			double const radius = std::sqrt( sample.normalisedFlux/options.discEdge );
			diagnostic.smallestSampledRadius =
				std::min( diagnostic.smallestSampledRadius, radius );
			diagnostic.largestSampledRadius =
				std::max( diagnostic.largestSampledRadius, radius );

			double const scale = std::sqrt( sample.weight );

			// modeValue() rather than evaluateMode(): the design matrix wants
			// values alone, and the derivatives cost three further special
			// function evaluations per entry on the Zernike branch.
			for ( std::size_t j = 0; j < modes; ++j )
				matrix[ i*modes + j ] = scale*modeValue( modeList[ j ], argument,
				                                         sample.theta );

			rhsR[ i ] = scale*sample.r;
			rhsZ[ i ] = scale*sample.z;
		}

		householderQr( matrix, rows, modes, rhsR, rhsZ );

		// The triangular factor, column-major, for the Jacobi sweep.
		std::vector<double> work( modes*modes, 0.0 );
		for ( std::size_t j = 0; j < modes; ++j )
			for ( std::size_t i = 0; i <= j; ++i )
				work[ j*modes + i ] = matrix[ i*modes + j ];

		std::vector<double> rightFactor;
		std::vector<double> singularValues;
		oneSidedJacobi( work, rightFactor, modes, singularValues );

		double largest = 0.0;
		double smallest = std::numeric_limits<double>::infinity();
		for ( double value : singularValues )
		{
			largest = std::max( largest, value );
			smallest = std::min( smallest, value );
		}

		if ( !( largest > 0.0 ) )
			throw std::runtime_error(
				"SurfaceFit: every singular value of the design matrix is zero, so"
				" the samples determine nothing at all" );

		diagnostic.largestSingularValue = largest;
		diagnostic.smallestSingularValue = smallest;
		diagnostic.conditionNumber = ( smallest > 0.0 )
			? largest/smallest : std::numeric_limits<double>::infinity();

		// c = V Sigma^-1 U^T y, with y the first `modes` entries of Q^T b and
		// U_j = work_j / sigma_j -- so U_j^T y / sigma_j is work_j^T y / sigma_j^2
		// and the columns never have to be normalised.
		double const floor = options.singularValueFloor*largest;
		coefficientR.assign( modes, 0.0 );
		coefficientZ.assign( modes, 0.0 );

		for ( std::size_t j = 0; j < modes; ++j )
		{
			if ( !( singularValues[ j ] > floor ) )
			{
				++diagnostic.discardedModes;
				continue;
			}

			double projectionR = 0.0;
			double projectionZ = 0.0;
			for ( std::size_t k = 0; k < modes; ++k )
			{
				projectionR += work[ j*modes + k ]*rhsR[ k ];
				projectionZ += work[ j*modes + k ]*rhsZ[ k ];
			}

			double const inverse = 1.0/( singularValues[ j ]*singularValues[ j ] );
			projectionR *= inverse;
			projectionZ *= inverse;

			for ( std::size_t k = 0; k < modes; ++k )
			{
				coefficientR[ k ] += rightFactor[ j*modes + k ]*projectionR;
				coefficientZ[ k ] += rightFactor[ j*modes + k ]*projectionZ;
			}
		}

		// The residual against the samples themselves, in the units of r and z
		// and WITHOUT the weights, because a residual is read as a distance.
		double sumR = 0.0;
		double sumZ = 0.0;
		for ( SurfaceSample const &sample : samples )
		{
			double r = 0.0;
			double z = 0.0;
			position( sample.normalisedFlux, sample.theta, r, z );

			double const errorR = r - sample.r;
			double const errorZ = z - sample.z;
			sumR += errorR*errorR;
			sumZ += errorZ*errorZ;
			diagnostic.worstR = std::max( diagnostic.worstR, std::abs( errorR ) );
			diagnostic.worstZ = std::max( diagnostic.worstZ, std::abs( errorZ ) );
		}

		diagnostic.residualR = std::sqrt( sumR/static_cast<double>( rows ) );
		diagnostic.residualZ = std::sqrt( sumZ/static_cast<double>( rows ) );
	}

	double SurfaceFit::argumentOf( double normalisedFlux ) const
	{
		double const fraction = normalisedFlux/option.discEdge;

		return ( option.coordinate == FitRadialCoordinate::DiscRadius )
			? std::sqrt( fraction ) : fraction;
	}

	double SurfaceFit::discRadiusOf( double normalisedFlux ) const
	{
		return std::sqrt( normalisedFlux/option.discEdge );
	}

	double SurfaceFit::normalisedFluxOf( double discRadius ) const
	{
		return discRadius*discRadius*option.discEdge;
	}

	double SurfaceFit::modeValue( FitMode mode, double argument,
	                              double theta ) const
	{
		if ( option.basis == FitBasis::Zernike )
			return zernike( mode.radial, mode.angular, argument, theta );

		int const order = ( mode.angular < 0 ) ? -mode.angular : mode.angular;
		double const radial = integerPower( argument, mode.radial );

		return ( mode.angular >= 0 ) ? radial*std::cos( order*theta )
		                             : radial*std::sin( order*theta );
	}

	void SurfaceFit::evaluateMode( FitMode mode, double argument, double theta,
	                               double &value, double &dArgument,
	                               double &dTheta ) const
	{
		if ( option.basis == FitBasis::Zernike )
		{
			value = zernike( mode.radial, mode.angular, argument, theta );
			dArgument = zernikeRadialDerivative( mode.radial, mode.angular,
			                                     argument, theta );
			dTheta = zernikeAngularDerivative( mode.radial, mode.angular,
			                                   argument, theta );
			return;
		}

		// The tensor-product control: argument^p times a bare Fourier mode, with
		// the ZernikeMode sign convention on m so that the two bases index the
		// same way. The p = 0 branch on the derivative is not an optimisation:
		// argument^-1 at argument = 0 is an infinity and 0 * infinity is a NaN,
		// which is the trap Zernike.cpp records at the same place.
		int const order = ( mode.angular < 0 ) ? -mode.angular : mode.angular;
		double const angle = order*theta;
		double const radial = integerPower( argument, mode.radial );
		double const radialPrime = ( mode.radial == 0 )
			? 0.0
			: mode.radial*integerPower( argument, mode.radial - 1 );

		if ( mode.angular >= 0 )
		{
			double const cosine = std::cos( angle );
			value = radial*cosine;
			dArgument = radialPrime*cosine;
			dTheta = -order*radial*std::sin( angle );
		}
		else
		{
			double const sine = std::sin( angle );
			value = radial*sine;
			dArgument = radialPrime*sine;
			dTheta = order*radial*std::cos( angle );
		}
	}

	void SurfaceFit::evaluateAll( double argument, double theta,
	                              double &r, double &z,
	                              double &dArgumentR, double &dArgumentZ,
	                              double &dThetaR, double &dThetaZ ) const
	{
		r = 0.0;
		z = 0.0;
		dArgumentR = 0.0;
		dArgumentZ = 0.0;
		dThetaR = 0.0;
		dThetaZ = 0.0;

		for ( std::size_t j = 0; j < modeList.size(); ++j )
		{
			double value = 0.0;
			double dArgument = 0.0;
			double dTheta = 0.0;
			evaluateMode( modeList[ j ], argument, theta, value, dArgument, dTheta );

			r += coefficientR[ j ]*value;
			z += coefficientZ[ j ]*value;
			dArgumentR += coefficientR[ j ]*dArgument;
			dArgumentZ += coefficientZ[ j ]*dArgument;
			dThetaR += coefficientR[ j ]*dTheta;
			dThetaZ += coefficientZ[ j ]*dTheta;
		}
	}

	void SurfaceFit::position( double normalisedFlux, double theta,
	                           double &r, double &z ) const
	{
		double dArgumentR = 0.0;
		double dArgumentZ = 0.0;
		double dThetaR = 0.0;
		double dThetaZ = 0.0;
		evaluateAll( argumentOf( normalisedFlux ), theta, r, z, dArgumentR,
		             dArgumentZ, dThetaR, dThetaZ );
	}

	void SurfaceFit::radialDerivative( double normalisedFlux, double theta,
	                                   double &r, double &z ) const
	{
		double value = 0.0;
		double height = 0.0;
		double dArgumentR = 0.0;
		double dArgumentZ = 0.0;
		double dThetaR = 0.0;
		double dThetaZ = 0.0;
		evaluateAll( argumentOf( normalisedFlux ), theta, value, height, dArgumentR,
		             dArgumentZ, dThetaR, dThetaZ );

		// d/d( discRadius ). For the disc coordinate the argument IS the disc
		// radius; for the control it is its square, so the chain factor is
		// 2 discRadius and applying it here is what makes the two settings
		// comparable in one column.
		double const chain = ( option.coordinate == FitRadialCoordinate::DiscRadius )
			? 1.0 : 2.0*discRadiusOf( normalisedFlux );

		r = chain*dArgumentR;
		z = chain*dArgumentZ;
	}

	void SurfaceFit::fluxDerivative( double normalisedFlux, double theta,
	                                 double &r, double &z ) const
	{
		double radialR = 0.0;
		double radialZ = 0.0;
		radialDerivative( normalisedFlux, theta, radialR, radialZ );

		// The 1/( 2 discRadius ) of Zernike.hpp, and then the discEdge that the
		// rescaling of the disc adds to it. No exception at the axis: this is
		// called per surface point and an infinity there is the coordinate, not
		// an error -- see the header.
		double const radius = discRadiusOf( normalisedFlux );
		r = fluxDerivativeFromRadial( radialR, radius )/option.discEdge;
		z = fluxDerivativeFromRadial( radialZ, radius )/option.discEdge;
	}

	void SurfaceFit::angularDerivative( double normalisedFlux, double theta,
	                                    double &r, double &z ) const
	{
		double value = 0.0;
		double height = 0.0;
		double dArgumentR = 0.0;
		double dArgumentZ = 0.0;
		evaluateAll( argumentOf( normalisedFlux ), theta, value, height, dArgumentR,
		             dArgumentZ, r, z );
	}

	double SurfaceFit::angularSpeed( double normalisedFlux, double theta ) const
	{
		double r = 0.0;
		double z = 0.0;
		angularDerivative( normalisedFlux, theta, r, z );

		return std::sqrt( r*r + z*z );
	}

	void SurfaceFit::axis( double &r, double &z ) const
	{
		axisAtAngle( 0.0, r, z );
	}

	void SurfaceFit::axisAtAngle( double theta, double &r, double &z ) const
	{
		position( 0.0, theta, r, z );
	}

	double SurfaceFit::coefficientEnvelope( int degreeWanted ) const
	{
		double envelope = 0.0;

		for ( std::size_t j = 0; j < modeList.size(); ++j )
		{
			int const order = ( modeList[ j ].angular < 0 )
				? -modeList[ j ].angular : modeList[ j ].angular;

			// The degree of a mode: l for Zernike, where the radial index IS the
			// degree, and p + |m| for the tensor product, where it is not.
			int const modeDegree = ( option.basis == FitBasis::Zernike )
				? modeList[ j ].radial : modeList[ j ].radial + order;

			if ( modeDegree != degreeWanted )
				continue;

			envelope = std::max( envelope, std::abs( coefficientR[ j ] ) );
			envelope = std::max( envelope, std::abs( coefficientZ[ j ] ) );
		}

		return envelope;
	}

	std::vector<double> const &SurfaceFit::majorRadiusCoefficients() const
	{
		return coefficientR;
	}

	std::vector<double> const &SurfaceFit::heightCoefficients() const
	{
		return coefficientZ;
	}

	void SurfaceFit::checkExpansionIsPlain( char const *where ) const
	{
		// A ZernikeExpansion is a degree and a coefficient vector and nothing
		// else, so it cannot say which basis, which radial coordinate or which
		// disc edge produced it. Handing one back from a fit that differs in any
		// of the three would hand back an object whose argument means something
		// other than sqrt( Psi_N ) with nothing anywhere to say so -- which is
		// the same species of silent mislabelling this file's header warns about
		// for the radial coordinate itself.
		if ( option.basis != FitBasis::Zernike
		     || option.coordinate != FitRadialCoordinate::DiscRadius
		     || option.discEdge != 1.0 )
			throw std::runtime_error(
				std::string( where ) + ": this fit is not a plain expansion in"
				" sqrt( Psi_N ) -- basis "
				+ fitBasisName( option.basis ) + ", coordinate "
				+ fitRadialCoordinateName( option.coordinate ) + ", disc edge "
				+ std::to_string( option.discEdge ) + ". A ZernikeExpansion"
				" records none of those, so handing one back from here would hand"
				" back a mislabelled object. Evaluate through this class instead" );
	}

	ZernikeExpansion SurfaceFit::majorRadiusExpansion() const
	{
		checkExpansionIsPlain( "SurfaceFit::majorRadiusExpansion" );

		return ZernikeExpansion( degree, coefficientR );
	}

	ZernikeExpansion SurfaceFit::heightExpansion() const
	{
		checkExpansionIsPlain( "SurfaceFit::heightExpansion" );

		return ZernikeExpansion( degree, coefficientZ );
	}

	int SurfaceFit::maxDegree() const
	{
		return degree;
	}

	SurfaceFitOptions const &SurfaceFit::options() const
	{
		return option;
	}

	SurfaceFitDiagnostics const &SurfaceFit::diagnostics() const
	{
		return diagnostic;
	}

}
