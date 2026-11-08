#include "meq/SafetyFactor.hpp"

#include "meq/StageChecksum.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace meq
{

	namespace
	{

		/// Same sign, treating either zero as "no".
		bool sameSign( double a, double b )
		{
			return ( a > 0.0 && b > 0.0 ) || ( a < 0.0 && b < 0.0 );
		}

		void requireAscending( std::vector<double> const &x, char const *what )
		{
			for ( std::size_t i = 1; i < x.size(); ++i )
			{
				if ( !( x[ i ] > x[ i - 1 ] ) )
				{
					std::ostringstream message;
					message << what << ": the labels must be strictly ascending, "
					           "and entry " << i << " is " << x[ i ]
					        << " after " << x[ i - 1 ]
					        << ". A repeated surface is what a family with two "
					           "traces at the same level gives, and the slope "
					           "rule below would divide by its zero interval.";
					throw std::invalid_argument( message.str() );
				}
			}
		}

	}

	/*
	 * FRITSCH-CARLSON, WITH THE STANDARD THREE-POINT ENDS.
	 *
	 * The interior rule is the weighted harmonic mean of the two neighbouring
	 * secants, which is zero whenever they disagree in sign -- that is exactly
	 * the no-overshoot property, and it is why a local wiggle in the extracted
	 * g^2 cannot become a sign change in gg'. The ends take the usual
	 * three-point formula and are then limited against the adjacent secant,
	 * because an unlimited endpoint is where a monotone rule most often is not.
	 */
	ToroidalFieldMap::ToroidalFieldMap( Solve solve, Options options )
		: solve( std::move( solve ) ), options( std::move( options ) )
	{
		if ( !this->solve )
			throw std::invalid_argument(
				"meq::ToroidalFieldMap: no solve was given, so there is nothing "
				"to take a Picard step with" );
		if ( !this->options.target )
			throw std::invalid_argument(
				"meq::ToroidalFieldMap: no target safety factor was given, so "
				"there is nothing to invert against" );
		if ( this->options.degree == 0 )
			throw std::invalid_argument(
				"meq::ToroidalFieldMap: a degree-zero g^2 is a constant, which "
				"makes gg' identically zero and the loop unable to move" );
	}

	std::vector<double> ToroidalFieldMap::refuse(
		std::vector<double> const &coefficients ) const
	{
		/*
		 * THE RESIDUAL AT AN INADMISSIBLE POINT, AND WHY IT IS THIS ONE.
		 *
		 * `G( c ) - c` has to be finite here and it has to point somewhere the
		 * line search can go. Returning the last admissible coefficients makes
		 * the residual `lastGood - c`, which is large wherever `c` has wandered
		 * and is exactly zero only at `lastGood` itself -- and `lastGood` is a
		 * point the map really was evaluated at, so it is not a root unless it
		 * is the answer. Before ANY step has succeeded there is no such point,
		 * and `c - 1` is used instead: a constant nonzero residual, which is
		 * the honest statement that nothing is known yet.
		 */
		std::vector<double> out( coefficients.size() );
		for ( std::size_t i = 0; i < coefficients.size(); ++i )
			out[ i ] = lastGood.size() == coefficients.size()
			           ? lastGood[ i ] : coefficients[ i ] - 1.0;
		return out;
	}

	std::vector<double> ToroidalFieldMap::operator()(
		std::vector<double> const &coefficients )
	{
		// ADMISSIBLE MEANS g^2 > 0 EVERYWHERE THE SOURCE WILL BE ASKED, which
		// is the whole of Psi and not the family's cut -- the knots span
		// [ 0, 1 ] and the profile is evaluated at every quadrature point of
		// the domain, including outside the plasma.
		std::size_t const checks = 64;
		for ( std::size_t i = 0; i <= checks; ++i )
		{
			double const psi = static_cast<double>( i )
			                   /static_cast<double>( checks );
			double value = 0.0, power = 1.0;
			for ( double c : coefficients )
			{
				value += c*power;
				power *= psi;
			}
			if ( !( value > 0.0 ) || !std::isfinite( value ) )
			{
				++refusalCount;
				return refuse( coefficients );
			}
		}

		++solveCount;
		// THE OUTERMOST SEAM OF THE q-DRIVEN LOOP, PROBED. `G( c ) - c` is what
		// the outer Newton roots, so a device run reporting `final norm 0` in
		// zero iterations is saying `G( c0 ) == c0` -- and these two lines say
		// whether that is true of the map or of something upstream of it.
		// MEASUREMENTS.md M-205.
		stageChecksum( "toroidalMap.in", coefficients );
		FluxSurfaceFamily const *family = solve( ggPrimeKnotsFromCoefficients(
			coefficients, options.knotSamples, options.extension ) );
		if ( family == nullptr || family->empty() )
		{
			++refusalCount;
			stageChecksum( "toroidalMap.refusedEmptyFamily", coefficients );
			return refuse( coefficients );
		}

		// THE REFLECTION, ONCE. The family labels by Psi_N and the target is
		// written in the source's Psi; see the header.
		auto const &target = options.target;
		ToroidalField const inverted = invertSafetyFactor(
			*family,
			[ &target ]( double normalisedFluxFamily )
			{
				return target( 1.0 - normalisedFluxFamily );
			} );

		lastGood = coefficients;
		labels = inverted.normalisedFlux;
		std::vector<double> const fitted =
			fitToroidalFieldSquared( inverted, options.degree );
		lastFloor = toroidalMapFloor( inverted, options.degree,
		                              options.solveAccuracy );
		stageChecksum( "toroidalMap.out", fitted );
		return fitted;
	}

	PinnedToroidalFieldMap::PinnedToroidalFieldMap( ToroidalFieldMap &map,
	                                                unsigned int degree, Pin pin )
		: map( map ), degree( degree ), pin( std::move( pin ) )
	{
		if ( degree == 0 )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap: a degree-zero g^2 has no q "
				"condition to give way" );
		if ( !std::isfinite( this->pin.value ) || this->pin.value == 0.0 )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap: the pinned value must be finite "
				"and non-zero" );
		if ( this->pin.kind == ScalePin::EdgeField && !( this->pin.value > 0.0 ) )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap: g on the boundary is pinned as a "
				"magnitude, so it must be positive; its sign is q's (see "
				"section 3)" );
		if ( this->pin.kind != ScalePin::EdgeField && !this->pin.measured )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap: a current or flux pin is a row of "
				"the outer system and needs the solved state's own value" );
		if ( !( this->pin.weight > 0.0 ) || !std::isfinite( this->pin.weight ) )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap: the pin row's weight must be "
				"positive" );
	}

	std::vector<double> PinnedToroidalFieldMap::unknowns(
		std::vector<double> const &coefficients ) const
	{
		if ( coefficients.size() != degree + 1 )
			throw std::invalid_argument(
				"meq::PinnedToroidalFieldMap::unknowns: the coefficients are not "
				"of this degree" );
		if ( pin.kind == ScalePin::EdgeField )
			return std::vector<double>( coefficients.begin() + 1, coefficients.end() );
		return coefficients;
	}

	std::vector<double> PinnedToroidalFieldMap::coefficients(
		std::vector<double> const &unknowns ) const
	{
		if ( pin.kind != ScalePin::EdgeField )
			return unknowns;
		std::vector<double> out;
		out.reserve( unknowns.size() + 1 );
		out.push_back( pin.value*pin.value );
		out.insert( out.end(), unknowns.begin(), unknowns.end() );
		return out;
	}

	std::vector<double> PinnedToroidalFieldMap::operator()(
		std::vector<double> const &unknowns )
	{
		/*
		 * A REFUSAL POINTS BACK AT THE LAST ADMISSIBLE POINT, as the free
		 * map's does: finite, and zero nowhere but at a point that really was
		 * evaluated.
		 *
		 * **AND IT IS NEVER SMALLER THAN ANY RESIDUAL THE MAP HAS RETURNED.**
		 * `lastGood - x` alone is large only where the iterate has WANDERED:
		 * when the last admissible point is a differencing column a step away,
		 * a refused trial point beside it reads a residual smaller than the
		 * iterate's own, the line search's Armijo test accepts it, and the
		 * loop then sits on a point the map cannot evaluate. Measured on
		 * QDrivenScalePin's base run in single: accepted at 0.16 from 0.50,
		 * stalled there at 0.14 (MEASUREMENTS.md M-358). Scaled to twice the
		 * largest admissible 2-norm seen, a refusal fails every sufficient-
		 * decrease test the loop can pose.
		 */
		auto refuse = [ & ]
		{
			lastWasRefused = true;
			std::vector<double> out( unknowns.size() );
			std::vector<double> back( unknowns.size() );
			double length = 0.0;
			for ( std::size_t i = 0; i < unknowns.size(); ++i )
			{
				back[ i ] = lastGood.size() == unknowns.size()
				            ? lastGood[ i ] - unknowns[ i ] : -1.0;
				length += back[ i ]*back[ i ];
			}
			length = std::sqrt( length );
			double const scale = length > 0.0
				? std::max( 1.0, 2.0*largestResidual/length ) : 1.0;
			for ( std::size_t i = 0; i < unknowns.size(); ++i )
				out[ i ] = unknowns[ i ] + scale*back[ i ];
			return out;
		};
		lastWasRefused = false;

		std::vector<double> const c = coefficients( unknowns );
		int const refusedBefore = map.refusals();
		std::vector<double> const image = map( c );
		if ( map.refusals() != refusedBefore )
			return refuse();

		double measured = 0.0;
		if ( pin.kind != ScalePin::EdgeField )
		{
			try
			{
				measured = pin.measured();
			}
			catch ( ... )
			{
				return refuse();
			}
			if ( !std::isfinite( measured ) )
				return refuse();
		}

		/*
		 * THE FIT'S ORTHONORMAL COORDINATES, on the labels this step was fitted
		 * on -- the same Vandermonde in the SOURCE's Psi that
		 * fitToroidalFieldSquared() solves, so `y = R c` is the coefficient
		 * vector's own discrete L2 coordinate over the surfaces.
		 */
		std::vector<double> const &label = map.lastLabels();
		Eigen::Index const terms = static_cast<Eigen::Index>( degree ) + 1;
		Eigen::MatrixXd design( static_cast<Eigen::Index>( label.size() ), terms );
		for ( std::size_t i = 0; i < label.size(); ++i )
		{
			double power = 1.0;
			for ( Eigen::Index j = 0; j < terms; ++j )
			{
				design( static_cast<Eigen::Index>( i ), j ) = power;
				power *= 1.0 - label[ i ];
			}
		}
		Eigen::HouseholderQR<Eigen::MatrixXd> const qr( design );
		Eigen::MatrixXd const r =
			Eigen::MatrixXd( qr.matrixQR().topRows( terms ).triangularView<Eigen::Upper>() )
			/std::sqrt( static_cast<double>( label.size() ) );

		Eigen::VectorXd difference( terms );
		for ( Eigen::Index j = 0; j < terms; ++j )
			difference( j ) = image[ static_cast<std::size_t>( j ) ]
			                  - c[ static_cast<std::size_t>( j ) ];
		Eigen::VectorXd const y = r*difference;

		// y_d GIVES WAY: rows y_0 .. y_{d-1}, and the pin's row after them.
		std::vector<double> residual( static_cast<std::size_t>( terms - 1 ) );
		for ( Eigen::Index j = 0; j + 1 < terms; ++j )
			residual[ static_cast<std::size_t>( j ) ] = y( j );
		if ( pin.kind != ScalePin::EdgeField )
			residual.push_back( pin.weight*( measured - pin.value )/pin.value );

		/*
		 * THE FLOOR. Each fit row is a combination of the map's coefficients
		 * with the weights of R's row, so its noise is at most the row sum
		 * times the map's floor; the pin row carries the solve's relative
		 * accuracy with toroidalMapFloor()'s four standard deviations of four.
		 */
		double rowSum = 0.0;
		for ( Eigen::Index i = 0; i + 1 < terms; ++i )
			rowSum = std::max( rowSum, r.row( i ).cwiseAbs().sum() );
		lastFloor = rowSum*map.floor();
		if ( pin.kind != ScalePin::EdgeField )
			lastFloor = std::max( lastFloor, pin.weight*16.0*map.solveAccuracy() );

		double norm = 0.0;
		for ( double value : residual )
			norm += value*value;
		largestResidual = std::max( largestResidual, std::sqrt( norm ) );
		lastGood = unknowns;
		measuredValue = measured;
		std::vector<double> out( unknowns.size() );
		for ( std::size_t i = 0; i < unknowns.size(); ++i )
			out[ i ] = unknowns[ i ] + residual[ i ];
		return out;
	}

	std::vector<double> monotoneSlopes( std::vector<double> const &x,
	                                    std::vector<double> const &y )
	{
		if ( x.size() != y.size() )
			throw std::invalid_argument(
				"meq::monotoneSlopes: the abscissae and the values differ in "
				"length" );
		if ( x.size() < 2 )
			throw std::invalid_argument(
				"meq::monotoneSlopes: at least two points are needed to have a "
				"slope at all" );
		requireAscending( x, "meq::monotoneSlopes" );

		std::size_t const n = x.size();
		std::vector<double> h( n - 1 ), delta( n - 1 );
		for ( std::size_t i = 0; i + 1 < n; ++i )
		{
			h[ i ] = x[ i + 1 ] - x[ i ];
			delta[ i ] = ( y[ i + 1 ] - y[ i ] )/h[ i ];
		}

		std::vector<double> slope( n, 0.0 );

		if ( n == 2 )
		{
			slope[ 0 ] = delta[ 0 ];
			slope[ 1 ] = delta[ 0 ];
			return slope;
		}

		for ( std::size_t i = 1; i + 1 < n; ++i )
		{
			if ( !sameSign( delta[ i - 1 ], delta[ i ] ) )
			{
				// A local extremum of the data. Zero here is what makes the
				// interpolant monotone on each side of it and is the whole rule.
				slope[ i ] = 0.0;
				continue;
			}
			double const w1 = 2.0*h[ i ] + h[ i - 1 ];
			double const w2 = h[ i ] + 2.0*h[ i - 1 ];
			slope[ i ] = ( w1 + w2 )/( w1/delta[ i - 1 ] + w2/delta[ i ] );
		}

		auto endSlope = []( double hNear, double hFar, double deltaNear,
		                    double deltaFar )
		{
			double value = ( ( 2.0*hNear + hFar )*deltaNear - hNear*deltaFar )
			               /( hNear + hFar );
			if ( !sameSign( value, deltaNear ) )
				value = 0.0;
			else if ( !sameSign( deltaNear, deltaFar )
			          && std::abs( value ) > 3.0*std::abs( deltaNear ) )
				value = 3.0*deltaNear;
			return value;
		};

		slope[ 0 ] = endSlope( h[ 0 ], h[ 1 ], delta[ 0 ], delta[ 1 ] );
		slope[ n - 1 ] = endSlope( h[ n - 2 ], h[ n - 3 ], delta[ n - 2 ],
		                           delta[ n - 3 ] );
		return slope;
	}

	ToroidalField invertSafetyFactor(
		FluxSurfaceFamily const &family,
		std::function<double( double normalisedFlux )> const &target )
	{
		if ( family.empty() )
			throw std::invalid_argument(
				"meq::invertSafetyFactor: the family has no surfaces, so there "
				"is no geometry to invert against" );
		if ( !target )
			throw std::invalid_argument(
				"meq::invertSafetyFactor: no target safety factor was given" );

		double const fourPiSquared = 4.0*M_PI*M_PI;

		ToroidalField out;
		out.normalisedFlux.reserve( family.surfaces.size() );
		out.safetyFactor.reserve( family.surfaces.size() );
		out.g.reserve( family.surfaces.size() );
		out.gSquared.reserve( family.surfaces.size() );

		for ( std::size_t i = 0; i < family.surfaces.size(); ++i )
		{
			FluxSurface const &surface = family.surfaces[ i ];

			// POSITIVE BY CONSTRUCTION, so a non-positive one is a corrupt
			// family and not a hard case. Dividing by it returns a signed
			// infinity, which the spline would carry into the source and the
			// solve would meet as a NaN several frames away from here.
			if ( !( surface.vPrime > 0.0 ) || !( surface.inverseRSquared > 0.0 ) )
			{
				std::ostringstream message;
				message << "meq::invertSafetyFactor: surface " << i
				        << " at Psi_N = " << surface.normalisedFlux
				        << " has V' = " << surface.vPrime
				        << " and < R^-2 > = " << surface.inverseRSquared
				        << ". Both are positive by construction -- V' is an "
				           "integral of a positive integrand and < R^-2 > is an "
				           "average of one -- so this is a corrupt family "
				           "rather than a difficult equilibrium.";
				throw std::invalid_argument( message.str() );
			}

			double const q = target( surface.normalisedFlux );
			if ( !std::isfinite( q ) )
			{
				std::ostringstream message;
				message << "meq::invertSafetyFactor: the target safety factor "
				           "at Psi_N = " << surface.normalisedFlux << " is "
				        << q << ", which is not finite";
				throw std::invalid_argument( message.str() );
			}

			double const g = fourPiSquared*q
			                 /( surface.vPrime*surface.inverseRSquared );

			out.normalisedFlux.push_back( surface.normalisedFlux );
			out.safetyFactor.push_back( q );
			out.g.push_back( g );
			out.gSquared.push_back( g*g );
		}

		return out;
	}

	std::vector<Knot> ggPrimeKnots( ToroidalField const &field,
	                                EdgeExtension extension )
	{
		if ( field.size() < 2 )
			throw std::invalid_argument(
				"meq::ggPrimeKnots: at least two surfaces are needed to "
				"differentiate g^2 at all" );
		requireAscending( field.normalisedFlux, "meq::ggPrimeKnots" );

		/*
		 * THE REFLECTION, AND IT IS THE WHOLE REASON THIS IS A FUNCTION.
		 *
		 * The field arrives against the family's Psi_N, zero on the axis; the
		 * source reads Psi = 1 - Psi_N, one on the axis. So
		 *
		 *     gg'( Psi ) = ( 1/2 ) d( g^2 )/dPsi
		 *                = -( 1/2 ) d( g^2 )/dPsi_N
		 *
		 * and the MINUS is the half that does not announce itself: without it
		 * the solve converges, at full order, to an equilibrium with its shear
		 * reversed -- which is a configuration a real machine can have, so
		 * nothing downstream looks wrong. See the header, section 2.
		 */
		std::vector<double> const slopeOfGSquared =
			monotoneSlopes( field.normalisedFlux, field.gSquared );

		std::vector<double> ggPrime( field.size() );
		for ( std::size_t i = 0; i < field.size(); ++i )
			ggPrime[ i ] = -0.5*slopeOfGSquared[ i ];

		// The knot derivative is d( gg' )/dPsi, so the rule is applied a second
		// time and reflected a second time.
		std::vector<double> const slopeOfGGPrime =
			monotoneSlopes( field.normalisedFlux, ggPrime );

		std::vector<Knot> knots;
		knots.reserve( field.size() + 1 );
		for ( std::size_t i = 0; i < field.size(); ++i )
		{
			Knot knot;
			knot.psi = 1.0 - field.normalisedFlux[ i ];
			knot.value = ggPrime[ i ];
			knot.derivative = -slopeOfGGPrime[ i ];
			knots.push_back( knot );
		}

		// Psi descends as Psi_N ascends, and SplineProfile wants ascending.
		std::reverse( knots.begin(), knots.end() );

		if ( extension == EdgeExtension::VacuumOutside )
		{
			// Psi = 0 is the plasma edge. Outside it the toroidal field is the
			// vacuum one, g = const, so gg' = 0 -- and a knot there is what
			// makes SplineProfile's clamp carry that rather than the innermost
			// value it happens to hold. Only added when the family does not
			// already reach the edge, which it does not by construction.
			if ( knots.front().psi > 0.0 )
			{
				Knot vacuum;
				vacuum.psi = 0.0;
				vacuum.value = 0.0;
				vacuum.derivative = 0.0;
				knots.insert( knots.begin(), vacuum );
			}
		}

		return knots;
	}


	/*
	 * LEAST SQUARES ON THE MONOMIALS, WHICH IS ADEQUATE HERE AND WOULD NOT BE
	 * AT HIGH DEGREE. Psi is in [ 0, 1 ] and the degree is two or three, so the
	 * Vandermonde's conditioning is a non-issue; a Jacobi or Chebyshev basis
	 * would be the answer at ten. Eigen's column-pivoting QR is used rather than
	 * the normal equations because it costs nothing here and squares nothing.
	 */
	std::vector<double> fitToroidalFieldSquared( ToroidalField const &field,
	                                             unsigned int degree )
	{
		if ( degree == 0 )
			throw std::invalid_argument(
				"meq::fitToroidalFieldSquared: a degree of zero is a constant "
				"g^2, i.e. gg' identically zero, which is a vacuum toroidal "
				"field and not a fit" );
		std::size_t const terms = degree + 1;
		if ( field.size() < terms )
		{
			std::ostringstream message;
			message << "meq::fitToroidalFieldSquared: degree " << degree
			        << " needs " << terms << " coefficients and the family has "
			        << field.size() << " surfaces";
			throw std::invalid_argument( message.str() );
		}
		requireAscending( field.normalisedFlux, "meq::fitToroidalFieldSquared" );

		Eigen::MatrixXd design( field.size(), terms );
		Eigen::VectorXd rhs( field.size() );
		for ( std::size_t i = 0; i < field.size(); ++i )
		{
			// THE SOURCE'S Psi, not the family's. See the header, section 2.
			double const psi = 1.0 - field.normalisedFlux[ i ];
			double power = 1.0;
			for ( std::size_t j = 0; j < terms; ++j )
			{
				design( static_cast<Eigen::Index>( i ),
				        static_cast<Eigen::Index>( j ) ) = power;
				power *= psi;
			}
			rhs( static_cast<Eigen::Index>( i ) ) = field.gSquared[ i ];
		}

		Eigen::VectorXd const solution =
			design.colPivHouseholderQr().solve( rhs );

		std::vector<double> out( terms );
		for ( std::size_t j = 0; j < terms; ++j )
			out[ j ] = solution( static_cast<Eigen::Index>( j ) );
		return out;
	}

	double toroidalMapFloor( ToroidalField const &field, unsigned int degree,
	                         double solveAccuracy )
	{
		if ( degree == 0 )
			throw std::invalid_argument(
				"meq::toroidalMapFloor: a degree of zero is not a fit" );
		std::size_t const terms = degree + 1;
		if ( field.size() < terms || field.gSquared.size() != field.size() )
			throw std::invalid_argument(
				"meq::toroidalMapFloor: the family has fewer surfaces than the "
				"fit has coefficients, or no g^2 on each of them" );
		requireAscending( field.normalisedFlux, "meq::toroidalMapFloor" );

		// THE SAME DESIGN MATRIX AS fitToroidalFieldSquared(), in the source's
		// Psi, so P below is that fit's own pseudo-inverse and not a model of
		// it.
		Eigen::MatrixXd design( field.size(), terms );
		for ( std::size_t i = 0; i < field.size(); ++i )
		{
			double const psi = 1.0 - field.normalisedFlux[ i ];
			double power = 1.0;
			for ( std::size_t j = 0; j < terms; ++j )
			{
				design( static_cast<Eigen::Index>( i ),
				        static_cast<Eigen::Index>( j ) ) = power;
				power *= psi;
			}
		}
		Eigen::MatrixXd const pseudoInverse =
			design.completeOrthogonalDecomposition().pseudoInverse();

		// Row k of P, each surface's entry weighted by that surface's g^2 --
		// the relative error is per surface, so the absolute one is g^2 times
		// it -- and summed in quadrature, the surfaces' errors being
		// independent. The worst coefficient is the one the max norm sees.
		double worst = 0.0;
		for ( Eigen::Index k = 0; k < pseudoInverse.rows(); ++k )
		{
			double sum = 0.0;
			for ( Eigen::Index i = 0; i < pseudoInverse.cols(); ++i )
			{
				double const term = pseudoInverse( k, i )
				                    *field.gSquared[ static_cast<std::size_t>( i ) ];
				sum += term*term;
			}
			worst = std::max( worst, std::sqrt( sum ) );
		}

		// 4 for g^2's two averages, each squared; 4 standard deviations.
		return 4.0*4.0*std::abs( solveAccuracy )*worst;
	}

	std::vector<Knot> ggPrimeKnotsFromFit( ToroidalField const &field,
	                                       unsigned int degree,
	                                       std::size_t samples,
	                                       EdgeExtension extension )
	{
		return ggPrimeKnotsFromCoefficients(
			fitToroidalFieldSquared( field, degree ), samples, extension );
	}

	std::vector<Knot> ggPrimeKnotsFromCoefficients(
		std::vector<double> const &c, std::size_t samples,
		EdgeExtension extension )
	{
		if ( samples < 2 )
			throw std::invalid_argument(
				"meq::ggPrimeKnotsFromCoefficients: at least two knots are "
				"needed" );
		if ( c.empty() )
			throw std::invalid_argument(
				"meq::ggPrimeKnotsFromCoefficients: no coefficients, so there is "
				"no g^2 to differentiate" );

		// gg' = ( 1/2 ) dG/dPsi and its derivative, both exact on the
		// polynomial -- which is the whole reason for fitting rather than
		// interpolating. No slope rule, no reflection: the fit was built in the
		// source's Psi already.
		auto ggPrime = [ &c ]( double psi )
		{
			double value = 0.0, power = 1.0;
			for ( std::size_t j = 1; j < c.size(); ++j )
			{
				value += 0.5*static_cast<double>( j )*c[ j ]*power;
				power *= psi;
			}
			return value;
		};
		auto ggPrimeSlope = [ &c ]( double psi )
		{
			double value = 0.0, power = 1.0;
			for ( std::size_t j = 2; j < c.size(); ++j )
			{
				value += 0.5*static_cast<double>( j )
				         *static_cast<double>( j - 1 )*c[ j ]*power;
				power *= psi;
			}
			return value;
		};

		/*
		 * THE KNOTS SPAN THE WHOLE OF Psi, NOT THE FAMILY'S RANGE, AND THE
		 * DISTINCTION IS THE DIFFERENCE BETWEEN A MODEL AND A MEASUREMENT.
		 *
		 * meq::FluxSurfaceFamily refuses to extrapolate, rightly: outside the
		 * cut a traced surface is made of something else and only the band mask
		 * can say so. A FIT is not a traced surface. It is a low-order model of
		 * g^2 and evaluating it at Psi = 0 or 1 is what a model is for.
		 *
		 * SAMPLING ONLY [ 0.05, 0.95 ] AND LETTING SplineProfile CLAMP IS A
		 * DIFFERENT PROFILE, AND IT COSTS THE FIXED POINT. Measured: with the
		 * knots on the family's range, the outer residual at the CLOSED-FORM
		 * answer is 1.3e-01 rather than zero -- because the clamped profile
		 * differs from the closed form exactly where it clamps, near the axis
		 * and near the edge, so the equilibrium it produces is a different one.
		 * Both outer methods then converge, correctly, to a fixed point that is
		 * not the answer, and no solver could have found the answer because it
		 * was not a root of the map they were given.
		 */
		double const lowest = 0.0;
		double const highest = 1.0;

		std::vector<Knot> knots;
		knots.reserve( samples + 1 );
		for ( std::size_t i = 0; i < samples; ++i )
		{
			double const t = static_cast<double>( i )
			                 /static_cast<double>( samples - 1 );
			double const psi = lowest + t*( highest - lowest );
			knots.push_back( Knot{ psi, ggPrime( psi ), ggPrimeSlope( psi ) } );
		}

		if ( extension == EdgeExtension::VacuumOutside && knots.front().psi > 0.0 )
			knots.insert( knots.begin(), Knot{ 0.0, 0.0, 0.0 } );

		return knots;
	}

}
