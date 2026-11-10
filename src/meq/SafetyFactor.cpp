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


	std::size_t toroidalKnotCount( double cells, int order )
	{
		if ( !( cells > 0.0 ) || !std::isfinite( cells ) )
			throw std::invalid_argument(
				"meq::toroidalKnotCount: the plasma must cover a positive number of cells" );
		if ( order < 1 )
			throw std::invalid_argument(
				"meq::toroidalKnotCount: the polynomial degree must be at least one" );
		// Elements across the minor radius of a disc of this many cells.
		double const across = std::sqrt( cells/M_PI );
		// The cubic's O( Delta^4 ) keeps pace with O( h^( k + 1 ) ) to k = 3;
		// above, Delta ~ h^( ( k + 1 )/4 ).
		double const exponent = std::max( 1.0, ( order + 1.0 )/4.0 );
		double const knots = std::ceil( 2.0*std::pow( across, exponent ) - 1e-9 );
		return static_cast<std::size_t>( std::clamp( knots, 6.0, 128.0 ) );
	}

	ToroidalKnotCut toroidalKnotCut( std::size_t knots )
	{
		if ( knots < 4 )
			throw std::invalid_argument(
				"meq::toroidalKnotCut: a not-a-knot cubic needs four knots at least" );
		double const spacing = 1.0/static_cast<double>( knots + 1 );
		ToroidalKnotCut cut;
		cut.surfaces = knots;
		cut.innerCut = spacing*spacing;
		cut.outerCut = ( 1.0 - spacing )*( 1.0 - spacing );
		return cut;
	}

	std::vector<double> toroidalKnotLabels( std::size_t knots )
	{
		// fluxSurfaceLabels()'s arithmetic, so that the labels are the
		// family's to the bit: rho uniform between the square roots of the cut.
		ToroidalKnotCut const cut = toroidalKnotCut( knots );
		double const inner = std::sqrt( cut.innerCut );
		double const outer = std::sqrt( cut.outerCut );
		double const last = static_cast<double>( knots - 1 );
		std::vector<double> labels;
		for ( std::size_t i = 0; i < knots; ++i )
		{
			double const rho = inner + ( static_cast<double>( i )/last )*( outer - inner );
			labels.push_back( rho*rho );
		}
		return labels;
	}

	namespace
	{
		/// gg' = ( g^2 )'/2 of a ToroidalFieldSpline evaluator.
		class SplineGGPrime : public Profile
		{
			public:
				explicit SplineGGPrime( std::function<double( double, int )> s )
					: spline( std::move( s ) ) {}
				double operator()( double psi ) const override { return 0.5*spline( psi, 1 ); }
				double prime( double psi ) const override { return 0.5*spline( psi, 2 ); }
				double doublePrime( double psi ) const override { return 0.5*spline( psi, 3 ); }
			private:
				std::function<double( double, int )> spline;
		};
	}

	/*
	 * THE NOT-A-KNOT CUBIC, AS A LINEAR MAP FROM VALUES TO SECOND DERIVATIVES.
	 *
	 * The second derivatives M_i at the nodes solve the classical tridiagonal
	 * system inside and the not-a-knot conditions -- the third derivative
	 * continuous across the second and the second-to-last node -- at the
	 * ends. That system depends on the nodes only, so its inverse applied to
	 * the right-hand side's operator is kept once, n x n: the cardinal
	 * functions are its columns, and a change of values is a matrix-vector
	 * product. n is the knot count, at most 128.
	 */
	ToroidalFieldSpline::ToroidalFieldSpline( std::vector<double> labels )
		: knotLabels( std::move( labels ) )
	{
		std::size_t const n = knotLabels.size();
		if ( n < 4 )
			throw std::invalid_argument(
				"meq::ToroidalFieldSpline: a not-a-knot cubic needs four knots at least" );
		requireAscending( knotLabels, "meq::ToroidalFieldSpline" );
		if ( !( knotLabels.front() > 0.0 ) || !( knotLabels.back() < 1.0 ) )
			throw std::invalid_argument(
				"meq::ToroidalFieldSpline: the knots are flux surfaces, inside ( 0, 1 ) in Psi_N" );

		psi.resize( n );
		for ( std::size_t i = 0; i < n; ++i )
			psi[ i ] = 1.0 - knotLabels[ n - 1 - i ];

		Eigen::MatrixXd a = Eigen::MatrixXd::Zero( static_cast<Eigen::Index>( n ),
		                                           static_cast<Eigen::Index>( n ) );
		Eigen::MatrixXd d = a;
		auto h = [ this ]( std::size_t i ) { return psi[ i + 1 ] - psi[ i ]; };
		auto at = []( std::size_t i ) { return static_cast<Eigen::Index>( i ); };
		a( 0, 0 ) = -1.0/h( 0 );
		a( 0, 1 ) = 1.0/h( 0 ) + 1.0/h( 1 );
		a( 0, 2 ) = -1.0/h( 1 );
		for ( std::size_t i = 1; i + 1 < n; ++i )
		{
			a( at( i ), at( i - 1 ) ) = h( i - 1 );
			a( at( i ), at( i ) ) = 2.0*( h( i - 1 ) + h( i ) );
			a( at( i ), at( i + 1 ) ) = h( i );
			d( at( i ), at( i + 1 ) ) += 6.0/h( i );
			d( at( i ), at( i ) ) += -6.0/h( i ) - 6.0/h( i - 1 );
			d( at( i ), at( i - 1 ) ) += 6.0/h( i - 1 );
		}
		std::size_t const l = n - 1;
		a( at( l ), at( l - 2 ) ) = -1.0/h( l - 2 );
		a( at( l ), at( l - 1 ) ) = 1.0/h( l - 2 ) + 1.0/h( l - 1 );
		a( at( l ), at( l ) ) = -1.0/h( l - 1 );
		Eigen::MatrixXd const map = a.colPivHouseholderQr().solve( d );
		valuesToSeconds.resize( n*n );
		for ( std::size_t i = 0; i < n; ++i )
			for ( std::size_t j = 0; j < n; ++j )
				valuesToSeconds[ i*n + j ] = map( at( i ), at( j ) );

		setValues( std::vector<double>( n, 1.0 ) );
	}

	void ToroidalFieldSpline::setValues( std::vector<double> values )
	{
		std::size_t const n = knotLabels.size();
		if ( values.size() != n )
			throw std::invalid_argument(
				"meq::ToroidalFieldSpline::setValues: one value per knot" );
		for ( double v : values )
			if ( !std::isfinite( v ) )
				throw std::invalid_argument(
					"meq::ToroidalFieldSpline::setValues: a value is not finite" );
		knotValues = std::move( values );
		sourceValues.resize( n );
		for ( std::size_t i = 0; i < n; ++i )
			sourceValues[ i ] = knotValues[ n - 1 - i ];
		seconds.assign( n, 0.0 );
		for ( std::size_t i = 0; i < n; ++i )
			for ( std::size_t k = 0; k < n; ++k )
				seconds[ i ] += valuesToSeconds[ i*n + k ]*sourceValues[ k ];
	}

	double ToroidalFieldSpline::evaluate( std::vector<double> const &y,
	                                      std::vector<double> const &m,
	                                      double t, int derivative ) const
	{
		std::size_t const n = psi.size();
		std::size_t i = 0;
		if ( t >= psi[ n - 2 ] )
			i = n - 2;
		else if ( t > psi[ 0 ] )
			i = static_cast<std::size_t>(
				    std::upper_bound( psi.begin(), psi.end(), t ) - psi.begin() ) - 1;
		double const h = psi[ i + 1 ] - psi[ i ];
		double const a = psi[ i + 1 ] - t, b = t - psi[ i ];
		double const ca = y[ i ]/h - m[ i ]*h/6.0, cb = y[ i + 1 ]/h - m[ i + 1 ]*h/6.0;
		switch ( derivative )
		{
			case 0:
				return m[ i ]*a*a*a/( 6.0*h ) + m[ i + 1 ]*b*b*b/( 6.0*h ) + ca*a + cb*b;
			case 1:
				return -m[ i ]*a*a/( 2.0*h ) + m[ i + 1 ]*b*b/( 2.0*h ) - ca + cb;
			case 2:
				return m[ i ]*a/h + m[ i + 1 ]*b/h;
			default:
				return ( m[ i + 1 ] - m[ i ] )/h;
		}
	}

	double ToroidalFieldSpline::operator()( double at, int derivative ) const
	{
		return evaluate( sourceValues, seconds, at, derivative );
	}

	double ToroidalFieldSpline::cardinal( std::size_t j, double at, int derivative ) const
	{
		std::size_t const n = psi.size();
		if ( j >= n )
			throw std::out_of_range( "meq::ToroidalFieldSpline::cardinal: no such knot" );
		// Knot j in the labels' order is node n - 1 - j in the source's Psi.
		std::size_t const node = n - 1 - j;
		std::vector<double> y( n, 0.0 ), m( n, 0.0 );
		y[ node ] = 1.0;
		for ( std::size_t i = 0; i < n; ++i )
			m[ i ] = valuesToSeconds[ i*n + node ];
		return evaluate( y, m, at, derivative );
	}

	std::shared_ptr<Profile const> ToroidalFieldSpline::ggPrime() const
	{
		auto copy = std::make_shared<ToroidalFieldSpline const>( *this );
		return std::make_shared<SplineGGPrime const>(
			[ copy ]( double at, int derivative ) { return ( *copy )( at, derivative ); } );
	}

	std::shared_ptr<Profile const> ToroidalFieldSpline::cardinalGGPrime( std::size_t j ) const
	{
		if ( j >= size() )
			throw std::out_of_range( "meq::ToroidalFieldSpline::cardinalGGPrime: no such knot" );
		auto copy = std::make_shared<ToroidalFieldSpline const>( *this );
		return std::make_shared<SplineGGPrime const>(
			[ copy, j ]( double at, int derivative ) { return copy->cardinal( j, at, derivative ); } );
	}

	bool ToroidalFieldSpline::positive( int samples ) const
	{
		for ( int k = 0; k <= samples; ++k )
		{
			double const value = ( *this )( static_cast<double>( k )/samples );
			if ( !( value > 0.0 ) || !std::isfinite( value ) )
				return false;
		}
		for ( double v : knotValues )
			if ( !( v > 0.0 ) )
				return false;
		return true;
	}

	SafetyFactorSystem::SafetyFactorSystem( std::vector<double> labels,
	                                        std::function<double( double psi )> const &targetIn,
	                                        Pin pinIn, Solve solveIn, Differentiate differentiateIn,
	                                        double solveAccuracy )
		: spline( std::move( labels ) ), pin( pinIn ), solve( std::move( solveIn ) ),
		  differentiate( std::move( differentiateIn ) ), accuracy( solveAccuracy )
	{
		if ( !solve || !differentiate )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem: a solve and its derivative are both needed" );
		if ( !targetIn )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem: no target safety factor was given" );
		if ( !std::isfinite( pin.value ) || pin.value == 0.0 )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem: the pinned value must be finite and non-zero" );
		if ( pin.kind == ScalePin::EdgeField && !( pin.value > 0.0 ) )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem: g on the boundary is pinned as a magnitude, so it "
				"must be positive; its sign is q's (SafetyFactor.hpp, section 3)" );
		if ( !( accuracy >= 0.0 ) )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem: the solve's accuracy must be non-negative" );
		// THE REFLECTION, ONCE: the target is in the source's Psi.
		for ( double label : spline.labels() )
		{
			double const q = targetIn( 1.0 - label );
			if ( !std::isfinite( q ) || q == 0.0 )
			{
				std::ostringstream message;
				message << "meq::SafetyFactorSystem: the target safety factor at Psi_N = "
				        << label << " is " << q << "; a row q/q* - 1 needs it finite and "
				           "non-zero";
				throw std::invalid_argument( message.str() );
			}
			target.push_back( q );
		}
	}

	bool SafetyFactorSystem::residual( std::vector<double> const &u, std::vector<double> &r )
	{
		std::size_t const m = spline.size();
		if ( u.size() != m )
			throw std::invalid_argument(
				"meq::SafetyFactorSystem::residual: one unknown per knot" );
		ToroidalFieldSpline trial( spline );
		bool finite = true;
		for ( double v : u )
			finite = finite && std::isfinite( v );
		if ( !finite )
		{
			++refusalCount;
			return false;
		}
		trial.setValues( u );
		// TOTAL: g^2 <= 0 is not an equilibrium, so it is refused without a
		// solve -- the step is rejected, not thrown.
		if ( !trial.positive() )
		{
			++refusalCount;
			return false;
		}
		stageChecksum( "toroidalSystem.in", u );
		++solveCount;
		SafetyFactorEvaluation evaluation;
		bool solved = false;
		try
		{
			solved = solve( trial, evaluation );
		}
		catch ( std::exception const & )
		{
			solved = false;
		}
		if ( !solved || evaluation.qOverG.size() != m )
		{
			++refusalCount;
			return false;
		}
		r.assign( m + 1, 0.0 );
		for ( std::size_t i = 0; i < m; ++i )
		{
			if ( !( evaluation.qOverG[ i ] > 0.0 ) || !std::isfinite( evaluation.qOverG[ i ] ) )
			{
				++refusalCount;
				return false;
			}
			r[ i ] = evaluation.qOverG[ i ]*std::sqrt( u[ i ] )*( target[ i ] > 0.0 ? 1.0 : -1.0 )
			         /target[ i ] - 1.0;
		}
		double measured = 0.0;
		if ( pin.kind == ScalePin::EdgeField )
		{
			double const edge = trial( 0.0 );
			measured = std::sqrt( edge );
			r[ m ] = ( edge - pin.value*pin.value )/( pin.value*pin.value );
		}
		else
		{
			measured = evaluation.measured;
			if ( !std::isfinite( measured ) )
			{
				++refusalCount;
				return false;
			}
			r[ m ] = ( measured - pin.value )/pin.value;
		}
		spline = trial;
		lastGood = u;
		qOverG = evaluation.qOverG;
		measuredValue = measured;
		stageChecksum( "toroidalSystem.out", r );
		return true;
	}

	bool SafetyFactorSystem::jacobian( std::vector<double> const &u, std::vector<double> &j )
	{
		std::size_t const m = spline.size();
		if ( u != lastGood )
			throw std::logic_error(
				"meq::SafetyFactorSystem::jacobian: asked at a point that is not the last "
				"admissible one; the derivative is formed at the state the solver holds" );
		SafetyFactorDerivative derivative;
		bool formed = false;
		try
		{
			formed = differentiate( spline, derivative );
		}
		catch ( std::exception const & )
		{
			formed = false;
		}
		if ( !formed || derivative.qOverG.size() != m*m
		     || ( pin.kind != ScalePin::EdgeField && derivative.measured.size() != m ) )
			return false;
		j.assign( ( m + 1 )*m, 0.0 );
		for ( std::size_t i = 0; i < m; ++i )
		{
			double const root = std::sqrt( u[ i ] );
			double const sign = target[ i ] > 0.0 ? 1.0 : -1.0;
			for ( std::size_t k = 0; k < m; ++k )
			{
				double value = root*derivative.qOverG[ i*m + k ];
				if ( i == k )
					value += qOverG[ i ]/( 2.0*root );
				j[ i*m + k ] = sign*value/target[ i ];
			}
		}
		for ( std::size_t k = 0; k < m; ++k )
			j[ m*m + k ] = pin.kind == ScalePin::EdgeField
			               ? spline.cardinal( k, 0.0 )/( pin.value*pin.value )
			               : derivative.measured[ k ]/pin.value;
		stageChecksum( "toroidalSystem.jacobian", j );
		return true;
	}

	std::vector<double> SafetyFactorSystem::floor() const
	{
		// See the header: 8 eps_solve on a q row, 4 eps_solve on a pin row
		// that reads the solve; R0B0's reads none.
		std::vector<double> out( spline.size() + 1, 8.0*accuracy );
		out.back() = pin.kind == ScalePin::EdgeField ? 0.0 : 4.0*accuracy;
		return out;
	}

	std::vector<double> SafetyFactorSystem::scaleGradient( std::vector<double> const &u )
	{
		std::vector<double> out( u.size() );
		for ( std::size_t j = 0; j < u.size(); ++j )
			out[ j ] = 1.0/u[ j ];
		return out;
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
			// vacuum one, g = const, so gg' = 0 -- and a knot there, value and
			// slope zero, is what makes SplineProfile's continuation carry that
			// rather than the innermost value and slope it happens to hold. Only added when the family does not
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

}
