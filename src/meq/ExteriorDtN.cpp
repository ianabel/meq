#include "ExteriorDtN.hpp"

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>

/*
 * HOW C_n IS EVALUATED, AND IT IS NOT THE FORMULA THE PLAN PRINTS.
 *
 * FREE-BOUNDARY-PLAN.md section 3 writes the angular functions as
 *
 *     C_n( mu ) = ( P_{n-2}( mu ) - P_n( mu ) ) / ( 2n - 1 ),
 *
 * which is correct and is a bad way to compute them. This file uses
 *
 *     C_n( mu ) = ( 1 - mu^2 ) P'_{n-1}( mu ) / ( n( n - 1 ) )
 *
 * instead. The two are the same function -- checked at 1.4e-42 in 40-digit
 * arithmetic over n = 2..10 and several mu -- by the standard Legendre identity
 *
 *     ( 2n - 1 )( 1 - mu^2 ) P'_{n-1}( mu ) = n( n - 1 )[ P_{n-2}( mu ) - P_n( mu ) ].
 *
 * TWO REASONS TO PREFER THE SECOND, AND THE FIRST IS ABOUT THE AXIS.
 *
 * The printed form is a DIFFERENCE OF TWO QUANTITIES THAT BOTH TEND TO 1 as
 * mu -> ±1, so it loses digits exactly where this basis is most delicate --
 * mu = ±1 IS the axis, where the semicircle's ends are. Measured, the ratio of
 * the larger operand to the result at n = 8:
 *
 *     mu       0.9      0.99     0.999    0.99999
 *     ratio    2.4      6.2      66       6666
 *
 * -- growing without bound, so about four digits are already gone at 1e-5 from
 * the axis and more beyond. The derivative form carries ( 1 - mu^2 ) as an
 * EXPLICIT FACTOR, so it loses nothing there and C_n( ±1 ) = 0 comes out
 * exactly zero rather than to round-off. Given that MEQ has just measured the
 * axis to cost one power of h in the trace conditioning (FB-A), a basis that
 * silently degrades there is not what this should be built on.
 *
 * AND HAVING THE FACTOR IS NOT ENOUGH -- HOW IT IS WRITTEN MATTERS TOO, by
 * another eight orders. See weightFactor() below: ( 1 - mu )( 1 + mu ) is not
 * the same computation as 1 - mu*mu near the axis, and only the first of them
 * keeps its digits there.
 *
 * The second reason is that it makes the ( 1 - mu^2 ) available for
 * cancellation. Every integral in this file carries the weight
 * dGamma/R = dmu/( 1 - mu^2 ), and dividing the derivative form by that weight
 * leaves P'_{n-1}( mu )/( n( n - 1 ) ) -- A POLYNOMIAL, with no 0/0 at the
 * endpoints and no singular quadrature anywhere. The printed form would have to
 * evaluate 0/0 at mu = ±1 and a near-cancellation beside them.
 *
 * THIS IS THE SAME FINDING Zernike.cpp RECORDS about the factorial sum for the
 * radial polynomials, from the other direction: there the textbook formula
 * cancels catastrophically for large degree, here for mu near the axis, and in
 * both cases an equivalent special-function form has no cancellation at all.
 * Both call Boost rather than carrying a recurrence, for the reason CLAUDE.md
 * gives under *Prefer a maintained library to a hand-rolled algorithm*.
 *
 * AND h_n FOLLOWS FROM IT WITHOUT BEING TRANSCRIBED. With the derivative form
 * the mass integral is the standard
 *
 *     integral_-1^1 ( 1 - mu^2 )[ P'_m( mu ) ]^2 dmu = 2 m( m + 1 )/( 2m + 1 ),
 *
 * at m = n - 1, divided by ( n( n - 1 ) )^2 -- which gives
 * 2/( n( n - 1 )( 2n - 1 ) ) exactly. So the closed form in the header is
 * DERIVED here rather than copied from the plan, and the two agreeing is a
 * check rather than a restatement. Confirmed numerically as well, at n = 2..8,
 * to every digit printed.
 */

namespace meq
{
	namespace
	{
		/// THE SMALLEST Gauss-Legendre rule coefficients() spends, and the one
		/// it spends up to fourteen modes.
		///
		/// **THE RULE SCALES WITH THE MODE COUNT, AS 2 lastMode() POINTS ABOVE
		/// THIS FLOOR.** A p-point rule is exact to degree 2p - 1 and the mass
		/// integrals have degree m + n - 2, so exactness of the basis alone
		/// needs only p >= lastMode(). A GENERAL trace needs the factor of two:
		/// its degree-d content aliases into the projection on mode n once
		/// d + n - 2 > 2p - 1, so at p = 2 lastMode() the first degree that
		/// can reach the last retained mode is 3 lastMode() + 1, a whole
		/// spectrum beyond it. A fixed thirty points reaches only m + n <= 61,
		/// and M-260 measured what that does past it: a copper trace's tail at
		/// 96 modes read a flat 0.82 out to N = 48 -- the quadrature, not the
		/// trace.
		///
		/// The floor keeps every count up to fourteen on Boost's tabulated
		/// thirty-point rule, summed in Boost's order, so those results are
		/// what they always were to the bit.
		constexpr int minimumGaussPoints = 30;

		using FloorRule = boost::math::quadrature::gauss<double, minimumGaussPoints>;

		int gaussPointsFor( int lastMode )
		{
			return std::max( minimumGaussPoints, 2*lastMode );
		}

		/// C_n( mu ) / ( 1 - mu^2 ), which is a POLYNOMIAL of degree n - 2 and
		/// is the shape every weighted integral in this file actually needs.
		///
		/// Having this rather than dividing C_n by the weight at the call site
		/// is the whole reason the singular weight never appears: the
		/// cancellation is done in closed form here, once.
		double basisOverWeight( int n, double mu )
		{
			return boost::math::legendre_p_prime( n - 1, mu )
			       /( static_cast<double>( n )*static_cast<double>( n - 1 ) );
		}

		/// 1 - mu^2, FACTORED AS ( 1 - mu )( 1 + mu ) AND NOT WRITTEN
		/// 1 - mu*mu, WHICH IS A REAL DIFFERENCE NEAR THE AXIS.
		///
		/// mu = ±1 is the axis, so this is evaluated at its worst argument
		/// exactly where MEQ is most delicate. Sterbenz's lemma makes ( 1 - mu )
		/// EXACT in floating point for mu in [ 1/2, 2 ], whereas mu*mu rounds
		/// first and the subtraction then cancels against a value that has
		/// already lost the low bits. Measured against a 60-digit reference at
		/// n = 8, absolute error in C_n:
		///
		///     1 - mu        1 - mu*mu      ( 1 - mu )( 1 + mu )
		///     7.5e-09       3.7e-09        5.0e-16
		///     9.1e-13       4.6e-13        5.2e-23
		///
		/// -- eight orders at the first of those, and the naive form has NO
		/// correct digits by 1e-16. Both forms still give exactly 0.0 AT the
		/// axis, so the difference is invisible to a test that only samples the
		/// endpoints; it is the approach that degrades.
		///
		/// Re-measured through THIS code against a long-double reference, at
		/// n = 8 and 1 - mu = 2^-7 .. 2^-53: the factored form is EXACT --
		/// 0.000e+00 absolute at every one of them -- while the naive form
		/// reaches 2.78e-17 absolute at 1 - mu = 7.45e-09, which against a value
		/// of 7.45e-09 is 3.7e-09 relative and matches the independent figure
		/// above. So this is not a transcribed measurement; it is the shipped
		/// function's.
		double weightFactor( double mu )
		{
			return ( 1.0 - mu )*( 1.0 + mu );
		}
	}

	ExteriorDtN::ExteriorDtN( double zCentreIn, double rhoGammaIn, int modesIn )
		: zCentreValue( zCentreIn ),
		  rhoGammaValue( rhoGammaIn ),
		  modeCountValue( modesIn ),
		  quadraturePointCount( gaussPointsFor( firstMode() + modesIn - 1 ) )
	{
		if ( !( rhoGammaIn > 0.0 ) )
			throw std::invalid_argument(
				"meq::ExteriorDtN: the radius of Gamma must be positive" );
		if ( modesIn < 1 )
			throw std::invalid_argument(
				"meq::ExteriorDtN: at least one mode is needed; modes are "
				"indexed by DEGREE from 2, so a count of k carries degrees "
				"2 .. k + 1" );

		/*
		 * THE RULE, FROM BOOST EITHER WAY. At the floor it is the tabulated
		 * one; above it the nodes are boost::math::legendre_p_zeros() -- the
		 * non-negative zeros of P_p, ascending, zero first when p is odd,
		 * which is the layout Boost's own tables use -- and the weights the
		 * standard 2/( ( 1 - x^2 ) P'_p( x )^2 ), with ( 1 - x )( 1 + x ) for
		 * weightFactor()'s reason.
		 */
		if ( quadraturePointCount == minimumGaussPoints )
		{
			ruleAbscissa.assign( FloorRule::abscissa().begin(),
			                     FloorRule::abscissa().end() );
			ruleWeights.assign( FloorRule::weights().begin(),
			                    FloorRule::weights().end() );
		}
		else
		{
			ruleAbscissa = boost::math::legendre_p_zeros<double>( quadraturePointCount );
			ruleWeights.resize( ruleAbscissa.size() );
			for ( std::size_t i = 0; i < ruleAbscissa.size(); ++i )
			{
				double const x = ruleAbscissa[ i ];
				double const slope = boost::math::legendre_p_prime( quadraturePointCount, x );
				ruleWeights[ i ] = 2.0/( weightFactor( x )*slope*slope );
			}
		}
	}

	int ExteriorDtN::firstMode()
	{
		return 2;
	}

	int ExteriorDtN::modeCount() const
	{
		return modeCountValue;
	}

	int ExteriorDtN::lastMode() const
	{
		return firstMode() + modeCountValue - 1;
	}

	double ExteriorDtN::zCentre() const
	{
		return zCentreValue;
	}

	double ExteriorDtN::rhoGamma() const
	{
		return rhoGammaValue;
	}

	int ExteriorDtN::quadraturePoints() const
	{
		return quadraturePointCount;
	}

	void ExteriorDtN::requireMode( int n ) const
	{
		if ( n < firstMode() || n > lastMode() )
			throw std::invalid_argument(
				"meq::ExteriorDtN: mode " + std::to_string( n ) + " is outside "
				"the range " + std::to_string( firstMode() ) + " .. "
				+ std::to_string( lastMode() ) + ". Modes are indexed by DEGREE "
				"and start at 2: degrees 0 and 1 are inadmissible, not merely "
				"omitted -- the mass divides by n( n - 1 ), and neither vanishes "
				"on the axis" );
	}

	void ExteriorDtN::direction( double radius, double z, double &mu, double &rho ) const
	{
		double const dz = z - zCentreValue;
		rho = std::hypot( radius, dz );

		// The centre itself has no direction. Reported as mu = 1 rather than
		// NaN because basis() is documented to depend on direction alone and a
		// caller sampling a grid should get a number; every admissible C_n
		// vanishes there anyway, so the value is the right one.
		mu = ( rho > 0.0 ) ? dz/rho : 1.0;
	}

	double ExteriorDtN::basis( int n, double radius, double z ) const
	{
		requireMode( n );
		double mu = 0.0;
		double rho = 0.0;
		direction( radius, z, mu, rho );

		// ( 1 - mu )( 1 + mu ) as an explicit factor, so the axis is exactly
		// zero and the approach to it does not lose digits.
		return weightFactor( mu )*basisOverWeight( n, mu );
	}

	double ExteriorDtN::basisDerivative( int n, double radius, double z ) const
	{
		requireMode( n );
		double mu = 0.0;
		double rho = 0.0;
		direction( radius, z, mu, rho );

		/*
		 * dC_n/dmu = -P_{n-1}( mu ), EXACTLY, which is one Legendre evaluation
		 * and no difference at all.
		 *
		 * Verified symbolically for n = 2..12. It follows from differentiating
		 * the printed form and applying the standard
		 * ( 2n - 1 ) P_{n-1} = P'_n - P'_{n-2}; but it is worth having as an
		 * identity in its own right, because the obvious implementation --
		 * differencing the two Legendre DERIVATIVES the way basis() would have
		 * to difference the two Legendre VALUES -- is both slower and
		 * needlessly indirect.
		 *
		 * Note the endpoints come out exactly: C_n'( +1 ) = -1 for every n, and
		 * C_n'( -1 ) = +1 for even n, -1 for odd. So the basis meets the axis at
		 * unit slope in mu whatever the mode, which is the derivative statement
		 * of the simple zero recorded above.
		 */
		return -boost::math::legendre_p( n - 1, mu );
	}

	double ExteriorDtN::symbol( int n ) const
	{
		requireMode( n );
		return ( 1.0 - static_cast<double>( n ) )/rhoGammaValue;
	}

	double ExteriorDtN::mass( int n ) const
	{
		requireMode( n );
		double const d = static_cast<double>( n );
		return 2.0/( d*( d - 1.0 )*( 2.0*d - 1.0 ) );
	}

	double ExteriorDtN::blockEntry( int n ) const
	{
		requireMode( n );
		return -symbol( n )*mass( n );
	}

	void ExteriorDtN::requireCoefficients( std::vector<double> const &a,
	                                       char const *method ) const
	{
		if ( static_cast<int>( a.size() ) != modeCountValue )
			throw std::invalid_argument(
				std::string( "meq::ExteriorDtN::" ) + method + ": expected "
				+ std::to_string( modeCountValue )
				+ " coefficients, one per mode, but was given "
				+ std::to_string( a.size() ) );
	}

	double ExteriorDtN::exterior( double radius, double z,
	                              std::vector<double> const &a ) const
	{
		requireCoefficients( a, "exterior" );

		double mu = 0.0;
		double rho = 0.0;
		direction( radius, z, mu, rho );

		/*
		 * REFUSED INSIDE Gamma, RATHER THAN EXTRAPOLATED.
		 *
		 * Every exterior mode is rho^( 1 - n ), which decays outward and GROWS
		 * WITHOUT BOUND inward -- at n = 20 a point at half the radius is
		 * 2^19 times the trace value. So an accidental interior call would not
		 * return a poor answer, it would return a huge one, and it would do so
		 * silently in a quantity a coupled residual is built from. The interior
		 * is what the mesh is for.
		 *
		 * The tolerance is relative and generous: a point ON Gamma is the
		 * ordinary case and must not be refused for round-off in a caller's own
		 * geometry.
		 */
		double const tolerance = 1.0e-12*rhoGammaValue;
		if ( rho < rhoGammaValue - tolerance )
			throw std::invalid_argument(
				"meq::ExteriorDtN::exterior: the point is INSIDE Gamma, where "
				"this expansion does not apply -- its modes decay outward and so "
				"grow without bound inward. The interior field is the mesh's" );

		double const scaled = rho/rhoGammaValue;
		double const weight = weightFactor( mu );

		double sum = 0.0;
		for ( int i = 0; i < modeCountValue; ++i )
		{
			int const n = firstMode() + i;
			double const radial = std::pow( scaled, 1.0 - static_cast<double>( n ) );
			sum += a[ static_cast<std::size_t>( i ) ]*radial*weight
			       *basisOverWeight( n, mu );
		}
		return sum;
	}

	std::vector<double> ExteriorDtN::coefficients(
		std::function<double( double, double )> const &trace ) const
	{
		/*
		 * a_n = ( 1/h_n ) integral_Gamma psi C_n dGamma/R, and the whole point
		 * of basisOverWeight() is that this integral has NO SINGULAR WEIGHT
		 * once C_n/( 1 - mu^2 ) is taken in closed form:
		 *
		 *     integral_-1^1 psi( mu ) * C_n( mu )/( 1 - mu^2 ) dmu.
		 *
		 * A plain Gauss-Legendre rule in mu is therefore the right instrument
		 * and not a compromise. On Gamma the point at parameter mu is
		 * R = rhoGamma sqrt( 1 - mu^2 ), z = zCentre + rhoGamma mu.
		 */
		std::vector<double> out( static_cast<std::size_t>( modeCountValue ), 0.0 );

		for ( int i = 0; i < modeCountValue; ++i )
		{
			int const n = firstMode() + i;
			auto integrand = [ & ]( double mu )
			{
				// The same factoring, for the same reason: this is R on Gamma and
				// the quadrature's outermost nodes sit near the axis.
				double const radius = rhoGammaValue
				                 *std::sqrt( std::max( 0.0, weightFactor( mu ) ) );
				double const z = zCentreValue + rhoGammaValue*mu;
				return trace( radius, z )*basisOverWeight( n, mu );
			};

			// boost::math::quadrature::gauss::integrate()'s own sum: the
			// centre once when the count is odd, then each symmetric pair.
			std::size_t start = 0;
			double sum = 0.0;
			if ( quadraturePointCount % 2 == 1 )
			{
				sum = integrand( 0.0 )*ruleWeights[ 0 ];
				start = 1;
			}
			for ( std::size_t q = start; q < ruleAbscissa.size(); ++q )
			{
				double const plus = integrand( ruleAbscissa[ q ] );
				double const minus = integrand( -ruleAbscissa[ q ] );
				sum += ( plus + minus )*ruleWeights[ q ];
			}
			out[ static_cast<std::size_t>( i ) ] = sum/mass( n );
		}
		return out;
	}

	/*
	 * THE DECAY DIAGNOSTIC, AND WHY A RAW COEFFICIENT IS THE WRONG THING TO
	 * PRINT.
	 *
	 * The modes are orthogonal in dGamma/R but not orthoNORMAL: mode n has norm
	 * sqrt( mass( n ) ) rather than 1, so the trace decomposes as
	 *
	 *     || sum a_n C_n ||^2 = sum a_n^2 mass( n )
	 *
	 * and what mode n actually contributes is | a_n | sqrt( mass( n ) ). With
	 * mass( n ) = 2/( n( n - 1 )( 2n - 1 ) ) ~ 1/n^3, that conversion factor is
	 * ~ n^( -3/2 ) -- a strong decay in its own right, and one a reader of the
	 * raw numbers has no way to see.
	 *
	 * SO THE RAW VIEW IS WRONG IN BOTH DIRECTIONS, AND IT IS REASSURING IN THE
	 * ONE THAT MATTERS. Raw coefficients that look FLAT -- the picture of a
	 * truncation going nowhere -- are an energy spectrum already falling like
	 * n^( -3/2 ); raw coefficients that look like they are GROWING at n^( 3/2 )
	 * are the flat-in-energy case, which is the genuinely unconverged one. A
	 * caller reading raw a_n therefore worries about the healthy case and is
	 * comforted by the sick one.
	 *
	 * The measured stake, from FREE-BOUNDARY-PLAN.md section 11.6: on section
	 * 7.16's machine case N = 6 costs 0.3% in psi_ax and N = 10 is converged.
	 * That is a small enough margin that it has to be read off an honest
	 * quantity rather than eyeballed.
	 */
	std::vector<double> ExteriorDtN::modeAmplitudes(
		std::vector<double> const &a ) const
	{
		requireCoefficients( a, "modeAmplitudes" );

		std::vector<double> out( static_cast<std::size_t>( modeCountValue ), 0.0 );
		for ( int i = 0; i < modeCountValue; ++i )
		{
			std::size_t const k = static_cast<std::size_t>( i );
			out[ k ] = std::fabs( a[ k ] )*std::sqrt( mass( firstMode() + i ) );
		}
		return out;
	}

	double ExteriorDtN::traceNorm( std::vector<double> const &a ) const
	{
		requireCoefficients( a, "traceNorm" );

		// The root-sum-square of the amplitudes, which by the orthogonality
		// above IS the L2( dGamma/R ) norm of the trace. Accumulated from the
		// amplitudes rather than from a_n^2 mass( n ) directly so that the two
		// cannot drift apart: this is the same arithmetic modeAmplitudes()
		// does, and it should stay that way.
		std::vector<double> const amplitudes = modeAmplitudes( a );

		double sum = 0.0;
		for ( double value : amplitudes )
			sum += value*value;
		return std::sqrt( sum );
	}

	/*
	 * THE TAIL IS THE LAST TWO MODES AND NOT THE LAST ONE, AND THAT IS A
	 * PARITY FACT RATHER THAN A HEDGE.
	 *
	 * C_n( -mu ) = ( -1 )^n C_n( mu ), so an UP-DOWN SYMMETRIC trace -- which
	 * is the ordinary tokamak case, and every current loop centred on
	 * z = zCentre -- has identically zero coefficients in every odd mode. Half
	 * the spectrum is exactly 0 by symmetry rather than by convergence.
	 *
	 * A summary reading amplitudes.back() alone therefore returns EXACTLY ZERO
	 * whenever the last retained degree happens to be the absent parity, which
	 * is half of all mode counts -- and zero reads as "perfectly converged".
	 * Measured on a unit loop at four radii with twelve modes (degrees 2..13,
	 * so the last is the odd 13): 0.000e+00 at every one of them, while the raw
	 * ratio | a_12/a_2 | moved over six orders. The diagnostic would have been
	 * silently useless on exactly the configurations MEQ is for.
	 *
	 * Two consecutive degrees always span both parities, so the largest of the
	 * last two is the smallest robust answer and needs no symmetry argument
	 * from the caller.
	 */
	double ExteriorDtN::truncationRatio( std::vector<double> const &a ) const
	{
		requireCoefficients( a, "truncationRatio" );

		std::vector<double> const amplitudes = modeAmplitudes( a );

		double largest = 0.0;
		for ( double value : amplitudes )
			largest = std::max( largest, value );

		// A zero trace has no unresolved tail, and it is the state a coupled
		// Newton starts from -- so this must be a number a driver can print
		// rather than a NaN. See the header.
		if ( !( largest > 0.0 ) )
			return 0.0;

		double tail = amplitudes.back();
		if ( amplitudes.size() > 1 )
			tail = std::max( tail, amplitudes[ amplitudes.size() - 2 ] );

		return tail/largest;
	}

	/*
	 * THE EXTERIOR'S SUPPORT, CONNECTED TO Gamma. See the header for what it is
	 * and why only the connected part; this is the fill and the quadrature.
	 *
	 * u = ln( rho/rhoGamma ) and theta are both cut at the same step
	 * delta = pi/angularCells, so a cell is square in the metric
	 * d( ln rho )^2 + dtheta^2 -- the conformal one, in which a support of a
	 * given shape costs the same number of cells at every radius. Ring 0 is the
	 * one whose inner edge is Gamma.
	 *
	 * EVERY CELL IS EVALUATED AT MOST ONCE, the visited mark doubling as the
	 * cache, and the fill is breadth-first from ring 0 in ( ring, angle )
	 * order, so the sum is taken in an order fixed by the geometry alone.
	 */
	ExteriorSupportIntegral exteriorSupportIntegral(
		ExteriorDtN const &exterior,
		std::function<double( double, double )> const &psi,
		std::function<double( double, double, double )> const &source,
		int angularCells, double radiusRatio )
	{
		if ( angularCells < 4 )
			throw std::invalid_argument(
				"meq::exteriorSupportIntegral: at least four angular cells are "
				"needed to resolve anything" );
		if ( !( radiusRatio > 1.0 ) || !std::isfinite( radiusRatio ) )
			throw std::invalid_argument(
				"meq::exteriorSupportIntegral: the fill's outer limit must be a "
				"finite radius beyond Gamma, rho/rhoGamma > 1" );

		double const pi = std::acos( -1.0 );
		double const delta = pi/static_cast<double>( angularCells );
		int const rings = static_cast<int>( std::ceil( std::log( radiusRatio )/delta ) );
		double const rhoGamma = exterior.rhoGamma();
		double const zCentre = exterior.zCentre();

		// Three-point Gauss-Legendre on [ 0, 1 ].
		double const offset = 0.5*std::sqrt( 0.6 );
		double const node[ 3 ] = { 0.5 - offset, 0.5, 0.5 + offset };
		double const weight[ 3 ] = { 5.0/18.0, 8.0/18.0, 5.0/18.0 };

		struct Cell
		{
			bool inSupport = false;
			bool finite = true;
			double sum = 0.0;
		};

		auto evaluate = [ & ]( int ring, int angle )
		{
			Cell cell;
			for ( int a = 0; a < 3; ++a )
			{
				double const rho = rhoGamma*std::exp( ( ring + node[ a ] )*delta );
				for ( int b = 0; b < 3; ++b )
				{
					double const theta = ( angle + node[ b ] )*delta;
					double const radius = rho*std::sin( theta );
					double const z = zCentre + rho*std::cos( theta );
					double const value = source( radius, z, psi( radius, z ) );
					if ( value == 0.0 )
						continue;
					cell.inSupport = true;
					if ( !std::isfinite( value ) )
					{
						cell.finite = false;
						continue;
					}
					// | F |/R dR dz = | F | rho^2/R du dtheta.
					cell.sum += weight[ a ]*weight[ b ]*delta*delta
					            *std::abs( value )*rho*rho/radius;
				}
			}
			return cell;
		};

		std::size_t const width = static_cast<std::size_t>( angularCells );
		std::vector<char> visited( static_cast<std::size_t>( rings )*width, 0 );
		auto index = [ & ]( int ring, int angle )
		{
			return static_cast<std::size_t>( ring )*width
			       + static_cast<std::size_t>( angle );
		};

		ExteriorSupportIntegral result;
		std::deque<std::pair<int, int>> queue;
		auto visit = [ & ]( int ring, int angle )
		{
			if ( angle < 0 || angle >= angularCells || ring < 0 )
				return;
			if ( ring >= rings )
			{
				result.bounded = false;
				return;
			}
			char &mark = visited[ index( ring, angle ) ];
			if ( mark )
				return;
			mark = 1;
			Cell const cell = evaluate( ring, angle );
			if ( !cell.inSupport )
				return;
			++result.cells;
			if ( ring == 0 )
				++result.touching;
			result.integral += cell.sum;
			if ( !cell.finite )
				result.finite = false;
			result.reach = std::max( result.reach,
			                         rhoGamma*std::exp( ( ring + 1 )*delta ) );
			queue.emplace_back( ring, angle );
		};

		for ( int angle = 0; angle < angularCells; ++angle )
			visit( 0, angle );
		while ( !queue.empty() )
		{
			auto const [ ring, angle ] = queue.front();
			queue.pop_front();
			visit( ring + 1, angle );
			visit( ring - 1, angle );
			visit( ring, angle + 1 );
			visit( ring, angle - 1 );
		}

		if ( !result.finite )
			result.integral = std::numeric_limits<double>::infinity();
		return result;
	}

}
