#include "ExteriorDtN.hpp"

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <algorithm>
#include <cmath>
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
 * dGamma/r = dmu/( 1 - mu^2 ), and dividing the derivative form by that weight
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
		/// The Gauss-Legendre rule coefficients() spends.
		///
		/// FIXED AT 30 POINTS RATHER THAN SCALED WITH THE MODE COUNT, because
		/// Boost's rule is selected by a template parameter and a runtime choice
		/// would mean a switch over a handful of instantiations for no gain: 30
		/// points is exact for polynomials up to degree 59, and the mass
		/// integrals have degree m + n - 2, so it is exact for every pair of
		/// modes with m + n <= 61. Nobody will use forty modes -- the spectrum
		/// falls geometrically at a rate set by rho_plasma/rho_Gamma, which is
		/// why the plan says N is small.
		///
		/// For a general trace it is a quadrature and not an identity, which is
		/// what quadraturePoints() exists to let a caller see.
		constexpr int gaussPoints = 30;

		using GaussRule = boost::math::quadrature::gauss<double, gaussPoints>;

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
		  quadraturePointCount( gaussPoints )
	{
		if ( !( rhoGammaIn > 0.0 ) )
			throw std::invalid_argument(
				"meq::ExteriorDtN: the radius of Gamma must be positive" );
		if ( modesIn < 1 )
			throw std::invalid_argument(
				"meq::ExteriorDtN: at least one mode is needed; modes are "
				"indexed by DEGREE from 2, so a count of k carries degrees "
				"2 .. k + 1" );
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

	void ExteriorDtN::direction( double r, double z, double &mu, double &rho ) const
	{
		double const dz = z - zCentreValue;
		rho = std::hypot( r, dz );

		// The centre itself has no direction. Reported as mu = 1 rather than
		// NaN because basis() is documented to depend on direction alone and a
		// caller sampling a grid should get a number; every admissible C_n
		// vanishes there anyway, so the value is the right one.
		mu = ( rho > 0.0 ) ? dz/rho : 1.0;
	}

	double ExteriorDtN::basis( int n, double r, double z ) const
	{
		requireMode( n );
		double mu = 0.0;
		double rho = 0.0;
		direction( r, z, mu, rho );

		// ( 1 - mu )( 1 + mu ) as an explicit factor, so the axis is exactly
		// zero and the approach to it does not lose digits.
		return weightFactor( mu )*basisOverWeight( n, mu );
	}

	double ExteriorDtN::basisDerivative( int n, double r, double z ) const
	{
		requireMode( n );
		double mu = 0.0;
		double rho = 0.0;
		direction( r, z, mu, rho );

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

	double ExteriorDtN::exterior( double r, double z,
	                              std::vector<double> const &a ) const
	{
		requireCoefficients( a, "exterior" );

		double mu = 0.0;
		double rho = 0.0;
		direction( r, z, mu, rho );

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
		 * a_n = ( 1/h_n ) integral_Gamma psi C_n dGamma/r, and the whole point
		 * of basisOverWeight() is that this integral has NO SINGULAR WEIGHT
		 * once C_n/( 1 - mu^2 ) is taken in closed form:
		 *
		 *     integral_-1^1 psi( mu ) * C_n( mu )/( 1 - mu^2 ) dmu.
		 *
		 * A plain Gauss-Legendre rule in mu is therefore the right instrument
		 * and not a compromise. On Gamma the point at parameter mu is
		 * r = rhoGamma sqrt( 1 - mu^2 ), z = zCentre + rhoGamma mu.
		 */
		std::vector<double> out( static_cast<std::size_t>( modeCountValue ), 0.0 );

		for ( int i = 0; i < modeCountValue; ++i )
		{
			int const n = firstMode() + i;
			auto integrand = [ & ]( double mu )
			{
				// The same factoring, for the same reason: this is r on Gamma and
				// the quadrature's outermost nodes sit near the axis.
				double const r = rhoGammaValue
				                 *std::sqrt( std::max( 0.0, weightFactor( mu ) ) );
				double const z = zCentreValue + rhoGammaValue*mu;
				return trace( r, z )*basisOverWeight( n, mu );
			};

			out[ static_cast<std::size_t>( i ) ] =
				GaussRule::integrate( integrand, -1.0, 1.0 )/mass( n );
		}
		return out;
	}

	/*
	 * THE DECAY DIAGNOSTIC, AND WHY A RAW COEFFICIENT IS THE WRONG THING TO
	 * PRINT.
	 *
	 * The modes are orthogonal in dGamma/r but not orthoNORMAL: mode n has norm
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
		// above IS the L2( dGamma/r ) norm of the trace. Accumulated from the
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

}
