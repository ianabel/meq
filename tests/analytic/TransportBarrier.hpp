#ifndef MEQ_TESTS_TRANSPORTBARRIER_HPP
#define MEQ_TESTS_TRANSPORTBARRIER_HPP

/*
 * An internal transport barrier, refs/HDG-GradShafranov-Adaptive.pdf section
 * 4.3, eq (25):
 *
 *     p( psi ) = [ ( 1 + H erf( s ( psi - psi_0 ) ) )/( 1 + H ) ]
 *                ( 1 - ( 1 - psi )^a )^b
 *
 * with H = 0.5, s = 40, psi_0 = 0.3, a = 4, b = 2 -- the steepest of the four
 * profiles drawn in that paper's Figure 10. Transcribed from the RENDERED page
 * 17, not from pdftotext; see CLAUDE.md on why that distinction is load bearing
 * for this pair of papers.
 *
 * H and s control the height and steepness of the barrier, psi_0 its radial
 * location, and a and b are natural numbers shaping the underlying pressure.
 * The paper notes that psi is taken normalised to [0,1] here; the profile is
 * still perfectly well defined outside that interval, but it grows like
 * psi^( a b ) = psi^8 for large psi, which is worth knowing when a Newton step
 * overshoots.
 *
 * THE CONVENTION MAP. The paper takes g = constant for this example -- the
 * caption of its Figure 10 says so in as many words -- so g dg/dpsi = 0 and
 *
 *     F( r, psi ) = mu0 r^2 dp/dpsi,     mu0 = 1,
 *
 * which is MEQ's F exactly. Unlike section 4.2, no closed form for F is printed:
 * only p is given, and the source has to be obtained by differentiating it. The
 * check that it was done right is Figure 10's centre panel, which plots f( psi )
 * -- the source without its r^2 -- peaking a little above 10 just below
 * psi = 0.4. pPrime( 0.3 ) here is 10.08.
 *
 * DIFFERENTIATED, NOT FINITE-DIFFERENCED. erf appears, and its derivative
 * ( 2/sqrt( pi ) ) exp( -x^2 ) is elementary, so there is no excuse for a
 * difference quotient: MEQ closes the semi-linear problem by Newton, and
 * CLAUDE.md records that a Jacobian wrong by 5 per cent leaves every error and
 * every convergence rate unchanged to six significant figures while quietly
 * dropping Newton to linear convergence. The test suite compares both
 * derivatives against central differences of the level below.
 *
 * NO EXACT SOLUTION, so there is no psi(), gradPsi(), flux() or deltaStarFD()
 * here -- see the header comment of PressurePedestal.hpp, which also records
 * the more awkward property this profile shares with those three: p'( 0 ) = 0,
 * so F( r, 0 ) = 0 and psi == 0 solves the homogeneous Dirichlet problem. That
 * p'( 0 ) vanishes is not obvious from eq (25) and is worth spelling out: with
 * b = 2 the factor ( 1 - ( 1 - psi )^a )^b is O( psi^2 ) at the origin, so both
 * it and its first derivative vanish there.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/// The transport-barrier pressure profile of eq (25) and the Grad-Shafranov
/// source it generates.
class TransportBarrier
{
	public:
		/// @param hIn      barrier height. 0.5 in the paper.
		/// @param sIn      barrier steepness. 40 in the paper.
		/// @param psi0In   barrier location in psi. 0.3 in the paper.
		/// @param aIn      exponent inside the bracket, a natural number. 4.
		/// @param bIn      exponent of the bracket, a natural number. 2.
		TransportBarrier( double hIn, double sIn, double psi0In, int aIn, int bIn )
			: hValue( hIn ), sValue( sIn ), psi0Value( psi0In ),
			  aValue( aIn ), bValue( bIn )
		{
		}

		/// Section 4.3 as published: the steepest barrier of its Figure 10.
		static TransportBarrier barrier()
		{
			return TransportBarrier( 0.5, 40.0, 0.3, 4, 2 );
		}

		/// The pressure, eq (25).
		double p( double psi ) const
		{
			return errorFactor( psi )*shape( psi );
		}

		/// dp/dpsi = A' B + A B', with A the erf factor and B the bracket.
		double pPrime( double psi ) const
		{
			return errorFactorPrime( psi )*shape( psi )
			       + errorFactor( psi )*shapePrime( psi );
		}

		/// d2p/dpsi2 = A'' B + 2 A' B' + A B''.
		double pDoublePrime( double psi ) const
		{
			return errorFactorDoublePrime( psi )*shape( psi )
			       + 2.0*errorFactorPrime( psi )*shapePrime( psi )
			       + errorFactor( psi )*shapeDoublePrime( psi );
		}

		/// F = mu0 r^2 p'( psi ) with mu0 = 1 and g constant. Returns F, not
		/// F/r, as meq::Source::f() is documented to.
		double f( double r, double /*z*/, double psi ) const
		{
			return r*r*pPrime( psi );
		}

		/// dF/dpsi = r^2 p''( psi ).
		double dFdPsi( double r, double /*z*/, double psi ) const
		{
			return r*r*pDoublePrime( psi );
		}

	private:
		/// A( psi ) = ( 1 + H erf( s ( psi - psi_0 ) ) )/( 1 + H ).
		double errorFactor( double psi ) const
		{
			return ( 1.0 + hValue*std::erf( sValue*( psi - psi0Value ) ) )
			       /( 1.0 + hValue );
		}

		/// A'( psi ). d/dx erf( x ) = ( 2/sqrt( pi ) ) exp( -x^2 ).
		double errorFactorPrime( double psi ) const
		{
			double const arg = sValue*( psi - psi0Value );
			return hValue*sValue*2.0*std::exp( -arg*arg )
			       /( std::sqrt( M_PI )*( 1.0 + hValue ) );
		}

		/// A''( psi ) = A'( psi ) * ( -2 s^2 ( psi - psi_0 ) ).
		double errorFactorDoublePrime( double psi ) const
		{
			return errorFactorPrime( psi )
			       *( -2.0*sValue*sValue*( psi - psi0Value ) );
		}

		/// w( psi ) = 1 - ( 1 - psi )^a, the bracket before it is raised to b.
		double bracket( double psi ) const
		{
			return 1.0 - std::pow( 1.0 - psi, static_cast<double>( aValue ) );
		}

		/// B( psi ) = w^b.
		double shape( double psi ) const
		{
			return std::pow( bracket( psi ), static_cast<double>( bValue ) );
		}

		/// B' = b w^( b-1 ) w', with w' = a ( 1 - psi )^( a-1 ).
		double shapePrime( double psi ) const
		{
			double const w = bracket( psi );
			double const wPrime = aValue*std::pow( 1.0 - psi,
			                                       static_cast<double>( aValue - 1 ) );
			return bValue*std::pow( w, static_cast<double>( bValue - 1 ) )*wPrime;
		}

		/// B'' = b ( b-1 ) w^( b-2 ) w'^2 + b w^( b-1 ) w'', with
		/// w'' = -a ( a-1 ) ( 1 - psi )^( a-2 ).
		double shapeDoublePrime( double psi ) const
		{
			double const w = bracket( psi );
			double const wPrime = aValue*std::pow( 1.0 - psi,
			                                       static_cast<double>( aValue - 1 ) );
			double const wDoublePrime = -aValue*( aValue - 1 )
			                            *std::pow( 1.0 - psi,
			                                       static_cast<double>( aValue - 2 ) );
			return bValue*( bValue - 1 )
			       *std::pow( w, static_cast<double>( bValue - 2 ) )*wPrime*wPrime
			       + bValue*std::pow( w, static_cast<double>( bValue - 1 ) )
			         *wDoublePrime;
		}

		double hValue, sValue, psi0Value;
		int aValue, bValue;
};

}
}

#endif // MEQ_TESTS_TRANSPORTBARRIER_HPP
