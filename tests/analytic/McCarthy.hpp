#ifndef MEQ_TESTS_MCCARTHY_HPP
#define MEQ_TESTS_MCCARTHY_HPP

/*
 * McCarthy equilibria: exact solutions of the Grad-Shafranov equation for
 * "dissimilar source functions", meaning a pressure and a toroidal field
 * function that are not proportional to one another.
 *
 * These sit between the two other fixtures in this directory and are the
 * reason to keep all three:
 *
 *   Soloviev.hpp               F constant in psi   -- linear, dF/dpsi = 0
 *   McCarthy.hpp               F linear in psi     -- dF/dpsi = T, constant
 *   ManufacturedNonlinear.hpp  F nonlinear in psi  -- psi^2 and exp(-psi)
 *
 * The middle case is the useful one for a Newton solve: the source genuinely
 * depends on psi, so the Jacobian picks up a real mass term, but dF/dpsi is a
 * constant that can be checked by inspection. A Jacobian bug that survives
 * Soloviev (where dF/dpsi vanishes) and is masked by the algebra of the
 * nonlinear case will show up here plainly.
 *
 * Sources:
 *   refs/Refs.md                doi 10.1063/1.873630 -- McCarthy, Phys.
 *                               Plasmas 6 (1999) 3554
 *   refs/HDG-GradShafranov.pdf  eqs (13)-(14) and Example 4, the ASDEX
 *                               Upgrade configuration whose coefficients
 *                               asdex() reproduces
 */

#include <array>
#include <cmath>
#include <stdexcept>

namespace meq
{
namespace analytic
{

/*
 * THE PARAMETRISATION, which is where the traps are.
 *
 * McCarthy takes
 *
 *     p(psi)   = (S/mu0) psi,        g(psi)^2 = T psi^2 + 2 U psi + g0^2,
 *
 * so that, with F := mu0 r^2 dp/dpsi + g dg/dpsi as everywhere else in MEQ,
 *
 *     mu0 r^2 dp/dpsi = S r^2,
 *     g dg/dpsi       = (1/2) d(g^2)/dpsi = T psi + U,
 *     ------------------------------------------------
 *     F(r, z, psi)    = T psi + S r^2 + U.                      (HDG-GS-1 (13))
 *
 * Linear in psi, so dF/dpsi = T identically.
 *
 * The eighteen coefficients below index the eigenfunction expansion of
 * HDG-GS-1 eq (14). The first two are not free: they carry the homogeneous
 * part, and HDG-GS-1 fixes them as
 *
 *     U = -c[0] T,      S = -c[1] T,
 *
 * (their c_1 and c_2, one-based) which gives
 *
 *     F = T ( psi - c[0] - c[1] r^2 ).
 *
 * A CHECK ON THE OLD CODE, which had this wrong. The pre-modernisation
 * fixture returned pPrime() = -c[0] T and ffPrime() = T psi - c[1]. Both are
 * wrong against its own source term: c[0] and c[1] are swapped, and pPrime was
 * additionally missing its 1/mu0. Its source function was right, which is
 * presumably why nothing ever noticed -- the solver only ever called that one.
 * The corrected forms are below, and mu0 is an explicit constructor argument
 * rather than an assumption.
 *
 * NOTE ON UNITS: f() returns F, NOT F/r. The 1/r belongs to the weak form and
 * GradShafranovSolver applies it. The old fixture returned F/r, matching a
 * different convention. Soloviev.hpp returns F too.
 */

/// An exact equilibrium with a source linear in psi, from the eigenfunction
/// families of McCarthy (1999). Lengths are normalised to the major radius;
/// r must be strictly positive, since the expansion contains Y_1(r).
class McCarthyEquilibrium
{
	public:
		/// @param tIn   the parameter T: F = T psi + S r^2 + U, so dF/dpsi = T.
		/// @param cIn   the eighteen expansion coefficients of HDG-GS-1 eq (14).
		///              The first two also fix U and S, see the note above.
		/// @param mu0In vacuum permeability. Defaults to 1, i.e. normalised
		///              units, which is what the benchmarks use.
		McCarthyEquilibrium( double tIn, std::array<double, 18> const &cIn,
		                     double mu0In = 1.0 )
			: t( tIn ), c( cIn ), mu0( mu0In )
		{
			if ( !( mu0 > 0.0 ) )
			{
				throw std::invalid_argument( "mu0 must be positive" );
			}
		}

		/// The ASDEX Upgrade configuration of refs/HDG-GradShafranov.pdf
		/// Example 4 — an upwards-oriented X-point.
		static McCarthyEquilibrium asdex()
		{
			return McCarthyEquilibrium( 17.8116, {
				 0.17795,   -0.03291,   1.4934,    -0.4818,
				-1.1759,    -0.162,     0.3722,     0.07697,
				 1.2959,     0.5881,    1.5820,    -0.009059,
				 2.2388,     0.4186,    1.195,     -0.4265,
				 0.8057,    -0.004804
			} );
		}

		/// The poloidal flux function, HDG-GS-1 eq (14).
		double psi( double r, double z ) const
		{
			double const p = std::sqrt( t );
			double const q = 0.5 * p;
			double const nu = std::sqrt( 0.75 ) * p;
			double const s = std::sqrt( r * r + z * z );

			return c[ 0 ]
			     + c[ 1 ] * r * r
			     + r * std::cyl_bessel_j( 1, p * r ) * ( c[ 2 ] + c[ 3 ] * z )
			     + c[ 4 ] * std::cos( p * z ) + c[ 5 ] * std::sin( p * z )
			     + r * r * ( c[ 6 ] * std::cos( p * z ) + c[ 7 ] * std::sin( p * z ) )
			     + c[ 8 ] * std::cos( p * s ) + c[ 9 ] * std::sin( p * s )
			     + r * std::cyl_bessel_j( 1, nu * r )
			           * ( c[ 10 ] * std::cos( q * z ) + c[ 11 ] * std::sin( q * z ) )
			     + r * std::cyl_bessel_j( 1, q * r )
			           * ( c[ 12 ] * std::cos( nu * z ) + c[ 13 ] * std::sin( nu * z ) )
			     + r * std::cyl_neumann( 1, nu * r )
			           * ( c[ 14 ] * std::cos( q * z ) + c[ 15 ] * std::sin( q * z ) )
			     + r * std::cyl_neumann( 1, q * r )
			           * ( c[ 16 ] * std::cos( nu * z ) + c[ 17 ] * std::sin( nu * z ) );
		}

		/// grad_bar(psi) = ( d_r psi, d_z psi ). Not the HDG flux; see flux().
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			double const p = std::sqrt( t );
			double const q = 0.5 * p;
			double const nu = std::sqrt( 0.75 ) * p;
			double const s = std::sqrt( r * r + z * z );

			// d( r C_1( a r ) )/dr for C = J and C = Y, using
			// C_1'(x) = ( C_0(x) - C_2(x) ) / 2.
			auto dRBesselJ = [ & ]( double a )
			{
				return std::cyl_bessel_j( 1, a * r )
				     + 0.5 * a * r * ( std::cyl_bessel_j( 0, a * r )
				                     - std::cyl_bessel_j( 2, a * r ) );
			};
			auto dRNeumann = [ & ]( double a )
			{
				return std::cyl_neumann( 1, a * r )
				     + 0.5 * a * r * ( std::cyl_neumann( 0, a * r )
				                     - std::cyl_neumann( 2, a * r ) );
			};

			dPsiDr = 2.0 * c[ 1 ] * r
			       + dRBesselJ( p ) * ( c[ 2 ] + c[ 3 ] * z )
			       + 2.0 * r * ( c[ 6 ] * std::cos( p * z ) + c[ 7 ] * std::sin( p * z ) )
			       - ( c[ 8 ] * p * r / s ) * std::sin( p * s )
			       + ( c[ 9 ] * p * r / s ) * std::cos( p * s )
			       + dRBesselJ( nu ) * ( c[ 10 ] * std::cos( q * z )
			                           + c[ 11 ] * std::sin( q * z ) )
			       + dRBesselJ( q ) * ( c[ 12 ] * std::cos( nu * z )
			                          + c[ 13 ] * std::sin( nu * z ) )
			       + dRNeumann( nu ) * ( c[ 14 ] * std::cos( q * z )
			                           + c[ 15 ] * std::sin( q * z ) )
			       + dRNeumann( q ) * ( c[ 16 ] * std::cos( nu * z )
			                          + c[ 17 ] * std::sin( nu * z ) );

			dPsiDz = r * std::cyl_bessel_j( 1, p * r ) * c[ 3 ]
			       + p * ( -c[ 4 ] * std::sin( p * z ) + c[ 5 ] * std::cos( p * z ) )
			       + r * r * p * ( -c[ 6 ] * std::sin( p * z )
			                     + c[ 7 ] * std::cos( p * z ) )
			       - ( c[ 8 ] * p * z / s ) * std::sin( p * s )
			       + ( c[ 9 ] * p * z / s ) * std::cos( p * s )
			       + r * std::cyl_bessel_j( 1, nu * r )
			           * q * ( -c[ 10 ] * std::sin( q * z ) + c[ 11 ] * std::cos( q * z ) )
			       + r * std::cyl_bessel_j( 1, q * r )
			           * nu * ( -c[ 12 ] * std::sin( nu * z ) + c[ 13 ] * std::cos( nu * z ) )
			       + r * std::cyl_neumann( 1, nu * r )
			           * q * ( -c[ 14 ] * std::sin( q * z ) + c[ 15 ] * std::cos( q * z ) )
			       + r * std::cyl_neumann( 1, q * r )
			           * nu * ( -c[ 16 ] * std::sin( nu * z ) + c[ 17 ] * std::cos( nu * z ) );
		}

		/// The HDG flux q = grad_bar(psi) / r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// F = T psi + S r^2 + U = T ( psi - c[0] - c[1] r^2 ).
		/// Returns F, not F/r.
		double f( double r, double /*z*/, double psiValue ) const
		{
			return t * ( psiValue - c[ 0 ] - c[ 1 ] * r * r );
		}

		/// dF/dpsi = T, constant. This is the property that makes the fixture
		/// worth having: a Newton Jacobian either has this mass term or it
		/// does not, and no algebra hides the difference.
		double dFdPsi( double, double, double ) const
		{
			return t;
		}

		/// dp/dpsi = S/mu0, with S = -c[1] T. Constant in psi.
		double pPrime( double ) const
		{
			return -c[ 1 ] * t / mu0;
		}

		/// g dg/dpsi = T psi + U, with U = -c[0] T.
		double ggPrime( double psiValue ) const
		{
			return t * psiValue - c[ 0 ] * t;
		}

		/// Delta*(psi), by central differences of psi(). The Grad-Shafranov
		/// equation is -Delta*(psi) = F, so this must equal -f(); the test
		/// suite asserts exactly that, which checks the eighteen-term
		/// transcription rather than trusting it.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) ) / ( 2.0 * h ) / rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) ) / ( 2.0 * h );
			double const dZZ = ( psi( r, z + h ) - 2.0 * psi( r, z ) + psi( r, z - h ) )
			                 / ( h * h );

			return r * dRInner + dZZ;
		}

		double getT() const
		{
			return t;
		}

	private:
		double t;
		std::array<double, 18> c;
		double mu0;
};

}
}

#endif // MEQ_TESTS_MCCARTHY_HPP
