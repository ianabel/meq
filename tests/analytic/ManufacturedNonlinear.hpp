#ifndef MEQ_TESTS_MANUFACTUREDNONLINEAR_HPP
#define MEQ_TESTS_MANUFACTUREDNONLINEAR_HPP

/*
 * A manufactured solution whose source depends genuinely on psi: Example 5 of
 * refs/HDG-GradShafranov.pdf (Sanchez-Vizuet & Solano, CPC 235 (2019) 120-132).
 *
 * The Solov'ev equilibria in Soloviev.hpp are exact, but their source is
 * independent of psi, so a Newton solve converges on them in a single step and
 * a wrong dF/dpsi cannot be detected. This is the case that exercises Newton:
 * the paper's own words are that "analytic solutions for pressure and current
 * profiles that result in non linear source terms as a function of psi are not
 * available", so it manufactures one, "including a linear, a quadratic and an
 * exponential term as functions of psi in the source".
 *
 * The imposed flux is
 *
 *     psi( r, z ) = sin( kr ( r + r0 ) ) cos( kz z ),
 *     r0 = -0.5,   kr = 1.15 pi,   kz = 1.15,
 *
 * and the source is whatever makes that psi solve the equation.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/*
 * SIGN AND SCALING CONVENTION -- read this before using the class.
 *
 * The equation, refs/HDG-GradShafranov.pdf eq (3a), is
 *
 *     -div_bar( ( 1/r ) grad_bar psi ) = F( r, z, psi ) / r,
 *
 * equivalently -Delta*( psi ) = F, with
 *
 *     Delta*( psi ) := r d_r( ( 1/r ) d_r psi ) + d_zz psi.
 *
 * f() BELOW RETURNS F, NOT F/r. That matches the paper -- its Example 5 prints
 * "the effective source term was F( r, z, psi ) := ..." for exactly the F of
 * eq (3a) -- and it matches meq::Source::f(), whose header says in as many
 * words that "the 1/r on the right hand side belongs to the weak form, not
 * here". meq::GradShafranovSolver applies the 1/r itself.
 *
 * The pre-port version of this file returned F/r from its operator(), because
 * the old driver fed the right hand side of (3a) directly to the assembly.
 * That factor is worth being loud about: an extra or missing r here converges,
 * at the full k+1 rate, to the wrong function, and no order-of-accuracy study
 * can see it. Measured, with the /r reinstated: the L2 error in psi sits at
 * 1.0e-1 and in q at 6.9e-1 through four dyadic refinements, rate 0.00.
 * deltaStarFD() below, and the first test case in NewtonConvergence.cpp that
 * uses it, check f() against Delta* by finite differences instead, which does
 * see it.
 *
 * The source is, with sc := sin( kr ( r + r0 ) ) cos( kz z ) the exact flux,
 *
 *     F( r, z, psi ) = ( kr^2 + kz^2 ) psi
 *                      + ( kr / r ) cos( kr ( r + r0 ) ) cos( kz z )
 *                      + r ( sc^2 - psi^2 + exp( -sc ) - exp( -psi ) ).
 *
 * The first two terms are exactly -Delta*( psi ) evaluated at the manufactured
 * flux; the bracket vanishes identically at psi = sc and so does not disturb
 * the solution, while contributing a quadratic and an exponential dependence on
 * psi that a Picard or a Newton iteration has to work through. Note that this
 * means F is NOT a function of ( r, z ) alone even where the solution is known:
 * asking for f( r, z, 0 ) gives something quite different from f( r, z, psi ),
 * which is the whole point.
 *
 * Differentiating,
 *
 *     dF/dpsi = ( kr^2 + kz^2 ) + r ( -2 psi + exp( -psi ) ),
 *
 * which is what the Newton Jacobian carries as a mass term -( dF/dpsi )/r on
 * the potential block.
 */

/// The manufactured non-linear benchmark of refs/HDG-GradShafranov.pdf
/// Example 5.
///
/// r must be strictly positive: the source carries a 1/r, as the operator
/// does. r0 is an offset of the radial coordinate that places the arch of the
/// sine inside the domain, and is NOT a major radius.
class ManufacturedNonlinear
{
	public:
		/// @param r0In  the radial offset r0 inside sin( kr ( r + r0 ) ).
		/// @param krIn  the radial wavenumber kr.
		/// @param kzIn  the axial wavenumber kz.
		ManufacturedNonlinear( double r0In, double krIn, double kzIn )
			: r0Value( r0In ), krValue( krIn ), kzValue( kzIn )
		{
		}

		/// The paper's parameters: r0 = -0.5, kr = 1.15 pi, kz = 1.15. These
		/// are the values examples/manufactured.toml carries, and with the
		/// box in that file the flux runs over roughly one arch of the sine
		/// with a critical point inside the domain, which is the qualitative
		/// behaviour the paper wanted from a flux function.
		static ManufacturedNonlinear example5()
		{
			double const pi = 3.14159265358979323846;
			return ManufacturedNonlinear( -0.5, 1.15*pi, 1.15 );
		}

		/// The poloidal flux function.
		double psi( double r, double z ) const
		{
			return std::sin( krValue*( r + r0Value ) )*std::cos( kzValue*z );
		}

		/// grad_bar( psi ) = ( d_r psi, d_z psi ). Not the HDG flux: that is
		/// this divided by r, see flux().
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			dPsiDr =  krValue*std::cos( krValue*( r + r0Value ) )*std::cos( kzValue*z );
			dPsiDz = -kzValue*std::sin( krValue*( r + r0Value ) )*std::sin( kzValue*z );
		}

		/// The HDG flux q = grad_bar( psi ) / r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// The Grad-Shafranov source F( r, z, psi ), NOT F/r. See the
		/// convention note above before using it.
		double f( double r, double z, double psiValue ) const
		{
			double const sc = psi( r, z );
			double const linear = ( krValue*krValue + kzValue*kzValue )*psiValue;
			double const geometric = ( krValue/r )
			                         *std::cos( krValue*( r + r0Value ) )
			                         *std::cos( kzValue*z );
			double const nonlinear = sc*sc - psiValue*psiValue
			                         + std::exp( -sc ) - std::exp( -psiValue );
			return linear + geometric + r*nonlinear;
		}

		/// dF/dpsi at fixed ( r, z ). Required by the Newton solve, and the
		/// only thing in this file that the Solov'ev benchmark cannot test:
		/// there dF/dpsi is identically zero.
		double dFdPsi( double r, double, double psiValue ) const
		{
			return krValue*krValue + kzValue*kzValue
			       + r*( -2.0*psiValue + std::exp( -psiValue ) );
		}

		/// Delta*( psi ), by central differences of psi(). Used to check the
		/// scaling convention documented above -- Delta*( psi ) must equal
		/// -f( r, z, psi( r, z ) ) -- without depending on the hand-derived
		/// gradients or on a solver.
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			// r d_r( ( 1/r ) d_r psi ) as a second difference of the inner
			// quantity, plus d_zz psi. The same arrangement as
			// Soloviev.hpp::deltaStarFD().
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) )/( 2.0*h );
			double const dZZ = ( psi( r, z + h ) - 2.0*psi( r, z ) + psi( r, z - h ) )
			                   /( h*h );

			return r*dRInner + dZZ;
		}

		/// The radial offset r0.
		double r0() const
		{
			return r0Value;
		}

		/// The radial wavenumber kr.
		double kr() const
		{
			return krValue;
		}

		/// The axial wavenumber kz.
		double kz() const
		{
			return kzValue;
		}

	private:
		double r0Value;
		double krValue;
		double kzValue;
};

}
}

#endif // MEQ_TESTS_MANUFACTUREDNONLINEAR_HPP
