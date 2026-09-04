#ifndef MEQ_TESTS_SIMILARITYEXPONENTIAL_HPP
#define MEQ_TESTS_SIMILARITYEXPONENTIAL_HPP

/*
 * An EXACT solution of a genuinely nonlinear Grad-Shafranov equation, from
 * Kaltsas & Throumoulopoulos, Phys. Lett. A 380 (2016) 3373, section 3.2,
 * eq (22) -- refs/GS-SimilarityReduction.pdf, doi 10.1016/j.physleta.2016.08.011.
 *
 * WHY THIS IS A DIFFERENT KIND OF TEST FROM ManufacturedNonlinear.hpp.
 *
 * Example 5 of HDG-GradShafranov.pdf is a MANUFACTURED solution: a convenient
 * psi is chosen and F is then constructed, by substitution, to make it solve the
 * equation. The resulting F is not of the form any physical profile would give.
 *
 * This one is the other way round. The free function is chosen -- an exponential
 * poloidal current profile, f(u) = f0 exp(n u) -- and the solution follows from
 * a similarity reduction. So it tests the solver against a nonlinear equation
 * somebody might actually pose, rather than against an equation reverse
 * engineered from an answer. Nothing else found in the literature offers an
 * exact solution, in elementary functions, of a nonlinear Grad-Shafranov
 * equation.
 *
 * THE CONVENTION MAP, which is where the care is needed.
 *
 * Their eq (1) is written
 *
 *     d_rr u - (1/r) d_r u + d_zz u + f(u) + g(u) r^2 = 0,
 *
 * and the first three terms ARE Delta*(u):
 *
 *     Delta* u = r d_r( (1/r) d_r u ) + d_zz u
 *              = d_rr u - (1/r) d_r u + d_zz u.
 *
 * So their equation is -Delta*(u) = f(u) + g(u) r^2, which is MEQ's
 * -Delta*(psi) = F with
 *
 *     F( r, z, psi ) = f(psi) ( 1 + epsilon r^2 ),     epsilon := a^2 / b^2.
 *
 * The reduction requires g(u) = epsilon f(u) with epsilon = a^2/b^2 (their
 * page 2, below eq (4)). NOTE: pdftotext renders that line as "g(u) = f(u)",
 * having dropped the epsilon -- the same extraction unreliability that dropped
 * every minus sign in the Solov'ev coefficients. Read the rendered page.
 *
 * Also note u is their eq (2) transform of the poloidal flux, which reduces to
 * u = psi when the poloidal Mach number vanishes. There is no flow here, so
 * u = psi throughout.
 *
 * THE SOLUTION, their eq (22), on the + branch of the reduction variable:
 *
 *     psi(r,z) = (1/n) ln[ K sech^2( m ( a r^2 + 2 b z + cTilde ) ) ]
 *     K := b^2 c n / ( 2 f0 ),      m := n sqrt(c) / 4
 *
 * Checked by hand against their eq (11), w''(x) + b^-2 f(w(x)) = 0 with
 * x = a r^2/2 + b z: with s := m( 2x + cTilde ),
 *
 *     w''  = -( 8 m^2 / n ) sech^2 s,     f(w) = ( b^2 c n / 2 ) sech^2 s,
 *
 * and 8m^2/n = n c / 2, so the two cancel exactly. deltaStarFD() below
 * re-checks it numerically and the test suite asserts that.
 *
 * NO CLOSED FLUX SURFACES, and it does not matter. The paper is explicit that
 * its nonlinear solutions "do not form closed surfaces" -- the level sets here
 * are parabolas in (r^2, z). For a fixed-boundary benchmark that is irrelevant:
 * the domain is a rectangle and the Dirichlet data is the exact solution
 * restricted to its boundary, exactly as ManufacturedNonlinear.hpp is used. It
 * would matter only for a psi = 0 plasma boundary, which this is not for.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/// The exponential-profile similarity solution of Kaltsas & Throumoulopoulos
/// eq (22). r must be positive; the solution itself is regular everywhere.
class SimilarityExponential
{
	public:
		/// @param nIn   exponent of the free function f(u) = f0 exp( n u ).
		///              dF/dpsi = n F, so n sets how nonlinear the problem is.
		/// @param f0In  amplitude of the free function.
		/// @param aIn   coefficient of r^2/2 in the reduction variable.
		/// @param bIn   coefficient of z in the reduction variable.
		/// @param cIn   the integration constant of their eq (5); must be > 0
		///              for the solution to be real on this branch.
		/// @param cTildeIn  the translation constant, shifting the solution in z.
		SimilarityExponential( double nIn, double f0In, double aIn, double bIn,
		                       double cIn, double cTildeIn )
			: nValue( nIn ), f0Value( f0In ), aValue( aIn ), bValue( bIn ),
			  cValue( cIn ), cTildeValue( cTildeIn )
		{
		}

		/// A well-conditioned choice for the benchmark rectangle
		/// [0.6, 1.4] x [-0.6, 0.6].
		///
		/// n = 3 was chosen by measurement over that box: psi ranges over
		/// [-0.757, +0.366], so exp( n psi ) varies by a factor of 29 and F
		/// over [0.15, 3.57]. That is a genuinely nonlinear problem rather than
		/// a perturbation of a linear one. Larger n keeps going -- 139x at
		/// n = 4, 675x at n = 5 -- but F falls to 0.011 at n = 5, which starts
		/// to be a conditioning question rather than a harder test of Newton.
		static SimilarityExponential benchmark()
		{
			// n = 3, f0 = 1/2, a = b = c = 1  =>  K = 3, m = 3/4.
			return SimilarityExponential( 3.0, 0.5, 1.0, 1.0, 1.0, 0.0 );
		}

		/// The poloidal flux.
		double psi( double r, double z ) const
		{
			return std::log( kConstant()*sech2( argument( r, z ) ) )/nValue;
		}

		/// grad_bar(psi) = ( d_r psi, d_z psi ).
		///
		/// psi = (1/n)( ln K + 2 ln sech s ) with s = m( a r^2 + 2 b z + cT ),
		/// so d psi/d s = -(2/n) tanh s, and the chain rule supplies
		/// ds/dr = 2 m a r, ds/dz = 2 m b.
		void gradPsi( double r, double z, double &dPsiDr, double &dPsiDz ) const
		{
			double const s = argument( r, z );
			double const dPsiDs = -2.0*std::tanh( s )/nValue;
			double const m = mConstant();
			dPsiDr = dPsiDs*2.0*m*aValue*r;
			dPsiDz = dPsiDs*2.0*m*bValue;
		}

		/// The HDG flux q = grad_bar(psi) / r.
		void flux( double r, double z, double &qR, double &qZ ) const
		{
			gradPsi( r, z, qR, qZ );
			qR /= r;
			qZ /= r;
		}

		/// F = f(psi) ( 1 + epsilon r^2 ) with f(psi) = f0 exp( n psi ) and
		/// epsilon = a^2/b^2. Returns F, not F/r, as meq::Source::f() does.
		double f( double r, double /*z*/, double psiValue ) const
		{
			return f0Value*std::exp( nValue*psiValue )*( 1.0 + epsilon()*r*r );
		}

		/// dF/dpsi = n F. The whole psi-dependence is a single exponential, so
		/// the Jacobian is the residual times a constant -- which makes this a
		/// clean test of the Newton path without the algebra of Example 5.
		double dFdPsi( double r, double z, double psiValue ) const
		{
			return nValue*f( r, z, psiValue );
		}

		/// Delta*(psi) by central differences. The equation is
		/// -Delta*(psi) = F, so this must equal -f( r, z, psi(r,z) ).
		double deltaStarFD( double r, double z, double h = 1.0e-4 ) const
		{
			auto innerR = [ & ]( double rr )
			{
				return ( psi( rr + h, z ) - psi( rr - h, z ) )/( 2.0*h )/rr;
			};

			double const dRInner = ( innerR( r + h ) - innerR( r - h ) )/( 2.0*h );
			double const dZZ = ( psi( r, z + h ) - 2.0*psi( r, z ) + psi( r, z - h ) )
			                   /( h*h );

			return r*dRInner + dZZ;
		}

		double n() const
		{
			return nValue;
		}

		/// epsilon = a^2/b^2, the ratio the reduction fixes between the
		/// pressure and poloidal-current free functions.
		double epsilon() const
		{
			return aValue*aValue/( bValue*bValue );
		}

	private:
		/// s = m ( a r^2 + 2 b z + cTilde ).
		double argument( double r, double z ) const
		{
			return mConstant()*( aValue*r*r + 2.0*bValue*z + cTildeValue );
		}

		double kConstant() const
		{
			return bValue*bValue*cValue*nValue/( 2.0*f0Value );
		}

		double mConstant() const
		{
			return nValue*std::sqrt( cValue )/4.0;
		}

		static double sech2( double s )
		{
			double const c = std::cosh( s );
			return 1.0/( c*c );
		}

		double nValue, f0Value, aValue, bValue, cValue, cTildeValue;
};

}
}

#endif // MEQ_TESTS_SIMILARITYEXPONENTIAL_HPP
