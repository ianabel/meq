#ifndef MEQ_TESTS_MILLERDSHAPE_HPP
#define MEQ_TESTS_MILLERDSHAPE_HPP

/*
 * Example 6 of refs/HDG-GradShafranov.pdf, section 4.2: a polynomial
 * non-linearity on a smooth D-shaped domain given by a Miller parametrisation.
 *
 * Transcribed from the RENDERED page 17 of that PDF. The paper prints, verbatim:
 *
 *     r( t ) = 1 + eps cos( t + arcsin( delta sin( t ) ) ),
 *     z( t ) = eps kappa sin( t ),                          t in [0, 2 pi)
 *
 * with eps = 0.32, delta = 0.33, kappa = 1.7, "corresponding to an ITER-like
 * configuration", and
 *
 *     p = ( psi/2 )( 1 + 2 psi^2/3 - psi^5/5 )   and   g = 0,
 *
 *     F( r, z, psi ) = r^2 ( 1 - ( 1 - psi^2 )^2/2 ).
 *
 * THE PAPER'S p AND ITS F ARE INCONSISTENT, and F is the one to keep.
 *
 *     printed p = psi/2 + psi^3/3 - psi^6/10     =>  p' = 1/2 + psi^2 - 3 psi^5/5
 *     printed F/r^2 = 1 - ( 1 - psi^2 )^2/2      =   1/2 + psi^2 - psi^4/2
 *
 * The two disagree in the last term, -0.6 psi^5 against -0.5 psi^4. Changing the
 * printed psi^5/5 to psi^4/5 removes the discrepancy exactly:
 *
 *     p = ( psi/2 )( 1 + 2 psi^2/3 - psi^4/5 ) = psi/2 + psi^3/3 - psi^5/10
 *     p' = 1/2 + psi^2 - psi^4/2                                    == F/r^2
 *
 * so the exponent in p is a typesetting slip and F is self-consistent. F is what
 * the solver is fed and F is what is implemented; p() below is the corrected
 * profile, and the test suite asserts F == r^2 p' so the correction is checked
 * rather than asserted. This is the third transcription defect found in this pair
 * of papers, after the Solov'ev source sign and the printed NSTX coefficients --
 * see CLAUDE.md.
 *
 * A SECOND, HARMLESS DIFFERENCE FROM CERFON & FREIDBERG. The Miller form in
 * refs/CerfonFreidberg.pdf eq (9) is
 *
 *     x = 1 + eps cos( tau + alpha sin( tau ) ),   sin( alpha ) = delta,
 *
 * that is arcsin( delta ) TIMES sin( tau ), where Example 6 prints
 * arcsin( delta sin( t ) ). The two agree exactly at t = 0, +-pi/2 and pi -- so
 * both give the same eps, delta and kappa -- and differ by at most 6e-4 in r
 * anywhere else, which is 0.2 per cent of the minor radius. What is implemented
 * is the form Example 6 prints. Both are available; see boundaryPoint() and
 * boundaryPointCerfonFreidberg().
 *
 * NO EXACT SOLUTION. The paper says so explicitly and measures convergence by
 * successive differences instead: "For this final example we do not impose a
 * manufactured solution; instead we set homogeneous boundary conditions and
 * estimate the convergence by measuring the maximum difference between two
 * consecutive levels of refinement, Delta^k_j( f ) := || f_h^k - f_h^(k-1) ||_j."
 * So there is no psi(), gradPsi(), flux() or deltaStarFD() here.
 *
 * F( r, 0 ) = r^2/2, WHICH IS NOT ZERO -- unlike every source in
 * PressurePedestal.hpp and TransportBarrier.hpp. So this one IS well posed with
 * the paper's own homogeneous Dirichlet data: psi == 0 does not solve it, and
 * Newton started from zero has somewhere to go. It is the only benchmark of the
 * five new non-linear ones that can be posed exactly as its paper poses it.
 *
 * AND THE NON-LINEARITY IS NUMERICALLY NEGLIGIBLE, which is worth knowing before
 * reading anything into the Newton counts. On this geometry the source is
 * O( 1/2 ) and the minor radius is 0.32, so psi comes out around 1.3e-2; the
 * psi^2 term of F/r^2 is then about 3e-4 of the 1/2, and dF/dpsi = 2 r^2 psi
 * ( 1 - psi^2 ) is about 3e-2. This is a very slightly perturbed LINEAR problem.
 * It exercises the geometry and the self-convergence measurement; it does not
 * exercise Newton, and NewtonConvergence.cpp's Example 5 and
 * PedestalConvergence.cpp's four sources are what do.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/// The Miller D-shaped boundary of refs/HDG-GradShafranov.pdf Example 6, and the
/// Grad-Shafranov source posed inside it.
class MillerDShape
{
	public:
		/// @param epsIn    inverse aspect ratio a/R0. 0.32 in the paper.
		/// @param deltaIn  triangularity. 0.33.
		/// @param kappaIn  elongation. 1.7.
		MillerDShape( double epsIn, double deltaIn, double kappaIn )
			: epsValue( epsIn ), deltaValue( deltaIn ), kappaValue( kappaIn )
		{
		}

		/// Example 6 as published: an ITER-like D shape.
		static MillerDShape example6()
		{
			return MillerDShape( 0.32, 0.33, 1.7 );
		}

		/*
		 * ---------------------------------------------------------------
		 * The source
		 * ---------------------------------------------------------------
		 */

		/// The corrected pressure profile, p = ( psi/2 )( 1 + 2 psi^2/3
		/// - psi^4/5 ). See the header comment: the paper prints psi^5/5, which
		/// does not differentiate to its own F.
		static double p( double psi )
		{
			return 0.5*psi*( 1.0 + 2.0*psi*psi/3.0
			                 - psi*psi*psi*psi/5.0 );
		}

		/// dp/dpsi = 1/2 + psi^2 - psi^4/2, which is F/r^2 exactly.
		static double pPrime( double psi )
		{
			double const psi2 = psi*psi;
			return 0.5 + psi2 - 0.5*psi2*psi2;
		}

		/// d2p/dpsi2 = 2 psi - 2 psi^3 = 2 psi ( 1 - psi^2 ).
		static double pDoublePrime( double psi )
		{
			return 2.0*psi*( 1.0 - psi*psi );
		}

		/// F = r^2 ( 1 - ( 1 - psi^2 )^2/2 ), written as the paper writes it
		/// rather than as the expanded polynomial, so that the two forms can be
		/// compared in the test.
		double f( double r, double /*z*/, double psi ) const
		{
			double const oneMinusPsiSquared = 1.0 - psi*psi;
			return r*r*( 1.0 - 0.5*oneMinusPsiSquared*oneMinusPsiSquared );
		}

		/// dF/dpsi = 2 r^2 psi ( 1 - psi^2 ) = r^2 p''.
		double dFdPsi( double r, double /*z*/, double psi ) const
		{
			return r*r*pDoublePrime( psi );
		}

		/*
		 * ---------------------------------------------------------------
		 * The boundary
		 * ---------------------------------------------------------------
		 */

		/// The boundary curve as Example 6 prints it:
		/// r = 1 + eps cos( t + arcsin( delta sin t ) ), z = eps kappa sin t.
		void boundaryPoint( double t, double &r, double &z ) const
		{
			r = 1.0 + epsValue*std::cos( t + std::asin( deltaValue*std::sin( t ) ) );
			z = epsValue*kappaValue*std::sin( t );
		}

		/// The same curve in Cerfon & Freidberg's eq (9) form, with
		/// alpha = arcsin( delta ) multiplying sin( tau ). Provided for
		/// comparison; the test asserts the two agree at t = 0, pi/2, pi, 3 pi/2
		/// and differ by less than 1e-3 elsewhere.
		void boundaryPointCerfonFreidberg( double t, double &r, double &z ) const
		{
			double const alpha = std::asin( deltaValue );
			r = 1.0 + epsValue*std::cos( t + alpha*std::sin( t ) );
			z = epsValue*kappaValue*std::sin( t );
		}

		/**
		 * A level set of the D-shaped region: negative inside, zero on the
		 * boundary curve, positive outside.
		 *
		 * The curve is star shaped about ( 1, 0 ), so for each polar angle about
		 * that point there is exactly one boundary point, and
		 *
		 *     levelSet( r, z ) = | ( r, z ) - ( 1, 0 ) | - R_b( angle )
		 *
		 * is continuous everywhere and vanishes exactly on the curve. R_b is
		 * found by bisecting on the curve parameter, which is legitimate because
		 * the polar angle of the boundary point is a strictly increasing function
		 * of t -- checked in the test, not assumed.
		 *
		 * This is NOT a signed distance function: it is the radial gap, which is
		 * larger. It is what mfem::MarkLevelSetSubdomain and
		 * mfem::VertexConePath need -- a sign at the vertices and a root along a
		 * ray -- and nothing here depends on it being a distance.
		 */
		double levelSet( double r, double z ) const
		{
			double const dr = r - 1.0;
			double const rho = std::hypot( dr, z );
			if ( rho < 1.0e-14 )
				return -epsValue;

			double angle = std::atan2( z, dr );
			if ( angle < 0.0 )
				angle += 2.0*M_PI;

			double const t = parameterAtAngle( angle );
			double br, bz;
			boundaryPoint( t, br, bz );
			return rho - std::hypot( br - 1.0, bz );
		}

		/// The polar angle about ( 1, 0 ) of the boundary point at parameter @a t,
		/// in [ 0, 2 pi ). Strictly increasing in t for a star-shaped curve.
		double angleAtParameter( double t ) const
		{
			double br, bz;
			boundaryPoint( t, br, bz );
			double angle = std::atan2( bz, br - 1.0 );
			if ( angle < 0.0 )
				angle += 2.0*M_PI;
			return angle;
		}

		double eps() const
		{
			return epsValue;
		}

		double delta() const
		{
			return deltaValue;
		}

		double kappa() const
		{
			return kappaValue;
		}

	private:
		/// Invert angleAtParameter() by bisection. Sixty halvings takes the
		/// bracket to below 1e-17 of 2 pi, so the result is exact to double
		/// precision and the cost is a fixed sixty evaluations rather than an
		/// iteration count that depends on the point.
		double parameterAtAngle( double angle ) const
		{
			double low = 0.0;
			double high = 2.0*M_PI;
			for ( int i = 0; i < 60; ++i )
			{
				double const middle = 0.5*( low + high );
				if ( angleAtParameter( middle ) < angle )
					low = middle;
				else
					high = middle;
			}
			return 0.5*( low + high );
		}

		double epsValue, deltaValue, kappaValue;
};

}
}

#endif // MEQ_TESTS_MILLERDSHAPE_HPP
