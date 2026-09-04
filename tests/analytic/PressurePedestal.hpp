#ifndef MEQ_TESTS_PRESSUREPEDESTAL_HPP
#define MEQ_TESTS_PRESSUREPEDESTAL_HPP

/*
 * The three stiff source terms built on a pressure pedestal, from
 * refs/HDG-GradShafranov-Adaptive.pdf sections 4.2, 4.4 and 4.5.
 *
 *   4.2  eq (23), (24)   a pressure pedestal
 *   4.4  eq (26)         a current hole      = the pedestal source plus a term
 *   4.5  eq (27)         an internal layer   = the pedestal source plus a term
 *
 * Section 4.3, the internal transport barrier, has a different pressure profile
 * entirely and lives in TransportBarrier.hpp.
 *
 * All three were transcribed from RENDERED pages 16, 18 and 19 of that PDF, not
 * from pdftotext -- which is recorded in CLAUDE.md as having silently dropped
 * every minus sign in one of these papers and an epsilon in another.
 *
 * NONE OF THESE HAS AN EXACT SOLUTION, so there is no psi(), no gradPsi(), no
 * flux() and no deltaStarFD() here: there is nothing for them to return. The
 * paper judges these four by its residual error estimator, which MEQ does not
 * have yet (stage 6). What they are used for instead is a SELF-convergence
 * study -- successive refinements compared against each other, as
 * refs/HDG-GradShafranov.pdf does for its own Example 6 -- and as a stress test
 * of the Newton path, which is the thing MEQ does differently from both papers.
 *
 * THE CONVENTION MAP.
 *
 * The paper's eq (1) is -Delta*(psi) = F with F := mu0 r^2 dp/dpsi + g dg/dpsi,
 * which is MEQ's convention exactly (see CLAUDE.md, "The equation being
 * solved"). Sections 4.2 and 4.4 take g = constant, so g dg/dpsi = 0 and
 * F = mu0 r^2 p'(psi) with mu0 = 1 in the paper's normalisation. Checked term by
 * term: differentiating eq (23) gives
 *
 *     p'(psi) = 2 c2 psi ( 1 - e ) + ( c1 + c2 psi^2 )( 2 psi/sigma^2 ) e
 *             = 2 psi [ c2 ( 1 - e ) + ( 1/sigma^2 )( c1 + c2 psi^2 ) e ],
 *             e := exp( -( psi/sigma )^2 )
 *
 * which is eq (24) divided by r^2. So the paper's F really is mu0 r^2 p', and
 * pPrime() and f() below are related by exactly that factor -- asserted in the
 * test suite rather than trusted here.
 *
 * Unlike the Solov'ev source, there is no sign disagreement between the two
 * papers to resolve here: eq (24) is consistent with eq (1) as printed.
 *
 * ALL OF THEM VANISH AT psi = 0, WHICH IS THE INTERESTING PART.
 *
 * F( r, 0 ) = 0 for every source in this file and for the transport barrier as
 * well. On the paper's own domain -- an ITER-like region with psi = 0 on the
 * plasma boundary -- that means psi == 0 is an exact solution of the boundary
 * value problem, and the physically interesting equilibrium is a SECOND branch.
 * refs/HDG-GradShafranov.pdf says as much, in one word: its Algorithm 2 begins
 * "psi^0 ;   // Non-trivial initial guess".
 *
 * MEQ's Newton iteration starts from the Dirichlet data and nothing else, so
 * with homogeneous data it starts exactly on the trivial branch, where the
 * residual is identically zero, and reports convergence in zero iterations
 * having produced psi == 0. That is asserted in
 * tests/convergence/PedestalConvergence.cpp -- deliberately, so the finding is
 * recorded in code and not just in a comment -- and it is why the benchmarks
 * there are posed with a non-homogeneous datum instead. See that file.
 *
 * WHAT MAKES THESE HARD. sigma^2 = 0.005, so sigma = 0.0707: the pedestal is a
 * layer of that width in psi, across which dF/dpsi swings by a factor of
 * c1/( sigma^2 c2 ) = 800. At psi = 0 the Jacobian's mass term is
 * 2 r^2 c1/sigma^2 = 320 r^2, far above the first Dirichlet eigenvalue of any
 * domain used here, so the Newton Jacobian is indefinite. That is exactly the
 * regime CLAUDE.md's note on KINSolver( KIN_LINESEARCH ) was written for.
 */

#include <cmath>

namespace meq
{
namespace analytic
{

/**
 * The pedestal core, eq (23) and eq (24) of
 * refs/HDG-GradShafranov-Adaptive.pdf section 4.2.
 *
 *     p( psi )    = ( c1 + c2 psi^2 )( 1 - exp( -( psi/sigma )^2 ) )
 *     F( r, psi ) = 2 r^2 psi ( c2 ( 1 - e ) + ( 1/sigma^2 )( c1 + c2 psi^2 ) e )
 *
 * Sections 4.4 and 4.5 reuse this with their own constants, which is why it is
 * a class of its own rather than a set of literals.
 */
class PressurePedestal
{
	public:
		/// @param c1In     the flat-top pressure. 0.8 in section 4.2.
		/// @param c2In     the quadratic correction. 0.2 in section 4.2.
		/// @param sigmaSquaredIn  the SQUARE of the pedestal width in psi. The
		///                        paper quotes sigma^2, not sigma, and 0.005 is
		///                        what makes this the stiffest case available.
		PressurePedestal( double c1In, double c2In, double sigmaSquaredIn )
			: c1Value( c1In ), c2Value( c2In ), sigmaSquaredValue( sigmaSquaredIn )
		{
		}

		/// Section 4.2 as published: c1 = 0.8, c2 = 0.2, sigma^2 = 0.005.
		static PressurePedestal pedestal()
		{
			return PressurePedestal( 0.8, 0.2, 0.005 );
		}

		/// The pressure, eq (23). Not used by the solver -- F is -- but it is
		/// what f() is checked against.
		double p( double psi ) const
		{
			return ( c1Value + c2Value*psi*psi )*( 1.0 - decay( psi ) );
		}

		/// dp/dpsi, differentiated in closed form. f() is r^2 times this.
		double pPrime( double psi ) const
		{
			double const e = decay( psi );
			return 2.0*psi*( c2Value*( 1.0 - e )
			                 + ( c1Value + c2Value*psi*psi )*e/sigmaSquaredValue );
		}

		/// d2p/dpsi2. dFdPsi() is r^2 times this.
		double pDoublePrime( double psi ) const
		{
			double const e = decay( psi );
			double const s = sigmaSquaredValue;
			double const psi2 = psi*psi;

			// d/dpsi of 2 c2 psi ( 1 - e ), using de/dpsi = -( 2 psi/s ) e.
			double const first = 2.0*c2Value*( 1.0 - e ) + 4.0*c2Value*psi2*e/s;

			// d/dpsi of ( 2/s ) psi ( c1 + c2 psi^2 ) e.
			double const second = 2.0*e*( c1Value + 3.0*c2Value*psi2
			                              - 2.0*psi2*( c1Value + c2Value*psi2 )/s )/s;

			return first + second;
		}

		/// F = mu0 r^2 p'( psi ), eq (24), with mu0 = 1. Returns F, not F/r:
		/// meq::Source::f() is documented to be F as eq (2) writes it.
		double f( double r, double /*z*/, double psi ) const
		{
			return r*r*pPrime( psi );
		}

		/// dF/dpsi = r^2 p''( psi ). Differentiated, not finite-differenced: a
		/// wrong Jacobian does not move the converged answer, it only wrecks the
		/// path to it, and no convergence rate can see the difference.
		double dFdPsi( double r, double /*z*/, double psi ) const
		{
			return r*r*pDoublePrime( psi );
		}

		double c1() const
		{
			return c1Value;
		}

		double c2() const
		{
			return c2Value;
		}

		double sigmaSquared() const
		{
			return sigmaSquaredValue;
		}

		/// exp( -( psi/sigma )^2 ), the factor every term is built from.
		double decay( double psi ) const
		{
			return std::exp( -psi*psi/sigmaSquaredValue );
		}

	private:
		double c1Value, c2Value, sigmaSquaredValue;
};

/**
 * The current hole, eq (26) of section 4.4:
 *
 *     F( r, psi ) = < the pedestal source, with c1, c2, sigma_1 >
 *                   + c3 ( 1 - exp( -( psi/sigma_2 )^2 ) ) cos( c4 psi )
 *
 * with c1 = 0.4, c2 = 0.1, c3 = -18, c4 = 10 pi, sigma_1^2 = 5e-3 and
 * sigma_2^2 = 3e-3. Note that the added term carries no r^2: it is a pure
 * function of psi, which is what mu0 r^2 p' would not be. It still fits the
 * canonical form, as a g dg/dpsi contribution.
 *
 * mu0 J_phi = F/r, so the oscillating cosine is what drives the toroidal
 * current close to zero over an extended core region -- the "hole" -- while
 * leaving sharp peaks near the boundary.
 */
class CurrentHole
{
	public:
		CurrentHole( PressurePedestal const &pedestalIn, double c3In, double c4In,
		             double sigma2SquaredIn )
			: pedestalPart( pedestalIn ), c3Value( c3In ), c4Value( c4In ),
			  sigma2SquaredValue( sigma2SquaredIn )
		{
		}

		/// Section 4.4 as published.
		static CurrentHole currentHole()
		{
			return CurrentHole( PressurePedestal( 0.4, 0.1, 5.0e-3 ),
			                    -18.0, 10.0*M_PI, 3.0e-3 );
		}

		double f( double r, double z, double psi ) const
		{
			return pedestalPart.f( r, z, psi )
			       + c3Value*( 1.0 - decay2( psi ) )*std::cos( c4Value*psi );
		}

		double dFdPsi( double r, double z, double psi ) const
		{
			double const e2 = decay2( psi );
			double const added = c3Value*( 2.0*psi*e2*std::cos( c4Value*psi )
			                               /sigma2SquaredValue
			                               - c4Value*( 1.0 - e2 )
			                                 *std::sin( c4Value*psi ) );
			return pedestalPart.dFdPsi( r, z, psi ) + added;
		}

		PressurePedestal const &pedestal() const
		{
			return pedestalPart;
		}

	private:
		double decay2( double psi ) const
		{
			return std::exp( -psi*psi/sigma2SquaredValue );
		}

		PressurePedestal pedestalPart;
		double c3Value, c4Value, sigma2SquaredValue;
};

/**
 * The internal layer, eq (27) of section 4.5:
 *
 *     F( r, psi ) = < the pedestal source, with c1, c2, sigma_1 >
 *                   + c3 ( 1 - exp( -( psi/sigma_1 )^2 ) )
 *                     exp( -( 1 - r - psi )^2/sigma_2^2 )
 *
 * with c1 = 0.8, c2 = 0.2, c3 = 15, sigma_1^2 = 5e-3, sigma_2^2 = 7.5e-4. Note
 * that BOTH factors of the added term use sigma_1 in the first exponential and
 * sigma_2 only in the second; that is what the rendered page 19 prints.
 *
 * THIS ONE IS NOT A PHYSICAL EQUILIBRIUM, and the paper says so: "This source
 * term is not physically relevant for magnetic confinement fusion applications,
 * because it cannot be cast in the canonical form of the source in (1) due to
 * the explicit appearance of the coordinate r in the argument of the last
 * exponential." There is no p( psi ) and no g( psi ) that produce it, so this
 * class has no p() -- only f() and dFdPsi(). It is here because it is a good
 * benchmark for detecting a localised internal layer: sigma_2 = 0.027, so the
 * added term is a ridge of that width along the curve r + psi = 1.
 *
 * For MEQ it is also the one case where dF/dpsi and dF/dr are genuinely
 * independent, since r enters the nonlinearity rather than multiplying it.
 * Nothing in the solver cares -- only dF/dpsi is ever asked for -- but it means
 * a fixture that got the r-dependence wrong would still pass a dF/dpsi check,
 * so f() is also compared against the printed expression term by term in the
 * test.
 */
class InternalLayer
{
	public:
		InternalLayer( PressurePedestal const &pedestalIn, double c3In,
		               double sigma2SquaredIn )
			: pedestalPart( pedestalIn ), c3Value( c3In ),
			  sigma2SquaredValue( sigma2SquaredIn )
		{
		}

		/// Section 4.5 as published.
		static InternalLayer internalLayer()
		{
			return InternalLayer( PressurePedestal( 0.8, 0.2, 5.0e-3 ),
			                      15.0, 7.5e-4 );
		}

		double f( double r, double z, double psi ) const
		{
			return pedestalPart.f( r, z, psi )
			       + c3Value*( 1.0 - pedestalPart.decay( psi ) )*ridge( r, psi );
		}

		double dFdPsi( double r, double z, double psi ) const
		{
			double const e1 = pedestalPart.decay( psi );
			double const g = ridge( r, psi );
			double const u = 1.0 - r - psi;

			// d/dpsi of ( 1 - e1 ) is ( 2 psi/sigma_1^2 ) e1, and d/dpsi of g is
			// ( 2 u/sigma_2^2 ) g -- the sign coming from du/dpsi = -1.
			double const added =
				c3Value*( 2.0*psi*e1*g/pedestalPart.sigmaSquared()
				          + ( 1.0 - e1 )*2.0*u*g/sigma2SquaredValue );

			return pedestalPart.dFdPsi( r, z, psi ) + added;
		}

		PressurePedestal const &pedestal() const
		{
			return pedestalPart;
		}

	private:
		/// exp( -( 1 - r - psi )^2/sigma_2^2 ). The r inside is what makes this
		/// source non-physical.
		double ridge( double r, double psi ) const
		{
			double const u = 1.0 - r - psi;
			return std::exp( -u*u/sigma2SquaredValue );
		}

		PressurePedestal pedestalPart;
		double c3Value, sigma2SquaredValue;
};

}
}

#endif // MEQ_TESTS_PRESSUREPEDESTAL_HPP
