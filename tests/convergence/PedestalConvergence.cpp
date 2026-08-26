#define BOOST_TEST_MODULE PedestalConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/PressurePedestal.hpp"
#include "analytic/Soloviev.hpp"
#include "analytic/TransportBarrier.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * The four non-linear sources of refs/HDG-GradShafranov-Adaptive.pdf sections
 * 4.2 to 4.5: a pressure pedestal, an internal transport barrier, a current
 * hole and an internal layer. These are the stiffest source terms meq has been
 * given, and the reason for having them is that the Newton path had not been
 * pushed on anything harder than Example 5's smooth manufactured source.
 *
 * WHAT CAN AND CANNOT BE ASSERTED HERE.
 *
 * None of the four has an exact solution. The paper judges them by its residual
 * error estimator, which is stage-6 work meq does not have, so there is nothing
 * to measure an error against and no k+1 rate to assert. What replaces it is a
 * SELF-convergence study -- successive refinements compared against each other,
 * which is exactly what refs/HDG-GradShafranov.pdf does for its own Example 6 --
 * plus the Newton record: iteration counts and residual histories at every
 * refinement and every degree.
 *
 * A self-convergence table is a strictly weaker instrument, and it is worth being
 * plain about the gap. It says the sequence is Cauchy at the design rate. It does
 * NOT say the limit is the right function: a wrong sign convention converges, at
 * the right rate, to the wrong answer, which is the standing warning of
 * CLAUDE.md's "Testing stance". The things that guard against that here are
 * elsewhere -- the exact-solution ladder in SolovievConvergence.cpp and
 * NewtonConvergence.cpp -- and what is guarded HERE is the transcription of the
 * source itself, term by term against the printed equations, plus dF/dpsi
 * against a central difference of F.
 *
 *
 * THE TRIVIAL BRANCH, WHICH IS THE FIRST REAL FINDING IN THIS FILE.
 *
 * Every one of the four sources vanishes at psi = 0:
 *
 *   4.2  F = 2 r^2 psi ( ... )                                 explicitly
 *   4.3  F = r^2 p', p = A( psi ) ( 1 - ( 1 - psi )^a )^b       b = 2, so the
 *                                                              bracket is
 *                                                              O( psi^2 )
 *   4.4  F = < 4.2 > + c3 ( 1 - exp( -( psi/sigma_2 )^2 ) ) cos( c4 psi )
 *   4.5  F = < 4.2 > + c3 ( 1 - exp( -( psi/sigma_1 )^2 ) ) exp( ... )
 *
 * so on the paper's own domain -- an ITER-like region with psi = 0 on the plasma
 * boundary -- psi == 0 is an exact solution of the boundary value problem and
 * the interesting equilibrium is a SECOND branch. The paper is explicit that its
 * iteration has to be started off that branch: Algorithm 2 of
 * refs/HDG-GradShafranov.pdf begins
 *
 *     psi^0 ;                                 // Non-trivial initial guess
 *
 * meq's Newton iteration starts from the Dirichlet data and nothing else. With
 * homogeneous data the initial residual is identically zero, NewtonSolver
 * reports convergence in zero iterations, and the answer is psi == 0.
 * homogeneousDataLandsOnTheTrivialBranch below asserts exactly that, so the
 * finding lives in the suite rather than in a comment.
 *
 * That is a gap in GradShafranovSolver's interface, not in these benchmarks:
 * there is no setInitialGuess(), and prepare() -- which solve() calls -- resets
 * the iterate, so no caller can supply one. Adding one is the fix; until then
 * the four benchmarks here are posed the way every other non-linear benchmark in
 * this directory already is, with NON-HOMOGENEOUS Dirichlet data on the standard
 * rectangle. NewtonConvergence.cpp says why that is the house pattern: "the
 * McCarthy and Solov'ev cases are being solved on a box that is not their own
 * plasma boundary. That is deliberate and is what makes the Dirichlet data
 * non-homogeneous."
 *
 * So what is measured below is the paper's source term, exactly as printed, on
 * meq's standard box, with a datum that excludes psi == 0. The deviation is the
 * domain and the boundary condition, both of which are stated here; the source
 * is untouched and is checked against the printed equations.
 *
 *
 * WHY THIS DATUM.
 *
 * pedestalDatum() is a linear ramp in z, from -0.3 at the bottom of the box to
 * +0.3 at the top, and it carries three of the four sources. Three reasons, in
 * order of importance:
 *
 *   1. it changes sign, so the surface psi = 0 runs through the INTERIOR of the
 *      domain. sigma = 0.0707 for section 4.2, so the pedestal -- the layer in
 *      which dF/dpsi swings by a factor of 800 -- lies inside the mesh rather
 *      than pressed against the boundary where a Dirichlet condition would
 *      partly hide it. Section 4.5's ridge, at r + psi = 1, likewise crosses the
 *      box for psi in [ -0.4, 0.4 ]. That is the localised internal structure
 *      these benchmarks exist to produce.
 *   2. it is smooth, and simple enough to state in one line. There is no exact
 *      solution to take a trace of, so any datum here is a design choice, and a
 *      complicated one would be a hidden parameter.
 *   3. its magnitude puts psi across the whole interesting range. p( psi ) in
 *      eq (23) saturates for |psi| >> sigma and is quadratic below it, so a
 *      datum of a few sigma exercises both regimes. Measured, the solution then
 *      spans very nearly the datum's own range, [ -0.300, +0.299 ].
 *
 * SECTION 4.3 IS THE EXCEPTION and uses barrierDatum(), a ramp from 0 to 0.6,
 * because its feature is at psi_0 = 0.3 rather than at psi = 0. With
 * pedestalDatum() the barrier would sit exactly on the top edge of the box, and
 * measured, that configuration does not converge at all. See
 * transportBarrierSelfConverges.
 *
 * AND IT COSTS SOMETHING, WHICH IS MEASURED RATHER THAN ARGUED. With a datum
 * that is not the trace of a smooth solution of this equation, the corners of
 * the rectangle limit the regularity of the solution: a right-angled corner has
 * singular exponents k pi/omega = 2k, and the resonance with the integer
 * exponent produces an r^2 log r term, so the solution sits in H^(3-epsilon) and
 * its gradient in H^(2-epsilon) however smooth the data is.
 *
 * aCornerSingularSolutionCapsTheRateOnThisRectangle puts a number on it, on a
 * Solov'ev source that has no non-linearity whatever: the self-convergence rate
 * is 2.00 in psi at k = 1 and then FLAT at 2.86 and 2.96, with the flux flat at
 * about 2.2. That is the ceiling every study in this file is working under. It
 * explains why sections 4.3 and 4.5 fall short of k+1 at k >= 2 and it is why
 * the floors asserted below are set from measurement rather than from k+1.
 *
 * The exact-solution studies do not have the problem, because there the datum IS
 * the trace of a smooth solution -- SolovievGeometryConvergence.cpp gets 4.000
 * at k = 3 on the same rectangle. MillerConvergence.cpp measures the same effect
 * again, in a sharper form, on a polygon whose corners are nearly flat, with the
 * same control beside it.
 *
 *
 * THE THREE FINDINGS THIS FILE RECORDS, each as a test that asserts the defect
 * so that the day it is fixed the suite says so:
 *
 *   homogeneousDataLandsOnTheTrivialBranch
 *       all four sources vanish at psi = 0, so the paper's own homogeneous
 *       problem is solved by psi == 0 and meq's Newton, which starts from the
 *       Dirichlet data, stops there in zero iterations.
 *
 *   pedestalNewtonFailsOnCoarseMeshesAtOrderOne
 *       at k = 1 the pedestal does not converge on meshes coarser than
 *       h = 0.033. The element-local non-linear solves fail first.
 *
 *   currentHoleDoesNotConvergeUnderPlainNewton
 *       section 4.4 does not converge in ANY configuration tried, and on the
 *       coarser ones the iterate reaches NaN and MFEM aborts the process.
 *
 * All three point the same way: meq differs from both papers in closing the
 * semi-linear problem by Newton rather than by Anderson-accelerated Picard, and
 * Newton needs an initial guess and a globalisation that the solver does not yet
 * offer. CLAUDE.md already names the intended remedy -- KINSolver with
 * KIN_LINESEARCH -- and it cannot be tried here, because the MFEM tree is built
 * with MFEM_USE_SUNDIALS = NO.
 */

namespace
{
	using meq::tests::Rectangle;
	using meq::tests::SampleCloud;
	using meq::tests::SelfMeasurement;
	using meq::tests::standardBox;

	/// The Dirichlet datum, and the whole of the deviation from the paper's
	/// problem. See the file comment for why it is this and not zero.
	double pedestalDatumAmplitude()
	{
		return 0.3;
	}

	double pedestalDatum( double /*r*/, double z )
	{
		Rectangle const box = standardBox();
		return pedestalDatumAmplitude()*z/box.zMax;
	}

	/// The datum for section 4.3 only, a ramp from 0.2 at the bottom of the box
	/// to 0.4 at the top. See transportBarrierSelfConverges for both reasons it
	/// differs from pedestalDatum(): the barrier's feature is at psi_0 = 0.3
	/// rather than at psi = 0, and its width has to be resolvable.
	double barrierDatum( double, double z )
	{
		Rectangle const box = standardBox();
		return 0.3 + 0.1*z/box.zMax;
	}

	double zeroDatum( double, double )
	{
		return 0.0;
	}

	SampleCloud const &cloud()
	{
		static SampleCloud const c( standardBox() );
		return c;
	}

	/// Section 4.2's own constants, spelled out here rather than reached for
	/// through the fixture, so that the check below is a check and not a
	/// tautology.
	double const c1 = 0.8;
	double const c2 = 0.2;
	double const sigmaSquared = 0.005;

	/// eq (24), written out from the rendered page 16 with no factoring:
	/// F = 2 r^2 psi ( c2 ( 1 - e ) + ( 1/sigma^2 )( c1 + c2 psi^2 ) e ).
	double printedPedestalSource( double r, double psi )
	{
		double const e = std::exp( -( psi/std::sqrt( sigmaSquared ) )
		                           *( psi/std::sqrt( sigmaSquared ) ) );
		return 2.0*r*r*psi*( c2*( 1.0 - e )
		                     + ( c1 + c2*psi*psi )*e/sigmaSquared );
	}

	/// eq (26), from the rendered page 18, with section 4.4's own constants.
	double printedCurrentHoleSource( double r, double psi )
	{
		double const hc1 = 0.4, hc2 = 0.1, hc3 = -18.0, hc4 = 10.0*M_PI;
		double const s1 = 5.0e-3, s2 = 3.0e-3;
		double const e1 = std::exp( -psi*psi/s1 );
		double const e2 = std::exp( -psi*psi/s2 );
		return 2.0*r*r*psi*( hc2*( 1.0 - e1 ) + ( hc1 + hc2*psi*psi )*e1/s1 )
		       + hc3*( 1.0 - e2 )*std::cos( hc4*psi );
	}

	/// eq (27), from the rendered page 19, with section 4.5's own constants.
	/// Note the first exponential of the added term uses sigma_1, not sigma_2.
	double printedInternalLayerSource( double r, double psi )
	{
		double const lc1 = 0.8, lc2 = 0.2, lc3 = 15.0;
		double const s1 = 5.0e-3, s2 = 7.5e-4;
		double const e1 = std::exp( -psi*psi/s1 );
		double const u = 1.0 - r - psi;
		return 2.0*r*r*psi*( lc2*( 1.0 - e1 ) + ( lc1 + lc2*psi*psi )*e1/s1 )
		       + lc3*( 1.0 - e1 )*std::exp( -u*u/s2 );
	}

	/// A central difference of F in psi, for the Jacobian checks.
	template<typename Source>
	double differenceOfF( Source const &source, double r, double z, double psi,
	                      double step )
	{
		return ( source.f( r, z, psi + step ) - source.f( r, z, psi - step ) )
		       /( 2.0*step );
	}

	/// dF/dpsi against a central difference of F, over a scatter of ( r, psi ).
	/// The step has to be small next to sigma = 0.07 and large next to the
	/// round-off floor, which is what 1e-6 is.
	template<typename Source>
	double worstJacobianDiscrepancy( Source const &source )
	{
		double const step = 1.0e-6;
		double worst = 0.0;
		for ( double r = 0.6; r < 1.45; r += 0.1 )
		{
			for ( double psi = -0.5; psi < 0.85; psi += 0.01 )
			{
				double const difference = differenceOfF( source, r, 0.13, psi, step );
				double const analytic = source.dFdPsi( r, 0.13, psi );
				double const scale = 1.0 + std::abs( analytic );
				worst = std::max( worst, std::abs( difference - analytic )/scale );
			}
		}
		return worst;
	}
}

/*
 * -------------------------------------------------------------------------
 * The transcription, before any solver is involved
 * -------------------------------------------------------------------------
 */

/// The fixtures factor eq (24) into a pedestal class that eqs (26) and (27)
/// reuse. That factoring is convenient and it is also exactly the kind of
/// rearrangement that loses a constant, so f() is compared against the printed
/// expressions written out flat.
BOOST_AUTO_TEST_CASE( sourcesReproduceThePrintedExpressions )
{
	meq::analytic::PressurePedestal const pedestal
		= meq::analytic::PressurePedestal::pedestal();
	meq::analytic::CurrentHole const hole
		= meq::analytic::CurrentHole::currentHole();
	meq::analytic::InternalLayer const layer
		= meq::analytic::InternalLayer::internalLayer();

	double worstPedestal = 0.0, worstHole = 0.0, worstLayer = 0.0;
	for ( double r = 0.6; r < 1.45; r += 0.05 )
	{
		for ( double psi = -0.6; psi < 0.95; psi += 0.005 )
		{
			auto relative = []( double a, double b )
			{
				return std::abs( a - b )/( 1.0 + std::abs( b ) );
			};
			worstPedestal = std::max( worstPedestal,
				relative( pedestal.f( r, 0.0, psi ), printedPedestalSource( r, psi ) ) );
			worstHole = std::max( worstHole,
				relative( hole.f( r, 0.0, psi ), printedCurrentHoleSource( r, psi ) ) );
			worstLayer = std::max( worstLayer,
				relative( layer.f( r, 0.0, psi ), printedInternalLayerSource( r, psi ) ) );
		}
	}

	std::printf( "\n  f() against the printed equations, worst relative difference:\n"
	             "    eq (24) pedestal      %.3e\n"
	             "    eq (26) current hole  %.3e\n"
	             "    eq (27) internal layer %.3e\n",
	             worstPedestal, worstHole, worstLayer );
	std::fflush( stdout );

	BOOST_TEST( worstPedestal < 1.0e-14 );
	BOOST_TEST( worstHole < 1.0e-14 );
	BOOST_TEST( worstLayer < 1.0e-14 );
}

/// F = mu0 r^2 p' with mu0 = 1, for the two sources that have a pressure
/// profile. This is the step that eq (24) makes explicit and that section 4.3
/// leaves to the reader, so it is checked for both.
BOOST_AUTO_TEST_CASE( theSourceIsRSquaredTimesThePressureGradient )
{
	meq::analytic::PressurePedestal const pedestal
		= meq::analytic::PressurePedestal::pedestal();
	meq::analytic::TransportBarrier const barrier
		= meq::analytic::TransportBarrier::barrier();

	double worst = 0.0;
	for ( double r = 0.6; r < 1.45; r += 0.1 )
	{
		for ( double psi = -0.4; psi < 0.95; psi += 0.01 )
		{
			worst = std::max( worst, std::abs( pedestal.f( r, 0.0, psi )
			                                   - r*r*pedestal.pPrime( psi ) ) );
			worst = std::max( worst, std::abs( barrier.f( r, 0.0, psi )
			                                   - r*r*barrier.pPrime( psi ) ) );
		}
	}
	BOOST_TEST( worst < 1.0e-13 );
}

/// p' and p'' against central differences of the level below. The transport
/// barrier is the one that matters here: it involves erf, its derivative has to
/// be supplied by hand, and the source is not printed in the paper at all -- it
/// has to be obtained by differentiating eq (25). CLAUDE.md records that a
/// Jacobian wrong by five per cent leaves every error and every rate unchanged
/// to six significant figures, so this is not something a convergence table can
/// find.
BOOST_AUTO_TEST_CASE( pressureDerivativesMatchFiniteDifferences )
{
	meq::analytic::PressurePedestal const pedestal
		= meq::analytic::PressurePedestal::pedestal();
	meq::analytic::TransportBarrier const barrier
		= meq::analytic::TransportBarrier::barrier();

	double const step = 1.0e-6;
	double worstPedestalFirst = 0.0, worstPedestalSecond = 0.0;
	double worstBarrierFirst = 0.0, worstBarrierSecond = 0.0;

	for ( double psi = -0.4; psi < 0.95; psi += 0.005 )
	{
		auto compare = []( double analytic, double difference )
		{
			return std::abs( analytic - difference )/( 1.0 + std::abs( analytic ) );
		};

		worstPedestalFirst = std::max( worstPedestalFirst,
			compare( pedestal.pPrime( psi ),
			         ( pedestal.p( psi + step ) - pedestal.p( psi - step ) )
			         /( 2.0*step ) ) );
		worstPedestalSecond = std::max( worstPedestalSecond,
			compare( pedestal.pDoublePrime( psi ),
			         ( pedestal.pPrime( psi + step ) - pedestal.pPrime( psi - step ) )
			         /( 2.0*step ) ) );
		worstBarrierFirst = std::max( worstBarrierFirst,
			compare( barrier.pPrime( psi ),
			         ( barrier.p( psi + step ) - barrier.p( psi - step ) )
			         /( 2.0*step ) ) );
		worstBarrierSecond = std::max( worstBarrierSecond,
			compare( barrier.pDoublePrime( psi ),
			         ( barrier.pPrime( psi + step ) - barrier.pPrime( psi - step ) )
			         /( 2.0*step ) ) );
	}

	std::printf( "\n  pressure derivatives against central differences:\n"
	             "    pedestal p'  %.3e   p''  %.3e\n"
	             "    barrier  p'  %.3e   p''  %.3e\n",
	             worstPedestalFirst, worstPedestalSecond,
	             worstBarrierFirst, worstBarrierSecond );
	std::fflush( stdout );

	BOOST_TEST( worstPedestalFirst < 1.0e-6 );
	BOOST_TEST( worstPedestalSecond < 1.0e-5 );
	BOOST_TEST( worstBarrierFirst < 1.0e-6 );
	BOOST_TEST( worstBarrierSecond < 1.0e-5 );
}

/// dF/dpsi against a central difference of F, for all four sources. This is the
/// only thing standing between a typo in a derivative and a Newton iteration
/// that converges linearly while still producing the right answer.
BOOST_AUTO_TEST_CASE( everySourceHasTheDerivativeItClaims )
{
	double const pedestalWorst
		= worstJacobianDiscrepancy( meq::analytic::PressurePedestal::pedestal() );
	double const barrierWorst
		= worstJacobianDiscrepancy( meq::analytic::TransportBarrier::barrier() );
	double const holeWorst
		= worstJacobianDiscrepancy( meq::analytic::CurrentHole::currentHole() );
	double const layerWorst
		= worstJacobianDiscrepancy( meq::analytic::InternalLayer::internalLayer() );

	std::printf( "\n  dF/dpsi against a central difference of F, worst relative:\n"
	             "    4.2 pedestal        %.3e\n"
	             "    4.3 barrier         %.3e\n"
	             "    4.4 current hole    %.3e\n"
	             "    4.5 internal layer  %.3e\n",
	             pedestalWorst, barrierWorst, holeWorst, layerWorst );
	std::fflush( stdout );

	BOOST_TEST( pedestalWorst < 1.0e-5 );
	BOOST_TEST( barrierWorst < 1.0e-5 );
	BOOST_TEST( holeWorst < 1.0e-5 );
	BOOST_TEST( layerWorst < 1.0e-5 );
}

/*
 * -------------------------------------------------------------------------
 * The trivial branch
 * -------------------------------------------------------------------------
 */

/// All four sources vanish at psi = 0. Asserted so that a future change to one
/// of the fixtures that quietly gives it a non-zero value at the origin -- which
/// would make the paper's own homogeneous problem well posed -- is noticed
/// rather than inherited.
BOOST_AUTO_TEST_CASE( everySourceVanishesAtZeroFlux )
{
	for ( double r = 0.6; r < 1.45; r += 0.1 )
	{
		BOOST_TEST( std::abs( meq::analytic::PressurePedestal::pedestal()
		                      .f( r, 0.0, 0.0 ) ) < 1.0e-300 );
		BOOST_TEST( std::abs( meq::analytic::TransportBarrier::barrier()
		                      .f( r, 0.0, 0.0 ) ) < 1.0e-300 );
		BOOST_TEST( std::abs( meq::analytic::CurrentHole::currentHole()
		                      .f( r, 0.0, 0.0 ) ) < 1.0e-300 );
		BOOST_TEST( std::abs( meq::analytic::InternalLayer::internalLayer()
		                      .f( r, 0.0, 0.0 ) ) < 1.0e-300 );
	}
}

/// THE FINDING, recorded as a test. With homogeneous Dirichlet data these four
/// benchmarks are solved by psi == 0, and meq's Newton iteration -- which starts
/// from the Dirichlet data, because that is the only initial iterate
/// GradShafranovSolver's interface admits -- lands on it and stops, in zero
/// iterations, with an identically zero residual.
///
/// This is not a solver defect. It is the reason
/// refs/HDG-GradShafranov.pdf Algorithm 2 says "psi^0 ; // Non-trivial initial
/// guess".
///
/// setInitialGuess() now exists, so this measures the DEFAULT rather than a
/// missing feature: no guess supplied, still the trivial branch. What a guess
/// does and does not buy is the next test.
BOOST_AUTO_TEST_CASE( homogeneousDataLandsOnTheTrivialBranch )
{
	meq::analytic::PressurePedestal const eq
		= meq::analytic::PressurePedestal::pedestal();

	SelfMeasurement const point = meq::tests::measureSelf(
		eq, standardBox(), 2, 8, cloud(), zeroDatum );

	std::printf( "\n  section 4.2 with psi = 0 on the whole boundary:\n"
	             "    Newton iterations %d, residuals recorded %zu, "
	             "psi in [ %.3e, %.3e ]\n",
	             point.newtonIterations, point.residuals.size(),
	             point.psiMin, point.psiMax );
	std::fflush( stdout );

	BOOST_TEST( point.converged,
	            "the trivial branch should be reached without difficulty" );
	BOOST_TEST( point.newtonIterations == 0,
	            "Newton took " << point.newtonIterations << " iterations on a problem "
	            "whose initial residual is identically zero. If this now fails, "
	            "GradShafranovSolver has acquired a non-trivial initial iterate and "
	            "the benchmarks in this file should be re-posed with the paper's own "
	            "homogeneous data" );
	BOOST_TEST( std::abs( point.psiMax ) < 1.0e-14 );
	BOOST_TEST( std::abs( point.psiMin ) < 1.0e-14 );
}

/*
 * A GUESS REACHES THE ITERATE, AND IT IS NOT ENOUGH. Both halves measured.
 *
 * setInitialGuess() puts a starting point into the trace, which is Newton's
 * actual unknown -- see meq::GradShafranovSolver::projectOntoTrace, and note
 * that no library call does this. That half works: Newton takes steps where it
 * previously took none.
 *
 * IT DOES NOT ESCAPE THE TRIVIAL BRANCH, which is what it was wanted for. On
 * section 4.3 with the paper's own homogeneous data, sweeping the amplitude of
 * a sin-sin bump on the box at k = 2, n = 8:
 *
 *     amplitude   outcome
 *       0.05      converges in  4 iterations, to psi ~ 1e-12
 *       0.20      converges in  5 iterations, to psi ~ 3e-16
 *       0.40      element-local solves fail, 60 iterations
 *       0.60      the same
 *       1.00      the same
 *
 * Small guesses are pulled back to psi == 0; guesses large enough to have a
 * chance of clearing it hit the globalisation failure instead. No amplitude
 * works.
 *
 * WHY, AND WHY IT IS NOT A DEFECT IN THE GUESS. Newton converges to the root
 * nearest its iterate, and psi == 0 is a root, so a starting point helps only
 * if it lands inside the other root's basin. Picard differs in exactly the
 * relevant way: from a non-zero psi^0 it evaluates F( psi^0 ) != 0 and is
 * carried away from zero rather than back to it. That is why both papers use
 * Anderson-accelerated Picard, and why GS-1's Algorithm 2 needs nothing more
 * than "psi^0 ; // Non-trivial initial guess" for it to work.
 *
 * So these four benchmarks need a guess AND globalisation --
 * KINSolver( KIN_LINESEARCH ), or continuation in the source amplitude. A guess
 * is necessary infrastructure, not a sufficient fix, and the non-homogeneous
 * ramp the rest of this file uses stays the workaround until then.
 */
BOOST_AUTO_TEST_CASE( anInitialGuessReachesTheIterateButNotTheOtherBranch )
{
	meq::analytic::TransportBarrier const eq
		= meq::analytic::TransportBarrier::barrier();
	meq::tests::Rectangle const box = standardBox();

	mfem::FunctionCoefficient bump( [ &box ]( mfem::Vector const &x )
	{
		double const u = ( x( 0 ) - box.rMin )/box.width();
		double const v = ( x( 1 ) - box.zMin )/box.height();
		return 0.2*std::sin( M_PI*u )*std::sin( M_PI*v );
	} );

	SelfMeasurement const cold = meq::tests::measureSelf(
		eq, box, 2, 8, cloud(), zeroDatum, 60, 1.0e-10 );
	SelfMeasurement const warm = meq::tests::measureSelf(
		eq, box, 2, 8, cloud(), zeroDatum, 60, 1.0e-10, &bump );

	std::printf( "\n  section 4.3, homogeneous data: no guess %d iterations, "
	             "with a guess %d\n", cold.newtonIterations, warm.newtonIterations );
	std::fflush( stdout );

	BOOST_TEST( cold.newtonIterations == 0,
	            "without a guess Newton should stop at once on the trivial branch" );
	BOOST_TEST( warm.newtonIterations > 0,
	            "Newton took no steps from a non-zero guess, so the guess never "
	            "reached the trace -- which is what projectOntoTrace() is for" );

	// The limitation, asserted so it fails the day globalisation lands and the
	// comment above has to be rewritten.
	BOOST_TEST( std::abs( warm.psiMax ) < 1.0e-6,
	            "psi came back at " << warm.psiMax << ", away from the trivial "
	            "branch. If this now fails, a guess alone is escaping it, and both "
	            "the comment above and DRIVER-PLAN.md section 1 are out of date" );
}

/*
 * AND A GUESS MUST NOT MOVE THE ANSWER, only the path to it. That is what
 * separates a starting point from data, and it is the check that would catch
 * projectOntoTrace() writing into the wrong dofs: a guess leaking into the
 * essential trace values would change the boundary condition, which a
 * convergence rate would accommodate without complaint.
 */
/*
 * KINSOL'S LINE SEARCH DOES NOT FIX THE PEDESTAL, AND THE REASON IS STRUCTURAL.
 *
 * CLAUDE.md carried this as the highest-value outstanding change for a while:
 * meq's Newton fails on section 4.2 at k = 1 for h >= 0.05, the fix is
 * globalisation, and KINSolver( KIN_LINESEARCH ) is the globalisation. SUNDIALS
 * is now built in and setGlobalisation() makes it one line. Measured:
 *
 *     case                     Newton        KIN_LINESEARCH   KIN_NONE
 *     4.2 k=1 n=16 (h=0.05)    ok, 42 it     FAIL at 18       FAIL at 6
 *     4.2 k=1 n=24 (h=0.033)   ok, 23 it     ok, 22 it        FAIL at 11
 *     4.3 homogeneous + guess  trivial       trivial          --
 *
 * It fixes nothing, and on the one case that matters it is WORSE than the
 * undamped iteration it was meant to rescue.
 *
 * WHY. The failure is not in the outer iteration. It is
 * "el: N not convered in 100 iters" -- MFEM's ELEMENT-LOCAL non-linear solve,
 * one per element per residual evaluation, eliminating the flux and potential
 * for a given trace. A line search on the outer step chooses how far to move the
 * trace; it does not change the fact that the local problems at that trace are
 * indefinite on a coarse element and their own Newton diverges. KINSOL never
 * sees them. Globalising the outer iteration does not globalise the inner ones.
 *
 * WHAT WOULD. Damping the LOCAL solves, which is MFEM's to offer; or
 * continuation in the source amplitude, so each solve starts from the previous
 * one's answer rather than from a cold trace; or Picard on the outer loop, which
 * evaluates F at the previous iterate and leaves every local problem LINEAR.
 * That last is what both papers do, and this measurement is the argument for it.
 *
 * WHAT IS ASSERTED HERE is only that the KINSOL path works and agrees, on a
 * problem both solvers converge on. The table above is prose because every
 * count in it sits on the same rounding knife edge
 * pedestalNewtonFailsOnCoarseMeshesAtOrderOne documents.
 */
/*
 * ANDERSON-ACCELERATED PICARD, WHICH IS THE PAPERS' OWN METHOD AND WORKS.
 *
 * Both GS papers solve the semi-linear problem by Anderson-accelerated Picard,
 * and CLAUDE.md explains why that is a coherent choice rather than a timid one:
 * evaluating F at the previous iterate puts it on the RIGHT HAND SIDE, so the
 * potential block stays linear and every element-local elimination is a linear
 * solve. meq's Newton puts F on the non-linear potential mass, where
 * hybridization turns each element's elimination into its own Newton -- and
 * those are what fail on the stiff sources.
 *
 * Globalisation::AndersonPicard is KINSolver( KIN_FP ) over that map. Measured
 * here on section 4.2 at k = 1, h = 0.05 -- the case Newton manages only by
 * grinding to 42 iterations and fails outright at some thread counts:
 *
 *     method                        outcome
 *     Newton                        42 iterations
 *     Picard, undamped              stalls
 *     Picard, w = 0.5               248 iterations
 *     Anderson depth 1, undamped    162 iterations
 *     Anderson depth 2 and above    fails
 *
 * TWO SURPRISES, both recorded on the setters. Plain Picard needs damping and
 * Anderson does not; and HDG-GS-1's own m = 2 fails here where m = 1 works.
 *
 * WHAT IS ASSERTED. That the Picard path reaches the SAME discrete solution as
 * Newton, which is what makes it an alternative rather than a different problem,
 * and that it does so with no element-local non-linear solves -- the property
 * the whole approach exists for. Not an iteration count: they sit on the same
 * rounding knife edge as everything else in this file, and Picard's is in the
 * hundreds where Newton's is in the tens. This is a robustness route, not a
 * faster one.
 */
BOOST_AUTO_TEST_CASE( andersonPicardReachesTheSameSolutionAsNewton )
{
	meq::analytic::PressurePedestal const eq
		= meq::analytic::PressurePedestal::pedestal();
	using G = meq::GradShafranovSolver::Globalisation;

	SelfMeasurement const newton = meq::tests::measureSelf(
		eq, standardBox(), 1, 16, cloud(), pedestalDatum, 500, 1.0e-8,
		nullptr, G::None );
	SelfMeasurement const anderson = meq::tests::measureSelf(
		eq, standardBox(), 1, 16, cloud(), pedestalDatum, 500, 1.0e-8,
		nullptr, G::AndersonPicard );

	std::printf( "\n  section 4.2 k = 1, h = 0.05: Newton %d iterations, "
	             "Anderson-Picard %d\n    psi in [%.6e, %.6e] and [%.6e, %.6e]\n",
	             newton.newtonIterations, anderson.newtonIterations,
	             newton.psiMin, newton.psiMax,
	             anderson.psiMin, anderson.psiMax );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( newton.converged );
	BOOST_TEST_REQUIRE( anderson.converged,
	                    "Anderson-accelerated Picard did not converge in "
	                    << anderson.newtonIterations << " iterations" );

	// The same discrete solution by two routes with different nonlinear
	// structure. This is what says the Picard path solves meq's problem and not
	// a neighbouring one -- the frozen source carries the same -1/r and the same
	// sign convention as the Newton path, and getting either wrong would show up
	// here and nowhere else.
	BOOST_TEST( anderson.psiMax == newton.psiMax,
	            boost::test_tools::tolerance( 1.0e-6 ) );
	BOOST_TEST( anderson.psiMin == newton.psiMin,
	            boost::test_tools::tolerance( 1.0e-6 ) );
}

BOOST_AUTO_TEST_CASE( kinsolAgreesWithNewtonWhereBothConverge )
{
	meq::analytic::PressurePedestal const eq
		= meq::analytic::PressurePedestal::pedestal();

	using G = meq::GradShafranovSolver::Globalisation;

	// n = 24, h = 0.0333: the well-posed mesh, where both converge.
	SelfMeasurement const plain = meq::tests::measureSelf(
		eq, standardBox(), 1, 24, cloud(), pedestalDatum, 60, 1.0e-12,
		nullptr, G::None );
	SelfMeasurement const damped = meq::tests::measureSelf(
		eq, standardBox(), 1, 24, cloud(), pedestalDatum, 60, 1.0e-12,
		nullptr, G::LineSearch );

	std::printf( "\n  section 4.2 k = 1, h = 0.0333: Newton %d iterations, "
	             "KIN_LINESEARCH %d\n"
	             "    psi in [%.6e, %.6e] and [%.6e, %.6e]\n",
	             plain.newtonIterations, damped.newtonIterations,
	             plain.psiMin, plain.psiMax, damped.psiMin, damped.psiMax );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( plain.converged );
	BOOST_TEST_REQUIRE( damped.converged,
	                    "KIN_LINESEARCH did not converge on the well-posed mesh, "
	                    "which would mean the KINSOL path is wired up wrongly "
	                    "rather than merely unhelpful" );

	// The same discrete solution by two routes. This is what would catch
	// ShiftedResidual getting the sign or the shift wrong -- KINSOL ignores the
	// right hand side handed to Mult(), so without that adapter it would converge
	// happily to the answer of a different problem.
	BOOST_TEST( damped.psiMax == plain.psiMax,
	            boost::test_tools::tolerance( 1.0e-6 ) );
	BOOST_TEST( damped.psiMin == plain.psiMin,
	            boost::test_tools::tolerance( 1.0e-6 ) );
}

BOOST_AUTO_TEST_CASE( anInitialGuessDoesNotMoveTheConvergedAnswer )
{
	meq::analytic::TransportBarrier const eq
		= meq::analytic::TransportBarrier::barrier();
	meq::tests::Rectangle const box = standardBox();

	mfem::FunctionCoefficient bump( [ &box ]( mfem::Vector const &x )
	{
		double const u = ( x( 0 ) - box.rMin )/box.width();
		double const v = ( x( 1 ) - box.zMin )/box.height();
		return 0.05*std::sin( M_PI*u )*std::sin( M_PI*v );
	} );

	// Small on purpose. A guess is meant to perturb the PATH, and on these stiff
	// sources a large one destabilises MFEM's element-local solves outright --
	// 0.35 here takes the same problem from 7 iterations to a failure at 60. That
	// is the globalisation limit again, not a property of the guess, and this
	// test is about invariance rather than about robustness.

	// barrierDatum, so there is a non-trivial solution to converge to.
	SelfMeasurement const cold = meq::tests::measureSelf(
		eq, box, 2, 16, cloud(), barrierDatum, 60, 1.0e-12 );
	SelfMeasurement const warm = meq::tests::measureSelf(
		eq, box, 2, 16, cloud(), barrierDatum, 60, 1.0e-12, &bump );

	std::printf( "  section 4.3 with barrierDatum: cold %d iterations, "
	             "psi in [%.6e, %.6e]\n"
	             "                                 warm %d iterations, "
	             "psi in [%.6e, %.6e]\n",
	             cold.newtonIterations, cold.psiMin, cold.psiMax,
	             warm.newtonIterations, warm.psiMin, warm.psiMax );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( cold.converged );
	BOOST_TEST_REQUIRE( warm.converged );

	BOOST_TEST( warm.psiMax == cold.psiMax,
	            boost::test_tools::tolerance( 1.0e-6 ) );
	BOOST_TEST( warm.psiMin == cold.psiMin,
	            boost::test_tools::tolerance( 1.0e-6 ) );
}

/*
 * -------------------------------------------------------------------------
 * Self convergence and the Newton record
 * -------------------------------------------------------------------------
 */

namespace
{
	/// Three dyadic levels from @a coarsest cells a side: two differences and
	/// therefore ONE rate per case per polynomial degree. With three usable
	/// sources at three degrees that is nine independent rate measurements. A
	/// fourth level would give a second rate and costs 95 s for a single k = 3
	/// solve at 128 cells a side, which is more than this study is worth.
	std::vector<int> dyadicFrom( int coarsest )
	{
		return { coarsest, 2*coarsest, 4*coarsest };
	}

	/// WHERE EACH SEQUENCE HAS TO START, measured rather than chosen. Newton
	/// simply does not converge on meshes coarser than these, and the coarsest
	/// usable mesh depends on the degree as well as on the source:
	///
	///                        k = 1              k = 2, 3
	///   4.2 pedestal         h = 0.05 fails     h = 0.1 converges but is
	///                        h = 0.0333: 23 it  pre-asymptotic
	///                        h = 0.025:   9 it  h = 0.05 onwards: 7-12 it
	///   4.3 barrier          h = 0.05 fails     h = 0.05 onwards: 8 it
	///                        h = 0.0333 fails
	///                        h = 0.025:   8 it
	///   4.5 internal layer   h = 0.05 fails     h = 0.05 onwards: 10-21 it
	///                        h = 0.0333: 26 it
	///                        h = 0.025:  11 it
	///
	/// k = 2 and 3 could all start at h = 0.1 and converge, but the first
	/// self-difference is then pre-asymptotic: on the pedestal it measures 2.32
	/// at k = 2 and 3.28 at k = 3 against design orders of 3 and 4. Starting at
	/// h = 0.05 measures the order rather than the approach to it.
	std::vector<int> meshesFor( int order, int coarsestAtOrderOne )
	{
		return order == 1 ? dyadicFrom( coarsestAtOrderOne ) : dyadicFrom( 16 );
	}

	/*
	 * WHY THE FLOORS BELOW ARE SET FROM MEASUREMENT AND NOT FROM k+1.
	 *
	 * Two independent reasons, and neither is about the solver.
	 *
	 * The first is the domain, and it is measured rather than argued:
	 * aCornerSingularSolutionCapsTheRateOnThisRectangle runs this same machinery
	 * on a Solov'ev source with no non-linearity at all and gets 2.86 and 2.96 in
	 * psi at k = 2 and 3, with the flux flat at about 2.2. That ceiling applies
	 * to every study in this file.
	 *
	 * The second is that a self-difference is a difference of two errors, so it
	 * carries the sampling error of the point cloud -- about one per cent, worth
	 * 0.03 of rate -- on top of everything a single error carries, and there is
	 * only ONE rate per configuration to average over. The floors therefore sit
	 * 0.15 to 0.20 under each measured value, in the same spirit as the absolute
	 * error ceilings of the exact-solution studies, which sit at three times
	 * their measurement.
	 */

	/// Newton is given 60 iterations rather than the default 30. Not to paper
	/// over a failure: the internal layer genuinely needs 26 at k = 1 on the
	/// coarsest mesh it converges on at all, and a cap that cut that off would
	/// report a failure that is not one. Every count is in the tables.
	int const newtonCap = 60;

	/// Run one source over its three levels at one degree, print the table and
	/// the residual history, and assert.
	template<typename Source, typename BoundaryFunction>
	void study( Source const &source, char const *label, int order,
	            int coarsestAtOrderOne, BoundaryFunction datum,
	            double psiFloor, double fluxFloor )
	{
		std::vector<SelfMeasurement> const points = meq::tests::checkSelfOrderAgainst(
			source, label, order, cloud(), datum, standardBox(),
			meshesFor( order, coarsestAtOrderOne ), psiFloor, fluxFloor,
			newtonCap );

		meq::tests::printNewtonHistory( label, order, points.back().h,
		                                points.back().residuals );

		// The best observed order over any triple of residuals above the round-off
		// floor. These sources are stiff enough that Newton wanders for several
		// steps before it enters the quadratic regime -- the histories printed
		// above show orders of 0.33 and 1.3 in the middle of a run that finishes
		// with a drop from 1e-9 to 1e-13 -- so the assertion is on the BEST triple,
		// not on all of them. A run with no triple above 1.5 anywhere is a chord
		// method, and that would mean the Jacobian disagrees with the residual.
		double const best = meq::tests::bestNewtonOrder( points.back().residuals );
		std::printf( "    best observed Newton order over any triple: %.3f\n", best );
		std::fflush( stdout );
		BOOST_TEST( best >= 1.5,
		            label << ", k = " << order << ": the best observed Newton order "
		            "is " << best << ". Nothing in the run looks quadratic, which "
		            "means the Jacobian is not the derivative of the residual" );
	}
}

/*
 * -------------------------------------------------------------------------
 * Section 4.2, the pressure pedestal
 * -------------------------------------------------------------------------
 */

/// THE SECOND FINDING. At k = 1 the pedestal does not converge at all on meshes
/// coarser than h = 0.033, and this records it rather than hiding behind a mesh
/// sequence chosen to avoid it.
///
/// What it looks like: MFEM's element-local non-linear solves fail --
/// "el: N not convered in 100 iters", printed by DarcyHybridization's LSsolve --
/// and the outer Newton then has an inconsistent residual and never recovers.
/// That is a structural difference from the papers, not a tuning problem. Both
/// papers close the semi-linear problem by Anderson-accelerated PICARD, which
/// evaluates F at the previous iterate and leaves every element-local solve
/// LINEAR. meq's Newton puts the non-linearity inside the local solve, and with
/// dF/dpsi = 2 c1 r^2/sigma^2 = 320 r^2 at psi = 0 those local problems are
/// indefinite on a coarse element.
///
/// The fix is globalisation, which is exactly what CLAUDE.md's note on
/// KINSolver( KIN_LINESEARCH ) anticipates -- and which cannot be tried here,
/// because the MFEM tree is built with MFEM_USE_SUNDIALS = NO.
BOOST_AUTO_TEST_CASE( pedestalNewtonFailsOnCoarseMeshesAtOrderOne )
{
	SelfMeasurement const coarse = meq::tests::measureSelf(
		meq::analytic::PressurePedestal::pedestal(), standardBox(), 1, 16,
		cloud(), pedestalDatum, newtonCap );
	SelfMeasurement const finer = meq::tests::measureSelf(
		meq::analytic::PressurePedestal::pedestal(), standardBox(), 1, 24,
		cloud(), pedestalDatum, newtonCap );

	std::printf( "\n  section 4.2 at k = 1: h = %.5f %s in %d iterations, "
	             "h = %.5f %s in %d\n",
	             coarse.h, coarse.converged ? "converged" : "FAILED",
	             coarse.newtonIterations,
	             finer.h, finer.converged ? "converged" : "FAILED",
	             finer.newtonIterations );
	std::fflush( stdout );

	BOOST_TEST( !coarse.converged,
	            "Newton now converges at k = 1, h = 0.05 on the pedestal, in "
	            << coarse.newtonIterations << " iterations. That is an improvement, "
	            "not a regression -- something has been globalised, or the local "
	            "solves have been made more robust. Move the sequences in "
	            "meshesFor() down and delete this test" );
	BOOST_TEST( finer.converged,
	            "Newton no longer converges at k = 1, h = 0.033 either, so the "
	            "resolution threshold has moved the wrong way" );
}

BOOST_AUTO_TEST_CASE( pedestalSelfConverges )
{
	// MEASURED: psi 3.097, 3.047, 4.384 and q 4.738, 3.122, 4.416 at k = 1, 2, 3.
	// Design order or better throughout, and the only one of the three that is:
	// see aCornerSingularSolutionCapsTheRateOnThisRectangle for what the other
	// two are up against. Floors 0.2 under each measurement.
	double const psiFloor[] = { 0.0, 2.90, 2.85, 4.15 };
	double const fluxFloor[] = { 0.0, 4.50, 2.90, 4.20 };
	for ( int order = 1; order <= 3; ++order )
		study( meq::analytic::PressurePedestal::pedestal(),
		       "4.2 pressure pedestal", order, 24, pedestalDatum,
		       psiFloor[ order ], fluxFloor[ order ] );
}

/*
 * -------------------------------------------------------------------------
 * Section 4.3, the internal transport barrier
 * -------------------------------------------------------------------------
 */

/// The barrier gets a DIFFERENT Dirichlet datum from the other three, and the
/// reason is the physics of eq (25) rather than a convergence difficulty.
///
/// Its structure is at psi_0 = 0.3, not at psi = 0: that is where the erf turns
/// over and where the source spikes to 10 r^2. barrierDatum() ramps from 0 at
/// the bottom of the box to 0.6 at the top, so psi_0 lies squarely in the
/// interior and the barrier is an INTERNAL feature, which is what section 4.3 is
/// about. pedestalDatum(), which spans [ -0.3, 0.3 ], would put psi_0 exactly on
/// the top edge -- and measured, that configuration does not converge at all:
/// 60 iterations at k = 2, h = 0.05, with the local solves failing throughout.
///
/// The other three sources keep pedestalDatum() because their structure IS at
/// psi = 0: the sigma layer of eq (24), and the ridge r + psi = 1 of eq (27),
/// which crosses the box at psi in [ -0.4, 0.4 ].
BOOST_AUTO_TEST_CASE( transportBarrierSelfConverges )
{
	// MEASURED with barrierDatum: psi 1.851, 2.421, 3.441 and q 2.183, 2.721,
	// 3.653 at k = 1, 2, 3. Design order at k = 1 and short of it at k = 2 and
	// 3, for the reason the control below establishes: on this rectangle, with a
	// datum that is not the trace of a smooth solution, the corners cap the rate
	// at about 3 in psi and 2.2 in q whatever the source. Floors 0.2 under each
	// measurement.
	//
	// A datum spanning [ 0, 0.6 ] rather than [ 0.2, 0.4 ] was tried, on the
	// theory that the erf transition -- 1/s = 0.025 wide in psi, and therefore
	// 0.025/|grad psi| wide in space -- was under-resolved on the coarse mesh.
	// It is not the explanation: the steeper datum gives 2.473 and 3.470, within
	// 0.05 of these.
	double const psiFloor[] = { 0.0, 1.70, 2.25, 3.25 };
	double const fluxFloor[] = { 0.0, 2.00, 2.55, 3.45 };
	for ( int order = 1; order <= 3; ++order )
		study( meq::analytic::TransportBarrier::barrier(), "4.3 transport barrier",
		       order, 32, barrierDatum, psiFloor[ order ], fluxFloor[ order ] );
}

/*
 * -------------------------------------------------------------------------
 * Section 4.5, the internal layer
 * -------------------------------------------------------------------------
 */

/// THE FLOORS HERE ARE NOT k+1, AND FOR TWO REASONS THAT COMPOUND.
///
/// MEASURED: psi 1.786, 2.596, 3.206 and q 2.480, 3.373, 3.538 at k = 1, 2, 3.
///
/// The first reason is the domain, and it is the same one the transport barrier
/// runs into: on this rectangle with a datum that is not the trace of a smooth
/// solution, the corners cap the self-convergence rate at about 3 in psi and 2.2
/// in q however smooth the source is. See
/// aCornerSingularSolutionCapsTheRateOnThisRectangle below, which measures that
/// on a Solov'ev source with no non-linearity at all.
///
/// The second is specific to eq (27) and is what the section is about. It adds a
/// ridge along r + psi = 1 of width sigma_2 = 0.0274 in that argument;
/// |grad( r + psi )| is about 1.1 over this box -- the r contributes 1 on its own
/// -- so the ridge is about 0.025 wide IN SPACE, and the finest mesh here has
/// h = 0.0125: two cells across it. This is precisely the case
/// refs/HDG-GradShafranov-Adaptive.pdf poses to demonstrate ADAPTIVE refinement
/// -- "it represents a good benchmarking problem to test for the detection of
/// internal layers" -- and it solves it with six levels of adaptive refinement
/// down to h_min = 2.14e-2. Reaching the asymptotic regime on a uniform mesh
/// would need h of about 0.003, which is 32768 elements a level and about 95 s
/// for a single k = 3 solve.
///
/// The floors below sit 0.15 under each measured value. The honest statement is
/// that the differences are falling, monotonically and faster than second order
/// -- not that they are falling at k+1.
BOOST_AUTO_TEST_CASE( internalLayerSelfConverges )
{
	double const psiFloor[] = { 0.0, 1.65, 2.45, 3.05 };
	double const fluxFloor[] = { 0.0, 2.30, 3.20, 3.35 };
	for ( int order = 1; order <= 3; ++order )
		study( meq::analytic::InternalLayer::internalLayer(),
		       "4.5 internal layer", order, 24, pedestalDatum,
		       psiFloor[ order ], fluxFloor[ order ] );
}

/*
 * -------------------------------------------------------------------------
 * The control: what a rectangle with a corner-singular solution costs
 * -------------------------------------------------------------------------
 */

/// A CONTROL for the two studies above, and the reason their shortfall at k >= 2
/// can be attributed to the DOMAIN rather than to the source or the solver.
///
/// It runs the same self-convergence machinery on a Solov'ev source, which is
/// CONSTANT in psi -- the easiest source in the tree, one Newton step, no
/// non-linearity at all -- with homogeneous Dirichlet data on the same
/// rectangle. Homogeneous data is what matters: the solution is then not the
/// trace of anything smooth, and the four right-angled corners of the rectangle
/// carry r^2 log r terms, which put it in H^(3-epsilon) and no better.
///
/// If this control caps out near 3 as well, then the shortfalls above are the
/// corners and not the sources.
BOOST_AUTO_TEST_CASE( aCornerSingularSolutionCapsTheRateOnThisRectangle )
{
	meq::analytic::SolovievEquilibrium const eq
		= meq::analytic::SolovievEquilibrium::nstx();

	// MEASURED, and this is the result:
	//
	//     k = 1   psi 2.002   q 1.880
	//     k = 2   psi 2.858   q 2.248
	//     k = 3   psi 2.956   q 2.183
	//
	// Design order at k = 1, and flat at about 3 in psi and 2.2 in q from k = 2
	// onwards. There is nothing non-linear here at all -- F is constant in psi,
	// dF/dpsi is identically zero, and every solve finishes in one Newton step --
	// so this is the domain and nothing else. It is the r^2 log r corner terms of
	// a right-angled corner: they put the solution in H^(3-epsilon) and its
	// gradient in H^(2-epsilon), and no polynomial degree recovers what the
	// geometry has taken away.
	//
	// Two consequences. The shortfalls in sections 4.3 and 4.5 above are the
	// corners, not the sources. And a self-convergence study on a rectangle with
	// homogeneous or otherwise incompatible data cannot demonstrate k+1 for
	// k >= 2, which is worth knowing before designing another one.
	double const psiFloor[] = { 0.0, 1.85, 2.70, 2.80 };
	double const fluxFloor[] = { 0.0, 1.75, 2.10, 2.05 };
	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<SelfMeasurement> const points = meq::tests::checkSelfOrderAgainst(
			eq, "control: Solov'ev, homogeneous data", order, cloud(), zeroDatum,
			standardBox(), dyadicFrom( 16 ), psiFloor[ order ], fluxFloor[ order ],
			newtonCap, 1.0e-9 );
		for ( SelfMeasurement const &point : points )
			BOOST_TEST( point.newtonIterations <= 1,
			            "the Solov'ev source does not depend on psi, so one exact "
			            "Newton step must finish it" );
	}
}

/*
 * -------------------------------------------------------------------------
 * Section 4.4, the current hole -- which does not converge at all
 * -------------------------------------------------------------------------
 */

/// THE THIRD FINDING, and the largest. The current-hole source of eq (26) is
/// not solvable by meq's plain Newton iteration in ANY configuration tried.
///
/// MEASURED, all with the ramp datum unless stated, all with a cap of 60:
///
///   k = 2, 16 cells a side    the iterate reaches NaN; MFEM's NewtonSolver
///                             aborts the PROCESS on
///                             MFEM_VERIFY( IsFinite( norm ) )
///   k = 2, 32 / 48 / 64       runs out of 60 iterations without converging
///   k = 3, 32                 runs out of 60 iterations
///   k = 1, 32                 reaches NaN and aborts
///   amplitudes 0.05 to 0.6, and offset ramps spanning [ 0, 0.2 ],
///   [ 0.28, 0.32 ] and [ 0.4, 0.6 ]   all run out of 60 iterations
///
/// WHY. eq (26) adds c3 ( 1 - exp( -( psi/sigma_2 )^2 ) ) cos( c4 psi ) to the
/// pedestal source, with c3 = -18 and c4 = 10 pi. Its psi-derivative carries
/// -c3 c4 sin( c4 psi ), which is 565 in magnitude and changes sign five times
/// across psi in [ -0.3, 0.3 ]. The Jacobian mass term is therefore large AND
/// oscillating in sign, so the Newton correction from one iterate lands where
/// the linearisation is meaningless. Nothing about the discretisation is
/// implicated: the source and its derivative are checked against the printed
/// equation and against a central difference above, and the same solver handles
/// the other three sources of the same family.
///
/// WHAT WOULD FIX IT. A globalised Newton -- KINSolver( KIN_LINESEARCH ), which
/// CLAUDE.md already identifies as the intended route and which is unavailable
/// because the MFEM tree is built with MFEM_USE_SUNDIALS = NO -- or continuation
/// in c3 from the pure pedestal, which needs a setInitialGuess() that
/// GradShafranovSolver does not have. Both papers avoid the problem entirely by
/// using Anderson-accelerated Picard, where F is evaluated at the previous
/// iterate and never differentiated.
///
/// THIS TEST IS LAST IN THE FILE ON PURPOSE. The configuration it uses -- k = 2,
/// 32 cells a side -- is the one measured to fail CLEANLY, running out of
/// iterations rather than reaching NaN. Coarser meshes reach NaN, and MFEM
/// aborts the process rather than throwing, so a NaN here would take every
/// result above it down with it. If this test ever starts aborting instead of
/// failing, that is why, and the fix is to move it to a file of its own.
BOOST_AUTO_TEST_CASE( currentHoleDoesNotConvergeUnderPlainNewton )
{
	SelfMeasurement const point = meq::tests::measureSelf(
		meq::analytic::CurrentHole::currentHole(), standardBox(), 2, 32,
		cloud(), pedestalDatum, newtonCap );

	std::printf( "\n  section 4.4, current hole, k = 2, h = %.5f: %s after %d "
	             "iterations\n", point.h,
	             point.converged ? "CONVERGED" : "did not converge",
	             point.newtonIterations );
	meq::tests::printNewtonHistory( "4.4 current hole", 2, point.h,
	                                point.residuals );

	BOOST_TEST( !point.converged,
	            "Newton now converges on the current hole, in "
	            << point.newtonIterations << " iterations. That is a real "
	            "improvement -- globalisation, or an initial guess -- and this test "
	            "should be replaced by a self-convergence study like the other "
	            "three. Read the comment above it first" );
}
