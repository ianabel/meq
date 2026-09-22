/*
 * XP-2: psi_bnd FROM THE LOCATED X-POINT, AS AN OUTER FIXED POINT.
 * FREE-BOUNDARY-PLAN.md section 10.6.
 *
 * FB-3's setBoundaryFluxPoint() pins psi_bnd = psi_h at a PRESCRIBED point,
 * which is right for a limiter -- the contact is a piece of hardware -- and
 * wrong for a divertor, whose X-point is a functional of the solution and moves
 * as Newton moves. XP-3 makes the X-point three more unknowns of the same
 * Newton; XP-2 is the halfway house that needs no new border at all:
 *
 *     solve with psi_bnd pinned at the CURRENT estimate of the X-point
 *     locate the saddle of the SOLVED q_h
 *     re-pin, re-solve
 *
 * an outer fixed point over ( R_X, z_X ), warm started from the previous sweep.
 *
 * THE DISCRETE CHOICE IS FROZEN WITHIN A SOLVE AND RE-DECIDED BETWEEN THEM,
 * which is section 10.5's answer to the one thing this stage cannot
 * differentiate: psi_bnd is a min over candidates and a min is not smooth, so
 * Newton is never asked to differentiate through it. Each solve is a smooth
 * problem with the bounding point FIXED.
 *
 *
 * WHY THIS RUNS ON THE MACHINE CASE AND NOT ON THE HALF-DISC FIXTURE, WHICH IS
 * A MEASURED FALSIFICATION AND NOT A PREFERENCE.
 *
 * FreeBoundaryCoupling.cpp's theTwoBorderSolveReportsATrueMagneticAxis is the
 * obvious host: it is already a physical free-boundary equilibrium with coils,
 * a prescribed current and ConfineToPlasma. It cannot be made diverted, and the
 * obstruction is structural rather than a matter of trying harder.
 *
 * psi vanishes IDENTICALLY on R = 0 -- the symmetry axis is fitted Dirichlet
 * boundary -- so the axis sits at normalised flux -psi_bnd/span. A NEGATIVE
 * psi_bnd therefore puts it at POSITIVE Psi, inside the plasma support, where
 * ConfineToPlasma does not switch gg' off; and F/R is mu_0 j_phi, so that is an
 * infinite toroidal current density on the axis. What the field then grows is
 * section 11.3's axis layer, and what CriticalPointFinder reports is the layer.
 * Measured on one such run: THIRTY critical points, every one at R < 0.1, a
 * ladder of maxima marching up the symmetry axis from z = -0.11 to z = +1.44,
 * and the "magnetic axis" at R = 0.022.
 *
 * That fixture reaches only psi_bnd = +9.5e-04 against psi_ax = 1.19e-02, a
 * margin of 8%, while a Shafranov vertical field alone contributes about
 * psi = B_v R^2/2 = -1.8e-02 at the limiter. The plasma's own flux outweighs
 * the coils' by a few per cent, so any conductor strong enough to make a null
 * tips psi_bnd negative first. Twelve configurations were run across three
 * routes -- an exterior filament swept in current, the same with the
 * equilibrium pair re-derived to hold the total vertical field, and a designed
 * three-conductor single null with the divertor meshed inside Omega -- and
 * every one of them ended at negative psi_bnd, at no off-axis saddle, or at a
 * saddle out in the vacuum ABOVE the plasma. Superposition predicted X-points
 * in three of six cases where the solve produced none, which is the other half
 * of it: the screen drops the plasma's response and the response is most of the
 * answer.
 *
 * examples/diverted-tokamak.toml IS THE MACHINE THIS RUNS ON INSTEAD, AND IT IS
 * BORROWED RATHER THAN DESIGNED. It is freegs4e's own A_testtokamak_classic:
 * an up-down asymmetric double null whose coil currents were SOLVED FOR by its
 * control system, which section 7.15 argues is what makes a free-boundary
 * Newton converge at all. Its psi_bndry IS the lower X-point's flux to every
 * digit and its psi_bnd/psi_ax is 0.39.
 *
 * A run of hand-designed single nulls on the LIMITED machine's geometry came
 * first and every one of them failed the same way -- negative psi_bnd and the
 * axis layer -- so the design machinery that used to live in this file is gone.
 * The finding it produced is kept: every design solve pinned psi_bnd at the
 * LIMITER, and a diverted machine's limiter is OUTSIDE its separatrix, so the
 * solve was being asked the wrong question rather than answering it badly.
 *
 *
 * THE INPUTS WERE CHECKED AND THE INNER SOLVE WAS THE DEFECT. WHAT FOLLOWS IS
 * THE DIAGNOSIS, AND THE TWO THINGS IT TOOK TO CLOSE IT.
 *
 * Both halves of the problem statement are verified independently of MEQ:
 *
 *   * the two profile tables reproduce the reference's own plasma current when
 *     integrated over the reference's own core -- 1.999667e+05 A against
 *     2.0e+05, so the span conversion is right and a healthy solve must report
 *     a profile scale of 1;
 *   * the Green's-function guess reproduces psi at the reference's magnetic
 *     axis to 4.1e-04 relative ( 8.268402e-02 against 8.271751e-02 ) and at its
 *     active X-point to 4.8e-04 ( 3.238847e-02 against 3.240413e-02 ), so the
 *     iteration starts on the physical branch.
 *
 * It does not stay there. Measured, one key changed at a time:
 *
 *   pin      k  guess    outcome
 *   X-point  2  33^2     99 steps to || R ||/|| r_0 || = 4.1e-11, psi_ax
 *                        8.532768e-02 -- and psi_h AT the reference axis reads
 *                        1.32e-02 against the reference's 8.27e-02. The located
 *                        axis is ( 2.2626, 0.5413 ), hard against Gamma
 *   X-point  3  33^2     14 steps, axis at ( 1.4691, 1.4401 )
 *   X-point  2  129^2    the residual falls to 3.4e-04 and STALLS there for 150
 *                        steps -- no runaway and no convergence -- with the
 *                        axis already at ( 2.2590, 0.3105 ) by step 13
 *   (1.68,0) 2  33^2     a LIMIT CYCLE: || R || descends to 7.1e-03 and jumps
 *                        back to 2.7e-02, six times over the 200-step cap. That
 *                        pin is inside the separatrix, where the reference
 *                        carries psi = 4.52e-02, so it is the easy question
 *
 * AND THE TRACE SAYS WHERE IT BREAKS. Printing the located axis at every
 * evaluation of the first run: the first fourteen are right -- ( 1.3636,
 * 0.0008 ) carrying psi = 8.077e-02, with the span climbing 4.14e-02, 6.20e-02,
 * 7.24e-02, ... , 8.268e-02 towards the reference's 8.2718e-02 -- and then
 * psi_ax collapses, 7.24e-02, 3.78e-02, 2.28e-02, with the axis following it
 * out to ( 2.2042, 0.2313 ). The solve is nearly there and then leaves.
 *
 * What that ruled out: the profile conversion, the guess, the polynomial
 * degree, and the axis competition -- a conductor's O-point is excluded on both
 * paths now ( meq::Source::conductors ) and the runaway axis is in no coil.
 * What was left is the pair section 10.5 names: a plasma support that MOVES
 * within the Newton, and a psi_bnd pinned exactly at a saddle, where the
 * enclosed area's derivative is not bounded.
 *
 *
 * IT TOOK TWO THINGS, NEITHER OF THEM SUFFICIENT ALONE, AND THE FIRST WAS A
 * DEFECT IN THIS FIXTURE.
 *
 * setPlasmaSupport() was called on the PLASMA source before the coil wrapper
 * existed, so the wrapper's own flag stayed false -- and the wrapper is what
 * the solver holds, so GradShafranovSolver::plasmaComponentWanted() read false
 * and XP-1's flood fill never ran. F was still confined, the inner source's
 * pointwise test being live, so nothing failed loudly. Measured once it was
 * fixed: the pointwise test carries 752 elements over THREE components and the
 * fill keeps the 419 holding the axis. The other 333 -- 44% of the candidates
 * -- are the private flux region, and they were getting a current channel
 * nobody asked for on the only diverted case in this tree. Section 10.3
 * predicted exactly that and XP-1 was built for it.
 *
 * The second is section 10.5's own prescription applied to the support as well
 * as to the bounding point: hold the edge insidePlasma() tests against fixed
 * for the whole of one solve and move it between solves --
 * meq::NormalisedSource::freezePlasmaEdge, with
 * GradShafranovSolver::setPlasmaSupportFrozen holding the fill's own mask the
 * same way. psi_bnd is pinned at psi_h at a SADDLE, so it is a non-local
 * functional of the iterate; letting the support edge follow it sweeps F's
 * on/off region across the plasma edge every residual, and the Jacobian carries
 * no surface term for that.
 *
 * A FULL CROSS AND NOT A SWEEP OF ONE KEY, WHICH IS THE POINT:
 *
 *     fill   frozen edge   outcome
 *     off    off           FAILED at the 200 cap, || R || stalled at 4.8e-04
 *     off    on            converged in 17, and to the WRONG BRANCH:
 *                          psi_ax -7.99e-02 against the reference's +8.27e-02
 *     on     off           FAILED at the 200 cap
 *     on     on            CONVERGED in 7, and XP-2 is met
 *
 * So neither is sufficient and both are necessary -- and every experiment that
 * changes ONE of them from the configuration this case shipped with lands in a
 * failing row, which is why a long run of one-key experiments narrowed the
 * suspects without ever reaching the cause. MEASUREMENTS.md M-82 has the
 * numbers.
 *
 * AND THE OLD "SUPPORT OFF CONVERGES" MEASUREMENT WAS TRUE AND ITS CONCLUSION
 * WAS WRONG. Turning ConfineToPlasma off did make the solve converge -- in 21
 * steps, to psi_ax = -8.09e-02, which is a DIFFERENT equilibrium with F switched
 * on in the vacuum. Reading that as "the moving support is the hazard" was an
 * inference from a real measurement to the wrong cause, because a third variable
 * -- the fill -- was silently off in both of its rows. It is the same shape as
 * this tree's other instrument-not-answer findings, and it is worth the space:
 * a one-key experiment only separates two hypotheses if everything else is
 * where you think it is.
 *
 * XP-2 IS MET. The outer fixed point contracts quadratically -- 1.909e-03,
 * 3.732e-06, 9.554e-10, 2.255e-13 over four sweeps, the inner solve falling to
 * two Newton steps as it does -- and it finds BOTH nulls of the double-null
 * machine, the active one at ( 1.093, -0.604 ) carrying psi 3.2379e-02 and the
 * upper at ( 1.109, 0.796 ) carrying 2.8931e-02, which is what
 * examples/diverted-tokamak.toml records freegs4e reporting for each. The
 * located X-point sits 4.378e-04 m from freegs4e's, psi_bnd agrees to 0.08% and
 * psi_ax to 0.07%, and psi_h AT the reference axis reads 8.265989e-02 against
 * this solve's own psi_ax of 8.266004e-02 -- the border closed on the field it
 * is a constraint on.
 *
 * The loop was never the thing that was wrong: it is the eight lines below, and
 * nothing in the diagnosis above is about it. What it lacked was an equilibrium
 * to iterate on.
 */
#define BOOST_TEST_MODULE XPointOuter
#include <boost/test/unit_test.hpp>

#include <cstdio>

// The machine, the reference values, the solve helper, the saddle filter and
// the support freeze. Shared with XPointBorder.cpp, which is XP-3 and whose
// acceptance is agreement with this case -- so the two must solve the same
// fixture and not a copy of it.
#include "convergence/DivertedMachine.hpp"

using namespace meqtest;

BOOST_AUTO_TEST_CASE( theBoundaryFluxConvergesToTheLocatedXPoint )
{
	/*
	 * THE BOOTSTRAP PIN IS DELIBERATELY NOT THE ANSWER, and 0.05 m is the
	 * distance a machine's own drawings would get you to: it is 8% of the
	 * null's height below the midplane and about two thirds of an element on
	 * this mesh, so the loop has somewhere to travel and the seeded search
	 * still reaches from one sweep to the next.
	 *
	 * Pinning AT the reference X-point would make the first sweep's step the
	 * only thing measured, and it would measure zero.
	 */
	double const bootstrapR = referenceXPointR + 0.05;
	double const bootstrapZ = referenceXPointZ + 0.05;

	Machine m = buildMachine();
	completeMachine( m );

	/*
	 * THE SUPPORT IS FROZEN BEFORE THE FIRST SOLVE AND RE-DECIDED IN THE LOOP,
	 * WHICH IS THE ONLY THING THAT DISTINGUISHES THIS FROM THE RUN M-82
	 * MEASURED FAILING.
	 *
	 * The bootstrap has no solved field to freeze at, so it freezes at the
	 * problem statement instead: the file's own PsiAxis as the axis and the
	 * transferred guess evaluated AT THE PIN as the edge, which is exactly what
	 * psi_bnd will mean once the border closes. The guess is freegs4e's
	 * equilibrium, so this is the reference's own topology and not a
	 * construction of this test.
	 */
	double bootstrapEdge = 0.0;
	BOOST_TEST_REQUIRE( valueAt( *m.carried, bootstrapR, bootstrapZ,
	                             bootstrapEdge ),
	                    "the transferred guess cannot be evaluated at the "
	                    "bootstrap pin, so there is no edge to freeze at" );
	freezeSupportAt( m, *m.carried,
	                 m.config->getSource().psiAxisGuess(), bootstrapEdge );
	std::printf( "\n  XP-2: the support is FROZEN within each solve and moved "
	             "between them\n    bootstrap edge %.6e at the pin, psi_ax guess "
	             "%.6e, %d of %d elements carry plasma over %d component(s)\n",
	             bootstrapEdge, m.config->getSource().psiAxisGuess(),
	             m.solver->plasmaComponentElements(),
	             m.solver->plasmaCandidateElements(),
	             m.solver->plasmaComponentCount() );
	std::fflush( stdout );

	// THE BOOTSTRAP: one solve pinned at the prior, and the saddle of its own
	// solved q_h.
	BOOST_TEST_REQUIRE( solveAt( m, bootstrapR, bootstrapZ ),
	                    "the diverted machine did not converge pinned at ( "
	                    << bootstrapR << ", " << bootstrapZ << " ) with the "
	                    "support frozen, so the outer loop has no starting "
	                    "point. M-82 measured the MOVING support as the cause "
	                    "of this failing; if it still fails frozen, that "
	                    "separation is wrong or the freeze does not cover the "
	                    "whole support" );

	meq::CriticalPoint x;
	{
		meq::CriticalPointFinder finder( *m.solver );
		BOOST_TEST_REQUIRE( findXPoint( finder, x ),
			"the limiter solve of the diverted machine carries no off-axis "
			"saddle, so this configuration is not diverted and XP-2 has nothing "
			"to iterate on" );
	}

	std::printf( "\n  XP-2: psi_bnd FROM THE LOCATED X-POINT, bootstrapped at "
	             "( %.3f, %.3f )\n", bootstrapR, bootstrapZ );
	std::printf( "    %-6s %5s %19s %11s %13s %13s %11s %7s %8s\n",
	             "sweep", "its", "X-point", "step", "psi_ax", "psi_bnd",
	             "psi_X - psi_h", "depth", "route" );
	std::printf( "    %-6s %5zu   (%6.3f,%7.3f) %11s %13.6e %13.6e %11s %7.3f "
	             "%8s\n", "boot", m.solver->newtonResiduals().size() - 1, x.radius,
	             x.z, "-", m.solver->psiAxis(), m.solver->psiBoundary(), "-",
	             elementDepth( x ), "SWEEP" );
	std::fflush( stdout );

	struct Sweep
	{
		double radius = 0.0;
		double z = 0.0;
		double step = 0.0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		double pinResidual = 0.0;
		double depth = 0.0;
		bool seeded = false;
		std::size_t iterations = 0;
	};
	std::vector<Sweep> sweeps;
	int seededSweeps = 0;
	// Section 10.5's own warning is that a discrete outer state can alternate
	// rather than settle, so how often the support changes size is counted and
	// printed rather than assumed away.
	int supportMoved = 0;

	// THE WARM START LIVES ACROSS THE WHOLE LOOP, the guess being BORROWED and
	// having to outlive the solve it seeds. Seeded from the bootstrap.
	mfem::GridFunction previous( m.solver->potential() );

	bool broke = false;
	for ( int outer = 0; outer < 8; ++outer )
	{
		double const pinnedR = x.radius;
		double const pinnedZ = x.z;
		// THE SUPPORT MOVES HERE AND NOWHERE ELSE, at the answer the previous
		// sweep reached -- the same outer fixed point the X-point itself is on,
		// and re-decided in the same place so the two cannot disagree about
		// which iterate they belong to.
		int const supportBefore = m.solver->plasmaComponentElements();
		freezeSupportAt( m, previous, m.solver->psiAxis(),
		                 m.solver->psiBoundary() );
		supportMoved += ( m.solver->plasmaComponentElements() != supportBefore );
		if ( !solveAt( m, pinnedR, pinnedZ, &previous ) )
		{
			broke = true;
			std::printf( "    %-6d %5s %19s\n", outer + 1, "-", "NO SOLVE" );
			break;
		}

		/*
		 * THE SEEDED SADDLE SEARCH, WHICH IS WHAT IT WAS BUILT FOR. The
		 * previous sweep's X-point is a prior worth having: the point moves by
		 * a fraction of an element, and CriticalPointFinder::sweep() costs one
		 * Newton per element where the seeded search costs the rings around
		 * one -- 32x to 205x on the XP-0 fixture, widening with refinement.
		 *
		 * ITS REACH IS ABOUT ONE AND A HALF ELEMENTS, so a step larger than
		 * that DECLINES rather than returning a wrong point, and the sweep is
		 * the fallback. How often that fires is a measurement rather than an
		 * implementation detail, so it is counted and printed.
		 */
		meq::CriticalPointFinder finder( *m.solver );
		meq::CriticalPoint next;
		bool const seeded = finder.tryFindCriticalPointFrom(
			pinnedR, pinnedZ, meq::AxisSense::Saddle, next );
		if ( seeded )
			++seededSweeps;
		if ( !seeded && !findXPoint( finder, next ) )
		{
			broke = true;
			std::printf( "    %-6d %5zu %19s\n", outer + 1,
			             m.solver->newtonResiduals().size() - 1, "NO SADDLE" );
			break;
		}

		Sweep sweep;
		sweep.radius = next.radius;
		sweep.z = next.z;
		sweep.step = std::hypot( next.radius - pinnedR, next.z - pinnedZ );
		sweep.psiAxis = m.solver->psiAxis();
		sweep.psiBoundary = m.solver->psiBoundary();
		// THE FIXED POINT'S OWN DEFINING PROPERTY. The border makes psi_bnd
		// equal psi_h at the PINNED point exactly, whatever that point is; what
		// is not automatic is that the pinned point IS the saddle.
		sweep.pinResidual = next.psi - m.solver->psiBoundary();
		sweep.depth = elementDepth( next );
		sweep.seeded = seeded;
		sweep.iterations = m.solver->newtonResiduals().size() - 1;
		sweeps.push_back( sweep );

		std::printf( "    %-6d %5zu   (%6.3f,%7.3f) %11.3e %13.6e %13.6e "
		             "%11.2e %7.3f %8s\n", outer + 1, sweep.iterations, sweep.radius,
		             sweep.z, sweep.step, sweep.psiAxis, sweep.psiBoundary,
		             sweep.pinResidual, sweep.depth,
		             sweep.seeded ? "seeded" : "SWEEP" );
		std::fflush( stdout );

		x = next;
		previous = m.solver->potential();
		if ( sweep.step < 1.0e-10 )
			break;
	}

	BOOST_TEST_REQUIRE( !broke, "the outer loop broke down" );
	BOOST_TEST_REQUIRE( sweeps.size() >= 3,
	                    "the outer loop ran " << sweeps.size() << " sweeps, "
	                    "too few to say anything about contraction" );

	Sweep const &last = sweeps.back();
	double const span = last.psiAxis - last.psiBoundary;

	/*
	 * THE PRECONDITION ON THE FIXTURE, CHECKED RATHER THAN TRUSTED. An X-point
	 * landing on a mesh line is held by two elements at once, each with its own
	 * root and both strictly contained, and which one is reported can flip for
	 * a sub-h^(k+1) move of the field -- so the border row would be built from a
	 * different element between sweeps and the loop would chatter rather than
	 * converge. elementDepth() says at length what was measured.
	 */
	BOOST_TEST( last.depth > 0.05,
		"the converged X-point sits " << last.depth << " into its element in "
		"reference coordinates, so it is on a mesh line and two elements hold "
		"it. That is a statement about where this machine's null falls on this "
		"mesh, not about the solver: move the designed target." );

	// AND THE SEEDED SEARCH IS WHAT THE LOOP RUNS ON, worth asserting rather
	// than only printing -- a loop silently falling back to a full sweep every
	// time is a different cost and would go unnoticed.
	std::printf( "    %d of %zu sweeps reached the X-point from the previous "
	             "one\n", seededSweeps, sweeps.size() );

	/*
	 * AND THE SUPPORT IS THE OTHER DISCRETE STATE ON THIS FIXED POINT, so it
	 * gets the same treatment as the X-point: reported, and asserted to settle.
	 *
	 * Section 10.5 warns that a combinatorial outer iteration can alternate
	 * rather than converge -- two element sets swapping on alternate sweeps --
	 * and says there is no globalisation for one. There is no cure to assert
	 * here, only the symptom: a support still changing size on the LAST sweep
	 * has not settled, whatever the X-point did, and a psi_bnd read off it is a
	 * statement about which of two states the loop stopped on.
	 */
	std::printf( "    the frozen support changed size on %d of %zu sweeps, and "
	             "holds %d of %d candidate elements at the end\n", supportMoved,
	             sweeps.size(), m.solver->plasmaComponentElements(),
	             m.solver->plasmaCandidateElements() );
	BOOST_TEST( supportMoved < static_cast<int>( sweeps.size() ),
		"the frozen support changed size on every one of " << sweeps.size()
		<< " sweeps, so the outer iteration has not settled on a topology -- "
		"section 10.5's alternating discrete state, which has no globalisation. "
		"The X-point's own convergence below says nothing about this: they are "
		"two states on one fixed point and either can chatter while the other "
		"contracts." );
	BOOST_TEST( seededSweeps + 1 >= static_cast<int>( sweeps.size() ),
		"only " << seededSweeps << " of " << sweeps.size() << " sweeps reached "
		"the X-point from the previous one. The seeded search reaches about one "
		"and a half elements, so repeated fallback is a statement about the "
		"outer map rather than about the search." );

	// ONE: IT IS A FIXED POINT, and it stops well inside an element -- a loop
	// merely wandering within one would look converged on a coarser measure.
	BOOST_TEST( last.step < 1.0e-8,
		"the outer iteration left the X-point moving by " << last.step
		<< " metres on its last sweep, so psi_bnd is whatever the loop happened "
		"to stop at rather than the flux at the saddle." );

	// TWO: IT CONTRACTS.
	BOOST_TEST( last.step < 1.0e-3*sweeps.front().step,
		"the outer step went from " << sweeps.front().step << " to "
		<< last.step << " over " << sweeps.size() << " sweeps. XP-2 is a "
		"fixed-point iteration and a slow one is a different algorithm from a "
		"fast one -- section 10.4's three-row border is the answer if this is "
		"what the map looks like." );

	// THREE: psi_bnd IS THE FLUX AT THE SADDLE.
	BOOST_TEST( std::abs( last.pinResidual ) < 1.0e-10*span,
		"psi at the located saddle differs from the psi_bnd the solve returned "
		"by " << last.pinResidual << ", which is "
		<< std::abs( last.pinResidual )/span << " of the span. The loop has "
		"converged to a point that is not the saddle." );

	/*
	 * FOUR: THE PLASMA IS DIVERTED AND THE NULL XP-2 FOUND IS THE ACTIVE ONE.
	 *
	 * This machine has TWO off-axis saddles -- freegs4e puts the lower at
	 * psi = 3.240412550738516e-02 and the upper at 2.891018784455272e-02, and
	 * that 3.5e-03 gap is the whole of what makes it a SINGLE null. The
	 * bounding surface is the null of LARGEST psi, psi having its maximum on
	 * the axis, so a loop that converged on the upper one would satisfy every
	 * assertion above and describe a different plasma.
	 *
	 * So the case sweeps the converged field for every off-axis saddle and
	 * asserts both halves: that there is more than one, which is the fixture's
	 * own premise, and that the one the loop stopped at carries the largest
	 * psi.
	 */
	std::vector<meq::CriticalPoint> nulls;
	{
		meq::CriticalPointFinder finder( *m.solver );
		for ( meq::CriticalPoint const &p : finder.sweep() )
			if ( p.type == meq::CriticalPointType::Saddle && p.radius > 0.30 )
				nulls.push_back( p );
	}
	std::printf( "    %zu off-axis saddles in the converged field:", nulls.size() );
	for ( meq::CriticalPoint const &p : nulls )
		std::printf( "  ( %6.3f, %7.3f ) psi %11.4e", p.radius, p.z, p.psi );
	std::printf( "\n" );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( nulls.size() >= 2,
		"the converged field carries " << nulls.size() << " off-axis saddle( s ). "
		"This machine is an up-down asymmetric DOUBLE null and the reference "
		"finds both, so one is a statement about the discretisation or about "
		"which branch the solve reached, and it makes the active-null check "
		"below vacuous." );

	double largest = -std::numeric_limits<double>::infinity();
	for ( meq::CriticalPoint const &p : nulls )
		largest = std::max( largest, p.psi );
	BOOST_TEST( std::abs( largest - last.psiBoundary ) < 1.0e-8*span,
		"the loop settled on a null carrying psi = " << last.psiBoundary
		<< " where the field's largest off-axis saddle carries " << largest
		<< ". The bounding surface is a MIN over candidates -- the null met "
		"first walking out from the axis -- and this one is not it." );

	/*
	 * AND AGAINST THE INDEPENDENT CODE, WHICH IS THE ONLY ASSERTION HERE THAT
	 * THE ANSWER IS RIGHT RATHER THAN SELF-CONSISTENT.
	 *
	 * The tolerance is per cent and the reason is recorded in
	 * examples/diverted-tokamak.toml: fgsref.py splines its own analytic
	 * profile before solving and the fit MOVES it -- 2.514e-05 of the amplitude
	 * in p' and 1.707e-02 in ff' -- while MEQ's tables are the analytic shape.
	 * So this is a cross-check at the level the two problem statements agree,
	 * not a convergence measurement.
	 */
	double const nullGap = std::hypot( last.radius - referenceXPointR,
	                                   last.z - referenceXPointZ );
	std::printf( "    against freegs4e: X-point ( %.6f, %.6f ) vs ( %.6f, "
	             "%.6f ), %.3e m apart;  psi_bnd %.6e vs %.6e;  psi_ax %.6e vs "
	             "%.6e\n", last.radius, last.z, referenceXPointR, referenceXPointZ,
	             nullGap, last.psiBoundary, referencePsiBoundary, last.psiAxis,
	             referencePsiAxis );
	std::fflush( stdout );

	BOOST_TEST( nullGap < 2.0e-2,
		"the located X-point is " << nullGap << " m from freegs4e's ( "
		<< referenceXPointR << ", " << referenceXPointZ << " )." );
	BOOST_TEST( std::abs( last.psiBoundary - referencePsiBoundary )
	            < 2.0e-2*std::abs( referencePsiBoundary ),
		"psi_bnd is " << last.psiBoundary << " against freegs4e's "
		<< referencePsiBoundary );
	BOOST_TEST( std::abs( last.psiAxis - referencePsiAxis )
	            < 2.0e-2*std::abs( referencePsiAxis ),
		"psi_ax is " << last.psiAxis << " against freegs4e's "
		<< referencePsiAxis );

	/*
	 * AND THE CORE IS WHERE THE REFERENCE PUTS IT, WHICH IS A DIFFERENT CLAIM
	 * FROM ANY OF THE ABOVE AND IS THE ONE THAT CATCHES A WRONG BRANCH.
	 *
	 * psi_ax, psi_bnd and I_p are all border unknowns, so a solve that has run
	 * away to a different equilibrium still reports three plausible numbers and
	 * a converged residual. What it cannot fake is psi_h AT the reference's own
	 * magnetic axis: on the branch this machine reached before XP-2 -- axis
	 * driven out to ( 2.2626, 0.5413 ), against Gamma -- psi at ( 1.3513,
	 * 0.0622 ) read 1.32e-02 where the reference carries 8.27e-02, while
	 * psi_ax itself read 8.53e-02 and looked healthy.
	 */
	double corePsi = 0.0;
	BOOST_TEST_REQUIRE( potentialAt( *m.solver, referenceAxisR, referenceAxisZ,
	                                 corePsi ),
	                    "the reference's magnetic axis is not in the mesh" );
	std::printf( "    psi_h at the reference axis ( %.4f, %.4f ) = %.6e "
	             "against psi_ax = %.6e\n", referenceAxisR, referenceAxisZ,
	             corePsi, last.psiAxis );
	std::fflush( stdout );

	BOOST_TEST( std::abs( corePsi - last.psiAxis ) < 5.0e-2*span,
		"psi_h at the reference equilibrium's magnetic axis is " << corePsi
		<< " where this solve's psi_ax is " << last.psiAxis << ", "
		<< std::abs( corePsi - last.psiAxis )/span << " of the span apart. The "
		"core of this solve is not where the reference's is, so the solve has "
		"converged to a different equilibrium and every border unknown above "
		"describes it rather than this machine." );

	// FIVE: AND THE EQUILIBRIUM IS THE ONE THE SHIPPED MACHINE CASE CALLS
	// HEALTHY, by the same three checks.
	meq::GradShafranovSolver::AxisSourceCheck const axisSource =
		m.solver->checkAxisSource();
	meq::CriticalPointFinder finder( *m.solver );
	meq::AxisAgreement const axis =
		finder.checkAxis( m.solver->psiAxis(), m.solver->psiBoundary() );

	std::printf( "    | F | on R = 0 is %.3e, psi_ax attained at ( %.3f, %.3f ),"
	             " Psi at the O-point %.4f over %d extrema and %d saddles, "
	             "mu0 I_p delivered %.6e\n", axisSource.worstOnAxis, axis.nodeR,
	             axis.nodeZ, axis.normalisedFlux, axis.extrema, axis.saddles,
	             m.solver->plasmaCurrent() );
	std::fflush( stdout );

	BOOST_TEST( axisSource.bounded,
		"| F | on the symmetry axis is " << axisSource.worstOnAxis
		<< ", so F/R = mu_0 j_phi is an unbounded toroidal current density on "
		"R = 0." );
	BOOST_TEST( axis.nodeR > 0.30,
		"psi_ax is attained at R = " << axis.nodeR << ", which is on or beside "
		"the symmetry axis -- section 11.3's axis layer." );
	BOOST_TEST( axis.agrees,
		"the O-point of q_h carries Psi = " << axis.normalisedFlux
		<< " against the 1 it must carry by definition." );

	double const wanted = m.mu0*m.config->getSource().getMHD().plasmaCurrent;
	BOOST_TEST( std::abs( m.solver->plasmaCurrent() - wanted )
	            < 1.0e-5*std::abs( wanted ),
		"the delivered current is " << m.solver->plasmaCurrent()
		<< " against the " << wanted << " prescribed" );
}
