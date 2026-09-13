/*
 * XP-3: THE X-POINT IS TWO MORE UNKNOWNS OF THE SAME NEWTON.
 * FREE-BOUNDARY-PLAN.md sections 10.4 and 10.6.
 *
 * XP-2 pins psi_bnd at the current estimate of the X-point, solves, locates the
 * saddle of the solved q_h, re-pins and re-solves -- an OUTER fixed point over
 * ( r_X, z_X ). This is the same statement made INSIDE the Newton:
 *
 *     q_r( r_X, z_X )             = 0
 *     q_z( r_X, z_X )             = 0
 *     psi_bnd - psi_h( r_X, z_X ) = 0
 *
 * three rows and three unknowns closing with everything else on one
 * factorisation per step. meq::GradShafranovSolver::setXPointBoundary is the
 * whole of the interface and its documentation carries the derivation.
 *
 *
 * WHAT THIS CASE ASSERTS, AND WHY IT IS NOT AN ERROR NORM.
 *
 * The plan's own acceptance is "agreement with XP-2 at round-off, and the
 * observed Newton order, which is the only thing that can see the inexact
 * grad q corner block". Both are here, and so is the reason the second matters:
 * CLAUDE.md's *A wrong Jacobian is invisible to a convergence table* -- Newton
 * converges to the same discrete solution whatever Jacobian carried it there,
 * so a corner block that is wrong, or absent, moves the ITERATION COUNT and the
 * observed order and moves no number in the answer. An equilibrium comparison
 * alone would pass on a border with the whole corner block zeroed.
 *
 *
 * AND THE PLAN'S PREDICTION ABOUT THAT CORNER BLOCK IS FALSIFIED HERE, WHICH IS
 * THE ONE FINDING OF THE STAGE.
 *
 * Section 10.4 expects to lose the quadratic rate: "the corner block
 * d( q_r, q_z )/d( r_X, z_X ) is grad q -- the Hessian of the potential -- and
 * there is no solved variable for it: differentiating an L2 field of degree k
 * leaves k - 1. This is the same wall recorded for the band continuation of B."
 * It then plans a fallback of differencing the two rows.
 *
 * The band analogy does not carry, and the difference is between approximating
 * a CONTINUOUS object and differentiating a DISCRETE one. The band continuation
 * needs grad q as an approximation of the true Hessian, and there it really is
 * an order down. Newton needs the derivative of the discrete residual with
 * respect to the discrete unknowns -- and q_h is a polynomial on its element,
 * so dq_h/dx there is exact arithmetic, not an approximation of anything. The
 * corner block is exact, nothing on this path is differenced, and the measured
 * order says so.
 *
 * WHAT IS GENUINELY NOT SMOOTH IS THE ELEMENT CHANGING. q_h is discontinuous
 * across a face, so carrying the point over a mesh line changes the polynomial
 * the rows evaluate, by O( h^( k+1 ) ). That is a floor on the order rather
 * than a cap on it, it only bites while the point is still travelling, and it
 * is why elementDepth() is a precondition on this fixture exactly as it is on
 * XP-2's.
 *
 *
 * THE SUPPORT IS STILL AN OUTER LOOP AND THAT IS THE POINT OF THE STAGE.
 *
 * Section 10.5 names two discrete states on XP-2's fixed point: which critical
 * point bounds the plasma, and which elements carry current. M-82 measured both
 * as necessary -- the outer loop fails at the 200-iteration cap with either one
 * moving inside the Newton. XP-3 removes the FIRST from the outer loop by
 * making the X-point's POSITION differentiable, which it is; it does not touch
 * the second, which is not. So this case still freezes the support and still
 * re-decides it between solves, and what it demonstrates is one outer state
 * where XP-2 has two.
 */
#define BOOST_TEST_MODULE XPointBorder
#include <boost/test/unit_test.hpp>

#include <cstdio>

// The machine, the reference values, the saddle filter and the support freeze,
// shared with XPointOuter.cpp. The SAME fixture and not a copy of it: this
// case's headline assertion is agreement with that one, and a second copy would
// make a disagreement ambiguous between the two borders and the two fixtures.
#include "convergence/DivertedMachine.hpp"
#include "convergence/ConvergenceHarness.hpp"

using namespace meqtest;

namespace
{
	/// THE BOOTSTRAP PIN, AND IT IS XP-2's OWN so that the two cases start from
	/// the same distance away. 0.05 m in each coordinate is what a machine's
	/// drawings would get you to: about two thirds of an element here, so the
	/// border has somewhere to travel and travels across at least one face
	/// getting there.
	double const bootstrapR = referenceXPointR + 0.05;
	double const bootstrapZ = referenceXPointZ + 0.05;

	/// One solve of the machine with the X-point in the border, with the plasma
	/// support frozen at @a state and re-decided here rather than inside Newton.
	///
	/// setXPointBoundary() is called ONCE, on the first sweep: the solver keeps
	/// the X-point where the last solve left it, so a re-solve continues from
	/// it. Calling it again is refused, which is the interface saying that the
	/// point is an unknown rather than a setting.
	bool solveWithBorder( Machine &m, bool first, mfem::GridFunction const &state,
	                      double axis, double boundary )
	{
		freezeSupportAt( m, state, axis, boundary );
		if ( first )
			m.solver->setXPointBoundary( bootstrapR, bootstrapZ );
		else
			m.solver->setInitialGuess( state );

		try
		{
			m.solver->solve();
			m.converged = !m.solver->newtonResiduals().empty()
			              && m.solver->newtonResiduals().back() < 1.0e-8;
		}
		catch ( std::exception const &error )
		{
			std::printf( "      solve threw: %s\n", error.what() );
			m.converged = false;
		}
		return m.converged;
	}

	void printHistory( std::vector<double> const &history )
	{
		for ( std::size_t i = 0; i < history.size(); ++i )
		{
			std::printf( "        %2zu  %14.6e", i, history[ i ] );
			if ( i >= 2 && history[ i - 1 ] > 0.0 && history[ i - 2 ] > 0.0 )
				std::printf( "   order %6.3f",
				             meq::tests::newtonOrder( history[ i - 2 ],
				                                      history[ i - 1 ],
				                                      history[ i ] ) );
			std::printf( "\n" );
		}
		std::fflush( stdout );
	}

	/// The last triple whose middle and last terms are above the round-off
	/// floor. The TAIL rather than the best triple anywhere, for LimiterCurve's
	/// reason: a best-of over a wandering history reports the one place it
	/// happened to look quadratic.
	double tailOrder( std::vector<double> const &history )
	{
		if ( history.size() < 3 )
			return 0.0;
		double const floorLevel = 1.0e-13*history.front();
		double tail = 0.0;
		for ( std::size_t i = 2; i < history.size(); ++i )
			if ( history[ i ] > floorLevel && history[ i - 1 ] > floorLevel )
				tail = meq::tests::newtonOrder( history[ i - 2 ], history[ i - 1 ],
				                                history[ i ] );
		return tail;
	}

	/// Everything one run of the machine leaves behind that the assertions read.
	struct Answer
	{
		double xR = 0.0;
		double xZ = 0.0;
		double psiAxis = 0.0;
		double psiBoundary = 0.0;
		/// The saddle an INDEPENDENT root find reports on the solved field, and
		/// how far into its element it sits.
		double saddleR = 0.0;
		double saddleZ = 0.0;
		double saddlePsi = 0.0;
		double depth = 0.0;
		bool saddleFound = false;
		int sweeps = 0;
		int supportMoved = 0;
		std::size_t iterations = 0;
		double order = 0.0;
		/// The order off the FIRST sweep, which is the only history long enough
		/// to carry one -- see theXPointBorderClosesOnTheSaddle's fifth
		/// assertion. A warm sweep reaches the floor in two steps.
		double bootstrapOrder = 0.0;
		std::size_t bootstrapIterations = 0;
	};

	/// XP-3, end to end: bootstrap, then the SUPPORT alone as an outer state.
	Answer runBorder( Machine &m )
	{
		Answer answer;

		double bootstrapEdge = 0.0;
		BOOST_TEST_REQUIRE( valueAt( *m.carried, bootstrapR, bootstrapZ,
		                             bootstrapEdge ),
		                    "the transferred guess cannot be evaluated at the "
		                    "bootstrap pin, so there is no edge to freeze at" );

		std::printf( "\n  XP-3: the X-POINT IS IN THE BORDER, bootstrapped at "
		             "( %.3f, %.3f )\n", bootstrapR, bootstrapZ );
		std::printf( "    %-6s %5s %19s %11s %13s %13s %7s %8s\n", "sweep",
		             "its", "X-point", "moved", "psi_ax", "psi_bnd", "order",
		             "support" );
		std::fflush( stdout );

		mfem::GridFunction previous( *m.carried );
		double axis = m.config->getSource().psiAxisGuess();
		double boundary = bootstrapEdge;

		double lastR = bootstrapR;
		double lastZ = bootstrapZ;

		for ( int sweep = 0; sweep < 6; ++sweep )
		{
			int const supportBefore = m.solver->plasmaComponentElements();
			if ( !solveWithBorder( m, sweep == 0, previous, axis, boundary ) )
			{
				std::printf( "    %-6d %5s %19s\n", sweep, "-", "NO SOLVE" );
				return answer;
			}
			answer.supportMoved +=
				( m.solver->plasmaComponentElements() != supportBefore );

			answer.xR = m.solver->xPointR();
			answer.xZ = m.solver->xPointZ();
			answer.psiAxis = m.solver->psiAxis();
			answer.psiBoundary = m.solver->psiBoundary();
			answer.iterations = m.solver->newtonResiduals().size() - 1;
			answer.order = tailOrder( m.solver->newtonResiduals() );
			if ( sweep == 0 )
			{
				answer.bootstrapOrder = answer.order;
				answer.bootstrapIterations = answer.iterations;
			}
			double const moved = std::hypot( answer.xR - lastR,
			                                 answer.xZ - lastZ );

			std::printf( "    %-6d %5zu   (%6.3f,%7.3f) %11.3e %13.6e %13.6e "
			             "%7.3f %4d/%d\n", sweep, answer.iterations, answer.xR,
			             answer.xZ, moved, answer.psiAxis, answer.psiBoundary,
			             answer.order, m.solver->plasmaComponentElements(),
			             m.solver->plasmaCandidateElements() );
			std::printf( "      Newton history of sweep %d\n", sweep );
			printHistory( m.solver->newtonResiduals() );

			previous = m.solver->potential();
			axis = answer.psiAxis;
			boundary = answer.psiBoundary;
			lastR = answer.xR;
			lastZ = answer.xZ;
			++answer.sweeps;

			// The support is the only outer state left, so the loop stops when
			// it stops moving -- and one more sweep after that, whose Newton
			// history is the one the order is read off, since by then nothing
			// outside the solve has changed under it.
			if ( sweep > 0 && moved < 1.0e-10
			     && m.solver->plasmaComponentElements() == supportBefore )
				break;
		}

		/*
		 * THE FIXED-POINT PROPERTY, AND IT IS WHAT XP-2 SPENDS FOUR SWEEPS
		 * REACHING. The border makes psi_bnd equal psi_h at the point it carries
		 * whatever that point is; what is NOT automatic is that the point is a
		 * saddle of the solved q_h. An independent root find is what says so,
		 * and it is the same one XP-2's outer loop is built out of.
		 */
		meq::CriticalPointFinder finder( *m.solver );
		meq::CriticalPoint saddle;
		answer.saddleFound = finder.tryFindCriticalPointFrom(
			answer.xR, answer.xZ, meq::AxisSense::Saddle, saddle );
		if ( !answer.saddleFound )
			answer.saddleFound = findXPoint( finder, saddle );
		if ( answer.saddleFound )
		{
			answer.saddleR = saddle.r;
			answer.saddleZ = saddle.z;
			answer.saddlePsi = saddle.psi;
			answer.depth = elementDepth( saddle );
		}

		return answer;
	}

	/// XP-2, end to end, for the comparison. It is XPointOuter's loop written
	/// again in twenty lines rather than called, because that case asserts as it
	/// goes and what is wanted here is only its answer.
	Answer runOuterLoop( Machine &m )
	{
		Answer answer;

		double bootstrapEdge = 0.0;
		BOOST_TEST_REQUIRE( valueAt( *m.carried, bootstrapR, bootstrapZ,
		                             bootstrapEdge ),
		                    "the transferred guess cannot be evaluated at the "
		                    "bootstrap pin" );
		freezeSupportAt( m, *m.carried, m.config->getSource().psiAxisGuess(),
		                 bootstrapEdge );

		std::printf( "\n  XP-2 ON THE SAME FIXTURE, FOR THE COMPARISON\n" );
		std::fflush( stdout );

		if ( !solveAt( m, bootstrapR, bootstrapZ ) )
			return answer;

		// THE COMPARABLE HISTORY IS THE BOOTSTRAP'S, not the loop's first sweep.
		// It starts from the transferred guess exactly as XP-3's first solve
		// does; every sweep after it is warm-started from a converged answer and
		// reaches the floor in two or three steps, where no order exists to be
		// read. So this is the control for XP-3's own bootstrap order.
		answer.bootstrapOrder = tailOrder( m.solver->newtonResiduals() );
		answer.bootstrapIterations = m.solver->newtonResiduals().size() - 1;
		std::printf( "    bootstrap: %zu iterations, tail order %.3f\n",
		             answer.bootstrapIterations, answer.bootstrapOrder );
		printHistory( m.solver->newtonResiduals() );

		meq::CriticalPointFinder first( *m.solver );
		meq::CriticalPoint x;
		if ( !findXPoint( first, x ) )
			return answer;

		mfem::GridFunction previous( m.solver->potential() );
		for ( int sweep = 0; sweep < 8; ++sweep )
		{
			double const pinnedR = x.r;
			double const pinnedZ = x.z;
			freezeSupportAt( m, previous, m.solver->psiAxis(),
			                 m.solver->psiBoundary() );
			if ( !solveAt( m, pinnedR, pinnedZ, &previous ) )
				return answer;

			meq::CriticalPointFinder finder( *m.solver );
			meq::CriticalPoint next;
			if ( !finder.tryFindCriticalPointFrom( pinnedR, pinnedZ,
			                                       meq::AxisSense::Saddle, next )
			     && !findXPoint( finder, next ) )
				return answer;

			double const step = std::hypot( next.r - pinnedR, next.z - pinnedZ );
			x = next;
			previous = m.solver->potential();
			++answer.sweeps;

			answer.xR = x.r;
			answer.xZ = x.z;
			answer.psiAxis = m.solver->psiAxis();
			answer.psiBoundary = m.solver->psiBoundary();
			// NOT bootstrapOrder, which the bootstrap solve above already set:
			// sweep 0 of THIS loop is already warm-started from it.
			answer.order = tailOrder( m.solver->newtonResiduals() );
			answer.iterations = m.solver->newtonResiduals().size() - 1;
			std::printf( "    sweep %d: %zu iterations, ( %.6f, %.6f ), step "
			             "%.3e, psi_ax %.6e, psi_bnd %.6e\n", sweep + 1,
			             answer.iterations, x.r, x.z, step, answer.psiAxis,
			             answer.psiBoundary );
			// THE CONTROL ON THE ORDER, and the reason it is printed rather than
			// only summarised: XP-2's inner solve carries the SAME border as
			// XP-3's less the two flux rows, so whatever its history does to the
			// order is the machine and the line search rather than anything this
			// stage added.
			printHistory( m.solver->newtonResiduals() );
			std::fflush( stdout );

			if ( step < 1.0e-10 )
				break;
		}

		answer.saddleFound = answer.sweeps > 0;
		answer.saddleR = answer.xR;
		answer.saddleZ = answer.xZ;
		return answer;
	}
}

/*
 * THE STAGE ITSELF: the border closes, on the machine, against freegs4e -- and
 * the two things only an inner border can say.
 */
BOOST_AUTO_TEST_CASE( theXPointBorderClosesOnTheSaddle )
{
	Machine m = buildMachine();
	completeMachine( m );

	Answer const border = runBorder( m );

	BOOST_TEST_REQUIRE( border.sweeps > 0,
		"the diverted machine did not converge with the X-point in the border. "
		"M-82 measured the two things this needs -- XP-1's fill running and the "
		"plasma edge frozen within each solve -- and the fixture supplies both; "
		"a failure here with those in place is the border and not the "
		"configuration." );
	BOOST_TEST_REQUIRE( border.saddleFound,
		"no off-axis saddle in the solved field, so this configuration is not "
		"diverted and the border closed on something else" );

	std::printf( "\n  WHERE IT CLOSED\n" );
	std::printf( "    border    ( %.6f, %.6f )\n", border.xR, border.xZ );
	std::printf( "    root find ( %.6f, %.6f ), %.3e apart, depth %.3f\n",
	             border.saddleR, border.saddleZ,
	             std::hypot( border.saddleR - border.xR,
	                         border.saddleZ - border.xZ ), border.depth );
	std::printf( "    freegs4e  ( %.6f, %.6f ), %.3e apart\n",
	             referenceXPointR, referenceXPointZ,
	             std::hypot( referenceXPointR - border.xR,
	                         referenceXPointZ - border.xZ ) );
	std::printf( "    psi_bnd %.6e against the saddle's own psi_h %.6e, "
	             "%.2e apart\n", border.psiBoundary, border.saddlePsi,
	             std::abs( border.saddlePsi - border.psiBoundary ) );
	std::printf( "    psi_ax  %.6e against freegs4e's %.6e\n", border.psiAxis,
	             referencePsiAxis );
	std::fflush( stdout );

	/*
	 * ONE: THE POINT THE BORDER CARRIES IS A SADDLE OF THE FIELD IT SOLVED.
	 *
	 * This is the whole content of the stage. XP-2 reaches it by iterating a
	 * root find against a solve; here the root find IS two rows of the solve, so
	 * the two coincide to whatever the root finder's own tolerance is rather
	 * than to whatever the outer loop was stopped at.
	 *
	 * The bound is a fraction of an element rather than round-off, and
	 * deliberately: CriticalPointFinder runs its own Newton to its own tolerance
	 * and, near a face, can report the neighbouring element's root -- which is
	 * a different polynomial's zero and legitimately O( h^( k+1 ) ) away. See
	 * elementDepth().
	 */
	double const apart = std::hypot( border.saddleR - border.xR,
	                                 border.saddleZ - border.xZ );
	BOOST_TEST( apart < 1.0e-6,
		"the border left the X-point at ( " << border.xR << ", " << border.xZ
		<< " ) and an independent root find on the SAME solved field puts the "
		"saddle at ( " << border.saddleR << ", " << border.saddleZ << " ), "
		<< apart << " m away. The two rows q_r = q_z = 0 are exactly the "
		"equations that root find solves, so a disagreement is the rows being "
		"built wrong -- suspect the sign that undoes DarcyForm's -q, which is "
		"right in magnitude and wrong in direction if it is missed." );

	/*
	 * TWO: psi_bnd IS THE FLUX THERE, which is the third row and is the reason
	 * the other two are worth having. A border that found the saddle and left
	 * psi_bnd somewhere else would report an equilibrium bounded by nothing.
	 */
	double const span = border.psiAxis - border.psiBoundary;
	BOOST_TEST( std::abs( border.saddlePsi - border.psiBoundary )
	            < 1.0e-6*std::abs( span ),
		"psi_bnd is " << border.psiBoundary << " and psi_h at the located "
		"saddle is " << border.saddlePsi << ", which is " << ( ( border.saddlePsi
		- border.psiBoundary )/span ) << " of the span apart. The third border "
		"row is psi_bnd - psi_h( x ) = 0 and it is solved to the same tolerance "
		"as the rest, so this is the row not reaching the same point the other "
		"two do." );

	/*
	 * THREE: THE ELEMENT IS NOT A MESH LINE. XPointOuter's own precondition, and
	 * it matters more here than there: the border rows are built from ONE
	 * element's polynomial, so a converged point sitting on a face has two of
	 * them and the Jacobian would be the neighbour's on alternate steps.
	 */
	BOOST_TEST( border.depth > 0.05,
		"the converged X-point sits " << border.depth << " into its element in "
		"reference coordinates, so it is on a mesh line and two elements hold "
		"it. That is a statement about where this machine's null falls on this "
		"mesh rather than about the border: move the designed target." );

	/*
	 * FOUR: AGAINST AN INDEPENDENT CODE. The same 4.4e-04 m XP-2 reaches, which
	 * is what says the inner border found the same physical null and not merely
	 * a self-consistent one. The bound is loose because it is a statement about
	 * two discretisations of one machine -- freegs4e's 129^2 finite differences
	 * against k = 2 HDG on 1333 elements -- and not about either code's own
	 * convergence.
	 */
	double const fromReference = std::hypot( referenceXPointR - border.xR,
	                                         referenceXPointZ - border.xZ );
	BOOST_TEST( fromReference < 5.0e-03,
		"the border's X-point is " << fromReference << " m from freegs4e's. "
		"XP-2 reaches 4.378e-04 on this fixture, so a figure much larger than "
		"that is a different null or a different equilibrium rather than a "
		"different discretisation of this one." );
	BOOST_TEST( std::abs( border.psiBoundary/referencePsiBoundary - 1.0 ) < 0.01,
		"psi_bnd is " << border.psiBoundary << " against freegs4e's "
		<< referencePsiBoundary );
	BOOST_TEST( std::abs( border.psiAxis/referencePsiAxis - 1.0 ) < 0.01,
		"psi_ax is " << border.psiAxis << " against freegs4e's "
		<< referencePsiAxis );

	/*
	 * FIVE: THE OBSERVED ORDER, WHICH IS THE ONLY ASSERTION HERE THAT CAN SEE
	 * THE CORNER BLOCK AT ALL -- AND IT IS READ OFF THE BOOTSTRAP.
	 *
	 * Everything above is satisfied by a border whose corner block is merely
	 * WRONG: Newton reaches the same discrete solution whatever Jacobian
	 * carried it there, which is CLAUDE.md's *A wrong Jacobian is invisible to
	 * a convergence table*. So this is the assertion the stage stands on.
	 *
	 * **IT IS THE FIRST SWEEP AND NOT THE LAST, AND THE REASON IS A PROPERTY OF
	 * THE FIXTURE RATHER THAN A CHOICE.** Every sweep after the first is
	 * warm-started from a converged answer and reaches the solve's own floor in
	 * TWO steps -- 8.26e-03, 1.96e-08, 7.15e-12 -- so its "tail triple" is the
	 * whole history and its last term is the floor. An order read off that
	 * measures where the floor is, and reads 0.61 for XP-2's identical
	 * two-step sweeps as readily as for XP-3's. The bootstrap is the only
	 * history on this fixture with an asymptotic regime in it.
	 *
	 * **AND THE BOUND IS 1.2 RATHER THAN 1.5, WHICH IS THIS MACHINE'S CEILING
	 * AND NOT THIS BORDER'S.** Measured, XP-2's own bootstrap -- the same
	 * bordered system less exactly these two rows -- runs in the same band, and
	 * theBorderAndTheOuterLoopReachTheSameEquilibrium asserts the comparison
	 * directly rather than leaving it to two absolute bounds. What caps it is
	 * not known and is not XP-3's: the candidates are the axis argmax hopping
	 * elements, an O( h^( k+1 ) ) kink CLAUDE.md already records, and the
	 * element-granular plasma support. A clean quadratic run on a bordered
	 * free-boundary solve is LimiterCurve's, on a fixture with none of that.
	 */
	std::printf( "    tail order: bootstrap %.3f over %zu iterations, last "
	             "sweep %.3f over %zu\n", border.bootstrapOrder,
	             border.bootstrapIterations, border.order, border.iterations );
	BOOST_TEST( border.bootstrapOrder > 1.2,
		"the bordered Newton's tail order on the bootstrap is "
		<< border.bootstrapOrder << ". A rate near 1 is a Jacobian statement, "
		"and the entries this stage owns are the two flux rows and the corner "
		"block grad q -- so the first thing to suspect is the corner block, "
		"which FREE-BOUNDARY-PLAN.md section 10.4 expects to be inexact and "
		"which this file argues is exact. If it really cannot be made exact, "
		"the plan's own fallback is to difference the two rows in "
		"( r_X, z_X ), which is two residual evaluations in a 2-vector. Read "
		"the printed history before believing this number: a run that damped "
		"its way through the first half has no asymptotic regime." );
	BOOST_TEST( border.bootstrapOrder < 3.5,
		"the tail order reads " << border.bootstrapOrder << ", which is not a "
		"rate -- it is a short or non-monotone history being read as one. The "
		"printed history above is the thing to look at." );

	/*
	 * SIX: ONE OUTER STATE WHERE XP-2 HAS TWO. The support is still re-decided
	 * between solves -- section 10.5's combinatorial state, which has no
	 * globalisation -- and it has to settle for psi_bnd to mean anything. The
	 * X-point is no longer on that list, which is the stage.
	 */
	std::printf( "    the frozen support changed size on %d of %d sweeps\n",
	             border.supportMoved, border.sweeps );
	BOOST_TEST( border.supportMoved < border.sweeps,
		"the frozen support changed size on every one of " << border.sweeps
		<< " sweeps, so the one outer state this stage leaves has not settled "
		"on a topology -- section 10.5's alternating discrete state. The "
		"X-point's own convergence says nothing about this: they are separate "
		"states and either can chatter while the other closes." );
}

/*
 * AGREEMENT WITH XP-2, WHICH IS THE PLAN'S OWN ACCEPTANCE.
 *
 * The two solve the same problem by different means: the outer loop alternates
 * a root find with a solve until they stop moving each other, and the border
 * makes the root find part of the solve. They must therefore reach the same
 * equilibrium, and the only thing that can separate them is the tolerance each
 * stops at -- which is why the bound here is far tighter than the one against
 * freegs4e above and far looser than round-off.
 *
 * IT RUNS BOTH IN ONE BINARY RATHER THAN QUOTING XPointOuter's PRINTED ANSWER.
 * A hardcoded number from another case is a number that goes stale silently:
 * the fixture is shared, so anything that moves the machine moves both, and
 * only a comparison taken in one process can say so.
 */
BOOST_AUTO_TEST_CASE( theBorderAndTheOuterLoopReachTheSameEquilibrium )
{
	Machine inner = buildMachine();
	completeMachine( inner );
	Answer const border = runBorder( inner );
	BOOST_TEST_REQUIRE( border.sweeps > 0, "the border solve did not converge" );

	Machine outer = buildMachine();
	completeMachine( outer );
	Answer const loop = runOuterLoop( outer );
	BOOST_TEST_REQUIRE( loop.sweeps > 0,
		"XP-2's outer loop did not converge on this fixture, so there is "
		"nothing to compare against. XPointOuter is the case that owns that "
		"failure." );

	double const apart = std::hypot( border.xR - loop.xR, border.xZ - loop.xZ );
	double const span = loop.psiAxis - loop.psiBoundary;

	std::printf( "\n  XP-3 AGAINST XP-2, SAME FIXTURE, SAME PROCESS\n" );
	std::printf( "    %-8s %19s %13s %13s %6s\n", "route", "X-point", "psi_ax",
	             "psi_bnd", "its" );
	std::printf( "    %-8s ( %.6f,%10.6f ) %13.6e %13.6e %6zu\n", "border",
	             border.xR, border.xZ, border.psiAxis, border.psiBoundary,
	             border.iterations );
	std::printf( "    %-8s ( %.6f,%10.6f ) %13.6e %13.6e %6zu\n", "outer",
	             loop.xR, loop.xZ, loop.psiAxis, loop.psiBoundary,
	             loop.iterations );
	std::printf( "    apart: %.3e m, psi_ax %.3e, psi_bnd %.3e of the span\n",
	             apart, std::abs( border.psiAxis - loop.psiAxis )/span,
	             std::abs( border.psiBoundary - loop.psiBoundary )/span );
	std::printf( "    sweeps: border %d, outer %d\n", border.sweeps,
	             loop.sweeps );
	std::printf( "    bootstrap: border %zu iterations at order %.3f, outer "
	             "%zu at %.3f\n", border.bootstrapIterations,
	             border.bootstrapOrder, loop.bootstrapIterations,
	             loop.bootstrapOrder );
	std::fflush( stdout );

	BOOST_TEST( apart < 1.0e-5,
		"the border puts the X-point at ( " << border.xR << ", " << border.xZ
		<< " ) and XP-2's outer loop at ( " << loop.xR << ", " << loop.xZ
		<< " ), " << apart << " m apart. They solve the same equations, so a "
		"disagreement is one of them not having reached its own fixed point -- "
		"the printed histories say which." );
	BOOST_TEST( std::abs( border.psiAxis - loop.psiAxis ) < 1.0e-6*std::abs( span ),
		"psi_ax disagrees between the two routes by "
		<< std::abs( border.psiAxis - loop.psiAxis )/span << " of the span" );
	BOOST_TEST( std::abs( border.psiBoundary - loop.psiBoundary )
	            < 1.0e-6*std::abs( span ),
		"psi_bnd disagrees between the two routes by "
		<< std::abs( border.psiBoundary - loop.psiBoundary )/span
		<< " of the span" );

	/*
	 * AND THE ORDER, AS A COMPARISON RATHER THAN AS A BOUND, WHICH IS THE ONLY
	 * HONEST WAY TO READ IT ON THIS FIXTURE.
	 *
	 * The two bootstraps solve the same problem from the same guess with the
	 * same border, differing in exactly the two flux rows and the corner block
	 * this stage adds. So whatever caps the observed order here is common to
	 * both, and the difference between them is XP-3's alone. A corner block
	 * that was wrong -- not missing, which makes the dense system singular,
	 * but wrong -- would show up here and nowhere else in this file.
	 */
	BOOST_TEST( border.bootstrapOrder > loop.bootstrapOrder - 0.2,
		"the border's bootstrap runs at observed order "
		<< border.bootstrapOrder << " against XP-2's " << loop.bootstrapOrder
		<< " on the same fixture from the same guess. The two differ by the "
		"two flux rows and the corner block grad q, so a LOWER order is that "
		"corner block being wrong -- the one thing an equilibrium comparison "
		"cannot see." );
}

/*
 * WHAT THE INTERFACE REFUSES, AND THE CONTROL THAT KEEPS THE REFUSALS FROM
 * BEING VACUOUS.
 *
 * All three of setBoundaryFluxPoint(), setLimiterSurface() and
 * setXPointBoundary() pin psi_bnd, and a diverted plasma's limiter is OUTSIDE
 * its separatrix -- so naming two of them is a choice between two equilibria
 * made by call order. XPointOuter's own header records what the wrong one costs:
 * every hand-designed single null in this tree pinned psi_bnd at the limiter and
 * every one of them was asked the wrong question rather than answering it badly.
 */
BOOST_AUTO_TEST_CASE( theXPointBorderRefusesASecondConstraintOnTheSameUnknown )
{
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 4, 4, mfem::Element::TRIANGLE,
	                                               false, 1.0, 1.0 );
	for ( int v = 0; v < mesh.GetNV(); ++v )
		mesh.GetVertex( v )[ 0 ] += 1.0;

	{
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setXPointBoundary( 1.5, 0.5 );
		BOOST_CHECK_THROW( solver.setBoundaryFluxPoint( 1.5, 0.5 ),
		                   std::logic_error );
		BOOST_CHECK_THROW( solver.setLimiterSurface( 1 ), std::logic_error );
		BOOST_CHECK_THROW( solver.setXPointBoundary( 1.5, 0.5 ),
		                   std::logic_error );
	}
	{
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setBoundaryFluxPoint( 1.5, 0.5 );
		BOOST_CHECK_THROW( solver.setXPointBoundary( 1.5, 0.5 ),
		                   std::logic_error );
	}

	// THE CONTROL. Without it every refusal above passes on a setter that
	// throws unconditionally, which is the failure mode a refusal test has.
	{
		meq::GradShafranovSolver solver( mesh, 1 );
		BOOST_CHECK_NO_THROW( solver.setXPointBoundary( 1.5, 0.5 ) );
		BOOST_TEST( solver.xPointIsAnUnknown() );
		BOOST_TEST( solver.xPointR() == 1.5 );
		BOOST_TEST( solver.xPointZ() == 0.5 );
		BOOST_TEST( !solver.xPointWasLocated(),
			"xPointWasLocated() is true before any solve, so it reports the "
			"SETTING rather than what a solve did with it" );
	}

	// AND THE AXIS IS NOT AN X-POINT: psi vanishes identically on r = 0, so q_h
	// is small along the whole of it and a sweep reports a ladder of
	// near-saddles that no divertor put there.
	{
		meq::GradShafranovSolver solver( mesh, 1 );
		BOOST_CHECK_THROW( solver.setXPointBoundary( 0.0, 0.5 ),
		                   std::invalid_argument );
	}
}
