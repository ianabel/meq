#define BOOST_TEST_MODULE WarmStartConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Sampler.hpp"
#include "meq/Source.hpp"
#include "meq/WarmStart.hpp"

#include "analytic/ManufacturedNonlinear.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * STAGE 7d: RESTARTING A SOLVE FROM ANOTHER ONE, AT THE ORDER IT WAS COMPUTED
 * WITH.
 *
 * DRIVER-PLAN.md section 4 argues that the full-order route is worth a
 * dependency and the structured-grid route is not a substitute for it. This file
 * is that argument as measurements, because the difference is exactly the sort
 * that a working restart hides: BOTH routes converge, and the cheap one simply
 * does less than it looks. A guess is only a guess, so nothing fails -- the
 * iteration count is a little higher and nobody investigates.
 *
 * The problem throughout is HDG-GradShafranov Example 5, whose closed form makes
 * "how good is the guess" answerable rather than relative.
 */

namespace
{
	using meq::tests::EquilibriumSource;

	double const rMin = 0.6;
	double const rMax = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ rMin, rMax, zMin, zMax };
	}

	mfem::Mesh makeMesh( int n )
	{
		return meq::tests::makeMesh( box(), n );
	}

	meq::analytic::ManufacturedNonlinear equilibrium()
	{
		return meq::analytic::ManufacturedNonlinear::example5();
	}

	/*
	 * THE INTERCHANGE ROUTE'S ARITHMETIC: psi on a structured ( R, Z ) grid, read
	 * back by bilinear interpolation. This is what a NetCDF file from another
	 * code offers and all it can offer -- no basis, no mesh, no elements.
	 *
	 * It is written here rather than in src/ because meq does not yet READ its
	 * own NetCDF files; what is under test is the arithmetic that route would
	 * use, which is what decides whether it is a substitute for the full-order
	 * one. Second order in the GRID spacing, whatever the solve's degree was.
	 */
	class BilinearGridCoefficient : public mfem::Coefficient
	{
		public:
			BilinearGridCoefficient( meq::GridSampler const &sampler,
			                         std::vector<double> values,
			                         double outsideIn )
				: nR( sampler.nodesR() ), nZ( sampler.nodesZ() ),
				  data( std::move( values ) ), outside( outsideIn )
			{
				dR = ( rMax - rMin )/( nR - 1 );
				dZ = ( zMax - zMin )/( nZ - 1 );
			}

			double Eval( mfem::ElementTransformation &tr, // NOLINT(readability-identifier-naming)
			             mfem::IntegrationPoint const &ip ) override
			{
				mfem::Vector x( 2 );
				tr.Transform( ip, x );

				double const fi = ( x( 0 ) - rMin )/dR;
				double const fj = ( x( 1 ) - zMin )/dZ;
				int i = static_cast<int>( std::floor( fi ) );
				int j = static_cast<int>( std::floor( fj ) );
				i = std::min( std::max( i, 0 ), nR - 2 );
				j = std::min( std::max( j, 0 ), nZ - 2 );

				double const s = std::min( std::max( fi - i, 0.0 ), 1.0 );
				double const t = std::min( std::max( fj - j, 0.0 ), 1.0 );

				double const v00 = at( i, j ), v10 = at( i + 1, j );
				double const v01 = at( i, j + 1 ), v11 = at( i + 1, j + 1 );
				if ( !std::isfinite( v00 ) || !std::isfinite( v10 )
				     || !std::isfinite( v01 ) || !std::isfinite( v11 ) )
					return outside;

				return ( 1.0 - s )*( 1.0 - t )*v00 + s*( 1.0 - t )*v10
				       + ( 1.0 - s )*t*v01 + s*t*v11;
			}

		private:
			double at( int i, int j ) const { return data[ j*nR + i ]; }

			int nR, nZ;
			std::vector<double> data;
			double outside;
			double dR, dZ;
	};

	/// A converged Example 5 solve, kept with everything a restart needs.
	struct Solved
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::GradShafranovSolver> solver;
	};

	/// A once-refined mesh, built fresh each time, so that no two live solvers
	/// share one. Nothing here has been shown to require that; it is the cheap
	/// side of a question not worth opening in a test about something else.
	mfem::Mesh refinedMesh( int n )
	{
		mfem::Mesh mesh = makeMesh( n );
		mesh.UniformRefinement();
		return mesh;
	}

	Solved solve( int order, int n, meq::Source const &source,
	              mfem::Coefficient &datum )
	{
		Solved out;
		out.mesh = std::make_unique<mfem::Mesh>( makeMesh( n ) );
		out.solver = std::make_unique<meq::GradShafranovSolver>( *out.mesh, order );
		out.solver->setSource( source );
		out.solver->setBoundaryData( datum );
		out.solver->solve();
		return out;
	}
}

/*
 * THE TRANSFER IS EXACT ONTO A REFINED MESH, AND THAT IS THE SHARPEST THING
 * THAT CAN BE ASSERTED ABOUT IT.
 *
 * Uniform refinement is nested, so every element of the fine mesh sits inside
 * one coarse element, and the coarse solution restricted to it is a polynomial
 * of degree k -- which the fine element's own P_k space represents exactly.
 * There is no approximation to make. So the transferred field must have the same
 * L2 error against the exact solution as the coarse field it came from, to
 * round-off, and any real number there is the transfer losing something.
 *
 * This is the check that a rate cannot make: an interpolation that quietly lost
 * an order would still produce a usable guess and a converging solve.
 */
BOOST_AUTO_TEST_CASE( theTransferIsExactOntoARefinedMesh )
{
	meq::analytic::ManufacturedNonlinear const eq = equilibrium();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	std::printf( "\n  transfer onto a refined mesh, Example 5\n" );
	std::printf( "    %2s %14s %14s %12s %8s\n",
	             "k", "L2 coarse", "L2 transferred", "relative", "missed" );

	for ( int order = 1; order <= 3; ++order )
	{
		Solved const coarse = solve( order, 8, source, exact );
		double const coarseError = coarse.solver->potentialError( exact );

		// The target: the same degree on a once-refined mesh, so the fine space
		// strictly contains the coarse function.
		mfem::Mesh fine = makeMesh( 8 );
		fine.UniformRefinement();
		meq::GradShafranovSolver fineSolver( fine, order );
		mfem::GridFunction transferred( &fineSolver.potentialSpace() );
		transferred = 0.0;

		meq::FieldTransfer transfer( *coarse.mesh );
		int const missed = transfer.transfer( coarse.solver->potential(), exact,
		                                      transferred );

		double const transferredError = transferred.ComputeL2Error( exact );
		double const relative =
			std::abs( transferredError - coarseError )/coarseError;

		std::printf( "    %2d %14.6e %14.6e %12.2e %8d\n",
		             order, coarseError, transferredError, relative, missed );
		std::fflush( stdout );

		BOOST_TEST( missed == 0,
		            "k = " << order << ": " << missed << " of "
		            << transfer.queried() << " nodes fell outside the coarse mesh, "
		            "which covers the same domain -- the point search is not "
		            "finding elements it should" );
		// MEASURED: 9.6e-7, 3.3e-9, 5.9e-8 at k = 1, 2, 3, through an
		// element-local L2 projection. Nodal interpolation instead gives 9%, 6%
		// and 11%, because meq's L2 spaces are on a GAUSS-LOBATTO basis whose dof
		// points include the element boundary -- where a discontinuous field has
		// two values and the search picks a side. Gauss quadrature points are
		// strictly interior, so the projection never evaluates on a jump. The
		// floor here is not gslib's search tolerance; see WarmStart.cpp.
		BOOST_TEST( relative < 1.0e-5,
		            "k = " << order << ": the transferred field has L2 error "
		            << transferredError << " against " << coarseError
		            << " for the field it came from, a relative difference of "
		            << relative << ". A refined space CONTAINS the coarse "
		            "solution, so this transfer has an exact answer and is not "
		            "producing it" );
	}
}

/*
 * THE SIMPLEST RESTART THERE IS: the same mesh, the same degree, and the answer
 * itself as the guess. DRIVER-PLAN.md section 4's first acceptance bullet, and
 * the one that isolates the seeding from the interpolation -- nothing is
 * interpolated here at all.
 */
BOOST_AUTO_TEST_CASE( anExactRestartFinishesImmediately )
{
	int const order = 3;
	meq::analytic::ManufacturedNonlinear const eq = equilibrium();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	Solved const first = solve( order, 8, source, exact );

	/*
	 * HOW MUCH OF THE ANSWER CAN BE HANDED BACK BEFORE IT BREAKS, which is what
	 * turns "the restart fails" into something actionable. Measured, k = 3,
	 * n = 8, the converged potential scaled and fed back as the guess:
	 *
	 *     0.0  ok    4 iterations      0.4  ok    4
	 *     0.1  ok    4                 0.8  ok    3   <- the guess IS helping
	 *     0.2  ok    4                 1.0  FAIL 30
	 *
	 * So it is not an amplitude limit and not a bad guess: everything up to 80%
	 * of the answer converges, and 80% converges FASTER than cold. Only the
	 * answer itself fails, and it fails at the iteration cap rather than
	 * diverging. That is the signature of a stopping test rather than of an
	 * iteration -- a relative tolerance measured against an initial residual that
	 * is already small has nothing left to ask for -- but it is a signature and
	 * not yet a diagnosis, and the residual histories are where to look next.
	 */
	for ( double scale : { 0.0, 0.1, 0.2, 0.4, 0.8, 1.0 } )
	{
		mfem::Mesh probeMesh = makeMesh( 8 );
		meq::GradShafranovSolver probe( probeMesh, order );
		mfem::GridFunction scaled( first.solver->potential() );
		scaled *= scale;
		probe.setSource( source );
		probe.setBoundaryData( exact );
		probe.setInitialGuess( scaled );
		bool ok = true;
		try { probe.solve(); } catch ( std::exception const & ) { ok = false; }
		std::printf( "    guess = %.1f x the answer: %s, %d iterations\n", scale,
		             ok ? "ok " : "FAIL", probe.newtonIterations() );
		std::fflush( stdout );
	}

	mfem::Mesh again = makeMesh( 8 );
	meq::GradShafranovSolver restarted( again, order );
	restarted.setSource( source );
	restarted.setBoundaryData( exact );
	restarted.setInitialGuess( first.solver->potential() );
	restarted.solve();

	std::printf( "\n  exact restart, k = %d: cold %d iterations, restarted %d\n",
	             order, first.solver->newtonIterations(),
	             restarted.newtonIterations() );
	std::fflush( stdout );

	BOOST_TEST( restarted.newtonIterations() <= 2,
	            "restarting from the converged answer took "
	            << restarted.newtonIterations() << " Newton iterations. Handed its "
	            "own solution the iteration should have nothing to do" );
	BOOST_TEST( restarted.newtonIterations() < first.solver->newtonIterations(),
	            "restarting from the converged answer took as many iterations as "
	            "starting cold, so the guess is not reaching the iterate" );
}

/*
 * AND IT IS WORTH SOMETHING -- IN ITERATIONS. THE FIRST RESIDUAL IS NO LONGER
 * THE INSTRUMENT, AND WHY NOT IS A PROPERTY OF NonlinearOrdering::NPC.
 *
 * This test used to assert || r_0 || fell by 10x, on the reasoning that Newton's
 * iteration count is coarse -- it moves in whole steps and a good guess often
 * buys less than one -- while || r_0 || is the quantity a guess actually
 * changes. That was true while the unknown was the TRACE ALONE: q and psi were
 * functions of it, so seeding psi and the trace seeded everything there was.
 *
 * UNDER NPC q IS AN UNKNOWN OF ITS OWN, AND meq DOES NOT SEED IT. prepare()
 * projects the guess onto the potential and the trace and leaves the flux block
 * at zero, deliberately: a guess for psi says nothing about q without
 * differentiating it, and the guess arrives as a bare mfem::Coefficient, which
 * cannot be differentiated. So the guessed state is INCONSISTENT in exactly the
 * row that couples them -- the flux equation q - (1/r) grad psi = 0 reads worst
 * when psi is the converged answer and q is zero -- and || r_0 || goes UP.
 * Measured here at k = 3: cold 1.771e-01, warm 2.638e-01, a factor of 1.5 the
 * wrong way.
 *
 * THE GUESS IS STILL DOING ITS JOB, WHICH IS WHY THIS IS A CHANGE OF INSTRUMENT
 * AND NOT A REGRESSION: the warm solve takes 2 Newton iterations against 4 cold
 * and lands on the same L2 error to every figure printed. The flux row is LINEAR
 * in q, so one Newton step recovers the q that belongs to the guessed psi and
 * the iteration proceeds from a genuinely warm state -- which is also why an
 * exact restart above still finishes in 1.
 *
 * THE FIX, IF THE STRONGER PROPERTY IS EVER WANTED, is to seed the flux:
 * darcyFlux = -(1/r) grad psi_guess, via GradientGridFunctionCoefficient, for
 * the setInitialGuess( GridFunction const & ) overload that has a differentiable
 * guess to work with. That is a small change with its own measurement to make,
 * and it is not done.
 *
 * The converged answer must not move. That is what separates a starting point
 * from data, and it is the same invariance the +5% Jacobian experiment relied
 * on: a guess that changed the answer would be entering as a boundary condition
 * or a source, not as a guess.
 */
BOOST_AUTO_TEST_CASE( aWarmStartCutsTheWorkAndNotTheAnswer )
{
	int const order = 3;
	meq::analytic::ManufacturedNonlinear const eq = equilibrium();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	Solved const coarse = solve( order, 8, source, exact );

	mfem::Mesh coldMesh = refinedMesh( 8 );
	mfem::Mesh warmMesh = refinedMesh( 8 );

	// Cold.
	meq::GradShafranovSolver cold( coldMesh, order );
	cold.setSource( source );
	cold.setBoundaryData( exact );
	cold.solve();

	// Warm, at full order.
	meq::GradShafranovSolver warm( warmMesh, order );
	mfem::GridFunction guess( &warm.potentialSpace() );
	guess = 0.0;
	meq::FieldTransfer transfer( *coarse.mesh );
	transfer.transfer( coarse.solver->potential(), exact, guess );

	warm.setSource( source );
	warm.setBoundaryData( exact );
	warm.setInitialGuess( guess );
	warm.solve();

	double const coldFirst = cold.newtonResiduals().front();
	double const warmFirst = warm.newtonResiduals().front();
	double const coldError = cold.potentialError( exact );
	double const warmError = warm.potentialError( exact );

	std::printf( "\n  warm start from a coarse solve, k = %d\n", order );
	std::printf( "    cold: || r_0 || = %.6e, %d iterations, L2 %.6e\n",
	             coldFirst, cold.newtonIterations(), coldError );
	std::printf( "    warm: || r_0 || = %.6e, %d iterations, L2 %.6e\n",
	             warmFirst, warm.newtonIterations(), warmError );
	std::printf( "    the guess moved the first residual by %.2fx "
	             "(NOT the instrument here -- see above)\n",
	             coldFirst/warmFirst );
	std::fflush( stdout );

	// STRICTLY fewer, not merely no more. Under NPC this is the whole of the
	// head start a psi-only guess can show, so a warm start that only matched
	// the cold iteration count would mean the guess was not reaching the solve
	// at all -- which is the failure this is here to catch, and it is the same
	// failure the old || r_0 || assertion was catching by another route.
	BOOST_TEST( warm.newtonIterations() < cold.newtonIterations(),
	            "the warm start took " << warm.newtonIterations()
	            << " Newton iterations against " << cold.newtonIterations()
	            << " cold, so the guess bought nothing. Under NPC the first "
	            "residual is not the instrument -- the flux block is unseeded, so "
	            "|| r_0 || can rise while the solve still starts warm -- and the "
	            "iteration count is what is left. If THAT has stopped moving, the "
	            "guess is not reaching the iterate" );

	// Six figures, per DRIVER-PLAN.md section 4's acceptance.
	BOOST_TEST( std::abs( warmError - coldError ) < 1.0e-6*coldError,
	            "the warm start converged to L2 error " << warmError
	            << " where the cold start reached " << coldError
	            << ". A starting point must move the path and not the answer" );
}

/*
 * THE MEASUREMENT THAT JUSTIFIES THE GSLIB DEPENDENCY RATHER THAN ASSERTING IT.
 *
 * DRIVER-PLAN.md section 4: the structured-grid route is SECOND ORDER in the
 * grid spacing whatever the solve's degree, so restarting a k = 3 solve through
 * it discards most of what the solve computed. Both routes converge, which is
 * exactly why this needs measuring -- the cheap one does not fail, it just
 * quietly does less.
 *
 * The comparison is deliberately generous to the grid: 257 x 257 nodes, which is
 * far finer than the k = 3 mesh the field came from and far finer than anything a
 * plotting file would carry. If the full-order route did not win here it would
 * not be worth a dependency.
 */
BOOST_AUTO_TEST_CASE( fullOrderCarriesMoreThanAStructuredGrid )
{
	int const order = 3;
	int const gridNodes = 257;

	meq::analytic::ManufacturedNonlinear const eq = equilibrium();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	Solved const coarse = solve( order, 8, source, exact );
	double const coarseError = coarse.solver->potentialError( exact );

	mfem::Mesh fullMesh = refinedMesh( 8 );
	mfem::Mesh griddedMesh = refinedMesh( 8 );

	// Route one: full order, never leaving the finite element representation.
	meq::GradShafranovSolver full( fullMesh, order );
	mfem::GridFunction fullGuess( &full.potentialSpace() );
	fullGuess = 0.0;
	{
		meq::FieldTransfer transfer( *coarse.mesh );
		transfer.transfer( coarse.solver->potential(), exact, fullGuess );
	}
	double const fullGuessError = fullGuess.ComputeL2Error( exact );

	// Route two: out to a structured grid and back by bilinear interpolation,
	// over a SEQUENCE of grids -- because a single grid size cannot say what the
	// route costs. Bilinear interpolation is second order in the GRID spacing,
	// so its accuracy is a property of the file, not of the solve that wrote it.
	std::printf( "\n  full order against a structured grid, k = %d, mesh h = %.4f\n",
	             order, ( rMax - rMin )/8.0 );
	std::printf( "    %10s %10s %14s %8s\n", "grid", "spacing", "L2 of guess", "rate" );
	std::printf( "    %10s %10s %14.6e %8s   <- the coarse solve itself\n",
	             "-", "-", coarseError, "-" );

	std::vector<double> gridErrors;
	std::vector<double> spacings;
	for ( int gridNodes : { 65, 129, 257 } )
	{
		meq::GridSampler sampler( *coarse.mesh, rMin, rMax, gridNodes,
		                          zMin, zMax, gridNodes );
		std::vector<double> values;
		sampler.sample( coarse.solver->potential(), values,
		                std::numeric_limits<double>::quiet_NaN() );
		BilinearGridCoefficient grid( sampler, values, 0.0 );

		mfem::Mesh probeMesh = refinedMesh( 8 );
		meq::GradShafranovSolver probe( probeMesh, order );
		mfem::GridFunction gridGuess( &probe.potentialSpace() );
		gridGuess.ProjectCoefficient( grid );

		double const spacing = ( rMax - rMin )/( gridNodes - 1 );
		double const error = gridGuess.ComputeL2Error( exact );
		spacings.push_back( spacing );
		gridErrors.push_back( error );

		std::printf( "    %10d %10.5f %14.6e", gridNodes, spacing, error );
		if ( gridErrors.size() > 1 )
		{
			std::size_t const i = gridErrors.size() - 1;
			std::printf( " %8.3f\n",
			             std::log( gridErrors[ i - 1 ]/gridErrors[ i ] )
			                 /std::log( spacings[ i - 1 ]/spacings[ i ] ) );
		}
		else
		{
			std::printf( " %8s\n", "-" );
		}
		std::fflush( stdout );
	}
	std::printf( "    %10s %10s %14.6e %8s   <- full order, no grid involved\n",
	             "-", "-", fullGuessError, "-" );
	std::fflush( stdout );

	// And what the full-order guess is worth as a starting point.
	full.setSource( source );
	full.setBoundaryData( exact );
	full.setInitialGuess( fullGuess );
	full.solve();

	std::printf( "    full order as a start: || r_0 || %.6e, %d Newton iterations\n",
	             full.newtonResiduals().front(), full.newtonIterations() );
	std::fflush( stdout );

	// THE FULL-ORDER GUESS IS THE COARSE SOLUTION, so its error against the exact
	// solution is the coarse discretisation error and nothing more. It has no
	// grid to depend on, which is the whole of its advantage.
	BOOST_TEST( std::abs( fullGuessError - coarseError ) < 1.0e-5*coarseError,
	            "the full-order guess has L2 error " << fullGuessError
	            << " where the solve it came from had " << coarseError
	            << ": the transfer is losing accuracy the representation already had" );

	// AND THE GRID ROUTE IS SECOND ORDER IN THE GRID SPACING, whatever degree the
	// solve was. That is DRIVER-PLAN.md section 4's claim and it is the whole
	// argument for the dependency, so it is measured rather than repeated.
	// MEASURED: 1.813 and 1.735 over 65 -> 129 -> 257 nodes, drifting DOWN from
	// two rather than up to it -- because by 257 nodes the grid error (1.01e-5)
	// is within a factor of two of the coarse solve's own (6.00e-6), which the
	// grid cannot represent better than. The sequence is leaving the asymptotic
	// regime from below, so the bound allows for that rather than demanding a
	// clean 2.
	for ( std::size_t i = 1; i < gridErrors.size(); ++i )
	{
		double const rate = std::log( gridErrors[ i - 1 ]/gridErrors[ i ] )
		                    /std::log( spacings[ i - 1 ]/spacings[ i ] );
		BOOST_TEST( rate > 1.7,
		            "the structured-grid guess converged at " << rate
		            << " in the grid spacing, not the second order bilinear "
		            "interpolation gives. If it is now higher, the interchange "
		            "route has stopped being the second-order thing "
		            "DRIVER-PLAN.md section 4 argues against" );
		BOOST_TEST( rate < 2.3,
		            "the structured-grid guess converged at " << rate
		            << ", above second order, which bilinear interpolation cannot "
		            "do -- check the sampler before believing it" );
	}

	// AND THAT IS WHY IT IS NOT A SUBSTITUTE: its accuracy is a property of the
	// FILE and the full-order route's is a property of the SOLVE. Here they are
	// within a factor of two, because a 257 x 257 grid over this box has spacing
	// 0.003 against a mesh h of 0.1 -- fine enough to keep up with k = 3 on eight
	// cells. Refine the mesh or raise the degree and the grid must be refined
	// with it, quadratically, to keep pace; the full-order route never has to be.
	// The claim under test is the SCALING above, not this ratio.
	BOOST_TEST( gridErrors.back() > fullGuessError,
	            "the structured-grid guess (" << gridErrors.back()
	            << ") is at least as accurate as the full-order one ("
	            << fullGuessError << ") even at 257 nodes, so nothing here "
	            "justifies the GSLIB dependency" );

	// The full-order guess must be worth what it claims as a STARTING POINT too,
	// which is what a restart is actually for. MEASURED: || r_0 || falls from
	// 1.57e+01 cold to 1.80e-03, and Newton finishes in one step.
	BOOST_TEST( full.newtonIterations() <= 2,
	            "restarting at full order from a converged coarse solve took "
	            << full.newtonIterations() << " Newton iterations. Moving to an "
	            "adjacent equilibrium is what this is for, and it should cost one "
	            "or two steps" );
}
