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
 * AND IT IS WORTH SOMETHING: THE FIRST RESIDUAL, WHICH IS WHERE A WARM START
 * SHOWS UP.
 *
 * Newton's iteration count is a coarse instrument -- it moves in whole steps and
 * a good guess often buys less than one. || r_0 || is the quantity a guess
 * actually changes, and it changes it by orders.
 *
 * The converged answer must not move. That is what separates a starting point
 * from data, and it is the same invariance the +5% Jacobian experiment relied
 * on: a guess that changed the answer would be entering as a boundary condition
 * or a source, not as a guess.
 */
BOOST_AUTO_TEST_CASE( aWarmStartCutsTheFirstResidualAndNotTheAnswer )
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
	std::printf( "    the guess cut the first residual by %.1fx\n",
	             coldFirst/warmFirst );
	std::fflush( stdout );

	BOOST_TEST( warmFirst < 0.1*coldFirst,
	            "the warm start's first residual is " << warmFirst
	            << " against " << coldFirst << " cold, which is not the head start "
	            "a full-order guess from a converged solve should give" );

	BOOST_TEST( warm.newtonIterations() <= cold.newtonIterations(),
	            "the warm start took " << warm.newtonIterations()
	            << " Newton iterations against " << cold.newtonIterations()
	            << " cold" );

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

	// Route two: out to a structured grid and back by bilinear interpolation.
	meq::GridSampler sampler( *coarse.mesh, rMin, rMax, gridNodes,
	                          zMin, zMax, gridNodes );
	std::vector<double> values;
	sampler.sample( coarse.solver->potential(), values,
	                std::numeric_limits<double>::quiet_NaN() );
	BilinearGridCoefficient grid( sampler, values, 0.0 );

	meq::GradShafranovSolver gridded( griddedMesh, order );
	mfem::GridFunction gridGuess( &gridded.potentialSpace() );
	gridGuess.ProjectCoefficient( grid );
	double const gridGuessError = gridGuess.ComputeL2Error( exact );

	// And what each is worth as a starting point.
	full.setSource( source );
	full.setBoundaryData( exact );
	full.setInitialGuess( fullGuess );
	full.solve();

	gridded.setSource( source );
	gridded.setBoundaryData( exact );
	gridded.setInitialGuess( gridGuess );
	gridded.solve();

	std::printf( "\n  full order against a %d x %d structured grid, k = %d\n",
	             gridNodes, gridNodes, order );
	std::printf( "    the coarse solve itself       L2 %.6e\n", coarseError );
	std::printf( "    full-order transfer           L2 %.6e   || r_0 || %.6e\n",
	             fullGuessError, full.newtonResiduals().front() );
	std::printf( "    through the grid, bilinear    L2 %.6e   || r_0 || %.6e\n",
	             gridGuessError, gridded.newtonResiduals().front() );
	std::printf( "    the grid route is %.1fx worse as a guess\n",
	             gridGuessError/fullGuessError );
	std::fflush( stdout );

	// THE FULL-ORDER GUESS IS THE COARSE SOLUTION, so its error against the exact
	// solution is the coarse discretisation error and nothing more.
	BOOST_TEST( std::abs( fullGuessError - coarseError ) < 1.0e-5*coarseError,
	            "the full-order guess has L2 error " << fullGuessError
	            << " where the solve it came from had " << coarseError
	            << ": the transfer is losing accuracy that the representation "
	            "already had" );

	// AND THE GRID ROUTE IS MEASURABLY WORSE, on a grid far finer than the mesh.
	BOOST_TEST( gridGuessError > 10.0*fullGuessError,
	            "the structured-grid guess has L2 error " << gridGuessError
	            << " against " << fullGuessError << " for the full-order one, so "
	            "on this problem the grid route loses little and the GSLIB "
	            "dependency is not buying what DRIVER-PLAN.md section 4 claims. "
	            "Check the grid resolution before concluding that -- this is "
	            "deliberately a generous grid" );
}

/*
 * A RESTART ONTO A DOMAIN THE STORED ONE DOES NOT COVER MUST SAY SO.
 *
 * The guess is only a guess, so the fallback need not be clever -- but a restart
 * that found no data for most of the domain has quietly become a cold start, and
 * the iteration count will not obviously say so. DRIVER-PLAN.md section 4 asks
 * for the count to be reported, and this is what makes the count true.
 */
BOOST_AUTO_TEST_CASE( nodesOutsideTheStoredMeshAreCountedAndFallBack )
{
	int const order = 2;
	meq::analytic::ManufacturedNonlinear const eq = equilibrium();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	Solved const stored = solve( order, 8, source, exact );

	// A target box reaching well beyond the stored one in r and z, so that a
	// known fraction of it has no data.
	meq::tests::Rectangle const wider{ rMin, rMax + 0.8, zMin, zMax + 0.6 };
	mfem::Mesh target = meq::tests::makeMesh( wider, 8 );

	meq::GradShafranovSolver solver( target, order );
	mfem::GridFunction guess( &solver.potentialSpace() );
	guess = 0.0;

	mfem::ConstantCoefficient fallback( -7.0 );
	meq::FieldTransfer transfer( *stored.mesh );
	int const missed = transfer.transfer( stored.solver->potential(), fallback,
	                                      guess );

	std::printf( "\n  a target reaching outside the stored mesh\n" );
	std::printf( "    %d of %d nodes had no data (%.1f%%), worst search "
	             "distance %.3e\n",
	             missed, transfer.queried(),
	             100.0*missed/transfer.queried(), transfer.worstDistance() );
	std::fflush( stdout );

	BOOST_TEST( missed > 0,
	            "no node fell outside the stored mesh, on a target box that "
	            "extends 0.8 beyond it in r and 0.6 in z -- the miss count is not "
	            "counting" );
	BOOST_TEST( missed < transfer.queried(),
	            "every node missed, so the transfer found nothing at all rather "
	            "than partially covering the target" );

	// The fallback really was used, and is not a zero that happens to look like
	// one. -7 is nowhere in this solution.
	BOOST_TEST( guess.Min() < -6.0,
	            "the fallback coefficient's value never appears in the guess, so "
	            "the missed nodes were left at whatever the target held rather "
	            "than filled from it" );
}
