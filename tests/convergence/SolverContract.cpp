#define BOOST_TEST_MODULE SolverContract
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/ManufacturedNonlinear.hpp"
#include "analytic/Soloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * THE SOLVER'S CONTRACT, WHICH THE RATE STUDIES DO NOT TOUCH.
 *
 * Everything else in tests/convergence/ drives GradShafranovSolver down its happy
 * path and measures the answer. That leaves most of the class untested: measured
 * with gcovr, GradShafranov.cpp sat at 60% of lines against 92-100% for the
 * configuration layer, and 22 of its 51 methods were called by nothing at all.
 * The best-covered code in meq parsed TOML and the least-covered was the part
 * every physical claim rests on.
 *
 * WHAT IS TESTED HERE IS NOT COVERAGE FOR ITS OWN SAKE. Three kinds of thing,
 * each of which has already gone wrong once in this project:
 *
 *   - SETTERS THAT DO NOT TAKE. setGlobalisation() reset `prepared` but not
 *     `built`, so switching a live solver between a Picard path and a Newton one
 *     silently reused the other path's blocks. That was found by accident while
 *     wiring PicardThenNewton. A round-trip through the getter would not have
 *     caught that particular one, but the class of bug is why the setters are
 *     exercised rather than assumed.
 *   - THE REFUSALS. Twenty-eight throws state what the class will not do -- set a
 *     source after the forms are built, hand out the reduced operator before
 *     prepare(), post-process before solving. An untested refusal is one that
 *     drifts into a segfault or, worse, into silently doing the thing.
 *   - THE QUANTITIES NOTHING READ. trace(), totalFlux() and postProcessedFlux()
 *     are computed on every solve and were checked by nothing. Two of them carry
 *     DarcyForm's sign convention, which is the single most expensive convention
 *     in this tree to get wrong.
 */

namespace
{
	using meq::tests::EquilibriumSource;

	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}

	mfem::Mesh makeMesh( int n )
	{
		return meq::tests::makeMesh( box(), n );
	}
}

/*
 * THE CONSTRUCTOR REFUSES WHAT IT CANNOT SOLVE, and reports what it was given.
 */
BOOST_AUTO_TEST_CASE( theConstructorRejectsWhatItCannotSolve )
{
	mfem::Mesh mesh = makeMesh( 4 );

	BOOST_CHECK_THROW( meq::GradShafranovSolver( mesh, -1 ), std::invalid_argument );

	// One dimension is not ( r, z ), and neither is three. A 3D mesh reaching the
	// assembly would produce a shape mismatch a long way from the cause.
	mfem::Mesh line = mfem::Mesh::MakeCartesian1D( 4 );
	BOOST_CHECK_THROW( meq::GradShafranovSolver( line, 1 ), std::invalid_argument );

	mfem::Mesh cube = mfem::Mesh::MakeCartesian3D( 2, 2, 2, mfem::Element::TETRAHEDRON );
	BOOST_CHECK_THROW( meq::GradShafranovSolver( cube, 1 ), std::invalid_argument );

	// Degree zero is legal -- P_0 flux and potential with a P_0 trace is a finite
	// volume scheme and hybridization does not mind -- so it must NOT throw.
	BOOST_CHECK_NO_THROW( meq::GradShafranovSolver( mesh, 0 ) );

	meq::GradShafranovSolver solver( mesh, 3, 2.5 );
	BOOST_TEST( solver.order() == 3, "order() reported " << solver.order() );
	BOOST_TEST( solver.tau() == 2.5, "tau() reported " << solver.tau() );

	// The default tau is 1, which both papers use and which CLAUDE.md records the
	// pre-port code having quietly set to 5 with no measurement behind it.
	meq::GradShafranovSolver defaulted( mesh, 1 );
	BOOST_TEST( defaulted.tau() == 1.0,
	            "the default tau is " << defaulted.tau() << ", not 1. The pre-port "
	            "code used 5 for no recorded reason; do not reinstate it without a "
	            "measurement" );
}

/*
 * THE CONFIGURATION ROUND-TRIPS, and the ones with a valid range enforce it.
 */
BOOST_AUTO_TEST_CASE( theConfigurationRoundTrips )
{
	using G = meq::GradShafranovSolver::Globalisation;
	using N = meq::GradShafranovSolver::NonlinearOrdering;
	using L = meq::GradShafranovSolver::LocalSolver;

	mfem::Mesh mesh = makeMesh( 4 );
	meq::GradShafranovSolver solver( mesh, 1 );

	// Compared as integers because Boost.Test insists on printing whatever it
	// compares and a scoped enum has no operator<<.
	BOOST_TEST( static_cast<int>( solver.globalisation() ) == static_cast<int>( G::None ),
	            "the default globalisation is not None" );
	BOOST_TEST( static_cast<int>( solver.nonlinearOrdering() )
	                == static_cast<int>( N::CondenseThenLinearise ),
	            "the default ordering is not CondenseThenLinearise" );

	for ( G choice : { G::PicardOnly, G::AndersonPicard, G::PicardThenNewton, G::None } )
	{
		solver.setGlobalisation( choice );
		BOOST_TEST( static_cast<int>( solver.globalisation() ) == static_cast<int>( choice ),
		            "setGlobalisation did not take for choice "
		            << static_cast<int>( choice ) );
	}

	for ( N choice : { N::LineariseThenCondense, N::CondenseThenLinearise } )
	{
		solver.setNonlinearOrdering( choice );
		BOOST_TEST( static_cast<int>( solver.nonlinearOrdering() )
		                == static_cast<int>( choice ),
		            "setNonlinearOrdering did not take for choice "
		            << static_cast<int>( choice ) );
	}

	// No getter for these two, so what is checked is the range they enforce.
	BOOST_CHECK_NO_THROW( solver.setLocalSolver( L::Lbfgs ) );
	BOOST_CHECK_NO_THROW( solver.setLocalSolver( L::Newton ) );

	// Damping in ( 0, 1 ]: one is undamped, which is the Anderson default, and
	// zero would not move at all.
	// invalid_argument, not the bare logic_error these two used to throw: the
	// constructor already distinguishes a bad VALUE from a bad CALL ORDER and
	// these did not. invalid_argument derives from logic_error, so a caller
	// catching the wider type is unaffected.
	BOOST_CHECK_NO_THROW( solver.setPicardDamping( 1.0 ) );
	BOOST_CHECK_NO_THROW( solver.setPicardDamping( 0.5 ) );
	BOOST_CHECK_THROW( solver.setPicardDamping( 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( solver.setPicardDamping( 1.5 ), std::invalid_argument );
	BOOST_CHECK_THROW( solver.setPicardDamping( -1.0 ), std::invalid_argument );

	// Depth zero is plain Picard with no acceleration, which is meaningful.
	BOOST_CHECK_NO_THROW( solver.setAndersonDepth( 0 ) );
	BOOST_CHECK_NO_THROW( solver.setAndersonDepth( 3 ) );
	BOOST_CHECK_THROW( solver.setAndersonDepth( -1 ), std::invalid_argument );
}

/*
 * THE GUESS IS REPORTED AND CAN BE TAKEN BACK.
 *
 * clearInitialGuess() matters more than it looks since solve() started scaling
 * its convergence target to the COLD residual when a guess is present -- see
 * CLAUDE.md's trap on a target scaled to || r_0 ||. A clear that left the guess
 * behind would leave the target scaled too, and the only symptom would be a
 * solve that stops slightly too early.
 */
BOOST_AUTO_TEST_CASE( theGuessCanBeSetAndCleared )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	mfem::Mesh mesh = makeMesh( 8 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( exact );

	BOOST_TEST( !solver.hasInitialGuess(), "a fresh solver reports a guess" );
	solver.setInitialGuess( exact );
	BOOST_TEST( solver.hasInitialGuess(), "setInitialGuess did not take" );
	solver.clearInitialGuess();
	BOOST_TEST( !solver.hasInitialGuess(), "clearInitialGuess did not take" );

	// And the cleared solver behaves as a cold one. Measured against a solver
	// that never had a guess, on the same problem: identical iteration counts and
	// the same answer.
	solver.solve();

	mfem::Mesh reference = makeMesh( 8 );
	meq::GradShafranovSolver cold( reference, 2 );
	cold.setSource( source );
	cold.setBoundaryData( exact );
	cold.solve();

	std::printf( "\n  cleared %d iterations, never-guessed %d iterations\n",
	             solver.newtonIterations(), cold.newtonIterations() );
	std::fflush( stdout );

	BOOST_TEST( solver.newtonIterations() == cold.newtonIterations(),
	            "a solver whose guess was cleared took " << solver.newtonIterations()
	            << " iterations against " << cold.newtonIterations()
	            << " for one that never had one, so clearing left something behind "
	            "-- most likely the convergence target still scaled to the cold "
	            "residual" );
	BOOST_TEST( std::abs( solver.potentialError( exact ) - cold.potentialError( exact ) )
	            < 1.0e-12*cold.potentialError( exact ),
	            "the two disagree on the answer" );
}

/*
 * THE REFUSALS. Each of these states something the class will not do, and an
 * untested refusal drifts into a segfault or into silently doing the thing.
 */
BOOST_AUTO_TEST_CASE( theSolverRefusesToWorkOutOfOrder )
{
	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	meq::SolovievSource const source( eq.getA() );
	mfem::ConstantCoefficient zero( 0.0 );

	// Nothing may be had from the reduced system before prepare().
	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		BOOST_CHECK_THROW( solver.reducedOperator(), std::logic_error );
		BOOST_CHECK_THROW( solver.reducedRhs(), std::logic_error );
		BOOST_CHECK_THROW( solver.reducedSolution(), std::logic_error );
		BOOST_CHECK_THROW( solver.essentialTraceDofs(), std::logic_error );
	}

	// prepare() needs a source and boundary data, and says which is missing.
	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		BOOST_CHECK_THROW( solver.prepare(), std::logic_error );

		solver.setSource( source );
		BOOST_CHECK_THROW( solver.prepare(), std::logic_error );

		solver.setBoundaryData( zero );
		BOOST_CHECK_NO_THROW( solver.prepare() );
	}

	// post-processing before a solve, and reading it before post-processing.
	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		solver.setSource( source );
		solver.setBoundaryData( zero );

		BOOST_CHECK_THROW( solver.postProcess(), std::logic_error );
		mfem::ConstantCoefficient nought( 0.0 );
		BOOST_CHECK_THROW( solver.postProcessedPotentialError( nought ), std::logic_error );

		solver.solve();
		BOOST_TEST( !solver.isPostProcessed(),
		            "a fresh solve reports itself post-processed" );
		BOOST_CHECK_NO_THROW( solver.postProcess() );
		BOOST_TEST( solver.isPostProcessed(),
		            "postProcess() ran and did not say so" );

		// A new solve invalidates the old post-processing rather than pairing the
		// two silently.
		solver.solve();
		BOOST_TEST( !solver.isPostProcessed(),
		            "a second solve left the previous post-processing marked valid, "
		            "so psi* and psi_h can be read as a matched pair when they are "
		            "not" );
	}
}

/*
 * ONE SOLVER HOLDS ONE PROBLEM, and says so rather than quietly taking the
 * second thing it was handed.
 */
BOOST_AUTO_TEST_CASE( theSolverRefusesASecondProblem )
{
	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	meq::SolovievSource const source( eq.getA() );
	mfem::ConstantCoefficient zero( 0.0 );
	mfem::ConstantCoefficient other( 1.0 );

	mfem::Mesh mesh = makeMesh( 4 );
	meq::GradShafranovSolver solver( mesh, 1 );
	solver.setSource( source );

	// A second source of either kind, before anything is built. The same-kind
	// case is the one that used to be allowed: the checks were symmetric across
	// the two overloads and silent within each, so setSource( a ); setSource( b )
	// replaced a with b and the solve answered a different question.
	BOOST_CHECK_THROW( solver.setSource( source ), std::logic_error );
	BOOST_CHECK_THROW( solver.setSource( other ), std::logic_error );

	solver.setBoundaryData( zero );
	solver.solve();

	// And after the forms are built, no source at all.
	BOOST_CHECK_THROW( solver.setSource( other ), std::logic_error );
}

/*
 * THE ESSENTIAL TRACE CONDITION REALLY IMPOSES THE DATUM, ON BOTH PATHS.
 *
 * CLAUDE.md flags this as the combination with no regression behind it:
 * DarcyHybridization::SetEssentialBC together with a NON-LINEAR reduced operator
 * was broken on the MFEM branch until recently -- EliminateTraceTrueDofsInRHS
 * returned early for non-linear problems and the essential condition was
 * silently ignored, so the solver answered a different problem -- and that file
 * says outright that if a converged answer ever looks wrong near Gamma, look
 * here before looking at meq's assembly. Nothing looked.
 *
 * A CONSTANT datum makes it exactly checkable. A constant lies in every
 * polynomial trace space, so its projection is itself and every essential trace
 * dof must equal it to round-off -- no discretisation error to allow for, no
 * tolerance to argue about. Anything else means the condition was not imposed.
 */
BOOST_AUTO_TEST_CASE( theEssentialTraceConditionImposesTheDatum )
{
	double const datumValue = 0.37;
	mfem::ConstantCoefficient datum( datumValue );

	meq::analytic::SolovievEquilibrium const soloviev
		= meq::analytic::SolovievEquilibrium::nstx();
	meq::SolovievSource const nonlinear( soloviev.getA() );
	mfem::FunctionCoefficient linear( [ &soloviev ]( mfem::Vector const &x )
	{
		return soloviev.f( x( 0 ), x( 1 ), 0.0 );
	} );

	for ( int path = 0; path < 2; ++path )
	{
		for ( int order = 1; order <= 3; ++order )
		{
			mfem::Mesh mesh = makeMesh( 6 );
			meq::GradShafranovSolver solver( mesh, order );
			if ( path == 0 )
				solver.setSource( linear );
			else
				solver.setSource( nonlinear );
			solver.setBoundaryData( datum );
			solver.solve();

			mfem::GridFunction const &trace = solver.trace();
			mfem::Array<int> const &essential = solver.essentialTraceDofs();

			BOOST_TEST_REQUIRE( essential.Size() > 0,
			                    "no trace dof is essential, so the Dirichlet "
			                    "condition is not being imposed anywhere at all" );

			double worst = 0.0;
			for ( int i = 0; i < essential.Size(); ++i )
				worst = std::max( worst,
				                  std::abs( trace( essential[ i ] ) - datumValue ) );

			// MEASURED: 0.000e+00 at every order on both paths. Not "small" --
			// exactly zero, which is what a constant projected onto a polynomial
			// trace space and then imposed as an essential condition has to give.
			std::printf( "    %s path, k = %d: %d essential trace dofs, worst "
			             "departure from the datum %.3e\n",
			             path == 0 ? "linear " : "Newton ", order,
			             essential.Size(), worst );
			std::fflush( stdout );

			BOOST_TEST( worst < 1.0e-12,
			            ( path == 0 ? "linear" : "Newton" ) << " path, k = " << order
			            << ": an essential trace dof sits " << worst
			            << " from the datum. A constant projects exactly onto the "
			            "trace space, so this should be round-off. See CLAUDE.md on "
			            "DarcyHybridization::SetEssentialBC with a non-linear "
			            "reduced operator" );
		}
	}
}

/*
 * AND THE QUANTITIES NOTHING READ. trace(), totalFlux() and postProcessedFlux()
 * are computed on every solve and post-process and were checked by nothing.
 *
 * Two of them carry DarcyForm's sign convention, which CLAUDE.md records as the
 * most expensive convention in this tree to get wrong -- flux() negates what
 * DarcyForm holds and totalFlux() deliberately does NOT, because it is defined
 * by the constraint equation which is written in DarcyForm's convention
 * throughout. A test that never looks at totalFlux() cannot notice that drifting.
 */
BOOST_AUTO_TEST_CASE( theReconstructedQuantitiesAreThereAndSigned )
{
	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );
	mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
	                                                       mfem::Vector &value )
	{
		eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
	} );

	std::printf( "\n  the reconstructed quantities, Solov'ev on the linear path\n" );
	std::printf( "    %2s %14s %14s %14s\n", "k", "L2( q )", "L2( q* )", "|q + total|" );

	for ( int order = 1; order <= 3; ++order )
	{
		mfem::Mesh mesh = makeMesh( 8 );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( sourceCoeff );
		solver.setBoundaryData( psiCoeff );
		solver.solve();
		solver.postProcess();

		double const fluxError = solver.fluxError( fluxCoeff );
		double const enrichedError =
			solver.postProcessedFlux().ComputeL2Error( fluxCoeff );

		// totalFlux() is in DarcyForm's convention and flux() has been negated out
		// of it, so the two should be opposites. They live on different spaces --
		// the total flux is enriched -- so this is a norm comparison rather than a
		// dof one: || q + q_total || against || q ||.
		mfem::VectorGridFunctionCoefficient totalCoeff( &solver.totalFlux() );
		mfem::GridFunction sum( solver.flux() );
		mfem::Vector zeroes( 2 );
		zeroes = 0.0;
		mfem::VectorConstantCoefficient zeroVector( zeroes );
		double const fluxNorm = solver.flux().ComputeL2Error( zeroVector );

		// q evaluated against -q_total: if the conventions agree these cancel.
		mfem::GridFunction negatedTotal( solver.totalFlux() );
		negatedTotal.Neg();
		mfem::VectorGridFunctionCoefficient negatedCoeff( &negatedTotal );
		double const mismatch = solver.flux().ComputeL2Error( negatedCoeff )/fluxNorm;

		std::printf( "    %2d %14.6e %14.6e %14.6e\n",
		             order, fluxError, enrichedError, mismatch );
		std::fflush( stdout );

		BOOST_TEST( std::isfinite( enrichedError ),
		            "k = " << order << ": the post-processed flux is not finite" );
		// MEASURED: 6.92e-4 against 9.58e-4, 8.90e-6 against 1.58e-5, 1.43e-7
		// against 2.22e-7 at k = 1, 2, 3 -- about 1.5x better at every order. So
		// "better than q" is the assertion, not "not much worse".
		BOOST_TEST( enrichedError < fluxError,
		            "k = " << order << ": the post-processed flux has L2 error "
		            << enrichedError << " against " << fluxError << " for q itself, "
		            "so the enriched reconstruction is not adding anything" );
		// MEASURED: 2.53e-3, 4.92e-5, 5.38e-7 -- converging, which is what says the
		// two conventions agree rather than happening to be close on one mesh.
		BOOST_TEST( mismatch < 0.05,
		            "k = " << order << ": || q - ( -q_total ) || / || q || is "
		            << mismatch << ". flux() negates DarcyForm's convention and "
		            "totalFlux() deliberately does not, so the two must be "
		            "opposites; if one of those negations has moved, this is where "
		            "it shows" );
	}
}
