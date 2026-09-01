#define BOOST_TEST_MODULE SolverContract
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
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
	// NPC is the default, and CondenseThenLinearise is kept as the backup rather
	// than as an equal. NPC is the canonical Nguyen-Peraire-Cockburn method --
	// Newton on the full ( q, psi, psihat ) system, hybridized elimination for
	// the Jacobian -- which is what the GS papers' discretisation actually is,
	// and it leaves every element-local operation a LINEAR solve against one
	// factorisation. Asserted here because it is a decision rather than a
	// preference: switching it back would put a non-linear solve inside every
	// element of every residual evaluation and change what the whole suite is
	// measuring.
	//
	// A THIRD VALUE USED TO BE HERE. MFEM's NLOrdering::LineariseThenCondense
	// claimed to be this method and was a condensation with the linearisation
	// kept as hidden state; upstream deleted it, and meq's Relinearised -- the
	// wrapper that paired the residual with the gradient to work around exactly
	// that -- went with it.
	BOOST_TEST( static_cast<int>( solver.nonlinearOrdering() )
	                == static_cast<int>( N::NPC ),
	            "the default ordering is not NPC" );

	for ( G choice : { G::PicardOnly, G::AndersonPicard, G::PicardThenNewton, G::None } )
	{
		solver.setGlobalisation( choice );
		BOOST_TEST( static_cast<int>( solver.globalisation() ) == static_cast<int>( choice ),
		            "setGlobalisation did not take for choice "
		            << static_cast<int>( choice ) );
	}

	for ( N choice : { N::NPC, N::CondenseThenLinearise } )
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
 * THE NORMALISED PATH'S REFUSALS, WHICH ARE FIVE AND ARE ALL LOAD BEARING.
 *
 * setSource( NormalisedSource &, double ) turns psi_ax into an unknown and closes
 * the system by a bordered Newton, and each thing it refuses is refused because
 * the alternative is a solve that answers a different question without saying so:
 *
 *   - a normalisation of zero, which is where the degenerate branch lives;
 *   - a second source, for the reason above;
 *   - NOTHING about the ordering. Both surviving orderings are accepted, and
 *     they carry psi_ax differently rather than one of them badly -- see the
 *     block below and solveWithNormalisation();
 *   - any Globalisation but None: KINSOL drives a residual of its own and the
 *     Picard paths build no Jacobian to put a border on. Refused at solve()
 *     rather than at the setter, because the two can be set in either order;
 *   - axisFlux() before prepare(), which would recover from a system that does
 *     not exist.
 */
BOOST_AUTO_TEST_CASE( theNormalisedPathRefusesWhatItCannotDo )
{
	auto profile = std::make_shared<meq::ConstantProfile>( -0.5 );
	mfem::ConstantCoefficient zero( 0.0 );

	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		meq::NormalisedMHDSource source( profile, profile, 1.0, 1.0 );

		BOOST_CHECK_THROW( solver.setSource( source, 0.0 ), std::invalid_argument );
		BOOST_CHECK( !solver.normalisationIsUnknown() );

		solver.setSource( source, 0.5 );
		BOOST_CHECK( solver.normalisationIsUnknown() );
		BOOST_CHECK_EQUAL( solver.psiAxis(), 0.5 );

		// A second source of any kind, through the same one-solver-one-source
		// check the other overloads use.
		BOOST_CHECK_THROW( solver.setSource( source, 0.5 ), std::logic_error );
		BOOST_CHECK_THROW( solver.setSource( zero ), std::logic_error );
	}

	/*
	 * AND IT REFUSES NEITHER ORDERING, WHICH IS A CHANGE TWICE OVER.
	 *
	 * A refusal of MFEM's LineariseThenCondense stood here once, written from
	 * darcyhybridization.hpp's summary rather than from the code under it, and
	 * was removed when both orderings were measured to the same psi_ax. That
	 * mode is now deleted upstream and the value is gone; what is left is NPC,
	 * under which the question does not arise at all.
	 *
	 * Under NPC psi is an unknown of the system rather than a function of the
	 * trace, so the border row is EXACTLY the unit vector -e_j and the corner is
	 * EXACTLY one -- neither is differenced, so there is no residual whose
	 * linearisation history could make a difference mean something other than a
	 * derivative. Only c = dR/ds is still differenced, and it is differenced in
	 * a scalar on both paths.
	 *
	 * So this asserts both orderings are ACCEPTED. A throw here would mean a
	 * refusal came back without the measurement being redone.
	 */
	for ( auto ordering : { meq::GradShafranovSolver::NonlinearOrdering::NPC,
	                        meq::GradShafranovSolver::NonlinearOrdering::CondenseThenLinearise } )
	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		meq::NormalisedMHDSource source( profile, profile, 1.0, 1.0 );

		solver.setNonlinearOrdering( ordering );
		BOOST_CHECK_NO_THROW( solver.setSource( source, 0.5 ) );
		BOOST_CHECK( solver.normalisationIsUnknown() );
	}

	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		meq::NormalisedMHDSource source( profile, profile, 1.0, 1.0 );

		solver.setSource( source, 0.5 );
		solver.setBoundaryData( zero );

		// Before prepare() there is nothing to recover from.
		mfem::Vector trace( solver.traceSpace().GetVSize() );
		trace = 0.0;
		BOOST_CHECK_THROW( solver.axisFlux( trace ), std::logic_error );

		// And the coupling round-trips, since it decides which terms are built.
		BOOST_CHECK( solver.normalisationCoupling()
		             == meq::GradShafranovSolver::Normalisation::Coupled );
		solver.setNormalisationCoupling(
			meq::GradShafranovSolver::Normalisation::Decoupled );
		BOOST_CHECK( solver.normalisationCoupling()
		             == meq::GradShafranovSolver::Normalisation::Decoupled );
		solver.setNormalisationCoupling(
			meq::GradShafranovSolver::Normalisation::Coupled );

		solver.setGlobalisation(
			meq::GradShafranovSolver::Globalisation::PicardThenNewton );
		BOOST_CHECK_THROW( solver.solve(), std::logic_error );
	}
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

/*
 * Threading the element-local assembly must not change the answer.
 *
 * MFEM's DarcyHybridization::SetAssemblyMode() parallelises the element half of
 * ComputeH() and documents that the two modes agree BIT FOR BIT -- not to a
 * tolerance. That is a strong claim and it rests on something specific: the
 * element-local arithmetic is per element and so reassociates nothing, and the
 * scatter into the trace matrix stays serial and in element order. It is worth
 * meq asserting rather than trusting, for two reasons.
 *
 * It is the property that makes the option usable at all. A threaded assembly
 * that agreed only to round-off would put a run-to-run difference into every
 * convergence table in this suite, and the tables assert rates to two decimal
 * places. Bit-for-bit is what lets setAssemblyMode() be a pure performance knob
 * that no other test has to know about.
 *
 * And it is exactly the kind of property that decays quietly. A future
 * `static` reinstated on some element-local scratch would show up here as a
 * handful of wrong entries on the smallest mesh -- MFEM's own test found
 * precisely that, 143 stored entries against 144 -- and would show up nowhere
 * else in meq until a rate table moved for no visible reason.
 *
 * Equality is asserted on the SPARSITY as well as the values, and separately,
 * because the two failures mean different things: different values is a
 * reduction-order question, a different structure is not.
 */
BOOST_AUTO_TEST_CASE( threadedAssemblyReproducesSerialAssemblyExactly )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	std::printf( "\n  threaded assembly: skipped, this MFEM has no OpenMP or no "
	             "thread safety\n" );
#else
	using AM = meq::GradShafranovSolver::AssemblyMode;

	meq::analytic::SolovievEquilibrium const eq =
		meq::analytic::SolovievEquilibrium::nstx();
	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psi( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	std::printf( "\n  threaded assembly against serial, entry by entry\n" );

	// Two degrees and two meshes: the smallest case is the one MFEM's own test
	// found the scratch hazard on, and the largest is where there is enough work
	// per thread for the scheduler to interleave differently between runs.
	for ( int order : { 1, 2 } )
	{
		for ( int n : { 4, 16 } )
		{
			mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), n );

			// Both solvers are kept alive at once, deliberately: each owns its
			// own reduced operator, so comparing them entry by entry needs both
			// still standing. A helper returning one matrix would leave a dangling
			// pointer the moment its solver went out of scope.
			auto serialSolver = std::make_unique<meq::GradShafranovSolver>( mesh, order );
			serialSolver->setAssemblyMode( AM::Serial );
			serialSolver->setSource( source );
			serialSolver->setBoundaryData( psi );
			serialSolver->prepare();

			auto threadedSolver = std::make_unique<meq::GradShafranovSolver>( mesh, order );
			threadedSolver->setAssemblyMode( AM::Threaded );
			threadedSolver->setSource( source );
			threadedSolver->setBoundaryData( psi );
			threadedSolver->prepare();

			auto *a = dynamic_cast<mfem::SparseMatrix *>( &serialSolver->reducedOperator() );
			auto *b = dynamic_cast<mfem::SparseMatrix *>( &threadedSolver->reducedOperator() );

			BOOST_TEST_REQUIRE( ( a != nullptr && b != nullptr ),
			                    "the reduced operator is not a SparseMatrix, so the "
			                    "two assembly modes cannot be compared entry by entry" );

			BOOST_TEST_REQUIRE( a->Height() == b->Height(),
			                    "k = " << order << ", n = " << n << ": threaded "
			                    "assembly produced a trace system of a different SIZE" );
			BOOST_TEST_REQUIRE( a->NumNonZeroElems() == b->NumNonZeroElems(),
			                    "k = " << order << ", n = " << n << ": threaded "
			                    "assembly stored " << b->NumNonZeroElems()
			                    << " entries against serial's " << a->NumNonZeroElems()
			                    << ". That is the element-local scratch being shared "
			                    "between threads, not a rounding question -- MFEM's "
			                    "own test catches the same fault the same way" );

			int const nnz = a->NumNonZeroElems();
			double worst = 0.0;
			int columnMismatches = 0;
			for ( int i = 0; i < nnz; ++i )
			{
				if ( a->GetJ()[ i ] != b->GetJ()[ i ] )
					++columnMismatches;
				worst = std::max( worst,
				                  std::fabs( a->GetData()[ i ] - b->GetData()[ i ] ) );
			}

			std::printf( "    k = %d, n = %2d : %6d entries, worst difference %.3e%s\n",
			             order, n, nnz, worst,
			             worst == 0.0 ? "  (exact)" : "  *** NOT EXACT ***" );
			std::fflush( stdout );

			BOOST_TEST( columnMismatches == 0,
			            "k = " << order << ", n = " << n << ": " << columnMismatches
			            << " entries sit in different COLUMNS between the two assembly "
			            "modes. The structures differ, which is worse than the values "
			            "differing" );

			// NOT a tolerance. MFEM documents exactness and the mechanism for it is
			// specific; if this ever needs slack, the mechanism has changed and the
			// right response is to find out how, not to widen the gate.
			BOOST_TEST( worst == 0.0,
			            "k = " << order << ", n = " << n << ": threaded and serial "
			            "assembly differ by " << worst << " in the worst entry. MFEM "
			            "documents these as bit for bit, so this is a real change "
			            "rather than a tolerance to widen -- the element-local "
			            "arithmetic has started reassociating, or the scatter is no "
			            "longer in element order" );
		}
	}
#endif
}

/*
 * Every trace solver this build has must reach the same equilibrium.
 *
 * setTraceSolver() is a PERFORMANCE choice and this is the test that entitles it
 * to be one. UMFPack, PARDISO and cuDSS are three unrelated sparse LU
 * implementations with three different orderings and three different pivoting
 * strategies; nothing but measurement says they agree, and if they ever stop
 * agreeing the whole idea of a selectable solver is unsound rather than merely
 * slower.
 *
 * **The gate is 1e-10 and the measured difference is around 1e-14**, so there is
 * four orders of slack. That is deliberate: the two are solving the same matrix
 * to different roundings, not approximating each other, and a gate at 1e-13
 * would be a test of the machine's arithmetic rather than of meq. If this ever
 * fails it will fail by a mile -- a wrong matrix type, a wrong ordering, a
 * solver handed a matrix it misdescribes -- not by a factor of three.
 *
 * cuDSS is NOT exercised here and its absence is not an oversight: it reads its
 * matrix and vectors through MFEM's device-aware accessors, so it needs an
 * mfem::Device configured for CUDA, and mfem::Device is global state that must
 * be set before any other MFEM object exists. Configuring one inside a test case
 * would change how every other test in this binary allocates. cuDSS's agreement
 * is checked instead by tests/performance/TraceSolverScaling, which configures
 * the device at program start and returns non-zero if any solver disagrees --
 * registered as the ctest `cuDSSTraceSolver`.
 */
BOOST_AUTO_TEST_CASE( theTraceSolversAgree )
{
	using TS = meq::GradShafranovSolver::TraceSolver;

	meq::analytic::SolovievEquilibrium const eq =
		meq::analytic::SolovievEquilibrium::nstx();
	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psi( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	struct Named { TS choice; char const *name; };
	std::vector<Named> const all = {
		{ TS::UMFPack, "UMFPack" },
		{ TS::Pardiso, "Pardiso" },
		{ TS::cuDSS,   "cuDSS"   } };

	std::printf( "\n  the trace solvers, against each other\n" );

	int const order = 2;
	int const n = 16;
	mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), n );

	std::unique_ptr<mfem::GridFunction> reference;
	char const *referenceName = nullptr;
	int available = 0;

	for ( Named const &entry : all )
	{
		if ( !meq::GradShafranovSolver::traceSolverAvailable( entry.choice ) )
		{
			std::printf( "    %-8s : not in this build\n", entry.name );
			continue;
		}
		if ( entry.choice == TS::cuDSS )
		{
			// See the note above. Available, deliberately not driven from here.
			std::printf( "    %-8s : available, checked by the cuDSSTraceSolver "
			             "ctest instead (it needs an mfem::Device)\n", entry.name );
			continue;
		}

		++available;
		meq::GradShafranovSolver solver( mesh, order );
		solver.setTraceSolver( entry.choice );
		BOOST_TEST( static_cast<int>( solver.traceSolver() )
		                == static_cast<int>( entry.choice ),
		            "setTraceSolver did not take" );
		solver.setSource( source );
		solver.setBoundaryData( psi );
		solver.solve();

		if ( !reference )
		{
			reference = std::make_unique<mfem::GridFunction>( solver.potential() );
			referenceName = entry.name;
			std::printf( "    %-8s : reference, %d trace dofs\n",
			             entry.name, solver.trace().Size() );
			continue;
		}

		mfem::GridFunction difference( solver.potential() );
		difference -= *reference;
		double const relative =
			difference.Norml2()/std::max( 1.0e-300, reference->Norml2() );

		std::printf( "    %-8s : agrees with %s to %.3e\n",
		             entry.name, referenceName, relative );
		std::fflush( stdout );

		BOOST_TEST( relative < 1.0e-10,
		            entry.name << " and " << referenceName << " reach solutions "
		            "differing by " << relative << " relative. Two direct solvers "
		            "on the same matrix must agree to round-off; check the matrix "
		            "type each is given before anything else, since meq's trace "
		            "matrix is structurally symmetric and, on the extension path, "
		            "not symmetric in its values" );
	}

	BOOST_TEST( available >= 1,
	            "no trace solver at all is available in this build, so the "
	            "convergence suite cannot have run either" );
}

/*
 * THE TWO ORDERINGS REACH THE SAME DISCRETE SOLUTION, AND EXACTLY ONE OF THEM
 * ITERATES INSIDE AN ELEMENT.
 *
 * This is the test that says CondenseThenLinearise is a BACKUP rather than a
 * different problem, and it is the only thing in the suite that can tell the
 * two apart from outside. Two claims, and they are different in kind.
 *
 * THE SOLUTIONS AGREE. NPC and the condensation are different methods -- Newton
 * on the full ( q, psi, psihat ) system against Newton on a residual whose
 * element-local eliminations are themselves non-linear solves -- and they
 * converge to the same discrete solution because they discretise the same
 * equations. A wrong sign, a dropped 1/r or a load left out of one of them
 * would show O( 1 ) here, not round-off. Both are also checked against the
 * exact solution, so "they agree" cannot be satisfied by both being wrong the
 * same way.
 *
 * ONLY ONE ITERATES LOCALLY, and that is what says the NPC path is NPC rather
 * than a condensation wearing the name -- which is not a hypothetical failure:
 * MFEM's NLOrdering::LineariseThenCondense claimed exactly that and was
 * deleted for not being it. localNonlinearIterations() must read EXACTLY zero
 * under NPC, because every element-local operation there is one linear solve
 * against one factorisation, and must read a large number under the
 * condensation, where it is one non-linear solve per element per residual
 * evaluation. An equality against zero rather than a tolerance, because the
 * quantity is a count.
 */
BOOST_AUTO_TEST_CASE( theOrderingsAgreeAndOnlyOneIteratesLocally )
{
	using N = meq::GradShafranovSolver::NonlinearOrdering;

	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	std::printf( "\n  the two non-linear orderings, Example 5, against each other\n" );
	std::printf( "    %2s %3s %-22s %7s %14s %14s %14s\n",
	             "k", "n", "ordering", "Newton", "local NL its",
	             "L2 vs exact", "rel. in psi" );

	for ( int order = 1; order <= 3; ++order )
	{
		int const n = 8;

		std::unique_ptr<mfem::GridFunction> reference;

		for ( N ordering : { N::NPC, N::CondenseThenLinearise } )
		{
			bool const npc = ordering == N::NPC;
			char const *name = npc ? "NPC" : "CondenseThenLinearise";

			mfem::Mesh mesh = makeMesh( n );
			meq::GradShafranovSolver solver( mesh, order );
			solver.setSource( source );
			solver.setBoundaryData( exact );
			solver.setNonlinearOrdering( ordering );
			solver.solve();

			long const local = solver.localNonlinearIterations();
			double const error = solver.potentialError( exact );

			double relative = 0.0;
			bool const isReference = !reference;
			if ( isReference )
			{
				reference = std::make_unique<mfem::GridFunction>( solver.potential() );
			}
			else
			{
				mfem::GridFunction difference( solver.potential() );
				difference -= *reference;
				relative = difference.Norml2()
				           /std::max( 1.0e-300, reference->Norml2() );
			}

			std::printf( "    %2d %3d %-22s %7d %14ld %14.6e %14.6e\n",
			             order, n, name, solver.newtonIterations(), local,
			             error, relative );

			if ( npc )
				BOOST_TEST( local == 0L,
				            "k = " << order << ": NonlinearOrdering::NPC ran "
				            << local << " element-local NON-LINEAR iterations, and "
				            "NPC has none by construction. Either the solve did not "
				            "take the NPC path, or DarcyNPCOperator is reaching a "
				            "condensation somewhere. This is the acceptance signal "
				            "for the ordering and it is an equality, not a "
				            "tolerance" );
			else
				BOOST_TEST( local > 0L,
				            "k = " << order << ": CondenseThenLinearise ran NO "
				            "element-local non-linear iterations, which it cannot "
				            "do on a non-linear source -- eliminating flux and "
				            "potential on an element IS a non-linear solve there. "
				            "Either the ordering did not take or the source stopped "
				            "being non-linear" );

			// Loose because it is a floor, not a rate: Example 5 at n = 8 is a
			// coarse mesh and this is only here so that "the two agree" cannot be
			// satisfied by both of them being wrong.
			BOOST_TEST( error < 1.0e-1,
			            "k = " << order << ", " << name << ": L2 against the exact "
			            "solution is " << error << ", which is not a solution of "
			            "Example 5 at all" );

			// THE AGREEMENT ITSELF. Two iterations with different non-linear
			// structure, converged to their own tolerances, landing on the same
			// discrete solution. The bound is set from the measurement printed
			// above rather than chosen: see the table in CLAUDE.md.
			if ( !isReference )
				BOOST_TEST( relative < 1.0e-9,
				            "k = " << order << ": the two orderings differ by "
				            << relative << " relative in psi, which is far above "
				            "round-off. They discretise the same equations and must "
				            "reach the same discrete solution -- a difference this "
				            "size is one of them solving a different problem, not "
				            "two iterations stopping at different places. Note that "
				            "a COARSE discretisation can carry more than one "
				            "solution and two methods can land on different ones; "
				            "check the L2 column before assuming an assembly fault" );
		}
	}
}

/*
 * A solver this build does not have must be refused, not substituted.
 */
BOOST_AUTO_TEST_CASE( anUnavailableTraceSolverIsRefused )
{
	using TS = meq::GradShafranovSolver::TraceSolver;

	meq::tests::Rectangle const box = meq::tests::standardBox();
	mfem::Mesh mesh = meq::tests::makeMesh( box, 4 );
	meq::GradShafranovSolver solver( mesh, 1 );

	for ( TS choice : { TS::UMFPack, TS::Pardiso, TS::cuDSS } )
	{
		if ( meq::GradShafranovSolver::traceSolverAvailable( choice ) )
		{
			BOOST_CHECK_NO_THROW( solver.setTraceSolver( choice ) );
		}
		else
		{
			// invalid_argument and not logic_error: it is the ARGUMENT that this
			// build cannot honour, which is the same distinction the constructor
			// draws for a bad degree.
			BOOST_CHECK_THROW( solver.setTraceSolver( choice ),
			                   std::invalid_argument );
		}
	}
}
