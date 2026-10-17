#define BOOST_TEST_MODULE SolverContract
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>

#include "mfem.hpp"

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
	// For omp_get_max_threads() alone, so that the nonlinear threaded-assembly
	// case can say whether the run it is reporting could have exhibited a race.
	#include <omp.h>
#endif

#include "RoundOffFloor.hpp"

#include "meq/Config.hpp"
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

	// One dimension is not ( R, z ), and neither is three. A 3D mesh reaching the
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
 * tau's SCALING IS SET BEFORE THE FORMS, AND REFUSED WHERE IT HAS NO MEANING.
 *
 * Constant is the default, so every caller written before the choice existed
 * solves what it solved. Diffusion leaves tau() reading the constructor's
 * number -- it is tau0 now, not a different tau -- and it TAKES: the same
 * problem solved both ways is a different discrete answer. After a solve the
 * setter refuses, like setSource(), rather than rebuilding under the caller.
 * A mesh reaching R < 0 is refused for Diffusion, where 1/R is negative, and
 * accepted for Constant, which never reads R.
 *
 * THE FLOOR RADIUS R_f in tau0/max( R, R_f ) defaults to
 * DiffusionScaledStabilization::defaultFloorRadius -- the same number the
 * configuration's TauFloorRadius defaults to, asserted here because the
 * MFEM-free header cannot name the solver's -- and the two-argument form sets
 * it. It refuses, and leaves the solver as it was: a floor beside Constant,
 * which would ignore it; R_f <= 0 or not finite, R_f = 0 being the unbounded
 * tau0/R measured defective at the axis (MEASUREMENTS.md M-290); and any call
 * once the forms are built.
 */
BOOST_AUTO_TEST_CASE( theTauScalingIsSetBeforeTheFormsAndRefusedWhereItIsMeaningless )
{
	using S = meq::GradShafranovSolver::TauScaling;

	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	EquilibriumSource<meq::analytic::SolovievEquilibrium> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	mfem::Mesh mesh = makeMesh( 4 );
	meq::GradShafranovSolver constant( mesh, 2, 1.5 );
	BOOST_TEST( static_cast<int>( constant.tauScaling() ) == static_cast<int>( S::Constant ),
	            "the default tau scaling is not Constant" );
	constant.setSource( source );
	constant.setBoundaryData( exact );
	constant.solve();
	BOOST_CHECK_THROW( constant.setTauScaling( S::Diffusion ), std::logic_error );
	BOOST_CHECK_THROW( constant.setTauScaling( S::Diffusion, 0.2 ), std::logic_error );

	mfem::Mesh scaledMesh = makeMesh( 4 );
	meq::GradShafranovSolver scaled( scaledMesh, 2, 1.5 );
	scaled.setTauScaling( S::Diffusion );
	BOOST_TEST( static_cast<int>( scaled.tauScaling() ) == static_cast<int>( S::Diffusion ),
	            "setTauScaling( Diffusion ) did not take" );
	BOOST_TEST( scaled.tau() == 1.5, "tau() under Diffusion reported " << scaled.tau()
	            << ", not the tau0 the constructor was given" );
	scaled.setSource( source );
	scaled.setBoundaryData( exact );
	scaled.solve();

	mfem::Vector difference( scaled.potential() );
	difference -= constant.potential();
	BOOST_TEST( difference.Normlinf() > 1.0e-8,
	            "the two scalings gave the same psi_h to " << difference.Normlinf()
	            << ", so the scaled stabilisation was never installed" );
	BOOST_TEST( scaled.potentialError( exact ) < 3.0*constant.potentialError( exact ),
	            "Diffusion's psi error " << scaled.potentialError( exact )
	            << " is not of the size of Constant's " << constant.potentialError( exact ) );

	mfem::Mesh across = meq::tests::makeMesh( meq::tests::Rectangle{ -0.2, 0.6, -0.4, 0.4 }, 4 );
	meq::GradShafranovSolver negative( across, 1 );
	BOOST_CHECK_THROW( negative.setTauScaling( S::Diffusion ), std::invalid_argument );
	BOOST_CHECK_THROW( negative.setTauScaling( S::Diffusion, 0.2 ), std::invalid_argument );
	BOOST_CHECK_NO_THROW( negative.setTauScaling( S::Constant ) );

	// The floor radius.
	double const fallback = meq::DiffusionScaledStabilization::defaultFloorRadius;
	BOOST_TEST( meq::DiscretisationConfig().tauFloorRadius == fallback,
	            "[discretisation] TauFloorRadius defaults to " << meq::DiscretisationConfig().tauFloorRadius
	            << " m and the solver to " << fallback << " m; a run through the driver and one through "
	            "the library would solve different problems" );
	BOOST_TEST( scaled.tauFloorRadius() == fallback );

	mfem::Mesh flooredMesh = makeMesh( 4 );
	meq::GradShafranovSolver floored( flooredMesh, 2, 1.5 );
	BOOST_TEST( floored.tauFloorRadius() == fallback );
	BOOST_CHECK_THROW( floored.setTauScaling( S::Constant, 0.2 ), std::invalid_argument );
	for ( double bad : { 0.0, -0.1, std::numeric_limits<double>::infinity(),
	                     std::numeric_limits<double>::quiet_NaN() } )
	{
		BOOST_CHECK_THROW( floored.setTauScaling( S::Diffusion, bad ), std::invalid_argument );
		BOOST_CHECK_THROW( meq::DiffusionScaledStabilization( 1.5, bad ), std::invalid_argument );
	}
	BOOST_TEST( static_cast<int>( floored.tauScaling() ) == static_cast<int>( S::Constant ),
	            "a refused floor radius changed the scaling" );
	BOOST_TEST( floored.tauFloorRadius() == fallback, "a refused floor radius was kept" );
	floored.setTauScaling( S::Diffusion, 0.25 );
	BOOST_TEST( static_cast<int>( floored.tauScaling() ) == static_cast<int>( S::Diffusion ) );
	BOOST_TEST( floored.tauFloorRadius() == 0.25 );
	BOOST_TEST( floored.tau() == 1.5, "setting the floor moved tau0 to " << floored.tau() );
}

/*
 * THE ASSEMBLED FACE TERM IS < ( tau0/max( R, R_f ) ) psi, w > AT EVERY
 * QUADRATURE POINT.
 *
 * Eval() itself first: at every face quadrature point of a mesh REACHING THE
 * AXIS, from both sides, it returns tau0/max( R, R_f ) EXACTLY, R read off
 * the side's own transformation -- with R_f = 0.3 on a 0.2 m mesh, so the
 * floor binds on the first two columns of elements and not beyond, both
 * counts asserted. Then MFEM's own HDGDiffusionIntegrator, handed
 * meq::DiffusionScaledStabilization, assembles every face, and each block is
 * compared with the same rule summed here at the face's physical point --
 * which is what would catch Eval() reading the reference point, the wrong
 * side's point or a stale one. Faces lying on R = 0 take tau0/R_f with no
 * rule of their own, and are checked on faces that exist: the count of axis
 * faces is asserted, not assumed.
 *
 * And it is assembled as a BILINEAR form, which MFEM refuses for a
 * stabilisation reporting IsConstant() false -- "a state dependent
 * stabilization makes the face term nonlinear" -- so getting here at all is
 * the statement that position dependence is not state dependence and no
 * EvalGrad() is owed. IsConstant() is asserted directly as well.
 */
BOOST_AUTO_TEST_CASE( theDiffusionScaledFaceTermIsTau0OverTheFlooredRadiusAtEveryQuadraturePoint )
{
	double const tau0 = 2.5;
	double const floorRadius = 0.3;
	int const order = 2;
	mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::Rectangle{ 0.0, 0.8, -0.4, 0.4 }, 4 );
	mfem::L2_FECollection potentialColl( order, 2, mfem::BasisType::GaussLobatto );
	mfem::DG_Interface_FECollection traceColl( order, 2 );
	mfem::FiniteElementSpace potentialSpace( &mesh, &potentialColl );
	mfem::FiniteElementSpace traceSpace( &mesh, &traceColl );

	meq::DiffusionScaledStabilization const stab( tau0, floorRadius );
	BOOST_TEST( stab.IsConstant(), "a position-dependent tau reported itself state dependent" );
	BOOST_TEST( stab.tau0() == tau0 );
	BOOST_TEST( stab.floorRadius() == floorRadius );
	BOOST_TEST( stab.at( 0.5 ) == tau0/0.5 );
	BOOST_TEST( stab.at( 0.1 ) == tau0/floorRadius );
	BOOST_TEST( stab.at( 0.0 ) == tau0/floorRadius );
	BOOST_TEST( meq::DiffusionScaledStabilization( tau0 ).floorRadius()
	            == meq::DiffusionScaledStabilization::defaultFloorRadius );

	mfem::ConstantCoefficient one( 1.0 );
	mfem::HDGDiffusionIntegrator integrator( one, 1.0 );
	integrator.SetStabilization( stab );
	mfem::IntegrationRule const &rule = mfem::IntRules.Get( mfem::Geometry::SEGMENT, 2*order + 2 );
	integrator.SetIntRule( &rule );

	// Eval(), sampled: every face point, from each side's own transformation.
	int flooredPoints = 0, scaledPoints = 0;
	for ( int f = 0; f < mesh.GetNumFaces(); ++f )
	{
		mfem::FaceElementTransformations *ft = mesh.GetFaceElementTransformations( f );
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			ft->SetAllIntPoints( &rule.IntPoint( q ) );
			for ( int side = 0; side < ( ft->Elem2No >= 0 ? 2 : 1 ); ++side )
			{
				mfem::ElementTransformation &tr = side == 0 ? *ft->Elem1 : *ft->Elem2;
				mfem::Vector x;
				tr.Transform( tr.GetIntPoint(), x );
				double const wanted = tau0/std::max( static_cast<double>( x( 0 ) ), floorRadius );
				double const got = stab.Eval( 0.0, 0.0, 0.0, 0.0, tr );
				BOOST_TEST( got == static_cast<mfem::real_t>( wanted ),
				            "face " << f << ", point " << q << " at R = " << x( 0 ) << ": tau is "
				            << got << ", not tau0/max( R, R_f ) = " << wanted );
				if ( x( 0 ) < floorRadius )
					++flooredPoints;
				else
					++scaledPoints;
			}
		}
	}
	BOOST_TEST( flooredPoints > 0, "no face point sits below R_f, so the floor was not sampled" );
	BOOST_TEST( scaledPoints > 0, "every face point sits below R_f, so tau0/R was not sampled" );

	int axisFaces = 0, interiorFaces = 0;
	double worst = 0.0;
	for ( int f = 0; f < mesh.GetNumFaces(); ++f )
	{
		mfem::FaceElementTransformations *ft = mesh.GetFaceElementTransformations( f );
		mfem::FiniteElement const &el1 = *potentialSpace.GetFE( ft->Elem1No );
		bool const interior = ft->Elem2No >= 0;
		mfem::FiniteElement const &el2 = interior ? *potentialSpace.GetFE( ft->Elem2No ) : el1;
		mfem::FiniteElement const &traceEl = *traceSpace.GetFaceElement( f );
		mfem::DenseMatrix elmat;
		integrator.AssembleHDGFaceMatrix( traceEl, el1, el2, *ft, elmat );

		int const n1 = el1.GetDof(), n2 = interior ? el2.GetDof() : 0, nt = traceEl.GetDof();
		mfem::DenseMatrix expected( n1 + n2 + nt );
		expected = 0.0;
		mfem::Vector shape1( n1 ), shape2( n2 ), traceShape( nt ), x;
		bool onAxis = true;
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			ft->SetAllIntPoints( &ip );
			ft->Transform( ip, x );
			onAxis = onAxis && x( 0 ) == 0.0;
			double const faceWeight = ip.weight*ft->Weight();
			el1.CalcPhysShape( *ft->Elem1, shape1 );
			traceEl.CalcShape( ip, traceShape );
			double const tau1 = tau0/std::max( static_cast<double>( x( 0 ) ), floorRadius );
			double tauSum = tau1;
			for ( int i = 0; i < n1; ++i )
			{
				for ( int j = 0; j < n1; ++j )
					expected( i, j ) += faceWeight*tau1*shape1( i )*shape1( j );
				for ( int j = 0; j < nt; ++j )
					expected( i, n1 + n2 + j ) -= faceWeight*tau1*shape1( i )*traceShape( j );
			}
			if ( interior )
			{
				el2.CalcPhysShape( *ft->Elem2, shape2 );
				double const tau2 = tau0/std::max( static_cast<double>( x( 0 ) ), floorRadius );
				tauSum += tau2;
				for ( int i = 0; i < n2; ++i )
				{
					for ( int j = 0; j < n2; ++j )
						expected( n1 + i, n1 + j ) += faceWeight*tau2*shape2( i )*shape2( j );
					for ( int j = 0; j < nt; ++j )
						expected( n1 + i, n1 + n2 + j ) -= faceWeight*tau2*shape2( i )*traceShape( j );
				}
			}
			for ( int i = 0; i < nt; ++i )
				for ( int j = 0; j < nt; ++j )
					expected( n1 + n2 + i, n1 + n2 + j ) -= faceWeight*tauSum*traceShape( i )*traceShape( j );
		}
		axisFaces += onAxis;
		interiorFaces += interior;

		// The potential-potential and potential-trace blocks and the trace
		// block; the trace-potential block is MFEM's negated transpose and is
		// not a separate claim.
		double scale = 0.0, error = 0.0;
		for ( int i = 0; i < n1 + n2 + nt; ++i )
			for ( int j = 0; j < n1 + n2 + nt; ++j )
			{
				if ( i >= n1 + n2 && j < n1 + n2 )
					continue;
				scale = std::max( scale, static_cast<double>( std::abs( expected( i, j ) ) ) );
				error = std::max( error, static_cast<double>( std::abs( elmat( i, j ) - expected( i, j ) ) ) );
			}
		// The same sum in another order, so in single the format's floor is
		// the whole difference; a tau read off the wrong point is orders above it.
		BOOST_TEST( ( error <= 1.0e-13*scale || !meq::tests::aboveTheFormatsFloor( error, scale ) ),
		            "face " << f << ( onAxis ? " (on the axis)" : "" ) << ": the assembled face term is "
		            << error << " from tau0/max( R, R_f ) summed at its points, against " << scale );
		worst = std::max( worst, error/scale );
	}
	BOOST_TEST( axisFaces == 4, axisFaces << " faces on R = 0, not the 4 a 4 x 4 grid has there" );
	BOOST_TEST( interiorFaces > 0 );
	BOOST_TEST_MESSAGE( "worst relative face-term difference " << worst << " over "
	                    << mesh.GetNumFaces() << " faces, " << axisFaces << " on the axis" );
}

/*
 * THE CONFIGURATION ROUND-TRIPS, and the ones with a valid range enforce it.
 */
BOOST_AUTO_TEST_CASE( theConfigurationRoundTrips )
{
	using G = meq::GradShafranovSolver::Globalisation;
	using L = meq::GradShafranovSolver::LocalSolver;

	mfem::Mesh mesh = makeMesh( 4 );
	meq::GradShafranovSolver solver( mesh, 1 );

	// Compared as integers because Boost.Test insists on printing whatever it
	// compares and a scoped enum has no operator<<.
	BOOST_TEST( static_cast<int>( solver.globalisation() ) == static_cast<int>( G::None ),
	            "the default globalisation is not None" );
	for ( G choice : { G::PicardOnly, G::AndersonPicard, G::PicardThenNewton, G::None } )
	{
		solver.setGlobalisation( choice );
		BOOST_TEST( static_cast<int>( solver.globalisation() ) == static_cast<int>( choice ),
		            "setGlobalisation did not take for choice "
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
	 * A NORMALISED SOURCE IS ACCEPTED, AND THERE IS NO LONGER AN ORDERING AXIS
	 * TO SWEEP IT OVER.
	 *
	 * This looped over both NonlinearOrdering values and asserted that neither
	 * was refused. NPC is now the only ordering MEQ has -- see the note on the
	 * header where the enum was -- so what is left is the acceptance itself.
	 *
	 * Under NPC psi is an unknown of the system rather than a function of the
	 * trace, so the border row is EXACTLY the unit vector -e_j and the corner is
	 * EXACTLY one -- neither is differenced, so there is no residual whose
	 * linearisation history could make a difference mean something other than a
	 * derivative. Only c = dR/ds is still differenced, and it is a scalar.
	 */
	{
		mfem::Mesh mesh = makeMesh( 4 );
		meq::GradShafranovSolver solver( mesh, 1 );
		meq::NormalisedMHDSource source( profile, profile, 1.0, 1.0 );

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

			// FLOORED AT ROUND-OFF ON THE DATUM'S OWN SIZE. In double this is
			// the 1.0e-12 it always was; in single a dof sits 0.11 of one eps
			// of a 0.37 datum away, which no scale-free bound can express.
			// M-219.
			BOOST_TEST( worst < meq::tests::scalarAgreementFloor( 1.0e-12, datumValue ),
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
		/*
		 * THE MFEM-FREE SEAM, WRITTEN OUT. `meq::analytic`'s `flux()` takes
		 * `double &` out-parameters because CLAUDE.md's *Layout* keeps the
		 * Source and Profile hierarchies free of MFEM -- so they stay `double`
		 * whatever `MFEM_PRECISION` says -- while `mfem::Vector::operator()`
		 * yields `real_t &`. In the double build those bind; in single they do
		 * not, and no cast can fix a REFERENCE. MEASUREMENTS.md M-217.
		 */
		double fluxR = 0.0, fluxZ = 0.0;
		eq.flux( x( 0 ), x( 1 ), fluxR, fluxZ );
		value( 0 ) = fluxR;
		value( 1 ) = fluxZ;
	} );

	std::printf( "\n  the reconstructed quantities, Solov'ev on the linear path\n" );
	std::printf( "    %2s %14s %14s %14s %14s\n", "k", "L2( q )", "L2( q* )",
	             "|q + total|", "|| q ||" );

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

		std::printf( "    %2d %14.6e %14.6e %14.6e %14.6e\n",
		             order, fluxError, enrichedError, mismatch, fluxNorm );
		std::fflush( stdout );

		BOOST_TEST( std::isfinite( enrichedError ),
		            "k = " << order << ": the post-processed flux is not finite" );
		// MEASURED: 6.92e-4 against 9.58e-4, 8.90e-6 against 1.58e-5, 1.43e-7
		// against 2.22e-7 at k = 1, 2, 3 -- about 1.5x better at every order. So
		// "better than q" is the assertion, not "not much worse".
		/*
		 * AND THE COMPARISON IS ONLY AVAILABLE WHILE THERE IS A DISCRETISATION
		 * ERROR TO COMPARE. "q* is better than q" is a statement about
		 * TRUNCATION, and at k = 3 on this mesh the double build reads
		 * 1.43e-07 against 2.22e-07 -- both below what an L2 norm of a field
		 * of size || q || can carry in single, where the same two read
		 * 8.80e-06 and 7.90e-06 and their ORDER is a coin toss. Skipped
		 * loudly rather than loosened: widening the bound would leave a case
		 * that cannot tell an enriched reconstruction from an absent one and
		 * reports green. The finiteness check above still runs, and so does
		 * every degree whose truncation is still above the floor. M-220.
		 */
		double const floor = meq::tests::roundOffFloor( fluxNorm );
		if ( fluxError < floor && enrichedError < floor )
			BOOST_TEST_MESSAGE(
				"k = " << order << ": SKIPPED, L2( q ) = " << fluxError
				<< " and L2( q* ) = " << enrichedError << " are both under "
				<< floor << ", which is round-off on a flux of norm "
				<< fluxNorm << ". There is no truncation left here to compare." );
		else
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
				worst = std::max(
					worst,
					static_cast<double>( std::fabs( a->GetData()[ i ]
					                                - b->GetData()[ i ] ) ) );
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
 * AND THE SAME PROPERTY ON A NONLINEAR SOURCE, WHICH IS A DIFFERENT LOOP.
 *
 * The case above cannot see the defect this one is for, and the reason is
 * structural rather than a matter of coverage. It builds its source as an
 * mfem::FunctionCoefficient, which takes meq's LINEAR path:
 * usesNonlinearForms() is false, meq::SourceIntegrator is never installed, and
 * DarcyHybridization::MultNL() is never called. So it exercises ComputeH()'s
 * element loop and nothing else -- which was the only threaded loop in the
 * class when it was written.
 *
 * It is not any more. MFEM now threads MultNL() as well -- the residual and the
 * Jacobian assembly, and so NPCResidual() and NPCGradient(), which is to say
 * every NPC step -- walking the elements in colour order under OpenMP. That
 * loop calls the caller's own integrators, and MFEM says plainly that it cannot
 * check them: "Any integrator the caller installs ... sits on this loop and
 * must be thread-safe too. An integrator holding per-point scratch as a plain
 * member will race, silently."
 *
 * meq::SourceIntegrator held exactly that -- an mfem::Vector shape member,
 * resized and refilled per quadrature point -- and it is now guarded the way
 * MFEM guards its own. THIS TEST IS WHAT SAYS SO. Nothing else in the suite
 * can: the source integrator IS the whole semi-linear term, and CLAUDE.md's
 * *A wrong Jacobian is invisible to a convergence table* records that a
 * degraded dFdPsi leaves every error and every rate unchanged to six
 * significant figures while costing Newton its order. So a corrupted Jacobian
 * here would show up as neither a wrong answer nor a wrong rate.
 *
 * THREE PROPERTIES, AND THE THIRD IS THE ONE WITH TEETH.
 *
 *   - the recovered psi, bit for bit. Colouring changes the ORDER in which a
 *     face's two elements accumulate into the trace row, and a + b == b + a
 *     exactly, so MFEM's exactness claim covers this loop too.
 *   - the flux, likewise, since it is recovered from the same local solves.
 *   - THE NEWTON ITERATION COUNT. This is the assertion that reaches the
 *     Jacobian. A racing residual would move the answer and be caught by the
 *     first two; a racing GRADIENT converges to the same discrete solution by a
 *     different path, so it moves the count and nothing else. That is precisely
 *     the failure a rate table cannot see, and it is why the count is asserted
 *     rather than printed.
 */
BOOST_AUTO_TEST_CASE( threadedAssemblyReproducesSerialAssemblyOnANonlinearSource )
{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
	std::printf( "\n  threaded assembly, nonlinear: skipped, this MFEM has no "
	             "OpenMP or no thread safety\n" );
#else
	using AM = meq::GradShafranovSolver::AssemblyMode;

	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient psi( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	std::printf( "\n  threaded assembly against serial, on a nonlinear source\n" );
	std::printf( "    this is the case that reaches MultNL, and so the case that\n"
	             "    reaches meq::SourceIntegrator\n" );

	// omp_get_max_threads() is printed rather than asserted: a one-thread run
	// cannot demonstrate the absence of a race, and saying so is better than a
	// green tick that means nothing. The gate stays exactness either way.
	std::printf( "    omp_get_max_threads() = %d%s\n", omp_get_max_threads(),
	             omp_get_max_threads() > 1 ? ""
	             : "  <- ONE THREAD: this run cannot exhibit a race" );

	// Both stabilisations: meq::DiffusionScaledStabilization evaluates on the
	// transformation each thread is handed, which is a different claim from a
	// constant's holding one double.
	for ( bool const scaled : { false, true } )
	for ( int order : { 1, 2 } )
	{
		for ( int n : { 4, 16 } )
		{
			mfem::Mesh mesh = meq::tests::makeMesh( meq::tests::standardBox(), n );
			double const tau = scaled ? 3.0 : 1.0;
			auto const scaling = scaled ? meq::GradShafranovSolver::TauScaling::Diffusion
			                            : meq::GradShafranovSolver::TauScaling::Constant;
			if ( order == 1 && n == 4 )
				std::printf( "    tau %s\n", scaled ? "= 3/R" : "= 1" );

			auto serialSolver = std::make_unique<meq::GradShafranovSolver>( mesh, order, tau );
			serialSolver->setTauScaling( scaling );
			serialSolver->setAssemblyMode( AM::Serial );
			serialSolver->setSource( source );
			serialSolver->setBoundaryData( psi );
			serialSolver->solve();

			auto threadedSolver = std::make_unique<meq::GradShafranovSolver>( mesh, order, tau );
			threadedSolver->setTauScaling( scaling );
			threadedSolver->setAssemblyMode( AM::Threaded );
			threadedSolver->setSource( source );
			threadedSolver->setBoundaryData( psi );
			threadedSolver->solve();

			mfem::GridFunction const &a = serialSolver->potential();
			mfem::GridFunction const &b = threadedSolver->potential();
			mfem::GridFunction const &qa = serialSolver->flux();
			mfem::GridFunction const &qb = threadedSolver->flux();

			BOOST_TEST_REQUIRE( a.Size() == b.Size(),
			                    "k = " << order << ", n = " << n << ": the two "
			                    "assembly modes produced potentials of different size" );
			BOOST_TEST_REQUIRE( qa.Size() == qb.Size(),
			                    "k = " << order << ", n = " << n << ": the two "
			                    "assembly modes produced fluxes of different size" );

			double worstPsi = 0.0;
			for ( int i = 0; i < a.Size(); ++i )
				worstPsi = std::max(
			worstPsi, static_cast<double>( std::fabs( a( i ) - b( i ) ) ) );

			double worstFlux = 0.0;
			for ( int i = 0; i < qa.Size(); ++i )
				worstFlux = std::max(
			worstFlux, static_cast<double>( std::fabs( qa( i ) - qb( i ) ) ) );

			int const itsSerial = serialSolver->newtonIterations();
			int const itsThreaded = threadedSolver->newtonIterations();

			std::printf( "    k = %d, n = %2d : psi %.3e, q %.3e, Newton %d against %d%s\n",
			             order, n, worstPsi, worstFlux, itsThreaded, itsSerial,
			             ( worstPsi == 0.0 && worstFlux == 0.0
			               && itsSerial == itsThreaded ) ? "  (exact)"
			             : "  *** DIFFERS ***" );
			std::fflush( stdout );

			BOOST_TEST( worstPsi == 0.0,
			            "k = " << order << ", n = " << n << ": the threaded and "
			            "serial assembly modes recovered potentials differing by "
			            << worstPsi << ". On a nonlinear source this is MultNL's "
			            "loop, so the first thing to check is whether an integrator "
			            "meq installs holds per-point scratch as a member -- "
			            "meq::SourceIntegrator's `shape` is guarded on "
			            "MFEM_THREAD_SAFE for exactly this reason" );

			BOOST_TEST( worstFlux == 0.0,
			            "k = " << order << ", n = " << n << ": the recovered FLUXES "
			            "differ by " << worstFlux << ". The flux comes out of the "
			            "same element-local solves as the potential, so this and the "
			            "potential should fail together; one without the other is "
			            "the more interesting finding" );

			// The Jacobian assertion. A racing gradient reaches the same discrete
			// solution -- Newton converges to the root whatever carried it there --
			// so it shows up here and in no error norm and no rate.
			BOOST_TEST( itsSerial == itsThreaded,
			            "k = " << order << ", n = " << n << ": threaded assembly "
			            "took " << itsThreaded << " Newton iterations against "
			            "serial's " << itsSerial << ", having reached the same "
			            "answer. That is the signature of a corrupted JACOBIAN "
			            "rather than a corrupted residual: Newton converges to the "
			            "same root by a longer path, so no error norm and no "
			            "convergence rate can see it. AssembleElementGrad is where "
			            "to look" );
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
			difference.Norml2()
		       /std::max( 1.0e-300,
		                  static_cast<double>( reference->Norml2() ) );

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
 * NOTHING ITERATES INSIDE AN ELEMENT, WHICH IS THE ACCEPTANCE SIGNAL THAT NPC
 * IS NPC.
 *
 * This case used to sweep both NonlinearOrdering values and assert that they
 * reach the same discrete solution while exactly one of them iterates locally.
 * CondenseThenLinearise is gone -- nothing in MEQ reached it, and this case was
 * one of the three places keeping it alive -- so the comparison has no second
 * arm. What survives is the half that was never about the comparison.
 *
 * localNonlinearIterations() must read EXACTLY zero, because every
 * element-local operation under NPC is one linear solve against one
 * factorisation. That is not a hypothetical failure to guard against: MFEM's
 * NLOrdering::LineariseThenCondense claimed to be the NPC method, was not, and
 * was deleted for it. An equality against zero rather than a tolerance, because
 * the quantity is a count.
 *
 * The L2 against the exact solution is checked beside it so that "no local
 * iterations" cannot be satisfied by a solve that did not solve anything.
 */
BOOST_AUTO_TEST_CASE( theOrderingNeverIteratesInsideAnElement )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	std::printf( "\n  NPC's element-local work, Example 5\n" );
	std::printf( "    %2s %3s %7s %14s %14s\n",
	             "k", "n", "Newton", "local NL its", "L2 vs exact" );

	for ( int order = 1; order <= 3; ++order )
	{
		int const n = 8;

		mfem::Mesh mesh = makeMesh( n );
		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( exact );
		solver.solve();

		long const local = solver.localNonlinearIterations();
		double const error = solver.potentialError( exact );

		std::printf( "    %2d %3d %7d %14ld %14.6e\n",
		             order, n, solver.newtonIterations(), local, error );

		BOOST_TEST( local == 0L,
		            "k = " << order << ": the solve ran " << local
		            << " element-local NON-LINEAR iterations, and NPC has none by "
		            "construction. Either the solve did not take the NPC path, or "
		            "DarcyNPCOperator is reaching a condensation somewhere. This "
		            "is the acceptance signal for the ordering and it is an "
		            "equality, not a tolerance" );

		// Loose because it is a floor, not a rate: Example 5 at n = 8 is a
		// coarse mesh and this is only here so that "no local iterations"
		// cannot be satisfied by a solve that did not solve anything.
		BOOST_TEST( error < 1.0e-1,
		            "k = " << order << ": L2 against the exact solution is "
		            << error << ", which is not a solution of Example 5 at all" );
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

/*
 * AxisConstraint::PrescribedPoint's REFUSALS, which are what keep a prescribed
 * point from being DEFAULTED. MIRROR-PLAN.md MR-2; the acceptance, which needs
 * a solve, is PrescribedAxisPoint.cpp.
 *
 * Every candidate default chooses an equilibrium -- the symmetry axis carries
 * psi = 0, the domain's centre is nobody's plasma, and the last located axis is
 * another constraint's answer -- so choosing the mode without a point throws,
 * and a point that cannot be one throws at the setter rather than at the solve.
 * A refused setter leaves the solver as it was.
 */
BOOST_AUTO_TEST_CASE( aPrescribedAxisPointIsRefusedRatherThanDefaulted )
{
	using A = meq::GradShafranovSolver::AxisConstraint;

	mfem::Mesh mesh = makeMesh( 4 );
	meq::GradShafranovSolver solver( mesh, 1 );

	BOOST_CHECK( solver.axisConstraint() == A::LocatedAxis );
	BOOST_CHECK_THROW( solver.setAxisConstraint( A::PrescribedPoint ),
	                   std::logic_error );
	BOOST_CHECK( solver.axisConstraint() == A::LocatedAxis );

	double const nan = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_THROW( solver.setAxisFluxPoint( nan, 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( solver.setAxisFluxPoint( 1.0, nan ), std::invalid_argument );

	// ON A MESH THAT REACHES THE AXIS, so that R = 0 is a point of Omega and it
	// is the radius check refusing it rather than the mesh. On box() R = 0 is
	// outside the mesh as well, so with the radius check deleted the setter
	// still throws -- a runtime_error from the location, which this case would
	// read as an unexpected exception rather than as the missing check.
	{
		mfem::Mesh axisMesh =
			meq::tests::makeMesh( meq::tests::Rectangle{ 0.0, 1.0, -0.5, 0.5 }, 4 );
		meq::GradShafranovSolver onAxis( axisMesh, 1 );
		BOOST_CHECK_THROW( onAxis.setAxisFluxPoint( 0.0, 0.0 ),
		                   std::invalid_argument );
		BOOST_CHECK_THROW( onAxis.setAxisFluxPoint( 0.0, 0.25 ),
		                   std::invalid_argument );
		BOOST_CHECK( onAxis.axisConstraint() == A::LocatedAxis );
		BOOST_CHECK_NO_THROW( onAxis.setAxisFluxPoint( 0.5, 0.0 ) );
	}
	BOOST_CHECK_THROW( solver.setAxisFluxPoint( -1.0, 0.0 ), std::invalid_argument );

	// Outside box() [ 0.6, 1.4 ] x [ -0.6, 0.6 ], at a positive radius, so it
	// is the mesh refusing it and not the radius check.
	BOOST_CHECK_THROW( solver.setAxisFluxPoint( 2.0, 0.0 ), std::runtime_error );
	BOOST_CHECK_THROW( solver.setAxisFluxPoint( 1.0, 0.9 ), std::runtime_error );
	BOOST_CHECK( solver.axisConstraint() == A::LocatedAxis );
	// And still refused afterwards: a failed setter named no point.
	BOOST_CHECK_THROW( solver.setAxisConstraint( A::PrescribedPoint ),
	                   std::logic_error );

	// A point inside SELECTS the mode, and the point survives a switch away.
	BOOST_CHECK_NO_THROW( solver.setAxisFluxPoint( 1.0, 0.0 ) );
	BOOST_CHECK( solver.axisConstraint() == A::PrescribedPoint );
	solver.setAxisConstraint( A::NodalMaximum );
	BOOST_CHECK( solver.axisConstraint() == A::NodalMaximum );
	BOOST_CHECK_NO_THROW( solver.setAxisConstraint( A::PrescribedPoint ) );
	BOOST_CHECK( solver.axisConstraint() == A::PrescribedPoint );
}
