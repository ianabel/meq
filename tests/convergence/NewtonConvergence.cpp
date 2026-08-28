#define BOOST_TEST_MODULE NewtonConvergence
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/ManufacturedNonlinear.hpp"
#include "analytic/SimilarityExponential.hpp"
#include "analytic/McCarthy.hpp"
#include "analytic/Soloviev.hpp"
#include "convergence/ConvergenceHarness.hpp"

/*
 * The stage-4 acceptance test: the semi-linear HDG Grad-Shafranov solver,
 * closed by Newton, measured against the manufactured solution of
 * refs/HDG-GradShafranov.pdf Example 5 and against the two exact equilibria
 * either side of it on the ladder.
 *
 * Two things have to be true, and they fail independently:
 *
 *   1. the discretisation converges at k+1 in psi and in q, which says the
 *      right equation is being solved;
 *   2. Newton converges quadratically, which says the Jacobian is the
 *      derivative of the residual the solver actually assembles.
 *
 * The second is the one this stage exists to be able to detect. A Jacobian that
 * disagrees with its residual still converges -- to the same answer, at the same
 * spatial rate -- it just grinds down linearly getting there. Nothing in a
 * convergence table sees that, which is why the residual history is asserted on
 * and not merely printed. The shape to expect is CEDRES++ Table 2 on a 577k
 * unknown production problem: relative residual 2.7e0 -> 9.2e-2 -> 1.8e-3 ->
 * 5.3e-6 -> 3.9e-12 in five iterations, an observed order that climbs towards 2
 * as the iterate enters the quadratic regime.
 *
 * Before either of those, two checks that need no solver at all: that each
 * benchmark's source really is the one that makes its psi solve the equation
 * (finite differences of Delta*), and that Example 5's dF/dpsi really is the
 * derivative of its F. Both are cheap and both catch a transcription error that
 * would otherwise present as a solver bug.
 *
 * Then the check that separates the two failure modes above: the assembled
 * Jacobian against a finite difference of the assembled residual, on a small
 * mesh. tests/unit/SourceTests.cpp does the analogous thing to dF/dpsi alone;
 * this does it to the whole reduced operator, so it also covers the hybridized
 * elimination, the local non-linear solves and the essential trace condition.
 *
 * ALL THREE RUNGS OF THE LADDER ARE HERE, and the middle one does most of the
 * work. CLAUDE.md sets out why tests/analytic carries three fixtures rather than
 * one; this is the test that spends them.
 *
 *   Soloviev.hpp     F constant in psi     dF/dpsi = 0    the same problem the
 *                                                         linear path solves, so
 *                                                         the two setSource()
 *                                                         overloads must agree
 *                                                         to round-off
 *   McCarthy.hpp     F affine in psi       dF/dpsi = T    exact, so the whole
 *                                                         discrete residual is
 *                                                         affine and ONE Newton
 *                                                         step must finish it
 *   Manufactured-    psi^2 and e^-psi      varies         the quadratic
 *   Nonlinear.hpp                                         convergence itself
 *
 * The first rung cannot test dF/dpsi at all, since it is zero; the third can, but
 * its algebra is messy enough to hide a factor. The second is the sharp one: T
 * is a single constant, so the Jacobian's mass term is either there or it is
 * not, and the difference is one iteration against eighty.
 *
 * The domain is the same rectangle for all three, which also means the McCarthy
 * and Solov'ev cases are being solved on a box that is not their own plasma
 * boundary. That is deliberate and is what makes the Dirichlet data
 * non-homogeneous.
 *
 * The domain is the rectangle of examples/manufactured.toml, not the paper's
 * ITER-like double-null boundary: the curved boundary is stage 5, and a fitted
 * polygon keeps Gamma_h == Gamma. r is bounded away from zero because the
 * operator and the source both carry a 1/r. The manufactured psi does not
 * vanish on that rectangle, so the Dirichlet data is non-homogeneous, which on
 * the non-linear path exercises a route through DarcyHybridization that had no
 * MFEM regression covering it (see CLAUDE.md).
 *
 * WHAT GOES WRONG WHEN A CONVENTION IS WRONG, measured by introducing each
 * defect and reading the output, as SolovievConvergence.cpp did for the linear
 * ones. At k = 2 unless stated:
 *
 *   SourceIntegrator's sign flipped, in     psi and q flat at 4.7e-1 and 2.8e0,
 *   both the residual and the Jacobian      rate 0.00. Newton still converges
 *                                           quadratically and the Jacobian check
 *                                           still passes -- it is the residual
 *                                           that is wrong, not its derivative
 *
 *   ManufacturedNonlinear::f() returning    psi and q flat at 1.0e-1 and 6.9e-1,
 *   F/r rather than F, which is what the    rate 0.00, and Newton needs 10 to 12
 *   pre-port version of that file did       iterations. Caught first, and much
 *                                           more legibly, by
 *                                           manufacturedSourceMatchesTheOperator
 *
 *   the HDG face stabilisation left on      MFEM aborts in
 *   the linear potential mass form while    SetPotMassNonlinearIntegrator with
 *   the source sits on the non-linear one   "Non-linear mass cannot work with a
 *                                           linear constraint"
 *
 *   dF/dpsi too large by 5%                 every error in the table unchanged
 *                                           to six digits, every rate unchanged,
 *                                           and Newton's observed order exactly
 *                                           1.000 over eight consecutive steps,
 *                                           taking 10 iterations instead of 4.
 *                                           The Jacobian check reports 1.7e-3
 *
 *   dF/dpsi too large by 50%                Newton diverges to NaN and the run
 *                                           aborts
 *
 *   the Jacobian's mass term dropped        on McCarthy, 87/84/80/76 iterations
 *   entirely, that is SourceIntegrator::    down the four meshes instead of one,
 *   AssembleElementGrad returning zero      the residual settling to a constant
 *                                           contraction of 0.7711 a step, and
 *                                           every error unchanged to six digits
 *
 * The last two rows are what this whole stage is arranged around: the answer is
 * right, the table is right to six digits, and the only things that notice are
 * the iteration counts, the residual history and the finite-difference Jacobian
 * check.
 */

/*
 * Everything this file used to define for itself -- the rectangle mesh, the
 * meq::Source adapter, one measured solve, the rate arithmetic, the tables and
 * checkOrder() -- now lives in convergence/ConvergenceHarness.hpp, so that the
 * five benchmarks added since can use the same machinery rather than a second
 * copy of it. What is left here is the box this study is posed on and thin
 * wrappers that supply it.
 */
namespace
{
	using meq::tests::EquilibriumSource;
	using meq::tests::Measurement;
	using meq::tests::checkOrder;
	using meq::tests::newtonOrder;
	using meq::tests::rate;

	// The box of examples/manufactured.toml, which encloses the ITER-like
	// double-null boundary the paper poses Example 5 on.
	//
	// The domain is that rectangle, not the paper's ITER-like double-null
	// boundary: the curved boundary is stage 5, and a fitted polygon keeps
	// Gamma_h == Gamma. r is bounded away from zero because the operator and the
	// source both carry a 1/r. The manufactured psi does not vanish on that
	// rectangle, so the Dirichlet data is non-homogeneous, which on the
	// non-linear path exercises a route through DarcyHybridization that had no
	// MFEM regression covering it (see CLAUDE.md).
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

	template<typename Equilibrium>
	Measurement measure( Equilibrium const &eq, int order, int n,
	                     std::vector<double> *residualHistory = nullptr )
	{
		return meq::tests::measure( eq, box(), order, n, residualHistory );
	}

	void printTable( char const *label, int order,
	                 std::vector<Measurement> const &points )
	{
		meq::tests::printTable( label, order, box(), points );
	}

	/// Four dyadic refinements from 4 cells a side, so three measured rates per
	/// quantity per order -- the same ladder SolovievConvergence.cpp uses.
	std::vector<int> const meshSizes = meq::tests::dyadicMeshes();

	/// k+1, less the slack allowed for a two-mesh rate estimate.
	double const rateSlack = meq::tests::rateSlack;
}

/// The benchmark before the solver, part one: -Delta*( psi ) must equal
/// F( r, z, psi ) evaluated at the manufactured psi, or everything measured
/// against it is measured against the wrong thing. This is what catches the
/// factor of r that the pre-port version of ManufacturedNonlinear.hpp carried
/// -- its operator() returned F/r -- and which no convergence rate can see.
BOOST_AUTO_TEST_CASE( manufacturedSourceMatchesTheOperator )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();

	for ( double r = rMin; r <= rMax + 1.0e-12; r += 0.1 )
	{
		for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.15 )
		{
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, eq.psi( r, z ) );
			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		}
	}
}

/// The benchmark before the solver, part two: dF/dpsi must be the derivative of
/// F. The Newton iteration is the only thing that reads it, and a wrong
/// derivative does not move the converged answer, so this is checked here
/// rather than inferred from the answer.
BOOST_AUTO_TEST_CASE( manufacturedDerivativeMatchesAFiniteDifference )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();

	double const h = 1.0e-6;

	for ( double r = rMin; r <= rMax + 1.0e-12; r += 0.2 )
	{
		for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
		{
			for ( double psi : { -1.0, -0.3, 0.0, 0.4, 1.0 } )
			{
				double const difference = ( eq.f( r, z, psi + h ) - eq.f( r, z, psi - h ) )
				                          /( 2.0*h );
				double const analytic = eq.dFdPsi( r, z, psi );
				BOOST_TEST( std::abs( difference - analytic )
				            < 1.0e-6*( 1.0 + std::abs( analytic ) ),
				            "at ( " << r << ", " << z << ", psi = " << psi
				            << " ): dF/dpsi = " << analytic
				            << " but a difference of F gives " << difference );
			}
		}
	}
}

/// The single highest-value check in the stage: the Jacobian the Newton solve
/// uses, against a finite difference of the residual it is meant to
/// differentiate.
///
/// This is what separates "my Jacobian is wrong" from "my discretisation is
/// wrong", and those two look identical from a convergence table -- both give
/// k+1 in psi. It is done on the assembled reduced operator rather than on
/// dF/dpsi alone, so it covers the hybridized elimination, the element-local
/// non-linear solves and the essential trace condition as well as the source.
///
/// The perturbation is zero on the essential trace dofs deliberately.
/// DarcyHybridization carries the Dirichlet condition the way a NonlinearForm
/// does: the values ride in the iterate, the residual is masked to zero on
/// those rows and the Jacobian is given a unit row. A difference quotient of a
/// masked row is zero while the Jacobian row is the identity, so probing those
/// directions would report a discrepancy that is not one.
BOOST_AUTO_TEST_CASE( jacobianMatchesAFiniteDifferenceOfTheResidual )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> source( eq );

	mfem::Mesh mesh = makeMesh( 3 );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( psiCoeff );
	solver.prepare();

	mfem::Operator &residual = solver.reducedOperator();
	int const size = solver.reducedSolution().Size();

	std::vector<bool> isEssential( size, false );
	mfem::Array<int> const &essential = solver.essentialTraceDofs();
	for ( int i = 0; i < essential.Size(); ++i )
		isEssential[ essential[ i ] ] = true;

	BOOST_TEST( essential.Size() > 0,
	            "no essential trace dofs: the Dirichlet condition is not being imposed, "
	            "which would make this check vacuous and the solver wrong" );

	// Evaluate away from the solution and away from zero, so that the quadratic
	// and exponential terms in F are genuinely being exercised. The wiggle is a
	// deterministic function of the dof index rather than a random draw: a test
	// that fails only sometimes is worse than no test.
	mfem::Vector state( solver.reducedSolution() );
	for ( int i = 0; i < size; ++i )
		if ( !isEssential[ i ] )
			state( i ) += 0.3*std::sin( 1.7*i + 0.4 );

	// Every residual evaluation happens before the gradient is asked for.
	// GetGradient() writes through the same element-block storage the residual
	// path reads, and while the configuration here does not read what it writes,
	// ordering the calls this way means the check does not depend on that.
	double const step = 1.0e-5;
	int const numDirections = 4;

	std::vector<mfem::Vector> directions( numDirections );
	std::vector<mfem::Vector> differences( numDirections );

	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector &direction = directions[ trial ];
		direction.SetSize( size );
		for ( int i = 0; i < size; ++i )
			direction( i ) = isEssential[ i ]
			                 ? 0.0
			                 : std::cos( 2.3*i + 1.1*trial ) + 0.5*std::sin( 0.7*i - 0.3*trial );
		direction /= direction.Norml2();

		mfem::Vector forward( size );
		mfem::Vector backward( size );
		mfem::Vector shifted( state );

		shifted.Add( step, direction );
		residual.Mult( shifted, forward );

		shifted = state;
		shifted.Add( -step, direction );
		residual.Mult( shifted, backward );

		// A central difference: its truncation error is O(step^2), which keeps
		// the comparison well clear of the discrepancy a wrong Jacobian gives.
		mfem::Vector &difference = differences[ trial ];
		difference.SetSize( size );
		difference = forward;
		difference -= backward;
		difference /= 2.0*step;
	}

	mfem::Operator &jacobian = residual.GetGradient( state );

	double worst = 0.0;
	for ( int trial = 0; trial < numDirections; ++trial )
	{
		mfem::Vector applied( size );
		jacobian.Mult( directions[ trial ], applied );

		mfem::Vector error( differences[ trial ] );
		error -= applied;

		double const relative = error.Norml2()/applied.Norml2();
		std::printf( "  Jacobian vs central difference, direction %d: "
		             "||J d|| = %.6e, relative error %.3e\n",
		             trial, applied.Norml2(), relative );
		worst = std::max( worst, relative );
	}
	std::fflush( stdout );

	BOOST_TEST( worst < 1.0e-8,
	            "the assembled Jacobian differs from a central difference of the "
	            "assembled residual by " << worst << " relative; Newton will still "
	            "converge, but linearly" );
}

/// The two setSource() overloads must agree wherever both apply, and the
/// Solov'ev source is where they do: F is constant in psi, so it can be handed
/// over either as a coefficient on the right hand side or as a meq::Source on
/// the non-linear potential mass form, and the two must assemble the same
/// numbers.
///
/// That makes this the tightest available check of SourceIntegrator's sign and
/// its 1/r, because the linear path is not merely self-consistent: it has been
/// measured against a closed-form equilibrium to k+1 in both psi and q, in
/// tests/convergence/SolovievConvergence.cpp. A sign error in the integrator
/// would still let Newton converge quadratically -- measured, it does -- and
/// would still leave the Jacobian consistent with the residual. It just answers
/// a different question. This is what notices.
///
/// It also pins the other end of the Newton path: with dF/dpsi identically zero
/// the Jacobian is the linear operator and one step must finish the job, which
/// is exactly why the Solov'ev benchmark cannot be used to test dF/dpsi and
/// Example 5 has to be.
///
/// Agreement to round-off needs the two paths to integrate F on the same
/// quadrature rule, which they do here and which is why GradShafranov.cpp asks
/// DomainLFIntegrator for order 2k+4 rather than its default 2k. With the
/// default the two differ by 4.3e-7 relative -- not a bug, just two rules on a
/// rational integrand, and well below the 1.2e-5 discretisation error -- but a
/// tolerance that loose would no longer notice a small error in the integrator.
/// The rules coincide because these are affine triangles; on a curved mesh
/// SourceIntegrator's extra ElementTransformation::OrderW() would separate them
/// again, and this tolerance would have to be revisited with it.
BOOST_AUTO_TEST_CASE( theNewtonPathReproducesTheLinearPathOnSoloviev )
{
	meq::analytic::SolovievEquilibrium const eq
		= meq::analytic::SolovievEquilibrium::nstx();
	meq::SolovievSource const source( eq.getA() );

	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );
	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );

	int const order = 2;
	int const cells = 8;

	mfem::Mesh linearMesh = makeMesh( cells );
	meq::GradShafranovSolver linearSolver( linearMesh, order );
	linearSolver.setSource( sourceCoeff );
	linearSolver.setBoundaryData( psiCoeff );
	linearSolver.solve();

	mfem::Mesh newtonMesh = makeMesh( cells );
	meq::GradShafranovSolver newtonSolver( newtonMesh, order );
	newtonSolver.setSource( source );
	newtonSolver.setBoundaryData( psiCoeff );
	newtonSolver.solve();

	BOOST_TEST( linearSolver.isNonlinear() == false,
	            "the coefficient overload should not have selected the Newton path" );
	BOOST_TEST( newtonSolver.isNonlinear() == true,
	            "the meq::Source overload should have selected the Newton path" );

	mfem::Vector difference( newtonSolver.potential() );
	difference -= linearSolver.potential();
	double const relative = difference.Norml2()/linearSolver.potential().Norml2();

	std::printf( "\n  Solov'ev through both paths, k = %d, %d cells a side:\n"
	             "    L2(psi) linear %.6e, Newton %.6e, coefficients differ by %.3e "
	             "relative, %d Newton iteration(s)\n",
	             order, cells,
	             linearSolver.potentialError( psiCoeff ),
	             newtonSolver.potentialError( psiCoeff ),
	             relative, newtonSolver.newtonIterations() );
	std::fflush( stdout );

	BOOST_TEST( relative < 1.0e-10,
	            "the two setSource() overloads disagree by " << relative
	            << " relative on a source that does not depend on psi, so the "
	            "sign or the 1/r in SourceIntegrator does not match the right hand "
	            "side the linear path was measured with" );

	BOOST_TEST( newtonSolver.newtonIterations() <= 2,
	            "Newton took " << newtonSolver.newtonIterations() << " iterations on "
	            "a problem whose dF/dpsi is identically zero; one step should solve "
	            "it exactly" );
}

/*
 * THE RECONSTRUCTION DEFECT IS CONFINED TO THE RECONSTRUCTION, AND THIS IS THE
 * MEASUREMENT THAT SAYS SO.
 *
 * theReconstructionIsWrongWhereTheJacobianVanishes records that psi* comes back at
 * 8.98, 5.07, 6.90e11 and 4.77 when dF/dpsi vanishes identically. A number like
 * 6.9e11 is alarming enough that the first question to ask is whether the SOLVE
 * is also wrong somewhere less visible -- whether the same degeneracy reaches
 * the element-local elimination that produces psi_h and q_h, and is merely
 * better hidden there.
 *
 * IT DOES NOT, and the reason is structural before it is measured. The forward
 * solve's local problem takes its potential block from the HDG stabilisation,
 * which lives on Mnl_p's INTERIOR FACE integrators and is a fixed bilinear form
 * whatever the source does. Reconstruct() builds a DIFFERENT local problem, and
 * for its potential mass it lifts Mnl_p's DOMAIN integrators -- which carry the
 * reaction term -dF/dpsi/r and nothing else, so they vanish with it. Two
 * problems, two operators, one of which degenerates. And Reconstruct() writes
 * only into the post-processing GridFunctions: potentialGf and fluxGf are
 * already final when postProcess() is entered, and are not arguments to it.
 *
 * That is an argument. What follows is the measurement, over four degrees and
 * three meshes: the SAME source, once as a coefficient on the linear path and
 * once as a meq::Source on the Newton path, compared in BOTH psi_h and q_h. If
 * the degeneracy reached the forward solve at all, the two columns would part.
 *
 * q_h is in here deliberately and is not decoration. It is the physically
 * interesting output -- B_poloidal is a relabelling of it -- it is what the
 * mixed method exists to deliver at k+1, and it comes out of the same local
 * elimination as psi_h. A defect that moved q while leaving psi alone is exactly
 * the sort of thing a psi-only check would wave through.
 *
 * AND THE THIRD CHECK IS THAT A BROKEN postProcess() LEAVES THE ANSWER ALONE.
 * The driver post-processes and then writes potential() and flux(); if the
 * reconstruction had scribbled on either of them, it would be putting a
 * corrupted equilibrium on disk while reporting success.
 */
BOOST_AUTO_TEST_CASE( theReconstructionDefectDoesNotReachTheSolve )
{
	meq::analytic::SolovievEquilibrium const eq = meq::analytic::SolovievEquilibrium::nstx();
	meq::SolovievSource const source( eq.getA() );

	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );
	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );

	std::printf( "\n  the same Solov'ev source down both paths, psi_h and q_h\n" );
	std::printf( "    %2s %4s %14s %14s %14s\n",
	             "k", "n", "rel psi", "rel q", "L2(psi) error" );

	double worstPsi = 0.0;
	double worstFlux = 0.0;

	for ( int order = 1; order <= 4; ++order )
	{
		for ( int n : { 4, 8, 16 } )
		{
			mfem::Mesh linearMesh = makeMesh( n );
			meq::GradShafranovSolver linear( linearMesh, order );
			linear.setSource( sourceCoeff );
			linear.setBoundaryData( psiCoeff );
			linear.solve();

			mfem::Mesh newtonMesh = makeMesh( n );
			meq::GradShafranovSolver newton( newtonMesh, order );
			newton.setSource( source );
			newton.setBoundaryData( psiCoeff );
			newton.solve();

			BOOST_TEST_REQUIRE( linear.isNonlinear() == false );
			BOOST_TEST_REQUIRE( newton.isNonlinear() == true );

			mfem::Vector psiDifference( newton.potential() );
			psiDifference -= linear.potential();
			double const relativePsi =
				psiDifference.Norml2()/std::max( 1.0e-300, linear.potential().Norml2() );

			mfem::Vector fluxDifference( newton.flux() );
			fluxDifference -= linear.flux();
			double const relativeFlux =
				fluxDifference.Norml2()/std::max( 1.0e-300, linear.flux().Norml2() );

			std::printf( "    %2d %4d %14.3e %14.3e %14.4e\n",
			             order, n, relativePsi, relativeFlux,
			             linear.potentialError( psiCoeff ) );
			std::fflush( stdout );

			worstPsi = std::max( worstPsi, relativePsi );
			worstFlux = std::max( worstFlux, relativeFlux );

			// The Newton path is the one whose reconstruction is degenerate. Its
			// psi_h and q_h must still be the linear path's, which have been
			// measured against the closed form at k+1 in SolovievConvergence.cpp.
			BOOST_TEST( relativePsi < 1.0e-9,
			            "k = " << order << ", n = " << n << ": psi_h differs by "
			            << relativePsi << " between the paths" );
			BOOST_TEST( relativeFlux < 1.0e-9,
			            "k = " << order << ", n = " << n << ": q_h differs by "
			            << relativeFlux << " between the paths -- the flux is the "
			            "physical output and it has moved where psi may not have" );
		}
	}

	std::printf( "    worst over the sweep: psi %.3e, q %.3e\n", worstPsi, worstFlux );
	std::fflush( stdout );

	// AND THE BROKEN RECONSTRUCTION LEAVES BOTH OF THEM WHERE THEY WERE.
	{
		mfem::Mesh mesh = makeMesh( 8 );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( psiCoeff );
		solver.solve();

		mfem::Vector const potentialBefore( solver.potential() );
		mfem::Vector const fluxBefore( solver.flux() );

		// This is the degenerate case -- Solov'ev has dF/dpsi == 0 -- so what
		// comes back is not a post-processed potential. It must still not have
		// touched what was already computed.
		solver.postProcess();

		mfem::Vector potentialAfter( solver.potential() );
		potentialAfter -= potentialBefore;
		mfem::Vector fluxAfter( solver.flux() );
		fluxAfter -= fluxBefore;

		BOOST_TEST( potentialAfter.Normlinf() == 0.0,
		            "postProcess() moved psi_h by " << potentialAfter.Normlinf()
		            << ". The driver writes the potential after post-processing "
		            "has run, so a reconstruction that scribbles on it would put a "
		            "corrupted equilibrium on disk while reporting success" );
		BOOST_TEST( fluxAfter.Normlinf() == 0.0,
		            "postProcess() moved q_h by " << fluxAfter.Normlinf() );
	}
}

/// The middle rung of the analytic ladder, and the sharpest single statement
/// stage 4 can make.
///
/// McCarthy's source is F = T ( psi - c1 - c2 r^2 ): linear in psi, so
/// dF/dpsi = T identically and the whole discrete residual is affine in the
/// unknowns. An exact Jacobian therefore solves it in ONE Newton step, exactly,
/// and this asserts that.
///
/// It is what CLAUDE.md's ladder puts this fixture in the tree for. Solov'ev has
/// dF/dpsi = 0, so it cannot tell a correct Jacobian from an absent one; Example
/// 5's dF/dpsi is messy enough to hide a factor. Here the Jacobian's mass term
/// is a single constant, and it is either present or it is not: drop it and the
/// iteration becomes a chord method that still lands on the same answer, but
/// takes many steps to do it.
///
/// Measured, with meq::SourceIntegrator::AssembleElementGrad() returning zero:
/// 87, 84, 80 and 76 iterations down the four meshes instead of one, with the
/// residual settling to a constant contraction of 0.7711 per step -- linear
/// convergence, exactly the shape CLAUDE.md warns about -- and every error in
/// the table below unchanged to six significant figures. That is the whole
/// argument for asserting on the iteration and not only on the answer.
///
/// The spatial rates are checked at the same time, because an exact solution
/// with a psi-dependent source is not otherwise available: this is the only
/// benchmark here that is both non-trivial for Newton and closed form.
BOOST_AUTO_TEST_CASE( mccarthyNeedsExactlyOneNewtonStep )
{
	meq::analytic::McCarthyEquilibrium const eq
		= meq::analytic::McCarthyEquilibrium::asdex();

	// The fixture's own transcription first, on this box. Eighteen Bessel and
	// Neumann terms is a lot to get right by eye, and a mistyped one would
	// present as a solver that converges to the wrong equilibrium.
	for ( double r = rMin; r <= rMax + 1.0e-12; r += 0.2 )
	{
		for ( double z = zMin; z <= zMax + 1.0e-12; z += 0.3 )
		{
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, eq.psi( r, z ) );
			BOOST_TEST( std::abs( deltaStar - minusF ) < 1.0e-5,
			            "at ( " << r << ", " << z << " ): Delta*(psi) = " << deltaStar
			            << " but -F = " << minusF );
		}
	}

	BOOST_TEST( eq.dFdPsi( 1.0, 0.0, 0.0 ) != 0.0,
	            "McCarthy's dF/dpsi is zero, so this fixture has become another "
	            "Solov'ev and tests nothing Solov'ev does not" );

	int const order = 2;
	std::vector<int> const sizes = { 4, 8, 16, 32 };

	std::vector<Measurement> points;
	points.reserve( sizes.size() );
	for ( int n : sizes )
		points.push_back( measure( eq, order, n ) );

	printTable( "McCarthy ASDEX Upgrade", order, points );

	double const expected = order + 1.0 - rateSlack;

	for ( std::size_t i = 0; i < points.size(); ++i )
	{
		BOOST_TEST( points[ i ].newtonIterations == 1,
		            "h = " << points[ i ].h << ": Newton took "
		            << points[ i ].newtonIterations << " iterations on a source that "
		            "is affine in psi. One exact step should finish it; more than one "
		            "means the Jacobian is not the derivative of the residual" );

		if ( i == 0 )
			continue;

		double const ratio = points[ i - 1 ].h/points[ i ].h;
		double const ratePsi = rate( points[ i - 1 ].errorPsi, points[ i ].errorPsi, ratio );
		double const rateFlux = rate( points[ i - 1 ].errorFlux, points[ i ].errorFlux, ratio );

		BOOST_TEST( ratePsi >= expected,
		            "h = " << points[ i ].h << ": psi converged at " << ratePsi
		            << ", wanted " << expected );
		BOOST_TEST( rateFlux >= expected,
		            "h = " << points[ i ].h << ": q converged at " << rateFlux
		            << ", wanted " << expected );
	}
}

/// Newton's residual history, asserted rather than admired.
///
/// The assertion is on the observed order log( r2/r1 )/log( r1/r0 ), which is 2
/// for a quadratically convergent iteration in its asymptotic regime and 1 for a
/// linearly convergent one anywhere. Only triples above the round-off floor are
/// considered: once the residual reaches machine precision relative to the first
/// one, the ratios are noise and the order estimated from them is meaningless.
/*
 * EVERY NON-LINEAR PATH REACHES THE SAME EXACT SOLUTION.
 *
 * meq now has four ways to close the semi-linear problem, and they have
 * genuinely different nonlinear structure -- Newton puts F on the non-linear
 * potential mass, where hybridization makes each element's elimination its own
 * Newton; the two Picard paths freeze F at the previous iterate and put it on
 * the right hand side, so the potential block is linear and no element-local
 * non-linear solve exists at all. PicardThenNewton is BOTH, in sequence, and it
 * rebuilds the forms in the middle of one solve() to change from one to the
 * other -- which is the reason it is pinned here rather than trusted to be a
 * composition of two paths that are pinned separately. See CLAUDE.md.
 *
 * Agreeing with each other would not be enough -- three routes through the same
 * wrong sign would agree perfectly. Example 5 has a closed form, so this
 * measures all three against it: same mesh, same degree, same Dirichlet data,
 * and the L2 error must come out at the DISCRETISATION error, identically,
 * because the discretisation is identical and only the path to it differs.
 *
 * This is what would catch FrozenSource getting the -1/r or the sign convention
 * wrong. Those are applied in setSource( Coefficient & ) for the linear path and
 * re-applied for the frozen one; a mismatch would leave the Picard paths
 * converging beautifully to the wrong function, which is precisely the failure
 * mode this project's testing stance exists for.
 */
BOOST_AUTO_TEST_CASE( everyNonlinearPathReachesTheSameExactSolution )
{
	using G = meq::GradShafranovSolver::Globalisation;
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );

	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	struct Run { char const *name; G globalisation; double damping; };
	Run const runs[] = {
		{ "Newton          ", G::None,             1.0 },
		{ "Anderson-Picard ", G::AndersonPicard,   1.0 },
		{ "Picard, w = 0.5 ", G::PicardOnly,       0.5 },
		// The path the driver ships on an observed failure, and the one that
		// switches forms mid-solve: setGlobalisation() reset `prepared` but not
		// `built` until stage 4 found it, and this is the configuration that
		// would have caught it -- stage 2 would have reused stage 1's blocks and
		// converged to something else. Reported iterations are stage 2's.
		{ "Picard->Newton  ", G::PicardThenNewton, 1.0 },
	};

	double reference = -1.0;
	std::printf( "\n  Example 5 at k = 2, n = 16, by four non-linear paths\n" );

	for ( Run const &run : runs )
	{
		mfem::Mesh mesh = makeMesh( 16 );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( exact );
		solver.setGlobalisation( run.globalisation );
		solver.setPicardDamping( run.damping );
		// Generous: the Picard paths need hundreds where Newton needs four, and
		// the point here is where they land, not how fast.
		solver.setNewtonControl( 1.0e-10, 1.0e-14, 500 );

		BOOST_REQUIRE_NO_THROW( solver.solve() );

		double const error = solver.potentialError( exact );
		std::printf( "    %s %3d iterations, L2 error %.6e\n",
		             run.name, solver.newtonIterations(), error );
		std::fflush( stdout );

		// MEASURED, and it is the reason this run is noisy on stdout: stage 2 of
		// PicardThenNewton emits 2164 "el: N not convered in 100 iters" warnings
		// from MFEM's element-local non-linear solves, where plain Newton on the
		// same problem emits NONE. Starting Newton from Picard's iterate drives
		// the local problems harder than starting it from the Dirichlet datum
		// does -- the same effect CLAUDE.md measures for every KINSOL strategy
		// under "Why meq's Newton struggles", where it is fatal. Here it is not:
		// the outer iteration still finishes in two steps at the same L2 error to
		// seven figures, which is what the assertion below holds it to. So a
		// local solve giving up is not by itself a failed solve, and the warnings
		// are not a reason to reach for a tolerance.

		if ( reference < 0.0 )
			reference = error;
		else
			BOOST_TEST( std::abs( error - reference ) <= 1.0e-6*reference,
			            "path " << run.name << " reached L2 error " << error
			            << " where Newton reached " << reference
			            << " -- the same discretisation should give the same error "
			            "whatever iteration found it" );
	}

	// And it must be the discretisation error, not some larger number all four
	// happen to share. k = 2 on this mesh: measured 2.983495e-5, by all four
	// paths, to every digit printed.
	BOOST_TEST( reference < 1.0e-4,
	            "the common error is " << reference << ", too large to be the "
	            "k = 2 discretisation error on this mesh" );
}

/*
 * PSI* THROUGH NEWTON. An MFEM defect used to make this impossible, and the
 * measurement is what says it does not any more.
 *
 * DarcyForm::ReconstructFluxAndPot() consulted only the LINEAR potential mass
 * form M_p. meq's Newton path puts the whole potential block on Mnl_p -- see
 * buildForms() -- so the local reconstruction problem got no potential mass and
 * no face constraint, was singular, and was factored and solved anyway.
 * Measured, it returned psi* of 9.9e14, 8.4e15 and 3.9e14 at k = 1, 2, 3 where
 * the same problem solved linearly gives 3.8e-6, 2.4e-7 and 1.5e-8 -- and it
 * returned them WITHOUT A WORD, with psi_h agreeing to six figures either way.
 * postProcess() therefore threw rather than passing that on.
 *
 * The fix is on the branch meq builds from: darcyform.cpp now lifts the
 * non-linear potential integrators as a Jacobian frozen at the computed
 * potential, and verifies that the face constraint's gradient does not depend on
 * the trace -- which meq's constant tau satisfies.
 *
 * SO THE REFUSAL IS RETIRED, AND THIS IS WHAT RETIRED IT. Not reading the fix:
 * the failure mode was a silent 1e15, so only a rate can tell the difference
 * between post-processing and a singular local solve that happened to return
 * something finite. k+2 is asserted, and separately that psi* actually beats
 * psi_h -- a plausible-looking psi* that is no better than the potential it came
 * from would be the quiet version of the same defect.
 *
 * IT MATTERS BEYOND TIDINESS: eq (20) builds FOUR of its five terms on psi*, so
 * with the refusal in place the adaptive loop was linear-only, and the driver
 * would have had to run the estimator on ResidualEstimator::Potential::Raw --
 * which CLAUDE.md measures losing exactly one order at every k.
 */
BOOST_AUTO_TEST_CASE( thePostProcessedPotentialSurvivesNewton )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> const source( eq );

	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	struct Point { double h; double errorPsi; double errorStar; int iterations; };

	for ( int order = 1; order <= 3; ++order )
	{
		std::vector<Point> points;
		for ( int n : meshSizes )
		{
			mfem::Mesh mesh = makeMesh( n );
			meq::GradShafranovSolver solver( mesh, order );
			solver.setSource( source );
			solver.setBoundaryData( exact );
			solver.solve();
			// The whole point: this used to throw, and before it threw it
			// returned 1e15.
			solver.postProcess();

			points.push_back( Point{ ( rMax - rMin )/n,
			                         solver.potentialError( exact ),
			                         solver.postProcessedPotentialError( exact ),
			                         solver.newtonIterations() } );
		}

		std::printf( "\n  psi* through Newton, Example 5, k = %d\n", order );
		std::printf( "    %8s %5s %12s %6s %12s %6s\n",
		             "h", "it", "L2(psi)", "rate", "L2(psi*)", "rate" );
		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			Point const &p = points[ i ];
			std::printf( "    %8.5f %5d %12.4e", p.h, p.iterations, p.errorPsi );
			if ( i == 0 )
			{
				std::printf( " %6s %12.4e %6s\n", "-", p.errorStar, "-" );
			}
			else
			{
				double const ratio = points[ i - 1 ].h/p.h;
				std::printf( " %6.3f %12.4e %6.3f\n",
				             rate( points[ i - 1 ].errorPsi, p.errorPsi, ratio ),
				             p.errorStar,
				             rate( points[ i - 1 ].errorStar, p.errorStar, ratio ) );
			}
		}
		std::fflush( stdout );

		for ( std::size_t i = 1; i < points.size(); ++i )
		{
			double const ratio = points[ i - 1 ].h/points[ i ].h;
			double const measured = rate( points[ i - 1 ].errorStar,
			                              points[ i ].errorStar, ratio );
			BOOST_TEST( measured >= order + 2.0 - rateSlack,
			            "k = " << order << ", h = " << points[ i ].h
			            << ": psi* converged at " << measured << " on the NEWTON "
			            "path, wanting " << order + 2.0 - rateSlack
			            << ". If this has gone, DarcyForm::Reconstruct() is back to "
			            "reading only the linear potential mass and postProcess() "
			            "must refuse the non-linear path again" );
		}

		// A rate alone would not catch a psi* that converges at k+2 to something
		// no better than psi_h, and the defect this replaces was silent.
		BOOST_TEST( points.back().errorStar < 0.1*points.back().errorPsi,
		            "k = " << order << ": psi* is " << points.back().errorStar
		            << " against " << points.back().errorPsi << " for psi_h, so the "
		            "post-processing is not buying the extra order it converges at" );
	}
}

/*
 * AND WHERE dF/dpsi VANISHES, psi* IS STILL WRONG -- WHICH IS WHY THE DRIVER'S
 * ADAPTIVE LOOP DOES NOT USE IT.
 *
 * THIS TEST ASSERTS A DEFECT. It is written to FAIL the day MFEM fixes it, and
 * the failure message says what to do; do not weaken it to keep it green.
 *
 * THE MECHANISM, which was worth locating exactly because the first two guesses
 * at it were wrong. DarcyForm::ReconstructFluxAndPot() closes the local
 * post-processing problem with a mean-value constraint -- the branch commented
 * "adjust the element average of potential" -- because that problem is a pure
 * NEUMANN problem and is singular by construction, determined only up to a
 * constant. Correct, and exactly what the method needs. But that branch is
 * skipped whenever nl_src is set, and nl_src is set on the mere PRESENCE of a
 * non-convection non-linear integrator:
 *
 *     nlfi->AssembleElementGrad( *fe_p, *Tr, p_lift, Mp_k );  // zero if dF/du == 0
 *     Mp_z += Mp_k;
 *     if ( !convection-like ) { Mp_k.AddMult( p_lift, rhs_p_nl ); nl_src = true; }
 *
 * So where dF/dpsi vanishes, Mp_z and rhs_p_nl are both zero, the regularisation
 * is skipped anyway, and inv.Factor() factors a singular matrix -- with
 * DenseMatrixInverse::Factor( m, TOL = 0.0 ) testing abs( pivot ) <= 0.0 and the
 * return value discarded. It is NOT that a zero right hand side is
 * unhandleable: the singularity is by design and the handling exists. It is
 * unreachable.
 *
 * IT IS PER ELEMENT, because nl_src is per element. That is the part that makes
 * this worth a test rather than a footnote: a profile with a flat segment gives
 * dF/dpsi == 0 on part of the domain and nowhere else, which is ordinary rather
 * than pathological, and it corrupts psi* on exactly those elements -- which are
 * then exactly the elements an estimator built on psi* would mark.
 *
 * meq carries no runtime check for this. A solver should not stand permanently
 * on guard against its own dependency, and the condition is one flag test away
 * from never arising. This test is the record of its state instead. Written up
 * as ../mfem-hdg-dev/doc/HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md.
 */
BOOST_AUTO_TEST_CASE( theReconstructionIsWrongWhereTheJacobianVanishes )
{
	// F depends on psi only where z > threshold. Below it dF/dpsi is exactly
	// zero -- what a pressure profile with a flat segment does. `floor` adds a
	// constant to dF/dpsi and is the experiment: it changes nothing physical and
	// everything about whether a matrix is invertible.
	struct Layered : public meq::Source
	{
		Layered( double t, double e ) : threshold( t ), floor( e ) {}
		double f( double r, double z, double psi ) const override
		{
			return std::max( 0.0, z - threshold )*psi*psi + floor*psi + r*r;
		}
		double dFdPsi( double, double z, double psi ) const override
		{
			return 2.0*std::max( 0.0, z - threshold )*psi + floor;
		}
		double threshold;
		double floor;
	};

	mfem::FunctionCoefficient datum( []( mfem::Vector const &x )
	{
		return 0.2 + 0.1*x( 1 );
	} );

	/// max | psi* | against max | psi_h | on one element, from the dof values.
	/// GridFunction::ComputeElementL2Errors would be the natural call and is a
	/// trap: it MFEM_ASSERTs that the result vector is already sized, so in a
	/// Release build it writes off the end of an empty one.
	auto elementRatios = []( meq::GradShafranovSolver &solver, mfem::Mesh &mesh,
	                         double threshold, double &worstDead, double &worstLive,
	                         int &dead )
	{
		worstDead = 0.0;
		worstLive = 0.0;
		dead = 0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::Vector starDofs, rawDofs, centre;
			solver.postProcessedPotential().GetElementDofValues( e, starDofs );
			solver.potential().GetElementDofValues( e, rawDofs );
			mesh.GetElementCenter( e, centre );

			double const ratio = starDofs.Normlinf()
			                     /std::max( 1.0e-300, rawDofs.Normlinf() );
			if ( centre( 1 ) < threshold ) { ++dead; worstDead = std::max( worstDead, ratio ); }
			else                            worstLive = std::max( worstLive, ratio );
		}
	};

	// ---- one: the corruption is local to the vanishing region ------------
	std::printf( "\n  psi* where dF/dpsi vanishes on part of the domain, k = 2, n = 8\n" );
	std::printf( "    %8s %10s %14s %14s\n",
	             "z <", "dead frac", "worst K dead", "worst K live" );

	for ( double threshold : { -0.4, -0.2, 0.0 } )
	{
		Layered const source( threshold, 0.0 );
		mfem::Mesh mesh = makeMesh( 8 );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( datum );
		solver.solve();
		solver.postProcess();

		double worstDead = 0.0, worstLive = 0.0;
		int dead = 0;
		elementRatios( solver, mesh, threshold, worstDead, worstLive, dead );

		std::printf( "    %8.2f %10.3f %14.4e %14.4e\n", threshold,
		             dead/static_cast<double>( mesh.GetNE() ), worstDead, worstLive );
		std::fflush( stdout );

		// MEASURED: 20.3, 64.1, 61.6 on the dead elements against 1.0049, 1.0024,
		// 1.0021 on the live ones. Post-processing perturbs the potential by
		// O( h^(k+1) ), so a live element sits within a per cent of 1 and a dead
		// one is off by more than an order of magnitude.
		BOOST_TEST( worstDead > 5.0,
		            "z < " << threshold << ": the worst element where dF/dpsi "
		            "vanishes has || psi* || / || psi_h || = " << worstDead
		            << ", which is no longer the corruption this test records. If "
		            "MFEM's ReconstructFluxAndPot() now reaches its mean-value "
		            "branch when the non-linear Jacobian is zero, then DELETE THIS "
		            "TEST and take setPotential( Potential::Raw ) out of "
		            "apps/meq.cpp -- the driver can use the published estimator "
		            "again" );
		BOOST_TEST( worstLive < 1.1,
		            "z < " << threshold << ": elements where dF/dpsi does NOT "
		            "vanish give " << worstLive << ", so the corruption is not "
		            "confined to the vanishing region and this test's account of "
		            "the mechanism is wrong" );
	}

	// ---- two: a floor on dF/dpsi, and nothing else, repairs it ------------
	//
	// This is what identifies the cause as a singular matrix rather than anything
	// about the physics. 1e-12 is twelve orders below every other quantity in the
	// problem: it cannot move a solution, and it moves this one from 7.57 to
	// 0.9985.
	std::printf( "\n  the same source with a floor on dF/dpsi\n" );
	std::printf( "    %12s %16s\n", "floor", "|psi*| / |psi|" );

	double ratioAtZero = 0.0;
	for ( double floor : { 0.0, 1.0e-12, 1.0e-8 } )
	{
		Layered const source( 0.0, floor );
		mfem::Mesh mesh = makeMesh( 8 );
		meq::GradShafranovSolver solver( mesh, 2 );
		solver.setSource( source );
		solver.setBoundaryData( datum );
		solver.solve();
		solver.postProcess();

		mfem::ConstantCoefficient nought( 0.0 );
		double const ratio = solver.postProcessedPotential().ComputeL2Error( nought )
		                     /solver.potential().ComputeL2Error( nought );
		std::printf( "    %12.0e %16.6f\n", floor, ratio );
		std::fflush( stdout );

		if ( floor == 0.0 )
		{
			ratioAtZero = ratio;
		}
		else
		{
			BOOST_TEST( std::abs( ratio - 1.0 ) < 0.01,
			            "a floor of " << floor << " on dF/dpsi gives "
			            << ratio << ", so an infinitesimal perturbation no longer "
			            "repairs psi* and the singular-matrix account is wrong" );
		}
	}

	BOOST_TEST( ratioAtZero > 2.0,
	            "with dF/dpsi vanishing over half the domain psi* now gives "
	            << ratioAtZero << ", which is not the corruption this test "
	            "records -- see the deletion instructions above" );
}

BOOST_AUTO_TEST_CASE( newtonConvergesQuadratically )
{
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();

	std::vector<double> history;
	Measurement const point = measure( eq, 3, 8, &history );

	BOOST_TEST( history.size() >= 3u,
	            "Newton produced only " << history.size() << " residuals; there is "
	            "nothing to estimate an order from. A problem this non-linear "
	            "should not be solved in one step" );

	std::printf( "\n  Newton residual history, k = 3, h = %.5f, %d trace dofs\n",
	             point.h, point.traceDofs );
	std::printf( "  %5s %16s %16s %8s\n", "it", "||r||", "||r||/||r_0||", "order" );
	for ( std::size_t i = 0; i < history.size(); ++i )
	{
		if ( i >= 2 )
		{
			std::printf( "  %5zu %16.6e %16.6e %8.3f\n", i, history[ i ],
			             history[ i ]/history[ 0 ],
			             newtonOrder( history[ i - 2 ], history[ i - 1 ], history[ i ] ) );
		}
		else
		{
			std::printf( "  %5zu %16.6e %16.6e %8s\n", i, history[ i ],
			             history[ i ]/history[ 0 ], "-" );
		}
	}
	std::fflush( stdout );

	// Round-off floor: below this the residual is not converging any more, it is
	// bouncing off the conditioning of the trace system.
	double const floor = 1.0e-12*history.front();

	double best = 0.0;
	for ( std::size_t i = 2; i < history.size(); ++i )
	{
		if ( history[ i ] < floor || history[ i - 1 ] < floor )
			continue;
		best = std::max( best, newtonOrder( history[ i - 2 ], history[ i - 1 ],
		                                    history[ i ] ) );
	}

	BOOST_TEST( best >= 1.8,
	            "the best observed Newton order is " << best << ", not the ~2 a "
	            "Jacobian consistent with its residual gives. A run that grinds "
	            "down linearly means the two disagree" );

	// A quadratically convergent iteration reaching 1e-12 relative takes a
	// handful of steps, not a dozen. Asserted separately from the order because
	// an iteration can show a good order over one triple and still be limping.
	BOOST_TEST( point.newtonIterations <= 8,
	            "Newton took " << point.newtonIterations << " iterations to reach "
	            "the tolerance; the reference shape is five" );

	// And the problem must actually be non-linear here, or the whole test is
	// measuring nothing. Two residuals would mean a single Newton step sufficed,
	// which is what the Solov'ev source gives and is why it cannot be used.
	BOOST_TEST( point.newtonIterations >= 2,
	            "Newton converged in " << point.newtonIterations << " iterations, so "
	            "the source is behaving linearly and dF/dpsi is untested" );
}

BOOST_AUTO_TEST_CASE( orderOneConvergesAtTwo )
{
	checkOrder( meq::analytic::ManufacturedNonlinear::example5(),
	            "Example 5 manufactured", 1, 1.5e-3, 2.5e-3 );
}

BOOST_AUTO_TEST_CASE( orderTwoConvergesAtThree )
{
	checkOrder( meq::analytic::ManufacturedNonlinear::example5(),
	            "Example 5 manufactured", 2, 1.2e-5, 2.0e-5 );
}

BOOST_AUTO_TEST_CASE( orderThreeConvergesAtFour )
{
	checkOrder( meq::analytic::ManufacturedNonlinear::example5(),
	            "Example 5 manufactured", 3, 7.0e-8, 1.2e-7 );
}


/*
 * THE SIMILARITY-REDUCTION BENCHMARK, and why it is worth having beside
 * Example 5.
 *
 * Example 5 is manufactured: a psi was chosen and F built to fit it, so its F
 * is not of a form any physical profile produces. The solution below is the
 * other way round -- the free function is chosen, f(u) = f0 exp(n u), and the
 * solution follows from a similarity reduction of the Grad-Shafranov equation
 * (Kaltsas & Throumoulopoulos, Phys. Lett. A 380 (2016) 3373, eq (22);
 * refs/GS-SimilarityReduction.pdf). So this tests the solver against a
 * nonlinear equation somebody might pose rather than one reverse engineered
 * from an answer.
 *
 * It is also a different shape of nonlinearity. Example 5 mixes psi^2 with
 * exp(-psi) and a geometric term; here the entire psi-dependence is one
 * exponential, so dF/dpsi = n F exactly. Over the benchmark rectangle at n = 3
 * that exponential varies by a factor of 29, which is a real nonlinearity
 * rather than a perturbation.
 */

/// The same guard as manufacturedSourceMatchesTheOperator, on the new fixture:
/// if -Delta*(psi) does not equal F at the exact solution, everything measured
/// against it is measured against the wrong problem, and no convergence rate
/// can see it.
BOOST_AUTO_TEST_CASE( similaritySourceMatchesTheOperator )
{
	meq::analytic::SimilarityExponential const eq
		= meq::analytic::SimilarityExponential::benchmark();

	double worst = 0.0;
	for ( double r = rMin + 0.05; r < rMax; r += 0.05 )
	{
		for ( double z = zMin + 0.05; z < zMax; z += 0.05 )
		{
			double const deltaStar = eq.deltaStarFD( r, z );
			double const minusF = -eq.f( r, z, eq.psi( r, z ) );
			worst = std::max( worst, std::abs( deltaStar - minusF )/std::abs( minusF ) );
		}
	}

	BOOST_TEST_MESSAGE( "  similarity solution: -Delta*(psi) vs F, worst relative "
	                    "difference " << worst );
	BOOST_TEST( worst < 1.0e-5,
	            "the similarity solution does not solve the equation: worst "
	            "relative difference " << worst );
}

/// dF/dpsi = n F for an exponential free function. Checked against a difference
/// of f() rather than assumed, because a Jacobian error is invisible to every
/// convergence rate -- see the +5% experiment recorded in CLAUDE.md.
BOOST_AUTO_TEST_CASE( similarityDerivativeMatchesAFiniteDifference )
{
	meq::analytic::SimilarityExponential const eq
		= meq::analytic::SimilarityExponential::benchmark();

	double const step = 1.0e-7;
	double worst = 0.0;
	for ( double r = rMin + 0.07; r < rMax; r += 0.07 )
	{
		for ( double psi = -0.8; psi < 0.45; psi += 0.15 )
		{
			double const difference =
				( eq.f( r, 0.11, psi + step ) - eq.f( r, 0.11, psi - step ) )/( 2.0*step );
			worst = std::max( worst,
				std::abs( eq.dFdPsi( r, 0.11, psi ) - difference )/std::abs( difference ) );
		}
	}

	BOOST_TEST_MESSAGE( "  similarity dF/dpsi vs central difference: worst relative "
	                    "difference " << worst );
	BOOST_TEST( worst < 1.0e-6 );
}

BOOST_AUTO_TEST_CASE( similarityOrderOneConvergesAtTwo )
{
	// Ceilings at 3x measured: psi 1.1024e-04, q 2.3063e-04 at h = 0.025.
	checkOrder( meq::analytic::SimilarityExponential::benchmark(),
	            "Kaltsas-Throumoulopoulos exponential", 1, 3.4e-4, 7.0e-4 );
}

BOOST_AUTO_TEST_CASE( similarityOrderTwoConvergesAtThree )
{
	// Measured: psi 8.6282e-07, q 2.0018e-06.
	checkOrder( meq::analytic::SimilarityExponential::benchmark(),
	            "Kaltsas-Throumoulopoulos exponential", 2, 2.6e-6, 6.1e-6 );
}

BOOST_AUTO_TEST_CASE( similarityOrderThreeConvergesAtFour )
{
	// Measured: psi 1.0006e-08, q 2.4542e-08.
	checkOrder( meq::analytic::SimilarityExponential::benchmark(),
	            "Kaltsas-Throumoulopoulos exponential", 3, 3.1e-8, 7.4e-8 );
}

/*
 * THE SYMBOLIC ANALYSIS IS COMPUTED ONCE, NOT ONCE PER NEWTON STEP.
 *
 * UMFPACK splits a solve into a symbolic analysis, which depends only on the
 * sparsity pattern, and a numeric factorisation, which depends on the values.
 * NewtonSolver::Mult calls prec->SetOperator( *grad ) once per iteration, and
 * the hybridized trace system's pattern never changes between them -- so
 * without SetReuseSymbolic() the analysis is recomputed and discarded every
 * step. Measured before it was enabled, that was 22 to 24% of each step, and
 * meq asks for METIS ordering, which makes it dearer than the default.
 *
 * This asserts the RATIO rather than a time, which is the house rule: a timing
 * would be a measurement about this machine, where the count is a measurement
 * about the code. One analysis, and one factorisation per iteration.
 *
 * It is also the only thing that would notice the reuse silently lapsing. The
 * pattern check inside SetReuseSymbolic() is exact and re-analyses whenever it
 * fails, so a lapse costs performance and NOTHING ELSE -- no wrong answer, no
 * failed convergence, nothing a rate table or an error norm could see.
 */
BOOST_AUTO_TEST_CASE( theSymbolicAnalysisIsReusedAcrossNewtonSteps )
{
#ifndef MFEM_USE_SUITESPARSE
	BOOST_TEST_MESSAGE( "no SuiteSparse, so there is no UMFPACK analysis to reuse" );
#else
	meq::analytic::ManufacturedNonlinear const eq
		= meq::analytic::ManufacturedNonlinear::example5();
	EquilibriumSource<meq::analytic::ManufacturedNonlinear> source( eq );

	mfem::Mesh mesh = makeMesh( 8 );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( psiCoeff );
	solver.setNewtonControl( 1.0e-10, 1.0e-14, 60 );
	solver.solve();

	std::printf( "\n  Newton took %d iterations: %ld symbolic analyses, "
	             "%ld numeric factorisations\n",
	             solver.newtonIterations(), solver.symbolicFactorisations(),
	             solver.numericFactorisations() );
	std::fflush( stdout );

	// More than one step, or the test cannot distinguish reuse from there being
	// nothing to reuse.
	BOOST_TEST_REQUIRE( solver.newtonIterations() >= 2 );

	BOOST_TEST( solver.symbolicFactorisations() == 1,
	            "the sparsity was analysed " << solver.symbolicFactorisations()
	            << " times across " << solver.newtonIterations()
	            << " Newton iterations, and once is the point of "
	            "SetReuseSymbolic(). This costs only speed, which is why nothing "
	            "else in the suite would catch it" );

	// Reuse must not have skipped a refactorisation: the values change every
	// step even though the pattern does not, so each iteration needs its own.
	BOOST_TEST( solver.numericFactorisations() > solver.symbolicFactorisations(),
	            "only " << solver.numericFactorisations() << " numeric "
	            "factorisations for " << solver.newtonIterations()
	            << " iterations -- reuse must never skip one of these" );
#endif
}
