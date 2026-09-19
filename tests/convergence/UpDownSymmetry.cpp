#define BOOST_TEST_MODULE UpDownSymmetry
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

#include "convergence/ConvergenceHarness.hpp"

/*
 * UP-DOWN SYMMETRY AS A CONSTRAINT, AND THE MESH IT NEEDS.
 *
 * `[solver] UpDownSymmetry` projects the iterate onto the subspace of fields
 * even in z, so a DOUBLE NULL's two saddles are exactly degenerate when the
 * X-point search looks for them. It exists so that a double-null machine need
 * not be told WHICH null to follow, which is a choice nothing physical makes.
 *
 * **IT IS A DOF-FOR-DOF AVERAGE AND SO IT NEEDS A MIRROR-SYMMETRIC MESH.**
 * Every dof is paired with the dof at its own reflection and the two are
 * averaged; a dof with no partner has nothing to average against.
 * setUpDownSymmetry() therefore REFUSES such a mesh by name rather than
 * projecting onto something that is not a reflection, and the second case here
 * is that refusal.
 *
 * **AND THAT REFUSAL IS NOT HYPOTHETICAL: IT IS WHAT MAST-U MEETS.** Measured
 * on `examples/mastu-nke.msh`, **3624 of 4735 vertices have no mirror partner**
 * -- while the machine's own [[coils]] are mirror-paired to the last digit, one
 * pair excepted whose currents are equal and opposite at 1.9e-04 A against a
 * total of 2.3e+06 A. The GEOMETRY is symmetric and gmsh's triangulation of it
 * is not, which is a property of unstructured meshing rather than of the
 * machine. So the constraint is blocked on the MESHER, and the two cases here
 * are what says the projection itself is sound while that stands.
 *
 * MakeCartesian2D's TRIANGLES cannot be the symmetric mesh, which is the trap
 * to know before reaching for one: it splits every cell along ONE diagonal, so
 * a box symmetric in z has a triangulation that is not. LimiterCurve.cpp
 * records the same thing from the other side, its discrete limiter contact
 * converging to the midplane rather than sitting on it. The QUADRILATERAL
 * variant has no diagonal to choose and is symmetric whenever its box is.
 */

namespace
{
	/// The standard box, which is symmetric about z = 0, in QUADRILATERALS.
	mfem::Mesh symmetricMesh( int n )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::QUADRILATERAL, false, box.width(),
			box.height() );
		double const rMin = box.rMin;
		double const zMin = box.zMin;
		mesh.Transform( [ rMin, zMin ]( mfem::Vector const &in,
		                                mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	/// A symmetric bump, peaking on the midplane. It is the guess AND the
	/// thing that makes this problem have the symmetry being imposed: a bump
	/// off the midplane would break it in the one place the constraint acts.
	mfem::FunctionCoefficient evenBump( double height )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		double const rMin = box.rMin;
		double const width = box.width();
		double const depth = box.height();
		return mfem::FunctionCoefficient(
			[ height, rMin, width, depth ]( mfem::Vector const &x )
			{
				return height*std::sin( M_PI*( x( 0 ) - rMin )/width )
				       *std::cos( M_PI*x( 1 )/depth );
			} );
	}

	/// The bordered solve of HighBetaConvergence, on a mesh this file chooses
	/// and with every ingredient EVEN in z.
	///
	/// THE BORDERED PATH IS THE POINT AND NOT AN INCIDENTAL CHOICE. The
	/// projection is wired into solveWithNormalisation() and nowhere else,
	/// because psi_ax being an unknown is what makes the branch a question at
	/// all -- so a plain linear solve() does not touch it. THE FIRST VERSION
	/// OF THIS CASE USED ONE and read a step of exactly 0.000e+00 in every
	/// block at every ( k, n ), which is what a comparison of two identical
	/// runs reads. An exact zero from a routine that reassociates sums is not
	/// a pass, it is a routine that did not run.
	std::unique_ptr<meq::GradShafranovSolver> solveOn(
		mfem::Mesh &mesh, int order, bool symmetric,
		mfem::Coefficient &datum, mfem::Coefficient &guess,
		meq::NormalisedMHDSource &source )
	{
		meq::tests::Rectangle const box = meq::tests::standardBox();
		auto solver =
			std::make_unique<meq::GradShafranovSolver>( mesh, order );
		// ON THE MIDPLANE EXACTLY, so the boundary-flux constraint does not
		// itself break the symmetry it is being asked to keep.
		solver->setBoundaryFluxPoint( box.rMin + 0.68*box.width(), 0.0 );
		solver->setSource( source, 0.30 );
		solver->setBoundaryData( datum );
		solver->setInitialGuess( guess );
		solver->setUpDownSymmetry( symmetric );
		solver->setNewtonControl( 1.0e-12, 1.0e-14, 40 );
		solver->solve();
		return solver;
	}
}

/*
 * THE PROJECTION IS THE IDENTITY ON A PROBLEM THAT ALREADY HAS THE SYMMETRY,
 * AND THAT IS THE ASSERTION WITH TEETH.
 *
 * Both blocks are compared, and the FLUX is the one that can fail: psi is even
 * in z and its map carries no sign, but q = grad_bar( psi )/r has q_r EVEN and
 * q_z ODD, so the flux map has to carry a sign per COMPONENT. Get it backwards
 * and the projection lands on the ANTIsymmetric subspace, where the only even
 * field is zero -- a trivial branch reached silently, which is this tree's
 * most-recorded failure shape. Comparing whole solution blocks catches that in
 * one number where sampling psi alone would not catch it at all.
 */
BOOST_AUTO_TEST_CASE( theSymmetryProjectionLeavesASymmetricSolutionAlone )
{
	std::printf( "\n  THE PROJECTION AGAINST THE SAME BORDERED SOLVE WITHOUT IT\n" );
	std::printf( "    %2s %4s %9s %13s %13s %13s %13s\n",
	             "k", "n", "dofs", "psi_ax", "d psi_ax", "d psi", "d flux" );

	double worstPotential = 0.0;
	double worstFlux = 0.0;

	for ( int order : { 1, 2 } )
		for ( int n : { 8, 12 } )
		{
			auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
			auto ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
			meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );

			mfem::ConstantCoefficient zero( 0.0 );
			mfem::FunctionCoefficient guess = evenBump( 0.30 );

			mfem::Mesh plainMesh = symmetricMesh( n );
			mfem::Mesh projectedMesh = symmetricMesh( n );
			auto plain =
				solveOn( plainMesh, order, false, zero, guess, source );
			auto projected =
				solveOn( projectedMesh, order, true, zero, guess, source );

			mfem::Vector potentialStep( plain->potential() );
			potentialStep -= projected->potential();
			mfem::Vector fluxStep( plain->flux() );
			fluxStep -= projected->flux();

			double const scale = std::max( plain->potential().Normlinf(),
			                               1.0e-300 );
			double const dPotential = potentialStep.Normlinf()/scale;
			double const dFlux =
				fluxStep.Normlinf()/std::max( plain->flux().Normlinf(),
				                              1.0e-300 );
			double const dAxis =
				std::abs( plain->psiAxis() - projected->psiAxis() )
				/std::abs( plain->psiAxis() );

			worstPotential = std::max( worstPotential, dPotential );
			worstFlux = std::max( worstFlux, dFlux );

			std::printf( "    %2d %4d %9d %13.6e %13.3e %13.3e %13.3e\n",
			             order, n, plain->potential().Size(),
			             plain->psiAxis(), dAxis, dPotential, dFlux );
		}
	std::fflush( stdout );

	/*
	 * 1e-10 RATHER THAN ROUND-OFF. The two runs are not the same arithmetic:
	 * the projected one averages every dof with its partner before each
	 * residual, so it reassociates sums that the plain one does not, and it
	 * starts Newton from a projected guess. What is asserted is that the
	 * projection changes no digit that matters, not that it changes none.
	 */
	BOOST_TEST( worstPotential < 1.0e-10,
	            "the symmetry projection moved psi on a problem that is "
	            "already symmetric -- worst relative " << worstPotential
	            << ". On such a problem it must be the identity, so this is "
	            "the potential map pairing dofs that are not partners" );

	BOOST_TEST( worstFlux < 1.0e-10,
	            "the symmetry projection moved the FLUX on a problem that is "
	            "already symmetric -- worst relative " << worstFlux
	            << ". q_r is even in z and q_z is odd, so this is the SIGN the "
	            "flux map carries per component; with it backwards the "
	            "projection lands on the antisymmetric subspace, where the "
	            "only symmetric field is zero" );
}

/*
 * AND IT REFUSES A MESH THAT IS NOT MIRROR-SYMMETRIC, NAMING A NODE.
 *
 * This is the case MAST-U meets. The refusal is the feature: a dof with no
 * partner cannot be averaged, and averaging it against something else -- or
 * quietly leaving it alone -- would impose a constraint that is not a
 * reflection and report an equilibrium nobody asked for.
 */
BOOST_AUTO_TEST_CASE( theSymmetryProjectionRefusesAnAsymmetricMesh )
{
	meq::tests::Rectangle const box = meq::tests::standardBox();

	// The SAME box, which is symmetric about z = 0 -- so what is asymmetric
	// here is the triangulation and nothing else.
	mfem::Mesh mesh = meq::tests::makeMesh( box, 8 );

	auto pPrime = std::make_shared<meq::ConstantProfile const>( 0.45 );
	auto ggPrime = std::make_shared<meq::ConstantProfile const>( 0.30 );
	meq::NormalisedMHDSource source( pPrime, ggPrime, 1.0, 1.0 );
	mfem::ConstantCoefficient zero( 0.0 );
	mfem::FunctionCoefficient guess = evenBump( 0.30 );

	bool refused = false;
	std::string message;
	try
	{
		solveOn( mesh, 2, true, zero, guess, source );
	}
	catch ( std::exception const &error )
	{
		refused = true;
		message = error.what();
	}

	std::printf( "\n  A TRIANGULATED SYMMETRIC BOX: %s\n",
	             refused ? message.c_str() : "ACCEPTED" );
	std::fflush( stdout );

	BOOST_TEST( refused,
	            "MakeCartesian2D's triangulation splits every cell along one "
	            "diagonal, so this mesh is NOT symmetric in z even though its "
	            "box is -- and the projection accepted it. A dof with no "
	            "mirror partner has nothing to average against" );

	BOOST_TEST( message.find( "not up-down symmetric" ) != std::string::npos,
	            "the refusal does not say what is wrong with the mesh: \""
	            << message << "\". A caller who meets this needs to be told it "
	            "is the MESH and not the machine -- MAST-U's own conductors "
	            "are mirror-paired and its mesh is not" );
}
