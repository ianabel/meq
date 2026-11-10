#define BOOST_TEST_MODULE InnerSolveStall
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorTangent.hpp"
#include "meq/Source.hpp"

#include "SourceTree.hpp"

/*
 * RED, DELIBERATELY: THE BORDERED INNER NEWTON STALLS ON A g g' THAT
 * CONTINUATION REACHES (TODO, "The bordered inner Newton stalls on a sharp
 * g g' near the axis"; MEASUREMENTS.md M-398).
 *
 * The q-driven box at level 0, k = 2, g^2 the closed form of
 * QDrivenConsistency at the loop's own knots (meq::toroidalKnotCount()). The
 * innermost knot, Psi_N = 0.0083, alone raised by 0.3%: the end cubic carries
 * that onto the axis as g g' from 1.16 to 1.94. From the closed form's own
 * equilibrium the warm bordered Newton plateaus -- 60 iterations, the residual
 * stuck between 1e-2 and 4e-2 -- and so does a cold one; four quarter steps
 * of continuation reach it, in 28 iterations. An equilibrium exists and Newton
 * from its neighbour does not find it: a crawl, which is a defect.
 *
 * It is not the axis on a face row, which was suspected (M-392) and is
 * falsified: smooth changes of g g' that move the axis across the z = 0 row
 * converge in 5 or 6 iterations on this mesh as on one moved off it
 * (M-398). The q-driven loop meets the case as refused trial steps, which its
 * line search absorbs; the defect is the inner solve's.
 */
namespace
{
	double closedGSquared( double psi )
	{
		return 2.2*2.2 + ( std::exp( 2.0*psi ) - 1.0 )/( std::exp( 2.0 ) - 1.0 );
	}

	struct Box
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<meq::NormalisedMHDSource> source;
		std::unique_ptr<meq::GradShafranovSolver> solver;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<mfem::FunctionCoefficient> guess;
		mfem::Vector warm;

		explicit Box( std::shared_ptr<meq::Profile const> ggPrime )
		{
			mesh = std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian2D(
				8, 8, mfem::Element::TRIANGLE, false, 0.8, 1.2 ) );
			mesh->Transform( []( mfem::Vector const &in, mfem::Vector &out )
			{
				out = in;
				out( 0 ) += 0.6;
				out( 1 ) -= 0.6;
			} );
			auto pPrime = std::make_shared<meq::SplineProfile const>(
				meq::SplineProfile::fromFile( meq::test::sourcePath( "examples/q-driven-pprime.dat" ) ) );
			source = std::make_unique<meq::NormalisedMHDSource>( pPrime, ggPrime, 0.3, 1.0 );
			solver = std::make_unique<meq::GradShafranovSolver>( *mesh, 2 );
			guess = std::make_unique<mfem::FunctionCoefficient>( []( mfem::Vector const &x )
			{
				return 0.3*std::sin( M_PI*( x( 0 ) - 0.6 )/0.8 )*std::sin( M_PI*( x( 1 ) + 0.6 )/1.2 );
			} );
			solver->setSource( *source, 0.3 );
			solver->setBoundaryData( zero );
			solver->setInitialGuess( *guess );
			solver->setNewtonControl( 1.0e-12, 1.0e-14, 60 );
		}

		bool solve( std::shared_ptr<meq::Profile const> ggPrime )
		{
			source->setGGPrime( std::move( ggPrime ) );
			if ( warm.Size() > 0 )
				solver->setInitialState( warm, {} );
			try
			{
				solver->solve();
			}
			catch ( std::exception const & )
			{
				return false;
			}
			warm = solver->state();
			return true;
		}
	};
}

BOOST_AUTO_TEST_CASE( aWarmSolveReachesWhatContinuationReaches )
{
	mfem::Mesh probe = mfem::Mesh::MakeCartesian2D( 8, 8, mfem::Element::TRIANGLE, false, 0.8, 1.2 );
	std::size_t const knots = meq::toroidalKnotCount( meq::plasmaCells( probe ), 2 );
	meq::ToroidalFieldSpline spline( meq::toroidalKnotLabels( knots ) );
	std::vector<double> base;
	for ( double label : spline.labels() )
		base.push_back( closedGSquared( 1.0 - label ) );
	auto with = [ & ]( double factor )
	{
		std::vector<double> values( base );
		values[ 0 ] *= factor;
		meq::ToroidalFieldSpline s( spline );
		s.setValues( values );
		return s.ggPrime();
	};
	double const factor = 1.003;

	// The precondition: the equilibrium exists, reached by continuation.
	Box path( with( 1.0 ) );
	BOOST_TEST_REQUIRE( path.solve( with( 1.0 ) ) );
	bool reached = true;
	for ( int k = 1; k <= 4 && reached; ++k )
		reached = path.solve( with( 1.0 + ( factor - 1.0 )*k/4.0 ) );
	BOOST_TEST_REQUIRE( reached, "continuation does not reach the profile either, so this case "
	                             "is not the stall it is about" );

	Box warm( with( 1.0 ) );
	BOOST_TEST_REQUIRE( warm.solve( with( 1.0 ) ) );
	bool const solved = warm.solve( with( factor ) );
	int const iterations = warm.solver->newtonIterations();
	std::printf( "\n  THE INNERMOST KNOT x %.3f, g g' on the axis %.3f: warm Newton %s in %d iterations\n",
	             factor, ( *with( factor ) )( 1.0 ), solved ? "converged" : "did not converge", iterations );
	BOOST_TEST( ( solved && iterations <= 15 ),
		"the warm bordered Newton " << ( solved ? "converged" : "did not converge" ) << " in "
		<< iterations << " iterations from its neighbour's equilibrium, where four quarter steps "
		"of continuation reach the answer: a stall on a sharp g g' near the axis. TODO, \"The "
		"bordered inner Newton stalls on a sharp g g' near the axis\"; MEASUREMENTS.md M-398" );
}
