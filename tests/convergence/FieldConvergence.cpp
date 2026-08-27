#define BOOST_TEST_MODULE MeqFieldConvergence

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>

#include "mfem.hpp"

#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "ConvergenceHarness.hpp"
#include "analytic/Soloviev.hpp"

/*
 * B_poloidal, and the sign convention that no convergence rate can see.
 *
 * src/meq/Field.hpp derives B_R = -q_z, B_Z = +q_r from refs/Miller.pdf eq (1).
 * That is a derivation, and this project's record on convention questions
 * settled by derivation is one in four -- the Solov'ev source sign, tau in
 * eq (8e) and DarcyForm holding -q all went the other way. So the derivation is
 * not trusted here; it is measured, three ways:
 *
 *   1. against the analytic flux of the exact Solov'ev solution, which tests
 *      the relabelling and meq's flux() sign together;
 *   2. against CENTRAL DIFFERENCES of the exact psi, which tests the same thing
 *      without going through the fixture's own gradient -- so a transcription
 *      error in Soloviev.hpp's gradPsi() cannot make this pass;
 *   3. by convergence rate, because a field that is right at one mesh and
 *      converging at the wrong order is a different bug from a field with a
 *      wrong sign.
 *
 * WHY A WRONG SIGN WOULD OTHERWISE SURVIVE. Nothing in the suite reads B. The
 * solve is in psi and q, every rate is measured on psi and q, and B is produced
 * only for output. Flip a sign in Field.cpp and every existing test still
 * passes, the driver still writes a file, and every field line in it points the
 * wrong way. That is the failure this file exists for.
 */

namespace
{
	using meq::tests::SampleCloud;

	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq =
			meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	/// The benchmark rectangle of the Solov'ev studies.
	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}

	/// Solve, convert, and report the worst error in B over a sample cloud.
	struct FieldMeasurement
	{
		double h;
		double worstAgainstAnalytic;
		double worstAgainstDifference;
		double magnitude;
	};

	FieldMeasurement measure( int order, int n )
	{
		meq::analytic::SolovievEquilibrium const &eq = equilibrium();
		mfem::Mesh mesh = meq::tests::makeMesh( box(), n );

		mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{
			return eq.f( x( 0 ), x( 1 ), 0.0 );
		} );
		mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{
			return eq.psi( x( 0 ), x( 1 ) );
		} );

		meq::GradShafranovSolver solver( mesh, order );
		solver.setSource( source );
		solver.setBoundaryData( exact );
		solver.solve();

		mfem::GridFunction field;
		meq::poloidalField( solver.flux(), field );

		// Sampled well inside the box: B is compared pointwise and the boundary
		// layer of a DG field is where the jumps live.
		FieldMeasurement out{ box().width()/static_cast<double>( n ), 0.0, 0.0, 0.0 };
		double const inset = 0.08;
		double const step = 1.0e-5;

		for ( double r = box().rMin + inset; r <= box().rMax - inset; r += 0.05 )
			for ( double z = box().zMin + inset; z <= box().zMax - inset; z += 0.05 )
			{
				mfem::Vector point( 2 );
				point( 0 ) = r;
				point( 1 ) = z;

				mfem::DenseMatrix points( 2, 1 );
				points( 0, 0 ) = r;
				points( 1, 0 ) = z;
				mfem::Array<int> elements;
				mfem::Array<mfem::IntegrationPoint> ips;
				if ( mesh.FindPoints( points, elements, ips ) < 1 || elements[ 0 ] < 0 )
					continue;

				mfem::Vector computed( 2 );
				field.GetVectorValue( elements[ 0 ], ips[ 0 ], computed );

				// (1) the fixture's analytic flux, relabelled the same way.
				double qR = 0.0, qZ = 0.0;
				eq.flux( r, z, qR, qZ );
				double const analyticBr = -qZ;
				double const analyticBz = qR;

				// (2) central differences of the exact psi, which does not use
				// the fixture's gradient at all.
				double const dPsiDr = ( eq.psi( r + step, z ) - eq.psi( r - step, z ) )
				                      /( 2.0*step );
				double const dPsiDz = ( eq.psi( r, z + step ) - eq.psi( r, z - step ) )
				                      /( 2.0*step );
				double const differenceBr = -dPsiDz/r;
				double const differenceBz = dPsiDr/r;

				out.worstAgainstAnalytic = std::max( out.worstAgainstAnalytic,
					std::hypot( computed( 0 ) - analyticBr, computed( 1 ) - analyticBz ) );
				out.worstAgainstDifference = std::max( out.worstAgainstDifference,
					std::hypot( computed( 0 ) - differenceBr, computed( 1 ) - differenceBz ) );
				out.magnitude = std::max( out.magnitude,
					std::hypot( analyticBr, analyticBz ) );
			}

		return out;
	}
}

BOOST_AUTO_TEST_SUITE( field_convergence )

/// The relabelling itself, with no solve involved: B from a flux GridFunction
/// that was set by hand. This is what fails first if the component swap or the
/// sign is wrong, and it fails unambiguously.
BOOST_AUTO_TEST_CASE( theRelabellingSwapsComponentsAndFlipsOneSign )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 2 );
	mfem::L2_FECollection collection( 1, mesh.Dimension() );
	mfem::FiniteElementSpace space( &mesh, &collection, 2 );

	mfem::GridFunction q( &space );
	// q_r = 3 everywhere, q_z = -7 everywhere, whatever the ordering.
	mfem::Vector componentValue( 2 );
	componentValue( 0 ) = 3.0;
	componentValue( 1 ) = -7.0;
	mfem::VectorConstantCoefficient constant( componentValue );
	q.ProjectCoefficient( constant );

	mfem::GridFunction field;
	meq::poloidalField( q, field );

	// B_R = -q_z = +7, B_Z = +q_r = +3.
	mfem::Vector sampled( 2 );
	mfem::IntegrationPoint centre;
	centre.Set2( 1.0/3.0, 1.0/3.0 );
	field.GetVectorValue( 0, centre, sampled );

	BOOST_TEST( sampled( 0 ) == 7.0, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( sampled( 1 ) == 3.0, boost::test_tools::tolerance( 1.0e-12 ) );
}

/*
 * AND THE FIELD OF AN ACTUAL SOLVE, against the exact one two ways.
 *
 * The error must be the DISCRETISATION error of q -- there is no extra
 * approximation in the relabelling -- so it converges at k+1 and its size is
 * that of solver.fluxError(). A sign error would show as an error of the size of
 * B itself, which the printed magnitude column is there to compare against.
 */
BOOST_AUTO_TEST_CASE( theFieldMatchesTheExactSolutionAndConvergesAtKPlusOne )
{
	for ( int order = 1; order <= 3; ++order )
	{
		std::printf( "\n  B_poloidal, Solov'ev NSTX, k = %d\n", order );
		std::printf( "  %8s %14s %14s %12s %7s\n", "h", "vs analytic",
		             "vs differences", "|B|", "rate" );

		std::vector<FieldMeasurement> points;
		for ( int n : { 8, 16, 32 } )
			points.push_back( measure( order, n ) );

		for ( std::size_t i = 0; i < points.size(); ++i )
		{
			FieldMeasurement const &p = points[ i ];
			double rate = 0.0;
			if ( i > 0 )
				rate = std::log( points[ i - 1 ].worstAgainstAnalytic
				                 /p.worstAgainstAnalytic )/std::log( 2.0 );
			std::printf( "  %8.5f %14.6e %14.6e %12.4e %7s",
			             p.h, p.worstAgainstAnalytic, p.worstAgainstDifference,
			             p.magnitude, i > 0 ? "" : "-" );
			if ( i > 0 ) std::printf( "%7.3f", rate );
			std::printf( "\n" );
		}
		std::fflush( stdout );

		FieldMeasurement const &finest = points.back();

		// The two references must agree with each other to the finite-difference
		// floor: if they do not, the fixture's gradPsi() and its psi() disagree
		// and this file is measuring the wrong thing.
		BOOST_TEST( std::abs( finest.worstAgainstAnalytic
		                      - finest.worstAgainstDifference )
		            < 1.0e-6*std::max( 1.0, finest.magnitude ),
		            "the analytic flux and the central differences of psi "
		            "disagree, so Soloviev.hpp is inconsistent with itself" );

		// A WRONG SIGN would leave an error of order |B|. This is the assertion
		// that catches it, and it is deliberately crude -- it does not need to
		// be tight to be decisive.
		BOOST_TEST( finest.worstAgainstAnalytic < 0.05*finest.magnitude,
		            "k = " << order << ": the field is wrong by "
		            << finest.worstAgainstAnalytic << " against a field of size "
		            << finest.magnitude << ", which is the scale of a flipped "
		            "sign or a swapped component rather than of a discretisation "
		            "error" );

		// And the rate, so a field that is merely inaccurate is caught too.
		//
		// Slack of 0.5, wider than the 0.15 the L2 studies use, and deliberately:
		// this is a MAX norm over a cloud of FIXED points on a discontinuous
		// field. Which element a point lands in changes with the mesh, so the
		// measured maximum jumps about -- the per-pair rates printed above run
		// 3.215 then 0.827 at k = 1 on a sequence whose overall rate is 2.02.
		// Measured overall: 2.02, 2.72, 4.29 for k = 1, 2, 3, so 2.72 against a
		// 2.65 threshold would have been a test that flakes rather than one that
		// measures.
		//
		// The rate is the secondary assertion here in any case. The sign check
		// above is what this file exists for and it has four orders of margin.
		double const rate = std::log( points.front().worstAgainstAnalytic
		                              /points.back().worstAgainstAnalytic )
		                    /std::log( 4.0 );
		BOOST_TEST( rate > order + 1.0 - 0.5,
		            "k = " << order << ": B converged at " << rate
		            << ", short of the k+1 that q converges at -- the relabelling "
		            "should add no error of its own" );
	}
}

BOOST_AUTO_TEST_SUITE_END()
