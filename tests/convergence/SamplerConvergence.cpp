#define BOOST_TEST_MODULE MeqSamplerConvergence

#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/Field.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Sampler.hpp"

#include "ConvergenceHarness.hpp"
#include "analytic/Soloviev.hpp"

/*
 * The ( R, Z ) sampler: does it find the right value, and does it scale?
 *
 * Two claims to check and they fail differently. A sampler that locates points
 * in the wrong element gives wrong numbers, which an accuracy test catches. A
 * sampler that locates them by brute force gives the RIGHT numbers slowly, which
 * only a timing catches -- and DRIVER-PLAN section 3 says to time it rather than
 * assert linearity in a comment, because complexity claims rot.
 */

namespace
{
	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq =
			meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	meq::tests::Rectangle box()
	{
		return meq::tests::Rectangle{ 0.6, 1.4, -0.6, 0.6 };
	}
}

BOOST_AUTO_TEST_SUITE( sampler_convergence )

/// Every interior node of a grid over the mesh must be found, and found in an
/// element that really contains it. Checked by transforming the reference point
/// forward again and comparing with the node it came from -- which catches a
/// point recorded against the wrong element, the failure an accuracy test on a
/// smooth field would hide.
BOOST_AUTO_TEST_CASE( everyNodeIsFoundInAnElementThatContainsIt )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 8 );

	// Inset by a hair so that nodes on the boundary are unambiguous.
	double const inset = 1.0e-9;
	meq::GridSampler sampler( mesh,
		box().rMin + inset, box().rMax - inset, 41,
		box().zMin + inset, box().zMax - inset, 41 );

	BOOST_TEST( sampler.locatedCount() == 41*41,
	            "only " << sampler.locatedCount() << " of " << 41*41
	            << " nodes were located, on a grid entirely inside the mesh" );

	// The round trip: sampling the identity coefficients gives back the physical
	// coordinates of whatever element and reference point each node was recorded
	// against. If a node was filed under the wrong element, they will not match.
	mfem::FunctionCoefficient rCoefficient(
		[]( mfem::Vector const &x ) { return x( 0 ); } );
	mfem::FunctionCoefficient zCoefficient(
		[]( mfem::Vector const &x ) { return x( 1 ); } );

	std::vector<double> rs, zs;
	sampler.sampleCoefficient( rCoefficient, rs, 0.0 );
	sampler.sampleCoefficient( zCoefficient, zs, 0.0 );

	double worst = 0.0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			BOOST_TEST_REQUIRE( sampler.located( i, j ) );
			std::size_t const at = static_cast<std::size_t>( j )*sampler.nodesR() + i;
			worst = std::max( worst, std::hypot( rs[ at ] - sampler.rAt( i ),
			                                     zs[ at ] - sampler.zAt( j ) ) );
		}

	BOOST_TEST( worst < 1.0e-10,
	            "a node was located in an element that does not contain it: the "
	            "recorded reference point maps back " << worst << " away" );
}

/// Outside is not an error, it is the mask. A grid larger than the mesh must
/// find the interior and miss the exterior, which is what the output uses to
/// decide where psi is defined.
BOOST_AUTO_TEST_CASE( nodesOutsideTheMeshAreNotLocated )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 8 );

	// Half again as wide as the mesh in each direction, centred on it.
	double const padR = 0.5*box().width();
	double const padZ = 0.5*box().height();
	meq::GridSampler sampler( mesh,
		box().rMin - padR, box().rMax + padR, 61,
		box().zMin - padZ, box().zMax + padZ, 61 );

	int inside = 0, outside = 0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			bool const within = sampler.rAt( i ) > box().rMin + 1.0e-9
			                 && sampler.rAt( i ) < box().rMax - 1.0e-9
			                 && sampler.zAt( j ) > box().zMin + 1.0e-9
			                 && sampler.zAt( j ) < box().zMax - 1.0e-9;
			if ( within )
			{
				++inside;
				BOOST_TEST( sampler.located( i, j ),
				            "a node strictly inside the mesh was not located" );
			}
			else if ( sampler.rAt( i ) < box().rMin - 1.0e-9
			       || sampler.rAt( i ) > box().rMax + 1.0e-9
			       || sampler.zAt( j ) < box().zMin - 1.0e-9
			       || sampler.zAt( j ) > box().zMax + 1.0e-9 )
			{
				++outside;
				BOOST_TEST( !sampler.located( i, j ),
				            "a node outside the mesh was located anyway" );
			}
		}

	BOOST_TEST( inside > 0 );
	BOOST_TEST( outside > 0 );
	std::printf( "\n  mask: %d nodes inside, %d outside, %d located of %d\n",
	             inside, outside, sampler.locatedCount(),
	             sampler.nodesR()*sampler.nodesZ() );
	std::fflush( stdout );
}

/// The sampled field must be the field. psi and both components of B, against
/// the exact solution, at the discretisation error -- so the sampler adds
/// nothing of its own.
BOOST_AUTO_TEST_CASE( theSampledFieldsMatchTheExactSolution )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 32 );

	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{ return eq.f( x( 0 ), x( 1 ), 0.0 ); } );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{ return eq.psi( x( 0 ), x( 1 ) ); } );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( exact );
	solver.solve();

	mfem::GridFunction field;
	meq::poloidalField( solver.flux(), field );

	double const inset = 0.05;
	meq::GridSampler sampler( mesh,
		box().rMin + inset, box().rMax - inset, 65,
		box().zMin + inset, box().zMax - inset, 65 );

	std::vector<double> psi, bR, bZ;
	sampler.sample( solver.potential(), psi, std::nan( "" ) );
	sampler.sampleComponent( field, 0, bR, std::nan( "" ) );
	sampler.sampleComponent( field, 1, bZ, std::nan( "" ) );

	double worstPsi = 0.0, worstB = 0.0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			std::size_t const at = static_cast<std::size_t>( j )*sampler.nodesR() + i;
			BOOST_TEST_REQUIRE( sampler.located( i, j ) );

			double const r = sampler.rAt( i ), z = sampler.zAt( j );
			worstPsi = std::max( worstPsi, std::abs( psi[ at ] - eq.psi( r, z ) ) );

			double qR = 0.0, qZ = 0.0;
			eq.flux( r, z, qR, qZ );
			worstB = std::max( worstB, std::hypot( bR[ at ] + qZ, bZ[ at ] - qR ) );
		}

	std::printf( "  sampled on 65x65 at k = 2, h = %.4f: worst psi %.4e, worst B %.4e\n",
	             box().width()/32.0, worstPsi, worstB );
	std::fflush( stdout );

	// Ceilings set from the measurement at about three times it, as the
	// exact-solution studies in this directory set theirs. These are MAX norms
	// at sample points, so they sit an order above the L2 errors the rate tables
	// record -- 1.9e-7 for psi at k = 2 on this mesh. Measured: psi 1.88e-6,
	// B 1.19e-6.
	BOOST_TEST( worstPsi < 6.0e-6 );
	BOOST_TEST( worstB < 4.0e-6 );
}

/*
 * AND THE COMPLEXITY CLAIM, timed rather than asserted in a comment.
 *
 * The inverted loop is O( elements x points per element ). Doubling the grid in
 * each direction quadruples the node count and should roughly quadruple the
 * cost; a brute-force locator would be quadratic in the same sweep. The bound
 * below is deliberately loose -- this measures wall clock on a shared machine --
 * but it separates linear from quadratic by a wide margin, which is all it needs
 * to do.
 */
BOOST_AUTO_TEST_CASE( locatingIsLinearInTheNodeCount )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 32 );

	auto time = [ &mesh ]( int nodes )
	{
		auto const start = std::chrono::steady_clock::now();
		meq::GridSampler sampler( mesh,
			box().rMin + 0.01, box().rMax - 0.01, nodes,
			box().zMin + 0.01, box().zMax - 0.01, nodes );
		auto const stop = std::chrono::steady_clock::now();
		BOOST_TEST( sampler.locatedCount() > 0 );
		return std::chrono::duration<double, std::milli>( stop - start ).count();
	};

	double const small = time( 65 );
	double const large = time( 129 );      // four times the nodes

	std::printf( "  locating 65x65 %.1f ms, 129x129 %.1f ms, ratio %.2f "
	             "(4 is linear, 16 is quadratic)\n",
	             small, large, large/std::max( 1.0e-6, small ) );
	std::fflush( stdout );

	BOOST_TEST( large/std::max( 1.0e-6, small ) < 8.0,
	            "locating cost grew by " << large/std::max( 1.0e-6, small )
	            << " for four times the nodes, which is closer to quadratic than "
	            "to linear -- the element loop is not being inverted" );
}

BOOST_AUTO_TEST_SUITE_END()
