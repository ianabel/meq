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
 * only a timing catches -- so it is timed rather than asserted to be linear in
 * a comment, because complexity claims rot.
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
		box().minRadius + inset, box().maxRadius - inset, 41,
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

/**
 * The affine inverse must locate exactly what the Newton inverse locates.
 *
 * A STRAIGHT-SIDED TRIANGLE'S MAP IS AFFINE, so meq::GridSampler inverts it
 * with a 2x2 solve rather than through ElementTransformation::TransformBack(),
 * which constructs an InverseElementTransformation and iterates. That is worth
 * about seventy times the cost of the constructor and it is only allowed to be
 * if it decides the SAME THING -- both which element owns each node and where
 * in it, since a node on an inter-element face is a legitimate candidate for
 * two elements and the class's contract is that the first one to reach it keeps
 * it.
 *
 * SO THE CONTROL IS THE NEWTON ROUTE ITSELF, run here over every element in
 * index order, which is what the constructor would have done. The grid is
 * chosen so that every fifth line falls exactly on a mesh line -- 41 nodes
 * across 8 cells -- because a node in the middle of an element is the case that
 * cannot distinguish the two and a node on a face is the case that can.
 *
 * THE ELEMENT IS READ OUT THROUGH A PIECEWISE-CONSTANT FIELD WHOSE VALUE IS THE
 * ELEMENT INDEX, which is exact in P_0 and so reports ownership rather than
 * approximating it. A smooth field is sampled beside it to pin the reference
 * POINT as well: the same element with a different point would pass the first
 * assertion and fail the second.
 */
BOOST_AUTO_TEST_CASE( theAffineInverseLocatesWhatTheNewtonInverseDoes )
{
	mfem::Mesh mesh = meq::tests::makeMesh( box(), 8 );

	// The precondition, asserted rather than assumed: a mesh carrying nodes
	// would take the Newton route in both arms and the case would pass
	// vacuously.
	BOOST_TEST_REQUIRE( mesh.GetNodes() == nullptr,
	                    "this mesh is curved, so the affine path is not the one "
	                    "under test" );

	int const nodes = 41;
	meq::GridSampler sampler( mesh,
		box().minRadius, box().maxRadius, nodes,
		box().zMin, box().zMax, nodes );

	// P_0: one dof per element, set to the element's own index.
	mfem::L2_FECollection constants( 0, mesh.Dimension() );
	mfem::FiniteElementSpace constantSpace( &mesh, &constants );
	mfem::GridFunction owner( &constantSpace );
	{
		mfem::Array<int> dofs;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			constantSpace.GetElementDofs( e, dofs );
			owner( dofs[ 0 ] ) = e;
		}
	}

	// Something with structure in it, so that a wrong reference point in the
	// right element shows up.
	mfem::L2_FECollection smoothColl( 3, mesh.Dimension() );
	mfem::FiniteElementSpace smoothSpace( &mesh, &smoothColl );
	mfem::GridFunction smooth( &smoothSpace );
	mfem::FunctionCoefficient smoothCoeff( []( mfem::Vector const &x )
	{
		return std::sin( 7.0*x( 0 ) )*std::cos( 5.0*x( 1 ) ) + 2.0*x( 0 )*x( 1 );
	} );
	smooth.ProjectCoefficient( smoothCoeff );

	std::vector<double> sampledOwner, sampledSmooth;
	sampler.sample( owner, sampledOwner, -1.0 );
	sampler.sample( smooth, sampledSmooth, 0.0 );

	// The control: the element loop the constructor would have run, with the
	// Newton inverse deciding.
	int checked = 0, onAFace = 0;
	double worstValue = 0.0;
	mfem::Vector physical( 2 );
	for ( int j = 0; j < nodes; ++j )
		for ( int i = 0; i < nodes; ++i )
		{
			physical( 0 ) = sampler.rAt( i );
			physical( 1 ) = sampler.zAt( j );

			int expected = -1;
			int claims = 0;
			mfem::IntegrationPoint reference;
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				mfem::IsoparametricTransformation transformation;
				mesh.GetElementTransformation( e, &transformation );
				mfem::IntegrationPoint here;
				if ( transformation.TransformBack( physical, here )
				     == mfem::InverseElementTransformation::Inside )
				{
					++claims;
					if ( expected < 0 )
					{
						expected = e;
						reference = here;
					}
				}
			}
			if ( claims > 1 )
				++onAFace;

			std::size_t const at = static_cast<std::size_t>( j )*nodes + i;
			BOOST_TEST_REQUIRE( sampler.located( i, j ) == ( expected >= 0 ),
			                    "node ( " << i << ", " << j << " ) is located by "
			                    "one inverse and not the other" );
			if ( expected < 0 )
				continue;

			++checked;
			BOOST_TEST_REQUIRE( static_cast<int>( sampledOwner[ at ] + 0.5 ) == expected,
			                    "node ( " << i << ", " << j << " ) was filed under "
			                    "element " << static_cast<int>( sampledOwner[ at ] + 0.5 )
			                    << " and the Newton inverse gives " << expected );

			worstValue = std::max( worstValue,
				std::abs( sampledSmooth[ at ]
				          - smooth.GetValue( expected, reference ) ) );
		}

	// A grid line on a mesh line is the whole point of the geometry above, so
	// if none of them landed there the case is weaker than it reads.
	BOOST_TEST( onAFace > 0,
	            "no grid node fell on an inter-element face, so this case never "
	            "exercised the ownership rule it exists for" );

	std::printf( "\n  affine vs Newton: %d nodes located identically, %d of them "
	             "claimed by more than one element, worst value difference %.3e\n",
	             checked, onAFace, worstValue );
	std::fflush( stdout );

	BOOST_TEST( worstValue < 1.0e-12,
	            "the two inverses agree on the element and not on the point in "
	            "it: the sampled field differs by " << worstValue );
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
		box().minRadius - padR, box().maxRadius + padR, 61,
		box().zMin - padZ, box().zMax + padZ, 61 );

	int inside = 0, outside = 0;
	for ( int j = 0; j < sampler.nodesZ(); ++j )
		for ( int i = 0; i < sampler.nodesR(); ++i )
		{
			bool const within = sampler.rAt( i ) > box().minRadius + 1.0e-9
			                 && sampler.rAt( i ) < box().maxRadius - 1.0e-9
			                 && sampler.zAt( j ) > box().zMin + 1.0e-9
			                 && sampler.zAt( j ) < box().zMax - 1.0e-9;
			if ( within )
			{
				++inside;
				BOOST_TEST( sampler.located( i, j ),
				            "a node strictly inside the mesh was not located" );
			}
			else if ( sampler.rAt( i ) < box().minRadius - 1.0e-9
			       || sampler.rAt( i ) > box().maxRadius + 1.0e-9
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
		box().minRadius + inset, box().maxRadius - inset, 65,
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

			double const radius = sampler.rAt( i ), z = sampler.zAt( j );
			worstPsi = std::max( worstPsi, std::abs( psi[ at ] - eq.psi( radius, z ) ) );

			double qR = 0.0, qZ = 0.0;
			eq.flux( radius, z, qR, qZ );
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
			box().minRadius + 0.01, box().maxRadius - 0.01, nodes,
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

/*
 * THE TWO LOCATORS, ASKED THE SAME QUESTION AT THE SAME POINTS.
 *
 * MEQ finds points in a mesh by two unrelated routes and neither had ever been
 * checked against the other. meq::GridSampler inverts the loop -- element
 * bounding box to grid index range, then an affine or Newton TransformBack --
 * and meq::FieldTransfer goes through mfem::FindPointsGSLIB, whose simplex path
 * splits triangles into quads internally. MEQ's meshes are triangles, so that
 * is a code path a quad-based user never exercises. WarmStartConvergence does
 * run it on triangles, but against the EXACT SOLUTION, which cannot separate
 * "both right" from "both wrong the same way".
 *
 * IT IS gslib THAT IS COMPARED HERE AND NOT meq::FieldTransfer, and the
 * distinction is why this case is shaped as it is. FieldTransfer::transfer is
 * an element-local L2 PROJECTION onto the target space, sampling at the
 * target's Gauss points and never at its dofs -- a deliberate choice its own
 * comment defends. So a transfer's output differs from a point sample by its
 * own projection error, which is a property of the target space and says
 * nothing about locating. Comparing the LOCATORS means comparing what each
 * READS at one set of points.
 *
 * THE NODES SPLIT IN TWO AND ONLY ONE HALF CARRIES A CLAIM. psi_h is L2 and so
 * discontinuous across an element face: a node ON a face has two right answers
 * and the two routes need not choose the same one. A node strictly inside an
 * element has exactly one, and there the two must agree to round-off or one of
 * them is in the wrong element -- the failure no comparison against an exact
 * solution can see, because a smooth field is nearly right in the neighbour
 * too.
 *
 * gslib's own code 1 is what classifies a node as on a border, because
 * GridSampler has no notion of one -- it picks whichever element claims the
 * node first and says nothing about the choice. That asymmetry is what TODO's
 * first sampling question proposed to close with SetL2AvgType, and the
 * measurement below closes the question the other way: see the comment on the
 * NONE setting a few lines down. What GridSampler is missing is a way to
 * REPORT that a node was ambiguous, not a different value to return for it.
 */
BOOST_AUTO_TEST_CASE( theTwoLocatorsAgreeAtTheSamePoints )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	int const n = 32;
	mfem::Mesh mesh = meq::tests::makeMesh( box(), n );

	mfem::FunctionCoefficient source( [ &eq ]( mfem::Vector const &x )
		{ return eq.f( x( 0 ), x( 1 ), 0.0 ); } );
	mfem::FunctionCoefficient exact( [ &eq ]( mfem::Vector const &x )
		{ return eq.psi( x( 0 ), x( 1 ) ); } );

	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( source );
	solver.setBoundaryData( exact );
	solver.solve();

	/*
	 * gslib REQUIRES A NODAL MESH and Setup() installs one -- the same call
	 * meq::FieldTransfer's constructor makes, for the same reason. That matters
	 * here beyond satisfying gslib: GridSampler's affine fast path asks whether
	 * the GEOMETRY is affine, and a mesh that has just acquired an ORDER-1
	 * nodal field still is. Every sampler below is therefore built AFTER this,
	 * deliberately, so the comparison runs on a mesh in the state a warm start
	 * leaves it in -- which is the state in which testing GetNodes() for null
	 * would silently cost two orders of magnitude.
	 */
	mesh.EnsureNodes();

	mfem::FindPointsGSLIB finder;
	finder.Setup( mesh );

	/*
	 * NONE, AND IT IS NOT THE DEFAULT, WHICH IS THE FIRST THING THIS CASE
	 * FOUND. FindPointsGSLIB::Interpolate() on an L2 field goes back over every
	 * point it classified as ON A BORDER, projects the whole field to H1 and
	 * re-interpolates those -- so with the default averaging left in place the
	 * value gslib returns at a border node belongs to NEITHER element's
	 * polynomial and the comparison below stops being about locating at all.
	 * meq::FieldTransfer sets NONE in its own constructor for the same reason
	 * and records what the default costs it: 28%, 10% and 12% in L2 at
	 * k = 1, 2, 3 on a transfer whose exact answer exists.
	 *
	 * Measured here, leaving it at the default makes EVERY border node differ
	 * -- 577 of 577 on the generic grid, worst 1.93e-06 -- because the value
	 * returned is no longer any element's polynomial. With NONE only the 60
	 * that the two routes actually resolve to different elements differ, and
	 * the other 517 agree to round-off. So the default does not measure
	 * locating at all.
	 *
	 * AND THAT SETTLES TODO's FIRST SAMPLING QUESTION RATHER THAN LEAVING IT
	 * OPEN. The suggestion there was that SetL2AvgType is "a principled answer
	 * where meq::GridSampler has a shrug". The rule exists and is the WRONG
	 * one for an L2 field: it returns a smoothed value belonging to neither
	 * side of the jump, which is why meq::FieldTransfer measured it 28%, 10%
	 * and 12% wrong and turned it off. Picking a side -- which is what
	 * GridSampler already does -- is the better rule, and the shrug was right.
	 */
	finder.SetL2AvgType( mfem::FindPointsGSLIB::NONE );

	/// A P0 field whose value IS the element index, so that "which element did
	/// you choose" is reported exactly rather than inferred from a value.
	mfem::L2_FECollection p0( 0, 2 );
	mfem::FiniteElementSpace ownerSpace( &mesh, &p0 );
	mfem::GridFunction owner( &ownerSpace );
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		mfem::Array<int> dofs;
		ownerSpace.GetElementDofs( e, dofs );
		owner( dofs[ 0 ] ) = e;
	}

	struct Outcome
	{
		int interior = 0, border = 0, misses = 0;
		int interiorDiffering = 0, borderDiffering = 0, differentElement = 0;
		double worstInterior = 0.0, worstBorder = 0.0;
	};

	auto compare = [ & ]( char const *what,
	                      double minRadius, double maxRadius, int nR,
	                      double zMin, double zMax, int nZ )
	{
		meq::GridSampler sampler( mesh, minRadius, maxRadius, nR, zMin, zMax, nZ );

		std::vector<double> grid, mine;
		sampler.sample( solver.potential(), grid, std::nan( "" ) );
		sampler.sample( owner, mine, -1.0 );

		// The SAME points, in gslib's byNODES ordering: every R, then every z.
		// Passed explicitly because the two orderings differ silently and a
		// transposed cloud gives a plausible wrong answer rather than an error.
		int const total = nR*nZ;
		mfem::Vector positions( 2*total );
		for ( int j = 0; j < nZ; ++j )
			for ( int i = 0; i < nR; ++i )
			{
				positions( j*nR + i ) = sampler.rAt( i );
				positions( total + j*nR + i ) = sampler.zAt( j );
			}

		finder.FindPoints( positions, mfem::Ordering::byNODES );

		mfem::Vector values;
		finder.Interpolate( solver.potential(), values );
		mfem::Array<unsigned int> const &code = finder.GetCode();
		mfem::Array<unsigned int> const &element = finder.GetElem();

		Outcome out;
		for ( int at = 0; at < total; ++at )
		{
			int const i = at%nR, j = at/nR;
			if ( code[ at ] == 2 )
			{
				++out.misses;
				continue;
			}
			if ( !sampler.located( i, j ) )
				continue;

			double const difference =
				std::abs( grid[ static_cast<std::size_t>( at ) ] - values( at ) );
			// Round-off, not a tolerance on physics: two evaluations of one
			// polynomial at one point in one element.
			bool const differs = difference > 1.0e-12;

			if ( code[ at ] == 1 )
			{
				++out.border;
				out.worstBorder = std::max( out.worstBorder, difference );
				if ( differs )
					++out.borderDiffering;
				if ( static_cast<int>( element[ at ] )
				     != static_cast<int>( std::lround( mine[ at ] ) ) )
					++out.differentElement;
			}
			else
			{
				++out.interior;
				out.worstInterior = std::max( out.worstInterior, difference );
				if ( differs )
					++out.interiorDiffering;
			}
		}

		std::printf( "  %-14s interior %5d (%d differing, worst %.3e)   "
		             "border %5d (%d differing, %d in a different element, "
		             "worst %.3e)   misses %d\n",
		             what, out.interior, out.interiorDiffering, out.worstInterior,
		             out.border, out.borderDiffering, out.differentElement,
		             out.worstBorder, out.misses );
		std::fflush( stdout );
		return out;
	};

	std::printf( "\n  GridSampler against FindPointsGSLIB, psi_h at k = 2 on %d triangles\n",
	             mesh.GetNE() );

	double const inset = 0.05;
	Outcome const generic = compare( "generic:", box().minRadius + inset, box().maxRadius - inset, 65,
	                                 box().zMin + inset, box().zMax - inset, 65 );
	// Every node on a mesh line: n + 1 nodes across the full box is the mesh's
	// own spacing, and the triangulation's diagonals cross those nodes too.
	Outcome const aligned = compare( "face-aligned:", box().minRadius, box().maxRadius, n + 1,
	                                 box().zMin, box().zMax, n + 1 );

	BOOST_TEST_REQUIRE( generic.interior > 3000 );

	/*
	 * THE CLAIM. Where the point is unambiguous the two locators are reading one
	 * polynomial on one element, and they share no locating machinery whatever,
	 * so they must agree to round-off. Measured: 3648 of 3648 interior nodes,
	 * worst difference below 1e-12.
	 */
	BOOST_TEST( generic.interiorDiffering == 0,
	            generic.interiorDiffering << " of " << generic.interior
	            << " nodes strictly inside an element disagree, worst "
	            << generic.worstInterior << " -- one locator is placing points in "
	            "the wrong element" );
	BOOST_TEST( aligned.interiorDiffering == 0 );

	/*
	 * AND ON A BORDER BOTH ARE RIGHT, so what is asserted there is a BOUND. The
	 * difference is the face jump of an L2 field, O( h^(k+1) ), which is bounded
	 * by the solution's own error up to a constant; far above it would mean a
	 * non-adjacent element rather than the ambiguity.
	 */
	double worstAgainstExact = 0.0;
	{
		meq::GridSampler check( mesh, box().minRadius, box().maxRadius, n + 1,
		                              box().zMin, box().zMax, n + 1 );
		std::vector<double> psi;
		check.sample( solver.potential(), psi, std::nan( "" ) );
		for ( int j = 0; j < check.nodesZ(); ++j )
			for ( int i = 0; i < check.nodesR(); ++i )
				if ( check.located( i, j ) )
					worstAgainstExact = std::max( worstAgainstExact,
						std::abs( psi[ static_cast<std::size_t>( j )*check.nodesR() + i ]
						          - eq.psi( check.rAt( i ), check.zAt( j ) ) ) );
	}

	std::printf( "  worst border difference %.3e generic, %.3e aligned, against "
	             "%.3e of error in psi_h itself\n",
	             generic.worstBorder, aligned.worstBorder, worstAgainstExact );
	std::fflush( stdout );

	BOOST_TEST( aligned.worstBorder < 20.0*worstAgainstExact,
	            "the two locators differ by " << aligned.worstBorder << " on border "
	            "nodes against " << worstAgainstExact << " of error in psi_h. A face "
	            "jump is bounded by the solution error; this much larger is a "
	            "non-adjacent element, not the ambiguity" );

	/*
	 * AND THE AMBIGUITY IS COMMON RATHER THAN A CURIOSITY. On a grid with no
	 * alignment intended, 577 of 4225 nodes -- ONE IN SEVEN -- land on a
	 * border, because a uniform grid over a uniform mesh is commensurate with
	 * it far more often than a generic point would be. On the aligned grid it
	 * is all of them, and there the two routes agree on every one.
	 *
	 * THE DIFFERENT-ELEMENT COUNT IS PRINTED AND NEVER ASSERTED, because it is
	 * not a stable quantity: it is decided by which side of a face a node falls
	 * on at the last bit. Measured, building the same mesh with an extent of
	 * 0.8 rather than 1.4 - 0.6 -- one ulp, 0.7999999999999999 against 0.8 --
	 * moves it from 202 to 331 on an otherwise identical grid. A gate on it
	 * would be a gate on round-off in the mesh generator.
	 */
	BOOST_TEST( generic.border > 0,
	            "no node landed on a border, so this grid exercises none of the "
	            "ambiguity and the bound above is vacuous" );
	std::printf( "  border nodes: %d of %d generic (%.1f%%), %d of %d aligned\n\n",
	             generic.border, generic.border + generic.interior,
	             100.0*generic.border/( generic.border + generic.interior ),
	             aligned.border, aligned.border + aligned.interior );
	std::fflush( stdout );
}

BOOST_AUTO_TEST_SUITE_END()
