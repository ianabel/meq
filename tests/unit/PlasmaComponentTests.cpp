// Unit tests for the plasma connectivity fill in src/meq/PlasmaComponent.
//
// FREE-BOUNDARY-PLAN.md stage XP-1. meq::NormalisedSource::insidePlasma() is a
// POINTWISE test on the value, so `{ Psi > 0 }` is whatever the level set
// happens to be and not necessarily the plasma. This is the connectivity the
// pointwise test is missing, as a graph algorithm on plain data.
//
// MFEM-FREE, like the unit under test, so CI can run it -- and the
// configurations that matter are far easier to write down as graphs than as
// meshes:
//
//   * TWO LOBES TOUCHING AT A VERTEX. `../freegs4e` has to block a
//     neighbourhood of each X-point before its grid fill, because on a uniform
//     grid two lobes meeting at a point are eight-neighbours. Over FACE
//     neighbours they are not adjacent at all, and section 10.3 predicts MEQ
//     therefore needs no blocking. Half of that is right and this file pins the
//     half that is: a vertex-touching pair really is separated.
//   * TWO LOBES BRIDGED BY ONE ELEMENT, which is the other half and is what the
//     X-point's own element does -- it is cut by both branches and is a face
//     neighbour of both lobes. The fill goes straight through it, and once it
//     is on the far side every further element is an ordinary face neighbour.
//     So the leak is a REGION and not one element, and blocking that one
//     element is what stops it.
//   * A SEED IN THE WRONG LOBE, because the fill names whichever component the
//     seed is in and nothing else can tell it which one is the plasma.
//
// tests/convergence/PlasmaConnectivity.cpp is the same three questions on a
// real mesh with a real field, where the graph comes from
// mfem::Mesh::ElementToElementTable() and the X-point is Soloviev::nstx()'s
// closed-form one.

#define BOOST_TEST_MODULE PlasmaComponentTests
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "meq/PlasmaComponent.hpp"

using meq::PlasmaComponent;

namespace
{
	/**
	 * A 4-connected `rows` x `columns` grid of quadrilaterals, as CSR.
	 *
	 * FACE neighbours only -- left, right, up, down -- which is the whole point
	 * of the comparison with a grid code: the diagonal neighbours a uniform-grid
	 * flood fill would take are deliberately absent, because a mesh element
	 * shares a FACE with four of its neighbours and only a VERTEX with the
	 * corner four.
	 */
	PlasmaComponent grid( int rows, int columns )
	{
		std::vector< int > offsets{ 0 };
		std::vector< int > neighbours;

		auto index = [ & ]( int i, int j ) { return i*columns + j; };

		for ( int i = 0; i < rows; ++i )
			for ( int j = 0; j < columns; ++j )
			{
				if ( i > 0 )            neighbours.push_back( index( i - 1, j ) );
				if ( i + 1 < rows )     neighbours.push_back( index( i + 1, j ) );
				if ( j > 0 )            neighbours.push_back( index( i, j - 1 ) );
				if ( j + 1 < columns )  neighbours.push_back( index( i, j + 1 ) );
				offsets.push_back( static_cast< int >( neighbours.size() ) );
			}

		PlasmaComponent component;
		component.setAdjacency( std::move( offsets ), std::move( neighbours ) );
		return component;
	}
}


/*
 * THE STATE BEFORE fill() IS THE CONSTANT TRUE, AND THAT IS NOT A DETAIL.
 *
 * Every solve in this tree that does not ask for a connectivity test must be
 * BIT-unchanged by this class existing, and the way that is guaranteed is that
 * an unfilled PlasmaComponent answers true everywhere -- so the caller tests
 * unconditionally and there is no second code path to keep in step. A version
 * that answered false, or that the caller had to remember to skip, would be the
 * same shape of trap as a mask nobody consults.
 */
BOOST_AUTO_TEST_CASE( an_unfilled_component_holds_everywhere )
{
	PlasmaComponent empty;
	BOOST_TEST( !empty.active() );
	BOOST_TEST( empty.holds( 0 ) );
	BOOST_TEST( empty.holds( 12345 ) );
	BOOST_TEST( empty.holds( -1 ) );

	PlasmaComponent g = grid( 3, 3 );
	BOOST_TEST( g.nodeCount() == 9 );
	BOOST_TEST( !g.active() );
	for ( int e = 0; e < 9; ++e )
		BOOST_TEST( g.holds( e ) );

	// And clear() puts it back, because the solver clears the mask whenever the
	// support is switched off and the next solve must not inherit one.
	g.fill( std::vector< char >( 9, 1 ), 0 );
	BOOST_TEST( g.active() );
	g.clear();
	BOOST_TEST( !g.active() );
	BOOST_TEST( g.holds( 4 ) );
}


/*
 * TWO LOBES TOUCHING AT A VERTEX ARE SEPARATE, WHICH IS THE PREDICTION.
 *
 * Section 10.3: "A fill over FACE neighbours cannot cross a shared vertex, so
 * it should be blocked at the saddle automatically". Here is that statement in
 * its purest form -- a 3x3 grid whose two diagonal corners carry plasma and
 * whose centre does not:
 *
 *     X . .          a grid code's eight-neighbour fill joins these two,
 *     . . .          and has to block the centre to stop it. A face
 *     . . X          neighbourhood never had them adjacent.
 *
 * The CONTROL is the same two lobes with the centre carrying plasma too, which
 * is the bridged case below and reads one component -- so this is not passing
 * because the fill fails to move at all.
 */
BOOST_AUTO_TEST_CASE( two_lobes_meeting_at_a_vertex_are_not_connected )
{
	PlasmaComponent g = grid( 3, 3 );

	std::vector< char > mask( 9, 0 );
	mask[ 0 ] = 1;   // ( 0, 0 )
	mask[ 8 ] = 1;   // ( 2, 2 ), diagonally opposite and vertex-adjacent to it
	                 // only through ( 1, 1 ), which is empty.

	g.fill( mask, 0 );

	BOOST_TEST( g.componentCount() == 2,
	            "two diagonally opposite cells came back as "
	            << g.componentCount() << " component(s): a face-neighbour fill "
	            "must not cross a shared vertex, which is the property that "
	            "lets MEQ skip the X-point blocking freegs4e needs" );
	BOOST_TEST( g.candidateNodes() == 2 );
	BOOST_TEST( g.componentNodes() == 1 );
	BOOST_TEST( g.holds( 0 ) );
	BOOST_TEST( !g.holds( 8 ) );
}


/*
 * AND ONE ELEMENT IS ENOUGH TO JOIN THEM, WHICH IS THE OTHER HALF.
 *
 * The element that CONTAINS an X-point is cut by both branches of the
 * separatrix, so it carries Psi > 0 and is a face neighbour of both lobes.
 * Section 10.3 predicts the resulting leak is "one element wide"; it is not.
 * Once the fill is through the bridge, every further cell of the far lobe is an
 * ordinary face neighbour and the whole region is taken -- which this case
 * measures directly by making the far lobe large.
 *
 * blockNode is the cure and it is ONE node, which is what the face-neighbour
 * graph buys over a grid: freegs4e blocks a NEIGHBOURHOOD of the saddle.
 */
BOOST_AUTO_TEST_CASE( one_bridging_element_leaks_the_whole_far_lobe )
{
	// A 3 x 9 grid. The middle row carries plasma from end to end; the bridge
	// is the single cell ( 1, 4 ) and the far lobe is everything past it.
	int const rows = 3;
	int const columns = 9;
	PlasmaComponent g = grid( rows, columns );

	std::vector< char > mask( static_cast< std::size_t >( rows*columns ), 0 );
	auto at = [ & ]( int i, int j ) { return static_cast< std::size_t >( i*columns + j ); };
	for ( int j = 0; j < columns; ++j )
		mask[ at( 1, j ) ] = 1;

	int const bridge = 1*columns + 4;

	g.fill( mask, 1*columns + 0 );
	BOOST_TEST( g.componentCount() == 1 );
	BOOST_TEST( g.componentNodes() == columns,
	            "the fill reached " << g.componentNodes() << " of " << columns
	            << " cells through the bridge, so the leak past a single "
	            "bridging element is the WHOLE far lobe and not one element" );

	// One blocked node, and the far lobe is gone.
	g.fill( mask, 1*columns + 0, bridge );
	BOOST_TEST( g.componentCount() == 2 );
	BOOST_TEST( g.componentNodes() == 4,
	            "blocking the bridge left " << g.componentNodes()
	            << " cells where 4 were expected" );
	BOOST_TEST( g.holds( 1*columns + 3 ) );
	BOOST_TEST( !g.holds( bridge ) );
	BOOST_TEST( !g.holds( 1*columns + 5 ) );
}


/*
 * THE SEED DECIDES WHICH LOBE IS THE PLASMA, AND NOTHING ELSE CAN.
 *
 * A component is not intrinsically "the plasma" -- the largest one need not be,
 * and on a diverted machine the private flux region can be substantial. The
 * seed is the element holding psi_ax, which is the only thing in the problem
 * that says where the plasma is.
 */
BOOST_AUTO_TEST_CASE( the_seed_names_the_component )
{
	PlasmaComponent g = grid( 1, 7 );

	// Two lobes: cells 0..1 and cells 4..6. Cell 2 and 3 are empty.
	std::vector< char > mask{ 1, 1, 0, 0, 1, 1, 1 };

	g.fill( mask, 0 );
	BOOST_TEST( g.componentCount() == 2 );
	BOOST_TEST( g.componentNodes() == 2 );
	BOOST_TEST( g.candidateNodes() == 5 );
	BOOST_TEST( g.holds( 1 ) );
	BOOST_TEST( !g.holds( 4 ) );

	// The SAME mask, seeded in the other lobe, gives the other answer. The
	// larger component is not privileged.
	g.fill( mask, 5 );
	BOOST_TEST( g.componentCount() == 2 );
	BOOST_TEST( g.componentNodes() == 3 );
	BOOST_TEST( !g.holds( 1 ) );
	BOOST_TEST( g.holds( 4 ) );
}


/*
 * A SEED OUTSIDE THE PLASMA GIVES AN EMPTY MASK, LOUDLY.
 *
 * If the element holding psi_ax does not carry Psi > 0 then the normalisation
 * and the field disagree about where the plasma is, and the honest answer is
 * that the confined source is switched off everywhere -- a vacuum field, which
 * is visible in one glance. Quietly falling back to the pointwise test would
 * hide a real inconsistency behind a plausible answer.
 */
BOOST_AUTO_TEST_CASE( a_seed_outside_the_candidates_holds_nothing )
{
	PlasmaComponent g = grid( 1, 5 );
	std::vector< char > mask{ 0, 1, 1, 0, 0 };

	g.fill( mask, 0 );
	BOOST_TEST( g.active() );
	BOOST_TEST( g.seedLabel() == -1 );
	BOOST_TEST( g.componentNodes() == 0 );
	BOOST_TEST( g.candidateNodes() == 2 );
	BOOST_TEST( g.componentCount() == 1 );
	for ( int e = 0; e < 5; ++e )
		BOOST_TEST( !g.holds( e ) );
}


/*
 * ONE COMPONENT MEANS THE CONNECTIVITY TEST CHANGED NOTHING, AND IT MUST SAY SO.
 *
 * This is the control every measurement of the fill needs: on a configuration
 * whose level set really is connected, holds() must agree with the candidate
 * mask exactly. Without it a fill that quietly dropped elements would look like
 * a fix.
 */
BOOST_AUTO_TEST_CASE( a_connected_candidate_set_is_unchanged )
{
	PlasmaComponent g = grid( 5, 5 );

	// A plus sign: connected through the centre, and not convex, so a fill that
	// only walked a bounding box would get it wrong.
	std::vector< char > mask( 25, 0 );
	auto at = [ & ]( int i, int j ) { return static_cast< std::size_t >( i*5 + j ); };
	for ( int j = 0; j < 5; ++j ) { mask[ at( 2, j ) ] = 1; mask[ at( j, 2 ) ] = 1; }

	g.fill( mask, at( 2, 2 ) );
	BOOST_TEST( g.componentCount() == 1 );
	BOOST_TEST( g.candidateNodes() == 9 );
	BOOST_TEST( g.componentNodes() == 9 );
	for ( int e = 0; e < 25; ++e )
		BOOST_TEST( g.holds( e ) == ( mask[ static_cast< std::size_t >( e ) ] != 0 ) );
}


/*
 * THE TWO-RULE FILL: TRAVERSE THE INTERIOR, SHARE OUT THE BAND.
 *
 * The one-rule fill above is bridged by a single element, and on a real mesh
 * that bridge is the band STRADDLING the level set near a saddle -- which is
 * several elements wide, so the leak is the whole far lobe rather than the one
 * element section 10.3 predicted. The two-rule fill walks only elements that
 * are unambiguously inside, where the two lobes are not adjacent at all, and
 * then gives each straddling element to whichever interior component reaches it
 * first, so that the plasma EDGE is not eaten.
 *
 * Here that is a 1 x 11 strip: interior at both ends, a five-cell straddling
 * band in the middle, and nothing to tell the two ends apart except that the
 * fill may not walk the band.
 */
BOOST_AUTO_TEST_CASE( the_two_rule_fill_walks_the_interior_and_shares_the_band )
{
	PlasmaComponent g = grid( 1, 11 );

	//               0  1  2  3  4  5  6  7  8  9 10
	// interior      I  I  I  .  .  .  .  .  I  I  I
	// candidates    C  C  C  C  C  C  C  C  C  C  C
	std::vector< char > interior{ 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1 };
	std::vector< char > candidates( 11, 1 );

	g.fill( interior, candidates, 0 );

	// TWO components of the interior, and the seed names one of them.
	BOOST_TEST( g.componentCount() == 2 );

	// AND candidateNodes() REPORTS THE INCLUSIVE RULE, not the one that was
	// walked: it is what the pointwise support test would have selected, and it
	// is what the fill has to be compared against for the comparison to mean
	// anything.
	BOOST_TEST( g.candidateNodes() == 11 );

	// THE BAND SPLITS DOWN THE MIDDLE. Five cells, so the near lobe takes 3, 4
	// and 5 -- the wave from cell 2 reaches 5 at three hops and the wave from
	// cell 8 reaches it at three too, and the tie goes to the lower component.
	BOOST_TEST( g.holds( 3 ) );
	BOOST_TEST( g.holds( 4 ) );
	BOOST_TEST( !g.holds( 6 ) );
	BOOST_TEST( !g.holds( 7 ) );
	for ( int e = 8; e <= 10; ++e )
		BOOST_TEST( !g.holds( e ),
		            "cell " << e << " -- the far lobe's own interior -- is held. "
		            "The watershed never re-assigns an interior element, so this "
		            "means the traversal itself crossed the band" );
}


/*
 * ONE INTERIOR COMPONENT MEANS THE WHOLE BAND, WHICH IS THE OTHER HALF AND IS
 * WHY THIS IS A WATERSHED RATHER THAN A FIXED NUMBER OF RINGS.
 *
 * A confined plasma on a fitted domain has psi = 0 on the whole boundary, so
 * every element touching it straddles and NONE of them is interior. That band
 * is three elements thick at a corner, and a rule that took a fixed two rings
 * would drop the third -- measured on a real confined rectangle, 4 of 512
 * elements, moving psi_ax by 1.5e-04 and costing seven Newton steps. With one
 * component there is nothing to share with, so the wave takes the band entire.
 */
BOOST_AUTO_TEST_CASE( a_single_interior_component_keeps_the_whole_band )
{
	PlasmaComponent g = grid( 1, 11 );

	// Interior only in the middle; a five-cell band on the left and a
	// three-cell band on the right, both attached to it and neither reachable
	// in the same number of hops.
	std::vector< char > interior{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
	std::vector< char > candidates( 11, 1 );

	g.fill( interior, candidates, 5 );

	BOOST_TEST( g.componentCount() == 1 );
	BOOST_TEST( g.ringNodes() == 10 );
	BOOST_TEST( g.componentNodes() == 11 );
	for ( int e = 0; e <= 10; ++e )
		BOOST_TEST( g.holds( e ),
		            "cell " << e << " was dropped from a band with only ONE "
		            "interior component to give it to. There is nothing to share "
		            "with, so the whole band belongs to the plasma however far "
		            "from the interior it is" );
}


/*
 * AND A STRADDLING ELEMENT ATTACHED TO NO INTERIOR AT ALL IS DROPPED, WHICH IS
 * THE ONE THING THE WATERSHED CANNOT KEEP AND IS WORTH KNOWING.
 *
 * A pocket thinner than one element everywhere has no interior for a wave to
 * start from, so it is left out entirely. That is a mesh too coarse to resolve
 * the plasma rather than a case to accommodate -- and
 * GradShafranovSolver::refreshPlasmaComponent falls back to the one-rule fill
 * when the SEED itself has no interior, so a plasma that is thin everywhere
 * gets the inclusive support MEQ had before XP-1 rather than nothing at all.
 */
BOOST_AUTO_TEST_CASE( a_band_touching_no_interior_is_dropped )
{
	PlasmaComponent g = grid( 1, 9 );

	// Interior at 0..2; a detached candidate pocket at 6..7 with a gap at 4..5.
	std::vector< char > interior{ 1, 1, 1, 0, 0, 0, 0, 0, 0 };
	std::vector< char > candidates{ 1, 1, 1, 1, 0, 0, 1, 1, 0 };

	g.fill( interior, candidates, 0 );

	BOOST_TEST( g.holds( 3 ) );
	BOOST_TEST( !g.holds( 6 ) );
	BOOST_TEST( !g.holds( 7 ) );
	BOOST_TEST( g.componentNodes() == 4 );
	BOOST_TEST( g.candidateNodes() == 6 );
}


/*
 * A RING CANNOT BE TAKEN FROM A SET THAT DOES NOT CONTAIN THE INTERIOR./*
 * A RING CANNOT BE TAKEN FROM A SET THAT DOES NOT CONTAIN THE INTERIOR.
 *
 * Two rules that disagree about which is stricter would give a mask with HOLES
 * in it -- kept elements with dropped elements inside them -- and a mask with
 * holes is a plausible answer that nothing downstream can check. Refused here,
 * where the message can name what is wrong.
 */
BOOST_AUTO_TEST_CASE( the_two_rules_must_agree_about_which_is_stricter )
{
	PlasmaComponent g = grid( 1, 5 );

	std::vector< char > interior{ 1, 1, 1, 0, 0 };
	std::vector< char > candidates{ 1, 1, 0, 0, 0 };

	BOOST_CHECK_THROW( g.fill( interior, candidates, 0 ), std::invalid_argument );
	BOOST_CHECK_THROW( g.fill( interior, std::vector< char >( 3, 1 ), 0 ),
	                   std::invalid_argument );
}


/*
 * A MALFORMED GRAPH THROWS RATHER THAN PRODUCING A PLAUSIBLE MASK.
 *
 * The output of this class is a list of element numbers, and nothing
 * downstream -- not a residual, not a convergence rate, not a picture -- can
 * tell a wrong one from a right one. So the checks are here, where they are
 * cheap and where the message can name what is wrong.
 */
BOOST_AUTO_TEST_CASE( a_malformed_graph_is_refused )
{
	PlasmaComponent g;

	BOOST_CHECK_THROW( g.setAdjacency( {}, {} ), std::invalid_argument );
	BOOST_CHECK_THROW( g.setAdjacency( { 1, 2 }, { 0, 0 } ), std::invalid_argument );
	BOOST_CHECK_THROW( g.setAdjacency( { 0, 2, 1 }, { 0, 1 } ), std::invalid_argument );
	BOOST_CHECK_THROW( g.setAdjacency( { 0, 1 }, {} ), std::invalid_argument );
	BOOST_CHECK_THROW( g.setAdjacency( { 0, 1, 2 }, { 1, 7 } ), std::invalid_argument );

	PlasmaComponent good = grid( 2, 2 );
	BOOST_CHECK_THROW( good.fill( std::vector< char >( 3, 1 ), 0 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( good.fill( std::vector< char >( 4, 1 ), 4 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( good.fill( std::vector< char >( 4, 1 ), -1 ),
	                   std::invalid_argument );
}
