#define BOOST_TEST_MODULE AdaptiveRefinement
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "mfem.hpp"

#include "meq/Estimator.hpp"
#include "meq/GradShafranov.hpp"
#include "meq/Source.hpp"

#include "analytic/Soloviev.hpp"

/*
 * The other half of stage 6: the solve -> estimate -> mark -> refine loop of
 * refs/HDG-GradShafranov-Adaptive.pdf section 3, and the computational-domain
 * update of its section 3.3.
 *
 * Two claims are checked here, and they are deliberately checked apart, because
 * they can fail independently and because only one of them is about the
 * estimator.
 *
 * 1. THE DOMAIN UPDATE IS NECESSARY AND IT WORKS. Section 3.3's whole argument
 *    is that refining T_h alone lets dist( Gamma_h, Gamma ) drift away from
 *    O( h_loc ): the children of an element inside Omega are still inside Omega,
 *    so Gamma_h does not move while h_loc halves. The companion mesh is what
 *    fixes it. Measured below on geometry alone -- no solver, no estimator -- so
 *    this test says something even if the estimator is wrong. Three treatments,
 *    marking the elements along Gamma_h on an eight-cell background, giving
 *    max_y d( y, Gamma ) / h_loc( y ):
 *
 *      cycle              0       1       2       3
 *      companion       0.98    1.02    1.39    1.36     bounded
 *      background      0.98    1.40    1.99    1.99     half the effect, by luck
 *      frozen          0.98    1.97    3.94    7.88     doubling, i.e. Figure 3
 *
 *    with the gap itself constant at 0.29595 to every digit in the frozen case
 *    and falling 0.296 -> 0.154 -> 0.104 -> 0.051 with the companion update.
 *
 *    The middle row is the one worth knowing about, because it is the version
 *    somebody would ship by mistake. "Background" is steps 4 and 5 with step 3's
 *    second half -- mark the band Gamma cuts as well -- left out. It is NOT the
 *    same as doing nothing: MFEM's conforming bisection propagates into a marked
 *    element's neighbours to keep the mesh conforming, some of those neighbours
 *    are in the band, and re-selection then admits a few of their children. So it
 *    recovers about half the effect for free and looks as though it works. The
 *    frozen row -- the computational mesh refined as a standalone object, which is
 *    what Figure 3 actually draws -- is the honest control, and it fails exactly
 *    as the paper says, at a clean factor of two per cycle.
 *
 * 2. THE LOOP RUNS, AND THE ESTIMATOR STAYS HONEST ON A GRADED MESH. eta was
 *    verified in tests/convergence/EstimatorConvergence.cpp on uniform meshes
 *    only, where every h_K and h_e is the same number. A residual estimator can
 *    be right there and wrong on a mesh with a spread of element sizes, so the
 *    effectivity index is measured again on the adapted meshes. It is, and that
 *    is the one thing here that a plain "the adaptive loop converges" run would
 *    not have said.
 *
 * 3. AND THE TWO TOGETHER, on the curved boundary, which is the configuration the
 *    paper is actually about. That one is asserted more weakly -- eta and the
 *    true error come down every cycle, the proximity condition holds, and the
 *    transfer paths stay admissible -- because the sequence of computational
 *    domains is not a dyadic refinement of anything, so there is no h to take a
 *    rate against.
 *
 * WHAT MADE 3 POSSIBLE, AND WHAT IT WOULD HAVE LOOKED LIKE OTHERWISE.
 *
 * eta_5 on Gamma_h compares psi* against a trace value that was pinned to zero
 * rather than being the phi_h actually imposed, so the difference there is
 * O( dist( Gamma_h, Gamma ) ) = O( h ) and not O( h^(k+2) ). At k = 2 that makes
 * eta 4.09e-1 where eta_1 is 2.12e-3, converging at about a half: the estimator
 * becomes that one term, and the marking is driven by it. A loop built on that
 * would have run, produced pictures, and refined the wrong elements -- which is
 * precisely the failure mode a "the adaptive loop converges" acceptance test
 * cannot see. meq::ResidualEstimator::setTransferredBoundary() leaves those faces
 * out, and eta and all five components then converge at k+1 on the extension
 * path too. The other worry, that psi* itself would be wrong there because
 * DarcyForm::Reconstruct() drops the boundary-face integrator carrying the datum,
 * turned out to be groundless: psi* still converges at k+2. Both were measured.
 *
 * MAXIMUM MARKING AT GAMMA = 0.3 IS UNIFORM REFINEMENT ON THIS PROBLEM, and that
 * is not a broken marker. Every element's eta_K is within a factor of 0.3 of the
 * largest, because the Solov'ev solution is smooth and this rectangle cuts no
 * feature out of it, so the criterion marks everything at every cycle. It is
 * section 3.2's own statement -- maximum marking "becomes uniform as the
 * parameter approaches zero" -- meeting an indicator that was nearly uniform
 * already. GS-2 used 0.3 on domains where an unfitted boundary and internal
 * layers spread eta_K over orders of magnitude. The loop below runs it anyway,
 * because it is the paper's choice, and runs 0.8 beside it for the graded case.
 *
 * THE WARNING AT THE END OF SECTION 3.3, WHICH APPLIES TO THIS FILE'S OUTPUT.
 * The refinement will look as though it is crowding Gamma. That is the local
 * proximity condition at work: step 3 of the update marks the band Gamma cuts
 * whatever the indicator says, so new and necessarily smaller elements appear
 * there by construction. It is not the error estimator concentrating there, and
 * reading it as the latter sends you hunting a bug that does not exist. The
 * "prox" column of the first table is that rule's contribution as a number --
 * 29, 64 and 134 elements on the three cycles, against a marked set of about the
 * same size -- so the distinction is measurable rather than a matter of
 * impression.
 */

namespace
{

	meq::analytic::SolovievEquilibrium const &equilibrium()
	{
		static meq::analytic::SolovievEquilibrium const eq
			= meq::analytic::SolovievEquilibrium::nstx();
		return eq;
	}

	// ---------------------------------------------------------------------
	// The curved geometry, as tests/convergence/ExtensionConvergence.cpp sets
	// it up: Gamma is the closed flux surface psi_nstx = -0.03, written as the
	// zero set of psi := psi_nstx + 0.03. The separatrix itself passes through
	// an X-point, which is a corner of Gamma where the transfer paths give out;
	// that file argues it out in full.
	// ---------------------------------------------------------------------
	double const psiOffset = 0.03;

	double psiExact( double r, double z )
	{
		return equilibrium().psi( r, z ) + psiOffset;
	}

	double levelSet( mfem::Vector const &x )
	{
		return psiExact( x( 0 ), x( 1 ) );
	}

	double const boxRMin = 0.25;
	double const boxRMax = 1.95;
	double const boxZMin = -1.75;
	double const boxZMax = 1.65;

	mfem::Mesh makeBackground( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			n, 2*n, mfem::Element::TRIANGLE, false,
			boxRMax - boxRMin, boxZMax - boxZMin );
		mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + boxRMin;
			out( 1 ) = in( 1 ) + boxZMin;
		} );
		return mesh;
	}

	/// The distance from a point of Gamma_h to Gamma, along the outward normal of
	/// the level set, by marching and then bisecting.
	///
	/// This is what section 2.2's d( y, y-bar ) is, computed from the level set
	/// alone. It deliberately does NOT go through mfem::TransferPath: the claim
	/// being tested is about where Gamma_h sits relative to Gamma, and measuring
	/// it with the same machinery that is supposed to bridge the gap would make
	/// the measurement depend on the thing it is meant to be independent of.
	double distanceToGamma( double r, double z, double searchLength )
	{
		double const step = 1.0e-5;
		double const dr = ( psiExact( r + step, z ) - psiExact( r - step, z ) )/( 2.0*step );
		double const dz = ( psiExact( r, z + step ) - psiExact( r, z - step ) )/( 2.0*step );
		double const norm = std::sqrt( dr*dr + dz*dz );
		if ( norm <= 0.0 )
			return searchLength;

		// psi increases outward across Gamma, so +grad points out of Omega.
		double const ur = dr/norm;
		double const uz = dz/norm;

		double lo = 0.0;
		double hi = -1.0;
		int const marches = 400;
		for ( int i = 1; i <= marches; ++i )
		{
			double const s = searchLength*static_cast<double>( i )/marches;
			if ( psiExact( r + s*ur, z + s*uz ) > 0.0 )
			{
				hi = s;
				lo = searchLength*static_cast<double>( i - 1 )/marches;
				break;
			}
		}
		if ( hi < 0.0 )
			return searchLength;

		for ( int i = 0; i < 60; ++i )
		{
			double const mid = 0.5*( lo + hi );
			if ( psiExact( r + mid*ur, z + mid*uz ) > 0.0 )
				hi = mid;
			else
				lo = mid;
		}
		return 0.5*( lo + hi );
	}

	/// The two halves of the proximity condition, kept apart: the largest gap
	/// d( y, Gamma ) over the vertices y of Gamma_h, and the largest ratio
	/// d( y, Gamma ) / h_loc( y ). Section 2.2 requires the ratio to stay O( 1 );
	/// section 3.3 says plain refinement of T_h does not manage it.
	///
	/// Reported apart because they fail differently. The ratio is the condition;
	/// the gap says WHY it fails, since h_loc always shrinks and the question is
	/// whether Gamma_h follows it in.
	///
	/// Takes a plain mesh and the attribute of Gamma_h on it, not an
	/// AdaptiveDomain, so that the same measurement can be made on a
	/// computational mesh refined as a standalone object -- which is the control
	/// the paper's Figure 3 actually draws.
	struct Proximity
	{
		double gap;
		double ratio;
		double smallest;
	};

	Proximity worstProximity( mfem::Mesh &mesh, int gammaH, double searchLength )
	{
		// h_loc( y ), the smallest diameter over the elements containing y, in one
		// pass over the elements rather than one pass per vertex.
		std::vector<double> local( mesh.GetNV(), 0.0 );
		mfem::Array<int> vertices;
		double smallest = 0.0;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			double const h = meq::elementDiameter( mesh, e );
			smallest = ( e == 0 ) ? h : std::min( smallest, h );
			mesh.GetElementVertices( e, vertices );
			for ( int v = 0; v < vertices.Size(); ++v )
			{
				double &slot = local[ vertices[ v ] ];
				slot = ( slot == 0.0 ) ? h : std::min( slot, h );
			}
		}

		std::vector<char> onGammaH( mesh.GetNV(), 0 );
		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			if ( mesh.GetBdrAttribute( be ) != gammaH )
				continue;
			mesh.GetBdrElementVertices( be, vertices );
			for ( int v = 0; v < vertices.Size(); ++v )
				onGammaH[ vertices[ v ] ] = 1;
		}

		Proximity worst;
		worst.gap = 0.0;
		worst.ratio = 0.0;
		worst.smallest = smallest;
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			if ( !onGammaH[ v ] || local[ v ] <= 0.0 )
				continue;
			double const *x = mesh.GetVertex( v );
			double const gap = distanceToGamma( x[ 0 ], x[ 1 ], searchLength );
			worst.gap = std::max( worst.gap, gap );
			worst.ratio = std::max( worst.ratio, gap/local[ v ] );
		}
		return worst;
	}

	/// The elements with a face on Gamma_h.
	///
	/// A geometric marking, not an error indicator, and that is the point: the
	/// claim of section 3.3 is about what refinement does to the geometry, so the
	/// marked set is chosen to isolate it. It is also the worst case, since it is
	/// exactly the elements whose refinement leaves Gamma_h where it was.
	void markAlongGammaH( mfem::Mesh &mesh, int gammaH, mfem::Array<int> &marked )
	{
		std::vector<char> flag( mesh.GetNE(), 0 );
		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			if ( mesh.GetBdrAttribute( be ) != gammaH )
				continue;
			int const face = mesh.GetBdrElementFaceIndex( be );
			mfem::FaceElementTransformations *ftr = mesh.GetFaceElementTransformations( face );
			if ( ftr )
				flag[ ftr->Elem1No ] = 1;
		}

		marked.SetSize( 0 );
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			if ( flag[ e ] )
				marked.Append( e );
		}
	}

	// ---------------------------------------------------------------------
	// The fitted rectangle, for the loop itself. The same one
	// SolovievConvergence.cpp and EstimatorConvergence.cpp use, so the
	// estimator running on it has been measured against a closed form.
	// ---------------------------------------------------------------------
	double const rMin = 0.6;
	double const rMax = 1.4;
	double const zMin = -0.6;
	double const zMax = 0.6;

	mfem::Mesh makeFittedMesh( int n )
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( n, n, mfem::Element::TRIANGLE, false,
		                                               rMax - rMin, zMax - zMin );
		mesh.Transform( []( mfem::Vector const &in, mfem::Vector &out )
		{
			out( 0 ) = in( 0 ) + rMin;
			out( 1 ) = in( 1 ) + zMin;
		} );
		return mesh;
	}

	enum class Marking
	{
		Doerfler,
		Maximum
	};

	struct Cycle
	{
		int elements;
		int traceDofs;
		int marked;
		double eta;
		double errorFlux;
		double errorPsi;
		double largest;
		double smallest;
	};

	/// One turn of solve -> estimate -> mark -> refine on the fitted mesh.
	///
	/// A fresh solver every cycle. GradShafranovSolver builds its three spaces in
	/// its constructor and DarcyForm's hybridization takes what it finds when
	/// EnableHybridization() runs, so a refined mesh needs a new solver rather
	/// than an Update() -- which is what CLAUDE.md records the deleted
	/// Solution::Prolong() and Update() as having been for. Nothing here needs the
	/// previous iterate: the linear path has no initial guess to carry, so the
	/// cheap route is also the correct one. A semi-linear adaptive run would want
	/// the prolongation, and that is a separate piece of work.
	std::vector<Cycle> loop( int order, int start, int cycles, Marking marking,
	                         double gamma )
	{
		meq::analytic::SolovievEquilibrium const &eq = equilibrium();
		meq::SolovievSource const source( eq.getA() );

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

		mfem::Mesh mesh = makeFittedMesh( start );
		std::vector<Cycle> history;

		for ( int c = 0; c < cycles; ++c )
		{
			meq::GradShafranovSolver solver( mesh, order );
			solver.setSource( sourceCoeff );
			solver.setBoundaryData( psiCoeff );
			solver.solve();
			solver.postProcess();

			meq::ResidualEstimator estimator( solver, source );
			mfem::Vector const &local = estimator.GetLocalErrors();

			mfem::Array<int> marked;
			if ( marking == Marking::Doerfler )
				meq::markDoerfler( local, gamma, marked );
			else
				meq::markMaximum( local, gamma, marked );

			Cycle entry;
			entry.elements = mesh.GetNE();
			entry.traceDofs = solver.numTraceDofs();
			entry.marked = marked.Size();
			entry.eta = estimator.GetTotalError();
			entry.errorFlux = solver.fluxError( fluxCoeff );
			entry.errorPsi = solver.potentialError( psiCoeff );
			entry.largest = 0.0;
			entry.smallest = 0.0;
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				double const h = meq::elementDiameter( mesh, e );
				entry.largest = std::max( entry.largest, h );
				entry.smallest = ( e == 0 ) ? h : std::min( entry.smallest, h );
			}
			history.push_back( entry );

			if ( c + 1 == cycles )
				break;
			BOOST_TEST_REQUIRE( marked.Size() > 0,
			                    "nothing was marked at cycle " << c
			                    << ", so the loop cannot continue" );
			mesh.GeneralRefinement( marked );
		}

		return history;
	}

}

/// Section 3.3, measured both ways. Refining T_h alone lets d( y, Gamma )/h_loc
/// grow without bound; the companion update holds it.
///
/// Marked geometrically, along Gamma_h, so that the only thing under test is what
/// refinement does to the domain. No solver and no estimator appear.
BOOST_AUTO_TEST_CASE( theCompanionUpdateHoldsTheProximityCondition )
{
	int const cells = 8;
	int const cycles = 3;
	double const searchLength = 12.0*( boxRMax - boxRMin )/cells;

	struct Row
	{
		int elements;
		int companion;
		int proximity;
		Proximity worst;
	};

	// Three treatments, because two of them turn out to be different and the
	// difference is worth knowing:
	//
	//   Companion  section 3.3 in full -- mark the band Gamma cuts as well, refine
	//              the background, re-select T_h and T_c^h.
	//   Background steps 4 and 5 with step 3's second half left out. NOT the same
	//              as doing nothing: MFEM's conforming bisection propagates into
	//              the neighbours of a marked element to keep the mesh conforming,
	//              and some of those neighbours are in the cut band, so a little
	//              of the companion update happens by accident. Measured, it
	//              recovers about half of the effect.
	//   Frozen     the computational mesh refined as a standalone object, with no
	//              re-selection at all. This is what Figure 3 draws, and Omega_h
	//              is then genuinely fixed.
	enum class Treatment { Companion, Background, Frozen };
	Treatment const treatments[ 3 ] =
		{ Treatment::Companion, Treatment::Background, Treatment::Frozen };
	char const *names[ 3 ] = { "companion", "background", "frozen" };

	std::vector<Row> table[ 3 ];

	for ( int variant = 0; variant < 3; ++variant )
	{
		Treatment const treatment = treatments[ variant ];

		mfem::Mesh background = makeBackground( cells );
		meq::AdaptiveDomain domain( background, levelSet );

		// The frozen case detaches a copy of T_h and never looks at the
		// background again.
		mfem::Mesh frozen( domain.computational(), true );
		int const frozenGammaH = domain.gammaHAttribute();

		for ( int c = 0; c <= cycles; ++c )
		{
			mfem::Mesh &current = ( treatment == Treatment::Frozen )
				? frozen : static_cast<mfem::Mesh &>( domain.computational() );
			int const gammaH = ( treatment == Treatment::Frozen )
				? frozenGammaH : domain.gammaHAttribute();

			Row row;
			row.elements = current.GetNE();
			row.companion = domain.numCompanion();
			row.proximity = domain.lastProximityAdditions();
			row.worst = worstProximity( current, gammaH, searchLength );
			table[ variant ].push_back( row );

			if ( c == cycles )
				break;

			mfem::Array<int> marked;
			markAlongGammaH( current, gammaH, marked );
			BOOST_TEST_REQUIRE( marked.Size() > 0,
			                    names[ variant ] << ": no element has a face on Gamma_h at "
			                    "cycle " << c );

			switch ( treatment )
			{
				case Treatment::Companion:  domain.refine( marked ); break;
				case Treatment::Background: domain.refineWithoutCompanion( marked ); break;
				case Treatment::Frozen:     frozen.GeneralRefinement( marked ); break;
			}
		}
	}

	std::printf( "\n  section 3.3, three ways. Gamma_h against Gamma over %d cycles of "
	             "refinement marked along Gamma_h.\n", cycles );
	for ( int variant = 0; variant < 3; ++variant )
	{
		std::printf( "    %-11s %6s %8s %6s %9s %9s %8s\n",
		             names[ variant ], "cycle", "elem", "prox", "h_min", "max gap",
		             "gap/h" );
		for ( std::size_t c = 0; c < table[ variant ].size(); ++c )
		{
			Row const &row = table[ variant ][ c ];
			std::printf( "    %-11s %6zu %8d %6d %9.5f %9.5f %8.4f\n",
			             "", c, row.elements, row.proximity, row.worst.smallest,
			             row.worst.gap, row.worst.ratio );
		}
	}
	std::fflush( stdout );

	Row const &firstCompanion = table[ 0 ].front();
	Row const &lastCompanion = table[ 0 ].back();
	Row const &firstFrozen = table[ 2 ].front();
	Row const &lastFrozen = table[ 2 ].back();
	Row const &lastBackground = table[ 1 ].back();

	// Figure 3, exactly: with Omega_h held fixed, Gamma_h does not move at all,
	// so the gap is unchanged to the last bit while h_loc halves and the ratio
	// grows with it. This is the assertion that needs no asymptotics and no
	// tolerance.
	BOOST_TEST( std::abs( lastFrozen.worst.gap - firstFrozen.worst.gap )
	            < 1.0e-12*firstFrozen.worst.gap,
	            "with Omega_h frozen the largest gap moved from " << firstFrozen.worst.gap
	            << " to " << lastFrozen.worst.gap
	            << " -- refining a mesh cannot move its own boundary, so the "
	            "measurement is wrong" );
	BOOST_TEST( lastFrozen.worst.ratio > 3.0*firstFrozen.worst.ratio,
	            "with Omega_h frozen gap/h_loc went from " << firstFrozen.worst.ratio
	            << " to " << lastFrozen.worst.ratio << " over " << cycles
	            << " cycles, so h_loc is not shrinking and this control has stopped "
	            "being a control" );

	// And with the companion update it stays bounded. O( h_loc ) is what the
	// analysis asks for, not a particular constant: the staircase boundary of a
	// subdomain of a triangulation has faces whose normal is nearly tangent to
	// Gamma, which is what makes the constant bigger than one.
	// ExtensionConvergence.cpp measures 1.3 h on a uniform mesh by a different
	// route and bounds it at 3.
	for ( std::size_t c = 0; c < table[ 0 ].size(); ++c )
		BOOST_TEST( table[ 0 ][ c ].worst.ratio < 3.0,
		            "with the companion update, cycle " << c << " has gap/h_loc = "
		            << table[ 0 ][ c ].worst.ratio
		            << ", so the proximity condition has been lost" );
	BOOST_TEST( lastCompanion.worst.ratio < 1.6*firstCompanion.worst.ratio,
	            "with the companion update gap/h_loc went from "
	            << firstCompanion.worst.ratio << " to " << lastCompanion.worst.ratio
	            << " over " << cycles << " cycles, which is a drift rather than a bound" );

	// The mechanism: Gamma_h closes on Gamma at the rate the mesh refines,
	// because the cut band subdivides and its children are admitted to T_h.
	BOOST_TEST( lastCompanion.worst.gap < 0.35*firstCompanion.worst.gap,
	            "with the companion update the largest gap went from "
	            << firstCompanion.worst.gap << " to " << lastCompanion.worst.gap
	            << " over " << cycles << " cycles" );
	BOOST_TEST( lastCompanion.elements > 5*firstCompanion.elements,
	            "T_h has " << lastCompanion.elements << " elements after " << cycles
	            << " cycles against " << firstCompanion.elements
	            << ", so the companion refinement is not feeding new elements into "
	            "Omega" );

	// The companion update beats the frozen control outright, and it also beats
	// leaving step 3's second half out -- which is the one that would be easy to
	// ship by mistake, because conforming bisection's propagation makes it look
	// as though it works.
	BOOST_TEST( lastFrozen.worst.ratio > 2.5*lastCompanion.worst.ratio,
	            "the companion update and the frozen domain end at "
	            << lastCompanion.worst.ratio << " and " << lastFrozen.worst.ratio
	            << ", which is not the difference section 3.3 exists to make" );
	BOOST_TEST( lastBackground.worst.ratio > 1.3*lastCompanion.worst.ratio,
	            "leaving the proximity rule out ends at " << lastBackground.worst.ratio
	            << " against " << lastCompanion.worst.ratio << " with it, so step 3's "
	            "second half is not earning its place -- check whether the bisection "
	            "propagation is doing all the work" );

	// The warning at the end of section 3.3, as a number. A good part of what
	// gets refined near Gamma is put there by the proximity rule and not by the
	// marked set, so a picture of the output will look as though the indicator is
	// concentrating at the boundary when it is not.
	BOOST_TEST( lastCompanion.proximity > 0,
	            "the proximity rule added no elements on the last cycle, so step 3's "
	            "second half is doing nothing and refine() is plain AMR" );
}

/// The loop, on the geometry where the estimator has been measured. eta and the
/// true error must both come down every cycle, and the effectivity index must
/// stay where EstimatorConvergence.cpp found it -- which is a new statement,
/// because that file only ever put the estimator on a uniform mesh.
BOOST_AUTO_TEST_CASE( theAdaptiveLoopRunsAndTheEstimatorStaysHonest )
{
	int const order = 2;
	int const cycles = 5;

	struct Variant
	{
		char const *name;
		Marking marking;
		double gamma;
		bool expectGraded;
	};

	// GS-2's own experiments used maximum marking at gamma = 0.3; the convergence
	// proof of Cockburn, Nochetto and Zhang assumes Doerfler. Both are run, and so
	// is maximum marking at a parameter that actually grades this problem.
	//
	// MEASURED, and worth recording because it looks at first like a broken
	// marker: maximum marking at gamma = 0.3 marks EVERY element here, at every
	// cycle, and the loop degenerates to uniform refinement. That is correct
	// behaviour and it is section 3.2's own statement -- maximum marking "becomes
	// uniform as the parameter approaches zero" -- applied to a problem where the
	// indicator is nearly uniform to begin with. The Solov'ev solution is smooth
	// and this rectangle cuts no interesting feature out of it, so no element's
	// eta_K falls below 0.3 of the largest. GS-2 used 0.3 on domains where the
	// unfitted boundary and internal layers spread eta_K over orders of magnitude.
	// So 0.3 is kept as the paper's own choice, with the grading checks turned
	// off, and 0.8 is added as the one that grades.
	std::vector<Variant> const variants = {
		{ "Doerfler, gamma = 0.6", Marking::Doerfler, 0.6, true },
		{ "maximum,  gamma = 0.8", Marking::Maximum, 0.8, true },
		{ "maximum,  gamma = 0.3", Marking::Maximum, 0.3, false }
	};

	for ( Variant const &variant : variants )
	{
		std::vector<Cycle> const history = loop( order, 4, cycles, variant.marking,
		                                         variant.gamma );

		std::printf( "\n  the adaptive loop on the fitted rectangle, k = %d, %s\n",
		             order, variant.name );
		std::printf( "  %6s %7s %8s %7s %11s %11s %11s %8s %8s %8s\n",
		             "cycle", "elem", "trace", "marked", "eta", "L2(q)", "L2(psi)",
		             "eta/L2q", "h_max", "h_min" );
		for ( std::size_t c = 0; c < history.size(); ++c )
		{
			Cycle const &e = history[ c ];
			std::printf( "  %6zu %7d %8d %7d %11.4e %11.4e %11.4e %8.3f %8.5f %8.5f\n",
			             c, e.elements, e.traceDofs, e.marked, e.eta, e.errorFlux,
			             e.errorPsi, e.eta/e.errorFlux, e.largest, e.smallest );
		}
		std::fflush( stdout );

		for ( std::size_t c = 1; c < history.size(); ++c )
		{
			BOOST_TEST( history[ c ].eta < history[ c - 1 ].eta,
			            variant.name << ": eta went from " << history[ c - 1 ].eta
			            << " to " << history[ c ].eta << " at cycle " << c );
			BOOST_TEST( history[ c ].errorFlux < history[ c - 1 ].errorFlux,
			            variant.name << ": the flux error went from "
			            << history[ c - 1 ].errorFlux << " to " << history[ c ].errorFlux
			            << " at cycle " << c
			            << " -- eta came down and the true error did not, which is what a "
			            "misleading indicator looks like" );
		}

		// The effectivity index on a graded mesh. EstimatorConvergence.cpp
		// measures 15.3 for k = 2 on uniform meshes; a residual estimator whose
		// h_K or h_e weight were wrong would be right there and wrong here, since
		// on a uniform mesh every h is the same number and a misplaced power of h
		// is only a constant.
		for ( Cycle const &e : history )
		{
			double const index = e.eta/e.errorFlux;
			BOOST_TEST( ( index > 5.0 && index < 40.0 ),
			            variant.name << ": the effectivity index is " << index
			            << " on a mesh whose elements span " << e.smallest << " to "
			            << e.largest << ", against 15.3 on a uniform mesh at this k" );
		}

		if ( !variant.expectGraded )
			continue;

		// The mesh really is graded by the end, or the effectivity paragraph above
		// is checking nothing new.
		Cycle const &last = history.back();
		BOOST_TEST( last.largest > 1.5*last.smallest,
		            variant.name << ": the final mesh spans only " << last.smallest
		            << " to " << last.largest
		            << ", so it is close to uniform and the graded-mesh check above is "
		            "not exercising anything" );

		// And the loop is adaptive rather than uniform: a cycle that marked
		// everything would be a dyadic refinement in disguise.
		for ( std::size_t c = 0; c + 1 < history.size(); ++c )
			BOOST_TEST( history[ c ].marked < history[ c ].elements,
			            variant.name << ": cycle " << c << " marked all "
			            << history[ c ].elements << " elements, which is uniform "
			            "refinement" );
	}
}

/// The whole of stage 6 at once: solve -> estimate -> mark -> refine on the
/// CURVED boundary, with the computational-domain update of section 3.3 in the
/// refine step and the estimator driving the marking.
///
/// This is the configuration the paper is about, and it only became runnable once
/// eta_5 was told to leave Gamma_h alone -- see
/// meq::ResidualEstimator::setTransferredBoundary(). With those faces in, eta is
/// 4.09e-1 on the coarsest mesh against eta_1 = 2.12e-3, converging at about a
/// half, and the marking is then driven entirely by a term comparing psi* against
/// a trace value that was pinned rather than imposed. A loop built on that would
/// have run, produced pictures, and refined the wrong elements.
///
/// What is asserted is what can be: eta and the true error both come down every
/// cycle, the proximity condition holds throughout, and the transfer paths stay
/// admissible. Not a rate -- the sequence of computational domains is not a
/// dyadic refinement of a fixed domain, so there is no h to take a rate against.
BOOST_AUTO_TEST_CASE( theAdaptiveLoopRunsOnTheCurvedBoundary )
{
	int const order = 2;
	int const cells = 8;
	int const cycles = 4;
	double const gamma = 0.6;

	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	meq::SolovievSource const source( eq.getA() );

	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psiCoeff( []( mfem::Vector const &x )
	{
		return psiExact( x( 0 ), x( 1 ) );
	} );
	mfem::VectorFunctionCoefficient fluxCoeff( 2, [ &eq ]( mfem::Vector const &x,
	                                                       mfem::Vector &value )
	{
		eq.flux( x( 0 ), x( 1 ), value( 0 ), value( 1 ) );
	} );
	mfem::ConstantCoefficient zero( 0.0 );

	struct Step
	{
		int elements;
		int traceDofs;
		int marked;
		int widened;
		double eta;
		double errorPsi;
		double errorFlux;
		double ratio;
	};

	mfem::Mesh background = makeBackground( cells );
	meq::AdaptiveDomain domain( background, levelSet );

	std::vector<Step> history;

	for ( int c = 0; c < cycles; ++c )
	{
		mfem::Array<int> marked;
		Step step;

		{
			// Scoped, because domain.refine() replaces the SubMesh these are all
			// built on. Everything wanted afterwards is copied out into step first.
			mfem::SubMesh &sub = domain.computational();
			int const gammaH = domain.gammaHAttribute();

			// Twelve times the LARGEST element, not the mesh parameter: on a graded
			// mesh the coarse part needs the long search and the fine part is not
			// harmed by being given one. ExtensionConvergence.cpp uses six h on a
			// uniform mesh, where the two are the same number.
			mfem::VertexConePath path( sub, gammaH, levelSet,
			                           12.0*domain.largestElement() );

			meq::GradShafranovSolver solver( sub, order );
			solver.setSource( sourceCoeff );
			solver.setBoundaryData( zero );
			solver.setExtension( path, domain.gammaHMarker() );
			solver.solve();
			solver.postProcess();

			meq::ResidualEstimator estimator( solver, source );
			estimator.setTransferredBoundary( domain.gammaHMarker() );
			mfem::Vector const &local = estimator.GetLocalErrors();

			meq::markDoerfler( local, gamma, marked );

			step.elements = sub.GetNE();
			step.traceDofs = solver.numTraceDofs();
			step.marked = marked.Size();
			step.widened = path.NumWidened();
			step.eta = estimator.GetTotalError();
			step.errorPsi = solver.potentialError( psiCoeff );
			step.errorFlux = solver.fluxError( fluxCoeff );
			step.ratio = worstProximity( sub, gammaH,
			                             12.0*domain.largestElement() ).ratio;
		}

		history.push_back( step );

		if ( c + 1 == cycles )
			break;
		BOOST_TEST_REQUIRE( marked.Size() > 0,
		                    "nothing was marked at cycle " << c );
		domain.refine( marked );
	}

	std::printf( "\n  stage 6 end to end: the curved boundary, k = %d, Doerfler "
	             "gamma = %.1f\n", order, gamma );
	std::printf( "  %6s %7s %8s %7s %5s %11s %11s %11s %8s\n",
	             "cycle", "elem", "trace", "marked", "wide", "eta", "L2(psi)",
	             "L2(q)", "gap/h" );
	for ( std::size_t c = 0; c < history.size(); ++c )
	{
		Step const &s = history[ c ];
		std::printf( "  %6zu %7d %8d %7d %5d %11.4e %11.4e %11.4e %8.4f\n",
		             c, s.elements, s.traceDofs, s.marked, s.widened, s.eta,
		             s.errorPsi, s.errorFlux, s.ratio );
	}
	std::fflush( stdout );

	for ( std::size_t c = 1; c < history.size(); ++c )
	{
		BOOST_TEST( history[ c ].eta < history[ c - 1 ].eta,
		            "eta went from " << history[ c - 1 ].eta << " to " << history[ c ].eta
		            << " at cycle " << c );
		BOOST_TEST( history[ c ].errorPsi < history[ c - 1 ].errorPsi,
		            "the L2 error in psi went from " << history[ c - 1 ].errorPsi
		            << " to " << history[ c ].errorPsi << " at cycle " << c
		            << " -- eta came down and the true error did not" );
		BOOST_TEST( history[ c ].errorFlux < history[ c - 1 ].errorFlux,
		            "the L2 error in q went from " << history[ c - 1 ].errorFlux
		            << " to " << history[ c ].errorFlux << " at cycle " << c );
	}

	for ( std::size_t c = 0; c < history.size(); ++c )
	{
		// The domain update is doing its job throughout the run, not just in the
		// geometry-only test above.
		BOOST_TEST( history[ c ].ratio < 3.0,
		            "cycle " << c << " has gap/h_loc = " << history[ c ].ratio
		            << ", so the proximity condition was lost during the loop" );

		// Assumption P.1 of the Cockburn-Solano analysis: every vertex of Gamma_h
		// needs a path leaving D_h through both faces meeting there. VertexConePath
		// widens the fan when it cannot and says so; the method may still run, but
		// the estimate no longer covers it. A graded Gamma_h is where this is most
		// likely to fail, so it is checked here rather than assumed from the
		// uniform-mesh study.
		BOOST_TEST( history[ c ].widened == 0,
		            "cycle " << c << ": " << history[ c ].widened
		            << " vertices of Gamma_h needed a widened fan, so assumption P.1 "
		            "has been given up on a graded boundary" );
	}

	// And the loop is adaptive, not uniform.
	for ( std::size_t c = 0; c + 1 < history.size(); ++c )
		BOOST_TEST( history[ c ].marked < history[ c ].elements,
		            "cycle " << c << " marked all " << history[ c ].elements
		            << " elements, which is uniform refinement" );
}

/// The two marking strategies pull in opposite directions as gamma grows, which
/// is section 3.2's own description and is easy to get backwards. Asserted on the
/// real indicator rather than a made-up one, because the statement is about how
/// the marked fraction responds to gamma on a realistic distribution of eta_K.
BOOST_AUTO_TEST_CASE( theMarkingParametersActInOppositeDirections )
{
	meq::analytic::SolovievEquilibrium const &eq = equilibrium();
	meq::SolovievSource const source( eq.getA() );

	mfem::FunctionCoefficient sourceCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.f( x( 0 ), x( 1 ), 0.0 );
	} );
	mfem::FunctionCoefficient psiCoeff( [ &eq ]( mfem::Vector const &x )
	{
		return eq.psi( x( 0 ), x( 1 ) );
	} );

	mfem::Mesh mesh = makeFittedMesh( 8 );
	meq::GradShafranovSolver solver( mesh, 2 );
	solver.setSource( sourceCoeff );
	solver.setBoundaryData( psiCoeff );
	solver.solve();
	solver.postProcess();

	meq::ResidualEstimator estimator( solver, source );
	mfem::Vector const &local = estimator.GetLocalErrors();

	mfem::Array<int> few;
	mfem::Array<int> many;

	// Doerfler: small gamma marks few, large gamma marks many.
	meq::markDoerfler( local, 0.2, few );
	meq::markDoerfler( local, 0.9, many );
	std::printf( "\n  on %d elements: Doerfler 0.2 marks %d, 0.9 marks %d;",
	             mesh.GetNE(), few.Size(), many.Size() );
	BOOST_TEST( few.Size() < many.Size(),
	            "Doerfler marked " << few.Size() << " at gamma = 0.2 and " << many.Size()
	            << " at 0.9, so the parameter is acting backwards" );

	// Maximum: the reverse.
	meq::markMaximum( local, 0.8, few );
	meq::markMaximum( local, 0.1, many );
	std::printf( " maximum 0.8 marks %d, 0.1 marks %d\n", few.Size(), many.Size() );
	std::fflush( stdout );
	BOOST_TEST( few.Size() < many.Size(),
	            "maximum marking marked " << few.Size() << " at gamma = 0.8 and "
	            << many.Size() << " at 0.1, so the parameter is acting backwards" );
}
