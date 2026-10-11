/*
 * A MIRROR'S PLASMA SUPPORT, ON THE SOLVER'S SIDE -- MIRROR-PLAN.md MR-1.
 *
 * meq::PlasmaSupportShape::Box makes the support `0 < Psi <= 1` and
 * `|z| <= z_edge`, and RotatingSourceTests asserts that the source vanishes on
 * all three sides of it. What it cannot assert is the SOLVER's half, which is
 * two statements and both need a mesh:
 *
 *   XP-1's ELEMENT FILL IS OFF ON THIS PATH, AS A STATED BRANCH. The fill seeds
 *   at the magnetic axis and keeps the connected component around it; a mirror
 *   has no axis, and its annulus is connected through the tube to the ends of
 *   the machine. MIRROR-PLAN.md §3 takes the decision that the box makes the
 *   fill unnecessary -- and CLAUDE_FB.md records this fill once being silently
 *   off as a real defect, found only by a 2 x 2 cross, so "off" has to be
 *   asserted against a control in which it is ON on the same field.
 *
 *   THE CUTOFF IS REPORTED. `z_edge` is prescribed, so a configuration whose
 *   plasma still carries current at the throat is cut through it, and
 *   GradShafranovSolver::checkAxialCutoff() is what makes that visible.
 *
 * THE FIRST TWO CASES RUN NO SOLVE. The field is a closed form projected onto
 * the potential space, which the diagnostic and the fill then read as they
 * would read a converged answer -- so every number they check is checked
 * against arithmetic on the closed form and not against another run of the
 * same code. The third SOLVES, because what it asks about -- the solver's own
 * labels, and its axis diagnostic -- exists only after a solve.
 *
 * THE FIELD is `psi = R^2 + z^2` on R in [ 0.2, 1.4 ], z in [ -1.5, 1.5 ],
 * labelled as a mirror: psi_ax = 0.25 is the INNER surface and psi_bnd = 1 the
 * OUTER, so the span is negative. Its surfaces close on the axis at
 * `|z| = sqrt( psi )`, so the outer one has left the domain by `|z| = 1` --
 * which is what lets one field carry a cutoff through current and a cutoff
 * past the plasma's own end. Quadratic, so the order-2 potential space holds
 * it exactly.
 */

#define BOOST_TEST_MODULE MirrorSupport
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"

using meq::GradShafranovSolver;

namespace
{
	double const innerLabel = 0.25;   // psi_ax
	double const outerLabel = 1.0;    // psi_bnd
	double const minRadius = 0.2;
	double const maxRadius = 1.4;
	double const halfHeight = 1.5;

	/// p'( Psi ) = 4 Psi ( 1 - Psi ): zero at both labels, so the box's two
	/// flux edges are where the source vanishes anyway and the AXIAL edge is
	/// the one that can cut through current. Exact at every derivative level.
	class BumpProfile : public meq::Profile
	{
		public:
			double operator()( double psi ) const override
			{
				return 4.0*psi*( 1.0 - psi );
			}
			double prime( double psi ) const override
			{
				return 4.0 - 8.0*psi;
			}
			double doublePrime( double ) const override
			{
				return -8.0;
			}
	};

	double exactPsi( double radius, double z )
	{
		return radius*radius + z*z;
	}

	/// `|F|/R`, from the closed form and nothing of MEQ's: mu0 = 1, gg' = 0,
	/// so `F = R^2 p'( Psi )/span` inside the box and zero outside it.
	double exactDensity( double radius, double z, double zEdge )
	{
		double const span = innerLabel - outerLabel;
		double const normalised = ( exactPsi( radius, z ) - outerLabel )/span;
		if ( std::abs( z ) > zEdge || !( normalised > 0.0 ) || normalised > 1.0 )
			return 0.0;
		return std::abs( radius*radius*4.0*normalised*( 1.0 - normalised )/span )
		       /radius;
	}

	struct Fixture
	{
		mfem::Mesh mesh;
		meq::NormalisedMHDSource source;
		mfem::ConstantCoefficient zero{ 0.0 };
		std::unique_ptr<GradShafranovSolver> solver;

		Fixture()
			: mesh( mfem::Mesh::MakeCartesian2D( 30, 75, mfem::Element::TRIANGLE,
			                                     false, maxRadius - minRadius,
			                                     2.0*halfHeight ) ),
			  source( std::make_shared<BumpProfile const>(),
			          std::make_shared<meq::ConstantProfile const>( 0.0 ),
			          innerLabel, 1.0 )
		{
			for ( int v = 0; v < mesh.GetNV(); ++v )
			{
				mesh.GetVertex( v )[ 0 ] += minRadius;
				mesh.GetVertex( v )[ 1 ] -= halfHeight;
			}

			solver = std::make_unique<GradShafranovSolver>( mesh, 2 );
			solver->setSource( source, innerLabel );
			solver->setBoundaryData( zero );
			solver->prepare();

			// The mirror's labels, and the field, AFTER prepare(): nothing
			// solves here, so nothing will overwrite either.
			source.setNormalisation( innerLabel, outerLabel );
			mfem::FunctionCoefficient field( [ ]( mfem::Vector const &x )
			{
				return exactPsi( x( 0 ), x( 1 ) );
			} );
			solver->potential().ProjectCoefficient( field );
		}
	};
}

/*
 * THE FILL IS OFF ON THE BOX AND ON FOR THE TOKAMAK'S SHAPE, ON ONE FIELD.
 *
 * The one-sided arm is the CONTROL and comes first: under it the fill must
 * run and reach elements, or "zero under the box" would be the zero of a fill
 * that never runs here at all. The box arm then refreshes the SAME solver at
 * the SAME field, which also asserts that a mask filled before the shape
 * changed does not survive into the box -- holds() keeps restricting F to a
 * stale component otherwise.
 */
BOOST_AUTO_TEST_CASE( theElementFillIsOffOnAMirrorsBox )
{
	Fixture m;
	GradShafranovSolver &solver = *m.solver;

	BOOST_TEST_REQUIRE( ( solver.plasmaConnectivity()
	                      == GradShafranovSolver::PlasmaConnectivity::Component ),
	                    "the fill is the default and this case is about the "
	                    "default" );

	m.source.setPlasmaSupport( true );
	BOOST_TEST_REQUIRE( ( m.source.plasmaSupportShape()
	                      == meq::PlasmaSupportShape::OneSided ) );
	solver.refreshPlasmaComponent( solver.potential() );
	int const oneSidedReached = solver.plasmaComponentElements();
	int const oneSidedCandidates = solver.plasmaCandidateElements();

	m.source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, 0.5 );
	solver.refreshPlasmaComponent( solver.potential() );
	int const boxReached = solver.plasmaComponentElements();
	int const boxCandidates = solver.plasmaCandidateElements();

	int masked = 0;
	for ( int e = 0; e < m.mesh.GetNE(); ++e )
		if ( !solver.elementInPlasma( e ) )
			++masked;

	std::printf( "\n  XP-1's FILL ON A MIRROR-LABELLED FIELD, %d elements\n",
	             m.mesh.GetNE() );
	std::printf( "    %-12s %10s %10s\n", "shape", "reached", "candidates" );
	std::printf( "    %-12s %10d %10d\n", "one-sided", oneSidedReached,
	             oneSidedCandidates );
	std::printf( "    %-12s %10d %10d   elements masked off: %d\n", "box",
	             boxReached, boxCandidates, masked );
	std::fflush( stdout );

	BOOST_TEST_REQUIRE( oneSidedReached > 0,
	                    "the CONTROL: under the tokamak's one-sided support the "
	                    "fill must run on this field, or the box's zero below "
	                    "says nothing" );
	BOOST_TEST( boxReached == 0,
	            "the fill reached " << boxReached << " elements under "
	            "PlasmaSupportShape::Box. A mirror has no axis to seed it and "
	            "its annulus is connected through the tube, so "
	            "GradShafranovSolver::plasmaComponentWanted() must decline it" );
	BOOST_TEST( boxCandidates == 0,
	            "a mask is still active under the box -- the one-sided fill's, "
	            "left standing" );
	BOOST_TEST( masked == 0,
	            masked << " elements are masked off under the box. The support "
	            "is pointwise there and no element may be switched off whole" );
}

/*
 * THE CUTOFF DIAGNOSTIC READS THE CURRENT THE BOX CUTS THROUGH, AND READS ZERO
 * WHERE THE PLASMA HAS ENDED FIRST.
 *
 * Three cutoffs on one field, each a different statement:
 *
 *   z_edge = 0.5   the surfaces are still open at the cutoff, so it cuts
 *                  through current. The reading must be the closed form's own
 *                  maximum of |F|/R on |z| = 0.5, and so must the peak over
 *                  the box.
 *   z_edge = 1.2   the outer surface has met the axis by |z| = 1, so Psi < 0
 *                  everywhere on the cutoff: the reading is EXACTLY zero, from
 *                  samples that WERE located -- a zero with nothing behind it
 *                  would be the next case.
 *   z_edge = 2.0   beyond the mesh: no sample is located, and the report says
 *                  so rather than reading zero as "no current".
 *
 * and the tokamak's shape, where it does not apply at all.
 *
 * The closed-form maxima are taken over a sweep far finer than the sampler's,
 * so the solver's readings are bounded ABOVE by them to within that sweep's
 * own resolution -- every value it reads is an exact sample of the closed
 * form, the field being held exactly -- and below by how close 513 samples
 * and the nodes come to the maximum.
 */
BOOST_AUTO_TEST_CASE( theCutoffDiagnosticReadsTheCurrentTheBoxCutsThrough )
{
	Fixture m;
	GradShafranovSolver &solver = *m.solver;

	auto lineMaximum = [ ]( double z, double zEdge )
	{
		double best = 0.0;
		for ( int i = 0; i <= 200000; ++i )
		{
			double const radius = minRadius + ( maxRadius - minRadius )*i/200000.0;
			best = std::max( best, exactDensity( radius, z, zEdge ) );
		}
		return best;
	};

	auto boxMaximum = [ & ]( double zEdge )
	{
		double best = 0.0;
		for ( int j = 0; j <= 600; ++j )
		{
			double const z = -zEdge + 2.0*zEdge*j/600.0;
			for ( int i = 0; i <= 2000; ++i )
			{
				double const radius = minRadius + ( maxRadius - minRadius )*i/2000.0;
				best = std::max( best, exactDensity( radius, z, zEdge ) );
			}
		}
		return best;
	};

	// THE TOKAMAK'S SHAPE: nothing to report.
	m.source.setPlasmaSupport( true );
	GradShafranovSolver::AxialCutoffCheck const tokamak = solver.checkAxialCutoff();
	BOOST_TEST( tokamak.applies == false );
	BOOST_TEST( tokamak.atCutoff == 0.0 );

	std::printf( "\n  THE AXIAL CUTOFF ON psi = R^2 + z^2, labelled 0.25 inner / "
	             "1.0 outer\n" );
	std::printf( "    %7s %12s %12s %12s %12s %10s %8s\n", "z_edge", "peak",
	             "closed form", "at cutoff", "closed form", "relative",
	             "samples" );

	// THROUGH CURRENT.
	{
		double const zEdge = 0.5;
		m.source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, zEdge );
		GradShafranovSolver::AxialCutoffCheck const r = solver.checkAxialCutoff();
		double const wantPeak = boxMaximum( zEdge );
		double const wantCutoff = lineMaximum( zEdge, zEdge );

		std::printf( "    %7.2f %12.6e %12.6e %12.6e %12.6e %10.4f %8d\n", zEdge,
		             r.peak, wantPeak, r.atCutoff, wantCutoff, r.relative,
		             r.cutoffSamples );

		BOOST_TEST_REQUIRE( r.applies );
		BOOST_TEST( r.axialEdge == zEdge );
		BOOST_TEST_REQUIRE( wantCutoff > 0.0 );
		BOOST_TEST( r.cutoffSamples > 0 );

		BOOST_TEST( r.atCutoff <= wantCutoff*( 1.0 + 1.0e-6 ),
		            "the cutoff reads " << r.atCutoff << ", ABOVE the closed "
		            "form's maximum " << wantCutoff << " on |z| = " << zEdge );
		BOOST_TEST( r.atCutoff >= wantCutoff*( 1.0 - 1.0e-4 ),
		            "the cutoff reads " << r.atCutoff << " against the closed "
		            "form's maximum " << wantCutoff << " on |z| = " << zEdge );
		BOOST_TEST( std::abs( std::abs( r.cutoffZ ) - zEdge ) == 0.0,
		            "the cutoff value must be read ON the cutoff, at |z| = "
		            << zEdge << ", not at " << r.cutoffZ );

		BOOST_TEST( r.peak <= wantPeak*( 1.0 + 1.0e-6 ) );
		BOOST_TEST( r.peak >= wantPeak*( 1.0 - 1.0e-2 ),
		            "the peak reads " << r.peak << " against the closed form's "
		            << wantPeak );
		BOOST_TEST( std::abs( r.peakZ ) <= zEdge );
		BOOST_TEST( r.relative == r.atCutoff/r.peak );
		BOOST_TEST( r.relative > 0.1,
		            "a cutoff through open surfaces must read a current of the "
		            "plasma's own order" );
	}

	// PAST THE PLASMA'S OWN END.
	{
		double const zEdge = 1.2;
		m.source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, zEdge );
		GradShafranovSolver::AxialCutoffCheck const r = solver.checkAxialCutoff();

		std::printf( "    %7.2f %12.6e %12.6e %12.6e %12.6e %10.4f %8d\n", zEdge,
		             r.peak, boxMaximum( zEdge ), r.atCutoff,
		             lineMaximum( zEdge, zEdge ), r.relative, r.cutoffSamples );

		BOOST_TEST_REQUIRE( r.applies );
		BOOST_TEST( r.cutoffSamples > 0,
		            "the cutoff is inside the mesh, so samples must be located "
		            "-- a zero below would otherwise be a zero with nothing "
		            "behind it" );
		BOOST_TEST( r.atCutoff == 0.0 );
		BOOST_TEST( r.relative == 0.0 );
		BOOST_TEST( r.peak > 0.0 );
	}

	// BEYOND THE MESH.
	{
		double const zEdge = 2.0;
		m.source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, zEdge );
		GradShafranovSolver::AxialCutoffCheck const r = solver.checkAxialCutoff();

		std::printf( "    %7.2f %12.6e %12s %12.6e %12s %10.4f %8d\n", zEdge,
		             r.peak, "", r.atCutoff, "", r.relative, r.cutoffSamples );

		BOOST_TEST_REQUIRE( r.applies );
		BOOST_TEST( r.cutoffSamples == 0,
		            "the cutoff lies outside the mesh, so nothing on it can "
		            "have been located" );
		BOOST_TEST( r.atCutoff == 0.0 );
	}
	std::fflush( stdout );

	BOOST_CHECK_THROW( solver.checkAxialCutoff( 1 ), std::invalid_argument );
}

/*
 * A MIRROR'S SYMMETRY AXIS IS OUTSIDE ITS BOX, AND checkAxisSource() SAYS SO.
 *
 * psi( 0, z ) = 0 exactly, so Psi on the axis is -psi_bnd/span, and on a
 * mirror -- inner label the smaller flux, span negative -- that is ABOVE 1.
 * The one-sided rule reads only its sign and puts the axis inside; the box's
 * Psi <= 1 excludes it, and the gate is the box. So the diagnostic must ask the
 * gate's own rule, NormalisedSource::supportAdmits(), or it reports the axis
 * inside a plasma no source evaluation will ever place there -- harmless at
 * g g' = 0, where F vanishes on the axis anyway, and a false refusal of the
 * driver's at g g' != 0.
 *
 * A REAL SOLVE, because the labels checkAxisSource() reads are the solver's
 * converged ones. The smallest mirror there is: psi = R^2 on the boundary of
 * [ 0, 1.4 ] x [ -1.5, 1.5 ] -- Delta*-harmonic, so it is the vacuum field --
 * with both labels pinned at prescribed points on the midplane, the box
 * support, and BumpProfile, which vanishes at both labels as the box's
 * precondition asks. mu0 = 0.1 keeps the plasma a perturbation of R^2, so the
 * labels sit near 0.25 and 1 and Psi on the axis near 4/3.
 *
 * THE TOKAMAK'S READING IS THE CONTROL: the same solver, the same labels, the
 * shape put back to OneSided, must still read the axis INSIDE -- which is the
 * one-sided rule unchanged, and what makes the box's "outside" a statement
 * about the box rather than about the field.
 *
 * AND plasmaSupportNodes(), the count the bordered solve's failure report
 * reads in place of the span's sign: non-empty here with a negative span, and
 * empty at a frozen edge no node reaches.
 */
BOOST_AUTO_TEST_CASE( aMirrorsAxisIsOutsideItsBoxAndTheDiagnosticAgrees )
{
	double const width = 1.4;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 28, 60, mfem::Element::TRIANGLE,
	                                               false, width, 2.0*halfHeight );
	for ( int v = 0; v < mesh.GetNV(); ++v )
		mesh.GetVertex( v )[ 1 ] -= halfHeight;

	meq::NormalisedMHDSource source( std::make_shared<BumpProfile const>(),
	                                 std::make_shared<meq::ConstantProfile const>( 0.0 ),
	                                 innerLabel, 0.1 );
	source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, 1.0 );

	mfem::FunctionCoefficient vacuum( [ ]( mfem::Vector const &x )
	{
		return x( 0 )*x( 0 );
	} );

	GradShafranovSolver solver( mesh, 2 );
	solver.setAxisFluxPoint( std::sqrt( innerLabel ), 0.0 );
	solver.setBoundaryFluxPoint( std::sqrt( outerLabel ), 0.0 );
	// MANDATORY on a mirror: the default 0 makes the first span positive and
	// puts the box around the axis.
	solver.setBoundaryFluxInitial( outerLabel );
	solver.setSource( source, innerLabel );
	solver.setBoundaryData( vacuum );
	solver.setInitialGuess( vacuum );
	solver.setNewtonControl( 1.0e-12, 1.0e-14, 40 );
	solver.solve();

	double const span = solver.psiAxis() - solver.psiBoundary();
	GradShafranovSolver::AxisSourceCheck const box = solver.checkAxisSource();
	GradShafranovSolver::PlasmaSupportNodes const counted = solver.plasmaSupportNodes();

	std::printf( "\n  A MIRROR'S AXIS AGAINST ITS BOX: psi_ax %.10f, psi_bnd %.10f, "
	             "span %.4e, %d Newton steps\n", solver.psiAxis(),
	             solver.psiBoundary(), span, solver.newtonIterations() );
	std::printf( "    Psi on R = 0 %.6f; box: axis inside %d; support holds %d of "
	             "%d nodes\n", box.normalisedFluxOnAxis, box.axisInsidePlasma,
	             counted.inside, counted.total );

	BOOST_TEST_REQUIRE( span < 0.0,
	                    "the fixture is a mirror only if the inner label is the "
	                    "smaller flux; span " << span );
	BOOST_TEST_REQUIRE( box.reachesAxis );
	BOOST_TEST_REQUIRE( box.normalisedFluxOnAxis > 1.0,
	                    "Psi on the axis is " << box.normalisedFluxOnAxis
	                    << ", so the box's upper bound is not what this case "
	                    "exercises" );
	BOOST_TEST( !box.axisInsidePlasma,
	            "checkAxisSource() puts a mirror's axis inside its plasma at "
	            "Psi = " << box.normalisedFluxOnAxis << ", which the box's "
	            "Psi <= 1 excludes. It must ask the gate's own rule, "
	            "NormalisedSource::supportAdmits(), not the one-sided sign" );
	BOOST_TEST( box.sourceVanishesOnAxis );

	BOOST_TEST( counted.inside > 0,
	            "a converged mirror's support holds no node: the count the "
	            "border's failure report reads would call it EMPTY" );
	BOOST_TEST( counted.inside < counted.total );

	// THE CONTROL: the tokamak's one-sided rule at the same labels.
	source.setPlasmaSupportShape( meq::PlasmaSupportShape::OneSided );
	GradShafranovSolver::AxisSourceCheck const oneSided = solver.checkAxisSource();
	std::printf( "    one-sided at the same labels: axis inside %d\n",
	             oneSided.axisInsidePlasma );
	BOOST_TEST( oneSided.normalisedFluxOnAxis == box.normalisedFluxOnAxis );
	BOOST_TEST( oneSided.axisInsidePlasma,
	            "the one-sided rule reads Psi = " << oneSided.normalisedFluxOnAxis
	            << " > 0 as inside, and must go on doing so: that is the "
	            "tokamak's reading, which the driver refuses on" );

	// AND AN EMPTY SUPPORT IS COUNTED AS ONE: a frozen edge above every psi.
	source.setPlasmaSupportShape( meq::PlasmaSupportShape::Box, 1.0 );
	source.freezePlasmaEdge( 10.0, 9.0 );
	GradShafranovSolver::PlasmaSupportNodes const empty = solver.plasmaSupportNodes();
	source.thawPlasmaEdge();
	std::printf( "    frozen at ( 10, 9 ): support holds %d of %d nodes\n",
	             empty.inside, empty.total );
	std::fflush( stdout );
	BOOST_TEST( empty.inside == 0 );
	BOOST_TEST( empty.total == counted.total );
}
