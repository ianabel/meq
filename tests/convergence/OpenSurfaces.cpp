#define BOOST_TEST_MODULE OpenSurfaces
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mfem.hpp"

#include "meq/FluxSurfaces.hpp"
#include "meq/SurfaceFit.hpp"

/*
 * OPEN SURFACES -- INVERSION-PLAN.md stage IN-5.
 *
 * Every surface IN-0 to IN-6 handles is a loop: the tracer closes it, the disc
 * chart parametrises it by an angle about the axis, and the averages integrate
 * round it. **A level outside psi_bnd is none of those things.** It runs to the
 * wall and stops, so it has two genuine endpoints, no period, and no enclosed
 * area -- and the machinery that serves a closed surface says nothing useful
 * about it.
 *
 * v0-legacy:FluxSurfaces.cpp printed "Terminating because curve left domain"
 * and returned the arc it had reached labelled as a contour, which is the
 * failure this stage exists to replace: half a curve, presented as a whole one.
 *
 *
 * THE FIXTURE IS ANALYTIC AND THE TRACER IS POINTED AT IT DIRECTLY, which is
 * what separates "the open machinery is wrong" from "the discretisation is
 * coarse". psi = z - a ( r - r0 )^2 has level sets z = c + a ( r - r0 )^2 --
 * parabolas, open across the box, and known in closed form. Interpolated into
 * an H1 space of degree 2 the quadratic is represented EXACTLY, so the traced
 * points sit on the true curve to the corrector's tolerance and every number
 * below is a property of the tracing and the fit alone.
 *
 * That is meq::ContourTracer's second constructor being used for the reason its
 * own comment gives.
 */
namespace
{
	double const a = 2.0;
	double const r0 = 1.0;
	double const rMin = 0.6, rMax = 1.4, zMin = -0.8, zMax = 0.8;
	double const level = 0.0;

	double exactZ( double r ) { return level + a*( r - r0 )*( r - r0 ); }

	/// Distance from a point to the exact parabola, to first order in the
	/// offset: the vertical gap over the slope's hypotenuse. Exact enough at
	/// the sizes below -- the fits reach 1e-9 m against a curve of length 1 m --
	/// and a genuine foot-of-perpendicular would be a root find whose own
	/// tolerance would then be in the budget.
	double distanceToExact( double r, double z )
	{
		double const gap = z - exactZ( r );
		double const slope = 2.0*a*( r - r0 );
		return std::abs( gap )/std::sqrt( 1.0 + slope*slope );
	}

	struct Field
	{
		std::unique_ptr<mfem::Mesh> mesh;
		std::unique_ptr<mfem::H1_FECollection> scalar;
		std::unique_ptr<mfem::H1_FECollection> vector;
		std::unique_ptr<mfem::FiniteElementSpace> scalarSpace;
		std::unique_ptr<mfem::FiniteElementSpace> vectorSpace;
		std::unique_ptr<mfem::GridFunction> potential;
		std::unique_ptr<mfem::GridFunction> flux;
	};

	/// psi and grad psi on an n x n mesh. The tracer is handed grad psi rather
	/// than q, which its own header says "changes the corrector's scaling by r
	/// and nothing else" -- and here there is no r weight to respect, the field
	/// being a manufactured one rather than an equilibrium.
	Field buildField( int n )
	{
		Field out;
		out.mesh = std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian2D(
			n, n, mfem::Element::TRIANGLE, false, rMax - rMin, zMax - zMin ) );
		for ( int i = 0; i < out.mesh->GetNV(); ++i )
		{
			double *v = out.mesh->GetVertex( i );
			v[ 0 ] += rMin;
			v[ 1 ] += zMin;
		}

		out.scalar = std::make_unique<mfem::H1_FECollection>( 2, 2 );
		out.vector = std::make_unique<mfem::H1_FECollection>( 2, 2 );
		out.scalarSpace = std::make_unique<mfem::FiniteElementSpace>(
			out.mesh.get(), out.scalar.get() );
		out.vectorSpace = std::make_unique<mfem::FiniteElementSpace>(
			out.mesh.get(), out.vector.get(), 2 );

		out.potential = std::make_unique<mfem::GridFunction>( out.scalarSpace.get() );
		out.flux = std::make_unique<mfem::GridFunction>( out.vectorSpace.get() );

		mfem::FunctionCoefficient psi( []( mfem::Vector const &x )
		{
			return x( 1 ) - a*( x( 0 ) - r0 )*( x( 0 ) - r0 );
		} );
		mfem::VectorFunctionCoefficient grad( 2,
			[]( mfem::Vector const &x, mfem::Vector &value )
			{
				value( 0 ) = -2.0*a*( x( 0 ) - r0 );
				value( 1 ) = 1.0;
			} );

		out.potential->ProjectCoefficient( psi );
		out.flux->ProjectCoefficient( grad );
		return out;
	}
}

/*
 * THE TRACER RETURNS THE WHOLE CURVE, WITH BOTH ENDS ON THE WALL.
 *
 * trace() follows ONE direction from its seed, so on an open level it comes
 * back with half the curve and the status LeftMesh. traceOpen() traces the
 * other half as well -- the same march with the tangent reversed, which is
 * exact rather than approximate -- and joins them.
 */
BOOST_AUTO_TEST_CASE( anOpenLevelIsTracedFromWallToWall )
{
	Field const field = buildField( 40 );
	meq::ContourTracer tracer( *field.potential, *field.flux );

	meq::Contour const half = tracer.trace( level, r0, exactZ( r0 ) );
	meq::Contour const whole = tracer.traceOpen( level, r0, exactZ( r0 ) );

	std::printf( "\n  ONE DIRECTION AGAINST BOTH\n" );
	std::printf( "    trace()      %s, %zu points, r from %.4f to %.4f\n",
	             meq::contourStatusName( half.status ), half.points.size(),
	             half.points.front().r, half.points.back().r );
	std::printf( "    traceOpen()  %s, %zu points, r from %.4f to %.4f\n",
	             meq::contourStatusName( whole.status ), whole.points.size(),
	             whole.points.front().r, whole.points.back().r );

	BOOST_TEST( ( half.status == meq::ContourStatus::LeftMesh ),
		"trace() on an open level should reach the edge of the field and say "
		"so; this fixture's level does not close, but it returned "
		<< meq::contourStatusName( half.status ) );

	BOOST_TEST_REQUIRE( ( whole.status == meq::ContourStatus::Open ),
		"traceOpen() did not reach the wall at both ends: "
		<< meq::contourStatusName( whole.status ) );

	// BOTH ENDS ARE AT THE WALL, and at OPPOSITE walls -- which is the property
	// that says the two halves were traced in opposite directions rather than
	// the same one twice.
	double const lowR = std::min( whole.points.front().r, whole.points.back().r );
	double const highR = std::max( whole.points.front().r, whole.points.back().r );
	BOOST_TEST( lowR < rMin + 0.05,
		"the low end stopped at r = " << lowR << " rather than at the wall" );
	BOOST_TEST( highR > rMax - 0.05,
		"the high end stopped at r = " << highR << " rather than at the wall" );

	// THE JOIN DOES NOT DUPLICATE THE SEED and does not fold: consecutive
	// points advance monotonically in r on this fixture, which a reversed
	// second half would break at the seam while every point stayed on the
	// level set. That is the quiet failure the join has to avoid.
	BOOST_TEST( whole.points.size() > half.points.size() );
	std::size_t reversals = 0;
	for ( std::size_t i = 1; i < whole.points.size(); ++i )
		if ( ( whole.points[ i ].r - whole.points[ i - 1 ].r )
		     *( whole.points[ 1 ].r - whole.points[ 0 ].r ) <= 0.0 )
			++reversals;
	BOOST_TEST( reversals == 0u,
		"the joined curve turns back on itself " << reversals << " times, so "
		"the reversed half was appended without its orientation" );

	// ARC LENGTH IS CUMULATIVE OVER THE WHOLE CURVE, from one endpoint. Each
	// half measured its own from the shared seed, so the seam would otherwise
	// carry two origins.
	BOOST_TEST( whole.points.front().arcLength == 0.0 );
	for ( std::size_t i = 1; i < whole.points.size(); ++i )
		BOOST_TEST( whole.points[ i ].arcLength
		            > whole.points[ i - 1 ].arcLength );

	// AND EVERY POINT IS ON THE EXACT CURVE. The space represents this
	// quadratic exactly, so what is left is the corrector's tolerance and
	// nothing about the mesh.
	double worst = 0.0;
	for ( meq::ContourPoint const &p : whole.points )
		worst = std::max( worst, distanceToExact( p.r, p.z ) );

	/*
	 * THE BOUND IS THE RESIDUAL THE CORRECTOR ACHIEVED, NOT THE ONE IT WAS
	 * ASKED FOR, and the difference is a real property of the corrector rather
	 * than slack: it keeps its best iterate and accepts after four
	 * non-improving steps, so the target is a goal and worstResidual is the
	 * outcome. Asserting against the target would be asserting that the
	 * corrector always meets it, which it documents that it does not.
	 *
	 * AND THE RELATION IS EXACT ON THIS FIXTURE. Distance to a level set is
	 * | psi - level |/| grad psi | to first order, and | grad psi | >= 1
	 * everywhere here because d psi/dz is exactly 1. So the distance cannot
	 * EXCEED the residual, whatever the residual turns out to be -- which makes
	 * this a statement about the tracer with no chosen number in it. The mesh
	 * cannot contribute, the space representing this quadratic exactly.
	 */
	std::printf( "    worst distance to the exact parabola  %.3e m\n", worst );
	std::printf( "    corrector target %.3e, worst residual achieved %.3e\n",
	             whole.correctorTarget, whole.worstResidual );
	BOOST_TEST( worst <= whole.worstResidual,
		"the traced points are " << worst << " m from the exact level set "
		"while the potential residual there is only " << whole.worstResidual
		<< ". Since | grad psi | >= 1 on this fixture the distance cannot "
		"exceed the residual, so the two disagree about where the curve is" );

	// AND THE RESIDUAL ITSELF IS SMALL, which is the half the relation above
	// does not cover: it would hold just as well for a trace that wandered off
	// the level set by a millimetre.
	BOOST_TEST( whole.worstResidual < 1.0e-8,
		"the trace did not stay on the level set: worst residual "
		<< whole.worstResidual );
}

/*
 * THE CHEBYSHEV FIT CONVERGES AND A PERIODIC ONE CANNOT -- section 4.2's claim,
 * measured.
 *
 * The curve is analytic, so a Chebyshev series in normalised arc length should
 * converge GEOMETRICALLY: each added pair of modes multiplies the error by a
 * constant factor rather than adding an algebraic order.
 *
 * **THE CONTROL IS THE POINT.** A periodic basis is not merely a worse choice
 * here, it is an inadmissible one: it forces R( -1 ) = R( +1 ) on a curve whose
 * ends are 0.8 m apart, so no number of modes can represent it and the error
 * FLOORS. Given the same mode count and the same least squares -- the two
 * differ in their basis and in nothing else -- that floor is what says the
 * choice of basis was a decision rather than a preference.
 *
 * The error is measured AWAY FROM THE SAMPLES, at points the fit did not see,
 * against the closed form. A residual at its own samples is what a least
 * squares minimises and would flatter both columns.
 */
BOOST_AUTO_TEST_CASE( theChebyshevFitConvergesWhereAPeriodicOneFloors )
{
	Field const field = buildField( 60 );
	meq::ContourTracer tracer( *field.potential, *field.flux );
	meq::Contour const curve = tracer.traceOpen( level, r0, exactZ( r0 ) );
	BOOST_TEST_REQUIRE( ( curve.status == meq::ContourStatus::Open ),
		"the fixture's level is not open: "
		<< meq::contourStatusName( curve.status ) );

	std::vector<double> arc, r, z;
	for ( meq::ContourPoint const &p : curve.points )
	{
		arc.push_back( p.arcLength );
		r.push_back( p.r );
		z.push_back( p.z );
	}
	std::printf( "\n  AN OPEN SURFACE FITTED IN NORMALISED ARC LENGTH\n" );
	std::printf( "    %zu samples over %.4f m\n", arc.size(), arc.back() );
	std::printf( "    modes   Chebyshev, off-sample   periodic control\n" );

	auto offSample = []( meq::OpenSurfaceFit const &fit, bool periodic )
	{
		double worst = 0.0;
		// 401 points, deliberately not the sample count and not a divisor of
		// it, so none of them lands on a sample by construction.
		for ( int i = 0; i <= 400; ++i )
		{
			double const t = -1.0 + 2.0*static_cast<double>( i )/400.0;
			double fittedR = 0.0, fittedZ = 0.0;
			if ( periodic )
				meq::evaluateOpenSurfacePeriodic( fit, t, fittedR, fittedZ );
			else
				meq::evaluateOpenSurface( fit, t, fittedR, fittedZ );
			worst = std::max( worst, distanceToExact( fittedR, fittedZ ) );
		}
		return worst;
	};

	std::vector<std::size_t> const counts = { 4, 8, 12, 16, 20, 24, 28, 32 };
	std::vector<double> chebyshev, periodic;
	for ( std::size_t modes : counts )
	{
		meq::OpenSurfaceFit const fit = meq::fitOpenSurface( arc, r, z, modes );
		meq::OpenSurfaceFit const control =
			meq::fitOpenSurfacePeriodic( arc, r, z, modes );
		chebyshev.push_back( offSample( fit, false ) );
		periodic.push_back( offSample( control, true ) );
		std::printf( "     %2zu           %.3e            %.3e",
		             modes, chebyshev.back(), periodic.back() );
		if ( chebyshev.size() > 1 )
			std::printf( "     x%.2f",
			             chebyshev[ chebyshev.size() - 2 ]/chebyshev.back() );
		std::printf( "\n" );
	}

	/*
	 * IT CONVERGES GEOMETRICALLY, AND THE RATIO IS WHAT SAYS SO RATHER THAN THE
	 * TOTAL.
	 *
	 * The steps here are EQUAL AND ADDITIVE -- four modes each time -- so the
	 * two decays are told apart by the shape of the ratio column and not by its
	 * size. A geometric error `C q^m` gives the SAME ratio at every step. An
	 * algebraic `C m^-p` gives `( 16/12 )^p`, then `( 20/16 )^p`, then
	 * `( 24/20 )^p`, which FALLS steadily because the multiplier does. Measured,
	 * the column is flat at about 6.4.
	 */
	double const gain = chebyshev.front()/chebyshev.back();
	std::printf( "    Chebyshev improved by %.3e over %zu -> %zu modes\n",
	             gain, counts.front(), counts.back() );
	BOOST_TEST( gain > 1.0e4,
		"the Chebyshev fit improved by only " << gain << " over the sweep" );

	// THE FLAT RANGE STOPS BEFORE THE FLOOR. By 32 modes the fit has reached
	// about 6e-8 m, which is the traced points' own accuracy carried through a
	// least squares of 32 columns -- so the last ratio is the floor arriving
	// and not the decay changing, and asserting on it would be asserting on
	// round-off.
	std::vector<double> ratios;
	for ( std::size_t i = 2; i + 2 < counts.size(); ++i )
		ratios.push_back( chebyshev[ i - 1 ]/chebyshev[ i ] );
	BOOST_TEST_REQUIRE( ratios.size() >= 3u );

	double lowest = ratios.front(), highest = ratios.front();
	for ( double value : ratios )
	{
		lowest = std::min( lowest, value );
		highest = std::max( highest, value );
	}
	std::printf( "    the ratio over equal 4-mode steps runs %.2f to %.2f\n",
	             lowest, highest );
	BOOST_TEST( lowest > 3.0,
		"the fit is barely improving per mode, which is not a geometric decay" );
	BOOST_TEST( highest/lowest < 1.5,
		"the ratio per equal step ran from " << lowest << " to " << highest
		<< ", which falls the way an ALGEBRAIC decay does. A geometric one "
		"holds its ratio over equal additive steps" );

	// MONOTONE, which an algebraic basis fighting the parametrisation would not
	// be -- IN-3's prescribed-angle column decayed at L^-1.2 and wandered.
	for ( std::size_t i = 1; i < chebyshev.size(); ++i )
		BOOST_TEST( chebyshev[ i ] < chebyshev[ i - 1 ],
			"the Chebyshev error rose from " << counts[ i - 1 ] << " to "
			<< counts[ i ] << " modes" );

	/*
	 * AND THE CONTROL DOES NOT CONVERGE AT ALL -- it gets WORSE, which is
	 * sharper than merely stalling and is exactly what an inadmissible basis
	 * should do. Given more modes it fits the samples better in the
	 * least-squares sense while the periodicity it cannot escape pushes the
	 * curve further from the truth between them.
	 */
	double const controlGain = periodic.front()/periodic.back();
	std::printf( "    the periodic control improved by %.3e over the same "
	             "modes\n", controlGain );
	BOOST_TEST( controlGain < 2.0,
		"the periodic control improved by " << controlGain << ", so it is "
		"converging -- and if a periodic basis can represent this curve then "
		"section 4.2's argument for Chebyshev is not what these numbers say" );
	BOOST_TEST( periodic.back() > 1.0e5*chebyshev.back(),
		"the two bases came out comparable, so this control is not "
		"discriminating and the choice between them is not being tested" );
}
