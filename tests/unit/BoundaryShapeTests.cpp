#define BOOST_TEST_MODULE MeqBoundaryShapeTests

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <vector>

#include "meq/BoundaryShape.hpp"

#include "analytic/MillerDShape.hpp"

namespace
{
	double const pi = M_PI;

	/// An ITER-like D-shape: the parameters HDG-GS-1 Example 6 uses, in this
	/// file's convention. eps = 0.32 there is r/R0, so r = 0.32 with R0 = 1.
	meq::BoundaryShape iterLike()
	{
		return meq::BoundaryShape::miller( 1.0, 0.0, 0.32, 1.7, 0.33 );
	}
}

BOOST_AUTO_TEST_SUITE( boundary_shape_tests )

/*
 * MILLER IS MXH TRUNCATED, WHICH IS THE REASON THERE IS ONE EVALUATOR.
 *
 * refs/MXH.pdf eq (4) says Turnbull-Miller is recovered by keeping s_1 and s_2
 * alone, with s_1 = arcsin( delta ). If that is true then miller() and the
 * general constructor given { arcsin( delta ) } must trace the same curve, and
 * if it is not then this file has two parametrisations pretending to be one.
 */
BOOST_AUTO_TEST_CASE( millerIsMxhTruncatedAtOneHarmonic )
{
	double const delta = 0.33;
	meq::BoundaryShape const sugar = iterLike();
	meq::BoundaryShape const general( 1.0, 0.0, 0.32, 1.7, {},
	                                  { std::asin( delta ) } );

	double worst = 0.0;
	for ( int i = 0; i < 720; ++i )
	{
		double const theta = 2.0*pi*i/720.0;
		double ar = 0.0, az = 0.0, br = 0.0, bz = 0.0;
		sugar.point( theta, ar, az );
		general.point( theta, br, bz );
		worst = std::max( worst, std::hypot( ar - br, az - bz ) );
	}
	BOOST_TEST( worst < 1.0e-15 );
}

/*
 * THE PARAMETERS MEAN WHAT THEY SAY.
 *
 * A shape class that silently redefines elongation is worse than no shape class,
 * because every downstream number stays plausible. So the geometric definitions
 * are checked against the curve rather than trusted: kappa is the height-to-width
 * ratio of the bounding box, and delta places the top of the surface inboard of
 * the centre by delta times the minor radius.
 *
 * MXH section 3 gives exactly these as the way to recover { R0, Z0, kappa, r }
 * from a flux surface, so this is the paper's own inverse applied to its forward
 * map.
 */
BOOST_AUTO_TEST_CASE( theBoundingBoxRecoversTheParameters )
{
	meq::BoundaryShape const shape = iterLike();

	double rMin = 0.0, rMax = 0.0, zMin = 0.0, zMax = 0.0;
	shape.boundingBox( rMin, rMax, zMin, zMax );

	// MXH section 3: 2r = max R - min R, 2 kappa r = max Z - min Z,
	// 2 R0 = max R + min R, 2 Z0 = max Z + min Z.
	BOOST_TEST( 0.5*( rMax - rMin ) == 0.32, boost::test_tools::tolerance( 1.0e-6 ) );
	BOOST_TEST( 0.5*( rMax + rMin ) == 1.0, boost::test_tools::tolerance( 1.0e-6 ) );
	BOOST_TEST( 0.5*( zMax - zMin )/0.32 == 1.7, boost::test_tools::tolerance( 1.0e-9 ) );
	// Absolute, not relative: the centre height is zero here and a relative
	// tolerance against zero tests nothing.
	BOOST_TEST( std::abs( 0.5*( zMax + zMin ) ) < 1.0e-9 );
}

BOOST_AUTO_TEST_CASE( triangularityPlacesTheTopOfTheSurface )
{
	double const delta = 0.33;
	double const minor = 0.32;
	meq::BoundaryShape const shape = iterLike();

	// The top of the surface is theta = pi/2, where Z is largest. There
	// R = R0 + r cos( pi/2 + arcsin delta ) = R0 - r sin( arcsin delta )
	//   = R0 - r delta,
	// which is the standard definition of triangularity: the top sits inboard of
	// the centre by delta minor radii. That identity is why s_1 is arcsin( delta )
	// and not delta.
	double r = 0.0, z = 0.0;
	shape.point( 0.5*pi, r, z );

	BOOST_TEST( r == 1.0 - minor*delta, boost::test_tools::tolerance( 1.0e-12 ) );
	BOOST_TEST( z == 1.7*minor, boost::test_tools::tolerance( 1.0e-12 ) );

	// And the delta-versus-arcsin-delta confusion, asserted so a silent revert
	// fails: had s_1 been delta itself, the top would sit at
	// R0 - r sin( delta ) = R0 - 0.3239 r, which differs in the third figure.
	double const wrong = 1.0 - minor*std::sin( delta );
	BOOST_TEST( std::abs( r - wrong ) > 1.0e-4,
	            "the arcsin and the raw delta give the same point, so this test "
	            "cannot tell the two conventions apart" );
}

/*
 * THE LEVEL SET IS THE POINT OF THE CLASS, so it is checked three ways: it
 * vanishes on the curve, it has the right sign either side, and it is
 * consistent with the parametrisation it came from.
 */
BOOST_AUTO_TEST_CASE( theLevelSetVanishesOnTheCurve )
{
	meq::BoundaryShape const shape = iterLike();

	double worst = 0.0;
	for ( int i = 0; i < 360; ++i )
	{
		double const theta = 2.0*pi*i/360.0;
		double r = 0.0, z = 0.0;
		shape.point( theta, r, z );
		worst = std::max( worst, std::abs( shape.levelSet( r, z ) ) );
	}
	BOOST_TEST( worst < 1.0e-12 );
}

BOOST_AUTO_TEST_CASE( theLevelSetIsNegativeInsideAndPositiveOutside )
{
	meq::BoundaryShape const shape = iterLike();

	BOOST_TEST( shape.levelSet( 1.0, 0.0 ) < 0.0 );      // the centre
	BOOST_TEST( shape.levelSet( 1.0, 0.3 ) < 0.0 );      // inside, above it
	BOOST_TEST( shape.levelSet( 1.5, 0.0 ) > 0.0 );      // outboard
	BOOST_TEST( shape.levelSet( 1.0, 1.0 ) > 0.0 );      // above the top

	// Monotone along a ray, which the radial-gap construction guarantees and a
	// signed distance would not: the gap grows exactly as fast as the radius.
	double previous = -1.0e300;
	for ( double radius = 0.05; radius < 0.9; radius += 0.05 )
	{
		double const value = shape.levelSet( 1.0 + radius, 0.0 );
		BOOST_TEST( value > previous );
		previous = value;
	}
}

/*
 * AND THE VALIDATION THAT MATTERS, because it is the one a user will hit.
 *
 * levelSet() inverts the polar angle by bisection, which is single valued only
 * when the curve is star shaped about its centre. A large enough harmonic folds
 * it, bisection then has several roots to choose between, and it would pick one
 * and produce a domain that is quietly the wrong shape. So the constructor
 * refuses, and this checks that it does -- and that the refusal is not so eager
 * that ordinary shapes trip it.
 */
/*
 * THE LIBRARY AGAINST THE TEST FIXTURE, which are independent implementations.
 *
 * tests/analytic/MillerDShape.hpp was written for HDG-GS-1 Example 6 and carries
 * BOTH readings of the triangularity: boundaryPoint() uses that paper's printed
 * arcsin( delta sin t ), and boundaryPointCerfonFreidberg() the
 * arcsin( delta ) sin t that Miller eq (34) and Cerfon & Freidberg eq (9) print.
 * This class implements the second, so it must match the second to round-off and
 * the first only to the 6e-4 the fixture's own header quotes.
 *
 * Checking both is the point. Matching the CF form says the parametrisation is
 * right; DIFFERING from Example 6's form by the documented amount says the two
 * conventions have not been quietly conflated somewhere -- which, if they had
 * been, would leave every number here plausible.
 */
BOOST_AUTO_TEST_CASE( theShapeMatchesTheAnalyticFixture )
{
	double const eps = 0.32, delta = 0.33, kappa = 1.7;
	meq::analytic::MillerDShape const fixture( eps, delta, kappa );
	// The fixture works at R0 = 1 with eps as the minor radius, which is this
	// class's r when R0 = 1.
	meq::BoundaryShape const shape =
		meq::BoundaryShape::miller( 1.0, 0.0, eps, kappa, delta );

	double worstAgainstCf = 0.0;
	double worstAgainstExample6 = 0.0;

	for ( int i = 0; i < 720; ++i )
	{
		double const t = 2.0*pi*i/720.0;

		double fr = 0.0, fz = 0.0;
		fixture.boundaryPointCerfonFreidberg( t, fr, fz );
		double sr = 0.0, sz = 0.0;
		shape.point( t, sr, sz );
		worstAgainstCf = std::max( worstAgainstCf, std::hypot( sr - fr, sz - fz ) );

		double er = 0.0, ez = 0.0;
		fixture.boundaryPoint( t, er, ez );
		worstAgainstExample6 = std::max( worstAgainstExample6,
		                                 std::hypot( sr - er, sz - ez ) );
	}

	BOOST_TEST_MESSAGE( "  against Cerfon-Freidberg " << worstAgainstCf
	                    << ", against Example 6 " << worstAgainstExample6 );

	BOOST_TEST( worstAgainstCf < 1.0e-14,
	            "this class and the fixture's Cerfon-Freidberg form should be the "
	            "same curve, and differ by " << worstAgainstCf );

	// The fixture's header quotes at most 6e-4 between the two conventions, on a
	// minor radius of 0.32. Bounded above so a conflation is caught, and below so
	// that this test fails if the two forms ever silently become one.
	BOOST_TEST( worstAgainstExample6 < 1.0e-3 );
	BOOST_TEST( worstAgainstExample6 > 1.0e-5,
	            "Example 6's arcsin( delta sin t ) and Cerfon-Freidberg's "
	            "arcsin( delta ) sin t now agree to " << worstAgainstExample6
	            << ", so one of them has been changed into the other" );
}

BOOST_AUTO_TEST_CASE( aFoldedSurfaceIsRefused )
{
	// c_1 = 2.0 is far past the fold: theta_R doubles back on itself.
	BOOST_CHECK_THROW( meq::BoundaryShape( 1.0, 0.0, 0.32, 1.7, { 0.0, 2.0 }, {} ),
	                   meq::ShapeError );

	// And an extreme triangularity, which folds through the same mechanism.
	BOOST_CHECK_THROW( meq::BoundaryShape( 1.0, 0.0, 0.32, 1.7, {}, { 2.0 } ),
	                   meq::ShapeError );
}

BOOST_AUTO_TEST_CASE( ordinaryShapesAreAccepted )
{
	// Strongly shaped but still single valued: DIII-D-like, with the harmonics
	// MXH figure 2 plots for r/a = 0.95.
	BOOST_CHECK_NO_THROW(
		meq::BoundaryShape( 1.7, 0.0, 0.6, 1.8, { -0.1, 0.02, 0.01 },
		                    { 0.35, 0.08, -0.01 } ) );

	// A tilted surface, c_0 non-zero, which is the parameter most likely to be
	// mistaken for something harmless.
	BOOST_CHECK_NO_THROW(
		meq::BoundaryShape( 1.0, 0.0, 0.32, 1.7, { 0.15 }, { 0.34 } ) );
}

BOOST_AUTO_TEST_CASE( aSurfaceReachingTheAxisIsRefused )
{
	// r0 = 0.3 with minor 0.32 puts the inboard edge at r < 0, where the
	// operator's 1/r is not integrable. A shape class that let this through
	// would hand the solver an unsolvable problem with no diagnostic.
	BOOST_CHECK_THROW( meq::BoundaryShape( 0.3, 0.0, 0.32, 1.7 ),
	                   meq::ShapeError );
}

BOOST_AUTO_TEST_CASE( degenerateParametersAreRefused )
{
	BOOST_CHECK_THROW( meq::BoundaryShape( 1.0, 0.0, 0.0, 1.7 ), meq::ShapeError );
	BOOST_CHECK_THROW( meq::BoundaryShape( 1.0, 0.0, -0.3, 1.7 ), meq::ShapeError );
	BOOST_CHECK_THROW( meq::BoundaryShape( 1.0, 0.0, 0.32, 0.0 ), meq::ShapeError );
	// |delta| >= 1 has no arcsin.
	BOOST_CHECK_THROW( meq::BoundaryShape::miller( 1.0, 0.0, 0.32, 1.7, 1.0 ),
	                   meq::ShapeError );
}

BOOST_AUTO_TEST_SUITE_END()
