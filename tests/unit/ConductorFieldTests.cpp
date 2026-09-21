// Unit tests for the conductor field in src/meq/ConductorField.
//
// COIL-SUBTRACTION-PLAN.md CS-1 and CS-1b. The split is psi = psi_c + psi_p,
// with psi_c the conductors' own field evaluated exactly and psi_p what MEQ
// solves for, so that the conductors need not be in the mesh at all. This file
// is the unit under that: psi_c over FILAMENTS and RECTANGLES both, and the
// refusal that decides which meshes may carry it.
//
// FOUR THINGS HERE COULD BE SILENTLY WRONG AND NO CONVERGENCE RATE WOULD SEE
// ANY OF THEM:
//
//   * SUPERPOSITION. psi_c is a sum over conductors, and the whole plan rests
//     on Delta* being linear -- so a set of N must give exactly what the N
//     single-conductor fields give added up. A factor applied per conductor
//     rather than once converges at full rate to the wrong equilibrium.
//   * WHICH KIND THE COINCIDENCE TEST IS ABOUT. A filament has a line
//     singularity and a rectangle does not -- meq::coilPsi() is "Valid
//     EVERYWHERE, including inside the coil" -- so refusing a node inside a
//     coil would reject the ordinary configuration this plan exists to make
//     cheap, while NOT refusing one on a filament is an infinity in a solve.
//   * THE AXIS. psi_c must be EXACTLY zero at r = 0, bit for bit, because that
//     is the boundary condition the free-boundary problem imposes there and the
//     split must not disturb it. An epsilon there is a boundary condition
//     quietly becoming inhomogeneous.
//   * THE COINCIDENCE TEST. It is what stands between a user's TOML and an
//     infinity, and both of its failure directions are silent: too tight and a
//     coincident node reaches filamentPsi() and throws from inside a solve; too
//     loose and it rejects meshes that are fine.
//
// MFEM-FREE, like the unit under test, so CI can run it.

#define BOOST_TEST_MODULE ConductorFieldTests
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <stdexcept>

#include "meq/ConductorField.hpp"

namespace
{
	// A three-filament set with currents of both signs and different radii, so
	// that a per-filament error cannot cancel against another.
	meq::ConductorField threeFilaments()
	{
		meq::ConductorField field;
		field.add( meq::CurrentFilament( 1.20, -0.80, +3.0e5 ) );
		field.add( meq::CurrentFilament( 2.05, +0.35, -1.7e5 ) );
		field.add( meq::CurrentFilament( 0.75, +1.10, +9.0e4 ) );
		return field;
	}
}

// SUPERPOSITION, AND IT IS ASSERTED AS AN IDENTITY RATHER THAN A TOLERANCE.
// The sum is taken in insertion order in both arms, so the two are the same
// floating-point operations in the same order and anything but bit equality is
// a real difference rather than reassociation.
BOOST_AUTO_TEST_CASE( the_field_of_a_set_is_the_sum_of_its_filaments )
{
	meq::ConductorField const field = threeFilaments();
	BOOST_TEST_REQUIRE( field.filamentCount() == 3u );
	BOOST_TEST_REQUIRE( field.size() == 3u );
	BOOST_TEST_REQUIRE( field.coilCount() == 0u );

	double const points[][ 2 ] = { { 1.70, 0.00 }, { 0.40, -1.90 },
	                               { 2.60, +0.90 }, { 1.05, +1.05 } };

	for ( auto const &p : points )
	{
		double byHand = 0.0;
		for ( std::size_t i = 0; i < field.filamentCount(); ++i )
			byHand += meq::filamentPsi( field.filament( i ), p[ 0 ], p[ 1 ],
			                            field.mu0() );

		BOOST_TEST( field.psi( p[ 0 ], p[ 1 ] ) == byHand,
		            boost::test_tools::tolerance( 0.0 ) );

		double gr = 0.0;
		double gz = 0.0;
		field.gradPsi( p[ 0 ], p[ 1 ], gr, gz );

		double handR = 0.0;
		double handZ = 0.0;
		for ( std::size_t i = 0; i < field.filamentCount(); ++i )
		{
			double dr = 0.0;
			double dz = 0.0;
			meq::filamentGradPsi( field.filament( i ), p[ 0 ], p[ 1 ], dr, dz,
			                      field.mu0() );
			handR += dr;
			handZ += dz;
		}
		BOOST_TEST( gr == handR, boost::test_tools::tolerance( 0.0 ) );
		BOOST_TEST( gz == handZ, boost::test_tools::tolerance( 0.0 ) );
	}
}

// THE AXIS IS EXACTLY ZERO, AND "EXACTLY" IS THE ASSERTION.
// k^2 is an exact factor of the kernel, so every filament contributes a bit
// exact 0.0 and so does their sum. This is the one place in this file where a
// tolerance would destroy the point of the test.
BOOST_AUTO_TEST_CASE( psi_c_vanishes_on_the_axis_bit_for_bit )
{
	meq::ConductorField const field = threeFilaments();

	double const heights[] = { -2.0, -0.5, 0.0, +0.35, +1.10, +2.0 };
	for ( double z : heights )
	{
		BOOST_TEST( field.psi( 0.0, z ) == 0.0,
		            boost::test_tools::tolerance( 0.0 ) );
	}

	// And an EMPTY set is zero everywhere, which is what makes the split inert
	// for a machine with no filaments -- the property CS-2 will lean on to keep
	// a fixed-boundary case unchanged.
	meq::ConductorField const none;
	BOOST_TEST( none.empty() );
	BOOST_TEST( none.psi( 1.7, 0.3 ) == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );
	BOOST_TEST( none.totalCurrent() == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );
}

// THE FLUX IS THE GRADIENT OF THE SUM DIVIDED ONCE, NOT A SUM OF QUOTIENTS.
// Both are the same number in exact arithmetic; only one of them is a single
// division, and only one of them is NaN once rather than per filament on the
// axis. The test pins the identity q = ( 1/r ) grad_bar psi against gradPsi()
// so that a future rewrite cannot quietly go back to summing filamentFlux().
BOOST_AUTO_TEST_CASE( the_flux_is_one_over_r_times_the_gradient_of_the_sum )
{
	meq::ConductorField const field = threeFilaments();

	double const r = 1.55;
	double const z = -0.25;

	double gr = 0.0;
	double gz = 0.0;
	field.gradPsi( r, z, gr, gz );

	double qr = 0.0;
	double qz = 0.0;
	field.flux( r, z, qr, qz );

	BOOST_TEST( qr == gr/r, boost::test_tools::tolerance( 0.0 ) );
	BOOST_TEST( qz == gz/r, boost::test_tools::tolerance( 0.0 ) );
}

// THE TOTAL CURRENT IS THE NUMBER A BOUNDARY INTEGRAL IS CHECKED AGAINST.
BOOST_AUTO_TEST_CASE( the_total_current_is_the_signed_sum )
{
	meq::ConductorField const field = threeFilaments();
	BOOST_TEST( field.totalCurrent() == 3.0e5 - 1.7e5 + 9.0e4,
	            boost::test_tools::tolerance( 1.0e-15 ) );
}

// THE COINCIDENCE TEST CATCHES A NODE ON A FILAMENT AND NOTHING ELSE.
//
// This is CS-1's policy, and it is a COINCIDENCE test rather than a CLEARANCE
// test on purpose: Coils.hpp measures filamentPsi() still correct at eps = 1e-13
// from the conductor, where the textbook form has been NaN since 1e-09, so a
// point NEAR a filament is not a problem and only a point ON one is. The two
// halves below are what say the tolerance is in the right place -- an exact hit
// and a round-tripped hit are caught, and a point a millimetre away is not.
BOOST_AUTO_TEST_CASE( a_node_on_a_filament_is_caught_and_a_node_near_one_is_not )
{
	meq::ConductorField const field = threeFilaments();

	// Exactly on each of them, in insertion order.
	BOOST_TEST( field.indexAt( 1.20, -0.80 ) == 0 );
	BOOST_TEST( field.indexAt( 2.05, +0.35 ) == 1 );
	BOOST_TEST( field.indexAt( 0.75, +1.10 ) == 2 );
	BOOST_TEST( field.coincides( 1.20, -0.80 ) );

	// A COORDINATE THAT HAS BEEN THROUGH A TEXT ROUND TRIP, which is what a
	// .msh actually delivers: the same point, re-read from 17 significant
	// digits. This is the case the tolerance exists for, and an exact-equality
	// test would miss it.
	double const roundTripped = std::stod( "1.2000000000000002" );
	BOOST_TEST( roundTripped != 1.20 );
	BOOST_TEST( field.coincides( roundTripped, -0.80 ) );

	// AND A POINT THAT IS MERELY CLOSE IS FINE, at a separation far below any
	// mesh size and far above the tolerance. psi is large there and that is
	// correct rather than wrong -- nothing approximates psi_c by a polynomial.
	BOOST_TEST( !field.coincides( 1.20 + 1.0e-6, -0.80 ) );
	BOOST_TEST( field.indexAt( 1.20 + 1.0e-6, -0.80 ) == -1 );
	BOOST_TEST( std::isfinite( field.psi( 1.20 + 1.0e-9, -0.80 ) ) );

	// A point nowhere near anything.
	BOOST_TEST( field.indexAt( 1.70, 0.00 ) == -1 );
}

// AND THE REFUSAL BEHIND IT IS STILL THERE, which is what makes the
// coincidence test an EARLY warning rather than the only one. A caller who
// skips indexAt() does not get a wrong answer; it gets filamentPsi()'s throw.
BOOST_AUTO_TEST_CASE( evaluating_on_a_filament_still_throws )
{
	meq::ConductorField const field = threeFilaments();
	BOOST_CHECK_THROW( field.psi( 2.05, +0.35 ), std::invalid_argument );

	double a = 0.0;
	double b = 0.0;
	BOOST_CHECK_THROW( field.gradPsi( 2.05, +0.35, a, b ),
	                   std::invalid_argument );
}

// CS-1b: RECTANGLES SUM IN BESIDE THE FILAMENTS, AND AGAINST meq::CoilSet
// TERM FOR TERM. Superposition again, and again as an identity rather than a
// tolerance: the two arms are the same operations in the same order, so
// anything but bit equality is a real difference.
BOOST_AUTO_TEST_CASE( rectanglesSumInBesideTheFilaments )
{
	meq::ConductorField field = threeFilaments();
	field.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );
	field.add( meq::Coil( 0.90, +0.70, 0.05, 0.05, +2.1e5 ) );

	BOOST_TEST( field.filamentCount() == 3u );
	BOOST_TEST( field.coilCount() == 2u );
	BOOST_TEST( field.size() == 5u );
	BOOST_TEST( !field.empty() );

	// The rectangles' own set answers exactly what the pair added.
	meq::CoilSet bare;
	bare.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );
	bare.add( meq::Coil( 0.90, +0.70, 0.05, 0.05, +2.1e5 ) );

	double const points[][ 2 ] = { { 1.70, 0.00 }, { 0.40, -1.90 },
	                               { 2.60, +0.90 } };
	for ( auto const &p : points )
	{
		double filaments = 0.0;
		for ( std::size_t i = 0; i < field.filamentCount(); ++i )
			filaments += meq::filamentPsi( field.filament( i ), p[ 0 ], p[ 1 ],
			                               field.mu0() );

		BOOST_TEST( field.psi( p[ 0 ], p[ 1 ] )
		            == filaments + bare.psi( p[ 0 ], p[ 1 ] ),
		            boost::test_tools::tolerance( 0.0 ) );
	}

	// And the axis is still exactly zero with rectangles present, which is the
	// free-boundary boundary condition and must survive every conductor kind.
	BOOST_TEST( field.psi( 0.0, 0.25 ) == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );

	BOOST_TEST( field.totalCurrent()
	            == 3.0e5 - 1.7e5 + 9.0e4 - 4.4e5 + 2.1e5,
	            boost::test_tools::tolerance( 1.0e-15 ) );
}

// A MESH POINT INSIDE A RECTANGLE IS AN ORDINARY POINT, AND THAT IS A CONTRACT.
//
// coincides() is about filaments only. A coil's field is a quadrature of the
// filament kernel over its cross-section and meq::coilPsi() is documented
// "Valid EVERYWHERE, including inside the coil" -- so there is no line
// singularity for a node to land on, and refusing one would reject exactly the
// configuration this plan exists to make cheap.
BOOST_AUTO_TEST_CASE( aPointInsideARectangleIsNotACoincidence )
{
	meq::ConductorField field;
	field.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );

	// Dead centre of the coil, where a filament would be infinite.
	BOOST_TEST( !field.coincides( 1.60, -0.30 ) );
	BOOST_TEST( field.indexAt( 1.60, -0.30 ) == -1 );
	BOOST_TEST( std::isfinite( field.psi( 1.60, -0.30 ) ) );

	double gr = 0.0;
	double gz = 0.0;
	BOOST_CHECK_NO_THROW( field.gradPsi( 1.60, -0.30, gr, gz ) );
	BOOST_TEST( std::isfinite( gr ) );
	BOOST_TEST( std::isfinite( gz ) );

	// A corner, and a point on the edge: still ordinary.
	BOOST_TEST( !field.coincides( 1.60 + 0.08, -0.30 + 0.12 ) );
	BOOST_TEST( std::isfinite( field.psi( 1.60 + 0.08, -0.30 ) ) );
}

// THE QUADRATURE ORDER FORWARDS, which is what COIL-SUBTRACTION-PLAN.md §7.2's
// replacement for CS-5 needs: a reference field built far beyond what a solve
// would use, on the same conductors and the same code.
BOOST_AUTO_TEST_CASE( theQuadratureOrderForwardsAndChangesOnlyTheRectangles )
{
	meq::ConductorField field;
	field.add( meq::CurrentFilament( 1.20, -0.80, +3.0e5 ) );

	int const shipped = field.quadratureOrder();
	double const filamentOnly = field.psi( 1.70, 0.0 );

	field.setQuadratureOrder( shipped + 8 );
	BOOST_TEST( field.quadratureOrder() == shipped + 8 );

	// A filament has no quadrature, so raising the order must not move it by a
	// single bit -- which is what says the order reaches the rectangles alone.
	BOOST_TEST( field.psi( 1.70, 0.0 ) == filamentOnly,
	            boost::test_tools::tolerance( 0.0 ) );

	// And on a rectangle it DOES move, or the forwarding is inert and the test
	// above proves nothing.
	meq::ConductorField coarse;
	coarse.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );
	meq::ConductorField fine;
	fine.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );
	fine.setQuadratureOrder( 4 );

	BOOST_TEST( coarse.quadratureOrder() != fine.quadratureOrder() );
	BOOST_TEST( coarse.psi( 1.63, -0.28 ) != fine.psi( 1.63, -0.28 ),
	            boost::test_tools::tolerance( 0.0 ) );
}

// THE REFUSALS ARE THE CONTRACT.
BOOST_AUTO_TEST_CASE( the_refusals_are_the_contract )
{
	BOOST_CHECK_THROW( meq::ConductorField( 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::ConductorField( -1.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::ConductorField(
		std::numeric_limits<double>::quiet_NaN() ), std::invalid_argument );

	meq::ConductorField field = threeFilaments();
	BOOST_CHECK_THROW( field.filament( 3 ), std::out_of_range );
	BOOST_CHECK_THROW( field.setCoincidenceTolerance( -1.0e-12 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( field.indexAt(
		std::numeric_limits<double>::infinity(), 0.0 ),
		std::invalid_argument );

	// Zero is ALLOWED and means exact equality -- a defensible choice for a
	// caller whose coordinates never round-trip, and the setting under which
	// the round-tripped point above is correctly NOT a coincidence.
	BOOST_CHECK_NO_THROW( field.setCoincidenceTolerance( 0.0 ) );
	BOOST_TEST( field.coincides( 1.20, -0.80 ) );
	BOOST_TEST( !field.coincides( std::stod( "1.2000000000000002" ), -0.80 ) );
}

// THE AXIS IS THE ONE PLACE q HAS A LIMIT RATHER THAN A VALUE, AND
// poloidalField() IS THE ONLY ENTRY POINT THAT TAKES IT.
//
// flux() is 0/0 at r = 0 and reports NaN deliberately -- a caller who has not
// thought about the axis learns something from a NaN and nothing from a
// plausible number. But an output grid on a half-disc machine has its WHOLE
// FIRST COLUMN on r = 0, so somebody has to take the limit, and a limit taken
// wrongly is invisible: it would be one column of a 129 x 129 file, finite,
// smooth against its neighbours, and wrong.
//
// THREE THINGS ARE ASSERTED AND THEY FAIL IN DIFFERENT DIRECTIONS.
BOOST_AUTO_TEST_CASE( the_axis_limit_is_the_limit_and_not_a_substitute )
{
	double const a = 1.5;
	double const h = 0.3;
	double const current = 1.0e6;

	meq::ConductorField field;
	field.add( meq::CurrentFilament( a, h, current ) );

	// (1) THE CLOSED FORM. On the axis every point of the ring is the same
	// distance away, which is the whole reason this case needs no elliptic
	// integral: q_r -> mu0 I a^2 / ( 2 d^3 ).
	double const z = -0.4;
	double const d = std::hypot( a, z - h );
	double const expected =
		0.5*meq::vacuumPermeability*current*a*a/( d*d*d );

	double bR = 0.0;
	double bZ = 0.0;
	field.poloidalField( 0.0, z, bR, bZ );
	BOOST_TEST( bZ == expected, boost::test_tools::tolerance( 1.0e-15 ) );

	// (2) B_R IS EXACTLY ZERO AND NOT MERELY SMALL. psi ~ c( z ) r^2 near the
	// axis, so d_z psi ~ c'( z ) r^2 and q_z ~ c'( z ) r -- the limit is zero
	// identically, and a tolerance here would accept a limit taken by
	// evaluating at some small r instead of by algebra.
	BOOST_TEST( bR == 0.0, boost::test_tools::tolerance( 0.0 ) );

	// (3) AND IT IS THE LIMIT OF THE THING IT REPLACES, which is what says the
	// two kernels describe one field rather than two. flux() off the axis
	// approaches it from the side, at a rate that halves the error as r halves
	// -- so the check is that the approach happens at all and lands where the
	// closed form says, not that any one r is close.
	double previous = 0.0;
	for ( double r : { 1.0e-2, 1.0e-3, 1.0e-4 } )
	{
		double qR = 0.0;
		double qZ = 0.0;
		field.flux( r, z, qR, qZ );
		double const error = std::abs( qR - expected );
		if ( previous > 0.0 )
			BOOST_TEST( error < previous );
		previous = error;
	}
	BOOST_TEST( previous < 1.0e-6*std::abs( expected ) );

	// AND flux() ITSELF STILL REFUSES TO GUESS, which is the contract the
	// header states and the reason poloidalField() had to be written at all.
	double qR = 0.0;
	double qZ = 0.0;
	field.flux( 0.0, z, qR, qZ );
	BOOST_TEST( std::isnan( qR ) );
	BOOST_TEST( std::isnan( qZ ) );
}

// A RECTANGLE'S AXIS LIMIT IS THE FILAMENT'S IN THE THIN LIMIT, AND THE
// QUADRATURE IS WHAT IS BEING CHECKED.
//
// meq::coilAxisFlux() integrates a^2/( 2 d^3 ) over the cross-section on a
// PLAIN tensor Gauss rule -- no panelling and no grading, because the axis is
// outside every meq::Coil by that class's own refusal and the integrand is
// analytic there. If that reasoning were wrong the rule would be integrating a
// near-singular function with no grading, and the symptom would be a value that
// is merely a few per cent off rather than a failure.
BOOST_AUTO_TEST_CASE( a_thin_rectangle_reaches_its_filament_on_the_axis )
{
	double const a = 1.5;
	double const h = 0.3;
	double const current = -7.5e5;
	double const z = 0.9;

	meq::ConductorField filament;
	filament.add( meq::CurrentFilament( a, h, current ) );
	double bR = 0.0;
	double reference = 0.0;
	filament.poloidalField( 0.0, z, bR, reference );

	double previous = 0.0;
	for ( double halfWidth : { 1.0e-1, 1.0e-2, 1.0e-3 } )
	{
		meq::ConductorField rectangle;
		rectangle.add( meq::Coil( a, h, halfWidth, halfWidth, current ) );

		double thisBR = 0.0;
		double thisBZ = 0.0;
		rectangle.poloidalField( 0.0, z, thisBR, thisBZ );

		// The rectangle's own B_R on the axis is zero for the same reason the
		// filament's is, and it is the sum of per-conductor zeros rather than a
		// cancellation.
		BOOST_TEST( thisBR == 0.0, boost::test_tools::tolerance( 0.0 ) );

		double const error = std::abs( thisBZ - reference );
		if ( previous > 0.0 )
			BOOST_TEST( error < previous );
		previous = error;
	}
	BOOST_TEST( previous < 1.0e-5*std::abs( reference ) );

	// AND IT SUPERPOSES, which is the property the whole split rests on and is
	// as easy to break on this path as on psi's: a factor applied per conductor
	// rather than once gives the right answer for one and the wrong one for two.
	meq::ConductorField pair;
	pair.add( meq::CurrentFilament( a, h, current ) );
	pair.add( meq::CurrentFilament( 0.8, -0.2, 3.0e5 ) );

	meq::ConductorField second;
	second.add( meq::CurrentFilament( 0.8, -0.2, 3.0e5 ) );

	double pairR = 0.0;
	double pairZ = 0.0;
	double secondR = 0.0;
	double secondZ = 0.0;
	pair.poloidalField( 0.0, z, pairR, pairZ );
	second.poloidalField( 0.0, z, secondR, secondZ );

	BOOST_TEST( pairZ == reference + secondZ,
	            boost::test_tools::tolerance( 1.0e-15 ) );
}

// PSI_C IS EVEN IN r, AND THAT IS WHAT LETS AN EXTRAPOLATING EVALUATION OFF
// THE HALF-PLANE BE ANSWERED INSTEAD OF KILLING A SOLVE.
//
// Two of MEQ's evaluations of psi_c reach points with r < 0 by design and
// neither is a geometry error. The exterior datum is assembled at transfer
// path TARGETS on Gamma, and on a half-disc machine Gamma is a semicircle
// whose two endpoints lie exactly ON the axis -- so a target near an endpoint
// lands either side of r = 0 by an amount that is a property of the path map
// rather than of the mesh's validity. Measured on the DIII-D filament case at
// k = 3, that point is ( r, z ) = ( -1.0984e-03, 3.4 ), which is Gamma's own
// upper endpoint and about one per cent of h. The other is
// meq::CriticalPointFinder, whose element-local Newton is allowed to leave its
// element by up to 2 in reference coordinates on purpose.
//
// psi = r A_phi, and under r -> -r at fixed z the point is the same physical
// point rotated by pi in phi: phi-hat reverses, A_phi changes sign, and the
// product does not. So the reflection is the ANALYTIC CONTINUATION and the
// equality below is exact rather than a tolerance -- which is the assertion
// that separates it from a clamp to r = 0, a substitution that would also
// "work" and would be wrong by O( r^2 ).
//
// AND THE VECTOR ENTRY POINTS MUST STILL REFUSE, which is the other half and
// the one a careless widening of this fix would break: d_r psi is ODD where
// psi is even, so the three of them continue with a sign that differs BETWEEN
// THE TWO ENTRIES of one vector, and there is no single rule a caller holding
// one can apply. meq::CriticalPointFinder::totalFlux is the seam that meets
// this for q and it abandons the evaluation rather than continuing.
BOOST_AUTO_TEST_CASE( psiIsEvenInRAndTheVectorEntryPointsStillRefuse )
{
	meq::ConductorField field = threeFilaments();
	field.add( meq::Coil( 1.60, -0.30, 0.08, 0.12, -4.4e5 ) );
	field.add( meq::Coil( 0.90, +0.70, 0.05, 0.05, +2.1e5 ) );

	// The last of these is the measured point from the DIII-D case, to the
	// digits the throw reported it at.
	double const points[][ 2 ] = { { 1.70, 0.00 }, { 0.40, -1.90 },
	                               { 2.60, +0.90 }, { 1.0984e-03, 3.40 } };
	for ( auto const &p : points )
	{
		BOOST_TEST( field.psi( -p[ 0 ], p[ 1 ] ) == field.psi( p[ 0 ], p[ 1 ] ),
		            boost::test_tools::tolerance( 0.0 ) );
	}

	// AND IT IS NOT A CLAMP. psi_c ~ c( z ) r^2 near the axis, so reflecting
	// and clamping differ at second order -- small, and not zero. A point at
	// r = 1.0984e-03 must give the SAME number as its reflection and a
	// DIFFERENT one from the axis, or the continuation has been replaced by a
	// substitution that this test would otherwise pass.
	double const near = field.psi( -1.0984e-03, 3.40 );
	BOOST_TEST( field.psi( 0.0, 3.40 ) == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );
	BOOST_TEST( near != 0.0 );

	// -0.0 reflects to +0.0 and the axis value is still exactly zero, so the
	// bit-for-bit boundary condition psi_c_vanishes_on_the_axis_bit_for_bit
	// asserts survives the reflection.
	BOOST_TEST( field.psi( -0.0, 3.40 ) == 0.0,
	            boost::test_tools::tolerance( 0.0 ) );

	double a = 0.0;
	double b = 0.0;
	BOOST_CHECK_THROW( field.gradPsi( -1.0984e-03, 3.40, a, b ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( field.flux( -1.0984e-03, 3.40, a, b ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( field.poloidalField( -1.0984e-03, 3.40, a, b ),
	                   std::invalid_argument );
}
