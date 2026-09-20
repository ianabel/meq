// Unit tests for the conductor field in src/meq/ConductorField.
//
// COIL-SUBTRACTION-PLAN.md CS-1. The split is psi = psi_c + psi_p, with psi_c
// the conductors' own field evaluated exactly and psi_p what MEQ solves for, so
// that the conductors need not be in the mesh at all. This file is the unit
// under that: psi_c, and the refusal that decides which meshes may carry it.
//
// THREE THINGS HERE COULD BE SILENTLY WRONG AND NO CONVERGENCE RATE WOULD SEE
// ANY OF THEM:
//
//   * SUPERPOSITION. psi_c is a sum over filaments, and the whole plan rests on
//     Delta* being linear -- so a set of N filaments must give exactly what the
//     N single-filament fields give added up. A factor applied per filament
//     rather than once converges at full rate to the wrong equilibrium.
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
	BOOST_TEST_REQUIRE( field.size() == 3u );

	double const points[][ 2 ] = { { 1.70, 0.00 }, { 0.40, -1.90 },
	                               { 2.60, +0.90 }, { 1.05, +1.05 } };

	for ( auto const &p : points )
	{
		double byHand = 0.0;
		for ( std::size_t i = 0; i < field.size(); ++i )
			byHand += meq::filamentPsi( field.filament( i ), p[ 0 ], p[ 1 ],
			                            field.mu0() );

		BOOST_TEST( field.psi( p[ 0 ], p[ 1 ] ) == byHand,
		            boost::test_tools::tolerance( 0.0 ) );

		double gr = 0.0;
		double gz = 0.0;
		field.gradPsi( p[ 0 ], p[ 1 ], gr, gz );

		double handR = 0.0;
		double handZ = 0.0;
		for ( std::size_t i = 0; i < field.size(); ++i )
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
