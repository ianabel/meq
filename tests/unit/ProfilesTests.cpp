// Unit tests for meq::Profile and the profile implementations in src/meq/Profiles.
//
// A profile is a user-supplied function of the normalised flux, and it is half of
// the Grad-Shafranov right hand side; the other half is geometry. Two properties
// of this file matter more than the rest:
//
//   * A Hermite cubic interpolates a value *and* a slope at each knot, so it
//     reproduces any cubic exactly. That is a sharp, exact test -- not a
//     convergence rate -- and it is the first thing here.
//   * prime() is the analytic derivative of operator(), and meq's Newton solve
//     depends on the two agreeing. The finite-difference tests below check that
//     for the profile; SourceTests.cpp checks it again for the full source term.

#define BOOST_TEST_MODULE ProfilesTests
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#include <boost/test/unit_test.hpp>

#include "meq/Profiles.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace
{

	// f( x ) = 2x^3 - x^2 + 3x - 1, with its exact derivative. A Hermite cubic
	// must reproduce this to round-off; note that it has a root near x = 0.34, so
	// comparisons use an absolute floor rather than a pure relative tolerance.
	double cubic( double x )
	{
		return 2.0*x*x*x - x*x + 3.0*x - 1.0;
	}

	double cubicPrime( double x )
	{
		return 6.0*x*x - 2.0*x + 3.0;
	}

	// Something a cubic cannot be: smooth, non-polynomial, and with a derivative
	// that changes sign several times over [ 0, 1 ].
	double wiggle( double x )
	{
		return std::exp( std::sin( 5.0*x ) ) + 0.3*x*x;
	}

	double wigglePrime( double x )
	{
		return 5.0*std::cos( 5.0*x )*std::exp( std::sin( 5.0*x ) ) + 0.6*x;
	}

	// Absolute comparison with a floor of one, so that a value passing through
	// zero does not make a relative tolerance meaningless.
	void checkClose( double actual, double expected, double tolerance, char const * what, double where )
	{
		double const scale = std::max( 1.0, std::fabs( expected ) );
		BOOST_CHECK_MESSAGE( std::fabs( actual - expected ) <= tolerance*scale,
			what << " at " << where << ": got " << actual << ", expected " << expected
			<< " (difference " << actual - expected << ", allowed " << tolerance*scale << ")" );
	}

	template<typename Function>
	double centralDifference( Function const & f, double x, double h )
	{
		return ( f( x + h ) - f( x - h ) )/( 2.0*h );
	}

	meq::SplineProfile makeCubicSpline( unsigned int intervals )
	{
		return meq::SplineProfile( cubic, cubicPrime, intervals );
	}

}

BOOST_AUTO_TEST_SUITE( constant_profile_tests )

BOOST_AUTO_TEST_CASE( returns_its_value_everywhere_and_a_zero_derivative )
{
	meq::ConstantProfile const profile( 2.5 );

	BOOST_CHECK_EQUAL( profile.value(), 2.5 );

	for ( double psi : { -3.0, -1e-12, 0.0, 0.25, 1.0, 1.0 + 1e-12, 7.5 } )
	{
		BOOST_CHECK_EQUAL( profile( psi ), 2.5 );
		// Exactly zero, not merely small: a Solov'ev-like profile must contribute
		// nothing at all to the Newton Jacobian.
		BOOST_CHECK_EQUAL( profile.prime( psi ), 0.0 );
	}
}

BOOST_AUTO_TEST_CASE( is_usable_through_the_base_class )
{
	std::shared_ptr<meq::Profile> const profile = std::make_shared<meq::ConstantProfile>( -1.25 );
	meq::Profile const & base = *profile;

	BOOST_CHECK_EQUAL( base( 0.4 ), -1.25 );
	BOOST_CHECK_EQUAL( base.prime( 0.4 ), 0.0 );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( hermite_cubic_tests )

BOOST_AUTO_TEST_CASE( reproduces_a_cubic_exactly_on_one_interval )
{
	double const lower = 0.25, upper = 0.75;
	meq::HermiteCubicSpline const spline( lower, upper, cubic( lower ), cubic( upper ), cubicPrime( lower ), cubicPrime( upper ) );

	for ( int i = 0; i <= 100; ++i )
	{
		double const x = lower + ( upper - lower )*i/100.0;
		checkClose( spline( x ), cubic( x ), 1e-14, "value", x );
		checkClose( spline.prime( x ), cubicPrime( x ), 1e-13, "derivative", x );
	}
}

BOOST_AUTO_TEST_CASE( matches_the_data_at_both_endpoints )
{
	meq::HermiteCubicSpline const spline( -1.0, 2.0, 3.0, -4.0, 0.5, 7.0 );

	checkClose( spline( -1.0 ), 3.0, 1e-15, "value", -1.0 );
	checkClose( spline( 2.0 ), -4.0, 1e-15, "value", 2.0 );
	checkClose( spline.prime( -1.0 ), 0.5, 1e-13, "derivative", -1.0 );
	checkClose( spline.prime( 2.0 ), 7.0, 1e-13, "derivative", 2.0 );

	BOOST_CHECK_EQUAL( spline.interval().first, -1.0 );
	BOOST_CHECK_EQUAL( spline.interval().second, 2.0 );
	BOOST_CHECK_EQUAL( spline.values().first, 3.0 );
	BOOST_CHECK_EQUAL( spline.values().second, -4.0 );
	BOOST_CHECK_EQUAL( spline.derivatives().first, 0.5 );
	BOOST_CHECK_EQUAL( spline.derivatives().second, 7.0 );
}

BOOST_AUTO_TEST_CASE( clamps_to_a_constant_outside_its_interval )
{
	// The documented policy: outside the interval the interpolant is extended by
	// a constant, so the value is the nearer endpoint's and the derivative is
	// zero. The class this replaces returned f( x_l ) -- a *value* -- from
	// prime(), which is both the wrong number and the wrong units.
	meq::HermiteCubicSpline const spline( 0.0, 1.0, 3.0, -4.0, 0.5, 7.0 );

	BOOST_CHECK_EQUAL( spline( -1.0 ), 3.0 );
	BOOST_CHECK_EQUAL( spline( 2.0 ), -4.0 );
	BOOST_CHECK_EQUAL( spline.prime( -1.0 ), 0.0 );
	BOOST_CHECK_EQUAL( spline.prime( 2.0 ), 0.0 );
}

BOOST_AUTO_TEST_CASE( rejects_a_degenerate_interval )
{
	BOOST_CHECK_THROW( meq::HermiteCubicSpline( 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::HermiteCubicSpline( 1.0, 0.5, 0.0, 0.0, 0.0, 0.0 ), std::invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( spline_profile_tests )

BOOST_AUTO_TEST_CASE( reproduces_a_cubic_exactly )
{
	// Seven intervals, so most sample points are interior to an interval and some
	// land on knots. A piecewise Hermite cubic through exact cubic data is that
	// cubic, everywhere, to round-off.
	meq::SplineProfile const profile = makeCubicSpline( 7 );

	for ( int i = 0; i <= 500; ++i )
	{
		double const psi = i/500.0;
		checkClose( profile( psi ), cubic( psi ), 1e-14, "value", psi );
		checkClose( profile.prime( psi ), cubicPrime( psi ), 1e-12, "derivative", psi );
	}
}

BOOST_AUTO_TEST_CASE( interpolates_the_data_at_every_knot )
{
	std::vector<meq::Knot> data;
	for ( int i = 0; i <= 11; ++i )
	{
		double const psi = i/11.0;
		data.push_back( meq::Knot{ psi, wiggle( psi ), wigglePrime( psi ) } );
	}

	meq::SplineProfile const profile( data );

	BOOST_REQUIRE_EQUAL( profile.numIntervals(), data.size() - 1 );
	BOOST_CHECK_EQUAL( profile.domain().first, data.front().psi );
	BOOST_CHECK_EQUAL( profile.domain().second, data.back().psi );

	for ( auto const & knot : data )
	{
		// Both endpoints of every interior knot's two intervals; the interpolant
		// has to agree with the data, not merely be close to it.
		checkClose( profile( knot.psi ), knot.value, 1e-14, "value at knot", knot.psi );
		checkClose( profile.prime( knot.psi ), knot.derivative, 1e-12, "derivative at knot", knot.psi );
	}
}

BOOST_AUTO_TEST_CASE( is_continuous_across_interior_knots )
{
	meq::SplineProfile const profile( wiggle, wigglePrime, 11 );
	double const eps = 1e-9;

	for ( std::size_t i = 1; i < profile.knots().size() - 1; ++i )
	{
		double const psi = profile.knots()[ i ].psi;

		// C^0 across the knot...
		checkClose( profile( psi - eps ), profile( psi + eps ), 1e-8, "one-sided values", psi );
		// ...and C^1, which is the point of tabulating the slope as well as the
		// value.
		checkClose( profile.prime( psi - eps ), profile.prime( psi + eps ), 1e-7, "one-sided derivatives", psi );
	}
}

BOOST_AUTO_TEST_CASE( evaluates_on_knots_and_at_both_endpoints )
{
	// The interval lookup this replaces was a hand-rolled bisection that could
	// only be shown correct by an invariant argument about its up-front test of
	// interval zero. These are the cases that argument turned on: psi exactly on
	// an interior knot, and psi exactly at each end of the table.
	meq::SplineProfile const profile = makeCubicSpline( 8 );

	for ( auto const & knot : profile.knots() )
	{
		BOOST_CHECK_NO_THROW( profile( knot.psi ) );
		checkClose( profile( knot.psi ), cubic( knot.psi ), 1e-14, "value on knot", knot.psi );
		checkClose( profile.prime( knot.psi ), cubicPrime( knot.psi ), 1e-12, "derivative on knot", knot.psi );
	}

	checkClose( profile( 0.0 ), cubic( 0.0 ), 1e-14, "value at the lower end", 0.0 );
	checkClose( profile( 1.0 ), cubic( 1.0 ), 1e-14, "value at the upper end", 1.0 );
	checkClose( profile.prime( 0.0 ), cubicPrime( 0.0 ), 1e-12, "derivative at the lower end", 0.0 );
	checkClose( profile.prime( 1.0 ), cubicPrime( 1.0 ), 1e-12, "derivative at the upper end", 1.0 );

	// And a psi just inside each end, which the lookup must not push out of range.
	checkClose( profile( 1e-15 ), cubic( 0.0 ), 1e-13, "value just above the lower end", 1e-15 );
	checkClose( profile( 1.0 - 1e-15 ), cubic( 1.0 ), 1e-13, "value just below the upper end", 1.0 );
}

BOOST_AUTO_TEST_CASE( prime_agrees_with_a_central_difference )
{
	// prime() is the derivative of the *interpolant*, so the comparison is against
	// a difference of the interpolant, not of wiggle. Sample away from the knots:
	// the interpolant is only C^1 there, so a stencil straddling a knot carries an
	// O( h ) term from the jump in the second derivative.
	meq::SplineProfile const profile( wiggle, wigglePrime, 40 );
	double const h = 1e-6;

	for ( int i = 0; i < 40; ++i )
	{
		double const psi = ( i + 0.5 )/40.0;
		double const difference = centralDifference( profile, psi, h );
		checkClose( profile.prime( psi ), difference, 1e-7, "prime against a central difference", psi );
	}

	// Independently of the derivative check: the interpolant should also be a good
	// approximation to the function it was sampled from. Hermite cubic
	// interpolation is O( h^4 ) in the value and O( h^3 ) in the derivative, so
	// halving h twice from the 40 intervals above buys a factor of 256 and 64
	// respectively -- these tolerances are set from that, not from watching the
	// test pass.
	meq::SplineProfile const fine( wiggle, wigglePrime, 160 );

	for ( int i = 0; i < 160; ++i )
	{
		double const psi = ( i + 0.5 )/160.0;
		checkClose( fine( psi ), wiggle( psi ), 1e-7, "spline against the sampled function", psi );
		checkClose( fine.prime( psi ), wigglePrime( psi ), 1e-5, "spline derivative against the exact one", psi );
	}
}

BOOST_AUTO_TEST_CASE( clamps_outside_the_table_rather_than_throwing )
{
	// The documented policy, and it matters: a Newton iterate overshoots [ 0, 1 ]
	// routinely, and an exception out of a quadrature loop would abandon a solve
	// that was going to converge.
	meq::SplineProfile const profile( wiggle, wigglePrime, 10 );

	for ( double psi : { -5.0, -1.0, -1e-12 } )
	{
		BOOST_CHECK_EQUAL( profile( psi ), profile.knots().front().value );
		BOOST_CHECK_EQUAL( profile.prime( psi ), 0.0 );
	}

	for ( double psi : { 1.0 + 1e-12, 2.0, 17.0 } )
	{
		BOOST_CHECK_EQUAL( profile( psi ), profile.knots().back().value );
		BOOST_CHECK_EQUAL( profile.prime( psi ), 0.0 );
	}

	// At the end knots themselves the tabulated derivative is still returned.
	checkClose( profile.prime( 0.0 ), wigglePrime( 0.0 ), 1e-13, "derivative at psi = 0", 0.0 );
	checkClose( profile.prime( 1.0 ), wigglePrime( 1.0 ), 1e-13, "derivative at psi = 1", 1.0 );
}

BOOST_AUTO_TEST_CASE( write_then_read_round_trips )
{
	meq::SplineProfile const original( wiggle, wigglePrime, 9 );

	std::ostringstream out;
	original.write( out );

	std::istringstream in( out.str() );
	meq::SplineProfile const copy = meq::SplineProfile::fromStream( in );

	BOOST_REQUIRE_EQUAL( copy.knots().size(), original.knots().size() );

	for ( std::size_t i = 0; i < original.knots().size(); ++i )
	{
		// write() emits max_digits10 significant figures, so this is exact.
		BOOST_CHECK_EQUAL( copy.knots()[ i ].psi, original.knots()[ i ].psi );
		BOOST_CHECK_EQUAL( copy.knots()[ i ].value, original.knots()[ i ].value );
		BOOST_CHECK_EQUAL( copy.knots()[ i ].derivative, original.knots()[ i ].derivative );
	}

	for ( int i = 0; i <= 100; ++i )
	{
		double const psi = i/100.0;
		BOOST_CHECK_EQUAL( copy( psi ), original( psi ) );
		BOOST_CHECK_EQUAL( copy.prime( psi ), original.prime( psi ) );
	}

	// operator<< is write().
	std::ostringstream streamed;
	streamed << original;
	BOOST_CHECK_EQUAL( streamed.str(), out.str() );
}

BOOST_AUTO_TEST_CASE( reads_comments_and_stops_at_a_blank_line )
{
	// The blank-line terminator is what lets several profiles share one stream.
	std::istringstream in(
		"# a comment\n"
		"0.0\t1.0\t0.0\n"
		"   # an indented comment\n"
		"0.5\t2.0\t1.0\n"
		"1.0\t0.0\t-1.0\n"
		"\n"
		"0.0 -1.0 0.0\n"
		"1.0 -2.0 0.0\n" );

	meq::SplineProfile const first = meq::SplineProfile::fromStream( in );
	BOOST_REQUIRE_EQUAL( first.knots().size(), 3u );
	BOOST_CHECK_EQUAL( first( 0.5 ), 2.0 );
	BOOST_CHECK_EQUAL( first.prime( 0.5 ), 1.0 );

	meq::SplineProfile const second = meq::SplineProfile::fromStream( in );
	BOOST_REQUIRE_EQUAL( second.knots().size(), 2u );
	BOOST_CHECK_EQUAL( second( 0.0 ), -1.0 );
	BOOST_CHECK_EQUAL( second( 1.0 ), -2.0 );
}

BOOST_AUTO_TEST_CASE( rejects_input_it_cannot_turn_into_a_spline )
{
	// Not two knots.
	{
		std::istringstream in( "0.0 1.0 0.0\n" );
		BOOST_CHECK_THROW( meq::SplineProfile::fromStream( in ), std::runtime_error );
	}

	// Nothing at all.
	{
		std::istringstream in( "" );
		BOOST_CHECK_THROW( meq::SplineProfile::fromStream( in ), std::runtime_error );
	}

	// A line that is not three numbers.
	{
		std::istringstream in( "0.0 1.0 0.0\nthis is not data\n" );
		BOOST_CHECK_THROW( meq::SplineProfile::fromStream( in ), std::runtime_error );
	}

	// Two numbers where three are needed: the derivative is data, not optional.
	{
		std::istringstream in( "0.0 1.0 0.0\n1.0 2.0\n" );
		BOOST_CHECK_THROW( meq::SplineProfile::fromStream( in ), std::runtime_error );
	}

	// Knots out of order, or repeated, would give an interval of zero or negative
	// width.
	{
		std::vector<meq::Knot> const backwards{ { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } };
		BOOST_CHECK_THROW( meq::SplineProfile{ backwards }, std::invalid_argument );

		std::vector<meq::Knot> const repeated{ { 0.0, 0.0, 0.0 }, { 0.5, 1.0, 0.0 }, { 0.5, 2.0, 0.0 } };
		BOOST_CHECK_THROW( meq::SplineProfile{ repeated }, std::invalid_argument );
	}

	// One knot, or none, is not a spline.
	{
		std::vector<meq::Knot> const single{ { 0.0, 1.0, 0.0 } };
		BOOST_CHECK_THROW( meq::SplineProfile{ single }, std::invalid_argument );
		BOOST_CHECK_THROW( meq::SplineProfile( std::vector<meq::Knot>{} ), std::invalid_argument );
	}

	// Non-finite data.
	{
		std::vector<meq::Knot> const notFinite{ { 0.0, 1.0, 0.0 }, { 1.0, std::nan( "" ), 0.0 } };
		BOOST_CHECK_THROW( meq::SplineProfile{ notFinite }, std::invalid_argument );
	}

	// And the sampling constructor needs at least one interval.
	BOOST_CHECK_THROW( meq::SplineProfile( cubic, cubicPrime, 0 ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::SplineProfile( meq::SplineProfile::RealFunction(), cubicPrime, 4 ), std::invalid_argument );

	// intervalAt is bounds checked.
	meq::SplineProfile const profile = makeCubicSpline( 3 );
	BOOST_CHECK_THROW( profile.intervalAt( profile.numIntervals() ), std::out_of_range );
	BOOST_CHECK_NO_THROW( profile.intervalAt( profile.numIntervals() - 1 ) );
}

BOOST_AUTO_TEST_CASE( from_file_reads_a_table_and_fails_loudly_when_it_cannot )
{
	// The constructor this replaces had a body of `return;`: it silently produced
	// an empty spline table, and every configured profile then indexed out of
	// bounds on its first evaluation.
	BOOST_CHECK_THROW( meq::SplineProfile::fromFile( "/no/such/directory/no-such-profile.dat" ), std::runtime_error );

	std::filesystem::path const path = std::filesystem::temp_directory_path() / "meq-profile-round-trip.dat";
	meq::SplineProfile const original( wiggle, wigglePrime, 6 );

	{
		std::ofstream file( path );
		BOOST_REQUIRE( file.good() );
		original.write( file );
	}

	meq::SplineProfile const fromFile = meq::SplineProfile::fromFile( path.string() );
	BOOST_CHECK_EQUAL( fromFile.knots().size(), original.knots().size() );
	for ( int i = 0; i <= 50; ++i )
	{
		double const psi = i/50.0;
		BOOST_CHECK_EQUAL( fromFile( psi ), original( psi ) );
		BOOST_CHECK_EQUAL( fromFile.prime( psi ), original.prime( psi ) );
	}

	// A file that exists but holds nothing usable must throw, not produce an
	// empty profile.
	{
		std::ofstream file( path );
		file << "# only a comment\n";
	}
	BOOST_CHECK_THROW( meq::SplineProfile::fromFile( path.string() ), std::runtime_error );

	std::filesystem::remove( path );
}

BOOST_AUTO_TEST_CASE( copies_and_moves_are_value_semantics )
{
	// The special members are compiler generated -- the hand-written "move"
	// constructor this replaces copied its argument, which a nothrow move
	// constructor cannot do for a vector member.
	static_assert( std::is_copy_constructible_v<meq::SplineProfile>, "a profile should be copyable" );
	static_assert( std::is_nothrow_move_constructible_v<meq::SplineProfile>, "the move constructor should not be a copy" );

	meq::SplineProfile original = makeCubicSpline( 5 );
	meq::SplineProfile const copy = original;
	meq::SplineProfile const moved = std::move( original );

	for ( int i = 0; i <= 50; ++i )
	{
		double const psi = i/50.0;
		checkClose( copy( psi ), cubic( psi ), 1e-14, "copied profile", psi );
		checkClose( moved( psi ), cubic( psi ), 1e-14, "moved profile", psi );
		checkClose( moved.prime( psi ), cubicPrime( psi ), 1e-12, "moved profile derivative", psi );
	}
}

BOOST_AUTO_TEST_CASE( is_usable_through_the_base_class )
{
	std::vector<std::shared_ptr<meq::Profile>> profiles;
	profiles.push_back( std::make_shared<meq::ConstantProfile>( 4.0 ) );
	profiles.push_back( std::make_shared<meq::SplineProfile>( makeCubicSpline( 4 ) ) );

	BOOST_CHECK_EQUAL( ( *profiles[ 0 ] )( 0.3 ), 4.0 );
	checkClose( ( *profiles[ 1 ] )( 0.3 ), cubic( 0.3 ), 1e-14, "spline through Profile&", 0.3 );
	checkClose( profiles[ 1 ]->prime( 0.3 ), cubicPrime( 0.3 ), 1e-12, "spline derivative through Profile&", 0.3 );
}

BOOST_AUTO_TEST_SUITE_END()
