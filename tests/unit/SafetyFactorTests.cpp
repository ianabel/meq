#define BOOST_TEST_MODULE SafetyFactorTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "meq/SafetyFactor.hpp"

/*
 * DRIVING AN EQUILIBRIUM BY q( psi ), the inversion half.
 *
 * MFEM-free and so gated by CI, which is the whole reason the arithmetic lives
 * apart from the extraction: `q = V' g < R^-2 >/4 pi^2` is one division per
 * surface, and everything that can go wrong in it -- the constant, the sign,
 * and above all WHICH normalised flux -- is checkable without a mesh.
 */
namespace
{
	/// A family carrying only what invertSafetyFactor() reads. The geometry is
	/// not exercised here and is deliberately absent: a fixture that traced
	/// contours would be testing the extraction as well, and the extraction has
	/// its own acceptance in FluxGridConvergence.
	meq::FluxSurfaceFamily familyOf( std::vector<double> const &normalisedFlux,
	                                 std::vector<double> const &vPrime,
	                                 std::vector<double> const &inverseRSquared )
	{
		meq::FluxSurfaceFamily family;
		for ( std::size_t i = 0; i < normalisedFlux.size(); ++i )
		{
			meq::FluxSurface surface;
			surface.normalisedFlux = normalisedFlux[ i ];
			surface.radial = std::sqrt( normalisedFlux[ i ] );
			surface.vPrime = vPrime[ i ];
			surface.inverseRSquared = inverseRSquared[ i ];
			family.surfaces.push_back( surface );
		}
		return family;
	}

	double const fourPiSquared = 4.0*M_PI*M_PI;
}

/*
 * THE INVERSION ROUND-TRIPS, WHICH PINS THE CONSTANT AND THE DIRECTION.
 *
 * Choose a g, compute the q it gives through the SAME formula
 * meq::SurfaceAverages::safetyFactor publishes, invert, and get the g back.
 * That catches a wrong 4 pi^2, a reciprocal taken the wrong way round, and
 * V' and < R^-2 > swapped -- none of which any rate table anywhere could see,
 * because all three converge beautifully to the wrong equilibrium.
 */
BOOST_AUTO_TEST_CASE( the_inversion_recovers_the_field_it_was_built_from )
{
	std::vector<double> const label = { 0.05, 0.20, 0.45, 0.70, 0.95 };
	std::vector<double> const vPrime = { 41.3, 43.0, 45.1, 47.4, 49.7 };
	std::vector<double> const inverse = { 0.404, 0.421, 0.457, 0.510, 0.581 };
	std::vector<double> const g = { 3.1, 3.05, 2.95, 2.80, 2.60 };

	std::vector<double> q( label.size() );
	for ( std::size_t i = 0; i < label.size(); ++i )
		q[ i ] = vPrime[ i ]*g[ i ]*inverse[ i ]/fourPiSquared;

	meq::FluxSurfaceFamily const family = familyOf( label, vPrime, inverse );
	meq::ToroidalField const field = meq::invertSafetyFactor(
		family, [ & ]( double psiN )
		{
			for ( std::size_t i = 0; i < label.size(); ++i )
				if ( std::abs( psiN - label[ i ] ) < 1.0e-12 )
					return q[ i ];
			throw std::logic_error( "the target was asked at a label that is "
			                        "not one of the family's own surfaces" );
		} );

	BOOST_TEST_REQUIRE( field.size() == label.size() );
	for ( std::size_t i = 0; i < label.size(); ++i )
	{
		BOOST_TEST( std::abs( field.g[ i ] - g[ i ] ) < 1.0e-13*std::abs( g[ i ] ),
			"surface " << i << ": g came back " << field.g[ i ]
			<< " against the " << g[ i ] << " the target was built from" );
		BOOST_TEST( std::abs( field.gSquared[ i ] - g[ i ]*g[ i ] )
		            < 1.0e-13*g[ i ]*g[ i ] );
	}
}

/*
 * THE TWO NORMALISED FLUXES RUN IN OPPOSITE DIRECTIONS.
 *
 * THE ONE THAT DOES NOT ANNOUNCE ITSELF. meq::FluxSurfaceFamily's Psi_N is zero
 * on the axis; meq::NormalisedSource's Psi is ONE there. So a gg' differentiated
 * against the family's label and handed to the source without the reflection
 * has the wrong sign, and the solve does not fail -- it converges, at full
 * order, to an equilibrium with its shear reversed, which is a configuration a
 * real machine can have.
 *
 * EXACTLY LINEAR DATA MAKES IT AN EQUALITY RATHER THAN A TOLERANCE. With
 * g^2 = a + b Psi_N the derivative is b everywhere, so gg' against the SOURCE's
 * Psi is -b/2 at every knot, and the whole test is one number repeated with the
 * sign the reflection puts on it.
 */
BOOST_AUTO_TEST_CASE( the_two_normalised_fluxes_run_in_opposite_directions )
{
	std::vector<double> const label = { 0.05, 0.30, 0.55, 0.80, 0.95 };
	std::vector<double> const vPrime( label.size(), 40.0 );
	std::vector<double> const inverse( label.size(), 0.5 );

	// g^2 = a + b Psi_N, RISING toward the plasma edge.
	double const a = 4.0;
	double const b = 3.0;

	meq::ToroidalField field;
	field.normalisedFlux = label;
	for ( double x : label )
	{
		field.safetyFactor.push_back( 0.0 );
		field.gSquared.push_back( a + b*x );
		field.g.push_back( std::sqrt( a + b*x ) );
	}

	std::vector<meq::Knot> const knots = meq::ggPrimeKnots( field );
	BOOST_TEST_REQUIRE( knots.size() == label.size() );

	// ASCENDING IN THE SOURCE'S Psi, which is the reflected order.
	for ( std::size_t i = 1; i < knots.size(); ++i )
		BOOST_TEST( knots[ i ].psi > knots[ i - 1 ].psi,
			"the knots are not ascending in the source's Psi, so SplineProfile "
			"would refuse them" );

	for ( std::size_t i = 0; i < knots.size(); ++i )
	{
		double const expectedPsi = 1.0 - label[ label.size() - 1 - i ];
		BOOST_TEST( std::abs( knots[ i ].psi - expectedPsi ) < 1.0e-15,
			"knot " << i << " sits at Psi = " << knots[ i ].psi
			<< " where the reflection of the family's label puts it at "
			<< expectedPsi );

		// -b/2, AND THE MINUS IS THE POINT. Against the family's own label the
		// slope is +b; the source reads the other direction.
		BOOST_TEST( std::abs( knots[ i ].value - ( -0.5*b ) ) < 1.0e-13,
			"knot " << i << " carries gg' = " << knots[ i ].value
			<< " where the reflection gives " << -0.5*b
			<< ". A value of " << +0.5*b << " is the reflection omitted, which "
			"converges to an equilibrium with its shear reversed." );

		// g^2 is linear, so gg' is constant and its own derivative is zero.
		BOOST_TEST( std::abs( knots[ i ].derivative ) < 1.0e-13 );
	}
}

/*
 * THE SLOPE RULE CANNOT OVERSHOOT, WHICH IS WHY IT IS THIS RULE.
 *
 * g^2 arrives as a table extracted from a solved field and carries that
 * extraction's noise. An unlimited cubic through noisy data overshoots, and an
 * overshoot in gg' is a source term with a sign it should not have, fed back
 * into the next solve -- which is how a fixed-point iteration that ought to
 * converge instead oscillates. The defining property is asserted directly: at a
 * local extremum of the data the slope is EXACTLY zero, and nowhere does a
 * slope disagree in sign with both of its neighbouring secants.
 */
BOOST_AUTO_TEST_CASE( the_slope_rule_is_monotone_and_exact_on_a_line )
{
	{
		// Exactly linear data: every slope is the line's, to the last bit.
		std::vector<double> x = { 0.0, 0.1, 0.35, 0.8, 1.0 };
		std::vector<double> y;
		for ( double t : x )
			y.push_back( 2.0 - 1.75*t );
		std::vector<double> const slope = meq::monotoneSlopes( x, y );
		for ( std::size_t i = 0; i < slope.size(); ++i )
			BOOST_TEST( std::abs( slope[ i ] + 1.75 ) < 1.0e-14,
				"slope " << i << " is " << slope[ i ] << " on exactly linear "
				"data, where the line's own slope is -1.75" );
	}

	{
		// A local maximum in the middle. The rule must flatten there, which is
		// what stops the interpolant rising past it.
		std::vector<double> const x = { 0.0, 0.25, 0.5, 0.75, 1.0 };
		std::vector<double> const y = { 0.0, 1.0, 2.0, 1.0, 0.0 };
		std::vector<double> const slope = meq::monotoneSlopes( x, y );
		BOOST_TEST( slope[ 2 ] == 0.0,
			"the slope at the data's own maximum is " << slope[ 2 ]
			<< " and must be exactly zero, or the interpolant overshoots it" );

		for ( std::size_t i = 0; i + 1 < x.size(); ++i )
		{
			double const secant = ( y[ i + 1 ] - y[ i ] )/( x[ i + 1 ] - x[ i ] );
			for ( std::size_t k : { i, i + 1 } )
				BOOST_TEST( slope[ k ]*secant >= -1.0e-15,
					"slope " << k << " (" << slope[ k ] << ") opposes the "
					"secant it bounds (" << secant << "), which is an "
					"overshoot" );
		}
	}
}

/*
 * WHAT IT REFUSES, AND THE CONTROL THAT KEEPS THE REFUSALS FROM BEING VACUOUS.
 *
 * A non-positive V' or < R^-2 > is a CORRUPT family rather than a hard case --
 * both are positive by construction -- and dividing by one returns a signed
 * infinity that the spline carries into the source, where the solve meets it as
 * a NaN several frames from here. Refusing at the division is the only place
 * the message can name the surface.
 */
BOOST_AUTO_TEST_CASE( a_corrupt_family_is_refused_at_the_division )
{
	auto unitTarget = []( double ) { return 1.0; };

	BOOST_CHECK_THROW( meq::invertSafetyFactor( meq::FluxSurfaceFamily{},
	                                            unitTarget ),
	                   std::invalid_argument );

	BOOST_CHECK_THROW(
		meq::invertSafetyFactor( familyOf( { 0.1, 0.5 }, { 40.0, 0.0 },
		                                   { 0.5, 0.5 } ), unitTarget ),
		std::invalid_argument );

	BOOST_CHECK_THROW(
		meq::invertSafetyFactor( familyOf( { 0.1, 0.5 }, { 40.0, 40.0 },
		                                   { 0.5, -0.5 } ), unitTarget ),
		std::invalid_argument );

	BOOST_CHECK_THROW(
		meq::invertSafetyFactor( familyOf( { 0.1, 0.5 }, { 40.0, 40.0 },
		                                   { 0.5, 0.5 } ),
		                         []( double ) { return std::nan( "" ); } ),
		std::invalid_argument );

	// THE CONTROL. The same call on a sound family must be accepted, or the
	// four refusals above are compatible with a function that refuses
	// everything.
	BOOST_CHECK_NO_THROW(
		meq::invertSafetyFactor( familyOf( { 0.1, 0.5 }, { 40.0, 40.0 },
		                                   { 0.5, 0.5 } ), unitTarget ) );

	// And ggPrimeKnots needs two surfaces to differentiate between.
	meq::ToroidalField one;
	one.normalisedFlux = { 0.5 };
	one.gSquared = { 9.0 };
	one.g = { 3.0 };
	one.safetyFactor = { 1.0 };
	BOOST_CHECK_THROW( meq::ggPrimeKnots( one ), std::invalid_argument );
}

/*
 * THE EDGE EXTENSION IS A STATEMENT ABOUT THE PHYSICS AND THE CALLER MAKES IT.
 *
 * A family stops short of the plasma edge -- Psi_N = 0.95 as shipped, i.e.
 * Psi = 0.05 -- and SplineProfile CLAMPS outside its knots. Clamping carries
 * the innermost gg' out past the edge, which keeps g^2 growing where the field
 * should be the vacuum one. VacuumOutside puts a knot at Psi = 0 carrying zero
 * instead. Neither is right for every run: with [source] ConfineToPlasma the
 * source is switched off out there and the question does not arise at all.
 */
BOOST_AUTO_TEST_CASE( the_vacuum_extension_adds_a_zero_at_the_plasma_edge )
{
	meq::ToroidalField field;
	field.normalisedFlux = { 0.05, 0.50, 0.95 };
	field.gSquared = { 9.0, 8.0, 6.0 };
	for ( double v : field.gSquared )
		field.g.push_back( std::sqrt( v ) );
	field.safetyFactor = { 1.0, 1.0, 1.0 };

	std::vector<meq::Knot> const clamped =
		meq::ggPrimeKnots( field, meq::EdgeExtension::Clamp );
	std::vector<meq::Knot> const vacuum =
		meq::ggPrimeKnots( field, meq::EdgeExtension::VacuumOutside );

	BOOST_TEST( clamped.size() == 3u );
	BOOST_TEST( vacuum.size() == 4u );
	BOOST_TEST( vacuum.front().psi == 0.0 );
	BOOST_TEST( vacuum.front().value == 0.0 );
	BOOST_TEST( vacuum.front().derivative == 0.0 );

	// AND IT CHANGES NOTHING INSIDE THE FAMILY'S OWN RANGE, which is what makes
	// it an extension rather than a different profile.
	for ( std::size_t i = 0; i < clamped.size(); ++i )
	{
		BOOST_TEST( vacuum[ i + 1 ].psi == clamped[ i ].psi );
		BOOST_TEST( vacuum[ i + 1 ].value == clamped[ i ].value );
	}
}

/*
 * THE MAP IS TOTAL, AND WHAT IT RETURNS WHERE IT CANNOT BE EVALUATED MATTERS
 * MORE THAN WHERE IT CAN.
 *
 * meq::solveForToroidalField drives this through KINSOL, which is C: an
 * exception raised inside it unwinds through its frames and denies the line
 * search the finite value it needs to reject a step with. So the contract is
 * that no reachable input throws, and the residual at an inadmissible point is
 * finite, large, and not a root.
 */
BOOST_AUTO_TEST_CASE( the_map_is_total_where_g_squared_goes_negative )
{
	std::vector<double> const label = { 0.10, 0.30, 0.50, 0.70, 0.90 };

	// A family that inverts to a CONSTANT g of 2, whatever is asked of it, so
	// the map's own arithmetic is the only thing under test.
	std::vector<double> vPrime( label.size(), 3.0 );
	std::vector<double> inverseRSquared( label.size(), 0.5 );
	meq::FluxSurfaceFamily const family =
		familyOf( label, vPrime, inverseRSquared );

	int solves = 0;
	auto solve = [ & ]( std::vector<meq::Knot> const &knots )
	{
		++solves;
		BOOST_TEST( !knots.empty() );
		return &family;
	};

	meq::ToroidalFieldMap::Options options;
	options.degree = 2;
	options.target = []( double ) { return 4.0*M_PI*M_PI*2.0/( 3.0*0.5 ); };

	meq::ToroidalFieldMap map( solve, options );

	// A g^2 that is positive everywhere on [ 0, 1 ] is admissible and solves.
	std::vector<double> const good = { 4.0, 0.0, 0.0 };
	std::vector<double> const first = map( good );
	BOOST_TEST( map.solves() == 1 );
	BOOST_TEST( map.refusals() == 0 );
	BOOST_TEST( first.size() == good.size() );

	// g^2 = 1 - 4 Psi is negative over most of the range. It must NOT reach the
	// solve, must not throw, and must come back finite.
	std::vector<double> const bad = { 1.0, -4.0, 0.0 };
	std::vector<double> hopeless;
	BOOST_CHECK_NO_THROW( hopeless = map( bad ) );
	BOOST_TEST( map.solves() == 1,
		"an inadmissible g^2 reached the solve, which is a wasted equilibrium "
		"and, with a moving support, one that may not converge at all" );
	BOOST_TEST( map.refusals() == 1 );
	for ( double value : hopeless )
		BOOST_TEST( std::isfinite( value ) );

	// AND THE REFUSAL POINTS BACK AT THE LAST ADMISSIBLE POINT, so the residual
	// there is `lastGood - c`, which is large and is a root only at a place the
	// map really was evaluated.
	BOOST_TEST( map.lastAdmissible() == good, boost::test_tools::per_element() );
	BOOST_TEST( hopeless == good, boost::test_tools::per_element() );

	std::printf( "\n  the map refused g^2 = 1 - 4 Psi without solving, and "
	             "returned the last admissible c\n" );
}

/*
 * THE TARGET IS ASKED IN THE SOURCE'S Psi AND THE FAMILY IS LABELLED IN Psi_N,
 * AND GETTING IT BACKWARDS CONVERGES TO A REVERSED SHEAR.
 *
 * This is the same trap the_two_normalised_fluxes_run_in_opposite_directions
 * pins one level down, met where a CALLER meets it. The map owns the
 * reflection, so the assertion is that the target sees 1 - Psi_N and not
 * Psi_N -- checked by handing it a target that is not symmetric about 0.5, for
 * the reason that a symmetric one cannot tell the two apart.
 */
BOOST_AUTO_TEST_CASE( the_map_asks_its_target_in_the_sources_normalised_flux )
{
	// NOT SYMMETRIC ABOUT 0.5, WHICH IS THE WHOLE POINT. Labels of 0.10 and
	// 0.90 reflect to 0.90 and 0.10 -- the same SET -- so a check on the set
	// would pass either way round and test nothing. 0.10 and 0.40 reflect to
	// 0.90 and 0.60, which share no element with them.
	std::vector<double> const label = { 0.10, 0.40 };
	meq::FluxSurfaceFamily const family =
		familyOf( label, { 3.0, 3.0 }, { 0.5, 0.5 } );

	std::vector<double> asked;
	meq::ToroidalFieldMap::Options options;
	options.degree = 1;
	options.target = [ &asked ]( double psi )
	{
		asked.push_back( psi );
		return 1.0;
	};

	meq::ToroidalFieldMap map(
		[ & ]( std::vector<meq::Knot> const & ) { return &family; }, options );
	map( { 4.0, 0.0 } );

	BOOST_TEST_REQUIRE( asked.size() == label.size() );

	// THE PAIRING, IN ORDER, so a reflection applied to the wrong surface is
	// caught as well as no reflection at all.
	for ( std::size_t i = 0; i < label.size(); ++i )
		BOOST_TEST( asked[ i ] == 1.0 - label[ i ],
			boost::test_tools::tolerance( 1.0e-14 ) );

	// And the control: the un-reflected labels are NOT what was asked, which is
	// what says the assertion above has teeth on this fixture.
	for ( std::size_t i = 0; i < label.size(); ++i )
		BOOST_TEST( std::abs( asked[ i ] - label[ i ] ) > 0.1,
			"the target was asked at the FAMILY's label rather than the "
			"source's, so gg' comes back with its shear reversed" );

	std::printf( "  the target was asked at Psi = %.2f, %.2f for surfaces at "
	             "Psi_N = %.2f, %.2f\n",
	             asked[ 0 ], asked[ 1 ], label[ 0 ], label[ 1 ] );
}

/*
 * A MAP WITH NOTHING TO SOLVE WITH IS REFUSED AT CONSTRUCTION, because the
 * alternative is a null callable met for the first time inside KINSOL.
 */
BOOST_AUTO_TEST_CASE( a_map_without_a_solve_or_a_target_is_refused )
{
	meq::FluxSurfaceFamily const family =
		familyOf( { 0.5 }, { 3.0 }, { 0.5 } );
	auto solve = [ & ]( std::vector<meq::Knot> const & ) { return &family; };

	meq::ToroidalFieldMap::Options complete;
	complete.target = []( double ) { return 1.0; };

	BOOST_CHECK_THROW(
		meq::ToroidalFieldMap( nullptr, complete ), std::invalid_argument );

	meq::ToroidalFieldMap::Options noTarget;
	BOOST_CHECK_THROW(
		meq::ToroidalFieldMap( solve, noTarget ), std::invalid_argument );

	// A degree-zero g^2 is a constant, so gg' is identically zero and no
	// coefficient the loop moves can change the equilibrium.
	meq::ToroidalFieldMap::Options flat = complete;
	flat.degree = 0;
	BOOST_CHECK_THROW(
		meq::ToroidalFieldMap( solve, flat ), std::invalid_argument );

	// And a solve that cannot produce an equilibrium is a refusal rather than
	// a throw, since it is reachable from inside the Newton.
	meq::ToroidalFieldMap failing(
		[]( std::vector<meq::Knot> const & )
		{
			return static_cast<meq::FluxSurfaceFamily const *>( nullptr );
		},
		complete );
	std::vector<double> out;
	BOOST_CHECK_NO_THROW( out = failing( { 4.0, 0.0, 0.0 } ) );
	BOOST_TEST( failing.refusals() == 1 );
	for ( double value : out )
		BOOST_TEST( std::isfinite( value ) );
}
