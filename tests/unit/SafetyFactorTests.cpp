#define BOOST_TEST_MODULE SafetyFactorTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "meq/SafetyFactor.hpp"
#include "meq/SafetyFactorSolve.hpp"

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
 * Psi = 0.05 -- and SplineProfile continues linearly outside its knots, which
 * carries the innermost gg' and its slope out past the edge and keeps g^2
 * growing where the field should be the vacuum one. VacuumOutside puts a knot at Psi = 0 carrying zero
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
 * THE KNOTS FOLLOW THE MESH: about two per element across the minor radius,
 * a knot spacing from the axis and from the edge, and the labels are the
 * family's own to the bit.
 */
BOOST_AUTO_TEST_CASE( the_knots_follow_the_mesh_and_their_layers_shrink_with_it )
{
	// The q-driven box: 8 x 8 cells at level 0, four times that per level.
	std::size_t const level0 = meq::toroidalKnotCount( 64.0, 2 );
	std::size_t const level1 = meq::toroidalKnotCount( 256.0, 2 );
	std::size_t const level2 = meq::toroidalKnotCount( 1024.0, 2 );
	BOOST_TEST( level0 == 10u );
	BOOST_TEST( level1 == 19u );
	BOOST_TEST( level2 == 37u );
	// k = 3 is the cubic's own order: the same count.
	BOOST_TEST( meq::toroidalKnotCount( 256.0, 3 ) == level1 );
	// Above it the knots grow faster than 1/h.
	BOOST_TEST( meq::toroidalKnotCount( 256.0, 5 ) > level1 );
	BOOST_TEST( meq::toroidalKnotCount( 1.0, 2 ) == 6u );
	BOOST_CHECK_THROW( meq::toroidalKnotCount( 0.0, 2 ), std::invalid_argument );

	for ( std::size_t knots : { level0, level1, level2 } )
	{
		meq::ToroidalKnotCut const cut = meq::toroidalKnotCut( knots );
		std::vector<double> const labels = meq::toroidalKnotLabels( knots );
		double const spacing = 1.0/static_cast<double>( knots + 1 );
		BOOST_TEST( labels.size() == knots );
		BOOST_TEST( std::abs( std::sqrt( cut.innerCut ) - spacing ) < 1e-15 );
		BOOST_TEST( std::abs( 1.0 - std::sqrt( cut.outerCut ) - spacing ) < 1e-15,
			"the layer outside the outermost knot is not one knot spacing, so the "
			"part of the request no knot sees does not shrink with the mesh" );
		for ( std::size_t j = 0; j < knots; ++j )
			BOOST_TEST( std::abs( std::sqrt( labels[ j ] ) - ( j + 1.0 )*spacing ) < 1e-14 );
	}
}

/*
 * THE SPLINE: exact on a cubic in the source's Psi, gg' half its derivative,
 * linear in its values -- the cardinal functions sum to the spline -- and
 * reflected once: knots given in Psi_N, evaluated in Psi = 1 - Psi_N.
 */
BOOST_AUTO_TEST_CASE( the_spline_is_a_not_a_knot_cubic_in_the_sources_psi )
{
	std::vector<double> const labels = meq::toroidalKnotLabels( 9 );
	meq::ToroidalFieldSpline spline( labels );
	auto cubic = []( double psi ) { return 4.84 + 2.42*psi + 0.3025*psi*psi - 0.2*psi*psi*psi; };
	auto slope = []( double psi ) { return 2.42 + 0.605*psi - 0.6*psi*psi; };
	std::vector<double> values;
	for ( double label : labels )
		values.push_back( cubic( 1.0 - label ) );
	spline.setValues( values );

	double worst = 0.0, worstSlope = 0.0;
	for ( double psi = 0.0; psi <= 1.0; psi += 0.01 )
	{
		worst = std::max( worst, std::abs( spline( psi ) - cubic( psi ) ) );
		worstSlope = std::max( worstSlope, std::abs( spline( psi, 1 ) - slope( psi ) ) );
	}
	BOOST_TEST( worst < 1e-12, "a not-a-knot cubic reproduces a cubic, the end pieces included" );
	BOOST_TEST( worstSlope < 1e-11 );

	std::shared_ptr<meq::Profile const> const ggPrime = spline.ggPrime();
	for ( double psi : { 0.0, 0.3, 0.77, 1.0 } )
	{
		BOOST_TEST( std::abs( ( *ggPrime )( psi ) - 0.5*slope( psi ) ) < 1e-11,
			"gg' is half of d( g^2 )/dPsi in the SOURCE's Psi; the wrong sign is a "
			"reversed shear (section 2)" );
		double const h = 1e-6;
		BOOST_TEST( std::abs( ggPrime->prime( psi ) - ( ( *ggPrime )( psi + h ) - ( *ggPrime )( psi - h ) )/( 2*h ) ) < 1e-6 );
	}

	double sum = 0.0;
	for ( std::size_t j = 0; j < spline.size(); ++j )
		sum += values[ j ]*spline.cardinal( j, 0.4 );
	BOOST_TEST( std::abs( sum - spline( 0.4 ) ) < 1e-12 );
	BOOST_TEST( std::abs( spline.cardinal( 2, 1.0 - labels[ 2 ] ) - 1.0 ) < 1e-12 );
	BOOST_TEST( std::abs( spline.cardinal( 2, 1.0 - labels[ 3 ] ) ) < 1e-12 );

	BOOST_TEST( spline.positive() );
	values[ 0 ] = -1.0;
	spline.setValues( values );
	BOOST_TEST( !spline.positive() );
	BOOST_CHECK_THROW( meq::ToroidalFieldSpline( std::vector<double>{ 0.1, 0.2, 0.3 } ), std::invalid_argument );
	BOOST_CHECK_THROW( meq::ToroidalFieldSpline( std::vector<double>{ 0.1, 0.3, 0.2, 0.4 } ), std::invalid_argument );
}

namespace
{
	/*
	 * A SYNTHETIC q-DRIVEN SYSTEM: q/g at knot i is a fixed geometry times a
	 * response to the spline, q/g = a_i ( 1 + b sum_j u_j )^-1, standing in
	 * for the solve. Enough to check the rows, the pin and the totality.
	 */
	struct Synthetic
	{
		std::vector<double> labels = meq::toroidalKnotLabels( 6 );
		std::vector<double> a = { 1.0, 0.9, 0.8, 0.75, 0.7, 0.65 };
		double b = 0.01;
		int solves = 0;
		double askedPsi = -1.0;

		bool solve( meq::ToroidalFieldSpline const &s, meq::SafetyFactorEvaluation &out )
		{
			++solves;
			double sum = 0.0;
			for ( double v : s.values() )
				sum += v;
			out.qOverG.clear();
			for ( double ai : a )
				out.qOverG.push_back( ai/( 1.0 + b*sum ) );
			out.measured = sum;
			return true;
		}

		bool differentiate( meq::ToroidalFieldSpline const &s, meq::SafetyFactorDerivative &out )
		{
			std::size_t const m = s.size();
			double sum = 0.0;
			for ( double v : s.values() )
				sum += v;
			out.qOverG.assign( m*m, 0.0 );
			for ( std::size_t i = 0; i < m; ++i )
				for ( std::size_t j = 0; j < m; ++j )
					out.qOverG[ i*m + j ] = -a[ i ]*b/( ( 1.0 + b*sum )*( 1.0 + b*sum ) );
			out.measured.assign( m, 1.0 );
			return true;
		}

		meq::SafetyFactorSystem system( meq::ScalePin kind, double value, std::function<double( double )> target )
		{
			meq::SafetyFactorSystem::Pin pin;
			pin.kind = kind;
			pin.value = value;
			return meq::SafetyFactorSystem(
				labels, target, pin,
				[ this ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorEvaluation &o ) { return solve( s, o ); },
				[ this ]( meq::ToroidalFieldSpline const &s, meq::SafetyFactorDerivative &o ) { return differentiate( s, o ); },
				1.0e-12 );
		}
	};
}

/*
 * THE SYSTEM IS TOTAL: a spline with g^2 <= 0 anywhere on [ 0, 1 ] is refused
 * without a solve -- it is not an equilibrium -- and so is a solve that fails
 * or throws. Refused means residual() returns false and the least squares
 * rejects the step; nothing throws for a reason the loop can meet.
 */
BOOST_AUTO_TEST_CASE( the_system_is_total_where_g_squared_goes_negative )
{
	Synthetic synthetic;
	meq::SafetyFactorSystem system = synthetic.system( meq::ScalePin::EdgeField, 2.2,
	                                                   []( double ) { return 1.5; } );
	std::vector<double> r;
	std::vector<double> u( 6, 4.84 );
	BOOST_TEST( system.residual( u, r ) );
	BOOST_TEST( r.size() == 7u );

	std::vector<double> negative( u );
	negative[ 3 ] = -0.5;
	int const before = synthetic.solves;
	BOOST_TEST( !system.residual( negative, r ) );
	BOOST_TEST( synthetic.solves == before, "g^2 < 0 reached the solve" );
	BOOST_TEST( system.refusals() == 1 );

	// A spline positive at every knot that dips below zero between them, on
	// the continued end piece, is refused too: the source reads all of [ 0, 1 ].
	std::vector<double> dipping = { 5.0, 0.2, 5.0, 0.2, 5.0, 0.05 };
	BOOST_TEST( !system.residual( dipping, r ) );

	meq::SafetyFactorSystem failing(
		synthetic.labels, []( double ) { return 1.5; }, { meq::ScalePin::EdgeField, 2.2 },
		[]( meq::ToroidalFieldSpline const &, meq::SafetyFactorEvaluation & ) -> bool
		{
			throw std::runtime_error( "the inner solve did not converge" );
		},
		[]( meq::ToroidalFieldSpline const &, meq::SafetyFactorDerivative & ) { return false; },
		1.0e-12 );
	BOOST_TEST( !failing.residual( u, r ), "a throwing solve must be a refusal, not an exception" );
}

/*
 * THE TARGET IS ASKED IN THE SOURCE'S Psi: the knots are labelled by the
 * family's Psi_N, zero on the axis, and the target is read at 1 - Psi_N. And
 * the rows are what the header says: q/q* - 1, and the pin row weighted.
 */
BOOST_AUTO_TEST_CASE( the_rows_ask_the_target_in_the_sources_psi )
{
	Synthetic synthetic;
	auto target = []( double psi ) { return 1.0 + psi; };
	meq::SafetyFactorSystem system = synthetic.system( meq::ScalePin::EdgeField, 2.2, target );
	for ( std::size_t i = 0; i < synthetic.labels.size(); ++i )
		BOOST_TEST( std::abs( system.targets()[ i ] - ( 2.0 - synthetic.labels[ i ] ) ) < 1e-15,
			"the target at Psi_N = " << synthetic.labels[ i ] << " was read at the wrong "
			"end of the profile; a reflection missed converges to a reversed shear" );

	std::vector<double> u( 6, 4.84 ), r;
	BOOST_TEST_REQUIRE( system.residual( u, r ) );
	double const qOverG0 = synthetic.a[ 0 ]/( 1.0 + synthetic.b*6.0*4.84 );
	BOOST_TEST( std::abs( r[ 0 ] - ( qOverG0*2.2/system.targets()[ 0 ] - 1.0 ) ) < 1e-14 );
	// R0B0's row: the spline continued to Psi = 0, here the constant 4.84.
	BOOST_TEST( std::abs( r[ 6 ] ) < 1e-14 );
	u.assign( 6, 4.84*1.01 );
	BOOST_TEST_REQUIRE( system.residual( u, r ) );
	BOOST_TEST( std::abs( r[ 6 ] - 0.01 ) < 1e-9 );

	std::vector<double> j;
	BOOST_TEST_REQUIRE( system.jacobian( u, j ) );
	BOOST_TEST( j.size() == 7u*6u );
	// The R0B0 row is the cardinal functions at the edge; they sum to one.
	double sum = 0.0;
	for ( std::size_t k = 0; k < 6; ++k )
		sum += j[ 6*6 + k ];
	BOOST_TEST( std::abs( sum - 1.0/( 2.2*2.2 ) ) < 1e-9 );
	// Asked anywhere but the last admitted point, the Jacobian refuses loudly.
	std::vector<double> elsewhere( 6, 5.0 );
	BOOST_CHECK_THROW( system.jacobian( elsewhere, j ), std::logic_error );

	// A pin row that reads the solve: ( m - m* )/m*, weighted.
	meq::SafetyFactorSystem current = synthetic.system( meq::ScalePin::PlasmaCurrent, 30.0, target );
	BOOST_TEST_REQUIRE( current.residual( u, r ) );
	BOOST_TEST( std::abs( r[ 6 ] - ( 6.0*4.84*1.01 - 30.0 )/30.0 ) < 1e-9 );
	BOOST_TEST( current.floor().size() == 7u );
	BOOST_CHECK_THROW( synthetic.system( meq::ScalePin::EdgeField, -2.2, target ), std::invalid_argument );
	BOOST_CHECK_THROW( synthetic.system( meq::ScalePin::EdgeField, 2.2, []( double ) { return 0.0; } ),
	                   std::invalid_argument );
}

namespace
{
	using Residual = std::function<bool( std::vector<double> const &, std::vector<double> & )>;

	/// An affine least squares with a known minimiser: two exact rows and one
	/// inconsistent third, r = A x - b.
	struct Affine
	{
		std::vector<double> a = { 2.0, 0.5, -1.0, 3.0, 1.0, 1.0 };   // 3 x 2
		std::vector<double> b = { 1.0, 2.0, 2.5 };
		int calls = 0;
		bool operator()( std::vector<double> const &x, std::vector<double> &r )
		{
			++calls;
			r.resize( 3 );
			for ( int i = 0; i < 3; ++i )
				r[ i ] = a[ 2*i ]*x[ 0 ] + a[ 2*i + 1 ]*x[ 1 ] - b[ i ];
			return true;
		}
	};
}

/*
 * GAUSS-NEWTON IS EXACT ON AN AFFINE PROBLEM: one step to the least-squares
 * minimiser, whose misfit is not zero -- and the loop is converged there,
 * because what a step can remove is zero, not because the residual is.
 */
BOOST_AUTO_TEST_CASE( gauss_newton_finishes_an_affine_least_squares_in_one_step )
{
	auto affine = std::make_shared<Affine>();
	Residual residual = [ affine ]( std::vector<double> const &x, std::vector<double> &r ) { return ( *affine )( x, r ); };
	auto jacobian = [ affine ]( std::vector<double> const &, std::vector<double> &j )
	{
		j = affine->a;
		return true;
	};
	meq::LeastSquaresResult const result = meq::solveLeastSquares( { 10.0, -7.0 }, residual, jacobian );
	BOOST_TEST( result.converged, result.status );
	BOOST_TEST( result.iterations == 1 );
	BOOST_TEST( result.misfit > 0.1, "the third row is inconsistent; its misfit stays" );
	// The normal equations' answer, checked independently.
	double const ata[ 4 ] = { 4 + 1 + 1, 1 - 3 + 1, 1 - 3 + 1, 0.25 + 9 + 1 };
	double const atb[ 2 ] = { 2 - 2 + 2.5, 0.5 + 6 + 2.5 };
	double const det = ata[ 0 ]*ata[ 3 ] - ata[ 1 ]*ata[ 2 ];
	double const x0 = ( atb[ 0 ]*ata[ 3 ] - ata[ 1 ]*atb[ 1 ] )/det;
	double const x1 = ( ata[ 0 ]*atb[ 1 ] - ata[ 2 ]*atb[ 0 ] )/det;
	BOOST_TEST( std::abs( result.x[ 0 ] - x0 ) < 1e-12 );
	BOOST_TEST( std::abs( result.x[ 1 ] - x1 ) < 1e-12 );
}

/*
 * WITH A CONSTRAINT THE SEARCH JUDGES WHAT THE STEP CAN REMOVE, NOT || r ||^2
 * (M-QS-2). The answer is least squares within a level set of ln u_1 + ln u_2,
 * and the pin u_1 = 1.5 picks the level. It is not the minimum of || r ||^2 on
 * the pin's line. Started at that minimum, (1.5, 2.25), every step toward the
 * answer, u_2 = ( 3 + sqrt 3 )/2, raises || r ||^2. A search on || r ||^2
 * refuses them all and stalls at the start.
 */
BOOST_AUTO_TEST_CASE( a_constrained_step_is_judged_on_what_it_can_remove )
{
	Residual residual = []( std::vector<double> const &u, std::vector<double> &r )
	{
		r = { u[ 0 ] - 1.0, u[ 1 ] - 2.0, u[ 0 ] + u[ 1 ] - 4.0, u[ 0 ] - 1.5 };
		return true;
	};
	auto jacobian = []( std::vector<double> const &, std::vector<double> &j )
	{
		j = { 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0 };
		return true;
	};
	meq::LeastSquaresOptions options;
	options.tolerance = 1.0e-12;
	options.constraints = 1;
	options.level = []( std::vector<double> const &u )
	{
		return std::vector<double>{ 1.0/u[ 0 ], 1.0/u[ 1 ] };
	};
	meq::LeastSquaresResult const result = meq::solveLeastSquares( { 1.5, 2.25 }, residual, jacobian, options );
	BOOST_TEST( result.converged, result.status );
	BOOST_TEST( std::abs( result.x[ 0 ] - 1.5 ) < 1e-12 );
	BOOST_TEST( std::abs( result.x[ 1 ] - 0.5*( 3.0 + std::sqrt( 3.0 ) ) ) < 1e-10,
		"u_2 = " << result.x[ 1 ] << ": not least squares within the pin's level set" );
}

/*
 * A RANK-DEFICIENT JACOBIAN IS REPORTED, NOT STEPPED THROUGH; A THROW ENDS THE
 * LOOP AND IS REPORTED, NOT RETHROWN; A REFUSED TRIAL IS A REJECTED STEP.
 */
BOOST_AUTO_TEST_CASE( an_undetermined_jacobian_a_throw_and_a_refusal )
{
	Residual flat = []( std::vector<double> const &x, std::vector<double> &r )
	{
		r = { x[ 0 ] + x[ 1 ] - 1.0, 2.0*( x[ 0 ] + x[ 1 ] ) + 3.0 };
		return true;
	};
	auto singular = []( std::vector<double> const &, std::vector<double> &j )
	{
		j = { 1.0, 1.0, 2.0, 2.0 };
		return true;
	};
	meq::LeastSquaresResult const undetermined = meq::solveLeastSquares( { 0.0, 0.0 }, flat, singular );
	BOOST_TEST( !undetermined.converged );
	BOOST_TEST( undetermined.status.find( "RANK DEFICIENT" ) != std::string::npos, undetermined.status );

	Residual throwing = []( std::vector<double> const &x, std::vector<double> &r ) -> bool
	{
		if ( x[ 0 ] > 0.5 )
			throw std::runtime_error( "an inner solve exploded" );
		r = { x[ 0 ] - 1.0 };
		return true;
	};
	auto unit = []( std::vector<double> const &, std::vector<double> &j )
	{
		j = { 1.0 };
		return true;
	};
	meq::LeastSquaresResult const thrown = meq::solveLeastSquares( { 0.0 }, throwing, unit );
	BOOST_TEST( !thrown.converged );
	BOOST_TEST( thrown.status.find( "an inner solve exploded" ) != std::string::npos, thrown.status );

	// A residual refused beyond x = 0.7: the full step to 1 is rejected and
	// the search halves into the admissible side, then onward.
	int refusals = 0;
	Residual walled = [ &refusals ]( std::vector<double> const &x, std::vector<double> &r )
	{
		if ( x[ 0 ] > 0.7 && x[ 0 ] < 0.99 )
		{
			++refusals;
			return false;
		}
		r = { std::atan( x[ 0 ] - 1.0 ) };
		return true;
	};
	auto slope = []( std::vector<double> const &x, std::vector<double> &j )
	{
		j = { 1.0/( 1.0 + ( x[ 0 ] - 1.0 )*( x[ 0 ] - 1.0 ) ) };
		return true;
	};
	meq::LeastSquaresResult const around = meq::solveLeastSquares( { 0.0 }, walled, slope );
	BOOST_TEST( around.rejected >= 0 );
	BOOST_TEST( around.converged, around.status );
}

namespace
{
	/// A quadratic problem whose residual carries a fresh draw of noise
	/// every evaluation, as a warm-started solve's does.
	struct Noisy
	{
		double noise;
		int calls = 0;
		bool operator()( std::vector<double> const &x, std::vector<double> &r )
		{
			++calls;
			r.resize( x.size() );
			for ( std::size_t i = 0; i < x.size(); ++i )
			{
				double const away = x[ i ] - 1.0 - 0.1*static_cast<double>( i );
				r[ i ] = 1.8*away + 0.3*away*away
				         + noise*std::sin( 12.9898*calls + 78.233*static_cast<double>( i ) );
			}
			return true;
		}
	};
}

/*
 * CONVERGED MEANS THE REMOVABLE RESIDUAL REACHED max( TOLERANCE, FLOOR ), and
 * the floor is measured where it governs (M-363): a stall above a floor
 * understated a hundredfold is measured and then converged at the noise; a
 * stall above the noise is not converged.
 */
BOOST_AUTO_TEST_CASE( a_stall_at_the_measured_noise_is_converged_and_above_it_is_not )
{
	double const noise = 1.0e-4;
	auto noisy = std::make_shared<Noisy>( Noisy{ noise } );
	Residual residual = [ noisy ]( std::vector<double> const &x, std::vector<double> &r ) { return ( *noisy )( x, r ); };
	auto jacobian = []( std::vector<double> const &x, std::vector<double> &j )
	{
		std::size_t const n = x.size();
		j.assign( n*n, 0.0 );
		for ( std::size_t i = 0; i < n; ++i )
			j[ i*n + i ] = 1.8 + 0.6*( x[ i ] - 1.0 - 0.1*static_cast<double>( i ) );
		return true;
	};
	meq::LeastSquaresOptions understated;
	understated.tolerance = 1.0e-12;
	understated.floor = [ noise ]() { return std::vector<double>( 3, noise/100.0 ); };
	meq::LeastSquaresResult const measured = meq::solveLeastSquares( { 3.0, 2.0, 0.0 }, residual, jacobian, understated );
	std::printf( "\n  A NOISY RESIDUAL, FLOOR UNDERSTATED: %s, floor %.2e\n", measured.status.c_str(), measured.floor );
	BOOST_TEST( measured.converged, measured.status );
	BOOST_TEST( measured.floor >= noise, "the floor judged against is under the residual's own noise" );
	for ( std::size_t i = 0; i < 3; ++i )
		BOOST_TEST( std::abs( measured.x[ i ] - 1.0 - 0.1*static_cast<double>( i ) ) < 20.0*noise );

	// The control: noise far above any floor and no measurement can excuse a
	// stall at 1e-2 -- a residual that does not fall below its noise is a
	// stall, and the tolerance is far under both.
	auto loud = std::make_shared<Noisy>( Noisy{ 1.0e-2 } );
	Residual loudResidual = [ loud ]( std::vector<double> const &x, std::vector<double> &r ) { return ( *loud )( x, r ); };
	meq::LeastSquaresOptions honest;
	honest.tolerance = 1.0e-12;
	honest.maxIterations = 4;
	meq::LeastSquaresResult const stalled = meq::solveLeastSquares( { 3.0, 2.0, 0.0 }, loudResidual, jacobian, honest );
	BOOST_TEST( !stalled.converged, "a residual still moving at 1e-2 after four steps reads converged: "
	                                << stalled.status );
}
