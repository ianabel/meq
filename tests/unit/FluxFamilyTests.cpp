// Unit tests for the flux-surface family, its interpolation, and the per-psi
// geometry cache: src/meq/FluxFamily.
//
// INVERSION-PLAN.md stage IN-6. MFEM-FREE, like the unit under test, and that
// is the whole reason the split exists: the cache's invalidation rule IS
// MANTA-COUPLING.md section 8 -- a served answer must equal a cold one BIT FOR
// BIT and a stale cache must be impossible rather than unlikely -- and that is
// the part of the stage with a stated obligation attached. CI cannot obtain the
// MFEM branch meq needs, so this is where the contract is gated.
//
// THE EXTRACTOR HERE IS SYNTHETIC AND THAT IS THE POINT. A real extraction
// traces contours through a solved field and takes a quarter of a second on an
// ordinary mesh; what is being tested is not the geometry but WHEN the geometry
// is recomputed, so the extractor is a counter that manufactures a family from
// the psi vector it is handed. That makes "was this served or recomputed?" an
// exact question with an exact answer instead of a timing.
//
// tests/convergence/FluxGridConvergence.cpp is the other half: the same cache
// driven by the real extraction over a real solve, plus the geometry itself.

#define BOOST_TEST_MODULE FluxFamilyTests
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "meq/FluxFamily.hpp"

namespace
{
	double const twoPi = 6.283185307179586476925286766559;

	/// A family whose surfaces are concentric circles of minor radius
	/// a( rho ) = 0.3 rho about ( 1, 0 ), carrying quantities that are simple
	/// analytic functions of the label. Enough structure to interpolate and to
	/// compare, and no field anywhere.
	///
	/// @param mark a value folded into every quantity so that two families built
	///        from the same psi but a DIFFERENT extractor state are tellable
	///        apart. That is what the resetForRun() case turns on.
	meq::FluxSurfaceFamily circles( std::size_t count, std::size_t angles,
	                                double innerCut, double outerCut,
	                                double mark = 0.0 )
	{
		meq::FluxSurfaceFamily family;
		family.axisR = 1.0;
		family.axisZ = 0.0;
		family.psiAxis = 1.0;
		family.psiBoundary = 0.0;
		family.innerCut = innerCut;
		family.outerCut = outerCut;
		family.angles = angles;

		double const inner = std::sqrt( innerCut );
		double const outer = std::sqrt( outerCut );
		double const last = static_cast<double>( count - 1 );

		for ( std::size_t i = 0; i < count; ++i )
		{
			double const label = inner
				+ ( outer - inner )*static_cast<double>( i )/last;

			meq::FluxSurface surface;
			surface.radial = label;
			surface.normalisedFlux = label*label;
			surface.level = family.psiAxis
				- surface.normalisedFlux*( family.psiAxis - family.psiBoundary );

			double const minor = 0.3*label;
			surface.radius.resize( angles );
			surface.z.resize( angles );
			surface.extended.assign( angles, 0 );
			for ( std::size_t j = 0; j < angles; ++j )
			{
				double const theta = twoPi*static_cast<double>( j )
					/static_cast<double>( angles );
				surface.radius[ j ] = family.axisR + minor*std::cos( theta );
				surface.z[ j ] = family.axisZ + minor*std::sin( theta );
			}

			// Chosen so that every one is a DIFFERENT function of the label: a
			// bug that served the wrong surface would otherwise be invisible in
			// a quantity that happened to be flat.
			surface.vPrime = ( 1.0 + mark )*twoPi*twoPi*family.axisR*minor;
			surface.volume = ( 1.0 + mark )*twoPi*family.axisR
				*3.141592653589793*minor*minor;
			surface.arcLength = ( 1.0 + mark )*twoPi*minor;
			surface.inverseRSquared = ( 1.0 + mark )/( family.axisR*family.axisR );
			surface.absGradPsi = ( 1.0 + mark )*( 0.5 + label );
			surface.gradPsiSquared = ( 1.0 + mark )*( 0.5 + label )
				*( 0.5 + label );
			family.surfaces.push_back( surface );
		}

		return family;
	}

	/// An extractor that counts its calls and folds the psi it is handed into
	/// the answer, so that "was this recomputed?" and "was it recomputed from
	/// THIS psi?" are separate, checkable questions.
	class CountingExtractor
	{
		public:
			meq::FluxSurfaceFamily operator()( std::vector<double> const &psi )
			{
				++calls;
				if ( refuse )
					throw std::runtime_error(
						"CountingExtractor: refusing, as asked" );

				double total = 0.0;
				for ( double value : psi )
					total += value;

				meq::FluxSurfaceFamily family =
					circles( 9, 32, 0.05, 0.95, generation );
				for ( meq::FluxSurface &surface : family.surfaces )
					surface.vPrime *= total;
				return family;
			}

			int calls = 0;
			bool refuse = false;
			double generation = 0.0;
	};

	std::vector<double> ramp( std::size_t n, double scale = 1.0 )
	{
		std::vector<double> psi( n, 0.0 );
		for ( std::size_t i = 0; i < n; ++i )
			psi[ i ] = scale*( 1.0 + 0.01*static_cast<double>( i ) );
		return psi;
	}

	meq::SurfaceQuantity vPrimeOf()
	{
		return []( meq::FluxSurface const &s ) { return s.vPrime; };
	}
}

// ---------------------------------------------------------------------------
// The label
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( the_normalised_flux_convention_is_zero_on_the_axis )
{
	double const psiAxis = 0.35;
	double const psiBoundary = -0.05;

	BOOST_TEST( meq::normalisedFlux( psiAxis, psiAxis, psiBoundary ) == 0.0 );
	BOOST_TEST( meq::normalisedFlux( psiBoundary, psiAxis, psiBoundary ) == 1.0 );

	// Round trip, at a value that is neither end.
	double const level = meq::fluxAtNormalised( 0.37, psiAxis, psiBoundary );
	BOOST_TEST( meq::normalisedFlux( level, psiAxis, psiBoundary )
	            == 0.37, boost::test_tools::tolerance( 1.0e-15 ) );

	// A zero span is not a normalisation and is refused rather than dividing.
	BOOST_CHECK_THROW( meq::normalisedFlux( 0.1, 0.2, 0.2 ),
	                   std::invalid_argument );
	BOOST_CHECK_THROW( meq::fluxAtNormalised( 0.1, 0.2, 0.2 ),
	                   std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( the_label_derivative_is_analytic_and_refuses_the_axis )
{
	meq::FluxSurfaceFamily const family = circles( 5, 16, 0.05, 0.95 );

	// dPsi_N/dpsi = -1/( psi_ax - psi_bnd ), exactly.
	BOOST_TEST( family.normalisedFluxDerivative() == -1.0 );

	// drho/dpsi = ( dPsi_N/dpsi )/( 2 rho ), exactly, and NOT a difference.
	double const label = 0.5;
	BOOST_TEST( family.radialDerivative( label ) == -1.0/( 2.0*label ) );

	// rho = 0 is the magnetic axis, where the coordinate itself is singular.
	// A refusal, not an infinity handed to a consumer.
	BOOST_CHECK_THROW( family.radialDerivative( 0.0 ), std::runtime_error );
	BOOST_CHECK_THROW( family.radialDerivative( -0.1 ), std::runtime_error );
}

// ---------------------------------------------------------------------------
// The interpolation
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( the_interpolation_is_exact_at_its_own_nodes )
{
	meq::FluxSurfaceFamily const family = circles( 11, 16, 0.04, 0.96 );

	// AT THE NODES IT MUST BE THE NODE, TO THE LAST BIT. An interpolant that
	// merely comes close at its own data is not an interpolant, and a consumer
	// whose grid happens to coincide with meq's surfaces would be reading a
	// smoothed version of the answer without anything saying so.
	double worst = 0.0;
	for ( meq::FluxSurface const &surface : family.surfaces )
	{
		double const served = family.at( surface.radial, vPrimeOf() );
		worst = std::max( worst, std::abs( served - surface.vPrime ) );
	}
	std::printf( "  interpolation at its own nodes: worst %.3e\n", worst );
	BOOST_TEST( worst == 0.0,
	            "the monotone cubic does not reproduce its own data" );
}

BOOST_AUTO_TEST_CASE( the_interpolation_reproduces_a_straight_line )
{
	meq::FluxSurfaceFamily family = circles( 9, 16, 0.05, 0.95 );
	for ( meq::FluxSurface &surface : family.surfaces )
		surface.vPrime = 3.0 + 2.0*surface.radial;

	double worst = 0.0;
	for ( int i = 0; i <= 200; ++i )
	{
		double const t = static_cast<double>( i )/200.0;
		double const label = family.innerLabel()
			+ t*( family.outerLabel() - family.innerLabel() );
		double const exact = 3.0 + 2.0*label;
		worst = std::max( worst, std::abs( family.at( label, vPrimeOf() )
		                                   - exact ) );
	}

	std::printf( "  a straight line through the monotone cubic: worst %.3e\n",
	             worst );
	BOOST_TEST( worst < 1.0e-14,
	            "the monotone cubic does not reproduce a linear function" );
}

BOOST_AUTO_TEST_CASE( the_interpolation_cannot_overshoot_into_a_negative_vprime )
{
	// THE PROPERTY THE MONOTONE LIMITER EXISTS FOR, and it is not a nicety.
	// V' -> 0 at the magnetic axis, so an ordinary cubic that overshoots below
	// the innermost node hands a consumer a NEGATIVE volume derivative -- a
	// nonsense value arriving through exactly the route MANTA-COUPLING.md
	// section 8 forbids, and one no transport code would flag.
	//
	// The data here is deliberately nastier than a real family: V' rising
	// steeply from near zero, with a flat shoulder, which is where an unlimited
	// cubic rings.
	meq::FluxSurfaceFamily family = circles( 7, 16, 0.01, 0.99 );
	double const shoulder[ 7 ] = { 1.0e-6, 2.0e-3, 4.0e-1, 1.0, 1.02, 1.03,
	                               1.031 };
	for ( std::size_t i = 0; i < family.size(); ++i )
		family.surfaces[ i ].vPrime = shoulder[ i ];

	double lowest = shoulder[ 0 ];
	double highest = shoulder[ 6 ];
	double worstUnder = 0.0;
	double worstOver = 0.0;
	for ( int i = 0; i <= 4000; ++i )
	{
		double const t = static_cast<double>( i )/4000.0;
		double const label = family.innerLabel()
			+ t*( family.outerLabel() - family.innerLabel() );
		double const value = family.at( label, vPrimeOf() );
		worstUnder = std::max( worstUnder, lowest - value );
		worstOver = std::max( worstOver, value - highest );
	}

	std::printf( "  overshoot below the data %.3e, above it %.3e\n",
	             worstUnder, worstOver );
	BOOST_TEST( worstUnder <= 0.0,
	            "the interpolant went below its own smallest datum, which for V' "
	            "near the axis is a negative volume derivative" );
	BOOST_TEST( worstOver <= 0.0,
	            "the interpolant went above its own largest datum" );

	// And nothing anywhere in range is negative, which is the statement a
	// consumer actually depends on.
	for ( int i = 0; i <= 4000; ++i )
	{
		double const t = static_cast<double>( i )/4000.0;
		double const label = family.innerLabel()
			+ t*( family.outerLabel() - family.innerLabel() );
		BOOST_TEST_REQUIRE( family.at( label, vPrimeOf() ) > 0.0 );
	}
}

BOOST_AUTO_TEST_CASE( a_query_outside_the_cut_is_refused_and_not_extrapolated )
{
	// THE SEPARATRIX CUT, AND THE DECISION IT ENCODES. FreeGS extrapolates
	// outside [ 0.01, 0.99 ]; meq refuses. MANTA-COUPLING.md section 8: a state
	// meq cannot evaluate at must THROW, so the consumer's integrator treats it
	// as recoverable and retries with a smaller step. An extrapolated V' past
	// the separatrix is precisely the plausible-looking nonsense value that
	// converts a recoverable step into a wrong answer.
	meq::FluxSurfaceFamily const family = circles( 9, 16, 0.05, 0.95 );

	BOOST_TEST( family.covers( family.innerLabel() ) );
	BOOST_TEST( family.covers( family.outerLabel() ) );
	BOOST_TEST( !family.covers( 0.0 ) );
	BOOST_TEST( !family.covers( 1.0 ) );

	BOOST_CHECK_THROW( family.at( 0.0, vPrimeOf() ), std::domain_error );
	BOOST_CHECK_THROW( family.at( 1.0, vPrimeOf() ), std::domain_error );
	BOOST_CHECK_THROW( family.atNormalisedFlux( 0.999, vPrimeOf() ),
	                   std::domain_error );

	// A negative normalised flux has no square root and names no surface.
	BOOST_CHECK_THROW( family.atNormalisedFlux( -0.1, vPrimeOf() ),
	                   std::domain_error );

	// And the message has to say what the cut IS, or the caller cannot act on
	// it. Checked rather than assumed, because an unactionable refusal is only
	// marginally better than a wrong number.
	bool named = false;
	try
	{
		family.at( 1.0, vPrimeOf() );
	}
	catch ( std::domain_error const &error )
	{
		std::string const text( error.what() );
		named = text.find( "0.95" ) != std::string::npos
			&& text.find( "Psi_N" ) != std::string::npos;
		std::printf( "  refusal: %s\n", text.c_str() );
	}
	BOOST_TEST( named,
	            "the refusal does not say what the cut is, so a caller cannot "
	            "act on it" );
}

BOOST_AUTO_TEST_CASE( a_family_too_small_to_interpolate_refuses )
{
	meq::FluxSurfaceFamily empty;
	BOOST_CHECK_THROW( empty.innerLabel(), std::runtime_error );
	BOOST_CHECK_THROW( empty.outerLabel(), std::runtime_error );
	BOOST_TEST( !empty.covers( 0.5 ) );
	BOOST_CHECK_THROW( empty.at( 0.5, vPrimeOf() ), std::runtime_error );

	meq::FluxSurfaceFamily one = circles( 9, 16, 0.05, 0.95 );
	one.surfaces.resize( 1 );
	BOOST_CHECK_THROW( one.at( one.innerLabel(), vPrimeOf() ),
	                   std::runtime_error );

	// An empty quantity is a caller error and is named as one rather than
	// crashing inside std::function.
	meq::FluxSurfaceFamily const family = circles( 9, 16, 0.05, 0.95 );
	BOOST_CHECK_THROW( family.at( 0.5, meq::SurfaceQuantity() ),
	                   std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( the_label_metric_uses_the_analytic_chain_factor )
{
	meq::FluxSurfaceFamily const family = circles( 9, 16, 0.05, 0.95 );

	double const label = 0.6;
	double const chain = family.radialDerivative( label );

	BOOST_TEST( family.absGradLabelAt( label )
	            == std::abs( chain )*family.at( label,
	                []( meq::FluxSurface const &s ) { return s.absGradPsi; } ),
	            boost::test_tools::tolerance( 1.0e-15 ) );

	BOOST_TEST( family.gradLabelSquaredAt( label )
	            == chain*chain*family.at( label,
	                []( meq::FluxSurface const &s )
	                { return s.gradPsiSquared; } ),
	            boost::test_tools::tolerance( 1.0e-15 ) );
}

// ---------------------------------------------------------------------------
// The cache
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( a_served_answer_equals_a_cold_one_bit_for_bit )
{
	// MANTA-COUPLING.md section 8, and the tolerance is deliberately zero: MaNTA
	// pins the analogous property with a test asserting that a reused solver
	// matches a fresh one BIT FOR BIT, because the last defect of this kind left
	// the second run completing, plausible, and wrong in the eleventh digit.
	CountingExtractor cold;
	CountingExtractor warm;

	std::vector<double> const psi = ramp( 500 );

	meq::GeometryCache coldCache( std::ref( cold ) );
	meq::GeometryCache warmCache( std::ref( warm ) );

	std::size_t const nodes = 60;
	double worst = 0.0;
	for ( std::size_t residual = 0; residual < 4; ++residual )
	{
		for ( std::size_t node = 0; node < nodes; ++node )
		{
			double const label = 0.25 + 0.7*static_cast<double>( node )
				/static_cast<double>( nodes );

			// Cold: a fresh cache every single query, which is the naive
			// pointwise pattern MANTA-COUPLING.md section 5 describes.
			meq::GeometryCache once( std::ref( cold ) );
			double const fresh = once.at( psi, label, vPrimeOf() );

			double const served = warmCache.at( psi, label, vPrimeOf() );
			worst = std::max( worst, std::abs( served - fresh ) );
		}
	}

	std::printf( "  cold extractions %d, warm extractions %d over %zu queries\n",
	             cold.calls, warm.calls, warmCache.queries() );
	std::printf( "  served minus cold: %.3e\n", worst );

	BOOST_TEST( worst == 0.0,
	            "a served answer differs from a cold one, so the cache is not "
	            "transparent" );

	// THE COUNT IS THE PROMISE, AND IT IS A COUNT RATHER THAN A TIMING. What the
	// cache guarantees is one extraction per distinct psi; IN-P's 5.1x is what
	// that is worth in seconds on one machine, and belongs in
	// tests/performance.
	BOOST_TEST( warm.calls == 1 );
	BOOST_TEST( warmCache.queries() == 4*nodes );
	BOOST_TEST( warmCache.hits() == 4*nodes - 1 );
	BOOST_TEST( cold.calls == static_cast<int>( 4*nodes ) );
}

BOOST_AUTO_TEST_CASE( one_bit_of_psi_invalidates_the_cache )
{
	// The key is the psi vector compared BITWISE, not a hash: section 8 asks for
	// a stale cache to be impossible and a hash makes it unlikely. A single ulp
	// in one entry of five hundred has to miss.
	CountingExtractor extractor;
	meq::GeometryCache cache( std::ref( extractor ) );

	std::vector<double> psi = ramp( 500 );
	double const first = cache.at( psi, 0.5, vPrimeOf() );
	BOOST_TEST( extractor.calls == 1 );

	std::vector<double> nudged = psi;
	nudged[ 371 ] = std::nextafter( nudged[ 371 ],
	                                std::numeric_limits<double>::infinity() );
	double const second = cache.at( nudged, 0.5, vPrimeOf() );

	std::printf( "  one ulp in entry 371 of 500: extractions %d, "
	             "answer moved %.3e\n",
	             extractor.calls, std::abs( second - first ) );

	BOOST_TEST( extractor.calls == 2,
	            "a one-ulp change in psi was served from the cache" );

	// THE ANSWER DOES NOT MOVE, AND ASSERTING THAT IT WOULD WAS THE FIRST
	// VERSION OF THIS CASE. A one-ulp change in one entry of five hundred is
	// below the resolution of anything the extractor computes from them, so the
	// family it builds is bit-identical -- which is a statement about the
	// extractor and not about the cache. What the cache promises is that it
	// LOOKED AGAIN, and the count is that promise.
	//
	// The guard against a vacuous test is separate and is below: a materially
	// different psi has to move the answer, or the count would be measuring
	// nothing.
	BOOST_TEST( second == first );

	std::vector<double> const elsewhere = ramp( 500, 2.0 );
	double const third = cache.at( elsewhere, 0.5, vPrimeOf() );
	BOOST_TEST( extractor.calls == 3 );
	BOOST_TEST( third != first,
	            "the extractor's answer does not depend on psi at all, so the "
	            "extraction counts above are measuring nothing" );

	// And the size alone is enough to miss, which is the cheap half of the key.
	nudged.push_back( 1.0 );
	cache.at( nudged, 0.5, vPrimeOf() );
	BOOST_TEST( extractor.calls == 4 );
}

BOOST_AUTO_TEST_CASE( reset_for_run_forces_a_recomputation_the_key_cannot_see )
{
	// THE SUBTLE HALF OF SECTION 8. The bitwise key catches a changed psi. It
	// does NOT catch a changed EXTRACTOR -- a second run whose psi happens to
	// start where the first ended, against a different mesh or a different
	// solver, matches the key and is served the first run's geometry. That is
	// "a cache keyed on the object rather than the run", and resetForRun() is
	// the only thing that can see it.
	CountingExtractor extractor;
	meq::GeometryCache cache( std::ref( extractor ) );

	std::vector<double> const psi = ramp( 64 );

	double const runOne = cache.at( psi, 0.5, vPrimeOf() );
	BOOST_TEST( extractor.calls == 1 );
	BOOST_TEST( cache.held() );

	// The second run: same psi bits, different everything else.
	extractor.generation = 1.0;

	// Without the reset the cache serves, and it is RIGHT to -- it caches on
	// psi and psi has not moved. Asserted because it is what makes the reset
	// necessary rather than decorative.
	BOOST_TEST( cache.at( psi, 0.5, vPrimeOf() ) == runOne );
	BOOST_TEST( extractor.calls == 1 );

	cache.resetForRun();
	BOOST_TEST( !cache.held() );

	double const runTwo = cache.at( psi, 0.5, vPrimeOf() );
	std::printf( "  run 1 %.6e, run 2 after resetForRun %.6e, extractions %d\n",
	             runOne, runTwo, extractor.calls );

	BOOST_TEST( extractor.calls == 2,
	            "resetForRun() did not force a recomputation" );
	BOOST_TEST( runTwo != runOne,
	            "the second run reproduced the first exactly, so the extractor's "
	            "state change is not visible and this test cannot fail" );

	// Idempotent, and safe on an empty cache.
	cache.resetForRun();
	cache.resetForRun();
	BOOST_TEST( !cache.held() );
}

BOOST_AUTO_TEST_CASE( a_failed_extraction_commits_nothing )
{
	// If the extractor throws, the cache must be left holding NOTHING. Leaving
	// the previous family under the previous key serves a different state's
	// geometry the moment the caller retries at the old psi; leaving it under
	// the new key serves it immediately. Both are the fudge section 8 forbids.
	CountingExtractor extractor;
	meq::GeometryCache cache( std::ref( extractor ) );

	std::vector<double> const good = ramp( 32 );
	std::vector<double> const bad = ramp( 32, 2.0 );

	double const before = cache.at( good, 0.5, vPrimeOf() );
	BOOST_TEST( cache.held() );
	BOOST_TEST( extractor.calls == 1 );

	extractor.refuse = true;
	BOOST_CHECK_THROW( cache.at( bad, 0.5, vPrimeOf() ), std::runtime_error );

	BOOST_TEST( !cache.held(),
	            "a failed extraction left a family in the cache" );
	BOOST_TEST( cache.extractions() == 2u,
	            "extractions() counts attempts, so the one that threw is in it" );

	// The old psi must NOT be served from what is left behind.
	extractor.refuse = false;
	double const after = cache.at( good, 0.5, vPrimeOf() );
	std::printf( "  before %.6e, after a refusal %.6e, extractions %zu\n",
	             before, after, cache.extractions() );
	BOOST_TEST( extractor.calls == 3,
	            "the psi that succeeded before the refusal was served from a "
	            "family the refusal should have dropped" );
	BOOST_TEST( after == before );
}

BOOST_AUTO_TEST_CASE( a_state_meq_cannot_evaluate_at_is_refused )
{
	// MANTA-COUPLING.md section 1 guarantees meq is called at states far from
	// equilibrium as a normal part of the consumer's Newton iteration. A psi
	// carrying a NaN must come back as a refusal the integrator can retry from,
	// not as a geometry -- and it must be refused BEFORE the extractor is
	// entered, since a NaN reaching the tracer produces an unbounded search
	// rather than a message.
	CountingExtractor extractor;
	meq::GeometryCache cache( std::ref( extractor ) );

	std::vector<double> psi = ramp( 16 );
	BOOST_CHECK_NO_THROW( cache.family( psi ) );
	BOOST_TEST( extractor.calls == 1 );

	psi[ 7 ] = std::numeric_limits<double>::quiet_NaN();
	BOOST_CHECK_THROW( cache.family( psi ), std::invalid_argument );
	BOOST_TEST( extractor.calls == 1,
	            "the extractor was entered with a NaN in psi" );

	psi[ 7 ] = std::numeric_limits<double>::infinity();
	BOOST_CHECK_THROW( cache.family( psi ), std::invalid_argument );
	BOOST_TEST( extractor.calls == 1 );

	std::vector<double> const nothing;
	BOOST_CHECK_THROW( cache.family( nothing ), std::invalid_argument );

	// A cache that cannot fill itself would report every query as a miss and
	// then fail with something unrelated to the cause.
	BOOST_CHECK_THROW( meq::GeometryCache( meq::GeometryCache::Extractor() ),
	                   std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( the_pointwise_call_pattern_costs_one_extraction_per_state )
{
	// MANTA-COUPLING.md section 5's pattern, simulated: Geometry is called once
	// per physics node, per residual evaluation, handed the whole psi vector
	// each time. A Newton step visits a handful of distinct states and each one
	// is walked over every node several times.
	CountingExtractor extractor;
	meq::GeometryCache cache( std::ref( extractor ) );

	std::size_t const nodes = 60;
	std::size_t const residualsPerState = 3;
	std::size_t const states = 5;

	for ( std::size_t state = 0; state < states; ++state )
	{
		std::vector<double> const psi =
			ramp( 400, 1.0 + 0.1*static_cast<double>( state ) );

		for ( std::size_t residual = 0; residual < residualsPerState; ++residual )
			for ( std::size_t node = 0; node < nodes; ++node )
			{
				double const label = 0.25 + 0.7*static_cast<double>( node )
					/static_cast<double>( nodes );
				cache.at( psi, label, vPrimeOf() );
			}
	}

	std::size_t const expected = states*residualsPerState*nodes;
	std::printf( "  %zu queries over %zu states: %zu extractions, %zu hits\n",
	             cache.queries(), states, cache.extractions(), cache.hits() );

	BOOST_TEST( cache.queries() == expected );
	BOOST_TEST( cache.extractions() == states,
	            "the family was extracted more than once per distinct psi" );
	BOOST_TEST( cache.hits() == expected - states );

	// The key comparison is what a hit costs, and it is O( nDOF ) against an
	// extraction that traces and integrates a whole family. PRINTED AND NOT
	// ASSERTED ON: a timing is a measurement about this machine, which is what
	// tests/performance exists for.
	std::vector<double> const psi = ramp( 400 );
	auto const start = std::chrono::steady_clock::now();
	for ( int i = 0; i < 100000; ++i )
		cache.at( psi, 0.5, vPrimeOf() );
	auto const stop = std::chrono::steady_clock::now();
	std::printf( "  a served query over 400 dofs: %.3f us\n",
	             std::chrono::duration<double, std::micro>( stop - start ).count()
	             /100000.0 );
	std::fflush( stdout );
}
