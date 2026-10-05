// Unit tests for meq::AxisSearch -- the magnetic axis found without bringing
// the iterate to the host.
//
// WHAT THIS ASSERTS IS THE CANDIDATE LIST AND NOT THE ANSWER, WHICH IS A
// STRONGER THING AND IS THE ONLY ONE THAT CAN FAIL USEFULLY.
// meq::axisCandidatesOnDevice() is a port of
// meq::CriticalPointFinder::sweep()'s root finding; what the solver then does
// with the list is host code shared by both routes. So the place a port can go
// wrong is the list, and comparing the final psi_ax -- which
// tests/convergence/DeviceParity.cpp does end to end -- would pass a port that
// found a different set of critical points and happened to pick the same one.
//
// AND IT CONFIGURES AN mfem::Device, which is what makes any of it meaningful:
// the routine REFUSES without one, so a binary with no Device would assert on
// a refusal and call it a pass. "debug" by default, since it has device memory
// semantics with host arithmetic and needs no card; MEQ_TEST_DEVICE=cuda puts
// the same cases through nvcc's build of the same header.
//
// THE TOLERANCE IS A POSITION AND IT IS MEASURED, NOT CHOSEN. The two routes
// evaluate one polynomial through two bases -- MFEM's CalcShape and
// meq::DeviceBasis' change of basis -- so they agree to round-off and not to
// the bit, and what they are rooting is a STATIONARY point, where M-220
// measured that a field perturbation of d moves the located position by about
// sqrt( d ) and its value by d. The cases print both so the assertion can be
// read against what it is asserting on.

#define BOOST_TEST_MODULE AxisSearchTests
#include <boost/test/included/unit_test.hpp>

#include "meq/AxisSearch.hpp"
#include "meq/CriticalPoints.hpp"
#include "meq/Crossings.hpp"
#include "meq/DeviceBasis.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace
{
	/// One Device for the binary, function-local for the
	/// static-initialisation-order reason DeviceRowsTests.cpp sets out.
	mfem::Device &deviceOnce()
	{
		char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
		static mfem::Device device( requested ? requested : "debug" );
		return device;
	}

	/**
	 * A field with critical points in known places.
	 *
	 * `psi = cos( a R ) cos( b z )` has a lattice of maxima, minima and
	 * saddles, and `q = -grad psi` vanishes at every one of them -- so the
	 * flux interpolant has roots the finder must locate and classify, several
	 * of each kind, rather than the single O-point a plasma-shaped fixture
	 * would give. A degree-2 interpolant of it is not that function, which is
	 * the point: both routes root the SAME interpolant, and what is being
	 * compared is the two evaluations of it.
	 */
	struct Fixture
	{
		mfem::Mesh mesh;
		mfem::L2_FECollection collection;
		mfem::FiniteElementSpace fluxSpace;
		mfem::FiniteElementSpace potentialSpace;
		mfem::Vector state;
		mfem::GridFunction flux;
		mfem::GridFunction potential;
		mfem::Array<int> offsets;

		explicit Fixture( int degree, int cells )
			: mesh( mfem::Mesh::MakeCartesian2D(
				  cells, cells, mfem::Element::TRIANGLE, false, 2.0, 2.0 ) ),
			  collection( degree, 2, mfem::BasisType::GaussLobatto ),
			  fluxSpace( &mesh, &collection, 2 ),
			  potentialSpace( &mesh, &collection ),
			  flux( &fluxSpace ),
			  potential( &potentialSpace )
		{
			offsets.SetSize( 3 );
			offsets[ 0 ] = 0;
			offsets[ 1 ] = fluxSpace.GetVSize();
			offsets[ 2 ] = offsets[ 1 ] + potentialSpace.GetVSize();
			state.SetSize( offsets[ 2 ] );

			double const a = 2.3;
			double const b = 1.9;

			mfem::FunctionCoefficient psi(
				[ = ]( mfem::Vector const &x )
				{
					return std::cos( a*x( 0 ) )*std::cos( b*x( 1 ) );
				} );
			mfem::VectorFunctionCoefficient minusQ(
				2, [ = ]( mfem::Vector const &x, mfem::Vector &out )
				{
					// The state's flux block holds -q, as mfem::DarcyForm
					// does; q = -grad psi, so the block holds +grad psi.
					out( 0 ) = -a*std::sin( a*x( 0 ) )*std::cos( b*x( 1 ) );
					out( 1 ) = -b*std::cos( a*x( 0 ) )*std::sin( b*x( 1 ) );
				} );

			mfem::GridFunction blockFlux( &fluxSpace );
			blockFlux.ProjectCoefficient( minusQ );
			potential.ProjectCoefficient( psi );

			for ( int i = 0; i < blockFlux.Size(); ++i )
				state( i ) = blockFlux( i );
			for ( int i = 0; i < potential.Size(); ++i )
				state( offsets[ 1 ] + i ) = potential( i );

			// What the host finder is given: the NEGATED flux block, which is
			// the one thing locateAxisPoint() does to the iterate before
			// handing it over. CriticalPoints.hpp records that a finder handed
			// the raw block puts every Maximum where a Minimum is.
			for ( int i = 0; i < flux.Size(); ++i )
				flux( i ) = -state( i );
		}
	};

	/// The device route's inputs, bound exactly as the solver binds them.
	struct Route
	{
		/// Every element, ascending -- which is what a cold search passes and
		/// is the selection these cases compare against sweep() with.
		mfem::Array<int> everyElement;
		meq::DeviceTriangleBasis fluxBasis;
		meq::DeviceTriangleBasis potentialBasis;
		mfem::Array<meq::AffineTriangle> maps;
		mfem::Vector target;
		mfem::Array<meq::AxisCandidate> raw;

		bool run( Fixture &fixture, meq::CriticalPointFinder const &defaults,
		          std::vector<meq::CriticalPoint> &out, int &count )
		{
			out.clear();
			if ( !fluxBasis.build( fixture.fluxSpace )
			     || !potentialBasis.build( fixture.potentialSpace )
			     || !meq::affineTrianglesOf( fixture.mesh, maps ) )
				return false;

			if ( !meq::axisFluxScaleOnDevice(
				     fixture.state, fixture.fluxSpace.GetNDofs(),
				     static_cast<mfem::real_t>( defaults.newtonTolerance() ),
				     target ) )
				return false;

			everyElement.SetSize( fixture.mesh.GetNE() );
			everyElement.UseDevice( true );
			{
				int *const host = everyElement.HostWrite();
				for ( int e = 0; e < fixture.mesh.GetNE(); ++e )
					host[ e ] = e;
			}

			mfem::IntegrationPoint const centre =
				mfem::Geometries.GetCenter( mfem::Geometry::TRIANGLE );

			meq::AxisSearchInputs inputs{};
			inputs.flux = fluxBasis.view();
			inputs.potential = potentialBasis.view();
			inputs.maps = maps.Read();
			inputs.state = meq::countedDeviceRead( fixture.state, __LINE__ );
			inputs.target = target.Read();
			inputs.centreXi = static_cast<mfem::real_t>( centre.x );
			inputs.centreEta = static_cast<mfem::real_t>( centre.y );
			inputs.selected = meq::countedDeviceRead( everyElement, __LINE__ );
			inputs.selectedCount = everyElement.Size();
			inputs.elements = fixture.mesh.GetNE();
			inputs.fluxScalarDofs = fixture.fluxSpace.GetNDofs();
			inputs.potentialOffset = fixture.offsets[ 1 ];
			inputs.jacobianStep = static_cast<mfem::real_t>(
				defaults.newtonJacobianStep() );
			inputs.containment = static_cast<mfem::real_t>(
				defaults.rootContainment() );
			inputs.maxIterations = defaults.newtonMaxIterations();

			if ( !meq::axisCandidatesOnDevice( inputs, raw, count ) )
				return false;

			rawCount = count;
			for ( int i = 0; i < count; ++i )
			{
				meq::CriticalPoint point;
				point.radius = raw[ i ].radius;
				point.z = raw[ i ].z;
				point.psi = raw[ i ].psi;
				point.type =
					static_cast<meq::CriticalPointType>( raw[ i ].type );
				point.element = raw[ i ].element;
				point.overshoot = raw[ i ].overshoot;
				out.push_back( point );
			}

			// THE HOST'S OWN MERGE, WHICH THE SOLVER ALSO CALLS. The device
			// returns every root it found, in element order; sweep() merges
			// two candidates of the same type within half an element, keeping
			// the one least outside its own. Comparing an unmerged list
			// against a merged one is comparing two different things, and
			// that is what the first version of this file did -- it read
			// 7.2e-03 and called it a failure of the port.
			out = defaults.mergeCandidates( out );
			count = static_cast<int>( out.size() );
			return true;
		}

		/// How many roots the kernel found before the merge.
		int rawCount = 0;
	};

	/// The nearest host candidate of the same type, and how far away it is.
	double nearestOfType( std::vector<meq::CriticalPoint> const &all,
	                      meq::CriticalPoint const &wanted )
	{
		double best = 1.0e300;
		for ( auto const &point : all )
		{
			if ( point.type != wanted.type )
				continue;
			double const dr = point.radius - wanted.radius;
			double const dz = point.z - wanted.z;
			best = std::min( best, std::sqrt( dr*dr + dz*dz ) );
		}
		return best;
	}
}

// THE DEVICE FINDS THE SAME CRITICAL POINTS AS THE HOST SWEEP.
//
// Both directions are asserted, which is the half a one-sided check misses: a
// kernel that found a SUBSET would pass "every device point is a host point",
// and a kernel that accepted spurious roots would pass the converse. A missing
// candidate is how a coil's O-point wins a competition it should have lost.
BOOST_AUTO_TEST_CASE( theDeviceRootsWhatTheHostSweepRoots )
{
	deviceOnce();

	for ( int degree = 1; degree <= 3; ++degree )
	{
		Fixture fixture( degree, 6 );
		meq::CriticalPointFinder const finder( fixture.flux,
		                                       fixture.potential );

		std::vector<meq::CriticalPoint> const host = finder.sweep();

		Route route;
		std::vector<meq::CriticalPoint> device;
		int count = 0;
		BOOST_TEST_REQUIRE( route.run( fixture, finder, device, count ) );
		BOOST_TEST_REQUIRE( count == static_cast<int>( device.size() ) );

		// BOTH LISTS ARE MERGED BY THE SAME ROUTINE, so this is a comparison
		// of two candidate sets and not of a set against a superset. The
		// counts are asserted equal below for the same reason: once the merge
		// is shared, a difference in COUNT is a real finding rather than an
		// artefact of one side not having been tidied.
		double worstDeviceToHost = 0.0;
		for ( auto const &point : device )
			worstDeviceToHost = std::max( worstDeviceToHost,
			                              nearestOfType( host, point ) );

		double worstHostToDevice = 0.0;
		for ( auto const &point : host )
			worstHostToDevice = std::max( worstHostToDevice,
			                              nearestOfType( device, point ) );

		std::printf( "  degree %d: host %zu candidates, device %d "
		             "(%d before the merge); device->host %.3e, "
		             "host->device %.3e\n",
		             degree, host.size(), count, route.rawCount,
		             worstDeviceToHost, worstHostToDevice );

		BOOST_TEST_REQUIRE( !host.empty() );
		BOOST_TEST_REQUIRE( count > 0 );
		BOOST_TEST( count == static_cast<int>( host.size() ) );
		// sqrt( eps ) is the floor a STATIONARY point's position can be
		// located to at all (M-220), so this is one decade above the best
		// either route could do, not a convenience.
		BOOST_TEST( worstDeviceToHost < 1.0e-6 );
		BOOST_TEST( worstHostToDevice < 1.0e-6 );
	}
}

// AND IT CLASSIFIES THEM THE SAME WAY.
//
// The type decides which candidates even enter the competition -- `wanted` is
// Maximum or Minimum and a Saddle is dropped -- so a port that located every
// point correctly and classified them by the wrong sign would find no axis at
// all on half the iterates and the right one on the other half. The
// determinant-then-trace ladder is three branches and it is worth pinning that
// all three agree.
BOOST_AUTO_TEST_CASE( theDeviceClassifiesThemAsTheHostDoes )
{
	deviceOnce();

	Fixture fixture( 2, 6 );
	meq::CriticalPointFinder const finder( fixture.flux, fixture.potential );
	std::vector<meq::CriticalPoint> const host = finder.sweep();

	Route route;
	std::vector<meq::CriticalPoint> device;
	int count = 0;
	BOOST_TEST_REQUIRE( route.run( fixture, finder, device, count ) );

	int maxima = 0, minima = 0, saddles = 0;
	for ( auto const &point : host )
	{
		if ( point.type == meq::CriticalPointType::Maximum ) { ++maxima; }
		else if ( point.type == meq::CriticalPointType::Minimum ) { ++minima; }
		else if ( point.type == meq::CriticalPointType::Saddle ) { ++saddles; }
	}
	std::printf( "  host: %d maxima, %d minima, %d saddles\n",
	             maxima, minima, saddles );

	// The fixture must actually present all three, or this case asserts
	// nothing about two of the branches.
	BOOST_TEST_REQUIRE( maxima > 0 );
	BOOST_TEST_REQUIRE( minima > 0 );
	BOOST_TEST_REQUIRE( saddles > 0 );

	// Every host point has a device point of the SAME type within the floor.
	for ( auto const &point : host )
		BOOST_TEST( nearestOfType( device, point ) < 1.0e-6 );

	// And psi at a located point agrees to round-off, which is the quantity
	// the border row actually carries -- the position costs half the digits
	// and the value costs none.
	double worstPsi = 0.0;
	for ( auto const &point : host )
	{
		double best = 1.0e300;
		double value = 0.0;
		for ( auto const &other : device )
		{
			if ( other.type != point.type )
				continue;
			double const dr = other.radius - point.radius;
			double const dz = other.z - point.z;
			double const distance = std::sqrt( dr*dr + dz*dz );
			if ( distance < best ) { best = distance; value = other.psi; }
		}
		if ( best < 1.0e-6 )
			worstPsi = std::max( worstPsi, std::fabs( value - point.psi ) );
	}
	double const eps = static_cast<double>(
		std::numeric_limits<mfem::real_t>::epsilon() );
	std::printf( "  worst psi disagreement at a matched point: %.3e "
	             "(%.1f eps)\n", worstPsi, worstPsi/eps );
	// In epsilons, because the claim is about the two evaluations of one
	// polynomial and not about `double`. M-219's repair: 1.0e-12 is true in
	// double and an artefact of it, and a basis three decades worse would
	// have passed.
	BOOST_TEST( worstPsi < 100.0*eps );
}

// AND THE STOPPING SCALE IS THE SAME NUMBER, COMPUTED WHERE THE DATA IS.
//
// meq::axisFluxScaleOnDevice() reproduces CriticalPointFinder::fluxScale()
// times the tolerance, and it exists so that the scale never comes to the
// host: meq::AxisSearchInputs::target points at the device buffer. A scale
// that disagreed would give the Newton a different stopping rule and a
// different set of accepted roots, which is a difference no comparison of
// positions would separate from a bad Jacobian.
BOOST_AUTO_TEST_CASE( theStoppingScaleAgreesAndStaysOnTheDevice )
{
	deviceOnce();

	Fixture fixture( 2, 5 );
	meq::CriticalPointFinder const finder( fixture.flux, fixture.potential );

	mfem::Vector target;
	BOOST_TEST_REQUIRE( meq::axisFluxScaleOnDevice(
		fixture.state, fixture.fluxSpace.GetNDofs(),
		static_cast<mfem::real_t>( finder.newtonTolerance() ), target ) );

	// The host's own expression, written out here rather than called:
	// fluxScale() is private, and reproducing it is what is being checked.
	int const scalarDofs = fixture.fluxSpace.GetNDofs();
	double worst = 0.0;
	for ( int i = 0; i < scalarDofs; ++i )
	{
		double const a = fixture.flux( i );
		double const b = fixture.flux( scalarDofs + i );
		worst = std::max( worst, std::sqrt( a*a + b*b ) );
	}
	double const expected =
		finder.newtonTolerance()*( worst > 0.0 ? worst : 1.0 );

	// AND THE SAME SUM IN real_t THROUGHOUT, WHICH IS WHAT THE KERNEL DOES.
	//
	// THE FUNNEL IS NOT DECORATION: axisFluxScaleOnDevice() reads the state
	// through countedDeviceRead(), which leaves it device-valid -- and under
	// mfem::Device( "debug" ) that protects the host page, so the raw loop
	// below faults. It did, on the first run of this case. The same class as
	// everything else in this file's subject, met in the test rather than in
	// the code under test.
	meq::countedHostRead( fixture.state, __LINE__ );
	mfem::real_t narrow = 0.0;
	for ( int i = 0; i < scalarDofs; ++i )
	{
		mfem::real_t const a = fixture.state( i );
		mfem::real_t const b = fixture.state( scalarDofs + i );
		mfem::real_t const size = std::sqrt( a*a + b*b );
		narrow = size > narrow ? size : narrow;
	}
	mfem::real_t const narrowExpected =
		static_cast<mfem::real_t>( finder.newtonTolerance() )
		*( narrow > mfem::real_t( 0.0 ) ? narrow : mfem::real_t( 1.0 ) );

	meq::countedHostRead( target, __LINE__ );
	double const eps = static_cast<double>(
		std::numeric_limits<mfem::real_t>::epsilon() );
	std::printf( "  stopping scale: device %.17g, host(real_t) %.17g, "
	             "host(fluxScale) %.17g  [%.1f eps apart]\n",
	             static_cast<double>( target( 0 ) ),
	             static_cast<double>( narrowExpected ), expected,
	             std::fabs( target( 0 ) - expected )/( expected*eps ) );

	/*
	 * TWO CLAIMS, AND THE FIRST VERSION OF THIS CASE CONFLATED THEM.
	 *
	 * It asserted `target( 0 ) == expected` with the reason *"a max is
	 * exactly associative, so this one IS exact"*. The reduction is exactly
	 * associative and that half is true; what is not is the TERM. The host's
	 * `fluxScale()` reads each pair into `double` locals, takes
	 * `sqrt( a*a + b*b )` in double and rounds the result to `real_t`, while
	 * the kernel does the whole term in `real_t`. In a double build those are
	 * the same operation, so `==` held and read as a general claim. Under
	 * `MFEM_PRECISION = single` they are one rounding apart per term.
	 *
	 * **So the exactness belongs to the REDUCTION and not to the function**,
	 * and the two are asserted separately: bit-equality against a `real_t`
	 * reduction, which is what the kernel is, and a few epsilons against
	 * `fluxScale()`'s own arithmetic, which is what it must agree with.
	 * *A bit-exactness claim that is only true in one build is a claim about
	 * that build.*
	 */
	BOOST_TEST( target( 0 ) == narrowExpected );
	BOOST_TEST( std::fabs( target( 0 ) - expected ) <= 4.0*eps*expected );
}

// THE MERGE IS THE QUADRATIC LOOP'S, TO THE BIT, ON A LIST LONG ENOUGH TO MATTER.
//
// mergeCandidates() walks a spatial grid rather than the whole kept list, and
// the claim that makes that safe is that it returns EXACTLY what the walk did:
// same points, same order, same replacements. So the oracle here is the walk
// itself, verbatim, and the comparison is `==` on every field. The list is
// built to reach every branch of the rule -- coincident points (the
// separation test), same-type neighbours within half an element (the merge),
// a later candidate less outside its element (the replacement, which MOVES a
// kept point), and different types side by side (never merged) -- and it is
// two candidates per element, which is the shape a flux-free guess hands the
// sweep and the shape that made the walk 22 seconds on fixed-f-diiid.
BOOST_AUTO_TEST_CASE( theMergeIsTheWalkItReplaces )
{
	Fixture fixture( 1, 40 );
	meq::CriticalPointFinder const finder( fixture.flux, fixture.potential );
	mfem::Mesh &mesh = fixture.mesh;

	// A deterministic scatter: an LCG, so the list is the same on every run.
	unsigned long long state = 12345;
	auto next = [ & ]()
	{
		state = state*6364136223846793005ULL + 1442695040888963407ULL;
		return static_cast<double>( state >> 11 )/9007199254740992.0;
	};

	std::vector<meq::CriticalPoint> ordered;
	mfem::IsoparametricTransformation scratch;
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		mesh.GetElementTransformation( e, &scratch );
		double const h = mesh.GetElementSize( &scratch );
		mfem::IntegrationPoint centre;
		centre.Set2( 1.0/3.0, 1.0/3.0 );
		mfem::Vector x( 2 );
		scratch.Transform( centre, x );
		for ( int seed = 0; seed < 2; ++seed )
		{
			meq::CriticalPoint point;
			point.element = e;
			point.radius = x( 0 ) + ( next() - 0.5 )*h;
			point.z = x( 1 ) + ( next() - 0.5 )*h;
			// One in eight sits exactly on the previous candidate, which is
			// the separation branch and no other.
			if ( !ordered.empty() && next() < 0.125 )
			{
				point.radius = ordered.back().radius;
				point.z = ordered.back().z;
			}
			point.psi = next();
			point.overshoot = next() - 0.5;
			int const pick = static_cast<int>( 3.0*next() );
			point.type = pick == 0 ? meq::CriticalPointType::Maximum
			             : pick == 1 ? meq::CriticalPointType::Minimum
			                         : meq::CriticalPointType::Saddle;
			ordered.push_back( point );
		}
	}

	// THE ORACLE, which is the walk mergeCandidates() replaced, verbatim.
	double diameter = 0.0;
	{
		mfem::real_t lo[ 2 ] = { mesh.GetVertex( 0 )[ 0 ], mesh.GetVertex( 0 )[ 1 ] };
		mfem::real_t hi[ 2 ] = { lo[ 0 ], lo[ 1 ] };
		for ( int v = 1; v < mesh.GetNV(); ++v )
			for ( int d = 0; d < 2; ++d )
			{
				lo[ d ] = std::min( lo[ d ], mesh.GetVertex( v )[ d ] );
				hi[ d ] = std::max( hi[ d ], mesh.GetVertex( v )[ d ] );
			}
		diameter = std::sqrt( ( hi[ 0 ] - lo[ 0 ] )*( hi[ 0 ] - lo[ 0 ] )
		                      + ( hi[ 1 ] - lo[ 1 ] )*( hi[ 1 ] - lo[ 1 ] ) );
	}
	double const separation = finder.rootSeparation();
	std::vector<meq::CriticalPoint> expected;
	int replaced = 0;
	int separated = 0;
	for ( meq::CriticalPoint const &point : ordered )
	{
		mesh.GetElementTransformation( point.element, &scratch );
		double const reach = 0.5*mesh.GetElementSize( &scratch );
		bool duplicate = false;
		for ( std::size_t j = 0; j < expected.size(); ++j )
		{
			double const dr = expected[ j ].radius - point.radius;
			double const dz = expected[ j ].z - point.z;
			double const distance = std::sqrt( dr*dr + dz*dz );
			if ( distance < separation*diameter )
			{
				duplicate = true;
				++separated;
				break;
			}
			if ( expected[ j ].type == point.type && distance < reach )
			{
				duplicate = true;
				if ( point.overshoot < expected[ j ].overshoot )
				{
					expected[ j ] = point;
					++replaced;
				}
				break;
			}
		}
		if ( !duplicate )
			expected.push_back( point );
	}

	std::vector<meq::CriticalPoint> const got = finder.mergeCandidates( ordered );

	std::printf( "  %zu candidates over %d elements: %zu kept, %d separated, "
	             "%d replaced\n", ordered.size(), mesh.GetNE(), expected.size(),
	             separated, replaced );

	// Every branch was reached, or the comparison below proves less than it says.
	BOOST_TEST_REQUIRE( separated > 0 );
	BOOST_TEST_REQUIRE( replaced > 0 );
	BOOST_TEST_REQUIRE( expected.size() < ordered.size() );

	BOOST_TEST_REQUIRE( got.size() == expected.size() );
	int differing = 0;
	for ( std::size_t i = 0; i < got.size(); ++i )
		if ( got[ i ].radius != expected[ i ].radius
		     || got[ i ].z != expected[ i ].z
		     || got[ i ].psi != expected[ i ].psi
		     || got[ i ].type != expected[ i ].type
		     || got[ i ].element != expected[ i ].element
		     || got[ i ].overshoot != expected[ i ].overshoot )
			++differing;
	BOOST_TEST( differing == 0 );
}

// AND IT RETURNS THE WHOLE LIST, HOWEVER LONG, AND IT IS THE HOST'S.
//
// The route has no ceiling: the list is at most two candidates per element,
// and the long case is real -- a flux-free guess makes every element report a
// root. A truncated list would be a competition missing an entrant, which is
// how a conductor's O-point would win one, so the claim is the whole list or
// nothing, and the whole list here is the host sweep's. The fixture is that
// guess: the cosine potential with the flux block ZEROED, as `Type = "ramp"`
// leaves it, on 5000 elements -- past the 4096 a ceiling once refused.
BOOST_AUTO_TEST_CASE( itReturnsTheWholeCandidateListHoweverLong )
{
	deviceOnce();

	Fixture fixture( 1, 50 );
	{
		mfem::real_t *const host = fixture.state.HostReadWrite();
		for ( int i = 0; i < fixture.offsets[ 1 ]; ++i )
			host[ i ] = 0.0;
		fixture.flux = 0.0;
		fixture.flux.HostRead();
	}
	meq::CriticalPointFinder const finder( fixture.flux, fixture.potential );
	std::vector<meq::CriticalPoint> const host = finder.sweep();

	Route route;
	std::vector<meq::CriticalPoint> device;
	int count = 0;
	bool const took = route.run( fixture, finder, device, count );

	std::printf( "  %d elements: %d roots before the merge, %d after, "
	             "host %zu, route %s\n",
	             fixture.mesh.GetNE(), route.rawCount, count, host.size(),
	             took ? "taken" : "REFUSED" );

	BOOST_TEST_REQUIRE( took );
	BOOST_TEST( route.rawCount > 4096 );
	BOOST_TEST( count == static_cast<int>( host.size() ) );
	// Merged by one routine from two lists of the same roots, so the merged
	// lists are the same points; how close is theDeviceRootsWhatTheHostSweepRoots's
	// question, asked here only of the element each came from.
	int differing = 0;
	for ( std::size_t i = 0; i < host.size() && i < device.size(); ++i )
		if ( host[ i ].element != device[ i ].element
		     || host[ i ].type != device[ i ].type )
			++differing;
	BOOST_TEST( differing == 0 );
}
