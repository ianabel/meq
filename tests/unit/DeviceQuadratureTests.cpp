// Unit tests for meq::plasmaCurrentPerElementOnDevice() -- the plasma-current
// quadrature as a device element loop.
//
// WHY THIS ONE IS A TOLERANCE WHERE DeviceRowsTests IS AN EQUALITY, WHICH IS A
// DECISION AND NOT A WEAKENING. M-190 established that a device kernel doing a
// multiply-accumulate differs from the host in its last bits unless the FMA
// contraction is disabled, and that the trade is per kernel:
//
//   * meq::DeviceRows exists to REPRODUCE meq::CompressedRows::dot(), whose
//     bit-identity with a dense loop is a documented contract. There the FMA is
//     given up (-fmad=false) and the test asserts ==.
//   * this loop exists to COMPUTE an integral. The arithmetic is the work, an
//     FMA is the more accurate of the two roundings, and demanding the host's
//     exact rounding would make the kernel worse to make a test simpler.
//
// So the acceptance here is agreement to round-off against the host loop, and
// the case prints the worst relative difference so that a real change shows as
// a jump rather than as a pass.
//
// AND THE REFERENCE IS THE HOST ARITHMETIC ITSELF, not a closed form: the same
// rule, the same shape functions, the same source object. That is deliberate --
// what is under test is the device loop's INDEXING, which is where the
// layouts of GeometricFactors ( NQ x SDIM x NE ), detJ ( NQ x NE ), DofToQuad::B
// ( NQ x ND ) and the E-vector ( ND x NE ) can each be got wrong
// independently, and every one of those errors gives a plausible number.

#define BOOST_TEST_MODULE DeviceQuadratureTests
#include <boost/test/included/unit_test.hpp>

#include "meq/DeviceQuadrature.hpp"
#include "meq/Profiles.hpp"
#include "meq/Source.hpp"
#include "meq/SourceView.hpp"

#include "mfem.hpp"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <vector>

namespace
{
	/// See DeviceRowsTests.cpp: a function-local static, because a
	/// namespace-scope mfem::Device may be constructed before MFEM's own
	/// singleton and then reports itself already configured.
	mfem::Device &deviceOnce()
	{
		char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
		static mfem::Device device( requested ? requested : "debug" );
		return device;
	}

	std::shared_ptr<meq::Profile const> rampProfile()
	{
		std::vector<meq::Knot> knots;
		for ( int i = 0; i <= 8; ++i )
		{
			double const x = 0.125*i;
			knots.push_back( { x, 1.0 + 0.5*x, 0.5 } );
		}
		return std::make_shared<meq::SplineProfile>( std::move( knots ) );
	}

	/// The mesh, shifted off the axis so that F/R is finite everywhere.
	void shiftOffAxis( mfem::Vector const &in, mfem::Vector &out )
	{
		out.SetSize( in.Size() );
		out( 0 ) = in( 0 ) + 1.0;
		out( 1 ) = in( 1 ) - 0.5;
	}
}


BOOST_AUTO_TEST_CASE( theDeviceElementLoopReproducesTheHostQuadrature )
{
	deviceOnce();

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		6, 5, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	mesh.Transform( shiftOffAxis );

	int const order = 2;
	mfem::L2_FECollection collection( order, mesh.Dimension() );
	mfem::FiniteElementSpace space( &mesh, &collection );

	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( rampProfile(), rampProfile(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );

	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( source, view, knots ) );

	// A potential with structure, projected so the dofs are what the element
	// loop will gather rather than something written by hand.
	mfem::GridFunction potential( &space );
	mfem::FunctionCoefficient field( []( mfem::Vector const &x )
	{
		return 0.15 + 0.10*std::sin( 2.0*x( 0 ) )*std::cos( 1.5*x( 1 ) );
	} );
	potential.ProjectCoefficient( field );

	mfem::Array<int> inPlasma( mesh.GetNE() );
	inPlasma = 1;

	int const extraOrder = 4;

	mfem::Vector got;
	BOOST_TEST_REQUIRE(
		meq::plasmaCurrentPerElementOnDevice( space, potential, view, knots,
		                                      inPlasma, {}, {}, extraOrder,
		                                      got ),
		"the device loop refused: either no Device is configured, or one of "
		"its documented refusals fired on a mesh that should not trigger any "
		"of them" );
	BOOST_TEST_REQUIRE( got.Size() == mesh.GetNE() );
	got.HostRead();   // device-valid on return; the reader funnels

	/*
	 * THE HOST REFERENCE, ARRANGED EXACTLY AS assemblePlasmaCurrent() IS.
	 */
	double worst = 0.0;
	double largest = 0.0;
	potential.HostRead();
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		mfem::FiniteElement const &el = *space.GetFE( e );
		mfem::IsoparametricTransformation transformation;
		mesh.GetElementTransformation( e, &transformation );

		mfem::Array<int> dofs;
		space.GetElementDofs( e, dofs );
		int const dof = el.GetDof();
		mfem::Vector shape( dof );
		mfem::Vector point;

		mfem::IntegrationRule const &rule = mfem::IntRules.Get(
			el.GetGeomType(),
			2*el.GetOrder() + transformation.OrderW() + extraOrder );

		double elementSum = 0.0;
		for ( int q = 0; q < rule.GetNPoints(); ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			transformation.SetIntPoint( &ip );
			el.CalcShape( ip, shape );
			transformation.Transform( ip, point );

			double psi = 0.0;
			for ( int j = 0; j < dof; ++j )
				psi += shape( j )*potential( dofs[ j ] );

			elementSum += ip.weight*transformation.Weight()
			              *source.scaledF( point( 0 ), point( 1 ), psi )
			              /point( 0 );
		}

		largest = std::max( largest, std::abs( elementSum ) );
		double const scale = std::max( std::abs( elementSum ), 1.0e-14 );
		worst = std::max( worst, std::abs( got( e ) - elementSum )/scale );
	}

	std::printf( "    %d elements, worst relative difference %.3e\n",
	             mesh.GetNE(), worst );

	// A reference that is identically zero would make the comparison vacuous,
	// which is the failure DeviceRowsTests guards the same way.
	BOOST_TEST( largest > 1.0e-8,
	            "every reference element integral is at round-off, so the "
	            "agreement below means nothing" );

	BOOST_TEST( worst < 1.0e-12,
	            "the device element loop disagrees with the host by " << worst
	            << " relative. At this size that is an INDEXING error rather "
	            "than round-off: check the four layouts this file's opening "
	            "comment lists, each of which is independently wrong-able and "
	            "each of which gives a plausible number." );
}


/*
 * THE REFUSALS, which are the difference between a device path that declines
 * and one that integrates on a rule the caller did not ask for.
 */
BOOST_AUTO_TEST_CASE( theDeviceLoopRefusesWhatItCannotDo )
{
	deviceOnce();

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		3, 3, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	mesh.Transform( shiftOffAxis );

	mfem::L2_FECollection collection( 1, mesh.Dimension() );

	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( rampProfile(), rampProfile(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );
	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( source, view, knots ) );

	// A VECTOR space: the restriction, the shape layout and the gather all
	// assume a scalar one, so this must decline rather than index wrongly.
	{
		mfem::FiniteElementSpace vectorSpace( &mesh, &collection, 2 );
		mfem::GridFunction potential( &vectorSpace );
		potential = 0.2;
		mfem::Array<int> inPlasma( mesh.GetNE() );
		inPlasma = 1;
		mfem::Vector out;
		BOOST_TEST( !meq::plasmaCurrentPerElementOnDevice(
		                vectorSpace, potential, view, knots, inPlasma,
		                {}, {}, 4, out ),
		            "a vector space must be refused" );
	}

	// A per-element flag array of the wrong length: rows built for one mesh
	// used against another.
	{
		mfem::FiniteElementSpace space( &mesh, &collection );
		mfem::GridFunction potential( &space );
		potential = 0.2;
		mfem::Array<int> wrong( mesh.GetNE() + 1 );
		wrong = 1;
		mfem::Vector out;
		BOOST_TEST( !meq::plasmaCurrentPerElementOnDevice(
		                space, potential, view, knots, wrong, {}, {}, 4, out ),
		            "a flag array that is not one per element must be refused" );
	}
}


/*
 * THE CONDUCTOR SHIFT, WHICH IS THE ONE REFUSAL THAT STOOD BETWEEN THIS LOOP
 * AND THE PRODUCTION CONFIGURATION.
 *
 * COIL-SUBTRACTION-PLAN.md CS-2: the source must see the TOTAL flux while the
 * solver holds only the remainder, so `psi_c` is added before F is evaluated
 * and never afterwards. meq::SourceIntegrator caches it per quadrature point
 * as `values[ offsets[ e ] + q ]`, and this case is about that INDEX: a psi_c
 * read at the wrong point is a plausible number and a different equilibrium,
 * which no convergence table can see.
 *
 * The control is that the shift must actually MOVE the answer. Without it the
 * case would pass equally against a kernel that ignored the cache, which is
 * precisely the defect it exists to catch.
 */
BOOST_AUTO_TEST_CASE( theConductorShiftIsAddedAtTheRightQuadraturePoint )
{
	deviceOnce();

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	mesh.Transform( shiftOffAxis );

	int const order = 2;
	mfem::L2_FECollection collection( order, mesh.Dimension() );
	mfem::FiniteElementSpace space( &mesh, &collection );

	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( rampProfile(), rampProfile(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );
	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( source, view, knots ) );

	mfem::GridFunction potential( &space );
	potential = 0.18;

	mfem::Array<int> inPlasma( mesh.GetNE() );
	inPlasma = 1;
	int const extraOrder = 4;

	// The rule this loop will use, so the cache is built on its points.
	mfem::FiniteElement const &el0 = *space.GetFE( 0 );
	mfem::IsoparametricTransformation tr0;
	mesh.GetElementTransformation( 0, &tr0 );
	mfem::IntegrationRule const &rule = mfem::IntRules.Get(
		el0.GetGeomType(), 2*el0.GetOrder() + tr0.OrderW() + extraOrder );
	int const nq = rule.GetNPoints();

	// A shift that VARIES with the point, so an index error cannot cancel.
	std::vector<double> values;
	std::vector<int> offsets;
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		offsets.push_back( e*nq );
		for ( int q = 0; q < nq; ++q )
			values.push_back( 0.01 + 0.002*q + 0.0005*e );
	}

	mfem::Vector without;
	BOOST_TEST_REQUIRE( meq::plasmaCurrentPerElementOnDevice(
		space, potential, view, knots, inPlasma, {}, {}, extraOrder, without ) );

	mfem::Vector with;
	BOOST_TEST_REQUIRE( meq::plasmaCurrentPerElementOnDevice(
		space, potential, view, knots, inPlasma, values, offsets, extraOrder,
		with ) );
	without.HostRead();
	with.HostRead();

	// The host reference, with psi_c added at the same point.
	double worst = 0.0;
	int moved = 0;
	potential.HostRead();
	for ( int e = 0; e < mesh.GetNE(); ++e )
	{
		mfem::FiniteElement const &el = *space.GetFE( e );
		mfem::IsoparametricTransformation transformation;
		mesh.GetElementTransformation( e, &transformation );
		mfem::Array<int> dofs;
		space.GetElementDofs( e, dofs );
		mfem::Vector shape( el.GetDof() ), point;

		double elementSum = 0.0;
		for ( int q = 0; q < nq; ++q )
		{
			mfem::IntegrationPoint const &ip = rule.IntPoint( q );
			transformation.SetIntPoint( &ip );
			el.CalcShape( ip, shape );
			transformation.Transform( ip, point );

			double psi = 0.0;
			for ( int j = 0; j < el.GetDof(); ++j )
				psi += shape( j )*potential( dofs[ j ] );
			psi += values[ static_cast<std::size_t>( offsets[ e ] + q ) ];

			elementSum += ip.weight*transformation.Weight()
			              *source.scaledF( point( 0 ), point( 1 ), psi )
			              /point( 0 );
		}

		double const scale = std::max( std::abs( elementSum ), 1.0e-14 );
		worst = std::max( worst, std::abs( with( e ) - elementSum )/scale );
		if ( with( e ) != without( e ) )
			++moved;
	}

	std::printf( "    shift: worst relative difference %.3e, %d of %d elements moved\n",
	             worst, moved, mesh.GetNE() );

	BOOST_TEST( moved == mesh.GetNE(),
	            "the conductor shift did not move every element, so this case "
	            "would pass against a kernel that ignored the cache" );
	BOOST_TEST( worst < 1.0e-12,
	            "the shifted device loop disagrees with the host by " << worst
	            << " relative -- suspect the cache index, which is "
	            "values[ offsets[ e ] + q ] and not anything else" );
}


/*
 * AND A CACHE BUILT ON A DIFFERENT RULE IS REFUSED RATHER THAN INDEXED ANYWAY.
 */
BOOST_AUTO_TEST_CASE( theDeviceLoopRefusesAConductorCacheFromAnotherRule )
{
	deviceOnce();

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		3, 3, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	mesh.Transform( shiftOffAxis );
	mfem::L2_FECollection collection( 1, mesh.Dimension() );
	mfem::FiniteElementSpace space( &mesh, &collection );

	double const psiAxis = 0.4;
	meq::NormalisedMHDSource source( rampProfile(), rampProfile(), psiAxis );
	source.setNormalisation( psiAxis, 0.05 );
	meq::SourceView view{};
	std::vector<meq::KnotView> knots;
	BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( source, view, knots ) );

	mfem::GridFunction potential( &space );
	potential = 0.2;
	mfem::Array<int> inPlasma( mesh.GetNE() );
	inPlasma = 1;

	// Offsets that are not e*nq: a cache from a rule with a different point
	// count, whose values are not these points.
	std::vector<int> offsets;
	for ( int e = 0; e < mesh.GetNE(); ++e )
		offsets.push_back( e*3 );
	std::vector<double> values( static_cast<std::size_t>( mesh.GetNE() )*3, 0.01 );

	mfem::Vector out;
	BOOST_TEST( !meq::plasmaCurrentPerElementOnDevice(
	                space, potential, view, knots, inPlasma, values, offsets,
	                4, out ),
	            "a conductor cache whose offsets do not match this rule must "
	            "be refused: its values are not these quadrature points" );
}
