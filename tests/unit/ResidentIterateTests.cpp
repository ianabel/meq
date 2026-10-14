// Unit tests for meq::DevicePlasmaLoops and the three small kernels beside it
// -- the device routes of the last host readers of the bordered iterate.
//
// AN EQUALITY, NOT A TOLERANCE, AND THAT IS THE OTHER SIDE OF THE TRADE
// DeviceQuadratureTests DESCRIBES. That loop exists to compute an integral and
// is accepted to round-off. These exist to REPLACE host loops the bordered
// Newton already runs, so that a resident iterate is read where it lives
// without changing a digit the host prints -- and they are written for it:
// the radii and weights are the host's own Transform() and Weight() values,
// uploaded once; every expression is in the host's order; and the device file
// is compiled without FMA contraction. Measured on machine-f-diiid with every
// route forced, 126 evaluations of eight routes agreed with their host loops
// to the bit, after one Hermite derivative was rewritten in the host's form.
// This pins that.
//
// THE REFERENCE IS EACH HOST LOOP AS GradShafranov.cpp WRITES IT, inline here
// because the solver's loops are private. A drift between this copy and the
// solver's is the failure mode, so each is a transcription of one function and
// says which.
//
// MEQ_TEST_DEVICE=cuda puts the routes on a card; the default is the debug
// backend, which has device memory semantics with host arithmetic and so
// exercises the indexing and the validity flags without a GPU.

#define BOOST_TEST_MODULE ResidentIterateTests
#include <boost/test/included/unit_test.hpp>

#include "meq/DeviceNPC.hpp"
#include "meq/Profiles.hpp"
#include "meq/ResidentIterate.hpp"
#include "meq/Source.hpp"
#include "meq/SourceView.hpp"

#include "mfem.hpp"

#include "convergence/RoundOffFloor.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
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

	/// A spline whose second derivative is not zero, so the derivative routes
	/// read every Hermite basis function including the one that was not the
	/// host's double.
	std::shared_ptr<meq::Profile const> curvedProfile()
	{
		std::vector<meq::Knot> knots;
		for ( int i = 0; i <= 8; ++i )
		{
			double const x = 0.125*i;
			knots.push_back( { x, 1.0 + 0.5*x - 0.7*x*x, 0.5 - 1.4*x } );
		}
		return std::make_shared<meq::SplineProfile>( std::move( knots ) );
	}

	void shiftOffAxis( mfem::Vector const &in, mfem::Vector &out )
	{
		out.SetSize( in.Size() );
		out( 0 ) = in( 0 ) + 1.0;
		out( 1 ) = in( 1 ) - 0.5;
	}

	/// Everything the routes and the host references share.
	struct Fixture
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			6, 5, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
		int order = 2;
		int extraOrder = 4;
		std::unique_ptr<mfem::L2_FECollection> collection;
		std::unique_ptr<mfem::FiniteElementSpace> space;

		/// The plasma source behind a ScaledProfile, so the profile-scale
		/// path -- whose association against the division moved -- is read.
		meq::NormalisedMHDSource source{
			curvedProfile(),
			std::make_shared<meq::ScaledProfile>( curvedProfile(), 2.5 ),
			0.4 };

		/// The bordered unknown's layout: something before the potential block
		/// and something after it, so an offset error reads the wrong entries.
		int potentialStart = 37;
		mfem::Vector state;

		std::vector<int> flags;
		std::vector<double> shiftValues;
		std::vector<int> shiftOffsets;
		int quadPoints = 0;

		meq::DeviceConstants constants;
		meq::SourceView view{};

		Fixture()
		{
			deviceOnce();
			mesh.Transform( shiftOffAxis );
			collection = std::make_unique<mfem::L2_FECollection>(
				order, mesh.Dimension() );
			space = std::make_unique<mfem::FiniteElementSpace>(
				&mesh, collection.get() );

			source.setNormalisation( 0.4, 0.05 );
			source.setCurrentScale( 0.9 );

			int const potentialSize = space->GetVSize();
			state.SetSize( potentialStart + potentialSize + 11 );
			state.UseDevice( true );
			{
				mfem::real_t *const host = state.HostWrite();
				for ( int i = 0; i < state.Size(); ++i )
					host[ i ] = 0.12 + 0.2*std::sin( 0.37*i )*std::cos( 0.11*i );
			}

			// Every fifth element off, so the per-element skip is read.
			int const elements = mesh.GetNE();
			flags.assign( static_cast<std::size_t>( elements ), 1 );
			for ( int e = 0; e < elements; e += 5 )
				flags[ static_cast<std::size_t>( e ) ] = 0;

			// A conductor shift on the source's own points, as
			// meq::SourceIntegrator's cache is laid out.
			quadPoints = mfem::IntRules.Get( mfem::Geometry::TRIANGLE,
			                                 2*order + extraOrder ).GetNPoints();
			for ( int e = 0; e < elements; ++e )
			{
				shiftOffsets.push_back( e*quadPoints );
				for ( int q = 0; q < quadPoints; ++q )
					shiftValues.push_back( 0.01*std::cos( 0.3*e + 0.7*q ) );
			}

			meq::DeviceIntegrandInputs inputs;
			inputs.source = &source;
			inputs.elementInPlasma = flags.data();
			inputs.elementCount = elements;
			inputs.conductorValues = shiftValues.data();
			inputs.conductorOffsets = shiftOffsets.data();
			inputs.conductorValueCount = static_cast<int>( shiftValues.size() );

			std::vector<meq::KnotView> knots;
			BOOST_TEST_REQUIRE( meq::makePlasmaSourceView( source, view, knots ) );
			constants.refresh( inputs, knots, {} );
		}

		/// The host's per-point preamble, as all five loops in GradShafranov.cpp
		/// write it: SetIntPoint, CalcShape, Transform, then psi from zero in
		/// local-dof order and the conductor shift added.
		template <typename Body>
		void hostLoop( Body &&body ) const
		{
			state.HostRead();
			for ( int e = 0; e < mesh.GetNE(); ++e )
			{
				if ( flags[ static_cast<std::size_t>( e ) ] == 0 )
					continue;
				mfem::FiniteElement const &el = *space->GetFE( e );
				mfem::IsoparametricTransformation tr;
				const_cast<mfem::Mesh &>( mesh ).GetElementTransformation( e, &tr );
				mfem::Array<int> dofs;
				space->GetElementDofs( e, dofs );
				int const dof = el.GetDof();
				mfem::Vector shape( dof );
				mfem::Vector point;
				mfem::IntegrationRule const &ir = mfem::IntRules.Get(
					el.GetGeomType(), 2*el.GetOrder() + tr.OrderW() + extraOrder );
				for ( int i = 0; i < ir.GetNPoints(); ++i )
				{
					mfem::IntegrationPoint const &ip = ir.IntPoint( i );
					tr.SetIntPoint( &ip );
					el.CalcShape( ip, shape );
					tr.Transform( ip, point );
					double psi = 0.0;
					for ( int j = 0; j < dof; ++j )
						psi += shape( j )*state( potentialStart + dofs[ j ] );
					psi += shiftValues[ static_cast<std::size_t>(
						shiftOffsets[ static_cast<std::size_t>( e ) ] + i ) ];
					body( e, ip, tr, point, psi, dofs, shape );
				}
			}
		}
	};

	/**
	 * The count of entries that differ, and the largest difference: BITWISE IN
	 * DOUBLE, AND TO THE FORMAT IN SINGLE.
	 *
	 * The host loops evaluate the source in MEQ's MFEM-free double half and
	 * accumulate in double; the device routes are real_t throughout. In single
	 * those are different arithmetic, which no edit to a kernel could make
	 * agree, and they part at about one eps of the vector (measured 1.5e-08 to
	 * 1.6e-07). So an entry counts as differing when it is further than
	 * exactAgreementFloor() of the vector's largest entry, which is EXACTLY ZERO
	 * in double, where this stays the bit equality it always was. What single
	 * gives up is the reassociation check; what it keeps is the layout -- a
	 * wrong element, a dropped point, a transposed scatter is O( 1 ). DeviceRowsTests
	 * makes the same trade. M-273.
	 */
	int differing( mfem::Vector const &got, mfem::Vector const &want,
	               double &largest )
	{
		got.HostRead();
		want.HostRead();
		double scale = 0.0;
		for ( int i = 0; i < want.Size(); ++i )
			scale = std::max( scale, double( std::fabs( want( i ) ) ) );
		double const floor = meq::tests::exactAgreementFloor( scale );
		largest = 0.0;
		int count = 0;
		for ( int i = 0; i < want.Size(); ++i )
		{
			double const gap = std::fabs( double( got( i ) ) - double( want( i ) ) );
			if ( floor == 0.0 ? got( i ) != want( i ) : gap > floor )
			{
				++count;
				largest = std::max( largest, gap );
			}
		}
		return count;
	}
}


BOOST_AUTO_TEST_CASE( theTwoIntegralsMatchTheirHostLoopsToTheBit )
{
	Fixture f;
	meq::DevicePlasmaLoops loops;
	BOOST_TEST_REQUIRE( loops.prepare( *f.space, f.extraOrder ),
	                    "phase B refused an affine L2 triangle mesh, which is "
	                    "the one configuration it exists for" );

	int const elements = f.mesh.GetNE();

	// assemblePlasmaCurrent().
	mfem::Vector current( elements );
	current = 0.0;
	current.HostReadWrite();
	f.hostLoop( [ & ]( int e, mfem::IntegrationPoint const &ip,
	                   mfem::ElementTransformation &tr, mfem::Vector const &point,
	                   double psi, mfem::Array<int> const &,
	                   mfem::Vector const & )
	{
		current( e ) += ip.weight*tr.Weight()
		                *f.source.scaledF( point( 0 ), point( 1 ), psi )
		                /point( 0 );
	} );

	mfem::Vector got;
	BOOST_TEST_REQUIRE( loops.integrate( meq::PlasmaQuantity::Current, f.view,
	                                     f.constants, f.state, f.potentialStart,
	                                     got, nullptr ) );
	double largest = 0.0;
	int const currentOff = differing( got, current, largest );
	BOOST_TEST( currentOff == 0,
	            currentOff << " of " << elements << " per-element currents differ "
	            "from assemblePlasmaCurrent()'s loop, by up to " << largest );

	// assembleCurrentNormalisationCorner().
	mfem::Vector axis( elements ), boundary( elements );
	axis = 0.0;
	boundary = 0.0;
	axis.HostReadWrite();
	boundary.HostReadWrite();
	f.hostLoop( [ & ]( int e, mfem::IntegrationPoint const &ip,
	                   mfem::ElementTransformation &tr, mfem::Vector const &point,
	                   double psi, mfem::Array<int> const &,
	                   mfem::Vector const & )
	{
		double dAxis = 0.0, dBoundary = 0.0;
		BOOST_TEST_REQUIRE( f.source.normalisationDerivatives(
			point( 0 ), point( 1 ), psi, dAxis, dBoundary ) );
		double const w = ip.weight*tr.Weight()/point( 0 );
		axis( e ) += w*dAxis;
		boundary( e ) += w*dBoundary;
	} );

	mfem::Vector second;
	BOOST_TEST_REQUIRE( loops.integrate( meq::PlasmaQuantity::NormalisationCorner,
	                                     f.view, f.constants, f.state,
	                                     f.potentialStart, got, &second ) );
	int const axisOff = differing( got, axis, largest );
	BOOST_TEST( axisOff == 0, axisOff << " per-element d/dpsi_ax partials "
	            "differ from the corner's loop, by up to " << largest );
	int const boundaryOff = differing( second, boundary, largest );
	BOOST_TEST( boundaryOff == 0, boundaryOff << " per-element d/dpsi_bnd "
	            "partials differ from the corner's loop, by up to " << largest );
}

BOOST_AUTO_TEST_CASE( theFourScattersMatchTheirHostLoopsToTheBit )
{
	Fixture f;
	meq::DevicePlasmaLoops loops;
	BOOST_TEST_REQUIRE( loops.prepare( *f.space, f.extraOrder ) );

	double const scale = f.source.currentScale();
	int const ps = f.potentialStart;

	struct Scatter
	{
		meq::PlasmaQuantity quantity;
		char const *host;
	};
	Scatter const scatters[] = {
		{ meq::PlasmaQuantity::CurrentColumn, "assembleCurrentColumn()" },
		{ meq::PlasmaQuantity::CurrentRow, "assembleCurrentRow()" },
		{ meq::PlasmaQuantity::AxisColumn, "assembleNormalisationColumn( axis )" },
		{ meq::PlasmaQuantity::BoundaryColumn,
		  "assembleNormalisationColumn( boundary )" } };

	for ( Scatter const &which : scatters )
	{
		mfem::Vector want( f.state.Size() );
		want = 0.0;
		want.HostReadWrite();
		f.hostLoop( [ & ]( int, mfem::IntegrationPoint const &ip,
		                   mfem::ElementTransformation &tr,
		                   mfem::Vector const &point, double psi,
		                   mfem::Array<int> const &dofs,
		                   mfem::Vector const &shape )
		{
			double factor = 0.0;
			switch ( which.quantity )
			{
				case meq::PlasmaQuantity::CurrentColumn:
				{
					double const derivative =
						f.source.scaledF( point( 0 ), point( 1 ), psi )/scale;
					factor = -ip.weight*tr.Weight()*derivative/point( 0 );
					break;
				}
				case meq::PlasmaQuantity::CurrentRow:
					factor = ip.weight*tr.Weight()
					         *f.source.scaledDFdPsi( point( 0 ), point( 1 ), psi )
					         /point( 0 );
					break;
				default:
				{
					double dAxis = 0.0, dBoundary = 0.0;
					f.source.normalisationDerivatives( point( 0 ), point( 1 ),
					                                   psi, dAxis, dBoundary );
					double const derivative =
						which.quantity == meq::PlasmaQuantity::AxisColumn
						? dAxis : dBoundary;
					double const weight = ip.weight*tr.Weight();
					factor = -weight*derivative/point( 0 );
					break;
				}
			}
			for ( int j = 0; j < dofs.Size(); ++j )
				want( ps + dofs[ j ] ) += factor*shape( j );
		} );

		mfem::Vector got;
		BOOST_TEST_REQUIRE( loops.scatter( which.quantity, f.view, f.constants,
		                                   f.state, ps,
		                                   static_cast<mfem::real_t>( scale ),
		                                   got ) );
		BOOST_TEST_REQUIRE( got.Size() == want.Size() );
		double largest = 0.0;
		int const off = differing( got, want, largest );
		BOOST_TEST( off == 0, off << " of " << want.Size() << " entries differ "
		            "from " << which.host << "'s loop, by up to " << largest );
	}
}

BOOST_AUTO_TEST_CASE( theSmallKernelsAreTheHostLoops )
{
	deviceOnce();

	mfem::Vector source( 40 );
	source.UseDevice( true );
	{
		mfem::real_t *const host = source.HostWrite();
		for ( int i = 0; i < 40; ++i )
			host[ i ] = 0.3 + std::sin( 1.7*i );
	}

	mfem::Array<int> indices( { 3, 17, 4, 39, 0, 22 } );
	indices.UseDevice( true );
	mfem::Vector weights( { 0.25, -1.5, 3.0e-3, 0.7, 2.0, -0.125 } );
	weights.UseDevice( true );

	// gatherOnDevice(): out[ k ] = source[ indices[ k ] ].
	mfem::Vector gathered;
	BOOST_TEST_REQUIRE( meq::gatherOnDevice( source, indices, gathered ) );
	// EVERY BUFFER A KERNEL READ IS FUNNELLED BEFORE THE HOST READS IT, the
	// index list included: under the debug backend a device read leaves the
	// host page protected, and `indices[ k ]` is as raw as `source( i )`.
	gathered.HostRead();
	source.HostRead();
	indices.HostRead();
	weights.HostRead();
	for ( int k = 0; k < indices.Size(); ++k )
		BOOST_TEST( gathered( k ) == source( indices[ k ] ) );

	// orderedDotOnDevice(): limiterValue()'s loop, from zero, in k order.
	mfem::Vector dot;
	BOOST_TEST_REQUIRE( meq::orderedDotOnDevice( source, indices, weights, dot ) );
	source.HostRead();
	indices.HostRead();
	weights.HostRead();
	mfem::real_t want = 0.0;
	for ( int k = 0; k < indices.Size(); ++k )
		want += weights( k )*source( indices[ k ] );
	dot.HostRead();
	BOOST_TEST( dot( 0 ) == want );

	// copyEntriesOnDevice(): the cold state's mask, over a zeroed target.
	mfem::Vector target( 40 );
	target.UseDevice( true );
	target = 0.0;
	BOOST_TEST_REQUIRE( meq::copyEntriesOnDevice( source, indices, target ) );
	target.HostRead();
	source.HostRead();
	indices.HostRead();
	for ( int i = 0; i < 40; ++i )
	{
		bool kept = false;
		for ( int k = 0; k < indices.Size(); ++k )
			kept = kept || indices[ k ] == i;
		BOOST_TEST( target( i ) == ( kept ? source( i ) : 0.0 ) );
	}
}

// THE NODAL FALLBACK'S PEAK IS THE WALK'S: VALUE, INDEX AND SIGN OF ZERO.
//
// peakAt()'s fallback walks the potential block from -inf with a strict `>`,
// so it keeps the FIRST of equal maxima, skips NaN, and reports -1 when
// nothing beats -inf. A parallel reduction gets none of those for free --
// MFEM's ArgMaxReducer keeps whichever equal joined last -- so each is a case
// here, against the walk written out. `+ 0.0` is the walk's own `+
// conductorPsiAtDof()` with no table, and is what makes a -0 peak +0.
BOOST_AUTO_TEST_CASE( theNodalPeakIsTheWalk )
{
	deviceOnce();

	double const inf = std::numeric_limits<double>::infinity();
	double const nan = std::numeric_limits<double>::quiet_NaN();
	struct Case
	{
		char const *name;
		std::vector<double> values;
		int begin;
		int end;
	};
	std::vector<double> ramp( 5000 );
	for ( std::size_t i = 0; i < ramp.size(); ++i )
		ramp[ i ] = std::sin( 0.01*static_cast<double>( i ) );
	std::vector<double> ties = ramp;
	ties[ 4000 ] = ties[ 1200 ] = ties[ 2500 ] = 3.0;
	std::vector<Case> const cases = {
		{ "ramp", ramp, 0, 5000 },
		{ "three equal maxima", ties, 0, 5000 },
		{ "a sub-range missing the first maximum", ties, 1500, 5000 },
		{ "NaN ahead of the maximum", { nan, 1.0, nan, 2.0, nan }, 0, 5 },
		{ "signed zeros, -0 first", { -1.0, -0.0, 0.0, -2.0 }, 0, 4 },
		{ "nothing beats -inf", { -inf, nan, -inf }, 0, 3 },
		{ "an empty range", { 1.0, 2.0 }, 1, 1 },
	};

	for ( Case const &c : cases )
	{
		mfem::Vector source( static_cast<int>( c.values.size() ) );
		source.UseDevice( true );
		{
			mfem::real_t *const host = source.HostWrite();
			for ( std::size_t i = 0; i < c.values.size(); ++i )
				host[ i ] = static_cast<mfem::real_t>( c.values[ i ] );
		}
		source.Read();

		mfem::real_t got = 0.0;
		int gotIndex = -2;
		BOOST_TEST_REQUIRE( meq::nodalPeakOnDevice( source, c.begin, c.end,
		                                            got, gotIndex ) );
		double const reported = gotIndex < 0 ? -inf : static_cast<double>( got ) + 0.0;

		source.HostRead();
		double want = -inf;
		int wantIndex = -1;
		for ( int i = c.begin; i < c.end; ++i )
		{
			double const value = source( i ) + 0.0;
			if ( value > want )
			{
				want = value;
				wantIndex = i;
			}
		}

		std::printf( "  %-40s walk %d %+.17g, device %d %+.17g\n", c.name,
		             wantIndex, want, gotIndex, reported );
		BOOST_TEST( gotIndex == wantIndex );
		BOOST_TEST( std::signbit( reported ) == std::signbit( want ) );
		BOOST_TEST( ( reported == want ) );
	}
}

BOOST_AUTO_TEST_CASE( phaseBRefusesWhatTheScatterCannotHonour )
{
	deviceOnce();

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0 );
	mesh.Transform( shiftOffAxis );

	// H1 shares dofs between elements, so one thread per element writing its
	// own dofs would race. The routes must refuse it and leave the host loop.
	mfem::H1_FECollection continuous( 2, mesh.Dimension() );
	mfem::FiniteElementSpace shared( &mesh, &continuous );
	meq::DevicePlasmaLoops loops;
	BOOST_TEST( !loops.prepare( shared, 4 ) );

	// A vector space likewise: the shape table is scalar.
	mfem::L2_FECollection discontinuous( 2, mesh.Dimension() );
	mfem::FiniteElementSpace vector( &mesh, &discontinuous, 2 );
	meq::DevicePlasmaLoops vectorLoops;
	BOOST_TEST( !vectorLoops.prepare( vector, 4 ) );

	// And the scalar L2 space is accepted -- the refusals are not blanket.
	mfem::FiniteElementSpace scalar( &mesh, &discontinuous );
	meq::DevicePlasmaLoops scalarLoops;
	BOOST_TEST( scalarLoops.prepare( scalar, 4 ) );
}
