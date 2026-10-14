// Unit tests for meq::DeviceRows -- meq::CompressedRows' dot products, on the
// device, in one launch.
//
// THIS FILE ASSERTS EQUALITY TO THE BIT AND THAT IS THE WHOLE POINT OF IT.
// meq::CompressedRows documents a bit-exactness contract in terms: "summing the
// surviving terms in ASCENDING INDEX ORDER reproduces the dense loop's every
// partial sum". meq::DeviceRows exists to do that sum somewhere else, and the
// only thing that makes the two interchangeable is that they agree exactly. A
// tolerance here would hide precisely the edit this file exists to catch -- a
// tree reduction dropped into the kernel for occupancy, which is a reasonable
// thing to want and which would quietly cost MEQ the digit-for-digit comparison
// between --device cuda and --device cpu that root-caused M-79, M-130 and
// M-187. If these BOOST_TEST( a == b ) ever become BOOST_TEST( a == b,
// tolerance( ... ) ), meq::DeviceRows has stopped being what its header claims.
//
// AND IT RUNS ON A MACHINE WITH NO GPU, which is what makes it a ctest rather
// than a thing somebody runs by hand. mfem::Device( "debug" ) has device memory
// semantics with host arithmetic, so Device::Allows( DEVICE_MASK ) is true, the
// kernel is dispatched, the memory is moved and protected exactly as CUDA would
// -- and none of it needs the card. DEVICE-PLAN.md section 4 prescribes that
// instrument and this is the first test to take it.

#define BOOST_TEST_MODULE DeviceRowsTests
#include <boost/test/included/unit_test.hpp>

#include "convergence/RoundOffFloor.hpp"
#include "meq/Crossings.hpp"
#include "meq/DeviceRows.hpp"
#include "meq/GradShafranov.hpp"

#include "mfem.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace
{
	/*
	 * ONE Device FOR THE BINARY, AND IT MUST BE A FUNCTION-LOCAL STATIC RATHER
	 * THAN A NAMESPACE-SCOPE ONE.
	 *
	 * mfem::Device::Get() returns `device_singleton`, a static member defined
	 * in general/device.cpp, and Setup() opens MFEM_VERIFY( ngpu == -1 ). A
	 * namespace-scope Device in THIS translation unit may be constructed before
	 * that singleton's own dynamic initialisation runs, and then `ngpu` is
	 * still its zero-initialised 0 rather than the -1 the constructor writes --
	 * so MFEM reports "the mfem::Device is already configured!" on the FIRST
	 * and only Device in the process. Measured here, as an abort before the
	 * first test case.
	 *
	 * A function-local static is initialised on first call, which is after all
	 * static initialisation has finished, so the order cannot bite. This is the
	 * ordinary static-initialisation-order fiasco wearing a device costume, and
	 * it is worth the comment because the error message points at the wrong
	 * thing entirely.
	 */
	mfem::Device &deviceOnce()
	{
		/*
		 * "debug" BY DEFAULT AND CUDA ON REQUEST, because the two answer
		 * different questions and only one of them can run in CI. `debug` has
		 * device memory semantics with HOST arithmetic, so it proves the
		 * bookkeeping and needs no card; `cuda` proves the kernel. The
		 * environment variable is how one binary does both, since mfem::Device
		 * is process-wide and a test binary gets one:
		 *
		 *     MEQ_TEST_DEVICE=cuda build/tests/DeviceRowsTests
		 *
		 * The ctest registers the default. **A bit-exact claim has to be made
		 * on the card as well as on the debug backend** -- host arithmetic
		 * cannot tell you whether an FMA contracted differently in a kernel --
		 * so the cuda run is part of this file's acceptance even though CI
		 * cannot take it.
		 */
		char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
		static mfem::Device device( requested ? requested : "debug" );
		return device;
	}

	/// A support and a row set with structure rather than noise in them: the
	/// support is sparse and ascending as the class requires, and the values
	/// span several orders so that a reassociated sum would show.
	struct Fixture
	{
		std::vector<int> support;
		std::vector<double> values;
		int rows = 0;
		int length = 0;

		Fixture( int rowCount, int width, int stride, int vectorLength )
			: rows( rowCount ), length( vectorLength )
		{
			for ( int k = 0; k < width; ++k )
				support.push_back( k*stride );

			values.reserve( static_cast<std::size_t>( rowCount )*width );
			for ( int r = 0; r < rowCount; ++r )
				for ( int k = 0; k < width; ++k )
					values.push_back(
						std::pow( 10.0, ( ( r + k ) % 9 ) - 4 )
						* ( ( ( r*width + k ) % 7 ) - 3 ) );
		}

		/// The dense rows meq::CompressedRows::compress() wants.
		std::vector<mfem::Vector> dense() const
		{
			std::vector<mfem::Vector> out;
			for ( int r = 0; r < rows; ++r )
			{
				mfem::Vector row( length );
				row = 0.0;
				for ( std::size_t k = 0; k < support.size(); ++k )
					row( support[ k ] ) =
						values[ static_cast<std::size_t>( r )*support.size() + k ];
				out.push_back( row );
			}
			return out;
		}
	};

	/// A vector with no exact zeros and a wide dynamic range, so that a
	/// difference in summation order has something to bite on.
	mfem::Vector probe( int n )
	{
		mfem::Vector v( n );
		v.HostWrite();
		for ( int i = 0; i < n; ++i )
			v( i ) = std::sin( 0.7*i + 0.3 )*std::pow( 10.0, ( i % 11 ) - 5 );
		return v;
	}
}


/*
 * THE CENTRAL CLAIM. Every row, both routes, exact equality.
 */
BOOST_AUTO_TEST_CASE( theDeviceRowsReproduceTheHostDotToTheBit )
{
	deviceOnce();
	for ( auto const &shape : { Fixture( 14, 97, 13, 2048 ),
	                            Fixture( 1, 1, 1, 8 ),
	                            Fixture( 3, 512, 4, 4096 ) } )
	{
		meq::CompressedRows host;
		host.compress( shape.dense(), shape.support );

		meq::DeviceRows device;
		device.upload( host.support(), host.values(), host.rowCount() );

		BOOST_TEST_REQUIRE( device.valid() );
		BOOST_TEST( device.rowCount() == host.rowCount() );
		BOOST_TEST( device.width() == static_cast<int>( host.support().size() ) );

		mfem::Vector const v = probe( shape.length );

		mfem::Vector got;
		BOOST_TEST_REQUIRE( device.dotAll( v, got ),
		                    "the kernel refused under mfem::Device( \"debug\" ), "
		                    "which means Device::Allows( DEVICE_MASK ) is false "
		                    "there -- the test's instrument is gone, not the "
		                    "code under test" );
		BOOST_TEST_REQUIRE( got.Size() == host.rowCount() );
		// The kernel leaves its result device-valid; the reader funnels.
		got.HostRead();

		std::vector<double> scratch;
		for ( int r = 0; r < host.rowCount(); ++r )
		{
			double const expected = host.dot( r, v, scratch );
			/*
			 * EXACTLY EQUAL IN DOUBLE, AND IN SINGLE THE TWO ROUTES ARE NOT
			 * THE SAME COMPUTATION AT ALL. `meq::CompressedRows` stores
			 * `std::vector<double>` and accumulates in `double`;
			 * `meq::DeviceRows` stores `mfem::real_t` and accumulates in
			 * `float`. Under `MFEM_PRECISION = single` those differ by the
			 * format before any question of ORDER arises, so the bit-exactness
			 * contract is not available and no edit to the kernel could make
			 * it hold.
			 *
			 * `exactAgreementFloor()` is therefore EXACTLY ZERO in double,
			 * where this remains the bit-equality it has always been, and
			 * 64 eps in single. **What that gives up is the reassociation
			 * check** -- a reassociated sum and a float-vs-double sum are the
			 * same order -- and what it keeps is the layout: a transposed
			 * `r + nr*c`, a wrong support index or a dropped row is `O( 1 )`.
			 *
			 * The reassociation half is a DOUBLE-BUILD acceptance now, and that
			 * is a real narrowing rather than a tidy one. Making it hold in
			 * both would mean giving `meq::CompressedRows` the same `real_t`
			 * its mirror carries -- a production change, bit-identical in
			 * double, which nobody has measured in single. M-219.
			 */
			BOOST_TEST( std::abs( got( r ) - expected )
			            <= meq::tests::exactAgreementFloor( expected ),
			            "row " << r << " of " << host.rowCount()
			            << ": device " << got( r ) << " against host "
			            << expected << ". These must be EQUAL, not close -- "
			            "see this file's opening comment. A near miss is a "
			            "reassociated sum, which is the one change the kernel "
			            "may not make." );
		}
	}
}


/*
 * AND THE REFERENCE MUST NOT BE TRIVIAL, or the case above passes on zeros.
 * The same shape as meq::CompressedRows' own acceptance, which checks its
 * reference integral is not identically zero for exactly this reason.
 */
BOOST_AUTO_TEST_CASE( theReferenceDotsAreNotAllZero )
{
	deviceOnce();
	Fixture const shape( 14, 97, 13, 2048 );
	meq::CompressedRows host;
	host.compress( shape.dense(), shape.support );

	mfem::Vector const v = probe( shape.length );
	std::vector<double> scratch;

	double worst = 0.0;
	for ( int r = 0; r < host.rowCount(); ++r )
		worst = std::max( worst, std::abs( host.dot( r, v, scratch ) ) );

	BOOST_TEST( worst > 1.0e-6,
	            "every reference dot is at round-off (" << worst << "), so the "
	            "equality above is vacuous" );
}


/*
 * THE REFUSALS. Each is a read off the end of a device buffer taking a
 * different route there, and under CUDA an off-the-end read is not a fault --
 * it is whatever was in that memory, arriving as a plausible border row.
 */
BOOST_AUTO_TEST_CASE( theMirrorRefusesWhatItCannotDescribe )
{
	deviceOnce();
	meq::DeviceRows device;

	// Values that are not rowCount * width long.
	BOOST_CHECK_THROW( device.upload( { 0, 2, 4 }, { 1.0, 2.0 }, 1 ),
	                   std::invalid_argument );

	// A negative row count.
	BOOST_CHECK_THROW( device.upload( { 0, 2, 4 }, {}, -1 ),
	                   std::invalid_argument );

	// An empty row set is not an error -- it is the free-boundary case with no
	// exterior modes -- and it must refuse the launch rather than run an empty
	// one, so the caller keeps the host path.
	device.upload( {}, {}, 0 );
	BOOST_TEST( !device.valid() );
	mfem::Vector out;
	mfem::Vector const v = probe( 16 );
	BOOST_TEST( !device.dotAll( v, out ),
	            "an empty mirror must refuse rather than report success over "
	            "no rows" );

	// A support that reaches past the vector it is contracted against: rows
	// compressed on one space, used on another.
	device.upload( { 0, 4, 40 }, { 1.0, 2.0, 3.0 }, 1 );
	BOOST_TEST_REQUIRE( device.valid() );
	mfem::Vector const tooShort = probe( 8 );
	BOOST_CHECK_THROW( device.dotAll( tooShort, out ), std::invalid_argument );
}


/*
 * THE BORDERED STEP'S ( N + 4 )^2 FILL, IN ONE LAUNCH.
 *
 * The host writes `dense( i, j ) = cornerEntry( i, j ) - rowDot( i, columnZ( j ) )`
 * with one gather per column and a dot per pair. This is the dot per pair, and
 * the case asserts the same thing the single-column one does -- equality to the
 * bit against meq::CompressedRows::dot() -- plus the one property that only
 * exists here: the ORIENTATION. A transposed bordered matrix factors and solves
 * perfectly well, to the wrong step, so a test that only checked the set of
 * entries would pass against it.
 */
BOOST_AUTO_TEST_CASE( theColumnFillMatchesTheHostAndIsNotTransposed )
{
	deviceOnce();

	Fixture const shape( 14, 97, 13, 2048 );
	int const columnCount = 18;   // N + 4 on the shipped machine case

	meq::CompressedRows host;
	host.compress( shape.dense(), shape.support );

	meq::DeviceRows device;
	device.upload( host.support(), host.values(), host.rowCount() );
	BOOST_TEST_REQUIRE( device.valid() );

	// The columns as the bordered step produces them: one batch, contiguous,
	// column c at offset c*length. See M-98 -- ArrayMult already hands them
	// over this way.
	mfem::Vector packed( columnCount*shape.length );
	packed.HostWrite();
	std::vector<mfem::Vector> separate;
	for ( int c = 0; c < columnCount; ++c )
	{
		mfem::Vector column( shape.length );
		column.HostWrite();
		for ( int i = 0; i < shape.length; ++i )
		{
			double const v = std::cos( 0.37*i + 1.1*c )
			                 *std::pow( 10.0, ( ( i + c ) % 9 ) - 4 );
			column( i ) = v;
			packed( c*shape.length + i ) = v;
		}
		separate.push_back( column );
	}

	mfem::DenseMatrix got;
	BOOST_TEST_REQUIRE( device.dotColumns( packed, columnCount, shape.length,
	                                       got ) );
	BOOST_TEST_REQUIRE( got.Height() == host.rowCount() );
	BOOST_TEST_REQUIRE( got.Width() == columnCount );
	// The kernel leaves its result device-valid; the reader funnels, through
	// the matrix's OWN memory rather than a Vector wrapped round its pointer.
	got.HostRead();

	std::vector<double> scratch;
	double largest = 0.0;
	for ( int r = 0; r < host.rowCount(); ++r )
		for ( int c = 0; c < columnCount; ++c )
		{
			double const expected = host.dot( r, separate[ c ], scratch );
			largest = std::max( largest, std::abs( expected ) );
			// Exactly zero in double; 64 eps in single, where the host sums in
			// `double` and the device in `float`. See the note above. M-219.
			BOOST_TEST( std::abs( got( r, c ) - expected )
			            <= meq::tests::exactAgreementFloor( expected ),
			            "entry ( " << r << ", " << c << " ) : device "
			            << got( r, c ) << " against host " << expected
			            << ". Equal, not close -- and if the whole matrix looks "
			            "transposed, that is the r + nr*c layout, not the "
			            "arithmetic." );
		}

	BOOST_TEST( largest > 1.0e-6,
	            "every entry is at round-off, so the equality is vacuous" );

	/*
	 * THE ORIENTATION, ASSERTED SEPARATELY. The fixture's rows differ from each
	 * other and the columns differ from each other, so a transpose would show
	 * -- but only if something looks for it. This asks that entry ( 0, 1 ) is
	 * row 0 against column 1 and NOT row 1 against column 0.
	 */
	double const zeroOne = host.dot( 0, separate[ 1 ], scratch );
	double const oneZero = host.dot( 1, separate[ 0 ], scratch );
	BOOST_TEST_REQUIRE( zeroOne != oneZero,
	                    "the fixture cannot distinguish a transpose" );
	// The transpose this separates is O( 1 ) -- zeroOne and oneZero are
	// required above to differ -- so round-off cannot blur it. M-219.
	BOOST_TEST( std::abs( got( 0, 1 ) - zeroOne )
	            <= meq::tests::exactAgreementFloor( zeroOne ),
	            "entry ( 0, 1 ) must be row 0 against column 1" );
	BOOST_TEST( std::abs( got( 1, 0 ) - oneZero )
	            <= meq::tests::exactAgreementFloor( oneZero ),
	            "entry ( 1, 0 ) must be row 1 against column 0" );
}


/*
 * THE CROSSING COUNTER COUNTS TRANSFERS AND NOT CALLS, AND THAT DISTINCTION IS
 * THE WHOLE INSTRUMENT -- DEVICE-PLAN.md §4 stage 4b.
 *
 * A funnel call moves data only when the device holds the live copy; once any
 * read has brought a buffer to the host it stays valid until the device writes
 * again. A counter that could not tell those apart would report a number 16x
 * too large and name it a crossing count, which is exactly what §3.3's
 * leg-call proxy did.
 */
BOOST_AUTO_TEST_CASE( theCounterSeparatesAReadFromATransfer )
{
	deviceOnce();

	meq::crossings().clear();

	mfem::Vector v( 64 );
	v.UseDevice( true );
	v.HostWrite();
	for ( int i = 0; i < v.Size(); ++i )
		v( i ) = 0.5*i;

	// Host-valid: a read moves nothing.
	meq::countedHostRead( v );
	BOOST_TEST( meq::crossings().reads == 1 );
	BOOST_TEST( meq::crossings().transfers == 0,
	            "a read of a host-valid buffer moved data" );

	// Make the device the live copy, the way a kernel would.
	v.Write();

	meq::countedHostRead( v );
	BOOST_TEST( meq::crossings().reads == 2 );
	BOOST_TEST( meq::crossings().transfers == 1,
	            "a read of a device-valid buffer did not count as a transfer" );
	BOOST_TEST( meq::crossings().entries == 64,
	            "the transfer did not account for the buffer's size" );

	// And it is now host-valid again, so a second read is free. THIS is the
	// property that makes a call count an overestimate: without it the two
	// reads above would look like two crossings.
	meq::countedHostRead( v );
	BOOST_TEST( meq::crossings().reads == 3 );
	BOOST_TEST( meq::crossings().transfers == 1,
	            "a repeat read counted as a second transfer, so the counter is "
	            "measuring calls rather than data movement" );

	meq::crossings().clear();
	BOOST_TEST( meq::crossings().reads == 0 );
}

/*
 * AND IT SEES THE OTHER DIRECTION, WHICH IT COULD NOT UNTIL PHASE B NEEDED IT.
 *
 * DEVICE-PLAN.md section 3.3 defines phase B as *one upload*, and the counter
 * above measures host reads -- downloads. **An instrument that can see one
 * direction of the bus cannot measure a plan whose first phase is an upload**,
 * which is the same species of coverage defect as M-198's first crossing
 * counter covering 17 funnels of 36 and reporting a number 18x too low.
 *
 * The condition is the MIRROR and not the same one: a device read moves data
 * exactly when the host holds the live copy and the device does not.
 */
BOOST_AUTO_TEST_CASE( theCounterSeesTheUploadDirection )
{
	deviceOnce();

	meq::crossings().clear();

	mfem::Vector v( 64 );
	v.UseDevice( true );
	v.HostWrite();
	for ( int i = 0; i < v.Size(); ++i )
		v( i ) = 0.5*i;

	// Host-valid and device-stale: a device read has to move it.
	meq::countedDeviceRead( v );
	BOOST_TEST( meq::crossings().deviceReads == 1 );
	BOOST_TEST( meq::crossings().uploads == 1,
	            "a device read of a host-only buffer did not count as an upload" );
	BOOST_TEST( meq::crossings().uploadEntries == 64,
	            "the upload did not account for the buffer's size" );

	// And now both copies are valid, so a second device read is free. THIS is
	// what makes phase B's "one upload" an assertable property rather than a
	// hope: a constant that stays up is re-read for nothing.
	meq::countedDeviceRead( v );
	BOOST_TEST( meq::crossings().deviceReads == 2 );
	BOOST_TEST( meq::crossings().uploads == 1,
	            "a repeat device read counted as a second upload, so the "
	            "counter is measuring calls rather than data movement" );

	// A host WRITE takes the device copy away again, which is what a phase-C
	// host leg does to phase B's payload -- and the next upload is real.
	v.HostWrite();
	meq::countedDeviceRead( v );
	BOOST_TEST( meq::crossings().uploads == 2,
	            "a device read after a host write did not re-upload" );

	// The two directions are counted separately and neither touches the other.
	BOOST_TEST( meq::crossings().reads == 0,
	            "a device read was counted as a host read" );
	BOOST_TEST( meq::crossings().transfers == 0 );

	meq::crossings().clear();
}

/*
 * THE PHASE MARKER RESTORES WHAT IT FOUND, WHICH IS WHAT MAKES IT SAFE UNDER
 * RE-ENTRY.
 *
 * meq::GradShafranovSolver::solve() is re-entrant -- solveByPicardThenNewton()
 * calls it for its second stage -- and prepare(), which is phase B, is
 * re-entered from inside it by picardStep(). A marker that cleared to a default
 * on exit would have the inner scope leave the outer one mislabelled for the
 * rest of its run, so every crossing after it would be attributed to the wrong
 * phase and the headline would silently move. Same hazard solve() records about
 * `profileStart` being a local rather than a member, one level up.
 */
BOOST_AUTO_TEST_CASE( thePhaseScopeRestoresWhatItFound )
{
	BOOST_TEST_REQUIRE( ( meq::residencyPhase() == meq::Phase::Setup ) );

	{
		meq::PhaseScope const outer( meq::Phase::Solve );
		BOOST_TEST( ( meq::residencyPhase() == meq::Phase::Solve ) );

		{
			meq::PhaseScope const inner( meq::Phase::Upload );
			BOOST_TEST( ( meq::residencyPhase() == meq::Phase::Upload ) );
		}

		// The point of the case: Solve, not Setup.
		BOOST_TEST( ( meq::residencyPhase() == meq::Phase::Solve ),
		            "a nested phase scope cleared to the default instead of "
		            "restoring the phase it found, so every crossing after an "
		            "inner scope would be attributed to the wrong phase" );
	}

	BOOST_TEST( ( meq::residencyPhase() == meq::Phase::Setup ) );
}

/*
 * AND A CROSSING LANDS IN THE PHASE IT HAPPENED IN, which is the whole reason
 * the phases exist: a transfer in phase B is the design and the same transfer
 * in phase C is the defect, and one total cannot tell them apart.
 */
BOOST_AUTO_TEST_CASE( aCrossingIsAttributedToItsPhase )
{
	deviceOnce();

	meq::crossings().clear();
	for ( int k = 0; k < meq::phaseCount; ++k )
		meq::crossingPhases()[ k ].clear();

	mfem::Vector v( 32 );
	v.UseDevice( true );
	v.HostWrite();

	// FILLED IN A LOOP AND NOT WITH `v = 1.0`, AND THE FIRST DRAFT OF THIS CASE
	// USED THE ASSIGNMENT AND FAILED. mfem::Vector::operator=( double ) is
	// DEVICE-AWARE -- Write() then a forall -- so it leaves the vector
	// device-valid, and the countedDeviceRead() below then correctly counts no
	// upload because there is nothing to upload. That is M-200's finding, which
	// cost MEQ 12x its transfer count in production code, reproduced here by a
	// test written to measure uploads.
	for ( int i = 0; i < v.Size(); ++i )
		v( i ) = 1.0;

	{
		meq::PhaseScope const phase( meq::Phase::Upload );
		meq::countedDeviceRead( v );
	}

	v.Write();

	{
		meq::PhaseScope const phase( meq::Phase::Download );
		meq::countedHostRead( v );
	}

	int const upload = static_cast<int>( meq::Phase::Upload );
	int const download = static_cast<int>( meq::Phase::Download );
	int const solve = static_cast<int>( meq::Phase::Solve );

	BOOST_TEST( meq::crossingPhases()[ upload ].uploads == 1 );
	BOOST_TEST( meq::crossingPhases()[ download ].transfers == 1 );
	BOOST_TEST( meq::crossingPhases()[ solve ].uploads == 0,
	            "a crossing outside phase C was attributed to it" );
	BOOST_TEST( meq::crossingPhases()[ solve ].transfers == 0 );

	meq::crossings().clear();
	for ( int k = 0; k < meq::phaseCount; ++k )
		meq::crossingPhases()[ k ].clear();
}

// THE MIRROR IS REFRESHED ONCE PER NEWTON STEP, SO upload() IS CALLED MANY
// TIMES ON ONE OBJECT -- AND EVERY CALL AFTER THE FIRST WRITES INTO BUFFERS ITS
// OWN PREVIOUS CALL LEFT DEVICE-VALID.
//
// The first upload is the easy case and the only one that survives a missing
// HostWrite(): a freshly constructed mfem::Vector has UseDevice() false, so its
// fill is an ordinary host loop. After the closing Read() both members are
// device-valid, SetSize() at an unchanged size reallocates nothing, and the
// next fill is a RAW HOST WRITE INTO A PROTECTED PAGE. Under a debug device
// that aborts; under CUDA it writes the host copy while the kernel goes on
// reading the stale device one, which is a wrong border row arriving as a
// plausible number.
//
// It is exactly the shape of MEASUREMENTS.md M-187 one level down -- "mode 0
// survives, every mode after it reads zero" -- and this case is here because
// the first version of upload() had it and no single-upload test could tell.
BOOST_AUTO_TEST_CASE( theMirrorSurvivesBeingRefreshed )
{
	deviceOnce();

	Fixture const shape( 6, 50, 9, 512 );   // ( width - 1 )*stride = 441 < 512
	meq::CompressedRows host;
	host.compress( shape.dense(), shape.support );

	mfem::Vector const v = probe( shape.length );

	meq::DeviceRows device;
	mfem::Vector first;
	device.upload( host.support(), host.values(), host.rowCount() );
	BOOST_TEST_REQUIRE( device.dotAll( v, first ) );
	first.HostRead();

	// FOUR MORE, as a bordered solve would do, with the SAME rows: the answer
	// must not move, and the abort must not happen.
	for ( int again = 0; again < 4; ++again )
	{
		device.upload( host.support(), host.values(), host.rowCount() );

		mfem::Vector repeat;
		BOOST_TEST_REQUIRE( device.dotAll( v, repeat ),
		                    "the kernel refused on refresh " 
		                        + std::to_string( again + 1 ) );
		repeat.HostRead();

		BOOST_TEST_REQUIRE( repeat.Size() == first.Size() );
		for ( int r = 0; r < repeat.Size(); ++r )
			BOOST_TEST( repeat( r ) == first( r ),
			            "row " + std::to_string( r ) + " moved on refresh "
			                + std::to_string( again + 1 ) + ": "
			                + std::to_string( repeat( r ) ) + " against "
			                + std::to_string( first( r ) )
			                + " -- the mirror was refreshed into a buffer the "
			                  "previous upload left device-valid" );
	}
}

// EVERY BORDER ROW AGAINST EVERY COLUMN IN ONE LAUNCH, EQUAL TO THE HOST
// rowDot()'S ARMS TO THE BIT -- the bordered fill's device route, M-238.
//
// The arms are written here as GradShafranov.cpp's rowDot() writes them: a
// compressed row's sum from zero in support order; psi_ax's in dof order;
// psi_bnd's NEGATED after the sum; XP-3's by `total -=`. XP-3's sign folds into
// its weights exactly; psi_bnd's cannot, because `-( +0 )` is `-0` and a sum of
// negated terms is `+0` -- so it is meq::RowView's flag, and the zero column
// below is what tells the two apart.
BOOST_AUTO_TEST_CASE( theBatchedFillIsEveryHostArmToTheBit )
{
	deviceOnce();

	int const n = 61;
	int const columnCount = 4;
	std::vector<mfem::Vector> columns( columnCount );
	for ( int c = 0; c < columnCount; ++c )
	{
		columns[ c ].SetSize( n );
		columns[ c ].UseDevice( true );
		mfem::real_t *const host = columns[ c ].HostWrite();
		for ( int i = 0; i < n; ++i )
			host[ i ] = c == columnCount - 1
			            ? 0.0   // a zero column, as XP-3's z_j is
			            : std::sin( 0.7*i + 1.3*c ) + 0.25*c;
		columns[ c ].Read();
	}

	// A compressed row set on a shared ascending support, through the mirror.
	std::vector<int> support = { 2, 5, 9, 14, 20, 33, 47, 58 };
	std::vector<double> values;
	for ( int r = 0; r < 2; ++r )
		for ( std::size_t k = 0; k < support.size(); ++k )
			values.push_back( std::cos( 0.9*k + 2.1*r ) - 0.3 );
	meq::DeviceRows mirror;
	mirror.upload( support, values, 2 );

	// The short rows, in their own dof orders, with rowDot()'s signs folded.
	std::vector<int> shortIndices = { 40, 3, 17, 11,   7, 52, 30 };
	std::vector<double> shortWeights = { 0.5, -1.25, 3.0e-3, 0.7,   2.0, -0.125, 1.5 };
	mfem::Array<int> shortIndexTable( static_cast<int>( shortIndices.size() ) );
	mfem::Vector shortWeightTable( static_cast<int>( shortWeights.size() ) );
	shortIndexTable.UseDevice( true );
	shortWeightTable.UseDevice( true );
	for ( std::size_t k = 0; k < shortIndices.size(); ++k )
	{
		shortIndexTable[ static_cast<int>( k ) ] = shortIndices[ k ];
		shortWeightTable( static_cast<int>( k ) ) = shortWeights[ k ];
	}
	// Row 0 is psi_ax's, the first four; row 1 psi_bnd's, the last three with
	// the finished sum negated -- the flag, not the weights. See meq::RowView.
	int const *const dIndices = shortIndexTable.Read();
	mfem::real_t const *const dWeights = shortWeightTable.Read();

	mfem::Array<meq::RowView> rows( 4 );
	rows.UseDevice( true );
	rows[ 0 ] = meq::RowView{ dIndices, dWeights, 4, 0 };
	rows[ 1 ] = meq::RowView{ dIndices + 4, dWeights + 4, 3, 1 };
	rows[ 2 ] = mirror.row( 0 );
	rows[ 3 ] = mirror.row( 1 );

	mfem::Array<mfem::real_t const *> columnPointers;
	columnPointers.UseDevice( true );
	for ( mfem::Vector const &c : columns )
		columnPointers.Append( c.Read() );

	mfem::Vector dots;
	BOOST_TEST_REQUIRE( meq::dotRowsAgainstColumnsOnDevice( rows, columnPointers,
	                                                        dots ) );
	dots.HostRead();
	BOOST_TEST_REQUIRE( dots.Size() == 4*columnCount );

	for ( int c = 0; c < columnCount; ++c )
	{
		mfem::Vector const &v = columns[ c ];
		v.HostRead();

		// psi_ax's arm. Each arm also keeps the sum of its terms' magnitudes,
		// which is what its round-off scales with when the terms cancel.
		double axis = 0.0;
		double magnitude[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
		for ( int k = 0; k < 4; ++k )
		{
			axis += shortWeights[ k ]*v( shortIndices[ k ] );
			magnitude[ 0 ] += std::abs( shortWeights[ k ]*v( shortIndices[ k ] ) );
		}
		// psi_bnd's: limiterValue(), then negated.
		mfem::real_t limiter = 0.0;
		for ( int k = 4; k < 7; ++k )
		{
			limiter += shortWeights[ k ]*v( shortIndices[ k ] );
			magnitude[ 1 ] += std::abs( shortWeights[ k ]*v( shortIndices[ k ] ) );
		}
		double const boundary = -limiter;
		// The compressed rows: CompressedRows::dot() on the gathered column.
		double compressed[ 2 ] = { 0.0, 0.0 };
		for ( int r = 0; r < 2; ++r )
			for ( std::size_t k = 0; k < support.size(); ++k )
			{
				compressed[ r ] += values[ r*support.size() + k ]*v( support[ k ] );
				magnitude[ 2 + r ] += std::abs( values[ r*support.size() + k ]*v( support[ k ] ) );
			}

		/*
		 * TO THE BIT IN DOUBLE, AND TO THE FORMAT IN SINGLE, for the reason
		 * this file's first case gives: the host arms accumulate in double, and
		 * psi_ax's and the compressed rows' weights are double, where the
		 * launch is real_t throughout. In single the two are different
		 * arithmetic and part by about one eps of the row (measured 1.5e-08 to
		 * 3.2e-08 on entries near 0.25 to 0.8). exactAgreementFloor() is
		 * exactly zero in double, so there the check, sign bit included, is
		 * the bit equality it always was. M-273.
		 */
		double const want[ 4 ] = { axis, boundary, compressed[ 0 ], compressed[ 1 ] };
		for ( int r = 0; r < 4; ++r )
		{
			double const got = dots( r + 4*c );
			double const floor = meq::tests::exactAgreementFloor( magnitude[ r ] );
			BOOST_TEST( ( floor == 0.0
			              ? ( got == want[ r ]
			                  && std::signbit( got ) == std::signbit( want[ r ] ) )
			              : std::abs( got - want[ r ] ) <= floor ),
			            "row " << r << " against column " << c << ": the launch "
			            "gives " << got << " and the host arm " << want[ r ] );
		}
	}

	// AND THE XP-3 ARM, `total -=`, against the same table read as +w: one
	// more launch with the X-point's own sign convention folded.
	mfem::Vector xpointWeights( 3 );
	xpointWeights.UseDevice( true );
	for ( int k = 0; k < 3; ++k )
		xpointWeights( k ) = -shortWeights[ 4 + k ];
	mfem::Array<meq::RowView> xpointRow( 1 );
	xpointRow.UseDevice( true );
	xpointRow[ 0 ] = meq::RowView{ shortIndexTable.Read() + 4,
	                               xpointWeights.Read(), 3, 0 };
	// FRESH POINTERS, TAKEN NOW: the loop above read every column on the
	// host, and a device pointer is good for the launch it was taken for --
	// under the debug backend a stale one is a fault, which is how this case
	// first failed. The solver takes its own immediately before its launch.
	mfem::Array<mfem::real_t const *> freshPointers;
	freshPointers.UseDevice( true );
	for ( mfem::Vector const &c : columns )
		freshPointers.Append( c.Read() );
	BOOST_TEST_REQUIRE( meq::dotRowsAgainstColumnsOnDevice( xpointRow,
	                                                        freshPointers,
	                                                        dots ) );
	dots.HostRead();
	for ( int c = 0; c < columnCount; ++c )
	{
		mfem::Vector const &v = columns[ c ];
		v.HostRead();
		double total = 0.0;
		double magnitude = 0.0;
		for ( int k = 4; k < 7; ++k )
		{
			total -= shortWeights[ k ]*v( shortIndices[ k ] );
			magnitude += std::abs( shortWeights[ k ]*v( shortIndices[ k ] ) );
		}
		// The host's double against the launch's real_t, as above. M-273.
		double const floor = meq::tests::exactAgreementFloor( magnitude );
		BOOST_TEST( ( floor == 0.0
		              ? ( dots( c ) == total
		                  && std::signbit( dots( c ) ) == std::signbit( total ) )
		              : std::abs( dots( c ) - total ) <= floor ),
		            "XP-3's row against column " << c << ": " << dots( c )
		            << " against the host's " << total );
	}
}
