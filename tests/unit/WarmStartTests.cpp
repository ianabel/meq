// Unit tests for the warm start's threading under an mfem::Device:
// meq::setupLoopMayThread() and the two loops it gates, meq::FieldTransfer's
// projection and meq::PointLocator's host search.
//
// A warm start is setup, host work on every arm, so under a Device it threads
// where it provably cannot reach MFEM's memory manager, which is not
// reentrant. This file holds both halves of that: the predicate refuses
// exactly what it should, and where it allows threads the transfer is the
// one-thread transfer to the bit.
//
// "debug" by default and cuda on request, as tests/unit/DeviceBasisTests.cpp:
//
//     MEQ_TEST_DEVICE=cuda build/tests/WarmStartTests
//
// The two arms answer different questions. `debug`'s host memory type is
// HOST_DEBUG, which registers every allocation, so the predicate must refuse
// and the transfer runs on one thread. Only `cuda` (plain HOST) threads, so
// only the cuda arm exercises the threaded loops under a Device, and ctest
// registers it as WarmStartTests_cuda wherever MFEM has CUDA.

#define BOOST_TEST_MODULE WarmStart
#include <boost/test/unit_test.hpp>

#include "meq/Threading.hpp"
#include "meq/WarmStart.hpp"

#include "mfem.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>

#ifdef MFEM_USE_OPENMP
	#include <omp.h>
#endif

namespace
{
	/// One Device for the binary, a function-local static for the reason
	/// tests/unit/DeviceRowsTests.cpp gives.
	mfem::Device &deviceOnce()
	{
		char const *const requested = std::getenv( "MEQ_TEST_DEVICE" );
		static mfem::Device device( requested ? requested : "debug" );
		return device;
	}

	/// Whether this build threads MEQ's loops at all.
	constexpr bool threadingBuilt()
	{
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
		return true;
#else
		return false;
#endif
	}

	/// Whether the Device's host memory is plain new[], which is the one case
	/// meq::setupLoopMayThread() lets through.
	bool hostMemoryIsPlain()
	{
		return mfem::MemoryManager::GetHostMemoryType() == mfem::MemoryType::HOST;
	}

	/// A guess mesh of QUADRILATERALS with the interior vertices moved, so every
	/// element is bilinear and meq::PointLocator takes its host search: the
	/// device search takes straight-sided triangles only.
	mfem::Mesh guessMesh()
	{
		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			8, 6, mfem::Element::QUADRILATERAL, false, 1.7, 0.9 );
		for ( int v = 0; v < mesh.GetNV(); ++v )
		{
			mfem::real_t *const point = mesh.GetVertex( v );
			double const x = point[ 0 ];
			double const y = point[ 1 ];
			// Zero on the boundary of the box, so the two meshes cover the
			// same domain and no point misses.
			double const bump = std::sin( M_PI*x/1.7 )*std::sin( M_PI*y/0.9 );
			point[ 0 ] = static_cast<mfem::real_t>( x + 0.04*bump );
			point[ 1 ] = static_cast<mfem::real_t>( y - 0.03*bump );
		}
		return mesh;
	}

	double guessValue( mfem::Vector const &x )
	{
		return std::sin( 2.0*x( 0 ) )*std::cos( 3.0*x( 1 ) ) + 0.5*x( 0 )*x( 1 );
	}
}

// THE PREDICATE REFUSES A REGISTERED OBJECT AND A REGISTERING HOST MEMORY.
//
// A vector the device has written is registered with the memory manager, so a
// loop reading it through an MFEM accessor would reach the maps from every
// thread. One that never left the host is not.
BOOST_AUTO_TEST_CASE( theSetupPredicateRefusesWhatReachesTheMemoryManager )
{
	deviceOnce();
	BOOST_TEST_MESSAGE( "device: " << ( hostMemoryIsPlain() ? "plain HOST" : "registering host" )
	                    << " memory, threading " << ( threadingBuilt() ? "built" : "not built" ) );

	mfem::Vector hostOnly( 16 );
	hostOnly = 1.0;

	mfem::Vector onDevice( 16 );
	onDevice.UseDevice( true );
	onDevice = 2.0;
	onDevice.Write();

	bool const plain = threadingBuilt() && hostMemoryIsPlain();
	BOOST_TEST( meq::setupLoopMayThread( {} ) == plain );
	BOOST_TEST( meq::setupLoopMayThread( { nullptr } ) == plain );
	BOOST_TEST( meq::setupLoopMayThread( { hostOnly.GetData() } ) == plain );
	BOOST_TEST( !meq::setupLoopMayThread( { hostOnly.GetData(), onDevice.GetData() } ) );
}

// A THREADED WARM START UNDER THE DEVICE IS THE ONE-THREAD WARM START, TO THE
// BIT.
//
// Each element's projection and each point's search is the same arithmetic in
// the same order whichever thread runs it, so equality is the claim, not a
// tolerance. On the cuda arm this is the threaded loops under a Device; on the
// debug arm both runs are one thread, and the case checks the predicate said
// so.
BOOST_AUTO_TEST_CASE( aThreadedWarmStartIsTheSerialOneBitForBit )
{
	deviceOnce();

	mfem::Mesh source = guessMesh();
	mfem::L2_FECollection sourceCollection( 2, 2, mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace sourceSpace( &source, &sourceCollection );
	mfem::GridFunction guess( &sourceSpace );
	mfem::FunctionCoefficient guessCoefficient( guessValue );
	guess.ProjectCoefficient( guessCoefficient );

	mfem::Mesh target = mfem::Mesh::MakeCartesian2D(
		11, 7, mfem::Element::TRIANGLE, false, 1.7, 0.9 );
	mfem::L2_FECollection targetCollection( 3, 2, mfem::BasisType::GaussLobatto );
	mfem::FiniteElementSpace targetSpace( &target, &targetCollection );

	meq::FieldTransfer transfer( source );
	mfem::ConstantCoefficient fallback( 0.0 );

	bool const mayThread = meq::setupLoopMayThread(
		{ source.GetNodes()->GetData() } );
	BOOST_TEST( mayThread == ( threadingBuilt() && hostMemoryIsPlain() ) );
	BOOST_TEST_MESSAGE( "the warm start's loops " << ( mayThread ? "thread" : "run on one thread" ) );

#ifdef MFEM_USE_OPENMP
	int const threadsBefore = omp_get_max_threads();
	omp_set_num_threads( 1 );
#endif
	mfem::GridFunction serial( &targetSpace );
	serial = 0.0;
	BOOST_TEST( transfer.transfer( guess, fallback, serial ) == 0 );

#ifdef MFEM_USE_OPENMP
	omp_set_num_threads( 4 );
#endif
	mfem::GridFunction threaded( &targetSpace );
	threaded = 0.0;
	BOOST_TEST( transfer.transfer( guess, fallback, threaded ) == 0 );
#ifdef MFEM_USE_OPENMP
	omp_set_num_threads( threadsBefore );
#endif

	serial.HostRead();
	threaded.HostRead();
	BOOST_TEST( std::memcmp( serial.GetData(), threaded.GetData(),
	                         sizeof( mfem::real_t )*static_cast<std::size_t>( serial.Size() ) ) == 0 );

	// And it is a transfer, not two copies of the same wrong answer. Measured
	// 4.53e-4 on both arms, mostly the degree-2 guess's own error; the bound is
	// a decade above it.
	double const error = threaded.ComputeL2Error( guessCoefficient );
	BOOST_TEST_MESSAGE( "L2 error of the transferred guess: " << error );
	BOOST_TEST( error < 5.0e-3 );
}
