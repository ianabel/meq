#ifndef MEQ_TESTS_DEVICEKERNELS_HPP
#define MEQ_TESTS_DEVICEKERNELS_HPP

/*
 * A CASE THAT ASSERTS A KERNEL RAN, ON A BUILD THAT HAS NO KERNELS.
 *
 * MEQ's seven `*Device.cpp` files are host stubs on a build against an MFEM
 * without CUDA: each route refuses before it reads anything, and every caller
 * keeps its host loop. `mfem::Device( "debug" )` is still configurable there
 * and still allows `DEVICE_MASK`, so a case that requires a kernel to have run
 * under it reads the stub's refusal as a red -- about code that was never
 * compiled. Such a case is a precondition-skip here, and says so on stdout, the
 * way `DeviceParity`'s cuda arm says "SKIPPED -- this build has no CUDA".
 *
 * Only the cases that REQUIRE a kernel take it. A refusal case, a counter case
 * or a host-only case in the same file runs on every build.
 *
 * Usage:
 *
 *     BOOST_AUTO_TEST_CASE( theKernelDoesX,
 *                           *meq::tests::requiresDeviceKernels() )
 */

#include "meq/SourceView.hpp"

#include <boost/test/unit_test.hpp>

#include <cstdio>

namespace meq::tests
{
	inline boost::unit_test::precondition requiresDeviceKernels()
	{
		return boost::unit_test::precondition(
			[]( boost::unit_test::test_unit_id id )
			{
				if ( meq::deviceKernelsBuilt() )
					return boost::test_tools::assertion_result( true );
				std::printf( "\n  %s: SKIPPED -- this build has no device "
				             "kernels (MEQ's *Device.cpp files are host stubs "
				             "without CUDA), so the kernel this case asserts "
				             "cannot run. A CUDA build runs it under "
				             "mfem::Device( \"debug\" ) and MEQ_TEST_DEVICE=cuda.\n",
				             boost::unit_test::framework::get( id,
				                 boost::unit_test::TUT_ANY ).p_name.get().c_str() );
				boost::test_tools::assertion_result skipped( false );
				skipped.message() << "this build has no device kernels";
				return skipped;
			} );
	}
}

#endif
