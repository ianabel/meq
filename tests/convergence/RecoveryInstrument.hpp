#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>

#include "mfem.hpp"

#include "meq/GradShafranov.hpp"

/*
 * THE TWO PRECISION INSTRUMENTS' SETTINGS, FROM THE ENVIRONMENT --
 * RecoveryPrecision and ExtractionNoise, both built and not registered, so no
 * ctest ever reads these:
 *
 *   MEQ_TEST_DEVICE        an mfem::Device configured before any case runs
 *                          ("cuda", "debug"); unset, none. Under a device
 *                          the solver takes the driver's device routes:
 *                          AssemblyMode, LocalFactorMode and
 *                          TraceAssemblyMode batched, cuDSS under CUDA.
 *   MEQ_TEST_BATCHED       "1": the same batched keys with no device, so
 *                          the batched local routes can be read on the host.
 *   MEQ_LOCAL_REFINEMENT   "off" or "on": GradShafranovSolver::
 *                          setLocalRefinement(); unset, the solver's default
 *                          (on).
 *
 * Each case prints which routes were taken, since a batched key that
 * silently fell back would put the host's numbers under the device's name.
 */
namespace meq::tests::recovery
{
	inline char const *environment( char const *name )
	{
		char const *const value = std::getenv( name );
		return value != nullptr && *value != '\0' ? value : nullptr;
	}

	/// A function-local static, as in PlasmaConnectivity.cpp: a namespace-scope
	/// mfem::Device may be constructed before MFEM's own singleton.
	struct DeviceFromEnvironment
	{
		DeviceFromEnvironment()
		{
			if ( char const *const requested = environment( "MEQ_TEST_DEVICE" ) )
			{
				static mfem::Device device( requested );
				device.Print();
			}
		}
	};

	inline bool onDevice()
	{
		return mfem::Device::Allows( mfem::Backend::DEVICE_MASK );
	}

	inline bool batched()
	{
		char const *const forced = environment( "MEQ_TEST_BATCHED" );
		return onDevice() || ( forced != nullptr && std::string( forced ) == "1" );
	}

	inline void configure( meq::GradShafranovSolver &solver )
	{
		if ( char const *const refinement = environment( "MEQ_LOCAL_REFINEMENT" ) )
		{
			std::string const r( refinement );
			if ( r != "on" && r != "off" )
				throw std::invalid_argument( "MEQ_LOCAL_REFINEMENT is \"on\" or \"off\"" );
			solver.setLocalRefinement( r == "on" );
		}
		if ( batched() )
		{
			solver.setAssemblyMode( meq::GradShafranovSolver::AssemblyMode::Batched );
			solver.setLocalFactorMode( meq::GradShafranovSolver::LocalFactorMode::Batched );
			solver.setTraceAssemblyMode( meq::GradShafranovSolver::TraceAssemblyMode::Batched );
		}
#ifdef MFEM_USE_CUDSS
		if ( mfem::Device::Allows( mfem::Backend::CUDA_MASK ) )
			solver.setTraceSolver( meq::GradShafranovSolver::TraceSolver::cuDSS );
#endif
	}

	/// After a solve and its postProcess().
	inline void reportRoutes( meq::GradShafranovSolver const &solver )
	{
		// batchedLocalFactorTaken() is upstream's CanBatchLocalFactor(), which
		// asks the block sizes and not the mode, so the mode is printed beside
		// it; CanBatchLocalSolve() asks both.
		std::printf( "    routes: %s, local refinement %s; local factor mode %s (uniform blocks %s), local "
		             "solve batched %s, reconstruction batched %s\n",
		             onDevice() ? "device" : "host", solver.localRefinement() ? "ON" : "off",
		             solver.localFactorMode() == meq::GradShafranovSolver::LocalFactorMode::Batched
		                 ? "batched" : "serial",
		             solver.batchedLocalFactorTaken() ? "yes" : "no",
		             solver.batchedLocalSolveTaken() ? "yes" : "no",
		             solver.batchedReconstructionTaken() ? "yes" : "no" );
	}
}
