#include "AxisSearch.hpp"

#include "Crossings.hpp"

#include <cstdlib>
#include <cstring>

/*
 * THE HOST HALF OF meq::AxisSearch: THE REFUSALS, THE MARSHALLING AND THE ONE
 * DOWNLOAD. THE LAUNCHES ARE IN AxisSearchDevice.cpp.
 *
 * The split follows meq::DeviceQuadrature's, for the reason that file's own
 * header gives: a translation unit holding an `mfem::forall` must be compiled
 * by nvcc, nvcc's flags do not parse under clang, and CMake therefore excludes
 * every such file from the `naming` ctest. Everything a reader looks up --
 * every refusal, every bound, every name -- stays here, where the convention
 * check can see it, and only the kernels sit outside.
 */

namespace meq
{

	bool axisDeviceRouteForced()
	{
		// A function-local static: read once, on first use, after all static
		// initialisation has finished. The same arrangement the device
		// singletons in the unit tests use and for the same ordering reason.
		static bool const forced = []
		{
			char const *const value = std::getenv( "MEQ_AXIS_DEVICE" );
			return value != nullptr && std::strcmp( value, "force" ) == 0;
		}();
		return forced;
	}

	bool axisFluxScaleOnDevice( mfem::Vector const &state, int scalarDofs,
	                            mfem::real_t tolerance, mfem::Vector &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK )
		     || scalarDofs <= 0 || state.Size() < 2*scalarDofs )
			return false;

		return axisFluxScaleKernel( state, scalarDofs, tolerance, out );
	}

	bool axisCandidatesOnDevice( AxisSearchInputs const &inputs,
	                             mfem::Array<AxisCandidate> &out, int &count )
	{
		count = 0;
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK )
		     || inputs.elements <= 0 || inputs.flux.change == nullptr
		     || inputs.potential.change == nullptr || inputs.maps == nullptr
		     || inputs.state == nullptr || inputs.target == nullptr )
			return false;

		/*
		 * TWO CANDIDATES AN ELEMENT, BECAUSE `elementSeeds()` TRIES TWO SEEDS.
		 * The host roots from the element centre and from the node carrying
		 * the smallest `| q |`, keeps whichever are accepted, and merges them
		 * afterwards. Keeping both here is what makes the candidate list the
		 * host's, and the merge is the host's to skip -- see the header on why
		 * the selection does not need it.
		 *
		 * **FRESH ARRAYS AND NOT MEMBERS, WHICH IS M-226 APPLIED RATHER THAN
		 * RE-LEARNED.** A persistent buffer refilled in place after a
		 * `SetSize()` at an unchanged size reallocates nothing, so the fill
		 * writes into a page the manager still believes the device owns. A
		 * freshly constructed `mfem::Array` has `UseDevice()` false and cannot
		 * be in that state; these are per call, and the only persistent one is
		 * `out`, which the caller owns and which is written by `Write()`
		 * rather than by a host fill.
		 */
		mfem::Array<AxisCandidate> perElement( 2*inputs.elements );
		mfem::Array<int> offsets;
		mfem::Array<int> total( 1 );

		if ( !axisCandidatesCountKernel( inputs, perElement, offsets, total ) )
			return false;

		/*
		 * THE COUNT COMES BACK FIRST, AND THE COMPACTED ARRAY IS SIZED TO IT.
		 *
		 * **FOUR BYTES, SO THAT THE SECOND DOWNLOAD IS PROPORTIONAL TO THE
		 * CANDIDATES THAT EXIST RATHER THAN TO A CEILING.** The first draft
		 * downloaded a fixed `axisCandidateCap` records, which made the
		 * ceiling the cost and forced it down to 256 -- and on
		 * `fixed-f-diiid`'s 102,840 elements a transient iterate exceeded
		 * that, the route declined, and the host sweep that replaced it cost
		 * **22.4 seconds**, the whole leg. *A refusal that costs more than
		 * the thing it refuses is not a safe default.*
		 *
		 * Two downloads where there was one, and the pair is smaller than the
		 * single one was: a handful of records at convergence against 22 KB
		 * of mostly-empty array, and against the 830 KB state either replaces.
		 */
		countedHostRead( total, __LINE__ );
		count = total[ 0 ];

		// A TRUNCATED competition is one missing an entrant, which is exactly
		// how a conductor's O-point would win one. The count is reported
		// truthfully either way so a caller can say what happened.
		if ( count > axisCandidateCap )
			return false;

		out.SetSize( count > 0 ? count : 1 );
		axisCandidatesScatterKernel( inputs, perElement, offsets, out );
		countedHostRead( out, __LINE__ );
		return true;
	}

	bool xPointFieldOnDevice( XPointFieldInputs const &inputs,
	                          XPointFieldValues &out )
	{
		out = XPointFieldValues{};
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK )
		     || inputs.element < 0 || inputs.flux.change == nullptr
		     || inputs.potential.change == nullptr || inputs.maps == nullptr
		     || inputs.state == nullptr )
			return false;

		/*
		 * A FRESH ARRAY OF ONE, NOT A MEMBER. M-226 is a persistent buffer
		 * refilled in place after a `SetSize()` at an unchanged size, which
		 * reallocates nothing and leaves the fill writing into a page the
		 * manager still believes the device owns. A freshly constructed
		 * `mfem::Array` has `UseDevice()` false and cannot be in that state,
		 * and one record a call is not worth a lifetime to get wrong.
		 */
		mfem::Array<XPointFieldValues> result( 1 );
		if ( !xPointFieldKernel( inputs, result ) )
			return false;

		// The one download: eight reals and a flag, against the 830 KB the
		// host route brings over to compute the same eight.
		countedHostRead( result, __LINE__ );
		out = result[ 0 ];

		// The element was not affine. The caller keeps the host path, which
		// handles a curved element through mfem::ElementTransformation.
		return out.ok != 0;
	}

}
