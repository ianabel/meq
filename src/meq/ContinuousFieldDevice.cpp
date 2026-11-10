#include "ContinuousField.hpp"

/*
 * THE LAUNCH FOR meq::ContinuousInterpolant, AND NOTHING ELSE. The blocks, the
 * host route this reproduces and the account of both are in
 * ContinuousField.{hpp,cpp}, which stay under the naming check; this file is
 * handed to the device compiler and cannot be.
 *
 * COMPILED WITHOUT FMA CONTRACTION, for meq::DeviceRows' recorded reason: the
 * kernel is the host route's loop, operation for operation, so that the
 * `debug` backend reproduces the host to the bit and the card differs from it
 * only where its own rounding does. CMakeLists.txt carries the flag.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	bool applyLinearBlocksOnDevice( LinearBlocks const &blocks,
	                                mfem::Vector const &first,
	                                mfem::Vector const *second, int split,
	                                mfem::Vector &out, int outShift )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		int const count = blocks.count;
		if ( count == 0 )
			return true;

		auto const *d_outOffset = blocks.outOffset.Read();
		auto const *d_inOffset = blocks.inOffset.Read();
		auto const *d_knownOffset = blocks.knownOffset.Read();
		auto const *d_matrixOffset = blocks.matrixOffset.Read();
		auto const *d_outIndex = blocks.outIndex.Read();
		auto const *d_inIndex = blocks.inIndex.Read();
		auto const *d_knownIndex = blocks.knownIndex.Read();
		auto const *d_matrix = blocks.matrix.Read();
		auto const *d_first = first.Read();
		auto const *d_second = second != nullptr ? second->Read() : d_first;
		// ReadWrite and not Write: a stage reads what an earlier stage wrote.
		auto *d_out = out.ReadWrite();
		int const cut = split;
		int const shift = outShift;

		mfem::forall( count, [ = ] MFEM_HOST_DEVICE ( int b )
		{
			int const outBegin = d_outOffset[ b ];
			int const rows = d_outOffset[ b + 1 ] - outBegin;
			int const inBegin = d_inOffset[ b ];
			int const inputs = d_inOffset[ b + 1 ] - inBegin;
			int const knownBegin = d_knownOffset[ b ];
			int const knowns = d_knownOffset[ b + 1 ] - knownBegin;
			int const columns = inputs + knowns;
			mfem::real_t const *m = d_matrix + d_matrixOffset[ b ];
			for ( int r = 0; r < rows; ++r )
			{
				mfem::real_t sum = mfem::real_t( 0 );
				for ( int c = 0; c < inputs; ++c )
				{
					int const at = d_inIndex[ inBegin + c ];
					mfem::real_t const value = at < cut ? d_first[ at ] : d_second[ at - cut ];
					sum += m[ r*columns + c ]*value;
				}
				for ( int c = 0; c < knowns; ++c )
					sum += m[ r*columns + inputs + c ]*d_out[ d_knownIndex[ knownBegin + c ] + shift ];
				d_out[ d_outIndex[ outBegin + r ] + shift ] = sum;
			}
		} );
		return true;
	}

#else

	bool applyLinearBlocksOnDevice( LinearBlocks const &, mfem::Vector const &,
	                                mfem::Vector const *, int, mfem::Vector &, int )
	{
		return false;
	}

#endif

}
