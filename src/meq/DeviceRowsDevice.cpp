#include "DeviceRows.hpp"

/*
 * THE DEVICE TRANSLATION UNIT FOR meq::DeviceRows.
 *
 * One file per feature is compiled by the device compiler, because
 * `mfem::forall` static_asserts under g++ -- PointLocatorDevice.cpp carries the
 * full account and SourceViewDevice.cpp is the other one. CMakeLists.txt marks
 * these and only these with LANGUAGE CUDA, and only where a device language was
 * actually enabled; on a host build this compiles to the stub below and the
 * `naming` check reads it like any other file.
 */

namespace meq
{

#if defined( MFEM_USE_CUDA_OR_HIP_LANG )

	bool dotCompressedRowsOnDevice( mfem::Vector const &values,
	                                mfem::Array<int> const &indices,
	                                int rowCount, int width,
	                                mfem::Vector const &v,
	                                mfem::Vector &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( rowCount <= 0 || width <= 0 )
			return false;

		auto const *d_values  = values.Read();
		auto const *d_indices = indices.Read();
		auto const *d_v       = v.Read();
		auto *d_out           = out.Write();

		/*
		 * ONE THREAD PER ROW, AND THE SUM INSIDE IT IS SEQUENTIAL AND
		 * ASCENDING. That is the whole design and the header says why: it
		 * reproduces meq::CompressedRows::dot()'s every partial sum, so the
		 * answer is bit-identical rather than merely close, and MEQ keeps the
		 * digit-for-digit comparison between --device cuda and --device cpu
		 * that root-caused M-79, M-130 and M-187.
		 *
		 * A tree reduction within a row would reassociate the sum. It would
		 * also be the wrong trade: this launch exists to avoid bringing 698 KB
		 * of state to the host, not to make a few thousand multiply-adds fast,
		 * and the rows -- N + 4 of them, independent -- are the parallelism
		 * that matters.
		 */
		mfem::forall( rowCount, [ = ] MFEM_HOST_DEVICE ( int r )
		{
			mfem::real_t total = mfem::real_t( 0 );
			int const base = r*width;
			for ( int k = 0; k < width; ++k )
				total += d_values[ base + k ]*d_v[ d_indices[ k ] ];
			d_out[ r ] = total;
		} );

		/*
		 * LEFT DEVICE-VALID, AND THE CALLER FUNNELS. An earlier draft did
		 * `out.HostRead()` here, on the reasoning that the caller wants N + 4
		 * scalars and the bordered step's arithmetic around them is host
		 * arithmetic today. That is true today and is the WRONG place to
		 * decide it: a kernel that forces the read makes it impossible for a
		 * resident border ever to keep the result, which is the whole of
		 * DEVICE-PLAN.md's phase C. Whether the host needs these scalars is
		 * the CALLER's question, and MEQ's convention everywhere else (§1.2) is
		 * that the reader funnels.
		 */
		return true;
	}

	bool dotCompressedColumnsOnDevice( mfem::Vector const &values,
	                                   mfem::Array<int> const &indices,
	                                   int rowCount, int width,
	                                   mfem::Vector const &columns,
	                                   int columnCount, int length,
	                                   mfem::DenseMatrix &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( rowCount <= 0 || width <= 0 || columnCount <= 0 || length <= 0 )
			return false;

		auto const *d_values  = values.Read();
		auto const *d_indices = indices.Read();
		auto const *d_cols    = columns.Read();
		auto *d_out           = out.Write();

		int const nr = rowCount;
		int const nw = width;
		int const nl = length;

		/*
		 * ONE THREAD PER ( ROW, COLUMN ) PAIR, AND THE SUM INSIDE IT IS
		 * SEQUENTIAL AND ASCENDING -- the same contract dotCompressedRowsOnDevice()
		 * keeps and for the same reason. There are ( N + 4 )^2 of them, 196 on
		 * the shipped machine case, which is more parallelism than the
		 * per-row launch has and still not much; the launch exists to keep the
		 * COLUMNS on the device, not to be fast.
		 *
		 * THE OUTPUT IS COLUMN MAJOR, `r + nr*c`, because an mfem::DenseMatrix
		 * is. Stated rather than assumed: a transposed bordered matrix factors
		 * and solves perfectly well, to the wrong step.
		 */
		mfem::forall( rowCount*columnCount, [ = ] MFEM_HOST_DEVICE ( int t )
		{
			int const r = t % nr;
			int const c = t / nr;

			mfem::real_t total = mfem::real_t( 0 );
			int const base = r*nw;
			int const cbase = c*nl;
			for ( int k = 0; k < nw; ++k )
				total += d_values[ base + k ]*d_cols[ cbase + d_indices[ k ] ];
			d_out[ r + nr*c ] = total;
		} );

		return true;
	}

	bool dotRowsAgainstColumnsKernel( mfem::Array<RowView> const &rows,
	                                  mfem::Array<mfem::real_t const *> const &columns,
	                                  mfem::Vector &out )
	{
		auto const *d_rows = rows.Read();
		auto const *d_columns = columns.Read();
		auto *d_out = out.Write();
		int const nr = rows.Size();

		// One thread per ( row, column ), summing from zero in k order: the
		// host rowDot()'s loop exactly, which is what -fmad=false is for.
		mfem::forall( rows.Size()*columns.Size(), [ = ] MFEM_HOST_DEVICE ( int t )
		{
			int const r = t % nr;
			int const c = t / nr;
			RowView const row = d_rows[ r ];
			mfem::real_t const *const column = d_columns[ c ];
			mfem::real_t total = mfem::real_t( 0 );
			for ( int k = 0; k < row.count; ++k )
				total += row.weights[ k ]*column[ row.indices[ k ] ];
			d_out[ t ] = row.negate ? -total : total;
		} );
		return true;
	}

#else

	bool dotRowsAgainstColumnsKernel( mfem::Array<RowView> const &,
	                                  mfem::Array<mfem::real_t const *> const &,
	                                  mfem::Vector & )
	{
		return false;
	}

	bool dotCompressedRowsOnDevice( mfem::Vector const &, mfem::Array<int> const &,
	                                int, int, mfem::Vector const &,
	                                mfem::Vector & )
	{
		return false;
	}

	bool dotCompressedColumnsOnDevice( mfem::Vector const &,
	                                   mfem::Array<int> const &, int, int,
	                                   mfem::Vector const &, int, int,
	                                   mfem::DenseMatrix & )
	{
		return false;
	}

#endif

}
