#include "DeviceRows.hpp"

#include "Crossings.hpp"

#include <sstream>
#include <stdexcept>

namespace meq
{

	void DeviceRows::upload( std::vector<int> const &support,
	                         std::vector<double> const &valuesIn,
	                         int rowCountIn )
	{
		if ( rowCountIn < 0 )
			throw std::invalid_argument(
				"meq::DeviceRows::upload: the row count must not be negative" );

		int const w = static_cast<int>( support.size() );
		std::size_t const wanted =
			static_cast<std::size_t>( rowCountIn )*support.size();

		/*
		 * CHECKED RATHER THAN TRUSTED, because the failure it guards is a read
		 * off the end of a device buffer -- which under CUDA is not a fault,
		 * it is whatever was in that memory, arriving as a plausible border
		 * row. meq::CompressedRows makes the same argument for
		 * firstNonzeroOutside(): a row wrong in magnitude does not diverge,
		 * it reaches a different equilibrium.
		 */
		if ( valuesIn.size() != wanted )
		{
			std::ostringstream message;
			message << "meq::DeviceRows::upload: " << valuesIn.size()
			        << " values for " << rowCountIn << " rows of width " << w
			        << ", which wants " << wanted
			        << ". The values are row major and share one support, "
			           "exactly as meq::CompressedRows stores them.";
			throw std::invalid_argument( message.str() );
		}

		if ( rowCountIn == 0 || w == 0 )
		{
			clear();
			return;
		}

		rows = rowCountIn;
		columns = w;

		/*
		 * `HostWrite()` BEFORE EVERY FILL, AND THE FIRST CALL IS THE ONLY ONE
		 * THAT WOULD SURVIVE WITHOUT IT.
		 *
		 * The mirror is refreshed once per Newton step, and the `Read()` at
		 * the end of this function leaves both members DEVICE-valid. On the
		 * next call `SetSize()` at an unchanged size reallocates nothing, so
		 * the flags persist -- and `indices[ k ] = ...` is then a raw host
		 * write into a page the memory manager has protected. Under
		 * `Device( "debug" )` that is a named fault; under CUDA it writes the
		 * host copy and the kernel goes on reading the STALE device one, which
		 * is a wrong border row arriving as a plausible number. This file's
		 * own comment above makes that argument about reading off the end of
		 * the buffer; it applies here for a different reason.
		 *
		 * Found by `tests/convergence/DeviceParity.cpp` on its first run,
		 * which is what that test exists for: no host build can see this, and
		 * no CUDA run has yet, because `columnIsResident()` keeps the mirror
		 * unread today. MEASUREMENTS.md M-211.
		 */
		indices.SetSize( w );
		indices.HostWrite();
		for ( int k = 0; k < w; ++k )
			indices[ k ] = support[ static_cast<std::size_t>( k ) ];

		// Taken HERE, while the array is still host memory. See the member's
		// comment: after the Read() below it may not be looked at on the host.
		largestIndex = indices[ w - 1 ];

		values.SetSize( static_cast<int>( wanted ) );
		values.HostWrite();
		for ( std::size_t k = 0; k < wanted; ++k )
			values( static_cast<int>( k ) ) = valuesIn[ k ];

		/*
		 * THE MIRROR IS DEVICE MEMORY FROM HERE ON, AND SAYING SO IS WHAT MAKES
		 * IT SO. A freshly constructed mfem::Vector has UseDevice() FALSE --
		 * which is the exact property that let M-187's defect survive mode 0
		 * and fail every mode after it -- so a mirror that is never marked
		 * would be written on the host, read on the host, and quietly never
		 * reach the card at all. The indices are an mfem::Array<int> for the
		 * same reason: it is device state on the same footing as a Vector, as
		 * meq::CompressedRows' own comment about rowDot() records.
		 */
		values.UseDevice( true );
		values.Read();
		indices.Read();
	}

	void DeviceRows::clear()
	{
		indices.SetSize( 0 );
		values.SetSize( 0 );
		rows = 0;
		columns = 0;
		largestIndex = -1;
	}

	RowView DeviceRows::row( int r ) const
	{
		RowView view{};
		view.indices = indices.Read();
		view.weights = values.Read() + static_cast<std::size_t>( r )*columns;
		view.count = columns;
		view.negate = 0;
		return view;
	}

	bool dotRowsAgainstColumnsOnDevice( mfem::Array<RowView> const &rows,
	                                    mfem::Array<mfem::real_t const *> const &columns,
	                                    mfem::Vector &out )
	{
		if ( !mfem::Device::Allows( mfem::Backend::DEVICE_MASK ) )
			return false;
		if ( rows.Size() == 0 || columns.Size() == 0 )
			return false;
		out.SetSize( rows.Size()*columns.Size() );
		out.UseDevice( true );
		// The two descriptor tables are rebuilt per fill and small -- a few
		// dozen pointers -- and counted, so they show as the uploads they are.
		countedDeviceRead( rows, __LINE__ );
		countedDeviceRead( columns, __LINE__ );
		return dotRowsAgainstColumnsKernel( rows, columns, out );
	}

	bool DeviceRows::valid() const
	{
		return rows > 0 && columns > 0;
	}

	int DeviceRows::rowCount() const
	{
		return rows;
	}

	int DeviceRows::width() const
	{
		return columns;
	}

	bool DeviceRows::dotAll( mfem::Vector const &v, mfem::Vector &out ) const
	{
		if ( !valid() )
			return false;

		/*
		 * THE ONE PRECONDITION WORTH CHECKING HERE RATHER THAN IN THE KERNEL.
		 * A support index past the end of v is the same off-the-end read the
		 * upload guards, arriving by a different route: rows compressed
		 * against one space and contracted against a vector from another. It
		 * is cheap to check because the support is sorted, so the largest
		 * index is the last one, and it needs no host read of v -- Size() is a
		 * property of the container and not of its data.
		 */
		if ( largestIndex >= v.Size() )
		{
			std::ostringstream message;
			message << "meq::DeviceRows::dotAll: the support reaches index "
			        << largestIndex << " and the vector has "
			        << v.Size() << " entries. These rows were compressed "
			           "against a different space.";
			throw std::invalid_argument( message.str() );
		}

		out.SetSize( rows );
		out.UseDevice( true );

		return dotCompressedRowsOnDevice( values, indices, rows, columns,
		                                  v, out );
	}

	bool DeviceRows::dotColumns( mfem::Vector const &columnData,
	                             int columnCount, int length,
	                             mfem::DenseMatrix &out ) const
	{
		if ( !valid() || columnCount <= 0 || length <= 0 )
			return false;

		// The same bounds check dotAll() makes, against the COLUMN length --
		// and on `largestIndex` rather than on the array, which after upload()
		// is device memory. See its member comment.
		if ( largestIndex >= length )
			return false;
		if ( columnData.Size() != columnCount*length )
			return false;

		/*
		 * THE KERNEL WRITES THROUGH THE DenseMatrix's OWN MEMORY, NOT THROUGH A
		 * Vector WRAPPED ROUND ITS POINTER.
		 *
		 * An earlier draft did `mfem::Vector view( out.GetData(), n )` and
		 * handed that across. That is an ALIAS with its own flags: the kernel
		 * left the view device-valid, the view was destroyed, and the caller's
		 * matrix carried a host buffer nothing had synced -- so reading it
		 * faulted under Device( "debug" ) with "invalid permissions". **That is
		 * M-187's defect class, committed in the file whose header argues about
		 * it**, and mfem::DenseMatrix has had Read()/Write()/HostRead() of its
		 * own the whole time.
		 *
		 * DenseMatrix stores column major, so entry ( r, c ) is at
		 * r + rowCount*c -- the layout the launch writes, stated there rather
		 * than assumed, because a transposed bordered matrix factors and solves
		 * perfectly well to the wrong step.
		 */
		out.SetSize( rows, columnCount );

		return dotCompressedColumnsOnDevice( values, indices, rows, columns,
		                                     columnData, columnCount, length,
		                                     out );
	}

}
