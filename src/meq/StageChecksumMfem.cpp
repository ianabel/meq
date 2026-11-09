#include "meq/StageChecksumMfem.hpp"

#include <cmath>
#include <cstddef>

namespace meq
{

	void stageChecksum( char const *stage, mfem::Vector const &v )
	{
		if ( !stageChecksumWanted() ) { return; }

		StageReading r;
		r.stage = stage;
		r.size = v.Size();

		// BEFORE THE READ. This is the state a raw accessor here would meet,
		// and it is the half of the diagnostic that survives the instrument
		// being a sync -- see meq/StageChecksum.hpp.
		mfem::Memory<mfem::real_t> const &memory = v.GetMemory();
		r.hostValid = memory.HostIsValid();
		r.deviceValid = memory.DeviceIsValid();

		mfem::real_t const *const p = v.HostRead();
		r.hash = stageHash( p, static_cast< std::size_t >( v.Size() )
		                       *sizeof( mfem::real_t ) );
		double sum = 0.0;
		for ( int i = 0; i < v.Size(); ++i )
		{
			double const x = static_cast< double >( p[ i ] );
			sum += x*x;
		}
		r.norm = std::sqrt( sum );
		recordStageReading( r );
	}

	void stageChecksum( char const *stage, mfem::Array<int> const &a )
	{
		if ( !stageChecksumWanted() ) { return; }

		StageReading r;
		r.stage = stage;
		r.size = a.Size();

		mfem::Memory<int> const &memory = a.GetMemory();
		r.hostValid = memory.HostIsValid();
		r.deviceValid = memory.DeviceIsValid();

		int const *const p = a.HostRead();
		r.hash = stageHash( p, static_cast< std::size_t >( a.Size() )
		                       *sizeof( int ) );
		double sum = 0.0;
		for ( int i = 0; i < a.Size(); ++i )
		{
			double const x = static_cast< double >( p[ i ] );
			sum += x*x;
		}
		r.norm = std::sqrt( sum );
		recordStageReading( r );
	}

}
