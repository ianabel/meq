#include "meq/StageChecksum.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace meq
{

	namespace
	{

		std::vector<StageReading> &readings()
		{
			static std::vector<StageReading> theReadings;
			return theReadings;
		}

	}

	bool stageChecksumWanted()
	{
		static bool const wanted = std::getenv( "MEQ_STAGE_CHECKSUM" ) != nullptr;
		return wanted;
	}

	std::uint64_t stageHash( void const *data, std::size_t bytes )
	{
		auto const *p = static_cast< unsigned char const * >( data );
		std::uint64_t h = 1469598103934665603ULL;
		for ( std::size_t i = 0; i < bytes; ++i )
		{
			h ^= p[ i ];
			h *= 1099511628211ULL;
		}
		return h;
	}

	void recordStageReading( StageReading const &r )
	{
		readings().push_back( r );
		std::printf( "MEQ-STAGE %-38s n=%-8d %s%s hash=%016llx "
		             "norm=%.17e\n",
		             r.stage.c_str(), r.size,
		             r.hostValid ? "H+" : "H-",
		             r.deviceValid ? "D+" : "D-",
		             static_cast< unsigned long long >( r.hash ),
		             r.norm );
		// FLUSHED PER LINE ON PURPOSE. A walk's whole value is the FIRST
		// stage that disagrees, and a run that aborts -- which is what
		// `Device( "debug" )` does when it finds anything -- must not lose
		// the lines leading up to it to a buffer.
		std::fflush( stdout );
	}

	void stageChecksum( char const *stage, std::vector<double> const &v )
	{
		if ( !stageChecksumWanted() ) { return; }

		StageReading r;
		r.stage = stage;
		r.size = static_cast< int >( v.size() );
		// A std::vector has no memory manager behind it, so the host copy is
		// the only copy and saying so keeps the columns aligned with the
		// mfem::Vector lines of StageChecksumMfem.cpp rather than leaving a
		// reader to wonder.
		r.hostValid = true;
		r.deviceValid = false;
		r.hash = v.empty() ? 0
		                   : stageHash( v.data(), v.size()*sizeof( double ) );
		double sum = 0.0;
		for ( double const x : v ) { sum += x*x; }
		r.norm = std::sqrt( sum );
		recordStageReading( r );
	}

	std::vector<StageReading> const &stageReadings()
	{
		return readings();
	}

	void clearStageReadings()
	{
		readings().clear();
	}

}
