#include "meq/StageChecksum.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace meq
{

	namespace
	{

		/**
		 * FNV-1a OVER THE RAW BYTES, AND THE CHOICE IS ABOUT WHAT IT MUST
		 * CATCH.
		 *
		 * The question a walk asks is *did this buffer change at all*, so the
		 * hash has to see the last bit of every mantissa -- a checksum built
		 * from sums or norms cancels, and two vectors differing by a sign flip
		 * in two entries share a norm exactly. Bytes of the IEEE
		 * representation have no such coincidence, and FNV-1a is four lines.
		 *
		 * It is NOT a claim about collision resistance and nothing here needs
		 * one: the comparison is between two runs of the same code on the same
		 * problem, where the alternative to "the hashes match" is a difference
		 * large enough that the L2 beside it shows it.
		 */
		std::uint64_t fnv1a( void const *data, std::size_t bytes )
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

		std::vector<StageReading> &readings()
		{
			static std::vector<StageReading> theReadings;
			return theReadings;
		}

		void report( StageReading const &r )
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

	}

	bool stageChecksumWanted()
	{
		static bool const wanted = std::getenv( "MEQ_STAGE_CHECKSUM" ) != nullptr;
		return wanted;
	}

	void stageChecksum( char const *stage, mfem::Vector const &v )
	{
		if ( !stageChecksumWanted() ) { return; }

		StageReading r;
		r.stage = stage;
		r.size = v.Size();

		// BEFORE THE READ. This is the state a raw accessor here would meet,
		// and it is the half of the diagnostic that survives the instrument
		// being a sync -- see the header.
		mfem::Memory<mfem::real_t> const &memory = v.GetMemory();
		r.hostValid = memory.HostIsValid();
		r.deviceValid = memory.DeviceIsValid();

		mfem::real_t const *const p = v.HostRead();
		r.hash = fnv1a( p, static_cast< std::size_t >( v.Size() )
		                   *sizeof( mfem::real_t ) );
		double sum = 0.0;
		for ( int i = 0; i < v.Size(); ++i )
		{
			double const x = static_cast< double >( p[ i ] );
			sum += x*x;
		}
		r.norm = std::sqrt( sum );
		report( r );
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
		r.hash = fnv1a( p, static_cast< std::size_t >( a.Size() )
		                   *sizeof( int ) );
		double sum = 0.0;
		for ( int i = 0; i < a.Size(); ++i )
		{
			double const x = static_cast< double >( p[ i ] );
			sum += x*x;
		}
		r.norm = std::sqrt( sum );
		report( r );
	}

	void stageChecksum( char const *stage, std::vector<double> const &v )
	{
		if ( !stageChecksumWanted() ) { return; }

		StageReading r;
		r.stage = stage;
		r.size = static_cast< int >( v.size() );
		// A std::vector has no memory manager behind it, so the host copy is
		// the only copy and saying so keeps the columns aligned with the
		// mfem::Vector lines above rather than leaving a reader to wonder.
		r.hostValid = true;
		r.deviceValid = false;
		r.hash = v.empty() ? 0
		                   : fnv1a( v.data(), v.size()*sizeof( double ) );
		double sum = 0.0;
		for ( double const x : v ) { sum += x*x; }
		r.norm = std::sqrt( sum );
		report( r );
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
