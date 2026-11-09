#include "Digest.hpp"

#include <cstring>

namespace meq
{
	namespace
	{
		constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ULL;
		constexpr std::uint64_t fnvPrime = 1099511628211ULL;

		std::uint64_t mixBytes( std::uint64_t digest, std::uint64_t bits )
		{
			for ( int byte = 0; byte < 8; ++byte )
			{
				digest ^= ( bits >> ( 8*byte ) ) & 0xffULL;
				digest *= fnvPrime;
			}
			return digest;
		}
	}

	std::uint64_t digestSeed()
	{
		return fnvOffsetBasis;
	}

	std::uint64_t digestAppend( std::uint64_t digest, double value )
	{
		std::uint64_t bits = 0;
		std::memcpy( &bits, &value, sizeof( bits ) );
		return mixBytes( digest, bits );
	}

	std::uint64_t digestAppend( std::uint64_t digest, std::int64_t value )
	{
		return mixBytes( digest, static_cast< std::uint64_t >( value ) );
	}

}
