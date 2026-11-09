#ifndef MEQ_DIGEST_HPP
#define MEQ_DIGEST_HPP

#include <cstdint>

/**
 * @file
 *
 * THE DIGEST A CONDUCTOR CACHE IS KEYED ON, AND NOTHING ELSE: FNV-1a over the
 * IEEE bytes of what is mixed in.
 *
 * Exposed because the mesh and conductor digests are built by callers that own
 * those objects -- the solver knows its mesh, meq::ConductorField knows its
 * conductors -- and a digest is only comparable if both sides mix in the same
 * order. `meq/ConductorStore.hpp` includes this and is where the digest is
 * written and compared.
 *
 * **ITS OWN TRANSLATION UNIT, MFEM-FREE AND NetCDF-FREE**, because
 * `src/meq/ConductorField.cpp` is in both of those halves and computes a
 * digest, while `src/meq/ConductorStore.cpp`, the reader and writer, needs
 * netcdf-cxx4 and is left out of a build without it -- which CI's is.
 *
 * A digest is only ever compared with another digest written by this same
 * code, so the IEEE bytes are the right thing to mix -- and -0.0 and +0.0
 * hashing differently is CORRECT here, because a vertex at -0.0 and one at
 * +0.0 are the same point and `psi_c` is even in `R`, so the values agree and
 * a spurious refusal costs a recompute rather than an answer.
 */

namespace meq
{

	std::uint64_t digestSeed();
	std::uint64_t digestAppend( std::uint64_t digest, double value );
	std::uint64_t digestAppend( std::uint64_t digest, std::int64_t value );

}

#endif
