#ifndef MEQ_STAGE_CHECKSUM_MFEM_HPP
#define MEQ_STAGE_CHECKSUM_MFEM_HPP

#include "meq/StageChecksum.hpp"

#include "mfem.hpp"

/**
 * @file
 *
 * THE STAGE-CHECKSUM PROBES OVER MFEM'S CONTAINERS, which are the ones that
 * carry a memory manager and therefore validity flags worth reporting.
 *
 * The instrument itself -- the switch, the readings, the printed line and why
 * the flags are asked before the read -- is `meq/StageChecksum.hpp`, which is
 * MFEM-free so that the half of MEQ CI builds can probe the seam back into
 * plain doubles. These overloads record into the same list.
 */

namespace meq
{

	/**
	 * Checksum @a v under @a stage, print it when enabled, and record it.
	 *
	 * A no-op when `stageChecksumWanted()` is false, including the read -- so
	 * an instrumented build measured with the variable unset is the
	 * uninstrumented build.
	 */
	void stageChecksum( char const *stage, mfem::Vector const &v );

	/// The same for an integer array, which `Array<int>::operator[]` reads as
	/// rawly as `Vector::operator()` does.
	void stageChecksum( char const *stage, mfem::Array<int> const &a );

}

#endif
