#ifndef MEQ_STAGE_CHECKSUM_HPP
#define MEQ_STAGE_CHECKSUM_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/**
 * @file
 *
 * WHERE A DEVICE RUN FIRST STOPS AGREEING WITH A HOST RUN, PRINTED PER STAGE.
 *
 * **THIS IS THE INSTRUMENT [M-187](../MEASUREMENTS.md) PRESCRIBES AND IT HAD
 * NOT BEEN KEPT.** `mfem::Device( "debug" )` catches a raw access to the
 * **wrong location** -- a protected page -- and is blind to honest flags over
 * a **stale copy**, which is the defect class M-187 and
 * [M-205](../MEASUREMENTS.md) both ended in. What separates those is a
 * host-side checksum taken at every stage of a step under `--device cpu` and
 * again under `--device cuda`, diffed: *the first stage that disagrees names
 * the defect*. M-187's walk built one by hand and threw it away; this is the
 * same thing, kept.
 *
 * **IT REPORTS THE VALIDITY FLAGS BEFORE IT READS, WHICH IS THE HALF A PLAIN
 * CHECKSUM LACKS.** `Memory<T>::HostIsValid()` and `DeviceIsValid()` are asked
 * *first*, so a stage whose host copy is not live is visible as such even
 * though the checksum that follows is correct -- and a raw reader at that
 * stage is precisely the bug. A line reading `H- D+` is a stage where anything
 * not going through a funnel gets stale data.
 *
 * **AND THE INSTRUMENT IS ITSELF A SYNC, WHICH IS A LIMITATION AND NOT A
 * DEFECT.** `HostRead()` migrates, so a checksum placed *before* a raw read
 * repairs the very staleness that read would have suffered: the run then gives
 * the right answer and the walk shows nothing. **So a stage whose flags read
 * `H- D+` is the finding, whether or not its checksum differs**, and probes
 * are placed after the writes rather than before the reads. Where a defect
 * vanishes once probes are inserted, that is evidence about the class, not a
 * failure of the walk.
 *
 * **OFF UNLESS `MEQ_STAGE_CHECKSUM` IS SET IN THE ENVIRONMENT**, and an
 * environment variable rather than a key because this is a debugging
 * instrument and not a property of a problem: it must be reachable on a
 * shipped `examples` TOML file without editing it, which is the whole point
 * of walking a case that already misbehaves. Reading it is cached, so a
 * disabled probe costs one predicted branch.
 *
 * **THIS HEADER IS MFEM-FREE, AND THE PROBES OVER MFEM'S CONTAINERS ARE IN
 * `meq/StageChecksumMfem.hpp`.** The walk has to cross the seam where MEQ
 * hands MFEM's world back to its own, so the `std::vector<double>` probe is
 * called from `src/meq/SafetyFactor.cpp`, which is in the MFEM-free half CI
 * builds. A probe over an `mfem::Vector` needs MFEM's header, and that file
 * must compile where there is none, so the two cannot share a header. The
 * readings, the switch and the printing are here, so both halves record into
 * one list and print one format; the `MfemFreeIncludes` test is what keeps
 * this header free of MFEM.
 */

namespace meq
{

	/// Whether the stage-checksum probes print anything. Reads
	/// `MEQ_STAGE_CHECKSUM` once and caches it.
	bool stageChecksumWanted();

	/**
	 * What one probe found, before and after the read it performs.
	 *
	 * `hostValid` and `deviceValid` are asked BEFORE the read, so they
	 * describe the state a raw accessor at this point would have met.
	 */
	struct StageReading
	{
		std::string stage;            ///< the caller's label
		int size = 0;                 ///< entries
		bool hostValid = false;       ///< before the read
		bool deviceValid = false;     ///< before the read
		std::uint64_t hash = 0;       ///< FNV-1a over the raw bytes
		double norm = 0.0;            ///< L2, so a difference has a size
	};

	/**
	 * Checksum a bare host-side vector of doubles under @a stage, print it
	 * when enabled, and record it.
	 *
	 * It has no memory manager behind it and therefore always reports
	 * `H+ D-`. It is here so that a walk can follow a quantity across the seam
	 * where MEQ hands MFEM's world back to its own -- `meq::SafetyFactorSystem`'s
	 * knot values and residuals are the case it is used for. The overloads over
	 * `mfem::Vector` and `mfem::Array<int>` are in `meq/StageChecksumMfem.hpp`.
	 *
	 * A no-op when `stageChecksumWanted()` is false, including the read -- so
	 * an instrumented build measured with the variable unset is the
	 * uninstrumented build.
	 */
	void stageChecksum( char const *stage, std::vector<double> const &v );

	/**
	 * FNV-1a over the raw bytes, the hash every probe records.
	 *
	 * The question a walk asks is *did this buffer change at all*, so the
	 * hash has to see the last bit of every mantissa -- a checksum built from
	 * sums or norms cancels, and two vectors differing by a sign flip in two
	 * entries share a norm exactly. Bytes of the IEEE representation have no
	 * such coincidence, and FNV-1a is four lines.
	 *
	 * It is NOT a claim about collision resistance and nothing here needs
	 * one: the comparison is between two runs of the same code on the same
	 * problem, where the alternative to "the hashes match" is a difference
	 * large enough that the L2 beside it shows it.
	 */
	std::uint64_t stageHash( void const *data, std::size_t bytes );

	/// Record @a reading and print it as one flushed `MEQ-STAGE` line: the
	/// step every probe ends in, whichever container it read. Called only
	/// when `stageChecksumWanted()` is true.
	void recordStageReading( StageReading const &reading );

	/// Every reading taken so far, in order.
	std::vector<StageReading> const &stageReadings();

	/// Forget them. Used between the arms of a comparison inside one process.
	void clearStageReadings();

}

#endif
