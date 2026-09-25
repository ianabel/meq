#ifndef MEQ_CROSSINGS_HPP
#define MEQ_CROSSINGS_HPP

#include "mfem.hpp"

#include <map>

/**
 * @file
 *
 * HOW OFTEN THE HOST TAKES A STATE THE DEVICE WROTE, COUNTED RATHER THAN
 * INFERRED.
 *
 * **THIS IS THE ONE ACCEPTANCE A MACHINE LIKE THIS CAN HONOUR, AND
 * DEVICE-PLAN.md §4 STAGE 4b IS WHERE IT IS ASKED FOR.** Correctness cannot
 * tell a resident step from a partial one: [M-187](../MEASUREMENTS.md) is a run
 * that agreed with the host to every printed digit while doing **4× the work**.
 * Nor can a clock here — `CLAUDE.md` forbids reading throughput off this card,
 * and a reduction that is correct and does not yet pay is invisible to a timer
 * and visible to this. So the number that says whether residency is arriving is
 * a **count of crossings per Newton step**, and it needs no GPU.
 *
 * **IT COUNTS TRANSFERS AND NOT CALLS, WHICH IS THE WHOLE DIFFERENCE BETWEEN AN
 * INSTRUMENT AND A PROXY.** `DEVICE-PLAN.md` §3.3's table derives its 21.1 a
 * step from leg CALL counts, and says in terms that a count of reads is not a
 * count of transfers: once any read has brought the state to the host it stays
 * valid until the device writes again. MFEM exposes the distinction --
 * `mfem::Memory<T>::DeviceIsValid()` and `HostIsValid()` are public, and a host
 * read moves data exactly when the device copy is valid and the host copy is
 * not. So both are counted and both are reported, and their RATIO is the
 * interesting number: reads far above transfers means the funnels are
 * redundant, transfers approaching reads means the device is writing between
 * them, which is what a resident step looks like.
 *
 * **THE COUNTER IS PROCESS-WIDE AND THAT IS A DECISION.** The property it
 * measures is process-wide -- the memory manager's validity flags -- and the
 * alternative is threading a profile pointer through `CompressedRows::gather()`,
 * a free function whose whole reason to exist is being cheap. A `StepProfile`
 * member would also miss every crossing outside a solve, which is exactly where
 * the output stage's are.
 *
 * **AND A PLAIN `long` IS SAFE HERE, WHICH IS AN ARGUMENT RATHER THAN AN
 * OVERSIGHT.** Two arms and neither races: with no `mfem::Device` configured
 * nothing is ever device-valid, so `transfers` never increments and `reads` is
 * touched only from loops MEQ already runs serially; with one configured,
 * `meq::hostLoopsAreSafe()` gates every MEQ element loop to a single thread
 * ([M-188](../MEASUREMENTS.md)). **If that gate is ever removed, this needs an
 * atomic** -- and the gate's own removal is what would make the counter wrong
 * rather than anything here.
 */

namespace meq
{

	/// Host reads of a buffer the device may have written, and how many of
	/// them actually moved data.
	struct CrossingCount
	{
		long reads = 0;      ///< funnel calls
		long transfers = 0;  ///< those where the device copy was the live one
		long entries = 0;    ///< doubles moved, summed over the transfers

		void clear()
		{
			reads = 0;
			transfers = 0;
			entries = 0;
		}
	};

	/// The process-wide counter. See the file comment for why it is not on
	/// meq::GradShafranovSolver::StepProfile.
	CrossingCount &crossings();

	/**
	 * Per-site tallies, keyed by the line the funnel is on.
	 *
	 * **A SINGLE NUMBER SAYS RESIDENCY IS OR IS NOT ARRIVING AND NAMES NOTHING
	 * TO DO ABOUT IT**, which is the difference between a score and an
	 * instrument. The key is `__LINE__`, generated at compile time, so it
	 * cannot rot the way a line number written into prose does -- `CLAUDE.md`
	 * records a table whose line numbers were eight hundred lines wrong.
	 */
	std::map<int, CrossingCount> &crossingSites();

	/// Count a host read of @a v, attributing it to @a line, then perform it.
	void countedHostReadAt( mfem::Vector const &v, int line );

	/**
	 * Count a host read of @a v and then perform it.
	 *
	 * **THE COUNT AND THE READ ARE ONE CALL SO THEY CANNOT DRIFT APART**, which
	 * is `meq::limiterPolygonMaximum()`'s lesson and
	 * [M-189](../MEASUREMENTS.md)'s: two expressions of one quantity, kept in
	 * step by hand, diverge the moment one of them moves. A funnel that counted
	 * separately from the read it describes is the same shape.
	 *
	 * Replaces a bare `v.HostRead()` at a funnel and is otherwise identical to
	 * it, including being valid on a `const` vector.
	 */
	inline void countedHostRead( mfem::Vector const &v, int line = 0 )
	{
		CrossingCount &count = crossings();
		CrossingCount &site = crossingSites()[ line ];
		++count.reads;
		++site.reads;

		// A host read moves data exactly when the device holds the live copy.
		// Asked BEFORE the read, because the read is what changes the answer.
		mfem::Memory<mfem::real_t> const &memory = v.GetMemory();
		if ( memory.DeviceIsValid() && !memory.HostIsValid() )
		{
			++count.transfers;
			count.entries += v.Size();
			++site.transfers;
			site.entries += v.Size();
		}

		v.HostRead();
	}

	/**
	 * The same for a read-WRITE, which moves data on the same condition and
	 * additionally invalidates the device copy.
	 *
	 * Counted separately in neither field: a funnel is a funnel, and what the
	 * instrument is for is how often the host takes a state the device wrote.
	 * That a read-write also *keeps* it on the host is a property of the call
	 * site and shows up as the NEXT read being free.
	 */
	inline void countedHostReadWrite( mfem::Vector &v, int line = 0 )
	{
		CrossingCount &count = crossings();
		CrossingCount &site = crossingSites()[ line ];
		++count.reads;
		++site.reads;

		mfem::Memory<mfem::real_t> const &memory = v.GetMemory();
		if ( memory.DeviceIsValid() && !memory.HostIsValid() )
		{
			++count.transfers;
			count.entries += v.Size();
			++site.transfers;
			site.entries += v.Size();
		}

		v.HostReadWrite();
	}

	/**
	 * And for an `mfem::Array<T>`.
	 *
	 * **AN ARRAY IS DEVICE STATE ON THE SAME FOOTING AS A VECTOR**, which is
	 * the thing that is easy to miss and which `meq::CompressedRows`' own
	 * comment about `rowDot()` records: MFEM's device-aware operations take dof
	 * arrays by `Array::Read( use_dev )`, so one handed to a kernel comes back
	 * device-valid, and `Array::operator[]` is as raw as `Vector::operator()`.
	 * A crossing counter that saw only Vectors would miss them.
	 */
	template <typename T>
	inline void countedHostRead( mfem::Array<T> const &a, int line = 0 )
	{
		CrossingCount &count = crossings();
		CrossingCount &site = crossingSites()[ line ];
		++count.reads;
		++site.reads;

		mfem::Memory<T> const &memory = a.GetMemory();
		if ( memory.DeviceIsValid() && !memory.HostIsValid() )
		{
			++count.transfers;
			count.entries += a.Size();
			++site.transfers;
			site.entries += a.Size();
		}

		a.HostRead();
	}

}

#endif // MEQ_CROSSINGS_HPP
