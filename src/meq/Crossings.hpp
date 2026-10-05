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
 * **AND IT COUNTS BOTH DIRECTIONS, WHICH IT DID NOT UNTIL PHASE B NEEDED IT.**
 * DEVICE-PLAN.md section 3.3 divides a resident solve into *setup · one upload ·
 * device only · one download · output*, and everything above counts host
 * READS -- downloads. Phase B is an UPLOAD. *An instrument that can see one
 * direction of the bus cannot measure a plan whose first phase is an upload*,
 * and stage 4b reported a figure in the other direction and called it the
 * residency acceptance for as long as it existed.
 * **THE PER-SITE TABLE IS KEYED BY A LINE NUMBER ALONE, AND THAT IS A KNOWN
 * LIMITATION RATHER THAN AN OVERSIGHT.** Call sites pass `__LINE__` and not
 * `__FILE__`, so two funnels at the same line in different files are summed
 * into one row -- and funnels are registered from GradShafranov.cpp,
 * CriticalPoints.cpp, PointLocator.cpp, DeviceBasis.cpp and
 * AxisSearchDevice.cpp. **The totals, the two directions and the phase split
 * are exact**; only the attribution of one ROW can merge, and a reader chasing
 * a row should check the line in more than one file. Fixing it means a
 * `( file, line )` key and `__FILE__` at every call site, which is forty
 * edits for a table that is read by hand.
 *
 * `meq::countedDeviceRead()` is the mirror: a device read moves data exactly
 * when the HOST holds the live copy and the device does not.
 *
 * **AND EVERY CROSSING IS ATTRIBUTED TO A PHASE, BECAUSE THE SAME TRANSFER
 * MEANS OPPOSITE THINGS IN TWO OF THEM.** A transfer in phase B is the design
 * and a transfer in phase C is the defect, so a single total is an average over
 * populations with opposite meanings -- which is how the published
 * `0.2 transfers a step` turned out to be three transfers in phase A, divided
 * by a Newton step count. **The two directions are needed TOGETHER**: phase C's
 * transfer count reads the same zero for a perfectly resident solve and for one
 * where nothing reaches the device, and only the upload count beside it tells
 * those apart. M-226.
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

	/// Funnelled accesses of a buffer the other side may have written, and how
	/// many of them actually moved data, in each direction.
	struct CrossingCount
	{
		long reads = 0;      ///< host-read funnel calls
		long transfers = 0;  ///< those where the device copy was the live one
		long entries = 0;    ///< doubles moved, summed over the transfers

		long deviceReads = 0;   ///< device-read funnel calls
		long uploads = 0;       ///< those where the HOST copy was the live one
		long uploadEntries = 0; ///< doubles moved, summed over the uploads

		void clear()
		{
			reads = 0;
			transfers = 0;
			entries = 0;
			deviceReads = 0;
			uploads = 0;
			uploadEntries = 0;
		}

		void add( CrossingCount const &other )
		{
			reads += other.reads;
			transfers += other.transfers;
			entries += other.entries;
			deviceReads += other.deviceReads;
			uploads += other.uploads;
			uploadEntries += other.uploadEntries;
		}
	};

	/**
	 * The residency phase a crossing belongs to -- DEVICE-PLAN.md section 3.3's
	 * own table, which is the shape every item of that plan is ordered against.
	 *
	 * **A TRANSFER IN PHASE B IS THE DESIGN AND THE SAME TRANSFER IN PHASE C IS
	 * THE DEFECT, AND ONE TOTAL CANNOT TELL THEM APART.** The plan says *one
	 * upload*, *device only*, *one download* -- so the number that decides
	 * whether residency is arriving is a count taken BETWEEN the upload and the
	 * download, and a headline averaged over all five phases is an average over
	 * populations with opposite meanings. Stage 4b's published
	 * `0.2 transfers a step` is such an average.
	 */
	enum class Phase
	{
		Setup,     ///< A -- host, and it may take as long as it likes
		Upload,    ///< B -- configure the device problem
		Solve,     ///< C -- device only
		Download,  ///< D -- the answer returns, once
		Output     ///< E -- host
	};

	/// How many there are, for the per-phase tallies.
	int constexpr phaseCount = 5;

	/// The phase crossings are currently attributed to. Setup until something
	/// says otherwise, which is the truthful answer for a process that has not
	/// reached a solve.
	Phase &residencyPhase();

	/// Per-phase tallies, indexed by `static_cast<int>( Phase )`.
	CrossingCount *crossingPhases();

	/// The name of @a phase, for a report.
	char const *phaseName( Phase phase );

	/**
	 * Attribute every crossing to @a phase for as long as this object lives,
	 * and restore the previous phase afterwards.
	 *
	 * **RESTORING THE PREVIOUS PHASE RATHER THAN RETURNING TO A DEFAULT IS WHAT
	 * MAKES THIS SAFE UNDER RE-ENTRY**, and the solver re-enters:
	 * `solveByPicardThenNewton()` calls `solve()` for its second stage, so a
	 * plain set-then-clear would have the inner call's exit put the outer one
	 * into the wrong phase for the rest of its run. That is the same hazard
	 * `GradShafranovSolver::solve()` records about `profileStart` being a local
	 * rather than a member, met one level up.
	 */
	class PhaseScope
	{
		public:
			explicit PhaseScope( Phase phase )
				: previous( residencyPhase() )
			{
				residencyPhase() = phase;
			}

			~PhaseScope()
			{
				residencyPhase() = previous;
			}

			PhaseScope( PhaseScope const & ) = delete;
			PhaseScope &operator=( PhaseScope const & ) = delete;

		private:
			Phase previous;
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
	 * The shared body of every host funnel: tally in all three places.
	 *
	 * **THREE TALLIES AND NOT ONE, BECAUSE EACH ANSWERS A QUESTION THE OTHERS
	 * CANNOT.** The process-wide total is the headline; the per-site tally
	 * names what to do about it, which is why `crossingSites()` exists at all;
	 * and the per-phase tally says whether a transfer was the DESIGN or the
	 * DEFECT, which neither of the others can distinguish. They are updated
	 * together, in one function, for `meq::limiterPolygonMaximum()`'s reason:
	 * two expressions of one quantity kept in step by hand diverge the moment
	 * one of them moves.
	 */
	template <typename T>
	inline void recordHostRead( mfem::Memory<T> const &memory, int size, int line )
	{
		CrossingCount &count = crossings();
		CrossingCount &site = crossingSites()[ line ];
		CrossingCount &phase =
			crossingPhases()[ static_cast<int>( residencyPhase() ) ];
		++count.reads;
		++site.reads;
		++phase.reads;

		// A host read moves data exactly when the device holds the live copy.
		// Asked BEFORE the read, because the read is what changes the answer.
		if ( memory.DeviceIsValid() && !memory.HostIsValid() )
		{
			++count.transfers;
			count.entries += size;
			++site.transfers;
			site.entries += size;
			++phase.transfers;
			phase.entries += size;
		}
	}

	/**
	 * And its mirror, for a read ON the device.
	 *
	 * **THE CONDITION IS THE OTHER WAY ROUND AND THAT IS THE WHOLE OF THE
	 * DIFFERENCE**: a device read moves data exactly when the HOST holds the
	 * live copy and the device does not. Asked before the read, for the same
	 * reason -- `Read()` is what makes the answer stop being true.
	 */
	template <typename T>
	inline void recordDeviceRead( mfem::Memory<T> const &memory, int size, int line )
	{
		CrossingCount &count = crossings();
		CrossingCount &site = crossingSites()[ line ];
		CrossingCount &phase =
			crossingPhases()[ static_cast<int>( residencyPhase() ) ];
		++count.deviceReads;
		++site.deviceReads;
		++phase.deviceReads;

		if ( memory.HostIsValid() && !memory.DeviceIsValid() )
		{
			++count.uploads;
			count.uploadEntries += size;
			++site.uploads;
			site.uploadEntries += size;
			++phase.uploads;
			phase.uploadEntries += size;
		}
	}

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
		recordHostRead( v.GetMemory(), v.Size(), line );
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
		recordHostRead( v.GetMemory(), v.Size(), line );
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
		recordHostRead( a.GetMemory(), a.Size(), line );
		a.HostRead();
	}

	/**
	 * Count a device read of @a v and then perform it, returning the device
	 * pointer exactly as `Read()` would.
	 *
	 * **THIS DIRECTION HAD NO COUNTER AT ALL, AND PHASE B IS DEFINED ENTIRELY
	 * IN TERMS OF IT.** DEVICE-PLAN.md section 3.3 asks for *one upload*, and
	 * stage 4b's acceptance -- the one number this project trusts about
	 * residency -- counted host reads only. *An instrument that can see one
	 * direction of the bus cannot measure a plan whose first phase is an
	 * upload*, which is the same species of coverage defect as
	 * [M-198](../MEASUREMENTS.md)'s first crossing counter covering 17 funnels
	 * of 36 and reporting a number 18x too low.
	 *
	 * Replaces a bare `v.Read()` at an upload site and is otherwise identical.
	 */
	inline mfem::real_t const *countedDeviceRead( mfem::Vector const &v,
	                                              int line = 0 )
	{
		recordDeviceRead( v.GetMemory(), v.Size(), line );
		return v.Read();
	}

	/**
	 * The same for a read-WRITE: uploads on the same condition and
	 * additionally invalidates the HOST copy, which is what makes a vector
	 * `H- D+` -- the state every residency gate in `GradShafranov.cpp` asks
	 * for. A `Read()` leaves both copies valid and opens none of them.
	 */
	inline mfem::real_t *countedDeviceReadWrite( mfem::Vector &v, int line = 0 )
	{
		recordDeviceRead( v.GetMemory(), v.Size(), line );
		return v.ReadWrite();
	}

	/// The same for an `mfem::Array<T>`, for the reason the host overload gives.
	template <typename T>
	inline T const *countedDeviceRead( mfem::Array<T> const &a, int line = 0 )
	{
		recordDeviceRead( a.GetMemory(), a.Size(), line );
		return a.Read();
	}

}

#endif // MEQ_CROSSINGS_HPP
