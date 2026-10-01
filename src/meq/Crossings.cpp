#include "Crossings.hpp"

namespace meq
{

	CrossingCount &crossings()
	{
		// A function-local static rather than a namespace-scope object, for the
		// reason tests/unit/DeviceRowsTests.cpp records about mfem::Device: an
		// object at namespace scope may be constructed before something it
		// depends on, and a counter that is read during static initialisation
		// would be reading a zero it had not earned.
		static CrossingCount count;
		return count;
	}

	std::map<int, CrossingCount> &crossingSites()
	{
		static std::map<int, CrossingCount> sites;
		return sites;
	}

	Phase &residencyPhase()
	{
		// Setup until something says otherwise, which is the truthful answer
		// for a process that has not reached a solve -- and the right default
		// for a test binary that never will.
		static Phase phase = Phase::Setup;
		return phase;
	}

	CrossingCount *crossingPhases()
	{
		static CrossingCount phases[ phaseCount ];
		return phases;
	}

	char const *phaseName( Phase phase )
	{
		switch ( phase )
		{
			case Phase::Setup:    return "A setup";
			case Phase::Upload:   return "B upload";
			case Phase::Solve:    return "C solve";
			case Phase::Download: return "D download";
			case Phase::Output:   return "E output";
		}
		// Unreachable for any value of the enum, and a default: would stop
		// -Wswitch warning when a phase is added -- which CLAUDE.md records
		// costing a silent wrong trace solver at apps/meq.cpp's own switch.
		return "unknown";
	}

}
