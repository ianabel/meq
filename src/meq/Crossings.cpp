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

}
