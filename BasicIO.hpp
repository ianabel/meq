
#include "meq.hpp"

namespace meq {
	void WriteAsciiPSI( std::string const& filename, mfem::GridFunction const& psi, double R_min, double R_max, double Z_min, double Z_max, unsigned int NR, unsigned int NZ );
}
