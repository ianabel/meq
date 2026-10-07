#include "Inverse.hpp"

#include <vector>

#include "SourceFactory.hpp"

namespace meq
{
	void poseInverse( GradShafranovSolver &solver, std::shared_ptr< CoilSet > coils,
	                  CoilConfig const &config, InverseConfig const &inverse )
	{
		std::vector< FreeCircuit > circuits = makeFreeCircuits( config );
		std::vector< double > gamma( circuits.size(), inverse.gamma );
		solver.clearInverse();
		solver.setFreeCurrents( std::move( coils ), std::move( circuits ),
		                        std::move( gamma ) );
		for ( IsofluxTargetConfig const &t : inverse.isoflux )
			solver.addIsofluxTarget( t.r1, t.z1, t.r2, t.z2 );
		for ( XPointTargetConfig const &t : inverse.xpoint )
		{
			if ( t.track )
				solver.addXPointTarget( t.radius, t.z );
			else
				solver.addCriticalPointTarget( t.radius, t.z );
		}
	}
}
