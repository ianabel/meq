
#include "Refiners.hpp"

namespace mfem {

	MaximumMarkingRefiner::MaximumMarkingRefiner(ErrorEstimator &est) :
		estimator( est )
	{
		gamma = .5;
		max_elements = std::numeric_limits<long>::max();

		threshold = 0.0;
		num_marked_elements = 0L;
		current_sequence = -1;

		non_conforming = -1;
		nc_limit = 0;
	};

	int MaximumMarkingRefiner::ApplyImpl(Mesh &mesh)
	{
		threshold = 0.0;
		num_marked_elements = 0;
		marked_elements.SetSize(0);
		current_sequence = mesh.GetSequence();

		const long num_elements = mesh.GetGlobalNE();
		if (num_elements >= max_elements) { return STOP; }

		const int NE = mesh.GetNE();
		const Vector &local_err = estimator.GetLocalErrors();
		MFEM_ASSERT(local_err.Size() == NE, "invalid size of local_err");

		threshold = gamma * local_err.Max();	

		for (int el = 0; el < NE; el++)
		{
			if (local_err(el) > threshold)
			{
				marked_elements.Append(Refinement(el));
			}
		}

		num_marked_elements = mesh.ReduceInt(marked_elements.Size());
		if (num_marked_elements == 0) { return STOP; }

		mesh.GeneralRefinement(marked_elements, non_conforming, nc_limit);
		return CONTINUE + REFINED;
	}

	std::pair<int,int> Partition( mfem::Vector const& x, std::vector<int>& pi, int l, int u, int p )
	{
	}

	int QuickMark( mfem::Vector const& x, std::vector<int> & pi, int l, int u, double v )
	{
		double x_old_p = x( pi[ p ] ); 
		// We never need pi_old & pi_new at the same time. 
		// so we modify in place

		int p = Pivot( x, pi, l, u ); // Leaves pi fixed
		int g,s;
		[ g, s ] = Partition( x, pi, l, u, p ); // Updates pi, so we no longer have pi_old
		double sigma_g = 0;
		for ( int j = l; j <= g; ++j )
			sigma_g += x( pi[ j ] );
		if ( sigma_g >= v )
			return QuickMark( x, pi, l, g, v );
		if ( ( sigma_g + ( s - g - 1 )*x_old_p ) >= v )
			// We've updated pi already, so can just return the index.
			return g + std::ceil( ( v - sigma_g )/x_old_p );

		return QuickMark( x, pi, s, u, v - sigma_g - ( s - g - 1 )*x_old_p );
	}

	int DoerflerMarking::ApplyImpl(Mesh &mesh)
	{
		threshold = 0.0;
		num_marked_elements = 0;
		marked_elements.SetSize(0);
		current_sequence = mesh.GetSequence();

		const long num_elements = mesh.GetGlobalNE();
		if (num_elements >= max_elements) { return STOP; }

		const int NE = mesh.GetNE();
		const Vector &local_err = estimator.GetLocalErrors();
		MFEM_ASSERT(local_err.Size() == NE, "invalid size of local_err");

		threshold = gamma * local_err.Max();	

		for (int el = 0; el < NE; el++)
		{
			if (local_err(el) > threshold)
			{
				marked_elements.Append(Refinement(el));
			}
		}

		num_marked_elements = mesh.ReduceInt(marked_elements.Size());
		if (num_marked_elements == 0) { return STOP; }

		mesh.GeneralRefinement(marked_elements, non_conforming, nc_limit);
		return CONTINUE + REFINED;
	}
}
