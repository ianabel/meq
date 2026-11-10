#pragma once

#include <functional>
#include <memory>

#include "mfem.hpp"

#include "meq/FluxExtraction.hpp"
#include "meq/FluxFamily.hpp"
#include "meq/FluxSurfaces.hpp"
#include "meq/GradShafranov.hpp"

#include "ConvergenceHarness.hpp"

/*
 * The fixture ExtractionContinuity asserts on and ExtractionJumps measures
 * with: a family of linear equilibria on the standard box, psi = 0 on its
 * edge, under the Solov'ev-type source F = -( ( 1 - A ) R^2 + A ), A the one
 * parameter. Moving A moves every flux surface smoothly, so the nodes of a
 * fixed family slide across element faces as A does, which is what makes a
 * broken extraction jump (MEASUREMENTS.md M-364).
 */
namespace meq::tests::crossing
{
	/// One solve at source mix @a mix on n x n cells at degree @a order, kept
	/// alive for its fields.
	class Solved
	{
		public:
			/// @a configure, when given, is applied to the solver before the
			/// solve: the precision instruments' device routes and refinement
			/// switch (RecoveryInstrument.hpp).
			Solved( double mix, int n, int order,
			        std::function<void( meq::GradShafranovSolver & )> const &configure = {} )
				: a( mix ),
				  mesh( makeMesh( standardBox(), n ) ),
				  source( [ this ]( mfem::Vector const &x )
				  {
					  return -( ( 1.0 - a )*x( 0 )*x( 0 ) + a );
				  } ),
				  zero( 0.0 ),
				  solver( mesh, order )
			{
				solver.setSource( source );
				solver.setBoundaryData( zero );
				namedRoutesFromEnvironment( solver );
				if ( configure )
					configure( solver );
				solver.solve();
				solver.postProcess();
			}

			Solved( Solved const & ) = delete;
			Solved &operator=( Solved const & ) = delete;

			/// Twelve surfaces of 96 nodes from Psi_N = 0.05 to 0.9, rooted on
			/// @a which.
			meq::FluxSurfaceFamily family( meq::Potential which ) const
			{
				meq::FluxFamilyOptions options;
				options.surfaces = 12;
				options.angles = 96;
				options.innerCut = 0.05;
				options.outerCut = 0.90;
				return meq::extractFluxSurfaces( solver, options, which );
			}

			meq::GradShafranovSolver const &theSolver() const
			{
				return solver;
			}

			mfem::Mesh &theMesh()
			{
				return mesh;
			}

		private:
			double a;
			mfem::Mesh mesh;
			mfem::FunctionCoefficient source;
			mfem::ConstantCoefficient zero;
			meq::GradShafranovSolver solver;
	};

	inline meq::FluxSurfaceFamily familyAt( double mix, int n, int order, meq::Potential which )
	{
		return Solved( mix, n, order ).family( which );
	}
}
