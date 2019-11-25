#ifndef NONLINEARSOLVER_HPP
#define NONLINEARSOLVER_HPP

#include "mfem.hpp"
#include "Solution.hpp"
#include "GSSolver.hpp"

namespace meq {
	class NonlinearGSSolver {
		protected:
			GSSolver solver;
			mfem::KINSolver *nonlinearSolver;
		public:
			std::shared_ptr<DGSpace> getSolutionSpace() { return solver.SolutionSpace; };
			NonlinearGSSolver( std::shared_ptr<mfem::Mesh> mesh, int Order, GSSolver::Func PlasmaRHS ) 
				: solver( mesh, Order, PlasmaRHS )
			{
				nonlinearSolver = new mfem::KINSolver( KIN_FP, false );
				nonlinearSolver->SetMaxIter( 1000 );
				nonlinearSolver->SetAbsTol( 1e-4 );
				nonlinearSolver->iterative_mode = false;
				nonlinearSolver->SetOperator( solver );
			};

			~NonlinearGSSolver() 
			{
				delete nonlinearSolver;
			}

			void ApplyAMR( Solution &soln )
			{
				solver.ApplyAdaptiveRefinement( soln );
				soln.Reset();
				nonlinearSolver->SetOperator( solver );
			};

			void SetBCs( mfem::Coefficient& coeff ) {
				solver.SetBCs( coeff );
			};

			void Solve( Solution &solution )
			{
				nonlinearSolver->Mult( solution.qu, solution.qu );
				solver.Postprocess( solution );
			};
	};

};

#endif // NONLINEARSOLVER_HPP

