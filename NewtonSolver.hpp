#ifndef PICARDSOLVER_HPP
#define PICARDSOLVER_HPP

#include "mfem.hpp"
#include "Solution.hpp"
#include "GSSolver.hpp"


namespace meq {
	class NonlinearGSPicardSolver {
		protected:
			GSSolver solver;
			mfem::KINSolver *nonlinearSolver;
			int m_MaxIter;
			int m_AndersonLevels;
			int m_PrintLevel;
			double m_tol;
		public:
			std::shared_ptr<DGSpace> getSolutionSpace() { return solver.SolutionSpace; };
			NonlinearGSPicardSolver( std::shared_ptr<mfem::Mesh> mesh, int Order, GSSolver::Func RHS, double MaxIter = 1000, double AbsTol = 1e-3, int AndersonLevels = 0, int PrintLevel = 0 ) 
				: solver( mesh, Order, RHS ),
				m_MaxIter( MaxIter ),
				m_AndersonLevels( AndersonLevels ),
				m_PrintLevel( PrintLevel ),
				m_tol( AbsTol )
			{
				nonlinearSolver = new mfem::KINSolver( KIN_PICARD, false );
				nonlinearSolver->SetMaxIter( m_MaxIter );
				nonlinearSolver->SetAbsTol( m_tol );
				if ( m_AndersonLevels > 0 )
					nonlinearSolver->SetMAA( m_AndersonLevels );
				nonlinearSolver->iterative_mode = true;
				nonlinearSolver->SetPrintLevel( m_PrintLevel );
				nonlinearSolver->SetOperator( solver );
			};

			void Update()
			{
				delete nonlinearSolver;
				solver.Update();

				nonlinearSolver = new mfem::KINSolver( KIN_PICARD, false );
				nonlinearSolver->SetMaxIter( m_MaxIter );
				nonlinearSolver->SetAbsTol( m_tol );
				if ( m_AndersonLevels > 0 )
					nonlinearSolver->SetMAA( m_AndersonLevels );
				nonlinearSolver->iterative_mode = true;
				nonlinearSolver->SetPrintLevel( m_PrintLevel );
				nonlinearSolver->SetOperator( solver );
			}

			~NonlinearGSPicardSolver() 
			{
				delete nonlinearSolver;
			}

			void SetBCs( mfem::Coefficient& coeff ) {
				solver.SetBCs( coeff );
			};

			void Solve( Solution &solution )
			{
				nonlinearSolver->Mult( solution.qu, solution.qu );
				solver.Postprocess( solution );
			};
	};
}
#endif // PICARDSOLVER_HPP
