#ifndef NONLINEARSOLVER_HPP
#define NONLINEARSOLVER_HPP

#include "mfem.hpp"
#include "PoissonSolver.hpp"

class NonlinearPoissonSolver {
	protected:
		PoissonSolver solver;
		mfem::KINSolver *nonlinearSolver;
		int m_MaxIter;
		int m_AndersonLevels;
		int m_PrintLevel;
		double m_tol;
	public:
		std::shared_ptr<DGSpace> getSolutionSpace() { return solver.SolutionSpace; };
		NonlinearPoissonSolver( std::shared_ptr<mfem::Mesh> mesh, int Order, PoissonSolver::Func RHS, double MaxIter = 1000, double AbsTol = 1e-3, int AndersonLevels = 0, int PrintLevel = 0 ) 
			: solver( mesh, Order, RHS ),
			  m_MaxIter( MaxIter ),
			  m_AndersonLevels( AndersonLevels ),
			  m_PrintLevel( PrintLevel ),
			  m_tol( AbsTol )
		{
			nonlinearSolver = new mfem::KINSolver( KIN_FP, false );
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

			nonlinearSolver = new mfem::KINSolver( KIN_FP, false );
			nonlinearSolver->SetMaxIter( m_MaxIter );
			nonlinearSolver->SetAbsTol( m_tol );
			if ( m_AndersonLevels > 0 )
				nonlinearSolver->SetMAA( m_AndersonLevels );
			nonlinearSolver->iterative_mode = true;
			nonlinearSolver->SetPrintLevel( m_PrintLevel );
			nonlinearSolver->SetOperator( solver );
		}

		~NonlinearPoissonSolver() 
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

#endif // NONLINEARSOLVER_HPP

