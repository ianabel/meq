
/*
 * An MFEM-based nonlinear poisson solver, derived from the HDG Example Poisson Solver
 * of T. Horvath, S. Rhebergen, and A. Sivas.
 *
 * The discretization is based on the paper:
 * N.C. Nguyen, J. Peraire, B. Cockburn, An implicit high-order hybridizable
 * discontinuous Galerkin method for linear convection–diffusion equations,
 * J. Comput. Phys., 2009, 228:9, 3232--3254.
 * Contributed by: T. Horvath, S. Rhebergen, A. Sivas
 *                University of Waterloo
 */

#include "meq.hpp"

#include "SolovievEquilibrium.hpp"
#include "GSSolver.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace mfem;
using namespace meq;

// Define the analytical solution and forcing terms / boundary conditions

static const double pi = 3.14159265358979323844;

// Solov'ev from Tonatiuh
const double A = -0.52;
const double C = 1.52;
const double c1 = -0.001479661575325;
const double c2 = -0.366568333204813;
const double c3 =  0.002409406149732;
const double c4 = -0.023957517168316;

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

int main(int argc, char *argv[])
{
	StopWatch chrono;

	// 1. Parse command-line options.
	int order = 3;
	int N_AMR=0;
	
	OptionsParser args(argc, argv);
	args.AddOption(&order, "-o", "--order",
		"Finite element order (polynomial degree).");
	args.AddOption(&N_AMR, "-amr", "--mesh-refinement-levels",
		"The number of loops of AMR to perform" );
	args.Parse();

	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);


	// Mesh generation
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(3, 3, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.5;
		constexpr double R_max = 1.5;
		constexpr double Z_min = -1;
		constexpr double Z_max =  1;

		out( 0 ) = R_min + in( 0 )*( R_max - R_min );
		out( 1 ) = Z_min + in( 1 )*( Z_max - Z_min );
	};

	mesh->Transform( xform );

	SolovievEquilibrium TestEq( A, C, c1, c2, c3, c4 );

	std::function<double( const mfem::Vector &, double )> J_Plasma = [ &TestEq ]( const mfem::Vector& pt, double psi ) {
		double R = pt( 0 );
		return -( R*R*TestEq.Pprime( psi ) + TestEq.FFprime( psi ) )/R;
	};

	NonlinearGSSolver Solver( mesh, order, J_Plasma );

	std::function<double( const mfem::Vector & )> uFun_ex = std::bind( &SolovievEquilibrium::psi, &TestEq, std::placeholders::_1 );
	std::function<void( const mfem::Vector &, mfem::Vector & )> qFun_ex = std::bind( &SolovievEquilibrium::q, &TestEq, std::placeholders::_1, std::placeholders::_2 );

	StdFunctionCoefficient bcFunCoeff( uFun_ex );
	Solver.SetBCs( bcFunCoeff );

	int i_amr;
	Solution soln( Solver.getSolutionSpace() );

	for (i_amr = 0; i_amr< N_AMR; i_amr++ )
	{
		Solver.Solve( soln );
		std::cout << "After " << i_amr << " levels of refinement:" << std::endl;
		auto [ err_u, err_q, err_u_star ] = soln.l2_errors( uFun_ex, qFun_ex );

		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;

		Solver.ApplyAMR( soln );

		std::cerr << "Refined Mesh" << std::endl;

	}

	{
		Solver.Solve( soln );
		std::cout << "After " << i_amr << " levels of refinement:" << std::endl;
		auto [ err_u, err_q, err_u_star ] = soln.l2_errors( uFun_ex, qFun_ex );

		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;
	}


	// 13. Save the mesh and the solution.
	{
		ofstream mesh_ofs("ex_hdg.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		soln.WriteOutputMFEM( "sol_q.gf", "sol_u.gf" );

	}

	// 14. Send the solution by socket to a GLVis server.
	{
		char vishost[] = "localhost";
		int  visport   = 19916;
		socketstream u_sock(vishost, visport);
		u_sock.precision(8);
		u_sock << "solution\n" << *mesh << soln.u_variable << "window_title 'Solution u'" <<
			endl;

		socketstream q_sock(vishost, visport);
		q_sock.precision(8);
		q_sock << "solution\n" << *mesh << soln.q_variable << "window_title 'Solution q'" <<
			endl;
	}

	return 0;
}


