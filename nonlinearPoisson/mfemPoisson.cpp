
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

#include "mfem.hpp"
#include "NonlinearSolver.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/math/constants/constants.hpp>


using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions

static constexpr double pi = boost::math::double_constants::pi;

double uFun_ex(const Vector & x);
void qFun_ex(const Vector & x, Vector & q);
double fFun(const Vector & x);

int N_ANDERSON = 2;
int N_MAX_ITER = 1500;
double tol = 1e-4;

void TestSolution( std::shared_ptr<mfem::Mesh> mesh, int order )
{
	auto RHS = [ & ]( const mfem::Vector &pt, double u ) {
		return ( u*u - uFun_ex( pt )*uFun_ex( pt ) ) + fFun( pt );
	};

	NonlinearPoissonSolver Solver( mesh, order, RHS, N_MAX_ITER, tol, N_ANDERSON, 0 );

	FunctionCoefficient bcFunCoeff( uFun_ex );

	Solver.SetBCs( bcFunCoeff );

	Solution soln( Solver.getSolutionSpace() );
	soln.Zero(); // Initial guess is 0

	{
		Solver.Solve( soln );
		auto [ err_u, err_q, err_u_star ] = soln.l2_errors( uFun_ex, qFun_ex );

		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;
	}

	std::cout << " Refining Grid and Prolonging Solution " << std::endl;

	mesh->UniformRefinement();
	mesh->Finalize( true, true );
	Solver.Update();
	soln.Prolong();

	{
		double refined_tol = ::pow( tol, 1.5 );
		if ( refined_tol < 1e-12 )
			refined_tol = 1e-12;

		NonlinearPoissonSolver refined_solver( mesh, order, RHS, N_MAX_ITER, refined_tol, N_ANDERSON, 0 );
		Solution refined_soln( refined_solver.getSolutionSpace(), soln.qu );

		refined_solver.SetBCs( bcFunCoeff );
		refined_solver.Solve( refined_soln );

		auto [ err_u, err_q, err_u_star ] = refined_soln.l2_errors( uFun_ex, qFun_ex );

		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;
	}

	{
		std::cout << " Visualising" << std::endl;
		char vishost[] = "localhost";
		int  visport   = 19916;
		socketstream u_star_sock(vishost, visport);
		u_star_sock.precision(8);
		u_star_sock << "solution\n" << *mesh << soln.u_star_variable <<
			"window_title 'Solution u_star'" << endl;
	}
}

int main(int argc, char *argv[])
{
	// 1. Parse command-line options.
	int order = 3;
	int N_REFINE=0;
   const char *mesh_file = "./square.msh";
	
	OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
	args.AddOption(&order, "-o", "--order",
		"Finite element order (polynomial degree).");
	args.AddOption(&N_REFINE, "-r", "--refine",
		"Initial uniform refinements" );
	args.AddOption(&N_ANDERSON, "-a", "--maa",
		"Anderson Acceleration Level" );
	args.AddOption(&tol, "-t", "--tolerance",
		"Fixed-point Acceleration Tolerance" );
	args.Parse();

	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// Mesh generation
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>( mesh_file, 1.0, 1.0 );

	for ( int i=0; i<N_REFINE; i++ )
		mesh->UniformRefinement();

	TestSolution( mesh, order );

	return 0;
}

static const int n = 4;
static const double diff = 1.0;

double uFun_ex(const Vector & x)
{
   double xi(x(0));
   double yi(x(1));

	return 1.0 + xi + sin(n*pi*xi)*sin(n*pi*yi);
}

void qFun_ex(const Vector & x, Vector & q)
{
	double xi(x(0));
	double yi(x(1));
	{
		q(0) = -diff*1.0 - diff*n*pi*cos(n*pi*xi)*sin(n*pi*yi);
		q(1) =           - diff*n*pi*sin(n*pi*xi)*cos(n*pi*yi);
	}
}


double fFun(const Vector & x)
{
	double xi(x(0));
	double yi(x(1));
	return diff * ( 2.0 * n * n * pi * pi )*sin(n * pi * xi)*sin(n * pi * yi);
}
