
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
#include "PoissonSolver.hpp"

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

double tol = 1e-4;

void TestSolution( std::shared_ptr<mfem::Mesh> mesh, int order )
{
	auto RHS = [ & ]( const mfem::Vector &pt, double u ) {
		return fFun( pt );
	};

	PoissonSolver Solver( mesh, order, RHS );

	FunctionCoefficient bcFunCoeff( uFun_ex );

	Solver.SetBCs( bcFunCoeff );

	Solution soln( Solver.SolutionSpace );
	soln.Zero(); // Initial guess is 0

	{
		Solver.Solve( soln );
		Solver.Postprocess( soln );
		auto [ err_u, err_q, err_u_star ] = soln.l2_errors( uFun_ex, qFun_ex );

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
	
	OptionsParser args(argc, argv);
	args.AddOption(&order, "-o", "--order",
		"Finite element order (polynomial degree).");
	args.AddOption(&N_REFINE, "-r", "--refine",
		"Initial uniform refinements" );
	args.Parse();

	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// Mesh generation
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(4, 4, Element::Type::TRIANGLE, false, 1.0, 1.0, true );

	for ( int i=0; i<N_REFINE; i++ )
		mesh->UniformRefinement();

	TestSolution( mesh, order );

	return 0;
}

static const int n = 2;
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
