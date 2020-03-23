
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
#include <boost/math/constants/constants.hpp>


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

template<typename TestEqClass> 
void TestSolution(  TestEqClass const& TestEq, std::shared_ptr<mfem::Mesh> mesh, int order )
{
	GSSolver Solver( mesh, order, TestEq );

	std::function<double( const mfem::Vector & )> uFun_ex = [ & ]( const mfem::Vector& pt ){ return TestEq.Psi( pt );};
	std::function<void( const mfem::Vector &, mfem::Vector & )> qFun_ex = [ & ]( const mfem::Vector& pt, mfem::Vector &q_out ) { TestEq.GradPsi( pt, q_out ); q_out *= -1.0/pt( 0 );};

	StdFunctionCoefficient bcFunCoeff( uFun_ex );

	Solver.SetBCs( bcFunCoeff );

	Solution soln( Solver.SolutionSpace );
	soln.Zero();
	Solver.Solve( soln );
	Solver.Postprocess( soln );
	auto [ err_u, err_q, err_u_star ] = soln.l2_errors( uFun_ex, qFun_ex );

	std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
	std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
	std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
	std::cout << std::endl;

}



int main(int argc, char *argv[])
{
	// 1. Parse command-line options.
	int order = 2;
	int N_REFINE = 4;
	
	OptionsParser args(argc, argv);
	args.AddOption(&order, "-o", "--order",
		"Finite element order (polynomial degree).");
	args.AddOption(&N_REFINE, "-r", "--refine",
		"Total uniform refinements" );
	args.Parse();

	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);


	SolovievEquilibrium TestEq( A, C, c1, c2, c3, c4 );
	
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(4, 4, Element::Type::TRIANGLE, false, 1.0, 1.0, true );

	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.25;
		constexpr double R_max = 1.0;
		constexpr double Z_min = -1;
		constexpr double Z_max = +1;
		out( 0 ) = R_min + in( 0 )*( R_max - R_min );
		out( 1 ) = Z_min + in( 1 )*( Z_max - Z_min );
	};

	mesh->Transform( xform );
	mesh->Finalize( true, true );

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " Soloviev (Polynomial Order " << order << ") : " << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	for ( unsigned int i_ref=0; i_ref < N_REFINE; i_ref++ )
	{
		std::cout <<"Ref " << i_ref << std::endl;
		TestSolution( TestEq, mesh, order );
		mesh->UniformRefinement();
	}

	return 0;
}


