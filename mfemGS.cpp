
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
#include "TonatiuhManufacturedEq1.hpp"
#include "AnalyticTestEq2.hpp"
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

int N_ANDERSON = 2;
int N_MAX_ITER = 4500;
double tol = 1e-4;

template<typename TestEqClass> 
void TestSolution(  TestEqClass const& TestEq, std::shared_ptr<mfem::Mesh> mesh, int order )
{
	NonlinearGSSolver Solver( mesh, order, TestEq, N_MAX_ITER, tol, N_ANDERSON );

	std::function<double( const mfem::Vector & )> uFun_ex = [ & ]( const mfem::Vector& pt ){ return TestEq.Psi( pt );};
	std::function<void( const mfem::Vector &, mfem::Vector & )> qFun_ex = [ & ]( const mfem::Vector& pt, mfem::Vector &q_out ) { TestEq.GradPsi( pt, q_out ); q_out *= -1.0/pt( 0 );};

	StdFunctionCoefficient bcFunCoeff( uFun_ex );

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
		double refined_tol = tol * .01;
		// refined_tol = tol;

		NonlinearGSSolver refined_solver( mesh, order, TestEq, N_MAX_ITER, refined_tol, N_ANDERSON, 0 );
		Solution refined_soln( refined_solver.getSolutionSpace(), soln.qu );

		refined_solver.SetBCs( bcFunCoeff );
		refined_solver.Solve( refined_soln );

		auto [ err_u, err_q, err_u_star ] = refined_soln.l2_errors( uFun_ex, qFun_ex );

		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;
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
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(3, 3, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.25;
		constexpr double R_max = 1.;
		constexpr double Z_min = -1;
		constexpr double Z_max =  1;

		out( 0 ) = R_min + in( 0 )*( R_max - R_min );
		out( 1 ) = Z_min + in( 1 )*( Z_max - Z_min );
	};

	mesh->Transform( xform );

	for ( int i=0; i<N_REFINE; i++ )
		mesh->UniformRefinement();

	mesh->Finalize(true, true);

	SolovievEquilibrium TestEq( A, C, c1, c2, c3, c4 );

	constexpr double pi = boost::math::double_constants::pi;
	TSVSoln1 TestEq2( -0.5, 1.15*pi, 1.15 );

	McCarthyEquilibrium MCE( 17.8116, { 0.17795, -0.03291, 1.4934, -0.4818, -1.1759, -0.162, 0.3722, 0.07697, 1.2959, 0.5881, 1.5820, -0.009059, 2.2388, 0.4186, 1.195, -0.4265, 0.8057, -0.004804} );

	

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " Soloviev (No Flow) : " << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	TestSolution( TestEq, mesh, order );

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " Tonatiuh Manufactured Solution: " << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	TestSolution( TestEq2, mesh, order );

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << " McCarthy Equilibrium: " << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;
	TestSolution( MCE, mesh, order );
	return 0;
}


