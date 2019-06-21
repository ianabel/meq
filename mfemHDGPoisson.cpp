
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

#include "GSInverter.hpp"
#include "CockburnEstimator.hpp"
#include "FreeBoundary.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>


using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions

static const double pi = 3.14159265358979323844;

double zeroFun( const Vector & )
{
	return 0;
}

double uFun_ex(const Vector & pt)
{
   double x(pt(0));
   double y(pt(1));
	double a = .5;
	double b = .5;
	return x*x*::exp( -4.0 * ( x-a )*( x-a ) - 4.0*( y-b )*( y-b ) );
}

void qFun_ex(const Vector & pt, Vector & q)
{
   double x(pt(0));
   double y(pt(1));
	double a = .5;
	double b = .5;
	double e = ::exp( -4.0 * ( x-a )*( x-a ) - 4.0*( y-b )*( y-b ) );

	q( 0 ) = -e*( 2. + 8. * a * x - 8. *x *x );
	q( 1 ) = -8.*e*x*( b - y );

	return;
}

double bcFun( const Vector& pt )
{
	return uFun_ex( pt );
}


double fFun(const Vector & pt)
{
   double x(pt(0));
   double y(pt(1));

	double a = .5;
	double b = .5;
	double e = ::exp( -4.0 * ( x-a )*( x-a ) - 4.0*( y-b )*( y-b ) );
	// with nu = 1/x
	// f = -div( 1/x grad(u)) 
	//
	return -8.0*e*( 8. *a*a*x + 3.*a - 16*a*x*x + x *( -5. + 8.*b*b + 8.*x*x - 16.*b*y + 8.*y*y ) );
}

int main(int argc, char *argv[])
{
	StopWatch chrono;

	// 1. Parse command-line options.
	const char *mesh_file = "grid.mesh";
	int order = 3;
	int initial_ref_levels = 0;
	bool visualization = true;
	bool save = true;

	OptionsParser args(argc, argv);
	args.AddOption(&order, "-o", "--order",
			"Finite element order (polynomial degree).");
	args.AddOption(&initial_ref_levels, "-mr", "--mesh-refinement-levels",
			"The number of levels of uniform refinement to apply to the grid.");

	args.Parse();
	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// Mesh up (R,z) in [0.1,1] x [0,1]
	Mesh *mesh = new Mesh(10, 10, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.1;
		out( 1 ) = in( 1 );
		out( 0 ) = R_min + in( 0 )*( 1 - R_min );
	};
	mesh->Transform( xform );

	int dim = mesh->Dimension();

	for (int ii=0; ii<initial_ref_levels; ii++)
	{
		mesh->UniformRefinement();
	}

	
	GSSolver solver( mesh, order, fFun );

	ConstantCoefficient zero( 0.0 );
	FunctionCoefficient bcCoeff( bcFun );
	solver.SetBCs( zero );

	
	GridFunction q_variable,u_variable,u_hat_variable;

	mfem::Vector qu;
	solver.Solve( qu );

	q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
	u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );
	u_hat_variable.MakeRef( solver.GetMSpace(), qu, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );

	// Perform the Lackner trick to compute the boundary condition
	mfem::Vector zeroSolution = qu;
	GreensFunctionBoundaryCoefficient Lackner( mesh, solver.GetQSpace(), zeroSolution );

	mfem::Vector test_pt( 2 );
	double b_val;
	test_pt( 0 ) = 0.1;
	test_pt( 1 ) = 0.5;
	b_val = BoundaryPsi( solver.GetQSpace(), zeroSolution, test_pt );
	std::cout << "Boundary value u(0.1,0.5) is computed to be " << b_val << " but is " << uFun_ex( test_pt ) << std::endl;

	test_pt( 0 ) = 0.1;
	test_pt( 1 ) = 0.2;
	b_val = BoundaryPsi( solver.GetQSpace(), zeroSolution, test_pt );
	std::cout << "Boundary value u(0.1,0.2) is computed to be " << b_val << " but is " << uFun_ex( test_pt ) << std::endl;

	test_pt( 0 ) = 1;
	test_pt( 1 ) = 0.5;
	b_val = BoundaryPsi( solver.GetQSpace(), zeroSolution, test_pt );
	std::cout << "Boundary value u(1,0.5) is computed to be " << b_val << " but is " << uFun_ex( test_pt ) << std::endl;




	return 0;

	solver.SetBCs( Lackner );
	solver.Solve( qu );

	// 12. Compute the discretization error
	int order_quad = max(2, 3*order+2);
	const IntegrationRule *irs[Geometry::NumGeom];
	for (int i=0; i < Geometry::NumGeom; ++i)
	{
		irs[i] = &(IntRules.Get(i, order_quad));
	}
   FunctionCoefficient ucoeff(uFun_ex);
   VectorFunctionCoefficient qcoeff(dim, qFun_ex);
	{
		double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
		double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
		double err_mean  = u_variable.ComputeMeanLpError(2.0, ucoeff, irs);

		std::cout << "|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "|| mean(u_h) - mean(u_ex) || = " << err_mean << "\n";
	}

	StdFunctionCoefficient fFunCoeff( fFun );
	auto kappaF = []( const Vector& pt ) { return 1.0 / pt( 0 ); };
	StdFunctionCoefficient kappa( kappaF );
	CockburnZhangEstimator errorEstimator( q_variable, u_variable, u_hat_variable, kappa, fFunCoeff );


	ThresholdRefiner refiner( errorEstimator );
	refiner.SetTotalErrorFraction( 0.2 );
	refiner.Apply( *mesh );

	std::cout << std::endl << " Adaptive Refinement Applied " << std::endl << std::endl;

	solver.Update();
	qu = 0.0;
	solver.Solve( qu );

	q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
	u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );



	{
		double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
		double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
		double err_mean  = u_variable.ComputeMeanLpError(2.0, ucoeff, irs);

		std::cout << "|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "|| mean(u_h) - mean(u_ex) || = " << err_mean << "\n";
	}

	GridFunction ustar( solver.GetUStarSpace() );
	solver.Postprocess( ustar, qu );
	double err_ustar    = ustar.ComputeL2Error(ucoeff, irs);
	std::cout << "|| u^*_h - u_ex || = " << err_ustar << std::endl;

	// 14. Send the solution by socket to a GLVis server.
	if (visualization)
	{
		char vishost[] = "localhost";
		int  visport   = 19916;
		socketstream u_sock(vishost, visport);
		u_sock.precision(8);
		u_sock << "solution\n" << *mesh << u_variable << "window_title 'Solution u'" <<
			endl;

		socketstream q_sock(vishost, visport);
		q_sock.precision(8);
		q_sock << "solution\n" << *mesh << q_variable << "window_title 'Solution q'" <<
			endl;
	}

	// 13. Save the mesh and the solution.
	if (save)
	{
		ofstream q_variable_ofs("sol_q.gf");
		q_variable_ofs.precision(8);
		q_variable.Save(q_variable_ofs);

		ofstream u_variable_ofs("sol_u.gf");
		u_variable_ofs.precision(8);
		u_variable.Save(u_variable_ofs);

		ofstream rmesh_ofs("refined.mesh");
		rmesh_ofs.precision(8);
		mesh->Print(rmesh_ofs);
	}

	delete mesh;
	return 0;
}


