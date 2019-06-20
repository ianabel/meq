
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
	return ::exp( 2*x*y )*::sin( pi*x )*::cos( pi*y );
}

void qFun_ex(const Vector & pt, Vector & q)
{
   double x(pt(0));
   double y(pt(1));

	q( 0 ) = -::exp( 2*x*y )*::cos( pi*y )*(  pi*::cos( pi*x ) + 2*y*::sin( pi*x ) ) / x;
	q( 1 ) = -::exp( 2*x*y )*::sin( pi*x )*( -pi*::sin( pi*y ) + 2*x*::cos( pi*y ) ) / x;

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
	// with nu = 1/x
	double t = ( pi*pi*x - 2*x*x*x + y - 2*x*y*y )*::cos( pi*y ) + 2*pi*x*x*::sin( pi*y );
	return ( ::exp( 2*x*y ) / ( x*x ) )*( 
			  pi*( 1.0 - 4.0 *  x  * y ) * ::cos( pi*x ) * ::cos( pi*y ) + ( 2 * t ) * ::sin( pi*x ) 
			);

	/*
	// with nu = 1
	double t = ( pi*pi - 2*x*x - 2*y*y )*::cos( pi*y ) + 2*pi*x*::sin( pi*y );
	return ( 2.0 * ::exp( 2*x*y ) )*( -2.0 * pi * y * ::cos( pi*x ) * ::cos( pi*y ) + t * ::sin( pi*x ) );
	*/

}

int main(int argc, char *argv[])
{
	StopWatch chrono;

	// 1. Parse command-line options.
	const char *mesh_file = "grid.mesh";
	int order = 3;
	int initial_ref_levels = 0;
	bool visualization = true;
	bool post = true;
	bool save = true;
	double memA = 0.0;
	double memB = 0.0;

	OptionsParser args(argc, argv);
	args.AddOption(&mesh_file, "-m", "--mesh",
			"Mesh file to use.");
	args.AddOption(&order, "-o", "--order",
			"Finite element order (polynomial degree).");
	args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
			"--no-visualization",
			"Enable or disable GLVis visualization.");
	args.AddOption(&post, "-post", "--postprocessing",
			"-no-post", "--no-postprocessing",
			"Enable or disable postprocessing.");
	args.AddOption(&save, "-save", "--save-files", "-no-save",
			"--no-save-files",
			"Enable or disable file saving.");
	args.AddOption(&initial_ref_levels, "-mr", "--mesh-refinement-levels",
			"The number of levels of uniform refinement to apply to the grid.");
	args.AddOption(&memA, "-memA", "--memoryA",
			"Storage of A.");
	args.AddOption(&memB, "-memB", "--memoryB",
			"Storage of B.");

	args.Parse();
	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// memA, memB \in [0,1], memB <= memA
	if (memB > memA)
	{
		std::cout << "memB cannot be more than memA. Resetting to be equal" << std::endl
			<< std::flush;
		memA = memB;
	}
	if (memA > 1.0)
	{
		std::cout << "memA cannot be more than 1. Resetting to 1" << std::endl <<
			std::flush;
		memA = 1.0;
	}
	else if (memA < 0.0)
	{
		std::cout << "memA cannot be less than 0. Resetting to 0." << std::endl <<
			std::flush;
		memA = 0.0;
	}
	if (memB > 1.0)
	{
		std::cout << "memB cannot be more than 1. Resetting to 1" << std::endl <<
			std::flush;
		memB = 1.0;
	}
	else if (memB < 0.0)
	{
		std::cout << "memB cannot be less than 0. Resetting to 0." << std::endl <<
			std::flush;
		memB = 0.0;
	}

	// Mesh up [0,1] x [0,1]
	Mesh *mesh = new Mesh(8,8,Element::Type::TRIANGLE, false, 1.0, 1.0, true );
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

	FunctionCoefficient bcFunCoeff( bcFun );
	solver.SetBCs( bcFunCoeff );

	
	GridFunction q_variable,u_variable,u_hat_variable;

	mfem::Vector qu;
	solver.Solve( qu );

	q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
	u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );
	u_hat_variable.MakeRef( solver.GetMSpace(), qu, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );

	StdFunctionCoefficient fFunCoeff( fFun );
	auto kappaF = []( const Vector& pt ) { return 1.0 / pt( 0 ); };
	StdFunctionCoefficient kappa( kappaF );
	CockburnZhangEstimator errorEstimator( q_variable, u_variable, u_hat_variable, kappa, fFunCoeff );

	ThresholdRefiner refiner( errorEstimator );
	refiner.SetTotalErrorFraction( 0.7 );
	refiner.Apply( *mesh );



	// 12. Compute the discretization error
	int order_quad = max(2, 3*order+2);
	const IntegrationRule *irs[Geometry::NumGeom];
	for (int i=0; i < Geometry::NumGeom; ++i)
	{
		irs[i] = &(IntRules.Get(i, order_quad));
	}
   FunctionCoefficient ucoeff(uFun_ex);
   VectorFunctionCoefficient qcoeff(dim, qFun_ex);
	double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
	double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
	double err_mean  = u_variable.ComputeMeanLpError(2.0, ucoeff, irs);

	std::cout << "|| u_h - u_ex || = " << err_u << "\n";
	std::cout << "|| q_h - q_ex || = " << err_q << "\n";
	std::cout << "|| mean(u_h) - mean(u_ex) || = " << err_mean << "\n";

	// 13. Save the mesh and the solution.
	if (save)
	{
		ofstream mesh_ofs("refined.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		ofstream q_variable_ofs("sol_q.gf");
		q_variable_ofs.precision(8);
		q_variable.Save(q_variable_ofs);

		ofstream u_variable_ofs("sol_u.gf");
		u_variable_ofs.precision(8);
		u_variable.Save(u_variable_ofs);

	}

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

	GridFunction ustar( solver.GetUStarSpace() );
	solver.Postprocess( ustar, qu );
	double err_ustar    = ustar.ComputeL2Error(ucoeff, irs);
	std::cout << "|| u^*_h - u_ex || = " << err_ustar << std::endl;


	delete mesh;
	return 0;
}


