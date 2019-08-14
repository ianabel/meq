
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

// Solov'ev from Tonatiuh
const double A = -0.52;
const double c1 = -0.001479661575325;
const double c2 = -0.366568333204813;
const double c3 =  0.002409406149732;
const double c4 = -0.023957517168316;

double SolovevRHSF( const mfem::Vector &pt )
{
	double R = pt( 0 );
	return -( R*R*( 1.0 - A ) + A )/R;
}

// Psi Solov'ev
double uFun_ex(const Vector & pt)
{
   double R(pt(0));
   double Z(pt(1));
	double psi_1 = 1.0;
	double psi_2 = R*R;
	double psi_3 = Z*Z - R*R*::log( R );
	double psi_4 = R*R*R*R - 4.0 * Z*Z * R*R;

	double psi_0 = ( 1.0 - A ) * R*R*R*R/8 + A * R * R * ::log( R ) / 2.0;

	return psi_0 + c1 * psi_1 + c2 * psi_2 + c3 * psi_3 + c4 * psi_4;
}


void qFun_ex(const Vector & pt, Vector & q)
{
   double R(pt(0));
   double Z(pt(1));

	double d_psi_0_dR = ( 1.0 - A ) * R*R*R/2 + A * R * ::log( R ) + A * R / 2.0;

	double d_psi_2_dR = 2.0 * R;
	double d_psi_3_dR = - 2.0 * R * ::log( R ) - R;
	double d_psi_4_dR = 4.0 * R*R*R - 8.0 * Z*Z * R;

	double d_psi_3_dZ = 2.0 * Z;
	double d_psi_4_dZ = - 8.0 * Z * R * R;

	q( 0 ) = ( d_psi_0_dR + c2 * d_psi_2_dR + c3 * d_psi_3_dR + c4 * d_psi_4_dR ) / R;
	q( 1 ) = ( c3 * d_psi_3_dZ + c4 * d_psi_4_dZ ) / R;

	return;
}

double bcFun( const Vector& pt )
{
	return uFun_ex( pt );
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
	Mesh *mesh = new Mesh(3, 3, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	auto xform = []( const Vector& in, Vector& out ) { 
		constexpr double R_min = 0.5;
		constexpr double R_max = 1.5;
		constexpr double Z_min = -1;
		constexpr double Z_max =  1;

		out( 0 ) = R_min + in( 0 )*( R_max - R_min );
		out( 1 ) = Z_min + in( 1 )*( Z_max - Z_min );
	};

	mesh->Transform( xform );

	int dim = mesh->Dimension();

	
	GSSolver solver( mesh, order, SolovevRHSF );

	FunctionCoefficient bcCoeff( bcFun );
	solver.SetBCs( bcCoeff );

	

	mfem::Vector qu;
	solver.Solve( qu );

	GridFunction q_variable,u_variable,u_hat_variable;
	q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
	u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );

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

		std::cout << "|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "|| q_h - q_ex || = " << err_q << "\n";
	}

	/*
	solver.ApplyAdaptiveRefinement( qu );
	qu = 0.0;
	solver.Solve( qu );
	solver.ApplyAdaptiveRefinement( qu );

	std::cout << std::endl << " Adaptive Refinement Applied " << std::endl << std::endl;

	qu = 0.0;
	solver.Solve( qu );

	q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
	u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );



	{
		double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
		double err_q    = q_variable.ComputeL2Error(qcoeff, irs);

		std::cout << "|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "|| q_h - q_ex || = " << err_q << "\n";
	}
	*/

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
	}

	delete mesh;
	return 0;
}


