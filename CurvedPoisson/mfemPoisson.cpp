
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
double qNormalFun(const Vector & x);
void qFun_ex(const Vector & x, Vector & q);
double fFun(const Vector & x);

double tol = 1e-4;

double qNormalFun( const Vector &pt )
{
	double x( pt[ 0 ] );
	double y( pt[ 1 ] );

	mfem::Vector q( 2 );

	qFun_ex( pt, q );
	
	const double eps_tol = 1e-8;

	if ( ::fabs( x - 0 ) < eps_tol )
		return -q[ 0 ];
	if ( ::fabs( y - 0 ) < eps_tol )
		return -q[ 1 ];
	if ( ::fabs( x - 1 ) < eps_tol )
		return q[ 0 ];
	if ( ::fabs( y - 1 ) < eps_tol )
		return q[ 1 ];

	std::cerr << "being evaluated at " << x << ", " << y << std::endl;
	
	throw std::logic_error( "FTAGN!" );

	return std::nan( "" );
}

void TestSolution( std::shared_ptr<mfem::Mesh> mesh, int order )
{
	auto RHS = [ & ]( const mfem::Vector &pt, double u ) {
		return fFun( pt );
	};

	PoissonSolver Solver( mesh, order, RHS );

	FunctionCoefficient bcFunCoeff( uFun_ex );
	FunctionCoefficient bcNeumannCoeff( qNormalFun );

	Solver.SetBCs( &bcFunCoeff, 1, &bcNeumannCoeff );

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
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(3, 3, Element::Type::TRIANGLE, false, 1.0, 1.0, true );
	for ( int i=0; i < mesh->GetNBE(); i++ )
	{
		Array<int> verts( 2 );
		mesh->GetBdrElementVertices( i, verts );
		double *vert1 = mesh->GetVertex( verts[ 0 ] );
		double *vert2 = mesh->GetVertex( verts[ 1 ] );
		if ( vert1[ 0 ] == 0 && vert2[ 0 ] == 0 )
		{
			mesh->GetBdrElement( i )->SetAttribute( 1 );
		}
		else
		{
			mesh->GetBdrElement( i )->SetAttribute( 2 );
		}
	}
	mesh->SetAttributes();

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

	return sin(n*pi*xi)*sin(n*pi*yi);
}

void qFun_ex(const Vector & x, Vector & q)
{
	double xi(x(0));
	double yi(x(1));
	{
		q(0) = - diff*n*pi*cos(n*pi*xi)*sin(n*pi*yi);
		q(1) = - diff*n*pi*sin(n*pi*xi)*cos(n*pi*yi);
	}
}


double fFun(const Vector & x)
{
	double xi(x(0));
	double yi(x(1));
	return diff * ( 2.0 * n * n * pi * pi )*sin(n * pi * xi)*sin(n * pi * yi);
}
