
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
#include "NLRHSIntegrator.hpp"
#include "CockburnEstimator.hpp"
#include "StdFnCoeffs.hpp"
#include "SolovievEquilibrium.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>


using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions

static const double pi = 3.14159265358979323844;

// Solov'ev from Tonatiuh
const double A = -0.52;
const double C = 1.52;
const double c1 = -0.001479661575325;
const double c2 = -0.366568333204813;
const double c3 =  0.002409406149732;
const double c4 = -0.023957517168316;

class NLGSSolver : public mfem::Operator
{
	public:
		using Func = std::function< double( double ) >;
	protected:
		GSInverter solver;
		Func Pprime,FFprime;
	public:
		NLGSSolver(mfem::Mesh *meshPtr, unsigned int order, Func Pprime_orig, Func FFprime_orig )
			: solver( meshPtr, order ), Pprime( Pprime_orig ), FFprime( FFprime_orig )
		{
			height = solver.NumRows();
			width = solver.NumRows();
		};

		void SetBCs( mfem::Coefficient& coeff ) {
			solver.SetBCs( coeff );
		};

		void Update() { 
			solver.Update(); 
			height = solver.NumRows();
			width = solver.NumRows();
		};

		virtual void Mult( mfem::Vector const& qu_in, mfem::Vector & qu_out ) const
		{
			mfem::Vector rhs_F( solver.NumCols() );
			qu_out.SetSize( solver.NumRows() );
			// Assemble the RHS and the Schur complement
			mfem::LinearForm *fform = new mfem::LinearForm;

			mfem::GridFunction u;
			u.MakeRef( const_cast<mfem::FiniteElementSpace* >( solver.GetUSpace() ), static_cast<double*>( qu_in.GetData() + solver.GetQSpace()->GetVSize() ) );

			auto PlasmaRHS = [ this ]( const mfem::Vector &pt, double psi ){
				double R = pt( 0 );
				return -( R*R*Pprime( psi ) + FFprime( psi ) )/R;
			};

			fform->AddDomainIntegrator( new NonlinearDomainLFIntegrator( u, PlasmaRHS, 4, 2 ) );
			fform->Update(const_cast<mfem::FiniteElementSpace* >( solver.GetUSpace() ), rhs_F, 0);
			fform->Assemble();

			solver.Mult( rhs_F, qu_out );
			delete fform;
		};

		void Postprocess( mfem::GridFunction &u_out, mfem::Vector & qu_in )
		{
			solver.Postprocess( u_out, qu_in );
		};

		mfem::FiniteElementSpace const * GetMSpace() const { return solver.GetMSpace(); };
		mfem::FiniteElementSpace const * GetQSpace() const { return solver.GetQSpace(); };
		mfem::FiniteElementSpace const * GetUSpace() const { return solver.GetUSpace(); };
		mfem::FiniteElementSpace const * GetUStarSpace() const { return solver.GetUStarSpace(); };
		mfem::FiniteElementSpace * GetMSpace() { return solver.GetMSpace(); };
		mfem::FiniteElementSpace * GetQSpace() { return solver.GetQSpace(); };
		mfem::FiniteElementSpace * GetUSpace() { return solver.GetUSpace(); };
		mfem::FiniteElementSpace * GetUStarSpace() { return solver.GetUStarSpace(); };

		void Prolong( mfem::Vector const& qu_old, mfem::Vector & qu_new ) const
		{
			solver.Prolong( qu_old, qu_new );
		};

		void ApplyAdaptiveRefinement(  mfem::Vector & soln_vector )
		{
			mfem::GridFunction q_variable,u_variable,u_hat_variable;

			q_variable.MakeRef( solver.GetQSpace(), soln_vector, 0 );
			u_variable.MakeRef( solver.GetUSpace(), soln_vector, solver.GetQSpace()->GetVSize() );
			u_hat_variable.MakeRef( solver.GetMSpace(), soln_vector, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );


			mfem::GridFunction u_star( solver.GetUStarSpace() );
			solver.Postprocess( u_star, soln_vector );

			auto PlasmaRHS = [ this ]( const mfem::Vector &pt, double psi ){
				double R = pt( 0 );
				return -( R*R*Pprime( psi ) + FFprime( psi ) )/R;
			};
			mfem::GradShafranovEstimator errorEstimator( q_variable, u_star, u_hat_variable, PlasmaRHS );
			mfem::ThresholdRefiner refiner( errorEstimator );

			refiner.SetTotalErrorFraction( 0.2 );
			refiner.Apply( *( solver.GetMesh() ) );
			Update();
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


	SolovievEquilibrium TestEq( A, C, c1, c2, c3, c4 );

	std::function<double( double )> Pprime  = std::bind(  &SolovievEquilibrium::Pprime, &TestEq, std::placeholders::_1 );
	std::function<double( double )> FFprime = std::bind( &SolovievEquilibrium::FFprime, &TestEq, std::placeholders::_1 );
	NLGSSolver solver( mesh, order, Pprime, FFprime );

	std::function<double( const mfem::Vector & )> uFun_ex = std::bind( &SolovievEquilibrium::psi, &TestEq, std::placeholders::_1 );
	std::function<void( const mfem::Vector &, mfem::Vector & )> qFun_ex = std::bind( &SolovievEquilibrium::q, &TestEq, std::placeholders::_1, std::placeholders::_2 );

	StdFunctionCoefficient bcFunCoeff( uFun_ex );
	solver.SetBCs( bcFunCoeff );

	GridFunction q_variable,u_variable,u_hat_variable;

	mfem::Vector qu( solver.NumRows() );

	qu = 0.0;

	int i_amr = 0;
	KinSolver *nonlinearPoissonSolver = nullptr;
	mfem::Vector qu_refined;
	for ( ; i_amr< N_AMR; i_amr++ )
	{
		nonlinearPoissonSolver = new KinSolver( KIN_FP, false );
		nonlinearPoissonSolver->SetMaxIter( 1000 );
		nonlinearPoissonSolver->SetFuncNormTol( 1e-4 );
		nonlinearPoissonSolver->iterative_mode = false;
		nonlinearPoissonSolver->SetOperator( solver );
		qu.SetSize( solver.Height() );
		nonlinearPoissonSolver->Mult( qu, qu );

		q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
		u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );
		u_hat_variable.MakeRef( solver.GetMSpace(), qu, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );

		std::cout << "After " << i_amr << " levels of refinement:" << std::endl;

		int order_quad = max(2, 2*order+2);
		const IntegrationRule *irs[Geometry::NumGeom];
		for (int i=0; i < Geometry::NumGeom; ++i)
		{
			irs[i] = &(IntRules.Get(i, order_quad));
		}
		StdFunctionCoefficient ucoeff(uFun_ex);
		VectorStdFunctionCoefficient qcoeff(dim, qFun_ex);
		double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
		double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
		mfem::GridFunction u_star( solver.GetUStarSpace() );
		solver.Postprocess( u_star, qu );
		double err_u_star    = u_star.ComputeL2Error(ucoeff, irs);


		std::cout << "\t|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "\t|| q_h - q_ex || = " << err_q << "\n";
		std::cout << "\t|| u*_h - u_ex || = " << err_u_star << "\n";
		std::cout << std::endl;

		solver.ApplyAdaptiveRefinement( qu );

		std::cerr << "Refined Mesh" << std::endl;

		delete nonlinearPoissonSolver;
	}

	{
		nonlinearPoissonSolver = new KinSolver( KIN_FP, false );
		nonlinearPoissonSolver->SetMaxIter( 1000 );
		nonlinearPoissonSolver->SetFuncNormTol( 1e-4 );
		nonlinearPoissonSolver->iterative_mode = false;
		nonlinearPoissonSolver->SetOperator( solver );
		qu.SetSize( solver.Height() );
		nonlinearPoissonSolver->Mult( qu, qu );


		q_variable.MakeRef( solver.GetQSpace(), qu, 0 );
		u_variable.MakeRef( solver.GetUSpace(), qu, solver.GetQSpace()->GetVSize() );
		u_hat_variable.MakeRef( solver.GetMSpace(), qu, solver.GetQSpace()->GetVSize() + solver.GetUSpace()->GetVSize() );

		std::cout << "After " << i_amr << " levels of refinement:" << std::endl;

		int order_quad = max(2, 2*order+2);
		const IntegrationRule *irs[Geometry::NumGeom];
		for (int i=0; i < Geometry::NumGeom; ++i)
		{
			irs[i] = &(IntRules.Get(i, order_quad));
		}
		StdFunctionCoefficient ucoeff(uFun_ex);
		VectorStdFunctionCoefficient qcoeff(dim, qFun_ex);
		double err_u    = u_variable.ComputeL2Error(ucoeff, irs);
		double err_q    = q_variable.ComputeL2Error(qcoeff, irs);
		mfem::GridFunction u_star( solver.GetUStarSpace() );
		solver.Postprocess( u_star, qu );
		double err_u_star    = u_star.ComputeL2Error(ucoeff, irs);


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

		ofstream q_variable_ofs("sol_q.gf");
		q_variable_ofs.precision(8);
		q_variable.Save(q_variable_ofs);

		ofstream u_variable_ofs("sol_u.gf");
		u_variable_ofs.precision(8);
		u_variable.Save(u_variable_ofs);

	}

	// 14. Send the solution by socket to a GLVis server.
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


	delete mesh;
	return 0;
}


