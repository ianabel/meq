
/*
 * An MFEM-based nonlinear poisson solver, derived from the HDG Example Poisson Solver
 * of T. Horvath, S. Rhebergen, and A. Sivas.
 *
 * The discretization is based on the paper:
 * N.C. Nguyen, J. Peraire, B. Cockburn, An implicit high-order hybridizable
 * discontinuous Galerkin method for linear convection–diffusion equations,
 * J. Comput. Phys., 2009, 228:9, 3232--3254.
 *
 * Modifications are by I. G. Abel, University of Maryland
 */


#include "GSInverter.cpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>


using namespace std;
using namespace mfem;

// Define the analytical solution and forcing terms / boundary conditions

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

class DummyOperator : public mfem::Operator {
	public:
	DummyOperator() {};
	~DummyOperator() {};
	virtual void Mult( const Vector &, Vector & ) const {};
};

class HDGPostProcessing
{
private:
   GridFunction *q, *u;

   FiniteElementSpace *fes;

   Coefficient *diffcoeff;

protected:
   const IntegrationRule *IntRule;

public:
   HDGPostProcessing(FiniteElementSpace *f, GridFunction &_q, GridFunction &_u,
                     Coefficient &_diffcoeff)
      : q(&_q), u(&_u), fes(f), diffcoeff(&_diffcoeff) {IntRule = NULL; }

   void Postprocessing(GridFunction &u_postprocessed) ;
};

// Returns (-u^2 + f) as a suitable RHS for the poissson
// solver class. It also maintains the copy of the PoissonSolver.
class NLPoissonOperator : public Operator
{
	public:
		using RealFunc = double (*)( const mfem::Vector& );
	protected:
		FiniteElementCollection *dg_2_coll;
		FiniteElementSpace *W2_space;
		DummyOperator dummy;
	public:
		PoissonSolver pSolver;
		RealFunc RHS_ref;
		NLPoissonOperator(Mesh *meshPtr, unsigned int order, RealFunc fRHS ) : 
			RHS_ref( fRHS ),pSolver( meshPtr, order )
		{
			// More accurate space for computing the nonlinear RHS before projecting back.
			dg_2_coll = new DG_FECollection( 2*order, 2 );
			W2_space = new FiniteElementSpace( meshPtr, dg_2_coll);
			
			height = pSolver.NumCols();
			width = pSolver.NumRows();

		}

		~NLPoissonOperator()
		{
			delete W2_space;
			delete dg_2_coll;
		}

		virtual Operator& GetGradient( const Vector & ) const
		{
			return dynamic_cast< mfem::Operator& >( const_cast<DummyOperator&>( dummy ) );
		}

		void Update() {
			pSolver.Update();
			W2_space->Update(true);
			height = pSolver.NumCols();
			width = pSolver.NumRows();
		}


		virtual void Mult( Vector const& qu_in, Vector & qu_out ) const
		{
			GridFunction u_old( pSolver.GetUSpace(), static_cast<double*>( qu_in.GetData() + pSolver.bOffsets[ 1 ] ) );

			GridFunction u_fine( W2_space );
			u_fine.ProjectGridFunction( u_old );

			GridFunction u2_gf( W2_space );
			for ( int i=0; i<u2_gf.Size(); i++ )
				u2_gf[ i ] = u_fine[ i ]*u_fine[ i ];

			GridFunctionCoefficient u2_coeff( &u2_gf ); 

			// Assemble the RHS and the Schur complement
			LinearForm *fform = new LinearForm;

			qu_out.SetSize( height );
			qu_out = 0.0;
			FunctionCoefficient fcoeff( RHS_ref );
			fform->AddDomainIntegrator( new DomainLFIntegrator( fcoeff ) );
			fform->AddDomainIntegrator( new DomainLFIntegrator( u2_coeff, 3, 3 ) );
			fform->Update(pSolver.GetUSpace(), qu_out, pSolver.bOffsets[ 1 ] );
			fform->Assemble();
		}

			
};

int main(int argc, char *argv[])
{

	// 1. Parse command-line options.
	const char *mesh_file = "grid.mesh";
	int order = 3;
	int initial_ref_levels = 0;
	bool visualization = true;
	bool post = true;
	bool save = true;

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
	args.Parse();
	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// 2. Read the mesh from the given mesh file. Refine it up to the initial_ref_levels.
	Mesh *mesh = new Mesh(mesh_file, 1, 1);
	int dim = mesh->Dimension();

	for (int ii=0; ii<initial_ref_levels; ii++)
	{
		mesh->UniformRefinement();
	}
	
	Vector qu_data;
	Vector qu_refined;	

	KinSolver *nonlinearPoissonSolver = new KinSolver( KIN_PICARD, false );
	NLPoissonOperator nlPoissonEq( mesh, order, fFun );

	qu_data.SetSize( nlPoissonEq.Height() );
	qu_data = 0.0; // Initial Guess is zero

	// Quick check on coarse grid

	nonlinearPoissonSolver->SetMaxIter( 20 );
	nonlinearPoissonSolver->SetFuncNormTol( 1e-3 );
	nonlinearPoissonSolver->SetAndersonAcceleration( 3 );
	nonlinearPoissonSolver->iterative_mode = false;
	nonlinearPoissonSolver->SetSolver( dynamic_cast<mfem::Solver&>( nlPoissonEq.pSolver ) );
	nonlinearPoissonSolver->SetOperator( nlPoissonEq );
	nonlinearPoissonSolver->Mult( qu_data, qu_data );

	/*
	// Update mesh
	mesh->UniformRefinement();
	nlPoissonEq.Update();
	nlPoissonEq.pSolver.QUUpdate( qu_data, qu_refined );
	qu_data.SetSize( nlPoissonEq.Height() );
	qu_data = 0.0;

	std::cout << "Mesh Has been Refined" << std::endl;

	delete nonlinearPoissonSolver;
	nonlinearPoissonSolver = new KinSolver( KIN_FP, false );
	
	nonlinearPoissonSolver->SetMaxIter( 50 );
	nonlinearPoissonSolver->SetFuncNormTol( 1e-5 );
	nonlinearPoissonSolver->iterative_mode = true;
	nonlinearPoissonSolver->SetAndersonAcceleration( 3 );
	nonlinearPoissonSolver->SetOperator( nlPoissonEq );
	nonlinearPoissonSolver->Mult( qu_refined, qu_data );
	*/
	// Turn the vector of nodal values into Grid funcions
	Array<int> const &offsets = nlPoissonEq.pSolver.GetOffsets();
	GridFunction q_solution,u_solution;
	q_solution.MakeRef( nlPoissonEq.pSolver.GetQSpace(), qu_data, offsets[ 0 ] );
	u_solution.MakeRef( nlPoissonEq.pSolver.GetUSpace(), qu_data, offsets[ 1 ] );


	// Exact Solutions
	FunctionCoefficient ucoeff(uFun_ex);
	VectorFunctionCoefficient qcoeff(dim, qFun_ex);

	// 12. Compute the error compared to the exact solution
	int order_quad = max(2, 4*order+4);
	const IntegrationRule *irs[Geometry::NumGeom];
	for (int i=0; i < Geometry::NumGeom; ++i)
	{
		irs[i] = &(IntRules.Get(i, order_quad));
	}
	double err_u    = u_solution.ComputeL2Error(ucoeff, irs);
	double err_q    = q_solution.ComputeL2Error(qcoeff, irs);
	double err_mean  = u_solution.ComputeMeanLpError(2.0, ucoeff, irs);

	std::cout << "|| u_h - u_ex || = " << err_u << "\n";
	std::cout << "|| q_h - q_ex || = " << err_q << "\n";
	std::cout << "|| mean(u_h) - mean(u_ex) || = " << err_mean << "\n";

	// 13. Save the mesh and the solution.
	if (save)
	{
		ofstream mesh_ofs("ex_hdg.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		ofstream q_solution_ofs("sol_q.gf");
		q_solution_ofs.precision(8);
		q_solution.Save(q_solution_ofs);

		ofstream u_solution_ofs("sol_u.gf");
		u_solution_ofs.precision(8);
		u_solution.Save(u_solution_ofs);

	}

	// 14. Send the solution by socket to a GLVis server.
	if (visualization)
	{
		char vishost[] = "localhost";
		int  visport   = 19916;
		socketstream u_sock(vishost, visport);
		u_sock.precision(8);
		u_sock << "solution\n" << *mesh << u_solution << "window_title 'Solution u'" <<
			endl;

		socketstream q_sock(vishost, visport);
		q_sock.precision(8);
		q_sock << "solution\n" << *mesh << q_solution << "window_title 'Solution q'" <<
			endl;
	}

	if (post)
	{
		FiniteElementCollection *dg_coll_pstar(new DG_FECollection(order+1, 2));
		FiniteElementSpace *Vstar_space = new FiniteElementSpace(mesh, dg_coll_pstar);

		GridFunction u_post(Vstar_space);

		ConstantCoefficient diffusion( 1.0 );
		HDGPostProcessing *hdgpost(new HDGPostProcessing(Vstar_space, q_solution,
					u_solution, diffusion));

		hdgpost->Postprocessing(u_post);

		order_quad = max(2, 2*order+5);
		for (int i=0; i < Geometry::NumGeom; ++i)
		{
			irs[i] = &(IntRules.Get(i, order_quad));
		}
		double err_u_post   = u_post.ComputeL2Error(ucoeff, irs);

		std::cout << "|| u^*_h - u_ex || = " << err_u_post << "\n";

		if (save)
		{
			ofstream u_post_ofs("sol_u_star.gf");
			u_post_ofs.precision(8);
			u_post.Save(u_post_ofs);
		}

		if (visualization)
		{
			char vishost[] = "localhost";
			int  visport   = 19916;
			socketstream u_star_sock(vishost, visport);
			u_star_sock.precision(8);
			u_star_sock << "solution\n" << *mesh << u_post <<
				"window_title 'Solution u_star'" << endl;
		}
	}
	// 18. Free the used memory.
	delete mesh;

	return 0;

}

// Postprocessing
void HDGPostProcessing::Postprocessing(GridFunction &u_postprocessed)
{
	Mesh *mesh = fes->GetMesh();
	Array<int>  vdofs;
	Vector      elmat2, shape, RHS, to_RHS, vals, uvals;
	double      RHS2;
	DenseMatrix elmat, invdfdx, dshape, dshapedxt, qvals;

	int  ndofs;
	const FiniteElement *fe_elem;
	ElementTransformation *Trans;

	for (int i = 0; i < fes->GetNE(); i++)
	{
		fes->GetElementVDofs(i, vdofs);
		ndofs = vdofs.Size();
		vals.SetSize(ndofs);
		// elmat is the matrix for the -(nabla w_h, q_h) term
		elmat.SetSize(ndofs);
		// elmat 1 is the vector for the (1, u_h^*) term
		elmat2.SetSize(ndofs);
		shape.SetSize(ndofs);

		RHS.SetSize(ndofs);
		to_RHS.SetSize(ndofs);

		elmat = 0.0;
		elmat2 = 0.0;
		RHS = 0.0;
		RHS2 = 0.0;

		fe_elem = fes->GetFE(i);
		int dim = fe_elem->GetDim();
		int spaceDim = dim;
		invdfdx.SetSize(dim, spaceDim);
		dshape.SetSize(ndofs, spaceDim);
		dshapedxt.SetSize(ndofs, spaceDim);

		Trans = mesh->GetElementTransformation(i);

		const IntegrationRule *ir = IntRule;
		if (ir == NULL)
		{
			int order = 3*fe_elem->GetOrder() + 3;
			ir = &IntRules.Get(fe_elem->GetGeomType(), order);
		}

		// Get the values of u_h and q_h
		u->GetValues(i, *ir, uvals);
		q->GetVectorValues(*Trans, *ir, qvals);

		for (int j = 0; j < ir->GetNPoints(); j++)
		{
			const IntegrationPoint &ip = ir->IntPoint(j);

			fe_elem->CalcDShape(ip, dshape);
			fe_elem->CalcShape(ip, shape);

			Trans->SetIntPoint(&ip);
			// Compute invdfdx = / adj(J),       if J is square
			//               \ adj(J^t.J).J^t, otherwise
			CalcAdjugate(Trans->Jacobian(), invdfdx);
			double w = Trans->Weight();
			w = ip.weight / w;
			w *= diffcoeff->Eval(*Trans, ip);
			Mult(dshape, invdfdx, dshapedxt);

			// compute the (nabla w_h, \nu \nabla u_h^*) term
			AddMult_a_AAt(w, dshapedxt, elmat);

			dshapedxt *= ip.weight ;

			Vector qval_col;
			qvals.GetColumn(j, qval_col);

			// compute (nabla w_h, q_h)
			dshapedxt.Mult(qval_col, to_RHS);

			// subtract it from the rhs
			RHS -= to_RHS;

			// compute (1, u_h^*)
			shape *= (Trans->Weight() * ip.weight);
			elmat2 += shape;

			// compute (1, u_h)
			double rhs_weight = (Trans->Weight() * ip.weight);
			RHS2  += (uvals(j)*rhs_weight);

		}

		// changing the last row and the last entry
		for (int j = 0; j < ndofs; j++)
		{
			elmat(ndofs-1,j) = elmat2(j);
		}
		RHS(ndofs-1) = RHS2;

		// solve the local problem
		elmat.Invert();
		elmat.Mult(RHS, vals);
		u_postprocessed.SetSubVector(vdofs, vals);

	}
}
