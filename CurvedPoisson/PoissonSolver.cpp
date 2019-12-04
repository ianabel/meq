
#include "PoissonSolver.hpp"

using namespace mfem;


PoissonSolver::PoissonSolver(std::shared_ptr<mfem::Mesh> meshPtr, unsigned int order, Func RHS_f ) : 
	Order( order ),
	Dim( 2 ),
	boundary_conditions( nullptr ),
	diff_c( 1.0 ),
	tau_D( 5.0 ),
	diffusion( 1.0 ),
	RHS( RHS_f )
{
	SolutionSpace = std::make_shared<DGSpace>( meshPtr, Order );

	// Define the different forms, and initialise them with the linearised problem
	AVarf = new HDGBilinearForm( SolutionSpace->QSpace(), SolutionSpace->USpace(), SolutionSpace->MSpace() );

	AVarf->AddHDGDomainIntegrator( new HDGDomainIntegratorDiffusion( diff_c ) );
	AVarf->AddHDGFaceIntegrator( new HDGFaceIntegratorDiffusion(tau_D) );

	height = SolutionSpace->GetOffsets()[ 3 ];
	width = height;
};

void PoissonSolver::SetBCs( Coefficient& coeff )
{
	boundary_conditions = &coeff;
}

void PoissonSolver::Solve( Solution &soln )
{
	mfem::Vector qu_old( soln.qu );
	Mult( qu_old, soln.qu );
}

void PoissonSolver::Mult( const Vector& qu_in , Vector& qu_out ) const
{
	mfem::Vector rhs_F( SolutionSpace->USpace()->GetVSize() );
	qu_out.SetSize( height );
	Solution old_soln( SolutionSpace, const_cast<double*>( qu_in.GetData() ) );
	Solution new_soln( SolutionSpace, const_cast<double*>( qu_out.GetData() ) );

	// Assemble the RHS and the Schur complement
	mfem::LinearForm *fform = new mfem::LinearForm;

	this->Postprocess( old_soln );
	
	fform->AddDomainIntegrator( new NonlinearDomainLFIntegrator( old_soln.u_star_variable, RHS ) );
	fform->Update(const_cast<mfem::FiniteElementSpace* >( SolutionSpace->USpace() ), rhs_F, 0);
	fform->Assemble();

	Vector rhs_R( SolutionSpace->QSpace()->GetVSize() );

	rhs_R = 0.0;

	// To eliminate the boundary conditions we project the BC to a grid function
	// defined for the facet unknowns.
	Array<int> ess_bdr(SolutionSpace->Mesh()->bdr_attributes.Max());
	ess_bdr = 1;

	new_soln.u_hat_variable.ProjectBdrCoefficient( *boundary_conditions, ess_bdr );

	GridFunction R(SolutionSpace->QSpace(), rhs_R);
	GridFunction F(SolutionSpace->USpace(), rhs_F);
	mfem::Array<mfem::GridFunction*> F_arr( 2 );
	F_arr[ 0 ] = &R;
	F_arr[ 1 ] = &F;
	AVarf->AssembleSC(F_arr, ess_bdr, new_soln.u_hat_variable, 1.0, 1.0, 1);
	AVarf->Finalize();

	SparseMatrix* SC = AVarf->SpMatSC();
	Vector* SC_RHS = AVarf->VectorSC();
	// AVarf->VectorSC() provides -C*A^{-1} RF, the RHS for the
	// Schur complement is  L - C*A^{-1} RF, but L is zero for this case.

	// Solve the Schur complement system
	int maxIter(4000);
	double rtol(1.e-6);
	double atol(1.e-12);
	GSSmoother M(*SC);
	BiCGSTABSolver solver;
	solver.SetAbsTol(atol);
	solver.SetRelTol(rtol);
	solver.SetMaxIter(maxIter);
	solver.SetOperator(*SC);
	solver.SetPrintLevel(-1);
	solver.SetPreconditioner(M);
	solver.Mult(*SC_RHS, new_soln.u_hat_variable);

	// Reconstruct the solution u and q from the facet solution lambda
	
	mfem::Array<mfem::GridFunction*> soln_arr( 2 );
	soln_arr[ 0 ] = &new_soln.q_variable;
	soln_arr[ 1 ] = &new_soln.u_variable;

	AVarf->Reconstruct(F_arr, &new_soln.u_hat_variable, soln_arr );


	if (!solver.GetConverged())
	{
		std::cout << "Iterative method failed to converge!" << std::endl;
	}
};

void PoissonSolver::Update() {
	SolutionSpace->Update();
	AVarf->Update();

	height = SolutionSpace->GetOffsets()[ 3 ];
	width = height;
};


// Postprocessing
void PoissonSolver::Postprocess( Solution &soln ) const
{
	soln.AllocateUStar();
	GridFunction &u = soln.u_variable;
	GridFunction &q = soln.q_variable;
	GridFunction &u_postprocessed = soln.u_star_variable;

	Array<int>  vdofs;
	Vector      elmat2, shape, RHS, to_RHS, vals, uvals;
	double      RHS2;
	DenseMatrix elmat, invdfdx, dshape, dshapedxt, qvals;

	int  ndofs;
	const FiniteElement *fe_elem;
	ElementTransformation *Trans;

	FiniteElementSpace *postproc_space = SolutionSpace->UStarSpace();

	for (int i = 0; i < postproc_space->GetNE(); i++)
	{
		postproc_space->GetElementVDofs(i, vdofs);
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

		fe_elem = postproc_space->GetFE(i);
		int dim = fe_elem->GetDim();
		int spaceDim = dim;
		invdfdx.SetSize(dim, spaceDim);
		dshape.SetSize(ndofs, spaceDim);
		dshapedxt.SetSize(ndofs, spaceDim);

		Trans = SolutionSpace->Mesh()->GetElementTransformation(i);

		const mfem::IntegrationRule *ir;
		{
			int order = 3*fe_elem->GetOrder() + 4;
			ir = &IntRules.Get(fe_elem->GetGeomType(), order);
		}

		// Get the values of u_h and q_h
		u.GetValues(i, *ir, uvals);
		q.GetVectorValues(*Trans, *ir, qvals);

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

			mfem::Mult(dshape, invdfdx, dshapedxt);

			// compute the (nabla w_h, \nu \nabla u_h^*) term
			AddMult_a_AAt(w, dshapedxt, elmat);

			dshapedxt *= ip.weight ;

			Vector qval_col;
			qvals.GetColumn(j, qval_col);
			// Multiply q_h by inverse diffusion coefficient ( R in our case )
			Vector pt( 3 );
			Trans->Transform( ip, pt );
			qval_col *= 1.0/diffusion;

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

void NonlinearDomainLFIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                                mfem::ElementTransformation &Tr,
                                                mfem::Vector &elvect)
{
   int dof = el.GetDof();
	mfem::Vector shape( dof );

   elvect.SetSize(dof);
   elvect = 0.0;

   const mfem::IntegrationRule *ir;
   ir = &mfem::IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);

	mfem::Vector u_vals( ir->GetNPoints() );

	u.GetValues( Tr.ElementNo, *ir, u_vals );

	mfem::Vector point( 2 );


   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(i);

      Tr.SetIntPoint (&ip);
		Tr.Transform( ip, point );

      double val = Tr.Weight() * F( point, u_vals( i ) );

		if ( Q != nullptr )
			val *= Q->Eval( Tr, ip );

      el.CalcShape(ip, shape);
		shape *= ( ip.weight * val );

		elvect += shape;
   }
}
