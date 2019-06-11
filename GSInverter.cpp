#include "GSInverter.hpp"

using namespace mfem;

GSInverter::GSInverter(Mesh *meshPtr, unsigned int order, GSInverter::RealFunc fRHS ) : 
	mesh( meshPtr ),
	Order( order ),
	Dim( 2 ),
	tau_D( 5.0 ),
	RHS_ref( fRHS )
{
	// Define finite element collections and spaces on the mesh.
	dg_coll   = new DG_FECollection(Order, Dim);
	face_coll = new DG_Interface_FECollection(Order, Dim);

	// Finite element spaces:
	// V_space is the scalar DG space on elements for u_h
	// W_space is the vector valued DG space on elements for q_h
	// M_space is the DG space on faces for lambda_h
	V_space = new FiniteElementSpace(mesh, dg_coll, Dim);
	W_space = new FiniteElementSpace(mesh, dg_coll);
	M_space = new FiniteElementSpace(mesh, face_coll);
	dimV = V_space->GetVSize();
	dimW = W_space->GetVSize();
	dimM = M_space->GetVSize();

	// Define the different forms, and initialise them with the linearised problem
	AVarf = new HDGBilinearForm3(V_space, W_space, M_space);

	AVarf->AddHDGDomainIntegrator(new HDGDomainIntegratorGS());
	AVarf->AddHDGBdrIntegrator(new HDGFaceIntegratorGS(tau_D));


	bOffsets.SetSize( 3 );
	bOffsets[ 0 ] = 0;
	bOffsets[ 1 ] = dimV;
	bOffsets[ 2 ] = dimV + dimW;
	height = dimV + dimW;
	width = dimV + dimW;
};


/* 
 * Prolongs a vector containing q & u from the old mesh to the 
 * new. This will not handle increasing polynomial order.
 */
void GSInverter::QUUpdate( Vector const& qu_old, Vector &qu_new ) const
{
	const Operator* U_update = W_space->GetUpdateOperator();
	const Operator* Q_update = V_space->GetUpdateOperator();
	int U_old_dim = U_update->Width();
	int U_new_dim = U_update->Height();
	int Q_old_dim = Q_update->Width();
	int Q_new_dim = Q_update->Height();
	if ( U_old_dim + Q_old_dim != qu_old.Size() )
		throw new std::logic_error( "Stop trying to multiply badgers by goats!" );
	qu_new.SetSize( U_new_dim + Q_new_dim );
	qu_new = 0.;

	Array<int> oldOffsets; oldOffsets.SetSize( 3 );
	oldOffsets[ 0 ] = 0; oldOffsets[ 1 ] = Q_old_dim; oldOffsets[ 2 ] = U_old_dim + Q_old_dim;
	BlockVector QU_old_blk( qu_old.GetData(), oldOffsets );
	BlockVector QU_new_blk( qu_new.GetData(), bOffsets );

	Q_update->Mult( QU_old_blk.GetBlock( 0 ), QU_new_blk.GetBlock( 0 ) );
	U_update->Mult( QU_old_blk.GetBlock( 1 ), QU_new_blk.GetBlock( 1 ) );

}

void GSInverter::SetBCs( Coefficient& coeff )
{
	BoundaryConditions.SetSpace( M_space );
	Array<int> boundary( mesh->bdr_attributes.Max() );
	boundary = 1; // The entire boundary is dirichlet
	BoundaryConditions.ProjectBdrCoefficient( coeff, boundary );
}

// Actually solve the problem:
// which in this case doesn't depend on the input vector
// and store in the Vector y
void GSInverter::Mult( const Vector& qu_in , Vector& qu_out ) const
{
	StopWatch chrono;
	Vector rhs_R(dimV);
	Vector rhs_F(dimW);
	Vector V_aux(dimV);
	Vector W_aux(dimW);

	V_aux = 0.0;
	W_aux = 0.0;
	rhs_R = 0.0;

	// To eliminate the boundary conditions we project the BC to a grid function
	// defined for the facet unknowns.
	GridFunction lambda_variable( BoundaryConditions );


	Array<int> ess_bdr(mesh->bdr_attributes.Max());
	ess_bdr = 1;

	// Assemble the RHS and the Schur complement
	LinearForm *fform = new LinearForm;
	StdFunctionCoefficient fcoeff( RHS_ref );
	fform->AddDomainIntegrator( new DomainLFIntegrator( fcoeff ) );
	fform->Update(W_space, rhs_F, 0);
	fform->Assemble();

	GridFunction R(V_space, rhs_R);
	GridFunction F(W_space, rhs_F);
	AVarf->AssembleSC(R, F, ess_bdr, lambda_variable, 1, 1);
	AVarf->Finalize();

	SparseMatrix* SC = AVarf->SpMatSC();
	Vector* SC_RHS = AVarf->VectorSC();
	// AVarf->VectorSC() provides -C*A^{-1} RF, the RHS for the
	// Schur complement is  L - C*A^{-1} RF, but L is zero for this case.

	// Solve the Schur complement system
	int maxIter(4000);
	double rtol(1.e-13);
	double atol(0.0);
	GSSmoother M(*SC);
	BiCGSTABSolver solver;
	solver.SetAbsTol(atol);
	solver.SetRelTol(rtol);
	solver.SetMaxIter(maxIter);
	solver.SetOperator(*SC);
	solver.SetPrintLevel(-1);
	solver.SetPreconditioner(M);
	chrono.Clear();
	chrono.Start();
	solver.Mult(*SC_RHS, lambda_variable);
	chrono.Stop();

	if (!solver.GetConverged())
	{
		std::cout << "Iterative method failed to converge!" << std::endl;
	}
	else
	{
		std::cout << "Linear Solve took " << chrono.RealTime() << " seconds." << std::endl;
	}

	// Reconstruct the solution u and q from the facet solution lambda
	GridFunction q_variable,u_variable;
	qu_out.SetSize( dimV + dimW );
	q_variable.MakeRef( V_space, qu_out, 0 );
	u_variable.MakeRef( W_space, qu_out, dimV );
	AVarf->Reconstruct(&R, &F, lambda_variable, &q_variable, &u_variable);
};

void GSInverter::Update() {
	V_space->Update(true);
	W_space->Update(true);
	M_space->Update(false);

	AVarf->Update();

	dimV = V_space->GetVSize();
	dimW = W_space->GetVSize();
	dimM = M_space->GetVSize();


	bOffsets.SetSize( 3 );
	bOffsets[ 0 ] = 0;
	bOffsets[ 1 ] = dimV;
	bOffsets[ 2 ] = dimV + dimW;
	height = dimV + dimW;
	width = dimV + dimW;
};

