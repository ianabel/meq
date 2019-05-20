
using namespace mfem;

GSInverter(Mesh *meshPtr, unsigned int order, RealFunc fRHS ) : 
			mesh( meshPtr ),
			Order( order ),
			Dim( 2 ),
			diffusion( 1.0 ),
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
		AVarf->AddHDGDomainIntegrator(new HDGDomainIntegratorDiffusion(diffusion));
		AVarf->AddHDGBdrIntegrator(new HDGFaceIntegratorDiffusion(tau_D));
		

		bOffsets.SetSize( 3 );
		bOffsets[ 0 ] = 0;
		bOffsets[ 1 ] = dimV;
		bOffsets[ 2 ] = dimV + dimW;
		height = dimV + dimW;
		width = dimV + dimW;
	};

	FiniteElementSpace* GetQSpace() { return V_space; };
	FiniteElementSpace* GetUSpace() { return W_space; };

	/* 
	 * Prolongs a vector containing q & u from the old mesh to the 
	 * new. This will not handle increasing polynomial order.
	 */
	void QUUpdate( Vector const& qu_old, Vector &qu_new ) const
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


	Array<int> const & GetOffsets() { return bOffsets; };

	// Actually solve the problem:
	// which in this case doesn't depend on the input vector
	// and store in the Vector y
	virtual void Mult( const Vector& qu_in , Vector& qu_out ) const
	{
		StopWatch chrono;
		Vector rhs_R(dimV);
		Vector rhs_F(dimW);
		Vector V_aux(dimV);
		Vector W_aux(dimW);

		V_aux = 0.0;
		W_aux = 0.0;
		rhs_R = 0.0;

		GridFunction lambda_variable( M_space );

		// To eliminate the boundary conditions we project the BC to a grid function
		// defined for the facet unknowns.
		FunctionCoefficient lambda_coeff(zeroFun);
		lambda_variable.ProjectCoefficientSkeletonDG(lambda_coeff);

		Array<int> ess_bdr(mesh->bdr_attributes.Max());
		ess_bdr = 1;

		Vector u_squared( dimW );
		for ( int i=0; i< dimW; i++ )
		{
			double u_val = qu_in[ i + bOffsets[ 1 ] ];
			u_squared[ i ] = 2*u_val * u_val;
		}

		GridFunction u2_gf;
		u2_gf.MakeRef( W_space, u_squared, 0 );
		GridFunctionCoefficient u2_coeff( &u2_gf ); 

		// Assemble the RHS and the Schur complement
		LinearForm *fform = new LinearForm;
		FunctionCoefficient fcoeff( RHS_ref );
		fform->AddDomainIntegrator( new DomainLFIntegrator( fcoeff ) );
		fform->AddDomainIntegrator( new DomainLFIntegrator( u2_coeff ) );
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
		q_variable.MakeRef( V_space, qu_out, 0 );
		u_variable.MakeRef( W_space, qu_out, dimV );
		AVarf->Reconstruct(&R, &F, lambda_variable, &q_variable, &u_variable);
	};

	void Update() {
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

	~NLPoissonOperator() 
	{
		delete V_space;
		delete W_space;
		delete M_space;
		delete AVarf;
		delete dg_coll;
		delete face_coll;
	};
		
};
