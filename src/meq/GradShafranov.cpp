#include "GradShafranov.hpp"

#include <stdexcept>

namespace meq
{

	ConstantStabilization::ConstantStabilization( double tauIn )
		: tauValue( tauIn )
	{
	}

	bool ConstantStabilization::IsConstant() const
	{
		return true;
	}

	mfem::real_t ConstantStabilization::Eval( mfem::real_t, mfem::real_t, mfem::real_t,
	                                          mfem::real_t,
	                                          mfem::ElementTransformation & ) const
	{
		return tauValue;
	}

	double ConstantStabilization::tau() const
	{
		return tauValue;
	}

	GradShafranovSolver::GradShafranovSolver( mfem::Mesh &meshIn, int orderIn, double tauIn )
		: mesh( meshIn ),
		  orderValue( orderIn ),
		  stabilization( tauIn ),
		  radius( []( mfem::Vector const &x ) { return x( 0 ); } ),
		  negativeInverseRadius( []( mfem::Vector const &x ) { return -1.0/x( 0 ); } ),
		  source( nullptr ),
		  boundaryData( nullptr ),
		  assembled( false )
	{
		if ( orderValue < 0 )
			throw std::invalid_argument( "meq::GradShafranovSolver: the polynomial order must not be negative" );
		if ( mesh.Dimension() != 2 )
			throw std::invalid_argument( "meq::GradShafranovSolver: the mesh must be two dimensional ( r, z )" );

		int const dim = mesh.Dimension();

		// The closed (Gauss-Lobatto) basis, as miniapps/hdg/convdiff.cpp puts it,
		// "as it is customary for HDG to match trace DOFs". All three spaces carry
		// the same degree; hybridization is what makes that legal.
		fluxColl = std::make_unique<mfem::L2_FECollection>( orderValue, dim,
		                                                    mfem::BasisType::GaussLobatto );
		potentialColl = std::make_unique<mfem::L2_FECollection>( orderValue, dim,
		                                                         mfem::BasisType::GaussLobatto );
		traceColl = std::make_unique<mfem::DG_Interface_FECollection>( orderValue, dim );

		fluxFes = std::make_unique<mfem::FiniteElementSpace>( &mesh, fluxColl.get(), dim );
		potentialFes = std::make_unique<mfem::FiniteElementSpace>( &mesh, potentialColl.get() );
		traceFes = std::make_unique<mfem::FiniteElementSpace>( &mesh, traceColl.get() );

		darcy = std::make_unique<mfem::DarcyForm>( fluxFes.get(), potentialFes.get() );

		dirichletMarker.SetSize( mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0 );
		dirichletMarker = 1;

		// ( r q, v ). DarcyForm's flux mass form holds the INVERSE of the diffusion
		// coefficient -- convdiff puts 1/k there -- and the coefficient here is
		// 1/r, so this is r. Measured rather than assumed: putting 1/r here instead
		// still converges, to a different function, with the L2 error against the
		// exact Solov'ev solution flat at 1.9e-2 through four refinements.
		mfem::BilinearForm *fluxMass = darcy->GetFluxMassForm();
		fluxMass->AddDomainIntegrator( new mfem::VectorMassIntegrator( radius ) );

		// < tau( psi_h - psihat_h ), w > on every face of every element, interior
		// and boundary alike. The coefficient handed to HDGDiffusionIntegrator is
		// dead weight once SetStabilization() is called -- the hook replaces the
		// built-in expression entirely -- but it is the diffusion coefficient the
		// integrator is documented to take, so it is the honest thing to pass.
		mfem::BilinearForm *potentialMass = darcy->GetPotentialMassForm();
		{
			auto *interior = new mfem::HDGDiffusionIntegrator( radius, stabilization.tau() );
			auto *boundary = new mfem::HDGDiffusionIntegrator( radius, stabilization.tau() );
			interior->SetStabilization( stabilization );
			boundary->SetStabilization( stabilization );
			potentialMass->AddInteriorFaceIntegrator( interior );
			potentialMass->AddBdrFaceIntegrator( boundary, dirichletMarker );
		}

		// ( div_bar q, w ) and, by transposition, ( psi, div_bar v ).
		//
		// The two face integrators need a word, because under hybridization their
		// numerical values are never used. DarcyForm::Assemble() builds B from
		// ComputeElementMatrix(), which sums domain integrators only; the flux-trace
		// coupling < psihat, v.n > comes instead from the transpose of the
		// constraint operator below. What the boundary face integrator does do is
		// carry a marker: EnableHybridization() reads B's boundary-face markers and
		// registers a boundary flux constraint on exactly those attributes. Drop it
		// and the Dirichlet boundary faces get no constraint at all -- measured, the
		// L2 error then sits at 1.5e-1 and does not move under refinement. Changing
		// its coefficient from -2 to anything else, or removing the interior one,
		// changes the answer in not one digit. The values are convdiff's, kept so
		// that the form is still right if hybridization is ever switched off.
		mfem::MixedBilinearForm *fluxDiv = darcy->GetFluxDivForm();
		fluxDiv->AddDomainIntegrator( new mfem::VectorDivergenceIntegrator() );
		fluxDiv->AddInteriorFaceIntegrator(
			new mfem::TransposeIntegrator( new mfem::DGNormalTraceIntegrator( -1.0 ) ) );
		fluxDiv->AddBdrFaceIntegrator(
			new mfem::TransposeIntegrator( new mfem::DGNormalTraceIntegrator( -2.0 ) ),
			dirichletMarker );

		// < qhat_h.n, mu > = 0. This must come after every AddIntegrator above:
		// EnableHybridization() reaches into the potential mass and flux divergence
		// forms as it runs and takes what is there at that moment.
		darcy->EnableHybridization( traceFes.get(), new mfem::NormalTraceJumpIntegrator(),
		                            essentialFluxTdofs );
		darcy->GetHybridization()->SetEssentialBC( dirichletMarker );

		blockOffsets.SetSize( 4 );
		blockOffsets[ 0 ] = 0;
		blockOffsets[ 1 ] = fluxFes->GetVSize();
		blockOffsets[ 2 ] = potentialFes->GetVSize();
		blockOffsets[ 3 ] = traceFes->GetVSize();
		blockOffsets.PartialSum();

		solution.Update( blockOffsets );
		rhs.Update( blockOffsets );
		solution = 0.0;
		rhs = 0.0;

		darcyFlux.MakeRef( fluxFes.get(), solution.GetBlock( 0 ), 0 );
		potentialGf.MakeRef( potentialFes.get(), solution.GetBlock( 1 ), 0 );
		traceGf.MakeRef( traceFes.get(), solution.GetBlock( 2 ), 0 );

		fluxGf.SetSpace( fluxFes.get() );
		fluxGf = 0.0;
	}

	void GradShafranovSolver::setSource( mfem::Coefficient &fIn )
	{
		source = &fIn;

		// The potential right hand side is -( F/r, w ), and both signs in that are
		// real. The 1/r is the equation's: the right hand side is F/r, not F. The
		// minus is DarcyForm's: constructed with its default bsymmetrize = true it
		// assembles the second block row as -B q - Mp psi = bp, so the datum handed
		// to it is the negative of the source of div q. Measured: with +F/r the L2
		// error against the exact solution is flat at 7.3e-2 through four
		// refinements, with -F/r it converges at k+1.
		potentialRhsCoeff = std::make_unique<mfem::ProductCoefficient>( negativeInverseRadius,
		                                                                *source );
	}

	void GradShafranovSolver::setBoundaryData( mfem::Coefficient &boundaryIn )
	{
		boundaryData = &boundaryIn;
	}

	void GradShafranovSolver::assembleForms()
	{
		if ( assembled )
			return;
		darcy->Assemble();
		assembled = true;
	}

	void GradShafranovSolver::solve()
	{
		if ( !source )
			throw std::logic_error( "meq::GradShafranovSolver::solve: no source has been set" );
		if ( !boundaryData )
			throw std::logic_error( "meq::GradShafranovSolver::solve: no boundary data has been set" );

		assembleForms();

		solution = 0.0;
		rhs = 0.0;

		// The Dirichlet datum lives on the trace, and only on the trace: the flux
		// right hand side stays zero because the essential trace condition, not a
		// boundary linear form, is what imposes psi = g_D here.
		traceGf.ProjectBdrCoefficient( *boundaryData, dirichletMarker );

		mfem::LinearForm potentialRhs;
		potentialRhs.Update( potentialFes.get(), rhs.GetBlock( 1 ), 0 );
		potentialRhs.AddDomainIntegrator( new mfem::DomainLFIntegrator( *potentialRhsCoeff ) );
		potentialRhs.Assemble();

		// X and B alias the trace block of the solution and right hand side. That
		// aliasing is not cosmetic: FormLinearSystem() only calls
		// EliminateTraceTrueDofsInRHS() -- the step that actually moves the
		// essential trace values into the right hand side -- when X arrives already
		// sized to the reduced system and copy_interior is set.
		mfem::Vector traceX;
		mfem::Vector traceB;
		traceX.MakeRef( solution, blockOffsets[ 2 ], traceFes->GetVSize() );
		traceB.MakeRef( rhs, blockOffsets[ 2 ], traceFes->GetVSize() );

		mfem::OperatorHandle reduced;
		darcy->FormLinearSystem( essentialFluxTdofs, solution, rhs, reduced, traceX, traceB,
		                         true );

#ifdef MFEM_USE_SUITESPARSE
		mfem::UMFPackSolver solver;
		solver.Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
		solver.SetOperator( *reduced.Ptr() );
		solver.Mult( traceB, traceX );
#else
		// Only a fallback. The hybridized trace system is small but not symmetric
		// positive definite in this sign convention, so GMRES rather than CG.
		mfem::SparseMatrix &reducedMatrix = *reduced.As<mfem::SparseMatrix>();
		mfem::GSSmoother preconditioner( reducedMatrix );
		mfem::GMRESSolver solver;
		solver.SetOperator( reducedMatrix );
		solver.SetPreconditioner( preconditioner );
		solver.SetRelTol( 1.0e-12 );
		solver.SetAbsTol( 0.0 );
		solver.SetMaxIter( 5000 );
		solver.SetPrintLevel( -1 );
		solver.Mult( traceB, traceX );
		if ( !solver.GetConverged() )
			throw std::runtime_error( "meq::GradShafranovSolver::solve: the trace solve did not converge" );
#endif

		darcy->RecoverFEMSolution( traceX, rhs, solution );

		// The one place the sign convention is undone. See the file comment.
		fluxGf = darcyFlux;
		fluxGf.Neg();
	}

	mfem::GridFunction &GradShafranovSolver::potential()
	{
		return potentialGf;
	}

	mfem::GridFunction const &GradShafranovSolver::potential() const
	{
		return potentialGf;
	}

	mfem::GridFunction &GradShafranovSolver::flux()
	{
		return fluxGf;
	}

	mfem::GridFunction const &GradShafranovSolver::flux() const
	{
		return fluxGf;
	}

	mfem::GridFunction &GradShafranovSolver::trace()
	{
		return traceGf;
	}

	mfem::GridFunction const &GradShafranovSolver::trace() const
	{
		return traceGf;
	}

	mfem::FiniteElementSpace &GradShafranovSolver::fluxSpace()
	{
		return *fluxFes;
	}

	mfem::FiniteElementSpace &GradShafranovSolver::potentialSpace()
	{
		return *potentialFes;
	}

	mfem::FiniteElementSpace &GradShafranovSolver::traceSpace()
	{
		return *traceFes;
	}

	namespace
	{
		/// A quadrature rule of order 2k+4 on every geometry, so that the error
		/// integral is not itself what limits a measured rate.
		void errorRules( int order, mfem::IntegrationRule const **irs )
		{
			int const quadratureOrder = 2*order + 4;
			for ( int i = 0; i < mfem::Geometry::NumGeom; ++i )
				irs[ i ] = &( mfem::IntRules.Get( i, quadratureOrder ) );
		}
	}

	double GradShafranovSolver::potentialError( mfem::Coefficient &exact ) const
	{
		mfem::IntegrationRule const *irs[ mfem::Geometry::NumGeom ];
		errorRules( orderValue, irs );
		return potentialGf.ComputeL2Error( exact, irs );
	}

	double GradShafranovSolver::fluxError( mfem::VectorCoefficient &exact ) const
	{
		mfem::IntegrationRule const *irs[ mfem::Geometry::NumGeom ];
		errorRules( orderValue, irs );
		return fluxGf.ComputeL2Error( exact, irs );
	}

	int GradShafranovSolver::numTraceDofs() const
	{
		return traceFes->GetVSize();
	}

	int GradShafranovSolver::order() const
	{
		return orderValue;
	}

	double GradShafranovSolver::tau() const
	{
		return stabilization.tau();
	}

}
