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

	SourceIntegrator::SourceIntegrator( Source const &sourceIn, int extraOrderIn )
		: source( &sourceIn ), extraOrder( extraOrderIn )
	{
	}

	mfem::IntegrationRule const &SourceIntegrator::rule( mfem::FiniteElement const &el,
	                                                     mfem::ElementTransformation &tr ) const
	{
		// 2k for the product of the shape function with a field of the same
		// degree, plus the transformation's own weight, plus room for the fact
		// that F is not a polynomial in psi at all: Example 5's source carries
		// exp( -psi ). The residual and the Jacobian use this same rule, so
		// whatever it costs in accuracy it costs consistently and Newton still
		// converges quadratically to the answer the rule defines.
		int const quadratureOrder = 2*el.GetOrder() + tr.OrderW() + extraOrder;
		return mfem::IntRules.Get( el.GetGeomType(), quadratureOrder );
	}

	void SourceIntegrator::AssembleElementVector( mfem::FiniteElement const &el,
	                                              mfem::ElementTransformation &tr,
	                                              mfem::Vector const &elfun,
	                                              mfem::Vector &elvect )
	{
		int const dof = el.GetDof();
		shape.SetSize( dof );
		elvect.SetSize( dof );
		elvect = 0.0;

		mfem::IntegrationRule const &ir = rule( el, tr );
		mfem::Vector point;

		for ( int i = 0; i < ir.GetNPoints(); ++i )
		{
			mfem::IntegrationPoint const &ip = ir.IntPoint( i );
			tr.SetIntPoint( &ip );
			el.CalcShape( ip, shape );
			tr.Transform( ip, point );

			double const r = point( 0 );
			double const z = point( 1 );
			double const psi = shape*elfun;
			double const weight = ip.weight*tr.Weight();

			elvect.Add( -weight*source->f( r, z, psi )/r, shape );
		}
	}

	void SourceIntegrator::AssembleElementGrad( mfem::FiniteElement const &el,
	                                            mfem::ElementTransformation &tr,
	                                            mfem::Vector const &elfun,
	                                            mfem::DenseMatrix &elmat )
	{
		int const dof = el.GetDof();
		shape.SetSize( dof );
		elmat.SetSize( dof );
		elmat = 0.0;

		mfem::IntegrationRule const &ir = rule( el, tr );
		mfem::Vector point;

		for ( int i = 0; i < ir.GetNPoints(); ++i )
		{
			mfem::IntegrationPoint const &ip = ir.IntPoint( i );
			tr.SetIntPoint( &ip );
			el.CalcShape( ip, shape );
			tr.Transform( ip, point );

			double const r = point( 0 );
			double const z = point( 1 );
			double const psi = shape*elfun;
			double const weight = ip.weight*tr.Weight();

			mfem::AddMult_a_VVt( -weight*source->dFdPsi( r, z, psi )/r, shape, elmat );
		}
	}

	namespace
	{
		/// Copies the l2 norm of the non-linear residual out of NewtonSolver at
		/// every iteration. MFEM's own print level would put the same numbers on
		/// stdout, but stage 4's acceptance criterion is an assertion on the
		/// *order* of the convergence, so the history has to be a value rather
		/// than a log line.
		class ResidualRecorder : public mfem::IterativeSolverMonitor
		{
			public:
				explicit ResidualRecorder( std::vector<double> &historyIn )
					: history( historyIn )
				{
				}

				/// MFEM's spelling, from IterativeSolverController. Called on the
				/// first iteration, which is the natural place to start a fresh
				/// history and keeps a second solve from appending to the first.
				void Reset() override // NOLINT(readability-identifier-naming)
				{
					mfem::IterativeSolverMonitor::Reset();
					history.clear();
				}

				/// MFEM's spelling, from IterativeSolverController.
				///
				/// The final callback is skipped deliberately. NewtonSolver reports
				/// every iterate from inside the loop and then reports final_norm
				/// once more after it, which on every exit path is the norm it has
				/// just reported. Keeping the copy would put a ratio of exactly one
				/// on the end of the history and an observed order of zero with it,
				/// which is a convergence failure that did not happen.
				void MonitorResidual( int, mfem::real_t norm, // NOLINT(readability-identifier-naming)
				                      mfem::Vector const &, bool isFinal ) override
				{
					if ( !isFinal )
						history.push_back( norm );
				}

			private:
				std::vector<double> &history;
		};
	}

	GradShafranovSolver::GradShafranovSolver( mfem::Mesh &meshIn, int orderIn, double tauIn )
		: mesh( meshIn ),
		  orderValue( orderIn ),
		  stabilization( tauIn ),
		  radius( []( mfem::Vector const &x ) { return x( 0 ); } ),
		  negativeInverseRadius( []( mfem::Vector const &x ) { return -1.0/x( 0 ); } ),
		  linearSource( nullptr ),
		  nonlinearSource( nullptr ),
		  boundaryData( nullptr ),
		  transferPath( nullptr ),
		  extensionLineOrder( -1 ),
		  newtonRelativeTolerance( 1.0e-12 ),
		  newtonAbsoluteTolerance( 1.0e-14 ),
		  newtonMaxIterations( 30 ),
		  newtonIterationCount( 0 ),
		  built( false ),
		  prepared( false )
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

		dirichletMarker.SetSize( mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0 );
		dirichletMarker = 1;

		gammaHMarker.SetSize( dirichletMarker.Size() );
		gammaHMarker = 0;
		fittedMarker.SetSize( dirichletMarker.Size() );
		fittedMarker = 1;

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
		if ( built )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: the forms are already built; the source has to be set before the first solve" );
		if ( nonlinearSource )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: a psi-dependent source is already set; one solver holds one source" );

		linearSource = &fIn;

		// The potential right hand side is -( F/r, w ), and both signs in that are
		// real. The 1/r is the equation's: the right hand side is F/r, not F. The
		// minus is DarcyForm's: constructed with its default bsymmetrize = true it
		// assembles the second block row as -B q - Mp psi = bp, so the datum handed
		// to it is the negative of the source of div q. Measured: with +F/r the L2
		// error against the exact solution is flat at 7.3e-2 through four
		// refinements, with -F/r it converges at k+1.
		potentialRhsCoeff = std::make_unique<mfem::ProductCoefficient>( negativeInverseRadius,
		                                                                *linearSource );
	}

	void GradShafranovSolver::setSource( Source const &fIn )
	{
		if ( built )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: the forms are already built; the source has to be set before the first solve" );
		if ( linearSource )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: a psi-independent source is already set; one solver holds one source" );

		// Nothing else to do here. Unlike the linear case there is no coefficient
		// to build: the source does not reach a right hand side at all, it becomes
		// a SourceIntegrator on the non-linear potential mass form, where Newton
		// can differentiate it. See buildForms().
		nonlinearSource = &fIn;
	}

	void GradShafranovSolver::setBoundaryData( mfem::Coefficient &boundaryIn )
	{
		boundaryData = &boundaryIn;
	}

	void GradShafranovSolver::setExtension( mfem::TransferPath &pathIn,
	                                        mfem::Array<int> const &gammaHMarkerIn,
	                                        int lineOrderIn )
	{
		if ( built )
			throw std::logic_error( "meq::GradShafranovSolver::setExtension: the forms are already built; the extension has to be set before the first solve" );
		if ( gammaHMarkerIn.Size() != dirichletMarker.Size() )
			throw std::invalid_argument( "meq::GradShafranovSolver::setExtension: the Gamma_h marker must be sized by the largest boundary attribute of the mesh" );

		transferPath = &pathIn;
		extensionLineOrder = lineOrderIn;

		gammaHMarkerIn.Copy( gammaHMarker );
		for ( int i = 0; i < fittedMarker.Size(); ++i )
			fittedMarker[ i ] = gammaHMarker[ i ] ? 0 : 1;

		bool any = false;
		for ( int i = 0; i < gammaHMarker.Size(); ++i )
			any = any || gammaHMarker[ i ];
		if ( !any )
			throw std::invalid_argument( "meq::GradShafranovSolver::setExtension: the Gamma_h marker selects no boundary attribute" );
	}

	bool GradShafranovSolver::isExtended() const
	{
		return transferPath != nullptr;
	}

	void GradShafranovSolver::setNewtonControl( double relativeToleranceIn,
	                                            double absoluteToleranceIn,
	                                            int maxIterationsIn )
	{
		newtonRelativeTolerance = relativeToleranceIn;
		newtonAbsoluteTolerance = absoluteToleranceIn;
		newtonMaxIterations = maxIterationsIn;
	}

	bool GradShafranovSolver::isNonlinear() const
	{
		return nonlinearSource != nullptr;
	}

	void GradShafranovSolver::buildForms()
	{
		if ( built )
			return;
		if ( !linearSource && !nonlinearSource )
			throw std::logic_error( "meq::GradShafranovSolver: no source has been set" );

		darcy = std::make_unique<mfem::DarcyForm>( fluxFes.get(), potentialFes.get() );

		// ( r q, v ). DarcyForm's flux mass form holds the INVERSE of the diffusion
		// coefficient -- convdiff puts 1/k there -- and the coefficient here is
		// 1/r, so this is r. Measured rather than assumed: putting 1/r here instead
		// still converges, to a different function, with the L2 error against the
		// exact Solov'ev solution flat at 1.9e-2 through four refinements.
		mfem::BilinearForm *fluxMass = darcy->GetFluxMassForm();
		fluxMass->AddDomainIntegrator( new mfem::VectorMassIntegrator( radius ) );

		if ( transferPath )
		{
			// < L_e( q_h ), v.n > on Gamma_h: the solution-dependent half of the
			// transferred datum, which is the whole of it for a homogeneous g.
			// Two arguments here were measured rather than argued.
			//
			// The coefficient is radius, the same r the flux mass form carries.
			// HDGExtensionIntegrator documents C as "the same coefficient the
			// flux mass form carries", and CLAUDE.md's mapping table says the
			// same thing from meq's side, but it was checked: with 1/r here the
			// error is flat under refinement, exactly as it is when the flux mass
			// form itself is given 1/r.
			//
			// The sign is +1, HDGExtensionIntegrator's own default, and it is the
			// default for the same reason it is right here: DarcyForm's flux block
			// holds -q, which is precisely the u = -K grad p of the Darcy problem
			// the extension was written for, so meq's convention and the
			// integrator's coincide. Measured: with -1 the rates collapse. See
			// tests/convergence/ExtensionConvergence.cpp for the numbers.
			fluxMass->AddBdrFaceIntegrator(
				new mfem::HDGExtensionIntegrator( *transferPath, radius, +1.0,
				                                  extensionLineOrder ),
				gammaHMarker );
		}

		// < tau( psi_h - psihat_h ), w > on every face of every element, interior
		// and boundary alike. The coefficient handed to HDGDiffusionIntegrator is
		// dead weight once SetStabilization() is called -- the hook replaces the
		// built-in expression entirely -- but it is the diffusion coefficient the
		// integrator is documented to take, so it is the honest thing to pass.
		auto *interior = new mfem::HDGDiffusionIntegrator( radius, stabilization.tau() );
		auto *boundary = new mfem::HDGDiffusionIntegrator( radius, stabilization.tau() );
		interior->SetStabilization( stabilization );
		boundary->SetStabilization( stabilization );

		if ( nonlinearSource )
		{
			// THE WHOLE POTENTIAL BLOCK GOES ON THE NON-LINEAR FORM, not just the
			// source. This is not tidiness, it is what DarcyHybridization requires,
			// and it is the one structural decision in this file that is not
			// obvious from the weak form.
			//
			// The two forms cannot be mixed. Leaving the HDG face stabilisation on
			// the linear potential mass form while the source sits on the
			// non-linear one is refused outright -- MFEM aborts in
			// SetPotMassNonlinearIntegrator() with "Non-linear mass cannot work
			// with a linear constraint", because M_p's face integrators become the
			// hybridization's linear potential constraint. Measured, not reasoned.
			//
			// A domain integrator left on M_p would be worse, because it would not
			// abort: LocalNLOperator::AddMultDE() reads
			//
			//     if ( m_nlfi_p ) { ... }  else if ( !D_empty ) { ... linear D ... }
			//
			// and ConstructGrad() has the same shape, so with a non-linear
			// potential mass present the assembled linear D is simply not visited.
			// miniapps/hdg/convdiff.cpp puts every potential term -- domain terms,
			// face stabilisation, convection -- on whichever of the two forms is in
			// use, and so does this.
			mfem::NonlinearForm *potentialMass = darcy->GetPotentialMassNonlinearForm();
			potentialMass->AddDomainIntegrator( new SourceIntegrator( *nonlinearSource ) );
			potentialMass->AddInteriorFaceIntegrator( interior );
			potentialMass->AddBdrFaceIntegrator( boundary, fittedMarker );
		}
		else
		{
			mfem::BilinearForm *potentialMass = darcy->GetPotentialMassForm();
			potentialMass->AddInteriorFaceIntegrator( interior );
			potentialMass->AddBdrFaceIntegrator( boundary, fittedMarker );
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
			fittedMarker );

		// < qhat_h.n, mu > = 0. This must come after every AddIntegrator above:
		// EnableHybridization() reaches into the potential mass and flux divergence
		// forms as it runs and takes what is there at that moment.
		darcy->EnableHybridization( traceFes.get(), new mfem::NormalTraceJumpIntegrator(),
		                            essentialFluxTdofs );
		darcy->GetHybridization()->SetEssentialBC( dirichletMarker );

		if ( nonlinearSource )
		{
			// Hybridization eliminates the flux and the potential element by
			// element, and when the potential block is non-linear that elimination
			// is itself a small non-linear solve, one per element per residual
			// evaluation. Two things about it are load bearing.
			//
			// Its tolerance is part of the global residual. Solve the local
			// problems loosely and the outer Newton is differentiating a function
			// it is not quite evaluating, which does not raise an error -- it
			// shows up as an outer residual history that stalls short of round-off
			// and an observed order below two. 1e-12 relative is tight enough that
			// the outer iteration reaches 4e-14 in four steps; the elements are
			// small and the extra inner steps cost nothing measurable.
			//
			// And the inner solver is Newton as well, with a dense LU rather than
			// the default GMRES on the local Jacobian, so that neither the inner
			// rate nor an inner linear tolerance can be what a stalled outer
			// history is blamed on.
			darcy->GetHybridization()->SetLocalNLSolver(
				mfem::DarcyHybridization::LSsolveType::Newton, 100, 1.0e-12, 1.0e-16, -1 );
			darcy->GetHybridization()->SetLocalNLPreconditioner(
				mfem::DarcyHybridization::LPrecType::LU );
		}

		darcy->Assemble();
		built = true;
	}

	void GradShafranovSolver::prepare()
	{
		if ( !linearSource && !nonlinearSource )
			throw std::logic_error( "meq::GradShafranovSolver::prepare: no source has been set" );
		// On the extension path a Gamma_h attribute carries no datum of its own --
		// what is imposed there is phi_h -- so boundary data is needed only if some
		// attribute is still fitted.
		bool anyFitted = false;
		for ( int i = 0; i < fittedMarker.Size(); ++i )
			anyFitted = anyFitted || fittedMarker[ i ];
		if ( !boundaryData && anyFitted )
			throw std::logic_error( "meq::GradShafranovSolver::prepare: no boundary data has been set" );

		buildForms();

		solution = 0.0;
		rhs = 0.0;

		// The Dirichlet datum lives on the trace, and only on the trace: the flux
		// right hand side stays zero because the essential trace condition, not a
		// boundary linear form, is what imposes psi = g_D here. On the non-linear
		// path this is doubly load bearing -- the reduced operator masks its
		// residual to zero on the essential trace dofs and puts a unit row in the
		// Jacobian there, so the value Newton starts from on those dofs is the
		// value it finishes with. That path through DarcyHybridization was broken
		// until recently (EliminateTraceTrueDofsInRHS returned early for non-linear
		// problems and the condition was ignored outright) and no MFEM regression
		// covers the combination; if a converged answer ever looks wrong near
		// Gamma, look there before looking here.
		// fittedMarker, not dirichletMarker: the trace dofs of Gamma_h are pinned
		// to zero rather than to a datum, since nothing references them. See
		// setExtension().
		if ( boundaryData && anyFitted )
			traceGf.ProjectBdrCoefficient( *boundaryData, fittedMarker );

		if ( linearSource )
		{
			// Order 2k+4, not DomainLFIntegrator's default 2k. F/r is not a
			// polynomial -- for Solov'ev it is a rational function, and for a
			// tabulated profile it is a spline in psi -- so the default rule
			// integrates it to its own accuracy rather than to the solution's.
			// The number matters for a second reason: SourceIntegrator uses
			// 2k+4 as well, so on the affine meshes stage 4 uses, the same F
			// handed to either overload of setSource() is integrated by the
			// same rule and the two paths agree to round-off. They do not
			// otherwise -- measured, the Solov'ev problem solved both ways
			// differs by 4.3e-7 relative with the default rule here and by
			// 3.2e-14 with this one, and
			// NewtonConvergence.cpp's theNewtonPathReproducesTheLinearPathOnSoloviev
			// is what watches it. On the Solov'ev convergence table this moved
			// no rate and no more than the fourth significant figure of any
			// error.
			mfem::LinearForm potentialRhs;
			potentialRhs.Update( potentialFes.get(), rhs.GetBlock( 1 ), 0 );
			potentialRhs.AddDomainIntegrator(
				new mfem::DomainLFIntegrator( *potentialRhsCoeff, 2, 4 ) );
			potentialRhs.Assemble();
		}

		// traceX and traceB alias the trace blocks of the solution and right hand
		// side. That aliasing is not cosmetic: FormLinearSystem() only calls
		// EliminateTraceTrueDofsInRHS() -- the step that moves the essential trace
		// values into the reduced problem -- when X arrives already sized to the
		// reduced system and copy_interior is set.
		traceX.MakeRef( solution, blockOffsets[ 2 ], traceFes->GetVSize() );
		traceB.MakeRef( rhs, blockOffsets[ 2 ], traceFes->GetVSize() );

		// DarcyForm knows about two blocks, not three; see the header.
		darcySolution.Update( solution, darcy->GetOffsets() );
		darcyRhs.Update( rhs, darcy->GetOffsets() );

		darcy->FormLinearSystem( essentialFluxTdofs, darcySolution, darcyRhs, reduced,
		                         traceX, traceB, true );

		prepared = true;
	}

	void GradShafranovSolver::solve()
	{
		prepare();

		if ( nonlinearSource )
		{
			ResidualRecorder recorder( newtonResidualHistory );

#ifdef MFEM_USE_SUITESPARSE
			mfem::UMFPackSolver linear;
			linear.Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
#else
			mfem::GMRESSolver linear;
			linear.SetRelTol( 1.0e-14 );
			linear.SetAbsTol( 0.0 );
			linear.SetMaxIter( 5000 );
			linear.SetPrintLevel( -1 );
#endif

			// The reduced operator is DarcyHybridization itself, whose GetGradient()
			// differentiates the assembled residual rather than the continuous
			// equation. That is the point of doing it this way: CEDRES++ rejected
			// the continuous-level derivative because it "seems to blow up if psi
			// reaches a critical point", and a Jacobian obtained by differentiating
			// the discrete residual cannot disagree with the residual by
			// construction. What it can disagree with is the physics, and only a
			// convergence study against a closed form catches that.
			mfem::NewtonSolver newton;
			newton.SetOperator( *reduced.Ptr() );
			newton.SetSolver( linear );
			newton.SetMonitor( recorder );
			newton.SetRelTol( newtonRelativeTolerance );
			newton.SetAbsTol( newtonAbsoluteTolerance );
			newton.SetMaxIter( newtonMaxIterations );
			newton.SetPrintLevel( -1 );

			// The Dirichlet data rides in traceX, so the iteration must start from
			// it rather than from zero. With iterative_mode false NewtonSolver
			// zeroes x on entry and the boundary condition disappears without a
			// word -- the residual is masked on those rows, so nothing complains.
			newton.iterative_mode = true;

			newton.Mult( traceB, traceX );
			newtonIterationCount = newton.GetNumIterations();

			// Loudly, and without a recovered solution. A Newton iteration that
			// ran out of iterations has produced a vector, not an equilibrium, and
			// on this problem a failure to converge means something is wrong with
			// the Jacobian rather than with the initial guess -- there is no
			// globalisation here to have been the thing that failed.
			if ( !newton.GetConverged() )
				throw std::runtime_error( "meq::GradShafranovSolver::solve: the Newton iteration did not converge" );
		}
		else
		{
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
		}

		darcy->RecoverFEMSolution( traceX, darcyRhs, darcySolution );

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

	mfem::Operator &GradShafranovSolver::reducedOperator()
	{
		if ( !prepared )
			throw std::logic_error( "meq::GradShafranovSolver::reducedOperator: prepare() has not been called" );
		return *reduced.Ptr();
	}

	mfem::Vector &GradShafranovSolver::reducedRhs()
	{
		if ( !prepared )
			throw std::logic_error( "meq::GradShafranovSolver::reducedRhs: prepare() has not been called" );
		return traceB;
	}

	mfem::Vector &GradShafranovSolver::reducedSolution()
	{
		if ( !prepared )
			throw std::logic_error( "meq::GradShafranovSolver::reducedSolution: prepare() has not been called" );
		return traceX;
	}

	mfem::Array<int> const &GradShafranovSolver::essentialTraceDofs() const
	{
		if ( !built )
			throw std::logic_error( "meq::GradShafranovSolver::essentialTraceDofs: the forms have not been built" );
		return darcy->GetHybridization()->GetEssentialTrueDofs();
	}

	std::vector<double> const &GradShafranovSolver::newtonResiduals() const
	{
		return newtonResidualHistory;
	}

	int GradShafranovSolver::newtonIterations() const
	{
		return newtonIterationCount;
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
