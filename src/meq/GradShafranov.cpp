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
		  initialGuess( nullptr ),
		  globalisationChoice( Globalisation::None ),
		  localSolverChoice( LocalSolver::Newton ),
		  orderingChoice( NonlinearOrdering::CondenseThenLinearise ),
		  andersonDepth( 1 ),
		  picardDamping( 1.0 ),
		  transferPath( nullptr ),
		  extensionLineOrder( -1 ),
		  newtonRelativeTolerance( 1.0e-12 ),
		  newtonAbsoluteTolerance( 1.0e-14 ),
		  newtonMaxIterations( 30 ),
		  newtonIterationCount( 0 ),
		  built( false ),
		  prepared( false ),
		  postProcessed( false )
	{
		if ( orderValue < 0 )
			throw std::invalid_argument( "meq::GradShafranovSolver: the polynomial order must not be negative" );
		if ( mesh.Dimension() != 2 )
			throw std::invalid_argument( "meq::GradShafranovSolver: the mesh must be two dimensional ( r, z )" );

		int const dim = mesh.Dimension();

		/*
		 * The closed (Gauss-Lobatto) basis, as miniapps/hdg/convdiff.cpp puts it,
		 * "as it is customary for HDG to match trace DOFs". All three spaces carry
		 * the same degree; hybridization is what makes that legal.
		 *
		 * IT IS A CONVENTION AND NOT A REQUIREMENT, and that was an inherited
		 * quotation until it was measured. A nodal basis does not change the
		 * SPACE -- Gauss-Lobatto and Gauss-Legendre both span P_k( K ), and only
		 * the shape functions differ -- so the discretisation cannot see the
		 * choice. Measured, with both volume spaces switched to GaussLegendre and
		 * nothing else touched: SolovievConvergence gives 1.996/3.001/3.998 in psi
		 * at k = 1, 2, 3, ExtensionConvergence's curved benchmark gives L2
		 * 2.742813e-05 down both paths agreeing to 6.6e-15, and the k = 3 Newton
		 * history reads 1.121994e+01, 4.085930e-02, 5.744224e-04, 1.411244e-07 --
		 * the same digits this file records for Lobatto, differing only in the
		 * last place. Nothing in the hybridization needs the volume dofs to sit on
		 * the faces: every coupling is a face INTEGRAL, computed by quadrature
		 * against both bases, and integrals do not care where dofs live.
		 *
		 * So it is kept for alignment with the miniapp meq was ported from, which
		 * is worth something when debugging against MFEM, and for nothing else.
		 * WHAT IT COSTS is that a dof is a point value ON the element boundary,
		 * where an L2 field is discontinuous -- so reading this space by nodal
		 * interpolation at another mesh's dof points is ambiguous, and measured,
		 * 9% to 28% wrong. meq::FieldTransfer therefore projects rather than
		 * interpolates, which is basis-agnostic and is the right thing for
		 * non-nested meshes anyway. Switching the basis would remove that
		 * ambiguity and buy nothing else, which is why it was not switched.
		 */
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

	/*
	 * F( r, z, psi^k( r, z ) ): the source frozen at the previous iterate.
	 *
	 * This is what makes a Picard path linear. setSource( Source const & ) puts
	 * the source on the NON-LINEAR potential mass form, where hybridization
	 * turns each element's elimination into its own Newton. Handing the same
	 * Source through this coefficient instead puts it on the right hand side,
	 * the potential block stays linear, and every local elimination is a linear
	 * solve -- which is the ordering Nguyen, Peraire & Cockburn use and meq's
	 * Newton path does not. See CLAUDE.md.
	 */
	class FrozenSource : public mfem::Coefficient
	{
		public:
			FrozenSource( Source const &sourceIn, mfem::GridFunction const &psiIn )
				: source( sourceIn ), psi( psiIn )
			{
			}

			/// MFEM's spelling, from mfem::Coefficient.
			double Eval( mfem::ElementTransformation &tr, // NOLINT(readability-identifier-naming)
			             mfem::IntegrationPoint const &ip ) override
			{
				mfem::Vector x;
				tr.Transform( ip, x );
				return source.f( x( 0 ), x( 1 ), psi.GetValue( tr, ip ) );
			}

		private:
			Source const &source;
			mfem::GridFunction const &psi;
	};

	/*
	 * The Picard map as an Operator, so KINSOL can iterate it.
	 *
	 * KINSOL reads Mult() differently by strategy: for KIN_NONE and
	 * KIN_LINESEARCH it is the residual F( u ) and the target is F = 0; for
	 * KIN_FP it is the fixed point map G( u ) and the target is u = G( u ). This
	 * is the latter -- one frozen-source assembly and one linear solve per call.
	 */
	class PicardMap : public mfem::Operator
	{
		public:
			using Step = std::function<void( mfem::Vector const &, mfem::Vector & )>;

			PicardMap( int size, Step stepIn )
				: mfem::Operator( size ), step( std::move( stepIn ) )
			{
			}

			/// MFEM's spelling, from mfem::Operator.
			void Mult( mfem::Vector const &x, // NOLINT(readability-identifier-naming)
			           mfem::Vector &y ) const override
			{
				step( x, y );
			}

		private:
			Step step;
	};

	/*
	 * R( x ) - b as an Operator, because KINSOL will not take a right hand side.
	 *
	 * mfem::KINSolver derives from mfem::NewtonSolver, which reads as though the
	 * two were interchangeable, and CLAUDE.md said so. They are not, and the
	 * difference is silent: NewtonSolver::Mult( b, x ) forms r = oper( x ) - b,
	 * while KINSolver::Mult declares its first argument WITHOUT A NAME and solves
	 * oper( x ) = 0. Hand KINSOL a problem with a non-zero right hand side and it
	 * converges -- to the solution of a different problem.
	 *
	 * This reproduces NewtonSolver's residual exactly, which is what keeps the
	 * comparison between the two honest: whatever b holds on the essential trace
	 * rows, both paths subtract the same thing.
	 *
	 * GetGradient forwards untouched. The shift is constant, so it contributes
	 * nothing to the Jacobian.
	 */
	class ShiftedResidual : public mfem::Operator
	{
		public:
			ShiftedResidual( mfem::Operator &operatorIn, mfem::Vector const &rhsIn )
				: mfem::Operator( operatorIn.Height(), operatorIn.Width() ),
				  residual( operatorIn ), shift( rhsIn )
			{
			}

			/// MFEM's spelling, from mfem::Operator.
			void Mult( mfem::Vector const &x, // NOLINT(readability-identifier-naming)
			           mfem::Vector &y ) const override
			{
				residual.Mult( x, y );
				y -= shift;
			}

			/// MFEM's spelling, from mfem::Operator.
			mfem::Operator &GetGradient( // NOLINT(readability-identifier-naming)
				mfem::Vector const &x ) const override
			{
				return residual.GetGradient( x );
			}

		private:
			mfem::Operator &residual;
			mfem::Vector const &shift;
	};

	void GradShafranovSolver::setBoundaryData( mfem::Coefficient &boundaryIn )
	{
		boundaryData = &boundaryIn;
	}

	void GradShafranovSolver::setInitialGuess( mfem::Coefficient &psiGuess )
	{
		ownedInitialGuess.reset();
		initialGuess = &psiGuess;
		prepared = false;
	}

	void GradShafranovSolver::setInitialGuess( mfem::GridFunction const &psiGuess )
	{
		// GridFunctionCoefficient takes a non-const pointer but only reads, which
		// is why the public interface can promise const.
		ownedInitialGuess = std::make_unique<mfem::GridFunctionCoefficient>(
			const_cast<mfem::GridFunction *>( &psiGuess ) );
		initialGuess = ownedInitialGuess.get();
		prepared = false;
	}

	void GradShafranovSolver::setGlobalisation( Globalisation choice )
	{
#ifndef MFEM_USE_SUNDIALS
		if ( choice != Globalisation::None )
			throw std::logic_error(
				"meq::GradShafranovSolver::setGlobalisation: MFEM was built without "
				"MFEM_USE_SUNDIALS, so no KINSOL strategy is available" );
#endif
		globalisationChoice = choice;
		// built, not just prepared. usesNonlinearForms() reads globalisationChoice,
		// and buildForms() branches on it to decide whether the potential block
		// goes on the LINEAR or the NON-LINEAR form -- so switching a live solver
		// between a Picard path and a Newton one while built stayed true would
		// silently reuse the other path's blocks. Found while wiring
		// PicardThenNewton, which switches twice in one solve(); setNonlinearOrdering
		// and setLocalSolver beside it always did reset built, and this did not.
		built = false;
		prepared = false;
	}

	void GradShafranovSolver::setPicardDamping( double damping )
	{
		if ( !( damping > 0.0 ) || damping > 1.0 )
			throw std::logic_error( "meq::GradShafranovSolver::setPicardDamping: the damping must be in ( 0, 1 ]" );
		picardDamping = damping;
	}

	void GradShafranovSolver::setAndersonDepth( int depth )
	{
		if ( depth < 0 )
			throw std::logic_error( "meq::GradShafranovSolver::setAndersonDepth: the depth cannot be negative" );
		andersonDepth = depth;
	}

	/*
	 * One Picard step: freeze F at @a in, reassemble, solve, hand back the new
	 * potential. This is the fixed point map, and KINSOL calls it once per
	 * iteration.
	 *
	 * prepare() is re-entered deliberately. buildForms() is guarded by `built`,
	 * so the bilinear forms are assembled once and only the right hand side and
	 * the hybridization's reduction are redone -- which they must be, the source
	 * having changed.
	 */
	void GradShafranovSolver::picardStep( mfem::Vector const &in, mfem::Vector &out )
	{
		*picardIterate = in;

		prepared = false;
		prepare();

#ifdef MFEM_USE_SUITESPARSE
		// Held across calls rather than built per iteration, which is the whole
		// point: Picard runs 122 to 290 of these, each a full factorisation of a
		// matrix whose sparsity never changes. prepare() rebuilds `reduced` every
		// iteration, but SetReuseSymbolic() compares the pattern rather than the
		// object -- it documents accepting "a matrix rebuilt into a fresh object
		// with the same structure" -- so the analysis survives that.
		if ( !picardSolver )
		{
			picardSolver = std::make_unique<mfem::UMFPackSolver>();
			picardSolver->Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
			picardSolver->SetReuseSymbolic();
		}
		picardSolver->SetOperator( *reduced.Ptr() );
		picardSolver->Mult( traceB, traceX );
#else
		mfem::SparseMatrix &matrix = *reduced.As<mfem::SparseMatrix>();
		mfem::GSSmoother preconditioner( matrix );
		mfem::GMRESSolver step;
		step.SetOperator( matrix );
		step.SetPreconditioner( preconditioner );
		step.SetRelTol( 1.0e-12 );
		step.SetAbsTol( 0.0 );
		step.SetMaxIter( 5000 );
		step.SetPrintLevel( -1 );
		step.Mult( traceB, traceX );
		if ( !step.GetConverged() )
			throw std::runtime_error( "meq::GradShafranovSolver::picardStep: the trace solve did not converge" );
#endif

		darcy->RecoverFEMSolution( traceX, darcyRhs, darcySolution );
		out = potentialGf;
	}

	/*
	 * Anderson-accelerated Picard, which is the GS papers' own method.
	 *
	 * KINSOL's KIN_FP reads the operator as the fixed point map G( u ) and drives
	 * u = G( u ); KINSetMAA turns on Anderson with the given subspace depth.
	 * MFEM's KINSolver::Mult forwards to oper->Mult() and lets KINSOL interpret
	 * it, so the same wrapper serves both readings -- which is convenient and is
	 * also the trap recorded in setGlobalisation().
	 *
	 * The unknown is psi_h in W_h, not the trace. That is a different iteration
	 * from every other path in this file, and it is why this is a separate
	 * function rather than another case in solve()'s switch.
	 */
	void GradShafranovSolver::solveByPicard()
	{
#ifndef MFEM_USE_SUNDIALS
		throw std::logic_error(
			"meq::GradShafranovSolver::solve: a Picard globalisation was asked for "
			"but MFEM was built without MFEM_USE_SUNDIALS" );
#else
		// One assembly to size everything and to give the fixed point its start:
		// the Dirichlet data extended inward, which is what Newton starts from too.
		prepared = false;
		prepare();

		mfem::Vector iterate( potentialGf.Size() );
		iterate = potentialGf;
		if ( initialGuess )
		{
			mfem::GridFunction seeded( potentialFes.get() );
			seeded.ProjectCoefficient( *initialGuess );
			iterate = seeded;
		}

		PicardMap map( iterate.Size(),
			[ this ]( mfem::Vector const &in, mfem::Vector &out )
			{
				picardStep( in, out );
			} );

		mfem::KINSolver fixedPoint( KIN_FP, false );
		// One damping, not two. KINSetDampingAA damps the Anderson combination and
		// KINSetDamping the underlying fixed point; setting both compounds them
		// and is worse than either -- measured, Anderson with both fails at 500
		// iterations where plain damped Picard converges in 194.
		if ( globalisationChoice == Globalisation::AndersonPicard )
			fixedPoint.EnableAndersonAcc( andersonDepth, KIN_ORTH_MGS, 0, picardDamping );
		else
			fixedPoint.SetDamping( picardDamping );
		fixedPoint.SetOperator( map );
		fixedPoint.SetRelTol( newtonRelativeTolerance );
		fixedPoint.SetAbsTol( newtonAbsoluteTolerance );
		fixedPoint.SetMaxIter( newtonMaxIterations );
		fixedPoint.SetPrintLevel( -1 );
		fixedPoint.iterative_mode = true;

		mfem::Vector unused( iterate.Size() );
		unused = 0.0;
		fixedPoint.Mult( unused, iterate );
		newtonIterationCount = fixedPoint.GetNumIterations();

		if ( !fixedPoint.GetConverged() )
			throw std::runtime_error( "meq::GradShafranovSolver::solve: the Picard iteration did not converge" );

		// One last step at the converged iterate, so that the flux, the trace and
		// the potential all come from the SAME assembly. Without it the flux would
		// be one iteration stale, which no rate would notice and every field plot
		// would.
		mfem::Vector settled( iterate.Size() );
		picardStep( iterate, settled );

		fluxGf = darcyFlux;
		fluxGf.Neg();
		postProcessed = false;
#endif
	}

	void GradShafranovSolver::setNonlinearOrdering( NonlinearOrdering choice )
	{
		orderingChoice = choice;
		built = false;
		prepared = false;
	}

	GradShafranovSolver::NonlinearOrdering
	GradShafranovSolver::nonlinearOrdering() const
	{
		return orderingChoice;
	}

	void GradShafranovSolver::setLocalSolver( LocalSolver choice )
	{
		localSolverChoice = choice;
		built = false;
		prepared = false;
	}

	GradShafranovSolver::Globalisation GradShafranovSolver::globalisation() const
	{
		return globalisationChoice;
	}

	void GradShafranovSolver::clearInitialGuess()
	{
		ownedInitialGuess.reset();
		initialGuess = nullptr;
		prepared = false;
	}

	bool GradShafranovSolver::hasInitialGuess() const
	{
		return initialGuess != nullptr;
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

	bool GradShafranovSolver::usesNonlinearForms() const
	{
		if ( !nonlinearSource )
			return false;
		return globalisationChoice != Globalisation::AndersonPicard
		    && globalisationChoice != Globalisation::PicardOnly;
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

		if ( usesNonlinearForms() )
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

		if ( usesNonlinearForms() )
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
			mfem::DarcyHybridization::LSsolveType localType =
				mfem::DarcyHybridization::LSsolveType::Newton;
			switch ( localSolverChoice )
			{
				case LocalSolver::Newton:
					localType = mfem::DarcyHybridization::LSsolveType::Newton;
					break;
				case LocalSolver::Lbfgs:
					localType = mfem::DarcyHybridization::LSsolveType::LBFGS;
					break;
				case LocalSolver::Lbb:
					localType = mfem::DarcyHybridization::LSsolveType::LBB;
					break;
			}
			darcy->GetHybridization()->SetLocalNLSolver(
				localType, 100, 1.0e-12, 1.0e-16, -1 );
			darcy->GetHybridization()->SetLocalNLPreconditioner(
				mfem::DarcyHybridization::LPrecType::LU );

			// And the ordering, which decides whether either of the two above
			// means anything: under LineariseThenCondense the local problem is a
			// linear solve and there is no local iteration to configure.
			darcy->GetHybridization()->SetNonlinearOrdering(
				orderingChoice == NonlinearOrdering::LineariseThenCondense
					? mfem::DarcyHybridization::NLOrdering::LineariseThenCondense
					: mfem::DarcyHybridization::NLOrdering::CondenseThenLinearise );
		}

		darcy->Assemble();
		built = true;
	}

	/*
	 * Interpolate a coefficient onto the hybrid trace space.
	 *
	 * There is no library call for this. GridFunction::ProjectCoefficient loops
	 * over fes->GetNE() -- VOLUME elements -- so on a trace space it never
	 * reaches a face dof; and ProjectBdrCoefficient reaches boundary faces only,
	 * which is where almost none of the dofs are. Checked in
	 * mfem/linalg/../fem/gridfunc.cpp; this is the one place in this file where
	 * the obvious call is the wrong one.
	 *
	 * The pattern is the one Estimator.cpp uses to read a trace value back:
	 * GetFaceElement gives the face's own element, GetFaceVDofs its dofs, and
	 * the two agree on orientation because every face carries its own dofs. The
	 * coefficient is evaluated through the VOLUME transformation of the owning
	 * element rather than the face's own, so that a GridFunctionCoefficient --
	 * which needs an element to look a value up in -- works as well as a
	 * FunctionCoefficient does.
	 *
	 * This is nodal interpolation, not an L2 projection. For a starting point
	 * the difference is immaterial, and DG_Interface_FECollection is nodal with
	 * VALUE map type, so the nodes are where the dofs live.
	 */
	void GradShafranovSolver::projectOntoTrace( mfem::Coefficient &coeff,
	                                            mfem::GridFunction &target ) const
	{
		mfem::Mesh &mesh = *traceFes->GetMesh();
		mfem::Array<int> vdofs;
		mfem::Vector values;

		for ( int f = 0; f < mesh.GetNumFaces(); ++f )
		{
			mfem::FaceElementTransformations *ftr =
				mesh.GetFaceElementTransformations( f );
			if ( !ftr )
				continue;

			mfem::FiniteElement const *faceFe = traceFes->GetFaceElement( f );
			if ( !faceFe )
				continue;

			traceFes->GetFaceVDofs( f, vdofs );
			int const dof = faceFe->GetDof();
			values.SetSize( dof );

			mfem::IntegrationRule const &nodes = faceFe->GetNodes();
			for ( int i = 0; i < dof; ++i )
			{
				ftr->SetAllIntPoints( &nodes.IntPoint( i ) );
				values( i ) = coeff.Eval( *ftr->Elem1, ftr->GetElement1IntPoint() );
			}

			target.SetSubVector( vdofs, values );
		}
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

		// The starting point, before the Dirichlet datum and after the zeroing,
		// so that g_D wins on the essential dofs and the guess supplies the rest.
		// See setInitialGuess() for why this is needed at all: with a source that
		// vanishes at psi = 0, starting from the Dirichlet data alone lands Newton
		// on psi == 0 and it stops there, converged.
		//
		// The trace is the half that matters, being Newton's actual unknown. The
		// potential block is seeded as well, for MFEM's element-local non-linear
		// solves, which start from whatever the block vector holds -- see the
		// header. The flux block is left at zero: nothing iterates from it, and a
		// guess for psi says nothing about q without differentiating it.
		if ( initialGuess && nonlinearSource )
		{
			projectOntoTrace( *initialGuess, traceGf );
			potentialGf.ProjectCoefficient( *initialGuess );
		}

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

		// On a Picard path the source is the frozen coefficient, so the linear
		// right hand side below is what assembles it. Built here rather than in
		// buildForms() because it reads picardIterate, which changes every step.
		mfem::Coefficient *rhsSource = linearSource;
		if ( nonlinearSource && !usesNonlinearForms() )
		{
			if ( !picardIterate )
			{
				picardIterate = std::make_unique<mfem::GridFunction>( potentialFes.get() );
				*picardIterate = 0.0;
			}
			if ( !frozenSource )
				frozenSource = std::make_unique<FrozenSource>( *nonlinearSource,
				                                              *picardIterate );
			rhsSource = frozenSource.get();

			// The same -1/r that setSource( Coefficient & ) applies, and for the
			// same two reasons -- the equation's 1/r and DarcyForm's sign. Built
			// once; ProductCoefficient holds a reference, so the frozen source
			// re-reads picardIterate on every assembly without rebuilding this.
			if ( !potentialRhsCoeff )
				potentialRhsCoeff = std::make_unique<mfem::ProductCoefficient>(
					negativeInverseRadius, *rhsSource );
		}

		if ( rhsSource )
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

	void GradShafranovSolver::solveByPicardThenNewton()
	{
#ifndef MFEM_USE_SUNDIALS
		throw std::logic_error(
			"meq::GradShafranovSolver::solve: Globalisation::PicardThenNewton needs "
			"the Picard stage, and MFEM was built without MFEM_USE_SUNDIALS" );
#else
		// Restored however this exits, so a throw out of either stage does not
		// leave the solver reporting a globalisation it is no longer set up for.
		// built goes with it, for the reason setGlobalisation() records.
		struct Restore
		{
			GradShafranovSolver *solver;
			Globalisation choice;
			~Restore()
			{
				solver->globalisationChoice = choice;
				solver->built = false;
			}
		} restore { this, globalisationChoice };

		// ---- stage 1: Anderson-accelerated Picard, to reach Newton's basin ----
		//
		// Not to solve the problem. GS-2 section 4.5 converges at both orders from
		// a Picard state that never met its own tolerance, so stage 1 stopping
		// short is an expected outcome and not an error. What it leaves behind is
		// the last iterate it evaluated, which is what stage 2 wants.
		setGlobalisation( Globalisation::AndersonPicard );
		try
		{
			solve();
		}
		catch ( std::runtime_error const & )
		{
			// Ran out of iterations. Deliberately swallowed; see above.
		}
		// Held in a local and published only once stage 2 is done, because stage
		// 2 re-enters solve(), which zeroes the member on the way in.
		int const stageOneIterations = newtonIterationCount;

		// A COPY, because setInitialGuess() only references what it is given and
		// stage 2 writes through potentialGf.
		picardSeed = std::make_unique<mfem::GridFunction>( potentialGf );

		// ---- stage 2: plain Newton from there ----
		setGlobalisation( Globalisation::None );
		setInitialGuess( *picardSeed );
		solve();

		picardIterationCount = stageOneIterations;
#endif
	}

	void GradShafranovSolver::solve()
	{
		// The Picard paths iterate a fixed point on the POTENTIAL, not a residual
		// on the trace, so they do not go through prepare()-then-Newton at all --
		// picardStep() re-enters prepare() itself, once per iteration.
		// Stale on a solver reused across globalisations otherwise; the handoff
		// republishes it after stage 2.
		picardIterationCount = 0;

		if ( nonlinearSource && globalisationChoice == Globalisation::PicardThenNewton )
		{
			solveByPicardThenNewton();
			return;
		}

		if ( nonlinearSource && !usesNonlinearForms() )
		{
			solveByPicard();
			return;
		}

		prepare();

		if ( nonlinearSource )
		{
			ResidualRecorder recorder( newtonResidualHistory );

#ifdef MFEM_USE_SUITESPARSE
			mfem::UMFPackSolver linear;
			linear.Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
			// NewtonSolver::Mult calls prec->SetOperator( *grad ) on THIS object
			// once per iteration, and the trace system's sparsity does not change
			// between them -- so without this the METIS analysis is recomputed and
			// thrown away every step, at a fifth to a quarter of the factorisation.
			// The pattern is compared entry by entry rather than assumed, so a
			// pattern that did change is re-analysed and the answer is unaffected.
			linear.SetReuseSymbolic();
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
			// KINSOL ignores the right hand side handed to Mult(), so the shift
			// has to be in the operator. Built unconditionally and used only on
			// the KINSOL paths: it costs one subtraction per residual evaluation
			// and keeps the two paths reading the same residual, which is what
			// makes a difference between them attributable to the line search.
			ShiftedResidual shifted( *reduced.Ptr(), traceB );

			std::unique_ptr<mfem::NewtonSolver> nonlinear;
			bool kinsol = false;

			switch ( globalisationChoice )
			{
				case Globalisation::None:
					nonlinear = std::make_unique<mfem::NewtonSolver>();
					break;
#ifdef MFEM_USE_SUNDIALS
				case Globalisation::LineSearch:
					nonlinear = std::make_unique<mfem::KINSolver>( KIN_LINESEARCH, true );
					kinsol = true;
					break;
				case Globalisation::KinsolNoLineSearch:
					nonlinear = std::make_unique<mfem::KINSolver>( KIN_NONE, true );
					kinsol = true;
					break;
#else
				case Globalisation::LineSearch:
				case Globalisation::KinsolNoLineSearch:
					throw std::logic_error(
						"meq::GradShafranovSolver::solve: a KINSOL globalisation was "
						"asked for but MFEM was built without MFEM_USE_SUNDIALS" );
#endif
			}

			// The reduced operator is DarcyHybridization itself, whose GetGradient()
			// differentiates the assembled residual rather than the continuous
			// equation. That is the point of doing it this way: CEDRES++ rejected
			// the continuous-level derivative because it "seems to blow up if psi
			// reaches a critical point", and a Jacobian obtained by differentiating
			// the discrete residual cannot disagree with the residual by
			// construction. What it can disagree with is the physics, and only a
			// convergence study against a closed form catches that.
			//
			// SetOperator before SetSolver on both paths, because KINSolver
			// requires that order -- its SetSolver documents "must be called after
			// SetOperator()" -- and NewtonSolver does not care.
			if ( kinsol )
				nonlinear->SetOperator( shifted );
			else
				nonlinear->SetOperator( *reduced.Ptr() );
			nonlinear->SetSolver( linear );
			nonlinear->SetMonitor( recorder );
			nonlinear->SetRelTol( newtonRelativeTolerance );
			nonlinear->SetAbsTol( newtonAbsoluteTolerance );
			nonlinear->SetMaxIter( newtonMaxIterations );
			nonlinear->SetPrintLevel( -1 );

			// The Dirichlet data rides in traceX, so the iteration must start from
			// it rather than from zero. With iterative_mode false NewtonSolver
			// zeroes x on entry and the boundary condition disappears without a
			// word -- the residual is masked on those rows, so nothing complains.
			// KINSolver::Mult has the same line and the same consequence.
			nonlinear->iterative_mode = true;

			// traceB on the Newton path, where Mult() subtracts it; ignored on the
			// KINSOL paths, where ShiftedResidual has already done so.
			nonlinear->Mult( traceB, traceX );
			newtonIterationCount = nonlinear->GetNumIterations();

#ifdef MFEM_USE_SUITESPARSE
			// Recorded so the reuse can be asserted on rather than timed: a Newton
			// solve refactorises once per iteration and must analyse ONCE.
			symbolicFactorisationCount = linear.GetNumSymbolicFactorizations();
			numericFactorisationCount = linear.GetNumNumericFactorizations();
#endif

			// Loudly, and without a recovered solution: an iteration that ran out
			// of steps has produced a vector, not an equilibrium.
			if ( !nonlinear->GetConverged() )
				throw std::runtime_error( "meq::GradShafranovSolver::solve: the non-linear iteration did not converge" );
		}
		else
		{
#ifdef MFEM_USE_SUITESPARSE
			// No SetReuseSymbolic() here, deliberately: the linear path factorises
			// once and this object is destroyed straight after, so retaining the
			// analysis would buy nothing and cost a copy of the pattern.
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

		// A new solution invalidates the old post-processing rather than being
		// silently paired with it.
		postProcessed = false;
	}

	void GradShafranovSolver::postProcess()
	{
		if ( !prepared )
			throw std::logic_error( "meq::GradShafranovSolver::postProcess: solve() has not been called" );

		// THIS USED TO REFUSE THE NON-LINEAR PATH, AND NO LONGER DOES.
		//
		// DarcyForm::ReconstructFluxAndPot() consulted only the LINEAR potential
		// mass form M_p, and meq's Newton path puts the whole potential block on
		// Mnl_p -- see buildForms() -- so the local problem got no potential mass
		// and no face constraint, was singular, and was factored and solved
		// anyway. What came back was not a degraded psi* but 1e15, without a word,
		// with psi_h agreeing to six figures either way.
		//
		// MFEM now lifts the non-linear potential integrators as a Jacobian frozen
		// at the computed potential, checking that the face constraint's gradient
		// does not depend on the trace -- which meq's constant tau satisfies, and
		// which is one more reason to keep it constant.
		//
		// MEASURED before this line was deleted, because a silent 1e15 is exactly
		// the failure a code read cannot detect. Example 5 on the Newton path,
		// L2( psi* ) over four dyadic meshes:
		//
		//     k = 1    rates 3.240, 3.101, 3.049      47x smaller than psi_h
		//     k = 2    rates 4.289, 4.114, 4.049     113x
		//     k = 3    rates 5.150, 5.062, 5.025     125x
		//
		// k+2 at every order, and buying the order it converges at rather than
		// merely converging. NewtonConvergence.cpp's
		// thePostProcessedPotentialSurvivesNewton is that measurement as an
		// assertion, and it is what would notice the defect returning.

		// DarcyForm::Reconstruct() builds the spaces on first use, so the second
		// call reuses them. The block vector handed in is the TWO-block view over
		// DarcyForm's own offsets, not the three-block solution: the trace arrives
		// separately as sol_r, and a three-block vector is the size mismatch
		// CLAUDE.md records against ReduceRHS().
		darcy->Reconstruct( darcySolution, traceX, totalFluxGf, enrichedFluxGf,
		                    postProcessedGf, enrichedTraceGf );

		// The same negation solve() applies to fluxGf, for the same reason: what
		// Reconstruct() writes is DarcyForm's -q. The total flux is deliberately
		// left alone -- it is defined by the constraint equation, which is written
		// in DarcyForm's convention throughout.
		enrichedFluxGf.Neg();

		/*
		 * AND IT IS NOT CHECKED HERE, WHICH IS A DECISION RATHER THAN AN
		 * OVERSIGHT.
		 *
		 * psi* is wrong on any element where dF/dpsi vanishes -- see the header,
		 * and NewtonConvergence.cpp's theReconstructionIsWrongWhereTheJacobianVanishes
		 * for the measurement. A version of this function did detect that, by
		 * comparing || psi* || against || psi_h ||, and it has been taken out: the
		 * condition it detects is a defect in a library meq does not own and is one
		 * flag test away from never arising, and a solver should not carry a
		 * standing defence against its dependency. The suite establishes the state
		 * of that defect; this code assumes the library works.
		 */
		postProcessed = true;
	}

	bool GradShafranovSolver::isPostProcessed() const
	{
		return postProcessed;
	}

	mfem::GridFunction &GradShafranovSolver::postProcessedPotential()
	{
		return postProcessedGf;
	}

	mfem::GridFunction const &GradShafranovSolver::postProcessedPotential() const
	{
		return postProcessedGf;
	}

	mfem::GridFunction &GradShafranovSolver::postProcessedFlux()
	{
		return enrichedFluxGf;
	}

	mfem::GridFunction const &GradShafranovSolver::postProcessedFlux() const
	{
		return enrichedFluxGf;
	}

	mfem::GridFunction &GradShafranovSolver::totalFlux()
	{
		return totalFluxGf;
	}

	mfem::GridFunction const &GradShafranovSolver::totalFlux() const
	{
		return totalFluxGf;
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

	int GradShafranovSolver::picardIterations() const
	{
		return picardIterationCount;
	}

	long GradShafranovSolver::symbolicFactorisations() const
	{
		return symbolicFactorisationCount;
	}

	long GradShafranovSolver::numericFactorisations() const
	{
		return numericFactorisationCount;
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

	double GradShafranovSolver::postProcessedPotentialError( mfem::Coefficient &exact ) const
	{
		if ( !postProcessed )
			throw std::logic_error( "meq::GradShafranovSolver::postProcessedPotentialError: postProcess() has not been called" );

		// orderValue + 1, because psi*_h lives one degree up and a rule chosen for
		// psi_h would be what limited the k+2 rate this exists to measure.
		mfem::IntegrationRule const *irs[ mfem::Geometry::NumGeom ];
		errorRules( orderValue + 1, irs );
		return postProcessedGf.ComputeL2Error( exact, irs );
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
