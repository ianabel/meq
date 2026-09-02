#include "GradShafranov.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
#include <omp.h>
#endif

// Whether this build has ANY direct trace solver. The fallback paths below are
// guarded on this rather than on MFEM_USE_SUITESPARSE alone, which was the same
// question only while UMFPack was the only choice.
#if defined( MFEM_USE_SUITESPARSE ) || defined( MFEM_USE_MKL_PARDISO ) \
	|| defined( MFEM_USE_CUDSS )
#define MEQ_HAVE_DIRECT_TRACE_SOLVER 1
#endif

namespace meq
{
namespace
{
	/*
	 * The assembly mode a fresh solver starts in: SERIAL, ALWAYS.
	 *
	 * A gate on omp_get_max_threads() was written, measured and REMOVED, and the
	 * measurement is why this function still exists rather than being a constant.
	 *
	 * In isolation the threaded assembly wins everywhere worth caring about --
	 * 1.15x to 1.33x from four threads up, and still 1.19x on a 512-element mesh.
	 * On that evidence defaulting to Threaded whenever more than one thread is
	 * available looks free. It is not. `HighBetaConvergence`, which assembles a
	 * small system many times over inside a bordered Newton, goes from 21.5 s to
	 * 39 s -- **1.8x SLOWER**, reproducibly, on the same binary.
	 *
	 * The distinguishing variable is not mesh size, which is what makes this
	 * undecidable from here: HighBeta's meshes are 128 and 512 elements, and 512
	 * is exactly where the isolated benchmark still showed a win. What differs is
	 * how OFTEN assembly is called relative to everything else. MFEM's threaded
	 * path forks a team and buffers element blocks per call, and a caller that
	 * assembles repeatedly pays that per call while a caller that assembles once
	 * amortises it. The solver cannot know which caller it has.
	 *
	 * So there is no safe automatic default, and the honest default is MFEM's
	 * own. setAssemblyMode( Threaded ) is an informed choice: worth taking on a
	 * large mesh assembled a few times, worth avoiding on a small one assembled
	 * hundreds of times.
	 *
	 * This is the same lesson as the trace matrix's symmetry and the threaded-MKL
	 * collapse: a property measured on the easy configuration is not a property
	 * of the code.
	 */
	GradShafranovSolver::AssemblyMode defaultAssemblyMode()
	{
		return GradShafranovSolver::AssemblyMode::Serial;
	}

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
	/*
	 * Build the chosen direct solver, configured as meq wants it.
	 *
	 * One place rather than four, which is the point: the four call sites --
	 * Picard, Newton, the bordered Newton and the linear path -- had four copies
	 * of the same UMFPack configuration, and the ONLY difference between them
	 * that ever mattered was whether the symbolic analysis is retained.
	 *
	 * `reuseSymbolic` is that difference and it is not cosmetic. Every path that
	 * re-solves with the same sparsity wants it: NewtonSolver::Mult calls
	 * SetOperator on this object once per iteration, and Picard runs 122 to 290
	 * full factorisations. The linear path does not, because it factorises once
	 * and destroys the object, so retaining the analysis would buy nothing and
	 * cost a copy of the pattern. All three packages compare the PATTERN rather
	 * than the object, so a matrix rebuilt into a fresh object with the same
	 * structure still hits the reuse.
	 */
	std::unique_ptr<mfem::Solver>
	makeTraceSolver( GradShafranovSolver::TraceSolver choice, bool reuseSymbolic )
	{
		switch ( choice )
		{
			case GradShafranovSolver::TraceSolver::UMFPack:
			{
#ifdef MFEM_USE_SUITESPARSE
				auto solver = std::make_unique<mfem::UMFPackSolver>();
				solver->Control[ UMFPACK_ORDERING ] = UMFPACK_ORDERING_METIS;
				if ( reuseSymbolic )
					solver->SetReuseSymbolic();
				return solver;
#else
				break;
#endif
			}
			case GradShafranovSolver::TraceSolver::Pardiso:
			{
#ifdef MFEM_USE_MKL_PARDISO
				auto solver = std::make_unique<mfem::PardisoSolver>();
				// STRUCTURE symmetric, not symmetric. The trace matrix is
				// symmetric to 2e-16 on a fitted mesh and asymmetric at 5.4e-1 on
				// the extension path, where HDGExtensionIntegrator deposits an
				// outer product into the flux block; the SPARSITY is symmetric on
				// both. So this is the type that is right for meq's headline
				// configuration as well as for the easy one.
				solver->SetMatrixType( mfem::PardisoSolver::REAL_STRUCTURE_SYMMETRIC );
				solver->SetPrintLevel( 0 );
				if ( reuseSymbolic )
					solver->SetReuseSymbolic();
				return solver;
#else
				break;
#endif
			}
			case GradShafranovSolver::TraceSolver::cuDSS:
			{
#ifdef MFEM_USE_CUDSS
				auto solver = std::make_unique<mfem::CuDSSSolver>();
				// NONSYMMETRIC + FULL for the reason above, and because cuDSS's
				// symmetric modes would be wrong on the extension path rather than
				// merely slower.
				solver->SetMatrixSymType( mfem::CuDSSSolver::NONSYMMETRIC );
				solver->SetMatrixViewType( mfem::CuDSSSolver::FULL );
				// cuDSS spells the same idea differently and requires it BEFORE
				// the first SetOperator -- it verifies against its own null handle.
				if ( reuseSymbolic )
					solver->SetReorderingReuse( true );
				return solver;
#else
				break;
#endif
			}
		}

		// Unreachable through setTraceSolver(), which refuses an unavailable
		// choice. Reachable only if that check and this switch disagree, which is
		// worth saying out loud rather than returning null into a dereference.
		throw std::logic_error(
			"meq: the requested trace solver is not available in this build, and "
			"setTraceSolver() should already have refused it" );
	}

	/// The factorisation counters, where the package keeps any. UMFPack and
	/// PARDISO both do and cuDSS does not, so a cuDSS solve reports zero rather
	/// than a wrong number.
	void readFactorisationCounts( mfem::Solver const &solver,
	                              long &symbolic, long &numeric )
	{
#ifdef MFEM_USE_SUITESPARSE
		if ( auto const *umf = dynamic_cast<mfem::UMFPackSolver const *>( &solver ) )
		{
			symbolic = umf->GetNumSymbolicFactorizations();
			numeric = umf->GetNumNumericFactorizations();
			return;
		}
#endif
#ifdef MFEM_USE_MKL_PARDISO
		if ( auto const *par = dynamic_cast<mfem::PardisoSolver const *>( &solver ) )
		{
			symbolic = par->GetNumSymbolicFactorizations();
			numeric = par->GetNumNumericFactorizations();
			return;
		}
#endif
		symbolic = 0;
		numeric = 0;
	}
#endif // MEQ_HAVE_DIRECT_TRACE_SOLVER
}


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
		  normalisedSource( nullptr ),
		  psiAxisValue( 0.0 ),
		  normalisationResidualValue( 0.0 ),
		  normalisationChoice( Normalisation::Coupled ),
		  boundaryData( nullptr ),
		  initialGuess( nullptr ),
		  globalisationChoice( Globalisation::None ),
		  localSolverChoice( LocalSolver::Newton ),
		  orderingChoice( NonlinearOrdering::NPC ),
		  assemblyModeChoice( defaultAssemblyMode() ),
		  traceSolverChoice( TraceSolver::UMFPack ),
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
		if ( linearSource )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: a source is already set; one solver holds one source" );

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
		// AND A SECOND SOURCE OF THE SAME KIND, which until SolverContract.cpp went
		// looking was the one way round "one solver holds one source": the checks
		// were symmetric ACROSS the two overloads and silent WITHIN each, so
		// setSource( a ); setSource( b ) replaced a with b without a word and the
		// solve answered a different question than the caller had asked.
		if ( nonlinearSource )
			throw std::logic_error( "meq::GradShafranovSolver::setSource: a source is already set; one solver holds one source" );

		// Nothing else to do here. Unlike the linear case there is no coefficient
		// to build: the source does not reach a right hand side at all, it becomes
		// a SourceIntegrator on the non-linear potential mass form, where Newton
		// can differentiate it. See buildForms().
		nonlinearSource = &fIn;
	}

	void GradShafranovSolver::setSource( NormalisedSource &fIn, double psiAxisGuessIn )
	{
		if ( !std::isfinite( psiAxisGuessIn ) || psiAxisGuessIn == 0.0 )
			throw std::invalid_argument( "meq::GradShafranovSolver::setSource: the psi_ax guess must be finite and non-zero" );
		/*
		 * NO ORDERING IS REFUSED HERE, AND THAT HAS BEEN SETTLED TWICE.
		 *
		 * A guard once refused MFEM's NLOrdering::LineariseThenCondense, on the
		 * reasoning that the reduced operator was a linearised residual between
		 * GetGradient() calls, so psi_ax -- which enters only through the source
		 * -- would be invisible to a finite difference of it and the border would
		 * come back zero. That was read out of a header summary rather than out
		 * of the code under it, and the code read the source afresh on every
		 * residual evaluation. Measured, both orderings reached the same psi_ax
		 * to every digit printed. The guard went, and then the mode did: upstream
		 * deleted it as a condensation in disguise.
		 *
		 * WHAT REPLACED IT MAKES THE QUESTION MOOT. Under NonlinearOrdering::NPC
		 * psi is an unknown of the system rather than a function of the trace, so
		 * two of the three bordered quantities stop being differences at all --
		 * the border row is exactly -e_j and the corner is exactly 1. Only
		 * c = dR/ds is differenced, in a SCALAR, and s reaches every element's
		 * source under either ordering. There is nothing left for a linearisation
		 * history to hide.
		 *
		 * The moral is the one this tree keeps relearning, and it is why the
		 * guard is recorded rather than merely deleted: it was written from a
		 * header comment rather than from the code under it.
		 */

		// The ordinary checks first, and through the ordinary overload, so that
		// "one solver holds one source" is enforced in exactly one place.
		setSource( static_cast<Source const &>( fIn ) );

		normalisedSource = &fIn;
		psiAxisValue = psiAxisGuessIn;
		normalisationResidualValue = 0.0;
	}

	double GradShafranovSolver::axisFlux( mfem::Vector const &trace, int *element )
	{
		if ( !normalisedSource )
			throw std::logic_error( "meq::GradShafranovSolver::axisFlux: psi_ax is not an unknown of this solver" );
		if ( !prepared )
			throw std::logic_error( "meq::GradShafranovSolver::axisFlux: prepare() has not been called" );
		return recoverPeak( trace, psiAxisValue, element );
	}

	void GradShafranovSolver::setNormalisationCoupling( Normalisation choice )
	{
		normalisationChoice = choice;
	}

	GradShafranovSolver::Normalisation GradShafranovSolver::normalisationCoupling() const
	{
		return normalisationChoice;
	}

	bool GradShafranovSolver::normalisationIsUnknown() const
	{
		return normalisedSource != nullptr;
	}

	double GradShafranovSolver::psiAxis() const
	{
		return psiAxisValue;
	}

	double GradShafranovSolver::normalisationResidual() const
	{
		return normalisationResidualValue;
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
			throw std::invalid_argument( "meq::GradShafranovSolver::setPicardDamping: the damping must be in ( 0, 1 ]" );
		picardDamping = damping;
	}

	void GradShafranovSolver::setAndersonDepth( int depth )
	{
		if ( depth < 0 )
			throw std::invalid_argument( "meq::GradShafranovSolver::setAndersonDepth: the depth cannot be negative" );
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

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
		// Held across calls rather than built per iteration, which is the whole
		// point: Picard runs 122 to 290 of these, each a full factorisation of a
		// matrix whose sparsity never changes. prepare() rebuilds `reduced` every
		// iteration, but the reuse compares the pattern rather than the object --
		// it documents accepting "a matrix rebuilt into a fresh object with the
		// same structure" -- so the analysis survives that.
		if ( !picardSolver )
			picardSolver = makeTraceSolver( traceSolverChoice, true );
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

	/*
	 * Threading the element-local assembly.
	 *
	 * MFEM aborts the process if asked for AssemblyMode::Threaded in a build
	 * without MFEM_USE_OPENMP or without MFEM_THREAD_SAFE, and it is right to --
	 * a caller asking for threads is asking a performance question, and quietly
	 * running the serial loop would report a speedup nobody got. But aborting is
	 * not a library's decision to impose on meq's callers, so the build is
	 * checked here and the refusal comes back as an exception like every other
	 * value fault in this class.
	 */
	void GradShafranovSolver::setAssemblyMode( AssemblyMode choice )
	{
		if ( choice == AssemblyMode::Threaded )
		{
#if !defined( MFEM_USE_OPENMP ) || !defined( MFEM_THREAD_SAFE )
			throw std::invalid_argument(
				"AssemblyMode::Threaded needs an MFEM built with both "
				"MFEM_USE_OPENMP and MFEM_THREAD_SAFE; this one has at least one "
				"of them off, and MFEM would abort the process rather than fall "
				"back to the serial loop" );
#endif
		}

		assemblyModeChoice = choice;
		built = false;
		prepared = false;
	}

	GradShafranovSolver::AssemblyMode
	GradShafranovSolver::assemblyMode() const
	{
		return assemblyModeChoice;
	}

	/*
	 * The trace solver, and which of them this build actually has.
	 *
	 * Kept as a compile-time question rather than a runtime one because that is
	 * what it is: MFEM either wraps the package or it does not, and there is no
	 * state in which asking would give a different answer later.
	 */
	bool GradShafranovSolver::traceSolverAvailable( TraceSolver choice )
	{
		switch ( choice )
		{
			case TraceSolver::UMFPack:
#ifdef MFEM_USE_SUITESPARSE
				return true;
#else
				return false;
#endif
			case TraceSolver::Pardiso:
#ifdef MFEM_USE_MKL_PARDISO
				return true;
#else
				return false;
#endif
			case TraceSolver::cuDSS:
#ifdef MFEM_USE_CUDSS
				return true;
#else
				return false;
#endif
		}
		return false;
	}

	void GradShafranovSolver::setTraceSolver( TraceSolver choice )
	{
		if ( !traceSolverAvailable( choice ) )
			throw std::invalid_argument(
				"meq::GradShafranovSolver::setTraceSolver: this MFEM was built "
				"without the package that solver needs -- UMFPack wants "
				"MFEM_USE_SUITESPARSE, Pardiso wants MFEM_USE_MKL_PARDISO, cuDSS "
				"wants MFEM_USE_CUDSS. Refused rather than silently substituted, "
				"because a caller naming a solver has a reason for naming it; "
				"traceSolverAvailable() answers the question without throwing" );

		traceSolverChoice = choice;
		// NOT a rebuild. The trace solver is chosen when the reduced system is
		// solved, not when the forms are assembled, so `built` and `prepared`
		// both stay valid -- unlike every other setter in this block.
	}

	GradShafranovSolver::TraceSolver
	GradShafranovSolver::traceSolver() const
	{
		return traceSolverChoice;
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

	std::unique_ptr<mfem::Coefficient>
		GradShafranovSolver::transferredDatum( mfem::PositionFunction g )
	{
		if ( !transferPath )
			throw std::logic_error(
				"GradShafranovSolver::transferredDatum: there is no transferred "
				"datum on the fitted path -- setExtension() was never called, so "
				"the trace unknown on the boundary IS the condition imposed" );

		// The fixed-boundary problem puts psi = 0 on Gamma, and the extension
		// benchmark shifts its solution so that its datum is homogeneous too, so
		// this default is every case in the suite. It is still a parameter,
		// because g is a property of the problem and not of the technique.
		if ( !g )
			g = []( mfem::Vector const & ) { return 0.0; };

		// darcyFlux, NOT flux(): the same convention HDGExtensionIntegrator was
		// assembled against. radius and extensionLineOrder likewise have to be
		// the ones buildForms() gave the integrator -- a different rule along the
		// path is a different lifting, which the MFEM header says in as many
		// words.
		return std::make_unique<mfem::TransferredDatumCoefficient>(
			*transferPath, std::move( g ), darcyFlux, radius, extensionLineOrder );
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

		// Who runs the element loop. OUTSIDE the branch below, deliberately:
		// ComputeH() factors A, forms and factors the Schur complement and does
		// one local back-substitution per trace dof on BOTH paths, so the linear
		// solve has exactly as much element-local work to thread as the Newton
		// one. setAssemblyMode() has already refused Threaded if the build cannot
		// honour it, so this cannot reach MFEM's abort.
		darcy->GetHybridization()->SetAssemblyMode(
			assemblyModeChoice == AssemblyMode::Threaded
				? mfem::DarcyHybridization::AssemblyMode::Threaded
				: mfem::DarcyHybridization::AssemblyMode::Serial );

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

			// AND NEITHER OF THE TWO ABOVE MEANS ANYTHING UNDER
			// NonlinearOrdering::NPC, which is meq's default. They are set
			// regardless, because they are properties of the hybridization rather
			// than of the solve, and because CondenseThenLinearise is one
			// setNonlinearOrdering() call away and must find them configured.
			//
			// There is nothing to select here: NPC is NOT A MODE OF THIS OBJECT.
			// MFEM used to offer SetNonlinearOrdering() with a third value that
			// claimed to be NPC and was a condensation in disguise; it is deleted,
			// and NPC is a separate mfem::DarcyNPCOperator over the full
			// ( q, psi, psihat ) vector, built in solve(). So the ordering choice
			// changes which OPERATOR the outer iteration drives, not how this one
			// is configured.
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

		formSystem();
	}

	void GradShafranovSolver::formSystem()
	{
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

	double GradShafranovSolver::recoverPeak( mfem::Vector const &trace, double psiAxisIn,
	                                         int *element, int *dof )
	{
		normalisedSource->setNormalisation( psiAxisIn );

		// Into scratch, not into the solution blocks: this is called from inside
		// finite differences, and leaving the caller's psi_h perturbed would be a
		// silent corruption of exactly the field the answer is read from.
		//
		// What the element-local Newtons inside ComputeSolution() start from is
		// NOT this vector -- DarcyHybridization keeps a copy taken at
		// FormLinearSystem() time and every local solve begins there. That is why
		// solveWithNormalisation() re-forms the system once per Newton step; see
		// the comment on the seed there, which is the single thing that made
		// these differences mean anything.
		darcy->RecoverFEMSolution( trace, darcyRhs, recoveryScratch );

		/*
		 * psi_ax IS THE LARGEST NODAL VALUE, not the largest value of the
		 * polynomial, and that is a definition rather than an approximation.
		 *
		 * The two differ by O( h^(k+1) ) and both converge to max psi, so either
		 * would do for the physics. The nodal one is chosen because it is what
		 * makes the constraint DIFFERENTIABLE in a form the border can use: it is
		 * one entry of the recovered potential, so d psi_ax/d lambda is supported
		 * on the trace dofs of a single element. The maximum over an element's
		 * interior would move that support to wherever the interior maximum is,
		 * and the derivative would acquire the argmax's own sensitivity.
		 *
		 * Which element attained it is returned rather than searched for again,
		 * because the border needs it and a second search over a different
		 * recovery could find a different element.
		 */
		mfem::Vector const &potential = recoveryScratch.GetBlock( 1 );
		double best = -std::numeric_limits<double>::infinity();
		int bestElement = -1;
		int bestDof = -1;

		mfem::Array<int> dofs;
		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			potentialFes->GetElementDofs( e, dofs );
			for ( int i = 0; i < dofs.Size(); ++i )
			{
				int const index = dofs[ i ];
				if ( potential( index ) > best )
				{
					best = potential( index );
					bestElement = e;
					bestDof = index;
				}
			}
		}

		if ( element )
			*element = bestElement;
		if ( dof )
			*dof = bestDof;
		return best;
	}

	void GradShafranovSolver::traceDofsOfElement( int element, mfem::Array<int> &dofs ) const
	{
		dofs.SetSize( 0 );
		if ( element < 0 )
			return;

		mfem::Array<int> faces, orientations, faceDofs;
		// Two dimensional throughout -- the constructor refuses anything else --
		// so an element's faces are its edges.
		mesh.GetElementEdges( element, faces, orientations );

		for ( int i = 0; i < faces.Size(); ++i )
		{
			traceFes->GetFaceVDofs( faces[ i ], faceDofs );
			for ( int j = 0; j < faceDofs.Size(); ++j )
			{
				// A negative index is MFEM's sign-carrying encoding; the trace
				// space here is scalar, so only the index is wanted.
				int const d = faceDofs[ j ] >= 0 ? faceDofs[ j ] : -1 - faceDofs[ j ];
				dofs.Append( d );
			}
		}

		dofs.Sort();
		dofs.Unique();
	}

	/*
	 * THE BORDERED NEWTON: the trace and psi_ax solved together.
	 *
	 * The system is
	 *
	 *     R( lambda, s ) = 0        the hybridized trace residual, with the
	 *                               source normalised by s
	 *     G( lambda, s ) = s - max psi_h( lambda, s ) = 0
	 *
	 * and the step comes from
	 *
	 *     [  A   c  ] [ dlambda ]     [ R ]
	 *     [  b^T d  ] [   ds    ]  = -[ G ]
	 *
	 * solved by block elimination: A y = R, A z = c, then
	 *
	 *     ds      = ( b.y - G ) / ( d - b.z )
	 *     dlambda = -y - z ds.
	 *
	 * ONE FACTORISATION AND TWO BACKSOLVES. That is the whole cost of the extra
	 * unknown on the linear-algebra side, and it is why the border is not
	 * assembled into an ( n + 1 ) matrix: a dense row and a dense column would
	 * cost fill in the factorisation for no gain.
	 *
	 * WHERE c AND b COME FROM, and why they are differenced rather than
	 * assembled. Both are derivatives of the CONDENSED residual. Assembling them
	 * would need the sensitivity of the element-local eliminations -- for c the
	 * derivative of each local solve with respect to a parameter of its own
	 * source, for b the derivative of the recovered potential with respect to the
	 * trace -- and DarcyHybridization exposes neither. So they are obtained by
	 * differencing the assembled residual, which is the same principle CEDRES++
	 * states for the local term (refs/CEDRES.pdf): differentiate the DISCRETE
	 * residual, never the continuous equation, because the continuous formula
	 * "seems to blow up if psi reaches a critical point" -- which is precisely
	 * the point psi_ax is defined at.
	 *
	 * c IS DENSE AND b IS NOT, and the asymmetry is structural rather than a
	 * saving. s enters every element's source, so dR/ds has an entry on every
	 * trace dof: one central difference in a scalar, two residual evaluations,
	 * done. max psi_h is one nodal value of one element, and under hybridization
	 * that element's recovered potential depends only on the trace dofs of its
	 * own faces -- so b has 3( k + 1 ) entries at most and the rest are exactly
	 * zero. That claim is asserted rather than assumed, in
	 * HighBetaConvergence.cpp's theAxisSensitivityIsLocalToItsElement.
	 *
	 * THE ORDER OF OPERATIONS IS LOAD BEARING. Every finite difference here runs
	 * ComputeSolution() through recoverPeak(), which refreshes the factored local
	 * Jacobians DarcyHybridization keeps; GetGradient() must therefore be the
	 * LAST thing called before the linear solve, or the matrix handed to UMFPACK
	 * belongs to a trace that has since been perturbed and put back.
	 *
	 * THE PRINTED RESIDUAL IS ||( R, gamma G )||, with gamma frozen at ||c|| from
	 * the first iterate. G is a flux and R is a trace residual, so the two cannot
	 * simply be concatenated; ||c|| is the factor that converts a perturbation of
	 * psi_ax into the units R is measured in, which is exactly the conversion
	 * wanted. Freezing it keeps the history a comparison of like with like -- a
	 * gamma recomputed each step would put the Jacobian's own variation into the
	 * convergence history and manufacture orders out of it.
	 */
	void GradShafranovSolver::solveWithNormalisation()
	{
		if ( globalisationChoice != Globalisation::None )
			throw std::logic_error( "meq::GradShafranovSolver::solve: psi_ax as an unknown is implemented for Globalisation::None only -- the KINSOL paths drive a residual of their own and the Picard ones do not build a Jacobian at all" );
		// No ordering guard: both surviving orderings carry psi_ax. They carry it
		// DIFFERENTLY, and the difference is the whole content of this function --
		// see the border and the corner below.

		bool const npcOrdering = orderingChoice == NonlinearOrdering::NPC;

		/*
		 * THE UNKNOWN, AND WHY ITS LENGTH DECIDES EVERYTHING ELSE HERE.
		 *
		 * Under NonlinearOrdering::NPC it is the whole ( q, psi, psihat ) vector
		 * -- `solution` itself -- and psi is therefore an INDEPENDENT unknown.
		 * Under CondenseThenLinearise it is the trace alone and psi is a function
		 * of it, recovered by an element-local non-linear solve.
		 *
		 * That one difference collapses two of the three differenced quantities
		 * in this bordered system:
		 *
		 *   b = -d( max psi_h )/d( unknown )   NPC: EXACTLY -e_j, a unit vector,
		 *                                      because max psi_h is literally one
		 *                                      entry of the unknown. Condensation:
		 *                                      3( k + 1 ) central differences over
		 *                                      the trace dofs of one element.
		 *
		 *   d = dG/ds                          NPC: EXACTLY 1, because
		 *                                      G = s - max psi( x ) and s does not
		 *                                      appear in psi. Condensation: 1
		 *                                      minus a central difference, because
		 *                                      s enters every element's source and
		 *                                      so moves the recovered psi.
		 *
		 * Only c = dR/ds is still differenced, and it has to be on both paths: s
		 * enters every element's source, so dR/ds is dense and there is no
		 * assembled route to it. Two residual evaluations, one central difference
		 * in a SCALAR.
		 *
		 * AND THE WHOLE FROZEN-SEED APPARATUS GOES AWAY UNDER NPC. The re-forming
		 * of the system after every accepted step -- the single thing that made
		 * these differences derivatives rather than noise under the condensation,
		 * see the note on it below -- exists because DarcyHybridization freezes
		 * the element-local Newton's initial guess at FormLinearSystem() time.
		 * NPC has no element-local non-linear solve at all, so there is no seed,
		 * nothing to go stale, and nothing to re-form.
		 */
		mfem::Vector &unknown = npcOrdering ? solution : traceX;
		int const n = unknown.Size();

		recoveryScratch.Update( darcy->GetOffsets() );

		newtonResidualHistory.clear();
		newtonIterationCount = 0;
		symbolicFactorisationCount = 0;
		numericFactorisationCount = 0;

		std::unique_ptr<mfem::DarcyNPCOperator> npc;
		if ( npcOrdering )
			npc = std::make_unique<mfem::DarcyNPCOperator>(
				*darcy->GetHybridization(), blockOffsets, darcyRhs );

		mfem::Vector residual( n ), column( n ), y( n ), z( n ), scratch( n );

		auto fieldResidual = [ & ]( mfem::Vector const &state, double normalisation,
		                            mfem::Vector &out )
		{
			normalisedSource->setNormalisation( normalisation );
			if ( npcOrdering )
			{
				npc->Mult( state, out );
				return;
			}
			// Refetched rather than held: formSystem() replaces the handle every
			// time the local seed is refreshed, and a reference taken before the
			// loop would outlive the operator it names.
			reduced.Ptr()->Mult( state, out );
			out -= traceB;
		};

		/*
		 * max psi_h, AND WHERE IT IS ATTAINED.
		 *
		 * Under NPC this is a READ rather than a recovery: the potential block of
		 * the unknown holds every nodal value, so the peak is one scan and the
		 * index it is attained at is the whole of the border row. Under the
		 * condensation it is recoverPeak(), which runs ComputeSolution() to
		 * rebuild psi from the trace and reports which ELEMENT attained it,
		 * because the border there is supported on that element's trace dofs.
		 *
		 * The normalisation is set on both paths even though NPC does not need it
		 * for the read, so that the two agree on the source's state afterwards --
		 * and so that a psi_ax the source refuses still throws from here, which is
		 * what the line search below is catching.
		 */
		auto peakAt = [ & ]( mfem::Vector const &state, double normalisation,
		                     int *element, int *dof )
		{
			if ( !npcOrdering )
				return recoverPeak( state, normalisation, element, dof );

			normalisedSource->setNormalisation( normalisation );

			double best = -std::numeric_limits<double>::infinity();
			int bestIndex = -1;
			for ( int i = blockOffsets[ 1 ]; i < blockOffsets[ 2 ]; ++i )
			{
				if ( state( i ) > best )
				{
					best = state( i );
					bestIndex = i;
				}
			}

			// There is no element to report: the index is into the FULL vector
			// and is all the border needs.
			if ( element )
				*element = -1;
			if ( dof )
				*dof = bestIndex;
			return best;
		};

		// A central difference is worth the second evaluation here. The residual
		// carries the element-local solves' own stopping tolerance as noise, so a
		// forward difference would leave an O( h ) truncation error on top of it
		// and the border would be the thing that limits Newton's order. Under NPC
		// there are no element-local solves and so no such noise, and the central
		// difference is kept anyway: it is the same quantity, measured the same
		// way, which is what lets the two paths be compared.
		auto normalisationStep = []( double value )
		{
			return 1.0e-5*std::max( std::abs( value ), 1.0e-10 );
		};

		auto sourceColumn = [ & ]( mfem::Vector const &state, double normalisation,
		                           mfem::Vector &out )
		{
			if ( normalisationChoice == Normalisation::Decoupled )
			{
				out = 0.0;
				return;
			}
			double const h = normalisationStep( normalisation );
			fieldResidual( state, normalisation + h, out );
			fieldResidual( state, normalisation - h, scratch );
			out -= scratch;
			out /= 2.0*h;
		};

		// Sized by the TRACE on both paths, because that is what it indexes;
		// under NPC the trace block starts at blockOffsets[ 2 ] of the unknown.
		mfem::Array<int> const &essentialTrace =
			darcy->GetHybridization()->GetEssentialTrueDofs();
		std::vector<bool> essential( traceX.Size(), false );
		for ( int i = 0; i < essentialTrace.Size(); ++i )
			essential[ essentialTrace[ i ] ] = true;

		double s = psiAxisValue;
		int argElement = -1;
		int argDof = -1;
		double peak = peakAt( unknown, s, &argElement, &argDof );
		double constraint = s - peak;
		fieldResidual( unknown, s, residual );

		// c at the starting iterate, computed whatever the coupling: it is both
		// the first step's column and the scale gamma. gamma converts a
		// perturbation of psi_ax into the units the trace residual is measured
		// in, and || c || is exactly that conversion. Normalisation::Decoupled
		// does not use the column but is given the same gamma, or the two
		// convergence histories would not be comparable -- which is the whole
		// point of having a control.
		mfem::Vector initialColumn( n );
		{
			double const h = normalisationStep( s );
			fieldResidual( unknown, s + h, initialColumn );
			fieldResidual( unknown, s - h, scratch );
			initialColumn -= scratch;
			initialColumn /= 2.0*h;
		}

		double gamma = initialColumn.Norml2();
		if ( !( gamma > 0.0 ) || !std::isfinite( gamma ) )
			gamma = 1.0;

		if ( normalisationChoice == Normalisation::Coupled )
			column = initialColumn;
		else
			column = 0.0;

		/*
		 * The convergence target is taken from the COLD iterate, for the reason
		 * the plain Newton path takes it from there: MFEM's rule scales the target
		 * by the residual at the iterate it was handed, so a good guess shrinks
		 * the target with it and past a point drives it under the round-off floor.
		 * On this path a guess is not an optimisation but part of the problem
		 * statement -- see setSource( NormalisedSource &, double ) -- so the
		 * effect would be systematic rather than occasional.
		 *
		 * COLD means the Dirichlet datum and nothing else, which under NPC is the
		 * flux and potential blocks zeroed as well as the free trace dofs.
		 */
		double reference = 0.0;
		{
			mfem::Vector coldState( unknown );
			int const traceBase = npcOrdering ? blockOffsets[ 2 ] : 0;
			for ( int i = 0; i < traceBase; ++i )
				coldState( i ) = 0.0;
			for ( int i = traceBase; i < n; ++i )
				if ( !essential[ i - traceBase ] )
					coldState( i ) = 0.0;

			mfem::Vector coldResidual( n );
			fieldResidual( coldState, s, coldResidual );
			double const coldPeak = peakAt( coldState, s, nullptr, nullptr );
			reference = std::hypot( coldResidual.Norml2(), gamma*( s - coldPeak ) );
		}
		double const target = std::max( newtonAbsoluteTolerance,
		                                newtonRelativeTolerance*reference );

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
		// The sparsity of the trace system does not change between Newton steps
		// and this object outlives the loop, so the analysis is done once.
		std::unique_ptr<mfem::Solver> const linearOwned =
			makeTraceSolver( traceSolverChoice, true );
		mfem::Solver &linear = *linearOwned;
#else
		mfem::GMRESSolver linear;
		linear.SetRelTol( 1.0e-14 );
		linear.SetAbsTol( 0.0 );
		linear.SetMaxIter( 5000 );
		linear.SetPrintLevel( -1 );
#endif

		// Under NPC the border is solved against the FULL Jacobian, and A^-1 is
		// the hybridized elimination -- reduce to the trace, solve there with
		// `linear`, recover the local increments. So the two backsolves the
		// border costs are still two TRACE solves against one factorisation, and
		// the extra unknown still costs one factorisation and two backsolves.
		mfem::DarcyNPCSolver npcLinear( linear );

		bool converged = false;
		for ( int iteration = 0; iteration <= newtonMaxIterations; ++iteration )
		{
			double const norm = std::hypot( residual.Norml2(), gamma*constraint );
			newtonResidualHistory.push_back( norm );
			newtonIterationCount = iteration;

			if ( norm <= target )
			{
				converged = true;
				break;
			}
			if ( !std::isfinite( norm ) || iteration == newtonMaxIterations )
				break;

			bool const coupled = normalisationChoice == Normalisation::Coupled;

			// d = 1 - d( max psi_h )/d psi_ax, the corner of the border.
			//
			// EXACTLY 1 UNDER NPC and not differenced: psi is an unknown of the
			// system, so moving s moves the RESIDUAL and not psi, and G = s - psi_j
			// differentiates to 1 in s. Under the condensation psi is a function of
			// s through every element's source, so the derivative is real and has
			// to be measured.
			double corner = 1.0;
			if ( coupled && !npcOrdering )
			{
				double const h = normalisationStep( s );
				double const peakUp = recoverPeak( traceX, s + h );
				double const peakDown = recoverPeak( traceX, s - h );
				corner = 1.0 - ( peakUp - peakDown )/( 2.0*h );
			}

			// b = -d( max psi_h )/d( unknown ). Essential dofs are left at zero:
			// the step does not move them, so what they would contribute is
			// multiplied by an increment that is identically zero.
			mfem::Array<int> borderDofs;
			mfem::Vector border;

			if ( npcOrdering )
			{
				// ONE ENTRY, EXACT, NOT DIFFERENCED. max psi_h is the argDof'th
				// entry of the unknown, so d( max psi_h )/d( unknown ) is the unit
				// vector e_argDof and b is its negation. Nothing is measured, so
				// nothing here carries a truncation error, and the 3( k + 1 )
				// recoveries the condensation spends on this are not spent.
				borderDofs.SetSize( 1 );
				borderDofs[ 0 ] = argDof;
				border.SetSize( 1 );
				border( 0 ) = coupled ? -1.0 : 0.0;
			}
			else
			{
				// On the trace dofs of the element that attained the maximum and
				// nowhere else -- under hybridization that element's recovered
				// potential depends on no others. Asserted rather than assumed, in
				// HighBetaConvergence.cpp's theAxisSensitivityIsLocalToItsElement.
				traceDofsOfElement( argElement, borderDofs );
				border.SetSize( borderDofs.Size() );
				border = 0.0;

				double const traceStep = 1.0e-6*std::max( traceX.Normlinf(), 1.0 );
				for ( int i = 0; coupled && i < borderDofs.Size(); ++i )
				{
					int const dof = borderDofs[ i ];
					if ( essential[ dof ] )
						continue;

					double const saved = traceX( dof );
					traceX( dof ) = saved + traceStep;
					double const up = recoverPeak( traceX, s );
					traceX( dof ) = saved - traceStep;
					double const down = recoverPeak( traceX, s );
					traceX( dof ) = saved;
					border( i ) = -( up - down )/( 2.0*traceStep );
				}
			}

			// LAST, and after the source is put back to s: every difference above
			// ran a residual evaluation of its own, and the gradient is what leaves
			// the factored local blocks the backsolves below are entitled to. That
			// ordering is load bearing under the condensation, where the
			// differences run local solves that overwrite exactly those blocks; it
			// is kept under NPC because it is the right thing to write either way.
			normalisedSource->setNormalisation( s );

			if ( npcOrdering )
			{
				mfem::Operator &jacobian = npc->GetGradient( unknown );
				npcLinear.SetOperator( jacobian );
				npcLinear.Mult( residual, y );
				npcLinear.Mult( column, z );
			}
			else
			{
				mfem::Operator &gradient = reduced.Ptr()->GetGradient( traceX );
				linear.SetOperator( gradient );
				linear.Mult( residual, y );
				linear.Mult( column, z );
			}

			double borderDotY = 0.0;
			double borderDotZ = 0.0;
			for ( int i = 0; i < borderDofs.Size(); ++i )
			{
				borderDotY += border( i )*y( borderDofs[ i ] );
				borderDotZ += border( i )*z( borderDofs[ i ] );
			}

			double const denominator = corner - borderDotZ;
			if ( denominator == 0.0 || !std::isfinite( denominator ) )
				throw std::runtime_error( "meq::GradShafranovSolver::solve: the bordered Jacobian is singular in psi_ax -- the normalisation has no influence on the solution it normalises" );

			double const deltaS = ( borderDotY - constraint )/denominator;

			/*
			 * BACKTRACKING, AND IT IS NOT OPTIONAL HERE.
			 *
			 * Measured, at nu = 4 and a pressure amplitude of 10 on the standard
			 * box: the full step converges for the mild profiles and wanders for
			 * this one, the augmented residual reading 2.7e-1, 6.6e-2, 2.0e0,
			 * 7.7e-1, 2.5e2, 6.7e4, 7.7e6, 3.8e8 before psi_ax crosses zero and
			 * the source refuses the normalisation. That is not the Jacobian
			 * being wrong -- the same Jacobian finishes the milder cases at
			 * observed order 2 -- it is the equilibrium being a MOUNTAIN-PASS
			 * solution of a superlinear problem, where the linearised operator is
			 * indefinite and an undamped step happily leaves the basin.
			 *
			 * Halving on the augmented norm, and the best trial kept if none of
			 * them improves it, so a bad step costs iterations rather than the
			 * solve. The trial evaluation is not wasted: whichever step is
			 * accepted, its residual and constraint are the ones the next
			 * iteration starts from.
			 *
			 * Under NPC the damping scales the fields and the trace TOGETHER,
			 * because they are one vector -- which is the half of a line search a
			 * trace-only operator cannot express.
			 */
			mfem::Vector const savedState( unknown );
			double const savedS = s;

			double bestNorm = std::numeric_limits<double>::infinity();
			double bestDamping = 0.0;
			double damping = 1.0;
			bool accepted = false;

			for ( int trial = 0; trial < 12 && !accepted; ++trial, damping *= 0.5 )
			{
				double trialNorm = std::numeric_limits<double>::infinity();
				try
				{
					unknown = savedState;
					unknown.Add( -damping, y );
					unknown.Add( -damping*deltaS, z );
					s = savedS + damping*deltaS;

					peak = peakAt( unknown, s, &argElement, &argDof );
					constraint = s - peak;
					fieldResidual( unknown, s, residual );
					trialNorm = std::hypot( residual.Norml2(), gamma*constraint );
				}
				catch ( std::exception const & )
				{
					// A normalisation the source will not accept -- psi_ax through
					// zero, most often -- is a step that left the branch. Reject
					// it like any other non-improving step rather than letting it
					// end the solve.
					trialNorm = std::numeric_limits<double>::infinity();
				}

				if ( std::isfinite( trialNorm ) && trialNorm < bestNorm )
				{
					bestNorm = trialNorm;
					bestDamping = damping;
				}
				// Armijo, with the mildest useful constant: what is wanted is a
				// step that does not make things worse, not an optimal one.
				accepted = std::isfinite( trialNorm )
				           && trialNorm < ( 1.0 - 1.0e-4*damping )*norm;
			}

			if ( !accepted )
			{
				if ( bestDamping == 0.0 )
					throw std::runtime_error( "meq::GradShafranovSolver::solve: no damping of the bordered Newton step gave a finite residual -- psi_ax through zero, most often, which is the branch leaving the physical one" );
				unknown = savedState;
				unknown.Add( -bestDamping, y );
				unknown.Add( -bestDamping*deltaS, z );
				s = savedS + bestDamping*deltaS;
				peak = peakAt( unknown, s, &argElement, &argDof );
				constraint = s - peak;
				fieldResidual( unknown, s, residual );
			}

			/*
			 * THE STEP IS SETTLED, SO THE ELEMENT-LOCAL SOLVES ARE GIVEN A FRESH
			 * STARTING POINT -- and under the condensation this is the single
			 * thing that makes the border a derivative rather than noise.
			 *
			 * DarcyHybridization takes its local initial guess from the solution
			 * blocks at FormLinearSystem() time and keeps it for the life of the
			 * reduced system. Left alone, every local Newton in every residual
			 * evaluation restarts from the ORIGINAL guess, however far the trace
			 * has since travelled. Measured on nu = 4 at amplitude 10, n = 16,
			 * k = 2: 40,000 to 60,000 element-local iterations per outer step,
			 * most of them hitting the cap of 100 -- and a local solve that ran
			 * out of iterations returns whatever it had reached, which is not a
			 * function of anything. Differencing psi_ax by 9e-6 then moved
			 * max psi_h from 0.8961 to 2.04 and on the next step to 3.84. The
			 * corner of the border read 1.6e5 where it should read about 1, the
			 * step in psi_ax collapsed to 1e-8 against a constraint residual of
			 * 3e-3, and the iteration stalled -- looking exactly like a singular
			 * border and being nothing of the kind.
			 *
			 * Re-forming the system from the recovered state puts every local
			 * solve within one or two iterations of its answer, so it converges,
			 * so it is continuous in psi_ax, so the difference is a derivative.
			 * It also removes the cost: those tens of thousands of local
			 * iterations were the run time.
			 *
			 * The right hand side is rebuilt from zero exactly as prepare() does,
			 * and traceB comes back the same -- the essential trace values have
			 * not moved and the source is on the operator, not the right hand
			 * side -- so the residual being differenced does not drift.
			 *
			 * NONE OF WHICH APPLIES UNDER NPC, and that is the point rather than
			 * an omission: there is no element-local non-linear solve, so there is
			 * no seed to freeze, no local iteration count to blow up, and nothing
			 * for a stale linearisation to corrupt. q and psi are Newton state and
			 * are already in `unknown`.
			 */
			if ( !npcOrdering )
			{
				darcySolution = recoveryScratch;
				rhs = 0.0;
				formSystem();
			}

			sourceColumn( unknown, s, column );
		}

		psiAxisValue = s;
		normalisationResidualValue = constraint;
		normalisedSource->setNormalisation( s );

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
		readFactorisationCounts( linear, symbolicFactorisationCount,
		                         numericFactorisationCount );
#endif

		if ( !converged )
			throw std::runtime_error( "meq::GradShafranovSolver::solve: the non-linear iteration did not converge" );
	}

	void GradShafranovSolver::solve()
	{
		// The Picard paths iterate a fixed point on the POTENTIAL, not a residual
		// on the trace, so they do not go through prepare()-then-Newton at all --
		// picardStep() re-enters prepare() itself, once per iteration.
		// Stale on a solver reused across globalisations otherwise; the handoff
		// republishes it after stage 2.
		picardIterationCount = 0;

		// Before the Picard dispatches below, not after: a normalised source with
		// a Picard globalisation would otherwise be routed into a fixed point on
		// the potential that has no idea psi_ax is an unknown, and would converge
		// to the solution of a different problem.
		if ( normalisedSource && globalisationChoice != Globalisation::None )
			throw std::logic_error( "meq::GradShafranovSolver::solve: psi_ax as an unknown is implemented for Globalisation::None only -- the KINSOL paths drive a residual of their own and the Picard ones build no Jacobian to border" );

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

		/*
		 * WHETHER q AND psi COME BACK IN `solution` OR HAVE TO BE REBUILT.
		 *
		 * Under NonlinearOrdering::NPC they are Newton state: the outer unknown
		 * is the whole ( q, psi, psihat ) vector, which IS `solution`, so they
		 * are already there when the iteration returns and there is nothing to
		 * recover. Under the condensation the unknown is the trace alone and the
		 * fields are a function of it, so RecoverFEMSolution() is what produces
		 * them; the linear path is the same.
		 *
		 * Skipping it under NPC is a correctness decision and not a saving.
		 * RecoverFEMSolution() runs ComputeSolution(), whose element-local solves
		 * are NON-LINEAR under this discretisation, and upstream records that
		 * function as never having been exercised against an NPC solution. Asking
		 * an unchecked back-substitution to reproduce fields meq already holds
		 * exactly is a way to lose them, not a way to confirm them.
		 */
		bool const fieldsAreState = usesNonlinearForms()
		                            && orderingChoice == NonlinearOrdering::NPC;

		// psi_ax has to be in the source before anything assembles, so that
		// prepare() and the first residual evaluation see the same normalisation.
		if ( normalisedSource )
			normalisedSource->setNormalisation( psiAxisValue );

		prepare();

		if ( normalisedSource )
		{
			solveWithNormalisation();
		}
		else if ( nonlinearSource )
		{
			bool const npcOrdering = orderingChoice == NonlinearOrdering::NPC;

			/*
			 * THE COLD REFERENCE IS TAKEN FIRST, AND BY PREPARING WITHOUT THE
			 * GUESS RATHER THAN BY EDITING THE TRACE. Both halves of that are
			 * bug fixes; see the note where the tolerance is set for the
			 * measurement, and note that it has to happen HERE, before anything
			 * binds a reference to the reduced operator, because the second
			 * prepare() replaces it.
			 *
			 * IT IS A NORM OF A DIFFERENT VECTOR UNDER THE TWO ORDERINGS, and
			 * that is correct rather than something to reconcile. NPC's residual
			 * is the full ( q, psi, psihat ) system and the condensation's is the
			 * trace alone, so the two reference values are not comparable -- but
			 * neither is either one's residual history, and the reference exists
			 * only to set a target in the same units as the history it gates.
			 * Comparing the two paths means comparing the SOLUTIONS, which is
			 * what the suite does.
			 */
			auto coldNorm = [ & ]()
			{
				if ( npcOrdering )
				{
					mfem::DarcyNPCOperator cold( *darcy->GetHybridization(),
					                             blockOffsets, darcyRhs );
					mfem::Vector coldResidual( cold.Height() );
					cold.Mult( solution, coldResidual );
					return coldResidual.Norml2();
				}

				mfem::Vector coldResidual( traceX.Size() );
				reduced.Ptr()->Mult( traceX, coldResidual );
				coldResidual -= traceB;
				return coldResidual.Norml2();
			};

			double coldReference = -1.0;
			if ( initialGuess )
			{
				mfem::Coefficient *const guess = initialGuess;
				initialGuess = nullptr;
				prepare();

				coldReference = coldNorm();

				initialGuess = guess;
				prepare();
			}

			ResidualRecorder recorder( newtonResidualHistory );

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
			// NewtonSolver::Mult calls prec->SetOperator( *grad ) on THIS object
			// once per iteration, and the trace system's sparsity does not change
			// between them -- so without the reuse the ordering is recomputed and
			// thrown away every step, at a fifth to a quarter of the factorisation.
			// The pattern is compared entry by entry rather than assumed, so a
			// pattern that did change is re-analysed and the answer is unaffected.
			std::unique_ptr<mfem::Solver> const linearOwned =
				makeTraceSolver( traceSolverChoice, true );
			mfem::Solver &linear = *linearOwned;
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
			/*
			 * NPC: THE UNKNOWN IS THE WHOLE SYSTEM, AND meq ALREADY HAD IT.
			 *
			 * mfem::DarcyNPCOperator is an Operator on ( q, psi, psihat )
			 * together -- Nguyen, Peraire & Cockburn eqs (14)-(18) -- with the
			 * Jacobian solved by hybridized elimination inside
			 * mfem::DarcyNPCSolver: reduce to the trace, solve there once, recover
			 * the local increments. Every element-local operation is ONE linear
			 * solve against ONE factorisation, and GetNumLocalNLIterations()
			 * staying at zero is the acceptance signal that it really is NPC.
			 *
			 * It costs meq almost nothing to give the fields to Newton, because
			 * `solution` has always been a three-block vector on `blockOffsets`
			 * with darcyFlux, potentialGf and traceGf MakeRef'd into it. That IS
			 * the NPC unknown, block for block. So the fields are already in
			 * place when Mult() returns and RecoverFEMSolution() drops out of this
			 * path entirely -- see the guard on it below, which is not an
			 * optimisation: upstream records ComputeSolution() as never having
			 * been exercised against an NPC solution, and running the
			 * condensation's element-local NON-LINEAR back-substitution over an
			 * answer that already satisfies the full system would be asking a
			 * question nobody has checked for an answer meq already has.
			 *
			 * `darcyRhs` is the ( flux, potential ) load and is held BY REFERENCE
			 * by the operator, so it must outlive it; it is a member, and
			 * formSystem() has already pointed it at rhs's first two blocks.
			 */
			std::unique_ptr<mfem::DarcyNPCOperator> npc;
			std::unique_ptr<mfem::DarcyNPCSolver> npcLinear;
			if ( npcOrdering )
			{
				npc = std::make_unique<mfem::DarcyNPCOperator>(
					*darcy->GetHybridization(), blockOffsets, darcyRhs );
				npcLinear = std::make_unique<mfem::DarcyNPCSolver>( linear );
			}

			mfem::Operator &residualOperator =
				npcOrdering ? static_cast<mfem::Operator &>( *npc )
				            : static_cast<mfem::Operator &>( *reduced.Ptr() );

			// The unknown, and the right hand side Newton subtracts from the
			// residual. BOTH RIGHT HAND SIDES ARE ZERO AND THEY ARE ZERO FOR
			// DIFFERENT REASONS, which is worth not conflating: under the
			// condensation ReduceRHS() zeroes traceB for a non-linear problem and
			// puts the load inside the operator, while under NPC the load is the
			// ( flux, potential ) pair passed to the operator and the trace row
			// carries none -- meq imposes psi = g_D as an ESSENTIAL condition,
			// not as a Neumann datum, and a Neumann datum is the one thing that
			// would have to ride here instead.
			mfem::Vector &unknown = npcOrdering ? solution : traceX;
			mfem::Vector newtonRhs( residualOperator.Height() );
			if ( npcOrdering )
				newtonRhs = 0.0;
			else
				newtonRhs = traceB;

			// KINSOL ignores the right hand side handed to Mult(), so the shift
			// has to be in the operator. Built unconditionally and used only on
			// the KINSOL paths: it costs one subtraction per residual evaluation
			// and keeps the two paths reading the same residual, which is what
			// makes a difference between them attributable to the line search.
			ShiftedResidual shifted( residualOperator, newtonRhs );

			std::unique_ptr<mfem::NewtonSolver> nonlinear;
			bool kinsol = false;

			switch ( globalisationChoice )
			{
				case Globalisation::None:
					// NO LINE SEARCH HERE, AND THAT IS MEASURED RATHER THAN
					// ASSUMED. Backtracking on the full residual is the
					// globalisation upstream recommends for NPC; implemented as a
					// NewtonSolver::ComputeScalingFactor subclass and measured, it
					// made EVERY case worse, including the five that converge
					// undamped. The reason is structural and is written up in
					// CLAUDE.md under *Why a line search on the full residual does
					// not work here*: the flux and trace rows of the NPC residual
					// are LINEAR, so a full step annihilates them exactly, and any
					// damping restores ( 1 - alpha ) of them. An l2 merit function
					// over the whole residual therefore rewards the very step that
					// ruins the potential block, alpha collapses to about 1e-2 and
					// the iteration creeps. KINSOL's line search fails identically,
					// on the same merit function. What works is
					// Globalisation::PicardThenNewton, which fixes WHERE the
					// iterate is rather than how far it steps.
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
				nonlinear->SetOperator( residualOperator );

			// The Jacobian solve. Under NPC it is the hybridized elimination
			// rather than a trace solve, and the difference is not optional: the
			// handle DarcyNPCOperator::GetGradient() returns is SOLVE-ONLY -- the
			// local blocks are factored in place, so J cannot be applied out of
			// them and its Mult() aborts. Handing the trace solver straight to
			// Newton would fail loudly, which is upstream's design and better
			// than the alternative.
			if ( npcOrdering )
				nonlinear->SetSolver( *npcLinear );
			else
				nonlinear->SetSolver( linear );
			nonlinear->SetMonitor( recorder );

#ifdef MFEM_USE_SUNDIALS
			// KINSOL calls LinSysSetup every tenth step by default, which under
			// NPC makes this a LAGGED-JACOBIAN Newton: legitimate, and
			// self-consistent because the reduction and the recovery both
			// eliminate with whatever factorisation is currently held, but it
			// costs iterations -- upstream measures 12 against 4 on one case,
			// both converged to round-off. meq asks for a Jacobian per step so
			// that a KINSOL run and a plain Newton run differ in the LINE SEARCH
			// and in nothing else, which is what makes
			// kinsolAgreesWithNewtonWhereBothConverge worth asserting.
			if ( kinsol && npcOrdering )
				static_cast<mfem::KINSolver &>( *nonlinear ).SetMaxSetupCalls( 1 );
#endif

			/*
			 * THE CONVERGENCE TARGET MUST NOT DEPEND ON WHERE THE ITERATION
			 * STARTED, AND MFEM'S DOES.
			 *
			 * NewtonSolver stops at || r || <= max( rel_tol * || r_0 ||, abs_tol )
			 * with || r_0 || measured at the iterate it was handed. A warm start
			 * makes || r_0 || small, so the target shrinks with it -- and a good
			 * enough guess drives the target below the round-off floor, where it
			 * can never be met. THE BETTER THE GUESS, THE MORE CERTAIN THE
			 * FAILURE, which is the exact opposite of what a restart is for.
			 *
			 * MEASURED on Example 5 at k = 3, n = 8, restarting from the converged
			 * answer. The solve reaches the floor in TWO iterations and is then
			 * reported as a failure at the thirtieth:
			 *
			 *     it 0   1.180694e-03
			 *     it 1   8.782203e-13
			 *     it 2   4.151551e-14      <- converged; the floor
			 *     ...    ~4e-14 for 28 more iterations, then FAIL
			 *
			 * because the target was max( 1e-12 * 1.18e-3, 1e-14 ) = 1e-14, under
			 * the 3.7e-14 this problem can actually reach. The same solve started
			 * cold has || r_0 || = 11.2, a target of 1.1e-11, and converges.
			 *
			 * THE FIX IS A REFERENCE THAT THE GUESS CANNOT MOVE: the residual at
			 * the COLD iterate -- the Dirichlet datum alone, which is where this
			 * solve would have started with no guess. rel_tol then keeps exactly
			 * the meaning it has always had, and a cold solve is bit-identical,
			 * because there the reference IS || r_0 ||. Only a warm one changes,
			 * and it changes from failing to converging in one step.
			 *
			 * It costs one extra residual evaluation, and only when a guess was
			 * set. On this path that is a full set of element-local solves, which
			 * is the price of the criterion meaning something.
			 */
			if ( coldReference >= 0.0 )
			{
				// A pure absolute target at the right scale. SetRelTol( 0 ) rather
				// than leaving it, because MFEM takes the LARGER of the two and a
				// live rel_tol would put || r_0 || back into the test.
				nonlinear->SetRelTol( 0.0 );
				nonlinear->SetAbsTol( std::max( newtonAbsoluteTolerance,
				                                newtonRelativeTolerance*coldReference ) );
			}
			else
			{
				nonlinear->SetRelTol( newtonRelativeTolerance );
				nonlinear->SetAbsTol( newtonAbsoluteTolerance );
			}
			nonlinear->SetMaxIter( newtonMaxIterations );
			nonlinear->SetPrintLevel( -1 );

			// The Dirichlet data rides in traceX, so the iteration must start from
			// it rather than from zero. With iterative_mode false NewtonSolver
			// zeroes x on entry and the boundary condition disappears without a
			// word -- the residual is masked on those rows, so nothing complains.
			// KINSolver::Mult has the same line and the same consequence.
			nonlinear->iterative_mode = true;

			// newtonRhs on the Newton path, where Mult() subtracts it; ignored on
			// the KINSOL paths, where ShiftedResidual has already done so. The
			// unknown is `solution` entire under NPC and the trace alone under
			// the condensation.
			nonlinear->Mult( newtonRhs, unknown );
			newtonIterationCount = nonlinear->GetNumIterations();

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
			// Recorded so the reuse can be asserted on rather than timed: a Newton
			// solve refactorises once per iteration and must analyse ONCE.
			readFactorisationCounts( linear, symbolicFactorisationCount,
			                         numericFactorisationCount );
#endif

			// Loudly, and without a recovered solution: an iteration that ran out
			// of steps has produced a vector, not an equilibrium.
			if ( !nonlinear->GetConverged() )
				throw std::runtime_error( "meq::GradShafranovSolver::solve: the non-linear iteration did not converge" );
		}
		else
		{
#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
			// No reuse here, deliberately: the linear path factorises once and
			// this object is destroyed straight after, so retaining the analysis
			// would buy nothing and cost a copy of the pattern.
			std::unique_ptr<mfem::Solver> const solver =
				makeTraceSolver( traceSolverChoice, false );
			solver->SetOperator( *reduced.Ptr() );
			solver->Mult( traceB, traceX );
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

		if ( !fieldsAreState )
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

	long GradShafranovSolver::localNonlinearIterations() const
	{
		if ( !built )
			throw std::logic_error( "meq::GradShafranovSolver::localNonlinearIterations: the forms have not been built" );
		return darcy->GetHybridization()->GetNumLocalNLIterations();
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
