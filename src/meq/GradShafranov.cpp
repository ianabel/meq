#include "GradShafranov.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
	 * THE ASSEMBLY MODE A FRESH SOLVER STARTS IN, AND IT CHANGED ON 2026-09-04
	 * FROM Serial TO Threaded. THE MEASUREMENT THAT SETTLED IT BEFORE WAS TAKEN
	 * AGAINST A DIFFERENT OPTION.
	 *
	 * What this used to say, and it was right at the time: a gate on
	 * omp_get_max_threads() was written and REMOVED because
	 * `HighBetaConvergence` went 21.5 s to 39 s under it -- 1.8x SLOWER,
	 * reproducibly -- since MFEM forks a team and buffers element blocks PER
	 * CALL, so a caller that assembles hundreds of times inside a bordered
	 * Newton pays that every time while one that assembles once amortises it.
	 * Mesh size did not separate the two cases -- HighBeta's meshes are 128 and
	 * 512 elements and 512 is where the isolated benchmark still showed a win --
	 * so the solver could not know which caller it had, and the honest default
	 * was MFEM's own.
	 *
	 * THAT ARGUMENT WAS ABOUT ComputeH() AND ONLY ComputeH(), because in August
	 * that was the only loop AssemblyMode touched. MFEM has since threaded
	 * MultNL() as well -- the residual and the Jacobian assembly, and so
	 * NPCResidual() and NPCGradient(), which is every NPC step. A bordered
	 * Newton is dominated by residual evaluations rather than by assembly, so
	 * the case that most OPPOSED the flag is now the case that most favours it.
	 *
	 * RE-MEASURED, SAME TEST, SAME MACHINE, MKL_NUM_THREADS=1:
	 *
	 *     HighBetaConvergence    Serial 3.22, 3.19 s     Threaded 1.34, 1.33 s
	 *
	 * 2.4x FASTER where it was 1.8x slower. That inversion is the whole
	 * justification, and it makes this a change of measurement rather than a
	 * change of mind.
	 *
	 * On a whole nonlinear solve the flag is worth 2.8x-3.0x at eight threads:
	 * example5 at k = 2, n = 24 goes 0.369 s to 0.131 s, and at k = 3, n = 16
	 * 0.280 s to 0.093 s.
	 *
	 * AND AT ONE THREAD IT IS A WASH, WHICH IS WHAT MAKES IT SAFE AS A DEFAULT.
	 * The old note recorded 0.86x at one thread, the chunk buffering costing
	 * more than the serial loop it imitates. Measured now: 0.3664 s threaded
	 * against 0.3686 s serial, a ratio of 1.01. A build with OpenMP that happens
	 * to run at OMP_NUM_THREADS=1 loses nothing.
	 *
	 * THE LARGER REASON IS NOT A SPEEDUP AT ALL. Threaded assembly is what makes
	 * MKL_NUM_THREADS > 1 survivable. MKL suppresses its own threading inside an
	 * ACTIVE OpenMP region, so the element-local dense work is nested and pays
	 * nothing for MKL threads; the same work in a SERIAL element loop pays for
	 * them on every call. Measured on a whole nonlinear solve, k = 3, n = 16,
	 * OMP_NUM_THREADS=8:
	 *
	 *                       MKL=1         MKL=8
	 *      Serial          0.2797 s     107.19 s     <- 383x
	 *      Threaded        0.0927 s       0.0834 s   <- immune
	 *
	 * 1285x between the two modes at MKL=8. That is the answer to CLAUDE.md's
	 * *What to do* item 0 -- "get ComputeH()'s element-local dense LU off
	 * threaded MKL" -- reached without writing any of the code that item
	 * proposes, and it is what makes PARDISO's MKL threads spendable: the trace
	 * solve runs on the master thread OUTSIDE any parallel region and takes all
	 * of them.
	 *
	 * WHAT IS GIVEN UP, AND IT IS REAL. The two modes agree BIT FOR BIT only at
	 * MKL_NUM_THREADS=1. Above that they differ at round-off -- measured 1.3e-15
	 * in psi and 1.3e-13 in the flux -- because the serial path hands MKL eight
	 * threads and the threaded path hands it one, and a blocked BLAS-3 sums in a
	 * different order from an unblocked loop. That is arithmetic reassociation
	 * inside MKL, not a race in MFEM or in MEQ, and the suite pins exactness at
	 * MKL_NUM_THREADS=1 which every registered test sets.
	 */
	GradShafranovSolver::AssemblyMode defaultAssemblyMode()
	{
		/*
		 * BUILD-CONDITIONAL, AND IT HAS TO BE. buildForms() hands this straight
		 * to DarcyHybridization::SetAssemblyMode(), and the comment there says
		 * setAssemblyMode() has already refused Threaded if the build cannot
		 * honour it -- a guarantee that holds for the SETTER only. The
		 * constructor assigns this value directly, so a Threaded default on a
		 * build without OpenMP or without thread safety would reach MFEM
		 * unchecked and ABORT THE PROCESS rather than throw.
		 *
		 * So the capability question is asked here too. It is the same question
		 * assemblyModeAvailable() answers, asked at the one place that bypasses
		 * the setter.
		 */
#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
		return GradShafranovSolver::AssemblyMode::Threaded;
#else
		return GradShafranovSolver::AssemblyMode::Serial;
#endif
	}

#ifdef MEQ_HAVE_DIRECT_TRACE_SOLVER
	/*
	 * Build the chosen direct solver, configured as MEQ wants it.
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
				// both. So this is the type that is right for MEQ's headline
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
			"MEQ: the requested trace solver is not available in this build, and "
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
#ifdef MFEM_THREAD_SAFE
		// Local, because this runs on DarcyHybridization::MultNL()'s threaded
		// element loop. See the declaration in the header.
		mfem::Vector shape( dof );
#else
		shape.SetSize( dof );
#endif
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
#ifdef MFEM_THREAD_SAFE
		mfem::Vector shape( dof );
#else
		shape.SetSize( dof );
#endif
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
		  // In DECLARATION order, which is the order these are actually
		  // constructed in whatever this list says. -Wreorder had transferPath
		  // and extensionLineOrder written after picardDamping and initialised
		  // before it, which costs nothing while every entry is a scalar or a
		  // null pointer and is a trap the moment one of them reads another.
		  transferPath( nullptr ),
		  extensionLineOrder( -1 ),
		  /*
		   * 40, NOT 12, AND THE DIFFERENCE WAS MEASURED RATHER THAN CHOSEN.
		   *
		   * This is a rule ACROSS a face of Gamma_h, integrating a foot map
		   * xi -> a( x( xi ) ) whose smoothness is the path family's business
		   * and not MEQ's. mfem::VertexConePath's cone drives the two
		   * interpolated vertex directions apart and roughens that map, and at
		   * order 12 a quadrature of it is short by O( h^2 ):
		   * FreeBoundaryCoupling's tiling sweep reads 1.01e-04 at order 12 and
		   * 6.68e-08 at 40 on the same geometry, against a coverage floor near
		   * 1e-10.
		   *
		   * A caller cannot know whether the path handed in cones -- HasCone()
		   * is on the concrete class, not on mfem::TransferPath -- so the
		   * default has to be adequate for one that does. 40 is; 12 is not.
		   * Raise it with setTransmissionQuadratureOrder() and check the answer
		   * stops moving, which is the only way to know it is enough.
		   *
		   * Still a setup cost paid once per mesh rather than once per Newton
		   * step, so the order is cheap.
		   */
		  transmissionQuadratureOrder( 40 ),
		  globalisationChoice( Globalisation::None ),
		  localSolverChoice( LocalSolver::Newton ),
		  exteriorCoupling( nullptr ),
		  sourceQuadratureExtra( 4 ),
		  orderingChoice( NonlinearOrdering::NPC ),
		  assemblyModeChoice( defaultAssemblyMode() ),
		  traceSolverChoice( TraceSolver::UMFPack ),
		  andersonDepth( 1 ),
		  picardDamping( 1.0 ),
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
		 * So it is kept for alignment with the miniapp MEQ was ported from, which
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
		return recoverPeak( trace, psiAxisValue, psiBoundaryValue, element );
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
	 * solve -- which is the ordering Nguyen, Peraire & Cockburn use and MEQ's
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
	 * not a library's decision to impose on MEQ's callers, so the build is
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

	bool GradShafranovSolver::assemblyModeAvailable( AssemblyMode choice )
	{
		if ( choice == AssemblyMode::Serial )
			return true;

#if defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )
		return true;
#else
		return false;
#endif
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

	void GradShafranovSolver::setSourceQuadratureOrder( int extraOrder )
	{
		if ( extraOrder < 0 )
			throw std::invalid_argument(
				"meq::GradShafranovSolver::setSourceQuadratureOrder: the extra "
				"order is added to 2k and cannot be negative" );
		sourceQuadratureExtra = extraOrder;
		built = false;
		prepared = false;
	}

	int GradShafranovSolver::sourceQuadratureOrder() const
	{
		return sourceQuadratureExtra;
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

	double GradShafranovSolver::starShapedMargin( double centreR,
	                                              double centreZ ) const
	{
		mfem::Mesh &mesh = *traceFes->GetMesh();
		bool const extended = ( transferPath != nullptr );

		double margin = 1.0;
		bool sawAny = false;

		mfem::Vector x( 2 );
		mfem::Vector normal( 2 );

		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			if ( extended )
			{
				int const attribute = mesh.GetBdrAttribute( be );
				if ( attribute < 1 || attribute > gammaHMarker.Size()
				     || !gammaHMarker[ attribute - 1 ] )
					continue;
			}

			mfem::FaceElementTransformations *ftr =
				mesh.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;

			/*
			 * Both ENDS of the face, and that is enough on a straight one: the
			 * quantity ( x - c ).n is affine in x along the face because n is
			 * constant there, so its minimum is at an endpoint. A curved face
			 * would need the interior too, and MEQ's Gamma_h is straight --
			 * it is the union of background element faces.
			 */
			for ( double xi : { 0.0, 1.0 } )
			{
				mfem::IntegrationPoint ip;
				ip.Set1w( xi, 1.0 );
				ftr->SetAllIntPoints( &ip );
				ftr->Transform( ip, x );

				mfem::CalcOrtho( ftr->Jacobian(), normal );
				double const length = normal.Norml2();
				if ( !( length > 0.0 ) )
					continue;
				normal /= length;

				double const dr = x( 0 ) - centreR;
				double const dz = x( 1 ) - centreZ;
				double const distance = std::hypot( dr, dz );
				if ( !( distance > 0.0 ) )
					continue;

				/*
				 * CalcOrtho's sign follows the face's own parametrisation and
				 * says nothing about which side the domain is on, so it is
				 * oriented here the same way ExtensionBoundaryQuadrature
				 * orients its own: outward means agreeing with the direction
				 * away from the interior. For a boundary face, Elem1 is the
				 * interior element, so pointing away from its centre is the
				 * test -- and using the CENTRE rather than the candidate star
				 * centre is deliberate, since the latter is what is under test.
				 */
				mfem::Vector elementCentre( 2 );
				mesh.GetElementCenter( ftr->Elem1No, elementCentre );
				double const outward = normal( 0 )*( x( 0 ) - elementCentre( 0 ) )
				                     + normal( 1 )*( x( 1 ) - elementCentre( 1 ) );
				double const sign = ( outward < 0.0 ) ? -1.0 : 1.0;

				double const cosine =
					sign*( normal( 0 )*dr + normal( 1 )*dz )/distance;

				margin = std::min( margin, cosine );
				sawAny = true;
			}
		}

		if ( !sawAny )
			throw std::logic_error(
				"meq::GradShafranovSolver::starShapedMargin: no boundary faces "
				"were examined; on the extension path that means the Gamma_h "
				"marker selects nothing" );

		return margin;
	}

	/*
	 * PROJECTING A PATH COEFFICIENT ONTO Gamma_h'S TRACE DOFS, AND WHY
	 * GridFunction::ProjectBdrCoefficient CANNOT DO IT.
	 *
	 * FREE-BOUNDARY-PLAN.md section 4.3 says each column of P is "one call to
	 * ProjectBdrCoefficient against mfem::PathTraceCoefficient". IT IS NOT, and
	 * the attempt aborts rather than misbehaving quietly:
	 *
	 *     PathTraceCoefficient must be evaluated on a face: the path family may
	 *     need the outward normal of Gamma_h
	 *
	 * ProjectBdrCoefficient evaluates through the ELEMENT transformation, and a
	 * path coefficient needs the FACE one -- mfem::TransferredDatumCoefficient's
	 * own header says the same thing about itself, and this is the same
	 * requirement arriving one class over. So the projection is written out
	 * here.
	 *
	 * IT IS projectOntoTrace() WITH TWO CHANGES, and both are the point. The
	 * coefficient is evaluated on `*ftr` rather than on `*ftr->Elem1`, which is
	 * what supplies the face and its normal; and the loop is over BOUNDARY
	 * ELEMENTS carrying a Gamma_h attribute rather than over every face, since
	 * a mode has nothing to say anywhere else. Everything else -- nodal
	 * interpolation at the face element's own nodes, DG_Interface being nodal
	 * with VALUE map type so the nodes are where the dofs live -- is that
	 * function's reasoning and is not repeated.
	 *
	 * Serial, and it uses Mesh's shared face transformation deliberately: this
	 * runs once per mode at setup, never inside an element loop, so the
	 * reentrancy hazard CLAUDE.md records for GetBdrFaceTransformations does not
	 * arise. If it is ever threaded, the caller-allocated overload is the fix.
	 */
	void GradShafranovSolver::projectPathTraceOntoGammaH(
		mfem::Coefficient &coeff, mfem::Vector &target ) const
	{
		mfem::Mesh &mesh = *traceFes->GetMesh();
		mfem::Array<int> vdofs;
		mfem::Vector values;

		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			int const attribute = mesh.GetBdrAttribute( be );
			if ( attribute < 1 || attribute > gammaHMarker.Size()
			     || !gammaHMarker[ attribute - 1 ] )
				continue;

			mfem::FaceElementTransformations *ftr =
				mesh.GetBdrFaceTransformations( be );
			if ( !ftr )
				continue;

			int const face = mesh.GetBdrElementFaceIndex( be );
			mfem::FiniteElement const *faceFe = traceFes->GetFaceElement( face );
			if ( !faceFe )
				continue;

			traceFes->GetFaceVDofs( face, vdofs );
			int const dof = faceFe->GetDof();
			values.SetSize( dof );

			mfem::IntegrationRule const &nodes = faceFe->GetNodes();
			for ( int i = 0; i < dof; ++i )
			{
				ftr->SetAllIntPoints( &nodes.IntPoint( i ) );
				// ON THE FACE, not on Elem1. That is the whole difference.
				values( i ) = coeff.Eval( *ftr, nodes.IntPoint( i ) );
			}

			target.SetSubVector( vdofs, values );
		}
	}

	/*
	 * THE COLUMNS OF P, AND WHY EACH ONE IS A PROJECTION AT THE FOOT OF A PATH
	 * RATHER THAN A PROJECTION ON Gamma_h.
	 *
	 * The exterior expansion lives on Gamma, the TRUE boundary. Gamma_h is the
	 * inscribed polygon the mesh actually has, and the whole of stage 5 is the
	 * machinery for carrying a datum between them: for a point x on Gamma_h,
	 * a( x ) is its foot on Gamma, and mfem::PathTraceCoefficient( path, g )
	 * evaluates g there. So column n is the projection onto Gamma_h's trace dofs
	 * of C_n evaluated ON Gamma -- not of C_n evaluated on Gamma_h, which would
	 * be a different function by O( dist( Gamma_h, Gamma ) ) and would throw the
	 * whole transfer technique's accuracy away at the one place it is needed.
	 *
	 * This is the same call the fixed-boundary datum already makes, with g = 0
	 * replaced by a mode. That is section 4.3's point: the coupling is stage 5
	 * with the datum unknown instead of zero, so the machinery is already here.
	 *
	 * ONE PROJECTION PER MODE, AND THE MARKER IS gammaHMarker. Note the contrast
	 * with prepare(), which projects the fixed-boundary datum against
	 * fittedMarker precisely because Gamma_h's dofs are pinned to zero there.
	 * Free boundary un-pins exactly those, and leaves the fitted ones alone.
	 */
	std::vector<mfem::Vector>
		GradShafranovSolver::exteriorTraceColumns( ExteriorDtN const &exterior ) const
	{
		if ( !transferPath )
			throw std::logic_error(
				"meq::GradShafranovSolver::exteriorTraceColumns: there is no "
				"Gamma_h on the fitted path -- setExtension() was never called, "
				"so the trace unknown on the boundary IS the condition imposed "
				"and there is nothing for an exterior expansion to drive" );

		std::vector<mfem::Vector> columns;
		columns.reserve( static_cast<std::size_t>( exterior.modeCount() ) );

		for ( int i = 0; i < exterior.modeCount(); ++i )
		{
			int const n = ExteriorDtN::firstMode() + i;

			/*
			 * The mode as a PositionFunction. Captured by value except for the
			 * ExteriorDtN, which the caller owns and which must outlive the
			 * projection -- it does, since the projection happens inside this
			 * loop and nothing escapes.
			 *
			 * basis() takes ( r, z ) and depends on the DIRECTION alone, so it
			 * is well defined at any point of the plane and in particular at a
			 * foot on Gamma, which is where PathTraceCoefficient evaluates it.
			 */
			mfem::PositionFunction mode =
				[ &exterior, n ]( mfem::Vector const &x )
			{
				return exterior.basis( n, x( 0 ), x( 1 ) );
			};

			mfem::PathTraceCoefficient traceOfMode( *transferPath, mode );

			mfem::Vector column( traceFes->GetVSize() );
			column = 0.0;
			projectPathTraceOntoGammaH( traceOfMode, column );
			columns.push_back( std::move( column ) );
		}

		return columns;
	}

	void GradShafranovSolver::setExteriorDatum( mfem::PositionFunction g )
	{
		if ( !transferPath )
			throw std::logic_error(
				"meq::GradShafranovSolver::setExteriorDatum: there is no Gamma_h "
				"on the fitted path -- setBoundaryData() is the datum there, and "
				"it reaches the boundary through essential trace dofs rather than "
				"through a load term" );

		exteriorDatumFunction = std::move( g );

		// prepare() is where it reaches the right hand side, and solve()
		// re-prepares, so a datum set after a solve would otherwise be ignored
		// until something else invalidated the state.
		prepared = false;
	}

	void GradShafranovSolver::setExteriorCoupling( ExteriorDtN const &exterior )
	{
		if ( !transferPath )
			throw std::logic_error(
				"meq::GradShafranovSolver::setExteriorCoupling: there is no "
				"Gamma_h to transfer the exterior datum from -- setExtension() "
				"comes first, and the fitted path cannot carry this coupling at "
				"all because its datum reaches the boundary through essential "
				"trace dofs rather than through a load term" );

		exteriorCoupling = &exterior;
		exteriorCoefficientValues.assign(
			static_cast<std::size_t>( exterior.modeCount() ), 0.0 );

		// THE DATUM READS THE COEFFICIENT VECTOR RATHER THAN A COPY OF IT, so
		// that a Newton step which moves `a` moves the transferred boundary
		// condition with it. The vector is a member and outlives every solve;
		// capturing it by reference is what makes the coupling live rather than
		// a snapshot taken at setup.
		ExteriorDtN const *dtn = &exterior;
		std::vector<double> const *coefficients = &exteriorCoefficientValues;
		setExteriorDatum( [ dtn, coefficients ]( mfem::Vector const &x )
		{
			double total = 0.0;
			int const first = ExteriorDtN::firstMode();
			for ( std::size_t i = 0; i < coefficients->size(); ++i )
				total += ( *coefficients )[ i ]
				         *dtn->basis( first + static_cast<int>( i ), x( 0 ), x( 1 ) );
			return total;
		} );
	}

	std::vector<double> const &GradShafranovSolver::exteriorCoefficients() const
	{
		return exteriorCoefficientValues;
	}

	double GradShafranovSolver::outwardFlux() const
	{
		mfem::Mesh &mesh = *traceFes->GetMesh();
		mfem::Array<int> vdofs;
		double total = 0.0;

		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			int const attribute = mesh.GetBdrAttribute( be );
			if ( attribute < 1 || attribute > gammaHMarker.Size() )
				continue;
			bool const transferred = transferPath
			                         && gammaHMarker[ attribute - 1 ];

			thread_local mfem::FaceElementTransformations faceScratch;
			thread_local mfem::IsoparametricTransformation faceElem1;
			thread_local mfem::IsoparametricTransformation faceElem2;
			mesh.GetBdrFaceTransformations( be, faceScratch, faceElem1, faceElem2 );
			if ( faceScratch.GetGeometryType() == mfem::Geometry::INVALID )
				continue;

			int const element = faceScratch.Elem1No;
			mfem::FiniteElement const *fluxFe = fluxFes->GetFE( element );
			if ( !fluxFe )
				continue;

			fluxFes->GetElementVDofs( element, vdofs );
			int const dof = fluxFe->GetDof();
			int const dim = mesh.Dimension();

			thread_local mfem::IsoparametricTransformation elementScratch;
			mesh.GetElementTransformation( element, &elementScratch );
			mfem::ElementExtension extender;
			extender.SetElement( elementScratch );

			mfem::IntegrationRule const &faceRule =
				mfem::IntRules.Get( faceScratch.GetGeometryType(),
				                    transmissionQuadratureOrder );

			mfem::Vector shape( dof );
			bool reached = true;

			// q.nu from the element's flux dofs at a reference point. The minus
			// undoes DarcyForm's convention, exactly as the transmission rows
			// do: the block holds -q and the identity is written for q.
			auto normalFlux = [ & ]( mfem::IntegrationPoint const &eip,
			                         mfem::Vector const &nu )
			{
				fluxFe->CalcShape( eip, shape );
				double normalComponent = 0.0;
				for ( int d = 0; d < dim; ++d )
				{
					double component = 0.0;
					for ( int j = 0; j < dof; ++j )
						component += shape( j )*solution( vdofs[ dof*d + j ] );
					normalComponent += component*nu( d );
				}
				return -normalComponent;
			};

			if ( transferred )
			{
				mfem::ExtensionBoundaryQuadrature( faceScratch, *transferPath,
					faceRule,
					[ & ]( mfem::ExtensionBoundaryPoint const &pt )
				{
					if ( !reached )
						return;

					mfem::IntegrationPoint eip;
					if ( !extender.TransformBack( pt.y, eip ) )
					{
						reached = false;
						return;
					}
					// pt.weight is SIGNED and is used as it stands.
					total += pt.weight*normalFlux( eip, pt.nu );
				} );
			}
			else
			{
				/*
				 * A FITTED FACE, AND ON THE HALF-DISC THOSE ARE THE AXIS.
				 *
				 * LEAVING THEM OUT IS WHAT MADE THIS 93% WRONG AND FLAT.
				 * Gamma is a SEMICIRCLE, so the boundary enclosing the current
				 * is the arc PLUS the axis segment, and the divergence theorem
				 * wants all of it. The axis does not contribute zero: psi ~
				 * c( z ) r^2 there, so q_r = 2c is finite and generally
				 * non-zero even though psi itself vanishes.
				 *
				 * The symptom was diagnostic once seen -- an error that does
				 * not move under either h or k is not a discretisation error,
				 * and this one sat at 93% across three meshes and three
				 * degrees. It was a missing PIECE OF THE CONTOUR.
				 *
				 * No extension here: a fitted face is on Gamma already.
				 */
				mfem::Vector nu( dim );
				for ( int i = 0; i < faceRule.GetNPoints(); ++i )
				{
					mfem::IntegrationPoint const &ip = faceRule.IntPoint( i );
					faceScratch.SetAllIntPoints( &ip );
					mfem::CalcOrtho( faceScratch.Jacobian(), nu );
					double const measure = nu.Norml2();
					if ( !( measure > 0.0 ) )
						continue;
					nu /= measure;

					// CalcOrtho follows the face's parametrisation, so orient
					// it away from the interior element, as
					// ExtensionBoundaryQuadrature orients its own.
					mfem::Vector centre( dim ), here( dim );
					mesh.GetElementCenter( faceScratch.Elem1No, centre );
					faceScratch.Transform( ip, here );
					double outward = 0.0;
					for ( int d = 0; d < dim; ++d )
						outward += nu( d )*( here( d ) - centre( d ) );
					if ( outward < 0.0 )
						nu.Neg();

					total += ip.weight*measure
					         *normalFlux( faceScratch.GetElement1IntPoint(), nu );
				}
			}

			if ( !reached )
				throw std::runtime_error(
					"meq::GradShafranovSolver::outwardFlux: the extension of an "
					"element of Gamma_h did not reach its foot on Gamma" );
		}

		return total;
	}

	void GradShafranovSolver::setBoundaryFluxPoint( double r, double z )
	{
		if ( !std::isfinite( r ) || !std::isfinite( z ) )
			throw std::invalid_argument(
				"meq::GradShafranovSolver::setBoundaryFluxPoint: the limiter "
				"contact must be finite" );
		if ( orderingChoice != NonlinearOrdering::NPC )
			throw std::logic_error(
				"meq::GradShafranovSolver::setBoundaryFluxPoint: psi_bnd as an "
				"unknown is implemented for NonlinearOrdering::NPC only -- under "
				"the condensation psi is a function of the trace through every "
				"element's source, so both the border row and its corner would "
				"have to be differenced rather than being exact" );

		boundaryFluxIsUnknown = true;
		boundaryFluxR = r;
		boundaryFluxZ = z;
		prepared = false;
	}

	double GradShafranovSolver::psiBoundary() const
	{
		return psiBoundaryValue;
	}

	/*
	 * THE POTENTIAL DOF NEAREST A POINT.
	 *
	 * psi_bnd is pinned to ONE nodal value, exactly as psi_ax is pinned to the
	 * largest one, because that is what makes the constraint differentiable in a
	 * form the border can use: the row is then a unit vector and nothing about it
	 * is measured. Choosing the nearest dof rather than interpolating is the same
	 * trade CLAUDE.md records for psi_ax -- it differs from psi at the point by
	 * O( h^{k+1} ) and both converge to it.
	 *
	 * Walked once at setup over every element's nodes, which is why it is not
	 * worth an octree.
	 */
	int GradShafranovSolver::nearestPotentialDof( double r, double z ) const
	{
		mfem::Mesh &mesh = *potentialFes->GetMesh();
		mfem::Array<int> dofs;
		mfem::Vector point( 2 );

		double best = std::numeric_limits<double>::infinity();
		int bestDof = -1;

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::FiniteElement const *fe = potentialFes->GetFE( e );
			if ( !fe )
				continue;
			potentialFes->GetElementDofs( e, dofs );

			thread_local mfem::IsoparametricTransformation scratch;
			mesh.GetElementTransformation( e, &scratch );

			mfem::IntegrationRule const &nodes = fe->GetNodes();
			for ( int i = 0; i < fe->GetDof(); ++i )
			{
				scratch.Transform( nodes.IntPoint( i ), point );
				double const d = std::hypot( point( 0 ) - r, point( 1 ) - z );
				if ( d < best )
				{
					best = d;
					bestDof = dofs[ i ];
				}
			}
		}

		if ( bestDof < 0 )
			throw std::runtime_error(
				"meq::GradShafranovSolver::nearestPotentialDof: the mesh has no "
				"potential dofs" );
		return bestDof;
	}

	void GradShafranovSolver::setTransmissionQuadratureOrder( int order )
	{
		if ( order < 0 )
			throw std::invalid_argument(
				"meq::GradShafranovSolver::setTransmissionQuadratureOrder: the "
				"rule order must not be negative" );

		transmissionQuadratureOrder = order;
	}

	/*
	 * THE ROWS OF T: THE TRANSMISSION CONDITION, TESTED AGAINST EACH MODE.
	 *
	 * This is the Neumann half of the coupling and the piece FREE-BOUNDARY-PLAN.md
	 * section 7.5 called "genuinely new code". Its MFEM half is not: the sweep of
	 * Gamma is mfem::ExtensionBoundaryQuadrature(), which was written for this and
	 * merged upstream. What is here is the contraction.
	 *
	 * WHAT THE QUADRATURE HANDS BACK, AND THE TWO THINGS IN IT THAT ARE EASY TO
	 * GET WRONG.
	 *
	 *   pt.y       the point ON GAMMA -- the foot a( x ), not the face point x.
	 *   pt.nu      the OUTWARD unit normal of Gamma, oriented by the PATHS rather
	 *              than by the face. Those differ wherever Gamma and Gamma_h are
	 *              not parallel, which is everywhere that matters.
	 *   pt.weight  SIGNED, and the sign is the point. A staircase Gamma_h has
	 *              faces whose foot map reverses along Gamma; summing the weights
	 *              integrates over Gamma, while summing their absolute values
	 *              integrates the length the map traverses, which is larger and is
	 *              a different quantity. Writing std::abs here would silently
	 *              reintroduce exactly the defect the signed weight exists to fix.
	 *
	 * A degenerate face -- one whose image on Gamma is a single point -- is skipped
	 * by the routine rather than refused, so visit() is simply not called there.
	 * That is correct: it covers none of Gamma and its neighbours cover Gamma
	 * between them.
	 *
	 * THE MEASURE IS PLAIN dGamma AND THE 1/r IS ALREADY IN q. The exterior block
	 * is diagonal in the weight dGamma/r -- that is section 3.2, and it is what
	 * makes the whole method cheap -- so the interior term must be tested in the
	 * same weight. It is, without dividing by anything: the condition matches
	 * ( 1/r ) dpsi/dnu across Gamma, and MEQ's q IS ( 1/r ) grad_bar( psi ), so
	 * q.nu tested in the plain measure already carries the radius the exterior side
	 * carries in its weight. Dividing by pt.y( 0 ) here would do it twice. The
	 * header says this at more length; it is repeated because the wrong version
	 * converges.
	 *
	 * AND IT IS LINEAR IN q, WHICH IS WHY A ROW EXISTS AT ALL. E_h is a linear
	 * operator on the flux dofs and nu, C_m and the measure are geometry, so the
	 * integral is a fixed covector applied to the flux block. It is built by
	 * pushing each flux basis function through the extension in turn, which is
	 * what the inner loop over element dofs is.
	 */
	std::vector<mfem::Vector>
		GradShafranovSolver::exteriorTransmissionRows( ExteriorDtN const &exterior ) const
	{
		if ( !transferPath )
			throw std::logic_error(
				"meq::GradShafranovSolver::exteriorTransmissionRows: there is no "
				"Gamma_h on the fitted path -- there is no band between Gamma_h "
				"and Gamma for the extension to cross, and so no transmission "
				"condition to impose" );

		int const modes = exterior.modeCount();
		int const fluxSize = fluxFes->GetVSize();

		std::vector<mfem::Vector> rows;
		rows.reserve( static_cast<std::size_t>( modes ) );
		for ( int i = 0; i < modes; ++i )
		{
			rows.emplace_back( solution.Size() );
			rows.back() = 0.0;
		}

		mfem::Mesh &mesh = *traceFes->GetMesh();
		mfem::Array<int> vdofs;

		for ( int be = 0; be < mesh.GetNBE(); ++be )
		{
			int const attribute = mesh.GetBdrAttribute( be );
			if ( attribute < 1 || attribute > gammaHMarker.Size()
			     || !gammaHMarker[ attribute - 1 ] )
				continue;

			/*
			 * The caller-allocated variants throughout, per CLAUDE.md: both
			 * GetBdrFaceTransformations( int ) and GetElementTransformation( int )
			 * hand out the Mesh's own shared scratch, and this loop is exactly the
			 * shape that would be threaded one day. The face transformation has to
			 * outlive the sweep, since ExtensionBoundaryQuadrature reads it at
			 * every quadrature point.
			 */
			thread_local mfem::FaceElementTransformations faceScratch;
			thread_local mfem::IsoparametricTransformation faceElem1;
			thread_local mfem::IsoparametricTransformation faceElem2;
			mesh.GetBdrFaceTransformations( be, faceScratch, faceElem1, faceElem2 );
			if ( faceScratch.GetGeometryType() == mfem::Geometry::INVALID )
				continue;

			int const element = faceScratch.Elem1No;
			mfem::FiniteElement const *fluxFe = fluxFes->GetFE( element );
			if ( !fluxFe )
				continue;

			fluxFes->GetElementVDofs( element, vdofs );
			int const dof = fluxFe->GetDof();
			int const dim = mesh.Dimension();

			/*
			 * A SECOND element transformation, and it must not be the mesh's:
			 * ElementExtension::TransformBack runs a Newton solve that moves the
			 * transformation's own integration point, and faceScratch.Elem1 is
			 * being used by the sweep at the same time. Two separate objects is
			 * the whole fix.
			 */
			thread_local mfem::IsoparametricTransformation elementScratch;
			mesh.GetElementTransformation( element, &elementScratch );
			mfem::ElementExtension extender;
			extender.SetElement( elementScratch );

			mfem::IntegrationRule const &faceRule =
				mfem::IntRules.Get( faceScratch.GetGeometryType(),
				                    transmissionQuadratureOrder );

			/*
			 * A SCALAR shape vector, not a DenseMatrix, because the flux space is
			 * an L2_FECollection with vdim 2 rather than a vector FE: GetFE()
			 * returns the scalar element and CalcVShape() would refuse it. The
			 * ordering is byNODES, so component d of basis function j is vdof
			 * dof*d + j -- and getting THAT wrong swaps the radial and vertical
			 * components of the normal trace, which is a rotation of the field
			 * and converges to a plausible wrong answer.
			 */
			mfem::Vector shape( dof );
			bool reached = true;

			mfem::ExtensionBoundaryQuadrature( faceScratch, *transferPath, faceRule,
				[ & ]( mfem::ExtensionBoundaryPoint const &pt )
			{
				if ( !reached )
					return;

				// E_h evaluated at the foot: the element's own polynomial, read
				// outside it. mfem::ElementExtension is what does not clamp the
				// reference point back into the element -- an ordinary
				// TransformBack does, and a clamped point turns the extension
				// into a constant without saying so.
				mfem::IntegrationPoint eip;
				if ( !extender.TransformBack( pt.y, eip ) )
				{
					reached = false;
					return;
				}

				fluxFe->CalcShape( eip, shape );

				for ( int m = 0; m < modes; ++m )
				{
					int const n = ExteriorDtN::firstMode() + m;
					double const mode =
						exterior.basis( n, pt.y( 0 ), pt.y( 1 ) );

					// pt.weight is SIGNED and is used as it stands. The minus
					// undoes DarcyForm's convention: the flux block holds -q, so
					// the row that multiplies it must carry the sign that turns
					// it back into q. See the file comment.
					double const factor = -pt.weight*mode;
					mfem::Vector &row = rows[ static_cast<std::size_t>( m ) ];

					for ( int d = 0; d < dim; ++d )
					{
						double const weighted = factor*pt.nu( d );
						for ( int j = 0; j < dof; ++j )
							row( vdofs[ dof*d + j ] ) += weighted*shape( j );
					}
				}
			} );

			if ( !reached )
				throw std::runtime_error(
					"meq::GradShafranovSolver::exteriorTransmissionRows: the "
					"extension of an element of Gamma_h did not reach its foot on "
					"Gamma -- the inverse element map failed to converge. That is "
					"assumption P.1 giving way: dist( Gamma_h, Gamma ) has grown "
					"large against the local mesh size. meq::AdaptiveDomain is "
					"what keeps it bounded through refinement" );
		}

		// Sanity: a row must live on the flux block alone. Anything outside it is
		// an indexing error, and a silent one -- the bordered solve would simply
		// couple the exterior to a potential or trace dof and converge to
		// something.
		for ( auto const &row : rows )
			for ( int i = fluxSize; i < row.Size(); ++i )
				if ( row( i ) != 0.0 )
					throw std::logic_error(
						"meq::GradShafranovSolver::exteriorTransmissionRows: a "
						"transmission row has an entry outside the flux block" );

		return rows;
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
			// same thing from MEQ's side, but it was checked: with 1/r here the
			// error is flat under refinement, exactly as it is when the flux mass
			// form itself is given 1/r.
			//
			// The sign is +1, HDGExtensionIntegrator's own default, and it is the
			// default for the same reason it is right here: DarcyForm's flux block
			// holds -q, which is precisely the u = -K grad p of the Darcy problem
			// the extension was written for, so MEQ's convention and the
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
			potentialMass->AddDomainIntegrator(
				new SourceIntegrator( *nonlinearSource, sourceQuadratureExtra ) );
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
			// NonlinearOrdering::NPC, which is MEQ's default. They are set
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

		/*
		 * THE DATUM ON Gamma_h, AND IT IS A LOAD TERM RATHER THAN AN ESSENTIAL
		 * VALUE. THE FIRST ATTEMPT SET THE TRACE AND DID NOTHING AT ALL.
		 *
		 * Setting Gamma_h's trace dofs and letting FormLinearSystem eliminate
		 * them is how the FITTED datum is imposed, four lines above, and it is
		 * inert on Gamma_h. The reason is on the flux divergence form: its
		 * boundary face integrator carries `fittedMarker`, and
		 * EnableHybridization registers a boundary flux constraint on exactly
		 * the attributes it finds marked there. Gamma_h is NOT among them, so
		 * its trace dofs are essential in name while nothing couples to them --
		 * eliminating them at any value whatever changes not one digit, which is
		 * what was measured before this was understood.
		 *
		 * MFEM'S OWN miniapps/hdg/extension.cpp IS THE WORKED EXAMPLE and it does
		 * it the other way, which is the way the weak form asks for:
		 *
		 *     fform->AddBdrFaceIntegrator(
		 *        new VectorBoundaryFluxLFIntegrator( datum ), bdr_gamma_h );
		 *
		 * with `datum` a PathTraceCoefficient. That is < psihat, v.n > of (8a) as
		 * a LOAD on the flux equation, and it is where a non-homogeneous g
		 * belongs. HDGExtensionIntegrator supplies the other half, the path
		 * integral of the flux, which is the whole of it only when g vanishes.
		 *
		 * THE SIGN IS THE MINIAPP'S, AND IT NEGATES: pNatural = -pExact, "the
		 * datum as the flux equation takes it". MEQ's flux block holds -q for the
		 * same reason the Darcy problem's does, so the two conventions coincide
		 * and the negation carries over. It is asserted by a rate rather than by
		 * this paragraph -- see FreeBoundaryCoupling's half-disc case, where the
		 * opposite sign converges to a different function.
		 */
		if ( exteriorDatumFunction && transferPath )
		{
			mfem::PositionFunction g = exteriorDatumFunction;
			exteriorDatumCoefficient =
				std::make_unique<mfem::PathTraceCoefficient>(
					*transferPath,
					[ g ]( mfem::Vector const &x ) { return -g( x ); } );

			// Whole, every time. See the declaration.
			fluxRhs = std::make_unique<mfem::LinearForm>();
			fluxRhs->Update( fluxFes.get(), rhs.GetBlock( 0 ), 0 );
			fluxRhs->AddBdrFaceIntegrator(
				new mfem::VectorBoundaryFluxLFIntegrator(
					*exteriorDatumCoefficient ),
				gammaHMarker );
			fluxRhs->Assemble();
		}

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
	                                         double psiBoundaryIn,
	                                         int *element, int *dof )
	{
		// BOTH, and psi_bnd used not to be here either. It is inert today --
		// setBoundaryFluxPoint() is NPC-only and this is the condensation's peak,
		// so psi_bnd is zero whenever this runs -- but it is the same latent trap
		// that cost the Jacobian its sign one screen up, and a signature that
		// cannot express psi_bnd is how that trap gets set again.
		normalisedSource->setNormalisation( psiAxisIn, psiBoundaryIn );

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
	void GradShafranovSolver::setPlasmaCurrent( double muZeroCurrent )
	{
		if ( !std::isfinite( muZeroCurrent ) || muZeroCurrent == 0.0 )
			throw std::invalid_argument(
				"meq::GradShafranovSolver::setPlasmaCurrent: mu0 * I_p must be "
				"finite and non-zero -- a zero target is the trivial branch "
				"asked for by name" );
		if ( orderingChoice != NonlinearOrdering::NPC )
			throw std::logic_error(
				"meq::GradShafranovSolver::setPlasmaCurrent: the current "
				"constraint is implemented for NonlinearOrdering::NPC only -- "
				"its row is a covector on the POTENTIAL, which is an unknown of "
				"the system only under NPC" );

		currentIsUnknown = true;
		targetMuZeroCurrent = muZeroCurrent;
		prepared = false;
	}

	double GradShafranovSolver::plasmaCurrentScale() const
	{
		return currentScaleValue;
	}

	double GradShafranovSolver::plasmaCurrent() const
	{
		return plasmaCurrentValue;
	}

	/*
	 * THE THREE LOOPS BELOW ARE ONE LOOP WITH THREE INTEGRANDS, and they are
	 * written out rather than shared because what differs is not only the
	 * integrand but where the answer goes: two are covectors on the potential
	 * block and one is a scalar.
	 *
	 * ALL THREE USE meq::SourceIntegrator'S OWN QUADRATURE RULE, and they have
	 * to: these are derivatives of the ASSEMBLED residual, not of the continuous
	 * one, so a different rule would differentiate a different function.
	 */
	namespace
	{
		/// The rule SourceIntegrator uses, so that every derivative of the
		/// assembled source term is taken on the same points it was assembled on.
		mfem::IntegrationRule const &sourceRule( mfem::FiniteElement const &el,
		                                         mfem::ElementTransformation &tr,
		                                         int extra )
		{
			return mfem::IntRules.Get( el.GetGeomType(),
			                           2*el.GetOrder() + tr.OrderW() + extra );
		}
	}

	double GradShafranovSolver::assemblePlasmaCurrent( mfem::Vector const &state ) const
	{
		if ( !normalisedSource )
			return 0.0;

		mfem::Mesh &mesh = *potentialFes->GetMesh();
		mfem::Array<int> dofs;
		mfem::Vector shape;
		mfem::Vector point;
		int const potentialStart = blockOffsets[ 1 ];
		double total = 0.0;

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::FiniteElement const &el = *potentialFes->GetFE( e );
			thread_local mfem::IsoparametricTransformation scratch;
			mesh.GetElementTransformation( e, &scratch );
			potentialFes->GetElementDofs( e, dofs );
			int const dof = el.GetDof();
			shape.SetSize( dof );

			mfem::IntegrationRule const &ir =
				sourceRule( el, scratch, sourceQuadratureExtra );
			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				scratch.SetIntPoint( &ip );
				el.CalcShape( ip, shape );
				scratch.Transform( ip, point );

				double psi = 0.0;
				for ( int j = 0; j < dof; ++j )
					psi += shape( j )*state( potentialStart + dofs[ j ] );

				// scaledF, not f: with coils present f() is the SUM and the
				// prescribed current is the plasma's alone.
				total += ip.weight*scratch.Weight()
				         *normalisedSource->scaledF( point( 0 ), point( 1 ), psi )
				         /point( 0 );
			}
		}
		return total;
	}

	void GradShafranovSolver::assembleCurrentColumn( mfem::Vector const &state,
	                                                 mfem::Vector &out ) const
	{
		out.SetSize( state.Size() );
		out = 0.0;
		if ( !normalisedSource )
			return;

		double const scale = normalisedSource->currentScale();
		if ( scale == 0.0 )
			return;

		mfem::Mesh &mesh = *potentialFes->GetMesh();
		mfem::Array<int> dofs;
		mfem::Vector shape;
		mfem::Vector point;
		int const potentialStart = blockOffsets[ 1 ];

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::FiniteElement const &el = *potentialFes->GetFE( e );
			thread_local mfem::IsoparametricTransformation scratch;
			mesh.GetElementTransformation( e, &scratch );
			potentialFes->GetElementDofs( e, dofs );
			int const dof = el.GetDof();
			shape.SetSize( dof );

			mfem::IntegrationRule const &ir =
				sourceRule( el, scratch, sourceQuadratureExtra );
			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				scratch.SetIntPoint( &ip );
				el.CalcShape( ip, shape );
				scratch.Transform( ip, point );

				double psi = 0.0;
				for ( int j = 0; j < dof; ++j )
					psi += shape( j )*state( potentialStart + dofs[ j ] );

				// F carries the scale linearly, so dF/d(scale) is F/scale --
				// and the residual's source term is -w F/r, so this is that
				// term divided by the scale. Exact, and one loop.
				double const derivative =
					normalisedSource->scaledF( point( 0 ), point( 1 ), psi )/scale;
				double const factor =
					-ip.weight*scratch.Weight()*derivative/point( 0 );
				for ( int j = 0; j < dof; ++j )
					out( potentialStart + dofs[ j ] ) += factor*shape( j );
			}
		}
	}

	void GradShafranovSolver::assembleCurrentNormalisationCorner(
		mfem::Vector const &state, double &againstAxis,
		double &againstBoundary ) const
	{
		againstAxis = 0.0;
		againstBoundary = 0.0;
		if ( !normalisedSource )
			return;

		mfem::Mesh &mesh = *potentialFes->GetMesh();
		mfem::Array<int> dofs;
		mfem::Vector shape;
		mfem::Vector point;
		int const potentialStart = blockOffsets[ 1 ];

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::FiniteElement const &el = *potentialFes->GetFE( e );
			thread_local mfem::IsoparametricTransformation scratch;
			mesh.GetElementTransformation( e, &scratch );
			potentialFes->GetElementDofs( e, dofs );
			int const dof = el.GetDof();
			shape.SetSize( dof );

			mfem::IntegrationRule const &ir =
				sourceRule( el, scratch, sourceQuadratureExtra );
			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				scratch.SetIntPoint( &ip );
				el.CalcShape( ip, shape );
				scratch.Transform( ip, point );

				double psi = 0.0;
				for ( int j = 0; j < dof; ++j )
					psi += shape( j )*state( potentialStart + dofs[ j ] );

				double dAxis = 0.0;
				double dBoundary = 0.0;
				if ( !normalisedSource->normalisationDerivatives(
					     point( 0 ), point( 1 ), psi, dAxis, dBoundary ) )
					return;

				double const w = ip.weight*scratch.Weight()/point( 0 );
				againstAxis += w*dAxis;
				againstBoundary += w*dBoundary;
			}
		}
	}

	void GradShafranovSolver::assembleCurrentRow( mfem::Vector const &state,
	                                              mfem::Vector &out ) const
	{
		out.SetSize( state.Size() );
		out = 0.0;
		if ( !normalisedSource )
			return;

		mfem::Mesh &mesh = *potentialFes->GetMesh();
		mfem::Array<int> dofs;
		mfem::Vector shape;
		mfem::Vector point;
		int const potentialStart = blockOffsets[ 1 ];

		for ( int e = 0; e < mesh.GetNE(); ++e )
		{
			mfem::FiniteElement const &el = *potentialFes->GetFE( e );
			thread_local mfem::IsoparametricTransformation scratch;
			mesh.GetElementTransformation( e, &scratch );
			potentialFes->GetElementDofs( e, dofs );
			int const dof = el.GetDof();
			shape.SetSize( dof );

			mfem::IntegrationRule const &ir =
				sourceRule( el, scratch, sourceQuadratureExtra );
			for ( int i = 0; i < ir.GetNPoints(); ++i )
			{
				mfem::IntegrationPoint const &ip = ir.IntPoint( i );
				scratch.SetIntPoint( &ip );
				el.CalcShape( ip, shape );
				scratch.Transform( ip, point );

				double psi = 0.0;
				for ( int j = 0; j < dof; ++j )
					psi += shape( j )*state( potentialStart + dofs[ j ] );

				// d/dx of int F/r: the plasma's own dF/dpsi against the shape
				// functions. NOT negated -- this is the constraint's gradient,
				// not a residual contribution.
				double const factor =
					ip.weight*scratch.Weight()
					*normalisedSource->scaledDFdPsi( point( 0 ), point( 1 ), psi )
					/point( 0 );
				for ( int j = 0; j < dof; ++j )
					out( potentialStart + dofs[ j ] ) += factor*shape( j );
			}
		}
	}

	void GradShafranovSolver::solveWithNormalisation()
	{
		if ( globalisationChoice != Globalisation::None )
			throw std::logic_error( "meq::GradShafranovSolver::solve: psi_ax as an unknown is implemented for Globalisation::None only -- the KINSOL paths drive a residual of their own and the Picard ones do not build a Jacobian at all" );
		// No ordering guard: both surviving orderings carry psi_ax. They carry it
		// DIFFERENTLY, and the difference is the whole content of this function --
		// see the border and the corner below.

		bool const npcOrdering = orderingChoice == NonlinearOrdering::NPC;

		/*
		 * FB-5. THIS FUNCTION NOW CARRIES THREE KINDS OF BORDER AND NOT ONE,
		 * and the arithmetic below is the general bordered elimination rather
		 * than the scalar it used to be:
		 *
		 *   psi_ax    one row, when a NormalisedSource is set
		 *   psi_bnd   one more, when setBoundaryFluxPoint() named a point
		 *   a_n       N more, when setExteriorCoupling() coupled the exterior
		 *
		 * They are not alike, and the differences are the whole content of the
		 * code: psi_ax's row needs an argmax and its column is DIFFERENCED;
		 * psi_bnd's dof is fixed at setup and its column is differenced too;
		 * the exterior rows are the transmission integrals, their corner block
		 * is DIAGONAL, and their columns are CONSTANT in the iterate because
		 * `a` reaches the residual only through a load term.
		 *
		 * With no normalised source the psi_ax and psi_bnd rows are absent and
		 * this is a pure exterior-coupled Newton, which is what FB-1's vacuum
		 * problem wants.
		 */
		bool const hasNormalisation = normalisedSource != nullptr;
		int const nModes = exteriorCoupling ? exteriorCoupling->modeCount() : 0;

		if ( nModes > 0 && !npcOrdering )
			throw std::logic_error(
				"meq::GradShafranovSolver::solve: the exterior coupling needs "
				"NonlinearOrdering::NPC. Its rows are the transmission "
				"integrals, which are a covector on the FLUX, and only under NPC "
				"is the flux an unknown of the system -- under the condensation "
				"it is recovered from the trace and the row would need the "
				"recovery's derivative, which DarcyHybridization does not expose" );
		if ( nModes > 0 && !usesNonlinearForms() )
			throw std::logic_error(
				"meq::GradShafranovSolver::solve: the exterior coupling needs a "
				"meq::Source rather than a coefficient, because NPC needs a "
				"non-linear form to build its operator on. A source whose f() "
				"ignores psi and whose dFdPsi() is zero is the vacuum case and "
				"costs one Newton step" );

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

		/*
		 * RE-ASSEMBLE THE RIGHT HAND SIDE AND REBUILD THE NPC OPERATOR.
		 *
		 * The transferred exterior datum reaches the system as a LOAD on the
		 * flux equation, and a load is assembled in prepare(). So a Newton step
		 * that moves `a` and does not come back here evaluates its next residual
		 * against the datum of the PREVIOUS step -- measured, and it does not
		 * diverge: the transmission constraints sit at 1e-17 while || R || falls
		 * by a factor of 2/3 per iteration, which is `a` and `x` chasing each
		 * other rather than a Newton converging. On an AFFINE problem.
		 *
		 * prepare() re-runs FormLinearSystem(), which re-finalises the
		 * hybridization, and mfem::DarcyNPCOperator caches what it found there
		 * when it was built -- so the operator is rebuilt with it. Without that
		 * the second residual evaluation reads through a stale handle and the
		 * process dies inside this function with nothing in the trace to say so.
		 *
		 * prepare() also zeroes the iterate, so every caller puts it back.
		 */
		auto reprepare = [ & ]()
		{
			prepare();
			if ( npcOrdering )
			{
				npc = std::make_unique<mfem::DarcyNPCOperator>(
					*darcy->GetHybridization(), blockOffsets, darcyRhs );
			}
		};

		mfem::Vector residual( n ), column( n ), y( n ), z( n ), scratch( n );
		mfem::Vector columnL( n ), zL( n ), currentRow( n );
		double currentIntegral = 0.0;
		double constraintL = 0.0;
		double currentAgainstAxis = 0.0;
		double currentAgainstBoundary = 0.0;

		double sB = 0.0;
		int boundaryDof = -1;
		if ( boundaryFluxIsUnknown )
		{
			// INTO THE FULL VECTOR, not into the potential space. peakAt() scans
			// [ blockOffsets[1], blockOffsets[2] ) and reports an index into the
			// unknown, so the second border has to be shifted the same way or it
			// reads the FLUX block instead. Measured before it was: psi_bnd came
			// back 7.6e-02 from the field at the limiter and FLAT under
			// refinement, which is the signature -- an O( h^{k+1} ) nodal
			// difference would have fallen by a factor of eight.
			boundaryDof = blockOffsets[ 1 ]
			              + nearestPotentialDof( boundaryFluxR, boundaryFluxZ );
			sB = psiBoundaryValue;
		}
		// psi_ax is border 0 always; psi_bnd and the current scale take the next
		// slots when they are unknowns, and the exterior modes follow them all.
		int const boundaryIndex = 1;
		int const currentIndex = boundaryFluxIsUnknown ? 2 : 1;
		int const nBorders = 1 + ( boundaryFluxIsUnknown ? 1 : 0 )
		                     + ( currentIsUnknown ? 1 : 0 );
		(void)boundaryIndex;

		// The scale starts where it was left, so a re-solve continues rather
		// than restarting, and the source is told before any residual is taken.
		double sL = currentIsUnknown ? currentScaleValue : 1.0;
		if ( normalisedSource && currentIsUnknown )
			normalisedSource->setCurrentScale( sL );

		auto fieldResidual = [ & ]( mfem::Vector const &state, double normalisation,
		                            mfem::Vector &out )
		{
			if ( normalisedSource )
				normalisedSource->setNormalisation( normalisation, sB );
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
				return recoverPeak( state, normalisation, sB, element, dof );

			if ( normalisedSource )
				normalisedSource->setNormalisation( normalisation, sB );

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
			if ( normalisationChoice == Normalisation::Decoupled || !normalisedSource )
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
		double peak = hasNormalisation ? peakAt( unknown, s, &argElement, &argDof )
		                               : 0.0;
		double constraint = hasNormalisation ? s - peak : 0.0;
		double constraintB = boundaryFluxIsUnknown
		                     ? sB - unknown( boundaryDof ) : 0.0;
		fieldResidual( unknown, s, residual );

		/*
		 * THE EXTERIOR BORDERS, ASSEMBLED ONCE.
		 *
		 * The rows are the transmission integrals, which depend on the geometry
		 * and the mode alone -- exteriorTransmissionRows() sweeps Gamma with
		 * mfem::ExtensionBoundaryQuadrature and contracts the extended flux
		 * against each mode. They are a covector on the FLUX BLOCK of the
		 * unknown. DarcyForm's flux block holds -q and the rows are built for
		 * flux(), which undoes that sign, so contracting the row against the
		 * unknown directly is the right thing and no negation belongs here.
		 *
		 * THE COLUMNS ARE CONSTANT AND THAT IS WHY THEY ARE OUT HERE. `a`
		 * reaches the residual only through the transferred datum, which
		 * arrives as a LOAD TERM on the flux equation and is linear in `a`; the
		 * operator never sees it. So one difference per mode is EXACT and it is
		 * taken once per mesh rather than once per Newton step, which is what
		 * makes the border cost N backsolves and nothing else.
		 *
		 * THEY ARE MEASURED RATHER THAN ASSEMBLED BY HAND, and the reason is
		 * the elimination. prepare() runs FormLinearSystem(), which moves the
		 * essential trace values into the right hand side, so the load vector as
		 * VectorBoundaryFluxLFIntegrator deposits it is NOT the vector the
		 * residual is measured against. Differencing the residual asks the
		 * question the residual answers. It costs one prepare() per mode, at
		 * setup; prepare() zeroes the iterate, so it is saved and put back.
		 */
		std::vector<mfem::Vector> exteriorRows, exteriorColumns, exteriorZ;
		if ( nModes > 0 )
		{
			std::vector<mfem::Vector> const rows =
				exteriorTransmissionRows( *exteriorCoupling );

			for ( int mode = 0; mode < nModes; ++mode )
			{
				// NO NEGATION, AND THE FIRST VERSION HAD ONE. The rows are built
				// to be contracted against flux(), which is the SIGN-CORRECTED
				// field; the unknown carries DarcyForm's raw block, which is -q.
				// So contracting the row against the unknown directly is already
				// what FB-1b writes as rows[ m ] . ( -flux() ), and negating
				// again gives a border that is right in magnitude and wrong in
				// sign -- which does not diverge, it just fails to converge.
				mfem::Vector row( n );
				row = 0.0;
				for ( int i = 0; i < rows[ static_cast<std::size_t>( mode ) ].Size(); ++i )
					row( i ) = rows[ static_cast<std::size_t>( mode ) ]( i );
				exteriorRows.push_back( row );
				exteriorZ.emplace_back( n );
			}

			mfem::Vector const savedIterate( unknown );
			std::vector<double> const savedCoefficients = exteriorCoefficientValues;

			// The baseline, a = 0, and one unit response per mode.
			auto residualAt = [ & ]( mfem::Vector &out )
			{
				reprepare();
				unknown = savedIterate;
				fieldResidual( unknown, s, out );
			};

			std::fill( exteriorCoefficientValues.begin(),
			           exteriorCoefficientValues.end(), 0.0 );
			mfem::Vector baseline( n );
			residualAt( baseline );

			for ( int mode = 0; mode < nModes; ++mode )
			{
				std::fill( exteriorCoefficientValues.begin(),
				           exteriorCoefficientValues.end(), 0.0 );
				exteriorCoefficientValues[ static_cast<std::size_t>( mode ) ] = 1.0;
				mfem::Vector response( n );
				residualAt( response );
				response -= baseline;
				exteriorColumns.push_back( response );
			}

			exteriorCoefficientValues = savedCoefficients;
			reprepare();
			unknown = savedIterate;
			fieldResidual( unknown, s, residual );
		}

		/// T_m = ( transmission integral of x )_m + blockEntry( m ) a_m.
		auto transmissionConstraint = [ & ]( mfem::Vector const &state, int mode )
		{
			double total = 0.0;
			mfem::Vector const &row = exteriorRows[ static_cast<std::size_t>( mode ) ];
			for ( int i = 0; i < n; ++i )
				total += row( i )*state( i );
			return total
			       + exteriorCoupling->blockEntry( ExteriorDtN::firstMode() + mode )
			         *exteriorCoefficientValues[ static_cast<std::size_t>( mode ) ];
		};

		std::vector<double> transmission( static_cast<std::size_t>( nModes ), 0.0 );
		for ( int mode = 0; mode < nModes; ++mode )
			transmission[ static_cast<std::size_t>( mode ) ] =
				transmissionConstraint( unknown, mode );

		// c at the starting iterate, computed whatever the coupling: it is both
		// the first step's column and the scale gamma. gamma converts a
		// perturbation of psi_ax into the units the trace residual is measured
		// in, and || c || is exactly that conversion. Normalisation::Decoupled
		// does not use the column but is given the same gamma, or the two
		// convergence histories would not be comparable -- which is the whole
		// point of having a control.
		mfem::Vector columnB( n );
		mfem::Vector zB( n );
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

		/*
		 * THE AUGMENTED NORM, over every constraint there is.
		 *
		 * gamma converts a perturbation of a BORDER unknown into the units the
		 * field residual is measured in, and it is frozen at the first iterate
		 * so that the printed history compares like with like -- a gamma
		 * recomputed each step would put the Jacobian's own variation into the
		 * convergence order and manufacture one. The exterior constraints get
		 * the same scale: T_m is a flux integral like the rest of the border and
		 * there is no second natural scale to give it.
		 */
		auto augmentedNorm = [ & ]( double fieldNorm, double cAx, double cBnd,
		                            double cCurrent,
		                            std::vector<double> const &cModes )
		{
			double total = fieldNorm*fieldNorm
			             + gamma*gamma*( cAx*cAx + cBnd*cBnd + cCurrent*cCurrent );
			for ( double v : cModes )
				total += gamma*gamma*v*v;
			return std::sqrt( total );
		};

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
			double const coldPeak = hasNormalisation
			                        ? peakAt( coldState, s, nullptr, nullptr ) : 0.0;
			std::vector<double> coldModes( static_cast<std::size_t>( nModes ), 0.0 );
			for ( int mode = 0; mode < nModes; ++mode )
				coldModes[ static_cast<std::size_t>( mode ) ] =
					transmissionConstraint( coldState, mode );
			reference = augmentedNorm( coldResidual.Norml2(),
			                           hasNormalisation ? s - coldPeak : 0.0,
			                           boundaryFluxIsUnknown
			                             ? sB - coldState( boundaryDof ) : 0.0,
			                           currentIsUnknown
			                             ? assemblePlasmaCurrent( coldState )
			                               - targetMuZeroCurrent : 0.0,
			                           coldModes );
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
			double const norm = augmentedNorm( residual.Norml2(), constraint,
			                                   constraintB, constraintL,
			                                   transmission );
			newtonResidualHistory.push_back( norm );
			newtonIterationCount = iteration;

			if ( norm <= target )
			{
				converged = true;
				break;
			}
			if ( !std::isfinite( norm ) || iteration == newtonMaxIterations )
				break;

			bool const coupled = hasNormalisation
			                     && normalisationChoice == Normalisation::Coupled;

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
				double const peakUp = recoverPeak( traceX, s + h, sB );
				double const peakDown = recoverPeak( traceX, s - h, sB );
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
					double const up = recoverPeak( traceX, s, sB );
					traceX( dof ) = saved - traceStep;
					double const down = recoverPeak( traceX, s, sB );
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
			// BOTH NORMALISATIONS, AND THE SECOND ONE WAS MISSING UNTIL
			// 2026-09-06. The one-argument overload forwards to
			// setNormalisation( psi_ax, 0 ), so this silently ZEROED psi_bnd for
			// the rest of the iteration -- and the NPC Jacobian immediately
			// below, and every plasma-current assembly further down, were built
			// against Psi = psi/psi_ax on a support of { psi > 0 } instead of
			// the iterate's own. Measured before the repair: the Newton
			// direction failed to solve its own linearised system by 109%, and
			// the current constraint came out with the OPPOSITE SIGN, so the
			// right-hand side drove the profile scale the wrong way and no
			// damping was a descent direction. Inert wherever psi_bnd = 0, which
			// is every fixed-boundary case in this tree -- which is why it
			// survived FB-3.
			if ( normalisedSource )
				normalisedSource->setNormalisation( s, sB );

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

			for ( int mode = 0; mode < nModes; ++mode )
				npcLinear.Mult( exteriorColumns[ static_cast<std::size_t>( mode ) ],
				                exteriorZ[ static_cast<std::size_t>( mode ) ] );

			/*
			 * THE BORDERED ELIMINATION, IN ITS GENERAL FORM.
			 *
			 * The system Newton solves is
			 *
			 *     [ J   C ] [ dx ]     [ -R ]
			 *     [ B   D ] [ dp ]  =  [ -G ]
			 *
			 * with p the border unknowns -- psi_ax, psi_bnd and the N exterior
			 * coefficients -- C their columns, B their rows and D the corner.
			 * Eliminating dx = -( y + sum_j dp_j z_j ) with y = J^-1 R and
			 * z_j = J^-1 c_j leaves a dense system of the border's own size:
			 *
			 *     M dp = f,    M_ij = D_ij - b_i . z_j,   f_i = b_i . y - G_i
			 *
			 * ONE FACTORISATION AND ONE BACKSOLVE PER BORDER, which is the whole
			 * argument for putting the coupling here rather than solving N + 1
			 * separate problems and adding them: superposition is exact only
			 * while the problem is linear, and a plasma source is not.
			 *
			 * THE SCALAR CASE IS KEPT SEPARATE AND THAT IS DELIBERATE. With one
			 * border the dense route would compute the same quotient by a
			 * different sequence of roundings, and HighBetaConvergence asserts
			 * psi_ax to every digit it printed before this generalisation
			 * existed. A refactor that moves the last bit of a published number
			 * is a refactor that has to be argued about; keeping the division is
			 * cheaper than the argument.
			 */
			int const nBorderTotal = nBorders + nModes;
			std::vector<double> step( static_cast<std::size_t>( nBorderTotal ), 0.0 );

			// b_i . v, for each border row.
			auto rowDot = [ & ]( int i, mfem::Vector const &v )
			{
				if ( i == 0 )
				{
					double total = 0.0;
					for ( int j = 0; j < borderDofs.Size(); ++j )
						total += border( j )*v( borderDofs[ j ] );
					return total;
				}
				if ( boundaryFluxIsUnknown && i == boundaryIndex )
					return coupled ? -v( boundaryDof ) : 0.0;
				if ( currentIsUnknown && i == currentIndex )
				{
					double total = 0.0;
					for ( int j = 0; j < n; ++j )
						total += currentRow( j )*v( j );
					return total;
				}

				mfem::Vector const &row =
					exteriorRows[ static_cast<std::size_t>( i - nBorders ) ];
				double total = 0.0;
				for ( int j = 0; j < n; ++j )
					total += row( j )*v( j );
				return total;
			};

			// The corner block, row by row.
			auto cornerEntry = [ & ]( int i, int j )
			{
				if ( i < nBorders || j < nBorders )
				{
					// THE CURRENT ROW IS NOT DIAGONAL. int F/r depends on
					// psi_ax and psi_bnd explicitly, through the normalisation
					// the profiles are evaluated at, so it has entries against
					// both. Leaving them out does not move the answer -- it
					// costs the quadratic rate, and it was measured doing
					// exactly that: a clean geometric contraction of about 0.8
					// per step where Newton should be quadratic.
					if ( currentIsUnknown && i == currentIndex && i != j )
					{
						if ( j == 0 )
							return currentAgainstAxis;
						if ( boundaryFluxIsUnknown && j == boundaryIndex )
							return currentAgainstBoundary;
						return 0.0;
					}
					if ( i != j )
						return 0.0;
					if ( i == 0 )
						return corner;
					// d( int F/r )/d( scale ) = ( int F/r )/scale, F being
					// linear in it. Not 1, which is psi_bnd's corner.
					if ( currentIsUnknown && i == currentIndex )
						return sL != 0.0 ? currentIntegral/sL : 0.0;
					return 1.0;
				}
				return i == j
				       ? exteriorCoupling->blockEntry(
				             ExteriorDtN::firstMode() + i - nBorders )
				       : 0.0;
			};

			// The constraint residual, row by row.
			auto constraintAt = [ & ]( int i )
			{
				if ( i == 0 ) return constraint;
				if ( boundaryFluxIsUnknown && i == boundaryIndex ) return constraintB;
				if ( currentIsUnknown && i == currentIndex ) return constraintL;
				return transmission[ static_cast<std::size_t>( i - nBorders ) ];
			};

			// z_j, the backsolved column.
			auto columnZ = [ & ]( int j ) -> mfem::Vector const &
			{
				if ( j == 0 ) return z;
				if ( boundaryFluxIsUnknown && j == boundaryIndex ) return zB;
				if ( currentIsUnknown && j == currentIndex ) return zL;
				return exteriorZ[ static_cast<std::size_t>( j - nBorders ) ];
			};

			if ( currentIsUnknown )
			{
				/*
				 * THE CURRENT CONSTRAINT, ALL THREE PIECES ANALYTIC.
				 *
				 *   G_I    = int F/r - mu0 I_p        the constraint
				 *   dR/dL  = ( source term )/L        F is LINEAR in the scale
				 *   dG/dx  = int ( dF/dpsi )/r phi_j  a covector on the potential
				 *   dG/dL  = ( int F/r )/L            for the same linearity
				 *
				 * None of it is differenced, which is the whole reason this is
				 * cheap: prescribing the current costs one more backsolve and
				 * three element loops, not a second factorisation.
				 */
				currentIntegral = assemblePlasmaCurrent( unknown );
				constraintL = currentIntegral - targetMuZeroCurrent;
				assembleCurrentColumn( unknown, columnL );
				assembleCurrentRow( unknown, currentRow );
				assembleCurrentNormalisationCorner( unknown, currentAgainstAxis,
				                                    currentAgainstBoundary );
				npcLinear.Mult( columnL, zL );
			}

			if ( boundaryFluxIsUnknown )
			{
				// dR/d psi_bnd, the second column, by the same central difference
				// the first one uses.
				double const hB = normalisationStep( sB );
				double const columnBase = sB;
				sB = columnBase + hB;
				fieldResidual( unknown, s, columnB );
				sB = columnBase - hB;
				fieldResidual( unknown, s, scratch );
				sB = columnBase;
				columnB -= scratch;
				columnB /= 2.0*hB;
				npcLinear.Mult( columnB, zB );
			}

			if ( nBorderTotal == 1 )
			{
				double const borderDotY = rowDot( 0, y );
				double const borderDotZ = rowDot( 0, z );
				double const denominator = corner - borderDotZ;
				if ( denominator == 0.0 || !std::isfinite( denominator ) )
					throw std::runtime_error( "meq::GradShafranovSolver::solve: the bordered Jacobian is singular in psi_ax" );
				step[ 0 ] = ( borderDotY - constraint )/denominator;
			}
			else
			{
				mfem::DenseMatrix dense( nBorderTotal );
				mfem::Vector right( nBorderTotal ), solved( nBorderTotal );
				for ( int i = 0; i < nBorderTotal; ++i )
				{
					right( i ) = rowDot( i, y ) - constraintAt( i );
					for ( int j = 0; j < nBorderTotal; ++j )
						dense( i, j ) = cornerEntry( i, j ) - rowDot( i, columnZ( j ) );
				}

				mfem::DenseMatrixInverse inverse( dense );
				inverse.Mult( right, solved );
				for ( int i = 0; i < nBorderTotal; ++i )
				{
					if ( !std::isfinite( solved( i ) ) )
						throw std::runtime_error( "meq::GradShafranovSolver::solve: the bordered Jacobian is singular in ( psi_ax, psi_bnd, a )" );
					step[ static_cast<std::size_t>( i ) ] = solved( i );
				}
			}

			double const deltaS = step[ 0 ];
			double const deltaB = boundaryFluxIsUnknown
			                      ? step[ static_cast<std::size_t>( boundaryIndex ) ] : 0.0;
			double const deltaL = currentIsUnknown
			                      ? step[ static_cast<std::size_t>( currentIndex ) ] : 0.0;

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
			double const savedB = sB;
			double const savedL = sL;
			std::vector<double> const savedCoefficients = exteriorCoefficientValues;

			double bestNorm = std::numeric_limits<double>::infinity();
			double bestDamping = 0.0;
			double damping = 1.0;
			bool accepted = false;

			for ( int trial = 0; trial < 12 && !accepted; ++trial, damping *= 0.5 )
			{
				double trialNorm = std::numeric_limits<double>::infinity();
				try
				{
					// `a` FIRST, then the right hand side it changes, and only
					// then the state -- prepare() zeroes the iterate, so the
					// order is not a preference.
					for ( int mode = 0; mode < nModes; ++mode )
						exteriorCoefficientValues[ static_cast<std::size_t>( mode ) ] =
							savedCoefficients[ static_cast<std::size_t>( mode ) ]
							+ damping*step[ static_cast<std::size_t>( nBorders + mode ) ];
					if ( nModes > 0 )
						reprepare();

					unknown = savedState;
					unknown.Add( -damping, y );
					unknown.Add( -damping*deltaS, z );
					s = savedS + damping*deltaS;
					if ( boundaryFluxIsUnknown )
					{
						unknown.Add( -damping*deltaB, zB );
						sB = savedB + damping*deltaB;
					}
					if ( currentIsUnknown )
					{
						unknown.Add( -damping*deltaL, zL );
						sL = savedL + damping*deltaL;
						if ( normalisedSource )
							normalisedSource->setCurrentScale( sL );
					}
					for ( int mode = 0; mode < nModes; ++mode )
						unknown.Add( -damping*step[ static_cast<std::size_t>( nBorders + mode ) ],
						             exteriorZ[ static_cast<std::size_t>( mode ) ] );

					peak = hasNormalisation
					       ? peakAt( unknown, s, &argElement, &argDof ) : 0.0;
					constraint = hasNormalisation ? s - peak : 0.0;
					constraintB = boundaryFluxIsUnknown ? sB - unknown( boundaryDof ) : 0.0;
					if ( currentIsUnknown )
						constraintL = assemblePlasmaCurrent( unknown ) - targetMuZeroCurrent;
					fieldResidual( unknown, s, residual );
					for ( int mode = 0; mode < nModes; ++mode )
						transmission[ static_cast<std::size_t>( mode ) ] =
							transmissionConstraint( unknown, mode );
					// BOTH constraints, or the line search is blind to the one it
					// is not told about and will happily accept a step that has
					// wrecked psi_bnd to improve psi_ax.
					// EVERY constraint, or the line search is blind to the one it
					// is not told about. constraintL was omitted when the current
					// border was added, which is the same mistake this comment
					// already warned against for psi_bnd.
					trialNorm = augmentedNorm( residual.Norml2(), constraint,
					                           constraintB, constraintL,
					                           transmission );
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
				// EVERY border, or the fallback step puts the fields back on a
				// damping the exterior coefficients and psi_bnd never received,
				// which is a state no equation in this system describes.
				for ( int mode = 0; mode < nModes; ++mode )
					exteriorCoefficientValues[ static_cast<std::size_t>( mode ) ] =
						savedCoefficients[ static_cast<std::size_t>( mode ) ]
						+ bestDamping*step[ static_cast<std::size_t>( nBorders + mode ) ];
				if ( nModes > 0 )
					reprepare();

				unknown = savedState;
				unknown.Add( -bestDamping, y );
				unknown.Add( -bestDamping*deltaS, z );
				s = savedS + bestDamping*deltaS;
				if ( boundaryFluxIsUnknown )
				{
					unknown.Add( -bestDamping*deltaB, zB );
					sB = savedB + bestDamping*deltaB;
				}
				if ( currentIsUnknown )
				{
					unknown.Add( -bestDamping*deltaL, zL );
					sL = savedL + bestDamping*deltaL;
					if ( normalisedSource )
						normalisedSource->setCurrentScale( sL );
				}
				for ( int mode = 0; mode < nModes; ++mode )
					unknown.Add( -bestDamping*step[ static_cast<std::size_t>( nBorders + mode ) ],
					             exteriorZ[ static_cast<std::size_t>( mode ) ] );
				peak = hasNormalisation
				       ? peakAt( unknown, s, &argElement, &argDof ) : 0.0;
				constraint = hasNormalisation ? s - peak : 0.0;
				constraintB = boundaryFluxIsUnknown ? sB - unknown( boundaryDof ) : 0.0;
					if ( currentIsUnknown )
						constraintL = assemblePlasmaCurrent( unknown ) - targetMuZeroCurrent;
				fieldResidual( unknown, s, residual );
				for ( int mode = 0; mode < nModes; ++mode )
					transmission[ static_cast<std::size_t>( mode ) ] =
						transmissionConstraint( unknown, mode );
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
		psiBoundaryValue = sB;
		if ( currentIsUnknown )
		{
			currentScaleValue = sL;
			plasmaCurrentValue = currentIntegral;
			if ( normalisedSource )
				normalisedSource->setCurrentScale( sL );
		}
		normalisationResidualValue = constraint;
		// Both, for the reason above: leaving the source at psi_bnd = 0 after the
		// solve makes every later evaluation -- postProcess(), the estimator, the
		// sampler -- read a different source from the one that was solved.
		if ( normalisedSource )
			normalisedSource->setNormalisation( s, sB );

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
		 * an unchecked back-substitution to reproduce fields MEQ already holds
		 * exactly is a way to lose them, not a way to confirm them.
		 */
		bool const fieldsAreState = usesNonlinearForms()
		                            && orderingChoice == NonlinearOrdering::NPC;

		// psi_ax has to be in the source before anything assembles, so that
		// prepare() and the first residual evaluation see the same normalisation.
		if ( normalisedSource )
			normalisedSource->setNormalisation( psiAxisValue );

		prepare();

		if ( normalisedSource || exteriorCoupling )
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
			 * NPC: THE UNKNOWN IS THE WHOLE SYSTEM, AND MEQ ALREADY HAD IT.
			 *
			 * mfem::DarcyNPCOperator is an Operator on ( q, psi, psihat )
			 * together -- Nguyen, Peraire & Cockburn eqs (14)-(18) -- with the
			 * Jacobian solved by hybridized elimination inside
			 * mfem::DarcyNPCSolver: reduce to the trace, solve there once, recover
			 * the local increments. Every element-local operation is ONE linear
			 * solve against ONE factorisation, and GetNumLocalNLIterations()
			 * staying at zero is the acceptance signal that it really is NPC.
			 *
			 * It costs MEQ almost nothing to give the fields to Newton, because
			 * `solution` has always been a three-block vector on `blockOffsets`
			 * with darcyFlux, potentialGf and traceGf MakeRef'd into it. That IS
			 * the NPC unknown, block for block. So the fields are already in
			 * place when Mult() returns and RecoverFEMSolution() drops out of this
			 * path entirely -- see the guard on it below, which is not an
			 * optimisation: upstream records ComputeSolution() as never having
			 * been exercised against an NPC solution, and running the
			 * condensation's element-local NON-LINEAR back-substitution over an
			 * answer that already satisfies the full system would be asking a
			 * question nobody has checked for an answer MEQ already has.
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
			// carries none -- MEQ imposes psi = g_D as an ESSENTIAL condition,
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
				// NAMED RATHER THAN LEFT TO A default:, AND THAT IS THE WHOLE
				// VALUE OF WRITING THEM OUT. These three never arrive here --
				// solve() sends PicardThenNewton to solveByPicardThenNewton() and
				// the other two to solveByPicard() before this block is reached --
				// but the switch had no arm for them and no default either, so
				// falling through left `nonlinear` null and the SetOperator()
				// below dereferenced it. The invariant lived in two places that
				// did not reference each other.
				//
				// A default: would silence -Wswitch permanently, which is exactly
				// what should not happen: naming them keeps the warning live, so
				// a SEVENTH globalisation added tomorrow fails to compile here
				// rather than reaching a null pointer.
				case Globalisation::AndersonPicard:
				case Globalisation::PicardOnly:
				case Globalisation::PicardThenNewton:
					throw std::logic_error(
						"meq::GradShafranovSolver::solve: a Picard globalisation reached "
						"the Newton block, which means the dispatch at the top of solve() "
						"no longer routes it -- these paths iterate a fixed point on the "
						"potential and build no Newton solver" );
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
			// both converged to round-off. MEQ asks for a Jacobian per step so
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
		// mass form M_p, and MEQ's Newton path puts the whole potential block on
		// Mnl_p -- see buildForms() -- so the local problem got no potential mass
		// and no face constraint, was singular, and was factored and solved
		// anyway. What came back was not a degraded psi* but 1e15, without a word,
		// with psi_h agreeing to six figures either way.
		//
		// MFEM now lifts the non-linear potential integrators as a Jacobian frozen
		// at the computed potential, checking that the face constraint's gradient
		// does not depend on the trace -- which MEQ's constant tau satisfies, and
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
		 * condition it detects is a defect in a library MEQ does not own and is one
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
