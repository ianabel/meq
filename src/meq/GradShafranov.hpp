#ifndef MEQ_GRADSHAFRANOV_HPP
#define MEQ_GRADSHAFRANOV_HPP

#include <memory>
#include <vector>

#include "mfem.hpp"

#include "Source.hpp"

/*
 * The HDG discretisation of the Grad-Shafranov operator, on MFEM's DarcyForm.
 *
 * meq solves the fixed-boundary problem
 *
 *     -div_bar( ( 1/r ) grad_bar( psi ) ) = F( r, z, psi ) / r    in Omega
 *                                     psi = g_D                   on Gamma
 *
 * with the flux q := ( 1/r ) grad_bar( psi ) carried as an unknown in its own
 * right, so that the magnetic field is obtained at the same order as psi rather
 * than one order down. In first-order form,
 *
 *     r q - grad_bar( psi ) = 0,     -div_bar( q ) = F/r.
 *
 * This is the LDG-H method of refs/HDG-GradShafranov.pdf eq (8), restated with
 * the block structure explicit as refs/HDG-GradShafranov-Adaptive.pdf eq (13):
 *
 *     ( r q_h, v )_Th + ( psi_h, div_bar v )_Th - < psihat_h, v.n >_dTh  = 0
 *     ( q_h, grad_bar w )_Th - < qhat_h.n, w >_dTh                       = ( F/r, w )_Th
 *     < qhat_h.n, mu >_dTh\Gamma_h                                       = 0
 *     psihat_h = g_D  on Gamma_h
 *
 * with the numerical flux qhat_h.n := q_h.n -+ tau( psi_h - psihat_h ); the sign
 * there is the subject of a long comment in the .cpp, and it is not the one the
 * papers print.
 *
 * TWO PATHS THROUGH THIS CLASS, chosen by which setSource() is called.
 *
 * When F does not depend on psi the problem is linear: F goes to the right hand
 * side as a coefficient, the trace system is a matrix, and one direct solve
 * finishes it. When F does depend on psi the problem is semi-linear and is
 * closed by Newton, with the source moved off the right hand side and into the
 * operator as a non-linear potential mass term. Both papers use an
 * Anderson-accelerated Picard iteration instead; the reasons for departing are
 * in CLAUDE.md, and the price is that every Source must supply dF/dpsi.
 *
 * TWO BOUNDARY REGIMES, chosen by whether setExtension() is called.
 *
 * Without it the domain is polygonal and fitted, Gamma_h == Gamma, and psi = g_D
 * is an essential condition on the trace: setBoundaryData() supplies g_D,
 * DarcyHybridization::SetEssentialBC eliminates those trace dofs, and that is the
 * whole of it. That is stages 2 and 4, and it stays the simpler configuration.
 *
 * With it the true boundary Gamma is a curved level set that the mesh does not
 * follow, the mesh is a polygonal subdomain D_h of the region Gamma encloses, and
 * the datum is carried from Gamma to Gamma_h = dD_h along transferring paths, by
 * the technique of Cockburn and Solano -- refs/HDG-GradShafranov-Adaptive.pdf
 * sections 2.1-2.2. The transferred datum
 *
 *     phi_h( x ) = g( a( x ) ) + int_sigma C E_h( u_h ) . m ds
 *
 * -- with C = r, u_h the flux as DarcyForm holds it, which is -q, and m the unit
 * tangent of the path from x on Gamma_h to a( x ) on Gamma -- depends on the
 * unknown flux through the second term, so it is not data to be eliminated; it
 * is a coupling into the operator. See setExtension() for how that
 * changes the assembly, and note that the change is structural rather than an
 * extra term.
 *
 * There is no local post-processing (stage 3) and no adaptivity (stage 6).
 *
 * WHAT IS AND IS NOT MEQ'S SIGN CONVENTION FOR q -- read before using flux().
 *
 * DarcyForm is built for the mixed system
 *
 *     ┌         ┐ ┌   ┐   ┌    ┐
 *     | Mu  -Bt | | u |   | bu |
 *     | -B  -Mp | | p | = | bp |
 *     └         ┘ └   ┘   └    ┘
 *
 * whose flux obeys u = -k grad p, the opposite sign to the Grad-Shafranov q.
 * The MFEM integrators that make the hybridization consistent -- the constraint
 * NormalTraceJumpIntegrator and the trace rows of HDGDiffusionIntegrator -- have
 * that sign baked in and take no scaling argument, so the flux block of the
 * assembled system necessarily holds -q, not q. That negation is undone once,
 * in solve(), and only in the separate GridFunction that flux() returns; the
 * block vector itself stays in DarcyForm's convention throughout, because the
 * Newton residual is assembled by DarcyForm and expects it there. Do not "fix"
 * one without the other.
 */

namespace meq
{

	/**
	 * A constant HDG stabilisation parameter tau.
	 *
	 * Both papers set tau = 1 and note that optimal order needs only tau = O(1).
	 * MFEM's HDGDiffusionIntegrator does not do that by default: its built-in
	 * stabilisation is ( beta +- ( alpha/2 )( u.n )/|u.n| ){ h^-1 Q }, scaled by
	 * the inverse local mesh size and by the diffusion coefficient. That is the
	 * LDG choice, and it is not harmless here -- measured on the Solov'ev
	 * benchmark it costs a full order in the flux, k instead of k+1, while psi
	 * still converges at k+1. See tests/convergence/SolovievConvergence.cpp,
	 * which carries the numbers.
	 *
	 * HDGDiffusionIntegrator::SetStabilization() is the designed way out.
	 * HDGStabilization::StabValue() divides the quadrature weight out before
	 * calling Eval() and multiplies it back afterwards, so returning a bare
	 * constant tau assembles exactly < tau psi, w >.
	 *
	 * IsConstant() is true, which means EvalGrad() is never called. That matters
	 * more than it looks: fem/darcy/bilininteg_hdg.hpp warns that omitting
	 * EvalGrad for a state-dependent stabilisation gives "no wrong answer, only
	 * slow Newton convergence -- a failure that survives a passing regression
	 * suite". A constant tau cannot fall into that hole, which is a further
	 * reason to prefer it now that the problem can be non-linear. Do not make
	 * tau solution dependent without supplying EvalGrad.
	 */
	class ConstantStabilization : public mfem::HDGStabilization
	{
		public:
			explicit ConstantStabilization( double tauIn );

			/// MFEM's spelling, from HDGStabilization.
			bool IsConstant() const override; // NOLINT(readability-identifier-naming)

			/// MFEM's spelling, from HDGStabilization. Every argument is ignored:
			/// that is the whole point of a constant tau.
			mfem::real_t Eval( mfem::real_t sDiff, mfem::real_t un, mfem::real_t u, // NOLINT(readability-identifier-naming)
			                   mfem::real_t uhat,
			                   mfem::ElementTransformation &tr ) const override;

			double tau() const;

		private:
			double tauValue;
	};

	/**
	 * The psi-dependent source, as a domain integrator on the potential space.
	 *
	 * This is the whole of what makes the problem semi-linear, and the only
	 * place dF/dpsi is used. It contributes
	 *
	 *     residual: -( F( r, z, psi_h ), w )/r
	 *     Jacobian: -( ( dF/dpsi )( r, z, psi_h ) w, v )/r
	 *
	 * on each element. The minus and the 1/r are the same pair that
	 * setSource( mfem::Coefficient & ) applies to the linear right hand side,
	 * and for a Source whose f() does not depend on psi the two paths assemble
	 * exactly the same numbers -- which is worth knowing, because it makes the
	 * Solov'ev benchmark usable as a cross-check of the Newton path.
	 *
	 * The sign is not a free choice. DarcyForm assembles the potential row as
	 * -B q - Mp psi = bp, and under hybridization the local solve is handed the
	 * negated datum, so the local potential residual reads
	 *
	 *     B u + D psi + E psihat = -bp.
	 *
	 * Moving -( F/r, w ) from bp into D is what puts the source under the
	 * Newton iteration, and it arrives with the sign it had on the right hand
	 * side. See the note in the .cpp for the arithmetic.
	 *
	 * The Jacobian is the exact derivative of the residual this same class
	 * assembles, evaluated on the same quadrature rule, so the two cannot drift
	 * apart through a rule change. What they can drift apart through is a wrong
	 * Source::dFdPsi, which is what tests/unit/SourceTests.cpp checks against a
	 * finite difference of f(), and what the finite-difference Jacobian check in
	 * tests/convergence/NewtonConvergence.cpp checks once more at the level of
	 * the assembled reduced operator.
	 */
	class SourceIntegrator : public mfem::NonlinearFormIntegrator
	{
		public:
			/// @param sourceIn  F and dF/dpsi. Borrowed; it must outlive the
			///                  integrator, which means outliving the solver.
			/// @param extraOrderIn  quadrature order added to 2k, to keep the
			///                  integration of an exponential in psi from being
			///                  what limits a measured rate.
			explicit SourceIntegrator( Source const &sourceIn, int extraOrderIn = 4 );

			/// MFEM's spelling, from NonlinearFormIntegrator.
			void AssembleElementVector( mfem::FiniteElement const &el, // NOLINT(readability-identifier-naming)
			                            mfem::ElementTransformation &tr,
			                            mfem::Vector const &elfun,
			                            mfem::Vector &elvect ) override;

			/// MFEM's spelling, from NonlinearFormIntegrator.
			void AssembleElementGrad( mfem::FiniteElement const &el, // NOLINT(readability-identifier-naming)
			                          mfem::ElementTransformation &tr,
			                          mfem::Vector const &elfun,
			                          mfem::DenseMatrix &elmat ) override;

		private:
			mfem::IntegrationRule const &rule( mfem::FiniteElement const &el,
			                                   mfem::ElementTransformation &tr ) const;

			Source const *source;
			int extraOrder;

			/// Scratch. Not thread safe, in the manner of every MFEM integrator.
			mfem::Vector shape;
	};

	/**
	 * The hybridized HDG Grad-Shafranov solver, linear or semi-linear.
	 *
	 * Owns the three finite element spaces, the DarcyForm built on them and the
	 * solution block vector. The mesh is borrowed and must outlive the solver;
	 * so must the source and boundary-data coefficients handed to setSource()
	 * and setBoundaryData().
	 *
	 * Usage, linear:
	 *
	 *     meq::GradShafranovSolver solver( mesh, order );
	 *     solver.setSource( fCoefficient );          // F( r, z ), not F/r
	 *     solver.setBoundaryData( psiCoefficient );  // psi on Gamma
	 *     solver.solve();
	 *     solver.potential();  solver.flux();
	 *
	 * Usage, semi-linear -- the only difference is the argument to setSource():
	 *
	 *     solver.setSource( source );                // a meq::Source
	 *     solver.solve();
	 *     solver.newtonResiduals();                  // the convergence history
	 *
	 * Usage on a curved Gamma, either of the above plus:
	 *
	 *     solver.setExtension( path, gammaHMarker ); // psi = 0 carried from Gamma
	 *
	 * The forms are built on the first call to solve() or prepare(), not in the
	 * constructor, because the two paths need different forms -- a linear
	 * potential mass and a source on the right hand side, or a non-linear
	 * potential mass carrying the source -- and EnableHybridization() takes what
	 * it finds at the moment it runs. So setSource() must be called before
	 * anything that assembles, and changing which overload is used afterwards is
	 * refused rather than silently ignored.
	 *
	 * Every boundary attribute of the mesh is Dirichlet: the fixed-boundary
	 * problem is an interior Dirichlet problem by construction, so there is no
	 * knob for that and no Neumann path to get wrong. setExtension() chooses HOW
	 * the condition is imposed on an attribute, never whether.
	 */
	class GradShafranovSolver
	{
		public:
			/// @param meshIn   the computational domain. Borrowed, and r > 0 is
			///                 required everywhere on it -- the operator carries a
			///                 1/r and the Solov'ev expansion carries a log r.
			/// @param orderIn  the polynomial degree k, used for all three spaces.
			///                 Hybridization removes the inf-sup condition that
			///                 would otherwise force them apart.
			/// @param tauIn    the HDG stabilisation, constant. tau = 1 is both
			///                 papers' choice and the measured optimum here; the
			///                 pre-port code used 5.0 with no recorded reason.
			GradShafranovSolver( mfem::Mesh &meshIn, int orderIn, double tauIn = 1.0 );

			GradShafranovSolver( GradShafranovSolver const & ) = delete;
			GradShafranovSolver &operator=( GradShafranovSolver const & ) = delete;

			/// The right hand side F( r, z ) of a source that does not depend on
			/// psi -- NOT F/r. The 1/r belongs to the weak form and is applied
			/// here, which keeps meq::Source free of it too (see Source.hpp).
			/// The problem is then linear and solve() does one direct solve.
			void setSource( mfem::Coefficient &fIn );

			/// The source F( r, z, psi ) of a semi-linear problem, with its
			/// derivative. Borrowed. The problem is closed by Newton, and
			/// Source::dFdPsi is what the Jacobian is built from -- an error
			/// there does not move the converged answer, it only wrecks, or
			/// silently slows, the convergence to it.
			void setSource( Source const &fIn );

			/// The Dirichlet datum g_D for psi on Gamma. Non-homogeneous data is
			/// the normal case: the level set psi = 0 is the plasma boundary, but a
			/// benchmark on a rectangle cut out of an exact equilibrium is not.
			void setBoundaryData( mfem::Coefficient &boundaryIn );

			/**
			 * Start Newton from @a psiGuess rather than from the Dirichlet data
			 * alone. Borrowed; it must outlive the next solve().
			 *
			 * WHY THIS EXISTS. Every source of
			 * refs/HDG-GradShafranov-Adaptive.pdf sections 4.2 to 4.5 satisfies
			 * F( r, 0 ) = 0, so with homogeneous Dirichlet data psi == 0 SOLVES
			 * the discrete problem. Newton starts from the Dirichlet data, lands
			 * on that branch, and stops in zero iterations with an identically
			 * zero residual -- which looks exactly like success. GS-1's
			 * Algorithm 2 opens `psi^0 ; // Non-trivial initial guess` for this
			 * reason and no other.
			 *
			 * WHERE THE GUESS ACTUALLY GOES, which is not where it looks.
			 * Newton's unknown is the TRACE: solve() runs it on the condensed
			 * system and the volume unknowns are recovered afterwards. So a guess
			 * written as psi( r, z ) has to reach M_h, and it does through an
			 * L2( e ) projection onto each face -- see projectOntoTrace() in the
			 * .cpp, and note that GridFunction::ProjectCoefficient does NOT do
			 * this: it loops over volume elements and never touches a face dof.
			 *
			 * The potential block is seeded too. That is not for Newton, which
			 * never reads it, but for MFEM's ELEMENT-LOCAL non-linear solves,
			 * which iterate from whatever the block vector holds. Whether it
			 * helps them is an open question -- see CLAUDE.md on the pressure
			 * pedestal -- and it is free to try.
			 *
			 * ORDER MATTERS. The Dirichlet datum is applied AFTER the guess, so
			 * the boundary condition always wins on essential dofs. A guess that
			 * disagrees with g_D there is simply overwritten rather than fought
			 * over.
			 *
			 * Ignored on the linear path, where the solve is direct and a
			 * starting point means nothing.
			 */
			void setInitialGuess( mfem::Coefficient &psiGuess );

			/// The same, from a potential computed elsewhere -- a previous solve
			/// on this or another mesh. Borrowed; it must outlive the next
			/// solve(). Evaluated through a GridFunctionCoefficient, so a guess
			/// on a DIFFERENT mesh needs the caller to have transferred it first.
			void setInitialGuess( mfem::GridFunction const &psiGuess );

			/// Forget the guess; the next solve() starts from the Dirichlet data
			/// alone, which is the default.
			void clearInitialGuess();

			/// True once setInitialGuess() has been called and not cleared.
			bool hasInitialGuess() const;

			/**
			 * Carry a homogeneous Dirichlet datum from the curved Gamma to the
			 * polygonal Gamma_h, by extension from the subdomain.
			 *
			 * @param pathIn          the transferring paths. Borrowed; it must
			 *                        outlive the solver. Any mfem::TransferPath
			 *                        will do -- LevelSetPath is the family the a
			 *                        priori analysis is written for, and
			 *                        VertexConePath is the general one meq's
			 *                        benchmark uses, because a flux surface has no
			 *                        closest-point map in closed form.
			 * @param gammaHMarkerIn  the boundary attributes of Gamma_h, the part
			 *                        of the mesh boundary that SubMesh had to
			 *                        generate. Attributes NOT marked here are
			 *                        treated as fitted and keep the essential
			 *                        trace condition, so a domain with both kinds
			 *                        of boundary works.
			 * @param lineOrderIn     order of the quadrature along a path;
			 *                        negative takes twice the element order plus
			 *                        two, which is HDGExtensionIntegrator's own
			 *                        default.
			 *
			 * WHAT THIS CHANGES, which is more than it looks.
			 *
			 * On a marked face psihat is no longer an unknown with an essential
			 * value. It is phi_h, and phi_h depends on the flux, so the datum
			 * splits: g( a( x ) ) is data and the line integral is an operator.
			 * For meq g == 0 on Gamma -- the plasma boundary IS the level set
			 * psi = 0 -- so the data half vanishes identically and only the
			 * operator half is left, which is
			 * refs/HDG-GradShafranov-Adaptive.pdf eq (14). That is why this takes
			 * no datum: a non-homogeneous g would need
			 * < g o a, v.n > on the flux right hand side as well, which
			 * mfem::PathTraceCoefficient supplies and nothing here calls for yet.
			 *
			 * The operator half goes on the flux mass form as an
			 * HDGExtensionIntegrator, where it is local to the element owning the
			 * face and the hybridization never sees it -- which is why the weak
			 * route costs nothing structural and the essential-trace route would.
			 *
			 * Two things come OFF a marked attribute in exchange. The HDG
			 * stabilisation does, so tau is zero on Gamma_h: it is
			 * < tau( psi_h - psihat_h ), w >, and psihat_h there is not the datum
			 * any more. Leave it on and the method loses an order at k = 1 and two
			 * at k = 2 -- measured, see
			 * tests/convergence/ExtensionConvergence.cpp. And the flux constraint
			 * does, so B's boundary face integrator is registered on the fitted
			 * attributes only: there is no trace unknown on Gamma_h to constrain.
			 *
			 * The trace dofs of those faces stay in the essential list all the
			 * same, at zero, because SetEssentialBC is still given every attribute.
			 * Measured, that is inert -- the answer does not change in any digit
			 * without it, and the reduced matrix has no zero row either way -- but
			 * it is what keeps the one combination that does NOT work from arising:
			 * with the flux constraint registered on Gamma_h AND nothing pinning
			 * the dofs it constrains, psihat_h becomes a free unknown answering
			 * < qhat_h.n, mu > = 0, a natural condition, and the error reaches
			 * 5e13.
			 *
			 * So trace() is NOT psihat on Gamma_h. The datum actually imposed
			 * there is phi_h, and it is never stored.
			 *
			 * miniapps/hdg/extension.cpp in the MFEM tree is the worked driver
			 * this follows.
			 */
			void setExtension( mfem::TransferPath &pathIn,
			                   mfem::Array<int> const &gammaHMarkerIn,
			                   int lineOrderIn = -1 );

			/// True once setExtension() has been called.
			bool isExtended() const;

			/// Which non-linear solver drives the trace system.
			enum class Globalisation
			{
				/// mfem::NewtonSolver: full steps, no line search. Quadratic when
				/// it works, and it is what every rate in the suite was measured
				/// with.
				None,
				/// KINSolver( KIN_LINESEARCH ): Newton with SUNDIALS' line search
				/// backtracking. For the stiff GS-2 sources whose element-local
				/// solves an undamped step drives out of their basin.
				LineSearch,
				/// KINSolver( KIN_NONE ): the same machinery taking full steps, so
				/// that a difference between None and LineSearch can be attributed
				/// to the line search rather than to SUNDIALS.
				KinsolNoLineSearch,
				/// KINSolver( KIN_FP ) with Anderson acceleration: the papers'
				/// own method. F is evaluated at the previous iterate, so the
				/// potential block is LINEAR and every element-local elimination
				/// is a linear solve -- which is the whole point, see CLAUDE.md
				/// under "Why meq's Newton struggles". Depth from
				/// setAndersonDepth().
				AndersonPicard,
				/// The same fixed point without acceleration, for the comparison
				/// that says whether Anderson is doing the work. Measured
				/// undamped Picard stalls on the pedestal, so expect this to.
				PicardOnly,
				/// Anderson-accelerated Picard to walk the iterate into Newton's
				/// basin, then plain Newton from there for the quadratic endgame.
				///
				/// **This is the route for a COARSE MESH.** On GS-2 section 4.3 at
				/// k = 1, h = 0.05 Newton from the Dirichlet ramp fails at 60,
				/// while from the converged Picard state it finishes in FOUR
				/// iterations at observed order 2.01, agreeing with Picard's own
				/// answer to 4.7e-10. That agreement is what makes it a handoff
				/// rather than a change of problem, and
				/// picardThenNewtonRecoversQuadraticOrder asserts it.
				///
				/// **It is not the only route to those cases, and usually not the
				/// one to prefer.** Sections 4.2, 4.3 and 4.5 are under-resolved
				/// rather than stiff: raw Newton solves all three in 7 to 17
				/// iterations once resolved, and both refining h and raising k cure
				/// them independently. Reach for refinement when it is available.
				/// What this is for is the initial solve of an adaptive run, which
				/// must happen before there is an estimator to refine with.
				///
				/// It does not rescue section 4.4, the current hole, which fails at
				/// every order and mesh tried up to k = 3, n = 48. That one's
				/// problem is the trivial branch, not the iteration.
				///
				/// Picard's job here is NOT to solve the problem. Section 4.5
				/// converges at both orders from a Picard state that never met its
				/// own tolerance, so this is a globalisation, not a two-solver
				/// pipeline. Stage 1 failing to converge is therefore not an error
				/// and does not stop stage 2.
				///
				/// It is not a cheap option: stage 1 spent 122 to 290 iterations on
				/// the cases above, each a full linear solve. Reach for it when
				/// Globalisation::None fails, not before.
				///
				/// **Do not replace the tolerance with an iteration budget.** The
				/// handoff is not monotone in Picard effort -- on 4.5 at k = 1,
				/// budgets of 400 and 3 converge while 40 and 10 fail, and on 4.3
				/// at k = 1 a budget of 3 diverges to 1e4. A budget tuned on one
				/// mesh will betray you on the next; Picard's own tolerance is the
				/// trigger that worked wherever it was reached.
				PicardThenNewton
			};

			/// True when the potential block is assembled non-linearly, which is
			/// the Newton paths and not the Picard ones.
			bool usesNonlinearForms() const;

			/**
			 * Choose the non-linear solver. Globalisation::None is the default and
			 * is what every convergence rate in the suite was measured with.
			 *
			 * KINSOL IS NOT QUITE A DROP-IN, despite mfem::KINSolver deriving from
			 * mfem::NewtonSolver, and the difference is silent. NewtonSolver::Mult(
			 * b, x ) solves oper( x ) - b = 0; KINSolver::Mult ignores its first
			 * argument entirely -- it is unnamed in the signature -- and solves
			 * oper( x ) = 0. meq's trace right hand side is not zero, so the
			 * residual has to be shifted before KINSOL sees it. solve() wraps it;
			 * see ShiftedResidual in the .cpp.
			 *
			 * @throws std::logic_error if a KINSOL option is asked for and MFEM was
			 *         built without MFEM_USE_SUNDIALS, rather than silently falling
			 *         back to an undamped Newton and reporting rates for a solver
			 *         nobody asked for.
			 */
			void setGlobalisation( Globalisation choice );

			/// Which solver solve() will use.
			Globalisation globalisation() const;

			/// Which solver eliminates the flux and potential on each element.
			/// This is a DIFFERENT iteration from the one setGlobalisation()
			/// chooses, and on a stiff source it is the one that fails: see
			/// CLAUDE.md, "The nonlinearity is inside the local solve".
			enum class LocalSolver
			{
				Newton,   ///< undamped, MFEM's LSsolveType::Newton
				Lbfgs,    ///< limited-memory BFGS, which line searches
				Lbb       ///< Barzilai-Borwein
			};

			/// Choose it. Newton is the default and is what every rate in the
			/// suite was measured with. **Inert under
			/// NonlinearOrdering::LineariseThenCondense**, where there is no local
			/// non-linear solve to configure.
			void setLocalSolver( LocalSolver choice );

			/// In which order the hybridization and the linearisation are applied.
			/// A DIFFERENT axis from setGlobalisation(): that picks the outer
			/// iteration, this decides what the outer iteration's residual costs.
			enum class NonlinearOrdering
			{
				/// Condense first. Eliminating flux and potential on an element is
				/// then itself a non-linear solve, one per element per residual
				/// evaluation -- and on the stiff GS-2 sources those are what
				/// fail. MFEM's default, and meq's until measured otherwise.
				CondenseThenLinearise,
				/// Linearise first: Newton on the full ( q, psi, psihat ) system,
				/// hybridizing the linear system that results. Every local
				/// operation is then a linear solve. This is how the method is
				/// defined -- Nguyen, Peraire & Cockburn, refs/HDG-NPC-2.pdf
				/// section 2.6, eqs (14)-(18) -- and the ordering meq wants, on
				/// the argument in CLAUDE.md under "Why meq's Newton struggles".
				LineariseThenCondense
			};

			/// Choose it. Needs an MFEM carrying
			/// DarcyHybridization::SetNonlinearOrdering; see CLAUDE.md on the
			/// meq-integration branch.
			void setNonlinearOrdering( NonlinearOrdering choice );

			/// The ordering solve() will use.
			NonlinearOrdering nonlinearOrdering() const;

			/// Damping for the Picard paths, in ( 0, 1 ]. Which knob it reaches
			/// depends on the path: KINSetDampingAA for AndersonPicard,
			/// KINSetDamping for PicardOnly. Setting both compounds them and is
			/// worse than either, so solveByPicard() sets exactly one.
			///
			/// **The default of 1.0 -- undamped -- is right for Anderson and wrong
			/// for plain Picard**, and that asymmetry is measured rather than
			/// chosen. On the section 4.2 pedestal at k = 1, h = 0.05, 500
			/// iterations allowed:
			///
			///     depth   w = 1.0        w = 0.5
			///       0     fails          converges, 248
			///       1     converges, 162 converges, 358
			///       2+    fails          fails
			///
			/// So PicardOnly needs w around 0.5 and will stall undamped; Anderson
			/// is best left alone. See setAndersonDepth() for the other surprise.
			void setPicardDamping( double damping );

			/// Anderson subspace depth for Globalisation::AndersonPicard.
			///
			/// **The default is 1, and HDG-GS-1's m = 2 does not work here.** That
			/// paper takes m = 2 on Toth & Kelley's evidence that there is no gain
			/// beyond m = 3, and this file originally defaulted to 2 for the same
			/// reason. Measured on the pedestal above, depth 1 converges in 162
			/// iterations and every depth from 2 to 10 fails at 500, damped or
			/// not. Whether that is this fixed point's conditioning or something
			/// about KINSOL's implementation is not established, and until it is,
			/// do not raise this expecting the papers' behaviour.
			///
			/// Ignored by every other Globalisation.
			void setAndersonDepth( int depth );

			/// Newton's stopping rule and iteration cap. Ignored on the linear
			/// path. The defaults are tight on purpose: a Newton iteration that
			/// stops early looks exactly like one that converges slowly, and this
			/// solver is measured on the shape of its residual history.
			void setNewtonControl( double relativeToleranceIn,
			                       double absoluteToleranceIn,
			                       int maxIterationsIn );

			/// True once setSource( Source const & ) has been called.
			bool isNonlinear() const;

			/// Assemble the forms and reduce to the trace system, without solving
			/// it. solve() calls this first; it is public so that the Jacobian can
			/// be checked against a finite difference of the residual it claims to
			/// differentiate, which is the one check that separates a wrong
			/// Jacobian from a wrong discretisation.
			///
			/// Afterwards reducedOperator() is the residual operator R, where the
			/// Newton residual is R.Mult( lambda ) - reducedRhs(), and
			/// R.GetGradient( lambda ) is the Jacobian the solve uses.
			void prepare();

			/// Assemble and solve. Both a source and boundary data must have been
			/// set. On the semi-linear path this runs Newton and fills
			/// newtonResiduals().
			void solve();

			/// psi_h in W_h. Valid after solve().
			mfem::GridFunction &potential();
			mfem::GridFunction const &potential() const;

			/// q_h = ( 1/r ) grad_bar( psi ) in V_h, in meq's sign convention.
			/// Valid after solve(); see the sign note at the top of this file.
			mfem::GridFunction &flux();
			mfem::GridFunction const &flux() const;

			/// psihat_h in M_h, the hybrid unknown. Valid after solve().
			mfem::GridFunction &trace();
			mfem::GridFunction const &trace() const;

			/**
			 * Build the post-processed potential psi*_h in P_(k+1), and the
			 * enriched flux and total flux that come with it.
			 *
			 * This is the local post-processing of Stenberg that
			 * refs/HDG-GradShafranov-Adaptive.pdf section 2.7 uses, and meq does
			 * not implement it: DarcyForm does. Reconstruct() takes the
			 * hybridized solution, projects the normally continuous total flux
			 * onto RT_k through the constraint equation, and then solves one
			 * small mixed problem per element on spaces one order higher, with
			 * the element average of psi*_h pinned to that of psi_h. That last
			 * part is eq (19b), and the average is imposed by
			 * ReconstructFluxAndPot() replacing one potential equation with it --
			 * which it does only when the potential mass form carries no domain
			 * integrator to use as a source instead. It does not on the linear
			 * path, so the plain Stenberg constraint is what is applied.
			 *
			 * WHAT THIS IS FOR. psi_h converges at k+1 and so does q, so an extra
			 * order in psi buys a magnetic-confinement calculation nothing --
			 * which is why CLAUDE.md records post-processing as dropped, and why
			 * neither GS paper implemented it for accuracy. It is needed for the
			 * residual estimator of eq (20), whose eta_1, eta_2, eta_4 and eta_5
			 * are all built on psi*_h rather than psi_h -- four of the five, not
			 * the three CLAUDE.md records; eta_4 is [[ psi*_h ]]. eta_2 differentiates the
			 * potential, and on psi_h that costs an order; measured in
			 * tests/convergence/EstimatorConvergence.cpp, which reports both.
			 *
			 * MEASURED, on the fitted Solov'ev benchmark: psi*_h converges at k+2
			 * for k = 1, 2, 3 -- 3.03, 4.03, 5.00 across the sequence -- so the
			 * library route delivers the superconvergence the paper wants and no
			 * hand-written local solve is needed. See EstimatorConvergence.cpp.
			 *
			 * IT USED TO BE REFUSED ON THE SEMI-LINEAR PATH, and the history is
			 * kept because the failure was silent and could return the same way.
			 * ReconstructFluxAndPot() read the LINEAR potential mass form and never
			 * looked at the non-linear one, so on meq's Newton path -- where the
			 * whole potential block including the HDG stabilisation lives on the
			 * non-linear form, of necessity, see buildForms() -- the local problem
			 * it built had no potential mass and no potential constraint. It did
			 * not abort. It returned numbers: psi* of 9.9e14, 8.4e15 and 3.9e14 on
			 * three successive meshes against 3.8e-6, 2.4e-7 and 1.5e-8 for the
			 * same problem solved linearly, with psi_h agreeing to six figures
			 * either way.
			 *
			 * MFEM fixed it, and MEASURING IT IS WHAT RETIRED THE REFUSAL -- not
			 * reading the fix, because a silent 1e15 is precisely what a code read
			 * cannot detect. Example 5 through Newton, L2( psi* ) over four dyadic
			 * meshes: rates 3.05 at k = 1, 4.05 at k = 2, 5.03 at k = 3, and 47x,
			 * 113x, 125x smaller than psi_h on the finest. k+2 on the non-linear
			 * path, in other words, exactly as on the linear one.
			 * NewtonConvergence.cpp's thePostProcessedPotentialSurvivesNewton
			 * asserts it, and eq (20) needs it: FOUR of the estimator's five terms
			 * are built on psi*, so while the refusal stood the adaptive loop was
			 * linear-only.
			 *
			 * IT DOES SURVIVE THE EXTENSION PATH, which was not expected.
			 * ReconstructFluxAndPot() lifts only the DOMAIN integrators of the flux
			 * mass form onto the enriched space, so the boundary-face
			 * HDGExtensionIntegrator that carries the transferred datum is dropped
			 * -- and yet psi* still converges at k+2 on the stage-5 benchmark: 2.62
			 * and 3.00 at k = 1, 3.46 and 3.90 at k = 2. The local problem is driven
			 * by the reconstructed total flux and the element average of psi_h, both
			 * of which already know about the extension, and that is apparently
			 * enough. Measured rather than assumed, in both directions: an earlier
			 * version of this comment asserted the opposite.
			 *
			 * What does NOT survive the extension path is eta_5, and for an
			 * unrelated reason -- psihat_h on Gamma_h is pinned rather than being
			 * the phi_h actually imposed. See
			 * meq::ResidualEstimator::setTransferredBoundary().
			 */
			void postProcess();

			/// True once postProcess() has been called since the last solve().
			bool isPostProcessed() const;

			/// psi*_h in P_(k+1). Valid after postProcess().
			mfem::GridFunction &postProcessedPotential();
			mfem::GridFunction const &postProcessedPotential() const;

			/// q*_h in [P_(k+1)]^2, in meq's sign convention -- the negation
			/// flux() applies is applied here too. Valid after postProcess().
			mfem::GridFunction &postProcessedFlux();
			mfem::GridFunction const &postProcessedFlux() const;

			/// The normally continuous total flux qhat_h in RT_k, in DarcyForm's
			/// sign convention because that is the one the constraint equation it
			/// is projected through is written in. Valid after postProcess().
			mfem::GridFunction &totalFlux();
			mfem::GridFunction const &totalFlux() const;

			mfem::FiniteElementSpace &fluxSpace();
			mfem::FiniteElementSpace &potentialSpace();
			mfem::FiniteElementSpace &traceSpace();

			/// The reduced trace system, valid after prepare(). On the linear path
			/// the operator is the assembled matrix; on the semi-linear path it is
			/// the non-linear DarcyHybridization operator itself, whose Mult() is
			/// the residual and whose GetGradient() differentiates that residual
			/// rather than the continuous equation. CEDRES++ rejected the
			/// continuous derivative deliberately -- see CLAUDE.md -- and this is
			/// how meq gets the discrete one for free.
			mfem::Operator &reducedOperator();

			/// The reduced right hand side and the reduced unknown, valid after
			/// prepare(). The unknown aliases the trace block of the solution and
			/// arrives carrying the Dirichlet data on the essential trace dofs.
			mfem::Vector &reducedRhs();
			mfem::Vector &reducedSolution();

			/// The trace dofs the Dirichlet condition is imposed on. On the
			/// semi-linear path the residual is masked to zero on these and the
			/// Jacobian carries a unit row, so a finite-difference check of the
			/// Jacobian must perturb only their complement.
			mfem::Array<int> const &essentialTraceDofs() const;

			/// The l2 norm of the non-linear residual at the start of every Newton
			/// iteration, plus the final one. Empty on the linear path. This is
			/// the thing to look at when a semi-linear run misbehaves: a history
			/// that grinds down linearly means the Jacobian disagrees with the
			/// residual, which no amount of mesh refinement will fix.
			std::vector<double> const &newtonResiduals() const;

			/// Iterations spent in stage 1 of Globalisation::PicardThenNewton.
			/// Zero on every other path. newtonIterations() and newtonResiduals()
			/// report stage 2, which is what an order assertion wants.
			int picardIterations() const;

			/// UMFPACK symbolic analyses performed during the last Newton solve,
			/// and numeric factorisations. The point of the pair is the RATIO:
			/// symbolic reuse is on, so a converged Newton run analyses the
			/// sparsity ONCE while refactorising once per iteration. Both are zero
			/// without SuiteSparse and on the linear and Picard paths.
			long symbolicFactorisations() const;
			long numericFactorisations() const;

			/// The number of Newton iterations the last solve took. Zero on the
			/// linear path. One fewer than newtonResiduals().size(), since that
			/// counts the residual at the initial guess too.
			int newtonIterations() const;

			/// L2 errors against a closed form, on a quadrature rule generous
			/// enough that it does not itself limit the measured rate.
			double potentialError( mfem::Coefficient &exact ) const;
			double fluxError( mfem::VectorCoefficient &exact ) const;

			/// The same for psi*_h, on a rule scaled to its own higher degree.
			/// Valid after postProcess().
			double postProcessedPotentialError( mfem::Coefficient &exact ) const;

			/// The number of globally coupled unknowns, that is the size of the
			/// system actually solved. For an HDG method that is the trace space,
			/// not the sum of the three.
			int numTraceDofs() const;

			int order() const;
			double tau() const;

		private:
			void buildForms();

			mfem::Mesh &mesh;
			int orderValue;

			std::unique_ptr<mfem::FiniteElementCollection> fluxColl;
			std::unique_ptr<mfem::FiniteElementCollection> potentialColl;
			std::unique_ptr<mfem::FiniteElementCollection> traceColl;

			std::unique_ptr<mfem::FiniteElementSpace> fluxFes;
			std::unique_ptr<mfem::FiniteElementSpace> potentialFes;
			std::unique_ptr<mfem::FiniteElementSpace> traceFes;

			std::unique_ptr<mfem::DarcyForm> darcy;

			ConstantStabilization stabilization;

			/// r, the weight of the flux mass form. The diffusion coefficient of
			/// the operator is 1/r, and DarcyForm's flux mass form wants its
			/// inverse -- verified by experiment, see the .cpp.
			mfem::FunctionCoefficient radius;

			/// -1/r, the factor the linear source is multiplied by to form the
			/// potential right hand side. Both halves of that are load bearing;
			/// see the .cpp.
			mfem::FunctionCoefficient negativeInverseRadius;

			mfem::Coefficient *linearSource;
			Source const *nonlinearSource;
			mfem::Coefficient *boundaryData;
			std::unique_ptr<mfem::Coefficient> potentialRhsCoeff;

			/// The Newton starting point, or null. Borrowed when it came in as a
			/// Coefficient; owned when setInitialGuess( GridFunction ) had to
			/// wrap one. Only one of the two is ever live.
			mfem::Coefficient *initialGuess;
			std::unique_ptr<mfem::Coefficient> ownedInitialGuess;

			/// Interpolate a coefficient onto the trace space, face by face.
			/// GridFunction::ProjectCoefficient cannot: it loops volume elements.
			/// See the .cpp.
			void projectOntoTrace( mfem::Coefficient &coeff,
			                       mfem::GridFunction &target ) const;

			/// The transferring paths, or null on the fitted path. Borrowed.
			mfem::TransferPath *transferPath;
			int extensionLineOrder;

			Globalisation globalisationChoice;
			LocalSolver localSolverChoice;
			NonlinearOrdering orderingChoice;
			int andersonDepth;
			double picardDamping;

			/// The iterate the frozen source reads, on the Picard paths. Lives in
			/// the potential space, which is also the fixed point's unknown.
			std::unique_ptr<mfem::GridFunction> picardIterate;
#ifdef MFEM_USE_SUITESPARSE
			/// The Picard path's trace solver, held across iterations so that its
			/// retained symbolic analysis has something to be reused by. Built on
			/// first use; see picardStep().
			std::unique_ptr<mfem::UMFPackSolver> picardSolver;
#endif

			/// F( r, z, picardIterate ), as the right hand side coefficient the
			/// linear path takes. Rebuilt with the spaces.
			std::unique_ptr<mfem::Coefficient> frozenSource;

			/// One Picard step: freeze F at @a in, assemble, solve, return the new
			/// potential in @a out. The fixed point map KINSOL iterates.
			void picardStep( mfem::Vector const &in, mfem::Vector &out );

			/// The Picard paths' whole solve: a KINSOL fixed point on the
			/// potential, optionally Anderson accelerated.
			void solveByPicard();

			/// Globalisation::PicardThenNewton: stage 1 then stage 2, re-entering
			/// solve() for each so that neither stage duplicates its body.
			void solveByPicardThenNewton();

			double newtonRelativeTolerance;
			double newtonAbsoluteTolerance;
			int newtonMaxIterations;
			int newtonIterationCount;
			/// Stage 1's count under Globalisation::PicardThenNewton, zero elsewhere.
			int picardIterationCount = 0;
			long symbolicFactorisationCount = 0;
			long numericFactorisationCount = 0;
			std::vector<double> newtonResidualHistory;
			/// The Picard iterate that seeds stage 2. It must be a COPY: the
			/// GridFunction overload of setInitialGuess() keeps a coefficient that
			/// only references its argument, and stage 2 overwrites potentialGf.
			std::unique_ptr<mfem::GridFunction> picardSeed;

			/// Every boundary attribute, marked. See the class comment. This is
			/// what the essential trace condition is imposed on, on both paths:
			/// on Gamma_h it pins dofs nothing references, on a fitted attribute
			/// it imposes g_D.
			mfem::Array<int> dirichletMarker;

			/// The attributes of Gamma_h, and its complement in dirichletMarker.
			/// Empty and equal to dirichletMarker respectively until
			/// setExtension() is called.
			mfem::Array<int> gammaHMarker;
			mfem::Array<int> fittedMarker;

			/// Empty, and a member rather than a local so that it outlives the
			/// DarcyHybridization that was handed a reference to it. A discontinuous
			/// flux space has no essential flux dofs -- the Dirichlet condition on
			/// psi is carried entirely by the trace.
			mfem::Array<int> essentialFluxTdofs;

			/// Four entries: flux, potential, trace, end. DarcyForm::GetOffsets()
			/// gives only the first three -- it does not know about the trace space
			/// until hybridization is enabled, and never grows to include it. The
			/// miniapp spells this DarcyOperator::ConstructOffsets(), which lives
			/// in miniapps/hdg and is not part of libmfem.
			mfem::Array<int> blockOffsets;

			mfem::BlockVector solution;
			mfem::BlockVector rhs;

			/// Two-block views -- flux and potential only -- over the first two
			/// blocks of solution and rhs. DarcyForm's own offsets stop at the
			/// potential, and BlockVector::operator= checks the block count, so
			/// handing it the three-block vector aborts inside ReduceRHS() on the
			/// semi-linear path ("Number of Blocks don't match"). The linear path
			/// survived it only because the corresponding checks there are
			/// MFEM_ASSERTs, compiled out of a release build. miniapps/hdg's
			/// DarcyOperator::ImplicitSolve() takes the same view for the same
			/// reason.
			mfem::BlockVector darcySolution;
			mfem::BlockVector darcyRhs;

			/// Aliases into solution. darcyFlux holds -q, see the file comment.
			mfem::GridFunction darcyFlux;
			mfem::GridFunction potentialGf;
			mfem::GridFunction traceGf;

			/// A copy of darcyFlux with the sign corrected, filled by solve().
			mfem::GridFunction fluxGf;

			/// The post-processed quantities, filled by postProcess(). Each owns
			/// the collection and space DarcyForm::Reconstruct() builds for it on
			/// first use -- RT_k for the total flux, and the primary collections
			/// cloned at k+1 for the rest -- which is why they are default
			/// constructed here and not given a space.
			mfem::GridFunction totalFluxGf;
			mfem::GridFunction enrichedFluxGf;
			mfem::GridFunction postProcessedGf;
			mfem::GridFunction enrichedTraceGf;

			/// The reduced trace system. traceX and traceB alias the trace blocks
			/// of solution and rhs; that aliasing is what makes FormLinearSystem()
			/// carry the essential trace values into the reduced problem, so it is
			/// not cosmetic. See the .cpp.
			mfem::OperatorHandle reduced;
			mfem::Vector traceX;
			mfem::Vector traceB;

			bool built;
			bool prepared;
			bool postProcessed;
	};

}

#endif // MEQ_GRADSHAFRANOV_HPP
