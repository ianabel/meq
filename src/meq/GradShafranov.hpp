#ifndef MEQ_GRADSHAFRANOV_HPP
#define MEQ_GRADSHAFRANOV_HPP

#include <memory>
#include <vector>

#include "mfem.hpp"

#include "ExteriorDtN.hpp"
#include "PlasmaComponent.hpp"
#include "Source.hpp"

/*
 * The HDG discretisation of the Grad-Shafranov operator, on MFEM's DarcyForm.
 *
 * MEQ solves the fixed-boundary problem
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
 * THREE PATHS THROUGH THIS CLASS, chosen by which setSource() is called.
 *
 * When F does not depend on psi the problem is linear: F goes to the right hand
 * side as a coefficient, the trace system is a matrix, and one direct solve
 * finishes it. When F does depend on psi the problem is semi-linear and is
 * closed by Newton, with the source moved off the right hand side and into the
 * operator as a non-linear potential mass term. Both papers use an
 * Anderson-accelerated Picard iteration instead; the reasons for departing are
 * in CLAUDE.md, and the price is that every Source must supply dF/dpsi.
 *
 * And when the source's profiles are functions of NORMALISED flux -- which is
 * how an equilibrium is actually specified -- psi on the magnetic axis is a
 * functional of the solution rather than data, and becomes an UNKNOWN of the
 * non-linear system: setSource( NormalisedSource &, double ) closes the trace
 * and psi_ax together by a bordered Newton. That is a structural change and not
 * an extra term, for the reason recorded on that overload.
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

			/**
			 * THE PLASMA'S CONNECTED COMPONENT, as a per-element test.
			 *
			 * XP-1. meq::NormalisedSource::insidePlasma() is pointwise on the
			 * value and has no connectivity, so `{ Psi > 0 }` is a level set
			 * rather than a plasma. The fill that fixes it needs element
			 * adjacency, which is the mesh's, so it is computed by the solver
			 * and handed here: on an element the fill did not reach, the plasma
			 * term is off and NormalisedSource::fOutsidePlasma() is what
			 * remains -- the coil term, or zero.
			 *
			 * @param component borrowed and may be null. Null, or a component
			 *        that has not been filled, means no connectivity test and
			 *        this class is bit-unchanged.
			 *
			 * The element number comes from ElementTransformation::ElementNo,
			 * which mfem::Mesh::GetElementTransformation sets and
			 * DarcyHybridization's own element workspaces carry through
			 * unchanged -- checked rather than assumed, in
			 * tests/convergence/PlasmaConnectivity.cpp.
			 */
			void setPlasmaComponent( PlasmaComponent const *component );

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

			/// Whether the fill reached @a element, and the constant true when
			/// no fill is live.
			bool elementCarriesPlasma( int element ) const;

			/// `F` at a point in @a element: the source's own `f()` where the
			/// fill reached, and what survives outside the plasma where it did
			/// not.
			double sourceValue( double r, double z, double psi, int element ) const;

			Source const *source;

			/// The same object as `source` when it is a NormalisedSource, and
			/// null otherwise. Held so that a masked-out element can ask what
			/// survives there without a dynamic_cast per element.
			NormalisedSource const *normalised;

			/// setPlasmaComponent(). Borrowed; null means no connectivity test.
			PlasmaComponent const *plasmaComponent = nullptr;

			int extraOrder;

			/// Per-point scratch, and a MEMBER only in a build that cannot
			/// thread. This used to read "not thread safe, in the manner of
			/// every MFEM integrator", which was true when it was written and
			/// is now false: MFEM converted its own -- HDGDiffusionIntegrator's
			/// dozen scratch members are guarded exactly like this
			/// (fem/darcy/bilininteg_hdg.hpp), and the installed library has
			/// MFEM_THREAD_SAFE = YES, so in MEQ's build they do not exist.
			///
			/// **AND THIS INTEGRATOR IS ON THE LOOP MFEM HAS JUST THREADED.**
			/// DarcyHybridization::MultNL() -- the residual and the Jacobian
			/// assembly, and so NPCResidual() and NPCGradient() -- now walks its
			/// elements in colour order under OpenMP when
			/// AssemblyMode::Threaded is asked for, and calls
			/// AssembleElementVector() and AssembleElementGrad() from inside it.
			/// Two elements of one colour share no face, but they shared THIS
			/// OBJECT: `shape` is resized and refilled per quadrature point, so
			/// a thread would compute `shape*elfun` against another element's
			/// basis. No crash, no error -- a wrong residual and a wrong
			/// Jacobian.
			///
			/// The Jacobian half is the worse one, and it is invisible to the
			/// obvious test: SourceIntegrator IS the whole semi-linear term, and
			/// CLAUDE.md's *A wrong Jacobian is invisible to a convergence
			/// table* records that a degraded dFdPsi leaves every error and
			/// every rate unchanged to six figures while costing Newton its
			/// order. A rate table cannot see this.
#ifndef MFEM_THREAD_SAFE
			mfem::Vector shape;
#endif
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

			/**
			 * A source specified in NORMALISED flux, with psi on the magnetic
			 * axis carried as an UNKNOWN of the non-linear system rather than as
			 * data.
			 *
			 * @param fIn            the source. Borrowed and MUTATED: the solver
			 *                       calls setNormalisation() on it before every
			 *                       residual evaluation, so it must outlive the
			 *                       solve and must not be shared with anything
			 *                       that reads it meanwhile.
			 * @param psiAxisGuessIn the starting value of psi_ax. It is a guess in
			 *                       the Newton sense -- see below, it matters.
			 *
			 * WHAT THIS CHANGES, which is structural rather than an extra term.
			 * The unknown becomes the pair ( lambda, psi_ax ): the trace, and one
			 * scalar. The system Newton closes is
			 *
			 *     R( lambda, psi_ax ) = 0        the hybridized trace residual
			 *     psi_ax - max psi_h  = 0        the normalisation, as an equation
			 *
			 * and the Jacobian is bordered,
			 *
			 *     [  A   c  ]        A = dR/dlambda,  c = dR/dpsi_ax
			 *     [  b^T d  ]        b = -d( max psi_h )/dlambda,  d = 1 - d( max psi_h )/dpsi_ax
			 *
			 * with c and b THE NON-LOCAL TERMS. They are non-local in the precise
			 * sense CEDRES++ (refs/CEDRES.pdf) means when it says a normalised
			 * profile "leads to non-local entries in the stiffness matrix": psi_ax
			 * is a functional of the whole solution, so a perturbation of the
			 * trace anywhere near the magnetic axis moves the source EVERYWHERE.
			 *
			 * WHY IT CANNOT BE A RANK-ONE UPDATE INSIDE THE ELEMENT BLOCKS, which
			 * is what a CG code would do. In an H^1 discretisation psi_ax is one
			 * entry of the global unknown and the Jacobian simply acquires a
			 * rank-one term. Hybridization eliminates flux and potential ELEMENT BY
			 * ELEMENT, and a term coupling every element to the one element holding
			 * the axis is exactly what that elimination cannot represent. The
			 * border is where it goes instead, and the bordered system is solved by
			 * block elimination -- two backsolves against one factorisation of A,
			 * so the extra unknown costs a backsolve and not a second matrix.
			 *
			 * WHY THE GUESS IS NOT OPTIONAL. With psi_ax held fixed the equation
			 * has a small solution -- the one Newton finds from zero, where the
			 * profile is never sampled beyond Psi ~ 1e-2 and is inert -- and a
			 * large one, which is the equilibrium. Only the large one can satisfy
			 * max psi_h = psi_ax. Adding the constraint removes the small branch
			 * from the SOLUTION SET but not from the iteration's reach, so the
			 * starting point still has to be on the right side of it:
			 * setInitialGuess() with a bump of about the right height, and a
			 * psiAxisGuessIn of about the right size. Dimensionally, for a
			 * pressure p = A Psi^nu on a box whose first Dirichlet eigenvalue is
			 * lambda_1, psi_ax is around sqrt( nu A / lambda_1 ).
			 *
			 * @throws std::logic_error if a source is already set or if the forms
			 *         are built. **No ordering is refused**, and the guard that
			 *         used to be here named an ordering MFEM has since deleted.
			 *         Both surviving orderings carry psi_ax, and they carry it
			 *         differently rather than one of them badly: under
			 *         NonlinearOrdering::NPC psi is an unknown, so the border row
			 *         is a UNIT VECTOR and the corner is exactly one, where under
			 *         CondenseThenLinearise psi is a function of the trace and
			 *         both have to be differenced. See solveWithNormalisation().
			 * @throws std::invalid_argument if psiAxisGuessIn is not finite or is
			 *         zero.
			 */
			void setSource( NormalisedSource &fIn, double psiAxisGuessIn );

			/**
			 * psi on the magnetic axis implied by @a trace: the potential
			 * recovered from it at the current normalisation, and the largest
			 * nodal value of that. This is exactly the functional the
			 * normalisation constraint is built on.
			 *
			 * It is public for the reason prepare() and reducedOperator() are:
			 * the border of the bordered Jacobian is d( this )/d lambda, and its
			 * being supported on the trace dofs of a single element is a claim
			 * about the hybridization that ought to be measured rather than
			 * argued. See theAxisSensitivityIsLocalToItsElement.
			 *
			 * @param element if not null, receives the element that attained it.
			 *
			 * Valid after prepare(). Leaves the solution blocks alone.
			 *
			 * @throws std::logic_error if psi_ax is not an unknown of this solver,
			 *         or if nothing has been prepared.
			 */
			double axisFlux( mfem::Vector const &trace, int *element = nullptr );

			/// How psi_ax is coupled to the field when it is an unknown.
			enum class Normalisation
			{
				/// The bordered Newton. psi_ax is a genuine unknown and the
				/// Jacobian carries the two non-local terms. The default, and the
				/// only one to run a calculation with.
				Coupled,
				/// The same iteration with the border DROPPED: c, b and d - 1 all
				/// taken to be zero. The step in psi_ax then reduces to
				/// psi_ax <- max psi_h and the trace step never sees the
				/// normalisation move -- which is psi_ax OUTSIDE the residual,
				/// done inside the same loop with everything else held fixed.
				///
				/// It exists so that what the non-local terms buy can be MEASURED
				/// rather than asserted, in the manner of
				/// ResidualEstimator's TraceComparison::Literal. Nothing else
				/// should use it: see theNonLocalTermsAreWhatMakeItConverge for
				/// what it does to the convergence.
				Decoupled
			};

			/// Choose it. Normalisation::Coupled is the default and is what every
			/// number in the suite was measured with. Ignored unless psi_ax is an
			/// unknown.
			void setNormalisationCoupling( Normalisation choice );

			/**
			 * WHAT `psi_ax` IS CONSTRAINED TO BE.
			 *
			 * `psi_ax` is what the profiles are NORMALISED by, so what it is
			 * defined as decides which equilibrium is reported. The two answers
			 * differ by `O( h )` in position and `O( h^2 )` in value, both
			 * independent of `k`, so on a refined high-order mesh they SEPARATE
			 * rather than converge -- measured on the finest Solov'ev mesh, the
			 * gap is 202x `psi_h`'s own L2 error at `k = 2` and 4204x at `k = 3`.
			 */
			enum class AxisConstraint
			{
				/**
				 * `psi_ax = max_j psi_h( x_j )` over the potential NODES.
				 *
				 * Cheapest by a wide margin: under NPC `psi` is part of the
				 * unknown, so the largest nodal value is literally one entry of
				 * it and the border row is exactly `-e_j` -- one entry, nothing
				 * assembled and nothing differenced.
				 *
				 * **AND NOTHING IN IT SAYS THE LARGEST NODAL VALUE IS A MAGNETIC
				 * AXIS**, which is the whole of FREE-BOUNDARY-PLAN.md section 11:
				 * `G = psi_ax - max psi_h = 0` is satisfied at machine zero by a
				 * spurious nodal spike exactly as it is by an axis, and three
				 * separate configurations were found reporting one.
				 */
				NodalMaximum,

				/**
				 * `psi_ax = psi_h( x* )` at the located magnetic axis, `x*` being
				 * a zero of `q_h`. The DEFAULT since 2026-09-07, and the
				 * principled one: it is what `psi_ax` MEANS.
				 *
				 * **IT COSTS NOTHING IN THE JACOBIAN, WHICH IS NOT OBVIOUS AND IS
				 * THE ENVELOPE THEOREM.** `x*` moves with the solution, so the
				 * chain rule gives
				 *
				 *     dG/dlambda = -[ dpsi_h/dlambda |_x*
				 *                     + grad( psi_h )( x* ) . dx* / dlambda ]
				 *
				 * and `grad_bar( psi ) = r q`, so `grad( psi_h )( x* ) = 0` at a
				 * zero of `q_h` **identically**. The position term vanishes, no
				 * sensitivity of the root find is needed, and the row is the
				 * potential shape functions of `x*`'s element evaluated at `x*` --
				 * `( k+1 )( k+2 )/2` entries, exact, one element. `-e_j` is the
				 * special case where `x*` lands on a node.
				 *
				 * **WHAT IT DOES COST** is a seeded root find per Jacobian, warm
				 * started from the previous iterate's axis, and the possibility
				 * that WHICH O-point is followed changes between iterations. The
				 * warm start is what keeps both small: it follows one axis rather
				 * than re-competing a field of them every step.
				 *
				 * **AND IT WEAKENS ONE DIAGNOSTIC, DELIBERATELY AND WITH A
				 * REPLACEMENT.** `CriticalPointFinder::checkAxis()` compares the
				 * normalised flux at a located O-point against 1, which under this
				 * choice is very nearly true by construction. It is not vacuous --
				 * it still catches a run with NO O-point, and one where the
				 * largest-`Psi` O-point is not the one the constraint followed --
				 * but the guard that carries the weight is now
				 * checkAxisSource(), which asks whether the toroidal current
				 * density is bounded on the symmetry axis and is untouched by
				 * this.
				 *
				 * **NPC ONLY.** Under the condensation `psi` is a function of the
				 * trace through every element's source, so the row would have to
				 * be differenced against `3( k+1 )` trace dofs with a root find
				 * inside each difference. It is REFUSED there rather than
				 * silently downgraded.
				 */
				LocatedAxis
			};

			/// Choose it. AxisConstraint::LocatedAxis is the default.
			///
			/// AxisConstraint::NodalMaximum is kept as the CONTROL rather than as
			/// a fallback -- every number in this tree published before
			/// 2026-09-07 was measured with it, so a table that moves is a
			/// statement about the definition rather than about the solver, and
			/// there has to be a way to take both readings on one problem.
			void setAxisConstraint( AxisConstraint choice );

			/// The value setAxisConstraint() last set.
			AxisConstraint axisConstraint() const;

			/// Where `psi_ax` was attained on the last solve, and whether it was
			/// a located axis or the fallback. Valid after solve().
			///
			/// `locatedAxis` is FALSE when AxisConstraint::NodalMaximum was asked
			/// for, and also when LocatedAxis was asked for and no O-point could
			/// be reached -- the annulus branch, where there is no closed surface
			/// to be an axis of. The two are distinguished by axisConstraint().
			bool axisWasLocated() const;
			double axisR() const;
			double axisZ() const;

			/// How the bordered Newton obtains its column `dR/ds`.
			enum class BorderColumn
			{
				/// Assembled from the source's own `dF/ds`, when it has one.
				/// Falls back to Differenced otherwise, so this is always safe.
				Analytic,
				/// A central difference of two full residual evaluations. What
				/// MEQ always did, kept so the two can be measured against each
				/// other -- and the only route under the condensation, whose
				/// residual is the reduced trace one rather than this assembly.
				Differenced
			};

			/**
			 * Choose it. The default is Analytic.
			 *
			 * IT IS NOT A PERFORMANCE KNOB, and the difference is largest
			 * exactly where it matters most. A differenced column costs two
			 * residual evaluations and so does the analytic one's fallback, so
			 * there is no speed in it. What there is:
			 *
			 *   - a differenced column FLOORS the iteration at the difference's
			 *     own accuracy, measured at about 3e-09 on a coupled half-disc
			 *     solve, where the residual then sits for as many iterations as
			 *     it is given;
			 *   - and with setPlasmaSupport() on it is worse than a floor,
			 *     because perturbing the normalisation MOVES THE PLASMA EDGE and
			 *     the two evaluations then have different supports.
			 *
			 * Differenced is kept because a control that can be switched off is
			 * how this project tells a repair from a coincidence.
			 */
			void setBorderColumn( BorderColumn choice );

			/// Which coupling solve() will use.
			Normalisation normalisationCoupling() const;

			/// True once setSource( NormalisedSource &, double ) has been called:
			/// psi_ax is an unknown and solve() runs the bordered Newton.
			bool normalisationIsUnknown() const;

			/// psi on the magnetic axis: the guess before solve(), the converged
			/// value after it. Zero unless normalisationIsUnknown().
			double psiAxis() const;

			/// The constraint residual psi_ax - max psi_h at the end of the last
			/// solve. It is the half of the augmented residual that says whether
			/// the normalisation is self consistent, and it is worth reporting
			/// separately from the trace residual because the two have different
			/// units. Zero unless normalisationIsUnknown().
			double normalisationResidual() const;

			/**
			 * WHETHER THE CONFINED SOURCE IS THE LEVEL SET OR THE CONNECTED
			 * PLASMA: stage XP-1 of FREE-BOUNDARY-PLAN.md section 10.6.
			 *
			 * meq::NormalisedSource::setPlasmaSupport() switches F off wherever
			 * `Psi <= 0`, and that test is POINTWISE -- it has no connectivity
			 * in it, so what it selects is a level set and not a plasma. Section
			 * 10.3 records the direction it is wrong in for a DIVERTED plasma,
			 * where the private flux region carries `Psi > 0` and gets a second
			 * current channel nobody asked for; MEASURED, it is wrong on
			 * ordinary LIMITER cases too, because a level of a field that is not
			 * monotone in radius cuts the domain into as many lobes as it likes.
			 * tests/convergence/FreeBoundaryCoupling.cpp's own half-disc, with
			 * no X-point anywhere, already converges to a `{ psi > psi_bnd }`
			 * with more than one component.
			 *
			 * Component is the default, because a connected plasma is what
			 * `ConfineToPlasma` MEANS. Pointwise is kept as the CONTROL, in the
			 * manner of Normalisation::Decoupled and ResidualEstimator's
			 * TraceComparison::Literal: a measurement of what the fill buys
			 * needs the configuration it replaces, and a measurement that cannot
			 * distinguish "the fix worked" from "there was nothing to fix" is
			 * not a measurement.
			 *
			 * **INERT UNLESS setPlasmaSupport() IS ON.** With the support off
			 * the domain IS the plasma and there is nothing to disconnect, so no
			 * mask is built and every existing solve is bit-unchanged.
			 */
			enum class PlasmaConnectivity
			{
				/// The pointwise test alone: `{ Psi > 0 }`, whatever its
				/// topology. The control, and what MEQ did before 2026-09-07.
				Pointwise,
				/// A face-neighbour flood fill over the elements carrying
				/// `Psi > 0`, seeded at the element holding psi_ax. The default.
				Component
			};

			/// Choose it. PlasmaConnectivity::Component is the default.
			void setPlasmaConnectivity( PlasmaConnectivity choice );

			/// Which one the next solve will use.
			PlasmaConnectivity plasmaConnectivity() const;

			/**
			 * Recompute the plasma's connected component from @a state.
			 *
			 * Called by solve() before every residual and every Jacobian
			 * assembly, so a caller normally never needs it; it is public
			 * because the MEASUREMENT of what the fill costs Newton has to drive
			 * it at a state of its own choosing, and because reporting on a
			 * converged answer means refreshing at that answer.
			 *
			 * @param state the full ( flux, potential, trace ) block vector on
			 *        this solver's offsets -- the NPC unknown -- or the
			 *        potential block alone, which is what the condensation path
			 *        has.
			 *
			 * Does nothing unless a NormalisedSource is set, its plasma support
			 * is on, and PlasmaConnectivity::Component is chosen.
			 *
			 * @throws std::invalid_argument if @a state is neither of those two
			 *         sizes. A vector of the wrong length would otherwise be
			 *         read as a potential and produce a plausible mask.
			 *
			 * NOT THREAD SAFE, and it does not need to be: it runs on the master
			 * thread before the element loop rather than inside it, so the mask
			 * meq::SourceIntegrator reads under AssemblyMode::Threaded is
			 * finished being written before any thread sees it.
			 */
			void refreshPlasmaComponent( mfem::Vector const &state );

			/// Elements the fill reached: the plasma. Zero when no fill is live.
			int plasmaComponentElements() const;

			/// Elements carrying `Psi > 0` at any of their potential dofs, over
			/// every component. The difference from plasmaComponentElements() is
			/// what the connectivity test removed, and it is the number that
			/// says whether the fill did anything.
			int plasmaCandidateElements() const;

			/// Components among those candidates. ONE means the pointwise test
			/// was already right on this configuration, which is the control
			/// every measurement of the fill needs.
			int plasmaComponentCount() const;

			/// True where the confined source is switched on, at element
			/// granularity. Always true when no fill is live, so the assembly
			/// tests it unconditionally.
			bool elementInPlasma( int element ) const;

			/// Which component an element is in, or -1 where it carries no
			/// plasma. For reporting on a fill rather than for the assembly.
			int plasmaComponentLabel( int element ) const;

			/// The seed's component number. For the same reason.
			int plasmaComponentSeedLabel() const;

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
			///
			/// AND THIS OVERLOAD SEEDS THE FLUX AS WELL, WHICH THE COEFFICIENT
			/// ONE STRUCTURALLY CANNOT. Under NonlinearOrdering::NPC `q` is an
			/// unknown, so a state carrying the right `psi` and `q = 0` is
			/// inconsistent in exactly the row that couples them: the flux row
			/// reads `( r q, v ) + ( psi, div v ) - < psihat, v.n >`, which at
			/// `q = 0` is the whole of `( grad psi, v )` and dominates the
			/// initial residual. A GridFunction can be differentiated and a bare
			/// Coefficient cannot, which is why the seed lives here.
			///
			/// It is a WEIGHTED projection, `( r q_h, v ) = ( grad psi_g, v )`
			/// element by element, and not `q = ( 1/r ) grad psi` interpolated at
			/// the nodes. Two reasons, and the second is the load-bearing one:
			/// the weighted form IS the flux row of the residual, so it makes the
			/// state consistent rather than merely close; and `1/r` is singular
			/// on the axis, which is where FB-A's domain reaches and where a
			/// nodal interpolation would divide a numerical zero by zero. The
			/// weight `r` removes the singularity instead of guarding it.
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
			 *                        VertexConePath is the general one MEQ's
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
			 * For MEQ g == 0 on Gamma -- the plasma boundary IS the level set
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

			/**
			 * The Dirichlet datum `phi_h` ACTUALLY IMPOSED on `Gamma_h`, as a
			 * Coefficient to be evaluated on a face of it.
			 *
			 * THE TRACE UNKNOWN IS NOT THIS. On `Gamma_h` the trace dofs are
			 * pinned to zero because nothing references them, and what the solve
			 * really imposes is `phi_h = g( a( x ) ) + L_e( q_h )` -- a value
			 * transferred along a path from the true boundary, which is never
			 * stored anywhere. Anything that wants to compare against the
			 * imposed condition, rather than against the zero standing in for
			 * it, has to rebuild it, and this is that. meq::ResidualEstimator's
			 * eta_5 is the caller that needs it.
			 *
			 * IT IS BUILT FROM THE RAW FLUX BLOCK AND NOT FROM flux(), WHICH IS
			 * THE ONE THING TO GET RIGHT HERE. `mfem::PathLiftCoefficient`
			 * evaluates the same `HDGExtensionIntegrator` the assembly used, so
			 * it wants the same convention that integrator was assembled
			 * against: `DarcyForm`'s block, which holds `-q`. `flux()` is a copy
			 * with the sign undone and would give a lifting of the wrong sign --
			 * a silent error, since the result stays smooth and merely converges
			 * to the wrong datum.
			 *
			 * @param g  the datum on the TRUE boundary as a function of position,
			 *           NOT its negation, and NOT the boundary-data Coefficient:
			 *           `mfem::PositionFunction` is evaluated at a bare point,
			 *           which a Coefficient cannot be. Defaults to zero, which is
			 *           MEQ's fixed-boundary problem and every case in the suite.
			 *
			 * @throws std::logic_error on the fitted path, where there is no path
			 *         family and no datum to transfer.
			 */
			std::unique_ptr<mfem::Coefficient> transferredDatum(
				mfem::PositionFunction g = mfem::PositionFunction() );

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
				/// under "Why MEQ's Newton struggles". Depth from
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
			 * oper( x ) = 0. MEQ's trace right hand side is not zero, so the
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
			/// suite was measured with.
			///
			/// **THE WHOLE SETTING IS INERT UNDER NonlinearOrdering::NPC**, which
			/// is MEQ's default -- solver type, preconditioner, cap and tolerance
			/// alike. NPC has no element-local non-linear solve to configure: the
			/// local work is one linear solve against one factorisation, and
			/// `GetNumLocalNLIterations()` staying at zero is how you check that.
			/// It bites only under CondenseThenLinearise.
			///
			/// **There it is load bearing, and this file used to say the
			/// tolerances did not matter.** MFEM's doxygen is blunt about why:
			/// `GetGradient()` is the Schur complement of the Jacobian at the
			/// fields the local solve reached, so it is the derivative of
			/// `Mult()` only as far as those fields solve the local problem. A
			/// fixed budget of two corrections shipped for a while and put the
			/// gradient 3e-04 out on a stiff source **that still converged** --
			/// which is the shape of every Jacobian defect this project has met:
			/// no wrong answer, only a wrong path to it. MEQ asks for 100
			/// iterations at rtol 1e-12, so it is bounded by the tolerance rather
			/// than by the cap.
			void setLocalSolver( LocalSolver choice );

			/**
			 * Quadrature order ADDED TO 2k for the semi-linear source term, on
			 * the meq::Source path. The default is 4 and every rate in this
			 * suite was measured with it.
			 *
			 * IT IS AN INSTRUMENT KNOB AND NOT A MODELLING ONE. Raising it
			 * cannot change the continuous problem, so a rate that MOVES when it
			 * moves is a rate limited by the rule rather than by the
			 * discretisation -- which is how this project tells the two apart
			 * everywhere else, by refining the instrument at fixed geometry.
			 *
			 * The one place it is not decoration is a source that STOPS inside
			 * an element. Free boundary's chi_{Omega_p} does exactly that, and a
			 * Gauss rule cannot see a discontinuity between its points however
			 * many of them there are: PlasmaEdgeConvergence measures psi_h's
			 * rate pinned to THREE FIGURES across a sweep of 4 to 20, which is
			 * what says the loss is the rule's BLINDNESS and not its resolution.
			 * What the sweep does move is psi*, 2.87 to 3.57 against its own
			 * bound of 3.5 -- so raising this is the cheap half of what a cut
			 * quadrature would have bought, and MEQ builds no cut quadrature for
			 * that reason.
			 *
			 * THERE IS A CEILING ON TRIANGLES AND NOTHING WARNS YOU. MFEM's
			 * symmetric triangle rules are exact and positive-weighted up to
			 * order 25 and switch at 26 to a construction whose least weight is
			 * -3.6e+01, reaching -1.9e+07 by order 64; a solve at 2k + 30
			 * returns errors of order 1e+3. Keep 2k + extraOrder <= 25, which at
			 * k = 3 means extraOrder <= 19. It is documented rather than
			 * enforced because a quadrilateral mesh has different rules, and
			 * theLossIsTheRulesBlindnessAndNotItsResolution pins the boundary so
			 * a change upstream fails loudly.
			 */
			void setSourceQuadratureOrder( int extraOrder );

			/// The value setSourceQuadratureOrder() last set.
			int sourceQuadratureOrder() const;

			/// Which non-linear method the hybridization is asked for. A
			/// DIFFERENT axis from setGlobalisation(): that picks the outer
			/// iteration, this decides what the outer iteration's unknown IS,
			/// and what one residual evaluation costs.
			///
			/// **A THIRD VALUE USED TO BE HERE AND MFEM DELETED IT.**
			/// `LineariseThenCondense` was an operator on the trace alone whose
			/// local blocks were eliminated against a retained linearisation,
			/// and it claimed to be the NPC method. It was not -- NPC's fields
			/// are Newton state, and a trace-only operator has nowhere to keep
			/// them, which is why that mode needed `lin_u`, `lin_p` and
			/// `lin_trace` as hidden state and why MEQ needed `Relinearised` to
			/// pair the residual with the gradient. Upstream measured it slower
			/// than the plain condensation on stiff problems and failing four
			/// configurations that one solves, and removed it. MEQ's
			/// `Relinearised` went with it. See
			/// ../mfem-hdg-dev/doc/HDG-ORDERING-API.md.
			enum class NonlinearOrdering
			{
				/// Condense first. Eliminating flux and potential on an element
				/// is then itself a non-linear solve, one per element per
				/// residual evaluation, and the outer unknown is the trace
				/// alone. MFEM's own default, and **MEQ's backup rather than
				/// MEQ's choice**.
				///
				/// It is kept, and is not merely legacy: it is the only route
				/// that is parallel, the only one that accepts an H(div) flux,
				/// and the only one whose reduced residual is an exact function
				/// of the trace -- which is what setLocalSolver()'s tolerance
				/// buys and what a differenced border needs when the fields are
				/// not state. PedestalConvergence measures the two against each
				/// other and needs this one for that.
				CondenseThenLinearise,
				/// Newton on the FULL ( q, psi, psihat ) system, with the
				/// Jacobian solved by hybridized elimination -- Nguyen, Peraire
				/// & Cockburn, refs/HDG-NPC-2.pdf section 2.6, eqs (14)-(18).
				/// `mfem::DarcyNPCOperator` and `mfem::DarcyNPCSolver`.
				///
				/// **THIS IS MEQ'S DEFAULT.** It is how the method is defined,
				/// and no paper in refs/ runs the other one: GS-1 and GS-2 avoid
				/// the question with Anderson-accelerated Picard, NPC linearises
				/// first.
				///
				/// **What it buys, and none of it is speed.** Every
				/// element-local operation is ONE linear solve against ONE
				/// factorisation, so `GetNumLocalNLIterations()` stays at zero
				/// -- which is the acceptance signal that this really is NPC and
				/// not a condensation wearing its name. The convergence test is
				/// on the full residual rather than on the trace alone, and a
				/// line search scales the fields and the trace together because
				/// they are one vector. Upstream's own caveat is worth
				/// repeating: **NPC is not automatically faster.** Its advantage
				/// is the UNIFORMITY of the local work, which is also what makes
				/// it the better batched or threaded workload, not fewer
				/// floating-point operations.
				///
				/// **What it costs MEQ is that the unknown is the whole
				/// system.** MEQ pays almost nothing for that, because
				/// `solution` was already a three-block
				/// { flux, potential, trace } vector on `blockOffsets` with
				/// every GridFunction MakeRef'd into it -- so the NPC unknown IS
				/// MEQ's solution vector, and `RecoverFEMSolution()` leaves the
				/// Newton path entirely rather than needing rework. The fields
				/// are already there when the solve returns.
				///
				/// **And it removes a trap rather than working around one.**
				/// `DarcyHybridization` freezes the element-local Newton's
				/// initial guess at `FormLinearSystem()` time, which cost the
				/// bordered Newton its correctness until `formSystem()` was
				/// factored out to re-form once per accepted step. NPC has no
				/// element-local non-linear solve, so there is no seed to go
				/// stale and no re-forming to do; see solveWithNormalisation().
				///
				/// **Two hard refusals**, both `MFEM_VERIFY` in `NPCCheck()`:
				/// an H(div) flux space, and `LocalOpType::FluxNL`. MEQ meets
				/// neither -- its flux space is L2 and its non-linearity is on
				/// the potential mass.
				NPC
			};

			/// Choose it. Needs an MFEM carrying `mfem::DarcyNPCOperator`; see
			/// CLAUDE.md on the MEQ-integration branch.
			void setNonlinearOrdering( NonlinearOrdering choice );

			/// The ordering solve() will use.
			NonlinearOrdering nonlinearOrdering() const;

			/// How the element loop that builds the reduced system runs.
			/// A THIRD axis, orthogonal to the two above: those decide what is
			/// computed, this decides who computes it. Purely a performance
			/// choice -- MFEM guarantees the two modes agree **bit for bit**,
			/// because the element-local arithmetic is per element and so
			/// reassociates nothing, and the scatter stays serial and in element
			/// order.
			enum class AssemblyMode
			{
				/// One thread. MFEM's default and **MEQ's**, unconditionally --
				/// see setAssemblyMode() for why an automatic gate was tried and
				/// removed.
				Serial,
				/// Thread EVERY element-local loop in DarcyHybridization, which is
				/// more than this option used to mean and is the reason to re-read it.
				///
				/// It once covered ComputeH() alone -- the factorisation of A, the
				/// Schur complement and its factorisation, and one local
				/// back-substitution per trace dof -- so it touched assembly and
				/// nothing else. MFEM has since threaded MultNL() as well: the
				/// residual and the Jacobian assembly, and therefore NPCResidual() and
				/// NPCGradient(), which is to say **every NPC step**. Plus ReduceRHS(),
				/// ComputeSolution() and the two RHS eliminations.
				///
				/// The two kinds of loop are threaded differently, and the difference
				/// is the scatter TARGET rather than the loop. ComputeH() scatters into
				/// an unfinalized mfem::SparseMatrix, which carries one current_row,
				/// one column-pointer scratch and one RowNode allocator for the whole
				/// matrix -- so two threads writing rows that are disjoint by
				/// construction still collide, and the failure is a hang rather than a
				/// wrong answer. That scatter stays serial and in element order, and it
				/// is 40-47% of NPCGradient(). MultNL() scatters into a Vector, where
				/// disjoint indices really are independent, so it is walked in element
				/// COLOUR order instead -- two elements of a colour share no face, so
				/// no trace dof takes two writes at once.
				///
				/// **Still bit for bit on both**, and on the colouring for a reason
				/// worth knowing: it changes the ORDER in which a face's two elements
				/// accumulate, and a + b == b + a exactly. MFEM notes this would NOT
				/// survive a trace space whose dofs are shared between faces -- an
				/// H1_Trace (EDG) one -- where a dof takes more than two contributions
				/// and associativity would start to matter. MEQ's is
				/// DG_Interface_FECollection, so it holds.
				///
				/// MFEM measures NPCResidual at 5.6-6.1x, NPCGradient at 2.6-3.3x and a
				/// whole NPC step at **1.9-2.1x** on eight threads. MEQ's own numbers,
				/// and what they mean for the default, are under setAssemblyMode().
				///
				/// **Requires MFEM_USE_OPENMP and MFEM_THREAD_SAFE, and MFEM ABORTS
				/// rather than falling back without either** -- deliberately, since a
				/// caller asking for this is asking a performance question and a silent
				/// serial loop is not an answer to it. MEQ therefore checks the build
				/// before passing it on.
				Threaded
			};

			/// Choose it. **Serial is the default and there is no automatic
			/// gate**, which was decided by measurement rather than caution.
			///
			/// Threaded is worth **1.15x to 1.33x** on a mesh assembled a few
			/// times, and it is **0.86x** -- a loss -- at one thread, where the
			/// chunk buffering pays for parallelism nobody asked for. A gate on
			/// `omp_get_max_threads() > 1` was therefore written, and then
			/// removed: `HighBetaConvergence`, which assembles a *small* system
			/// many times inside a bordered Newton, went from 21.5 s to 39 s
			/// under it. **1.8x slower, reproducibly.**
			///
			/// Mesh size does not separate the two cases -- HighBeta's meshes are
			/// 128 and 512 elements and 512 is where the isolated benchmark still
			/// showed a win. What separates them is how often assembly is called
			/// relative to everything else, which the solver cannot know. So the
			/// choice is the caller's: **take Threaded for a large mesh assembled
			/// a few times, leave it alone for a small one assembled hundreds of
			/// times.**
			///
			/// **A CALLER OBLIGATION COMES WITH IT, AND MFEM CANNOT CHECK IT.**
			/// Its own words: "Any integrator the caller installs -- a source term
			/// in the potential mass, a constraint integrator -- sits on this loop
			/// and must be thread-safe too. An integrator holding per-point scratch
			/// as a plain member will race, silently."
			///
			/// MEQ's side of that is discharged and asserted.
			/// meq::SourceIntegrator and meq::PoloidalFieldCoefficient guard their
			/// scratch on MFEM_THREAD_SAFE, exactly as MFEM's own HDG integrators
			/// now do; meq::ConstantStabilization holds one double and reads it;
			/// and every meq::Source and meq::Profile behind them is immutable
			/// during evaluation -- no interval memo in SplineProfile, no cached
			/// previous root in RotatingSource, which is why meq::maxSpecies is a
			/// compile-time cap.
			/// `threadedAssemblyReproducesSerialAssemblyOnANonlinearSource` is what
			/// says so, and it is the NONLINEAR case for a reason: the linear one
			/// takes meq's linear path, never installs meq::SourceIntegrator and so
			/// never reaches MultNL at all.
			///
			/// **The one hazard NOT closed is an exception leaving the loop.**
			/// meq::RotatingSource throws from f() and dFdPsi() when a species
			/// temperature goes non-positive or the quasineutrality root find fails
			/// -- reachable from a Newton iterate that overshoots. An exception
			/// escaping an OpenMP structured block is undefined behaviour, so on a
			/// rotating source that diagnostic is a clean throw under Serial and a
			/// terminate under Threaded. It is recorded rather than repaired
			/// because the loop is MFEM's: MFEM has the same exposure through its
			/// own MFEM_VERIFY under MFEM_USE_EXCEPTIONS, and inventing an error
			/// path the library does not support would replace a crash with a
			/// silent NaN. **Prefer Serial for a rotating source until the iterate
			/// is known good.**
			///
			/// Throws std::invalid_argument rather than letting MFEM abort the
			/// process when the library was built without OpenMP or without
			/// thread safety.
			void setAssemblyMode( AssemblyMode choice );

			/// The mode buildForms() will ask for.
			AssemblyMode assemblyMode() const;

			/// Whether this build can honour a given mode, so a caller can offer
			/// only what is there instead of catching to find out. The symmetric
			/// companion of traceSolverAvailable(), and it exists for the same
			/// reason: a driver reading a mode out of a configuration file wants to
			/// refuse it with a message about the build, not to relay an exception.
			///
			/// Serial is always available. Threaded needs MFEM_USE_OPENMP and
			/// MFEM_THREAD_SAFE, and this reports the BUILD rather than the thread
			/// count -- a build that can thread but is running at
			/// OMP_NUM_THREADS=1 answers true, because that is a run-time
			/// configuration and not a capability.
			static bool assemblyModeAvailable( AssemblyMode choice );

			/// Which direct solver factorises the hybridized trace system.
			///
			/// All three reach the same answer -- asserted, not assumed:
			/// `theTraceSolversAgree` pins them against each other to 1e-10 and
			/// they measure 1e-14 or better. So this is a **performance** choice
			/// and a licence choice, never a numerical one.
			enum class TraceSolver
			{
				/// SuiteSparse, `mfem::UMFPackSolver`, METIS ordering. **The
				/// default**, because it is the one every rate in the suite was
				/// measured with and the only one present in every build.
				UMFPack,
				/// oneMKL, `mfem::PardisoSolver`, `REAL_STRUCTURE_SYMMETRIC` --
				/// structurally symmetric on both paths, symmetric in value only
				/// on the fitted one. Needs `MFEM_USE_MKL_PARDISO`.
				///
				/// **Faster than UMFPack even single-threaded** -- 1.50x on the
				/// factorisation and 1.41x on the backsolve at 37,248 trace dofs
				/// -- and it is NOT the default anyway, because oneMKL's licence
				/// is not everybody's to accept and most builds do not have it.
				/// It scales a further 1.9x on MKL threads, which MEQ cannot
				/// currently spend: see CLAUDE.md, *Threaded MKL is a
				/// catastrophe*.
				Pardiso,
				/// NVIDIA cuDSS, `mfem::CuDSSSolver`, `NONSYMMETRIC` + `FULL`.
				/// Needs `MFEM_USE_CUDSS` **and an `mfem::Device` configured for
				/// CUDA before the solver is built** -- it reads its matrix and
				/// vectors through the device-aware accessors, which hand back
				/// host pointers otherwise.
				///
				/// **AND IT IS NOT WORTH TAKING ON ITS OWN, WHICH IS A STRONGER
				/// STATEMENT THAN THE TIMING ONE BELOW.** A device solver pays
				/// only if the data STAYS on the device, and MEQ's integrators
				/// and its scatter have no device kernels -- 58% to 70% of an NPC
				/// step by MFEM's own measurement. Its
				/// `doc/HDG-DEVICE-OFFLOAD.md`, under construction, gates the
				/// whole device path on the integrators for exactly this reason:
				/// doing the trace solve alone means copying the system across
				/// the bus once per iteration to accelerate a quarter of it,
				/// which is "plausibly slower than staying on the host
				/// throughout". `apps/meq.cpp` therefore refuses this choice
				/// from a config file; the library keeps it so that correctness
				/// on the device path can be checked at all.
				///
				/// **Correct, and not recommended on the strength of any timing
				/// taken here.** It agrees with UMFPack to 3.5e-14 from 9,408 to
				/// 148,224 trace dofs; its wall time on this machine varies by a
				/// factor of THIRTY between runs of the same binary on the same
				/// problem, so no ranking against the other two is possible
				/// locally. It is selectable so that correctness can be checked
				/// and so that somebody with a datacentre part can answer the
				/// performance question MEQ cannot.
				///
				/// Spelled the way NVIDIA spells it, which is why it breaks the
				/// UpperCamelCase rule for enum values: an external name keeps its
				/// author's capitalisation, as UMFPack and Pardiso do above. The
				/// rule this suppresses is house style, and house style does not
				/// get to rename other people's products.
				cuDSS // NOLINT(readability-identifier-naming)
			};

			/// Choose it. Throws std::invalid_argument when the library was built
			/// without the backing package, rather than silently falling back to
			/// a different solver -- a caller naming a solver has a reason.
			void setTraceSolver( TraceSolver choice );

			/// The trace solver solve() will use.
			TraceSolver traceSolver() const;

			/// Whether this build can honour a given choice. Lets a caller offer
			/// only what is there instead of catching to find out.
			static bool traceSolverAvailable( TraceSolver choice );

			/// The columns of `P`: the trace projection of each exterior mode onto
			/// `Gamma_h`, for the free-boundary coupling of
			/// FREE-BOUNDARY-PLAN.md section 4.3.
			///
			/// **WHAT THIS IS FOR.** Free boundary replaces the fixed-boundary datum
			/// on `Gamma_h` with `psihat|_Gamma_h = P a`, where `a` is the vector of
			/// Gegenbauer coefficients of the exterior expansion. `P` is
			/// `n_trace x N` and column `n` is one projection of
			/// `mfem::PathTraceCoefficient( path, C_n )` -- the mode evaluated at the
			/// FOOT of the transfer path, which is where `Gamma` actually is. Each
			/// returned vector is full trace length and is zero off `Gamma_h`.
			///
			/// **AND THE COLUMN IS CONSTANT IN THE ITERATE**, which is the claim
			/// section 4.3 makes and calls "an argument and not a measurement".
			/// `psihat` enters the flux row as `<psihat, v.n>`, the potential row as
			/// `<tau psihat, w>` and the trace row as `<tau psihat, mu>` -- linearly
			/// in all three -- while every non-linearity is `F( r, z, psi )`, which
			/// depends on `psi` and not on `psihat`. So `dF/da = ( dF/dpsihat ) P`
			/// cannot move, and `P` is built once per mesh rather than once per
			/// Newton step. This method takes no iterate for exactly that reason:
			/// **if it ever needs one, the claim has failed** and the caller should
			/// find out why rather than passing one in.
			///
			/// **THE MODES ARE INDEXED BY DEGREE FROM 2**, so column `i` carries
			/// degree `ExteriorDtN::firstMode() + i`. That is an off-by-TWO waiting
			/// to happen; see meq::ExteriorDtN.
			///
			/// Throws std::logic_error on the fitted path, where there is no
			/// `Gamma_h` and no transfer path to evaluate a mode at the foot of.
			std::vector<mfem::Vector>
			exteriorTraceColumns( ExteriorDtN const &exterior ) const;

			/// Impose a NON-ZERO datum on `Gamma_h`'s trace dofs, which is what
			/// `psihat|_Gamma_h = P a` of FREE-BOUNDARY-PLAN.md section 4.3
			/// means concretely.
			///
			/// **WITHOUT THIS THERE IS NO WAY TO PUT A DATUM ON `Gamma`, AND THE
			/// PLAN ASSUMED THERE WAS.** `setBoundaryData()` is projected against
			/// `fittedMarker` alone -- see `prepare()` -- because until free
			/// boundary nothing referenced `Gamma_h`'s trace values, so they were
			/// pinned to zero. `HDGExtensionIntegrator` supplies only the
			/// SOLUTION-dependent half of the transferred datum, the path integral
			/// of the flux, which is the whole of it exactly when `g` is
			/// homogeneous. So a fixed-boundary solve on a curved `Gamma` could
			/// only ever impose `psi = 0` there, and every extension study in this
			/// tree happens to want that.
			///
			/// Hand it `sum_m a_m P_m`, built from exteriorTraceColumns(). The
			/// vector is trace length and is ADDED into the trace, which is safe
			/// because a column is zero off `Gamma_h` by construction -- so this
			/// cannot disturb a fitted datum sitting on the same boundary, and on
			/// the half-disc the axis is exactly such a boundary.
			///
			/// The imposed condition is then `g( a( x ) )` plus the flux lifting,
			/// which is Cockburn-Solano's form with `g` no longer zero.
			///
			/// Passing an empty vector clears it. Throws std::invalid_argument on
			/// a length that is neither.
			void setExteriorDatum( mfem::PositionFunction g );

			/// The rows of `T`: the transmission condition of
			/// FREE-BOUNDARY-PLAN.md section 4.2, tested against each exterior
			/// mode. The Neumann half of the coupling, and the other half of
			/// exteriorTraceColumns().
			///
			/// **WHAT THIS IS FOR.** `P` says what the exterior expansion imposes
			/// on `Gamma_h`; on its own that leaves `a` undetermined, since
			/// nothing yet says the two fields agree in their normal derivative.
			/// Row `m` is that statement:
			///
			///     T_m( q, a ) = INT_Gamma E_h( q_h ).nu C_m dGamma
			///                 + a_m * exterior.blockEntry( m )
			///
			/// and the returned vector is the first term's derivative with respect
			/// to the flux unknowns -- which, the term being LINEAR in `q_h`, is
			/// also the term itself contracted with the iterate. Each vector is
			/// full SOLUTION length, not trace length, and is non-zero only on the
			/// flux dofs of the elements owning a `Gamma_h` face.
			///
			/// **THE INTEGRAL CARRIES NO 1/r AND THAT IS NOT AN OMISSION.** The
			/// exterior block is diagonal in the weight `dGamma/r`, which is what
			/// makes section 3.2 work at all, so the interior term has to be
			/// tested in that same weight or the two do not meet. It is: the
			/// transmission condition equates `(1/r) dpsi/dnu` across `Gamma`, and
			/// MEQ's `q` IS `(1/r) grad_bar( psi )` -- so testing `q.nu` against
			/// `C_m` in the PLAIN measure `dGamma` already carries the `1/r` the
			/// exterior side carries in its weight. Writing `dGamma/r` here would
			/// divide by the radius twice. The flux is the asset again, for the
			/// fourth time in this tree.
			///
			/// **AND THE RADIUS THAT WOULD BE WRONG IS THE FOOT'S, NOT THE NODE'S**
			/// -- see meq::GridSampler, where reading a band quantity at the foot
			/// rather than at the point cost a factor of 1.7e5. Here the question
			/// does not arise, because no radius is read.
			///
			/// **CONSTANT IN THE ITERATE, LIKE `P`, AND FOR A DIFFERENT REASON.**
			/// `P` is constant because `psihat` enters every row linearly; this row
			/// is constant because the extension `E_h` is a linear operator on the
			/// flux and `nu`, `C_m` and the measure are geometry. So it is built
			/// once per mesh, and this method takes no iterate. If it ever needs
			/// one, the same claim has failed.
			///
			/// **THE SIGN IS CHOSEN HERE RATHER THAN INHERITED.** The assembled
			/// flux block holds `-q` (see the file comment), so the extension of it
			/// is negated on the way out. Compare transferredDatum(), which must
			/// hand mfem::PathLiftCoefficient the RAW block precisely because that
			/// class re-runs the integrator the assembly used; nothing is re-run
			/// here, so the convention is written out.
			///
			/// **THE MODES ARE INDEXED BY DEGREE FROM 2**, as in
			/// exteriorTraceColumns(), and the same off-by-TWO is waiting.
			///
			/// Throws std::logic_error on the fitted path, for the same reason
			/// exteriorTraceColumns() does.
			std::vector<mfem::Vector>
			exteriorTransmissionRows( ExteriorDtN const &exterior ) const;

			/// Make `psi_bnd` an unknown too, pinned by `psi` at a prescribed
			/// point — FREE-BOUNDARY-PLAN.md's FB-3.
			///
			/// **THE PLASMA EDGE IS WHERE THE PROFILES STOP, AND IN A FREE
			/// BOUNDARY IT IS NOT KNOWN IN ADVANCE.** The profiles are functions
			/// of `Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd )`, and `psi_bnd`
			/// is the flux at the limiter contact or the X-point — a functional
			/// of the solution, exactly as `psi_ax` is. So it gets a border row
			/// of its own and the bordered Newton becomes 2x2.
			///
			/// `( r, z )` is the limiter contact. The constraint is `psi_bnd =
			/// psi_h` at the NEAREST POTENTIAL DOF to it, which is a definition
			/// rather than an approximation and is the same choice `psi_ax` makes
			/// in taking the largest nodal value: it is what makes the constraint
			/// differentiable in a form the border can use, and under NPC the row
			/// is then exactly `-e_j` and the corner exactly 1.
			///
			/// **NPC ONLY**, and refused otherwise. Under the condensation `psi`
			/// is a function of the trace through every element's source, so both
			/// the row and the corner would have to be differenced — which is
			/// possible and is not built, because §4.5 records that two of the
			/// three borders are exact under NPC and this is one of them.
			///
			/// Call with no point set — the default — and `psi_bnd` stays zero
			/// and the solve is the 1x1 it always was, arithmetically unchanged.
			void setBoundaryFluxPoint( double r, double z );

			/**
			 * PRESCRIBE THE PLASMA CURRENT AND SOLVE FOR THE PROFILE SCALE,
			 * which is how every production free-boundary code poses this.
			 *
			 * The profiles then give the SHAPE of the current and this gives its
			 * SIZE. meq::NormalisedSource::setCurrentScale explains why the
			 * alternative -- fixing the amplitude -- is a non-linear eigenvalue
			 * problem whose eigenvalue depends on a plasma region that is itself
			 * unknown, and why scaling the amplitude cannot fix it.
			 *
			 * @param muZeroCurrent  `mu0 * I_p`, NOT `I_p`.
			 *
			 * **THE ARGUMENT IS `mu0 I_p` DELIBERATELY, AND IT IS THE QUANTITY
			 * EVERYTHING HERE ALREADY SPEAKS IN.** Ampere's law through this
			 * solver reads `oint q.nu dGamma = -mu0 I_p`, outwardFlux() returns
			 * the left side, and the constraint below is assembled as
			 * `int F/r dOmega`, which IS `mu0 I_p`. Taking a current in amperes
			 * would mean knowing `mu0` here, and a `mu0` that disagreed with the
			 * source's own would scale two terms of one equation differently and
			 * converge, at full order, to a machine nobody described -- the trap
			 * CLAUDE.md records for why a coil block carries no `Mu0` key.
			 *
			 * **NPC ONLY**, and refused otherwise, for setBoundaryFluxPoint()'s
			 * reason: the row is a covector on the POTENTIAL, which is an
			 * unknown only under NPC.
			 *
			 * @throws std::invalid_argument if not finite or zero.
			 */
			void setPlasmaCurrent( double muZeroCurrent );

			/// The converged scale. One unless setPlasmaCurrent() was called.
			double plasmaCurrentScale() const;

			/// `int F_plasma/r` over the domain at the converged state, which is
			/// `mu0 I_p`. Zero unless setPlasmaCurrent() was called.
			double plasmaCurrent() const;

			/**
			 * COUPLE THE EXTERIOR TO THE SOLVE, so that the Gegenbauer
			 * coefficients on `Gamma` become unknowns of the same Newton rather
			 * than being recovered afterwards by superposition.
			 *
			 * FB-5. FB-1b already solves for them, and does it by running one
			 * full solve per mode and adding the answers -- which is exact and
			 * is available only because that problem is LINEAR. A plasma source
			 * is not, and superposition stops meaning anything the moment
			 * `F` depends on `psi`. This is the same system solved as one
			 * bordered Newton instead:
			 *
			 *     T_m( x, a ) = ( transmission integral of x )_m
			 *                   + blockEntry( m ) a_m = 0
			 *
			 * with the `psi_ax` and `psi_bnd` rows beside it, N + 2 borders
			 * against ONE factorisation.
			 *
			 * WHAT IT COSTS IS N + 2 EXTRA BACKSOLVES AND NO EXTRA
			 * FACTORISATION, which is the whole reason the border is the right
			 * structure: `DarcyNPCSolver::Mult()` is reduce, backsolve, recover,
			 * and `SetReuseSymbolic()` keeps the symbolic analysis across the
			 * Newton run.
			 *
			 * AND EVERY EXTERIOR COLUMN IS CONSTANT, which is worth more than it
			 * looks. `a` reaches the residual only through the transferred
			 * datum, which `setExteriorDatum()` deposits as a LOAD TERM on the
			 * flux equation -- linear in `a`, and independent of the iterate.
			 * So the N columns are assembled once per mesh, not once per Newton
			 * step, and what the border costs per step is the backsolves alone.
			 *
			 * @param exterior  the DtN. Borrowed; it must outlive the solve.
			 *
			 * Requires setExtension(): there is no `Gamma_h` to transfer from
			 * without it, and the fitted path imposes its datum through
			 * essential trace dofs where this needs a load term.
			 */
			void setExteriorCoupling( ExteriorDtN const &exterior );

			/// The converged Gegenbauer coefficients, in degree order from
			/// n = 2. Empty unless setExteriorCoupling() was called.
			std::vector<double> const &exteriorCoefficients() const;

			/**
			 * THE TRANSMISSION RESIDUAL PER ELEMENT: a BOUNDARY error indicator,
			 * and the one thing `meq::ResidualEstimator` structurally cannot see.
			 *
			 * FB-5 measured the gap rather than predicting it. Driving the
			 * coupled solve through the adaptive loop, `eta` fell 2.5109e-01 ->
			 * 3.2079e-02 over four cycles while the exterior coefficients sat at
			 * **1.3194e-03 at every cycle, to five digits**, and `Gamma_h` kept
			 * its 34 faces while the element count doubled. Under UNIFORM
			 * refinement the same quantity converges at 3.32, so it is not a
			 * mode-truncation floor -- the loop simply never refines there,
			 * because those 34 elements are 11% of the mesh and carry 0.00% of
			 * `eta^2`.
			 *
			 * AND `eta` IS RIGHT, WHICH IS WHY THIS IS A SEPARATE QUANTITY RATHER
			 * THAN A FIX TO IT. `eta` estimates the INTERIOR discretisation
			 * error, and `eta_5` on `Gamma_h` compares `psi*` against the datum
			 * actually imposed -- so it correctly reports the boundary as well
			 * resolved FOR THE INTERIOR PROBLEM. The coefficients are a BOUNDARY
			 * FUNCTIONAL: a transmission integral over `Gamma` reached by
			 * extension from `Gamma_h`. Nothing in `eta` measures it. Both are
			 * true at once.
			 *
			 * WHAT IS COMPUTED. The transmission condition is that the interior
			 * and exterior normal derivatives agree on `Gamma`, and
			 * setExteriorCoupling() imposes its PROJECTION onto the retained
			 * modes: `int_Gamma ( q_h.nu ) C_m dGamma + blockEntry( m ) a_m = 0`
			 * for each `m`. So at convergence the mismatch
			 *
			 *     d( x ) = q_h.nu( x ) - ( 1/r ) sum_n a_n symbol( n ) C_n( mu )
			 *
			 * is orthogonal to `C_2 .. C_{N+1}` and is NOT zero: what survives is
			 * the modes past the truncation and the discretisation error, which
			 * is exactly the quantity that froze. This returns
			 * `h_e int |d|^2 dGamma` per element -- the scaling `eta_3` uses for
			 * the flux jump, so that the numbers are commensurate with the other
			 * terms and Doerfler marking over their sum means something.
			 *
			 * @param exterior  the same DtN the solve was coupled to.
			 * @param out       sized to the element count; zero on every element
			 *                  not touching `Gamma_h`.
			 *
			 * @throws std::logic_error if there is no coupling, no `Gamma_h`, or
			 *         the solve has not run.
			 */
			void exteriorTransmissionResidual( ExteriorDtN const &exterior,
			                                   mfem::Vector &out ) const;

			/// The converged `psi_bnd`. Zero unless setBoundaryFluxPoint() was
			/// called. Valid after solve().
			double psiBoundary() const;

			/**
			 * IS THE TOROIDAL CURRENT DENSITY BOUNDED ON THE SYMMETRY AXIS?
			 *
			 * `meq::SourceIntegrator` assembles the load `-( F/r, w )`, and
			 * **`F/r` IS `mu_0 j_phi`**:
			 *
			 *     j_phi  =  r p'( Psi )  +  g g'( Psi ) / ( mu_0 r )
			 *
			 * so a finite current density on the axis REQUIRES `F( 0, z ) = 0`.
			 * `F = mu_0 r^2 p' + g g'` leaves only `g g'` there: `p'` is
			 * protected by its own `r^2` and `g g'` is not.
			 *
			 * **AND THE DISCRETE HALF IS WHY THIS NEEDS A CHECK RATHER THAN A
			 * COMMENT.** The CONTINUOUS problem is well posed whatever `F` does:
			 * the energy `int ( 1/r )|grad_bar psi|^2` forces its members to
			 * vanish faster than `r` at the axis -- the physical `psi ~ r^2` --
			 * and against such test functions `int ( F/r ) w` converges. The
			 * DISCRETE space is `L2` polynomials, free to be nonzero at `r = 0`,
			 * and against those the load functional is **unbounded**. The
			 * quadrature is the only thing making it finite, so the answer
			 * depends on the RULE and not on the mesh -- measured, sweeping
			 * setSourceQuadratureOrder() at fixed `h` moves `psi_h` on the axis
			 * while `h`-refinement does not fix it, and `psi_h` develops an
			 * `O( 1 )` layer along the WHOLE axis at a value that can exceed the
			 * true peak.
			 *
			 * **WHICH `Psi` THE AXIS SITS AT IS THE WHOLE OF IT, AND IT IS FB-3
			 * THAT OPENS THE TRAP.** `psi( 0, z ) = 0` exactly -- `psi` is the
			 * poloidal flux through a circle of radius `r`, which vanishes with
			 * the area -- so `Psi_axis = -psi_bnd/span`. On a FIXED boundary
			 * `psi_bnd = 0`, the axis sits at `Psi = 0`, and every profile in
			 * this tree vanishes there; that is why nothing met this before.
			 * setBoundaryFluxPoint() makes `psi_bnd` an unknown, it comes out
			 * POSITIVE, and the axis is then at NEGATIVE `Psi` -- in the VACUUM,
			 * where the physics is `g = const` so `g g' = 0`, and where an
			 * unconfined profile EXTRAPOLATES and returns a current instead.
			 *
			 * **THE REPAIR IS `NormalisedSource::setPlasmaSupport()`**, which
			 * sets `F = 0` wherever `Psi <= 0` -- the statement that the vacuum
			 * carries no current. For a domain reaching `r = 0` with `psi_bnd`
			 * free it is a PRECONDITION rather than an option, in the same sense
			 * as `j >= 1` at the plasma edge.
			 *
			 * A mesh that does not reach the axis is unaffected and reports
			 * `reachesAxis == false`: there is no `1/r` to be unbounded.
			 */
			struct AxisSourceCheck
			{
				/// Whether any potential node sits at `r = 0`. FALSE is not a
				/// failure -- it is the ordinary case for a fitted rectangle,
				/// and every other field is then meaningless.
				bool reachesAxis = false;

				/// The largest `| F |` over the nodes AT `r = 0`, evaluated at
				/// `psi = 0` -- which `psi( 0, z )` is EXACTLY, for any
				/// axisymmetric field with bounded `B`. See the implementation
				/// for why the iterate's own `psi_h` is the wrong thing to ask:
				/// it measures the layer rather than its cause, and cannot
				/// separate a healthy run from a failing one.
				double worstOnAxis = 0.0;
				double worstR = 0.0;
				double worstZ = 0.0;

				/// The largest `| F |` anywhere, as the scale the one above is
				/// judged against. A bare tolerance on `F` would be a statement
				/// about the units the source is written in.
				double sourceScale = 0.0;

				/// worstOnAxis/sourceScale, and whether it is within tolerance.
				double relative = 0.0;
				bool bounded = true;

				/**
				 * `Psi` ON THE SYMMETRY AXIS, AND WHETHER THE PLASMA CONTAINS IT.
				 *
				 * `psi( 0, z ) = 0` exactly, so
				 * `Psi_axis = -psi_bnd/( psi_ax - psi_bnd )`, and
				 * `insidePlasma()`'s own test makes the axis part of the plasma
				 * when that is POSITIVE -- which needs `psi_bnd` and the span to
				 * carry opposite signs, i.e. a NEGATIVE `psi_bnd` on an ordinary
				 * positive span.
				 *
				 * **FOR A TOKAMAK THAT IS NOT A LARGE ERROR, IT IS THE WRONG
				 * TOPOLOGY.** A tokamak plasma is a torus about `R_0 > 0` and its
				 * symmetry axis is in the vacuum, so a support containing `r = 0`
				 * describes current threading the machine's own centre line, and
				 * no refinement makes it into the equilibrium that was asked for.
				 * It is a separate and worse statement than `bounded` above: that
				 * one says the load carries a `1/r` the quadrature is papering
				 * over, this one says the answer is not an equilibrium of the
				 * intended kind at all.
				 *
				 * **BUT IT IS NOT WRONG FOR EVERY DEVICE, AND THAT IS WHY THE
				 * ESCAPE CLAUSE IS PHYSICS RATHER THAN A HEDGE.** A **levitated
				 * dipole** has plasma right up to the axis, and so does a
				 * **magnetic mirror** -- and NEITHER HAS A TOROIDAL FIELD, so
				 * `g` vanishes identically in both. That is not a coincidence:
				 * it is the same fact twice. `B_phi = g/r` has to be finite on
				 * the axis, so a device whose plasma reaches `r = 0` cannot carry
				 * a toroidal field there, and `g == 0` is what makes the
				 * configuration admissible in the first place.
				 *
				 * **SO THE TEST IS `g g' == 0` AND NOT "IS THIS A TOKAMAK".**
				 * `F( 0, z, . )` is `g g'` and nothing else -- `p'` is killed by
				 * its own `r^2` -- so a source with `g g' == 0` identically puts
				 * no current on the axis whatever `Psi` reads there, and
				 * `j_phi = r p'` vanishes with `r` regardless. A dipole or a
				 * mirror passes; a tokamak whose `psi_bnd` has gone negative does
				 * not. `../geq`, the rotating-mirror wrapper this tree already
				 * compares against, sets `g g' == 0` unconditionally for exactly
				 * this reason.
				 *
				 * `sourceVanishesOnAxis` is taken over a SPREAD of `psi` rather
				 * than at one value: asking only at `psi = 0` cannot tell
				 * `g g' == 0` from `g g'( Psi_axis ) == 0` by luck, and under a
				 * moving support it would read zero for any profile at all.
				 */
				double normalisedFluxOnAxis = 0.0;
				bool axisInsidePlasma = false;
				bool sourceVanishesOnAxis = false;
			};

			/**
			 * Evaluate the above at the current iterate. Valid after solve(),
			 * and cheap: one pass over the potential dofs.
			 *
			 * @param tolerance how large `| F |` on the axis may be relative to
			 *        `| F |`'s own scale. Default 1e-6, which is not tuned: a
			 *        source that vanishes on the axis does so to ROUND-OFF -- the
			 *        clean case measures 1e-15 relative -- and one that does not
			 *        measures 1.5e-01. There is nothing in between to calibrate
			 *        against, so the default sits six orders below the failure
			 *        and nine above the floor.
			 */
			AxisSourceCheck checkAxisSource( double tolerance = 1.0e-6 ) const;

			/// The outward flux of `q` through the true boundary `Gamma`:
			/// `oint_Gamma q.nu dGamma`.
			///
			/// **THIS IS AMPERE'S LAW AND IT IS THE SHARPEST WHOLE-ASSEMBLY
			/// CHECK MEQ HAS.** Since `q = ( 1/r ) grad_bar( psi )`, the
			/// integrand is `( 1/r ) dpsi/dn`, and integrating the equation over
			/// the enclosed region gives
			///
			///     oint_Gamma ( 1/r ) dpsi/dn dl = -mu0 * I_enclosed
			///
			/// exactly, with no discretisation anywhere in the statement. So it
			/// ties the assembled operator, the source, the boundary condition
			/// and the transfer to ONE number that is known in advance --
			/// `meq::CoilSet::totalCurrent()` for a prescribed current, and
			/// `tests/unit/CoilsTests.cpp` already pins the identity itself on
			/// the exact field to 3.3e-11, so a discrepancy here is the SOLVE.
			///
			/// Swept over `Gamma` with `mfem::ExtensionBoundaryQuadrature`, the
			/// same routine and the same signed weight the transmission rows
			/// use, with the mode dropped. So it inherits their caveats: the
			/// weight is signed on purpose, and the rule must resolve the foot
			/// map -- see setTransmissionQuadratureOrder().
			///
			/// Valid after solve(). Throws std::logic_error on the fitted path,
			/// where there is no band and `Gamma` is `Gamma_h`.
			double outwardFlux() const;

			/// The quadrature rule order used along each `Gamma_h` face by
			/// exteriorTransmissionRows(). Defaults to a rule generous enough that
			/// it is not what limits the row; raise it to check that.
			void setTransmissionQuadratureOrder( int order );

			/// **THIS IS THE HYPOTHESIS RADIAL TRANSFER PATHS REST ON, AND IT IS
			/// MEASURED RATHER THAN ASSUMED.** Free boundary puts `Gamma` on a
			/// semicircle centred on the axis, so the natural path family is rays
			/// from that centre -- `mfem::ClosestPointPath::Sphere`, whose foot map
			/// is exact, monotone and needs no search. What it needs in return is
			/// that every ray from the centre meets `Gamma_h` exactly once. If
			/// `Gamma_h` is star-shaped about the centre it does; if it is not, a ray
			/// meets it twice, two faces claim the same piece of `Gamma`, and the
			/// transfer is double valued there.
			///
			/// The condition is local and linear on a polygon: every boundary face
			/// must satisfy `( x - c ) . n > 0` for the outward normal `n`, which is
			/// visibility of that face from the centre. This returns the minimum of
			/// `( x - c ) . n / |x - c|` over the faces of `Gamma_h` -- a cosine, so
			/// it is dimensionless and comparable between meshes, exactly as
			/// `meq::AngleParametrisation`'s `min |u x t|` is for the contour tracer.
			///
			/// **AND IT MUST BE RE-CHECKED AFTER EVERY REFINEMENT.** A domain that
			/// starts star-shaped need not stay so: `meq::AdaptiveDomain` re-cuts
			/// `D_h` from the background mesh each cycle, and the elements it admits
			/// change. This is cheap -- one dot product per boundary face -- so the
			/// adaptive loop can afford to assert it every cycle, and should.
			///
			/// Returns the margin over `Gamma_h` alone on the extension path. On the
			/// fitted path there is no `Gamma_h`, so it measures the whole boundary
			/// and answers the same question about it.
			double starShapedMargin( double centreR, double centreZ ) const;

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

			/// q_h = ( 1/r ) grad_bar( psi ) in V_h, in MEQ's sign convention.
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
			 * refs/HDG-GradShafranov-Adaptive.pdf section 2.7 uses, and MEQ does
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
			 * which is why stage 3 was dropped as MEQ code, and why neither GS
			 * paper implemented it for accuracy. It is needed for the residual
			 * estimator of eq (20), whose eta_1, eta_2, eta_4 and eta_5 are all
			 * built on psi*_h rather than psi_h. That is FOUR of the five, and
			 * eta_4 is the one people miss: it is [[ psi*_h ]], not [[ psi_h ]].
			 * eta_2 differentiates the potential, and on psi_h that costs an
			 * order; measured in tests/convergence/EstimatorConvergence.cpp,
			 * which reports both.
			 *
			 * MEASURED, on the fitted Solov'ev benchmark: psi*_h converges at k+2
			 * for k = 1, 2, 3 -- 3.03, 4.03, 5.00 across the sequence -- so the
			 * library route delivers the superconvergence the paper wants and no
			 * hand-written local solve is needed. See EstimatorConvergence.cpp.
			 *
			 * IT USED TO BE REFUSED ON THE SEMI-LINEAR PATH, and the history is
			 * kept because the failure was silent and could return the same way.
			 * ReconstructFluxAndPot() read the LINEAR potential mass form and never
			 * looked at the non-linear one, so on MEQ's Newton path -- where the
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

			/// q*_h in [P_(k+1)]^2, in MEQ's sign convention -- the negation
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
			/// how MEQ gets the discrete one for free.
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
			///
			/// **When psi_ax is an unknown this is the AUGMENTED residual**,
			/// || ( R, gamma G ) ||, where G = psi_ax - max psi_h and gamma is
			/// || dR/dpsi_ax || frozen at the first iterate -- the factor that
			/// puts a perturbation of psi_ax into the units R is measured in. It
			/// is frozen rather than recomputed so that the history compares like
			/// with like; normalisationResidual() reports G on its own.
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

			/// Element-local NON-LINEAR iterations, summed over elements and over
			/// every residual and gradient evaluation since the forms were built.
			///
			/// **THIS IS THE ACCEPTANCE SIGNAL FOR NonlinearOrdering::NPC, and it
			/// is the only way to tell the two orderings apart from outside.**
			/// NPC linearises the full ( q, psi, psihat ) system and hybridizes
			/// the linear system that results, so every element-local operation
			/// is ONE linear solve against one factorisation and this reads
			/// EXACTLY ZERO. CondenseThenLinearise eliminates first, which makes
			/// each element's elimination its own non-linear solve, one per
			/// element per residual evaluation, and this reads in the thousands.
			///
			/// A non-zero count under NPC would mean the solve was not NPC --
			/// which is precisely the failure MFEM's deleted
			/// NLOrdering::LineariseThenCondense was, a condensation wearing the
			/// name. Zero on the linear path, where there is nothing to iterate.
			///
			/// @throws std::logic_error if the forms have not been built.
			long localNonlinearIterations() const;

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

			/// The same object as nonlinearSource when psi_ax is an unknown, and
			/// null otherwise. Held non-const because the solver writes the
			/// normalisation into it once per residual evaluation.
			NormalisedSource *normalisedSource;
			double psiAxisValue;
			double normalisationResidualValue;
			Normalisation normalisationChoice;
			BorderColumn borderColumnChoice;
			mfem::Coefficient *boundaryData;
			std::unique_ptr<mfem::Coefficient> potentialRhsCoeff;

			/// The Newton starting point, or null. Borrowed when it came in as a
			/// Coefficient; owned when setInitialGuess( GridFunction ) had to
			/// wrap one. Only one of the two is ever live.
			mfem::Coefficient *initialGuess;
			std::unique_ptr<mfem::Coefficient> ownedInitialGuess;

			/// The same guess as a FIELD, when it arrived as one, so that its
			/// gradient is available to seed the flux block. Null for the
			/// Coefficient overload, which has nothing to differentiate.
			mfem::GridFunction const *initialGuessField = nullptr;

			/// Solve `( r q_h, v ) = ( grad psi_g, v )` on each element and put
			/// `-q_h` -- DarcyForm's convention -- into the flux block.
			void seedFluxFromGuess( mfem::GridFunction const &psiGuess );

			/// Interpolate a coefficient onto the trace space, face by face.
			/// GridFunction::ProjectCoefficient cannot: it loops volume elements.
			/// See the .cpp.
			/// How nearly `Gamma_h` fails to be star-shaped about a centre, as a
			/// cosine in [ -1, 1 ]. Positive means it IS star-shaped; the value is
			/// the margin.
			///
			/// Project a PATH coefficient onto Gamma_h's trace dofs.
			///
			/// Separate from projectOntoTrace() because a path coefficient
			/// must be evaluated on the FACE transformation rather than the
			/// element one, and MFEM aborts rather than coping. See the
			/// definition.
			/// The potential dof nearest a point; FB-3's border pins psi_bnd
			/// to one nodal value, as psi_ax is pinned to the largest.
			int nearestPotentialDof( double r, double z ) const;

			void projectPathTraceOntoGammaH( mfem::Coefficient &coeff,
			                                 mfem::Vector &target ) const;

			void projectOntoTrace( mfem::Coefficient &coeff,
			                       mfem::GridFunction &target ) const;

			/// The transferring paths, or null on the fitted path. Borrowed.
			mfem::TransferPath *transferPath;
			int extensionLineOrder;

			/// `g` on Gamma, or empty. Wrapped in a PathTraceCoefficient and
			/// added to the FLUX right hand side in prepare(), which is where a
			/// non-homogeneous datum belongs; see the long note there.
			mfem::PositionFunction exteriorDatumFunction;
			std::unique_ptr<mfem::PathTraceCoefficient> exteriorDatumCoefficient;
			/// The flux equation's load, which carries the transferred exterior
			/// datum. A POINTER, and not a value, because prepare() is re-entrant
			/// on the coupled path: setExteriorCoupling() moves `a` and asks for
			/// the right hand side again, and mfem::LinearForm owns its
			/// integrators with no way to drop them. Re-adding one per prepare()
			/// while destroying the coefficient it references is a dangling read
			/// on the second Assemble() -- measured, as a segfault inside
			/// VectorBoundaryFluxLFIntegrator with nothing in the trace naming
			/// this file. Rebuilt whole instead.
			std::unique_ptr<mfem::LinearForm> fluxRhs;

			/// Face rule order for exteriorTransmissionRows(). SEPARATE from
			/// extensionLineOrder, which is a rule ALONG a path and has to match
			/// what buildForms() gave HDGExtensionIntegrator; this one is a rule
			/// ACROSS the face and is nobody else's business.
			int transmissionQuadratureOrder;

			/// FB-3's limiter contact, and whether one was given.
			bool boundaryFluxIsUnknown = false;
			bool currentIsUnknown = false;
			double targetMuZeroCurrent = 0.0;
			double currentScaleValue = 1.0;
			double plasmaCurrentValue = 0.0;
			double boundaryFluxR = 0.0;
			double boundaryFluxZ = 0.0;
			double psiBoundaryValue = 0.0;

			/// setAxisConstraint(), and where the last solve put the axis.
			AxisConstraint axisConstraintChoice = AxisConstraint::LocatedAxis;
			bool axisLocatedValue = false;
			double axisRValue = 0.0;
			double axisZValue = 0.0;

			Globalisation globalisationChoice;
			LocalSolver localSolverChoice;

			/// The exterior coupling of setExteriorCoupling(), borrowed, and the
			/// coefficients it solves for. The datum function reads the vector,
			/// so the transferred boundary condition follows the iterate.
			ExteriorDtN const *exteriorCoupling;
			std::vector<double> exteriorCoefficientValues;

			/// Extra quadrature order for meq::SourceIntegrator; see
			/// setSourceQuadratureOrder().
			int sourceQuadratureExtra;
			NonlinearOrdering orderingChoice;
			AssemblyMode assemblyModeChoice;
			TraceSolver traceSolverChoice;
			int andersonDepth;
			double picardDamping;

			/// The iterate the frozen source reads, on the Picard paths. Lives in
			/// the potential space, which is also the fixed point's unknown.
			std::unique_ptr<mfem::GridFunction> picardIterate;
#ifdef MFEM_USE_SUITESPARSE
			/// The Picard path's trace solver, held across iterations so that its
			/// retained symbolic analysis has something to be reused by. Built on
			/// first use; see picardStep().
			/// Hoisted out of picardStep() so the symbolic analysis survives
			/// between iterations. Typed as the base class since
			/// setTraceSolver() decides which one it is.
			std::unique_ptr<mfem::Solver> picardSolver;
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

			/// The bordered Newton of setSource( NormalisedSource &, double ):
			/// the trace and psi_ax solved together. See the .cpp.
			void solveWithNormalisation();

			/**
			 * The bordered Newton's COLUMN, `dR/ds`, assembled rather than
			 * differenced.
			 *
			 * `s` -- either normalisation -- reaches the residual ONLY through
			 * the source, so `dR/ds` is the assembly of `dF/ds` by exactly the
			 * loop meq::SourceIntegrator runs on `F`: same quadrature rule, same
			 * `-w F/r` sign, into the potential block and nowhere else. The flux
			 * and trace rows carry no `F` and are left at zero.
			 *
			 * @param axis  true for `dR/d(psi_ax)`, false for `dR/d(psi_bnd)`.
			 *
			 * @return false when it cannot be done -- no normalised source, not
			 *         NonlinearOrdering::NPC, or a source that does not supply
			 *         meq::NormalisedSource::normalisationDerivatives() -- and
			 *         then the caller differences as it always did. Under the
			 *         condensation the residual is the REDUCED trace one, which
			 *         is not this assembly at all.
			 */
			bool assembleNormalisationColumn( mfem::Vector const &state,
			                                  bool axis,
			                                  mfem::Vector &out ) const;

			/// `int F_plasma/r` over the domain at @a state, i.e. `mu0 I_p`.
			double assemblePlasmaCurrent( mfem::Vector const &state ) const;

			/// The current constraint's COLUMN, `dR/d(scale)`. `F` is linear in
			/// the scale, so this is the residual's own source term divided by
			/// it -- assembled by the same loop, with the same sign.
			void assembleCurrentColumn( mfem::Vector const &state,
			                            mfem::Vector &out ) const;

			/// Its ROW, `d( int F/r )/dx`: the plasma's own `dF/dpsi` integrated
			/// against the potential shape functions. A covector on the
			/// potential block and zero everywhere else.
			void assembleCurrentRow( mfem::Vector const &state,
			                         mfem::Vector &out ) const;

			/// AND ITS OFF-DIAGONAL CORNER ENTRIES, which are not zero and whose
			/// absence costs the quadratic rate rather than the answer.
			/// `int F/r` depends on `psi_ax` and `psi_bnd` EXPLICITLY, through
			/// the normalisation the profiles are evaluated at, so the current
			/// row of the corner block has entries against both of them.
			void assembleCurrentNormalisationCorner( mfem::Vector const &state,
			                                         double &againstAxis,
			                                         double &againstBoundary ) const;

			/// psi_h recovered from @a trace at normalisation @a psiAxisIn, and
			/// its largest nodal value -- which is the discrete psi_ax. Writes
			/// recoveryScratch and leaves the source's normalisation at
			/// @a psiAxisIn.
			///
			/// @param element  if not null, receives the element that attained it.
			/// @param dof      if not null, receives the potential dof that did.
			double recoverPeak( mfem::Vector const &trace, double psiAxisIn,
			                    double psiBoundaryIn,
			                    int *element = nullptr, int *dof = nullptr );

			/// The trace dofs of the faces of @a element: the only trace dofs the
			/// recovered potential on that element can depend on, and therefore
			/// the support of d( max psi_h )/dlambda. Measured rather than
			/// assumed -- see theAxisSensitivityIsLocalToItsElement.
			void traceDofsOfElement( int element, mfem::Array<int> &dofs ) const;

			/// Scratch for a trial recovery, so that a finite difference does not
			/// disturb the solution blocks the caller is going to read.
			mfem::BlockVector recoveryScratch;

			/// setPlasmaConnectivity(). Component by default: a connected
			/// plasma is what a confined source MEANS, and the pointwise test is
			/// kept as the control. Inert unless the source is confined.
			PlasmaConnectivity connectivityChoice = PlasmaConnectivity::Component;

			/// The fill. Unfilled -- so holds() is the constant true -- unless a
			/// confined source and PlasmaConnectivity::Component ask for it.
			PlasmaComponent plasmaComponentMask;

			/// The element adjacency, built once per mesh and reused: it is
			/// geometry, and refreshPlasmaComponent() runs once per residual.
			bool plasmaAdjacencyBuilt = false;

			/// Whether the confined source needs a fill at all.
			bool plasmaComponentWanted() const;

			/// The element adjacency of the solve mesh, as CSR, into
			/// plasmaComponentMask. Idempotent.
			void buildPlasmaAdjacency();

			/// Re-form the reduced system from whatever the solution blocks hold,
			/// which is how the element-local non-linear solves are given a fresh
			/// starting point: DarcyHybridization captures one at
			/// FormLinearSystem() time and keeps it. See solveWithNormalisation().
			void formSystem();

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
			/// MFEM_ASSERTs, which are dead in this build -- gated on
			/// MFEM_DEBUG, which the installed MFEM sets to NO, and NOT on
			/// NDEBUG, so building MEQ itself in Debug does not revive them.
			/// miniapps/hdg's
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
