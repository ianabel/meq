#ifndef MEQ_GRADSHAFRANOV_HPP
#define MEQ_GRADSHAFRANOV_HPP

#include <memory>
#include <vector>

#include "mfem.hpp"

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
				/// Thread the element-local half of ComputeH() -- the factorisation
				/// of A, the Schur complement and its factorisation, and one local
				/// back-substitution per trace dof. The scatter into the trace
				/// matrix is NOT threaded and cannot be: an unfinalized
				/// mfem::SparseMatrix carries one current_row and one column-pointer
				/// scratch for the whole matrix, so two threads writing rows that
				/// are disjoint by construction still collide, and the failure is a
				/// hang rather than a wrong answer. That serial scatter is the
				/// Amdahl ceiling on this option.
				///
				/// **Requires MFEM_USE_OPENMP and MFEM_THREAD_SAFE, and MFEM
				/// ABORTS rather than falling back if the build lacks either** --
				/// deliberately, since a caller asking for this is asking a
				/// performance question and a silent serial loop is not an answer
				/// to it. MEQ therefore checks the build before passing it on.
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
			/// Throws std::invalid_argument rather than letting MFEM abort the
			/// process when the library was built without OpenMP or without
			/// thread safety.
			void setAssemblyMode( AssemblyMode choice );

			/// The mode buildForms() will ask for.
			AssemblyMode assemblyMode() const;

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

			/// psi_h recovered from @a trace at normalisation @a psiAxisIn, and
			/// its largest nodal value -- which is the discrete psi_ax. Writes
			/// recoveryScratch and leaves the source's normalisation at
			/// @a psiAxisIn.
			///
			/// @param element  if not null, receives the element that attained it.
			/// @param dof      if not null, receives the potential dof that did.
			double recoverPeak( mfem::Vector const &trace, double psiAxisIn,
			                    int *element = nullptr, int *dof = nullptr );

			/// The trace dofs of the faces of @a element: the only trace dofs the
			/// recovered potential on that element can depend on, and therefore
			/// the support of d( max psi_h )/dlambda. Measured rather than
			/// assumed -- see theAxisSensitivityIsLocalToItsElement.
			void traceDofsOfElement( int element, mfem::Array<int> &dofs ) const;

			/// Scratch for a trial recovery, so that a finite difference does not
			/// disturb the solution blocks the caller is going to read.
			mfem::BlockVector recoveryScratch;

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
