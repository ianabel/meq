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

			/// The number of Newton iterations the last solve took. Zero on the
			/// linear path. One fewer than newtonResiduals().size(), since that
			/// counts the residual at the initial guess too.
			int newtonIterations() const;

			/// L2 errors against a closed form, on a quadrature rule generous
			/// enough that it does not itself limit the measured rate.
			double potentialError( mfem::Coefficient &exact ) const;
			double fluxError( mfem::VectorCoefficient &exact ) const;

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

			/// The transferring paths, or null on the fitted path. Borrowed.
			mfem::TransferPath *transferPath;
			int extensionLineOrder;

			double newtonRelativeTolerance;
			double newtonAbsoluteTolerance;
			int newtonMaxIterations;
			int newtonIterationCount;
			std::vector<double> newtonResidualHistory;

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

			/// The reduced trace system. traceX and traceB alias the trace blocks
			/// of solution and rhs; that aliasing is what makes FormLinearSystem()
			/// carry the essential trace values into the reduced problem, so it is
			/// not cosmetic. See the .cpp.
			mfem::OperatorHandle reduced;
			mfem::Vector traceX;
			mfem::Vector traceB;

			bool built;
			bool prepared;
	};

}

#endif // MEQ_GRADSHAFRANOV_HPP
