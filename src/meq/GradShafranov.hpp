#ifndef MEQ_GRADSHAFRANOV_HPP
#define MEQ_GRADSHAFRANOV_HPP

#include <memory>

#include "mfem.hpp"

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
 * Stage 2 scope: the linear operator. F does not depend on psi, the domain is
 * polygonal and fitted so Gamma_h == Gamma, and there is no Newton iteration,
 * no local post-processing and no adaptivity. See CLAUDE.md for the stage plan.
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
 * block vector itself stays in DarcyForm's convention throughout, so that a
 * later Newton residual assembled by DarcyForm sees what it expects. Do not
 * "fix" one without the other.
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
	 * reason to prefer it once stage 4 turns the problem nonlinear.
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
	 * The hybridized HDG Grad-Shafranov solver.
	 *
	 * Owns the three finite element spaces, the DarcyForm built on them and the
	 * solution block vector. The mesh is borrowed and must outlive the solver;
	 * so must the source and boundary-data coefficients handed to setSource()
	 * and setBoundaryData().
	 *
	 * Usage:
	 *
	 *     meq::GradShafranovSolver solver( mesh, order );
	 *     solver.setSource( fCoefficient );          // F, not F/r
	 *     solver.setBoundaryData( psiCoefficient );  // psi on Gamma
	 *     solver.solve();
	 *     solver.potential();  solver.flux();
	 *
	 * Every boundary attribute of the mesh is Dirichlet: the fixed-boundary
	 * problem is an interior Dirichlet problem by construction, so there is no
	 * knob for that and no Neumann path to get wrong.
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

			/// The right hand side F( r, z ) of the equation as written above --
			/// NOT F/r. The 1/r belongs to the weak form and is applied here, which
			/// keeps meq::Source free of it too (see Source.hpp).
			void setSource( mfem::Coefficient &fIn );

			/// The Dirichlet datum g_D for psi on Gamma. Non-homogeneous data is
			/// the normal case: the level set psi = 0 is the plasma boundary, but a
			/// benchmark on a rectangle cut out of an exact equilibrium is not.
			void setBoundaryData( mfem::Coefficient &boundaryIn );

			/// Assemble and solve. Both coefficients must have been set.
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
			void assembleForms();

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

			/// -1/r, the factor the source is multiplied by to form the potential
			/// right hand side. Both halves of that are load bearing; see the .cpp.
			mfem::FunctionCoefficient negativeInverseRadius;

			mfem::Coefficient *source;
			mfem::Coefficient *boundaryData;
			std::unique_ptr<mfem::Coefficient> potentialRhsCoeff;

			/// Every boundary attribute, marked. See the class comment.
			mfem::Array<int> dirichletMarker;

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

			/// Aliases into solution. darcyFlux holds -q, see the file comment.
			mfem::GridFunction darcyFlux;
			mfem::GridFunction potentialGf;
			mfem::GridFunction traceGf;

			/// A copy of darcyFlux with the sign corrected, filled by solve().
			mfem::GridFunction fluxGf;

			bool assembled;
	};

}

#endif // MEQ_GRADSHAFRANOV_HPP
