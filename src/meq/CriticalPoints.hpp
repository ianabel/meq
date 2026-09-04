#ifndef MEQ_CRITICALPOINTS_HPP
#define MEQ_CRITICALPOINTS_HPP

#include <vector>

#include "mfem.hpp"

#include "GradShafranov.hpp"

/*
 * Critical points of psi_h as objects: the magnetic axis, any X-point, and a
 * global audit of how many of each there can be. INVERSION-PLAN.md section 5
 * and stage IN-A.
 *
 * The magnetic axis and an X-point are where grad_bar(psi) vanishes. Every
 * contour method wants them located FIRST, because they are exactly where a
 * tracer's tangent is undefined and its corrector divides by zero -- so they
 * are not hard cases to be survived, they are objects to be found before the
 * tracing starts. This unit finds them and nothing else; IN-0 is the tracer.
 *
 * THE RESIDUAL IS A SOLVED FIELD, AND THAT IS THE WHOLE REASON THIS IS CHEAP
 * AND ACCURATE.
 *
 * MEQ's mixed formulation carries the flux q as an unknown of the same degree
 * as the potential, with r q = grad_bar(psi). So the equation to solve for a
 * critical point is
 *
 *     q_h( r, z ) = 0,
 *
 * a 2x2 system whose residual is a SOLVED variable converging at the
 * potential's own order, not a derivative of the potential converging one order
 * down. Compare CEDRES++, which records as an open problem that in P1
 * continuous Galerkin the axis and the X-point are confined to mesh vertices,
 * and TokaMaker, which notes for Lagrange order >= 2 that saddles "can exist
 * anywhere within the mesh". MEQ resolves them sub-element because q is a
 * polynomial inside each element and its zero set is found by root finding
 * rather than by looking at nodes.
 *
 * ONE LEVEL OF DIFFERENCING, NOT TWO -- AND THE ACCURACY OF THAT LEVEL DOES
 * NOT REACH THE ANSWER.
 *
 * Newton on q = 0 needs dq/dx, which is one differentiation of a degree-k L2
 * field and is therefore O(h^k) at best however it is obtained. It is taken
 * here by central differences of q in the element's own reference coordinates.
 * That looks like a shortcut and is not one, for a reason worth stating
 * plainly: the located root is where q_h vanishes, and where q_h vanishes is a
 * property of q_h alone. The Jacobian steers the iteration and appears nowhere
 * in its fixed point. A wrong Jacobian costs Newton steps and buys no error --
 * which is the same observation CLAUDE.md records for the Grad-Shafranov Newton
 * itself, from the other side: a Jacobian perturbed by 5% leaves every error
 * and every convergence rate unchanged to six figures.
 *
 * What the Jacobian IS load bearing for is the classification below, and there
 * only its two signs are used.
 *
 * THE HESSIAN OF PSI AT A ZERO OF q IS r TIMES THE JACOBIAN OF q, EXACTLY.
 *
 * Differentiating r q = grad_bar(psi) gives
 *
 *     Hess( psi ) = q (x) e_r + r dq/dx,
 *
 * and at a point where q = 0 the first term is identically zero. So
 * Hess( psi ) = r dq/dx there, with r > 0 throughout an axisymmetric domain.
 * The determinant scales by r^2 and the trace by r, both positive, so the SIGNS
 * that classify the point -- and therefore its Poincare-Hopf index -- can be
 * read off dq/dx without ever forming the Hessian. No second derivative of
 * psi_h is taken anywhere in this file.
 *
 * THIS IS NOT GradShafranovSolver::psiAxis(), AND THE TWO MUST NOT BE
 * RECONCILED.
 *
 * The solver's psi_ax is "the largest NODAL value of psi_h". It is deliberately
 * that and not the maximum of the polynomial, because the bordered Newton of
 * the normalised-profile path needs a constraint it can differentiate: one
 * nodal value is one entry of the discrete unknown, so the border row is sparse
 * -- exactly -e_j under NPC -- while the maximum of a polynomial over an
 * element is not a differentiable function of the coefficients at all. See
 * GradShafranovSolver::setSource( NormalisedSource &, double ) and CLAUDE.md's
 * "What works: psi_ax inside the residual".
 *
 * IN-A's axis is the critical point: the place where q_h vanishes. The two
 * quantities are different, both correct, and they differ by O(h) in POSITION
 * -- the distance from the polynomial's extremum to the nearest nodal point --
 * and by O(h^2) in VALUE, since psi_h is smooth and quadratic about its own
 * extremum. BOTH ORDERS ARE INDEPENDENT OF k, where psi_h's own error is k+1,
 * so the two readings separate rather than converge: measured on the finest
 * mesh of the Solov'ev benchmark the gap between them is 202 times psi_h's own
 * L2 error at k = 2 and 4204 times at k = 3.
 *
 * Neither is a defective version of the other and neither should be changed to
 * match. tests/convergence/CriticalPointConvergence.cpp measures the gap so
 * that the claim is a number rather than an assertion.
 *
 * THE AUDIT: A DEGREE IS A SUM OF INDICES AND NEVER A COUNT.
 *
 * audit() walks the boundary of the mesh and accumulates the turning of q,
 * which is the topological degree of q on that boundary. By the Poincare index
 * theorem that equals the SUM of the indices of the zeros of q inside -- +1 for
 * a maximum, +1 for a minimum, -1 for a saddle. It is a one-dimensional
 * integral: no subdivision, no root finding, and its cost is the boundary
 * rather than the domain.
 *
 * It is a CERTIFICATION, never an EXCLUSION. Degree zero does not imply no
 * root: a maximum and a saddle inside sum to zero and the boundary cannot tell
 * the difference from an empty domain. INVERSION-PLAN.md section 5 says this in
 * capitals because it was got wrong once during the survey that produced it,
 * and tests/convergence/CriticalPointConvergence.cpp keeps a domain containing
 * exactly one maximum and one saddle as a live demonstration: it reads a
 * winding number of zero with two critical points inside. Anything that uses
 * this class to decide "there is nothing here" is wrong.
 *
 * AND IT IS BLIND TO SPURIOUS PAIRS, WHICH IS WHY IT IS NOT THE ONLY CHECK
 * WANTED.
 *
 * Numerical noise in q_h creates critical points strictly in pairs -- a
 * spurious maximum next to a spurious saddle -- because a small perturbation of
 * a field cannot change its degree. So the pair sums to zero and the audit
 * passes with the pair present. The complementary test is a persistence
 * threshold, which does not need tuning: the stability theorem
 * (Cohen-Steiner, Edelsbrunner & Harer, 10.1007/s00454-006-1276-5) bounds
 * every spurious feature's persistence by 2 || psi_h - psi ||_inf, a quantity
 * MEQ already measures to convergence-rate precision. The index check and
 * persistence are COMPLEMENTARY, not redundant, and neither subsumes the
 * other. Persistence is not implemented here; this comment is where the reader
 * finds out that the audit alone does not cover it.
 *
 * WHEN THE DEGREE EQUALS THE EULER CHARACTERISTIC, AND WHEN IT MERELY HAPPENS
 * TO.
 *
 * Poincare-Hopf says the sum of the indices equals chi( Omega ) when the field
 * is TRANSVERSE to the boundary -- pointing outward (or inward) everywhere on
 * it and vanishing nowhere on it. That is a condition on q . n and NOT on the
 * boundary being a level set, which is worth separating because it is easy to
 * assume otherwise. A level set of psi with grad(psi) non-zero on it is one way
 * to get transversality and is what MEQ's own fixed-boundary Gamma gives; a
 * boundary that merely happens to lie outside every critical point, with the
 * flux pointing consistently outward across it, is another and is just as good.
 *
 * Both occur in the tests. On the standard benchmark rectangle
 * [0.6,1.4]x[-0.6,0.6] -- which is not a level set of anything -- q . n keeps
 * one sign the whole way round with min |q . n|/|q| = 0.15, so the hypothesis
 * holds and winding == chi == 1 is a theorem there. On a box drawn wide enough
 * to enclose an X-point it fails outright, measured at 0.00 with q . n changing
 * sign, and the degree reads 0 against chi = 1 -- no contradiction, because the
 * hypothesis is not satisfied.
 *
 * What survives in every case is degree == sum of the interior indices, which
 * needs no transversality at all. IndexAudit::transverse records which
 * situation the caller is in, so that "winding == chi" is not read as a theorem
 * where it is a coincidence, or as a defect where the hypothesis simply does
 * not hold.
 *
 * WHAT THIS IS NOT: AN EXHAUSTIVE SEARCH.
 *
 * findAxis() and sweep() are seeded Newton. Newton certifies the root it
 * converges to and says nothing whatever about the roots it does not.
 * INVERSION-PLAN.md section 5 specifies the exhaustive construction --
 * subdivision in the barycentric Bernstein basis with the convex-hull test of
 * Reuter et al., 10.1007/s00371-007-0184-x, where a sub-triangle all of whose
 * Bernstein coefficients share a sign provably contains no zero and can be
 * discarded -- and it is deliberately not built yet, because IN-A's acceptance
 * needs the axis and the audit and neither needs exhaustiveness. Do not read a
 * sweep() result as "these are all of them".
 *
 * SIGN CONVENTIONS, BOTH OF WHICH BITE DIFFERENTLY.
 *
 * GradShafranovSolver::flux() is +q, the sign flip out of DarcyForm's own
 * convention having already been undone. Handing this class the raw flux block
 * instead, which holds -q, leaves every winding number unchanged -- in even
 * dimension index( -v ) = index( v ) -- and silently swaps every Maximum for a
 * Minimum. That is the worse of the two failures, because the audit still
 * passes.
 *
 * And MEQ's psi is not sign-normalised across sources. With F single-signed
 * negative -- which is what the Solov'ev benchmarks have, F = -((1-A) r^2 + A)
 * -- psi is a subsolution, its maximum is on the boundary and the magnetic axis
 * is an interior MINIMUM. With F positive it is an interior maximum, which is
 * the case the high-beta source and INVERSION-PLAN.md section 6's maximum
 * principle argument are written for. So "seed from the largest nodal value" is
 * right for one sign of F and finds a corner of the mesh for the other.
 * findAxis() therefore seeds from BOTH nodal extremes and returns whichever
 * yields a genuine interior extremum, and refuses rather than guesses if both
 * do. AxisSense is there for a caller who knows which they want.
 */

namespace meq
{

	/// What the Hessian says a critical point is. Degenerate means the
	/// determinant is at round-off, where no classification is entitled -- it is
	/// reported rather than resolved, and carries index zero, which is not an
	/// index but an admission.
	enum class CriticalPointType
	{
		Maximum,
		Minimum,
		Saddle,
		Degenerate
	};

	/// "maximum", "minimum", "saddle", "degenerate". For printing.
	char const *criticalPointName( CriticalPointType type );

	/// The Poincare-Hopf index: +1 for either extremum, -1 for a saddle, 0 for
	/// a degenerate point. Note the first of those: a maximum and a minimum are
	/// NOT distinguished by their index in two dimensions, which is why
	/// INVERSION-PLAN.md section 6 needs the maximum principle and not just the
	/// topology to conclude that there is exactly one axis and no saddle.
	int criticalPointIndex( CriticalPointType type );

	/// One located zero of q_h.
	struct CriticalPoint
	{
		double r = 0.0;
		double z = 0.0;

		/// psi_h at the located point, from the same element.
		double psi = 0.0;

		CriticalPointType type = CriticalPointType::Degenerate;

		/// criticalPointIndex( type ), cached so that a sum over a vector of
		/// these does not have to re-derive it.
		int index = 0;

		/// The element whose polynomial was rooted. A zero lying near a face can
		/// be reached from either side, and the two answers differ by the
		/// O( h^(k+1) ) jump in q_h across it; sweep() merges them and keeps the
		/// one least outside its own element, so this names that one.
		int element = -1;

		/// | q_h | at the returned point, in the units of q. The Newton residual,
		/// kept because a root reported at 1e-3 is not a root.
		double fluxResidual = 0.0;

		/// det and trace of dq/dx at the point. Hess( psi ) is r times this, so
		/// the determinant of the Hessian is r^2 times determinant and its trace
		/// is r times trace -- both positive multiples, so the classification is
		/// the same either way. See the header comment.
		double determinant = 0.0;
		double trace = 0.0;

		/// How far outside its element the root lies, in reference-element
		/// units; zero when it is strictly inside, which is the ordinary case.
		///
		/// A NON-ZERO VALUE HERE IS NOT AN ERROR AND IS WORTH REPORTING RATHER
		/// THAN HIDING. q_h is discontinuous, so a zero lying within the jump of
		/// a face belongs to NEITHER of the two elements strictly: each side's
		/// polynomial puts its own zero a little way into the other's territory,
		/// and refusing both would make the search fail whenever the axis
		/// happens to land on a mesh line. Measured on the Solov'ev benchmark at
		/// k = 1, n = 4 -- where the axis at z = 0.0111 sits beside the mesh line
		/// z = 0 -- the two candidates are 6.6e-4 and 8.9e-2 outside their
		/// elements, and with no allowance at all the axis is not found. It is
		/// the ONLY point of the twelve in
		/// tests/convergence/CriticalPointConvergence.cpp's k = 1, 2, 3 by
		/// n = 4, 8, 16, 32 sweep that needs one. See setContainment().
		double overshoot = 0.0;
	};

	/// Which extremum findAxis() should accept.
	enum class AxisSense
	{
		/// Take whichever of the two is found, and throw if both are. The
		/// default, because MEQ's psi is not sign-normalised and a caller
		/// usually does not want to have to know which way round F points.
		Either,
		Maximum,
		Minimum
	};

	/// The result of the boundary audit. Everything in it is a measurement of
	/// one walk around the boundary of the mesh.
	struct IndexAudit
	{
		/// The topological degree of q on the boundary, that is the total
		/// turning divided by 2 pi and rounded. Equals the sum of the indices of
		/// the interior zeros. NEVER a count of them: see the header.
		int windingNumber = 0;

		/// The unrounded total turning, in units of 2 pi.
		double turning = 0.0;

		/// | turning - windingNumber |. A degree is an integer, so this is a
		/// direct measurement of whether the boundary was sampled finely enough:
		/// it is not a discretisation error that shrinks with h, it is either at
		/// round-off or the walk missed a rotation.
		double windingDefect = 0.0;

		/// chi( Omega_h ) = V - E + F of the mesh. 1 for a disc, 0 for an
		/// annulus. On a non-conforming mesh the vertex and edge counts include
		/// hanging entities, and this is then a count of the refined mesh rather
		/// than of the domain -- which is still the right answer for a
		/// conforming refinement of a disc, and is not checked here.
		int eulerCharacteristic = 0;

		/// How many closed loops the boundary faces form. For a connected planar
		/// domain chi = 2 - boundaryLoops, which is an independent route to the
		/// line above and disagrees with it if the mesh is not what it is taken
		/// to be.
		int boundaryLoops = 0;

		/// The largest single angular increment in the walk, in radians. A value
		/// approaching pi means the walk is under-sampled and the winding number
		/// is not to be believed; setBoundarySamples() is the remedy.
		double worstTurn = 0.0;

		/// The smallest | q | seen on the boundary, in the units of q. A zero of
		/// q ON the boundary makes the degree undefined, and this is how close
		/// the walk came to one.
		double smallestFlux = 0.0;

		/// min | q . n | / | q | over the walk, with n the outward normal, and
		/// whether q . n kept one sign throughout. Together they say whether q is
		/// transverse to the boundary and therefore whether comparing
		/// windingNumber against eulerCharacteristic is Poincare-Hopf or a
		/// coincidence. See the header.
		double transversality = 0.0;
		bool transverse = false;

		/// windingNumber == eulerCharacteristic. Meaningful as a theorem only
		/// when transverse; otherwise it is a statement about what is inside.
		bool consistent() const
		{
			return windingNumber == eulerCharacteristic;
		}
	};

	/// V - E + F for a 2D mesh. Free function because it is a property of the
	/// mesh and nothing to do with q.
	int eulerCharacteristic( mfem::Mesh &mesh );

	/**
	 * Zeros of q_h, and the boundary audit over them.
	 *
	 * Borrows the flux and the potential; both must outlive it, and the solver
	 * they came from must have been solved. Nothing is cached: every call walks
	 * the mesh again, which costs at worst two Newtons per element -- one from
	 * the element centre and one from its quietest flux node -- and is
	 * negligible beside the solve that produced the field.
	 *
	 * The potential is used for two things only -- seeding findAxis() from its
	 * extreme nodal values, and reporting psi at a located point. The roots
	 * themselves are a property of the flux alone.
	 */
	class CriticalPointFinder
	{
		public:
			/// The ordinary way in. Takes solver.flux() and solver.potential().
			explicit CriticalPointFinder( GradShafranovSolver const &solverIn );

			/// The same over bare fields, so that the finder can be pointed at an
			/// interpolated exact q -- which is how a test separates "the root
			/// finder is wrong" from "the discretisation is coarse".
			///
			/// @param fluxIn      q in MEQ's sign convention, vdim 2. See the
			///                    header on what handing it -q does.
			/// @param potentialIn psi_h, on the same mesh.
			CriticalPointFinder( mfem::GridFunction const &fluxIn,
			                     mfem::GridFunction const &potentialIn );

			/**
			 * The magnetic axis: the interior extremum of psi_h, as a zero of q_h.
			 *
			 * Seeded from the elements holding the largest and the smallest nodal
			 * values of psi_h and from two rings of face neighbours around each,
			 * then Newton on q_h = 0 in each seed element using that element's own
			 * polynomial. A root is accepted only if it lies in the element whose
			 * polynomial produced it, up to setContainment()'s allowance.
			 *
			 * THAT IS THE FAST PATH AND NOT THE ONLY ONE. Where it does not
			 * produce exactly one extremum strictly inside an element, the search
			 * falls back to sweep() and takes its answer -- which costs one Newton
			 * per element and is negligible beside the solve that produced the
			 * field. Both halves are measured in the implementation; the seed is a
			 * heuristic and the answer is not allowed to depend on it.
			 *
			 * @throws std::runtime_error if no interior extremum is found, or if
			 *         @a sense is Either and both a maximum and a minimum are.
			 *         The second is a genuine ambiguity rather than a failure --
			 *         INVERSION-PLAN.md section 6 argues it cannot happen for the
			 *         fixed-boundary problem with single-signed F -- and guessing
			 *         would be worse than refusing.
			 */
			CriticalPoint findAxis( AxisSense sense = AxisSense::Either ) const;

			/// The same without the throw. Returns false where findAxis() would
			/// throw, and leaves @a found untouched.
			bool tryFindAxis( CriticalPoint &found,
			                  AxisSense sense = AxisSense::Either ) const;

			/**
			 * Newton from every element, deduplicated: every zero of q_h that a
			 * seeded search happens to reach.
			 *
			 * NOT EXHAUSTIVE, and the header says why at length. Its value is that
			 * the sum of the indices it returns can be compared against audit()'s
			 * winding number, and a disagreement is then positive evidence that
			 * something was missed -- which is a use of the degree as a
			 * certification, the only use it has.
			 */
			std::vector<CriticalPoint> sweep() const;

			/// Walk the boundary and accumulate the turning of q.
			IndexAudit audit() const;

			/// Newton stops when | q | falls below this times the largest | q | on
			/// the mesh. Default 1e-13: q_h is a polynomial and Newton on it is
			/// quadratic, so this is reached in a handful of steps or not at all.
			void setTolerance( double toleranceIn );

			/// Default 50. A cap, not a target.
			void setMaxIterations( int maxIterationsIn );

			/// The central-difference step for dq/dx, in reference-element units.
			/// Default 1e-4, which puts the truncation error at 1e-8 relative and
			/// the round-off at about the same -- and neither reaches the answer,
			/// per the header.
			void setJacobianStep( double stepIn );

			/// Points sampled per boundary face by audit(). Default 16. Raise it
			/// if IndexAudit::worstTurn approaches pi.
			void setBoundarySamples( int samplesIn );

			/// Two roots closer than this, relative to the diameter of the mesh,
			/// are one root. Default 1e-8. Only sweep() uses it.
			void setSeparation( double separationIn );

			/// How far outside its own reference element a root is still
			/// attributed to that element, in reference-element units. Default
			/// 0.10, against a worst measured need of 4.0e-2 on the benchmarks in
			/// tests/convergence/CriticalPointConvergence.cpp.
			///
			/// It exists because q_h is discontinuous and CriticalPoint::overshoot
			/// says what it costs. The units are the point: a FIXED tolerance in
			/// reference space is a shrinking one in physical space, so the
			/// allowance vanishes with the mesh exactly as the ambiguity it covers
			/// does. Set it to zero to refuse anything not strictly inside, which
			/// is the pedantically correct behaviour and which fails on coarse
			/// meshes whose axis lands on a mesh line.
			///
			/// IT IS NOT A TUNING PARAMETER, and that was checked rather than
			/// hoped. Swept over 0.001, 0.01, 0.05, 0.10 and 0.20 on the whole
			/// k = 1, 2, 3 by n = 4, 8, 16, 32 benchmark, the located axis is
			/// IDENTICAL at every point to every digit printed. Two rules are what
			/// make it so: a candidate outside its element never stops the search
			/// early -- findAxis() spends a full sweep whenever the seeded path
			/// returns one -- and where two candidates are merged the one least
			/// outside its own element wins, rather than the one with the lower
			/// element index.
			void setContainment( double containmentIn );

		private:
			/// Newton on q_h = 0 inside one element, from one reference-space
			/// seed. Returns false unless it converged to a point inside that
			/// element.
			///
			/// @a target is the absolute residual to stop at, which the callers
			/// compute once from fluxScale(). It is a parameter rather than a
			/// lookup because fluxScale() costs a pass over every flux dof, and
			/// calling it once per element would make a sweep quadratic in the
			/// mesh for no gain whatever -- the scale is a property of the field,
			/// not of the element.
			bool rootInElement( int element, mfem::IntegrationPoint const &seed,
			                    double target, CriticalPoint &found ) const;

			/// dq/dxi by central differences, with the stencil slid so that both
			/// samples stay inside the reference element.
			void referenceJacobian( int element, mfem::IntegrationPoint const &ip,
			                        double jacobian[ 2 ][ 2 ] ) const;

			/// The largest | q | over the flux dofs, as the scale the Newton
			/// tolerance is relative to.
			double fluxScale() const;

			/// Element indices to seed findAxis() from: the elements holding the
			/// extreme nodal values of psi_h, and two rings of face neighbours
			/// around each.
			///
			/// TWO RINGS RATHER THAN ONE, AND THAT IS A MEASUREMENT. The extreme
			/// nodal value can be further from the critical point than it looks:
			/// where the point sits on a mesh line the extreme node is the shared
			/// vertex, and which element is credited with it is decided by the L2
			/// jump at that vertex -- a difference of 1e-8 picking between
			/// elements one apart. Measured on iterExample2 at k = 2, n = 24, the
			/// minimum nodal value is in element 694 and the root is in element
			/// 696, which is not a face neighbour of it.
			void axisSeeds( std::vector<int> &elements ) const;

			/// The extrema among the roots reachable from a list of seed
			/// elements, merged and filtered by @a sense.
			std::vector<CriticalPoint> extremaFrom( std::vector<int> const &elements,
			                                        AxisSense sense ) const;

			/// The seeds this element offers: its centre, and the node where
			/// | q_h | is smallest.
			void elementSeeds( int element,
			                   std::vector<mfem::IntegrationPoint> &seeds ) const;

			mfem::GridFunction const &fluxField;
			mfem::GridFunction const &potentialField;
			mfem::Mesh &meshRef;

			double tolerance = 1.0e-13;
			int maxIterations = 50;
			double jacobianStep = 1.0e-4;
			int boundarySamples = 16;
			double separation = 1.0e-8;
			double containment = 0.10;
	};

}

#endif
