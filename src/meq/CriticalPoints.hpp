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

		/// WHERE THE ROOT SITS IN @a element's REFERENCE COORDINATES, which is
		/// what a caller needs to evaluate that element's shape functions there.
		///
		/// IT IS PLUMBED OUT RATHER THAN RECOVERABLE, AND THAT IS THE POINT. The
		/// Newton runs in reference space, so this is the iterate it converged
		/// to and costs nothing to report. The alternative -- handing back only
		/// ( r, z ) and letting the caller invert the element map with
		/// TransformBack -- re-solves a problem that was already solved, and
		/// CLAUDE.md records the failure that invites: a CLAMPED inverse map
		/// returns a point on the element boundary instead of failing, and a
		/// field that is constant over the element cannot tell the difference.
		/// The one consumer that would notice is exactly the one this exists
		/// for, FREE-BOUNDARY-PLAN.md section 11.5's option 3, where the shape
		/// functions at this point ARE a row of the bordered Jacobian.
		///
		/// TWO DOUBLES RATHER THAN AN mfem::IntegrationPoint, because that class
		/// has no default member initialiser -- a defaulted CriticalPoint would
		/// carry an uninitialised one, and the rest of this struct is
		/// zero-initialised. referencePoint() builds the IntegrationPoint a
		/// caller wants, so the convenience is kept without the hazard.
		double referenceX = 0.0;
		double referenceY = 0.0;

		/// ( referenceX, referenceY ) as MFEM wants it, for CalcShape(),
		/// GetValue() and the transformations. The weight is meaningless here
		/// and is left at whatever Set2() leaves it: this is a POSITION, not a
		/// quadrature point.
		mfem::IntegrationPoint referencePoint() const
		{
			mfem::IntegrationPoint point;
			point.Set2( referenceX, referenceY );
			return point;
		}
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

	/**
	 * Whether the psi_ax a solve REPORTS is the flux at a magnetic axis.
	 *
	 * THE DEFINITION SAYS NOTHING ABOUT AN AXIS, AND THAT IS THE WHOLE OF THE
	 * PROBLEM THIS ANSWERS. GradShafranovSolver::psiAxis() is the largest NODAL
	 * value of psi_h, chosen for the reason the header gives at length: one
	 * nodal value is one entry of the discrete unknown, so the bordered Newton's
	 * row is exactly -e_j. Nothing in "the largest number in the potential
	 * vector" says that number sits at a magnetic axis, and the bordered Newton
	 * cannot notice: G( lambda, s ) = s - max psi_h is satisfied at machine zero
	 * by a spurious nodal spike exactly as it is by an axis.
	 *
	 * MEASURED, and this is why the check exists rather than being a precaution.
	 * On a free-boundary machine case at k = 2 -- the SAME mesh, source and guess
	 * on which k = 3 reproduces its reference to four digits -- the solve
	 * converged in 17 Newton steps with psi_ax - max psi_h reading 0.000e+00, the
	 * prescribed plasma current delivered to seven figures and every border at
	 * machine zero, and reported psi_ax = 2.734289e+00 against a peak of
	 * 8.64e-02: twenty-nine times too large, and a SINGLE dof of one element.
	 * The runaway is self consistent -- a spurious psi_ax inflates the span, Psi
	 * collapses, and the current border raises the profile scale by 980 to hold
	 * the current -- so every diagnostic the run prints is green.
	 *
	 * WHAT IS COMPARED, AND WHY IT IS THE NORMALISED FLUX RATHER THAN A RATIO OF
	 * FLUXES. The profiles are functions of
	 *
	 *     Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ),
	 *
	 * and Psi at the magnetic axis is 1 BY DEFINITION when psi_ax is the axis
	 * flux. So evaluating Psi at a genuine O-point -- a zero of q_h, found by
	 * this class rather than read off the nodes -- is a direct statement about
	 * the quantity the source actually consumes, and it reduces to psi_O/psi_ax
	 * on every fixed-boundary run, where psi_bnd is zero.
	 *
	 * IT IS ONE SIDED, AND THE SIDE IS FORCED. psi_ax is the largest nodal value,
	 * and the peak of a polynomial over a closed element is at least its largest
	 * nodal value, so on a healthy field Psi at the axis is 1 from ABOVE, short
	 * only of the O( h^(k+1) ) by which q_h disagrees with itself across the face
	 * that carries the extreme node. The failure drives it DOWN -- 0.032 on the
	 * machine case above, 0.31 and 0.36 on a half-disc that latched onto the
	 * corner where Gamma meets the axis. So the guard is on the low side and the
	 * high side is reported rather than tested: a Psi above 1 is a field whose
	 * element interior peaks well above its own nodes, which is coarseness rather
	 * than a wrong equilibrium.
	 *
	 * WHAT IT CANNOT DO. sweep() is seeded Newton and is not exhaustive, so a
	 * spurious extremum it happens to reach with a HIGHER Psi than the true axis
	 * makes this pass -- the largest Psi is taken deliberately, so that noise
	 * costs a missed detection rather than a false alarm. extrema and saddles are
	 * reported for that reason: an axis found among nine other maxima is not the
	 * same evidence as an axis found alone. And separation is corroboration and
	 * not a second trigger, because on a graded mesh the axis element's own
	 * diameter is not the scale the spike lives on.
	 */
	struct AxisAgreement
	{
		/// What the solve reported, and the boundary flux it is measured against.
		double psiAxis = 0.0;
		double psiBoundary = 0.0;

		/// The extreme nodal value of psi_h recomputed here, and where it sits.
		/// It is the same argmax GradShafranovSolver takes -- same field, same
		/// rule -- so a disagreement with psiAxis above is itself a finding, and
		/// is the reason this is recomputed rather than only passed in.
		double nodalExtreme = 0.0;
		double nodeR = 0.0;
		double nodeZ = 0.0;
		int nodeElement = -1;

		/// Whether an extremum of the sense the span asks for was reached at all.
		/// FALSE IS A RESULT, not an absence of one: a normalised solve whose flux
		/// has no interior extremum anywhere is the wall-hugging annulus branch,
		/// where psi rises monotonically to Gamma and there is no closed surface
		/// to be an axis of.
		bool located = false;

		/// The O-point carrying the largest normalised flux. Valid if located.
		CriticalPoint axis;

		/// How many extrema of that sense the sweep reached, and how many saddles
		/// beside them. See the header on why noise arrives in pairs.
		int extrema = 0;
		int saddles = 0;

		/// Psi at @a axis. IT MUST BE 1.
		double normalisedFlux = 0.0;

		/// | x_axis - x_node |, in metres and in diameters of the axis element.
		/// Corroboration for a reader, never the trigger.
		double separation = 0.0;
		double separationInElements = 0.0;

		/// located, and normalisedFlux is within the tolerance below 1.
		bool agrees = false;
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
			 * The axis NEAR a point already believed to be close to it: Newton on
			 * `q_h = 0` seeded from the element nearest @a r, @a z and a couple of
			 * rings of face neighbours around it, and nothing else.
			 *
			 * THIS IS THE WARM-START ENTRY POINT AND ITS WHOLE PURPOSE IS COST.
			 * findAxis() and tryFindAxis() are written for a caller with no prior:
			 * they seed from the extreme NODAL values -- the quantity
			 * FREE-BOUNDARY-PLAN.md section 11 is about -- and they fall back to a
			 * full sweep() whenever the seeded path is not unambiguous, which is
			 * the right trade when the answer is wanted once after a solve. It is
			 * the wrong trade INSIDE a Newton loop, where the axis is wanted once
			 * per Jacobian and the previous iterate's axis is a seed that is
			 * already correct to the size of the last step. This pays a search
			 * over a handful of elements instead of two Newtons over every one.
			 *
			 * WHAT IT COSTS, MEASURED IN ELEMENTS ROOTED RATHER THAN IN SECONDS
			 * -- a timing on one machine is a measurement about that machine, and
			 * this is a property of the algorithm. Solov'ev at k = 2, seeded a
			 * given fraction of an element from the answer:
			 *
			 *     offset      n = 16, of 512      n = 32, of 2048
			 *     0            1   ( 0.20% )       1   ( 0.05% )
			 *     0.5 h       14   ( 2.73% )      19   ( 0.93% )
			 *     1.0 h       14   ( 2.73% )      19   ( 0.93% )
			 *     2.0 h       37   ( 7.23% )      48   ( 2.34% )
			 *
			 * against sweep(), which roots every element. So a seed that is still
			 * on the answer costs ONE element, and the ratio IMPROVES with
			 * refinement -- the search is bounded by rings, which is an absolute
			 * element count, while a sweep is the whole mesh. That is the property
			 * a Newton loop needs: the per-Jacobian cost does not grow with the
			 * problem.
			 *
			 * THERE IS NO SWEEP FALLBACK, DELIBERATELY. Adding one would restore
			 * exactly the cost this exists to avoid, and would do it on the
			 * iterations where the seed is worst -- which in a continuation are
			 * the ones where the answer matters least, because a later iterate
			 * will correct it. A caller who needs the guaranteed answer should
			 * call findAxis(), and the natural pattern is findAxis() once to start
			 * and this thereafter.
			 *
			 * AND IT BREAKS TIES BY DISTANCE, WHERE tryFindAxis() REFUSES THEM.
			 * That is the other difference and it follows from the same premise:
			 * a caller with a prior is FOLLOWING one critical point, so when the
			 * seed region offers more than one extremum of the requested sense the
			 * nearest to ( @a r, @a z ) is the continuation of the one being
			 * followed. tryFindAxis() has no prior and so cannot prefer one, and
			 * refuses instead. Do not use this to DISCOVER an axis: seeded far
			 * from one it will return whatever extremum happens to lie in reach,
			 * which is a different question from "where is the axis".
			 *
			 * @param r,z   where to start looking. Need not be inside the mesh and
			 *              need not be near an element boundary; the nearest
			 *              element CENTRE is what is used, which costs one pass
			 *              over the elements with no field evaluation in it.
			 * @param sense which extremum is wanted. The same caution applies as
			 *              everywhere else in this class: pass +q, never the raw
			 *              flux block, or every Maximum silently becomes a
			 *              Minimum. See the file header.
			 * @param found written only on success.
			 *
			 * @return false, rather than throwing, when no extremum of that sense
			 *         is reachable from the seed region. There is nothing
			 *         exceptional about that -- an iterate whose axis has left the
			 *         search radius is an ordinary event in a continuation -- and
			 *         the caller is expected to have a fallback.
			 */
			bool tryFindAxisFrom( double r, double z, AxisSense sense,
			                      CriticalPoint &found ) const;

			/// The MOST rings of face neighbours tryFindAxisFrom() will grow
			/// around its seed element before giving up. Default 6.
			///
			/// IT IS A CAP AND NOT A COST, because the search stops at the first
			/// ring that yields an extremum of the requested sense. A warm start
			/// whose axis is still in the seed element roots ONE element; only a
			/// seed that has fallen behind pays for the outer rings. That is why
			/// the default can be generous.
			///
			/// AND SIX RATHER THAN axisSeeds()' TWO, WHICH IS MEASURED. A ring is
			/// a hop across a face, not a distance: on a diagonally split
			/// Cartesian mesh each hop advances about a quarter of an element, so
			/// two rings reach barely half of one. Measured on the Solov'ev
			/// benchmark at k = 1, n = 16, a seed half an element away diagonally
			/// needs THREE rings and two finds nothing. axisSeeds() gets away with
			/// two because it grows them from the element holding the extreme
			/// NODAL value, which is already within a node of the answer; a
			/// general seed is not.
			///
			/// It is a RADIUS rather than a tolerance: raising it lets the search
			/// follow an axis that has moved further between calls, and cannot
			/// change the answer where the axis has not moved.
			void setSeedRings( int ringsIn );

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

			/**
			 * Is @a psiAxisIn the flux at a magnetic axis, or merely the largest
			 * number in the potential vector? See AxisAgreement, which is where
			 * the reasoning and the measurements are.
			 *
			 * @param psiAxisIn     what the solve reported, ordinarily
			 *                      GradShafranovSolver::psiAxis().
			 * @param psiBoundaryIn psi_bnd, ordinarily
			 *                      GradShafranovSolver::psiBoundary(). Zero on
			 *                      every fixed-boundary run, where the profiles
			 *                      take Psi = psi/psi_ax.
			 * @param toleranceIn   how far below 1 the normalised flux at the
			 *                      located axis may fall before AxisAgreement
			 *                      says the two disagree. Default 0.10, which is
			 *                      generous rather than tuned: a healthy field
			 *                      reads 1 from above, and the failures this
			 *                      exists for read 0.03 and 0.31.
			 *
			 * The sense follows the SIGN OF THE SPAN and is not guessed: the
			 * plasma is where ( psi - psi_bnd ) has the span's sign, so a
			 * positive span puts the axis at a maximum and a negative one at a
			 * minimum. That is meq::NormalisedSource::insidePlasma()'s own test,
			 * and it is why this does not meet AxisSense::Either's ambiguity.
			 *
			 * @throws std::invalid_argument if the span is zero, where Psi is
			 *         undefined and there is nothing to compare.
			 */
			AxisAgreement checkAxis( double psiAxisIn, double psiBoundaryIn = 0.0,
			                         double toleranceIn = 0.10 ) const;

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

			/// The extreme nodal value of the potential and where it sits: the
			/// same search GradShafranovSolver runs to produce psi_ax, so that
			/// checkAxis() compares against the point psi_ax is actually AT
			/// rather than against one recovered some other way.
			///
			/// The position is the node's own, from the element's nodal
			/// IntegrationRule. A basis whose node count does not match its dof
			/// count -- which no space MEQ builds has -- falls back to the
			/// element centre, so that a diagnostic never throws over a basis
			/// choice.
			void nodalExtreme( bool wantMaximum, double &value, int &element,
			                   double &r, double &z ) const;

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

			/// The element whose CENTRE is nearest ( r, z ), or -1 on an empty
			/// mesh.
			///
			/// NEAREST CENTRE AND NOT THE CONTAINING ELEMENT, WHICH IS A COST
			/// DECISION AND NOT AN APPROXIMATION THAT COULD BE WRONG. Locating
			/// the containing element properly means inverting element maps --
			/// Mesh::FindPoints, which CLAUDE.md records as O( elements x points )
			/// and whose shared scratch makes it non-reentrant besides. This is
			/// one distance per element and no field evaluation at all, which is
			/// about a per cent of what a sweep costs, and it feeds a search that
			/// grows by face neighbours anyway: a seed one element out is
			/// absorbed by the first ring.
			int nearestElementCentre( double r, double z ) const;

			mfem::GridFunction const &fluxField;
			mfem::GridFunction const &potentialField;
			mfem::Mesh &meshRef;

			double tolerance = 1.0e-13;
			int maxIterations = 50;
			double jacobianStep = 1.0e-4;
			int boundarySamples = 16;
			double separation = 1.0e-8;
			double containment = 0.10;
			int seedRings = 6;
	};

}

#endif
