#ifndef MEQ_SURFACEAVERAGE_HPP
#define MEQ_SURFACEAVERAGE_HPP

#include <cstddef>
#include <functional>
#include <vector>

#include "CriticalPoints.hpp"
#include "FluxSurfaces.hpp"

/*
 * Flux-surface averages: INVERSION-PLAN.md stage IN-2, on the fitted path.
 *
 * THE CONVENTION IS PART OF THE DEFINITION, so it is first and it is stated
 * once:
 *
 *     V'( psi ) = closed-integral 2 pi R dl / | grad psi |
 *     < X >_psi = ( 1 / V' ) closed-integral 2 pi R X dl / | grad psi |
 *
 * The 2 pi R is IN. A per-unit-length V' differs by exactly that factor, and
 * the flux-surface-averaged Grad-Shafranov identity below -- which is this
 * file's reference-free acceptance -- is stated for this one and not for the
 * other. Anything reading a V' out of meq is reading the volume derivative,
 * dV/dpsi with V the volume enclosed by the surface.
 *
 * ONE FACILITY, NOT A LIST OF QUANTITIES, AND THAT IS A RESPONSE TO THE
 * CONSUMER RATHER THAN A PREFERENCE.
 *
 * MANTA-COUPLING.md section 5 says the geometry slots are "a set of scalar
 * functions of the flux label, each differentiable with respect to every psi
 * DOF", and says of the concrete list that it "is negotiated with the transport
 * physics case, not fixed by MaNTA" -- and that the illustrative tokamak set it
 * prints is known to be not entirely correct and is being revised. A list that
 * is still moving must not be baked in as a list of functions.
 *
 * So the primitive is: given a surface, take a callable integrand in
 * ( R, z, psi, q ) and return < f > and V'. < R^-2 >,
 * < | grad psi |^2 / R^2 >, the arc length and the safety factor are then one
 * line each. Adding or removing a slot is a line; the convergence machinery is
 * written once; and the dGeometry_dpsi derivative MANTA-COUPLING.md section 6
 * needs is taken of the facility rather than of eight separate expressions.
 * tests/analytic/FluxSurfaceReference.hpp is the same shape on the EXACT field
 * and is the prototype this is the discrete counterpart of.
 *
 * A NAME COLLISION THAT IS RESPECTED THROUGHOUT. In this project q is the
 * FLUX, q = ( 1/r ) grad_bar( psi ), a solved unknown of the discretisation and
 * the asset the whole inversion item is built on. The safety factor is also
 * universally written q. IN THIS FILE AND IN src/meq THE SAFETY FACTOR IS
 * safetyFactor AND NEVER q: a reader who meets q here is entitled to assume the
 * flux, and one silent conflation would be very hard to see afterwards.
 *
 * | grad psi | COMES FROM q, POINTWISE, AND THAT IS THE POINT OF THE WHOLE
 * ITEM.
 *
 * grad_bar( psi ) = r q, so | grad psi | = r | q | with q read at the node from
 * the solved flux field -- converging at the potential's own order rather than
 * one order down, which is what a differentiated L2 potential would give. Every
 * weight in this file is 2 pi R | dx/ds | ds / ( r | q | ) for whichever
 * parametrisation variable s is in use. Nothing is differentiated and nothing
 * is differenced.
 *
 * A pleasant consequence worth stating because it looks like a coincidence:
 * < | grad psi |^2 / R^2 > is exactly < | q |^2 >. The wrapper below is
 * nonetheless written as gradient^2 / r^2, because it is the QUANTITY that is
 * being named and a reader checking it against the identity should not have to
 * re-derive the cancellation.
 *
 * THE ORDER OF EVERY QUANTITY HERE IS THE FLUX'S AND NOT THE POTENTIAL'S, AND
 * THAT IS THE ONE RESULT OF THIS STAGE THAT WAS NOT EXPECTED.
 *
 * ContourTracer's default pairing is psi* with q*, and psi* converges at k+2
 * where psi_h converges at k+1 -- so the natural expectation is that the
 * post-processed pairing buys an order in the averages too. IT DOES NOT. The
 * level set improves by an order and the WEIGHT does not, because the weight
 * divides by | grad psi | = r | q | and q* converges at k+1 like q_h: the local
 * post-processing buys its extra order in the POTENTIAL, and there is no k+2
 * flux to be had. An average built on both inherits the worse of the two.
 *
 * So every quantity in this file converges at k+1 with either pairing, and what
 * psi* buys is a CONSTANT -- a factor of a few. A constant is not an order and
 * the two must not be confused. This is the same shape as the band continuation
 * of B recorded in CLAUDE.md, where psi is continued by a SOLVED q at the
 * potential's own order and B cannot be, because differentiating an L2 field of
 * degree k leaves k-1 and no solved variable exists for it.
 *
 * The table is in tests/convergence/SurfaceAverageConvergence.cpp and the
 * numbers are not repeated here, because a measurement in a header goes stale
 * silently -- which is the stance FluxSurfaces.hpp takes for the same reason.
 *
 * AND THE AVERAGES ARE IMMUNE TO q's SIGN, WHERE THE TRACER IS NOT. Only | q |
 * enters, so handing this facility -q -- the raw DarcyForm block, which holds
 * exactly that -- changes nothing at all. FluxSurfaces.hpp records that the same
 * substitution traces the same curve backwards. Neither is a defect; the point
 * is that a sign error which the tracer survives visibly is invisible here, so
 * do not use these numbers as evidence about a sign.
 *
 * THE METRIC TRAP, AND WHAT IT DOES TO AN AVERAGE AS OPPOSED TO A LENGTH.
 *
 * INVERSION-PLAN.md section 3.2 warns that it is easy to build a spectrally
 * accurate rule and then feed it a SECOND-ORDER Jacobian obtained by
 * differencing neighbouring node positions, at which point the whole scheme is
 * second order and nothing in its output says so. IN-1 measured that on an arc
 * length -- 7.03 against 1.97 -- and tests/analytic/FluxSurfaceReference.hpp
 * measures it again on the same fixture, 1.08e-02 -> 3.11e-15 against
 * 1.37e-01 -> 6.84e-04.
 *
 * AN AVERAGE IS A DIFFERENT INTEGRAND FROM AN ARC LENGTH AND IS ENTITLED TO ITS
 * OWN MEASUREMENT, WHICH IS WHY THE CONTROL IS KEPT LIVE HERE TOO RATHER THAN
 * INHERITED. The obvious objection is that < X > is a RATIO of two integrals
 * over the same metric, so a bad metric ought to cancel. It does not: measured
 * on the analytic surface at Psi_N = 0.25 of SolovievEquilibrium::nstx(), the
 * differenced metric gives a relative error in < R^-2 > falling by a factor of
 * four per doubling of the node count from N = 32 on -- 2.6e-03 at N = 16 down
 * to 3.9e-06 at N = 512 -- against a pointwise metric that is at round-off from
 * N = 64 onwards. What the ratio buys is a smaller CONSTANT, about forty times
 * smaller than the same control does to V', and nothing whatever in the order.
 * A constant is not an order and the two must not be confused.
 *
 * differencedWeight and averageDifferenced() are therefore first-class members
 * and not test scaffolding. They are filled only by the equispaced-angle
 * builder, because that is the parametrisation the trap lives in;
 * differencedAvailable() says so and averageDifferenced() REFUSES rather than
 * returning a number when it is not.
 *
 * TWO EXTRACTIONS, WHICH IS TWO OF THE THREE LEGS SECTION 3.3 ASKS FOR.
 *
 * There are two builders and they share the field and the level and nothing
 * else:
 *
 *   surfaceAverages( tracer, fit )       equispaced in poloidal angle about the
 *                                        axis, each node a 1-D ray Newton, the
 *                                        metric | dx/dtheta | from IN-1's
 *                                        pointwise identity; periodic
 *                                        trapezoid, which for a smooth periodic
 *                                        integrand converges geometrically
 *   surfaceAverages( tracer, contour )   Gauss-Legendre on the cubic Hermite
 *                                        segments of the TRACED curve, whose
 *                                        nodes are wherever the tracer's step
 *                                        controller put them rather than
 *                                        equispaced in anything, and whose
 *                                        metric is the interpolant's own
 *                                        | dx/dt |
 *
 * One is a ray parametrisation and the other is a predictor-corrector trace
 * with an interpolant; they agree or one of them is wrong, and agreement is
 * worth more than either being plausible. There is a THIRD leg and IT IS NOT
 * BUILT: section 3.3's implicit quadrature, a rule constructed on the level set
 * with no curve extracted at all. It is deferred deliberately rather than
 * forgotten -- Saye's construction is for hyperrectangles and meq is on
 * triangles, so it would be Fries-Omerovic or moment fitting, which is a piece
 * of work rather than an afternoon. Read the two-route agreement below as two
 * legs of three.
 *
 * A THIRD, FREE CHECK THAT IS NOT A THIRD EXTRACTION. The angle builder's nodes
 * are laid out about a ray ORIGIN, and the integral cannot depend on it: a
 * reparametrisation of the same curve integrates to the same number. Fitting
 * the same contour about a deliberately displaced origin and getting the same
 * V' is a real check on the metric -- it is the invariance the identity
 * | dx/dtheta | = sqrt( rho'^2 + rho^2 ) exists to satisfy -- and it costs one
 * more fit. Measured at 4e-10 and better over displacements up to a third of
 * the minor radius. It does not make a third leg, because the field, the curve
 * and the rule are all shared.
 *
 * THE BAND: THE FLAG IS CARRIED, IT IS FALSE ON THE FITTED PATH, AND WHAT
 * CHANGES IS NAMED.
 *
 * On the curved path Omega_h is inscribed in Gamma, so there is a band O( h )
 * wide that is inside the plasma and outside the mesh, and the outermost
 * surfaces -- exactly the ones a q( psi ) profile most wants -- cross it.
 * INVERSION-PLAN.md section 4.3 is explicit that a COUNT is not enough and that
 * a consumer must be able to tell WHICH data came from an extension, the way
 * the .nc carries a byte extrapolated( Z, R ) beside inside; CLAUDE.md records
 * that the count-not-mask version of that was half of a real defect.
 *
 * So SurfaceNode::extended is read from the producing object per node --
 * AngleParametrisation::extended for the angle builder and
 * ContourPoint::extended for the contour one -- and SurfaceAverages::extended
 * is what a consumer drops a surface on before computing an error norm or
 * differencing two runs. A surface with extended true is limited by the
 * extension rather than by the discretisation and must be reported SEPARATELY:
 * quoting one rate over both populations would hide precisely what the band
 * stage exists to measure.
 *
 * EVERY ONE OF THEM READS FALSE ON THE FITTED PATH, which is not the same
 * statement as "the flag is not wired" and the difference is worth keeping. The
 * fitted path has no band -- Gamma_h IS Gamma there -- so a zero is the correct
 * answer rather than an absent one, and
 * tests/convergence/SurfaceAverageConvergence.cpp asserts the zero rather than
 * assuming it.
 *
 * THE CONTOUR BUILDER MARKS EACH GAUSS NODE FOR ITSELF, THROUGH
 * ContourTracer::sampleAt()'s seven-argument overload, which exists for this
 * caller. It used to mark a whole segment whenever either endpoint was
 * extended, which was conservative in one direction and wrong in the other: the
 * band does not respect the segment a node sits in, so a segment can have both
 * endpoints inside Omega_h and still cross Gamma_h in between. That
 * over-reported at the ends of a band excursion and UNDER-reported in the
 * middle of one, and under-reporting is the failure that matters -- it is a
 * band quantity presented as a solved one.
 *
 * The fit builder does not sample at all. AngleParametrisation keeps the
 * potential, the flux and the band flag its own ray Newton found at each node,
 * so this reads them rather than re-deriving them: cheaper, and in agreement
 * with the fit by construction rather than by coincidence, which on a
 * discontinuous field is not the same thing.
 *
 * THE REFERENCE-FREE ACCEPTANCE: THE AVERAGED EQUATION ITSELF.
 *
 * From < div G > = ( 1 / V' ) d/dpsi ( V' < G . grad psi > ) with
 * G = grad psi / R^2, using div( grad psi / R^2 ) = Delta* psi / R^2 and
 * Delta* psi = -F,
 *
 *     ( 1 / V' ) d/dpsi ( V' < | grad psi |^2 / R^2 > )  =  - < F / R^2 >
 *
 * -- three averages checking each other with nothing but the equation and no
 * reference value at all. THE RIGHT-HAND SIDE IS WRITTEN WITH THE F THE SOLVER
 * IS ACTUALLY FED and never as -mu_0 p' - g g' < R^-2 >: the second is
 * Solov'ev-specific and is re-derived by hand, so it is not independent of the
 * hand that derived it. That is the same discipline
 * SolovievEquilibrium::deltaStarFD() follows, and for the same reason.
 *
 * THE d/dpsi MUST BE RICHARDSON-EXTRAPOLATED, ( 4 D( h/2 ) - D( h ) ) / 3. A
 * plain central difference carries its own O( step^2 ) truncation and the
 * comparison is then floored by the INSTRUMENT rather than by the identity.
 * Measured twice. On the ANALYTIC nstx surface at Psi_N = 0.25, where the
 * identity is exact and only the difference is approximate, the extrapolation
 * converges at 4.000 in the step and the plain difference at 2.003, and the two
 * columns separate by a factor of 590 at 5% of | psi_ax | and 37,000 at 0.6%. On
 * the DISCRETE field at 2% of | psi_ax | the plain column is FLAT to three
 * figures -- 9.01e-06 at every mesh from n = 12 to n = 96 at k = 3 -- while the
 * extrapolated one falls to 2.80e-09. An identity checked with that instrument
 * would pass at 1e-5 with a real defect underneath it.
 * FluxDerivative::CentralDifference is kept so that the control is a live
 * column and not a claim.
 *
 * AND THE STEP IS A REAL CHOICE, NOT A DETAIL. Too large and the extrapolation's
 * own truncation dominates; too small and the difference of two nearly equal
 * integrals divides the surfaces' own DG-jump noise by the step and amplifies
 * it. Measured at k = 2, n = 96: 1.5e-08 at 2% of | psi_ax | against 3.4e-07 at
 * 1%, the SMALLER step being twenty times worse.
 *
 * WHAT THIS IS NOT.
 *
 * It is not valid on an open surface. Every rule here is periodic or closed and
 * assumes the contour returns to its start; an open surface terminating on the
 * domain boundary has endpoints and wants Chebyshev, which is IN-5 and is
 * deferred with free boundary.
 *
 * It does not check star-shapedness and does not need to. ContourTracer::
 * fitByAngle() refuses when its rays stop being transverse and reports
 * transversality when they do not; both are carried through onto
 * SurfaceAverages so that a caller reading an average never has to go back for
 * the hypothesis. The contour builder does not need the hypothesis at all.
 *
 * It is not an X-point-aware facility. Approaching the separatrix 1/| grad psi |
 * diverges, the surface develops a corner and every quantity here diverges
 * logarithmically with it. That is the physics, it is INVERSION-PLAN.md section
 * 8's second risk, and where meq cuts is a decision to be recorded rather than
 * discovered.
 *
 * And safetyFactor() carries no evidence of its own. It is
 * V' g < R^-2 > / 4 pi^2, an algebraic combination of two quantities that ARE
 * measured, times a g( psi ) the caller supplies from its own profile, so it
 * inherits their rates and nothing more. The suite pins its algebra and its
 * 4 pi^2 -- which is worth doing, since a constant is exactly as easy to get
 * wrong here as anywhere and no rate table would see it -- and that is all. It
 * is here because it is what ROADMAP.md item 10 reads, not because it is
 * independently checked.
 */

namespace meq
{

	/// One node of a surface quadrature. An integrand is a function of this and
	/// of nothing else, deliberately: no integrand should have to call back into
	/// a field and risk evaluating it at a slightly different point from the one
	/// the weight was built at.
	struct SurfaceNode
	{
		/// The parametrisation variable: the poloidal angle theta for the
		/// equispaced-angle builder, the accumulated segment parameter for the
		/// contour builder. Provenance for printing, not something an integrand
		/// should need.
		double parameter = 0.0;

		double r = 0.0;
		double z = 0.0;

		/// psi_h at the node, from the same evaluation the weight was built
		/// from, and | psi_h - level | with it. The second is
		/// INVERSION-PLAN.md section 2's error ( b ) for this node.
		double psi = 0.0;
		double residual = 0.0;

		/// q at the node, in meq's sign convention, and | grad psi | = r | q |
		/// with it. POINTWISE FROM THE SOLVED FLUX: see the header.
		double qR = 0.0;
		double qZ = 0.0;
		double gradient = 0.0;

		/// | dx / d(parameter) | at the node, pointwise. For the angle builder
		/// this is IN-1's sqrt( rho'^2 + rho^2 ); for the contour builder it is
		/// the cubic Hermite's own | dx/dt |.
		double metric = 0.0;

		/// 2 pi R | dx/ds | ds / | grad psi |, the whole quadrature weight, so
		/// that V' is a sum of these and < f > is a weighted mean with them.
		double weight = 0.0;

		/// The same two with the metric obtained by CENTRAL DIFFERENCING the
		/// neighbouring node positions -- the second-order control of section
		/// 3.2, filled only by the equispaced-angle builder. Never used by
		/// average(); a caller asks for it deliberately.
		double differencedMetric = 0.0;
		double differencedWeight = 0.0;

		/// Whether the field at this node came from the band extension outside
		/// the mesh rather than from an element. Read from the producing
		/// object -- AngleParametrisation::extended per node, ContourPoint::
		/// extended per segment endpoint -- and false throughout on the fitted
		/// path, where there is no band. See the header on why this is per node
		/// and not a count.
		bool extended = false;
	};

	/// An integrand. Takes a node and returns the value of the quantity being
	/// averaged there.
	using SurfaceIntegrand = std::function<double( SurfaceNode const & )>;

	/// A closed flux surface, sampled, with the weights that make it a
	/// quadrature. This is the discrete counterpart of
	/// meq::analytic::SurfaceQuadrature, which does the same on the exact field
	/// and is what these numbers are measured against.
	struct SurfaceAverages
	{
		/// The level psi was rooted at, carried through from the contour.
		double level = 0.0;

		std::vector<SurfaceNode> nodes;

		/// V' = closed-integral 2 pi R dl / | grad psi |, which is the sum of
		/// the node weights. Computed at build because everything else is
		/// relative to it.
		double vPrime = 0.0;

		/// The same from the differenced metric. THE CONTROL: second order, and
		/// nothing about the rule it came from says so. Zero, and not to be
		/// read, when differencedAvailable() is false.
		double vPrimeDifferenced = 0.0;

		/// How many nodes carried a field value from an extension, whether any
		/// did, and how far outside Gamma_h the deepest of them sat. A surface
		/// with extended true is limited by the extension rather than by the
		/// discretisation and must be reported separately; the depth is the
		/// length the extension's own remainder is taken over, so it is what a
		/// band error is read against.
		int extendedNodes = 0;
		bool extended = false;
		double deepestBandNode = 0.0;

		/// Which extension answered those calls, carried through from the
		/// tracer so that a stored surface says what its band nodes were
		/// computed with rather than leaving it to be remembered.
		/// BandExtension::None on the fitted path.
		BandExtension bandExtension = BandExtension::None;

		/// Carried through from ContourTracer::fitByAngle(): min | u x t | over
		/// the fit, and whether it cleared the floor. Zero and true for the
		/// contour builder, which needs no such hypothesis.
		double transversality = 0.0;
		bool transverse = true;

		/// The worst | psi_h( x ) - level | over the nodes. Error ( b ) for the
		/// whole surface, and the thing to read before believing any of the
		/// rest.
		double worstResidual = 0.0;

		/// Times the element walk fell back on Mesh::FindPoints. A performance
		/// statement and nothing else -- CLAUDE.md records FindPoints as
		/// O( elements x points ) -- and nothing else would say it.
		///
		/// IT IS THE PRODUCING ROUTINE'S COUNT AND NOT THIS QUADRATURE'S, which
		/// is a limitation rather than a choice: ContourTracer::sampleAt() does
		/// not report whether its walk fell back, so the sampling loops here
		/// cannot tally their own. What is carried through is the fit's count
		/// for the angle builder and the trace's for the contour builder.
		int fallbackLocations = 0;

		/// Whether the differenced control was filled. True for the
		/// equispaced-angle builder and false for the contour builder, where
		/// there is no equispaced neighbour to difference against. A flag rather
		/// than a test on vPrimeDifferenced, because a control that reports
		/// itself absent by being zero is a control that reports itself present
		/// the moment somebody writes a zero into it.
		bool differencedControl = false;

		bool differencedAvailable() const
		{
			return differencedControl;
		}

		/// closed-integral 2 pi R f dl / | grad psi |. What an integrand that is
		/// already a density wants, and what V' d/dpsi terms are built from.
		double integrate( SurfaceIntegrand const &f ) const;

		/// < f > = integrate( f ) / V'.
		double average( SurfaceIntegrand const &f ) const;

		/// The same two through the differenced metric. THE CONTROL.
		/// @throws std::runtime_error when differencedAvailable() is false --
		///         refusing rather than returning a number silently built on the
		///         pointwise metric, which would make the control agree with the
		///         answer for the worst possible reason.
		double integrateDifferenced( SurfaceIntegrand const &f ) const;
		double averageDifferenced( SurfaceIntegrand const &f ) const;

		/// < R^-2 >. The geometric one, and the factor RoPP (142)'s safety
		/// factor is built from.
		double inverseRSquared() const;

		/// < | grad psi |^2 / R^2 >, which is exactly < | q |^2 > and is written
		/// the long way for the reason in the header. This is the one the flux
		/// earns its keep on, and it is the left-hand side of the averaged
		/// Grad-Shafranov identity.
		double gradPsiSquaredOverRSquared() const;

		/// The arc length of the surface, closed-integral dl, through the same
		/// facility.
		///
		/// BE PRECISE ABOUT WHAT THIS CHECKS, BECAUSE IT LOOKS LIKE MORE THAN IT
		/// IS. Since every weight is 2 pi R | dx/ds | ds / | grad psi |, this
		/// integrand recovers | dx/ds | ds exactly and the result is
		/// ALGEBRAICALLY equal to the producing object's own length --
		/// AngleParametrisation::length() for the angle route and
		/// Contour::hermiteLength() for the contour one. So it is not an
		/// independent measurement of the metric and must not be quoted as one.
		/// What it does check, exactly, is that the weights are built the way
		/// this header says they are: a missing r, a wrong dtheta or a gradient
		/// on the wrong side of the division all break the identity.
		double arcLength() const;

		/// The safety factor, RoPP (142):
		/// safetyFactor = V' g < R^-2 > / 4 pi^2, with g( psi ) = R B_toroidal
		/// supplied by the caller from its own profile because a Source does not
		/// carry it. NOT q, and never named q -- see the header. It carries no
		/// evidence of its own: it inherits the rates of V' and < R^-2 >, and
		/// what the suite pins is its algebra and its 4 pi^2.
		double safetyFactor( double g ) const;
	};

	/**
	 * Build a quadrature from an equispaced-in-poloidal-angle fit.
	 *
	 * THE PRIMARY ROUTE. The nodes are @a fit's, the metric is IN-1's pointwise
	 * sqrt( rho'^2 + rho^2 ), and the rule is the periodic trapezoid -- which
	 * for a smooth periodic integrand converges geometrically, and on a DISCRETE
	 * contour converges very fast and then floors at the DG jump, because
	 * psi_h is discontinuous across faces and no quadrature is geometric on a
	 * function with jumps. FluxSurfaces.hpp records that floor and the control
	 * that identifies it as the field rather than the rule.
	 *
	 * psi and q are RE-READ at each node through @a tracer's own seam, so that
	 * the field behind the weights is the same field the fit rooted, located by
	 * the same element walk. That is one extra evaluation per node and it is
	 * avoidable: AngleParametrisation computes q at every node to build rho' and
	 * does not store it. The band flag does NOT come through this call -- it is
	 * read from AngleParametrisation::extended, because sampleAt() answers a
	 * band point without saying that it did.
	 *
	 * @throws std::invalid_argument if the fit is empty.
	 */
	SurfaceAverages surfaceAverages( ContourTracer const &tracer,
	                                 AngleParametrisation const &fit );

	/**
	 * Build a quadrature from the traced contour itself, by Gauss-Legendre on
	 * its cubic Hermite segments.
	 *
	 * THE INDEPENDENT ROUTE. It shares the field and the level with the fit
	 * above and nothing else: the nodes are wherever the tracer's step
	 * controller put them rather than equispaced in anything, the metric is the
	 * interpolant's own | dx/dt | rather than a ray identity, and no
	 * star-shapedness is assumed anywhere. Agreement
	 * between the two is the acceptance; disagreement would not say which is
	 * wrong, which is why both are kept.
	 *
	 * It fills no differenced control: there is no equispaced neighbour to
	 * difference against, and a chord-based substitute would be a different
	 * control answering a different question.
	 *
	 * @param gaussPoints per segment. Eight is exact for a degree-15 integrand
	 *        and the cost is one field evaluation each.
	 * @throws std::invalid_argument if the contour is not closed or has no
	 *         segments, or @a gaussPoints is below one.
	 */
	SurfaceAverages surfaceAverages( ContourTracer const &tracer,
	                                 Contour const &contour,
	                                 int gaussPoints = 8 );

	/// How the d/dpsi of the averaged equation is taken.
	enum class FluxDerivative
	{
		/// ( 4 D( step/2 ) - D( step ) ) / 3. THE DEFAULT, and see the header:
		/// without it the identity is floored by the instrument.
		Richardson,

		/// A plain central difference. THE CONTROL, kept live so that the
		/// paragraph above is a measurement rather than a claim. Never an
		/// answer.
		CentralDifference
	};

	/// The result of checking the flux-surface-averaged Grad-Shafranov equation
	/// on the discrete field. Everything needed to print the comparison, not
	/// only the residual, because a residual with no scale beside it is not a
	/// reading.
	struct AveragedEquationResidual
	{
		/// ( 1 / V' ) d/dpsi ( V' < | grad psi |^2 / R^2 > ).
		double leftHandSide = 0.0;

		/// - < F / R^2 >, from the F the solver is fed.
		double rightHandSide = 0.0;

		/// leftHandSide - rightHandSide, in the units of F / R^2.
		double residual = 0.0;

		/// d/dpsi ( V' < | grad psi |^2 / R^2 > ) itself, and V' at the level,
		/// so that a reader can see which of the two the residual is a small
		/// difference of.
		double derivative = 0.0;
		double vPrime = 0.0;

		/// The worst | psi_h - level | over the five surfaces the difference
		/// touched, and whether any of them crossed the band.
		double worstResidual = 0.0;
		bool extended = false;
	};

	/**
	 * The flux-surface average of the Grad-Shafranov equation itself, as a
	 * residual that a correct set of averages drives to zero.
	 *
	 * It needs no reference value: three averages check each other with nothing
	 * but the equation. On the DISCRETE field it is not expected to reach
	 * round-off -- the identity holds exactly for psi and only to the
	 * discretisation error for psi_h -- so the property to assert on is that the
	 * residual CONVERGES with h, not that it is small.
	 *
	 * Five surfaces are traced: the level itself and the four the Richardson
	 * difference needs. Three with FluxDerivative::CentralDifference, which is
	 * cheaper and is not the answer.
	 *
	 * @a step is in the units of psi and it is a real choice, NOT MONOTONE in
	 * either direction: too large and the extrapolation's own truncation
	 * dominates, too small and the difference of two nearly equal integrals
	 * divides the surfaces' DG-jump noise by the step and amplifies it. Two
	 * percent of | psi_ax - psi_bnd | is what was measured to be near the
	 * optimum on the fixture in the suite; one percent is twenty times worse.
	 *
	 * @param f the source, F( r, z, psi ) -- THE F THE SOLVER IS FED, per the
	 *          header. Not a re-derived right-hand side.
	 * @throws whatever ContourTracer::traceFromAxis() and fitByAngle() throw,
	 *         which is how a level too close to the separatrix or a surface that
	 *         is not star-shaped announces itself.
	 */
	AveragedEquationResidual averagedGradShafranovResidual(
		ContourTracer const &tracer, CriticalPoint const &axis, double level,
		std::size_t angles, double step,
		std::function<double( double r, double z, double psi )> const &f,
		FluxDerivative derivative = FluxDerivative::Richardson );

}

#endif
