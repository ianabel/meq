#ifndef MEQ_FLUXSURFACES_HPP
#define MEQ_FLUXSURFACES_HPP

#include <cstddef>
#include <vector>

#include "mfem.hpp"

#include "CriticalPoints.hpp"
#include "GradShafranov.hpp"

/*
 * Flux surfaces as curves: the predictor-corrector tracer of INVERSION-PLAN.md
 * stage IN-0, and the poloidal-angle parametrisation and metric of stage IN-1.
 *
 * THE BAND BETWEEN Gamma_h AND Gamma IS THE SECOND HALF OF IN-0, AND IT IS THE
 * INTERESTING PART RATHER THAN A CORNER CASE.
 *
 * On the curved path Omega_h is the union of background elements lying INSIDE
 * Gamma, so Gamma_h is inscribed and there is a band O( h ) wide that is inside
 * the plasma and outside the mesh. The outermost closed surface is psi = 0,
 * which IS Gamma, so the surfaces the band affects are exactly those with
 * Psi_N -> 1 -- which is where q( psi ) and the flux-surface averages are most
 * wanted and least forgiving.
 *
 * setBandExtension() says which field answers a call there, and
 * BandExtension::None -- THE DEFAULT -- is the fitted path unchanged: a trace
 * that reaches the edge of the mesh returns ContourStatus::LeftMesh and says
 * so. That is deliberate rather than timid. v0-legacy:FluxSurfaces.cpp printed
 * "Terminating because curve left domain" and returned a partial curve, which
 * for an outer surface is an arc labelled as a closed contour, and a
 * flux-surface average over it is wrong by an amount nothing reports.
 *
 * TWO EXTENSIONS, MEASURED SIDE BY SIDE, AND THE ORDERS ARE STRUCTURAL RATHER
 * THAN A PROPERTY OF THE BENCHMARK.
 *
 *   BandExtension::FluxTaylor    psi( p ) = psi( x0 ) + r0 q( x0 ) . ( p - x0 )
 *                                from the foot x0 of p on Gamma_h, which is
 *                                what GridSampler::samplePotentialWithFlux()
 *                                does for the .nc file. NOTHING IS EVER
 *                                EVALUATED OUTSIDE AN ELEMENT, which is the
 *                                property that made it right there. Its error
 *                                is | psi_h - psi | at the foot, O( h^(k+1) ),
 *                                PLUS the Taylor remainder over a band of width
 *                                O( h ), which is O( h^2 ) whatever k is. So it
 *                                is exact in order at k = 1 and caps the band
 *                                at second order from k = 2 on. And because
 *                                grad psi is frozen at the foot the extended
 *                                field is AFFINE: contours in the band are
 *                                straight lines, at every k.
 *
 *   BandExtension::TransferLift  psi( p ) = g( a( x0 ) )
 *                                           + integral over the segment from p
 *                                             to a( x0 ) of ( -r q ) . m ds,
 *                                with q outside the mesh supplied by the
 *                                method's OWN extension operator E_h -- the
 *                                element's polynomial evaluated outside it,
 *                                mfem::ElementExtension. This is the extension
 *                                technique's own construction, and it is the
 *                                identity the whole curved path rests on:
 *                                C u = -grad_bar( psi ) is a gradient, so the
 *                                line integral of it from p to any point of
 *                                Gamma is psi( p ) - g there. Its error is the
 *                                error of q integrated over a path of length
 *                                O( h ) -- so O( h^(k+2) ) in principle, one
 *                                order BETTER than psi_h's own, and the band
 *                                does not limit the contour at all. Measured on
 *                                Gamma itself it reaches k+2 at k = 2 and 3 and
 *                                about k+1 at k = 1, where every term in sight
 *                                is second order anyway.
 *
 * BOTH LIFT WHICHEVER PAIR THE TRACER WAS GIVEN, so Potential::PostProcessed
 * integrates q* and steps with it while Potential::Raw integrates q. The rates
 * are measured on the RAW pair, because k+1 and k+2 are the two orders being
 * told apart and a study of the band has to say which field it is the band of.
 *
 * IT IS A REAL DIFFERENCE AND NOT A REFINEMENT OF THE SAME IDEA. The Taylor
 * step is affine and never leaves an element; the lift is curved and leaves one
 * deliberately, through the operator the analysis of the method is written for.
 * Which wins, by how much, and at what order each converges, is measured in
 * tests/convergence/FluxSurfaceConvergence.cpp with the two run over the SAME
 * contours -- and the Taylor column is kept as the control in the way
 * ExtensionConvergence.cpp keeps its pinned-zero column: if it ever converges
 * as fast as the lift, the comparison is empty and the test is worthless.
 *
 * THE EXACT BOUNDARY IS AN ANCHOR AND NEVER A CORRECTION, WHICH IS WHY THE
 * LIFT IS ANCHORED ON Gamma AND NOT BLENDED TOWARDS IT. psi = 0 exactly on
 * Gamma and Gamma is known analytically, so the outermost contour is known for
 * free with no tracing at all. What does NOT work is scaling a bad
 * extrapolation towards it: CLAUDE.md records that ( 1 - t ) v scales a
 * positive value down and never changes its sign, so all seventeen offending
 * nodes of the .nc survived the blend. The error was in WHERE the field was
 * evaluated, not in how it was weighted. TransferLift fixes the evaluation --
 * g( a( x0 ) ) is the lower limit of an integral, not a weight.
 *
 * A POINT IN THE BAND SAYS SO, PER POINT, AND A COUNT WOULD NOT BE ENOUGH.
 * ContourPoint::extended is the flag, Contour::extendedPoints the count and
 * Contour::crossesBand() the predicate; AngleParametrisation carries the same
 * pair for its own nodes. That is the .nc's byte extrapolated( Z, R ) beside
 * inside, and it is here for the reason it is there: a band point holds real
 * data and is therefore indistinguishable from a solved one without being
 * told. CLAUDE.md records extrapolated_nodes having been a COUNT and not a
 * MASK as the other half of a real defect -- nothing downstream could tell
 * WHICH. IN-2 must be able to report band-crossing surfaces separately rather
 * than averaging them in with the others, and that needs the mask.
 *
 * WHAT A TRACED CONTOUR'S ERROR IS MADE OF, WHICH IS THE FRAME FOR EVERY
 * MEASUREMENT IN tests/convergence/FluxSurfaceConvergence.cpp.
 *
 * INVERSION-PLAN.md section 2 splits it three ways and the split is the whole
 * reason this is worth doing carefully:
 *
 *   (a) FIELD.          the level set of the field being rooted, against
 *                       { psi = c }. O( h^(k+1) / |grad psi| ) for psi_h and
 *                       O( h^(k+2) / |grad psi| ) for psi*, by the implicit
 *                       function theorem. Not ours to improve here -- it is the
 *                       discretisation, and WHICH POTENTIAL TO ROOT below is
 *                       the one lever this file has on it.
 *   (b) POINT LOCATION. how far an accepted point is from { psi_h = c }. This
 *                       is what the corrector buys: at the corrector tolerance
 *                       and, the property that matters, INDEPENDENT OF PATH
 *                       LENGTH. Measured flat over ten circuits.
 *   (c) REPRESENTATION. how far the interpolant BETWEEN accepted points is from
 *                       the curve. Governed by the point spacing and the
 *                       interpolation order, and DECOUPLED FROM h AND k
 *                       ENTIRELY.
 *
 * Conflating them is how this problem gets done badly, so they are measured as
 * three separate numbers rather than as one.
 *
 * WHY PREDICTOR-CORRECTOR, AND WHY THE FIRST-ORDER PREDICTOR IS NOT A DEFECT.
 *
 * The corrector is a root find, so every accepted point sits on the discrete
 * level set to tolerance whatever the predictor did. There is nothing for a
 * higher-order predictor to buy in ACCURACY; what it buys is STEP LENGTH. That
 * is why explicit Euler plus Newton is the norm, and why the first-order
 * predictor in v0-legacy:FluxSurfaces.cpp is not the defect an earlier reading
 * of that file called it.
 *
 * What the corrector buys is that error off the contour does not ACCUMULATE.
 * Integrate the tangent field without projection and the curve drifts onto a
 * neighbouring contour and does not close; with the projection each point is
 * independently on the curve and there is no accumulation at all. Allgower &
 * Georg state the same preference from the other side: "the general opinion is
 * that it is preferable to exploit the contractive properties of the zero set
 * H^-1( 0 ) relative to such iterative methods as those of Newton type."
 *
 * It also makes the DG jump a near-non-issue. psi_h disagrees with itself by
 * O( h^(k+1) ) across a face; a tracer that steps across one lands in whichever
 * element it lands in and the corrector re-anchors it there. Compare an ODE
 * integrator, for which a face is a discontinuity violating the smoothness
 * every Runge-Kutta order derives from.
 *
 * "NEAR-NON-ISSUE" RATHER THAN "NON-ISSUE", AND THE REMAINDER IS A REAL ONE
 * THAT HAD TO BE MET BEFORE IT WAS BELIEVED. { psi_h = c } is a union of
 * per-element arcs offset from each other by the jump, so a point landing
 * within jump/|grad psi| of a face may be on NEITHER arc: the Newton step
 * computed in element A pushes it across into B, the step computed in B pushes
 * it back, and the residual alternates between two values without ever falling
 * below a tolerance tighter than the jump. It is rare -- it needs a point to
 * land inside a band about 2e-5 wide on a contour of length 1.7 -- and it gets
 * commoner as Delta_s falls, simply because more points are placed. Left alone
 * it ends a trace, and it ended several before it was diagnosed.
 *
 * The corrector therefore keeps its BEST iterate and accepts it when four steps
 * running fail to improve on it, refusing to let the iterate travel further
 * than one predictor step from where it started. Contour::stalledCorrections
 * counts those points and ContourPoint::residual says how close each got, to be
 * read against Contour::correctorTarget. That is honest rather than convenient:
 * on a discontinuous field a point cannot be closer to "the" level set than the
 * field's own ambiguity at a face.
 *
 * NO Mesh::FindPoints IN THE INNER LOOP, AND THE LAST-RESORT CALLS ARE COUNTED.
 *
 * CLAUDE.md records mfem::Mesh::FindPoints as O( elements x points ): a
 * brute-force scan over element centres. A tracer calls for a point location
 * once per corrector iteration, so using it there makes a trace quadratic in
 * the mesh. The walk instead is ElementTransformation::TransformBack into the
 * element the previous evaluation used, then its face neighbours from
 * Mesh::ElementToElementTable(), then a second ring, and only then FindPoints.
 *
 * Contour::fallbackLocations COUNTS the last-resort calls, and a caller who
 * sees it climbing is being told the walk is failing -- which is a performance
 * statement, not a correctness one, since FindPoints returns the same element.
 * On the benchmarks it is zero.
 *
 * THE STEP IS CURVATURE CONTROLLED, AND THE STANDARD STRATEGIES OPTIMISE THE
 * WRONG OBJECTIVE FOR US.
 *
 * Allgower & Georg ch. 6 details two steplength strategies and both are keyed
 * on CORRECTOR EFFORT -- Georg (1983) and Den Heijer & Rheinboldt (1981).
 * Continuation software wants to traverse a path cheaply, so equidistributing
 * corrector work is right for it. We want a well-shaped CURVE, so we want to
 * equidistribute INTERPOLATION error, which means curvature control. Take a
 * constant turning angle per step, Delta_s = targetTurn / kappa, and the
 * product kappa * Delta_s is then constant -- which by the Hermite error below
 * makes the deviation per segment proportional to Delta_s and the deviation
 * relative to the segment's own length constant. Equidistributing the deviation
 * itself would want Delta_s ~ kappa^(-3/4); the plan asks for 1/kappa, which is
 * the robust geometric variant and is what is implemented.
 *
 * kappa IS A DIFFERENCE -- the turning of the unit tangent between the last two
 * accepted points, divided by the step between them -- AND THAT IS FINE, for
 * exactly the reason CriticalPoints.hpp gives for its differenced Jacobian: it
 * STEERS the iteration and appears nowhere in its fixed point. Every accepted
 * point is corrected onto the level set to tolerance whatever the step
 * controller did, so a wrong kappa costs a badly shaped point distribution and
 * cannot move a single accepted point. The step is bounded from both sides -- a
 * floor, so a momentary large kappa cannot stall the trace, and a ceiling that
 * is a fraction of the local element size, so a step cannot leap a whole
 * element and hand the walk a point it has to search for.
 *
 * CUBIC HERMITE BETWEEN ACCEPTED POINTS, WITH THE TANGENTS FROM q.
 *
 * INVERSION-PLAN.md section 3.2( i ), and the highest-value single use of q in
 * the whole item. q gives the unit tangent at every accepted point at the
 * potential's own order and costs nothing, the tracer having evaluated it
 * already for the predictor. Two points and two unit tangents, with the tangent
 * magnitudes taken as the CHORD LENGTH, give
 *
 *     x( t ) = h00( t ) x_i + h10( t ) L T_i + h01( t ) x_(i+1) + h11( t ) L T_(i+1)
 *
 * and that is fourth order in the step. Chord scaling is not the optimal
 * scaling and the difference is worth recording, because the arithmetic looks
 * as though it should cost an order: on a unit-radius arc of angle theta the
 * exact optimal Bezier control offset is ( 4/3 ) tan( theta/4 ) and the chord
 * scaling gives ( 4/3 ) sin( theta/4 ) cos( theta/4 ), which differ by a
 * relative O( theta^2 ). It does not cost an order. The midpoint deviation
 * works out at
 *
 *     1 - cos( theta/2 ) - ( 1/2 ) sin^2( theta/2 ) = theta^4 / 128 + O( theta^6 ),
 *
 * so the deviation is kappa^3 Delta_s^4 / 128 -- fourth order in the spacing,
 * against second order for the straight chords the same points would otherwise
 * be joined by. Both are measured, side by side, in one table, with the chord
 * column kept as the control in the way ExtensionConvergence.cpp keeps its
 * pinned-zero column.
 *
 * BUT THE TANGENT FROM q IS NOT THE TANGENT OF THE CURVE BEING INTERPOLATED,
 * AND THAT COSTS THE FOURTH ORDER ON THE WRONG PAIRING. The tangent of
 * { psi_h = c } is built from grad( psi_h ), and q_h is a separate unknown that
 * agrees with grad( psi_h )/r only to O( h^k ). So an interpolant carrying q's
 * tangents passes through points of { psi_h = c } with the tangents of a
 * slightly different curve, and measured against { psi_h = c } it is fourth
 * order until that tilt takes over and second order afterwards. With psi* the
 * tilt is O( h^(k+1) ) and the fourth order survives; see WHICH POTENTIAL TO
 * ROOT below, which is where the two facts meet. Against the TRUE contour q's
 * tangent is the better of the two by a full order either way, which is why it
 * is what the tracer stores.
 *
 * CLOSURE: "AT MACHINE PRECISION" NEEDS A REFINEMENT, AND HERE IT IS.
 *
 * INVERSION-PLAN.md asks for "closure error over a full circuit at machine
 * precision". Taken literally that is not available and it is not a defect that
 * it is not. psi_h is discontinuous across faces, so { psi_h = c } is a union of
 * per-element arcs with O( h^(k+1) ) jumps between them; the curve is closed as
 * a set, but it is not one analytic curve. What IS true, and what is measured:
 *
 *   * TANGENTIALLY the error is whatever the final step happens to leave, which
 *     is O( Delta_s ) if nothing is done about it. That is why the final step is
 *     SHORTENED to aim at the start point: the predictor is given the length
 *     ( x_start - x_last ) . t rather than Delta_s, so the remaining tangential
 *     gap g_t comes out at O( kappa^2 Delta_s^3 ) rather than O( Delta_s ). It
 *     is reported as Contour::closureTangential;
 *   * NORMALLY the gap is bounded by the curve BENDING over that tangential
 *     gap, kappa g_t^2 / 2, and is independent of path length. It is reported
 *     as Contour::closureNormal.
 *
 * THE SECOND OF THOSE WAS FIRST WRITTEN AS "at the corrector tolerance" AND
 * THAT WAS WRONG, in the same way and for the same kind of reason as the
 * original "at machine precision" it was replacing. Both the final corrected
 * point and the start point are on { psi_h = c } to the corrector's 1e-13; but
 * they are separated ALONG the curve by g_t, and an arc departs from its own
 * tangent line by kappa g_t^2 / 2 over that distance. So the normal gap is
 * geometry, not drift. Measured on the benchmark it sits at about a quarter of
 * that bound and does not move when the path length is doubled, and the test
 * asserts against the bound rather than against a tolerance -- because a tracer
 * that was integrating rather than projecting would show a normal gap growing
 * with path length and unrelated to g_t, which is exactly the failure the
 * corrector exists to prevent.
 *
 * The measurement is taken BEFORE the curve is closed, and the closing point is
 * then set to the start point exactly, so that the polyline a consumer receives
 * has no stub and no duplicate-but-not-quite node.
 *
 * FACE CROSSINGS ARE COUNTED AND THEIR JUMPS MEASURED, AND THE JUMP IS ERROR
 * ( a ) RATHER THAN A TRACER DEFECT.
 *
 * Where two consecutive accepted points are in different elements the segment
 * between them crosses a face, and psi_h evaluated from the two sides at the
 * crossing point disagrees by the DG jump. Contour::faceCrossings and
 * Contour::worstFaceJump report it. It is a measurement of the DISCRETISATION
 * -- it converges at the order of whichever field is being rooted, k+1 for
 * psi_h and k+2 for psi*, and is a component of error ( a ) -- and nothing the
 * tracer does can reduce it. It is reported because it is the floor under
 * measurement ( c ): a Hermite segment whose endpoints are corrected onto two
 * different arcs cannot be closer to either than the two arcs are to each
 * other, so refining Delta_s on a fixed mesh drives the representation error
 * down at fourth order until it meets that floor and then stops.
 *
 * THE SEAM: EVERY EVALUATION OF psi AND OF q GOES THROUGH ONE FUNCTION.
 *
 * ContourTracer::sampleField() is the only place in this file where psi_h or
 * q_h is read at a physical point. The predictor, the corrector, the tangent,
 * the curvature and the ray Newton of the angle parametrisation all call it.
 * That is deliberate, and it is what made the band a change to ONE function:
 * both extensions above are answered there and nowhere else, and every counter,
 * every status and every flag downstream of them came along for free.
 *
 * The seam held. What it did NOT anticipate is that a band point has no
 * element containing it, so FieldSample::element now carries the element whose
 * polynomial was extended -- the one owning the nearest face of Gamma_h --
 * rather than the element the point is in. Consecutive points in that element
 * therefore look like a face crossing when they are not, so the face-jump
 * measurement is skipped wherever either endpoint is extended. It is the DG
 * jump of psi_h that is being measured there and a band point has no second
 * side to disagree with.
 *
 * MFEM's PathLiftCoefficient IS NOT THE PRIMITIVE THE LIFT IS BUILT ON, AND
 * INVERSION-PLAN.md section 4.3 NAMES IT AS THOUGH IT WERE. It evaluates on a
 * FaceElementTransformations of a face of Gamma_h, at that face's own
 * integration point, and returns the lifting to the foot -- it answers "what is
 * phi_h on Gamma_h", which is the question the estimator's eta_5 asks and is
 * not this one. The usable primitive is one level down and is public:
 * mfem::PathIntegral( Cu, x, xbar, line_ir ) takes ARBITRARY endpoints, so it
 * integrates from a band point rather than from a face point, and
 * mfem::ElementExtension supplies the flux along the way. So the answer to
 * "is the transfer-path lift usable at an arbitrary band point" is yes, but
 * not through the class the plan names.
 *
 * WHICH POTENTIAL TO ROOT: MEASURED, AND THE MEASUREMENT MOVED THE DEFAULT AND
 * CHANGED WHAT THE TANGENT FROM q IS WORTH.
 *
 * The corrector roots a scalar field and takes its direction from a flux, and
 * MEQ has two candidate pairs: psi_h with q_h, which converge at k+1, and
 * psi*_h with q*_h from DarcyForm::Reconstruct(), where psi* converges at k+2
 * and is what MEQ reports everywhere else. Three things were measured and all
 * three point the same way; tests/convergence/FluxSurfaceConvergence.cpp
 * carries the tables, and the numbers are not repeated here because a
 * measurement in a header goes stale silently.
 *
 * 1. psi* CARRIES ITS OWN ORDER INTO THE CONTOUR. Every accepted point is on
 *    { psi* = c } to tolerance, so by the implicit function theorem the traced
 *    curve is k+2 from the true one rather than k+1 -- measured, one to two
 *    orders of magnitude closer on the same mesh, and at fewer corrector
 *    iterations per point rather than more.
 *
 * 2. THE CONSISTENCY ARGUMENT FOR IT IS TRUE BY ONE ORDER AND NOT EXACTLY, AND
 *    IT IS NOT TRUE FOR THE REASON USUALLY GIVEN. The reason to expect psi* to
 *    pair better is that the local post-processing is built so that its
 *    gradient matches the flux. Checked rather than repeated: Stenberg's local
 *    problem is driven by the reconstructed TOTAL flux qhat_h in RT_k -- the
 *    normally continuous field the constraint equation projects onto -- and not
 *    by q*_h, so grad( psi* ) and r q* are NOT the same object. What is
 *    measured is that they agree to O( h^(k+1) ) where grad( psi_h ) and r q_h
 *    agree only to O( h^k ), which is the whole point of a mixed method read
 *    from an unexpected direction: q converges at k+1 and a differentiated L2
 *    potential at k.
 *
 * 3. AND THAT ORDER IS WHAT THE CUBIC HERMITE NEEDS, WHICH IS THE FINDING THIS
 *    STAGE DID NOT EXPECT. The interpolant is built on tangents from the flux
 *    and is measured against the level set of the potential. Those are the same
 *    curve only to the extent that the two fields agree. Paired with psi_h the
 *    Hermite is fourth order in Delta_s until the O( h^k ) tangent tilt takes
 *    over and second order afterwards; paired with psi* the tilt is a full
 *    order smaller and the fourth order survives the whole sweep. So the
 *    default is not merely more accurate, it is what makes the section 3.2( i )
 *    claim about q hold at all on this discretisation.
 *
 * Potential::PostProcessed is therefore the default and the constructor refuses
 * a solver that has not been post-processed, because postProcessedPotential()
 * returns an empty field rather than complaining. Potential::Raw is the other
 * door, needs nothing but a solve, and is what the field-error rate of IN-0 is
 * pinned on -- k+1 being psi_h's order and k+2 being psi*'s, a test that wants
 * to check the implicit function theorem has to say which.
 *
 * WHAT THIS IS NOT.
 *
 * It is not an exhaustive extraction of the level set. trace() follows the
 * connected component through the point it is given and says nothing about any
 * other component -- a second, disjoint island at the same level is not found,
 * not reported and not looked for. That is the same disclaimer
 * CriticalPoints.hpp makes about seeded Newton, and for the same reason: the
 * exhaustive construction is subdivision with a convex-hull test and it is not
 * built.
 *
 * It is not an X-point-aware tracer. The level set through a saddle is not a
 * 1-manifold and the tangent is not defined there; a trace at that level will
 * stall or turn a corner arbitrarily. INVERSION-PLAN.md section 5 is why
 * CriticalPointFinder exists, and a caller is expected to have used it.
 *
 * ARC LENGTH AND THE METRIC -- IN-1, AND THE TRAP MADE INTO A MEASUREMENT.
 *
 * INVERSION-PLAN.md section 3.2 warns that it is easy to build a spectrally
 * accurate quadrature rule and then feed it a SECOND-ORDER Jacobian obtained by
 * differencing neighbouring node positions, at which point the whole scheme is
 * second order and nothing in the output says so. AngleParametrisation -- the
 * struct below -- exists to measure exactly that, with the differenced version
 * kept as the control.
 *
 * The parametrisation is by poloidal angle about the magnetic axis a. With
 *
 *     x( theta ) = a + rho( theta ) u( theta ),   u = ( cos theta, sin theta ),
 *
 * so that u' = ( -sin theta, cos theta ) and u . u' = 0, the curve tangent t is
 * known POINTWISE from q -- t = ( -q_z, +q_r ) / |q| -- and dx/dtheta must be
 * parallel to it. Writing a x b := a_r b_z - a_z b_r for the scalar cross
 * product in the plane,
 *
 *     dx/dtheta = rho' u + rho u',
 *     ( rho' u + rho u' ) x t = 0,
 *     rho' ( u x t ) = -rho ( u' x t ),
 *
 * and since u x t = u' . t and u' x t = -u . t this is
 *
 *     rho'( theta ) = rho ( u . t ) / ( u' . t ),
 *     | dx/dtheta |  = sqrt( rho'^2 + rho^2 ).
 *
 * POINTWISE, FROM q, WITH NOTHING DIFFERENCED. Two checks that it is right,
 * both of which a transcription error fails: on a circle about a the tangent is
 * perpendicular to the radius, u . t = 0 and rho' = 0; on the straight line
 * z = 1 with a at the origin, rho = 1 / sin theta and the identity returns
 * -cos theta / sin^2 theta, which is rho' differentiated by hand. It is also
 * invariant under t -> -t, so it does not care which way round the tracer went.
 *
 * STAR-SHAPEDNESS IS A HYPOTHESIS AND IS REPORTED, NEVER ASSUMED. The
 * denominator u' . t = u x t vanishes exactly when the ray is TANGENT to the
 * curve, which is where a ray parametrisation stops being one.
 * AngleParametrisation::transversality is min | u x t | over the fit and
 * transverse is whether it cleared the floor, mirroring
 * IndexAudit::transversality and IndexAudit::transverse and for the same
 * reason. INVERSION-PLAN.md section 3.4 records that ray methods fail on
 * indented cross-sections; this is that hypothesis made measurable. The
 * parametrisation REFUSES rather than returning a number when it degenerates,
 * and the traced curve's own polar angle is checked for monotonicity as a
 * second, independent statement of the same hypothesis -- one that does not
 * depend on the fit succeeding.
 *
 * The rays are seeded and BRACKETED from the traced curve rather than from an
 * assumption: for each target angle the two traced points straddling it give a
 * seed by interpolation and an interval to keep Newton inside, so a Newton step
 * that leaves the bracket falls back to bisection instead of walking to another
 * branch.
 */

namespace meq
{

	/// Which potential the corrector roots. See the header on which won and by
	/// how much; the tracer's constructor takes the fields, so this is a
	/// convenience over the solver rather than a mode.
	enum class Potential
	{
		/// psi_h in W_h, paired with q_h. Needs nothing but a solve.
		Raw,

		/// psi*_h in P_(k+1), paired with the enriched flux q*_h. THE DEFAULT,
		/// on the measurement recorded in the header. Requires
		/// GradShafranovSolver::postProcess() to have been called, and the
		/// constructor refuses rather than rooting an empty field if it has not.
		PostProcessed
	};

	/// Which field answers a call at a point OUTSIDE the mesh. See the header:
	/// the two are a real choice and the choice was made by measurement, with
	/// the table in tests/convergence/FluxSurfaceConvergence.cpp.
	enum class BandExtension
	{
		/// There is no field outside the mesh. A trace that reaches the edge
		/// returns ContourStatus::LeftMesh. THE DEFAULT, and the whole of the
		/// fitted path, where Gamma IS the mesh boundary and there is no band.
		None,

		/// psi( p ) = psi( x0 ) + r0 q( x0 ) . ( p - x0 ) from the foot x0 of p
		/// on Gamma_h, with q frozen at the foot. Nothing is evaluated outside
		/// an element; the extended field is AFFINE, so contours in the band
		/// are straight lines at every k, and the band error caps at O( h^2 ).
		/// It is GridSampler::samplePotentialWithFlux() made pointwise, and it
		/// is kept as the CONTROL rather than as a fallback.
		FluxTaylor,

		/// The extension technique's own construction: the datum on Gamma plus
		/// the line integral of -grad_bar( psi ) = -r q back from it, with q
		/// outside the mesh supplied by mfem::ElementExtension. Needs the same
		/// mfem::TransferPath the solve was given. THE ONE TO USE; see the
		/// header for what it costs and what it buys.
		TransferLift
	};

	/// "none", "flux Taylor", "transfer lift". For printing.
	char const *bandExtensionName( BandExtension which );

	/// How a trace ended. Every one of these is a distinct outcome and a
	/// consumer must not treat the last three as "a shorter curve":
	/// v0-legacy:FluxSurfaces.cpp printed "Terminating because curve left
	/// domain" and returned a partial curve, which for an outer surface is an
	/// arc labelled as a closed contour.
	enum class ContourStatus
	{
		/// A complete circuit, or the requested number of them.
		Closed,

		/// The curve reached the edge of the field. On the fitted path that is
		/// the plasma boundary. On the curved path with BandExtension::None it
		/// is Gamma_h, and the answer is to configure an extension rather than
		/// to accept the arc; with one configured it is the far edge of the
		/// band -- the point sat further outside Gamma_h than
		/// setBandExtension()'s reach allows, which means either the reach is
		/// too small for dist( Gamma_h, Gamma ) or the level is outside Gamma
		/// altogether.
		LeftMesh,

		/// The corrector failed to reach its tolerance at some point, or the
		/// flux vanished so that there is no tangent. The latter is a critical
		/// point and means the level was chosen at or beyond a saddle.
		Stalled,

		/// The point budget ran out before the curve closed. Raise
		/// setMaxPoints(), or the step floor is too small for the contour.
		TooLong
	};

	/// "closed", "left mesh", "stalled", "too long". For printing.
	char const *contourStatusName( ContourStatus status );

	/// One accepted point of a traced contour. Everything in it is either free
	/// -- the tracer evaluated it anyway -- or is provenance the caller cannot
	/// reconstruct afterwards.
	struct ContourPoint
	{
		double r = 0.0;
		double z = 0.0;

		/// The unit tangent from q: ( -q_z, +q_r ) / |q|. Free, and it is what
		/// makes the cubic Hermite of Contour::pointOnSegment() fourth order.
		double tangentR = 0.0;
		double tangentZ = 0.0;

		/// | q | at the point, in the units of q. Free for the same reason, and
		/// it is what IN-2's flux-surface averages weight by: dl / |grad psi|
		/// is dl / ( r |q| ).
		double fluxMagnitude = 0.0;

		/// Accumulated POLYLINE length from the first point. Not the arc length
		/// of the curve -- that is Contour::hermiteLength(), which is the same
		/// quantity computed to fourth order rather than second.
		double arcLength = 0.0;

		/// The element the point was corrected in. Consecutive points with
		/// different elements are a face crossing; see Contour::faceCrossings.
		int element = -1;

		/// Corrector iterations spent on this point, and | psi_h( x ) - c |
		/// where it stopped. The second is measurement ( b ) of the header, per
		/// point, and the property worth checking is that it does not grow with
		/// the index.
		int correctorIterations = 0;
		double residual = 0.0;

		/// WHETHER THIS POINT IS IN THE BAND: the field was answered by an
		/// extension outside the mesh rather than read from the element the
		/// point sits in, because it sits in no element. False on the fitted
		/// path and false for every point of a contour that stays inside
		/// Omega_h.
		///
		/// PER POINT, AND A COUNT WOULD NOT DO. A consumer computing an error
		/// norm, or differencing two runs, has to be able to drop these; and a
		/// flux-surface quantity taken over a surface with any of them is known
		/// only to the extension's order, not to the discretisation's. See the
		/// header, and Contour::crossesBand().
		bool extended = false;

		/// For an extended point, how far outside Gamma_h it sits -- the
		/// distance to its foot. Zero for every point inside the mesh. It is
		/// the length the extension's Taylor remainder is taken over, so it is
		/// the number to read a band error against.
		double bandDepth = 0.0;
	};

	/// A traced contour: the points, how the trace ended, and what it cost.
	///
	/// The Hermite evaluators are on the struct rather than on the tracer
	/// because the interpolant is a property of the accepted points and their
	/// tangents alone. A test can therefore build a Contour by hand from an
	/// exact curve and measure the interpolation in isolation, which is what
	/// separates "the Hermite is wrong" from "the trace is wrong".
	struct Contour
	{
		/// The level c that psi was rooted at.
		double level = 0.0;

		ContourStatus status = ContourStatus::Stalled;

		/// The accepted points, in trace order. For a closed contour the last
		/// point IS the first point, exactly -- see the header on why the final
		/// step is shortened rather than left as a stub.
		std::vector<ContourPoint> points;

		/// Accumulated turning of the unit tangent, in radians, signed. A
		/// closed circuit is +/- 2 pi per circuit and the sign says which way
		/// round the trace went.
		double turning = 0.0;

		/// Completed circuits. More than one is how the "( b ) does not
		/// accumulate" property is measured: trace ten circuits and the
		/// residuals and the closure error must not know about it.
		int circuits = 0;

		/// The gap between the final corrected point and the start point,
		/// resolved along the start point's normal and tangent, and measured
		/// BEFORE the curve is closed onto the start point. See the header:
		/// the normal component is the one that says the corrector is doing its
		/// job, and the tangential one is why the final step is shortened.
		double closureNormal = 0.0;
		double closureTangential = 0.0;

		/// Consecutive accepted points in different elements, and the worst and
		/// mean disagreement of psi_h across the face between them, measured at
		/// the crossing point from both sides. Zero unless
		/// ContourTracer::setMeasureFaceJumps() is on, which it is by default.
		/// THE JUMP IS ERROR ( a ), not a tracer defect; see the header.
		int faceCrossings = 0;
		double worstFaceJump = 0.0;
		double meanFaceJump = 0.0;

		/// The ABSOLUTE residual the corrector was asked to reach:
		/// setTolerance() times the scale of psi_h over the mesh. Every
		/// ContourPoint::residual is to be read against it, and a point whose
		/// residual exceeds it is one of the stalledCorrections below.
		double correctorTarget = 0.0;

		/// Points accepted at the best residual the corrector could reach rather
		/// than at its tolerance. That happens where a point lands within the DG
		/// jump of a face, where { psi_h = c } is genuinely ambiguous and no
		/// tolerance tighter than the jump is attainable; ContourPoint::residual
		/// says how close each got and worstResidual is the worst of them. It is
		/// rare and it becomes commoner as Delta_s falls, simply because more
		/// points are placed. A LARGE count is a different matter and means the
		/// step is far below the mesh size.
		int stalledCorrections = 0;

		/// Times the element walk gave up and fell back on Mesh::FindPoints.
		/// Should be zero or near it; a climbing count says the walk is failing
		/// and the trace is going quadratic in the mesh.
		int fallbackLocations = 0;

		/// Accepted points whose field came from the band extension, and the
		/// deepest any of them sat outside Gamma_h. Zero on the fitted path and
		/// zero for a surface wholly inside Omega_h.
		///
		/// THE COUNT IS A SUMMARY OF ContourPoint::extended AND NOT A
		/// SUBSTITUTE FOR IT. It says how much of the curve is in the band; the
		/// per-point flag says which parts, which is what a consumer dropping
		/// them needs.
		int extendedPoints = 0;
		double deepestBandPoint = 0.0;

		/// Which extension answered the calls, so that a stored contour says
		/// what its band points were computed with rather than leaving it to be
		/// remembered.
		BandExtension bandExtension = BandExtension::None;

		/// Corrector cost, and the worst residual over the whole curve. The
		/// last is measurement ( b ) as one number.
		int correctorIterationsTotal = 0;
		int worstCorrectorIterations = 0;
		double worstResidual = 0.0;

		/// The step the tracer actually used, and the extremes it took where
		/// curvature control was on.
		double nominalStep = 0.0;
		double shortestStep = 0.0;
		double longestStep = 0.0;

		bool closed() const
		{
			return status == ContourStatus::Closed;
		}

		/// Whether any of the curve lies outside the mesh. INVERSION-PLAN.md
		/// section 4.3 requires that every flux-surface quantity computed on
		/// such a surface be flaggable as such, so that IN-2 reports the two
		/// populations separately rather than quoting one rate over both --
		/// which would hide exactly the thing that stage exists to measure.
		bool crossesBand() const
		{
			return extendedPoints > 0;
		}

		/// Segments between accepted points. One fewer than the points.
		std::size_t segments() const;

		/// The polyline length: points.back().arcLength. Second order in the
		/// spacing, and kept because it is what the chord control measures.
		double length() const;

		/// The arc length of the cubic Hermite interpolant, by Gauss-Legendre
		/// per segment. Fourth order in the spacing, and the independent route
		/// to the same quantity AngleParametrisation::length() computes from
		/// the metric.
		double hermiteLength( int gaussPoints = 8 ) const;

		/// The cubic Hermite point on segment @a i at parameter @a t in [0,1],
		/// with the tangents from q scaled by the chord length. This is the
		/// representation: it, and not the polyline, is what the contour IS.
		void pointOnSegment( std::size_t i, double t, double &r, double &z ) const;

		/// dx/dt of the same, which is the tangent up to a positive factor.
		void tangentOnSegment( std::size_t i, double t, double &r, double &z ) const;

		/// The straight chord on segment @a i. THE CONTROL, kept so that the
		/// fourth-order claim above is measured against something rather than
		/// asserted.
		void chordOnSegment( std::size_t i, double t, double &r, double &z ) const;

		/// The Hermite point at accumulated polyline length @a s, wrapped into
		/// [ 0, length() ]. The segment parameter is taken linearly in the
		/// polyline length, which is not the arc-length parametrisation and does
		/// not claim to be -- it is a way of asking for a point, not a
		/// reparametrisation.
		void pointAtArcLength( double s, double &r, double &z ) const;
	};

	/// Declared here and defined below: ContourTracer::fitByAngle() returns one,
	/// and the fit needs the tracer's own seam to read psi and q.
	struct AngleParametrisation;

	/**
	 * The predictor-corrector tracer.
	 *
	 * Borrows the potential and the flux; both must outlive it, and the solver
	 * they came from must have been solved. Nothing is cached between traces
	 * except the mesh's element-to-element table, which is the mesh's own, and
	 * -- once setBandExtension() has been called -- the faces of Gamma_h. THOSE
	 * ARE A SNAPSHOT: a mesh refined under the tracer, which is what an adaptive
	 * cycle does, leaves them naming faces that no longer exist. Build a new
	 * tracer after a refinement, as a caller borrowing the solver's fields has
	 * to anyway.
	 *
	 * const throughout: a trace reads the fields and returns a Contour, and the
	 * counters that would otherwise be state live on the Contour instead.
	 */
	class ContourTracer
	{
		public:
			/// The ordinary way in. Potential::PostProcessed -- the default --
			/// takes solver.postProcessedPotential() and
			/// solver.postProcessedFlux() and requires postProcess() to have
			/// been called; Potential::Raw takes solver.potential() and
			/// solver.flux() and needs nothing but a solve.
			/// @throws std::invalid_argument if @a which is PostProcessed and
			///         GradShafranovSolver::postProcess() has not been called.
			explicit ContourTracer( GradShafranovSolver const &solverIn,
			                        Potential which = Potential::PostProcessed );

			/// The same over bare fields, so that the tracer can be pointed at
			/// an interpolated exact psi and q -- which is how a test separates
			/// "the tracer is wrong" from "the discretisation is coarse".
			///
			/// @param potentialIn the scalar field to root.
			/// @param fluxIn      q in MEQ's sign convention, vdim 2, on the same
			///                    mesh. The tangent and the corrector direction
			///                    both come from it; handing it -q traces the
			///                    same curve backwards, which is harmless, and
			///                    handing it grad_bar( psi ) rather than q
			///                    changes the corrector's scaling by r and
			///                    nothing else.
			ContourTracer( mfem::GridFunction const &potentialIn,
			               mfem::GridFunction const &fluxIn );

			/// Trace the connected component of { psi_h = level } through the
			/// point nearest ( startR, startZ ). The start point is corrected
			/// onto the level set before the first step, so it need only be
			/// close.
			///
			/// @throws std::runtime_error if the start point is not in the mesh,
			///         or if the corrector cannot reach the level from it --
			///         which usually means the level is not attained near there.
			Contour trace( double level, double startR, double startZ ) const;

			/// Trace outward from a known axis: walk out along a ray until the
			/// potential brackets @a level, bisect, and trace from there. This
			/// is the ordinary entry point, IN-A's findAxis() being the
			/// prerequisite INVERSION-PLAN.md section 5 says it is.
			///
			/// EIGHT RAYS RATHER THAN ONE, AND THAT IS A MEASUREMENT. +r alone
			/// is the obvious choice and fails on an ordinary benchmark; the
			/// implementation records which one and why. It is not an assumption
			/// of star-shapedness -- one crossing on one ray is all it needs,
			/// and it is a seed rather than a parametrisation.
			///
			/// @throws std::runtime_error if the level is not bracketed on any
			///         ray before the mesh runs out.
			Contour traceFromAxis( double level, CriticalPoint const &axis ) const;

			/// psi_h and q at a physical point, located by the same walk the
			/// tracer uses. Returns false if the point is not in the mesh.
			///
			/// Public because a caller measuring a contour -- which is what
			/// tests/convergence/FluxSurfaceConvergence.cpp does at every
			/// segment midpoint -- needs the same field the tracer rooted,
			/// evaluated the same way. @a hint is IN AND OUT: pass -1 the first
			/// time and the element it leaves behind on every call after, or
			/// every call falls back on Mesh::FindPoints and the measurement
			/// goes quadratic in the mesh.
			bool sampleAt( double r, double z, double &psi, double &qR,
			               double &qZ, int &hint ) const;

			/// The same for a caller with no hint to give.
			bool sampleAt( double r, double z, double &psi, double &qR,
			               double &qZ ) const;

			/// The same, and it also says whether the value came from the band
			/// extension outside the mesh rather than from an element.
			///
			/// A SEPARATE OVERLOAD RATHER THAN AN EXTRA DEFAULTED ARGUMENT,
			/// because the flag is the thing INVERSION-PLAN.md section 4.3
			/// insists a consumer must be able to read PER POINT, and a defaulted
			/// out-parameter is easy to not notice. Without it a caller placing
			/// its own quadrature points -- which is what
			/// meq::surfaceAverages( tracer, contour ) does -- cannot tell
			/// whether a point it sampled is band data, and has to mark whole
			/// segments conservatively from their endpoints instead.
			bool sampleAt( double r, double z, double &psi, double &qR,
			               double &qZ, int &hint, bool &extended ) const;

			/**
			 * Re-parametrise a closed @a contour by poloidal angle about
			 * @a axis, at @a count equispaced angles, and build the metric.
			 * INVERSION-PLAN.md stage IN-1; see AngleParametrisation.
			 *
			 * A member rather than a free function because every ray Newton
			 * reads psi and q through the same seam the tracer does, and threads
			 * the same element hint through it -- a free function would have to
			 * locate each sample from scratch, which is the Mesh::FindPoints
			 * trap this file exists partly to avoid.
			 *
			 * @throws std::invalid_argument if the contour is not closed or
			 *         @a count is below four.
			 * @throws std::runtime_error if the transversality falls below
			 *         @a floor, if the traced angle is not monotone about the
			 *         axis, or if a ray fails to find the level within its
			 *         bracket. Refusing is the point: INVERSION-PLAN.md section
			 *         3.4 records that ray methods fail on indented
			 *         cross-sections, and a number returned from a degenerate
			 *         fit is worse than no number.
			 */
			AngleParametrisation fitByAngle( Contour const &contour,
			                                 CriticalPoint const &axis,
			                                 std::size_t count,
			                                 double floor = 1.0e-6 ) const;

			/// The nominal step, which is the MAXIMUM step where curvature
			/// control is on. Default 0, meaning "half the size of the element
			/// the trace starts in".
			void setStep( double stepIn );

			/// Radians of turning per step. Default 0.15, about forty points to
			/// a circuit on a circle. ZERO DISABLES CURVATURE CONTROL, which is
			/// what a convergence study in Delta_s wants: the step is then
			/// exactly setStep() everywhere except the shortened final one.
			void setTargetTurn( double turnIn );

			/// The step floor, as a fraction of the nominal step. Default 0.02.
			/// It stops a momentary large curvature -- which is a DIFFERENCE and
			/// can be wrong -- from stalling the trace.
			void setMinStepFraction( double fractionIn );

			/// The step ceiling, as a fraction of the current element's size.
			/// Default 1.0. Zero or less disables it, which a fixed-Delta_s
			/// study wants, since otherwise the mesh silently overrides the step
			/// the study is sweeping.
			void setLocalStepCeiling( double fractionIn );

			/// The corrector stops when | psi_h( x ) - c | falls below this
			/// times the scale of psi_h over the mesh. Default 1e-12: the
			/// corrector is Newton on a polynomial along a line and is
			/// quadratic, so this is reached in two or three steps or not at
			/// all.
			void setTolerance( double toleranceIn );

			/// Default 30. A cap, not a target.
			void setMaxCorrectorIterations( int maxIterationsIn );

			/// Default 200000. A trace that reaches it returns
			/// ContourStatus::TooLong rather than a curve.
			void setMaxPoints( int maxPointsIn );

			/// How many circuits to trace before closing. Default 1. More than
			/// one is how "( b ) does not accumulate" is measured, and the
			/// points of the extra circuits are kept rather than discarded.
			void setCircuits( int circuitsIn );

			/// Whether to bisect for the crossing point at each face and
			/// measure the jump of psi_h there. Default true. It costs about
			/// sixty extra locations per crossing and nothing else, and it is on
			/// by default because the jump is the floor under the representation
			/// error and a reader who does not know it is there will read the
			/// floor as a defect in the interpolant.
			void setMeasureFaceJumps( bool measureIn );

			/// How many rings of face neighbours the walk widens over before
			/// falling back on Mesh::FindPoints. Default 12.
			///
			/// IT WAS 4, AND 4 IS ENOUGH FOR trace() AND NOT FOR ANYTHING ELSE.
			/// A trace steps a fraction of an element at a time, so the next
			/// point is almost always in the element it just left or a face
			/// neighbour of it, and Contour::fallbackLocations reads zero at
			/// either depth. The RAYS of an AngleParametrisation do not: they
			/// are placed by angle rather than by arc length, so consecutive
			/// nodes sit one to two cells apart, and four rings of FACE
			/// neighbours reach about four triangles in a straight line and
			/// fewer diagonally. Measured in tests/performance, depth 4 took the
			/// FindPoints fallback on 183 of 576 rays and depth 12 on none.
			///
			/// AND THE FALLBACK IS 73% OF THE COST OF SAMPLING. CLAUDE.md
			/// records Mesh::FindPoints as O( elements x points ) -- it is a
			/// brute-force scan over element centres -- so a miss is enormously
			/// more expensive than the walk it replaces. Raising the depth is
			/// worth about 2x on fitByAngle and surfaceAverages and about 1.57x
			/// on a whole extraction, and the answers are BIT-IDENTICAL at every
			/// depth, which is what makes the change free rather than a
			/// trade: the walk decides how a point is FOUND and not where it is.
			///
			/// It is also what makes threading possible at all. FindPoints is
			/// not reentrant -- it loops over every element through the mesh's
			/// own shared ElementTransformation -- so a shared tracer aborts the
			/// moment any thread takes the fallback. See CLAUDE.md's Traps.
			///
			/// A deeper walk costs more only when it is going to fail anyway,
			/// and it is bounded by the rings it visits where FindPoints is
			/// bounded by the whole mesh, so the deeper default is also the
			/// safer one as the mesh grows.
			void setWalkDepth( int depthIn );

			/**
			 * Answer calls in the band between Gamma_h and Gamma with an
			 * extension, so that a contour crossing it COMPLETES rather than
			 * stopping at the edge of the mesh. INVERSION-PLAN.md section 4.3.
			 *
			 * @param which          which extension. BandExtension::None
			 *                       restores the fitted behaviour and ignores
			 *                       every other argument.
			 * @param gammaHMarkerIn the boundary attributes of Gamma_h, sized by
			 *                       the largest boundary attribute of the mesh --
			 *                       the same array setExtension() was given. An
			 *                       attribute NOT marked is FITTED, and a point
			 *                       outside across one of its faces is outside
			 *                       the plasma and gets no extension, which is
			 *                       why this is not simply "the mesh boundary".
			 * @param pathIn         the transferring paths, for TransferLift and
			 *                       for nothing else. Borrowed, and it must be
			 *                       the family the solve used: a different one is
			 *                       a different lifting. Must have been built on
			 *                       THIS mesh, since it indexes its own faces.
			 * @param g              the datum on the TRUE boundary, as a
			 *                       function of position. Defaults to zero, which
			 *                       is MEQ's fixed-boundary problem -- Gamma IS
			 *                       the level set psi = 0.
			 * @param reach          how far outside Gamma_h a point may sit and
			 *                       still be answered, as a multiple of its
			 *                       nearest face's own length. Default 2.
			 *
			 *                       TWO IS SLACK AND THE SLACK WAS MEASURED
			 *                       RATHER THAN CHOSEN.
			 *                       GridSampler::extendOutward() defaults near
			 *                       one because it is filling a picture; here the
			 *                       band has to be CROSSED, and a reach too small
			 *                       truncates the outermost surfaces on the very
			 *                       faces where the band is widest. Measured on
			 *                       the benchmark by sweeping it: the deepest
			 *                       point of Gamma itself sits at 0.92 to 0.98 of
			 *                       a face length outside Gamma_h, so one works
			 *                       with a margin of two per cent and 0.75
			 *                       truncates at every mesh. Two is a factor of
			 *                       two over what is needed there.
			 *
			 *                       AND A REACH TOO SMALL IS LOUD. The trace ends
			 *                       as ContourStatus::LeftMesh after a handful of
			 *                       points rather than returning a short arc
			 *                       labelled as a closed contour, which is the
			 *                       failure the prior art had.
			 * @param lineOrderIn    order of the quadrature along a path;
			 *                       negative takes twice the flux order plus two,
			 *                       which is HDGExtensionIntegrator's own
			 *                       default. Pass what setExtension() was given.
			 *
			 * @throws std::invalid_argument if TransferLift is asked for without
			 *         a path, if the marker is the wrong size, if it selects no
			 *         attribute, or if @a reach is not positive.
			 * @throws std::runtime_error if the marked attributes carry no
			 *         boundary face -- which means the marker names the wrong
			 *         attribute, and silently extending nothing would read as
			 *         the band being absent.
			 */
			void setBandExtension( BandExtension which,
			                       mfem::Array<int> const &gammaHMarkerIn,
			                       mfem::TransferPath const *pathIn = nullptr,
			                       mfem::PositionFunction g
			                           = mfem::PositionFunction(),
			                       double reach = 2.0,
			                       int lineOrderIn = -1 );

			/// Forget it; a trace then stops at the edge of the mesh again.
			void clearBandExtension();

			/// Which extension is configured, BandExtension::None if none is.
			BandExtension bandExtension() const;

			/// How many faces of Gamma_h the extension was built over. Zero
			/// unless setBandExtension() has been called, and a caller checking
			/// that the marker named the right attribute wants this rather than
			/// the mesh's own boundary count.
			std::size_t bandFaceCount() const;

		private:
			/// trace(), given the element the start point is in. traceFromAxis()
			/// already knows it from its own bracket search, so this is what
			/// keeps a trace free of Mesh::FindPoints entirely.
			Contour traceFrom( double level, double startR, double startZ,
			                   int seed ) const;

			/// THE SEAM. The only place psi and q are read at a physical point.
			/// See the header: the band extension of IN-0's second half is a
			/// change to this function and to nothing else.
			struct FieldSample
			{
				int element = -1;
				mfem::IntegrationPoint ip;
				double psi = 0.0;
				double qR = 0.0;
				double qZ = 0.0;

				/// Whether the value came from an extension outside the mesh
				/// rather than from an element. See the header: it is carried
				/// through to ContourPoint::extended, because a count is not
				/// enough and a consumer has to be able to tell WHICH points are
				/// in the band.
				bool extended = false;

				/// For an extended sample, how far outside Gamma_h the point sat.
				/// Zero otherwise, and it is the length the extension's Taylor
				/// remainder is taken over -- the number a band error is read
				/// against.
				double depth = 0.0;

				/// For an extended sample the element above is the one whose
				/// polynomial was extended and it does NOT contain the point; this
				/// is the foot on its face. Kept because the Taylor step's own r is
				/// the FOOT's radius and not the point's, and getting that wrong is
				/// the factor of 1.7e5 CLAUDE.md records against sampleCoefficient().
				double footR = 0.0;
				double footZ = 0.0;
			};

			/// Locate ( r, z ) and evaluate both fields there. @a hint is the
			/// element the previous call used, or -1. Returns false when the point
			/// is outside the mesh and no band extension answers for it, which on
			/// the fitted path means outside the plasma.
			bool sampleField( double r, double z, int hint, FieldSample &sample,
			                  int &fallbacks ) const;

			/// The element walk: TransformBack into @a hint, then its face
			/// neighbours, then a second ring, then Mesh::FindPoints. Increments
			/// @a fallbacks when it reaches the last resort.
			///
			/// @a allowFallback EXISTS FOR THE BAND AND FOR NOTHING ELSE. With an
			/// extension configured most failed walks are band points, and paying an
			/// O( elements ) FindPoints scan for each of them would make an outer
			/// contour quadratic in the mesh -- which is the cost this walk was
			/// written to avoid. sampleField() therefore tries the walk without the
			/// last resort, then the band, and only then the last resort, so a point
			/// that is genuinely inside and merely lost is still found and still
			/// counted.
			bool locate( double r, double z, int hint, int &element,
			             mfem::IntegrationPoint &ip, int &fallbacks,
			             bool allowFallback = true ) const;

			/// One face of Gamma_h, as the band search needs it. The endpoints are
			/// in the FACE transformation's own parametrisation -- xi = 0 at
			/// ( r0, z0 ) -- because that is what mfem::TransferPath::Endpoint
			/// interpolates its vertex directions along, and the mesh's own vertex
			/// order is not required to agree with it.
			struct BoundaryFace
			{
				double r0 = 0.0, z0 = 0.0;
				double r1 = 0.0, z1 = 0.0;

				/// The OUTWARD unit normal, oriented away from the element's
				/// centroid rather than read off CalcOrtho: one fewer convention to
				/// be wrong about, and it costs three subtractions.
				double normalR = 0.0, normalZ = 0.0;

				double length = 0.0;
				int element = -1;
				int boundaryElement = -1;
			};

			/// The nearest face of Gamma_h to ( r, z ), its foot, the face parameter
			/// of the foot and the distance to it. Returns -1 if no face is within
			/// reach, or if ( r, z ) is on the INSIDE of the nearest one -- in which
			/// case the point is not in the band at all and the caller still owes it
			/// a proper location.
			int nearestBandFace( double r, double z, double &footR, double &footZ,
			                     double &parameter, double &depth ) const;

			/// Answer psi and q at a point outside the mesh. Returns false if the
			/// point is not in the band. Fills @a sample with extended = true and
			/// with the element whose polynomial was extended -- which does NOT
			/// contain the point.
			bool extendField( double r, double z, FieldSample &sample ) const;

			/// The minimum-norm Newton corrector: x <- x + grad( psi )( c - psi )
			/// / | grad( psi ) |^2, iterated. Fills @a sample with the accepted
			/// point. Returns false if it left the mesh, ran out of iterations,
			/// or met a vanishing flux.
			/// @a target is the absolute residual to stop at, which the callers
			/// compute once from potentialScale(); @a maxMove is how far the
			/// iterate may travel from where the predictor put it before the
			/// correction is refused. @a stalled says the point was accepted at
			/// the best residual reached rather than at the tolerance, which is
			/// what happens where a point lands within the DG jump of a face.
			/// See the implementation: both are consequences of psi_h being
			/// discontinuous and neither is a defect of the iteration.
			bool correct( double level, double target, double maxMove, double &r,
			              double &z, int hint, FieldSample &sample,
			              int &iterations, bool &stalled, int &fallbacks ) const;

			/// sqrt of the element's own area measure, as the length scale the
			/// step ceiling is a fraction of.
			double elementSize( int element ) const;

			/// The scale | psi_h | is measured against by the corrector's
			/// stopping rule: the largest | psi_h | over the potential's dofs.
			double potentialScale() const;

			/// Bisect along the segment for the point where the located element
			/// changes, and return the disagreement of psi_h across it. Zero if
			/// the two elements turn out to be the same after all.
			double faceJump( double r0, double z0, int element0,
			                 double r1, double z1, int element1,
			                 int &fallbacks ) const;

			mfem::GridFunction const &potentialField;
			mfem::GridFunction const &fluxField;
			mfem::Mesh &meshRef;
			mfem::Table const &neighbours;

			double step = 0.0;
			double targetTurn = 0.15;
			double minStepFraction = 0.02;
			double localStepCeiling = 1.0;
			double tolerance = 1.0e-12;
			int maxCorrectorIterations = 30;
			int maxPoints = 200000;
			int circuits = 1;
			int walkDepth = 12;
			bool measureFaceJumps = true;

			/// The band. Empty and None on the fitted path, which is the default
			/// and costs a trace there one comparison per failed walk.
			BandExtension bandMethod = BandExtension::None;
			std::vector<BoundaryFace> bandFaces;
			mfem::TransferPath const *bandPath = nullptr;
			mfem::PositionFunction bandDatum;
			double bandReach = 2.0;
			int bandLineOrder = -1;
	};

	/**
	 * A closed contour re-parametrised by poloidal angle about the axis, and
	 * the metric that comes with it. INVERSION-PLAN.md stage IN-1.
	 *
	 * The points are at equispaced theta and each is found by a 1-D Newton along
	 * its ray, using q for the derivative -- no differencing anywhere in the
	 * search. rho' comes from the identity in the header, pointwise from q, and
	 * | dx/dtheta | = sqrt( rho'^2 + rho^2 ) follows because u and u' are
	 * orthonormal.
	 *
	 * THE POINT OF THE CLASS IS THAT IT ALSO COMPUTES THE WRONG ANSWER. Three
	 * lengths are returned from one fit:
	 *
	 *   length()             periodic trapezoid, rho' from q          the answer
	 *   differencedLength()  periodic trapezoid, rho' by central       THE TRAP
	 *                        difference of the rho_j
	 *   chordLength()        sum of chords                             the naive
	 *                                                                  control
	 *
	 * The first is the rule the plan wants and the second is the same rule fed a
	 * second-order Jacobian -- a spectrally accurate quadrature reduced to
	 * second order by its metric, with nothing in its own output saying so. If
	 * either control ever converges as fast as the first the comparison is empty
	 * and the test built on it is worthless, which is what
	 * tests/convergence/FluxSurfaceConvergence.cpp says in its failure message.
	 *
	 * THE FIRST COLUMN IS NOT SPECTRAL ON A DISCRETE CONTOUR AND THE PLAN SAID
	 * IT WOULD BE. It cannot be: psi_h is discontinuous across faces, so
	 * rho( theta ) for a level set of it is piecewise analytic with a jump at
	 * every face crossing, and no quadrature rule is geometric on a function
	 * with jumps. What is measured is that it falls very much faster than
	 * second order and then FLOORS at the jump -- and the floor is quantitative:
	 * on the benchmark the DG jump of psi* converts to a distance of 6.8e-10
	 * and the length error stops at about 1.2e-9. The control that says the
	 * floor is the FIELD and not the RULE is
	 * theMetricIdentityIsSpectralOnTheClosedForm, which runs the identical rule
	 * on the analytic contour, where rho really is analytic, and reaches 4e-15.
	 *
	 * Which does not weaken the comparison at all: five orders of magnitude
	 * separate the first column from the second by the end of the sweep, and
	 * every one of them is the Jacobian.
	 */
	struct AngleParametrisation
	{
		/// The axis the rays are drawn from.
		double axisR = 0.0;
		double axisZ = 0.0;

		/// The level, carried through from the contour.
		double level = 0.0;

		/// theta_j = 2 pi j / N, j = 0 .. N-1. Not stored: it is implied by the
		/// index and by size().
		std::size_t count() const
		{
			return radius.size();
		}

		/// rho( theta_j ), and the point it puts on the curve.
		std::vector<double> radius;
		std::vector<double> pointR;
		std::vector<double> pointZ;

		/// q at the node, in MEQ's sign convention -- the very flux the ray
		/// Newton used for its derivative and the identity below used for rho'.
		///
		/// STORED BECAUSE IT WAS COMPUTED AND THROWN AWAY. Every node of a fit
		/// evaluates q to build rho' and, before 2026-09-03, discarded it, so
		/// meq::surfaceAverages() re-read the field at every node through
		/// sampleAt() to get the same numbers back. That is one avoidable
		/// evaluation per node in the one place CLAUDE.md records
		/// Mesh::FindPoints as O( elements x points ). Keeping it costs two
		/// doubles a node.
		std::vector<double> fluxR;
		std::vector<double> fluxZ;

		/// psi_h at each node, from the same evaluation that placed it.
		///
		/// STORED FOR THE SAME REASON AS fluxR AND fluxZ, AND THE TWO ARE ONLY
		/// USEFUL TOGETHER. With the potential and the flux both here, a
		/// consumer integrating over these nodes -- meq::surfaceAverages() is
		/// the one this exists for -- needs no field evaluation of its own at
		/// all. Re-sampling would cost one location per node in the one place
		/// CLAUDE.md records Mesh::FindPoints as O( elements x points ), and it
		/// would risk landing on the other side of a face from the fit, which
		/// on a discontinuous field is a different number.
		///
		/// NOT THE SURFACE'S LEVEL, AND THE DIFFERENCE IS THE POINT. On a
		/// stalled ray the accepted node sits as close to the level as the
		/// field's own jump allows and no closer -- see stalledRays -- so this
		/// records where the node actually is rather than where it was asked to
		/// be.
		std::vector<double> potential;

		/// rho'( theta_j ), POINTWISE FROM q by the identity in the header.
		std::vector<double> radiusPrime;

		/// | dx/dtheta |_j = sqrt( rho'^2 + rho^2 ), the arc-length element.
		std::vector<double> speed;

		/// | u x t |_j at each node. The denominator of the identity, so this is
		/// how close the ray came to being tangent to the curve.
		std::vector<double> crossing;

		/// 1 where the node's field came from the band extension rather than
		/// from an element, 0 otherwise. THE SAME OBLIGATION ContourPoint
		/// carries, one level up: a fit is a set of points and the band does not
		/// respect its angular grid, so a surface can be inside Omega_h at one
		/// theta and outside it at the next. IN-2 integrates over these nodes,
		/// so it needs to know which of them the extension answered for -- a
		/// count over the whole surface would say a quantity is affected and not
		/// where.
		///
		/// char rather than bool because std::vector<bool> is a bitset whose
		/// elements have no address, and this is data a consumer will want to
		/// pass around like the columns beside it.
		std::vector<char> extended;

		/// How many of those there are, and the deepest. Zero on the fitted path
		/// and zero for a surface wholly inside Omega_h.
		int extendedNodes = 0;
		double deepestBandNode = 0.0;

		/// Whether any node is in the band. The predicate IN-2 reports the two
		/// populations by.
		bool crossesBand() const
		{
			return extendedNodes > 0;
		}

		/// min over j of crossing, and whether it cleared the floor. See the
		/// header: star-shapedness is a hypothesis and this is it, measured.
		double transversality = 0.0;
		bool transverse = false;

		/// Whether the polar angle of the TRACED curve was monotone about the
		/// axis, and the worst backward step in radians if it was not. An
		/// independent statement of the same hypothesis, and one that does not
		/// depend on the fit having succeeded.
		bool angleMonotone = true;
		double worstBacktrack = 0.0;

		/// The worst | psi_h( x_j ) - c | over the fit, and the worst and total
		/// Newton iterations. The first is measurement ( b ) for the rays.
		double worstResidual = 0.0;
		int worstIterations = 0;
		int totalIterations = 0;

		/// How many rays needed the bisection fallback because Newton left its
		/// bracket. Zero on a well-behaved surface, and a climbing count is the
		/// star-shapedness hypothesis fraying before it fails.
		int bisections = 0;

		/// Rays accepted at the closest point they could reach rather than at
		/// the tolerance. THE EXACT COUNTERPART OF Contour::stalledCorrections,
		/// and it is here for the same reason: { psi_h = c } is a union of
		/// per-element arcs offset by the DG jump, so a ray crossing a face where
		/// c falls inside the jump has NO point on it with psi_h = c, and no
		/// tolerance tighter than the jump can be met. worstResidual says how
		/// close the worst of them got and is to be read against
		/// setTolerance() times the scale of psi_h.
		///
		/// IT IS RARE AND IT GETS COMMONER WITH THE ANGLE COUNT, more rays being
		/// more chances -- the same shape as stalledCorrections getting commoner
		/// as Delta_s falls. A LARGE count is a different matter and says the
		/// tolerance asked for is below the field's own ambiguity everywhere,
		/// which is a statement about the mesh and the degree rather than about
		/// this fit.
		int stalledRays = 0;

		/// Times the element walk fell back on Mesh::FindPoints during the fit.
		/// The same measurement Contour::fallbackLocations is, for the same
		/// reason: it is a performance statement and nothing else would say it.
		int fallbackLocations = 0;

		/// The periodic trapezoid rule on sqrt( rho'^2 + rho^2 ), with rho' from
		/// q. THE ANSWER.
		double length() const;

		/// The same rule with rho' replaced by a central difference of the
		/// rho_j. THE TRAP: second order, and nothing about the rule says so.
		double differencedLength() const;

		/// The sum of the chords between consecutive points. The naive control.
		double chordLength() const;

		/// rho' by Fourier differentiation of the sampled rho_j -- the spectral
		/// derivative of the same periodic sequence. An INDEPENDENT route to
		/// rho', used to check the identity rather than to compute anything: it
		/// needs only the rho_j and knows nothing of q. O( N^2 ), which is
		/// nothing at these sizes.
		std::vector<double> spectralRadiusPrime() const;
	};

}

#endif
