Flux surfaces
=============

A Grad–Shafranov solve produces :math:`\psi(r, z)` on a mesh. Almost everything
downstream — a transport code, a gyrokinetic code, a plot of the safety factor
— wants the *inverse* of that: the flux surfaces themselves, as curves, indexed
by the flux label rather than by position. meq extracts them, and this chapter
is how.

Three things are built on top of one another, and each is usable on its own:

#. **The critical points.** The magnetic axis, and any X-point, as roots of the
   flux :math:`q`. A tracer needs the axis before it can start, and it needs to
   know where the saddles are because that is where a level set stops being a
   curve.
#. **The tracer.** A predictor–corrector walk along :math:`\{\psi = c\}`,
   returning a closed contour with a cubic interpolant between its points.
#. **The parametrisation.** The same contour re-expressed as a radius against
   poloidal angle, with the arc-length element that turns it into a quadrature
   rule.

:doc:`surface_geometry` is what those three are for: flux-surface averages, and
the surfaces as a smooth map from a disc.

.. note::

   **This is library-only.** None of it is reachable from a configuration file
   and none of it appears in the output files described in :doc:`output`. The
   headers are ``meq/CriticalPoints.hpp``, ``meq/FluxSurfaces.hpp``,
   ``meq/SurfaceAverage.hpp``, ``meq/Zernike.hpp`` and ``meq/SurfaceFit.hpp``,
   and none of them is in the ``meq/meq.hpp`` umbrella; include them directly.

.. code-block:: cpp

   #include "meq/CriticalPoints.hpp"
   #include "meq/FluxSurfaces.hpp"

   solver.solve();
   solver.postProcess();

   meq::CriticalPointFinder finder( solver );
   meq::CriticalPoint axis = finder.findAxis();

   meq::ContourTracer tracer( solver );
   meq::Contour surface = tracer.traceFromAxis( level, axis );
   meq::AngleParametrisation fit = tracer.fitByAngle( surface, axis, 128 );

Everything borrows: the finder and the tracer hold the solver's fields by
reference, and the solver must have been solved. See :ref:`api-borrowing`.

.. _flux-surfaces-q:

Why the flux is the asset
-------------------------

Every step below is cheaper and more accurate than it would otherwise be for
one reason, and it is the same reason meq is a mixed method at all: **the
gradient of** :math:`\psi` **is a solved unknown, not a derivative of one**.

Writing :math:`\gradbar\psi = r q`:

* A critical point is where :math:`q = 0`. That is a residual converging at the
  potential's own order, rather than a differentiated potential converging an
  order down.
* The tangent to a contour is :math:`(-q_z, +q_r)/|q|`, available at every
  point the tracer has already visited, at no cost — the tracer evaluated
  :math:`q` for the predictor anyway.
* The weight in every flux-surface average is :math:`1/|\gradbar\psi| = 1/(r
  |q|)`, again pointwise and again at the flux's own order.
* The gauge-free fit of :doc:`surface_geometry` needs
  :math:`\nabla\Psi_{\mathrm N}`, which is the same field once more.

Contrast the continuous Galerkin codes. :cite:t:`Heumann2015` records as an
open problem that in a P1 discretisation the axis and the X-point are confined
to mesh vertices, because that is where the discrete gradient can be said to
vanish. meq resolves both **sub-element**, by root finding inside an element's
own polynomial.

Critical points
---------------

The magnetic axis and any X-point are the zeros of :math:`q_h`, and
`meq::CriticalPointFinder` finds them by Newton inside one element
at a time.

Classification is free and needs no second derivative of :math:`\psi`.
Differentiating :math:`r q = \gradbar\psi` gives

.. math::

   \mathrm{Hess}(\psi) = q \otimes e_r + r \, \frac{\partial q}{\partial x},

whose first term vanishes identically at a point where :math:`q = 0`. So the
Hessian there is :math:`r` times the Jacobian of :math:`q`, and :math:`r > 0`
throughout an axisymmetric domain — the determinant scales by :math:`r^2` and
the trace by :math:`r`, both positive, so the *signs* that decide maximum,
minimum or saddle can be read straight off :math:`\partial q/\partial x`.

.. note::

   The Jacobian is taken by central differences in the element's reference
   coordinates, and that is not a shortcut. The located root is where
   :math:`q_h` vanishes, which is a property of :math:`q_h` alone; the Jacobian
   steers the iteration and appears nowhere in its fixed point. A wrong
   Jacobian costs Newton steps and buys no error. It is the same observation
   :ref:`sources-jacobian` makes about the Grad–Shafranov Newton itself, from
   the other side.

.. important::

   **This is not** :cpp:func:`meq::GradShafranovSolver::psiAxis`, **and the two
   must not be reconciled.**

   The solver's :math:`\psiax` is *the largest nodal value* of
   :math:`\psi_h`. It is deliberately that rather than the maximum of the
   polynomial, because the bordered Newton of :doc:`normalised_flux` needs a
   constraint it can differentiate, and one nodal value is one entry of the
   discrete unknown. The critical point is *the place where* :math:`q_h`
   *vanishes*.

   They differ by :math:`O(h)` in position and :math:`O(h^2)` in value, and
   **both orders are independent of** :math:`k` — where :math:`\psi_h`'s own
   error is :math:`k+1`. So on a refined high-order mesh the two readings
   **separate** rather than converge. Neither is a defective version of the
   other. Use the solver's for the normalisation and the finder's for the
   geometry.

Which extremum is the axis
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   **meq's** :math:`\psi` **is not sign-normalised across sources, so the
   magnetic axis is not always a maximum.** With :math:`F` single-signed
   *negative* — which is what the Solov'ev benchmarks have — :math:`\psi` is a
   subsolution, its maximum is on :math:`\Gamma`, and the axis is an interior
   **minimum**. With :math:`F` positive it is a maximum.

   ``findAxis()`` therefore seeds from *both* nodal extremes and returns
   whichever yields a genuine interior extremum, and **refuses rather than
   guessing** if both are present. ``AxisSense`` is there for a caller who knows
   which they want. A search seeded only from the largest nodal value finds a
   corner of the mesh for one sign of :math:`F`.

.. warning::

   **Hand the finder** ``solver.flux()``, **not the raw block.** The
   library's own convention holds :math:`-q`; in two dimensions the index of a
   vector field is unchanged by negating it, so every winding number below comes
   out the same and **every maximum silently becomes a minimum**. The audit
   still passes. See :ref:`sign-conventions`.

.. _flux-surfaces-audit:

The index audit, and what a degree is
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``audit()`` walks the boundary of the mesh and accumulates the turning of
:math:`q`. That is the topological degree of :math:`q` on the boundary, which by
the Poincaré index theorem equals the **sum of the indices** of the zeros
inside: :math:`+1` for either extremum, :math:`-1` for a saddle. It is a
one-dimensional integral — no subdivision, no root finding, and its cost is the
boundary rather than the domain.

.. danger::

   **A degree is a sum of indices and never a count.** Degree zero does **not**
   mean there is nothing there: a maximum and a saddle inside sum to zero, and
   the boundary cannot tell that from an empty domain. The test suite keeps a
   box drawn round both an axis and an X-point as a live demonstration of
   exactly that. Anything that reads a zero winding number as "no critical
   points here" is wrong.

The audit is a **certification** and never an exclusion. Its value is that the
sum of the indices of whatever ``sweep()`` found can be compared against it, and
a disagreement is then positive evidence that something was missed.

.. note::

   It is also blind to *spurious pairs*. Noise in :math:`q_h` creates critical
   points strictly in twos — a spurious maximum beside a spurious saddle —
   because a small perturbation of a field cannot change its degree, so the pair
   sums to zero and the audit passes with it present. The complementary check is
   a persistence threshold, which is not implemented; this note is where a
   reader finds out that the audit alone does not cover it.

``windingNumber == eulerCharacteristic`` is Poincaré–Hopf only when :math:`q` is
**transverse** to the boundary — pointing consistently outward (or inward) and
vanishing nowhere on it. That is a condition on :math:`q \cdot n` and **not** on
the boundary being a level set, which is easy to assume otherwise. Both
situations occur in practice: the standard benchmark rectangle is a level set of
nothing yet :math:`q\cdot n` keeps one sign all the way round, so the equality is
a theorem there; a box drawn wide enough to enclose an X-point fails the
hypothesis outright, and the degree then reads zero against an Euler
characteristic of one with no contradiction whatever.
``IndexAudit::transversality`` and ``IndexAudit::transverse`` record which
situation the caller is in, so the comparison is not read as a theorem where it
is a coincidence.

.. warning::

   **A search is not exhaustive.** ``findAxis()`` and ``sweep()`` are seeded
   Newton, and Newton certifies the root it converges to while saying nothing
   about the roots it does not. Do not read a ``sweep()`` result as "these are
   all of them". The exhaustive construction — subdivision in a barycentric
   Bernstein basis, discarding sub-triangles whose coefficients provably share
   a sign — is deliberately not built.

.. note::

   A discontinuous :math:`q_h` can carry boundary degree one with no zero
   strictly inside **any** element. Each element's polynomial puts its zero just
   into a neighbour's territory, and with a face jump much smaller than an
   element there is a window where the zero belongs to neither. Poincaré–Hopf is
   a theorem about continuous fields; this is the discontinuous Galerkin jump
   meeting it head on. The window closes with refinement, so it is a property of
   the mesh rather than a defect — and the practical rule is that **where the
   audit and the search disagree, believe the audit**.

   ``CriticalPoint::overshoot`` is the same phenomenon at the level of one root,
   and ``setContainment()`` is the allowance for it, in *reference*-element units
   so that it shrinks with the mesh exactly as the ambiguity it covers does. It
   is not a tuning parameter, and that was checked by sweeping it: the located
   axis does not move.

.. _flux-surfaces-tracer:

Tracing a surface
-----------------

`meq::ContourTracer` follows the connected component of
:math:`\{\psi = c\}` through a starting point. ``traceFromAxis()`` is the
ordinary entry: it walks out along rays from a known axis until the potential
brackets the level, bisects, and traces from there.

Each step is an explicit Euler predictor along the tangent, followed by a
**Newton corrector onto the level set**. That the predictor is only first order
is not a defect, and the reason is worth stating because it settles a natural
objection: the corrector is a root find, so every accepted point sits on the
discrete level set to tolerance *whatever the predictor did*. There is nothing
for a higher-order predictor to buy in accuracy. What it buys is step length.
:cite:t:`AllgowerGeorg2003` is the reference for the whole construction, and
states the same preference from the other side.

What the corrector buys is that error off the contour does not **accumulate**.
Integrate the tangent field without projection and the curve drifts onto a
neighbouring contour and does not close; with the projection each point is
independently on the curve, and there is no accumulation at all. It also makes
the face jump a near-non-issue: a step that crosses a face lands in whichever
element it lands in, and the corrector re-anchors it there. An ODE integrator
would instead be meeting a discontinuity that violates the smoothness every
Runge–Kutta order is derived from.

Three errors, measured separately
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A traced contour is wrong in three independent ways, and conflating them is how
this problem gets done badly.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Error
     - What it is
   * - **Field**
     - How far :math:`\{\psi_h = c\}` is from :math:`\{\psi = c\}`. This is the
       discretisation, divided by :math:`|\gradbar\psi|`. Not the tracer's to
       improve — except through *which potential it roots*, below.
   * - **Point location**
     - How far an accepted point is from :math:`\{\psi_h = c\}`. This is what
       the corrector buys, at its tolerance, and the property that matters is
       that it is **independent of path length**. ``setCircuits()`` traces more
       than one lap for exactly that measurement, and it does not move.
   * - **Representation**
     - How far the interpolant *between* accepted points is from the curve.
       Governed by the point spacing and the interpolation order, and
       **decoupled from** :math:`h` **and** :math:`k` **entirely**.

The last of those is the one a caller controls directly, and it is why the point
spacing is a free parameter rather than something tied to the mesh.

.. _flux-surfaces-potential:

Which potential to root
~~~~~~~~~~~~~~~~~~~~~~~

meq has two candidate pairs: :math:`\psi_h` with :math:`q_h`, both converging at
:math:`k+1`; and the post-processed :math:`\psi^\star` with the enriched
:math:`q^\star`, where :math:`\psi^\star` converges at :math:`k+2` (see
:doc:`postprocessing`). ``Potential::PostProcessed`` is **the default**, and the
constructor refuses a solver on which ``postProcess()`` has not been called.

Two things follow from rooting the better field, and the second was not
expected.

**The contour inherits the potential's order.** Every accepted point is on
:math:`\{\psi^\star = c\}` to tolerance, so by the implicit function theorem the
traced curve is :math:`k+2` from the true one rather than :math:`k+1` —
measured, orders of magnitude closer on the same mesh, and at *fewer* corrector
iterations per point rather than more.

.. important::

   **And it is what makes the cubic interpolant work at all on this
   discretisation.** The interpolant is built on tangents from the *flux* and is
   measured against the level set of the *potential*, and those are the same
   curve only so far as the two fields agree. :math:`\gradbar\psi_h / r` agrees
   with :math:`q_h` only to :math:`O(h^k)` — differentiating a discontinuous
   Galerkin potential of degree :math:`k` loses an order, while :math:`q_h`
   keeps :math:`k+1`. So paired with :math:`\psi_h` the interpolant is fourth
   order in the spacing **until that tangent tilt takes over, and second order
   afterwards**. Paired with :math:`\psi^\star` the tilt is a full order smaller
   and the fourth order survives.

.. warning::

   **The usual reason given for preferring** :math:`\psi^\star` **is wrong**,
   and the difference matters. It is *not* that the local post-processing is
   built so that :math:`\gradbar\psi^\star` matches :math:`r q^\star`: the
   constraint equation projects the **total** flux onto a face restriction of a
   Raviart–Thomas space, and that field is what drives Stenberg's local problem
   :cite:p:`Stenberg1991` — so the two are different objects. What is true, and
   measured, is that they agree an order better than the raw pair does. The
   wrong mechanism predicts exactness, and would have made the measured tilt
   look like a defect.

``Potential::Raw`` is the other door. It needs nothing but a solve, and it is
what a study of the field error has to use, :math:`k+1` and :math:`k+2` being
the two orders such a study is telling apart.

The interpolant, and the step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Between two accepted points the contour **is** a cubic Hermite, with the unit
tangents from :math:`q` scaled by the chord length. That is fourth order in the
spacing, against second order for the straight chords the same points would
otherwise be joined by; ``Contour::pointOnSegment`` is the interpolant and
``Contour::chordOnSegment`` is the control kept beside it.

.. note::

   Chord scaling is not the optimal tangent scaling, and the arithmetic looks as
   though it should cost an order. It does not — the two differ by a relative
   :math:`O(\theta^2)` on an arc of turning angle :math:`\theta`, which leaves
   the fourth order intact.

The step is **curvature controlled**: a constant turning angle per step, so that
the interpolation error per segment is equidistributed. That is a deliberate
departure from the standard continuation strategies, which key the step on
*corrector effort* :cite:p:`AllgowerGeorg2003` — right for software whose object
is to traverse a path cheaply, and wrong for one whose object is a well-shaped
curve. ``setTargetTurn( 0 )`` disables it, which is what a convergence study in
the spacing wants.

The curvature is itself a difference, and that is fine for the same reason the
finder's Jacobian is: it steers the step controller and appears nowhere in any
accepted point. A wrong curvature costs a badly shaped point distribution and
cannot move a single point of the curve. The step is bounded from both sides —
a floor, so a momentary large curvature cannot stall the trace, and a ceiling
that is a fraction of the local element size, so a step cannot leap a whole
element.

Closure
~~~~~~~

.. note::

   **"Closes at machine precision" is not the right statement, and the correct
   one is geometric.** :math:`\psi_h` is discontinuous across faces, so
   :math:`\{\psi_h = c\}` is a union of per-element arcs: the curve is closed as
   a set, but it is not one analytic curve.

   *Tangentially*, the gap is whatever the final step leaves, so the final step
   is **shortened** to aim at the start point rather than taking a full step.
   *Normally*, the gap is then bounded by the curve **bending** over that
   tangential offset, and is independent of path length. Both are reported, as
   ``closureTangential`` and ``closureNormal``, and the second is the one that
   says the corrector is doing its job — a tracer that integrated rather than
   projected would show a normal gap growing with path length and unrelated to
   the tangential one.

   The measurement is taken *before* the curve is closed, and the closing point
   is then set to the start point exactly, so the polyline a consumer receives
   has no stub and no almost-duplicate node.

What the tracer reports, and why
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every counter on a ``Contour`` is either free — the tracer computed it anyway —
or is provenance a caller cannot reconstruct afterwards.

``status``
   ``Closed``, ``LeftMesh``, ``Stalled`` or ``TooLong``. **These are four
   distinct outcomes and the last three are not "a shorter curve".** The prior
   art in this tree printed a message and returned a partial curve, which for an
   outer surface is an arc labelled as a closed contour — and a flux-surface
   average over it is wrong by an amount nothing reports.

``faceCrossings``, ``worstFaceJump``
   Where two consecutive points are in different elements, the segment crosses a
   face and :math:`\psi_h` disagrees with itself there by the jump. **That is
   the field error, not a tracer defect** — it converges at the order of
   whichever field is being rooted. It is reported because it is the *floor*
   under the representation error: a segment whose endpoints are corrected onto
   two different arcs cannot be closer to either than the arcs are to each
   other, so refining the spacing on a fixed mesh drives the interpolation error
   down at fourth order until it meets that floor and then stops. A reader who
   does not know the floor is there will read it as a defect in the interpolant.

``stalledCorrections``, ``worstResidual``
   .. warning::

      **A point can land where no tolerance is attainable, and it gets commoner
      as the spacing falls.** Because :math:`\{\psi_h = c\}` is a union of arcs
      offset by the face jump, a point landing within the jump of a face is on
      **neither**: the Newton step computed in element A pushes it into B, B's
      pushes it back, and the residual alternates without ever meeting a
      tolerance tighter than the jump. It is rare, it becomes *more* likely the
      finer the spacing simply because more points are placed, and left alone it
      ends the trace.

      The corrector therefore keeps its **best** iterate, accepts after four
      non-improving steps, and refuses to travel more than one predictor step
      from where it started. ``stalledCorrections`` counts those points and
      ``ContourPoint::residual`` says how close each got, to be read against
      ``Contour::correctorTarget``. That is honest rather than convenient: on a
      discontinuous field, a point cannot be closer to "the" level set than the
      field's own ambiguity at a face. A **large** count is a different matter,
      and means the spacing is far below the mesh size.

``fallbackLocations``
   Times the element walk gave up and fell back on ``mfem::Mesh::FindPoints``,
   which is a brute-force scan over element centres. A corrector asks for a
   point location once per iteration, so using it there makes a trace quadratic
   in the mesh; the walk instead transforms back into the previous element, then
   widens over rings of face neighbours. This is a performance statement and
   nothing else — the fallback returns the same element — but a climbing count
   says the walk is failing.

What the tracer is not
~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   * **It is not an exhaustive extraction of the level set.** ``trace()``
     follows the connected component through the point it is given. A second,
     disjoint island at the same level is not found, not reported, and not
     looked for.
   * **It is not X-point aware.** The level set through a saddle is not a
     one-dimensional manifold and the tangent is undefined there, so a trace at
     that level will stall or turn a corner arbitrarily. Locate the critical
     points first; that is what they are for.
   * ``pointAtArcLength()`` **parametrises linearly in polyline length**, not in
     true arc length. It is a way of asking for a point, not a
     reparametrisation.

.. _flux-surfaces-band:

The band, on a curved boundary
------------------------------

On the curved path :math:`\Omega_h` is the union of background elements lying
*inside* :math:`\Gamma`, so :math:`\Gamma_h` is inscribed and there is a band
:math:`O(h)` wide that is inside the plasma and outside the mesh — see
:ref:`curved-band`. The surfaces that band affects are exactly those with
:math:`\Psi_{\mathrm N} \to 1`, since the outermost closed surface *is*
:math:`\Gamma` — which is where the safety factor and the flux-surface averages
are most wanted and least forgiving.

``setBandExtension()`` says which field answers a call out there. There are
three settings, and on a curved boundary only one of them is an answer.

``BandExtension::None``
   **The default.** There is no field outside the mesh, and a trace that reaches
   the edge returns ``LeftMesh``. This is the fitted path unchanged, where
   :math:`\Gamma_h` *is* :math:`\Gamma` and there is no band. It is deliberate
   rather than timid: an arc labelled as a closed contour is worse than a
   refusal.

``BandExtension::FluxTaylor``
   A Taylor step from the point's foot on :math:`\Gamma_h` with :math:`q` frozen
   at the foot — which is exactly what the gridded output does for its band
   nodes (:ref:`output-band`). **Kept as the control, not as a fallback.** Its
   remainder is second order over a band of width :math:`O(h)` however good
   :math:`q` is, so it cannot improve with the polynomial degree and does not.
   Because the gradient is frozen, the extended field is *affine*: contours in
   the band are straight lines, at every :math:`k`.

``BandExtension::TransferLift``
   **The one to use.** The datum on :math:`\Gamma`, plus the line integral of
   :math:`-\gradbar\psi = -r q` back from it, with :math:`q` outside the mesh
   supplied by the method's own extension operator. This is the extension
   technique's own construction — the identity the whole curved path rests on
   :cite:p:`CockburnSolano2012` — so its error is the error of :math:`q`
   integrated along a path of length :math:`O(h)`, which is one order *better*
   than :math:`\psi_h`'s own. It needs the same transfer path the solve was
   given: a different family is a different lifting.

Measured against :math:`\Gamma` itself — which lies entirely in the band at
every mesh, so every point of it is answered by the extension and the exact
answer is known — the lift reaches :math:`k+2` from :math:`k = 2` upwards while
the Taylor step is stuck at second order, and the gap between them widens by
orders with refinement. At :math:`k = 1` the two are comparable, every term in
sight being second order anyway.

.. important::

   **The exact boundary is an anchor and never a correction.** :math:`\psi = 0`
   exactly on :math:`\Gamma`, and :math:`\Gamma` is known analytically, so the
   outermost contour is available for free with no tracing at all. What does
   *not* work is scaling a bad extrapolation towards it — see the same warning
   under :ref:`output-band`, where blending failed for the same reason. The
   error is in *where the field is evaluated*, not in how it is weighted, and
   the lift fixes the evaluation.

.. warning::

   **A band point says so, per point.** ``ContourPoint::extended`` is the flag,
   ``Contour::extendedPoints`` the count and ``Contour::crossesBand()`` the
   predicate, with ``bandDepth`` saying how far outside :math:`\Gamma_h` each
   one sits.

   A count is **not** enough, and that is not a hypothetical: the gridded output
   carried a count rather than a mask, and nothing downstream could tell *which*
   nodes were affected. A band point holds real data and is therefore
   indistinguishable from a solved one without being told. A quantity taken over
   a surface with any band points is known only to the extension's order, and
   must be reported separately from one that is not.

.. note::

   The reach — how far outside :math:`\Gamma_h` a point may sit and still be
   answered — defaults to twice the nearest face's own length, which is
   comfortably more than the deepest excursion measured on the benchmark. It is
   slack on purpose: here the band has to be **crossed**, and a reach too small
   truncates the outermost surfaces on the very faces where the band is widest.
   A reach too small is at least loud: the trace ends as ``LeftMesh`` after a
   handful of points rather than returning a short arc labelled as closed.

.. _flux-surfaces-angle:

Parametrising by poloidal angle
-------------------------------

``fitByAngle()`` re-expresses a closed contour as a radius :math:`\rho(\theta)`
about the magnetic axis, at equispaced angles, each node found by a
one-dimensional Newton along its own ray. For a closed surface a periodic
parameter plus the equispaced trapezoidal rule is the right quadrature
:cite:p:`TrefethenWeideman2014`, and this is what supplies it.

The metric comes from :math:`q`, pointwise. With :math:`x(\theta) = a +
\rho(\theta) u(\theta)` and :math:`u = (\cos\theta, \sin\theta)`, so that
:math:`u \cdot u' = 0`, the curve tangent :math:`t` is known at every node and
:math:`\mathrm{d}x/\mathrm{d}\theta` must be parallel to it. That gives

.. math::

   \rho'(\theta) = \rho \, \frac{u \cdot t}{u' \cdot t},
   \qquad
   \left| \frac{\mathrm{d}x}{\mathrm{d}\theta} \right|
       = \sqrt{ \rho'^2 + \rho^2 },

with nothing differenced anywhere. Two checks that the identity is right, both
of which a transcription error fails: on a circle about the axis the tangent is
perpendicular to the radius, so :math:`u\cdot t = 0` and :math:`\rho' = 0`; and
on a straight line it reproduces the hand-differentiated answer. It is also
invariant under :math:`t \to -t`, so it does not care which way round the tracer
went.

.. danger::

   **The trap this class exists to measure: a spectral rule fed a second-order
   Jacobian is a second-order scheme, and nothing in its output says so.**

   It is easy to build an exponentially convergent quadrature rule and then
   obtain :math:`|\mathrm{d}x/\mathrm{d}\theta|` by central-differencing
   neighbouring node positions. The rule is unchanged, the two versions converge
   to the **same limit**, and the differenced one is orders worse at any finite
   node count.

   ``AngleParametrisation`` therefore returns **three** lengths from one fit:
   ``length()``, the answer, with :math:`\rho'` from :math:`q`;
   ``differencedLength()``, the identical rule fed the differenced metric — the
   trap; and ``chordLength()``, the naive control. If either control ever
   converges as fast as the first, the comparison is empty and the test built on
   it is worthless.

.. note::

   **The first column is not spectral on a discrete contour, and that is not a
   defect.** :math:`\psi_h` jumps across faces, so :math:`\rho(\theta)` is
   piecewise analytic with jumps, and no quadrature rule is geometric on a
   function with jumps. What is measured is that it falls very much faster than
   second order and then **floors**, at the point where the field's own face jump
   converts into a distance. The control that identifies the floor as the *field*
   and not the *rule* is the identical rule run on an analytic contour, where it
   reaches round-off.

.. warning::

   **Star-shapedness is a hypothesis, and it is measured and refused rather than
   assumed.** The denominator :math:`u' \cdot t` vanishes exactly when the ray is
   tangent to the curve, which is where a ray parametrisation stops being one —
   an indented cross-section is the standard failure.
   ``AngleParametrisation::transversality`` is the minimum over the fit and
   ``transverse`` is whether it cleared the floor, mirroring
   ``IndexAudit::transversality`` one stage earlier. The fit **throws** when it
   degenerates, and separately checks the traced curve's own polar angle for
   monotonicity — an independent statement of the same hypothesis, one that does
   not depend on the fit succeeding. A number returned from a degenerate fit is
   worse than no number.

Rays are seeded and **bracketed** from the traced curve rather than from an
assumption: the two traced points straddling each target angle give both a seed
and an interval, so a Newton step that leaves the bracket falls back on
bisection instead of walking off to another branch. ``bisections`` counts those,
and a climbing count is the star-shapedness hypothesis fraying before it fails.

.. note::

   ``stalledRays`` is the exact counterpart of the tracer's
   ``stalledCorrections``, and is there for the same reason: a ray crossing a
   face where the level falls inside the jump has **no** point on it with
   :math:`\psi_h = c` at all, so no tolerance tighter than the jump can be met.
   The fit accepts its best iterate there and reports it. It gets commoner with
   the angle count, more rays being more chances.
