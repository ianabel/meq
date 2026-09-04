Flux-surface geometry
=====================

Once the surfaces exist as curves (:doc:`flux_surfaces`), two things are wanted
from them, and meq supplies both.

**Flux-surface averages** are what a transport code reads: the volume
derivative :math:`V'`, the geometric coefficients that appear in the averaged
transport equations, the safety factor. They are integrals over one surface.

**A representation** is what a code that needs the geometry as a *function* of
the flux label reads: :math:`R` and :math:`z` as smooth functions on a disc
whose centre is the magnetic axis, differentiable with respect to the label.
That is a fit over all the surfaces at once.

.. note::

   As with :doc:`flux_surfaces`, this is library-only: ``meq/SurfaceAverage.hpp``,
   ``meq/Zernike.hpp`` and ``meq/SurfaceFit.hpp``, none of them in the
   ``meq/meq.hpp`` umbrella, and none of it written to any output file.

.. _geometry-averages:

Flux-surface averages
---------------------

.. code-block:: cpp

   #include "meq/SurfaceAverage.hpp"

   meq::AngleParametrisation fit = tracer.fitByAngle( surface, axis, 128 );
   meq::SurfaceAverages avg = meq::surfaceAverages( tracer, fit );

   double vPrime = avg.vPrime;
   double invR2  = avg.inverseRSquared();
   double q      = avg.safetyFactor( g );

The convention is part of the definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   V'(\psi) = \oint \frac{2\pi R \, \mathrm{d}l}{|\nabla\psi|},
   \qquad
   \langle X \rangle_\psi
       = \frac{1}{V'} \oint \frac{2\pi R \, X \, \mathrm{d}l}{|\nabla\psi|}.

**The** :math:`2\pi R` **is in.** A per-unit-length :math:`V'` differs by
exactly that factor, so anything reading a :math:`V'` out of meq is reading the
volume derivative, :math:`\mathrm{d}V/\mathrm{d}\psi` with :math:`V` the volume
the surface encloses. The averaged Grad–Shafranov identity below — which is this
facility's reference-free acceptance — is stated for this convention and not for
the other.

.. warning::

   **The safety factor is** ``safetyFactor()`` **and is never called** ``q``.
   In meq, :math:`q` is the **flux**, :math:`q = \gradbar\psi / r`, a solved
   unknown of the discretisation. The safety factor is also universally written
   :math:`q`. A reader who meets ``q`` anywhere in ``src/meq`` is entitled to
   assume the flux, and one silent conflation would be very hard to find
   afterwards.

One facility, not a list of quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The primitive is: given a surface, take a callable integrand in
:math:`(R, z, \psi, q)` and return :math:`\langle f \rangle` and :math:`V'`.
Everything named — :math:`\langle R^{-2} \rangle`, :math:`\langle
|\nabla\psi|^2/R^2 \rangle`, the arc length, the safety factor — is a one-line
wrapper on it.

That shape is a response to the consumer rather than a preference. The list of
geometry slots a coupled transport code wants is negotiated with the physics
case rather than fixed, so a list that is still moving must not be baked in as a
list of functions. Adding or removing one is a line; the convergence machinery is
written once; and the derivative with respect to :math:`\psi` that a coupling
needs is taken of the facility rather than of eight separate expressions.

.. note::

   :math:`\langle |\nabla\psi|^2/R^2 \rangle` is exactly :math:`\langle |q|^2
   \rangle`. The wrapper is nonetheless written the long way, because it is the
   *quantity* that is being named and a reader checking it against the identity
   below should not have to re-derive the cancellation.

The weight comes from the flux, pointwise
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\gradbar\psi = r q`, so :math:`|\nabla\psi| = r|q|` with :math:`q` read
at the node from the solved flux field. Every weight here is :math:`2\pi R
\,|\mathrm{d}x/\mathrm{d}s|\,\mathrm{d}s / (r|q|)`. **Nothing is differentiated
and nothing is differenced.** This is :ref:`flux-surfaces-q` paying off again.

.. important::

   **The post-processed potential does not buy an order in an average, and the
   reason is structural.**

   The tracer's default pairing roots :math:`\psi^\star`, which converges at
   :math:`k+2` where :math:`\psi_h` converges at :math:`k+1` (see
   :ref:`flux-surfaces-potential`), so the natural expectation is that the
   averages inherit the better order. **They do not.** The level set improves and
   the *weight* does not: the weight divides by :math:`|\nabla\psi| = r|q|`, and
   the enriched flux converges at :math:`k+1` like the raw one. The local
   post-processing buys its extra order in the **potential**, and there is no
   :math:`k+2` flux to be had. An average built on both inherits the worse.

   So every quantity here converges at :math:`k+1` with either pairing, and what
   :math:`\psi^\star` buys is a **constant** — a factor of a few. A constant is
   not an order.

   This is the same shape as the band continuation of :math:`\mathbf{B}` in
   :ref:`output-field`: a quantity limited by the one factor with no solved
   variable behind it.

.. note::

   **The averages are immune to the flux's sign, where the tracer is not.** Only
   :math:`|q|` enters, so handing this facility :math:`-q` changes nothing at
   all — while the same substitution traces the same curve backwards. Neither is
   a defect; the point is that a sign error the tracer survives visibly is
   invisible here, so do not use these numbers as evidence about a sign.

Two extractions, and why both are kept
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two builders. They share the field and the level and nothing else.

``surfaceAverages( tracer, fit )``
   **The primary route.** Equispaced in poloidal angle about the axis, each node
   a one-dimensional ray Newton, the metric from the pointwise identity of
   :ref:`flux-surfaces-angle`, integrated by the periodic trapezoidal rule
   :cite:p:`TrefethenWeideman2014`.

``surfaceAverages( tracer, contour )``
   **The independent route.** Gauss–Legendre on the cubic Hermite segments of
   the traced curve, whose nodes are wherever the step controller put them rather
   than equispaced in anything, whose metric is the interpolant's own, and which
   assumes no star-shapedness at all.

One is a ray parametrisation and the other a predictor–corrector trace with an
interpolant. **They agree or one of them is wrong**, and agreement is worth more
than either being plausible on its own. A missing :math:`2\pi R`, a metric taken
about the wrong point, or a gradient on the wrong side of the division would
each break the agreement, and no single-route table could see any of them.

.. note::

   The two agree only to about their own error, and **cannot** do better:
   :math:`\{\psi_h = c\}` is a union of per-element arcs offset by the face jump,
   so two routes placing nodes differently sample different arcs. The property to
   assert on is therefore the **rate** at which the gap closes, not that the gap
   is small.

   There is a third leg and it is **not built**: a quadrature rule constructed on
   the level set with no curve extracted at all. Read the two-route agreement as
   two legs of three.

.. note::

   A third check is free and is not a third extraction. The angle builder's nodes
   are laid out about a ray *origin*, and the integral cannot depend on it — a
   reparametrisation of the same curve integrates to the same number. Fitting the
   same contour about a deliberately displaced origin and getting the same
   :math:`V'` is a real check on the metric, and it costs one more fit.

.. danger::

   **The metric trap of** :ref:`flux-surfaces-angle` **does not go away because a
   quantity is an average, and the plausible argument that it does is wrong.**

   :math:`\langle X \rangle` is a *ratio* of two integrals over the same metric,
   so a bad metric looks as though it ought to cancel. It cancels a **constant**
   and **nothing in the order**: with the metric differenced rather than taken
   pointwise from :math:`q`, the averages revert to second order in the node
   count exactly as an arc length does, and the two columns separate by orders
   by the time the node count is large.

   ``differencedWeight`` and ``averageDifferenced()`` are therefore first-class
   members rather than test scaffolding. They are filled only by the
   equispaced-angle builder, because that is the parametrisation the trap lives
   in; ``differencedAvailable()`` says so, and ``averageDifferenced()`` **refuses**
   rather than quietly returning a number built on the pointwise metric — which
   would make the control agree with the answer for the worst possible reason.

.. _geometry-identity:

The reference-free acceptance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are **no closed-form Solov'ev flux-surface averages**. :math:`\psi` is
elementary, but :math:`V'`, :math:`\langle R^{-2}\rangle` and the safety factor
are integrals over a *contour* of it, and those contours have no elementary arc
length. What the test fixtures supply is a converged **reference value**, and
that distinction should be made wherever the number is printed.

So the sharpest check is one that needs no reference value at all. Taking
:math:`\langle \nabla\cdot G \rangle = (1/V')\,\mathrm{d}(V'\langle G \cdot
\nabla\psi\rangle)/\mathrm{d}\psi` with :math:`G = \nabla\psi / R^2`, and using
:math:`\dstar\psi = -F`,

.. math::

   \frac{1}{V'} \frac{\mathrm{d}}{\mathrm{d}\psi}
       \left( V' \left\langle \frac{|\nabla\psi|^2}{R^2} \right\rangle \right)
   = - \left\langle \frac{F}{R^2} \right\rangle.

Three averages checking each other with nothing but the equation.
``averagedGradShafranovResidual()`` computes it.

.. important::

   **The right-hand side is written with the** :math:`F` **the solver is
   actually fed**, and never as :math:`-\mu_0 p' - g g' \langle R^{-2}\rangle`.
   The second is Solov'ev-specific and is re-derived by hand, so it is not
   independent of the hand that derived it. That is the same discipline the
   analytic fixtures follow when they recompute :math:`\dstar\psi` by
   differences rather than trusting a transcription (see :ref:`testing-fixtures`).

.. danger::

   **The** :math:`\mathrm{d}/\mathrm{d}\psi` **must be Richardson-extrapolated**,
   :math:`(4D(h/2) - D(h))/3`, and ``FluxDerivative::Richardson`` is the default
   for that reason. A plain central difference carries its own second-order
   truncation, so the comparison is floored by the **instrument** rather than by
   the identity.

   Measured on a discrete field over a sixteenfold refinement, the plain column
   **stops moving at all** while the extrapolated one continues to fall at
   fourth order. A column that does not converge under mesh
   refinement is measuring the instrument, not the answer — and an identity
   checked with that instrument would pass with a real defect underneath it.
   ``FluxDerivative::CentralDifference`` is kept live so that this is a
   measurement and not a claim.

.. warning::

   **The step is a real choice and it is not monotone in either direction.** Too
   large and the extrapolation's own truncation dominates; too small and the
   difference of two nearly equal integrals divides the surfaces' own face-jump
   noise by the step and amplifies it. Measured on the suite's fixture, halving
   the step from the value near the optimum made the residual **worse**, not
   better. There is an optimum; find it on your own problem rather than
   assuming smaller is better.

Surfaces that cross the band
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SurfaceNode::extended`` is read from the producing object per node, and
``SurfaceAverages::extended`` is what a consumer drops a surface on before
computing an error norm or differencing two runs. A surface with any band node is
limited by the **extension's** order rather than by the discretisation, and must
be reported separately; quoting one rate over both populations hides precisely
what :ref:`flux-surfaces-band` exists to measure.

Every one of those flags reads false on the fitted path, which is not the same
statement as "the flag is not wired" — there is no band there, so zero is the
correct answer rather than an absent one, and the suite asserts the zero rather
than assuming it.

.. note::

   The contour builder marks a whole **segment**, deliberately. Its nodes are
   Gauss points strictly between two accepted points, so a node cannot say for
   itself; one is therefore marked extended when *either* endpoint of its segment
   is. That over-reports by at most one segment at each end of a band excursion
   and never under-reports. Under-reporting is the failure that matters: it is
   band data presented as solved data.

What the averages are not
~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   * **Not valid on an open surface.** Every rule here is periodic or closed and
     assumes the contour returns to its start. An open surface terminating on the
     domain boundary has endpoints and wants a different rule; that is deferred
     with free boundary.
   * **Not X-point aware.** Approaching the separatrix, :math:`1/|\nabla\psi|`
     diverges, the surface develops a corner, and every quantity here diverges
     logarithmically with it. That is the physics. Where to cut is a decision to
     be recorded rather than discovered.
   * ``safetyFactor()`` **carries no evidence of its own.** It is an algebraic
     combination of :math:`V'` and :math:`\langle R^{-2}\rangle` — both of which
     *are* measured — times a :math:`g(\psi)` the caller supplies from its own
     profile, so it inherits their rates and nothing more. What the suite pins is
     its algebra and its :math:`4\pi^2`, which is worth doing since a constant is
     exactly as easy to get wrong here as anywhere and no rate table would see
     it.
   * ``arcLength()`` **is not an independent measurement of the metric.** Since
     every weight already contains the arc-length element, this integrand
     recovers it exactly and the answer is *algebraically* equal to the producing
     object's own length. What it checks is that the weights are built the way
     the convention above says they are.

.. _geometry-disc:

The surfaces as a map from a disc
---------------------------------

The second deliverable is :math:`R` and :math:`z` as functions of a flux label
and an in-surface angle — a map from a disc whose centre is the magnetic axis and
whose edge is the outermost surface. ``meq::SurfaceFit`` is that map, fitted to
the traced points.

.. code-block:: cpp

   #include "meq/SurfaceFit.hpp"

   std::vector<meq::SurfaceSample> samples;   // one per traced point
   // ... fill from each surface's AngleParametrisation ...

   samples = meq::relabelByAxisShape(
       samples, meq::axisShapeFromSamples( samples, axis.r, axis.z ) );

   meq::SurfaceFit fit( 12, samples );

   meq::GaugeFreeFitReport report;
   meq::SurfaceFit better = meq::gaugeFreeFit(
       fit, field, meq::discNodesFrom( samples ),
       meq::GaugeFreeFitOptions(), report );

The relabelling line is **not optional**, and the gauge-free refit is what makes
the whole thing converge near the axis. Both are explained below.

.. note::

   ``SurfaceFit`` and ``Zernike`` are **MFEM-free**: plain doubles in,
   coefficients out. That is what makes them unit-testable without the finite
   element library and buildable in continuous integration, which cannot obtain
   the MFEM branch meq needs (see :ref:`organization-mfem-free`). The price is
   that a caller writes the two-line loop that turns an ``AngleParametrisation``
   into samples, rather than the fit taking a dependency on the tracer for it.

.. _geometry-radius:

The radial coordinate is :math:`\sqrt{\Psi_{\mathrm N}}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. danger::

   **Getting this wrong is silent.**

:math:`\psi` has a *quadratic* extremum at the magnetic axis — a smooth function
with a non-degenerate critical point — so :math:`\psi - \psiax` falls off like
(distance)\ :sup:`2` and the normalised flux :math:`\Psi_{\mathrm N}` behaves
like (distance)\ :sup:`2` near the axis. The geometry is smooth in the
**distance**, so it is smooth in :math:`\sqrt{\Psi_{\mathrm N}}` and carries a
square-root branch point in :math:`\Psi_{\mathrm N}` itself.

Parametrise by :math:`\Psi_{\mathrm N}` directly and *every* basis — Zernike,
Chebyshev, anything — converges algebraically against that branch point, worst
near the axis, **with nothing in a convergence table to say why**. The rate
simply comes out below the design order and looks like a discretisation problem.
It is the same species of defect as a wrong Jacobian (:ref:`sources-jacobian`):
the answer converges, to the wrong rate, for a reason no assertion in the suite
is looking at.

So the disc radius is :math:`\varrho = \sqrt{\Psi_{\mathrm N}}`, and
``radiusFromNormalisedFlux()``, ``normalisedFluxFromRadius()`` and
``fluxDerivativeFromRadial()`` are the conversions. Use them rather than writing
a square root and a chain rule at each call site: the chain factor is
:math:`1/(2\varrho)`, and it is the thing that gets dropped.

.. note::

   That factor is **unbounded at the axis, and that is the coordinate rather
   than a defect.** Even for a representation perfectly smooth in Cartesian
   coordinates, :math:`\mathrm{d}/\mathrm{d}\Psi_{\mathrm N}` at fixed angle
   diverges at the axis whenever the :math:`m = 1` content is non-zero — those
   modes are the rigid shift of a surface, and shifting the axis sideways at a
   rate :math:`1/(2\varrho)` is exactly what a Shafranov shift does as the
   surfaces shrink onto a point. A caller wanting a quantity that stays finite
   there wants ``radialDerivative()``, which is in the disc radius. Its growth as
   the innermost surface approaches the axis is a **conditioning number for a
   coupling to read**, not a defect to fix.

.. _geometry-zernike:

Why Zernike
~~~~~~~~~~~

The Zernike polynomials on the unit disc are

.. math::

   Z_l^m(\varrho, \theta) = R_l^m(\varrho) \times
       \begin{cases} \cos(m\theta) & m \ge 0 \\
                     \sin(|m|\theta) & m < 0 \end{cases}

with :math:`R_l^m` containing only the powers :math:`\varrho^{|m|},
\varrho^{|m|+2}, \ldots, \varrho^{l}`, admissible exactly when :math:`l \ge |m|`
and :math:`l - |m|` is **even**.

**The index constraint is the entire point of choosing the basis.** A naive
tensor product — any polynomial in the radius times any Fourier mode in the angle
— admits terms such as :math:`\varrho^2\cos\theta`, which in Cartesian
coordinates is :math:`x\sqrt{x^2+y^2}`, whose second derivative has no limit at
the origin: it reads a different number down every ray. Such a basis makes the
centre of the disc a coordinate singularity that every consumer has to
special-case, and the special case is where the bugs live.

The parity-and-minimum-power constraint excludes exactly those terms, and what
survives is the statement worth having:

.. important::

   **Every admissible** :math:`Z_l^m` **is a bivariate polynomial in**
   :math:`(x, y)` **of degree** :math:`l`.

   So a Zernike expansion is a polynomial in Cartesian coordinates, **the
   magnetic axis is an ordinary interior point of it**, and derivatives there are
   whatever the polynomial says they are. No limit to take, no ray to choose, no
   branch — at precisely the place where the surfaces shrink to a point,
   :math:`1/|\nabla\psi|` diverges, and the in-surface angle is undefined. That
   is a better argument for the basis than "DESC does it".

A consequence that comes free and exactly: every mode with :math:`m \ne 0`
carries :math:`\varrho^{|m|}` and above, so all of them **vanish** at the centre
and what is left there is one number with no angle in it. A Zernike fit
therefore puts every surface's collapse point at the same place *by
construction*, whatever the coefficients are and however badly the fit was
conditioned. ``axis()`` is that point; ``axisAtAngle()`` exists to measure that
it does not depend on the angle, and to measure what the tensor-product control
does instead, where it is a **curve** rather than a point.

.. note::

   What is *not* free is whether that point is the magnetic axis. A sample set
   has a hole in the middle — nobody traces the surface at
   :math:`\Psi_{\mathrm N} = 0`, because it is a point — so ``axis()`` is an
   extrapolation into the hole, and the suite measures how far it lands from the
   axis the critical-point finder reports.

.. note::

   The radial polynomial is a Jacobi polynomial under a change of variable, and
   that is how it is evaluated. **The explicit factorial sum every reference
   prints is not used**: its terms are binomial-sized while its answer is order
   one, so it costs roughly one digit per decade of the largest term to
   cancellation, and a fit wanting twenty or thirty modes would be reading noise.
   It is kept in the tests as a *control* rather than as an implementation.

.. _geometry-angle:

The angle is a label, and the obvious choice is the wrong one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. danger::

   **The geometric poloidal angle is inadmissible, and everything else in this
   section is downstream of it.**

Labelling surfaces by the geometric angle about the axis is the obvious choice:
it is what ``fitByAngle()`` produces, it is consistent across surfaces because
the axis is one point, and it needs no extra information. It is also, for any
surface that is not a **circle**, a parametrisation in which the disc map is not
smooth at the axis — so no basis that is smooth there can converge against it.

The argument is exact and takes one line. A function smooth at the origin has an
expansion :math:`f(0) + ax + by + O(\varrho^2)`, so its coefficient of
:math:`\varrho^1` carries exactly **one** angular harmonic. Take nested similar
ellipses of semi-axes :math:`A` and :math:`B`, labelled by geometric angle:

.. math::

   R - R_{\mathrm{axis}} = \varrho \,
       \frac{A B \cos\theta}{\sqrt{A^2 \sin^2\theta + B^2\cos^2\theta}},

whose coefficient of :math:`\varrho^1` carries :math:`\cos\theta`,
:math:`\cos 3\theta`, :math:`\cos 5\theta` and so on for ever. Every one of
those beyond the first is a mode the Zernike index constraint **excludes** —
precisely because it is not smooth there. The parametrisation puts content
exactly where the basis refuses to look.

Measured on the ellipses themselves, where the answer is known exactly: under the
ellipse's own parameter the family is a polynomial of degree one and the fit
reproduces it to round-off at degree two; the **identical points** relabelled by
their geometric angle are wrong at degree two and still wrong at degree twenty,
decaying algebraically. Same points, same basis, same solver — only the label
differs.

**The repair needs no field.** Near the axis every equilibrium's surfaces are
ellipses, since :math:`\psi` has a non-degenerate critical point there and its
level sets are those of a quadratic form. Two numbers describe such an ellipse up
to scale — how it is turned and how flat it is — and both are recoverable from
the innermost traced surface alone. ``AxisShape`` is that ellipse,
``axisShapeFromSamples()`` and ``axisShapeFromHessian()`` are the two routes to
it, and ``relabelByAxisShape()`` turns geometric angles into the label that makes
the leading term smooth. It is worth orders of magnitude in the worst fit
error, at no cost in conditioning.

.. important::

   **It is an exact reparametrisation, not an approximation.** The map from the
   geometric angle to the label is a fixed bijection of the circle, so a caller
   asking "where is this surface at geometric angle :math:`\theta`" evaluates the
   fit at ``shapedPoloidalAngle( shape, theta )`` and gets the right point, at
   every radius. What is approximate is only the *claim* that this label makes
   the map smooth, and that claim is exact at leading order and no better.

.. warning::

   ``SurfaceFit::position()`` **takes the label the samples carried.** A caller
   who fitted relabelled samples and then asks at a geometric angle gets the
   wrong point on the right surface. Ask at
   ``shapedPoloidalAngle( shape, theta )``, with the same shape the samples were
   relabelled by.

And a relabelling that does not depend on the radius has exactly one function's
worth of freedom, so it can fix exactly **one order**: it makes the
:math:`\varrho^1` coefficient a single harmonic and leaves the rest carrying
whatever harmonics the surface shaping puts there. Pull the innermost surface in
towards the axis and an algebraic tail appears that no degree removes. That is
the tension the next section resolves.

Conditioning: the hole and the disc edge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Zernike is orthogonal on the disc under the uniform measure, so a
well-distributed sample gives a well-conditioned problem — but the sample is
whatever the caller traced. The conditioning is therefore a measurement and not
an assumption, and ``SurfaceFitDiagnostics::conditionNumber`` is it, **of the
design matrix** rather than of the normal matrix; the two differ by a square, and
quoting the wrong one understates the difficulty by half its digits.

.. warning::

   **The plausible advice about where to put the surfaces is wrong.** The disc
   measure is :math:`\varrho \,\mathrm{d}\varrho = \mathrm{d}\Psi_{\mathrm N}/2`,
   so Gauss–Legendre levels in :math:`\Psi_{\mathrm N}` with equispaced angles
   *ought* to make the discrete inner product the continuous one and the fit a
   projection. Measured over equispaced-in-flux, equispaced-in-radius and
   Gauss-in-flux layouts on the same surfaces, the three condition numbers agree
   closely at every inner limit tried, and Gauss is sometimes the **worst** of
   the three.

   The orthogonality argument needs nodes spanning the whole disc, and a sample
   set has a **hole** in the middle.

What decides the conditioning, by orders, is that hole and the disc edge. Pulling
the innermost surface out from the axis costs orders of magnitude. And the disc
edge is a genuine choice:

``discEdge = 1``
   :math:`\varrho = \sqrt{\Psi_{\mathrm N}}` literally, so a coefficient means
   what it says and two fits with different outermost surfaces are directly
   comparable. **The default**, on that ground and no other.

``discEdge`` set to the outermost fitted flux
   :math:`\varrho = \sqrt{\Psi_{\mathrm N}/\Psi_{\max}}`, so the samples reach the
   edge of the disc. Worth orders of magnitude in conditioning.

.. important::

   **The two fit the same function.** A Zernike expansion of degree :math:`L`
   spans exactly the polynomials of degree :math:`L` in :math:`(x, y)`, and that
   space is closed under scaling — so rescaling the disc is a change of **basis**
   and not of model, and the two fits agree to every digit. The choice is purely
   conditioning.

   **Rescale unless you need coefficients comparable across fits.** The default is
   what it is because a silent change of coordinate is worse than a conditioning
   number a caller can read, and ``majorRadiusExpansion()`` refuses to hand back
   a bare ``ZernikeExpansion`` unless the basis, the coordinate *and* the edge are
   all the plain ones — such an object carries no record of any of the three.

.. note::

   The solve is a column-pivoted QR followed by a singular value decomposition,
   **never normal equations**, which would square the condition number: a design
   matrix an ill-placed sample set reaches easily would then lose every digit
   through a route whose only symptom is a fit that is merely poor. Directions
   whose singular value falls below a relative floor are **discarded** rather than
   inverted, and the count is reported.

Controls, kept live
~~~~~~~~~~~~~~~~~~~

Two of the decisions above change the convergence *rate* with nothing in the
output to say why, so each keeps its losing alternative buildable — in the same
way :ref:`flux-surfaces-band` keeps its Taylor step and the averages keep their
differenced metric.

``FitRadialCoordinate::NormalisedFlux``
   Feeds the basis :math:`\Psi_{\mathrm N}` instead of its square root. This is
   what turns :ref:`geometry-radius`'s argument into a measurement, and it loses
   by orders in the worst fit error while leaving the **conditioning untouched**
   — which is what says the effect is the branch point and not a conditioning
   artefact.

``FitBasis::TensorProduct``
   Any power of the radius times any Fourier mode, with no parity constraint. It
   is a **generous** control, carrying more modes than Zernike at the same degree,
   so anything it does worse it does worse with more freedom rather than less.
   Measured, it fits the sample cloud considerably better and is useless: the
   condition number is at the edge of double precision, the extrapolated axis is
   orders worse, and the axis **moves depending which angle you approach it
   along**. A better residual on the data you fitted and nothing anywhere else.

.. note::

   A caveat for anyone rebuilding that second control: a tensor product that
   keeps :math:`l \ge |m|` *also* gives an exactly angle-independent axis,
   because those modes vanish at the centre too. Only admitting modes with **no
   radial factor** breaks it, and a control built the obvious way demonstrates
   nothing.

Neither is an answer. Both are returned by nothing except a caller who asked for
them by name.

.. _geometry-gauge-free:

The gauge-free fit
~~~~~~~~~~~~~~~~~~

Everything above pins the angle as an **input**: a ``SurfaceSample`` carries a
:math:`\theta`, so the fit is asked to put a particular point at a particular
label. That is what leaves the algebraic tail near the axis, and it is a
structural limit rather than a tuning problem — a parametrisation smooth to all
orders has to vary with the surface, which is what a code that *solves* for its
coefficients gets for free and a post-hoc fit does not. The DESC authors say as
much of their own solver: it enforces no poloidal angle constraint and ends up
finding a good representation through the course of the optimisation.

**So ask for less.** Require each disc node only to *land on the right surface*:

.. math::

   \min_{c_R,\, c_z} \; \sum_j
       \left[ \Psi_{\mathrm N}\big( x(\varrho_j, \theta_j) \big) - \Psi_j
       \right]^2

with :math:`x` the expansion being solved for and :math:`\Psi_j` the surface
node :math:`j` belongs on. Nothing says where *along* its surface a node must
sit, so the angle is free and the truncated basis may choose
it. ``meq::DiscNode`` is the input, and it deliberately **carries no position** —
that is the entire difference from a ``SurfaceSample``.

It needs no force balance and no second physics solver, because meq already has
:math:`\psi`; and the Jacobian needs :math:`\nabla\Psi_{\mathrm N}`, which is the
**solved flux** once more. Measured on nested ellipses, where the answer is known
exactly, the gauge-free solve reaches round-off at low degree where the
prescribed-angle fit never converges at all — and that is the theoretically right
answer rather than a lucky one, since for similar nested ellipses the family
genuinely *is* a degree-one map under the correct angle. On a discrete field it
removes the algebraic tail and the coefficient envelope resumes geometric decay.

.. important::

   **The warm start is not an optimisation.** The surface residual is minimised
   by *any* map onto the right surfaces, including folded ones and ones that
   traverse a surface twice, so the basin matters. The linear fit is what puts
   the iterate in the right one; starting from zero coefficients is not supported
   and would not be meaningful.

.. note::

   **The residual is scaled into metres, and that is not cosmetic.**
   :math:`\Psi_{\mathrm N}` is dimensionless and its gradient varies by orders
   across the disc, so a least-squares problem in the flux residual weights outer
   surfaces against inner ones by whatever :math:`|\nabla\Psi_{\mathrm N}|`
   happens to be. Dividing each row by it makes the residual the **normal
   distance from the node to its surface, in metres** — which is the quantity a
   reader cares about, and is itself gauge invariant.

Fixing the gauge
~~~~~~~~~~~~~~~~

The freedom is a tangential slide, one function's worth, and the linearised
system cannot see it: moving a node along its own surface changes no residual.
The expectation going in was therefore that the Gauss–Newton matrix is **rank
deficient** along the gauge, by about half its columns.

.. warning::

   **Measured, it is not, and the answer depends on the field.** On nested
   similar ellipses the exactly-null directions are a handful at every degree,
   and on a Solov'ev equilibrium there are **none at all**. What both have
   instead is a long **soft tail with no gap in it**, running smoothly down to a
   tiny fraction of the largest singular value.

   The reason is that a tangential slide is exactly null only when the displaced
   map is still *in* the truncated space, and generally it is not. So the gauge
   is **not a subspace to be projected out**; it is a direction in which the
   problem is merely very soft — and the floor deciding which directions to
   invert therefore has to be *chosen* rather than read off a gap.

   **And the soft tail is enough on its own.** Removing the gauge on the field
   where nothing is exactly singular still leaves the map folded and the fit
   orders worse. A field with no exact null space is not a field that can do
   without a gauge.

``SurfaceGauge::MinimumNorm``
   **The default.** Pseudo-inverse steps in which directions below a relative
   singular-value floor are simply not moved: among the ways to satisfy the
   surface constraint it prefers the smallest step, which is a spectral
   condensation preference arrived at for free and with no weights to choose.

   It is min-norm **relative to the warm start**, not absolutely: each *step*
   avoids the null directions, so the answer stays as close as it can to the
   coefficients the linear fit arrived with. That is a property to like — the
   linear fit's label is a sensible place to be pinned — but it is not "the
   smallest coefficients", and a reader who assumed it was would be wrong.

``SurfaceGauge::SpectralWidth``
   An explicit quadratic penalty on the angular content, which is Hirshman and
   Breslau's spectral condensation and the constraint VMEC uses where DESC lets
   the solve find the angle. It fixes the gauge outright rather than by
   preference, and it **biases** the surface residual, so its weight is a real
   tuning parameter where the floor above is a threshold on a cliff.

   .. note::

      **It loses on its own metric as well as on accuracy, which was not
      expected.** The spectral width is a *ratio* of two weighted sums of the
      same coefficients, and a quadratic penalty shrinks both — so driving the
      coefficients down does not drive the ratio down. Measured, sweeping the
      weight across its whole useful range moves the width **the wrong way** and
      costs more than an order of magnitude in surface error, while the
      minimum-norm step, which asks nothing about the spectrum at all, already
      sits below every penalised value. Hirshman and Breslau minimise the width
      *itself*, which is not a quadratic problem, and a quadratic surrogate for
      it is not the same thing. Kept because a comparison with no losing column
      is not a comparison.

``SurfaceGauge::None``
   **The control.** Every direction inverted however small its singular value —
   and **no trust region either**.

.. important::

   **The trust region is itself a gauge, which the control has to know.** What
   actually makes this robust is adaptive Levenberg–Marquardt damping: as the
   damping falls the step tends to the plain pseudo-inverse, so on an easy problem
   it is inert, and where it is not it is what stops a soft direction being
   inverted into an enormous step. Measured, the undamped pseudo-inverse reaches
   round-off at moderate degree and **fails** at high degree, while the damped one
   converges at every degree tried.

   Damping the step *is* choosing among the directions the constraint does not
   determine, so ``SurfaceGauge::None`` disables the damping as well as the floor.
   A "no gauge" control that kept the trust region would be a control of nothing,
   and it would quietly pass.

.. danger::

   **Check the map for folding. A beautiful residual over a folded map is the
   quiet wrong answer this class is most exposed to.**

   Requiring only that nodes *land* on their surfaces says nothing about their
   order along one, so nodes may bunch, cross, and turn the map over while every
   residual stays perfect. ``SurfaceFit::mapJacobian()`` is the safeguard: sample
   it over the fitted annulus and check that it keeps **one sign**.

   *Which* sign is not fixed and must not be assumed — it is positive for an
   angle running one way round and negative for the other, and both are
   legitimate labels. It vanishes at the centre for the ordinary reason polar
   coordinates do, so a check that includes the centre measures that and not a
   fold. ``GaugeFreeFitReport`` reports the minimum and maximum; the ungauged
   control produces a negative minimum on every field tried.

.. note::

   ``GaugeFreeFitReport::worstExcursion`` is the worst surface error at *any*
   accepted iterate, and it is there because a final error alone cannot tell a
   converged run from one that wandered and came back.

A free angle changes what can be asserted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any check that compares the fit's position *at a given angle* against a surface
traced at that angle is measuring precisely the freedom being granted. So the
gauge-free fit is accepted on **gauge-invariant** properties instead: the
distance from each node to its surface, the fitted perimeter against the exact
one, and the sign of the map Jacobian.

.. important::

   **This costs the consumer nothing, and that is the point.** A coupled
   transport code reads flux-surface *averages*, and an average does not know how
   its surface was parametrised. The deliverable was gauge invariant all along.

What the representation is not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   * **Not a fit of open surfaces.** Every mode is periodic in the angle and
     regular at the disc centre. An open surface terminating on the domain
     boundary has endpoints and wants a different basis; deferred with free
     boundary.
   * **Not an equilibrium solve.** Nothing here enforces force balance.
     :math:`\psi` is an input, and the only thing solved for is where the disc
     coordinates put their points.
   * **It does not constrain the axis.** A node at :math:`\Psi_{\mathrm N} = 0`
     would ask for the flux to vanish, which is true at the axis and would pin it
     for free — except that :math:`\nabla\Psi_{\mathrm N}` vanishes there, so the
     row scaling divides by zero and the linearisation carries no information
     about where to move. ``discNodesFrom()`` therefore refuses such a node, and
     the axis stays an extrapolation into the hole whose accuracy is measured
     rather than imposed.
   * **It does not know where its samples came from**, so it cannot flag a
     surface as band-limited. Mixing surfaces that crossed the band with surfaces
     that did not fits two different accuracies into one expansion, which is a
     decision to be taken deliberately. ``Contour::crossesBand()`` and
     ``AngleParametrisation::crossesBand()`` are where a caller finds out.
