Output
======

Every run writes the same equilibrium **three times**, in three formats. That is
not redundancy: no single format is simultaneously exact, portable and
convenient, and each of the three gives up a different one of those. On request
there is a fourth file, which is not the equilibrium in another resolution but a
**reduction** of it: the flux surfaces and the integrals over them.

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Files
     - Read by
     - What it is
   * - ``<stem>.mesh``, ``<stem>_psi.gf``, ``<stem>_grad_psi.gf``
     - GLVis; MEQ itself
     - **Exact.** Every polynomial coefficient of the finite element solution.
       This is the restart format, and the only one from which MEQ can resume a
       run. It carries :math:`\psi_h`.
   * - ``<stem>_psistar.gf``
     - GLVis
     - The **post-processed** potential :math:`\psi^\star`, in
       :math:`P_{k+1}`. See :ref:`output-which-potential`.
   * - ``<stem>/<name>.pvd`` and ``<stem>/Cycle000000/``
     - ParaView, VisIt
     - **The picture.** VTK Lagrange cells at degree :math:`k+1`, with the mesh
       boundary bent onto the true :math:`\Gamma`.
   * - ``<stem>.nc``
     - anything that reads NetCDF
     - **The interchange format.** :math:`\psi` and :math:`\mathbf{B}` sampled on
       a uniform :math:`(R, Z)` grid. Lossy, and portable.
   * - ``<stem>_surfaces.nc``
     - a 1-D transport code; ``plot_equilibrium.py``
     - **The reduction.** The flux surfaces themselves and the flux-surface
       averages over them, against a flux label. Written only when
       ``[output] FluxSurfaces`` is set. See :ref:`output-flux-surfaces`.

``<stem>`` is the output directory and prefix from the ``[output]`` table. MEQ
does **not** create the output directory; if it does not exist the run exits 3
(see :ref:`running-exit-codes`).

.. note::

   **None of the first three carries the flux surfaces**, and that is what the
   fourth file is for. A consumer wanting :math:`V'`,
   :math:`\langle R^{-2}\rangle` or the shape of a surface reads
   ``<stem>_surfaces.nc``; one wanting the family as a differentiable map from a
   disc, or wanting to drive the extraction itself, links against the library —
   see :doc:`flux_surfaces` and :doc:`surface_geometry`.

.. _output-which-potential:

Which potential each file carries
---------------------------------

.. important::

   **Everything MEQ draws or exports carries** :math:`\psi^\star`, **the
   post-processed potential** — the ``.vtu``, the ``.nc``, and
   ``<stem>_psistar.gf``. It converges one order faster than the solved
   :math:`\psi_h` at essentially no cost, so there is no reason to report the
   worse field.

   **The one exception is** ``<stem>_psi.gf``, **which keeps** :math:`\psi_h`
   deliberately. That file is the restart format and is read back into a
   degree-\ :math:`k` potential space; :math:`\psi^\star` lives in
   :math:`P_{k+1}` and would not fit. Making the two consistent would break
   restart, which is why the asymmetry is commented at the write site.

Two consequences for anyone reading MEQ's output:

* **A** ``.nc`` **differenced against a run from before this change measures the
  post-processing, not the physics.** The variable has the same name, units,
  layout and mask; only its meaning moved. The file records which potential it
  holds in a global attribute ``potential``, and the *absence* of that attribute
  is what identifies an older file.
* **VTK point counts rose**, because the Lagrange cells went from degree
  :math:`k` to :math:`k+1` to match the field they now draw — a factor of 1.5 on
  triangles at :math:`k = 3`.

.. warning::

   Making :math:`\psi^\star` the reported potential puts every output behind
   MFEM's local reconstruction, which had a defect that corrupted it **per
   element** wherever :math:`\partial F/\partial\psi` vanishes — invisible to
   any whole-domain norm. Until this change that defect could only reach the
   error estimator; it now reaches the primary outputs. The fix is on one MFEM
   branch and no other, and ``INSTALL.md`` in the source tree says which. It is
   regression-tested rather than latent: see :ref:`postprocessing-singular`.

The exact format
----------------

MFEM's own mesh and grid-function files, written at full precision. ``_psi.gf``
is the solved potential :math:`\psi_h`; ``_psistar.gf`` is the post-processed
:math:`\psi^\star`; ``_grad_psi.gf`` is the HDG flux
:math:`q = \gradbar\psi / R` — the solved unknown, **not** the magnetic field,
which is a relabelling of it (see :ref:`output-field`).

.. code-block:: sh

   glvis -m run.mesh -g run_psi.gf

On a curved run the mesh written here is :math:`\Omega_h`, the polygonal
computational domain, not the true plasma region: it is the mesh actually solved
on, which is what a restart needs.

.. _output-vtk:

The VTK format
--------------

Written through ParaView's data collection format, as **Lagrange cells at the
degree of the field being drawn** — which is :math:`k+1`, since what is drawn is
:math:`\psi^\star`.

That is a deliberate choice with a quiet failure mode behind it. VTK's native
cells are linear, so the default path draws a cubic solution as though it were
linear — and the result does not look like a bug. It looks like a coarse mesh.
A picture that is silently one to three orders less accurate than the
computation it came from is worse than no picture, so the point count is
asserted against the vertex count in the test suite for exactly this reason.

**The mesh boundary is bent onto** :math:`\Gamma`. On the curved path the
computational domain :math:`\Omega_h` is inscribed in the plasma, so its
boundary :math:`\Gamma_h` is a chord of the true boundary and there is a band
between them (see :ref:`output-band`). For the picture, MEQ installs a curvature
on the mesh and moves each boundary face out onto :math:`\Gamma`. Since the VTK
was already high-order cells, this needed nothing further from the format; the
two features compose.

The bending is per-node and backs off where a node cannot reach without tangling
its element. Two things about that were got wrong on the way and are worth not
repeating:

* **Smoothing the displacement into the interior on the degree-of-freedom graph
  couples the two coordinates**, because they share one index range — so a
  radial displacement gets averaged against a vertical one. Measured, that made
  tangling *worse* than moving the boundary alone.
* **Backing the displacement off globally costs every face the worst face's
  limit.** Per-node backoff reaches a far larger fraction of the boundary.

The gap can exceed an element's own size, so some faces genuinely cannot reach.
The driver reports the fraction that did, and warns when it is not all of them.

An adaptive run also writes ``<stem>_cycles/``, one VTK frame per cycle, which
ParaView scrubs through as a time series. It is a **separate collection** from
the answer, because ``<stem>`` gets its boundary bent onto :math:`\Gamma`, and
doing that mid-loop would hand the next refinement a geometry the estimator
never saw.

.. note::

   The frames carry :math:`\psi^\star`, as the answer does — *which field is
   drawn* and *which geometry it is drawn on* are separate questions. The
   geometry stays as solved, faceted and unbent; the field is the better one
   either way, and it costs nothing, since the post-processing has already been
   done for the estimator by the time each frame is written.

.. note::

   One frame per cycle is easy to get *almost* right, and the near miss is
   invisible. Rebuilding the collection object for each frame puts every cycle
   directory on disk with the correct refined mesh, and leaves the ``.pvd``
   index listing only the last of them — so ParaView opens the file, shows a
   single frame, reports no error, and all the data is present on disk. The
   collection appends to its index and never scans the directory, so it has to
   survive between frames and be rebound to each new mesh.

.. _output-netcdf:

The gridded format
------------------

:math:`\psi` and :math:`\mathbf{B}` sampled onto a uniform :math:`(R, Z)` grid,
whose extent is :math:`\Gamma`'s bounding box when there is a shape and the
``[mesh]`` box otherwise. This is the format to give to anything that is not a
finite element code.

Variables: the grid coordinates; ``psi`` (which is :math:`\psi^\star` — see
:ref:`output-which-potential`); ``B_R`` and ``B_Z``; a byte mask ``inside``; a
byte mask ``extrapolated``; and the boundary curve as a pair of coordinate
arrays. A rotating run additionally carries a density for each species
and the electrostatic potential. Nodes outside the plasma are ``NaN`` *and*
carry ``inside = 0``, so a reader that checks either one is safe.

The global attributes carry the run's provenance — the configuration file name,
the MEQ and MFEM versions, the polynomial degree, the element count, the Newton
iteration count and final residual, the boundary treatment, and for adaptive
runs the cycle count, the final estimator and the marking strategy. This is what
makes a ``.nc`` file self-describing enough to plot six months later.

**The two-dimensional variables are indexed** ``(Z, R)``, **with** ``R``
**varying fastest.** That is C row-major order, and it is the layout a transport
code coupling to MEQ along :math:`(t, x)` already expects.

.. _output-restart-variables:

The same file also carries the exact solution, and the conductors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beside the rasterization above, ``<stem>.nc`` carries a second representation
of the same answer:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable or attribute
     - What it is
   * - ``psi_coefficients(psi_dof)``
     - every coefficient of :math:`\psi_h` in the finite-element space — the
       solved potential, not the post-processed :math:`\psi^*`
   * - ``flux_coefficients(flux_dof)``
     - the same for the solved flux :math:`q`, with
       :math:`\bar\nabla\psi = R\,q`
   * - ``fe_collection``, ``fe_order``, ``flux_collection``, ``flux_vdim``,
       ``flux_ordering``
     - the spaces those coefficients live in, spelled as MFEM's
       ``FiniteElementCollection::New()`` takes them back
   * - ``mesh_file``
     - the ``.mesh`` written beside this file. The mesh is **named, not
       embedded**
   * - ``content``
     - ``"total (psi)"`` or ``"remainder (psi - psi_c)"``
   * - ``conductor_kind``, ``conductor_R``, ``conductor_Z``,
       ``conductor_half_width``, ``conductor_half_height``,
       ``conductor_current``, ``mu0``
     - every conductor, enough to rebuild :math:`\psi_c` exactly at any point.
       ``kind`` is 0 for a rectangle of uniform current density and 1 for a
       point filament, whose two half-extents are zero by definition rather
       than by omission
   * - ``input_toml``, ``input_file``
     - the configuration that produced the run, verbatim

.. important::

   **Under a subtracting** ``[conductors] Model`` **the lossy format carries the
   physically exact field and the exact format carries a difference.** That
   inversion is the opposite of what the names suggest and is worth reading
   twice. The gridded variables ``psi``, ``B_R`` and ``B_Z`` *sample*, so they
   can add :math:`\psi_c` back at every node and do; the ``.gf`` files
   *represent*, and :math:`\psi_c` has no representation in the space — for a
   filament it is logarithmic at the conductor.

   This is why the conductor table is here and not a flag. ``content =
   "remainder"`` would tell you that you hold the wrong field without letting
   you fix it: :math:`\psi_c` is recoverable from neither the mesh, nor the
   space, nor the coefficients — only from the conductors.

   A run under a split also writes ``<stem>_psi_total.gf``, holding
   :math:`\psi_h + I_h(\psi_c)` — exact at the nodes and interpolated between
   them. It is for **looking at**, and its name says so.

.. _output-grid-cost:

What the grid format costs
~~~~~~~~~~~~~~~~~~~~~~~~~~

This file is **lossy, by construction**, and it is worth being precise about how,
because the loss does not shrink when the solve improves.

.. important::

   **A structured grid is second order however fine it is.** :math:`\psi_h` is a
   piecewise polynomial of degree :math:`k`; sampling it onto a grid and reading
   it back by bilinear interpolation gives an error that is
   :math:`O(h_{\text{grid}}^2)` — in the **grid** spacing, not the mesh's, and
   with no dependence on :math:`k` at all.

   So the accuracy of a ``.nc`` file is a property of **the file**. The accuracy
   of the ``.gf`` output is a property of **the solve**. Refine the mesh or raise
   the degree and the grid has to be refined quadratically to keep pace; the
   exact route never has to be.

That is the whole reason MEQ writes three formats rather than one, and it is why
a warm start reads ``.gf`` and not ``.nc`` — see :ref:`running-warm-start`. It is
not an argument against the gridded file, which is the right thing for its job:
an outside code has to produce or consume :math:`\psi` on a rectangle and nothing
else, with no finite element library, no mesh format and no agreement about
element types.

.. note::

   A grid node landing exactly on an inter-element face is **ambiguous**, because
   :math:`\psi_h` is discontinuous there and the two elements disagree. Whichever
   element claims it first wins. The disagreement is the size of the face jump,
   :math:`O(h^{k+1})`, so it converges away with everything else — but it is a
   genuine arbitrary choice rather than an averaging rule.

.. _output-band:

The band between :math:`\Gamma_h` and :math:`\Gamma`
----------------------------------------------------

On the curved path the computational domain is the union of background elements
lying *inside* :math:`\Gamma`. So :math:`\Gamma_h` is inscribed, and there is a
band, :math:`O(h)` wide, which is inside the plasma and outside the mesh. A
uniform grid over the plasma has nodes there — a noticeable fraction of them —
and something has to be said about those nodes.

Both output formats deal with it, differently, because they have to.

**The gridded output continues into the band using the flux**, which is the
mixed method paying off somewhere nobody expected. Because :math:`q` is a solved
unknown carried at the same order as :math:`\psi`, and :math:`\gradbar\psi = R
q`, a node :math:`p` outside the mesh can be reached from its foot :math:`x_0` on
:math:`\Gamma_h` by a Taylor step

.. math::

   \psi(p) \;\approx\; \psi(x_0) + r_0 \, q(x_0) \cdot (p - x_0),

in which **nothing is ever evaluated outside an element**. Every quantity on the
right is read at :math:`x_0`, which is on the mesh.

.. warning::

   **The obvious alternative was implemented first and is bounded by nothing.**
   Continuing :math:`\psi_h`'s own polynomial past the edge of its element put
   band nodes at *positive* :math:`\psi` on a test case where the boundary
   condition makes :math:`\psi` exactly zero on :math:`\Gamma` and strictly
   negative inside — so the :math:`\psi = 0` contours in the band were visibly
   wrong. The flux version puts none there, and its band error converges;
   extrapolation does not converge in the band at all.

   **Blending the extrapolation toward the known** :math:`\psi = 0` **on**
   :math:`\Gamma` **does not fix it**, which is worth recording because it looks
   as though it should. Scaling a positive value down by :math:`(1-t)` never
   changes its sign. The error was in *where the field was evaluated*, not in how
   it was weighted.

.. _output-field:

:math:`\mathbf{B}` in the band, and the order it gets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The poloidal field is a pure relabelling of the flux,
:math:`\mathbf{B} = (-q_z, +q_r)`, so it inherits :math:`q`'s order in the
interior. In the band it takes the same Taylor step from the same foot, using the
field's own gradient — which is read *inside* the element, so the guarantee that
nothing is evaluated outside one survives.

.. important::

   **That step does not reach** :math:`\psi`'s **order, and the gap is
   structural rather than a shortcut.** :math:`\psi` is continued with
   :math:`q`, a *solved* variable carrying the potential's own order — that is
   the mixed method paying off. There is no solved variable for
   :math:`\nabla q`: differentiating a discontinuous Galerkin field of degree
   :math:`k` leaves degree :math:`k-1`, so the band continuation of
   :math:`\mathbf{B}` is second order at every :math:`k`.

   The route to better is known and was not taken. The equation being solved
   pins two of :math:`\nabla q`'s four entries exactly, and
   :math:`R q = \gradbar\psi` pins a relation among the others — but what
   remains still has to be differentiated, so they buy *structure* rather than
   an order, at the cost of plumbing the source term into the sampler.

.. warning::

   **Use the** ``extrapolated`` **mask.** It is a strict subset of ``inside``,
   and it names exactly the nodes that were reached by a Taylor step rather than
   evaluated in an element. Those nodes carry real data — they are in the plasma
   — which is precisely what makes them look trustworthy: a check that the mask
   agrees with the non-``NaN`` data passes on them.

   Drop them before computing an error norm, before differencing two runs, and
   before anything else quantitative. Keep them for a picture.

   Any field derived from the geometry must also use **the node's own radius**,
   not the foot's. A quantity that depends on :math:`R` and is evaluated at the
   foot is wrong by the width of the band; for a quantity whose exponent carries
   :math:`R^2`, it is wrong twice over. This was measured, and the difference is
   several orders of magnitude, which is why it is stated as a rule rather than
   left as a detail.

.. _output-flux-surfaces:

The flux-surface format
-----------------------

``[output] FluxSurfaces = true`` writes ``<stem>_surfaces.nc``: the flux
surfaces, sampled at equispaced poloidal angle about the magnetic axis, and the
flux-surface averages over them, against a flux label. It is what a 1-D
transport code reads, and it is the one output that is a *reduction* rather than
a resolution of the answer.

**Off by default**, for two reasons. It costs a contour trace and an angle fit
per surface, which on a coarse mesh is comparable with the solve; and it is the
one output that can be impossible on a run that solved perfectly well — a level
whose surface is not closed, or not star-shaped about the axis, has no
flux-surface average, and MEQ refuses rather than inventing one. A failure there
is reported on standard error and does not change the exit code: the equilibrium
has already been written and is unaffected.

**Under** ``[conductors] Model`` **these are surfaces of the physical flux, not
of the solved remainder.** A subtracting model leaves the solver holding
:math:`\psi_p`, and a flux surface is a level set of :math:`\psi_c + \psi_p`;
the tracer evaluates :math:`\psi_c` in closed form at each point it visits, so
the surfaces, :math:`V'`, the safety factor and the metric are all quantities of
the machine rather than of the remainder. The file records which conductor model
produced it in its ``conductor_model`` attribute, beside ``coils`` and
``coil_current``, for the same reason the gridded format does: a filament set and
a rectangle set of the same currents are different machines, and two files that
do not say which cannot be differenced.

The layout is ``flux × theta`` with ``theta`` fastest, which is C row-major and
the same convention the gridded format uses:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Variable
     - What it is
   * - ``rho(flux)``, ``normalised_flux(flux)``, ``psi(flux)``
     - The label :math:`\rho = \sqrt{\Psi_N}`, the normalised flux
       :math:`\Psi_N`, and the level :math:`\psi` itself.
   * - ``theta(theta)``
     - Geometric poloidal angle about the magnetic axis, radians.
   * - ``R(flux, theta)``, ``Z(flux, theta)``
     - The surface.
   * - ``extrapolated(flux, theta)``
     - 1 where the node's field came from the band extension rather than from an
       element.
   * - ``V_prime(flux)``
     - :math:`V' = \oint 2\pi R \, \mathrm{d}l / |\nabla\psi|`. The
       :math:`2\pi R` is part of the definition, so this is
       :math:`|\mathrm{d}V/\mathrm{d}\psi|` with :math:`V` the enclosed volume.
   * - ``volume(flux)``, ``cross_section_area(flux)``
     - :math:`\oint \pi R^2 \, \mathrm{d}z` and :math:`\oint R \, \mathrm{d}z`,
       by Green's theorem on the same nodes.
   * - ``arc_length(flux)``, ``surface_area(flux)``
     - :math:`\oint \mathrm{d}l` and :math:`\oint 2\pi R \, \mathrm{d}l`.
   * - ``inverse_R_squared(flux)``
     - :math:`\langle R^{-2} \rangle`.
   * - ``grad_psi_squared_over_R_squared(flux)``
     - :math:`\langle |\nabla\psi|^2 / R^2 \rangle`.
   * - ``abs_grad_psi(flux)``, ``grad_psi_squared(flux)``
     - :math:`\langle |\nabla\psi| \rangle` and
       :math:`\langle |\nabla\psi|^2 \rangle`, which are what
       :math:`\langle |\nabla\rho| \rangle` and
       :math:`\langle |\nabla\rho|^2 \rangle` are, times the analytic
       :math:`\mathrm{d}\rho/\mathrm{d}\psi`.
   * - ``safety_factor(flux)``
     - :math:`V' g \langle R^{-2}\rangle / 4\pi^2`. **Present only when
       :math:`g(\psi) = R B_\phi` is known**, which on the driver is exactly
       the :math:`q`-driven route: a ``meq::Source`` carries :math:`g g'` and
       not :math:`g`, so recovering
       :math:`g = \sqrt{g_{\text{edge}}^2 + 2\int g g'\,\mathrm{d}\Psi}`
       from a *prescribed* field needs a constant of integration no key
       supplies, and there the variable is **absent rather than zero** —
       because zero is a real machine. Under ``[source] SafetyFactorFile`` the
       outer Newton solves for the coefficients of :math:`g^2` itself, so the
       column is written and the run can be checked against the target it was
       given. The file then also carries ``toroidal_field_driven``,
       ``g_squared_coefficients`` and ``safety_factor_target``, so a reader can
       tell a solved :math:`g` from a supplied one.
   * - ``band(flux)``, ``worst_residual(flux)``, ``transversality(flux)``
     - Per-surface diagnostics: whether any node is band data, the worst
       :math:`|\psi_h - c|` over the surface, and how close a ray came to being
       tangent to it.

Global attributes carry the magnetic axis, :math:`\psi_{\mathrm{ax}}`,
:math:`\psi_{\mathrm{bnd}}`, the cut, the total number of band nodes, and the
provenance the gridded file carries.

The label is :math:`\rho`, not :math:`\Psi_N`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file's ``flux_label`` attribute says so. :math:`\psi` has a quadratic
maximum at the magnetic axis, so :math:`\Psi_N` behaves like
:math:`(\text{distance})^2` there and a surface's minor radius grows like
:math:`\sqrt{\Psi_N}`. Parametrising by :math:`\Psi_N` therefore puts a
square-root branch point on the axis, and every representation converges
algebraically against it with nothing in the numbers to say why. Measured on a
Solov'ev equilibrium, the same fit against :math:`\rho` rather than
:math:`\Psi_N` is a factor of 204 better in the worst error, with the
conditioning untouched — so it is the branch point and not the algebra.

:math:`\Psi_N` is written beside it, so a consumer with its own normalised-flux
grid need not square anything, and MEQ interpolates the family in :math:`\rho`.

Both ends are cut, for different reasons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``FluxInnerCut`` and ``FluxOuterCut`` default to :math:`\Psi_N \in [0.05,
0.95]`.

At the **inner** end a surface shrinks to a point, :math:`V' \to 0`, and
:math:`\mathrm{d}\rho/\mathrm{d}\psi` diverges like
:math:`1/(2\sqrt{\Psi_N})`. That divergence belongs to the *coordinate* rather
than to the extraction — the product with :math:`2\sqrt{\Psi_N}` settles at
1.148 down to :math:`\Psi_N = 0.005` — but the trace gives out first: MEQ
brackets a level by walking outward along a ray from the axis, and near the axis
the bracket is a fraction of an element wide.

At the **outer** end nothing fails, and that is the point. On a curved boundary
:math:`\Omega_h` is the union of background elements lying *inside*
:math:`\Gamma`, so the outer surfaces cross the band between :math:`\Gamma_h`
and :math:`\Gamma` and the extension answers for those nodes. Measured on a
Miller boundary at two resolutions, the share of each surface that is band data:

.. list-table::
   :header-rows: 1

   * - :math:`\Psi_N`
     - 0.50
     - 0.80
     - 0.90
     - **0.95**
     - 0.98
     - 0.99
     - 0.995
   * - :math:`h = 0.0425`
     - 0%
     - 1%
     - 15%
     - **41%**
     - 80%
     - 94%
     - 98%
   * - :math:`h = 0.0212`
     - 0%
     - 0%
     - 2%
     - **28%**
     - 46%
     - 78%
     - 91%

Every one of those traces closed, every fit converged, no ray stalled, and
:math:`|\psi_h - c|` sat at :math:`2\times10^{-13}` throughout. **The residual
does not distinguish an element from an extension**, because the extension
answers as confidently as an element does. So the cut cannot be discovered by
pushing outward until something breaks; the per-node mask is the only signal
there is, and where to stop is a decision. MEQ's default is close to what
production codes use, and is where the outermost surface is still mostly solved
data at production resolutions. The band is :math:`O(h)`, so the right cut moves
with the mesh — read ``extrapolated`` rather than trusting the default.

.. warning::

   **A query outside the cut is refused, not extrapolated.** Some codes
   extrapolate past the last surface they can trace; MEQ does not. A plausible
   :math:`V'` beyond the plasma boundary is worse than no answer, because
   nothing downstream can tell it from a real one.

``extrapolated`` is a mask, ``extrapolated_nodes`` is a count
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same obligation the gridded format carries, one dimension up, and sharper: a
flux surface can be inside :math:`\Omega_h` at one :math:`\theta` and outside it
at the next, so the per-surface ``band(flux)`` flag under-reports in the middle
of a band excursion — which is exactly where a :math:`q(\psi)` profile is being
read. Drop nodes with ``extrapolated = 1`` before computing an error norm or
differencing two runs; keep them for a picture.

Looking at one
~~~~~~~~~~~~~~

:ref:`plot_equilibrium.py <tools-plot>` reads this file as well as the gridded
one, and tells them apart by their dimensions rather than by the name:

.. code-block:: sh

   tools/plot_equilibrium.py run_surfaces.nc -o family.png

It draws the traced surfaces in the poloidal plane beside a panel for every
average the file carries — which is how a family with a ``safety_factor`` gets a
panel for it and one without gets no empty axes — with the band drawn **per
node** over the surfaces it affects, the affected surfaces ringed on the profile
panels, and the two cut regions shaded so that a family over :math:`\Psi_N \in
[0.05, 0.95]` reads as a cut rather than as an answer that starts somewhere.
``tools/README.md`` has the options.
