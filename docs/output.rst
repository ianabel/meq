Output
======

Every run writes the same equilibrium **three times**, in three formats. That is
not redundancy: no single format is simultaneously exact, portable and
convenient, and each of the three gives up a different one of those.

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

``<stem>`` is the output directory and prefix from the ``[output]`` table. MEQ
does **not** create the output directory; if it does not exist the run exits 3
(see :ref:`running-exit-codes`).

.. note::

   **None of the three carries the flux surfaces.** MEQ can locate the magnetic
   axis, trace the surfaces, take flux-surface averages over them and fit the
   whole family as a map from a disc — see :doc:`flux_surfaces` and
   :doc:`surface_geometry` — but all of that is reached through the C++ library
   and none of it is written to a file or driven from a configuration key. A
   consumer that wants :math:`V'`, :math:`\langle R^{-2}\rangle` or a safety
   factor links against MEQ rather than reading its output.

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
:math:`q = \gradbar\psi / r` — the solved unknown, **not** the magnetic field,
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
unknown carried at the same order as :math:`\psi`, and :math:`\gradbar\psi = r
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
   :math:`r q = \gradbar\psi` pins a relation among the others — but what
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
   not the foot's. A quantity that depends on :math:`r` and is evaluated at the
   foot is wrong by the width of the band; for a quantity whose exponent carries
   :math:`r^2`, it is wrong twice over. This was measured, and the difference is
   several orders of magnitude, which is why it is stated as a rule rather than
   left as a detail.
