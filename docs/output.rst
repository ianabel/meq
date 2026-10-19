Output
======

Every run writes **one file**, ``<stem>.nc``, and it is the whole answer:
:math:`\psi^\star` and :math:`\mathbf{B}` sampled on a uniform
:math:`(R, Z)` grid for anything that reads NetCDF, **and** the solve's mesh,
every coefficient of :math:`\psi_h`, :math:`q` and :math:`\psi^\star`, and the
boundary curve, exactly. So it is also the **restart**, and the GLVis and
ParaView files are made from it — by the run when asked, or afterwards by
``meq-convert``, to the same bytes. On request there is one more file, which is
not the equilibrium in another resolution but a **reduction** of it: the flux
surfaces and the integrals over them.

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Files
     - Written
     - What it is
   * - ``<stem>.nc``
     - always
     - **The equilibrium.** The interchange grid, read by anything that reads
       NetCDF, and beside it the exact solution and the mesh it lives on — the
       restart. See :ref:`output-gridded`.
   * - ``<stem>.mesh``, ``<stem>_psi.gf``, ``<stem>_grad_psi.gf``
     - ``--gf``, or ``[output] GridFunctions``
     - **Exact,** in MFEM's own format, for GLVis. Every polynomial coefficient
       of the finite element solution; carries :math:`\psi_h`, and is a
       restart too.
   * - ``<stem>_psistar.gf``
     - the same
     - The **post-processed** potential :math:`\psi^\star`, in
       :math:`P_{k+1}`. See :ref:`output-which-potential`.
   * - ``<stem>/<name>.pvd`` and ``<stem>/Cycle000000/``
     - ``--vtk``, or ``[output] VTK``
     - **The picture,** for ParaView and VisIt. VTK Lagrange cells at degree
       :math:`k+1`, with the mesh boundary bent onto the true :math:`\Gamma`,
       and what each cell is made of.
   * - ``<stem>_surfaces.nc``
     - ``[output] FluxSurfaces``
     - **The reduction.** The flux surfaces themselves and the flux-surface
       averages over them, against a flux label. See
       :ref:`output-flux-surfaces`.

.. important::

   **The GLVis and ParaView files are views of the** ``.nc`` **and are not
   written unless asked for.** On a large mesh they were most of a run's output
   time — the VTK bends the mesh onto :math:`\Gamma` before it draws it — for
   files most runs are never opened in. Ask the run for them with ``--gf`` and
   ``--vtk`` (see :ref:`running-options`), or make them from the ``.nc`` at any
   time:

   .. code-block:: sh

      meq-convert run.nc               # both
      meq-convert --vtk run.nc         # the VTK only
      meq-convert --gf --output elsewhere/run run.nc

   ``meq-convert`` writes exactly the files the run would have written, byte
   for byte — the curved VTK included — with two exceptions it cannot make
   from the ``.nc``: an adaptive run's per-cycle series ``<stem>_cycles/``,
   whose meshes only the run sees, and a subtracting run's
   ``<stem>_psi_total.gf``, which evaluates the conductor field on the mesh.
   Ask the run for those with ``--vtk`` and ``--gf``.

``<stem>`` is the output directory and prefix from the ``[output]`` table. MEQ
does **not** create the output directory; if it does not exist the run exits 3
(see :ref:`running-exit-codes`).

.. note::

   **No equilibrium file carries the flux surfaces**, and that is what
   ``<stem>_surfaces.nc`` is for. A consumer wanting :math:`V'`,
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

   **The exception is the restart field**, ``<stem>_psi.gf`` and the
   ``.nc``'s ``psi_coefficients``, **which keep** :math:`\psi_h` deliberately.
   A restart reads it back into a degree-\ :math:`k` potential space;
   :math:`\psi^\star` lives in :math:`P_{k+1}` and would not fit. Making the
   two consistent would break restart, which is why the asymmetry is commented
   at the write site. The ``.nc`` carries :math:`\psi^\star`'s coefficients
   too, as ``psistar_coefficients``.

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

Written with ``--gf`` or ``[output] GridFunctions = true``, or by
``meq-convert`` from the ``.nc``. MFEM's own mesh and grid-function files,
written at full precision. ``_psi.gf``
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

Written with ``--vtk`` or ``[output] VTK = true``, or by ``meq-convert`` from
the ``.nc``, through ParaView's data collection format, as **Lagrange cells at
the degree of the field being drawn** — which is :math:`k+1`, since what is drawn is
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

   Each frame also carries the cell arrays of :ref:`output-vtk-cells`,
   classified on that cycle's own solve, so the plasma boundary can be watched
   as the mesh refines. The last frame's ``region`` is the answer's
   ``element_region``.

.. note::

   One frame per cycle is easy to get *almost* right, and the near miss is
   invisible. Rebuilding the collection object for each frame puts every cycle
   directory on disk with the correct refined mesh, and leaves the ``.pvd``
   index listing only the last of them — so ParaView opens the file, shows a
   single frame, reports no error, and all the data is present on disk. The
   collection appends to its index and never scans the directory, so it has to
   survive between frames and be rebound to each new mesh.

.. _output-vtk-cells:

What each cell is made of
~~~~~~~~~~~~~~~~~~~~~~~~~

Beside ``psi`` and ``B``, every cell carries the summary of the ``.nc``'s
element record (see :ref:`output-element-materials`):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Array
     - What it is
   * - ``region``
     - the ``.nc``'s ``element_region``: 0 ``vacuum``, 1 ``plasma``, 2
       ``coil``, 3 ``passive``, 4 ``ferromagnetic``, 5 ``plasma_boundary``, 6
       ``conductor_boundary``, 7 ``mixed``
   * - ``plasma_fraction``
     - the plasma's share of the element
   * - ``conductor_material``
     - the material-table row of the conductor holding the largest share of
       the element, or −1 where there is none. A tie goes to the lower row
   * - ``conductor_fraction``
     - that conductor's share, or 0
   * - ``attribute``
     - the element's attribute in the ``.mesh``. ParaView's data collection
       writes it as true cell data

The first four are written as point data, at each cell's own Lagrange points.
**No two cells share a point**, so each value is constant on its cell and
nothing is averaged across a face: colour by ``region`` and the cells come out
flat, and a threshold on ``region = 5`` picks out the cells the plasma boundary
cuts. The test suite reads the files back to check this, on curved cells as
well as straight ones. VTK has no ragged cell arrays, so a cell holding
two conductors names only the larger. The full record is in the ``.nc``.

.. _output-netcdf:

.. _output-gridded:

The gridded format
------------------

:math:`\psi` and :math:`\mathbf{B}` sampled onto a uniform :math:`(R, Z)` grid,
whose extent is :math:`\Gamma`'s bounding box when there is a shape and the
mesh's box otherwise, unless ``[output] GridRMin`` to ``GridZMax`` name it. Node
:math:`i` is at :math:`R_{\min} + (R_{\max} - R_{\min})\,i/(\texttt{GridNR} - 1)`,
so naming another code's grid puts MEQ's samples on that code's nodes, read
from the exact field rather than interpolated. This is the format to give to
anything that is not a finite element code.

Variables: the grid coordinates; ``psi`` (which is :math:`\psi^\star` — see
:ref:`output-which-potential`); ``B_R`` and ``B_Z``; a byte mask ``inside``; a
byte mask ``extrapolated``; ``node_material``, what each node is made of (see
:ref:`output-node-material`); and the boundary curve as a pair of coordinate
arrays. A rotating run additionally carries a density for each species
and the electrostatic potential. Nodes outside the plasma are ``NaN`` *and*
carry ``inside = 0``, so a reader that checks either one is safe.

The global attributes carry the run's provenance — the configuration file name,
the MEQ and MFEM versions, the polynomial degree, the element count, the Newton
iteration count and final residual, the boundary treatment, and for adaptive
runs the cycle count, the final estimator and the marking strategy. This is what
makes a ``.nc`` file self-describing enough to plot six months later.

``meq_version`` names the commit the binary was built from, and it is read at
**build** time rather than at configure time, so a long-lived build directory
cannot stamp a fresh binary with a stale commit. A ``-dirty`` suffix means the
tracked files had uncommitted edits.

Two more say which *build* wrote the file, which is not the same question as
which commit:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - What it is
   * - ``meq_precision``
     - ``"double"`` or ``"single"`` — what ``mfem::real_t`` is in the MFEM this
       binary linked. Not a MEQ option and not a key in the configuration file:
       ``MFEM_PRECISION`` is a cache variable of MFEM's own build, and MEQ
       inherits it through ``MFEM_DIR``. It belongs in the file because the two
       builds answer to floors three orders apart, and an equilibrium's
       accuracy is qualified per problem and per format
   * - ``meq_device_backends``
     - the device backends this build *could* configure — ``"cuda"``,
       ``"hip"``, ``"cuda, hip"`` or ``"none"``. Available is not configured:
       whether a particular run used one is ``--device``, which is process-wide
       runtime state

Both are on ``<stem>_surfaces.nc`` as well.

A normalised run carries ``psi_axis``, and says what it is. With
:math:`\psi_{\mathrm{ax}}` at the located magnetic axis — the default —
``axis_normalised_flux``, ``axis_r`` and ``axis_z`` give the axis the check
found, and :math:`\Psi` there must read 1. With ``[source.axis]``, or a
mirror's ``[boundary.electrode]``, it is a **label** at a prescribed point — a
mirror's :math:`\psi_{\mathrm{inner}}` — and the file carries
``axis_constraint = "prescribed point"`` with ``axis_point_r`` and
``axis_point_z`` instead: no magnetic axis was asked for, so none is reported.
A magnetic mirror's run (``[source] SupportShape = "box"``) also carries:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - What it is
   * - ``plasma_support``
     - ``"box (...)"``: :math:`F = 0` outside :math:`0 < \Psi \le 1` and
       :math:`|z| \le` ``plasma_support_axial_edge``
   * - ``plasma_support_axial_edge``
     - ``[source] AxialEdge``, the prescribed throat
   * - ``plasma_support_nodes``, ``plasma_support_nodes_total``
     - how many potential nodes the converged support admits, of those visited.
       Zero is a converged vacuum, which the run also warns of
   * - ``axial_cutoff_peak``, ``axial_cutoff_value``, ``axial_cutoff_relative``,
       ``axial_cutoff_samples``
     - the cutoff's report: the peak :math:`|F|/R` (:math:`\mu_0 j_\phi`), its
       largest value on :math:`|z| =` ``AxialEdge``, their ratio, and how many
       samples along the cutoff lay in the mesh. A report of how much current
       the prescribed cutoff cuts through, not a verdict

.. _output-midplane:

**A rotating mirror's run** — ``[source] Type = "rotating"`` with
``SupportShape = "box"``, in any gauge — **also carries its true midplane**, on
a dimension ``midplane`` of its own: where each *solved* surface crosses
:math:`z = 0`, and each species' density there. Every gauge states its
``Density`` tables somewhere else — on ``ReferenceRadius``, on the vacuum's
crossing, or on the solved crossing itself — and this is the one statement of
the plasma that does not depend on which.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Variable
     - What it is
   * - ``midplane_Psi(midplane)``
     - :math:`\Psi` at each knot, uniform from 0 (the outer label) to 1 (the
       inner)
   * - ``midplane_R``, ``midplane_dR_dPsi``
     - :math:`R_0(\Psi)`, where the solved surface crosses :math:`z = 0`, and
       its slope :math:`1/(\partial\Psi/\partial R)` — the crossing the
       midplane gauge reads, from the solved :math:`q` integrated along the
       midplane
   * - ``midplane_n_<species>``, ``midplane_dn_<species>_dPsi``
     - :math:`n_s(R_0(\Psi), \Psi)` in :math:`\mathrm{m}^{-3}` and its
       :math:`\Psi`-derivative along the crossing, both from the closure at that
       radius and exact on the solved state

The three columns ``midplane_Psi``, ``midplane_n_<species>`` and
``midplane_dn_<species>_dPsi`` are a ``DensityFile`` table as they stand
(``Psi f df/dPsi``): under ``Gauge = "midplane"`` they reproduce the run's own
tables at the gauge's knots, and between them up to the small shift the knots'
interpolant leaves (:ref:`the midplane gauge <rotation-gauge>`); from a run in
any other gauge they pose **the same plasma** in the midplane gauge, to the same
shift — which is how to state, from a run in the vacuum-midplane or
fixed-radius gauge, the densities at the true midplane. A solved midplane on
which :math:`\psi` is not monotone has no crossing; the run then warns and the
columns are left out.

.. note::

   **A rotating run's** ``n_<species>`` **and** ``e_phi_0`` **are evaluated at
   every located node**, inside the support or not, because they are algebraic
   in :math:`(R, \psi)`. Outside the plasma they are what the closure would
   give there, not a density anyone has. Keep the nodes where
   ``node_material`` is 1, the plasma row, before reading them as a plasma
   (see :ref:`output-node-material`).

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
   * - ``psistar_coefficients(psistar_dof)``, ``psistar_collection``
     - the post-processed :math:`\psi^\star` and its space, which is what the
       grid above samples and what the VTK draws
   * - ``mesh``
     - the mesh the solve ran on, as one string in **MFEM's own mesh format**:
       write it to a file and MFEM or GLVis reads it as a ``.mesh``. It is
       byte for byte the ``<stem>.mesh`` the run writes under ``--gf``, and
       the element index of every ``element_*`` variable below is its
   * - ``shape_major_radius``, ``shape_centre_height``, ``shape_minor_radius``,
       ``shape_elongation``, ``shape_cosines``, ``shape_sines``
     - the ``[boundary.shape]`` curve, as MXH coefficients, when the run had
       one — what the VTK is bent onto
   * - ``mesh_file``
     - the ``.mesh`` written beside this file, present only when the run wrote
       one
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

.. _output-element-materials:

What each element is made of
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``<stem>.nc`` also records, for every element of ``<stem>.mesh``, the fraction
of it occupied by each **material**: vacuum, plasma, each ``[[coils]]`` block,
and two kinds reserved for later, passive structures and iron. The element
index is the element's index in the ``.mesh`` that ``mesh_file`` names.

The classification is the solver's own. Plasma is decided by the test the
source integrates with: the connectivity fill, ``[source] ExcludeAttributes``,
and the support (one-sided or a mirror's box) on :math:`\psi_h + \psi_c`.
Copper is decided by each coil's rectangle, whichever ``[conductors] Model``
carried its current. Both are asked at the solve's own quadrature points, and
a fraction is the quadrature-weighted share of the element's area. The run's
``regions:`` line reports the same record, the VTK carries a summary of it per
cell (:ref:`output-vtk-cells`), and the grid carries the same test at each node
(:ref:`output-node-material`).

**The fractions are a CF contiguous ragged array**, because almost every
element holds one material:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - What it is
   * - ``element_material_count(element)``
     - how many materials the element holds. Its ``sample_dimension`` attribute
       names ``element_entry``
   * - ``element_material(element_entry)``
     - a row of the material table
   * - ``element_fraction(element_entry)``
     - that material's share of the element. An element's fractions sum to 1

Element :math:`e`'s entries are the ``element_material_count[e]`` that follow
those of the elements before it, so a reader takes a running sum of the counts
to find where each element starts. Beside them, one value per element:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - What it is
   * - ``element_region(element)``
     - a byte for plotting, derived from the fractions, with CF
       ``flag_values`` and ``flag_meanings``: ``vacuum``, ``plasma``, ``coil``,
       ``passive`` and ``ferromagnetic`` where one material fills the element;
       ``plasma_boundary`` where plasma and vacuum share it;
       ``conductor_boundary`` where conductors and vacuum share it; ``mixed``
       for anything else
   * - ``element_plasma_fraction(element)``
     - the plasma's share, which is also in the ragged array
   * - ``element_excluded(element)``
     - 1 where ``[source] ExcludeAttributes`` names the element's attribute
   * - ``element_attribute(element)``
     - the element's attribute in the ``.mesh``

**The material table**, on the dimension ``material``. Row 0 is vacuum, row 1
plasma, and the coils follow in ``[[coils]]`` order. A column that does not
apply to a row holds its fill value:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - What it is
   * - ``material_name``
     - the name: ``vacuum``, ``plasma``, or the ``[[coils]]`` block's ``Name``
   * - ``material_kind``
     - 0 vacuum, 1 plasma, 2 coil, 3 passive, 4 ferromagnetic, with CF flags
   * - ``material_circuit``
     - a row of the circuit table, or −1 for none
   * - ``material_current``
     - a coil's resolved current, in amperes. After an inverse solve this is
       the solved current
   * - ``material_resistivity``
     - a passive structure's, in Ω m
   * - ``material_relative_permeability``
     - 1, except for iron. It is the fill (``NaN``) for an iron whose
       permeability is a curve
   * - ``material_permeability_curve``
     - a row of the permeability curve table, or −1 for none

**The circuit table**, on ``circuit``, holds ``circuit_name``,
``circuit_current`` and ``circuit_supplied``. A ``[[circuits]]`` block is one
supplied row. Every coil outside such a block belongs to a single row named
``default``. That row is a grouping and not a supply: its members carry
independent currents, so ``circuit_supplied`` is 0 and ``circuit_current`` is
the fill. Each coil's own current is its ``material_current``. The ``default``
row appears only when some coil is unwired, which is why ``[[circuits]] Name =
"default"`` is refused (:ref:`config-circuits`). A run with no coils has no
circuit table.

**The permeability curve table** describes non-linear iron, which MEQ does not
produce yet. It is on ``permeability_curve``, as a second ragged array:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - What it is
   * - ``permeability_curve_material``
     - the ``material_name`` of the row the curve belongs to
   * - ``permeability_curve_point_count``
     - how many points the curve has, with ``sample_dimension =
       "permeability_curve_point"``
   * - ``permeability_curve_h``, ``permeability_curve_b``
     - :math:`H` in A/m and :math:`B` in T, increasing in :math:`H` from
       :math:`(0, 0)`. The ``interpolation`` attribute on
       ``permeability_curve_b`` names the rule used between the points

.. important::

   **Check that a curve's name matches the row that points to it.**
   ``permeability_curve_material[material_permeability_curve[m]]`` must equal
   ``material_name[m]``. The name is repeated so that an index error can be
   caught. A curve belongs to one material, so two regions of the same steel
   each hold their own copy.

**Global attributes** give the totals: ``region_elements_<region>`` for every
region, ``material_area_<name>`` for every row in m² (the name made into a legal
attribute name, with ``_<row>`` appended where two rows would collide), and
``region_faults``.

.. note::

   **A fault is plasma current inside a conductor** at some quadrature point.
   The run warns, ``region_faults`` counts the elements, and the fractions keep
   it: the point is counted as the conductor. The usual cause is an unconfined
   source whose plasma term reaches the coils. A ``mixed`` element without a
   fault holds plasma and copper at different points. ``[mesh.generate]
   AlignConductors`` splits such cells along the copper's edges. Under a
   subtracting ``[conductors] Model``, ``conductor_boundary`` cells are
   expected: the mesh is not aligned to the copper, which is the point of
   subtracting it.

.. _output-node-material:

And at each grid node
^^^^^^^^^^^^^^^^^^^^^

``node_material(Z, R)`` is an ``int`` on the grid: the row of the material
table at each node, with ``-1`` (its ``_FillValue``) where ``inside`` is 0.
It is the mask to read a rotating run's densities through, and the support a
reader would otherwise have to rebuild from ``psi`` and the labels.

It is decided by the same test as the element record, asked at the node.
Copper is decided by each coil's rectangle. Plasma is decided by the
connectivity fill and ``[source] ExcludeAttributes`` of the element holding
the node, and by the support on ``psi`` as this file writes it
(:math:`\psi^\star + \psi_c`). A node in the band of a curved run is judged
in the element its value was continued from. A point has no fraction, so
nothing here is a boundary or a fault: a node in copper is that coil's row.

On a mirror, the plasma rows agree with :math:`0 < \Psi \le 1` and
:math:`|z| \le` ``plasma_support_axial_edge``, computed from ``psi``,
``psi_axis`` and ``psi_boundary``. ``DriverAcceptance`` checks this at every
node more than one grid spacing from an edge.

The variable is ``node_material`` and not ``material``, because ``material``
is the material table's dimension. In NetCDF, a variable named after a
dimension is read as that dimension's coordinate variable.

The ``.gf`` files carry none of this. They are the exact restart format and are
unchanged.

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

That is the whole reason the ``.nc`` carries the coefficients beside the grid,
and why a warm start reads those coefficients — from the ``.nc`` or a ``.gf`` —
and never the grid; see :ref:`running-warm-start`. It is
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

.. _psi-coil-group:

The ``psi_coil`` group
----------------------

A run that takes its conductors out of the mesh — ``[conductors] Model`` set to
``"subtracted"`` or ``"filament"``, see :ref:`caching-psi-c` — writes one more
thing into the gridded ``.nc``: a netCDF **group** named ``psi_coil`` holding
:math:`\psi_c`, the flux those conductors contribute, and
:math:`q_c = (1/R)\,\bar\nabla\psi_c`.

.. code-block:: text

   group: psi_coil {
     dimensions:
       potential_dof = 18230 ;
       quadrature_point = 45575 ;
       element_offset = 1824 ;
     variables:
       double psi_c_dof( potential_dof ) ;
       double psi_c_quadrature( quadrature_point ) ;
       int quadrature_offset( element_offset ) ;
       double critical_q_c( critical_flux_node, component ) ;
       byte critical_q_c_usable( critical_flux_node ) ;
       double critical_psi_c( critical_potential_node ) ;
     // :format = "meq-conductor-cache" ;
     // :mesh_digest, :conductor_digest, :polynomial_degree, ... ;
   }

**It is not a field on a grid and must not be read as one.** Everything else in
this file is sampled onto a uniform :math:`(R, Z)` lattice and is meant to be
interpolated; this is the opposite. The values sit at exactly the points MEQ
evaluated them at — one per potential degree of freedom, one per source
quadrature point — in MEQ's own ordering, which is meaningless without the mesh
that produced it. There is deliberately no ``R`` or ``Z`` beside them.

``critical_q_c_usable`` is 0 where :math:`q_c` is NaN, which is the symmetry
axis; the screen travels with the values rather than being re-derived from them,
since deciding whether to read a NaN by reading it does not work.

**What it is for** is handing the next run the expensive half of its setup. The
attributes are the identity it is valid for, and MEQ refuses the data outright
if any of them has moved, so the group is safe to carry around and impossible to
apply to the wrong machine. ``[conductors] CacheFile`` writes the same group as
the root of a file of its own when you want it separately from the equilibrium.

If you want :math:`\psi_c` as a *field*, evaluate it: it is a sum of elliptic
integrals over conductors whose geometry is in the configuration file, and
:math:`\psi = \psi_c + \psi_p` holds everywhere, so the gridded ``psi`` in this
file already includes it.
