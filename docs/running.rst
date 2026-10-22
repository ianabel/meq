Running MEQ
===========

.. code-block:: text

   meq <config.toml>    solve the equilibrium the file describes
   meq --help           usage
   meq --version        the build this is

There are no subcommands, and the flags there are describe the *process*
rather than the equilibrium: ``--device``, ``--mesh-command``, ``--mesh-ready``,
``--version``, the two that print diagnostics, ``--profile`` and
``--trace-solver-stats``, one that sets a diagnostic of the trace solve,
``--trace-refinement-cap`` (see :doc:`linear_solvers`), one that bounds a
:math:`q`-driven run's outer loop, ``--q-iteration-cap`` (see
:ref:`configuration-q-driven`), and the two that ask
for files beside the ``.nc``, ``--gf`` and ``--vtk`` (see
:ref:`running-options`). Everything about a run is
in the configuration file, which
:doc:`configuration` documents key by key, and which ``examples/`` has worked
instances of — see :doc:`examples`.

.. code-block:: sh

   mkdir -p runs/soloviev && cd runs/soloviev
   ../../install/double/bin/meq ../../source/examples/soloviev-nstx.toml

A run is made from a directory under the project's ``runs/`` (see
:ref:`install-layout`), because that is where it writes: every output file, and
a mesh ``meq-run`` generates, lands relative to the working directory. A path
the configuration *reads* is relative to the configuration file instead, so a
shipped example finds its profile tables wherever it is run from — see
:ref:`configuration-paths`.

``--version`` reports four things, and two of them are not MEQ's own choices:

.. code-block:: text

   MEQ 0.1.0 (6a078e2) (MFEM 4.10.1, double, cuda)

the MEQ version and the commit the binary was built from, then the MFEM it
linked, the floating-point precision that MFEM was configured at, and the device
backends this build could be given. **Precision is not a key in the
configuration file and cannot be**: ``MFEM_PRECISION`` decides ``real_t``
throughout the library, so one source tree built against two MFEM installs
produces two binaries with the same name that answer to different round-off
floors. ``--version`` and the output files' ``meq_precision`` attribute are how
a run says which one it was. See :doc:`install`.

.. note::

   The absence of a command line is a deliberate choice, not an unfinished one.
   A configuration file is a record of what was run; a command line is not. When
   somebody asks six months later which equilibrium a plot came from, the answer
   should be a file that can be re-run, and every run's output carries the name
   of the file that produced it as an attribute.

Making the mesh from the same file
----------------------------------

A free-boundary machine needs a mesh ``mfem::Mesh::MakeCartesian2D`` cannot
build: a semicircle reaching :math:`R = 0` exactly, with the conductors
fragmented in. That comes from gmsh, through ``tools/mesh/halfdisc.py``.

``[mesh.generate]`` puts the geometry of that mesh in the same file as the
equilibrium, and ``meq-run`` is the one command that makes it and then solves:

.. code-block:: sh

   ../../install/double/bin/meq-run ../../source/examples/diverted-tokamak-generated.toml

``meq-run`` is a script, configured into the build tree beside the binary it
runs and installed beside it in ``install/<name>/bin``, with the generator under
``install/<name>/share/meq/tools``; the installed copy finds both relative to
itself, so it keeps working after the build tree is deleted. That is a decision rather than a stage: ``meq`` links MFEM and not gmsh —
``tools/mesh/README.md`` has the reasoning — so *something* outside the solver
has to run the generator, and this is that something. It is deliberately thin,
and in particular it does **not** parse the configuration. There is one reader
of MEQ's schema and it is ``meq``:

.. code-block:: text

   meq --mesh-command <config>   the generator's argument list, shell-quoted —
                                 or nothing at all, and exit 0, when the file
                                 names a mesh that already exists
   <the generator>               writes the .msh
   meq --mesh-ready <config>     the solve, exec'd, so the exit code is meq's

Every number in that command line came out of the same parser the solve will
use.

.. code-block:: text

   meq-run [options] <config.toml> [meq options ...]

   --remesh       generate even when the stamp says the mesh is current
   --mesh-only    stop after the mesh; do not solve
   --dry-run      print what would be run, and run nothing
   -v, --verbose  echo each command before running it

Anything else is passed through, so ``meq-run f.toml --device debug`` works, and
a file with no ``[mesh.generate]`` is simply solved.

.. note::

   **The generator command is written beside the mesh, as** ``<mesh>.meq-mesh``,
   and a later run producing the same command skips the meshing. So editing a
   coil re-meshes and editing ``PolynomialDegree`` does not — which is why the
   stamp holds the *command* rather than the configuration file's modification
   time.

.. warning::

   **``meq`` refuses a configuration carrying ``[mesh.generate]`` unless it is
   given ``--mesh-ready``.** MEQ cannot make that mesh, so it cannot check it
   either, and the failure it would otherwise walk into is the quiet kind: edit
   a coil, forget to re-mesh, and the run converges at full order to the machine
   the *previous* mesh described, with every printed number looking exactly as
   it should. ``meq-run`` passes the flag; a person who has just run the
   generator by hand can pass it too, and is then making the claim themselves.

What a run does
---------------

In order:

#. **Reads the configuration.** Every key is validated up front — unknown keys,
   misspelt keys (with a suggestion), out-of-range values and contradictory
   combinations are all reported before any work is done.
#. **Builds the mesh**, either by reading a mesh file or by making a Cartesian
   triangulation of the ``[mesh]`` box, then applying the requested uniform
   refinements.
#. **Builds the source term** from the ``[source]`` table — see :doc:`sources`.
#. **Builds the computational domain.** With no ``[boundary.shape]`` this is the
   mesh itself and :math:`\psi = 0` is imposed on its boundary: the *fitted*
   path. With a shape, MEQ marks the background elements lying inside
   :math:`\Gamma`, extracts them as a submesh, and builds the transfer paths that
   carry the boundary condition from :math:`\Gamma_h` out to the true
   :math:`\Gamma`: the *curved* path, which is what :doc:`curved_boundary` is
   about.
#. **Solves**, by Newton — possibly inside an adaptive loop, and possibly with
   :math:`\psiax` as an additional unknown. See :doc:`nonlinear`.
#. **Post-processes**, building :math:`\psi^\star` — which is the potential
   every output but the restart file carries.
#. **Writes the answer** to ``<stem>.nc`` — and the GLVis and ParaView files
   beside it when ``--gf``/``--vtk`` or ``[output] GridFunctions``/``VTK`` ask,
   and, if ``[output] FluxSurfaces`` is set, reduces it to flux surfaces and
   writes those too. See :doc:`output`.

.. _running-options:

The files beside the ``.nc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   meq --gf  <config.toml>     also <stem>.mesh, _psi.gf, _grad_psi.gf, _psistar.gf
   meq --vtk <config.toml>     also <stem>/<name>.pvd, bent onto Gamma
   meq-convert <stem>.nc       either or both, afterwards, to the same bytes

The ``.nc`` alone is what a run writes unless asked: it carries the mesh and
every coefficient, so it is the restart, and ``meq-convert`` makes the other
formats from it whenever they are wanted. ``--gf`` and ``--vtk`` say the same
as ``[output] GridFunctions = true`` and ``[output] VTK = true``, and either
asking is enough; ``meq-run f.toml --vtk`` passes them through. Ask the *run*
for them when you need what ``meq-convert`` cannot make — an adaptive run's
per-cycle VTK series, or a subtracting run's ``_psi_total.gf``. See
:doc:`output`.

What it prints
--------------

Progress and diagnostics go to standard output; errors and warnings go to
standard error, every one of them prefixed ``meq:`` and a space.

The line that matters most on an ordinary run is the convergence report:

.. code-block:: text

   MEQ: converged in 4 Newton iterations on 512 elements, degree 3

followed by the residual history, which is printed on **every** run and not only
on failure:

.. code-block:: text

      it            ||R||    ||R||/||r_0||    order

The ``order`` column is the observed order of convergence, computed once three
residuals exist and while they are above the round-off floor. It is there
because it is the one diagnostic that can see a wrong Jacobian: a Jacobian with
an error in it still converges to the *right answer* at the *right rate* in
:math:`h`, and the only thing that changes is the path — the observed order
drops from 2 to 1. :doc:`nonlinear` explains why that matters so much here.

An adaptive run additionally prints a table with one row per cycle — elements,
trace degrees of freedom, elements marked, transfer paths widened, the estimator
:math:`\eta` and the Newton iteration count — and says which stopping condition
ended it.

A run with :math:`\psiax` as an unknown prints the converged axis flux and the
residual of the constraint that defines it.

.. _running-exit-codes:

Exit codes
----------

.. list-table::
   :header-rows: 1
   :widths: 10 24 66

   * - Code
     - Name
     - Meaning
   * - ``0``
     - solved
     - The equilibrium was computed and everything asked for was written.
       ``--help`` and ``--version`` also exit 0.
   * - ``1``
     - configuration error
     - The file was unreadable or invalid, a key was wrong, the requested
       geometry could not be built, or the configuration asked for something
       MEQ refuses (see below). Usually nothing was attempted; two of the
       refusals can only be recognised *after* a solve, and those say so and
       write nothing.
   * - ``2``
     - solve failed
     - The configuration was good and the nonlinear iteration did not converge,
       including after the fallback. The residual history is printed, so the
       failure is diagnosable rather than merely reported. A
       :math:`q`-driven run whose outer loop did not converge exits 2 too,
       with no fallback and nothing written (:ref:`configuration-q-driven`).
   * - ``3``
     - output failed
     - The equilibrium was computed and something went wrong writing it — most
       often an output directory that does not exist. MEQ does **not** create
       the output directory.
   * - ``4``
     - targets not met
     - An inverse solve (:ref:`config-inverse`) converged and everything was
       written, but its targets were **not** met: the answer is the best fit
       the free currents allow, with the equilibrium itself exact. It is what
       the machine can do, not the shape that was asked for, and the residual
       of every target is in the output.

The split between 1 and 2 is the useful one for a script: a 1 means the input is
wrong and re-running will not help, a 2 means the input is well formed and the
solver could not do it, which is a case for refining the mesh, raising the
polynomial degree, or supplying an initial guess.

.. important::

   Exit code 2 exists only because MFEM is built with ``MFEM_USE_EXCEPTIONS``.
   Without it, an MFEM error terminates the process, and a failed solve is a
   ``SIGABRT`` rather than a report. See :ref:`install-mfem`.

.. _running-refusals:

What MEQ refuses rather than approximates
-----------------------------------------

Four things the driver could do approximately and declines to. Each exits 1
with an explanation.

``[boundary] Type = "exact"``
   asks for the boundary datum to be the exact solution of the problem. That
   requires a closed form, and the :cpp:class:`meq::Source` interface does not
   carry one — a source is :math:`F` and :math:`\partial F/\partial\psi` and
   nothing else. Closed-form solutions exist in MEQ, but they live in the test
   fixtures, where they are the thing being converged *against*.

a :math:`\psi_{\mathrm{ax}}` that is not the flux at a magnetic axis
   :math:`\psi_{\mathrm{ax}}` is the largest *nodal* value of
   :math:`\psi_h`, and nothing in that definition makes it an axis. Because it
   is what the profiles are normalised by, a wrong one is not a bad number in
   one field — it is a different equilibrium, with the current, the geometry and
   every profile-derived quantity downstream of it. After a converged
   normalised solve the driver locates the axis as a zero of :math:`q_h` and
   checks the normalised flux there, which must be 1. Only a *positive*
   detection refuses: the search is seeded rather than exhaustive and misses
   rather than false-alarms, so finding no interior extremum at all is a
   warning, not a refusal. On a diverted free-boundary run — an X-point in the
   border and ``[boundary.exterior]`` set — it is first treated as an observed
   failure and retried, as described under :ref:`running-failure`; the warning
   is what is left if the retry does no better.

a plasma that contains the symmetry axis
   :math:`\psi(0, z) = 0` exactly, so the normalised flux on the axis is
   :math:`-\psi_{\mathrm{bnd}}/(\psi_{\mathrm{ax}} - \psi_{\mathrm{bnd}})`,
   and the plasma is wherever that is positive. A **tokamak** is a torus about
   :math:`R_0 > 0` and its symmetry axis is in the vacuum, so a support
   containing :math:`R = 0` is toroidal current threading the machine's own
   centre line — not a large error but the wrong topology, and no refinement
   turns it into the equilibrium ``[source]`` asks for.

   **A levitated dipole or a magnetic mirror genuinely reaches the axis, and
   this does not refuse those.** Neither has a toroidal field, so :math:`g`
   vanishes identically in both — which is the same fact twice, since
   :math:`B_\phi = g/R` must be finite on the axis and a plasma reaching
   :math:`R = 0` therefore cannot carry a toroidal field there. So the test is
   :math:`g g' \equiv 0`, checked across a spread of :math:`\Psi` rather than at
   one value, and not a device name. If a dipole or a mirror is what you meant,
   set ``GGPrime`` to zero and the run proceeds.

a source that does not vanish on the symmetry axis
   :math:`F/R` is :math:`\mu_0 j_\phi`, so a domain reaching :math:`R = 0`
   needs :math:`F(0, z) = 0` or the toroidal current density is infinite there.
   :math:`F = \mu_0 R^2 p' + g g'` leaves only :math:`g g'`, and
   :math:`\psi(0, z) = 0` exactly, so the profiles are evaluated at
   :math:`\Psi = -\psi_{\mathrm{bnd}}/(\psi_{\mathrm{ax}} -
   \psi_{\mathrm{bnd}})`. With a fixed boundary that is zero and every profile
   vanishes there; with ``[boundary.limiter]`` it is *negative* — in the vacuum,
   where :math:`g` is constant and :math:`g g'` must be zero, and where an
   unconfined profile extrapolates instead. Set
   ``[source] ConfineToPlasma = true``, which is what says the vacuum carries no
   current, or give a ``GGPrime`` that vanishes for :math:`\Psi \le 0`. This
   one cannot be checked before the solve, because
   :math:`\psi_{\mathrm{bnd}}` is an unknown of it.

**A prescribed** :math:`\psi_{\mathrm{ax}}` **is not checked as an axis.** With
``[source.axis]``, or a mirror's ``[boundary.electrode]``, the label is the flux
at the point the file names, so the
magnetic-axis check above is not asked and the run prints that it was not; a
magnetic mirror has no axis to locate. A mirror's own refusals —
``[source] SupportShape = "box"`` beside a key it has no state for — are made
at parse time and are listed under :ref:`configuration-mirror-support`.

Everything else the driver rejects is a plain configuration error: a shape that
encloses no element, a shape that touches the edge of the mesh box, an
unreadable initial guess. The library adds one refusal of the same character as
the above — handing a normalised source to the ordinary construction path, see
:doc:`normalised_flux`.

In each case the alternative to refusing is not an approximation but a
*different problem* that would converge perfectly well to the wrong
equilibrium, which is the failure mode this codebase guards against hardest.

.. _running-failure:

When the solve does not converge
--------------------------------

The driver runs a **reactive ladder**: plain Newton first, and on *observed*
failure it rebuilds the solver and tries again with Anderson–Picard followed by
Newton, which is a globalisation rather than a different solver. It says so when
it does:

.. code-block:: text

   MEQ: Newton did not converge: ... Retrying with Anderson-Picard ...

**A free-boundary run with a normalised source has a ladder of its own**,
because its unknowns include :math:`\psi_{\mathrm{ax}}`, :math:`\psi_{\mathrm{bnd}}`
and the rest of the border, which Anderson–Picard has nowhere to put. Its second
rung keeps the border and runs Picard in the field block before Newton finishes,
starting from the same state the first attempt started from. Where the X-point is
in the border it also **holds the X-point within a twentieth of**
``[boundary.exterior] Radius`` **of the position the file gives**: a step that
would carry it further is shortened. The rung is reached in two ways — the first
attempt failing, or the first attempt converging to a flux with **no magnetic
axis at all**, an annulus — and in the second case the first answer stands if the
retry does not find an axis either. Where ``[solver] PlasmaSupportSweeps`` is
set, the retry is the first sweep of the support loop and the support is
re-decided after it exactly as on the first attempt, so a retried run reports
``plasma_support_sweeps`` and ``plasma_support_settled`` like any other.

.. code-block:: text

   MEQ: the bordered Newton did not converge: ... Retrying with PICARD in the field block ...
   MEQ: the retry holds the X-point within 0.12 m of its seed ( 1.1000, -0.6000 ): 0.05 of Gamma's radius

Holding the X-point is a choice about *which* equilibrium is found, so MEQ makes
it only on a run that has already failed: applied from the start, the same
restriction changes the answer on machines that converge without it.

**A mirror has no second rung.** Beside ``[source] SupportShape = "box"`` a
failed first attempt is reported as it is and the run exits 2. On the shipped
mirror the rung was measured to converge from starts the first attempt cannot,
but to a different equilibrium from some of them, so a retry there would turn a
failure into a wrong answer. A mirror's labels start at the starting state's
flux at ``[boundary.electrode]`` and ``[boundary.limiter]``, so the start is
moved only by an ``[initialguess]``; ``[source] PsiAxis`` is refused beside the
box.

**Nor does a run driven by** :math:`q(\Psi)` **whose outer loop failed.** The rung would
re-solve one equilibrium at the toroidal field the loop last tried, which is not
the field the run was asked to find, so the run exits 2 and writes nothing. See
:ref:`configuration-q-driven`.

The ladder is deliberately never *predictive*. Nothing is inferred from the
source about which strategy to use, and two plausible cheap indicators — the
size of the reaction term relative to the operator's first eigenvalue, and
whether the first Newton step makes the residual worse — have both been measured
and found not to separate the cases that work from the cases that do not. One of
them turned out to be *anti*-correlated. See :ref:`nonlinear-no-discriminator`.

If a run exits 2, the things worth trying, in the order they are most likely to
help:

#. **Refine.** Most of the hard cases in the literature benchmarks are not stiff
   problems; they are under-resolved ones, and both refinement paths cure them
   independently — halve :math:`h`, or raise ``PolynomialDegree``. Turning on the
   adaptive loop is the systematic version of this.
#. **Give it an initial guess.** Several standard test sources vanish at
   :math:`\psi = 0`, so :math:`\psi \equiv 0` solves the problem exactly and
   Newton, which starts from the boundary data, stops there in zero iterations.
   ``[initialguess]`` is what gets you off that branch. See
   :ref:`sources-trivial-branch`.
#. **Warm start from a nearby answer.** A converged run at another resolution is
   a much better starting point than the boundary datum, and MEQ will
   interpolate it onto the new mesh for you.
#. **Change the nonlinear ordering.** The default is not the more robust of the
   two on a coarse mesh. See :ref:`nonlinear-ordering`.

.. _running-warm-start:

Warm starts
-----------

A run can start from a previously computed answer, and there are two ways it
happens:

* **Explicitly**, from ``[initialguess]``, naming a previous run's ``.nc`` —
  ``Type = "gridfunction"``, ``File = "<stem>.nc"`` — which carries its own
  mesh and says what field it holds, or a stored grid function and the mesh it
  lives on, ``File`` and ``MeshFile``. The two read the same coefficients and
  restart from the same potential, but a ``.nc`` carries the rest of the
  solver's state too. On the mesh it was solved on it restarts from the whole
  iterate — the flux and the trace as well as :math:`\psi_h` — and from the
  border's unknowns as the run left them: :math:`\psiax`, :math:`\psi_{\rm
  bnd}`, the X-point, the plasma-current scale, the free circuit currents and a
  :math:`q`-driven run's :math:`g^2`. So a converged run restarted from its own
  ``.nc`` starts at its answer, and says what it took; a ``.gf`` starts from
  :math:`\psi_h` and the TOML's seeds. If the stored mesh matches the one being solved on, the
  restart is *exact* — every polynomial coefficient, not an interpolation. If
  it does not, the guess is interpolated onto the new mesh, which is the
  ordinary way to restart a run at a different resolution.
* **Automatically**, between cycles of an adaptive run. Each cycle after the
  first begins from the previous cycle's converged answer, interpolated onto the
  refined mesh, with the Dirichlet datum as the fallback at any node the coarse
  mesh does not cover. On a curved boundary there really are such nodes, because
  the computational domain *grows* as it refines.

What a warm start buys depends entirely on the source. For a source whose
:math:`\partial F/\partial\psi` is zero the problem is affine, Newton takes one
step from anywhere, and a warm start is unmeasurable by construction. For a
genuinely nonlinear source it cuts the number of Newton iterations
substantially. It must not change the answer, and the test suite asserts both
halves of that: strictly less work, and an :math:`L^2` error that does not move.

.. _running-restart-routes:

Three routes, and they are not the same problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 34 46

   * - Route
     - When
     - What it carries
   * - **Exact**
     - Same mesh, same degree
     - Every polynomial coefficient. Bitwise resumable, and it reads the stored
       grid function directly.
   * - **Full order**
     - A different mesh or degree, same code
     - The full :math:`k+1` accuracy of the solve it came from. This is
       :cpp:class:`meq::FieldTransfer`, and it is why MEQ takes a dependency on
       GSLIB.
   * - **Interchange**
     - Another code entirely
     - :math:`\psi(R, Z)` on a structured grid, read back by interpolation. That
       is the ``.nc`` file of :ref:`output-netcdf`.

The first two are what ``[initialguess]`` selects between, and it selects
automatically: a stored guess whose mesh matches is taken exactly, and one whose
mesh does not is interpolated.

.. important::

   **The third is not a substitute for the second, and the gap does not close
   with a finer grid.** Bilinear interpolation on a structured grid is second
   order however fine that grid is (:ref:`output-grid-cost`), so restarting a
   :math:`k = 3` solve through a gridded file throws away most of what the
   previous solve computed. It still converges — a guess is only a guess — but
   "warm" is then doing less work than it looks.

   Staying inside the finite element representation is what carries the full
   order, and that is the entire argument for the extra dependency. The gridded
   route keeps its own job, which is genuinely the right one for it: foreign
   input, where no shared mesh format exists to interpolate through.

.. warning::

   **A restart can quietly become a cold start.** A refined domain can reach
   beyond the stored one, and on the curved path it does so by construction,
   since the computational domain *grows* as it refines. Nodes with no data fall
   back to the Dirichlet datum.

   That fallback is deliberate and does not need to be clever — but the count is
   **reported**, because a restart that found no data for most of the domain
   converges like a cold start and nothing else in the output would say so.

.. warning::

   **A warm start interacts badly with a relative convergence tolerance, and the
   better the guess the worse it gets.** A stopping rule of the form
   :math:`\|R\| \le \texttt{rel\_tol} \cdot \|r_0\|` measures :math:`\|r_0\|` at
   the iterate it was handed; start from a converged answer and that target
   shrinks with it, until it falls below the round-off floor and *nothing* can
   meet it. The symptom is a solve that reaches machine precision in two
   iterations and is then reported as a failure at the thirtieth.

   MEQ handles this: it takes the reference residual from the **cold** iterate —
   where the solve would have started with no guess — and converts it into an
   absolute target. ``RelativeTolerance`` keeps exactly the meaning it always
   had, a cold solve is bit-identical, and only a warm one changes: from failing
   to converging. It is recorded here because anyone writing an outer loop
   around a Newton solver will meet the same trap.

.. _tools-plot:

Looking at the answer
---------------------

``tools/plot_equilibrium.py`` reads either of the NetCDF files a run writes. On
the gridded one it draws flux surfaces, the magnetic field, or both:

.. code-block:: sh

   tools/plot_equilibrium.py run.nc                            # both, to a window
   tools/plot_equilibrium.py run.nc --what surfaces -o psi.png
   tools/plot_equilibrium.py run.nc --what field --levels 40 -o b.png

and on ``<stem>_surfaces.nc``, when ``[output] FluxSurfaces`` asked for one, the
traced surfaces in the poloidal plane and the flux-surface averages against
:math:`\rho = \sqrt{\Psi_N}`:

.. code-block:: sh

   tools/plot_equilibrium.py run_surfaces.nc -o family.png
   tools/plot_equilibrium.py run_surfaces.nc --what profiles -o averages.png

The two are told apart by the file's own dimensions rather than by its name, and
a ``--what`` belonging to the other one is refused rather than ignored. See
:ref:`output-flux-surfaces` for what the second file carries, and how the band
and the cut are drawn on it.

It selects a non-interactive backend automatically when writing to a file, and
it captions the figure with the run's provenance — polynomial degree, element
count, Newton iterations, adaptive cycles, final residual — read from the
attributes in the file rather than supplied on the command line.

For anything else, use the format that suits it: GLVis for the exact solution,
ParaView or VisIt for the high-order VTK — both made from the ``.nc`` by
``meq-convert``, or by the run under ``--gf`` and ``--vtk`` — and any NetCDF
reader for the gridded interchange file. :doc:`output` says which is which and
what each one loses.
