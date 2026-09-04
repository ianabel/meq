Running MEQ
===========

.. code-block:: text

   meq <config.toml>    solve the equilibrium the file describes
   meq --help           usage
   meq --version        the build this is

There are no other flags and no subcommands. Everything about a run is in the
configuration file, which :doc:`configuration` documents key by key, and which
``examples/`` has worked instances of — see :doc:`examples`.

.. code-block:: sh

   ./build/meq examples/soloviev-nstx.toml

.. note::

   The absence of a command line is a deliberate choice, not an unfinished one.
   A configuration file is a record of what was run; a command line is not. When
   somebody asks six months later which equilibrium a plot came from, the answer
   should be a file that can be re-run, and every run's output carries the name
   of the file that produced it as an attribute.

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
#. **Writes the answer**, in three formats. See :doc:`output`.

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

      it            ||r||    ||r||/||r_0||    order

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
     - Nothing was attempted. The file was unreadable or invalid, a key was
       wrong, the requested geometry could not be built, or the configuration
       asked for something MEQ refuses (see below).
   * - ``2``
     - solve failed
     - The configuration was good and the nonlinear iteration did not converge,
       including after the fallback. The residual history is printed, so the
       failure is diagnosable rather than merely reported.
   * - ``3``
     - output failed
     - The equilibrium was computed and something went wrong writing it — most
       often an output directory that does not exist. MEQ does **not** create
       the output directory.

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

There is exactly one thing the driver could do approximately and declines to,
and it exits 1 with an explanation:

``[boundary] Type = "exact"``
   asks for the boundary datum to be the exact solution of the problem. That
   requires a closed form, and the :cpp:class:`meq::Source` interface does not
   carry one — a source is :math:`F` and :math:`\partial F/\partial\psi` and
   nothing else. Closed-form solutions exist in MEQ, but they live in the test
   fixtures, where they are the thing being converged *against*.

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

* **Explicitly**, from ``[initialguess]``, naming a stored mesh and grid
  function. If the stored mesh matches the one being solved on, the restart is
  *exact* — every polynomial coefficient, not an interpolation. If it does not,
  the guess is interpolated onto the new mesh, which is the ordinary way to
  restart a run at a different resolution.
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

.. warning::

   **A warm start interacts badly with a relative convergence tolerance, and the
   better the guess the worse it gets.** A stopping rule of the form
   :math:`\|r\| \le \texttt{rel\_tol} \cdot \|r_0\|` measures :math:`\|r_0\|` at
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

``tools/plot_equilibrium.py`` reads the NetCDF output and draws flux surfaces,
the magnetic field, or both:

.. code-block:: sh

   tools/plot_equilibrium.py run.nc                            # both, to a window
   tools/plot_equilibrium.py run.nc --what surfaces -o psi.png
   tools/plot_equilibrium.py run.nc --what field --levels 40 -o b.png

It selects a non-interactive backend automatically when writing to a file, and
it captions the figure with the run's provenance — polynomial degree, element
count, Newton iterations, adaptive cycles, final residual — read from the
attributes in the file rather than supplied on the command line.

For anything else, use the format that suits it: GLVis for the exact solution,
ParaView or VisIt for the high-order VTK, and any NetCDF reader for the gridded
interchange file. :doc:`output` says which is which and what each one loses.
