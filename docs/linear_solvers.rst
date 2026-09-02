Linear solvers and threading
============================

A hybridized HDG scheme needs exactly two linear solvers: one for the global
face-coupled trace system, one for the small dense per-element systems. meq has
a **third**, and it is not an oversight — it is the price of using Newton.

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - System
     - What meq uses
     - Notes
   * - Global trace
     - A sparse direct solver, selectable at run time
     - See below. Falls back to GMRES if the build has no direct solver.
   * - Per-element dense
     - Partial-pivot LU, as a :math:`2\times 2` block
     - LU on the flux block, a local Schur complement, LU on that.
   * - Per-element **nonlinear**
     - An element-local Newton
     - **Only under** ``CondenseThenLinearise``. meq's default ordering removes
       it entirely; see :ref:`nonlinear-ordering`.

And a fourth thing that is not a fourth solver: when :math:`\psiax` is an
unknown, the bordered system is solved by block elimination against the *same*
factorisation of the trace Jacobian, costing one extra backsolve per Newton step
and nothing else. See :doc:`normalised_flux`.

Choosing a trace solver
-----------------------

.. code-block:: cpp

   if ( meq::GradShafranovSolver::traceSolverAvailable(
            meq::GradShafranovSolver::TraceSolver::Pardiso ) )
       solver.setTraceSolver( meq::GradShafranovSolver::TraceSolver::Pardiso );

.. list-table::
   :header-rows: 1
   :widths: 18 26 56

   * - Backend
     - Needs
     - Notes
   * - **``UMFPack``**
     - SuiteSparse
     - **The default.** The only one present in every build, and what every
       convergence rate in the suite was measured with.
   * - ``Pardiso``
     - oneMKL
     - Measurably faster than UMFPACK on both the factorisation and the
       backsolve, even single-threaded, and it scales with MKL threads. Not the
       default because most builds do not have it and its licence is not
       everybody's to accept.
   * - ``cuDSS``
     - CUDA and cuDSS, and an ``mfem::Device`` configured **before** the solver
       is built
     - Correct — it agrees with UMFPACK to round-off — but see the warning
       below about timing it.

.. important::

   **Choosing among them is never a numerical decision.** All three drive the
   whole solver to the same recovered :math:`\psi` to round-off, and the test
   suite pins them against each other. It is a performance and licensing choice.

Why a direct solver at all
--------------------------

Because at these sizes it wins. A two-dimensional serial problem with a few tens
of thousands of trace degrees of freedom factorises in a fraction of a second.
Algebraic multigrid earns its place in three dimensions or in parallel, and
serial MFEM has none in any case.

The trace system is also **not symmetric** on meq's headline configuration —
see :ref:`formulation-symmetry` — so an unsymmetric LU is the right kind of
factorisation. On a *fitted* mesh the matrix is symmetric and negative definite,
so a Cholesky factorisation or a symmetric Krylov method on :math:`-A` would
apply there; it would be wrong wherever :math:`\Gamma` is curved, which is what
meq is for.

.. note::

   **The symbolic factorisation is reused across Newton steps.** The sparsity
   pattern does not change between iterations, so re-analysing it every time
   throws away a substantial fraction of each step for no numerical gain. It is
   switched on for the Newton and Picard paths and deliberately off for the
   linear path, where the object is factorised once and destroyed.

   meq verifies this **by count rather than by clock** —
   :cpp:func:`meq::GradShafranovSolver::symbolicFactorisations` against
   :cpp:func:`meq::GradShafranovSolver::numericFactorisations` — on purpose: a
   timing would be a measurement about the machine, and this is a measurement
   about the code. It is also the only thing that could notice the reuse
   lapsing, since a lapse costs speed and nothing an error norm could ever see.

.. _linear-threading:

Threading
---------

There are **two** independent knobs, and sweeping them together is actively
misleading.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Drives
   * - ``OMP_NUM_THREADS``
     - The threaded element assembly, if enabled.
   * - ``MKL_NUM_THREADS``
     - The trace solver's BLAS, PARDISO's internals, **and the element-local
       dense factorisations**.

That third entry is the surprise, and it is why the axes must be separated.

.. warning::

   **Set** ``MKL_NUM_THREADS=1``. meq's inner loop factorises a great many
   *small* dense matrices — one or two per element, per residual evaluation —
   and above a block size that depends on your BLAS, each of those calls pays for
   a thread fork and a barrier that dwarfs the arithmetic. Measured on the
   development machine the effect was dramatic and grew rapidly with polynomial
   degree, to the point where a production run at high degree would have been
   unusable and would have looked like a solver problem rather than a threading
   one.

   Where the threshold sits is a property of your BLAS, not of meq, so **measure
   it on your machine**. Every registered test sets the variable for itself, and
   it is a no-op on a machine without MKL.

.. important::

   **This is in genuine tension with the threaded trace solvers**, which want
   MKL threads and are measurably faster with them. ``MKL_NUM_THREADS`` is
   process-wide, so meq cannot have both, and the setting that makes the trace
   solve fast is the setting that makes the element-local factorisations slow.
   Resolving it needs the element-local factorisation to stop going through
   threaded MKL — either a thread-count scope around the trace solve, or a
   batched local factorisation. Neither is done, and it is the largest single
   performance item outstanding.

Threaded assembly
~~~~~~~~~~~~~~~~~

:cpp:func:`meq::GradShafranovSolver::setAssemblyMode` threads the element-local
half of the hybridization's assembly. It requires MFEM built with both OpenMP
and thread safety, and it *aborts* rather than falling back without them — so
meq refuses it at the setter with a clear message instead.

The two modes agree **bit for bit**, which meq asserts as an exact equality
rather than a tolerance: the mechanism is specific, in that element-local
arithmetic reassociates nothing and the scatter stays serial and in element
order. The scatter is also the ceiling, and it cannot be threaded — an
unfinalized sparse matrix carries one insertion cursor for the whole matrix, so
two threads writing provably disjoint rows still collide, and the failure is a
hang rather than a wrong answer.

.. warning::

   **The default is** ``Serial``, **unconditionally, and an automatic gate on
   thread availability was written and removed.** This is where an isolated
   benchmark most misled: the library forks a team and buffers *per call*, so a
   caller that assembles once amortises that and a caller that assembles
   hundreds of times inside a bordered Newton pays it every time. Under the
   gate, a test that assembles repeatedly on small meshes got substantially
   *slower*.

   **Mesh size does not separate the two cases** — the meshes in question were
   in the range where the isolated benchmark still showed a win. The solver
   cannot know which kind of caller it has.

   So ``Threaded`` is an informed opt-in: **take it for a large mesh assembled a
   few times; leave it for a small one assembled hundreds of times.** Measure it
   on your own workload.

Measuring performance
---------------------

``tests/performance/`` holds the harness and ``scan.sh`` drives it, one process
per point — because MKL fixes its threading at first use, so an in-process sweep
would measure the first setting several times over.

It is deliberately **not** a registered test. Everything it reports is a timing,
and a threaded timing on one machine is a measurement about that machine.

.. warning::

   **A GPU device timing must synchronise inside the timing loop.** A device
   solver queues work on a stream and returns; timed without a synchronise, a
   large factorisation reads several orders of magnitude too fast — and
   *plausibly* so. meq's harness synchronises for **every** solver, CPU ones
   included, so that the synchronise can never be the thing that was forgotten.

.. note::

   **Local GPU results here are correctness evidence, not performance
   evidence.** meq is built for double precision, and consumer graphics cards
   run double precision at a small fraction of their single-precision rate where
   datacentre parts run it at about half — so a device timing on a workstation
   can invert the conclusion a production part would give.
