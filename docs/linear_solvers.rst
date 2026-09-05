Linear solvers and threading
============================

A hybridized HDG scheme needs exactly two linear solvers: one for the global
face-coupled trace system, one for the small dense per-element systems. MEQ has
a **third**, and it is not an oversight — it is the price of using Newton.

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - System
     - What MEQ uses
     - Notes
   * - Global trace
     - A sparse direct solver, selectable at run time
     - See below. Falls back to GMRES if the build has no direct solver.
   * - Per-element dense
     - Partial-pivot LU, as a :math:`2\times 2` block
     - LU on the flux block, a local Schur complement, LU on that.
   * - Per-element **nonlinear**
     - An element-local Newton
     - **Only under** ``CondenseThenLinearise``. MEQ's default ordering removes
       it entirely; see :ref:`nonlinear-ordering`.

And a fourth thing that is not a fourth solver: when :math:`\psiax` is an
unknown, the bordered system is solved by block elimination against the *same*
factorisation of the trace Jacobian, costing one extra backsolve per Newton step
and nothing else. See :doc:`normalised_flux`.

.. _linear-trace-solver:

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

The trace system is also **not symmetric** on MEQ's headline configuration —
see :ref:`formulation-symmetry` — so an unsymmetric LU is the right kind of
factorisation. On a *fitted* mesh the matrix is symmetric and negative definite,
so a Cholesky factorisation or a symmetric Krylov method on :math:`-A` would
apply there; it would be wrong wherever :math:`\Gamma` is curved, which is what
MEQ is for.

.. note::

   **The symbolic factorisation is reused across Newton steps.** The sparsity
   pattern does not change between iterations, so re-analysing it every time
   throws away a substantial fraction of each step for no numerical gain. It is
   switched on for the Newton and Picard paths and deliberately off for the
   linear path, where the object is factorised once and destroyed.

   MEQ verifies this **by count rather than by clock** —
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

   **With serial assembly, set** ``MKL_NUM_THREADS=1``. MEQ's inner loop
   factorises a great many *small* dense matrices — one or two per element, per
   residual evaluation — and above a block size that depends on your BLAS each of
   those calls pays for a thread fork and a barrier that dwarfs the arithmetic.
   Measured on the development machine, a whole nonlinear solve at :math:`k = 3`
   went from 0.28 s to **107 s** when MKL threads were turned on under serial
   assembly.

   Where the threshold sits is a property of your BLAS, not of MEQ, so **measure
   it on your machine**. Every registered test sets the variable for itself, and
   it is a no-op on a machine without MKL.

.. important::

   **Threaded assembly removes this constraint, and that is the main reason to
   take it.** MKL suppresses its own threading inside an *active* OpenMP parallel
   region. Once the element loop is itself such a region, the element-local dense
   work is nested, MKL runs it sequentially, and ``MKL_NUM_THREADS`` costs it
   nothing — while the trace solve, which runs on the master thread outside any
   parallel region, still gets all of them.

   Measured on a whole nonlinear solve at :math:`k = 3`, ``OMP_NUM_THREADS=8``:

   .. list-table::
      :header-rows: 1
      :widths: 34 33 33

      * -
        - ``MKL_NUM_THREADS=1``
        - ``MKL_NUM_THREADS=8``
      * - Serial assembly
        - 0.2797 s
        - **107.19 s**
      * - Threaded assembly
        - 0.0927 s
        - **0.0834 s**

   A factor of **1285** between the two modes with MKL threads on, and threaded
   assembly is simply immune. This is what makes PARDISO's MKL threads spendable
   at all: earlier versions of this page said MEQ "cannot have both" and that
   resolving it needed either a thread-count scope around the trace solve or a
   batched local factorisation. Neither turned out to be necessary.

   **The recipe, if you want threaded PARDISO**: ``AssemblyMode = "threaded"``,
   ``TraceSolver = "pardiso"``, and ``OMP_NUM_THREADS`` and ``MKL_NUM_THREADS``
   set to the *same* value greater than one.

.. warning::

   **UMFPACK can never take MKL threads**, whatever the assembly mode. Its BLAS
   calls happen in the trace solve, on the master thread, outside any parallel
   region — so they are not nested and get the full thread count, which is
   exactly the fork-and-barrier collapse described above. If
   ``MKL_NUM_THREADS`` is greater than one, use PARDISO.

   And **do not** combine threaded assembly with ``OMP_NUM_THREADS=1`` and MKL
   threads. A one-thread team does not get the nested-region suppression, and the
   element-local work then pays full MKL threading on every call. MEQ warns about
   this combination at startup.

Threaded assembly
~~~~~~~~~~~~~~~~~

:cpp:func:`meq::GradShafranovSolver::setAssemblyMode` threads **every**
element-local loop in the hybridization: the assembly, and — since MFEM threaded
its ``MultNL`` — the residual and the Jacobian, which is to say every step of
the nonlinear solve. It requires MFEM built with both OpenMP and thread safety,
and it *aborts* rather than falling back without them, so MEQ refuses it at the
setter with a clear message instead and falls back to serial when the mode was
merely inherited rather than asked for.

The two modes agree **bit for bit at** ``MKL_NUM_THREADS=1``, which MEQ asserts
as an exact equality rather than a tolerance, on a linear source *and* on a
nonlinear one. The mechanism is specific: element-local arithmetic reassociates
nothing, the scatter into the trace matrix stays serial and in element order, and
where a colouring is used instead it only changes the order in which a face's two
elements accumulate — and :math:`a + b = b + a` exactly.

.. note::

   Above ``MKL_NUM_THREADS=1`` the two modes differ at round-off — about
   1e-15 in :math:`\psi` — because the serial path hands MKL many threads and
   the threaded path hands it one, and a blocked BLAS-3 sums in a different order
   from an unblocked loop. That is arithmetic reassociation inside MKL, not a
   race, and it is why the exactness assertions run at one MKL thread.

.. note::

   **The default is** ``Threaded`` **on a build that supports it, and it was**
   ``Serial`` **until 2026-09-04.** The change is a change of measurement rather
   than of opinion, and the earlier reasoning is worth knowing because it was
   correct at the time.

   An automatic gate on thread availability had been written and removed: a test
   that assembles a small system hundreds of times inside a bordered Newton got
   **1.8× slower** under it, because the library forks a team and buffers *per
   call*, so a caller that assembles once amortises that and a caller that
   assembles repeatedly pays it every time. Mesh size did not separate the two
   cases, so the solver could not know which caller it had.

   **That argument was about the assembly loop, which was then the only loop the
   flag touched.** Now that the residual and the Jacobian are threaded too, a
   bordered Newton is dominated by work the flag *does* parallelise. The same
   test, same machine, re-measured: **3.2 s serial against 1.33 s threaded** —
   2.4× faster where it was 1.8× slower. The heaviest test in the suite went from
   233 s to 138 s.

   At one thread the flag is a wash (1.01), which is what makes it safe as a
   default; the older note recording 0.86× there no longer holds.

There is still a ceiling, and it is the scatter into the trace matrix. That
cannot be threaded — an unfinalized sparse matrix carries one insertion cursor
for the whole matrix, so two threads writing provably disjoint rows still
collide, and the failure is a hang rather than a wrong answer.

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
   *plausibly* so. MEQ's harness synchronises for **every** solver, CPU ones
   included, so that the synchronise can never be the thing that was forgotten.

.. note::

   **Local GPU results here are correctness evidence, not performance
   evidence.** MEQ is built for double precision, and consumer graphics cards
   run double precision at a small fraction of their single-precision rate where
   datacentre parts run it at about half — so a device timing on a workstation
   can invert the conclusion a production part would give.
