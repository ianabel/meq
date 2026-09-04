Testing and verification
========================

.. code-block:: sh

   cd build && ctest --output-on-failure
   ctest -R Soloviev                       # one, or a family
   MKL_NUM_THREADS=1 ./tests/SolovievConvergence     # a binary directly

Every registered test sets ``MKL_NUM_THREADS=1`` for itself and runs with the
source tree as its working directory, because several of them open
``examples/*.toml`` by relative path. Running a binary by hand needs the same
variable set — see the warning in :doc:`install`.

.. _testing-stance:

The testing stance
------------------

Two rules, and the second is the unusual one.

**A test asserts the behaviour that is wanted, and fails until it is there.**
Never the reverse. A test that asserts a known defect *passes while the defect
stands*, and that makes a green suite compatible with a broken solver — which
destroys the only property a suite is for.

**So the suite is expected to be red while a known defect stands**, and that is
the intended signal rather than a broken build. The failing test's message is
the record of the defect, and it must say what would fix it, because a red suite
nobody can action is one people learn to ignore.

This was got wrong once and then vindicated. A test of the post-processed
solution was first written to assert the *corruption* that was then present, and
it passed. Rewritten to assert the correct behaviour it went red, stayed red
until the underlying library was fixed, and then flipped green on its own — and
that flip is what said the solver could go back to using the quantity. A test
asserting the defect would have said nothing at all.

.. _testing-ladder:

Why convergence rates, and not just unit tests
----------------------------------------------

.. warning::

   **A wrong sign convention converges, at the right rate, to the wrong
   function.** An order-of-accuracy study cannot detect it. Only a comparison
   against a closed form can.

That is the specific hazard this discretisation presents, and it is why the
central benchmark is pinned to an exact Solov'ev equilibrium with published
coefficients rather than to a self-convergence study. There is a second,
sharper version of the same problem for the nonlinear solver:

.. important::

   **A wrong Jacobian is invisible to a convergence table.** Newton converges to
   the same *discrete* solution whatever Jacobian carried it there, so every
   error and every convergence rate is unchanged by an error in
   :math:`\partial F/\partial\psi`. What changes is only the path: the observed
   order of the Newton iteration drops from 2 to 1, and it takes more steps.

   Three things can see it, and all three are in the suite: a finite-difference
   check of the *assembled Jacobian* against the assembled residual; an
   assertion on the observed Newton order; and an affine test source, for which
   Newton must finish in exactly one step.

So the ladder is, in order of what it can catch:

#. **Unit tests** (``tests/unit/``, Boost.Test). Configuration parsing, profile
   values *and* derivatives against closed forms,
   :math:`\partial F/\partial\psi` against a finite difference of :math:`F`.
#. **An exact solution**, as an absolute-error regression — the only thing that
   catches a converged answer to the wrong problem.
#. **Convergence rates** (``tests/convergence/``), asserting a *rate* rather than
   a tolerance: :math:`k+1` for :math:`\psi` and for :math:`q`, :math:`k+2` for
   the post-processed :math:`\psi^\star`, over a range of degrees and several
   dyadic refinements.

Each stage of MEQ's development ended at a measured convergence rate, not at "it
runs".

.. _testing-fixtures:

The analytic solutions
----------------------

``tests/analytic/`` holds closed-form solutions, and they form a deliberate
ladder in how the source depends on :math:`\psi` — which is to say, in what each
one demands of a Newton Jacobian.

.. list-table::
   :header-rows: 1
   :widths: 22 22 16 40

   * - Fixture
     - Source
     - :math:`\partial F/\partial\psi`
     - What it catches
   * - Solov'ev
     - constant in :math:`\psi`
     - :math:`0`
     - The discretisation alone. Newton must converge in one step.
   * - McCarthy
     - linear in :math:`\psi`
     - a nonzero constant
     - A Jacobian missing its mass term — either the term is there or it is not.
   * - Manufactured nonlinear
     - :math:`\psi^2`, :math:`e^{-\psi}`
     - varies
     - A Jacobian that is present but wrong.
   * - Similarity exponential
     - :math:`f_0 e^{n\psi}(1 + \varepsilon r^2)`
     - :math:`nF`
     - The same, against an *exact* rather than manufactured solution.

The last two are both nonlinear and are not redundant. The manufactured one
picks a convenient :math:`\psi` and builds :math:`F` to fit it, so its :math:`F`
is not of a form any physical profile produces. The similarity solution
:cite:p:`KaltsasThroumoulopoulos2016` goes the other way: the free function is
chosen and the solution *follows*, so it tests the solver against a nonlinear
equation somebody might actually pose.

Each fixture checks its own transcription rather than trusting it — every one
recomputes :math:`\dstar\psi` by central differences and asserts it against the
source it hands the solver, so a typo in a twelve-term expansion fails the test
rather than converging beautifully to something else.

.. note::

   **The published Solov'ev coefficients are not self-checking**, and were
   wrong twice in this project's own transcription before being settled. Every
   basis function of that solution is :math:`\dstar`-harmonic, so *any*
   coefficients whatsoever leave :math:`F`, :math:`\dstar\psi` and every
   convergence rate exact. Nothing in a rate table can see a wrong one.

   What settles it is an **independent** quantity — here the curvature of the
   model boundary surface the coefficients are supposed to reproduce
   :cite:p:`CerfonFreidberg2010`. The moral generalises past this project:
   *checking a solve against the formula it used cannot detect a misread
   formula*. The first correction verified the constraints and passed, because
   it verified them against its own wrong constants.

Reading a convergence table
---------------------------

Three things about these tables are not obvious and cause false alarms.

**Rates stop improving at high degree or fine mesh.** The source papers report
round-off dominating beyond roughly degree 5 or eight refinements, so a table
that flattens there is behaving correctly.

**A self-convergence study on a rectangle cannot demonstrate** :math:`k+1`, and
the cap has nothing to do with the physics. A right-angled corner puts the
solution in :math:`H^{3-\epsilon}` and its gradient in :math:`H^{2-\epsilon}`,
and no polynomial degree recovers that. It is worse on a polygon approximating a
curved boundary, where the interior angle sets a singular exponent close to 1.
This is precisely the difficulty the curved-boundary technique exists to remove,
and it is worth knowing before anyone designs another fitted-polygon study. The
exact-solution studies are immune, because there the boundary datum *is* the
trace of a smooth solution.

**Unfitted convergence needs a two-tier rate assertion.** :math:`\Omega_h` is
the union of background elements lying inside :math:`\Gamma`, and *which*
elements those are is not a smooth function of :math:`h` — so the error need not
even be monotone in the mesh count, and one dyadic sequence can be worse than
another. That is geometry, not the transfer. The extension tests allow more
slack per pair than across the sequence as a whole.

**Newton is not monotone on a stiff source**, so an order claim is made on the
best triple above the round-off floor, not on every triple. A history that
wanders for several steps and *then* enters its quadratic endgame is the normal
shape here. Taking the best triple of a short or non-monotone history
manufactures orders out of runs that are not converging at all, so the printed
order is bounded on both sides: too low means a broken Jacobian, too high means
an artefact.

Mutation testing
----------------

**Mutation-test a suite you are relying on.** MEQ's profile and source tests
were checked by deliberately introducing defects — a dropped :math:`r^2`,
:math:`p'` replaced by :math:`p`, :math:`\mu_0` on the wrong term, a flipped
sign, an off-by-one in an interval lookup at a knot, a truncated file
round-trip — and confirming each was caught. That is a cheap way to find out
whether a green suite means anything, and it is worth repeating for the
convergence tests, where the risk of a test that passes regardless is higher.

The same technique is what established the claim above that a rate table cannot
see a Jacobian error: perturbing :math:`\partial F/\partial\psi` by a few
percent and re-running an entire convergence study leaves every error and every
rate unchanged to the digits printed.

.. _testing-coverage:

Coverage
--------

.. code-block:: sh

   cmake -B build-cov -DMEQ_ENABLE_COVERAGE=ON && cmake --build build-cov -j4
   cd build-cov && ctest
   gcovr --root .. --filter 'src/meq/' --print-summary

``MEQ_ENABLE_COVERAGE`` is an option rather than a build type on purpose: a
coverage build wants ``-O0`` and a release build wants ``-O3``, and that
difference should be chosen rather than inherited from whatever
``CMAKE_BUILD_TYPE`` happens to hold.

.. note::

   A previous version of this page gave a different reason — that routing
   ``-O0`` through ``CMAKE_BUILD_TYPE`` would drop ``NDEBUG`` and change which
   ``MFEM_ASSERT``\ s are live. **That was wrong**, and is corrected here rather
   than quietly removed, because it is the sort of claim that sends someone to
   build in Debug for nothing. ``MFEM_ASSERT`` is gated on ``MFEM_DEBUG``, not
   on ``NDEBUG``; see the note in :doc:`install`.

.. note::

   **A line-coverage percentage on this codebase is partly a measurement of
   comment density.** This project comments heavily, and different tools
   disagree substantially about which of those lines are executable at all; the
   same file can read very differently to ``gcov`` and to ``gcovr`` on the same
   run. Use one tool's number for the gate and the other's for "is this code
   tested", and do not compare across them.

   The figure with real headroom in it is **branch** coverage, and it is the one
   nobody has looked at.

What CI can and cannot check
----------------------------

**Continuous integration cannot build the solver, and this is structural.** The
MFEM branch MEQ needs is not published on any remote, so a hosted runner cannot
obtain it and caching does not help — there is nothing to fetch. ``find_package``
is therefore not ``REQUIRED``, and CI builds the MFEM-free half: the
configuration parser, the profiles, the sources, the boundary shapes, plus the
naming check and a line-coverage gate.

.. important::

   **Treat a green CI badge as evidence about the configuration and profile
   layers and about nothing else.** The HDG assembly, Newton, the extension
   technique, the estimator and every convergence claim in this documentation
   are exercised only by a local build. A full job exists in the workflow file
   and is disabled until either the branch is published or a self-hosted runner
   exists.

Performance measurements are not tests
--------------------------------------

``tests/performance/`` is built but is deliberately **not** registered with
CTest, and ``tests/performance/scan.sh`` drives it. Everything it reports is a
timing, and this project's standing rule is that a threaded timing on one
machine is a measurement about that machine.

The one thing in it that *is* a test — and it returns non-zero — is the pair of
correctness properties that make its timings mean anything at all: that threaded
assembly reproduces serial assembly bit for bit, and that the alternative trace
solvers reproduce the default one. Those are in the ordinary suite too.

.. warning::

   **Any device timing must synchronise inside the timing loop.** A GPU solver
   queues work on a stream and returns, so a factorisation timed without a
   device synchronise reads several orders of magnitude too fast — and
   *plausibly* so, because "the reordering was reused, of course it is fast" is
   a story that fits. MEQ's harness synchronises for **every** solver, CPU ones
   included, so that the sync can never be the thing that was forgotten.

The naming check
----------------

``ctest -R naming`` runs ``clang-tidy``'s identifier-naming check over the
library as a warnings-as-errors gate. The conventions it enforces are in
:ref:`organization-naming`. It is skipped with a message if ``clang-tidy`` is
not installed.
