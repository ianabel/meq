Post-processing: :math:`\psi^\star`
===================================

HDG methods admit an element-local post-processing that produces a potential
one polynomial degree richer than the one the solver carries, and which
converges one order faster :cite:p:`Stenberg1991`. MEQ computes it on request,
and its error estimator is built on it.

.. code-block:: cpp

   solver.solve();
   solver.postProcess();
   mfem::GridFunction const & psiStar = solver.postProcessedPotential();

What it is
----------

On each element, :math:`\psi^\star \in P_{k+1}(K)` is reconstructed from
:math:`(\psi_h, q_h)` by matching its gradient to the computed flux and pinning
its mean to that of :math:`\psi_h`:

.. math::

   \left(\gradbar\psi^\star, \gradbar v\right)_K
       &= \left(r\, q_h, \gradbar v\right)_K
       \quad \forall\, v \in P_{k+1}(K) \\
   \left(\psi^\star, 1\right)_K &= \left(\psi_h, 1\right)_K

The first line is a **pure Neumann problem**, so it determines
:math:`\psi^\star` only up to a constant on each element. The second line is
what closes it, and it is unconditional — the total flux driving the problem is
already in :math:`H(\mathrm{div})`, so the element's flux balance is satisfied
and the constant is the only thing left undetermined.

Measured, :math:`\psi^\star` converges at :math:`k+2`, and it survives both the
Newton path and the extension path.

.. important::

   **It is what MEQ reports.** Every output except the restart file carries
   :math:`\psi^\star` rather than :math:`\psi_h`, so the driver post-processes
   on every run and not only on adaptive ones. :ref:`output-which-potential`
   says which file is which, and why the restart file is the exception.

Why it is needed rather than merely nice
----------------------------------------

The residual error estimator of :cite:t:`SanchezVizuet2020adaptive` uses
:math:`\psi^\star` in **four of its five terms**. Building those terms on the
raw :math:`\psi_h` instead loses exactly one order at every :math:`k` — which
was measured, and is the defect the pre-modernisation estimator in this project
had. :doc:`adaptivity` says what that costs an adaptive run.

MEQ keeps the degraded version reachable, as
:cpp:enumerator:`meq::ResidualEstimator::Potential::Raw`, precisely so that the
difference stays measured rather than argued about.

.. note::

   The published method post-processes *nonlinearly*, carrying the source terms
   into the local problem. :cite:t:`SanchezVizuet2020adaptive` records that
   dropping them and using plain linear post-processing is also effective in
   the semi-linear case, and that is what MEQ does — the underlying library
   lifts the nonlinear potential integrators as a Jacobian frozen at the
   computed potential, and the result converges at :math:`k+2` regardless.

.. _postprocessing-singular:

A failure mode worth knowing about
----------------------------------

.. warning::

   **A local post-processing can be corrupted per element, silently, and no
   whole-domain norm can see it.**

The local problem above is singular by construction — it is Neumann — and the
mean-value constraint is what regularises it. Treat that constraint as
*optional*, as one implementation once did on the strength of a flag set by the
mere presence of a nonlinear integrator, and the local matrix is factored
singular wherever the flag is set.

The case where that bites is not exotic. Where :math:`\partial F/\partial\psi`
vanishes — which a tabulated profile with a flat segment produces routinely —
the nonlinear term's Jacobian *is* the zero matrix, so the local problem really
is singular there and nowhere else. The corruption is therefore **per element**,
which is exactly what makes it invisible: with a fraction of the domain
affected, the elements there can be badly wrong while the global norm looks
almost right.

The diagnosis that settled it is worth reusing: putting a floor of
:math:`10^{-12}` on :math:`\partial F/\partial\psi` — twelve orders below
anything else in the problem, incapable of moving a solution — restored the
correct answer. That is a singular matrix and nothing else.

.. important::

   MEQ deliberately carries **no runtime check** for this. A solver should not
   stand permanently on guard against its dependency. The state of the defect
   lived in the test suite instead, as a test asserting the *correct* behaviour
   and failing until the library was fixed — at which point it flipped green on
   its own and became the regression that says the driver can rely on
   :math:`\psi^\star`. See :ref:`testing-stance`.

.. note::

   **The defect did not reach the solve, and that was measured rather than
   argued.** The forward solve's local problem takes its potential block from
   the HDG stabilisation on interior faces, a fixed bilinear form whatever the
   source does; the reconstruction builds a *different* local problem, whose
   regularisation is what failed. Both paths were compared in :math:`\psi_h`
   and :math:`q_h` over a range of degrees and meshes and agreed to round-off,
   far below the discretisation error at every point.

The other outputs
-----------------

``postProcess()`` also produces two flux quantities:

:cpp:func:`meq::GradShafranovSolver::postProcessedFlux`
   the enriched flux, in MEQ's sign convention.

:cpp:func:`meq::GradShafranovSolver::totalFlux`
   the total flux including the stabilisation contribution, in **the library's**
   sign convention — deliberately, because the constraint equation it is
   projected through is written that way. See :ref:`sign-conventions`.

A new ``solve()`` invalidates all three;
:cpp:func:`meq::GradShafranovSolver::isPostProcessed` reports whether they are
current.
