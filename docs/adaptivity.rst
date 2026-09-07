Adaptivity
==========

.. code-block:: toml

   [adaptivity]
   Enabled = true
   MaxIterations = 4
   Strategy = "doerfler"
   Theta = 0.6
   TargetError = 1.0e-8

One turn of the loop is **solve → post-process → estimate → mark → refine**. It
stops when the estimated error falls below ``TargetError``, when
``MaxIterations`` solves have been done, or when marking selects nothing — and
it says which.

.. note::

   ``MaxIterations`` counts **solves**, not refinements. Four means at most
   three refinements, and one means a single solve with an estimate printed —
   which is a useful thing to ask for on its own.

The estimator
-------------

MEQ implements the residual estimator of
:cite:t:`SanchezVizuet2020adaptive`, whose analysis for the semi-linear
unfitted case is :cite:t:`SanchezSanchezVizuetSolano2021`. It has five terms,
exposed individually by :cpp:func:`meq::ResidualEstimator::component`:

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Term
     - Symbol
     - Measures
   * - ``Divergence``
     - :math:`\eta_1`
     - The residual of :math:`-\gradbar\cdot q = F/r` on each element.
   * - ``Constitutive``
     - :math:`\eta_2`
     - The residual of :math:`q = \gradbar\psi^\star / r`.
   * - ``FluxJump``
     - :math:`\eta_3`
     - The jump in :math:`q_h` across interior edges.
   * - ``PotentialJump``
     - :math:`\eta_4`
     - The jump in :math:`\psi^\star` across interior edges.
   * - ``TraceMismatch``
     - :math:`\eta_5`
     - :math:`\hat\psi_h` against :math:`\psi^\star` on element boundaries.

**Four of the five are built on the post-processed potential**, which is why
:doc:`postprocessing` is a prerequisite rather than a refinement. The estimator
throws if it is asked to compute before ``postProcess()`` has been called.

.. warning::

   **A ``TargetError`` chosen against a run from before the estimator moved
   onto** :math:`\psi^\star` **is not comparable.** It is a different and
   smaller quantity on the same mesh. ``TargetError`` is absolute, in the
   estimator's own norm, so there is no way to make it portable across such a
   change; re-calibrate it.

.. note::

   The estimator **caches**, keyed on the mesh sequence. It recomputes when the
   mesh changes or when one of its setters is called, and **not** when the
   solver it borrows is solved again on the same mesh — which the mesh sequence
   cannot see. Call ``Reset()`` after a second solve, or build a fresh
   estimator.

.. _adaptivity-boundary:

:math:`\eta_6`, the boundary indicator
--------------------------------------

The five terms above estimate the **interior** discretisation error. On a
free-boundary run — ``[boundary.exterior]`` in :doc:`configuration` — there is a
second error they structurally cannot see: the Gegenbauer coefficients on
:math:`\Gamma` are a **boundary functional**, a transmission integral reached by
extension from :math:`\Gamma_h`, and nothing in :math:`\eta_1 \ldots \eta_5`
measures it.

That is not a defect in :math:`\eta`. :math:`\eta_5` on :math:`\Gamma_h` compares
:math:`\psi^*` against the datum actually imposed, so it correctly reports the
boundary as well resolved *for the interior problem*. Both are true at once, and
left alone the loop refines the interior of a problem whose answer has stopped
moving on its boundary: measured on the half-disc, the elements touching
:math:`\Gamma_h` are 11 % of the mesh and carry **0.00 %** of :math:`\eta^2`.

:math:`\eta_6` is the transmission residual. The coupling imposes only the
**projection onto the retained modes**, so at convergence the pointwise mismatch

.. math::

   d(x) = q_h\cdot\nu(x)
          - \frac{1}{r}\sum_n a_n\,\sigma_n\,C_n(\mu),
   \qquad \sigma_n = \frac{1-n}{\rho_\Gamma}

is orthogonal to :math:`C_2 \ldots C_{N+1}` and is **not zero** — what survives is
the modes past the truncation and the discretisation error. It is accumulated as
:math:`h_e \int |d|^2\,\mathrm{d}\Gamma` per element, the scaling :math:`\eta_3`
uses for a flux jump.

.. important::

   **Summing it into** :math:`\eta` **is not enough on its own.** :math:`\eta_6`
   is small against the interior terms — measured, :math:`8.6\times10^{-4}`
   against an :math:`\eta` of :math:`2.5\times10^{-1}`, a share of
   :math:`\eta^2` around :math:`10^{-5}` — so a Dörfler competition never
   reaches it, and :math:`\Gamma_h` keeps every face it started with.

   This is **not a threshold to lower.** An interior discretisation error and a
   boundary functional are different quantities in different units, and one sum
   over both is a comparison with no meaning however it is weighted. The driver
   therefore marks the boundary term on **its own** distribution and unions the
   two sets, so the loop drives both errors down; :math:`\eta_6` stays *in*
   :math:`\eta` because the **stopping** rule does have to see it.

With the second pass in place, :math:`\Gamma_h` refines 34 → 46 → 57 → 64 faces
over four cycles and the coefficients' error falls
:math:`1.32\times10^{-3} \rightarrow 1.53\times10^{-4}`, against a control that
does not move at all. It is automatic: a run with ``[boundary.exterior]`` and
``[adaptivity]`` gets it, and a run without a coupling has :math:`\eta_6`
identically zero and is unchanged.

.. _adaptivity-eta5:

Two problems with :math:`\eta_5`, and both are recorded because they converge
-----------------------------------------------------------------------------

**The first is an erratum in the published term.** As printed, :math:`\eta_5`
compares :math:`\hat\psi_h \in P_k(e)` against :math:`\psi^\star`, whose trace
has degree :math:`k+1` — which no element of :math:`M_h` can represent.
Orthogonality splits it:

.. math::

   \|\hat\psi - \psi^\star\|^2_e
     = \|\hat\psi - P_M \psi^\star\|^2_e
     + \|(I - P_M)\psi^\star\|^2_e .

The second piece survives even when the *exact* solution and the best possible
trace are substituted. Scaled by the :math:`h_e^{-1}` the estimator carries, it
becomes an :math:`O(h^k)` floor, and measured, the printed term **is** that
floor — so it converges at :math:`k` rather than :math:`k+1` and drags the
total down with it.

Taking the difference *inside* :math:`M_h` restores :math:`k+1` and reproduces
the rates the paper's own table reports. That is
:cpp:enumerator:`meq::ResidualEstimator::TraceComparison::Projected`, and it is
the default; ``Literal`` is kept so the suite keeps measuring the difference.

**The second is specific to the curved path and was nearly fatal.** On
:math:`\Gamma_h` the term was comparing :math:`\psi^\star` against a trace
pinned to *zero*, rather than against the :math:`\varphi_h` actually imposed
there — see :doc:`curved_boundary`. The difference is then the geometric gap,
:math:`O(h)`. Unmitigated, :math:`\eta_5` was orders larger than every other
term and converged at about half an order: the loop would have run, produced
plausible pictures, and **refined the wrong elements**.

The repair is to give the estimator the datum:

.. code-block:: cpp

   auto datum = solver.transferredDatum();
   estimator.setTransferredBoundary( domain.gammaHMarker(), datum.get() );

so that :math:`\eta_5` compares against the condition actually imposed. The
term then converges at :math:`k+1` and is a small fraction of :math:`\eta_1`
rather than orders larger than it. The driver does this automatically.

.. note::

   **The total estimator barely moves under that repair, and that is the point
   rather than a disappointment.** Passing ``nullptr`` instead simply *excludes*
   those faces, which also restores the rate — by deleting the term exactly
   where the geometry error lives. The two agree closely, because a correctly
   evaluated :math:`\eta_5` on :math:`\Gamma_h` is small.

   What changes is that the boundary elements now *have* an :math:`\eta_5`
   contribution instead of none, so the marking can see them. That the adaptive
   loop is otherwise unchanged is what says the repair did not perturb what
   already worked.

Marking
-------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - ``Strategy``
     - ``Theta``
     - Behaviour
   * - ``"doerfler"``
     - fraction in :math:`(0, 1]`
     - Marks the **smallest** set of elements carrying that share of the total
       estimated error. Small ``Theta`` refines very locally; ``Theta``
       approaching 1 approaches uniform refinement.
   * - ``"maximum"``
     - fraction in :math:`[0, 1]`
     - Marks every element with :math:`\eta_K \ge \texttt{Theta} \cdot \max_K
       \eta_K`. **The parameter runs the opposite way**: large ``Theta``
       refines locally, small ``Theta`` approaches uniform.

Dörfler marking is the default and is what the contraction analysis assumes
:cite:p:`CockburnNochettoZhang2015`. The published Grad–Shafranov experiments
used *maximum* marking, for which the analysis is open — so MEQ offers both, and
means it.

.. note::

   Both marking functions take :math:`\eta_K` itself, **not** the squares, so
   that they take the same argument and neither can be fed the wrong one.
   ``markDoerfler`` squares them itself. ``localSquares()`` is available
   separately for anyone assembling their own criterion, because the squares are
   the additive quantity.

.. _adaptivity-companion:

The companion mesh, and why the curved path needs one
-----------------------------------------------------

.. warning::

   **Refining the computational mesh does not move** :math:`\Gamma_h`.

Refine an element of :math:`\Omega_h` and its children are still inside
:math:`\Gamma`, so :math:`\Gamma_h` stays exactly where it was while the local
element size halves. The ratio :math:`\mathrm{dist}(\Gamma_h, \Gamma)/h_{loc}`
therefore **doubles every cycle** — and the transfer of the Dirichlet datum is
only analysed, and only optimal, while that ratio is :math:`O(1)`. A naive
adaptive loop on the curved path silently walks out of the regime the method is
proved in.

:cpp:class:`meq::AdaptiveDomain` implements the companion-mesh update of
:cite:t:`SanchezVizuet2020adaptive`: it keeps the full background mesh beside
the computational one, and when refining it marks not only the elements the
estimator chose but also **every companion element that** :math:`\Gamma`
**cuts and that shares an edge with a marked one**. Some of those children then
fall inside :math:`\Gamma`, so :math:`\Omega_h` *grows* towards it and the ratio
stays bounded. Measured, it does; without the update it doubles per cycle as
predicted.

``refineWithoutCompanion()`` exists to do the wrong thing deliberately — "it is
what the paper says does not work, and it is here to be measured failing rather
than argued about".

.. important::

   **The refinement will look as though it is crowding** :math:`\Gamma`, **and
   that is the companion update doing its job** rather than the estimator
   concentrating there. Reading it the other way sends you looking for a bug in
   the estimator that is not there.
   :cpp:func:`meq::AdaptiveDomain::lastProximityAdditions` reports how many
   elements came from the proximity rule rather than from the indicator, which
   is how to tell the two apart. The driver prints it as the ``wide`` column.

.. note::

   :cpp:func:`meq::AdaptiveDomain::refine` **invalidates**
   :cpp:func:`meq::AdaptiveDomain::computational` — it builds a new submesh — so
   the solver, the transfer path and every grid function living on it must be
   rebuilt. The driver tears down in order for exactly this reason, and a
   warm-started next cycle transfers the previous answer across the boundary
   first.

Warm starting between cycles
----------------------------

Each cycle after the first begins from the previous cycle's converged answer,
interpolated onto the refined mesh by :cpp:class:`meq::FieldTransfer`, with the
Dirichlet datum as the fallback at nodes the coarse mesh does not cover. On the
curved path there really are such nodes, since :math:`\Omega_h` grows as it
refines.

For a source whose :math:`\partial F/\partial\psi` vanishes this changes
nothing — Newton takes one step from anywhere — so the shipped adaptive example
cannot demonstrate it. For a genuinely nonlinear source it cuts the total Newton
work substantially while leaving the answer unmoved, and the test suite asserts
both halves of that. See :ref:`running-warm-start` for the trap that comes with
it.
