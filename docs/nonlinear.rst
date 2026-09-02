Solving the nonlinear system
============================

The Grad–Shafranov equation is semi-linear, so something has to iterate. This
page is about what meq does, what the options are, and — since this is where
runs actually fail — how to get a difficult one to converge.

.. _nonlinear-newton:

Newton, and the obligation it creates
-------------------------------------

**meq uses Newton. The papers it implements use Anderson-accelerated Picard.**
This is a deliberate departure with a cost that must be respected everywhere a
source is written.

:cite:t:`SanchezVizuetSolano2019` state the design meq is reversing: they keep
:math:`F` as opaque problem data so that the solver "relies only on the
discretization of the toroidal operator", and pay for it by iterating on *every*
source — even one linear in :math:`\psi`. Newton makes the opposite trade,
putting :math:`\partial F/\partial\psi` into the operator as a mass term on the
potential block. Hence :ref:`sources-jacobian`: **every source must supply its
own derivative, and every profile must be differentiable twice.**

There is independent support for the choice from the free-boundary literature.
:cite:t:`Heumann2015` report that fixed-point iterations "usually suffer from
very slow convergence or even fail to converge, which made researchers move
towards Newton-type methods", and name vertically unstable plasmas as a case
where Picard does not converge at all. :cite:t:`Lackner1976` says the same from
the other end: plain Picard "will converge to the physically trivial solution
:math:`\psi \equiv 0` if admitted by the formulation of the problem" — which is
:ref:`sources-trivial-branch`.

.. code-block:: toml

   [solver]
   NewtonMaxIterations = 20
   NewtonRelativeTolerance = 1.0e-8
   NewtonAbsoluteTolerance = 1.0e-12

What working Newton looks like: a handful of iterations with the residual
roughly squaring each step once it gets going. :cite:t:`Heumann2015` publish a
five-iteration run over twelve orders of magnitude on more than half a million
unknowns; that is the shape to expect. **A history that grinds down linearly
means the Jacobian disagrees with the residual.**

.. _nonlinear-ordering:

The nonlinear ordering
----------------------

This is the largest structural choice in the solver, and it is worth
understanding even if you never change it.

Hybridization eliminates the flux and potential element by element. When
:math:`F` depends on :math:`\psi`, **there are two orders in which linearisation
and elimination can be done**, and they are genuinely different algorithms:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Ordering
     - What it does
   * - ``CondenseThenLinearise``
     - Condense first. Each element's elimination is then *itself a nonlinear
       solve*, one per element per residual evaluation, and the outer unknown is
       the trace alone.
   * - **``NPC``**
     - Linearise the full :math:`(q, \psi, \hat\psi)` system first, then
       hybridize the linear system that results
       :cite:p:`NguyenPeraireCockburn2009nonlinear`. **Every element-local
       operation is one linear solve. This is meq's default.**

The acceptance signal for NPC is that
:cpp:func:`meq::GradShafranovSolver::localNonlinearIterations` reads **exactly
zero** — where under the condensation it reads in the thousands. That is
asserted rather than assumed, and it is the only way to tell the two apart from
outside.

What NPC buys
~~~~~~~~~~~~~

**An NPC iteration is several times cheaper than a condensation iteration**, at
the same or fewer iterations, because the element-local nonlinear solves are
simply gone. Iteration counts understate this badly: the honest currency is wall
clock, and comparing counts across the two orderings compares unlike things.

**When it fails, it fails cheaply.** On a case that defeats both, NPC reaches its
iteration cap far faster, because it is not spending millions of element-local
nonlinear iterations to find out. That makes it the better thing to *try first*
even on cases it loses.

**It simplifies the bordered Newton of** :doc:`normalised_flux`. With
:math:`\psi` an unknown of the system rather than a function of the trace, two
of the three bordered quantities stop being finite differences — one becomes a
single unit entry and one becomes exactly 1 — and the whole apparatus for
keeping element-local solves consistent across a difference goes away.

What NPC costs
~~~~~~~~~~~~~~

.. warning::

   **NPC is not uniformly better, and the expectation that it would be was
   wrong.** On coarse, under-resolved meshes ``CondenseThenLinearise``
   converges on cases where NPC does not. The pattern in the measurements is
   **resolution, not stiffness**: NPC loses where the mesh is too coarse for the
   source and wins where it is not.

   This is not a defect in either. What the condensation is doing is *nonlinear
   elimination*, which is a genuine and well-understood preconditioning
   mechanism: its element-local solves take :math:`(q, \psi)` to the exact
   solution of each element's nonlinear problem at the current trace, so the
   fields never take an unphysical linear excursion. That robustness is exactly
   what NPC gives up in exchange for uniform local work — a trade the method's
   own authors are explicit about.

   MFEM's own test suite reproduces this gap on its own fixtures, which is
   independent evidence that it is a property of the two algorithms rather than
   of meq's discretisation.

``CondenseThenLinearise`` is therefore kept as a **backup**, and it is the one
to reach for on a stiff, under-resolved mesh. It is also the only ordering that
is parallel, and the only one that accepts an :math:`H(\mathrm{div})` flux
space.

.. note::

   **Compare the two orderings on a resolved mesh, or not at all.** A
   disagreement between them at coarse resolution is :ref:`nonlinear-multiple`,
   not a defect. meq's own contract test runs them against each other only where
   both converge in a handful of steps and the discretisation has one solution.

.. _nonlinear-why-hard:

Why the local structure matters, and where the difficulty actually is
----------------------------------------------------------------------

Three Grad–Shafranov codes solve this equation by Newton and report it robust
:cite:p:`Serino2024,Heumann2015`. In a continuous-Galerkin or finite-difference
discretisation :math:`\psi` is one global unknown vector and :math:`F` enters
only the global residual, so Newton linearises once. Hybridization changes that,
under the condensation ordering: the elimination is itself nonlinear, none of
those local solves is globalised, and any one of them failing poisons the whole
residual.

That structural story led to a prediction — that NPC would fix the hard cases —
**and the prediction was falsified.** NPC has no element-local nonlinear
iteration at all, and it loses on precisely those cases. The element-local
solves were never the cause.

.. important::

   **The cause, where it has been identified, is resolution.** Most of the
   published "difficult" benchmark sources are ordinary problems on a resolved
   mesh — a handful of Newton iterations — and *both* refinement paths cure them
   independently: halve :math:`h`, or raise the polynomial degree without
   touching the mesh. That two independent cures work is what identifies
   under-resolution as the cause.

   So the remedy for a production run is **resolution**, which means the
   adaptive loop, rather than a cleverer globalisation. What globalisation buys
   is the **coarse start** an adaptive run necessarily begins from, which is
   exactly where raw Newton fails.

One benchmark family resists this. Where the Jacobian's reaction term has swept
well past the first several eigenvalues of the operator it is added to, the
linearised operator is strongly indefinite and the *continuous* problem is
multi-valued. Refinement is powerless there because there is no discretisation
error to remove. Continuation in the source amplitude does reach it — but see
the note under :ref:`nonlinear-no-discriminator` for why that cannot be offered
through the ``Source`` interface.

.. _nonlinear-globalisation:

The globalisation ladder
------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - ``Globalisation``
     - 
   * - **``None``**
     - Plain Newton, full steps. **The default**, and what every convergence
       rate in the test suite was measured with.
   * - ``PicardThenNewton``
     - Anderson-accelerated Picard walks the iterate into Newton's basin; plain
       Newton takes it from there and supplies the quadratic endgame Picard
       structurally cannot. **The route for a coarse mesh**, and what the driver
       falls back to.
   * - ``AndersonPicard``
     - Anderson-accelerated fixed point — the source papers' own method.
       :math:`F` is evaluated at the previous iterate, so the potential block is
       **linear** and every element-local problem is a linear solve.
   * - ``PicardOnly``
     - The same fixed point without acceleration, for comparison.
   * - ``LineSearch``, ``KinsolNoLineSearch``
     - KINSOL with and without its backtracking line search. Having both means a
       difference between them is attributable to the line search rather than to
       SUNDIALS.

``PicardThenNewton`` reaches cases nothing else in the solver touches *at that
resolution*, finishing in a few Newton steps and agreeing with Picard's own
answer — which is what makes it a handoff rather than a change of problem.

.. note::

   **Picard's job is not to solve the problem.** Stage one stopping short of its
   own tolerance is an expected outcome, not an error, and the handoff swallows
   that deliberately. This is a globalisation, not a two-solver pipeline.

.. warning::

   **Do not replace Picard's tolerance with an iteration budget.** The handoff
   is *not monotone* in Picard effort: budgets that are too small and budgets
   that are large both converge on cases where intermediate ones fail, and one
   tuned on a given mesh will betray you on the next. Picard's own tolerance is
   the trigger that works wherever it is reached.

Two defaults here were set from measurement rather than from the papers, and
both are surprising. **Plain Picard needs damping and Anderson does not**, so
the damping default is 1.0 — right for one path and wrong for the other. And
**the published Anderson depth of 2 fails here** where depth 1 converges, so the
depth default is 1. Whether that is this fixed point's conditioning or the
implementation's is not established, and raising it expecting the published
behaviour will not work.

Should ``PicardThenNewton`` simply be the default? No
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three reasons, in increasing order of how decisive they are.

#. **It costs a factor of a few on everything that already works** — real, but
   the reason you would trade away for robustness if the other two did not
   exist.
#. **It does not cover everything.** The hardest benchmark fails under plain
   Newton and under the handoff alike. A default that is slower and still needs
   a fallback is not a default.
#. **It silently changes which discrete solution you get.** See below.

.. _nonlinear-multiple:

Coarse discretisations of these sources carry more than one solution
--------------------------------------------------------------------

.. warning::

   **This is the finding that most affects how you should read a converged
   run.** On an under-resolved mesh, the *same* discrete system solved by
   different routes — different ordering, different globalisation — converges
   fully, to a tight tolerance, to answers that differ by as much as several per
   cent in the peak flux.

Both iterations are genuinely converged: tightening the tolerance by orders of
magnitude leaves the disagreement bit-identical. They are different fixed points
of the same discrete equations.

Three consequences:

* **Choosing a globalisation is not a performance decision** on a coarse mesh.
  It changes the equilibrium meq reports — and a coarse mesh is precisely what
  an adaptive run starts from. That is why the default stays plain Newton and
  the ladder stays *reactive*: what the ladder is for is a solve that would
  otherwise not happen at all, and paying for it only then keeps the answer on
  the branch plain Newton finds.
* **Do not gate a cross-solver agreement test at a fixed tolerance on one mesh.**
  meq's own test asserts the two things actually entailed: that the *best*
  agreement across a sweep is at round-off, which says the alternative path
  solves meq's problem rather than a neighbouring one, and that the *worst* is
  bounded well below a different problem.
* **Refine before concluding anything from a disagreement.** At adequate
  resolution the routes agree to round-off.

.. _nonlinear-no-discriminator:

There is no cheap discriminator, and the obvious one is anti-correlated
-----------------------------------------------------------------------

The reactive ladder pays for a failure by running Newton to its cap first. The
obvious improvement is to notice earlier. **Two candidate indicators have been
measured and both are useless**, which is why meq has no predictive trigger.

**The size of the reaction term** relative to the operator's first eigenvalue is
computable from a black-box source, since ``dFdPsi`` is mandatory. But cases
converge at values where others fail, giving no threshold — and the ratio needs
the range of :math:`\psi`, which is not known before solving.

**Whether the first Newton step makes the residual worse** is worse than
uninformative: it is *anti*-correlated. The largest first-step blow-up in the
measured set converges comfortably; the cases where the residual came *down* at
step one both run to the cap. And the initial residual itself carries almost no
information, because at the cold iterate it is dominated by the boundary datum
entering through the trace couplings — it is a measurement of the boundary data
and the mesh, not of the difficulty.

.. important::

   **Do not add a predictive trigger without a measurement that separates the
   cases.** Note also that failure is cheap to *observe* under the default
   ordering — see "when it fails, it fails cheaply" above — so the reactive
   ladder is not paying much for its information.

   Continuation in the source amplitude *does* reach the hardest case, with
   adaptive step control. It cannot be offered through the ``Source`` interface,
   and this is a matter of principle rather than of effort: ``Source`` exposes
   :math:`F` and :math:`\partial F/\partial\psi` and nothing else, so **there is
   no amplitude parameter in it and no way to recover one from a black-box**
   :math:`F`. Offering continuation generally would need a source that can
   parameterise itself, which is an interface change to argue on its own merits.

Why a line search does not rescue NPC
-------------------------------------

Worth recording, because backtracking on the full residual is the recommended
globalisation for this method and it was implemented, measured, and removed.

Under NPC the flux and trace rows of the residual are **linear** in
:math:`(q, \psi, \hat\psi)` — those equations carry no :math:`F` — so a Newton
step annihilates them exactly, and *all* of the Newton remainder lives in the
potential row. That structure is common to any problem of this form. What
differs between problems where the line search works and meq's is not the
structure but the **size of the first correction**: where a full step improves
the potential residual, backtracking helps; where it multiplies it by a large
factor, no step length recovers.

Two further findings from that exercise are worth having:

* **A monotone acceptance test is not a line search.** A backtracking rule that
  accepts on any decrease will always find *some* step that passes, because for
  Newton on an :math:`\ell^2` merit the direction is always a descent direction.
  The result is not failure but a **creep** — one per cent a step, forever —
  which disguises its own failure. An Armijo sufficient-decrease condition
  rejects those steps and fails honestly. KINSOL's line search uses one, and
  fails on the same cases, which is how the two were reconciled.
* **What works is fixing where the iterate is, not how far it steps** — which is
  the ladder meq already has. Picard hands NPC a physically sensible state
  instead of a cold one, and the huge first correction never arises.

.. note::

   **Globalising the outer iteration does not globalise the inner ones.** Under
   ``CondenseThenLinearise`` the failure is often an *element-local* solve
   running out of iterations, and a line search chooses how far to move the
   trace — it cannot make the local problems at that trace well posed. NPC
   removes the question entirely by having no local nonlinear solve.

The stopping rule
-----------------

MFEM's Newton solver stops at :math:`\|r\| \le \max(\texttt{rel\_tol}\cdot
\|r_0\|, \texttt{abs\_tol})`, with :math:`\|r_0\|` measured at the iterate it
was handed. That interacts badly with warm starts, and meq works around it — see
the warning in :ref:`running-warm-start`.

.. note::

   **Assert on the residual drop, not on the iteration count.** A test that an
   affine problem converges "in one iteration" is measuring the *stopping rule*,
   not the solver: whether an exact step is also the *last* step depends on where
   round-off over the residual's scale falls relative to the tolerance, which
   moves when anything about the residual's definition changes. Asserting that
   the first step drops the residual by many orders of magnitude is the property
   actually meant, and it is stable.

.. warning::

   **An iteration that ran out of steps has produced a vector, not an
   equilibrium.** :cpp:func:`meq::GradShafranovSolver::solve` throws rather than
   returning a partial answer, and does not recover a solution from a failed
   iteration.

   Relatedly: a ``GradShafranovSolver`` should be considered **unusable after a
   caught exception from the library**. The throw unwinds out of the middle of
   MFEM and leaves objects as the throw found them. meq's own paths construct a
   fresh solver per solve, and the driver rebuilds before retrying.
