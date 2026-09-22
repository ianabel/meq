Normalised flux, and :math:`\psiax` as an unknown
=================================================

Equilibrium codes pose profiles against the **normalised** flux

.. math::

   \Psi = \frac{\psi - \psibnd}{\psiax - \psibnd},

which runs over :math:`[0, 1]` from the magnetic axis to the boundary
:cite:p:`GourdainLeboeuf2004`. MEQ solves the fixed-boundary problem with
:math:`\psi = 0` on :math:`\Gamma`, so :math:`\psibnd = 0` and
:math:`\Psi = \psi/\psiax`.

.. code-block:: toml

   [source]
   Type = "rotating"
   Normalised = true
   PsiAxis = 0.08          # a Newton GUESS, not a scale factor

.. important::

   **This makes** :math:`\psiax` **an unknown of the nonlinear system**, not an
   input. It is a global functional of the solution, so fixing it does not
   approximate the problem — it replaces it with a different one.

Why a fixed :math:`\psiax` is a different problem
-------------------------------------------------

Three cheaper answers were tried before the one that works, and each failed in
an instructive way.

**Fixing** :math:`\psiax` **leaves the profile inert.** With the wrong axis
flux, the solution never reaches the part of :math:`[0, 1]` the profile is
shaped over. Measured on a peaked pressure profile, the computed
:math:`\psi` was identical to every printed digit across a factor of five
hundred in pressure amplitude — the profile was contributing essentially
nothing. A fixed :math:`\psiax` is not a simplification of a normalised
profile; unless the value happens to be right it is a **different problem**, and
MEQ keeps a test asserting exactly that as the control.

**A consequence: "is this case stiff?" becomes unanswerable.** Every
configuration converged in one or two Newton steps across an enormous range of
reaction strength, because no nonlinearity was ever switched on.

**Closing the loop outside the solver does not rescue it, and the reason is a
pole.** Iterating :math:`\psiax \leftarrow \max\psi` with relaxation *does*
converge — to a degenerate fixed point where :math:`\psi` and :math:`\psiax`
shrink together, :math:`\Psi` stays :math:`O(1)`, and the pressure gradient runs
away as the solution it drives goes to zero. Mapped out afterwards, the cause is
plain: the outer map has a **pole immediately beside its own fixed point**.
Relaxing harder is not a fix for a pole.

The bordered Newton
-------------------

What works is putting :math:`\psiax` inside the residual. The system closed by
Newton is

.. math::

   R(\lambda, s) &= 0 && \text{the hybridized trace residual, source normalised by } s \\
   G(\lambda, s) &= s - \max \psi_h(\lambda, s) = 0 && \text{the normalisation, as an equation}

in the pair :math:`(\lambda, s)` — the trace, and one scalar. The Jacobian is
**bordered**:

.. math::

   \begin{bmatrix} A & c \\ b^{\mathsf T} & d \end{bmatrix},
   \qquad
   \begin{aligned}
   A &= \partial R/\partial\lambda && \text{the existing hybridized Jacobian} \\
   c &= \partial R/\partial s && \text{dense} \\
   b &= -\partial(\max\psi_h)/\partial\lambda && \text{sparse} \\
   d &= 1 - \partial(\max\psi_h)/\partial s
   \end{aligned}

:math:`c` and :math:`b` are **the non-local terms**, in exactly the sense
:cite:t:`Heumann2015` mean when they warn that a normalised profile "leads to
non-local entries in the stiffness matrix": :math:`\psiax` is a functional of
the whole solution, so perturbing the trace near the magnetic axis moves the
source *everywhere*.

.. note::

   **Why it cannot be a rank-one update inside the element blocks**, which is
   what a continuous-Galerkin code would do. In an :math:`H^1` discretisation
   :math:`\psiax` is one entry of the global unknown and the Jacobian simply
   acquires a rank-one term. Hybridization eliminates flux and potential
   **element by element**, and a term coupling every element to the one element
   holding the axis is precisely what that elimination cannot represent.

   The border is where it goes instead, and it costs **one factorisation and
   two backsolves** rather than a second matrix: solve :math:`A y = R` and
   :math:`A z = c`, then :math:`\delta s = (b\cdot y - G)/(d - b\cdot z)` and
   :math:`\delta\lambda = -y - z\,\delta s`. Assembling the border into an
   :math:`(n+1)`-square matrix would put a dense row and column into the
   factorisation for no gain.

The asymmetry between :math:`c` and :math:`b` is structural: :math:`s` enters
every element's source, so :math:`\partial R/\partial s` has an entry on every
trace degree of freedom — but it is one central difference in a *scalar*.
:math:`\max\psi_h` is one nodal value in one element, and under hybridization
that element's recovered potential depends only on the trace degrees of freedom
of its own faces, so :math:`b` has a handful of entries and the rest are exactly
zero. That locality is measured, not assumed.

.. note::

   :math:`\psiax` **is the flux at the located magnetic axis** — the point where
   the flux :math:`q_h` vanishes — and that is a definition rather than an
   approximation. :cpp:func:`meq::GradShafranovSolver::setAxisConstraint`
   selects it; ``AxisConstraint::LocatedAxis`` is the default.
   ``AxisConstraint::NodalMaximum`` takes the largest nodal value of
   :math:`\psi_h` instead, and is kept as a control: the two differ by
   :math:`O(h^2)` in value and :math:`O(h)` in position, **both independent of**
   :math:`k`, so on a refined high-order mesh they separate rather than converge.

**The located axis costs nothing in the Jacobian, which is the envelope
theorem.** :math:`x^*` moves with the solution, so the chain rule gives

.. math::

   \frac{dG}{d\lambda}
     = -\left[ \left.\frac{\partial\psi_h}{\partial\lambda}\right|_{x^*}
       + \nabla\psi_h(x^*)\cdot\frac{\partial x^*}{\partial\lambda} \right]

and :math:`\bar\nabla\psi = R\,q`, so :math:`\nabla\psi_h(x^*) = 0` at a zero
of :math:`q_h` **identically**. The position term vanishes, no sensitivity of the
root find is needed, and the row is the potential shape functions of
:math:`x^*`'s element evaluated at :math:`x^*` — exact, undifferenced, one
element.

Under the default nonlinear ordering — see :ref:`nonlinear-ordering` — none of
the bordered quantities is a finite difference: with :math:`\psi` an unknown of
the system the row is that shape-function row (:math:`-e_j` in the special case
where the axis lands on a node) and :math:`d` is exactly 1. The located-axis
constraint **requires** that ordering and is refused under the condensation,
where :math:`\psi` is a function of the trace and the row would have to be
differenced with a root find inside every difference.

The search is **warm started** from the previous iterate's axis, which is what
keeps it affordable and what keeps the iteration following one axis rather than
re-running a competition between every extremum at every step. Seeded on the
answer it roots a single element; a whole element away, nineteen of two thousand.
The cost is bounded by a ring count rather than by the mesh.

.. _normalised-analytic-column:

The column in closed form
~~~~~~~~~~~~~~~~~~~~~~~~~

That leaves :math:`c = \partial R/\partial s` as the only differenced quantity
in the border, and it is available analytically. :math:`s` reaches the residual
only through the source, and for
:math:`F = g(\Psi)/\sigma` with :math:`\sigma = \psiax - \psi_{\text{bnd}}`
and :math:`\Psi = (\psi - \psi_{\text{bnd}})/\sigma`,

.. math::

   \frac{\partial F}{\partial \psiax}
       &= -\frac{g'(\Psi)\,\Psi + g(\Psi)}{\sigma^2}, \\
   \frac{\partial F}{\partial \psi_{\text{bnd}}}
       &= \frac{g'(\Psi)\,(\Psi - 1) + g(\Psi)}{\sigma^2},

where :math:`g(\Psi) = \mu_0 R^2 p'(\Psi) + gg'(\Psi)` is what the source
already evaluates and :math:`g'(\Psi)` is one further derivative of each stored
profile — the level :cpp:func:`meq::Profile::doublePrime` supplies.
:cpp:func:`meq::NormalisedSource::normalisationDerivatives` is the interface and
:cpp:func:`meq::GradShafranovSolver::setBorderColumn` chooses between the
assembled column and the differenced one, defaulting to assembled and falling
back silently for a source that does not supply the derivatives.

The formulae agree with a Richardson-extrapolated difference of the source to
:math:`10^{-12}`, and end to end the two routes reach the same
:math:`\psiax` and :math:`\psi_{\text{bnd}}` to ten digits while the assembled
one finishes at a residual some **300 times lower** — the difference's own
accuracy was the floor.

.. note::

   **It is not a universal speed-up, and the differenced route is kept rather
   than removed.** On a well-conditioned border the two are indistinguishable —
   same iteration count, same residual to every digit. The assembled column
   earns its place where the *difference* is poor: a stiff border, and
   structurally wherever a moving support makes the two evaluations straddle the
   plasma edge. Keeping both is how the difference between them can be measured
   rather than assumed.

.. important::

   **Outside the plasma both derivatives are exactly zero, and that is the
   reason to prefer the closed form over a difference rather than a mere
   efficiency.** When the support moves with the solution
   (``ConfineToPlasma``, see :doc:`configuration`), perturbing the normalisation
   *moves the edge*, so a difference evaluates the two sides at supports that do not
   coincide and straddles a kink instead of measuring a derivative. A closed
   form is evaluated pointwise at the current state and knows the point is
   outside.

   The precondition is the one :cpp:func:`meq::NormalisedSource::setPlasmaSupport`
   already carries: with a moving support the true derivative picks up a surface
   term where the edge sweeps, and that term vanishes exactly when the profiles
   vanish at the edge.

.. note::

   **The magnetic axis needs no position derivative, and that is a theorem
   rather than an approximation.** :math:`\psiax` is a *stationary* value, so if
   the axis position :math:`(r_\ast, z_\ast)` is treated as moving with the
   solution, the chain rule gives

   .. math::

      \frac{\mathrm{d}}{\mathrm{d}\lambda}
        \psi\big(\lambda; r_\ast(\lambda), z_\ast(\lambda)\big)
        = \frac{\partial \psi}{\partial \lambda}
        + \nabla\psi \cdot \frac{\partial (r_\ast, z_\ast)}{\partial\lambda},

   and :math:`\nabla\psi = 0` at an interior extremum. The second term vanishes
   identically, so taking the largest *nodal* value loses nothing the Jacobian
   would have used.

   That argument does **not** extend to an X-point, where the constraint is
   :math:`q = 0` rather than a stationary value, so the corner block
   :math:`\nabla q` has to be written down.  It does **not** have to be
   approximated: ``[boundary.xpoint]`` writes it in closed form, and it is
   exact.  :math:`\nabla q` is an order down as an approximation of the
   *continuous* Hessian, which is a real wall elsewhere — see the band
   continuation in :doc:`output` — and what a Jacobian needs is the derivative
   of the *discrete* residual, where :math:`q_h` is a polynomial on its element.
   Approximating a continuous object and differentiating a discrete one are
   different questions.

.. _normalised-guess:

The guess is part of the problem statement
------------------------------------------

.. warning::

   ``PsiAxis`` **is where the Newton iteration starts, not what the answer will
   be.** Reading it as the axis flux of the computed equilibrium is reading the
   starting point of an iteration. The converged value is printed and written to
   the output file as an attribute.

And it is not optional, because **at a fixed** :math:`\psiax` **this equation
has two solutions**: a small positive one and a large one. Only the large one
can satisfy :math:`\max\psi = \psiax`, so the constraint takes the small branch
out of the *solution set* — but not out of the *iteration's reach*, and Newton
starting from the Dirichlet datum walks straight onto it.

A dimensional estimate is enough to pick a starting value of the right size: for
:math:`p \sim A\Psi^\nu` on a box with first Dirichlet eigenvalue
:math:`\lambda_1`, :math:`\psiax \sim \sqrt{\nu A/\lambda_1}`.

.. note::

   This is **not** the trivial-branch trap of :ref:`sources-trivial-branch`. A
   normalised source of this kind does not vanish at :math:`\psi = 0` at all;
   the small solution is a genuine second root, not the zero function.

Globalisation is refused on this path
-------------------------------------

Loudly, and at the driver level too. Every globalisation MEQ has either drives a
residual of its own (the KINSOL ones) or builds no Jacobian to border (the
Picard ones), so there is nothing for the border to attach to. If a normalised
run does not converge, the levers are:

#. **A better** ``PsiAxis``. See the dimensional estimate above.
#. **An** ``[initialguess]`` that puts :math:`\psi` near the right size.
#. **Resolution.**

.. warning::

   **A backtracking line search on the border itself is not optional**, and MEQ
   applies one internally. The full step converges for mild profiles and wanders
   for peaked ones, with the augmented residual climbing through several orders
   of magnitude before :math:`\psiax` crosses zero and the source refuses the
   normalisation. That is not the Jacobian being wrong — the same Jacobian
   finishes the milder cases at observed order 2. It is the equilibrium being a
   **mountain-pass solution of a superlinear problem**, where the linearised
   operator is indefinite and an undamped step leaves the basin.

Reading the residual history
----------------------------

On this path the printed residual is :math:`\|(R, \gamma G)\|` over the full
augmented system, not the trace residual alone. :math:`G` is a flux and
:math:`R` is a trace residual, so the two cannot simply be concatenated;
:math:`\gamma` is the factor that converts a perturbation of :math:`\psiax` into
the units :math:`R` is measured in, and it is **frozen at the first iterate**.

Freezing it keeps the history a comparison of like with like. A :math:`\gamma`
recomputed each step would put the Jacobian's own variation into the convergence
history and manufacture orders out of it.

:cpp:func:`meq::GradShafranovSolver::normalisationResidual` reports :math:`G` on
its own.

.. _normalised-axis-check:

What the axis check is for
--------------------------

Constraining :math:`\psiax` at a located axis removes the failure that
motivated this check: a border on the largest nodal value is satisfied at machine
zero by a spurious nodal spike exactly as it is by an axis, so a run could
converge, report its constraint at ``0.000e+00``, deliver a prescribed current to
seven figures, and describe an equilibrium nobody asked for.

MEQ still locates the axis independently after the solve, as a zero of the
**flux** :math:`q_h` — a solved variable carrying the potential's own order
rather than a derivative of one — and reports the normalised flux there:

.. code-block:: text

   psi_ax = 1.039325e-01 Wb/rad, constraint psi_ax - max psi_h = -5.551e-17
   psi_ax is constrained at the located magnetic axis, a zero of q_h
   the axis, as a zero of q_h: psi = 1.039269e-01 at ( 1.0919, -0.0000 ),
        normalised flux 0.9999

**That last number must be 1.** :math:`\Psi` at the magnetic axis is 1 by
definition when :math:`\psiax` is the axis flux, so the reading is a direct
statement about the quantity the profiles consume. Under the default constraint
it is very nearly 1 by construction, and the check is correspondingly weaker —
what it still catches is a run with **no** interior extremum at all, and one
where the extremum carrying the largest :math:`\Psi` is not the one the
constraint followed. A reading materially below 1 means the profiles were
evaluated over a range the plasma never reaches, and MEQ **refuses**: on a
normalised run :math:`\psiax` is what the profiles are divided by, so a wrong one
is a different equilibrium rather than a bad number.

.. important::

   The guard that carries the weight on a domain reaching :math:`R = 0` is a
   different one — see :ref:`running-refusals`. :math:`F/R` is
   :math:`\mu_0 j_\phi`, so a source that does not vanish on the symmetry axis
   is an infinite current density there, and that is a statement about the
   **field** which no definition of :math:`\psiax` can repair.

A run whose flux carries no interior extremum at all is warned about separately:
that is a plasma with no closed surface around an axis, and :math:`\psiax` is
then the edge of nothing.

The ``.nc`` file carries ``axis_normalised_flux``, ``axis_r`` and ``axis_z``
beside ``psi_axis``, so a consumer differencing two runs can see the same thing
without re-deriving it. Their **absence** is informative too: it means no axis
was located.

:cpp:func:`meq::CriticalPointFinder::checkAxis` is the check, and
:cpp:class:`meq::AxisAgreement` is what it returns. The search it runs is seeded
Newton rather than an exhaustive one, so a clean reading is evidence and not
proof; the cost is linear in the mesh and small — 0.04 s over 768 elements and
0.21 s over 12,288, against solves of 1 s and 36 s.

.. _normalised-decoupled:

The control that makes the measurement mean anything
----------------------------------------------------

:cpp:enumerator:`meq::GradShafranovSolver::Normalisation::Decoupled` is the same
solver, mesh, guess, line search and stopping rule with exactly three quantities
zeroed — :math:`c`, :math:`b` and :math:`d - 1`. The step in :math:`\psiax` then
reduces to :math:`\psiax \leftarrow \max\psi_h`, and the trace step is a Newton
step that does not know :math:`\psiax` is about to move. That is
":math:`\psiax` outside the residual, done as favourably as possible".

Measured, the coupled iteration converges and the decoupled one **does not move
at all** — not slowly, not to a different answer: the residual sits where it
started.

That control exists because it is the *only* thing that can see the missing
terms. The finite-difference check on ``dFdPsi`` structurally cannot: ``f`` and
``dFdPsi`` are both evaluated at whatever normalisation is set and agree with
each other however wrong it is. Never run a calculation with ``Decoupled``; it
is there to be measured failing.

What is not done
----------------

:math:`\psibnd` is zero, because MEQ solves the fixed-boundary problem. Free
boundary makes it an unknown as well, which is a **second border row and column
of the same shape** — :cpp:class:`meq::NormalisedSource` is where it would go.
