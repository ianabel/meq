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

   :math:`\psiax` **is the largest NODAL value, and that is a definition rather
   than an approximation.** It differs from the maximum of the polynomial by
   :math:`O(h^{k+1})` and both converge to :math:`\max\psi`. The nodal one is
   chosen because it is what makes the constraint differentiable in a form the
   border can use.

Under the default nonlinear ordering — see :ref:`nonlinear-ordering` — two of
the three bordered quantities are not finite differences at all: with
:math:`\psi` an unknown of the system, :math:`b` is exactly :math:`-e_j` (one
entry) and :math:`d` is exactly 1.

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
