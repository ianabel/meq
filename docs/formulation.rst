The discretisation
==================

This page describes what the solver actually does: the equation it forms, the
spaces it forms it in, and the conventions those spaces carry. The parts you
cannot skip if you are extending meq are :ref:`the two sign conventions
<sign-conventions>`.

The equation
------------

meq solves the fixed-boundary Grad–Shafranov problem
:cite:p:`SanchezVizuetSolano2019`,

.. math::

   -\gradbar \cdot \left( \frac{1}{r} \gradbar \psi \right)
       = \frac{F(r, z, \psi)}{r} \quad \text{in } \Omega \subset \mathbb{R}^2,
   \qquad \psi = 0 \ \text{ on } \Gamma = \partial\Omega,

with

.. math::

   F(r, z, \psi) := \mu_0 r^2 \frac{\mathrm{d}p}{\mathrm{d}\psi}
                  + g \frac{\mathrm{d}g}{\mathrm{d}\psi}.

.. important::

   :math:`\gradbar := (\partial_r, \partial_z)` **acts formally like a vector
   of partial derivatives and is not the cylindrical gradient.** The
   distinction is the whole content of the :math:`1/r` and :math:`r` weights
   above. Writing the operator as a cylindrical divergence of a cylindrical
   gradient gives a different equation.

   The operator on the left is :math:`-\dstar`, the Grad–Shafranov operator.

:math:`p(\psi)` is the plasma pressure and :math:`g(\psi)/r` the toroidal field
function; both are user input, and it is their :math:`\psi`-dependence that
makes the problem semi-linear. Everything meq knows about the physics arrives
through :math:`F` and :math:`\partial F/\partial\psi` — see :doc:`sources`.

The first-order system
----------------------

Introduce the **flux**

.. math::

   q := \frac{1}{r}\gradbar\psi,

and the problem becomes a first-order system

.. math::

   q - \frac{1}{r}\gradbar\psi = 0, \qquad
   -\gradbar\cdot q = \frac{F}{r}, \qquad
   \psi = 0 \ \text{ on } \Gamma.

Introducing :math:`q` is not a numerical convenience. The physically
interesting output is the magnetic field, which is built from
:math:`\gradbar\psi`; discretising :math:`q` directly gives the derivative at
the *same* order as the potential rather than one order down. That is the whole
reason to prefer a mixed method for this equation.

The LDG-H method
----------------

On a triangulation :math:`\Th` of the domain, the hybridizable discontinuous
Galerkin method in LDG-H form :cite:p:`SanchezVizuetSolano2019,SanchezVizuet2020adaptive`
seeks :math:`(q_h, \psi_h, \hat\psi_h)` satisfying

.. math::

   (r\, q_h, v)_{\Th} + (\psi_h, \gradbar\cdot v)_{\Th}
       - \langle \hat\psi_h, v\cdot n \rangle_{\partial\Th} &= 0 \\
   (q_h, \gradbar w)_{\Th}
       - \langle \hat{q}_h\cdot n, w \rangle_{\partial\Th}
       &= \left( \frac{F}{r}, w \right)_{\Th} \\
   \langle \hat{q}_h\cdot n, \mu \rangle_{\partial\Th \setminus \Gamma_h} &= 0 \\
   \hat\psi_h &= \varphi_h \ \text{ on } \Gamma_h

for all test functions, with the **numerical flux**

.. math::

   \hat{q}_h\cdot n := q_h\cdot n + \tau\,(\psi_h - \hat\psi_h)
   \qquad \text{on } \partial\Th .

The third equation is the transmission condition: it says the normal numerical
flux is single valued across every interior face, and it is what couples the
elements to one another.

.. _formulation-spaces:

The spaces
----------

.. list-table::
   :header-rows: 1
   :widths: 14 20 26 40

   * - Space
     - Holds
     - MFEM
     - Notes
   * - :math:`V_h = [P_k(K)]^2`
     - the flux :math:`q_h`
     - ``L2_FECollection``, ``vdim 2``
     - Discontinuous, element by element.
   * - :math:`W_h = P_k(K)`
     - the potential :math:`\psi_h`
     - ``L2_FECollection``
     - Discontinuous.
   * - :math:`M_h = P_k(e)`
     - the trace :math:`\hat\psi_h`
     - ``DG_Interface_FECollection``
     - Lives on faces only. **This is the globally coupled unknown.**

All three are of the **same degree** :math:`k`. That is not an oversight or a
simplification: hybridization removes the inf-sup compatibility condition that
a classical mixed method would impose between the flux and potential spaces
:cite:p:`CockburnGopalakrishnanLazarov2009`, so equal order is admissible and
gives :math:`k+1` convergence in both :math:`\psi` and :math:`q`.

.. note::

   **The volume spaces use a closed nodal (Gauss–Lobatto) basis, and that is a
   convention rather than a requirement.** A nodal basis does not change the
   *space*: Lobatto and Legendre both span :math:`P_k(K)`, so the
   discretisation cannot see the choice, and nothing in the hybridization needs
   volume degrees of freedom to sit on faces — every coupling is a face
   *integral*, computed by quadrature against both bases. This was measured
   both ways and the answers agree to the last place.

   It is kept for alignment with the MFEM miniapp meq was ported from. **What
   it costs** is that a degree of freedom is a point value *on* the element
   boundary, where a discontinuous field is ambiguous — so reading
   :math:`W_h` by nodal interpolation at another mesh's node points is not
   well defined. :cpp:class:`meq::FieldTransfer` projects rather than
   interpolating, which is basis-agnostic and is what a non-nested transfer
   needs in any case.

Hybridization
-------------

The point of the method is that the flux and potential can be eliminated
**element by element**. Every element's contribution to the global system is
expressed in terms of the trace on its own faces, leaving one globally coupled
system in :math:`\hat\psi_h` alone, whose size is the number of face degrees of
freedom — independent of how many volume unknowns each element carries.

That elimination is a small dense solve per element, and it is where meq's
Newton iteration meets its most interesting structural question, because when
:math:`F` depends on :math:`\psi` the elimination can itself be nonlinear. That
is the subject of :ref:`nonlinear-ordering`, and it is the single largest
performance and robustness lever in the code.

.. _formulation-tau:

The stabilisation :math:`\tau`
------------------------------

:math:`\tau` is a constant, and its default is **1**. Both source papers set it
there and note that optimal order requires only :math:`\tau = O(1)`. It is
configurable as ``[discretisation] Tau``.

.. warning::

   **MFEM will not give you a constant** :math:`\tau` **unless you ask.** The
   HDG diffusion integrator's built-in stabilisation is the LDG choice, scaled
   by the inverse local mesh size *and* by the diffusion coefficient — which is
   not what the papers use, and which was measured to cost **a full order in
   the flux**: :math:`q` converges at :math:`k` rather than :math:`k+1`, while
   :math:`\psi` still converges at :math:`k+1`.

   **A convergence study of** :math:`\psi` **alone would have passed it.** That
   is the argument for measuring both norms in every rate table, and it is why
   :cpp:class:`meq::ConstantStabilization` exists.

.. important::

   **Keep** :math:`\tau` **constant unless something measured says otherwise.**
   A solution-dependent stabilisation must also supply its own derivative for
   the Newton Jacobian, and omitting that gives "no wrong answer, only slow
   Newton convergence" — a failure that survives a passing regression suite. A
   constant :math:`\tau` never reaches that code path at all.

   Deriving :math:`\tau` from the local coefficient is specifically not
   advised.

.. _sign-conventions:

The two sign conventions
------------------------

.. warning::

   **The assembled flux block holds** :math:`-q`, **not** :math:`q`.

MFEM's ``DarcyForm`` is built for :math:`u = -k\nabla p`, the opposite sign to
:math:`q = \gradbar\psi/r`, and the integrators that make the hybridization
consistent have that sign baked in and take no scaling argument — so there is
no way to flip it during assembly.

:cpp:func:`meq::GradShafranovSolver::flux` undoes it once, into a separate grid
function. The block vector inside the solver stays in the library's convention,
because a Newton residual assembled by the library expects it there.

The consequences, all of which have bitten:

* :cpp:func:`meq::GradShafranovSolver::flux` and
  :cpp:func:`meq::GradShafranovSolver::postProcessedFlux` return meq's
  :math:`+q`. :cpp:func:`meq::GradShafranovSolver::totalFlux` is in the
  library's convention, deliberately, because the constraint equation it is
  projected through is written that way.
* :cpp:func:`meq::poloidalField` and
  :cpp:func:`meq::GridSampler::samplePotentialWithFlux` want ``flux()``.
* :cpp:func:`meq::GradShafranovSolver::transferredDatum` internally wants the
  **raw** block. Feeding it ``flux()`` gives :math:`-\psi` where :math:`\psi`
  was wanted — the answer with its sign reversed, not a small bias.
* The potential right-hand side is assembled as :math:`-(F/r, w)` to match.

The second convention concerns :math:`\tau` itself:

.. warning::

   :math:`\tau` **carries the opposite sign to the published numerical flux.**

The papers print :math:`\hat q\cdot n := q\cdot n + \tau(\psi - \hat\psi)`.
Testing the system against itself with that sign gives

.. math::

   (r q, q) - \tau \|\psi - \hat\psi\|^2_{\partial\Th} = 0,

which is indefinite, so the element-local solves are not guaranteed invertible.
The stable sign is :math:`-\tau`, which is exactly what assembling in the
library's convention with a positive :math:`\tau` produces.

This is a well-posedness argument rather than an observed failure — both signs
converge at :math:`k+1` on the benchmark. It is the second sign slip in the
same pair of papers as the one below.

.. _formulation-soloviev-sign:

The papers disagree about the sign of the Solov'ev source
---------------------------------------------------------

Checked, resolved, and recorded here so that nobody rediscovers it.
:cite:t:`SanchezVizuetSolano2019` gives :math:`F = -((1-A)r^2 + A)`;
:cite:t:`SanchezVizuet2020adaptive` gives :math:`F = +((1-A)r^2 + A)`.

**The first is right**, and the second contradicts its own statement of the
equation. Applying :math:`\dstar` to the particular solution *both* papers
publish settles it analytically:

.. math::

   \dstar\!\left(\tfrac{1}{8} r^4\right) &= r^2 \\
   \dstar\!\left(\tfrac{A}{2} r^2 \ln r\right) &= A \\
   \dstar\!\left(-\tfrac{A}{8} r^4\right) &= -A r^2

so that, summing,

.. math::

   \dstar \psi_P = (1-A) r^2 + A

and since both papers define :math:`-\dstar\psi = F`, the source is
:math:`F = -((1-A)r^2 + A)`. The twelve homogeneous terms contribute nothing,
being :math:`\dstar`-harmonic.

meq does not take this on trust. The Solov'ev fixture recomputes
:math:`\dstar\psi` by central differences and the test suite asserts it against
the very :math:`F` the solver is fed, so the whole twelve-term transcription is
checked rather than believed.

.. warning::

   Take this as the standing caution about analytic benchmarks: **the published
   coefficients are not self-checking**, and a sign error here shows up as a
   solver that converges beautifully to the wrong equilibrium. See
   :ref:`testing-fixtures`.

Boundary conditions
-------------------

On the **fitted** path, :math:`\Gamma` is the mesh boundary and the trace
degrees of freedom there are essential: they carry the Dirichlet datum
directly, and ``[boundary] Type = "zero"`` sets it to zero.

On the **extension** path, :math:`\Gamma_h` is the boundary of a polygonal
subdomain strictly inside the plasma, and the datum :math:`\varphi_h` imposed
there is *transferred* from the true :math:`\Gamma` along paths. Two things come
off a face marked that way: the HDG stabilisation (:math:`\tau` is zero on
:math:`\Gamma_h`, and leaving it on loses one order at :math:`k = 1` and two at
:math:`k = 2`, measured) and the flux constraint. See :doc:`curved_boundary`.

.. note::

   **On the extension path** :cpp:func:`meq::GradShafranovSolver::trace` **is
   not** :math:`\hat\psi` **on** :math:`\Gamma_h`. Those trace degrees of
   freedom remain in the essential list at zero; the datum actually imposed is
   :math:`\varphi_h`, and it is never stored anywhere.
   :cpp:func:`meq::GradShafranovSolver::transferredDatum` rebuilds it on
   request, which is what the error estimator needs.

.. _formulation-symmetry:

Is the trace matrix symmetric?
------------------------------

On a fitted mesh, **yes**, to round-off, and negative definite besides — so
:math:`-A` is symmetric positive definite and every symmetric Krylov method and
a Cholesky factorisation would apply. The Newton Jacobian is the same, for
:math:`\partial F/\partial\psi` of either sign.

On the **extension path, no**, and that is not a defect — it is what the
transfer technique *is*.

The extension integrator's element matrix is an outer product of the normal
trace of one basis function against the *path lifting* of another: two
unrelated vectors, so the contribution is structurally unsymmetric, and it is
substantially larger than the mass term it sits beside. It goes straight into
the hybridization's per-element flux block, so it reaches both the local
factorisation and, through the Schur complement, the global trace matrix.

.. important::

   **The transfer technique costs self-adjointness.** The continuous
   Grad–Shafranov operator is self-adjoint; :math:`\dstar` with a *transferred*
   Dirichlet datum is not, because the datum on a face depends on the flux
   along a path leaving the element. That is a property of the method
   :cite:p:`CockburnSolano2012`, not of this implementation — and it means
   meq's headline configuration is genuinely non-symmetric, which is why an
   unsymmetric sparse LU is the right solver for it. See
   :doc:`linear_solvers`.

What survives on both paths is that the **symmetric part is negative
definite**, which is what gives a Krylov method a convergence bound and makes
preconditioning meaningful.

.. warning::

   **Measure the free–free block, not the assembled matrix.** Imposing
   essential boundary conditions zeroes rows without eliminating the matching
   columns, which makes even the fitted-path Jacobian look wildly unsymmetric
   when the whole matrix is examined. That is a boundary-condition artefact and
   says nothing about the operator.

   More generally: this is the place where measuring only the easy
   configuration was most misleading. Symmetry held to round-off on a fitted
   rectangle and failed outright on the geometry meq is actually for.
