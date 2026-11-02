Free boundary
=============

On a fixed boundary :math:`\Gamma` is given and :math:`\psi = 0` on it. On a
**free** boundary the plasma edge is an unknown: the field reaches to infinity
through a vacuum containing the coils, and what is given instead is the coil
currents and the machine. MEQ solves

.. math::

   -\gradbar\cdot\left(\frac{1}{R}\gradbar\psi\right) =
   \begin{cases}
     F(R, z, \psi)/R & \text{in the plasma } \Omega_p(\psi), \\
     \mu_0 j_{\mathrm{coil}} & \text{in the conductors}, \\
     0 & \text{everywhere else, out to infinity,}
   \end{cases}
   \qquad \psi(0, z) = 0, \qquad \psi \to 0 \text{ at infinity},

with :math:`F = [\mu_0 R^2 p'(\Psi) + gg'(\Psi)]` on the plasma, the profiles
functions of :math:`\Psi = (\psi - \psibnd)/(\psiax - \psibnd)`, and the plasma
:math:`\Omega_p(\psi)` the region inside the last closed flux surface — the one
through the limiter contact or the X-point.

Four things are new against a fixed boundary, and they are of very different
sizes: the **unbounded domain**, :math:`\psibnd` as **a second unknown**
beside :math:`\psiax`, the **plasma support** :math:`\Omega_p(\psi)` moving with
the solution, and the **conductors**. This page is how MEQ does each.
:doc:`configuration` has the keys — ``[boundary.exterior]``,
``[boundary.limiter]``, ``[boundary.xpoint]``, ``[source] PlasmaCurrent``,
``[source] ConfineToPlasma`` and ``[[coils]]`` — and :doc:`validation` the
comparison against an independent free-boundary code.

.. code-block:: toml

   [boundary.exterior]
   Radius = 2.4          # rho_Gamma: the semicircle the exterior map lives on
   Modes = 10

   [boundary.limiter]
   R = 1.3375            # psi_bnd is psi_h here, an unknown of the same Newton
   Z = 0.0

   [source]
   Normalised = true
   ConfineToPlasma = true
   PlasmaCurrent = 3.0e5 # amperes: the profile scale becomes an unknown

Coupling at a distance
----------------------

MEQ computes on a polygonal subdomain :math:`\Omega_h` with HDG, puts a smooth
artificial boundary :math:`\Gamma` outside it, represents everything beyond
:math:`\Gamma` by an exact exterior operator, and couples the two **at a
distance** by the same transfer technique that carries a curved boundary
condition (:doc:`curved_boundary`). This is :cite:t:`CockburnSayasSolano2012`,
and the choice is structural rather than a preference: their extension and
lifting are, term for term, the machinery MEQ's curved path already runs, so the
free-boundary coupling is that path with the datum on :math:`\Gamma` an unknown
instead of zero. The coupling is the transmission condition

.. math::

   E_h(q_h)\cdot\nu + \lambda = 0 \quad \text{on } \Gamma,

with :math:`E_h(q_h)` the element extension of the solved flux and
:math:`\lambda` the exterior normal derivative. :math:`\Gamma_h` need not fit
:math:`\Gamma`: optimal order holds with
:math:`\mathrm{dist}(\Gamma_h, \Gamma) = O(h)`, which is the assumption the
adaptive loop maintains.

The classical alternative is the Green's-function route of :cite:t:`Lackner1976`
and von Hagenow, which FreeGS-family codes use. It amortises a precomputed
boundary operator over a fixed grid, and MEQ's adaptive loop changes the domain
every cycle. A boundary-element method with a kernel — :cite:t:`Heumann2015`'s
boundary form — is the right operator assembled the expensive way, as the next
section shows. And an outer fixed point between the interior and the exterior
sits outside the residual, so Newton cannot see it.

.. _free-boundary-semicircle:

The exterior on a semicircle
----------------------------

Take :math:`\Gamma` to be a **semicircle of radius** :math:`\rho_\Gamma`
**centred on the axis**, so the computational domain is a half-disc whose flat
side is :math:`R = 0`. In three dimensions that semicircle is a sphere, and the
exterior operator on a sphere is diagonal.

**The separation.** In spherical coordinates :math:`R = \rho\sin\theta`,
:math:`z = z_c + \rho\cos\theta`, :math:`\mu = \cos\theta`,

.. math::

   \dstar = \partial_{\rho\rho}
          + \frac{1}{\rho^2}\left(\partial_{\theta\theta}
          - \cot\theta\,\partial_\theta\right):

the :math:`-(1/R)\partial_r` that distinguishes :math:`\dstar` from the
axisymmetric Laplacian cancels the Laplacian's :math:`(2/\rho)\partial_\rho`
exactly. Separating :math:`\psi = \rho^\alpha f(\mu)` gives

.. math::

   (1 - \mu^2) f'' + \alpha(\alpha - 1) f = 0,

the **Gegenbauer equation of order** :math:`-1/2`. So :math:`\alpha = n` or
:math:`\alpha = 1 - n`, and the angular functions are

.. math::

   C_n(\mu) = \frac{P_{n-2}(\mu) - P_n(\mu)}{2n - 1}
            = \frac{(1 - \mu^2)\,P'_{n-1}(\mu)}{n(n-1)},
   \qquad n \ge 2,

with :math:`P_n` the Legendre polynomials. Three consequences, and each removes
something:

* :math:`C_n(\pm 1) = 0` **for every** :math:`n \ge 2`, so every mode vanishes on
  the axis identically and the flat side needs no treatment in the exterior.
* **The exterior modes** :math:`\rho^{1-n}` **all decay, and no admissible mode
  is constant.** The degenerate :math:`n = 0, 1` sector — :math:`1`,
  :math:`z/\rho`, :math:`\rho` and :math:`z`, all :math:`\dstar`-harmonic — does
  not vanish on the axis, and the two of those that decay have infinite norm in
  the weight below. So the undetermined constant and the compatibility
  condition of the planar problem do not arise; :math:`\{C_n\}_{n \ge 2}` is
  complete.
* **The Dirichlet-to-Neumann map is diagonal.** If
  :math:`\psi|_\Gamma = \sum_n a_n C_n(\mu)`, then
  :math:`\psi_{\mathrm{ext}} = \sum_n a_n (\rho/\rho_\Gamma)^{1-n} C_n(\mu)`
  and

  .. math::

     \left.\frac{\partial\psi}{\partial\rho}\right|_\Gamma
       = \sum_n a_n\,\frac{1 - n}{\rho_\Gamma}\,C_n(\mu).

**And the weight is MEQ's own.** The Gegenbauer functions of order
:math:`-1/2` are orthogonal in :math:`(1 - \mu^2)^{-1}`, and on a semicircle
centred on the axis

.. math::

   \frac{\mathrm{d}\Gamma}{R} = \frac{\rho_\Gamma\,\mathrm{d}\theta}
                                     {\rho_\Gamma \sin\theta}
                              = \frac{\mathrm{d}\mu}{1 - \mu^2},

the weight the Grad–Shafranov weak form carries anyway, with :math:`\rho_\Gamma`
cancelling out of it. So

.. math::

   \int_\Gamma C_m C_n \frac{\mathrm{d}\Gamma}{R} = \delta_{mn} h_n,
   \qquad h_n = \frac{2}{n(n-1)(2n-1)},

and the exterior adds to the coupled system a **diagonal** block with entry
:math:`(n-1) h_n/\rho_\Gamma = 2/\big(n(2n-1)\rho_\Gamma\big)` — no layer
potentials, no singular quadrature, no elliptic integrals.

.. note::

   **The second form of** :math:`C_n` **is the one to compute with.** The
   printed difference of Legendre polynomials cancels exactly at the axis, where
   both values tend to 1, and loses digits without bound as :math:`\mu \to \pm
   1`. The form carrying :math:`(1 - \mu^2)` as a factor is exactly zero on the
   axis, makes :math:`C_n/(1 - \mu^2)` a polynomial so a plain Gauss–Legendre
   rule integrates the mass exactly, and pairs with
   :math:`\mathrm{d}C_n/\mathrm{d}\mu = -P_{n-1}(\mu)`, which needs no
   difference at all. The printed form is also valid for :math:`n \ge 2` only:
   below it it returns :math:`1 - \mu` for both :math:`n = 0` and :math:`n =
   1`, the wrong function rather than a division by zero.

:cpp:class:`meq::ExteriorDtN` is that statement — the basis, the symbol
:math:`(1 - n)/\rho_\Gamma`, the mass :math:`h_n` and the field at an exterior
point — and it is MFEM-free, so continuous integration tests it. Its modes are
**indexed by degree from 2**, so a loop over modes is a loop over degrees, and an
error there is an off-by-two. It is checked against two references that share
the equation and nothing else: the field of a circular current loop by complete
elliptic integrals, expanded on :math:`\Gamma` and pushed through the symbol;
and :cite:t:`Heumann2015`'s boundary form, a double integral with a
hypersingular elliptic-integral kernel, assembled against the Gegenbauer basis —
which comes out diagonal, with the diagonal equal to the block entry above.

.. important::

   **The diagonality is a property of a semicircle centred on the axis and of no
   other curve.** That is why ``[boundary.exterior]`` requires
   ``[mesh] RMin = 0``: on a domain starting at :math:`R = 0.05` the modes do not
   span the exterior at all.

.. _free-boundary-coupled-system:

The coupled system
------------------

The unknowns are the whole HDG state :math:`x = (q, \psi, \hat\psi)`, the
Gegenbauer coefficients :math:`a = (a_2, \ldots, a_{N+1})`, and the scalars of
the normalisation and its constraints:

.. math::

   \begin{aligned}
   \mathcal{F}(x, a, \psiax, \psibnd) &= 0
     && \text{the HDG residual, with the transferred datum on } \Gamma_h \\
   T_m(x, a) = \int_\Gamma E_h(q_h)\cdot\nu\, C_m\,\mathrm{d}\Gamma
       - a_m\,\frac{(1 - m)\,h_m}{\rho_\Gamma} &= 0
     && m = 2 \ldots N+1 \\
   G_{\mathrm{ax}} = \psiax - \psi_h(x^\ast) &= 0
     && \text{at the located axis, see } \textit{Normalised flux} \\
   G_{\mathrm{bnd}} = \psibnd - (\text{limiter or X-point functional}) &= 0
   \end{aligned}

with one more row for a prescribed current and two more for a located X-point,
below.

**Where** :math:`a` **enters.** The datum on :math:`\Gamma` is
:math:`\sum_n a_n C_n`, carried back to :math:`\Gamma_h` along the transfer
paths and imposed as a load on the flux equation, :math:`\langle \hat\psi,
v\cdot n\rangle`. Each mode is evaluated at the **foot** of its path, on
:math:`\Gamma`, and not at the face of :math:`\Gamma_h` it is imposed on; the two
differ by :math:`O(\mathrm{dist}(\Gamma_h, \Gamma)) = O(h)`, which is exactly
the error the transfer exists to remove. :math:`\hat\psi` enters every HDG row
linearly and the non-linearity is :math:`F(R, z, \psi)` alone, so
:math:`\partial\mathcal{F}/\partial a` **does not depend on the iterate**: the
columns are a function of the geometry and are built once per mesh.

**The transmission row** integrates the extended flux over the true
:math:`\Gamma`, swept by MFEM's ``ExtensionBoundaryQuadrature``, in the plain
measure :math:`\mathrm{d}\Gamma`: :math:`q = \gradbar\psi/R` already carries the
:math:`1/R` the orthogonality wants, and :math:`\mathrm{d}\Gamma/R` would divide
by it twice. On a semicircle the normal is :math:`\hat\rho` analytically.

**The bordered solve.** Newton on the whole system is the field's Jacobian with
a dense border, eliminated against **one** factorisation per step:

.. math::

   \begin{bmatrix} J & C \\ B & D \end{bmatrix}
   \begin{bmatrix} \delta x \\ \delta s \end{bmatrix}
   = -\begin{bmatrix} \mathcal{F} \\ G \end{bmatrix},
   \qquad
   \delta x = -y - Z\,\delta s,\quad Jy = \mathcal{F},\quad JZ = C,

with :math:`s` every scalar unknown and :math:`D` the small dense corner. The
cost is one backsolve per border column with a non-zero field part — the
:math:`N` exterior columns, :math:`\psiax`'s and :math:`\psibnd`'s, and the
profile scale's — applied together through MFEM's batched entry point, and a
dense solve of the corner. It is :doc:`normalised_flux`'s bordered Newton with
the scalar corner replaced by a matrix. Under the default ordering, where
:math:`\psi` and :math:`q` are unknowns, the entries are exact rather than
differenced:

* the exterior columns are :math:`-P`, the projections above, exactly;
* :math:`D`'s exterior block is diagonal, :math:`(n-1)h_n/\rho_\Gamma`;
* the :math:`\psiax` and :math:`\psibnd` rows are the potential shape functions
  of one element at one point (see below), and their diagonal entries are 1;
* the :math:`\psiax` and :math:`\psibnd` columns are
  :math:`\partial F/\partial\psiax` and :math:`\partial F/\partial\psibnd` in
  closed form (:ref:`normalised-analytic-column`).

**What the border does cost is one re-assembly per accepted step**, because the
datum reaches the system as a load. With the exterior coupling alone the
residual is affine in :math:`(x, a)` and an exact Jacobian finishes it in one
step, which is how the coupling is tested.

.. note::

   **The coupled Jacobian can be inverted but not applied.** MFEM's
   ``DarcyNPCOperator`` gradient is solve-only — after its local factorisations
   its ``Mult()`` aborts — so any scheme wanting a product with the full coupled
   operator, a Jacobian-free Krylov method over the border included, has that
   constraint to design around.

.. _free-boundary-axis:

The domain reaches the axis
---------------------------

The half-disc's flat side is :math:`R = 0`, where the flux mass form
:math:`(Rq, v)` degenerates. It is survivable, and what it costs is measured on
harmonic fixtures that vanish on the axis:

* :math:`\psi` **keeps** :math:`k+1` on a mesh reaching the axis, at every
  degree;
* :math:`q` **loses about half an order** at low degree, because the element
  touching the axis carries a weight of order :math:`h` and an unweighted norm
  counts it in full;
* **the trace system's conditioning grows as** :math:`O(1/h)`, **not**
  :math:`O(1/h^2)`: that same element has the smallest diagonal in the system.
  A direct trace solver is nearly insensitive to it at MEQ's sizes, and it would
  matter to an iterative one.

The corner where :math:`\Gamma` meets the axis is benign: it is a corner between
the arc and a fitted straight side, and the lifting's weight vanishes there, so
:math:`\psi` converges at :math:`k+1` across it.

.. important::

   **A domain reaching the axis needs** :math:`F(0, z) = 0`. :math:`F/R` is
   :math:`\mu_0 j_\phi`, and :math:`F = \mu_0 R^2 p' + gg'` leaves
   :math:`gg'` on the axis. Since :math:`\psi(0, z) = 0` exactly, the axis sits
   at :math:`\Psi = -\psibnd/(\psiax - \psibnd)`, which a limiter makes
   *negative* — the vacuum, where an unconfined profile extrapolates and hands
   back a current density growing like :math:`1/R`. The discrete space does not
   forbid a non-zero value at :math:`R = 0` the way the continuous one does, so
   the load is unbounded there and :math:`\psi_h` grows a layer along the whole
   axis. ``[source] ConfineToPlasma`` is the statement that the vacuum carries no
   current, and the driver refuses the case otherwise
   (:ref:`running-refusals`).

.. _free-boundary-borders:

The plasma boundary flux, and the other borders
-----------------------------------------------

:math:`\psibnd` enters only through :math:`\Psi`, so it is a functional of the
solution exactly as :math:`\psiax` is, and it is closed the same way: one more
row and column of the border. Which row depends on what bounds the plasma.

**A prescribed limiter contact** — ``[boundary.limiter] R, Z`` — makes
:math:`\psibnd = \psi_h(R, Z)` inside the element containing the point. The row
is that element's potential shape functions there, exact. The point does not
move with the solution, so there is no position term to argue away.

.. warning::

   **Taking the nearest potential degree of freedom instead is an**
   :math:`O(h)` **error at every degree**, not an :math:`O(h^{k+1})` one: the
   nearest node is up to half a node spacing from the contact, so the error is
   that distance times :math:`|\nabla\psi|`, and raising the degree does not
   shrink it. It is kept as ``LimiterConstraint::NearestDof``, the control, and
   it is the instructive kind of wrong: a sweep of the requested contact is a
   staircase of bit-identical solves, and a fixture whose contact happens to sit
   beside a node converges cleanly to the wrong number. A clean convergence
   table is not evidence that a boundary condition is right.

**A limiter surface** — ``[boundary.limiter] SurfaceAttribute`` — makes
:math:`\psibnd = \max\psi_h` over a polygon of mesh faces, with the contact
**found**. The row is again the shape functions at the contact, for both kinds
of contact a polygon admits: on an edge the tangential derivative vanishes at a
maximum while the contact moves tangentially, and at a vertex it does not move,
so the envelope theorem removes the position term. Newton's observed order is
what checks that, since no error norm can.

**An X-point** — ``[boundary.xpoint]`` — is a functional of the solution that
moves as Newton moves, so its position becomes two more unknowns:

.. math::

   q_r(R_X, z_X) = 0, \qquad q_z(R_X, z_X) = 0, \qquad
   \psibnd - \psi_h(R_X, z_X) = 0.

:math:`\psi` is analytic at an X-point — :math:`\psi - \psi_X \approx
\tfrac12(a\xi^2 - b\eta^2)` — and neither the equation nor the flux mass
degenerates there; only the level curve through it has a corner. The corner
block :math:`\partial(q_r, q_z)/\partial(R_X, z_X)` is exact, because
:math:`q_h` is a polynomial on its element and Newton differentiates the
*discrete* residual. The two new unknowns cost no backsolve: the field residual
reaches them only through :math:`\psibnd`, so their columns are zero. What is
not smooth is the point crossing into another element, which moves the residual
by :math:`O(h^{k+1})`; the element is frozen within a Jacobian and re-decided
at each accepted step, as for the axis.

**A prescribed plasma current** — ``[source] PlasmaCurrent`` — makes the
profile scale :math:`\lambda` an unknown and adds the row
:math:`\int F/R - \mu_0 I_p = 0`. Every entry is analytic: the column is the
source term over :math:`\lambda`, since :math:`F` is linear in it; the row is
:math:`\int (\partial F/\partial\psi)/R\,\phi_j` on the potential block; and the
corner entries against :math:`\psiax` and :math:`\psibnd` are **not zero**,
because :math:`\int F/R` depends on both through :math:`\Psi`. Leave them out
and the solve still converges, linearly — the answer is right and the rate is
gone.

.. important::

   **A free scale is what makes a free boundary posable.** With the amplitude
   fixed and the support moving, the problem is a non-linear eigenvalue problem:
   with :math:`p' \sim \Psi^j`, :math:`\Lambda = A/(\psiax - \psibnd)^2` must be
   an eigenvalue *on the plasma region*, which is itself unknown, so scaling the
   amplitude changes nothing. A fixed amplitude closes only when it came from an
   equilibrium of these very inputs — consistent inputs, not a guess. The
   current border is how production codes pose it, and how MEQ's machine cases
   do.

.. _free-boundary-support:

The plasma support
------------------

With ``ConfineToPlasma`` the source is zero outside the plasma, so the region
carrying current moves with the iterate. Three things decide whether that is
solvable and at what order.

**The support is connected, not pointwise.** :math:`\Psi > 0` is not the plasma:
a level of :math:`\psi` can cut the domain into several pieces, and beyond an
X-point the private flux region carries :math:`\Psi > 0` too. MEQ takes the
connected component containing the magnetic axis, by a flood fill on the
element adjacency graph — ``meq::PlasmaComponent``, MFEM-free — refreshed
before every residual and every Jacobian so both see the same support. The fill
traverses only elements strictly inside the level set; the band of elements
straddling the separatrix, every one of which carries :math:`\Psi > 0` at some
vertex, is shared out between components by a **watershed**, one breadth-first
wave from each at once. An inclusive candidate rule would join the lobes through
that band all along the divertor legs, not merely at the saddle, and blocking
the saddle's own element is resolution-dependent. The watershed needs no
X-point and no parameter. And it does not make the residual jump: a lobe never
changes hands whole, only a straddling element where :math:`\Psi \approx 0` and
the profile with it, so the step between neighbouring states falls like
:math:`h^2` and the fill stays inside Newton. The seed is the element holding
the located axis, excluding conductors and the symmetry axis.
``GradShafranovSolver::setPlasmaConnectivity()`` keeps the pointwise test as
``PlasmaConnectivity::Pointwise``, the control.

**The profile's edge order caps the order and decides solvability.** With
:math:`p' \sim \Psi^j` at the edge the exact :math:`\psi` carries
:math:`|d|^{j+2}` across it, so the best approximation by any method caps
:math:`\psi^\star` at :math:`\min(k + 2, j + 2.5)`: :math:`k+2` holds exactly
when :math:`k \le j`. That cap is the profile's, and it is the same for an exact
cut quadrature rule as for MEQ's plain one, which is why MEQ adopts no cut rule.
With fixed quadrature points the assembled Jacobian is the exact derivative of
the assembled residual. At :math:`j = 0` that residual is **discontinuous** in
the unknowns, a quadrature point crossing the edge making :math:`F` jump, and the
moving support does not converge at any degree, on any mesh, from the exact
solution or under ``PicardThenNewton``. **So** :math:`j \ge 1` **is a
precondition of the free-boundary path**, not a convention. ``freegs4e``'s
:math:`(1 - \Psi_n^{\alpha})^{\beta}` has :math:`j = \beta`, and a published
profile with :math:`\beta < 2` also has an unbounded :math:`\partial
F/\partial\psi` at the edge.

**An X-point on the edge costs nothing.** :math:`\Psi` vanishes quadratically
at the null, so :math:`F \sim \Psi^j` vanishes to order :math:`2j` there, and
the corner the support acquires is measured to cost the best approximation no
order.

**On a diverted plasma the support is frozen within a solve.** The fill's
choice of elements is not differentiable, and on a diverted machine letting it
move inside Newton stalls the solve. ``[solver] PlasmaSupportSweeps`` fixes it
within a solve and re-decides it between solves, stopping when a sweep changes
nothing; see that key in :doc:`configuration`. A near-double-null makes the
choice of bounding saddle arbitrarily sensitive, and that is ill-conditioning
in the problem rather than a defect.

.. _free-boundary-conductors:

Conductors
----------

A conductor carrying a uniform current density contributes
:math:`F = \mu_0 R\,j_\phi` on its cross-section — :math:`\mu_0 R I/|\Omega_c|`,
proportional to :math:`R` because it is the density that is uniform — in the
same units as the plasma term, so the two simply add. Its current is data and
contributes nothing to the Jacobian. :doc:`configuration`'s ``[[coils]]`` and
``[conductors]`` have the models MEQ offers for a conductor inside the domain.

**Inside the domain, mesh to it.** :math:`F` jumps at a uniform conductor's
edge, and where that edge cuts an element the quadrature integrates a
discontinuous integrand with a rule assuming smoothness. The signature is a rate
that *falls* with :math:`k` — a regularity limit is flat in :math:`k`, a
quadrature error worsens as the rest improves — and aligning the mesh to the
conductor restores :math:`k+1` up to the cap a rectangle's corners impose. A
conductor is prescribed input and can always be meshed to; the plasma edge moves
and cannot.

**Outside** :math:`\Gamma`, **through the coupling instead.** On the interior
route a conductor the mesh does not reach is never sampled: the run converges
and describes a machine with it switched off. The exterior problem is linear, so
split

.. math::

   \psi = \psi_{\mathrm{coil}} + \tilde\psi,

with :math:`\psi_{\mathrm{coil}}` the conductor's own field in closed form. The
conductor is outside :math:`\Gamma \supseteq \partial\Omega`, so
:math:`\dstar\psi_{\mathrm{coil}} = 0` inside :math:`\Omega` and the interior
equation is untouched; and :math:`\tilde\psi` is :math:`\dstar`-harmonic in the
whole exterior and decays, so the Gegenbauer expansion is valid for it. The
conductor therefore enters only through :math:`\Gamma`, as known terms on both
halves of the transmission condition:

.. math::

   \psi|_\Gamma = \sum_n a_n C_n + \psi_{\mathrm{coil}}|_\Gamma, \qquad
   q\cdot\nu|_\Gamma = \frac{1}{R}\sum_n a_n\frac{1-n}{\rho_\Gamma} C_n
                       + q_{\mathrm{coil}}\cdot\nu|_\Gamma.

No new unknowns and no change to the border. Both halves land together or not
at all: the datum without the transmission term is an inconsistent pair that
still converges. ``GradShafranovSolver::setExteriorConductors()`` takes a
``meq::ExteriorCoilSet`` — rectangles and point filaments, with the field and
its gradient in closed form and **no source term**, since an exterior conductor
contributes nothing to the interior equation; handing the solver interior coils
there is a compile error. A conductor reaching inside :math:`\Gamma` is refused,
measured from the nearest point of each member rather than its centre. It is a
library capability; the driver carries every ``[[coils]]`` conductor inside the
domain.

.. note::

   **The two routes answer different questions.** The coefficients :math:`a_n`
   *are* the plasma's field outside :math:`\Gamma`, so the plasma's share of
   :math:`B` at an exterior conductor is a modal evaluation there. But the mode
   count is chosen against :math:`\psiax`, not against :math:`B` near
   :math:`\Gamma`; a multipole expansion is worst just outside :math:`\Gamma`
   and improves outward; and differentiating it costs a mode. A filament has no
   self-field at all. So mesh a conductor whose field you need, and take the
   exterior route for one far enough out that meshing to it is wasted elements.

.. _free-boundary-modes:

Choosing the radius and the mode count
--------------------------------------

The coefficients fall geometrically, at a rate set by the ratio of the plasma's
and the conductors' extent to :math:`\rho_\Gamma`, so the mode count is small and
:math:`\rho_\Gamma` trades modes against elements: a larger semicircle needs
fewer modes and more mesh. The mesh between the plasma and :math:`\Gamma` is the
source-free, :math:`\dstar`-harmonic vacuum, the smoothest part of the problem,
so it is where a coarse mesh is cheapest and an elongated plasma does not force
a large mesh.

**Read the spectrum in energy, not in raw coefficients.** :math:`a_n` carries
the mode's own normalisation, :math:`h_n \sim n^{-3}`, so raw coefficients of a
converged trace need not decay at all.
``ExteriorDtN::modeAmplitudes()`` returns :math:`|a_n|\sqrt{h_n}`, and
``ExteriorDtN::truncationRatio()`` — the larger of the **last two** of those
over the largest — is what the driver prints on every coupled run, advising more
modes when it is large. The last two because :math:`C_n(-\mu) = (-1)^n
C_n(\mu)`: an up-down symmetric trace has identically zero odd modes, and the
last mode alone reads exactly zero on half of all mode counts. The advice
threshold is calibrated against the move in :math:`\psiax`, about a per cent
there; sweep ``Modes`` on your own case before trusting a value.

The truncation is also what the adaptive loop cannot see from the interior error
alone; :ref:`adaptivity-boundary` is the boundary indicator that refines
:math:`\Gamma_h` for it.

Which equilibrium is reached
----------------------------

**Once the boundary is free, the starting state chooses which equilibrium is
reported.** The same file reaches different discrete equilibria from different
guesses, each satisfying every border to round-off. The initial
:math:`\psibnd` is such a lever, with no setting right for every machine:
started near its answer on a limited machine it puts the edge on the limiter
from the first step and stalls, while started at zero it lands cleanly — and on
a diverted machine the two starts reach two equilibria that refinement does not
merge. So it is not a configuration key. A spurious equilibrium has tells: a
profile scale far from :math:`O(1)` under ``PlasmaCurrent``, a located axis
inside a conductor, or a negative scale. The driver warns on the second and
checks the axis (:ref:`normalised-axis-check`); :ref:`normalised-guess` is the
general point, and :doc:`validation` shows what consistent inputs reach.
