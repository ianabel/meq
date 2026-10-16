Toroidal rotation
=================

A tokamak plasma rotating at sonic speed is not described by the ordinary
Grad–Shafranov equation, because **the density is no longer a flux function**.
Centrifugal force sweeps the heavy species to the outboard side, and an
electrostatic potential arises to stop that separating the charges.

MEQ solves the generalised equation of :cite:t:`Abel2013` (their eq. 136),
closed by their expressions for the poloidal density variation and for the
potential, for an arbitrary number of species in a local gauge.

.. code-block:: toml

   [source]
   Type = "rotating"
   ReferenceRadius = 1.0
   Omega = 4.0e5
   GGPrime = 0.8

   [[source.species]]
   Name = "D"
   Mass = 3.3435837768e-27
   Charge = 1.0
   Temperature = 1.0
   TemperatureScale = 1.602176634e-16
   DensityFile = "rotating-density.dat"
   DensityScale = 1.0e20

   [[source.species]]
   Name = "e"
   Mass = 9.1093837015e-31
   Charge = -1.0
   Temperature = 0.8
   TemperatureScale = 1.602176634e-16
   Neutralising = true

.. _rotation-equation:

The equation, in MEQ's convention
---------------------------------

MEQ solves :math:`-\gradbar\cdot\left(\gradbar\psi/R\right) = F/R`, so what
rotation changes is :math:`F`. Written out, with :math:`g \equiv rB_\phi` and
:math:`e` the elementary charge,

.. math::

   F(R, z, \psi) = \mu_0 R^2 \sum_s n_s
       \Big\{ T_s (\ln N_s)'
         + \big[ Z_s e \phi_0 - \tfrac{1}{2} m_s \omega^2 R^2 + T_s \big]
           (\ln T_s)' \Big\}
     + \mu_0 R^4 \omega \omega' \sum_s m_s n_s
     + g g',

closed pointwise by the poloidal density variation and by quasineutrality,

.. math::

   n_s(R, \psi) &= N_s(\psi)\,
       \exp\!\left[\frac{m_s \omega^2(\psi) R^2}{2 T_s(\psi)}
                   - \frac{Z_s e \phi_0}{T_s(\psi)}\right], \\
   0 &= \sum_s Z_s n_s(R, \psi) \qquad \text{which determines } \phi_0(R,\psi).

Here :math:`N_s`, :math:`T_s`, :math:`\omega` and :math:`g` are flux functions;
:math:`n_s` and :math:`\phi_0` are **not**, and that is the whole of what makes
this a different equation.

That is the form as the source paper writes it, in *its* gauge. MEQ fixes the
gauge differently, which replaces :math:`R^2` in the exponent by
:math:`R^2 - r_{\text{ref}}^2` and gives each :math:`N_s` a physical meaning —
see :ref:`the next section <rotation-changes>` for the form MEQ evaluates and
`The gauge`_ for why.

.. _rotation-changes:

.. note::

   **The brace collapses, and that is the form MEQ implements.** Differentiating
   :math:`p = \sum_s n_s T_s` at fixed :math:`R`, the
   :math:`\partial\phi_0/\partial\psi` terms collect into
   :math:`-e\,(\partial\phi_0/\partial\psi)\sum_s Z_s n_s`, which vanishes
   *identically* by quasineutrality. What is left is

   .. math::

      F(R, z, \psi) = \mu_0 R^2
          \left.\frac{\partial p}{\partial\psi}\right|_r + g g',
      \qquad p(R, \psi) = \sum_s n_s(R, \psi)\, T_s(\psi),

   which is :cpp:class:`meq::MHDSource`'s shape with an :math:`R`-dependent
   :math:`p`. MEQ codes :math:`p` and differentiates it, rather than coding the
   brace term by term: there are fewer places to drop a factor, and the brace
   then becomes a *check* on the derivative rather than the thing being checked.

Two independent confirmations, both worth repeating before anything new rests on
this form. It is :cite:t:`Abel2013`'s own force balance projected on
:math:`\nabla\psi`; and at :math:`\omega \to 0` it gives
:math:`F \to \mu_0 R^2 \sum_s p_s' + gg'`, which is that paper's low-Mach
result and is the ordinary Grad–Shafranov source MEQ already solves.

Conventions, and where they differ from the source paper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every one of these has to be right for the equation above to be the equation in
the paper, and none of them is visible in a convergence rate.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * -
     -
   * - **Units**
     - :cite:t:`Abel2013` is Gaussian — :math:`4\pi` where MEQ has
       :math:`\mu_0`, and a :math:`c` MEQ does not carry. Their :math:`I` is
       MEQ's :math:`g`.
   * - **The sign**
     - Their current expression makes :math:`\Delta^\star\psi` *negative* in
       the brace, so :math:`F = -\Delta^\star\psi` is **positive** — matching
       :cpp:class:`meq::MHDSource`'s :math:`F = \mu_0 R^2 p' + gg'`.
   * - :math:`2\pi` **in** :math:`\psi`
     - :math:`\psi` is poloidal flux **per radian**, in Wb/rad. That is what
       MEQ's :math:`gg'` already assumes and what EQDSK tabulates as ``FF'``.
       There is no :math:`2\pi` to insert.
   * - :math:`\Delta^\star`
     - :math:`\partial_{rr} - R^{-1}\partial_r + \partial_{zz}`, on which the
       source paper, :cite:t:`LiZhu2021` and MEQ all agree.

What changes, and what does not
-------------------------------

Exactly one thing changes: the pressure. Each species' density acquires a
poloidal variation,

.. math::

   n_s(R, \psi) = n_{s0}(\psi)\,
       \exp\!\left[\frac{m_s \omega^2 (R^2 - r_{\text{ref}}^2)}{2 T_s}
                   - \frac{Z_s e \phi_0}{T_s}\right],

with :math:`\phi_0` determined by quasineutrality :math:`\sum_s Z_s n_s = 0`.
Everything else about MEQ — the operator, the discretisation, :math:`\tau`, the
hybridization, the estimator, the adaptive loop, the curved boundary — is
untouched.

.. note::

   **The equation collapses to** :math:`F = \mu_0 R^2\,
   \partial p/\partial\psi|_r + gg'` **with** :math:`p = \sum_s n_s T_s`,
   because the :math:`\partial\phi_0/\partial\psi` terms cancel identically
   against quasineutrality. So the residual needs :math:`\phi_0` and never its
   derivative; only the **Jacobian** does.

With :math:`\omega = 0` the source reduces to :cpp:class:`meq::MHDSource`'s
equation exactly, and MEQ checks that it does — both pointwise against the
static source and end to end through the solver, where it reproduces the static
benchmark's errors to every printed digit.

.. _rotation-gauge:

The gauge
---------

:math:`\phi_0` is determined by quasineutrality only up to an additive function
of :math:`\psi`. :cite:t:`Abel2013` resolve this with a flux-surface average,
which MEQ deliberately does not have; MEQ instead pins

.. math::

   \phi_0(r_{\text{ref}}, \psi) = 0

on the curve :math:`R = \texttt{ReferenceRadius}`, a constant radius rather than
the magnetic axis and not a flux-surface quantity at all.

.. important::

   **The consequence is the useful part**: each ``Density`` is then the
   *physical* density of that species on :math:`R = r_{\text{ref}}` — a number a
   user can state and another code can be checked against.

   **The consequence is also the trap**: two sets of densities differing by the
   gauge factor :math:`\exp(Z_s e\,\delta/T_s)` describe *exactly the same
   plasma*. So "our density disagrees with yours" is a statement about
   ``ReferenceRadius`` until proved otherwise, and two codes must agree on it
   before their profiles can be compared at all.

The midplane gauge
~~~~~~~~~~~~~~~~~~

A mirror takes the other local gauge, ``[source] Gauge = "midplane"``:

.. math::

   \phi_0(R_0(\psi), \psi) = 0,

with :math:`R_0(\psi)` the radius at which the surface :math:`\psi` crosses the
midplane :math:`z = 0`. Each ``Density`` is then the physical density **where its
surface crosses the midplane**.

**Why a mirror needs it.** The fixed-radius gauge carries every surface's
densities to one radius through
:math:`\exp\bigl(m_s\omega^2(r_{\text{ref}}^2 - R_0^2)/2T_s\bigr)`, an exponent
*linear in the species mass*. On a mirror at Mach 6 that factor spans hundreds of
decades for carbon, past a double's :math:`e^{709.8}`, so a three-species mirror
at its operating point cannot be written down in the fixed-radius gauge at all.
Under the midplane gauge the factor is one on every surface by construction, and
each ``Density`` table spans only the physical density range.

**Why it is a mirror's only.** It is one condition per surface only where
:math:`\psi` is monotone along :math:`z = 0`. A tokamak's midplane crosses every
surface twice, inboard and outboard of the magnetic axis, so the parse refuses
``Gauge = "midplane"`` without ``SupportShape = "box"``, and the solver checks
the monotonicity itself, on every iterate, rather than assuming it.

**What it costs.** :math:`R_0(\psi)` is a functional of the solution, so the
source is no longer a pointwise function of :math:`(R, z, \psi)`. MEQ tabulates
:math:`R_0(\Psi)` by a root find along the midplane once per Newton iteration of
the bordered solve and *freezes* it within the iteration; its slope
:math:`R_0' = 1/(R_0 q_R)` is read off the solved flux :math:`q`, and the flux
along the midplane is :math:`q` integrated rather than :math:`\psi_h` read
pointwise, since :math:`\psi_h` is discontinuous between elements. The source's
Jacobian is exact for the frozen table, and the table's own dependence on the
field is taken up by the next iteration — so Newton converges **linearly**
rather than quadratically under this gauge. An iterate whose midplane flux is
not monotone keeps the previous table for one more iteration; a solve may not
*converge* on a table its own answer cannot give, and is refused instead.

**The crossings as unknowns, in the library.** The frozen table is what makes
Newton linear, and the library therefore also solves the
midplane gauge with :math:`R_0` at a set of knots in :math:`\Psi` as unknowns of
the bordered Newton, beside :math:`\psi_{\text{ax}}`, :math:`\psi_{\text{bnd}}` and
the exterior modes —
:cpp:func:`meq::GradShafranovSolver::setMidplaneGaugeKnots`. Each knot carries
one equation, that the flux along the midplane at :math:`R_0` be its surface's;
the source's table at every residual is the :math:`C^2` cubic spline of
:math:`R_0^2/2` through the knots, so the source is again a pointwise function
of :math:`(R, z, \psi)` and the knots, and the Jacobian is the whole one: Newton
is quadratic. It costs one more backsolve per knot per step against the same
factorisation. ``meq`` itself solves with the frozen table, and there is no key
for the bordered one: it moves the answer by the interpolant's error between the
knots, which is a choice of discrete problem rather than of speed. A run can
start from another solve's whole state and crossings —
:cpp:func:`meq::GradShafranovSolver::setInitialState` and
:cpp:func:`meq::GradShafranovSolver::setMidplaneCrossingInitial` — so a cheaper
solve can supply the start. Measure the knot count against your own mesh's
error before relying on one.

**And it changes what "flat at a label" means.** The support's precondition is
that :math:`F` vanish at both labels. Under the midplane gauge a ``Density``
flat at a label still carries :math:`-m_s\omega^2 R_0 R_0'/T_s` in its
:math:`\psi`-derivative at fixed :math:`R` wherever :math:`\omega` does not
vanish there, so it is the densities themselves, or the rotation, that must
fall away at the labels — or the tables must be transferred from a
fixed-radius statement of the same plasma, as the test suite's mirror cases are.

``examples/mirror-three-species.toml`` is the worked example: electrons,
deuterons and carbon at Mach 6 on the mirror of :ref:`rotation-mirror`, with
densities that fall to a millionth of their peak at the labels.

.. code-block:: sh

   ../../install/double/bin/meq-run ../../source/examples/mirror-three-species.toml

The vacuum-midplane gauge
~~~~~~~~~~~~~~~~~~~~~~~~~

``[source] Gauge = "vacuum-midplane"`` pins the potential on a **fixed** curve
instead,

.. math::

   \phi_0(R_g(\Psi), \Psi) = 0, \qquad R_g = R_{0,\text{vac}},

where :math:`R_{0,\text{vac}}(\Psi)` is the radius at which the *vacuum's*
surface of normalised flux :math:`\Psi` crosses :math:`z = 0` — the conductors'
own field, normalised by its own flux at the two label points, so that the
curve passes through both. Each ``Density`` is then the density **where the
vacuum's surface crossed the midplane**.

**Why it is a gauge at all.** A gauge pin does not have to lie on the surface
it labels: quasineutrality fixes :math:`\phi_0(R, \psi)` pointwise, at any
:math:`(R, \psi)`, which is also why ``ReferenceRadius`` is a gauge on a mirror
whose surfaces never reach it. Any curve known before the solve is one.

**What it buys.** The curve is fixed, so the source is again a pointwise
function of :math:`(R, z, \psi)` and the Newton Jacobian is exact, border
columns included: Newton is **quadratic**, where the midplane gauge's frozen
crossings make it linear. And it keeps what the midplane gauge is for. The
transfer from a ``Density`` to the density where the *solved* surface crosses
:math:`z = 0` is
:math:`\exp\bigl(m_s\omega^2(R_{0,\text{vac}}^2 - R_0^2)/2T_s\bigr)` with
:math:`\phi_0` aside — an exponent set by how far the plasma moves its own
surfaces, not by the machine's span of radii.

**What it costs.** The tables are stated at the vacuum's crossing, which is not
where the solved surface crosses, so the same tables pose a slightly different
plasma in the two midplane gauges. The ``.nc`` carries every rotating mirror's
densities at its solved crossing (:ref:`the true midplane <output-midplane>`),
and those columns, written as ``DensityFile`` tables, pose the same plasma in
``Gauge = "midplane"``.

**What it needs.** The vacuum is read as the subtracted conductors' own field
:math:`\psi_c`, which is known before the solve only when ``[conductors]``
subtracts them and :math:`\Gamma`'s data is zero, so the gauge is refused
otherwise. With the coils meshed the vacuum would be a solve of its own, and a
curve read off anything else — the starting guess, a previous answer — would
make what each ``Density`` means depend on it. For the same reason there is no
key for a curve of the user's own.

**And the box's precondition is checked before the solve.** The support needs
:math:`F` to vanish at both labels: a jump there, across a contour that cuts
elements, caps the convergence order as a plasma edge does, and the Jacobian
misses its edge term — the lesser cost, which Newton barely notices. On the
vacuum's curve a ``Density`` flat at a label still carries
:math:`-n_s m_s\omega^2 R_g R_g'/T_s` in its :math:`\psi`-derivative wherever
:math:`\omega` does not vanish there; but the curve is known in advance, so
whether the profiles meet the precondition is too. A run whose :math:`F` at a
label exceeds a millionth of its peak, at any radius that label's surface can
reach, is refused before it solves, exit code 1, naming the label. The midplane
gauge cannot be checked the same way: its :math:`R_0'` at a label is the
answer's, and the vacuum's is not close enough to it to refuse on — a profile
set that meets the precondition on its own solved crossings can fail it on the
vacuum's.

.. code-block:: toml

   [source]
   Gauge = "vacuum-midplane"
   # no ReferenceRadius; [conductors] Model = "filament" or "subtracted"

Species
-------

Between two and :cpp:var:`meq::maxSpecies` of them, and the two-species minimum
is not arbitrary: quasineutrality is what determines :math:`\phi_0`, and it
needs charges of both signs.

.. warning::

   **Charge is** :math:`Z_s` — **signed and dimensionless.** Not coulombs.
   :math:`+1` for a proton, :math:`-1` for an electron, :math:`+6` for stripped
   carbon.

   **Temperature is in JOULES.** Not eV, not keV. There is no hidden conversion
   anywhere in MEQ. ``TemperatureScale`` is how a file stays readable — write
   ``1.0`` keV with the conversion in the scale, rather than a raw number with a
   comment.

**Exactly one species must set** ``Neutralising = true``, and it must not carry
a density of its own. Fixing the gauge removes one function's worth of freedom
from the set of densities, so for :math:`n` species there are :math:`n-1`
independent ones. The marked species' density is derived as
:math:`-(1/Z_s)\sum_{\text{others}} Z_{s'} n_{s'0}`, **exact at every derivative
level** rather than differenced. Asking an author for two profiles that happen
to balance is asking for two that do not.

The upper limit is a fixed cap rather than a vector so that the per-quadrature-
point work allocates nothing and needs no mutable scratch — which a source
evaluated from a threaded assembly must not have.

Closing the potential
---------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``Closure``
     - 
   * - ``ClosedForm``
     - **Two species need no root find at all.** After taking logarithms the
       quasineutrality condition is linear in :math:`\phi_0`, giving a single
       exponent that both species share, exact with the electron mass kept.
   * - ``RootFind``
     - Three or more: a safeguarded scalar Newton on the quasineutrality
       condition, with :math:`\phi_0`'s two :math:`\psi`-derivatives obtained by
       implicit differentiation rather than by differencing.
   * - ``Automatic``
     - The default: closed form at two species, root find above.

The two agree to round-off at two species, which is what makes ``RootFind``
usable as a cross-check rather than merely as a fallback.

At two species the condition is linear in :math:`\phi_0` after taking
logarithms, and both species end up sharing one exponent:

.. math::

   e\phi_0 = \frac{\omega^2 (R^2 - r_{\text{ref}}^2)}{2}\,
             \frac{m_1 T_2 - m_2 T_1}{Z_1 T_2 - Z_2 T_1},
   \qquad
   C = \omega^2\,\frac{Z_1 m_2 - Z_2 m_1}{Z_1 T_2 - Z_2 T_1}.

The numerator is a mass-weighted temperature difference, so two species with
equal :math:`m/T` leave nothing for the field to separate and :math:`\phi_0`
vanishes identically. MEQ keeps :math:`m_e/T_e` rather than dropping it: it costs
one term and removes a question.

Above two species the condition is transcendental and is solved at each
evaluation point. **It is as well behaved as such a thing gets**, and that is a
property of the equation rather than of the solver:

.. math::

   \frac{\partial}{\partial\phi_0}\sum_s Z_s n_s
       = -e \sum_s \frac{Z_s^2}{T_s} n_s \;<\; 0
   \qquad\text{strictly,}

because every term carries :math:`Z_s^2`. With at least one charge of each sign
the left-hand side runs monotonically from :math:`+\infty` to :math:`-\infty`,
so the root **exists, is unique, and can be bracketed** — a safeguarded Newton
cannot fail on it. The bracket comes from the two-species formula above, which is
exact at :math:`n = 2` and close when the impurity fraction is small.

.. warning::

   :math:`\phi_0`'s :math:`\psi`-derivatives, which the Jacobian needs, are
   obtained by **implicit differentiation** of the quasineutrality condition and
   never by differencing the root find. Differencing an inner solve from outside
   gives a derivative whose accuracy is the inner tolerance, and Newton then
   degrades from quadratic to linear with no wrong answer and no failing test.

   For the same reason MEQ does **not** cache :math:`\phi_0` across Newton steps,
   although the previous iterate at the same quadrature point would be a good
   guess: anything that makes an evaluation depend on history stops it being a
   function of its arguments, and an assembled Jacobian is then differentiating
   something that is not one.

The cost in derivatives
-----------------------

This is the one structural change rotation forced on the rest of MEQ. A
rotating :math:`F` is *already* :math:`\partial p/\partial\psi` of something
built from flux functions, so the Jacobian spends a **second** derivative of
every input — where :cpp:class:`meq::MHDSource` stores the products :math:`p'`
and :math:`gg'` and needs only one. :cpp:func:`meq::Profile::doublePrime` exists
for this, and no reparametrisation avoids it. See :doc:`profiles`, including the
caveat that a tabulated profile's second derivative jumps at every knot.

.. _rotation-errata:

Three published errors, all of which converge
---------------------------------------------

Recorded because each is exactly the kind that produces a beautiful convergence
table for the wrong equation.

**A Mach number that is really its square.** :cite:t:`LiZhu2021` write
:math:`M_0^2` in their solutions where their prose defines :math:`M_0` as a
group without a square root. The group is an energy ratio, so it *is* a Mach
number squared. MEQ's fixture names its member accordingly and cites the
exponent rather than the symbol.

**Two reversed signs in a derivative.** The same paper's expression for
:math:`\partial p/\partial\psi` carries the wrong sign on both the
:math:`\mathrm{d}\Omega/\mathrm{d}\psi` and the
:math:`\mathrm{d}T/\mathrm{d}\psi` corrections, relative to differentiating its
own definition of :math:`p`. Found independently three times — by
transcription, by an unrelated numerical check, and by MEQ's own derivation
agreeing with the corrected form.

.. important::

   **Neither of that paper's own benchmarks can see it**, because both have
   constant temperature and constant rotation, so the offending terms are
   identically zero. That is the general shape of the hazard here: a benchmark
   that does not exercise a term cannot validate it.

   The same gap is why **no published rotating benchmark exercises the term
   MEQ's Jacobian is most at risk in** — the derivative of the shared density
   exponent with respect to :math:`\psi`. Every exact rotating solution holds
   that exponent constant, because holding it constant is what makes the
   equation solvable in closed form. Only MEQ's own finite-difference sweep
   touches it.

**A citation that is wrong nearly everywhere.** :cite:t:`MaschkePerrin1980` is
a second exact rotating benchmark, and almost every citation of the result names
a different, later paper by the same authors. It is section 4 of the 1980 paper
— temperature a surface quantity — that is the isothermal closure MEQ solves;
that paper's section 3 is a genuine polytrope and is a *different* equation.

That section 4 constrains only the *ratio* of the rotation frequency squared to
the temperature, leaving both free functions of :math:`\psi`. It is therefore
the exact solution MEQ measures :cpp:class:`meq::RotatingSource` against with
**every profile varying** — two temperatures, two densities and the rotation —
whose variations have to cancel to a source independent of :math:`\psi`. Its
equilibrium equation is otherwise the same one :cite:t:`LiZhu2021` solve, so it
is a second test of the source rather than of the discretisation.

.. warning::

   **Rotating-equilibrium papers use at least three different closures that look
   alike on the page**: isothermal on a flux surface, adiabatic/polytropic, and
   variants with the density as a flux function. Check any borrowed closed form
   against the operator by finite differences before trusting it. MEQ's fixtures
   all do.

Normalised flux
---------------

:cpp:class:`meq::NormalisedRotatingSource` puts the profiles in normalised flux,
where :math:`\psiax` becomes an unknown and the bordered Newton of
:doc:`normalised_flux` closes it. ``examples/rotating-normalised.toml`` is the
worked example, and it is the same plasma as ``examples/rotating-rectangle.toml``
re-expressed — which is a good way to see what the two formulations do and do
not have in common.

.. _rotation-mirror:

A magnetic mirror
-----------------

``examples/mirror-rotating.toml`` solves a rotating mirror free boundary from a
file:

.. code-block:: sh

   ../../install/double/bin/meq-run ../../source/examples/mirror-rotating.toml

It is ``geq``'s machine and ``geq``'s problem — two coil packs carried as twenty
point filaments, electrons and deuterons in rigid rotation — which the test
suite also poses through the library and compares with the binary's answer
degree of freedom for degree of freedom.

**A mirror has no magnetic axis.** Its plasma is an annulus between two
prescribed flux surfaces, so both labels are read at points:
:math:`\psi_{\mathrm{ax}}` at ``[source.axis] R, Z``, the inner one, and
:math:`\psi_{\mathrm{bnd}}` at ``[boundary.limiter] R, Z``, the outer one.
Both float with the solve, as ``geq``'s do. The inner label is the *smaller*
flux, so the span :math:`\psi_{\mathrm{ax}} - \psi_{\mathrm{bnd}}` is
**negative**, and nothing needs saying about that: :math:`\Psi` is written in
the span.

**The profiles are in** :math:`\Psi`, **which is** :math:`1 - \psi_N` in
``geq``'s convention: :math:`\Psi = 1` at the inner label and 0 at the outer.
A table carried over from a code that labels the other way round has to have
its abscissa reversed and its slopes negated. Their slopes, and so the source
:math:`F`, must vanish at **both** labels, because :math:`\Psi = 1` is an edge
of the support as well as :math:`\Psi = 0`; the example's slopes fall as
:math:`(1 - \Psi)^2` and :math:`\Psi^3` at the two.

**The support is a box**, ``[source] SupportShape = "box"``:
:math:`0 < \Psi \le 1` and :math:`|z| \le` ``AxialEdge``. The axial end is a
prescribed modelling cutoff at the throat — beyond it the same :math:`\Psi`
recurs in the expander, where the centrifugal potential no longer confines —
so the run reports the current density it cuts through (``axial cutoff`` on the
console, ``axial_cutoff_relative`` in the ``.nc``). ``[mesh.generate]
AlignPlasma`` puts the cutoff on element faces, since a jump across a line that
cuts elements costs the order. :ref:`configuration-mirror-support` has the keys
and what is refused beside them.

``GGPrime = 0``: a mirror has no toroidal field, and that is also what lets the
plasma's flux reach :math:`R = 0` without an infinite current density there
(:ref:`running-refusals`).

Output
------

A rotating run writes a density field per species and the electrostatic
potential beside :math:`\psi` and :math:`\mathbf{B}` in the NetCDF file, each
field's description saying which gauge pinned it. The
potential is written as :math:`e\phi_0` in **joules** rather than volts, because
that is the combination every exponent above contains. A rotating *mirror*, in
any gauge, also writes each species' density where each solved surface crosses
the midplane, as a table that can be fed back as a ``DensityFile``
(:ref:`the true midplane <output-midplane>`).

.. warning::

   Any field derived from the geometry must be evaluated at **the node's own
   radius**. In the band between :math:`\Gamma_h` and :math:`\Gamma` (see
   :ref:`output-band`) it is tempting to reuse the value at the foot on
   :math:`\Gamma_h`, and for a density whose exponent carries :math:`R^2` that
   is wrong by orders of magnitude rather than by the width of the band. This
   was measured, deliberately, as a controlled experiment — the same closed form
   evaluated both ways over the same band nodes — because the trap had just been
   met and dodged for one field while a neighbouring one was left in it.
