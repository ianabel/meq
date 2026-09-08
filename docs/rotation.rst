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
   DensityFile = "examples/rotating-density.dat"
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

MEQ solves :math:`-\gradbar\cdot\left(\gradbar\psi/r\right) = F/r`, so what
rotation changes is :math:`F`. Written out, with :math:`g \equiv rB_\phi` and
:math:`e` the elementary charge,

.. math::

   F(r, z, \psi) = \mu_0 r^2 \sum_s n_s
       \Big\{ T_s (\ln N_s)'
         + \big[ Z_s e \phi_0 - \tfrac{1}{2} m_s \omega^2 r^2 + T_s \big]
           (\ln T_s)' \Big\}
     + \mu_0 r^4 \omega \omega' \sum_s m_s n_s
     + g g',

closed pointwise by the poloidal density variation and by quasineutrality,

.. math::

   n_s(r, \psi) &= N_s(\psi)\,
       \exp\!\left[\frac{m_s \omega^2(\psi) r^2}{2 T_s(\psi)}
                   - \frac{Z_s e \phi_0}{T_s(\psi)}\right], \\
   0 &= \sum_s Z_s n_s(r, \psi) \qquad \text{which determines } \phi_0(r,\psi).

Here :math:`N_s`, :math:`T_s`, :math:`\omega` and :math:`g` are flux functions;
:math:`n_s` and :math:`\phi_0` are **not**, and that is the whole of what makes
this a different equation.

That is the form as the source paper writes it, in *its* gauge. MEQ fixes the
gauge differently, which replaces :math:`r^2` in the exponent by
:math:`r^2 - r_{\text{ref}}^2` and gives each :math:`N_s` a physical meaning —
see :ref:`the next section <rotation-changes>` for the form MEQ evaluates and
`The gauge`_ for why.

.. _rotation-changes:

.. note::

   **The brace collapses, and that is the form MEQ implements.** Differentiating
   :math:`p = \sum_s n_s T_s` at fixed :math:`r`, the
   :math:`\partial\phi_0/\partial\psi` terms collect into
   :math:`-e\,(\partial\phi_0/\partial\psi)\sum_s Z_s n_s`, which vanishes
   *identically* by quasineutrality. What is left is

   .. math::

      F(r, z, \psi) = \mu_0 r^2
          \left.\frac{\partial p}{\partial\psi}\right|_r + g g',
      \qquad p(r, \psi) = \sum_s n_s(r, \psi)\, T_s(\psi),

   which is :cpp:class:`meq::MHDSource`'s shape with an :math:`r`-dependent
   :math:`p`. MEQ codes :math:`p` and differentiates it, rather than coding the
   brace term by term: there are fewer places to drop a factor, and the brace
   then becomes a *check* on the derivative rather than the thing being checked.

Two independent confirmations, both worth repeating before anything new rests on
this form. It is :cite:t:`Abel2013`'s own force balance projected on
:math:`\nabla\psi`; and at :math:`\omega \to 0` it gives
:math:`F \to \mu_0 r^2 \sum_s p_s' + gg'`, which is that paper's low-Mach
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
       :cpp:class:`meq::MHDSource`'s :math:`F = \mu_0 r^2 p' + gg'`.
   * - :math:`2\pi` **in** :math:`\psi`
     - :math:`\psi` is poloidal flux **per radian**, in Wb/rad. That is what
       MEQ's :math:`gg'` already assumes and what EQDSK tabulates as ``FF'``.
       There is no :math:`2\pi` to insert.
   * - :math:`\Delta^\star`
     - :math:`\partial_{rr} - r^{-1}\partial_r + \partial_{zz}`, on which the
       source paper, :cite:t:`LiZhu2021` and MEQ all agree.

What changes, and what does not
-------------------------------

Exactly one thing changes: the pressure. Each species' density acquires a
poloidal variation,

.. math::

   n_s(r, \psi) = n_{s0}(\psi)\,
       \exp\!\left[\frac{m_s \omega^2 (r^2 - r_{\text{ref}}^2)}{2 T_s}
                   - \frac{Z_s e \phi_0}{T_s}\right],

with :math:`\phi_0` determined by quasineutrality :math:`\sum_s Z_s n_s = 0`.
Everything else about MEQ — the operator, the discretisation, :math:`\tau`, the
hybridization, the estimator, the adaptive loop, the curved boundary — is
untouched.

.. note::

   **The equation collapses to** :math:`F = \mu_0 r^2\,
   \partial p/\partial\psi|_r + gg'` **with** :math:`p = \sum_s n_s T_s`,
   because the :math:`\partial\phi_0/\partial\psi` terms cancel identically
   against quasineutrality. So the residual needs :math:`\phi_0` and never its
   derivative; only the **Jacobian** does.

With :math:`\omega = 0` the source reduces to :cpp:class:`meq::MHDSource`'s
equation exactly, and MEQ checks that it does — both pointwise against the
static source and end to end through the solver, where it reproduces the static
benchmark's errors to every printed digit.

The gauge
---------

:math:`\phi_0` is determined by quasineutrality only up to an additive function
of :math:`\psi`. :cite:t:`Abel2013` resolve this with a flux-surface average,
which MEQ deliberately does not have; MEQ instead pins

.. math::

   \phi_0(r_{\text{ref}}, \psi) = 0

on the curve :math:`r = \texttt{ReferenceRadius}`, a constant radius rather than
the magnetic axis and not a flux-surface quantity at all.

.. important::

   **The consequence is the useful part**: each ``Density`` is then the
   *physical* density of that species on :math:`r = r_{\text{ref}}` — a number a
   user can state and another code can be checked against.

   **The consequence is also the trap**: two sets of densities differing by the
   gauge factor :math:`\exp(Z_s e\,\delta/T_s)` describe *exactly the same
   plasma*. So "our density disagrees with yours" is a statement about
   ``ReferenceRadius`` until proved otherwise, and two codes must agree on it
   before their profiles can be compared at all.

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

   e\phi_0 = \frac{\omega^2 (r^2 - r_{\text{ref}}^2)}{2}\,
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

Output
------

A rotating run writes a density field per species and the electrostatic
potential beside :math:`\psi` and :math:`\mathbf{B}` in the NetCDF file. The
potential is written as :math:`e\phi_0` in **joules** rather than volts, because
that is the combination every exponent above contains.

.. warning::

   Any field derived from the geometry must be evaluated at **the node's own
   radius**. In the band between :math:`\Gamma_h` and :math:`\Gamma` (see
   :ref:`output-band`) it is tempting to reuse the value at the foot on
   :math:`\Gamma_h`, and for a density whose exponent carries :math:`r^2` that
   is wrong by orders of magnitude rather than by the width of the band. This
   was measured, deliberately, as a controlled experiment — the same closed form
   evaluated both ways over the same band nodes — because the trap had just been
   met and dodged for one field while a neighbouring one was left in it.
