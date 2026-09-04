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
   MEQ's Jacobian is most at risk in**. Only MEQ's own finite-difference sweep
   over profiles with genuine :math:`\psi`-dependence in temperature and
   rotation touches it.

**A citation that is wrong nearly everywhere.** :cite:t:`MaschkePerrin1980` is
a second exact rotating benchmark, and almost every citation of the result names
a different, later paper by the same authors. It is section 4 of the 1980 paper
— temperature a surface quantity — that is the isothermal closure MEQ solves;
that paper's section 3 is a genuine polytrope and is a *different* equation.

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
