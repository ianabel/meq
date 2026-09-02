Sources
=======

``#include "meq/Source.hpp"``, ``"meq/RotatingSource.hpp"``,
``"meq/SourceFactory.hpp"``

.. cpp:namespace:: meq

.. cpp:var:: inline constexpr double vacuumPermeability

   :math:`\mu_0` in SI units. Sources take it as a constructor argument so that
   a run in normalised units can set it to 1.

The interface
-------------

.. cpp:class:: Source

   .. cpp:function:: virtual double f( double r, double z, double psi ) const = 0

      :math:`F(r, z, \psi)`: the full right-hand-side numerator, with **no**
      :math:`1/r` applied — that weight belongs to the weak form.

   .. cpp:function:: virtual double dFdPsi( double r, double z, double psi ) const = 0

      The **exact** derivative of :cpp:func:`f`. See :ref:`sources-jacobian` for
      why an error here is invisible to every convergence table.

   Neither may throw: a Newton iterate routinely overshoots, and an exception
   out of a quadrature loop kills a solve that would otherwise have converged.

.. cpp:class:: NormalisedSource : public Source

   A source whose profiles are functions of :math:`\Psi = \psi/\psiax`, so that
   :math:`\psiax` is an unknown of the system. Necessarily of the form
   :math:`F(r, z, \psi) = H(r, z, \psi/\psiax)/\psiax`.

   .. cpp:function:: virtual void setNormalisation( double psiAxis ) = 0

      Called by the solver **once per residual evaluation**, so it must be cheap
      and must not allocate. Throws ``std::invalid_argument`` for a
      :math:`\psiax` that is not finite or is zero.

   .. cpp:function:: virtual double normalisation() const = 0

Implementations
---------------

.. cpp:class:: SolovievSource : public Source

   :math:`F = -\left((1-A)r^2 + A\right)`, with :math:`A + C = 1` baked in.
   :cpp:func:`Source::dFdPsi` is exactly zero, so the problem is linear and
   Newton converges in one step.

   .. cpp:function:: explicit SolovievSource( double a )
   .. cpp:function:: double a() const
   .. cpp:function:: double c() const

.. cpp:class:: MHDSource : public Source

   :math:`F = \mu_0 r^2 p'(\psi) + (gg')(\psi)`.

   .. cpp:function:: MHDSource( std::shared_ptr<Profile const> pPrime, std::shared_ptr<Profile const> ggPrime, double mu0 = vacuumPermeability )

      The profiles are the **derivative** quantities. No :math:`\mu_0` is
      applied to the second term and no sign is applied to either. Throws
      ``std::invalid_argument`` for a null profile or a non-finite
      :math:`\mu_0`.

   .. cpp:function:: Profile const & pPrime() const
   .. cpp:function:: Profile const & ggPrime() const
   .. cpp:function:: double mu0() const

.. cpp:class:: NormalisedMHDSource : public NormalisedSource

   The same thing with profiles in normalised flux.

   .. cpp:function:: NormalisedMHDSource( std::shared_ptr<Profile const> pPrime, std::shared_ptr<Profile const> ggPrime, double psiAxis, double mu0 = vacuumPermeability )

   .. note::

      :cpp:func:`Source::f` divides by :math:`\psiax` **once** and
      :cpp:func:`Source::dFdPsi` divides by it **twice**. Dropping the second is
      the classic error here.

Rotating sources
----------------

.. cpp:var:: inline constexpr std::size_t maxSpecies

   The most species a rotating source will carry. A fixed cap rather than a
   dynamic bound so that the per-quadrature-point work allocates nothing and
   needs no mutable scratch — which a source evaluated from a threaded assembly
   must not have.

.. cpp:struct:: Species

   .. cpp:var:: double mass

      In kilogrammes; positive.

   .. cpp:var:: double charge

      :math:`Z_s`: **signed and dimensionless**, not coulombs. Non-zero.

   .. cpp:var:: std::shared_ptr<Profile const> temperature

      :math:`T_s(\psi)` **in joules**.

   .. cpp:var:: std::shared_ptr<Profile const> density

      :math:`n_{s0}(\psi)` in :math:`\mathrm{m}^{-3}`, on the curve
      :math:`r = r_{\text{ref}}`.

.. cpp:function:: double chargeNeutralityResidual( std::vector<Species> const & species, double psi )

   :math:`\sum_s Z_s n_{s0}(\psi)`.

.. cpp:function:: std::shared_ptr<Profile const> neutralisingDensity( std::vector<Species> const & species, std::size_t index )

   The density the nominated species must have for charge neutrality, **exact at
   every derivative level** rather than differenced. The nominated species' own
   density is never read, so it may be null — which is what lets a caller fill
   the set in one pass.

.. cpp:class:: RotatingSource : public Source

   The generalised Grad–Shafranov source of :cite:t:`Abel2013`. See
   :doc:`../rotation`.

   .. cpp:enum-class:: Closure

      .. cpp:enumerator:: Automatic

         Closed form at two species, root find above. **The default**, and it is
         resolved in the constructor, so :cpp:func:`closure` never returns it.

      .. cpp:enumerator:: ClosedForm

         Two species only; refused for more.

      .. cpp:enumerator:: RootFind

         Any number, and the cross-check at two.

   .. cpp:function:: RotatingSource( std::vector<Species> species, std::shared_ptr<Profile const> omega, std::shared_ptr<Profile const> ggPrime, double referenceRadius, double mu0 = vacuumPermeability, Closure closure = Closure::Automatic )

      ``omega`` **may be null**, meaning no rotation, in which case the source
      reduces to :cpp:class:`MHDSource`'s equation.

      Throws ``std::invalid_argument`` for: fewer than two species; more than
      :cpp:var:`maxSpecies`; a closed form asked for with more than two; a null
      required profile; a non-positive reference radius, mass or temperature; a
      zero charge; charges all of one sign; or charge neutrality violated on the
      reference curve.

   The closure is exposed for testing and for output:

   .. cpp:function:: double potential( double r, double psi ) const

      :math:`e\phi_0` in **joules**, exactly zero at :math:`r = r_{\text{ref}}`
      by construction.

   .. cpp:function:: double density( std::size_t index, double r, double psi ) const
   .. cpp:function:: double dPotentialDPsi( double r, double psi ) const
   .. cpp:function:: double pressure( double r, double psi ) const
   .. cpp:function:: double dPressureDPsi( double r, double psi ) const
   .. cpp:function:: double densityExponent( std::size_t index, double r, double psi ) const

   .. cpp:function:: std::vector<Species> const & species() const
   .. cpp:function:: Closure closure() const
   .. cpp:function:: Profile const * omega() const

      A **null pointer** when there is no rotation.

   .. cpp:function:: Profile const & ggPrime() const
   .. cpp:function:: double referenceRadius() const
   .. cpp:function:: double mu0() const

.. cpp:class:: NormalisedRotatingSource : public NormalisedSource

   A wrapper rather than a reimplementation: it holds a
   :cpp:class:`RotatingSource` by value.

   .. cpp:function:: NormalisedRotatingSource( std::vector<Species> species, std::shared_ptr<Profile const> omega, std::shared_ptr<Profile const> ggPrime, double referenceRadius, double psiAxis, double mu0 = vacuumPermeability, RotatingSource::Closure closure = RotatingSource::Closure::Automatic )

   .. cpp:function:: double potential( double r, double psi ) const
   .. cpp:function:: double density( std::size_t index, double r, double psi ) const
   .. cpp:function:: double pressure( double r, double psi ) const

      These three take **physical** :math:`\psi` and convert internally.

   .. cpp:function:: RotatingSource const & unnormalised() const

      Whose arguments are :math:`\Psi`, not :math:`\psi`.

The factories
-------------

``#include "meq/SourceFactory.hpp"``

.. cpp:function:: std::shared_ptr<Source const> makeSource( SourceConfig const & config, std::string const & configFileName = std::string() )

.. cpp:function:: std::shared_ptr<NormalisedSource> makeNormalisedSource( SourceConfig const & config, std::string const & configFileName = std::string() )

Both throw :cpp:class:`ConfigError` and nothing else — library exceptions from
the sources are caught and translated, so a bad profile file names the key it
came from.

.. important::

   **Each refuses the other's case**, on
   :cpp:func:`SourceConfig::isNormalised`. That refusal is deliberate: a
   normalised source *is-a* ``Source``, so returning one from the ordinary
   factory would compile and solve — with :math:`\psiax` frozen at the guess for
   ever.

   ``makeNormalisedSource`` returns a **non-const** pointer by design, because
   the solver mutates it. The object must outlive the solve.
