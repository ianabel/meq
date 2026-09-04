Source terms
============

Everything MEQ knows about the physics arrives through one small interface:

.. code-block:: cpp

   class Source
   {
   public:
       virtual double f( double r, double z, double psi ) const = 0;
       virtual double dFdPsi( double r, double z, double psi ) const = 0;
   };

:cpp:func:`meq::Source::f` returns :math:`F(r, z, \psi)` — the full numerator of
the right-hand side, with **no** :math:`1/r` **applied**, because that weight
belongs to the weak form. Coordinates are cylindrical, in metres, with
:math:`r > 0`.

.. _sources-jacobian:

Every source must differentiate itself
--------------------------------------

.. important::

   :cpp:func:`meq::Source::dFdPsi` **must be the exact derivative of**
   :cpp:func:`meq::Source::f`. **A source that cannot differentiate itself
   cannot be used.**

This is the price of using Newton where the source papers use Picard, and it is
not negotiable. :cite:t:`SanchezVizuetSolano2019` state the design MEQ is
reversing: they keep :math:`F` as opaque problem data so that the solver
"relies only on the discretization of the toroidal operator :math:`\dstar`", and
pay for it by iterating on *every* source — even one linear in :math:`\psi` that
could have been folded into the bilinear form. Newton makes the opposite trade,
putting :math:`\partial F/\partial\psi` into the operator as a mass term on the
potential block. :doc:`nonlinear` argues the case.

.. warning::

   **An error in** ``dFdPsi`` **does not change the answer.** Newton converges
   to the same discrete solution whatever Jacobian carried it there, so every
   error norm and every convergence rate is unchanged. What changes is the path:
   the observed order of the Newton iteration drops from 2 to 1, and it takes
   more steps. Past a certain size of error it diverges instead.

   This is why ``dFdPsi`` is checked against a finite difference of ``f`` in the
   unit tests, why the *assembled* Jacobian is checked against a finite
   difference of the assembled residual, and why the driver prints an observed
   order. See :ref:`testing-ladder`.

.. note::

   **Differentiate the discrete residual, not the continuous equation.** MEQ
   gets this for free — the library differentiates the assembled operator — but
   it is worth knowing that it is a deliberate property. A continuous-level
   Newton derivative for the plasma-current term exists, and
   :cite:t:`Heumann2015` explicitly decline to use it, on the grounds that
   there is no evidence it holds for equilibria whose boundary contains an
   X-point and that one of its terms appears to blow up where :math:`\psi`
   reaches a critical point. The shortcut fails exactly where the physics is
   interesting.

The sources that ship
---------------------

:cpp:class:`meq::SolovievSource` — ``Type = "soloviev"``
   :math:`F = -\left((1-A)r^2 + A\right)`, with the flux normalised so that
   :math:`A + C = 1`. It does not depend on :math:`\psi` at all, so
   :math:`\partial F/\partial\psi \equiv 0`, the problem is **linear**, and
   Newton converges in one step. That is exactly why it is the first case to
   run and the fixture that isolates the discretisation from the solver. The
   sign is the subject of :ref:`formulation-soloviev-sign`.

:cpp:class:`meq::MHDSource` — ``Type = "mhd"``
   :math:`F = \mu_0 r^2 p'(\psi) + (g g')(\psi)`, built from two tabulated
   profiles. Note what the profiles are: the **derivative** quantities
   :math:`\mathrm{d}p/\mathrm{d}\psi` and :math:`g\,\mathrm{d}g/\mathrm{d}\psi`
   — the latter is what EQDSK calls ``FF'``. No :math:`\mu_0` is applied to the
   second term and no sign is applied to either.

   :math:`\mu_0` is a constructor argument (``Mu0``) so that a run in
   normalised units can set it to 1.

:cpp:class:`meq::RotatingSource` — ``Type = "rotating"``
   A rotating multi-species plasma. See :doc:`rotation`.

A *manufactured* source is also reachable as ``Type = "manufactured"``: the
nonlinear test case of :cite:t:`SanchezVizuetSolano2019`, whose
:math:`\psi`-dependence is linear, quadratic *and* exponential, and which is
what exercises the Newton path end to end.

Normalised sources are a different interface
--------------------------------------------

.. warning::

   **A source whose profiles are functions of normalised flux must be a**
   :cpp:class:`meq::NormalisedSource`, **not a** :cpp:class:`meq::Source`.

When profiles are given against :math:`\Psi = \psi/\psiax`, the axis flux
:math:`\psiax` is a *functional of the solution*. The Jacobian then acquires
non-local terms through it, which ``dFdPsi`` structurally cannot carry — and
which **the finite-difference check on** ``dFdPsi`` **cannot see either**,
because ``f`` and ``dFdPsi`` are both evaluated at whatever normalisation
happens to be set and agree with each other however wrong it is.

:doc:`normalised_flux` is the account of what MEQ does instead. The practical
consequence here is that the two construction routes are separate and each
refuses the other's case:

.. code-block:: cpp

   auto source = meq::makeSource( config.getSource() );            // throws if normalised
   auto normalised = meq::makeNormalisedSource( config.getSource() ); // throws if not

That refusal is deliberate rather than defensive. A normalised source *is-a*
``Source``, so returning one from the ordinary factory would compile and solve
— with :math:`\psiax` frozen at the initial guess for ever, converging
beautifully to a different equilibrium.

.. _sources-trivial-branch:

The trivial branch
------------------

.. warning::

   **Every source in the standard benchmark family vanishes at**
   :math:`\psi = 0`, **so with homogeneous Dirichlet data** :math:`\psi \equiv 0`
   **solves the problem** — and Newton, which starts from the Dirichlet data,
   lands on it and stops in *zero* iterations with an identically zero residual.

This is not a hypothesis about the benchmarks; the published algorithm for them
literally opens by demanding a non-trivial initial guess. ``[initialguess] Type
= "ramp"`` is what MEQ offers: a linear ramp in :math:`z` from
:math:`-\texttt{Amplitude}` to :math:`+\texttt{Amplitude}`, so that
:math:`\psi = 0` falls in the *interior* rather than being the whole answer.
The amplitude must be positive, and the error message says why.

A source that is nowhere zero on the domain — which a physical configuration
with a real pressure gradient generally is — cannot land on the trivial branch
and needs no guess.

.. note::

   The trivial branch is a *different* phenomenon from the multiple solutions
   discussed in :ref:`nonlinear-multiple`, and confusing the two wastes time.
   The trivial branch is an exact zero solution admitted by the boundary
   conditions; the other is two genuinely distinct non-zero solutions of an
   under-resolved discrete system.

Writing your own
----------------

Derive from :cpp:class:`meq::Source`, implement both methods, and hand the
object to :cpp:func:`meq::GradShafranovSolver::setSource`. The source is
**borrowed**: it must outlive the solver.

.. code-block:: cpp

   class MySource : public meq::Source
   {
   public:
       double f( double r, double /*z*/, double psi ) const override
       {
           return alpha*r*r*psi + beta;
       }
       double dFdPsi( double r, double /*z*/, double /*psi*/ ) const override
       {
           return alpha*r*r;
       }
   private:
       double alpha, beta;
   };

Three rules the shipped sources all follow, and which are easy to break:

* **Do not throw from** ``f`` **or** ``dFdPsi``. A Newton iterate routinely
  overshoots into unphysical values, and an exception out of a quadrature loop
  kills a solve that would otherwise have converged. Clamp instead — which is
  what :doc:`profiles` do.
* **Be cheap and allocate nothing.** These are called at every quadrature point
  of every element on every residual evaluation, and potentially from a threaded
  assembly.
* **Keep MFEM out of it.** ``Source`` and ``Profile`` take plain ``double``
  arguments deliberately, so that both are unit-testable without the finite
  element library — and they are tested that way, in the half of the build that
  continuous integration can actually run.
