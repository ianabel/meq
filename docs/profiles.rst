Profiles
========

A **profile** is a one-dimensional function of flux that a source is built from
— a pressure gradient, a temperature, a density, a rotation frequency.
:cpp:class:`meq::Profile` is the interface:

.. code-block:: cpp

   class Profile
   {
   public:
       virtual double operator()( double psi ) const = 0;
       virtual double prime( double psi ) const = 0;
       virtual double doublePrime( double psi ) const = 0;
   };

Three derivative levels, and all three are pure virtual
-------------------------------------------------------

``prime()`` must be the **exact** derivative of ``operator()``, not a finite
difference of it, and ``doublePrime()`` must be the exact derivative of
``prime()``. The reason for the first is :ref:`sources-jacobian`: a Newton
Jacobian assembled from an approximate derivative converges linearly and still
reaches the right answer, which makes the defect invisible to every error norm.

The reason for the **second** is toroidal rotation, and it is the one structural
change that feature forced. :cpp:class:`meq::MHDSource` stores the *products*
:math:`p'` and :math:`gg'`, so :math:`F` is one profile evaluation and
:math:`\partial F/\partial\psi` is one ``prime()`` — two levels, which is all
the interface used to have. A rotating :math:`F` is *already*
:math:`\partial p/\partial\psi` of something built from flux functions, so its
Jacobian spends a second derivative of every input. No reparametrisation avoids
it. See :doc:`rotation`.

``doublePrime()`` is pure virtual rather than defaulting to zero so that **every**
implementation had to answer, rather than one silently returning the wrong thing.

.. warning::

   **A tabulated profile is** :math:`C^1` **and no more, so its second
   derivative jumps at every interior knot.** MEQ's spline is a piecewise
   Hermite cubic, which matches value and slope across a knot but not curvature.
   The test suite asserts that jump rather than pretending otherwise. Nothing
   has measured what it costs a Newton iteration; if you are chasing an
   iteration count on a tabulated profile, this is a candidate.

Out of range, nothing throws
----------------------------

.. important::

   **Implementations must not throw when handed a** :math:`\psi` **outside their
   natural range.** A Newton iterate routinely overshoots, and an exception
   thrown out of a quadrature loop kills a solve that would otherwise have
   converged.

The convention is a **constant extension**: below the first knot the value is
the first knot's value, above the last it is the last knot's, and the first and
second derivatives are **zero** outside the table — which is the exact
derivative of that constant extension.

.. warning::

   **Zero derivative outside the table means the source term switches itself
   off there, silently.** For a pressure profile that is a plasma with no
   pressure gradient outside the tabulated range. Both shipped example tables
   extend well past the flux range their run actually visits, for exactly this
   reason, and a table given in normalised flux needs to extend past
   :math:`[0, 1]` because early Newton iterates put :math:`\Psi` anywhere at
   all.

The implementations
-------------------

:cpp:class:`meq::ConstantProfile`
   A constant. Both derivatives are exactly zero, so a source built from these
   contributes nothing to the Newton Jacobian — which is the Solov'ev case, and
   why it converges in one step.

:cpp:class:`meq::SplineProfile`
   A piecewise Hermite cubic through tabulated ``(psi, value, slope)`` knots.
   The slope is **data**, not inferred, which is what lets ``prime()`` be the
   analytic derivative of ``operator()`` rather than an approximation of it —
   and that is what keeps the Jacobian consistent with the residual. A
   values-only table would have to invent slopes.

:cpp:class:`meq::ScaledProfile`
   Multiplies another profile at all three derivative levels. This is what
   ``*Scale`` keys in a configuration file become, and it is how a table stays
   readable — a temperature written as ``1.0`` keV with the conversion in a
   scale factor, rather than as a raw number in joules with a comment.

   A negative scale is allowed: a sign convention is a scale like any other.

.. _profiles-file-format:

The tabulated file format
-------------------------

One knot per line, whitespace separated:

.. code-block:: text

   # psi [Wb/rad]     n_D0 [1e20 m^-3]     dn_D0/dpsi [1e20 m^-3 per Wb/rad]
     -0.02             0.1                  20.0
      0.00             0.5                  20.0
      0.05             1.5                  20.0
      0.10             2.5                  20.0
      0.15             3.5                  20.0

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Rule
     - 
   * - Three columns
     - :math:`\psi`, :math:`f(\psi)`, :math:`f'(\psi)`. All three are required.
   * - Comments
     - A line whose first non-blank character is ``#``.
   * - **A blank line ends the table**
     - So several profiles can be concatenated in one stream.
   * - At least two knots
     - Fewer is an error.
   * - Strictly increasing in :math:`\psi`
     - And every value must be finite.

``write()`` emits exactly this format, terminated by a blank line, at full
precision — so writing and reading back is the identity.

Relative paths in a configuration file are resolved against the **working
directory of the run**, not against the directory the configuration file is in.

.. warning::

   **Re-expressing a profile against normalised flux divides the abscissa by the
   axis flux and MULTIPLIES the derivative column by it.** Since
   :math:`\mathrm{d}/\mathrm{d}\Psi = \psiax\,\mathrm{d}/\mathrm{d}\psi`, the
   values are unchanged and the slopes are not.

   Hand a table to the wrong configuration and it parses, solves and converges —
   to a plasma whose gradient is wrong by a factor of :math:`\psiax`. MEQ ships
   the same profile written both ways, in ``examples/rotating-density.dat`` and
   ``examples/rotating-density-normalised.dat``, and both headers say which is
   which. There is no way for the code to detect the mistake; the units are in
   the comment or nowhere.

A profile that is linear in flux
--------------------------------

Worth calling out, because it is a common way to write a test case that quietly
tests nothing. If a density or pressure profile is **linear** in :math:`\psi`,
then :math:`\mathrm{d}p/\mathrm{d}\psi` does not depend on :math:`\psi`, so
:math:`F` does not either, so :math:`\partial F/\partial\psi` vanishes
identically and Newton finishes in one step. The problem is affine.

That is a perfectly good thing to want — it isolates whatever else is being
tested from the nonlinear solve — but it means the run says nothing whatever
about the Jacobian. Two of MEQ's own shipped examples are like this deliberately,
and their headers say so; ``examples/mhd-rectangle.toml`` is the counterexample,
with a pressure gradient quadratic in :math:`\psi`. See :doc:`examples`.
