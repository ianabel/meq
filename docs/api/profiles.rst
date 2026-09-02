Profiles
========

``#include "meq/Profiles.hpp"``

.. cpp:namespace:: meq

.. cpp:class:: Profile

   A one-dimensional function of flux. See :doc:`../profiles` for the
   conventions, which are load-bearing.

   .. cpp:function:: virtual double operator()( double psi ) const = 0
   .. cpp:function:: virtual double prime( double psi ) const = 0

      The **exact** derivative of :cpp:func:`operator()`, not a finite
      difference of it.

   .. cpp:function:: virtual double doublePrime( double psi ) const = 0

      The exact derivative of :cpp:func:`prime`. Pure virtual so that every
      implementation had to answer rather than one silently returning zero; it
      exists because a rotating source's Jacobian spends a second derivative of
      every input.

   None of the three may throw for a :math:`\psi` outside the profile's natural
   range. The convention is a **constant extension**: the value clamps and both
   derivatives are zero.

.. cpp:class:: ConstantProfile : public Profile

   .. cpp:function:: explicit ConstantProfile( double value )
   .. cpp:function:: double value() const

.. cpp:struct:: Knot

   .. cpp:var:: double psi
   .. cpp:var:: double value
   .. cpp:var:: double derivative

   The slope is **data**, not inferred — which is what lets
   :cpp:func:`Profile::prime` be an analytic derivative rather than an
   approximation.

.. cpp:class:: SplineProfile : public Profile

   A piecewise Hermite cubic through tabulated knots.

   .. cpp:type:: RealFunction = std::function<double( double )>

   .. cpp:function:: explicit SplineProfile( std::vector<Knot> data )

      Throws ``std::invalid_argument`` for fewer than two knots, non-finite
      data, or knots not strictly increasing in :math:`\psi`.

   .. cpp:function:: SplineProfile( RealFunction f, RealFunction fPrime, unsigned int intervals )

      Samples the two functions on :math:`[0, 1]`.

   .. cpp:function:: static SplineProfile fromStream( std::istream & is )
   .. cpp:function:: static SplineProfile fromFile( std::string const & fileName )

      See :ref:`profiles-file-format`. Throws ``std::runtime_error`` on an
      unreadable or malformed file; never returns an empty profile.

   .. cpp:function:: void write( std::ostream & os ) const

      Emits the same format at full precision, so writing and reading back is
      the identity.

   .. cpp:function:: std::vector<Knot> const & knots() const
   .. cpp:function:: std::pair<double, double> domain() const
   .. cpp:function:: std::size_t numIntervals() const
   .. cpp:function:: HermiteCubicSpline intervalAt( std::size_t i ) const

.. cpp:class:: ScaledProfile : public Profile

   .. cpp:function:: ScaledProfile( std::shared_ptr<Profile const> inner, double scale )

      Multiplies at all three derivative levels. A negative scale is allowed —
      a sign convention is a scale like any other. Throws
      ``std::invalid_argument`` for a null inner profile or a non-finite scale.

   .. cpp:function:: Profile const & unscaled() const
   .. cpp:function:: double scale() const

.. cpp:class:: HermiteCubicSpline

   A single interval, not itself a :cpp:class:`Profile`.

   .. cpp:function:: HermiteCubicSpline( double lower, double upper, double fLower, double fUpper, double fPrimeLower, double fPrimeUpper )
   .. cpp:function:: HermiteCubicSpline( Knot const & lower, Knot const & upper )

      Both throw ``std::invalid_argument`` unless the interval has positive
      width.

   .. cpp:function:: double operator()( double x ) const
   .. cpp:function:: double prime( double x ) const
   .. cpp:function:: double doublePrime( double x ) const
   .. cpp:function:: std::pair<double, double> interval() const
   .. cpp:function:: std::pair<double, double> values() const
   .. cpp:function:: std::pair<double, double> derivatives() const

   .. note::

      Out of range, the value is clamped and **both derivatives are zero** —
      the exact derivative of the constant extension. The code this replaced
      returned a function *value* from ``prime()``, which is a units error.

.. cpp:function:: std::ostream & operator<<( std::ostream & os, SplineProfile const & profile )
