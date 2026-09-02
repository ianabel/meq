Boundary shapes
===============

``#include "meq/BoundaryShape.hpp"``

.. cpp:namespace:: meq

.. cpp:class:: ShapeError : public std::runtime_error

   Thrown when a shape cannot be used: degenerate parameters, a curve that
   crosses the axis, or one that is not star-shaped about its own centre.

.. cpp:class:: BoundaryShape

   A closed plasma boundary in the Miller extended harmonic parametrisation
   :cite:p:`ArbonCandyBelli2021`. See :doc:`../curved_boundary`.

   .. cpp:function:: BoundaryShape( double r0In, double z0In, double minorIn, double elongationIn, std::vector<double> cosIn = {}, std::vector<double> sinIn = {} )

      ``cosIn`` is :math:`c_0, c_1, \ldots` and **starts at** :math:`c_0`, the
      tilt; ``sinIn`` is :math:`s_1, s_2, \ldots` and starts at :math:`s_1` —
      there is no :math:`s_0`. Both may be empty.

   .. cpp:function:: static BoundaryShape miller( double r0In, double z0In, double minorIn, double elongationIn, double deltaIn, double squarenessIn = 0.0 )

      The Miller D-shape :cite:p:`Miller1998`, mapped onto the harmonic form
      with :math:`s_1 = \arcsin\delta` — **arcsin, not** :math:`\delta`.
      Throws :cpp:class:`ShapeError` unless :math:`|\delta| < 1`.

   .. cpp:function:: void point( double theta, double & r, double & z ) const

      ``theta`` is the parametrisation's poloidal angle, **not** the polar angle
      about the centre.

   .. cpp:function:: double levelSet( double r, double z ) const

      Negative inside, zero on the curve, positive outside.

      .. note::

         This is the **radial gap** — the distance from the centre minus the
         distance from the centre to the curve at the same polar angle — and
         **not a signed distance**, which is larger. The curve radius at a given
         polar angle is found by bisecting on ``theta``, which is legitimate
         exactly when the polar angle is strictly increasing in it.

   .. cpp:function:: double polarAngle( double theta ) const
   .. cpp:function:: void boundingBox( double & rMin, double & rMax, double & zMin, double & zMax ) const
   .. cpp:function:: double majorRadius() const
   .. cpp:function:: double centreHeight() const
   .. cpp:function:: double minorRadius() const
   .. cpp:function:: double elongation() const

   The constructor throws :cpp:class:`ShapeError` for a non-positive minor
   radius, elongation or major radius; for a curve whose **sampled bounding
   box** reaches the axis, since the operator's :math:`1/r` is not integrable
   there — the check is on the box rather than on :math:`R_0 - a`, because the
   harmonics move the innermost point; and for a curve that is not star-shaped,
   naming the angle at which monotonicity failed.
