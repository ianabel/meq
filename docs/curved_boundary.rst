Curved boundaries
=================

The plasma boundary :math:`\Gamma` is a smooth curve. A mesh is not. meq does
not attempt to reconcile the two by fitting the mesh: it solves on a polygonal
subdomain lying strictly inside the plasma and **transfers** the boundary
condition outward, which is the technique of :cite:t:`CockburnSolano2012`
applied to Grad–Shafranov by :cite:t:`SanchezVizuet2020adaptive`.

.. code-block:: toml

   [boundary.shape]
   Type = "miller"
   R0 = 1.5
   MinorRadius = 0.5
   Elongation = 1.6
   Triangularity = 0.35

Why not just mesh the shape
---------------------------

Because a polygonal approximation to a curved boundary **caps the convergence
rate whatever the polynomial degree**, and the cap is severe.

At an interior angle :math:`\pi - 2\pi/N` of an :math:`N`-gon, the solution
acquires a corner singularity whose exponent approaches 1 as :math:`N` grows,
and interior pollution spreads it. Measured on a fixed 40-gon at :math:`k = 3`,
meq converged at close to *second* order where a Solov'ev control on the very
same meshes converged at fourth. No amount of accuracy inside the elements
recovers that, because there is nothing wrong with the elements — the *domain*
is wrong.

This is precisely the difficulty the transfer technique exists to remove, and
it is worth knowing before anyone designs another fitted-polygon study.

How it works
------------

Given a level set for :math:`\Gamma`:

#. **Mark the background elements** lying inside :math:`\Gamma` and extract them
   as a submesh. That union is :math:`\Omega_h`, and its boundary is
   :math:`\Gamma_h` — inscribed in the plasma, since every element of
   :math:`\Omega_h` is entirely inside.
#. **Build transfer paths.** From each point on :math:`\Gamma_h`, a short path
   runs out to :math:`\Gamma`. meq uses a vertex-cone family, which is one of
   several the underlying library provides.
#. **Impose the transferred datum.** The value on :math:`\Gamma_h` is not the
   boundary condition itself; it is the boundary condition on :math:`\Gamma`
   *carried back along the path*, using the computed flux to integrate along
   it. Since :math:`q` is a solved unknown, that line integral is available at
   the flux's own order — this is the mixed method paying off again.

The essential property, and the reason the method is more than a heuristic, is
that **optimal order is retained even when** :math:`\mathrm{dist}(\Gamma_h,
\Gamma)` **is only** :math:`O(h)` :cite:p:`CockburnSolano2012`. Earlier
approaches of this kind needed the gap to shrink like :math:`h^{k+1}`, which in
practice means fitting the mesh after all. That assumption — usually called
P.1 — is what :doc:`adaptivity` has to work to preserve.

What comes off a transferred face
---------------------------------

Two things, and both matter:

* **The HDG stabilisation.** :math:`\tau` is zero on :math:`\Gamma_h`. Leaving
  it on was measured to lose one order at :math:`k = 1` and two at
  :math:`k = 2`.
* **The flux constraint.** The transmission condition is not imposed there;
  the transferred datum is.

The trace degrees of freedom on those faces stay in the essential list at zero,
which is why ``trace()`` does not report :math:`\varphi_h` — see the note at the
end of :doc:`formulation`.

.. warning::

   **The box must contain the shape with room to spare.** If :math:`\Omega_h`
   reaches the edge of the background mesh, part of :math:`\Gamma_h` is an
   inherited mesh boundary carrying no transferred datum, and zero would be
   imposed there instead — quietly, and the run would converge. meq checks for
   this and refuses rather than solving it, and it also refuses a shape that
   encloses no background element at all.

.. important::

   **A silently absent datum is the failure mode to guard against, because it
   converges.** :math:`\Gamma_h` carrying no transferred value simply imposes
   zero, which is a well-posed problem — just not the one asked for. The
   acceptance test for this path therefore pins the driver against the library
   on the same configuration **and** against a deliberately zero-datum solve,
   asserting that the two differ substantially. Agreement with the first alone
   would pass with the transfer doing nothing.

.. _curved-band:

The band, and what it means for output
--------------------------------------

Because :math:`\Omega_h` is the union of elements lying *inside* :math:`\Gamma`,
there is a band :math:`O(h)` wide that is inside the plasma and outside the
mesh. Nothing was solved there. Both output formats have to say something about
it, and they do so differently; see :ref:`output-band`.

.. note::

   The error estimator has to deal with the band too, and getting it wrong was
   nearly fatal to the adaptive loop. The term comparing the trace against
   :math:`\psi^\star` on :math:`\Gamma_h` was originally comparing against a
   trace pinned to *zero* rather than against the :math:`\varphi_h` actually
   imposed, so its size was the geometric gap — orders larger than every other
   term in the estimator, and converging at about half an order. The loop would
   have run, produced plausible pictures, and refined the wrong elements. See
   :ref:`adaptivity-eta5`.

Shapes
------

Two parametrisations ship, both configured under ``[boundary.shape]`` and both
implemented by :cpp:class:`meq::BoundaryShape`:

``Type = "miller"``
   The Miller D-shape :cite:p:`Miller1998` in physical parameters: major
   radius, minor radius, elongation, triangularity and squareness.

   .. warning::

      **Triangularity enters as** :math:`\arcsin\delta`, **not as**
      :math:`\delta`. The two differ by about one per cent at
      :math:`\delta = 0.35`, which is enough to be visible in a convergence
      study and not enough to look like a bug.

``Type = "mxh"``
   The Miller extended harmonic parametrisation
   :cite:p:`ArbonCandyBelli2021`, taking cosine and sine coefficient arrays.
   The cosine array **starts at** :math:`c_0`, the tilt; the sine array starts
   at :math:`s_1`, and there is no :math:`s_0`. Physically, :math:`s_1` is
   triangularity, :math:`-s_2` is squareness, :math:`c_0` is tilt and
   :math:`c_1` is ovality; the :math:`n`-th harmonic induces an
   :math:`(n-2)`-sided polygonal deformation.

Both refuse rather than approximate. A surface is rejected if it reaches or
crosses the axis — the operator's :math:`1/r` is not integrable through
:math:`r = 0` — and if it is **not star-shaped** about its own centre, since
the level set is evaluated by bisecting on the poloidal angle and that is only
legitimate when the polar angle increases monotonically. The message names the
angle at which monotonicity failed.

.. note::

   :cpp:func:`meq::BoundaryShape::levelSet` returns the **radial gap** — the
   distance from the centre minus the distance from the centre to the curve at
   the same polar angle — and *not* a signed distance, which is larger. It is
   negative inside, zero on the curve and positive outside, which is what the
   subdomain marking needs.
