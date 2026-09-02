Worked examples
===============

``examples/`` holds configurations that are run by the acceptance tests, so they
are known to work and known to stay working. Each carries a long header
explaining what it demonstrates and — usually more usefully — what it does
*not*.

Run them from the repository root, since two of them name profile files by
relative path:

.. code-block:: sh

   ./build/meq examples/soloviev-nstx.toml

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - File
     - Demonstrates
   * - ``soloviev-nstx.toml``
     - The headline Solov'ev benchmark on a rectangle. Linear source, one Newton
       step, fitted boundary. Every default written out explicitly.
   * - ``manufactured-driver.toml``
     - The Newton path end to end: a source that depends on :math:`\psi`
       linearly, quadratically *and* exponentially. Expect a short quadratic
       run.
   * - ``manufactured.toml``
     - The same source with ``[boundary] Type = "exact"`` — a **study**
       configuration the driver refuses. Kept as the definition of what the
       convergence tests measure.
   * - ``mhd-rectangle.toml``
     - The general tabulated-profile source: :math:`F` built from two files of
       numbers with nothing about the plasma hard-coded. **The only example in
       which Newton has real work to do** — its pressure gradient is quadratic
       in :math:`\psi`, so the problem is genuinely semi-linear.
   * - ``mhd-pprime.dat``, ``mhd-ggprime.dat``
     - Its two profiles, and the second full worked description of the
       tabulated file format.
   * - ``miller-curved.toml``
     - The curved boundary by transfer, on a background rectangle that knows
       nothing about the shape. :doc:`curved_boundary`.
   * - ``miller-adaptive.toml``
     - The full adaptive loop on the curved path, including the companion-mesh
       update. :doc:`adaptivity`.
   * - ``rotating-rectangle.toml``
     - Sonic toroidal rotation, two species, a tabulated density. The simplest
       complete rotating run. :doc:`rotation`.
   * - ``rotating-normalised.toml``
     - The same plasma with its profiles against normalised flux, so
       :math:`\psiax` is an unknown. :doc:`normalised_flux`.
   * - ``rotating-density.dat``, ``rotating-density-normalised.dat``
     - The same profile written both ways. The **only** documentation in
       ``examples/`` of the tabulated file format — and a demonstration of the
       trap in :ref:`profiles-file-format`.

.. _examples-benchmark-caveat:

The Solov'ev example is not the NSTX equilibrium
------------------------------------------------

Worth stating plainly, because the file is named after it and the temptation to
compare against the published closed form is strong.

For the published coefficients, :math:`\psi = 0` **is the separatrix** — a
curved contour through an X-point. ``soloviev-nstx.toml`` imposes :math:`\psi =
0` on a *rectangle* that circumscribes it, which is a different, well-posed
problem, and measured, the difference is the same order as :math:`\psi` itself.

Reproducing the published equilibrium needs :math:`\Gamma` to *be* that contour,
which is ``[boundary.shape]`` and the extension path. What this example
demonstrates is the driver and the discretisation on a well-posed problem.

.. note::

   This is also why meq's driver acceptance test pins the driver against the
   **library** on the same configuration rather than against a closed form: the
   question it answers is "does the driver reproduce what the library does",
   and a comparison against an analytic solution would answer a different
   question badly. The closed-form comparisons live in the convergence tests,
   where the boundary condition is the exact trace.

Two examples that deliberately test nothing about the Jacobian
--------------------------------------------------------------

``soloviev-nstx.toml`` and ``rotating-rectangle.toml`` both have
:math:`\partial F/\partial\psi \equiv 0` — the first because a Solov'ev source
is constant in :math:`\psi`, the second because its density profile is
**linear** in flux, so :math:`\mathrm{d}p/\mathrm{d}\psi` does not depend on
:math:`\psi`. Both converge in one Newton step, and both say nothing whatever
about the nonlinear solve.

That is deliberate in each case, and their headers say so. It is also a good
illustration of the warning in :doc:`profiles`: a profile linear in flux makes an
affine problem, and it is an easy way to write a test case that quietly measures
nothing.

``mhd-rectangle.toml`` is the counterexample, and it is worth reading beside
them. Its pressure gradient is *quadratic* in :math:`\psi`, so
:math:`\partial F/\partial\psi` genuinely depends on :math:`\psi` and Newton
converges quadratically over several steps rather than finishing in one.

.. note::

   That example carries a one-key control worth knowing about as a technique.
   Setting ``PPrimeScale = 0.0`` and changing nothing else collapses the source
   to its affine half, and Newton then finishes in a single step — **with a
   bit-identical initial residual**, because the pressure gradient vanishes at
   the boundary datum Newton starts from. The two problems are
   indistinguishable until the first step has been taken, which is as clean a
   demonstration as a configuration file can give of which profile supplies the
   nonlinearity.

   It also illustrates :ref:`sources-trivial-branch` from the other side: it is
   the *toroidal field* profile, not the pressure, that keeps that run off the
   trivial branch, since the pressure gradient vanishes exactly where the
   boundary condition puts :math:`\psi`.

Coverage gaps
-------------

For completeness, since a reader looking for an example of these will not find
one. No shipped example uses ``[mesh] File``, ``[boundary.shape] Type = "mxh"``,
``[initialguess]`` in any form, ``Strategy = "maximum"``, a tabulated ``Omega``
or ``Temperature``, or more than two species.

Those spellings come from :doc:`configuration` and from
``tests/unit/ConfigTests.cpp``, which does exercise all of them and is a good
secondary source when writing a new configuration.
