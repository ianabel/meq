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
   * - ``flux-surfaces.toml``
     - The fourth output file: the equilibrium reduced to flux surfaces and
       flux-surface averages against a flux label. Curved on purpose, so the
       band mask is not identically zero. :ref:`output-flux-surfaces`.
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
   * - ``free-boundary-halfdisc.toml``
     - The exterior coupling on its own: :math:`\Gamma` an artificial boundary
       in the vacuum, no limiter, no coils. The smallest free-boundary run.
   * - ``limiter-halfdisc.toml``
     - The limiter as a **curve** rather than a point — the contact is found on
       a meshed limiter surface rather than prescribed.
   * - ``limited-tokamak.toml``
     - A whole machine, free boundary: four coils, a prescribed :math:`I_p`, a
       limiter contact and the exterior map, reproducing ``freegs4e``'s own
       converged equilibrium. :doc:`validation`.
   * - ``diverted-tokamak-xpoint.toml``
     - The same, **diverted**: the X-point is two more unknowns of the Newton
       rather than a point the file names, and the plasma support gets an outer
       loop of its own. :doc:`validation`.
   * - ``diverted-tokamak.toml``
     - The same machine with the bounding point *prescribed* at the reference's
       X-point instead. Kept as the comparison that isolates the border.
   * - ``diverted-tokamak-generated.toml``
     - The same machine again, **meshing itself**: ``[mesh.generate]`` describes
       the half-disc and the conductors' rectangles are derived from the
       ``[[coils]]`` blocks, so each coil is written once. Run it with
       ``meq-run`` — see :doc:`running`.
   * - ``fixed-h-circular.toml``, ``fixed-a-testtokamak.toml``,
       ``fixed-e-diamagnetic.toml``, ``fixed-d-tcv.toml``,
       ``fixed-g-mastu.toml``, ``fixed-f-diiid.toml``
     - Six **machines posed fixed boundary**, from :math:`10^3` to
       :math:`10^5` elements, each with a real tabulated profile pair and an
       independent code's answer to check against. See below.

.. _examples-machine-fixed-boundary:

Six machines, posed fixed boundary
----------------------------------

**A real machine's equilibrium is the same equilibrium whichever boundary you
pose it on, and that is what these are.** ``freegs4e`` solves each machine
*free* boundary on a :math:`257^2` grid — its own conductors, its own inverse
solve for their currents, its own algorithm — and its
:math:`\psi_{\mathrm N} = 0.95` surface is fitted to an MXH parametrisation.
:math:`\psi` is defined only up to an additive constant, so subtracting that
surface's flux gives the **same** equilibrium with :math:`\psi = 0` on
:math:`\Gamma`. Inside it, that is the fixed-boundary Grad–Shafranov equation
with the machine's own :math:`p'` and :math:`gg'`, exactly — and the free
solution, shifted, is a reference for it.

``tools/freegs4e-benchmark/make_case.py`` is the generator and
``make_fixed.sh`` the driver; both are re-runnable and the TOML headers say so.

**Why** :math:`\psi_{\mathrm N} = 0.95` **and not the separatrix** is the one
decision in the recipe that is not mechanical, and it buys two things at once.
The last closed surface of a diverted machine passes through an X-point and has
a **corner**, which a truncated Fourier series cannot turn; and the source *at*
that surface is the plasma edge, where a profile exponent of 1.2 makes
:math:`p'` a fractional power and caps any convergence rate at about 1.2
whatever the polynomial degree. At 0.95 the boundary is smooth and closed and
the profiles are analytic across it. **These are therefore the only
machine-shaped problems in the tree on which a high-order claim can be made**,
and the claim is measured: :math:`k+2` in :math:`\psi^{*}`.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - File
     - Shape
     - What it is there for
   * - ``fixed-h-circular.toml``
     - :math:`\kappa` 1.01, :math:`\delta` 0.05
     - The nearly-circular **control**, and **limited**, so no X-point exists
       anywhere in it. The smallest, about 1100 elements.
   * - ``fixed-a-testtokamak.toml``
     - :math:`\kappa` 1.19, :math:`\delta` 0.32
     - The plain diverted case — and the same machine that
       ``machine-a-testtokamak.toml`` solves *free* boundary, so the two are a
       pair.
   * - ``fixed-e-diamagnetic.toml``
     - :math:`\kappa` 1.31, :math:`\delta` 0.35
     - Exactly up–down symmetric, and :math:`gg'` of the other sign.
   * - ``fixed-d-tcv.toml``
     - :math:`\kappa` 1.75, :math:`\delta` 0.23
     - The most elongated, and the smallest minor radius.
   * - ``fixed-g-mastu.toml``
     - aspect 1.56
     - Spherical — the lowest aspect ratio in the set.
   * - ``fixed-f-diiid.toml``
     - :math:`\kappa` 1.47, :math:`a` 0.59 m
     - The large conventional machine, near double null, about
       :math:`10^5` elements.

.. note::

   **Every one uses** ``[source] Normalised = true``, **and that is not a
   preference.** Against :math:`\psi` in Wb/rad the profile table necessarily
   stops at the axis, because ``freegs4e``'s profiles are defined on
   :math:`\psi_{\mathrm N} \in [0, 1]` and there is no :math:`\psi_{\mathrm
   N} < 0`; MEQ's documented out-of-range policy extends a table by a constant,
   so above the axis flux the source becomes a plateau — and that plateau
   carries a solution of its own. Sweeping only the initial guess on one
   machine reached 0.0105, then 1.2604, then 0.9996 of the reference's axis
   height as the ramp grew. With ``Normalised = true``, :math:`\Psi` is
   confined to :math:`[0, 1]` by the constraint :math:`\psi_{\rm
   ax} = \max\psi`, so the region that branch lived in does not exist.

.. note::

   **They carry no conductors, and their absence is a property of the problem.**
   :math:`\Gamma` is the plasma edge, so every coil of the real machine is
   outside the computational domain and its whole influence is in a Dirichlet
   datum that is identically zero. ``[conductors] Model`` would have nothing to
   act on. A fixed-boundary case that *can* exercise coil subtraction needs a
   domain larger than the plasma; ``coils-rectangle.toml`` is the small one that
   does.

**What they cost and what they agree to** is in ``MEASUREMENTS.md`` M-147, and
``tests/convergence/MachineFixedBoundary.cpp`` is the acceptance — it runs the
shipped files through the shipped binary, because what rots about an example is
the *file*.

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

   This is also why MEQ's driver acceptance test pins the driver against the
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
one. No shipped example uses ``[mesh] File``, ``Strategy = "maximum"``, a
tabulated ``Omega`` or ``Temperature``, or more than two species.

``[boundary.shape] Type = "mxh"`` and ``[initialguess] Type = "ramp"`` were on
that list until the six machine cases above were added, and both are now
exercised by all six.

Those spellings come from :doc:`configuration` and from
``tests/unit/ConfigTests.cpp``, which does exercise all of them and is a good
secondary source when writing a new configuration.
