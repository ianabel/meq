Validation against an independent code
======================================

Everything in :doc:`testing` checks MEQ against closed forms, manufactured
solutions and its own refinement. Those catch a great deal, but they share one
blind spot: a convention MEQ has misread stays misread in the fixture that
checks it. The Solov'ev coefficients, the sign of :math:`\tau`, and
``DarcyForm``'s :math:`-q` were each settled by measurement only after an
argument said otherwise, and a self-consistent tree cannot find that class of
error at all.

**So MEQ is also checked against a code that shares the equation and essentially
no code.** ``freegs4e`` — a FreeGS derivative — solves the same
Grad–Shafranov equation by an entirely different route: free boundary through
von Hagenow Green's functions, second- and fourth-order finite differences on a
uniform :math:`(R, Z)` grid, and Picard iteration with adaptive blending. MEQ is
hybridized discontinuous Galerkin with Newton. They agree on
:math:`\Delta^* \psi = -F`, on :math:`F = \mu_0 r^2 p' + g g'`, and on
:math:`\psi` in Wb/rad, and on very little else.

The harness is ``tools/freegs4e-benchmark/``; its ``README.md`` carries the
per-experiment numbers and is the authoritative record.

What it establishes
-------------------

Seven tokamak configurations — a test tokamak in three profile shapes, plus
MAST, MAST-U, TCV and DIII-D — spanning aspect ratios 1.4 to 3.4 and
elongations 1.3 to 2.1, including one case with :math:`g g' < 0` throughout.

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - case
     - machine
     - relative :math:`L^2`
     - relative :math:`L^\infty`
   * - A
     - test tokamak, classic
     - 1.5e-04
     - 5.4e-04
   * - G
     - MAST-U
     - 6.4e-04
     - 1.6e-03
   * - B
     - test tokamak, :math:`ff'`-dominated
     - 7.1e-04
     - 1.5e-03
   * - E
     - test tokamak, diamagnetic
     - 2.7e-03
     - 5.8e-03
   * - F
     - DIII-D
     - 2.9e-03
     - 5.2e-03
   * - D
     - TCV
     - 4.4e-03
     - 7.9e-03
   * - C
     - MAST
     - 6.8e-03
     - 1.3e-02

Over roughly 11,000 grid nodes per case, worst case 0.7 %.

.. note::

   These numbers date from when the page was written and will drift. The
   harness regenerates them in about ten minutes; the reason to re-run it is
   any change to the source term, the boundary shape machinery or the profile
   reader, since those are what it actually exercises end to end.

How the comparison is posed
---------------------------

MEQ solves the **fixed**-boundary problem and ``freegs4e`` solves the **free**
one, so they cannot be handed the same input. The comparison instead runs in
one direction:

#. ``freegs4e`` converges a free-boundary equilibrium, with coils, an X-point
   and a plasma current constraint.
#. A closed interior flux surface of that answer is extracted and fitted to the
   Miller eXtended Harmonic form, which is what :doc:`configuration`'s
   ``[boundary.shape] Type = "mxh"`` takes.
#. Its profiles are converted onto MEQ's :math:`\psi` and written as MEQ
   profile tables.
#. MEQ solves the fixed-boundary problem on that surface, and the two fields
   are differenced on a common grid.

So MEQ is given ``freegs4e``'s own boundary and its own source, and asked to
reproduce the field inside. Nothing about MEQ's discretisation is shared with
the reference, and neither is the solver, the mesh, or the language.

.. warning::

   **The band must be dropped before any error norm.** MEQ's NetCDF output
   continues :math:`\psi` past :math:`\Gamma_h` by a Taylor step at nodes that
   were never solved for, and marks them in ``extrapolated``. They are not the
   same quantity another code computes there. See :doc:`output`.

Four conversions, each of which converges to a wrong answer
-----------------------------------------------------------

None of the following fails loudly. Each produces a solve that converges
cleanly to an equilibrium that is not the one intended, which is why they are
written down rather than left to the next person.

**The reference's profiles are** :math:`dp/d\psi`, **not** :math:`dp/d\psi_N`.
``freegs4e``'s own docstring says normalised flux and is wrong; the arrays are
against unnormalised :math:`\psi`. Dividing by
:math:`\psi_\mathrm{bnd} - \psi_\mathrm{ax}` on the docstring's word costs a
factor that for a tokamak is of order one — it does not look wrong. Only the
*derivative* column of a MEQ profile table carries a chain-rule factor, because
:doc:`profiles` uses a Hermite cubic and each row holds a value and its
derivative.

**Not the separatrix.** Every reference case is diverted, so its last closed
surface passes through an X-point and has a *corner* there. MXH is a truncated
Fourier series and cannot represent one. MEQ's own curved-boundary studies avoid
the separatrix for the same reason — see :doc:`curved_boundary`. An interior
surface fits two orders of magnitude better.

**Zero must land on the surface actually used.** MEQ's :math:`\psi` vanishes on
:math:`\Gamma`, so the profile abscissa has to be shifted by the flux of the
surface handed over, not by the flux of the separatrix. Getting this wrong puts
MEQ's boundary where the source has nearly died, and the solve converges in
three Newton iterations to a field more than twenty times too small.

**The initial guess selects which equilibrium is reported**, which is the next
section.

Root selection is the real difficulty
-------------------------------------

The fixed-boundary problem with :math:`p'` and :math:`g g'` given as functions
of :math:`\psi` is nonlinear, and it has **more than one solution**. On one of
the seven cases three were found:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - initial guess
     - :math:`\max \psi`
     - against the reference
   * - ramp, small amplitude
     - 1.93e-03
     - −95.7 %
   * - **seeded from the reference**
     - **4.525e-02**
     - **−0.067 %**
   * - ramp, larger amplitude
     - 5.225e-02
     - +15.4 %
   * - ``freegs4e``
     - 4.528e-02
     - —

The physical solution lies **between** the two a ramp can reach. No ramp
amplitude finds it: the sweep steps from the lower root to the upper one
without landing on it, which is the signature of the unstable middle branch of
an S-curve — Newton slides off it from either side. Seeded from the reference
through ``[initialguess] Type = "gridfunction"``, MEQ converges to it in **two
Newton iterations** and then agrees better than any other case.

.. important::

   A converged solve is not by itself evidence of the right equilibrium. Where
   a problem admits several, the initial guess is part of the problem statement
   rather than a convenience, and MEQ will report whichever one it is steered
   toward without complaint.

**Distinguishing a second solution from an error takes one measurement**:
refine. The upper root sat at +15.39 % to +15.41 % across two refinement levels
and two polynomial degrees, converged to six digits. A discretisation error
falls under refinement; a different solution does not.

What limits the agreement
-------------------------

The residual disagreement is most likely geometric rather than either solver's.
The boundary handed to MEQ is a *fit* to a contour that was itself extracted
from the reference's finite-difference grid, and the fit stops improving once
enough harmonics are used — which is the signature of the contour being the
limit rather than the fit. Refining the reference's grid is the way to test
that, and it has not been done.

The source conversion is *not* a candidate, and that is checked rather than
assumed: evaluating MEQ's profile tables on the reference's own :math:`\psi`
reproduces the reference's own current density to about 2e-05. Both codes are
solving the same equation.

Which field is compared
-----------------------

The NetCDF output carries the **post-processed** potential :math:`\psi^*`, not
the raw :math:`\psi_h` — see :doc:`postprocessing` and :doc:`output`. The
comparison reads that file, so every number on this page is :math:`\psi^*`.

That is the field a consumer of MEQ actually receives, which is why it is the
one compared. It is worth stating because :math:`\psi^*` converges at
:math:`k+2` rather than :math:`k+1`, so the cost table below should not be read
as a statement about the raw discretisation. The convergence studies in the test
suite measure :math:`\psi_h` instead, and those are the :math:`k+1` rates.

Cost
----

The two codes do not solve the same problem, so a wall-clock ratio is not a
statement about either: ``freegs4e`` converges a *free*-boundary equilibrium
with coils, an X-point and a control system, while MEQ solves the
*fixed*-boundary problem inside a surface it is handed. The comparison that is
language-neutral, and that actually asks something about the discretisation, is
**accuracy per unknown**.

Measured on one case, MEQ reaches the benchmark's accuracy floor with about
**five times fewer degrees of freedom at** :math:`k = 3` **than at**
:math:`k = 1`, and in about a fifth of the time — roughly 400 elements against
7,000. That is the high-order argument made concrete, and it is the reason
:doc:`formulation` puts the flux in the discretisation rather than
differentiating the potential.

.. note::

   **MEQ saturates this benchmark**, which is the most useful thing the cost
   table says. The error stops improving once it reaches the *reference's* own
   accuracy — set by the boundary fit and by the contour extracted from
   ``freegs4e``'s finite-difference grid — so refining MEQ further buys nothing.
   Future work here should refine the reference, not MEQ.

Read the timings as an order of magnitude only. This project's standing rule is
that a timing is a measurement about the machine it was taken on; the
reproducible column is the dof count.

**One convergence failure is worth knowing**, because it is a useful piece of
advice rather than a defect: at :math:`k = 2` on the coarsest mesh two of the
cases do not converge at all — Newton fails, the reactive ladder's
Picard-then-Newton fails too, and the driver says to raise the degree.
:math:`k = 3` on the *same mesh* then converges. See :doc:`nonlinear`;
`p`-refinement reaches cases that neither refinement in `h` nor a globalisation
does.

The free-boundary comparison
----------------------------

Everything above is the **fixed**-boundary rehearsal, and it contains a boundary
fit, so it cannot reach round-off however fine either grid is. The
free-boundary comparison removes the fit: MEQ is given the same four coils, the
same profile shapes, the same prescribed :math:`I_p` and the same limiter
contact, **and no plasma boundary at all**. It ships as
``examples/limited-tokamak.toml``.

The reference is ``freegs4e``'s ``H_limited_circular`` — vertical-field coils
only, so no X-point exists in range and the plasma boundary is the flux surface
through the limiter contact. That is the one topology MEQ's pointwise
plasma-support test can represent; see :doc:`normalised_flux`.

Against the reference converged on a :math:`513^2` grid, with MEQ at
:math:`k = 3` on **26,375 elements** and ten Gegenbauer modes:

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - quantity
     - ``freegs4e``, :math:`513^2`
     - MEQ
     - relative
   * - :math:`\psi_{\mathrm{ax}}`
     - 9.308752e-02
     - 9.307342e-02
     - 1.5e-04
   * - :math:`\psi_{\mathrm{bnd}}`
     - 2.622462e-02
     - 2.622580e-02
     - 4.5e-05
   * - profile scale
     - 0.939672
     - 0.9400254
     - 3.8e-04
   * - :math:`I_p`
     - 3.0e+05 A
     - 3.000003e+05 A
     - it is the constraint

Three points about reading that table.

**The reference is not exact either, and its own grid error is the floor.** At
:math:`513^2` it sits **7.92e-04** from its own Richardson limit, so agreement
closer than that is measuring the reference rather than MEQ. All three
quantities above are inside it.

**MEQ converges in its own mesh on this case.** Uniform refinement of the
shipped mesh at :math:`k = 3` gives :math:`\psi_{\mathrm{ax}}` = 9.484400e-02,
9.482652e-02, 9.482423e-02 on 1601, 6527 and 26,375 elements — order **2.93**,
with the finest 3.2e-05 from the extrapolated limit. (Those are against the
:math:`129^2` reference's coil currents, which the shipped file carries.)

**The profile scale is predicted, not fitted.** ``freegs4e`` rescales both
profiles by a factor :math:`L` to deliver the requested current, and MEQ's
tables carry the :math:`129^2` run's :math:`L`. The scale MEQ's own current
border should then report is
:math:`(L_{513}/L_{129})(\mathrm{span}_{513}/\mathrm{span}_{129}) = 0.939672`,
which it reproduces to 3.8e-04. So the scale is a third independent check and
not a free parameter.

Two inputs have to be taken from the same reference run being compared against,
and neither is obvious. ``freegs4e``'s **coil currents** are an output of its
control system, not an input, so they converge in the grid like everything else
— one of them moves 14 % between :math:`129^2` and :math:`513^2`. And its
**limiter contact** is the maximum over the innermost ring of grid *cells*
inside the wall, so it too moves with the grid: on this case it sits on the
outboard midplane at :math:`129^2` and on the inboard shoulder, 0.6 m away, at
:math:`513^2`. Mixing a reference's currents with another's contact compares two
different machines.

The diverted free-boundary comparison
-------------------------------------

The comparison above hands **both** codes the same limiter contact, on purpose:
``freegs4e``'s own contact is the maximum over a ring of grid cells and moves
0.6 m between :math:`129^2` and :math:`513^2`, so pinning it takes the
contact-finding out of the comparison.

The diverted case leaves it in. ``examples/diverted-tokamak-xpoint.toml`` is
``freegs4e``'s ``A_testtokamak_classic`` — the FreeGS worked example's geometry,
an up-down asymmetric double null whose four coil currents its control system
solved for — and MEQ is told **nothing** about where the null is beyond a seed
deliberately placed 7.1 cm away. It solves for the X-point as two more unknowns
of its Newton (``[boundary.xpoint]``, see :doc:`configuration`); ``freegs4e``
finds it by a critical-point search on a finite-difference field. So the
agreement in its **position** is a result rather than a precondition.

.. list-table::
   :header-rows: 1
   :widths: 24 28 28 20

   * - quantity
     - ``freegs4e``, :math:`129^2`
     - MEQ, :math:`k = 2`, 2642 el
     - apart
   * - X-point
     - ( 1.093144118, −0.603965084 )
     - ( 1.093103369, −0.603529197 )
     - **4.4e-04 m**
   * - :math:`\psi_{\mathrm{ax}}`
     - 8.271751445e-02
     - 8.266003630e-02
     - 6.9e-04
   * - :math:`\psi_{\mathrm{bnd}}`
     - 3.240412551e-02
     - 3.237931762e-02
     - 7.7e-04
   * - profile scale
     - 1, by construction
     - 9.985753e-01
     - 1.4e-03
   * - :math:`I_p`
     - 2.0e+05 A
     - 2.000000e+05 A
     - it is the constraint

On a diverted plasma :math:`\psi_{\mathrm{bnd}}` **is** the active null's flux,
so the second and third rows are two views of one disagreement rather than two
independent checks. The upper saddle carries 2.891019e-02, 3.5e-03 further out,
which is what makes this a single-active-null equilibrium.

Two things about reading that table.

**The floor is the profile fit, not either discretisation.** ``fgsref.py`` fits a
spline to its own analytic profile shape before solving, and the fit moves it —
1.7e-02 of the amplitude in :math:`ff'` — while MEQ's tables are the analytic
shape, because that extends below :math:`\Psi = 0` into the vacuum where a
tabulation on :math:`\psi_n \in [0, 1]` cannot. The two codes are therefore
solving sources that differ at the per-cent level, and agreement much tighter
than what is measured would be evidence of a shared mistake.

**The pointwise field comparison is dominated by the conductor model, and only
by it.** Over every comparable node the relative :math:`L^\infty` in :math:`\psi`
is 0.30, which read alone is a failure — and all of it is inside two coils.
``freegs4e``'s default ``Coil`` is an exact **filament**, a point source with a
logarithmic singularity; MEQ's conductors are rectangles meshed into the domain
carrying a uniform current density. Excluded with a 5 cm collar, 5478 nodes agree
at **5.3e-04** relative :math:`L^2` and 1.8e-03 :math:`L^\infty`, worst on the
outboard midplane near the plasma edge rather than anywhere near a conductor. The
benchmark's ``compare.py --exclude-box`` reports the excluded region as a row of
its own rather than dropping it.

What it does not establish
--------------------------

**Nothing about a found limiter contact under a free boundary.** In the limited
comparison the contact is prescribed to MEQ as a point rather than searched for
on the limiter curve, so the surface MEQ reports is the one through the point it
is given. Finding it on a meshed curve is a capability
(``[boundary.limiter] SurfaceAttribute``) and is checked against its own
reference, not against ``freegs4e``.

**Nothing about the conductor model.** The diverted comparison above measures
it and cannot resolve it: two of that machine's four coils are filaments in the
reference and meshed rectangles in MEQ, which is a modelling difference no
refinement of either code closes.

**Nothing about the flux-surface machinery.** Only :math:`\psi` on a grid is
differenced. The quantities in :doc:`surface_geometry` are checked against
their own references, not against ``freegs4e``, which computes no flux-surface
averages at all.

**Nothing about a real reconstruction.** Both codes are given analytic profile
shapes. Agreement here says the two implementations of one equation agree; it
says nothing about whether either matches a machine.

Running it
----------

``freegs4e`` is not part of MEQ and is not required to build or test it. The
harness needs a Python environment with it importable, and its ``README.md``
records the dependency pins that do and do not work, along with two defects in
that library worth knowing before trusting a result from it.

.. code-block:: sh

   cd tools/freegs4e-benchmark
   python3 -m venv venv
   venv/bin/pip install numpy scipy matplotlib h5py shapely numba netCDF4 contourpy
   PYTHONPATH=/path/to/freegs4e venv/bin/python fgsref.py     # the references
   venv/bin/python make_case.py out/*.npz                     # TOML + profiles
   for f in out/*.toml; do meq "$f"; done
   venv/bin/python compare.py out/*.npz

The six **shipped** fixed-boundary machine cases are built by the same two
scripts, from references refined a level further:

.. code-block:: sh

   cd tools/freegs4e-benchmark
   FGSREF_OUT=$PWD PYTHONPATH=/path/to/freegs4e \
       venv/bin/python fgsref.py --shaped --nx=257 --seed-from=auto
   sh make_fixed.sh                                    # -> examples/fixed-*

``--seed-from`` is what makes the refinement affordable: the grids are
:math:`2^n + 1` and therefore nest, so a :math:`257^2` run starts from the
committed :math:`129^2` answer and converges in a handful of Picard passes
rather than the twenty to a hundred a cold start costs. All nine machines take
about six minutes. The references at that resolution are **not** committed —
they are 1.6 MB each and nothing but the comparison needs them; the shipped
TOML carries the MXH coefficients and the profile tables, which is everything a
*run* needs. See :ref:`examples-machine-fixed-boundary`.

**The resolution is not a detail.** What floors the comparison is the MXH fit of
:math:`\Gamma`, and that fit is extracted from ``freegs4e``'s grid: at
:math:`129^2` with ten harmonics it reads 2.7e-04 to 7.6e-04 m, and at
:math:`257^2` with twenty it reads 4.4e-05 to 1.8e-04. Above twenty harmonics
nothing improves, which is what says the remaining residual is the grid rather
than the truncation.

Two pieces of the harness self-test, and both should be run first if anything
looks wrong: ``mxh.py`` fits a shape whose coefficients it knows and demands
them back, and ``convert.py`` checks its derivative column against a closed
form and carries a control that deliberately drops the chain rule.
