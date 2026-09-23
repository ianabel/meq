Configuration reference
=======================

A run is described entirely by a TOML file. This page documents every table and
every key; :doc:`examples` walks through the shipped configurations.

.. code-block:: toml

   [mesh]
   RMin = 0.6
   RMax = 1.4
   ZMin = -0.6
   ZMax = 0.6
   NR = 8
   NZ = 8
   RefinementLevels = 1

   [discretisation]
   PolynomialDegree = 3

   [source]
   Type = "soloviev"
   A = -0.52

Three tables are **required**: ``[mesh]``, ``[discretisation]``, ``[source]``.
The rest — ``[boundary]``, ``[solver]``, ``[output]``, ``[initialguess]``,
``[adaptivity]``, ``[conductors]`` — are optional, and an absent table behaves
exactly like a present but empty one. ``[[coils]]`` is an *array* of tables
rather than a table, so "absent" there means no conductors at all.

.. note::

   **Configuration keys are** ``UpperCamelCase``, **and table names are
   lowercase.** That is a deliberate mismatch with the C++ naming convention,
   inherited from a sibling project, and it is not an oversight to be
   reconciled.

How the parser behaves
----------------------

.. important::

   **An unknown key is an error, and so is an unknown table.** The message names
   the nearest accepted spelling:

   .. code-block:: text

      MEQ: configuration error in 'run.toml', key 'mesh.RefinmentLevels':
      is not a key of [mesh]; did you mean 'RefinementLevels'? accepted keys
      are: RMin, RMax, ZMin, ZMax, NR, NZ, RefinementLevels, File

   This is the whole point of having a schema. A key that is silently ignored is
   a run that quietly did something other than what the file says.

Numbers accept either TOML spelling: ``RMin = 0`` and ``RMin = 0.0`` both work.
Counts do **not** — ``NR = 4.0`` is refused rather than truncated.

.. warning::

   That distinction cost a real bug once, and the obvious fix was worse than the
   bug. Reading every number as a float rejects ``RMin = 0``; but the natural
   remedy in the TOML library used here **returns the default** when a
   conversion fails, because a failed conversion is indistinguishable from a
   missing key. That turns a loud failure into a silent wrong answer in the
   configuration, which is the last place anybody wants one. MEQ therefore
   checks the node's type explicitly.

Every error is a :cpp:class:`meq::ConfigError`, carrying the file name and the
fully-qualified key — including array elements, as
``source.species[2].Mass``. **Construction either succeeds and leaves every
accessor meaningful, or throws**; there is no partially valid configuration.

``[mesh]``
----------

The background mesh: either a box that MEQ triangulates, or a file.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``RMin``
     - *required*
     - Lower radial bound. Must not be negative — :math:`R` is a cylindrical
       radius, and a box reaching :math:`R = 0` contains the operator's
       singularity.
   * - ``RMax``
     - *required*
     - Upper radial bound; must exceed ``RMin``.
   * - ``ZMin``, ``ZMax``
     - *required*
     - Vertical bounds; ``ZMax`` must exceed ``ZMin``.
   * - ``NR``, ``NZ``
     - ``1``
     - Cells across each direction **before** refinement. Cells are split into
       triangles, so the initial element diameter is the cell diagonal.
   * - ``RefinementLevels``
     - ``0``
     - Levels of uniform refinement; each halves :math:`h`. Applies whether the
       mesh came from a box or a file.
   * - ``File``
     - ``""``
     - Any mesh format MFEM reads. Non-empty selects file mode, in which the box
       keys and ``NR``/``NZ`` are not used.

.. note::

   **The gridded output takes its extent from the mesh itself when the mesh
   comes from a file**, since ``RMin``…``ZMax`` describe a box that was not
   built. That is also the right answer: the mesh's extent is exactly the
   region the solve claims anything about.

   It used to be degenerate, and a run that had already solved was then lost at
   the output stage — which is not a corner case for free boundary, where the
   half-disc reaching the axis cannot come from ``MakeCartesian2D`` and the mesh
   is *always* a file. See ``tools/mesh/README.md``.

.. warning::

   ``File`` **still gives up something on the curved path.** The box bounds are
   what set the search length for the transfer paths of ``[boundary.shape]``,
   and a mesh file supplies none — so a file plus a shape is not a configuration
   to reach for. A file with no shape, which is the free-boundary case, is fine.

``[mesh.generate]``
-------------------

The mesh as **this file's own build product** rather than as an input. With
this block present, ``File`` is where a generator *writes* and where the solve
then reads, and ``meq-run`` — :doc:`running` — is the one command that does
both.

MEQ links MFEM and not gmsh, deliberately (``tools/mesh/README.md`` has the
reasoning; the short form is that the coupling between MEQ and a mesher is a
*file*). So ``meq`` does not run the generator. What it does is print the
command the block describes:

.. code-block:: sh

   ./build/meq --mesh-command examples/diverted-tokamak-generated.toml

and **refuse to solve** such a configuration unless the caller passes
``--mesh-ready`` to say the mesh has been made from it. Without that refusal an
edited geometry is answered from the previous geometry's mesh, at full order,
with every printed number looking exactly as it should.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Tool``
     - *required*
     - Which generator. ``"halfdisc"`` is the only one — ``tools/mesh/halfdisc.py``,
       a semicircle reaching :math:`R = 0` exactly with the conductors
       fragmented in. Naming it is what makes the block recognisable at all.
   * - ``Radius``
     - *required*
     - The **disc's** radius, metres. Not :math:`\Gamma` — see the warning below.
   * - ``Size``
     - *required*
     - Background element size, metres.
   * - ``Order``
     - ``1``
     - Geometric order. Above 1 the arc's mid-edge nodes are placed on the true
       circle rather than on the chord, and the axis stays exact.
   * - ``CoilSize``
     - *the background* ``Size``
     - Element size inside the conductors.
   * - ``PlasmaRMin``, ``PlasmaRMax``, ``PlasmaZMin``, ``PlasmaZMax``
     - *none*
     - A box to refine inside, in the same four-bounds form ``[mesh]``'s own box
       uses. All four go together with ``PlasmaSize``.
   * - ``PlasmaSize``
     - *none*
     - Element size in that box.
   * - ``LimiterR``, ``LimiterZ``, ``LimiterRadius``
     - *none*
     - A circular limiter **fragmented into** the geometry, so its edges are mesh
       faces, written as element attribute 20 — which is what
       ``[boundary.limiter] SurfaceAttribute`` reads.
   * - ``Vessel``
     - *none*
     - A closed vessel polygon, alternating :math:`R` and :math:`Z` in metres,
       **fragmented into** the geometry so its edges are mesh faces. Everything
       inside :math:`\Gamma`, outside this and not a conductor takes element
       attribute **30**, which ``[source] ExcludeAttributes`` can then name as a
       region that can never be plasma. At least three points, every
       :math:`R \ge 0`, and a non-zero area. It is **not** confinement:
       ``psi_bnd`` is what confines the plasma, and this covers the one case
       connectivity cannot settle — several O-points across a saddle.
   * - ``Transition``
     - *four background sizes*
     - Width of the graded transition out of a refined region, metres.
   * - ``Check``
     - ``true``
     - Re-read the written file and assert MEQ's preconditions on it: that
       :math:`R` reaches 0 **exactly**, that :math:`\Gamma` and the axis are the
       outer boundary and nothing else, and that each coil attribute covers its
       rectangle.
   * - ``Symmetric``
     - ``false``
     - Mesh :math:`z \ge 0` and **reflect it**, so the mesh is exactly
       mirror-symmetric about :math:`z = 0`. This is what makes
       ``[solver] UpDownSymmetry`` usable: that projection averages every
       degree of freedom with the one at its own reflection, and gmsh's
       triangulation of a symmetric geometry is *not* symmetric — it picks a
       diagonal and picks freely. The generator **refuses** this on a geometry
       that is not itself mirror-symmetric — an unpaired conductor, a limiter
       off the midplane, a vessel outline that is not its own image — rather
       than reflecting your machine into a different one. ``Plasma*`` is the
       exception and may be asymmetric: it is a size field, not geometry, so
       the meshed half is refined as asked and the other half gets the
       reflection.

.. note::

   **The conductors are not repeated here, and that is the point.** The
   generator's ``--coil`` rectangles are derived from the ``[[coils]]`` blocks,
   in file order — which is the order ``halfdisc.py`` assigns its ``10 + i``
   element attributes in. So a machine's conductors are written once.

   That is not tidiness. FB-2 measured what meshing *to* a conductor is worth
   against cutting through it — rates of 1.99 / 2.88 / 3.01 aligned, against
   1.33 / 1.27 / 1.09 cut — so a coil the mesh is not aligned to costs a full
   order, silently. Deriving one list from the other makes that unreachable rather than
   unlikely, and it aligns the mesh to the rectangle ``meq::Coil``'s quadrature
   uses *to the last bit*, which a hand-copied command line cannot be.

   The two conventions differ and MEQ's is what the file carries: ``[[coils]]``
   is a centre and half-extents, ``halfdisc.py`` takes a corner and two extents,
   and the refined box is four bounds here and a corner plus extents there. The
   driver converts. One file does not get to hold two meanings of four numbers.

.. warning::

   **``Radius`` is the disc and not** :math:`\Gamma`, and it is easy to read
   them as one number. :math:`\Omega_h` is cut *from* the generated mesh as the
   elements lying inside ``[boundary.exterior] Radius``, so the arc gmsh draws
   is the background's outer edge and :math:`\Gamma` is the smaller semicircle
   inside it. The two are checked against each other at parse time, before gmsh
   runs.

.. note::

   Two more things become parse errors once the mesh's geometry is in the file,
   each of which otherwise costs a mesh generation to discover: a
   ``[boundary.limiter] SurfaceAttribute`` naming a region this mesh will not
   carry — which at the solve reads *"psi_bnd = max psi_h over the empty set"* —
   and a limiter circle through the axis.

``[discretisation]``
--------------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``PolynomialDegree``
     - *required*
     - The degree :math:`k` of all three HDG spaces. The source papers report
       :math:`k = 1` to :math:`5`; beyond about 5, or eight refinements,
       round-off dominates and a table that flattens there is behaving
       correctly.
   * - ``Tau``
     - ``1.0``
     - The HDG stabilisation, dimensionless and positive. Optimal order needs
       only :math:`\tau = O(1)`. See :ref:`formulation-tau` before changing it.

``[source]``
------------

``Type`` is required, and the accepted keys depend on it — so ``Normalised``
under a Solov'ev source is an unknown key, not an ignored one.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - ``Type``
     - 
   * - ``"soloviev"``
     - :math:`F = -\left((1-A)R^2 + A\right)`, independent of :math:`\psi`. The
       problem is linear and Newton converges in one step.
   * - ``"mhd"``
     - :math:`F = \mu_0 R^2 p'(\psi) + (gg')(\psi)`, from two tabulated
       profiles.
   * - ``"manufactured"``
     - The nonlinear manufactured solution of
       :cite:t:`SanchezVizuetSolano2019`, whose :math:`\psi`-dependence is
       linear, quadratic and exponential.
   * - ``"rotating"``
     - A rotating multi-species plasma. See :doc:`rotation`.

``Type = "soloviev"``
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``A``
     - *required*
     - Dimensionless, with the flux normalised so that :math:`A + C = 1`. The
       NSTX benchmark uses ``-0.52``.

``Type = "mhd"``
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``PPrimeFile``
     - *required*
     - Path to tabulated :math:`\mathrm{d}p/\mathrm{d}\psi`. See
       :ref:`profiles-file-format`.
   * - ``GGPrimeFile``
     - *required*
     - Path to tabulated :math:`g\,\mathrm{d}g/\mathrm{d}\psi` — what EQDSK
       calls ``FF'``. **An alternative to** ``SafetyFactorFile``; naming both is
       refused.
   * - ``SafetyFactorFile``
     - —
     - Path to a tabulated **target** :math:`q(\Psi)`, which makes
       :math:`g\,\mathrm{d}g/\mathrm{d}\psi` an **output** of the run rather
       than an input. Requires ``Normalised = true`` and ``ToroidalFieldGuess``.
       See :ref:`configuration-q-driven`.
   * - ``SafetyFactorDegree``
     - ``2``
     - Degree of the :math:`g^2` polynomial the outer loop solves for. Only with
       ``SafetyFactorFile``.
   * - ``ToroidalFieldGuess``
     - —
     - **Required with** ``SafetyFactorFile``. The constant
       :math:`g = R B_\phi` the loop opens at — a machine's vacuum
       :math:`R_0 B_0`.
   * - ``PPrimeScale``, ``GGPrimeScale``
     - ``1.0``
     - Constant multiplying the table as read, so unit conversion needs no edit
       to the file.
   * - ``Mu0``
     - SI :math:`\mu_0`
     - Vacuum permeability on the :math:`R^2 p'` term. Set to 1 for a run in
       normalised units.
   * - ``Normalised``
     - ``false``
     - The tables are functions of :math:`\Psi = \psi/\psiax` rather than of
       :math:`\psi`. See :doc:`normalised_flux`.
   * - ``PsiAxis``
     - —
     - **Required when** ``Normalised = true``, **refused otherwise.** A
       starting guess for the Newton iteration, not a scale factor. Must be
       finite and non-zero.
   * - ``ConfineToPlasma``
     - ``false``
     - :math:`F` and :math:`\partial F/\partial\psi` are **zero outside the
       plasma**, so the plasma's *support* moves with the solution instead of
       being the whole domain. The plasma is the **connected** region of
       :math:`\Psi > 0` containing the magnetic axis, found by a flood fill over
       the mesh's element adjacency: a level of :math:`\psi` can cut the domain
       into several pieces, and the ones that do not contain the axis are
       pockets rather than plasma. **Refused unless** ``Normalised = true``,
       since the test is on :math:`\Psi`. See the warning below before setting
       it.
   * - ``ExcludeAttributes``
     - ``[]``
     - Mesh **element attributes** that can never be plasma, whatever the flux
       says there — the far side of a vessel wall, a port, a pocket the mesh
       carries for the coils' sake. ``psi_bnd`` already confines the plasma and
       the connectivity fill already separates the lobes of :math:`\Psi > 0`
       that a level set leaves joined; what neither can do is know that a lobe
       is *behind a wall*, which matters where several O-points sit across a
       saddle. An attribute rather than a polygon because the support is
       re-decided on every residual evaluation and a vessel does not move.
       ``[mesh.generate] Vessel`` is how to produce one, as attribute ``30``.
       **Refused unless** ``Normalised = true``, and an attribute this mesh does
       not carry is **refused** rather than ignored: a silently empty exclusion
       would converge, report nothing unusual, and describe a machine with a
       current channel behind its own wall.

.. note::

   Profile paths are resolved against the **working directory of the run**, not
   against the directory the configuration file lives in.

.. warning::

   ``ConfineToPlasma = true`` **requires a profile that vanishes at the plasma
   edge**, :math:`p'(0) = 0`, and nothing can check that for you — the profile
   is a table.

   With :math:`p'(0) \ne 0` the source *jumps* across the plasma boundary, so
   the assembled residual is discontinuous in the unknowns and there is no
   Jacobian to iterate with. Measured, it does not converge at any degree, on
   any mesh, from the exact solution, or under ``PicardThenNewton``.

   Note also that :math:`\psi_{\mathrm{bnd}}` is zero unless a limiter is
   given.  With ``[boundary.limiter]`` present it is an unknown of the same
   bordered Newton, pinned by :math:`\psi_h` **at the contact point given**,
   and the plasma edge is then *found* rather than pinned at
   :math:`\psi = 0`.

.. note::

   ``ProfileFile``, and ``PPrimeVariable``/``PPrimeFit``/``GGPrimeVariable``/
   ``GGPrimeFit`` beside it, are **reserved** for reading several profiles out
   of one NetCDF file. That is not implemented, and all five are *refused* with
   a message saying so rather than accepted and ignored. The same holds of the
   corresponding keys on every other source type.

``Type = "manufactured"``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``R0``
     - *required*
     - The radial **offset** in :math:`\sin(K_r(R + R_0))`. **Not a major
       radius.**
   * - ``Kr``, ``Kz``
     - *required*
     - Radial and vertical wavenumbers.

``Type = "rotating"``
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``species``
     - *required*
     - An array of tables, written ``[[source.species]]``. Between 2 and
       :cpp:var:`meq::maxSpecies` of them.
   * - ``Omega`` / ``OmegaFile``
     - ``0``
     - Rotation frequency :math:`\omega(\psi)`, as a constant or a table.
       **Giving neither means no rotation**, and the source reduces to the
       static equation.
   * - ``GGPrime`` / ``GGPrimeFile``
     - *required*
     - :math:`g\,\mathrm{d}g/\mathrm{d}\psi`, as a constant or a table.
   * - ``OmegaScale``, ``GGPrimeScale``
     - ``1.0``
     - Scales for whichever form was given.
   * - ``ReferenceRadius``
     - ``1.0``
     - **The gauge**: the radius at which the electrostatic potential vanishes,
       and therefore at which each ``Density`` is the physical density. A
       constant radius — not the magnetic axis, not a flux-surface average. See
       :doc:`rotation`.
   * - ``Mu0``, ``Normalised``, ``PsiAxis``, ``ConfineToPlasma``
     - as ``"mhd"``
     - 

``[[source.species]]``
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Name``
     - ``"species<i>"``
     - Used in output variable names.
   * - ``Mass``
     - *required*
     - Particle mass in kilogrammes.
   * - ``Charge``
     - *required*
     - :math:`Z_s`: **signed and dimensionless**, not coulombs.
   * - ``Temperature`` / ``TemperatureFile``
     - *required*
     - :math:`T_s` **in joules**, as a constant or a table.
   * - ``TemperatureScale``
     - ``1.0``
     - Scale on whichever was given — this is how a file stays readable, with
       ``1.0`` keV and the conversion in the scale.
   * - ``Density`` / ``DensityFile``
     - *required unless* ``Neutralising``
     - :math:`n_{s0}` in :math:`\mathrm{m}^{-3}`, **on the curve**
       :math:`R = \texttt{ReferenceRadius}`.
   * - ``DensityScale``
     - ``1.0``
     - 
   * - ``Neutralising``
     - ``false``
     - This species' density is *derived* from the others by charge neutrality.
       **Exactly one species must set it**, and that species must not carry a
       density.

Cross-key rules, each with its own message: a constant and a file for the same
profile is an error; a ``*Scale`` with nothing to scale is an error; the species
must carry charges of both signs; and ``*Variable``/``*Fit`` keys are reserved
for a facility that does not exist yet and are refused rather than ignored, on
this source type and on every other.

``PlasmaCurrent`` — the third border
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``[source] PlasmaCurrent`` makes the profile **scale** an unknown of the bordered
Newton and prescribes the current instead: the profiles then give the current's
**shape** and this gives its **size**. It is how every production free-boundary
code poses the problem. Absent, the amplitude is fixed and the current is
whatever it comes out as.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``PlasmaCurrent``
     - *absent*
     - :math:`I_p` in **amperes**, signed. Requires ``Normalised = true``.

.. note::

   **It is in amperes here and** :math:`\mu_0 I_p` **at the library, and the
   difference is deliberate.** :cpp:func:`meq::GradShafranovSolver::setPlasmaCurrent`
   takes :math:`\mu_0 I_p` because everything inside the solver already does —
   Ampère's law reads the flux integral as :math:`-\mu_0 I_p` and the constraint
   is assembled as :math:`\int F/R`, which *is* :math:`\mu_0 I_p` — and a solver
   taking amperes would need a :math:`\mu_0` of its own, which could disagree
   with the source's and scale two terms of one equation differently.

   The configuration layer has no such problem: the file names exactly one
   :math:`\mu_0`, under ``[source]``, and it is the same one ``[[coils]]`` uses
   for its ``Current``. The conversion happens once, in the driver.

An explicit ``PlasmaCurrent = 0.0`` is **refused** rather than read as the
default: a prescribed current of zero is satisfied by driving the profile scale
to zero, which is a converged solve describing no plasma.

.. warning::

   **With one constraint and one scale the counting is right and the solution is
   not unique.** ``freegs4e`` fixes *two* quantities against two parameters — the
   axis pressure **and** :math:`I_p`. MEQ's border fixes one, so Newton has no
   reason to prefer the physical branch over a larger, flatter plasma carrying
   the same current. The initial guess is what selects it; see
   ``[initialguess]``.

``[boundary]``
--------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Type``
     - ``"zero"``
     - ``"zero"`` imposes :math:`\psi = 0` — the fixed-boundary problem proper.
       ``"exact"`` imposes the source's known exact solution and is a
       convergence-study device.

.. warning::

   ``Type = "exact"`` **is refused by the driver**, with an explanation and exit
   code 1. :cpp:class:`meq::Source` carries :math:`F` and
   :math:`\partial F/\partial\psi` and no closed-form solution — those live in
   the test fixtures, where they are the thing being converged against. It
   remains usable from the library. See :ref:`running-refusals`.

``[boundary.shape]``
~~~~~~~~~~~~~~~~~~~~

Present, this selects the **curved** path: see :doc:`curved_boundary`.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Type``
     - ``"none"``
     - ``"none"``, ``"miller"`` or ``"mxh"``.
   * - ``R0``
     - *required unless* ``"none"``
     - Major radius.
   * - ``Z0``
     - ``0.0``
     - Centre height.
   * - ``MinorRadius``
     - *required unless* ``"none"``
     - 
   * - ``Elongation``
     - ``1.0``
     - :math:`\kappa`.
   * - ``Triangularity``
     - ``0.0``
     - ``"miller"`` only. :math:`\delta`, which **enters as**
       :math:`\arcsin\delta`.
   * - ``Squareness``
     - ``0.0``
     - ``"miller"`` only. Zero gives the original three-parameter Miller shape.
   * - ``CosCoefficients``
     - ``[]``
     - ``"mxh"`` only. :math:`c_0, c_1, \ldots` — **starts at** :math:`c_0`,
       the tilt.
   * - ``SinCoefficients``
     - ``[]``
     - ``"mxh"`` only. :math:`s_1, s_2, \ldots` — **starts at** :math:`s_1`;
       there is no :math:`s_0`.

Each shape refuses the other's keys rather than ignoring them, and ``"mxh"``
with no harmonics at all is refused with a message suggesting ``"miller"``,
since that is an ellipse.

``[solver]``
------------

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key
     - Default
     - Meaning
   * - ``NewtonMaxIterations``
     - ``20``
     - 
   * - ``NewtonRelativeTolerance``
     - ``1.0e-8``
     - Relative to the residual at the **cold** iterate — see
       :ref:`running-warm-start` for why that matters.
   * - ``NewtonAbsoluteTolerance``
     - ``1.0e-12``
     - 
   * - ``PlasmaSupportSweeps``
     - ``0``
     - Freeze the plasma support within each solve and re-decide it between
       them, at most this many times. ``0`` leaves it moving inside Newton.
       Needs ``[source] ConfineToPlasma = true``; see the note below.
   * - ``XPointMeritWeight``
     - ``1.0``
     - How heavily the X-point's two rows count in the **line search's merit**,
       as a multiplier on the length :math:`R h` that converts :math:`q` into a
       flux. It changes the merit and nothing else — the border still solves
       :math:`q_r = q_z = 0`, so a converged answer is the same answer at any
       weight, and what moves is how many iterations it costs. **There is no
       good universal value and the sensitivity inverts between machines**: one
       benchmark case goes from 56 iterations to 35 at weight 20 and fails
       outright at 30, while another goes from 14 to 82 at weight 4. Measure
       your own case and do not carry a value between them. Refused without
       ``[boundary.xpoint]``, and refused if non-positive.
   * - ``AssemblyMode``
     - ``"threaded"``
     - ``"serial"`` or ``"threaded"``. Who computes the element-local work.
       See :ref:`linear-threading`.
   * - ``TraceSolver``
     - ``"pardiso"``
     - ``"umfpack"`` or ``"pardiso"``. Which direct solver factorises the
       hybridized trace system. The default falls back to ``"umfpack"`` on a
       build without oneMKL — a value the file *states* is refused there, one it
       merely inherits is downgraded. See :ref:`linear-trace-solver`.
       ``"cudss"`` parses but the driver refuses it — see the note below.

.. note::

   **What** ``PlasmaSupportSweeps`` **is for.** With ``ConfineToPlasma = true``
   the set of elements carrying :math:`F` is a functional of the iterate: the
   pointwise :math:`\Psi > 0` test moves with :math:`\psi`, and so does the
   connected component the flood fill reaches. Its derivative is a surface term
   on a moving edge and the Jacobian does not carry it.

   On a **limited** machine that costs nothing measurable and the default of
   ``0`` is right. On a **diverted** one it is the difference between a solve
   and a failure: measured with one key changed and nothing else, the first
   solve of ``examples/diverted-tokamak-xpoint.toml`` stalls at the 200-iteration
   cap with the support moving, and converges in 14 steps with it held. So the
   support is fixed within a solve and re-decided between solves — freeze,
   solve, refreeze at the answer, solve again — and the loop stops when a sweep
   changes nothing.

   The run reports how many sweeps it took and whether it **settled**, and
   writes both into the ``.nc`` as ``plasma_support_sweeps`` and
   ``plasma_support_settled``. A loop that hit the cap has written the solution
   of the problem its last sweep posed, which is not the same object as a
   converged equilibrium; treat it as a failure to converge rather than as a
   loose tolerance. Alternation between two supports is possible and the cap is
   what bounds it.

   The key is **refused without** ``ConfineToPlasma``, because there is then no
   support to freeze and the loop would re-solve the identical problem until the
   cap.

.. note::

   **Neither of the last two may change the answer, and that is what makes them
   safe to expose.** The two assembly modes are asserted **bit for bit** against
   each other, on a linear source and on a nonlinear one; the three trace
   solvers agree to about 1e-14. So the same configuration file with a different
   value for either must produce the same equilibrium, and if it does not, that
   is a defect rather than a tuning outcome.

   They are also the only keys in the file whose validity depends on how MFEM
   was **built**. A misspelling is refused when the file is parsed; a value
   naming a solver or a threading mode this binary does not have is refused at
   startup, before any mesh is built, with a message naming the CMake option.
   Two different faults, two different messages — MEQ never substitutes a
   solver you did not ask for, precisely because all three reach the same answer
   and the substitution would be invisible.

.. note::

   ``TraceSolver = "cudss"`` is a valid spelling and the **library** supports it,
   but the **driver refuses it** even on a build that has cuDSS — and it is
   withheld rather than merely unimplemented.

   A device solver is only worth having if the data **stays** on the device.
   MEQ's element-local integrators and its scatter into the trace matrix — most
   of a Newton step between them — have no device kernels yet, so a device trace
   solve would copy the system across the bus once per iteration in order to
   accelerate one part of it. MFEM's own HDG device-offload plan reaches the same
   conclusion about doing this group on its own: *"doing only those is worse than
   doing nothing … plausibly slower than staying on the host throughout"*. That
   work is under construction upstream, and this key will open when it lands.

   There is an immediate failure too, which is what made the refusal urgent
   rather than only principled: without an ``mfem::Device`` configured, cuDSS
   does not fall back. It reads host pointers as device pointers and aborts
   inside CUDA, with a message that says nothing about the key that caused it.

.. warning::

   ``AssemblyMode = "threaded"`` with ``OMP_NUM_THREADS=1`` and
   ``MKL_NUM_THREADS`` greater than one is a **pathological** combination and
   MEQ warns about it at startup. A one-thread OpenMP team does not get MKL's
   nested-region suppression, so the element-local dense work pays full MKL
   threading on every call: measured 12.4 s against 0.069 s at :math:`k = 3`.
   Either raise ``OMP_NUM_THREADS`` or set ``MKL_NUM_THREADS=1``.

.. note::

   ``LinearMaxIterations`` and ``LinearTolerance`` were once accepted here and
   are now **refused**, with a message saying why. They configure an *iterative*
   inner solve, and MEQ solves the hybridized trace system with a **direct**
   solver — UMFPACK, PARDISO or cuDSS — which has no iteration count and no
   tolerance to set.

   They were parsed and validated and read by nothing, which is worse than
   either accepting or rejecting them: a key that validates is a key its author
   believes is doing something. If you have a configuration file carrying
   either, delete the line.

``[[coils]]``
-------------

Poloidal field coils of rectangular cross-section, each carrying a **uniform**
current density. Written as an *array of tables* — one ``[[coils]]`` block per
conductor, in any order — and the set may be empty, which is what every
fixed-boundary configuration has.

A coil adds to the Grad–Shafranov source

.. math::

   F_{\mathrm{coil}}(R, z) = \mu_0\, R\, \frac{I_k}{|\Omega_{c_k}|}

inside coil :math:`k` and exactly zero outside it, summed over the coils
containing the point. The coil current is **data**: it does not depend on
:math:`\psi`, it contributes nothing to the Jacobian, and a solve driven by
coils alone is affine and finishes in one Newton step.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Name``
     - ``"coil<i>"``
     - Diagnostics only. It is what a refusal or a mesh-coverage warning quotes
       back, so it is worth writing.
   * - ``CentreR``, ``CentreZ``
     - *required*
     - The centre, in metres.
   * - ``HalfWidth``, ``HalfHeight``
     - *required*
     - Half-extents in :math:`R` and :math:`z`, metres. Strictly positive, and
       ``CentreR - HalfWidth`` must be strictly positive too — a coil reaching
       the axis is refused, because the operator's :math:`1/R` is not integrable
       through :math:`R = 0`.
   * - ``Current``
     - *one of the two*
     - The **total** current through the cross-section, amperes. Signed, and
       zero is allowed. A real winding's turns are not modelled, so this is
       turns × amps per turn.
   * - ``CurrentDensity``
     - *one of the two*
     - The uniform :math:`j_\phi` in A/m², multiplied by the area
       :math:`4\,\texttt{HalfWidth}\,\texttt{HalfHeight}`.
   * - ``FilamentsR``, ``FilamentsZ``
     - *from* ``FilamentSize``
     - This block's own filament stack, under ``[conductors] Model =
       "filament"`` **only**. Both or neither. See
       :ref:`dividing-a-rectangle`.

.. important::

   **Give exactly one of** ``Current`` **and** ``CurrentDensity``. Naming both is
   refused rather than resolved by precedence: an author who writes both has two
   numbers in mind, and silently honouring one of them is how a coil set ends up
   carrying a current nobody chose. Naming neither is refused for the same
   reason.

.. note::

   **There is no** ``Mu0`` **key on a coil block, deliberately.** The coil term
   and the plasma term are *added*, so they must share a permeability, and the
   coils take ``[source] Mu0``. A run in normalised units setting
   ``[source] Mu0 = 1`` while the coils kept the SI value would be summing two
   terms scaled a million-fold apart — and would converge, at full order, to a
   machine nobody described. Two keys that must agree are a way of writing down
   a disagreement.

.. warning::

   **A coil the mesh does not reach contributes nothing, silently.** :math:`F`
   is assembled by quadrature over the elements, so a coil outside the ``[mesh]``
   box is never sampled: the run converges, writes its files, and describes a
   machine with that conductor switched off. MEQ prints a warning naming the
   coil and both boxes, because nothing in the output could otherwise show it —
   the current appears in the ``.nc`` file's ``coil_current`` attribute whether
   or not it did any work.

The source is **discontinuous at every coil edge**, which is physics rather than
a modelling shortcut: a conductor has a boundary. It costs the convergence
*rate* wherever an element straddles that edge, exactly as a re-entrant corner
does, and the remedy is the same — put element edges on the coil edges, or
accept the rate. A study measuring an order should align the mesh; nothing else
needs to.

``examples/coils-rectangle.toml`` is the worked example, and it writes its two
currents one each way so that both spellings are exercised.

``[conductors]``
----------------

*How* the ``[[coils]]`` blocks enter the equation. One key decides it, and the
default is what every file written before this table existed already means.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Model``
     - ``"meshed"``
     - ``"meshed"``, ``"subtracted"`` or ``"filament"``. See below.
   * - ``QuadratureOrder``
     - the library's
     - Gauss points per direction per rectangle, for ``"subtracted"`` **only**.
   * - ``FilamentSize``
     - one per block
     - Target filament **cell** size in metres, for ``"filament"`` **only**.
       See below.
   * - ``CacheFile``
     - none
     - A netCDF file holding :math:`\psi_c` at this run's own points, read if
       it is there and matches and written if it is not. Refused under
       ``"meshed"``, which subtracts nothing. See below.
   * - ``InterpolatedFlux``
     - ``false``
     - Root the critical points on :math:`q_c`'s interpolant rather than on
       :math:`q_c`. A large speedup on a **warm** run and not safe on a cold
       one. See below.

``"meshed"``
   Each rectangle carries a uniform current density and enters as a **domain
   source** on the mesh, which is the route described above. The mesh must
   resolve the conductor, and ``[mesh.generate] CoilSize`` is what grades it.
   This is the only route that puts the conductor's own current inside the
   interior equation, which is where a force calculation would have to read it.

``"subtracted"``
   The same rectangles, integrated analytically and **taken out of the mesh**.
   MEQ writes :math:`\psi = \psi_c + \psi_p`, computes :math:`\psi_c` from the
   conductors exactly, and solves

   .. math::

      \Delta^{*}\psi_p = -\mu_0\, R\, J_{\mathrm{plasma}}(\psi_c + \psi_p),

   which is legitimate because :math:`\Delta^{*}` is linear and the conductor
   currents are prescribed inputs of a forward solve: the conductor's own
   :math:`\delta`-like source cancels exactly, so the mesh need not carry it.
   The conductors are then resolved to machine precision rather than to the
   mesh's order, and the mesh loses both the conductor elements and the grading
   around them.

``"filament"``
   A **point filament** at each rectangle's centre carrying the same total
   current, subtracted the same way. This is the only way MEQ can carry a
   filament at all: a point source has no finite-element representation as a
   current density, and :math:`\psi` near it is logarithmic, so the subtraction
   *is* the representation. It is also the cheaper field — one closed-form
   evaluation per point per conductor against a tensor quadrature per rectangle.

.. important::

   **The three are not three ways of writing one machine.** ``"meshed"`` and
   ``"subtracted"`` describe the same conductors and must agree to the
   discretisation; ``"filament"`` describes **different** ones. Measured on
   DIII-D's eighteen conductors, a filament model differs from the rectangles by
   :math:`6.1\times10^{-3}` globally and :math:`4.1\times10^{-2}` *inside* the
   conductors, against :math:`6.7\times10^{-5}` off them for the matched model
   — a per-cent-level answer near the coils and far better away from them. That
   is a mode to choose deliberately, and the ``.nc`` file records which one
   produced it in its ``conductor_model`` attribute.

   The difference is concentrated **at** the conductor, where :math:`\psi` of a
   line current diverges. Measured through a whole MEQ solve on
   ``examples/coils-rectangle.toml``, filament against rectangle: 1.23 relative
   everywhere, :math:`3.8\times10^{-3}` beyond 0.10 m of a coil centre — i.e.
   outside the conductor — and :math:`6.9\times10^{-4}` beyond 0.40 m. Choose
   ``"filament"`` when the field near the conductors is not what you are asking
   about.

.. _caching-psi-c:

Reusing :math:`\psi_c`: ``CacheFile``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The subtracted conductors' flux is built once per mesh, at every potential
degree of freedom and every source quadrature point. That is a precompute rather
than a per-iteration cost, but it is paid again by every *run* — and a transport
code driving MEQ changes the profiles while leaving the machine alone, so each
of those runs rebuilds a :math:`\psi_c` that could not have changed.

``CacheFile`` names a netCDF file to keep it in:

.. code-block:: toml

   [conductors]
   Model = "subtracted"
   CacheFile = "psi_c.nc"

The first run writes it and later runs read it. The file holds two things:
:math:`\psi_c` at every potential degree of freedom and every source quadrature
point, and :math:`q_c = (1/R)\,\bar\nabla\psi_c` at every flux-space element
node — the second being what the critical-point search needs, and the more
expensive of the two.

**What that is worth depends entirely on the conductor model**: for rectangles,
where one field point costs a whole cross-section quadrature, it is a measured
**8.4×** on the wall clock; for filaments, where a field point is one elliptic
integral, **1.05×** — there is simply nothing to save.

**The answer does not move.** Both fields are stored at exactly the points MEQ
evaluates them at, not sampled onto a grid and interpolated back, so a reload
*is* the recompute — every printed number is identical between a cold run and a
warm one.

**One thing is not cached**, and cannot be: writing the ``.nc`` grid adds
:math:`\psi_c` and :math:`B_c` back at every grid node, and the grid is not the
mesh. At ``GridNR = GridNZ = 129`` that is about a fifth of a warm run. Lower
the grid if you are running many warm solves and do not need it.

**A cache can only cost you time, never correctness.** It is bound to the mesh,
the polynomial degree and the conductors it was built for, and any of those
moving means the file is refused and rebuilt, with the reason printed. A file
that is absent, stale or corrupt is likewise a rebuild and a message. If you
change a coil current and forget to delete the file, MEQ notices.

The same data is written into the equilibrium ``.nc`` as a ``psi_coil`` group
whenever the split is in use, so a file that carries the answer also carries the
means to warm start from it. See :doc:`output`.

.. _interpolated-flux:

Searching faster: ``InterpolatedFlux``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finding the magnetic axis is an element-local Newton, and it evaluates
:math:`q_c` at every **iterate**. For subtracted rectangles each of those is a
cross-section quadrature of elliptic integrals, which makes the axis search the
largest item in a warm run — larger than the solve it serves.

MEQ already tabulates :math:`q_c` at the flux space's element nodes, and that
space is nodal, so the table is exactly a finite-element field. Setting

.. code-block:: toml

   [conductors]
   Model = "subtracted"
   InterpolatedFlux = true

roots that field instead: one polynomial evaluation per iterate. On a warm run
of a rectangle machine it is worth about **1.5×** on the wall clock, and takes
the axis leg from 9.5 s to 0.4 s.

.. warning::

   **This can change which equilibrium you get, and it is off by default for
   that reason.** The interpolant differs from :math:`q_c` by the interpolation
   error between nodes, which perturbs *which* candidate roots the search finds
   — and a normalised solve selects among discrete equilibria. Measured on a
   cold solve it moved :math:`\psi_{ax}` by up to 7.5e-02 where the exact field
   is reproducible to 5.6e-17, and the discrepancy does **not** fall under mesh
   refinement.

   **Use it warm.** Starting from a guess near the answer there is no branch
   left to select: the search refines one root rather than choosing among
   several, which is the case for a transport code re-solving a machine whose
   coils have not moved. In that regime it is measured to leave every printed
   digit of :math:`\psi_{ax}` unchanged.

.. _dividing-a-rectangle:

Dividing a rectangle: ``FilamentSize``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One point at a rectangle's centre is exact in the **far** field — the collapse
conserves the total current — and it is not a small approximation anywhere
else. On MAST-U's solenoid, 0.012 m by 3.180 m sitting beside the plasma rather
than beyond it, one filament at its centre gives a flux at the plasma that is a
**factor of 5.3** from the rectangle's own.

``FilamentSize`` divides each block instead. It is a target **cell** size in
metres, so a block of half-extents :math:`(a, b)` becomes a
:math:`\lceil 2a/s \rceil \times \lceil 2b/s \rceil` stack of filaments at the
cell centres, each carrying its share of the current:

.. code-block:: toml

   [conductors]
   Model = "filament"
   FilamentSize = 0.05      # metres

.. list-table::
   :header-rows: 1
   :widths: 22 26 18 34

   * - MAST-U block
     - cross-section
     - aspect
     - stack at 0.05 m
   * - ``Solenoid``
     - 0.012 × 3.180 m
     - 265
     - 1 × 64
   * - ``PX1``
     - 0.025 × 0.403 m
     - 16
     - 1 × 9
   * - ``D11``
     - 0.086 × 0.086 m
     - 1
     - 2 × 2

A **length** rather than a count, because the accuracy criterion is one: a
cell's error at a field point goes as :math:`(h/d)^2` in its standoff
:math:`d`, so one number buys a predictable accuracy across conductors of every
shape — and no single *count* is right for a machine whose blocks range over an
aspect ratio of 265.

A block may also name its own stack, which overrides the size for that block
alone:

.. code-block:: toml

   [[coils]]
   Name = "Solenoid"
   FilamentsR = 2
   FilamentsZ = 128

Both counts or neither: naming one is a parse error rather than a block whose
stack is half explicit and half derived. ``1`` is a legal count and is the
single filament in that direction.

.. important::

   **A divided block is a different conductor model from an undivided one, not
   a finer one, and the default is undivided for a reason.** The stack is the
   midpoint rule for the integral the rectangle's own field is, so it converges
   to the **rectangle** at second order in the cell size — measured, order
   2.000 over four refinements, taking an aspect-8 coil from
   :math:`7.4\times10^{-2}` at one filament to :math:`1.5\times10^{-6}` at
   16 × 128.

   That is towards the rectangle and therefore **away from a point-filament
   reference.** MEQ reproduces ``freegs4e``'s ``H_limited_circular`` — which
   carries point filaments at exactly the centres MEQ puts its own — to
   :math:`3.5\times10^{-5}` in :math:`\psi_{\mathrm{ax}}` undivided, and to
   :math:`1.2\times10^{-4}` with ``FilamentSize = 0.02``: **3.5 times worse**,
   and 10.8 times worse in :math:`\psi_{\mathrm{bnd}}`. Neither setting is the
   safe one. Divide when you are modelling the winding; leave it alone when you
   are comparing against a code that models filaments. A run says which it took,
   on stdout and in the ``.nc`` conductor table, which carries every filament
   individually.

   Refinement does **not** help *inside* the winding, and no cell size closes
   it: a stack has :math:`n_R n_Z` logarithmic singularities where a rectangle
   has none. Use ``"subtracted"`` if the field inside the conductor is what you
   are asking about.

.. warning::

   **Under a subtracting model the warm start changes meaning, and so does the**
   ``.gf`` **MEQ writes.** The solved field is :math:`\psi_p`, so
   ``[initialguess] File`` is read as :math:`\psi_p` and ``<stem>_psi.gf`` holds
   :math:`\psi_p`. A ``.gf`` carries no record of which it is — the format has
   no slot for one — so a guess written by a run *without* ``[conductors]`` is
   read as though :math:`\psi_c` had already been taken out of it. MEQ warns on
   standard output whenever both are set.

   **Both halves of that now have an answer.** ``[initialguess] Content =
   "total"`` says the stored file holds :math:`\psi` and MEQ takes
   :math:`\psi_c` off it as it reads it; and a subtracting run writes a fourth
   grid function, ``<stem>_psi_total.gf``, holding
   :math:`\psi_h + I_h( \psi_c )` for looking at. That one is exact at the
   nodes and an interpolation between them, which is the right trade for a
   picture and the wrong one for a restart — **the restart is the** ``.nc``,
   which carries the coefficients, the space, ``content``, and the conductor
   table that says what the remainder is a remainder *from*.

   **The two sampled outputs are unaffected and always carry the physical**
   :math:`\psi`. The gridded ``.nc`` adds :math:`\psi_c` back at every located
   node, and :math:`B` with it; ``<stem>_surfaces.nc`` traces level sets of
   :math:`\psi_c + \psi_p` and averages over them. Both can, because they
   *evaluate* :math:`\psi_c` at a point rather than representing it in the
   finite-element space — which is the same reason the ``.gf`` cannot, and for
   a filament could not even in principle.

.. note::

   **The table is refused when there is nothing to subtract.** ``Model`` naming
   anything but ``"meshed"`` with no ``[[coils]]`` blocks is a parse error, not a
   no-op: a run that asked to take its conductors out of the mesh and had none
   would converge, report every diagnostic it always did, and be slower than it
   should be, with nothing anywhere to look at. ``QuadratureOrder`` is refused
   on the other two models for the same reason — a filament has no cross-section
   to integrate over and a meshed coil is integrated by the mesh — and
   ``FilamentSize`` and ``[[coils]] FilamentsR``/``FilamentsZ`` are refused on
   the other two, which divide nothing. And
   ``[mesh.generate] CoilSize`` is refused beside a subtracting model: it grades
   the mesh around conductors the generator is no longer told about, so it names
   an element size for a region that does not exist.

.. note::

   **A subtracting model changes the mesh, and that is the point.** With
   ``[mesh.generate]``, ``meq --mesh-command`` emits no ``--coil`` at all under
   ``"subtracted"`` or ``"filament"``, so the generator neither fragments nor
   grades the conductors. ``meq-run`` stamps the mesh with the command it used,
   so changing ``Model`` changes the command and the mesh is rebuilt with
   nothing else to do. Measured: ``examples/mastu-nke``'s 23 conductors take its
   mesh from 9361 triangles to 4716, and
   ``examples/diverted-tokamak-generated``'s four from 2854 to 2322.

.. note::

   **This example is not free boundary.** With ``[boundary] Type = "zero"`` on a
   mesh that is the plasma, the conductors sit *inside* the plasma and the
   boundary is where the file says it is rather than where the coil currents put
   it. A machine case additionally needs a domain with a vacuum region, the
   plasma boundary flux as an unknown, and the exterior coupling that makes
   :math:`\psi` on the outer boundary the field the currents produce. The last
   two are ``[boundary.limiter]`` and ``[boundary.exterior]`` above.

``[boundary.limiter]``
~~~~~~~~~~~~~~~~~~~~~~

What pins :math:`\psi_{\mathrm{bnd}}`, making it an unknown of the bordered
Newton beside :math:`\psi_{\mathrm{ax}}`.  Give **either** a contact point or a
meshed limiter surface to find the contact on — not both.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``R``, ``Z``
     - *both required if either is given*
     - The limiter contact, in metres, prescribed. ``R`` must be strictly
       positive.
   * - ``SurfaceAttribute``
     - *none*
     - The **element** attribute of the region the limiter encloses.  The
       limiter is that region's boundary and the contact is **found** on it.

The profiles are functions of
:math:`\Psi = (\psi - \psi_{\mathrm{bnd}})/(\psi_{\mathrm{ax}} - \psi_{\mathrm{bnd}})`,
so :math:`\psi_{\mathrm{bnd}}` is a functional of the solution exactly as
:math:`\psi_{\mathrm{ax}}` is.

With ``R`` and ``Z`` the constraint is that it equals :math:`\psi_h` **at the
point given**, evaluated inside the element containing it.  With
``SurfaceAttribute`` it is
:math:`\psi_{\mathrm{bnd}} = \max \psi_h` **over the limiter surface**, which is
what a machine actually does: the plasma edge is the flux surface that touches
the limiter, and where it touches is an output of the solve.  Either way the
border row is the contact element's potential shape functions there, exact and
undifferenced — for the found contact because at a maximum along the surface the
*tangential* derivative vanishes while the contact moves tangentially, so the
envelope theorem removes the position term, as it does for
:math:`\psi_{\mathrm{ax}}` at the located axis.

.. note::

   **The mesh must be fitted to the limiter**, which is what
   ``SurfaceAttribute`` means: the limiter is the boundary of a region of
   elements, so it is a union of mesh faces and the contact lies on one.
   :program:`tools/mesh/halfdisc.py --limiter` fragments the limiter circle into
   the geometry and writes the enclosed region as attribute ``20``.  The limiter
   MEQ then uses is the **polygon** those faces make; fitting the mesh to the
   curve puts its vertices on the true curve, so the polygon inscribes rather
   than approximates.

.. note::

   Prescribing the contact is an :math:`O(h)` choice and finding it is not, so
   prefer ``SurfaceAttribute`` when the limiter can be meshed.  A contact
   named in the wrong place is wrong by
   :math:`\mathrm{dist} \times |\nabla\psi|` however fine the mesh — measured,
   a contact 0.15 m along the same curve is nearly thirty thousand times further
   from the answer than the found one, and its error does not fall under
   refinement at all.  If you must prescribe a point, give it the point the
   machine actually touches, and be aware that a contact read off another
   code's grid is a maximum over *cells*.

   Where the contact was actually found is reported on the run and written into
   the ``.nc`` as ``limiter_r``, ``limiter_z`` and ``limiter_contact_located``.

.. warning::

   A limiter is **refused without** ``[source] Normalised = true``.
   :math:`\psi_{\mathrm{bnd}}` enters only through :math:`\Psi`, so on a source
   that does not read it the border would be solved and its answer discarded —
   and the run would converge, at full order, to the equilibrium the file did
   not describe.

``[boundary.xpoint]``
~~~~~~~~~~~~~~~~~~~~~

The same unknown as ``[boundary.limiter]`` — :math:`\psi_{\mathrm{bnd}}` — for a
**diverted** plasma, where the bounding point is an X-point rather than a piece
of hardware.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``R``, ``Z``
     - *both required if either is given*
     - Where the X-point is **believed** to be, in metres.  A seed for an
       unknown, not a prescription.  ``R`` must be strictly positive.

A limiter contact is where the drawings say it is.  An X-point is not: it is a
functional of the solution and moves as Newton moves.  So this block does not
pin :math:`\psi_{\mathrm{bnd}}` at a point — it makes the point itself two more
unknowns of the same bordered Newton, closing

.. math::

   q_r(R_X, z_X) = 0, \qquad
   q_z(R_X, z_X) = 0, \qquad
   \psi_{\mathrm{bnd}} - \psi_h(R_X, z_X) = 0

on the same factorisation as everything else.  The two new unknowns cost **no
extra backsolve**: the field residual does not contain :math:`(R_X, z_X)` at
all — they reach it only through :math:`\psi_{\mathrm{bnd}}`, which has a column
already — so both new columns are exactly zero and the system grows only in its
dense corner.

.. note::

   **The seed is an answer's starting value, not the answer.**  A file whose
   numbers are a few centimetres out describes the same equilibrium as one whose
   numbers are exact, and the run reports where the solve actually put the
   null — on the terminal, and in the ``.nc`` as ``xpoint_r``, ``xpoint_z`` and
   ``xpoint_located``.  ``examples/diverted-tokamak-xpoint.toml`` deliberately
   seeds 7 cm away and reports the distance it travelled.

   What the seed **does** decide is **which** null.  This follows one saddle,
   exactly as ``[source] PsiAxis`` follows one O-point, and a machine with two
   nulls has two answers: the shipped example's is an up-down asymmetric double
   null whose upper saddle sits 1.4 m away and 3.5e-03 further out in
   :math:`\psi`.  Get the seed into the right half of the machine.

.. warning::

   ``[boundary.limiter]`` and ``[boundary.xpoint]`` are **alternatives** and
   naming both is refused: all three routes to :math:`\psi_{\mathrm{bnd}}` pin
   one unknown, and a precedence rule would decide which equilibrium the run
   reports on the strength of key order.  Like a limiter, an X-point is refused
   without ``[source] Normalised = true``.

.. note::

   **A diverted run wants** ``[solver] PlasmaSupportSweeps`` **as well.**  The
   border makes the X-point's *position* an unknown, which it can be because it
   is differentiable.  *Which elements carry current* is not, and on a diverted
   plasma the level set is disconnected across the null, so that choice is a
   real one.  See that key, and ``[source] ConfineToPlasma``.

``examples/diverted-tokamak-xpoint.toml`` is the worked example, and
``examples/diverted-tokamak.toml`` is the same machine with the contact
prescribed instead — the comparison between the two is what
``tests/convergence/XPointBorder.cpp`` measures.

``[boundary.exterior]``
~~~~~~~~~~~~~~~~~~~~~~~

Present, this makes the run **free boundary**: :math:`\Gamma` is an artificial
boundary in the vacuum carrying no prescribed datum, and the exact exterior
Dirichlet-to-Neumann map stands in for everything outside it.

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Radius``
     - *required*
     - :math:`\rho_\Gamma`, the radius of the semicircle, in metres. Must fit
       strictly inside the ``[mesh]`` box.
   * - ``CentreZ``
     - ``0.0``
     - The axial position of its centre.
   * - ``Modes``
     - ``4``
     - How many Gegenbauer modes, degrees :math:`2 \ldots \texttt{Modes}+1`.

The map is **diagonal**, which is why this block is three numbers rather than a
boundary-element solver: separating :math:`\Delta^*` in spherical coordinates
gives the Gegenbauer equation, whose modes
:math:`\rho^{1-n} C_n(\cos\theta)` each decay independently. The truncation at
``Modes`` is the only approximation in the exterior — the map is exact mode by
mode — and its error falls spectrally in the smoothness of the trace on
:math:`\Gamma`. The converged coefficients are printed by the driver and written
into the ``.nc`` as ``exterior_a2``, ``exterior_a3``, …; they *are* the exterior
solution, since :math:`\psi` outside the mesh is their sum against the basis.

.. warning::

   **The mesh must reach the axis exactly**: ``[mesh] RMin = 0``, and the driver
   refuses anything else. The separation above holds on a semicircle centred on
   the axis and nowhere else, so a box starting at :math:`R = 0.05` does not
   give a slightly worse version of the problem — the modes do not span its
   exterior at all. Every mode vanishes on the axis identically, so the flat
   side needs no treatment in the exterior and is ordinary fitted boundary.

   ``[boundary.exterior]`` and ``[boundary.shape]`` are **alternatives**, and
   naming both is refused: one makes :math:`\Gamma` a semicircle about the axis
   and the other a closed surface that may not reach :math:`R = 0`.
   ``[boundary] Type`` must be ``"zero"`` beside an exterior block, since
   :math:`\Gamma` carries the transmission condition rather than data.

.. note::

   **A limiter and an exterior coupling together are the machine case**, and
   they converge: :math:`\psi_{\mathrm{bnd}}` is one more border and the
   Gegenbauer coefficients are :math:`N` more, all on one factorisation per
   Newton step.  ``examples/limited-tokamak.toml`` is that case and reaches
   freegs4e's own answer to 1.3e-04.  With the limiter absent
   :math:`\psi_{\mathrm{bnd}}` stays at zero, so the plasma edge is the
   :math:`\psi = 0` contour rather than a contact.

``examples/free-boundary-halfdisc.toml`` is the worked example without a
limiter; ``examples/limited-tokamak.toml`` and
``examples/diverted-tokamak-xpoint.toml`` are the two machine cases.

``[output]``
------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Directory``
     - ``"."``
     - **Not created by MEQ.** A missing directory is exit code 3.
   * - ``Prefix``
     - ``"meq"``
     - The stem of every output name; must not be empty.
   * - ``GridNR``, ``GridNZ``
     - ``129``
     - Sampling **nodes** for the gridded output, so the spacing is
       :math:`(R_{\max} - R_{\min})/(\texttt{GridNR} - 1)`. Nothing to do with
       ``[mesh] NR``.
   * - ``FluxSurfaces``
     - ``false``
     - Write ``<stem>_surfaces.nc``, the flux-surface geometry and the
       flux-surface averages against a flux label. Off by default: it costs a
       contour trace and an angle fit per surface, and it is the one output that
       can be impossible on a run that solved perfectly well.
   * - ``FluxSurfaceCount``
     - ``24``
     - Surfaces in the family. At least 2.
   * - ``FluxAngleCount``
     - ``128``
     - Nodes on each surface, equispaced in the poloidal angle about the
       magnetic axis. At least 3.
   * - ``FluxInnerCut``, ``FluxOuterCut``
     - ``0.05``, ``0.95``
     - The range of normalised flux the family covers, strictly inside
       :math:`(0, 1)`. Both ends are cut and for different reasons; see
       :ref:`output-flux-surfaces`. A value at or past either end is **refused**
       rather than clamped, whether or not the file is being written.

See :doc:`output` for what gets written.

``[initialguess]``
------------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Type``
     - ``"none"``
     - ``"none"`` starts from the Dirichlet datum — a cold start. ``"ramp"``
       makes :math:`\psi` run from :math:`-\texttt{Amplitude}` to
       :math:`+\texttt{Amplitude}` across :math:`z`. ``"conductors"`` builds a
       guess out of this file's own ``[[coils]]`` — see below. ``"bump"`` is a
       paraboloid core. ``"gridfunction"`` reads a stored answer.
   * - ``Amplitude``
     - ``0.3``
     - For ``"ramp"`` and ``"bump"``; must be positive. For ``"bump"`` it is the
       peak value at the centre.
   * - ``CentreR``, ``CentreZ``
     - *``CentreR`` required for* ``"bump"``
     - The centre of the paraboloid, in metres. For ``"conductors"`` it is the
       guessed magnetic axis and is optional; leaving it out gives the
       conductors' field with no plasma column in it.
   * - ``RadiusR``, ``RadiusZ``
     - *``RadiusR`` required for* ``"bump"``
     - Its extent. ``RadiusZ`` defaults to ``RadiusR``. ``CentreR - RadiusR``
       must be strictly positive: a bump reaching the axis describes no plasma.
       For ``"conductors"`` these are the semi-axes of the plasma column and
       ``RadiusR`` defaults to half of ``CentreR``, which cannot reach the axis
       whatever ``CentreR`` is.

   * - ``File``, ``MeshFile``
     - *required for* ``"gridfunction"``
     - The stored grid function and the mesh it lives on — a grid function
       cannot be read without its mesh.
   * - ``Content``
     - ``"remainder"``
     - What that file *holds*: ``"total"`` for :math:`\psi` or ``"remainder"``
       for :math:`\psi_p = \psi - \psi_c`. Read only with
       ``Type = "gridfunction"``, and it matters only when this run subtracts —
       see below.

``Content`` says what a stored guess holds, because the file cannot
--------------------------------------------------------------------

A ``.gf`` carries the finite-element space and the coefficients and nothing
else, so a file written under ``[conductors] Model`` holds
:math:`\psi_p = \psi - \psi_c` and looks exactly like one that holds
:math:`\psi`. ``Content`` is where the file's author says which it is.

``"remainder"``, the default, is what a subtracting run writes and what every
file written before this key existed means. ``"total"`` says the file holds the
physical flux — a *meshed* run's answer, say — and MEQ takes :math:`\psi_c` off
it at every node as it reads it, reporting the largest shift on standard output.

**Without it one direction is unreachable, which is why this is not a second
spelling of the same thing.** Restarting a subtracted run from its own output
wants ``"remainder"``. Restarting from a meshed answer — the whole of the
coarse-then-fine route, where a cheap filament solve hands a meshed one a
starting point — wants ``"total"``, and there was no way to ask: the guess
arrived one entire conductor field away from where the solver believed it was.

The key describes the *file*, which does not change because the run reading it
did, so it is accepted and ignored on a run that does not subtract rather than
refused. :math:`\psi_c` is then zero and the two spellings mean the same file.

``Type = "conductors"`` builds the guess from what the file already says
-------------------------------------------------------------------------

A free-boundary problem has more than one converged solution, so something has
to say which is wanted — but that is a physical statement, not a flux field, and
MEQ has everything it needs to build one.

The guess is the sum of two things. The **conductors' own vacuum field**, over
the ``[[coils]]`` in this same file, computed by Carlson's elliptic integrals
and :math:`\Delta^*`-harmonic off the conductors at measured rate 2.00. And,
when ``CentreR`` is given, **the plasma's own current** — ``[source]
PlasmaCurrent`` spread uniformly over an ellipse about that guessed magnetic
axis.

The plasma half is what selects the branch. The conductors alone have the right
scale and topology outside the plasma and say nothing about the core; on
freegsnke's MAST-U they converge in nine iterations to a magnetic axis outside
the machine. The column is an **ellipse** rather than a point or a rectangle for
two measured reasons: the flux of a filament diverges logarithmically at the
filament, which puts an unbounded spike exactly where the axis search has to
look, and a uniform current density over a rectangle carries a logarithm in its
second derivatives at each of its four corners.

How well the ellipse matches the real plasma does not matter much. On MAST-U a
column tuned to the machine and a round one at the default ``0.5*CentreR`` reach
the same :math:`\psi_\mathrm{ax}`, the same :math:`\psi_\mathrm{bnd}` and the
same X-point to every printed digit, in the same number of Newton iterations.
``"bump"`` gives
:math:`\psi = \texttt{Amplitude}\,(1 - (\Delta R/\texttt{RadiusR})^2 - (\Delta z/\texttt{RadiusZ})^2)`
where that is positive and zero elsewhere.

.. important::

   **Once the boundary is free the guess is part of the problem statement.** A
   ramp is antisymmetric in :math:`z` and describes no plasma; it exists so that
   :math:`\psi = 0` falls in the interior and the trivial branch is not a fixed
   point. A free-boundary problem has **more than one converged solution**, and
   the guess is what says which of them is reported — so a core is what selects
   the branch that is an equilibrium. See :doc:`validation`, where three solve
   routes reach discrete solutions 9.4 % apart and the physical root lies
   *between* two that an amplitude sweep reaches.

The ramp is not a nicety: see :ref:`sources-trivial-branch`. A stored guess on
the *same* mesh is an exact restart; on a different one it is interpolated. It
is applied on the first cycle only — later adaptive cycles warm-start from the
previous cycle.

``[adaptivity]``
----------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``Enabled``
     - ``false``
     - 
   * - ``MaxIterations``
     - ``10``
     - Counts **solves**, not refinements: 4 means at most 3 refinements.
   * - ``Strategy``
     - ``"doerfler"``
     - ``"doerfler"`` or ``"maximum"``.
   * - ``Theta``
     - ``0.6``
     - In :math:`(0, 1]`. **The two strategies read it oppositely** — see
       :doc:`adaptivity`.
   * - ``TargetError``
     - ``1.0e-6``
     - Absolute, in the estimator's own norm.

.. _configuration-q-driven:

Driving the equilibrium by :math:`q(\Psi)`
------------------------------------------

Every other configuration in :doc:`examples` **prescribes** the toroidal field
and reports the safety factor as an output. A transport code hands an
equilibrium code the other way round: :math:`q(\Psi)` is the target and
:math:`g(\Psi)` is what has to be found. ``[source] SafetyFactorFile`` is that
run, and ``examples/q-driven.toml`` is the worked example.

The algebra is a division. RoPP (142) gives

.. math::

   q = \frac{V' \, g \, \langle R^{-2}\rangle}{4\pi^2},
   \qquad\text{so}\qquad
   g = \frac{4\pi^2 q}{V' \, \langle R^{-2}\rangle},

one division per flux surface at **fixed geometry**. What makes it a solver is
that :math:`V'` and :math:`\langle R^{-2}\rangle` are functionals of the
solution, and the solution depends on :math:`g`. So the loop is

.. code-block:: text

   solve  ->  extract the surfaces  ->  invert  ->  rebuild gg'  ->  solve

and its fixed point is an equilibrium whose own :math:`q` is the one asked for.

What it costs
~~~~~~~~~~~~~

**Every map evaluation is a whole equilibrium** — one per residual, and
:math:`2(d+1)` more per Jacobian column. The shipped example converges in 12
outer iterations and **43 equilibria**. It is affordable only because the
profile is *fitted*: the unknown is the handful of coefficients of
:math:`g^2 = \sum_j c_j \Psi^j` rather than every degree of freedom of the
field.

One solver serves all of them. Only :math:`g\,\mathrm{d}g/\mathrm{d}\Psi`
changes between steps, so the mesh, the finite element spaces, the assembled
forms and the trace solver's symbolic factorisation all survive, and each solve
is warm-started from the previous one's field.

A relaxed Picard iteration cannot replace the outer Newton, and that is a
theorem rather than a preference: :math:`c \leftarrow c + \omega(G(c) - c)` has
derivative :math:`1 + \omega(G' - 1)` at the fixed point, which exceeds one for
**every** :math:`\omega > 0` when :math:`G' > 1`. Under-relaxation stabilises a
map that oscillates and does nothing at all for one that runs away.

Two things to get right
~~~~~~~~~~~~~~~~~~~~~~~

**The table is in the source's** :math:`\Psi`, one on the axis — like every
other profile table in ``examples/``, and unlike the flux-surface family's own
label, which is zero there. MEQ owns the reflection between them. A table
written the other way round does **not** fail: it converges, at full order, to
an equilibrium with its **shear reversed**, which is a configuration a real
machine can have and which nothing downstream will look at twice.

**The degree is a modelling choice and higher is not better.** A degree the
surface family does not determine leaves the outer Jacobian rank deficient,
which MEQ reports — ``did not converge … the outer Jacobian is RANK
DEFICIENT`` — rather than solving. Two or three is the usual.

What the run writes
~~~~~~~~~~~~~~~~~~~

The answer is :math:`g`, so it goes into the interchange file: ``<stem>.nc``
carries ``toroidal_field_driven``, ``g_squared_coefficients`` — ascending in
:math:`\Psi` — and ``safety_factor_target``. There is no ``GGPrimeFile`` beside
the output for a reader to look the field up in, which is why it is recorded
rather than left to the configuration.

**And with** ``[output] FluxSurfaces`` **this is the one route whose surfaces
file reports** :math:`q` **itself.** Everywhere else a ``meq::Source`` carries
:math:`g g'` and not :math:`g`, so there is no :math:`g` to build
:math:`V' g \langle R^{-2}\rangle / 4\pi^2` from and the ``safety_factor``
column is absent rather than zero — see :doc:`output`. Here the loop solved for
:math:`g`, so the column is written, and the round trip can be read straight off
the file: the reported :math:`q` against the table the run was given. On
``examples/q-driven.toml`` the two agree to **4.5e-07** over all 24 surfaces.

``BorderRegularisation``, ``BorderCollinearityRegularisation``
   Levenberg damping on the dense border solve, so a near-singular Schur
   complement gives a damped step instead of failing. Both default to ``0.0``,
   which is off and leaves the solve bit-identical. The damping is scaled by the
   relative merit, so it vanishes as the solve converges and cannot cost the
   quadratic endgame. ``BorderCollinearityRegularisation`` weights a per-row
   penalty on border rows that have gone parallel to one another, and does
   nothing unless ``BorderRegularisation`` is also set.

``TopologyRetry``
   How many gentle (0.75) reductions a trial step that could not be evaluated at
   all -- a normalisation the source refuses, an X-point outside the mesh -- gets
   before the ordinary halving line search takes over. Default ``0``, which is
   off. A step that breaks the plasma's topology and a step that merely made the
   residual worse are different failures; this is what lets the first be treated
   as one.

``UpDownSymmetry``
   Project every Newton iterate onto the subspace of fields even in :math:`z`.
   Default ``false``. A double null has two saddles at identical flux and
   ``[boundary.xpoint]`` follows whichever its seed is nearer; with the
   constraint on, the two are exactly degenerate and the choice stops mattering.

   **It needs a mirror-symmetric mesh and it refuses one that is not.** The
   projection averages each degree of freedom with the one at its own
   reflection, so a degree of freedom with no partner has nothing to average
   against; MEQ names the first such element rather than projecting away an
   asymmetry the mesh describes. That is a real restriction rather than a
   formality: an unstructured mesher given a perfectly symmetric machine does
   not generally return a symmetric mesh, and neither does a triangulated
   Cartesian grid, which splits every cell along one diagonal. A quadrilateral
   grid symmetric about :math:`z = 0` does.

   **Where MEQ makes the mesh, say so there too.** Set
   ``[mesh.generate] Symmetric = true``, which meshes :math:`z \ge 0` and
   reflects it. Asking for ``UpDownSymmetry`` on a generated mesh *without*
   that key is a **parse error**, because the mirror maps are built inside the
   nonlinear solve rather than at setup — the refusal would otherwise arrive
   minutes in, after the mesh, the spaces, the assembly and the initial guess.
   A mesh you made yourself is left alone: MEQ cannot know whether it is
   symmetric without looking, and finds out the slow way.
