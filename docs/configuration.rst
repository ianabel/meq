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
``[adaptivity]`` — are optional, and an absent table behaves exactly like a
present but empty one. ``[[coils]]`` is an *array* of tables rather than a
table, so "absent" there means no conductors at all.

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
     - Lower radial bound. Must not be negative — :math:`r` is a cylindrical
       radius, and a box reaching :math:`r = 0` contains the operator's
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
     - :math:`F = -\left((1-A)r^2 + A\right)`, independent of :math:`\psi`. The
       problem is linear and Newton converges in one step.
   * - ``"mhd"``
     - :math:`F = \mu_0 r^2 p'(\psi) + (gg')(\psi)`, from two tabulated
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
     - Vacuum permeability on the :math:`r^2 p'` term. Set to 1 for a run in
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
     - The radial **offset** in :math:`\sin(K_r(r + R_0))`. **Not a major
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
       :math:`r = \texttt{ReferenceRadius}`.
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
   is assembled as :math:`\int F/r`, which *is* :math:`\mu_0 I_p` — and a solver
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

   F_{\mathrm{coil}}(r, z) = \mu_0\, r\, \frac{I_k}{|\Omega_{c_k}|}

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
     - Half-extents in :math:`r` and :math:`z`, metres. Strictly positive, and
       ``CentreR - HalfWidth`` must be strictly positive too — a coil reaching
       the axis is refused, because the operator's :math:`1/r` is not integrable
       through :math:`r = 0`.
   * - ``Current``
     - *one of the two*
     - The **total** current through the cross-section, amperes. Signed, and
       zero is allowed. A real winding's turns are not modelled, so this is
       turns × amps per turn.
   * - ``CurrentDensity``
     - *one of the two*
     - The uniform :math:`j_\phi` in A/m², multiplied by the area
       :math:`4\,\texttt{HalfWidth}\,\texttt{HalfHeight}`.

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
   the axis and nowhere else, so a box starting at :math:`r = 0.05` does not
   give a slightly worse version of the problem — the modes do not span its
   exterior at all. Every mode vanishes on the axis identically, so the flat
   side needs no treatment in the exterior and is ordinary fitted boundary.

   ``[boundary.exterior]`` and ``[boundary.shape]`` are **alternatives**, and
   naming both is refused: one makes :math:`\Gamma` a semicircle about the axis
   and the other a closed surface that may not reach :math:`r = 0`.
   ``[boundary] Type`` must be ``"zero"`` beside an exterior block, since
   :math:`\Gamma` carries the transmission condition rather than data.

.. note::

   **A limiter and an exterior coupling together do not yet converge.** Each
   works on its own — :math:`\psi_{\mathrm{bnd}}` as a second border, and the
   Gegenbauer coefficients as :math:`N` more — and the combination is the one
   still open. With the limiter absent :math:`\psi_{\mathrm{bnd}}` stays at
   zero, so the plasma edge is the :math:`\psi = 0` contour rather than a
   limiter contact.

``examples/free-boundary-halfdisc.toml`` is the worked example.

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
       :math:`+\texttt{Amplitude}` across :math:`z`. ``"bump"`` is a paraboloid
       core. ``"gridfunction"`` reads a stored answer.
   * - ``Amplitude``
     - ``0.3``
     - For ``"ramp"`` and ``"bump"``; must be positive. For ``"bump"`` it is the
       peak value at the centre.
   * - ``CentreR``, ``CentreZ``
     - *``CentreR`` required for* ``"bump"``
     - The centre of the paraboloid, in metres.
   * - ``RadiusR``, ``RadiusZ``
     - *``RadiusR`` required for* ``"bump"``
     - Its extent. ``RadiusZ`` defaults to ``RadiusR``. ``CentreR - RadiusR``
       must be strictly positive: a bump reaching the axis describes no plasma.
   * - ``File``, ``MeshFile``
     - *required for* ``"gridfunction"``
     - The stored grid function and the mesh it lives on — a grid function
       cannot be read without its mesh.

``"bump"`` gives
:math:`\psi = \texttt{Amplitude}\,(1 - (\Delta r/\texttt{RadiusR})^2 - (\Delta z/\texttt{RadiusZ})^2)`
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
