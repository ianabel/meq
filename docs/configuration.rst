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
present but empty one.

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

      meq: configuration error in 'run.toml', key 'mesh.RefinmentLevels':
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

.. warning::

   ``File`` **is effectively supported on the fitted path only.** The box bounds
   are what the gridded output samples over and what sets the search length for
   the curved boundary's transfer paths; with a mesh file and no
   ``[boundary.shape]`` to supply a bounding box, both are degenerate.

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
       calls ``FF'``.
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

.. note::

   Profile paths are resolved against the **working directory of the run**, not
   against the directory the configuration file lives in.

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
   * - ``Mu0``, ``Normalised``, ``PsiAxis``
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
for a facility that does not exist yet and are refused rather than ignored.

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
       :math:`+\texttt{Amplitude}` across :math:`z`. ``"gridfunction"`` reads a
       stored answer.
   * - ``Amplitude``
     - ``0.3``
     - For ``"ramp"``; must be positive.
   * - ``File``, ``MeshFile``
     - *required for* ``"gridfunction"``
     - The stored grid function and the mesh it lives on — a grid function
       cannot be read without its mesh.

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
