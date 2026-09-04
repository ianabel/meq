Code organisation
=================

.. code-block:: text

   src/meq/     the library
   apps/        drivers -- only meq.cpp
   tests/       unit/       Boost.Test, no MFEM needed
                convergence/ rate assertions
                analytic/   closed-form solutions used by both
                performance/ built, deliberately NOT a registered test
   tools/       plotting and visualisation
   examples/    TOML run configurations
   refs/        Refs.md is tracked; the PDFs are gitignored, fetch by DOI
   docs/        this documentation; docs/manual/ is the pre-Sphinx LaTeX manual
   attic/       free-boundary work, not ported and not built, kept visible

The library
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - Header
     - Needs MFEM
     - What it is
   * - ``Config``
     - no
     - The TOML schema and its parser. Throws on anything it does not
       understand.
   * - ``Profiles``
     - no
     - One-dimensional functions of flux, with two exact derivative levels.
   * - ``Source``, ``RotatingSource``
     - no
     - :math:`F` and :math:`\partial F/\partial\psi`.
   * - ``SourceFactory``
     - no
     - Configuration to source, including the refusal that keeps the normalised
       path separate.
   * - ``BoundaryShape``
     - no
     - Miller and MXH plasma boundaries, and their level sets.
   * - ``GradShafranov``
     - yes
     - The solver. Everything else in this list serves it.
   * - ``Estimator``
     - yes
     - The residual estimator, the marking strategies, and the companion mesh.
   * - ``Field``
     - yes
     - :math:`q \to \mathbf{B}`, a relabelling.
   * - ``Sampler``
     - yes
     - Finite element fields onto a uniform grid, including the band.
   * - ``WarmStart``
     - yes
     - Interpolating one mesh's answer onto another.
   * - ``Output``
     - yes
     - The three output formats.
   * - ``CriticalPoints``
     - yes
     - The magnetic axis and any X-point, as roots of the flux, plus the
       boundary index audit over them.
   * - ``FluxSurfaces``
     - yes
     - The contour tracer, the band extension, and the poloidal-angle
       parametrisation.
   * - ``SurfaceAverage``
     - yes
     - Flux-surface averages over a callable integrand.
   * - ``Zernike``
     - no
     - The disc basis, and the conversions between the flux label and the disc
       radius.
   * - ``SurfaceFit``
     - no
     - The surfaces as a map from a disc: the linear fit, and the gauge-free
       refit.

.. _organization-mfem-free:

Why half of it does not include MFEM
------------------------------------

``Config``, ``Profiles``, ``Source``, ``RotatingSource``, ``SourceFactory``,
``BoundaryShape``, ``Zernike`` and ``SurfaceFit`` take plain ``double``
arguments and know nothing about finite elements. The ``mfem::Coefficient``
adapters live with the assembly that needs them, and the two geometry headers
take their field as a callable rather than as a grid function — which is why a
caller writes the short loop that turns a traced surface into samples. See
:ref:`geometry-disc`.

Two things follow, and the second is what makes it worth the discipline:

* Both layers are **unit-testable without the library**, and they are tested
  that way, exhaustively.
* **Continuous integration can build them.** The MFEM branch meq needs is
  published nowhere a hosted runner can fetch it, so ``find_package`` is
  deliberately not ``REQUIRED`` and CI builds this half. See :doc:`testing` for
  what that does and does not establish.

The data flow of a run
----------------------

.. code-block:: text

   TOML file
     -> meq::Configuration                     validated, all-or-nothing
     -> meq::makeSource / makeNormalisedSource  -> meq::Source
     -> mesh                                    box or file, then refinement
     -> meq::BoundaryShape                      optional; the curved path
        -> subdomain + transfer paths           -> meq::AdaptiveDomain if adaptive
     -> meq::GradShafranovSolver
        setSource / setExtension / setBoundaryData / setInitialGuess
        solve()                                 -> potential, flux, trace
        postProcess()                           -> psi*, enriched flux
     -> meq::ResidualEstimator                  if adaptive: estimate, mark, refine
     -> meq::poloidalField                      q -> B
     -> meq::GridSampler + meq::NetCDFWriter    the gridded file
     -> meq::writeMfem                          the exact restart
     -> meq::curveBoundaryOnto + meq::writeVtu  the picture, LAST

The ordering at the end is not incidental: ``curveBoundaryOnto`` mutates the
mesh geometry, so everything that reads the solved geometry must run before it.

.. note::

   The flux-surface machinery is **not** in that flow. ``CriticalPoints``,
   ``FluxSurfaces``, ``SurfaceAverage``, ``Zernike`` and ``SurfaceFit`` are
   library-only: nothing in the TOML schema reaches them and nothing they
   produce is written to an output file. They are a second consumer of a solved
   ``GradShafranovSolver``, alongside the estimator and the sampler. See
   :doc:`flux_surfaces` and :doc:`surface_geometry`.

Relationship to MFEM
--------------------

meq is a thin layer over ``DarcyForm`` with hybridization enabled. What meq
supplies on top:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - meq
     - Because
   * - :cpp:class:`meq::ConstantStabilization`
     - The library's default stabilisation is the LDG choice, which costs an
       order in the flux. See :ref:`formulation-tau`.
   * - :cpp:class:`meq::SourceIntegrator`
     - :math:`F` and :math:`\partial F/\partial\psi` in the weak form's
       :math:`1/r` weighting.
   * - The offsets for a three-block vector
     - The library's own offsets stop at the potential; the trace is meq's to
       track.
   * - The bordered Newton
     - :math:`\psiax` as an unknown is not something a hybridized solver can
       express internally. See :doc:`normalised_flux`.
   * - The residual estimator
     - The library's built-in error estimator computes a different quantity and
       does not consult an installed stabilisation.

.. note::

   **The old HDG interface in MFEM is gone**, which is why meq's relationship to
   it is a port rather than a recompile. The hand-written bilinear form, its
   domain and face integrators, and its explicit condense-and-reconstruct calls
   have no counterparts; the replacement is ``DarcyForm`` with hybridization,
   which owns the block structure and the spaces meq used to wrap itself.

.. warning::

   **One integrator in the assembly is load-bearing for a reason unrelated to
   the integral it appears to compute.** Under hybridization the potential-flux
   coupling comes from the transpose of the trace-jump integrator supplied to
   the hybridization, *not* from the face integrators on the mixed block —
   which are never assembled at all, since the assembly sums domain integrators
   only. Changing their coefficients does not move a single digit.

   But the boundary-face one still has to be there, as a **marker**: enabling
   hybridization reads the boundary-face markers to decide where to register the
   flux constraint. Remove it and the Dirichlet faces get no constraint, and the
   error goes flat. This is exactly the kind of thing to leave a note about.

.. _organization-naming:

Naming conventions
------------------

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - Kind
     - Convention
   * - Types — class, struct, enum, alias
     - ``UpperCamelCase``
   * - Enum values
     - ``UpperCamelCase``
   * - Functions and methods
     - ``lowerCamelCase``
   * - Variables, parameters, members
     - ``lowerCamelCase``
   * - **TOML configuration keys**
     - **``UpperCamelCase``**

That last row is a deliberate mismatch with the C++ rule, inherited from a
sibling project's convention set. **Do not reconcile the two.**

The style is tab indentation, Allman braces, C++17.

.. note::

   **External names keep their author's capitalisation**, whatever the table
   says. ``TraceSolver::cuDSS`` is spelled the way its vendor spells it and
   carries a suppression comment saying so; the other backends happen to need
   none. The house rule governs meq's own identifiers and does not extend to
   renaming other people's products — the same exemption the linter
   configuration already records for interface methods MFEM imposes.

This is **enforced**, by a ``clang-tidy`` identifier-naming check registered as
the ``naming`` test over every source in the library. See :doc:`testing`.
