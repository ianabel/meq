C++ API reference
=================

MEQ is a library first and a program second. Everything the ``meq`` executable
does is done through the public interface documented here, and
``apps/meq.cpp`` is a worked example of using it.

.. code-block:: cpp

   #include "meq/meq.hpp"

   mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D( 8, 8, mfem::Element::TRIANGLE,
                                                  false, 0.8, 1.2 );
   meq::SolovievSource source( -0.52 );
   mfem::ConstantCoefficient zero( 0.0 );

   meq::GradShafranovSolver solver( mesh, /*order=*/3 );
   solver.setSource( source );
   solver.setBoundaryData( zero );
   solver.solve();

   mfem::GridFunction const & psi = solver.potential();
   mfem::GridFunction const & q   = solver.flux();

.. note::

   ``meq/meq.hpp`` is an umbrella header covering the configuration, the
   solver, the profiles and the sources. The geometry, estimator, sampler,
   warm-start, output and flux-surface headers are **not** in it and must be
   included directly — ``meq/Estimator.hpp``, ``meq/BoundaryShape.hpp``,
   ``meq/Sampler.hpp``, ``meq/WarmStart.hpp``, ``meq/Output.hpp``,
   ``meq/Field.hpp``, ``meq/SourceFactory.hpp``, ``meq/CriticalPoints.hpp``,
   ``meq/FluxSurfaces.hpp``, ``meq/SurfaceAverage.hpp``, ``meq/Zernike.hpp``,
   ``meq/SurfaceFit.hpp``.

   The last five have no hand-written reference page here yet; their headers
   carry the reasoning, and :doc:`../flux_surfaces` and
   :doc:`../surface_geometry` are the prose.

Everything is in ``namespace meq``. There are no nested namespaces.

.. _api-borrowing:

Borrowing, and object lifetimes
-------------------------------

.. warning::

   **meq borrows almost everything and owns almost nothing.** A mesh, a source,
   a boundary coefficient, a transfer path, an initial guess — all are held by
   reference and must outlive the object that was given them.

   :cpp:class:`meq::GradShafranovSolver` is the one to watch: it holds the mesh,
   the source, the boundary data, the initial guess and the transfer path, and a
   solve touches all of them.

A :cpp:class:`meq::NormalisedSource` is not merely borrowed but **mutated**: the
solver calls ``setNormalisation()`` on it before every residual evaluation. It
must not be shared with anything that reads it meanwhile.

.. _api-mfem-free:

The MFEM-free half
------------------

``Config``, ``Profiles``, ``Source``, ``RotatingSource``, ``SourceFactory`` and
``BoundaryShape`` deliberately do not include MFEM. They take plain ``double``
arguments, and the ``mfem::Coefficient`` adapters live with the assembly that
needs them.

That is what makes the configuration and physics layers unit-testable without
the finite element library — and, in practice, what makes continuous
integration possible at all. See :doc:`../install`.

.. toctree::
   :maxdepth: 2

   solver
   sources
   profiles
   config
   geometry
   estimator
   output
