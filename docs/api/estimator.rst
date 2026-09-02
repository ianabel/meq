Error estimation and adaptivity
===============================

``#include "meq/Estimator.hpp"``

.. cpp:namespace:: meq

.. cpp:class:: ResidualEstimator

   The residual error estimator of :cite:t:`SanchezVizuet2020adaptive`. Derives
   from MFEM's error estimator interface. Non-copyable. See :doc:`../adaptivity`.

   .. cpp:function:: ResidualEstimator( GradShafranovSolver & solverIn, Source const & sourceIn )
   .. cpp:function:: ResidualEstimator( GradShafranovSolver & solverIn, mfem::Coefficient & sourceIn )

      Both arguments are **borrowed** and must outlive the estimator.

   .. cpp:enum-class:: Term

      .. cpp:enumerator:: Divergence

         :math:`\eta_1`, the residual of :math:`-\gradbar\cdot q = F/r`.

      .. cpp:enumerator:: Constitutive

         :math:`\eta_2`, the residual of :math:`q = \gradbar\psi^\star / r`.

      .. cpp:enumerator:: FluxJump

         :math:`\eta_3`, the jump in :math:`q_h` across interior edges.

      .. cpp:enumerator:: PotentialJump

         :math:`\eta_4`, the jump in :math:`\psi^\star`.

      .. cpp:enumerator:: TraceMismatch

         :math:`\eta_5`, the trace against :math:`\psi^\star`.

   .. cpp:enum-class:: Potential

      .. cpp:enumerator:: PostProcessed

         The published estimator, built on :math:`\psi^\star`. **The default.**

      .. cpp:enumerator:: Raw

         Built on :math:`\psi_h`, which loses exactly one order at every degree.
         Kept so that the difference stays measured; see
         :doc:`../postprocessing`.

   .. cpp:enum-class:: TraceComparison

      .. cpp:enumerator:: Projected

         Compares inside :math:`M_h`. **The default**, and what restores the
         published rates; see :ref:`adaptivity-eta5`.

      .. cpp:enumerator:: Literal

         The term exactly as printed in the paper, kept as the control.

   .. cpp:function:: void setPotential( Potential potentialIn )
   .. cpp:function:: Potential potential() const
   .. cpp:function:: void setTraceComparison( TraceComparison comparisonIn )
   .. cpp:function:: TraceComparison traceComparison() const
   .. cpp:function:: void setExtraQuadratureOrder( int extraIn )

      Added to twice the degree of the potential in use; the default is 4,
      because neither :math:`F/r` nor :math:`1/r` is a polynomial.

   .. cpp:function:: void setTransferredBoundary( mfem::Array<int> const & markerIn, mfem::Coefficient * datumIn = nullptr )

      Pass the solver's own :math:`\Gamma_h` marker, which is **copied**, and
      the datum from
      :cpp:func:`GradShafranovSolver::transferredDatum`, which is
      **borrowed**. With a null datum those faces are simply excluded; with one
      supplied they are included and compared against the condition actually
      imposed. See :ref:`adaptivity-eta5`.

      The datum reads the solver's flux, so it goes stale the moment the solver
      is solved again.

   .. cpp:function:: mfem::Vector const & GetLocalErrors()
   .. cpp:function:: mfem::real_t GetTotalError() const
   .. cpp:function:: void Reset()

   .. cpp:function:: double component( Term term ) const
   .. cpp:function:: mfem::Vector const & localSquares( Term term ) const

      The per-element **squares**, which is the additive quantity a marking
      criterion needs.

   .. cpp:function:: static char const * name( Term term )

   Every computing entry point throws ``std::logic_error`` if the default
   post-processed potential is in use and the solver has not been
   post-processed.

   .. warning::

      **The estimator caches, keyed on the mesh sequence.** It recomputes when
      the mesh changes or a setter is called, and **not** when the borrowed
      solver is solved again on the same mesh — which the mesh sequence cannot
      see. Call ``Reset()`` after a second solve, or build a fresh estimator.

Marking
-------

.. cpp:function:: void markDoerfler( mfem::Vector const & localErrors, double gamma, mfem::Array<int> & marked )

   Marks a minimal set carrying ``gamma`` of the total estimated error.
   ``gamma`` in :math:`(0, 1]`. ``marked`` comes back in **decreasing order of
   local error**.

.. cpp:function:: void markMaximum( mfem::Vector const & localErrors, double gamma, mfem::Array<int> & marked )

   Marks every element with :math:`\eta_K \ge \gamma \max_K \eta_K`. ``gamma``
   in :math:`[0, 1]`, and it runs the **opposite way** to Dörfler's. ``marked``
   comes back in increasing order of element index.

.. note::

   Both take :math:`\eta_K` itself, **not** the squares, so that they take the
   same argument and neither can be fed the wrong one. ``markDoerfler`` squares
   them internally.

The companion mesh
------------------

.. cpp:class:: AdaptiveDomain

   The companion-mesh construction of :cite:t:`SanchezVizuet2020adaptive`,
   without which the transfer silently leaves the regime it is analysed in. See
   :ref:`adaptivity-companion`. Non-copyable.

   .. cpp:function:: AdaptiveDomain( mfem::Mesh const & backgroundIn, mfem::PositionFunction levelSetIn, int extraRefineIn = 1 )

      The background mesh is **copied**, because refining mutates it. The level
      set is negative inside.

   .. cpp:function:: mfem::SubMesh & computational()

      .. warning::

         **Invalidated by** :cpp:func:`refine`, which builds a new submesh. A
         solver, transfer path or grid function built on it must be rebuilt.

   .. cpp:function:: int gammaHAttribute() const
   .. cpp:function:: mfem::Array<int> const & gammaHMarker() const

      Ready to pass straight to :cpp:func:`GradShafranovSolver::setExtension`.

   .. cpp:function:: void refine( mfem::Array<int> const & marked )

      Refines the marked elements **and** every companion element that
      :math:`\Gamma` cuts and that shares an edge with a marked one, which is
      what makes the computational domain grow towards :math:`\Gamma`.

   .. cpp:function:: void refineWithoutCompanion( mfem::Array<int> const & marked )

      Deliberately omits that second half. It is what the paper says does not
      work, and it is here to be measured failing rather than argued about.

   .. cpp:function:: int lastProximityAdditions() const

      How many elements the last refinement added through the proximity rule
      rather than through the indicator. **Read this before concluding that the
      estimator is concentrating on the boundary.**

   .. cpp:function:: int numComputational() const
   .. cpp:function:: int numCompanion() const
   .. cpp:function:: int numBackground() const
   .. cpp:function:: int refinements() const
   .. cpp:function:: double largestElement() const
   .. cpp:function:: double smallestElement() const

Mesh measures
-------------

.. cpp:function:: double elementDiameter( mfem::Mesh & mesh, int element )
.. cpp:function:: double faceDiameter( mfem::Mesh & mesh, int face )

   The largest distance between two vertices — the :math:`h_K` and :math:`h_e`
   the estimator is written in terms of.

   .. note::

      These are deliberately **not** MFEM's own element size, which is a
      Jacobian-derived length scale and a different number on an anisotropic
      triangle.
