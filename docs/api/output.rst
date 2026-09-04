Fields, sampling, restart and output
====================================

.. cpp:namespace:: meq

The magnetic field
------------------

``#include "meq/Field.hpp"``

.. cpp:function:: void poloidalField( mfem::GridFunction const & q, mfem::GridFunction & field )

   :math:`B_R = -q_z`, :math:`B_Z = +q_r` :cite:p:`Miller1998`, applied degree
   of freedom by degree of freedom — no interpolation, no quadrature — so
   :math:`\mathbf{B}` carries :math:`q`'s own convergence order.

   .. warning::

      ``q`` **must be** :cpp:func:`GradShafranovSolver::flux`, in MEQ's sign
      convention — **not** the raw block, which holds :math:`-q`. See
      :ref:`sign-conventions`.

   ``field`` is resized and given ``q``'s space if it has none; the spaces must
   match if it has one. Throws ``std::logic_error`` if ``q`` has no space, if
   that space is not two-component, or if ``field`` already has a different one.

.. cpp:class:: PoloidalFieldCoefficient : public mfem::Coefficient

   .. cpp:function:: PoloidalFieldCoefficient( mfem::GridFunction const & qIn, int componentIn )

      ``componentIn`` is 0 for :math:`B_R` and 1 for :math:`B_Z`; anything else
      throws ``std::logic_error``. Holds ``q`` by reference, and carries scratch,
      so it is not thread safe.

Sampling onto a grid
--------------------

``#include "meq/Sampler.hpp"``

.. cpp:class:: GridSampler

   .. cpp:function:: GridSampler( mfem::Mesh & mesh, double rMinIn, double rMaxIn, int nRIn, double zMinIn, double zMaxIn, int nZIn )

      The counts are **nodes**, at least 2, so the spacing is
      :math:`(r_{\max} - r_{\min})/(n_R - 1)`. The mesh is borrowed.

      Locating is done once, at construction, by inverting the loop — element
      bounding box to grid index range by arithmetic — rather than by searching
      per point. **Being outside is not an error**: a node in the box but
      outside the mesh is simply not found, and that is the mask the output
      wants.

   Index order is ``j*nR + i``, **R fastest**, which is the ordering the NetCDF
   file uses.

   .. cpp:function:: int nodesR() const
   .. cpp:function:: int nodesZ() const
   .. cpp:function:: double rAt( int i ) const
   .. cpp:function:: double zAt( int j ) const
   .. cpp:function:: bool located( int i, int j ) const
   .. cpp:function:: int locatedCount() const

   .. cpp:function:: void sample( mfem::GridFunction const & field, std::vector<double> & values, double fill ) const
   .. cpp:function:: void sampleComponent( mfem::GridFunction const & field, int component, std::vector<double> & values, double fill ) const
   .. cpp:function:: void sampleCoefficient( mfem::Coefficient & coefficient, std::vector<double> & values, double fill ) const

      ``values`` is resized and set to ``fill`` wherever the node was not
      located.

      .. warning::

         For a node in the band beyond :math:`\Gamma_h`, these read the value at
         the **foot** on :math:`\Gamma_h` — piecewise constant, and first order.
         For anything written to an interchange file, use the two continuing
         samplers below instead. See :ref:`output-band`.

   .. cpp:function:: int extendOutward( double reach, std::function<bool( double, double )> const & accept = std::function<bool( double, double )>(), std::function<double( double, double )> const & gapToBoundary = std::function<double( double, double )>() )

      Fills nodes just outside the mesh by extrapolating from the element owning
      the nearest boundary face, and records the foot and the step. Returns how
      many were newly filled.

      ``reach`` is a multiple of the boundary face's own length — **a hard limit,
      not a tuning knob**; values much above 1 are outside anything anybody has
      analysed. ``accept`` is an optional predicate on position, which the
      curved path uses to mean "inside :math:`\Gamma`".

   .. cpp:function:: void samplePotentialWithFlux( mfem::GridFunction const & potential, mfem::GridFunction const & flux, std::vector<double> & values, double fill ) const

      Continues :math:`\psi` across the band by a Taylor step from the recorded
      foot, using :math:`\gradbar\psi = r q`. ``flux`` must be
      :cpp:func:`GradShafranovSolver::flux` — not the raw block and not the
      poloidal field.

   .. cpp:function:: void sampleComponentWithGradient( mfem::GridFunction const & field, int component, std::vector<double> & values, double fill ) const

      The same step for a vector component, using the field's own gradient.
      **It does not reach the flux's own order, and that is structural** — see
      :ref:`output-field`.

   .. cpp:function:: bool wasExtended( int i, int j ) const
   .. cpp:function:: int extendedCount() const

      Public because a reader of the output file is entitled to it: the
      ``inside`` mask says 1 in the band, so nothing else distinguishes a solved
      node from a continued one.

   .. cpp:function:: double blendWeight( int i, int j ) const

      0 on :math:`\Gamma_h`, 1 on :math:`\Gamma`, and 0 for every node found
      inside an element. Zero throughout unless ``extendOutward`` was given a
      gap function.

Warm starts
-----------

``#include "meq/WarmStart.hpp"``

.. cpp:class:: FieldTransfer

   Interpolates a field from one mesh onto another. Non-copyable.

   .. cpp:function:: explicit FieldTransfer( mfem::Mesh & sourceMeshIn )

      The mesh is borrowed and must outlive the object; it is non-const because
      the underlying point-location setup is. **Setup is the expensive part and
      is done once**, so a caller transferring several fields from the same mesh
      should keep one of these.

   .. cpp:function:: int transfer( mfem::GridFunction const & source, mfem::Coefficient & fallback, mfem::GridFunction & target )

      Fills ``target`` with ``source`` evaluated at ``target``'s nodes and
      returns how many nodes fell outside the source mesh. ``target``'s space
      decides where the source is sampled, so it must **already be sized**.
      ``fallback`` is evaluated at any uncovered node — the Dirichlet datum is
      the sensible choice.

   .. cpp:function:: int missed() const
   .. cpp:function:: int queried() const
   .. cpp:function:: double worstDistance() const

      The largest distance a point had to be moved to land in an element.
      **Large values mean the two meshes disagree about where the domain is,
      which a miss count of zero will not tell you.**

Writing files
-------------

``#include "meq/Output.hpp"``. See :doc:`../output` for what each format is for.

.. cpp:function:: bool hasNetCDF()

   Whether MEQ was built with netcdf-cxx4. Everything in
   :cpp:class:`NetCDFWriter` throws without it.

.. cpp:function:: void writeMfem( std::string const & stem, mfem::Mesh & mesh, mfem::GridFunction const & potential, mfem::GridFunction const & flux )

   Writes the exact restart triple at full precision — not the library default,
   because these are read back. ``flux`` is MEQ's :math:`+q`.

   .. important::

      Pass :math:`\psi_h` here, **not** :math:`\psi^\star`. This is the restart
      format and is read back into a degree-\ :math:`k` space. See
      :ref:`output-which-potential`.

.. cpp:function:: void writePostProcessed( std::string const & stem, mfem::GridFunction const & postProcessed )

   Writes ``<stem>_psistar.gf``. Separate from :cpp:func:`writeMfem` rather than
   folded into it, because :math:`\psi^\star` exists only after
   :cpp:func:`GradShafranovSolver::postProcess`, and the mesh and restart pair
   should not depend on a post-processing they have nothing to do with.

.. cpp:function:: void writeVtu( std::string const & stem, mfem::Mesh & mesh, mfem::GridFunction const & potential, mfem::GridFunction const & field, int levelsOfDetail )

   Produces a **directory**, not a file, with the index inside it rather than
   beside it. ``field`` is the poloidal field :math:`\mathbf{B}`, not the flux —
   use :cpp:func:`poloidalField`.

   ``levelsOfDetail`` is the subdivision of the VTK Lagrange cells and should be
   **the degree of the field being passed**, not the degree of the solve. The
   driver draws :math:`\psi^\star` and so passes :math:`k+1`; see
   :ref:`output-vtk`.

.. cpp:function:: void boundaryPolyline( mfem::Mesh & mesh, std::vector<double> & r, std::vector<double> & z, int & unreached )

   The domain boundary as an **ordered** closed polyline, walked by vertex
   adjacency, since boundary elements come out in no particular order. The loop
   is **not** repeated at the end. Only the loop containing the first boundary
   element is returned; ``unreached`` counts the boundary vertices not on it, and
   a non-zero value means the domain is not simply connected and the answer is
   partial.

.. cpp:function:: int curveBoundaryOnto( mfem::Mesh & mesh, int order, std::function<void( double, double, double &, double & )> const & project, double & applied )

   Bends the mesh boundary out onto the true :math:`\Gamma`, **for output
   only**. Returns the number of nodes moved and sets ``applied`` to the
   fraction of the displacement that survived the tangling check.

   .. warning::

      **This changes the geometry and must be the last thing done.** Call it
      after :cpp:func:`writeMfem` and after the grid sampling, immediately
      before :cpp:func:`writeVtu`.

   Moving a boundary by an element's own size can turn an element inside out, so
   the displacement is applied whole, every Jacobian is checked, and on failure
   the displacement is halved and rechecked; if no fraction works the mesh is
   left exactly as found.

.. cpp:class:: VtuSeries

   One frame per adaptive cycle, as a time series. Non-copyable.

   .. cpp:function:: VtuSeries( std::string const & stem, int levelsOfDetail )
   .. cpp:function:: void append( mfem::Mesh & mesh, mfem::GridFunction const & potential, mfem::GridFunction const & field, int cycle, double time )

      Writes immediately and retains nothing, so the caller may destroy the
      solver on the next line.

   .. cpp:function:: int frames() const

   .. note::

      This is a **separate collection** from the answer, and its frames are
      *uncurved* — :math:`\Gamma_h` as solved. The converged output has its
      boundary bent onto :math:`\Gamma`, and doing that mid-loop would hand the
      next refinement a geometry the estimator never saw.

.. cpp:class:: NetCDFWriter

   The gridded interchange file. Non-copyable.

   .. cpp:function:: NetCDFWriter( std::string const & path, GridSampler const & sampler )

      Writes the coordinate variables immediately, so the file is
      self-describing even if a later write fails.

   .. cpp:function:: void attribute( std::string const & name, std::string const & value )
   .. cpp:function:: void attribute( std::string const & name, double value )
   .. cpp:function:: void attribute( std::string const & name, int value )
   .. cpp:function:: void field( std::string const & name, std::vector<double> const & values, std::string const & longName, std::string const & units )
   .. cpp:function:: void boundary( std::vector<double> const & r, std::vector<double> const & z )
   .. cpp:function:: void close()

      The destructor calls ``close()`` and **swallows any error on purpose**,
      since a throw from a destructor terminates. Call it explicitly to see the
      error.

   A node outside the domain gets **both** a fill value of NaN and a zero in the
   ``inside`` mask, deliberately, because some tools honour the fill attribute
   and some do not.

   .. important::

      ``inside`` **and** ``extrapolated`` **answer different questions**, and a
      reader that conflates them will be misled. ``inside`` is "is there data
      here", and a band node scores 1 — correctly. ``extrapolated`` is "was this
      solved for", and the same node scores 1 there too. Drop those nodes before
      computing an error norm or differencing two runs.
