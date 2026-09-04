The solver
==========

``#include "meq/GradShafranov.hpp"``

.. cpp:namespace:: meq

.. cpp:class:: GradShafranovSolver

   The Grad–Shafranov solver. Non-copyable.

   .. cpp:function:: GradShafranovSolver( mfem::Mesh & meshIn, int orderIn, double tauIn = 1.0 )

      Builds the three finite element spaces of :ref:`formulation-spaces` on
      ``meshIn``, all of degree ``orderIn``.

      The mesh is **borrowed** and must outlive the solver, and :math:`r > 0`
      is required everywhere on it. Throws ``std::invalid_argument`` for a
      negative order or a mesh that is not two dimensional.

Setting the source
------------------

Exactly one source per solver, set before anything is assembled.

.. cpp:namespace-push:: GradShafranovSolver

.. cpp:function:: void setSource( mfem::Coefficient & fIn )

   The **linear** path: :math:`F(r, z)`, independent of :math:`\psi`. Note that
   this is :math:`F`, not :math:`F/r`.

.. cpp:function:: void setSource( Source const & fIn )

   The **semi-linear** path, solved by Newton. See :doc:`../sources`.

.. cpp:function:: void setSource( NormalisedSource & fIn, double psiAxisGuessIn )

   The **bordered** path, in which :math:`\psiax` is an unknown of the system
   and ``psiAxisGuessIn`` is where Newton starts looking for it. See
   :doc:`../normalised_flux`. The source is borrowed **and mutated**.

All three throw ``std::logic_error`` if the forms are already built, or if any
source is already set — including a second source of the same kind. The
normalised overload additionally throws ``std::invalid_argument`` for a guess
that is not finite or is zero.

Boundary data and initial guess
-------------------------------

.. cpp:function:: void setBoundaryData( mfem::Coefficient & boundaryIn )

   The Dirichlet datum. Required unless every boundary attribute is on
   :math:`\Gamma_h`.

.. cpp:function:: void setInitialGuess( mfem::Coefficient & psiGuess )
.. cpp:function:: void setInitialGuess( mfem::GridFunction const & psiGuess )
.. cpp:function:: void clearInitialGuess()
.. cpp:function:: bool hasInitialGuess() const

   The guess reaches the trace space by an :math:`L^2` projection onto each
   face — note that projecting a coefficient onto a grid function does **not**
   do this, since it loops over volume elements and never touches a face degree
   of freedom.

   **Order matters**: the Dirichlet datum is applied *after* the guess, so the
   boundary condition always wins on essential degrees of freedom. The guess is
   ignored on the linear path, where the solve is direct.

   The grid function overload requires the guess to be on the solver's own mesh;
   use :cpp:class:`meq::FieldTransfer` first if it is not.

The curved boundary
-------------------

.. cpp:function:: void setExtension( mfem::TransferPath & pathIn, mfem::Array<int> const & gammaHMarkerIn, int lineOrderIn = -1 )

   Marks the boundary attributes on which the Dirichlet datum is *transferred*
   rather than imposed. See :doc:`../curved_boundary`. The path is borrowed. A
   negative ``lineOrderIn`` takes the integrator's own default.

   Throws ``std::logic_error`` if the forms are built, and
   ``std::invalid_argument`` if the marker is not sized by the mesh's largest
   boundary attribute or selects nothing.

.. cpp:function:: bool isExtended() const

.. cpp:function:: std::unique_ptr<mfem::Coefficient> transferredDatum( mfem::PositionFunction g = mfem::PositionFunction() )

   Rebuilds :math:`\varphi_h`, the datum actually imposed on :math:`\Gamma_h`,
   which is not otherwise stored anywhere. ``g`` is the datum on the **true**
   boundary as a function of position; an empty ``g`` means zero. Throws
   ``std::logic_error`` on the fitted path.

   This is what the error estimator needs; see :ref:`adaptivity-eta5`.

Normalisation
-------------

.. cpp:function:: bool normalisationIsUnknown() const
.. cpp:function:: double psiAxis() const
.. cpp:function:: double normalisationResidual() const

   ``psiAxis()`` is the guess before ``solve()`` and the converged value after.
   ``normalisationResidual()`` is the constraint residual
   :math:`\psiax - \max\psi_h`. Both are zero unless :math:`\psiax` is an
   unknown.

.. cpp:function:: double axisFlux( mfem::Vector const & trace, int * element = nullptr )

   The axis flux implied by a given trace. Non-const, and valid only after
   ``prepare()``; it leaves the solution blocks alone. ``element`` receives the
   element that attained the maximum.

.. cpp:enum-class:: Normalisation

   .. cpp:enumerator:: Coupled

      The bordered Newton with its two non-local terms. **The default, and the
      only one to run a calculation with.**

   .. cpp:enumerator:: Decoupled

      The same iteration with the border dropped. It exists so the value of the
      non-local terms can be *measured*; see :ref:`normalised-decoupled`.

.. cpp:function:: void setNormalisationCoupling( Normalisation choice )
.. cpp:function:: Normalisation normalisationCoupling() const

Choosing algorithms
-------------------

.. cpp:enum-class:: Globalisation

   .. cpp:enumerator:: None

      Plain Newton, full steps. **The default.**

   .. cpp:enumerator:: LineSearch

      KINSOL with backtracking.

   .. cpp:enumerator:: KinsolNoLineSearch

      KINSOL with full steps, so that a difference against ``LineSearch`` is
      attributable to the line search rather than to SUNDIALS.

   .. cpp:enumerator:: AndersonPicard

      Anderson-accelerated fixed point — the source papers' own method. The
      potential block is assembled linearly.

   .. cpp:enumerator:: PicardOnly

      The same fixed point without acceleration.

   .. cpp:enumerator:: PicardThenNewton

      Picard into Newton's basin, then plain Newton. See
      :ref:`nonlinear-globalisation`.

.. cpp:enum-class:: NonlinearOrdering

   .. cpp:enumerator:: CondenseThenLinearise

      Condense first, so each element's elimination is itself a nonlinear
      solve. MEQ's **backup**, and the more robust of the two on a coarse mesh.

   .. cpp:enumerator:: NPC

      Linearise the full system, then hybridize. **MEQ's default.** See
      :ref:`nonlinear-ordering`.

.. cpp:enum-class:: LocalSolver

   .. cpp:enumerator:: Newton
   .. cpp:enumerator:: Lbfgs
   .. cpp:enumerator:: Lbb

   The element-local nonlinear solver. **Entirely inert under**
   :cpp:enumerator:`NonlinearOrdering::NPC`, which has no local nonlinear solve
   to configure.

.. cpp:enum-class:: AssemblyMode

   .. cpp:enumerator:: Serial

      **The default, unconditionally.**

   .. cpp:enumerator:: Threaded

      Threads the element-local half of the assembly. Requires MFEM built with
      OpenMP and thread safety. See :ref:`linear-threading` before taking it.

.. cpp:enum-class:: TraceSolver

   .. cpp:enumerator:: UMFPack

      **The default**, and the only backend present in every build.

   .. cpp:enumerator:: Pardiso
   .. cpp:enumerator:: cuDSS

      Spelled the way its vendor spells it, deliberately breaking the house
      naming rule. ``cuDSS`` needs an ``mfem::Device`` configured for CUDA
      **before** the solver is built.

.. cpp:function:: void setGlobalisation( Globalisation choice )
.. cpp:function:: Globalisation globalisation() const
.. cpp:function:: void setNonlinearOrdering( NonlinearOrdering choice )
.. cpp:function:: NonlinearOrdering nonlinearOrdering() const
.. cpp:function:: void setLocalSolver( LocalSolver choice )
.. cpp:function:: void setAssemblyMode( AssemblyMode choice )
.. cpp:function:: AssemblyMode assemblyMode() const
.. cpp:function:: void setTraceSolver( TraceSolver choice )
.. cpp:function:: TraceSolver traceSolver() const
.. cpp:function:: static bool traceSolverAvailable( TraceSolver choice )

   The setters throw rather than falling back silently: ``setGlobalisation``
   throws ``std::logic_error`` for a KINSOL option in a build without SUNDIALS,
   and ``setAssemblyMode`` and ``setTraceSolver`` throw
   ``std::invalid_argument`` when the backing facility is absent.
   ``traceSolverAvailable`` answers the question without throwing.

.. cpp:function:: void setPicardDamping( double damping )
.. cpp:function:: void setAndersonDepth( int depth )

   Damping must be in :math:`(0, 1]` and defaults to 1.0 — right for Anderson
   and wrong for plain Picard. The Anderson depth defaults to **1**, not to the
   published 2, which does not work here. See
   :ref:`nonlinear-globalisation`.

.. cpp:function:: void setNewtonControl( double relativeToleranceIn, double absoluteToleranceIn, int maxIterationsIn )
.. cpp:function:: bool isNonlinear() const
.. cpp:function:: bool usesNonlinearForms() const

Solving
-------

.. cpp:function:: void prepare()

   Assembles and reduces to the trace system without solving. Public so that
   the Jacobian can be checked against a finite difference of the residual it
   claims to differentiate. Throws ``std::logic_error`` if no source is set, or
   if boundary data is missing while some attribute is still fitted.

.. cpp:function:: void solve()

   Calls ``prepare()`` first. Throws ``std::runtime_error`` if the iteration
   does not converge — **no solution is recovered from a failed iteration** —
   and ``std::logic_error`` if :math:`\psiax` is an unknown and the
   globalisation is anything but ``None``.

.. cpp:function:: void postProcess()
.. cpp:function:: bool isPostProcessed() const

   Builds :math:`\psi^\star` and the two post-processed flux quantities. See
   :doc:`../postprocessing`. A new ``solve()`` invalidates them.

Results
-------

.. cpp:function:: mfem::GridFunction & potential()
.. cpp:function:: mfem::GridFunction & flux()
.. cpp:function:: mfem::GridFunction & trace()
.. cpp:function:: mfem::GridFunction & postProcessedPotential()
.. cpp:function:: mfem::GridFunction & postProcessedFlux()
.. cpp:function:: mfem::GridFunction & totalFlux()

   Each has a ``const`` overload. The first three are valid after ``solve()``,
   the last three after ``postProcess()``.

.. warning::

   ``flux()`` **and** ``postProcessedFlux()`` **return MEQ's** :math:`+q`;
   ``totalFlux()`` **is in the library's convention**, which is :math:`-q`.
   See :ref:`sign-conventions`. On the extension path ``trace()`` is *not*
   :math:`\varphi_h` on :math:`\Gamma_h`; use ``transferredDatum()``.

.. cpp:function:: mfem::FiniteElementSpace & fluxSpace()
.. cpp:function:: mfem::FiniteElementSpace & potentialSpace()
.. cpp:function:: mfem::FiniteElementSpace & traceSpace()

.. cpp:function:: double potentialError( mfem::Coefficient & exact ) const
.. cpp:function:: double fluxError( mfem::VectorCoefficient & exact ) const
.. cpp:function:: double postProcessedPotentialError( mfem::Coefficient & exact ) const

The reduced system, and diagnostics
-----------------------------------

.. cpp:function:: mfem::Operator & reducedOperator()
.. cpp:function:: mfem::Vector & reducedRhs()
.. cpp:function:: mfem::Vector & reducedSolution()
.. cpp:function:: mfem::Array<int> const & essentialTraceDofs() const

   Available after ``prepare()``. On the semi-linear path the operator's
   gradient differentiates the **discrete** residual. The residual is masked to
   zero on the essential trace degrees of freedom and the Jacobian carries a
   unit row there, so a finite-difference Jacobian check must perturb only their
   complement.

.. cpp:function:: std::vector<double> const & newtonResiduals() const
.. cpp:function:: int newtonIterations() const
.. cpp:function:: int picardIterations() const
.. cpp:function:: long localNonlinearIterations() const
.. cpp:function:: long symbolicFactorisations() const
.. cpp:function:: long numericFactorisations() const
.. cpp:function:: int numTraceDofs() const
.. cpp:function:: int order() const
.. cpp:function:: double tau() const

   ``localNonlinearIterations()`` is **the acceptance signal for**
   :cpp:enumerator:`NonlinearOrdering::NPC` — it reads exactly zero there and in
   the thousands under the condensation, and is the only way to tell the two
   apart from outside.

   The two factorisation counts matter as a **ratio**: symbolic reuse is on, so
   a converged Newton run should analyse the sparsity once and refactorise once
   per step.

   When :math:`\psiax` is an unknown, ``newtonResiduals()`` is the *augmented*
   residual over the full bordered system; see :doc:`../normalised_flux`.

.. cpp:namespace-pop::

Supporting classes
------------------

.. cpp:class:: ConstantStabilization

   A constant HDG stabilisation :math:`\tau`, which is what the method requires
   and what MFEM does not supply by default. See :ref:`formulation-tau`.

   .. cpp:function:: explicit ConstantStabilization( double tauIn )
   .. cpp:function:: double tau() const

.. cpp:class:: SourceIntegrator

   Contributes :math:`-(F/r, w)` to the residual and
   :math:`-((\partial F/\partial\psi)\,w, v)/r` to the Jacobian.

   .. cpp:function:: explicit SourceIntegrator( Source const & sourceIn, int extraOrderIn = 4 )

      The source is **borrowed** and must outlive the integrator, which means
      outliving the solver. Residual and Jacobian use the same quadrature rule,
      so they cannot drift apart through a change to one of them.

   Like every MFEM integrator it carries scratch and is **not thread safe**.
