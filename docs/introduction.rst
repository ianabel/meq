Introduction
============

MEQ solves the **fixed-boundary Grad–Shafranov problem**: given a plasma
boundary and the profiles that drive the plasma, find the poloidal flux
function whose level sets are the flux surfaces.

The problem
-----------

An axisymmetric magnetohydrodynamic equilibrium in cylindrical coordinates
:math:`(r, \varphi, z)` is described entirely by a scalar, the poloidal flux per
radian :math:`\psi(r, z)`. Force balance reduces to one semi-linear elliptic
equation for it,

.. math::

   -\gradbar \cdot \left( \frac{1}{r} \gradbar \psi \right)
       = \frac{F(r, z, \psi)}{r},
   \qquad
   F(r, z, \psi) := \mu_0 r^2 \frac{\mathrm{d}p}{\mathrm{d}\psi}
                  + g \frac{\mathrm{d}g}{\mathrm{d}\psi},

with :math:`p(\psi)` the plasma pressure and :math:`g(\psi)/r` the toroidal
field function. The operator on the left is the Grad–Shafranov operator
:math:`-\dstar`, and it is *not* the Laplacian: the weights :math:`1/r` and
:math:`r` are what distinguish it, and :math:`\gradbar := (\partial_r,
\partial_z)` acts formally like a vector of partial derivatives rather than
being the cylindrical gradient. :doc:`formulation` is precise about this,
because the distinction is the entire content of those weights.

*Fixed boundary* means the plasma boundary :math:`\Gamma` is given rather than
found. Taking it to be the level set :math:`\psi = 0` makes this an interior
Dirichlet problem. The complementary *free boundary* problem — in which the
boundary is determined by external coils and the plasma finds its own shape —
is a harder problem that MEQ does not currently solve;
:cite:t:`Serino2024` note that the fixed problem is "significantly easier",
which is worth keeping in mind when comparing against a free-boundary code.

What makes it a numerical problem
---------------------------------

Three things, and MEQ's design is a response to each.

**The physically interesting output is a derivative.** What experiments and
downstream codes want is the magnetic field, which is built from
:math:`\gradbar\psi`. A method that computes :math:`\psi` to order :math:`k+1`
and then differentiates it delivers the field at order :math:`k` — one order
worse, in the quantity people actually use. MEQ therefore uses a **mixed**
method, in which the scaled gradient :math:`q = \gradbar\psi / r` is an
independent unknown solved for directly, at the same order as :math:`\psi`.
That is the reason to prefer HDG here, and it keeps paying off in places nobody
designed for: it is what continues the solution into the band outside a curved
mesh (:ref:`output-band`), and it is what makes the flux surfaces themselves
cheap to extract accurately — the magnetic axis is a *root* of it rather than a
turning point of something differentiated. See :ref:`flux-surfaces-q`.

**The equation is semi-linear**, because :math:`p` and :math:`g` are functions
of :math:`\psi`. Somebody has to iterate. The published HDG Grad–Shafranov
solvers use Anderson-accelerated Picard iteration, keeping :math:`F` as opaque
problem data; MEQ uses **Newton**, which is faster and more robust but requires
every source to be able to differentiate itself. :doc:`nonlinear` sets out that
trade and its consequences, which reach into every part of the code.

**The plasma boundary is curved and the mesh is not.** A polygonal
approximation to a smooth boundary caps the convergence rate whatever the
polynomial degree, because the interior angle of the polygon introduces a
singularity that no amount of accuracy in the elements removes. MEQ does not
mesh the boundary at all: it solves on a polygonal subdomain strictly inside
:math:`\Gamma` and *transfers* the boundary condition outward along short
paths, which is :cite:t:`CockburnSolano2012`'s technique. See
:doc:`curved_boundary`.

What MEQ is built on
--------------------

The discretisation is the LDG-H hybridizable discontinuous Galerkin method of
:cite:t:`SanchezVizuetSolano2019`, with the adaptivity, error estimator and
curved-boundary treatment of :cite:t:`SanchezVizuet2020adaptive`. Those two
papers are not background: they *are* the method, and ``src/meq`` is an
implementation of them.

The finite element machinery is `MFEM <https://mfem.org>`_
:cite:p:`Anderson2021`, and specifically its ``DarcyForm`` and the HDG
integrators that go with it. MEQ is an early user of that work, and it requires
a development branch rather than a release; :ref:`install-mfem` says why.

.. _introduction-related:

Related codes, and where MEQ sits among them
--------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - Code
     - Discretisation
     - Relation to MEQ
   * - :cite:t:`SanchezVizuetSolano2019`, :cite:t:`SanchezVizuet2020adaptive`
     - HDG, LDG-H
     - The method MEQ implements. MEQ departs from them in using Newton rather
       than Picard, and in a handful of sign conventions that follow from the
       library it is built on.
   * - CEDRES++ :cite:p:`Heumann2015`
     - Continuous Galerkin, free boundary
     - The source of three design rules MEQ follows: differentiate the
       *discrete* residual rather than the continuous equation; expect
       normalised-flux profiles to make the Jacobian non-local; and expect
       fixed-point iteration to fail where Newton does not.
   * - :cite:t:`Serino2024`
     - Continuous Galerkin, free boundary, MFEM
     - Newton on Grad–Shafranov, in the same library. Reports Newton succeeding
       where "conventional Picard-based solvers fail to converge".
   * - :cite:t:`Pataki2013`, :cite:t:`PalhaKorenFelici2016`
     - Integral equation; mimetic spectral element
     - Other high-order approaches to the same equation.
   * - :cite:t:`TakedaTokuda1991`
     - —
     - The standing review of tokamak equilibrium computation.

The distinctive combination in MEQ is **a mixed high-order method, on a mesh
that does not conform to the plasma boundary, driven by Newton**. Each of those
three is well established on its own; taking all three at once is what produces
most of the interesting behaviour documented here, and most of the traps.

Reading this documentation
--------------------------

:doc:`formulation` is the one chapter to read before extending MEQ — it
contains the two sign conventions, which are the most common way to get an
extension of it subtly wrong.

The remaining "how it works" chapters can be read on demand:
:doc:`postprocessing` for the superconvergent reconstruction the error
estimator is built on, :doc:`curved_boundary` for the transfer technique, and
:doc:`adaptivity` for the refinement loop. :doc:`nonlinear` is about solving
the system and is the place to go when a run does not converge.

:doc:`flux_surfaces` and :doc:`surface_geometry` are about what happens *after*
a solve: turning :math:`\psi(r, z)` into the flux surfaces, the flux-surface
averages and the geometry a transport code reads. Both are library-only — none
of it is driven from a configuration file or written to an output file.

For the physics input — what a source term is, what a profile must be able to
do, and what changes when profiles are given in normalised flux or the plasma
is rotating — see :doc:`sources`, :doc:`profiles`, :doc:`normalised_flux` and
:doc:`rotation`.

.. note::

   **This documentation explains why MEQ is the way it is, and many of those
   reasons are measurements.** Where a choice could only be settled by
   measurement, these pages say so and say which way it came out; they do not
   reproduce the numbers, which belong with the tests that produce them and go
   stale as soon as anything changes. Where the choice is exposed as an option,
   the recommendation is to measure it on your own problem — the settings that
   are right here were established on one machine, on one family of test cases.
