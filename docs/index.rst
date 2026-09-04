meq
===

**meq** — the Maryland Equilibrium Solver — computes axisymmetric plasma
equilibria by solving the **Grad–Shafranov** equation with a hybridizable
discontinuous Galerkin (HDG) discretisation built on `MFEM
<https://mfem.org>`_.

For a poloidal flux :math:`\psi(r, z)` on a domain :math:`\Omega \subset
\mathbb{R}^2` whose boundary :math:`\Gamma` is the level set :math:`\psi = 0`,
it solves

.. math::

   -\gradbar \cdot \left( \frac{1}{r} \gradbar \psi \right)
       = \frac{F(r, z, \psi)}{r} \quad \text{in } \Omega,
   \qquad \psi = 0 \ \text{ on } \Gamma,

where :math:`\gradbar := (\partial_r, \partial_z)` and the source

.. math::

   F(r, z, \psi) := \mu_0 r^2 \frac{\mathrm{d}p}{\mathrm{d}\psi}
                  + g \frac{\mathrm{d}g}{\mathrm{d}\psi}

is built from the plasma pressure :math:`p(\psi)` and the toroidal field
function :math:`g(\psi)`. Because those are functions of :math:`\psi`, the
problem is semi-linear, and meq solves it by **Newton** rather than by the
Picard iteration the source papers use — a choice with consequences that reach
into every part of the code, and which :doc:`nonlinear` explains.

What meq does
-------------

* **Solves the fixed-boundary problem to optimal order.** Both the flux
  :math:`\psi` and its scaled gradient :math:`q = \gradbar\psi / r` converge at
  :math:`k+1` in the polynomial degree :math:`k`. The gradient is a *solved
  unknown*, not a post-processed derivative, which is the reason to prefer a
  mixed method here: the magnetic field is what most users want out of an
  equilibrium code, and it arrives at the same order as the potential.
* **Handles curved plasma boundaries** without meshing them, by solving on a
  polygonal subdomain and *transferring* the boundary condition outward along
  paths — the technique of :cite:t:`CockburnSolano2012`, as applied to
  Grad–Shafranov by :cite:t:`SanchezVizuet2020adaptive`. See
  :doc:`curved_boundary`.
* **Refines the mesh adaptively** against a residual error estimator, in the
  usual solve → estimate → mark → refine loop. See :doc:`adaptivity`.
* **Takes profiles in normalised flux**, in which case the value of
  :math:`\psi` on the magnetic axis is itself an unknown of the nonlinear
  system rather than an input. See :doc:`normalised_flux`.
* **Includes toroidal rotation**, solving the generalised Grad–Shafranov
  equation for a rotating multi-species plasma. See :doc:`rotation`.
* **Extracts the flux surfaces.** The magnetic axis and any X-point as roots of
  the solved flux, the surfaces themselves as traced curves, flux-surface
  averages over them, and the whole family as a smooth map from a disc. See
  :doc:`flux_surfaces` and :doc:`surface_geometry`.
* **Writes three output formats**, deliberately, because no one of them is both
  exact and portable. See :doc:`output`.

meq is usable two ways: as the ``meq`` program, driven by a TOML configuration
file, and as a C++ library whose central class is
:cpp:class:`meq::GradShafranovSolver`.

.. note::

   meq is **serial** throughout, and 2D. There is no MPI in it. The problem
   sizes it is built for — a few tens of thousands of trace degrees of freedom
   — are ones a sparse direct solver handles comfortably on one core, and
   :doc:`linear_solvers` explains why that is the right trade at this size
   rather than a limitation waiting to be lifted.

Where to start
--------------

Building it and running something is :doc:`install` then :doc:`running`.
:doc:`configuration` is the reference for every key in a configuration file.

If you want to understand what the solver is doing — and in particular the two
sign conventions, which are the most common way to get an extension of meq
subtly wrong — read :doc:`formulation` first. If you are choosing a nonlinear
solver strategy for a problem that is not converging, go to :doc:`nonlinear`.

.. admonition:: A note on numbers in these pages

   Many of meq's design choices could only be settled by measurement, and this
   documentation says so where that is the case. It does **not** reproduce the
   measurements. Where a choice is exposed as a configurable option — the trace
   solver, the assembly mode, the nonlinear ordering, the globalisation — the
   right setting depends on your problem, your mesh and your machine, and the
   only way to find it is to measure it there. Treat the recommendations here
   as the defaults they are, and measure before departing from them.

.. toctree::
   :maxdepth: 1

   introduction

.. toctree::
   :maxdepth: 2
   :caption: Using meq

   install
   running
   configuration
   output
   examples

.. toctree::
   :maxdepth: 2
   :caption: How it works

   formulation
   postprocessing
   curved_boundary
   adaptivity
   flux_surfaces
   surface_geometry

.. toctree::
   :maxdepth: 2
   :caption: The physics input

   sources
   profiles
   normalised_flux
   rotation

.. toctree::
   :maxdepth: 2
   :caption: Solving the system

   nonlinear
   linear_solvers

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/index
   organization
   bibliography

.. toctree::
   :maxdepth: 2
   :caption: Development

   testing
