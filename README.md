# MEQ — the Maryland Equilibrium Solver

MEQ computes axisymmetric plasma equilibria by solving the **Grad–Shafranov**
equation with a high-order hybridizable discontinuous Galerkin (HDG) method,
following the algorithms of Sánchez-Vizuet, Solano and Cerfon. It is built on the
[MFEM](https://mfem.org) finite element library.

The equation, in the fixed-boundary form MEQ solves, is

$$-\bar\nabla\cdot\left(\frac{1}{r}\bar\nabla\psi\right) = \frac{F(r,z,\psi)}{r}
\quad\text{in }\Omega,\qquad \psi = 0 \ \text{ on } \partial\Omega,$$

with

$$F(r,z,\psi) := \mu_0 r^2 \frac{\mathrm{d}p}{\mathrm{d}\psi} + g\frac{\mathrm{d}g}{\mathrm{d}\psi},$$

where $\psi$ is the poloidal flux function, $p(\psi)$ the plasma pressure and
$g(\psi)/r$ the toroidal field function. The dependence of $p$ and $g$ on $\psi$
makes the problem semi-linear.

MEQ solves it as a first-order system in $\psi$ and the flux
$\boldsymbol{q} = \bar\nabla\psi / r$. That is deliberate: the physically
interesting output is the magnetic field, which is built from $\bar\nabla\psi$,
and a mixed method resolves it to the same order as $\psi$ itself rather than one
order lower.

## Status

**MEQ solves the linear Grad–Shafranov equation**, reproducing an exact Solov'ev
equilibrium at order $k+1$ in both $\psi$ and the flux $\boldsymbol{q}$, for
$k = 1,2,3$. The semi-linear problem — a source depending on $\psi$, solved by
Newton — is in progress, and the command-line driver is not yet ported, so there
is currently no way to run MEQ except through its test suite.

`CLAUDE.md` has the full picture, including a stage-by-stage account of what is
done and what is next. Read it before believing anything below about behaviour.

Two things worth knowing up front, because earlier revisions of this file claimed
otherwise:

* MEQ is now a **fixed-boundary** solver. The free-boundary work — the von
  Hagenow / Lackner Green's-function scheme — is unported and sits in
  `attic/free-boundary/`, with a README explaining why it will not simply be
  revived.
* **The pre-port code never solved a Grad–Shafranov equation.** It computed the
  vacuum field from coil currents; the right-hand side took $\psi$ and ignored
  it, and the plasma model was never connected. Treat any description of MEQ
  predating this note with suspicion.

## Building

Requirements:

| | |
|---|---|
| C++17 compiler | g++ 15 is what it is developed against |
| CMake | ≥ 3.20 |
| [MFEM](https://mfem.org) | **4.9.1, the HDG branch** — see below |
| [toml11](https://github.com/ToruNiina/toml11) | vendored as a submodule at `extern/toml11` |
| SuiteSparse | UMFPACK and friends, pulled in through MFEM |
| Boost | `unit_test_framework`, for the test suite only |

MEQ needs a specific MFEM: the branch carrying `DarcyForm` and the HDG
integrators under `fem/darcy/`. `mfem/master` will not work — the older
`HDGBilinearForm` API that MEQ used to depend on has been replaced.

```sh
git submodule update --init --recursive
cmake -B build -DMFEM_DIR=/path/to/mfem-hdg-dev
cmake --build build -j4
```

toml11 is a pinned submodule, so `--recursive` is not optional. MFEM is not a
submodule — its history is very large and it needs its own out-of-tree build —
so point `MFEM_DIR` at a checkout you have built yourself.

`MFEM_DIR` defaults to `../mfem-hdg-dev` and also reads the environment.

## Running

The driver is not yet ported, so this does not work today:

```sh
./build/apps/meq examples/soloviev-nstx.toml
```

**The environment variable is not optional on a machine whose BLAS resolves to
MKL** — without it, MKL silently returns wrong results from UMFPACK's BLAS-3
calls. Silently: you get numbers, and they are wrong. The test suite sets it for
you; a direct invocation does not.

Run configurations are TOML. See `examples/` for worked cases, including the
Solov'ev benchmark whose exact solution is published.

## Testing

```sh
cd build && ctest --output-on-failure
```

The suite is in three layers, described in full in `CLAUDE.md`:

* **unit** — configuration parsing, spline interpolation and its derivative,
  and the source term's $\partial F/\partial\psi$ against a finite difference.
* **exact solution** — absolute error against a published closed-form Solov'ev
  equilibrium.
* **convergence** — assertions on the observed *rate*, not on a tolerance:
  $k+1$ for both $\psi$ and $\bar\nabla\psi$, over several dyadic refinements.

That last layer is the one that matters, and the ordering is deliberate. A wrong
sign convention converges at the right rate to the wrong function, so a
convergence study alone cannot catch it — only comparison against a closed form
can.

## Numerical method

The discretisation is the LDG-H formulation of Sánchez-Vizuet & Solano, with
equal-degree polynomial spaces for the flux, the potential and the hybrid trace,
and stabilisation $\tau = 1$. Hybridization condenses the element interiors away,
leaving a global system in the trace unknown alone; the interior solution is then
recovered element by element, in parallel.

Curved plasma boundaries are handled by **extension from polygonal subdomains**:
the computational domain is the union of background mesh elements lying inside
$\Omega$, and boundary data is carried from the true boundary to the mesh
boundary along transfer paths. Consequently MEQ needs only a uniform,
shape-regular background mesh of a box — there is no geometry-conforming meshing
step, and no mesh generator dependency.

`refs/Refs.md` indexes the papers, with a doi for each. The PDFs themselves are
not in the repository.

## Parallelism

None yet. The method is well suited to it — the local solves are independent by
construction — but MEQ builds against a serial MFEM.

## Built with

* [MFEM](https://github.com/mfem/mfem) — the finite element framework MEQ is
  built on, from Lawrence Livermore National Laboratory
* [toml11](https://github.com/ToruNiina/toml11) — for parsing
  [TOML](https://github.com/toml-lang/toml) configuration files
* [Boost](https://boost.org) — the unit testing framework

## Authors

* **Ian Abel** — *Initial work* — [Ian Abel at UMD](https://ireap.umd.edu/faculty/abel)

For full copyright attribution see [COPYRIGHT](COPYRIGHT); for contributors, the
[contributors](https://github.com/ianabel/meq/contributors) page.

## Acknowledgments

* The [MFEM](https://github.com/mfem) team at
  [LLNL](https://computing.llnl.gov) for a great C++ FEM library.
* The HDG algorithm for the Grad–Shafranov equation is due to T. Sánchez-Vizuet
  (New York University) and M. Solano (Universidad de Concepción), with the
  adaptive extension developed together with A. Cerfon.

## Licence

3-Clause BSD — see [LICENCE.md](LICENCE.md).
