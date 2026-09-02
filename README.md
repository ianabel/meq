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

## Documentation

Full documentation — installation, the configuration reference, the C++ API, and
an account of the numerics and why they are the way they are — is at
**<https://meq.readthedocs.io>**, built from `docs/`.

To build it locally:

```sh
cmake --build build --target docs        # or: make -C docs html
```

`CLAUDE.md` is a separate thing and stays: it is the working technical record,
including measurements, dead ends and the reasoning behind decisions. The
documentation is what a user needs; `CLAUDE.md` is what a maintainer needs.

## Status

**MEQ solves the semi-linear Grad–Shafranov equation by Newton**, reproducing an
exact Solov'ev equilibrium and a manufactured non-linear solution at order
$k+1$ in both $\psi$ and the flux $\boldsymbol{q}$, with Newton converging
quadratically. It is a program as well as a library: `meq config.toml` builds
the mesh and the source, solves, and writes the equilibrium in three formats.

What is present, in addition to the linear solve:

* **Newton on a source depending on $\psi$**, with every source supplying its own
  $\partial F/\partial\psi$.
* **Curved plasma boundaries**, by extension from polygonal subdomains, driven
  from a configuration file.
* **Adaptive refinement** against a residual error estimator, on both the fitted
  and the curved path.
* **Profiles in normalised flux**, in which the flux on the magnetic axis is an
  unknown of the system rather than an input, closed by a bordered Newton.
* **Toroidal rotation**, for an arbitrary number of species.

Not present: free boundary, and any parallelism.

`CLAUDE.md` has the full picture, including a stage-by-stage account of what is
done and what is next, and the measurements behind every claim above.

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
| [MFEM](https://mfem.org) | **the HDG branch** — see below. The one hard dependency |
| [toml11](https://github.com/ToruNiina/toml11) | vendored as a submodule at `extern/toml11` |
| netcdf-cxx4 | optional; without it the gridded output format is dropped |
| Boost | `unit_test_framework`, for the test suite only |
| clang-tidy | optional; gates the `naming` test only |

Reached *through* MFEM rather than found separately: SuiteSparse (UMFPACK, the
default direct solver), SUNDIALS (KINSOL, which backs the Picard globalisations),
GSLIB (point location, for interpolating a warm start between meshes), and
optionally MKL PARDISO and cuDSS. MEQ reads which of them are present out of
MFEM's own configuration and adapts.

MEQ needs a specific MFEM: the branch carrying `DarcyForm` and the HDG
integrators under `fem/darcy/`. `mfem/master` will not work — the older
`HDGBilinearForm` API that MEQ used to depend on has been replaced.

```sh
git submodule update --init --recursive
cmake -B build -DMFEM_DIR=/path/to/mfem/install
cmake --build build -j4
```

Never a bare `make -j` or `cmake --build -j` with no number: unbounded, the job
count goes to the host's core count, which on a container or a WSL2 virtual
machine is not the same thing as the memory behind it.

toml11 is a pinned submodule, so `--recursive` is not optional. MFEM is not a
submodule — its history is very large and it needs its own out-of-tree build —
so point `MFEM_DIR` at a checkout you have built yourself.

`MFEM_DIR` defaults to `../mfem/install` and also reads the environment, so on
a machine laid out that way the `-D` is unnecessary. The finder accepts either
layout an MFEM build leaves behind — an installed tree or an in-source make
build — and reads the enabled features out of MFEM's own `config.mk`.

`find_package(MFEM)` is deliberately **not** `REQUIRED`. Without it CMake builds
the MFEM-free half of the library — the configuration parser, the profiles, the
sources, the boundary shapes — and skips every test that needs a solver. That is
what continuous integration builds, since the MFEM branch MEQ needs is published
nowhere a hosted runner can fetch it; **read a green CI badge as evidence about
those layers and nothing else.**

## Running

```sh
./build/meq examples/soloviev-nstx.toml
```

One positional argument and no other flags, apart from `--help` and `--version`:
everything about a run is in the configuration file, which is a record of what
was run in a way a command line is not. Exit codes are 0 solved, 1 configuration,
2 did not converge, 3 could not write the output.

Run configurations are TOML. See `examples/` for worked cases — the Solov'ev
benchmark, a genuinely non-linear source, a curved boundary, an adaptive run, and
two rotating equilibria — and the [configuration
reference](https://meq.readthedocs.io/en/latest/configuration.html) for every key.

**Set `MKL_NUM_THREADS=1` on a machine whose BLAS is MKL.** MEQ's inner loop
factorises a great many *small* dense matrices, and above a block size that
depends on your BLAS each of those calls pays for a thread fork and a barrier
that dwarfs the arithmetic — measured here, the effect was dramatic and grew
rapidly with polynomial degree. Every registered test sets the variable for
itself; a direct invocation does not.

## Testing

```sh
cd build && ctest --output-on-failure
```

The suite is in three layers, described in full in the
[testing chapter](https://meq.readthedocs.io/en/latest/testing.html):

* **unit** — configuration parsing, spline interpolation and its derivative,
  and the source term's $\partial F/\partial\psi$ against a finite difference.
* **exact solution** — absolute error against a published closed-form Solov'ev
  equilibrium.
* **convergence** — assertions on the observed *rate*, not on a tolerance:
  $k+1$ for both $\psi$ and $\bar\nabla\psi$, over several dyadic refinements.

That last layer is the one that matters, and the ordering is deliberate. A wrong
sign convention converges at the right rate to the wrong function, so a
convergence study alone cannot catch it — only comparison against a closed form
can. A wrong *Jacobian* is worse: it changes neither the answer nor the rate,
only the observed order of the Newton iteration, which is why that order is
asserted and printed.

One consequence of the stance is worth stating here rather than being met by
surprise: **a test asserts the behaviour that is wanted and fails until it is
there**, never the reverse, so the suite is expected to be red while a known
defect stands. That is the intended signal, and the failing test's message is
the record.

## Numerical method

The discretisation is the LDG-H formulation of Sánchez-Vizuet & Solano, with
equal-degree polynomial spaces for the flux, the potential and the hybrid trace,
and stabilisation $\tau = 1$. Hybridization condenses the element interiors away,
leaving a global system in the trace unknown alone; the interior solution is then
recovered element by element, independently.

Curved plasma boundaries are handled by **extension from polygonal subdomains**:
the computational domain is the union of background mesh elements lying inside
$\Omega$, and boundary data is carried from the true boundary to the mesh
boundary along transfer paths. Consequently MEQ needs only a uniform,
shape-regular background mesh of a box — there is no geometry-conforming meshing
step, and no mesh generator dependency.

Curved boundaries and adaptive refinement interact, and the interaction is the
subtle part: refining the computational mesh does not move its boundary, so the
gap to the true boundary must be closed by a companion mesh or the transfer
silently leaves the regime it is analysed in. The
[documentation](https://meq.readthedocs.io/en/latest/adaptivity.html) covers it.

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
