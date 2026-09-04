# Installing MEQ

**MEQ needs a development branch of MFEM that is not published anywhere.** That
is the single thing standing between a fresh clone and a working solver, and it
is why this file exists separately from `README.md`.

Everything else — toml11, Boost, Eigen, netCDF — is either vendored or a package
your distribution has.

## The short version

```sh
git clone --recursive github:ianabel/meq.git
cd meq
cmake -B build -DMFEM_DIR=/path/to/mfem/install
cmake --build build -j4
cd build && ctest --output-on-failure
```

If you do not have the right MFEM, that still works — it just builds and tests
about half of MEQ. See [Without MFEM](#without-mfem).

> **Never run a bare `make -j` or `cmake --build -j` with no number.** Unbounded,
> the job count goes to the host's core count, which on a container or a WSL2
> virtual machine is not the same thing as the memory behind it. Use `-j4` for
> MEQ, and `-j4` for MFEM, whose translation units are large.

## Why a custom MFEM

MEQ is built on `DarcyForm` and the HDG integrators under `fem/darcy/`, and on
the curved-boundary machinery in `fem/darcy/extension_hdg.*`. **None of this is
in a released MFEM, and `mfem/master` will not work** — the older
`HDGBilinearForm` API that MEQ was originally written against has been removed,
and its replacement is the branch work below.

MEQ needs pieces from **four** branches:

| branch | what MEQ needs from it |
|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — the curved boundary by extension from polygonal subdomains, and `TransferredDatumCoefficient`. **The merge base.** |
| `direct-solver-symbolic-reuse` | `UMFPackSolver` and `PardisoSolver` keeping their symbolic factorisation across Newton steps |
| `gf-hdg-linearise-first` | `DarcyNPCOperator` / `DarcyNPCSolver` — the NPC method, which is MEQ's default nonlinear ordering — plus `SetAssemblyMode`, `SetGradientMode`, `SetLocalFactorMode` |
| `gf-hdg-dev` | The post-processing fix: *"the postprocessing closes on the element average, always"* |

**That last one is not optional and is the easiest to lose**, because it is
currently an *ancestor* of two of the others and so arrives without a merge of
its own. Without it, the post-processed potential `ψ*` is silently wrong on
every element where `∂F/∂ψ` vanishes — per element, so no global norm can see
it — which degrades the error estimator, the adaptive loop, and (since
2026-09-02) the potential written to the `.vtu` and `.nc` outputs.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` is the regression
that catches it; if that test is green, your MFEM has the fix.

### The integration branch

The four are combined into a local branch called `meq-integration`, which
exists **only to be built against**. It is not a development branch, it is
published on no remote, and it must never be pushed: it is re-created whenever
any of the four moves, so anything committed directly to it is lost.

The recipe, the current merge topology, and the three conflict resolutions
(which are decided rather than fiddly, but include one seam where "keep both
sides" leaves the file a brace short) are in **`CLAUDE.md`**, under
*`meq-integration`: the branch MEQ builds from*. Read that before rebuilding.

Two checks from it are worth repeating here because they have both fired:

* **Verify the branches are actually contained** in your merge, rather than
  assuming — the topology belongs to the upstream tree and changes without MEQ
  doing anything.
* **Verify the *install*, not just the branch.** They are different questions,
  and a check that passed on the branch while the installed library was a day
  out of date has already cost this project a set of measurements:

  ```sh
  diff <( git -C ../mfem/mfem-src show meq-integration:fem/darcy/darcyhybridization.hpp ) \
       ../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp \
    && echo "install is current"
  ```

## Building MFEM

Out of tree, installed to a prefix MEQ can point at. MFEM is deliberately **not**
a submodule: its history is enormous, it needs its own configure-and-build, and
the development tree it comes from has its own active work that a submodule pin
would fight.

The options MEQ reads, and what each buys:

| option | required? | what MEQ does with it |
|---|---|---|
| `MFEM_USE_SUITESPARSE` | **effectively yes** | UMFPACK, the default trace solver. Without any direct solver MEQ falls back to GMRES. |
| `MFEM_USE_EXCEPTIONS` | **yes, for the driver** | Makes an MFEM error a C++ exception rather than an abort. Without it a failed solve is a `SIGABRT` instead of exit code 2. |
| `MFEM_USE_SUNDIALS` | for globalisation | KINSOL, behind the Anderson–Picard and Picard-then-Newton paths |
| `MFEM_USE_GSLIB` | for warm starts | `FindPointsGSLIB`, used to interpolate a stored answer onto a different mesh |
| `MFEM_USE_LAPACK` | recommended | Also gates MFEM's cut-element quadrature |
| `MFEM_USE_MKL_PARDISO` | optional | A second trace solver, selectable at run time |
| `MFEM_USE_CUDA`, `MFEM_USE_CUDSS` | optional | A third trace solver. A prerequisite, not a speed-up: `fem/darcy/` has no partial-assembly kernels and MEQ issues no device kernels of its own. |
| `MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE` | optional, **both or neither** | Threaded element assembly, which MEQ exposes as an opt-in. It *aborts* rather than falling back if only one is set. |
| `MFEM_USE_MPI` | **no** | MEQ is serial. Detected only so that a parallel MFEM fails intelligibly. |

MEQ reads all of these out of MFEM's own `share/mfem/config.mk` rather than
probing, and prints what it found at configure time. It accepts either layout an
MFEM build leaves behind — an installed tree, or an in-source `make` build with
`config/config.mk`.

> **`make clean` after editing any MFEM header.** MFEM's makefiles have no
> header dependency tracking, and a stale object has produced heap corruption in
> unrelated functions and "unimplemented" aborts for methods that had just been
> added. MEQ's own CMake build tracks headers properly; this applies only to the
> MFEM tree.

## Building MEQ

`MFEM_DIR` is a cache variable, also reads the environment, and defaults to
`../mfem/install` relative to the source tree — so on a machine laid out that
way the `-D` is unnecessary.

The rest of the dependencies:

| dependency | required? | notes |
|---|---|---|
| C++17 compiler, CMake ≥ 3.20 | yes | Boost is found in `CONFIG` mode, because CMake 4 removed `FindBoost` |
| toml11 | yes | Vendored as `extern/toml11`. **A clone without `--recursive` fails at configure time**, and the error says so — `git submodule update --init --recursive` fixes it. |
| netcdf-cxx4 | optional | Found with `pkg-config`; Debian calls it `libnetcdf-c++4-dev`. Without it the gridded `.nc` output is dropped from the library and the driver cannot write it. |
| Boost | **yes** | **Boost.Math is header-only and is needed to build the LIBRARY**, not merely the tests: `src/meq/Zernike.cpp` gets its Jacobi polynomials from it, the Zernike radial polynomial being a Jacobi polynomial under a coordinate change. No link dependency is added. |
| Eigen | **yes** | **Header-only, and needed to build the LIBRARY.** `src/meq/SurfaceFit.cpp` takes its least-squares solve from Eigen's `JacobiSVD` — a column-pivoted QR followed by a one-sided Jacobi sweep, which is what that file used to carry by hand. Nothing is linked and no MEQ header mentions it, so a consumer of `meq_core` takes on nothing. Debian calls it `libeigen3-dev`. |
| Boost.Test | for the tests | `unit_test_framework`, the shared-library build |
| clang-tidy | optional | Gates the `naming` test only; skipped with a message if absent |

SuiteSparse, SUNDIALS, GSLIB, LAPACK, MKL and CUDA are all reached *through*
MFEM rather than found separately.

### Running the tests

```sh
cd build && ctest --output-on-failure
```

Nothing needs to be set in the environment: the test registration puts
`MKL_NUM_THREADS=1` on every test itself. Running a test binary **by hand** needs
that same one variable:

```sh
MKL_NUM_THREADS=1 ./tests/SolovievConvergence
```

> **That variable is not cosmetic.** MEQ's inner loop factorises a great many
> *small* dense matrices, and on a threaded MKL each call pays for a thread fork
> and a barrier that dwarfs the arithmetic. Measured on the development machine
> the effect was dramatic and grew rapidly with polynomial degree. Where the
> threshold sits is a property of your BLAS, so measure it on yours; the safe
> setting is 1.

Expect the suite to be **red while a known defect stands** — MEQ's tests assert
the behaviour that is wanted and fail until it is there, never the reverse. A
failing test's message is the record of what is outstanding.

## Without MFEM

`find_package(MFEM)` is deliberately **not** `REQUIRED`. Without it, CMake builds
the MFEM-free half of the library — the configuration parser, the profiles, the
source terms, the boundary shapes — skips every file that includes `mfem.hpp`,
and skips every test that needs a solver, printing a status line for each thing
it left out.

```sh
cmake -B build -DMFEM_DIR=/nonexistent
```

This is not a convenience. It is what makes continuous integration possible at
all, since the branch above cannot be fetched by a hosted runner. **Read a green
CI badge as evidence about the configuration and profile layers and about
nothing else**: the HDG assembly, Newton, the extension technique, the estimator
and every convergence claim are exercised only by a local build.

## Where the MFEM branch should live — unresolved

**Today, an outside user cannot build MEQ's solver at all.** `meq-integration`
is a local merge of four branches from a working tree on one machine. There is
nothing to clone, so caching does not help and neither does a submodule: there
is nothing to point one at.

That is a real limitation rather than an oversight waiting to be tidied, and it
is the reason `find_package(MFEM)` is optional and the reason CI covers what it
covers. The options, none of which has been chosen:

* **Publish the merge**, as a tag or branch on a fork of MFEM, and pin MEQ to
  the commit. Cheapest to consume; the cost is that it must be re-published
  whenever any of the four moves, and a stale pin is worse than none.
* **Upstream the four branches**, which is the real fix and is not MEQ's to do.
* **Ship a build script** that fetches the four and performs the merge, with the
  conflict resolutions encoded. Honest about the situation, but it bakes in a
  topology that belongs to somebody else's tree and has already changed once.
* **A container or binary MFEM install** published alongside MEQ.

Until one is chosen, the practical answer for a new machine is to follow
`CLAUDE.md`'s recipe by hand and check the install afterwards.

---

The [documentation](https://meq.readthedocs.io) covers everything past the
build: `README.md` is the overview, and
[the install chapter](https://meq.readthedocs.io/en/latest/install.html) is the
same material as this file in a form that assumes you already have MFEM.
