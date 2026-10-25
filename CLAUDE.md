# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Do not append session logs, temporary notes, or todo lists to this file

MEQ (the Maryland Equilibrium Solver) computes axisymmetric plasma equilibria by
solving the **Grad–Shafranov** equation with a hybridizable discontinuous
Galerkin (HDG) discretisation built on MFEM.

`refs/Refs.md` indexes the papers behind the numerics, with the doi to fetch each
from. The two marked ✔ there are not background reading — they *are* the method,
and `src/meq` is an implementation of them.

## This file is an index

It carries what a maintainer reads every session: the status, the commands, the
traps as a checklist, the testing stance, the layout. Everything else is in a
companion, **read when working on that subject and not before**:

| | |
|---|---|
| **`CLAUDE_HDGGS.md`** | the equation, the discretisation and its sign conventions (`DarcyForm` holds `−q`; the papers' `τ` carries the wrong sign), post-processing, the NPC port, Newton and the bordered Newton, the linear solves and threading, `MKL_NUM_THREADS`, the solver's own traps, and **the fixture ladder and how to assert a rate** |
| **`CLAUDE_FB.md`** | free boundary, FB-A to FB-7, the inverse solve, and the `freegs4e` benchmark |
| **`CLAUDE_INVERSION.md`** | `ψ(R, z)` to `R(Ψ, l)` — critical points, the tracer, the averages, the disc basis, the `(Ψ, θ)` output, the `q(ψ)`-driven solve |
| **`CLAUDE_FLOW.md`** | toroidal flow, FL-0 to FL-8 — RoPP (136), the species closure, and Maschke & Perrin |
| **`CLAUDE_DRIVER.md`** | `apps/meq.cpp` and `meq-run`, the output files and the band between `Γ_h` and `Γ`, the adaptive loop and the reactive ladder through the driver, what the driver refuses, the two performance keys, the warm start, the configuration parser's traps |
| **`CLAUDE_DEVICE.md`** | running under an `mfem::Device` — every defect found, the instruments, residency, cuDSS, the threaded-assembly promise, the integrator inventory and the batched routes, and **single precision** |
| **`CLAUDE_MFEM.md`** | which MFEM MEQ links and how it is configured, `meq-integration` and how to refresh it, `../../mfem-hdg-dev` and what is open upstream, the old-to-new HDG API mapping, and MFEM's API traps |
| **`CLAUDE_TRAPS.md`** | the long form of the *Traps* checklist below — build, provenance, code hygiene, shell tooling, and reading a red test |
| **`CLAUDE_DECLINED.md`** | what was considered and decided against, with its measurement. **Read it before proposing a route**, since most obvious ones are in it |

**THREE RULES APPLY IN EVERY ONE OF THESE FILES.**

* **`MEASUREMENTS.md` holds the published tables** under stable anchors, and a
  `→ **[M-nn](MEASUREMENTS.md#m-nn)**` line stands wherever one is cited.
  **Do not renumber the anchors.** The argument a table supports stays in the
  companion; the digits live there. Reference tables — output formats, stages,
  MFEM options, the branch list, the fixture ladder — stay inline.
* **A hypothesis that was MEASURED AND FALSIFIED stays**, because recording it
  is what stops it being re-derived.
* **The history of the files does not.** Present tense; no "used to", no
  "an earlier version", no dated "now fixed".

The split is by subject and not by size, so a claim touching two companions is
written in the one that owns the subject and cross-referenced from the other.
Plans live in `*-PLAN.md`; **`ROADMAP.md` is the priority order** and says which
items are MEQ's and which are requests on upstream. **`docs/` is the user-facing
account** (see *Layout* for the rule on what may go there).

## Status

**The solver works.** MEQ solves the semi-linear Grad–Shafranov equation by
Newton at `k+1` in both `ψ` and `q`, against an exact Solov'ev equilibrium on the
linear path and HDG-GS-1's Example 5 on the non-linear one, with Newton
quadratic; `SolovievConvergence` and `NewtonConvergence` are the acceptance.
`src/meq/GradShafranov.{hpp,cpp}` is the solver — `DarcyForm` with
`EnableHybridization`, ordering `NonlinearOrdering::NPC`. **`ψ_ax` is an unknown
of the non-linear system**: `setSource( NormalisedSource &, double )` closes it
by a bordered Newton (`HighBetaConvergence`), and `meq::NormalisedMHDSource` is
the production source, reachable as `[source] Type = "mhd"`, `Normalised = true`.
`v0-legacy` has the pre-port tree.

**MEQ is a program.** `meq config.toml` parses, meshes, solves — adaptively if
asked — and writes the equilibrium to **one file by default, the `.nc`**, plus
an optional reduction of it; exit codes 0–4 are `docs/running.rst`'s. `meq-run config.toml` runs the
`[mesh.generate]` generator first and then execs `meq`; it is a script because
`meq` links MFEM and not gmsh, and it **does not parse the configuration** — it
asks `meq --mesh-command`, so there is one reader of the schema.

| | |
|---|---|
| `.nc` | **always.** `ψ*` and **B** on a uniform `(R, Z)` grid, with `inside` and `extrapolated` masks — **drop `extrapolated` nodes before an error norm or a diff** — and beside it the mesh as MFEM's own text and every coefficient of `ψ_h`, `q` and `ψ*`. The interchange format **and the restart** (`[initialguess] File = "<stem>.nc"`) |
| `.mesh`, `_psi.gf`, `_grad_psi.gf`, `_psistar.gf` | `--gf` / `[output] GridFunctions`. Exact, every coefficient. GLVis, and the `.gf` restart |
| `<stem>/<stem>.pvd` | `--vtk` / `[output] VTK`. VTK at the degree of the field it draws (`k+1`, since `ψ*` is that field), boundary bent onto `Γ`. ParaView |
| `_surfaces.nc` | flux surfaces and averages against `ρ = √Ψ_N`, only with `[output] FluxSurfaces` |

**`meq-convert <stem>.nc` makes the `--gf` and `--vtk` files from the `.nc`,
byte for byte what the run would have written** — `OutputConversion` asserts it
— except an adaptive run's per-cycle series and a split run's `_psi_total.gf`,
which only the run can make (→ **[M-297](MEASUREMENTS.md#m-297)**).

`CLAUDE_DRIVER.md` has all of it: the band and why `B` is continued with its
gradient, the curved and adaptive paths through the driver, the reactive ladder
(**never predictive** — nothing may be inferred from `F` about which solver to
run), the bordered retry rung that is MEQ's one branch control
(→ **[M-247](MEASUREMENTS.md#m-247)**), and the refusals.

**THE TWO PERFORMANCE KEYS ARE THE ONLY SOLVER KNOBS IN THE TOML, AND THE LINE
IS DRAWN ON PURPOSE.** `AssemblyMode` and `TraceSolver` cannot change the answer
where both converge. **`Globalisation` is not exposed and must not
be** (and there is no `NonlinearOrdering` to expose: NPC is the only
ordering): three solve routes reached discrete solutions **9.4% apart** on an
under-resolved mesh, which is exactly where an adaptive run starts
(→ **[M-26](MEASUREMENTS.md#m-26)**; `CLAUDE_HDGGS.md`, *Should
`PicardThenNewton` simply be the default?*). That is a dated measurement: its
third route, condense-first, is not in MEQ, and its two NPC routes, plain and
`PicardThenNewton`, read `max ψ_h` 1.0% apart; no test prints the comparison
on the current Picard map (M-303). A key
that silently changes which equilibrium is reported is not a performance key.

**The device path reproduces the host.** Under `--device cuda` the production
cases match `--device cpu` to every printed digit, three of upstream's four
element-local groups run on the device, and the remaining cost of a device run is
host threading it must give up (M-188, → **[M-223](MEASUREMENTS.md#m-223)**).
The driver's device defaults are deliberate — cuDSS as the trace solver,
`AssemblyMode = "serial"` to keep bit-reproducibility, threaded assembly refused
above one thread, `TraceSolver = "cudss"` refused without a Device.
**Single precision** builds and runs at parity with double, qualified per
problem. `CLAUDE_DEVICE.md` is the record. **This machine's GPU is for
correctness, not timing** — consumer FP64 is 1/32 to 1/64 of FP32.

**The warm start is MEQ's own `meq::PointLocator`**, so `MFEM_USE_GSLIB` is
optional; it carries a stored guess across meshes and each adaptive cycle into
the next, for a third of the Newton work (→ **[M-01](MEASUREMENTS.md#m-01)**).

| stage | | |
|---|---|---|
| 0 | Git reconciliation, tag `v0-legacy` | **done** |
| 1 | Tree, CMake, `Config` / `Profiles` / `Source` | **done** |
| 2 | Linear `Δ*` on `DarcyForm`, fitted polygonal domain | **done** |
| 3 | Local post-processing `ψ*_h` | **dropped as MEQ code** — `DarcyForm::Reconstruct()` supplies it, at `k+2` |
| 4 | Newton on the semi-linear source | **done** |
| 5 | Curved `Γ` by extension from subdomains | **done** |
| 6 | Adaptivity: the residual estimator and mesh update | **done** |

| campaign | | |
|---|---|---|
| toroidal flow | FL-0 to FL-8, **done** | `docs/rotation.rst`, `CLAUDE_FLOW.md` |
| solution inversion | IN-A to IN-P, **done** | `INVERSION-PLAN.md`, `CLAUDE_INVERSION.md` |
| free boundary | FB-A to FB-7 done or answered, FB-6 met, XP-0 to XP-4 done — MEQ reproduces `freegs4e`'s diverted equilibrium | `FREE-BOUNDARY-PLAN.md`, `CLAUDE_FB.md` |
| rotating mirror | MR-0 to MR-5 and §5's midplane gauge done — every free-boundary mirror run reports the fraction of `∫|F|/R` its answer puts beyond `Γ` and warns above 1e-4, exactly zero on both shipped examples (M-310); MEQ solves geq's rotating mirror and agrees with it to geq's own resolution in either gauge (M-261 to M-263, M-278); `examples/mirror-rotating.toml` runs it through the driver (M-264) with the inner label `ψ_inner` at the limiting electrode, `[boundary.electrode]`, started like `ψ_bnd` at the starting state's flux at its point (so `PsiAxis` and `[source.axis]` are refused beside the box), adaptively under a judged estimator (M-279), on the card (M-280), and with no retry rung beside the box (M-277); `examples/mirror-three-species.toml` solves e, D and C⁶⁺ at `M = 6` in the midplane gauge (M-278), quadratic with `R₀` at 32 knots as border unknowns (M-283, M-284); `Gauge = "vacuum-midplane"` pins the potential on the vacuum's crossing, quadratic with an exact Jacobian (M-281), and the `.nc` carries the true-midplane densities (M-282). A box support defaults to `τ = 3/max(R, 0.1 m)` (`meq::MirrorStabilisation`, M-293): the three-species ceiling at constant `τ = 1` was the stabilisation's fold inside one element (M-291), and under the default 2e19 converges on level 0 in both midplane gauges and 5e19 in the midplane gauge, asserted; only the vacuum-midplane gauge has a fold of its own, near 4.16e19 (M-294). Open: a second discrete root reached from `ψ_inner = 0.15` under the default, which only a far `[initialguess]` can now pose; a meshed-pack mirror through the driver in single, whose floor is above the solver's target; `q_h` reference-limited on the rotating mirror's finest rungs; `ψ*` gaining nothing over `ψ_h` on the vacuum mirror; both midplane gauges untried under `--device` and adaptivity | `MIRROR-PLAN.md`, `CLAUDE_FLOW.md` |
| element regions | ER-0 to ER-3 done — each element's material fractions, the material, circuit and permeability curve tables in the `.nc`, VTK cell data and the `node_material` grid mask, classified by the solver's own source predicate (M-266 to M-270) | `ELEMENT-REGIONS-PLAN.md`, `CLAUDE_DRIVER.md` |
| inverse solve | IV-0 to IV-6 done — coil currents are border unknowns, `[inverse]` runs through the driver and exits 4 on a best fit (M-250 to M-253) | `INVERSE-AND-RECONSTRUCTION-PLAN.md`, `CLAUDE_FB.md` |

Each stage ends at a **measured convergence rate**, not at "it runs".

## Where things go

**THIS REPOSITORY IS `source/`, AND NOTHING WRITES INTO IT.** The project
directory holds four things, and only `source/` is tracked:

| | |
|---|---|
| `source/` | this repository — the session's working directory |
| `builds/<name>/` | one CMake build tree per configuration (`double`, `single`, `coverage`); deletable |
| `install/<name>/` | `cmake --install builds/<name>`: `bin/meq`, `bin/meq-run`, `share/meq/tools/`; what runs and tools execute |
| `runs/` | everything a run writes, a tool's under `runs/<tool>/`; deletable |

CMake **refuses** a build tree inside `source/`. Tests run in their build tree
and read the source through `MEQ_SOURCE_DIR` (`tests/SourceTree.hpp`); tools run
the installed binary and write under `runs/` through `tools/meqpaths.py`. **A
path a configuration READS is relative to the TOML file; a path the run WRITES
— `[output] Directory`, `CacheFile`, a `[mesh.generate]` mesh — is relative to
the working directory**, so a run is made from a directory under `runs/`:
`../../install/double/bin/meq ../../source/examples/x.toml`. **Untracked
output in `git status` is a bug in whatever wrote it**; `.gitignore` does not
hide run products. Sibling checkouts (`mfem`, `mfem-hdg-dev`, `freegs4e`,
`geq`) sit beside the project directory, so `../../<name>` from here.

## Commands

From the project directory, `..`:

```sh
git -C source submodule update --init --recursive     # extern/toml11
cmake -S source -B builds/double
cmake --build builds/double -j6
cmake --install builds/double                         # -> install/double
cd builds/double && OMP_NUM_THREADS=4 ctest -j4
```

**`-j4` with `OMP_NUM_THREADS=4`** holds the product at the core count, since
every solver test already threads its assembly → **[M-14](MEASUREMENTS.md#m-14)**.
**`OMP_NUM_THREADS` is deliberately not pinned in `tests/CMakeLists.txt`**, so a
plain `ctest` exercises threaded assembly as it ships. The one red is
`PlasmaEdgeConvergence`, deliberate while the staircase `Γ_{p,h}` stands. The
wall is dominated by `FreeBoundaryCoupling` (forty solves reproducing M-134, kept
whole on purpose) and then `naming` (one single-threaded clang-tidy run).
**A suite time is a measurement only on an idle machine**; identical code has
read 1.16× to 1.6× apart here, so read ratios.

**ctest sets `MKL_NUM_THREADS=1` on every test**; run a binary directly the same
way, from the build tree, since a test writes into its working directory —
`cd builds/double/tests && MKL_NUM_THREADS=1 ./SolovievConvergence`. Without it
the suite takes over an hour (see *Traps*).

`DeviceParity` runs three arms, `cpu`, `debug` and `cuda`, and skips loudly
where no device is configured (→ **[M-228](MEASUREMENTS.md#m-228)**).

The performance harness is not a ctest, since every number in it is a timing:

```sh
tests/performance/scan.sh build 3            # LINEAR path, one thread axis at a time
tests/performance/npc-scan.sh build 3        # NONLINEAR path, OMP = MKL = N together
tests/performance/newton-step-profile.sh build 5   # where a Newton step's time goes, by leg
```

`npc-scan.sh` sweeps both axes together because under threaded assembly MKL
suppresses its own threading inside the OpenMP region (`CLAUDE_HDGGS.md`,
*Threading, measured*).

`MFEM_DIR` defaults to `../../mfem/install` from `source/` and reads the
environment; the single build takes `-DMFEM_DIR=` its own install.

### Coverage and CI

```sh
cmake -S source -B builds/coverage -DMEQ_ENABLE_COVERAGE=ON
cmake --build builds/coverage -j6 && (cd builds/coverage && ctest)
cd source && gcovr --root . ../builds/coverage --filter 'src/meq/' --print-summary
```

`MEQ_ENABLE_COVERAGE` is an option and not a build type so that `-O0` is chosen
rather than inherited. **A debug build revives no assertion**: `MFEM_ASSERT` is
gated on `MFEM_DEBUG`, which the install sets `NO`, and MEQ has no `assert()`
(`CLAUDE_MFEM.md`). `gcovr` is not installed here — use a venv, with
`--merge-mode-functions=separate`, from `source/`.

**CI CANNOT BUILD THE SOLVER, AND THIS IS STRUCTURAL.** `meq-integration` is a
local merge published nowhere, so `find_package(MFEM)` is not `REQUIRED` and CI
builds the MFEM-free half: four unit suites, `naming` and a 90% `gcovr` line
gate. **A green CI badge is evidence about the configuration and profile layers
and nothing else.** The `full` job is `if: false` until the branch is published.

**No coverage percentage is quoted, deliberately.** `gcovr` counts comment-heavy
lines `gcov` does not (92.7% against 64% on one file, one run), so the gate's
number partly measures comment density. Branch coverage is the figure with real
headroom. What a percentage cannot show is uncalled API: `SolverContract.cpp` exists
because 22 of `GradShafranovSolver`'s 51 methods were called by nothing.

### Docs

```sh
.venv-docs/bin/sphinx-build -W -b html docs ../builds/docs/html
make -C docs/manual                       # the LaTeX manual, into ../builds/manual
```

`.venv-docs` carries exactly `docs/requirements.txt` — sphinx 7.4.7 and
sphinx-material 0.0.36 — matching Read the Docs' `fail_on_warning`. **Never use
the system `sphinx-build` (8.2.3, above the pin) and never create a second
venv**; `ls -d .venv*` before concluding a tool is absent.
`sphinx-build -W -b dummy` checks prose and cross-references without a theme.

### Dependencies, and preferring a maintained library

**toml11 is `extern/toml11`**, pinned (`b32a2ff`, v4.4.0-31, as MaNTA pins it),
because its conversion behaviour is version-dependent (see *Traps*). A configure
failing on toml11 is a clone without `--recursive`. **MFEM stays out of tree**,
found through `MFEM_DIR` (`CLAUDE_DECLINED.md`, *MFEM as a submodule*).

**PREFER A MAINTAINED LIBRARY TO A HAND-ROLLED ALGORITHM**, for maintenance
rather than speed: a performance tie is not an argument for bespoke code, and the
bar for a header-only, well-maintained library is that it breaks nothing.
**For device code the order is `mfem` abstraction > CUB > a manual kernel** —
CUB before a kernel, after MFEM, because MFEM's abstractions also dispatch to
HIP. **A contract MEQ has already made outranks the order**: where a sum is
documented bit-identical to a dense loop (`meq::CompressedRows`,
`assemblePlasmaCurrent()`), a reassociating `mfem::reduce` is the wrong
primitive; `meq::DeviceRows` is the worked example
(→ **[M-190](MEASUREMENTS.md#m-190)**). Taken on these grounds: Boost.Math for
the Zernike radials, Eigen for `SurfaceFit` — whose truncation threshold and LM
damping *are* the gauge, so after a swap the ellipse result, the `nstx` tail and
the zero axis spread must not move.

## MFEM → `CLAUDE_MFEM.md`

**MEQ links `../../mfem/install`**: MFEM 4.10.1 built from branch
**`meq-integration`** in `../../mfem/mfem-src`, a local merge of seventeen branches
that exists only to be built against — never pushed, re-created rather than
merged into, tagged before each re-creation. `mfem/master` will not do.
`CLAUDE_MFEM.md` has the options, the CUDA obstacles, the branch table, the
recipe with its conflict resolutions, and the checks.

**After any refresh**: run the containment loop **and** the unfiltered
merge-subject listing (the loop only checks branches it names); check the
install's timestamp against the branch tip **and** diff an installed header
(they fail in different directions); rebuild MFEM with `make -j4` at most, touch
the installed headers, and rebuild MEQ `--clean-first`. A class-layout change
behind `MFEM_THREAD_SAFE` is otherwise a silent ODR violation.

**`../../mfem-hdg-dev` IS RECEIVE-ONLY.** Write MFEM development *requests* into
its `doc/` and nothing else — no branches, commits, checkouts or builds; another
agent works there. Upstream does not track `*-FROM-*.md` in git, so a request
appears as an untracked file by design, and an absent document means nothing.
**The only test of "landed" is the symbol in `../../mfem/install/include`.** Do not
file findings against unfinished work.

## Traps

The checklist; the measurement and the story of each is in the file named.

**Environment and machine**

* **`MKL_NUM_THREADS=1`** on every run of the serial element loop: MFEM's
  element-local dense LU at `k ≥ 3` is 40× slower with MKL threading, which looks
  like a solver problem. Threaded assembly suppresses it, which is why PARDISO
  is the default trace solver (`CLAUDE_HDGGS.md`, M-59, M-78).
* **Never a bare `make -j`.** WSL2 falls over. `-j6` for MEQ, **`make -j4` at
  most for MFEM**, and `make clean` after editing an MFEM header (no dependency
  tracking).
* **The GPU is an RTX 2070 SUPER**: correctness evidence, never a timing.

**Build and provenance** (`CLAUDE_TRAPS.md`, `CLAUDE_MFEM.md`)

* **`cmake --install` keeps source mtimes**, so MEQ looks up to date across an
  MFEM layout change: `find ../../mfem/install/include -name '*.hpp' -exec touch {} +`.
* **A fresh configure re-searches cached results** (`BLAS_LIBRARIES`,
  `MKL_PARDISO_*`, MEQ's own copied link line): copy them from the reference
  cache and diff `MFEM_EXT_LIBS` token by token before comparing two installs.
* **Never `cmake --build` while `ctest` runs** — the result mixes binaries.
* **`meq --version` and the `.nc` carry commit, precision and backends**,
  written at build time by `cmake/GitCommit.cmake`; anything read by
  `execute_process()` is read once, at configure.
* **`#ifdef MFEM_USE_SINGLE` in a translation unit with no MFEM header is
  silently false**; include `mfem/config/config.hpp`.

**Code hygiene no compiler checks** (`CLAUDE_TRAPS.md`)

* **A braceless `for`/`if`/`while` whose next code line is a bare `{` after a
  gap adopts that block** — no flag in g++ or clang warns. Grep for the shape
  after any edit to a control-flow head.
* **A `switch` over a `Config` enum with no `default` after an initialised
  variable keeps the initialiser**, and `-Wswitch` did not fire here. Find a new
  enum value's switches by hand.
* **A mechanical rename is not verified by a build and a lint.** Build a
  separate tree clean with `-Wshadow -Winit-self` and read every warning for the
  renamed name — a reason that makes one site safe is not a reason about the
  category. `python3 -m py_compile` every tracked `.py`. Fenced blocks of
  captured output are data and must not change.
* **A verification that compares printed values names which binaries it ran**,
  and a red test cannot double as a regression for something else. To separate
  *already failing* from *I broke it*: `git worktree add <dir> HEAD --detach`,
  configure with `-DMFEM_DIR=` and `-DTOML11_DIR=`, and diff the failure
  inventories.
* **Never commit a file a subagent owns on the strength of its compiling** —
  scaffolding compiles. Read the diff.

**MFEM's API** (`CLAUDE_MFEM.md`)

* **`Mesh::GetElementTransformation( int )` and
  `GetBdrFaceTransformations( int )` return shared scratch**; threaded, that is a
  silent wrong answer. Use the caller-allocated overloads into a `thread_local`.
* **`Mesh::FindPoints` is not reentrant and is `O(elements × points)`.**
* **`TransformBack()` clamps an outside point onto the element** (its default
  solver is `NewtonElementProject`), so an overshoot read off the reference
  coordinates is zero. Use an unprojected `InverseElementTransformation`.
* **A `SubMesh` keeps a pointer to its parent**; anything returning one must
  outlive the parent.
* **A decorator that overrides `Mult()` and not `ArrayMult()` silently
  un-batches** the solver it wraps (→ **[M-98](MEASUREMENTS.md#m-98)**).
* **`EnsureNodes()` makes an affine mesh look curved** to anything keying on
  `GetNodes() == nullptr`; ask the nodal degree instead.
* **A nonlinear potential mass and a linear one do not mix**, loudly for face
  integrators and silently for domain ones. `DarcyHybridization` freezes the
  element-local initial guess at `FormLinearSystem()`. `Reconstruct()` is right
  where `∂F/∂ψ` vanishes only with `gf-hdg-dev`'s fix in the merge.

**The solve** (`CLAUDE_HDGGS.md`, *Traps in the solve*)

* **A target scaled to `‖r₀‖` makes a good guess fail**; `solve()` takes its
  reference from the cold iterate. `NewtonSolver` needs `iterative_mode = true`,
  or the Dirichlet data is silently discarded.
* **GS-2's §4.2–4.5 sources have a trivial branch at `ψ = 0`** and Newton finds
  it in zero iterations; `setInitialGuess()` is the way out.
* **No residual evaluation between `NPCGradient()` and the solves it
  entitles.**
* **The suite's fixtures and bounds were measured at `τ = 1`**: no constant `τ`
  is best for both `ψ_h` and `q_h`, `τ = 0.3` breaks Newton across the
  non-linear suite, and from `τ = 10` `CriticalPointFinder` finds no zero of
  `q_h`, so a user's `Tau = 10` loses the axis search (M-286). **A box support
  defaults to `τ = 3/max(R, 0.1 m)`, every other run to constant `τ = 1`**;
  beside a box a `Tau` with no `TauScaling` is refused, since it would be
  `τ₀` in m⁻¹ there and a constant in m⁻² elsewhere (M-294).
* **`TauScaling = "diffusion"` is `τ₀/max(R, R_f)`, and `R_f ≤ 0` is
  refused**: an unbounded `τ₀/R` costs `q` half an order on any mesh touching
  `R = 0`, since the face rule makes `τ h` constant there, and gives `τ ~ 1e18`
  at the gmsh half-disc's `R ~ 1e-18` quadrature points, which stalls Newton.
  A tiny `R_f` (1e-6) still brings the half order back (M-292). In single the
  floor of a solve rises with `τ`'s size, so a single-precision bound set at
  `τ = 1` is re-measured, not assumed (M-294).

**Configuration and seams** (`CLAUDE_DRIVER.md`, `CLAUDE_FB.md`)

* **toml11's `find_or<double>` returns the default for an integer node**;
  `Config.cpp` reads numbers through its own `asFloat()`.
* **A reserved key is only reserved on the paths that call the refusal**, so a
  test asserts the refusal on every source type, not on one.
* **A bound on how wrong an estimate may be says nothing about an estimate of a
  different quantity**; `meq::limiterPolygonMaximum()` is the single
  implementation of that constraint.
* **A term added to an extrapolating evaluation must be defined wherever the
  extrapolation goes, or the evaluation must be able to refuse** — and which,
  at a seam, is decided by parity (→ **[M-156](MEASUREMENTS.md#m-156)**).

**Shell tooling** (`CLAUDE_TRAPS.md`)

* `cd` persists between Bash calls; shell **variables** do not, and `grep -q`
  with an empty filename reads stdin forever.
* `pgrep -f`/`pkill -f` match the invoking shell; `-x` silently matches nothing
  past 15 characters, which most test binaries are. **Wait on an artefact, with a
  bound**, match the marker the writer actually wrote, unanchored, and kill a
  run's waiters when the run is restarted. `grep -c` exits 1 on zero matches.

## Testing stance

**A test asserts the behaviour that is wanted, and fails until it is there.**
Never the reverse: a test asserting a known defect passes while the defect
stands, and **100% green must mean there is no lurking defect.** So a defect
gets a *failing* test naming it, whose message says what fixes it. The suite is
expected to be red while a known defect stands.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` is the worked
example: written to assert `ψ*` is a post-processing of `ψ_h` everywhere, red
until MFEM fixed the reconstruction, then green on its own — which is what put
the driver back on `ψ*`.

**A wrong sign convention converges, at the right rate, to the wrong
function.** Only a closed form catches it, which is why stage 3 is pinned to
`refs/HDG-GradShafranov-Adaptive.pdf` §4.1's exact NSTX Solov'ev solution and
not to a self-convergence study. The ladder:

1. **Unit**, `tests/unit/`, Boost.Test — config parsing, splines and `dFdPsi`
   against closed forms, and the *assembled* Jacobian against a central
   difference of the assembled residual.
2. **An exact solution**, as an absolute-error regression.
3. **Convergence rates**, `tests/convergence/` — a *rate*, not a tolerance:
   `k+1` for `ψ` and `∇ψ`, `k+2` for `ψ*`, over `k = 1…4` and dyadic meshes.
   Rates flatten past `k ≈ 5` or `h/8`, which is round-off and not a bug.

The analytic fixtures (`Soloviev`, `McCarthy`, `ManufacturedNonlinear`,
`SimilarityExponential` — a ladder in `∂F/∂ψ`), the six fixed-boundary machines,
the two-tier rate assertion on unfitted meshes, asserting Newton on its best
triple, and why a rectangle cannot show `k+1` are in `CLAUDE_HDGGS.md`, *The
fixture ladder, and how to assert a rate*. **Mutation-test a suite you rely
on.**

**WHEN A CASE GOES RED** (`CLAUDE_TRAPS.md`, *Reading a red test*):

* **Suspect the fixture before the solver.** Twice a long-red case was a fixture
  missing something production has (coils and a prescribed current; calling
  `setPlasmaSupport()` on the wrapper). **A one-key experiment separates two
  hypotheses only if everything else is where you think it is** — it took a
  `2 × 2` cross both times (→ **[M-134](MEASUREMENTS.md#m-134)**,
  **[M-82](MEASUREMENTS.md#m-82)**).
* **Check it is not a coin toss.** `MKL_NUM_THREADS` is a numerically neutral
  perturbation; apply it to **both** binaries before calling a red a regression
  (→ **[M-186](MEASUREMENTS.md#m-186)**). ctest's pin at 1 is one draw.
* **Relax a precondition, never a tolerance.** A case on a knife edge sweeps
  several meshes or weights, skips one that does not solve, and requires a
  majority — measured to be reachable before it is asserted.
* **A test that cannot fail is worse than none** — e.g. `checkAxis().agrees`
  under the located-axis constraint, which is 1 by construction.
* **When a plan's stages compose, the acceptance that matters is the one on the
  composition** (→ **[M-148](MEASUREMENTS.md#m-148)**).

## Layout

```
src/meq/     the library: Config, Profiles, Source, SourceFactory,
             GradShafranov, BoundaryShape, Estimator, Field, Sampler,
             WarmStart, Output, StoredEquilibrium (the .nc read back as
             mesh and fields: the restart and meq-convert's input),
             Materials (the element material table and
             its record, MFEM-free); CriticalPoints, FluxSurfaces, SurfaceAverage,
             SurfaceFit, Zernike (the inversion); ConductorStore (psi_c on
             disk, MFEM-free, a CACHE at the solver's own points -- CS-6);
             DeviceBasis, AxisSearch, ResidentIterate (device routes);
             Crossings, StageChecksum (device instruments, not part of the
             solve). See CLAUDE_DEVICE.md for the device files.
apps/        meq.cpp, meq-convert.cpp (the GLVis and ParaView files from a
             .nc), and meq-run.in -- configured twice: into the build
             directory as `meq-run`, and for install/<name>/bin, paths
             relative to itself. Not runnable from apps/.
tests/       SourceTree.hpp (how a test names a source file), unit/ (Boost.Test), convergence/ (rate assertions), analytic/
             (closed forms), performance/ (timings, none a ctest:
             TraceSolverScaling, NpcThreadScaling, NewtonStepProfile,
             InversionScaling).
tools/       meqpaths.py (the layout, for every tool), plot_equilibrium.py (both NetCDF files; tools/README.md says which
             format goes with which reader), mesh/halfdisc.py (FB-6's gmsh
             geometry; --symmetric meshes z >= 0 and reflects it for
             [solver] UpDownSymmetry, --symmetry-check asks it of any mesh,
             M-127), and three independent-code benchmarks:
             freegs4e-benchmark/ (free boundary); desc-benchmark/ (DESC, fixed
             boundary, races examples/fixed-*.toml; profiles are labelled by
             toroidal flux, so each case is posed self-consistently UNTIMED
             and only the forward cold solve is clocked, M-157);
             chease-benchmark/ (CHEASE, same equation and free functions;
             converge.py refines each code separately, as a MATRIX, because
             only a difference that falls rules out a conversion error, M-158).
examples/    TOML run configurations.
docs/        the Sphinx manual for Read the Docs; docs/manual/ is the older
             LaTeX one.
refs/        Refs.md is tracked; the PDFs are gitignored, fetch by doi.
```

**SEVEN FILES IN `src/meq` ARE COMPILED BY A DEVICE COMPILER** —
`*Device.cpp`, holding every `mfem::forall` in the tree, because `mfem::forall`
`static_assert`s under `g++`. Without CUDA they are host stubs and the build is
host-only and correct, which is what CI builds.

**`Profiles` and `Source` are MFEM-free on purpose** — plain `double`
arguments — so they are unit-testable without the library; the
`mfem::Coefficient` adapters live with the assembly.

**WHAT MAY GO IN `docs/`: documentation reflects the state of the CODE, not a
historical record.** A page says what MEQ does and, where a choice is exposed as
an option, which way it came out and that it was measured — then tells the
reader to measure it themselves. The **numbers** stay in `MEASUREMENTS.md`,
beside the tests that produce them, because a measurement in a manual goes stale
silently. A live trap a user meets ("freegs4e's docstring says normalised flux
and is wrong") belongs in both.

**Style** follows MaNTA: **tabs**, Allman braces, C++17.

| | |
|---|---|
| Types — class, struct, enum, alias | `UpperCamelCase` |
| Enum values | `UpperCamelCase` |
| Functions and methods | `lowerCamelCase` |
| Variables, parameters, members | `lowerCamelCase` |
| **TOML configuration keys** | **`UpperCamelCase`** — a deliberate mismatch; do not reconcile |

**The project is MEQ in prose and `meq` in code** — `meq::`, `src/meq/`,
`meq_core`, `meq config.toml`. **External names keep their author's
capitalisation** (`cuDSS`, with a `NOLINT` saying so). Enforced by
`.clang-tidy`'s `readability-identifier-naming` through the `naming` ctest over
every file in `src/meq` — except, on a CUDA build only, the seven device files,
because clang-tidy cannot read nvcc's flags; on a host build they are checked
like everything else.

## Git

Branch **`main`**; remote **`origin` = `github:ianabel/meq.git`** (**public**;
`github` is an `~/.ssh/config` alias), so a push or a GitHub release publishes.
**Tag `v0-legacy` is the pre-modernisation tree** — reach for it before
concluding something was lost. **History is not rewritten.** Commit messages end
with the `Co-Authored-By` trailer.

**The version is set once, in `project(meq VERSION …)` in `CMakeLists.txt`**;
`meq --version`, the outputs' `meq_version` and `docs/conf.py`'s `release` all
read it. A release bumps that line, commits, and tags `vX.Y.Z` annotated, with a
GitHub release on the tag.
