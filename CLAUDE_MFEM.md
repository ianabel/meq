# MFEM: which build, the integration branch, upstream, and MFEM's API traps

**`CLAUDE.md`'s MFEM half** — the install MEQ links and its options, turning CUDA on, `meq-integration` and how to refresh it, the receive-only rule for `../../mfem-hdg-dev` and what is open there, the old-to-new HDG API mapping, and the traps in MFEM's API and build that have cost MEQ time.

**The standing rules are `CLAUDE.md`'s.** Measurement tables live in `MEASUREMENTS.md` under stable `M-nn` anchors, which are not renumbered. A hypothesis that was measured and falsified stays; the history of the file does not — present tense, no "used to", no dated "now fixed". What was decided against is `CLAUDE_DECLINED.md`.
## Which MFEM, and why not master

**MEQ builds against `../../mfem/install`** — MFEM **4.10.1** on branch
**`meq-integration`**, CMake-built in `../../mfem/build` from sources in
`../../mfem/mfem-src`. The options actually set, read out of
`share/mfem/config.mk` and re-checked 2026-09-01 rather than remembered:

| | |
|---|---|
| `MFEM_USE_SUNDIALS` | KINSol, so `KINSolver(KIN_LINESEARCH)` is reachable. **Now `../../sundials/cuda-install`, not `../../sundials/install`** — see the CUDA row |
| `MFEM_USE_GSLIB` | `FindPointsGSLIB`; gslib v1.0.9 built alongside at `../../mfem/gslib` |
| `MFEM_USE_SUITESPARSE` | UMFPACK, the direct solver MEQ's own solver runs on. **`../../suitesparse/install`, not Debian's** — v7.12.2, the same version Debian ships, built against oneAPI MKL. Debian's carried a `NEEDED` on `libblas.so.3` -> `libmkl_rt.so`; see `CLAUDE_HDGGS.md`, *PARDISO and the MKL link line* |
| `MFEM_USE_MKL_PARDISO` | `PardisoSolver`, oneAPI MKL 2026.1, **on the threaded `mkl_gnu_thread` layer**. Resolving it to Debian's `mkl_sequential` instead is a bigger deal than it sounds — see *Rebuilding MFEM* below |
| `MFEM_USE_EXCEPTIONS` | `MFEM_ERROR_THROW` by default, which is what makes the driver's exit code 2 reachable |
| `MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE` | **both, and both load bearing.** `DarcyHybridization::SetAssemblyMode( Threaded )` needs both and **aborts** rather than falling back without them. MEQ reaches it through `setAssemblyMode()` |
| `MFEM_USE_LAPACK`, `MFEM_USE_ZLIB`, `MFEM_USE_DOUBLE`, `MFEM_USE_MEMALLOC` | |
| `MFEM_USE_MPI = NO` | MEQ is serial throughout; `../../mfem-hdg-dev-par` exists for parallel work |
| **`MFEM_USE_CUDA = YES`**, `CUDA_ARCH = sm_75` | CUDA 13.3. **`sm_75` because the card is compute capability 7.5**; MFEM's cache defaulted to `sm_60`, which is wrong here and nothing warns you |
| **`MFEM_USE_CUDSS = YES`** | cuDSS 0.8.0 from Debian (`/usr/include/cudss.h`), found with `-DCUDSS_DIR=/usr`. `mfem::CuDSSSolver` agrees with UMFPACK to 3.5e-14; its **timings on this machine are not reproducible** and MEQ's solver does not use it. See `CLAUDE_HDGGS.md`, *Threading, measured* |
| `RAJA`, `OCCA`, `CEED`, `SIMD`, `AMGX` all `NO` | |

**Turning CUDA on is not one flag, and two of the three obstacles are silent.**

* **SUNDIALS must itself be a CUDA build.** `linalg/sundials.hpp` hard-errors
  with `MFEM_USE_CUDA=TRUE requires SUNDIALS to be built with CUDA support`,
  testing `SUNDIALS_NVECTOR_CUDA`. MEQ cannot drop SUNDIALS — `KINSolver` backs
  `AndersonPicard` and `PicardThenNewton` — so the answer is a second SUNDIALS.
  `../../sundials/cuda-install` is 7.5.0 configured `-DENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=75`, otherwise matching the old one.
* **`MFEM_ENABLE_TESTING` is ON, and it needed an upstream fix to be.** A
  serial CUDA build would not generate at all: the bug was that a device test loop merged the serial and parallel example lists and
  never re-split them, registering a test against `sundials_ex9p` in a build
  that never created it. Both loops — `examples/CMakeLists.txt` and
  `examples/sundials/CMakeLists.txt` — now guard the parallel branch with
  `elseif (MFEM_USE_MPI)`, so a serial build skips it. **Verified**:
  reconfigure, full rebuild and install all returned 0, and ctest registers
  **131 tests**. `MFEM_ENABLE_TESTING` does not appear in the installed
  `config.mk`, so it changes nothing MEQ's finder reads.

  **These notes claimed the bug was "filed as
  `../../mfem-hdg-dev/doc/CMAKE-SERIAL-CUDA-SUNDIALS.md`". IT WAS NOT.** That path
  has never existed in that tree's history and exists nowhere on disk — the
  report was never written. Treat a "filed as" claim in these notes as needing a
  check, not as evidence.
* **Stale derived cache variables survive a `-DSUNDIALS_DIR=` change.**
  `SUNDIALS_INCLUDE_DIR` and friends are separate cache entries, so repointing
  the directory leaves the old include path on `TPL_INCLUDE_DIRS` *in front of
  or behind* the new one, and which wins is an ordering accident. Delete every
  `SUNDIALS_*` entry but `SUNDIALS_DIR` and reconfigure.
* **AND THE SAME CACHE SURVIVES A NEW SUNDIALS *COMPONENT*, WHICH IS WORSE
  BECAUSE THE SYMPTOM NAMES NOTHING.** `meq-integration` was
  re-created with a `sundials-ida-integration` branch merged in, which adds
  `mfem::IDASolver` and adds `IDAS` to `SUNDIALS_COMPONENTS` in MFEM's
  `CMakeLists.txt`. `SUNDIALS_LIBRARIES` is a cached STRING, so
  `find_package(SUNDIALS)` short-circuited on the old list and
  `-lsundials_idas` never reached the link line — while `libmfem.a` now carries
  `sundials.cpp.o` referencing `IDASVtolerances`, `IDAPrintAllStats` and eight
  more. **Every MEQ test binary failed to link**, with an error naming IDA and
  nothing naming a cache. MFEM's source was correct throughout. Same cure:
  delete every `SUNDIALS_*` entry but `SUNDIALS_DIR`, reconfigure, reinstall.
  **Not filed upstream**, per the standing rule about findings against work that
  has just landed — and there was nothing to file, since the defect was entirely
  in MEQ's own build directory.

**And a CUDA build broke MEQ's own `FindMFEM.cmake`, which is worth knowing
because every serial build hides it.** MFEM puts the toolkit headers on the line
as `-isystem <dir>` inside `MFEM_CXXFLAGS`, not as `-I<dir>` in
`MFEM_TPLFLAGS`, while `general/backends.hpp` includes `<cuda_runtime.h>`
unconditionally once `MFEM_USE_CUDA` is defined. MEQ's finder took only
`TPLFLAGS`, so the first `#include <mfem.hpp>` failed on a missing CUDA header —
a confusing place to meet a build-configuration problem. The finder now also
parses `-isystem` out of `MFEM_CXXFLAGS`.

**WHAT CUDA BUYS IS THREE THINGS AND `fem/darcy/`'s MISSING PARTIAL ASSEMBLY IS
NOT ONE OF THEM.** Partial assembly is still absent and still irrelevant — its
economic case is sum factorisation and MEQ is triangles. What the flag enables
is `BatchedLinAlg`'s `gpu_blas` backend and `CuDSSSolver`; **MFEM's batched
element and face kernels**, which are reached through `AssemblyMode`,
`LocalFactorMode` and `TraceAssemblyMode` and **are** the device route; and
**MEQ's own device translation units** — `DeviceNPCDevice.cpp`,
`DeviceQuadratureDevice.cpp`, `DeviceRowsDevice.cpp`, `PointLocatorDevice.cpp`,
`SourceViewDevice.cpp`, `AxisSearchDevice.cpp` and `ResidentIterateDevice.cpp`,
which CMake hands to `nvcc` at the architecture it reads out of `config.mk`, and
which hold every `mfem::forall` launch in `src/meq`.
Without CUDA those seven compile to host stubs and the build is host-only and
correct, which is what CI builds. → **[M-223](MEASUREMENTS.md#m-223)**.

**AND IT IS NOT FREE — IT COSTS 7% OF THE WALL CLOCK ON THE CPU PATH, MEASURED
AGAINST AN OTHERWISE IDENTICAL INSTALL.** `../../mfem/install-nocuda` differs from
`../../mfem/install` in `MFEM_USE_CUDA` and the `MFEM_USE_CUDSS` that must follow
it, and in nothing else; the two produce **byte-identical output** on the machine
case, so this is cost and not answers.

→ **[M-77](MEASUREMENTS.md#m-77)** — wall · user · allocation count

The mechanism is visible and only partly the obvious one: `mfem::forall` on the
host path builds an `__nv_hdl_wrapper_t`, nvcc's host-lambda wrapper, whose
`manager::do_call` shows in `perf` at 1.9–2.9% — and interposing `operator new`
at its call site counts **226.7 million allocations in a 30-second run**, of
which **55.4 million (24%) go away with CUDA off**. The other 171 million are
MFEM's ordinary per-element temporaries and are not CUDA's doing. So the flag is
a prerequisite that currently charges rent, and a CPU-only production build has
a measured reason to be a separate install.

**`../../mfem-hdg-dev` is now the development tree and not what MEQ links.** That
split exists because the alternative was demonstrated: mid-session that tree was
switched to `gf-hdg-dev`, which has no `fem/darcy/extension_hdg.{hpp,cpp}`, and
MEQ's `naming` check failed on a tree MEQ had not touched. A library MEQ depends
on should not move under it while somebody is working on it. Point `MFEM_DIR` at
the dev tree deliberately when testing a fix; otherwise leave it alone.

## `meq-integration`: the branch MEQ builds from, and why it is local only

MEQ needs work from **seventeen branches**, all fetched from `hdgdev`. Several
are fixes branched off MFEM master rather than the HDG line, and two are merged
under their local names rather than as `hdgdev/<name>`. See the warning under the
containment loop for how branches missing from this table were found and why no
check built from the table could find them:

| branch | what MEQ needs from it | |
|---|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — stage 5's curved `Γ`, and `TransferredDatumCoefficient`. **AND IT IS THE TARGET FOR ANYTHING MEQ SENDS BACK ABOUT THE EXTENSION MACHINERY**, free boundary's transmission quadrature included: it carries `TransferPath`, `ElementExtension`, `PathTraceCoefficient` and `ExtensionRegionQuadrature`, and carries **no NPC at all** — verified 2026-09-05 — so a patch against it provably does not depend on the ordering | merge base |
| `direct-solver-symbolic-reuse` | `UMFPackSolver` and `PardisoSolver` keeping their symbolic factorisation across Newton steps | merge 1 |
| `gf-hdg-linearise-first` | **`DarcyNPCOperator` / `DarcyNPCSolver`** — the NPC method, which is what MEQ's default ordering now is — plus `SetAssemblyMode`, `SetGradientMode`, `SetLocalFactorMode`. **`SetNonlinearOrdering` was on this list and is DELETED**; see `CLAUDE_HDGGS.md`, *The NPC port* | merge 2 |
| `gf-hdg-dev` | ***"The postprocessing closes on the element average, always"*** — the reconstruction fix, without which `ψ*` is a different function wherever `∂F/∂ψ` vanishes | **merge 3, again** |
| **`pardiso-multi-rhs`** | **`PardisoSolver::ArrayMult`** — the multi-RHS trace solve. MEQ's bordered step depends on it for **1.34×**, → **[M-98](MEASUREMENTS.md#m-98)**, and **losing it does not break the build**: `Operator::ArrayMult`'s base loops `Mult()`, so a recipe that drops this branch compiles, answers identically and is silently slower | merge 4 |
| **`pardiso-device-host-sync`** | ***"PardisoSolver: synchronise the vectors, not just the matrix"*** — the other half of M-79's device chain, on the default trace solver | merge 5 |
| **`sundials-ida-integration`** | `mfem::IDASolver`, and the `IDAS` component that broke every MEQ test binary's link when the cache went stale — see the CUDA section's fourth bullet | merge 6 |
| **`gf-hdg-p-adaptivity`** | upstream's own §2 work; MEQ calls nothing from it, and it is in this table because it is in the merge, which is the point of the table | merge 7 |
| **`memory-manager-registration-refcount`** | **Off MFEM master, not the HDG line**, though `hdgdev` carries it, which is why a recipe built from the HDG list cannot see it. `internal::Memory` gains the reference count `internal::Alias` always had, so registering one address twice no longer lets the first holder's `Delete` evict the second's entry. Four routines register lazily and reach it: `MakeAlias` on an unregistered base, the device-class arms of `Read`/`Write`/`ReadWrite`, and `SetDeviceMemoryType`. MEQ proposed this shape after both projects' first mechanism proved wrong; see the `Reconstruct()` row above | merge 9 |
| **`cmake-find-package-component-cache`** | the fix to the defect MEQ reported when `SUNDIALS_LIBRARIES` short-circuited on a stale cached list and `-lsundials_idas` never reached the link line — see the CUDA section's fourth bullet. **It is in the merge because MEQ asked for it**, and its absence is a link failure rather than a silent one | merge 10 |
| **`cudss-solver-configuration`** | `CuDSSSolver`'s setters, including the one that refuses `CUDSS_MTYPE_SYMMETRIC` by name rather than returning an exactly zero solution every phase reports as success. MEQ does not call cuDSS from a config file and still wants this: `TraceSolverScaling` constructs one | merge 11 |
| **`intrules-triangle-high-order`** | `c6693d5f09`, the simplex-rule ceiling. Triangle rules are tabulated to **order 25** and above that `GrundmannMollerSimplexRule()` returns partly NEGATIVE weights — correct in exact arithmetic, destructive in double, and returned with a one-time warning rather than refused. **MEQ is clear and it was checked rather than assumed**: `SourceIntegrator::rule()` asks for `2k + OrderW + 4`, which is 12 at the suite's `k = 4`, and `PlasmaEdgeConvergence`'s triangle rules are `2k + 14` over `k <= 3`, so 20. `setSourceQuadratureOrder()`'s doxygen already stated the bound | merge 12 |
| **`gf-interp-hdg-dev`** | **CCSZ INTERPOLATORY HDG** — `fem/darcy/reaction_hdg.{hpp,cpp}` with `HDGInterpolatoryReactionIntegrator`, `NodalReactionFunction` and a quadrature control, plus `HDGPostprocessBlocks` in `postprocess_hdg.*` and `DarcyHybridization::Bg_data`. **15,024 insertions across `fem/darcy/`**, verified in `libmfem.a` by `nm`. MEQ does not call any of it yet; `INTERPOLATORY-HDG-PLAN.md` is the plan | merge 8 |
| **`sundials-real-t`** | **`linalg/sundials.{hpp,cpp}` in `real_t` rather than `double`**, plus a `static_assert` that `real_t` and `sunrealtype` agree and a makefile that stops refusing SUNDIALS in a single-precision build. **MEQ reported it and MEQ needs it**: without it MFEM does not compile at all in the `MFEM_PRECISION = single` + `MFEM_USE_SUNDIALS` pairing, which is `install-single`, and KINSOL is what backs `AndersonPicard` and `PicardThenNewton`. MEQ's own six sites were the ones that stop the compile; upstream found **73**. It is a branch off `origin/master` rather than one of the HDG line, merged as `hdgdev/sundials-real-t` | merge 13 |
| **`pardiso-refinement-control`** | `PardisoSolver::SetMaxRefinementSteps()`, `GetNumRefinementSteps()` and `GetNumPerturbedPivots()` — `iparm[7]` settable, `iparm[6]` and `iparm[13]` readable. MEQ asked for it ([M-300](MEASUREMENTS.md#m-300)): `setTraceRefinementCap()` and `StepProfile`'s refinement counts call it, so **its absence is a compile error**. The same refresh carries `DarcyHybridization::GetNPCReduceTime()`, `GetNPCRecoverTime()` and `GetNumElementColours()`, on `gf-hdg-linearise-first` and `gf-interp-hdg-dev` | merge 14 |
| **`transformback-single-tolerance`** | `d10e8c48e7`: `TransformBack()`'s default `tol_0` is `1e-6` in single, the value `meq::inverseTransformTolerance()` passes at MEQ's own sites, so the default and MEQ's call are the same arithmetic ([M-271](MEASUREMENTS.md#m-271); `CLAUDE_DEVICE.md`, *What was not a tolerance*). MEQ passes its own tolerance, so its absence moves nothing MEQ calls. Merged as `transformback-single-tolerance`, not `hdgdev/…`, so the `hdgdev/`-anchored half of the pattern below cannot see it | merge 15 |
| **`densematrixsvd-ldvt`** | `cee6a378af`: `DenseMatrixSVD` with `jobvt = 'S'` hands LAPACK an `LDVT` equal to `Vt`'s height, where it wrote past its own buffer for a wide matrix. MEQ reported it ([M-249](MEASUREMENTS.md#m-249)) and calls nothing from it, the inverse corner being Eigen's `JacobiSVD`. Merged under its local name, like merge 15 | merge 16 |

**THE TOPOLOGY HAS NOW CHANGED THREE TIMES AND THIS ROW HAS BEEN WRONG TWICE.**
It has read "pairwise independent", then "an ANCESTOR of both, so it needs no
merge of its own" — verified 2026-08-30 and true then — and as of **2026-09-05
`gf-hdg-dev` is once more ahead of the merge and needs a merge of its own**, by
a single commit (*"Four Darcy test files the trunk's CMake build never ran"*),
which conflicts in `tests/unit/CMakeLists.txt` alone. **The standing instruction
is therefore to run the containment loop and believe it over this table**, which
is a claim about someone else's branch and goes stale without MEQ doing
anything. The reconstruction fix itself has been contained throughout; a "not
contained" here has always meant *the branch has moved*, never *the fix is
missing*.

What has not changed is *why* `gf-hdg-dev` is in the list. Dropping it gives a
silently degraded `ψ*` on any element where the Jacobian's reaction term
vanishes, an adaptive loop a full order worse, and a suite that stays green
until `thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` runs. It is
free today rather than unnecessary.

`meq-integration` is their merge, and it exists **only** to be built against. It
is not a development branch and must not be pushed anywhere: it will be
re-created whenever any of them moves, so anything committed directly to it is
lost. Do the work on the branch it belongs to and re-merge.

**It lives in `../../mfem/mfem-src`, which is MEQ's own clone — never in
`../../mfem-hdg-dev`.** See the rule below.

**Expect the fetch to be a forced update.** Two of the four were rebased between
2026-08-29 and 08-30, so a merge *into* the existing `meq-integration` is not
the operation you want; re-create it. Tag the old one first — the tree carries
`meq-integration-before-2026-08-30` for exactly that reason — because a bad
re-merge is otherwise unrecoverable without redoing the conflict resolution.

Refreshing it, from `../../mfem/mfem-src`:

```sh
git fetch hdgdev                                 # every branch below, not a subset
git tag -f meq-integration-before-$(date +%F) meq-integration
git checkout -B meq-integration hdgdev/gf-hdg-subdomains-dev
git merge hdgdev/direct-solver-symbolic-reuse     # 1; 3 conflicts, see below
git merge hdgdev/gf-hdg-linearise-first           # 2; 3 conflicts, see below
git merge hdgdev/gf-hdg-dev                       # 3; 1 conflict as of 2026-09-05
git merge hdgdev/pardiso-multi-rhs                # 4; M-98 depends on this
git merge hdgdev/pardiso-device-host-sync         # 5
git merge hdgdev/sundials-ida-integration         # 6
git merge hdgdev/gf-hdg-p-adaptivity              # 7
git merge hdgdev/gf-interp-hdg-dev                # 8; the interpolatory half
git merge hdgdev/memory-manager-registration-refcount   # 9; off master
git merge hdgdev/cmake-find-package-component-cache     # 10; else a link failure
git merge hdgdev/cudss-solver-configuration       # 11; TraceSolverScaling
git merge hdgdev/intrules-triangle-high-order     # 12
git merge hdgdev/sundials-real-t                  # 13; single precision needs it
git merge hdgdev/pardiso-refinement-control       # 14; MEQ's setTraceRefinementCap()
git merge hdgdev/transformback-single-tolerance   # 15
git merge hdgdev/densematrixsvd-ldvt              # 16
# ALWAYS finish with this, and believe it over the table above:
for b in gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
         gf-hdg-linearise-first gf-hdg-dev pardiso-multi-rhs \
         pardiso-device-host-sync sundials-ida-integration \
         gf-hdg-p-adaptivity gf-interp-hdg-dev \
         memory-manager-registration-refcount \
         cmake-find-package-component-cache cudss-solver-configuration \
         intrules-triangle-high-order sundials-real-t \
         pardiso-refinement-control transformback-single-tolerance \
         densematrixsvd-ldvt; do
  git merge-base --is-ancestor hdgdev/$b HEAD && echo "$b contained" \
                                              || echo "$b NOT CONTAINED"
done
```

**Run that last loop.** It is seventeen cheap commands and it is the only thing that
catches a branch having moved out from under the recipe — which is precisely
what happened to the three-branch version of this section.

**AND IT CANNOT CATCH A BRANCH THE RECIPE NEVER NAMED, WHICH IS THE WORSE
FAILURE AND IS THE ONE THAT HAPPENED.** The loop asks *is each branch I listed
contained*, so it is only ever as complete as the list — and the list was
missing **five of nine**. `gf-interp-hdg-dev` alone is 15,024 lines of
`fem/darcy/`, the whole interpolatory-HDG half; `pardiso-multi-rhs` is what
M-98's 1.34× rests on. Re-creating `meq-integration` from the recipe as it stood
would have deleted all five, and **four of them without a compile error**:
`Operator::ArrayMult` has a working base implementation, and MEQ calls nothing
in `reaction_hdg.hpp` yet. The build stays green, the answers do not move, and
the solver gets slower.

**Ask the tree the other question instead.** It needs no list, which is the
whole point. It is what found `memory-manager-registration-refcount`, and it
finds `transformback-single-tolerance` and `densematrixsvd-ldvt` (merges 15 and
16), which are merged under their local names and which an `hdgdev/`-anchored
pattern misses. It also shows one merge of `origin/master` itself
(`463593cbb3`, master at `45799ebea8`), which is the trunk and no branch of the
table, and **the recipe needs no step for it**: it brought two files MEQ never
builds, `.github/workflows/builds-and-tests.yml` and
`miniapps/nurbs/makefile`, and that master commit is an ancestor of merges 11,
13, 14, 15 and 16, so re-creating the branch from the recipe contains it.
`git merge-base --is-ancestor 45799ebea8 meq-integration` is the check if the
recipe changes. Match every form, over the merges this branch made rather
than MFEM's whole history:

```sh
git -C ../../mfem/mfem-src log --first-parent --merges --format=%s \
    hdgdev/gf-hdg-subdomains-dev..meq-integration \
  | grep -oE "hdgdev/[a-z0-9-]+|branch '[a-z0-9/-]+'|^Merge [a-z0-9-]+( into|:)" \
  | sed -E "s/^branch '(.*)'$/\1/; s/^Merge ([a-z0-9-]+)( into|:)$/\1/; s|^hdgdev/||" \
  | sort -u
```

**The pattern is as much a list as the loop was**, which is the lesson repeating
one level down: an `hdgdev/`-anchored grep encodes the assumption that every
merge comes from that remote, and the first branch that does not is invisible to
it. Widen the pattern or read the merge subjects unfiltered.

**BUT THE LOOP CHECKS THE BRANCH, NOT THE INSTALL, AND THAT GAP HAS ALREADY
BITTEN.** Every command in it runs in `../../mfem/mfem-src` and says nothing
whatever about `../../mfem/install`, which is what MEQ actually links. On
2026-09-01 the loop reported **all four contained** while the installed library
was `fa65a2f932` and `meq-integration` had moved to `d244329d8c` — **514 changed
lines in `fem/darcy/` across seven files**, including *"The HDG face quadrature
never saw the trace element"* and *"DarcyOperator dereferenced a null
prolongation on a hanging-node-free NC mesh"*, the second in exactly the
non-conforming meshes the adaptive loop makes. So the one check this project
trusts passed on a tree whose measurements were all taken against a library a
day out of date.

**A REFRESH IS A CLEAN REBUILD OF BOTH TREES, AND THE REASON IS THE CLASS
LAYOUT.** Anything that moves a member behind `#ifndef MFEM_THREAD_SAFE` — which
MEQ's build defines, and MEQ constructs such objects itself — changes an object's
size without changing a signature, so an incremental link is an ODR violation
with no diagnostic. **MEQ must be fully rebuilt too**, `cmake --build ../builds/double -j6
--clean-first` and the same for every other tree under `../builds`, and the absolute timing tables are historical afterwards:
[M-138](MEASUREMENTS.md#m-138)'s budget is the one that wants re-taking, since
it is absolute where [M-137](MEASUREMENTS.md#m-137)'s payoff is a ratio between
two MEQ binaries on one library and survives.

**BOTH INSTALLS ARE CURRENT.** Both checks below pass on `../../mfem/install` and
on `../../mfem/install-single` against `meq-integration` at `6e1cc66618`:

* the batched `Reconstruct()`, `ConstructC()` and assembly predicates;
* setup's threaded loops;
* from `9e5ec65b50`, `ComputeH()` leaving the batched factors on the
  device ([M-320](MEASUREMENTS.md#m-320), [M-321](MEASUREMENTS.md#m-321));
* from `6ef93b0e28`, `HDGElementGeometry`, `HDGStabilizationRule` and
  the two weight predicates ([M-322](MEASUREMENTS.md#m-322));
* and, from `9d57974336`, `GPUBlasBatchedLinAlg`'s persistent pointer arrays,
  `HDGElementPointsAndWeights()` under `EnsureElementQuadratureTables()`, and
  the extension's bisection to round-off with `SetLevelSetGradient()`
  ([M-324](MEASUREMENTS.md#m-324)). Three path classes gained a member, so
  that refresh was a layout change;
* and, from `6e1cc66618`, a `SundialsNVector` view under CUDA returning the
  memory-manager registration it takes. Every KINSOL callback left one count
  behind, so a freed vector's entry outlived its buffer and a later, larger
  allocation at the same address inherited the smaller device buffer:
  `NewtonConvergencePaths_cuda`'s illegal address in double and its "invalid
  alias" in single, 0 runs in 30 in each build since. Two `.cpp` files, so
  an incremental rebuild.

The DIII-D machine reads `psi_ax = 3.759850e-01` and a constraint of
`-6.466e-12` under `--device cuda`, and `DeviceParity` passes in both builds.

* **Double passes 88 of 89.** The red is the deliberate
  `PlasmaEdgeConvergence`.
* **Single passes 87 of 89**, adding `DriverAcceptance`'s adaptive cases
  (MKL 2026.2, [M-303](MEASUREMENTS.md#m-303)), which fail on `cb96cf3` too. Every
  binary takes
MKL through the RUNPATH `mkl/latest/lib`, which moved from 2026.1 to 2026.2 on
2026-10-20 and moved single's knife edges ([M-302](MEASUREMENTS.md#m-302)).
`LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2026.1/lib` overrides it, which is how
an MKL is ruled in or out. **So an MKL upgrade is a library refresh too.** The single install is a
**second** tree with the same obligation and its own reason to move first:
`sundials-real-t` is in the merge because MFEM does not compile **at all** in
the `MFEM_PRECISION = single` + `MFEM_USE_SUNDIALS` pairing without it.

**AND REFRESHING IT IS NOT FREE: THE LAST REFRESH TOOK A CASE RED BY SEVEN
ORDERS, AND THE ANSWER IT MEASURES HAD NOT MOVED AT ALL**
→ **[M-189](MEASUREMENTS.md#m-189)**. `../../mfem/install` rebuilt from
`meq-integration` at `6d05189467` takes `FreeBoundaryCoupling` red at
`theTransmissionRowIsTheBoundaryIntegralItClaims`, which has read 7.8e-16 to
1.4e-13 since it was written and now reads **1.27e-06**. Bisected to one default
argument: `e38906d853` — upstream's fix to MEQ's own M-183 — changes
`ExtensionBoundaryQuadrature`'s `fd_step` default from the literal `1e-6` to
`ExtensionFDStep()`, which is `cbrt( eps )` and therefore **6.06e-06 in
double**. **MEQ's four production sweeps pass their own step and were never on
that default**; what took it was the reference the TEST builds, so the case
contracted a row at `1e-6` against a reference at `6.06e-6`. Give both halves
the same step and it returns to 1.4e-16.

**THE LESSON IS ABOUT THE TEST AND NOT ABOUT THE LIBRARY.** A central
difference's error is a deterministic function of its step, so it is common to
two sweeps that share one and **cancels exactly** — which is the only way that
case could ever have read 7.8e-16, orders below what a difference at `1e-6` can
resolve. **Where a caller computes one quantity two ways and compares them, both
sides must be given the same step**, and `meq::extensionFdStep()` is now
`inline` in `GradShafranov.hpp` so they can share one expression rather than two
constants kept in step by hand — the `meq::limiterPolygonMaximum()` precedent,
M-165. And MEQ's first report of this blamed the production Jacobian: it was one
`grep` from being checked and was not checked.

**Check the install too, and it is one command:**

```sh
diff <( git -C ../../mfem/mfem-src show meq-integration:fem/darcy/darcyhybridization.hpp ) \
     ../../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp \
  && echo "install is current"
```

**AND THAT CHECK IS BLIND TO A `.cpp`-ONLY FIX, WHICH IS MOST OF THEM.** It
diffs an installed HEADER, and only headers are installed — so a library that
is a bug fix behind reports *current*. Measured: `meq-integration` moved from
`93952626e2` to `dfb5525a35` carrying `d720f02d62`, the fix to MEQ's own
M-214, and every installed header in `fem/darcy/` was byte-identical across it.
The check passed on a library missing the fix MEQ had asked for.

**Ask the TIMESTAMPS instead, which cannot miss a translation unit:**

```sh
test ../../mfem/install/lib/libmfem.a -nt \
     <( git -C ../../mfem/mfem-src log -1 --format=%cd --date=iso meq-integration >/dev/null; \
        echo ) 2>/dev/null
# or, plainly:
[ "$( date -r ../../mfem/install/lib/libmfem.a +%s )" -ge \
  "$( git -C ../../mfem/mfem-src log -1 --format=%ct meq-integration )" ] \
  && echo "install is newer than the branch tip" || echo "INSTALL IS STALE"
```

The header diff still earns its place — it catches a *layout* change, which is
the case that silently produces an ODR violation rather than a missing fix — so
run both. They fail in different directions.

A branch that another agent re-creates is a branch that moves without MEQ doing
anything, so this belongs beside the containment loop rather than in a rebuild
checklist nobody reads.

**And it has already fired again. Run 2026-09-01 against `meq-integration` at
`fa65a2f932`: three of the four report contained and `gf-hdg-linearise-first`
does not**, because that branch has advanced **17 commits** since the merge.
**Read a "not contained" as *the branch has moved*, not as *the merge failed*** —
what MEQ needs from it is in the tree (`DarcyNPCOperator` is in the installed
`fem/darcy/darcyhybridization.hpp`, and NPC is what MEQ runs on). It says a
re-merge is available, not that the build is wrong. Re-create rather than merge
into the existing branch, per the recipe above.

**The second merge conflicts, and WHICH THREE FILES IT IS HAS ALREADY CHANGED
ONCE.** These notes recorded `gf-hdg-linearise-first` as merging clean, then as
bringing `HDG-DEFECTS-FROM-MEQ.md`, `HDG-ROADMAP.md` and the makefile. As of
**2026-09-05 the defects report no longer conflicts** — it is gone from both
sides — and `tests/unit/CMakeLists.txt` has taken its place:

| file | resolution |
|---|---|
| `doc/HDG-ROADMAP.md` | **Not a disjoint pair.** Four regions, two of them large *divergent rewrites of the same prose*, so keeping both sides duplicates whole sections and produces a worse file than either. Take `--ours`, the subdomains side: it is the merge base, it is the newer of the two, and each dev branch keeps its own copy of that roadmap anyway. MEQ neither builds nor reads it |
| `miniapps/hdg/makefile` | a union: `extension` from one side, `navierstokes` from the other. Merge the `SEQ_MINIAPPS` and `PAR_MINIAPPS` lines into one, keep the `ifeq (MFEM_USE_SUNDIALS,NO)` filter, and keep both test targets |
| `tests/unit/CMakeLists.txt` | a union of a sorted source list — `test_darcy_extension.cpp` from subdomains, `test_darcy_npc.cpp` from linearise-first. **`test_darcy_system.cpp` appears INSIDE one hunk and again in the common region below it**, so a mechanical keep-both duplicates it and CMake then compiles it twice. Take the union and assert no entry repeats |

Only the makefile and the CMakeLists affect a build. **The `HDG-ROADMAP.md`
marker count is not scriptable**; count them each time, for the `CHANGELOG`
reason below.

**AND THE SAME UNION TRAP APPEARS AGAIN IN MERGE 3.** `gf-hdg-dev` conflicts in
`tests/unit/CMakeLists.txt` alone and by then the merged file already carries
the union, so `--ours` is right there — but verify it rather than assuming, by
checking every darcy entry from all four parents is present exactly once.

**The first merge conflicts and the resolution is always the same.** The
symbolic-reuse branch is based on a much later upstream (767 commits past where
the subdomains branch left), so `CHANGELOG`, `doc/CodeDocumentation.dox` and
`fem/nonlininteg.cpp` collide where both sides appended. In every case the two
sides are **disjoint** — `nonlininteg.cpp` gets `SumNLFIntegrator` and
`SumBlockNLFIntegrator` from one side and the Navier–Stokes PA kernels
`ConvectiveVectorConvectionNLFIntegrator` and
`SkewSymmetricVectorConvectionNLFIntegrator` from the other — so **keep both
sides** of all three. Check the function names are still disjoint before assuming
that; it is true today and is a property of the branches rather than a rule.

**And "keep both sides" is not quite enough for `nonlininteg.cpp`.** Git leaves
the file's trailing

```
}

}
```

as *common* content after the `>>>>>>>` marker — one brace closing
`SumBlockNLFIntegrator::~SumBlockNLFIntegrator`, one closing `namespace mfem`.
That is correct for whichever single side you keep, and one brace short when you
keep both: the first side's destructor is then left open and the second side's
functions sit at namespace depth 2. **Add the destructor's closing brace at the
seam.** It fails to compile rather than failing silently, but the error appears
forty lines later on the *other* side's first function, which is a misleading
place to start looking. `git show <branch>:fem/nonlininteg.cpp | tail` on both
parents shows what the ending should be, and counting braces across the three
versions finds it in seconds — both parents balance and a bad merge does not.

*A trap when resolving `CHANGELOG` programmatically*: it contains long `======`
heading underlines, so a script that strips lines *beginning* with `=======`
destroys the file. Match the conflict marker exactly, alone on its line.

Then rebuild:

```sh
cd /home/ian/projects/mfem && cmake --build build -j4 && cmake --install build
```

## `../../mfem-hdg-dev` is receive-only

**Write MFEM development *requests* into `../../mfem-hdg-dev/doc/` and nothing
else.** No branches, no commits, no checkouts, no builds. Another agent works in
that tree and it is theirs; MEQ consumes it through `git fetch` from
`../../mfem/mfem-src` and in no other way.

This is written down because it was got wrong: `meq-integration` was first
created *in* the dev tree, given two merge commits and a fix commit, and that
tree was left checked out on it — so the next person to touch it would have been
on someone else's branch. Fetching from the dev tree is reading; everything else
is not.

**WHAT IS OPEN WITH UPSTREAM, AND NOTHING ELSE ABOUT IT.**

| | |
|---|---|
| ~~**The form build is one host core under a Device**~~ | **CLOSED** at `06d6680729` (→ **[M-327](MEASUREMENTS.md#m-327)**): the local blocks are zeroed and kept on the card, `Finalize()`'s backup copied there, the extension's crossing search is Brent's, and `AssembleFluxMassBdrMatricesBatched()`'s face loop threads on `SetIntegratorsThreadSafe()`, which now reaches MEQ's level set and its gradient from several threads, so **a level set must stay reentrant**. The face scatter's `cudaFree` was the kernel's own time. `ff3`'s re-assembly under CUDA 1.20 s to 0.7 s here. Left, not asked: 278 MB of blocks still downloaded for host readers |
| ~~**`R` and the face and element weights, evaluated per point on the host inside the batched routines**~~ | **CLOSED** at `6ef93b0e28` (→ **[M-322](MEASUREMENTS.md#m-322)**): weights and `R` in kernels, with `HDGStabilizationRule` describing MEQ's two stabilisations and `R` read from a `GridFunctionCoefficient` of `meq::CoordinateField`; `HDGFaceWeightsOnDevice()` and `HDGElementMassOnDevice()` say yes on `ff3` and the mirror. `HDGElementGeometry` serves the batched routes without giving the mesh nodes, which closed the mirror's batched regression |
| ~~**Nothing on one host core under a Device: `Reconstruct()`, `ConstructC()`, the batched packing, setup's loops**~~ | **CLOSED** at `4a87f9280f` (→ **[M-321](MEASUREMENTS.md#m-321)**): a batched `Reconstruct()` (both halves, the enriched trace written by each face's highest-numbered element, the serial replay kept as the exact mode); a batched `ConstructC()`, which MEQ takes by `SetDefaultAssemblyMode(Batched)` around `EnableHybridization()` because C is built inside it; the batched assembly's predicates answered once per mesh; and `MarkLevelSetSubdomain()` and `VertexConePath` threaded on a caller's `phi_thread_safe`, effective only with no Device configured. MEQ asks `CanBatchConstraintAssembly()`, `CanBatchReconstructTotalFlux()` and `DarcyForm::CanBatchReconstruction()` and prints the answers under `--profile` |
| ~~**A per-element correction to the local gradient before it is factored, on the batched `NPCGradient()` route**~~ | **CLOSED.** Upstream built `DarcyHybridization::AddNextGradientCorrection( int el, DenseMatrix const &dD )`, with `ClearNextGradientCorrections()` and `NumNextGradientCorrections()`, at `9d57561f69`, in `meq-integration` at `1a15a28363`. `dD` is added to element `el`'s potential-potential block after assembly and before the local factorisation, by the next `NPCGradient()` only, on the element loop and the batched route, on host and device. It aborts under `LocalOpType::FluxNL`, on a non-NPC `GetGradient()` with corrections queued, and on a size mismatch. MEQ's border fold (`GradShafranovSolver::setBorderFolding()`, [M-303](MEASUREMENTS.md#m-303)) is queued there and no longer rides on `SourceIntegrator::AssembleElementGrad()`, so it folds under a Device too. In single on the card the confined control loses the shelf: the folded solve meets its target in four steps where the unfolded one lands 7.5 times over its bound ([M-316](MEASUREMENTS.md#m-316)). `PlasmaConnectivityFoldOn_debug` and `_cuda` are the ctests |
| **Auxiliary globally-coupled unknowns** | `SetNumAuxiliaryUnknowns()` and the two per-element assemble hooks. **The one unbuilt piece, and no longer worth asking for**: the entry points MEQ named are public already as `NPCReduce()` and `NPCRecover()`, and `DarcyNPCSolver::ArrayMult` applies them to several right-hand sides in one pass, so a differenced border is `K` applications of a routine that blocks them. **MEQ CALLS IT** — the bordered step queues every column and flushes once, worth **1.34×** on the DIII-D solve leg at 14 columns, → **[M-98](MEASUREMENTS.md#m-98)** |
| ~~**`DarcyNPCSolver::ArrayMult()`'s per-column staging**~~ | **CLOSED.** MEQ measured it as the NPC leg less upstream's two halves, 32 ms a call at every thread count ([M-301](MEASUREMENTS.md#m-301)). Upstream timed it (`DarcyNPCSolver::GetStagingTime()`, which `--profile` prints) and replaced the copies with views of the caller's columns, at `f3f4fe1ea7`. It is now about 7 ms a call, the three timers close the leg, and the leg at eight threads on TCV is 44% shorter. The views give the same answer on a real card ([M-304](MEASUREMENTS.md#m-304)). The refinement cap is MEQ's `--trace-refinement-cap`, 0 by default, and the colour-free reduce was not wanted. M-304's faster recover comes from first-touch page faults in the old fresh output columns, which upstream measured. So the two halves now scale about alike, and if the leg matters again, both should be measured before choosing one. For the reduce, what remains to try is fewer colours and fewer barriers |
| **`SetIntegratorsThreadSafe()`'s SCOPE, which is the durable half and is now upstream's doxygen** | **CLOSED.** MEQ reported that `mfem::NormalTraceJumpIntegrator` carries `Vector face_shape, normal, shape1_n, shape2_n` and two `DenseMatrix` unguarded (`fem/bilininteg.hpp:4872`), that MEQ installs exactly one as the object handed to `EnableHybridization()`, and that `DarcyForm::ReconstructFluxAndPot()` takes it per face inside its element loop. **Guarded upstream, on the trunk**, so all five branches take it. The class layout moves in a `MFEM_THREAD_SAFE` build, so it needs a clean MFEM rebuild, and the `f3f4fe1ea7` installs are that rebuild, verified (*Both installs are current*); the installed `bilininteg.hpp` carries the guard. **MEQ's report overstated one thing and the correction is the tense**: it said the class is reached by "the loop you are threading", and there is no such loop — `darcyform.cpp`'s three `#pragma omp parallel` regions are all assembly and both reconstruction entry points are serial. MEQ inferred that from being told the work was in progress and wrote the inference as an observation. What upstream took as the real defect was **theirs and not the class**: `SetIntegratorsThreadSafe()` made a promise whose scope was written down nowhere, so threading one more loop silently widens what a caller who set it months earlier has undertaken. Its doxygen now says the promise is about **every integrator installed on the hybridization** rather than about whichever loops happen to be threaded |
| ~~**`mfem::HDGExtensionIntegrator` carries unguarded scratch**~~ | **CLOSED, and MEQ's own files went on asserting it afterwards.** Upstream's `26fccc9f3d` put `shape, shape_ext, nor, x, xbar, m, y, CTm`, `Cmat`, `L` and `ext` behind `#ifndef MFEM_THREAD_SAFE` (`extension_hdg.hpp:24`), and `git merge-base --is-ancestor` confirms it is in the `meq-integration` MEQ links. The row that stood here, and the matching comment at MEQ's `SetIntegratorsThreadSafe()` call, both described the window between the eight-class commit and that one — **a true conclusion resting on a premise that had expired**. Neither was re-read against the tree. *A claim about somebody else's code is only as fresh as the last time it was checked* |
| ~~**`Reconstruct()` aborts under `Device("debug")`**~~ — **CLOSED, and neither project's mechanism was right** | Fixed at `86b0b0fd13`, **contained and installed**, and `DeviceParity` is green on it. The abort is `CheckHostMemoryType_` at `mem_manager.cpp:1778`, from `Vector b_ze( b_z.GetData() + e*nd_ut, nd_ut )` with `e == 0` aliasing `b_z`'s base pointer exactly. **MEQ's mechanism was wrong and so was upstream's first one, in the same way.** Both said `Wrap( ptr, n, own = false )` takes `GetHostMemoryType()` so `std_delete` is false and `Delete_` evicts — but `Delete_` opens `if (!mm.exists \|\| !registered) { return; }`, so an unregistered wrap is **inert** and `h_mt` decides nothing. What actually registers the view is **`Memory<T>::MakeAlias()`**, whose "Register 'base'" branch fires whenever `IsDeviceMemory( GetDeviceMemoryType() )` — true under `debug` and `cuda` — setting `Registered \| OWNS_INTERNAL` on the **view**; `MemoryManager::Insert()`'s `emplace` then keeps the existing entry **silently**, and the view's destructor takes the `Known` branch and erases the **owner's**. In `ReconstructTotalFlux` the trigger is `b_zi.MakeRef( b_ze, ... )` and needs no device flag. **Both of MEQ's candidate fixes fail**, measured: keying on `flags & Registered` is a no-op because `registered` IS true by then, and making `Wrap` take `HOST` when `!own` aborts on a different assertion. What landed is MEQ's **third** option, the one it did not recommend — build the sub-views with `MakeRef` into the owner rather than wrapping a raw `GetData()` — which also syncs, where a raw `GetData()` is a host read of a buffer whose live copy may be on the device. **Scope: 2 of 5 branches**, `gf-hdg-linearise-first` and `gf-interp-hdg-dev`, both of which `meq-integration` merges. **"Why `cuda` never met this" has now been wrong three times and this is the version with evidence**: not that the destructor is inert (`h_mt` decides nothing), and not that the eviction "erases a key that is not there" (`Erase` does `if (mem_map_iter == maps->memories.end()) { mfem_error( "Unknown pointer!" ); }`, so a missing key **aborts**). Under `cuda`, `GetHostMemoryType() == HOST` leaves the **owner** unregistered, so the view's own lazy registration is what **creates** the entry and erasing it on destruction is **balanced**. The immunity is a property of the host memory type, and it is the same observation that makes a reference count the right shape: one registration one erase is already correct, and it is two registrations of one address where today's first erase takes the second's entry. **A second question is live and is upstream's user's**: whether to close the trap in `Memory` rather than per caller. MEQ said yes, in `../../mfem-hdg-dev/doc/HDG-MEMORY-OWNERSHIP-FROM-MEQ.md`, and its addendum argues the change is **smaller than either project first thought** — `maps->aliases` is already reference-counted against exactly this hazard (`InsertAlias` does `alias.counter++` on a collision, `EraseAlias` does `if (--alias.counter) return;`) and `maps->memories` is not (`Insert` discards the `emplace` result, and the `auto res =` is itself inside `#ifdef MFEM_DEBUG`; `Erase` erases unconditionally). The lazy registrations the counter must be correct over are **four routes across five call sites** — `MakeAlias`'s "Register 'base'", `SetDeviceMemoryType()`, and `Read`/`Write`/`ReadWrite`; the two `Wrap` sites need a caller ownership error first. So a counter on `internal::Memory` closes it **without touching ownership semantics or any call site**, where the proposed `own_internal` flag cannot: the call site cannot know whether the address was already registered, so no constant value of it is right. Upstream verified the shape on their trunk, confirmed `mem_manager.cpp` is unmodified from MFEM master in that family, and **withdrew the flag**. **The counter's own cost is bounded rather than hypothetical**: an unbalanced register/erase already leaks a `memories` entry today — `Memory<T>::Reset()` sets `flags = 0` without going near `Delete_`, and `Vector::StealData` is `*p = data; data.Reset(); size = 0;` on a `Memory` that may be registered, with `linalg/hypre.cpp` and `linalg/petsc.cpp` doing the same — so a counter **extends** that class rather than creating it, and turns an unbalanced pair from *an entry erased out from under its owner* into *an entry that stays*. Strictly milder and still silent |

**AND THE BATCHED LOCAL RESIDUAL THAT TAKES A CALLER'S INTEGRAND IS BUILT ON
BOTH SIDES AND CLOSED** → **[M-197](MEASUREMENTS.md#m-197)**,
**[M-223](MEASUREMENTS.md#m-223)**. `machine-f-diiid` under `--device cuda`
reports `batched nonlinear residual yes, gradient yes, device integrand yes`,
so **all eight `CanBatch*` predicates** read yes on MEQ's production bordered
configuration, the whole NPC traversal included. **A default device run prints
`face assembly NO` and that is the KEY and not the predicate**: the driver
leaves `AssemblyMode` Serial under a Device to keep bit-reproducibility
([M-215](MEASUREMENTS.md#m-215)), and asking for `AssemblyMode = "batched"` on
the same file reads `face assembly yes` — with the constraint moving from
`-6.465e-12` to `-6.466e-12`, which is exactly the last-bit difference that
default is protecting. *A predicate printed under a mode nobody asked for
answers about the mode.* What it was asked for
is below, because the SHAPE of the seam is the transferable part.
`HDG-DEVICE-RESIDENCY-FROM-MEQ.md`, request 3.
`CanBatchLinearResidual()` refuses MEQ on two clauses — `m_nlfi_p` not being
linear, and `lop_type == LocalOpType::PotNL` — and **both refusals are
correct**, so what MEQ asked for is a *sibling* rather than a loosening: a batched
residual for the `PotNL` case in which the element-local nonlinearity is
supplied by the caller. The shape proposed is a **header template over a
callable**, instantiated in MEQ's own device translation unit, so that upstream
owns the element loop and the block layout while MEQ owns only the integrand
and neither header names the other's types — the same seam
`meq::PointLocator` already uses one level down. MEQ's half is built and tested
→ **[M-177](MEASUREMENTS.md#m-177)**. **It was the last route MEQ's
configuration needed**: upstream's own audit put MEQ at YES on the other seven
`CanBatch*` predicates, and with it all eight read yes (above).

**THE SEAM, AS SENT AND AS BUILT.** Upstream's reply put it on an interface
rather than on their request 1 — *"it is not blocked on request 1; request 1 is
done. It is blocked on an interface, and the interface is ours"* — and asked
MEQ for the one thing only MEQ could supply: *"Send us the POD and the
signature you want to write `F` and `F'` in and we will build the other side."*
What MEQ sent, `../../mfem-hdg-dev/doc/HDG-BATCHED-RESIDUAL-SEAM-FROM-MEQ.md`,
is what upstream built the other side against: `meq::SourceView`
and its two POD companions, the `MFEM_HOST_DEVICE` evaluators, the
element-local `sourceResidualAtPoint()` / `sourceGradientAtPoint()` pair, and the
four constraints that are not decoration — `radius` unpre-divided because `F/R`
is `μ₀ j_φ` and MEQ's meshes reach `R = 0` exactly, the `psi_c` shift, the
per-ELEMENT connectivity flag, and the view being a **per-call** argument because
`psi_ax`, `psi_bnd` and the current scale move between residual evaluations
inside one Newton step.

**AND ONE HAS CLOSED, WHICH IS THE ROW THAT USED TO SIT HERE.** Threading
`NPCReduce()` and `NPCRecover()` — MEQ's largest measured leg, 0.212 s at one
thread and 0.212 s at eight, **30.6% of a threaded step** →
**[M-126](MEASUREMENTS.md#m-126)** — is **done upstream**, at 4.35× on eight
threads, bit for bit. Two things about how it closed are worth keeping.

**The reason it had been declined was a share measured on the wrong problem.**
Upstream's doxygen declined both loops at *"under 6% of the step"*, taken on
four **fixed-boundary, single-right-hand-side** cases. A bordered Newton applies
`J⁻¹` to `N + 4` columns, so the traversal is `O( elements × columns )` where
the integrator loops are `O( elements )` — the ratio inverts, and the same two
legs that are 6% there are 30.6% here. *A share is a property of a problem, not
of a loop.*

**AND IT IS NOW IN THE LIBRARY MEQ LINKS, AND MEASURED** →
**[M-135](MEASUREMENTS.md#m-135)**. `../../mfem/install` is rebuilt from
`meq-integration` at `3424f33bd5`, the install is verified against the branch,
and an interleaved pair on a quiet machine — a pre-upgrade binary kept aside,
`libmfem.a` being static, so the library is the only variable — reads the leg at
**1.315 s at 1.00 cores against 0.503 s at 3.32 cores, 2.61×**, with every other
leg flat within 2% and the **whole run 1.25×**. Upstream predicted 1.24× from
MEQ's own arithmetic. `psi_ax` is `3.759851e-01` in all twelve runs.

**THE PREDICTION LANDED FOR A REASON WORTH KEEPING: THE SHARE WAS RIGHT EVEN
THOUGH THE FACTOR WAS NOT.** The leg factor here is 2.61× where upstream
measured 4.35× on their own fixture — MEQ's traversal is twelve calls of
fourteen columns and theirs is not — and the *run* number still came out where
the arithmetic said, because what that arithmetic rested on was M-126's
**share**. A leg speedup quoted from somebody else's problem does not transfer;
a share measured on yours does.

**AND IT IS A 5% LOSS AT ONE THREAD**, which a speedup table hides: 1.375 →
1.446 s on the leg and 6.983 → 7.121 s on the run at `OMP = MKL = 1`, the
OpenMP region and the colouring being paid for and unusable. Upstream's own
table has the same sign at the same size. A serial deployment of MEQ is
very slightly worse off for this change.

**MEQ'S OWN HALF OF THE THREADING IS ALSO MEASURED**, and its rules are
`CLAUDE_HDGGS.md`'s *MEQ's own threaded loops*. The three items — the sparse border rows and the five plasma
element loops and the critical-point sweep — are worth **`1.16×` on the whole
DIII-D run at eight threads**, against the `1.14×` predicted, and
`1.03×` at one thread where only the sparsity item has anything to give
→ **[M-137](MEASUREMENTS.md#m-137)**. The `before` arm is the pre-threading
commit carrying **one** cherry-picked line, `SetIntegratorsThreadSafe()`,
without which it exits 2 against today's library — both arms on the shipping
`libmfem.a`, so MEQ's own source is the only variable.

**AND THE BUDGET EVERYTHING WAS SIZED AGAINST HAS INVERTED** →
**[M-138](MEASUREMENTS.md#m-138)**. About **59% of a MEQ run now threads**
against the ~30% M-126 measured, and the largest single-threaded item was
`DarcyForm::Reconstruct()` at 18.5% — MFEM's, unmoved by MEQ's own threading,
and ahead of the NPC traversal it used to sit behind. **Upstream has since
threaded it**, worth `1.177×` on the whole run
→ **[M-172](MEASUREMENTS.md#m-172)**, so the largest single-threaded item MEQ
pays for is `re-assembly` at about 6%, which MEQ owns: the line search's
right-hand side, linear in the exterior coefficients, which a trial now forms in
place (`CLAUDE_FB.md`, FB-5, → **[M-317](MEASUREMENTS.md#m-317)**). **Size anything on this path
against M-172, then M-138, and not against M-126.**

Everything else MEQ has sent is closed. **A CLOSED REPORT NEEDS NO ENTRY HERE**:
what it changed is in the code with a test on it, or it is a measurement under an
`M-nn` anchor, and either way these notes are the wrong place to re-tell it. Three
operational facts are worth the space and the rest is not:

* **UPSTREAM DOES NOT TRACK INTER-PROJECT CORRESPONDENCE UNDER `doc/` AT ALL, SO
  AN ABSENT DOCUMENT MEANS NOTHING.** `46ac3d3bc7` deletes every `*-FROM-MEQ.md`
  in one commit, on the principle that working notes exchanged with a consumer
  are not the library's documentation and that what they establish belongs in
  doxygen on the thing it is about. **That is a policy and not a verdict** —
  reading an empty `doc/` as "every request closed" is the obvious inference and
  is wrong. What a request established survives in upstream's own plan documents
  and in the doxygen; the correspondence does not.
* **SO THE ONLY TEST OF "LANDED" IS THE CODE MEQ BUILDS AGAINST.** Look for the
  symbol in `../../mfem/install/include`, not for the document. That was always the
  only test MEQ could apply and it is now the only one that exists.
* **Which documents exist is a question for `git`, not for `ls`** — a listing
  reports whatever branch that tree is checked out on, which is not MEQ's to
  control and has changed under these notes more than once. **And a request MEQ
  writes will never appear in `git` at all**: that tree's `.gitignore` carries
  `doc/*-FROM-*.md`, so a report lands as an untracked working-tree file by
  design and is read there. `git status` not showing it is the policy working,
  not the write having failed — check with `ls`, which is the one question `ls`
  answers better.

**And do not file findings against unfinished work.** A branch that exists is not
a branch that is done. Measure it if it is useful to know, keep the numbers in
MEQ, and wait to be told it is ready before writing anything into that tree about
how it behaves.

Rebuilding the install, after fetching whatever is wanted into `../../mfem/mfem-src`:

```sh
cd /home/ian/projects/mfem
cmake --build build -j4 && cmake --install build
```

The finder reads `share/mfem/config.mk` from an install and `config/config.mk`
from an in-source make build, so either layout works and neither needs a flag.

`mfem/master` will not do. MEQ needs `DarcyForm` and the HDG integrators in
`fem/darcy/`, and for stage 5 it needs the transfer-path machinery in
`fem/darcy/extension_hdg.{hpp,cpp}` — which is an implementation of exactly the
curved-boundary technique the GS papers use.

## The old HDG API, and how it maps onto `DarcyForm`

**The old MFEM's HDG API is gone, and that is why this is a port and not a
recompile.** `HDGBilinearForm`, `HDGDomainIntegratorGS`, `HDGFaceIntegratorGS`,
`AssembleSC` and `Reconstruct` have no counterparts; the replacement is
`DarcyForm` with `EnableHybridization`. The mapping, term for term:

| Weak term | Old | New |
|---|---|---|
| `(R q_h, v)` | `HDGDomainIntegratorGS` | flux mass: `VectorMassIntegrator` with an `R` coefficient |
| `(ψ_h, ∇̄·v)`, `(q_h, ∇̄w)` | `HDGDomainIntegratorGS` | flux divergence: `VectorDivergenceIntegrator` |
| `⟨ψ̂_h, v·n⟩` | `HDGFaceIntegratorGS` | the transpose of `NormalTraceJumpIntegrator` — **not** the face integrators on `B`, see below |
| `⟨τ(ψ_h − ψ̂_h), w⟩` | `HDGFaceIntegratorGS` | potential mass: `HDGDiffusionIntegrator` |
| `⟨q̂_h·n, μ⟩ = 0` | `AssembleSC` | `EnableHybridization(trace_space, new NormalTraceJumpIntegrator(), ess_list)` |
| condense / reconstruct | `AssembleSC` + `Reconstruct` | `FormLinearSystem` / `RecoverFEMSolution` |
| `φ_h` on `Γ_h` (stage 5) | — | `HDGExtensionIntegrator` on the **flux mass form**, `C = R`, sign `+1` |

`miniapps/hdg/convdiff.cpp` in that tree is the worked example to copy from;
lines 445–469 build exactly the three spaces above.

**A trap in that table, measured.** The `⟨ψ̂_h, v·n⟩` coupling does **not** come
from `TransposeIntegrator(DGNormalTraceIntegrator)` added to `B` on interior and
boundary faces, which is the natural reading. Under hybridization
`DarcyForm::Assemble()` builds `B` through `ComputeElementMatrix()`, which sums
**domain integrators only**, so those face integrators are never assembled —
changing the boundary one's coefficient from `−2` to `−0.7`, and deleting the
interior one outright, does not move a single digit of the answer. The coupling
comes from the transpose of `NormalTraceJumpIntegrator`, supplied to
`EnableHybridization`.

The boundary-face integrator on `B` still has to be there, but as a **marker**:
`EnableHybridization()` reads `B`'s boundary-face markers to decide where to
register the flux constraint. Remove it and the Dirichlet faces get no
constraint at all — the error goes flat at 1.5e-1. So it is load-bearing for a
reason unrelated to the integral it appears to compute, which is exactly the kind
of thing to leave a note about.

**And `DarcyForm::GetOffsets()` returns three entries, not four** — it never
learns about the trace space. The miniapps get the 4-entry version from
`DarcyOperator::ConstructOffsets()`, which is not in `libmfem.a`, so
`GradShafranovSolver` builds its own.

## `MFEM_ASSERT` and `NDEBUG`

**AND THE TEMPTING REASON IS WRONG IN A WAY THAT COSTS SOMEBODY AN AFTERNOON.**
Routing `-O0` through `CMAKE_BUILD_TYPE` does **not** drop `NDEBUG` in a way that
changes which `MFEM_ASSERT`s are live. Audited:

* `MFEM_ASSERT` is gated on **`MFEM_DEBUG`**, never on `NDEBUG`
  (`general/error.hpp`), and nothing anywhere in the installed tree derives one
  from the other.
* The installed MFEM has **`MFEM_DEBUG = NO`**, so **every `MFEM_ASSERT` is
  already dead in every MEQ build**, release or debug.
* MEQ's own sources contain **no `assert()` and no `<cassert>`**, so `-DNDEBUG`
  disables nothing of MEQ's either. It is inert here.

So `cmake -DCMAKE_BUILD_TYPE=Debug` on MEQ buys the debugger and nothing else,
and the advice to take a debug build "for this reason
alone" **did not work as written** and has been corrected. Recovering those
assertions needs a *second MFEM install* configured `MFEM_DEBUG=YES` — a
different and much larger undertaking, and one nobody has done.

## Rebuilding MFEM: headers, mtimes and cached link lines

**And `make clean` after editing any MFEM header.** MFEM's
makefiles have no `.d` files and no header dependency tracking. That trap has
produced heap corruption in unrelated functions and "unimplemented" aborts for
methods that had just been added. MEQ's own CMake build tracks headers properly;
this applies only to the MFEM tree.

**AND `cmake --install` PRESERVES THE SOURCE FILE'S MTIME, SO A REFRESHED MFEM
INSTALL LEAVES MEQ SILENTLY UP TO DATE ACROSS A CLASS-LAYOUT CHANGE.** The
installed header carries the mtime of the file in `../../mfem/mfem-src`, which is
when that tree was last checked out — routinely EARLIER than MEQ's existing
objects. Measured after one ordinary reinstall: **all 317 installed headers read
older than MEQ's objects**, `cmake --build` of MEQ reported nothing to do, and
`DarcyHybridization` had just gained a member. That links MEQ objects built
against the old layout to a `libmfem.a` built against the new one — an ODR
violation with no error anywhere, which is the failure MEQ's own header tracking
is supposed to make impossible. **Touch the installed headers after every
install** — `find ../../mfem/install/include -name '*.hpp' -exec touch {} +` — and
check the dry run reports work: `make -n` in MEQ's build tree going from 0 to 57 is the
confirmation. This is why "MEQ must then be fully rebuilt too" appears wherever
the install is refreshed; the mechanism is the mtime, not the flags.

**AND A FRESH CONFIGURE RE-SEARCHES FOR WHAT A LONG-LIVED BUILD DIRECTORY WAS
GIVEN BY HAND, WHICH IS HOW A BASELINE STOPS BEING AN ARM.** Configuring a
second MFEM build from the `MFEM_USE_*` flags and TPL directories alone — the
list a reader would write down — resolved **PARDISO to Debian's
`libmkl_sequential`** and **BLAS/LAPACK to Debian's MKL with `libiomp5`**, where
the reference build has oneAPI and `libgomp` throughout. The first fails to
link, which is the good case; the second links two OpenMP runtimes into one
process and runs. `BLAS_LIBRARIES`, `LAPACK_LIBRARIES`, `MKL_PARDISO_LIBRARIES`
and `MKL_PARDISO_INCLUDE_DIRS` are cached **results**, not inputs, and the same
is true one level up: **MEQ's own build directory caches the link line it read
out of `config.mk`**, so repairing the library and rebuilding MEQ in place left
Debian's MKL on MEQ's `link.txt`. Copy the resolved values out of the reference
cache, wipe rather than repoint, and **diff `MFEM_EXT_LIBS` token by token
before believing any comparison between two installs**.
→ **[M-172](MEASUREMENTS.md#m-172)**.

## Traps in MFEM's API

**`ElementTransformation::TransformBack()` CLAMPS AN OUTSIDE POINT ONTO THE
ELEMENT.** Its default solver is `NewtonElementProject`, so a point outside the
element comes back as a reference point *on* its boundary. Any test that reads
"how far outside" off the reference coordinates then reads zero, and accepts a
point outside the mesh. Use an `mfem::InverseElementTransformation` with
`SetSolverType( Newton )` wherever the reference coordinates are examined rather
than the verdict alone. Every site does: `locateFieldPoint`, the X-point's own
search, and `locatePotentialPoint` — which places the limiter contact, `ψ_bnd`'s
point and the prescribed axis, and whose test is `SolverContract`'s
outside-the-mesh refusals for `setAxisFluxPoint()`.

**A `BlockVector` size mismatch that Release builds do not catch.** Stage 2
passed a three-block vector (flux, potential, trace) where `DarcyHybridization`
expects two: `GetOffsets()` stops at the potential, and `ReduceRHS` does
`darcy_rhs = b_t`, whose `BlockVector::operator=` calls
`mfem_error("Number of Blocks don't match")`. Found in stage 4 and
fixed by passing views over `darcy->GetOffsets()`, as
`DarcyOperator::ImplicitSolve` does — **a latent defect in the linear path, not
something Newton introduced.**

**THE EXPLANATION THIS ENTRY USED TO CARRY WAS WRONG.** It said the linear path
survived because its checks are `MFEM_ASSERT`, "compiled out with `NDEBUG`", and
concluded that a debug build was "worth doing now and then for this reason
alone". Neither half holds: `MFEM_ASSERT` is gated on `MFEM_DEBUG`, not on
`NDEBUG`, and the installed MFEM sets `MFEM_DEBUG = NO` — so those checks are
dead in **every** MEQ build and **a MEQ debug build revives none of them**. What
actually caught the semi-linear path was `mfem_error`, which is unconditional.
See *Coverage* above for the audit. The failure mode is still real; only the
advice was.

**A nonlinear potential mass and a linear one do not mix, and how it fails
depends on where the integrators sit.** With *face* integrators on the linear
form (they become `c_bfi_p`) MFEM aborts loudly: "Non-linear mass cannot work
with a linear constraint". With *domain* integrators it is silent —
`LocalNLOperator::AddMultDE` and `ConstructGrad` simply drop them. So leaving the
HDG stabilisation on `M_p` beside a nonlinear source aborts, which is the good
case; the silent one is the reason `buildForms()` documents both.

**`DarcyHybridization` freezes the element-local Newton's initial guess at
`FormLinearSystem()` time, and `ComputeSolution()` will not take another**, so
seeding the vector you hand `RecoverFEMSolution()` is inert. On an ordinary
Newton path that costs iterations; on anything that **differentiates the
residual by differencing it** it costs correctness — measured, it moved
`max ψ_h` from 0.896 to 3.84 for a perturbation of 9e-6. `formSystem()` is the
fix. Full account under `CLAUDE_HDGGS.md`, *The trap that cost the most*.

**`DarcyForm::Reconstruct()` RETURNS A DIFFERENT FUNCTION WHERE `∂F/∂ψ` VANISHES
UNLESS THE FIX IS IN.** It fails silently and per element, so a profile with a
flat segment corrupts part of the domain and leaves the rest exact, and a
whole-domain check misses it. Full account and the transferable lesson under
`CLAUDE_HDGGS.md`, *Post-processing is back*; the part that matters here is that **the fix is on `gf-hdg-dev` and on no
other branch**, so a `meq-integration` rebuilt without it silently loses this.

**`mfem::Mesh::GetElementTransformation( int )` HANDS OUT SHARED SCRATCH, AND
THREADING IT IS A SILENT WRONG ANSWER.** MFEM's own comment: *"The returned
object is owned by the class and is shared, i.e., calling this function resets
pointers obtained from previous calls."* It is one `IsoparametricTransformation`
member of the `Mesh`. Two threads evaluating in different elements overwrite
each other's transformation — no crash, no error, a point transformed by the
wrong element. The reentrant route is the `( i, IsoparametricTransformation * )`
overload into a thread-local.

**SIX CALL SITES IN MEQ TOOK THE SHARED OVERLOAD, WHERE A COUNT MADE BY EYE
FOUND THREE — THE UNDERCOUNT IS THE POINT.** Five in `src/meq/FluxSurfaces.cpp` — `setBandExtension`'s face
loop, `elementSize()`, `locate()`, and **both** branches of `extendField()` —
and one in `src/meq/CriticalPoints.cpp`. All six are **fixed**, each into a
function-local `thread_local` so a transformation held live across a call cannot
be reset underneath it, and every printed number in the four affected
convergence tests is **byte-identical** afterwards.

**AND THERE IS A FOURTH KIND OF SHARED SCRATCH THAT NOTHING HAD NAMED.**
`Mesh::GetBdrFaceTransformations( int )` returns the mesh's own
`FaceElementTransformations`, with the same hazard. The caller-allocated variant
signals failure by `GetGeometryType() == Geometry::INVALID` where the pointer
version returns `nullptr`.

**A `SubMesh` KEEPS A POINTER TO ITS PARENT, AND THREE FIXTURES LET THE PARENT
DIE. LATENT FOR MONTHS, THEN A SEGFAULT IN SOMEBODY ELSE'S CONSTRUCTOR.**
`mfem::SubMesh::CreateFromDomain( background, ... )` stores `&background`.
`ExtensionConvergence`, `FluxSurfaceConvergence` and `FreeBoundaryCoupling` each
built the background as a **local** and returned the `SubMesh` out of the
function, so the parent was dangling the moment they returned. That is undefined
behaviour from the first day and cost nothing for as long as nothing asked a
`SubMesh` where it came from.

`mfem::VertexConePath` now does. Its cone `C(x)` reads the **parent's** edges at
each vertex of `Γ_h` — `HasCone()` is documented as *"whether the mesh handed to
the constructor was a SubMesh with a parent to read edges from"* — so the
constructor walks freed memory and dies in `Mesh::GetVertexToVertexTable`. **It
presents as an MFEM regression and is MEQ's own bug**, which is the reason to
write it down: the backtrace names three MFEM frames and no MEQ frame at all.

The production code was never exposed — `Estimator` holds its background as a
member and `apps/meq.cpp` takes it by reference — so this is a fixture defect,
and the fix is a `static` pool that lives as long as the process. **Anything
that returns a `SubMesh` must return, or otherwise outlive, its parent.**

**AND THE TRIGGER IS NOW GATED, WHICH CHANGES WHAT IS LATENT AND NOT WHAT IS
OWED.** The constructor opens
`const SubMesh *sub = use_cone_ ? dynamic_cast<const SubMesh *>( mesh ) : NULL;`
— so it reaches the parent **only when the cone is asked for**, and `use_cone_`
is `false` by default with MEQ passing it at none of its construction sites. So
MEQ has never used `C(x)`, and this backtrace is not reachable on MEQ's own
paths today. **The cone being off is upstream's measured decision rather than an
oversight, and half the measurement is MEQ's**: it restricts every vertex at
every refinement, `π/2` becoming `π/8`, and changes neither the tiling residual
nor the flux rates — while it costs the far face's quadrature, which is the half
MEQ found. The cost is not a coverage failure but the SMOOTHNESS of
`ξ ↦ a( x( ξ ) )` along a face, and refining the RULE rather than the mesh is
what separates those two: the same cone-on sums fall to the instrument's own
floor at a 40th-order rule. *A caller who wants the cone must raise the face
rule with it.* The lifetime rule above stands on its own — a dangling parent is
undefined behaviour whether or not anything currently dereferences it, and what
the gate bought is time, not safety.

**AND IT IS THE SECOND TIME AN MFEM UPDATE HAS TURNED A LATENT MEQ DEFECT INTO A
CRASH IN MFEM'S OWN CODE.** The first was `SourceIntegrator`'s shared scratch
meeting the threaded element loop. Both were dormant contracts — *the parent
outlives the child*, *an integrator on a threaded loop is reentrant* — that MEQ
had never had to honour. Expect the next library update to find a third.

**AND A MESH ACQUIRES NODES WITHOUT BECOMING CURVED, WHICH COSTS `GridSampler`
TWO ORDERS OF MAGNITUDE AND CHANGES NO ANSWER.** The affine fast path keyed on
`Mesh::GetNodes() == nullptr`; `meq::FieldTransfer`'s constructor calls
`Mesh::EnsureNodes()`, because `FindPointsGSLIB` refuses a mesh without one. So
constructing a `FieldTransfer` dropped every later `GridSampler` on that mesh
onto the Newton route for the life of the process — **159× at 129² and 83× at
513²**, with every located node and every sampled value identical, so only a
clock can see it. The fix is to ask the geometry's **degree**: an order-1 nodal
field on a triangle is the same affine map the vertices are, checked against the
vertex array rather than assumed equal to it, since MFEM lets a caller deform a
mesh through its nodes while leaving the vertices behind.

**It was latent rather than live, and nothing was keeping it that way.** Both
`FieldTransfer`s in `apps/meq.cpp` are built on the stored guess's mesh or on
`previousMesh`, which is a **copy** — so the mesh the driver samples was never
the mesh gslib was set up on. That is an accident of how the adaptive loop
happens to be written, not a rule anybody stated, and it is the second entry in
this section where one component's requirement silently disabled another's
optimisation. → **[M-93](MEASUREMENTS.md#m-93)**

**AND `Mesh::GetGeometricFactors()` IS ONE OF THE THINGS THAT GIVES A MESH
NODES.** On a nodeless mesh it calls `EnsureNodes()`, so the first call
partway through a run changes the mesh under everything set up before it.
Every element and face quantity is bit-identical with and without the nodes;
`Mesh::FindPoints()` is not, by 1.8e-15, and the rotating mirror's
free-boundary solve turned that into 9.06e-10 of `psi_h` (`TODO`). Upstream's
batched routes did it, and so did three of MEQ's device routes. **Take the
factors from a field instead**: `GeometricFactors( const GridFunction &, ...)`
reads any coordinate field and touches nothing, and `meq::CoordinateField`
is that field for MEQ (→ **[M-322](MEASUREMENTS.md#m-322)**). A breakpoint on
`Mesh::EnsureNodes` across a `--device cuda` run is the check; it should hit
only `FieldTransfer`, on the guess's mesh.

**AND THE THIRD IS A DECORATOR, WHICH IS PASS-THROUGH ONLY FOR THE METHODS IT
NAMES.** `mfem::Operator::ArrayMult()` has a base implementation that loops
`Mult()`, so it is never missing and never errors — it is merely slow. A
decorator that overrides `Mult()` and not `ArrayMult()` therefore un-blocks
whatever it wraps: `PardisoSolver`, `CuDSSSolver`, `MUMPSSolver`,
`SuperLUSolver` and `STRUMPACKSolver` all override `ArrayMult` to walk their
factors once for every column, and the base loop sends them `K` separate walks
instead, **with identical answers and identical call counts**. MEQ wraps its
trace solver in `TimedSolver` on every bordered path, so the blocked trace solve
would have been bought and thrown away at the wrapper. The rule generalises past
this one class: **when a library adds a batched entry point with a working
default, every decorator in the chain has to be revisited, and nothing will tell
you.** → **[M-98](MEASUREMENTS.md#m-98)**

**THE ONE THAT ACTUALLY BLOCKS THREADING IS `Mesh::FindPoints`, AND IT CANNOT BE
FIXED LOCALLY.** It loops over every element through that same shared
transformation *and* builds a vertex-to-element table on the way, so it is not
reentrant — an attempt to share one `ContourTracer` across threads aborts with
*"the axis is not in the mesh"*. Every entry point can reach it. The
**unconditional** per-surface call is gone: `traceFromAxis()` seeds the walk with
`CriticalPoint::element`, taking `theTracerClosesAndTheElementWalkDoesNotFallBack`
from 6 calls to 0 over its six surfaces and `SurfaceAverageConvergence` from 196
to 1 over the whole binary. What is left is the data-dependent
last resort inside `fitByAngle()` and `faceJump()`, so a shared tracer still
needs those unreachable or serialised.

**AND "THE TRACER REPORTS ZERO FALLBACKS" COULD NOT HAVE SEEN ANY OF THAT.**
`Contour::fallbackLocations` never counts a **seed** call —  `sampleAt()` hands
`sampleField()` a local counter and discards it on both overloads — so a
per-surface full-mesh scan sat under a test asserting that count is zero,
reading 0 before and 0 after while the real count moved. The numbers above are
from a breakpoint on `mfem::Mesh::FindPoints`, which is the only instrument that
can see them. `INVERSION-PLAN.md` §11.3 has the detail.

**`mfem::Mesh::FindPoints` is `O(elements × points)`** — a brute-force scan over
element centres. It caps sample-cloud sizes in any off-grid error measure.

**MEQ will be an early user of a thinly-tested MFEM combination.** Fixed-boundary
Grad–Shafranov is a Dirichlet problem, so the trace carries an essential BC; and
Newton makes it nonlinear. That pairing —
`DarcyHybridization::SetEssentialBC` together with a nonlinear reduced operator —
was broken on the MFEM branch until recently: `EliminateTraceTrueDofsInRHS`
returned early for nonlinear problems and the essential BC was silently ignored,
so the solver was solving a different problem than asked. **It is fixed**
(`fem/darcy/darcyhybridization.cpp`, and the reasoning is commented there), but
`../../mfem-hdg-dev/CLAUDE.md` records that **no regression covers the combination**.
If stage 4 produces a converged-but-wrong answer near the boundary, look here
first rather than at MEQ's assembly.

## The four defects MEQ reported against the HDG branch

**The four defects MEQ reported to MFEM are closed**, though
`HDG-DEFECTS-FROM-MEQ.md` itself is **not** gone — it is alive on `gf-hdg-dev`
and `gf-hdg-subdomains-dev` and deleted only on the symbolic-reuse line, which
is the modify/delete conflict in the merge recipe — a working-tree listing taken
while that tree sits on the branch that deletes it will say it is gone. Checked
one at a time, because "closed" arrived in four different ways: `Reconstruct()`
on a singular local matrix was **fixed** (MEQ's test flipped red to green on
it); `φ_h` unreachable after a solve was **fixed**, as
`mfem::TransferredDatumCoefficient` in `extension_hdg.hpp` — MEQ still calls
`setTransferredBoundary()` — and rebuilding `η₅` on it, which this sentence used
to say was "MEQ's work and is not done", **is done**: see *A separate `η₅`
problem on the extension path*, where it converges at 2.78 against the pinned
version's 0.40; `ReconstructFluxAndPot()` lifting only domain integrators was
**withdrawn as not a defect**, and MEQ had measured it harmless, which was the
right answer for the wrong reason; and `ComputeHDGFaceEnergy()` ignoring an
installed `HDGStabilization` is **FIXED — `StabValue()` is called in the
function body**, with a comment saying what it returns when no hook is
installed.

**THAT LAST ROW SAID THE OPPOSITE UNTIL RECENTLY AND THE REASON IS WORTH
KEEPING.** It read *"a code read today still shows it computing the `{h⁻¹Q}`
form with no call into an installed hook"*, and that code read used a **fixed
line range** which stopped short of the call — the function is 194 lines and the
call is at 153. Re-read by extracting the function body,
`awk '/ComputeHDGFaceEnergy/,/^}/'`, it is there. A line range is an instrument
and the function is the answer, which is the same species as every other
instrument-not-answer finding in these notes — in the one place where it produced
a false claim about somebody else's code. **All four defects are now closed and
the report is deletable** — and upstream has since deleted it from every working
branch, along with the three other closed ones.
