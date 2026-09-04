# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

meq (the Maryland Equilibrium Solver) computes axisymmetric plasma equilibria by
solving the **Grad–Shafranov** equation with a hybridizable discontinuous
Galerkin (HDG) discretisation built on MFEM.

`refs/Refs.md` indexes the papers behind the numerics, with the doi to fetch each
from. The two marked ✔ there are not background reading — they *are* the method,
and `src/meq` is an implementation of them.

## Status: the solver works, and meq is a program

**meq solves the semi-linear Grad–Shafranov equation by Newton**, which before
stage 2 it had never done at all. `k+1` in both `ψ` and `q` for `k = 1,2,3` over
four dyadic meshes, against an exact Solov'ev equilibrium on the linear path and
HDG-GS-1's Example 5 manufactured solution on the non-linear one, with Newton
converging quadratically. `tests/convergence/SolovievConvergence.cpp` and
`NewtonConvergence.cpp` are the acceptance criteria and print the tables.

* `src/meq/GradShafranov.{hpp,cpp}` is the solver — `DarcyForm` with
  `EnableHybridization`. `Config`, `Profiles`, `Source`, `BoundaryShape`,
  `Estimator`, `Field`, `Sampler`, `WarmStart` and `Output` are all in
  `meq_core` and all under the `naming` check.
* `src/meq/Solution.hpp` is **gone**, not ported: it wrapped hand-rolled block
  offsets and a second set of FE spaces around what `DarcyForm` now owns.
  `v0-legacy` has the original.
* `Estimator.hpp` called `GridFunction::GetValueFacet`, which 4.9.1 does not
  have; it is replaced by `traceFes->GetFaceElement(f)` + `GetFaceVDofs()` +
  `CalcShape()`, exploiting `DG_Interface_FECollection`'s `VALUE` map type — the
  same pattern `estimators_hdg.cpp` uses.
* **`ψ_ax` is an unknown of the non-linear system**, not an input.
  `meq::NormalisedSource` is the interface and
  `setSource( NormalisedSource &, double )` closes it by a bordered Newton, with
  the two non-local terms in the border. `HighBetaConvergence` is the acceptance
  criterion and it is green. `meq::NormalisedMHDSource` is the production source
  built on two `Profile`s and **is not yet reachable from a TOML file**.

**meq is runnable.** `apps/meq.cpp` is the driver, `MEQ_BUILD_APP` defaults
`ON`, and `meq config.toml` parses, builds the mesh and source, solves — with
the adaptive loop if asked — and writes the same equilibrium **three times**,
with exit codes 0/1/2/3 as `DRIVER-PLAN.md` §5 specifies:

| | |
|---|---|
| `.mesh`, `_psi.gf`, `_grad_psi.gf` | exact — every P_k coefficient. GLVis, and meq's own restart |
| `<stem>/<stem>.pvd` + `Cycle000000/` | VTK, **at the degree of the field it draws**, which is `k+1` since `ψ*` became that field — `apps/meq.cpp` passes `polynomialDegree + 1` at four sites. This row used to say "at the polynomial degree of the solve" and was one commit stale. ParaView, VisIt |
| `.nc` | ψ and **B** on a uniform `(R, Z)` grid. Lossy, and the interchange format |

**The VTK is written high-order deliberately**, and it is the one thing in that
row that can be silently wrong: VTK's native cells are linear, so the default
path draws a `k = 3` solution as though it were `k = 1` — a picture that looks
like a coarse mesh rather than like a bug.
`OutputConvergence::theVtkFilesCarryTheHighOrderSolution` asserts the point
count against the vertex count for exactly that reason, and reads **320 points
against 25 vertices** at `k = 3`.

**THE CURVED PATH LEAVES A BAND BETWEEN `Γ_h` AND `Γ`, AND BOTH OUTPUT FORMATS
NOW DEAL WITH IT — DIFFERENTLY, BECAUSE THEY HAVE TO.** `Ω_h` is the union of
background elements lying *inside* `Γ`, so `Γ_h` is inscribed and there is a
band `O(h)` wide that is inside the plasma and outside the mesh.

* **The `.nc` grid continues into it USING THE FLUX**, which is the mixed
  method paying off somewhere nobody expected. `q` is computed at the *same*
  order as ψ and `∇̄ψ = r q`, so a node `p` outside the mesh is reached from its
  foot `x₀` on `Γ_h` as `ψ(x₀) + r₀ q(x₀)·(p − x₀)` — **nothing is ever
  evaluated outside an element**.

  **The obvious alternative was implemented first and is bounded by nothing.**
  Continuing `ψ_h`'s own polynomial past its element put **17 nodes at positive
  ψ** on `miller-curved`, where `[boundary] Type = "zero"` makes ψ exactly zero
  on `Γ` and strictly negative inside — worst +1.06e-02 against a peak of
  2.5e-01, and the ψ = 0 contours in that band were visibly wrong. The flux
  version puts **none** there, with a maximum of −5.9e-05, and its band error
  converges at 3.75 where an extrapolation does not converge in the band at all.
  `theBandIsContinuedByTheFluxAtTheFluxesOwnOrder` asserts the rate.

  **Blending the extrapolation toward the known ψ = 0 on `Γ` does NOT fix it**,
  which is worth recording because it looks like it should: `(1−t)·v` scales a
  positive value down and never changes its sign, so all 17 nodes survived it.
  The error was in *where the field was evaluated*, not in how it was weighted.

  **`B` GETS THE BAND TOO, SINCE 2026-09-02, AND IT USED TO BE ZEROTH ORDER
  THERE BEHIND AN `inside = 1` MASK.** Only `samplePotentialWithFlux()` applied
  the Taylor step; `sample()`, `sampleComponent()` and `sampleCoefficient()`
  never read `offsetR`/`offsetZ` at all, so for a band node they return the
  value **at the foot on `Γ_h`** — and `apps/meq.cpp` sent `ψ` through the first
  and `B` through `sampleComponent`. Measured on `miller-curved.nc`:
  `extrapolated_nodes = 1667` of 129×129 = 16,641, so about **one node in ten
  carried a piecewise-constant `B`**, `O(h)`, in a field `q` itself resolves at
  `k+1` and that `ψ` receives at `O(h²)` in the same band. `B = (−q_z, +q_r)` is
  a pure relabelling, so nothing softened it.

  **`sampleComponentWithGradient()` is the fix**, and it takes the same Taylor
  step from the same foot, using the field's own `∇u` — which is read *inside*
  the element, so the guarantee that nothing is evaluated outside one survives.
  `theBandVectorContinuesAtItsGradientsOrder` measures both sides on a quadratic
  its own space represents exactly, so the band error is the truncation alone:
  **rate 2.20 with the gradient against 0.92 reading the foot, and 124× smaller
  at `n = 16`**, interior nodes exact to 2.2e-15.

  **IT DOES NOT REACH `ψ`'s ORDER AND THE GAP IS STRUCTURAL, NOT A SHORTCUT.**
  `ψ` is continued with `q`, a **solved** variable carrying the potential's own
  order — that is the mixed method paying off. There is no solved variable for
  `∇q`: differentiating an L2 field of degree `k` leaves `k−1`, so this is
  `O(h²)` at every `k`. A full order better than what it replaced, and not the
  same thing. **The route to better is known and was not taken**: `div q = −F/r`
  is the equation being solved and `∂_r q_z − ∂_z q_r = −q_z/r` follows from
  `r q = ∇̄ψ`, which pin two of `∇q`'s four entries exactly — but they leave the
  symmetric traceless part still differentiated, so they buy *structure* rather
  than an order, at the cost of plumbing the source into `GridSampler`.

  **AND THE COST IS MEASURED, NOT ESTIMATED.** Wiring the rotating output gave
  a free controlled experiment: `n_s` is algebraic in `(r, ψ)`, so the same
  closed form can be evaluated both ways over the same 356 band nodes. Against
  the node's own `r` it is **5.13e-07** wrong; against the foot's `r`, which is
  what `sampleCoefficient` hands you, it is **8.57e-02** — a factor of
  **1.7e5**. `n_s` gets the radius wrong twice over, since the exponent of (96)
  carries `r²`; `B` gets it wrong once, being a relabelling of `q`. The
  experiment is `theRotatingFieldsUseTheNodesOwnRadius` in
  `OutputConvergence.cpp`, and it exists because the trap was met and dodged for
  the new fields while `B` was left in it two lines above.

  **AND THE FILE NOW SAYS WHICH NODES THOSE ARE**, which was the other half of
  the defect and is not fixed by making `B` better. `extendOutward()` counts band
  nodes as found, so `located()` is true and `inside` says 1 there — correctly,
  since the node is in the plasma and carries real data — and
  `theMaskAgreesWithTheData` asserts mask-matches-non-NaN, which is precisely the
  property that made the band look trustworthy. `extrapolated_nodes` was a
  **count, not a mask**, so nothing downstream could tell *which*. The `.nc`
  carries a second `byte extrapolated( Z, R )` beside `inside`, from
  `GridSampler::wasExtended()`; on `miller-curved` it sums to 1667, agreeing with
  the attribute, and is a strict subset of `inside`. Drop those nodes before
  computing an error norm or differencing two runs. **This is the interchange
  format, and `B` is what most readers open it for**, which is why it was an item
  rather than a footnote.
* **The `.vtu` bends the mesh onto `Γ`** — a curvature is installed and each
  boundary face is moved out. Since the VTK is already Lagrange cells this
  needed **nothing further from the format**; the two features composed.
  `theBoundaryBendsOntoTheTrueGamma` asserts the nodes land on the target to
  **1.1e-16**.

**Two things about that bending were got wrong on the way and are worth not
repeating.** Smoothing the displacement into the interior on the *vdof* graph
couples R to Z — they share one index range — so a radial displacement gets
averaged against a vertical one; measured, that made tangling **worse**, 25%
surviving where moving the boundary alone managed 50%. And backing the
displacement off **globally** costs every face the worst face's limit: per-node
backoff takes `miller-curved` from 50% to **96%** of the boundary reaching `Γ`.
The gap can exceed an element's own size, so some faces genuinely cannot reach,
and the driver reports the fraction that did.

**AN ADAPTIVE RUN ALSO WRITES `<stem>_cycles/`, ONE VTK FRAME PER CYCLE**, which
ParaView scrubs through as a time series — 97 → 254 → 342 → 449 elements over
`miller-adaptive`'s four cycles. It is a separate collection from the answer
because `<stem>` gets its boundary bent onto `Γ`, and doing that mid-loop would
hand the next refinement a geometry the estimator never saw.

**One frame per cycle is easy to get ALMOST right, and the near miss is
invisible.** Rebuilding the collection per frame puts every `Cycle` directory on
disk with the correct refined mesh, and leaves the `.pvd` index listing only the
last of them — so ParaView opens the file and shows a single frame, with no
error and all the data present. `ParaViewDataCollection` appends to its `.pvd`
and does not scan the directory, so the collection has to survive between frames
and be rebound with `SetMesh()`. `theAdaptiveSeriesIndexesEveryFrame` asserts on
the **index**, not the pieces, because the pieces were never what broke.

`tools/README.md` is the guide to which format goes with which reader.
`DriverAcceptance.cpp` asserts the driver reproduces the *library* on the same
configuration — 1.189e-16 over 15,360 dofs — rather than comparing against a
closed form, for the reason recorded beside `examples/soloviev-nstx.toml`.

**The curved boundary works through the driver.** `[boundary.shape]` builds the
shape, marks `D_h`, locates `Γ_h` and hands a `VertexConePath` to
`setExtension()` — GS-2's technique end to end from a config file;
`examples/miller-curved.toml` is the worked example. **Two checks guard it,
because the failure mode is quiet**: `Γ_h` carrying no transferred datum
silently imposes zero and still converges, so `theDriverSolvesOnACurvedBoundary`
pins the driver against the library (1.6e-16) *and* against that zero-datum
solve (2.7e-1, i.e. the transfer is doing real work).

**The adaptive loop works through the driver, on both paths.** solve →
post-process → estimate → mark → refine, stopping at `TargetError` or
`MaxIterations` and saying which. On the curved path it uses
`meq::AdaptiveDomain` — the companion mesh of GS-2 §3.3, without which
`dist(Γ_h, Γ)/h_loc` doubles every cycle and the transfer silently leaves the
regime it is analysed in — and calls `setTransferredBoundary()` automatically.
`examples/miller-adaptive.toml`: η monotone 4.73e-4 → 6.87e-5 over four cycles,
97 → 449 elements, 0 transfer paths widened. `theDriverRunsTheAdaptiveLoop` pins
it against the same loop driven through the library at **1.666e-16 over 2694
dofs**, and separately asserts that it refined, that it refined *adaptively*,
that η came down, and that assumption P.1 survived the graded `Γ_h`.

**A `TargetError` chosen against a pre-2026-08-29 run is not comparable**: η is
built on `ψ*` now rather than `Potential::Raw`, which is a different and smaller
quantity on the same mesh. See *Post-processing is back*.

**The driver's non-linear path is a reactive ladder**: Newton, and on *observed*
failure `PicardThenNewton`, rebuilding the solver because a caught
`ErrorException` leaves one unusable. **Never predictive** — nothing may be
inferred from `F` about which solver to run, and two candidate detectors have
now been measured and killed. See *Why meq's Newton struggles* and *There is no
cheap discriminator*.

**What the driver refuses rather than approximates**: `[boundary] Type =
"exact"` needs a closed form `meq::Source` does not carry, and exits 1 with an
explanation. **That is now the only thing it refuses.**

**THE INTERPOLATING WARM START IS WIRED, 2026-09-02, AND THIS PARAGRAPH USED TO
SAY IT "NEEDS GSLIB AND IS NOT WRITTEN".** Both halves of that were false and had
been for months: `MFEM_USE_GSLIB` is `YES`, `meq::FieldTransfer` is
`FindPointsGSLIB` in `WarmStart.{hpp,cpp}`, and `WarmStartConvergence` was
already a registered ctest. Only the driver-side wiring was missing, in two
places, and both are done:

* **A stored guess on a DIFFERENT mesh** used to throw. It now transfers, so
  restarting from a run at another resolution — the ordinary way to use a stored
  answer — works. The exact restart is still taken when the meshes match, since
  it is every coefficient rather than an interpolation.
* **Adaptive cycles after the first used to start COLD**, throwing away a
  converged answer at every refinement. They now start from the previous cycle
  interpolated onto the refined mesh, with the Dirichlet datum as the fallback
  at nodes the coarse mesh does not cover — and on the curved path there really
  are such nodes, since the computational domain grows as it refines: measured,
  48 and 360 of them on `miller-adaptive`.

**WHAT IT BUYS, AND THE MEASUREMENT NEEDED A NONLINEAR SOURCE TO EXIST AT ALL.**
meq's adaptive examples run a Solov'ev source, whose `∂F/∂ψ` is zero, so Newton
takes one step whatever it starts from — `miller-adaptive` shows `1` in every
cycle warm or cold, and the warm start is unmeasurable there by construction. On
`ManufacturedNonlinear` over three cycles at `k = 2`,
`carryingTheAnswerAcrossCyclesCutsTheWork` reads:

| | per cycle | total | final L2 |
|---|---|---|---|
| cold | 4, 4, 4 | 12 | 3.652757e-06 |
| **warm** | **4, 2, 2** | **8** | 3.652757e-06 |

**A third of the Newton work, and the same answer to 2.5e-12.** Both halves are
asserted: strictly fewer iterations, and an L2 that does not move — a warm start
must change the work and not the answer.

**`theDriverRunsTheAdaptiveLoop` now reads 4.4e-14 rather than 1.7e-16** for
exactly this reason, and that is a stronger assertion rather than a weaker one:
the driver warm-starts each cycle and the library loop it is pinned against
deliberately does not, so the number is now a statement that the warm start left
the equilibrium alone.

**The non-linear ordering is `NonlinearOrdering::NPC`**, since MFEM deleted the
mode meq's default used to be. `CondenseThenLinearise` is kept as the backup and
is still the one that converges on stiff under-resolved meshes, at three to four
times the wall clock everywhere else. See *The NPC port*.

**~~`README.md` overstates what the old code did and has not been rewritten~~ —
REWRITTEN 2026-09-02, and it now says this itself.** The fact is worth keeping
because it is the reason to distrust anything written about meq before the port:
meq then solved the *vacuum* coil field and nothing else — `meq.cpp`'s
right-hand side took `psi` and ignored it, `Configuration::plasma` was hardcoded
`nullptr`, and the profile loader was a `{ return; }` stub that silently
produced an empty spline. The four test files that existed were empty Boost
stubs, and one — asserting `foo == bar` with neither declared — could not
compile. **`docs/` is now the user-facing account** and `README.md` points at
it; this file stays the maintainer's one. See *Layout* for why the numbers live
here and not there.

**`ROADMAP.md` is the priority order**, and says which items belong to meq and
which are requests on `../mfem-hdg-dev`, where another agent is working. This
file stays the technical record; that one is only about what to do first.

### The stages, and where we are

| | | |
|---|---|---|
| 0 | Git reconciliation, tag `v0-legacy` | **done** |
| 1 | Tree, CMake, `Config` / `Profiles` / `Source` | **done** |
| 2 | Linear `Δ*` on `DarcyForm`, fitted polygonal domain | **done** |
| 3 | Local post-processing `ψ*_h` | **dropped as meq code** — `DarcyForm::Reconstruct()` supplies it, at `k+2`. See below |
| 4 | Newton on the semi-linear source | **done** |
| 5 | Curved `Γ` by extension from subdomains | **done** |
| 6 | Adaptivity: the residual estimator and mesh update | **done** |

Each stage ends at a **measured convergence rate**, not at "it runs". See
*Testing stance* below for why that is the acceptance criterion.

## Commands

```sh
git submodule update --init --recursive     # extern/toml11
cmake -B build
cmake --build build -j4
cd build && ctest --output-on-failure       # ~490-680 s, 33/33
```

**ctest needs no environment set by hand.** `tests/CMakeLists.txt` puts
`MKL_NUM_THREADS=1` on every registered test, without which the suite takes well
over an hour instead of nine. It is the only variable meq needs; for the one
this file used to insist on and no longer does, see *Traps*. Running a test
binary **directly** needs the same one thing:

```sh
MKL_NUM_THREADS=1 ./tests/SolovievConvergence
```

The performance harness is separate and is not a ctest:

```sh
tests/performance/scan.sh build 3            # both thread axes, ~40 min
```

`MFEM_DIR` defaults to `../mfem/install` and also reads the environment, so the
`-D` is usually unnecessary. Never a bare `make -j` or `cmake --build -j`; see
*Traps*.

### Coverage, and what CI can and cannot check

```sh
cmake -B build-cov -DMEQ_ENABLE_COVERAGE=ON && cmake --build build-cov -j4
cd build-cov && ctest
gcovr --root .. --filter 'src/meq/' --print-summary        # or --html-details
```

`MEQ_ENABLE_COVERAGE` is an option, not a build type, on purpose: a coverage
build wants `-O0` and a release build wants `-O3`, and that difference should be
chosen rather than inherited from whatever `CMAKE_BUILD_TYPE` happens to hold.

**THE REASON THIS FILE USED TO GIVE WAS WRONG, AND IT IS THE KIND OF WRONG THAT
COSTS SOMEBODY AN AFTERNOON.** It said routing `-O0` through
`CMAKE_BUILD_TYPE` would drop `NDEBUG` and change which `MFEM_ASSERT`s are
live. It would not. Audited 2026-09-02:

* `MFEM_ASSERT` is gated on **`MFEM_DEBUG`**, never on `NDEBUG`
  (`general/error.hpp`), and nothing anywhere in the installed tree derives one
  from the other.
* The installed MFEM has **`MFEM_DEBUG = NO`**, so **every `MFEM_ASSERT` is
  already dead in every meq build**, release or debug.
* meq's own sources contain **no `assert()` and no `<cassert>`**, so `-DNDEBUG`
  disables nothing of meq's either. It is inert here.

So `cmake -DCMAKE_BUILD_TYPE=Debug` on meq buys the debugger and nothing else,
and the standing advice under *Traps* to take a debug build "for this reason
alone" **did not work as written** and has been corrected. Recovering those
assertions needs a *second MFEM install* configured `MFEM_DEBUG=YES` — a
different and much larger undertaking, and one nobody has done.

**CI CANNOT BUILD THE SOLVER, AND THIS IS STRUCTURAL.** The MFEM meq needs is
`meq-integration`, a local merge published on no remote, so a hosted runner
cannot obtain it and caching does not help — there is nothing to fetch.
`find_package(MFEM)` is therefore **not `REQUIRED`**: without it `CMakeLists.txt`
builds the MFEM-free half — `Config`, `Profiles`, `Source`, `SourceFactory` — and
`tests/CMakeLists.txt`'s existing "unit not present" guard skips every
convergence test without needing to know why. So `ci.yml` runs four unit suites,
the `naming` check and a coverage gate at 90% lines. **Treat a green CI badge as
evidence about the configuration and profile layers and about nothing else** —
the `full` job is written but `if: false` until the branch is published or a
self-hosted runner exists.

**EVERY COVERAGE FIGURE THIS FILE USED TO CARRY IS DELETED RATHER THAN
UPDATED.** They predated the bordered Newton and the NPC port —
`GradShafranov.cpp` has gained about 950 lines and two whole solve paths since —
and quoting a stale percentage is worse than quoting none. Two things about the
*measurement* survive, and they are the transferable part.

**A line-coverage percentage on this codebase is partly a measurement of comment
density.** On the same run, `GradShafranov.cpp` read **92.7%** to `gcov`
(executable lines only, 586 of them) and **64%** to `gcovr` (its own line set,
838). This project comments heavily, and `gcovr` counts roughly 250 lines that
`gcov` does not consider executable at all — then reports most of them as
uncovered. The CI gate uses `gcovr`, so that is the operative number for the
gate; `gcov`'s is the operative number for "is this code tested".

**A claim this file carried, and which was wrong, is worth leaving recorded**:
that the solver sat at 60% against 92–100% for the configuration layer, so "the
best-covered code in meq parses TOML and the least-covered is the part every
physical claim rests on". The 60% was `gcovr`'s number compared against
`gcov`-flavoured intuitions. On one denominator the solver and `Config.cpp` were
the same. **A measurement about the tool, not about the code** — the same trap
as the threaded-MKL one.

What is true regardless of the denominator is that **22 of
`GradShafranovSolver`'s 51 methods were called by nothing**, and
`SolverContract.cpp` was written for that rather than for a percentage — it
found `setSource()` silently accepting a second source of the same kind, and
value faults throwing `logic_error` where the constructor throws
`invalid_argument`. **Branch coverage, last seen at 48%, is the figure with real
headroom in it and the one nobody has looked at.**

**`gcovr` is not installed on this machine**, so the recipe above fails with
`command not found`; a venv is the way round it, and it needs
`--merge-mode-functions=separate` or it aborts on `Config.hpp`'s inline
accessors appearing at two line numbers across translation units. Run it from
the repo root: `--root .` after a `cd` into the build directory silently matches
nothing and reports 0%.

### Why toml11 is a submodule and MFEM is not

**toml11 is `extern/toml11`**, pinned to a commit, as MaNTA does it. It is header
only, small, and its version is something meq should control rather than inherit
from whatever is installed — the `find_or<double>` behaviour recorded under
*Traps* is version-dependent, and a silent config bug is exactly what a floating
dependency buys you. Currently `b32a2ff`, v4.4.0-31, the same commit MaNTA pins.

**MFEM stays out of tree**, hand-built, found through `MFEM_DIR`. Its history is
enormous, so a submodule would make every clone of meq expensive; it needs an
out-of-tree configure-and-build of its own anyway; and `../mfem-hdg-dev` is a
working tree with its own active development on `gf-hdg-subdomains-dev`, which a
submodule pin would fight rather than help.

If configure fails with toml11 not found, the cause is almost always a clone
without `--recursive`; the error message says so.

### Prefer a maintained library to a hand-rolled algorithm

**Standing preference, stated 2026-09-03, and it is about maintenance rather
than speed.** Where a well-known algorithm has a well-maintained
implementation — a Householder QR, a rank-revealing least squares, a special
function, a spline — **take the library**, and take it even when a measurement
says the hand-rolled version is no slower. *"The chance that you found the
perfect Householder implementation and are able to maintain it indefinitely
into the future is low."*

A performance tie is therefore **not** an argument for keeping bespoke code. The
bar for bringing in a header-only, well-known, well-maintained library is that
it does not break something else — not that it wins a benchmark.

What this does **not** license: a dependency that is heavy, obscure, unmaintained
or hard to obtain, and it does not override the reasons recorded above for why
MFEM is out of tree and toml11 is pinned. And it is not a reason to rewrite
working code on sight — it decides which way to go when the question is live.

Already taken on these grounds: **Boost.Math** for the Zernike radial
polynomials, where the Jacobi route is both more accurate than the textbook
factorial sum and somebody else's to maintain — see *The disc basis*. **Eigen**
(3.4.0, `/usr/include/eigen3`, header only) is the natural home for the dense
linear algebra in `src/meq/SurfaceFit.cpp`, which currently hand-rolls a
Householder QR and a one-sided Jacobi SVD in about 150 lines.

**The one thing to hold fixed when swapping in a library**: `SurfaceFit`'s
truncation threshold and its Levenberg–Marquardt damping together *are* the
gauge that makes the gauge-free fit well posed, and a different rank-revealing
decomposition can truncate different directions. The ellipse result, the `nstx`
tail, the exactly-zero axis spread and the positive minimum Jacobian are the
numbers that say the swap was clean. If one moves, that is a finding, not a new
baseline.

## The equation being solved

Fixed boundary: the plasma boundary `Γ` is known and taken to be the level set
`ψ = 0`, so this is an interior Dirichlet problem. From
`refs/HDG-GradShafranov.pdf` eqs (1)–(4):

```
-∇̄·( (1/r) ∇̄ψ ) = F(r,z,ψ) / r     in Ω ⊂ R²
ψ = 0                               on Γ = ∂Ω

F(r,z,ψ) := μ₀ r² dp/dψ + g dg/dψ
```

`∇̄ := (∂_r, ∂_z)` acts formally like a vector of partial derivatives independent
of the coordinate system — it is *not* the cylindrical gradient, and the
distinction is the whole content of the `1/r` and `r` weights below. `p(ψ)` is
the plasma pressure and `g(ψ)/r` the toroidal field function; both are user
input, and it is their `ψ`-dependence that makes the problem semi-linear.

Recast as a first-order system by introducing the **flux** `q = (1/r)∇̄ψ`:

```
q − (1/r)∇̄ψ = 0,      −∇̄·q = F/r,      ψ = 0 on Γ
```

Introducing `q` is not a numerical convenience. The physically interesting output
is the magnetic field, which is built from `∇ψ`, so discretising `q` directly is
the reason to prefer HDG here — it gives the derivative at the same order as the
potential rather than one order down.

### The two papers disagree about the sign of the Solov'ev source

Checked, resolved, and worth not rediscovering.

`refs/HDG-GradShafranov.pdf` eq. (10) gives `F = −((1−A)r² + A)`.
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (21) gives `F = +((1−A)r² + A)`.

**The first is right**, and the second contradicts its own eq. (1). Applying
`Δ*` to the particular solution *both* papers publish settles it analytically:

```
Δ*( r⁴/8 )            =  r²
Δ*( (A/2) r² ln r )   =  A
Δ*( −A r⁴/8 )         = −A r²
                        ─────────────────
Δ*( ψ_P )             =  (1−A) r² + A
```

and since both papers define `−Δ*ψ = F`, `F = −((1−A)r² + A)`. The twelve
homogeneous terms contribute nothing, being `Δ*`-harmonic. This is confirmed
numerically: `Soloviev.hpp`'s `deltaStarFD()` recomputes `Δ*ψ` by central
differences, and `SolovievConvergence.cpp`'s `solovievSourceMatchesTheOperator`
asserts it against `−f()` over the benchmark rectangle, so the whole twelve-term
transcription is checked rather than trusted, and checked against the same `f()`
the solver is actually fed.

Take this as the standing warning about the benchmark: **the published
coefficients are not self-checking**, and a sign error here would show up as a
solver that converges beautifully to the wrong equilibrium.

### A third erratum: eq. (20)'s `η₅` does not vanish on the exact solution

`ψ̂_h ∈ P_k(e)`, but `ψ*`'s trace has degree `k+1`, which no element of `M_h` can
represent. Orthogonality splits the term:

```
‖ψ̂ − ψ*‖²_e = ‖ψ̂ − P_M ψ*‖²_e + ‖(I − P_M) ψ*‖²_e
```

The second piece survives even when you substitute the *exact* solution and the
best possible trace, at `O(h^{k+3/2})`, which `h_e^{-1}` and `O(h^{-2})` edges
turn into an `O(h^k)` floor. Measured, the printed term **is** that floor —
agreeing with it to between 0.05% and 1.8% at every `k` and every mesh — so it
converges at `k`, not `k+1`, and drags the total to 1.44/2.32/3.55 instead of
1.99/2.99/3.97.

Taking the difference *inside* `M_h` restores `k+1` and reproduces the rates
GS-2's own Table 1 reports. `TraceComparison::Projected` is the default;
`Literal` is kept so the suite keeps measuring the difference.

**A separate `η₅` problem on the extension path**, and this one was nearly fatal
to the adaptive loop: on `Γ_h` the term compares `ψ*` against a trace pinned to
*zero* rather than the `φ_h` actually imposed, so the difference is
`O(dist(Γ_h, Γ)) = O(h)`. Unmitigated, `η = 4.09e-1` where `η₁ = 2.12e-3`,
converging at ~0.5 — the loop would have run, produced plausible pictures, and
refined the wrong elements.

**THAT WAS AN OMISSION AND IT IS NOW A REPAIR, 2026-09-02.**
`setTransferredBoundary()` used to *exclude* those faces, which restored `k+1`
by deleting the term exactly where the geometry error lives. It now takes the
datum as a second argument — `GradShafranovSolver::transferredDatum()`, built on
`mfem::TransferredDatumCoefficient` — so `η₅` compares `ψ*` against the `φ_h`
actually imposed and the faces stay in. Measured on the extension benchmark at
`k = 2` over `h = 0.213 / 0.106 / 0.053`, all three treatments on one solve:

| | coarsest | finest | rate |
|---|---|---|---|
| `η₁`, for scale | 2.1161e-03 | 4.7256e-05 | — |
| **`η₅` on the datum** | **9.5889e-05** | **2.0281e-06** | **2.78** |
| `η₅` on the pinned zero | 4.0688e-01 | 2.3226e-01 | 0.40 |
| `η` with the datum | 2.1746e-03 | 4.8375e-05 | 2.75 |
| `η` with the faces excluded | 2.1745e-03 | 4.8375e-05 | 2.75 |

So the term is now **two percent of `η₁`** rather than four orders larger than
it, and it converges at `k+1` rather than at a half.
`theTransferredDatumRestoresEtaFive` in `ExtensionConvergence.cpp` is that table,
and it keeps the pinned column as its control: if that ever converges, the
comparison is empty.

**η BARELY MOVES, AND THAT IS THE POINT RATHER THAN A DISAPPOINTMENT.** The
excluded and datum columns agree to four figures, because a correctly evaluated
`η₅` on `Γ_h` is small. What changes is that the boundary elements now *have* an
`η₅` contribution instead of none, so the marking sees them. The adaptive loop is
unchanged — 97 → 254 → 342 → 449, `η` 4.7352e-04 → 6.8668e-05 — which is what
says the repair did not perturb what already worked.

**THE SIGN IS THE THING TO GET RIGHT, AND IT IS TESTED SEPARATELY.**
`transferredDatum()` must hand `mfem::PathLiftCoefficient` the **raw flux
block**, which holds `−q`, and not `flux()`, whose sign is undone. Feed it the
wrong one and `φ_h` comes back as `−ψ` rather than `ψ`, because the `g(a(x))`
that would otherwise survive is zero on `Γ`. So the failure is the answer with
its sign reversed, not a small bias.
`theTransferredDatumReproducesTheImposedCondition` checks `φ_h` against the
exact `ψ` it transfers — 1.65e-3 → 4.0e-5 relative, converging at 3.25 then 3.63
— and a sign error reads about 2.

## The discretisation

HDG, in the LDG-H form of `refs/HDG-GradShafranov.pdf` eq. (8) — restated in
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (13) with the block structure made
explicit, which is the more useful version when writing assembly.

```
(r q_h, v)_Th + (ψ_h, ∇̄·v)_Th − ⟨ψ̂_h, v·n⟩_∂Th   = 0
(q_h, ∇̄w)_Th − ⟨q̂_h·n, w⟩_∂Th                    = (F/r, w)_Th
⟨q̂_h·n, μ⟩_∂Th\Γh                                 = 0
ψ̂_h = φ_h on Γh

q̂_h·n := q_h·n + τ(ψ_h − ψ̂_h)        on ∂Th
```

Spaces, all of the **same degree** `k` — hybridization removes the inf-sup
compatibility condition that a classical mixed method would impose:

```
V_h = [P_k(K)]²    flux            L2_FECollection, vdim 2
W_h = P_k(K)       potential       L2_FECollection
M_h = P_k(e)       hybrid trace    DG_Interface_FECollection
```

**The volume spaces are on the closed Gauss-Lobatto basis, and that is a
convention rather than a requirement** — `convdiff.cpp`'s "as it is customary for
HDG to match trace DOFs", quoted and inherited until it was measured. A nodal
basis does not change the *space*: Lobatto and Legendre both span `P_k(K)`, so
the discretisation cannot see the choice, and nothing in the hybridization needs
volume dofs to sit on faces — every coupling is a face **integral**, computed by
quadrature against both bases. Measured, with both volume spaces switched to
`GaussLegendre` and nothing else touched: rates 1.996/3.001/3.998 in `ψ`, the
curved benchmark's `L2` at 2.742813e-05 down both paths agreeing to 6.6e-15, and
the `k = 3` Newton history reproducing the one recorded below to seven digits.
The finest Solov'ev error reads 5.510484e-10 against 5.510486e-10 — the last
place.

So it is kept for alignment with the miniapp meq was ported from, and for nothing
else. **What it costs** is that a dof is a point value *on* the element boundary,
where an L2 field is discontinuous, so reading `W_h` by nodal interpolation at
another mesh's dof points is ambiguous — measured, 9% to 28% wrong.
`meq::FieldTransfer` projects rather than interpolates, which is basis-agnostic
and is what a non-nested transfer needs anyway.

### The assembled flux is −q, and the papers' τ carries the wrong sign

Two conventions settled by measurement rather than argument, both of which cost
time if rediscovered. `tests/convergence/SolovievConvergence.cpp` records what
every alternative produces.

**`DarcyForm` holds `−q`, not `q`.** It is built for `u = −k ∇p`, the opposite
sign to `q = (1/r)∇̄ψ`, and the two integrators that make the hybridization
consistent — `NormalTraceJumpIntegrator` and the trace rows of
`HDGDiffusionIntegrator` — have that sign baked in and take no scaling argument,
so there is no way to flip it in the assembly. `GradShafranovSolver::flux()`
undoes it once, into a separate GridFunction; the block vector stays in
`DarcyForm`'s convention, because a stage-4 Newton residual assembled by
`DarcyForm` will expect it there. Do not change one without the other.

The same convention makes the potential right-hand side `−(F/r, w)`: with the
default `bsymmetrize = true` the second block row is assembled as
`−B q − M_p ψ = b_p`. `bsymmetrize = false` is not an escape — it does not work
here at all, giving an error flat at 3e-1 for every combination of signs.

**And `τ` carries the opposite sign to `refs/HDG-GradShafranov.pdf` eq (8e)**,
which prints `q̂·n := q·n + τ(ψ − ψ̂)`. With that paper's own `q = (1/r)∇̄ψ` and
`−∇̄·q = F/r`, testing (8a)–(8d) against `(v,w,μ) = (q_h, ψ_h, ψ̂_h)` gives

```
(r q, q) − τ‖ψ − ψ̂‖²_∂Th = 0
```

which is indefinite, and the local solves are not guaranteed invertible. The
stable sign is `−τ`, which is exactly what assembling in `DarcyForm`'s convention
with a positive `τ` produces. This is a second sign slip in the same pair of
papers as the Solov'ev source one above. It is a well-posedness argument, not an
observed failure: **both** signs converge at `k+1` on this benchmark, with
`τ > 0` in `DarcyForm`'s convention giving slightly the smaller error at `k = 1`.

**`τ = 1`.** Both papers set it there and note that optimal order needs only
`τ = O(1)`. The pre-port code used `τ = 5.0` with no recorded reason; do not
reinstate that without a measurement to justify it.

**And MFEM will not give you a constant `τ` unless you ask.**
`HDGDiffusionIntegrator`'s built-in stabilisation is

```
τ± = ( β ± (α/2)(u·n)/|u·n| ) { h⁻¹ Q }
```

— the LDG choice, scaled by the inverse local mesh size *and* by the diffusion
coefficient, and not what the papers use. Measured on the Solov'ev benchmark it
costs a full order in the flux: `q` converges at `k`, not `k+1`, while `ψ` still
converges at `k+1`, so **a study of `ψ` alone would have passed it**.

The designed way out is the `HDGStabilization` hook in
`fem/darcy/bilininteg_hdg.hpp`: subclass it, return the constant from `Eval()`,
install it with `SetStabilization()`. `meq::ConstantStabilization` is that
subclass. Read `StabValue()` there — with a hook installed it divides the
quadrature weight out before calling and multiplies it back after, so returning
a bare `τ` yields exactly `⟨τψ, w⟩`. The MFEM tree has no ready-made constant
implementation; the only subclasses are the test fixtures in
`tests/unit/fem/test_bilininteg_hdg.cpp` and `test_darcy_degenerate.cpp`, which
are the idiom to copy.

**Keeping `IsConstant()` true also keeps meq out of an `EvalGrad` trap.** That
header warns that omitting `EvalGrad` for a *non-constant* stabilisation gives
"no wrong answer, only slow Newton convergence — a failure that survives a
passing regression suite". A constant `τ` never calls it. One more reason to
keep `τ` constant unless something measured says otherwise, which is also the
sibling project's standing advice: **do not derive `τ` from the local
coefficient.**

### Which MFEM, and why not master

**meq builds against `../mfem/install`** — MFEM **4.9.1** on branch
**`meq-integration`**, CMake-built in `../mfem/build` from sources in
`../mfem/mfem-src`. The options actually set, read out of
`share/mfem/config.mk` and re-checked 2026-09-01 rather than remembered:

| | |
|---|---|
| `MFEM_USE_SUNDIALS` | KINSol, so `KINSolver(KIN_LINESEARCH)` is reachable. **Now `../sundials/cuda-install`, not `../sundials/install`** — see the CUDA row |
| `MFEM_USE_GSLIB` | `FindPointsGSLIB`; gslib v1.0.9 built alongside at `../mfem/gslib` |
| `MFEM_USE_SUITESPARSE` | UMFPACK, the direct solver meq's own solver runs on. **`../suitesparse/install` since 2026-09-01, not Debian's** — v7.12.2, the same version Debian ships, built against oneAPI MKL. Debian's carried a `NEEDED` on `libblas.so.3` -> `libmkl_rt.so`; see *PARDISO and the MKL link line* |
| `MFEM_USE_MKL_PARDISO` | `PardisoSolver`, oneAPI MKL 2026.1, **on the threaded `mkl_gnu_thread` layer since 2026-08-30**. It used to be resolved to Debian's `mkl_sequential`, which is a bigger deal than it sounds — see *Traps* |
| `MFEM_USE_EXCEPTIONS` | `MFEM_ERROR_THROW` by default, which is what makes the driver's exit code 2 reachable |
| `MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE` | **both, and both are now load bearing.** `DarcyHybridization::SetAssemblyMode( Threaded )` needs both and **aborts** rather than falling back without them. meq reaches it through `setAssemblyMode()`. This row used to read "not exploited"; that is no longer true |
| `MFEM_USE_LAPACK`, `MFEM_USE_ZLIB`, `MFEM_USE_DOUBLE`, `MFEM_USE_MEMALLOC` | |
| `MFEM_USE_MPI = NO` | meq is serial throughout; `../mfem-hdg-dev-par` exists for parallel work |
| **`MFEM_USE_CUDA = YES`**, `CUDA_ARCH = sm_75` | CUDA 13.3. **`sm_75` because the card is compute capability 7.5**; MFEM's cache defaulted to `sm_60`, which is wrong here and nothing warns you |
| **`MFEM_USE_CUDSS = YES`** | cuDSS 0.8.0 from Debian (`/usr/include/cudss.h`), found with `-DCUDSS_DIR=/usr`. `mfem::CuDSSSolver` agrees with UMFPACK to 3.5e-14; its **timings on this machine are not reproducible** and meq's solver does not use it. See *Threading, measured* |
| `RAJA`, `OCCA`, `CEED`, `SIMD`, `AMGX` all `NO` | |

**Turning CUDA on is not one flag, and two of the three obstacles are silent.**

* **SUNDIALS must itself be a CUDA build.** `linalg/sundials.hpp` hard-errors
  with `MFEM_USE_CUDA=TRUE requires SUNDIALS to be built with CUDA support`,
  testing `SUNDIALS_NVECTOR_CUDA`. meq cannot drop SUNDIALS — `KINSolver` backs
  `AndersonPicard` and `PicardThenNewton` — so the answer is a second SUNDIALS.
  `../sundials/cuda-install` is 7.5.0 configured `-DENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=75`, otherwise matching the old one.
* **~~CMake will not generate at all with `MFEM_ENABLE_TESTING=ON` in a serial
  CUDA build~~ — FIXED UPSTREAM, and `MFEM_ENABLE_TESTING` is now ON.** The bug
  was that a device test loop merged the serial and parallel example lists and
  never re-split them, registering a test against `sundials_ex9p` in a build
  that never created it. Both loops — `examples/CMakeLists.txt` and
  `examples/sundials/CMakeLists.txt` — now guard the parallel branch with
  `elseif (MFEM_USE_MPI)`, so a serial build skips it. **Verified 2026-08-31**:
  reconfigure, full rebuild and install all returned 0, and ctest registers
  **131 tests**. `MFEM_ENABLE_TESTING` does not appear in the installed
  `config.mk`, so it changes nothing meq's finder reads.

  **This file claimed the bug was "filed as
  `../mfem-hdg-dev/doc/CMAKE-SERIAL-CUDA-SUNDIALS.md`". IT WAS NOT.** That path
  has never existed in that tree's history and exists nowhere on disk — the
  report was never written. Treat a "filed as" claim in this file as needing a
  check, not as evidence.
* **Stale derived cache variables survive a `-DSUNDIALS_DIR=` change.**
  `SUNDIALS_INCLUDE_DIR` and friends are separate cache entries, so repointing
  the directory leaves the old include path on `TPL_INCLUDE_DIRS` *in front of
  or behind* the new one, and which wins is an ordering accident. Delete every
  `SUNDIALS_*` entry but `SUNDIALS_DIR` and reconfigure.

**And a CUDA build broke meq's own `FindMFEM.cmake`, which is worth knowing
because every serial build hides it.** MFEM puts the toolkit headers on the line
as `-isystem <dir>` inside `MFEM_CXXFLAGS`, not as `-I<dir>` in
`MFEM_TPLFLAGS`, while `general/backends.hpp` includes `<cuda_runtime.h>`
unconditionally once `MFEM_USE_CUDA` is defined. meq's finder took only
`TPLFLAGS`, so the first `#include <mfem.hpp>` failed on a missing CUDA header —
a confusing place to meet a build-configuration problem. The finder now also
parses `-isystem` out of `MFEM_CXXFLAGS`.

**What CUDA does NOT buy meq today**, so nobody reads the flag as a capability:
`fem/darcy/` still has no partial-assembly kernels, meq uses no `mfem::forall`,
and meq's own sources still compile with `g++` rather than `nvcc`. What it
enables is `BatchedLinAlg`'s `gpu_blas`/`magma` backends and `CuDSSSolver`. It
is a prerequisite, not a speedup.

**`../mfem-hdg-dev` is now the development tree and not what meq links.** That
split exists because the alternative was demonstrated: mid-session that tree was
switched to `gf-hdg-dev`, which has no `fem/darcy/extension_hdg.{hpp,cpp}`, and
meq's `naming` check failed on a tree meq had not touched. A library meq depends
on should not move under it while somebody is working on it. Point `MFEM_DIR` at
the dev tree deliberately when testing a fix; otherwise leave it alone.

### `meq-integration`: the branch meq builds from, and why it is local only

meq needs work from **four MFEM branches**, and as of 2026-08-30 it takes
**two merges** rather than three to get them, because the topology changed:

| branch | what meq needs from it | |
|---|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — stage 5's curved `Γ`, and `TransferredDatumCoefficient` | merge base |
| `direct-solver-symbolic-reuse` | `UMFPackSolver` and `PardisoSolver` keeping their symbolic factorisation across Newton steps | merge 1 |
| `gf-hdg-linearise-first` | **`DarcyNPCOperator` / `DarcyNPCSolver`** — the NPC method, which is what meq's default ordering now is — plus `SetAssemblyMode`, `SetGradientMode`, `SetLocalFactorMode`. **`SetNonlinearOrdering` was on this list and is DELETED**; see *The NPC port* | merge 2 |
| `gf-hdg-dev` | ***"The postprocessing closes on the element average, always"*** — the reconstruction fix, without which `ψ*` is a different function wherever `∂F/∂ψ` vanishes | **now free** |

**`gf-hdg-dev` is now an ANCESTOR of both `gf-hdg-subdomains-dev` and
`gf-hdg-linearise-first`** — verified pairwise 2026-08-30 — so it arrives with
the base and needs no merge of its own. **This file previously said the four
branches were pairwise independent and that is no longer true.** Re-verify
before trusting either statement: both were true when written, and the topology
is the other tree's to change.

What has not changed is *why* `gf-hdg-dev` is in the list. Dropping it gives a
silently degraded `ψ*` on any element where the Jacobian's reaction term
vanishes, an adaptive loop a full order worse, and a suite that stays green
until `thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` runs. It is
free today rather than unnecessary.

`meq-integration` is their merge, and it exists **only** to be built against. It
is not a development branch and must not be pushed anywhere: it will be
re-created whenever any of them moves, so anything committed directly to it is
lost. Do the work on the branch it belongs to and re-merge.

**It lives in `../mfem/mfem-src`, which is meq's own clone — never in
`../mfem-hdg-dev`.** See the rule below.

**Expect the fetch to be a forced update.** Two of the four were rebased between
2026-08-29 and 08-30, so a merge *into* the existing `meq-integration` is not
the operation you want; re-create it. Tag the old one first — the tree carries
`meq-integration-before-2026-08-30` for exactly that reason — because a bad
re-merge is otherwise unrecoverable without redoing the conflict resolution.

Refreshing it, from `../mfem/mfem-src`:

```sh
git fetch hdgdev gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
          gf-hdg-linearise-first gf-hdg-dev
git tag -f meq-integration-before-$(date +%F) meq-integration
git checkout -B meq-integration hdgdev/gf-hdg-subdomains-dev
git merge hdgdev/direct-solver-symbolic-reuse     # 3 conflicts; see below
git merge hdgdev/gf-hdg-linearise-first           # 3 conflicts; see below
# gf-hdg-dev needs no merge -- verify with:
for b in gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
         gf-hdg-linearise-first gf-hdg-dev; do
  git merge-base --is-ancestor hdgdev/$b HEAD && echo "$b contained"
done
```

**Run that last loop.** It is four cheap commands and it is the only thing that
catches a branch having moved out from under the recipe — which is precisely
what happened to the three-branch version of this section.

**BUT THE LOOP CHECKS THE BRANCH, NOT THE INSTALL, AND THAT GAP HAS ALREADY
BITTEN.** Every command in it runs in `../mfem/mfem-src` and says nothing
whatever about `../mfem/install`, which is what meq actually links. On
2026-09-01 the loop reported **all four contained** while the installed library
was `fa65a2f932` and `meq-integration` had moved to `d244329d8c` — **514 changed
lines in `fem/darcy/` across seven files**, including *"The HDG face quadrature
never saw the trace element"* and *"DarcyOperator dereferenced a null
prolongation on a hanging-node-free NC mesh"*, the second in exactly the
non-conforming meshes the adaptive loop makes. So the one check this project
trusts passed on a tree whose measurements were all taken against a library a
day out of date.

**Check the install too, and it is one command:**

```sh
diff <( git -C ../mfem/mfem-src show meq-integration:fem/darcy/darcyhybridization.hpp ) \
     ../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp \
  && echo "install is current"
```

A branch that another agent re-creates is a branch that moves without meq doing
anything, so this belongs beside the containment loop rather than in a rebuild
checklist nobody reads.

**And it has already fired again. Run 2026-09-01 against `meq-integration` at
`fa65a2f932`: three of the four report contained and `gf-hdg-linearise-first`
does not**, because that branch has advanced **17 commits** since the merge.
**Read a "not contained" as *the branch has moved*, not as *the merge failed*** —
what meq needs from it is in the tree (`DarcyNPCOperator` is in the installed
`fem/darcy/darcyhybridization.hpp`, and NPC is what meq runs on). It says a
re-merge is available, not that the build is wrong. Re-create rather than merge
into the existing branch, per the recipe above.

**The second merge now conflicts too, and it did not used to.** This file used
to record `gf-hdg-linearise-first` as merging clean; after the 08-30 rebase it
brings three, and all three are decided rather than fiddly:

| file | resolution |
|---|---|
| `doc/HDG-DEFECTS-FROM-MEQ.md` | modify/delete. **Take the delete** — it is meq's own retired bug report, retired deliberately by *"Retire two bug reports whose findings are all fixed and covered"* |
| `doc/HDG-ROADMAP.md` | two regions, large and disjoint. Keep both sides |
| `miniapps/hdg/makefile` | a union: `extension` from one side, `navierstokes` from the other. Merge the `SEQ_MINIAPPS` line into one and keep both test targets |

Only the makefile affects a build, and only of the miniapps. **`HDG-ROADMAP.md`
has exactly two lines that are `=======` and exactly two conflict starts**, so
the marker is unambiguous there — but check that before scripting it, for the
`CHANGELOG` reason below.

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

### `../mfem-hdg-dev` is receive-only

**Write MFEM development *requests* into `../mfem-hdg-dev/doc/` and nothing
else.** No branches, no commits, no checkouts, no builds. Another agent works in
that tree and it is theirs; meq consumes it through `git fetch` from
`../mfem/mfem-src` and in no other way.

This is written down because it was got wrong: `meq-integration` was first
created *in* the dev tree, given two merge commits and a fix commit, and that
tree was left checked out on it — so the next person to touch it would have been
on someone else's branch. Fetching from the dev tree is reading; everything else
is not.

**WHICH REQUESTS EXIST IS A QUESTION FOR `git`, NOT FOR `ls`, AND THIS FILE GOT
IT WRONG EARLIER THE SAME DAY.** A listing of `../mfem-hdg-dev/doc/` reports
whatever branch that tree happens to be checked out on, which is not meq's to
control and changed twice on 2026-09-01. Asked properly —
`git cat-file -e <branch>:doc/<file>` across every branch — the answer is:

| document | lives on | |
|---|---|---|
| `HDG-ELEMENT-LOCAL-PARALLELISM.md` | `gf-hdg-linearise-first` | **open** |
| `HDG-BEM-COUPLING-FROM-MEQ.md` | `gf-hdg-linearise-first` | **open** |
| `HDG-NPC-GLOBALISATION-FROM-MEQ.md` | `gf-hdg-linearise-first` | **open**, and answered in place |
| `HDG-DEFECTS-FROM-MEQ.md` | **`gf-hdg-dev` and `gf-hdg-subdomains-dev`** | **NOT retired** — deleted only on the symbolic-reuse line, which is exactly why that merge conflicts modify/delete |
| `HDG-LINEARISE-THEN-CONDENSE.md` | backup refs only | retired with the mode |
| `DIRECT-SOLVER-SYMBOLIC-REUSE.md` | no branch at all | retired |

All three open ones sit on **one** branch, so with the tree on `gf-hdg-dev` they
are invisible and `doc/` looks nearly empty. An earlier version of this
paragraph read that emptiness as three retirements, and said
`HDG-DEFECTS-FROM-MEQ.md` had been retired in `2a50119ba1` — that commit deletes
it on one line of history and the file is alive on two of the four branches meq
merges. **The merge table below already knew this**, and recording the two
against each other is the point: a claim about someone else's working tree is a
claim about their checkout.

**`HDG-NPC-GLOBALISATION-FROM-MEQ.md`, 2026-08-31 — FILED, AND ANSWERED THE
SAME DAY** (`af82d42b14` in that tree). Deliberately **not** a defect report:
it disputed §6 of `HDG-ORDERING-API.md`, whose recommended backtracking line
search meq had implemented from `miniapps/hdg/navierstokes.cpp` and measured
making every case worse, and it asked what configuration produced §6's numbers
rather than asserting they were wrong. Upstream withdrew two §6 claims on meq's
evidence, corrected meq's own account of the mechanism, and found a defect in
their reference implementation that meq had inherited by copying it faithfully.
All of it is recorded under *Why it fails, measured*. **That two of meq's claims
were corrected by the exchange rather than by meq is the argument for writing
these notes at all.**

**And do not file findings against unfinished work.** A branch that exists is not
a branch that is done. Measure it if it is useful to know, keep the numbers in
meq, and wait to be told it is ready before writing anything into that tree about
how it behaves.

Rebuilding the install, after fetching whatever is wanted into `../mfem/mfem-src`:

```sh
cd /home/ian/projects/mfem
cmake --build build -j4 && cmake --install build
```

The finder reads `share/mfem/config.mk` from an install and `config/config.mk`
from an in-source make build, so either layout works and neither needs a flag.

`mfem/master` will not do. meq needs `DarcyForm` and the HDG integrators in
`fem/darcy/`, and for stage 5 it needs the transfer-path machinery in
`fem/darcy/extension_hdg.{hpp,cpp}` — which is an implementation of exactly the
curved-boundary technique the GS papers use.

**The old MFEM's HDG API is gone, and that is why this is a port and not a
recompile.** `HDGBilinearForm`, `HDGDomainIntegratorGS`, `HDGFaceIntegratorGS`,
`AssembleSC` and `Reconstruct` have no counterparts; the replacement is
`DarcyForm` with `EnableHybridization`. The mapping, term for term:

| Weak term | Old | New |
|---|---|---|
| `(r q_h, v)` | `HDGDomainIntegratorGS` | flux mass: `VectorMassIntegrator` with an `r` coefficient |
| `(ψ_h, ∇̄·v)`, `(q_h, ∇̄w)` | `HDGDomainIntegratorGS` | flux divergence: `VectorDivergenceIntegrator` |
| `⟨ψ̂_h, v·n⟩` | `HDGFaceIntegratorGS` | the transpose of `NormalTraceJumpIntegrator` — **not** the face integrators on `B`, see below |
| `⟨τ(ψ_h − ψ̂_h), w⟩` | `HDGFaceIntegratorGS` | potential mass: `HDGDiffusionIntegrator` |
| `⟨q̂_h·n, μ⟩ = 0` | `AssembleSC` | `EnableHybridization(trace_space, new NormalTraceJumpIntegrator(), ess_list)` |
| condense / reconstruct | `AssembleSC` + `Reconstruct` | `FormLinearSystem` / `RecoverFEMSolution` |
| `φ_h` on `Γ_h` (stage 5) | — | `HDGExtensionIntegrator` on the **flux mass form**, `C = r`, sign `+1` |

`miniapps/hdg/convdiff.cpp` in that tree is the worked example to copy from;
lines 445–469 build exactly the three spaces above.

**A correction to that table, measured in stage 2.** An earlier version said the
`⟨ψ̂_h, v·n⟩` coupling came from `TransposeIntegrator(DGNormalTraceIntegrator)`
added to `B` on interior and boundary faces. It does not. Under hybridization
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


### Post-processing is back, and it was free

Stage 3 dropped `ψ*` on the grounds that `q` already gives the physical quantity
at `k+1`, and noted that stage 6's estimator would need it again. It does — eq.
(20) uses `ψ*` in **four of its five terms** (`η₁`, `η₂`, `η₄`, `η₅`).

**`DarcyForm::Reconstruct()` supplies it, measured at `k+2`**: rates 3.03, 4.03,
5.00 for `k = 1, 2, 3`, and it survives the extension path and Newton. No
hand-written local solve was needed and the old `GSSolver::Postprocess()` stays
deleted. `thePostProcessedPotentialSurvivesNewton` is that measurement: MFEM
lifts the non-linear potential integrators as a Jacobian frozen at the computed
potential, which on Example 5 delivers 3.05, 4.05, 5.03 and is 47×, 113×, 125×
smaller than `ψ_h` on the finest mesh.

The measurement that justifies *requiring* `ψ*` rather than quoting it: building
`η₂` on raw `ψ_h` loses **exactly one order at every `k`** — 2.002 vs 0.998,
3.000 vs 2.002, 3.992 vs 2.981 — and is 124× to 407× larger on the finest mesh.
That is the defect the pre-modernisation estimator had.

**IT WAS BROKEN, SILENTLY AND PER ELEMENT, AND THE SHAPE OF THAT IS THE PART
WORTH KEEPING.** The local post-processing problem is a pure **Neumann** problem
by construction — the total flux driving it is in H(div), so the element's flux
balance is already satisfied and the potential is determined only up to a
constant. The element average is what closes it, and that is *unconditional*.
MFEM had been treating it as one of three choices: a flag set on the mere
presence of a non-linear integrator skipped the mean-value branch, and a
singular matrix was then factored anyway with `Factor( m, TOL = 0.0 )` testing
`|pivot| <= 0.0` and the return value discarded. Where `∂F/∂ψ` vanishes the
integrator's Jacobian *is* the zero matrix, so the term **is** singular.

**The flag is per element, so the corruption was too**, and that is what made it
invisible. With `∂F/∂ψ` vanishing on an eighth of the domain — a tabulated
profile with a flat segment, which is ordinary — elements were **20× wrong**
while the whole-domain norm read 1.87; at half the domain, 61× against 7.57.
**Any global check misses it.** And a floor of **1e-12** on `∂F/∂ψ` — twelve
orders below everything else in the problem, incapable of moving a solution —
took the ratio from 7.565317 to 0.998525. That is a singular matrix and nothing
else.

The fix is `gf-hdg-dev`'s *"The postprocessing closes on the element average,
always"*. Measured from meq's side on the case that used to read 20.3, 64.1 and
61.6 on dead elements: **1.0069, against 1.0048 where `∂F/∂ψ` does not vanish at
all.** The driver's adaptive loop is back on the published estimator — `η`
reaching 6.87e-5 on **449** elements where the degraded one needed 1069 to reach
4.77e-4.

**meq carried no runtime check for this, deliberately.** A solver should not
stand permanently on guard against its dependency. The state of the defect lived
in the suite instead, as a test asserting the *correct* behaviour and failing.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` still runs, now as
a regression rather than a tripwire.

**And the defect did not reach the solve — measured, not argued.** The forward
solve's local problem takes its potential block from the HDG stabilisation on
`Mnl_p`'s **interior faces**, a fixed bilinear form whatever the source does,
while `Reconstruct()` builds a different local problem whose regularisation is
what failed. `theReconstructionDefectDoesNotReachTheSolve` compares both paths
in `ψ_h` **and** `q_h` over `k = 1…4` and three meshes: worst relative
difference **1.6e-13** and **2.6e-13**, against discretisation errors from 2.2e-3
down to 7.3e-11. It also checks that post-processing leaves both bit-identical,
since the driver writes them afterwards.

**Two things to keep.** The fix is on `gf-hdg-dev` and on no other branch, so a
`meq-integration` rebuilt without it silently loses this — see the merge recipe.
And the failure mode is the transferable part: **a per-element defect that a
global norm cannot see, in a quantity four of the estimator's five terms are
built on.**

## The NPC port, and the parity gap that came back the other way

**MFEM deleted `NLOrdering::LineariseThenCondense` on 2026-08-31**, as *"a
condensation in disguise"*: it was an `Operator` on the trace alone that kept the
linearisation as hidden state, and NPC's fields are Newton state, which a
trace-only operator has nowhere to put. `SetNonlinearOrdering()` went with it.
That was a breaking change for meq, whose default ordering it was.

**meq's default is now `NonlinearOrdering::NPC`** — `mfem::DarcyNPCOperator`
with `mfem::DarcyNPCSolver`, Newton on the full `(q, ψ, ψ̂)` system with the
Jacobian solved by hybridized elimination, Nguyen, Peraire & Cockburn eqs
(14)–(18). `CondenseThenLinearise` is kept as the **backup**, and the reference
is `../mfem-hdg-dev/doc/HDG-ORDERING-API.md`.

**The port cost meq almost nothing structurally, and the reason is an accident
worth knowing.** `solution` has always been a three-block
`{ flux, potential, trace }` `BlockVector` on a four-entry `blockOffsets`, with
`darcyFlux`, `potentialGf` and `traceGf` all `MakeRef`'d into it. **That is the
NPC unknown, block for block.** So the fields are already in place when the
iteration returns, and `RecoverFEMSolution()` leaves the Newton path rather than
needing rework. `meq::Relinearised` — the wrapper that existed only to pair the
residual with the gradient under the deleted mode — was deleted with it.

### What NPC bought, measured

**`GetNumLocalNLIterations()` is exactly zero**, which is the acceptance signal
that this is NPC and not a condensation wearing the name.
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally` asserts it on
Example 5 and prints the comparison:

| `k` | NPC local NL its | condense local NL its | L2 vs exact | relative in `ψ` |
|---|---|---|---|---|
| 1 | **0** | 3644 | 7.614902e-03 both | 2.9e-15 |
| 2 | **0** | 3560 | 2.494695e-04 both | 1.6e-11 |
| 3 | **0** | 3412 | 6.000146e-06 both | 9.0e-12 |

Same L2 to seven figures, same 4 Newton iterations, ψ agreeing to round-off.
That is what says the backup is a backup and not a different problem.

**AND AN NPC ITERATION IS THREE TO FOUR TIMES CHEAPER THAN A CONDENSATION
ONE**, which is what the iteration-count tables everywhere else in this file do
*not* show and is where the suite's 872 s → 499 s actually came from. Measured
2026-08-31 on two cases where both orderings converge and agree:

| | NPC | CondenseThenLinearise | |
|---|---|---|---|
| §4.2 pedestal `k = 2, n = 24` | **9 its, 0.92 s**, 0 local NL | 10 its, **3.93 s**, 75,490 local NL | agree to 4.1e-11 |
| similarity `k = 2, n = 16` | **4 its, 0.19 s**, 0 local NL | 4 its, **0.70 s**, 10,898 local NL | agree to 3.6e-11 |

**4.3× and 3.7× on wall clock at the same or fewer iterations.** So reading the
parity table as "NPC needs 17 where the condensation needs 11" understates NPC
by most of a factor of four; the honest currency here is seconds, and the local
non-linear iteration counts in the third column are what NPC deletes.

**The bordered Newton got simpler and stronger, and this is the real prize.**
With `ψ` an unknown rather than a function of the trace, two of the three
bordered quantities stop being finite differences:

| | condensation | **NPC** |
|---|---|---|
| `b = −∂(max ψ_h)/∂x` | `3(k+1)` central differences over one element's trace dofs | **exactly `−e_j`**, one entry, not differenced |
| `d = ∂G/∂s` | `1 −` a central difference | **exactly `1`** |
| `c = ∂R/∂s` | one central difference in a scalar | the same |

and the whole frozen-seed apparatus goes away — `formSystem()` is no longer
re-called once per accepted step, because there is no element-local non-linear
solve whose initial guess could go stale. `HighBetaConvergence` reproduces
**every `ψ_ax` in the table under *The measurement* to every digit printed**, and
the constraint is now satisfied to machine zero:

| `ν` | `A` | `ψ_ax` | `ψ_ax − max ψ_h`, was | now | Newton, was | now |
|---|---|---|---|---|---|---|
| 2 | 1 | 3.058984e-01 | 1.7e-16 | 0.0 | 4 | 4 |
| 2 | 10 | 9.607537e-01 | −4.4e-16 | 0.0 | 4 | 4 |
| 2 | 100 | 3.036075e+00 | 4.4e-16 | 0.0 | 4 | 4 |
| 4 | 1 | 2.834510e-01 | **3.1e-12** | −5.6e-17 | 8 | **6** |
| 4 | 10 | 8.643745e-01 | 8.9e-16 | 0.0 | 11 | **7** |
| 4 | 100 | 2.720222e+00 | **1.0e-10** | 0.0 | 10 | **7** |

`HighBetaConvergence` went **17 s → 3.3 s**. The printed residual is a different
quantity — `‖(R, γG)‖` over the full system rather than the trace — so the
histories are not comparable across the port and `γ` is a different number; the
converged `ψ_ax` is what is comparable, and it is unchanged.

### What NPC cost, and it is the parity gap in the other direction

**THIS IS THE FINDING TO READ BEFORE ASSUMING THE PORT WAS FREE.** The
expectation going in was that NPC would remove the parity gap the deleted mode
had. It does not. It replaces it with a mirror image, and a wider one. Measured
2026-08-31, raw undamped Newton, cap 60:

**WIDENED 2026-08-31 (later) FROM 7 CASES TO 14, AND ONE SENTENCE OF IT WAS
WRONG.** The seven-case version said NPC was *never* better on any case
measured. It is better on two, and they are the same two corner: well resolved,
higher order.

| case | **NPC** | CondenseThenLinearise |
|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | **fails at 60** | ok, 14 |
| §4.2 pedestal `k = 1, n = 24` | ok, 26 | ok, 23 |
| §4.2 pedestal `k = 1, n = 32` | ok, 10 | ok, 9 |
| §4.2 pedestal `k = 2, n = 16` | ok, 17 | ok, 11 |
| §4.2 pedestal `k = 2, n = 24` | ok, **9** | ok, 10 |
| §4.2 pedestal `k = 3, n = 16` | ok, **8** | ok, 12 |
| §4.5 layer `k = 1, n = 30` | **fails at 60** | ok, 13 |
| §4.5 layer `k = 1, n = 60` | ok, 10 | ok, 10 |
| §4.5 layer `k = 2, n = 16` | ok, **51** | ok, 10 |
| §4.5 layer `k = 2, n = 24` | ok, 13 | ok, 12 |
| §4.3 barrier `k = 1, n = 16` | fails at 60 | **fails at 60** |
| §4.3 barrier `k = 2, n = 16` | fails at 60 | **fails at 60** |
| similarity solution `k = 2, n = 16` | ok, 4 | ok, 4 |
| similarity solution `k = 3, n = 16` | ok, 4 | ok, 4 |

**The pattern is resolution, not stiffness.** NPC loses where the mesh is too
coarse for the source and wins where it is not — `k = 2, n = 24` and
`k = 3, n = 16` go to NPC. And since an NPC step carries **no element-local
non-linear solve at all**, an equal iteration count is already an NPC win on
wall clock.

**The §4.3 rows are new and they matter, because they are the control.** The
barrier defeats *both* orderings at `k = 1` and `k = 2` on `n = 16`, so this is
not a table about one ordering being sound and the other not. (These rows are
run **without** `setInitialGuess()`, unlike the §4.3 rows under *Picard, then
Newton*, which seed the ramp — that is why `k = 2` reads a failure here and
`ok, 10` there. The guess is doing the work in that table, and it is worth not
confusing the two.)

**AND WHEN NPC FAILS IT FAILS CHEAPLY, WHICH IS NOT A SMALL THING.** Same
barrier case, `k = 1, n = 16`, both reaching the cap of 60:

| | wall clock | element-local non-linear iterations |
|---|---|---|
| **NPC** | **1.6 s** | **0** |
| CondenseThenLinearise | **53.5 s** | **2,027,130** |

A factor of 33 in the cost of finding out that neither works. That is exactly
the uniformity §6 promises, and it holds on an iteration going nowhere as much
as on one that is working — which makes NPC the better thing to *try first*
even on the cases it loses.

That NPC is not a free win is not a surprise to upstream, whose own §6
says *"Reach for [the reduced trace operator] unless you have a reason not to"*
and *"**NPC is not automatically faster.** Its advantage is uniformity of the
local work … not fewer floating-point operations."*

### MFEM's own test suite reproduces the parity gap, and two of its tests are red

**`MFEM_ENABLE_TESTING` IS NOW ON, and turning it on found this in the first
run.** meq had been building with `-DMFEM_ENABLE_TESTING=OFF` because of a CMake
bug in serial CUDA builds; that bug is fixed (see *Which MFEM*), so MFEM's own
131 registered tests are available. Run over the Darcy and HDG tags on
`meq-integration` at `974871c456`:

```
MKL_NUM_THREADS=1 \
  ./tests/unit/unit_tests "[DarcyForm],[NPC],[HDG],[DarcyHybridization],[NonlinearDarcy]"
```

**72 of 74 test cases pass; 98,434 of 98,436 assertions.** Both failures are in
`tests/unit/fem/test_darcy_npc.cpp`, and both are NPC on a stiff pedestal:

| test | what fails |
|---|---|
| *"A stiff source converges by condensation and by NPC alike"* | the condensation section **passes**, the NPC section **fails** |
| *"NPC solves stiff problems LineariseThenCondense cannot"* | fails outright |

Both at `n = 32`, order 1, `σ² = 0.003`, stalling at **1.1590e-03** after 41
iterations with `local_nl_iters = 0` — so it really is NPC, and the residual is
sitting at 1e-3 exactly as meq's §4.5 does at its own `n = 32`. **And the call is
`RunNPC( P, 40, true, GM::Assembled )`, where the `true` is the backtracking line
search** — upstream's recommended globalisation, on, and still failing.

**This is meq's parity gap reproduced inside MFEM's own suite**, on their fixture
rather than meq's, in two tests whose *names* encode the claim it contradicts. It
is independent evidence for everything in the two sections above, and it was
free: it took one flag and one run.

**Do not read it as a regression report.** meq has never run these tests before —
testing was off until 2026-08-31 — so there is no earlier green to compare
against, and `974871c456` predates upstream's own response commit. What it
establishes is that the finding is not about meq's discretisation, since it
reproduces on a fixture meq did not write.

### Why it fails, measured — and why a line search does not fix it

**IT IS A PROPERTY OF THE TWO ALGORITHMS, NOT A DEFECT IN EITHER, AND THE
MECHANISM IS MEASURED RATHER THAN ARGUED.**

Both orderings take a bad first step — §4.2 at `k = 1, n = 16` goes 1.14e-01 →
6.25e+00 under NPC and 7.78e+00 → 5.57e+01 under the condensation. The
condensation recovers and reaches its quadratic endgame in 14; NPC wanders
non-monotonically between 1e-2 and 1e+1 for sixty steps and never enters a
basin. Not a stall, not divergence to NaN — undamped Newton outside its basin.

**Where the Newton remainder goes, per block.** At the cold iterate
`x = [0, 0, 2.771]` — flux and potential zero, trace carrying the Dirichlet
datum — one full step gives

| | flux | potential | trace |
|---|---|---|---|
| residual at `x` | 7.93e-02 | 8.13e-02 | **0.00e+00** |
| residual after a full step | **1.14e-16** | **6.25e+00** | 1.82e-13 |
| the correction `c` | **3.50e+02** | **1.51e+02** | 1.48e+01 |

**The flux and trace rows of the NPC residual are LINEAR in `(q, ψ, ψ̂)`** — the
flux equation and the flux constraint carry no `F` — so a Newton step annihilates
them *exactly*, to 1e-16 and 1e-13. **All the non-linearity is in the potential
row**, and so is the entire Newton remainder: 8.1e-02 → 6.25e+00. And the
correction is `O(10²)` against a solution where `ψ ~ 0.3`, because at `ψ ≡ 0` the
linearised operator `−∇̄·((1/r)∇̄·) − (∂F/∂ψ)/r` is at its most indefinite — the
pedestal sits at `max|∂F/∂ψ|/λ₁ ≈ 7`.

**What the condensation does instead is nonlinear elimination, and that is a
preconditioner.** Its element-local solves take `(q, ψ)` to the *exact* solution
of each element's non-linear problem at the current trace, so the fields never
take an unphysical linear excursion and the outer residual is a tamer function of
a much smaller unknown. That is a real and well-understood robustness mechanism,
and it is precisely what NPC gives up in exchange for uniform local work.

**A LINE SEARCH ON THE FULL RESIDUAL DOES NOT FIX IT, AND THE TABLE ABOVE SAYS
WHY.** Upstream recommends backtracking on the full residual as NPC's
globalisation — `NewtonSolver::ComputeScalingFactor`, a dozen lines,
`miniapps/hdg/navierstokes.cpp` carries one. **It was implemented in meq and
measured, and it made every case worse**, including the five that converge
undamped: all of them went to sixty iterations, creeping by about 1% a step.

**THE MECHANISM OF THE CREEP IS A DEFECT IN THE REFERENCE IMPLEMENTATION, AND
UPSTREAM FOUND IT AFTER MEQ REPORTED THE SYMPTOM.** meq's first reading — that
`α` collapses because the `ℓ²` merit charges for restoring `(1 − α)` of the flux
and trace residuals — describes the pressure correctly and gets the mechanism
wrong. `NSBacktrackingNewton` accepts on `Norm(rt) < n0`, **a monotone test with
no sufficient-decrease constant**. For Newton on an `ℓ²` merit the direction is
always a descent direction, so `merit(α) ≈ merit(0)(1 − α)` and *any* small
enough `α` passes: upstream swept it and `α = 1.2e-4` still "succeeds",
improving the merit by 1e-4 relative. **So where `α = 1` is rejected the search
does not fail — it creeps**, which is exactly the 1%-a-step crawl meq measured.
An Armijo test would reject those steps and fail honestly. Fixing it is open
upstream, and **it would not make meq's problem converge**; it would stop the
globalisation disguising its own failure. **KINSOL's `KIN_LINESEARCH` does use
an Armijo condition and fails on the same cases**, which is why the two rows
agree.

**AND THE BLOCK STRUCTURE IS NOT WHAT SEPARATES MEQ'S PROBLEM FROM UPSTREAM'S —
SEVERITY IS.** This file claimed the failure follows from two rows being linear
and one non-linear. Upstream measured the same shape on their own Darcy pedestal
— a full step takes the flux row to 4e-17 and the trace row to 3e-14 with the
whole remainder in the potential row — and there the line search *works*: 3 of 4
configurations fail undamped and converge with backtracking. **What differs is
that a full step IMPROVES their potential residual, 4.7e-01 → 1.7e-01, and
MULTIPLIES meq's BY 77.** The block structure is common to both; the size of the
first correction is not. So the recommendation in their §6 is supported there
and wrong here, and meq's original explanation was over-general.

meq **asked whether §6's baseline was the deleted mode, and it was not** — it
was plain undamped `NewtonSolver` on NPC, so the disagreement is real rather
than an artefact of what was being compared against. On the strength of meq's
measurements upstream has withdrawn two §6 claims: *"NPC is not automatically
faster … not fewer floating-point operations"*, and *"reach for [the reduced
trace operator] unless you have a reason not to"* as a general default. §6 now
says the choice turns on what the element-local non-linear solve costs on your
problem.

So the line search was **not kept**, and the comment at the `Globalisation::None`
switch in `GradShafranov.cpp` records why. What a step-length rule would need
here is to damp the non-linear block without un-doing the exact annihilation of
the linear ones — different, and unattempted.

**What DOES work is fixing WHERE the iterate is rather than how far it steps**,
which is the reactive ladder meq already has, and this also explains why: Picard
hands NPC a physically sensible `(q, ψ, ψ̂)` instead of `(0, 0, g_D)`, and the
huge first correction never arises. It is the same observation as the warm-start
one below — **under NPC the fields are state, so what they start at matters**,
and `(0, 0)` is the worst place to start.

**Both failures are cured by the reactive ladder, which is what the driver
runs**, so nothing in production is affected:

| | NPC, plain | NPC, `PicardThenNewton` |
|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | fails at 60 | **ok, 5** |
| §4.5 layer `k = 1, n = 30` | fails at 60 | **ok, 5** |

and both are cured by refinement as well, which is the standing reading of these
sources: *a difficulty measured at one resolution is not a property of the
problem*.

### There is no cheap discriminator, and the obvious one is ANTI-correlated

**ASKED AND ANSWERED BY MEASUREMENT, 2026-08-31.** The reactive ladder pays for
a failure by running Newton to its cap before handing off. The obvious
improvement is to notice earlier — and the obvious signal, the first Newton
step making the residual worse, **is not a signal at all**. NPC, plain Newton,
cap 60, cold start:

| case | out | its | `‖r₀‖` | `‖r₁‖/‖r₀‖` |
|---|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | **FAIL** | 60 | 1.136e-01 | 5.50e+01 |
| §4.2 pedestal `k = 1, n = 24` | ok | 26 | 9.242e-02 | **1.17e+03** |
| §4.2 pedestal `k = 1, n = 32` | ok | 10 | 7.988e-02 | 2.28e+02 |
| §4.2 pedestal `k = 2, n = 16` | ok | 17 | 1.124e-01 | 2.91e+02 |
| §4.2 pedestal `k = 2, n = 24` | ok | 9 | 9.172e-02 | 5.04e+01 |
| §4.2 pedestal `k = 3, n = 16` | ok | 8 | 9.542e-02 | 5.72e+01 |
| §4.5 layer `k = 1, n = 30` | **FAIL** | 60 | 8.253e-02 | 3.79e+02 |
| §4.5 layer `k = 1, n = 60` | ok | 10 | 5.817e-02 | 3.55e+01 |
| §4.5 layer `k = 2, n = 16` | ok | 51 | 1.124e-01 | 2.91e+02 |
| §4.5 layer `k = 2, n = 24` | ok | 13 | 9.172e-02 | 5.04e+01 |
| §4.3 barrier `k = 1, n = 16` | **FAIL** | 60 | 1.136e-01 | **8.56e-01** |
| §4.3 barrier `k = 2, n = 16` | **FAIL** | 60 | 1.124e-01 | **8.67e-01** |
| similarity `k = 2, n = 16` | ok | 4 | 1.152e-01 | 5.18e-02 |
| similarity `k = 3, n = 16` | ok | 4 | 9.784e-02 | 5.23e-02 |

**The largest first-step blow-up in the table converges, and the two cases
where the residual came DOWN at step one both fail.** 1169× at
`k = 1, n = 24` finishes in 26; the barrier's 0.86 and 0.87 run to the cap.
Any threshold that catches the pedestal at `n = 16` throws away `n = 24`, and
no threshold whatever catches the barrier. Under the condensation it is the
same story from the other side — the barrier drops to 4.3e-02 and 1.1e-02 at
step one and still fails.

**And `‖r₀‖` carries no information either, because under NPC it is very nearly
source-INDEPENDENT.** Look down that column: the pedestal and the barrier at
`k = 1, n = 16` both read 1.136e-01, and the pedestal and the layer at
`k = 2, n = 16` both read 1.124e-01 — to four figures, on different sources. At
the cold iterate `(0, 0, g_D)` the residual is dominated by the Dirichlet datum
entering the flux row through `⟨ψ̂, v·n⟩` and the potential row through the
stabilisation, and `(F/r, w)` is a small correction to both. `‖r₀‖` is a
measurement of the boundary data and the mesh, not of the difficulty.

**So the ladder stays reactive and stays triggered by observed failure**, which
is what *Why meq's Newton struggles* already concluded from `max|∂F/∂ψ|/λ₁`.
Two independent candidate detectors have now been measured and both are
useless; the difference is that this one is *anti*-correlated, which is worse
than uninformative. **Do not add a predictive trigger without a measurement
that separates these fourteen rows**, and note that the failure is cheap to
observe anyway — see the 1.6 s against 53.5 s above.

### Should `PicardThenNewton` simply be the default? No, and cost is the weakest of the three reasons

**MEASURED 2026-08-31.** Plain Newton against `PicardThenNewton`, both under
NPC, Picard capped at the same 60 iterations, wall clock on an otherwise idle
machine:

| case | plain | | `PicardThenNewton` | | slowdown | rel. `max ψ_h` |
|---|---|---|---|---|---|---|
| | out, its | s | out, picard+newton | s | | |
| §4.2 pedestal `k = 1, n = 16` | FAIL 60 | 1.40 | **ok** 60+5 | 0.67 | **0.48** | — |
| §4.2 pedestal `k = 1, n = 24` | ok 26 | 1.37 | ok 60+11 | 1.85 | 1.35 | **6.1e-06** |
| §4.2 pedestal `k = 1, n = 32` | ok 10 | 0.98 | ok 60+3 | 2.60 | 2.65 | 3.7e-16 |
| §4.2 pedestal `k = 2, n = 16` | ok 17 | 0.61 | ok 60+2 | 0.88 | 1.43 | 0.0 |
| §4.2 pedestal `k = 2, n = 24` | ok 9 | 0.78 | ok 60+1 | 1.98 | 2.52 | 5.6e-16 |
| §4.2 pedestal `k = 3, n = 16` | ok 8 | 0.52 | ok 60+2 | 1.61 | 3.10 | 5.6e-16 |
| §4.5 layer `k = 1, n = 30` | FAIL 60 | 5.02 | **ok** 60+5 | 2.61 | **0.52** | — |
| §4.5 layer `k = 1, n = 60` | ok 10 | 3.72 | ok 60+3 | 11.41 | 3.07 | 8.9e-16 |
| §4.5 layer `k = 2, n = 16` | ok 51 | 1.81 | ok 60+4 | 0.93 | 0.52 | **1.0e-02** |
| §4.5 layer `k = 2, n = 24` | ok 13 | 1.09 | ok 60+3 | 2.15 | 1.97 | **1.6e-05** |
| §4.3 barrier `k = 1, n = 16` | FAIL 60 | 1.52 | **FAIL** 0+60 | 2.22 | 1.46 | — |
| §4.3 barrier `k = 2, n = 16` | FAIL 60 | 2.28 | **FAIL** 0+60 | 3.14 | 1.38 | — |
| similarity `k = 2, n = 16` | ok 4 | 0.16 | ok 13+1 | 0.30 | 1.87 | 1.8e-16 |
| similarity `k = 3, n = 16` | ok 4 | 0.28 | ok 13+1 | 0.54 | 1.93 | 1.8e-16 |

**Three reasons, in increasing order of how decisive they are.**

**1. It costs 1.4× to 3.1× on everything that already works**, which is real
but is the reason you would trade away for robustness if the other two did not
exist. Note also that this is with Picard *capped at 60*; the cap is a tuning
parameter and *Picard, then Newton* has already measured the handoff to be
**non-monotone** in Picard effort, so 60 working on these fourteen rows is not
a property to build a default on.

**2. It does not cover everything.** The barrier fails under plain Newton and
fails under `PicardThenNewton`, and the `picard 0` in those rows is Picard
itself throwing on the first iteration. A default that is 2× slower and still
needs a fallback is not a default.

**3. IT SILENTLY CHANGES WHICH DISCRETE SOLUTION YOU GET, AND THIS IS THE ONE
THAT SETTLES IT.** Look at the last column on the rows that converge both ways:
6.1e-06, 1.6e-05, and **1.0e-02** on §4.5 at `k = 2, n = 16`. These are
*converged* solves of the same discrete system to `rel_tol = 1e-12`, differing
by up to one percent in `max ψ_h`. Run three ways rather than two, the same
case gives three answers:

| route | Newton steps | `max ψ_h` |
|---|---|---|
| NPC, plain | 51 | 3.1831e-01 |
| CondenseThenLinearise, plain | 10 | 3.4779e-01 |
| NPC, `PicardThenNewton` | 4 | 3.1514e-01 |

**A spread of 9.4%.** This is the multiple-solution finding recorded under
*The suite is green* — coarse discretisations of these sources
carry more than one solution — meeting the solver from a third direction. It
means changing the default globalisation is **not** a performance decision:
it changes the equilibrium meq reports on an under-resolved mesh, which is
precisely the mesh an adaptive run starts from.

So `Globalisation::None` stays the default and the ladder stays reactive. What
the ladder is for is a solve that would otherwise not happen at all, and paying
for it only then is what keeps the answer on the branch plain Newton finds.

**AND THE OTHER ORDERING IS NOT A CHEAPER RUNG THAN PICARD — MEASURED, BECAUSE
IT LOOKED LIKE IT WOULD BE.** `CondenseThenLinearise` converges both cases NPC
fails on, in 13–14 iterations against Picard's ~200, so switching ordering on
failure looks like the obvious cheap rung. It is not: an iteration is not the
currency, and a condensation iteration carries thousands of element-local
non-linear solves.

| case | NPC plain | `CondenseThenLinearise` | NPC `PicardThenNewton` |
|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | FAIL 60, 1.48 s | ok 14, 1.91 s | **ok 5, 0.74 s** |
| §4.5 layer `k = 1, n = 30` | FAIL 60, 5.31 s | ok 13, 5.75 s | **ok 5, 2.84 s** |
| §4.5 layer `k = 2, n = 16` | ok 51, 1.87 s | ok 10, 1.99 s | **ok 4, 1.04 s** |
| §4.3 barrier `k = 1, n = 16` | FAIL 60, 1.60 s | FAIL 60, **53.5 s** | FAIL 60, 2.27 s |
| §4.3 barrier `k = 2, n = 16` | FAIL 60, 2.36 s | FAIL 60, **96.6 s** | FAIL 60, 3.77 s |

**`PicardThenNewton` is between 2.0× and 2.6× cheaper than the ordering switch
on every case either of them cures**, and on the cases neither cures the
ordering switch costs 25× to 30× more to find that out. The ladder meq already
has is the right one and nothing needs adding to it.

The two cures also do **not** agree at these resolutions — 2.5e-13 on the
pedestal, but 1.0e-03 on §4.5 at `n = 30` and 9.4e-02 at `k = 2, n = 16` — which
is reason 3 again, and a second argument against reaching for whichever cure is
convenient.

**A corollary for anyone comparing the two orderings**: do it on a *resolved*
mesh. A disagreement at `k = 2, n = 16` is this, not a defect.
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally` runs on Example 5
at `n = 8` where both converge in 4 steps, and asserts 1e-9 — that bound is
only reachable because the discretisation there has one solution.

### KINSOL: the two orderings are exact opposites

Measured over §4.2 at `k = 1, n = 32` and `n = 40`, `k = 2, n = 24`, `k = 3,
n = 32`:

| | NPC | CondenseThenLinearise |
|---|---|---|
| plain Newton | ok, 8–10 | ok, 7–10 |
| `KIN_NONE` | **ok, 8–10** | **fails at 60, every case** |
| `KIN_LINESEARCH` | **fails at 60, every case** | ok, 24–31 |

**The wiring is ruled out, and sharply.** Under NPC `KIN_NONE` goes through the
same `ShiftedResidual`, the same `DarcyNPCOperator` and the same
`DarcyNPCSolver` and reproduces plain Newton's iteration count **exactly** — 9
against 9, 8 against 8 — so the residual, the Jacobian and the sign convention
are all correct and the fault is in the globalisation alone.
`kinsolAgreesWithNewtonWhereBothConverge` now picks the strategy that converges
for the ordering under test, which is a *sharper* adapter check than before: a
line search may take a different path to the same answer, an undamped KINSOL
may not.

**`KIN_LINESEARCH` and meq's own backtracking fail on the same cases**, both
driving an `ℓ²` merit over the whole residual — KINSOL's is `½‖fscale·F‖²` with
`fscale = 1`. **This is one finding, not two.** They fail *differently*, though,
and the difference is instructive: KINSOL applies an Armijo sufficient-decrease
condition and so gives up honestly, while the hand-written one has only a
monotone test and creeps instead. See the mechanism under *Why it fails*.

### A warm start no longer shows up in `‖r₀‖`, and that is structural

`prepare()` projects the guess onto the potential and the trace and leaves the
**flux block at zero**, deliberately — a guess for `ψ` says nothing about `q`
without differentiating it, and the guess arrives as a bare `mfem::Coefficient`,
which cannot be differentiated. While the unknown was the trace alone that cost
nothing, because `q` was a function of it. **Under NPC `q` is an unknown**, so
the guessed state is inconsistent in exactly the row that couples them and
`‖r₀‖` goes *up*: measured at `k = 3`, cold 1.771e-01 against warm 2.638e-01.

The guess still works — 2 Newton iterations against 4 cold, same L2 to every
figure, and an exact restart still finishes in 1 — because the flux row is
**linear** in `q`, so one Newton step recovers the `q` belonging to the guessed
`ψ`. `aWarmStartCutsTheWorkAndNotTheAnswer` (renamed from
`…CutsTheFirstResidual…`) asserts strictly fewer iterations and an unmoved
answer. **The fix, if the stronger property is wanted, is to seed the flux** —
`darcyFlux = −(1/r) ∇̄ψ_guess` via `GradientGridFunctionCoefficient`, for the
`setInitialGuess( GridFunction const & )` overload that has something
differentiable to work with. Not done, and it wants its own measurement.

### One more test moved from the stopping rule to the property

`MillerConvergence::diagnosticExactSolutionOnThePolygon` asserted
`newtonIterations() <= 1` on the Solov'ev source, whose `∂F/∂ψ` is zero, so the
system is affine and step one must be exact. **It was measuring the stopping
rule.** MFEM stops at `max(rel_tol·‖r₀‖, abs_tol)`, so whether an exact step is
also the *last* step depends on where `‖r₁‖/‖r₀‖` — round-off over the residual's
scale — falls relative to `rel_tol`. NPC's residual is the full system rather
than the trace, the ratio moved from just under 1e-12 to just over, and a solve
that had always had a step to spare started taking it:

```
1.466e-02  ->  1.567e-14  ->  2.708e-17
```

Twelve orders in step one, against a target of 1.466e-14 that 1.567e-14 misses by
7%. It now asserts the **drop**, `‖r₁‖/‖r₀‖ < 1e-8`, which the sweep meets
between 4.8e-13 and 3.1e-11 — the floor *degrades* with refinement, because
`‖r₁‖` sits at round-off while `‖r₀‖` shrinks with the mesh — against the
`O(1)` an inexact step on an affine system would give. Rates were unmoved
throughout: 1.995 / 2.999 / 3.999.

## Newton, and the obligation it creates

meq uses **Newton**. Both papers use Anderson-accelerated Picard. This is a
deliberate departure, and it has a cost that must be respected everywhere a
source term is written.

`refs/HDG-GradShafranov.pdf` §4 states the design meq is reversing: the authors
keep `F` as opaque problem data so the solver "relies only on the discretization
of the toroidal operator `Δ*`", and pay for it by iterating on *every* source,
even one linear in `ψ` that could have been folded into the bilinear form.

Newton makes the opposite trade. It puts `∂F/∂ψ` into the operator — a mass term
`−(∂F/∂ψ)/r` on the potential block, carried through hybridization by
`DarcyForm::GetGradient`. So:

**Every `Source` must supply `dFdPsi`, and every `Profile` must supply `Prime`.**
A profile that cannot differentiate itself cannot be used. This is why
`src/meq/Profiles.hpp` keeps the Hermite cubic — it gives `f` and `f'` from the
same data — and why `SourceTests.cpp` checks `dFdPsi` against a finite difference
of `F`. That test is not box-ticking; it is what stands between a typo in a
derivative and a Newton iteration that quietly fails to converge quadratically
while still converging.

**And a source whose profiles are functions of NORMALISED flux must be a
`meq::NormalisedSource`, not a `meq::Source`** — because `ψ_ax` is then a
functional of the solution and the Jacobian acquires non-local terms through it,
which `dFdPsi` cannot carry and the finite-difference check on `dFdPsi` cannot
see. That is a whole section of its own, below.

There is independent support for the choice, from the free-boundary literature
rather than the HDG one. CEDRES++ (`refs/CEDRES.pdf`) reports that fixed-point
iterations "usually suffer from very slow convergence or even fail to converge,
which made researchers move towards Newton-type methods", and names **vertically
unstable plasmas** as a case where Picard does not converge at all. Lackner's own
review says the same from the other end: his plain Picard, eq. (3), "will
converge to the physically trivial solution ψ ≡ 0 if admitted by the formulation
of the problem".

Three further things from CEDRES++, all of which bite before free boundary does.

**Differentiate the discrete residual, not the continuous one.** A continuous-
level Newton derivative for the plasma-current term exists (Blum 1989) and
CEDRES++ deliberately does not use it: *"there is no theoretical evidence that
this formula holds also for plasma equilibria with boundaries that contain
X-points. In particular the second term on the right-hand side seems to blow up
if ψ reaches a critical point."* They differentiate the Galerkin form instead.
meq gets this right for free, since `DarcyForm::GetGradient` differentiates the
assembled operator — but know that the shortcut fails exactly where the physics
is interesting.

**Normalised flux will make the Jacobian non-local.** Their profiles, like
`meq::Profile`, are functions of `ψ_N = (ψ − ψ_ax)/(ψ_bnd − ψ_ax)` on `[0,1]`.
`ψ_ax` and `ψ_bnd` are global functionals of the solution, so `∂F/∂ψ` acquires
terms through them, and those "lead to non-local entries in the stiffness
matrix". Fixed boundary with `ψ = 0` on a known `Γ` does not need the
normalisation and so does not have the problem — but the day the profiles are
driven by normalised flux, `MHDSource::dFdPsi` is incomplete, **and the existing
finite-difference test will not catch it**, because `f()` and `dFdPsi()` would be
missing the same terms.

**That day arrived, and it is done: `ψ_ax` is now an unknown of the non-linear
system.** `HighBetaConvergence` was the standing red test that said it was not,
and it is green. This is the section to read before touching normalised
profiles or starting free boundary, because three cheaper answers were measured
and killed on the way to the one that works.

`refs/GourdainContour.pdf` §V eq (39) specifies profiles the way equilibrium
codes actually do — `p = Σ aᵢ Ψⁱ`, `F² = F²_ax − (F²_ax − F²_bnd) Σ bᵢ Ψⁱ`, with
`Ψ` the normalised flux — and
`refs/CowleyHighBeta.pdf` / `refs/HsuCowleyHighBeta.pdf` construct high-`β_p`
equilibria of that kind asymptotically. `tests/analytic/HighBetaPoloidal.hpp` is
that source.

### What did not work, in the order it was tried

1. **With `ψ_ax` fixed, the profile is inert.** The solution reaches
   `Ψ = 0.0013`, so a peaked `p ~ Ψ^ν` contributes `Ψ^(ν−1) ~ 1e-9` of itself.
   Pressure amplitudes of 1 and 512 gave `ψ ∈ [2.570e-07, 1.259e-03]` — identical
   to every digit. A fixed `ψ_ax` is not a simplification of a normalised
   profile; unless the value happens to be right it is a **different problem**.
   `theFixedNormalisationIsADifferentProblem` still asserts that, as the control.
2. **So "is it a stiff case?" was unanswerable, and the answer looked like no
   for the wrong reason.** Every configuration converged in 1–2 Newton steps at
   reaction ratios `max|∂F/∂ψ|/λ₁` from 0.09 to **2523** — the current hole is
   unsolvable at 26. Nothing was stiff because no non-linearity was switched on.
   It also showed the ratio was being sampled over `[0, ψ_ax]`, a range the
   solve never visits: **that diagnostic needs the solution's actual range.** It
   now gets it, and the answer is **1.9 at `ν = 2` and 13–14 at `ν = 4`** — real
   but not extreme, between the pedestal's 7 and the hole's 26.
3. **Closing the loop outside the solver does not rescue it, and the reason is a
   pole.** Iterating `ψ_ax ← max ψ` with relaxation converges — to
   `ψ_ax ~ 1e-12`, a degenerate fixed point where `ψ` and `ψ_ax` shrink together,
   `Ψ` stays `O(1)`, and the pressure gradient `νA/ψ_ax` runs to 1e12 while the
   solution it drives sits at 1e-12. Mapped afterwards, the cause is plain: the
   outer map `ψ_ax ↦ max ψ` has a **pole beside its own fixed point**. At `ν = 2`,
   `A = 1` the fixed point is 0.3059 and the pole is at 0.2996, six parts in a
   thousand away. Relaxing harder is not a fix for a pole.

**And there is a second solution that is not the equilibrium, which is separate
from the normalisation and cost its own afternoon.** At a *fixed* `ψ_ax` this
equation has a small positive solution and a large one — measured at `ν = 4`,
`A = 1`, `ψ_ax = 0.42` they are 3.0e-3 and 6.2e-1. Only the large one can satisfy
`max ψ = ψ_ax`, so the constraint takes the small branch out of the solution set
— **but not out of the iteration's reach**, and Newton from the Dirichlet datum
walks straight onto it. So this path needs `setInitialGuess()` with a bump of
about the right height, and that is part of the problem statement rather than an
optimisation. It is *not* the trivial-branch trap: unlike every GS-2 §4.2–4.5
source this one does not vanish at `ψ = 0`, the smallest `|F(r,z,0)|` over the box
being 0.2.

### What works: `ψ_ax` inside the residual, as a bordered Newton

`meq::NormalisedSource` is the interface — a `Source` that also has
`setNormalisation()` — and
`GradShafranovSolver::setSource( NormalisedSource &, double )` takes `ψ_ax` as an
unknown. The system Newton closes is

```
R( λ, s ) = 0                     the hybridized trace residual, source normalised by s
G( λ, s ) = s − max ψ_h( λ, s )   = 0
```

with the Jacobian bordered,

```
[  A    c  ]     A = ∂R/∂λ            the existing hybridized Jacobian
[  bᵀ   d  ]     c = ∂R/∂s            dense
                 b = −∂(max ψ_h)/∂λ   sparse
                 d = 1 − ∂(max ψ_h)/∂s
```

**`c` and `b` are the non-local terms**, in exactly the sense CEDRES++ means when
it warns that a normalised profile "leads to non-local entries in the stiffness
matrix": `ψ_ax` is a functional of the whole solution, so a perturbation of the
trace near the magnetic axis moves the source *everywhere*.

**Why it cannot be a rank-one update inside the element blocks**, which is what a
CG code would do. In an `H¹` discretisation `ψ_ax` is one entry of the global
unknown and the Jacobian simply acquires a rank-one term. Hybridization
eliminates flux and potential **element by element**, and a term coupling every
element to the one element holding the axis is precisely what that elimination
cannot represent. The border is where it goes instead — and it costs **one
factorisation and two backsolves**, not a second matrix: solve `A y = R` and
`A z = c`, then `δs = (b·y − G)/(d − b·z)` and `δλ = −y − z δs`.

**`c` is dense and `b` is not, and the asymmetry is structural.** `s` enters
every element's source, so `∂R/∂s` has an entry on every trace dof — one central
difference in a *scalar*, two residual evaluations. `max ψ_h` is one nodal value
of one element, and under hybridization that element's recovered potential
depends only on the trace dofs of its own faces, so `b` has at most `3(k+1)`
entries and the rest are **exactly zero**. Measured, not assumed:
`theAxisSensitivityIsLocalToItsElement` perturbs trace dofs off that element and
`max ψ_h` moves by **0.0000e+00**, against 7.0e-5 for the element's own dofs.

**Both borders are differenced rather than assembled, and there is no choice
about it.** They are derivatives of the *condensed* residual. Assembling them
would need the sensitivity of the element-local eliminations — for `c` the
derivative of each local solve with respect to a parameter of its own source, for
`b` the derivative of the recovered potential with respect to the trace — and
`DarcyHybridization` exposes neither. Differencing the assembled residual is the
same principle CEDRES++ states for the local term: differentiate the **discrete**
residual, never the continuous equation.

**`ψ_ax` is the largest NODAL value, and that is a definition rather than an
approximation.** It differs from the maximum of the polynomial by `O(h^{k+1})`
and both converge to `max ψ`. The nodal one is chosen because it is what makes
the constraint differentiable in a form the border can use.

**A backtracking line search is not optional here.** The full step converges for
the mild profiles and wanders for `ν = 4` at amplitude 10, the augmented residual
reading 2.7e-1, 6.6e-2, 2.0e0, 7.7e-1, 2.5e2, 6.7e4, 7.7e6, 3.8e8 before `ψ_ax`
crosses zero and the source refuses the normalisation. That is not the Jacobian
being wrong — the same Jacobian finishes the milder cases at observed order 2 —
it is the equilibrium being a **mountain-pass solution of a superlinear problem**,
where the linearised operator is indefinite and an undamped step leaves the
basin.

**The printed residual is `‖(R, γG)‖` with `γ` frozen at `‖c‖` from the first
iterate.** `G` is a flux and `R` is a trace residual, so the two cannot simply be
concatenated; `‖c‖` is the factor that converts a perturbation of `ψ_ax` into the
units `R` is measured in. Freezing it keeps the history a comparison of like with
like — a `γ` recomputed each step would put the Jacobian's own variation into the
convergence history and manufacture orders out of it. The **control** path gets
the same `γ`, computed the same way, or the two histories would not be
comparable.

### The measurement, and the control that makes it mean something

`k = 2`, `n = 8`, converged `ψ_ax` against the dimensional estimate
`√(νA/λ₁)`:

| `ν` | `A` | `ψ_ax` | `ψ_ax − max ψ_h` | ratio to estimate | `max\|∂F/∂ψ\|/λ₁` | Newton |
|---|---|---|---|---|---|---|
| 2 | 1 | 3.058984e-01 | 1.7e-16 | 1.021 | 1.88 | 4 |
| 2 | 10 | 9.607537e-01 | −4.4e-16 | 1.014 | 1.91 | 4 |
| 2 | 100 | 3.036075e+00 | 4.4e-16 | 1.013 | 1.91 | 4 |
| 4 | 1 | 2.834510e-01 | 3.1e-12 | 0.669 | 13.14 | 8 |
| 4 | 10 | 8.643745e-01 | 8.9e-16 | 0.645 | 14.13 | 11 |
| 4 | 100 | 2.720222e+00 | 1.0e-10 | 0.642 | 14.27 | 10 |

At `ν = 2` the source is linear in `ψ` and the problem is a linear eigenvalue
problem, which is why `ψ_ax` is converged in the mesh to six figures between
`n = 8` and `n = 16` — 3.058984e-01 against 3.058988e-01. At `ν = 4` it is
genuinely non-linear and the same pair reads 0.28345 and 0.28608.

**`Normalisation::Decoupled` is the control, and it is what a test of this can
actually assert on.** It is the same solver, mesh, guess, line search and
stopping rule with exactly three quantities zeroed — `c`, `b` and `d − 1` — so
the step in `ψ_ax` reduces to `ψ_ax ← max ψ_h` and the trace step is a Newton
step that does not know `ψ_ax` is about to move. That is `ψ_ax` outside the
residual, done as favourably as possible. Measured at `ν = 2`, `A = 1`:

| | residual | iterations |
|---|---|---|
| coupled | 8.310e-02 → **4.447e-15** | 4 |
| decoupled | 8.310e-02 → 8.239e-02 | 15, not converged |

**It does not converge slowly. It does not move.** That is the test the ROADMAP
asked for when it said the work "needs a test that can *see* the missing terms",
and it is needed for the reason recorded under *A wrong Jacobian is invisible to
a convergence table*: `SourceTests`' finite-difference check on `dFdPsi`
structurally cannot see this one, because `f()` and `dFdPsi()` are both evaluated
at whatever normalisation is set and agree with each other however wrong it is.

### The trap that cost the most, and it is in MFEM's local solves

**`DarcyHybridization` captures the element-local Newton's initial guess at
`FormLinearSystem()` time and keeps it for the life of the reduced system.**
`EliminateVDofsInRHS` copies the flux and potential blocks into `darcy_u` and
`darcy_p`, and every local solve in every subsequent residual evaluation,
gradient and recovery starts there — however far the trace has since travelled.
`ComputeSolution()` does **not** use the output vector as a guess, so seeding
that is inert.

Left alone this is a performance problem on an ordinary Newton path and a
**correctness** problem on this one. Measured at `ν = 4`, amplitude 10, `n = 16`,
`k = 2`: 40,000 to 60,000 element-local iterations per outer step, most of them
hitting the cap of 100 — and a local solve that ran out of iterations returns
whatever it had reached, which is not a function of anything. Differencing
`ψ_ax` by 9e-6 then moved `max ψ_h` from 0.8961 to **2.04**, and on the next step
to **3.84**. The corner of the border read 1.6e5 where it should read about 1,
the step in `ψ_ax` collapsed to 1e-8 against a constraint residual of 3e-3, and
the iteration stalled — looking exactly like a singular border and being nothing
of the kind.

`solveWithNormalisation()` therefore re-forms the system from the recovered state
once per accepted step (`formSystem()`, factored out of `prepare()` for this).
Every local solve is then within an iteration or two of its answer, so it
converges, so it is continuous in `ψ_ax`, so the difference is a derivative. The
same run then finishes in five steps.

**Anything else that differentiates a hybridized residual by differencing it will
hit this.** It is the reason a solver-level finite difference is not simply "one
more residual evaluation".

### What is not done

`ψ_bnd` is zero, because meq solves the fixed-boundary problem with `ψ = 0` on
`Γ`. Free boundary makes it an unknown as well, which is a **second** border row
and column of the same shape — `meq::NormalisedSource` is where it goes.

`Globalisation` other than `None` is refused on this path, loudly: the KINSOL
paths drive a residual of their own and the Picard ones build no Jacobian to
border. **No ordering is refused, and that has been settled twice.** A guard once
refused `NonlinearOrdering::LineariseThenCondense`, written from
`darcyhybridization.hpp`'s summary rather than from the code; it was removed when
both orderings were measured to the same `ψ_ax`, and MFEM has since deleted that
mode outright. Under `NonlinearOrdering::NPC` the question does not arise at all
— `ψ` is an unknown, so **the border row is exactly `−e_j` and the corner is
exactly `1`**, neither of them differenced. See *The NPC port*.

`meq::NormalisedMHDSource` is the production source built on two `meq::Profile`s
in normalised flux. It is unit tested and **not yet wired into `Config` or
`SourceFactory`** — that is driver work and belongs with the rest of it.

**What working Newton looks like.** CEDRES++ Table 2, on 577k unknowns: relative
residual `2.7e0 → 9.2e-2 → 1.8e-3 → 5.3e-6 → 3.9e-12` in five iterations. That is
the shape stage 4 should produce. A run that grinds down linearly means the
Jacobian disagrees with the residual.

meq's own, measured at `k = 3`, `h = 0.1`, 832 trace dofs:

```
   it            ||r||    ||r||/||r_0||    order
    0     1.121994e+01     1.000000e+00        -
    1     4.085930e-02     3.641669e-03        -
    2     5.744224e-04     5.119658e-05    0.759
    3     1.411243e-07     1.257800e-08    1.949
    4     4.142741e-14     3.692303e-15    1.810
```

### The Solov'ev coefficients were wrong twice, and are now checked

`Soloviev.hpp`'s `nstx()` does **not** use the coefficients printed in
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (22c); `nstxAsPublished()` keeps those.
It also does not use this file's *first* correction of them. Both errors are
recorded in the fixture's header, and the shape of the second is the more
instructive.

**Error one: the published set satisfies none of Cerfon & Freidberg's twelve
constraints.** `refs/CerfonFreidberg.pdf` §IX eq. (28) gives them for an up-down
asymmetric single null. At the four points where `ψ` must vanish the printed
coefficients give `−7.5e-3`, `+3.9e-4`, `+2.9e-2`, `−9.8e-3` — the third being
11% of the axis flux.

**Error two: `α` is not `δ`.** The re-solve substituted `sin α = δ` into eq.
(11)'s `N₁ = −(1+α)²/(εκ²)` and `N₂ = (1−α)²/(εκ²)`. But those mean `α` itself,
which for `δ = 0.35` is `arcsin(0.35) = 0.3576`. `N₃ = −κ/(ε cos²α)` is immune,
since `cos²α = 1 − δ²` either way — so the error hit two of three conditions and
nothing else at all.

Settled by differentiating C&F's model surface eq. (9) directly at `τ = 0, π,
π/2`: `N₁ = −0.5907049043`, `N₂ = +0.1322804125`, `N₃ = −2.9220542041`, each
matching the `α` form exactly and the `δ` form by 1.1% and 2.4% at the two
equatorial points.

**The moral, which is the transferable part: checking a solve against the
formula it used cannot detect a misread formula.** The first correction verified
the constraints, and passed, because it verified them against its own wrong `N`.
Only an independent quantity catches that — here the curvature of the surface the
coefficients are supposed to reproduce. The corrected set matches it to 2e-11,
3e-10 and 2e-12.

**And they are no longer asserted nowhere.**
`tests/convergence/SolovievGeometryConvergence.cpp` evaluates all of C&F's
conditions on every set in the fixture, and `nstxUsesTheCorrectAlpha` checks both
readings of `α` side by side so a silent revert fails rather than passing both
ways. Any statement elsewhere in the tree that the coefficients are unchecked is
out of date.

The reason this needed a dedicated test at all: every `ψ_i` is `Δ*`-harmonic, so
**any** coefficients leave `F`, `Δ*ψ` and every convergence rate exact. Nothing
else in the suite can see a wrong one. The absolute-error ceilings move, which is
why they are recorded beside each `checkOrder()` call with the value they were
set from.

`ExtensionConvergence` still takes `Γ` to be the interior surface `ψ = −0.03`
rather than `ψ = 0`: with correct coefficients `ψ = 0` *is* the separatrix, which
passes through an X-point — a **corner** of `Γ`, where both transfer-path
families give out and the Cockburn–Solano analysis does not reach.

**A tooling warning that has now cost time twice.** `pdftotext` silently drops
this paper's minus signs, and an `ε` in another. Read the rendered page.

### A wrong Jacobian is invisible to a convergence table

The single most useful thing measured in stage 4, and the justification for
asserting on Newton's *order* rather than merely on the rates.

Perturbing `∂F/∂ψ` by **+5%** and re-running the whole Example 5 study leaves
**every error and every convergence rate unchanged to six significant figures**.
The discretisation is untouched, because Newton converges to the same discrete
solution whatever Jacobian carried it there. What changes is only the path:
observed order drops to exactly `1.000`, and it takes 10 iterations instead of 4.
At +50% Newton diverges to NaN; drop the Jacobian mass term entirely and McCarthy
takes 87 iterations instead of 1, at a constant contraction of 0.7711 per step,
with the errors again unchanged.

So a rate table cannot see a Jacobian error at all. Three things can, and all
three are in the suite: the finite-difference check on the assembled Jacobian
(4e-11 relative, at the `O(step²)` floor), the assertion on observed Newton
order, and the McCarthy rung, whose affine source must finish in exactly one
step.

### Why meq's Newton struggles where other Newton solvers do not

**READ THIS FIRST: THREE OF THE FOUR "STIFF" SOURCES WERE MERELY
UNDER-RESOLVED.** The rest of this section is kept because its falsified
hypotheses are worth not repeating, but its premise is largely wrong. Raw
Newton, undamped, cold start, iterations to converge over `n = 16, 24, 32, 48`:

| | `k = 1` | `k = 2` | `k = 3` |
|---|---|---|---|
| §4.2 pedestal | **42**, **22**, 9, 9 | 10, 9, 8, 6 | 11, 9, 7, **7** |
| §4.3 barrier | fail, fail, **7**, 8 | 10, 7, 7, 8 | —, 8, 8, 7 |
| §4.5 layer | fail, **17**, 11, 11 | —, 12, 14, 11 | —, 10, 21, 12 |
| §4.4 current hole | **abort, abort, fail** | abort, fail, fail | —, **fail, fail, fail** |

**§4.2, §4.3 and §4.5 are ordinary problems on a resolved mesh** — 7 to 17
iterations — and both refinement paths cure them independently: `h`-refinement
takes §4.2 at `k = 1` from 42 to 9, and `p`-refinement takes `n = 16` from 42 to
10 without touching the mesh. `k = 3, n = 48` on the pedestal takes **7**. The
"meq is doing the easier problem and finding it harder" red flag below was
raised from a benchmark run at a single under-resolved point, `k = 1, h = 0.05`.

**§4.4 is the exception and is genuinely unsolved.** It fails at *every* order
and every mesh tried, including `k = 3, n = 48`, spending 1.8M element-local
iterations to do it. Refinement does nothing, and neither does Picard.

**And the reason is measured, not guessed: the Jacobian's reaction term has swept
past ~26 eigenvalues of the operator it is added to.** `∂F/∂ψ` for eq. (26)
ranges over `[−579.5, +565.7]`, against a first Dirichlet eigenvalue of
`λ₁ = π²(1/w² + 1/h²) = 22.3` on the benchmark box. The linearised operator
`−∇̄·((1/r)∇̄·) − (∂F/∂ψ)/r` is therefore strongly indefinite, and the continuous
problem is **multi-valued**. The driver is the `c₃(1 − e₂)cos(c₄ψ)` term at
`c₃ = −18`, `c₄ = 10π` — the very feature that empties the core of current —
whose derivative carries `c₃c₄ ≈ 566`.

That is why refinement is powerless: there is no discretisation error to remove.
Ratio against `λ₁`, by amplitude:

| `c₃` | 0 | −2 | −4 | −9 | **−18** |
|---|---|---|---|---|---|
| `max\|∂F/∂ψ\| / λ₁` | 7 | 7 | 8 | 13 | **26** |

Note the pedestal itself already sits at 7 and converges, so exceeding `λ₁` is
not by itself fatal; 26 is. **This supersedes an earlier reading of §4.4 as a
trivial-branch problem.** `F(r, 0) = 0` is true and the trivial branch is real —
see *Traps* — but it is not what defeats the iteration here, since the runs above
carry non-homogeneous ramp data that keeps `ψ` away from zero and fail anyway.

Do not file §4.4 under stiffness.

**§4.4 does have a solution, and continuation reaches it.** Adaptive
continuation in the added term's amplitude — halve the step on failure, grow it
by 1.3 on success — walked `c₃` from 0 to −18 at `k = 2, n = 32` in **9 solves
with 2 retreats**, the last taking 8 Newton iterations. Uniform steps do not do
it (ten of them stall at `c₃ = −10.8`), so the step control is the point, and
there is **no limit point on the branch** — the step never had to fall below
1e-3, it simply had to shrink over the last sixth.

**THIS MUST NOT GO INTO THE DRIVER, AND THE REASON IS NOT TASTE.** What was
ramped is `c₃`, a constructor argument of the `CurrentHole` *test fixture*. The
`Source` interface exposes `f( r, z, ψ )` and `dFdPsi( r, z, ψ )` and nothing
else — **there is no amplitude parameter in it and no way to recover one from a
black-box `F`**. The continuation is unavailable to a driver on principle, not
merely inadvisable. The black-box analogue, `F_λ = λF`, *is* expressible but is
**untested and a different homotopy**: at `λ = 0` it degenerates to the harmonic
problem, where the path above starts from the converged *pedestal*.

**And it cannot be predicted that continuation is needed.** `max|∂F/∂ψ|/λ₁` is
computable black-box, `dFdPsi` being mandatory — but the pedestal converges at 7
and the hole fails at 26 (two points, no threshold), and the ratio needs the
range of `ψ`, which is not known before solving. A detector calibrated on that
would be fitting noise. **The driver gets a reactive ladder, never a predictive
one**; see also *There is no cheap discriminator*, where a second candidate was
measured and turned out to be anti-correlated. If continuation is ever wanted
generally the prerequisite is a `Source` that can parameterise itself, which is
an interface change to argue on its own merits.

So the remedy for a production run is **resolution** — the adaptive loop rather
than a globalisation. What globalisation buys is the **coarse start** an
adaptive loop necessarily begins from, which is exactly where raw Newton fails.

What survives below: hybridization really does put a non-linear solve inside
every element elimination under `CondenseThenLinearise`, and that is why an
under-resolved §4.2 grinds where a CG code would not. What does **not** survive
is the conclusion that this makes the problems unreachable — it makes an
under-resolved discretisation expensive, and NPC removes the local solves
without removing the difficulty.

**Two hypotheses below are falsified and one is now known to be incomplete.**

Three Grad–Shafranov codes solve this equation by Newton and report it robust:

| | discretisation | nonlinear systems per Newton step |
|---|---|---|
| Serino, Tang, Tang, Kolev & Lipnikov (`refs/MFEM-GS-Newton.pdf`) | CG, `ψ ∈ H¹`, MFEM | **1, global** |
| CEDRES++ (`refs/CEDRES.pdf`) | CG finite elements | **1, global** |
| FreeGSNKE | 4th-order finite differences, Newton–Krylov | **1, global** |
| **meq** | **hybridized HDG** | **1 global + one per element** |

**In a CG or finite-difference discretisation `ψ` is one global unknown vector
and `F(ψ)` enters only the global residual**, so Newton linearises once — Serino
et al. eq (3.1) and eq (3.7) are the whole nonlinear structure. **Hybridization
changes that.** Static condensation expresses `(q, ψ)` element by element in
terms of `ψ̂`, and when `F` depends on `ψ` that elimination is *itself a
nonlinear solve, once per element per residual evaluation* — none of them
globalised, and any one failing poisons the whole residual. That is what
`el: N not convered in 100 iters` is, and it is why globalising the outer
iteration does nothing; see *On SUNDIALS*.

**Two comfortable explanations are wrong, and the papers say so directly.**
Newton is not the problem: Serino et al. built their Newton solver precisely
because "conventional Picard-based solvers fail to converge" on the Taylor
state, and report the residual reaching 1e-6 "in a small handful of iterations".
Free boundary is not the problem either — they note the fixed problem is
"significantly easier", so **meq is doing the easier problem and finding it
harder**. It also explains the GS papers' choice: keeping `F` as opaque problem
data leaves **every local solve linear**, which in a hybridized method is
coherence rather than fastidiousness.

**AND THE STRUCTURAL STORY WAS ITSELF FALSIFIED.** This section predicted that
applying Newton to the full `(q, ψ, ψ̂)` system — `refs/HDG-NPC-2.pdf` §2.6,
which is what `NonlinearOrdering::NPC` is and what meq now defaults to — would
fix the stiff sources. It does not. NPC has **no element-local non-linear
iteration at all**, confirmed at `GetNumLocalNLIterations() == 0`, and it loses
on precisely the under-resolved cases this section is about. What it buys is
uniform local work and 3.7×–4.3× of wall clock, not robustness. **The
element-local non-linear solves were never the cause.** See *The NPC port*.

**What follows, in order of how targeted it is.**

1. **The local solver was never chosen.** `SetLocalNLSolver` offers `Newton`,
   `LBFGS` and `LBB`, and meq hardcoded `Newton` — undamped, on exactly the
   problems that were failing. Measured on §4.2 at `k = 1, h = 0.05`: Newton 42
   outer iterations, LBFGS 36, **LBB 25**. `setLocalSolver()` exposes it. **It
   is inert under NPC**, which has no local non-linear solve to configure.
2. **The control that pointed at the local solves, and why it did not prove what
   it looked like it proved.** The variable was isolated by **Picard on meq's
   own linear path** — identical mesh, spaces, hybridization and `τ`, but `F`
   evaluated at the previous iterate and handed to `setSource( Coefficient & )`,
   which makes the potential block linear and every local elimination a linear
   solve. On §4.2 at `k = 1, h = 0.05`, the case that fails: Newton **fails**,
   undamped Picard **stalls at 3.5e-1**, and Picard at `ω = 0.5` reaches
   **2.8e-8 and falling**. So the mesh is fine, the discretisation is fine and
   the problem *is* solvable there.

   **But that control differed in GLOBALISATION as well as in local linearity**,
   and the later globalisation cross says which half was doing the work — the
   damping, not the local linearity. It is why the `ω = 1` row matters: undamped
   Picard stalls, which is the weakness CEDRES++ and Serino et al. both report
   and the reason the GS papers use **Anderson**-accelerated Picard rather than
   plain. **A diagnosis, not a recommendation**: relaxed Picard took 200
   iterations to reach 2.8e-8 where Newton takes 42 on the mesh Newton manages.

**A LARGE BLOCK THAT STOOD HERE IS DELETED.** It was the investigation of
`NLOrdering::LineariseThenCondense`, a mode MFEM has since retired as *"a
condensation in disguise"*. It is not re-runnable, it describes nothing meq
does, and it is recoverable from git (`cbb210b`, `05c864d`). Four things in it
were worth keeping and are kept here.

**The prediction in this section was falsified.** It argued that condensing
before linearising was the *cause* of meq's stiff-source trouble and that
reversing the order would fix it. The reversed ordering landed, was confirmed
genuinely in effect by `GetNumLocalNLIterations() == 0`, and was **strictly
worse on every stiff case**. NPC — which is the canonical version of the same
idea and is meq's default today — reproduces that verdict; see *The NPC port*.
**Do not expect an ordering to fix a stiff source.** It is the standing warning
against reasoning from structure to robustness in this file.

**A measurement technique worth reusing: cross-linearisation residuals.**
Evaluate the reduced residual twice at one trace, relinearising *at a different
trace* in between. Under condense-first the difference is exactly 0.000e+00 at
every local stiffness — it must be, since the local problems are solved to
convergence and nothing about how they were reached survives. Any nonzero
answer says the operator carries hidden state and no Newton can converge on it.
That is what condemned the deleted mode, at 149% of the residual's own value.

**AND IT FOUND A LIVE DEFECT IN THE ORDERING MEQ STILL SHIPS AS THE BACKUP.**
Under condense-first the residual is a perfect function of the trace, but the
assembled gradient disagrees with a central difference of it **by a factor of
100** once the source's local width reaches `σ² ≤ 0.01`. That is the
element-local solves hitting their 100-iteration cap: a residual computed from
an unconverged local solve is not smooth, so the difference is meaningless and
the Jacobian may or may not be. **It is the measured explanation for the
pedestal's wandering iterations at `k = 1, n = 16` under
`CondenseThenLinearise`**, and it is the same defect that corrupts a differenced
border — see the seed note under *Traps*. Refinement cures it, because a
resolved local problem converges inside its cap.

**Globalising the outer iteration makes the local solves worse, not better**,
which is the other half of *On SUNDIALS*. On §4.2 at `k = 1, n = 16` under
condense-first: plain Newton converges in 31 with 133,168 element-local
iterations; `KIN_LINESEARCH` fails at 45 having spent **1,381,527**; `KIN_NONE`
fails at 19. A line search chooses how far to move the trace and cannot make
the local problems at that trace well posed.

3. **Picard, keeping the local problems linear — implemented, as a bridge.**
   `Globalisation::AndersonPicard` is `KINSolver(KIN_FP)` over the fixed point
   `ψ^{k+1} = G(ψ^k)`, where `G` freezes `F` at the previous iterate, puts it on
   the right-hand side and does one linear solve. The papers' own method, and it
   works where Newton grinds. On §4.2 at `k = 1, h = 0.05`:

   | method | outcome |
   |---|---|
   | Newton | 42 iterations |
   | Picard, undamped | stalls |
   | Picard, `ω = 0.5` | 248 iterations |
   | **Anderson, depth 1, undamped** | **162 iterations** |
   | Anderson, depth 2 and above | fails |

   **Two surprises, both defaults now set from measurement rather than from the
   papers.** Plain Picard *needs* damping and Anderson does not — so
   `setPicardDamping` defaults to 1.0, which is right for one path and wrong for
   the other. And **HDG-GS-1's `m = 2` fails here** where `m = 1` converges, so
   `setAndersonDepth` defaults to 1; whether that is this fixed point's
   conditioning or KINSOL's implementation is not established, and raising it
   expecting the papers' behaviour will not work.

   It is a **robustness route, not a faster one**: 162 iterations against
   Newton's 42, each one a full linear solve.
   `andersonPicardReachesTheSameSolutionAsNewton` pins that both reach the same
   discrete solution — `ψ` agreeing to seven figures — which is what makes it an
   alternative rather than a different problem.
4. **Continuation in the source amplitude**, which also addresses the trivial
   branch.

### Picard, then Newton — the route for a coarse mesh

**This is what reaches the hard cases without refining them, and it is
`Globalisation::PicardThenNewton`.**
Anderson-accelerated Picard walks the iterate into Newton's basin; plain Newton
takes it from there and supplies the quadratic endgame Picard structurally cannot.
Measured, with `setInitialGuess()` seeding the ramp on every row:

| case | Newton alone | Picard | **Newton from Picard** |
|---|---|---|---|
| §4.2 pedestal `k = 1, n = 16` | ok, 31 | ok, 197 | **ok, 4** |
| §4.2 pedestal `k = 1, n = 24` | ok, 7 | ok, 290 | **ok, 4** — order 2.02 |
| §4.3 barrier `k = 1, n = 16` | **fails at 60** | ok, 122 | **ok, 4** — order 2.01 |
| §4.3 barrier `k = 2, n = 16` | ok, 10 | ok, 120 | **ok, 3** — order 1.96 |
| §4.5 layer `k = 1, n = 16` | **fails at 60** | 400, not converged | **ok, 28** |
| §4.5 layer `k = 2, n = 16` | **fails at 60** | 400, not converged | **ok, 5** |
| §4.4 hole `k = 1, n = 16` | aborts | not converged | **aborts** |

Three of the four unreachable cases become reachable. §4.3 at `k = 1, h = 0.05`
is the one to quote: nothing else in the solver touches it **at that
resolution** — plain Newton fails at 60, a line search fails at 24 having spent
5.4M element-local iterations, and linearise-first aborts — and the handoff
finishes in four Newton steps,
`8.3e-01 → 1.5e-03 → 1.0e-04 → 7.4e-07 → 3.7e-11`, agreeing with Picard's own
answer to **4.7e-10**. That agreement is what makes it a handoff rather than a
change of problem, and `picardThenNewtonRecoversQuadraticOrder` asserts it.

**Be precise about what this is worth, because refinement reaches the same three
cases.** Raw Newton solves §4.2, §4.3 and §4.5 perfectly well once resolved — see
the table under *Why meq's Newton struggles* — so the handoff is **not** the only
route to them, and it is not the route to prefer when refining is available. What
it is for is the **coarse mesh**: an adaptive run must solve on its initial mesh
before it has an estimator to refine with, and that first solve is exactly the
under-resolved regime where raw Newton fails. That is a real job, and it is a
narrower one than "the route for stiff sources".

**It does not rescue §4.4.** The current hole fails under Picard, under the
handoff, and at every order and mesh up to `k = 3, n = 48`. Its problem is the
trivial branch, not the iteration.

**Picard's job is not to solve the problem.** §4.5 converges at both orders from
a Picard state that never met its own tolerance, so stage 1 stopping short is an
expected outcome, not an error, and `solveByPicardThenNewton()` swallows that
throw deliberately. This is a globalisation, not a two-solver pipeline.

**Do not replace the tolerance with an iteration budget.** The handoff is *not
monotone* in Picard effort — on §4.5 at `k = 1`, budgets of 400 and 3 converge
while 40 and 10 fail, and on §4.3 at `k = 1` a budget of 3 diverges to `1e4`. A
budget tuned on one mesh will betray you on the next. Picard's own tolerance is
the trigger that worked wherever it was reached.

It is not cheap: stage 1 spent 122 to 290 iterations, each a full linear solve.
Reach for it when `Globalisation::None` fails, not before.

**A caution on reading the printed order.** Taking the best observed order over
any triple of a short or non-monotone history manufactures values of 3.2, 3.6 and
9.36 out of runs that are not converging at all. Only a monotone tail supports an
order claim; the 2.02, 2.01 and 1.96 above are those.

**And it exposed a latent defect, now fixed.** `setGlobalisation()` reset
`prepared` but not `built`, while `buildForms()` branches on
`usesNonlinearForms()` — which reads `globalisationChoice` — to decide whether the
potential block goes on the linear or the nonlinear form. Switching a live solver
between a Picard path and a Newton one therefore reused the other path's blocks.
`setNonlinearOrdering()` and `setLocalSolver()` beside it always did reset
`built`. Latent until now because every caller built a fresh solver per
globalisation; `PicardThenNewton` switches twice inside one `solve()`.

### On SUNDIALS

`mfem::KINSolver` **derives from `mfem::NewtonSolver`** (`linalg/sundials.hpp`)
and is reachable through `setGlobalisation()`; `SetJFNK` and `EnableAndersonAcc`
come with it. SUNDIALS 7.5.0 is built in at `../sundials/cuda-install`.

**IT IS NOT A DROP-IN, AND THE DIFFERENCE IS SILENT.** `NewtonSolver::Mult(b,x)`
forms `r = oper(x) − b`, while **`KINSolver::Mult` declares its first argument
without a name** and solves `oper(x) = 0`. meq's trace right-hand side is not
zero, so handing it straight to KINSOL converges — to the solution of a
different problem. `ShiftedResidual` in `GradShafranov.cpp` is the adapter, and
it reproduces `NewtonSolver`'s residual exactly, which is what makes any
comparison between the two paths mean anything.
`kinsolAgreesWithNewtonWhereBothConverge` would catch it being dropped or its
sign flipped. For the NPC-era behaviour see *KINSOL: the two orderings are exact
opposites*.

**Globalising the outer iteration does not globalise the inner ones**, which is
the structural point and is why a line search rescues nothing under
`CondenseThenLinearise`. The failure there is
`el: N not convered in 100 iters` — MFEM's **element-local** non-linear solve,
one per element per residual evaluation. A line search chooses how far to move
the *trace*; it cannot make the local problems at that trace well posed, and
KINSOL never sees them. What would: damping the **local** solves, which is
MFEM's to offer; continuation, so each solve starts from the previous answer; or
**Picard on the outer loop**, which leaves every local problem linear. That last
is what both papers do. **NPC removes the whole question** by having no
element-local non-linear solve at all.

**`MFEM_USE_LAPACK` APPEARS TO FIX THE PEDESTAL AND DOES NOT — IT TIPS A
MARGINAL ITERATION BY CHANGING THE ROUNDING.** Against a LAPACK build §4.2 at
`k = 1, h = 0.05` converges in 42; against `install-nolapack`, identical but for
that flag, it fails at 60. **But at `h = 0.0333` both builds converge in 23, the
same number.** Only the marginal mesh moves. There is no mechanism for it to be
anything else: meq sets `LPrecType::LU`, both implementations partial-pivot on
the same rule, and `Factor( m, TOL = 0.0 )` is `dgetrf`'s singularity condition
exactly. What is left is blocked BLAS-3 against unblocked scalar loops —
identical arithmetic, different summation order, `O(1e-16)`. **And the BLAS here
is threaded MKL**, whose reduction order depends on the thread count, so whether
that local Newton converges is machine- and environment-dependent. Do not treat
42 iterations as reproducible. An intermediate claim in this file that LAPACK
had fixed a third of the globalisation problem was wrong.

**`MFEM_USE_EXCEPTIONS` is enabled** and does what `DRIVER-PLAN.md` §5 needs:
`MFEM_ERROR_THROW` is the default error action and `mfem::ErrorException`
derives from `std::exception`, so **no meq-side change was required** — an
existing `catch ( std::exception const & )` already catches it. §4.4 at
`k = 1, n = 16` used to take the process down with SIGABRT and now reports a
failure with a usable iteration count, so a driver can return exit code 2 rather
than dying.

**A caveat worth keeping in view**: the throw unwinds out of the middle of MFEM,
and the objects are left as the throw found them. meq's paths construct a fresh
solver per solve and so do not care, but **do not assume a
`GradShafranovSolver` is reusable after a caught `ErrorException`.**

### The suite is green, and the last red one was a mesh chosen for a dead solver

**GREEN, 2026-08-31 (later still).** The NPC port dropped the whole suite from
872 s to 499 s — the element-local non-linear solves that cost `PedestalConvergence`
its time are gone, 213 s against 572 s, and `HighBetaConvergence` 3.3 s against
17 s. It left one red assertion, and that is now fixed rather than tolerated.

**What was red**: `internalLayerSelfConverges`, GS-2 §4.5 at `k = 1`, on the
**coarsest** mesh of its self-convergence sweep, `n = 24`.

**WHAT IT ACTUALLY WAS: A MESH SEQUENCE CALIBRATED AGAINST AN ORDERING MFEM HAS
DELETED.** `meshesFor()`'s table recorded `h = 0.0333: 26 it` for that point,
measured under the old default. The discretisation reason was there all along and
nobody had connected it: §4.5's ridge is **about 0.025 wide in space**, and
`n = 24` gives `h = 0.0333` — **the feature was thinner than a cell**. The
self-difference at that point was measuring the approach to the asymptotic
regime, not the rate, which is the same objection the file already records
against starting `k = 2, 3` coarser than 16.

Moving the `k = 1` sequence to `{ 48, 96, 192 }` fixes both at once, and the
rates improve rather than being relaxed to fit:

| §4.5, `k = 1` | `ψ` | `q` |
|---|---|---|
| `{ 24, 48, 96 }` | 1.786 | 2.480 |
| **`{ 48, 96, 192 }`** | **2.160** | **2.184** |

Design order for `k = 1` is 2, and the corner control caps `q` at about 2.2, so
the finer sequence sits at design order in `ψ` and at the cap in `q` — which is
what the study is for. Plain NPC converges at all three points, 12/13/12.
`PedestalConvergence` went 216 s → 262 s, and the whole suite only 509 s →
**522 s**: the 60 wasted iterations at `n = 24` paid for most of `n = 192`.

**And the failure had two distinct characters, which is worth keeping**:

| `n` | `h`/0.025 | plain NPC |
|---|---|---|
| 24 | 1.33 | fails at 60, **wandering** to `ψ ∈ [−1.51, +1.44]` |
| 32 | 1.00 | fails at 60, in a **PERIOD-4 LIMIT CYCLE** at ~4e-3 |
| 48 | 0.67 | ok, 12 |

The `n = 32` row is not divergence: the residual descends 7.99e-02 → 2.69e-03 in
ten steps and then *orbits*, period 4 to a mean relative mismatch of 0.041
against 0.32–0.44 at every other period, eleven orders short of its target.
Newton with a stable periodic orbit, which is a different animal from the
wandering at `n = 24` and from anything else recorded in this file.

**NONE OF IT IS THE DISCRETISATION.** Every route that converges reaches the same
discrete solution — `CondenseThenLinearise` and NPC-with-`PicardThenNewton` agree
to between 2.9e-11 and 2.3e-09 in L2 at every mesh from 24 to 192, and their
self-differences are **bit-identical**, so the rate does not depend on the route.
Both failing points are reachable: `PicardThenNewton` takes 11 Newton steps at
`n = 24` and **4** at `n = 32`.

**AND THE BASIN FOLLOWS THE FEATURE, NOT THE DOF COUNT.** Marking on `|∂F/∂ψ|`
at the datum ramp — an a-priori marker using nothing the solve is not given — and
refining that band once from `n = 24`:

| | trace dofs | plain NPC |
|---|---|---|
| uniform `n = 32` | 6,272 | **fails at 60** |
| **band-refined from `n = 24`** | **7,588** | **ok, 22** |
| uniform `n = 48` | 14,016 | ok, 12 |

Inside the basin on **54%** of uniform `n = 48`'s dofs, and with 21% more than
the uniform mesh that fails. The band must be generous, though: marking at half
the peak takes only 102 of 1152 elements and still fails, because the datum is
only a proxy for where the ridge sits. This does not go into the
self-convergence study — a rate needs a uniform `h` — but it is the right move
for anyone who needs a coarse solve on a localised-feature source.

**THE TRANSFERABLE LESSON, AND IT IS THE ONE THIS FILE KEEPS RELEARNING**: a
"coarsest usable mesh" is a property of the **solver** as much as of the
benchmark. That table was calibrated against a solver that no longer exists, and
the red assertion was the calibration going stale rather than anything about
meq. Re-measure it after any change to the ordering or the globalisation.

**A LONG PRE-PORT ACCOUNT STOOD HERE AND IS DELETED.** It was the 872 s run
against the deleted ordering, its per-file timing tables and the history of a
parity gap MFEM has since closed by removing the mode. Recoverable from git
(`cbb210b`, `05c864d`). Three things in it are still live.

**`pedestalConvergenceIsAResolutionThreshold` asserts the CURES, not the knife
edge.** It used to assert that Newton fails at `k = 1, h = 0.05`, which it was
never entitled to do in either direction: 42 iterations against a LAPACK build
and a failure at 60 against `install-nolapack` is one flag's difference in
threaded-MKL reduction order, so an assertion on it is an assertion about this
machine today. It now asserts the two independent refinement cures, each with a
factor of four in it:

| | | |
|---|---|---|
| `k = 1, n = 16` | 42 iterations | the knife edge, recorded and *not* asserted on |
| `k = 1, n = 32` | 9 iterations | `h`-refinement cures it |
| `k = 2, n = 16` | 11 iterations | `p`-refinement cures it, mesh untouched |

That two independent paths both cure it is what identifies the cause as
under-resolution. §4.4 is the control: it fails at every order and mesh tried,
because its trouble is that the problem is multi-valued. The threshold itself
is still asserted, as `iterations(n = 16) ≥ 2 × iterations(n = 32)`, which holds
whichever side the rounding falls on.

**THESE COARSE DISCRETISATIONS CARRY MORE THAN ONE SOLUTION, and it is measured
from three directions now.** `andersonPicardReachesTheSameSolutionAsNewton`
used to gate one mesh at 1e-6; over a sweep it reads 9.1e-05, 1.2e-13, 3.3e-06
and 4.9e-13 at `n = 16, 24, 32, 48` — round-off on some meshes and 1e-5 on
others, **with no trend in the mesh**. It is not a stopping tolerance:
tightening rtol from 1e-8 to 1e-12 leaves the `n = 16` figures *bit identical*.
Both iterations are fully converged and their fixed points differ. The test now
sweeps three meshes and asserts the two things actually entailed — the **best**
agreement is at round-off, which is what says the Picard path solves meq's
problem rather than a neighbouring one, and the **worst** is bounded well below
a different problem. A per-mesh gate at 1e-6 is not reinstatable. See
*Should `PicardThenNewton` simply be the default?* for the same phenomenon at
9.4% across three solve routes.

**`NewtonConvergence`'s finite-difference Jacobian check reads 4e-11**, the
`O(step²)` floor, and meq keeps its hoisted `GetGradient()` regardless of
whether the library requires it: holding the linearisation fixed across a
difference is the right thing to write, and upstream's own finding is that
re-taking the gradient *after* a difference silently buys extra corrections and
can hide exactly this class of defect.

## Toroidal flow: FL-0 to FL-8 are done and green

**`meq::RotatingSource` solves the generalised Grad-Shafranov equation of
`refs/RotatingGK.pdf` (136)**, closed by its (96) and (97), for two species in
the local gauge `φ₀(r_ref) = 0`. `FLOW-PLAN.md` is the design and the staged
plan; this section is only what a reader of the code needs. **Three or more
species is `Closure::RootFind`** — a safeguarded scalar Newton on (97) with
`φ₀`'s two `ψ`-derivatives by implicit differentiation — and
`meq::NormalisedRotatingSource` puts the profiles in normalised flux, where
`ψ_ax` is an unknown and the existing bordered Newton closes it.

**FL-8 IS DONE: ROTATION IS REACHABLE FROM A TOML FILE.** `[source] Type =
"rotating"` takes an array of `[[source.species]]` tables — mass, charge and a
temperature profile each — plus `Omega`, one density profile and `GGPrime`.
That array is **meq's first array-of-tables in the schema**, read by
`Table::getTableArrayOr()`, which names its elements `source.species[i]` so a
fault in the third one says so. `examples/rotating-rectangle.toml` is the
worked example and `examples/rotating-normalised.toml` the same thing in
normalised flux, where `ψ_ax` is an unknown and the bordered Newton closes it —
`makeSource` **throws** on a normalised config and `makeNormalisedSource` is
the other door, because handing a normalised source to the plain path would
converge to a different equilibrium rather than fail.

**AND ROTATION REACHES `ψ` ITSELF THERE, WHICH IS A MEASUREMENT RATHER THAN A
CLAIM.** Against the same configuration with `Omega = 0` and nothing else
changed, `‖ψ(ω) − ψ(0)‖/‖ψ(0)‖ = 1.2683e-01` in L2 over 4608 dofs, and
`max ψ` moves 9.755876e-02 → 1.079860e-01, up 10.7%.
`theDriverSolvesARotatingEquilibrium` re-measures both every run and gates the
shift at `> 5e-2`, so an edit that shrinks the pressure term fails loudly
instead of quietly re-vacuating the example. **It is not free**: it needs
`GGPrime` down from 2.0 to 0.8 so the pressure term is a real part of `F`
rather than a correction on `g g′` — 15 / 29 / **47** / 66 / 81 % of `F` from
`RMin` to `RMax` at `Omega = 4.0e5`, which is `M = 1.33` and squarely the sonic
regime RoPP is about.

**TWO PROFILE TABLES SHIP, AND THE REASON IS A TRAP.** Re-expressing `f(ψ)`
against `Ψ = ψ/ψ_ax` divides the abscissa by the axis flux and **multiplies the
derivative column by it**. `examples/rotating-density.dat` and
`rotating-density-normalised.dat` are the same profile written both ways; hand
one to the wrong configuration and it parses, solves and converges to a plasma
whose density gradient is out by a factor of ten. Both headers say so, and they
are also the only place in `examples/` that documents the `SplineProfile` file
format at all.

**(136) collapses to `F = μ₀ r² ∂p/∂ψ|_r + g g′`** with `p = Σ_s n_s T_s`,
because the `∂φ₀/∂ψ` terms cancel identically against quasineutrality. So the
residual needs `φ₀` and never its derivative; only the Jacobian does. **Two
species need no root find at all** — (97) is linear in `φ₀` after logs, giving
`C = ω²(Z₁m₂ − Z₂m₁)/(Z₁T₂ − Z₂T₁)` as the exponent both species share, exact
with the electron mass kept. **An earlier version of this line ended "three or
more is FL-6 and the constructor refuses", which FL-6 made untrue** — the root
find handles them, `Closure::Automatic` picks between the two, and the ceiling
is `meq::maxSpecies = 8`, enforced at parse rather than at construction.

**IT COST `meq::Profile` A THIRD DERIVATIVE LEVEL, AND THAT IS THE ONE
STRUCTURAL CHANGE.** `MHDSource` stores the *products* `p′` and `g g′`, so `F`
is one evaluation and `∂F/∂ψ` is one `prime()` — two levels, which is all
`Profile` had. A rotating `F` is *already* `∂p/∂ψ` of something built from flux
functions, so the Jacobian spends a second derivative of every input. No
reparametrisation avoids it. `doublePrime()` is a **pure** virtual, so every
subclass had to answer rather than one silently returning zero. **The caveat is
real**: a Hermite cubic is `C¹`, so its second derivative jumps at every
interior knot — `the_second_derivative_jumps_at_an_interior_knot` asserts that
rather than pretending otherwise, and nothing has measured what it costs Newton.

### What is measured

| | |
|---|---|
| `φ₀` against a **brentq root of (97)**, independent Python | 1.9e-14 |
| `Σ_s Z_s n_s = 0` at every radius, over a Mach sweep | 1e-13 |
| `p` against Li & Zhu (8) | 1e-14 |
| `dFdPsi` against a central difference, five Mach scales, two steps | the `O(h²)` floor |
| `ω = 0` against `MHDSource`, pointwise | 1e-13 |
| **`ω = 0` through the solver** | reproduces `SolovievConvergence`'s errors **to every printed digit** |
| rotating Solov'ev, `k = 1,2,3` | 1.995 / 2.995 / 3.995 in `ψ`, 1.980 / 2.985 / 3.984 in `q`, Newton = 1 |
| source against fixture, `M² = 0, 1, 4` | asserted at 1e-12 — two independent implementations of the same physics |
| assembled Jacobian vs a difference of the assembled residual | 3.1e-11 |
| Newton order on the manufactured nonlinear case | 1.980, and `k+1` at 2.007 / 2.999 / 4.002 |
| root find vs closed form at two species | `φ₀` 1e-12, `∂φ₀/∂ψ` 1e-11, `F` 1e-11, `∂F/∂ψ` 1e-9 |
| three species (D, C⁶⁺, e), `Σ_s Z_s n_s = 0` | 1e-12 at every radius |
| bordered Newton, `ψ_ax − max ψ_h` | **0.000e+00** on three meshes; tail order 2.000 |

**THREE PAPER ERRORS WERE FOUND ON THE WAY AND ALL THREE ARE THE KIND THAT
CONVERGE BEAUTIFULLY.** Li & Zhu's (12)–(16) write `M₀²` where their prose
defines `M₀` as the group without a square root — the group is an energy ratio,
so it is a Mach number *squared*, and `RotatingSoloviev.hpp` names its member
`machSquared` and cites the exponent rather than their symbol. Their (9) carries
**two reversed signs**, on the `dΩ/dψ` and `dT/dψ` corrections, found
independently three times — by transcription, by an unrelated numerical check,
and by meq's own derivation agreeing with the corrected form. **Neither of their
own benchmarks can see it**, because both have `dC/dψ = 0`. And their (6) omits
the `μ₀` their (9) carries.

**THE GAP THAT LEAVES IS WORTH KNOWING**: no published rotating benchmark
exercises the `C′(ψ)` term — Li & Zhu's Solov'ev case has `T` and `Ω` constant,
and Maschke & Perrin's (4.7) *forces* `C` constant — and that is precisely the
term Li & Zhu got wrong. Only `RotatingSourceTests`' `dFdPsi` sweep, over
profiles with genuine `ψ`-dependence in `ω` and `T`, touches it.

**THE +5% MUTATION TEST WAS RE-RUN ON THE ROTATING CASE AND REPRODUCES THE
STATIC RESULT EXACTLY.** Perturbing `RotatingSource::dFdPsi` by 5% leaves every
L2 error and every convergence rate unchanged **to all seven digits printed**, at
`k = 1, 2, 3`; Newton goes 3 iterations to 6, observed order 1.980 to 1.055, and
the assembled-Jacobian check 3.5e-11 to 2.1e-04. So the rate tables FL-4 rests on
are blind to the defect, and the Jacobian check and the order are what see it.

**THREE TRAPS FROM FL-5 TO FL-7, ALL OF THE SAME SPECIES.**

**A control on the reference curve is blind to the entire rotation chain rule.**
The gauge pins `φ₀(r_ref) = 0`, so at `r = r_ref` the exponent and *both* its
`ψ`-derivatives vanish and `∂²p/∂ψ²` collapses to `P₀″(ψ)` — the answer a
non-rotating source gives. A check that `∂F/∂ψ` varies with `ψ`, placed there,
read 0.333 against 7.300 at the outboard edge. **The one radius where the gauge
is exact is the one radius where a rotating source is indistinguishable from a
static one.**

**The best observed Newton order is not the order.** The bordered history opens
6.24e-02 → 4.01e-02 → 7.44e-03, whose "order" reads **3.81** — the iterate
walking into the basin. Asserting on the best triple would pass on that and keep
passing with a Jacobian degraded enough to destroy the tail. The assertion is on
the last triple above the round-off floor and is bounded **both** sides: 1 is a
broken Jacobian, 3.8 is an artefact.

**A comparison against an exact zero has no relative tolerance.** `φ₀` vanishes
identically on `r = r_ref`, and identically everywhere at `ω = 0`, so the
root-find-against-closed-form check compared 1e-33 with 0.0 and failed at every
such point. The fix is a floor at the problem's own energy scale, not a
case-dependent one read off the configuration — which was the first attempt and
failed again at `ω = 0`, where the scale is itself zero.

**MASCHKE & PERRIN IS A SECOND EXACT BENCHMARK, AND THIS FILE'S PLAN SAID IT WAS
NOT.** `FLOW-PLAN.md` rejected it as an adiabatic closure on the strength of a
`γ` in the equations. Wrong section of the paper: `refs/MaschkePerrin.pdf` —
*Plasma Physics* **22** (1980) 579, not the Phys. Lett. A 102 (1984) everyone
cites — carries two solutions, and its **§4** takes the temperature as a surface
quantity and is (136)'s isothermal closure. `γ` appears there once, only inside
`γΩ²`, and cancels out of the solution: **every `γ` works, not just `γ = 1`**.
Verified by substitution rather than re-derivation — 8e-26 relative at 50 digits,
exactly zero symbolically, `μ₀` restored by `p_SI = p_M&P/μ₀`. Its §3, the actual
polytrope, is a power law in `r²` and is **not** ours. The mistake is left
recorded because it is this project's standing hazard — three closures that look
alike — biting the plan that warns about it.

## Solution inversion: IN-A to IN-4 and IN-P are done; IN-5 and IN-6 remain

**`INVERSION-PLAN.md` is the design and the staged plan** — `ψ(R, z)` to
`R(Ψ, l)`, `z(Ψ, l)`, which is what `MANTA-COUPLING.md` needs, what
`DRIVER-PLAN.md` §3's `(Ψ, θ)` grid needs, and what `ROADMAP.md` item 10 is.
This section is only what a reader of the code needs.

**Nothing about the solve changes.** This is post-processing, in the same way
`FLOW-PLAN.md` was a change to `F` alone: a new consumer of `ψ_h` and `q_h`.

**AND `q` IS THE ASSET AGAIN, FOR THE THIRD TIME.** The band continuation of `ψ`
uses it, the band continuation of `B` uses it, and now the inversion does: a
critical point is a root of `q_h = 0`, which is a **solved** field converging at
the potential's own order, not a derivative of one converging an order down.
Compare CEDRES++, which records as an open problem that in P1 continuous
Galerkin the axis and the X-point are confined to mesh vertices, and TokaMaker,
which notes for Lagrange order ≥ 2 that saddles "can exist anywhere within the
mesh". meq resolves both sub-element by root finding.

The two pieces that exist are `src/meq/CriticalPoints.{hpp,cpp}` (stage IN-A)
and `src/meq/Zernike.{hpp,cpp}` (the basis IN-3 will fit in). Both headers carry
the long-form reasoning; what is below is what a maintainer needs to not
misread them.

### IN-A: the axis as a root of `q`, and a degree is never a count

**Measured against the analytic Solov'ev axis**, `CriticalPointConvergence.cpp`,
`k = 1, 2, 3` over `n = 4, 8, 16, 32`:

| | `k = 1` | `k = 2` | `k = 3` |
|---|---|---|---|
| position of the axis, rate over the whole sequence | **2.340** | **3.484** | **4.447** |

against design orders of 2, 3 and 4 — clearing `k+1` itself rather than `k+1`
less slack. **The per-pair rate is not a rate here** and the test says so: the
error is *pointwise*, so it oscillates as the axis moves within its element —
4.17 / 1.54 / 1.31 at `k = 1` — and the two-tier assertion pattern
`ExtensionConvergence` needs for the same reason is what this uses.

**The sharpest assertion in the stage was not in the brief**: the ratio of the
position error to `|q_h − q|` evaluated at the *exact* axis is **0.77 to 3.37**,
against **0.77 to 3.27** predicted by linearising `q` about its root. That is a
direct statement that the root finder adds nothing to the error of the field it
is rooting — the same shape as *A wrong Jacobian is invisible to a convergence
table*, from the other side.

**A DEGREE IS A SUM OF INDICES AND NEVER A COUNT, and the suite demonstrates it
rather than asserting it.** `audit()` walks the mesh boundary and accumulates
the turning of `q`, which by the Poincaré index theorem is the sum of the
indices of the interior zeros — `+1` for either extremum, `−1` for a saddle. A
box drawn round `iterExample2`'s axis *and* its X-point therefore reads
**winding 0 with two critical points inside**, the saddle located to **4.5e-6**
of the published X-point. Anything that reads a zero degree as "there is nothing
here" is wrong, and `theWindingNumberIsASumOfIndicesAndNotACount` is the live
demonstration.

**The Poincaré–Hopf hypothesis is transversality of `q·n`, NOT that the boundary
is a level set**, and the two boxes in the suite are the two cases. On the
standard benchmark rectangle — which is a level set of nothing — `q·n` keeps one
sign the whole way round at `min |q·n|/|q| = 0.15`, so `winding == χ == 1` **is**
a theorem there. On a box reaching past the X-point it reads **0.00** with a
sign change, and the degree is 0 against `χ = 1` with no contradiction whatever.
`IndexAudit::transverse` records which situation the caller is in, so the
comparison is not read as a theorem where it is a coincidence.

**AND A DISCONTINUOUS `q_h` CAN CARRY BOUNDARY DEGREE 1 WITH NO ZERO IN ANY
ELEMENT.** At `h = 0.4, k = 1` the audit reads 1 and an element-by-element
search finds nothing: each element's polynomial puts its zero just inside a
neighbour's territory, and with a face jump of `O(h^{k+1})` against an element of
size `h` there is a window where the zero belongs to neither. **Poincaré–Hopf is
a theorem about continuous fields**, and this is the DG jump meeting it head on.
The window closes with refinement, so it is a property of `h` rather than a
defect — and the test pins it to that one coarse mesh, so a *finer* mesh losing
the axis fails rather than passing quietly. The practical rule: where the audit
and the search disagree, **believe the audit**.

**`CriticalPoint::overshoot` is the same phenomenon at the level of one root**,
and refusing it outright does not work: at `k = 1, n = 4`, where the axis sits
beside the mesh line `z = 0`, the two candidates are 6.6e-4 and 8.9e-2 outside
their own elements and **with no allowance at all the axis is not found**.
`setContainment()` is in *reference*-element units on purpose, so the allowance
shrinks with the mesh exactly as the ambiguity it covers does. **It is not a
tuning parameter and that was checked**: swept over 0.001, 0.01, 0.05, 0.10 and
0.20 across the whole `k × n` benchmark the located axis is identical to every
digit printed.

### `CriticalPointFinder`'s axis is NOT `GradShafranovSolver::psiAxis()`

**They are different quantities, both correct, and neither should be changed to
match the other.** `ψ_ax` is *the largest nodal value* — chosen because the
bordered Newton needs a constraint it can differentiate, which under NPC makes
the border row exactly `−e_j`. IN-A's axis is *the point where `q_h` vanishes*.
They differ by `O(h)` in position and `O(h²)` in value, **both independent of
`k`**, so on a refined high-order mesh the two readings *separate* rather than
converge: measured on the finest Solov'ev mesh the gap is **202×** `ψ_h`'s own L2
error at `k = 2` and **4204×** at `k = 3`.

**`findAxis()` seeds from BOTH nodal extremes, and the reason is a sign error
this file's plan carried.** meq's `ψ` is not sign-normalised across sources:
the Solov'ev fixtures have `F` single-signed **negative**, so `ψ` is a
*subsolution*, its maximum is on `Γ`, and **the magnetic axis is an interior
MINIMUM** (Hessian determinant +0.693, trace +2.121 on `nstx`). With `F` positive
— the high-beta source — it is a maximum. So "seed from the largest nodal value"
is right for one sign of `F` and finds **a corner of the benchmark rectangle**
for the other. `AxisSense` is there for a caller who knows which they want, and
`findAxis()` refuses rather than guesses when both are present.

**Two rings of face neighbours, not one, and that is a measurement.** Where the
critical point sits on a mesh line the extreme nodal value is the shared vertex,
and which element is credited with it is decided by an L2 jump of 1e-8: on
`iterExample2` at `k = 2, n = 24` the minimum nodal value is in element 694 and
the root is in element **696**, which is not a face neighbour of it. The seeded
path is a fast path only — where it does not produce exactly one interior
extremum the search spends a full `sweep()`, so **the answer is not allowed to
depend on the seed**.

**Handing it the raw flux block instead of `flux()` is the failure to watch
for.** The raw block holds `−q`; in even dimension `index(−v) = index(v)`, so
every winding number is unchanged and **every Maximum silently becomes a
Minimum**. The audit still passes. `sweep()` is seeded Newton and is **not
exhaustive** — the certified subdivision of `INVERSION-PLAN.md` §5 is
deliberately not built, because IN-A's acceptance needs the axis and the audit
and neither needs exhaustiveness.

### The disc basis, and why `ρ = √Ψ_N` rather than `Ψ_N`

`src/meq/Zernike.{hpp,cpp}` is the basis `IN-3` fits in, landed early because it
is **MFEM-free** — plain doubles, like `Profiles` and `Source`, so CI can build
and test it without the MFEM branch it cannot obtain.

**The index constraint is the entire point.** `l − |m|` even and `l ≥ |m|` is
exactly what excludes `ρ² cos θ` and friends, which are not smooth at the
origin; what survives is that **every admissible mode is a bivariate polynomial
in `(x, y)`**, so the centre of the disc — the magnetic axis — is an ordinary
interior point and needs no special case. That is a better argument for the
basis than "DESC does it".

**GETTING THE RADIAL COORDINATE WRONG IS SILENT.** `ψ` has a quadratic maximum
at the axis, so `Ψ_N` behaves like (distance)² there and the geometry is smooth
in the *distance*: parametrise by `Ψ_N` directly and every basis converges
algebraically against a square-root branch point, worst near the axis, **with
nothing in a convergence table to say why**. Same species as *A wrong Jacobian
is invisible to a convergence table*. `radiusFromNormalisedFlux()` and
`fluxDerivativeFromRadial()` exist so the `1/(2ρ)` chain factor is not written
by hand at each call site, because that factor is the thing that gets dropped.

**The radial polynomial is a Jacobi polynomial under a change of variable**, so
`Zernike.cpp` calls `boost::math::jacobi()` rather than carrying a recurrence,
and **the explicit factorial sum every reference prints is not used at all**:
its terms are binomial-sized while its answer is `O(1)` — the largest is about
1e10 at `l = 30` — so it costs about `log₁₀(largest term)` digits to
cancellation. Measured at `ρ = 0.83` it disagrees with the Jacobi route by
**4.5e-8 at `l = 30` and 1.3e-4 at `l = 40`**, and it is kept in the tests as a
*control* rather than as an implementation. A flux-surface fit wanting twenty or
thirty modes would be reading noise. `find_package(Boost CONFIG REQUIRED)` and
`Boost::headers` are new dependencies of `meq_core`; the include is confined to
the `.cpp`, so a consumer of the basis takes on nothing.

**A DERIVATIVE CHECKED AGAINST A CENTRAL DIFFERENCE IS FLOORED BY THE
INSTRUMENT, NOT BY THE DERIVATIVE.** An earlier draft of the plan implied the
tolerance would tighten once the derivative was exact. It does not: a central
difference carries its own `O(h²)` truncation, so the comparison sits at
**1.3e-07** however exact the derivative is. Richardson extrapolation,
`(4D(h/2) − D(h))/3`, reaches **1.4e-11**. **That applies wherever this suite
checks a derivative against a difference**, which is several places.

### IN-0: the tracer, and the pairing that decides whether `q`'s tangent is worth anything

`src/meq/FluxSurfaces.{hpp,cpp}` is the predictor–corrector tracer and the
poloidal-angle parametrisation; `FluxSurfaceConvergence.cpp` is the acceptance.
**Fitted path only** — the band between `Γ_h` and `Γ` is deliberately not in it,
and `sampleField()` is the single seam it will be added at.

**THE DEFAULT IS `Potential::PostProcessed`, AND THE REASON IS NOT THE ONE THAT
WAS EXPECTED.** meq has two candidate pairs: `ψ_h` with `q_h`, both `k+1`; and
`ψ*` with `q*` from `DarcyForm::Reconstruct()`, where `ψ*` is `k+2`. Rooting
`ψ*` puts the traced curve `k+2` from the true one instead of `k+1` — measured
**60×, 54× and 83×** closer at `k = 1, 2, 3` on the same mesh, at *fewer*
corrector iterations per point.

**But the finding that matters is about the Hermite, and it was not predicted.**
The interpolant is built on tangents from the **flux** and measured against the
level set of the **potential**, and those are the same curve only so far as the
two fields agree. `∇ψ_h/r` agrees with `q_h` only to `O(h^k)` — differentiating
an L2 potential of degree `k` loses an order while `q_h` keeps `k+1`. So down a
`Δs` sweep the `ψ_h`/`q_h` pair is fourth order **until the tangent tilt takes
over and second order afterwards**: 3.809 → 1.400 → 1.569, tilt 3.5e-5 to
7.3e-5. The `ψ*`/`q*` pair has a tilt a full order smaller, 3.1e-7 to 5.9e-7,
and holds 3.960 / 3.996 / 3.812 across the whole sweep.

**So the post-processed pairing is what makes the cubic-Hermite-from-`q` claim
true at all on this discretisation**, rather than merely making it more
accurate. The control that identifies the tilt as the cause is a third column
built on `∇ψ_h`, the *exact* tangent of the curve being measured against, which
stays fourth order in both pairings — the right tangent for the representation
error and the wrong one for the field error, which is why all three are printed.

**AND THE USUAL REASON FOR PREFERRING `ψ*` IS WRONG.** It is not that the local
post-processing is built so `∇ψ*` matches `r q*`. MFEM's own documentation is
explicit: the constraint equation projects the **total** flux onto the face
restriction of an `RT_k` space, and *that* field is the source term of
Stenberg's local problem — so `∇ψ*` and `r q*` are different objects. What is
true, and measured, is that they agree an order better than the raw pair does.
**Right answer, wrong mechanism**, and the difference matters because the wrong
mechanism predicts exactness and would have made the measured tilt look like a
defect.

**THE FAILURE MODE THAT ENDS A TRACE, AND IT GETS COMMONER AS `Δs` FALLS.**
`{ψ_h = c}` is a union of per-element arcs offset by the face jump, so a point
landing within `jump/|∇ψ|` of a face is on **neither**: the Newton step computed
in element A pushes it into B, B's pushes it back, and the residual alternates
without ever meeting a tolerance tighter than the jump. The band is about 2e-5
wide on a contour of length 1.7 — rare, and *more* likely the finer the spacing,
simply because more points are placed. Left alone it ends the trace, and it
ended several before it was diagnosed. The corrector keeps its best iterate,
accepts after four non-improving steps, refuses to travel more than one
predictor step, and reports `stalledCorrections` against `correctorTarget`.
**One endpoint accepted at the jump level poisons that segment's interpolation
error by a factor of a hundred on a coarse mesh**, so the (c) measurement
excludes those segments as well as the face-crossing ones.

**Closure is not "machine precision", and the correct statement is geometric.**
Both the returning point and the start sit on the level set to 1e-13, but they
are separated *along* the curve by the final step's tangential offset `g_t`, and
an arc departs from its own tangent line by `κ g_t²/2`. Measured **3.573e-10
against a bound of 7.2e-10** — which is the discriminating assertion, because a
drifting tracer's normal error would grow with path length and have nothing to
do with `g_t`. The **residual** is what is flat with path length: 1.559e-13,
1.632e-13, 1.632e-13 over 1, 5 and 10 circuits. Do not compare the closure error
at 1 circuit against 10 — the step controller starts cold, so those two numbers
differ by 36× through their final steps rather than through their lengths.

**The rest, briefly.** `(c)` at fixed `Δs` is flat over a 16× change in dofs —
5.436e-09, 5.441e-09, 5.332e-09 — which is what says (b) and (c) are separated.
`(a)` against the exact contour converges at `k+1` on `Potential::Raw`. The
face-crossing jump converges at `k+2`. And the element walk —
`TransformBack`, then a breadth-first widening over face neighbours — reads
**zero** `FindPoints` fallbacks on every trace, which matters because
`CLAUDE.md` records `FindPoints` as `O(elements × points)` and a corrector calls
for a location once per iteration.

### IN-1: a spectral rule fed a second-order Jacobian is a second-order scheme

The trap `INVERSION-PLAN.md` §3.2 warns about, made into a measurement. Arc
length of one contour, three ways, rates in the number of angles:

| | |
|---|---|
| **from `q`, pointwise** | **7.03** |
| the trap: the same trapezoid, `ρ′` by differencing positions | 1.97 |
| chord sum | 1.99 |

**The rule is identical in the first two rows.** What differs is that `ρ′` comes
from `q` pointwise in one and from a central difference of neighbouring radii in
the other — and they converge to the **same limit**, so nothing in the second's
own output says it is orders worse. `ρ′ = ρ (u·t)/(u′·t)`, derived from
`ψ(a + ρ(θ)u(θ)) = c` and checked before use.

**Star-shapedness is a hypothesis, is measured, and is refused when it fails** —
`min|u × t|`, the same denominator, reading 0.844 / 0.823 / 0.805 on the
benchmarks. This is `IndexAudit::transversality` again, one stage later.

> **"SPECTRAL IN `N`" IS NOT ATTAINABLE ON A DISCRETE CONTOUR AND THAT IS NOT A
> DEFECT.** `ψ_h` jumps across faces, so `ρ(θ)` is piecewise analytic with jumps
> and no quadrature is geometric on it. The column plunges — 7.24, then 14.8 —
> and floors at about 1.2e-9, which is where the DG jump of `ψ*` converts to a
> distance (6.80e-10). **The control that says the floor is the field and not
> the rule** is the identical rule on the *analytic* contour, reaching
> 3.775e-15.

**`tests/analytic/FluxSurfaceReference.hpp` is the same statement from the other
side**, and it exists because IN-2 needs something to compare against.
Flux-surface averages on an *exact* equilibrium by rays plus the periodic
trapezoid: it reaches 3.11e-15 on the analytic contour, agreeing with the
tracer's own control, and it reproduces the trap at 1.37e-01 → 6.84e-04 for the
differenced metric while the pointwise one is at round-off by 256 angles.

**There are no closed-form Solov'ev flux-surface averages** — `ψ` is elementary,
but `V′`, `⟨R^{-2}⟩` and the safety factor are integrals over a *contour* of it,
and the Cerfon–Freidberg contours have no elementary arc length. What that file
supplies is a converged **reference value**, and the distinction should be said
wherever the number is printed.

**And it carries an identity that needs no reference value at all**: the
flux-surface average of the Grad–Shafranov equation,
`(1/V′) d/dψ ( V′ ⟨|∇ψ|²/R²⟩ ) = −⟨F/R²⟩`. **Quote it with its step or not at
all**: the residual is a property of the differencing step as much as of the
averages, running 9.6e-08 at 5% of `|ψ_ax|` down to 2.3e-11 at 0.6% on the exact
field, and on a *discrete* field it is not even monotone in the step — 2% reads
1.5e-08 where 1% reads 3.4e-07, because the difference divides the surfaces' own
DG-jump noise by the step.
Write the right-hand side with the `F` the solver is fed and **not** as
`−μ₀p′ − g g′⟨R^{-2}⟩`: the second is Solov'ev-specific and is re-derived by
hand, so it is not independent of the hand that derived it. **The `d/dψ` must be
Richardson-extrapolated** — a plain central difference floors the agreement at
8.1e-07 where `(4D(h/2) − D(h))/3` reaches 2.7e-13, which is `Zernike`'s
derivative finding for the third time in this file.

### The band, and the extension that was chosen by measuring both

`ContourTracer::setBandExtension()`. **`BandExtension::None` is the default**, so
the fitted path is bit-unchanged; the curved path is an explicit opt-in.

**THE DECIDING TEST TRACES `Γ` ITSELF**, which lies entirely in the band at every
mesh — `ψ_h` is strictly negative inside `Ω_h`, so every point of `{ψ = 0}` is
answered by the extension, and the exact answer is the curve `D_h` was cut from.
Worst distance from the exact `Γ` over `n = 16/32/64`:

| `k` | flux Taylor | transfer lift | lift closer at `n = 64` |
|---|---|---|---|
| 1 | 1.797 | **2.298** | **40×** |
| 2 | 2.138 | **3.995** | **1,610×** |
| 3 | 2.138 | **4.848** | **84,695×** |

**The flux Taylor step is second order at every `k` and that is structural**: its
remainder is `O(h²)` over a band of width `O(h)` however good `q` is, so it
cannot improve with the polynomial degree, and it does not — 2.138 at `k = 2`
and at `k = 3` alike. The transfer lift is the error of `q` **integrated along a
path of length `O(h)`**, which is `k+2`. `BandExtension::FluxTaylor` is kept as
the **control**, not as a fallback, and the tests say so in their failure
messages.

**AND THE PRIMITIVE THE PLAN NAMED CANNOT DO IT.** `mfem::PathLiftCoefficient`
`dynamic_cast`s its `ElementTransformation` to `FaceElementTransformations` and
lifts from *that face's own* integration point — it answers "what is `φ_h` on
`Γ_h`", which is `η₅`'s question and not this one. The usable primitive is one
level down and public: **`mfem::PathIntegral( Cu, x, xbar, line_ir )` takes
arbitrary endpoints**, with `mfem::ElementExtension` supplying `E_h(q_h)`.
MFEM's own comment is the licence — fed the exact flux it must return
`p(x) − p(a(x))` *"whatever the path"*.

**THREE TRAPS IN THE MEASURING, AND THE FIRST NEARLY PRODUCED A WRONG HEADLINE.**

* **A contour at fixed `Ψ_N` cannot measure the extension's order.** As `h`
  falls, `Γ_h` climbs toward `Γ`, so a contour a fixed distance inside has its
  band excursion shrink *faster* than `h` — `deep/h` goes 0.92 → 0.66 → 0.26 and
  the band population 136 points → 9. Both columns then converge faster than the
  extension beneath them, the Taylor step reading 3.75 at `k = 3` against its
  true 2.1. Tracing `Γ` is what fixes it.
* **The nearest face of `Γ_h` is not the face you are outside of.** A staircase
  `D_h` cut from a diagonally split Cartesian mesh **pinches** — two triangles
  meeting at one vertex — so both lobes' faces are equidistant from a point just
  outside and the tie goes to loop order. Half the time that picks the lobe
  whose outward normal points the other way and a genuine band point is refused:
  a trace stopped after 85 points of 320. **The outward test must be a filter
  applied first**, with the nearest taken among the survivors.
* **`ψ_h` is not strictly negative inside `Ω_h` to machine precision**, so a
  handful of `{ψ = 0}` points land back inside — one of 322 at `k = 1, n = 64`.
  That is the discretisation creeping above its own imposed datum near `Γ_h`,
  and the acceptance asserts ≥ 98% band rather than 100%.

**The flag is per point and every mask is asserted against `FindPoints`**, in
both directions, because `CLAUDE.md` records that a *count* rather than a mask
was the other half of a real defect in the `.nc`. `ContourPoint::extended` and
`bandDepth`, `Contour::extendedPoints` / `deepestBandPoint` / `bandExtension`,
and `AngleParametrisation::extended` per node, so a consumer can report
band-crossing surfaces separately.

### IN-2: the averages, and `ψ*` does NOT buy an order here

`src/meq/SurfaceAverage.{hpp,cpp}`. **One primitive with two builders** — an
integrand in `(R, z, ψ, q)` over either the angle parametrisation or the traced
contour with Gauss points — and every named quantity a one-line wrapper. That
shape is deliberate: `MANTA-COUPLING.md` says the slot list "is negotiated with
the transport physics case", so a list that is still moving must not be a list
of functions. Conventions, which are part of the definitions:
`V′ = ∮2πR dl/|∇ψ|` and `⟨X⟩ = (1/V′)∮2πR X dl/|∇ψ|`. **The safety factor is
`safetyFactor(g)` and never `q`**, `q` being the flux.

**`ψ*` DOES NOT BUY `k+2` IN AN AVERAGE AND THE REASON IS STRUCTURAL.** Both
pairings converge at **`k+1`** — `V′` at 2.230 / 3.185 / 4.296 raw against
1.957 / 2.827 / 4.265 post-processed. The weight is `2πR dl/|∇ψ|` and
**`|∇ψ| = r|q|`**: the reconstruction buys its extra order in the *potential*,
and there is no `k+2` flux to divide by. The level set improves, the weight does
not, and the average inherits the worse. What `ψ*` buys is a **constant** —
×1.29, ×1.46, ×1.73 in `V′`.

**That is the same shape as the band continuation of `B`**: a quantity limited
by the one factor with no solved variable behind it. Third time in this item
that the answer turned on *which field the error divides by* rather than on
which field was rooted.

**AN AVERAGE DOES NOT ESCAPE THE METRIC TRAP, AND THE PLAUSIBLE ARGUMENT THAT IT
DOES IS WRONG.** A ratio looks as though it should cancel a bad metric, the same
weight appearing above and below. It cancels a **constant** — about 40× — and
**nothing in the order**: from `q` pointwise the sequence rates are 6.77 for `V′`
and 7.38 for `⟨R^{-2}⟩`; with the metric differenced, 1.92 and 1.80, a separation
of 1.9e+06 by 128 angles.

**THE IDENTITY'S CONTROL IS FLAT, WHICH IS THE SHARPEST FORM OF THE RICHARDSON
FINDING YET.** At `k = 2` over a *sixteenfold* refinement the plain central
difference reads 5.8e-05, 7.9e-06, 8.7e-06, 8.9e-06 — **it stops moving at three
figures** — while Richardson goes 8.9e-05 → 1.5e-08 at rate 4.185. A column that
does not converge under mesh refinement is measuring the instrument.

**And the step is not monotone on a discrete field.** At `k = 2, n = 96` a step
of 2% of `|ψ_ax|` reads 1.5e-08 where **1% reads 3.4e-07** — the smaller step 20×
worse, because the difference divides the surfaces' own DG-jump noise by the
step. There is an optimum; find it rather than assuming smaller is better.

**Two extractions agree to about their own error and CANNOT do better.** The
angle fit and the Hermite contour disagree by an amount converging at the
field's own order — 2.96/3.09 at `k = 2`, 4.14/4.48 at `k = 3` — because
`{ψ_h = c}` is a union of per-element arcs offset by the jump and two routes
placing nodes differently sample different arcs. **So the assertion is on the
rate**, not on the gap being small. A missing `2πR`, a metric about the wrong
point, or a gradient on the wrong side of the division would each break the
rate, and no single-route table could see any of them.

**The fixture needed its own box, and one surface is not measurable at all.**
`standardBox()` cannot hold these surfaces — `Ψ_N = 0.25` on `nstx()` already
spans `r ∈ [0.99, 1.57]` against a box ending at 1.4 — so the study runs on
`[0.60, 1.90] × [-1.10, 1.10]`. `Ψ_N = 0.75` leaves under one cell of margin at
the coarsest mesh of a dyadic sweep, so an `h`-study there measures the
contour's distance to the mesh boundary; its reference value is asserted
instead. **That is the fixture's elongation, and it is one more reason the
curved path is where this item actually lives.**

### IN-3: the fit, and the angle is not free

`src/meq/SurfaceFit.{hpp,cpp}` fits `R` and `z` as truncated Zernike expansions
to a point cloud of traced surfaces. **MFEM-free**, like `Zernike`, `Profiles`
and `Source` — plain doubles in, coefficients out — so CI can build and test it;
a caller writes a two-line loop to turn an `AngleParametrisation` into samples,
and `CMakeLists.txt` records that reason beside the source list.

**THE GEOMETRIC POLOIDAL ANGLE IS INADMISSIBLE, AND EVERYTHING ELSE IN THE STAGE
IS DOWNSTREAM OF IT.** Labelling surfaces by the geometric angle about the axis
— the obvious choice, and what `AngleParametrisation` produces — **makes the
disc map non-smooth at the axis for any non-circular surface**, so no basis that
is smooth there can converge against it.

The argument is exact, and it takes one line to check. A function smooth at the
origin has precisely **one** angular harmonic multiplying `ρ¹`. Take nested
ellipses of semi-axes `a`, `b` labelled by geometric angle; with `u = ρ cos θ`
and `v = ρ sin θ`,

```
x( u, v ) = a b · u · √(u² + v²) / √(b² u² + a² v²)
```

whose trailing factor is **homogeneous of degree zero** — a function of
direction alone, with no limit at the origin unless `a = b`. The `ρ¹`
coefficient therefore carries `cos θ, cos 3θ, cos 5θ, …`, and every harmonic
past the first is one the Zernike index constraint **excludes**, precisely
because it is not smooth there. The parametrisation puts content exactly where
the basis refuses to look. Measured on nested ellipses, where the answer is
known exactly, that column decays like `L^{-1.2}` and **never converges**:
1.53e-01 at `L = 2`, still 1.07e-02 at `L = 20`, against 1.8e-14 relabelled.

**The repair needs no field.** Near the axis every equilibrium's surfaces are
ellipses, so a three-parameter fit of a quadratic form to the *innermost traced
surface* recovers tilt and flattening, and the relabelling is an exact
reparametrisation of the circle. `meq::relabelByAxisShape`, worth **45× to
660×**. On `nstx()` it recovers short/long = 0.4943 from samples against 0.4845
from the Hessian — two independent routes to the same shape.

**`ρ = √Ψ_N` IS NOW MEASURED RATHER THAN ARGUED**, which is what the control was
for. Worst fit error **3.44e-05 against 7.01e-03** for the same points
parametrised by `Ψ_N`, a factor of 204; the envelope over `l = 10 → 20` falls by
32.6 against the control's 4.5. **Conditioning is untouched by the choice**
(1.87e3 against 1.44e1), so it is not a conditioning artefact — it is the
square-root branch point at the axis, exactly as `Zernike.hpp` claims.

**The parity control fits the sample cloud eight times better and is useless**:
condition number **9.15e+16**, axis error 1.23e-04 against 3.56e-07, and an axis
that moves by 4.8e-03 depending which `θ` you approach along. A better residual
on the data you fitted and nothing anywhere else.

**The axis comes out `θ`-independent for free and exactly** — spread
`0.000e+00` at every degree, on analytic and discrete data alike, because every
mode with `m ≠ 0` carries `ρ^|m|` and above. Asserted as an exact zero, not a
tolerance. **A caveat for anyone rebuilding that control**: a tensor product
that keeps `l ≥ |m|` *also* gives an exactly `θ`-independent axis, because those
modes vanish at `ρ = 0` too. Only admitting modes with **no radial factor**
breaks it, and a control built the obvious way demonstrates nothing.

**THE LEVER ON CONDITIONING IS THE HOLE, NOT THE SAMPLE LAYOUT.** A sample set
has a hole in the middle — no surface is traced at `Ψ_N = 0` — and the
orthogonality argument needs nodes spanning the whole disc. Condition number
against the inner limit: 7.78 at `Ψ_min = 0.02`, 3.19e+02 at 0.10, **7.30e+04 at
0.25**. Four orders. The three layouts — equispaced in `Ψ_N`, equispaced in `ρ`,
Gauss in `Ψ_N` — agree to within **25%** at every hole size, and Gauss is
sometimes the worst. The test asserts that wrong story dead, at
`layout spread < 2×`.

**And rescaling the disc edge is a change of BASIS, not of model.** A Zernike
expansion of degree `L` spans the polynomials of degree `L` in `(x, y)` and that
space is closed under scaling, so the two extents fit the *same function* —
measured, identical worst errors to better than 1e-6 relative at every inner
limit. The choice is purely conditioning and is worth up to **11,500×**. The
default stays `1.0` so that a change of coordinate is never silent, and
`majorRadiusExpansion()` refuses unless basis, coordinate and edge are all the
plain ones.

**The Richardson finding needed a step sweep to appear at all, and that is
itself the lesson.** At the natural step it is worth **1.6×**, not the 10²–10⁴
seen elsewhere in this tree, because here the *fit's own* derivative error is the
binding constraint rather than the instrument. Swept, it separates properly —
70.5× at a step of 0.16 — and the diagnostic is that **the plain column falls by
50.8 across an eightfold refinement while the extrapolated one moves 11%**: the
converging column is the instrument, the flat one is the answer.

**The number a coupling reads.** `∂(geometry)/∂Ψ_N` grows exactly as the
coordinate demands: the product with `2√Ψ_N` settles at **1.148** as the inner
limit falls to 0.005, so an innermost node at `Ψ₁` carries a geometry derivative
of about **0.574 / √Ψ₁**. The assertion is that the product is bounded and
settling, which is the statement that the growth belongs to the coordinate and
not to the fit. `MANTA-COUPLING.md`'s consumer can keep its nodes off `Ψ = 0`,
so this is a conditioning number for a coupling to read rather than a defect.

**`fitByAngle()` now accepts its best iterate** where the tolerance is
unattainable, with `AngleParametrisation::stalledRays` beside `worstResidual`,
and keeps the throw for a ray that never brackets. The acceptance asserts against
the **measured face jump** rather than a chosen tolerance: 11 of 1536 rays
accepted at their best, worst residual 6.97e-04 against a DG jump of 1.63e-03 —
a ratio of 0.43, i.e. as close as the field allows. `AngleParametrisation` also
keeps the `q` it was computing and discarding, and `sampleAt()` has an overload
reporting `extended`.

### IN-4: the answer is that the question was the wrong one

**The `ψ`-element decision is closed and no element was built.** `§4.4`'s three
candidates were three ways to buy a `ρ`-dependent poloidal angle; **solving for
the angle directly is cheaper than all of them**, and it is what DESC was
measured doing.

`meq::gaugeFreeFit()` requires each disc node only to **land on the right
surface** rather than to sit at a prescribed angle — a geometric Gauss–Newton on
`Ψ_N( x(ρ,θ) ) − Ψ`, warm-started from IN-3's linear fit, with `∇Ψ_N` from the
**solved flux**, `∇ψ = r q`. No force balance and no second solver: meq already
has `ψ`. Rows are scaled by `1/|∇Ψ_N|`, so the residual is a **distance in
metres** and the error measure is itself gauge invariant. `SurfaceFit` stays
**MFEM-free** — the field arrives as one callable returning value, gradient and
a refusal.

**On nested ellipses, where the answer is known exactly**, started from the bad
fit: **1.718e-01 → 8.31e-16** at `L = 2`, and 1.623e-02 → 1.86e-09 at `L = 16`.
**That is the theoretically right answer and not a lucky one**: for nested
similar ellipses `x = a u`, `z = b v` in the disc's own Cartesian coordinates, so
the family **is a degree-1 map** under the correct angle, and the solve finds it
in four iterations at round-off. The prescribed-angle fit decays at `L^{-1.2}`
and never converges.

**IN-3's algebraic tail is gone and it holds to `Ψ_N = 0.005`** — linear
1.67e-02 → 3.78e-04 against gauge-free 1.78e-03 → **4.16e-09** over
`L = 4 → 16`, the difference concentrated in the last leg (1.69× against 52×).
**No inner limit tried costs it more than a factor of 1.5.** That is the number
`MANTA-COUPLING.md` needed, and it removes the constraint that motivated
elements in the first place.

**Panici's Figure 5 reproduced on a fit, with no penalty at all**: spectral width
`M(2,2)` falls 1.607 → 1.366 under minimum-norm damping alone, and 1.481 → 1.265
on the discrete field. A solve slides the angle to what its basis represents
best without being asked.

**THE EXPLICIT SPECTRAL-WIDTH PENALTY LOSES ON ITS OWN METRIC.** `M(p,q)` is a
**ratio** of two weighted sums of the same coefficients, so a quadratic penalty
is not a surrogate for minimising it: twelve decades of `λ` move `M` by 1.6% *in
the wrong direction* and cost **44×** in surface error. Hirshman & Breslau
minimise `M` itself, which is not a quadratic problem. Kept as the losing column.

**THE GAUGE IS A SOFT TAIL WITH NO GAP, NOT A NULL SUBSPACE.** Measured, the
ellipse family has **exactly 3** null directions at every degree from 2 to 16 —
of 6 to 306 columns — and **`nstx` has none at all**. What both have is a smooth
tail running to 8e-08 of the largest singular value, with 58 of 306 directions
below `1e-4 σ_max` on the ellipses and 97 on `nstx`. So there is nothing to
project out, and **the floor is a threshold that has to be chosen rather than
read off a gap**. The no-gauge control fires on both fields — first step
`8.6e+10 ×` the coefficient norm, Jacobian `−2.2e+13`, i.e. **folded** — but the
mechanism and the magnitude differ by six orders between them, so a control
measured on one field would have reported whichever it happened to meet.

**AND THE TRUST REGION IS ITSELF A GAUGE, WHICH THE CONTROL HAS TO KNOW.** The
undamped pseudo-inverse reaches round-off at `L ≤ 12` and **fails at `L = 16` and
20**; what makes it robust is adaptive Levenberg–Marquardt damping, which is a
Tikhonov term in disguise. `SurfaceGauge::None` therefore disables the damping as
well as the floor — **a "no gauge" control that kept the trust region would pass
while testing nothing.**

**THE MAP IS CHECKED FOR FOLDING.** Minimum Jacobian **+5.2e-02 to +6.3e-02**
with the gauge on, negative on **every** ungauged run. A surface residual alone
admits a beautiful number over a folded map, which is exactly the class of quiet
wrong answer this file exists to catalogue.

**ONCE THE ANGLE IS FREE, ANY ACCEPTANCE WRITTEN AGAINST A PRESCRIBED ANGLE
MEASURES THE GAUGE RATHER THAN THE FIT.** IN-3's derivative and metric checks
compare the fit's position *at a given `θ`* against a surface traced at that `θ`
— precisely the freedom being granted — so they are **not** re-asserted for the
gauge-free fit. They are replaced by gauge-invariant properties: distance to the
surface, fitted perimeter against exact (3.5e-07 relative), and the sign of the
map Jacobian. **This costs the consumer nothing**, and that is the point:
`MANTA-COUPLING.md` reads flux-surface *averages*, and an average does not know
how its surface was parametrised. The deliverable was gauge-invariant all along.

The conversion product is unmoved at 1.156 against 1.148, the axis spread stays
an exact `0.000e+00` at every degree on both fields, and on the discrete field
the fit sits **5.71e-07** from the exact surfaces — the post-processed pairing's
own `O(h^{k+2})`, so what remains is the discretisation and not the
representation.

### IN-P: what the inversion actually costs, and the cost model was wrong

`tests/performance/InversionScaling.cpp` and `inversion-scan.sh`. **Not a
ctest**, per the standing rule that every number in `tests/performance/` is a
timing and a threaded timing on this machine is a measurement about the machine;
it exits non-zero only for the correctness properties that make its timings mean
anything. Note that **`naming` does not cover `tests/performance`** — it runs
over `MEQ_CORE_SOURCES_PRESENT`, which is `src/meq` alone.

**85% OF THE CHAIN IS `gaugeFreeFit`, AND `INVERSION-PLAN.md` §11 DOES NOT
MENTION IT** — that section was written before IN-4 existed and puts the weight
on the tracer's per-point corrector instead. Measured at `k = 2`, `n = 48`, 12
surfaces × 48 angles, `L = 10`, serial: `trace` 11.2%, `fitByAngle` 1.8%,
`surfaceAverages` 1.8%, the linear `SurfaceFit` 0.6%, `findAxis` 0.01%, and
`gaugeFreeFit` **84.6%**. For scale the solve itself is 0.37 s and
`postProcess()` 0.62 s.

**And inside `gaugeFreeFit`, 82.5% is field evaluation** — 21,888 `sampleAt`
calls — so **about 70% of the whole chain is `ContourTracer::sampleAt`** and the
linear algebra is nowhere near the cost. §11 predicted the quadrature would not
be where the time is, and that much is right: it is 1.8%.

**The corrector is not the problem either.** 97.3% of accepted points take
**exactly two** iterations, mean 1.994, worst 3, with `stalledCorrections` and
`fallbackLocations` both zero on every trace.

**Extraction cost is independent of `k`** — 0.140, 0.141, 0.140 s at `k = 1, 2,
3` — and grows like `1/h`, because the step ceiling is a fraction of the element
size. **In the angle count it is not monotone: more points can be cheaper**,
0.048 s at 24 angles against 0.023 s at 192, because closer spacing makes the
element walk hit instead of missing.

### The largest lever is one integer, and it is free

**`sampleAt` decomposes as 26% walk, 0.9% evaluation, and 73%
`Mesh::FindPoints`** — the last-resort fallback, taken on 184 of 576 points.
Four rings of *face* neighbours reach about four triangles in a straight line and
fewer diagonally, while consecutive ray nodes are one to two cells apart. Swept
on the real code path, with the answers **bit-identical at every depth**:

| `setWalkDepth()` | seconds | fallbacks |
|---|---|---|
| 2 | 0.0895 | 367 |
| **4, the default** | 0.0531 | 183 |
| 8 | 0.0482 | 30 |
| 12 | **0.0262** | **0** |
| 16 | 0.0255 | 0 |

Worth **2.0× on `fitByAngle`, 2.1× on `surfaceAverages`, 1.8× on
`gaugeFreeFit`**, and about **1.57× on the whole chain**, for no change to any
answer. The default is deliberately left at 4 pending a decision, because it is a
library default and this project protects answers rather than timings — but the
harness asserts the bit-identity, so the change is provably free.

**And it is not only a speed question.** With no fallbacks there is no call into
`FindPoints`, which is what makes threading impossible today. Fixing the depth
and giving `traceFromAxis()` its axis element as a hint — `findAxis()` already
knows it — would remove the last unconditional `FindPoints` call and unblock the
shared-tracer construction below.

### Threading: available, and blocked by one non-reentrant function

Given a mesh, a solve and a tracer per thread, the parallelism §11.2 predicted is
there and is exact:

| threads | over surfaces | over rays |
|---|---|---|
| 2 | 1.89× | 1.89× |
| 4 | 3.16× | 3.36× |
| 16 | 3.31× | **7.55×** |

**Every count reproduces serial at `0.000e+00`** — contour points, ray radii,
quadrature weights and `V′` — which is available exactly because independent
surfaces and independent rays reassociate nothing. Surfaces cap at about 3.3×
because there are only twelve of them and the outer ones are longer; rays scale
properly, which is §11.2's asymmetry confirmed: **the tracer's steps are
sequential and the rays are not.**

**What is not available is sharing one tracer**, for the `FindPoints` reason
under *Traps*. That is the single blocker, and the walk-depth item above is most
of its cure.

### Two levers from §11.4, one of which does not exist

**Continuation in the flux label is worth 1.004× — nothing — and corrector
iterations went UP.** §11.4 framed this as a genuine trade against parallelism
over surfaces and warned against reasoning from structure. **There is no trade to
resolve**: the predictor for every point after the first already comes from the
previous point of the *same* surface, so continuation can only save
`traceFromAxis`'s bracket and one corrector per surface. Take the parallelism.

**The per-`ψ` cache MaNTA needs is worth `nodes/surfaces` and nothing eats it**:
60 physics nodes over 3 residual evaluations cost 6.55 s naively and 1.29 s
served from one family per `ψ`, a factor of **5.1**.

**AND HERE IS THE NUMBER THAT DECIDES THE COUPLING'S DESIGN.**
`dGeometry_dpsi` by differencing costs one extraction per `ψ` degree of
freedom: **0.447 s × 46,080 dofs ≈ 5.7 hours for a single Jacobian, serial.**
That is what §11.4's shape derivative has to beat, and it is the only part of the
chain where a core count is worth a factor rather than a few percent.

### Eigen, and what the swap did and did not move

`SurfaceFit.cpp`'s hand-rolled Householder QR and one-sided Jacobi SVD are gone,
replaced by `Eigen::JacobiSVD` at both call sites. `find_package(Eigen3 3.3
REQUIRED NO_MODULE)`, **REQUIRED rather than optional-with-fallback**, because a
second numerical path is exactly the maintenance burden the swap exists to
remove; `PRIVATE` to `meq_core` with the include confined to the `.cpp`, so
`SurfaceFit` stays MFEM-free *and* Eigen-free to its consumers.

**`JacobiSVD` and not the faster `CompleteOrthogonalDecomposition`**, for two
reasons: `gaugeFreeFit` applies `σ/(σ² + μ)` and so needs the singular values and
`V` explicitly, and the diagnostics report the design matrix's spectrum — a
rank-revealing QR gives a rank and a solve and no spectrum. It is also the same
algorithm meq had, with the same high *relative* accuracy in the small singular
values, which IN-4's soft tail needs.

**`EIGEN_DONT_PARALLELIZE` is set, and for a better reason than the one it was
asked about.** Read from the source, `Parallelizer.h` bails to the sequential
path when `omp_get_num_threads() > 1`, so Eigen does **not** nest inside an
OpenMP region and `setNbThreads(1)` is unnecessary for that. It is pinned because
*outside* a parallel region Eigen would take `omp_get_max_threads()`, and a
threaded GEMM blocks differently at different thread counts — **meq's answers
would depend on `OMP_NUM_THREADS`.** `EIGEN_USE_BLAS` and `EIGEN_USE_LAPACKE`
are not set, per the MKL rule.

**IN-4's numbers survive and two moved, reported rather than re-baselined.** The
nstx tail reads 4.163719e-09 against 4.163718e-09, the axis spread is still an
exact `0.000e+00`, and the minimum Jacobian is unmoved to five figures. The
ellipse `L = 2` figure moved 8.31e-16 → 6.80e-16 — both round-off, the headline
ratio reading 2.53e+14 instead of 2.07e+14. **The no-gauge control moved by a
factor of ten** and that is expected: it inverts singular values down to exactly
zero, so its step is dominated by the smallest resolved `σ` and two SVDs report
that differently. Its verdict is unchanged — still folded, still diverging.

### Three things in the kernel survey worth not re-doing

| | measured | |
|---|---|---|
| Zernike by Kintner recurrence against per-mode evaluation | **6.3×**, agreeing to 4.4e-15 | **not worth it** — 11.5% of a 0.6% stage |
| a fit at 20,000 points by Vandermonde + GEMM against per-point | **34.7×**, agreeing to 4.4e-16 m | **worth it for IN-6**, which is the many-point case, and not for anything today |
| batched field evaluation by element | **0.90–1.14×** | **do not** |

**And the reason the third one fails is worth keeping**, because the brief
asserted the opposite: `GridFunction::GetValues( i, ir, vals )` is **a loop
calling `CalcShape` per point, not a GEMM**. It amortises the dof gather and the
temporaries and nothing else. The access pattern is 1.15 points per element in
any case, so there is nothing to batch.

### What is next

**IN-A, IN-0 (both halves), IN-1, IN-2, IN-3, IN-4 and IN-P are done.** What
remains is IN-5 and IN-6, and a short list of things the stages left behind.

**IN-5, open surfaces**, is deferred with free boundary per §6. The canonical
in-surface coordinate for it is **poloidal arc length normalised to `2π` from a
fixed-`z` reference ray**: a disc chart has no meaning through a separatrix and
an angle about the axis has none on an open line. Note that arc length does
**not** fix axis regularity — for similar surfaces it is a `ρ`-independent
relabelling and buys the same one order — which is why the two concerns are
handled by separate machinery there.

**IN-6, the output**, is `DRIVER-PLAN.md` §3's `(Ψ, θ)` grid plus the per-`ψ`
cache `MANTA-COUPLING.md` §5's pointwise call pattern requires. **Both of its
numbers are already measured**: the cache is worth `nodes/surfaces`, 5.1× on a
60-node case, and evaluating a fit at many points by Vandermonde-plus-GEMM is
worth **34.7×** — a stage that is 0.2% of the chain today and is exactly the
many-point case IN-6 creates.

**And the number that decides the coupling's design**: `dGeometry_dpsi` by
differencing is **about 5.7 hours for one Jacobian, serial**. The shape
derivative has to beat that, and it is the only part of the chain where a core
count buys a factor rather than a few percent.

**Open, small, and each found by the stage after the one that caused it:**

* ~~**`setWalkDepth()`'s default of 4**~~ — **RAISED TO 12, 2026-09-04.** Four
  rings is enough for `trace()`, which steps a fraction of an element, and not
  for the *rays* of a parametrisation, which are placed by angle and land one to
  two cells apart: depth 4 took the `FindPoints` fallback on 183 of 576 rays and
  depth 12 on none. Worth about **1.57× on a whole extraction**, and the answers
  are bit-identical at every depth — the walk decides how a point is *found*,
  not where it is. It is also most of the cure for the threading blocker, since
  a shared tracer aborts the moment any thread reaches `FindPoints`.
* **`fitByAngle()` throws where the corrector would accept.** Its ray Newton
  demands a tolerance that on a discontinuous field is sometimes *unattainable*
  — a ray crossing a face where `c` falls inside the jump has no point on it
  with `ψ_h = c` at all — while the tracer's own corrector already handles that
  by keeping its best iterate. At `k = 1` on the raw pairing it fails on every
  mesh from `n = 12` to 32.
* ~~**Two headers describe a caller that was never moved.**~~ — **FIXED,
  2026-09-04, and the fix was a correctness one as well as a tidying.** The
  contour builder now marks **each Gauss node for itself** through the
  seven-argument `sampleAt()`. Marking a whole segment from its endpoints was
  conservative in one direction and **wrong in the other**: the band does not
  respect the segment a node sits in, so a segment can have both endpoints
  inside `Ω_h` and still cross `Γ_h` in between, which **under-reported** — a
  band quantity presented as a solved one. And the fit builder now samples
  **nothing at all**: `AngleParametrisation` keeps the potential, the flux and
  the band flag its own ray Newton found, so the averages read them instead of
  re-deriving them. `potential` was added beside `fluxR`/`fluxZ` for this, and
  it is deliberately *not* the surface's level — on a stalled ray the node sits
  as close to the level as the field's jump allows and no closer, and this
  records where it actually is.
* **§3.3's implicit quadrature is the missing third leg** of IN-2's
  cross-check. Its acceptance said "all three agreeing is worth more than any one
  being plausible" and two were delivered.
* **Maschke & Perrin is a verified exact rotating benchmark** — 8e-26 by
  substitution — and is still not in `tests/analytic/`.

**What is deliberately absent from the tracer**, so nobody reads more into it
than is there: it follows **one connected component** and neither finds nor
reports a disjoint island at the same level — the same disclaimer
`CriticalPoints.hpp` makes about seeded Newton, and for the same reason. It is
**not X-point aware**: the level set through a saddle is not a 1-manifold and
the tangent is undefined there, so a trace at that level will stall or turn a
corner arbitrarily. `pointAtArcLength()` parametrises segments linearly in
*polyline* length rather than true arc length, and says so.

## The linear solves, and what they should be

A hybridized HDG scheme needs exactly two linear solvers: one for the global
face-coupled trace system, one for the small dense per-cell systems. meq has a
**third**, and it is not an oversight — it is the price of Newton.

| | what meq uses | |
|---|---|---|
| global trace | **selectable**: `UMFPackSolver` (METIS ordering, the default), `PardisoSolver` (`REAL_STRUCTURE_SYMMETRIC`), `CuDSSSolver` (`NONSYMMETRIC` + `FULL`) | three unsymmetric sparse LUs, agreeing to 5e-14; `setTraceSolver()` |
| *fallback, no direct solver at all* | GMRES | **unpreconditioned** on the Newton path, `GSSmoother` on the linear path |
| per-cell dense | MFEM `LUFactors`, partial-pivot LU | as a 2×2 block: LU on the flux block `A`, local Schur `S`, LU on `S` |
| **per-cell nonlinear** | element-local `NewtonSolver`, 100 iters, rtol 1e-12, `LPrecType::LU` | one iteration *per element per residual evaluation* — **and meq's default ordering no longer reaches it**, see below |

And a **fourth** when `ψ_ax` is an unknown, which is not a fourth solver: the
bordered system is solved by block elimination against the same factorisation of
the trace Jacobian, so it costs one extra backsolve per Newton step and nothing
else. Assembling the border into an `(n+1)` matrix would put a dense row and a
dense column into the factorisation for no gain, which is why it is not done.

**The third one exists because meq uses Newton rather than Picard, and
`NonlinearOrdering::NPC` removes it.** Picard evaluates `F` at the previous
iterate and leaves every local problem linear — one dense factorisation and done.
Condense-first Newton puts `∂F/∂ψ` inside the local problem and makes it
nonlinear. NPC, which is **meq's default**, differentiates the full
`(q, ψ, ψ̂)` system first and hybridizes the linear system that results, so every
element-local operation is a linear solve and `GetNumLocalNLIterations()` stays
at **zero** — asserted, not assumed, by
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally`, which reads 0
against the condensation's 3644/3560/3412 at `k = 1, 2, 3`.

So the row above is what meq *has*, and it is one `setNonlinearOrdering()` call
away rather than unreachable — `CondenseThenLinearise` is the **backup**, and on
stiff under-resolved meshes it is the one that works. See *The NPC port*.

**And there is a fifth thing that is not a solver at all but decides how fast
the third and fourth run**: `setAssemblyMode()`, which threads the element loop
that builds all of them. Measured under *Threading, measured*.

**The trace solver became a run-time choice on 2026-08-31**, which it had not
been: `setTraceSolver()` picks among whichever of the three the build has, and
`traceSolverAvailable()` says which those are without throwing. All four
construction sites — Picard, Newton, the bordered Newton and the linear path —
now go through one factory, so the *only* thing that differs between them is
whether the symbolic analysis is retained. It is, everywhere except the linear
path, which factorises once and destroys the object.

**Choosing one is never a numerical decision.** `theTraceSolversAgree` drives
each of them through the whole solver and pins the recovered `ψ` against the
first — 5.0e-14 for PARDISO against UMFPack — and cuDSS is checked the same way
by the `TraceSolverCuDSS` ctest, which is separate only because cuDSS needs an
`mfem::Device` and that is global state a Boost test case should not be
configuring.

### The trace matrix is symmetric on the fitted path and is not on the extension path

**And that is not a defect — it is what the extension technique is.** This is the
single most important thing in this section, and it was nearly missed by
measuring only the easy configuration.

`GradShafranov.cpp` justifies GMRES over CG with "the hybridized trace system is
small but not symmetric positive definite in this sign convention". Measured on
the **fitted** benchmark that is true as written and misleading:

| `k` | `n` | nnz/row | rel `\|A_ij − A_ji\|` | Rayleigh quotients, free subspace |
|---|---|---|---|---|
| 1 | 416 | 9.4 | 2.2e-16 | 200/200 negative, largest −1.44 |
| 2 | 624 | 14.1 | 3.6e-16 | 200/200 negative, largest −1.30 |
| 3 | 832 | 18.8 | 6.9e-16 | 200/200 negative, largest −1.28 |

**Symmetric to round-off and negative definite** — so on a fitted mesh `−A` is
SPD and every symmetric method applies. The only positive diagonals are exactly
the essential trace dofs (64/96/128 of them, all exactly `1.0`), which is
`DIAG_ONE` putting a unit row in.

The **Newton Jacobian is the same** on that path: free–free block symmetric to
5e-16, every Rayleigh quotient in `[−1.91, −1.30]`, for `∂F/∂ψ` of either sign.
Its whole-matrix asymmetry reads 0.52, but that is entirely `GetGradient`
zeroing the essential rows without eliminating the matching columns — a
boundary-condition artefact, not a property of the operator. **Measure the
free–free block, not the assembled matrix**, or this looks like a wildly
unsymmetric problem.

**On the extension path none of that holds.** Same measurement, same code, a
curved `Γ`:

| | fitted | extension |
|---|---|---|
| free–free rel `\|A_ij − A_ji\|` | 2e-16 | **5.4e-1** |
| Rayleigh quotients | 200/200 negative | 200/200 negative, largest −0.898 |

The cause is `HDGExtensionIntegrator`, and it is structural rather than a bug.
Its element matrix is

```
elmat( dof*di + i, dof*dj + j ) += w * nor(di) * shape(i) * L(j, dj)
```

— an outer product of the normal trace of basis `i` against the **path lifting**
of basis `j`, two unrelated vectors. Measured on `Γ_h` faces its relative
asymmetry is exactly **1.0**, and it is **16.8× larger** than the `(r q, v)`
domain term it sits beside, which is itself symmetric to the last bit
(`max|A_ij − A_ji| = 0`). `DarcyForm::AssembleFluxMassBdrFaces()` deposits it
straight into the hybridization's per-element flux block, so it reaches both the
local factorisation and, through the Schur complement, the global trace matrix.

**So the transfer technique costs self-adjointness.** The continuous
Grad–Shafranov operator is self-adjoint; `Δ*` with a transferred Dirichlet datum
is not, because the datum on a face depends on the flux along a path leaving the
element. That is a property of Cockburn–Solano transfer, not of this
implementation, and it means **meq's headline configuration is genuinely
non-symmetric**. UMFPACK's unsymmetric LU is the right solver for it.

What survives on both paths is that the **symmetric part is negative definite** —
every Rayleigh quotient measured is negative on either path. So `−A` has positive
definite symmetric part, which is what gives GMRES a convergence bound and what
makes a preconditioned Krylov method reasonable at all.

### A quarter of every Newton step was thrown away — fixed, and switched on

`UMFPackSolver::SetOperator` used to declare `void *Symbolic` as a **local
variable**, while `NewtonSolver::Mult` calls `prec->SetOperator( *grad )` every
iteration. The sparsity pattern does not change between Newton steps, so the
symbolic analysis was recomputed and discarded each time — and meq asks for
METIS ordering, which makes it dearer than the default.

**`Symbolic` is now a member and meq calls `SetReuseSymbolic()`.** Verified by
count rather than by clock — `theSymbolicAnalysisIsReusedAcrossNewtonSteps`
prints and asserts it:

```
  Newton took 4 iterations: 1 symbolic analyses, 4 numeric factorisations
```

One analysis, one factorisation per step. **A count, not a timing, on purpose**:
a timing here would be a measurement about this machine, and this is a
measurement about the code. It is also the *only* thing that could notice the
reuse lapsing — the pattern check is exact and re-analyses whenever it fails, so
a lapse costs speed and nothing else. No wrong answer, no failed convergence,
nothing a rate table or an error norm could ever see.

**Where it is switched on, and where deliberately not:**

| site | reuse | why |
|---|---|---|
| Newton path | **yes** | `NewtonSolver::Mult` re-`SetOperator`s the same object every iteration |
| Picard path | **yes**, and the solver is **hoisted to a member** | it was built inside `picardStep()`, so it had nothing to reuse across calls; Picard runs 122 to 290 full factorisations, a larger absolute win than Newton's |
| linear path | no | one factorisation, object destroyed immediately: retaining the analysis buys nothing and costs a copy of the pattern |

The Picard hoist is safe across `prepare()` rebuilding `reduced` every iteration,
because `SetReuseSymbolic()` documents comparing the **pattern**, entry by entry,
rather than the object — it "accepts a matrix reassembled in place and a matrix
rebuilt into a fresh object with the same structure alike", and re-analyses
whenever the check fails. Reuse is a request, never an assumption.

The measurement that motivated it:

| `n` | symbolic | numeric | backsolve | symbolic share |
|---|---|---|---|---|
| 12,544 | 20.8 ms | 72.1 ms | 10.9 ms | **22%** |
| 49,664 | 104.2 ms | 325.6 ms | 54.6 ms | **24%** |

### Threading, measured — and the two axes are not the same axis

**`tests/performance/` is the harness and is deliberately NOT a ctest**:
everything it reports is a timing, and this project's standing rule is that a
threaded timing on this machine is a measurement about the machine.
`TraceSolverScaling` returns non-zero only for the two *correctness* properties
that make its timings mean anything — threaded assembly reproducing serial
assembly bit for bit, and PARDISO reproducing UMFPACK. `scan.sh` drives it one
process per point, because MKL fixes its threading at first use and an
in-process sweep would measure the first setting five times.

**Two knobs, and sweeping them together is actively misleading:**

| | drives |
|---|---|
| `OMP_NUM_THREADS` | `DarcyHybridization`'s threaded element assembly |
| `MKL_NUM_THREADS` | UMFPACK's BLAS, PARDISO's internals, **and MFEM's element-local dense LU** |

That third entry is the surprise, and it is why the axes must be separated: swept
together, the `k = 3` assembly blow-up swamps every other column. Full account
under *Traps*; the short of it is `MKL_NUM_THREADS=1`, non-negotiable, and the
culprit is `ComputeH()`'s dense LU rather than UMFPACK.

**Threaded assembly is bit-exact and worth about 1.2x.**
`SetAssemblyMode( Threaded )` threads the element-local half of `ComputeH()`;
`SolverContract::threadedAssemblyReproducesSerialAssemblyExactly` requires
`0.000e+00` over two degrees and two meshes — not a tolerance, because the
mechanism is specific (element-local arithmetic reassociates nothing, the
scatter stays serial and in element order). Measured 1.15x–1.33x from four
threads up, and **0.86x at one thread**, the buffering costing more than the
serial loop it imitates. MFEM measures 2.09x for the hybridized assembly
*alone*; meq times all of `prepare()`, of which only that element loop threads.
The scatter is the ceiling and cannot be threaded — an unfinalized
`SparseMatrix` carries one `current_row` for the whole matrix, so two threads
writing provably disjoint rows still collide, and the failure is a hang.

**THE ASSEMBLY DEFAULT IS `Serial`, AND A GATE ON `omp_get_max_threads() > 1`
WAS TRIED AND REMOVED.** This is the place in this campaign where the isolated
benchmark actively misled. Under that gate `HighBetaConvergence` went from
**21.5 s to 39 s** — 1.8x slower, reproducibly — because MFEM forks a team and
buffers *per call*, so a caller that assembles once amortises it and one that
assembles hundreds of times inside a bordered Newton pays it every time. **Mesh
size does not separate the two cases**: HighBeta's meshes are 128 and 512
elements, and 512 is exactly where the isolated benchmark still showed a win.
The solver cannot know which caller it has. So `Threaded` is an informed
opt-in: **take it for a large mesh assembled a few times, leave it for a small
one assembled hundreds of times.**

**PARDISO beats UMFPACK on the trace solve, and beats it even sequentially** —
1.50x on analyse-plus-factor and 1.41x on the backsolve at 37,248 trace dofs
with `MKL_NUM_THREADS=1` on both sides, agreeing to 1.0e-14 or better at every
point. It scales to about 8 threads (1.87x setup, 1.96x solve; 16 buys nothing
more), which would make the end-to-end gap 2.83x and 3.13x. **meq cannot have
that**, because `MKL_NUM_THREADS` is process-wide and the setting that makes
PARDISO fast is the setting that makes `ComputeH()` forty times slower at
`k = 3`. It stays unreachable until the element-local factorisation stops going
through threaded MKL — `LocalFactorMode::Batched`, or an
`mkl_set_num_threads_local()` around the trace solve. Neither is done.

**cuDSS is correct here and not measurable here.** It agrees with UMFPACK to
**3.5e-14 or better** from 9,408 to 148,224 trace dofs, which is the question
that needed answering locally. The timings are irreproducible **by a factor of
thirty at fixed size** — three runs of the same binary on the same 9,408-dof
problem gave 1.05 s, 1.13 s and 2.81 s — because WSL2 shares the GPU with the
Windows host (`nvidia-smi` reports 7.5 of 8 GB used with no compute processes)
and because a consumer card runs FP64 at 1/32 of FP32 where a datacentre part
runs it at about 1/2. It won nothing on this hardware. **Do not re-time it
here**; it needs a different machine, not another afternoon.

**AND ONE TRAP THAT WOULD HAVE POISONED ALL OF IT.** cuDSS queues work on a
stream and returns. Timed without a device synchronise, the warm setup of a
**148,224-unknown factorisation reads 2.0e-04 s** — four orders of magnitude
out, and *plausible enough to publish*, because "the reordering was reused, so
of course it is fast" is a story that fits. `TraceSolverScaling` calls
`cudaDeviceSynchronize()` inside the timing loop for **every** solver, CPU ones
included, so the sync can never be the thing that was forgotten. **Any future
device timing in this project must do the same.**

**What a fresh `GradShafranovSolver` does:** `setTraceSolver()` is `UMFPack`
(the only backend present in every build, and what every rate in the suite was
measured with); `setAssemblyMode()` is `Serial`, unconditionally, per the gate
above; `MKL_NUM_THREADS` is 1, set on every ctest.

**The fastest reachable set is PARDISO plus threaded assembly at 8**, and it is
worth **1.24x** end to end on one linear solve at `k = 2, n = 64` — 0.590 s
against the default's 0.732 s. It is not the default because
`MFEM_USE_MKL_PARDISO` is off in most builds and oneMKL's terms are not
everybody's to accept; `setTraceSolver()` and `traceSolverAvailable()` are how
a caller who has it takes it. **Read that number for its size**: the whole
ladder is 1.24x, dwarfed by the `MKL_NUM_THREADS=1` fix that preceded it (3x on
ordinary tests, 140x at `k ≥ 3`). The large win here has already been taken.

### What to do, in order of value

**Reordered 2026-08-30 by the threading measurements, which put a new item at
the top and demoted the one that used to be there.**

0. **Get `ComputeH()`'s element-local dense LU off threaded MKL.** The largest
   single item, and it was invisible until the link line was fixed. Forty times
   at `k = 3`, and it is *also* what makes PARDISO's 2.8x unreachable, since
   `MKL_NUM_THREADS` is process-wide and meq must keep it at 1. Two routes and
   neither is done: `mkl_set_num_threads_local()` around the trace solve, or
   `LocalFactorMode::Batched` for the local factorisations. Anything else on
   this list is smaller.

**And the older list below, which the asymmetry finding had already reordered.**
An earlier version ranked Cholesky second and a per-cell Cholesky fourth, on
measurements taken only on the fitted path. Both are wrong wherever `Γ` is
curved.

1. ~~**Reuse the symbolic factorisation across Newton steps**~~ — **DONE**, and
   `theSymbolicAnalysisIsReusedAcrossNewtonSteps` asserts it by count. About 23%
   off each step for no numerical change, valid on both paths.
2. **Precondition the Newton path's GMRES fallback.** It has none, while the
   linear path's fallback has a `GSSmoother`. That inconsistency is a plain bug
   and is worth fixing whatever else happens; without SuiteSparse the Newton
   path is currently running unpreconditioned GMRES on every step.
3. **Symmetric methods — fitted path only.** CG or MINRES on `−A`, and a CHOLMOD
   Cholesky, are all correct on a fitted mesh and all invalid on the extension
   path. CHOLMOD is already linked (`-lcholmod` via SuiteSparse) but MFEM wraps
   only UMFPACK and KLU. Worth having only if fitted-domain runs become a
   workload in their own right; the driver's target is curved boundaries, where
   none of it applies.
4. **~~Cholesky on the per-cell flux block~~ — do not.** `(r q, v)` is SPD, so
   this looked like a free 2×. It is wrong on the extension path for exactly the
   reason above: `AssembleFluxMassBdrFaces()` puts a term with relative
   asymmetry 1.0, sixteen times larger than the mass term, into the very block
   `LU_A` factors. A correct version would have to branch on whether an
   extension is configured, for at best a small gain on 9×9 to 30×30 blocks.
   Not worth a conditional in that inner loop. **`LUFactors` is the right
   choice and should stay.**
5. **Do not reach for AMG.** At 2D serial sizes direct wins: 50k trace dofs
   factorise in 0.4 s. AMG earns its place in 3D or in parallel, and serial MFEM
   has none anyway — `HypreBoomerAMG` needs MPI.
6. ~~**Make the trace solver selectable at run time**~~ — **DONE**,
   `setTraceSolver()`. All three backends agree to 5e-14 through the whole
   solver. **The default stays UMFPack**, because PARDISO needs a licence most
   builds do not have; a caller who has it gets 1.24x for one line, without
   waiting on item 0 — PARDISO's *sequential* advantage is 1.50x on the
   factorisation and needs no MKL threads at all. Item 0 is what unlocks the
   further 1.9x on top.
7. **Do not re-time cuDSS here.** Correctness is established, the timings are
   irreproducible on this machine by a factor of thirty, and the reason is the
   machine. It needs a datacentre part, not another afternoon.

### PARDISO and the MKL link line: what is still true

**Two sections stood here and are deleted**, both superseded by *Threading,
measured* which re-measured every number in them. They were taken while a stray
`libmkl_sequential` on the link line made MKL sequential, so they were
sequential-against-sequential timings that were never labelled as such.
Recoverable from git (`cbb210b`). Three things survive.

**oneAPI MKL 2026.1 fixed a real defect, and it was packaging.** Debian's
`intel-mkl` 2020.4.304 returned error `-3` from PARDISO at 12,544 and 28,032
trace dofs — the sizes meq actually runs. Against oneAPI at
`/opt/intel/oneapi/mkl/latest` it runs at every size and agrees with UMFPACK to
round-off. **Both paths must keep shipping**: oneMKL's licence is not
everybody's to accept and `MFEM_USE_MKL_PARDISO` is off in most builds, so
`tests/convergence/TraceSolverComparison.cpp` compiles and passes either way,
skipping the PARDISO columns when they are absent. It asserts only the
**agreement**; the timings are printed, per the standing rule.

**`SetReuseSymbolic()` exists on `PardisoSolver` too**, with the same
pattern-comparison contract as `UMFPackSolver`'s, so *A quarter of every Newton
step* carries over without re-arguing.

**~~THE LINK LINE STRADDLES TWO MKL INSTALLATIONS AND TWO THREADING LAYERS~~ —
FIXED 2026-09-01, AND THE LINK LINE WAS ONLY HALF OF IT.** It used to carry
oneAPI's `mkl_gnu_thread` from the explicit `BLAS_LIBRARIES`, and Debian's
`mkl_intel_thread` plus `iomp5` behind SuiteSparse — two MKL *versions*, two MKL
threading layers, and two OpenMP runtimes, working only because oneAPI came
first in the link order.

**THE HALF THAT NO CMAKE VARIABLE COULD FIX**: Debian's `libumfpack.so` and
`libcholmod.so` carry a hard `NEEDED` on `libblas.so.3`, which on this machine
is a Debian alternatives symlink to **`libmkl_rt.so`**. So Debian's MKL 2020
loaded at runtime whatever the link line said. **Editing the link line alone
would have looked like a fix and left the real one in place.**

**The fix is meq's own SuiteSparse**, at `../suitesparse`, exactly as
`../sundials/cuda-install` is meq's own SUNDIALS. It is built from **v7.12.2 —
the same version Debian ships** — so the only variable that changes is the BLAS,
not the numerics:

```sh
MKL="-L/opt/intel/oneapi/mkl/latest/lib;-Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib;\
-lmkl_intel_lp64;-lmkl_gnu_thread;-lmkl_core;-lgomp;-lpthread;-lm;-ldl"
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/home/ian/projects/suitesparse/install \
  -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;btf;camd;ccolamd;colamd;cholmod;klu;umfpack" \
  -DBLAS_LIBRARIES="$MKL" -DLAPACK_LIBRARIES="$MKL" -DBLA_VENDOR=Intel10_64lp \
  -DBLAS_INCLUDE_DIRS=/opt/intel/oneapi/mkl/latest/include \
  -DCMAKE_C_FLAGS=-I/opt/intel/oneapi/mkl/latest/include \
  -DBUILD_STATIC_LIBS=OFF -DSUITESPARSE_USE_CUDA=OFF -DSUITESPARSE_USE_FORTRAN=OFF
```

**Three snags, none of them obvious, all cost a configure round trip.**

* SuiteSparse takes `BLAS_LIBRARIES` as-is if you set it — but its threading
  probe then does `string(REGEX MATCH "^Intel" ... ${BLA_VENDOR})` on a
  `BLA_VENDOR` that is **empty**, because supplying the libraries skips
  `find_package(BLAS)`. Set `BLA_VENDOR` anyway; it is informational.
* That probe's `try_run` passes `LINK_LIBRARIES` but **never
  `BLAS_INCLUDE_DIRS`**, so it fails on a missing `mkl.h`. The include path has
  to go in `CMAKE_C_FLAGS` as well.
* Only the nine projects MFEM needs are built. `SUITESPARSE_ENABLE_PROJECTS`
  keeps GraphBLAS — which dwarfs everything else — out of the build entirely.

**AND ON THE MFEM SIDE THE STALE CACHE HAD TO GO.** `mfem_find_package` resolves
SuiteSparse's `BLAS` requirement through **CMake's own `FindBLAS`**, which had
cached Debian's paths in `BLAS_mkl_*_LIBRARY`. Deleting those
(`cmake -U "SuiteSparse_*" -U "BLAS_mkl_*" -DSuiteSparse_DIR=...`) is what lets
MFEM's explicit oneAPI `BLAS_LIBRARIES` win.

**Why `FindBLAS` chose the WRONG threading layer, which is worth knowing because
it will do it again**: it selects `mkl_gnu_thread` + `gomp` only when a
**Fortran** compiler is loaded and is GNU; otherwise it takes
`mkl_intel_thread` + `iomp5`. SuiteSparse and MFEM are C/C++ builds, so the
default is Intel threading — mismatched with everything else here. That is
`FindBLAS.cmake:505-513`, and it is why the libraries are handed over
explicitly rather than discovered.

**VERIFIED THREE WAYS, because "one MKL" is a claim about the loader and not
about a text file.** The installed `config.mk` carries oneAPI and nothing else;
`ldd` on a test binary resolves `libmkl_intel_lp64`, `libmkl_gnu_thread` and
`libmkl_core` to `/opt/intel/oneapi/...` with **no `libmkl_rt`, no `libblas.so.3`,
no `libiomp5`, and `libgomp` as the only OpenMP runtime**; and `LD_DEBUG=libs`
shows the only MKL objects *initialised* are oneAPI's. Debian's `libmetis` is
still on the line and is **not even loaded** — CHOLMOD now bundles its own,
prefixed `SuiteSparse_*`/`cholmod_*`, so there is no symbol collision and
UMFPACK keeps its METIS ordering. Suite green twice — **23/23 as it then was**, 490 s and 540 s on
the same code, which is the run-to-run spread on this machine and is worth
knowing before anyone reads a 10% change in a suite time as a result.

**`MKL_THREADING_LAYER=GNU` is inert and has been dropped everywhere.** It only
ever configured the `libmkl_rt` **dispatcher**, and there is no dispatcher any
more: `SolovievConvergence`, `NewtonConvergence` and `FieldConvergence` produce
**bit-identical output** with it set and with it unset. It was briefly kept as a
free guard, and that argument was rejected for a reason that generalises —
**most machines will not have MKL at all**, so it named a library the majority
of readers do not have, about a failure they cannot suffer. A guard that is
usually inapplicable is a false instruction to whoever reads it next, and it
spends the credibility of the ones that are real. The driver's matching runtime
warning went for the same reason. **Choosing a threading layer is CMake's job**,
which can look at what the build actually links.

**The hazard went off once, and in the opposite direction to the one predicted.**
The prediction was wrong answers from mixed threading layers. What actually
happened is that a stray `libmkl_sequential` was **load bearing**: it made MKL
resolve sequential for everything, which is why every timing in this file before
2026-08-30 was fast and why nobody noticed threaded MKL is ruinous here.
Removing it — the *correct* thing to do — exposed a 140x regression.
**The suite was never deliberately sequential; it was accidentally so.** Same
species of finding as the threaded-BLAS one: a property of the link line
masquerading as a property of the code.

**And the honest caveat on all of it**: on a hard case the dominant cost is not
the global solve. Globalisation is a bigger lever than anything in this
section, and under `CondenseThenLinearise` the element-local *non-linear*
iteration dominates outright.

**The transferable lesson**, which is the same one the Solov'ev coefficients
taught: a property measured on the easy configuration is not a property of the
code. Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on the
geometry meq is actually for.

## Traps

**Every run needs `MKL_NUM_THREADS=1`**, which every registered ctest sets. It
is a no-op wherever MKL is absent, so it costs a non-MKL machine nothing.

**`MKL_THREADING_LAYER=GNU` was the other half of this trap and is now gone**,
from the ctests, the performance harness and the driver alike. It guarded the
`libmkl_rt` **dispatcher** — what `/usr/lib/x86_64-linux-gnu/libblas.so.3`
resolves to on this machine — whose default threading layer silently corrupted
UMFPACK's BLAS-3: you got numbers, and they were wrong. meq no longer loads it,
building against its own SuiteSparse; see *PARDISO and the MKL link line* for
why the variable was dropped rather than kept. **The trap is still live for
anyone using Debian's SuiteSparse**, and `../mfem-hdg-dev/CLAUDE.md` records it.

`MKL_NUM_THREADS=1` is the **factor of 140** one. Measured 2026-08-30, one
process, nothing else running: `SolovievConvergence` **1.24 s** against
**177.83 s** at `=16`, and `FieldConvergence` 1.11 s against 198.35 s — at 273%
of *sixteen* cores, which is the shape of the answer: the threads are spinning
on barriers, not working. `SamplerConvergence`, which never reaches a direct
solver, does not move.

**THE CULPRIT IS NOT UMFPACK, AND THE FIRST VERSION OF THIS SECTION SAID IT
WAS.** The plausible story — UMFPACK calling threaded BLAS on the small dense
frontal matrices of a sparse LU, thousands of times, each call costing a fork
and a barrier — fits every whole-test fact and is **wrong**. Separated out,
`UMFPackSolver::SetOperator` degrades by about **40%** across the whole thread
range, never more.

**It is MFEM's element-local dense LU.** `ComputeH()` factors `A`, forms the
Schur complement and factors that, on blocks of order 10–30, once per element,
through `LUFactors` → LAPACK → MKL:

| assembly + reduction | `MKL=1` | `MKL=2` |
|---|---|---|
| `k = 2, n = 64` | 0.469 s | 0.440 s — untouched |
| `k = 3, n = 16` | 0.056 s | **2.306 s** |
| `k = 3, n = 32` | 0.216 s | **8.973 s** |

**`k = 2` does not move and `k = 3` degrades by a factor of forty**, putting
MKL's threading threshold between those two block sizes. The suite runs `k` up
to 4 and every element of every mesh pays it, which is why the whole-test
factors dwarf anything in the solver columns. **A production run at `k ≥ 3` with
MKL threading left on would be unusable and would look like a solver problem.**

**It is MKL and not OpenMP**, which the obvious cross settles: `FieldConvergence`
at `OMP=16 MKL=1` is **1.19 s at 99% CPU**, and at `OMP=1 MKL=16` is 183.54 s at
266%.

**Keep the wrong answer in view, because of HOW it was wrong.** It fitted every
fact then available — the collapse tracked exactly the tests that reach a direct
solver — and **no whole-test timing could have refuted it**, because in those
tests UMFPACK and `ComputeH()` always appear together. Only separating the
columns did. Same lesson as the Solov'ev coefficients: checking a story against
the evidence that suggested it proves nothing.

**Why it never bit before**: the link line used to carry `libmkl_sequential`
behind SuiteSparse, so MKL resolved sequential whatever it was asked for. See
*PARDISO and the MKL link line*. **`MKL_NUM_THREADS=1` is right for meq's
solver and wrong for PARDISO**, which is a genuine tension rather than an
oversight and is what item 0 of *What to do* is about.

**This machine's GPU is for development, not for performance conclusions.**
There is an RTX 2070 SUPER with CUDA 13.3, which is enough to write and debug a
device path and to check that it gives the same answers. It is **not** a
representative part for meq: MFEM is built `MFEM_USE_DOUBLE`, and consumer
NVIDIA cards run FP64 at 1/32 to 1/64 of their FP32 rate where datacentre parts
run it at about 1/2. A GPU timing taken here can invert the conclusion a
production part would give, so treat local device runs as correctness evidence
and take timings elsewhere. `TODO` carries the detail. This is the same species
of warning as the threaded-MKL note under *On SUNDIALS* — a measurement on this
machine that is not a measurement about the code.

**Never run a bare `make -j`, anywhere.** With no argument it is unbounded, and
this machine is WSL2 — the job count goes to the host's core count with a fraction
of the host's memory behind it, and the whole VM falls over rather than the build
merely failing. Always give a number: **4 to 8**, and for the MFEM tree
specifically **`make -j4`, never more**, since its translation units are large
enough that even 8 exhausts memory here. `cmake --build ... -j4` likewise.

**And `make clean` after editing any MFEM header.** MFEM's
makefiles have no `.d` files and no header dependency tracking. That trap has
produced heap corruption in unrelated functions and "unimplemented" aborts for
methods that had just been added. meq's own CMake build tracks headers properly;
this applies only to the MFEM tree.

**A convergence target scaled to `‖r₀‖` makes a good initial guess FAIL, and
the better the guess the more certain the failure.** MFEM's `NewtonSolver` stops
at `‖r‖ ≤ max(rel_tol·‖r₀‖, abs_tol)` with `‖r₀‖` measured at the iterate it was
handed. Warm-start from a converged answer and `‖r₀‖` is already small, so the
target shrinks with it — and past a point it falls below the round-off floor,
where nothing can meet it. Measured on Example 5 at `k = 3, n = 8`, restarting
from the exact answer: the residual reaches the floor in **two** iterations and
is then reported as a failure at the thirtieth.

```
   it 0   1.180694e-03
   it 1   8.782203e-13
   it 2   4.151551e-14      <- converged; the floor
   ...    ~4e-14 for 28 more iterations, then FAIL
```

The target was `max(1e-12 × 1.18e-3, 1e-14) = 1e-14`, under the 3.7e-14 this
problem can reach. Cold, `‖r₀‖ = 11.2` gives a target of 1.1e-11 and it converges
in four. **This is the failure mode that matters most for how meq will be used** —
moving to an adjacent equilibrium should cost one or two Newton steps, and it
instead threw away a converged answer.

`solve()` therefore takes the reference from the **cold** iterate — the Dirichlet
datum alone, where this solve would have started with no guess — and sets a pure
absolute target from it. `rel_tol` keeps exactly the meaning it always had, a
cold solve is bit-identical because there the reference *is* `‖r₀‖`, and only a
warm one changes: from failing to converging in one step. It costs one extra
residual evaluation, and only when a guess was set.

**`NewtonSolver` needs `iterative_mode = true`, and the failure is silent.** The
Dirichlet values ride in the iterate, so with `iterative_mode` left false
`NewtonSolver` zeroes `x` on entry and throws the boundary data away — silently,
because the residual is masked on exactly those rows, so nothing complains and
the answer is merely wrong.

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
dead in **every** meq build and **a meq debug build revives none of them**. What
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

**Every GS-2 §4.2–4.5 source vanishes at `ψ = 0`, so the paper's own problem has
a trivial branch — and meq falls into it.** `F(r, 0) = 0` for eqs (24), (25),
(26) and (27); for (25) because `b = 2` makes the bracket `O(ψ²)`. With
homogeneous Dirichlet data `ψ ≡ 0` therefore *solves* the problem, and Newton —
which starts from the Dirichlet data — lands on it and stops in **zero
iterations** with an identically zero residual. Not a hypothesis: GS-1's
Algorithm 2 literally opens `ψ⁰ ; // Non-trivial initial guess`.

**`setInitialGuess()` exists** — in `Coefficient` and `GridFunction` overloads —
and it is what those four cases need; without it `prepare()` resets the iterate
and `solve()` calls `prepare()`, so no caller could work round the trivial
branch. They are posed with a non-homogeneous ramp that puts `ψ = 0` in the
interior. **Under NPC the guess seeds the potential and the trace and leaves the
flux block at zero**, which is a real difference; see *A warm start no longer
shows up in `‖r₀‖`*.

**`DarcyHybridization` freezes the element-local Newton's initial guess at
`FormLinearSystem()` time, and `ComputeSolution()` will not take another**, so
seeding the vector you hand `RecoverFEMSolution()` is inert. On an ordinary
Newton path that costs iterations; on anything that **differentiates the
residual by differencing it** it costs correctness — measured, it moved
`max ψ_h` from 0.896 to 3.84 for a perturbation of 9e-6. `formSystem()` is the
fix. Full account under *The trap that cost the most*.

**~~`DarcyForm::Reconstruct()` returns a different function where `∂F/∂ψ`
vanishes~~ — FIXED.** Silently, per element, so a profile with a flat segment
corrupted part of the domain and left the rest exact, and a whole-domain check
missed it. Full account and the transferable lesson under *Post-processing is
back*; the part that matters here is that **the fix is on `gf-hdg-dev` and on no
other branch**, so a `meq-integration` rebuilt without it silently loses this.

**`mfem::Mesh::GetElementTransformation( int )` HANDS OUT SHARED SCRATCH, AND
THREADING IT IS A SILENT WRONG ANSWER.** MFEM's own comment: *"The returned
object is owned by the class and is shared, i.e., calling this function resets
pointers obtained from previous calls."* It is one `IsoparametricTransformation`
member of the `Mesh`. Two threads evaluating in different elements overwrite
each other's transformation — no crash, no error, a point transformed by the
wrong element. The reentrant route is the `( i, IsoparametricTransformation * )`
overload into a thread-local.

**~~Three call sites in meq use the shared overload~~ — SIX, AND THE UNDERCOUNT
IS THE POINT.** Five in `src/meq/FluxSurfaces.cpp` — `setBandExtension`'s face
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

**THE ONE THAT ACTUALLY BLOCKS THREADING IS `Mesh::FindPoints`, AND IT CANNOT BE
FIXED LOCALLY.** It loops over every element through that same shared
transformation *and* builds a vertex-to-element table on the way, so it is not
reentrant — an attempt to share one `ContourTracer` across threads aborts with
*"the axis is not in the mesh"*. Every entry point can reach it, and
`traceFromAxis()` takes it **once per surface unconditionally**, because it
samples the axis with no element hint. So "the tracer reports zero fallbacks" —
which `INVERSION-PLAN.md` §11.3 offered as the reason honouring this rule is
free — **is true of `trace()` and false of the entry points**.

**`mfem::Mesh::FindPoints` is `O(elements × points)`** — a brute-force scan over
element centres. It caps sample-cloud sizes in any off-grid error measure.

**meq will be an early user of a thinly-tested MFEM combination.** Fixed-boundary
Grad–Shafranov is a Dirichlet problem, so the trace carries an essential BC; and
Newton makes it nonlinear. That pairing —
`DarcyHybridization::SetEssentialBC` together with a nonlinear reduced operator —
was broken on the MFEM branch until recently: `EliminateTraceTrueDofsInRHS`
returned early for nonlinear problems and the essential BC was silently ignored,
so the solver was solving a different problem than asked. **It is fixed**
(`fem/darcy/darcyhybridization.cpp`, and the reasoning is commented there), but
`../mfem-hdg-dev/CLAUDE.md` records that **no regression covers the combination**.
If stage 4 produces a converged-but-wrong answer near the boundary, look here
first rather than at meq's assembly.

**toml11's `find_or<double>` silently returns the default when the node is an
integer.** Paid for once. The original bug was that everything was read through
`.as_floating()`, so `RMin = 0` threw instead of converting — and the obvious fix
is worse than the bug. In toml11 4.4.0, `toml::find<double>` on an integer node
throws just the same, and `toml::find_or<double>(v, "RMin", 1.0)` **returns
1.0**, because a failed conversion is indistinguishable from a missing key. That
turns a loud failure into a silent wrong answer in the configuration, which is
the last place you want one. `Config.cpp` therefore reads numbers through its own
`asFloat()`, accepting `is_floating()` or `is_integer()` explicitly.

**Four tooling traps that cost one session six idle polling loops**, and the
reason they are worth writing down is that a polling loop which can never
terminate is *invisible* when the thing it polls for is also reported some other
way.

* **`cd` persists between Bash tool calls**, so a `cd x && ...` that assumes the
  repo root silently runs somewhere else. Use absolute paths.
* **Shell VARIABLES do not persist**, only the working directory does — the same
  trap read the other way round, and the more dangerous half. `L=...` in one call
  is gone by the next, so `until grep -q '^EXIT=' "$LOG"` with an unset `LOG`
  becomes `grep -q '^EXIT=' ""`, and **`grep -q` with no file argument reads
  stdin**. In a detached background process stdin never reaches EOF, so the loop
  blocks silently for the life of the session.
* **A completion marker appended to a log is not necessarily at the start of a
  line.** `echo "EXIT=$?" >> log` after a Boost test binary lands on the end of
  its final ANSI reset sequence, because Boost's coloured output does not
  terminate with a newline — so `grep -q '^EXIT='` never matches. ctest's output
  *does* end with a newline, so the same waiter works on `ctest` and hangs on a
  bare test binary, which gets diagnosed as "the run is slow". Drop the `^`.
* **`grep -c` exits 1 on zero matches**, so a background check reports failure
  spuriously.
* **`pgrep -f` AND `pkill -f` MATCH THE INVOKING SHELL'S OWN COMMAND LINE.**
  This is the one that has now fired three times in a day, twice on the same
  afternoon, and it fails in two different directions. A waiter built on
  `pgrep -f "cmake --install build"` matches *itself*, so it reports the install
  running for as long as the session lasts and the thing it waits for is never
  seen — the install had in fact finished minutes earlier. And `pkill -f ctest`
  in a shell whose own command line contains `ctest` **kills that shell**, which
  surfaces as an unexplained exit code 144 (128 + SIGTERM) rather than as
  anything mentioning the pattern. Use `pgrep -x`/`pkill -x`, which match the
  process *name*, or check for the artefact — a file's timestamp, an exit
  marker — rather than for a process at all.

**If a waiter is used, check it actually fired.**

**AND NEVER `cmake --build` WHILE `ctest` IS RUNNING.** Done here on 2026-09-02
for a comment-only header change, which relinked ten test binaries *including
the one ctest was starting*. Nothing errored — replacing a running executable's
file leaves the running process on its old inode — so the run continued and
would have reported a pass over a mix of old and new binaries. A suite result is
only an acceptance measurement if every binary in it is the one being committed;
that run was killed and re-run rather than believed.

**The four defects meq reported to MFEM are closed**, though
`HDG-DEFECTS-FROM-MEQ.md` itself is **not** gone — it is alive on `gf-hdg-dev`
and `gf-hdg-subdomains-dev` and deleted only on the symbolic-reuse line, which
is the modify/delete conflict in the merge recipe. An earlier version of this
sentence said it was gone, from a working-tree listing taken while that tree was
on the branch that deletes it. Checked 2026-08-29
one at a time, because "closed" arrived in four different ways: `Reconstruct()`
on a singular local matrix was **fixed** (meq's test flipped red to green on
it); `φ_h` unreachable after a solve was **fixed**, as
`mfem::TransferredDatumCoefficient` in `extension_hdg.hpp` — meq still calls
`setTransferredBoundary()`, because rebuilding `η₅` on it is meq's work and is
not done; `ReconstructFluxAndPot()` lifting only domain integrators was
**withdrawn as not a defect**, and meq had measured it harmless, which was the
right answer for the wrong reason; and `ComputeHDGFaceEnergy()` ignoring an
installed `HDGStabilization` is **not re-measured** — a code read today still
shows it computing the `{h⁻¹Q}` form with no call into an installed hook, but
meq uses `meq::ResidualEstimator` rather than `HDGErrorEstimator` so it costs
meq nothing either way. **That last row is a code read and not a measurement.**
Anyone who needs `Energy` mode should measure before trusting it.

## Testing stance

**A test asserts the behaviour that is wanted, and fails until it is there.**
Never the reverse. A test that asserts a known defect passes while the defect
stands, and that makes a green suite compatible with a broken solver — which
destroys the only property the suite is for: **100% green must mean there is no
lurking defect.** So a defect gets a *failing* test naming it, not a passing test
recording it.

**This was got wrong, corrected, and then vindicated.**
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` was first written
to assert the *corruption*, and passed. Rewritten to assert that `ψ*` is a
post-processing of `ψ_h` on every element, it went **red** — and stayed red
until MFEM fixed the reconstruction, at which point it flipped green on its own
and became the regression. **That flip is what said to put the driver back on
`ψ*`.** A test asserting the defect would have said nothing.

Two consequences. The suite is expected to be **red while a known defect
stands**, and that is the intended signal rather than a broken build — the
failing test's message is the record. And a failing test must say what fixes it,
because a red suite nobody can action is one people learn to ignore.

**A convergence table finds what unit tests do not.** This is the lesson from the
sibling MFEM branch, where reproducing a published table turned up three library
defects a passing unit suite had walked past for years — two of them silent.

The specific hazard here is worth stating plainly, because it defeats the obvious
test design: **a wrong sign convention converges, at the right rate, to the wrong
function.** An order-of-accuracy study cannot catch it. Only a comparison against
a closed form can. That is why stage 3 is pinned to
`refs/HDG-GradShafranov-Adaptive.pdf` §4.1 — an exact NSTX Solov'ev solution with
`A = −0.52` whose twelve coefficients are published to 15 digits — and not to a
self-convergence study.

So the test ladder is, in order:

1. **Unit**, in `tests/unit/` — Boost.Test. Config parsing, spline value *and*
   derivative against closed forms, `dFdPsi` against a finite difference.
   And, in `NewtonConvergence.cpp`, the *assembled* Jacobian against a central
   difference of the assembled residual — which is a different and stronger
   check than the one on `dFdPsi` alone.
2. **An exact solution**, as an absolute-error regression.
3. **Convergence rates**, in `tests/convergence/` — asserting a *rate*, not a
   tolerance: `k+1` for `ψ` and `∇ψ`, `k+2` for `ψ*` at `k ≥ 1`, over `k = 1…4`
   and several dyadic refinements, against Tables 1–5 of the first paper and
   Table 1 of the second.

Expect rates to stop improving beyond `k ≈ 5` or `h/8`: the papers report
round-off dominating there, so a table that flattens at that point is behaving
correctly and is not a bug to chase.

**Unfitted convergence needs a two-tier rate assertion.** On a fitted mesh a
single-pair rate is tight enough to assert at 0.15 of slack. On the extension
path it is not: `Ω^h` is the union of background elements lying inside `Ω`, and
*which* elements those are is not a smooth function of `h`. At `k = 3` the error
is not even monotone in the mesh count — 2.02e-5, 3.03e-5, 3.88e-5 at
`n = 6, 7, 8` — and `{12,24,48,96}` is a worse sequence than `{8,16,32,64}`.
That is geometry, not the transfer: all four path families give the same
numbers. So `ExtensionConvergence` allows 0.30 per pair and asserts 0.15 on the
rate across the whole sequence.

**Newton is not monotone on a stiff source, so assert on the best triple, not on
every one.** Example 5 gives a clean four-step quadratic run; the GS-2 sources do
not. A typical history wanders for four to eight steps with observed orders of
0.02, 0.33, −0.37, −5.9 and *then* enters the quadratic regime and finishes
1e-9 → 1e-13 in two. §4.5 at `k = 2` even rises mid-run, 1.7e-2 → 5.3e-1 at
iteration 3, before recovering. Print the whole history; assert on the best
triple above the round-off floor.

**A self-convergence study on a rectangle cannot demonstrate `k+1`, and the cap
has nothing to do with the physics.** Measured on a *Solov'ev* source — constant
in `ψ`, `dF/dψ ≡ 0`, one Newton step — with homogeneous data on the benchmark
box, the self-difference rate is 2.00/2.86/2.96 in `ψ` and 1.88/2.25/2.18 in `q`:
flat at about 3 and 2.2 from `k = 2` on, with nothing nonlinear anywhere. That is
the `r² log r` corner term of a right angle — the solution sits in `H^{3−ε}` and
its gradient in `H^{2−ε}`, and no polynomial degree recovers it. The
exact-solution studies are immune because there the datum *is* the trace of a
smooth solution.

It is worse on a polygon approximating a curved boundary. GS-1 Example 6 on a
fixed 40-gon gives 2.12 in `ψ` at `k = 3`, against 1.995/3.000/4.000 for a
Solov'ev control **on exactly the same meshes** — the singular exponent at an
interior angle of `π − 2π/40` is `40/38 = 1.0526`, and interior pollution goes at
about twice that. **This is precisely the difficulty GS-2's curved-boundary
technique exists to remove**, and it is worth knowing before anyone designs
another fitted-polygon study.

**Mutation-test a suite you are relying on.** `ProfilesTests` and `SourceTests`
were checked by deliberately introducing fifteen defects — a dropped `r²`, `p'`
replaced by `p`, `μ₀` on the wrong term, the Solov'ev sign flipped, an
off-by-one in the interval lookup at a knot, `write()` truncated to six
significant figures — and confirming each was caught. That is a cheap way to
find out whether a green suite means anything, and worth repeating for the
convergence tests, where the risk of a test that passes regardless is higher.

Analytic solutions live in `tests/analytic/`, and they form a deliberate ladder
in how the source depends on `ψ` — which is to say, in what they demand of a
Newton Jacobian:

| fixture | source | `∂F/∂ψ` | what it catches |
|---|---|---|---|
| `Soloviev.hpp` | constant in `ψ` | `0` | the discretisation alone; Newton must converge in one step |
| `McCarthy.hpp` | linear in `ψ` | `T`, a nonzero constant | a Jacobian missing its mass term |
| `ManufacturedNonlinear.hpp` | `ψ²` and `e^(−ψ)` | varies | a Jacobian that is present but wrong |
| `SimilarityExponential.hpp` | `f₀e^{nψ}(1+εr²)` | `nF` | the same, against an *exact* rather than manufactured solution |

The last two are both nonlinear but are not redundant, and the difference is
worth keeping straight. **`ManufacturedNonlinear` is manufactured**: a convenient
`ψ` was chosen and `F` built to fit it, so its `F` is not of a form any physical
profile produces. **`SimilarityExponential` is exact**: the free function is
chosen — `f(u) = f₀e^{nu}` — and the solution follows from a similarity
reduction (Kaltsas & Throumoulopoulos; `refs/Refs.md`). So it tests the solver
against a nonlinear equation somebody might pose, rather than one reverse
engineered from an answer. It is also a cleaner shape: the whole `ψ`-dependence
is one exponential, so `∂F/∂ψ = nF` exactly, and at `n = 3` that exponential
varies by a factor of 29 across the benchmark rectangle.

The middle rung earns its place: with `∂F/∂ψ = 0` a Solov'ev run cannot tell a
correct Jacobian from an absent one, and in the nonlinear case the algebra is
messy enough to hide a factor. McCarthy's `∂F/∂ψ` is a single constant, so the
mass term is either there or it is not.

Each fixture checks its own transcription rather than trusting it: `deltaStarFD()`
recomputes `Δ*ψ` by central differences and the suite asserts it against `−f()`.
For McCarthy's eighteen-term Bessel and Neumann expansion that agrees to 3.6e-6,
and its analytic gradients match finite differences to 5e-9.

## Layout

```
src/meq/     the library. Config, Profiles, Source, SourceFactory,
             GradShafranov, BoundaryShape, Estimator, Field, Sampler,
             WarmStart, Output -- all ported and all under the naming check.
             CriticalPoints, FluxSurfaces, SurfaceAverage, SurfaceFit and
             Zernike are the solution-inversion work and were written here
             rather than ported; see that section.
apps/        drivers. Only meq.cpp, and MEQ_BUILD_APP defaults ON.
tests/       unit/ (Boost.Test), convergence/ (rate assertions),
             analytic/ (closed-form solutions used by both),
             performance/ (TraceSolverScaling + scan.sh -- built, NOT a ctest,
             because every number in it is a timing)
tools/       plotting and visualisation. plot_equilibrium.py reads the
             NetCDF; tools/README.md says which of the three output formats
             goes with which reader, and why they are not interchangeable
examples/    TOML run configurations
refs/        Refs.md is tracked; the PDFs are gitignored, fetch by doi
attic/       free-boundary/ -- not ported, not built, kept visible; its own
             README says why
docs/        the Sphinx manual, published to Read the Docs. Built by
             `make -C docs html` or the `docs` CMake target, both under
             -W to match .readthedocs.yaml's fail_on_warning. Citations
             are sphinxcontrib-bibtex over docs/references.bib, which is
             written by hand from refs/Refs.md and is meant to agree with
             it. docs/manual/ is the pre-Sphinx LaTeX manual, moved intact
             and keeping its own Makefile.

             THE SPLIT IS DELIBERATE. These pages are what a USER needs;
             this file is what a maintainer needs. Almost every choice in
             meq was settled by measurement, so the docs say which way it
             came out and that it was measured, and then tell the reader
             to measure it themselves wherever the choice is exposed as an
             option -- the trace solver, the assembly mode, the ordering,
             the globalisation, MKL_NUM_THREADS. The NUMBERS stay here,
             beside the tests that produce them, because a measurement in
             a manual goes stale silently.
```

Code style follows the sibling project MaNTA: **tab indentation**, Allman braces,
C++17. Naming, as decided for `../gffp` and carried here:

| | |
|---|---|
| Types — class, struct, enum, alias | `UpperCamelCase` |
| Enum values | `UpperCamelCase` |
| Functions and methods | `lowerCamelCase` |
| Variables, parameters, members | `lowerCamelCase` |
| **TOML configuration keys** | **`UpperCamelCase`** |

That last row is a deliberate mismatch with the C++ rule, not an oversight — it
comes from the same MaNTA convention set and applies to the key names in
`examples/*.toml` and to the string literals the parser looks them up by. Do not
reconcile the two.

**And external names keep their author's capitalisation, whatever the table
says.** `TraceSolver::cuDSS` is spelled the way NVIDIA spells it and carries a
`// NOLINT(readability-identifier-naming)` saying so; `UMFPack` and `Pardiso`
are spelled as MFEM's wrappers spell them and happen to need no suppression.
The house rule governs meq's own identifiers and does not extend to renaming
other people's products. This is the same exemption `.clang-tidy` already
records for MFEM-imposed overrides like `Eval` and `Mult`, and it is written in
both places because the naming check is what people meet first.

**This is enforced**, by `.clang-tidy`'s `readability-identifier-naming` and a
ctest named `naming`. It runs over `MEQ_CORE_SOURCES_PRESENT`, which is now
**every** file in `src/meq` — the exclusions this file used to describe were for
unported legacy sources and there are none left. `attic/` stays out permanently;
see `attic/free-boundary/README.md`.

`src/meq` deliberately keeps MFEM out of `Profiles` and `Source` — plain `double`
arguments, no `mfem::Vector` — so both are unit-testable without the library, and
the `mfem::Coefficient` adapters live with the assembly that needs them.

## Git

Branch **`main`**. The authoritative remote is **`origin` =
`github:ianabel/meq.git`** (private; `github` is a `Host` alias in
`~/.ssh/config`). A second remote `tantalum` points at the old
`tantalum:~/git/meq`, which is **unreachable from this machine** — do not expect
a fetch from it to work.

Those two had diverged, with neither a superset of the other, and were reconciled
by a real merge rather than by forcing one over the other. If old commits look
duplicated in the log, that is why.

**Tag `v0-legacy` is the pre-modernisation tree**, and it is where everything the
restructure deleted still lives: the one-off drivers (`mfemGS`, `mfemProjector`,
`mfemCheck`, `mfemLinearConvergence`, `FluxSurfaces`, `meshgen`, `KIntegrator`,
`BasicIO`, `NetCDFIO`), the `CurvedPoisson/` and `nonlinearPoisson/` experiments,
`output/`, `tools/plasma-output/`, and `HDGGSIntegrator.*`. Reach for it before
concluding something was lost.

**History is not rewritten**, deliberately: the deletions above are only
recoverable because it is not.

Commit messages end with the `Co-Authored-By` trailer.
