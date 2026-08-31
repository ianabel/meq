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
HDG-GS-1's Example 5 manufactured solution on the nonlinear one, with Newton
converging quadratically. `tests/convergence/SolovievConvergence.cpp` and
`NewtonConvergence.cpp` are the acceptance criteria and print the tables.
Concretely:

* `src/meq/GradShafranov.{hpp,cpp}` is ported — `DarcyForm` with
  `EnableHybridization` — is in `meq_core`, and is covered by the `naming` check
  like everything else.
* `src/meq/Solution.hpp` is **gone**, not ported. It wrapped hand-rolled block
  offsets and a second set of finite element spaces around what `DarcyForm` now
  owns, and `GradShafranovSolver` carries what was left of its job. Its
  `Prolong()` and `Update()` are stage-6 work against a mesh update that does not
  exist yet, and its `WriteOutputMFEM()` went with it. `src/meq/Output.{hpp,cpp}`
  replaced that job in stage 7c. `v0-legacy` has the original.
* `src/meq/Estimator.{hpp,cpp}` is written, compiles, is in `meq_core` and is
  under the `naming` check. `GridFunction::GetValueFacet`, which the old header
  called and which 4.9.1 does not have, is replaced by
  `traceFes->GetFaceElement(f)` + `GetFaceVDofs()` + `CalcShape()`, exploiting
  `DG_Interface_FECollection`'s `VALUE` map type — the same pattern
  `estimators_hdg.cpp` uses.
* `apps/meq.cpp` is the **driver**, and `MEQ_BUILD_APP` now defaults `ON`.
  `meq config.toml` parses, builds the mesh and source, solves — with the
  adaptive loop if it is asked for — prints a residual history and writes
  `.mesh`, `_psi.gf`, `_grad_psi.gf` and a `(R, Z)` NetCDF file. Exit codes
  0/1/2/3 as `DRIVER-PLAN.md` §5 specifies, and exit 2 is only reachable because
  `MFEM_USE_EXCEPTIONS` is now on.
* `Config`, `Profiles`, `Source` compile and test.
* **`ψ_ax` is an unknown of the non-linear system**, not an input.
  `meq::NormalisedSource` is the interface for a profile posed in normalised
  flux and `setSource( NormalisedSource &, double )` closes it by a bordered
  Newton — the trace and `ψ_ax` together, with the two non-local terms in the
  border. `HighBetaConvergence` is the acceptance criterion and it is green.
  `meq::NormalisedMHDSource` is the production source built on two `Profile`s;
  it is not yet reachable from a TOML file.

**meq is runnable.** `tests/convergence/DriverAcceptance.cpp` is the test that
says so, and it asserts the driver reproduces the *library* on the same
configuration — 1.189e-16 relative over 15,360 dofs — rather than comparing
against a closed form, for the reason recorded beside
`examples/soloviev-nstx.toml`.

**The curved boundary works through the driver.** `[boundary.shape]` builds the
shape, marks `D_h` out of the background mesh, locates `Γ_h` and hands a
`VertexConePath` to `setExtension()` — GS-2's technique, end to end, from a
config file. `examples/miller-curved.toml` is the worked example. Two checks
guard it, because the failure mode here is quiet: `Γ_h` carrying no transferred
datum silently imposes zero and still converges, so
`theDriverSolvesOnACurvedBoundary` pins the driver against the library (1.6e-16)
*and* against that zero-datum solve (2.7e-1, i.e. the transfer is doing real
work).

**The adaptive loop works through the driver, on both paths.**
`[adaptivity] Enabled = true` runs solve → post-process → estimate → mark →
refine, stopping at `TargetError` or `MaxIterations` and saying which, and prints
a per-cycle table. On the curved path it uses `meq::AdaptiveDomain` — the
companion mesh of GS-2 §3.3, without which `dist(Γ_h, Γ)/h_loc` doubles every
cycle and the transfer silently leaves the regime it is analysed in — and it
calls `setTransferredBoundary()` automatically, saying so in the log.
`examples/miller-adaptive.toml` is the worked example: η monotone
4.73e-4 → 2.47e-4 → 1.07e-4 → 6.87e-5 over four cycles, 97 → 449 elements, 0
transfer paths widened at any cycle — re-measured 2026-08-29, on `ψ*`. The
figures this file used to carry, 2.75e-3 → 4.77e-4 over 97 → 1069 elements, were
the same four cycles with `η` built on the raw potential, which is what the loop
had to do while MFEM's reconstruction was broken.
`theDriverRunsTheAdaptiveLoop` pins it
against the same loop driven through the library at **1.653e-16 over 6414 dofs**,
and separately asserts that it refined at all, that it refined *adaptively*, that
η came down, and that assumption P.1 survived the graded `Γ_h`.

**The driver's non-linear path is a reactive ladder**: Newton, and on *observed*
failure `PicardThenNewton`, rebuilding the solver because a caught
`ErrorException` leaves one unusable. Never predictive — nothing may be inferred
from `F` about which solver to run, for the reason recorded under *Why meq's
Newton struggles*.

**What the driver still does not do, and refuses rather than approximates**:
`[boundary] Type = "exact"` (needs a closed form `meq::Source` does not carry),
which exits 1 with an explanation.
The interpolating warm start of `DRIVER-PLAN.md` §4 needs GSLIB and is not
written; the exact restart — same mesh, same degree — is, and it is a
*first-cycle* guess only, since carrying an answer across a refinement is that
same missing interpolation.

**The loop builds `η` on `ψ*`, which it could not do until MFEM was fixed.** Four
of eq. (20)'s five terms are built on the post-processed potential, and while
`ReconstructFluxAndPot()` skipped its mean-value regularisation the driver ran on
`Potential::Raw` instead — one order down, correct but blunt. The fix landed; see
*Post-processing is back*. **A `TargetError` chosen against a pre-fix run is not
comparable**, because η is now a different and smaller quantity on the same mesh.

Stage 3's local post-processing was dropped as **meq code** and is not absent
from the calculation: `DarcyForm::Reconstruct()` supplies `ψ*` at `k+2` and the
estimator is built on it. See *Post-processing is back, and it was free*.

**And it is worth being clear about what the old code did, because the README
overstates it considerably.** Before the port, `meq` solved the *vacuum* coil
field and nothing else: `meq.cpp`'s right-hand side took `psi` and ignored it,
`Configuration::plasma` was hardcoded `nullptr`, and the profile loader was a
`{ return; }` stub that silently produced an empty spline. **meq had never solved
a Grad–Shafranov equation** until stage 2, and the semi-linear one still waits on
stage 4.

The four test files that existed were empty Boost stubs, and one of them —
asserting `foo == bar` with neither declared — could not compile. There was no
`test` target. Treat any claim in `README.md` about testing as aspirational until
it is checked; that file has not yet been rewritten.

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
cd build && ctest --output-on-failure       # ~870 s, 22/23
```

**ctest no longer needs any environment set by hand.** `tests/CMakeLists.txt`
puts both `MKL_THREADING_LAYER=GNU` and `MKL_NUM_THREADS=1` on every registered
test — the first for correctness, the second because without it the suite takes
well over an hour instead of thirteen minutes. Running a test binary **directly**
still needs both:

```sh
MKL_THREADING_LAYER=GNU MKL_NUM_THREADS=1 ./tests/SolovievConvergence
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
cd build-cov && MKL_THREADING_LAYER=GNU ctest
gcovr --root .. --filter 'src/meq/' --print-summary        # or --html-details
```

`MEQ_ENABLE_COVERAGE` is an option, not a build type, on purpose: a coverage
build wants `-O0`, and routing that through `CMAKE_BUILD_TYPE` would silently
drop `NDEBUG` for anyone who had set it, changing which `MFEM_ASSERT`s are live.
Stage 4 found a `BlockVector` size mismatch that *only* an assert catches, so
that difference is worth choosing rather than inheriting.

**CI cannot build the solver, and this is structural rather than an oversight.**
The MFEM meq needs is `meq-integration`, a local merge of four local branches,
published on no remote — so a hosted runner cannot obtain it and caching does not
help, there is nothing to fetch. (`ci.yml`'s own comment still names
`gf-hdg-subdomains-dev`, which was true when it was one branch.) `find_package(MFEM)` is therefore **not `REQUIRED`**: without
it, `CMakeLists.txt` builds the MFEM-free half — `Config`, `Profiles`, `Source`,
`SourceFactory` — and `tests/CMakeLists.txt`'s existing "unit not present"
guard already skips every convergence test without needing to know why.

So `.github/workflows/ci.yml` runs four unit suites, the `naming` check and a
coverage gate at 90% lines (94.8% measured). **That 94.8% is over four files,
not over the library** — reading it as a statement about meq is a mistake this
file used to invite.

**And the whole-library figure depends on which tool you ask, by a lot.** Two
denominators, on the same run:

| `src/meq/GradShafranov.cpp` | coverable lines | covered |
|---|---|---|
| `gcov`, executable lines only | 586 | **92.7%** |
| `gcovr`, its own line set | 838 | **64%** |

The file is 1554 lines for about 893 of code — this project comments heavily —
and `gcovr` counts roughly 250 lines that `gcov` does not consider executable at
all, most of which it then reports as uncovered. **So a line-coverage percentage
on this codebase is partly a measurement of comment density.** The CI gate uses
`gcovr`, so that is the operative number for the gate; `gcov`'s is the operative
number for the question "is this code tested".

On the honest denominator, measured over the whole suite:

| | |
|---|---|
| `Source.cpp` 100% · `Profiles.cpp` 99.3% · `BoundaryShape.cpp` 98.1% | |
| `Sampler.cpp` 97.3% · `SourceFactory.cpp` 95.2% | |
| **`GradShafranov.cpp` 92.7%** · `Config.cpp` 92.2% · `Estimator.cpp` 89.1% | |
| `WarmStart.cpp` 88.6% · `Output.cpp` 73.1% · `Field.cpp` 61.3% | |

`gcovr` over everything reads 82.2% of lines, 83.0% of functions and 48.4% of
branches.

**EVERY FIGURE IN THIS SUBSECTION PREDATES THE BORDERED NEWTON AND HAS NOT BEEN
RE-MEASURED.** `GradShafranov.cpp` gained about 350 lines and a whole solve path
in it; `Source.cpp` gained `NormalisedMHDSource`, which the CI gate's
`--filter src/meq/Source` does cover and which has unit tests of its own. So the
shape of the argument here — two denominators, comment density, branch coverage
being the figure with headroom — still stands, and the specific percentages
should be re-run before anyone quotes them. `gcovr` is still not installed on
this machine; see the recipe below.

**A claim this file carried briefly, and which was wrong, is worth leaving
recorded**: that the solver sat at 60% against 92–100% for the configuration
layer, so "the best-covered code in meq parses TOML and the least-covered is the
part every physical claim rests on". The 60% was `gcovr`'s number and the
comparison was against `gcov`-flavoured intuitions. On one denominator the solver
is 92.7% and `Config.cpp` is 92.2% — they are the same. The two weakest files are
`Field.cpp` and `Output.cpp`, both small. It is the same trap as the threaded-MKL
one: **a measurement about the tool, not about the code.**

What is true regardless is that 22 of `GradShafranovSolver`'s 51 methods were
called by nothing, and `SolverContract.cpp` was written for that rather than for
a percentage — it found `setSource()` silently accepting a second source of the
same kind, and value faults throwing `logic_error` where the constructor throws
`invalid_argument`. Branch coverage at 48% is the figure with real headroom in
it, and it is the one nobody has looked at.

**And `gcovr` is not installed on this machine**, so the recipe above fails with
`command not found`; a venv is the way round it, and it needs
`--merge-mode-functions=separate` or it aborts on `Config.hpp`'s inline
accessors appearing at two line numbers across translation units. Run it from the
repo root: `--root .` after a `cd` into the build directory silently matches
nothing and reports 0%.

**What it does not run is every
claim about rates** — the HDG assembly, Newton, the extension technique and the
estimator are all local-only, and the `full` job in that workflow is written but
`if: false` until either the branch is published or a self-hosted runner exists.
Treat a green CI badge as evidence about the configuration and profile layers
and about nothing else.

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
refined the wrong elements. `setTransferredBoundary()` excludes those faces and
`η` then converges at `k+1` on the extension path too. **That is an omission, not
a repair**: the proper fix is to evaluate `φ_h`, and **MFEM now offers the route
— `mfem::TransferredDatumCoefficient` in `fem/darcy/extension_hdg.hpp`.** So the
blocker is gone and the omission is meq's own; rebuilding `η₅` on it wants its
own convergence measurement rather than a switch, which is the only reason it has
not been done. ROADMAP item 3.

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
`HDGDiffusionIntegrator`'s built-in stabilisation is `{h⁻¹Q}`-scaled — the LDG
choice, not the papers'. Measured on the Solov'ev benchmark it costs a full order
in the flux: `q` converges at `k`, not `k+1`, while `ψ` still converges at `k+1`.
A study of `ψ` alone would have passed it. `meq::ConstantStabilization` is the
`HDGStabilization` hook that fixes it, installed with
`HDGDiffusionIntegrator::SetStabilization()`. Keeping `IsConstant()` true also
keeps meq out of the `EvalGrad` trap that header warns about, where a missing
derivative gives "no wrong answer, only slow Newton convergence".

**Getting a constant `τ` takes deliberate work — `HDGDiffusionIntegrator` will
not give you one by default.** Its documented stabilisation is

```
τ± = ( β ± (α/2)(u·n)/|u·n| ) { h⁻¹ Q }
```

— scaled by the inverse local mesh size *and* by the diffusion coefficient. That
is a reasonable default and it is not what the papers use. The designed way out
is the `HDGStabilization` hook in `fem/darcy/bilininteg_hdg.hpp`: subclass it,
return the constant from `Eval()`, install it with `SetStabilization()`. Read
`StabValue()` there — with a hook installed it divides the quadrature weight out
before calling and multiplies it back after, so returning a bare `τ` yields
exactly `⟨τψ, w⟩`. The MFEM tree has no ready-made constant implementation; the
only subclasses are test fixtures in `tests/unit/fem/test_bilininteg_hdg.cpp`
and `test_darcy_degenerate.cpp`, which are the idiom to copy.

This also matters for Newton. That header warns that omitting `EvalGrad` for a
*non-constant* stabilisation gives "no wrong answer, only slow Newton
convergence — a failure that survives a passing regression suite". A constant
`τ` has `IsConstant() == true`, so `EvalGrad` is never called and the trap does
not arise. One more reason to keep `τ` constant unless something measured says
otherwise, which is also the sibling project's standing advice: **do not derive
`τ` from the local coefficient.**

### Which MFEM, and why not master

**meq builds against `../mfem/install`** — MFEM **4.9.1** on branch
**`meq-integration`**, CMake-built in `../mfem/build` from sources in
`../mfem/mfem-src`. The options actually set, read out of
`share/mfem/config.mk` 2026-08-30 rather than remembered:

| | |
|---|---|
| `MFEM_USE_SUNDIALS` | KINSol, so `KINSolver(KIN_LINESEARCH)` is reachable. **Now `../sundials/cuda-install`, not `../sundials/install`** — see the CUDA row |
| `MFEM_USE_GSLIB` | `FindPointsGSLIB`; gslib v1.0.9 built alongside at `../mfem/gslib` |
| `MFEM_USE_SUITESPARSE` | UMFPACK, the direct solver meq's own solver runs on |
| `MFEM_USE_MKL_PARDISO` | `PardisoSolver`, oneAPI MKL 2026.1, **on the threaded `mkl_gnu_thread` layer since 2026-08-30**. It used to be resolved to Debian's `mkl_sequential`, which is a bigger deal than it sounds — see *Threaded MKL is a catastrophe here* |
| `MFEM_USE_EXCEPTIONS` | `MFEM_ERROR_THROW` by default, which is what makes the driver's exit code 2 reachable |
| `MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE` | **both, and both are now load bearing.** `DarcyHybridization::SetAssemblyMode( Threaded )` needs both and **aborts** rather than falling back without them. meq reaches it through `setAssemblyMode()`. This row used to read "not exploited"; that is no longer true |
| `MFEM_USE_LAPACK`, `MFEM_USE_ZLIB`, `MFEM_USE_DOUBLE`, `MFEM_USE_MEMALLOC` | |
| `MFEM_USE_MPI = NO` | meq is serial throughout; `../mfem-hdg-dev-par` exists for parallel work |
| **`MFEM_USE_CUDA = YES`**, `CUDA_ARCH = sm_75` | CUDA 13.3. **`sm_75` because the card is compute capability 7.5**; MFEM's cache defaulted to `sm_60`, which is wrong here and nothing warns you |
| **`MFEM_USE_CUDSS = YES`** | cuDSS 0.8.0 from Debian (`/usr/include/cudss.h`), found with `-DCUDSS_DIR=/usr`. `mfem::CuDSSSolver` agrees with UMFPACK to 3.5e-14; its **timings on this machine are not reproducible** and meq's solver does not use it. See *cuDSS* below |
| `RAJA`, `OCCA`, `CEED`, `SIMD`, `AMGX` all `NO` | |

**Turning CUDA on is not one flag, and two of the three obstacles are silent.**

* **SUNDIALS must itself be a CUDA build.** `linalg/sundials.hpp` hard-errors
  with `MFEM_USE_CUDA=TRUE requires SUNDIALS to be built with CUDA support`,
  testing `SUNDIALS_NVECTOR_CUDA`. meq cannot drop SUNDIALS — `KINSolver` backs
  `AndersonPicard` and `PicardThenNewton` — so the answer is a second SUNDIALS.
  `../sundials/cuda-install` is 7.5.0 configured `-DENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=75`, otherwise matching the old one.
* **CMake will not generate at all** with `MFEM_ENABLE_TESTING=ON` in a serial
  CUDA build. That is an upstream MFEM bug, present unchanged in `origin/master`
  and filed as `../mfem-hdg-dev/doc/CMAKE-SERIAL-CUDA-SUNDIALS.md`: the device
  test loop merges the serial and parallel example lists and never re-splits
  them, so it registers a test against `sundials_ex9p` in a build that never
  created it. **`-DMFEM_ENABLE_TESTING=OFF` is the workaround** and costs meq
  nothing, since meq does not run MFEM's ctest suite.
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
| `gf-hdg-linearise-first` | `SetNonlinearOrdering`, `SetAssemblyMode`, `SetGradientMode`, `SetLocalFactorMode` | merge 2 |
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

`doc/HDG-DEFECTS-FROM-MEQ.md`, `HDG-LINEARISE-THEN-CONDENSE.md`,
`HDG-ELEMENT-LOCAL-PARALLELISM.md` and `DIRECT-SOLVER-SYMBOLIC-REUSE.md` are what
meq puts there, and they are requests, not work.

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

**Also gone: `GridFunction::GetValueFacet`**, which `Estimator.hpp` calls. It was
a patch to the old branch and does not exist in 4.9.1.

### Post-processing is back, and it was free

Stage 3 dropped `ψ*` on the grounds that `q` already gives the physical quantity
at `k+1`, and noted that stage 6's estimator would need it again. It does — eq.
(20) uses `ψ*` in **four of its five terms** (`η₁`, `η₂`, `η₄`, `η₅`), not three
as an earlier version of this file said.

**`DarcyForm::Reconstruct()` supplies it, measured at `k+2`**: rates 3.03, 4.03,
5.00 for `k = 1, 2, 3`. No hand-written local solve was needed, and it survives
the extension path too. The old `GSSolver::Postprocess()` stays deleted.

**Through Newton it works, everywhere, since MFEM's reconstruction stopped
skipping its own close.** What follows is the history, kept because the failure
was silent and the shape of it is the useful part.

The fix is `gf-hdg-dev`'s *"The postprocessing closes on the element average,
always"*, merged into `meq-integration`. The local post-processing problem is a
pure Neumann one by construction — the total flux driving it is in H(div), so the
element's flux balance is already satisfied and the potential is determined only
up to a constant — and the element average is what closes it. That is
unconditional, and the code had been treating it as one of three choices.
Measured from meq's side on the case that used to read 20.3, 64.1 and 61.6 on
dead elements: **1.0069, against 1.0048 for elements where `∂F/∂ψ` does not
vanish at all.** `ψ*` is a post-processing everywhere, and the driver's adaptive
loop is back on the published estimator — `η` reaching 6.87e-5 on **449**
elements where the degraded one needed 1069 to reach 4.77e-4.

**What it was, and what it cost to find.** That is a narrowing of what this file used to say, arrived at
through two wrong diagnoses, and the mechanism is worth knowing exactly because
it is one condition.

`ReconstructFluxAndPot()` used to consult only the linear `M_p`. MFEM now lifts
the non-linear potential integrators as a Jacobian frozen at the computed
potential, and **measured on Example 5 that delivers `k+2`**: rates 3.05, 4.05,
5.03 at `k = 1, 2, 3`, and 47×, 113×, 125× smaller than `ψ_h` on the finest mesh.
`thePostProcessedPotentialSurvivesNewton` is that measurement.

**But the local post-processing problem is a pure Neumann problem — singular by
construction, determined only up to a constant — and MFEM's fix for that is
skipped by a flag set on the wrong condition.** `darcyform.cpp:1464` sets
`nl_src` on the mere *presence* of a non-convection non-linear integrator, under
a comment reading *"use non-singular terms as a source"* that states the intended
test exactly. Where `∂F/∂ψ` vanishes the integrator's Jacobian is the zero
matrix, so the term **is** singular; `nl_src` claims otherwise; the mean-value
branch at `:1586` is skipped; and `:1629` factors a singular matrix with
`Factor( m, TOL = 0.0 )` testing `|pivot| <= 0.0` and the return value discarded.

**It is not "a zero right-hand side it cannot handle".** The singularity is by
design and the handling is written. It is unreachable.

**And `nl_src` is per element, so the corruption is too.** Measured, with `∂F/∂ψ`
vanishing where `z < threshold`, comparing `max|ψ*|` against `max|ψ_h|` per
element:

| dead fraction | global `‖ψ*‖/‖ψ_h‖` | worst **dead** element | worst **live** element |
|---|---|---|---|
| 0.125 | 1.87 | **20.3** | 1.0049 |
| 0.312 | 5.69 | **64.1** | 1.0024 |
| 0.500 | 7.57 | **61.6** | 1.0021 |

At an eighth of the domain dead — a tabulated profile with a flat segment, which
is ordinary — elements are 20× wrong while the whole-domain norm is 1.87. **Any
global check misses it.** And a floor of **1e-12** on `∂F/∂ψ` — twelve orders
below everything else in the problem, incapable of moving a solution — takes the
ratio from 7.565317 to 0.998525. That is a singular matrix and nothing else.

**meq carried no runtime check for this**, deliberately — a solver should not
stand permanently on guard against its dependency, and the condition was one flag
test away from never arising. The state of the defect lived in the suite instead,
as a test asserting the *correct* behaviour and failing. That is what flipped
green when the fix landed, and what said to put the driver back on `ψ*`.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` still runs, now as
a regression rather than a tripwire. The report is
`../mfem-hdg-dev/doc/HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md`, and the fix
commit credits it.

**And the defect does not reach the solve — measured, not argued.** A number like
6.9e11 raises the question of whether `ψ_h` or `q_h` are also corrupted somewhere
less visible. They are not, and the reason is structural: the forward solve's
local problem takes its potential block from the HDG stabilisation on `Mnl_p`'s
**interior faces**, a fixed bilinear form whatever the source does, while
`Reconstruct()` builds a different local problem whose regularisation is what
fails. `theReconstructionDefectDoesNotReachTheSolve` compares both paths in `ψ_h`
**and** `q_h` over `k = 1…4` and three meshes:

| worst over the sweep | `ψ_h` | `q_h` |
|---|---|---|
| relative difference between the paths | **1.6e-13** | **2.6e-13** |

against discretisation errors from 2.2e-3 down to 7.3e-11. It also checks that
post-processing leaves both **bit-identical**, since the driver writes them
afterwards.

The measurement that justifies the `ψ*` requirement, rather than quoting it:
building `η₂` on raw `ψ_h` loses **exactly one order at every `k`** — 2.002 vs
0.998, 3.000 vs 2.002, 3.992 vs 2.981 — and is 124× to 407× larger on the finest
mesh. That is the defect the pre-modernisation estimator had.

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

| `ν` | `A` | `ψ_ax` | `ψ_ax − max ψ_h` | ratio to estimate | `max|∂F/∂ψ|/λ₁` | Newton |
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
border. **`NonlinearOrdering::LineariseThenCondense` used to be refused too and
no longer is** — that guard was written from `darcyhybridization.hpp`'s summary
of the ordering rather than from the code, which re-evaluates the source on every
residual; both orderings reach the same `ψ_ax` to every digit printed. See the
note on the guard in `GradShafranov.cpp` and `NORMALISED-LINEARISE-FIRST.md`.

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
| `max|∂F/∂ψ| / λ₁` | 7 | 7 | 8 | 13 | **26** |

Note the pedestal itself already sits at 7 and converges, so exceeding `λ₁` is
not by itself fatal; 26 is. **This supersedes an earlier reading of §4.4 as a
trivial-branch problem.** `F(r, 0) = 0` is true and the trivial branch is real —
see *Traps* — but it is not what defeats the iteration here, since the runs above
carry non-homogeneous ramp data that keeps `ψ` away from zero and fail anyway.

Do not file §4.4 under stiffness.

**§4.4 does have a solution, and continuation reaches it.** Adaptive continuation
in the added term's amplitude — halve the step on failure, grow it by 1.3 on
success — walked `c₃` from 0 to −18 at `k = 2, n = 32` in **9 solves with 2
retreats**, the final one taking 8 Newton iterations:

```
  c3 =  -1.800 : ok     7        c3 = -18.000 : FAIL  60   (halving)
  c3 =  -4.140 : ok     6        c3 = -17.887 : FAIL  60   (halving)
  c3 =  -7.182 : ok    10        c3 = -16.762 : ok    13
  c3 = -11.137 : ok     9        c3 = -18.000 : ok     8
  c3 = -15.637 : ok    10
```

Uniform steps do not do it — ten of them stall at `c₃ = −10.8` — so the step
control is the point. There is **no limit point on the branch**: the step never
had to fall below 1e-3, it simply had to shrink over the last sixth.

**THIS MUST NOT GO INTO THE DRIVER, AND THE REASON IS NOT TASTE.** What was
ramped is `c₃`, a constructor argument of the `CurrentHole` *test fixture*. The
`Source` interface exposes `f( r, z, ψ )` and `dFdPsi( r, z, ψ )` and nothing
else — **there is no amplitude parameter in it and no way to recover one from a
black-box `F`**. This continuation is therefore unavailable to a driver on
principle, not merely inadvisable.

The black-box analogue is scaling the whole source, `F_λ = λF`, which *is*
expressible. It is **untested and is a different homotopy**: at `λ = 0` it
degenerates to the harmonic problem, where the path above starts from the
converged *pedestal*, already a nontrivial equilibrium. Do not assume it behaves
the same.

**And it cannot currently be predicted that continuation is needed.** The ratio
`max|∂F/∂ψ| / λ₁` is computable black-box, `dFdPsi` being mandatory. But the
pedestal converges at 7 and the hole fails at 26 — two points, no threshold —
and the ratio needs the range of `ψ`, which is not known before solving. A
detector calibrated on that would be fitting noise.

So: **the driver gets a reactive ladder, never a predictive one.** Newton, and on
failure `PicardThenNewton`. Failure is *observed*, and needs nothing from `F`
beyond the existing interface. Continuation stays a test-only investigation. If
it is ever wanted generally, the prerequisite is a `Source` that can
parameterise itself, which is an interface change to argue on its own merits.

So the remedy for a production run is **resolution**, which means stage 6's
adaptive loop rather than a globalisation. What globalisation buys is the
**coarse start** an adaptive loop necessarily begins from, which is exactly where
raw Newton fails — see *Picard, then Newton*.

What survives below: hybridization really does put a nonlinear solve inside every
element elimination, and that is why an under-resolved §4.2 grinds where a CG
code would not. What does not survive is the conclusion that this makes the
problems unreachable. It makes an under-resolved discretisation expensive.

**Two hypotheses below are falsified and one is now known to be incomplete.**

Three Grad–Shafranov codes solve this equation by Newton and report it robust:

| | discretisation | nonlinear systems per Newton step |
|---|---|---|
| Serino, Tang, Tang, Kolev & Lipnikov (`refs/MFEM-GS-Newton.pdf`) | CG, `ψ ∈ H¹`, MFEM | **1, global** |
| CEDRES++ (`refs/CEDRES.pdf`) | CG finite elements | **1, global** |
| FreeGSNKE | 4th-order finite differences, Newton–Krylov | **1, global** |
| **meq** | **hybridized HDG** | **1 global + one per element** |

**In a CG or finite-difference discretisation `ψ` is one global unknown vector
and `F(ψ)` enters only the global residual.** Newton linearises once. Serino et
al. eq (3.1) is `a(ψ, v_i) = l(I, v_i)` and eq (3.7) is one linearised system per
step — that is the whole nonlinear structure.

**Hybridization changes that structure.** Static condensation expresses `(q, ψ)`
element by element in terms of `ψ̂`, and when `F` depends on `ψ` that elimination
is *itself a nonlinear solve, once per element per residual evaluation*. meq runs
`N_elements` independent Newtons inside every outer residual, none of them
globalised, and **any one of them failing poisons the whole residual**. That is
what `el: N not convered in 100 iters` is, and it is why globalising the outer
iteration does nothing — see *On SUNDIALS*.

**Two comfortable explanations are wrong, and the paper says so directly.**
Newton is not the problem: Serino et al. built their Newton solver precisely
because "conventional Picard-based solvers fail to converge" on the Taylor state,
and report the residual reaching 1e-6 "in a small handful of iterations". And
free boundary is not the problem either — they note they had already built a
fixed-boundary adaptive solver and that the fixed problem is "significantly
easier". **meq is doing the easier problem and finding it harder**, which is the
red flag that this section exists to explain.

**It also explains the GS papers' choice.** Sánchez-Vizuet and co. use
Anderson-accelerated Picard and keep `F` as opaque problem data so the solver
"relies only on the discretization of the toroidal operator `Δ*`". In a
hybridized method that is not fastidiousness — Picard evaluates `F` at the
previous iterate, which leaves **every local solve linear** and the whole
difficulty disappears. Their design is coherent; meq's Newton bought
`∂F/∂ψ` in the global Jacobian and paid for it with `N_elements` nonlinear
subproblems, and nothing measured that price until now.

**And the canonical HDG paper says the ordering is wrong.** Nguyen, Peraire &
Cockburn wrote the method this branch implements, and `refs/HDG-NPC-2.pdf` is
their treatment of *nonlinear* problems. §2.6:

> Next, we apply the Newton–Raphson method to solve the above system … we then
> find an increment `(δq_h, δu_h, δû_h)` …

Newton is applied to the **full** `(q, u, û)` system, giving eq (14) — a *linear*
system in the increments. The hybridization is then applied to **that**,
eqs (16)–(18), producing `K δΛ = F` for the trace increment alone, and the local
elimination is a matrix inverse of block-diagonal `A`, `B`, `D`: "the inverse can
be computed on each element independently … it results from applying the LDG
method to solve the **linearized** PDE".

**Every local operation is a linear solve. The canonical method has no
element-local nonlinear iteration at all.** meq has one because it condenses
first and linearises second, which is the opposite order.

**So the fix has a name: linearise, then condense.** Whether MFEM can be asked to
do it that way is the open question — `DarcyHybridization` is built around the
other order (`LocalNLOperator`, `SetLocalNLSolver`, `LSsolveType`), and its
`GetGradient` differentiates a residual that already required the local
nonlinear solves to converge. That is an MFEM-side capability question and
belongs in `../mfem-hdg-dev/doc/`, not here.

**What follows, in order of how targeted it is.**

1. **The local solver was never chosen.** `SetLocalNLSolver` offers `Newton`,
   `LBFGS` and `LBB`, and meq hardcoded `Newton` — undamped, on exactly the
   problems that are failing. First measurement on §4.2 at `k = 1, h = 0.05`:
   Newton 42 outer iterations, LBFGS 36, **LBB 25**. `setLocalSolver()` now
   exposes it. This is the cheapest lever and it was sitting unused.
2. **The control has been run, and it confirms this section.** Not by dropping
   hybridization — that gives a different method, LDG with a weakly imposed
   datum, so it would confound the discretisation with the solve. The variable
   is isolated instead by **Picard on meq's own linear path**: identical mesh,
   spaces, hybridization and `τ`, but `F` evaluated at the previous iterate and
   handed to `setSource( Coefficient & )`, which puts it on the right-hand side
   and makes the potential block linear. **Every local elimination is then a
   linear solve and there are no local Newtons at all.** On §4.2 at
   `k = 1, h = 0.05`, the case that fails:

   | method | local solves | outcome |
   |---|---|---|
   | Newton, meq's | **nonlinear**, one per element per residual | **fails** |
   | Picard, `ω = 1` | linear | stalls at 3.5e-1 |
   | **Picard, `ω = 0.5`** | **linear** | **2.8e-8 and falling** |

   So the mesh is fine, the discretisation is fine, and the problem *is* solvable
   there. What fails is Newton-with-nonlinear-local-solves specifically. Note the
   `ω = 1` row: undamped Picard stalls, which is the weakness CEDRES++ and
   Serino et al. both report and the reason the GS papers use **Anderson**-
   accelerated Picard rather than plain. It is also why this control needed the
   relaxation to mean anything — plain Picard failing would have proved nothing.

   **This is a diagnosis, not a recommendation.** Relaxed Picard took 200
   iterations to reach 2.8e-8 where Newton takes 42 on the mesh Newton manages.
   The point is what it isolates, not that meq should adopt it.
**THE ORDERING WAS NOT THE FIX. IT IS IMPLEMENTED, IT IS CORRECT, AND IT DOES
NOT HELP.** This section predicted that linearise-then-condense would rescue the
stiff sources. `gf-hdg-linearise-first` landed, two bugs meq reported against it
were fixed, and the prediction is **falsified by measurement**. `local 0` from
`GetNumLocalNLIterations()` confirms the ordering genuinely takes effect:

| case | condense-first | linearise-first |
|---|---|---|
| Example 5, similarity, `k = 2, 3` | ok, 3 it | ok, 3 it |
| §4.2 pedestal `k = 1`, `n = 16/24/30` | ok — 31/7/5 it | **fails at 60, all three** |
| §4.2 pedestal `k = 2, n = 16` | ok, 5 it | **fails at 60** |
| §4.3 barrier `k = 1` / `k = 2` | fails / ok, 10 it | **aborts** |
| §4.4 hole `k = 1` | aborts | aborts |

Identical on the mild problems, strictly worse on every stiff one.

**RE-MEASURED 2026-08-30 against the current library, because that table was
taken before `gf-hdg-dev` was merged and before two fixes landed in exactly this
code path, and a stale verdict on an ordering meq keeps is worth more than it
costs to re-run.** The pedestal reproduces, at every resolution:

| §4.2 pedestal | condense-first | linearise-first |
|---|---|---|
| `k = 1, n = 16` | ok, **39** it | fails at 60 |
| `k = 1, n = 24` | ok, 7 it | fails at 60 |
| `k = 1, n = 30` | ok, 5 it | fails at 60 |
| `k = 2, n = 16` | ok, 5 it | fails at 60 |

39 against the 31 recorded above is the same knife edge `k = 1, n = 16` always
is — threaded-MKL reduction order, see the pedestal tripwire — and the other
three match to the iteration. (The §4.3 rows did **not** reproduce in that run:
the barrier failed under condense-first at both orders where this table says
`k = 2` takes 10. The scratch harness differed from `PedestalConvergence`'s in
some detail that was not chased, since the pedestal is what the verdict rests on.
`PedestalConvergence` itself passes.)

**So the verdict stands on current measurement and not on memory.** What is
missing is still an *explanation*, and there are now three dead hypotheses:

1. **The retained linearisation.** `NewtonSolver` evaluates the residual at the
   new iterate using the linearisation retained at the *old* one — but wrapping
   the operator so the residual relinearises at its own argument changes nothing.
2. **A frozen local state.** Plausible from the outside, and wrong: the branch
   keeps `lin_u`, `lin_p` and `lin_trace` separately from the
   `FormLinearSystem`-time seed, and `MultInvLin()` takes one local correction
   per residual and two per gradient precisely so that the linearisation point
   advances. Its own comment records a bug already fixed there, with the unit
   test that covers it.
3. **Outer globalisation.** A line search on the trace makes it *worse*, failing
   at 2 where the plain iteration reaches 60.

**THE MECHANISM, MEASURED: THE REDUCED RESIDUAL IS NOT A FUNCTION OF THE TRACE.**

An earlier version of this section said the local Newton correction overshoots.
That was loose — a step cannot overshoot from *wherever* it starts if the
function is smooth and the Jacobian consistent, which is the right objection and
is what sent this back to the measurements. The answer is sharper and it is a
property of the operator, not of a step.

Three quantities, at a wiggled ramp on §4.2 at `k = 1, n = 8`, as `σ²` sharpens
the local problem. **same-lin** is `‖R(x) − R(x)‖/‖R(x)‖` for two evaluations at
one trace with a relinearisation *at that same trace* between them. **cross-lin**
is the same thing with the relinearisation taken at a *different* trace — which
is precisely what `NewtonSolver` does, evaluating `R` at `x_k` while carrying the
linearisation from `x_{k−1}` and only then relinearising. **J vs central** is the
assembled gradient against a central difference, both at one linearisation.

| `σ²` | | same-lin | **cross-lin** | J vs central |
|---|---|---|---|---|
| 0.05 | condense-first | 0 | **0** | 2.7e−11 |
| 0.05 | linearise-first | 0 | **5.6e−04** | 2.3e−11 |
| 0.01 | condense-first | 0 | **0** | 1.02e+02 |
| 0.01 | linearise-first | 0 | **1.75e−01** | 2.5e−06 |
| 0.005 | condense-first | 0 | **0** | 1.26e+02 |
| 0.005 | linearise-first | 0 | **1.49e+00** | 5.0e−09 |

**Under condense-first the residual is exactly independent of the linearisation
history** — 0.000e+00, at every width. It must be: the local problems are solved
to convergence, so nothing about how they were reached survives.

**Under linearise-first it is not, and the dependence grows with local
stiffness** — to 149% at the published `σ² = 0.005`. The retained local fields
are a **hidden state variable of the operator**. `MultInvLin()` advances them by
two corrections when relinearising and one when evaluating, from wherever they
already were, so `R` is really `R( trace; fields )` and the fields carry history.
Newton is handed a function that changes under it by more than its own value
between iterations, and *no* Newton converges on that — which is why the failure
survives refinement, survives a line search, and survives Picard.

**So it is not the Jacobian.** At a fixed linearisation the gradient is the
derivative of the residual to 2.3e−11, 2.5e−6 and 5.0e−9 — cleaner, at the stiff
widths, than condense-first manages.

**And the 1.02e+02 in the condense-first column is a second finding, about the
ordering meq actually ships.** There the residual is a perfect function of the
trace but the gradient disagrees with a central difference of it by a factor of
100 once `σ² ≤ 0.01`. That is the element-local solves hitting their 100-iteration
cap: a residual computed from an unconverged local solve is not smooth, so the
difference is meaningless and the Jacobian may or may not be. **It is the
measured explanation for the pedestal's 39 wandering iterations at
`k = 1, n = 16`**, and it is the same defect that corrupts a differenced border —
see the seed note under *Traps*. Refinement cures it because a resolved local
problem converges inside its cap.

**What each ordering therefore needs, and they are different things.**
Linearise-first needs the fields to stop being hidden state: either the outer
solver owns the field increments too — which is what NPC's simultaneous Newton
actually is, and would let one line search damp the whole step — or the local
correction is iterated until the fields no longer move, which is condense-first
under another name. Condense-first needs its local solves to converge, which
`setLocalSolver()` and resolution both reach for and neither guarantees. **Both
are MFEM-side.**

**And Picard-then-Newton does not rescue it**, which is worth recording because
it is the obvious thing to try and it is the route that rescues condense-first on
these same sources. Measured on §4.2 at `k = 1, n = 16`: `PicardThenNewton` under
condense-first converges in 4 Newton steps from `‖r₀‖ = 1.88e+00`; under
linearise-first the same handoff starts at `‖r₀‖ = 4.89e+02` — from a state
Picard has already walked to the answer — and thrashes between 1e−1 and 1e4 for
60 iterations. Picard fixes where the iterate is; it does nothing about a local
correction that overshoots from wherever it starts.

**So the element-local nonlinear solves were never the cause.** The control that
pointed here — relaxed Picard converging where Newton fails — differed in
*globalisation* as well as in local linearity, and the globalisation cross says
which half was doing the work. Every KINSOL strategy that steers the outer trace
buys it by driving the local Newtons harder and still fails:

| §4.2 pedestal `k = 1, n = 16` | outcome | local nonlinear iterations |
|---|---|---|
| condense-first, plain Newton | **ok, 31 it** | 133,168 |
| condense-first + `KIN_LINESEARCH` | fails at 45 | **1,381,527** |
| condense-first + `KIN_NONE` | fails at 19 | 105,021 |
| linearise-first, plain | fails at 60 | 0 |
| linearise-first + line search | fails at **2** | 0 |

**What does work is Picard, then Newton** — see *Picard, then Newton* below. Keep
`setNonlinearOrdering()`: it costs nothing, it is the canonical NPC ordering, and
it is the only way to run this discretisation without element-local nonlinear
solves at all. Do not expect it to fix a stiff source.

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

`mfem::KINSolver` **derives from `mfem::NewtonSolver`** (`linalg/sundials.hpp`),
and is reachable through `setGlobalisation()`. `SetJFNK` and
`EnableAndersonAcc` come with it.

**But it is not a drop-in, and the difference is silent.** This file used to say
that code written against a `NewtonSolver&` takes either "with no abstraction
layer and no rewrite". Not so: `NewtonSolver::Mult( b, x )` forms
`r = oper(x) − b`, while **`KINSolver::Mult` declares its first argument without
a name** and solves `oper(x) = 0`. meq's trace right-hand side is not zero, so
handing it straight to KINSOL converges — to the solution of a different
problem. `ShiftedResidual` in `GradShafranov.cpp` is the adapter; it reproduces
`NewtonSolver`'s residual exactly, which is what makes a comparison between the
two paths mean anything. `kinsolAgreesWithNewtonWhereBothConverge` is what
watches it, and would catch the adapter being dropped or its sign flipped.

**And the line search does not fix the pedestal.** Measured:

| case | Newton | `KIN_LINESEARCH` | `KIN_NONE` |
|---|---|---|---|
| §4.2 `k=1, n=16` (`h = 0.05`) | ok, 42 it | **fails at 18** | fails at 6 |
| §4.2 `k=1, n=24` (`h = 0.033`) | ok, 23 it | ok, 22 it | fails at 11 |
| §4.3 homogeneous + a guess | trivial branch | trivial branch | — |

It rescues nothing, and on the one case that motivated it it is *worse* than the
undamped iteration.

**The reason is structural and worth internalising: globalising the outer
iteration does not globalise the inner ones.** The failure is
`el: N not convered in 100 iters` — MFEM's **element-local** nonlinear solve,
one per element per residual evaluation, eliminating flux and potential for a
given trace. A line search chooses how far to move the *trace*; it cannot make
the local problems at that trace well posed, and KINSOL never sees them.

What would: damping the **local** solves, which is MFEM's to offer;
continuation in the source amplitude, so each solve starts from the previous
answer rather than a cold trace; or **Picard on the outer loop**, which
evaluates `F` at the previous iterate and leaves every local problem *linear*.
That last is what both papers do, and this measurement is the argument for it —
see *Newton, and the obligation it creates*, which is the other side of the same
trade.

This needed `MFEM_USE_SUNDIALS`, which the tree did not have when the section was
written, and which turned out to be nothing worse than a flag: MFEM 4.9.1's
SUNDIALS interface is fully modernised (`sunrealtype`, `SUNContext`) and asks
only for SUNDIALS ≥ 5, and 7.5.0 is installed at
`/home/ian/projects/sundials/install`. That mismatch is precisely why the *old*
pinned MFEM 4.5.1 stopped compiling — its `sundials.hpp` still used `realtype`
and `booleantype`, which SUNDIALS 7 removed.

**SUNDIALS is built in** (`../mfem/install`), so `KINSolver` is reachable
without a rebuild. It is no longer the blocker it was, and one reason is worth
recording because it was measured rather than guessed.

**`MFEM_USE_LAPACK` appears to fix the pedestal, and it does not. It tips a
marginal iteration by changing the rounding.** This is worth reading before
anyone concludes the globalisation problem has gone away, because the surface
evidence says it has.

Stage 6 recorded three failures: the pressure pedestal not converging at
`k = 1` for `h ≥ 0.05` (`el: N not convered in 100 iters` from MFEM's
element-local nonlinear solves), and the current hole (§4.4) going NaN and
*aborting the process*. Against a LAPACK build the pedestal converges at
`k = 1, h = 0.05` in 42 iterations; against `../mfem/install-nolapack`, which
is identical but for that one flag, it fails at 60 as before.

**But the flag is not doing what that suggests.** At `h = 0.03333` — the
well-posed mesh — both builds converge in **23 iterations, the same number**.
Only the marginal mesh moves, and 42 iterations is itself grinding against the
five a healthy Newton takes. Two points on one knife edge, not a robustness
gain.

There is no mechanism for it to be anything else. meq sets `LPrecType::LU`, so
the local solve really is the dense LU; both implementations partial-pivot on
the largest `|a|` in the column, which is the same rule; and the singularity
test that looks like a difference is not one — `Factor( int m, real_t TOL = 0.0 )`
compares `abs( pivot ) <= TOL` against a default of **zero**, which is `dgetrf`'s
condition exactly, and `DenseMatrixInverse::Factor` discards the return value
either way. What is left is `dgetrf`/`dgetrs`/`dgemm` being blocked BLAS-3
where MFEM's fallbacks are unblocked scalar loops: identical arithmetic, different
summation order, `O(1e-16)`.

**And the BLAS here is threaded MKL** — `MFEM_EXT_LIBS` carries
`-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5`, the same library the
`MKL_THREADING_LAYER` trap is about. Reduction order in a threaded BLAS depends
on the thread count, so whether that local Newton converges is not merely
fragile but machine- and environment-dependent. Do not treat 42 iterations as
reproducible.

**So the earlier claim in this file that `MFEM_USE_SUNDIALS = NO` was the
blocker stands, and an intermediate claim that LAPACK had fixed a third of it
was wrong.** Both §4.2 at `h = 0.05` and §4.4 still need globalisation;
`KINSolver(KIN_LINESEARCH)` is now available to try on them. What LAPACK
bought is a marginal case that happens to fall the right side of the line on
this machine today.

**`MFEM_USE_EXCEPTIONS` is now enabled**, and it does what `DRIVER-PLAN.md` §5
needs. With it, `MFEM_ERROR_THROW` is the *default* error action and
`mfem::ErrorException` derives from `std::exception`, so **no meq-side change was
required** — an existing `catch ( std::exception const & )` already catches it.

Measured on §4.4 at `k = 1, n = 16`, which used to take the process down with
SIGABRT:

```
 --> norm = -nan
 ... in function: virtual void mfem::NewtonSolver::Mult( ... )
RESULT|4.4 hole|k=1|n=16|condense-first|none|seeded|FAIL|0|local 214583
```

Exit 0, a reported failure, a usable iteration count. A driver can now return
exit code 2 rather than dying.

**A caveat worth keeping in view**: the throw unwinds out of the middle of
MFEM — that NaN is detected inside `NewtonSolver::Mult`, and the element-local
failures throw from deeper still. The objects are left as the throw found them.
meq's paths construct a fresh solver per solve and so do not care, but do not
assume a `GradShafranovSolver` is reusable after a caught `ErrorException`.

### The suite is 22/23, and the one red assertion is characterised

Measured 2026-08-31 against `meq-integration` carrying MFEM's `ad04da3749`, with
`MKL_NUM_THREADS=1`. **872 s for the whole suite.**

**`PedestalConvergence` is the only failing file and it has ONE assertion in it,
down from seven.** MFEM fixed the cause of the parity gap: the linearisation
point took a fixed **two** frozen-Jacobian local corrections, and
`GetGradient()` is the Schur complement of the Jacobian at whatever fields that
left — so the correction count *was* the gradient's accuracy, and nothing said
so. It now iterates to `setLocalSolver()`'s tolerance. What that bought, on
meq's own cases:

| GS-2 case | | before | after |
|---|---|---|---|
| §4.2 pedestal | `k = 1, n = 16` | fail at 60 | **ok, 28** |
| §4.2 pedestal | `k = 2, n = 16` | fail at 60 | **ok, 8** (against 11 for condense-first) |
| §4.2 pedestal | `k = 1, h = 0.0333` | fail at 60 | **ok, 11** |
| §4.3 barrier | `k = 2, h = 0.05` | fail at **0** | **ok** |
| §4.5 layer | `k = 1, h = 0.0333` | fail at 60 | **ok** |
| **§4.5 layer** | **`k = 2, h = 0.05`** | fail at 60 | **still fails at 60** |

**And the `0 iterations` mystery resolved itself.** That signature — a throw out
of MFEM before completing one iteration, rather than wandering for sixty — was
the one thing nobody could explain and nobody could reproduce outside meq. It is
gone, as a side effect of the divergence guard added to `MultInvLin()`. Upstream
flagged that as unverified and asked meq to check; meq checked.

**The survivor is measured against both orderings rather than merely reported**,
which is what makes it a parity claim:

| §4.5 internal layer | condense-first | linearise-first |
|---|---|---|
| `k = 1, n = 24` | ok, 23 | ok, 42 |
| **`k = 2, n = 16`** | **ok, 10** | **fail at 60** |
| `k = 2, n = 24` | ok, 12 | ok, 11 |

Isolated: refining once cures it and so does dropping an order. It is one of
three surviving out of 144 configurations upstream, all at widths where the
frozen-Jacobian local correction cannot converge and the guard truncates.
Closing them needs the local step globalised or solved exactly, and **nothing
has been tried** — meq is not asking for it, because the ladder covers the
production path.

**THE FIX COSTS WALL CLOCK, AND THAT IS THE TRADE RATHER THAN A REGRESSION.**
The correction loop that bought the correctness runs 4.2 to 12.1 corrections per
element per linearisation against the fixed two. `PedestalConvergence` went
351 s → 572 s and the whole suite 682 s → 872 s with no other change. Upstream
measures linearise-first at 0.9 s against condense-first's 0.7 s on a stiff case
and says so in its own doxygen: **this ordering is not a wall-clock win on a
stiff problem.** What it buys is a local problem that is always a linear solve
against one factorisation, and a reduced operator whose gradient is the
assembled Schur complement.

**A test that was green went red, and chasing it was worth more than the fix
to it.** `andersonPicardReachesTheSameSolutionAsNewton` compared Newton against
Anderson-Picard on one mesh at a tolerance of 1e-6. It now reads 9.1e-05 there.
Over a sweep:

| n | relative in `max ψ` |
|---|---|
| 16 | 9.103e-05 |
| 24 | **1.160e-13** |
| 32 | 3.275e-06 |
| 48 | **4.944e-13** |

**It is not a stopping tolerance**, which was the first guess: tightening rtol
from 1e-8 to 1e-12 leaves the `n = 16` figures **bit identical**, at 230
Anderson iterations instead of 162. Both iterations are fully converged and
their fixed points differ. Round-off agreement on some meshes and 1e-5 on
others, with no trend in the mesh, says these coarse discretisations carry **more
than one solution** — which is what this file records elsewhere about the GS-2
sources — and that two iterations with different non-linear structure sometimes
land on different ones. The old single-mesh gate was green because the pre-fix
Newton took 134 iterations on a gradient that did not belong to its residual and
drifted onto Picard's branch. **A better Jacobian converging faster to a nearer
solution is the expected consequence of the fix.**

The test now sweeps three meshes and asserts the two things actually entailed:
the **best** agreement is at round-off, which is what says the Picard path
solves meq's problem rather than a neighbouring one — a wrong sign or a missing
`1/r` would show O(1) — and the **worst** is bounded well below a different
problem. A per-mesh gate at 1e-6 is not reinstatable and the test says why.

**What the earlier `Mult()` fix cleared** is the Jacobian-check ordering:
`NewtonConvergence`'s central difference against the assembled gradient reads
4e-11, the `O(step²)` floor. meq keeps the hoisted `GetGradient()` regardless,
because holding the linearisation fixed across a difference is the right thing
to write whether or not the library requires it — though note upstream's finding
that re-taking the gradient *after* a difference silently buys extra corrections
and can hide exactly this class of defect.


**Do not compare these timings against the ones this file used to carry.** The
old 2269 s baseline was taken with a *different default ordering*
(`CondenseThenLinearise`), against a different MFEM, and with test sources that
have since changed substantially — `HighBetaConvergence.cpp` alone by 703 lines.
Most of the speedup is `MKL_NUM_THREADS=1` and that is separately and cleanly
measured under *Traps*; attributing the per-file ratios to any one cause would
be wrong. The table below is context, not a measurement:

| | old, condense-first | 08-30, MKL=1 | 08-31, MFEM fixed |
|---|---|---|---|
| `PedestalConvergence` | 1334 s | 409 s (7 red) | 572 s (1 red) |
| `naming` | 127 s | 150 s | 128 s |
| `MillerConvergence` | 194 s | 77 s | 58 s |
| `HighBetaConvergence` | 224 s | 25 s | 17 s |
| `NewtonConvergence` | 117 s | 36 s | 30 s |
| `WarmStartConvergence` | 84 s | 5 s | 5 s |
| `ExtensionConvergence` | 75 s | 36 s | 28 s |
| **total** | **2269 s** | **798 s** | **872 s** |

The last column is *slower overall* and *faster in every file but one*: the
correction loop that fixed the parity gap costs `PedestalConvergence` 163 s, and
that one file outweighs the gains everywhere else.

`naming` is now the second most expensive item in the suite and is the only one
that went *up*; it is clang-tidy over `src/meq`, so it tracks the library
growing rather than anything about the solver.

The rest of this section predates that run and is kept for what it records about
individual regressions.

**`HighBetaConvergence` was the last red one and it is green.** It asserted the
physical answer for a profile posed in normalised flux and failed while `ψ_ax`
was outside the residual; `ψ_ax` is now an unknown of the non-linear system and
the six cases converge to a self-consistent equilibrium in 4 to 11 Newton steps.
See *Newton, and the obligation it creates*. The file runs in 183 s.

**The `ψ*`-on-a-curved-boundary regression is fixed.**
`ExtensionConvergence::thePostProcessedPotentialConvergesAtKPlusTwo` and
`AdaptiveRefinement` both pass. That regression was `ψ*` converging at 1.28 at
`k = 1` and 1.16 at `k = 2` — a rate independent of `k`, so a fixed `O(h)`
geometric error rather than a discretisation one — and being *larger* than `ψ_h`,
by 7× at `k = 1` and 1000× at `k = 2`, so the post-processing cost accuracy
rather than buying it. `AdaptiveRefinement` failed downstream of it, `η₂`, `η₄`
and `η₅` being built on `ψ*`, which is why it surfaced as `η` refusing to come
down four files from the cause. Both are green now.

**The pedestal tripwire is rewritten and now passes.** It used to assert that
Newton does not converge at `k = 1, h = 0.05`, told its reader to delete it when
that stopped holding, and blamed Picard-vs-Newton and SUNDIALS in its prose —
both of which are now known to be wrong. It was never entitled to assert that
point in either direction: 42 iterations against LAPACK and a failure at 60
against `install-nolapack` is the same library differing by one flag, decided by
threaded-MKL reduction order, so an assertion on it is an assertion about this
machine today.

`pedestalConvergenceIsAResolutionThreshold` asserts the two **cures** instead,
each with a factor of four in it and neither marginal:

| | | |
|---|---|---|
| `k = 1, n = 16` | 42 iterations | the knife edge, recorded and *not* asserted on |
| `k = 1, n = 32` | 9 iterations | `h`-refinement cures it |
| `k = 2, n = 16` | 11 iterations | `p`-refinement cures it, mesh untouched |

That two independent refinement paths both cure it is what identifies the cause
as under-resolution; §4.4 is the control, and it fails at every order and mesh
tried because its trouble is that the problem is multi-valued. The threshold
itself is still asserted, as `iterations(n = 16) ≥ 2 × iterations(n = 32)`, which
holds whichever side the rounding falls on — 42 or a failure at the cap both
clear it. §4.4, the current hole, still fails, and now reports rather than
aborting.

**That file is 758 s and is the largest single item in the suite by a wide
margin**, which is worth knowing before adding to it. The percentage this used to
quote is gone rather than updated: the denominator moved when
`HighBetaConvergence` went from 428 s to 183 s, and the only re-measurement
available was taken on a contended machine, which by this file's own standing
rule is not a measurement about the code.

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

**The third one existed because meq used Newton rather than Picard, and
`NonlinearOrdering::LineariseThenCondense` has since removed it.** Picard
evaluates `F` at the previous iterate and leaves every local problem linear —
one dense factorisation and done. Condense-first Newton puts `∂F/∂ψ` inside the
local problem and makes it nonlinear. Linearise-first, which is **meq's default
and only ordering**, differentiates the full `(q, ψ, ψ̂)` system first and
hybridizes the linear system that results, so every element-local operation is a
linear solve and `GetNumLocalNLIterations()` stays at **zero**.

So the row above is what meq *has*, not what it runs: reaching it means asking
for `CondenseThenLinearise` deliberately, which only `PedestalConvergence` does,
and only to measure the two against each other. The historical reading — `el: N
not convered in 100 iters` being where the pedestal fails — describes the
ordering meq no longer uses. What the pedestal fails on now is the parity gap;
see *The suite is 22/23*.

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

**`tests/performance/` is the harness, and it is new.** `TraceSolverScaling` is
built but deliberately **not registered as a ctest**: everything it reports is a
timing, and this project's standing rule is that a threaded timing on this
machine is a measurement about the machine. It returns non-zero only for the two
*correctness* properties that make its timings mean anything — threaded assembly
reproducing serial assembly bit for bit, and PARDISO reproducing UMFPACK.
`scan.sh` drives it, one process per point, because MKL fixes its threading at
first use and an in-process sweep would measure the first setting five times.

**Two knobs, swept separately, and sweeping them together is actively
misleading:**

| | drives |
|---|---|
| `OMP_NUM_THREADS` | `DarcyHybridization`'s threaded element assembly |
| `MKL_NUM_THREADS` | UMFPACK's BLAS, PARDISO's internals, **and MFEM's element-local dense LU** |

The third entry in that second row is the one that surprises, and it is why the
axes have to be separated: with both swept together, UMFPACK's collapse and the
`k = 3` assembly blow-up swamp every other column in the table.

#### Threaded assembly: real, bit-exact, and worth about 1.2x

`DarcyHybridization::SetAssemblyMode( Threaded )` threads the element-local half
of `ComputeH()`. meq reaches it through `setAssemblyMode()`, which refuses
rather than letting MFEM abort when the build lacks OpenMP or thread safety.

**It is bit for bit**, which is asserted rather than trusted:
`SolverContract::threadedAssemblyReproducesSerialAssemblyExactly` compares every
stored entry and every column index over two degrees and two meshes and requires
`0.000e+00`. Not a tolerance — MFEM documents exactness, the mechanism is
specific (element-local arithmetic reassociates nothing, the scatter stays
serial and in element order), and if it ever needs slack the mechanism has
changed.

Speedup on meq's own problem, `MKL_NUM_THREADS=1` throughout, best of three:

| OMP threads | `k=2, n=64` | `k=3, n=64` | best over all eight cases |
|---|---|---|---|
| 1 | 0.86x | 1.03x | 1.04x |
| 2 | 1.09x | 1.10x | 1.16x |
| 4 | 1.17x | — | 1.33x |
| 8 | 1.15x | 1.22x | 1.26x |

**MFEM measures 2.09x at 8 threads on a comparable mesh and both numbers are
right.** MFEM times the hybridized assembly *alone*; meq times all of
`prepare()`, of which only `ComputeH()`'s element loop threads — the integrator
assembly, the scatter, and the FE space construction are serial. 1.2x is Amdahl
on the quantity a caller actually waits for. The scatter is the ceiling and MFEM
says so: it cannot be threaded, because an unfinalized `SparseMatrix` carries
one `current_row` and one column-pointer scratch for the whole matrix, so two
threads writing provably disjoint rows still collide, and the failure is a hang.

**So it is worth switching on SOMETIMES, and it is not the win TODO expected.**
That item called threading the element-local loops "the largest structural win";
measured, it is 20% on a mesh assembled a few times — and a **1.8x loss** on a
small system assembled hundreds of times, which is why it is not the default.
See *The defaults* below before enabling it.

#### Threaded MKL: a catastrophe, and not where it looks

Full account under *Traps*, including the wrong first diagnosis and why no
whole-test timing could have caught it. The short of it:

* `MKL_NUM_THREADS=16` costs `SolovievConvergence` a factor of **143** and
  `FieldConvergence` a factor of **180**, at 273% CPU.
* **It is not UMFPACK.** Separated out, `UMFPackSolver::SetOperator` degrades by
  about 40% across the whole thread range, never more.
* **It is `ComputeH()`'s element-local dense LU**, on blocks of order 10-30,
  through LAPACK. `k = 2` is untouched; `k = 3` degrades **forty-fold**.

#### The trace solve: PARDISO wins, and wins more with threads

`k = 2`, best of three, one process per thread count, assembly timing off so the
solver columns are not buried by the effect above. Agreement with UMFPACK is
**1.0e-14 or better** at every point, which is the only thing asserted.

`SetOperator` (analyse + factor), seconds:

| trace dofs | UMF `MKL=1` | UMF `=16` | PAR `MKL=1` | PAR `=4` | PAR `=8` | PAR `=16` |
|---|---|---|---|---|---|---|
| 2,400 | 0.0094 | 0.0110 | 0.0069 | 0.0056 | 0.0038 | 0.0058 |
| 9,408 | 0.0408 | 0.0618 | 0.0298 | 0.0162 | 0.0220 | 0.0163 |
| 21,024 | 0.1142 | 0.1550 | 0.0661 | 0.0558 | 0.0411 | 0.0456 |
| 37,248 | **0.2091** | 0.2893 | 0.1395 | 0.0903 | 0.0745 | **0.0738** |

`Mult` (backsolve), seconds:

| trace dofs | UMF `MKL=1` | PAR `MKL=1` | PAR `=8` | PAR `=16` |
|---|---|---|---|---|
| 9,408 | 0.0050 | 0.0028 | 0.0019 | 0.0014 |
| 21,024 | 0.0147 | 0.0103 | 0.0052 | 0.0046 |
| 37,248 | **0.0266** | 0.0188 | 0.0096 | **0.0085** |

**Three things follow, and the third is the awkward one.**

**PARDISO beats UMFPACK even sequentially** — 1.50x on setup and 1.41x on the
solve at 37,248 dofs, both at `MKL=1`. That is a fair fight with no threading on
either side, and it reframes the "3.8x setup, 2-3x solve" this file used to
carry: those were taken when the stray `libmkl_sequential` made PARDISO
sequential too, so they were never a threading result and the ratio has since
come down as the measurement got cleaner.

**PARDISO does scale, to about 8 threads.** 1.87x on setup and 1.96x on the
solve from 1 to 8, then flat — 16 threads buys 1.89x, i.e. nothing more. Against
UMFPACK at its own best setting the end-to-end gap is **2.83x on setup and 3.13x
on the solve**.

**And meq cannot have both at once.** `MKL_NUM_THREADS` is process-wide. The
setting that makes PARDISO fast is the setting that makes `ComputeH()` forty
times slower at `k = 3`. In a solve where assembly is 0.47 s and the trace solve
is 0.21 s, threading MKL to speed up a 0.21 s factorisation while paying 40x on
the 0.47 s assembly is a catastrophic trade. **So `MKL_NUM_THREADS=1` stands,
and PARDISO's scaling is unreachable until the element-local factorisation stops
going through threaded MKL** — which is `LocalFactorMode::Batched`'s territory,
or an `mkl_set_num_threads_local()` around the trace solve. Neither is done.


#### cuDSS: correct here, and not measurable here

`MFEM_USE_CUDSS` with cuDSS 0.8.0, `mfem::CuDSSSolver`, `NONSYMMETRIC` +
`FULL` — the general setting, chosen for the same reason PARDISO gets
`REAL_STRUCTURE_SYMMETRIC`: the trace matrix is structurally symmetric on both
paths and symmetric in value only on the fitted one. Reached with
`--device cuda`, which constructs an `mfem::Device` before anything else.

**Correctness: established.** cuDSS agrees with UMFPACK to **3.5e-14 or better**
at every size tried, 9,408 to 148,224 trace dofs, cold and warm. That is the
question TODO asked to answer locally and it is answered.

**Performance: unanswerable on this machine, and that is a measurement rather
than a shrug.** Three runs of the *same binary* on the *same* 9,408-dof problem,
nothing else running:

| cuDSS `SetOperator` | run 1 | run 2 | run 3 |
|---|---|---|---|
| `n = 32`, 9,408 dofs | 1.0517 s | 1.1324 s | **2.8125 s** |
| `n = 96`, 83,520 dofs | 0.1774 s | **5.7505 s** | 2.5190 s |

**A factor of thirty at fixed size.** UMFPACK and PARDISO on the same runs
reproduce to a few percent. So the spread is the device path, not the machine
being busy, and no ranking against the CPU solvers can be built on it. For scale
only, and not to be quoted: cuDSS's *best* observed warm setup at 148,224 dofs
was 3.08 s against UMFPACK's 0.65 s and PARDISO's 0.32 s, and its solve 0.199 s
against 0.133 s and 0.079 s. It did not win anything on this hardware.

**Why the machine cannot answer it**, so the next person does not re-run it
expecting better:

* **WSL2 shares the GPU with the Windows host.** `nvidia-smi` reports 7.5 GB of
  8 GB used with **no compute processes at all** — that is host-wide accounting,
  so the memory actually available to cuDSS is small and varies with whatever
  the desktop is doing. At `n = 128` the runs became openly bimodal (12.05 s and
  0.36 s on consecutive attempts).
* **FP64 at 1/32.** The standing caveat under *Traps* — a consumer card runs
  double precision at 1/32 to 1/64 of single, where a datacentre part runs it at
  about 1/2. A sparse LU is FP64 throughout.

**So the verdict is exactly the one TODO specified in advance**: correctness
locally, timing on a datacentre part before any conclusion. Nothing here argues
against cuDSS; it argues that this machine cannot be asked.

**AND ONE TRAP THAT WOULD HAVE POISONED ALL OF IT.** cuDSS queues work on a
stream and returns. Timed without a device synchronise, the warm setup of a
**148,224-unknown factorisation reads 2.0e-04 s** — a millisecond answer to a
second-scale question, off by four orders of magnitude, and *plausible enough to
publish* because "the reordering was reused, so of course it is fast" is a
story that fits. `TraceSolverScaling` calls `cudaDeviceSynchronize()` inside the
timing loop for every solver, CPU ones included, so the sync can never be the
thing that was forgotten. **Any future device timing in this project must do the
same.** The CPU solvers are synchronous and pay nothing for it.

#### The defaults, and what the fastest reachable set actually is

**What a fresh `GradShafranovSolver` does, as of 2026-08-31:**

| | default | why |
|---|---|---|
| `setTraceSolver()` | `UMFPack` | the only one present in every build, and what every rate in the suite was measured with |
| `setAssemblyMode()` | **`Serial`, unconditionally** | an automatic gate was tried and removed; see below |
| `MKL_NUM_THREADS` | **1**, set on every ctest | *Traps*. Not negotiable |

**THE ASSEMBLY DEFAULT IS FLAT, AND A GATE WAS TRIED AND REMOVED.** This is the
one place in this campaign where the isolated benchmark actively misled, so the
history is worth keeping.

MFEM's threaded path buffers the element-local blocks and scatters them
afterwards. Measured in isolation the buffering costs **0.86x** at one thread —
slower than the serial loop it imitates — and pays **1.15x to 1.33x** from four
threads up, still 1.19x on a 512-element mesh. On that evidence a gate on
`omp_get_max_threads() > 1` looks free, so one was written.

It is not free. `HighBetaConvergence` — which assembles a *small* system many
times over inside a bordered Newton — went from **21.5 s to 39 s** under it.
1.8x slower, reproducibly, on the same binary:

| `HighBetaConvergence` | run 1 | run 2 | run 3 |
|---|---|---|---|
| `OMP_NUM_THREADS=1` → Serial | 20.96 s | 21.22 s | 22.00 s |
| `OMP_NUM_THREADS=16` → Threaded | 32.37 s | 38.84 s | 39.25 s |

**And mesh size does not separate the two cases**, which is what makes it
undecidable from inside the solver: HighBeta's meshes are 128 and 512 elements,
and 512 is exactly where the isolated benchmark still showed a win. What differs
is how *often* assembly is called relative to everything else — MFEM forks a
team and buffers per call, so a caller that assembles once amortises it and a
caller that assembles hundreds of times pays it every time. The solver cannot
know which caller it has.

So `Serial` stands, and `Threaded` is an informed opt-in: **take it for a large
mesh assembled a few times, leave it for a small one assembled hundreds of
times.** Same lesson as the trace matrix's symmetry and the threaded-MKL
collapse — a property measured on the easy configuration is not a property of
the code.

**And the fastest set is NOT the default**, because the fastest set needs a
package most builds do not have. Composed from the scan at `k = 2, n = 64`,
37,248 trace dofs, one linear solve:

| | assembly | trace setup | solve | total | |
|---|---|---|---|---|---|
| UMFPack, serial assembly — the old default | 0.470 | 0.229 | 0.033 | 0.732 | |
| UMFPack, threaded assembly @8 — opt-in, big meshes only | 0.410 | 0.229 | 0.033 | 0.671 | 1.09x |
| **PARDISO, threaded assembly @8** | 0.410 | 0.158 | 0.022 | **0.590** | **1.24x** |
| PARDISO threaded, if item 0 is ever done | 0.410 | 0.075 | 0.010 | 0.494 | 1.48x |

**PARDISO beats UMFPack even single-threaded** — 1.50x on the factorisation,
1.41x on the backsolve — so row three needs no MKL threads and costs nothing but
the licence. It is not the default because `MFEM_USE_MKL_PARDISO` is off in most
builds and oneMKL's terms are not everybody's to accept; `setTraceSolver()` and
`traceSolverAvailable()` are how a caller who has it takes it.

**Read that table for its size as well as its order.** The whole ladder is
**1.48x**, and it is dwarfed by the `MKL_NUM_THREADS=1` fix that preceded it —
3x on ordinary tests, 140x on `k >= 3`. Everything here is a single-digit-percent
improvement on a problem whose large win has already been taken. It is also
**one linear solve**, not a Newton solve: a real solve repeats `ComputeH` and the
factorisation per iteration but the integrator assembly less often, and neither
post-processing nor the estimator is counted. Indicative, not a profile.

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

### PARDISO: faster, and it runs now

**SUPERSEDED BY THE THREADING SECTION ABOVE, WHICH RE-MEASURED ALL OF IT.** The
verdicts here still hold in direction — PARDISO is faster, and it runs at sizes
where Debian's MKL returned error `-3` — but **every number below was taken
while a stray `libmkl_sequential` on the link line made MKL sequential**, so
they are sequential-PARDISO-against-sequential-UMFPACK timings that were never
labelled as such. The clean re-measurement is under *The trace solve*: 1.50x on
setup sequentially, 2.83x at PARDISO's best thread count, and a reason meq
cannot currently take the latter. Kept for the packaging argument in the second
half, which is unaffected.

The open question of the previous section is closed. That section's verdict —
"2.2× where it works, and it stops working just below the sizes meq actually
runs" — was a verdict on **Debian's `intel-mkl` 2020.4.304** and on nothing else,
and it said so. Against **oneAPI MKL 2026.1** at `/opt/intel/oneapi/mkl/latest`,
with `MFEM_USE_MKL_PARDISO=YES` and MFEM's own `PardisoSolver`:

| `k` | `n` | trace dofs | UMFPACK setup | solve | PARDISO setup | solve | agreement |
|---|---|---|---|---|---|---|---|
| 2 | 16 | 2,400 | 0.0176 | 0.0011 | 0.0233 | 0.0064 | 1.4e-15 |
| 2 | 32 | 9,408 | 0.0616 | 0.0060 | 0.0216 | 0.0029 | 3.5e-15 |
| 2 | 48 | 21,024 | 0.1608 | 0.0165 | 0.0453 | 0.0071 | 9.4e-15 |
| 3 | 16 | 3,200 | 0.0297 | 0.0019 | 0.0077 | 0.0021 | 2.3e-15 |
| 3 | 32 | 12,544 | 0.0820 | 0.0099 | 0.0216 | 0.0033 | 4.4e-15 |
| 3 | 48 | 28,032 | 0.2295 | 0.0277 | 0.0611 | 0.0120 | 1.2e-14 |

**About 3.8× on the setup and 2–3× on the solve from 9k dofs up**, agreeing with
UMFPACK to round-off throughout. It is *slower* at 2,400, which is thread startup
against a problem too small to need it.

**And it runs at 12,544 and 28,032**, which is the point: those are the sizes
where the old MKL returned error `-3`. The defect was the packaging, as that
section predicted, and the prediction is now confirmed rather than assumed.
`tests/convergence/TraceSolverComparison.cpp` is the measurement, and it asserts
only the **agreement** — the timings are printed, because PARDISO is threaded and
CLAUDE.md's standing rule is that a threaded measurement on this machine is not a
measurement about the code.

**Both paths must keep shipping.** oneMKL's licence is not everybody's to accept
and `MFEM_USE_MKL_PARDISO` is off in most builds, so the comparison compiles and
passes either way — without it the PARDISO columns are skipped and UMFPACK still
runs. meq's solver still uses UMFPACK; making PARDISO selectable at run time is a
small piece of work that has not been done.

**`SetReuseSymbolic()` exists on `PardisoSolver` too**, with the same
pattern-comparison contract as `UMFPackSolver`'s, so the one unqualified win of
the previous section carries over rather than having to be re-argued.

### The link line straddles two MKL installations, and it cost 140x

**THE HAZARD THIS SECTION WARNED ABOUT WENT OFF, AND IT WENT OFF IN THE
OPPOSITE DIRECTION TO THE ONE PREDICTED.** The prediction was wrong answers from
mixed threading layers. What actually happened is that the stray
`libmkl_sequential` was *load bearing*: it made MKL resolve sequential for
everything, which is why every timing in this file before 2026-08-30 was fast
and why nobody noticed MKL threading is ruinous here. Removing it — by pointing
MFEM's PARDISO at oneAPI's `mkl_gnu_thread` instead of Debian's
`mkl_sequential`, which was the *correct* thing to do — exposed a 140x
regression in the test suite. See *Threaded MKL is a catastrophe* under *Traps*.

**Current state, `MFEM_EXT_LIBS` 2026-08-30**: one threading layer fewer than
below. `mkl_sequential` is gone; oneAPI's `mkl_gnu_thread` appears twice (BLAS
and PARDISO), and Debian's `mkl_intel_thread` still arrives behind SuiteSparse.
So it is still two versions and two threading layers, and it is still an
ordering accident that oneAPI wins. **Repointing SuiteSparse's BLAS at oneAPI is
still the fix and is still not done.** The rest of this section describes the
state before that change and is kept because the argument is unchanged.

The previous section warned that the first thing likely to go wrong in a oneMKL
adoption is "meq's existing MKL BLAS/LAPACK link line ending up straddling two
MKL versions". It does. `MFEM_EXT_LIBS` reads, in order:

```
-L/opt/intel/oneapi/mkl/latest/lib -Wl,-rpath,... -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
...
-L/usr/lib/x86_64-linux-gnu -lumfpack ... -lcholmod ... -lmkl_intel_thread -lmkl_sequential
```

oneAPI 2026.1 first, from the explicit `BLAS_LIBRARIES`/`LAPACK_LIBRARIES`; then
Debian's `intel-mkl` 2020.4.304 behind it, pulled in as SuiteSparse's own BLAS
dependency, bringing **two more threading layers**. Linking more than one MKL
threading layer is unsupported by MKL, and mixing two MKL *versions* is worse.

**It works today** — PARDISO agrees with UMFPACK to 1e-15 and every convergence
test passes — because oneAPI comes first in the link order and its rpath is
baked in, so its symbols win at load time. That is an ordering accident, not a
configuration. Repointing SuiteSparse's BLAS at oneAPI would settle it and has
not been done.

**One consequence is already measured and is a simplification**: with oneAPI
linked explicitly as `mkl_gnu_thread` rather than resolved through
`libmkl_rt`, `SolovievConvergence` passes with `MKL_THREADING_LAYER` **unset**.
The trap under *Traps* is about the `libmkl_rt` dispatcher needing to be told
which layer to use; linking a layer directly removes the dispatcher. That is one
suite passing, not a proof about a failure mode whose signature is silence, so
the environment variable stays set everywhere until more of the suite has been
run without it.

**And the honest caveat on all of it****And the honest caveat on all of it**: on a hard case the dominant cost is the
element-local *nonlinear* iteration, not the global solve. 42 outer Newton steps
each running thousands of local Newtons is where the pedestal's time goes.
Globalisation is a bigger lever than anything above.

**The transferable lesson**, which is the same one the Solov'ev coefficients
taught: a property measured on the easy configuration is not a property of the
code. Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on the
geometry meq is actually for.

## Traps

**Every run needs `MKL_THREADING_LAYER=GNU` AND `MKL_NUM_THREADS=1`.** Two
variables, two entirely unrelated reasons, and both are now set on every
registered ctest.

`MKL_THREADING_LAYER=GNU` is the **correctness** one.
`/usr/lib/x86_64-linux-gnu/libblas.so.3` on this machine resolves to
`libmkl_rt.so`, which silently corrupts UMFPACK's BLAS-3 without it. *Silently* —
you get numbers, and they are wrong. Same trap as `../mfem-hdg-dev/CLAUDE.md`
records.

`MKL_NUM_THREADS=1` is the **factor of 140** one, and it is not a micro-
optimisation. Measured 2026-08-30, one process, nothing else running:

| | `MKL_NUM_THREADS=1` | `=16` |
|---|---|---|
| `SolovievConvergence` | **1.24 s** | 177.83 s, at 273% CPU |
| `FieldConvergence` | **1.11 s** | 198.35 s |
| `SamplerConvergence` (never reaches UMFPACK) | 0.64 s | 0.63 s |

273% of *sixteen* cores is the shape of the answer: the threads are spinning on
barriers, not working.

**THE CULPRIT IS NOT UMFPACK, AND THE FIRST VERSION OF THIS SECTION SAID IT
WAS.** The plausible story — UMFPACK calling BLAS on the small dense frontal
matrices of a sparse LU, thousands of times, each threaded call costing a fork
and a barrier — fits every whole-test fact and is **wrong**. Separated out in
the scan, UMFPACK barely moves:

| `UMFPackSolver::SetOperator` | `MKL=1` | `MKL=2` | `MKL=16` |
|---|---|---|---|
| `k = 3, n = 32` | 0.0503 s | 0.0454 s | — |
| `k = 2, n = 64` | 0.2091 s | 0.2406 s | 0.2893 s |

Worse with threads, but by 40%, not by 140x.

**It is MFEM's element-local dense LU.** `ComputeH()` factors `A`, forms the
Schur complement and factors that, on blocks of order 10-30, once per element,
through `LUFactors` -> LAPACK -> MKL:

| assembly + reduction | `MKL=1` | `MKL=2` |
|---|---|---|
| `k = 2, n = 64` | 0.469 s | 0.440 s — untouched |
| `k = 3, n = 16` | 0.056 s | **2.306 s** |
| `k = 3, n = 32` | 0.216 s | **8.973 s** |

**`k = 2` does not move and `k = 3` degrades by a factor of forty**, putting
MKL's threading threshold between those two block sizes. That is also why the
whole-test factors are so much larger than anything in the solver columns: the
suite runs `k` up to 4, and every element of every mesh pays it. A production
run at `k >= 3` with MKL threading left on would be unusable and would look
like a solver problem.

**It is MKL and not OpenMP**, which the obvious cross settles:

| `FieldConvergence` | |
|---|---|
| `OMP=16 MKL=1` | **1.19 s**, 99% CPU |
| `OMP=1 MKL=16` | 183.54 s, 266% CPU |

**Keep the wrong answer in view, because of how it was wrong.** It fitted every
fact then available — the collapse tracked exactly the tests that reach a direct
solver — and no whole-test timing could have refuted it, because in those tests
UMFPACK and `ComputeH()` always appear together. Only separating the columns
did. Same lesson as the Solov'ev coefficients: checking a story against the
evidence that suggested it proves nothing.


**Why this never bit before, which is the part worth internalising.** The link
line used to carry `libmkl_sequential` behind SuiteSparse, so MKL resolved to
the sequential layer whatever it was asked for. Pointing MFEM's PARDISO at
oneAPI's threaded layer removed it and revealed the real behaviour. **The suite
was never deliberately sequential; it was accidentally so**, and every timing
recorded in this file before 2026-08-30 was taken under that accident. This is
the same species of finding as the threaded-BLAS one above: a property of the
link line masquerading as a property of the code.

**`MKL_NUM_THREADS=1` is right for meq's solver and wrong for PARDISO**, which
is a genuine tension rather than an oversight — PARDISO is a real parallel
sparse solver and scales. See *The trace solve* for the numbers and for what
that implies about ever making PARDISO meq's default.

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
`mfem_error("Number of Blocks don't match")`. It survived because the checks on
that path are `MFEM_ASSERT`, compiled out with `NDEBUG`. Found in stage 4 and
fixed by passing views over `darcy->GetOffsets()`, as
`DarcyOperator::ImplicitSolve` does — **a latent defect in the linear path, not
something Newton introduced.** Worth a debug build now and then for this reason
alone.

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

**`GradShafranovSolver` has no `setInitialGuess()`, and no caller can work around
it** — `prepare()` resets the iterate and `solve()` calls `prepare()`. Until that
exists, those four cases are posed with a non-homogeneous ramp that puts `ψ = 0`
in the interior. This is the most actionable gap in the solver.

**`DarcyHybridization` freezes the element-local Newton's initial guess at
`FormLinearSystem()` time, and `ComputeSolution()` will not take another.**
`EliminateVDofsInRHS` copies the flux and potential blocks into `darcy_u` and
`darcy_p`, and every local solve afterwards — residual, gradient and recovery
alike — starts there, however far the trace has since travelled. Seeding the
vector you hand `RecoverFEMSolution()` is inert. On an ordinary Newton path that
costs iterations; on anything that **differentiates the residual by differencing
it** it costs correctness, because a local solve that hits its 100-iteration cap
returns whatever it had reached and is a function of nothing. Measured, it moved
`max ψ_h` from 0.896 to 3.84 for a perturbation of 9e-6. The fix is to re-form
the system from the recovered state; `GradShafranovSolver::formSystem()` exists
for it. Full account under *Newton, and the obligation it creates*.

**~~`DarcyForm::Reconstruct()` returns a different function where `∂F/∂ψ`
vanishes~~ — FIXED, and kept here because of the shape of it.** It returned a
different function, silently, on any element where `∂F/∂ψ` vanished: a singular
local matrix, factored anyway, because the mean-value branch that exists to
regularise it was skipped by a flag set on the wrong condition. Per element, so a
profile with a flat segment corrupted part of the domain and left the rest exact,
and a whole-domain check missed it — which is why the driver ran on
`Potential::Raw`.

The fix is `gf-hdg-dev`'s *"The postprocessing closes on the element average,
always"*: the close is unconditional, because the local post-processing problem
is a pure Neumann one by construction and there was never anything to decide.
`ψ*` is a post-processing on **every** element now and the driver is back on it.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` is the regression.
See *Post-processing is back*.

**Two things to keep from it.** The fix is on `gf-hdg-dev` and on no other
branch, so a `meq-integration` rebuilt without it silently loses this — see the
merge recipe. And the failure mode is the transferable part: a per-element defect
that a global norm cannot see, in a quantity four of the estimator's five terms
are built on.

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

**`cd` persists between Bash tool calls.** A `cd x && ...` that assumes the repo
root will silently run somewhere else. Use absolute paths.

**`grep -c` exits 1 on zero matches**, so a background check reports failure
spuriously.

**Shell VARIABLES do not persist between Bash tool calls, only the working
directory does** — which is the trap above read the other way round, and it is
the more dangerous half. `L=...` in one call is gone by the next, so
`until grep -q '^EXIT=' "$LOG"` with an unset `LOG` becomes
`grep -q '^EXIT=' ""`, and **`grep -q` with no file argument reads stdin**. In a
detached background process stdin never reaches EOF, so the loop blocks rather
than failing, and it blocks silently for the life of the session.

**And a completion marker appended to a log is not necessarily at the start of a
line.** `echo "EXIT=$?" >> log` after a Boost test binary lands on the end of
its final ANSI reset sequence, because Boost's coloured output does not
terminate with a newline:

```
^[[0;39;49mEXIT=201
```

`grep -q '^EXIT='` then never matches, however long you wait. ctest's own output
*does* end with a newline, so a waiter watching `ctest` works and the same
waiter watching a bare test binary does not — which is exactly the kind of
difference that gets diagnosed as "the run is slow". Drop the `^`, or `echo`
a newline first.

**Both of these cost a session six idle `sleep` loops and no wrong answers**,
which is the point worth keeping: a polling loop that can never terminate is
invisible when the thing it polls for is also reported some other way. If a
waiter is used, check it actually fired.

**The four defects meq reported to MFEM are closed, and
`HDG-DEFECTS-FROM-MEQ.md` is gone from that tree's `doc/`** — retired by
*"Retire two bug reports whose findings are all fixed and covered"*. Checked
2026-08-29, one at a time, because "closed" arrived in four different ways:

| | state |
|---|---|
| `Reconstruct()` on a singular local matrix | **fixed**, and meq's test flipped from red to green on it |
| `φ_h` unreachable after a solve | **fixed** — `mfem::TransferredDatumCoefficient` exists in `extension_hdg.hpp`. meq still calls `setTransferredBoundary()`, because rebuilding `η₅` on it is meq's work and is not done |
| `ReconstructFluxAndPot()` lifting only domain integrators | **withdrawn as not a defect** — *"the lift dropping boundary faces is required"*. meq had measured it harmless, which was the right answer for the wrong reason |
| `ComputeHDGFaceEnergy()` ignoring an installed `HDGStabilization` | **not re-measured.** A code read today still shows it computing the `{h⁻¹Q}` form with no call into an installed hook. meq does not use `HDGErrorEstimator` at all — `meq::ResidualEstimator` is its own — so this costs meq nothing either way, which is why nobody has looked again |

That last row is a code read and not a measurement, and is recorded as one.
Anyone who needs `Energy` mode should measure before trusting it.

## Testing stance

**A test asserts the behaviour that is wanted, and fails until it is there.**
Never the reverse. A test that asserts a known defect passes while the defect
stands, and that makes a green suite compatible with a broken solver — which
destroys the only property the suite is for: **100% green must mean there is no
lurking defect.** So a defect gets a *failing* test naming it, not a passing test
recording it.

This was got wrong and corrected.
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` was first written to
assert the corruption and passed; it now asserts that `ψ*` is a post-processing
of `ψ_h` on every element, and **fails**, because MFEM's reconstruction skips its
mean-value regularisation where `∂F/∂ψ` vanishes. Its message says what to delete
when it goes green.

Two consequences worth stating. The suite is expected to be **red while a known
defect stands**, and that is the intended signal rather than a broken build — the
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
apps/        drivers. Only meq.cpp, and MEQ_BUILD_APP defaults ON.
tests/       unit/ (Boost.Test), convergence/ (rate assertions),
             analytic/ (closed-form solutions used by both),
             performance/ (TraceSolverScaling + scan.sh -- built, NOT a ctest,
             because every number in it is a timing)
examples/    TOML run configurations
refs/        Refs.md is tracked; the PDFs are gitignored, fetch by doi
attic/       not ported, not built, kept visible. See its README.
docs/        the LaTeX manual, which predates the port
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
see its README.

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
