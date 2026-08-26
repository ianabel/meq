# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

meq (the Maryland Equilibrium Solver) computes axisymmetric plasma equilibria by
solving the **Grad–Shafranov** equation with a hybridizable discontinuous
Galerkin (HDG) discretisation built on MFEM.

`refs/Refs.md` indexes the papers behind the numerics, with the doi to fetch each
from. The two marked ✔ there are not background reading — they *are* the method,
and `src/meq` is an implementation of them.

## Status: the solver works, the driver does not

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
  exist yet, and its `WriteOutputMFEM()` went with it, so **nothing writes files
  today**. `v0-legacy` has the original.
* `src/meq/Estimator.{hpp,cpp}` is written, compiles, is in `meq_core` and is
  under the `naming` check. `GridFunction::GetValueFacet`, which the old header
  called and which 4.9.1 does not have, is replaced by
  `traceFes->GetFaceElement(f)` + `GetFaceVDofs()` + `CalcShape()`, exploiting
  `DG_Interface_FECollection`'s `VALUE` map type — the same pattern
  `estimators_hdg.cpp` uses.
* `apps/meq.cpp` likewise. Its target is behind `MEQ_BUILD_APP`, default `OFF`.
* `Config`, `Profiles`, `Source` compile and test.

Still missing from the solver: adaptivity (stage 6). And the driver — nothing
writes files, so meq is reachable only through its test suite. The local post-processing that was stage 3 has
been dropped — see *There is no separate post-processing stage* below.

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

### The stages, and where we are

| | | |
|---|---|---|
| 0 | Git reconciliation, tag `v0-legacy` | **done** |
| 1 | Tree, CMake, `Config` / `Profiles` / `Source` | **done** |
| 2 | Linear `Δ*` on `DarcyForm`, fitted polygonal domain | **done** |
| 3 | Local post-processing `ψ*_h` | **dropped** — see below |
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
cd build && MKL_THREADING_LAYER=GNU ctest --output-on-failure
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
The MFEM meq needs is a local branch (`gf-hdg-subdomains-dev`) published on no
remote, so a hosted runner cannot obtain it and caching does not help — there is
nothing to fetch. `find_package(MFEM)` is therefore **not `REQUIRED`**: without
it, `CMakeLists.txt` builds the MFEM-free half — `Config`, `Profiles`, `Source`,
`SourceFactory` — and `tests/CMakeLists.txt`'s existing "unit not present"
guard already skips every convergence test without needing to know why.

So `.github/workflows/ci.yml` runs four unit suites, the `naming` check and a
coverage gate at 90% lines (94.8% measured). **What it does not run is every
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
a repair**: the proper fix is to evaluate `φ_h`, which needs a route out of the
extension machinery that MFEM does not currently offer.

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
`gf-hdg-subdomains-dev`, CMake-built in `../mfem/build` from sources in
`../mfem/mfem-src`, and configured with everything meq wants:

| | |
|---|---|
| `MFEM_USE_SUNDIALS` | KINSol, so `KINSolver(KIN_LINESEARCH)` is reachable without another rebuild |
| `MFEM_USE_GSLIB` | `FindPointsGSLIB`; gslib v1.0.9 built alongside at `../mfem/gslib` |
| `MFEM_USE_SUITESPARSE` | UMFPACK, the direct solver the convergence suite runs on |
| `MFEM_USE_LAPACK`, `MFEM_USE_ZLIB` | |
| `MFEM_USE_MPI = NO` | meq is serial throughout; `../mfem-hdg-dev-par` exists for parallel work |

**`../mfem-hdg-dev` is now the development tree and not what meq links.** That
split exists because the alternative was demonstrated: mid-session that tree was
switched to `gf-hdg-dev`, which has no `fem/darcy/extension_hdg.{hpp,cpp}`, and
meq's `naming` check failed on a tree meq had not touched. A library meq depends
on should not move under it while somebody is working on it. Point `MFEM_DIR` at
the dev tree deliberately when testing a fix; otherwise leave it alone.

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

**But it is unusable through Newton, and fails silently.**
`ReconstructFluxAndPot()` reads only the *linear* `M_p`, and meq's Newton path
puts the whole potential block on `Mnl_p` — so the local problem gets no
potential mass and no constraint. Measured `ψ*` of **9.9e14, 8.4e15, 3.9e14**
against 3.8e-6, 2.4e-7, 1.5e-8 for the same problem solved linearly, with `ψ_h`
agreeing to six figures either way. `postProcess()` throws rather than returning
it. This is an MFEM-side defect, listed under *Traps*.

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

### On SUNDIALS

`mfem::KINSolver` **derives from `mfem::NewtonSolver`** (`linalg/sundials.hpp`).
So code written against a `NewtonSolver&` takes either, with no abstraction layer
and no rewrite. Start on the native solver; switch to
`KINSolver(KIN_LINESEARCH)` when globalisation is wanted, which a stiff pressure
pedestal will eventually want. `SetJFNK` and `EnableAndersonAcc` come with it.

The MFEM tree currently has `MFEM_USE_SUNDIALS = NO`, so this needs a rebuild of
*that* tree first. It should be nothing worse than a flag: MFEM 4.9.1's SUNDIALS
interface is fully modernised (`sunrealtype`, `SUNContext`) and asks only for
SUNDIALS ≥ 5, and 7.5.0 is installed at `/home/ian/projects/sundials/install`.
That mismatch is precisely why the *old* pinned MFEM 4.5.1 stopped compiling —
its `sundials.hpp` still used `realtype` and `booleantype`, which SUNDIALS 7
removed. `miniapps/hdg/darcyop.hpp` already offers `SolverType::KINSol`, so the
path is exercised in that branch.

**This is now the highest-value outstanding change, not a someday.** Stage 6's
benchmark work produced three concrete globalisation failures on the GS-2
sources: the pressure pedestal does not converge at `k = 1` for `h ≥ 0.05`
(MFEM's element-local nonlinear solves fail outright, `el: N not convered in 100
iters`), and the current hole (§4.4) does not converge under plain Newton in any
configuration tried — NaN at `k = 2, n = 16`, and MFEM *aborts the process* on
`MFEM_VERIFY(IsFinite(norm))` because the tree is built without
`MFEM_USE_EXCEPTIONS`. Both papers use Anderson-accelerated **Picard**, which
leaves every local solve linear; meq's Newton puts `∂F/∂ψ = 2c₁r²/σ² = 320r²`
inside the local solve, and those local problems are indefinite on a coarse
element. `MFEM_USE_SUNDIALS = NO` is the blocker.

**Confirm with the user before rebuilding `../mfem-hdg-dev`** — it has its own
active work on `gf-hdg-subdomains-dev`.

### The suite is 11/13 against the new library, and both failures are informative

**`ExtensionConvergence::thePostProcessedPotentialConvergesAtKPlusTwo` fails,
and it is a real regression in MFEM.** `ψ*` on the extension path converges at
**1.28 at `k = 1` and 1.16 at `k = 2`** — a rate independent of `k`, which is a
fixed `O(h)` geometric error rather than a discretisation one — where the
working configuration gave 2.62/3.00 and 3.46/3.90. It is also *larger than
`ψ_h`*, by 7× at `k = 1` and **1000× at `k = 2`**, so the post-processing now
costs accuracy rather than buying it.

The solve is untouched: `ψ_h` and `q_h` agree with the old library to every
digit printed. Only the four `ψ*`-dependent quantities moved. The suspect is
§3 of `../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md` — meq reported that
`ReconstructFluxAndPot()` dropped the boundary-face `HDGExtensionIntegrator`
and *measured that drop harmless*; the fix lifts it into the local flux block,
and that is exactly the term carrying the `Γ_h` stand-off.

`AdaptiveRefinement`'s curved loop fails for the same reason and was how this
was found — `η₂`, `η₄` and `η₅` are built on `ψ*`, so it surfaced as `η`
refusing to come down, four files away from the cause. **That is why the new
test exists**: nothing put `ψ*` and a curved boundary together, because
`ExtensionConvergence` never called `postProcess()` and `EstimatorConvergence`
works on a fitted rectangle.

**`PedestalConvergence::pedestalNewtonFailsOnCoarseMeshesAtOrderOne` fails
because the failure it asserts has gone away**, which is what that test is for
— it prints its own instructions. Newton now converges at `k = 1, h = 0.05` in
42 iterations where MFEM's element-local solves used to give up at
`el: N not convered in 100 iters`; §4.4, the current hole, no longer produces
NaN or aborts the process, though it still does not converge in 60. Do not
delete that test until the cause is settled, since its prose explains the
failure in terms that would then be wrong.

## Traps

**Every run needs `MKL_THREADING_LAYER=GNU`.**
`/usr/lib/x86_64-linux-gnu/libblas.so.3` on this machine resolves to
`libmkl_rt.so`, which silently corrupts UMFPACK's BLAS-3 without it. *Silently* —
you get numbers, and they are wrong. CMake sets it on every registered ctest;
you must set it by hand when running a binary directly. Same trap as
`../mfem-hdg-dev/CLAUDE.md` records.

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

**`DarcyForm::Reconstruct()` returns ~1e15 on the nonlinear path, silently.**
`ReconstructFluxAndPot()` consults only the linear `M_p`, never `Mnl_p` — while
`ReconstructTotalFlux()` *does* consult `Mnl_p`, so the asymmetry is internal to
MFEM. `postProcess()` throws rather than pass it on. MFEM-side, and written up
as §1 of `../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md`.

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

**Four defects found in MFEM are written up, and none is fixed.**
`../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md` carries them with the line numbers
and the measurements: the `Reconstruct()` one above; `ComputeHDGFaceEnergy()`
ignoring an installed `HDGStabilization`, so `HDGErrorEstimator`'s `Energy` mode
reports a different operator than was solved while its `Residual` mode does not;
`ReconstructFluxAndPot()` lifting only domain integrators, which drops
`HDGExtensionIntegrator` and is measured harmless; and `φ_h` being unreachable
after a solve, which is what forces `setTransferredBoundary()`. Nothing in that
tree was changed — meq works around all four. Read it before assuming a
surprise on the reconstruction or extension paths is meq's.

## Testing stance

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
src/meq/     the library. Config, Profiles, Source, GradShafranov are
             ported; Estimator is not (see Status).
apps/        drivers. Only meq.cpp, and it does not build yet.
tests/       unit/ (Boost.Test), convergence/ (rate assertions),
             analytic/ (closed-form solutions used by both)
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

**This is enforced**, by `.clang-tidy`'s `readability-identifier-naming` and a
ctest named `naming`. The unported legacy files are excluded, because they
violate it wholesale; they rejoin the check as they are ported.

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
