# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

meq (the Maryland Equilibrium Solver) computes axisymmetric plasma equilibria by
solving the **Grad–Shafranov** equation with a hybridizable discontinuous
Galerkin (HDG) discretisation built on MFEM.

`refs/Refs.md` indexes the papers behind the numerics, with the doi to fetch each
from. The two marked ✔ there are not background reading — they *are* the method,
and `src/meq` is an implementation of them.

## Status: the linear operator works, the driver does not

**meq solves a Grad–Shafranov equation**, which before stage 2 it never had.
The linear `Δ*` is on `DarcyForm` and reproduces an exact Solov'ev equilibrium at
`k+1` in both `ψ` and `q` for `k = 1,2,3` over four dyadic meshes;
`tests/convergence/SolovievConvergence.cpp` is the acceptance criterion and
prints the table. Concretely:

* `src/meq/GradShafranov.{hpp,cpp}` is ported — `DarcyForm` with
  `EnableHybridization` — is in `meq_core`, and is covered by the `naming` check
  like everything else.
* `src/meq/Solution.hpp` is **gone**, not ported. It wrapped hand-rolled block
  offsets and a second set of finite element spaces around what `DarcyForm` now
  owns, and `GradShafranovSolver` carries what was left of its job. Its
  `Prolong()` and `Update()` are stage-6 work against a mesh update that does not
  exist yet, and its `WriteOutputMFEM()` went with it, so **nothing writes files
  today**. `v0-legacy` has the original.
* `src/meq/Estimator.hpp` still targets the Waterloo `HDGBilinearForm` API from
  MFEM 4.5.1 and **does not compile**. It is stage-6 material and it also calls
  `GridFunction::GetValueFacet`, which does not exist in 4.9.1.
* `apps/meq.cpp` likewise. Its target is behind `MEQ_BUILD_APP`, default `OFF`.
* `Config`, `Profiles`, `Source` compile and test.

Still missing from the solver: Newton (stage 4), the curved boundary
(stage 5), adaptivity (stage 6). The local post-processing that was stage 3 has
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
| 4 | Newton on the semi-linear source | next |
| 5 | Curved `Γ` by extension from subdomains | |
| 6 | Adaptivity: the residual estimator and mesh update | |

Each stage ends at a **measured convergence rate**, not at "it runs". See
*Testing stance* below for why that is the acceptance criterion.

## Commands

```sh
cmake -B build -DMFEM_DIR=/home/ian/projects/mfem-hdg-dev
cmake --build build -j4
cd build && MKL_THREADING_LAYER=GNU ctest --output-on-failure
```

`MFEM_DIR` defaults to `../mfem-hdg-dev` and also reads the environment, so the
`-D` is usually unnecessary.

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

`../mfem-hdg-dev`, MFEM **4.9.1**, branch `gf-hdg-subdomains-dev`. It is an
in-source build: `libmfem.a` sits at the tree root and headers are rooted at the
same directory.

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

### There is no separate post-processing stage, and there should not be

Decided 2026-08-24. `ψ*_h`, the local post-processing that gains an order in the
scalar, is **not needed for accuracy**, and meq does not implement one.

**The papers agree.** `refs/HDG-GradShafranov.pdf` §3.2 describes the
post-processing and then says they did not implement it: the physical quantity is
`B ∝ ∇ψ`, which the mixed formulation already delivers as `q` at `k+1`. Adding an
order to `ψ` buys nothing a magnetic-confinement calculation uses. Stage 2's table
confirms `q` at `k+1` in practice, so there is nothing to recover.

**And where it *is* needed, the library already has it.**
`DarcyForm::ReconstructFluxAndPot()` builds its potential space as
`p_coll->Clone(p_coll->GetOrder() + 1)` — the `P_(k+1)` post-processing space,
off the hybridized solution, natively. `Reconstruct()` wraps it together with
`ReconstructTotalFlux()`. It requires a flux mass form that assembles, which meq
has. So if `ψ*` is ever wanted, call that; do not port the old
`GSSolver::Postprocess()`, which was written against a different flux convention
and whose `−(∇w, r q)` sign would need re-measuring against `DarcyForm`'s `−q`.

**The one place it will actually be needed is stage 6.** The residual estimator of
`refs/HDG-GradShafranov-Adaptive.pdf` eq. (20) uses `ψ*_h`, not `ψ_h`, in `η₁`,
`η₂` and `η₅` — and that is not decoration. The paper is explicit that `η₂` built
on the raw `ψ_h` converges at *reduced* order, because it differentiates the
approximation; substituting `ψ*_h` is what preserves `k+1`. meq's
pre-modernisation estimator used raw `ψ_h` in both places and was a degraded copy
of the published one.

So the sequencing is: nothing now, and when adaptivity arrives, **measure whether
`Reconstruct()` delivers `k+2`** before writing any local solve by hand.

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

**Confirm with the user before rebuilding `../mfem-hdg-dev`** — it has its own
active work on `gf-hdg-subdomains-dev`.

## Traps

**Every run needs `MKL_THREADING_LAYER=GNU`.**
`/usr/lib/x86_64-linux-gnu/libblas.so.3` on this machine resolves to
`libmkl_rt.so`, which silently corrupts UMFPACK's BLAS-3 without it. *Silently* —
you get numbers, and they are wrong. CMake sets it on every registered ctest;
you must set it by hand when running a binary directly. Same trap as
`../mfem-hdg-dev/CLAUDE.md` records.

**If you rebuild MFEM: `make -j4`, never more** — larger makes exhaust memory on
this machine — and **`make clean` after editing any MFEM header.** MFEM's
makefiles have no `.d` files and no header dependency tracking. That trap has
produced heap corruption in unrelated functions and "unimplemented" aborts for
methods that had just been added. meq's own CMake build tracks headers properly;
this applies only to the MFEM tree.

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
2. **An exact solution**, as an absolute-error regression.
3. **Convergence rates**, in `tests/convergence/` — asserting a *rate*, not a
   tolerance: `k+1` for `ψ` and `∇ψ`, `k+2` for `ψ*` at `k ≥ 1`, over `k = 1…4`
   and several dyadic refinements, against Tables 1–5 of the first paper and
   Table 1 of the second.

Expect rates to stop improving beyond `k ≈ 5` or `h/8`: the papers report
round-off dominating there, so a table that flattens at that point is behaving
correctly and is not a bug to chase.

**Mutation-test a suite you are relying on.** `ProfilesTests` and `SourceTests`
were checked by deliberately introducing fifteen defects — a dropped `r²`, `p'`
replaced by `p`, `μ₀` on the wrong term, the Solov'ev sign flipped, an
off-by-one in the interval lookup at a knot, `write()` truncated to six
significant figures — and confirming each was caught. That is a cheap way to
find out whether a green suite means anything, and worth repeating for the
convergence tests, where the risk of a test that passes regardless is higher.

Analytic solutions live in `tests/analytic/`. `ManufacturedNonlinear.hpp` is
`refs/HDG-GradShafranov.pdf` Example 5, `ψ = sin(k_r(r+r₀))cos(k_z z)` with a
source carrying `ψ²` and `e^{−ψ}` terms — the case that exercises Newton, and
whose `∂F/∂ψ` is analytic.

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
