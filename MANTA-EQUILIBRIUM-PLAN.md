# MEQ as MaNTA's equilibrium: tiers 0 and 1, and the response MEQ owes

**The goal.** MaNTA's transport solve takes its magnetic geometry from MEQ
in-process, through MaNTA's `EquilibriumModel` interface: a converged
equilibrium for each transport state (tier 0), and its **linearised response**
(tier 1), so that MaNTA's Newton carries the coupling exactly and its adjoint
can differentiate through the field. The response — `G'(p) dp` and its
transpose — is what this plan is mostly about: four derivatives MEQ does not
have today, each on its own feature branch.

**Where it comes from.** MaNTA's brief of 2026-11-06,
`MEQ-INTERFACE-BRIEF-FROM-MANTA.md` (inbox in the project root; its final
version, read at MaNTA `9f7aa12` for the C++ interface and `feat/field-python`
at `f38dc17`, `5b7e193` for the Python one — built and tested there, not yet
merged to MaNTA's `main`). What it establishes is restated in §1, so this file stands
without it. The four MEQ-side findings in §2 were read from MEQ at `c5a1ecd`.

**What it supersedes.** `MANTA-COUPLING.md` describes MaNTA's tier 2 — MEQ's
trace unknowns inside MaNTA's DAE — and stays as the account of that tier. **Tier
2 does not fit a 2-D Grad–Shafranov solve as MaNTA built it**: its field blocks
are dense `nField × nField` (17 GB each at MEQ's 46,000 trace DOFs), and nothing
in MaNTA's plans makes them sparse. Tiers 0 and 1 need no unknowns in MaNTA's
vector. `PYTHON-PLAN.md`'s "no derivatives" scope is the tier-0 half of this
plan; its PY-6 (one process, one MKL, one SUNDIALS) is a prerequisite.

**It lives on feature branches, never main** until the user says a piece lands.
§6 lists them.

## 1. What MaNTA asks

MaNTA solves 1-D transport on a **fixed** grid in a label `x` it does not
interpret, HDG of degree `k`, `k+1` Chebyshev nodes per cell. The model sees:

| call | in | out |
|---|---|---|
| `Solve(states, points, weights, t, warm)` | the transport state at every basis node — `u`, `q` (a solved slope, not a derivative of `u_h`), stored `sigma`, `aux` — and quadrature weights | an opaque `Equilibrium`, a **self-contained snapshot**; throw on failure |
| `Geometry(g, e, points)` | an `Equilibrium`, label points (physics nodes; the output grid and face τ **include both domain ends**) | `g`, `nGeometry × nPoints`; slot order is a contract with the transport case |
| `applyResponse(dg, e, states, dp, …)` | a direction `dp` shaped like `states` | `dg = G'(p) dp` at the physics nodes |
| `applyResponseTranspose(dp, e, states, dg, …)` | `dg` | `dp = G'(p)^T dg` |
| `setTolerances(rtol·1e-3, atol·1e-3)`, `resetForRun()` | | converge to round-off; drop every cache |

The contract that matters most for this plan (brief §3.2):

* **Memoised byte for byte** on the state: any other state is a new `Solve`.
* **`Geometry` and `applyResponse` are called on an equilibrium that is not the
  latest**: MaNTA keeps the Jacobian build's `Equilibrium` while residuals solve
  newer ones. A snapshot that is a handle onto one live solver returns the wrong
  answer silently.
* **`G` must be smooth** at relative perturbations of `1e-8` (MaNTA's
  differenced defaults step there) and **deterministic** given the input and the
  warm start; a reused MaNTA solver is pinned bit-for-bit to a fresh one.
* MaNTA's defaults already difference `Solve` for both the response and its
  assembled transpose. **Tier 1 runs with no MEQ derivatives at all** once tier 0
  is smooth; MEQ's own response is the speed-up: on MaNTA's native model it cuts
  a transient's equilibrium solves from 1,006–1,215 (differenced) to 221–235, and
  `FieldCoupling = "auto"` keeps steady solves reduced only for a model that
  declares `EquilibriumResponse::Model`.

### 1.1 What MaNTA evolves: RoPP's transport equations

The coupled transport case evolves **conserved** quantities, because MaNTA has no
`n_s ∂_t T_s` form on its left-hand side, and in a time-dependent equilibrium
**`V'` sits inside them**. The equations are Abel et al., *Multiscale
gyrokinetics for rotating tokamak plasmas*, Rep. Prog. Phys. 76 (2013) 116201
(RoPP), §8, in the poloidal flux `ψ` with `V' = dV/dψ` and `∂_t` at fixed `ψ`:

| RoPP | conserved variable | flux function MEQ's source needs |
|---|---|---|
| (166) particles | `Ñ_s = V' ⟨n_s⟩` | `N_s`, the Boltzmann prefactor: `N_s = ⟨n_s⟩ / ⟨exp(m_s ω² R²/2T_s − Z_s e φ₀/T_s)⟩` (172), with `φ₀` from quasineutrality (97) |
| (194) energy | `W_s = (3/2) V' ⟨n_s⟩ T_s` | `T_s = (2/3) W_s / Ñ_s` — **`V'` cancels; node-local** |
| (179) angular momentum | `L = V' J ω`, `J = Σ_s m_s ⟨n_s R²⟩` (180) | `ω = L / (V' J)` |
| (143) the field | `q(ψ)`: **fixed first**, given once per run; later evolved by (143) | `I(ψ)` from `q = V' I ⟨R⁻²⟩ / 4π²` (142) |

and `ψ(R, z)` follows **instantaneously** from the generalised Grad–Shafranov
equation (RoPP (136), MEQ's `RotatingSource`) — RoPP §7.3's closing sentence. So
the equilibrium is an instantaneous function of the transport state, which is
exactly tiers 0 and 1; nothing needs MaNTA's tier 2.

**Three consequences shape this plan.**

1. **Only `T_s` is a pointwise map of MaNTA's state.** `N_s`, `ω` and `I` are
   defined through flux-surface integrals of **MEQ's own solution** — `V'`,
   `⟨exp(…)⟩` with MEQ's `φ₀`, `⟨n_s R²⟩`, `⟨R⁻²⟩`. MEQ cannot be handed tables;
   it must find the flux functions that reproduce MaNTA's contents. That is a
   **content-driven equilibrium** (§4.8), the generalisation of MEQ's q-driven
   solve (ROADMAP item 10), which already finds `I` from `q` through `V'⟨R⁻²⟩`.
2. **The conserved contents enter only as right-hand sides.** With the flux
   functions' knots as border unknowns and "MEQ's content = MaNTA's content" as
   border rows, MaNTA's state appears in `R(x; p) = 0` linearly, in the rows'
   right-hand sides — and in `T_s`'s table. `∂R/∂p` is then trivial for the
   contents; the work moves into the rows' derivative in `ψ`, which is the
   shape derivative of §4.5 applied to content integrals.
3. **The energy equation needs the geometry's time derivative — only there.**
   The fluxes do not: in RoPP's evaluated forms (170), (196) the `n ∂_tψ` of
   (167) has cancelled against the inductive `E×B` drift. At `ω = 0` the one
   term is the compressional heating (200), in `ρ`
   `V'_ρ P_comp = −(2/3) W_s ∂_t ln V'_ρ` with `ψ_ax` held (the write-up's
   (ρ30)); with rotation, the potential-energy exchange (201)–(203) adds more. At tiers 0/1 the geometry is `G(p(t), t)`, so
   **`ġ = G'(p) ṗ + ∂_t G` — MEQ's own response along MaNTA's `ṗ`.** The response
   is therefore not only a speed-up: it is how these terms can be computed at
   all without differencing in time. MaNTA hands a model no `ṗ` and a case no
   geometry time derivative today (its time-derivative-sources spec has that
   phase unbuilt): **a request to MaNTA**, §5 D7.

## 2. What MEQ has, read on 2026-11-06

1. **Re-solves.** `meq.Equilibrium.solve()` (PY-3, `python/_meq.cpp`)
   writes a private restart `.nc` and builds a **new `meq::Run`** — config, mesh,
   forms, ψ_c unless `CacheFile` — then swaps it in on success. Transactional,
   but a whole setup per `Solve`. **`Run::resolve()` (§4.1, `manta/session`) is
   the same restart held in memory**, with no file and no second setup. The
   q-driven loop is the precedent for one live `GradShafranovSolver` swapping a
   profile (`setGGPrime`) and re-solving warm, but it is not transactional.
2. **The linear solve.** The bordered Newton (`solveWithNormalisation`,
   `GradShafranov.cpp:11626-18358`) eliminates the border — ψ_ax, ψ_bnd, the
   X-point, the current scale, DtN modes, gauge knots, free coil currents —
   against **one** trace factorisation, pushing the residual and every border
   column through one blocked `DarcyNPCSolver::ArrayMult`. The trace solver and
   its factorisation are **locals of that solve** (`:15450`): discarded on
   return, and stale anyway, since the convergence test precedes the last
   `GetGradient`. No caller-RHS entry point and **no transpose solve anywhere**;
   MFEM's `DarcyNPCSolver` has `Mult`/`ArrayMult` only. The field Jacobian and the
   border columns are each checked against differences (`NewtonConvergence.cpp:289`,
   `HighBetaConvergence.cpp:1023,1200`, `FreeBoundaryCoupling.cpp:5810`,
   `RotatingMirror.cpp:4722`); the whole bordered Jacobian has no single such test.
3. **Profiles → source.** Tables are `SplineProfile` Hermite cubics
   (`Profiles.cpp:728-829`); every abscissa map is **linear in the knot values
   and slopes**. `RotatingSource`'s `F` is linear in the density and `GGPrime`
   tables at two species and fixed `T_s, ω`; nonlinear in `T_s` (inside the
   exponent and φ₀), in `ω` (as `ω²`), and in the densities from three species
   (φ₀ by root find). The analytic Jacobian already differentiates `F` in ψ
   through `prime()`/`doublePrime()`. **No dF/d(table) exists**; the precedent is
   the midplane-gauge knot columns (`assembleGaugeColumns`,
   `GradShafranov.cpp:11175`), which are exactly dR/d(knot) for a spline table,
   threaded, with `RotatingSource::crossingDerivatives` on the source side.
4. **Geometry.** Extraction (`FluxExtraction.cpp:130-260`) roots ψ* along fixed
   rays from the located axis, weights by `R|q*|`, and averages by a periodic
   trapezoid; the family (`FluxSurfaceFamily`) is an MFEM-free **value type**, ~55 KB
   at the defaults, and `GeometryCache` keys it on the ψ vector by `memcmp`.
   ~0.14–0.45 s per family. It refuses labels outside `[innerCut, outerCut]`.
   **No shape derivative exists.** Twelve discrete decisions make ψ → averages
   piecewise; the ones that bite at `1e-8` are the ray Newton/bisection stopping
   (noise ~1e-12 relative, i.e. ~1e-4 in a derivative differenced at 1e-8), the
   best-iterate acceptance in a DG face jump (a true discontinuity), element
   location across DG faces, band-extension switching on the outer surfaces, and
   the located axis changing element.

**The trace solver's symbolic factorisation lives one solve.** Every solve
builds its trace solver (`makeTraceSolver()` inside `solve()`) and resets the
count, and `NewtonConvergence::theSymbolicAnalysisIsReusedAcrossNewtonSteps`
asserts one analysis per solve; the q-driven loop pays one reordering per map
evaluation, 37 for `q-driven.toml`'s 37 solves (→ M-352). A resolve pays
the same, plus the solver's construction (§4.1).

## 3. The map, and its derivative

With `p` MaNTA's nodal state — `Ñ_s`, `W_s`, `L`, `q` and their slopes — the
map `G` is a chain:

```
p  --Π-->  (T, C)   T: the T_s tables, (2/3) W_s/Ñ_s at the nodes        pointwise
                    C: the contents Ñ_s, L, q as row right-hand sides      linear
(T, C) -->  x       the converged bordered state: ψ fields, the labels,
                    and the knots of N_s, ω, I as border unknowns         R(x; T) = B C
x  --E-->  family  --I-->  g   (extraction, then interpolation at the label points)
```

so, with `J = ∂R/∂x` the bordered Newton matrix at the converged state (the
content rows and knot columns included),

```
G'  dp  =  I' E'  J^{-1} ( B Π_C' dp − ∂R/∂T · Π_T' dp )
G'^T dg =  Π'^T ( B^T z − (∂R/∂T)^T z ),     J^T z = E'^T I'^T dg
```

Each factor is one piece of work:

| | what | exists? | branch |
|---|---|---|---|
| `Π`, `Π'` | MaNTA rows → `T_s` knots and content right-hand sides, and its derivative | node-local; `Π_T'` a 2×2 per species per node (§4.7) | `manta/adapter` |
| content rows, knot columns | MEQ's `V'⟨n_s⟩`, `V'Jω`, `V'⟨R⁻²⟩I` against MaNTA's, with `N_s`, `ω`, `I` knots as border unknowns | **no**; the q-driven solve does `I` by an outer loop, the midplane gauge has knot columns | `manta/content-solve` (§4.8) |
| `∂R/∂T` | the residual's derivative in the table entries (values, slopes) — for `T_s`, and as the knot columns of §4.8 | **no**; gauge columns are the precedent | `manta/table-derivative` |
| `J^{-1}` | the bordered solve, caller RHS, at a **retained** linearisation | pieces only | `manta/linearised-solve` |
| `E'` | directional derivative of the extraction — of one continuous in the state (§4.5) | **no** | `manta/shape-derivative` |
| `I'` | interpolation in the label, linear given the family | trivial | `manta/shape-derivative` |
| `J^{-T}` | the transposed bordered solve | **no**, and not in MFEM | §4.4 |

`(∂R/∂T)^T`, `E'^T` and `I'^T` are sparse — each node or quadrature point
touches one element's DOFs and a few knots — so their transposes are scatters
written alongside the forward ones and checked by the dot-product test.
**`J^{-T}` is the only hard transpose**, and §4.4 defers it.

### Why the labels are in `J`

MaNTA's `x` is a label (`psi_n` or `sqrt_psi_n`), so a geometry value at fixed
`x` moves with ψ_ax and ψ_bnd as well as with the field. The border unknowns
must be in the linearisation and in `E'`: the level is
`c = ψ_ax − Ψ_N(ψ_ax − ψ_bnd)`, so `δc = (1−Ψ_N) δψ_ax + Ψ_N δψ_bnd`. The
bordered solve already carries them; the response just keeps them.

## 4. Design, piece by piece

### 4.1 `manta/session` — one live run, re-solved in memory (tier 0's cost) — built

**As built**: `int Run::resolve( ProfileTables, shared_ptr<EquilibriumSnapshot
const> warm = nullptr )` and `shared_ptr<EquilibriumSnapshot const>
Run::snapshot()`, a fourth phase of `meq::Run` and not a binding over the
solver (the PYTHON-PLAN rule: a front end, not a second driver). The tables
replace the slots they name through `buildSources()` — the factories and
refusals `configure()` uses, factored out so the two cannot drift — and the run
solves again through the whole solve phase (ladder, retry rung, support
sweeps), warm from `warm` or from its own state: the `.nc` restart's fields,
seeded in `makeSolverBody` ahead of every other start. Mesh, spaces, `ψ_c`
(adopted from the last solver) and the configuration survive.

**A fresh solver per resolve, not one surviving solver object** — a change
from the design, and deliberate: the answer is then a function of
`( tables, warm )` and of nothing the run holds, which is MaNTA's determinism
contract (a reused model pinned bit-for-bit to a fresh one), and the old
solver is what a failed resolve is restored to. The cost is the solver's
construction; a 0.1% `Omega` step on `rotating-normalised` resolves in
0.30 s against PY-3's 0.48 s (→ M-SESSION-3). **Transactional** by swapping
every member the solve phase assigns aside (`SolvedState`) and back on failure.
Refused with the run untouched: `[adaptivity]`, a box support (`R_midplane`
tables are iterated across solves), `[inverse]` (shared coil set), a Device
(D5).

**The snapshot** (`meq::EquilibriumSnapshot`, MFEM-free, a value): session
id and sequence, the tables it was solved at, `( q, ψ, ψ̂ )`, the exterior
coefficients, `ψ_ax`, `ψ_bnd` and whether it is an unknown, the located axis,
the X-point, the `I_p` scale, the free currents, the q-driven `g²`
coefficients, the Newton count, and the family at `[output]`'s options
(`snapshot()` extracts it once per state). It carries no extraction
decisions, and none is carried from one snapshot to the next: the geometry is a
function of the tables alone (§4.5). §4.3 needs from it the state, the border
values and the tables, all of which are there.

**Acceptance, met** (`RunSession`, M-SESSION-3): an unchanged resolve is the
answer in zero Newton steps and the same bits; a changed table resolves warm to
the cold solve of the same tables, 1.5e-14 of the state apart in 4 steps
against 6; a resolve is the file restart of the same tables **to the bit**; a
failed resolve leaves the restart file byte for byte as it was (mutation:
dropping the solver from `SolvedState` is caught); the same tables from the
same warm start give the same bits whatever the run held between.

### 4.2 `manta/table-derivative` — `∂R/∂T`

A source method `tableDerivatives(R, ψ, slot) → (∂F/∂g, ∂F/∂g′, ∂F/∂g″)` at a
point, from which the residual's column for each knot entry is the Hermite
basis contracted against it and scattered with `SourceIntegrator`'s own
quadrature — `assembleGaugeColumns`'s pattern. `g″` is in because `dFdPsi`
reads `doublePrime()`, and the Jacobian's consistency needs it once `J` is
differentiated again; for `∂R/∂T` itself `F` reads value and slope.

Applied, not assembled: `(∂R/∂T) δT` is one element loop for a direction, and
`(∂R/∂T)^T v` one loop accumulating into the knots. The abscissa maps
(`psi_n` reversal, `sqrt_psi_n`'s axis fit and continuation, `FluxTableProfile`'s
label dependence) are linear in the knots and differentiate with them; the
`psi` abscissa's dependence on ψ_ax, ψ_bnd is a border-column term and goes
into the label columns the solver already has.

Scope, first: `RotatingSource` (the slots MaNTA would feed: species
densities, temperatures, `Omega`, `GGPrime`), two species and more, and the
normalised wrapper's current scale. `R_midplane` tables are host-only and
iterated across solves; they are refused under a session.

Acceptance: every column against a central difference of the assembled
residual in that knot entry, as `theAnalyticColumnAgreesWithTheDifferencedOne`
does; the transpose by `⟨(∂R/∂T)δT, v⟩ = ⟨δT, (∂R/∂T)^T v⟩` to round-off.

### 4.3 `manta/linearised-solve` — `J^{-1}` at a retained linearisation

**The problem the brief makes sharp**: MaNTA asks for the response at the
Jacobian build's equilibrium while newer equilibria are being solved. MFEM's
hybridization keeps element-local factorisations as mutable state, and any
residual evaluation between `NPCGradient()` and its solves invalidates them
(CLAUDE.md, Traps). So the linearisation needs **its own** hybridized form:

* the session holds a second `GradShafranovSolver` on the same mesh and spaces,
  the **linearisation solver**, never used for Newton;
* `linearise(snapshot)` loads the snapshot's state, assembles the gradient and
  border columns once, factorises, and pins that factorisation to the
  snapshot's identity; a request on another snapshot re-linearises (one numeric
  factorisation, the cost of one Newton step);
* `solveLinearised(b) → dx`, batched (`ArrayMult` over a fixed blocking, so a
  column's bits do not depend on its batch), with the border eliminated exactly
  as the Newton step does — the same corner solve, the same columns.

This also fixes the stale linearisation for free: the Newton's last
factorisation belongs to the second-to-last iterate; this one is at the snapshot.

Acceptance: `J dx = b` to round-off by applying the assembled bordered Jacobian
(a test the solver lacks today and should have anyway); `dx` for
`b = −(∂R/∂T) δT` against a central difference of two converged solves in `δT`,
to the solve's convergence error over the step.

### 4.4 `J^{-T}` — deferred, with a route that needs none

**MaNTA's default transpose assembles `G'` column by column through
`applyResponse`** and keeps it on the `Equilibrium`. With MEQ's response that
is `(3 nVars + nAux) N(k+1)` *back-solves*, not equilibrium solves — and the
columns for rows MEQ never reads (`sigma`, unused `aux`) are zero without a
solve. At 2 species × (n, T) on 32 cells at `k = 3`, 512 read columns, batched
through `ArrayMult`: seconds, once per equilibrium. **That serves steady-state
adjoints with no transpose solve.** An exact `applyResponseTranspose` needs
`J^{-T}`: transposed element-local solves and reduce/recover with the roles
swapped in `DarcyHybridization` — an upstream request to MFEM HDG dev — plus
the transposed border corner. Only gradients through transients need it.
Requested when tier 1's forward response is measured, not before.

### 4.5 `manta/shape-derivative` — an extraction continuous in the state, then `E'`

**What makes `G` jump is measured** (→ M-364, `SessionJumps` and
`SessionSmoothness`'s decision report): **a flux-surface node crossing an
element face.** The extraction roots `ψ*` along fixed rays from the axis and
weights each node by `1/(R|q*|)`; both fields are element-wise polynomials, so
when a node's root slides across a face it lands on the neighbour's polynomial,
and its position (by the `ψ*` jump over `|∇ψ|`) and its weight (by the `q*`
jump) change discontinuously. Every jump above the noise in the scans coincides
with a node of that surface changing element at that step, or is the echo of
one at the next; in no scan did the axis element, a band face or a
band mask change, and the ray Newton's stopping test, which changes iteration
counts at dozens of nodes, makes no jump above the noise. Per crossing the
jump falls with the mesh — about `h^{3.4}` at `k = 2`, the `q*` jump — while
crossings come faster as `1/h`.

**THE FIX IS SELF-CONSISTENT, NOT FROZEN.** Nothing about one equilibrium's
extraction is carried to the next: the geometry is a function of the solved
state, and so of the tables, whatever warm start reached it. Which
continuous extraction was the user's choice (§5, D8: A). The candidates, as
they were costed against M-364:

| | `G` | band, axis, `ψ_bnd` | accuracy | cost to run | cost to build | `E'` |
|---|---|---|---|---|---|---|
| **A. A continuous reconstruction**: `ψ*` and `q*` averaged into `H¹` (nodal averaging, MFEM's `ProjectDiscCoefficient`, or a mass-matrix projection), and the existing rays run on it | **continuous**; a kink where a node crosses a face (`∇ψ̃` jumps), so a straddled difference is off by the kink, bounded and ~`1/N_θ` of the gradient jump | the band's nearest-face choice still switches between faces, so the extension must be blended across neighbouring faces; the axis as the root of the continuous `q̃` | the averaging moves a value by at most the face jumps, which are the fields' own orders (`k+2`, `k+1`); every extraction output moves at the discretisation level | one projection per state, linear and `O(dofs)`; extraction unchanged | small: two `H¹` fields after `postProcess()` and `ContourTracer( ψ̃, q̃ )`, which exists; the band blend is the larger part | the old ray-root formula on `ψ̃`, `q̃`, composed with the averaging (linear) and `Reconstruct` |
| **B. The trace `ψ̂`** | **still jumps**: `ψ̂` is per face and broken at vertices, so a level-set point jumps as it passes one; and there is no single-valued flux on the skeleton but the normal `q̂·n`, so the weight comes from `q_h` again | — | the curve between face points is interpolated, losing order | — | — | — **Not a fix; recorded so it is not re-proposed** |
| **C. Per-element level-set integrals** (the coarea formula on the broken field): `V′ = Σ_e ∫_{ψ*_e = c, x∈e} 2πR dl/|∇ψ*_e|`, every element's own arc, no choice at all — `V′(c)` is then `d/dc` of the volume `∫_{ψ*>c} 2πR dA`, which is continuous in the state | **continuous**; a kink only where an arc's end passes a mesh vertex | the band integrated per nearest-face region, whose boundaries are fixed geometry, keeps it continuous but needs arcs clipped at those boundaries; the axis potential as `max ψ*`, which is continuous | `k+1` from the weight, as now; Gauss on each smooth arc piece is not floored by the DG jump the periodic trapezoid meets (IN-1's 1.2e-9) | arcs in the ~`perimeter/h` crossed elements, comparable to the rays | the largest: a new extraction (edge roots, arc parametrisation, several arcs or a tangency in one element, band regions) on triangles, where MFEM's cut rules do not go (`CLAUDE_DECLINED.md`, *Cut quadrature*); `R(θ)`, `z(θ)` for the file still from rays | per-element: interior terms plus face-end terms; its transpose a scatter |
| **D. Accept it as discretisation error** | **jumps**, by a per-crossing amount converging at about `h^{3.4}` (`k = 2`) | unchanged | unchanged | a straddled column's median error at 1e-4 needs about five more uniform levels on TCV at `k = 2` (×700 elements) or two to three at `k = 3` (×20–60), while the straddle rate doubles each level and the worst crossing is 10–60× the median | none | — |
| **E. A smooth partition of unity across faces**: each point's value a blend of the neighbouring elements' polynomials (extended) with weights smooth in the distance to their faces | **smooth**, no kinks | as A, the band blended the same way; axis as the root of the blended `q` | the blend differs from either side by the jump; its gradient carries jump/width | 2–6 element evaluations near a face or vertex | moderate; only worth it if a kink matters, which a bounded difference error says it does not | the blend's derivative, linear |

**Two measures, and they are different.** MEQ owes `G` to about 1e-6, no
better. Every jump measured is at or below that on a mesh good to 1e-6 — the
box's 1.3e-5 is on a 768-element `k = 2` mesh that is not, and TCV's are
1e-8 and below — so as an *accuracy* the jumps are already paid for. They
matter only because MaNTA's difference divides them by its 1.5e-8 step. So
each candidate is read twice:

| | `G`'s accuracy against 1e-6 | what a 1.5e-8 difference amplifies |
|---|---|---|
| now | at the field's `k+1`; the jumps sit inside it | jumps: a straddled column off by 0.07–160× its derivative (median, M-364), worst 4–6000× |
| A | unchanged in order; values move by at most the face jumps (≤ 5e-8 on TCV), inside 1e-6 | kinks only: the straddled column off by the kink in `G'`, estimated at the `q*` jump (1e-4 relative, M-308b) over `N_θ = 128`, ~1e-6 of the derivative — **an estimate, not measured** |
| B | worse (interpolated between face points) | jumps at vertex passages; not a fix |
| C | `k+1`, as now | kinks at vertex passages, of the same estimated size as A's |
| D | unchanged | jumps, falling `h^{3.4}`; on TCV at `k = 3` still a median 0.07 and a worst 3.8 of the derivative |
| E | the blend moves values by the jump and gradients by jump/width, both inside 1e-6 | nothing: `G` smooth |

A crossing's straddle rate per MaNTA column, at its 1.5e-8 step, all 24
surfaces read and `Omega` moved, is 4e-4 to 1.3e-3 per column on these meshes:
at that rate, of 512 read columns (§4.4), a third to a half of all Jacobian
builds would carry one column wrong by O(0.1)–O(1000). The same mechanism is M-308b's step in the q-driven map; a
continuous extraction removes that too.

**Candidate A holds up** (→ M-365, measured on a prototype at
`aa5adcb`). The user asked that the continuous field use the trace `ψ̂` as
well as `ψ*` and `q*`; what is built:

* **`ψ̃` ∈ H¹ of degree `k+2`, the moment interpolant.** Vertex values are the
  mean of `ψ*` over the elements there. On each interior edge the `P_k` moments
  equal `ψ̂`'s — exactly the `k+1` edge dofs a degree-`k+2` edge has beside its
  two vertices, which is why the degree is one above `ψ*`'s: `ψ̃`'s trace then
  projects onto `ψ̂`, the HDG solve's single-valued face unknown, and on no
  lower degree can every interior edge carry all of `ψ̂`. On boundary edges
  `ψ*`'s own trace (one side, so no choice; on `Γ_h` the trace unknown is a
  pinned zero, not `ψ̂`). In each cell `ψ*`'s `P_{k−1}` moments, which is the
  interior dofs' count. Unisolvent, reproduces `P_{k+2}`, and continuous by
  construction; every choice is a mean or a moment, none a side.
* **`q̃` ∈ [H¹]² of `q*`'s degree `k+1`**, the same interpolant: vertex means,
  edge `P_{k−1}` moments of the two sides' mean trace, cell `P_{k−2}` moments.
  The weight needs a continuous flux as much as the root needs a continuous
  potential; `∇ψ̃` is not one (its normal component jumps).

Measured on M-364's cases: a 1.5e-8 central difference straddling a
bisected crossing is off by 2e-6 to 2.5e-5 of the derivative, the noise floor
there, where the shipped extraction is off by 7e-3 to 690; `G` by A and by the
shipped extraction differ by less than either differs from a finer mesh, and
A's error against that mesh matches the shipped one's within 1.4×. **So A's
kink is below what a 1.5e-8 difference can see, and E's smoothness buys
nothing measurable**: A over E. `ψ̂` against the control that takes the edges
from `ψ*`'s mean trace: the same `G` to the digits that matter.

**The production form, built** (→ M-367):

1. `src/meq/ContinuousField.{hpp,cpp}`: the moment interpolant and the two
   fields. `GradShafranovSolver::continuousPotential()` and
   `continuousFlux()` build them on first use after `postProcess()` and hold
   them until the next solve or post-processing. Vertices, then edges, then
   cells, each dof written by one entity and the vertex means summed in
   element order, so the threaded build is the serial one bit for bit.
2. `Potential::Continuous` is `ContourTracer`'s default and the only pair a
   run extracts with — `_surfaces.nc`, the q-driven loop (M-308b's step is a
   face crossing, which no longer moves its map) and the session.
   `PostProcessed` and `Raw` stay as the controls `FluxSurfaceConvergence`
   pairs.
3. **The axis and level stay the solve's located axis.** `q̃`'s root, with
   `ψ̃` there, moves the labels by the axis's own error and buys no
   smoothness on the measured cases (M-366): no axis element changes
   along any scan.
4. The band's nearest-face choice is left unblended; `SessionJumps` is the
   by-hand study that would catch it switching.
5. The regression values that moved are re-measured in M-367.
   `SessionSmoothness` is green: across the bisected crossing `G` moves by
   2.4e-14 where it moved 4.65e-7.
6. `Run::familyFrom()` is gone.
7. `E'` follows: `δψ̃` is the same interpolant of `δψ*` (`Reconstruct`,
   linear) and `δψ̂` (part of the state), so the ray-root formula applies with
   `ψ̃`, `q̃` and their exact gradients.

**Differentiate the algorithm that ships, not the continuum formula** — that
stands, whichever row is chosen. The node positions come from roots of the
potential while `|∇ψ|` comes from the flux, so the continuum shape derivative is
not the derivative of what MEQ computes; `E'` is the discrete one, built after
the choice, and the ray Newton (or the arc root) runs to round-off so that its
stopping test is not in `G`.

**The ends.** Axis limits from the axis Hessian (`axisShapeFromHessian`):
`V' → 4π² R_ax / √det H`, `⟨R⁻²⟩ → R_ax⁻²`, `⟨|∇ψ|⟩ ~ √Ψ_N → 0`,
`⟨|∇ψ|²⟩ ~ Ψ_N → 0`, `⟨|∇ρ|^n⟩` to their finite ellipse limits — served at
ρ = 0 by a fixed axis node in the label interpolant. The outer end stays a
refusal past `outerCut`; a coupled case sets MaNTA's `UpperBoundary` inside it
(§5, D1).

Acceptance: `SessionSmoothness` green — no jump above the noise — on the
continuous extraction, with no decision carried between states; `E'` against a
central difference of the same extraction to the difference's truncation, with
the step kept off a kink; `E'^T` by the dot-product test; each mutation-tested.

### 4.6 The smoothness measurement — first, before any derivative

Before building `∂R/∂T`, `J^{-1}` and `E'`, measure `G` as it is: perturb a
table entry at relative steps `10^{-4} … 10^{-10}`, difference the geometry, and
plot the error of the difference against the step. A smooth `G` gives the
`h` / `h²` and round-off V; a piecewise one plateaus. That decides whether
tier 0 needs a different extraction (§4.5), and it is the test every later
branch keeps. It belongs on `manta/session`, as a convergence test over the
session.

**Measured on `manta/session`** (`SessionSmoothness`, → M-SESSION-4): `Omega`
scaled, and one density knot moved, by `±h` through `resolve()` warm from one
snapshot, on `rotating-normalised` (768 elements) and on `fixed-d-tcv`'s MXH
boundary with the same rotating source (3993 elements, curved, the band);
`V′_ρ`, `⟨|∇ρ|²⟩`, `⟨R⁻²⟩`, `⟨|∇ψ|²/R²⟩` on the innermost, middle and
outermost surfaces, and `ψ_ax`. Three findings, each deciding something:

1. **A session must converge to round-off, not to a tolerance.** At the
   shipped `NewtonRelativeTolerance = 1e-10` the same tables from two warm
   starts differ by up to 2.8e-11 (TCV) — at MaNTA's 1e-8 step, 3e-3 of
   noise — and a perturbation below ~1e-9 takes no Newton step at all, so `G`
   is flat there and its difference is exactly zero. Below ten eps the solve
   stops at its round-off floor (M-306): the noise falls to ≤ 2.6e-14 and the
   difference at `h = 1e-8` is good to 6e-6 to 5e-5 of the derivative. The
   adapter's `setTolerances` should ask for that, whatever MaNTA passes.
2. **Away from a decision, `G` is smooth**: the central difference settles
   between `h = 1e-5` and `1e-7` to 3e-8 – 3e-7 and then grows as
   `ε/h` with `ε ≈ 1e-13`.
3. **`G` jumps, and a step that straddles a jump is wrong by the jump over the
   step.** On the rectangle the outermost surface's `⟨|∇ρ|²⟩` and `V′_ρ` jump
   by 1.3e-5 of their value within `±1e-3` of the base in both inputs (the
   difference at `h = 1e-3, 1e-4` is 12–16% off); on TCV the jumps are
   1.5e-8 of the value (the innermost surface's `⟨|∇ψ|²/R²⟩`), and the base sits within 1e-5 of one (the
   difference at `h = 1e-5` is 3.4e-3 off, at `1e-4` and `1e-6` it is not).
   At MaNTA's 1.5e-8 step a 1e-8 jump is an O(1) error in the column it
   straddles. **What makes them is measured** (M-364): each is a ray node
   crossing an element face onto the neighbour's `ψ*` and `q*` — M-308b's
   mechanism — the box's at surface 23's node 24 (elements 763 and 761), TCV's
   at surface 0's node 88 (1548 and 1547) in both inputs; the axis element, the
   band faces and the band masks do not change anywhere in the scans.
   `SessionSmoothness` asserts a jump below the noise across a crossing it
   finds, and is green on the continuous extraction (§4.5, D8; M-367).
   **Trimmed to the box**, because what it needs is a crossing it is sure to
   see: a scan that happens to straddle none passes whatever the extraction
   does. It steps `Omega` until a node changes element (the per-node elements
   the family records), bisects to a 1e-9 interval and asserts `G` continuous
   across it, beside the noise floor from two warm starts — 15–30 s. TCV's
   band and curved boundary, and the long sweep, are `SessionJumps`'.
   Fixed boundary has no support change and no X-point; the boundary-flux
   point is not measured here (tables feed only `RotatingSource`, and no
   free-boundary rotating case ships).

### 4.7 `manta/adapter` — the `EquilibriumModel`

**Python first** (recommended, D2, and the route MaNTA's final brief settles
on): MaNTA's `manta.EquilibriumModel` (`PyEquilibriumModel.hpp`, on
`feat/field-python`) has `Solve`, `Geometry` and the optional `applyResponse`,
`applyResponseTranspose`, `responseStep`, `setTolerances`, `resetForRun`, and
MEQ's module already shares a process with `_manta` (PY-6). A Python adapter
over a `meq.Session` binding (`resolve`, `snapshot`, `geometry`, `response`,
`responseTranspose`) is thin on both sides and dodges the C++ plugin's hazards
(C++23 adapter units, Eigen defines in one `.so`, MFEM PIC, `-march`). The
response arrays cross numpy once per call; against a back-solve that is noise.
A C++ plugin can follow if the GIL or marshalling ever shows in a profile.

What the Python surface fixes, from MaNTA's side:

* **Shapes.** `Solve(self, states, points, weights, t, warm)` gets the whole
  `GlobalState` dict, `(nPoints, nVars)` arrays under `"Variable"`,
  `"Derivative"`, `"Flux"`, `"Aux"`; `Geometry(self, e, points)` returns
  `(nGeometry, len(points))`. `applyResponse(e, states, dp, inputPoints, weights,
  t, geometryPoints) -> dg`; `applyResponseTranspose(…, dg, …) -> dict of dp rows`;
  both only with `response = manta.EquilibriumResponse.Model`.
* **Slots by name.** A Python transport case reads `state.geom["Vprime"]` by the
  model's slot names (`geometry = [manta.GeometrySlot(name, description,
  units), …]`), and a JAX case gets `d*_dGeometry` by `jax.grad`. **The names
  MEQ publishes are part of the contract** (D4).
* **Attaching.** `runner.set_equilibrium_model(model)`, or
  `manta.registerEquilibriumModel("MEQ", factory)` with `factory(config, grid)`,
  then `EquilibriumModel = "MEQ"` in a dict or in a config file the `manta`
  command runs (it imports the module first).
* **Snapshots are kept several at a time**: the last, the Jacobian build's,
  and perturbed ones while MaNTA differences. §4.1's ~1–2 MB snapshot is
  affordable; a snapshot holding the factorisation would not be, which is one
  more reason the linearisation lives in the session (§4.3), keyed to a
  snapshot, and not in the snapshot.
* **The template** is MaNTA's `python/Tests/theta_pinch_balance.py`: its native
  C++ model rewritten in numpy as a `manta.EquilibriumModel`, matching the C++
  one to 7e-11 at tiers 0 and 1. MEQ's adapter's first test is the same shape:
  the Python adapter against a fresh `meq` run on the same tables.

**`manta/adapter` waits on `feat/field-python` reaching MaNTA's `main`**; the
other branches do not.

`Π` lives here: the slot map from MaNTA rows to MEQ table slots is the model's
own TOML table (`[MEQEquilibrium]`), as the theta pinch's `PressureAux` is.

**`Π` is two maps** (§1.1). The temperature tables are node-local quotients,

```
T_s  = (2/3) W_s / Ñ_s,        T_s' = (2/3) ( W_s' − (3/2) T_s Ñ_s' ) / Ñ_s
```

(`V'` cancels between `W_s` and `Ñ_s`), so `Π_T'` is a 2×2 block per species per
node on `(value, slope)`, and its transpose the same blocks transposed. A node
with `Ñ_s ≤ 0` is a throw, never a clamped temperature. The contents `Ñ_s`, `L`,
`q` pass through unchanged as the right-hand sides of §4.8's rows.

### 4.8 `manta/content-solve` — the equilibrium that reproduces MaNTA's contents

**The problem.** Given MaNTA's contents on its label grid — `Ñ_s = V'⟨n_s⟩`,
`L = V'Jω`, and `q` — find `ψ` and the flux functions `N_s`, `ω`, `I` such that
the GS equation holds and MEQ's own flux-surface integrals reproduce the
contents. `T_s` is given (§4.7). This is one equilibrium problem, not a table
lookup followed by a solve.

**Two routes, and the recommendation.**

* *Outer fixed point*, as `SafetyFactorSolve` does for `q → I`: solve, extract,
  update the tables, repeat. It exists in shape, but it converges linearly
  (12 outer iterations and 43 equilibria for the q-driven example), its
  stopping test makes `G` piecewise at MaNTA's 1e-8 steps, and its response
  would need the implicit function theorem over the loop anyway. **Kept as the
  reference the monolithic solve is checked against.**
* *Monolithic, recommended*: the flux functions' knots are **border unknowns**
  and the content conditions are **border rows** of the bordered Newton — the
  midplane gauge's pattern (`R₀` at 32 knots, quadratic, M-283/M-284). Newton
  stays quadratic, `G` has no outer stopping test, and the response falls out
  of the same `J`.

**The rows.** State each content as the **amount inside a surface or a shell**,
not as a density at a node. MaNTA's own conserved unknowns are cell contents —
`∫_cell Ñ_s dx` is the number of particles of species `s` in the shell between
two faces — and MEQ computes the same number as a volume integral of `n_s(R, ψ)`
over the region between two level sets. A volume integral cut by a level set is
smoother in `ψ` than a surface average (its derivative is one surface integral,
§4.5's machinery applied once per face), and it is what is conserved. So: one
row per species per MaNTA cell for particles, one per cell for angular
momentum, and `q` at the faces for the field (fixed per run at first, D3); knots of `N_s`, `ω` on MaNTA's
faces — the knot layout MaNTA found necessary for the theta pinch (brief §7.4) —
with one slope condition per cell to close a C¹ table. **The knot count and its
closure are to be measured, not assumed** (first experiment on this branch).

**What it needs.** The knot columns are `∂R/∂(knot)` for `N_s`, `ω` and `I`
inside `RotatingSource` — `manta/table-derivative`'s work, with `ω` entering the
exponent and `φ₀` through quasineutrality. The rows' derivative in `ψ` and the
labels is the shape derivative of a cut volume integral —
`manta/shape-derivative`'s. So this branch merges both and adds the rows,
their right-hand sides and the corner.

**Cost to watch.** Each border column is one trace back-solve per Newton step;
two species and `ω` on 32 cells is ~100 knots plus `q`'s ~33, batched through
`ArrayMult`. The gauge's 32 knots are the measured precedent; this is 4× it.

## 5. Decisions for the user

* **D1. The label and the domain — decided 2026-11-06: `x = ρ = √Ψ_N`, fixed
  boundary first.** `ρ` is fixed under the solve, `∂g/∂ρ` is bounded at the axis
  where `∂g/∂Ψ_N ~ 1/√Ψ_N`, and MEQ's tables take it at third order. **Fixed
  boundary makes the normalisation a constant at the edge**: `ψ_bnd` is the
  Dirichlet value, so `ρ = 1` is `Γ` itself for all time, and the outer end is
  the boundary — no separatrix, no band. Extraction at `ρ = 1` traces `Γ`
  rather than refusing (the family's `outerCut < 1` is a free-boundary guard; a
  fixed-boundary run can lift it). **At the axis it is constant only to the
  extent `ψ_ax` is**: `ψ_ax` is an unknown of MEQ's solve (`PsiAxis` is a
  starting guess). **With `q(ρ)` and `Γ` fixed, the pressure does not set the
  scale**: in a cylinder the fixed-boundary q-driven equilibrium at the same
  `p(ψ_N)` and `q(ψ_N)` exists for a continuum of central fields (force balance
  is first order in `B_z²` with a free constant — the write-up's §7), so one
  scalar must be prescribed (Grad–Hogan). **Decided 2026-11-06: hold `ψ_ax`,
  with `q` fully fixed** (`α = 1`). The located-axis row pins `ψ` at the axis to
  the held value; `ρ`'s normalisation is then constant at both ends on a fixed
  boundary, and MaNTA's frame terms vanish by construction (D6). A floating
  scale `α` beside a held `ψ_ax` would make the bordered Jacobian singular or
  nearly so, and imply `⟨E·B⟩ ∝ α̇/α`, outside RoPP's ordering — so there is
  none. Free boundary, where `ψ_bnd` moves too, comes after.

  **The torus check, done on `manta/session` (→ M-353).** The continuum
  is there in the torus too — it is VMEC's statement that `p(s)`, `ι(s)` and
  the boundary need `PHIEDGE` beside them, and with `q(Ψ_N)` fixed
  `Φ_edge = 2π (ψ_ax − ψ_bnd) ∫q dΨ_N`, so the free scalar is `ψ_ax` — and
  **MEQ's q-driven solve today, `ψ_ax` free, is pinned on it only by the
  truncation of `g²` to a polynomial.** Measured on `SafetyFactorSolver`'s
  fixture (`SafetyFactorScale`):
  * `ψ_ax` free, `G′ − I` at the answer has exactly one near-null direction,
    and it is the scale: `δψ_ax/ψ_ax = 0.53 δc₀/c₀` (`c₀ = g²` at the edge), i.e. `g ∝ ψ_ax`;
    `σ_min/σ_max` (in an orthonormal basis of the fit's space) 9.6e-5, 4.2e-5,
    1.1e-6, ~1e-6 at degrees 2 to 5.
  * Holding `ψ_ax` 10% off the reference, a q-driven fixed point exists to a
    least-squares residual that falls with the degree of `g²` — 3.5e-3, 3.8e-3,
    4.3e-4, 1.2e-4 at degrees 2 to 5 — and does **not** fall with the mesh
    (1.85e-3, 1.83e-3, 1.85e-3 at 5% on `n = 12, 24, 48`). So neither the
    torus nor the discretisation pins the scale; the representation does.
  * `examples/q-driven.toml` converges, and to the same `ψ_ax` from
    `ToroidalFieldGuess` 1.6, 2.2 and 3.0 at every degree, only because its
    target was built from an exactly quadratic `g²`: the reference is then an
    exact fixed point at every degree ≥ 2 (degrees 2 and 4 agree to 1.6e-5 in
    `ψ_ax`). Posed on `fixed-d-tcv`'s boundary with a measured `q` (not
    polynomial in `g²`), the loop **fails**: 145 map evaluations, the line
    search stopped at 0.067, outer conditioning 3.4e-8.

  **So the decision to hold `ψ_ax` is right, and it does what it is for**: held,
  the smallest singular value of the (fit + axis row) system is 0.53 against
  2.4e-3 free at degree 2 — 220× — and a floating `α` beside it is
  near-singular as the plan says (`σ_min/σ_max` 4e-4 to 1e-6). **What it
  leaves open is the row count.** With `ψ_ax` held, the `q` rows alone are
  near-degenerate in exactly one direction (`σ_min/σ_max` 1e-3 to 3e-5 at
  degrees 2 to 5: the cylinder's free axis datum) and the axis row is what
  removes it, so in a square discretisation the axis row must *replace* one `q`
  condition, not join them: with `g` knots on MaNTA's faces and `q` rows on
  the same faces, one row too many. §4.8 must decide which — `q` imposed on
  the traced surfaces and not at the axis limit, with `g`'s axis value closed
  by the axis row, is the natural reading of the cylinder (where `B_z(0)` is
  exactly the free datum) — and measure it. The fitted loop sidesteps it in
  least squares, which is well conditioned at low degree (0.53/8.6 at degree 2)
  and less so at higher (1.6e-4 relative at degree 5). **And the q-driven
  driver itself should hold `ψ_ax`** (a `[source]` key, the answer then a
  well-posed function of it) rather than leave the scale to its degree: an
  open item for the user.
* **D2. Route.** Python adapter first (§4.7); the C++ plugin kept possible by making MEQ's build ABI-compatible with MaNTA's (C++23, a concrete `-march`, `EIGEN_USE_BLAS`), investigated on `manta/plugin-abi`.
* **D3. The input map — settled by RoPP.** MaNTA evolves the conserved `V'⟨n_s⟩`,
  `(3/2)V'⟨n_s⟩T_s` and `V'Jω` (§1.1), so MEQ takes `T_s` pointwise and the rest
  as contents through §4.8's content-driven solve. **`q` is fixed first**
  (decided 2026-11-06): one `q(ρ)` per run, from the run's configuration, so
  the `I` rows are the q-driven solve's with a constant right-hand side and `q`
  is not in `p` — no `dq` column in the response, nothing for MaNTA to hand
  over. **Later**, MaNTA evolves `q` by RoPP (143); the same rows then take their
  right-hand side from `p`, and the response gains that column with no change to
  `J`. Building the rows with `q` as a right-hand side from the start keeps
  that step small.
* **D4. The slot list.** Needs a tokamak transport case on MaNTA's side;
  proposed start: `V'`, `⟨|∇ρ|²⟩`, `⟨|∇ρ|⟩`, `⟨R⁻²⟩`, `⟨|∇ψ|²/R²⟩`, `q`.
  SI, as MEQ is; units go to MaNTA's output strings.
* **D5. Precision and device.** Double and host only: MaNTA is double, and
  bit-determinism is a host promise (device parity is to printed digits).
* **D6. MaNTA's label and RoPP's `∂_t|_ψ`.** RoPP's equations are at fixed `ψ`;
  MaNTA's grid is fixed in `ρ`. On a fixed boundary with `ψ_ax` held (D1) the
  map `ψ ↔ ρ` is constant in time and the frame terms vanish; in general they
  are in `ψ̇_ax` and, on a free boundary, `ψ̇_bnd`. The transformation of RoPP's
  equations to `ρ` with both varying, and the route from fixed to free boundary,
  is written up separately (an artifact, 2026-11-06).
* **D7. Geometry time derivatives — a request to MaNTA.** `P_comp`, `Γ_s`'s
  flux-surface motion and `P_pot` need `V_ψ`, i.e. `∂_t` of the geometry at fixed
  label (§1.1, consequence 3). At tiers 0/1 that is `G'(p) ṗ + ∂_t G`: MaNTA would
  hand the model `ṗ` (IDA has it) and let a case read time-derivative slots.
  MEQ's side is the response it builds anyway, plus `∂_t G` (zero unless coil
  currents are programmed in time).
* **D8. Which continuous extraction — decided: A, built.** `G` jumped where a
  surface node crossed an element face (M-364). The user chose A — a
  continuous `H¹` reconstruction from `ψ*`, `ψ̂` and `q*` under the existing
  rays — and it is the only extraction (`_surfaces.nc`, the q-driven loop, the
  session): it removes the jump to the noise and keeps `G`'s accuracy
  (M-365, M-367), so E is not needed. The axis stays the solve's
  (M-366). `E'` follows (§4.5).

## 6. The branches

All from `main`, kept off it; each merges `main` in when it needs to (no
rebases). `manta/design` holds this plan and is where the pieces meet.

| branch | content | needs |
|---|---|---|
| `manta/design` | this plan; `MANTA-COUPLING.md` marked as tier 2's; ROADMAP item 15 | — |
| `manta/session` | `Run::resolve`, `EquilibriumSnapshot`, transactional in memory; the smoothness measurement (§4.6); the stale factorisation claim corrected | — |
| `manta/table-derivative` | `∂R/∂T` applied and transposed, differenced-column tests | — |
| `manta/linearised-solve` | the linearisation solver, `solveLinearised`, the whole-bordered-Jacobian test | `session` (the snapshot) |
| `manta/shape-derivative` | what makes `G` jump (M-364); the continuous extraction, `ψ̃` and `q̃` (§4.5, D8, M-365 to M-367); then `E'`, `E'^T`, `I'`, axis limits | `session` |
| `manta/content-solve` | the content-driven equilibrium: `N_s`, `ω`, `I` knots as border unknowns, content rows, the fixed-point reference | `table-derivative`, `shape-derivative` |
| `manta/plugin-abi` | MEQ built C++23 with MaNTA's `-march` and `EIGEN_USE_BLAS`, so a C++ plugin can share MaNTA's process | — |
| `manta/adapter` | `meq.Session` binding, the Python `EquilibriumModel`, `Π`, the coupled acceptance runs | all |

**Order:** `session` (and its smoothness numbers) first; `table-derivative`
in parallel, since it needs nothing; then `shape-derivative` and
`linearised-solve`; then `content-solve`, which is where MaNTA's variables first
reach MEQ correctly; `adapter` last, with tier 0 runnable as soon as `session`
and the binding exist. The coupled acceptance is MaNTA's: tier 0 against
tier 1-differenced against tier 1-own on one tokamak case, the same answer to
the integration error, with the equilibrium-solve counts of brief §4.5.

## 7. Deliberately not in it

* **Tier 2** (MEQ's trace DOFs in MaNTA's DAE): dense in `nField` on MaNTA's
  side. Revisit only for a field with dynamics, or a reduced ψ.
* **The exact transposed solve `J^{-T}`** until tier 1's forward response is
  measured (§4.4).
* **Coil currents or any equilibrium parameter as an adjoint parameter**:
  MaNTA has no `dR/dparameter` hook at any tier.
* **A mirror slot list**: waits for PY-5's open-surface geometry.
* **The anisotropic operator** (`(1−Δ)/R`, ANISOTROPY-PLAN): `Δ` reads `B`, so
  the response would have to differentiate it; first versions refuse an
  anisotropic source in a session.
