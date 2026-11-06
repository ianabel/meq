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

## 2. What MEQ has, read on 2026-11-06

1. **Re-solves.** `meq.Equilibrium.solve()` (PY-3, `python/_meq.cpp:188-262`)
   writes a private restart `.nc` and builds a **new `meq::Run`** — config, mesh,
   forms, ψ_c unless `CacheFile` — then swaps it in on success. Transactional,
   but a whole setup per `Solve`. `meq::Run` is one-shot by design (`Run.hpp:136`).
   The q-driven loop (`Run.cpp:3876-3925`) is the precedent for one live
   `GradShafranovSolver` swapping a profile (`setGGPrime`) and re-solving warm,
   but it is not transactional.
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

**A stale claim found on the way.** `ROADMAP.md:437`, `CLAUDE_INVERSION.md:417`
and `Run.cpp:3921` say the q-driven loop's symbolic factorisation survives
between solves. It does not: every solve builds its trace solver (`:15450`) and
resets the count (`:11789`), and `NewtonConvergence.cpp:1572` asserts one
symbolic factorisation per solve. Corrected on `manta/session`, where it
starts to matter.

## 3. The map, and its derivative

With `p` MaNTA's nodal state, the map `G` is a chain:

```
p  --Π-->  T  (tables: knots at the nodes, values from u, slopes from q)       linear
T  ------>  x  (the converged bordered state: ψ fields and border b)           R(x; T) = 0
x  --E-->  family  --I-->  g  (extraction, then interpolation at the label points)
```

so, with `J = ∂R/∂x` the bordered Newton matrix at the converged state,

```
G'  dp  =  I' E'  J^{-1} ( −∂R/∂T · Π dp )
G'^T dg =  Π^T (∂R/∂T)^T  J^{-T} ( −E'^T I'^T dg )
```

Each factor is one piece of work:

| | what | exists? | branch |
|---|---|---|---|
| `Π` | MaNTA rows → table knots | trivial, linear | `manta/adapter` |
| `∂R/∂T` | the residual's derivative in the table entries (values, slopes) | **no**; gauge columns are the precedent | `manta/table-derivative` |
| `J^{-1}` | the bordered solve, caller RHS, at a **retained** linearisation | pieces only | `manta/linearised-solve` |
| `E'` | directional derivative of the extraction, decisions frozen | **no** | `manta/shape-derivative` |
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

### 4.1 `manta/session` — one live run, re-solved in memory (tier 0's cost)

A **fourth phase of `meq::Run`**, not a binding over the solver (the PYTHON-PLAN
rule: a front end, not a second driver): `resolve(profiles, warm)` on a Solved
run replaces the source's tables through the same refusals `makeSource` applies,
warm-starts from a **snapshot** (below), solves, and on failure restores the
snapshot it started from — transactional in memory, as PY-3 is through a file.
Mesh, spaces, forms, ψ_c and the solver object survive. Adaptivity, the retry
rung and the Picard ladder rebuild the solver; a session refuses `[adaptivity]`
(as PY-3 does) and accepts a rebuilt solver on a fallback rung.

**The snapshot** (`meq::EquilibriumSnapshot`, a value): the ψ/q/ψ̂ BlockVector,
the border values, the axis (point and element), the boundary-flux point, the
exterior coefficients, the extracted family, and the frozen extraction
decisions of §4.5. ~1–2 MB at 46,000 trace DOFs. It is what MaNTA's opaque
`Equilibrium` holds, and what `warm` hands back.

Acceptance: a `resolve` of unchanged tables is the answer in zero Newton steps
and the same bits; a `resolve` matches a fresh run on the same tables to the
convergence tolerance; a failed `resolve` leaves `dataset()` byte-identical;
cost per warm small-change `resolve` measured against PY-3's.

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

### 4.5 `manta/shape-derivative` — `E'`, with the decisions frozen

**Differentiate the algorithm that ships, not the continuum formula.** The
node positions come from roots of ψ* while `|∇ψ|` comes from `q*`, so the
continuum shape derivative is not the derivative of what MEQ computes. The
discrete one is local:

* ray root: `δρ_j = (δc − δψ*(x_j) − ∇ψ*·δx_ax) / (R q*·u_j)`;
* weights `w_j = 2π R_j s_j Δθ / (R_j |q*_j|)` through `δx_j`, `δq*(x_j)`,
  `∇q*(x_j)·δx_j` and the arc element;
* axis: `δx_ax = −(∇q_h)^{-1} δq_h(x_ax)`; `δc` from the border (§3);
* `δψ*`, `δq*` from `δx` through the post-processing (`Reconstruct`), which is
  linear in the HDG state.

**Freezing** is the same work and is also what makes tier 0 smooth: record,
per node, the element, the band/extension mask and its foot face, and per
equilibrium the axis element and the boundary-flux point; evaluate the frozen
element's polynomial (extended past its face where the node sits in a jump).
A `Solve` warm from a snapshot reuses that snapshot's decisions where they
still bracket the root and re-decides only where they do not — deterministic
given input and warm start, which is MaNTA's contract. The ray Newton runs to
round-off rather than to a tolerance, so its stopping test is not in `G`.

**The ends.** Axis limits from the axis Hessian (`axisShapeFromHessian`):
`V' → 4π² R_ax / √det H`, `⟨R⁻²⟩ → R_ax⁻²`, `⟨|∇ψ|⟩ ~ √Ψ_N → 0`,
`⟨|∇ψ|²⟩ ~ Ψ_N → 0`, `⟨|∇ρ|^n⟩` to their finite ellipse limits — served at
ρ = 0 by a fixed axis node in the label interpolant. The outer end stays a
refusal past `outerCut`; a coupled case sets MaNTA's `UpperBoundary` inside it
(§5, D1).

Acceptance: `E'` against a central difference of a **frozen** extraction to
the difference's truncation; `E'^T` by the dot-product test; frozen against
free extraction identical where no decision changes; the smoothness test
(§4.6) passing on the diverted and limited examples.

### 4.6 The smoothness measurement — first, before any derivative

Before building `∂R/∂T`, `J^{-1}` and `E'`, measure `G` as it is: perturb a
table entry at relative steps `10^{-4} … 10^{-10}`, difference the geometry, and
plot the error of the difference against the step. A smooth `G` gives the
`h` / `h²` and round-off V; a piecewise one plateaus. That decides how much of
§4.5's freezing tier 0 needs, and it is the test every later branch keeps. It
belongs on `manta/session`, as a convergence test over the session.

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

## 5. Decisions for the user

* **D1. The label and the domain.** Recommend `sqrt_psi_n` (MEQ's ρ): fixed
  under the solve, bounded `∂g/∂ρ` at the axis where `∂g/∂Ψ_N ~ 1/√Ψ_N`, and the
  abscissa MEQ's tables already take at third order. `ρ_tor` is usual for
  tokamak transport but moves with `q(ψ)` and is not one of MEQ's abscissae.
  Domain `[0, ρ_max]` with `ρ_max ≤ outerCut` — the axis served by its limits,
  the edge kept inside the extraction's range.
* **D2. Route.** Python adapter first (§4.7), or a C++ plugin.
* **D3. The input map.** Recommend species `n_s` and `T_s` as MaNTA *evolved*
  variables, so `q` supplies the slopes a Hermite table needs; knots at the
  nodes, one table per slot. A pressure carried as an aux has no slope and costs
  an order (MaNTA learnt this on the theta pinch).
* **D4. The slot list.** Needs a tokamak transport case on MaNTA's side;
  proposed start: `V'`, `⟨|∇ρ|²⟩`, `⟨|∇ρ|⟩`, `⟨R⁻²⟩`, `⟨|∇ψ|²/R²⟩`, `q`.
  SI, as MEQ is; units go to MaNTA's output strings.
* **D5. Precision and device.** Double and host only: MaNTA is double, and
  bit-determinism is a host promise (device parity is to printed digits).

## 6. The branches

All from `main`, kept off it; each merges `main` in when it needs to (no
rebases). `manta/design` holds this plan and is where the pieces meet.

| branch | content | needs |
|---|---|---|
| `manta/design` | this plan; `MANTA-COUPLING.md` marked as tier 2's; ROADMAP item 15 | — |
| `manta/session` | `Run::resolve`, `EquilibriumSnapshot`, transactional in memory; the smoothness measurement (§4.6); the stale factorisation claim corrected | — |
| `manta/table-derivative` | `∂R/∂T` applied and transposed, differenced-column tests | — |
| `manta/linearised-solve` | the linearisation solver, `solveLinearised`, the whole-bordered-Jacobian test | `session` (the snapshot) |
| `manta/shape-derivative` | frozen extraction, `E'`, `E'^T`, `I'`, axis limits | `session` |
| `manta/adapter` | `meq.Session` binding, the Python `EquilibriumModel`, `Π`, the coupled acceptance runs | all |

**Order:** `session` (and its smoothness numbers) first; `table-derivative`
in parallel, since it needs nothing; then `linearised-solve` and
`shape-derivative`; `adapter` last, with tier 0 runnable as soon as `session`
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
