# CLAUDE_DECLINED.md

**WHAT MEQ CONSIDERED AND DECIDED AGAINST, IN ONE PLACE.** Every list of
declined, falsified, reverted, withdrawn or out-of-scope ideas from the plan
files and the CLAUDE companions is gathered here. Each entry names the file it
came from, and that file carries a pointer back at the same spot:
`→ CLAUDE_DECLINED.md, *<entry title>*`. **Entry titles are the stable
references, so do not rename them.**

**The rule this file serves is `CLAUDE.md`'s: a hypothesis that was MEASURED
AND FALSIFIED STAYS.** Recording it is what stops it being re-derived, and the
measurement is the part worth keeping, so every entry keeps its numbers and its
`M-nn` anchors. What does NOT belong here is history of the files. An idea that
was later built is not declined; delete its entry.

**Read the verdict before the text.**

| verdict | meaning |
|---|---|
| **Falsified** | measured, and the claim does not hold |
| **Reverted** | built, measured, and taken out |
| **Declined** | not built, for a stated reason |
| **Withdrawn** | a claim or a request taken back |
| **Closed** | answered: not worth doing at the measured price |
| **Deferred** | not now, with the condition that would reopen it |
| **Out of scope** | not MEQ's problem as posed |

**Before reopening an entry, read its measurement and say what has changed
since it was taken.** A new library, a new mesh or a new case can move a cost;
it cannot move a falsification.

Sections follow the campaigns. Lists that are interleaved with live material,
or that code cites by section number, stay where they are and are indexed under
*Recorded in place* at the end.

## Contents

* **Newton, globalisation and solver selection**
  * *Border scaling: a per-border ruler for the bordered merit*
  * *A backtracking line search on the full residual*
  * *Predictive solver selection: the killed detectors*
  * *`PicardThenNewton` as the default, and `CondenseThenLinearise` as a cheaper rung*
  * *Fixing `ψ_ax`, or closing it outside the solver*
  * *Why MEQ's Newton failed: four wrong answers*
  * *Continuation in the driver*
  * *Methods that need Jᵀ: merit-gradient Armijo, dogleg, trust region, JFNK, full-space SQP*
  * *`Normalisation::Decoupled`: a Picard on the border round a Newton on the field*
  * *KINSOL and Anderson acceleration on the bordered fixed point*
  * *The bordered ladder's rejected-outright list*
* **Linear solves, threading and the build**
  * *UMFPACK's threaded BLAS as the `MKL_NUM_THREADS` culprit*
  * *MKL_THREADING_LAYER=GNU*
  * *Symmetric solvers, flux-block Cholesky, AMG, and re-timing cuDSS*
  * *An iterative trace solve: GMRES, and Ginkgo*
  * *MFEM as a submodule*
  * *PARDISO's old failure: the sequential timings and three dead hypotheses*
  * *Threading the output writers (E) and the Newton loop's remainder (F)*
* **Device and single precision**
  * *cuDSS determinism: SetDeterministic and per-column Mult*
  * *A per-pair rate criterion, or a wider `rateSlack`, in single precision*
  * *A GPU performance verdict from this machine*
  * *`TraceSolver = "cudss"` from a configuration file alone*
  * *Partial assembly*
  * *`AssemblyMode::Batched` as a host default*
  * *A fully device-resident solve as a commitment*
  * *Host threading as a device-plan item*
  * *Digit-for-digit agreement as the device instrument*
  * *Items 6 and 8 strictly last, as uniquely destructive*
  * *A missing device route as the gate*
  * *Device findings MEQ does not file upstream*
* **Free boundary: borders, limiter, X-point and the plasma support**
  * *Confinement, clamping and a prescribed current: the M-134 two-by-two*
  * *Conditioning, units and combinatorics as the cause of the bordered stall*
  * *Starting from vacuum, and amplitude continuation*
  * *PicardSweeps and BorderedPicardThenNewton*
  * *Globalisation, and a transferable merit weight, as the repair for the cold failures*
  * *Changing the fill, or blocking axis-touching elements, to keep the axis out of the support*
  * *Corner-only refinement for the ψ_ax layer*
  * *Narrowing ψ_ax's nodal maximum instead of locating the axis*
  * *Initialising ψ_bnd from its own border row*
  * *Auxiliary globally-coupled unknowns*
* **Free boundary: geometry, conductors, meshing and the benchmarks**
  * *The Lackner / von Hagenow Green's-function route*
  * *A directly coupled BEM, and an outer fixed point between interior and exterior*
  * *Lost coverage as the cause of the red tiling, and relaxing the gate*
  * *A closed-form conductor field as more accurate than one assembled on the mesh*
  * *Triangle, CGAL Mesh_2 and MMG as the mesh generator*
  * *Cold-start guesses: the vacuum field alone, a filament column, a rectangular column*
  * *η₆ in the marking*
* **Plasma edge**
  * *Cut quadrature, Algoim and moment fitting*
  * *Locating the plasma edge with `ψ*` rather than `ψ_h`*
  * *Fitting the mesh to `Γ_p`, or enriching the space (XFEM/GFEM)*
  * *The interior-interface plan's scope: `j` stays in the physics, and the fixed-boundary path gains nothing*
  * *A single global Fourier modal border on a cornered `Γ_p`*
* **Inverse solve**
  * *An outer loop round the forward solve for the coil currents*
  * *Machinery right only for reachable targets: one least squares over every row, row four with the iterate's `H`, and the augmented-norm merit*
* **Coil subtraction**
  * *CS-5's acceptance against `freegs4e`*
  * *A total `GridFunction` `psi_p + psi_c` for consumers to read*
* **Solution inversion**
  * *Open-surface averages and an X-point-aware tracer*
  * *Implicit quadrature on the level set as IN-2's third cross-check*
  * *Other ways to extract a contour: marching, subdivision, rational contours, unprojected ODEs, rays, Morse–Smale*
  * *A kernel survey: the Kintner recurrence, Vandermonde + GEMM, and batched field evaluation*
* **Toroidal flow and the mirror**
  * *XP-1's element flood fill on the mirror path*
  * *A plasma-current border for the mirror*
* **Interpolatory HDG**
  * *The interpolatory source: seven hypotheses rejected*
* **Output and sampling**
  * *Smoothing the boundary displacement on the vdof graph, and a global backoff*
  * *Rebuilding the ParaView collection per frame*
  * *Pinning ∇q in the band*
  * *Face averaging in the grid sampler: gslib's SetL2AvgType*
* **Out of scope**
  * *Iron and the evolution problem*
  * *What the X-point pathway does not reach: the plasma edge, inversion through a separatrix, a fixed boundary on one*
  * *The inversion's out-of-scope list: topology, 3D, straight-field-line angles, `p` contours and the inverse-equilibrium formulation*
  * *Anisotropic pressure, NetCDF profiles and MaNTA coupling*
* **Recorded in place**

## Newton, globalisation and solver selection

### Border scaling: a per-border ruler for the bordered merit

**Reverted, withdrawn and closed; one stage landed.** The plan held that
`augmentedNorm()` weights six structurally different constraints — `ψ_ax`,
`ψ_bnd`, the X-point, `μ₀ I_p` and the `N` exterior modes — with **one** scale
`gamma = ‖∂R/∂ψ_ax‖` taken from the first of them, plus a hand-derived,
mesh-dependent `xScale = |R_X|·h` for the X-point, and proposed a per-border
scale from each border's own column and corner. What came of each stage:

| | stage | outcome |
|---|---|---|
| BS-0 | print the border decomposition and classify the failures | **done**, [M-117](MEASUREMENTS.md#m-117). `directionFinite` is true on every step of all six exit-2 cases, so a broken elimination (H2) is excluded; the border's share of the merit is three orders of magnitude apart between MAST and TCV, so no single `gamma` is right for both. Machine B is an X-point excursion and not "one fatal step" |
| BS-1 | `gamma_i` = the column's norm over its corner diagonal `D_ii`, for `ψ_ax`, `ψ_bnd`, the current and the modes | **built, reverted**, [M-245](MEASUREMENTS.md#m-245). The scales differ as predicted — on DIII-D the current is overpriced tenfold and the modes underpriced up to twelvefold — and the ruler closes none of the six, loses MAST's correct solve, turns MAST shaped's failure into a wrong answer reported as converged, and moves TCV shaped to a neighbouring equilibrium at a hundred times the work. A control computing the same scales with the merit left on the shared `gamma` reproduces the baseline to every digit, so the ruler alone did it |
| BS-2 | choose the X-point's scale: `R·h` (X1) against the implied displacement `gamma_ax·‖(∇q)⁻¹q‖/L` (X3) | **declined.** It is a reweighting of the same merit, so it is a branch control with BS-1's hazard; a better-founded X-point scale is not thereby a neutral one |
| BS-3 | the four converging M-103 cases reach the same equilibrium | **closed**: a regression with nothing left to regress. It is what said BS-1 moves answers |
| BS-4 | re-take M-103 under the new ruler | **subsumed**: [M-244](MEASUREMENTS.md#m-244) is the re-take and M-245 the BS-1 arm against it |
| BS-5 | cap how far a line-search trial may move the X-point from its seed | **landed on the retry rung only**, [M-247](MEASUREMENTS.md#m-247), [M-248](MEASUREMENTS.md#m-248): `setXPointRelocationCap()` at `0.05 × [boundary.exterior] Radius` closes A, A shaped, E and E shaped and moves nothing. From the first solve it loses MAST and MAST-U, so it is legal only as a reactive rung |
| BS-6 | (a) does `checkAxisSource()` ask the plasma component? | **answered**, M-247: it does, and G shaped's refusal is right |

**Reading `gamma_ax` off iteration 0's assembled column**, which would save two
`fieldResidual()` calls a solve and make the yardstick consistent with the
Jacobian, is **closed**: it moves `HighBetaConvergence`'s history for a saving
too small to see, and with BS-1 reverted nothing else rests on it.

**The affine-covariant merit is declined.** Deuflhard's monotonicity test on the
simplified-Newton correction `‖J(x_k)⁻¹R(x_k + αδx)‖` is invariant under any
row scaling, needs no `Jᵀ`, and costs one more backsolve per trial against a
factorisation MEQ already holds. It is declined for three reasons: up to nine
extra backsolves a step on the logs' trial counts; it changes the merit on
**every** case including those that converge, which is BS-1's hazard at full
strength; and it cannot help B, which is an X-point excursion and not a merit
that mis-prices a good step.

**Continuation in `μ₀ I_p` (BG-7) is falsified.** A ramp over
`λ ∈ { 1.00, 0.90, 0.75, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05 }`: machine B
fails identically at all ten, *"no finite damped step"*, from
`I_p = 2.5e+05 A` to `1.25e+04 A`. Machine E fails the first six; at `λ = 0.30`
and `0.20` it converges in 13 and 12 iterations to `ψ_ax = 7.4e-04` and
`8.3e-04`, `ψ_bnd = −6.87e-02` and `−6.68e-02` — refused by the topology guard,
axis flux absent, X drift 1.093 m, a degenerate branch rather than the
equilibrium. The route's own falsifier, *if `θ` must fall below 1e-2 to take the
first step there is a limit point*, is met in the strongest form: at
`λ = 0.05` there is no first step. The failure does not depend on the
amplitude, so it is neither an eigenvalue balance nor stiffness. Ramping `I_p`
remains *legitimate* (`CLAUDE_FB.md`); it is measured useless on the two
machines it was aimed at.

**The lesson that survives the plan**: on these problems the merit decides
*which* step is accepted, and M-26's three routes reach solutions 9.4% apart,
so any reweighting of it is a branch control and must be judged as one — by
whether the converging cases keep their answers, not by how many failing cases
it closes.

*From `BORDER-SCALING-PLAN.md`, deleted; its status, §1.4, §3.1, §4, §5, §6 and
§7.*

### A backtracking line search on the full residual

**Reverted.** Built in MEQ, measured, made every case worse, and not kept.

**A LINE SEARCH ON THE FULL RESIDUAL DOES NOT FIX IT, AND [M-22](MEASUREMENTS.md#m-22) SAYS
WHY.** Upstream recommends backtracking on the full residual as NPC's
globalisation — `NewtonSolver::ComputeScalingFactor`, a dozen lines,
`miniapps/hdg/navierstokes.cpp` carries one. **It was implemented in MEQ and
measured, and it made every case worse**, including the five that converge
undamped: all of them went to sixty iterations, creeping by about 1% a step.

**THE MECHANISM OF THE CREEP IS A DEFECT IN THE REFERENCE IMPLEMENTATION, AND
UPSTREAM FOUND IT AFTER MEQ REPORTED THE SYMPTOM.** MEQ's first reading — that
`α` collapses because the `ℓ²` merit charges for restoring `(1 − α)` of the flux
and trace residuals — describes the pressure correctly and gets the mechanism
wrong. `NSBacktrackingNewton` accepts on `Norm(rt) < n0`, **a monotone test with
no sufficient-decrease constant**. For Newton on an `ℓ²` merit the direction is
always a descent direction, so `merit(α) ≈ merit(0)(1 − α)` and *any* small
enough `α` passes: upstream swept it and `α = 1.2e-4` still "succeeds",
improving the merit by 1e-4 relative. **So where `α = 1` is rejected the search
does not fail — it creeps**, which is exactly the 1%-a-step crawl MEQ measured.
An Armijo test would reject those steps and fail honestly. Fixing it is open
upstream, and **it would not make MEQ's problem converge**; it would stop the
globalisation disguising its own failure. **KINSOL's `KIN_LINESEARCH` does use
an Armijo condition and fails on the same cases**, which is why the two rows of
[M-28](MEASUREMENTS.md#m-28) agree.

**AND THE BLOCK STRUCTURE IS NOT WHAT SEPARATES MEQ'S PROBLEM FROM UPSTREAM'S —
SEVERITY IS.** The hypothesis was that the failure follows from two rows being linear
and one non-linear. Upstream measured the same shape on their own Darcy pedestal
— a full step takes the flux row to 4e-17 and the trace row to 3e-14 with the
whole remainder in the potential row — and there the line search *works*: 3 of 4
configurations fail undamped and converge with backtracking. **What differs is
that a full step IMPROVES their potential residual, 4.7e-01 → 1.7e-01, and
MULTIPLIES MEQ's BY 77.** The block structure is common to both; the size of the
first correction is not. So the recommendation in their §6 is supported there
and wrong here, and MEQ's original explanation was over-general.

MEQ **asked whether §6's baseline was the deleted mode, and it was not** — it
was plain undamped `NewtonSolver` on NPC, so the disagreement is real rather
than an artefact of what was being compared against. On the strength of MEQ's
measurements upstream has withdrawn two §6 claims: *"NPC is not automatically
faster … not fewer floating-point operations"*, and *"reach for [the reduced
trace operator] unless you have a reason not to"* as a general default. §6 now
says the choice turns on what the element-local non-linear solve costs on your
problem.

So the line search was **not kept**, and the comment at the `Globalisation::None`
switch in `GradShafranov.cpp` records why. What a step-length rule would need
here is to damp the non-linear block without un-doing the exact annihilation of
the linear ones — different, and unattempted.

*From `CLAUDE_HDGGS.md`, “Why it fails, measured — and why a line search does not fix it”.*

### Predictive solver selection: the killed detectors

**Falsified.**

**ASKED AND ANSWERED BY MEASUREMENT.** The reactive ladder pays for
a failure by running Newton to its cap before handing off. The obvious
improvement is to notice earlier — and the obvious signal, the first Newton
step making the residual worse, **is not a signal at all**. NPC, plain Newton,
cap 60, cold start:

→ **[M-24](MEASUREMENTS.md#m-24)** — case · out · its · `‖r₀‖` · `‖r₁‖/‖r₀‖`

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
stabilisation, and `(F/R, w)` is a small correction to both. `‖r₀‖` is a
measurement of the boundary data and the mesh, not of the difficulty.

**So the ladder stays reactive and stays triggered by observed failure**, which
is what `CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles* concludes from `max|∂F/∂ψ|/λ₁`.
Two independent candidate detectors have been measured and both are
useless; the difference is that this one is *anti*-correlated, which is worse
than uninformative. **Do not add a predictive trigger without a measurement
that separates these fourteen rows**, and note that the failure is cheap to
observe anyway — see the 1.6 s against 53.5 s of [M-21](MEASUREMENTS.md#m-21).

*From `CLAUDE_HDGGS.md`, “There is no cheap discriminator, and the obvious one is ANTI-correlated”.*

**And there is no predictive trigger to be had.** Two candidates are measured
and dead: `max|∂F/∂ψ|/λ₁` (`CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles*; two
points, no threshold, and it needs the range of `ψ` which is not known before
solving) and `‖r₁‖/‖r₀‖`, which [M-24](MEASUREMENTS.md#m-24) measures **anti**-correlated. **A third
observation, from this tree's own shipped outputs, says the iteration count
carries no signal either**: `machine-a-testtokamak.nc` and
`machine-a-testtokamak-shaped.nc` both record `newton_iterations = 2` and no
`axis_normalised_flux` at all, while `machine-f-diiid.nc` — the one that is
right — also records `newton_iterations = 2` with
`axis_normalised_flux = 1.00000000002119`. The wrong root and the right one
take the same number of steps. **That reading is an observation and not a
measurement** — those three files carry `meq_version` 969b626, 969b626 and
d4a203b, so they are not one run of one tree, and the count is the last support
sweep's rather than the solve's. BG-0's instrumentation ([M-117](MEASUREMENTS.md#m-117)) is where it
would be re-taken. Do not propose a fourth detector without
saying how it would be falsified on all fourteen of M-103's rows.

*From `BORDERED-GLOBALISATION-PLAN.md` §1.4.*

### `PicardThenNewton` as the default, and `CondenseThenLinearise` as a cheaper rung

**Declined.** It costs more, covers less, and changes the answer; and the other ordering is not a cheaper rung either.

Plain Newton against `PicardThenNewton`, both under
NPC, Picard capped at the same 60 iterations, wall clock on an otherwise idle
machine:

→ **[M-25](MEASUREMENTS.md#m-25)** — case · plain · `PicardThenNewton` · slowdown · rel. `max ψ_h`

**Three reasons, in increasing order of how decisive they are.**

**1. It costs 1.4× to 3.1× on everything that already works**, which is real
but is the reason you would trade away for robustness if the other two did not
exist. Note also that this is with Picard *capped at 60*; the cap is a tuning
parameter and `CLAUDE_HDGGS.md`, *Picard, then Newton*, measures the handoff to be
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

→ **[M-26](MEASUREMENTS.md#m-26)** — route · Newton steps · `max ψ_h`

**A spread of 9.4%.** This is the multiple-solution finding recorded under
`CLAUDE_HDGGS.md`, *The suite is green* — coarse discretisations of these sources
carry more than one solution — meeting the solver from a third direction. It
means changing the default globalisation is **not** a performance decision:
it changes the equilibrium MEQ reports on an under-resolved mesh, which is
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

→ **[M-27](MEASUREMENTS.md#m-27)** — case · NPC plain · `CondenseThenLinearise` · NPC `PicardThenNewton`

**`PicardThenNewton` is between 2.0× and 2.6× cheaper than the ordering switch
on every case either of them cures**, and on the cases neither cures the
ordering switch costs 25× to 30× more to find that out. The ladder MEQ already
has is the right one and nothing needs adding to it.

The two cures also do **not** agree at these resolutions — 2.5e-13 on the
pedestal, but 1.0e-03 on §4.5 at `n = 30` and 9.4e-02 at `k = 2, n = 16` — which
is reason 3 again, and a second argument against reaching for whichever cure is
convenient.

*From `CLAUDE_HDGGS.md`, “Should `PicardThenNewton` simply be the default? No, and cost is the weakest of the three reasons”.*

### Fixing `ψ_ax`, or closing it outside the solver

**Falsified.** Three routes tried, in order, before `ψ_ax` went inside the residual as a bordered Newton (`CLAUDE_HDGGS.md`, *What works*).

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

*From `CLAUDE_HDGGS.md`, “What did not work, in the order it was tried”.*

### Why MEQ's Newton failed: four wrong answers

**Falsified, three times, and the fourth is the finding.** MEQ's Newton was
thought to fail on the stiff GS-2 sources, and the question ran through four
answers.

* **Not Newton, and not the boundary condition.** `refs/MFEM-GS-Newton.pdf`
  reports Newton robust where "conventional Picard-based solvers fail to
  converge", and calls its own *fixed*-boundary solver "significantly easier"
  than the free one.
* **Not the hybridization ordering.** This one was believed and acted on.
  `NLOrdering::LineariseThenCondense` removes every element-local nonlinear
  solve (`GetNumLocalNLIterations()` reads 0) and the stiff cases get *worse*.
  **Confirmed twice**: that mode is deleted, and `NPC` — the canonical version
  of the same idea, and MEQ's default — reproduces the verdict.
* **Mostly it was under-resolution.** §4.2, §4.3 and §4.5 are not stiff at all:
  raw Newton solves all three in 7 to 17 iterations once resolved, and refining
  `h` or raising `k` cures them independently. §4.2 at `k = 3, n = 48` takes 7,
  and `k = 2` converges on the same `n = 16` mesh where `k = 1` takes 42. The
  whole "stiff sources" narrative came from benchmarking at one under-resolved
  point.
* **What is really left is the coarse mesh, and §4.4.** Newton must start inside
  its basin and on a coarse mesh the Dirichlet datum is not — hence
  `Globalisation::PicardThenNewton`, which matters because an adaptive run has to
  solve *before* it has an estimator to refine with. §4.4 is a different animal:
  `∂F/∂ψ` spans `[−580, +566]` against a first Dirichlet eigenvalue of 22.3, so
  the linearised operator has swept past ~26 eigenvalues and the problem is
  multi-valued. Refinement cannot touch that.

Numbers in `CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles* and *Picard, then
Newton*.

*From `ROADMAP.md` *The nonlinear question is settled*.*

### Continuation in the driver

**Declined for the driver; a test-only investigation.** §4.4, the current hole,
is solvable: adaptive continuation in the added term's amplitude walks `c₃` from
0 to −18 at `k = 2, n = 32` in **9 solves with 2 retreats**, final solve 8
Newton iterations, no limit point on the branch. Ten uniform steps stall at
−10.8, so the step control is the whole of it.

**It ramps a constructor argument of the _test fixture_.** `Source` exposes
`f( R, z, ψ )` and `dFdPsi( R, z, ψ )` and nothing else, so a black-box `F` has
no amplitude to continue in and this route is unavailable to a driver on
principle. The expressible analogue `F_λ = λF` is a **different** homotopy —
degenerate at `λ = 0`, where this path starts from the converged pedestal — and
is untested.

So this stays a **test-only investigation**, and its value is that §4.4 is
known to have a solution MEQ can reach. Making it general needs a
self-parameterising `Source`: an interface change to argue on its own merits,
not to smuggle in.

*From `ROADMAP.md` item 7.*

### Methods that need Jᵀ: merit-gradient Armijo, dogleg, trust region, JFNK, full-space SQP

**Declined.** `mfem::DarcyNPCOperator::Jacobian::Mult` **aborts**
(`../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp:3756-3761`):
after `ComputeH()` the local arrays hold factored blocks, so the Jacobian can be
solved with and cannot be applied. There is no `Jᵀ` either.

That removes, structurally:

* a **merit-function gradient** `∇‖R‖² = 2JᵀR`, hence any Armijo test with a
  real directional derivative rather than the `(1 − 1e-4α)` surrogate that is
  there;
* a **dogleg** or any trust region needing a Cauchy point;
* **Newton–Krylov / JFNK**, which needs `J·v`;
* and anything that wants to *compare* `J` between iterates.

*From `BORDERED-GLOBALISATION-PLAN.md` §2.2, whose closing paragraph — what IS
available, and why the families there change the operator or the problem rather
than the step rule — stays in that file.*

**CEDRES++'s architecture is therefore not available to MEQ, and this is not a
preference.** Its SQP iteration (3.33) is a 3×3 block KKT system containing
`D_y Bᵀ(yⁱ)` explicitly in the first block row. MEQ cannot form that block, cannot
apply it, and cannot hand it to a direct solver. A full-space KKT solve is
**structurally unbuildable** on the NPC path.

*From `INVERSE-AND-RECONSTRUCTION-PLAN.md` §3.4, which verified the abort against
the installed header and keeps what IS available — `J⁻¹` applied to a right-hand
side and to each border column, blocked through `DarcyNPCSolver::ArrayMult`.*

### `Normalisation::Decoupled`: a Picard on the border round a Newton on the field

**Falsified.** `Normalisation::Decoupled` (`setNormalisationCoupling()`) zeroes the column
(`sourceColumn`, `GradShafranov.cpp:5934-5940`) and the row (`coupled` false,
`GradShafranov.cpp:6425`), leaving `step[0] = −constraint` — a Picard update on
`ψ_ax` wrapped round a full Newton on the field. **M-31 measures it**: coupled
reaches `4.447e-15` in 4 iterations, decoupled sits at `8.239e-02` after 15 and
does not converge.

That is the *opposite* split from the bordered Picard of
`BORDERED-GLOBALISATION-PLAN.md` §3: Newton on the field, Picard on the border.
It is worth stating because it is the obvious first idea and it is already dead,
and because it is the control that says the coupling is doing real work. **Do
not re-derive it.**

*From `BORDERED-GLOBALISATION-PLAN.md` §2.4.*

### KINSOL and Anderson acceleration on the bordered fixed point

**Declined.** M-36 and `CLAUDE_HDGGS.md` both say plain Picard needs damping and Anderson at
depth 1 does not. MEQ gets Anderson from `mfem::KINSolver( KIN_FP )`
(`solveByPicard()`, `GradShafranov.cpp:1309-1370`).

**KINSOL cannot be handed this fixed point.** `KINSolver::Mult` drives
`oper( x ) = 0` or `u = G( u )` over one vector; the bordered iteration's state
is `(x, p)` with `p`'s equations being the constraints, and MEQ's residual is
not zero (hence `ShiftedResidual`). Wiring it would re-create exactly the
"drives a residual of its own" problem the throw complains about.

So a bordered Picard is a loop inside `solveWithNormalisation()`, and Anderson —
if wanted — is a small dense least squares over the last `m` increments of
`(x, p)`. **That collides with `CLAUDE.md`'s standing preference for a
maintained library over a hand-rolled algorithm**, and the collision should be
stated rather than resolved by taste: KINSOL is the maintained implementation
and it structurally cannot take this problem; Eigen is already a sanctioned
dependency for dense least squares and would carry the QR. **The plan's recommendation
is to build the damped fixed point first with no acceleration (BG-3), measure
it, and only add Anderson if the iteration count is the thing that kills it** —
M-36 reads 122 to 290 Picard iterations on the unbordered stiff cases, and at
M-90's 1.1 s per non-linear step on DIII-D that is two to five minutes for a
reactive rung, which is affordable. **BG-3 is built that way and is falsified
on other grounds** — 0 of 6 close, [M-119](MEASUREMENTS.md#m-119) — so the
iteration count is not what kills it and Anderson has nothing to add.

*From `BORDERED-GLOBALISATION-PLAN.md` §3.6.*

### The bordered ladder's rejected-outright list

**Declined.** **Rejected outright, with reasons**: KINSOL on the bordered
residual (*KINSOL and Anderson acceleration on the bordered fixed point* — it
cannot take the border); a dogleg, trust region or JFNK (*Methods that need Jᵀ*
— the NPC Jacobian is solve-only); `Normalisation::Decoupled` as a rung
(*`Normalisation::Decoupled`* — measured failing, M-31); any predictive trigger
(*Predictive solver selection: the killed detectors* — two signals measured
useless, one of them anti-correlated, and a third that looks the same way).

*From `BORDERED-GLOBALISATION-PLAN.md` §7, item 5.*

## Linear solves, threading and the build

### UMFPACK's threaded BLAS as the `MKL_NUM_THREADS` culprit

**Falsified.** The `MKL_NUM_THREADS` collapse is MFEM's element-local dense LU in `ComputeH()` ([M-59](MEASUREMENTS.md#m-59), `CLAUDE_HDGGS.md`, *`MKL_NUM_THREADS`, measured*), not UMFPACK.

**THE CULPRIT IS NOT UMFPACK.** The plausible story — UMFPACK calling threaded BLAS on the small dense
frontal matrices of a sparse LU, thousands of times, each call costing a fork
and a barrier — fits every whole-test fact and is **wrong**. Separated out,
`UMFPackSolver::SetOperator` degrades by about **40%** across the whole thread
range, never more.

**Keep the wrong answer in view, because of HOW it was wrong.** It fitted every
fact then available — the collapse tracked exactly the tests that reach a direct
solver — and **no whole-test timing could have refuted it**, because in those
tests UMFPACK and `ComputeH()` always appear together. Only separating the
columns did. Same lesson as `CLAUDE_HDGGS.md`, *The Solov'ev coefficients were wrong
twice*: checking a story against
the evidence that suggested it proves nothing.

*From `CLAUDE.md`, “Traps”.*

### MKL_THREADING_LAYER=GNU

**Declined.** Keeping it set as a free guard was rejected.

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

*From `CLAUDE_HDGGS.md`, “PARDISO and the MKL link line: what is still true”.*

### Symmetric solvers, flux-block Cholesky, AMG, and re-timing cuDSS

**Declined.** Items 2 to 5 of the ranking in `CLAUDE_HDGGS.md`, *What to do, in order of value*; item 1, preconditioning the Newton path's GMRES fallback, stays there. The ranking is the one the ASYMMETRY finding gives: ranking Cholesky highly is what measurements taken only on the fitted path suggest, and it is wrong wherever `Γ` is curved.

2. **Symmetric methods — fitted path only.** CG or MINRES on `−A`, and a CHOLMOD
   Cholesky, are all correct on a fitted mesh and all invalid on the extension
   path. CHOLMOD is already linked (`-lcholmod` via SuiteSparse) but MFEM wraps
   only UMFPACK and KLU. Worth having only if fitted-domain runs become a
   workload in their own right; the driver's target is curved boundaries, where
   none of it applies.
3. **Cholesky on the per-cell flux block — do not.** `(R q, v)` is SPD, so this
   looks like a free 2×. It is wrong on the extension path for the reason
   `CLAUDE_HDGGS.md`, *The trace matrix is symmetric on the fitted path and is
   not on the extension path*, gives: `AssembleFluxMassBdrFaces()` puts a term with relative
   asymmetry 1.0, sixteen times larger than the mass term, into the very block
   `LU_A` factors. A correct version would have to branch on whether an
   extension is configured, for at best a small gain on 9×9 to 30×30 blocks.
   Not worth a conditional in that inner loop. **`LUFactors` is the right
   choice and should stay.**
4. **Do not reach for AMG.** At 2D serial sizes direct wins: 50k trace dofs
   factorise in 0.4 s. AMG earns its place in 3D or in parallel, and serial MFEM
   has none anyway — `HypreBoomerAMG` needs MPI.
5. **Do not re-time cuDSS here.** Correctness is established, the timings are
   irreproducible on this machine by a factor of thirty, and the reason is the
   machine. It needs a datacentre part, not another afternoon.

*From `CLAUDE_HDGGS.md`, “What to do, in order of value”.*

### An iterative trace solve: GMRES, and Ginkgo

**Closed.**

**AND `TraceSolver = "gmres"` IS A FOURTH CHOICE, EXPOSED TO BE MEASURED AND
MEASURED AS UNVIABLE** → **[M-179](MEASUREMENTS.md#m-179)**. GMRES with a
Gauss-Seidel smoother needs **247 Krylov iterations** on the easiest case in the
tree and **4,179** on `machine-f-diiid`, where the run does not finish in 1,800 s
against PARDISO's 3.0 s. It also gives up the two things the direct solvers are
chosen for — one factorisation amortised over the bordered step's fourteen
columns (M-98's 1.34×) and the symbolic reuse kept across Newton steps. **So an
iterative trace solve is closed rather than pending**, and Ginkgo — which is
otherwise well fitted, being the one MFEM solver package that is not MPI-gated
and that takes an `mfem::Device` directly — would have to close that gap with a
preconditioner rather than with hardware. **And the driver silently ran UMFPack
for the first measurement**, because `apps/meq.cpp:1634` initialises
`TS traceSolver = TS::UMFPack` before a `switch` with no `default`; see `CLAUDE_DRIVER.md`,
*Configuration traps*.

*From `CLAUDE.md`, “Status: the solver works, and MEQ is a program”.*

### MFEM as a submodule

**Declined.** toml11 is a pinned submodule; MFEM is not.

**MFEM stays out of tree**, hand-built, found through `MFEM_DIR`. Its history is
enormous, so a submodule would make every clone of MEQ expensive; it needs an
out-of-tree configure-and-build of its own anyway; and `../mfem-hdg-dev` is a
working tree with its own active development on `gf-hdg-subdomains-dev`, which a
submodule pin would fight rather than help.

*From `CLAUDE.md`, “Why toml11 is a submodule and MFEM is not”.*

### PARDISO's old failure: the sequential timings and three dead hypotheses

**Withdrawn, and falsified.** Do not quote the PARDISO timings `TODO` carried
before 2026-08-30. They were taken when a stray `libmkl_sequential` made MKL
resolve sequential for everything, so they were sequential-against-sequential
and were never labelled as such. The live numbers are in `CLAUDE_HDGGS.md`,
*Threading, measured*.

**Three hypotheses were tested against the old failure and are all dead**, so
do not re-test them if a new one appears: a METIS symbol collision with MFEM's
`-lmetis`; an interaction with MFEM's link line; and something particular to
MEQ's matrix — a plain 5-point Laplacian failed identically.

*From `TODO` *PARDISO: settle reproducibility*.*

### Threading the output writers (E) and the Newton loop's remainder (F)

**Declined, and the decline is a decision rather than a deferral.** Re-checked
against [M-138](MEASUREMENTS.md#m-138) and unchanged. If either grows, that is a
new plan and not a reopening of `THREADING-PLAN.md`.

**E. The output writers, 3.8% — do not do.** `output` less `postProcess` is
0.269 s. [M-92](MEASUREMENTS.md#m-92) splits it on this same case: the `.vtu`
0.130 s, the `.mesh` and three `.gf` 0.053 s, the NetCDF write 0.016 s, and the
whole `GridSampler` path — constructor and three sampling passes — **0.010 s**.
That is serialization and file I/O, not compute.

The one loop in it with real arithmetic is the sampling pass, and M-92 observes
that grouping the located nodes by element makes it a `forall` with no
reduction and no contention, and records that it is **not done** because the
whole grid path is 27 ms. That judgement stands: threading 10 ms buys 0.14% of
a run.

Writing the three files concurrently would take 0.20 s to 0.13 s, about 1.0%.
MFEM's stream I/O is not documented reentrant and the three writers read the
same `GridFunction`s; 1% is not worth finding out. **Do not do.**

The part of this phase that *is* compute is `postProcess()` —
`DarcyForm::Reconstruct()`, four integrators re-assembled at the enriched order
per element. It is MFEM's and out of scope for MEQ's threading, and **it is the
one number in this item the budget moves**: 13.0% in the budget the plan was
written against, **18.5% at 1.00 cores** in [M-138](MEASUREMENTS.md#m-138).
Upstream has threaded it, worth `1.177×` on the whole run
→ [M-172](MEASUREMENTS.md#m-172). The writers around it are unchanged —
`output` less `postProcess` is 0.198 s, **5.8%**, and still serialization rather
than arithmetic — so **the decline stands.**

**F. `other (remainder)`, 2.7% — not addressable by a loop MEQ can write.** This
is MFEM `Vector` arithmetic and allocator churn inside the Newton loop:
`unknown = savedState`, the `Add()` chain in each line-search trial,
`Norml2()` in `augmentedNorm`, and the temporaries underneath them.

* **The `n`-length vectors are already hoisted.** `residual`, `column`, `y`, `z`,
  `scratch`, `columnL`, `zL` and `currentRow` are declared once per solve,
  outside the Newton loop. There is nothing to hoist.
* **`Vector::Add` and `Norml2` are MFEM's**, dispatched through `mfem::forall`.
  Threading them on the host means configuring an `mfem::Device`, which is a
  process-wide decision with its own hazards ([M-79](MEASUREMENTS.md#m-79)), and
  `Norml2` would then reassociate.
* **The measured MEQ-side lever is the build, not a loop.**
  [M-77](MEASUREMENTS.md#m-77) measures the CUDA-enabled install at **1.07× the
  wall and +55.4 M allocations** against an otherwise identical
  `../mfem/install-nocuda`, with byte-identical output, and attributes 24% of the
  allocation traffic to `mfem::forall`'s host-lambda wrapper. **7% of the whole
  run for one CMake variable is larger than every item `THREADING-PLAN.md` built
  except A**, and it costs nothing to take for a production CPU build.

**Do not thread this leg. Take the no-CUDA install instead, and say in the
release notes that a CPU production build is a separate install.** Re-checked
against [M-138](MEASUREMENTS.md#m-138): the leg is 0.116 s, **3.4%**, at 1.06
cores, and M-77's 7% is still one CMake variable.

**Nor the two Γ sweeps in `border assembly`**, `exteriorTransmissionRows()` and
`exteriorConductorMoments()`, which run three times each over 82 boundary faces
and are already below the noise.

**When E or F would be worth revisiting, and it is a new plan when it
happens.** A reader meeting a bigger number later needs to tell whether it means
anything, and `THREADING-PLAN.md` produced the counter-example itself:
`re-assembly` went from 4.2% to 5.6% of the run **without the leg moving at
all**, purely because the run got shorter around it.

**A SHARE IS NOT A TRIGGER. SECONDS ARE.** Every share in a leg budget has the
wall clock as its denominator, so threading anything at all inflates every share
that did not move. Re-deriving M-138 is what makes one number comparable with
another; comparing a share across budgets is not.

* **E — the output writers.** The leg is 0.198 s and **0.010 s of it is
  compute**; the rest is serialization and file I/O. So the trigger is not the
  leg, it is the `GridSampler` path inside it, and the two things that would
  move it are a much larger `[output] GridNR`/`GridNZ` — it is `O(nodes)` — or a
  case whose mesh makes the location work dominate. **Revisit at 0.1 s of
  sampling**, which is ten times what it costs now, and revisit the
  concurrent-write half only if somebody establishes MFEM's stream I/O is
  reentrant, which is the part that is not about size at all.
* **F — `other (remainder)`.** MEQ cannot thread it without configuring a
  process-wide `mfem::Device`, so its trigger is not a size at all: **it becomes
  reachable for free the day MEQ configures a Device on the host path for some
  other reason**, because `Vector::Add` and `Norml2` already dispatch through
  `mfem::forall`. Until then the lever is M-77's build variable and not a loop.
  If a Device does arrive, `Norml2` would reassociate and the augmented norm
  feeds the Armijo test, so `THREADING-PLAN.md`'s cross-cutting exactness rule
  applies to it.

**And either way it is a new plan.** `THREADING-PLAN.md` is closed, its budget
is superseded, and an item reopened against a stale denominator is how a closed
plan quietly becomes a wrong one.

*From `THREADING-PLAN.md` §E, §F, *When E or F would be worth revisiting*, and
*What to do first* (*Not to be done*).*

## Device and single precision

### cuDSS determinism: SetDeterministic and per-column Mult

**Declined.** cuDSS 0.8.0 leaves `CUDSS_CONFIG_DETERMINISTIC_MODE` at **0**: same matrix, same
right-hand sides, fresh handle, five repeats, **4 of 4 differ**, at 6.7e-16 /
1.8e-15 / 3.6e-15 for `K = 1 / 4 / 16`.

**`SetDeterministic( true )` is not the
answer and must not be reached for.** It is supported for **one right-hand side
only**; the solve phase refuses 4 and 16 columns with
`CUDSS_STATUS_NOT_SUPPORTED`, and MEQ's bordered step is `N + 4` columns through
`ArrayMult` — *the one configuration MEQ would want determinism in is the one
where it cannot have it*. With `CUDSS_MTYPE_SYMMETRIC` the switch is
additionally **silently wrong**, returning an exactly zero solution that every
phase reports as success. Buying determinism back with column-by-column `Mult()`
gives up M-98's 1.34× and validates a configuration that is not the one shipped.

*From `CLAUDE.md`, “Status: the solver works, and MEQ is a program”.*

What this copy adds: 0.8.0 is the newest cuDSS NVIDIA ships, so
this is the library as it stands rather than something about to be fixed. With
`CUDSS_MTYPE_SYMMETRIC` the switch is silently wrong, which is why upstream's
setter refuses that type by name. **MEQ is clear of both hazards and it is worth
saying why rather than assuming**: it builds `NONSYMMETRIC` / `FULL`, and it
calls `SetMatrixSymType()` before `SetMatrixViewType()` at both of its two call
sites — `GradShafranov.cpp:451` and `TraceSolverScaling.cpp:470` — which is the
order upstream's own order-dependence defect requires and which a caller had no
way to know was load bearing.

*From `DEVICE-PLAN.md` §3.5 (*cuDSS, and why its own determinism cannot be
bought back*).*

### A per-pair rate criterion, or a wider `rateSlack`, in single precision

**Falsified.** The per-pair criteria are falsified and widening `rateSlack` is declined; the skip is by degree.

**`k+1` ABOVE `k = 1` IS A DOUBLE-BUILD ACCEPTANCE AND A SINGLE BUILD CANNOT
CARRY IT.** Against the double build on the same tree, `Solov'ev NSTX` in the
flux is 1.09× contaminated at `k = 1` on the finest mesh, **203×** at `k = 2`,
and **1.5× on the COARSEST mesh** at `k = 3` — so coarsening leaves the
asymptotic regime and refining goes deeper into round-off, and there is no
window. **Three per-pair criteria were tried and each provably fails**: no
constant separates a 9%-contaminated pair reporting a good 1.873 from a
36%-contaminated one; *"the error must at least halve"* is necessary but lets
through a pair that halves and still reads **3.03** against 3.85. The tail
**bends before it flattens**. So the skip is by DEGREE, loudly, and the
**absolute ceiling is still asserted** — widening `rateSlack` until the tables
passed would weaken every mesh in every build to accommodate meshes that are
not measuring anything, which is what `CLAUDE.md`, *Testing stance*, forbids.

*From `CLAUDE.md`, “Which integrators need device offload, and the batched routes”.*

### A GPU performance verdict from this machine

**Deferred to other hardware.** The partial device path is built and it is
green — `DeviceParity` passes, and with the batched keys held equal
`--device cuda` reproduces `--device cpu` on every printed line but the final
residual, at 5.3e-17. What this machine cannot supply is the *verdict*: consumer
FP64 runs at 1/32 to 1/64 of FP32 where a datacentre part runs it at about 1/2,
so a local timing can invert the production conclusion. Take the numbers
elsewhere before concluding anything about whether the device is worth it.
**[M-223](MEASUREMENTS.md#m-223) is where they stand here** — the solve leg at
parity with the host-threaded default, and the whole 1.07–1.13× wall loss in the
**output** stage, which is `Reconstruct()` losing its host threading under a
Device rather than anything in the solve.

**Timing anything on this machine and calling it throughput** does not change
under residency either. [M-77](MEASUREMENTS.md#m-77) measures the CUDA build
costing 7% of the wall on the CPU path for byte-identical output. M-223 measures
the partial device path's SOLVE at parity with the host-threaded default, and
the residual 1.07–1.13× on the wall as one host leg. Neither is a measurement of
what a datacentre part would do, and `DEVICE-PLAN.md` §7 is the only thing that
makes a local device timing worth taking at all.

*From `ROADMAP.md` *Deliberately not yet*, and `DEVICE-PLAN.md` §3.5 (*What does
NOT change*); two copies, merged.*

### `TraceSolver = "cudss"` from a configuration file alone

**Declined.** The gate is the device, not the key: a run that configured none is
refused and `--device cuda` honours it. **This waits on nothing of upstream's** —
it rests on the trade, which is unchanged. M-126 puts the trace solve at 20.1%
of a threaded step, so configuring a device to reach cuDSS and nothing else buys
a share of a fifth while paying a bus crossing on the other four fifths; and
without a device cuDSS does not fall back but dies inside CUDA naming
`cudaMemcpy` and not the key. Under `--device cuda` cuDSS is the default trace
solver, `DEVICE-PLAN.md` §3.3 item 6, and the digit-for-digit instrument
survives it: cuDSS 0.8.0 is nondeterministic run to run, below the printed
precision, and over eight repeats every field `DeviceParity` compares is
identical → [M-229](MEASUREMENTS.md#m-229).

*From `ROADMAP.md` *Deliberately not yet*. The bullet's closing claim — that
cuDSS costs the instrument and is ordered last for that reason — is corrected
here against M-229, which took item 6 and measured the instrument surviving.*

### Partial assembly

**Declined: a discretisation decision, not a flag.** Not implemented for the HDG
integrators at all, and its case is sum factorisation, which wants
tensor-product elements where MEQ is triangles. **Not the blocker it reads
as**: the device route is MFEM's batched element and face kernels, reached
through `[solver] AssemblyMode`, `LocalFactorMode` and `TraceAssemblyMode`, and
MEQ is on it. Residency does not need it — batched dense kernels are the route
— so it stays out of `DEVICE-PLAN.md` entirely rather than moving up it
(`DEVICE-PLAN.md` §2.4).

*From `ROADMAP.md` *Deliberately not yet* and `DEVICE-PLAN.md` §3.5 (*What does
NOT change*); two copies, merged.*

### `AssemblyMode::Batched` as a host default

**Declined.** [M-99](MEASUREMENTS.md#m-99)'s 15–27% loss stands as a host
measurement and MEQ's host default stays `Serial`; residency does not change
that. `DEVICE-PLAN.md` §3.2 says why that is not a verdict on the device route.

*From `DEVICE-PLAN.md` §3.5 (*What does NOT change*).*

### A fully device-resident solve as a commitment

**Deferred.** `DEVICE-PLAN.md` is the audit and it is explicitly not a
commitment. Where it stands, **asked of the library with `--profile` rather than
inferred** ([M-223](MEASUREMENTS.md#m-223)): MEQ's own device integrand, the
batched nonlinear residual and gradient, the batched local solve, the whole NPC
traversal, the condensation cache and the batched trace assembly are all on
under `--device cuda`. Face assembly is serial **deliberately** — its kernels
accumulate `D` through `AtomicAdd`, which is what would destroy the
bit-reproducibility the device programme protects
([M-215](MEASUREMENTS.md#m-215)). The bordered columns stay on the device and
the dense fill contracts them in one launch ([M-238](MEASUREMENTS.md#m-238)), so
what is left is the output stage. Three MEQ-specific facts also keep upstream's
**1.33× at 4096 elements** from being MEQ's crossover: their fixture is quads and
MEQ is triangles everywhere, five of MEQ's eight shipped meshes are *below* 4096
elements, and MEQ's Newton is bordered with the border being host arithmetic
that would not go to a device.

*From `ROADMAP.md` *Deliberately not yet*. Its "the bordered row dots stay on the
host … what is left is the output stage and the border" is corrected here
against M-238.*

### Host threading as a device-plan item

**Closed.** Host threading of the NPC traversal
([M-135](MEASUREMENTS.md#m-135), 1.25×) and of `Reconstruct()`
([M-172](MEASUREMENTS.md#m-172), 1.177×) are both landed and both are **host**
parallelism that a resident run replaces rather than uses. Nobody should size
the device work against them.

*From `DEVICE-PLAN.md` §3.3 (*What is no longer on this list*).*

### Digit-for-digit agreement as the device instrument

**Falsified.** `DEVICE-PLAN.md` §3.5 spent its length protecting *"`--device
cuda` reproduces `--device cpu` to every printed digit"*, and the record does
not support the weight put on it. Every device-path defect either project has
found was enormous:

| | the disagreement |
|---|---|
| [M-79](MEASUREMENTS.md#m-79) | the potential identically **zero**, 0/0 Newton iterations |
| [M-130](MEASUREMENTS.md#m-130) | *"the bordered Jacobian is singular"* — a refusal, not a number |
| [M-183](MEASUREMENTS.md#m-183) | the exterior coupling dropped; `psi_ax` **10%** wrong |
| [M-205](MEASUREMENTS.md#m-205), [M-206](MEASUREMENTS.md#m-206) | **a factor of two** |
| [M-189](MEASUREMENTS.md#m-189) | **1.27e-06** against a 7.8e-16 history — the tightest ever recorded |

**A 1e-10 tolerance would have caught all of them, and the last twelve digits
have never once been what fired.** Worse for the claim: the one defect
digit-for-digit agreement demonstrably **could not** catch is
[M-187](MEASUREMENTS.md#m-187), where every printed number was identical and
the work was 4× — found by an iteration count, and root-caused by a probe on
the Newton increment ([M-207](MEASUREMENTS.md#m-207)).

**So what finds device defects is a reference arm you trust, and the tightness
of the comparison has been nearly irrelevant.** The plan conflated two things,
and only the first is load bearing for debugging; what tightness genuinely buys
— interpretable diffs, a tolerance nobody can renegotiate, re-derivable
published numbers — stays in `DEVICE-PLAN.md` §3.5.

*From `DEVICE-PLAN.md` §3.5 (*The instrument this section was ordered around is
weaker than it claimed*).*

### Items 6 and 8 strictly last, as uniquely destructive

**Withdrawn.** cuDSS as the device default (item 6) and single precision (item
8) were ordered strictly last on the premise that they alone take the
digit-for-digit check away; that premise is gone, because the check was never
the thing doing the work. **The constraint that replaces it is one sentence: a
reproducible arm must survive every step, and a step that would leave MEQ with
none of them does not go in.** cuDSS satisfies that — it is a *trace solver*
choice, so a PARDISO arm remains available on the same build — and so does
single precision, whose arm is the double install. What both owe is §7.2's
harder problem: an acceptance stated as a tolerance rather than as agreement,
and a tolerance loose enough to absorb a 3.6e-15 solve that differs every run
underneath a bordered Newton that differences its constraints through it. Both
items are taken: under cuDSS every field `DeviceParity` compares is identical
over eight runs ([M-229](MEASUREMENTS.md#m-229)), and the single build is at
parity with the double ([M-222](MEASUREMENTS.md#m-222)).

*From `DEVICE-PLAN.md` §3.5 (*What this changes about the ordering*). Its "the
default device configuration already has [taken the check away]" is dropped:
`AssemblyMode` is `Serial` under a Device, which restores bit-reproducibility
([M-215](MEASUREMENTS.md#m-215)).*

### A missing device route as the gate

**Falsified.** `AssemblyMode = "batched"` is the device route, and MEQ exposes it
as a TOML key. `DarcyForm::Assemble()` calls
`AssembleFluxMassMatricesBatched`, `AssembleDivMatricesBatched`,
`AssemblePotMassMatricesBatched` and `AssembleFluxMassBdrMatricesBatched`
**in front of** each of the three element loops; each returns false having done
nothing unless `AssemblyMode::Batched` is asked for, and only then does the
caller fall back to `BilinearForm::ComputeElementMatrix()`. So the sentence
`DEVICE-PLAN.md` quoted as the blocker — `darcyform.hpp`'s note that hybridized
assembly goes through `ComputeElementMatrix()` — describes **what happens when
the key is off**, which is a default and not a wall.

**The admission is narrow and it is exactly MEQ's.** `HDGElementMassCanBatch()`
takes `MassIntegrator` and `VectorMassIntegrator`; `HDGElementDivCanBatch()`
takes `VectorDivergenceIntegrator`. `buildForms()` installs
`VectorMassIntegrator( radius )` and `VectorDivergenceIntegrator`. Nothing had
to be written for MEQ to fit; it already did.

**Run on the production case rather than argued**
→ **[M-185](MEASUREMENTS.md#m-185)**. With all three batched keys set,
`machine-f-diiid` prints `face assembly yes, local factorisation yes, local
solve yes, trace assembly yes` and reaches `psi_ax = 3.759851e-01`, the host
default's answer to every digit. **Seven of the ten batched routes are taken and
the one refusal left is `CanBatchLinearResidual()`** — `meq::SourceIntegrator`
being genuinely nonlinear, which is Request 3 and not a route.

**How the plan got it wrong is the transferable part, and it is the plan's own
rule broken.** `DEVICE-PLAN.md` §6 says every claim about somebody else's code
is taken by extracting whole functions, *"never by reading a line range —
`CLAUDE.md` records a false claim produced exactly that way"*. This one was a
**doxygen comment on a member declaration**, read without reading the function
it describes or the call site in front of it. The comment was true; it was an
answer to a different question. **A comment is documentation of an intent, and
the call site is the fact** — and `grep`ping for the batched entry points from
`DarcyForm::Assemble()` outward is the ten seconds that separates them.

**And one thing that looked like a disagreement with upstream was a
configuration.** Upstream's predicate table reports `CanBatchLocalSolve` as
**YES**; MEQ's own diagnostic printed **NO** on the same problem. The predicate
opens `lfac_mode == LocalFactorMode::Batched`, upstream set all three keys and
MEQ set one, and both readings are correct. **Before reporting that a library
disagrees with its own table, check which keys the table was taken with.**

*From `DEVICE-PLAN.md` §3.6.*

### Device findings MEQ does not file upstream

**Declined.**

* **Anything measured against `gf-hdg-linearise-first`'s unfinished offload
  work.** The standing rule is not to report findings against work that has
  not landed.
* **`DEVICE-PLAN.md` §1.2's three candidate sites**, which are not upstream's at
  all. They are MEQ's own aliases, and **no unit test in MFEM's tree can see a
  caller's long-lived aliases onto the solution blocks** — upstream says so in
  its own words about [M-79](MEASUREMENTS.md#m-79).

*From `DEVICE-PLAN.md` §5. "Six of that plan's nine items landed only recently"
is dropped as dated.*

## Free boundary: borders, limiter, X-point and the plasma support

### Confinement, clamping and a prescribed current: the M-134 two-by-two

**Falsified.** Candidate repairs for `theTwoBorderSolveReportsATrueMagneticAxis`'s fixture, each run on its own and then as the full 2 x 2; the repair that worked is coils and a prescribed current with the confinement (`CLAUDE_TRAPS.md`, *Reading a red test*).

**AND "THE FIXTURE CANNOT BE REPAIRED IN PLACE" IS THE WRONG CONCLUSION TO DRAW
FROM IT BY ONE INGREDIENT.** `ConfineToPlasma` alone fails 4/4 and clamped
profiles alone fail 4/4; both measurements stand, and neither is a measurement
about the pair — the same run records "clamped WITH a prescribed current
converges in 22 steps". **A measurement that rules two things out is not a
measurement that rules everything out.**

**THE FOURTH CELL IS NOW RUN AND IT SAYS THE PAIR IS NOT ENOUGH EITHER** →
**[M-134](MEASUREMENTS.md#m-134)**. Confinement **with** a prescribed current,
swept over 250x in the target, reaches **2 of 4 at best** — 0 of 4 below
`μ₀I_p = 0.20`, 2 of 4 at 0.50 and again at 2.00, and back to 0 of 4 at 5.00, so
the maximum is interior and the sweep is not truncated. The clamped arm is cured
outright by that same lever. So a prescribed current moves confinement off zero
and does not repair it, and `theTwoBordersConvergeTogether` keeps
`FREE-BOUNDARY-PLAN.md` §11.2's retraction of its `ψ_ax` and `ψ_bnd` as physics.
**The lesson survives its own example**: it took the full 2 x 2 to say this, and
each of the four one-key runs on its own supported a different wrong reading.

*From `CLAUDE.md`, “Commands”.*

**Declined.** §7.12b's two-border fixture: four limiter radii, `ψ_ax` in every
row, and the column is an artefact of the `R = 0` layer. The convergence claims
in that section stand — the residuals, the iteration counts and `ψ_ax`'s own
constraint at 1e-17 are statements about the solve closing, and it closes. **The
`ψ_ax` and `ψ_bnd` values are not physics.**

**RETRACT, AND IT CANNOT BE RE-PUBLISHED.** §11.1 supplied the located-axis
numbers — 6.513594e-02, 6.420050e-02, 4.420624e-02 and 6.111920e-02 at the four
radii — but **replacing the column with those would be wrong too**: the O-point
flux is the axis flux of *the field that was computed*, and §11.3 shows that
field carries a `1/R` pole on the axis, so it is not the equilibrium `[source]`
describes either.

**AND THE FIXTURE CANNOT BE REPAIRED IN PLACE — WHICH IS A FIXTURE DEFECT AND NOT
A CAPABILITY GAP, AND THE DISTINCTION MATTERS.** §7.16's machine case has every
ingredient this one has — a domain reaching `R = 0`, a limiter, an exterior
coupling and `ConfineToPlasma` — and converges in 11 Newton steps. What it has
that this does not is **coils and a prescribed current**. Made physical while
keeping neither, this fixture becomes the amplitude-fixed moving-support problem
§7.14 records as a non-linear **eigenvalue** problem, which is ill posed rather
than merely hard — and measured, **clamped profiles WITH a prescribed current
converge in 22 steps** where clamping alone fails 4/4. The cure is §7.14's own.
Both physically-correct repairs were tried at all four limiter radii, 300 Newton
iterations, everything else held: `ConfineToPlasma` — **4 of 4 do not
converge** — and profiles clamped to vanish at `Ψ ≤ 0`, the same statement made
in the profile rather than by a support test — **4 of 4 do not converge**. So
`theTwoBordersConvergeTogether`'s result *depends on the unphysical
extrapolation*: with a source carrying no current in the vacuum, that
configuration does not close at all. FB-T is the fourth cell, confinement WITH a
prescribed current, and it reaches 2 of 4 at best →
[M-134](MEASUREMENTS.md#m-134).

**So the `ψ_ax` and `ψ_bnd` columns of §7.12b are retracted as physics and kept
as a record of the solve closing.** The convergence claim stands and is what that
section is about.

*From `FREE-BOUNDARY-PLAN.md` §11.2; section numbers are that file's. §11.2
keeps its heading and the verdict.*

### Conditioning, units and combinatorics as the cause of the bordered stall

**Falsified.** The stall that `setNormalisation( s )` — one argument where two
were meant, zeroing `ψ_bnd` inside the Jacobian window — produced was first read
three other ways. `FREE-BOUNDARY-PLAN.md` §7.16 named unscaled border rows in SI
as the prime suspect. It is wrong, and so are the two obvious alternatives:

* **Not conditioning.** The dense border matrix has `cond( M ) ≈ 1.0e3` — three
  digits of sixteen — and the corner solves its own system to **8.4e-16**.
  Equilibration is worth 36×, not eleven orders: the SI dynamic range never
  reaches `M`, because every border row is already a ratio.
* **Not units.** Non-dimensionalising the entire problem — `μ₀ = 1`, `ψ` scaled
  to `O(1)`, currents and amplitudes rescaled — reproduces the SI run **to ten
  digits**, floor and all. It neither fixes nor perturbs it.
* **Not combinatorial.** In the stalled phase the `ψ_ax` argmax is **frozen** for
  all 21 stalled steps, the support creeps one way only (39 points in, 0 out),
  and the line search returns the identical verdict every step. The control is
  sharper still: the configuration that **converges** is the combinatorially
  noisier one — the argmax jumping three times, hundreds of support points
  flipping **both ways**.

The cause was the one-argument call, measured against the Newton direction's own
linearised system → [M-03](MEASUREMENTS.md#m-03).

*From `FREE-BOUNDARY-PLAN.md` §7.17 and `CLAUDE_FB.md`, *One argument where two
were meant*; the two copies are merged.*

### Starting from vacuum, and amplitude continuation

**Declined.** Asked after `FREE-BOUNDARY-PLAN.md` §7.14's coil sweep failed, and
worth answering because the obvious homotopy is the wrong one.

**A VACUUM START DEGENERATES EVERY BORDER UNKNOWN AT ONCE, WHICH IS WHY IT IS
NOT THE HOMOTOPY TO PICK.** At `F_plasma = 0` there is no plasma, and therefore:

* `ψ_ax` — the largest nodal value — is a maximum of the **coil** field, not a
  magnetic axis, and on a vertical-field set it sits on the boundary;
* `ψ_bnd` — the flux at the limiter — is a number, but not the edge of anything;
* the support `{Ψ > 0}` is empty or arbitrary, and `span = ψ_ax − ψ_bnd` has no
  reason to be non-zero.

**The first two are CONFIRMED** (§7.18): the only critical points are the two
coil centres, minima at `( 1.1966, −0.6965 )` and `( 1.1912, +0.6954 )`, and the
grid maximum of `ψ_h` is `−6.06e-06`. **What is FALSE is that
`setNormalisation()` refuses it**: measured, the vacuum solve **converges in 2
Newton steps to 7.03e-17** with span `8.19e-03 − ( −3.01e-02 ) = 3.83e-02`, and
reports three unknowns that mean nothing. A solve that succeeds and describes
nothing is worse than one that refuses.

So the three unknowns the continuation exists to carry are exactly the three that
stop meaning anything at `λ = 0`. This is the same objection recorded against
`F_λ = λF` for the current hole — *"at `λ = 0` it degenerates to the harmonic
problem"* — and it is sharper here, because there the degeneracy cost a starting
point and here it costs the border its unknowns.

**AND CONTINUATION IN THE RHS HAS THE SAME DEFECT.** Ramping the profile
amplitude is the vacuum start in slow motion: the plasma vanishes at one end.
Ramping **`I_p`** does not — the plasma is present at every step, all three
unknowns keep their meaning, and the parameter is **prescribed input** rather
than something recovered from a black-box `F`, so it carries none of the
objection recorded against continuation in the driver; `setPlasmaCurrent()`
makes it a one-line loop. **But the `I_p` ramp is NOT BORNE OUT** (§7.18): at
coil `−0.05`, limiter 0.90 — the one place a cold solve fails — a ladder
converges at `μ₀I_p` = 0.01 … 0.04 and **fails at 0.05** of the 0.35 wanted, on
a coarse and a refined ladder alike, with every converged rung degenerate
(`ψ < 0` across the midplane, scale ≤ 5.7e-04).

**THE REAL LESSON IS THAT §7.14's FAILURE WAS NOT A CONTINUATION PROBLEM.**
`../freegs4e` never meets it, and checking why is decisive: its `constrain`
object is a **control system that solves for the coil currents inside every
Picard iteration** — `control.py` assembles a least-squares system over the coil
currents against shape constraints (`Br = 0` and `Bz = 0` at a prescribed
X-point, isoflux pairs) and applies it each step. **The coils adapt to the
plasma.** MEQ's §7.14 attempt did the opposite: it *prescribed* currents chosen
by eye and asked a cold Newton to find a core plasma consistent with them.

So there are two honest routes, and neither is a vacuum start:

1. **Choose currents that are consistent with the equilibrium wanted**, rather
   than continuing toward them. The vertical field a given `I_p` needs is a
   closed form — Shafranov's
   `B_v = μ₀I_p/(4πR)·[ ln(8R/a) + β_p + l_i/2 − 3/2 ]` — so the conductors and
   the initial guess can be made to agree **before** the solve rather than
   discovered to disagree during it. This is not continuation at all; it is
   picking inputs that describe one machine. §7.18 found that the cold failures
   this rested on do not occur once §7.17's defect is repaired, and §7.16
   supplies much better evidence for the same conclusion: with a genuine
   equilibrium's own coils and profiles the amplitude-fixed problem closes in
   five Newton steps.
2. **Solve for the currents**, as freegs4e does. That is the inverse solve, and
   it is built: IV-0 to IV-6, with coil currents as border unknowns of the same
   Newton ([M-250](MEASUREMENTS.md#m-250) to [M-253](MEASUREMENTS.md#m-253);
   `INVERSE-AND-RECONSTRUCTION-PLAN.md`).

**THE DIAGNOSTIC THAT TELLS WHICH BRANCH AN ITERATE IS ON** was proposed
alongside: on the annular branch `ψ_ax`'s argmax sits near `Γ`, on a core it is
interior — one integer per iteration, distinguishing "not converging" from
"converging to the wrong equilibrium", which a residual cannot. §7.18 found it
**CONFIRMED for the annulus and MISLEADING otherwise**: where the solve does not
find the wall-hugging O-point, `ψ_ax` is attained at a single node in the corner
element where `Γ` meets the axis. Read it against the field's own maximum or it
reports the corner element rather than the branch.

*From `FREE-BOUNDARY-PLAN.md` §7.15 and `CLAUDE_FB.md`, *The plasma current is a
border unknown*; section numbers are `FREE-BOUNDARY-PLAN.md`'s. The stale items
are corrected from §7.18's re-measurement: the vacuum start is not refused, the
`I_p` ramp is not borne out, and the inverse problem is no longer out of scope.*

### PicardSweeps and BorderedPicardThenNewton

**Falsified.** **`PicardSweeps` IS THE CHEAPEST UNTRIED THING AND IT IS A
TRAP.** It had never been set by anything in the tree. It fixes none of the six
cold `exit 2` cases and turns **three of them from `exit 2` into `exit 0` with a
wrong answer and no warning** — C shaped out by a factor of 3.9, both E rows by
23%. A failure is actionable; a green run reporting a different machine is not.
→ [M-116](MEASUREMENTS.md#m-116).

*From `CLAUDE_FB.md`, *The cold failures are three different things*.*

### Globalisation, and a transferable merit weight, as the repair for the cold failures

**Falsified.** **NOT ONE OF THE SIX COLD `exit 2` CASES HAS A NON-FINITE
DIRECTION**, which is the condition `BORDERED-GLOBALISATION-PLAN.md` §0.1 names
as the one that would make the whole globalisation campaign moot. They are three
distinct failures: two machines lose a SINGLE step of twenty-seven; one is
chronic, exhausting the halving ladder on 55 of 288; and three are not
step-length failures at all but **X-point EXCURSIONS** — TCV's null travels 1.6 m
in Z from its seed and E's flips sign in one step. **TCV's was the input**: its
filament arm carried its solenoid at 220 times its current, and corrected it is
a chronic plateau with the X-point within 8 cm of its seed →
[M-246](MEASUREMENTS.md#m-246). Two excursions, E and E shaped, remain.
→ [M-117](MEASUREMENTS.md#m-117) — the classification · the traces.

**AND THE SAME TRACES EXPLAIN WHY NO MERIT WEIGHT EVER TRANSFERRED BETWEEN
MACHINES.** The border constraints' share of `augmentedNorm` differs by three
orders of magnitude across two cases of one benchmark: on MAST `g·axis` alone is
**0.75** of a merit of 0.843, and on TCV the whole border is **~1e-03** against a
field residual of 2.953 — **a field residual measured on a TCV whose solenoid
carried 220 times its current** ([M-246](MEASUREMENTS.md#m-246)), so that end of
the range is not established; MAST's is. So MAST is the end where the border
already dominates — raising it cannot do anything, which is exactly what
`BorderMeritWeight`'s upward saturation measured — and TCV is the end where
nothing restrains the X-point at all. M-113's inversion, M-114's local optimum
and that saturation are one fact seen from opposite ends of a range, not three
findings ([M-113](MEASUREMENTS.md#m-113), [M-114](MEASUREMENTS.md#m-114)).

**AND E's DEFECT IS ROOT SELECTION, NOT GLOBALISATION.** The Picard pre-stage
and M-106's exact-seed run reach 1.0308e-01 and 1.0252e-01 — 0.5% apart, both
23% from the reference. Two routes sharing nothing find the same wrong branch,
so no step-length management is the repair there.

*From `CLAUDE_FB.md`, *The cold failures are three different things, and the
merit's balance is machine-dependent*. The `PicardSweeps` paragraph is under
*PicardSweeps and BorderedPicardThenNewton*.*

### Changing the fill, or blocking axis-touching elements, to keep the axis out of the support

**Declined.** On the iterates where an axis-touching element enters the support
transiently, **THE FILL IS NOT THE DEFECT AND CHANGING IT WOULD BE WRONG
TWICE.** At those iterates its answer is *correct* — `{Ψ > 0}` really is
connected and really does contain `R = 0` when `ψ_bnd < 0` — so a rule refusing
to say so would be a fill lying about the state it was handed. And blocking
axis-touching elements outright switches `p′` off in the one device class
`FREE-BOUNDARY-PLAN.md` §11.6 protects: a levitated dipole and a magnetic mirror
both carry plasma to `R = 0`, both have `gg′ ≡ 0` there, and the pole never
existed.

*From `CLAUDE_FB.md`, *§7.12b's fixture was the defect*.*

### Corner-only refinement for the ψ_ax layer

**Falsified.** The question was *"are the corner spike and the corner
FRAGILITY the same defect?"*, and its premise was that the spike sits in the
axis–`Γ` corner. **§11.1 measured that it does not**: the largest nodal values
are a flat layer along the **whole** `R = 0` boundary, agreeing to 3.5e-05 over
168 axis dofs, and which end of it wins the argmax changes with the polynomial
degree. So the corner is where the argmax lands, not where the phenomenon is,
and an experiment refining *only* the corner would have measured the wrong
thing.

What survived of it was the mechanism question: `PLASMA-EDGE-PLAN.md` §9.4
measured the extension losing order at a corner between two **transferred**
pieces; the axis–`Γ` corner is between a transferred piece and a **fitted** one,
and FB-1a reads **1.99 / 2.99 / 3.99** there — so on the face of it those are
different phenomena from this. What the layer *is* adjacent to is the thing
FB-A measured: the flux mass `( R q, v )` gives every element touching `R = 0` a
weight of order `h`, the smallest diagonal in the system, with an `O( 1/h )`
conditioning penalty.

**AND UNIFORM REFINEMENT DOES NOT REMOVE THE LAYER** — 1.0916e-01, 1.0982e-01,
1.0953e-01 at `n = 24, 32, 48`, against a true peak of 4.447e-02 and an imposed
datum of zero on that boundary — so it is not a discretisation error converging
away. FB-A's own study is the control that says this is not the axis *per se*:
there `ψ` reads 2.000 / 3.000 / 4.000 on a vacuum harmonic vanishing on the
axis, with no exterior coupling and no borders. The experiment it proposed —
drop the exterior coupling, then the limiter border, and see which the layer
follows — is answered by §11.3: a `1/R` pole in the load, which the limiter
border opens.

*From `FREE-BOUNDARY-PLAN.md` §11.3a; section numbers are that file's.*

### Narrowing ψ_ax's nodal maximum instead of locating the axis

**Declined.** Of three options for repairing what `ψ_ax` means, option 3 —
`AxisConstraint::LocatedAxis`, the located O-point — is built and is the
default. The other two, and the argument that nearly kept the nodal maximum:

1. **Keep the nodal max, exclude elements touching a fitted boundary.**
   Cheapest, would have prevented §7.12b, and is arbitrary — it names the
   symptom. **§11.1 sizes it**: on the half-disc that is the entire `R = 0`
   column, 168 dofs, not two corner elements — so the exclusion is a whole
   boundary layer and the arbitrariness is larger than the option made it sound.
2. **Restrict the argmax to the plasma component.** `meq::PlasmaComponent` and
   `refreshPlasmaComponent()`, which seeds off-axis, are the machinery. It
   couples `ψ_ax` to the support, which moves.
3. **Constrain `ψ_ax` at the located O-point.** The principled one. The
   objection to it was that it looks as though it costs the exact `−e_j` row,
   because the O-point moves with the solution — but the **envelope theorem**
   argument for why the nodal max needs no position derivative (`ψ_ax` is a
   *stationary* value, so `∇ψ = 0` at an interior extremum and the chain-rule
   term vanishes identically) transfers to a located axis, so option 3 is nearly
   free. It does **not** extend to an X-point, where the constraint is `q = 0`
   and the corner block is `∇q`; see §10.4.

**AND KEEPING THE NODAL MAXIMUM AND LETTING THE GUARD DO THE WORK WAS REVERSED.**
The resolution stood for part of a day as exactly that, on the argument that a
redefinition would prevent §7.16's defect and **hide** §7.12b's — the axis layer
would still be there and `ψ_ax` would simply stop reporting it. **That half is
true**, and it is why `checkAxisSource()` exists and why it refuses
independently. What the argument got wrong was treating the two as
alternatives: the redefinition is what makes `ψ_ax` mean what it says, and the
source guard is what catches a wrong field. Doing only the second leaves `ψ_ax`
defined as something that is not a magnetic axis, which is the defect §11.0
opens with.

*From `FREE-BOUNDARY-PLAN.md` §11.5, *Why this nearly was not done* and *The
three options, as originally written*; section numbers are that file's.*

### Initialising ψ_bnd from its own border row

**Declined.** `FREE-BOUNDARY-PLAN.md` §10.5's *initialise `ψ_bnd` from its own
border row* is **declined on its own measurement**, not pending: the limited
machine goes from 11 Newton steps to the 200 cap under it, and that item's own
conclusion (§10.8) is that this is a lever rather than a missing default.

*From `FREE-BOUNDARY-PLAN.md` §12's introduction.*

### Auxiliary globally-coupled unknowns

**Closed.** `FREE-BOUNDARY-PLAN.md` §6.4's FB-5 row — auxiliary unknowns carried
through the elimination — is closed because `NPCReduce()`/`NPCRecover()` are
public already and `DarcyNPCSolver::ArrayMult` blocks them →
[M-98](MEASUREMENTS.md#m-98).

*From `FREE-BOUNDARY-PLAN.md` §12's introduction.*

## Free boundary: geometry, conductors, meshing and the benchmarks

### The Lackner / von Hagenow Green's-function route

**Declined.** The classical route, and what MEQ's original free-boundary code
implemented; `attic/` held it and `git checkout 635aa3d -- attic/` brings it
back. Rejected for one reason, in Lackner's own words: the
precompute-and-reuse structure "will probably not be competitive with iteration
methods if the geometry of R is changed after each calculation". MEQ's stage 6
changes it every cycle.

**Anything from the original free-boundary code** is out of scope for the same
reason. The von Hagenow implementation there is the good algorithm badly
amortised, and the amortisation is what adaptivity destroys. Not worth
restoring — which is why `attic/` is gone rather than kept visible; it is at
`git checkout 635aa3d -- attic/` if the judgement is ever revisited, and
`refs/Refs.md` keeps the analysis of what it did and got wrong.

*From `FREE-BOUNDARY-PLAN.md` §2, *Why not the alternatives*, and §9.*

### A directly coupled BEM, and an outer fixed point between interior and exterior

**Declined.** The two other alternatives to coupling at a distance
(`refs/CouplingAtADistance.pdf`):

* **A directly coupled BEM with a kernel** — CEDRES++'s `c(·,·)`, eq (3.5), a
  double surface integral over `Γ` with a complete-elliptic-integral kernel. It
  is the right operator; `FREE-BOUNDARY-PLAN.md` §3 says why MEQ should not
  assemble it that way: the semicircle collapses it to a diagonal operator.
* **An outer fixed point** between interior and exterior — the paper's own §4.3,
  a Richardson iteration. Rejected for the reason MEQ rejects every outer fixed
  point: it is outside the residual, so the Jacobian cannot see it. The
  measurement that settled that argument is in `HighBetaConvergence` and it was
  not close — coupled 8.3e-02 → 4.4e-15 in four steps against decoupled
  8.3e-02 → 8.2e-02 in fifteen. **It does not converge slowly; it does not
  move.**

*From `FREE-BOUNDARY-PLAN.md` §2, *Why not the alternatives*.*

### Lost coverage as the cause of the red tiling, and relaxing the gate

**Withdrawn.** Filed and answered the same day. At a 12th-order face rule the
boundary weights summed to `|Γ|` only to 1.01e-04, 2.24e-05, 4.59e-06 at
`n = 12, 24, 48`, converging at about `O(h²)` where the expectation is a
mesh-independent floor near 1e-10. A one-variable control — the same `D_h` with
`VertexConePath`'s cone off — put it back at the floor, so the cone was
unambiguously the cause. MEQ filed that as **lost coverage**, arguing that *"a
quadrature residual that converges is measuring a geometry rather than an
instrument"*.

**THE ARGUMENT WAS WRONG AND THE HOLE IN IT IS WORTH REMEMBERING: TWO THINGS
CONVERGE.** The second is a curve the *rule* under-resolves, which straightens as
`h` falls. Refining the **rule** at fixed `h` separates them, because no
quadrature recovers coverage that is not there. Measured on MEQ's own circle at
`n = 12`, cone on ([M-12](MEASUREMENTS.md#m-12)):

| q8 | q12 | q20 | q40 | q80 |
|---|---|---|---|---|
| 2.61e-04 | 1.01e-04 | 1.34e-05 | 6.68e-08 | **5.40e-10** |

`5.40e-10` is the cone-off floor. **Coverage is exact either way.** What the cone
costs is the smoothness of `ξ ↦ a(x(ξ))` — it drives the two interpolated vertex
directions apart, the foot map roughens, and a fixed-order rule under-resolves
it.

**So the repair was the rule, not the gate.** The gate was never too tight. That
distinction is the whole finding: relaxing it would have made the suite green and
hidden a genuine under-resolution in MEQ's own transmission row, whose
`transmissionQuadratureOrder` was the same 12 and is 40.

**And the acceptance beside it could not have caught it**, which is worth stating
because it reads 7.8e-16.
`theTransmissionRowIsTheBoundaryIntegralItClaims` builds its reference by
sweeping `Γ` with the **same rule the row uses**, so the quadrature error is
common to both sides and cancels exactly. That is deliberate — it is what
isolates the contraction, the vdof ordering, the sign and the measure — but it
makes the case blind to whether the rule resolves the foot map at all.

**What upstream did with it**, all of it useful: reproduced it on their own disc,
made the cone a `use_cone` constructor flag and turned it **off by default**,
added the boundary-sweep tiling case MEQ asked for, and corrected their own
commit message's *"changes nothing"*. They also measured that the cone does not
do what it was added for — the aerofoil flux order it targeted is
**pre-asymptotic** and recovers on its own — and connected it to MEQ's earlier
signed-weight finding: **the cone is what makes the foot map backtrack**, so the
`O(h)` unsigned overcount and this `O(h²)` residual are two readings of one
thing.

*From `FREE-BOUNDARY-PLAN.md` §7.7, which keeps its heading and the verdict.
`CLAUDE_FB.md`, *At `j = 0` the question does not arise*, carries the same
account interleaved with what it moved in MEQ's code, and keeps it.*

### A closed-form conductor field as more accurate than one assembled on the mesh

**Withdrawn.** FB-7's first draft claimed a closed-form conductor field beats one
assembled on the mesh. **What FB-7 buys is DOMAIN REDUCTION, not fidelity**, and
the claim is wrong twice over.

**First, because a real coil has finite extent and finite current density, and
near the plasma that matters to the solution.** A filament is an *idealisation*:
evaluating it exactly still solves the wrong problem when the conductor is close
enough for its multipole content to reach the plasma. Exactness of evaluation
and adequacy of the model are different questions, and only the second decides
whether an answer is right.

**Second, because MEQ's coils were never filaments.** `meq::Coil` is *"one coil
of rectangular cross-section, carrying a prescribed total current uniformly
distributed over it"*, and `coilPsi()` integrates the Green's function over that
footprint with a graded Gauss rule — 32 points per direction by default, machine
precision outside the coil and about 1e-12 inside. So **finite extent and finite
current density are already modelled**, and the exterior route inherits that
fidelity rather than trading it away. The filament is the limiting case, not the
method.

*From `FREE-BOUNDARY-PLAN.md` §7.19, *IT IS NOT MORE ACCURATE, AND THE FIRST
DRAFT OF THIS SECTION SAID IT WAS*. `CLAUDE_FB.md`, *FB-7*, states the same
verdict in one paragraph and keeps it.*

### Triangle, CGAL Mesh_2 and MMG as the mesh generator

**Declined.** **What was NOT chosen, and why the survey was short.** Triangle is
2D-only and has no subdomain tagging worth the name; CGAL's `Mesh_2` is a large
dependency for a 2D job; MMG is a re-mesher rather than a generator and would sit
*after* gmsh rather than instead of it — worth revisiting if FB-5's adaptivity
ever wants to move the geometry between cycles. Gmsh reads a `.geo` script,
drives from Python or C++, and is what MFEM's own examples assume.

*From `CLAUDE_FB.md`, *Meshing beyond `MakeCartesian2D`: gmsh, and it is already
here*.*

### Cold-start guesses: the vacuum field alone, a filament column, a rectangular column

**Falsified** (the vacuum field alone, the filament) **and Declined** (the
rectangle). The driver's free-boundary guess is the conductors' own vacuum field
plus `I_p` over an ELLIPTICAL column about a guessed axis.

**THE VACUUM FIELD ALONE CONVERGES TO A BRANCH THAT IS NOT ONE**, which is the
measurement that makes the second part necessary rather than an improvement:
MAST-U converges in 9 iterations to `psi_ax` 3.5× the reference, with its axis
at ( 2.65, −0.91 ) — outside the machine.

**AND THE COLUMN IS AN ELLIPSE BECAUSE OF WHAT THE OTHER TWO SHAPES DO.** A
FILAMENT is singular at exactly the point `locateAxisPoint()` has to find —
measured, 6.667e-01 at 3 mm from the guessed axis against MAST-U's reference
`psi_axis` of 9.187e-02, and still climbing — and the solve fails. A RECTANGLE
is bounded, `meq::Coil` already is one, and a uniform current density over it
carries a logarithm in its second derivatives at each of four corners; those are
artefacts of the shape rather than anything the equilibrium puts there.
→ [M-124](MEASUREMENTS.md#m-124).

*From `CLAUDE_FB.md`, *MEQ solves MAST-U from a configuration file alone*.*

### η₆ in the marking

**Falsified.** What this copy adds is the measurement that raised it and the
proposal itself. The adaptive loop cannot see the exterior coupling, and `η` is
not wrong. Measured on the half-disc at `k = 2` over four cycles, with `Γ` fixed
and only `Γ_h` refining:

```
    cycle   elements   Gamma_h faces        eta      L2(psi)   |a - exact|
      0        314          34          2.5109e-01  3.5828e-03  1.3194e-03
      3        640          34          3.2079e-02  9.3347e-04  1.3194e-03
```

`η` falls eightfold, the interior error follows it, and the exterior
coefficients **do not move** — five digits, four cycles. Under **uniform**
refinement the same quantity converges at 3.32, so it is not a mode-truncation
floor.

**It is not a marking accident.** The 34 elements touching `Γ_h` are 11% of the
mesh and carry 0.00% of `η²`, rising only to 0.07% by cycle 3 as the interior
improves around them. Their indicator is essentially zero, so no threshold
whatever would mark them, and `Γ_h` keeps its 34 faces while the element count
doubles.

**And `η` is answering a different question.** It estimates the **interior**
discretisation error, and `η₅` on `Γ_h` compares `ψ*` against the datum actually
imposed — the repair `CLAUDE.md` records under *A separate `η₅` problem on the
extension path* — so it correctly reports the boundary as well resolved **for
the interior problem**. The exterior coefficients are a **boundary
functional**: a transmission integral over `Γ` reached by extension from `Γ_h`.
Nothing in `η` measures it. Both statements are true at once and neither is a
defect. Left alone, the second becomes the floor of the first within a few
cycles — the loop reporting a falling `η` while `ψ` stops improving.

**The proposal, which is the falsified half:** take the transmission residual
per face of `Γ_h` — the same integral `exteriorTransmissionRows()` sweeps, kept
per face instead of contracted and summed — and add it to the marking beside
`η`. **Added to `η` as a sixth term it does nothing**: `η₆` is 8.58e-04 against
an `η` of 2.51e-01, so its share of `η²` is about 1e-5 and a Dörfler competition
never reaches it — `Γ_h` kept all 34 faces for four cycles with the term summed
in. The two are different quantities in different units, so one sum over both
is a comparison with no meaning however it is weighted. **What works is a second
marking pass** on the boundary term's own distribution, unioned with the
interior one; `η₆` stays in `η` so the **stopping** rule sees it. Measured,
`Γ_h` then refines 34 → 46 → 57 → 64 and `|a − exact|` falls 1.32e-03 →
1.53e-04 against a control that does not move → [M-05](MEASUREMENTS.md#m-05).
`theCoupledSolveSurvivesTheAdaptiveLoop` is the control, and `CLAUDE_FB.md`,
FB-5's section, the write-up.

*From `TODO` *A boundary error estimator for the exterior coupling*.*

## Plasma edge

### Cut quadrature, Algoim and moment fitting

**Declined** (cut quadrature for the order, `MFEM_USE_ALGOIM`) **and Falsified**
(`MomentFittingIntRules` as robust). **NO CUT QUADRATURE WAS BUILT *FOR THE
ORDER*, AND THAT IS A DECISION RATHER THAN AN OMISSION.** An exact cut rule would
move `ψ_h` from `j+1.5` to `j+2.5` and `ψ*` by nothing that a higher Gauss order
does not already buy; it would buy `q_h` **nothing at all**, `q` being at its own
regularity bound already; and it would not move the `k ≤ j` threshold. **What it
WOULD buy is `j = 0` at all** (`CLAUDE_FB.md`, *At `j = 0` the question does not
arise*), so "a cut rule buys nothing" is too strong. `refs/CutElementQuadratureSurvey.pdf`
(Loibl et al., arXiv:2602.18130) is the survey of the ten open-source
implementations and it adds two practical reasons: every one of its benchmarks
is Cartesian, and **MFEM's two cut backends are quadrilateral and hexahedral
only** — the Algoim path aborts with *"supports only quads and hexes"* and
`MomentFittingIntRules` builds four-vertex local meshes — where MEQ's meshes are
triangles. Turning `MFEM_USE_ALGOIM` on would therefore cost a mesh-type change
as well, for a threshold that does not move.

**`MomentFittingIntRules` WAS TRIED AND IT WORKS, WITH A CAVEAT WORTH KEEPING.**
On quadrilaterals it is exact on a straight cut (9e-16) and reached 3.3e-06 on a
disc — but on one mesh in three it produced weights of **−30 and −63 on elements
of area 1e-4**, and the area came out 1.3e-02 wrong by cancellation. That is the
known conditioning fragility of moment fitting on a nearly degenerate cut, and
it is why the survey's other family — dimension reduction, which is Algoim — is
the one that reaches its design order.

**AND THE FIRST ATTEMPT AT IT READ ZERO FROM EVERY CUT ELEMENT**, because it
held `mesh.GetElementTransformation( int )` across the call and MFEM's own cut
code resets that shared scratch underneath it. Third time this trap has been
recorded in this tree — `CLAUDE_MFEM.md`, *Traps in MFEM's API*, is where it lives; the fix is a
local `IsoparametricTransformation`, which is what `ex38.cpp` does.

*From `CLAUDE_FB.md`, *The plasma edge caps the order, and it is the PROFILE
that sets the cap*.*

**Declined** for the plasma-edge plan. **It is not a cut quadrature and does not
become one.** If what is wanted is *solvability* at `j = 0` rather than high
order, a cut volume rule plus the matching cut surface rule for the Jacobian's
`∮ F φ/|∇ψ|` term is a much smaller piece of work with a much smaller prize —
about second order — and it should be costed separately rather than folded into
the interior-interface plan.

*From `PLASMA-EDGE-PLAN.md` §8.*

### Locating the plasma edge with `ψ*` rather than `ψ_h`

**Closed: measured and not worth building** → **[M-94](MEASUREMENTS.md#m-94)**.
It cost no new code because the oracle was already in the tree.
`tests/analytic/PlasmaEdge.hpp` carries two fixtures over one exact solution:
`MovingPlasmaEdge` reads its support off `ψ_h`, which is what MEQ does, and
`PlasmaEdge` cuts at the **exact** edge. Since `ψ*` is a better approximation to
`ψ` than `ψ_h` is, and the fixed cut *is* `ψ`, the second column is an **upper
bound on everything this item could buy**. `locatingTheEdgeExactlyBuysNoOrder`
is the case.

**No order is on the table.** If the level set limited the rate, the oracle's
advantage would grow as `h` fell; over a fourfold refinement it grows by at most
**1.018** in `ψ_h` and **1.121** in `ψ*`, and in the one row where the oracle is
clearly ahead it *falls*, 2.98 → 2.03.

**And in the regime MEQ runs in it buys nothing at all.** At `k ≤ j`, where `ψ*`
keeps `k+2`, cutting at the exact edge is worth between **0.95× and 1.28×** — the
moving edge is the better of the two in half those rows. The one corner with a
real factor is `j = 2, k = 3`, about **2×**, which is `k > j`: the regime where
the cut already caps the order, which is what `PLASMA-EDGE-PLAN.md` addresses by
a different route, and where `ψ*` would capture only the part of the gap lying
between `ψ_h` and `ψ`.

**The premise was wrong, and it is kept because it is the reasoning not to
repeat.** It read: *"The place `ψ`'s accuracy DOES bind is the level set.
`ConfineToPlasma` locates `{ Ψ > 0 }` as a level set of `ψ_h`, and
`PLASMA-EDGE-PLAN.md`'s rate cap is about the geometry of that set rather than
about the load. Locating the edge with `ψ*` while still assembling `F` at `ψ_h`
is independent of 12.1, needs no new analysis and no Jacobian change, and has
`PlasmaEdgeConvergence` waiting as its acceptance criterion. It is an
afternoon."* It does not bind there either. What limits the cut is the **rule's
blindness** — a quadrature rule cannot see a kink between its points, whatever
field decided where the kink is — which
`theLossIsTheRulesBlindnessAndNotItsResolution` measures separately and which no
improvement to the level set can touch.

**It was an afternoon, and it was spent on the measurement instead of on the
build — which is the cheaper order whenever an oracle exists.** Bounding what a
change can buy before writing it needed two fixtures that were already there and
one test case; writing it would have needed a second support field through
`meq::SourceIntegrator` and a source-side way to evaluate `F` at one `ψ` while
testing the support at another, which is an API change across all three
`NormalisedSource` subclasses.

**Nor `ψ*` in the source for its own sake.** Substituting `ψ*` into the source
with the quadrature left in place is not worth it: the consistency error of
`F( ψ_h )` is `O( h^(k+1) )`, the method's own rate, so `ψ*`'s `O( h^(k+2) )`
removes a term that limits nothing and buys a constant. Its value in
interpolatory HDG (`ROADMAP.md` item 12.1, `INTERPOLATORY-HDG-PLAN.md`) is
structural — it is what the interpolation needs to stay superconvergent — and
not that the source was inaccurate.

*From `ROADMAP.md` item 12.2.*

### Fitting the mesh to `Γ_p`, or enriching the space (XFEM/GFEM)

**Declined.** Three routes beat an unfitted polynomial cap, and all three work by
putting the interface into the discretisation. The plasma-edge plan takes the
third — two subdomains coupled at a distance — and not these two:

| | what it is | why not |
|---|---|---|
| **Fit the mesh** | remesh so element faces follow `Γ_p`, curved | `Γ_p` moves every iterate, so it is a remesh per iterate; the trace space changes shape; and a *straight-sided* fit caps at second order, which `refs/CutElementQuadratureSurvey.pdf` makes its central finding |
| **Enrich the space** | add `d_+^{j+2}` to `P_k` on cut elements (XFEM/GFEM) | cheapest in principle, and the exponent is known exactly. But HDG's `k+2` post-processing is a theorem about polynomial local problems, and it would have to be re-derived on an enriched space. A large piece of analysis MEQ cannot check by measurement alone |

*From `PLASMA-EDGE-PLAN.md` §2.*

### The interior-interface plan's scope: `j` stays in the physics, and the fixed-boundary path gains nothing

**Out of scope.** **It does not remove `j` from the physics.** A profile with
`p'(0) ≠ 0` still puts a jump in the current density at the plasma edge; the
plan represents that jump instead of smearing it. Whether a user wants such an
equilibrium is a modelling question and this changes nothing about it.

**It does not help the fixed-boundary path.** There `Γ` is given and stage 5
already does this.

*From `PLASMA-EDGE-PLAN.md` §8. That section's last item — not to start while
`j ≥ 1` is unfinished — is dropped: the precondition is met, and the file's own
header says so.*

### A single global Fourier modal border on a cornered `Γ_p`

**Declined.** **`PLASMA-EDGE-PLAN.md` §3 expands `λ( s )` in one curve parameter
and claims spectral convergence with `N ~ log( 1/h )`. That claim comes from
`CouplingAtADistance.pdf`, whose `Γ` is an ARTIFICIAL interface chosen smooth** —
a circle, in its §4.2 — and which says so: *"the BEM is defined on a suitably
chosen, **smooth** artificial interface"*, *"let `x : R → Γ` be a **smooth
regular**"* parametrisation. They pick a smooth `Γ` precisely to avoid this
question. It is the **extension** work, `HDG-CurvedExtensions.pdf` and the
Joukowsky aerofoil of its §3.4, that is built for a sharp `Γ` — and
`VertexConePath` exists because both closed-form families give out there.

**So the transfer half of the plan has a parent that meets corners and the
modal half does not**, and a global Fourier series on a closed curve with a
corner converges at `O( 1/m² )`, not spectrally: `N ~ h⁻¹`, and §3's whole cost
argument against a mortar space collapses. The per-arc basis that replaces it
stays in `PLASMA-EDGE-PLAN.md` §9.6.

*From `PLASMA-EDGE-PLAN.md` §9.6.*

## Inverse solve

### An outer loop round the forward solve for the coil currents

**Declined.** **MEQ has already built an outer loop around this solver and measured it.**
`CLAUDE_INVERSION.md`, *Driving by `q(ψ)`*: the `q(ψ)`-driven solve is a KINSOL
outer Newton whose unknowns are the coefficients of a fitted `g²`, wrapping the
bordered Newton. `DriverAcceptance::theDriverSolvesForTheToroidalField` reads
**42.5 s, 12 outer iterations, 43 equilibria** to recover a three-coefficient
closed form.

**43 equilibria for 3 outer unknowns.** A differenced outer Jacobian costs one
equilibrium per unknown per outer step, so at DIII-D's **18 coils** — and
`examples/machine-*.toml` carries 4, 4, 11, 24, 4, 18 and 26 — the same
architecture costs of order **19 equilibrium solves per outer step**, against
`N_free` **backsolves** per Newton step for the bordered route. The two are not
close, and the gap widens with the coil count.

**And that is only the cost argument.** The house objection is stronger and is
about the answer rather than the clock: `CLAUDE.md` records that
`PicardThenNewton` and `Globalisation` are not exposed as keys because three
solve routes reach discrete solutions **9.4% apart** on an under-resolved mesh,
and `CLAUDE_INVERSION.md` records that *"a converged outer solve is not evidence
the root is the answer"* — both failing configurations of the `q(ψ)` loop
converged, to 3.5e-07, on a fixed point that was not the answer. **An outer loop
adds a second convergence criterion, and a second thing that can converge to the
wrong place quietly.**

**The `q(ψ)` loop is an outer loop for a reason that does not apply here**, and
saying so is the honest version of this argument: `g²` enters the source through
`setGGPrime()`, a whole profile object, and there is no assembled route to
`∂R/∂(spline coefficient)`. **There is an assembled route to `∂R/∂I_c` — it is a
constant vector**, `INVERSE-AND-RECONSTRUCTION-PLAN.md` §3.3. That is precisely
the difference between a quantity that belongs in a border and one that does
not.

**What an outer loop would buy, and it is not nothing.** It would reuse
`freegs4e`'s algorithm unchanged, need no change to `solveWithNormalisation()`,
and let the regularisation be somebody else's well-tested least squares.
`INVERSE-AND-RECONSTRUCTION-PLAN.md`'s **IV-0** is exactly that, as a throwaway
instrument rather than as the product.

*From `INVERSE-AND-RECONSTRUCTION-PLAN.md` §4.2.*

### Machinery right only for reachable targets: one least squares over every row, row four with the iterate's `H`, and the augmented-norm merit

**Reverted**, [M-252](MEASUREMENTS.md#m-252). IV-5 falsified three things IV-3
and IV-4 built, each right only for reachable targets:

* one least squares over every row of the corner, which traded the
  equilibrium's equations against unreachable targets on the first, unsettled
  sweep;
* row four with the iterate's `H`, whose dropped term made a best fit creep
  uphill;
* the augmented norm as the inverse line search's merit, which refused every
  step when one current faced 23 unreachable rows.

The corner is an equality-constrained least squares, row four is the exact
displacement, and the inverse line search judges an exact-penalty merit.

*From `INVERSE-AND-RECONSTRUCTION-PLAN.md` IV-5.*

## Coil subtraction

### CS-5's acceptance against `freegs4e`

**Falsified**, [M-139](MEASUREMENTS.md#m-139). **CS-5 as staged cannot work and
M-139 is why.** The benchmark's floor after matching the conductor models is
`7.7e-04` and it does not refine, because `freegs4e.shaped_coil.ShapedCoil` caps
at 6 points per triangle. **A change whose entire claim is that the conductors
are resolved exactly cannot be accepted against a reference that resolves them
to 6 points.** The replacement is MEQ's own quadrature-based reference,
[M-149](MEASUREMENTS.md#m-149), and it stays in `COIL-SUBTRACTION-PLAN.md` §7.2.

*From `COIL-SUBTRACTION-PLAN.md` §7.2.*

### A total `GridFunction` `psi_p + psi_c` for consumers to read

**Declined.** **The rejected option is to form a total `GridFunction` once and
let consumers read that.** It is the smallest diff and it throws away the plan's
central claim: `psi_c` projected onto the discrete space is `psi_c` to the
mesh's order, which is exactly the error `COIL-SUBTRACTION-PLAN.md` §2(b) exists
to remove. It would also reintroduce the singularity into a finite-element
function, which for a filament is the thing that cannot be represented at all.

*From `COIL-SUBTRACTION-PLAN.md` §8.2.*

## Solution inversion

### Open-surface averages and an X-point-aware tracer

**Out of scope.** There are no flux-surface AVERAGES over an open surface, and
there should not be: `V′` and `⟨R^{-2}⟩` are integrals round a closed loop, and
the volume an open curve encloses is not defined. What IN-5 delivers is the
GEOMETRY — the curve, its arc length and a basis that represents it — which is
what a scrape-off-layer consumer asks for, and it is the scope
`INVERSION-PLAN.md` §6 describes. Nor is the tracer X-point aware: a level AT a
separatrix stalls at the saddle, where the level set is not a 1-manifold, and
`Stalled` is the honest answer rather than a defect.

**What is deliberately absent from the tracer**, so nobody reads more into it
than is there: it follows **one connected component** and neither finds nor
reports a disjoint island at the same level — the same disclaimer
`CriticalPoints.hpp` makes about seeded Newton, and for the same reason. It is
**not X-point aware**: the level set through a saddle is not a 1-manifold and
the tangent is undefined there, so a trace at that level will stall or turn a
corner arbitrarily. `pointAtArcLength()` parametrises segments linearly in
*polyline* length rather than true arc length, and says so.

*From `CLAUDE_INVERSION.md` IN-5 (WHAT IS DELIBERATELY ABSENT) and What is next
(What is deliberately absent from the tracer), and `INVERSION-PLAN.md` IN-5
(WHAT IS NOT BUILT); three copies of one list, merged.*

### Implicit quadrature on the level set as IN-2's third cross-check

**Declined.** Two of the three cross-checks were built and the third was not.
Ray bisection from the axis (`INVERSION-PLAN.md` §3.4) is kept as the second
extraction; `INVERSION-PLAN.md` §3.3's implicit quadrature — a rule on the level
set with no curve extracted at all — is the third leg, and it is declined rather
than outstanding, so nobody reads it as work waiting. Two measurements close it.
IN-P puts the quadrature at **1.8%** of the extraction chain against
`gaugeFreeFit`'s 84.6%, so an implicit rule saves essentially nothing on the leg
it replaces; and two independent extractions already disagree by less than the
discretisation — the ray-bisection control agrees with the tracer, and the
displaced-origin check moves `V′` by **3.6e-10**. A third route would be a third
implementation to maintain in exchange for a check the second one already
passes. If a reason to build it ever appears it will be a geometry the tracer
cannot close, not a cost.

*From `INVERSION-PLAN.md` IN-2.*

### Other ways to extract a contour: marching, subdivision, rational contours, unprojected ODEs, rays, Morse–Smale

**Declined**, all six; ray bisection is kept as a cross-check and not as the
primary route.

| | |
|---|---|
| **Marching squares / marching triangles** | Caps at `O(h_cell²)` regardless of `k` — Etiene et al. measured vertex position 1.94, normals 0.93, curvature *divergent* (`10.1109/tvcg.2009.194`). On a simplicial mesh it is at least unambiguous, the affine interpolant's level set being convex hence connected; but the order is the point and it throws away `k−1` of them. |
| **Subdivide, then march** — what ParaView, VisIt, VTK and MFEM's own exporters do | **Second order in the SUB-CELL size with no `k`-dependence at all**; Remacle et al.'s error falls as `2^{2r}` and the degree never enters. Fine for a picture, which is what it is for. Not fine for `q(ψ)`. |
| **Exact rational contours** | Wiley et al. give the contour of a bivariate quadratic over a triangle as an exact conic. **It caps at `k = 2` in principle**: a level set of a degree-`d` polynomial is an algebraic curve of degree `d`, and rational parametrisation needs genus 0 — a generic smooth plane cubic has genus 1. MEQ runs to `k = 4`. (Their *quartic* comes from curved element geometry; MEQ's solve mesh is straight-sided, so it would be a conic here.) |
| **Plain ODE tracing without projection** | Error accumulates monotonically; the curve drifts onto a neighbouring contour and does not close. With projection it *is* `INVERSION-PLAN.md` §3.1, the predictor–corrector continuation MEQ uses. |
| **Ray bisection from the axis** | Point accuracy is as good as a corrector's and it parallelises trivially — worth keeping as an **independent cross-check** and as the initialiser for a periodic parametrisation, since it produces points already indexed by angle. **Not the primary route**: it assumes star-shapedness, which fails on indented cross-sections, at and beyond the separatrix, and near the axis where the bracket degenerates. |
| **The Morse–Smale complex** | **It gives the wrong curves.** Its "separatrices" are integral curves of `∇ψ`, *orthogonal* to the level sets; the plasma separatrix is the level **set** through the X-point. Recorded because it is the natural thing to reach for by name. |

*From `INVERSION-PLAN.md` §3.4, What is rejected, and why (the heading stays
there as a stub, because `src/meq/FluxSurfaces.{hpp,cpp}`,
`tests/analytic/FluxSurfaceReference.hpp` and
`tests/convergence/FluxSurfaceConvergence.cpp` cite it for ray bisection).*

### A kernel survey: the Kintner recurrence, Vandermonde + GEMM, and batched field evaluation

**Closed**, three ways: two not worth it, one worth it only for a consumer that
does not exist yet.

→ **[M-51](MEASUREMENTS.md#m-51)** — measured: the Zernike basis by Kintner
recurrence against per-mode evaluation is **6.3×**, agreeing to 4.4e-15, and
**not worth it**, being 11.5% of a 0.6% stage; a fit at 20,000 points by
Vandermonde + GEMM against per-point is **34.7×**, agreeing to 4.4e-16 m, worth
it for the many-point case and for nothing today (see *Vandermonde + GEMM
extraction*); batched field evaluation by element is **0.90–1.14×**, and **do
not**.

**And the reason the third one fails is worth keeping**, because the brief
asserted the opposite: `GridFunction::GetValues( i, ir, vals )` is **a loop
calling `CalcShape` per point, not a GEMM**. It amortises the dof gather and the
temporaries and nothing else. The access pattern is 1.15 points per element in
any case, so there is nothing to batch.

*From `CLAUDE_INVERSION.md`, Three things in the kernel survey worth not
re-doing.*

## Toroidal flow and the mirror

### XP-1's element flood fill on the mirror path

**Declined**, dissolved by the axial cutoff (`MIRROR-PLAN.md` item 2). The
element-level half of the support takes the connected component containing the
magnetic axis. A mirror has no axis, and a flux-labelled annulus is connected
*through the tube* to the ends of the machine — so the fill neither starts nor
terminates where it should.

**THE AXIAL CUTOFF REMOVES THE QUESTION.** With `|z| ≤ z_edge` the support is a
bounded annular box, defined pointwise and needing no seed and no connected
component. So the fill is simply not wanted on this path — which is a decision
the plan can take rather than a `2 × 2` it has to run. What stays live in
`MIRROR-PLAN.md` is that turning it off must be a **stated** branch a test
asserts, not a flag that happens to be false.

*From `MIRROR-PLAN.md` §3.*

### A plasma-current border for the mirror

**Declined.** `CLAUDE_FB.md` records that with a fixed profile amplitude a
confined equilibrium is a non-linear *eigenvalue* problem and the coupled solve
fails from everywhere without `setPlasmaCurrent`. That argument turns on the
support being free. On the mirror **both** labels are pinned at prescribed
points, so the support is not free and the eigenvalue does not arise — geq
prescribes physical densities and temperatures and never a current. Expect to
need no current border, and treat needing one as a finding.

*From `MIRROR-PLAN.md`, What does NOT bite, which is worth saying.*

## Interpolatory HDG

### The interpolatory source: seven hypotheses rejected

**Declined.** Considered for the CCSZ interpolatory port and rejected, so they
are not re-derived. Section numbers are `INTERPOLATORY-HDG-PLAN.md`'s.

* **"Interpolate `F( ψ_h )` into `W_h`."** That is the 2019 predecessor, and it
  loses the superconvergence of the postprocessed solution — which for MEQ means
  `ψ*` at `k+1`, the estimator a full order worse, and the adaptive loop
  refining on a degraded indicator. `Z_h = P^{k+1}` is not a refinement of the
  idea, it is the theorem.
* **"Feed `DarcyForm::Reconstruct()`'s `ψ*` to the source, since MEQ already has
  it."** It lifts the non-linear potential integrators and takes their gradient
  at the computed potential, so it is not linear in the unknowns, `B11`/`B12`
  would not be constant, and each element's reconstruction would become an
  implicit local fixed point. §6.4.
* **"Interpolate `F/R`, since that is what the weak form carries."** A pole at
  `R = 0` on the half-disc, and it destroys the Solov'ev exactness. §1.6, §7.1.
* **"The saving is that quadrature of `F` disappears."** At MEQ's degrees
  `ns` is within a factor of two of `N_q`, and at `extraOrder = 0` it is
  **larger**. §3.2. The paper's framing is a time-stepping one — matrices
  assembled once before the integration — and MEQ already reassembles nothing
  but the source.
* **"It removes `ComputeH()`'s cost."** It does not; `ComputeH()` is 68–73% of
  the gradient leg on every row of [M-80](MEASUREMENTS.md#m-80) and is
  element-local dense algebra, not assembly. §3.3, §3.4.
* **"Paper I measures the speed-up."** It does not. Grepped: the paper reports
  Table 1's convergence history and pattern-formation figures, and **no timing
  of any kind**. Remark 2.3's claim of computational advantage is asserted. The
  only measurement of the gap this project has access to is MFEM's own, in
  `reaction_hdg.hpp`'s doxygen.
* **"MEQ's `AssemblyMode` or `TraceSolver` precedent means `SourceTerm` can be a
  TOML key."** Those cannot change the answer; this one changes it by
  construction. IH-8.

*From `INTERPOLATORY-HDG-PLAN.md` §10.*

## Output and sampling

### Smoothing the boundary displacement on the vdof graph, and a global backoff

**Falsified.**

**Two things about bending the `.vtu` boundary onto `Γ` were got wrong on the
way and are worth not repeating.** Smoothing the displacement into the interior on the *vdof* graph
couples R to Z — they share one index range — so a radial displacement gets
averaged against a vertical one; measured, that made tangling **worse**, 25%
surviving where moving the boundary alone managed 50%. And backing the
displacement off **globally** costs every face the worst face's limit: per-node
backoff takes `miller-curved` from 50% to **96%** of the boundary reaching `Γ`.

*From `CLAUDE.md`, “Status: the solver works, and MEQ is a program”.*

### Rebuilding the ParaView collection per frame

**Falsified.**

**One frame per cycle is easy to get ALMOST right, and the near miss is
invisible.** Rebuilding the collection per frame puts every `Cycle` directory on
disk with the correct refined mesh, and leaves the `.pvd` index listing only the
last of them — so ParaView opens the file and shows a single frame, with no
error and all the data present. `ParaViewDataCollection` appends to its `.pvd`
and does not scan the directory, so the collection has to survive between frames
and be rebound with `SetMesh()`. `theAdaptiveSeriesIndexesEveryFrame` asserts on
the **index**, not the pieces, because the pieces were never what broke.

*From `CLAUDE.md`, “Status: the solver works, and MEQ is a program”.*

### Pinning ∇q in the band

**Declined, recorded so nobody re-derives it.** `B` in the band is `O(h²)` at
every `k` and does not reach `ψ`'s order, because `ψ` is continued with a
*solved* variable and `∇q` is a *differentiated* one. `div q = −F/R` and
`∂_r q_z − ∂_z q_r = −q_z/R` pin two of `∇q`'s four entries exactly, but leave
the symmetric traceless part differentiated — structure rather than an order,
and it would need the source plumbed into `GridSampler`. `CLAUDE_DRIVER.md`,
*MEQ is runnable*, has the detail.

*From `ROADMAP.md` item 11.*

### Face averaging in the grid sampler: gslib's SetL2AvgType

**Falsified: it is a real rule and it is the wrong one** →
[M-93](MEASUREMENTS.md#m-93); `theTwoLocatorsAgreeAtTheSamePoints` in
`tests/convergence/SamplerConvergence.cpp` is the case.

**The proposal.** gslib has a principled answer where `meq::GridSampler` has a
shrug. A grid node landing exactly on an inter-element face is ambiguous, `ψ_h`
being discontinuous there: whichever element claims it first wins. The
disagreement is the face jump, `O( h^(k+1) )`, so it converges away and nothing
is wrong — but it is an arbitrary choice where a rule exists.
`FindPointsGSLIB::SetL2AvgType` does arithmetic or harmonic averaging over the
elements meeting at a point, for exactly this case, and both of MEQ's sampled
fields are L2. The route to it is not free: `meq::FieldTransfer` is
`meq::PointLocator`, MEQ's own, and `MFEM_USE_GSLIB` is optional — it has to be,
because gslib carries no floating-point abstraction and cannot be built beside
`MFEM_PRECISION = single` — so taking gslib's rule means writing it in
`meq::PointLocator` rather than calling a setter. The sampler does not otherwise
want gslib: the output grid is uniform, so locating a point is index arithmetic
and the search is already free.

**What it measured.** On `psi_h` at `k = 2`, with gslib's default averaging left
on, **all 577 border nodes** of a 65² grid disagree with `GridSampler`, worst
1.93e-06, because the averaged value belongs to **neither** element's
polynomial — which is the same finding `meq::FieldTransfer`'s constructor
records from the other direction, at 28%, 10% and 12% wrong in L2 at
`k = 1, 2, 3`. Picking a side, which is what `GridSampler` does, is the better
rule for a **discontinuous** field. The shrug was right and nothing is to be
wired.

**The comparison the entry also asked for agrees exactly where the question has
an answer.** The two samplers had never been compared to each other —
`WarmStartConvergence.cpp` measured `FieldTransfer` against the exact solution,
not against `GridSampler`. **3648 of 3648** nodes strictly inside an element
agree to 7.2e-16 — two implementations sharing no locating machinery, which is
the check an exact-solution comparison structurally cannot make. Every
disagreement is a border node, where both are right and the difference is the
face jump, 1.17e-06 against 2.03e-06 of error in `psi_h` itself. Border nodes are
**common** — 577 of 4225, one in seven, on a grid with no alignment intended,
because a uniform grid over a uniform mesh is commensurate with it far more
often than a generic point would be. And the count of nodes the two resolve to
**different** elements is not a stable number: one ulp in the mesh extent,
0.7999999999999999 against 0.8, moves it from 202 to 331. It is reported and
never gated.

**And the real defect was in neither question.** Setting gslib up on a mesh
calls `Mesh::EnsureNodes()`, and `GridSampler`'s affine fast path keyed on
`Mesh::GetNodes() == nullptr` — so constructing a `meq::FieldTransfer` on a mesh
put every later `GridSampler` on it onto the Newton route, 159× at 129² and 83×
at 513², with every answer unchanged. The fast path tests the geometry's
**degree**; `CLAUDE_MFEM.md`, *Traps in MFEM's API*, carries it.

*From `TODO` *Two sampling questions left behind by stage 7*.*

## Out of scope

### Iron and the evolution problem

**Out of scope.**

* **Iron.** Ferromagnetic structures make `K − I` nonlinear inside `Γ`. The
  method accommodates it — that is exactly the coupling-at-a-distance paper's
  hypothesis — but it is a second nonlinearity and belongs after FB-6.
* **The evolution problem.** Quasi-static resistive evolution, CEDRES++'s
  Problems 5 and 6.

*From `FREE-BOUNDARY-PLAN.md` §9, which keeps its heading. That list also named
the inverse problem — coil currents to achieve a shape — which is dropped here
because it is built: IV-0 to IV-6, [M-250](MEASUREMENTS.md#m-250) to
[M-253](MEASUREMENTS.md#m-253). Its fourth item, the original free-boundary
code, is under *The Lackner / von Hagenow Green's-function route*.*

### What the X-point pathway does not reach: the plasma edge, inversion through a separatrix, a fixed boundary on one

**Out of scope.** The X-point pathway (XP-0 to XP-4) makes a diverted
*free-boundary* solve reachable. It does **not** make a diverted plasma
reachable for:

* **`PLASMA-EDGE-PLAN.md`**, which needs transfer paths across the plasma edge
  and whose §7 already records that both path families give out at a corner.
  That plan is limiter-only and this item does not change it.
* **the flux-surface inversion.** `INVERSION-PLAN.md` is complete, IN-5
  included — but IN-5 delivers open surfaces as Chebyshev in arc length, and a
  **disc chart still has no meaning through a separatrix**. The tracer is not
  X-point aware: a level AT one stalls at the saddle, where the level set is not
  a 1-manifold, and `Stalled` is the honest answer.
* **a FIXED-boundary solve on a separatrix**, which is a domain with a genuine
  re-entrant corner and is why `ExtensionConvergence` takes `Γ` to be
  `ψ = −0.03` rather than `ψ = 0`. Nothing here rescues that, and nothing needs
  to: free boundary is what removes the need to mesh the separatrix at all.

*From `FREE-BOUNDARY-PLAN.md` §10.6, *what it does NOT cover*.*

### The inversion's out-of-scope list: topology, 3D, straight-field-line angles, `p` contours and the inverse-equilibrium formulation

**Out of scope.**

| | |
|---|---|
| **Topological classification** — contour trees, open/closed, saddle crossing | `INVERSION-PLAN.md` §6, where it is deferred. It belongs with free boundary, where `χ(annulus) = 0` forces X-points to exist. |
| **3D / stellarator geometry** | MEQ is axisymmetric. DESC's basis is borrowed; its toroidal Fourier factor is not. |
| **Straight-field-line, PEST, Boozer or Hamada angles** | Stability and gyrokinetic codes want these; flux-surface *averages* do not. Build them on top of IN-3 if a consumer appears. Note the practical warning that a straight-field-line angle resolves the low-field-side midplane poorly, which is where ballooning modes live. |
| **Contouring `p` in a rotating equilibrium** | A different curve family (`INVERSION-PLAN.md` §6). Not hard, but not this. |
| **The inverse-equilibrium formulation** | VMEC, DESC and CHEASE make the surfaces unknowns of the *solve*. That is a different code, not a post-processing step, and it hard-codes star-shapedness into the discretisation. (This is the flux-coordinate formulation, not MEQ's inverse solve for coil currents, IV-0 to IV-6, which is built.) |

*From `INVERSION-PLAN.md` §9.*

### Anisotropic pressure, NetCDF profiles and MaNTA coupling

**Deferred.** The rest of the physics is deliberately not yet on the roadmap:
anisotropic pressure, several profiles read out of one NetCDF file, and coupling
to MaNTA. `TODO` carries each with what has been established; anisotropic
pressure has **no reference pinned**, which is the first thing it needs before
any of it is worth building.

*From `ROADMAP.md`, Deliberately not yet.*

## Recorded in place

These lists are woven into live material or are cited by section number, so
they stay where they are; this index makes the catalogue complete.

- `BORDERED-GLOBALISATION-PLAN.md` *STATUS* table and §8's stages: BG-1 `PicardSweeps` falsified, three silent wrong answers (M-116); BG-2 `setFrozenPotential()` not built as written, replaced by `meq::FieldLinearisation`; BG-3 `BorderedPicardThenNewton` falsified, 0 of 6 (M-119); BG-7 continuation in `μ₀ I_p` dead, see *Border scaling: a per-border ruler for the bordered merit*. Interleaved with the DONE stages BG-0, BG-4, BG-5, BG-6.
- `BORDERED-GLOBALISATION-PLAN.md` §3.2: the obvious `buildForms()` route for a frozen source, rejected for four reasons in favour of `SourceIntegrator`, itself not built as written.
- `BORDERED-GLOBALISATION-PLAN.md` §4: family (b), freeze the borders and hand off; §4.3's basin question answered no by BG-1 (M-116).
- `BORDERED-GLOBALISATION-PLAN.md` §5: family (c), continuation in `μ₀ I_p`, dead — every `λ` down to 0.05 fails on machine B, and machine E reaches a degenerate branch (*Border scaling: a per-border ruler for the bordered merit*).
- `BORDERED-GLOBALISATION-PLAN.md` §10: the campaign's kill criteria, several of which fired (M-116, M-119).
- `BORDERED-GLOBALISATION-PLAN.md` §12: removing the line search's full re-assembly, deferred at a 1.05× ceiling (M-138).
- `CLAUDE_DEVICE.md` *Single precision*, single precision's setup loss: `extensionFdStep()` as the suspect (falsified); `meq::pathTolerance()` built, swept and reverted whole, upstream's bisection termination taken instead (M-242).
- `CLAUDE_DEVICE.md` *Which integrators need device offload, and the batched routes*: `AssemblyMode::Batched` a 15 to 27% host loss (closed, M-99); `LocalFactorMode::Batched` not the host default (declined); `AssemblyMode::Batched` under a Device (declined, M-215, M-209); the ordered list's items 2 and 3 (closed), "`HDGExtensionIntegrator` disqualifies the condensation cache" (falsified), nothing for the linear-form integrators (declined) (M-185).
- `CLAUDE_MFEM.md` *`../mfem-hdg-dev` is receive-only*, the `Reconstruct()` row: MEQ's two candidate fixes (falsified, measured); upstream's `own_internal` flag (withdrawn); three explanations of why `cuda` never met it (falsified).
- `CLAUDE_FB.md` *FB-6 is beaten* (the limiter pinned to a dof): `LimiterConstraint::NearestDof` replaced and kept as a control, and the analogy that excused it falsified (M-08).
- `CLAUDE_FB.md` *FB-6 is beaten*, and `FREE-BOUNDARY-PLAN.md` §7.20 *What is left*: *Painting the limiter by centroid, and a precedence rule* — painting falsified (an `O(h)` staircase), cutting elements for a smooth limiter declined, a `SurfaceAttribute`/`R`-`Z` precedence rule declined (both refused).
- `CLAUDE_FB.md` *At `j = 0` the question does not arise* (the tiling paragraphs): *Lost coverage as the cause of the red tiling, and relaxing the gate* — withdrawn and declined (M-12).
- `CLAUDE_FB.md` *At `j = 0` the question does not arise*: `PicardThenNewton` or an exact start for `j = 0` falsified (M-11); a cut rule to make `j = 0` solvable closed; the PLASMA-EDGE route at `j = 0` deferred; the `F·δ(Ψ)` surface term declined.
- `CLAUDE_FB.md` *FB-5* (the support fill): the face-neighbour leak estimate half falsified (2,275 leak), blocking the saddle element declined, fixed rings reverted, ring depth declined. `FREE-BOUNDARY-PLAN.md` §10.3 carries the same list plus "latent not broken" falsified, interior-only declined and freezing the fill not needed.
- `CLAUDE_FB.md` *PE-0*: PE-0's acceptance as written falsified (M-74, M-75), the X-point corner as the cause falsified (M-76), quadrature as the `j = 2, k = 3` anomaly falsified (M-84).
- `CLAUDE_FB.md` *§7.12b's fixture was the defect*: "machine C's guard refusal is a true property" withdrawn, a guard on `axisInsidePlasma` reverted (M-109); a point-in-polygon vessel declined (M-110).
- `CLAUDE_FB.md` *The bordered path has a globalisation*: *PicardSweeps and BorderedPicardThenNewton* — freezing the residual per `BORDERED-GLOBALISATION-PLAN.md` §3.2 declined, `BorderedPicardThenNewton`'s predicted value falsified (M-119, M-253), the handoff tolerance useless, BG-3's re-purposed `PicardOnly` and BG-4's two `solve()` calls declined.
- `CLAUDE_FB.md` *Meshing the half-disc*: linking gmsh's C++ API into `meq_core` declined, re-meshing during a solve not on the plan.
- `CLAUDE_FB.md` *MEQ against freegs4e*: the causes eliminated on machine A — `h`/`p` refinement, freegs4e's grid, the conductor model, mode truncation, the exterior radius, the `ff′` amplitude, the profile conversion, the support, the source arithmetic, `I_p` — and the forward/inverse mode mismatch as the explanation falsified (M-105, M-107).
- `CLAUDE_FB.md` *MEQ against freegs4e*: contour against fitter falsified (M-62), caching only `Jtor` columns not pursued, two 4.5e-08 gap hypotheses falsified, `--seed-from` 10× falsified (1.21×), the `n³` sparse-LU explanation falsified and multigrid/VCYCLE closed, `±eps` averaging falsified and a BEM declined, von Hagenow for FB-6 declined, pyMFEM not used.
- `CLAUDE_FB.md` *MEQ against freegs4e*, "Deliberately not established": what the fixed-boundary comparison does not establish — free boundary, the flux-surface machinery, a real reconstruction (out of scope for that comparison; its free-boundary reason is corrected in place).
- `CLAUDE_HDGGS.md` *Why MEQ's Newton struggles where other Newton solvers do not*: the "easier problem" red flag withdrawn (M-34); continuation in the source amplitude kept out of the driver (declined); `max|∂F/∂ψ|/λ₁` as a predictive detector (declined); "Newton is the problem" and "free boundary is the problem" (falsified); NPC, and any ordering, as a fix for a stiff source (falsified); the linear-Picard control reading (falsified); globalising the outer iteration with `KIN_LINESEARCH` (falsified); Anderson depth ≥ 2 (declined).
- `CLAUDE_HDGGS.md` *The flux is seeded, and it buys a diagnostic and not one iteration*: leaving the flux block at zero, and nodal interpolation of the guess, both wrong; seeding the exterior coefficients, no gain (falsified) (M-29).
- `CLAUDE_HDGGS.md` *The assembled flux is −q, and the papers' τ carries the wrong sign*: `bsymmetrize = false` (falsified); reinstating `τ = 5.0`, MFEM's built-in LDG `τ`, and deriving `τ` from the local coefficient (declined).
- `CLAUDE_INVERSION.md` *`CriticalPointFinder`'s axis and `GradShafranovSolver::psiAxis()` are THE SAME POINT*: "Three design decisions in it are worth not undoing" — `tryFindAxis()` declined for `sweep()`, a two-sided test declined for the one-sided one, and any rule but "largest `Ψ` wins" declined (M-46)
- `CLAUDE_INVERSION.md` *Driving by `q( ψ )`*, "Neither obvious repair reaches it": bounding the GMRES iterations (Falsified: back-substitution goes infinite) and a rank-revealing dense solve (does not rescue the line search; kept for another reason)
- `CLAUDE_INVERSION.md` *IN-3: the fit, and the angle is not free*: the geometric poloidal angle (Falsified, `L^-1.2`; *The geometric poloidal angle*), a parity control (useless, cond 9e16) and "the sample layout is the conditioning lever" (Falsified)
- `CLAUDE_INVERSION.md` *IN-4: the answer is that the question was the wrong one*: the `ψ`-varying element (Closed, none built), an explicit spectral-width penalty (Falsified, loses 44×; kept as the losing column) and an undamped pseudo-inverse (fails at `L = 16, 20`)
- `CLAUDE_INVERSION.md` *Eigen, and what the swap did and did not move*: Eigen as optional, `CompleteOrthogonalDecomposition` and `EIGEN_USE_BLAS`/`LAPACKE` (all Declined)
- `CLAUDE_INVERSION.md` *The cache is keyed bitwise, and a hash was rejected on §8's own wording*, and `INVERSION-PLAN.md` IN-6's cache box: *A hash as the cache key* — a hash (Declined: "unlikely" is not "impossible") and `==` (Declined: `memcmp`, for signed zero and NaN)
- `CLAUDE_INVERSION.md` *What is next*: *Vandermonde + GEMM extraction* (34.7×, not taken: the `(Ψ, θ)` file is written from traced nodes, not a `SurfaceFit`), a guard on the seed hint (Declined), the tolerance ladder in `SurfaceAverageConvergence` (Reverted) and segment-endpoint band marking (Reverted, per-node instead) (M-51, M-71)
- `COIL-SUBTRACTION-PLAN.md` §0c: the "shrink Gamma" argument, falsified both ways — Gamma ≈ 2.1 with a 2.6× saving, and `setExteriorConductors()` doing it (2 of 23 conductors clear).
- `COIL-SUBTRACTION-PLAN.md` §7.1: a clearance rule, perturbing a coincident node, and taking the finite part, all declined in favour of refusing coincidence only.
- `COIL-SUBTRACTION-PLAN.md` §13.5 and §14: subdivided filaments as the default (declined, 3.50× / 10.81× worse against a point-filament reference; M-154, M-166), `InterpolatedFlux` as the default (declined, 7.5e-02; M-171), root polishing and exact `ψ` reporting (did not help), and an mgrid-style interpolated grid (declined).
- `DEVICE-PLAN.md` *What is next* and §0: "porting the last reader opens every gate" (falsified, M-237); phase-C scalars measured and left (closed, M-240); Boost.Math to the device not needed, "interpolatory HDG dodges the device problem" (falsified)
- `DEVICE-PLAN.md` §1.2: `transmissionConstraint()` as a defect (withdrawn); two attributions of the `postProcess` abort, MEQ's aliases and a stack `Vector` (falsified, M-130)
- `DEVICE-PLAN.md` §3.0 table and §3.3 order table: M-199 "half taken" (falsified); the limiter not worth porting (M-201); the dense Schur port moves nothing (M-227); both cuDSS predictions wrong (M-229); the `UseDevice( true )` arm at 17× the bus (reverted, M-229); "the plan was a poor guide" (M-200, M-202); leg share as the yardstick (falsified, M-232)
- `DEVICE-PLAN.md` §3.1: the Amdahl prediction of 1.1–1.3× against a measured 2.14× slower (falsified, M-173)
- `DEVICE-PLAN.md` §3.4 and §3.4a: "none of it is add a kernel" (falsified, M-230, M-231); the limiter not worth porting, and `I_p` and the limiter misclassified as reductions (M-201)
- `DEVICE-PLAN.md` §3.5: "the batched LU is innocent" (M-200 falsified by M-209)
- `DEVICE-PLAN.md` §3.7 and §4: "the first 1.74× is not MEQ's" (M-185, M-188); stage 2 threaded-under-device refused (M-188); stage 4 `x_loc` against `b_tr` (falsified, M-207, M-223); stage 5 no speed on this machine (declined)
- `DEVICE-PLAN.md` §6: the stack-`Vector` attribution and "an inexact integrand raises the Newton count" (falsified; M-130, M-203, M-207, M-223)
- `DEVICE-PLAN.md` §7: "do not convert first" (declined); the §7.2 prediction inverted (M-217); per-pair rate criteria (M-218); porting gslib to single (declined; M-175, M-176)
- `FREE-BOUNDARY-PLAN.md` §4.1: `NonlinearOrdering::LineariseThenCondense`, deleted upstream, and `NORMALISED-LINEARISE-FIRST.md` with it.
- `FREE-BOUNDARY-PLAN.md` §5.3 (the blockquote), §7.10 and §8: *Cut quadrature, Algoim and moment fitting* — ignoring or differencing the cut-rule derivative, cut quadrature not built, MFEM's cut backends quad-only, the PLASMA-EDGE route deferred, a cut rule for `j = 0` not worth it.
- `FREE-BOUNDARY-PLAN.md` §6: §6.1's rectangular integrator not needed; §6.3/§6.4's *Auxiliary globally-coupled unknowns*, closed (M-98).
- `FREE-BOUNDARY-PLAN.md` §7.4a: a filament (log singular), a finite cross-section and a naive C¹ match as FB-1's exact answer — rate capped at 2.5; the match ships as the `p = 0` control.
- `FREE-BOUNDARY-PLAN.md` §7.11: CEDRES++'s ITER case, Serino's ITER case (a proprietary seed), building an ITER coil set, and `γ = 1.395` as the opener, all declined.
- `FREE-BOUNDARY-PLAN.md` §7.12 and §7.12a: *η₆ in the marking* — "boundary elements lose Dörfler" falsified, `η₆` summed into `η` reverted, one weighted sum declined, the growing-domain hypothesis ruled out (M-05).
- `FREE-BOUNDARY-PLAN.md` §7.12b and §7.16: §7.13's verdict stale, the "corner artefact" falsified; "`h` does not cure, `p` does" falsified, the amplitude-fixed formulation failing falsified, the conductor-shrink test establishing nothing.
- `FREE-BOUNDARY-PLAN.md` §7.18 (verdict table): the coil sweep FALSE, the vertical field NOT BORNE OUT, branch crossing FALSE, the corner entries COULD NOT REPRODUCE, coil continuation not needed, the argmax MISLEADING, "vacuum refused" FALSE, the `I_p` ramp NOT BORNE OUT, "prescribe consistent currents" evidence gone.
- `FREE-BOUNDARY-PLAN.md` §8 (Risks): the `F(0, z) = 0` freebie falsified, a cut rule for `j = 0` not worth it, *A backtracking line search on the full residual* making every case worse (M-22).
- `FREE-BOUNDARY-PLAN.md` §7.19: "a = 0 EXACTLY" wrong, the predicted rate half wrong, a cross-mesh projection declined.
- `FREE-BOUNDARY-PLAN.md` §7.20: the NearestDof analogy falsified, regenerating from 513² cosmetic, an additive `soft²` falsified (5× worse), a monotone `ψ_bnd` assertion falsified, "refinement cures branch sensitivity" falsified.
- `FREE-BOUNDARY-PLAN.md` §10.8: initialising `ψ_bnd` from its own border row declined (11 → 200 steps), the `hasInitialGuess()` gate falsified, the axis seed innocent, "make the argmax cleverer" declined.
- `FREE-BOUNDARY-PLAN.md` §11.6: the first `checkAxisSource` reverted, "the fill must not be changed", last-amplitude truncation declined, a 1e-2 threshold declined.
- `FREE-BOUNDARY-PLAN.md` §11.7: a fixed `a = 0.30`, the textbook on-axis loop formula, and a looser assertion past `R ≈ 1.20`, all declined (M-45).
- `FREE-BOUNDARY-PLAN.md` §12, FB-R and FB-T: the gap between the two discrete equilibria does not close (M-132); a prescribed current does not cure confinement (M-134); `setBoundaryFluxInitial` is never a TOML key.
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` IV-0 and §4.6: a flat coil list, and freeing every conductor, falsified (M-244).
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` §4.5: targets as hard equations and a square corner, declined — it fails on machine A.
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` §4.5, §5 and §16 items 7–8: inequality constraints on coil currents and forces, and gap constraints, out of scope for Part One.
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` §6.3: free exterior-conductor currents and the mixed case, declined (M-95).
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` §8: "the inverse solve is a superset of either globalisation plan", declined for three reasons.
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` IV-6: a zero-current cold start (0 of 14), the basin-finder claim carried over from `freegs4e` (one conversion), and a support re-decided inside Newton plus topology-holding continuation, named and not built (M-253).
- `INVERSE-AND-RECONSTRUCTION-PLAN.md` §13.4: modifying the Grad–Shafranov operator for an iron core (falsified) and porting Appel & Lupelli's factor of 2 (declined).
- `INVERSION-PLAN.md` §4.1: the geometric angle (Falsified) and the factorial-sum Zernike (Declined)
- `INVERSION-PLAN.md` §4.3: the flux Taylor step as a control and not a fallback, `PathLiftCoefficient` (cannot), blending toward `ψ = 0` (Falsified) and the "third argument for `ψ`-elements" (Withdrawn) (M-48)
- `INVERSION-PLAN.md` §4.4: the `ψ`-varying element (Closed, answered by IN-4) and two of its constraints struck
- `INVERSION-PLAN.md` §6: the global-structure work (Deferred) and a split tree / merge (not needed)
- `INVERSION-PLAN.md` IN-2: "two expectations died" — `ψ*` buying `k+2` in an average and an average escaping the metric trap (both Falsified)
- `INVERSION-PLAN.md` IN-3: a parity basis (Falsified, useless) and "the layout is the lever" (Falsified)
- `INVERSION-PLAN.md` IN-4: the width penalty (loses), the gauge as a null subspace (wrong), an undamped pseudo-inverse (fails) and the prescribed-angle acceptances (not re-asserted)
- `INVERSION-PLAN.md` §11: the cost model (Falsified, M-49), threaded BLAS for extraction (Declined) and automatic per-surface parallelism (Declined)
- `MIRROR-PLAN.md` item 6 and its fix: tuning `R_ref` (Falsified, M-167), a midplane gauge (Declined) and RoPP's average gauge (Rejected)
- `MIRROR-PLAN.md` *Risks*: risks 1 to 3 retired or dissolved (M-169)
- `PLASMA-EDGE-PLAN.md` header: the whole plan deferred on cost-benefit, `j ≥ 1` already giving `k+2` at `k ≤ j`; "the transfer costs an order" withdrawn (M-83); `k+1.5` as PE-0's acceptance falsified (M-81).
- `PLASMA-EDGE-PLAN.md` §6, notes under the staging table: quadrature and the fixture ruled out as the `j = 2, k = 3` anomaly's cause, the "270° staircase" withdrawn (M-84), and the corner as the cause of fragility falsified.
- `PLASMA-EDGE-PLAN.md` §7: the X-point diagnosis half wrong, and the evidence against the premise disposed of.
- `PLASMA-EDGE-PLAN.md` §9.2 and §9.4: `φ = ψ − ψ_X` as the defining function, falsified; and three causes ruled out — cone widening, a P.1 violation, the mesher / datum / path family.
- `PLASMA-EDGE-PLAN.md` §9.5 and §9.7: an interface at `ψ_X − ε` as a fallback rather than a design, and an `O(h)` disc exclusion, superseded by the closed-form local solution.
- `ROADMAP.md` plan-file table: the `THREADING-PLAN.md` row (E and F declined, M-137, M-138), the `BORDERED-GLOBALISATION-PLAN.md` row (the ladder closes none of six; do not start its §8; M-116, M-117, M-119), the `COIL-SUBTRACTION-PLAN.md` row (caching `psi_c` alone superseded, M-170; subdividing filaments moves MEQ away from `freegs4e`, M-166), the `FREE-BOUNDARY-PLAN.md` row (FB-T does not cure confinement, M-134) and the `PLASMA-EDGE-PLAN.md` row (not scheduled) (M-134, M-137, M-138, M-166, M-170)
- `ROADMAP.md` *So what is next* item 1 (l.80-117): `η₆` in the marking does nothing; coils the wrong instinct, a vertical field makes no core, vacuum start meaningless, continuation not needed; the §7.12 / §7.13 attempts (falsified / declined) (M-05)
- `ROADMAP.md` *So what is next* item 1, *The conductor model is what is left*: no `Γ` excludes the conductors on any of seven machines, and the mixed case cascades (falsified) (M-95, M-96)
- `ROADMAP.md` *Division of labour*: auxiliary globally-coupled unknowns no longer worth asking for (M-98); BEM coupling not wanted, ExteriorDtN instead; the local-solve seed not filed (deferred)
- `ROADMAP.md` item 1: predictive solver selection killed, and continuation must not go in (falsified / declined) (M-24)
- `ROADMAP.md` item 4: upstream's "under 6%" (withdrawn, M-126); host `LocalFactorMode::Batched` (declined); "`ψ_ax` borders need condense-first" (falsified); the `LineariseThenCondense` investigation deleted (falsified / reverted)
- `ROADMAP.md` item 6: `mkl_set_num_threads_local` and `LocalFactorMode::Batched` neither needed; two configurations to avoid; the PARDISO `n ≈ 3000` ceiling and run-time trace choice no longer worth an item (closed)
- `ROADMAP.md` item 8: FB-4 has no cut rule (declined); rectangular integrators not needed (closed); cut quadrature closed — MomentFitting −30/−63, Algoim rejected (declined)
- `ROADMAP.md` item 10: a non-local Jacobian not needed, an outer loop instead (declined); damped Picard (falsified); global-structure work deferred
- `ROADMAP.md` item 12.1: interpolating `F(ψ_h)`, the CSZ 2019 route, loses `k+2` (declined)
- `THREADING-PLAN.md` *STATUS* and §0: two of six candidates were not threading problems (the border dense solve is sparsity; the remainder is the build, M-77); the `I_p` half (0.5%) and the transmission sweep (~1%) not worth doing; payoff against install-before-npcthread (declined) (M-77, M-138)
- `THREADING-PLAN.md` *`border dense solve` is a sparsity problem* and *Cross-cutting*: BLAS-3 GEMM for the border solve (declined; reassociates); design rejections inside A–C — `reduction(+:)`, gating on source type, threading the dedup (declined)
- `TODO` *Read p'(psi) and FF'(psi) from NetCDF*: splitting the path at the last slash, and GDAL's `NETCDF:"f":var` / OPeNDAP's `?` separator, not adopted (declined); options (b) differentiate the spline again and (c) a values-only fit (declined)
- `TODO` *Report q( Psi )*: the prescribed-field `q` refusal stands (declined); the vacuum-field key (deferred)
- `TODO` *More nonlinear test cases*: the "stiff" premise wrong for three of four (falsified); the Weierstrass dependency undecided (deferred)
- `TODO` *Performance: profile the inside of a solve*: neither named route taken (closed); `LocalFactorMode::Batched` "is not the answer" (declined); "the culprit is the dense LU" (falsified); two configurations to avoid; perf / Caliper (closed); the InvertA/D corrections (withdrawn); `LSsolveType` (closed); and *Two arguments that decide the rest* — partial assembly needs quads, local device timings must not decide, AmgX (SPD; asymmetry 0.54), SuperLU_DIST / MUMPS / STRUMPACK (MPI) (declined)
- `TODO` *Half a dozen fixed-boundary cases* (deprioritised, M-148) and *Give up on a bordered Newton that has stopped moving* (deferred, M-236; `setLineSearchMerit( Field )` and `BorderedPicardThenNewton` do not rescue single, falsified)
