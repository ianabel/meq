# Where meq is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27 and again 2026-08-29.
`CLAUDE.md` is the operational record and is authoritative on anything technical;
`DRIVER-PLAN.md` is the stage-7 design; `TODO` holds work that is understood but
not scheduled. **This file is only about order** — what to do first, what waits
on what, and what is deliberately not being done yet.

## Nothing is red

**`HighBetaConvergence` was the one failing test and it is green.** Profiles
specified in NORMALISED flux — `Ψ = (ψ − ψ_bnd)/(ψ_ax − ψ_bnd)`, which is how
`refs/GourdainContour.pdf` §V eq (39) poses them and how every equilibrium code
does — need `ψ_ax` inside the residual, where the Jacobian can see the non-local
terms it contributes. It is now there: `meq::NormalisedSource` is the interface,
and `GradShafranovSolver::setSource( NormalisedSource &, double )` carries `ψ_ax`
as an unknown and closes the pair by a **bordered Newton**, one factorisation and
two backsolves per step.

Each of the cheaper options was measured and killed on the way, and each
measurement is worth not repeating. With `ψ_ax` **fixed** the profile is inert
(amplitudes 1 and 512 give an identical solution, because the solve only ever
reaches `Ψ = 0.0013`). With `ψ_ax` iterated **outside** the solver the fixed point
is degenerate, and mapping the outer iteration afterwards says why: it has a
**pole beside its own fixed point** — at `ν = 2`, `A = 1` the fixed point is
0.3059 and the pole is at 0.2996. And at a *fixed* `ψ_ax` the equation has a
small solution as well as the equilibrium, so the run needs a starting guess of
about the right height; that is part of the problem statement, not an
optimisation.

The test that can *see* the missing terms is `Normalisation::Decoupled`: the same
solver with `∂R/∂ψ_ax`, `∂(max ψ_h)/∂λ` and `∂(max ψ_h)/∂ψ_ax` zeroed and
everything else identical. Coupled reaches 4.4e-15 in four iterations; decoupled
is still at 8.24e-2 after fifteen, having started at 8.31e-2. Full detail in
`CLAUDE.md` under *Newton, and the obligation it creates*.

**What this unblocks is free boundary**, where `ψ_ax` and `ψ_bnd` are both
unknowns: `ψ_bnd` is a second border row and column of exactly the same shape,
and `meq::NormalisedSource` is where it goes.

**What is left over from it**, and neither is on the critical path:

* `meq::NormalisedMHDSource` — the production source built on two `meq::Profile`s
  in normalised flux — is written and unit tested but is **not reachable from a
  TOML file**. Wiring it through `Config` and `SourceFactory` is driver work and
  belongs with item 1.
* The bordered path is `Globalisation::None` and
  `NonlinearOrdering::CondenseThenLinearise` only, and refuses the others loudly
  rather than quietly doing something else. Both refusals have reasons recorded
  beside them; neither has been needed.

**And it turned up an MFEM finding that is not filed anywhere**, which is item 5
below: `DarcyHybridization` captures the element-local Newton's initial guess at
`FormLinearSystem()` time and never refreshes it, so every local solve in every
later residual evaluation restarts from it however far the trace has travelled.
That is a cost on an ordinary Newton path and a *correctness* problem for
anything that differentiates the residual by differencing it. meq works around
it; nobody has been told.

## So what is next

With nothing red, the order is what it was minus the blocker:

1. **Finish the driver** — item 1. Two things left in it: wire the normalised
   source through `Config`, and the GSLIB warm start.
2. **Then free boundary — and it now has a plan**, `FREE-BOUNDARY-PLAN.md`.
   HDG on a polygonal subdomain coupled at a distance to an exterior operator on
   a semicircular artificial boundary, by
   `refs/CouplingAtADistance.pdf` — whose own reference [5] is meq's stage 5, so
   the coupling is the extension technique already in the tree with the datum
   unknown instead of zero. Staged FB-0 to FB-6 with an acceptance measurement
   each. See item 8.
3. **Hygiene** — item 3 — alongside either. `README.md` is the oldest debt in
   the tree and still describes a project that did not exist.

Items 4 and 5 are the other tree's to act on; item 6 is meq's and is small and
not urgent; item 7 is a closed investigation kept for its answer; item 8 is free
boundary, which is step 2 above.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done, and **stage 7 is finished: meq is a program that solves on a
curved boundary, refines its own mesh, and restarts from a previous answer in one
Newton step.** The suite is **22 of 22**, and it is green rather than
green-with-a-known-red: the last standing failure was `HighBetaConvergence` and
the work it was asserting is done. The only thing the driver still refuses is
`[boundary] Type = "exact"`.

## The nonlinear question is settled

It ran through four wrong answers, which is why it is written down rather than
summarised away. meq's Newton was thought to fail on the stiff GS-2 sources.

* **Not Newton, and not the boundary condition.** `refs/MFEM-GS-Newton.pdf`
  reports Newton robust where "conventional Picard-based solvers fail to
  converge", and calls its own *fixed*-boundary solver "significantly easier"
  than the free one.
* **Not the hybridization ordering.** This one was believed and acted on.
  `NLOrdering::LineariseThenCondense` removes every element-local nonlinear
  solve (`GetNumLocalNLIterations()` reads 0) and the stiff cases get *worse*.
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

Numbers in `CLAUDE.md`, *Why meq's Newton struggles* and *Picard, then Newton*.
**Nothing on this path blocks anything else.** What follows is ordered on value.

## Division of labour

**The MFEM tree is worked by another agent** and is receive-only: fetch from it,
write requests into its `doc/`, never branch, commit, check out or build there.
meq owns `src/`, `tests/` and `apps/`.

**Every request meq has filed except one is landed in the code meq builds
against**, which is the only test of "landed" meq can apply:

| filed | state |
|---|---|
| `HDG-LINEARISE-THEN-CONDENSE.md` | landed — `setNonlinearOrdering()` uses it |
| `DIRECT-SOLVER-SYMBOLIC-REUSE.md` | landed — `SetReuseSymbolic()` is on |
| `HDG-DEFECTS-FROM-MEQ.md` | retired by that tree. Three of its four are verifiably closed — one fixed, one fixed, one withdrawn as not a defect; the fourth, `ComputeHDGFaceEnergy()` ignoring an installed stabilisation, meq has not re-measured and does not use. `CLAUDE.md` has the breakdown |
| `HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md` | **landed and retired** — the fix is *"The postprocessing closes on the element average, always"*, and the document is gone from that tree's `doc/` |
| `HDG-ELEMENT-LOCAL-PARALLELISM.md` | still there, and meq has seen no change |
| `HDG-BEM-COUPLING-FROM-MEQ.md` | **filed 2026-08-29**, for free boundary — see item 8 |

The reconstruction one is the interesting entry, because it changed what the
driver does: `ψ*` is a post-processing on **every** element now, so the adaptive
loop is back on the published estimator rather than one order down. See item 3.

**What is NOT filed is the local-solve seed** — item 5. It was found from meq's
side this week and nothing has been written into that tree about it, beyond one
paragraph inside the coupling request using it as the argument for why a caller
should not have to difference a condensed residual.

**And the coupling request is a request for a capability, not a finding.**
Nothing in that tree has been measured to be wrong by it; §2 and §3 of it are
things that do not exist rather than things that misbehave, and it says so.

---

## 0. FOUR MFEM branches, one integration branch — standing maintenance

**It is four now, not three, and this is the entry to read before rebuilding
`../mfem/install`.** meq needs work from

| branch | what meq needs from it |
|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — the curved `Γ`, and `TransferredDatumCoefficient` |
| `direct-solver-symbolic-reuse` | `SetReuseSymbolic()` on the direct solvers |
| `gf-hdg-linearise-first` | `NLOrdering::LineariseThenCondense` |
| **`gf-hdg-dev`** | **the reconstruction fix — *"The postprocessing closes on the element average, always"*** |

**`meq-integration`** is their merge in `../mfem/mfem-src` and is what
`../mfem/install` is built from. Local only, never pushed, re-created whenever any
of the four moves — so nothing may be committed directly to it.

**Verified 2026-08-29, and it is a foot-gun:** none of the four contains any
other, and the reconstruction fix (`25b7dae612`) is reachable from `gf-hdg-dev`
and from **no** other branch. Re-creating `meq-integration` from the three this
file used to name would silently drop it, `ψ*` would go back to being a different
function wherever `∂F/∂ψ` vanishes, and the adaptive loop would go back a full
order — with no error and a green suite until
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` caught it.
`gf-hdg-dev` cannot simply replace `gf-hdg-subdomains-dev` either: it still has
no `fem/darcy/extension_hdg.*` at all, so there is no shorter merge.

`CLAUDE.md` carries the recipe, now four-branch, including the one recurring
conflict on the first merge. **The fourth merge is not characterised** — the
`meq-integration` in the tree today was built with `gf-hdg-dev` merged last and
it works, but whether that merge conflicts, and how, is not written down.

A standing cost of the arrangement, and it falls on meq's side.

## 1. Finish the driver — meq

`DRIVER-PLAN.md` §3–5: **7c output done, 7e the driver done.** `meq config.toml`
parses, solves, prints a residual history and writes `.mesh`, `_psi.gf`,
`_grad_psi.gf` and the `(R, Z)` NetCDF file, with exit codes 0/1/2/3.

What is left, in the order it is worth doing:

* ~~**Wire the extension path.**~~ **Done.** `[boundary.shape]` builds a
  `BoundaryShape`, marks `D_h` out of the background mesh, finds `Γ_h` and
  hands a `VertexConePath` to `setExtension()`. `examples/miller-curved.toml`
  runs it: 450 of 2560 background elements inside, 0 paths widened.
  `theDriverSolvesOnACurvedBoundary` pins it against the library at 1.6e-16 —
  and separately against imposing *zero* on `Γ_h`, at 2.7e-1, because an
  extension that was never applied would converge quietly to the wrong answer.
* ~~**Wire the adaptive loop.**~~ **Done.** `[adaptivity] Enabled = true` runs
  solve → post-process → estimate → mark → refine, stopping at `TargetError` or
  `MaxIterations` and saying which. On the curved path it uses
  `meq::AdaptiveDomain`, so the companion mesh keeps `dist(Γ_h, Γ)/h_loc` from
  doubling every cycle, and it calls `setTransferredBoundary()` automatically and
  logs that it did. `theDriverRunsTheAdaptiveLoop` pins it against the library at
  1.653e-16 and asserts separately that it refined, refined *adaptively*, brought
  η down, and kept assumption P.1.

  **`examples/miller-adaptive.toml`, re-measured 2026-08-29** now that the loop
  is on `ψ*` rather than on the raw potential:

  | cycle | elements | trace | marked | wide | η |
  |---|---|---|---|---|---|
  | 0 | 97 | 480 | 37 | 0 | 4.7348e-04 |
  | 1 | 254 | 1194 | 18 | 0 | 2.4729e-04 |
  | 2 | 342 | 1593 | 18 | 0 | 1.0687e-04 |
  | 3 | 449 | 2097 | 0 | 0 | 6.8664e-05 |

  The numbers this file used to carry — η 2.75e-3 → 4.77e-4 over 97 → 1069
  elements — were the **degraded** estimator, built on `ψ_h` because MFEM's
  reconstruction could not be trusted. Same four cycles now reach 6.87e-5 on 449
  elements, which is what the extra order buys.

  **It cost a prerequisite that was filed as an independent item, and that is
  the transferable part.** The estimator needs `ψ*` — four of eq. (20)'s five
  terms — and `postProcess()` refused the semi-linear path, which is the *only*
  path the driver has. So item 3's `postProcess()` bullet was not optional and
  not separate. See below for what measuring it turned up.

* **The driver's non-linear ladder is in**: Newton, and on *observed* failure
  `PicardThenNewton`, rebuilding the solver because a caught `ErrorException`
  leaves one unusable. Reactive, never predictive.
* **Wire `meq::NormalisedMHDSource` through `Config` and `SourceFactory`.** The
  solver takes `ψ_ax` as an unknown and the production source that needs it
  exists and is unit tested, but a TOML file cannot ask for either. That is the
  gap between "meq can compute a profile equilibrium" and "meq can be *told* to",
  and it is small: a `[profiles] Normalised = true` and the `psi_ax` guess, which
  is a required input on that path and not an optional one.
* **7d, the interpolating warm start** — needs GSLIB, already enabled in
  `../mfem/install`. The exact restart (same mesh, same degree) is written. The
  acceptance measurement is the one `DRIVER-PLAN.md` §4 names: at `k = 3` the
  GSLIB route should start from a residual smaller than the NetCDF route's by
  about `h^{k+1}` against `h^2` — which is what justifies the dependency rather
  than asserting it.

**The nonlinear path the driver ships:** Newton, falling back to
`PicardThenNewton` on **observed** failure. A *reactive* ladder — never a
predictive one. Nothing may be inferred from `F` about which solver to run,
because nothing can be: the ratio `max|∂F/∂ψ|/λ₁` is computable black-box but the
pedestal converges at 7 where the hole fails at 26, which is two points and not a
threshold. **Continuation must not go in**, for the stronger reason that it has
no black-box form at all — see item 7.

**Prerequisite `MFEM_USE_EXCEPTIONS`: done.** `../mfem/install` is rebuilt with
it, and §4.4 — which used to take the process down with SIGABRT — now returns a
reported failure and an iteration count. Exit code 2 is implementable. No
meq-side change was needed: `ErrorException` derives from `std::exception`.

## 2. ~~Symbolic factorisation reuse~~ — **done**

`SetReuseSymbolic()` is on for the Newton path and for the Picard path, whose
solver was hoisted to a member so it had something to reuse across its 122 to 290
factorisations. Off deliberately on the linear path, which factorises once.
`theSymbolicAnalysisIsReusedAcrossNewtonSteps` asserts one analysis against one
factorisation per iteration — a count, not a timing, and the only thing that
could notice the reuse lapsing, since a lapse costs speed and nothing else.

## 3. Hygiene — meq, alongside the driver

* ~~**Rewrite the pedestal tripwire.**~~ **Done.**
  `pedestalConvergenceIsAResolutionThreshold` asserts the two cures — 9
  iterations at `k = 1, n = 32`, 11 at `k = 2, n = 16`, against 42 at
  `k = 1, n = 16` — and asserts the threshold as a *ratio*, which holds whichever
  side the threaded-MKL rounding falls on. The knife-edge point is recorded and
  deliberately not asserted on.
* ~~**Re-run `everyNonlinearPathReachesTheSameExactSolution`**~~ **Done.** Four
  paths now, all reaching L2 2.983495e-05 to every digit printed;
  `PicardThenNewton` gets there in 2 Newton steps after Picard's walk. It also
  turned up something worth keeping: that run emits 2164 element-local
  non-convergence warnings where plain Newton on the same problem emits none, and
  converges anyway.
* ~~**Retire `postProcess()`'s refusal — but measure first.**~~ **Done, and
  measuring first is what stopped it being wrong.** The rate on Example 5 is
  `k+2` — 3.05, 4.05, 5.03 — so the lift works. But the refusal was *narrower*
  than a blanket one should have been: with `∂F/∂ψ` identically zero the frozen
  Jacobian vanishes, the local problem is singular again, and `ψ*` comes back at
  8.98, 5.07, 6.90e11, 4.77 over four meshes against 3.8e-6 … 9.3e-10 for the
  same problem solved linearly. **Finite at three meshes in four**, so a code
  read would have passed it and so would a spot check.

  `postProcess()` now measures its own output — `‖ψ*‖/‖ψ_h‖` in `[0.5, 2.0]`,
  0.9945–1.000000 measured on every working case against 26 upwards on every
  broken one — and throws naming the cause. Filed for the MFEM tree as
  `doc/HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md`.

  **And then the mechanism turned out to be one flag test, which changed the
  answer twice over.** The local post-processing problem is a *pure Neumann*
  problem, singular by construction, and MFEM regularises it with a mean-value
  constraint — correctly. That branch is skipped whenever `nl_src` is set, and
  `nl_src` is set on the mere *presence* of a non-linear integrator rather than
  on whether it contributed a non-singular block. Where `∂F/∂ψ` vanishes it
  contributes nothing, the regularisation is skipped anyway, and a singular
  matrix is factored. **A 1e-12 floor on `∂F/∂ψ` repairs it completely** —
  7.565317 to 0.998525 — which is a perturbation incapable of moving a solution
  and capable only of making a matrix invertible.

  **`nl_src` is per element, so this is not the special case it looked like.** At
  an eighth of the domain dead — a tabulated profile with a flat segment —
  individual elements are 20× wrong while the whole-domain norm is 1.87. Two
  earlier meq-side responses were therefore both wrong and are reverted: a
  `Source::dependsOnPsi()` flag routing linear sources to the linear path, which
  dodges one instance and misses this one; and a runtime check inside
  `postProcess()`, which is a standing defence against a dependency and which a
  global norm cannot make reliably anyway.

  **And it has since been fixed on the MFEM side, so this bullet closes.** The
  request went in with line numbers as
  `HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md`; the fix is *"The postprocessing
  closes on the element average, always"* — unconditional, because the local
  post-processing problem is a pure Neumann one by construction and there was
  never anything to decide. Measured from meq's side, the case that read 20.3,
  64.1 and 61.6 on dead elements now reads **1.0069**, against 1.0048 for
  elements where `∂F/∂ψ` does not vanish at all.

  Three consequences, all of them already in the tree:

  * The driver's adaptive loop is **back on `ψ*`** rather than on
    `Potential::Raw`. `Potential::Raw` stays as a documented option and is no
    longer what the driver picks.
  * The test flipped from asserting the defect to asserting the behaviour:
    `theReconstructionIsWrongWhereTheJacobianVanishes` is gone and
    `thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` is what runs,
    as a regression rather than a tripwire. That is the testing stance working
    exactly as `CLAUDE.md` says it should — the red test was the record, and it
    went green on the day the defect went away.
  * The measurement that the defect never reached the solve stands and is worth
    keeping: `ψ_h` and `q_h` agreed between the two paths to 1.6e-13 and 2.6e-13
    over `k = 1…4` and three meshes, and post-processing left both
    bit-identical. The forward local problem takes its potential block from the
    HDG stabilisation on the interior faces, which never degenerates.
* **Rebuild `η₅` on `TransferredDatumCoefficient`** instead of excluding those
  faces. `mfem::TransferredDatumCoefficient` exists now
  (`fem/darcy/extension_hdg.hpp`), so the thing that blocked this is gone and the
  work is meq's: on `Γ_h` the term currently compares `ψ*` against a trace pinned
  to zero rather than the `φ_h` actually imposed, and `setTransferredBoundary()`
  excludes those faces — **an omission, not a repair**, and the driver says so in
  its own comment. It wants its own convergence measurement rather than a switch,
  which is the only reason it is not a five-line change.
* **Rewrite `README.md`**, which still describes a project that did not exist.
  The oldest debt in the tree, and the one a new reader hits first.

## 4. Element-local parallelism — MFEM, outstanding

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`, still present in that
tree's `doc/`. Probably the largest performance win available: twelve sequential
element loops in `darcyhybridization.cpp` over work independent by construction.
meq has no visibility into whether it is being worked on and should not guess;
what meq can say is that nothing it links has changed.

**`LineariseThenCondense` is what makes this tractable, and it now exists.** That
ordering was pursued as a convergence fix and failed as one — but it does exactly
what it says, running **zero** element-local nonlinear iterations, so under it the
batched work is fixed-size *linear* solves rather than nonlinear ones of
unpredictable iteration count. That is what `linalg/batched/` provides and what a
device wants. So the ordering is not a dead end; it is the enabler for this item,
and that is the reason to keep `setNonlinearOrdering()`.

**The borders are no longer an obstacle to it.** meq's `ψ_ax` bordered Newton was
believed to require condense-first and does not: measured 2026-08-30, it reaches
the same `ψ_ax` under both orderings to every digit printed, and the guard has
been lifted. So of the three things pulling towards this ordering — the batching,
the borders, and every reference in `refs/` — one is now settled in its favour.

**What still blocks it is diagnosed, and it is not what it first looked like.**
`LineariseThenCondense` fails at 60 on §4.2 at every resolution tried, including
ones condense-first solves in five iterations — and the cause is that **its
reduced residual is not a function of the trace**. The retained local fields are
hidden state: evaluated at one trace, `R` differs by **149%** at the published
pedestal width depending on which trace the linearisation was last taken at,
against exactly zero for condense-first. Newton is handed a function that moves
under it by more than its own value, which is why the failure survives
refinement, a line search and Picard alike. It is *not* the Jacobian — at a fixed
linearisation the gradient is the derivative to 5e−9.

**And the same measurement found something about the ordering meq ships.** Under
condense-first the residual is a perfect function of the trace, but the gradient
disagrees with a central difference of it by **100×** once `σ² ≤ 0.01`, because
the element-local solves are at their iteration cap. That is the measured
explanation for the pedestal's 39 wandering iterations. `CLAUDE.md` has both
tables.

**So the way out is a specific ask rather than an investigation**, and the two
orderings need different things: linearise-first needs the field increments to
belong to the outer solver rather than to the operator — which is what NPC's
simultaneous Newton is, and would let one line search damp the whole step —
and condense-first needs its local solves to converge. Both are MFEM's. Until
then the default stays condense-first, and `setNonlinearOrdering()` is a
per-problem choice rather than a project-wide one.

## 5. The element-local solve's frozen seed — MFEM, **found and not filed**

`DarcyHybridization::EliminateVDofsInRHS` copies the flux and potential blocks
into `darcy_u` and `darcy_p` at `FormLinearSystem()` time, and every
element-local Newton afterwards — in every residual evaluation, every
`GetGradient()` and every `ComputeSolution()` — starts from that same vector.
`ComputeSolution()` does not take the output vector as a guess, so there is no
way to supply a better one from outside.

**Two costs, and the second is the one that matters.** On an ordinary Newton path
it is iterations: measured on a high-β profile at `n = 16, k = 2`, 40,000 to
60,000 element-local iterations per outer step, most of them hitting the cap of
100. On anything that obtains a derivative by **differencing the residual** it is
correctness, because a local solve that ran out of iterations returns whatever it
had reached, which is not a function of anything: differencing `ψ_ax` by 9e-6
moved `max ψ_h` from 0.896 to 2.04 and then to 3.84.

meq works around it in `GradShafranovSolver::formSystem()`, re-forming the
reduced system from the recovered state once per Newton step, which puts every
local solve within an iteration or two of its answer. **That workaround is
cheap and meq is not blocked**, which is why this is item 5 and not item 0.

What a request would ask for is small and specific: a way to refresh the local
initial guess without re-forming the system — the obvious shape is
`ComputeSolution()` reading the vector it is handed. **It has not been written**,
and by this project's own rule it should be written only once somebody is sure
the behaviour is settled rather than mid-change in that tree.

## 6. PARDISO runs, and removing the link-line straddle cost 140x — meq

**RE-MEASURED 2026-08-30, AND THE HEADLINE CHANGED.** The straddle described
below was fixed on the PARDISO side — MFEM's `MKL_PARDISO_LIBRARIES` pointed at
Debian's `mkl_sequential` and now points at oneAPI's `mkl_gnu_thread` — and that
correct change **exposed a 140x regression in the test suite**. The stray
sequential layer had been making all of MKL sequential, which is why every
timing this project ever recorded was fast.

What threading MKL actually costs, measured with the columns separated:

* **Not UMFPACK.** `UMFPackSolver::SetOperator` degrades about 40% across the
  whole thread range, never more.
* **`ComputeH()`'s element-local dense LU**, through LAPACK, on blocks of order
  10-30. `k = 2` untouched; **`k = 3` forty-fold worse** at two threads.
* So `MKL_NUM_THREADS=1` is now set on every ctest, and the suite runs in 798 s
  against 2269 s.

**And that pins PARDISO's own scaling out of reach.** PARDISO beats UMFPACK
1.50x on setup sequentially and 2.83x at 8-16 threads — but `MKL_NUM_THREADS` is
process-wide, so buying that means paying 40x on assembly. **Making the trace
solver selectable is therefore no longer the whole job**: it needs either
`mkl_set_num_threads_local()` around the trace solve, or the element-local
factorisation off threaded MKL (`LocalFactorMode::Batched`). Neither is done.

`CLAUDE.md`'s *Threading, measured* has the tables. The rest of this item is the
2026-08-29 state.

## 6a. The 2026-08-29 measurement, superseded

**The "does not run above `n ≈ 3000`" this file used to record was a verdict on
Debian's `intel-mkl` 2020.4.304 and on nothing else.** Against oneAPI MKL 2026.1
with `MFEM_USE_MKL_PARDISO=YES`, re-measured 2026-08-29, it runs at every size
tried and agrees with UMFPACK throughout:

| `k` | trace dofs | UMF setup / solve | PARDISO setup / solve | agreement |
|---|---|---|---|---|
| 2 | 9,408 | 0.0471 / 0.0048 | 0.0141 / 0.0020 | 3.7e-15 |
| 3 | 28,032 | 0.2078 / 0.0254 | 0.0489 / 0.0088 | 1.2e-14 |

`TraceSolverComparison` asserts only the **agreement**; the timings are printed,
because PARDISO is threaded and a threaded timing on this machine is not a
measurement about the code. Two things follow, neither urgent:

* ~~**meq's solver still uses UMFPACK.**~~ **DONE 2026-08-31.**
  `setTraceSolver()` chooses `UMFPack`, `Pardiso` or `CuDSS`;
  `traceSolverAvailable()` reports what the build has; an unavailable choice
  throws rather than substituting. Both paths keep shipping and **UMFPack stays
  the default**, for exactly the reason this bullet gave — oneMKL's licence is
  not everybody's to accept. What changed is that a caller who does have it gets
  1.24x for one line.
* **`MFEM_EXT_LIBS` links oneAPI 2026.1 and Debian's 2020.4 MKL at once**, the
  second pulled in behind SuiteSparse, bringing two extra threading layers. It
  works today only because oneAPI comes first in the link order and its rpath is
  baked in — an ordering accident, not a configuration. Repointing SuiteSparse's
  BLAS at oneAPI would settle it. `CLAUDE.md` carries the full line.

## 7. §4.4, the current hole — a known answer that must stay out of the driver

Solvable: adaptive continuation in the added term's amplitude walks `c₃` from 0
to −18 at `k = 2, n = 32` in **9 solves with 2 retreats**, final solve 8 Newton
iterations, no limit point on the branch. Ten uniform steps stall at −10.8, so
the step control is the whole of it.

**It ramps a constructor argument of the _test fixture_.** `Source` exposes
`f( r, z, ψ )` and `dFdPsi( r, z, ψ )` and nothing else, so a black-box `F` has
no amplitude to continue in and this route is unavailable to a driver on
principle. The expressible analogue `F_λ = λF` is a **different** homotopy —
degenerate at `λ = 0`, where this path starts from the converged pedestal — and
is untested.

So this stays a **test-only investigation**, and its value is that §4.4 is now
known to have a solution meq can reach. Making it general needs a
self-parameterising `Source`: an interface change to argue on its own merits, not
to smuggle in.

---

## 8. Free boundary — meq, planned and not started

`FREE-BOUNDARY-PLAN.md` is the design, and it is a plan in the sense
`DRIVER-PLAN.md` was: nothing is built. The shape of it:

* **The exterior is exact, not a BEM.** With `Γ` a semicircle centred on the
  axis, the exterior Dirichlet-to-Neumann map for `Δ*` is **diagonal** in the
  Gegenbauer basis `C_n^{−1/2}(cosθ)`, with symbol `(1−n)/ρ_Γ` and mass
  `2/(n(n−1)(2n−1))` in the `dΓ/r` weight the weak form already carries. So
  CEDRES++'s double surface integral over an elliptic-integral kernel collapses
  to `N ≈ 20–40` numbers. Derived and measured — against an independent exact
  current-loop field, to 2.3e−9 — in that plan's §3, and **not taken from any of
  the papers**, which is why its §3.4 names the test that would falsify it.
* **The coupling is stage 5 with the datum unknown.**
  `refs/CouplingAtADistance.pdf`'s own reference [5] is Cockburn & Solano, which
  is what `setExtension()` already drives.
* **The augmented Newton is the bordered one, `N + 2` times over.** `ψ_ax` is one
  border column today; free boundary adds `ψ_bnd` and the `N` coefficients, and
  the borders stop being differenced and start being assembled.
* **The order of work is FB-0 to FB-6**, and FB-1 — vacuum with coils, whose
  answer is a sum of loop fields known to machine precision — is the stage that
  tests everything structural against an exact answer. The received wisdom that
  free boundary has no analytic solution is true of FB-4 upwards and false below
  it, and the exact answers available early are what should be spent first.

**The two things most likely to hurt**, both named in the plan rather than left
to be discovered: the axis, where the flux mass `(r q, v)` degenerates and
`BoundaryShape` currently refuses to go — which is measurable **today**, with no
free boundary at all, and is worth measuring before FB-1 — and cut-element
quadrature for the plasma support, which is the one place a published code says
it hit a wall.

What this needs from the other tree is `HDG-BEM-COUPLING-FROM-MEQ.md`. Its §2 —
two rectangular integrators beside the extension ones — is blocking. Its §3,
auxiliary globally-coupled unknowns, is an **improvement on a route that already
works**, and that is a correction: this file claimed on 2026-08-29 that a
differenced border returns exactly zero under `LineariseThenCondense` and that §3
was therefore what made the ordering a free choice. Measured 2026-08-30, it does
not — the ordering re-evaluates the source on every residual, meq's `ψ_ax` border
works under it to every digit, and the guard refusing it has been lifted. What §3
buys is an assembled border instead of a differenced one.
`NORMALISED-LINEARISE-FIRST.md` is the design for meq's half. Plan §4.1.

---

## Deliberately not yet

* **GPU and cuDSS.** Correctness-testable here; this card cannot say whether it
  is worth it. Consumer FP64 is 1/32 of FP32 where datacentre parts are 1/2, so a
  local timing can invert the production conclusion.
* **Partial assembly.** Not implemented for the HDG integrators at all, and its
  case is sum factorisation, which wants tensor-product elements where meq is
  triangles. A discretisation decision, not a flag.
* **The physics** — sonic rotation, anisotropic pressure, NetCDF profiles, MaNTA
  coupling. The actual science. This used to read "all downstream of a driver";
  the driver is now done bar two bullets, so what they are downstream of is free
  boundary and each other. `TODO` carries each with what has been established.
  **`TODO`'s "Try PARDISO again against a real oneMKL build" is answered** — see
  item 6 — and that entry can go.

## The standing rule

Every stage ends at a **measured** number, not at "it runs". That is what caught
the Solov'ev coefficients, the `τ` sign, the `ψ*` regression, the hybridization
ordering, the under-resolution and the inert normalisation — six times the code
or the story looked right and was not. Four corollaries this project has paid
for:

* **A property measured on the easy configuration is not a property of the
  code.** Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on
  the geometry meq exists for.
* **A difficulty measured at one resolution is not a property of the problem.**
  Three sources were called stiff for months on the strength of `k = 1,
  h = 0.05`, and a whole MFEM work item was requested to fix them. They were
  under-resolved.
* **A measurement taken on this machine may not be a measurement about the
  code.** Threaded MKL decides a marginal Newton; a consumer GPU decides nothing
  about a production one. And a *stale* one is worse than none: this file said
  PARDISO stopped working above `n ≈ 3000` for months after the MKL it was a
  verdict on had been replaced.
* **A derivative obtained by differencing is only as good as the thing being
  differenced is a function.** The bordered Newton's border read 1.6e5 where it
  should read 1, and looked exactly like a singular Jacobian, because an
  element-local solve underneath it was hitting an iteration cap and returning
  whatever it had reached. Nothing in the arithmetic was wrong. Before believing
  a finite difference of a hybridized residual, check that the local solves
  under it converged.
