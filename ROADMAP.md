# Where meq is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27. `CLAUDE.md` is the
operational record and is authoritative on anything technical; `DRIVER-PLAN.md`
is the stage-7 design; `TODO` holds work that is understood but not scheduled.
**This file is only about order** — what to do first, what waits on what, and
what is deliberately not being done yet.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done, and **stage 7 is all but finished: meq is a program that solves
on a curved boundary and refines its own mesh.** `meq config.toml` parses,
solves, adapts, and writes. The suite is **18/18** — the pedestal tripwire that
used to fail has been rewritten to assert what it is actually measuring. What is
left of stage 7 is one item: the warm start's interpolating route (§4, needs
GSLIB). The only thing the driver still refuses is `[boundary] Type = "exact"`.

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
Of the plans meq filed, `HDG-LINEARISE-THEN-CONDENSE.md` and
`DIRECT-SOLVER-SYMBOLIC-REUSE.md` have **landed**, and
`HDG-DEFECTS-FROM-MEQ.md` was retired by that tree once all four defects were
fixed and covered. `HDG-ELEMENT-LOCAL-PARALLELISM.md` remains outstanding, and
`HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md` is newly filed — a **narrowing**
of the retired document's §1 rather than a re-opening: the fix landed and works
wherever `∂F/∂ψ ≠ 0`, and returns a different function, silently, where it does
not. meq owns `src/`, `tests/` and `apps/`.

---

## 0. Three MFEM branches, one integration branch — standing maintenance

meq needs work from `gf-hdg-subdomains-dev`, `direct-solver-symbolic-reuse` and
`gf-hdg-linearise-first`, none of which contains the others. **`meq-integration`**
is their merge in `../mfem/mfem-src` and is what `../mfem/install` is built from.
Local only, never pushed, re-created whenever any of the three moves — so nothing
may be committed directly to it. `CLAUDE.md` carries the recipe and the one
recurring conflict.

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
  logs that it did. `examples/miller-adaptive.toml`:
  η 2.75e-3 → 1.44e-3 → 8.15e-4 → 4.77e-4, 97 → 1069 elements, 0 paths widened.
  `theDriverRunsTheAdaptiveLoop` pins it against the library at 1.653e-16 and
  asserts separately that it refined, refined *adaptively*, brought η down, and
  kept assumption P.1.

  **It cost a prerequisite that was filed as an independent item, and that is
  the transferable part.** The estimator needs `ψ*` — four of eq. (20)'s five
  terms — and `postProcess()` refused the semi-linear path, which is the *only*
  path the driver has. So item 3's `postProcess()` bullet was not optional and
  not separate. See below for what measuring it turned up.

* **The driver's non-linear ladder is in**: Newton, and on *observed* failure
  `PicardThenNewton`, rebuilding the solver because a caught `ErrorException`
  leaves one unusable. Reactive, never predictive.
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
no black-box form at all — see item 5.

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

  What is left is: the driver builds `η` on `Potential::Raw`, one order down and
  correct, as a standing decision rather than a runtime choice; the fix is filed
  with line numbers as
  `../mfem-hdg-dev/doc/HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md`; and
  `theReconstructionIsWrongWhereTheJacobianVanishes` **asserts the defect**, so
  the suite fails the day it is fixed and says to put the estimator back.

  **The defect does not reach the solve**, measured rather than argued: `ψ_h` and
  `q_h` agree between the two paths to 1.6e-13 and 2.6e-13 over `k = 1…4` and
  three meshes, and post-processing leaves both bit-identical. The forward local
  problem takes its potential block from the HDG stabilisation on the interior
  faces, which never degenerates.
* **Rebuild `η₅` on `TransferredDatumCoefficient`** instead of excluding those
  faces. Stage-6 work, unblocked by the MFEM fix; wants its own rate.
* **Rewrite `README.md`**, which still describes a project that did not exist.

## 4. Element-local parallelism — MFEM, **in flight**

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`. Probably the largest
performance win available: twelve sequential element loops in
`darcyhybridization.cpp` over work independent by construction.

**`LineariseThenCondense` is what makes this tractable, and it now exists.** That
ordering was pursued as a convergence fix and failed as one — but it does exactly
what it says, running **zero** element-local nonlinear iterations, so under it the
batched work is fixed-size *linear* solves rather than nonlinear ones of
unpredictable iteration count. That is what `linalg/batched/` provides and what a
device wants. So the ordering is not a dead end; it is the enabler for this item,
and that is the reason to keep `setNonlinearOrdering()`.

## 5. §4.4, the current hole — a known answer that must stay out of the driver

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

## Deliberately not yet

* **GPU and cuDSS.** Correctness-testable here; this card cannot say whether it
  is worth it. Consumer FP64 is 1/32 of FP32 where datacentre parts are 1/2, so a
  local timing can invert the production conclusion.
* **Partial assembly.** Not implemented for the HDG integrators at all, and its
  case is sum factorisation, which wants tensor-product elements where meq is
  triangles. A discretisation decision, not a flag.
* **PARDISO.** 2.2× where it runs, and it does not run above `n ≈ 3000` on this
  machine's MKL. Needs a real oneMKL, and brings a reproducibility question.
* **The physics** — sonic rotation, anisotropic pressure, NetCDF profiles, MaNTA
  coupling. The actual science, all downstream of a driver. `TODO` carries each
  with what has been established.

## The standing rule

Every stage ends at a **measured** number, not at "it runs". That is what caught
the Solov'ev coefficients, the `τ` sign, the `ψ*` regression, the hybridization
ordering and the under-resolution — five times the code or the story looked right
and was not. Three corollaries this project has paid for:

* **A property measured on the easy configuration is not a property of the
  code.** Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on
  the geometry meq exists for.
* **A difficulty measured at one resolution is not a property of the problem.**
  Three sources were called stiff for months on the strength of `k = 1,
  h = 0.05`, and a whole MFEM work item was requested to fix them. They were
  under-resolved.
* **A measurement taken on this machine may not be a measurement about the
  code.** Threaded MKL decides a marginal Newton; a consumer GPU decides nothing
  about a production one.
