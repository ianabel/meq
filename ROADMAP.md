# Where meq is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27. `CLAUDE.md` is the
operational record and is authoritative on anything technical; `DRIVER-PLAN.md`
is the stage-7 design; `TODO` holds work that is understood but not scheduled.
**This file is only about order** — what to do first, what waits on what, and
what is deliberately not being done yet.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done, and **stage 7e landed: meq is a program.** `meq config.toml`
solves and writes. The suite is **17/18**, the single failure being a tripwire
whose premise stopped holding. What is left of stage 7 is the warm start's
interpolating route (§4, needs GSLIB) and wiring the two things the driver
currently refuses — the extension path and the adaptive loop.

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
fixed and covered. `HDG-ELEMENT-LOCAL-PARALLELISM.md` remains outstanding.
meq owns `src/`, `tests/` and `apps/`.

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
* **Wire the adaptive loop.** Stage 6 exists and is tested; `[adaptivity]`
  parses; the driver refuses it. It must call `setTransferredBoundary()`
  automatically on the extension path and say so in the log, until `η₅` is
  rebuilt on `TransferredDatumCoefficient`.
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

* **Rewrite the pedestal tripwire.** `pedestalNewtonFailsOnCoarseMeshesAtOrderOne`
  asserts a failure that no longer reproduces, and **its stated diagnosis is now
  known to be wrong twice over** — the prose blames Picard-vs-Newton and names
  SUNDIALS, and a later note blamed BLAS thread count. It is measuring
  **under-resolution**: `k = 1` needs `n ≥ 32`, `k ≥ 2` converges everywhere.
  Assert that instead. Do not delete it.
* **Re-run `everyNonlinearPathReachesTheSameExactSolution`** with
  `PicardThenNewton` added to the paths it pins.
* **Retire `postProcess()`'s refusal — but measure first.**
  `ReconstructFluxAndPot()` now takes a nonlinear potential lift, so the defect
  that forced the refusal may be gone. It returned ~1e15 *silently* before, so
  this is verified by a rate, never by reading the source.
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
