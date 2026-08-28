# Where meq is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27. `CLAUDE.md` is the
operational record and is authoritative on anything technical; `DRIVER-PLAN.md`
is the stage-7 design; `TODO` holds work that is understood but not scheduled.
**This file is only about order** — what to do first, what waits on what, and
what is deliberately not being done yet.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done; stage 7 is half done — `setInitialGuess()`, `BoundaryShape` and
the 7c writers landed, the driver did not. The suite is **16/17**, the single
failure being a tripwire whose premise stopped holding. **meq still cannot be run
by anyone**: nothing reachable from a command line writes a file, so the solver
is exercised only through `ctest`. That has been true since stage 1 and is now
the only thing standing between meq and being useful.

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

## 1. The driver — meq, and nothing is in its way

`DRIVER-PLAN.md` §3–5: **7c output (done)**, 7d warm start, 7e the driver.

**This is the only item that changes what meq _is_.** Everything else improves a
thing nobody can execute. It was previously sequenced behind the nonlinear work;
that work is finished and the dependency is gone.

* **7d** — warm start, with 7c's NetCDF grid doubling as the interchange format.
* **7e** — `meq config.toml`, exit codes, residual history, the adaptive loop.

**The nonlinear path the driver ships:** Newton, falling back to
`PicardThenNewton` on **observed** failure. A *reactive* ladder — never a
predictive one. Nothing may be inferred from `F` about which solver to run,
because nothing can be: the ratio `max|∂F/∂ψ|/λ₁` is computable black-box but the
pedestal converges at 7 where the hole fails at 26, which is two points and not a
threshold. **Continuation must not go in**, for the stronger reason that it has
no black-box form at all — see item 5.

**Prerequisite: `MFEM_USE_EXCEPTIONS`** on the next `../mfem/install` rebuild.
§5's exit code 2 for a failed solve is unimplementable while MFEM aborts the
process, and §4.4 is exactly the case that aborts. It is a `config.hpp` change,
so a full rebuild.

## 2. Use the symbolic factorisation reuse that landed — meq, small

`UMFPackSolver::Symbolic` is now a **member**, with `SetReuseSymbolic()`,
`GetNumSymbolic()` and `GetNumNumeric()`. meq calls none of them, so it is still
throwing away **22–24% of every Newton step's linear cost** at all three
`UMFPackSolver` sites in `GradShafranov.cpp`.

Cheap, no design decisions, and it comes with its own acceptance criterion: assert
`GetNumSymbolic() == 1` across a multi-step Newton solve, which is a *measured*
number in the house style rather than a timing.

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

## 4. Element-local parallelism — MFEM, and now unblocked

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
