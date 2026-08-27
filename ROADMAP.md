# Where meq is, and what to do next

Written 2026-08-26. `CLAUDE.md` is the operational record and is authoritative
on
anything technical; `DRIVER-PLAN.md` is the stage-7 design; `TODO` holds work
that is understood but not scheduled. **This file is only about order** — what
to
do first, what is waiting on what, and what is deliberately not being done yet.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done; stage 7 is half done — `setInitialGuess()` and `BoundaryShape`
landed, the driver did not. The suite is **12/13**, the single failure being a
tripwire that fires on BLAS rounding rather than on a defect. **meq still cannot
be run by anyone**: nothing in `src/` writes a file, so it is reachable only
through `ctest`. That has been true since stage 1.

## The finding that reorders everything

meq's Newton fails on the stiff GS-2 sources. Three explanations were tested and
two are wrong:

* **Not Newton.** `refs/MFEM-GS-Newton.pdf` — Serino, Tang, Tang, Kolev &
  Lipnikov, built on MFEM — reports Newton robust on cases where "conventional
  Picard-based solvers fail to converge".
* **Not the boundary condition.** The same paper calls its earlier
  *fixed*-boundary solver "significantly easier" than the free one. meq does the
  easier problem and finds it harder.
* **It is the hybridization ordering.** meq condenses first and linearises
  second, so eliminating flux and potential on an element is itself a nonlinear
  solve — one per element per residual evaluation, none globalised, any one of
  which poisons the residual. Nguyen, Peraire & Cockburn, who defined the
  method,
  do it the other way: §2.6 of `refs/HDG-NPC-2.pdf` applies Newton to the full
  `(q, u, û)` system and hybridizes *that*, so every local operation is a linear
  solve.

Confirmed by control: Picard on meq's own linear path — same mesh, spaces and
`τ`, but `F` frozen so the local solves are linear — converges on the case
Newton
fails. Details and numbers in `CLAUDE.md`, *Why meq's Newton struggles*.

## Division of labour

**The MFEM tree is being worked by another agent.** Three plans sit in
`../mfem-hdg-dev/doc/`, all written from meq's measurements and none
implemented.
meq owns everything in `src/`, `tests/` and `apps/`; the items below say which
side each belongs to and what meq must do when an MFEM item lands.

---

## 0. Three MFEM branches, one integration branch — standing maintenance

meq needs work from `gf-hdg-subdomains-dev`, `direct-solver-symbolic-reuse` and
`gf-hdg-linearise-first`, none of which contains the others.
**`meq-integration`** is their merge and is what `../mfem/install` is built
from. It is local only, never pushed, and re-created from the three whenever any
of them moves — so nothing may be committed directly to it. `CLAUDE.md` carries
the recipe and the one recurring conflict.

That is now a standing cost of this arrangement rather than a one-off, and it
falls on meq's side.

## 1. Linearise, then condense — MFEM, **landed, and a bug reported against it**

`../mfem-hdg-dev/doc/HDG-LINEARISE-THEN-CONDENSE.md`, implemented on
`gf-hdg-linearise-first` as
`SetNonlinearOrdering( NLOrdering::LineariseThenCondense )` — **off by default**,
so meq must opt in. `setNonlinearOrdering()` does, and the default is unchanged.

**One bug reported and fixed, a second reported and open.** Both against the
nonlinearity on the potential mass form — `Mnl_p`, where a semi-linear source
goes — which the existing unit test does not cover, driving the *block*
nonlinear form instead.

1. **Fixed** by `c5cac09e2f`: the residual was not a function of the trace.
   `Mult(x)`, `GetGradient(x)`, `Mult(x)` gave two different answers for the
   same `x`. Verified — the demonstrator now returns exactly zero.
2. **Open**: `GetGradient` is not the derivative of `Mult`. The error scales as
   the square of the nonlinearity — 5e-9, 5e-7, 5e-5 for `c` = 1, 10, 100 on a
   `c ψ²` source — and is independent of the differencing step, 5.465e-05 at
   `c = 100` across four orders of `h`, where `CondenseThenLinearise` traces the
   round-off curve of an exact Jacobian throughout.

Both in `../mfem-hdg-dev/doc/LINEARISE-FIRST-RESIDUAL-BUG.md`, with
demonstrators beside it.

So this item is not closed and meq keeps the old ordering as its default. The
acceptance list below is what to work when it is.

**Why it is first.** meq chose Newton deliberately and pays for it in every
source and profile it accepts — `dFdPsi` is mandatory, `Prime` is mandatory, a
profile that cannot differentiate itself cannot be used, and the assembled
Jacobian is checked against a finite difference of the residual. That cost is
paid *in order to have* Newton. Under the present ordering the cost stands and
the benefit is absent on exactly the problems that motivated it.

The gap is not a preference between solvers. On Example 5, all three paths
reaching identical discretisation error:

| | iterations |
|---|---|
| Newton | 4 |
| Anderson-accelerated Picard | 19 |
| damped Picard | 97 |

Picard converges linearly; Anderson is superlinear at best and never quadratic.

**meq should release with the correct algorithm, not with it added afterwards.**
A code earns its use cases by being the right choice from the start.

**What meq does when it lands:**

* Re-run `everyNonlinearPathReachesTheSameExactSolution`, which pins Newton,
  Anderson-Picard and damped Picard to the same closed-form answer. It is the
  check that would catch the new ordering solving a subtly different problem.
* Re-measure the four GS-2 sources, and pose them with the papers' own
  homogeneous data if a guess now suffices — `setInitialGuess()` alone did not.
* Re-measure the pedestal and the current hole, and rewrite or retire
  `pedestalNewtonFailsOnCoarseMeshesAtOrderOne` on the evidence.
* Keep the Picard paths. Not as a fallback anyone should reach for, but because
  they are the independent check on Newton.

## 2. The driver — meq, and can start now

`DRIVER-PLAN.md` §3–5: stage 7c output, 7d warm start, 7e the driver itself.

Independent of item 1 — `setGlobalisation()` already makes the nonlinear path a
runtime choice, so the driver does not care which solver runs. What item 1
changes is what meq *ships with*, not whether the driver can be written.

This is the only item that changes what meq **is**. Everything else improves a
thing nobody can execute.

* **7c** — mesh and grid functions, plus a NetCDF file carrying `ψ` and both
  components of `B_poloidal` on a masked `(R, Z)` grid inside `Γ`. `B` is a
  relabelling of the solved flux, `B_R = −q_z`, `B_Z = +q_r`, with no derivative
  taken — the payoff for the mixed method. **Verify by finite difference against
  exact Solov'ev**; three of this project's four convention questions went
  against the derivation.
* **7d** — warm start, with the NetCDF grid doubling as the interchange format.
* **7e** — `meq config.toml`, exit codes, residual history, the adaptive loop.

**Prerequisite:** `MFEM_USE_EXCEPTIONS` on the next `../mfem/install` rebuild.
§5's exit code 2 for a failed solve is unimplementable while MFEM aborts the
process, and the current hole is exactly the case that aborts.

## 3. Element-local parallelism — MFEM, after item 1

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`.

Probably the largest performance win available: twelve element loops in
`darcyhybridization.cpp`, all sequential, over work that is independent by
construction.

**After item 1 deliberately, not merely conveniently.** Under the present
ordering it would batch *nonlinear* solves of unpredictable iteration count;
under NPC's it batches fixed-size *linear* ones, which is what `linalg/batched/`
provides and what a device wants. Building it first designs it around the
compromise.

## 4. Symbolic factorisation reuse — MFEM, whenever

`../mfem-hdg-dev/doc/DIRECT-SOLVER-SYMBOLIC-REUSE.md`. 22–24% of every Newton
step's linear cost, thrown away because `UMFPackSolver::SetOperator` declares
`Symbolic` as a local. Self-contained, no design decisions, in-tree precedent in
`CuDSSSolver::SetReorderingReuse`. Good work to slot between larger items.

## 5. Hygiene — meq, small and worth doing alongside the driver

* **Rewrite the pedestal tripwire.** It asserts a failure that now depends on
  BLAS thread count. Its prose blames Picard-vs-Newton and names SUNDIALS, which
  is the wrong diagnosis. It should assert something reproducible.
* **Retire `postProcess()`'s refusal**, now that `ψ*` is measured working
through
  Newton — in the same change that adds the measurement, not before.
* **Rebuild `η₅` on `TransferredDatumCoefficient`** instead of excluding those
  faces. Stage-6 work, unblocked by the MFEM §4 fix; wants its own rate.
* **Rewrite `README.md`**, which still describes a project that did not exist.

---

## Deliberately not yet

* **GPU and cuDSS.** Correctness-testable here; this card cannot say whether it
  is worth it. Consumer FP64 is 1/32 of FP32 where datacentre parts are 1/2, so
  a local timing can invert the production conclusion. Development platform, not
  a benchmark.
* **Partial assembly.** Not implemented for the HDG integrators at all, and its
  case is sum factorisation, which wants tensor-product elements where meq is
  triangles. A discretisation decision, not a flag.
* **PARDISO.** 2.2× where it runs, and it does not run above `n ≈ 3000` on this
  machine's MKL. Needs a real oneMKL, and brings a reproducibility question.
* **The physics work** — sonic rotation, anisotropic pressure, NetCDF profiles,
  MaNTA coupling. This is the actual science and all of it is downstream of a
  driver; MaNTA coupling additionally needs 7c's NetCDF layer to exist. `TODO`
  carries each with what has been established.

## The standing rule

Every stage ends at a **measured** number, not at "it runs". That is what caught
the Solov'ev coefficients, the `τ` sign, the `ψ*` regression, and the
hybridization ordering — four times where the code looked right and was not. Two
corollaries this project has paid for:

* **A property measured on the easy configuration is not a property of the
  code.** Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on
  the geometry meq exists for.
* **A measurement taken on this machine may not be a measurement about the
  code.** Threaded MKL decides a marginal Newton; a consumer GPU decides nothing
  about a production one.
