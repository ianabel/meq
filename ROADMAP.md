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
landed, the driver did not. The suite is **17/18**, the single failure being a
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
* **Not the hybridization ordering either — this one was believed, and it is
  now falsified.** The theory was that condensing first makes every element
  elimination its own nonlinear solve. `NLOrdering::LineariseThenCondense`
  removes those entirely (`GetNumLocalNLIterations()` reads 0) and the stiff
  cases get *worse*, not better. See item 1.
* **It was mostly under-resolution.** Three of the four sources are not stiff at
  all. Raw Newton solves §4.2, §4.3 and §4.5 in 7 to 17 iterations once
  resolved, and refining `h` or raising `k` cures them independently — §4.2 at
  `k = 3, n = 48` takes **7**, and `k = 2` converges on the same `n = 16` mesh
  where `k = 1` takes 42. The whole "stiff sources" narrative came from
  benchmarking at one under-resolved point.
* **What is left is the coarse mesh, and §4.4.** Newton still needs to begin
  inside its basin, and on a coarse mesh the Dirichlet datum is not; that is
  what `PicardThenNewton` is for, and it matters because an adaptive run has to
  solve *before* it has an estimator. §4.4, the current hole, fails at every
  order and mesh tried and is a trivial-branch problem, not a resolution or
  iteration one.

Details and numbers in `CLAUDE.md`, *Why meq's Newton struggles* and *Picard,
then Newton*.

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

## 1. Picard, then Newton — meq, **done**, and the ordering theory is dead

**The prediction this roadmap was built on is falsified.** It said meq's Newton
fails on the stiff GS-2 sources because meq condenses first and linearises
second, making every element elimination its own nonlinear solve, and that
`NLOrdering::LineariseThenCondense` would fix it. That ordering has now landed
on `gf-hdg-linearise-first`, the two bugs meq reported against it are fixed, and
`GetNumLocalNLIterations()` returns **0**, which is what says it genuinely takes
effect. It does not help:

| case | condense-first | linearise-first |
|---|---|---|
| Example 5, similarity | ok, 3 it | ok, 3 it |
| §4.2 pedestal `k = 1`, three meshes | ok — 31/7/5 it | **fails at 60, all three** |
| §4.3 barrier | fails / ok | **aborts** |

Identical where it is easy, worse or fatal where it is not. The obvious
mechanism — `NewtonSolver` testing a residual built one iterate stale — was
checked too, and forcing relinearisation at the residual's own argument changes
nothing. **The element-local nonlinear solves were never the cause.** The control
that pointed at them differed in globalisation as well as in local linearity, and
KINSOL on the outer loop only drives the local solves harder: on the pedestal a
line search spends 1.38M local iterations to fail where plain Newton spends 133k
to succeed.

**`Globalisation::PicardThenNewton` is implemented**, and it is what reaches the
hard cases *without refining them*: §4.3 at `k = 1, h = 0.05` finishes in **4
Newton iterations at observed order 2.01** where Newton alone fails at 60,
agreeing with Picard's own answer to 4.7e-10.
`picardThenNewtonRecoversQuadraticOrder` is the regression.

**But refinement reaches the same three cases**, so this is the coarse-mesh
route, not the general remedy. Full tables in `CLAUDE.md`, *Why meq's Newton
struggles* and *Picard, then Newton*.

So **meq releases with Newton, and with Picard as its globalisation** — which is
the answer to "release with the optimal algorithm". Newton is the algorithm; it
just needs to be started somewhere reachable, exactly as CEDRES++ and Serino et
al. start theirs. Picard is not pasted on afterwards, it is the globalisation
Newton was always missing.

**Keep `setNonlinearOrdering()` anyway.** It costs nothing, it is the canonical
NPC ordering, and it is the only way to run this discretisation with no
element-local nonlinear solves at all. It is no longer blocking anything.

**What is left on this path:**

* **§4.4, the current hole, is the one genuinely unsolved case.** It fails at
  every order and every mesh tried — including `k = 3, n = 48`, at 1.8M
  element-local iterations — and under Picard and the handoff too. Refinement
  does nothing, which points at the trivial branch (`F(r, 0) = 0`, so `ψ ≡ 0`
  solves the homogeneous problem) rather than at resolution or globalisation.
  Continuation in the source amplitude is the untried idea.
* **`MFEM_USE_EXCEPTIONS`** — needed before the hole can even be *reported* as a
  failure rather than killing the process, and a prerequisite for the driver's
  exit code 2. This is now the highest-value MFEM-side item.
* Re-run `everyNonlinearPathReachesTheSameExactSolution` with the new path added.
* Rewrite `pedestalNewtonFailsOnCoarseMeshesAtOrderOne` to say it is measuring a
  knife edge, per `CLAUDE.md`. Do not delete it.

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
