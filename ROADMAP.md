# Where meq is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27, 2026-08-29 and 2026-09-01.
`CLAUDE.md` is the operational record and is authoritative on anything technical;
`DRIVER-PLAN.md` is the stage-7 design, `FREE-BOUNDARY-PLAN.md` and
`FLOW-PLAN.md` are the two designs not yet started; `TODO` holds work that is
understood but not scheduled. **This file is only about order** — what to do
first, what waits on what, and what is deliberately not being done yet.

## Nothing is red

**`HighBetaConvergence` was the one failing test and it is green.** Profiles
specified in NORMALISED flux need `ψ_ax` inside the residual, where the Jacobian
can see the non-local terms it contributes;
`meq::NormalisedSource` is the interface and
`setSource( NormalisedSource &, double )` closes the pair by a **bordered
Newton**, one factorisation and two backsolves per step. Under NPC two of the
three border quantities are exact rather than differenced.

The account of it is in `CLAUDE.md` under *Newton, and the obligation it
creates*, including the three cheaper answers that were measured and killed on
the way — a fixed `ψ_ax` makes the profile inert, an outer iteration on
`ψ_ax ← max ψ` has a pole beside its own fixed point, and at fixed `ψ_ax` there
is a second solution that is not the equilibrium. The test that can *see* the
missing terms is `Normalisation::Decoupled`, which does not converge slowly so
much as not move at all: 8.31e-2 to 8.24e-2 in fifteen iterations, against
4.4e-15 in four.

**What this unblocks is free boundary**, where `ψ_ax` and `ψ_bnd` are both
unknowns: `ψ_bnd` is a second border row and column of exactly the same shape.

**What is left over from it**, and neither is on the critical path: the
production source `meq::NormalisedMHDSource` is **not reachable from a TOML
file** — item 1 — and the bordered path accepts `Globalisation::None` only,
refusing the others loudly rather than quietly doing something else.

**And it turned up an MFEM finding that is still not filed anywhere**, which is
item 5.
## So what is next

With nothing red, the order is what it was, plus a second planned physics item:

1. **Finish the driver** — item 1. Two things left in it: wire the normalised
   source through `Config`, and the GSLIB warm start. **The first is now a
   prerequisite for two things rather than one**: the normalised plumbing gets
   written once and both the high-β source and a rotating one use it.
2. **Then the two planned extensions, and they are independent of each other.**
   * **Toroidal flow**, `FLOW-PLAN.md` — item 9, and **done**, FL-0 to FL-8.
     `refs/RotatingGK.pdf` (136), closed by (96) and (97). A change to `F`
     alone: the operator, the discretisation and the whole curved-boundary and
     adaptive apparatus were untouched, and `Source`'s signature already
     carried the `r` it needs.
   * **Free boundary**, `FREE-BOUNDARY-PLAN.md` — item 8. HDG on a polygonal
     subdomain coupled at a distance to an exterior operator on a semicircular
     artificial boundary, by `refs/CouplingAtADistance.pdf` — whose own
     reference [5] is meq's stage 5, so the coupling is the extension technique
     already in the tree with the datum unknown instead of zero. Staged FB-A and
     FB-0 to FB-6.

   **Flow is the cheaper of the two** and shares more with what exists: it needs
   no new geometry, no cut quadrature and nothing from the other tree. Free
   boundary is the more valuable and the more structural. Neither blocks the
   other, and both want item 1's normalised plumbing first.
3. ~~**Hygiene** — item 3 — alongside either. `README.md` is the oldest debt in
   the tree and still describes a project that did not exist.~~ **Item 3 is
   closed, 2026-09-02**: the README is true, and `docs/` is a Sphinx manual on
   Read the Docs.

Items 4 and 5 are the other tree's to act on; item 6 is meq's and is small and
not urgent; item 7 is a closed investigation kept for its answer; item 10 — the
fixed-`q(ψ)` solver — is wanted but undesigned and should wait on item 9.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done, and **stage 7 is finished: meq is a program that solves on a
curved boundary, refines its own mesh, and restarts from a previous answer in one
Newton step.** The suite is **28 of 28**, and it is green rather than
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
  **Confirmed twice**: that mode is now deleted, and `NPC` — the canonical
  version of the same idea, and meq's default — reproduces the verdict.
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
against**, which is the only test of "landed" meq can apply. **Ask `git` which
documents exist, not `ls`**: all three open ones live on `gf-hdg-linearise-first`
alone, so a listing taken while that tree sits on `gf-hdg-dev` shows an almost
empty `doc/` and means nothing. `CLAUDE.md` has the branch-by-branch table and
the mistake that produced it.

| filed | state |
|---|---|
| `HDG-LINEARISE-THEN-CONDENSE.md` | landed, then **retired** with the mode itself — on backup refs only now. `setNonlinearOrdering()` is gone; meq's default is `NPC` |
| `DIRECT-SOLVER-SYMBOLIC-REUSE.md` | landed — `SetReuseSymbolic()` is on — and **retired**, on no branch at all |
| `HDG-NPC-GLOBALISATION-FROM-MEQ.md` | **filed 2026-08-31 and answered the same day** (`af82d42b14`). Not a defect report: two §6 claims withdrawn on meq's evidence, meq's own account of the mechanism corrected, and a defect found in the reference implementation meq had copied |
| `HDG-DEFECTS-FROM-MEQ.md` | **still on `gf-hdg-dev` and `gf-hdg-subdomains-dev`**, deleted only on the symbolic-reuse line — which is why that merge conflicts modify/delete. Three of its four are verifiably closed — one fixed, one fixed, one withdrawn as not a defect; the fourth, `ComputeHDGFaceEnergy()` ignoring an installed stabilisation, meq has not re-measured and does not use. `CLAUDE.md` has the breakdown |
| `HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md` | **landed and retired** — the fix is *"The postprocessing closes on the element average, always"*, and the document is on no branch |
| `HDG-ELEMENT-LOCAL-PARALLELISM.md` | open, on `gf-hdg-linearise-first`, and meq has seen no change |
| `HDG-BEM-COUPLING-FROM-MEQ.md` | **filed 2026-08-29**, for free boundary — open, on `gf-hdg-linearise-first`. See item 8 |

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
| `gf-hdg-linearise-first` | **`DarcyNPCOperator` / `DarcyNPCSolver`** — meq's ordering. `NLOrdering::LineariseThenCondense` was on this branch and is **deleted** |
| **`gf-hdg-dev`** | **the reconstruction fix — *"The postprocessing closes on the element average, always"*** |

**`meq-integration`** is their merge in `../mfem/mfem-src` and is what
`../mfem/install` is built from. Local only, never pushed, re-created whenever any
of the four moves — so nothing may be committed directly to it.

**Four branches, but TWO merges, and this file said three.** As of 2026-08-30
`gf-hdg-dev` is an **ancestor** of both `gf-hdg-subdomains-dev` and
`gf-hdg-linearise-first`, so it arrives with the base and needs no merge of its
own. The earlier record here — that none of the four contains any other — was
true when written and is not now. **The topology is the other tree's to change,
so re-verify rather than trusting either statement.**

What has not changed is *why* `gf-hdg-dev` is on the list. Dropping it silently
loses the reconstruction fix, `ψ*` goes back to being a different function
wherever `∂F/∂ψ` vanishes, and the adaptive loop goes back a full order — with
no error and a green suite until
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` catches it. It is
free today rather than unnecessary.

**Run the containment loop; it has already fired once.** `CLAUDE.md` carries it
with the recipe. Against `meq-integration` at `fa65a2f932` on 2026-09-01 three of
the four report contained and `gf-hdg-linearise-first` does not, because that
branch has advanced 17 commits since the merge. **Read that as "the branch has
moved", not "the merge failed"** — what meq needs from it is in the tree. It
says a re-merge is available.

A standing cost of the arrangement, and it falls on meq's side.

## 1. Finish the driver — meq

**`DRIVER-PLAN.md` §3–5 is done bar two bullets.** `meq config.toml` parses,
solves and writes the equilibrium three times over — `.mesh` plus `.gf`, VTK at
the polynomial degree of the solve, and `ψ` and **B** on a uniform `(R, Z)` grid
— with exit codes 0/1/2/3. The curved boundary works through it, the adaptive
loop works through it on both paths, and the non-linear ladder is in: Newton,
falling back to `PicardThenNewton` on **observed** failure. Every one of those
is pinned against the library rather than against a closed form, for the reason
recorded beside `examples/soloviev-nstx.toml`. `CLAUDE.md` has the numbers.

**What is left, in the order it is worth doing:**

* ~~**Wire `meq::NormalisedMHDSource` through `Config` and `SourceFactory`**~~ —
  **DONE**, by FL-8, and written once for both as this bullet predicted.
  `Normalised = true` plus a `PsiAxis` guess reaches it, `makeNormalisedSource()`
  is the door, and `makeSource` **throws** on a normalised config rather than
  quietly solving a different problem. `examples/rotating-normalised.toml`
  exercises the path.
* ~~**7d, the interpolating warm start**~~ — **DONE 2026-09-02**, and it was
  the wiring this bullet said it was. Cycles after the first now start from the
  previous one interpolated onto the refined mesh, and a stored guess on a
  different mesh transfers instead of throwing. Measured on a nonlinear source
  over three cycles: **4, 2, 2 Newton iterations against 4, 4, 4 cold — a third
  of the work — with the final L2 unmoved to 2.5e-12**. It is unmeasurable on
  meq's own adaptive examples, whose Solov'ev source takes one Newton step
  whatever it starts from, which is why the measurement is
  `carryingTheAnswerAcrossCyclesCutsTheWork` rather than an example. The record
  below is kept because the correction it makes is the transferable part:
  `src/meq/WarmStart.{hpp,cpp}` is `meq::FieldTransfer` on
  `mfem::FindPointsGSLIB`, it is in `MEQ_CORE_SOURCES`, `MFEM_USE_GSLIB` is
  `YES` in the install, and `WarmStartConvergence` is a registered ctest whose
  `fullOrderCarriesMoreThanAStructuredGrid` **is** the acceptance measurement
  `DRIVER-PLAN.md` §4 names, run over a 65/129/257 grid sequence. What is left
  is one wiring job: `apps/meq.cpp` still throws *"the interpolating one needs
  GSLIB"* on a mesh-count mismatch, and adaptive cycles after the first start
  cold for the same reason. Deleting that throw and carrying the answer across
  cycles is hours, not days, and it is the best value-per-hour item on this
  list. (`src/meq/Sampler.hpp` and `src/meq/Config.hpp` both asserted that GSLIB
  was off, which was false; **corrected 2026-09-02** while item 11 was in that
  file, so the only thing left here is the driver wiring itself.)

**The nonlinear path the driver ships is reactive and must stay so.** Nothing
may be inferred from `F` about which solver to run, because nothing can be: the
ratio `max|∂F/∂ψ|/λ₁` is computable black-box but the pedestal converges at 7
where the hole fails at 26, which is two points and not a threshold. A second
candidate detector — the first Newton step making the residual worse — has since
been measured and is *anti*-correlated. **Continuation must not go in** either,
for the stronger reason that it has no black-box form at all; see item 7.

**The only thing the driver still refuses is `[boundary] Type = "exact"`**,
which needs a closed form `meq::Source` does not carry. It exits 1 with an
explanation rather than approximating.
## 2. ~~Symbolic factorisation reuse~~ — **done**

`SetReuseSymbolic()` is on for the Newton and Picard paths and off deliberately
on the linear one, which factorises once.
`theSymbolicAnalysisIsReusedAcrossNewtonSteps` asserts one analysis against one
factorisation per iteration — a count and not a timing, and the only thing that
could notice the reuse lapsing, since a lapse costs speed and nothing else.

## 3. Hygiene — meq, alongside the driver

**Four of the six bullets that stood here are done**, and what they found is in
`CLAUDE.md` rather than here: the pedestal tripwire now asserts the two
refinement cures rather than a knife edge that threaded-MKL rounding decides;
`everyNonlinearPathReachesTheSameExactSolution` runs four paths to the same L2;
and `postProcess()`'s refusal is retired, which turned out to be the
prerequisite for the adaptive loop rather than a separate task.

**That last one is worth one paragraph, because it changed what the driver
does.** The estimator needs `ψ*` in four of eq. (20)'s five terms, and MFEM's
reconstruction was silently wrong **per element** wherever `∂F/∂ψ` vanishes — a
pure Neumann problem whose mean-value regularisation was being skipped, so a
singular matrix was factored. A whole-domain norm could not see it: at an
eighth of the domain dead, individual elements were 20× wrong while the norm
read 1.87. It is fixed upstream, unconditionally, as *"The postprocessing closes
on the element average, always"*; the driver's loop is back on `ψ*` and reaches
η = 6.87e-5 on 449 elements where the degraded estimator needed 1069 to reach
4.77e-4. The test flipped from asserting the defect to asserting the behaviour
and went green on its own, which is the testing stance working as intended.

~~**One bullet is still open.**~~ **None are, as of 2026-09-02.**

* ~~**Rebuild `η₅` on `TransferredDatumCoefficient`**~~ — **DONE 2026-09-02.**
  `setTransferredBoundary()` takes the datum as a second argument and the faces
  stay in, compared against the `φ_h` actually imposed rather than against the
  zero standing in for it. `η₅` goes from 4.07e-1 to **9.6e-5** on the coarsest
  extension mesh and converges at **2.78** against 0.40, so the term is two
  percent of `η₁` rather than four orders larger. The adaptive loop is unchanged
  — 97 → 254 → 342 → 449, `η` 4.7352e-04 → 6.8668e-05 — which is the point:
  boundary elements now have an `η₅` contribution instead of none, and nothing
  that already worked moved. `theTransferredDatumRestoresEtaFive` keeps the
  pinned column as its control, and
  `theTransferredDatumReproducesTheImposedCondition` checks `φ_h` against the
  exact `ψ` it transfers — which is what catches the one real hazard, that the
  lift must be given the raw `−q` block and not `flux()`, since the wrong one
  returns `−ψ` rather than `ψ`.
* ~~**Rewrite `README.md`**, which still describes a project that did not exist.
  The oldest debt in the tree, and the one a new reader hits first.~~ —
  **DONE 2026-09-02**, `57046c7`. Four things were wrong with it and none was a
  matter of taste: it said the driver was "not yet ported, so there is currently
  no way to run MEQ except through its test suite"; it gave the binary as
  `./build/apps/meq`, which is not where it is built; it gave `MFEM_DIR` as
  defaulting to `../mfem-hdg-dev`, which is the tree meq must **not** build
  against, rather than `../mfem/install`; and it kept a paragraph insisting that
  "the environment variable is not optional" **without naming a variable**,
  because the `MKL_THREADING_LAYER=GNU` it was about had been deleted from the
  ctests and the driver as inert.

  **That last one is the transferable part.** Removing a variable from the code
  and from every test left the sentence *about* it behind, pointing at nothing —
  the same species of residue as a stale measurement, and harder to notice
  because it reads as a warning rather than as a claim.

  **And `docs/` is a different thing now**: a Sphinx tree published to Read the
  Docs at <https://meq.readthedocs.io>, with the pre-Sphinx LaTeX manual moved
  intact to `docs/manual/` and keeping its own Makefile. `CLAUDE.md`'s Layout
  block still calls `docs/` "the LaTeX manual, which predates the port", which is
  now the stale line about documentation rather than this one.
## 4. Element-local parallelism — MFEM, outstanding

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`, on
`gf-hdg-linearise-first`. Twelve sequential element loops in
`darcyhybridization.cpp` over work independent by construction.

**IT IS NOT "PROBABLY THE LARGEST PERFORMANCE WIN AVAILABLE", WHICH IS WHAT
THIS FILE SAID, AND UPSTREAM HAS MEASURED IT.** `d7ea90a538`, *"Element-local
parallelism in NPC: the structural win is worth 6%"*, breaks an NPC step down
per phase over six steps at four `(n, k)`: integrator-bound loops **59–63%**,
integrator-free loops **5.4–6.2%**, trace solve **31–35%**. Their conclusion is
that the loops threadable *today* are under 6% of the step and flat in both mesh
size and order, so Amdahl caps any gain there until the integrators are made
thread-safe. meq should stop quoting this as the big win.

**And §1 of it has LANDED and meq is not using it**: `LocalFactorMode`,
`SetLocalFactorMode` and `CanBatchLocalFactor` are in the *installed*
`darcyhybridization.hpp`, and `LocalFactorMode` appears nowhere in `src/`,
`tests/` or `apps/`. See item 6 for what it does and does not buy.

**AND IT IS UNBLOCKED, MORE CLEANLY THAN WHEN IT WAS FILED.** meq's default is
now `NonlinearOrdering::NPC` — `mfem::DarcyNPCOperator` over the full
`(q, ψ, ψ̂)` system — and `GetNumLocalNLIterations()` reads exactly **0** under
it against 3644 for the condensation, asserted by
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally`. So the batched work
really is fixed-size **linear** solves, which is what `linalg/batched/` provides
and what a device wants.

**The `ψ_ax` borders are not an obstacle either.** They were believed to require
condense-first and do not: measured 2026-08-30, both orderings reach the same
`ψ_ax` to every digit printed, and under NPC the border row is exactly `−e_j` and
the corner exactly `1`, neither differenced.

**A LONG HISTORY STOOD HERE AND IS DELETED.** It was the investigation of
`NLOrdering::LineariseThenCondense` — a mode MFEM has since retired as *"a
condensation in disguise"* — including its frozen-correction diagnosis, its 149%
cross-linearisation measurement and the case-by-case parity table. It describes
nothing meq runs and is recoverable from git. Three things in it are still live
and are kept in `CLAUDE.md` rather than here: the **falsified prediction** that
an ordering would fix a stiff source, which NPC has now falsified a second time;
the **cross-linearisation technique** for detecting hidden state in a reduced
operator; and the **live defect in the ordering meq still ships as the backup**,
where the assembled gradient disagrees with a central difference by 100× once the
element-local solves hit their iteration cap. See *The NPC port* and *Why meq's
Newton struggles*.

**Condense-first is retained only so the two can be measured against each
other**, which is what turns "this fails" into "this fails and the other does
not" — and on stiff under-resolved meshes it is still the one that works.

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

## 6. Threaded MKL costs 140x, and PARDISO's scaling is out of reach — meq

**The link-line straddle this item was about is FIXED** (2026-09-01: meq builds
its own SuiteSparse, so exactly one MKL is loaded), **and fixing it exposed a
140x regression in the suite.** The stray sequential layer had been making all
of MKL sequential, which is why every timing this project recorded before
2026-08-30 was fast.

What threading MKL actually costs, with the columns separated:

* **Not UMFPACK.** `UMFPackSolver::SetOperator` degrades about 40% across the
  whole thread range, never more.
* **`ComputeH()`'s element-local dense LU**, through LAPACK, on blocks of order
  10–30. `k = 2` untouched; **`k = 3` forty-fold worse** at two threads.
* So `MKL_NUM_THREADS=1` is on every ctest.

**And that pins PARDISO's own scaling out of reach.** PARDISO beats UMFPACK
1.50x on setup even sequentially and about 1.9x more at 8 threads — but
`MKL_NUM_THREADS` is process-wide, so buying that means paying 40x on assembly.
`setTraceSolver()` landed and is not the whole job: what is left is either
`mkl_set_num_threads_local()` around the trace solve, or the element-local
factorisation off threaded MKL. Neither is done,
and it is item 0 of `CLAUDE.md`'s *What to do, in order of value*.

**Two things settled and no longer worth an item.** PARDISO's `n ≈ 3000`
ceiling was Debian's `intel-mkl` 2020.4.304 and nothing about the method; and
the trace solver is a run-time choice, `setTraceSolver()` picking among
`UMFPack`, `Pardiso` and `CuDSS`, with UMFPack still the default because
oneMKL's licence is not everybody's to accept. `TODO`'s PARDISO entry carries
the one question that is left, which is reproducibility under threading.

`CLAUDE.md`'s *Threading, measured* has the tables.
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
  border column today; free boundary adds `ψ_bnd` and the `N` coefficients.
  **Under NPC two of the three border quantities are exact rather than
  differenced** — `b = −e_j` and `d = 1`, because `ψ` is an unknown — and the
  coefficients' column is a raw block, since the NPC residual is unreduced and
  has no condensation for a rectangular block to survive.
* **The order of work is FB-A, then FB-0 to FB-6**, and FB-1 — vacuum with coils, whose
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

**WHAT THIS NEEDS FROM THE OTHER TREE IS NOW NOTHING, FOR FB-0 THROUGH FB-3, AND
THIS FILE SAID OTHERWISE.** It claimed §2 of `HDG-BEM-COUPLING-FROM-MEQ.md` — two
rectangular integrators — was blocking. It is not, and the reason is the NPC
port: the datum's data half is an **essential trace value**, not a weak form, so
its block is `ProjectBdrCoefficient` against the `PathTraceCoefficient` that
already exists; and the transmission block is reachable from `TransferPath::
Endpoint` and `ElementExtension::TransformBack`, both public, in about forty
lines of meq. §3, auxiliary globally-coupled unknowns, is an optimisation over
`N + 2` backsolves against one factorisation — the cost meq's `ψ_ax` border
already pays. Plan §6.4 is the per-stage table.

**The one real gap is in FB-4 and it is not the one the plan predicted**: MFEM
*does* have cut-element quadrature — `mfem::MomentFittingIntRules` in
`fem/intrules_cut.hpp`, `MFEM_USE_LAPACK`-gated and installed — and what it does
not have is the **sensitivity of a cut rule to the level set**, which is the half
CEDRES++ actually names. Not blocking either: difference it per cut element, or
accept an inconsistent Jacobian on cut elements and measure the cost in Newton's
observed order.

`NORMALISED-LINEARISE-FIRST.md` was the design for meq's half of this under
`LineariseThenCondense`. **Both the mode and the design are deleted** — the
mechanism it describes has no code path left to run on, and under NPC the
problem it solved does not arise. Plan §4.1; git has it if it is ever wanted.

---

## 9. Toroidal flow — meq, DONE, FL-0 to FL-8

**`src/meq/RotatingSource.{hpp,cpp}` solves `refs/RotatingGK.pdf` (136), closed
by its (96) and (97), and it is reachable from a TOML file.** Two species in
closed form, `n` species by a safeguarded root find, normalised flux through the
existing bordered Newton, and `[source] Type = "rotating"` with
`examples/rotating-rectangle.toml` and `rotating-normalised.toml` as the worked
examples. `FLOW-PLAN.md` is the design; `CLAUDE.md`'s *Toroidal flow* has every
measurement, the three errors found in Li & Zhu, and the Maschke–Perrin reading
**this file previously got wrong** — its §4 is (136)'s isothermal closure at
every `γ`, and the paragraph that stood here called it an adiabatic one.

**What is not done here is item 11**, and the rotating output is what measured
it: `B` in the band between `Γ_h` and `Γ` is evaluated at the foot on `Γ_h`, and
the same closed form evaluated both ways differs by a factor of **1.7e5**.

## 10. The fixed-`q(ψ)` solver — meq, wanted, not designed

**`INVERSION-PLAN.md` IS NOW THE DESIGN FOR THE MACHINERY THIS NEEDS**, written
2026-09-02 after a literature survey (`refs/Refs.md`, *Solution inversion*). It is
staged IN-A and IN-0 to IN-6, and item 10 becomes reachable at **IN-2**, where the
flux-surface averages `⟨r^{-2}⟩_ψ` and `V′(ψ)` are measured against a closed form.
Two things in it change this item's shape: `v0-legacy`'s `FluxSurfaces` **has never
compiled**, so there is an algorithm to reuse and no working tool to extend; and the
global-structure work this item was assumed to need is **deferred**, because a
fixed-boundary problem with one axis has no interior saddles at all.

**Take `q(ψ)` as input and find `I(ψ)` from it**, rather than taking `I(ψ)`
directly as items 1 and 9 both do. The terminology is what people call it; the
requirement is ordinary. It is how a transport code hands an equilibrium code its
target, and it is what a coupling to MaNTA will want.

RoPP (142) is the relation:

```
q(ψ) = V′(ψ) I(ψ) ⟨r^{-2}⟩_ψ / 4π²
```

**Which makes this the item that needs the machinery item 9 was designed to
avoid.** `V′` and `⟨r^{-2}⟩_ψ` are flux-surface quantities, so a fixed-`q` solver
needs flux-surface averaging over level sets of `ψ_h` — contour extraction on an
HDG solution, at the solution's own order if it is not to throw away `k+1` — plus
an inner iteration for `I(ψ)`, plus that whole non-local dependence inside
`∂F/∂ψ` if Newton is to stay quadratic.

Three things worth writing down before anyone starts:

* **`FluxSurfaces` was a `v0-legacy` driver and was not ported.** Reach for tag
  `v0-legacy` before concluding it has to be written from nothing.
* **The non-local Jacobian has a precedent in the tree.** `ψ_ax` is one border
  row and column today and `HighBetaConvergence` is the acceptance criterion for
  it; a `q`-driven `I(ψ)` is the same shape with a continuum of rows rather than
  one, which is the part that is genuinely new.
* **It is independent of rotation.** A fixed-`q` static solver is useful on its
  own and is the smaller problem; doing it first and then composing is likely
  cheaper than doing it inside item 9.

Not on the critical path for anything, and it should not be started until item 9
has settled what a rotating source's profile set looks like.

## 11. `B` in the band — meq, DONE 2026-09-02

**FOUND WHILE WIRING THE ROTATING OUTPUT AND FIXED THE SAME WEEK.** Only
`samplePotentialWithFlux()` applied the Taylor step into the band between `Γ_h`
and `Γ`, so `B` — which went through `sampleComponent()` — was read at the foot
on `Γ_h`: about one node in ten of the interchange file piecewise constant,
`O(h)`, behind a mask saying `inside = 1`.

**Both halves are done.** `GridSampler::sampleComponentWithGradient()` continues
a vector field with its own `∇u`, read inside the element so nothing is
evaluated outside one, and `theBandVectorContinuesAtItsGradientsOrder` measures
**rate 2.20 against 0.92 for the foot, 124× smaller at `n = 16`**, on a
quadratic its space represents exactly so the band error is the truncation
alone. And the `.nc` now carries `byte extrapolated( Z, R )` beside `inside`, so
a reader can drop continued nodes rather than being told only how many there
were — 1667 on `miller-curved`, agreeing with the attribute.

**What is deliberately NOT done**, recorded so nobody re-derives it: this is
`O(h²)` at every `k` and does not reach `ψ`'s order, because `ψ` is continued
with a *solved* variable and `∇q` is a *differentiated* one. `div q = −F/r` and
`∂_r q_z − ∂_z q_r = −q_z/r` pin two of `∇q`'s four entries exactly, but leave
the symmetric traceless part differentiated — structure rather than an order,
and it would need the source plumbed into `GridSampler`. `CLAUDE.md`'s *Status*
section has the detail.

## Deliberately not yet

* **GPU and cuDSS.** Correctness-testable here; this card cannot say whether it
  is worth it. Consumer FP64 is 1/32 of FP32 where datacentre parts are 1/2, so a
  local timing can invert the production conclusion.
* **Partial assembly.** Not implemented for the HDG integrators at all, and its
  case is sum factorisation, which wants tensor-product elements where meq is
  triangles. A discretisation decision, not a flag.
* **The rest of the physics** — anisotropic pressure, NetCDF profiles, MaNTA
  coupling. `TODO` carries each with what has been established, and anisotropic
  pressure still has **no reference pinned**, which is the first thing it needs.
  **Sonic rotation has left this list**: it is item 9, and `FLOW-PLAN.md` is the
  design, which is why `TODO`'s *Sonic toroidal rotation* entry is now a stub
  pointing at it. `TODO` was cut with this file on 2026-09-01 — its PARDISO
  entry is down to the one open question, reproducibility under threading, and
  its performance entry to what the threading campaign did not already answer.

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
