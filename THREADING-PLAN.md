# THREADING-PLAN.md

**WHAT MEQ'S OWN SERIAL CODE COSTS, AND WHICH PARTS OF IT ARE WORTH
PARALLELISING.** The measurement this plans against is
**[M-126](MEASUREMENTS.md#m-126)** and the whole-run leg budget `meq --profile`
prints: on `examples/machine-f-diiid.toml` at `OMP = MKL = 8` with
`OMP_WAIT_POLICY=passive`, roughly 30% of the run threads and 70% does not. Of
the serial 70%, `NPCReduce`/`NPCRecover` (26.8%) and `DarcyForm::Reconstruct()`
(13.0%) are **MFEM's**, filed as
`../mfem-hdg-dev/doc/HDG-NPC-TRAVERSAL-FROM-MEQ.md`, and are out of scope here.

What is left is **about 25% of the run, and MEQ owns all of it**.

**THE HEADLINE IS THAT TWO OF THE SIX CANDIDATES ARE NOT THREADING PROBLEMS.**
The largest single item in this plan — `border dense solve`, and the
transmission slice of `constraint location` with it — is a **sparsity** problem:
the exterior transmission rows carry at most 984 nonzeros in a vector of about
109,000, and every dot product against them walks the zeros. Fixing that is
worth more than eight cores would be, and it is **exact** where threading a sum
is not. The second, `other (remainder)`, is not addressable by any loop MEQ can
write; the measured lever there is the build (M-77).

---

## STATUS: A, B, C AND D ARE BUILT. WHAT IS LEFT IS THE CLOCK

Landed as `0108f49`, merged as `a226cd3`. Each item below carries its own
marker; this table is what a reader needs before any of the estimates.

| item | state |
|---|---|
| **A** sparse border rows | **DONE.** `meq::CompressedRows` and `meq::firstNonzeroOutside` in `GradShafranov.cpp`, with the support check the item asks for grown teeth — it now asserts the recorded support covers every nonzero rather than only that the potential and trace blocks are clean. `tests/convergence/BorderAssembly.cpp` is new and `theCompressedBorderRowsContractExactly` requires `0.000e+00` |
| **B** the four plasma element loops | **DONE**, and it is five of them — `assemblePlasmaCurrent()` with the four the item names — through `MEQ_OMP()` in the new `src/meq/Threading.hpp`, per-element partial sums rather than a `reduction(+:)`, and a per-thread `std::exception_ptr` so a throw inside a region is not undefined behaviour |
| **C** `CriticalPointFinder::sweep()` | **DONE.** Per-element candidate buffers and a serial dedup over an element-ordered concatenation, with `theAxisSweepDoesNotDependOnTheThreadCount` in `tests/convergence/CriticalPointConvergence.cpp` |
| **D** the unnamed 3.2% | **DONE.** Two legs, `of which driver prepare` and `of which sweep overhead`, printed by `--profile` |
| **E**, **F** | not to be done. Unchanged |
| **G** re-assembly | still a restructure of the line search rather than a threading item, and still wants its own plan |

**THE CROSS-CUTTING ACCEPTANCE IS MET AND IT IS THE ONE WITH TEETH.**
`machine-f-diiid` solved by the pre-threading baseline at `OMP=1`, by the merged
tree at `OMP=1` and by the merged tree at `OMP=8` gives **one md5 over all
three** — every field, including all ten Gegenbauer modes. `OMP_NUM_THREADS`
does not change a printed digit, which is what every item's design was chosen
for.

**AND THE `find_package( OpenMP )` THE CROSS-CUTTING SECTION ASKS FOR IS IN
`CMakeLists.txt:433`**, with `OpenMP::OpenMP_CXX` linked `PUBLIC` onto
`meq_core`, so MEQ's own parallelism no longer arrives by inheritance from
MFEM's compile flags.

**§0's `--profile` RUN IS NOW TAKEN, ON A QUIET MACHINE, AND IT CLOSES ONE OF
THE TWO HALVES IT WAS FOR.** `I_p` is 0.5% at 5.14 cores, so item B's second
half is closed by measurement and the answer is that it was not worth doing.
The `transmission` half is **not** closed and could not have been by that run:
the leg it names is never written by any code path. §0 has the detail.

**WHAT IS STILL NOT DONE IS THE PAYOFF.** The A+B+C figure against the
pre-threading baseline — predicted `1.14×`, or `1.106×` counting only the
measured parts — needs a MEQ binary from before item A, which this tree does not
have lying about the way it had a pre-upgrade MFEM. **Until it is taken those
estimates stay estimates**; the standing rule is five quiet minutes before a
timed run, and a figure taken under another agent's build is not a measurement
about this work. That rule earned itself again on the day §0 was taken: the
first pair was measured while a peer's build ran, and the peer then reported
that the job they believed they had killed had been running throughout — they
had killed the process group of a PID `setsid` had reparented, so an empty group
died and the launcher reported success. Both arms were binned. **`pkill` by PID
is the same class of instrument error this file's own §0 is about.**

## 0. Read the sub-slices before doing anything

`meq --profile` already prints five sub-slices under `constraint location` —
`axis`, `cold full sweep`, `X-point`, `limiter`, `I_p`, `transmission` — and the
budget this plan was handed omits them. Two of the estimates below are
**unverified for want of exactly those lines**, and they are one run away:

```sh
OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OMP_WAIT_POLICY=passive \
  build/meq --profile examples/machine-f-diiid.toml
```

Where this plan says *unverified*, it means that line was not read. Read it
before sizing item A's second half or item B's second half.

**THE RUN IS TAKEN AND ONE OF THE TWO LINES CANNOT EVER SAY ANYTHING.**
`machine-f-diiid`, `OMP = MKL = 8`, `OMP_WAIT_POLICY=passive`:

* **`of which I_p` reads 0.016 s, 0.5% of the run, at 5.14 cores** — item B
  threaded it and there is no headroom left in it. **Item B's second half is
  closed by measurement**, and the answer is that it was not worth doing.
  **And the `cores` column is only readable because of `passive`**: M-126
  records the same leg reading **8.95 cores** under the default wait policy
  while no thread but one entered it, the CPU being the other seven spinning.
  Here 5.14 is real parallelism, which is item B having landed.
* **`of which transmission` reads `0.000 0.0% — 0`, and it reads that on every
  case, because NOTHING IN `src/` EVER WRITES `StepProfile::transmissionSeconds`.**
  It is declared, summed in `add()`, printed by `apps/meq.cpp` — and never
  incremented at any site. So §0's instruction to read that line before sizing
  item A's second half is unsatisfiable as written, and a reader who followed it
  would have concluded the transmission sweep is free.

**IT IS NOT FREE AND IT IS NOT MISSING EITHER — IT IS FOLDED INTO ITS PARENT.**
`exteriorTransmissionRows()` carries the `border assembly` `LegTimer` and is
called from the bordered solve, so the Γ sweep is inside **`border assembly`,
0.149 s and 2.6%** — which is the number item A's second half has to be sized
against, and is an upper bound on it rather than the thing itself. The sub-slice
that would separate them is the dead one.

**THE FIX IS ONE `LegTimer` AND A MOVED PRINT**, and it has to move: the leg is
printed under `constraint location`, where the work is not, while
`exteriorTransmissionRows()` lives under `border assembly`. Until that is done,
**item A's second half stays unverified — but for a reason the plan can now
name**, which is the difference between an unread line and a line that reports
nothing. `meq::Estimator` also calls `exteriorTransmissionResidual()` on the
adaptive path, which is a *different* parent, and that is the reason the slice
needs a decision rather than a timer.

**THE TRANSFERABLE PART IS THIS FILE'S OWN**: an instrument that reports zero
because it never counts is indistinguishable from an instrument reporting a
true zero, and `CLAUDE_INVERSION.md` §11.3 records the same shape — the tracer's
`fallbackLocations` read 0 before and after while the real count moved, because
the seed call handed it a local counter and discarded it.

---

## The two things the brief asked to check rather than assume

### `border assembly`'s 66 calls are NOT 66 sweeps of Γ

`GradShafranovSolver::exteriorTransmissionRows()`
(`src/meq/GradShafranov.cpp:3712`) sweeps Γ **once** and builds **all `N` mode
rows in that one sweep** — the mode loop is innermost, inside the quadrature
callback at `:3816`. It is called from one site, `:6882`, which sits in the
per-solve setup block above the Newton loop, so it runs **three times** in this
run, once per plasma-support sweep, not 66.

`exteriorConductorMoments()` (`:3648`) is the same shape and runs once per solve
as well. `assembleExteriorColumns()` (`:5592`) is `N` `LinearForm::Assemble()`
calls over Γ_h's boundary faces, also once per solve.

The 66 calls are dominated by **four element loops over the plasma**, each run
once per Newton step:

| function | line | what it builds |
|---|---|---|
| `assembleCurrentColumn()` | `:5211` | `dR/dL`, a vector on the potential block |
| `assembleCurrentNormalisationCorner()` | `:5273` | two scalars; falls back to four `assemblePlasmaCurrent()` calls when the source has no analytic derivatives |
| `assembleCurrentRow()` | `:5413` | `dG_I/dx`, a covector on the potential block |
| `assembleNormalisationColumn()` | `:5644` | `dR/ds`, the analytic `psi_ax` column |

At 12 Newton steps that is 48 calls, plus 3 + 3 for the two per-solve sweeps and
12 for the boundary-flux column. **So this leg is four plasma element loops, not
a boundary sweep**, and item B threads exactly those.

### `border dense solve` is a sparsity problem before it is a BLAS-3 problem

`rowDot()` (`src/meq/GradShafranov.cpp:7819`) contracts a border row against a
column, and the dense elimination at `:8115` calls it `(N+4)²` times plus
`(N+4)` more. On this case `n = blockOffsets[3] ≈ 109,200` (flux
`4848 × 6 × 2 = 58,176`, potential `4848 × 6 = 29,088`, trace
`≈ 7313 × 3 = 21,939`) and `nBorderTotal = 15` — `psi_ax`, `psi_bnd`, XP-3's two,
the current scale, and ten exterior modes.

Of those fifteen rows, **eleven walk the whole vector**: the current row (`:7889`)
and the ten exterior rows (`:7898`). The other four are already short — `psi_ax`
contracts over one element's trace dofs (`:7853`), `psi_bnd` calls
`limiterValue()` on one element, and XP-3's two contract over one element's flux
dofs (`:7862`). So the leg is `11 × 15 × 109,200 ≈ 18 M` multiply-adds reading
`≈ 290 MB` per call, twelve times. At 0.311 s over 12 calls that is **10 GB/s**,
which is one core streaming — the leg is memory-bandwidth-bound and the
arithmetic is incidental.

**And most of what it streams is zero.** `exteriorTransmissionRows()` writes only
`row( vdofs[ dof*d + j ] )` for `Elem1` of each Γ_h boundary face (`:3833`), and
the mesh has **82 boundary elements in total**. At `k = 2` a flux element carries
`6 × 2 = 12` vdofs, so an exterior row has **at most 984 nonzeros in 109,200** —
a factor of **111**. The function's own closing check (`:3852`) already asserts
that every entry outside the flux block is zero. The current row is nonzero only
on potential dofs of elements the plasma fill reached (`:5435`), perhaps `12×`
sparser.

**So a `dgemm` is the wrong instrument.** Compressing the rows takes this leg
from `18 M` terms to about `0.26 M`; what is left is small enough that BLAS-3
buys nothing measurable, and a blocked GEMM reassociates the sum where dropping
exact zeros does not. The standing preference for a maintained library
(`CLAUDE.md`, *Prefer a maintained library to a hand-rolled algorithm*) is not in
play: the question here is *which terms to sum*, not *how to sum them*.

---

## The items, ordered by payoff per risk

Shares are of the whole run, read off the budget above. Where a threaded leg's
efficiency is needed, **6.5 cores of 8** is used — the residual leg's own
observed figure on this case in this process, which is the only honest number
available for a loop of this shape.

**The arithmetic is validated against a published number.** M-126 sizes the NPC
traversal at 1.24× whole-run from 26.8% at 4.5×; the same method gives
`1/(1 − 0.268×(1 − 1/4.5)) = 1.26×`. Agreement to the second digit, so the
method below is not inventing its own units.

| # | item | share | after | saves | run speedup | risk |
|---|---|---|---|---|---|---|
| **A** | sparse border rows | 4.4% + ~1.8%* | ~0.2% | ~6.1% | **1.065×** | low |
| **B** | thread the four plasma element loops | 4.4% + ~0.9%* | ~0.8% | ~4.7% | **1.049×** | medium |
| **C** | thread `CriticalPointFinder::sweep()` | 1.7% | 0.3% | ~1.5% | **1.015×** | medium |
| **D** | instrument `outside solve()`'s unnamed 3.2% | 3.2% | — | unknown | — | none |
| **E** | the output writers | 3.8% | — | — | **do not do** | — |
| **F** | `other (remainder)` | 2.7% | — | — | **do not do** | — |
| **G** | `re-assembly` | 4.1% | — | — | **not a threading item** | — |

\* the starred halves are the `transmission` and `I_p` sub-slices of
`constraint location`, unverified — see §0.

**A + B + C together remove about 12.3% of the run, which is `1.14×`.** Taking
only the parts that are measured rather than estimated — 4.3% + 3.8% + 1.5% =
9.6% — it is **`1.106×`**. Both figures are on this case, at eight cores.

**Whether the set is worth starting**: it is worth about half of what the NPC
traversal alone is worth (1.24×), and about the same as `Reconstruct()` (1.10×),
both of which are somebody else's to write. It is three self-contained pieces of
work in one translation unit against two requests on another project's
schedule. The four compose additively — they are all in legs that currently read
1.00 cores — so if upstream lands both, the run loses about 43% and the combined
figure is near **1.7×**.

---

### A. Compress the border rows. First, and by a distance — **DONE**

**1. Where the code is.**

| | |
|---|---|
| rows built | `src/meq/GradShafranov.cpp:3713` `exteriorTransmissionRows()`, scatter at `:3833` |
| rows consumed, per constraint | `:6972` the `transmissionConstraint` lambda, dot at `:6978` |
| rows consumed, per elimination | `:7819` `rowDot()`, exterior branch `:7898`, current branch `:7889` |
| the elimination that calls it | `:8115`–`:8120` |
| the caller's own `n`-length copy | `:6895`, a `mfem::Vector row( n )` per mode |

**2. What the parallelism is.** None. This is not a parallel loop; it is a loop
over data that is 99.1% zero. The change is a representation:

* Every mode row is built in the **same sweep over the same boundary faces**, so
  all `N` rows share **one index set**. Build it once, as a sorted
  `std::vector<int>` of the Γ_h flux vdofs, and store the rows as a dense
  `modes × m` array with `m ≈ 984`.
* `transmissionConstraint()` becomes a gather of `state` onto the support
  followed by `m` multiply-adds per mode — one gather serves all `N`.
* `rowDot()` gathers the column onto the support once and contracts against each
  row. In the elimination, hoist the gather out of the `j` loop: each
  `columnZ( j )` is gathered once, not eleven times.
* The current row gets the same treatment against the plasma potential dofs,
  which `assembleCurrentRow()` already knows (`:5435`).

**3. The correctness hazard.**

* **Dropping an exactly-zero term from a sum of finite values is exact.**
  `total += row(i)*state(i)` with `row(i) == 0.0` leaves `total` bit-unchanged,
  so keeping the surviving terms in ascending index order reproduces the dense
  loop's every partial sum. The one exception is the sign of a total that is
  exactly zero: `-0.0 + 0.0` is `+0.0`. No downstream arithmetic here can see
  that, but it is the honest caveat and the acceptance below is written so that
  it would show if it mattered.
* **The zeros must genuinely be zero, not small.** They are: the scatter writes
  only the listed vdofs and the vectors are zeroed at `:3731`. A row that ever
  acquires a nonzero outside its recorded support is a silent wrong answer, so
  the existing `:3852` check must be kept and extended to assert the support
  list covers every nonzero entry rather than only that the potential and trace
  blocks are clean.
* **Device state.** `rowDot()` carries three `HostRead()` calls and a long
  comment saying why (`:7821`–`:7850`), including that `borderDofs` is device
  state too. A new index array is device state on exactly the same footing and
  needs the same `HostRead()`.
* **No threading hazard at all**, which is the point: `Mesh::GetElementTransformation`,
  `Mesh::GetBdrFaceTransformations` and `Mesh::FindPoints` are nowhere near this
  change.

**4. The payoff.** `border dense solve` 4.4% → about 0.06%, saving **4.3%**.
The `transmission` sub-slice of `constraint location` — 10 dots of length
`n` per constraint evaluation, at roughly 150 evaluations over the run —
estimated at half that leg, so a further **~1.8%**, *unverified*. Together
**~6.1%**, whole-run **1.065×**. Measured alone: **1.045×**.

**5. The cost and the risk.** One new small struct and three call sites, all in
`GradShafranov.cpp`, no new dependency and no new environment sensitivity. The
risk is that the support list and the scatter fall out of step, which is why the
check at `:3852` grows teeth rather than being deleted.

**6. How it would be tested.** A new case in `tests/convergence/XPointBorder.cpp`
or `FreeBoundaryCoupling.cpp`:

* `theCompressedBorderRowsContractExactly` — build the rows on the DIII-D
  fixture, contract every row against a handful of columns both ways, and
  **require `0.000e+00`**, not a tolerance. Exactness is the claim, so exactness
  is the assertion.
* The whole-solve pin: `psi_ax` to every printed digit (`3.759851e-01`), the
  **Newton iteration count**, and the constraint residual (`−6.466e-12`). The
  iteration count is the one with teeth — `CLAUDE_HDGGS.md` records that a
  Jacobian perturbed without moving the answer reaches the same discrete
  solution by a longer path, which no error norm sees.
* Checked to **discriminate**: perturb one support index deliberately and
  confirm the case goes red. A test that cannot fail is worse than no test.

---

### B. Thread the four plasma element loops — **DONE**, and it is five

**1. Where the code is.** `src/meq/GradShafranov.cpp` — `assembleCurrentColumn()`
`:5211`, `assembleCurrentNormalisationCorner()` `:5273` (analytic branch from
`:5362`), `assembleCurrentRow()` `:5413`, `assembleNormalisationColumn()` `:5644`,
and `assemblePlasmaCurrent()` `:5157`, which is the `I_p` sub-slice of
`constraint location` and is the same loop with a scalar output.

**2. What the parallelism is.** One `omp parallel for` over
`e = 0 .. mesh.GetNE()`, skipping `!elementInPlasma( e )`.

* **Reads**: `plasmaComponentMask.holds( e )` (a plain array read,
  `PlasmaComponent.cpp:243`), the element transformation, the potential shape
  functions, `state` on this element's potential dofs, and the source's
  `scaledF` / `scaledDFdPsi` / `normalisationDerivatives`.
* **Writes, three of the five**: `out( potentialStart + dofs[ j ] )`.
  **`potentialFes` is an `L2_FECollection`** (`:644`), so no two elements share
  a potential dof and **the scatter is disjoint by construction**. No colouring,
  no atomics, nothing — which is a stronger position than MFEM's own threaded
  element loop is in, since that one writes face dofs shared by two elements and
  needs a colour order.
* **Writes, the other two**: `assemblePlasmaCurrent()` and the corner accumulate
  **scalars**. A `reduction(+:)` reassociates. **Do not use one.** Write one
  partial sum per element into a `std::vector<double>` sized by `GetNE()`
  (39 KB here) and sum it serially in element order afterwards — **bit-identical
  to the serial loop**, at the cost of one array.

**3. The correctness hazard.**

* **`Mesh::GetElementTransformation( int )` hands out shared scratch.** All five
  sites already take the caller-allocated overload into a function-scope
  `thread_local` (`:5182`, `:5241`, `:5373`, `:5439`, `:5688`), so this hazard is already
  paid. It must stay that way: a future edit that reaches for the one-argument
  overload here is a silent wrong answer, and `CLAUDE.md` records six call sites
  that took it.
* **The function-scope scratch declared above the loop must move inside it.**
  `mfem::Array<int> dofs; mfem::Vector shape; mfem::Vector point;` are declared
  before the element loop at every one of the five sites — `:5165`, `:5225`,
  `:5357`, `:5423`, `:5670`. Shared across threads they are **exactly the
  `meq::SourceIntegrator` defect**: a concurrent `Vector::SetSize` is a
  reallocation, so the failure mode is memory corruption rather than a stale
  read. Declare them inside the loop body, or `firstprivate`.
* **`GridFunction` is not read here**, so `GetValue` / `GetVectorValue` and their
  internal `GetElementTransformation` call are not in play — these loops read
  `state` by raw index.
* **`Mesh::FindPoints` is not reachable** from any of the five. Checked: none of
  them calls `locateFieldPoint()`, `ContourTracer` or `CriticalPointFinder`.
* **The source must be reentrant, and it already has to be.**
  `meq::SourceIntegrator` calls `source->f()` and `source->dFdPsi()` from MFEM's
  threaded element loop (`:377`, `:498`), and `scaledF()` defaults to `f()`
  (`Source.hpp:569`). So this adds **no new obligation** on `meq::Source`.
* **Except for the one already recorded and not closed.**
  `meq::RotatingSource` throws from `f()` and `dFdPsi()` when a species
  temperature goes non-positive or the quasineutrality root find fails, and an
  exception escaping an OpenMP structured block is undefined behaviour
  (`CLAUDE_HDGGS.md`, *Threading, measured*). These five loops are MEQ's own, so
  unlike MFEM's the remedy is local: catch inside the region into a
  per-thread flag plus a stored `std::exception_ptr`, and rethrow after it.
  Do that rather than gating on the source type — a gate has to be kept correct
  as source types are added.
* **No new MKL nesting.** None of these loops calls BLAS, so `MKL_NUM_THREADS`
  is not newly exposed, and none of them is entered from inside MFEM's element
  loop, so there is no new nested region.
* **Bit-exactness.** With disjoint scatter and deterministic partial sums the
  result is bit-identical at any thread count, which is what the suite's pins
  need: `threadedAssemblyReproducesSerialAssemblyExactly` and its nonlinear
  sibling require `0.000e+00`, and that standard applies to MEQ's own loops too.

**4. The payoff.** `border assembly` 4.4% → 0.68% at 6.5 cores, saving **3.8%**.
The `I_p` sub-slice of `constraint location` a further **~0.9%**, *unverified*.
Together **~4.7%**, whole-run **1.049×**; measured alone **1.040×**.

Efficiency is likely a little under the residual leg's 6.5: these loops run over
the plasma subset rather than the whole mesh, and `elementInPlasma` makes the
trip count uneven, so `schedule(dynamic)` or `schedule(guided)` is the right
clause. That does **not** change the answer — the partial-sum array is indexed by
element, not by arrival.

**5. The cost and the risk.** Five pragmas and five scratch declarations moved.
The risk is the scratch, which is the defect this project has already met once,
and the exception path, which is recorded and open. Medium rather than low
because both are the kind that fail silently or not at all depending on the
mesh.

**6. How it would be tested.**
* `theBorderAssemblersAreThreadSafe` — run each of the five on the same state
  under `omp_set_num_threads( 1 )` and again at the ambient count, and
  **require every entry to agree at `0.000e+00`**. Not a tolerance: the design
  above claims exactness and the assertion should say so.
* Checked to discriminate, the way the existing nonlinear case was: hoist one
  scratch vector back out of the loop, confirm it fails, put it back. A race
  that only shows on some meshes is not a reason to trust a green run — it is a
  reason to run the discrimination.
* The whole-solve pin as in A: `psi_ax`, Newton iteration count, constraint
  residual, at both thread counts.

---

### C. Thread `CriticalPointFinder::sweep()` — **DONE**

**1. Where the code is.** `src/meq/CriticalPoints.cpp:932` `sweep()`; the unit of
work is `rootInElement()` at `:324`; the driver reaches it from
`apps/meq.cpp:3919` through `checkAxis()` → `findAxis()` (`CriticalPoints.cpp:885`)
→ `tryFindAxis()` (`:654`), whose seeded search falls through to a full sweep at
`:697` whenever the seeds do not return exactly one clean root.

**2. What the parallelism is.** A loop over every element, each running up to
three seeded 2×2 Newton solves against `q_h`.

* **Reads**: `meshRef` geometry, `fluxField` and `potentialField`.
* **Writes**: a shared candidate list `points`, two `mutable` counters
  `elementCount` and `newtonSolveCount`, and a shared
  `std::vector<mfem::IntegrationPoint> seeds`.
* **The writes are not disjoint and must not be made so naively.** Collect
  candidates into a per-thread vector, concatenate **in element order**, and run
  the **existing** deduplication serially over the concatenation. The counters
  become per-thread partial counts summed afterwards.

**3. The correctness hazard.**

* **The deduplication is order-dependent and the order is the answer.** The merge
  rule at `:987` keeps *"the one least outside its own element"* and the comment
  above it records the measurement that made it necessary: on Solov'ev at
  `k = 1, n = 6`, the candidate strictly inside its element is 2.7e-3 from the
  true axis and a neighbour 8.5e-2 outside is 6.1e-3, and *"element order is not
  a tie break"*. Threading the dedup makes which candidate survives depend on
  thread arrival, so **the reported axis position moves at the 1e-3 level with
  `OMP_NUM_THREADS`**. Keeping the dedup serial over an element-ordered
  concatenation is what makes the threaded sweep bit-identical, and it is not
  optional.
* **`Mesh::GetElementTransformation( int )`** — `rootInElement` already takes
  the caller-allocated overload into a `thread_local` (`:426`), with the comment
  saying exactly why. Already paid.
* **`GridFunction::GetVectorValue( int, ip, val )` and `GetValue( int, ip, vdim )`
  are reentrant here, and only because of what MEQ's spaces are.** Both take a
  `CalcShape` branch and touch no transformation when the element has
  `RangeType::SCALAR` and `MapType::VALUE`; the `else` branch of each calls
  `fes->GetElementTransformation( i )`, which is the mesh's shared scratch.
  MEQ's flux and potential spaces are `L2_FECollection`
  (`GradShafranov.cpp:642`, `:644`) — scalar range, value map — so the safe
  branch is the one taken. **This is a property of MEQ's spaces, not of MFEM**,
  and a space whose map type were ever anything else would make this a silent
  wrong answer with no other symptom.
* **`FiniteElementSpace::GetElementDofs` builds its table lazily.** With
  `elem_dof` present it is a table row read and reentrant; without it, it walks
  `mesh->GetElementEdges` and can populate shared tables. Make one serial call
  before the region — a single `GetElementDofs( 0, … )` — so the first parallel
  call finds the table built.
* **`Mesh::FindPoints` is not reachable from `sweep()`.** It is reachable from
  `tryFindCriticalPointFrom()` → `nearestElementCentre()` (`:717`), which is a
  plain element-centre scan and not `FindPoints`; the genuinely non-reentrant
  `Mesh::FindPoints` sites are in `FluxSurfaces.cpp` and `SurfaceAverage.cpp`
  and are not on this path. That said, `sweep()` is called from
  `findAxis()`, `tryFindAxis()` and `checkAxis()`, so **every one of those entry
  points inherits whatever is decided here** and none may be called from inside
  another parallel region.
* **`L2_TriangleElement::CalcShape` holds mutable scratch unless
  `MFEM_THREAD_SAFE`.** The installed MFEM sets it
  (`install/include/mfem/config/_config.hpp:78`), and `MFEM_USE_OPENMP` with it,
  so those members do not exist here. A build without it must keep the sweep
  serial, on the same build-conditional footing as `defaultAssemblyMode()`.

**4. The payoff.** `axis checks` 1.7% → 0.26%, saving **1.5%**, whole-run
**1.015×**. Small, and the reason to do it anyway is that the same function is
the cold-sweep fallback inside `tryFindAxis()` — M-100 measured a single sweep at
**100 ms against a seeded search's 2.5 ms**, and the seeded route falls through
to it whenever the seeds are unclean. Threading it caps the worst case as well
as the common one.

**5. The cost and the risk.** One loop, one per-thread candidate buffer, one
serial dedup pass. Medium risk, entirely in the dedup: get it wrong and the
answer depends on an environment variable, quietly, at a magnitude that looks
like mesh noise.

**6. How it would be tested.** `theAxisSweepDoesNotDependOnTheThreadCount` in
`tests/convergence/` — run `sweep()` at one thread and at the ambient count on
the same field, and require the returned vectors to be **equal element for
element, in order, in every field of `CriticalPoint`**, at `0.000e+00`. Then the
driver-level pin: `checkAxis()` reporting the same `psi`, `r` and `z` to every
printed digit at both thread counts. Discriminate by removing the element-order
concatenation and confirming it goes red on a mesh where the sweep finds
duplicates — the Solov'ev `k = 1, n = 6` case named at `:987` is one.

---

### D. `outside solve()`'s unnamed 3.2% — instrument it, do not thread it — **DONE**

`outside solve()` is 6.0%, of which `makeSolver` is 1.1%, `support move` 0.0%
and `axis checks` 1.7%. **The remaining 3.2% has no timer at all.**

Reading `apps/meq.cpp`, the largest identifiable thing in it is
`solver->prepare()` at `:3074` — the driver's own call, made once per adaptive
cycle before the support sweeps, which is a full `FormLinearSystem()`. The
`re-assembly` leg times the *same function* at 18 ms a call, so one call is about
0.26% of the run and the other 2.9% is unaccounted. `edgeFluxOf( state )` per
sweep, the per-sweep `printf` and `fflush`, and `plasmaComponentElements()` are
the other candidates.

**Nothing here should be threaded before it is named.** This project's own
record is unambiguous on that: M-126 exists because `otherSeconds()` was 45% of
a step and named nothing, and M-100 exists because the work was in `const`
members the call-site timers could not see. Add a `LegTimer` around the driver's
`prepare()` and one around the per-sweep block, re-run, and then decide.

**Cost**: two timers. **Payoff**: unknown, and that is the finding.

---

### E. The output writers, 3.8% — do not do

`output` less `postProcess` is 0.269 s. M-92 splits it on this same case: the
`.vtu` 0.130 s, the `.mesh` and three `.gf` 0.053 s, the NetCDF write 0.016 s,
and the whole `GridSampler` path — constructor and three sampling passes —
**0.010 s**. That is serialization and file I/O, not compute.

The one loop in it with real arithmetic is the sampling pass, and M-92 already
observes that grouping the located nodes by element makes it a `forall` with no
reduction and no contention, and records that it was **not done** because the
whole grid path is 27 ms. That judgement stands and this plan does not reopen
it: threading 10 ms buys 0.14% of a run.

Writing the three files concurrently would take 0.20 s to 0.13 s, about 1.0%.
MFEM's stream I/O is not documented reentrant and the three writers read the
same `GridFunction`s; 1% is not worth finding out. **Do not do.**

The 13.0% in this phase that *is* compute is `postProcess()` —
`DarcyForm::Reconstruct()`, four integrators re-assembled at the enriched order
per element. It is MFEM's, it is sized at about 1.10× whole-run in M-126, and it
is out of scope here.

---

### F. `other (remainder)`, 2.7% — not addressable by a loop MEQ can write

This is MFEM `Vector` arithmetic and allocator churn inside the Newton loop:
`unknown = savedState`, the `Add()` chain in each line-search trial
(`GradShafranov.cpp:8459`–`:8476`), `Norml2()` in `augmentedNorm`, and the
temporaries underneath them.

* **The `n`-length vectors are already hoisted.** `residual`, `column`, `y`, `z`,
  `scratch`, `columnL`, `zL` and `currentRow` are declared once per solve at
  `:5876`–`:5877`, outside the Newton loop. There is nothing to hoist.
* **`Vector::Add` and `Norml2` are MFEM's**, dispatched through `mfem::forall`.
  Threading them on the host means configuring an `mfem::Device`, which is a
  process-wide decision with its own hazards (M-79), and `Norml2` would then
  reassociate.
* **The measured MEQ-side lever is the build, not a loop.** M-77 measures the
  CUDA-enabled install at **1.07× the wall and +55.4 M allocations** against an
  otherwise identical `../mfem/install-nocuda`, with byte-identical output, and
  attributes 24% of the allocation traffic to `mfem::forall`'s host-lambda
  wrapper. **7% of the whole run for one CMake variable is larger than every
  item in this plan except A**, and it costs nothing to take for a production
  CPU build.

**Do not thread this leg. Take the no-CUDA install instead, and say in the
release notes that a CPU production build is a separate install.**

---

### G. `re-assembly`, 4.1% — a structural item, not a threading one

The leg is `prepare( bool )` (`:4179`) and the `reprepare` lambda that wraps it
(`:5865`). On the NPC path `rhsSource` is null — `usesNonlinearForms()` is true,
so the domain load is the non-linear form's and not a `LinearForm` — and what is
left is the Γ_h flux load over 82 boundary faces (`:4309`) and
`DarcyForm::FormLinearSystem()` (`:4377`). The second is essentially all of it,
it is MFEM's, it reads 1.05 cores, and MEQ cannot thread it.

**What MEQ can do is call it less often.** `reprepare()` runs **once per
line-search trial** (`:8455`, `:8633`) and once per fallback step, for one
reason: the exterior coefficients `a` change and they reach the residual as a
boundary load on Γ_h. But MEQ already knows — and `assembleExteriorColumns()`
at `:5592` depends on knowing — that **that load is exactly linear in `a`**, and
it already computes the per-mode load vectors `L_m`. So the trial's right-hand
side is `b₀ + Σ a_m L_m`, a linear combination of vectors MEQ has in hand,
rather than a re-assembly.

Two things stop this being a one-line change and both need their own
measurement: `FormLinearSystem()` also performs the essential-trace elimination,
and `reprepare()` rebuilds the `mfem::DarcyNPCOperator`. The elimination is
linear in the right-hand side, so in principle it composes; whether MFEM's entry
points expose that is a question for `darcyhybridization.cpp` and not for this
plan.

**Named here rather than planned**, because it is worth up to 4% and is the only
item in the budget where MEQ is paying a full assembly to change ten numbers. It
belongs in `BORDERED-GLOBALISATION-PLAN.md` or a successor, not here.

---

## Cross-cutting, and it applies to every item above

**`OMP_NUM_THREADS` MUST NOT CHANGE A PRINTED DIGIT.** Every acceptance in this
plan is `0.000e+00` rather than a tolerance, and each item's design is chosen to
make that achievable: exact zeros dropped from a sum (A), disjoint L2 scatter
plus per-element partial sums (B), a serial dedup over an element-ordered
concatenation (C). **An item that cannot be made exact should be brought back
for an argument rather than taken with a tolerance** — the suite already asserts
`psi_ax` to every printed digit in several places, and a threading change that
moves the last bit turns those into flaky tests, which is how a green suite stops
meaning anything.

**The existing bit-exactness pins are about assembly modes and hold at
`MKL_NUM_THREADS=1`.** `CLAUDE_HDGGS.md` records that at `MKL_NUM_THREADS=8` the
two assembly modes differ by 1.3e-15 in `ψ` because a blocked BLAS-3
reassociates. None of the loops in this plan calls BLAS, so none of them widens
that exposure — but item A's rejected GEMM route would have, and that is a second
reason to take the sparse route rather than the library one.

**No new nested regions.** All five of item B's loops and item C's sweep are
entered from the master thread of the Newton loop or from the driver, never from
inside MFEM's threaded element loop. If any of them ever moves inside one, the
inner region must be gated: MEQ would then have two teams where MFEM has one,
which is the case `EIGEN_DONT_PARALLELIZE` at `CMakeLists.txt:400` already
declines for the same reason.

**The build already carries OpenMP and no CMake change is required to start.**
`-fopenmp` is on `meq_core`'s compile line
(`build/CMakeFiles/meq_core.dir/flags.make`), inherited from MFEM's
`MFEM_CXXFLAGS`, and `libgomp` is on the link line through MFEM's own libraries.
A `#pragma omp parallel for` in `src/meq/` compiles and links today. **It would
still be right to add `find_package( OpenMP )` and
`target_link_libraries( meq_core PUBLIC OpenMP::OpenMP_CXX )`**, so that MEQ's
own parallelism does not depend on a flag arriving by inheritance from a
dependency that could be rebuilt without it — and so that a build against an
MFEM without OpenMP fails at configure rather than silently serialising. Guard
every new region on `defined( MFEM_USE_OPENMP ) && defined( MFEM_THREAD_SAFE )`,
as `SolverContract.cpp:12` already does, since item C's `CalcShape` argument
depends on the second.

**Take the profile under `OMP_WAIT_POLICY=passive`.** M-126 measures two fifths
of the `%CPU` on this case as barrier spin; adding MEQ-side regions adds MEQ-side
barriers, so a `cores` column read without it will credit the new regions with
work they are not doing.

---

## What to do first

**A, B, C and D are built** — see the status table at the top. What survives of
this ordering is the two measurements nobody has been able to take, and they are
items 1 and 6 below. The rest is the record of why the set was taken in this
order.

1. **§0** — one `--profile` run, to read the `transmission` and `I_p` sub-slices.
   Two of the numbers above depend on them and they cost one command.
2. **A** — the sparse border rows. Biggest, exact, no threading hazard at all,
   one translation unit. 1.045× measured, ~1.065× including the unverified half.
3. **B** — the four plasma element loops. 1.049×. Second because the scratch and
   the exception path are real hazards that this project has already been bitten
   by once each, and because A's payoff is larger and its risk is lower.
4. **C** — the critical-point sweep. 1.015×, and it caps the cold-sweep worst
   case as well as the diagnostic.
5. **D** — two timers on `outside solve()`, at any point. It is not work, it is
   an instrument, and 3.2% of the run currently has no name.
6. **THE PAYOFF ITSELF**, against the pre-threading baseline binary, as
   interleaved pairs at `OMP = MKL = 8` under `OMP_WAIT_POLICY=passive`. The
   prediction is `1.14×`, or `1.106×` counting only the measured parts, and
   until a quiet machine produces it this plan has delivered code and not a
   number.

**Not to be done**: E, the output writers (1% at best, for I/O concurrency
nobody has established is safe); F, `other` (threading it is not available to
MEQ, and the 7% that *is* available is a build variable); and the two Γ sweeps
in `border assembly`, `exteriorTransmissionRows()` and
`exteriorConductorMoments()`, which run three times each over 82 boundary faces
and are already below the noise. G is real and worth up to 4% but is a
restructure of the line search, not a threading item, and belongs in its own
plan.
