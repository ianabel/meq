# DEVICE-PLAN.md — a fully device-resident MEQ solve

**THE TARGET IS RESIDENCY, AND THAT IS A DECISION RATHER THAN AN AMBITION.**
The whole calculation runs on the device: the iterate, the local blocks, the
trace system and the border all live in device memory from `prepare()` to the
last Newton step, and the host issues kernels and reads scalars. **Host
threading is then irrelevant** — not preserved, not traded against, simply not
part of the design. Anything that reads "we would lose the host threading" is
an argument about a *partial* offload and does not bear on this plan.

**THE PARTIAL PATH WAS MEASURED AS A LOSS AND IS NOW MEASURED AT PARITY ON THE
SOLVE, WHICH CHANGES THE ARGUMENT FOR RESIDENCY WITHOUT WEAKENING IT.**
`--device cuda` runs `machine-f-diiid` correctly, its **solve leg overlapping
the host-threaded default's** — 4.53/4.75/4.75 s against 4.46/4.85/4.94 — and
the whole 1.07–1.13x wall loss is the **output** stage
→ **[M-223](MEASUREMENTS.md#m-223)**. [M-173](MEASUREMENTS.md#m-173)'s 2.14x
was a measurement of a configuration and said so: its device arm ran
`AssemblyMode::Batched`, whose host cost alone is 1.74x
([M-185](MEASUREMENTS.md#m-185)), and since then the driver's device default
became `Serial` ([M-215](MEASUREMENTS.md#m-215)), the NPC traversal moved to
kernels ([M-209](MEASUREMENTS.md#m-209)) and MEQ's integrand reached the device
([M-197](MEASUREMENTS.md#m-197)).

**WHAT IS LEFT OF THE LOSS IS ONE HOST LEG AND IT IS ITEM 7.** The output leg is
`Reconstruct()` running serially because a Device and a host OpenMP region that
touches MFEM memory are mutually exclusive — measured by a thread count rather
than argued, 0.90 and 0.97 s at `OMP = 8` and the same at `OMP = 1`, against a
host leg that halves. So the partial path no longer argues for residency by
being slow; it argues for it by **still having a host leg to lose**, and that is
a narrower and more honest claim.

**AND FIVE SIXTHS OF THAT LOSS IS NOT THE DEVICE, WHICH ONE MISSING CONTROL
SETTLED** → **[M-185](MEASUREMENTS.md#m-185)**. M-173's two arms vary the
device *and* the assembly mode together; the cell it did not run is
`--device cpu` with `AssemblyMode = "batched"`, which needs no GPU. Three
interleaved rounds put it at **1.74x** slower than the threaded default **with
no device in the process at all**, and the device on top of it at **1.13x**.
So the thing to fix first is not a kernel: it is that a host loop cannot be
threaded while MFEM's memory bookkeeping is live, which is
`general/mem_manager.cpp` and is upstream's. The device's own share of the
penalty is 13%.

**AND THE PRECISION FOLLOWS FROM THE TARGET, BUT NOT FIRST.** This machine's
consumer part runs FP64 at 1/32 to 1/64 of FP32. On a partial offload that
hardly matters, because the bottleneck is host serialisation; **on a resident
run it is the whole economics**, which is why `MFEM_PRECISION = single` is a
section of this plan and not a footnote. **It is also the LAST item**, because
a double build is what validates every step of the residency work — bit
agreement between `--device cuda` and `--device cpu` — and converting first
destroys that instrument. §3.3 and §7. It settles the solver question on the
way: a single-precision build cannot carry UMFPACK, and does not need to.

**THE YARDSTICK IS A COUNT OF CROSSINGS, NOT A SHARE OF A PROFILE, AND THAT IS
THE ONE THING TO TAKE FROM THIS PAGE IF NOTHING ELSE.** A leg's share answers
*what is expensive on a clean host run*, which is the wrong question: the point
of a device run is that the whole calculation stays on the device and the host
never blocks it. Measured on `machine-f-diiid` over 12 Newton steps, MEQ's step
makes **893.2 host reads of a device-written state per step, of which 0.2 are
actual transfers** — 3 transfers and 0.888 MB over the whole run
→ [M-212](MEASUREMENTS.md#m-212), [M-223](MEASUREMENTS.md#m-223). **The two
columns answer different questions and conflating them is the mistake this plan
has already made**, and it made a second one on top of it: the *"3.0 transfers a
step"* this paragraph carried is [M-198](MEASUREMENTS.md#m-198)'s **pre-fix**
figure, 36 over 12 steps, kept after [M-200](MEASUREMENTS.md#m-200) and
[M-202](MEASUREMENTS.md#m-202) took the run to 3 — and kept beside `0.888 MB`,
which is a **total**. *A rate and a total printed side by side read as one
row.* transfers are what MEQ pays on today's partial offload, reads are what
it would pay once a resident device writes between them. The leg carrying most
of the READS is 0.9% of the profile and transfers nothing at all. §3.3 has the
table, §1.3 the inventory, and §3.4a the decomposition — which says the reads
are overwhelmingly contractions of a fixed operand against the state, and that
exactly one is a search.

**AND THERE ARE NOW TWO STEPS THAT DESTROY THE VALIDATION INSTRUMENT, NOT ONE.**
Bit agreement between `--device cuda` and `--device cpu` is what root-caused
M-79, M-130 and M-187. Single precision destroys it, and so does cuDSS: 0.8.0
leaves `DETERMINISTIC_MODE` off, so repeats of the same solve differ at up to
3.6e-15, and the switch that would fix it is unsupported at the `N + 4`
right-hand sides MEQ's bordered step uses. **Both go last, and neither is to be
worked around** — §3.5.

**This plan is an audit of what residency costs, not a commitment to pay it.**
It answers five questions separately: what MEQ's own code would have to change,
what is still missing in the MFEM it links, what order the work has to happen
in, what precision it runs at, and what would tell us it is working. It draws
no throughput conclusion from a timing taken on this machine, per `CLAUDE.md`'s
standing rule about the RTX 2070 SUPER.

**Every claim about somebody else's code here was taken by extracting whole
functions from the installed tree**, never by reading a line range —
`CLAUDE.md` records a false claim produced exactly that way. §6 says what could
not be verified.

## 0. Where this starts from, and one thing that has moved

`../mfem/install` is MFEM 4.10.1 on `meq-integration`, `MFEM_USE_CUDA = YES`,
`CUDA_ARCH = sm_75`, `MFEM_USE_CUDSS = YES`, `MFEM_USE_DOUBLE`,
`MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE`, and **`MFEM_USE_MAGMA = NO`** — so
`BatchedLinAlg`'s only device backend in this build is `GPU_BLAS`.

**The installed library is current with the merge, checked and not assumed.**
The nine-branch containment loop reports all nine contained; every `.hpp` under
`fem/darcy/` in the install is byte-identical to `meq-integration`'s; and
`nm` finds `LinearResidualBatched`, `HDGFaceBuildTables`,
`ComputeElementsHBatched`, `FactorElementsBatched`, `AssembleNLFaceGradBatched`,
`HDGElementMassBatched`, `CanBatchLinearResidual` and `CanCacheCondensation` in
`libmfem.a`. **So upstream's device machinery is linked into MEQ today.**

**AND UPSTREAM'S OWN POSITION HAS MOVED FURTHER THAN MEQ'S ACCOUNT OF IT.**
`doc/HDG-DEVICE-OFFLOAD.md` on `gf-hdg-linearise-first` now opens by withdrawing
the headline MEQ quotes it for:

> **The chain pays end to end on the device above ~4k elements**, which the
> headline of this file denied for the whole project.

with `hdgdevice -o 2 -cudss` against the same problem with every device setting
off, one build, one session:

| | device | control | |
|---|---|---|---|
| n=32, 1024 el | 0.1866 | 0.1688 | 1.11x **slower** |
| n=64, 4096 el | 0.4198 | 0.5589 | **1.33x faster** |
| n=128, 16384 el | 1.3059 | 2.1938 | **1.68x faster** |

Six of that plan's nine items are done; two remain (parallel shared faces,
non-NPC problems). **The gate MEQ's refusal of `TraceSolver = "cudss"` rests on
is still written there, unchanged**, and this plan does not propose lifting it —
see §3.3. But the sentence beside it in `apps/meq.cpp` and in `ROADMAP.md`, that
"group 2 needs a partial-assembly rewrite and is not built", is no longer the
whole truth: what is not built is a partial-assembly rewrite, and what *is*
built is a different thing — batched element and face kernels reached through
`CanBatch*` predicates — which MEQ already reaches on two of three axes.

**STATUS.** The audit stands; what has moved under it is §1.2 and §4's
stages 1 and 4 → **[M-130](MEASUREMENTS.md#m-130)**, which walked a bordered
free-boundary case under `--device debug`, found **two live sites this sweep
did not list**, root-caused the `postProcess()` abort as upstream's, and left
`machine-f-diiid` solving under `--device cuda` in 2 Newton steps to the host's
every digit. §5.1's five corrections are all taken. **Nothing about the
performance argument has changed** and §3.3's refusals stand unaltered.

**Three MEQ-specific facts make that crossover table not MEQ's**, and they are
the core finding of this audit:

1. **`hdgdevice` is `Mesh::MakeCartesian2D(n, n, Element::QUADRILATERAL)`**
   (`miniapps/hdg/hdgdevice.cpp:123`). **MEQ is triangles everywhere** — the
   driver's `MakeCartesian2D` passes `mfem::Element::TRIANGLE`
   (`apps/meq.cpp:311-312`), `tools/mesh/halfdisc.py` runs
   `gmsh.model.mesh.generate(2)` with no recombination, and
   `machine-f-diiid.mesh` carries geometry code 2 on all 4848 elements.
2. **MEQ's problems sit at or below that crossover.** Element counts of the
   shipped meshes: `miller-curved` 450, `soloviev-nstx` 1536, `limited-tokamak`
   1601, `diverted-tokamak` 2642, `free-boundary-halfdisc` 3866,
   **`machine-f-diiid` 4848**, `machine-d-tcv-shaped` 5093,
   `machine-g-mastu` 6707. **Five of eight are below 4096** — `miller-curved`,
   `soloviev-nstx`, `limited-tokamak`, `diverted-tokamak` and
   `free-boundary-halfdisc` at 3866, which this list first counted on the wrong
   side — and the largest is 1.6x the crossover where upstream measured 1.33x
   on quads.
3. **MEQ's Newton is a BORDERED Newton, and the border is host arithmetic
   today.** M-126 measures it at 30% of a threaded step (§1.4, §3.1). This list
   used to say it *"would not go to a device"*, which is an observation about
   how it is written rather than a property of the arithmetic — under residency
   it is **item 6**, work MEQ owns, not a share to be conceded.

## 1. What MEQ itself would have to do

### 1.1 `meq::SourceIntegrator` — the one integrator MEQ owns

`src/meq/GradShafranov.hpp:258` declares
`class SourceIntegrator : public mfem::NonlinearFormIntegrator`, and
`src/meq/GradShafranov.cpp:399` and `:445` define `AssembleElementVector` and
`AssembleElementGrad`. It is installed on the potential-mass **non-linear**
form's domain at `src/meq/GradShafranov.cpp:3965-3987`.

**It declares no device entry point and inherits the base's untouched.**
`mfem::NonlinearFormIntegrator` declares `AssemblePA`, `AssembleGradPA`,
`AddMultPA`, `AddMultGradPA` and `AssembleGradDiagonalPA`
(`fem/nonlininteg.hpp:99-144`); `SourceIntegrator` overrides none of them.

**A device path for it is blocked by the language, not by effort.** The
quadrature loop as written calls, per point:

```
tr.SetIntPoint( &ip );       // ElementTransformation, host virtual
el.CalcShape( ip, shape );   // FiniteElement,          host virtual
tr.Transform( ip, point );   // ElementTransformation,  host virtual
tr.Weight()                  // ElementTransformation,  host virtual
source->f( R, z, psi )       // meq::Source,            host virtual
source->dFdPsi( R, z, psi )  // meq::Source,            host virtual
```

`grep -c MFEM_HOST_DEVICE` over the installed `fem/eltrans.hpp`,
`fem/coefficient.hpp`, `fem/fe/fe_base.hpp` and `fem/nonlininteg.hpp` returns
**0, 0, 0, 0**. None of those calls can appear in a `[=] MFEM_HOST_DEVICE`
lambda. Upstream says the same thing about its own integrators in two places,
and both are doc-comments on the batched routes that exist *because* of it
(`darcyhybridization.hpp:1797`, `bilininteg_hdg.hpp:1335`).

**So a device `F( R, z, psi )` needs the profile machinery re-expressed as
plain data plus a device function, and MEQ's profiles cannot be taken as they
are.** What `meq::Profile` actually is:

* `meq::Profile::operator()`, `prime()` and `doublePrime()` are **pure
  virtuals** (`src/meq/Profiles.hpp:63, 68, 90`), so every evaluation is a
  vtable dispatch. A device lambda cannot call through a host vtable.
* `meq::SplineProfile` holds a `std::vector<Knot> knotData` and finds its
  interval by search (`src/meq/Profiles.hpp:329, 327`). The knots are a flat
  POD array and **would port**: an `mfem::Vector` of knots plus a branchless
  or binary interval search is ordinary device code. This is the easy case.
* `meq::SplineProfile`'s sampling constructor takes two
  `std::function<double(double)>` (`Profiles.hpp:258`) — used only at
  construction, so it does not reach the hot loop.
* **`meq::RotatingSource` does a safeguarded Newton root find per point**
  (`RotatingSource.cpp:276 rootFindState`, `:403` the Newton, `:419` a
  `throw std::runtime_error` when it does not converge). A `throw` from a
  device lambda is not expressible; the root find itself would have to become
  a fixed-iteration scheme with an error flag returned in an output array, and
  the non-convergence case becomes a host-side check after the kernel.
* `meq::CoilSet::f( R, z )` is **cheap and would port trivially** — a
  point-in-rectangle test per coil and a sum of current densities
  (`src/meq/Coils.cpp`). The expensive elliptic-integral work (`coilPsi`,
  `ellipsePsi`, Boost.Math) is in the **initial guess and the exterior
  datum**, not in `f()`, so it is not on the per-step path at all. That is
  worth knowing before anyone plans to port Boost.Math.
* The plasma-component mask is already flat: `PlasmaComponent::holds( element )`
  over a `std::vector<char>`, consulted once per point
  (`GradShafranov.cpp:365`). It ports as an `mfem::Array<char>`.

**The realistic shape of a device `SourceIntegrator`** is therefore: a
`NodalSourceFunction`-style POD description of the source — an enum tag, a
handful of scalars, and flat knot arrays — evaluated by a
`MFEM_HOST_DEVICE` free function, with the virtual `meq::Source` hierarchy kept
on the host for everything that is not the quadrature loop. That is a second
implementation of the physics, and **CLAUDE.md's standing preference against
hand-rolled duplicates applies to MEQ's own code too**: two expressions of
`F( R, z, psi )` that must agree is exactly the maintenance hazard upstream
names when it refuses to reproduce `EvalStabilization` in a device lambda —
*"reproducing that formula in a device lambda means duplicating an integrator's
internals and diverging from them silently"*.

**AND THE INTERPOLATORY ROUTE DOES NOT DODGE THIS.** `gf-interp-hdg-dev` is
merged and linked, and `mfem::NodalReactionFunction::Eval` /
`EvalJacobian` (`fem/darcy/reaction_hdg.hpp:36, 44` in the install) are **host
virtuals taking `mfem::Vector`**. Interpolatory HDG changes *how many* times
`F` is evaluated, not *where*. MEQ calls none of `reaction_hdg.hpp` today —
`grep` over `src/`, `apps/`, `tests/` returns nothing.

### 1.2 Every host-pointer read and every alias, swept

**Swept systematically, not by symptom**, over `src/meq`, `apps` and `tests`,
for `MakeRef`, `Update()`-as-alias, `GetData()`, the memory accessors, the
`Sync*` family, `GetSubVector`/`SetSubVector`, `GetBlock`, `GetTrueVector` and
raw `double *`. The totals:

| pattern | production (`src/meq`, `apps`) | tests |
|---|---|---|
| `MakeRef(` | **5**, all `GradShafranov.cpp` | 0 |
| `Update()` forming a view | **8**, all `GradShafranov.cpp` | — |
| `GetData()` | **0** | 10 |
| `HostRead`/`HostWrite`/`HostReadWrite` | **17** | 1 |
| bare `.Read(` / `.Write(` / `.ReadWrite(` | **0** | 0 |
| `SyncMemory` / `SyncToBlocks` | **4** | 0 |
| `SyncFromBlocks` / `SyncAliasMemory` / `UseDevice` | **0** | 0 |
| `GetSubVector` / `SetSubVector` | 10 | 0 |
| `GetBlock(` | 9 | 0 |
| `GetTrueVector` | **0** | 0 |
| raw `double *` from an MFEM container | **0** | 4, all `Mesh::GetVertex` |

**MEQ holds no raw pointer into an MFEM container in production code and never
calls the non-`Host` accessors.** That is the good half and it means the whole
question reduces to aliases and to raw `operator()`.

**The five aliases:**

| site | what it aliases | synced |
|---|---|---|
| `GradShafranov.cpp:672` `darcyFlux.MakeRef( fluxFes, solution.GetBlock( 0 ), 0 )` | solution block 0 | **yes**, `:9329` |
| `:673` `potentialGf.MakeRef( ..., solution.GetBlock( 1 ), 0 )` | solution block 1 | **yes**, `:9330` |
| `:674` `traceGf.MakeRef( ..., solution.GetBlock( 2 ), 0 )` | solution block 2 | **yes**, `:9331` |
| `:4371` `traceX.MakeRef( solution, blockOffsets[ 2 ], ... )` | solution's trace range | **no** |
| `:4372` `traceB.MakeRef( rhs, blockOffsets[ 2 ], ... )` | rhs's trace range | **no** |

**AND THE FIRST THREE ARE IN THE CONSTRUCTOR, NOT IN `buildForms()`.**
`GradShafranovSolver::GradShafranovSolver` runs `:546–679`; `buildForms()`
starts at `:3884`. `CLAUDE.md`'s M-79 account and its index both attribute them
to `buildForms()`. It changes nothing about the defect or the fix — the aliases
outlive both functions either way — but it sends a reader to the wrong function
to check them, and `buildForms()` is the one that gets re-entered.

**Eight more views, formed by `Update()` rather than `MakeRef`, and they are
the same thing**: `solution.Update( blockOffsets )` `:667`,
`rhs.Update( blockOffsets )` `:668`, `fluxRhs->Update( ..., rhs.GetBlock( 0 ), 0 )`
`:4297`, `potentialRhs.Update( ..., rhs.GetBlock( 1 ), 0 )` `:4355`,
`darcySolution.Update( solution, darcy->GetOffsets() )` `:4375`,
`darcyRhs.Update( rhs, darcy->GetOffsets() )` `:4376`,
`form.Update( fluxFes, load, 0 )` `:5621`,
`recoveryScratch.Update( darcy->GetOffsets() )` `:5834`. **None is synced**, and
`:4297`, `:4355`, `:4375` and `:4376` are re-formed on every `prepare()` — which
on the Picard path is once per iteration and on the bordered path once per
`reprepare()`.

**The explicit syncs MEQ carries today are 18 calls in six places**, and they
are more than the four sites `CLAUDE.md`'s index records:

| file:line | what | frequency |
|---|---|---|
| `GradShafranov.cpp:4764` | `refreshPlasmaComponent`, `state.HostRead()` | per residual |
| `:6150–6152` | `refreshXPoint`, `HostRead()` + two `HostReadWrite()` | per Newton step |
| `:7759` | the bordered `flush()`, `borderOut[ i ]->HostRead()` per column | per Newton step |
| `:7849–7851`, `:7882–7883`, `:7891`, `:7900` | `rowDot()` — seven, **including `borderDofs.HostRead()` on an `mfem::Array<int>`** | per border row per step |
| `:9328–9331` | `solve()`'s `SyncToBlocks()` + three `SyncMemory()` | per solve |
| `CriticalPoints.cpp:150–151` | the finder's constructor, both fields | per axis search |
| `WarmStart.cpp:179–181` | after `FindPointsGSLIB::Interpolate` | per warm start |
| `tests/performance/NpcThreadScaling.cpp:248` | `copyOf()`'s `g.HostRead()` | per solve |

The `Array<int>` one is the transferable half, and the comment at `:7840–7847`
says why: MFEM's device-aware vector operations take dof arrays by
`Array::Read( use_dev )`, so an array handed to one comes back device-valid and
`Array::operator[]` is as raw as `Vector::operator()`.

#### THIRTEEN UNGUARDED RAW READS OF THE NEWTON ITERATE, AND THEY ARE ON THE DEFAULT BORDERED PATH

The sync discipline covers `refreshPlasmaComponent`, `refreshXPoint`,
`rowDot()` and `flush()`. Every *sibling* that raw-indexes the iterate through
`operator()` was checked, and these have no `HostRead()` anywhere in the
function:

| function / lambda | line of the read | when it runs |
|---|---|---|
| `locateAxisPoint()` | `:6386–6389` | **every Newton step of every bordered solve** — this is `AxisConstraint::LocatedAxis`, the default |
| `peakAt()`'s nodal-maximum fallback | `:6772–6779` | per step, when the located axis is not in force |
| `limiterValue()` | `:5975–5985` | per step under a limiter; the `rowDot()` call site is incidentally safe, the direct one at `:8680` is not |
| `locateLimiterContact()`'s `valueAt` | `:3020` | per step under `LimiterConstraint::LocatedContact` |
| `assemblePlasmaCurrent()` | `:5199` | per step when the current is an unknown |
| `assembleCurrentColumn()` | `:5258`, `:5268` | " |
| `assembleCurrentNormalisationCorner()` | `:5390` | " |
| `assembleCurrentRow()` | `:5456`, `:5466` | " |
| `assembleNormalisationColumn()` | `:5717`, `:5736` | per step, `ψ_ax` being an unknown of every bordered solve |
| `transmissionConstraint()` | `:6979–6980` | per exterior mode per step |
| `projectUpDown()` | `:1955–1962` | per accepted step under `[solver] UpDownSymmetry` |
| `assembleExteriorColumns()` | `:5639` | once per mesh, a raw read of a `LinearForm` immediately after `Assemble()` |
| `recoveryScratch.GetBlock( 1 )` at `:4476` | — | wherever the recovered potential is read without going through `potentialGf` |

**THE FIVE `assemble*` HELPERS ARE THE ANALYTIC BORDER COLUMNS, WHICH IS THE
DEFAULT ROUTE — SO THE DEFAULT PATH IS LESS GUARDED THAN THE FALLBACK.** The
differenced route goes through `fieldResidual` and `reprepare()`, both of which
sit behind a sync; the analytic route it is tried *before* does not.

**AND YET M-91 MEASURES THE DIII-D CASE SOLVING UNDER `--device debug`, WITH
THE X-POINT, THE LOCATED AXIS AND TEN EXTERIOR MODES LIVE.** That is not a
contradiction and the reason is the rule that makes this whole class of defect
hard: **`HostRead()` leaves BOTH copies marked valid.** Once any funnel syncs
the iterate, every later raw read of the *same buffer* is correct until
something writes it on the device again. So these thirteen sites are protected
today **by ordering rather than by a guard** — `refreshXPoint` or
`refreshPlasmaComponent` happens to run on the same vector first on the
configuration that was measured.

#### AND THE SWEEP ASKED THE WRONG QUESTION AT A FOURTEENTH SITE

**A `HostRead()` ON THE BASE CANNOT SEE DATA THAT IS STILL ON THE ALIAS**, and
that is the one shape this whole sweep is blind to by construction. Every
funnel below syncs *the object about to be read*, which is right when that
object is the one that was written — and useless when the write went through an
`mfem::LinearForm` or a `GridFunction` viewing the same buffer.

`assembleExteriorColumns()` is the site. One `mfem::Vector load` is reused
across the modes; `load = 0.0` is device-aware, `LinearForm::Update()` aliases
it, `Assemble()` fills the **alias**, and `load.HostRead()` — a funnel this plan
asked for and got — reads the **base**, which still holds the zeros. Mode 0
survives because a fresh `Vector` has `UseDevice()` false; modes 1 to 9 read
exactly zero, so **nine of the ten exterior columns were missing from the
bordered Jacobian** under every `mfem::Device`. Newton falls from quadratic to
linear, the run takes 4× the iterations, and **the answer does not move** —
the exterior coupling is still exact in the residual.
→ **[M-187](MEASUREMENTS.md#m-187)**. Repaired with
`form.SyncAliasMemory( load )`.

**TWO CONSEQUENCES FOR THIS PLAN.** The §1.2 sweep needs a second pass asking a
different question — *for each device-aware write, who else views that buffer*
— rather than *is each read funnelled*. And **`mfem::Device( "debug" )` does not
close this class**: it catches an access to the wrong *location*, and here every
flag was honest about a location and wrong only about which copy was *newer*.
A 1170-second `debug` run reproduced the defect and faulted zero times. The
instrument that found it was a stage-by-stage host-side checksum diffed between
`--device cpu` and `--device cuda`, and that belongs in this plan's toolbox
beside `debug`.

#### THEY ARE FUNNELLED NOW, AND ONE OF THE THIRTEEN WAS NEVER A DEFECT

**TWELVE SITES TOOK A FUNNEL AND THE THIRTEENTH ALREADY HAD ONE.** Each
function now syncs once at its head rather than at each read, which is the
right granularity: `assemblePlasmaCurrent()` and the three `assembleCurrent*`
helpers read the iterate **inside an OpenMP region**, so a per-read guard would
be a per-thread call on a shared `Memory`.

Three of them wanted more than a `HostRead()`, and the distinction is the part
worth keeping. **`HostRead()` is for a vector that is only read.**
`projectUpDown()` reads `unknown` through a copy *and* writes it raw, so it
takes `HostReadWrite()` — a `HostRead()` there would leave the device copy
marked valid and the next sync would overwrite every write. The same applies to
the column `assembleExteriorColumns()` fills after `columns.back() = 0.0`,
which is a device-aware assignment. And the recovered potential is a
**`BlockVector` block**, i.e. an alias with its own flags, which is M-130's
family rather than this one.

**`transmissionConstraint()` WAS ON THE LIST AND IS NOT A DEFECT.** It reaches
the iterate only through `CompressedRows::dot( row, v, scratch )`, whose
`gather()` has called `v.HostRead()` since it was written. The row was inferred
from the *address* of the read rather than from the path to it. **A sweep that
lists call sites will over-report; only following each one to its accessor
settles it**, which is the same distinction between a comment and a call site
that section 3.6 records one level up.

**What this does NOT do is make any of it exercised.** These are unconditional
and inert without a Device, and the configuration that would falsify the old
ordering argument — free boundary, limited, `setPlasmaSupportFrozen()`, under
`--device debug` — is still not in the tree as a fixture. The guard is now
there; the test for it is not.

**That is exactly the property `rowDot()`'s own comment argues against relying
on**, and it names the measurement: with the producer syncs alone,
`Device( "debug" )` still faulted on the first line of that function's body. An
ordering property is not a guard, it is a coincidence that currently holds, and
the configurations where it plausibly does not are nameable rather than
hypothetical — a frozen plasma support skips `refreshPlasmaComponent`
(`:6245`), and a limited (non-diverted) machine never calls `refreshXPoint`. A
free-boundary limited case with `setPlasmaSupportFrozen()` reaches
`locateAxisPoint()` with nothing having synced the iterate.

#### AND THE SWEEP MISSED A FOURTH, WHICH THE DEVICE THEN FOUND

**THE TABLE OF FIVE ALIASES ABOVE READS THE FLOW ONE WAY ONLY.** It asks
whether each alias is synced *from the base*, and the first three are —
`solve()`'s `SyncToBlocks()` plus three `SyncMemory()` calls, `:9329–9331`.
What no column asks is whether anything carries a write **through** an alias
back **up** to the base, and nothing did. `prepare()` seeds the iterate by
assigning through `potentialGf` and `traceGf`; a device write through an alias
marks the **alias** device-valid, calls `AliasProtect` on the base's host range,
and leaves the base's flags exactly as they were. `formSystem()` then builds
`darcySolution` and `traceX` from those flags and inherits the lie.

→ **[M-130](MEASUREMENTS.md#m-130)**. Three `SyncAliasMemory()` calls at the top
of `formSystem()`, inert without a Device because `Memory::SyncAlias` returns on
its first line when the base is not registered. **Under `debug` it is a named
fault; under CUDA nothing is protected, the border columns are differenced
against a state that is not the iterate, and the run reports *"the bordered
Jacobian is singular in ( psi_ax, psi_bnd, a )"*** — one fault, two faces, and
the CUDA face names a border rather than a memory.

**AND A SECOND, WHICH IS THE `Array<int>` HAZARD THIS SECTION NAMES AND DOES NOT
SWEEP FOR.** `solveWithNormalisation()` reads `GetEssentialTrueDofs()` through
`Array<int>::operator[]` straight after `FormLinearSystem()` has eliminated on
it, so the array is device state and the accessor neither syncs nor
invalidates. One `HostRead()`. **The sweep above counts `Array<int>` only where
a `HostRead()` already exists** — `rowDot()`'s `borderDofs` — so the pattern was
named in the prose and not looked for in the tree.

#### THREE CANDIDATE DEFECTS THE ORDERING ARGUMENT DOES NOT COVER

These are candidates found by reading; no device run was taken to confirm any
of them (§6). Each is a *write* through one alias followed by a read through
another, which no amount of earlier `HostRead()` fixes.

1. **`picardStep()`, `:1393–1394`.**
   `darcy->RecoverFEMSolution( traceX, darcyRhs, darcySolution )` writes the
   field blocks through `darcySolution`; the next line is `out = potentialGf`,
   a **different** alias registration onto the same memory, with no sync
   between. This is M-79 site 3's exact mechanism.
2. **`solve()`, `:9319–9331`.** The `if ( !fieldsAreState )` branch runs the
   same `RecoverFEMSolution`, and the `solution.SyncToBlocks()` that follows
   carries flags **downward** from `solution` — whose own flags the write
   through `darcySolution` never touched. The order constraint is the one
   upstream's `DarcyNPCSolver::Mult` fix records: the blocks have to be carried
   **up** into `darcySolution` and `darcySolution` up into `solution` before
   `SyncToBlocks()` can mean anything.
3. **`postProcess()`, `:9379`.** `darcy->Reconstruct( darcySolution, traceX, ... )`
   hands in two alias registrations that `SyncToBlocks()` and the three
   `SyncMemory()` calls do not name. **M-91 records a remaining abort at exactly
   this call**, attributed to a stack-backed `Vector` inside MFEM on the
   strength of a `0x7fff...` host pointer.

   **WALKED, AND BOTH HYPOTHESES ARE WRONG** →
   **[M-130](MEASUREMENTS.md#m-130)**. It is neither MEQ's aliases nor a stack
   `Vector`: `Memory<T>::Wrap( ptr, n, own = false )` sets
   `h_mt = MemoryManager::GetHostMemoryType()`, which under a debug device is
   `HOST_DEBUG`, while `Memory<T>::Delete()` computes
   `std_delete = !registered && ( h_mt == MemoryType::HOST )` — so a
   **non-owning, unregistered** wrap calls `MemoryManager::Delete_` on a buffer
   it does not own. `Vector b_ze( b_z.GetData() + e*nd_ut, nd_ut )` at
   `darcyhybridization.cpp:10164` with `e == 0` aliases `b_z`'s base pointer
   exactly and de-registers it on destruction; the next element's `b_z = 0.` at
   `:10113` aborts. Evidence: `h_mt = 3`, `mc = 0`, `bytes = 192`, and a
   `Write_` breakpoint hit count of exactly **2**. **Upstream's, `debug`-only,
   and filed** as `../mfem-hdg-dev/doc/HDG-RECONSTRUCT-TOTAL-FLUX-FROM-MEQ.md`
   with the seven same-pattern sites beside it. Sites 1 and 2 above are still
   unwalked.

**Sites 1 and 2 are on the Picard fallback, which is the path a hard problem
lands on.** `usesNonlinearForms()` (`:3876`) returns **false** for
`Globalisation::AndersonPicard` and `PicardOnly`, and
`solveByPicardThenNewton()` (`:4384`) sets `AndersonPicard` at `:4411` — so the
driver's reactive ladder, on an *observed* Newton failure, moves the solve onto
exactly the two sites that carry no sync. **The device hardening covers the
Newton path and not the fallback from it.**

**AND THE TEST SUITE'S OWN HABIT IS NOT UNIFORM**, which matters for stage 0 of
§4: `NpcThreadScaling.cpp:241–253`'s `copyOf()` calls `HostRead()` and carries
the comment explaining why, while `SolverContract.cpp:643`,
`FluxGridConvergence.cpp:896–897`, `PlasmaConnectivity.cpp:1611–1612`, `:1686–1687`
and `TraceSolverScaling.cpp:140–141` reach raw buffers through `GetData()` with
no sync. Inert while nothing in the suite configures a Device — and a device
regression case would be built out of exactly these comparison helpers.

### 1.3 The host-only components — the phase-C barrier inventory, with counts

The expensive distinction is per-step against once. MEQ's own
`GradShafranovSolver::StepProfile` already names the per-step ones as legs
(`src/meq/GradShafranov.hpp:3411-3461`), which is how M-126 could measure them
— **and the call counts beside the seconds are what this section is now for.**
Under residency the question is not what a leg costs but how often it takes a
state the device has just written and brings it to the host. Counts below are
`machine-f-diiid`, 4848 elements, `k = 2`, **12 Newton steps**, §3.3's run.

**Per Newton step — each is a host read of a state the device produced.** ψ is
29,088 doubles on this case and `q` is 58,176, so the pair is **698 KB** and it
is `O( elements )`. **A count of reads is not a count of transfers and this
table does not claim it is**: once any one of these has brought the state to the
host it stays valid until the device writes again, so the transfers are bounded
below by the number of device-write/host-read **alternations**, not by the
column on the right. What the counts do establish is the *shape* — the step
alternates rather than separating into a device phase and a host phase, because
the border interleaves 14 backsolves with the functionals it differences.

| component | leg | calls | /step | what it is, and what it would take |
|---|---|---|---|---|
| `meq::CriticalPointFinder` | `axisSeconds`, `xPointSeconds` | 19 + 16 | 2.9 | roots of the **solved** `q`; the constructor takes BOTH fields to the host at `CriticalPoints.cpp:150-151`, under a comment saying so in terms. A device Newton per candidate element plus a reduction over candidates — the arithmetic is small and the *search* is the work |
| the exterior transmission **contraction** | `constraintSeconds`, unsliced | 190 | **15.8** | `T_m = transmissionRows.dot( m, state )` per exterior mode, so it scales with the mode count — 10 here. **The largest single count in the step.** It is a sparse gather of the state on the row's support and then a dot: entirely a reduction, and the most device-shaped item on this list |
| the exterior transmission **row assembly** | `transmissionSeconds` | 3 | — | `∫_Γ ( q·ν ) C_m dΓ`, the order-40 face quadrature at `GradShafranov.cpp:4571`. **Once per support sweep, not per step** — a different leg from the row above and easy to confuse with it |
| the plasma-current row | `currentSeconds` | 28 | 2.3 | a quadrature of `F/R` over the support. **Follows item 1**: it is the same integrand |
| the limiter row | `limiterSeconds` | 0 here | — | a max over a marked surface, or a point evaluation. A reduction; the X-point path does not use it |
| the border's dense algebra | `borderAssemblySeconds`, `borderSolveSeconds` | 66 + 12 | 6.5 | `rowDot()` and an `(N+4)²` dense solve, with the seven `HostRead()`s of §1.2. The dense solve is 14×14 — it belongs wherever the rows are, and the rows are the crossings |
| `meq::PlasmaComponent`'s flood fill | `componentSeconds` | 0 here | — | a graph traversal over element adjacency (`PlasmaComponent.hpp:157`). **The one item on this list with no device shape at all**, and it runs between sweeps rather than between steps, which is what makes that survivable |
| re-assembly | `prepareSeconds` | 16 | 1.3 | §1.4 — the DtN datum reaches the system through the right-hand side |

**READ THE TRANSMISSION ROW AGAIN BEFORE THE CRITICAL-POINT FINDER.** The
finder is the one this plan has always named — it has the explicit comment, it
is the obvious barrier, and it is **2.9 crossings a step**. The transmission
rows are **15.8**, have no such comment, and are five times the traffic. *The
component that documents its own transfer is not the one doing most of it*, and
nothing but a call count separates them.

**AND THE ORDER-40 FACE RULE IS NOT THE COST HERE.** `transmissionQuadratureOrder`
is 40 — 21 Gauss points on a SEGMENT, so positive-weighted and nowhere near the
order-25 simplex ceiling `GradShafranov.hpp:2043` documents — and it is measured
at **0.035 s, 0.5%** of the run. It is high because the integrand is not
polynomial: `VertexConePath`'s cone roughens the foot map `ξ ↦ a( x( ξ ) )`, so
the rule converges algebraically rather than spectrally — 2.61e-04 at q8,
1.01e-04 at q12, 1.34e-05 at q20, 6.68e-08 at q40, 5.40e-10 at q80, which is the
cone-off floor (`FREE-BOUNDARY-PLAN.md` §7.7). **What matters for residency is
the 15.8 crossings, not the 21 points inside each one.**

**Once per solve, or once per adaptive cycle — these cost one transfer at the
end and are not on the critical path:**

`DarcyForm::Reconstruct()` for `ψ*` (`postProcess()`, `GradShafranov.cpp:9379`
— once, but **0.60 s of a 5.15 s DIII-D run**, M-126's `perf` row);
`meq::ResidualEstimator` (once per adaptive cycle); `meq::GridSampler` and the
three output writers; `meq::FluxSurfaces`, `meq::SurfaceAverage`,
`meq::SurfaceFit`, `meq::Zernike` (the inversion, once); `meq::FieldTransfer`
(once per cycle, and already synced).

`meq::SafetyFactorSolve` is the exception that is neither: it is an **outer**
Newton wrapping whole solves, so its host work is once per outer step and the
solve inside it is the thing being discussed.

### 1.4 The free-boundary bordered path, specifically

**The answer changes, and it changes against the device.** Three reasons, all
measured:

* **`N + 4` right-hand sides against one factorisation.** M-126: 14 border
  columns on `machine-f-diiid.toml`. `NPCReduce` + `NPCRecover` is
  `O( elements × columns )` where the integrator-bound loops are
  `O( elements )`, and it reads **0.212 s at one thread and 0.212 s at eight**
  — **30.6% of the threaded step, the largest single leg**. Verified by
  extracting each whole function: none of `NPCReduce`'s two overloads,
  `NPCRecover`'s two, or `DarcyNPCSolver::ArrayMult` contains a `#pragma omp`
  or an `mfem::forall`.
* **Four or more residual evaluations per step.** M-80's high-beta row: the
  residual leg is **48.2% serial and 23.1% threaded**, against one gradient and
  one factorisation, because the border differences `ψ_ax`. Every one of them
  re-enters `meq::SourceIntegrator`.
* **The exterior datum forces a re-assembly per accepted step.** `prepare()`
  runs `FormLinearSystem()` again because the DtN datum reaches the system
  through the **right-hand side**. `HDG-NPC-TRAVERSAL-FROM-MEQ.md` measures
  **16 calls across a three-sweep, twelve-step DIII-D run** at **0.293 s and
  1.05 cores of a possible 8**. On a device that is 16 host→device uploads of
  a freshly assembled right-hand side, where upstream's fixed-boundary fixtures
  pay one.

**AND A RESIDUAL EVALUATION MAY NOT SIT BETWEEN `NPCGradient()` AND THE SOLVES
IT ENTITLES** — `CLAUDE.md`'s *Traps* entry — which is why the bordered step
flushes its queue before differencing. Under a moving support that constraint
is about *which problem the local blocks were factored for*, and it does not
relax on a device; it becomes a constraint on when the device-resident local
factors may be considered live.

## 2. What is still missing in MFEM

**The three questions kept apart.** `CLAUDE.md` records getting exactly this
wrong in MEQ's favour, by inferring a constraint's *route* from which form its
integrators were added to. The three are: **does a kernel exist**, **does a
route reach it**, and **does MEQ take that route**.

### 2.1 The integrator inventory, re-taken from the installed headers

| integrator | class decl | device entry points DECLARED | route under hybridization |
|---|---|---|---|
| `mfem::VectorMassIntegrator` | `fem/bilininteg.hpp:3468` | `AssemblePA`, `AssembleMF`, `AssembleDiagonalPA`, `AssembleDiagonalMF`, `AddMultPA`, `AddMultMF`, `SupportsCeed` | **not reached** — `ComputeElementMatrix()`; separately `HDGElementMassBatched()` under `AssemblyMode::Batched` |
| `mfem::VectorDivergenceIntegrator` | `fem/bilininteg.hpp:3941` | `AssemblePA`, `AddMultPA`, `AddMultTransposePA` | **not reached**; separately `HDGElementDivBatched()` |
| `mfem::TransposeIntegrator( mfem::DGNormalTraceIntegrator )` | `:431` / `:4481` | wrapper forwards eight; **wrapped class declares none** | never assembled at all — markers, per `DarcyForm::Assemble()`'s own comments |
| `mfem::NormalTraceJumpIntegrator` | `fem/bilininteg.hpp:4866` | `AssembleEAInteriorFaces` **only** | interior faces only, no boundary-face EA, no `AddMultPA` |
| **`mfem::HDGDiffusionIntegrator`** | `fem/darcy/bilininteg_hdg.hpp:844` | **NONE** | its face work goes through `HDGFaceScatterBatched()` / `HDGDiffusionFaceMatricesBatched()` instead |
| **`mfem::HDGExtensionIntegrator`** | `fem/darcy/extension_hdg.hpp:509` | **NONE** | curved and free-boundary only; a `BilinearFormIntegrator` on the flux mass form, so assembled once |
| **`meq::SourceIntegrator`** | `src/meq/GradShafranov.hpp:258` | **NONE** | §1.1 |
| `mfem::DomainLFIntegrator` | `fem/lininteg.hpp:107` | `SupportsDevice() → true` + `AssembleDevice()` — a **different API** | once per `prepare()` |
| `mfem::VectorBoundaryFluxLFIntegrator` | `fem/lininteg.hpp:505` | **NONE**, and no `SupportsDevice` override | once per `prepare()` |

**Calling a PA entry point on the `TransposeIntegrator` wrapper would abort.**
The wrapper's eight declared entry points are pass-throughs to `bfi`, and
`DGNormalTraceIntegrator` overrides none, so they reach
`BilinearFormIntegrator`'s base default, which is `MFEM_ABORT` —
`fem/bilininteg.cpp:24-58`. Moot in practice, since those two are markers.

### 2.2 The predicates, and which route MEQ takes

`DarcyHybridization` carries **thirteen** `CanBatch*` / `CanCache*` member
predicates and `bilininteg_hdg.cpp` **nine** free ones, each gating a `*Batched`
routine that self-refuses on entry. The ones MEQ's configuration meets or
misses:

| predicate | MEQ | consequence |
|---|---|---|
| `CanBatchPotFaceAssembly` | **passes** | `AssemblyMode::Batched` is reachable; MEQ's `ConstantStabilization` clears `HDGFaceScatterCanBatch()`'s refusal of a state-dependent hook |
| `CanBatchTraceAssembly` | **passes** | `TraceAssemblyMode::Batched` is MEQ's default — the one axis measured to pay, 7–9%, M-99 |
| `CanBatchLocalFactor` | **passes** | but see the next row |
| `CanCacheCondensation` | **passes as shipped**, and `LocalFactorMode::Batched` is one of its seven refusals | M-101: that mode is a **trade**, not a free knob |
| `CanBatchLinearResidual` | **refuses** | `Finalize()` puts MEQ on `LocalOpType::PotNL` and the predicate refuses `PotNL` outright. Upstream pins this with a `REQUIRE_FALSE` on a fixture built to reproduce MEQ's routing |

**The first four rows are MEQ's own instruments answering, not this audit
re-deriving the predicates**: `batchedPotFaceAssemblyTaken()`,
`batchedTraceAssemblyTaken()`, `batchedLocalFactorTaken()`,
`batchedLocalSolveTaken()`, `fluxMassIsPrefactored()` and
`condensationCacheTaken()` are declared at `GradShafranov.hpp:782–799` and
`meq --profile` prints them; M-99 and M-101 are the readings. **And
`batchedLocalFactorTaken()` answers `CanBatchLocalFactor()` — *can*, not
*did*** — M-99 records it reading `yes` even with `LocalFactorMode = "serial"`,
so that one row is a statement about the predicate and not about the run.

**`Finalize()` chooses `PotNL` for MEQ, and the condition is exact**
(`darcyhybridization.cpp`, whole function extracted; enum at
`darcyhybridization.hpp:726`):

```cpp
const bool block_pot_row_only = m_nlfi && m_nlfi->GetBlockRowMask() == (1 << 1);
if (IsNonlinear() && !m_nlfi_u && (!m_nlfi || block_pot_row_only) && !c_nlfi)
   lop_type = LocalOpType::PotNL;
```

MEQ's flux mass is a linear `VectorMassIntegrator`, so `m_nlfi_u` is null; MEQ
builds no block-nonlinear form, so `m_nlfi` and `c_nlfi` are null; `Mnl_p`
carries `meq::SourceIntegrator`. **`PotNL`.** That is the row that both buys MEQ
the condensation cache and refuses it `LinearResidualBatched()`.

**MEQ's face constraint takes the LINEAR route, confirmed by extraction.**
`DarcyForm::EnableHybridization()` falls past `if (M_p)` — null on the Newton
path, `GetPotentialMassForm()` being called only inside the `else` of
`usesNonlinearForms()` at `GradShafranov.cpp:3991` — and reaches

```cpp
else if (Mnl_p && fc_mode == FaceConstraintMode::Frozen
         && FaceIntegratorsAreLinear(Mnl_p, nl_elsewhere))
```

which fires: `nl_elsewhere` is true through `Mnl_p->GetDNFI()->Size() > 0`
(`meq::SourceIntegrator`), both `HDGDiffusionIntegrator`s are
`BilinearFormIntegrator`s, and `fc_mode` is its default. So the constraint is
`c_bfi_p`, **assembled once per solve rather than once per step**.

**`FaceConstraintMode` is a fourth axis MEQ does not name anywhere.**
`darcyform.hpp:129`, default `Frozen` at `:148`. MEQ never calls
`SetFaceConstraintMode()`, so it stays `Frozen` — which is what keeps it on the
linear route above. Worth recording precisely because nothing in MEQ mentions
it: a future upstream change of that default would move MEQ's constraint onto
the per-step route silently.

### 2.3 `DarcyHybridization`'s own loops

| | `mfem::forall` | `#pragma omp` |
|---|---|---|
| `fem/darcy/darcyhybridization.cpp` | 11 (+2 `forall_switch`) | 13 |
| `fem/darcy/bilininteg_hdg.cpp` | 12 | 0 |
| every other `.cpp`/`.hpp` in `fem/darcy/` | 0 | 0 |

Every `MFEM_HOST_DEVICE` in `fem/darcy/`'s **sources** is a `forall` lambda
capture. Every one in the installed **headers** — exactly two,
`darcyhybridization.hpp:1797` and `bilininteg_hdg.hpp:1335` — is prose saying
what does *not* carry it.

**`NPCReduce` and `NPCRecover` are neither threaded nor offloaded**, both
overloads of each, and `DarcyNPCSolver::ArrayMult` composes them around the
trace solve. The single-vector overloads reach `MultInvBatched()` when
`CanBatchLocalSolve()` is true, and that routine contains no `forall` and no
`omp` of its own — its parallelism, if any, is `BatchedLinAlg`'s backend, which
is `GPU_BLAS` or nothing in this build. **This is the leg M-126 measures at
30.6%**, and `../mfem-hdg-dev/doc/HDG-NPC-TRAVERSAL-FROM-MEQ.md` is MEQ's
standing request for the **host OpenMP** version of it.

### 2.4 The simplex wall, and it is a null dereference rather than a refusal

Upstream records that `Mesh::GetFaceGeometricFactors()` segfaults on simplex
meshes and that `HDGFaceTables` deliberately avoids it for that reason
(`bilininteg_hdg.cpp:2523-2531`). **Traced to the line, in the installed
sources**: `FaceGeometricFactors`'s constructor (`mesh/mesh.cpp:15504-15556`,
extracted whole) has **no geometry check of any kind** and calls
`fespace->GetFaceRestriction( LEXICOGRAPHIC, type, SingleValued )`
unconditionally; `FiniteElementSpace::GetFaceRestriction()`
(`fem/fespace.cpp:1509`) sends a conforming nodal space to
`ConformingFaceRestriction`, whose constructor (`fem/restriction.cpp:610-651`)
does

```cpp
const TensorBasisElement* el =
   dynamic_cast<const TensorBasisElement*>(fes.GetTypicalFE());
const Array<int> &dof_map_ = el->GetDofMap();   // no null check
```

A simplex Lagrange element is not a `TensorBasisElement`, so the cast is null
and the next line dereferences it.

**What this means for MEQ.** Every device face-geometry route in MFEM is
tensor-only today, and MEQ is triangles on every mesh it ships. It is not a
blocker for the batched face *scatter* — `HDGFaceTables` was built to avoid
that path — but it is a hard wall in front of anything that wants precomputed
face geometry, and it removes the whole sum-factorisation argument for partial
assembly, which `ROADMAP.md` already records: *"its case is sum factorisation,
which wants tensor-product elements where MEQ is triangles."*

**It is a crash where a refusal belongs, and upstream says so in its own
notes.** It is not MEQ's to file (§5).

### 2.5 What upstream still lists as open, in its own words

Two of nine items: **parallel shared faces** (the face kernels refuse
`ParallelC()`) and **non-NPC problems** (the kernel writes `H_data` and the
reduced route reads an assembled sparse `H`). Neither is MEQ's — MEQ is serial
and NPC. And the one thing upstream names as blocking the last step is not on
that list: the stabilization weight loop in front of the face kernel is a
**host virtual taking an `ElementTransformation`**, which is §1.1's wall met
from upstream's side.

**AND THE CORRESPONDENCE HAS MOVED SINCE, IN FOUR WAYS THAT CHANGE THIS PAGE.**
Read the replies in the repository root rather than this summary if anything
turns on the detail.

| | |
|---|---|
| **Request 1 — a device route** | **DONE, and it was never missing.** Upstream's own audit of MEQ's form configuration puts MEQ at **YES on seven of the ten** `CanBatch*` predicates, with the two NPC-only refusals correct for MEQ and `CanBatchPotMass` correctly `no` because `M_p` is null. §3.6 |
| **Request 2 — `HDGExtensionIntegrator` on a device** | **ANSWERED: nothing is owed.** It is on a bilinear form, assembled once, and its blocks already reach `Af`/`Ae` through `AssembleFluxMassBdrMatricesBatched()`, which batches the scatter and never dispatches on the integrator. Upstream corrected two of their own doxygen claims that said no such integrator existed |
| **Request 3 — a batched local residual taking a caller's integrand** | **CLOSED, BUILT ON BOTH SIDES, AND RUNNING** → [M-197](MEASUREMENTS.md#m-197). Upstream asked for the seam rather than the loosening — *"it is blocked on an interface, and the interface is ours … send us the POD and the signature you want to write `F` and `F'` in and we will build the other side"* — MEQ sent `HDG-BATCHED-RESIDUAL-SEAM-FROM-MEQ.md`, and `NPCResidual()` / `NPCGradient()` now take a caller's `f( element, point, x, y, u )`. MEQ's half is `meq::DeviceNPCOperator`, whose integrands are **functor structs** rather than lambdas because nvcc refuses an extended lambda in three of the shapes this wanted. **Two refusals had to fall first and only one was upstream's**: `CanBatchNonlinearResidual()` requires an explicit `IntegrationRule` and `meq::SourceIntegrator` chose its own per element ([M-193](MEASUREMENTS.md#m-193)), and the seam is reached only by a caller that calls `NPCResidual()` itself ([M-194](MEASUREMENTS.md#m-194)) |
| **The 2.14× of [M-173](MEASUREMENTS.md#m-173)** | **diagnosed as upstream's and half fixed, and the OTHER half went away with the driver's device default.** `AssemblyMode` is one enum, so asking for `Batched` asked for *serial* everywhere a kernel did not reach; `ThreadHostLoops()` now decides that and every host loop gates on it, worth 1.82× to 4.13× on their own fixture. **MEQ gains none of THAT**, because the same function returns false under a Device — §3.7. What MEQ gained instead is `AssemblyMode::Serial` under a device ([M-215](MEASUREMENTS.md#m-215)), removing the 1.74× M-173's arm was carrying: **the solve leg is now at parity and only the output stage loses** → [M-223](MEASUREMENTS.md#m-223) |

**AND ONE ITEM UPSTREAM NAMES AS NOT STARTED IS ONE THIS PAGE USED TO CARE
MOST ABOUT AND HAS SINCE DOWNGRADED** → [M-216](MEASUREMENTS.md#m-216):
*"`general/mem_manager.cpp`'s thread safety is now the cheapest large
win on the device path … We have not started it."* MEQ has since shown it is
reachable from a **consumer's** thread team as well as from MFEM's own loops —
`../mfem-hdg-dev/doc/HDG-MEMORY-MANAGER-CALLER-THREADS-FROM-MEQ.md`,
[M-188](MEASUREMENTS.md#m-188) — which widens the justification without changing
the ownership.

## 3. The ordered list, with the blocking dependencies

### 3.0 WHERE THIS ACTUALLY IS, so the rest of the section can be read as history

**THE TABLE IN §3.3 IS THE LIVE ONE AND EVERYTHING AROUND IT IS THE REASONING
THAT PRODUCED IT.** Sections 3.1 to 3.7 were written before most of the work and
several of them argue for orderings the measurements have since overturned; they
stay, because *a hypothesis that was measured and falsified is the half of this
file worth keeping*. What a reader needs first is the state:

| item | |
|---|---|
| 1 · the element-local residual **and gradient** | **DONE**, both halves, through upstream's caller-integrand seam → [M-197](MEASUREMENTS.md#m-197) |
| 2 · `NPCReduce` / `NPCRecover` | **DONE, MERGED, AND MEQ TAKES IT** → [M-209](MEASUREMENTS.md#m-209), [M-212](MEASUREMENTS.md#m-212). `CanBatchNPCTraversal()` reads **yes** on MEQ's production bordered case under `--device cuda`, so both legs and both overloads run as kernels including the blocked `ArrayMult()` path; upstream measure **4.0× and 4.5×** on the 14-column calls at MEQ's size. The crossing count did not move — those crossings were already free. **And M-199 was wrong that MEQ already took half of it**: `CanBatchLocalSolve()` reads yes about the SINGLE-column overloads |
| 3 · the constraint functionals | reductions **DONE** → [M-190](MEASUREMENTS.md#m-190); `I_p` **DONE** as an element loop → [M-192](MEASUREMENTS.md#m-192); the limiter is a search and is **measured not worth porting** → [M-201](MEASUREMENTS.md#m-201) |
| 4 · the bordered border | **DONE AND WIRED** → [M-208](MEASUREMENTS.md#m-208), and wiring it measured the route as an **11× LOSS in transfers** on today's layout, because `flush()` brings the columns to the host. Gated on `Memory::DeviceIsValid()` of the columns, so it turns itself on under residency and is reported by `borderDotsRanOnDevice()` |
| 5 · the host branches | **needed no kernel** → [M-195](MEASUREMENTS.md#m-195) |
| 6 · cuDSS as the device default | after everything, and §3.5's reason for that has itself been corrected — see below |
| 7 · `mem_manager` thread safety | **OPTIONAL for MEQ, and RE-SIZED UPWARDS ON THIS CASE** → [M-216](MEASUREMENTS.md#m-216), [M-223](MEASUREMENTS.md#m-223). M-216 put it under 2% of a device run and expected it to shrink as the routes wire up; the opposite happened, because the routes wiring up is what removed everything else. With the solve leg at parity, the serial `Reconstruct()` this gate forces is **the whole of the 1.07–1.13× a device run loses**, about 10% of the wall. *A transition cost is a property of which legs have kernels yet, so it grows as a share every time another leg stops being the problem.* MEQ's own half is done → [M-188](MEASUREMENTS.md#m-188) |
| 8 · single precision | last — and **the whole of it has now been run, not only the early experiment** → [M-204](MEASUREMENTS.md#m-204), [M-222](MEASUREMENTS.md#m-222). The bordered Newton survives FP32 on both bordered machines; the single suite reads **62 of 63** against the double control's 62/63, the one red in each being the deliberate `PlasmaEdgeConvergence`; and the cost was one solver constant, four laws about where a floor comes from, and a per-case accuracy that ranges over 560× |

**AND THE STATE OF THE PARTIAL PATH ITSELF HAS CHANGED, WHICH THIS TABLE'S ITEMS
1, 2 AND 7 ADD UP TO** → **[M-223](MEASUREMENTS.md#m-223)**. On
`machine-f-diiid` under `--device cuda`, asked of the library rather than
inferred: **device integrand yes, batched residual and gradient yes, batched
local solve yes, whole NPC traversal yes, condensation cache yes, trace assembly
batched yes**, face assembly deliberately serial, bordered row dots off behind
their data gate. Three of upstream's four groups reach the device, the solve leg
is at parity with the host-threaded default, and `DeviceParity` is green. **What
residency still has to buy is phases A, B, D and E and the border's
alternation** — not the legs.

**THE NUMBER THIS SECTION IS ORDERED BY HAS BEEN CORRECTED TWICE AND IS NOW
INSTRUMENTED.** §3.3's 21.1 host reads a step was a proxy off leg call counts;
[M-196](MEASUREMENTS.md#m-196) instrumented it at 333 and was itself wrong
(17 of 36 funnels); [M-198](MEASUREMENTS.md#m-198) is the complete one at
**891.8 reads and 3.0 transfers a step** — the transfer figure being a rate over
12 steps, and superseded by [M-200](MEASUREMENTS.md#m-200)'s repair to 3 over
the whole run. It reads **893.2 reads and 0.2 transfers a step** on the
refreshed library with the NPC traversal on kernels
→ [M-212](MEASUREMENTS.md#m-212). **Read the transfer column, not the
read column**: a read of a host-valid buffer is free, and sorting by reads names
the sites to leave alone — MEQ's busiest funnel is the bordered `( N + 4 )²` row
dot at 240 reads a step and it transfers **zero**.

**AND THE TWO LARGEST THINGS MEQ HAS FIXED ON THIS PATH WERE IN NO SECTION OF
THIS PLAN.** A device-aware `= 0.0` on host-only scratch — the exterior-column
load and the X-point fields — was **12× MEQ's transfer count and 17× its
bytes** → [M-200](MEASUREMENTS.md#m-200), [M-202](MEASUREMENTS.md#m-202).
Neither appeared in a leg profile, because both are once per support sweep.
*The plan was a poor guide to where the cost was; the transfer counter was a
good one.*

**FOUR ITEMS OF §3.3 CHANGED ON BEING READ RATHER THAN ON BEING WORKED** —
`I_p`, the limiter, item 5 and item 2. Three were classified from what they
PRODUCE rather than from what their body does. **A plan item is only work once
somebody has read the body, and measured it on the case where it runs.**


### 3.1 The old Amdahl ceiling answered a different question, and it is now the measurement of the partial path

**THE ARITHMETIC BELOW IS CORRECT AND IT BOUNDS THE WRONG THING.** It asks *how
much of a step could move to a device while the rest stays on the host*, which
is the partial-offload question. Under residency there is no "rest": the legs
marked **no** are not a ceiling to be divided into, they are **work items**. The
table is kept because it is the best inventory of which legs are hard, and
because it predicted M-173 correctly.

*Arithmetic on M-126's published shares, not a new measurement.*
`examples/machine-f-diiid.toml`, 4848 elements, `k = 2`, PARDISO,
`AssemblyMode::Threaded`, `OMP = MKL = 8`, `OMP_WAIT_POLICY=passive`, step
total 0.692 s:

| leg | s | share | reachable today? | under residency |
|---|---|---|---|---|
| residual | 0.086 | 12.4% | only with §1.1 | **must port** — `meq::SourceIntegrator` |
| gradient | 0.045 | 6.5% | partly — `ComputeH` yes | `ComputeH` is already `BatchedLinAlg`; the integrators need §1.1 |
| trace factorisation | 0.058 | 8.4% | yes, cuDSS | **free** |
| trace backsolve | 0.081 | 11.7% | yes, cuDSS | **free** — the one leg M-173 measured a gain on |
| `NPCReduce` + `NPCRecover` | 0.212 | **30.6%** | yes in principle | dense, integrator-free; the 25 `forall` sites already cover much of it |
| constraint location | 0.023 | 3.3% | **no** — `CriticalPoints` | **must port**, and it is a mesh search (§3.4) |
| border assembly | 0.047 | 6.8% | **no** | **must port** — §1.4 |
| border dense solve | 0.039 | 5.6% | **no** | small and dense; a single batched solve |
| re-assembly | 0.003 | 0.4% | no | follows the residual |
| other (remainder) | 0.097 | 14.0% | **effectively no** | scalar reductions feeding host control flow — **the hardest 14%** |

**Device-able 69.8% under partial offload.** The three figures that used to
close this section — 3.31x on the step, about 2x on a run, realistically
1.1–1.3x — are all statements about the partial case and all superseded by the
measurement: **the partial case was 2.14x SLOWER**, not 1.1–1.3x faster, because
the arithmetic credited the host-bound legs with keeping their threading and
they do not. → **[M-173](MEASUREMENTS.md#m-173)**.

**Under residency the ceiling is the card, and this machine cannot measure it
in double.** An RTX 2070 SUPER at 1/32 FP64 is slower than eight host cores on
dense double work before any kernel is written; at FP32 it is not. **So the
residency work is built and validated in double and only then measured in
single** — §3.3 and §7. A resident double build on this part is expected to be
unimpressive and is not the deliverable; it is the correctness step that earns
the right to convert.

### 3.2 What is already on the device, and it is more than this plan used to say

**TAKEN FROM THE SOURCES, NOT REMEMBERED.** `fem/darcy/` carries **25 `forall`
launch sites** — 11 in `darcyhybridization.cpp` and 14 in `bilininteg_hdg.cpp`
— and they are exactly the batched routines and their gather/scatter helpers:
`AssembleFluxMassMatricesBatched`, `AssemblePotMassMatricesBatched`,
`AssembleDivMatricesBatched`, `HDGElementMassBatched`, `HDGElementDivBatched`,
`HDGDiffusionFaceMatricesBatched`, `HDGFaceScatterBatched`,
`HDGBdrFaceScatterBatched`, `HDGNLFaceResidualBatched`,
`HDGNLFaceGradScatterBatched`, `HDGMixedConductionResidualBatched`, and the
`HDGGatherFaceCols` / `HDGPackElementH` / `HDGScatterTraceRows` family.

**And the dense local algebra is `BatchedLinAlg`, whose GPU backend is linked.**
`nm` finds `GPUBlasBatchedLinAlg::LUFactor`, `::LUSolve`, `::AddMult` and
`::Invert` in `libmfem.a`, and `darcyhybridization.cpp` calls
`BatchedLinAlg::LUFactor` / `LUSolve` / `AddMult` / `Mult` at about twenty
sites covering `ComputeH` (`:4436`–`:4717`) and the batched local residual
(`:5845`–`:5875`).

**SO THE BATCHED ROUTE IS THE DEVICE ROUTE, AND M-99 MUST NOT BE READ AS A
VERDICT ON IT.** M-99 measures `AssemblyMode::Batched` as a **15 to 27 per cent
loss** and that measurement is on the **host**, where batching pays a gather
and buys nothing because the host was already threaded. **A share is a property
of a problem, not of a loop** — `CLAUDE.md` records that lesson from the NPC
traversal, and this is the same shape one level over: *a cost is a property of
a target, not of a route*. On a device, batching is not a knob that lost a
benchmark; it is the only way the work executes at all.

### 3.3 The order, and it is by RESIDENCY PHASE rather than by leg share

**THIS SECTION ORDERED THE WORK BY EACH LEG'S SHARE OF A HOST PROFILE, AND THAT
IS THE WRONG YARDSTICK.** A share answers *which leg is expensive here*;
residency asks *what crosses the bus, and how often*. The two give different
orders, and only the second decides whether a device run is worth having: a leg
at 0.9% of a host profile that takes the whole state to the host to produce one
scalar is a **barrier**, and a leg at 20% that stays put is not. The shape every
item below is ordered against:

> setup, reading configs and so on — host · **configure the device problem** ·
> run the problem on the device, receiving only diagnostic data, the host never
> blocking device computation · the compute ends · **and only then** the full
> state returns to the host, to be written out.

| phase | where | what crosses the bus |
|---|---|---|
| **A** setup | host, and it may take as long as it likes | nothing |
| **B** configure the device problem | **one upload** | mesh, dof maps, the source POD, the coil table, the initial iterate |
| **C** the solve | **device only** | diagnostics — small, and the host never blocks on them |
| **D** the answer | **one download** | ψ, `q` and ψ̂, every coefficient, once |
| **E** output | host | nothing |

**SO WHAT ORDERS THIS PLAN IS A COUNT OF PHASE-C CROSSINGS, AND IT IS TAKEN
RATHER THAN ESTIMATED.** `machine-f-diiid`, 4848 elements, `k = 2`, all three
batched keys `"serial"`, `--device cpu`, 6.513 s of wall at 1.32 cores — the
host arm of [M-187](MEASUREMENTS.md#m-187)'s verification — over three support
sweeps and **12 Newton steps**, read off the leg profile's own call counts. **The
counts are a property of the algorithm and not of the assembly mode**; only the
percentage below is this run's:

| what the host does to a state the device produced | calls | per Newton step |
|---|---|---|
| **scalar constraint functionals** — `peakAt`, `refreshXPoint`, the limiter, `I_p`, and the exterior `transmissionConstraint` | **253** | **21.1** |
| trace backsolves — the 14 border columns | 168 | 14.0 |
| border row assembly | 66 | 5.5 |
| residual evaluations | 28 | 2.3 |
| right-hand-side re-assembly for the DtN datum | 16 | 1.3 |
| gradient · `NPCReduce` + `NPCRecover` · trace factorisation | 12 each | 1.0 each |

**Twenty-one host-side scalar functionals per Newton step, every one of them
reading a state the device produced.** The payload is exact rather than
estimated: ψ is **29,088** doubles on this case and `q` is **58,176**, so the
pair is **698 KB**, and it is `O( elements )` — the fixed-boundary DIII-D case
at 102,840 elements is **21×** that.

**A COUNT OF READS IS NOT A COUNT OF TRANSFERS, AND THIS TABLE DOES NOT CLAIM IT
IS.** Once any one of these reads has brought the state to the host it stays
valid until the device writes again, so the transfers are bounded below by the
number of **alternations** rather than by the counts above. What the counts
establish is the *shape*, and the shape is the finding: the step does not
separate into a device phase and a host phase, because the border interleaves
its 14 backsolves with the functionals it differences. **That is why ordering by
share was wrong** — `constraint location` is **0.9%** of this profile, 3.3% of
M-126's threaded one (§3.4), and is the largest single obstacle on the page
either way.

**AND THE BORDER IS WHAT MULTIPLIES THEM, WHICH RECLASSIFIES IT FROM SMALL TO
STRUCTURAL.** A bordered Newton differences its constraint functionals at
`N + 4` perturbed states per step, so a step is not *solve, then check* — it
**alternates**, device solve against host functional, fourteen times over. The
border's own arithmetic is 1.9% of the profile; its *shape* sets the crossing
count for everything else, and that is a fact about the algorithm rather than
about how the code happens to be written. **It cannot be ported late.**

**THE ORDER, AND WHO OWNS EACH ITEM.** Still in double, and the precision sweep
still last, for the reason that has not changed: `--device cuda` reproducing
`--device cpu` **to every printed digit** is the instrument that root-caused
[M-79](MEASUREMENTS.md#m-79), [M-130](MEASUREMENTS.md#m-130) and
[M-187](MEASUREMENTS.md#m-187), and converting precision first leaves every
later disagreement with two candidate causes and no way to separate them.

| # | item | phase | whose | blocked by |
|---|---|---|---|---|
| — | ~~**MEQ's aliases and raw-read sites (§1.2)**~~ | — | MEQ | **DONE**, and §1.2's closing note says the sweep asked the wrong question at a fourteenth site — [M-187](MEASUREMENTS.md#m-187) |
| — | ~~**A device route through `DarcyHybridization`**~~ | — | — | **DONE, AND NEVER MISSING. `AssemblyMode = "batched"` IS THE ROUTE** — §3.6 |
| — | ~~**`mfem::HDGExtensionIntegrator` on a device**~~ | — | — | **ANSWERED: NOTHING IS OWED** → [M-185](MEASUREMENTS.md#m-185) |
| **1** | ~~**The element-local nonlinear residual**~~ | **C** | MEQ + upstream | **DONE, BOTH HALVES** → [M-197](MEASUREMENTS.md#m-197). `meq::makeDeviceNPCOperator()` routes MEQ's POD integrand through upstream's batched residual AND gradient; `--device cuda` reports `residual yes, gradient yes, device integrand yes` and reproduces the host to every printed digit. The coil term turned out POD-able — a closed rectangle test and a constant — so the production configuration is covered |
| **2** | **`NPCReduce` / `NPCRecover` — the FACE loop** | **C** | upstream | **HALF DONE, AND MEQ TAKES IT** → [M-199](MEASUREMENTS.md#m-199). The element-local solves are already batched through `MultInvBatched()`, gated on `CanBatchLocalSolve()`, which reads **yes** under `--device cuda`. What remains is the FACE traversal — host dense work through `GetData()` — which is what forces `du_all.HostRead()`/`dp_all.HostRead()` per call. The "unowned, 1.00 cores, not started" this row carried was wrong on all three |
| **3** | ~~**The constraint functionals, device-side**~~ | **C** | MEQ | **THE REDUCTIONS ARE DONE** → [M-190](MEASUREMENTS.md#m-190), and `I_p` is built as an element loop → [M-192](MEASUREMENTS.md#m-192). The list had **two** real reductions and one kernel serves both. The limiter is a SEARCH and is **measured not worth porting** → [M-201](MEASUREMENTS.md#m-201): one shipped case, 0.7%, zero transfers |
| **4** | ~~**The bordered border** — `rowDot()` and the `(N+4)²` fill~~ | **C** | MEQ | **DONE AND WIRED** → [M-190](MEASUREMENTS.md#m-190), [M-208](MEASUREMENTS.md#m-208). `dotAll()` is called per column for both compressed row sets, the host gathers are skipped when it runs, and `psi_ax`, the constraint and the X-point are bit-identical between `cpu` and `cuda`. **Wiring it MEASURED the route as a loss** — 3.0 transfers a step becomes 17.7 with the dots in front of the gathers and **33.1** with them instead, because `flush()` brings every column to the host and `dotAll()` pushes each 698 KB one back. So the gate is `Memory::DeviceIsValid()` **on the columns**, a predicate about the data rather than a key about the build: false today, true the moment the columns stop coming down |
| **5** | ~~**The host branches**~~ | **C** | MEQ | **NEEDED NO KERNEL** → [M-195](MEASUREMENTS.md#m-195). `Norml2()` is already an `mfem::reduce`, so the convergence test and the line search are library primitives; the support decision is FROZEN within a solve by XP-2 and reads **0 calls**. What remained was §4 stage 4b's crossing instrument, now built → [M-196](MEASUREMENTS.md#m-196), [M-198](MEASUREMENTS.md#m-198) |
| **6** | **`TraceSolver = "cudss"` as the DEFAULT under a device** | **C** | MEQ | items 1–5, and **strictly after all of them**. §3.5 — the refusal inverts once the rest is resident, and cuDSS is nondeterministic run to run at 3.6e-15 with no usable switch at `N + 4` columns, so this step **destroys the digit-for-digit instrument** and must not precede anything validated with it |
| **7** | **`general/mem_manager.cpp` made thread-safe** | B, C | **upstream** | **OPTIONAL for MEQ** → [M-216](MEASUREMENTS.md#m-216). *A host loop worth threading under a device is a loop that should be a KERNEL*, so this is a TRANSITION cost rather than a requirement — worth under 2% of a device run and falling. MEQ's own gate and the threaded-assembly refusal stay either way, because they are what make a device run CORRECT rather than fast. The 1.74× of [M-185](MEASUREMENTS.md#m-185) is a different measurement and is upstream's own |
| **8** | **`MFEM_PRECISION = single`, MEQ swept `double` → `real_t`** | all | **MEQ** | **items 1–7, deliberately** — §7. The **second** instrument-destroying step, and it pairs with item 6 |
| — | **The simplex face-restriction crash** | — | upstream | blocks nothing here; it gates tensor-only device *geometry*, which residency does not route through (§2.4) |

**ITEM 7 IS NOT ONLY A SPEED ITEM, AND THE SECOND HALF IS MEASURED**
→ **[M-188](MEASUREMENTS.md#m-188)**. **Read that entry's refreshed-install
table before acting on anything below**: upstream's `ThreadHostLoops()` has
landed and is now in the library MEQ links, so the abort this paragraph
describes is history and what remains is MEQ's own half.
`--device cuda` with `AssemblyMode = "threaded"` **aborts** in the memory
manager — *"alias pointer is not registered"*, *"host pointer is not
registered"*, `CheckHostMemoryType_`. Discriminated with a thread count rather
than an argument, same binary and same config: **`OMP_NUM_THREADS=1` exits 0 and
reproduces the host's own 7 / 3 / 2 Newton counts; `OMP_NUM_THREADS=8` aborts.**
`--device debug` at 8 threads aborts identically **with no GPU arithmetic
anywhere**, which rules the CUDA driver out, and `general/mem_manager.cpp`
carries `MemoryMap` and `AliasMap` as bare `std::unordered_map`s with no mutex,
no critical section and no atomic. Under a Device the lazy registrations in
`MakeAlias`, `SetDeviceMemoryType` and the device arms of `Read`/`Write`/
`ReadWrite` fire from inside `DarcyForm`'s threaded element loops; without a
Device those arms are never taken, which is why host-threaded assembly has
always been safe.

**AND IT IS NOT ONLY MFEM'S ASSEMBLY: MEQ'S OWN ELEMENT LOOPS REACH THE SAME
MAPS.** With every batched key `"serial"` — nothing of MFEM's threaded at all —
`--device debug` at 8 threads **segfaults** in `MemoryManager::Insert()` beneath
`meq::CriticalPointFinder::sweep()`'s own `GOMP_parallel` frame, because each
thread's `GridFunction::GetVectorValue()` registers an `Array<int>` of vdofs.
`--device cuda` survives that same loop only because `GetHostMemoryType()` is
`HOST` there and leaves owners unregistered, which is a fact about how often the
map is touched and not about safety.

**BOTH HALVES ARE NOW CLOSED, AND THE REFRESHED INSTALL SEPARATES WHOSE EACH
WAS.** MEQ gates all of its own OpenMP through one macro in
`src/meq/Threading.hpp`, which now carries OpenMP's `if` clause on
`meq::hostLoopsAreSafe()` — `!Device::Allows( Backend::DEVICE_MASK )`, copied
from `DarcyHybridization::ThreadHostLoops()`. That repairs the serial-assembly
cell outright and leaves the host path bit-identical, **and it is still
necessary against the refreshed library**: mutated out, the new install
segfaults exactly as the old one did, because no gate of upstream's can see a
caller's OpenMP region. The threaded-assembly cell is upstream's and is now
fixed upstream. MEQ's driver additionally refuses that pair above one thread;
**that refusal is no longer load bearing** and M-188 says on what narrow ground
it is kept.

**AND THE INSTALL HAS BEEN REFRESHED, WHICH IS WHAT SETTLED IT.**
`../mfem/install` is rebuilt from `meq-integration` at `6d05189467` and matches
the branch file for file; `CLAUDE.md`'s install-currency one-liner is the check
and it now says *install is current*. **Every device measurement in this
document taken before that refresh was taken against a library predating
`ThreadHostLoops()`** — the sweep counts and the M-187 fix reproduce unchanged
on the new one, checked, but any timing is not comparable across the boundary.

**AND IT COSTS RESIDENCY NOTHING, WHICH IS WHY IT IS A REFUSAL RATHER THAN A
BLOCKER.** MFEM's own `ThreadHostLoops()` already returns false whenever
`Device::Allows( Backend::DEVICE_MASK )`, so the library turns its host threading
off under a Device by design; the collision is only with *MEQ's* threaded
assembly, which a resident solve does not want. Item 7 stays upstream's and stays
the largest single number on this page **for phase B**, where the assembly is.

**WHAT IS NO LONGER ON THIS LIST, AND WHY.** Host threading of the NPC traversal
([M-135](MEASUREMENTS.md#m-135), 1.25×) and of `Reconstruct()`
([M-172](MEASUREMENTS.md#m-172), 1.177×) are both landed and both are **host**
parallelism that a resident run replaces rather than uses. Nobody should size
the device work against them.

**AND THE ONE PIECE OF ITEM 8 THAT CAN RUN EARLY IS THE EXPERIMENT, NOT THE
SWEEP.** A single-precision **host** build costs a second MFEM install and a
suite run, needs no device and no `real_t` sweep, and answers what the whole of
§7 rests on: does a bordered Newton closing `psi_ax` against a constraint
currently reading `-6.466e-12` converge at all at a 1e-7 floor? **If the answer
is no, item 8 is not work to schedule — it is work not to do**, and the case for
this hardware collapses with it. Run it long before its turn comes.

### 3.4 The hardest part is not an integrator, and the plan used to hide it in a remainder

**`other (remainder)` at 14.0% plus `constraint location` at 3.3% is 17.3% of a
step that is neither dense algebra nor an integrator.** Those are M-126's
threaded shares; §3.3's table is a different, serial run and reads the same leg
at 0.9%. **Neither share is the reason this is item 3** — the crossing count is,
and it is 21.1 a step in both — a proxy since superseded by
[M-198](MEASUREMENTS.md#m-198)'s instrumented **893.2 reads and 0.2
transfers**. It is scalar reductions
feeding **host control flow**: `psi_ax` and `psi_bnd` closing a bordered Newton,
the X-point location, the plasma-support decision, the convergence test. Every
one of them is a device reduction followed by a host branch, and a host branch
on a device value is a synchronisation point.

**This is where residency is won or lost**, because a per-iteration
device→host→device round trip on a scalar costs more than the arithmetic it
guards, and MEQ's bordered Newton has several per step. The techniques are
known — fuse the reductions, keep the branch on the device where the control
flow allows, batch the columns, and accept a fixed iteration count where the
test is cheap to over-run — but none of it is "add a kernel", and none of it is
upstream's. **It is also the part no measurement on this machine can size**
until items 1 and 3 exist.

**`Mesh::FindPoints` sits inside this and is `O(elements x points)` and not
reentrant** (`CLAUDE.md`, *Traps*). The tracer already seeds from
`CriticalPoint::element` and takes the full-mesh scan to zero on the common
path, which is the right shape for a device too; what remains is the
data-dependent last resort in `fitByAngle()` and `faceJump()`.

### 3.4a The barrier is nine parts reduction and one part search, and that is the plan for item 3

**THE CROSSING INVENTORY DECOMPOSES, AND IT DOES NOT DECOMPOSE THE WAY THE
NAMES SUGGEST.** §1.3's counts per Newton step, sorted:

| what | /step | kernel shape |
|---|---|---|
| exterior transmission contraction | **15.8** | **reduction** — a fixed sparse row dotted against the state |
| border row assembly + dense solve | 6.5 | **reduction** — `rowDot()` is the same shape; the `(N+4)²` solve is 14×14 and belongs wherever the rows are |
| plasma current `I_p` | 2.3 | **NOT a reduction — an element-loop quadrature of the nonlinear source.** See the correction below |
| `CriticalPointFinder` axis + X-point | 2.9 | **search** — an element-local Newton per candidate, then an extremum over candidates |
| limiter row | 0 here | **NOT a reduction — a golden-section search per face, and MEASURED NOT WORTH PORTING** → [M-201](MEASUREMENTS.md#m-201): it runs on one shipped case, at 0.7%, causing zero transfers |

**So about 25 of the 27.5 crossings a step are contractions of a fixed operand
against the state, and one item is a search.** Every reduction is the same
kernel: gather the state on a support, multiply, sum. MEQ already owns the data
structure — `meq::CompressedRows`, which stores each row on its support alone
and gathers the state at the same indices in the same ascending order, *"so the
sum below has the same terms in the same sequence as the walk over all `n` it
replaces, and the answer is bit-identical"*. That comment was written for a host
optimisation and it is exactly the contract a device kernel needs.

**THE SEARCH IS THE ONLY GENUINELY NEW THING, AND IT IS THE SMALLEST COUNT.**
This plan has named `CriticalPointFinder` as the barrier since it was written,
because it carries the explicit *"THIS CLASS IS HOST CODE THROUGHOUT"* comment.
It is **2.9 crossings a step against the transmission contraction's 15.8.** *The
component that documents its own transfer is not the one doing most of it.*

**WHAT CAN BE WRITTEN BEFORE ITEM 1, AND WHAT ONLY PAYS AFTER — THE DISTINCTION
MATTERS FOR SCHEDULING AND THIS PLAN HAS NOT DRAWN IT.** A device reduction over
the state is correct as soon as it is written and can be validated against the
host digit for digit today. It does not **pay** until the state is already
device-resident, because under today's partial offload the trace solve and
`NPCRecover` leave the state host-valid, so a device contraction would *add* a
transfer rather than remove one. So:

* **writable and checkable now, against a host reference**: the transmission
  contraction, `rowDot()` and the limiter max. Each has an existing host
  implementation to be pinned against — the same arrangement `SourceViewTests`
  uses for the POD, and the reason that arrangement is trusted here.
* **pays only after item 1**: all of them.

**`I_p` IS NOW BUILT TOO, AS AN ELEMENT LOOP RATHER THAN AS A REDUCTION**
→ **[M-192](MEASUREMENTS.md#m-192)**. `meq::plasmaCurrentPerElementOnDevice()`
agrees with the host loop **exactly** on the debug backend and to **one ulp**
(2.915e-16) on the card, and the pair is the result: the exact zero says the
indexing is right — four column-major layouts, each independently transposable
into a plausible number — and the ulp is M-190's FMA, quantitatively. Built from
`GetElementRestriction()`, `GetGeometricFactors()`, `GetDofToQuad()` and
`mfem::forall`, with no hand-written kernel, so it serves HIP as well as CUDA.
**It is the first kernel for which the M-190 trade was settled the other way**:
contraction left ON, acceptance a round-off tolerance, because the arithmetic is
the work rather than a sum being reproduced.

**AND THE LIMITER IS NOT A REDUCTION EITHER, WHICH MAKES TWICE.**
`meq::limiterPolygonMaximum()` is not a max over a set of values: it runs a
**golden-section search along every limiter face**, evaluating `psi_h` at
arbitrary parameters through `Mesh::GetFaceElementTransformations()` — the
caller-allocated variant of MFEM's shared face scratch — and returns the
argmax's element, shape vector, dofs and physical point so the border row can be
built at it. A device port needs face geometry at arbitrary points rather than
at quadrature points, and sampling instead would choose a different contact and
therefore a different border row. It is the same species as `Mesh::FindPoints`
in §3.4: a **search**, and this section's own decomposition already said the
searches are the hard part.

**SO THE CLASSIFICATION ERROR IS THE FINDING, AND IT RECURRED.** `I_p` and the
limiter were both listed as reductions because their OUTPUT is a scalar. *A
scalar output does not make a reduction; the body does.* Written down once for
`I_p` and then repeated one row later in the same table. **The reduction list had
two entries, they were the same kernel, and they are done**; the other two rows
are an element loop (built, [M-192](MEASUREMENTS.md#m-192)) and a search (not
started).

**AND `I_p` WAS NOT A REDUCTION, WHICH THIS SECTION GOT WRONG ON ITS FIRST
WRITING.** It was listed above as a reduction because its *output* is a scalar.
Its body is not: `assemblePlasmaCurrent()` walks elements, builds the shape
functions, sums `ψ = shape · dofs` per quadrature point and evaluates
`normalisedSource->scaledF( R, z, ψ )` — **a virtual call into the polymorphic
source hierarchy**, which is item 1's wall met from a second direction. Putting
it on a device needs the element geometry, the dof maps and the shape tables
resident, which is item 1's element loop and not a reduction at all. *A scalar
output does not make a reduction; the body does.* The correct reading is that
the reduction list had **two** entries, they were the same kernel, and they are
done.

**That is a real scheduling result and it inverts the obvious order.** Item 3
looked like it had to wait for item 1; its *correctness* does not, only its
benefit does. Since item 1's other half is upstream's and MEQ is waiting on it,
**the reductions are the work MEQ can do while waiting**, and each lands with a
bit-exact host comparison that item 6 would later take away.

**AND THE ACCEPTANCE IS §4's STAGE 4b RATHER THAN A CLOCK.** A reduction that is
correct and does not yet pay is invisible to a timing on this machine and
visible to a crossing count, which is why that stage exists.

**AND ITS BLOCKER IS NOW REMOVED ANYWAY, WHICH IS WORTH HAVING FOR ITEM 1.**
`assemblePlasmaCurrent()` wants `scaledF()` — the PLASMA term alone — and
`makeSourceView()` refused every coil-augmented source, which is **every
production machine case in this tree**. The refusal's reasoning is right for
`f()`, whose value is the sum, and inverted for `scaledF()`, whose value is the
plasma term by definition. `meq::makePlasmaSourceView()` is therefore a second
entry point rather than a widened cast: it unwraps both coil wrappers and
delegates, so nothing is dropped silently because the caller asked for the
plasma term by name. Pinned to the bit against `scaledF()` and `scaledDFdPsi()`
through a real wrapper, with `makeSourceView()` asserted still to refuse, and
with a control that the coil term is non-zero at 13 of the 65 sweep points —
**which failed first time and caught a vacuous test**, the sweep radii having
all missed the coil's own footprint. So the POD now covers the configuration MEQ
exists to solve, and item 1's MEQ half is no longer fixture-only.

**THE FIRST REDUCTION IS BUILT** → **[M-190](MEASUREMENTS.md#m-190)**.
`meq::DeviceRows` mirrors a `meq::CompressedRows` and computes **every row's dot
against the state in one launch**. It covers BOTH consumers — the transmission
contraction and the border's `rowDot()` — which is **22.3 of the step's 27.5**
host reads. One thread per row, sum sequential and ascending inside it, so the
answer is bit-identical to `CompressedRows::dot()` rather than merely close;
`tests/unit/DeviceRowsTests.cpp` asserts `==` on the double, on the debug
backend by default and on the card under `MEQ_TEST_DEVICE=cuda`. Not wired into
the solve, for the reason above: it is correct now and pays after item 1.

**AND IT TURNED UP A CONSTRAINT ON THE WHOLE PLAN, NOT JUST ON ITSELF.**
Matching the summation order is **not sufficient** for bit-identity on a GPU:
nvcc contracts `total += a*b` into an FMA, rounding once where MEQ's host build
— which carries no `-march`, so baseline x86-64 has no FMA — rounds twice.
Measured, one row of many: `3037906131.2962995` on the card against
`...2962999` on the host, five rows failing, **while the debug backend passed**.
`-fmad=false` on that translation unit closes it and both backends then agree
exactly. Three things follow and they apply to every kernel after this one:

* **`Device( "debug" )` cannot see this class**, being host arithmetic. That is
  its second blind spot after M-187's, and it means **a bit-exactness claim
  about a kernel has to be made on the card**.
* ***"Reproduces `--device cpu` to every printed digit"* is weaker than it
  reads.** Printed is six or seven figures; an FMA difference is at 1e-16
  relative. M-79, M-130 and M-187 were all root-caused with that weaker
  instrument and it was enough for them, because the defects they found were
  order-unity. It would not be enough for a kernel that is subtly wrong.
* **So a kernel that wants the strong instrument has to be compiled to give up
  the FMA**, and that is a trade to make deliberately per kernel rather than
  globally. It is obviously right for a reduction that exists to avoid a
  transfer. It would be obviously wrong for the element-local residual of item
  1, where the arithmetic IS the work — so item 1's acceptance is the
  convergence rates of §4 stage 3, not bit-equality, and that distinction
  should be settled before the kernel is written rather than after.

### 3.5 What changes about the refusals, and what "bit-exactness" is actually worth

**`TraceSolver = "cudss"` — the refusal INVERTS under residency, and it is the
one item on this page whose sign this reframing flips.** The refusal rests on
upstream's group-4-alone gate: *"Groups 1 and 4 leave the integrators on the
host, so every iteration would copy the local blocks host↔device around
host-side integrator work."* That argument is **about a partial offload and is
correct about one** — M-173 measures exactly it. On a resident run the premise
is gone: the integrators are on the device, nothing is copied per iteration,
and cuDSS is the only trace solver that does not force a round trip. **cuDSS
becomes the default on a device rather than a key that is withheld.**

#### The instrument this section was ordered around is weaker than it claimed

**THIS SECTION SPENT ITS LENGTH PROTECTING *"`--device cuda` REPRODUCES
`--device cpu` TO EVERY PRINTED DIGIT"*, AND THE RECORD DOES NOT SUPPORT THE
WEIGHT PUT ON IT.** Every device-path defect either project has found was
enormous:

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

**SO WHAT FINDS DEVICE DEFECTS IS A REFERENCE ARM YOU TRUST, AND THE TIGHTNESS
OF THE COMPARISON HAS BEEN NEARLY IRRELEVANT.** This plan has been conflating
two different things, and only the first is load bearing for debugging.

#### What tightness genuinely buys, which is not defect-finding

1. **INTERPRETABLE DIFFS, AND THIS IS THE STRONGEST ARGUMENT.**
   [M-186](MEASUREMENTS.md#m-186) settled a red case by diffing two binaries'
   full transcripts: 192 lines of 377 differed and **every one was a residual
   at 1e-12 or below**, which is what said *no semantic cause to find, do not
   bisect*. That conclusion is only available when both arms are
   deterministic. Against a binary that varies run to run at 1e-15 the same
   diff says nothing and the session goes to a bisect instead.
2. **A TOLERANCE NOBODY CAN RENEGOTIATE.** `tests/unit/DeviceRowsTests.cpp`
   asserts `==` on the double and `DeviceRows.hpp` says that if it is ever
   relaxed to a tolerance the file has stopped being what it claims. A
   tolerance is a number a later session under pressure widens; `==` is not.
3. **PUBLISHED NUMBERS THAT ARE RE-DERIVABLE.** `assemblePlasmaCurrent()`
   refuses `reduction( + : )` so the answer does not depend on a thread count.
   That is about `MEASUREMENTS.md`'s digits being reproducible by somebody
   else, and it is a different argument from debugging.

**NONE OF THE THREE REQUIRES THE DEFAULT TO BE BIT-EXACT. ALL THREE REQUIRE
THAT A REPRODUCIBLE CONFIGURATION EXIST.** That is the operative rule, and it
is a much weaker constraint than the one this section has been enforcing.

#### Two reproducible configurations exist, and that is measured

→ **[M-200](MEASUREMENTS.md#m-200)**, **[M-209](MEASUREMENTS.md#m-209)**. Two
consecutive `--device cuda` runs of `machine-f-diiid` give `psi_ax - max psi_h`
as `-6.465e-12` then `-6.466e-12`. Upstream then swept the three batched keys
**one at a time**, bitwise over four runs:

| batched keys on | reproduces bitwise |
|---|---|
| none | **yes** |
| `LocalFactorMode::Batched`, including the new traversal kernels | **yes** |
| + `TraceAssemblyMode::Batched` | **yes** |
| + `AssemblyMode::Batched` | **no** — 3.5e-16 to 4.7e-16 |

**The batched LU is innocent** — M-200 named it and was wrong — and it is
`AssemblyMode::Batched`'s face kernels, which accumulate `D` through
`AtomicAdd` from three or four faces per element. So MEQ has an all-serial-keys
arm that exercises the whole device **memory** path, which is what root-caused
every defect in the table above, **and** a second arm that additionally
exercises batched **arithmetic**. Validate in either; compare the default
against them to a tolerance.

**AND THE DRIVER'S OWN DEVICE DEFAULT WAS WHAT LEFT THE REPRODUCIBLE SET —
CHANGED** → **[M-215](MEASUREMENTS.md#m-215)**. `apps/meq.cpp` turned on
`AssemblyMode::Batched` under a Device where the file said nothing, which is
the one key of the three that loses bit-reproducibility; it is `Serial` there
now. Two consecutive `--device cuda` runs agree exactly where they read
`-6.465e-12` then `-6.466e-12`, and against the host with the keys held equal
every printed line agrees but one residual at 5.3e-17. **What it gives up is
the face kernels alone**, whose host measurement is a 15 to 27 per cent loss.
*This section's own rule is that a reproducible arm must survive every step,
and the driver was quietly failing it by default.*

#### cuDSS, and why its own determinism cannot be bought back

**cuDSS 0.8.0 leaves `CUDSS_CONFIG_DETERMINISTIC_MODE` at 0**: same matrix,
same right-hand sides, fresh handle, five repeats, **4 of 4 differ**, at
6.7e-16 / 1.8e-15 / 3.6e-15 for `K = 1 / 4 / 16`.

**AND THE SWITCH THAT WOULD RESTORE IT IS UNUSABLE, FOR A REASON SPECIFIC TO
MEQ.** Deterministic mode is supported for **one right-hand side only** on
0.8.0; the solve phase refuses 4 and 16 columns with
`CUDSS_STATUS_NOT_SUPPORTED`. MEQ's bordered Newton applies one factorisation to
`N + 4` right-hand sides through `ArrayMult` — M-98's 1.34× — so **the one
configuration MEQ would want determinism in is the one where it cannot have
it.** 0.8.0 is the newest cuDSS NVIDIA ships, so this is the library as it
stands rather than something about to be fixed.

**DO NOT REACH FOR `SetDeterministic( true )`.** It does not work at `N + 4`
columns; with `CUDSS_MTYPE_SYMMETRIC` it is additionally *silently wrong*, an
exactly zero solution that every phase reports as success, which is why
upstream's setter refuses that type by name; and buying determinism with
column-by-column `Mult()` would mean validating a configuration that is not the
one being shipped. **MEQ is clear of both hazards today and it is worth saying
why rather than assuming**: it builds `NONSYMMETRIC` / `FULL`, and it calls
`SetMatrixSymType()` before `SetMatrixViewType()` at both of its two call sites
— `GradShafranov.cpp:451` and `TraceSolverScaling.cpp:470` — which is the order
upstream's own order-dependence defect requires and which a caller had no way to
know was load bearing.

#### What this changes about the ordering

**ITEMS 6 AND 8 ARE NO LONGER "STRICTLY LAST BECAUSE THEY ARE UNIQUELY
DESTRUCTIVE".** They were ordered that way on the premise that they alone take
the digit-for-digit check away; that premise is gone, because the default
device configuration already has and because the check was never the thing
doing the work. **The constraint that replaces it is one sentence: a
reproducible arm must survive every step, and a step that would leave MEQ with
none of them does not go in.** cuDSS satisfies that — it is a *trace solver*
choice, so a PARDISO arm remains available on the same build — and so does
single precision, whose arm is the double install. What both still owe is
§7.2's harder problem: an acceptance stated as a tolerance rather than as
agreement, and a tolerance loose enough to absorb a 3.6e-15 solve that differs
every run underneath a bordered Newton that differences its constraints through
it.

**The gate is already the right shape and needs no change today.** It is
*refuse without a device*, not *refuse always* (`apps/meq.cpp:1480`), so a
caller who has configured `--device cuda` already gets cuDSS and
`theDriverRefusesASolverItCannotHonour` keeps the config-file-alone route
closed. What changes at item 6 is the **default**, not the refusal.

**What does NOT change:**

* **Partial assembly for the HDG integrators.** Its economic case is sum
  factorisation and MEQ is triangles (§2.4). Residency does not need it —
  batched dense kernels are the route — so this stays out of the plan entirely
  rather than moving up it.
* **`AssemblyMode::Batched` as a HOST default.** M-99's 15–27% loss stands as a
  host measurement and MEQ's host default stays `Serial`. §3.2 says why that is
  not a verdict on the device route.
* **Timing anything on this machine and calling it throughput.** M-77 measures
  the CUDA build costing 7% of the wall on the CPU path for byte-identical
  output. M-223 measures the partial device path's SOLVE at parity with the
  host-threaded default, and the residual 1.07-1.13x on the wall as one host
  leg. Neither is a
  measurement of what a datacentre part would do, and §7 is the only thing that
  makes a local device timing worth taking at all.
* **A second expression of `F( R, z, psi )` before item 3 exists.** Two
  implementations that must agree, with nothing to run one on.

### 3.6 The route was never missing, and this plan spent a section calling it the gate

**`AssemblyMode = "batched"` IS THE DEVICE ROUTE, AND MEQ HAS EXPOSED IT AS A
TOML KEY THE WHOLE TIME.** `DarcyForm::Assemble()` calls
`AssembleFluxMassMatricesBatched`, `AssembleDivMatricesBatched`,
`AssemblePotMassMatricesBatched` and `AssembleFluxMassBdrMatricesBatched`
**in front of** each of the three element loops; each returns false having done
nothing unless `AssemblyMode::Batched` is asked for, and only then does the
caller fall back to `BilinearForm::ComputeElementMatrix()`. So the sentence this
plan quoted as the blocker — `darcyform.hpp`'s note that hybridized assembly
goes through `ComputeElementMatrix()` — describes **what happens when the key is
off**, which is a default and not a wall.

**THE ADMISSION IS NARROW AND IT IS EXACTLY MEQ'S.** `HDGElementMassCanBatch()`
takes `MassIntegrator` and `VectorMassIntegrator`; `HDGElementDivCanBatch()`
takes `VectorDivergenceIntegrator`. `buildForms()` installs
`VectorMassIntegrator( radius )` and `VectorDivergenceIntegrator`. Nothing had
to be written for MEQ to fit; it already did.

**RUN ON THE PRODUCTION CASE RATHER THAN ARGUED**
→ **[M-185](MEASUREMENTS.md#m-185)**. With all three batched keys set,
`machine-f-diiid` prints `face assembly yes, local factorisation yes, local
solve yes, trace assembly yes` and reaches `psi_ax = 3.759851e-01`, the host
default's answer to every digit. **Seven of the ten batched routes are taken and
the one refusal left is `CanBatchLinearResidual()`** — `meq::SourceIntegrator`
being genuinely nonlinear, which is Request 3 and not a route.

**HOW THIS PLAN GOT IT WRONG IS THE TRANSFERABLE PART, AND IT IS THIS FILE'S OWN
RULE BROKEN.** §6 says every claim about somebody else's code here was taken by
extracting whole functions, *"never by reading a line range — `CLAUDE.md`
records a false claim produced exactly that way"*. This one was a **doxygen
comment on a member declaration**, read without reading the function it
describes or the call site in front of it. The comment was true; it was an
answer to a different question. **A comment is documentation of an intent, and
the call site is the fact** — and `grep`ping for the batched entry points from
`DarcyForm::Assemble()` outward is the ten seconds that separates them.

**AND ONE THING THAT LOOKED LIKE A DISAGREEMENT WITH UPSTREAM WAS A
CONFIGURATION.** Upstream's predicate table reports `CanBatchLocalSolve` as
**YES**; MEQ's own diagnostic printed **NO** on the same problem. The predicate
opens `lfac_mode == LocalFactorMode::Batched`, upstream set all three keys and
MEQ set one, and both readings are correct. **Before reporting that a library
disagrees with its own table, check which keys the table was taken with.**

### 3.7 The one item that outranks everything else on this page is upstream's and is not a kernel

**1.74x OF THE DEVICE PATH'S 1.96x IS HOST LOOPS LOSING THEIR THREADS, AND THE
CAUSE IS `general/mem_manager.cpp`.** M-185's `--device cpu` + `batched` control
reproduces the collapse with no device in the process: `residual` and
`NPCReduce`/`NPCRecover` both fall to **1.00 cores**. Upstream has since gated
every host loop on `ThreadHostLoops()` and measured 1.82x to 4.13x back on
those legs — **and made it return false whenever a Device is configured**,
because MFEM's memory bookkeeping is live only then and is not thread-safe.
`AssemblyMode::Threaded` under `--device cuda` **used to abort** in
`MemoryManager::CheckHostMemoryType_`; with the install refreshed it **declines
with a one-time warning and runs serially**, correctly.
[M-188](MEASUREMENTS.md#m-188) has both halves, and the half upstream's gate
cannot reach: MEQ's own element loops hit the same maps through
`GridFunction::GetVectorValue()`, need only two threads to do it, and are gated
in `src/meq/Threading.hpp` on `meq::hostLoopsAreSafe()`.

**So MEQ gains nothing from that fix and loses nothing either**, and the item
that would pay is one upstream names and has not started. It is worth stating
plainly because it reorders this page: **no amount of MEQ's own residency work
recovers 1.74x, and the first 1.74x is not MEQ's to take.** What MEQ's work
buys is the *other* factor — the legs that never reach the device — and the two
compose rather than substitute.

## 4. What would tell us it is working

**`mfem::Device( "debug" )` is the instrument and it has already been used
twice, successfully, on exactly this problem.** It has device memory semantics
with host arithmetic and `mprotect`s the host page, so a raw host read of a
device-valid buffer is a **named fault with a backtrace** rather than a wrong
number. M-79's chain of four and M-91's pair of two were both walked that way,
one fix per rebuild, **each fix moving the fault forward to the next site** —
which is the only thing that distinguishes a fix from a coincidence. It costs
nothing: the driver takes `--device debug` (`apps/meq.cpp:97`, `:983`) and so
does `NpcThreadScaling` (`:604`).

**The acceptance at each stage, and each is a property rather than a
tolerance:**

| stage | acceptance |
|---|---|
| **0. Today's state, as a regression** | **BUILT AND GREEN** → [M-223](MEASUREMENTS.md#m-223). `DeviceParity` is the ctest M-130 argued for — one fixed-boundary case and one bordered free-boundary case, each run under `--device cpu` and `--device debug`, asserting the host's `psi_ax`, `psi_bnd`, X-point and **Newton count** to every printed digit — and it passes at 444.9 s. `debug` needs no GPU, so it runs in any build; the one thing it does NOT cover is real CUDA, and §3.5's *what tightness buys* says why that is a gap worth naming rather than closing with a second slow case |
| **1. §1.2's three candidates** | **the `postProcess()` one is DONE** → M-130: it faults, each fix moved the fault forward, and the root cause is `Memory<T>::Wrap`'s non-owning `Delete()` inside MFEM rather than either hypothesis this plan offered. The `PicardThenNewton` pair is still unrun; `--device debug` on such a fixture either faults — and each fix moves the fault forward — or does not, in which case the sites are safe and the reading was wrong. **Both outcomes are results** |
| **1b. §1.2's thirteen ordering-protected reads** | a **limited** free-boundary case with the plasma support frozen — no `refreshXPoint`, no `refreshPlasmaComponent` — under `--device debug`, which is the configuration where the ordering that currently protects them is absent. A clean run says the ordering holds more widely than the reading suggests; a fault names the first site |
| **2. Threaded assembly under a Device** | **ANSWERED, AND THE ANSWER IS THAT IT CANNOT BE MADE TO WORK TODAY** → [M-188](MEASUREMENTS.md#m-188). The acceptance stood as *the same digits at `OMP = 1` and `OMP = 8` under `debug`*, and that is now **MET for serial assembly**: `machine-f-diiid` reads 7 / 3 / 2 and `2.992753e-12` at both thread counts once MEQ's own element loops are gated on `meq::hostLoopsAreSafe()`. It is **unmeetable for threaded assembly**, because `mem_manager.cpp`'s maps are unguarded and upstream has not started on them; the driver refuses that pair above one thread instead. **A refusal is the resolution of this row, not a deferral of it** |
| **3. A device `SourceIntegrator`** | the analytic fixture ladder, unchanged: `k+1` in `ψ` and `q`, `k+2` in `ψ*`, on `ManufacturedNonlinear` and `SimilarityExponential` at `k = 1…3`, with the host and device rows agreeing to round-off. A device path that changes a rate is a wrong integrand, not a faster one |
| **4. The whole chain** | `psi_ax`, `psi_bnd`, the X-point and the Newton count identical to the host on `machine-f-diiid.toml`, under `debug` AND under `cuda`, at more than one thread count. M-79's own closing row is the template: cpu, cuda and debug agreeing to every printed digit. **MET at `OMP = MKL = 1`** → M-130 — `psi_ax = 3.759851e-01`, X-point `( 1.200929, −0.999491 )`, constraint `−6.466e-12`, 2 steps. The thread-count half is stage 2 and is still open, and **`mhd-rectangle` took 9 Newton steps under a device against 5 on the host** — **walked** → [M-203](MEASUREMENTS.md#m-203). It is a configured `mfem::Device` and nothing else: not GPU arithmetic (`debug` reproduces `cuda`), not `LocalFactorMode`, not `AssemblyMode`, and not the border, since this case has none. Every route predicate is character for character identical between an arm that takes 5 and an arm that takes 9, the answer agrees to **1.05e-14**, and the device run is the host run with a **no-op step interleaved after every contraction** that raises `‖r‖` by a constant **0.47%**. **WALKED, EXPLAINED AND NOW GREEN** → [M-207](MEASUREMENTS.md#m-207), [M-214](MEASUREMENTS.md#m-214). `meq::StageChecksum` is the instrument M-187 prescribed, kept this time, and `meq::ProbedSolver` sits where `mfem::NewtonSolver` keeps the increment to itself: on the same right-hand side the device returns an increment **36× too small**, which is a wrong JACOBIAN. MEQ's own operator had exactly that — `xLoc` written blockwise and never carried up — and with it fixed **M-187's 31 / 16 / 4 is now 7 / 3 / 2, the host's to every digit**. `mhd-rectangle` takes `mfem::DarcyNPCOperator` instead of MEQ's own subclass; **that half is upstream's and is CLOSED at `ce744b84de`**, and it reads **5 steps under `--device debug` against 5 on the host**, every printed line agreeing but the final residual at 5.3e-19 → [M-223](MEASUREMENTS.md#m-223). MEQ's report named `x_loc`, by analogy with its own defect, and the fault was **`b_tr`** — right about the symptom, the idiom and the file, wrong about the member, and what separated them was upstream running the same Newton two ways in one process. *A reproduction is not a diagnosis.* **AND THE ACCEPTANCE IS STATED ON ONE FILE AND IS NOT ENOUGH** → [M-205](MEASUREMENTS.md#m-205), [M-206](MEASUREMENTS.md#m-206): `examples/q-driven.toml` was wrong under any configured Device by a **factor of two** in `psi_ax` — **five** raw seams in MEQ's own outer Newton, each hiding behind the one in front of it, the fourth of which cost convergence RATE rather than correctness — while its inner bordered solve converged in the same 4 steps and satisfied its own constraint to 1e-16 throughout. *A bordered row certifies consistency, not correctness.* **`--device cuda` now reproduces `--device cpu` to every printed digit on it.** This row should name every nonlinear example, not one machine |
| **4b. The crossing count** | **NEW, and it is the acceptance §3.3's restructure implies.** Correctness alone cannot tell a resident solve from a partial one — M-187 is the proof, an answer right to every digit on a run doing 4× the work. So the acceptance is a **count of TRANSFERS of a device-written state per Newton step**, instrumented rather than inferred. **BUILT** → [M-196](MEASUREMENTS.md#m-196), [M-198](MEASUREMENTS.md#m-198): `meq::countedHostRead()` asks `Memory::DeviceIsValid()` before the funnel, so it separates a read from a data movement, and it names its sites by `__LINE__`. Today's reading is **893.2 reads and 0.2 transfers a step** — 3 transfers and 0.888 MB over a whole 12-step run — and every remaining transfer is a real one. **The `3.0 a step` this row carried was [M-198](MEASUREMENTS.md#m-198)'s pre-fix total read as a rate**, 36 over 12 steps, kept after [M-200](MEASUREMENTS.md#m-200) took the run to 3 and printed beside `0.888 MB`, which is a total. *An acceptance stated in the wrong units is an acceptance nobody can fail.* It needs no GPU and makes no speed claim, which is what makes it the one throughput-shaped acceptance this machine can honour. **Its own coverage is part of the measurement** — the first version covered 17 of 36 funnels and reported a number 18× too low — and it sees MEQ's funnels only, not MFEM's |
| **5. Anything about speed** | **not on this machine.** §3.3 |

**AND THE INSTRUMENT HAS AN EXPIRY DATE, WHICH IS NEW.** Every *"to every
printed digit"* in this table is available only while the trace solve is
PARDISO. cuDSS 0.8.0 is nondeterministic run to run at up to 3.6e-15 and its
deterministic mode does not work at the `N + 4` right-hand sides MEQ's bordered
step uses, so item 6 converts each of these rows to a tolerance — §3.5. **Take
every one of these acceptances before item 6, not after**, and do not try to buy
the instrument back with `SetDeterministic( true )`.

**Two things the instrument cannot see, and they are why stage 0 exists.** A
value whose correct answer is exactly zero is indistinguishable from a value
never read back — upstream adopted that criticism from MEQ and it cuts both
ways. And a host-side check is not enough on real hardware: if the solution is
written raw to a host page, a host read returns the right numbers while the
device copy is stale, so the output half of such a defect has to be read back
through a device-aware operation. Under `debug` the protected page turns either
into a fault; under CUDA only the device read-back does.

## 5. What belongs upstream, and what does not

**Already filed and open**: `HDG-NPC-TRAVERSAL-FROM-MEQ.md` — item 1 of §3.2,
the largest measured leg, host OpenMP, explicitly not a device request.

**Worth filing, once, and not before it is checked**: the
`ConformingFaceRestriction` null dereference of §2.4. Upstream's own notes call
it *"a crash where a refusal belongs — report upstream"*, so it is theirs to
file and it may already be filed; MEQ should check rather than duplicate. **A
"filed as" claim in MEQ's own files needs a check and not a citation** —
`CLAUDE.md` records one that was never written.

**Not to be filed**: anything measured against `gf-hdg-linearise-first`'s
unfinished offload work. The standing rule is not to report findings against
work that has not landed, and six of that plan's nine items landed only
recently.

**Not upstream's at all**: §1.2's three candidate sites. They are MEQ's own
aliases, and **no unit test in MFEM's tree can see a caller's long-lived
aliases onto the solution blocks** — upstream says so in its own words about
M-79.

### 5.1 Five corrections MEQ's own files should take — **ALL TAKEN**

Each of these is a claim in MEQ's tree that a reader meets before they meet this
one. All five are now corrected at their source, so the list below is the record
of what was wrong rather than a list of work.

**The line-number ones were taken a different way, and deliberately.**
`CLAUDE.md`'s integrator inventory no longer cites `GradShafranov.cpp` line
numbers at all — it names `buildForms()` and the integrator, which is the
durable half. The numbers this section offered as the correct ones were
themselves stale within three commits, which is the argument.

* ~~**`CLAUDE.md` attributes `darcyFlux` / `potentialGf` / `traceGf` to
  `buildForms()`.**~~ **TAKEN.** They are in the constructor; `buildForms()` is
  the one that gets re-entered, which is why sending a reader there to check
  them matters.
* ~~**`CLAUDE.md`'s integrator inventory cites `GradShafranov.cpp:3077`, `:3098`,
  `:3109`, `:3140`, `:3150–3152`, `:3176–3186`.**~~ **TAKEN, by removing the
  numbers.** Every one of the replacements offered here — `:3899`, `:3920`,
  `:3936`, `:3967` — was itself wrong within three commits; the sites are
  `:4076`, `:4097`, `:4113`, `:4144`, `:4188`, `:4198` today and will not be
  tomorrow. The inventory names `buildForms()` and the integrator instead.
* ~~**`apps/meq.cpp:954–956` has no `M-nn` anchor**~~ **TAKEN.** The `--device`
  comment table is re-taken against M-130 — soloviev 1→1, mhd-rectangle 5→5
  since upstream's `ce744b84de` ([M-223](MEASUREMENTS.md#m-223)),
  limited-tokamak 12→12, **machine-f 2→2 to the host's every digit** — and
  carries the rule the correction is really about: *a device row that looks like
  a cost is a suspect, not a datum*.
* ~~**`tests/performance/NpcThreadScaling.cpp:544–570` states that the solve does
  not survive a Device.**~~ **TAKEN** — it now points at M-130, which is the
  measurement that settles it.
* ~~**`CLAUDE.md` and `ROADMAP.md` both say group 2 "needs a partial-assembly
  rewrite and is not built".**~~ **TAKEN in `ROADMAP.md`.** Half of that is
  still exactly right — there is no
  partial-assembly route and MEQ is simplices, so its economic case does not
  apply to MEQ at all. The other half is stale: what upstream built instead is
  batched element and face kernels behind `CanBatch*` predicates, MEQ passes
  most of them, and upstream's own chain measurement now reads 1.33x at 4096
  elements rather than "worse than doing nothing". **The gate sentence MEQ
  quotes is still in that file, unchanged**, and §3.3 still concludes the same
  thing about `TraceSolver = "cudss"` — but on MEQ's own leg profile and element
  counts rather than on a quotation.

## 6. What could not be verified

Stated plainly, because an admitted gap is worth more than a confident guess.

* **Nothing in this plan was run when it was written.** No build, no ctest, no
  GPU; every number was from `MEASUREMENTS.md`, from upstream's own documents,
  or arithmetic on those, and each is labelled. **Since then M-130 has run the
  instrument §4 names**, on `machine-f-diiid` under `--device debug` and
  `--device cuda`, and the bullets below say which of these gaps it closed.
* ~~**§1.2's three candidate alias sites are a reading, not a reproduction.**~~
  **The third is reproduced and closed** → M-130, and it is upstream's rather
  than either hypothesis offered. The two on the Picard fallback are still a
  reading. **And the reading missed two live sites that a single `--device
  debug` run found in an afternoon**, both of them in §1.2's own subject matter
  — which is the argument for running the instrument before extending the
  sweep.
* **§1.2's thirteen ordering-protected reads are not claimed to be defects.**
  M-91 measures the DIII-D case solving clean under `--device debug` with the
  X-point, the located axis and ten exterior modes live, so on that
  configuration the ordering holds. Whether it holds on a limited, frozen-support
  free-boundary case is a one-command experiment that has not been run.
* ~~**M-91's remaining `postProcess()` abort is attributed to a stack-backed
  `Vector` inside MFEM.**~~ **SETTLED, and neither attribution was right** →
  M-130. It is `Memory<T>::Wrap`'s non-owning, unregistered path meeting
  `Memory<T>::Delete()`'s `std_delete` test under `HOST_DEBUG`, at
  `darcyhybridization.cpp:10164`. Filed upstream with the seven same-pattern
  sites.
* **The threaded-plus-device question is open**, §4 stage 2. M-79's `OMP=8`
  abort has no recorded resolution and the post-fix tables do not state a
  thread count for the runs that pass.
* ~~**`apps/meq.cpp:954-956` asserts that a device raises the Newton step count
  "because the device path's element-local evaluation is inexact where dF/dpsi
  is non-zero", at a 1.7x to 2.0x whole-solve loss on the fixtures.** That
  claim has **no `M-nn` anchor**, and the commit that introduced it reports the
  opposite direction on `examples/limited-tokamak.toml` — 8.96 s host against
  9.77 s device at `OMP=1`, *"no iteration penalty"*, on a case whose `dF/dpsi`
  is not zero. One of the two is describing a configuration the other is not.
  **It is load-bearing for item 7 of §3.2 and it should get an anchor or a
  correction before anything is built on it.**~~ **ANSWERED: it was two unsynced
  reads, not a device property** → M-130. `machine-f-diiid` reads 2 steps on
  both. `mhd-rectangle`'s 5→9 survived that repair, was walked
  ([M-203](MEASUREMENTS.md#m-203), [M-207](MEASUREMENTS.md#m-207)) and is
  **closed**: it takes `mfem::DarcyNPCOperator` rather than MEQ's own subclass,
  upstream fixed it at `ce744b84de`, and it reads **5→5** →
  [M-223](MEASUREMENTS.md#m-223). MEQ's report named `x_loc` by analogy with its
  own defect and the fault was `b_tr` — *a reproduction is not a diagnosis*.
* ~~**`tests/performance/NpcThreadScaling.cpp:544-570` states that "a whole NPC
  solve with an `mfem::Device` configured for CUDA does not compute the right
  answer".**~~ **CORRECTED** — it points at M-130.
* **Upstream's crossover table is quads on their fixture.** Whether MEQ's
  triangles at 450 to 6707 elements land anywhere near it is unknown and is the
  single measurement most worth taking first, on a machine that can produce an
  admissible FP64 number.

## 7. Single precision, and why it is LAST rather than first

**IT IS LAST FOR A REASON THAT IS NOT ABOUT EFFORT.** The payoff from residency
on this hardware is mostly a precision question (§3.1), which makes it tempting
to convert first and measure as you go. **Do not**: the double build is the
*instrument* for the seven items above it. Every one of them is validated by
`--device cuda` reproducing `--device cpu` to every printed digit, which is how
[M-79](MEASUREMENTS.md#m-79) and [M-130](MEASUREMENTS.md#m-130) were root-caused
and the only check that tells *this kernel is wrong* from *this format is
narrower*. Convert first and you debug a device port through a rising round-off
floor with two candidate causes for every disagreement.

**THE HARDWARE DECIDES THAT IT IS WORTH DOING AT ALL, AND NOTHING ELSE DOES.** `CLAUDE.md` records that
this machine's RTX 2070 SUPER runs FP64 at **1/32 to 1/64** of its FP32 rate
where a datacentre part runs it at about 1/2. Under a partial offload that is
nearly irrelevant, because the loss is host serialisation
(→ **[M-173](MEASUREMENTS.md#m-173)**). **Under residency it is the entire
economics**: a resident FP64 solve on this card is competing against eight host
cores on a part deliberately crippled for the format it would be using.

**MFEM SUPPORTS IT AND IT IS ONE CACHE VARIABLE.** `CMakeLists.txt:213-224`
maps `MFEM_PRECISION` to `MFEM_USE_SINGLE` / `MFEM_USE_DOUBLE`, defaulted to
`double` in `config/defaults.cmake:25`. `real_t` follows, throughout the
library.

**cuDSS ALREADY TRACKS `real_t`, SO THE TRACE SOLVER COSTS NOTHING.**
`linalg/cudss.cpp:20` selects `CUDSS_R_32F` under `MFEM_USE_SINGLE` and
`CUDSS_R_64F` otherwise, and every value pointer in that file is `real_t *`.
The device solver is precision-agnostic as shipped.

**AND THE HOST DIRECT SOLVERS DO NOT SURVIVE IT, WHICH IS THE RIGHT TRADE AND
NOT A CASUALTY.** MFEM's UMFPACK wrapper calls `umfpack_di_*` and
`umfpack_dl_*` throughout (`linalg/solvers.cpp:3511` onward) — the `d` is
double, and SuiteSparse ships no single-precision UMFPACK. So a
single-precision install loses `TraceSolver = "umfpack"`, and PARDISO wants
checking the same way. **On a resident run that is no loss at all**: UMFPACK is
the host fallback and the correctness reference, and a resident run's solver is
cuDSS. It does mean the two builds are not interchangeable, which is §7.2.

### 7.1 What MEQ would have to change, and the failure mode is loud

**MEQ IS WRITTEN IN `double` AND HAS BARELY MET `real_t`.** Counted in
`src/meq` alone: **3,913 uses of `double` against 13 of `real_t`**. That is not
13 places to fix; it is a codebase that has never had to care.

**BUT THE HARD-ERROR SURFACE IS MUCH SMALLER THAN THAT COUNT, AND THE
DIFFERENCE IS THE DESIGN PAYING OFF.** `src/meq` carries 33 `double`-returning
`override`s, and **none of them face MFEM**: they are `meq::Source::f`,
`meq::Profile::operator()`, `meq::Coil::scaledF` and their kind — MEQ's own
hierarchies, which `CLAUDE.md`'s *Layout* keeps MFEM-free **on purpose** so that
`Profiles` and `Source` are unit-testable without the library. Both sides of
those signatures are MEQ's, so they stay `double` and compile untouched.

**What actually couples to `real_t` is the twelve classes deriving from an MFEM
type** — `meq::PoloidalFieldCoefficient`, `meq::SourceIntegrator`,
`FrozenSource`, `PicardMap`, `ShiftedResidual`, `TimedOperator`, `TimedSolver`,
`ComponentRefreshed`, `RemainderDatumCoefficient`, `ResidualRecorder`,
`FixedPointResidual`, `DenseJacobianSolver` — and their `Eval` / `Mult`
signatures. A signature that no longer matches its virtual either fails to
compile or **silently stops overriding**, and `-Woverloaded-virtual` plus the
`override` keyword turns the second into the first. MEQ spells `override`
everywhere, so this is a compile-time list rather than a search.

**The number to quote is therefore the class count and not the `double` count**,
and the compiler is the instrument that settles it — §7.4.

**The quiet surface is narrowing arithmetic**: a `double` local fed from a
`real_t` field and back. It compiles, it runs, and it costs precision without
saying so. The mechanical part is a `double` → `real_t` sweep of everything
that touches an `mfem::` type; what must stay `double` is everything that does
not — `Config`, `Profiles`, `Source` and `SourceFactory` are MFEM-free **on
purpose** (`CLAUDE.md`, *Layout*), the NetCDF output is an interchange format,
and the TOML parse is `double` because the file is.

**AND A TREE-WIDE TYPE SWEEP IS EXACTLY THE OPERATION `CLAUDE.md` RECORDS
GOING WRONG.** The `r` → `radius` rename touched ~1,300 sites, built clean,
passed `naming` and `sphinx -W`, and **changed six answers**. The instruments
that found it are the ones to reach for here before starting: a second tree
configured `-DCMAKE_CXX_FLAGS="-Wshadow -Winit-self -Wconversion"` built
**clean**, `py_compile` over the tools, fenced-block diffs on the captured
transcripts, and a `git worktree` of `HEAD` to tell *this was already failing*
from *I broke it*. `-Wconversion` is the addition for this sweep and it will be
noisy; that is the point.

### 7.2 The convergence tests are the acceptance, and they will move

**RUN, AND THE PREDICTION BELOW IS INVERTED** → **[M-217](MEASUREMENTS.md#m-217)**.
The tree builds and runs against `install-single`: **38 files, 405 lines**, and
the suite reads **44 red of 63** against the double build's 62/63 on the same
working tree. What this section expected to dominate does not, and what does
dominate is not an assertion at all.

| class | count | |
|---|---|---|
| **A** the SOLVE did not converge | **5** (was 23) | one solver constant, not a tolerance |
| **B** names an absent library | 6 | SUNDIALS and the device, absent from the install |
| **C** a rate read as a failed rate | **8** | what this section predicted would dominate |
| **D** a value assertion below the floor | **25** | what this section called the biggest cost |

**THE BLOCKER WAS `newtonRelativeTolerance`, AND NO SWEEP OF THE ASSERTIONS
WOULD HAVE FOUND IT.** At `1.0e-12` it is seven orders below single's residual
floor, so `solve()` **refuses** and 23 of 45 tests failed before reaching the
assertion they exist for. `meq::newtonTolerance()` floors it and 18 tests moved
into a class that can be assessed — *the count moved by one and the inventory
moved by eighteen.* **A tolerance tuned against a solve that did not converge is
tuned against nothing**, so the sweep below could not have been done first.

**AND THE INSTALL WAS NOT THE DOUBLE ONE MINUS PRECISION. IT IS NOW** →
**[M-218](MEASUREMENTS.md#m-218)**. SuiteSparse and gslib are missing
**structurally** — neither has a single-precision build — and CUDA, cuDSS and
SUNDIALS were missing **incidentally**. All three incidental ones are restored,
which cost a single-precision SUNDIALS (MFEM has no guard that
`SUNDIALS_PRECISION` matches `real_t`) and **six `double`/`real_t` signatures in
`linalg/sundials.{hpp,cpp}` that stop MFEM compiling in that pairing at all**,
reported as `HDG-SUNDIALS-SINGLE-PRECISION-FROM-MEQ.md`.

**AND IT TURNED UP A SEVENTH DIFFERENCE THAT NO LIST OF `MFEM_USE_*` FLAGS
CONTAINS**: PARDISO was resolved against `libmkl_sequential` in
`/opt/intel/.../lib/intel64` where the reference uses `libmkl_gnu_thread` in
`.../lib`, and `CUDA_ARCH` sat at the cache default `sm_60` on a `sm_75` card.
PARDISO is the default trace solver, so that is the knob
[M-78](MEASUREMENTS.md#m-78) measured as deciding which problems converge.
*A configuration stated in flags is not a configuration.*

**AT THIS STEP THE SUITE READS 40 RED OF 63 AND THE DOUBLE CONTROL IS 62/63**
— it reaches **62/63 in both** at [M-222](MEASUREMENTS.md#m-222), and the
progression is kept because what each step cost is the finding.
Against the parity install it was 46, and the six that cleared are
`SolovievConvergence`, `FieldConvergence`, `RotatingSolovievConvergence`,
`RotatingSourceConvergence`, `MaschkePerrinConvergence` and `cuDSSTraceSolver`.
**The count is the wrong instrument here for the third time**: restoring the
install took the red count UP by two while every *"MFEM was built without
`MFEM_USE_SUNDIALS`"* and *"no Device is configured"* disappeared, because a
test that aborted at its first case now runs all of them —
`SourceViewTests` went from 15 error lines to 436, `DeviceRowsTests` from 4 to
272. Read the distinct failing CASES: **173 then, 157 now.**

**AND THE RATE SWEEP THIS SECTION DESCRIBES IS ANSWERED, IN THE NEGATIVE.**
Against the double build on the same tree, `Solov'ev NSTX` in the flux is
1.09× contaminated at `k = 1` on the finest mesh, **203×** at `k = 2`, and
**1.5× on the COARSEST mesh** at `k = 3`. So `k+1` above `k = 1` is a
double-build acceptance: there is no mesh sequence a single build can measure it
on, and three separate per-pair criteria were tried and each one provably fails.
**MEQ's rate acceptance is not a thing single precision verifies; agreement with
the double build is.** That is the finding this section was for.

**AND THE SWEEP IS NOW MOSTLY DONE, ON FOUR LAWS RATHER THAN A TABLE OF
CONSTANTS** → **[M-219](MEASUREMENTS.md#m-219)**, **[M-220](MEASUREMENTS.md#m-220)**.
What this section called *"a large, careful edit with a real risk of relaxing an
assertion into vacuity"* turned out to be four statements about where a floor
comes from, each measured, each expressed once in
`tests/convergence/RoundOffFloor.hpp` and reached from every caller:

* *a difference is not a derivative* — `eps` where two routes are exact,
  `eps^(2/3)` to `sqrt( eps )` where one is a finite difference;
* *a tolerance a solver was GIVEN is the floor on what it can achieve*, and a
  FIXED-POINT iteration cannot be given a Newton's, since it approaches the
  floor from above and never crosses it;
* *a stationary point's POSITION costs half the digits and its VALUE costs
  none* — `sqrt( eps )`, and the envelope argument that hides the one keeps the
  other;
* *a quantity obtained by differentiating rises as the mesh falls* — an
  estimator floors at `eps/h`, a flux-surface average at `eps/h²`, so the floor
  and the truncation CROSS and refining past the crossing makes the number
  worse. That is also why an adaptive loop in this format stops being a descent.

**AND THE SAFETY PROPERTY IS WHAT MAKES THE EDIT SAFE RATHER THAN CAREFUL
READING.** Every helper in that header either `std::max`es against the caller's
own historical value or is `#ifdef`'d to be inert in double, so a guard applied
in the wrong place cannot weaken the double build — and the double control,
62/63 at every step, says so. The risk this section named is closed by
construction rather than by attention.

**AND IT IS NOW DONE: THE SINGLE SUITE IS 1 RED OF 63 AND IT IS THE STANDING
DELIBERATE ONE** → **[M-222](MEASUREMENTS.md#m-222)**, with the double control
at 62/63 throughout. **What the last stretch cost was not tolerances**, which is
the finding this section did not predict: the largest item was
`ElementTransformation::TransformBack`'s default argument sitting at **0.84 of
one `float` epsilon**, so the inverse map returned `Unknown` for the element a
point was actually in and MEQ's element walk fell back on an `O(elements)` scan
108, 240 and 95 times over three surfaces. The rest was `mfem::VertexConePath`
differencing a displacement below one ulp, and **four cases of M-186's coin
toss** — measured with `MKL_NUM_THREADS` rather than inferred, and every one
repaired by relaxing a PRECONDITION rather than a tolerance.

Three laws join the four above, each a sibling of one of them: *a position
obtained by ROOTING carries `eps` where one obtained by DIFFERENCING carries
`eps/h`*; *an answer is worth the residual the solver STOPPED at, which is not
the one it was asked for*, since KINSOL rests on a step tolerance as readily as
on a function norm and reports `converged` for both; and *the floor and the
ceiling on one Newton order fail in opposite directions*, the floor being a
statement about the Jacobian and the ceiling about the fit.

**WHAT REMAINS IS THE SWEEP THIS SECTION DESCRIBES, AND IT IS NOW SIZED.**
Class C is eight rate assertions, all at `k = 2` and the finest `h`, which is
exactly the flattened tail described below. Class D is 25, of which several are
class A's residue reaching a different test. Two are structural and will not
yield to a tolerance: `SourceViewTests` pins the device POD to the **bit**
against MEQ's MFEM-free `Source`, which is `real_t` against `double` by design;
and `SolovievGeometryConvergence` says the floor is **per problem**, as
[M-204](MEASUREMENTS.md#m-204) already found over a range of 560×.

**THIS IS THE PART THAT IS NOT MECHANICAL.** MEQ's acceptance criteria are
**measured convergence rates** — `k+1` for `ψ`, `k+2` for `ψ*` — asserted over
dyadic refinements. Single precision has about **7 decimal digits**, so the
round-off floor rises from ~1e-15 to ~1e-7, and every rate table truncates
earlier: `CLAUDE.md` already warns that rates flatten at `k ≈ 5` or `h/8` from
round-off in **double**. In single that ceiling arrives several refinements
sooner, and a rate assertion that reads a flattened tail as a failed rate goes
red for a correct solver.

**So the suite cannot simply be run against a single-precision build and
believed.** What is wanted is a per-test floor that follows `real_t`'s epsilon,
and the honest statement that the single build is a **device production**
configuration validated against the double build, not a replacement for it.
Concretely: **two installs, and they answer different questions** — `install`
(double, host, UMFPACK/PARDISO) stays the correctness reference and the suite's
home, and a single-precision device install is measured against it.

**And the tolerance-based tests need a sweep of their own.** The tree is full
of assertions at 1e-14, 1e-16 and 2.5e-12 — `theDriverSolvesOnACurvedBoundary`
at 1.6e-16, `DriverAcceptance` at 1.189e-16, `theBoundaryBendsOntoTheTrueGamma`
at 1.1e-16. **Every one of those is below single precision's epsilon**, so they
do not "get looser", they become **unsatisfiable**. That is a large, careful
edit with a real risk of relaxing an assertion into vacuity, and `CLAUDE.md`'s
testing stance forbids exactly that. It is the single biggest cost in this
section and it is MEQ's alone.

### 7.3 RUN, and what it measured

**The experiment this section called for has been taken** →
**[M-174](MEASUREMENTS.md#m-174)**. `../mfem/install-single` exists —
`MFEM_PRECISION = single`, host, PARDISO, no SuiteSparse / SUNDIALS / GSLIB /
CUDA — MEQ is ported to it, and `soloviev-nstx` solves. Four results, in the
order they matter:

1. **MEQ runs at FP32.** `‖r₀‖` agrees with the double build to every printed
   digit; the residual floor rises from `7.5e-13` to `2.2e-06` relative.
2. **Out of the box nothing converges, and the cause is one constant.**
   `newtonRelativeTolerance = 1.0e-12` (`GradShafranov.cpp:896`) is seven orders
   below the FP32 floor, so even a Solov'ev case needing **one** Newton step
   fails. At `1.0e-5` it converges in two. **§7.2 was right that the tolerance
   sweep is the biggest single cost, and it is also the FIRST thing, not a
   clean-up afterwards** — before it, every FP32 number is a measurement of a
   default.
3. **The port is 45 sites in 11 files, not 3,913**, and it is **inert in
   double** — pre-port and post-port binaries give `9.674113e-14` byte-identical
   at `MKL = OMP = 1`. §7.1's revised reasoning holds.
4. **MFEM's own HDG code needed two fixes**, and one of them is `MFEM_VERIFY(err
   < 1e-10)` on an identity — unsatisfiable at FP32's `1.2e-07` epsilon. Both
   are reported upstream.

### 7.4 The blocker this section did not predict: gslib is double-only

> **RESOLVED, AND BY WORK DONE FOR A DIFFERENT REASON.**
> `meq::FieldTransfer` is `meq::PointLocator` — MEQ's own locator, reproducing
> gslib to every printed digit → [M-176](MEASUREMENTS.md#m-176) — so
> `MFEM_USE_GSLIB` is optional and `../mfem/install-single` is built with it
> **off**. The warm start survives single precision. **The paragraph below
> stands as the reasoning that made the locator worth writing**, and its
> conclusion that the bordered question was "still open" no longer follows from
> this cause.

**`fem/gslib.hpp` is raw `double` throughout**, so `MFEM_USE_GSLIB` does not
survive `MFEM_PRECISION = single` — and with it goes
`meq::FieldTransfer`, the interpolating warm start. That matters more than it
sounds, because **every bordered free-boundary case in the tree needs its
guess**: `machine-f-diiid` and `limited-tokamak` both fail cold **in double**,
reporting a singular bordered Jacobian, so the guess is load-bearing physics and
not a convenience.

**So the question this plan most wanted answered — does the bordered Newton
survive a 1e-7 floor — WAS still open, and was blocked on something concrete**
rather than on effort: either a single-precision point locator, or a bordered
fixture that converges from cold. **MEQ built the locator, for a different
reason, and the question is now ANSWERED — §7.5.**

**AND THE LOCATOR BLOCKER IS SMALLER AND STRANGER THAN "GSLIB IS DOUBLE"** →
**[M-175](MEASUREMENTS.md#m-175)**. Three facts change what to do about it:

* **The precision blocker is mostly MFEM's own code, not the C library.**
  `fem/gslib/` is **6,019 lines carrying `mfem::forall` device kernels, with 542
  `double` and zero `real_t`**. The external library (9,440 lines) has no `real`
  typedef at all, so it is a fork rather than a port — but MFEM's half is a
  mechanical `real_t` sweep and is upstream's.
* **MEQ can never reach those device kernels anyway.**
  `fem/gslib.cpp:1388` gates the device path on `tensor_product_only`, which
  requires `QUADRILATERAL` or `HEXAHEDRON`. **MEQ is triangles everywhere**, so
  the gate is false and the host path runs — splitting every triangle into
  **three** quads (`NEsplit = 3`) and resampling the field onto them.
* **What MEQ needs is point location plus FE evaluation, and MEQ already owns
  the expensive half.** `meq::GridSampler` has the affine triangle inverse
  (69–81x over `TransformBack`, [M-93](MEASUREMENTS.md#m-93)) and the
  element-grouped evaluation. The missing piece is a spatial index for
  arbitrary query points — a uniform bin grid — and the whole thing estimates at
  **300 to 500 lines**, device-shaped, sharing item 3's
  `MFEM_HOST_DEVICE` shape-evaluation blocker rather than competing with it.
  `SamplerConvergence` is already the harness that would accept it.

**This is a candidate, not a decision.** It is a replacement for MEQ's *use* of
gslib and not for gslib, the curved path still needs the Newton fallback, and
the warm start is load-bearing physics — so a locator that is subtly wrong is a
solver that converges quietly to the wrong equilibrium.

**Two MEQ defects fell out of building in an unusual configuration**, neither
about precision: `picardSolver` declared under `MFEM_USE_SUITESPARSE` but used
under `MEQ_HAVE_DIRECT_TRACE_SOLVER`, so **PARDISO-without-SuiteSparse does not
compile**; and `FieldTransfer` naming `FindPointsGSLIB` unguarded. *A build in a
configuration nobody uses is a cheap instrument for finding the guards that were
never exercised.*

### 7.4b Beyond CUDA: `rocsolver_csrrf_*` is the likely route, and the missing
first factorisation is the design

**MARKED AS THE MOST LIKELY AMD PATHWAY, NOT AS WORK.** Everything MEQ needs on
a device is already vendor-neutral *except the trace solve*:
`linalg/batched/gpu_blas.cpp` is one source behind a `MFEM_CU_or_HIP(stub)`
macro (cuBLAS or hipBLAS), `BatchedLinAlg` also offers MAGMA which is dual
vendor, `fem/darcy/`'s 25 `forall` sites are `MFEM_HOST_DEVICE`, and MEQ's own
two device translation units guard on `MFEM_USE_CUDA_OR_HIP_LANG` with CMake
reading the backend out of MFEM's `config.mk`. **There is no `hipDSS`**, so the
trace solve is the single CUDA-locked piece.

**`rocsolver_csrrf_*` is the candidate, and it has no initial factorisation ON
PURPOSE.** The family is analysis, `sumlu`, `splitlu`, `refactlu` /
`refactchol` and `solve` — refactorisation given an existing pattern. It mirrors
cuSOLVER's `cusolverRf`, whose design is the same: **the host performs the first
LU and the device refactorises every subsequent matrix of that sparsity.** So
the absent "factor from scratch" routine is not a hole to be worked around.

**AND THAT SHAPE IS ALREADY MEQ'S.** The trace system's sparsity is fixed across
Newton steps — which is exactly why `direct-solver-symbolic-reuse` is a branch of
`meq-integration`, and why both `UMFPackSolver` and `PardisoSolver` carry
`SetReuseSymbolic()` and a `RetainedSparsityPattern`. MEQ already does one
analysis and many numeric factorisations, and already owns a host solver that
could supply the first LU. A `csrrf` trace solver would slot behind the existing
`TraceSolver` interface rather than needing a new one.

**NOTHING HERE IS VERIFIED.** There is no ROCm on this machine — no
`/opt/rocm`, no `rocsolver` headers — so the API shape above is from
documentation and not from a compile, MFEM has no wrapper for any of it, and
MEQ has no AMD hardware to test on. **It is recorded as the direction to take
if extending beyond CUDA becomes a goal, and as nothing more.**

### 7.5 ANSWERED: the bordered Newton survives FP32, and the cost is one constant per file

**RUN** → **[M-204](MEASUREMENTS.md#m-204)**, on a single install rebuilt from
the current `meq-integration`. Both bordered machines, each against the **double
build stopped at the same tolerance**, which is the control that makes the
answer mean anything:

| | its | `psi_ax` single vs double | constraint, single |
|---|---|---|---|
| `limited-tokamak`, one border row | 10 | **1.50e-04** | −1.240e-07 |
| `machine-f-diiid`, four border rows | 1 | **2.7e-07** | −1.863e-07 |

**Three conclusions, and the third is the one that changes what §7.2 has to
do.**

1. **The border closes to the format's floor.** `−1.24e-07` and `−1.86e-07`
   against FP32's `eps = 1.19e-07` say the same thing double's `8.19e-16`
   against `2.22e-16` says. The bordered Newton is not disqualified by FP32, and
   its residual history is the double build's to three and four figures for ten
   steps.
2. **`NewtonRelativeTolerance` is the failure, and it is per FILE.**
   [M-174](MEASUREMENTS.md#m-174) found this on a Solov'ev case; on a bordered
   one at the shipped `1.0e-10` the run stalls at `3.3e-07`, falls into
   `PicardThenNewton`, and burns **400 iterations** before exiting 2. Every
   `examples/*.toml` carries such a constant.
3. **THE ACCURACY IS A PROPERTY OF THE CASE AND RANGES OVER 560×.** The harder
   machine is the more accurate one — four border rows, a diverted topology and
   a moving support at **2.2 × eps**, against the simple limited case at
   **1,260 × eps**. So §7.2's *"a per-test floor that follows `real_t`'s
   epsilon"* is not enough: a floor that admits `limited-tokamak` at 1.5e-04
   admits nothing about `machine-f-diiid` at 2.7e-07, and the honest form is a
   per-case qualification measured against the double build at the same
   stopping point.

**And the single build is a DRIVER, not a suite.** `src/meq` and `apps` compile
clean against `install-single`; `tests/convergence/ConductorSubtraction.cpp`
does not — eight errors, all `real_t *` against `double *` and
`std::max( double, float )`, in a test written after M-174's 45-site port. The
port is not finished, and finishing it is §7.2's cost rather than a new one.

### 7.6 What this still does not settle

**THE NUMERICAL QUESTION IS SETTLED AND THE PHYSICS QUESTION IS NOT.** §7.5
answers *does the border tolerate a 1e-7 floor* — it does, closing to `eps` on
both machines — and that was the piece this section said could start today and
would tell us whether the rest of the plan is worth writing. It is.

**Three things it does not settle, in the order they will cost.**

1. **THE RATE TABLES, WHICH ARE MEQ'S ACTUAL ACCEPTANCE.** §7.2 is unpaid:
   nothing in `tests/` has been run against `install-single`, and
   `tests/convergence/ConductorSubtraction.cpp` does not compile there. A rate
   study asserting `k+2` over dyadic refinements truncates several refinements
   sooner at a 1e-7 floor, so the question is not *does it pass* but *which
   pairs are still above the floor at each `k`* — and that has to be measured
   per test before any tolerance is touched, or the edit is the vacuity
   `CLAUDE.md`'s testing stance forbids.
2. **WHETHER 1.5e-04 ON `psi_ax` IS ACCEPTABLE PHYSICS**, which is a question
   for the person reading the equilibrium and not for this file. §7.5 measures
   it as a per-case number spanning 560×, so the answer cannot be a single
   build-wide claim.
3. **WHETHER FP32 PAYS ON A RESIDENT SOLVE**, which is what it was for. Nothing
   here is a timing: [M-173](MEASUREMENTS.md#m-173) and
   [M-184](MEASUREMENTS.md#m-184) both measure a **partial** offload, where the
   loss is host serialisation and the format is nearly irrelevant. The rate
   question above is the gate on that, not the other way round.
