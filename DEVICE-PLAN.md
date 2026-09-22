# DEVICE-PLAN.md — what a fully device-resident MEQ solve would need

**This is an audit, not a commitment.** It answers four questions separately:
what MEQ's own code would have to change, what is still missing in the MFEM it
links, what order the work has to happen in, and what would tell us it is
working. It draws no performance conclusion from a timing taken on this
machine, per `CLAUDE.md`'s standing rule about the RTX 2070 SUPER's FP64 rate.

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
3. **MEQ's Newton is a BORDERED Newton, and the border is host arithmetic that
   does not thread and would not go to a device.** M-126 measures it at 30% of
   a threaded step (§1.4, §3.1).

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

### 1.3 The host-only components, and which cost per step

The expensive distinction is per-step against once. MEQ's own
`GradShafranovSolver::StepProfile` already names the per-step ones as legs
(`src/meq/GradShafranov.hpp:3411-3461`), which is how M-126 could measure them.

**Per Newton step — these force a device→host transfer every step:**

| component | leg | what it is |
|---|---|---|
| `meq::CriticalPointFinder` | `axisSeconds`, `axisSweepSeconds`, `xPointSeconds` | the located-axis and X-point constraints are roots of the **solved** `q`; the constructor takes both fields to the host at `CriticalPoints.cpp:150-151` |
| the limiter row | `limiterSeconds` | a max over a marked surface or a point evaluation |
| the plasma-current row | `currentSeconds` | a quadrature of `F/R` over the support |
| `meq::ExteriorDtN`'s transmission rows | `transmissionSeconds` | a boundary-face quadrature at order 40 (`GradShafranov.cpp:3501`, `:3570`) |
| the border's dense algebra | `borderAssemblySeconds`, `borderSolveSeconds` | `rowDot()` and an `(N+4)²` dense solve, with the seven `HostRead()`s of §1.2 |
| `meq::PlasmaComponent`'s flood fill | `componentSeconds` | a graph traversal over element adjacency (`PlasmaComponent.hpp:157`) |
| re-assembly | `prepareSeconds` | §1.4 |

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

## 3. The ordered list, with the blocking dependencies

### 3.1 The Amdahl ceiling, and it is set by the legs that stay on the host

*Arithmetic on M-126's published shares, not a new measurement.*
`examples/machine-f-diiid.toml`, 4848 elements, `k = 2`, PARDISO,
`AssemblyMode::Threaded`, `OMP = MKL = 8`, `OMP_WAIT_POLICY=passive`, step
total 0.692 s:

| leg | s | share | could a device take it? |
|---|---|---|---|
| residual | 0.086 | 12.4% | only with §1.1 |
| gradient | 0.045 | 6.5% | partly — `ComputeH` yes, the integrators only with §1.1 |
| trace factorisation | 0.058 | 8.4% | yes, cuDSS |
| trace backsolve | 0.081 | 11.7% | yes, cuDSS |
| `NPCReduce` + `NPCRecover` | 0.212 | **30.6%** | yes in principle — dense, integrator-free |
| constraint location | 0.023 | 3.3% | **no** — `CriticalPoints` |
| border assembly | 0.047 | 6.8% | **no** |
| border dense solve | 0.039 | 5.6% | **no** |
| re-assembly | 0.003 | 0.4% | no |
| other (remainder) | 0.097 | 14.0% | **effectively no** — scalar reductions feeding host control flow |

**Device-able 69.8%, host-bound 30.2%. Taking the device-able legs to zero is
3.31x on the step and no more.** At upstream's own measured whole-chain gains
it is far less, and those are quad numbers on their fixture rather than MEQ's:

| gain on the device-able legs | step |
|---|---|
| 1.33x (upstream at 4096 el) | **1.21x** |
| 1.68x (upstream at 16384 el) | 1.39x |
| 3.00x | 1.87x |

**On the whole run it is weaker again.** M-102 splits a cold DIII-D run at
5.245 s into solve 3.733, gmsh 0.684, output 0.788, setup 0.019 — so **28.8% of
the run is not the solve at all** and is entirely host. Taking the solve's
device-able share to zero is **1.99x on the run**; at upstream's 1.33x it is
**1.14x**, and at 1.68x it is **1.25x**.

**Say that plainly: the ceiling on any device work MEQ could do is about 2x on
a whole run, and the realistic figure at upstream's own measured gains is 1.1
to 1.3x — before the triangle and crossover caveats of §0.**

### 3.2 The order, and who owns each item

| # | item | whose | blocked by |
|---|---|---|---|
| 1 | **`NPCRecover` and `NPCReduce` threaded on the HOST** | upstream | nothing. Already requested as `HDG-NPC-TRAVERSAL-FROM-MEQ.md`. Worth **1.31x on the step and 1.24x on the run** by MEQ's own arithmetic, and it is the same 30.6% leg any device path would target |
| 2 | **MEQ's three candidate alias sites, and a guard for the thirteen ordering-protected reads (§1.2)** | **MEQ** | **PART DONE, M-130.** The `postProcess()` candidate is walked and is upstream's, filed; the two on the Picard fallback are not. Two sites this audit did not list — `prepare()` writing through an alias, and `essentialTrace[ i ]` — are fixed, and `machine-f-diiid` now solves under `--device cuda` in 2 Newton steps to the host's every digit. The thirteen still want one funnel `HostRead()` each |
| 3 | **A device route through `DarcyHybridization` that reaches the batched kernels instead of `ComputeElementMatrix()`** | upstream | items 4–9 of upstream's own plan, six of which are done. This is the "route" half of §2 |
| 4 | **A POD source description and a `MFEM_HOST_DEVICE` `F( R, z, psi )`** | **MEQ** | item 3. Without a route it is unreachable; with one it is 12–23% of a bordered step |
| 5 | **`mfem::HDGExtensionIntegrator` on a device** | upstream | items 3 and 4. Needed only for the curved and free-boundary paths — which is every problem MEQ exists for, and none of M-80's fixtures. **It does not disqualify the condensation cache**: it is a `BilinearFormIntegrator` (`extension_hdg.hpp:509`) landing in `A`, assembled once; it makes `A` asymmetric and not solution dependent |
| 6 | **The simplex face-restriction crash** | upstream | nothing, but it gates any tensor-only device geometry ever reaching MEQ |
| 7 | **`TraceSolver = "cudss"` from a config file ALONE, with no `--device`** | MEQ | items 3 and 4, and it is the LAST thing to open rather than the first — §3.3. With `--device cuda` it is already open |

**Item 1 is first and it is not a device item.** It is the largest measured leg,
it is host OpenMP with no integrator and no `ElementTransformation` in it, it is
already asked for, and it buys more than the whole device chain buys upstream at
MEQ's element counts. Anything that starts with the device starts behind it.

**Item 2 is MEQ's only unconditional item.** It costs nothing and it is the
difference between "the device path is correct" and "the device path is correct
on the fixtures we happened to run".

### 3.3 What is not worth starting, and why

* **`TraceSolver = "cudss"` from a config file — the refusal stands.**
  Upstream's gate is unchanged and still written: *"Two of the four groups are
  nearly free and doing only those is worse than doing nothing. Groups 1 and 4
  leave the integrators on the host, so every iteration would copy the local
  blocks host↔device around host-side integrator work."* MEQ's own M-126 puts
  the trace solve at **20.1%** of a threaded step (factorisation + backsolve),
  so group-4-alone is buying a share of a fifth while paying a round trip on
  the other four fifths. **The gate is already the right shape**: it is
  *refuse without a device*, not *refuse always* (`apps/meq.cpp:1480`), so a
  caller who has deliberately configured one with `--device cuda` — and is
  therefore already paying for every `Vector` to live there — gets cuDSS, while
  the config file alone cannot reach it and
  `theDriverRefusesASolverItCannotHonour` keeps that closed. Nothing here
  proposes changing either half.
* **Partial assembly for the HDG integrators.** Its economic case is sum
  factorisation and MEQ is triangles (§2.4). This is a discretisation decision,
  not a flag, and nothing in MEQ should wait for it.
* **`AssemblyMode::Batched` as a host default.** Measured at a **15 to 27 per
  cent LOSS** on the host in all three rounds of M-99. It is already the
  default *on a device* (`apps/meq.cpp:1376-1378`) and that is the right place
  for it.
* **`LocalFactorMode::Batched`.** M-101: it buys a batched local factorisation
  and pays the condensation cache, since `CanCacheCondensation()` refuses
  outright under it. M-99's flat axis is two effects of opposite sign.
* **Timing anything on this machine.** M-77 already measures the CUDA build
  costing **7% of the wall clock on the CPU path** for byte-identical output,
  with 55.4 M of 226.7 M allocations attributable to `mfem::forall`'s host
  lambda wrapper. A consumer FP64 rate of 1/32 to 1/64 can invert a production
  conclusion.
* **A second expression of `F( R, z, psi )` before item 3 exists.** It is two
  implementations that must agree, and there is nothing to run it on.

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
| **0. Today's state, as a regression** | **STILL OPEN, and M-130 is the argument for it**: the two links it found were latent for as long as nothing ran a bordered free-boundary case under a Device. There is **no ctest that runs a whole solve under a Device** — `tests/CMakeLists.txt:323` registers exactly one, `cuDSSTraceSolver`, and it checks the trace solver's agreement alone. The first thing to build is a `--device debug` run of one free-boundary case and one fixed-boundary case asserting the CPU answer to every printed digit. `debug` needs no GPU, so it runs in any build |
| **1. §1.2's three candidates** | **the `postProcess()` one is DONE** → M-130: it faults, each fix moved the fault forward, and the root cause is `Memory<T>::Wrap`'s non-owning `Delete()` inside MFEM rather than either hypothesis this plan offered. The `PicardThenNewton` pair is still unrun; `--device debug` on such a fixture either faults — and each fix moves the fault forward — or does not, in which case the sites are safe and the reading was wrong. **Both outcomes are results** |
| **1b. §1.2's thirteen ordering-protected reads** | a **limited** free-boundary case with the plasma support frozen — no `refreshXPoint`, no `refreshPlasmaComponent` — under `--device debug`, which is the configuration where the ordering that currently protects them is absent. A clean run says the ordering holds more widely than the reading suggests; a fault names the first site |
| **2. Threaded assembly under a Device** | M-79's `OMP_NUM_THREADS=8` abort inside `MultNL`'s own OpenMP region — *"host pointer is not registered"* — has **no recorded resolution**, and every post-fix table in M-79 and M-91 is at one thread or does not state a count. The acceptance is the same digits at `OMP = 1` and `OMP = 8` under `debug` |
| **3. A device `SourceIntegrator`** | the analytic fixture ladder, unchanged: `k+1` in `ψ` and `q`, `k+2` in `ψ*`, on `ManufacturedNonlinear` and `SimilarityExponential` at `k = 1…3`, with the host and device rows agreeing to round-off. A device path that changes a rate is a wrong integrand, not a faster one |
| **4. The whole chain** | `psi_ax`, `psi_bnd`, the X-point and the Newton count identical to the host on `machine-f-diiid.toml`, under `debug` AND under `cuda`, at more than one thread count. M-79's own closing row is the template: cpu, cuda and debug agreeing to every printed digit. **MET at `OMP = MKL = 1`** → M-130 — `psi_ax = 3.759851e-01`, X-point `( 1.200929, −0.999491 )`, constraint `−6.466e-12`, 2 steps. The thread-count half is stage 2 and is still open, and **`mhd-rectangle` still takes 9 Newton steps under a device against 5 on the host, with no MMU trace — unexplained, and nobody has walked it** |
| **5. Anything about speed** | **not on this machine.** §3.3 |

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
  comment table is re-taken against M-130 — soloviev 1→1, mhd-rectangle 5→9,
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
  both. What survives is `mhd-rectangle`'s 5→9, which the syncs do **not** move
  and which nobody has walked.
* ~~**`tests/performance/NpcThreadScaling.cpp:544-570` states that "a whole NPC
  solve with an `mfem::Device` configured for CUDA does not compute the right
  answer".**~~ **CORRECTED** — it points at M-130.
* **Upstream's crossover table is quads on their fixture.** Whether MEQ's
  triangles at 450 to 6707 elements land anywhere near it is unknown and is the
  single measurement most worth taking first, on a machine that can produce an
  admissible FP64 number.
