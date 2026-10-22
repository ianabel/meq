# The device, threading contracts, and single precision

**`CLAUDE.md`'s device half**: every defect found while making MEQ run under an
`mfem::Device`, the instruments that found them, the residency work, cuDSS, the
threaded-assembly promise MEQ makes to MFEM, the integrator inventory and the
batched routes, and the single-precision build. `DEVICE-PLAN.md` is the plan;
these are the findings.

**The standing rules are `CLAUDE.md`'s.** Measurement tables live in
`MEASUREMENTS.md` under stable `M-nn` anchors, which are not renumbered. A
hypothesis that was measured and falsified stays; the history of the file does
not — present tense, no "used to", no dated "now fixed". What was decided
against is `CLAUDE_DECLINED.md`.

## Where it stands

* **`--device cuda` reproduces `--device cpu` to every printed digit** on the
  production cases. Of the 29 shipped examples M-238 ran, 27 agreed and the
  two exceptions are both resolved (M-239, M-299; below). The 29 are that
  sweep's set; `examples/` holds 44 `.toml` files now, and the sweep has not
  been re-taken over them. Three of upstream's four element-local groups run on
  the device — local factorisation, the NPC traversal, the trace assembly — and
  MEQ's integrand does too. The solve leg is at parity with the host; what a
  device run loses, 1.07–1.13× on the wall, is host threading it has to give up
  (→ **[M-223](MEASUREMENTS.md#m-223)**).
* **The driver's device defaults**, each applied only where the file said
  nothing: `LocalFactorMode` and `TraceAssemblyMode` batched, `AssemblyMode`
  **`Serial`** (bit-reproducibility, below), and cuDSS as the trace solver under
  CUDA. Threaded assembly is refused above one thread, and
  `TraceSolver = "cudss"` is refused without a Device.
* **Residency is built and gated on the data.** The iterate's readers, the axis
  search and the bordered columns have device routes that open by themselves
  once their input is device-valid. Over a `machine-f-diiid` run the bus carries
  508 MB up and 78.6 MB down (M-240).
* **Single precision builds and runs at parity with double** (83 of 84 in each
  build, the same deliberate red, single on MKL 2026.2 included, M-303), and its accuracy is qualified **per problem**.
* **This machine's GPU is an RTX 2070 SUPER** — evidence about correctness,
  never a timing (see the last section).

## The instruments

Every device defect below was found by one of these instruments, and most were
invisible to the others. Reach for them in this order.

1. **`--device cuda` reproducing `--device cpu` to every printed digit.** This
   is how M-79 and M-130 were root-caused, and every ordering decision in the
   device plan exists to protect it: residency is built in **double**,
   `AssemblyMode` is `Serial` under a Device, and `DeviceParity` pins an exact
   CUDA arm.
2. **`mfem::Device( "debug" )`** has device memory semantics with host
   arithmetic and `mprotect`s the host page. So **a raw host read of a
   device-valid buffer becomes a named fault with a backtrace** rather than a
   wrong number. It costs nothing (`NpcThreadScaling --device debug` takes the
   string). Each fix moves the fault **forward** to the next site, which is the
   only thing that distinguishes a fix from a coincidence. Reach for it before a
   debugger. **Its blind spot is structural**: it catches an access to the
   **wrong location** and cannot see honest flags over a **stale copy**. Three
   defects below (M-187, M-206, M-226) ran clean under it, one through a
   1170-second run with zero faults.
3. **A stage-by-stage checksum, diffed between `cpu` and `cuda`.** The first
   stage that disagrees names the defect. `meq::stageChecksum()`
   (`src/meq/StageChecksum.{hpp,cpp}`) reports each buffer's **validity flags
   before it reads**, so `H- D+` names a stage whose host copy is not live,
   whether or not its value happens to be wrong. It prints only under
   `MEQ_STAGE_CHECKSUM`, an environment variable rather than a key, so a walk
   can be run on a shipped `examples/*.toml` without editing it.
4. **A count of transfers, not of calls or seconds.** `meq::countedHostRead()`
   asks `Memory::DeviceIsValid()` before the read, and so separates a funnel
   call from a data movement. A count of *calls* points at the wrong sites:
   MEQ's busiest funnel, the bordered `( N + 4 )²` row dot, makes 240 reads a
   step and transfers **zero**. `meq::countedDeviceRead()` is the upload mirror
   (`HostIsValid() && !DeviceIsValid()` before the `Read()`). `meq::PhaseScope`
   attributes every crossing to a `DEVICE-PLAN.md` §3.3 residency phase, A to E,
   because the same transfer is the design in phase B and a defect in phase C.
   **Read the two directions together**: uploads in B, zero in C and one
   download in D is residency; zero everywhere is an empty device. Three limits
   of this instrument are measured:
   * A per-site table attributes a transfer to whoever asks **first**, not to
     whoever needs it (M-232).
   * It counts MEQ's funnels and nothing beyond them. It read 0.888 MB on a run
     that `nsys` puts near 950 MB.
   * An instrument that sees only one direction of the bus cannot measure a
     plan whose first phase is an upload (M-226).
5. **A `gdb` breakpoint on `mfem::CuMemcpyDtoH`**, for crossings inside
   upstream's code that the counters cannot see (M-237). **On `mprotect`, with
   the faulting page as the condition**, when `debug` faults somewhere
   unhelpful. Use a `RelWithDebInfo` build: without debug info, gdb collapses
   inlined template frames onto the wrong function.
6. **`OMP_NUM_THREADS` alone**, which separates a thread-safety defect from a
   memory one (M-188).
7. **`meq::ProbedSolver`**, which sits where `mfem::NewtonSolver` keeps
   `c = J⁻¹r` to itself. It tells a wrong Jacobian from a wrong residual: on
   `mhd-rectangle` it read **the same right-hand side to twelve digits returning
   an increment 36× too small** (M-210).

**`DeviceParity` runs three arms, `cpu`, `debug` and `cuda`**
→ **[M-228](MEASUREMENTS.md#m-228)**.
* The arms are not repeats of each other: `debug` and `cuda` catch disjoint
  classes, as item 2 above explains.
* `debug` is the expensive arm: 329 s of the test's roughly 360 s is
  `machine-f-diiid` under it. The card arm costs 35 s.
* It skips loudly where MFEM's own `Device configuration:` line reports no
  device.
* It runs a copy of `diverted-tokamak` with `PlasmaSupportSweeps` stripped, so
  that a **moving** support stays covered (M-239). *A device test's cases have to
  exercise every way host state can move, not only every leg.*
* **The mirror is `RotatingMirrorDevice`'s, not one of its cases**: the box,
  two labels at points and subtracted filaments, `cpu` against `cuda` on the
  exact `ψ_h` and `q_h` within 1e-12 of their largest values. A tolerance and
  not printed digits, because the card differs from itself between takes
  → **[M-280](MEASUREMENTS.md#m-280)**.

**THREE PATTERNS RECUR, AND RECOGNISING THEM IS MOST OF THE WORK.**

* **A partial fix looks like a new defect**, not like progress. The chains below
  are links of one fault, so fixing one link leaves a symptom standing, often a
  worse-looking one.
* ***A cure measured after changing four things attributes nothing*** (M-210).
  ***A reproduction is not a diagnosis.***
* **A bordered row certifies consistency and not correctness.** It is satisfied
  by the state it was differenced at, so a wrong solve can close its constraint
  at 1e-12 or 1e-16 with every printed line healthy (M-206, M-226).

## Where the device code lives

**Seven files in `src/meq` are compiled by `nvcc`** — `PointLocatorDevice.cpp`,
`SourceViewDevice.cpp`, `DeviceNPCDevice.cpp`, `DeviceQuadratureDevice.cpp`,
`DeviceRowsDevice.cpp`, `AxisSearchDevice.cpp` and `ResidentIterateDevice.cpp`.
They hold every `mfem::forall` launch in the tree, because `mfem::forall`
`static_assert`s under `g++`. CMake reads the architecture out of MFEM's
`config.mk`. Without CUDA each compiles to a host stub. The non-device halves
are these:

* **`meq::DeviceBasis`** evaluates an L2 grid function at an *arbitrary*
  reference point on a straight triangle, which MFEM cannot do on a device at
  all.
* **`meq::AxisSearch`** finds the magnetic axis by rooting elements on the
  device. The host keeps the selection, the conductor filter and the merge,
  which are the half carrying the physics.
* **`meq::ResidentIterate`** holds the iterate's other readers, each
  reproducing its host loop to the bit.
* **`meq::SourceView`** is the POD form of the source.
* **`meq::DeviceRows`** is the bordered rows' mirror.
* **`meq::DeviceConstants`** holds the per-mesh constants.
* **`Crossings`** and **`StageChecksum`** are the two instruments above, and
  neither is part of the solve.

## Memory-validity defects, by class

MFEM keeps a host copy and a device copy of each buffer, with validity flags
that say which copy is current. Every defect in this section is one way of
making those flags lie, or of being honest about the flags while reading the
wrong copy.

### An alias leaves its base's flags stale → [M-79](MEASUREMENTS.md#m-79), [M-130](MEASUREMENTS.md#m-130)

**The symptom**: a whole NPC solve mis-solved at `OMP_NUM_THREADS=1`.
* `example5` reported **0/0 Newton iterations on a case that needs four**, with
  the potential identically zero.
* All three trace solvers agreed to four figures, which said one common fault
  rather than three. **cuDSS was a bystander**: UMFPack and PARDISO never touch
  the device and failed identically.

**The mechanism**: `GetAliasDevicePtr()` ends in `AliasProtect()` on the
**base's** host range and never touches the base's flags. `Memory::SyncAlias`'s
comment says so and says what is owed.
* The base still reads `VALID_HOST`. The next `ReadWrite_` computes
  `copy = !( flags & VALID_DEVICE )` = true, and memcpys host-to-device out of
  the page it just protected.
* Under CUDA nothing is protected, so it silently copies the **stale host copy**
  over a device buffer the solve had written correctly.

**MEQ's sites, in the order the debug fault moved through them:**
* **The constructor** makes three long-lived `MakeRef` aliases (`darcyFlux`,
  `potentialGf`, `traceGf`) onto the blocks of one `BlockVector solution`. They
  are in the constructor, *not* in the re-entered `buildForms()` where a reader
  looks. `solve()` calls `solution.SyncToBlocks()` and a `SyncMemory()` per grid
  function before the first host read. Both calls are inert with no Device.
* **`prepare()`** seeds the iterate by writing *through* `potentialGf` and
  `traceGf`. That marks the alias device-valid and leaves the base's flags
  unchanged, so `formSystem()`'s `darcySolution` and `traceX` inherit the lie.
  `formSystem()` makes the `SyncAliasMemory` calls.
* **`solveWithNormalisation()`** reads `GetEssentialTrueDofs()` through
  `Array<int>::operator[]`, which is as raw an accessor as `Vector`'s. Under
  CUDA the border columns are then differenced against a state that is not the
  iterate, and the run reports *"the bordered Jacobian is singular"*. **One
  fault, two faces, and the CUDA face names a border rather than a memory.**
* **`assembleExteriorColumns()`** reuses one `mfem::Vector load` across the
  modes. It zeroes it device-aware (`load = 0.0` is `Write()` then a `forall`)
  and hands it to an `mfem::LinearForm` through `Update()`, which *aliases* it.
  It assembles on the host into the alias, then reads the **base** with
  `load.HostRead()`, which still holds the zeros.
  * **Mode 0 survives**: a fresh `Vector` has `UseDevice()` false, so its first
    `= 0.0` is a host loop. Every later mode reads zero.
  * So the bordered Jacobian's exterior block was nine tenths missing, while the
    exterior coupling stayed exact in the **residual**.
  * That is the M-187 signature: the right answer at 4× the work (**31 / 16 / 4**
    Newton iterations per support sweep against the host's **7 / 3 / 2**). It is
    linear convergence at a fixed factor of about ten, which is what an inexact
    Jacobian looks like.
  * The repair is `form.SyncAliasMemory( load )` before the read.
  * Neither `debug` nor any assertion could see it. The stage checksum found it,
    four stages past where anyone was looking.

**Upstream's half is `bf27f5a928`.** `DarcyNPCOperator::Mult` and
`DarcyNPCSolver::Mult` wrote a `BlockVector`'s blocks without syncing. In the
second, **order is a constraint**: `Neg()` acts on `xb`'s own `Memory`, so the
blocks are carried up to `xb` *before* it and `xb` up to `x` *after* it.
**Upstream's fix alone repairs the count (0/0 becomes 2/2) and leaves `max|ψ_h|`
identically zero.** A library-side test cannot see a caller's aliases.

**THE LESSON IS ABOUT THE FUNNELS.** Every sync in MEQ's raw-read sweep is a
`HostRead()` **on the object about to be read**. That is right when that object
is the one that was written, and useless when the write went through an alias.
***A sync on the base cannot see data that is still on the alias.***

### A view or a copy carries its own flags → [M-205](MEASUREMENTS.md#m-205), [M-238](MEASUREMENTS.md#m-238), [M-243](MEASUREMENTS.md#m-243)

* **A copy construction has no accessor in it.** `mfem::Vector v2( v1 )` is a
  raw `memcpy` of `v1`'s host pointer. A sweep asking *is this read funnelled*
  of every `operator()`, `GetData()` and `HostRead()` cannot see it, because a
  copy construction is a **declaration**. Three sites in
  `solveWithNormalisation()` copy the bordered `unknown`, and all three now
  funnel first. To find the rest,
  `grep -nE 'mfem::(Vector|Array<[^>]*>) (const )?[a-zA-Z]+\( *[a-zA-Z]'`.
  **The sweep's question was right and its surface was short.**
* ***Funnelling the source of a copy does not funnel the copy.***
  `Vector::Vector( Vector const & )` propagates `UseDevice()` and copies
  device-to-device. So `coldState( unknown )` is born device-valid, with its host
  page protected. This was unreachable while another funnel ran first on every
  path, and **deleting** that funnel exposed it.
* **A view is a copy of the `Memory` object, not an alias.** In
  `DarcyHybridization::NonlinearGradientFinish()`, `Df_data += Df_lin_data` read
  stale data. The caller-integrand residual reads `Df_lin_data` on the device
  through `DenseTensor::NewMemoryAndSize()`, and that device read is recorded in
  the view's flags alone. Reproducing it needs a residual *between* two
  gradients, which is a Newton loop. ***A sync on the object cannot see a read
  taken through a copy of its bookkeeping.***
* **`BlockVector::SyncFromBlocks()` makes each block match its BASE.** A kernel
  writing the blocks of a vector whose base says *host* therefore has its output
  downloaded on the spot.
  * Upstream's `DarcyNPCSolver::ArrayMult()` did this to every backsolved
    column: fourteen columns a step, 140 MB over `machine-f-diiid`.
  * MEQ's residual operator did it to its own output, while its comment said the
    call *carries the blocks up*.
  * The repair is `Write()` on the base before the blocks are written, plus
    `UseDevice( true )` on each block by name where `Update()` kept blocks made
    before the base was device-capable.
* **A view's validity flags belong to the view.** Uploading `darcyRhs` added
  11 MB and moved the count not at all. `ResidualRowsBatched()` reads the base's
  block **views**, and `darcyRhs.Update()` re-makes them host-only at every
  `prepare()` (16 calls of 28). *A fix that does not move its number tests the
  attribution only if it reached the object the reader holds.*

### A buffer refreshed in place is device-valid from its own previous use → [M-211](MEASUREMENTS.md#m-211), [M-226](MEASUREMENTS.md#m-226), [M-187](MEASUREMENTS.md#m-187), [M-207](MEASUREMENTS.md#m-207), [M-210](MEASUREMENTS.md#m-210)

The mechanism is the same at every site. A fresh `mfem::Vector` or
`mfem::Array` has `UseDevice()` false, so its first raw fill is safe. After its
first `Read()` it is device-valid. `SetSize()` at an unchanged size reallocates
nothing. So the next raw fill writes the host copy while MFEM still believes the
device copy is current, and the following `Read()` copies nothing. ***A fill
that is safe in a constructor is not safe in a refresh*** — and phase B is
exactly the work of turning the first into the second, so expect this defect
wherever work is hoisted out of a loop onto a device.

* **`meq::DeviceRows`' mirror**, re-uploaded once per Newton step: the first
  call survives and every later one writes into a protected page. The fix is
  `HostWrite()` before each fill, and `theMirrorSurvivesBeingRefreshed` uploads
  five times.
* **`meq::DeviceConstants`.**
  * Hoisting the five per-mesh constants out of `meq::DeviceNPCOperator`'s
    rebuild took **80 uploads and 1.22 MB to 6 and 0.113 MB**.
  * It shipped a **neighbouring answer**, `psi_ax = 3.759845e-01` against
    `3.759851e-01`, with the bordered constraint still satisfied at 1e-12.
  * The counter that found it, `countedDeviceRead()`, was written for it: its
    per-site table read **`2 device reads, 1 upload`** on the flag buffer, a
    fill that happened beside a transfer that did not.
* **Upstream's `DarcyNPCSolver` reused `b_tr`.**
  * It was written through an alias in `NPCReduce()`'s host route and never
    carried back, while the previous call's `b_tr -= r_tr` had left it
    device-valid.
  * This was M-187's 4× work, and the M-207 credit to MEQ's own
    `xLoc.SyncFromBlocks()` is falsified. Upstream added that exact sync and the
    count did not move, and MEQ's 2 × 2 says neither of its two changes holds it.
  * MEQ's report named `x_loc` by analogy with its own defect. Nothing in
    `DarcyHybridization` reads `x_loc` as a whole, because both entry points take
    its blocks, and that was checkable by reading two functions.
  * What separated the right member from the wrong one was upstream running the
    same Newton two ways in one process.
  * The `xLoc` sync stays because the flags should be honest, and its comment
    says it is not load-bearing.
  * **`mhd-rectangle` takes `mfem::DarcyNPCOperator` directly, and upstream
    closed that half at `ce744b84de`**: 5 Newton steps under `debug` against 5 on
    the host.

### An upload of host state is a snapshot → [M-239](MEASUREMENTS.md#m-239)

The device integrand's element flags were uploaded at the NPC operator's
rebuild. Under a **moving** plasma support, `refreshPlasmaComponent()` changes
the mask at every residual. The host integrand followed it and the device one
did not, so `diverted-tokamak` under any Device converged to a different
equilibrium.

**The signature tells semantics from arithmetic**: seven host perturbations
agreed exactly, and four device arms agreed exactly with each other but not with
the host. Read that before any code. The same entry records that this example
landed on M-132's other branch, 2.6% from its `freegs4e` reference, until
`PlasmaSupportSweeps = 4` put it on the benchmarked one.

### Host-filled vectors handed to SUNDIALS → [M-206](MEASUREMENTS.md#m-206)

**The symptom**: `examples/q-driven.toml` was wrong by a factor of two under any
configured Device, with nothing printed in red.
* `psi_ax` came out `1.723374e-01` against the host's `3.457380e-01`.
* KINSOL reported *"converged in 0 iterations … final norm 0"*.
* The inner bordered solve took the same 4 Newton steps on both sides and closed
  its constraint to 1e-16.

**`mfem::SundialsNVector` binds both of SUNDIALS' pointers into one `Memory`**
(`linalg/sundials.cpp:465`), and KINSOL's CUDA kernels read the device array.
Five seams were each hidden behind the one in front of it:
* `FixedPointResidual::Mult` filling `y` rawly.
* `::GetGradient`'s `plus( j ) += step` after a device-aware `plus = x`. This
  gave a zero finite-difference Jacobian, and `-nan`.
* `DenseJacobianSolver::Mult`'s `x = 0.0` followed by a raw fill.
* `::sameAs`'s raw comparison of two iterates. **This seam has no symptom of its
  own**: it decides whether the cached Jacobian is reused, so the iteration
  contracts linearly with every norm looking healthy (30 iterations against 12).
  *A raw read on a control decision costs convergence rate rather than
  correctness.*
* The starting guess and the answer in `solveForToroidalField`.

Fixing only the first seam replaced the wrong answer with 730,196 spinning NaN
iterations.

**`HostWrite()` is not enough on its own.** It makes the flags honest and moves
nothing. What is also needed is a `Read()` *after* the fill, to do the
host-to-device copy into the array the consumer actually reads. *A sync fixes a
lie about where the data is; it does not fix the data being in the wrong place.*

**It is not the NPC fault.** Forcing `makeDeviceNPCOperator()` to return null
removes that fault and leaves this answer wrong. With all five seams fixed, CUDA
matches the host:
* 5 KINSOL iterations and 36 map evaluations, on both.
* Conditioning `6.440e-03`, the `g²` coefficients and `psi_ax = 3.457383e-01`
  to every printed digit.
* Final max norm `3.28e-10` against the host's `3.54e-10`, both under the
  `1e-9` threshold — the last Newton step squares a residual whose rounding
  differs.

The outer Newton's own Jacobian-vector product for KINSOL's line search
(`OuterKINSolver::jacobianTimes`) is a sixth crossing of the same kind, and is
written to the same pattern: `countedHostRead` the direction, `HostWrite` the
product, `Read()` it back up.

`q-driven` is the only example that re-enters `solve()` from an outer Newton.
*A defect in re-entry needs a case that re-enters.*

### Element blocks read in a host loop → [M-214](MEASUREMENTS.md#m-214)

`LocalFactorMode::Batched` left `Ct_data` device-only, and
`DarcyHybridization::ReconstructTotalFlux()` reads it in a host `GOMP_parallel`
loop. One key isolated it: `Serial` passed, and `Batched` aborted in
post-processing with the solve already correct. Upstream's `d720f02d62` gives
the two batched NPC traversal kernels the `HostRead()` that every other kernel
reading those blocks already has.

**MEQ's claim that CUDA would ship a wrong `ψ*` is falsified.** `_psistar.gf` on
`cpu` against `cuda` agrees to 6.94e-15. The reasons:
* A device `Read()` only adds `VALID_DEVICE` and never clears `VALID_HOST`.
* No kernel writes `Ct_data`.
* Only the debug backend protects the page.

*An argument that a defect COULD corrupt something is not a measurement that it
did.*

### Pure waste: device-aware zeroing of host scratch → [M-200](MEASUREMENTS.md#m-200), [M-202](MEASUREMENTS.md#m-202)

`mfem::Vector::operator=( 0.0 )` is device-aware, so zeroing host-only scratch
that way gives it a device copy, and the next host read drags it back. Two sites
did this:
* `assembleExteriorColumns()`'s `load`.
* The X-point's `xFluxField` / `xPotentialField`.

Four `UseDevice( false )` calls took `machine-f-diiid` from 36 transfers and
14.9 MB to 3 and 0.888 MB. No answer changed, and no leg profile could see it,
because the waste happens once per support sweep.

## A Device and a host OpenMP region are mutually exclusive → [M-188](MEASUREMENTS.md#m-188)

Under any configured `mfem::Device`, MFEM registers buffers lazily into
`general/mem_manager.cpp`'s `maps->memories` and `maps->aliases`. The lazy
registrations are `MakeAlias()`, `SetDeviceMemoryType()`, and the device arms of
`Read`/`Write`/`ReadWrite`. Those maps are bare `std::unordered_map`s with **no
mutex, critical section or atomic anywhere in the file**, so any thread team that
reaches one corrupts it.

**The instrument is a thread count.** `machine-f-diiid` exits 0 at
`OMP_NUM_THREADS=1` and aborts at 8, under `cuda` and `debug` alike. **With every
batched key `"serial"`, `debug` still segfaults at 8**, in
`MemoryManager::Insert()` beneath `meq::CriticalPointFinder::sweep()`'s own
`GOMP_parallel` frame:
* `GridFunction::GetVectorValue()` allocates an `Array<int>` per thread, and
  under a Device that allocation registers.
* `cuda` survives the same loop by luck: its host memory type is `HOST`, which
  leaves owners unregistered, where `debug`'s `HOST_DEBUG` registers them.

**MEQ's half is one condition, stated once.** `meq::hostLoopsAreSafe()` in
`src/meq/Threading.hpp` copies `DarcyHybridization::ThreadHostLoops()`'s
`if ( Device::Allows( Backend::DEVICE_MASK ) ) { return false; }` and is carried
into MEQ's eight element loops by `MEQ_OMP()`'s `if` clause. **`--device cpu` is
not a device by that test**, so every host measurement stands and the host path
is bit-identical. **MFEM's half MEQ can only refuse**: with MEQ's loops gated,
`AssemblyMode = "threaded"` under CUDA at 8 threads still aborts. So the driver
exits 1 naming both the key and the device, and
`theDriverRefusesThreadedAssemblyOnADevice` asserts the **pairing**: refused at
two threads, honoured at one, host path untouched.

**How this survived M-79**: M-129 measured these aborts, reproduced them twice,
and wrote them down as a property of the device path. *A result that is stable
and reproducible is still a bug until something explains it.*

## Residency

**Residency is a property of a CYCLE and not of a leg.** Every host leg must
move with the others, or the one host leg turns every other leg's kernel into an
upload. Two measurements say so:

* **cuDSS bought no residency** (M-227, M-229). The bordered columns stayed
  **168 of 168 `H+ D-`**. They are declared `mfem::Vector v( n )`, and **a
  freshly constructed `mfem::Vector` has `UseDevice()` false**. So item 4's gate
  was a predicate the data could not satisfy under any trace solver. Setting
  `UseDevice( true )` on them cost **17× the bus traffic**: the crossing MOVED,
  because `unknown.Add( -damping, y )` makes the iterate device-valid and the
  next step's axis search drags it back.
* ***A leg ported in isolation does not remove a transfer, it RENAMES one***
  (M-232). Forced open, the axis constraint's 12.5 MB download reappears at
  `refreshXPoint()`, and with that ported too it reappears at `limiterValue()`.
  The count falls only when the **last** host reader goes.

**So every route is gated on the DATA** (`DeviceIsValid() && !HostIsValid()`),
the only form that is right both before and after residency arrives. A route
stays inert until its input lives on the device, and then opens by itself.
`meq::DeviceRows` has had this arrangement since M-208.

### The axis search on the device → [M-230](MEASUREMENTS.md#m-230), [M-231](MEASUREMENTS.md#m-231), [M-232](MEASUREMENTS.md#m-232), [M-233](MEASUREMENTS.md#m-233), [M-234](MEASUREMENTS.md#m-234)

`DEVICE-PLAN.md` §3.4 called the mesh search the hardest item, and was wrong
about why: it is *add a kernel*, once a kernel can evaluate `q_h`.
* **`FiniteElement::CalcShape` is a host virtual taking an
  `mfem::IntegrationPoint`**, so no device code could evaluate a field at an
  arbitrary point, and an element-local Newton is nothing but such evaluations.
* `PointLocatorDevice.cpp` had named the same blocker from the other side: *"IT
  LOCATES, IT DOES NOT INTERPOLATE"*.

* **`meq::DeviceBasis` borrows nothing from MFEM's internals.**
  * The host asks each element for its shape functions at its own nodes, and
    solves for the constant matrix that re-expresses them in barycentric
    monomials. That is phase B in its purest form.
  * A fixed quadrature rule could not serve: *a rule is a set of points; a
    search needs a function*.
  * **The ceiling is degree 4, a measured conditioning bound.** Outside the
    element, where the Newton deliberately spends time, the error grows twenty
    times per degree and reaches 4.8e-07 at degree 6.
* **The acceptance is the candidate list, not the answer.**
  `meq::axisCandidatesOnDevice()` agrees with `CriticalPointFinder::sweep()` to
  **2.2e-16 in both directions**. *A test on the final answer passes a port that
  found a different set and picked the same one*. `machine-f-diiid`'s `psi_ax`
  is the same whether or not the device list is merged, and the merge is the
  half that was nearly left out.
* **The search covers a neighbourhood.** The device roots every element in a
  sweep, because on a card M-100's 40× economy for a seeded search inverts — but
  a full sweep costs 36× the host's warm ring search, so the device roots a block around the previous axis.
  It finds the block through `meq::PointLocator`'s CSR bin grid, exposed as
  `PointLocator::grid()`: one grid, not two uploads. A warm attempt that finds
  nothing widens to the whole mesh once. *A neighbourhood that silently missed
  the axis would hand the selection a competition its winner was not entered
  in.*
* **The kernel was 65× off its own arithmetic.**
  * Both reductions were `mfem::forall( 1, … )`, so the flux scale walked 1.23
    million dofs on one CUDA core.
  * The loop bounds were runtime values, so the per-thread arrays went to local
    memory: **608 bytes of stack frame with zero register spills**. *A stack
    frame with no spills is not register pressure; it is an array whose bounds
    the compiler cannot see.*
  * Templating on the degree takes the frame to 0 bytes and the sweep from
    **297 ms to 11.9 ms** at 102,840 elements, with the lists still
    bit-identical.
  * It is not `mfem::reduce`, which takes a host reference and would put the
    kernel's one device scalar back on the bus.

**Measured end to end** → **[M-235](MEASUREMENTS.md#m-235)**: the device wins
the production bordered case in both formats.
* `mf`'s solve reads 0.81× and 0.85× of the host, and the fastest of sixteen arms
  is single + CUDA at 2.445 s against the double host's 3.727 s.
* The device loses the big fixed-boundary cases, and **not on the FP64 rate**.
  The axis route was correctly off (`0 taken, 15 declined`), so the host sweep
  ran with `hostLoopsAreSafe()` false: 1.00 cores, 18–20% of the wall, with the
  11.9 ms kernel built and idle beside it.

### The iterate's other readers, and the bordered columns → [M-237](MEASUREMENTS.md#m-237), [M-238](MEASUREMENTS.md#m-238)

`src/meq/ResidentIterate.hpp` gives device routes behind the same gate to:
* the five plasma element loops;
* the transmission gather;
* the `psi_bnd` row;
* the three copy constructions.

Forced open, 126 of 126 evaluations equal their host loops.

* ***A second expression of the physics is bit-identical only where it is the
  same expression.*** `hermite10Prime` factored where `Profiles.cpp` expands,
  costing an ulp in every derivative route. `SourceViewTests` could not see it,
  because its derivative envelope is four ulps of `psi` wide.
* **Removing the readers does not open the gate.** Unforced, the routes read
  `0 taken` until the bordered columns are device-capable, and then every route
  opens by itself. With the columns resident, device-to-host traffic over a
  `machine-f-diiid` run falls from 259 MB to 101 MB and host-to-device from
  692 MB to 543 MB.
* **The dense fill's device contraction is bit-identical to `rowDot()`** only
  because `psi_bnd`'s `-limiterValue()` is a flag on the finished sum and not a
  negated weight. *`-( +0 )` is `-0`, and a sum of negated terms from `+0` is
  `+0`.*
* **Of the 29 shipped examples M-238 ran, two differed under CUDA, and both
  are resolved.**
  * `diverted-tokamak` is the stale-snapshot defect above (M-239), now fixed.
  * `rotating-normalised` located a different axis on the host: its axis sits
    on the mesh line `z = 0`, and the host's warm ring walk followed one twin
    of it where the device and the sweep took the other. The ring walk now
    merges twins by the sweep's rule, and host and `--device debug` both read
    1.039297e-01 (→ **[M-299](MEASUREMENTS.md#m-299)**).

### The first evaluation of a solve → [M-240](MEASUREMENTS.md#m-240)

Uploading the iterate before its first reader made phase C *worse* on its own,
for this chain of reasons:
* A `Type = "ramp"` guess leaves the flux block zero, so every element reports a
  root.
* The device route refused that list as over its ceiling.
* The host fallback then downloaded the iterate it had just been sent.

`perf` put the 22 s in `CriticalPointFinder::mergeCandidates()`, a quadratic walk
of the kept list that the host paid too. The fixes:
* **The merge is now grid-indexed and returns the walk's list to the bit**, with
  the walk itself as the test's oracle.
* The ceiling is gone, the list being bounded by the mesh.
* The nodal fallback is two exact reductions on the device, because
  `mfem::ArgMaxReducer` gives a tie to whichever element joined *last*.

`ff3`'s axis leg reads 0.40 s against 22.1 s. *A cost attributed to a route is
the route's only if the other route does not pay it.*

Four more of MEQ's crossings went at the same time:
* Phase D downloaded the answer twice through its aliases.
* The output stage restored mesh nodes with a device-aware `=`.
* `fluxGf` was negated on the device.
* `rLoc` arrived host-valid.

`GradShafranovSolver::downloadSolution()` is phase D.

## cuDSS

**cuDSS 0.8.0 is nondeterministic run to run, with no usable switch at MEQ's
`N + 4` columns.** `CUDSS_CONFIG_DETERMINISTIC_MODE` defaults to 0. With the
same matrix, the same right-hand sides and a fresh handle, 4 of 4 repeats
differ: 6.7e-16, 1.8e-15 and 3.6e-15 at `K = 1, 4, 16`.
→ **[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md)**, *cuDSS determinism:
SetDeterministic and per-column Mult*.

**MEQ is clear of the two wrong-answer traps, and was checked**:
* `CUDSS_MTYPE_SYMMETRIC`, under which the deterministic switch returns an
  exactly zero solution.
* The setters' order dependence. MEQ builds `NONSYMMETRIC` / `FULL` and calls
  `SetMatrixSymType()` before `SetMatrixViewType()` at both call sites
  (`GradShafranov.cpp`, `TraceSolverScaling.cpp`).

**cuDSS is the trace solver under CUDA** → **[M-229](MEASUREMENTS.md#m-229)**,
whenever `Backend::CUDA_MASK` is allowed and the file named no trace solver.
**It is `CUDA_MASK` and not `DEVICE_MASK`**, because `DEVICE_MASK` includes
`DEBUG_DEVICE`, which has no GPU for a vendor library to run on. The batched
keys beside it rightly use `DEVICE_MASK`, being about memory and routes rather
than a vendor library.

**The parity instrument survived cuDSS, which §3.5 predicted it would not.**
* The bordered constraint is genuinely nondeterministic: `-6.465e-12`, then
  `-6.466e-12`.
* Over eight repeats, every field `DeviceParity` compares is identical: 2 Newton
  steps, `psi_ax = 3.759851e-01`, `psi_bnd = 7.081394e-02`, and the X-point at
  `( 1.200929, -0.999491 )`.
* The driver prints seven figures of a quantity perturbed at the fifteenth. So
  the exact assertion stays exact, resting on printed precision now rather than
  on determinism. **If it ever flickers, widen it to the measured spread**, not
  to a convenient tolerance.

**`TraceSolver = "cudss"` without a Device is refused, and the reason is
measured.**
* **The immediate failure**: without a Device, `CuDSSSolver` reads host pointers
  through device-aware accessors and aborts with a raw
  `CUDA error … cudaMemcpyDeviceToDevice … invalid argument`. That message names
  nothing about the key that caused it.
* **The trade**: a device solver pays only if the data stays on the device. M-126
  puts the trace solve at 20.1% of a threaded step. A config-file cuDSS on a host
  run would route every `Vector` through the Device and pay a round trip per
  Newton iteration for the one part that moved. That is upstream's
  "group 4 alone" case, which their own measurement calls plausibly slower than
  staying on the host.
* **The measurement moved under it**, which is why the gate is *refuse without a
  device* rather than *refuse always*. M-173 read `--device cuda` at 2.14× slower,
  but that arm ran `AssemblyMode::Batched`, which costs 1.74× on its own
  ([M-185](MEASUREMENTS.md#m-185)). With today's defaults the solve leg is at
  parity, and on a resident path cuDSS is the only trace solver that does not
  force a round trip.
* **The library route stays open**: `TraceSolverScaling` constructs the Device
  first, and the `cuDSSTraceSolver` ctest is how cuDSS's agreement with UMFPACK
  is checked. `theDriverRefusesASolverItCannotHonour` keeps the config-file route
  closed.

*Exposing a knob makes reachable, by somebody who has not read the library's
documentation, every precondition that documentation records.* The Device
requirement had been written in `GradShafranov.hpp` all along.

## `AssemblyMode::Threaded` is a PROMISE MEQ makes, and it is load bearing

**MFEM refuses threaded assembly on any problem whose threaded element loop
would evaluate an integrator.** `MultNL()` aborts, naming `SetAssemblyMode()`.
The library can see which integrator objects it will call, and nothing about
whether they are reentrant. MFEM's convention is `#ifndef MFEM_THREAD_SAFE`
around the scratch, which is a build switch rather than a property of the class,
and 67 of the 99 scratch-carrying classes do not use it.
`DarcyHybridization::SetIntegratorsThreadSafe()` is where a caller who has read
their integrators says so. It is off by default.

**MEQ calls it in `buildForms()`, beside `SetAssemblyMode()`, and the audit is
the justification.** `DarcyForm::CanThreadAssembly()` gates the assembly loops on
the same flag, so two sets of integrators are covered:
* `meq::SourceIntegrator`, whose per-point scratch is already under
  `#ifndef MFEM_THREAD_SAFE`.
* `VectorMassIntegrator` and `VectorDivergenceIntegrator`, which upstream guarded.

The `HDGExtensionIntegrator` and `NormalTraceJumpIntegrator` scratch are guarded
upstream too (`CLAUDE_MFEM.md`, *`../../mfem-hdg-dev` is receive-only*).

**It is not decorative, and that is falsified rather than argued.** Comment out
the one line and `machine-f-diiid` exits **2** with *"AssemblyMode::Threaded is
unsafe for this problem"*, on both the Newton and the Picard-then-Newton paths.

## Which integrators need device offload, and the batched routes

**The inventory, taken from the installed headers.** Sites are named by function
rather than line number, because line numbers in this table went stale by
hundreds of lines. `grep` for the integrator.

| integrator | where MEQ installs it | device path |
|---|---|---|
| `mfem::VectorMassIntegrator( radius )` | flux mass, domain — `buildForms()` | `AssemblePA`, `AssembleMF`, `AssembleDiagonalPA`, `AddMultPA`, and separately `HDGElementMassBatched()` under `AssemblyMode::Batched` |
| `mfem::VectorDivergenceIntegrator` | flux divergence, domain — `buildForms()` | `AssemblePA` |
| `mfem::TransposeIntegrator( mfem::DGNormalTraceIntegrator )` | flux divergence, interior + boundary faces — `buildForms()` | irrelevant: under hybridization these are never assembled. They are MARKERS (`CLAUDE_MFEM.md`, *A trap in that table*) |
| `mfem::NormalTraceJumpIntegrator` | the hybridization constraint, handed to `EnableHybridization()` | `AssembleEAInteriorFaces` **only**: no `AddMultPA`, nothing for boundary faces |
| **`mfem::HDGDiffusionIntegrator`** | potential mass, interior **and** boundary faces — `buildForms()` | **none** |
| **`meq::SourceIntegrator`** | potential mass NON-LINEAR form, domain — `buildForms()` | **none, and it is the one MEQ owns** |
| **`mfem::HDGExtensionIntegrator`** | flux mass, boundary faces, curved and free-boundary only — `buildForms()` | **needs no kernel of its own**: `AssembleFluxMassBdrMatricesBatched()` batches the scatter and only checks the element block's size. Agrees to **8.65e-15** → **[M-185](MEASUREMENTS.md#m-185)** |
| `mfem::DomainLFIntegrator`, `mfem::VectorBoundaryFluxLFIntegrator` | the two right-hand sides, in the `assemble*Rhs()` helpers | none, and it does not matter: once per mesh |

**The hybridized path never reaches a PA entry point, and PA is not the route.**
Hybridized assembly goes through `BilinearForm::ComputeElementMatrix()`. Partial
assembly stays out of the plan anyway, because its economic case is sum
factorisation and MEQ is triangles.

**`AssemblyMode = "batched"` IS the route.** `DarcyForm::Assemble()` calls the
four `Assemble*Batched()` routines **in front of** each element loop. Each
returns false having done nothing unless `Batched` was asked for, and only then
does `ComputeElementMatrix()` run. *The reading that a route was missing came
from a doxygen comment on a member declaration, believed without reading the
call site in front of it. A comment documents an intent; the call site is the
fact.* On `machine-f-diiid` with all three batched keys, face assembly, local
factorisation, local solve and trace assembly all read yes, with `psi_ax` to
every digit (M-185).

**MEQ's face constraint takes the LINEAR route**, and so is in
`CanBatchPotFaceAssembly()`'s domain. The obvious inference runs the other way —
both HDG face integrators sit on the potential-mass **non-linear** form, by the
decision in `buildForms()` — and it is wrong, as `darcyform.cpp` shows:
* `GetPotentialMassForm()` is called only on the Picard/Anderson path, so on the
  Newton path `M_p` is null.
* `EnableHybridization()` falls past `if (M_p)` to
  `else if (Mnl_p && FaceIntegratorsAreLinear(...))`.
* `FaceIntegratorsAreLinear()` holds, because `meq::SourceIntegrator` makes the
  problem non-linear and both face integrators are plain
  `BilinearFormIntegrator`s.
* So the constraint is `c_bfi_p`, assembled once per solve.
* The stabilization gate in `HDGFaceScatterCanBatch()` passes too, because MEQ's
  hook is `meq::ConstantStabilization`.

*"Which form an integrator is added to" and "which route its constraint takes"
are different questions, and the second is decided by a predicate in somebody
else's translation unit.*

**The three batched keys, and what each is worth.** All three exist on
`GradShafranovSolver` and are settable from a file.
* **`TraceAssemblyMode` batched** is worth 7–9% of the DIII-D solve and is the
  default everywhere.
* **`AssemblyMode::Batched` is a 15–27% LOSS on the host** (M-99, M-101) — a
  measurement of the host and not a verdict on the route. It is the device route,
  carrying the 25 `forall` launch sites in `fem/darcy/` and `BatchedLinAlg`.
  *A cost is a property of a target, not of a route.*
* **`LocalFactorMode::Batched`** now keeps the condensation cache
  (`condensationCacheTaken()` reads yes) and is a small consistent gain on the
  host, about 4% on `ComputeH`. The host default stays `Serial`, because 1% of
  the wall is inside the serial arm's spread, and that is not reason enough to
  move a default never asserted bit-for-bit against the other route. **Under a
  Device it is the default, as a route decision**: `CanBatchLocalSolve()` opens
  on it, so with it off the element-local dense algebra stays on the host.
  Worth 1.54× on `ComputeH` under CUDA (M-185).

**`AssemblyMode` is `Serial` under a Device** → **[M-215](MEASUREMENTS.md#m-215)**,
**[M-209](MEASUREMENTS.md#m-209)**.
* **Why not `Batched`**: swept bitwise, `LocalFactorMode` and
  `TraceAssemblyMode` batched reproduce the host exactly, while
  `AssemblyMode::Batched` differs at 3.5e-16 to 4.7e-16. Its face kernels
  accumulate `D` through `AtomicAdd` from three or four faces per element, which
  would destroy the parity instrument.
* **What `Serial` gives up** is the face kernels alone, which lose on the host
  anyway. Two consecutive CUDA runs of `machine-f-diiid` now agree exactly.
* **Why `Serial` and not simply inheriting the default**: the host default is
  `Threaded`, which the driver refuses under a Device above one thread (M-188).
  *An inherited default has to be honourable in the configuration that inherits
  it.* → **[M-136](MEASUREMENTS.md#m-136)**.

**The ordered list, weighted by the leg profile** → **[M-80](MEASUREMENTS.md#m-80)**.

1. **`meq::SourceIntegrator`**, the `F(R, z, ψ)` term, evaluated per quadrature
   point per residual **and** per Jacobian. On MEQ's bordered configuration it is
   the dominant leg — **48.2% serial and 23.1% threaded** on the high-beta case,
   because the border spends four residuals per step. **Both sides of the device
   interface are built**:
   * MEQ's half is `meq::SourceView` (the POD), `sourceF()`/`sourceDFdPsi()`,
     and `sourceResidualAtPoint()`/`sourceGradientAtPoint()`
     (→ **[M-177](MEASUREMENTS.md#m-177)**).
   * `meq::makePlasmaSourceView()` covers coil-augmented sources by viewing the
     plasma term alone (→ **[M-191](MEASUREMENTS.md#m-191)**).
   * It is the single `CanBatch*` predicate MEQ fails, `CanBatchLinearResidual()`,
     and that refusal is correct, because the integrator genuinely is
     non-linear. Upstream's half is a sibling for that case, the batched
     caller-integrand residual (`CLAUDE_MFEM.md`).
   * **The POD is a second expression of the physics, and `SourceViewTests` pins
     it to the bit**; that test is load-bearing.
   * `INTERPOLATORY-HDG-PLAN.md` is the plan for replacing this leg with CCSZ
     interpolatory HDG. Read its falsifying experiment first: it gives 1.6–2.0×
     fewer `F` evaluations, `ComputeH()` untouched, and a 23–28% Amdahl ceiling.
2. **Face assembly**: reached through `AssemblyMode = "batched"`, and a loss on
   the host (above).
3. **`mfem::HDGExtensionIntegrator`**: batched through the scatter (above). It
   does **not** disqualify the condensation cache.
   * `Bnl_data` is written only from `grad_Aup`, the potential-dependence of a
     **block non-linear** flux law.
   * This integrator is a `BilinearFormIntegrator` on the flux mass bilinear
     form, so it lands in `A` and is assembled once.
   * It makes `A` asymmetric (5.4e-01) but not solution-dependent. *Asymmetry
     was read as state dependence.*
4. **Nothing for the two linear-form integrators**, assembled once per mesh.

## Single precision

**It is built last, after residency, which is the ordering decision worth
keeping.** A single-precision build destroys the parity instrument: a
disagreement could be a wrong kernel or a narrower format, and nothing
distinguishes the two. So residency is built in double, where it can be checked,
and converted afterwards. The early experiment was a single-precision **host**
build, which needs no device.

**Why it is worth having at all**:
* This card's FP64 runs at 1/32 to 1/64 of FP32, where a datacentre part runs at
  about 1/2.
* Under a partial offload the format hardly matters, because the loss is host
  serialisation (M-173).
* A fully resident solve, though, competes against eight host cores in the one
  format the card is crippled for.

**The `override` surface is smaller than a count of `double`s suggests**:
* All 33 `double`-returning overrides in `src/meq` are on MEQ's own MFEM-free
  `Source` / `Profile` / `Coil` hierarchies, which stay `double`.
* What couples to `real_t` is the twelve classes deriving from an MFEM type.

### What it measures

* **It works** → **[M-184](MEASUREMENTS.md#m-184)**. `machine-f-diiid` in single
  precision, on the host and under CUDA, gives `psi_ax = 3.759862e-01` against
  `3.759851e-01`. That is 2.9e-6 relative, FP32's floor, with the axis at the
  same point to four decimals.
* **The accuracy belongs to the case, not the format, over a range of 560×**
  → **[M-204](MEASUREMENTS.md#m-204)**.
  * Against the double build stopped at the **same tolerance** (the control),
    the four-border diverted machine is 2.2 × eps in `psi_ax` and the one-border
    limited case is 1,260 × eps.
  * Borders close at FP32's `eps`.
  * A single-precision equilibrium is qualified per problem, not per build:
    across the examples, `psi_ax` agrees with double from 6.9e-07 on
    `fixed-f-diiid` to 1.1e-03 on the `q(psi)` outer loop
    (→ **[M-221](MEASUREMENTS.md#m-221)**).
  * `NewtonRelativeTolerance = 1.0e-10` cannot be met. The run stalls at 3.3e-07
    and burns 400 iterations.
* **Both builds read 83 of 84, the same deliberate red**
  → **[M-303](MEASUREMENTS.md#m-303)**; parity was first reached at 62 of 63
  each → **[M-222](MEASUREMENTS.md#m-222)**. MKL 2026.2 and PARDISO's
  refinement cap move single's knife edges both ways
  (→ **[M-302](MEASUREMENTS.md#m-302)**).
  * `k+1` above `k = 1` is a double-build acceptance: round-off contaminates the
    flux by 203× at `k = 2` and by 1.5× on the coarsest mesh at `k = 3`. So it is
    skipped by degree, loudly, with the absolute ceiling still asserted.
    → **[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md)**, *A per-pair rate criterion,
    or a wider `rateSlack`, in single precision*.
  * The port is 38 files and 405 lines (→ **[M-217](MEASUREMENTS.md#m-217)**).
    The double build was **bit-for-bit unchanged at every step**, which is what
    the `std::max`-against-the-historical-value construction is for. A guard in
    the wrong place cannot weaken double.

### Tolerances and floors → [M-174](MEASUREMENTS.md#m-174), [M-217](MEASUREMENTS.md#m-217), [M-219](MEASUREMENTS.md#m-219), [M-220](MEASUREMENTS.md#m-220)

***A tolerance that does not scale with `real_t` is a build-configuration bug
wearing a numerical costume.*** Six constants were right for double, wrong for
single, and silent:
* `PointLocator::insideTolerance`.
* MFEM's `MFEM_VERIFY( err < 1e-10 )`.
* `normalisationDifferenceStep`.
* MFEM's `ExtensionBoundaryQuadrature` `fd_step`.
* `CriticalPointFinder`'s `tolerance`.
* `CriticalPointFinder`'s `jacobianStep`.

Each is now `std::max`'d or `#ifdef`'d against its historical value. **The worst
was MFEM's** (→ **[M-183](MEASUREMENTS.md#m-183)**):
`ExtensionBoundaryQuadrature` yields 21 points per face in double and **0** in
single, because its rejection floor follows `eps` while its `fd_step` does not.
It returns normally, so MEQ's exterior DtN coupling dropped out of the bordered
Jacobian while staying in the residual.

* ***A tolerance tuned against a solve that did not converge is tuned against
  nothing.*** `newtonRelativeTolerance = 1.0e-12` sits seven orders below
  single's residual floor of 2.2e-06, so 23 of the first 45 red tests never
  reached an assertion. This is a solver constant, not an assertion, and no
  assertion sweep could have found it.
* ***A default is only a default on the paths that do not override it.***
  `setNewtonControl()` was called at 39 sites with their own literals. There is
  now one shared `meq::newtonTolerance()`, following the `limiterPolygonMaximum()`
  and `extensionFdStep()` precedent.
* ***A tolerance a solver was GIVEN is the floor on what it can achieve.*** Seven
  sites asserted `< 1e-8` on a solve whose target floors at ten `eps`.
* ***A fixed-point iteration's floor is its map's, not Picard's.*** Newton
  squares its residual across the target; Picard contracts by a fixed factor and
  approaches its floor from above. Solving the condensed system for the whole
  field each step, Anderson-Picard stalls at 5.80e-05, because the element-local
  eliminations' forward error is in every iterate at full size. Evaluated as a
  correction on the NPC residual it stalls at 1 to 2.5 `eps`, so
  `meq::picardTolerance()` is `newtonTolerance()`'s `10 eps` (`solveByPicard()`).
  SUNDIALS has no compensated-summation route into KINSOL: that option is ARKODE
  SPRK's alone, and checked.
* ***A difference is not a derivative.*** Two exact routes agree to `eps`; a
  finite difference agrees only to `eps^(2/3)` to `sqrt( eps )`.
* ***A stationary point's POSITION costs half the digits and its VALUE costs
  none***: `sqrt( eps )` = 3.45e-04, measured on the limiter contact and the
  located axis.
* ***An estimator differentiates, so its floor RISES as the mesh falls.*** `eta`
  reads 8.93e-05 rising to 3.93e-04 over four refinements, against double's
  1.927e-08. Past that crossing, refinement genuinely raises the error, so
  **an adaptive loop in a narrow format stops being a descent**.
* **Assert the claim, not double's artefact of it.** "Newton takes one step"
  becomes "one step reaches the round-off floor", which is often the *stronger*
  test. *The robust half of a rate claim survives where the `k+1` half does
  not.* And *a controlled experiment whose signal is below the floor is a
  double-build acceptance, not a tolerance to widen*: loosening one would leave
  a case that cannot tell a 5% Jacobian error from an exact one.
* ***A suite count is not a measurement of a suite.*** One repair moved the count
  by 1 and the inventory by 18. Restoring the install moved the count *up*,
  because a test that had aborted at its first case now ran all of them. Read the
  distinct failing cases.

### What was not a tolerance

* **`TransformBack()`'s default tolerance** is `1e-7` in single, 0.84 of one
  float epsilon.
  * It is tighter than the `phys_rtol` of the `InverseElementTransformation`
    that `TransformBack()` constructs.
  * So on a straight triangle the Newton reaches neither exit and returns
    **`Unknown`**, which a caller testing `== Inside` reads as *not this
    element*.
  * `meq::inverseTransformTolerance()` fixes MEQ's five verdict-consuming sites,
    and is `1e-15` in double.
  * **It was a wrong answer and not only a slow one**:
    * The tracer fell back on `Mesh::FindPoints` 108, 240 and 95 times over three
      surfaces, against none in double.
    * A band excursion read 0.628 h against 2.027 h.
    * The warm start's lost nodes were cured by *constructing* the inverse, not
      by the better initial guess they were credited to. *A fix that works is
      not an explanation of why.*
  * **Fixed upstream**: `tol_0` is `1e-6` in single from `meq-integration`
    `81ecf3281e`, which is MEQ's own value, so the explicit argument and the
    default agree. MEQ keeps the argument, so it does not depend on the install
    -> [M-271](MEASUREMENTS.md#m-271).
* **`mfem::VertexConePath` aborts on a vertex lying exactly on `Γ`**. It
  differences a mean of ray endpoints from the vertex, and gets the zero vector.
  `meq::markingOffset()` keeps `MarkLevelSetSubdomain()` off the set that
  round-off cannot decide. It is exactly zero in double. It also closed
  `DriverAcceptance` marking 451 elements against 450.
* **Structural absences.**
  * UMFPACK is double only, so an example naming `TraceSolver = "umfpack"` fails
    in single. `meq::defaultTraceSolver()` is build-conditional, and **an example
    must name no trace solver**.
  * gslib has no floating-point abstraction, so `MFEM_USE_GSLIB` cannot be
    built beside `MFEM_PRECISION = single`. Every bordered free-boundary case
    fails cold even in double, so the warm start is load-bearing physics, and
    `meq::PointLocator` replacing gslib is what gives a single build one
    (M-174, M-176).
  * `locatePotentialPoint()` keeps the best near miss within
    `meq::samplerInsideTolerance()`.
  * **`install-single` is `install` minus six things, not minus precision**
    (→ **[M-218](MEASUREMENTS.md#m-218)**). SuiteSparse and gslib have no single
    build. MFEM has no guard that `SUNDIALS_PRECISION` matches `real_t`, and
    does not compile in that pairing without `sundials-real-t`. PARDISO resolved
    against `libmkl_sequential`. *A configuration stated in flags is not a
    configuration*, and a red test there cannot say which difference it met
    (M-210's attribution rule, applied to a failure).
* **A single build type-checks a seam that the double build collapses.**
  * A function returning both `double` and `real_t` compiles only when the two
    are the same type.
  * `std::max`'s refusal to deduce a common type flagged all 13 of the tree's
    `std::max( 1.0e-300, ‖x‖ )` guards. As a float comparison, a `1e-300` guard
    is below float's smallest subnormal and does not loosen; it disappears.
  * `meq::KnotView::psi` is `real_t` where `meq::Profile`'s abscissae are
    `double`. So an interval boundary moves by an ulp, and at a coil edge that
    switches the whole coil term on or off. That is a boundary predicate, not a
    tolerance.
* **`#ifdef MFEM_USE_SINGLE` in a translation unit that includes no MFEM header
  is silently false.** A test that drives the shipped binary and includes
  nothing of the library is exactly that translation unit.
  `mfem/config/config.hpp` defines the macro.
* **Coin tosses**, measured with M-186's instrument. `pedestalSelfConverges`
  and `PlasmaConnectivity`'s confined half each fail at `MKL_NUM_THREADS=1` and
  converge at higher values, and both are repaired by relaxing a
  **precondition**. `fixed-h-circular` crawls from a ramp at every draw and
  converges in 4 steps from a core at every draw, so `MachineFixedBoundary`
  asserts it in single with no skip (M-296).
  `diverted-tokamak` lands on a branch the driver's axis-source guard refuses —
  *the safety net works; the branch does not* — while `diverted-tokamak-xpoint`,
  on the same machine and mesh, solves.
* ***A run that converges is not a run that converged to your fixed point.***
  `XPointBorder`'s weight 4 finished in 195 iterations against 11 for the others,
  in a neighbouring basin, reporting `converged`. The discriminator is the
  *work*.

### The setup loss was a bisection with no stopping rule → [M-242](MEASUREMENTS.md#m-242)

Single precision read as a 1.1× to 2.4× loss, while every threaded solve leg was
10–27% faster per call.
* The loss was two single-threaded host legs, `setup` and `re-assembly`, doing
  identical work.
* The suspect, `extensionFdStep()`, was wrong. `perf` named
  `meq::BoundaryShape::shiftedAngle`, which is a `double` in both builds, so its
  2.99× self time was a **call count**.
* `mfem::VertexConePath`'s `tol_` is `1e-13`, and `FirstCrossing()`'s bisection
  floors at `eps*length`. That is six orders above `1e-13` in single, so every
  ray ran to `max_iter = 100`.
* It misreports nothing, so only a clock could see it.

**The repair is the termination, not a tolerance.** `meq::pathTolerance()` was
built, swept and reverted whole: every value loose enough to help moved an
answer, taking `ConductorSubtraction`, `DriverAcceptance` and
`MachineFixedBoundary` red. Upstream's `238f719f50` stops the bisection when the
midpoint stops moving, `if ( tm <= ta || tm >= tb )`.
* Float then stops at about 20 halvings.
* Double still stops at 39, so double is bit-identical.
* All three path classes reach that one `FirstCrossing()`.

Upstream pins it **in double**, with a tolerance no double bracket can reach and
a count of calls: *a format-specific defect is pinned by putting the other format
in the same position*. With the fix, single is 0.81× to 0.90× of double per
Newton step, and 0.69× to 0.81× on setup.

## This machine's GPU

An RTX 2070 SUPER with CUDA 13.3, `sm_75`: enough to write and debug a device
path and check its answers. It is **not** representative. Consumer FP64 runs at
1/32 to 1/64 of FP32 where datacentre parts run at about 1/2, so a timing taken
here can invert the conclusion a production part would give. Treat local device
runs as correctness evidence and take timings elsewhere. `TODO` carries the
detail.
