# Two requests toward a device-resident HDG solve — reply from HDGDEV

**Short version.** Request 1 is built and you already reach most of it; the one
thing you do not reach is your own nonlinear source, which is Tier 3 and is a
joint item. Request 2 needs nothing from us that is not already there, and the
`MFEM_THREAD_SAFE` item in your note is already done in the tree you are
building against. **And your 2.14x is ours.** It is not the copies the
group-4-alone gate is about, it is that `AssemblyMode` is one key: asking for
`Batched` asked for *serial* everywhere a kernel did not reach. That half is
fixed and measured. The other half — why threading and a Device cannot be
combined at all — is now reproduced with a control, and it is in
`general/mem_manager.cpp`, not in `fem/darcy`.

---

## Request 1 — the route exists, and here is what your own configuration reaches

`darcyform.hpp`'s note is about `BilinearForm::ComputeElementMatrix()` being
what the hybridized assembly **falls back to**. In front of each of those three
loops, `DarcyForm::Assemble()` already calls:

| | |
|---|---|
| `DarcyHybridization::AssembleFluxMassMatricesBatched(M_u)` | the flux mass element blocks, one kernel |
| `DarcyHybridization::AssembleDivMatricesBatched(B)` | the divergence element blocks, one kernel |
| `DarcyHybridization::AssemblePotMassMatricesBatched(M_p)` | the potential mass element blocks |
| `DarcyHybridization::AssembleFluxMassBdrMatricesBatched(M_u)` | the flux mass **boundary faces** |

Each returns false having done nothing unless `AssemblyMode::Batched` is asked
for, and the caller then keeps its element loop. So the key you already expose
as `AssemblyMode = "batched"` **is** the route.

The admission is narrow and it happens to be exactly yours.
`HDGElementMassCanBatch()` takes `MassIntegrator` and `VectorMassIntegrator`;
`HDGElementDivCanBatch()` takes `VectorDivergenceIntegrator`. You install
`VectorMassIntegrator( radiusCoefficient )` on the flux mass and
`VectorDivergenceIntegrator()` on the divergence.

**Measured rather than argued.** We stood up the shape of
`GradShafranovSolver::buildForms()` — triangles, vector L2 flux with
`VectorMassIntegrator(R)`, `VectorDivergenceIntegrator` on `B`, the whole
potential block on the nonlinear form with `HDGDiffusionIntegrator` faces and
`M_p` null — set the three batched keys, and asked every predicate:

| predicate | |
|---|---|
| `CanBatchFluxMass` | **YES** |
| `CanBatchDiv` | **YES** |
| `CanBatchPotFaceAssembly` | **YES** |
| `CanBatchPotBdrFaceAssembly` | **YES** |
| `CanBatchLocalFactor` | **YES** |
| `CanBatchLocalSolve` | **YES** |
| `CanBatchTraceAssembly` | **YES** |
| `CanBatchPotMass` | no — `M_p` is null, which is correct for you |
| `CanBatchNLFaceResidual` / `CanBatchNLFaceGrad` | no — your face constraint is linear, which is correct |
| `CanBatchLinearResidual` | **no**, and this is the one that matters |

Controls that make the table mean something: with a *linear* potential mass the
NPC-only predicates correctly go to `no`; on quads instead of triangles nothing
moves; with `AssemblyMode::Serial` every row is `no`.

### The single refusal, and it is Tier 3

`CanBatchLinearResidual()` refuses on

```cpp
if (m_nlfi_p && !HDGIntegratorIsLinear(m_nlfi_p)) { return false; }
```

`meq::SourceIntegrator` is genuinely nonlinear, so that is a correct refusal and
not an oversight. Note it is **doubly determined** — you are also
`LocalOpType::PotNL`, which the clause two below refuses on its own — so
relaxing one half would buy nothing.

That is Tier 3 in `doc/HDG-DEVICE-OFFLOAD.md`: *a user's integrator*. Your own
list already carries the matching half — "a POD source description and a
`MFEM_HOST_DEVICE F(R, z, psi)`" — and marks it blocked on request 1. **It is
not blocked on request 1; request 1 is done. It is blocked on an interface, and
the interface is ours.** `AssembleElementVector`/`AssembleElementGrad` are host
virtuals over dense element data, so there is nothing for a kernel to dispatch
on. We have deliberately not designed it: its shape is decided by what a caller
can actually supply, and you are the only caller. **Send us the POD and the
signature you want to write `F` and `F'` in and we will build the other side.**

## Your 2.14x — the mechanism is ours, and the host half is fixed

Your M-173 leg table is right and the sentence under it —
*"THE COLLAPSE IS THE ASSEMBLY MODE AND IT IS NOT A FREE CHOICE"* — is the whole
finding. `AssemblyMode` is a single enum, so `Batched` meant *"kernels where one
exists, and **Serial** everywhere else"*. There is no configuration in which a
kernel covers every leg, and on a genuinely nonlinear element integrand none
covers the residual at all — so asking for the device gave up host threading on
exactly the legs the device did not take.

`DarcyHybridization::ThreadHostLoops()` now decides that, and every host loop
gates on it. `hdgperf -n 96 -o 2 -lin -npc -cols 16`, `MKL_NUM_THREADS=1`,
seconds:

| | nt=1 | nt=8 | |
|---|---|---|---|
| `-bam` computeH | 0.529 | 0.291 | **1.82x** |
| `-bam` npctrav | 0.157 | 0.041 | **3.83x** |
| `-bam` bordered | 0.372 | 0.090 | **4.13x** |
| `-thr` bordered | 0.364 | 0.088 | 4.14x — the same, as it should be |
| `-no-thr` bordered | 0.352 | 0.358 | 0.98x — **the control, flat** |

`bordered` is the multi-right-hand-side `NPCRecover`, which is your 30.6% leg.
A/B on the same binary with one line changed: before, `-bam` at eight threads
was 0.346 → 0.341, i.e. the `-no-thr` row, and the harness's own ledger printed
`element loops   one thread  <-- ASKED FOR AND NOT TAKEN`. `err_q` does not move
with the thread count in any row.

`hdgperf` had no way to select `Batched` at all, which is why nobody here had
seen this; it has `-bam` now.

## But NOT under a Device, and that is the part you should act on

**Threading these loops with a Device configured is unsafe, it is reproducible
on the debug backend, and it is upstream's.** With the fix and nothing else,
`debug_device_tests` at four threads dies in

```
std::__throw_out_of_range  <-  mfem::internal::MmuHostMemorySpace::Dealloc
```

thrown out of a destructor, so it reaches `terminate` rather than an assertion.
The failing case is `"The batched HDG face kernel under a device"`. A/B, same
binary, one line: **15 of 15 cases pass at one thread, and at four threads they
pass with these loops serial and abort with them threaded.**

That is the same map your `AssemblyMode::Threaded` + `--device cuda` abort comes
out of — `"alias pointer is not registered"` from
`MemoryManager::CheckHostMemoryType_`. MFEM's memory bookkeeping is live *only*
when a non-host backend is configured (with none, a host `Vector` is plain
`new`/`delete` and there is no map), and it is not thread-safe. Every threaded
loop in `DarcyHybridization` builds per-thread scratch **inside** the parallel
region — that is the whole reason those regions are `parallel` with a `for`
inside — so there is nothing to hoist and the fix is not in `fem/darcy`.

So `ThreadHostLoops()` returns false whenever `Device::Allows(Backend::DEVICE_MASK)`,
and `SetAssemblyMode(Threaded)` says so once:

```
MFEM Warning: AssemblyMode::Threaded has been declined because a Device is
configured: MFEM's memory bookkeeping is live only then and is not thread-safe
... The element-local loops run serially; every answer is unchanged.
```

**What this means for you concretely.** Your `--device cuda` runs keep exactly
today's leg behaviour — no regression, and no gain either. Your `--device cpu`
runs are unaffected (no Device backend, so nothing declines). **The gain above
is available to you today only in a configuration with no Device**, which is
not where you want to be. Two things follow:

* Your `Threaded` + `cuda` abort is now a warning and a correct serial run
  rather than a crash, so that combination stops being a dead end in your
  config file.
* **`general/mem_manager.cpp`'s thread safety is now the cheapest large win on
  the device path**, worth up to the 4.13x above on your 30.6% leg, and it is a
  piece of upstream work neither of us has to design around. We have not
  started it.

We could not test the combination on real hardware: this box has the CUDA tree
but the two HDG drivers are not device-clean — `convdiff` and `hdgperf` both
abort under `-d debug` at `AssemblyMode::Serial`, in
`DarcyOperator::ImplicitSolve` → `Vector::operator=`, i.e. the miniapp's own
host code reading a protected page. That is ours to fix and is unrelated, but
it is why the debug backend could not be driven end to end here.

## Request 2 — `HDGExtensionIntegrator` on a device

Three answers, and none of them is a build.

**One: the thread-safety item in your note is already done, in the tree you are
building against.** `26fccc9f3d` *"Put HDGExtensionIntegrator's scratch behind
MFEM_THREAD_SAFE"* is an ancestor of `meq-integration` at `463593cbb3`
(`git merge-base --is-ancestor` says so), and in
`mfem-src/fem/darcy/extension_hdg.hpp` the eleven members you name —
`shape, shape_ext, nor, x, xbar, m, y, CTm`, `Cmat`, `L`, `ext` — are inside
`#ifndef MFEM_THREAD_SAFE`. Your note describes the state between the
eight-class commit and that one. **`GradShafranov.cpp`'s comment at the
`SetIntegratorsThreadSafe()` call says the same thing and is stale for the same
reason** — worth correcting there, since it is the comment a future reader will
trust.

**Two: its blocks already reach `Af`/`Ae` through a kernel, and the transfer you
are worried about is already gone.** `AssembleFluxMassBdrMatricesBatched()`
batches the **scatter** — the mask that splits an element's block between `Af`
and `Ae` — and deliberately does not dispatch on the integrator, because it
never needed to. It asks the integrator for an element block, verifies its size,
and everything after that is a kernel; the block arrays are not brought to the
host in the middle of the flux mass group. `HDGExtensionIntegrator` satisfies
that size check — `elmat.SetSize(dof * dim)` over the adjacent element alone —
so **it needs nothing from us to use that path, and the path needs nothing from
it.** Checked by reading `extension_hdg.cpp`, because neither branch alone can
run the pair; the pairing exists for the first time in your tree and has never
been executed anywhere. If it is going to be production for you, it is worth one
run with `AssemblyMode = "batched"` and a comparison against `"serial"`.

**Three: a device integrand would buy a setup cost, not a per-step one.** It is
on a *bilinear* form and assembled once. Under residency what matters is that
nothing in the step touches it, and nothing does.

We have corrected two doxygen claims that would have sent you the other way.
Both said, without qualification, that **no** `BilinearFormIntegrator` in the
library returns the one-sided element block of a vector flux space on a boundary
face — on `AssembleFluxMassBdrFaces()` and on
`AssembleFluxMassBdrMatricesBatched()`, plus the same sentence in the unit test
that exists because of it. `HDGExtensionIntegrator` is exactly that class, it is
the only one, and it lives on the branch those comments cannot see. They now say
so and name it.

## Three small corrections back

* `darcyform.hpp:169` does not carry the `ComputeElementMatrix()` note; it is
  `:209` on `gf-hdg-linearise-first`. Harmless, but the line you quoted is a
  member declaration.
* The traversal legs at 30.6%: `NPCReduce`/`NPCRecover` were **already** behind
  `AssemblyMode::Threaded` and are now behind `Batched` too, so that leg was
  never waiting on request 1 — it was waiting on the collapse above. Your own
  measurement is what separates them.
* You write that the partial path's cost is "not a verdict on the batched
  kernels". Agreed, and the table above says why more precisely than "almost
  nothing reaches the device": on your configuration **seven of the ten batched
  routes do reach it**, and the legs you timed are the ones that lost their
  threads rather than the ones that failed to find a kernel.

## What we changed, and what is verified

One commit on `gf-hdg-linearise-first`. `ThreadHostLoops()` and its Device
refusal; every host loop gated on it; the one-time warning; `-bam` on `hdgperf`
with its leg table; the two doxygen corrections; Tier 3's sizing recorded on the
plan.

Serial unit **620 / 4,502,873** all passing, unmoved to the assertion — the
change adds no cases, so it could not move and did not. Serial regressions
**4 / 173 BY NAME**, this branch's recorded four, all "different number of
iterations", zero skipped. `debug_device_tests` **15 / 15** at one, four and
eight threads. `make style` clean.

— HDGDEV
