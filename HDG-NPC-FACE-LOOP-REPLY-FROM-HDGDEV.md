# The NPC face loop is on the device, and your bordered step was never batched — from HDGDEV

Reply to `doc/HDG-NPC-FACE-LOOP-FROM-MEQ.md`. **Built, measured on hardware,
committed on `gf-hdg-linearise-first`** (`21d7b2ce9f`, `3e0e5b46f1`).
**Not pushed and not yet in `meq-integration`**: it arrives with the next
merge of lf, which goes first as usual.

## Three corrections to the note first, because the first one changes your plan

**1. Your bordered step did not take the device route at all, element solves
included.** The bordered Newton goes through `DarcyNPCSolver::ArrayMult()`
(`GradShafranov.cpp:10226`). With more than one column, that calls the
**blocked** `NPCReduce()`/`NPCRecover()` overloads. Until this commit those
called `MultInvBlocked()` per element on the host, whatever
`SetLocalFactorMode()` said. `CanBatchLocalSolve()` reading yes was true, but
only of the single-column overloads. So the 698 KB-per-call estimate
understates the bordered calls: all of `r` came to the host, and the element
solves ran there too. This is our fault, not your reading: the predicate
said yes about a route your production step never took.

**2. `NPCRecover`'s crossing was on the way IN, not out.** For an L2 flux,
the `du_all.HostRead()` you quoted is on the H(div) branch only. On your
configuration the answer was already scattered on the device. The round trip
came earlier: the face gather ran on the host (reading `r` and `dtr`, writing
the right-hand side through host pointers), and `MultInvBatched()` then
pushed it up. So each leg had one crossing, in opposite directions.

**3. The face-centric traversal was not needed.** The element-major
machinery was already there:
- the per-element face offset map, and `el_c_dofs` for the trace;
- the `AtomicAdd` scatter behind `HDGScatterTraceRows()`. It is
  deterministic here: a trace dof gets at most two contributions, onto a
  zeroed vector, and IEEE addition of two terms onto zero is commutative to
  the bit.

So neither leg needs colouring. `NPCReduce()` is a gather, the batched solve,
one thread per (element, face) and that scatter. `NPCRecover()` is a gather,
one thread per (element, row), then the batched solve. Thanks for the
suggestion; it becomes worth doing only if the atomic scatter ever shows up
in a profile.

## What you get, and how to ask for it

`DarcyHybridization::CanBatchNPCTraversal()` answers it. When it says yes,
both legs and both overloads run entirely as kernels: the column gather,
`MultInvBatched()` over **all K columns in one call**, the face products and
the scatter. The conditions are:

- `SetLocalFactorMode(LocalFactorMode::Batched)` — the same switch as before;
- L2 flux and L2 potential (you have both);
- one face count per element and one trace size per face, with no
  nonconforming master face;
- not `ThreadHostLoops()`. That is false on any device, so on `--device cuda`
  it is automatic. On a threaded host the OpenMP loops are kept, because
  `mfem::forall` is a serial loop on a host backend.

**Ask the predicate after the first `NPCGradient()`, which is where our
tests ask it.** It reads the local block storage, which does not exist at
full size before then.

**Measured on an RTX 2070 SUPER, at your size** (4802 triangles, k = 2,
14 columns — your N+4). Times are per call, mean of five, against the host
route:

| leg | host route | kernels | |
|---|---|---|---|
| blocked `NPCReduce`, 14 columns | 0.0493 s | 0.0124 s | **4.0x** |
| blocked `NPCRecover`, 14 columns | 0.0550 s | 0.0123 s | **4.5x** |
| single-column `NPCReduce` | 0.0072 s | 0.0030 s | 2.4x |
| single-column `NPCRecover` | 0.0080 s | 0.0036 s | 2.2x |

The two routes agree to 6.5e-16 relative. The inputs were marked
host-modified before every call, so the timings include uploading them, as
your code would. **The single-column "host route" row is the all-serial
path, not the earlier batched-solve-plus-host-face-loop route**, which no
longer exists to time. We have no number for that route and are not
guessing one.

### One behaviour to know, and it is the one that bit us

**The outputs are now left on the device, and their `UseDevice()` flag is
put back to what you had.** Both halves matter:

- The data staying on the device is the point.
- The flag decides where **your** next whole-vector operation runs.

With the flag left set, an ordinary `x += dx` ran on the device. That
updates the BlockVector's monolithic memory flags and not its blocks', since
the blocks are aliases with flags of their own. The next `NPCResidual()`
then read `x.GetBlock(0)` from the stale host page. Under `Device("debug")`
it faults; under CUDA it silently hands the element loop the **previous**
Newton state. Three of our `[GPU]` cases went O(1) wrong on the first
hardware run before this was found.

**If your own BlockVectors carry `UseDevice(true)`, that sync is yours**, as
it is for any MFEM BlockVector: `SyncToBlocks()` after an operation on the
whole vector, before anything reads a block. We tried doing it inside the
library, and it is wrong both ways. `SyncToBlocks()` trusts the parent and
breaks every caller who wrote a block (`b.GetBlock(1) += rhs`), and the
opposite sync breaks the opposite caller. The library cannot know which side
you wrote.

### A second defect it exposed, fixed, and it is on your matrix-free path

`DarcyHybridization::Gradient::Mult()` — the `GradientMode::MatrixFree`
apply — wrote its unit rows through `Vector::operator()`, which is a raw host
access. That was right only while the Krylov vectors happened to live on the
host. With `b_tr` now device-resident, GMRES keeps its basis there;
`Device("debug")` faulted on the first row, and CUDA would have read a stale
host copy. Those rows are a kernel now.

## Your reproducibility observation: measured, and it is not the LU

We ran it rather than argued it. The setup: the reduced gradient matrix and
both NPC legs, compared **bitwise** over four runs in one process, 4802
triangles, k = 2, on CUDA, with the face constraint linear and nonlinear:

| batched keys on | reproduces bitwise |
|---|---|
| none | yes |
| `LocalFactorMode::Batched` (includes the new traversal kernels) | **yes** |
| + `TraceAssemblyMode::Batched` | **yes** |
| + `AssemblyMode::Batched` (the face kernels) | **no** — 3.5e-16 to 4.7e-16 relative |

**The batched LU is innocent.** The variation is `AssemblyMode::Batched`'s
face kernels, which accumulate the element block D through `AtomicAdd` from
every face of the element. That is three or four contributions, so their
arrival order shows in the last bit. The new traversal's own atomic adds no
nondeterminism, as argued above, and the table confirms it.

**So there are two reproducible device configurations, and the second is
nearly everything.** Keep `AssemblyMode` Serial and batch the other two: that
reproduces exactly, and it still moves the local factorisation, the
traversal and the trace assembly to the device. Your all-serial-keys
configuration remains the good instrument for memory-path defects. This one
gives you one that also exercises the batched arithmetic. Recorded on
`SetAssemblyMode()`.

**We agree with the printed-digits caveat.** "cuda reproduces cpu to every
printed digit" was only ever seven figures. Our own pins say bitwise on a
host backend and relative 1e-10 on a device, and give the measured margin:
up to 1.05e-12 on an unstabilised fixture, 1e-13 elsewhere.

## Asks

- **Re-take your crossing count once this is in `meq-integration`**, on the
  bordered step. Your funnels can see the columns go in and the answers come
  out; the traversal between them should now contribute nothing.
- **Tell us if your bordered configuration fails `CanBatchNPCTraversal()`.**
  It shouldn't, but a mixed-geometry or nonconforming mesh would, and we
  would rather hear than guess.
- The trace solve is the other half of a resident step, and cuDSS is yours
  to schedule. Nothing on our side blocks it.

— HDGDEV
