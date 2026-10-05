# Reply: the two uploads left after `969569440e`

Re `HDG-DEVICE-BUS-MEASUREMENTS-FROM-MEQ.md` (2026-10-04). Thank you for the
numbers. 543 → 102 MB up with every answer unchanged is the result we were
hoping for, and your by-request table is the attribution we could not make
from here.

Of your two new items, **item 2 is ours and is fixed**. **Item 1 is not ours**,
and below is how to find which of your vectors it is.

* `gf-hdg-linearise-first` at **`3c183cde03` and `35ab01a477`**;
* `gf-interp-hdg-dev` at **`0ce31ebcc2`** (merges lf);
* **`meq-integration` at `4929772e95`.** Verified there from clean: unit 715 / 4,524,232 and debug device 22 / 16,687, all passing. All seven HDG miniapps build, and `make style` is clean. On a CUDA build of lf the GPU cases pass 109 / 109. Rebuild your install from it.

## 2. `HDGTraceEliminateRows`'s 12 × 1,310,436 B: fixed

That buffer is the reduced gradient's column-index array, J. Every gradient
gets a fresh `SparseMatrix`, and `BuildTraceHPattern()` gave it a **host** copy
of the cached pattern. Every reader of a gradient's pattern is now a kernel:
`HDGTraceDiagIdentity()` reads I, and `HDGTraceEliminateRows()` reads J, which
is where your gdb caught it. So both went up again at every linearisation.
The comment on the cache said the pattern "never leaves the host because
every consumer wants it there". That was never true for cuDSS. It is
withdrawn in place.

For a gradient, the pattern is now copied **on the device** from the map's
device copy, which goes up once per mesh. The linear matrix (non-NPC
`Finalize()`) keeps its host pattern, since it is assembled once.

**Measured** with our probe on your route with
`TraceAssemblyMode::Batched` (4802 triangles, k = 2, 14 columns, 12 steps,
RTX 2070):

* the gradient phase's uploads fall from 41 copies / 19.1 MB to 18 / 3.8 MB;
* J goes from 12 uploads to 1, and I from 13 to 1;
* `r`, `r_tr` and every column of `Y` are bitwise unchanged.

Our probe solves the trace with UMFPack. Because the pattern now lives on the
device, it brings I down once per solve (12 × 87,616 B); it was already
downloading J. **cuDSS should download neither, but we have no cuDSS build
to confirm that.** Please check it in your nsys numbers: J should be gone from
the uploads, and no new 1.3 MB or 88 KB download should appear.

## 1. `ResidualRowsBatched()`'s 16 + 17 uploads: your load, we think

`ResidualRowsBatched()` gathers five vectors with `GetSubVector()`:
`darcy_u`, `darcy_p`, `bu`, `bp` and the trace. `bu` and `bp` are the
**blocks of the load** you pass to `NPCResidual()`. They have the same two
sizes as the state blocks, so the byte counts cannot tell them apart.

We separated them in our probe:

| iterate (blocks) | load (blocks) | residual uploads of those sizes |
|---|---|---|
| host-valid | host-valid | 13 + 14 |
| device-valid | host-valid, never rewritten | **none** after the first |
| host-valid | device-valid | 12 + 13 |
| device-valid | parent made device-valid only | **none** after the first |

(Answers bitwise identical in every arm.)

So in our hands:

* **The state copy follows your iterate's blocks.** Your `xLoc` blocks are
  `UseDevice(true)` and device-valid, so `darcy_u = x.GetBlock(0)` runs device
  to device and `GetSubVector()` finds `darcy_u` already there. If your
  iterate were the cause, the upload would be a frame earlier, inside
  `Vector::operator=` in `NPCResidualStart()` or in your own `xLoc` copy.
  It would not be in `GetSubVector()`.
* **A load that is not rewritten goes up once and stays valid.** That holds
  even when only its parent was marked device-valid: the blocks' first device
  read uploads them, and afterwards both copies are valid.

An upload on **every** residual from `GetSubVector()` therefore means the
load's blocks are host-valid again at each `NPCResidual()`. Either something
rewrites the load on the host between residuals (`FormLinearSystem()`, or a
right-hand side that depends on the state), or the `BlockVector` over it is
re-made, so its block aliases start again from host-only flags. Your control,
making `darcyRhs` device-valid once per build, would not survive either.

**One-line check:** at entry to `NPCResidual()`, print
`load->GetBlock(0).GetMemory().HostIsValid()` and `DeviceIsValid()`. Or, in
the gdb frame, compare the `GetSubVector` call's `this` against
`&load->GetBlock(0)`. If the load really changes every step, one upload per
change is the price of computing it on the host. If it does not change, find
what marks it host-valid.

## Something we found on the way, and fixed

On the **linear, non-NPC** route, `TraceAssemblyMode::Batched` left the
reduced matrix's values on the device. `EliminateTraceTrueDofs()` then ran
`SparseMatrix::EliminateRowCol()`, which reads and writes I, J and A raw on
the host:

* under `Device("debug")` it faulted;
* under CUDA the essential rows were silently never eliminated.

NPC never reaches that loop, so your route was not affected. It now syncs
first (`35ab01a477`), and it is pinned by a debug-device case.

## Your other points

* **The curved-path movement** (`miller-adaptive`'s first residual in the
  7th digit, and the two `1e-14`s) is consistent with `683646da4b`. Your old
  install was `5c2324e921`, so that commit and the bus change are the only
  two new to you. The bus change was bitwise identical in our probe. The
  vertex-on-`Gamma` fix changes a path only where a vertex of `Gamma_h`, or a
  quadrature point, sits within round-off of `Gamma`, and an adaptive mesh is
  where that happens.
* **`ComputeH()`'s read-only downloads** (`Df_data`, `Df_ipiv`, about 24 MB)
  are still open, as we said. Nothing on your route reads them on the host.
  Removing them means making their host readers sync for themselves, and
  that is not done.
* **`NPCResidualFinish()`'s trace-row downloads** come from your `rTr`, as you
  say. Doing for `rTr` what you did for `rLoc` should remove them.

## What we would like back

* Your nsys totals from `4929772e95` against the 102.3 MB / 51.4 MB column.
* Which vector the residual's `GetSubVector()` uploads turn out to be.
