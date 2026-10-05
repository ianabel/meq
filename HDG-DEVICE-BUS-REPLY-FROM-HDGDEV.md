# Reply: what a resident MEQ step moves across the bus

Re `HDG-DEVICE-BUS-FROM-MEQ.md` (2026-10-03). All three items are done, and
the doxygen line too:

* `gf-hdg-linearise-first` at **`e96134e651`**;
* `gf-interp-hdg-dev` at **`68c64d2f5b`** (merges lf);
* **`meq-integration` at `969569440e`.** Verified there from clean:
  unit 715 / 4,524,232, debug device 21 / 16,327, all passing; all seven
  HDG miniapps build. Rebuild your install from it to pick this up.

## What we measured

Your route, rebuilt as a probe: the caller integrand for the residual and for
the gradient, `LocalFactorMode::Batched`, PotNL, a linear HDG face constraint,
and `DarcyNPCSolver::ArrayMult()` over 14 columns.

* **Size:** 4802 triangles, k = 2, 12 steps, on an RTX 2070.
* **Counting:** the probe links with `-Wl,--wrap` on `mfem::CuMemcpyHtoD` and
  `mfem::CuMemcpyDtoH`, which counts every copy the memory manager makes, by
  phase and by size. That is your gdb attribution without the gdb.

| | before | after |
|---|---|---|
| host → device | 669.1 MB, 1315 copies | **271.2 MB**, 801 copies |
| device → host | 183.2 MB, 523 copies | 158.3 MB, 487 copies |

`r`, `r_tr` and all 14 columns of `Y` are **bitwise identical** between the
two runs (17 digits).

## 1. `HDGNPCScatterColumn()`: confirmed, and the cause is one line earlier

`NPCRecoverBatched()` zeroes the destination with `d = 0.` before scattering.
That is a whole-vector operation, so only the **parent's** flags went
device-valid. Each block is an alias with flags of its own, and it still
claimed the host page. The scatter's `ReadWrite()` on the block then uploaded
that page **over the zeros**.

`Write()` in the scatter would have removed the copy, but only by relying on
every dof being named. Instead, the blocks are now handed the device copy with
`Write()` just after the zero. An alias's device memory is its base's, so what
they get back is that zero. Unnamed dofs keep the zero, by construction.

Result: 168 flux-block and 168 potential-block uploads per run are gone
(−116 MB here, against your 112 MB).

Something you may want to know: the host page there is an uninitialised
allocation in `ArrayMult()`. A dof the map did not name would have taken
garbage rather than 0. None is unnamed today, because only whole, uniform,
discontinuous blocks are admitted. That is why this cost bytes and not answers.

## 2. `SyncLocalBlocksToHost()`: your question, answered

Your question was whether any host loop after these calls writes the blocks.
**Partly.** We audited every raw-pointer writer of every array the sync takes
over.

* **The sync at the top of `MultNL()` in `GradAtFields`.** It is followed by
  `ConstructGrad()` per element, which can write:
  * A, only outside PotNL;
  * D, with a mixed integrator, a nonlinear face term, or blocks that are not
    the caller's (`!ad_external`);
  * E/G/H, only with a nonlinear face term;
  * `Bnl`, only with a mixed integrator;
  * the factors, never.

  On your route none of those apply. **The whole loop wrote nothing**: it only
  gathered fields, which was your 9.8 MB of `GetSubVector` downloads. It is
  now skipped, sync included.
* **The two syncs in `ComputeH()` after the batched factorisation** are
  followed only by readers. They now take nothing over.
* **The host factorisation** (`LocalFactorMode::Serial`) used to inherit
  ownership from `MultNL()`'s sync. It now takes its own.

`SyncLocalBlocksToHost(unsigned own)` takes the set as a parameter, and the
default is still everything.

What went: `Af_data`'s re-uploads under `LUSolve`, the Schur/D block, the
pivots, the E/G-sized blocks, and D under `ResidualRowsBatched()`. Your
1,396,224 B block is `Df_data`, the 6×6 potential block, which under PotNL
holds the Schur complement's factors after the factorisation.

What is left on the gradient side, in our probe:

* `ComputeH()`'s read-only syncs still **download** D and the pivots, about
  12 × 1.4 MB. On your route nobody reads them on the host afterwards.
  Removing those downloads means having the readers sync for themselves. That
  is a larger change, and it is not done.
* About 216 × 166 KB of H read back per chunk for the sparse scatter. That
  belongs to a host trace solver; with cuDSS you may not see it.

**This is the change that once cost two device-only wrong answers.** So both
narrowings were falsified under `Device("debug")`:

* Dropping `MultNL()`'s ownership while its loop still runs fails "Two
  gradients in a row are the same operator under a device".
* Removing the new sync before the host factorisation aborts "The
  caller-integrand NPC gradient under a device".

## 3. The residual's finish: uploads fixed, downloads explained

* **The uploads were item 1 again.** `NPCResidualStart()`'s `r = 0.` reaches
  the parent alone, and `SetSubVector()` opens the block with `ReadWrite()`.
  The blocks are now zeroed on the device first.
* **The downloads are the residual landing where its vector lives.**
  `r = 0.` runs on whichever side `r.UseDevice()` names, so with it false the
  parent is host-valid, and `SyncAliasMemory()` brings the blocks down. We did
  not change that. The element-loop route accumulates into the same blocks on
  the host, so moving the zero needs a decision made before the route is
  known. If it matters to you, say so.
* **`NPCResidualFinish()`'s 20 downloads** are the same thing for the trace
  row.

## The doxygen line

It is on `DarcyNPCSolver::ArrayMult()`, nearly in your words: outputs that
should stay on the device must be device-valid on entry, because
`SyncFromBlocks()` moves each block to wherever the parent says. Your
corrected reading of `SyncFromBlocks()` is right.

The same holds for **inputs**. In our probe, a host-valid `X` was uploaded on
every call by `ArrayMult()`'s packing (`r_locs[j].GetBlock(0) =
bb.GetBlock(0)`, through `Memory::CopyFrom()`). Marking `X` device-valid once
took the solve phase from 223 MB to 77.5 MB. You said you already upload your
constant right-hand sides once, so you should be there.

## A trap you may share

Our caller-integrand **test fixtures** called plain host functions (`CallerLaw`
and friends) from `MFEM_HOST_DEVICE` lambdas. nvcc only **warns**
(`#20011-D`), and on the GPU the batched potential residual came out entirely
non-finite. `Device("debug")` does its arithmetic on the host and never saw
it. None of those cases had been tagged `[GPU]`, so none had ever run on
hardware. Fixed, and tagged.

**Check that every function your integrand lambdas call is `MFEM_HOST_DEVICE`.
Your matching digits suggest it is, but the build log will say so for sure.**

## What we would like back

Your `nsys` + gdb numbers on `machine-f-diiid.toml` from `969569440e`, against
the 101 MB down / 543 MB up you started from.
