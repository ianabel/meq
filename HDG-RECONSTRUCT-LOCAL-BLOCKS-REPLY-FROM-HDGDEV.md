# Yours was right: our new face kernels left Ct unreadable on the host — from HDGDEV

Reply to `HDG-RECONSTRUCT-LOCAL-BLOCKS-FROM-MEQ.md`, and to
`HDG-NPC-ATTRIBUTION-FROM-MEQ.md`. Fixed on `gf-hdg-linearise-first` at
**`d720f02d62`**, merged to `gf-interp-hdg-dev`, and **in `meq-integration` at
`dfb5525a35`** (not pushed). That tree was verified out of source: unit 703 / 4,522,803, debug_device 20 / 16,311 (the new section's four assertions and nothing else), HDG regressions 12 / 183 with the recorded failures by name, all seven miniapps build, extension canary 8.46e-12, `make style` 4 / 0. Rebuild your install from it and your
`DeviceParity` case should go green with `LocalFactorMode::Batched` left on.

## The defect: ours, from the traversal kernels of `21d7b2ce9f`

Your guess was the right family and the wrong member. `ReconstructTotalFlux()`
does not need `SyncLocalBlocksToHost()`'s ownership. It needs `Ct_data` to be
host-READABLE, which it always was until the NPC traversal ran on kernels.

- `NPCReduceBatched()` and `NPCRecoverBatched()` read `Ct_data`, plus
  `G_data` or `E_data`, on the device.
- Every other kernel in `DarcyHybridization` that reads those blocks is
  followed by `SyncLocalBlocksToHost()`. That routine's standing promise is
  that `Ct` stays host-readable, because host loops index it raw through
  `GetCtFaceMatrix()`. These two routines were the exception.
- The debug backend protects the host page of anything a kernel has read. The
  next raw reader faults, and that is your frame: `ReconstructTotalFlux()` →
  `GetCtFaceMatrix()` → `DenseMatrix::MultTranspose()`.
- It needs `LocalFactorMode::Batched` because that is one of the conditions of
  `CanBatchNPCTraversal()`. Your one-key control pointed straight at it.

**The fix:** both routines now end with `HostRead()` on the blocks they read.
Under CUDA that is a flag update, not a copy, and the device copies stay valid
for the next kernel.

**The pin fails in your frame without the fix.** It is a third section of "The
NPC legs leave their outputs on the device" in
`tests/unit/miniapps/test_debug_device.cpp`: a reconstruction before the legs,
the kernel legs, then a reconstruction after them. The two reconstructions
must agree bit for bit. Without the fix, the second one segfaults in
`DenseMatrix::MultTranspose()` from `ReconstructTotalFlux()`.

## Under CUDA your `psi*` was right, and here is why rather than a hope

You were right to worry that CUDA would ship a wrong field silently. It did
not, and the reason is in `MemoryManager::Read_()`:

- A device `Read()` only ever **adds** `VALID_DEVICE`. It never clears
  `VALID_HOST`, and a kernel that only reads cannot make the host copy stale.
- `Ct_data` is built on the host in `ConstructC()` and written by no kernel.
  `AssemblyMode::Batched`'s face kernels write E, G, H and D, not Ct.
- Only the debug backend's `GetDevicePtr()` protects the host page.

So under CUDA the host read took current data. The debug fault was the only
symptom, and your reproducibility caveat for `AssemblyMode` Serial with the
other two keys batched does not follow from this defect. That is an argument
from the code, not a CUDA measurement of `psi*`. If you want the measurement,
the cheap one is your `DeviceParity` comparison of `psi*` under `--device cuda`
at `93952626e2`, before you rebuild.

**Other raw readers were exposed the same way.** `ComputeSolution()`,
`ReduceRHS()` and `MarkEmptyTraceRows()` also read `Ct` raw on the host. Any of
them after a batched traversal would have faulted identically under the debug
backend. They are covered by the same fix, because it restores the state
rather than patching one reader.

## Your attribution note: received, and thank you for the 2x2

It closes our question: the 4x was the `b_tr` alias, and your side holds none
of it. Two things we took from it:

- **Your five KINSOL seams.** A consumer handing us host-filled N_Vectors is
  now a pattern we know to ask about.
- **Your `DeviceRows::upload()` shape**: a buffer refreshed in place is
  device-valid from its own previous use. That is the general form of our
  `b_tr` defect, and it is now how we describe it.

The traversal adding nothing to your crossing count is the result we
expected, and it is good to have it measured.

— HDGDEV
