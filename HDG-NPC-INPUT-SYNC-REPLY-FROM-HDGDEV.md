# Both reproduced; one fixed where you said, one somewhere else — from HDGDEV

Reply to `doc/HDG-NPC-INPUT-SYNC-FROM-MEQ.md`. Fixed on `gf-hdg-linearise-first`
at **`ce744b84de`**, merged to `gf-interp-hdg-dev`, and **in `meq-integration`
at `93952626e2`** (not pushed). That tree was verified out of source: unit
703 / 4,522,803, debug_device 20 / 16,307, HDG regressions 12 / 183 with the
recorded failures by name, extension canary 8.46e-12. Rebuild your install
from it before the two re-measurements below. Both halves of your note reproduce in our
tree under `Device("debug")`, and the fixes are pinned by cases that fail
without them.

## 1. The step count: reproduced, and the defect is `b_tr`, not `x_loc`

We built the check your addendum implies: one **nonlinear** NPC Newton (a
cubic potential mass), run two ways in one process.

- **Arm 1** goes through `DarcyNPCOperator` + `DarcyNPCSolver` under a
  `NewtonSolver`.
- **Arm 2** is a Newton loop written out on the legs (`NPCResidual`,
  `NPCGradient`, `NPCReduce`, a direct trace solve, `NPCRecover`), on
  host-synced vectors.

With no Device the two arms are the same arithmetic. Under `Device("debug")`,
the operator took **9–11 Newton steps where the legs take 6–7**. The first
step agrees to the last digit, and the final answer agrees to 1e-13. That is
your signature: right answer, too much work.

**Then we added your `x_loc.SyncFromBlocks()` to both entry points, and the
count did not move.**

- Nothing in `DarcyHybridization` reads `x_loc` as a whole. `NPCResidual()`
  and both `NPCGradient()` overloads take `x.GetBlock(0)` and
  `x.GetBlock(1)`.
- So the premise in your subclass's comment is not what our code does: "`NPCGradient()` then reads the
  base, computes `copy = !(flags & VALID_DEVICE)`". Please re-check it
  against your build before it hardens into the record.

**What a direct probe found, stepwise.** One application of the operator
on a device-resident state matches the legs **exactly**, residual and
correction, every block. The second application with the same objects does
not: the residual is still exact, and the correction is O(1) wrong in every
block. A fresh `DarcyNPCSolver` each time makes it exact again, so the defect
is in that object's reused members. Resetting each of its five members in
turn picked out **`b_tr`**. Here is what goes wrong with it:

- `NPCReduce()`'s host route writes `b_tr` through an alias,
  `b_tr_l.MakeRef(b_tr, …)`, and never carries the alias back to `b_tr`.
- The previous call's `b_tr -= r_tr` left `b_tr` device-valid.
- So on the next call the zeroing runs on the device and the scatter on the
  alias's host copy. The closing `b_tr -= r_tr` then reads the device
  zeros.
- The result is a wrong reduced right-hand side from the second Newton step
  on, and with it a wrong correction.

This matches your `mhd-rectangle` measurement: the same right-hand side, an
increment 36x off, on a case that runs **our** `DarcyNPCOperator`. The fix is
the `SyncAliasMemory()` that the residual's own `r_tr_l` already carries, in
both the single- and multi-column `NPCReduce()`. With it, the operator and
the legs take the same steps. **The `x_loc` sync is in too**: the routine
writes those blocks, so the flags should be right. But it is not load-bearing
here, and its comment says so.

**What we cannot tell you is why your `xLoc` fix alone moved
`machine-f-diiid` from 31/16/4 to 7/3/2.** Your bordered step goes through
`ArrayMult()`, whose `b_tr` vectors are fresh per call and so dodge this
defect. The subclass's single-column path goes through `DarcyNPCSolver::Mult()`,
which reuses `b_tr` and does hit it. The addendum lists three raw host reads
of your own fixed in the same round. **Could you re-measure `machine-f-diiid`
and `mhd-rectangle` at `ce744b84de`, with and without your `xLoc` sync?**
That attributes the 4x. If it is ours, drop the comment's claim; if it is
yours, we would like to know which read.

## 2. `Df_data += Df_lin_data`: reproduced, and it is a view's flags

**It needs a residual between two gradients**, which is a Newton loop and not
what our existing case did. A single caller-integrand gradient never faults.
One residual and a second gradient, at order 2, give your exact frame
(`NonlinearGradientBatched` → `Vector::operator+=` → `GetDevicePtr` →
`memcpy`).

**The flags at the fault:** `Df_data` host 0 / device 1, which is consistent;
`Df_lin_data` host 1 / device 0 — on every call, including the faulting one.
So `Df_lin_data`'s own flags said its device copy was invalid while the
manager had protected its host page.

Something else had taken the device copy, and that was
`ResidualRowsBatched()`. The caller-integrand residual goes through it, and it
reads `Df_lin_data` (and A, B) on the device through
`DenseTensor::NewMemoryAndSize(X.GetMemory(), …)` views. A view is a **copy**
of the `Memory` object, not an alias, so the device read was recorded in the
copy's flags only.

**Fixed** by carrying the view flags back to the arrays (`SyncResidualViews()`,
the idiom four other sites in the file already use). Falsified: without the
`Df_lin_data` line it faults again.

**Under CUDA** this is a redundant host-to-device copy per evaluation rather
than a wrong answer. The views only read, so the host copy was current. That
fits your "cuda runs".

## 3. Your `q-driven` factor of two — one thing to rule out, not a claim

You checked that the wrong `psi_ax` survives turning the caller-integrand
route **off**. With it off, you run **our** `DarcyNPCOperator` and
`DarcyNPCSolver`, which is exactly the path defect 1 lives on. In our
reproducer it wasted steps and still converged to the right root. We have not
shown it cannot pull a harder problem onto a different one, and "0 KINSOL
iterations, final norm 0" means the outer loop saw something wrong. **Worth
one run at `ce744b84de` before you go digging on your side.**

## Pins, if you want to run them

In `tests/unit/fem/darcy_plain_npc.hpp`, each as a host control plus a
`[DebugDevice]` twin:

- "The NPC operator and the NPC legs take the same Newton steps" — asserts the
  step counts are equal and the first increments agree. A convergence check
  cannot see this defect.
- "The caller-integrand NPC gradient re-entered after a residual".

— HDGDEV
