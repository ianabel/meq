# The condensation's per-face solves, blocked over an element's faces

From `hdgdev` (`/home/ian/projects/mfem-hdg-dev`, `gf-hdg-linearise-first`),
2026-09-15. **Not committed and not pushed** — this is a request for one
measurement we cannot take here, plus the reading of yours that made it worth
making.

## What changed

Two sites in `fem/darcy`, both the same defect: the same factorisation applied
to one right-hand side at a time when the right-hand sides were already
contiguous.

**(a) `DarcyHybridization::ComputeElementH()`.** The condensation formed
`A⁻¹Cᵀ` and `S⁻¹(B A⁻¹Cᵀ − E)` **one face at a time** — `nf` calls of
`LU_A.Solve(na, nc)`, `nf` of `mfem::Mult(B, ·)` and `nf` of
`LU_S.Solve(nd, nc)`, all against the same two factors. Each face's block
already sits contiguously with leading dimension `na` / `nd`, so the element's
faces together *are* one multi-column right-hand side, and
`LUFactors::Solve(m, n, X)`'s `n` has always been the column count. Now one
call each.

At your level-3 dimensions (`na = 12`, `nd = 6`, `nc = 3`, `nf = 3`, M-80) that
is **three `dgetrs` of 3 columns becoming one of 9, three times over**, per
element per gradient. Flops are unchanged; what falls is call count.

**The `S` solve is the one that runs on every gradient.** The condensation
cache holds the bracket `B A⁻¹Cᵀ − E`, which does not move across a Newton
loop, but `S` moves with `D` — so that solve is paid whether the cache is live
or not. The other two are paid when the cache is cold or refused.

**(b) `HDGPotentialPostprocessor::Compute()`** factored the local Neumann
matrix once per element and then solved it once per field. It is now one
`DenseMatrixInverse::Mult(DenseMatrix&)`. **This is a no-op for MEQ** — it
bites only at `vdim > 1` — and is listed for completeness.

## Why this is being sent to you rather than just done

**M-89 item 7 writes this work off, and on your own later evidence the reason
it gives no longer holds.** It names the element-local dense work as
`mkl_lapack__dgetrs_` 13.4%, `mfem::Mult(DenseMatrix, …)` 10.5% and `MultNL`
12.7%, and says reaching it means `AssemblyMode::Batched` — which M-99 then
measured as a **15 to 27 per cent LOSS** on a host, in all three rounds.

(a) reaches the first two of those with **no mode at all**. It is
unconditional, in the serial route you already run, and it is orthogonal to
`AssemblyMode`, `LocalFactorMode`, `TraceAssemblyMode` and the condensation
cache.

Your own numbers size the leg it acts on: `ComputeH()` is **68–73%** of the
gradient leg on every row of M-80, `ComputeElementH` is **72.5% measured / 82%
by flop count** of that, and the gradient leg is 33.3% / 33.4% of a threaded
unpinned step. So it acts inside roughly **17%** of a whole solve — the figure
M-80 already calls a floor.

## What we verified, and what we cannot

Verified here on `gf-hdg-linearise-first`, clean rebuild of both trees:

| | |
|---|---|
| serial unit | **593 cases / 4,500,565 assertions** — the baseline 592 / 4,500,487 plus one new case |
| `debug_device_tests` | **13 cases / 16,194 assertions**, `Device("debug")` |
| serial regressions | **4 / 160**, the same four BY NAME, all `DIFFERS` on iteration count |
| parallel unit, 2 ranks | **102 / 70,966** |
| parallel regressions | **17 / 129**, the recorded 8 DIFFERS + 9 FAILING |
| the guard fires, (a) | leaving one face's block unsolved fails **34 of 68** Darcy cases, 174 assertions |
| the guard fires, (b) | rotating the field columns fails **2 of 3** `[Postprocess]` cases |
| the guard fires, step 1 | disabling the replay fails the fill-count assertion in 4 places |
| the guard fires, step 2 | it failed for real, twice, before the memory-type fix above |

**We cannot measure the point of the change.** `MFEM_USE_LAPACK` is `NO` in
both HDG trees, so `LUFactors::Solve` is `LSolve()` then `USolve()` over
columns: the arithmetic is bitwise what the per-face loop gave, and the timing
cannot move. On your MKL build it is a genuine multi-RHS call.

**One caveat, and it is the same footnote `TraceAssemblyMode::Batched`
already carries.** On LAPACK a triangular solve and a GEMM may block
differently for nine columns than for three, so the reduction can reassociate
and the last bits may move. Anything of yours holding two configurations to a
*bitwise* comparison should expect that across this change and compare to a
tolerance instead.

## What we would like measured

On `examples/machine-f-diiid.toml`, `build-nocuda`, `OMP = MKL = 8`, the
configuration you ship (threaded / serial / batched):

1. `mkl_lapack__dgetrs_` and `mfem::Mult(DenseMatrix, …)` as shares, before and
   after. These are the two entries M-89 item 7 names.
2. `GetComputeHTime()` before and after — the stablest number in either of your
   tables, and the one this change should move if it moves anything.
3. Whether `psi_ax` still reads `3.759851e-01` to every digit. We expect the
   answer to move in the last bits only, per the caveat above.

The patch is `fem/darcy/darcyhybridization.{cpp,hpp}` and
`fem/darcy/postprocess_hdg.cpp` on `gf-hdg-linearise-first`, plus two test
files. The multi-RHS half is trunk material — the per-face pattern is on every
branch including `gf-hdg-dev` — and will move there before it reaches you
through `meq-integration`; the cache work builds on `CondensationCache`, which
is this branch's.

## M-101 is now stale in your favour: the trade is gone

**`LocalFactorMode = "batched"` no longer costs the condensation cache.** M-101
reads it as "a **trade** and not a knob that does nothing ... it buys a batched
local factorisation and pays the condensation cache", which was an accurate
reading of the code as it stood. It has been changed, in two steps.

The refusal was about STORAGE and never about arithmetic. Under
`LocalOpType::PotNL`:

* `A` is **not refactored on a gradient at all** — `FactorElementsBatched()`
  carries the same `lop_type != PotNL` guard the serial route does, so the
  factorisation was never what was being recomputed;
* the Schur complement moves with `D` and is cacheable in neither route;
* what that route recomputed needlessly, every gradient, is `A^-1(-/+B^T)` and
  the four constant face-pair products — which is what the cache holds.

**Step 1** hands `FactorElementsBatched()` the cache's own buffer as its
destination, so there is one owner and the second gradient replays instead of
solving. **Step 2** does the same for the batched face-pair kernel: the bracket
`B A^-1 C^T - E`, `C A^-1 B^T + G` and `-C A^-1 C^T` are kept, and a replayed
gradient is then one Schur solve, one product and the H diagonal. By MAC count
at your dimensions that second step is the larger half — 3888 of ~4770.

**What this changes for you:** `condensationCacheTaken()` now reports `yes` on
the `LocalFactorMode = "batched"` row of M-101's table, and M-99's flat
`LocalFactorMode` axis should stop being two effects of opposite sign. Whether
it becomes a win rather than a wash is a measurement only your build can take,
and M-99 is the natural place to re-run it.

**Two counters for exactly that**, because a cache cannot be told from a
recomputation by its answer:
`DarcyHybridization::GetBatchedAiBtSolves()` and `GetBatchedFaceFills()`, with
`ResetBatchedCacheCounts()`; statics over function-local accumulators like
`GetComputeHTime()`, always on. Over N gradients of one Newton loop the first
reads **1** and the second reads the chunk count **once**, not N times. If
either grows with N in your configuration the replay is not firing and we want
to know.

**And one device finding that is worth having whatever you decide about the
modes.** The first version of step 2 segfaulted under `Device("debug")` and,
at order 0, silently returned a trace operator with 24 nonzeros where the first
gradient had 132. The cause is not HDG-specific and is worth knowing if you
ever cache a buffer a kernel writes: **`Vector::SetSize(int)` preserves the
memory type it already has**, so `UseDevice(true)` followed by `SetSize` still
allocates an unregistered HOST block. An alias `Memory` into an unregistered
base is unregistered too, and `Memory::CopyFrom()` then takes its "neither is
Registered" branch — a straight `memcpy` between host pointers, one of which
the debug backend has protected. Naming the type, `SetSize(n,
Device::GetDeviceMemoryType())`, is the fix; with no Device configured that is
`HOST` and a host build is unchanged.

Also fixed on the way, and it was a pure loss rather than a trade:
`CanCacheCondensation()` refused on the *mode* while `FactorElementsBatched()`
refuses on the mode **plus five storage conditions**, so a caller with
non-uniform blocks — essential flux dofs, a mixed-element mesh, variable-order
elements — got the fallback loop AND no cache, in exchange for nothing. Both
now ask one shared predicate.

## What we read of yours and are NOT asking you to repeat

M-99 and M-101 retired four questions before they were asked, and the reading
below is what we took:

* **`TraceAssemblyMode = "batched"` is worth 7–9% and is now on.** We were
  about to recommend exactly that. Withdrawn as already done.
* **`LocalFactorMode = "batched"` turns off the condensation cache** (M-101).
  We had this from the other side — `CanCacheCondensation()`'s
  `lfac_mode == Batched` refusal — and your framing was the better one: a
  trade, not a neutral knob, and it explains M-99's flat axis as two effects of
  opposite sign. **This is the entry the section above supersedes.** The two
  had never been measured against each other in either tree, because the
  predicate made running both impossible; they can be now.
* **`AssemblyMode = "batched"` is a 15–27% loss on a host.** That also settles
  a question we had about the enum: `Threaded` and `Batched` are mutually
  exclusive values of one field, which is an API artefact rather than a
  constraint — but on your evidence choosing between them currently costs
  nothing, so it is not worth changing.
* **You are on `LocalOpType::PotNL` with the cache taken.** Confirmed by your
  own accessors; nothing left on the table there.

## One question that is still open, and it is one line

M-89 item 7's percentages — are they **self** time or inclusive of children?
Both entries are leaves so self time is the natural reading, but it decides
whether the dense work is ~24% of your run or a subset of the legs already
counted, and therefore what (a) can possibly be worth. `perf report`'s
`--no-children` is the flag.
