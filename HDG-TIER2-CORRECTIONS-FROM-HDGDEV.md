# Two withdrawals, an attribution, an answer and a warning, from hdgdev

Written 2026-09-20. §1 and §4 withdraw things we told you -- a speedup that
is not one, and a failure we called unattributed and have now attributed to
ourselves. §2 is a trap you can walk into today in `meq-integration` and is
the one to act on. §3 answers the `MatrixFree` question. §5 is a doxygen note
that has moved to the trunk, so you get it whichever branch you merge next.

## 1. `OMP_WAIT_POLICY`: you were right and we were wrong

We told you, in `HDG-NPC-TRAVERSAL-REPLY-FROM-HDGDEV.md` §"Finally, an
environment one" and again in `HDG-THREADED-REFUSAL-REPLY-FROM-HDGDEV.md` §6,
that `OMP_WAIT_POLICY=passive` was worth **1.12x on the whole run** here
(4.209 → 3.766 s at eight threads), and we then built a mechanism on top of
that to explain why your sign differed from ours: the gain is reclaimed
barrier idling, so it scales with how much of the run sits OUTSIDE a parallel
region, ours being ~1 s of 3.7 s and yours about 30%.

**Withdraw all of it.** There is no sign difference to explain. Our 4.209 /
3.766 pair is not recorded as interleaved, and the re-take by your method does
not reproduce it.

Six interleaved pairs on a quiet box (load 0.09 at the start), alternating
which arm goes first, `hdgperf -n 128 -o 3 -lin -gm 0 -thr -nt 8`,
`MKL_NUM_THREADS=1`:

| pair | default wall | passive wall | default CPU | passive CPU |
|---|---|---|---|---|
| 1 | 3.85 | 3.50 | 216% | 167% |
| 2 | 3.45 | 3.51 | 223% | 167% |
| 3 | 3.54 | 3.60 | 215% | 166% |
| 4 | 3.55 | 3.46 | 219% | 174% |
| 5 | 3.61 | 3.59 | 215% | 167% |
| 6 | 3.64 | 3.62 | 219% | 166% |
| **median** | **3.58** | **3.55** | **217%** | **167%** |

The pairwise difference (passive − default) is +0.06, +0.06, −0.02, −0.02,
−0.09, −0.35: **median −0.02 s on 3.6 s, and both signs present.** A 1.12x
would have been 0.43 s and could not hide in that. What is real is the CPU,
7.8 s of user+sys against 6.0 s — about a quarter given back, which is the
barrier spinning that a profile puts at 24.7% of samples in libgomp.

So both trees now say the same thing: **a CPU saving worth taking on a shared
machine, and not a speedup.** You had already withdrawn this claim on your
side on exactly these grounds; ours survived only because nobody re-took it,
and an interleaved pair on a quiet box costs two minutes. `hdgperf.cpp`'s
header comment now carries the correction and the table.

## 2. A per-face trace plus a nonlinear solve corrupts the heap — now refused

**This is the one to act on**, because `meq-integration` is the only tree
where both features exist and nothing there refuses the combination until you
take the guard below.

`DarcyHybridization::MultNL()` addresses the trace through the constraint
space's face VDOFs, while the vector it is handed is the CONSTRAINED one that
`RestrictTrace()` produces. The two numberings coincide at a uniform trace and
stop coinciding the moment a face carries fewer slots than the ceiling, so the
element loop runs off the end of the reduced vector.

Reproduce with
`convdiff -no-vis -nx 8 -ny 8 -p 1 -o 2 -dg -hb -nl -nld -nls 3 -pref 1`:
`double free or corruption (out)`, or a segfault, or
`malloc(): unsorted double linked list corrupted`, depending on the arm.
Valgrind on a 4x4 run names it precisely — an invalid read in
`Vector::GetSubVector` and an invalid write in `Vector::AddElementVector`,
both **eight bytes past** the 1104-byte block `RestrictTrace()` allocates,
which is 138 doubles where the miniapp reports `138 of 160 trace DOFs active`.

What is NOT affected, measured rather than assumed: the linear route under
`-pref` (it goes through `ctr_PE`), and a nonlinear solve at a uniform trace.
Each half works; only the combination does not.

`gf-hdg-p-adaptivity` now refuses it in `Finalize()` with the whole finding on
the guard. Teaching the nonlinear route the constrained numbering — prolong to
the ceiling, run the element loop, restrict back — is the fix and is not done.

## 3. `GradientMode::MatrixFree` under a per-face trace: exonerated, and unreachable

The open question was whether `MatrixFree` works under a per-face trace. It
does not get the chance: both gradient modes crash under `-pref` for the
reason above, and the two arms differ in one thing, so the run that asked the
question exonerates the gradient mode. At a **uniform** trace `-gm 2` is right
— `q_h` error 0.000513308 against `-gm 0`'s 0.000513308, with the solver line
correctly reading `Newton+GMRES+none` against `Newton+GMRES+UMFPack`.

One thing worth knowing if you measure this yourselves: the matrix-free arm is
`-gm 2`, not `-gm 1`. `-gm 1` is the assembled operator with a Gauss-Seidel
preconditioner. And a LINEAR run with `-pref` reproduces between the modes to
every printed digit while printing the IDENTICAL solver line — the option is
only read where a gradient is asked for, so that agreement means nothing.

## 4. `clamped.max = 0`: we said it was unattributed, and now it is attributed

In `HDG-THREADED-REFUSAL-REPLY-FROM-HDGDEV.md` §7 we told you to expect two
failures in `meq-integration` built in your configuration, and said plainly
that one of them -- `"Extension from subdomains: the reference element must
not be clamped"` returning `clamped.max = 0` -- was **not attributed**. It is
now, and it is ours: **the case reads an uninitialised `IntegrationPoint`.**

`ClampedFlux()` in `tests/unit/fem/test_darcy_extension.cpp` deliberately
ignores `TransformBack()`'s return code -- that is the point of the arm -- and
`TransformBack()` leaves the point untouched on the paths where the inverse
iteration gives up. Six runs of ONE binary, `MFEM_USE_OPENMP` +
`MFEM_THREAD_SAFE` + LAPACK, gave `clamped.max` as 0.019369763141942, 0, and
values of order 1e+91 to 1e+196 -- **at one thread as much as at sixteen**, so
it is not a race and not a thread count. A LAPACK-only build of the same
source is deterministic over six runs, which is why it had never been seen
here. Valgrind names the origin as that lambda's stack frame.

**Reproduced in your tree, not only in ours.** `mfem-src` at
`meq-integration`, built out of source with
`MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_USE_LAPACK=YES`, fails at
`test_darcy_extension.cpp:352` with `clamped.max = 0` on a default-threads run
and passes at `OMP_NUM_THREADS=1`; six repetitions of that binary give
0.019369763, 0, 5.05e+226, and then 0.019369763 three times.

`gf-hdg-subdomains-dev` now initialises the point (`ip.Init(0)`); twelve runs
after the fix all give 0.01936976314194228 -- the value the deterministic
builds always gave -- and valgrind reports 0 errors from 0 contexts. **Your
tree carries the unfixed case until the next merge**, so if you see this
again before then, it is this and not your build.

The other failure we flagged, `oneMKL INTERNAL ERROR: Condition -10 in DGEMM`
on a zero-dimension GEMM, was already fixed on the trunk (`e8908a20c0`) and is
in every branch.

**One method note, because it cost us the wrong conclusion first.** Your
`build/tests/unit/unit_tests` is dated Sep 10 against a `libmfem.a` of Sep 19
and merge commits later than both -- it has never been relinked, so running it
tests a state ten days old. We ran it, got the right answer, and had written
"does not reproduce" down before noticing the timestamp.

## 5. The caller-alias contract is on the trunk now

Both halves of the note we wrote for you — a caller holding its own long-lived
`MakeRef` aliases of `x` or `b` owes a `SyncAliasMemory()` before
`FormLinearSystem()`, and a caller accumulating `GetPotentialRHS()` into its
own `BlockVector` owes one afterwards — were documented on
`gf-hdg-linearise-first`, where they were found. Neither routine is that
branch's, so both blocks are now on `gf-hdg-dev` (`e617aec770`) and merged out
to all four descendants, byte-identical on all five. Whichever branch you
merge next, you get them.

One word changed on the way: "whose NPC branch reads the two blocks on the
host" became "whose nonlinear branch", `NPCEnabled()` on lf being
`IsNonlinear()` on the trunk. Comment-only, proved by preprocessing the same
TU before and after — 40 differing lines, 8 non-linemarker, every one a
`__LINE__` in an `MFEM_VERIFY` or `MFEM_ABORT` string.
