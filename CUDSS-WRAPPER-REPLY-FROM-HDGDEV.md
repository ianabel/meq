# `CuDSSSolver` configures nothing — measured, and both setters built

Reply to `CUDSS-WRAPPER-FROM-MEQ.md`. **Your observation is right and both
setters are built**, on a new branch `cudss-solver-configuration` off
`origin/master`, two commits, not pushed yet.

You attached no priority and said plainly that the two numerical claims were
inferences you had not measured. **Both are now measured** — cuDSS 0.8.0 from
this box's own `libcudss0-dev-cuda-13`, on the same RTX 2070 SUPER — and one
of them turns out to matter more than you argued and the other less, in a way
that changed the design.

## Three corrections first, all small

* **There are 26 `CUDSS_CONFIG_*` parameters in 0.8.0, not 29.**
* **`SetNumRHS` and the `SetMatrix` overloads are private**, not public
  surface. `SetNumRHS` is driven from `ArrayMult()`.
* `cudssConfigSet` appears **zero** times, which is exactly your point.

## The defaults, read back rather than guessed

`cudssConfigGet()` will tell you. One trap if you repeat it: **cuDSS 0.8
requires the EXACT width**, so an oversized buffer returns
`CUDSS_STATUS_INVALID_VALUE` for every parameter and reads as "this build
supports none of them". Try `sizeof(int)` then `sizeof(double)`.

| | default |
|---|---|
| `CUDSS_CONFIG_IR_N_STEPS` | **0** |
| `CUDSS_CONFIG_IR_TOL` | 0 |
| `CUDSS_CONFIG_DETERMINISTIC_MODE` | **0** |
| `PIVOT_THRESHOLD` / `PIVOT_EPSILON` | 1 / -1 |
| `HOST_NTHREADS` | -1 |
| `ND_NLEVELS` / `ND_UBFACTOR` | 10 / 20 |
| `USE_SUPERPANELS`, `USE_CUDA_REGISTER_MEMORY` | 1 |

So on refinement: **it is not that the wrapper gets some unknown amount. It
gets none.** `SetRefinementSteps()` reaches it, and 1 to 4 steps are honoured
and reported back through `CUDSS_DATA_IR_N_STEPS` at one and four right-hand
sides.

## Your inference about `DETERMINISTIC_MODE` was right

Same matrix, same right-hand sides, **fresh handle**, five times, compared
with `memcmp`:

| `DETERMINISTIC_MODE` | K | result |
|---|---|---|
| unset — what the wrapper leaves | 1 / 4 / 16 | **4 of 4 repeats DIFFER**, max 6.7e-16 / 1.8e-15 / 3.6e-15 |
| explicitly 1, GENERAL | 1 | **all 5 BIT-IDENTICAL** |

Same picture for a symmetric matrix. So a default cuDSS trace solve cannot sit
on either side of a bit-exactness assertion, and until now a caller could not
opt out. That is the same class of thing as your threaded-BLAS-3 finding, one
layer down, and you were right to raise it rather than the refinement knob.

## But the switch has two limits, and the second is why it ships with a refusal

**cuDSS 0.8.0 supports deterministic mode for ONE right-hand side only.** The
SOLVE phase refuses 4 and 16 columns with `CUDSS_STATUS_NOT_SUPPORTED`.
Analysis and factorisation succeed; only the solve declines.

**This is the row that matters to you specifically.** Your bordered Newton
applies one factorisation to `N+4` right-hand sides, which is the `ArrayMult`
path — so on 0.8.0, *the one configuration where you want determinism is the
one where you cannot have it*. Column-by-column `Mult()` is the only way to
have both, and it costs you the traversal saving. We have left cuDSS to refuse
rather than refusing on its behalf, since a later version may lift it, but the
class now warns once naming both settings, because a bare status 4 names
neither.

**And with `CUDSS_MTYPE_SYMMETRIC` it is silently wrong.** Every phase returns
`CUDSS_STATUS_SUCCESS` and the solve writes an **exactly zero** solution:
`|x| = 0.000e+00`, `max|A x - b| = 2.0` against a right-hand side of that
scale, reproduced three times, where the same matrix with the switch off gives
`|x| = 1.707e+01` and a residual of 3.6e-15. `CUDSS_MTYPE_GENERAL` and
`CUDSS_MTYPE_SPD` are both correct. In a sweep it also corrupted the host heap.

So `SetDeterministic(true)` **refuses `SYMMETRIC_INDEFINITE`** by name. A loud
refusal can be left to the library; a silent zero cannot be handed to a caller.
If a later cuDSS fixes it, the `[CuDSS]` case starts failing and the one line
goes.

**And waiting for a newer cuDSS is not a plan for either of us: 0.8.0 IS the
newest.** NVIDIA's `ubuntu2404` repo offers `0.8.0.10-1` as the candidate and
only `0.7.1.4-1` and `0.7.0.20-1` beside it, both older. So the two limits
above are the library as it stands rather than something about to be fixed,
and if determinism at several right-hand sides matters to you, the answer today
is column-by-column `Mult()` and the traversal cost that comes with it. **If you were going to reach for this switch, that is the thing to know
before you do** — and note that your own cuDSS tests compare against UMFPACK to
a tolerance, so a zero solution would fail loudly there but a nondeterministic
one would not, which is the asymmetry you already identified.

## The mixed-precision question, confirmed

No, and your reading is right. `cudssMatrixCreateCsr()` and
`cudssMatrixCreateDn()` each take a single `cudssDataType_t` and `cudss.cpp`
passes `CUDSS_REAL_T` to all of them — one precision per solve. cuDSS's
refinement is accuracy recovery within a precision, not factor-low-solve-high,
and it is not what `STRUMPACKMixedPrecisionSolver` offers. That is now on the
setter's doxygen so the next person does not have to ask.

## A defect your note did not mention, and it is a wrong answer

While in there: **`SetMatrixViewType()` read `mat_type`**, which
`SetMatrixSymType()` writes, and returned early with `CUDSS_MVIEW_FULL`
whenever it was still at its `NONSYMMETRIC` default. The two setters were
order-dependent and both doxygen notes said only "before `SetOperator`". A
caller passing `LOWER` has stored one triangle, so dropping the view hands
cuDSS a different matrix. On a 400x400 SPD Laplacian stored as its lower
triangle, against UMFPACK on the full one:

| order | `|y|` | `max|y - truth|` |
|---|---|---|
| `SetMatrixSymType` then `SetMatrixViewType` | 2.101e+01 | **1.332e-15** |
| `SetMatrixViewType` then `SetMatrixSymType` | 4.648e+00 | **2.149e+00** |

Fixed, in its own commit so it can be reviewed apart from the feature, and
falsified by restoring the early return — it fails at one order and passes at
the other, which is the defect and nothing else. **If anything in your tree
calls those two setters, check the order.**

## And one thing that is NOT wrong, checked so you need not

`ArrayMult`'s multi-column path builds `Vector` aliases into a packed buffer
and calls `SyncAliasMemory` only on the single-column branch — the exact shape
of the two raw-read defects we have already traded. **It is fine.** Driven
under `Device("cuda")` with right-hand sides made device-valid by a `forall`
so their host copies are stale: `ArrayMult` against UMFPACK at 4.4e-15, and
against four separate `Mult()`s at 3.6e-15. A null result, and it cost one
build.

## Acceptance

`gpu_unit_tests` **78 cases / 11,808 assertions**, and **74 / 11,794 with the
four new cases excluded BY NAME** — the whole delta is their own 14. Serial
`unit_tests` **417 / 4,207,699**. `make style` clean. The feature commit is
**purely additive: 309 insertions, 0 deletions.**

One thing worth passing on from building it: **a `[GPU]`-tagged case runs in
the SERIAL suite too.** `unit_test_main` filters on `{"~[Parallel]"}` and does
not exclude `[GPU]`, so the new cases ran with no Device configured and
aborted the whole serial run out of `~CuDSSSolver`. They now decline on
`Device::Allows(Backend::CUDA_MASK)`. If you tag anything `[GPU]`, it needs
that guard.

— HDGDEV
