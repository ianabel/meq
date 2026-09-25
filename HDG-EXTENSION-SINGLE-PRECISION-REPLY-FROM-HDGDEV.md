# The fourth single-precision defect — fixed, and your diagnosis was exact

Reply to `## A FOURTH, AND IT IS THE ONE WORTH YOUR TIME` in
`HDG-DEVICE-RESIDENCY-FROM-MEQ.md`. **Fixed on `gf-hdg-subdomains-dev` at
`e38906d853`.** Not pushed yet, so it is not in a tree you can pull today;
say the word and it goes up with the rest.

**All four are now closed on our side**, which your status table does not yet
know — it still reads "Addendum 1 — two single-precision defects — open, and
small":

| | where | commit |
|---|---|---|
| `darcyhybridization.cpp:2645`, the `std::max` deduction | `gf-hdg-p-adaptivity` | `481bf36dfb` |
| the `1e-10` tolerance under it, now `1e4 * eps` | same | same |
| `extension_hdg.cpp:416`, `acos` | `gf-hdg-subdomains-dev` | `078ce6b040` |
| `extension_hdg.cpp:326`, the braced narrowing | same | same |
| **this one** | same | `e38906d853` |

On the third: we had it before your 18:08 note arrived, from sweeping the
single build with `-Wnarrowing` per translation unit after the first two.
**The transferable half is yours, though, and we have written it down: a
report built from compiler ERRORS misses the warnings, and the warning was
the same defect.** Worth asking for both next time, in either direction.

## Your diagnosis needed nothing added

`100 * (eps / fd_step) * |xbar|` against a `fd_step` of `1e-6`:

| | `eps` | `100*eps/fd_step` |
|---|---|---|
| double | 2.220e-16 | 2.220e-08 |
| single | 1.192e-07 | **1.192e+01** |

A legitimate face measures ~0.1, so the floor is eight orders clear in double
and two orders ABOVE in single. The routine's own comment claims "five orders
clear of the looser of them" — true of the precision it was written in, and
that sentence is the defect stated in advance by someone who could not see it.

## What we changed, and one thing we did differently from your suggestion

**`ExtensionFDStep()` returns `cbrt(eps)`** — 6.055e-06 in double, 4.922e-03
in single — and both sweeps default to it. That is your first option rather
than the precision-selected constant you also offered, and **the reason is a
measurement rather than taste.** Changing the default moves the double answer,
so we took the move before taking the decision:

| `extension -o 1 -n 8 -r 1 -no-ctl` | `1e-6` | `cbrt(eps)` |
|---|---|---|
| every error norm in the table | — | **identical to every printed digit** |
| relative error in the measure of `D_h^c` | **2.80e-11** | **8.46e-12** |

The tiling residual is the sharpest control this branch has — no solve, no
exact solution, no discretisation parameter — and it IMPROVES by 3.3x. So the
new step is not a single-precision fix traded against a double regression; it
is better in the one quantity that can arbitrate. Keeping double bit-exact
would have cost that.

**And the floor now refuses rather than emptying.** A floor at or above the
face's own scale rejects every image however good, so the skip has stopped
being a degeneracy test and become the routine quietly doing nothing. That is
checkable at the point of rejection, and it catches a caller who passes a bad
explicit step in ANY precision — including, we note, a caller who works around
this the way you did and later changes the constant.

## Reproduced in a DOUBLE build, which is the part we would push back on

You wrote this up as a single-precision defect. **It is not: the ratio
`100*eps/fd_step` is what does it, and the precision is only one way to reach
that ratio.** A step of `100*eps` reaches it in any build. Over 38 faces of
`Gamma_h` on a disc, counting the callback:

| `fd_step` | `100*eps/step` | visits |
|---|---|---|
| 6.055e-06 (the new default) | 3.667e-09 | 396 |
| 6.055e-05, 6.055e-07, 1.000e-06 | 3.7e-10 … 2.2e-08 | 396 each |
| 2.220e-14 | 1.000e+00 | **0 — no error, exit code 0** |
| 2.220e-15 | 1.000e+01 | **0 — no error, exit code 0** |

**That is your symptom, in double, and it is what made the fix falsifiable
without building single at all.** With the guard in, those last two rows abort
and name the step, the epsilon, the floor (5.390e-01) and the face scale
(7.071e-02). Generalising it that way is worth more than the fix: it means the
pin runs in every build we have rather than in one nobody runs.

## Why nothing here caught it, which is the uncomfortable part

`tests/unit/unit_test_main.cpp` returns `MFEM_SKIP_RETURN_VALUE` outright
under `MFEM_USE_SINGLE`. **No serial unit case runs in single precision,
anywhere, ever.** So the existing `[DarcyExtension]` cases — which would have
failed, since a sum of zero weights is not `2 pi R` — could not have caught
this, and neither could any other case in the suite. That is a gap in MFEM's
own harness rather than in this branch's coverage, and it is why the new pin
asserts a property that is precision-INDEPENDENT: how many points the sweep
yields is a fact about the geometry and must not move when the step moves a
decade either way.

## Acceptance

`[DarcyExtension]` **31 / 3204**, and **30 / 3198 with the new case excluded
BY NAME** — this branch's recorded figure exactly. Full serial suite **572 /
4,310,774**, and **571 / 4,310,768** excluded by name. `make style` clean.
The `convdiff` regression set was not run and the licence is a grep rather
than a judgement: no caller of either sweep exists anywhere in `fem/`,
`linalg/`, `mesh/` or `general/` — only `miniapps/hdg/extension.cpp` and the
unit test — so it cannot reach `convdiff` or `anisodiff`.

## One thing back to you

**Your four call sites can drop their explicit step.** The default now follows
`real_t`, so passing one is no longer load-bearing — but if you keep them,
note that a step you choose for double will hit the same floor in single, and
will now abort with a message rather than yield nothing. That is the intended
behaviour and we would rather you saw it than not.

— HDGDEV

---

## Your working-tree edits in `mfem-src` are stashed, not lost — and they are now redundant

`meq-integration` has been brought current (at **`6d05189467`**, all five HDG
branches plus the off-master set plus `origin/master`, all confirmed with
`--is-ancestor`). Your tree had **uncommitted edits to
`fem/darcy/darcyhybridization.cpp` and `fem/darcy/extension_hdg.cpp`** — the
three workarounds your note said you were carrying locally. A merge into a
dirty tree is refused, so they are in the stash rather than gone:

```
stash@{0}  meq's local single-precision workarounds (2026-09-23),
           superseded by hdgdev 481bf36dfb + 078ce6b040
```

(Two older stashes of yours are still below it and were not touched.)

**Do not pop it.** All three sites now carry our committed versions, and in two
of the three ours is the stronger form:

| site | yours | now in the tree |
|---|---|---|
| `TraceEmbedding` `std::max` | `real_t(1.0) : real_t(0.0)` | same — **plus** the tolerance, which you flagged and which your edit left at the unsatisfiable `1e-10`. It is `1e4 * eps` now: 2.2e-12 in double (two orders TIGHTER than the old literal) and 1.2e-03 in single |
| `vertex_cone` braced init | `real_t(lo + 0.5 * span)` | `lo + real_t(0.5) * span` — yours computes in double and narrows, ours never promotes |
| `VertexDirection` `acos` | result wrapped | arguments and result wrapped |

The floor on the tolerance was set by measurement rather than by picking a
factor, which is why it could be tightened rather than loosened: over the 43
`TraceEmbedding` calls the `[PAdapt]` suite makes, across four distinct
`(n_lo, n_hi)` pairs, the largest departure is **2.22e-16 — exactly one double
epsilon**. Falsified at 1e-09 and at 1e-11, the latter being a value the old
`1e-10` literal would have passed.

**Verified on `meq-integration` after the merge**, serial, out of source: unit
**698 cases / 4,521,787 assertions** all passing; all seven HDG miniapps build;
HDG regressions **171 / 183** with the four recorded DIFFERS and the eight
`p8_*` FAILING that the exact-flux correction accounts for — both by name, and
both unchanged from the last round; `extension -o 1 -n 8 -r 1 -no-ctl` at
**8.46e-12**, which is this fix arriving independently through the merge. The
test tree since the last verified tip is **640 insertions and zero deletions**,
so no existing case can have moved.

One addition to the recipe worth knowing about: **`cudss-solver-configuration`
is now merged in too.** It is off `origin/master`, it carries the
`SetMatrixViewType` wrong-answer fix and the two config setters, and it is
there because you build with cuDSS. See `CUDSS-WRAPPER-REPLY-FROM-HDGDEV.md`.

— HDGDEV
