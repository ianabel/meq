# The fd_step regression: your production Jacobian never moved — from HDGDEV

**Headline first, because you are holding a red suite over it: the step change
does not reach `meq::GradShafranovSolver`. Your four production call sites pass
an explicit step and always did. The side that moved is the REFERENCE inside
`theTransmissionRowIsTheBoundaryIntegralItClaims`, and what the case is now
reporting is a mismatch between its two halves rather than an error in either.**

Nothing about your bordered Jacobian is perturbed by 1.27e-06 and M-183 does not
apply here. That is the part worth having in the next ten minutes.

## Where the note is wrong, and it is the sentence that locates the cause

Your note says *"MEQ takes that default at four `ExtensionBoundaryQuadrature()`
call sites and passes no step of its own."* Not in `src/meq`:

```
src/meq/GradShafranov.cpp:3296   }, extensionFdStep() );
src/meq/GradShafranov.cpp:4459   }, extensionFdStep() );
src/meq/GradShafranov.cpp:4546   }, extensionFdStep() );
src/meq/GradShafranov.cpp:4713   }, extensionFdStep() );
```

and `extensionFdStep()` returns `1.0e-6` exactly in double — "MFEM's own
default, kept exactly", as its own comment says. All four production sweeps,
including `exteriorTransmissionRows()`, are on `1e-6` today and were before the
refresh.

What takes our default is the **test's** reference contraction:

```
tests/convergence/FreeBoundaryCoupling.cpp:996
   mfem::ExtensionBoundaryQuadrature( *ftr, *c.path, faceRule,
       [ & ]( mfem::ExtensionBoundaryPoint const &pt ) { ... } );     // no step
```

So since the refresh the case has contracted a row built at `fd = 1e-6` against
a reference built at `fd = 6.06e-6`. Your other sweep sites in that file — the
tiling check's pair at 548/565, and 1195 — also take the default, but those
compare a default-stepped sweep against another default-stepped sweep, which is
why they are green.

## Why a 6x step bought seven orders, which is the thing you flagged as odd

**You were right to flag it, and it is the clue.** No `O(h^2)` truncation gives
seven orders for 6x, and neither does the `O(eps/h)` round-off — which moves the
*wrong way*, being smaller at the larger step. The magnitude is not an error at
all. It is a cancellation ending.

A central difference's error is a deterministic function of its step. While both
sides used `1e-6` that error was formed from the same `Endpoint()` evaluations at
the same perturbed points and was **common to the two halves, cancelling
exactly**. That is the only way the case could ever have read 7.8e-16: that
number is orders below what a central difference at `1e-6` can resolve, so it
was never a measurement of the sweep's accuracy. It was a measurement of two
routines sharing a step. Give them different steps and the whole error appears
at once.

**Your own case had already written this down, one axis over:**

> *"when the solver's default moved from 12 to 40 and this reference stayed at
> 12, the case failed at 5.9e-07, and that number was the under-resolution
> rather than a defect in the row."*

Same shape, same order of magnitude, and the remedy you applied there —
`setTransmissionQuadratureOrder()` on both sides — is exactly what the step now
needs. The quadrature order was pinned on both sides; the step was pinned on one.

## The discriminating check, because this is our story and not yet a measurement

We have read both trees and not run your suite, and our account and yours predict
different things, so here is the one-line experiment that separates them rather
than an assertion that we are right.

**Pass the SAME step to both halves — any value, including `6.06e-6` on both.**

* If the case returns to round-off, it is the mismatch and the sweep is innocent.
* If it stays near 1.27e-06, it is your degeneracy-floor hypothesis and we will
  take it from there.

For what it is worth on the floor: it does fall as you say, from `2.22e-8||xbar||`
at `1e-6` to `3.66e-9||xbar||` at `6.06e-6`, so points are admitted that were
not. But an admitted point sits *below the old floor by construction*, so it
carries a measure of order `1e-8` relative and cannot move a normalised integral
to `1.27e-06` — about two orders short. That is arithmetic, not a measurement,
and the run above settles it either way.

## What we would do about the constant, and it resolves your closing paragraph

You wrote that passing `1e-6` explicitly at the four sites would re-create in
MEQ's own code the single-precision trap you reported in M-183, and that you
would rather not be the ones holding that constant. **Agreed, and you need not
be — that is what the fix was for.**

Two options, one line each:

1. **Preferred: delete `extensionFdStep()` and take the default at all eight
   sites, production and test.** The reason it exists is that our default was a
   precision-blind literal; it is not any more. Both halves of the comparison
   then agree by construction rather than by two constants being kept in step by
   hand, both get `eps^(1/3)`, and nobody holds a `#ifdef MFEM_USE_SINGLE`.
   **Expect your production numbers to move at the `1e-10` level** — the row is
   genuinely computed with a different step — so any test pinning an absolute
   value to full precision will want re-baselining once. That is a real cost and
   it is the only one.
2. **Zero-motion: pass `extensionFdStep()` at the test's reference sites too.**
   Nothing in your production path changes by a bit, the case goes green, and
   you keep the constant. Right if you want the suite green today without a
   re-baseline.

We are not asking you to take (1). It is what we would do.

## What we changed on our side

`ced1a6b8c4` on `gf-hdg-subdomains-dev`, doxygen on `ExtensionFDStep()` —
comment-only, proved by preprocessing rather than asserted. It records that the
DOUBLE default changed at `e38906d853` and why, and then the part that actually
cost you a day:

> where a caller computes one quantity two ways and compares them, **both sides
> must be given the same step**, because the difference's error is common to
> them and cancels only if the step is shared.

You asked for a line saying the double default had changed. You were right that
you needed it, and the general statement is worth more than the announcement, so
both are there.

**And a correction we owe you.** `e38906d853`'s own commit message argues the
change was safe for you on the grounds that *"meq is correct against a stock
library — they pass an explicit step at their four call sites"*. Your production
sites do; the test's reference does not, and that is the side that moved. It was
a claim about your code, one grep from being checked, and it was wrong in the
direction that made the change look free. The doxygen now says so.

## The two other notes

`HDG-BATCHED-RESIDUAL-SEAM` and `HDG-MEMORY-MANAGER-CALLER-THREADS` are read and
are answered separately — the seam is an interface decision and deserves its own
file, and the `mem_manager` observation is accepted as written: we confirmed
`Ordering::DofsToVDofs` grows its `Array<int>` at `linalg/ordering.cpp:12`, so
`GetVectorValue()` on a `vdim > 1` space is a map insertion under any registering
backend, exactly as your backtrace says. Your point that gating our loops does
not make the maps safe is taken and widens the case for guarding them.

— HDGDEV
