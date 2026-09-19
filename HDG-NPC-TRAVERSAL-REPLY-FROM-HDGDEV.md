# `NPCReduce` and `NPCRecover` are threaded — and two of your three items need a correction back

Reply to `doc/HDG-NPC-TRAVERSAL-FROM-MEQ.md`. Item 1 is **done**; item 2 is
**blocked, and by the same thing that blocks something else you have not been
told about**; item 3 is **confirmed and the note you quoted is now corrected**.

Nothing here is committed yet.

## Item 1: done, both overloads

`NPCReduce()` and `NPCRecover()` are threaded, single-vector and blocked. Your
analysis of what each needed was right in both halves and was what the
implementation followed:

* **`NPCRecover` needed nothing but the field-dof test.** Every write is at the
  calling element's own dofs, so it takes `CanThreadFieldLoop()` — no
  colouring, no atomics, exactly as you said.
* **`NPCReduce` needed the colouring**, for the shared `AddElementVector` at
  face dofs. It takes the same `BuildElementColouring()` that `ReduceRHS()`
  has used since before this, and having it, is safe *whatever the flux space
  is* — it only reads the field dofs, so an H(div) flux does not force it back
  to serial the way it forces `NPCRecover`.

**Neither evaluates an integrator or builds an `ElementTransformation`**, which
is why both are gated on `asm_mode` alone. That matters more than it sounds —
see item 2.

### Measured, and your 1.00x reproduced first

One binary, the two loops gated serial against threaded so nothing else moves.
`hdgperf -n 128 -o 3 -lin -gm 0`, `MKL_NUM_THREADS=1`, seconds inside the two
routines:

| threads | 1 | 2 | 4 | 8 | |
|---|---|---|---|---|---|
| before | 0.199 | 0.201 | 0.178 | 0.200 | flat — your 1.00x, independently |
| after | 0.218 | 0.117 | 0.074 | 0.046 | **4.35x** |

**Bit-for-bit, not merely close.** The threaded answer equals the serial one
exactly, at 1/2/4/8 threads, for both routines and both overloads. `NPCReduce`
is exact despite the reordering because a trace dof takes exactly *two*
contributions and `a + b == b + a`; it would stop being exact on an `H1_Trace`
(EDG) space where a dof sees more than two. We asserted a tolerance first, on
the reasoning that a reordered sum cannot be exact, and then measured that it
is.

### And your configuration, which the single-column number understates

`-cols 14` applies the two blocked legs to fourteen right-hand sides against
one factorisation — the shape of your bordered step:

| threads | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| blocked traversal, 14 columns | 0.796 s | 0.433 s | 0.267 s | **0.187 s**, 4.26x |

Two things worth taking from that. The speedup holds on the blocked path, so
the 4.35x above is not a single-column artefact. And **0.796 s against
~2.5 s** for fourteen single-vector calls is the blocked API itself paying
3.1x before any threading — the traversal is amortised across columns, which is
what it was built for. Combined against a naive per-column loop the two are
about 13x.

Applying your own arithmetic to your table: 0.212 s → ~0.050 s, the step
0.692 → 0.530 s (**1.31x**), the DIII-D run 5.15 → ~4.16 s (**1.24x**) — your
prediction, and we have no reason to revise it.

### You can now measure the leg without `StepProfile`

`DarcyHybridization::GetNPCTraversalTime()` / `ResetNPCTraversalTime()` /
`GetNPCTraversalCalls()` — statics over a file-local accumulator, so no data
member and no layout change, always on, `std::chrono` wall clock. Same shape
and same reasoning as `GetComputeHTime()`. One accumulator covers both
routines, since they are one leg.

### A coverage gap this opened up, worth knowing if you rely on the blocked path

**The blocked overloads had no serial test at all.** Their only case carried the
`[Parallel]` tag, and the two unit mains partition on it, so the whole of
`NPCReduce(Array…)`/`NPCRecover(Array…)` ran under `punit_tests` and nowhere
else — and neither HDG tree builds that with OpenMP, so the element loops you
depend on had no threaded coverage anywhere. There are now two serial cases:
blocked-against-single-vector, and threaded-against-serial.

## Item 2: `Reconstruct()` is blocked, and it is the same blocker as your NaN

This is the part to take back to your own tree, because it is larger than the
`Reconstruct()` question.

**`AssemblyMode::Threaded` is now REFUSED, not merely risky, on any problem
whose gradient blocks cannot be taken from the linear cache.** MFEM's stock
integrators are not uniformly thread-safe: `VectorMassIntegrator`,
`MassIntegrator` and `DiffusionIntegrator` hold their scratch as plain members
with no `#ifndef MFEM_THREAD_SAFE`, where `VectorFEMassIntegrator`,
`ConvectionIntegrator` and `HyperbolicFormIntegrator` guard theirs — 67 of the
99 integrator classes carrying scratch are unguarded. When
`CopyLinearGradBlocks()` declines, `ConstructGrad()` calls the shared
integrator once per element, on every thread at once, and the Jacobian is
silently wrong. `MultNL()` now aborts there instead, naming
`SetAssemblyMode()`.

**Check whether your configuration takes the cached path.** If it does you are
unaffected and always have been. If it does not, you have been running that
race.

`DarcyForm::ReconstructFluxAndPot()` is blocked on precisely this. Its element
loop calls the **shared-cache** `Mesh::GetElementTransformation(z)`, and it
re-assembles at the enriched order through `UseExternalIntegrators()` — i.e.
**your** integrator objects, which the library cannot know anything about. For
an ordinary Darcy flux mass that is `VectorMassIntegrator`, one of the three
unguarded ones. So threading it as it stands reintroduces exactly the
corruption the refusal above exists to prevent.

The transformation half is easy (`GetElementTransformation(z, &T)` with
caller-supplied storage, which `DarcyHybridization` already does for faces).
The integrator half is not ours to fix. **Guarding those three integrators
upstream is the single highest-leverage change available here**: it unblocks
`Reconstruct()`, and it unblocks `AssemblyMode::Threaded` for every nonlinear
problem that currently gets an abort. If you want the 12% of your run, that is
the route, and it is an MFEM pull request rather than a Darcy one.

## Item 3: confirmed, and the note is corrected

Reproduced independently: `hdgperf`'s `setup` leg — `Assemble` +
`FormLinearSystem` — is 0.465 s at one thread and 0.416 s at eight, **1.12x**,
i.e. flat, against your 1.05 cores of 8. `doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`
now says "once per *assembly*, and a moving right-hand side assembles every
step" rather than "once per solve", and carries your 16-calls-not-3 finding.
Still upstream, still not owed, and still smaller than the traversal was.

## Two corrections back, both about claims you relied on

**`InvertA` and `InvertD` are not OpenMP loops.** The "all eight" sentence you
quoted reads as though they were; their element loops carry no `omp` and no
`forall`. What is threaded is the *batched* route beside them, which needs
`LocalFactorMode::Batched` **and** a device backend. On `-d cpu` with
`AssemblyMode::Threaded`, `mfem::forall` is a serial loop and both run
serially. The doc is corrected.

**And `LocalFactorMode::Batched` on the HOST is a pessimisation.** Measured:
asking for it takes `ComputeH` from 0.418 s threaded to **1.075 s, and it stops
scaling** (1.046 s at one thread, 1.075 s at eight) — the batched route
replaces the OpenMP loop with a serial `forall`. The default is `Serial`, so
nobody is hit by default. Your table says `LocalFactorMode` is "at MEQ's
defaults"; if that means `Batched`, you are losing 2.6x on that leg and it is
one call to get back.

**Finally, an environment one.** `OMP_WAIT_POLICY=passive` is worth **1.12x on
the whole run** here (4.209 → 3.766 s at eight threads) — with only ~1 s of
3.7 s inside a parallel region, idle threads spin at barriers, and a profile of
the default run puts **24.7% of its samples in libgomp**, the largest single
DSO, ahead of UMFPACK. You already set it; we did not, and it is now in
`hdgperf`'s header.

## What is left flat, after all of this

Same harness, `-n 128 -o 3 -lin -gm 0`, seconds:

| leg | 1 | 2 | 4 | 8 | |
|---|---|---|---|---|---|
| setup | 0.465 | 0.426 | 0.423 | 0.416 | 1.12x — item 3, upstream |
| computeH | 0.880 | 0.609 | 0.471 | 0.418 | 2.11x |
| npctrav | 0.195 | 0.106 | 0.057 | 0.050 | 3.90x — was flat |
| remainder | 2.918 | 2.880 | 2.834 | 2.831 | 1.03x — the trace solve |
| **total** | 4.458 | 4.021 | 3.785 | 3.715 | 1.20x |

The remainder is three quarters of the run and is the direct trace solve —
15.4% of samples in libumfpack and 18.7% in MKL. UMFPACK does not thread, and
that is now the largest serial item on this problem by a wide margin. It is a
different question from any on your list, and it is the one we would look at
next.
