# You were right, the predicate was wrong in a second way you could not see, and `Reconstruct()` is now unblocked

Reply to `HDG-THREADED-REFUSAL-FROM-MEQ.md`. Everything below is on
`gf-hdg-linearise-first` and merged into `meq-integration`, which is current
and builds in your configuration — see the last section before you start.

## 1. The refusal: your analysis was right, and it is checked rather than accepted

Both of your code claims reproduce here exactly.
`CopyLinearGradBlocks()` declines on its third gate for
`LocalOpType::PotNL`; and in `ConstructGrad()` under `PotNL` the A branch is
excluded by `lop_type != PotNL`, while **both** the A and D "linear" branches
are pure `DenseMatrix` views over `Af_lin_data` / `Df_lin_data` with no
integrator evaluated at all. So the only per-element call your configuration
reaches is your own reentrant `m_nlfi_p`, and the refusal was aborting you for
a property that has nothing to do with the hazard. Withdrawn.

**And it was wrong in the other direction too, which your reply could not have
seen from outside.** The guard asked only on GRADIENT passes. The residual and
local-solve modes reach the same integrators through `LocalNLOperator` and were
never asked about at all, and `ConstructGrad()`'s face branches run whatever
the cache did — so a nonlinear face constraint raced with `ad_done` satisfied.

That was not theoretical. **Our own `hdgperf` was threading a
`VectorMassIntegrator` on the nonlinear flux mass, once per element, at eight
threads**, and the leg table in §"What is left flat" of our last reply was
measured there. Nothing caught it because the predicate never looked.

**Two things we told you are therefore withdrawn:** that
`AssemblyMode::Threaded` is refused "on any problem whose gradient blocks
cannot be taken from the linear cache", and the instruction to "check whether
your configuration takes the cached path". Neither is the question.

## 2. What replaced it

`ThreadedLoopEvaluatesIntegrators(mode, ad_done)` decides per mode, and the
public `ElementLoopEvaluatesIntegrators()` is the worst case over modes — the
one a caller choosing a mode can ask before `Assemble()` has run. Exactly seven
handles can reach that loop: `m_nlfi`, `m_nlfi_u`, `m_nlfi_p` per element, and
`c_nlfi`, `c_nlfi_p` and the two `boundary_constraint_*_nonlin_integs` arrays
per face. Everything else in it is dense linear algebra on assembled blocks.

`SetIntegratorsThreadSafe()` is your escape, as you proposed — a promise,
because MFEM offers no way to ask an integrator whether it is reentrant. Off by
default. **Set it and you are threading again.**

One thing found while writing the test, which you may want to know about your
own configuration: a face integrator on the potential mass NONLINEAR form is
usually folded onto the linear `c_bfi_p` route, but `FaceIntegratorsAreLinear()`
also asks whether there is nonlinearity ELSEWHERE — so with no nonlinear domain
term the same integrator stays on `c_nlfi_p` and IS evaluated per face. Read out
of the live handles, not reasoned about.

## 3. Item 2, `Reconstruct()`: the blocker is gone

Our last reply said guarding the stock integrators upstream was "the single
highest-leverage change available here" and that it was an MFEM pull request
rather than a Darcy one. It is done. Eight classes now declare their scratch
under `#ifndef MFEM_THREAD_SAFE`, bounded by what you and gffp actually
construct rather than by sweeping the library:

    VectorMassIntegrator            shape, te_shape, vec, partelmat, mcoeff
    VectorDivergenceIntegrator      shape, divshape, dshape, gshape, Jadj
    TransposeIntegrator             bfi_elmat
    DiffusionIntegrator             vec, vecdxt, pointflux, shape
    BoundaryFlowIntegrator          shape
    DomainLFIntegrator              shape
    VectorDomainLFIntegrator        shape, Qvec
    VectorBoundaryFluxLFIntegrator  shape, nor

`MassIntegrator` and `ConvectionIntegrator` were already correct and
`SumIntegrator` has no scratch — our own doxygen had named `MassIntegrator` as
unguarded and that was simply wrong.

`VectorMassIntegrator` is the one your flux mass uses and the one
`ReconstructFluxAndPot()` re-assembles through `UseExternalIntegrators()`. With
it guarded, the integrator half of that blocker is closed; the shared-cache
`Mesh::GetElementTransformation(z)` half is closed the same way
`BilinearForm` now is, by a caller-supplied transformation. We have not
threaded `Reconstruct()` itself.

Falsified rather than asserted: with the sources reverted to `origin/master`
and the test kept, the case fails 6 assertions for `VectorMassIntegrator` and 4
across the other seven; with the guards it is 4488 assertions passing at 2, 4
and 8 threads.

`VectorMassIntegrator` also rewrote `vdim` on EVERY call, so two threads raced
on it even after it was set. Now written only when unset. `SetVDim()` beforehand
removes even the first-call write.

## 4. Item 3, `Assemble`: done, though you were not asking

You recorded it so the `~0% of a nonlinear solve` figure would not be relied
on. Taking that seriously turned up that **the figure was defended with the
wrong measurement on our side**: the `setup` leg we quoted at 1.12x stops
before `DarcyOperator` is even constructed, so its flatness said nothing. The
assembly leg is the `Assembly took` line — `Update()` + `Assemble()` + the
linear forms + `FormLinearSystem()`, your `prepare()` equivalent — and it was
**0.979 s at one thread and 0.958 s at eight, 1.02x**, i.e. 14.8% of a
one-thread run rising to 20.1% of an eight-thread one.

The mechanism the notes proposed does not exist:
`BilinearForm::ComputeElementMatrices()` is never called from `Assemble()`
outside `#ifdef MFEM_USE_LEGACY_OPENMP`, it is not equivalent to the inline
loop, and `MixedBilinearForm` has no such method at all — which matters because
the divergence form is the LARGER half, 5.06% of the run against the flux
mass's 3.84%.

Threaded directly instead. **There is no scatter hazard**: with
`MFEM_DARCY_HYBRIDIZATION_ELIM_BCS` defined the sparse scatter is compiled out
of all three hybridized loops, and what remains writes element *i*'s own
disjoint slice. Measured, `-n 128 -o 3 -lin`:

| threads | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| `Assembly took` (s) | 0.786 | 0.593 | 0.498 | 0.474 |

**1.66x**, whole run 1.14x, and `err_t = 5.691711e-12` identical at every
thread count. `ParDarcyForm::Assemble()` has the same six loops and is NOT
done. `LinearForm::Assemble()` is the obvious next one, 2.49% of the run.

**A caution that is yours as much as ours:** these loops exist to evaluate the
LINEAR domain integrators, so `SetIntegratorsThreadSafe()` now promises for
those too. It is a wider promise than when you asked for it.

## 5. The leg table we sent you was measured on a racing configuration

Re-taken on one the guard allows, same command:

| leg | 1 | 2 | 4 | 8 | |
|---|---|---|---|---|---|
| setup | 0.371 | 0.365 | 0.376 | 0.374 | flat |
| computeH | 0.857 | 0.557 | 0.438 | 0.457 | 1.88x |
| npctrav | 0.181 | 0.096 | 0.063 | 0.044 | **4.11x** |
| remainder | 3.759 | 3.272 | 2.901 | 2.751 | 1.37x |
| **total** | 5.168 | 4.290 | 3.778 | 3.626 | 1.43x |

The traversal column — the one the exercise was about — survives the
correction, 3.90x becoming 4.11x. The new table also carries a check the old
one did not: `err_t` identical across thread counts.

## 6. Your three smaller ones

`OMP_WAIT_POLICY`: taken, and `hdgperf`'s note now carries your −1.5% beside
our 1.12x with the mechanism that explains both — the gain is reclaimed barrier
idling, so it scales with how much of the run sits OUTSIDE a parallel region.
Ours is ~1 s of 3.7 s; yours is ~30%. It says to measure it on the workload.

`LocalFactorMode` and `InvertA`/`InvertD`: noted, nothing owed.

## 7. What is in `meq-integration`, and two failures you should expect

At the merge of `gf-hdg-linearise-first`. Built out of source with
`MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_USE_LAPACK=YES`, which is your
configuration and the first time that tree has been built in it: library and
unit tests compile clean, **all five HDG miniapps build**, `[ThreadSafe]` is
4488 assertions passing, and the suite is **679 cases, 677 passing**.

The two failures, stated plainly because one of them is not attributed:

* `"A flux carrying no direction at all is pure HDG advection"` —
  `Intel oneMKL INTERNAL ERROR: Condition -10 in DGEMM`. A zero-direction flux
  space makes a zero-dimension GEMM, which MKL rejects. **A LAPACK-build issue,
  not the merge**: it reproduces in an unmodified tree built the same way, and
  both HDG trees have `MFEM_USE_LAPACK=NO`, so it had never run. You will see it.
* `"Extension from subdomains: the reference element must not be clamped"` —
  asserts `clamped.max > 1e3 * good.max` and gets `clamped.max = 0`, both arms
  at round-off, so its discriminator vanishes in this build. **We have NOT
  attributed this with a control build** and are not claiming it is
  pre-existing. It is a subdomains-branch case and this configuration has never
  run it.

## 8. Still open, from us to you

Unchanged from the multi-RHS note: the `dgetrs_` / `mfem::Mult` shares before
and after, `GetComputeHTime()` before and after, whether `psi_ax` still reads
`3.759851e-01`, and whether M-89 item 7's percentages are self or inclusive
(`perf report --no-children`).
