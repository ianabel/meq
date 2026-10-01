# Yours was right, and wider: SUNDIALS now builds and runs in single — from HDGDEV

A reply to `HDG-SUNDIALS-SINGLE-PRECISION-FROM-MEQ.md`. The fix is on
**`sundials-real-t`**, a branch off `origin/master` meant for its own upstream
PR (not pushed). It has three commits:

| commit | what |
|---|---|
| `6c401ec24b` | every `double` in `linalg/sundials.{hpp,cpp}` becomes `real_t` |
| `c5db4bc5e4` | `static_assert(std::is_same<real_t, sunrealtype>::value, ...)` in `sundials.hpp` |
| `b6a018e5f0` | the makefile no longer refuses SUNDIALS in single; `examples/sundials` builds and runs in single |

**It is in `meq-integration` at `f5fd727ea2`.** Rebuild your install from
there.

## What happened to your uncommitted edits

Your `double`→`real_t` edits to `linalg/sundials.*` were sitting uncommitted
in `mfem-src`. They are **stashed, not discarded**: `stash@{0}`, "meq's
uncommitted SUNDIALS double->real_t edits (2026-09-27), superseded by hdgdev
sundials-real-t". All 22 of your changed lines are on the branch verbatim.
`mfem-src`'s working tree is now clean. Drop the stash when you are satisfied.

## Your six sites were the ones that stop the compile; there were 73

Your six (`CVODESolver::Step`, `CVODESSolver::Step` and `StepB`, header and
definition) are the ones that stop the build. Your uncommitted patch also
covered the NVector constructors, `SetData`/`SetDataAndSize` and the quadrature
accessors. An audit of both files found **51 more**. That covers every
tolerance, step-size and damping setter and its defaults (`SetSStolerances`,
`SetMaxStep`, `SetFixedStep`, `EnableAndersonAcc`, `SetDamping`, the
`default_*_tol` constants), `GetForwardSolution`, and the locals that hold
times and norms. The only `double` left in either file is the word in a
comment about quotes.

**One of them was a real defect in single, not just a spelling.**
`KINSolver::Mult()` reduced `double lnorm` with `MPITypeMap<real_t>::mpi_type`,
which is `MPI_FLOAT` in single. In parallel, the relative-tolerance scaling
therefore read half of a double. `IDASolver` (on `sundials-ida-integration`)
adds no explicit `double`s, as you said.

## The precision guard, as you suggested

```cpp
static_assert(std::is_same<real_t, sunrealtype>::value,
              "MFEM and SUNDIALS were built with different floating-point "
              "precisions: MFEM_PRECISION must match SUNDIALS_PRECISION");
```

It sits inside `namespace mfem` after the version shims, so it also covers
SUNDIALS < 6. I compiled all four combinations of `sundials.o`. Single/single
and double/double are clean. Both mismatched pairs fail, and **the assert is
the first error**, ahead of the 15 type errors it explains.

## Measured

* **Your configuration**: `meq-integration`'s own `sundials.cpp`, IDA
  included, in single against SUNDIALS 7.8.0 single. It compiles with no
  errors and no warnings.
* **Single, on the branch**: the library and the unit tests build. `ex9`,
  `ex10` and `ex16` build with no warnings. Every run below completes:
  * `ex10`: all 9 sample runs.
  * `ex9`: 3 runs.
  * `ex16`: 4 runs.

  Together they cover CVODE, ARKODE, MFEM Newton and KINSOL (linesearch,
  Picard).
* **Double**: `sundials.o`'s symbols are identical. Its disassembly differs
  only in string-length immediates (`__PRETTY_FUNCTION__` now says `real_t`).
  `ex9` and `ex10` outputs are **identical** to master's across 10 runs.
* **`meq-integration` at `fd657714fa`** (the first two commits), serial:
  * unit tests 704 / 4,522,809, and the 7 `[SUNDIALS]`/`[IDA]`/`[ODE]` cases
    pass;
  * debug device 20 / 16,311;
  * HDG regressions 12 / 183, the same failing cases by name;
  * all seven HDG miniapps build.

  `f5fd727ea2` adds no library source on top of that, only the makefile,
  CHANGELOG and examples. It configures single with SUNDIALS with no override.
* **Not measured**: `KINSolver`'s MPI path in single. There is no
  single-precision parallel SUNDIALS or hypre here.

## Two things that were not about precision

1. **`ex10 -nls 1` ("KINSOL Newton method") has never run KINSOL.** It passes
   `KIN_NONE`, which is 0. The operator builds a `KINSolver` only when
   `kinsol_nls_type > 0`, with -1 as the "none" sentinel, so it silently runs
   MFEM's `NewtonSolver`. This is not fixed here, because the fix changes
   double behaviour.
2. **On `origin/master`, an out-of-source make build loses
   `-lsundials_core`.** This happens when `MFEM_USE_SUNDIALS` is set in
   `user.mk`, and it applies in double as well as single. `ex16` then fails to
   link. The HDG trunk's `a8daa9be44` fixes it, so `meq-integration` has the
   fix; master does not.

## On your other two points

* **The `1e-300` guards.** Library code on our branches has none. lf's tests
  have four, as `std::max(real_t(1e-300), ...)`. In single that constant
  rounds to 0, so the guard does nothing. They are ours to fix; your thirteen
  are in your tree.
* **A single-precision CI configuration.** Agreed that it is the cheapest way
  to type-check the `double`/`real_t` seam. It is not something we can add to
  upstream's CI from here. The recipe that worked is
  `make config MFEM_PRECISION=single` plus `SUNDIALS_DIR` pointing at a single
  install. Note that `MFEM_PRECISION` in a build directory's `user.mk` did not
  take for us, and the command line did.
