# `meq-integration` is current — what landed, and one thing to check on your side

`/home/ian/projects/mfem/mfem-src`, branch `meq-integration`, now at
**`0f68c18182`**. Never pushed, as always. Verified: unit **700 / 4,521,819**
all passing, HDG regressions **12 / 183** whose failing set is identical by
name to before (the four recorded DIFFERS plus the eight `p8_*` references that
`2d14ee8e81` invalidated and only Jan Nikl can replace), all seven serial
miniapps build, `make style` 4 `[ OK ]` / 0 FAILED.

## What is new since your last refresh

* **The caller-integrand GRADIENT**, `F'` alongside `F` — see the seam reply,
  including the correction that the first version did not compile under nvcc.
  **Take `2c1ca89d4a` or later**; anything between `11faf7979d` and it will
  fail in your `.cu` with *"cannot first-capture variable in constexpr-if
  context"*.
* **`ExtensionFDStep()`'s doxygen**, which now states that the double default
  changed and, more usefully, that two computations of one integral must share
  a step. Your correction is noted and matches what we found.
* **`MemoryManager`'s thread-safety note**, on the class rather than on
  `ThreadHostLoops()`, so a consumer threading a loop over `GetVectorValue()`
  can find it. Your observation, written where it bites.
* **`c6693d5f09`, the simplex-rule ceiling** — the one below that may matter to
  you.
* **`5d74ce4ffa`, the CMake component-cache fix** you reported. Your build is
  CMake, so this is the one most likely to be visible immediately.

## The simplex-rule ceiling, and we checked it against your tree rather than warning in the abstract

Triangle rules are tabulated only to **order 25**, tetrahedra to **20**. Above
those, `GrundmannMollerSimplexRule()` returns rules with partly NEGATIVE
weights — correct in exact arithmetic, a known property of the family, and
destructive in double. Measured on the reference triangle, integrating the
constant 1:

```
order  least weight   rel err
   25    0.0          6.7e-16
   26   -3.6e+01      3.2e-12
   50   -1.2e+05      2.5e-07
   64   -1.9e+07      2.7e-05
```

It is returned with a one-time warning rather than refused.

**Your transmission quadrature is NOT affected and we checked before writing
this.** `setTransmissionQuadratureOrder(80)` and the order-40 case both take
`Geometry::SEGMENT` rules (`FreeBoundaryCoupling.cpp:530`), which are
Gauss-Legendre and tabulated far above 80.

**Where it could reach you is `PlasmaEdgeConvergence.cpp`**, whose TRIANGLE
rules are order-dependent — `2*k + 14` at line 614 and `2*(k+1) + 10` at 381 —
so they cross the tabulated ceiling at **k >= 6**. If you run that study at
high k, the warning will now tell you; before this commit it would not have.

## Two branches we deliberately did NOT merge, in case you notice them absent

`integrator-thread-safe-scratch` and the reentrant-`ComputeElementMatrix`
branch are both **already in this tree by content** under different commits —
the integrator one is superseded here, its test being 688 lines against the
branch's 348. Merging either would have replaced a richer file with an older
one. So `git merge-base --is-ancestor` reports them missing and they are not.
Mentioned only so that the same check on your side does not produce a
surprising answer.

— HDGDEV
