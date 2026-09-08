# MEQ against freegs4e

**An independent check by a different algorithm**, which is the thing
`FREE-BOUNDARY-PLAN.md` FB-6 wants and a fine-mesh self-comparison cannot be.
`../freegs4e` is FreeGS-derived: free boundary by von Hagenow Green's functions,
2nd/4th-order finite differences on a uniform `(R,Z)` grid, Picard with adaptive
blending. MEQ is HDG with Newton. They share the equation and essentially no
code.

**Result, seven configurations** — TestTokamak (three profile shapes), MAST,
MAST-U, TCV, DIII-D — each seeded from its own reference:

| case | rel L2 in ψ | rel L∞ |
|---|---|---|
| A TestTokamak, classic | **1.5e-04** | 5.4e-04 |
| G MAST-U | 6.4e-04 | 1.6e-03 |
| B TestTokamak, ff′-dominated | 7.1e-04 | 1.5e-03 |
| E TestTokamak, diamagnetic | 2.7e-03 | 5.8e-03 |
| F DIII-D | 2.9e-03 | 5.2e-03 |
| D TCV | 4.4e-03 | 7.9e-03 |
| C MAST | 6.8e-03 | 1.3e-02 |

**Two independent codes over ~11,000 nodes, worst case 0.7%.** The residual is
almost certainly geometric rather than either solver's: the MXH fit of the
boundary is 2–4e-04 m on a minor radius of 0.24–0.61 m, which is 5e-04 to 1e-03
relative, and the contour it fits is itself extracted from freegs4e's 129² grid.

## What it does

MEQ solves fixed boundary; freegs4e solves free boundary. So the comparison is:
freegs4e converges an equilibrium, and MEQ is then given **its own interior flux
surface as a boundary** and **its own profiles as a source**, and has to
reproduce the field inside.

1. `fgsref.py` — seven converged free-boundary equilibria, saved as `.npz`.
2. `surface.py` — extracts a closed `ψ_n = 0.9` contour.
3. `mxh.py` — fits it to the Miller eXtended Harmonic form MEQ's
   `[boundary.shape] Type = "mxh"` takes.
4. `convert.py` — freegs4e's profiles onto MEQ's `ψ`.
5. `make_case.py` — writes the TOML and the two profile tables.
6. `compare.py` — reads MEQ's `.nc` and diffs. **It refuses a run whose
   `ψ_ax` is not a magnetic axis**; see below.

```sh
python3 -m venv venv
venv/bin/pip install numpy scipy matplotlib h5py shapely numba netCDF4 contourpy
PYTHONPATH=/home/ian/projects/freegs4e venv/bin/python fgsref.py
venv/bin/python make_case.py out/*.npz
for f in out/*.toml; do meq "$f"; done
venv/bin/python compare.py out/*.npz
```

**Do not `pip install -e` freegs4e** — its `requirements.txt` pins `numpy<2`,
`numba~=0.60`, `Shapely~=2.0.6`, none of which have wheels for Python 3.14, and
Shapely then fails to build on a missing `geos_c.h`. Relaxed, `numpy 2.5.2` /
`shapely 2.1.2` / `numba 0.67` work. `numba` looks optional and is not:
`critical.py`'s fallback calls `warnings.warn` without importing `warnings`.

## Four traps, each of which produces a plausible wrong answer

**freegs4e's profiles are `dp/dψ`, NOT `dp/dψ_n`.** Its own docstring says
`dp/dpsi_n` and is **wrong** — established four ways, decisively by rebuilding
`Jtor` from the saved arrays and matching the solver's to 4e-16. Dividing by
`ψ_bnd − ψ_ax` on the strength of the docstring costs a factor that for a
tokamak is order 1 and does not look wrong. Only the *derivative* column carries
a chain-rule factor, because MEQ's `SplineProfile` is a Hermite cubic.

**Not the separatrix.** Every reference case is diverted, so its LCFS has an
X-point **corner**, and MXH is a truncated Fourier series that cannot represent
one: fitting it reads 2.6e-03 to 1.3e-02 m where the fitter does 1.9e-05 on a
smooth shape. MEQ's own `ExtensionConvergence` avoids the separatrix for exactly
this reason. An interior surface fits to 2–4e-04 m, which is then the floor of
the contour extracted from a 129² grid rather than of the fitter — it stops
improving past 10 harmonics, and that flatness is what identifies the floor.

**ψ = 0 must land on the surface actually used.** Mapping the profile with
`ψ = 0` at `ψ_n = 1` while handing MEQ the `ψ_n = 0.9` surface puts MEQ's
boundary where the source has nearly died — `p'` of 2.7 instead of 3.1e5 — and
the solve converges cleanly, in three Newton steps, to a field 23× too small.

**THE RAMP PICKS THE BRANCH, AND IT IS PART OF THE PROBLEM STATEMENT.** These
sources do not vanish at `ψ = 0`, so this is not the trivial branch; it is the
second solution beside it. Measured on case A, the amplitude sweep finds exactly
two roots — 1.93e-03 and 5.23e-02 — and lands on the weak one for any amplitude
below about 2× the axis height. `[initialguess] Type = "ramp"` at **6×** the
expected axis value clears it on all seven.

## compare.py refuses a run whose psi_ax is not a magnetic axis

**`ψ_ax` is the largest NODAL value of `ψ_h`, and nothing in that definition
makes it a magnetic axis.** The bordered Newton's constraint
`ψ_ax − max ψ_h = 0` is satisfied at machine zero by a spurious nodal value
exactly as it is by an axis, so nothing the solve reports can tell them apart —
and `ψ_ax` is what the profiles are normalised by, so a wrong one is not a bad
number, it is **a different equilibrium**. Differencing that against the
reference compares two different problems, at whatever precision.

So the driver writes `axis_normalised_flux` into the `.nc` — `Ψ` at the located
O-point of `q_h`, which **must be 1** — and `compare.py` reads it and **exits
non-zero without quoting any norm** when it falls below `0.90`. That threshold
is `meq::CriticalPointFinder::checkAxis`'s own default, and the test is
**one sided**: a healthy field approaches 1 from *above*, because `ψ_ax` is a
nodal value and a polynomial's peak over a closed element is at least that.

`--allow-bad-axis` prints the norms anyway, flagged `BAD AXIS`, for deliberately
inspecting a known-bad run.

**An absent attribute warns rather than refusing**, because three different
things produce one: a run that was not normalised (no `psi_axis` either, so
there is nothing to judge), the wall-hugging **annulus** branch where no O-point
exists anywhere, and a `.nc` written before the attribute existed. The last two
are indistinguishable from the file alone, and an old file must not become
unreadable.

## It is psi-star that is compared, not psi_h

**The `.nc` carries the POST-PROCESSED potential.** `apps/meq.cpp` samples
`solver->postProcessedPotential()`, and says so in its own comment: *"THE
POTENTIAL SAMPLED HERE IS psi*, NOT psi_h"*. `_psi.gf` keeps `psi_h`, which is
what a restart reads back — but `compare.py` reads the `.nc`, so every number on
this page is `psi*`.

**That is the right field to compare and it changes how the table reads.**
`psi*` comes from `DarcyForm::Reconstruct()` and converges at **k+2**, not k+1,
so `k = 3` on 397 elements is delivering a degree-4 field converging at 5. It is
also what a consumer of MEQ actually gets — the `.nc` is the interchange format
— so comparing it is honest rather than flattering. But attributing the
accuracy-per-dof result to `k+1` convergence would be wrong.

**The convergence tests in `tests/convergence/FreeBoundaryCoupling.cpp` use
`solver.potential()`, which is `psi_h`**, and their `k+1` rates are that field's.
The two are deliberately different measurements: those ask what the
discretisation does, this asks what a user receives.

If a comparison against `psi_h` is ever wanted, it needs the `.gf` route rather
than the `.nc` — which is a second reason to reach for pyMFEM.

## Performance, and the only column that means anything

`perf.py`. **The two codes do not solve the same problem**, so a wall-clock
ratio is not a statement about either: freegs4e converges a *free*-boundary
equilibrium — coil Green's functions, X-point finding, a control system, 23 to
103 Picard steps — while MEQ solves the *fixed*-boundary problem inside a
surface it is handed. MEQ is doing strictly less work, in C++ against Python.

**So the column that means something is accuracy per unknown**, which is
language-neutral and is the real question: does high-order hybridized DG reach a
given accuracy with fewer degrees of freedom than a 2nd/4th-order finite
difference grid? Case A, `MKL_NUM_THREADS=1`, `OMP_NUM_THREADS=1`:

| `k` | refine | elements | dofs | dofs / ref pt | wall | rel L2 |
|---|---|---|---|---|---|---|
| 1 | 2 | 1735 | 15,615 | 0.94 | 5.25 s | 3.15e-04 |
| 1 | 3 | 7185 | 64,665 | 3.89 | 13.58 s | 1.46e-04 |
| 2 | 1 | 397 | 7,146 | 0.43 | 2.69 s | 1.01e-03 |
| 2 | 2 | 1735 | 31,230 | 1.88 | 6.01 s | 1.52e-04 |
| **3** | **1** | **397** | **11,910** | **0.72** | **3.00 s** | **1.72e-04** |
| 3 | 2 | 1735 | 52,050 | 3.13 | 6.45 s | 1.44e-04 |

freegs4e's own reference: **84.7 s at 129² = 16,641 grid points**.

**High order is worth about 5× in unknowns.** `k = 3` on 397 elements reaches
1.7e-04 with 11,910 dofs and 3.0 s; `k = 1` needs 7,185 elements, 64,665 dofs
and 13.6 s for the same accuracy — 5.4× the unknowns and 4.5× the time.

**MEQ SATURATES THIS BENCHMARK, WHICH IS THE MOST USEFUL THING IN THE TABLE.**
The errors stop improving at 1.44e-04 for case A — and at 5e-03 and 2.4e-03 for
MAST and DIII-D — because that is the *reference's* own accuracy, set by the MXH
fit and the contour extracted from its 129² grid. Refining MEQ past the second
row buys nothing. Any future work on this benchmark should refine the
**reference**, not MEQ.

Timings include the whole driver — mesh generation, solve, post-processing and
all three output formats. The `.nc` sampling is ≤0.5 s of it, measured by
rerunning at a 17² output grid, so the wall clock is solve-dominated. Read it as
an order of magnitude: a timing here is a measurement about this machine.

**One failure worth knowing.** `k = 2, refine = 1` does not converge on MAST or
DIII-D — Newton fails, the reactive ladder's `PicardThenNewton` fails too, and
the driver says to raise the degree. `k = 3` on the *same mesh* then converges in
8 Newton steps. That is `p`-refinement curing what neither `h` nor a
globalisation reaches, which is the pedestal finding in `CLAUDE.md` met from a
new direction.

## Root selection is the whole difficulty, and case A proves it

**These problems have at least three solutions and only one is the equilibrium.**
Case A, measured:

| start | ψ max | vs reference |
|---|---|---|
| ramp, amplitude ≤ 0.1 | 1.931300e-02 → 1.9313e-03 | −95.7% |
| **seeded from the reference** | **4.525154e-02** | **−0.067%** |
| ramp, amplitude ≥ 0.2 | 5.225179e-02 | +15.4% |
| freegs4e | 4.528205e-02 | — |

The physical root sits **between** the two a ramp can reach, which is the
signature of the unstable middle branch of an S-curve: Newton slides off it from
either side. No ramp amplitude finds it — the sweep goes straight from the lower
root to the upper one between 0.1 and 0.2.

**Seeded, MEQ converges to it in TWO Newton iterations** — residual
4.5e-03 → 2.7e-09 → 7.1e-15 — and then agrees best of all seven cases. So the
reference is a genuine root of MEQ's discrete problem, MEQ finds it immediately
from a nearby start, and the earlier 15% was never a disagreement about the
equation.

**That the upper root was FLAT under refinement is what proved it was a root and
not an error**: +15.39% → +15.41% across refinement levels 2→4 and degrees 2→3,
converged to six digits. Discretisation error falls; a different solution does
not.

`mkguess.py` writes the reference as an MFEM mesh plus `H1_2D_P1` GridFunction —
MEQ's `[initialguess] Type = "gridfunction"` interpolates it through
`meq::FieldTransfer`, so the guess needs neither MEQ's mesh nor its degree.
**pyMFEM would be the better tool here** and is worth reaching for next time:
it would replace this file's hand-rolled ASCII writer, and would let `compare.py`
read MEQ's `.gf`/`.mesh` directly instead of going through the lossy `.nc` grid.

**The source conversion was ruled out before any of this**, directly rather than
by argument: evaluating MEQ's tables on the reference's own ψ reproduces
freegs4e's `μ₀ R J_φ` to **2.3e-05**. `CLAUDE.md` records the same multiplicity
from inside MEQ — three solve routes reaching discrete solutions 9.4% apart on
an under-resolved mesh.

## Refining the reference, 2026-09-06 — and what actually limits it

This file has said since it was written that **MEQ saturates this benchmark**:
its error stops falling at about 1.4e-04 because that is the *reference's* own
accuracy, and that the next work here is to refine the reference rather than
MEQ. That is now measured, and it is right — with numbers, and with three
findings about `freegs4e` that were not expected.

`fgsref.py` grew three flags for it: `--nx=N` (every case on an `N x N` grid,
rounded up to `2^k + 1` so the ladder nests point for point), `--hagenow` and
`--vcycle=L`. `ladder.sh` drives the sweep and records wall clock, peak resident
memory and the Picard count at each level.

### Refining the reference works, and it works through the boundary fit

Case A, the reference's grid refined and MEQ rerun against it unchanged:

| reference grid | MXH shape fit | MEQ rel `L2` | MEQ rel `Linf` |
|---|---|---|---|
| 129² | 2.060e-04 m | **1.519e-04** | 5.403e-04 |
| 257² | 1.063e-04 m | **6.751e-05** | 3.140e-04 |
| | **1.94x** | **2.25x** | 1.72x |

**MEQ's error tracks the shape fit**, 2.25x against 1.94x. So the saturation was
never MEQ's: the LCFS handed to MEQ is a fit to a contour extracted from the
reference's grid, the contour improves like `h`, and the fit and the answer
improve with it. Extrapolated, 513² is about 3e-05 and 2049² about 8e-06 —
**an order of magnitude, not round-off.** Grid refinement alone cannot reach
round-off here because the fit improves only linearly and is in the budget at
all *because this is the fixed-boundary rehearsal*. **FB-6 proper removes it**:
a free-boundary MEQ takes the same coils and profiles and never sees an LCFS.

### Caching the default boundary condition, which is the better lever

`boundary.freeBoundary` **rebuilds a geometry-only matrix on every Picard
step**. `Greens` depends only on the grid and `romb` is linear, so the whole
boundary condition is one fixed matrix applied to `Jtor`:

```
psi_bndry = M Jtor,      M[ b, ij ] = dR dZ w_i w_j G_b[ i, j ]
```

with `w` the Romberg weights and `G_b[ b ] = 0`, exactly as the loop zeroes its
self term. Building it costs about one pass of the loop; every step after that
is a matrix–vector product.

| `n` | loop / step | build once | matvec | boundary values agree to | speed-up |
|---|---|---|---|---|---|
| 65 | 126.7 ms | 98.6 ms | 2.34 ms | 4.09e-16 | 54× |
| 129 | 651.5 ms | 518.2 ms | **8.07 ms** | 5.46e-16 | **81×** |
| 257 | 4271.0 ms | 4234.9 ms | **24.02 ms** | 5.46e-16 | **178×** |

**End to end on `H_limited_circular`: 4.9 s against 29.2 s at 129², and 13.3 s
against 175.1 s at 257²** — 6.0× and **13.2×** on the whole run, the ratio
growing because the loop is `O(n³)` in elliptic-integral evaluations while the
matvec is one BLAS call.

**THE WEIGHTS COME FROM `scipy`, NOT FROM A REIMPLEMENTATION.** `romb` is linear,
so `romb( I )` **is** the weight vector: one call on an identity, and the rule
can never drift from the one the shipped loop uses. Writing Romberg out here
would be a second definition of the quadrature when the point is that there is
only one. Two details reproduced rather than tidied: the **corners appear twice**
in `freegs4e`'s boundary index list (harmless — both writes put the same value in
the same cell — but a different index set is a different matrix), and `M` is
cached **on the Equilibrium object**, so it cannot outlive the grid it was built
for. A case-keyed global would be one edit away from a stale matrix, which is a
wrong answer rather than a slow one.

**MEMORY IS THE CONSTRAINT.** `M` is `2( nx + ny ) × nx·ny` doubles: 69 MB at
129², 0.54 GB at 257², **4.3 GB at 513²**, 35 GB at 1025². Above
`CACHE_BOUNDARY_MAX_GB` the shipped loop runs instead, so too large a grid is
slow rather than an OOM.

#### The converged answers differ by 4.5e-08, and that is `freegs4e`'s floor

Not the reassociation, and establishing that took three tries — two hypotheses
were wrong and the third was the control that should have come first.

* **The boundary condition itself is exact.** On the loop run's own converged
  `Jtor`, the two agree to **6.8e-16** at 129² and 5.6e-16 at 257².
* **Not the stopping rule.** Tightening the Picard tolerance from 1e-9 to 1e-11
  to 1e-13 — 41, 60, then 401 iterations, the last hitting `maxits` — leaves the
  gap at **4.511e-08 to four digits, unchanged**.
* **Not amplification either.** Perturbing `M` deliberately by 1e-12, 1e-10 and
  1e-8 relative moves the converged field by 3.39e-08, 6.50e-08 and 7.71e-08:
  **four orders in the input, a factor of 2.3 in the output.** A floor, not a
  gain.
* **And the control that decides it**: each variant is **bit-reproducible against
  itself**, `0.000e+00` run to run, so the difference is real and deterministic
  rather than noise.

So `freegs4e`'s converged answer on this case is determinate to about **5e-08**
and no better, whatever the arithmetic does — a floor on any comparison against
it. It sits five orders below that grid's own discretisation error (2.0e-02 at
129², 7.9e-04 at 513²), so it does not matter for use; it matters for what may
be *claimed*. The cache is exact where exactness is checkable, and the run-level
agreement is the reference's own.

### `freegs4e`'s cost is its boundary condition, not its solve

| grid | default boundary | von Hagenow | peak RSS (Hagenow) |
|---|---|---|---|
| 129² | ~85 s | **29 s** | 293 MB |
| 257² | **656 s** | **71 s** | 503 MB |
| 513² | (hours) | **230 s** | 1564 MB |

`Equilibrium.__init__` takes `boundary=freeBoundary`, whose own docstring calls
it *"an integral over the area of the domain for each point"*: it loops the `4n`
boundary points and evaluates `Greens` over the whole `n²` grid for each. That
is **`O(n³)` per Picard step**, and it is the 7.7x for a fourfold rise in
unknowns. `boundary.freeBoundaryHagenow` is in the same file, is the method this
benchmark's own description credits `freegs4e` with, is **`O(n²)`**, and is not
the default.

**The linear solve is not the cost, and the multigrid is switched off anyway.**
`Equilibrium.__init__` calls `multigrid.createVcycle( ..., nlevels=1, ... )`,
and at one level a V-cycle is a direct sparse solve on the full grid. Two traps
beside it: `setSolverVcycle()` hard-codes `GSsparse`, the **second**-order
generator, and ignores `Equilibrium.order` — so calling it on a 4th-order
equilibrium silently solves a different problem; and a V-cycle built correctly
on the 4th-order generator **does not converge at all**, failing with
`ValueError: No opoints found!`. The multigrid path is second order only.

### THE TWO BOUNDARY CONDITIONS DISAGREED, FLAT UNDER REFINEMENT — AND THE CAUSE IS ONE CONSTANT

**RESOLVED 2026-09-07: the shipped von Hagenow branch is INCONSISTENT, and it
is usable once the observation-point offset is tied to the cell.** The section
below is what was measured before that, and it reads the flat gap as an open
question about which boundary condition to believe. It is not open: a boundary
condition whose disagreement with an exact integral does not shrink is wrong,
and this one is wrong for a reason that takes one line to fix. See
**Making von Hagenow usable** below.

Case A, the same equilibrium computed both ways:

| grid | rel. difference in `ψ_ax` | rel. difference in the field |
|---|---|---|
| 129² | 6.460e-04 | 3.201e-03 |
| 257² | 6.404e-04 | 3.147e-03 |

**It does not shrink.** A discretisation difference falls by 4x or 16x per
level; this falls by 1.7 per cent. Self-convergence says which one is right:

| | 129² → 257² | 257² → 513² |
|---|---|---|
| default boundary | **2.78e-05** | — |
| von Hagenow | 3.64e-04 | 1.22e-04 |

The default is **converged at 129²** on this case — `ψ_ax` stable to six figures
— while von Hagenow creeps toward a different value.

**AND "TOWARD A DIFFERENT VALUE" WAS THE ANSWER ALL ALONG, NOT A CURIOSITY.**
Two consistent discretisations of one problem must converge together; a gap that
does not shrink means one of them is not a discretisation of that problem at all.
That is what the next section establishes, and it makes the fast boundary
condition usable rather than merely suspect.

### Making von Hagenow usable

`boundary.freeBoundaryHagenow` displaces its observation point off the boundary
by a **hard-coded `eps = 1e-2` metres**, *"to avoid the singularity in `G(R,R')`
when `R'=R`"*. That is a fixed **physical** distance: it does not shrink with the
grid, so the `O(eps·∂ψ/∂n)` error it costs is a constant and the method converges
to the wrong answer.

Measured on one fixed `Jtor` with **no solve in the way** — so the direct
Green's integral is ground truth and what is printed is von Hagenow's own error,
on a `0.1 … 2.0` by `−1 … 1` box:

| `eps` | order in `h` at `n` = 65, 129, 257 |
|---|---|
| `1e-2`, as shipped | 0.63, 0.21, **0.07** — stalling |
| **`0.2 h`** | 0.90, **1.00, 0.99** — clean first order |
| `0.5 h` | 0.97, 0.99, 0.99 — the same, 1.8× larger |

**THE OPTIMUM IS NOT "AS SMALL AS POSSIBLE", WHICH IS WHY THE CONSTANT IS 0.2.**
Swept at a fixed 129² grid the gap reads 7.6e-03, **3.3e-03**, 8.9e-03, 1.6e-02,
2.3e-02 at `eps` = 1e-2, 3e-3, 1e-3, 3e-4, 1e-4. Below about `0.2 h` the log
spike is **narrower than the cell** and the Romberg rule misses it; above it the
displacement dominates. Tying `eps` to `h` is what keeps the two terms on the
same footing.

**AND ON THE REAL SOLVE**, `H_limited_circular`, against the direct integral:

| variant | `n` | rel `ψ_ax` | rel `ψ_bnd` | rel L2 field | rate |
|---|---|---|---|---|---|
| shipped | 129 | 5.310e-03 | 1.797e-02 | 1.008e-02 | — |
| shipped | 257 | 5.233e-03 | 1.836e-02 | 1.032e-02 | **−0.04** |
| **scaled** | 129 | 1.087e-03 | 3.679e-03 | 2.139e-03 | — |
| **scaled** | 257 | **5.378e-04** | **1.886e-03** | **1.093e-03** | **0.97** |

The shipped one gets slightly **worse** under refinement. The scaled one halves,
in all three measures at once.

**WHAT IT IS GOOD FOR, AND THE HONEST CEILING.** First order is as far as a
displaced observation point goes, so extrapolating: at 513² the boundary error is
about **2.7e-04** in `ψ_ax`, against that grid's own discretisation error of
**7.9e-04** on this case. So it does not dominate, and 513² costs **230 s**
instead of hours. That is the whole point — it makes the reference refinable.

**THE OBVIOUS HIGH-ORDER CURE DOES NOT WORK, AND MEASURING THAT IS WHAT SAYS
WHERE THE ERROR IS.** `∮G σ dl` is a **single-layer** potential — continuous
across the boundary, with its normal derivative jumping — so averaging the
evaluations at `+eps` and `−eps` ought to cancel the `O(eps)` term exactly. It
does not: 1.5e-02, 8.1e-03, 4.0e-03, 2.0e-03 against the one-sided 1.2e-02,
6.5e-03, 3.3e-03, 1.6e-03 — the same rate of 1.00 and slightly worse. So the
residual is the **quadrature of the near-singular kernel** rather than the
displacement, and getting past first order needs the log subtracted and
integrated in closed form. That is a boundary-element method and is deliberately
not attempted.

**`freegs4e` IS NOT EDITED.** The change is one constant, lifted out and tied to
the cell, reinstalled from `fgsref.py` as WORKAROUND 5 beside the four already
there. `--hagenow` gives the usable version; `--hagenow-eps=shipped` reproduces
`freegs4e` exactly, so both readings are available on one problem; and the
default boundary path is bit-unchanged — the 129² run still reports
`psi_axis 0.0948314088`.

### The converged 513² reference, and the baseline it sets

**Run 2026-09-07 with the exact boundary condition cached** — 129² → 257² →
513², each rung seeded from the one below, coincident points reproducing to
9.44e-16. `baseline.json` is the machine-readable record and `baseline.py`
regenerates it.

| `n` | points | `ψ_ax` | `ψ_bnd` | Picard | wall | peak RSS |
|---|---|---|---|---|---|---|
| 129 | 16,641 | 0.09483140885 | 0.02781828740 | 41 | 4.85 s | 340 MB |
| 257 | 66,049 | 0.09337971414 | 0.02649111261 | 31 | 12.84 s | 1000 MB |
| **513** | **263,169** | **0.09308752051** | **0.02622461849** | **27** | **72.75 s** | **5663 MB** |

**Successive differences fall by 4.97, so the order in `h` is 2.313** — slower
than the fourth-order generator because a LIMITED boundary is a maximum over the
limiter ring, a pointwise operation on a discrete set, where a diverted one is a
saddle located by interpolation. Richardson gives

```
psi_ax  -> 0.093013888     513^2 is 7.92e-04 out
psi_bnd -> 0.026157662     513^2 is 2.56e-03 out
```

**THAT 7.9e-04 IS THE FLOOR ON ANY COMPARISON AGAINST THIS REFERENCE.** A MEQ
result closer to the 513² numbers than that is measuring the reference's own
grid error, not MEQ. Beneath it sits `freegs4e`'s ~5e-08 indeterminacy, which is
four orders lower and therefore never the binding constraint.

**AND 72.75 s IS WHAT THE EXACT BOUNDARY CONDITION NOW COSTS AT 513²**, of which
**37.4 s is building the 4.02 GB matrix once**. For scale, the table below
estimates the shipped loop at *hours* for that rung and von Hagenow at 230 s —
so the cache is **3× faster than the approximate method and exact**, which is
why von Hagenow is not needed for refinement after all.

**The Richardson sign was got wrong on the first pass and CLAUDE.md caught it.**
The sequence DECREASES toward its limit, so the extrapolate is BELOW the finest
rung; adding the correction instead of subtracting it put it above, at
0.093161 rather than 0.093014. The tree's own independently recorded value is
what refuted it — which is the argument for writing a number down the first time
it is measured.

### How far the grid can go, and what sets it

**Memory, and this machine has 23 GB with 21 available.** Peak resident, measured:

| | 129² | 257² | 513² | 1025² |
|---|---|---|---|---|
| exact boundary, cached | 340 MB | 1.0 GB | **5.7 GB** | 35 GB — out |
| von Hagenow | 293 MB | 503 MB | 1.56 GB | ~6 GB |

**So 513² is the rung that is both exact and fast** — 72.75 s with the matrix
cached — and it is the one stored. **1025² is reachable for the RUN but not for
the cached matrix**: `CACHE_BOUNDARY_MAX_GB` falls back to the shipped loop
there, which makes it slow rather than impossible. 2049² is out on both counts,
and with the faithful boundary condition would additionally be about ninety
hours.

**Re-check `free` rather than trusting this table.** The WSL2 allocation was
raised on 2026-09-07 and every ceiling recorded before that was a claim about a
smaller machine. A memory ceiling is a property of the day.

**And the route past 513² is FB-6 rather than a finer grid** — a free-boundary
MEQ never sees an LCFS, so the fit that limits this comparison stops existing.


## The floors, peeled one at a time — 2026-09-06

The refinement above answered "does refining the reference help" with *yes, 2.25x
from 129² to 257²*. Carrying it to 513² and then chasing what was left turned
into a sequence of floors, each one only visible once the previous was removed.
**The end of it is a factor of 17.**

| what was changed | MXH fit | MEQ rel `L2` | what was limiting |
|---|---|---|---|
| reference 129², 10 harmonics, `k=2 r=2` | 2.060e-04 m | 1.519e-04 | the contour |
| reference 257² | 1.063e-04 m | 6.751e-05 | the contour |
| reference 513² | **1.032e-04 m** | 5.726e-05 | **the FITTER — 3 % for 4x the grid** |
| 16 harmonics | **2.206e-05 m** | 5.825e-05 | **MEQ — `L2` did not move at all** |
| MEQ `k=3` | 2.206e-05 m | **8.802e-06** | ? |
| MEQ `k=3 r=3` | | 8.744e-06 | not MEQ |
| output grid 257², 513² | | 8.797e-06, 8.800e-06 | **not the sampling either** |

**Each row is the previous row's diagnosis being wrong.** The fit stopped
improving at 513², so the floor was the fitter and not the contour — this file
previously said the opposite, and it was right *at 129²*, where ten harmonics
happen to be enough to reach the contour's own limit. Raising the harmonics then
bought a 4.7x better boundary and **`L2` did not move**, which is what identified
MEQ's own discretisation as the binding constraint — the first time in this
benchmark's history that has been true. Refining MEQ bought 6.7x. And the last
two rows are the controls: `k=3 r=3` moves the answer by 0.7 %, and a sixteenfold
denser comparison grid moves it by 0.03 %, so neither MEQ nor the sampling is
what is left.

**`Linf` behaves differently and says the same thing.** It *did* follow the
boundary fit — 2.900e-04 to 1.347e-04 when the harmonics went up, a factor of
2.15, while `L2` sat still. A boundary that is wrong by a shape error shows up
worst near the boundary, which is what a maximum norm sees and an `L2` over the
whole core does not.

**So the remaining 8.8e-06 is the boundary fit**, at 2.2e-05 m against a minor
radius of 0.42 m, and more harmonics do not help — the fit plateaus at 16, so it
is contour-limited again and the contour wants a finer reference than this
machine can hold.

**WHICH IS THE ANSWER TO "DO THE CODES AGREE TO ROUND-OFF".** They agree to
**8.8e-06** and cannot be made to agree better *through this comparison*,
because the comparison contains a boundary fit and the fit is now what limits
it. Round-off is not reachable here by construction, and no amount of grid is
going to change that. **FB-6 removes the fit** — a free-boundary MEQ takes the
same coils and the same profiles and never sees an LCFS — and that is the only
route to the question the refinement was asked to answer.


## Cost, measured on a quiet machine — 2026-09-06

Everything below is **single-threaded on both sides**, `MKL_NUM_THREADS=1
OMP_NUM_THREADS=1`, taken with nothing else running and re-taken after an
earlier set was polluted by another build on the same box. The load factor was
small for this workload -- `freegs4e` is single-threaded numpy, and its 129²
run went 85 s busy against 81.8 s quiet -- but the numbers here are the quiet
ones.

### The raw comparison, which is the unfair one

| | wall | reaches |
|---|---|---|
| MEQ, `k=3`, refine 2 | **7.72 s** | 8.80e-06 |
| `freegs4e` 129² | **81.8 s** | — |
| `freegs4e` 257² | **567 s** | — |
| `freegs4e` 513² | **~5300 s** | the reference |

**A ratio here is not a statement about either code and should not be quoted as
one.** `freegs4e` converges a FREE-boundary equilibrium -- coil Green's
functions, X-point finding, a control system, and the `O(n³)`-per-Picard-step
boundary condition documented above -- while MEQ solves the FIXED-boundary
problem inside a surface `freegs4e` had to find for it. MEQ is doing strictly
less work, in C++ against Python. The 690x against the 513² run is mostly a
statement about `freeBoundary` being an area integral per boundary point.

### At comparable accuracy, which is the one worth having

`freegs4e`'s own error curve comes from its self-convergence against its 513²
answer on the nested grid:

| | wall | rel `L2` vs its own 513² |
|---|---|---|
| `freegs4e` 129² | 81.8 s | 2.699e-05 |
| `freegs4e` 257² | 567 s | 2.642e-05 |

**IT IS FLAT**, which is the third time a number in this study has failed to
move when it should have. Sixteen times the unknowns buys two per cent, so
`freegs4e`'s answer stops converging in the grid at about **2.6e-05** relative
and the rest is something else -- the Picard tolerance, the LCFS determination,
or the boundary condition it does not agree with itself about. Note this is a
whole-grid RMS while MEQ's column drops the band, so the two are not the same
norm and should not be differenced; they are comparable as orders of magnitude.

Against MEQ's ladder at the same accuracy:

| target | MEQ | `freegs4e` |
|---|---|---|
| ~2.7e-05 | `k=1` refine 3, **17.1 s** | 129², **81.8 s** |
| ~2.7e-05 | `k=3` refine 2, **7.72 s**, and it overshoots to 8.8e-06 | |

so **roughly 5x to 11x at matched accuracy, single-threaded**, and the honest
reading of that is *not* "MEQ is 11x faster": it is that a high-order hybridized
DG reaches this accuracy on 11,302 comparison nodes where a 2nd/4th-order finite
difference grid needs 129² and a boundary condition that costs `O(n³)`. The
discretisation is doing the work; the language and the problem statement are
doing some of the rest.

### MEQ's own ladder

| degree | refine | wall | rel `L2` | rel `L∞` |
|---|---|---|---|---|
| 1 | 2 | 5.93 s | 2.7341e-04 | 5.447e-04 |
| 1 | 3 | 17.14 s | 2.2805e-05 | 7.183e-05 |
| 2 | 1 | 3.22 s | 1.0300e-03 | 1.954e-03 |
| 2 | 2 | 7.73 s | 5.8247e-05 | 1.347e-04 |
| 2 | 3 | 17.03 s | 9.4131e-06 | 3.532e-05 |
| 3 | 1 | 3.32 s | 1.0644e-04 | 2.182e-04 |
| **3** | **2** | **7.72 s** | **8.8020e-06** | 3.221e-05 |
| 3 | 3 | 19.55 s | 8.7442e-06 | 3.427e-05 |
| 4 | 2 | 12.03 s | 8.6725e-06 | 3.460e-05 |

**`k=3` refine 2 is the point to quote**: same wall clock as `k=2` refine 2 to
within one per cent, and 6.6x the accuracy. Going further buys nothing --
`k=3` refine 3 costs 2.5x the time for 0.7 %, and `k=4` costs 1.6x for 1.5 % --
because 8.7e-06 is the boundary fit and not MEQ.

## Refining the reference: run it NESTED, not cold

`fgsref.py --nx=N` keeps freegs4e's `2^n + 1` convention so that 129, 257, 513
and 1025 nest **point for point**, and `--seed-from` is what spends that: it
reads the coarse `plasma_psi`, lifts it onto the fine grid with a cubic
`RectBivariateSpline`, and hands it to Picard as its starting iterate.

**What that costs.** The dominant per-step term is the free-boundary condition,
measured at **0.512 s at 129², 3.540 s at 257²** — about 6.9x per doubling, so
roughly 24 s at 513² and 170 s at 1025². A hundred cold steps at 1025² is
therefore hours, and a 1025² run was in fact killed at **12h36m** for exactly
this reason.

**The recipe, coarsest first:**

```sh
export FGSREF_OUT=$SCRATCH/fgsref
python3 fgsref.py --nx=129                     # the coarse rung, cold
python3 fgsref.py --nx=257 --seed-from=auto    # from n129
python3 fgsref.py --nx=513 --seed-from=auto    # from n257
```

`--seed-from=auto` derives the coarse directory one doubling down by the same
naming rule `--nx` uses for its output, so a scan is a loop over `--nx` with
nothing else to keep in step. `--seed-from=DIR` names one explicitly.

**IT IS WORTH 8 PICARD STEPS, NOT AN ORDER OF MAGNITUDE, AND THIS SECTION
PREDICTED THE WRONG THING.** Measured on `H_limited_circular`, 129² → 257²:

| | Picard steps | wall | `ψ_ax` |
|---|---|---|---|
| cold | 39 | 152.7 s | 9.337971063593e-02 |
| **seeded from 129²** | **31** | **126.1 s** | 9.337971354137e-02 |

**1.21x, and the saving is a FIXED NUMBER OF STEPS rather than a fraction.**
The mechanism is that Picard's contraction here is geometric at about **0.60 a
step**, and a *converged* coarse answer is still **1.6 % wrong** on the fine grid
— that is the reference's own discretisation error, the thing the refinement
exists to measure. So the seed starts the run a factor of ~62 closer than cold
does and buys `log(62)/log(1/0.60) ≈ 8` steps at the top, and nothing after.
Eight steps out of 39 is 20 %; out of a longer run it is proportionally less.

**So it does not rescue 1025².** The earlier claim here — "a cold 1025² run is
hours where a seeded one should be a handful of steps" — was wrong in the way
worth recording: it assumed the cost of a Picard run is set by *where it starts*,
when it is set by *how far it has to go* and by a contraction the seed does not
change. A 1025² run seeded from 513² should still take about three quarters of
its cold step count.

**And the answers are NOT bit-identical: they agree to 3.1e-08 relative in
`ψ_ax` and 9.3e-08 in `ψ_bnd`.** That is far below the 1.6 % the grid refinement
is measuring, so the seed does not move the benchmark — but it is worth saying
plainly rather than claiming an identity that is not there. Picard converges to
the fine grid's own solution whatever it starts from; what differs at the eighth
digit is where its `rtol = 1e-9` stopping test happened to bite.

**The nesting itself is asserted, not assumed.** `seed_from_coarse()` refuses a
coarse grid that does not divide the fine one, a seed finer than the run, and a
coarse file whose extents differ — the last of which would otherwise interpolate,
converge, and describe the wrong machine. A *missing* seed is a warning and the
run continues cold, which is the right way round: seeding is an optimisation, so
its absence should cost time and never correctness. On the shared points the lift
reproduces the coarse data to **7.5e-16** of `|ψ|max`, which is the check that
the strides are right rather than a check on the spline.

**Do not go past 513².** The reference is not the binding constraint on the
benchmark — the **MXH boundary fit** is, and it improves only linearly in `h`,
so a finer reference buys progressively less. 1025² also carries about 6 GB,
which on this machine competes with whatever else is building.

**And do not reach for the two solver knobs first.** `--vcycle=` is a lever on
the wrong term: `multigrid.MGDirect` factorises **once** at construction and only
backsolves afterwards, so the linear solve is not a per-step cost. `--hagenow`
is genuinely `O(n²)` where the default is `O(n³)`, but freegs4e's two boundary
conditions **disagree by 3.2e-03 with a flat gap**, and self-convergence says the
default is the converged one — so it changes the answer at the level MEQ is being
compared at. Seeding costs nothing and changes nothing.


## The FREE-boundary case, and it is not driven from this directory

**Everything above is the fixed-boundary rehearsal**: MEQ is handed freegs4e's
own interior flux surface as a Dirichlet boundary and asked to reproduce the
field inside. That comparison contains an MXH fit and therefore cannot reach
round-off however fine either grid is.

**`H_limited_circular` is the case that removes the fit**, and MEQ solves it as a
FREE-boundary problem: same four coils, same profile shape, same prescribed
`I_p`, same limiter point, and **no LCFS handed over at all**. It ships as
`examples/limited-tokamak.toml` and its regression is
`theDriverSolvesALimitedTokamak` in `tests/convergence/DriverAcceptance.cpp`, so
it is driven from `examples/` rather than from here — but the reference, the
profile conversion and the initial guess are all this directory's.

| | freegs4e, 129² | MEQ, `k = 3`, 1601 elements | apart |
|---|---|---|---|
| `ψ_ax` | 9.483141e-02 | **9.484400e-02** | **1.33e-04** |
| `ψ_bnd` | 2.781829e-02 | **2.782291e-02** | **1.66e-04** |
| profile amplitude | 1, by construction | **1.001239** | 1.24e-03 |
| `I_p` | 3.0e+05 A | 3.000000e+05 A | it is the constraint |

**THOSE NUMBERS MOVED BY A FACTOR OF TWENTY ON 2026-09-07 AND THE MESH DID NOT
CHANGE.** They read 9.455354e-02, 2.774057e-02 and 0.994505 — 2.9e-03, 2.8e-03
and 5.5e-03 — while `ψ_bnd` was pinned at the potential dof NEAREST the limiter
point rather than at the point. Snapping is an `O( h )` error at every degree,
not the `O( h^{k+1} )` nodal choice it was described as; see
`FREE-BOUNDARY-PLAN.md` §7.20 and `GradShafranovSolver::LimiterConstraint`.

**AND MEQ CONVERGES IN ITS OWN MESH ON THIS CASE**, which could not be shown
before, because the boundary condition was quantised by the mesh. Uniform
refinement of the shipped `.msh`, `RefinementLevels = 0, 1, 2`:

| elements | 1601 | 6527 | 26375 | order | limit | vs the reference |
|---|---|---|---|---|---|---|
| `ψ_ax` | 9.484400e-02 | 9.482652e-02 | 9.482423e-02 | **2.93** | 9.482388e-02 | **7.9e-05** |

The regression ships the coarse mesh: it is on the same branch at a fraction of
the cost and is now within 1.4e-04 of that limit.

**WHY THE LIMITER IS PRESCRIBED AS A POINT AND NOT AS A CURVE.**
`FreeGSProfileMixin.attach_limiter` builds the innermost layer of grid CELLS
inside the wall polygon and takes `ψ_bndry` to be the maximum over that layer.
The wall here is `R0_LIM, A_LIM = 1.00, 0.35`, and on this case that maximum sits
at **( 1.3375, −0.0125 )** — one cell in `R` and one in `Z` inside the circle, `dR
= dZ = 0.0125`. Interpolating the saved `ψ` there returns **2.781829e-02**, the
reported `ψ_bndry` to every digit, which is the direct confirmation rather than
an inference from the code.

**So the reference's boundary flux carries an `O( h )` error and its plasma does
not touch the limiter anywhere.** The maximum of `ψ` on the TRUE circle is
**2.574498e-02 at θ = 121.8°, ( 0.8157, 0.2976 )** — on the INBOARD side, not
beside the reported contact at all — which is **7.5%** below the `ψ_bndry` the
run reports.

So the comparison hands BOTH codes that same point, as MEQ's
`[boundary.limiter] R = 1.3375, Z = 0`. It is a well-posed problem both codes
solve identically and it takes the contact-finding out of the comparison
entirely.

**THIS PARAGRAPH USED TO END "it also caps what the agreement can mean: MEQ pins
`ψ_bnd` at the nearest potential DOF … and is not done". IT IS DONE, 2026-09-07,
AND THE CAP WAS THE WHOLE ERROR BUDGET.** MEQ evaluates `ψ_h` AT the point now.
The size of what was being tolerated: sweeping the requested `R` with one key
changed, the entire solve is **bit-identical over `[ 1.3250, 1.3500 ]`** — a
plateau 0.025 m wide, 7% of the minor radius — and steps 5% in `ψ_ax` and 11% in
`ψ_bnd` at each end.

**AND THE CONTACT ITSELF IS A GRID ARTEFACT THAT MOVES TO THE OTHER SIDE OF THE
MACHINE.** At 513² the same ring maximum sits at **( 0.8250, −0.3000 )**, the
inboard shoulder, **0.6 m** from the 129² contact — and the maximum on the TRUE
circle is on the inboard shoulder, `θ = 121.8°`, on *both* grids. So the 129² run
reports a contact on the wrong side of the machine, and a comparison against a
finer reference must take the limiter point, **and the coil currents**, from that
same run.

**THE COIL CURRENTS ARE AN OUTPUT OF freegs4e's CONTROL SYSTEM, NOT AN INPUT**,
so they converge in the grid like everything else: P1 moves 1.4% between 129² and
513² and **P2 moves 14%**, −46186 to −52740. A MEQ run carrying the currents in
`examples/limited-tokamak.toml` and compared against `ref-n513/` is comparing two
different machines.

**AGAINST THE 513² REFERENCE, WITH ITS OWN CURRENTS, ITS OWN LIMITER POINT AND A
GUESS REBUILT FROM ITS OWN SOURCE:**

| | freegs4e, 513² | MEQ, `k = 3`, 26375 el | apart |
|---|---|---|---|
| `ψ_ax` | 9.308752e-02 | **9.307342e-02** | **1.5e-04** |
| `ψ_bnd` | 2.622462e-02 | **2.622580e-02** | **4.5e-05** |

both inside the **7.92e-04** the 513² run itself sits from its own Richardson
limit, so what is left is the reference's grid error rather than MEQ's. With the
dof-snapped constraint the same two runs read 1.0% and 1.5% out and did not
converge in the mesh at all.

**AND THE REFERENCE IS NOT CONVERGED AT 129² HERE, WHERE THE DIVERTED CASES
ARE** — see the grid scan above, 2.0% in `ψ_ax` and 6.3% in `ψ_bnd` from its own
Richardson limit. That is not a reason to compare against the limit instead: the
limit is a number neither code computed. It is a reason not to read the 2.9e-03
above as a statement about either code's discretisation.

### Regenerating the fixture

Three files in `examples/` come from this directory and none of them can be
rebuilt at test time, because the reference `.npz` they are built from is
gitignored and needs freegs4e, scipy and a Picard solve to recreate.

```sh
# 1. the reference, if it is not already here
python3 fgsref.py --case=H_limited_circular

# 2. the two profile tables.  freegs4e's arrays are dp/dpsi and MEQ's tables
#    are dp/dPsi, so the conversion MULTIPLIES by the span -- see below.
#    convert.py's header is the recipe; the tables carry it in their own.

# 3. the initial guess, over MEQ's half-disc rather than freegs4e's box
python3 mkexactguess.py H_limited_circular.npz \
        ../../examples/limited-tokamak-guess.mesh \
        ../../examples/limited-tokamak-guess.gf 2.60 32

# 4. the mesh -- BUT SEE BELOW BEFORE REGENERATING IT
python3 ../mesh/halfdisc.py --rho 2.6 --size 0.30 --coil-size 0.05 \
        --coil 1.70  0.85 0.10 0.10  --coil 1.70 -0.95 0.10 0.10 \
        --coil 0.50  1.05 0.10 0.10  --coil 0.50 -1.15 0.10 0.10 \
        --plasma 0.60 -0.45 0.85 0.90 --plasma-size 0.09 --check \
        -o ../../examples/limited-tokamak.msh
```

**THE SHIPPED `.msh` IS THE FIXTURE AND THAT COMMAND DOES NOT REPRODUCE IT.**
The four conductors and `--rho` are read straight off the file and are certain;
the three SIZE arguments are recovered from its element edges and are not. The
command above writes 1807 triangles against the shipped 1775, with the same 14
elements per coil and a plasma-region median edge of 0.090 against 0.090 — as
close as guessing gets, and **not close enough**:

| mesh | triangles | `ψ_ax` | from the reference |
|---|---|---|---|
| **shipped** | **1775** | **9.455354e-02** | **2.9e-03** |
| `--plasma-size 0.09` | 1807 | 9.317576e-02 | 1.7e-02 |
| `--plasma-size 0.10` | 1716 | 9.337595e-02 | 1.5e-02 |

**So MEQ's own mesh scatter at this resolution is 1.5–2%, and the shipped mesh
sits at 0.29% by fortune as much as by resolution.** Both regenerated meshes
converge cleanly, both satisfy every border at machine zero, and both are on the
physical branch — they are the same equilibrium at a different discretisation,
not a different root. But `theDriverSolvesALimitedTokamak` gates `ψ_ax` at 1e-02
and either of them would fail it.

**Regenerating the mesh is therefore a re-measurement and not a rebuild.** If it
has to happen, run the fine mesh too (3802 elements reads 9.484390e-02) to
separate scatter from a real change, and move the gate to what the new fixture
measures rather than relaxing it to whatever passes.

**THE GUESS IS PART OF THE PROBLEM STATEMENT, NOT AN OPTIMISATION.** A cold bump
start of `examples/free-boundary-halfdisc.toml`'s kind wanders for 200 iterations
around `‖r‖ = 1.3` and never converges on this problem. `mkexactguess.py` builds
it from the SOURCE rather than by interpolating the answer — `Jtor` on the 2199
core cells plus the four coil filaments, summed by Green's functions — because
freegs4e's `ψ` exists only on its own 1.6 × 1.6 box while MEQ's domain is a
half-disc of radius 2.6. That is also a free check on the whole conversion, since
the sum and the PDE solve share nothing but the source: the reconstruction
reproduces the reference's own `ψ_axis` to **2.5e-04**.

**AND THE ANSWER DOES NOT DEPEND ON THE GUESS'S RESOLUTION**, which is the
property a guess should have and is worth measuring rather than assuming on a
problem with more than one solution. One key changed, everything else fixed:

| `n` | mesh + gf | Newton | `ψ_ax` |
|---|---|---|---|
| 128 | 1.6 MB | 7 | 9.455354e-02 |
| 64 | 395 kB | 8 | 9.455354e-02 |
| 48 | 223 kB | 7 | 9.455354e-02 |
| **32** | **99 kB** | **7** | **9.455354e-02** |

Not one printed digit moves over a sixteenfold change in size, while the work
moves by one iteration. `n = 32` is what ships.

**THE PROFILE TABLES ARE `dp/dΨ` AND freegs4e's ARRAYS ARE `dp/dψ`, AND WITH A
PRESCRIBED CURRENT THE ERROR IS INVISIBLE.** `meq::NormalisedMHDSource::f`
evaluates `F = scale·( μ₀r²p′( Ψ ) + gg′( Ψ ) )/span`, so converting means
MULTIPLYING by `ψ_ax − ψ_bnd`. Getting it wrong is
`examples/rotating-density.dat`'s trap exactly — and with `[source]
PlasmaCurrent` set it produces **no symptom at all**: both profiles carry the
same wrong factor, the scale is an unknown, and the border absorbs it. Measured,
with the tables a factor of `span` too small the scale came back as **6.713e-02
against a span of 6.689e-02** and the equilibrium was right to every digit. **The
tell is a profile scale that is not `O(1)`**, and the regression asserts on it
for that reason.

### Degree 2 is not a cheaper degree 3 here

One key changed, same mesh and same guess:

| | `ψ_ax` | from the reference | profile scale |
|---|---|---|---|
| `k = 2` | 9.634577e-02 | 1.60e-02 | 1.0108 |
| **`k = 3`** | **9.455354e-02** | **2.93e-03** | **0.9945** |

`p`-refinement is worth **5.5×** in `ψ_ax`, 4.1× in `ψ_bnd` and 5.2× in the field
`L2` (3.7e-02 against 7.1e-03). That is the same shape as the fixed-boundary
finding above — *`k = 2, refine = 1` does not converge on MAST or DIII-D while
`k = 3` does on the same mesh* — and it is why the regression runs degree 2 as a
control rather than merely running degree 3.

**`FREE-BOUNDARY-PLAN.md` §7.16 records a WORSE degree-2 failure than this** — 17
Newton steps to `ψ_ax = 2.73e+00`, a single spiking nodal value with the profile
scale raised by 980 to compensate. **That is not reproducible on the code as it
stands**: the recorded configuration was rerun verbatim, same TOML and same 128²
guess, and converges in 8 steps to 9.676040e-02, on the physical branch. The
regression asserts that degree 2 stays on the branch, so it fails if that comes
back.
