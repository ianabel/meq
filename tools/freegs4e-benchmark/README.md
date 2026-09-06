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
6. `compare.py` — reads MEQ's `.nc` and diffs.

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

### AND THE TWO BOUNDARY CONDITIONS DISAGREE, FLAT UNDER REFINEMENT

This is the finding that matters, and it is why the fast one was not simply
adopted. Case A, the same equilibrium computed both ways:

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

The default is **converged at 129²** — `ψ_ax` stable to six figures — while von
Hagenow creeps at about `O(h^1.6)` toward a different value. So **the default is
the reference to use**, `freegs4e` does not agree with itself to better than
3.2e-03 across its own two boundary conditions, and MEQ's agreement with it is
already twenty times inside that spread.

### 2048² is not reachable on this machine, and why

Memory. The von Hagenow ladder runs 293, 503, 1564 MB at 129², 257², 513²;
1025² reached **6.0 GB** with 4 GB free and was stopped rather than risk an OOM
on a shared machine. 2049² needs roughly 24 GB against 15 GB total. With the
*faithful* boundary condition it would additionally be about ninety hours.

**So the honest ceiling here is 513², and the route to round-off is FB-6 rather
than a finer grid.**
