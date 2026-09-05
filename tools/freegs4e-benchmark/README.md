# MEQ against freegs4e

**An independent check by a different algorithm**, which is the thing
`FREE-BOUNDARY-PLAN.md` FB-6 wants and a fine-mesh self-comparison cannot be.
`../freegs4e` is FreeGS-derived: free boundary by von Hagenow Green's functions,
2nd/4th-order finite differences on a uniform `(R,Z)` grid, Picard with adaptive
blending. MEQ is HDG with Newton. They share the equation and essentially no
code.

**Result, seven configurations** — TestTokamak (three profile shapes), MAST,
MAST-U, TCV, DIII-D:

| case | rel L2 in ψ | axis value |
|---|---|---|
| B TestTokamak, ff′-dominated | **7.1e-04** | +0.06% |
| G MAST-U | **6.4e-04** | +0.04% |
| F DIII-D | 2.9e-03 | +0.46% |
| D TCV | 4.7e-03 | +0.75% |
| E TestTokamak, diamagnetic | 6.4e-03 | +1.09% |
| C MAST | 8.2e-03 | +1.36% |
| **A TestTokamak, classic** | 9.4e-02 | **+15.39%** |

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

## Case A, and why it is left disagreeing

A is on a **different root**, and that is measured rather than supposed: its
+15.39% is flat at +15.39 → +15.41% across refinement levels 2→4 and polynomial
degrees 2→3, converged to six digits. Discretisation error would fall; this does
not.

The reference sits at 4.53e-02, **between** the two roots MEQ finds (1.93e-03
and 5.23e-02), which is the signature of a third, middle root — the unstable
branch of an S-curve, which Newton passes over from either side. Reaching it
needs a start close to the answer rather than a ramp, i.e. MEQ's
`[initialguess] Type = "gridfunction"` restart seeded with the reference. That is
the obvious next step and is not done.

**The source conversion is not the cause and that is checked directly**:
evaluating MEQ's tables on the reference's own ψ reproduces freegs4e's `μ₀ R
J_φ` to **2.3e-05**, so both codes are solving the same equation and disagree
only about which of its solutions to report. `CLAUDE.md` records the same
phenomenon from inside MEQ — three solve routes reaching discrete solutions
9.4% apart on an under-resolved mesh.
