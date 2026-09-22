# How to drive freegs4e

`freegs4e` is a free-boundary tokamak equilibrium solver in Python — von Hagenow
Green's functions and a finite-difference Picard iteration on a uniform `(R, Z)`
grid. It is MEQ's oldest and most-used reference, and the source of every
`ref-n*/` case.

**Tree:** `/home/ian/projects/freegs4e`  ·  **Harness:** `tools/freegs4e-benchmark/`

## The one command

```sh
tools/freegs4e-benchmark/venv/bin/python tools/freegs4e-benchmark/fgsref.py H_limited_circular
tools/freegs4e-benchmark/race.py --grids 129,257,513
```

`fgsref.py` holds the case table (`CASES`); each entry names a machine, a grid,
the profile shape, the `isoflux` control constraints and a limiter. The runner
writes an `.npz` per case carrying `psi`, `psi_axis`, `psi_bndry`, the coil
geometry and currents, the LCFS and the wall.

It runs here only with **relaxed pins and numba 0.67**; the venv already has
them.

## The thing to know before comparing anything against it

**`psi_bndry` IS THE VALUE OF `psi` AT ONE GRID NODE, NOT THE FLUX AT A LIMITER
CONTACT.** `fgsref.py`'s `attach_limiter()` builds a `limiter_ring` — the grid
nodes inside the wall polygon having a neighbour outside it — and
`boundary_and_mask()` takes `psi_bndry = max psi over that ring`. On the limited
case at 513² that is a **632-node** ring and the maximum sits at
`( 0.82500, 0.30000 )`, radius 0.34731 against a wall radius of 0.35000 with
`h = 0.003125`. So:

* the boundary flux carries an **`O(h)` error**, `0.86 h` worth;
* the converged LCFS reaches only 0.34731 and **never touches its own limiter**;
* on the true circle the maximum of `freegs4e`'s own field is 2.588141e-02,
  **1.31% below** the 2.622462e-02 it reports as `psi_bndry`.

**AND THE NODE MOVES 0.6 m ACROSS THE MACHINE AS THE GRID REFINES** — it is at
`( 1.3375, -0.0125 )`, the outboard midplane, on the 129² reference and at
`( 0.8250, -0.3000 )`, the inboard shoulder, on the 513² one. MEQ ships one
config per node (`limited-tokamak.toml` and `limited-tokamak-filament.toml`),
each prescribing that reference's node, so **a MEQ run must take its limiter
point and its coil currents from the same reference it is scored against.**
Mixing them charges the other code for using the point it was correctly given.

**THIS IS ALSO WHY A CLOSE AGREEMENT WITH `freegs4e` ON THIS CASE IS NOT BY
ITSELF EVIDENCE.** MEQ is handed the number it is then scored on. The control
that separates a real agreement from a manufactured one is to **refine the
reference and see whether the agreement follows the shared input** — here it
does. M-164 section 3.

## Other live traps

**ITS DOCSTRING SAYS NORMALISED FLUX AND IS WRONG** in the profile interface —
read the arrays, not the prose.

**THE VON HAGENOW BOUNDARY CONDITION IS INCONSISTENT AS SHIPPED**, and the
default is a cacheable fixed matrix worth about **13×**. See
`freegs4e-boundary-condition-defects` in the project memory and
`CLAUDE_FB.md`.

**INTERPOLATE THE PROFILES WITH A CUBIC SPLINE**, which the `.npz` says in its
own `interpolation_note`: a cubic through `(psi_n, pprime)` reproduces the
solver's own `Jtor` to ~4e-16 where linear gives ~3e-5.

**COIL CURRENTS AND PROFILE AMPLITUDES ARE CONTROL-SYSTEM OUTPUTS**, not inputs
— they converge with the grid (P2 moves 14% and `p` on axis 6.03% between 129²
and 513²). Take them from the reference you are comparing against, never from a
sibling.
