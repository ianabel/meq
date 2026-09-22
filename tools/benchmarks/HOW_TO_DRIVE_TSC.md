# How to drive TSC

TSC is a vintage Fortran tokamak simulation code, finite differences on a
rectangular `(x, z)` grid, driven by **fixed-format punched-card input**. MEQ
uses it as a single free-boundary equilibrium with `NCYCLE = 0` — no time
evolution.

**Tree:** `/home/ian/projects/TSC`  ·  **Harness:** `meq-race/limited/`

## The one command

```sh
cd /home/ian/projects/TSC/meq-race/limited
./run.sh race_c --nx  61 --nz  69 --ppres 50235 --isurf 1 --neqmax 30000 --npsi 400 --tfmult 2
./run.sh race_m --nx 121 --nz 137 --ppres 48466 --isurf 1 --neqmax 30000 --npsi 400 --tfmult 2
./run.sh race_h --nx 181 --nz 205 --ppres 47905 --isurf 1 --neqmax 30000 --npsi 400 --tfmult 2
```

`mkdeck.py` writes the deck, `run.sh <tag> [options]` makes a clean run
directory and runs the binary in it. `shoot.py probe <tag>` pulls
`iter, gp1, gp2, psimin, psilim, span, xmag, zmag` out of the last cycle.

**The diagnostics are free**: `--isurf 1 --ipest 1 --npsi 400` adds the flux
surfaces and a geqdsk and reads **17.28 s** against the same grid's 17.28 s
without them. Turn them on.

## Things that will cost you an afternoon

**IT IS SERIAL AND CANNOT BE MADE OTHERWISE.** No `GOMP_*` symbols, no MPI, no
`-fopenmp` anywhere in `tsc.gfortran.mk`. It links `libmkl_rt` so MKL *could*
thread inside LAPACK, and does not: **17.49 s at `MKL_NUM_THREADS=1` against
17.89 s at 8, both at 99% CPU.** Do not give it threads and do not report it as
having declined them.

**`PPRES` IS AN INPUT WHOSE CORRECT VALUE IS DEFINED BY AN OUTPUT, AND IT IS
RESOLUTION-DEPENDENT.** It sets the `p'`/`gg'` amplitude ratio, the map between
the two depends on the discrete equilibrium, and it must be **re-shot at every
grid**: 50235 / 48466 / 47905 for the three above. Two secant steps each;
`shoot.py next` does the arithmetic. **Quote the ratio error beside any answer.**

The usable window is narrow and bounded above by an equilibrium limit: below
about `PPRES = 30000` (coarse) the set `{psi < psilim}` is not closed on the
inboard side and the current fills the card-03 `XLIM` box; above about 52000
(coarse) / 48500 (fine) the solve runs away or needs thousands of Picard
iterations. **`NEQMAX = 30000` is needed at 121 × 137 and finer** — at the
default 2000 the medium grid stops at `ineg = 12` with the answer already within
1.5%.

**THERE IS A SIZING BUG IN THE SHIPPED TREE AND IT FIRES ON ANY MACHINE WITH NO
PASSIVE STRUCTURE.** `modules/alloc_scr9.f90` allocates `ans1v / sns1v / ans2v /
sns2v` with third dimension `pnwire` (internal wires) while `bounda.f90` indexes
them by **coil number** up to `ncoilte`. MEQ's machine has four external coils
and no conducting structure, so `nwire = 0`, `pnwire = 2`, `ncoilte = 4`, and
coils 3 and 4 are written and read past the end of the array — returning almost
coil 1's field. **Fix: allocate with `pncoil`**, which is always `>= pnwire`.
It is latent in every deck TSC ships (ITER has 361 wires, C-Mod 569).

Diagnosed with `--lrswtch 5`, which names a coil group no coil belongs to and so
solves the **coils only, no plasma**. The four-coil vacuum field against an
exact elliptic-integral filament sum then reads 1.1589–1.2086 before the fix and
**1.0000–0.9999 after**. Moving coil 3 anywhere changed TSC's answer not at all
while its own echo showed the moved position — that is an aliased array, not a
physics choice.

**THE PROFILE RESIDUAL THROUGH THE G-FILE HAS A FLOOR OF ABOUT 1e-02 AND IT IS
NOT THE PROFILES.** `wrpest.f90` writes `pprime`/`ffprim` off the **flux-surface
spline**, not the analytic form the solver integrated; it also overwrites the
first two samples by linear extrapolation and forces the last to zero. Raising
`npsit` shrinks the residual and moves nothing else. The profiles the solver
actually uses are `Psihat²` exactly.

## What TSC is NOT doing that MEQ is

**TSC HAS NO EQUIVALENT OF MEQ'S FREEZE-AND-RE-DECIDE PLASMA SUPPORT.** Its
plasma is the pointwise set `{psi < psilim}`, bounded by card 03's
`XLIM`/`ZLIM`/`XLIM2` box.

**ITS LIMITER IS THE PRESCRIBED POINT ON CARD 05** — for MEQ's posing,
`( 1.3375, 0 )`, which is `examples/limited-tokamak.toml`'s point and therefore
the **129²** reference's staircase node. **Score TSC against the 129² reference,
not the 513² one**: against its own it reads +0.92% in `psi_ax` and +2.32% in
`psi_bnd` at 181 × 205, monotone in the grid; against 513² it appears to be 8.5%
wrong in `psi_bnd` and is being charged for using the point it was given.
M-164 section 3.

TSC's `psi` is `−psi_MEQ`; compare magnitudes.
