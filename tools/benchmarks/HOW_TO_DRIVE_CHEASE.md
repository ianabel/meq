# How to drive CHEASE

CHEASE is a fixed-boundary inverse equilibrium code, Fortran 90, straight
field-line coordinates with a Hermite bicubic representation. It solves the
**same equation from the same two free functions** MEQ does, which is what makes
its conversion the cheapest of the three: it needs nothing from anybody's
converged answer.

**Tree:** `/home/ian/projects/chease`  ·  **Harness:** `tools/chease-benchmark/`

## The one command

```sh
tools/freegs4e-benchmark/venv/bin/python tools/chease-benchmark/cheaserun.py \
    fixed-h-circular --ns 40 --nt 40
tools/freegs4e-benchmark/venv/bin/python tools/chease-benchmark/converge.py fixed-h-circular
```

`converge.py` is the acceptance and it is a **matrix, not a pair**: a conversion
error is a fixed offset, so only a difference that FALLS as each code is refined
*separately* rules one out. M-158.

## Things that will cost you an afternoon

**IT IS SERIAL AND HAS NOTHING TO OFFER EIGHT THREADS.** No `GOMP_*` symbols,
no MPI calls anywhere in `src-f90/*.f90`, and no external BLAS at all — LAPACK
is linked in. The `CHEASE_F90 = mpif90` lines in `Makefile.define_CHEASEF90`
and `Makefile.define_FLAGS` select a **compiler wrapper** for other machines and
do not make it a parallel code. Do not go looking for a thread knob.

**THERE IS NO ARGV.** `chease_prog.f90:116` and `iodisk.f90:361` open
`EXPEQ` and `chease_namelist` **by hard-coded name** in the working directory,
so every run needs its own directory. `cheaserun.py` makes one per case and
writes both files into it.

**THE AMPLITUDE IS SOLVED FOR, NOT SUPPLIED.** MEQ hides one scalar CHEASE
needs — `psi_ax` — and `convert_chease.py` closes it by a fixed point that is
**exact in one step**, because the map is linear. Do not hand-tune it.

**CHEASE's `psi` IS `−psi_MEQ`**, is zero on `Γ` and negative on the axis, so
the amplitude is `|SIMAG|` and `SIBRY` in its geqdsk is a literal zero.

## What it is good for

**CHEASE AND MEQ AGREE AN ORDER BETTER THAN EITHER AGREES WITH `freegs4e`**, and
that is the finding, not the timing. MEQ `k3r1` against CHEASE `NS80` is
**1.300e-06** in a relative L2 of `psi`, while each against the `freegs4e`
reference reads 1.1443e-04 and 1.1267e-04. **CHEASE's column is flat in `NS` at
1.126–1.135e-04**, which says that number is the REFERENCE's error and not
CHEASE's — and therefore that M-147's column is too. M-158, M-163 section 5.

Its `psi_ax` is converged by `NS = NT = 30`: 6.3545257e-02, 6.3545314e-02,
6.3545318e-02, 6.3545320e-02 at NS = 20, 30, 40, 60. Timings 139.7 / 92.9 /
98.6 / 101.4 s — note the coarsest is the *slowest*, because it takes 3 sweeps
of the amplitude fixed point where the others take 2. **Quote the sweep count
beside the seconds.**
