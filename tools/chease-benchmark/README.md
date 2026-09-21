# MEQ against CHEASE, on the six shipped fixed-boundary machine cases

> **NO TIMING IN THIS FILE IS PUBLISHABLE.** Every second here was taken on a
> machine carrying another agent's MFEM build and the user's own solves, at a
> one-minute load average of 5 to 20 against 8 physical cores. They are here to
> say what a cell costs and nothing else. `CLAUDE.md` records a 490 s / 540 s
> spread on identical code on an idle machine; this is worse than that. **The
> ERROR columns do not depend on load and are reproducible as they stand.**

`examples/fixed-*.toml` are six real machines posed as fixed-boundary
Grad–Shafranov problems on a smooth MXH curve —
[M-147](../../MEASUREMENTS.md#m-147) is what MEQ reports on them. This directory
races them against **CHEASE**, the EPFL/SPC fixed-boundary equilibrium code
(Lütjens, Bondeson & Sauter, *Comput. Phys. Commun.* **97** (1996) 219), at
`/home/ian/projects/chease`.

**CHEASE IS THE CLOSEST COMPARISON MEQ HAS, AND THAT IS THE POINT OF ADDING
IT.** `freegs4e` solves the same equation by a different algorithm but
*free* boundary, so the shipped cases are a re-posing of its answer. DESC does
not even carry the same unknown. CHEASE solves **the same equation, fixed
boundary, from the same two free functions, against the same normalised flux
label, on the same analytic curve** — it differs from MEQ only in the
discretisation: bicubic Hermite finite elements on a flux-aligned mesh that
CHEASE constructs for itself, against MEQ's hybridizable discontinuous Galerkin
on an unfitted triangulation with a transferred boundary datum.

**So there is almost no conversion, and what conversion there is has no free
parameters.** The DESC arm needs a flux-label map that has to come out of
somebody's converged equilibrium — `tools/desc-benchmark/README.md` spends a
section on how to stop that being circular. **This arm reads nothing but the
shipped `.toml`, `-pprime.dat`, `-ggprime.dat` and `-meta.json`.** Even
`psi_ax`, which MEQ makes an unknown of a bordered Newton and CHEASE's
formulation requires as an input, comes out of a two-run fixed point that is
exact in one step. See *The amplitude closure*.

---

## What was built, and the exact commands

`/home/ian/projects/chease` is a clean `git` checkout of
`https://gitlab.epfl.ch/spc/chease`, on `master` at `fb46366`, **with no
uncommitted changes before or after this work** — the build products are all
covered by its own `.gitignore` (`*.[oa]`, `*.mod`, `src-f90/chease`).

```sh
cd /home/ian/projects/chease/src-f90
make CHEASE_MACHINE=linux_nohdf5 CHEASE_F90=gfortran chease -j4
```

**It builds first time, zero errors and 98 warnings**, against **gfortran
15.2.0** — which is worth saying because CHEASE is old
Fortran and gfortran has been rejecting legacy argument mismatches as errors
since version 10. None of that is needed here: no `-fallow-argument-mismatch`,
no `-std=legacy`, no patches. The binary is `src-f90/chease`, 5.7 MB,
dynamically linked, and needs no external library at all — **no HDF5, no
NetCDF, no MKL, no LAPACK**: `linux_nohdf5` links `dpgbtrf_s.o`, CHEASE's own
banded factorisation, and `interpos_source.o`, its own spline package.

**Both `make` variables have to be given on the command line**, because
`Makefile.define_MACHINE` and `Makefile.define_CHEASEF90` select from
`$(HOSTNAME)` against a list of EPFL, ITER and CINECA machines, and on anything
else fall through to `CHEASE_MACHINE = none` with `CHEASE_F90 = ifort`. There is
no ifort here. `make` reads them from the environment too (`Makefile:64-69`
guards both `include`s with `ifndef`).

The flags that come out of that pair (`Makefile.define_FLAGS:235-245`):

```
gfortran -g -O2 -Wall -ffree-line-length-none -fdefault-real-8 -fdefault-double-8
```

**`-fdefault-real-8` is load bearing** — CHEASE declares its working precision
as plain `REAL` and relies on the promotion, so a build without it would
compile and give single-precision answers.

Not built, and not needed: `chease_hdf5` (wants `futils`), `chease_itm` /
`chease_imas` (want the EU-ITM or ITER access layer), `libchease_kepler`.
`make test_chease` needs `matlab -nodesktop`, which is not on this machine, so
CHEASE's own test target was **not** run; the validation here is the round trip
in *How each convention was established* instead.

The harness runs on **`tools/freegs4e-benchmark/venv`** — no new virtualenv. It
needs `numpy`, `scipy`, `netCDF4`, `skimage` and `freegs4e` (for
`_geqdsk.read`), and that venv is the only one on this machine with all five.

---

## What CHEASE solves, and what it is handed

**The problem.** Fixed-boundary axisymmetric Grad–Shafranov, in the form

```
Delta* psi = -mu0 r^2 p'( Flux ) - T T'( Flux ),      Flux = ( psi_edge - psi )/( psi_edge - psi_axis )
```

with `psi = 0` on a prescribed closed curve, `T = g = R B_phi`. `Flux` is 1 on
the magnetic axis and 0 on the boundary — **which is MEQ's `Psi` exactly**, and
is the one relabelling this comparison does *not* have to do. (`freegs4e`'s
`psi_n` runs the other way, and `examples/fixed-*-pprime.dat`'s own header is
about that.)

**The input.** Two files in the working directory, both by hard-coded name —
there is no argv.

* **`chease_namelist`**, a Fortran namelist `&EQDATA … /`
  (`src-f90/COMDAT.inc:95-120`), read on unit 13 at
  `src-f90/chease_prog.f90:144`. **Its first four lines are not comments**:
  `chease_prog.f90:133-141` reads them into `LABEL1..LABEL4` and only then
  rewinds, so a file that opens with `&EQDATA` loses its first four namelist
  lines. Defaults are set by `PRESET` (`src-f90/preset.f90`) before any input is
  read.
* **`EXPEQ`**, free-format ASCII, opened at `src-f90/iodisk.f90:361`, parsed at
  `:369-518`. Every read is list-directed, so the layout below is
  record-by-record and not column-by-column:

  ```
  ASPCT              inverse aspect ratio ( Rmax - Rmin )/( Rmax + Rmin )
  RZ0C               Z of the geometric centre, / R0EXP
  PREDGE             p on the boundary, CHEASE units
  NBPS               number of boundary points  ( <= NPBPS = 2300 )
  R  Z               NBPS pairs, both / R0EXP
  NPPF1  NPPFUN      profile table length, and 4 = "p' is given" ( 8 = "p" )
  NSTTP              1 = "the second table is T dT/dpsi"  ( 2 I*, 3 I||, 4 j||, 5 q )
  s        x NPPF1   the abscissa,  s = sqrt( psi_N ),  0 on axis, 1 on the edge
  p'       x NPPF1
  T dT/dpsi x NPPF1
  ```

  There is no `EXPEQ.IN`: CHEASE reads only `EXPEQ`. `EXPEQ_EXPEQ.IN` is an
  *echo* it writes.

**The output.** `EQDSK_COCOS_02.OUT` and `EQDSK_COCOS_02_POS.OUT`, standard
G-EQDSK, `NRBOX x NZBOX` — written **unconditionally** whenever `NVERBOSE >= 1`,
because the guard that used to depend on `NEQDSK` is commented out
(`src-f90/stepon.f90:255-259`). They are the only files carrying `psi( R, Z )`;
`EXPEQ.OUT`, `chease.dat` and `NUPLO` carry profiles and surfaces but no box.
**The EQDSK has a non-standard tail** — `iodisk.f90:2220` appends the same
human-readable summary table that goes on `EXPEQ.OUT` — so a reader that stops
after the limiter block is fine and one that keeps going is not.
`freegs4e._geqdsk.read` stops, which is why it is the reader used here.

---

## The conventions, and how each one was established

**Every row was checked twice: once by reading the source, and once by making
CHEASE state the number itself.** `CLAUDE.md` is explicit that a convention read
in a manual and not checked against the code is not evidence, and this tree has
been bitten by `freegs4e`'s own docstring being wrong about exactly this kind of
thing. The *measured* column below is the stronger half.

| convention | the code | measured |
|---|---|---|
| **`psi` is Wb/rad, zero on the boundary, NEGATIVE on the axis, increasing outward** | `norept.f90:245-258` sets the gauge; `iodisk.f90:2110` writes a literal zero for `SIBRY`; `psibox.f90:352` subtracts `CPSRF` | the shipped EQDSK reads `SIMAG = -4.222889774e-02`, `SIBRY = 0.0`, `CPASMA = +208580 A` |
| **the profile abscissa is `s = sqrt( psi_N )`, 0 on axis** | `ppspln.f90:76-79` forms `sqrt( 1 - psi/SPSIM )` before interpolating; `premap.f90:62` builds the iso-mesh from its square | interpolating CHEASE's own `EXPEQ.OUT` table against its own EQDSK `pprime` reads **3.3e-04** under `s = sqrt( psi_N )` and **2.5e-01** under `s = psi_N` — and the endpoints agree to 9 figures |
| `psi_CHEASE = psi_SI /( R0EXP^2 B0EXP )` | `iodisk.f90:1103`, `:2105`, `:2175` | printed: `-3.05238308e-02` against `-4.22288977e-02` SI, with `R0EXP = 0.99408`, `B0EXP = 1.4` |
| `p_CHEASE = p_SI mu0/B0EXP^2` | `iodisk.f90:1109`, `:2151` | `PREDGE = 2.3358e-10` ⇒ `3.6431903e-04 Pa`, against the EQDSK's `PRES[-1] = 3.643190338e-04` |
| `p'_CHEASE = p'_SI mu0 R0EXP^2/B0EXP` | `iodisk.f90:1108`, `:2169` | `EXPEQ.OUT`'s own axis value **−0.27972256** against **−0.2797225613** rebuilt from the EQDSK's SI `pprime[0]` |
| `TT'_CHEASE = TT'_SI/B0EXP` | `iodisk.f90:1107`, `:2160` | **−1.9426126** against **−1.9426125671**, same route |
| `I_CHEASE = I_SI mu0/( R0EXP B0EXP )` | `iodisk.f90:1105`, `:2119` | printed `1.88336373e-01` against `2.08580503e+05 A` |
| **`T( edge ) = R0EXP * B0EXP` identically, with `NTMF0 = 0`** | `isofun.f90:327-331` sets `TMF( KN ) = 0.5` and takes `sqrt( 2 x )` | the EQDSK's `FPOL[-1] = 1.391712` = `0.99408 x 1.4` exactly |
| boundary `R, Z` are metres **divided by `R0EXP`**, not by the minor radius | `iodisk.f90:1116-1120` does that division explicitly on the EQDSK path | `0.71314032 x 0.99408 = 0.70891853` against the run's printed `RMIN [m]` |
| **COCOS 2 in and out by default**; `EQDSK_COCOS_02` is the COCOS index, not a version | `globals.f90:94-95`; the filename is `I2.2` of `COCOS_OUT`, `iodisk.f90:2045` | with `SIGNIPXP = SIGNB0XP = +1` every sign factor is 1 and the `_POS` file's body is byte-identical |

### The one whole-code check: a round trip

The strongest single piece of evidence is not in that table. CHEASE writes
`EXPEQ.OUT` in the same format it reads, so **its own output is a valid input**.
Taking a run under `NCSCAL = 1` (rescale to `QSPEC = 0.9`) with `p` given
(`NPPFUN = 8`) and `I*` (`NSTTP = 2`), and feeding its `EXPEQ.OUT` straight back
under `NCSCAL = 4` with `p'` (`NPPFUN = 4`) and `TT'` (`NSTTP = 1`):

| | first run | round trip | relative |
|---|---|---|---|
| `psi_axis` | −4.22288977e-02 | −4.22289224e-02 | 5.8e-07 |
| `Ip` | 2.08580503e+05 | 2.08580721e+05 | 1.0e-06 |
| `q0` | 0.900000000 (imposed) | 0.900008137 (**not** imposed) | 9.0e-06 |
| `p` on axis | 6.61542092e+03 | 6.61544314e+03 | 3.4e-06 |

The residual is the 141-point table and the 8-figure text format, not a method
difference. **That single experiment establishes the sign, the normalisation,
the flux label and the format of all three input blocks at once**, and it
establishes that `NCSCAL = 4` really does take the profiles as given — `q0`
comes out at 0.900008 having been told nothing about `QSPEC`.

---

## What CHEASE would otherwise do to the data, which is the dangerous part

**Four defaults rewrite the input, and none of them says so at the default
verbosity.** Each is set the other way in `convert_chease.namelist_text`.

| key | default | what it does |
|---|---|---|
| **`NCSCAL`** | **2** (`preset.f90:91`) | rescales **both** profiles by `CURRT/CUROLD` so the total current hits `CURRT` (`norept.f90:123-157`). 1 and 3 do the same to hit `QSPEC`. **`NCSCAL = 4` is the no-rescale value** — `norept.f90`'s chain ends at `:224` with no branch for it, and `chease_input_choices_default.xml:48-52` names it `rescaling_options/no_rescaling` |
| **`TENSPROF`** | **−0.3** (`preset.f90:232`) | runs an `interpos` smoothing spline over **both** input profile tables (`iodisk.f90:490-495, 530-535`), with `dp'/ds = 0` forced on the axis |
| **`TENSBND`** | **−0.3** (`preset.f90:233`) | runs a periodic smoothing fit over the boundary **and doubles the point count** (`iodisk.f90:409-414`) |
| `NBLOPT`, `NBSOPT` | 0 | off by default, but if turned on they *rewrite `p'`* toward marginal ballooning stability or a bootstrap-consistent profile (`ballit.f90:179,199`) |

`ASPCT` and `RZ0C` are also recomputed from the boundary and silently replaced
if they disagree by more than 20% / `0.05 ASPCT` (`iodisk.f90:394-406`), so they
are diagnostics rather than inputs. They are written correctly anyway, so the
override never fires.

**And one combination is a silent multiply-by-zero.** `NCSCAL = 4` with
`NSTTP = 3` or `4`: `guess.f90:424-445` branches on `NCSCAL` 1, 3 and 2 and has
no arm for 4, so `SCALE` keeps the module initialiser `globals.f90:632
scale = 0` and the current-profile coefficients are multiplied by it. This
conversion uses `NSTTP = 1` and is not exposed, but the combination is
reachable from a namelist and nothing warns.

**The settings that make this MEQ's problem rather than a CHEASE one**, all in
`convert_chease.namelist_text` with the reasons beside them:

```
NSURF = 6, NEQDSK = 0     boundary point by point from EXPEQ
NPPFUN = 4                the first table is p', not p
NFUNC = 4, NSTTP = 1      the second is tabulated, and it is T dT/dpsi
NCSCAL = 4                do not rescale
NRSCAL = 0                do not move lengths to put the MAGNETIC axis at 1
NTMF0 = 0                 anchor T at the edge -- this is what makes B0EXP mean g_at_gamma/R0EXP
NBLOPT = 0, NBSOPT = 0    do not rewrite p'
CPRESS = 1, CFNRESS = 1   no one-shot multipliers
TENSPROF = 0, TENSBND = 0 do not smooth the inputs
```

With `NSTTP = 1` and `NFUNRHO = 0` CHEASE also skips its outer current-profile
iteration entirely (`itipr.f90:35`), so the map from the two tables to the
answer is a single Picard solve of the nonlinearity and nothing else.

---

## The conversion

### The two scalars Grad–Shafranov does not fix

`p_at_gamma` and `g_at_gamma`, **taken from `<stem>-meta.json`** (present as of
`b375fdc`) and not from the `freegs4e` reference. Neither touches `psi`:

* `g_at_gamma` enters as the **units choice**. With `NTMF0 = 0` CHEASE anchors
  `T = 1` at the edge, so `T_SI( Gamma ) = R0EXP B0EXP` identically, and setting
  `B0EXP = g_at_gamma/R0EXP` is what makes CHEASE's `F` and `q` physical. The
  Grad–Shafranov operator sees `g` only through `TT'`.
* `p_at_gamma` enters as `PREDGE`, record 3, and only sets the constant of
  integration on the reported pressure.

`R0EXP` is the boundary's **geometric** centre, `( Rmax + Rmin )/2`, and not
MXH's `R0` parameter — `NRSCAL = 0` normalises the boundary so that its
geometric major radius is 1 in CHEASE units, so a different `R0EXP` would put
the plasma somewhere else in metres. On `fixed-h-circular` the two agree to
1.1e-15 because the case is nearly circular; on a shaped one they will not.

### The sign, which is the one that converges to the wrong answer

CHEASE's `psi` is the negative of MEQ's. So

```
p'_CHEASE  = -( dp/dPsi )/psi_ax        TT'_CHEASE = -( g dg/dPsi )/psi_ax
psi_MEQ    = -psi_CHEASE
```

and the tables handed to CHEASE are **negative** where MEQ's are positive.
`convert_chease.py --check` prints both ends of both tables against the meta
file's own two scalars, because this is precisely the failure `CLAUDE.md` names.

**A trap inside that check, met and recorded**: `meta["pprime_at_gamma"]` is
`dp/dpsi` in **Pa per Wb/rad** — the `.toml` says so in its units string — while
the `.dat` table's value column is `dp/dPsi`, **per unit normalised flux**, which
the table's own header describes as carrying "ONE factor of `dpsi/dPsi =
psi_ax`". They differ by `psi_ax`, a factor of 15.7 on this case, and comparing
them without it reads a *correct* conversion as badly wrong. With the factor in,
the conversion reproduces the table to every digit:
`3.165642787599e+02` both ways.

### The amplitude closure, which is the only real arithmetic

MEQ and CHEASE are handed **different functions**, and the difference is exactly
`psi_ax`:

```
MEQ      Delta* psi = -[ mu0 r^2 P( Psi ) + G( Psi ) ]/psi_ax,  psi_ax = max psi, an UNKNOWN
CHEASE   Delta* psi = -mu0 r^2 p'( Flux ) - TT'( Flux )
```

Matching them needs `p' = -P/A` with `A` the axis flux **of the answer**, so a
trial has to be supplied and then made self-consistent. Using MEQ's `PsiAxis` or
the reference's `psi_ax` would work and would be circular.

**It is not necessary, because the map is exactly reciprocal.** Scaling both
profiles by `lambda` scales `psi` by `lambda` and leaves `Flux` invariant, so
with a trial `A0` CHEASE returns `A_C( A0 ) = K/A0`. Self-consistency is
`A0 = A_C( A0 )`, i.e.

```
A = sqrt( A0 * A_C( A0 ) )      for ANY A0
```

Measured on `fixed-h-circular` at `NS = NT = 40`, starting sixteen times away
from the answer:

| sweep | trial `A0` | returned `A_C` | residual |
|---|---|---|---|
| 0 | 1.0000000000e+00 | 4.0380075770e-03 | 2.466e+02 |
| 1 | 6.3545319080e-02 | 6.3545317030e-02 | **3.226e-08** |
| 2 | 6.3545318055e-02 | 6.3545318050e-02 | 7.998e-11 |

**One step, from a trial two and a half orders of magnitude out.** That is the
reciprocal law being exactly right rather than approximately, and it is itself a
check on the conversion: a wrong scaling anywhere upstream would leave a
fixed point that still exists but does not close in one sweep.

Contrast `tools/desc-benchmark`'s `--self-consistent`, which iterates a map that
is *not* exact, cycles with period 2 undamped, and needs eight sweeps at
`mix = 0.5`. The difference is a property of the two conversions, not of the two
codes.

**The second run is not optional even though `psi` could be rescaled by hand**,
because `T` does not obey the same law: `T^2 = T^2( edge ) + 2 int TT' dpsi`
picks up `lambda^2` where `T^2( edge )` is fixed, so `q`, `F` and `beta` from a
rescaled first run would be wrong while `psi` was merely scaled.

**So `psi_ax` is CHEASE's own independent measurement**, and comparing it with
the reference's is a result rather than an input.

### Four numbers the conversion never saw, and what CHEASE says about them

At `NS = NT = 40`, `NPSI = 100`:

| | CHEASE, from the tables alone | independently | relative |
|---|---|---|---|
| `psi_ax` | 6.3545318050e-02 Wb/rad | 6.354433603e-02, the reference's | **1.5e-05** |
| `Ip` through `Gamma` | 2.998698e+05 A | 2.9985773e+05 A, `convert_desc.py --refine 4`'s contour integral of the reference | **4.0e-06** |
| `g` on axis | 1.05435512 T m | 1.0543556, the reference's `fpol[0]` | **4.6e-07** |
| magnetic axis | ( 1.055133, 0.000000 ) m | ( 1.055122, 0.000000 ), the reference's | 1.1e-05 m |

The `Ip` row is the sharpest: a Fortran finite-element code integrating its own
solution, and a Python contour integral over a 257² finite-difference field,
agreeing to four parts in a million on a number neither was told.

### What is NOT exact, in one list

| | size on `fixed-h-circular` | who pays it |
|---|---|---|
| `Gamma` sampled at 2048 points instead of analytically | **1.5e-06** in enclosed area | CHEASE only |
| `Gamma` then re-represented on CHEASE's own `NT` poloidal mesh | not separately measured; folded into the `NS`/`NT` sweep | CHEASE only |
| the profile tables resampled onto 257 points even in `s` | not separately measured; the round trip bounds the same operation at ~3e-06 | CHEASE only |
| CHEASE's EQDSK box bilinearly interpolated onto MEQ's grid | measured by the `NRBOX x2` control below | the comparison |
| `MEQ`'s `Gamma` | exact — it **is** the MXH curve | — |
| the `freegs4e` reference as a yardstick | M-147 puts MEQ 5.0e-05 to 1.1e-04 from it | both, and only when it is an arm |

**Three of those six rows are CHEASE's alone and all three are small**, against
five of seven and much larger for the DESC arm. That difference is the whole
argument for having this one.

---

## Does the difference fall under refinement?

**THIS IS THE ACCEPTANCE.** `CLAUDE.md`'s *Testing stance* records why a
single-resolution agreement is worth nothing: *a wrong sign convention
converges, at the right rate, to the wrong function*. A conversion error is a
fixed offset and does not care about either mesh, so two codes whose difference
**falls** as each is refined cannot be agreeing for that reason.

`converge.py` runs a matrix rather than a pair: MEQ's ladder
(`PolynomialDegree` x `RefinementLevels`, taken through
`race_desc.meq_config` so it is the *same* ladder the DESC race uses) down one
axis, CHEASE's `NS x NT` along the other.

### The matrix

`fixed-h-circular`, MEQ's output grid 257², CHEASE's `NPSI = NCHI = NISO = 200`
and `NRBOX = NZBOX = 999`. Relative L2 of `psi` over MEQ's own nodes, band
dropped, through `compare.py`'s norm.

| | `NS=NT=20` | `NS=NT=40` | `NS=NT=80` |
|---|---|---|---|
| **MEQ `k=1 r=0`** | 8.925e-04 | 8.912e-04 | 8.911e-04 |
| **MEQ `k=2 r=0`** | 1.933e-04 | 1.945e-04 | 1.942e-04 |
| **MEQ `k=3 r=0`** | 1.446e-05 | 1.123e-05 | **2.696e-06** |
| **MEQ `k=1 r=1`** | 2.875e-04 | 2.886e-04 | 2.885e-04 |
| **MEQ `k=2 r=1`** | 1.999e-05 | 1.722e-05 | 1.340e-05 |
| **MEQ `k=3 r=1`** | 1.390e-05 | 1.085e-05 | **1.300e-06** |

**THE DIFFERENCE FALLS BY 686x FROM ONE CORNER TO THE OTHER, AND THE SHAPE OF
THE FALL IS THE RESULT RATHER THAN THE FACTOR.**

* **The first three rows are FLAT in `NS`** — a sixteen-fold change in CHEASE's
  element count moves them by under 1%. MEQ's own discretisation error dominates
  there, so refining the other code cannot be seen. That is what a row is
  supposed to do.
* **The last three rows FALL in `NS`**, by 5.4x, 1.5x and 10.7x. MEQ has passed
  CHEASE, so CHEASE's error becomes the visible one. The crossover is between
  `k = 2` and `k = 3`, which is a statement about where these two
  discretisations sit against each other at these sizes and about nothing else.
* **Every column falls with MEQ's ladder** — at `NS = 80`, 8.911e-04 →
  1.942e-04 → 2.696e-06 across `k = 1, 2, 3` at `r = 0`, a factor of 330.

**A conversion error is a fixed offset and cannot do any of that.** This is the
check `converge.py` exists for and it passes.

### The deepest cell is the INSTRUMENT, and finding that out changed the table

**The single largest item at the bottom of the matrix is not either code.** It
is the bilinear interpolation of CHEASE's EQDSK box onto MEQ's grid inside
`compare.py`. Measured by holding everything fixed — MEQ at `k=3 r=1`, CHEASE at
`NS=NT=80` — and changing only the grid CHEASE prints its answer on:

| `NRBOX = NZBOX` | 129 | 257 | 513 | 999 |
|---|---|---|---|---|
| MEQ vs CHEASE | 7.6057e-05 | 1.8991e-05 | 4.7723e-06 | **1.2999e-06** |
| ratio to the last | — | 4.005 | 3.980 | 3.671 |
| CHEASE's `psi_ax` | 6.3545319800e-02 | 6.3545319800e-02 | 6.3545319800e-02 | 6.3545319800e-02 |

**Exactly second order, and `psi_ax` identical to eleven digits throughout** — so
the solve did not move at all and the whole of that column is how the answer was
printed. Extrapolating the `h^2` law says the box error at 999 is about 1.26e-06
against a measured cell of 1.2999e-06, i.e. **the true MEQ-against-CHEASE
disagreement at `k=3 r=1` is below 1.3e-06 and this comparison cannot resolve
it.**

The whole matrix above was first taken at `NRBOX = 513` and read 4.77e-06 in
that cell. It was re-taken at 999 rather than quoted, because the control said
the number belonged to the grid. **This is the same lesson `tools/desc-benchmark`
records about contour chording**: the biggest item in a two-code comparison is
routinely the instrument, and it is found by changing the instrument and nothing
else.

The other control is clean: **`NPSI = NCHI = NISO` doubled from 200 to 400 moves
that cell from 4.7723e-06 to 4.7723e-06** — not a rounding of two nearby
numbers, the same five figures. The flux-surface mapping mesh is not what is
being measured.

### `psi_ax`, the scalar both codes solve FOR

Sharper than any norm, because it is one number with no interpolation anywhere
near it. CHEASE's own self-convergence first:

| | `NS=20` → `40` | `NS=40` → `80` |
|---|---|---|
| CHEASE `psi_ax` moves by | 9.57e-07 | **2.74e-08** |

so `6.354531980e-02` at `NS = 80` is converged to a few parts in `10^8`. Against
it:

| MEQ rung | `psi_ax` | vs CHEASE `NS=80` | vs the `freegs4e` reference |
|---|---|---|---|
| `k=1 r=0` | 6.359424741504e-02 | 7.700e-04 | 7.855e-04 |
| `k=1 r=1` | 6.353512042584e-02 | 1.605e-04 | 1.450e-04 |
| `k=2 r=0` | 6.353706768217e-02 | 1.299e-04 | 1.144e-04 |
| `k=2 r=1` | 6.354578707868e-02 | 7.353e-06 | 2.284e-05 |
| `k=3 r=0` | 6.354539589856e-02 | 1.198e-06 | 1.668e-05 |
| `k=3 r=1` | 6.354532156975e-02 | **2.785e-08** | 1.551e-05 |

**MONOTONE OVER FOUR AND A HALF ORDERS OF MAGNITUDE, ENDING AT 2.8e-08.** Two
codes with nothing in common but the equation and the input files agree on the
axis flux to eight significant figures.

**AND THE LAST COLUMN STOPS FALLING AT 1.55e-05, WHICH IS THE REFERENCE'S OWN
ERROR AND NOT MEQ'S.** The `freegs4e` reference sits 1.548e-05 from CHEASE's
converged `psi_ax`; MEQ's distance from the reference flattens at 2.28e-05,
1.67e-05, 1.55e-05 over the last three rungs while its distance from CHEASE
keeps falling by two more orders. **MEQ is converging; the yardstick is not.**

### The reference's error, measured rather than bounded

The same statement from the field norm rather than from one scalar:

| | `NS=20` | `NS=40` | `NS=80` |
|---|---|---|---|
| CHEASE vs `freegs4e` | 1.1237e-04 | 1.1319e-04 | 1.1267e-04 |

**Flat to 0.7% while CHEASE's own element count changes sixteen-fold**, and flat
while CHEASE-against-MEQ over the same sweep falls by a factor of ten. A number
that does not move when you refine the code producing it is a property of the
thing it is measured against. **So ~1.13e-04 is what the 257² `freegs4e`
free-boundary solve is worth on this case**, and M-147's 9.7e-05 for MEQ at
`k=2 r=0` is a difference between two comparable errors rather than a
measurement of MEQ's.

**That also makes one of M-147's rows legible in a way it was not.** MEQ at
`k=2 r=0` reads **9.7266e-05** against the reference — better than CHEASE's
1.1267e-04 — and **1.942e-04** against CHEASE. Both are true and the first is
partly cancellation: MEQ's error and the reference's happen to lie on the same
side. **Two codes agreeing better than either is accurate is exactly the failure
`compare_desc.py`'s header warns about, caught here in the wild.** The MEQ arm
here reproduces M-147's published 9.712e-05 as 9.7266e-05, so this is the same
measurement and not a different one.

### The floor

**NOT REACHED.** The diagonal is still falling in both directions at the bottom
right, and what stops it is the EQDSK box rather than the MXH curve, the profile
tables or either solver. So the floor of this comparison is **below 1.3e-06 and
unmeasured** — which is two orders of magnitude below the floor of the
`freegs4e` comparison and 90x below the DESC one's 1.18e-04.

Removing the instrument is cheap and was not done: CHEASE's `NRBOX` cannot go
past 999 (see below), so the route is either to put the EQDSK box exactly on
MEQ's output grid — `cheaserun.py --match-box` is written and **untested** — or
to compare on CHEASE's grid instead of MEQ's, which `compare.py` is not built
for.

### A CHEASE defect found on the way

**`NRBOX` or `NZBOX` at 1000 or more writes a G-EQDSK no reader can parse, and
CHEASE exits 0.** The header format is `(A28,I2.2,A10,A8,3I4)`
(`src-f90/iodisk.f90:2092`, and `9380 FORMAT(A40,A8,3I4)` at `:2667` for the
`_POS` file), so a four-digit box size fills its `I4` field and the three
integers run together:

```
      FROM CHEASE BUT COCOS=02 ,SI UNITS20260920   310251025
```

`freegs4e._geqdsk.read` dies on `int( 'UNITS20260920' )`.
`convert_chease.namelist_text` refuses above 999 rather than letting it be met
downstream. **Not reported upstream** — `CLAUDE.md`'s standing rule is to write
MFEM requests into `../mfem-hdg-dev/doc/`, and there is no equivalent channel
here; it is recorded for whoever wants to send it.

### Cost, ROUGH and on a loaded machine

`fixed-h-circular`, one solve including every sweep of the amplitude closure, at
`NPSI = 200` and `NRBOX = 999`: **CHEASE `NS=20` 99.9 s (2 sweeps), `NS=40`
52.1 s (1 sweep), `NS=80` 112.9 s (2 sweeps)**, against **MEQ 5.7 s at `k=1 r=0`
and 5.8 s at `k=2 r=0`**. Do not read anything into these. The load average was
5 to 20 against 8 physical cores, a CHEASE sweep at `NRBOX = 999` spends a large
and unmeasured share of itself writing a 999² ASCII box twice, and the two arms
are not matched in resolution by any definition. **A real timing needs an idle
machine and a stated matching**, exactly as `tools/desc-benchmark/race_desc.py`
argues; no such run has been made here.

### All six cases, and the conversion generalises

`NS = NT = 40`, `NPSI = 200`, `NRBOX = 999`, trial `A0 = 1.0` in every case —
i.e. **nothing case-specific was tuned, and the amplitude closure was started
sixteen to three hundred times away from the answer on every one.** Three sweeps
each, every one converged.

| case | CHEASE `psi_ax` | the reference's | relative | `Ip` [A] | `q0` | vs `freegs4e`, rel L2 |
|---|---|---|---|---|---|---|
| `fixed-f-diiid` | 2.8978148570e-01 | 2.8977947450e-01 | **6.94e-06** | 9.981580e+05 | 1.031 | 5.6104e-05 |
| `fixed-d-tcv` | 1.7556945580e-02 | 1.7556774814e-02 | **9.73e-06** | 1.197507e+05 | 2.733 | 7.7360e-05 |
| `fixed-e-diamagnetic` | 4.3545915020e-02 | 4.3545313386e-02 | **1.38e-05** | 1.996227e+05 | 0.479 | 1.1503e-04 |
| `fixed-a-testtokamak` | 4.7789133020e-02 | 4.7788462591e-02 | **1.40e-05** | 1.996712e+05 | 1.864 | 8.9321e-05 |
| `fixed-h-circular` | 6.3545318050e-02 | 6.3544336033e-02 | **1.55e-05** | 2.998698e+05 | 0.377 | 1.1319e-04 |
| `fixed-g-mastu` | 8.0433203720e-02 | 8.0431769680e-02 | **1.78e-05** | 5.989861e+05 | 1.024 | 1.6513e-04 |

**THE `psi_ax` COLUMN IS PROBABLY MEASURING THE REFERENCE AND NOT CHEASE, AND
SAYING SO IS THE POINT.** On `fixed-h-circular` — the one case where CHEASE's own
convergence was measured — CHEASE's `psi_ax` is settled to 2.7e-08 and the whole
of its 1.55e-05 from the reference belongs to the reference. The other five rows
sit in the same 7e-06 to 1.8e-05 band without that check having been run, so the
honest reading is that **all six are consistent with the reference carrying an
error of about 1e-05 in `psi_ax`**, and not that CHEASE is accurate to 1e-05.

**`fixed-d-tcv` is the row worth looking at twice.** It is the most elongated of
the six at `kappa = 1.75` and it is the case
`tools/desc-benchmark/README.md` has an open problem about: DESC reads 2.17e-03
there at `M = 12` and does not converge monotonically in `M`. CHEASE reads
9.73e-06 with no special handling, and neither does MEQ have any trouble with
it. **Elongation is hard for a truncated spectral boundary and is not hard for
either of the codes that take the curve as points.**

---

## What is not done

**The convergence matrix is one case.** `fixed-h-circular`, because it is the
easiest and because `tools/desc-benchmark` developed on it. The other five have
one CHEASE solve each and no ladder. Running the matrix on them is
`converge.py <stem>`; it was not timed, and on this machine the one that was
run took tens of minutes.

**The floor was not reached and the instrument that stops it is known.** See
*The floor*. `cheaserun.py --match-box` puts the EQDSK box exactly on the case's
own mesh rectangle so that no interpolation is needed at all; it is **written and
never run**, and it would want `NRBOX = GridNR` as well, which `converge.py` does
not wire up.

**No timing of any kind.** Everything here was measured under load. A race in
the shape of `race_desc.py` — cold wall clock at matched accuracy, with a stated
definition of "matched" — has not been written, let alone run. CHEASE's cost is
also badly distorted here by `NRBOX = 999`: it writes a 999² ASCII box **twice**
on every sweep, which is a large and unmeasured share of the numbers in *Cost*.

**CHEASE's own test suite was not run.** `make test_chease` wants
`matlab -nodesktop`.

**The `NSTTP = 5` route is unused.** CHEASE can be driven by a `q` profile
instead of `TT'`, which is what `CLAUDE_INVERSION.md`'s `q( psi )`-driven solve
does in MEQ, and that would be a second and independent comparison of the same
two codes on a different problem statement. Nothing here touches it.

**Nothing here writes a plot.** `tools/plot_equilibrium.py` reads MEQ's two
NetCDF formats; the CHEASE arm writes the reference's `.npz` layout and nothing
draws it — the same gap `tools/desc-benchmark` records.

---

## Running it

```sh
cd tools/chease-benchmark          # from the MEQ checkout
V=../freegs4e-benchmark/venv/bin/python

$V convert_chease.py fixed-h-circular --check         # the posing, audited
$V cheaserun.py      fixed-h-circular --ns 40 --nt 40 # one CHEASE solve
$V compare_chease.py fixed-h-circular \
      --chease runs/fixed-h-circular-chease-NS40NT40NPSI100.npz \
      --meq /tmp/x/fixed-h-circular.nc
$V converge.py       fixed-h-circular --mapping-control
```

`CHEASE_EXE` overrides the binary's path.

**CHEASE reads `EXPEQ` and `chease_namelist` from the working directory by
hard-coded name and writes twenty-odd files beside them**, so every run gets its
own scratch directory and `cwd` is how it is told anything. `cheaserun.py`
deletes the scratch unless given `--keep`.

## Being a good citizen of this machine

The CHEASE arm is `nice -n 10` by default (`cheaserun.run_once`). The MEQ arm
goes through `race_desc.run_meq`, which pins to the physical cores with
`taskset` — that bounds it but does not deprioritise it, so prefix the whole
thing with `nice -n 10` when anything else is running, and then do not read the
seconds at all.
