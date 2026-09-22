# MEQ against DESC, on the six shipped fixed-boundary machine cases

> **NO TIMING IN THIS FILE IS PUBLISHABLE.** Every second here was taken on a
> machine carrying another agent's `ctest` and its own solves, at a one-minute
> load average of 8 to 20 against 8 physical cores. They are here to fix the
> order of magnitude and to show the harness runs. The publishable numbers come
> from one command on an idle machine:
>
> ```sh
> ./venv/bin/python race_desc.py \
>     fixed-h-circular fixed-a-testtokamak fixed-e-diamagnetic \
>     fixed-g-mastu fixed-f-diiid \
>     --mode accuracy --self-consistent 12 \
>     --require-quiet --repeats 3 --out /tmp/race-desc
> ```
>
> **`--self-consistent 12` IS THE PROTOCOL AND NOT AN OPTION.** It splits each
> DESC arm in two. **Phase A, UNTIMED:** iterate the flux-label map to its fixed
> point — arriving at the pressure and the field as functions of the toroidal
> flux is *posing*, and posing is not solving, the same line `race.py` draws by
> leaving MEQ's meshing outside its clock. **Phase B, TIMED:** one cold forward
> solve, fresh process, on the converged profiles. That is the number.
>
> Without it the race is between three codes: posed from the reference, DESC's
> input carries freegs4e's own interior surfaces, and that flatters DESC by
> 2.9× — `psi_ax` 7.402e-05 from the reference against 2.129e-04
> self-consistent, on `fixed-h-circular` at `M = 10`. The posing is **refused
> rather than timed** if the fixed point is not reached in `N` sweeps;
> `race_desc.py` prints `POSING FAILED ... -- not timed`.
>
> `--require-quiet` refuses to time anything while the load is above 0.5. **It
> has not been run.** The ERROR columns below do not depend on load and are
> reproducible as they stand.

`examples/fixed-*.toml` are six real machines posed as fixed-boundary
Grad-Shafranov problems on a smooth MXH curve — `MEASUREMENTS.md`
[M-147](../../MEASUREMENTS.md#m-147) is what MEQ reports on them and how far it
sits from the `freegs4e` solve they were generated from. This directory races
them against **DESC**, a JAX-based inverse spectral equilibrium code from
Princeton (`/home/ian/projects/DESC`).

**DESC IS NOT ANOTHER GRAD-SHAFRANOV SOLVER AND THAT IS THE WHOLE INTEREST OF
THE COMPARISON.** `freegs4e` and `freegsnke` differ from MEQ in the
discretisation and in the algorithm; DESC differs in the *unknown*. It does not
solve for `psi( R, z )` at all — it solves for the SHAPES of the flux surfaces,
`R( rho, theta )` and `Z( rho, theta )`, with nestedness built into the
representation rather than emerging from the answer, and it drives the residual
of the force balance `J x B = grad p` rather than of `Delta* psi = -F`. So this
race is not "two ways of doing the same sum"; the two codes do not even carry
the same variables, and nothing about it can be read off from M-147.

---

## What is installed, and the exact commands

Python 3.14.4 is the only interpreter on this machine, and DESC's `setup.py`
claims 3.14 in its classifiers, so no separate interpreter was needed. Wheels
exist for every dependency including `jax-finufft`, so nothing was compiled.

```sh
# 1. The DESC checkout is READ-ONLY.  Copy what pip needs out of it rather than
#    installing -e, which would write an egg-info into somebody else's tree.
#    versioneer falls back to the PARENT DIRECTORY NAME for the version when
#    there is no .git, hence the directory name.
S=/tmp/desc-src
mkdir -p $S && cd /home/ian/projects/DESC
tar cf - desc setup.py setup.cfg pyproject.toml requirements.txt \
         MANIFEST.in README.rst LICENSE versioneer.py CITATION.cff \
  | ( mkdir -p $S/desc-0.17.3+8.g7582fb387 \
      && cd $S/desc-0.17.3+8.g7582fb387 && tar xf - )

# 2. Its own venv.  NOT tools/freegs4e-benchmark/venv, which is pinned tight
#    and which other work depends on.
python3 -m venv /home/ian/projects/meq/tools/desc-benchmark/venv
tools/desc-benchmark/venv/bin/python -m pip install --upgrade pip setuptools wheel
tools/desc-benchmark/venv/bin/python -m pip install $S/desc-0.17.3+8.g7582fb387
```

That installed `desc-opt 0.17.3+8.g7582fb387` with **jax 0.9.2 / jaxlib 0.9.2,
CPU only**, plus scipy 1.18.0, numpy 2.4.6, scikit-image 0.26.0 and the rest of
`requirements.txt`. No CUDA jaxlib: this machine's GPU is an RTX 2070 SUPER and
`CLAUDE.md`'s standing warning applies — consumer FP64 at 1/32 to 1/64 of FP32
makes a local device timing a measurement of the card and not of the code. DESC
prints `An NVIDIA GPU may be present ... Falling back to cpu` on every import,
which is the intended configuration here and not a warning to act on.

The version was **not** taken from PyPI on purpose: `/home/ian/projects/DESC` is
eight commits past `v0.17.3` and is the tree the comparison is about.

---

## The conversion, which is the hard part

### What the two codes are handed

MEQ takes the Grad-Shafranov free functions against its own **normalised
poloidal** flux — `src/meq/Source.hpp`, `meq::NormalisedMHDSource`:

```
F( r, z, psi ) = [ mu0 r^2 ( dp/dPsi )( Psi ) + ( g dg/dPsi )( Psi ) ] / psi_ax
Psi = psi / psi_ax,   Psi = 1 on the magnetic axis and 0 on Gamma
```

**`Psi` runs the opposite way to almost every other code's**, and the shipped
`-pprime.dat` header says so: `freegs4e`'s `psi_n` is 0 on the axis and 1 on the
separatrix, MEQ's `Psi` is 1 on the axis and 0 on `Gamma`, and
`Psi = 1 - psi_n/0.95` is the whole of the relabelling. Getting it backwards
converges perfectly well to a wrong equilibrium.

DESC takes a boundary surface, a total **toroidal** flux `Psi` in Webers, and
two profiles against `rho = sqrt( Phi / Phi_edge )`: the pressure `p( rho )` and
the **enclosed toroidal current** `I( rho )` in Amperes. Its poloidal flux comes
out as the compute quantity `chi`, in Wb per radian — the same units MEQ's `psi`
is in.

### The three things the conversion needs that MEQ's input does not carry

Each is a real asymmetry rather than a detail, and each is measured.

**1. The flux-label map.** `Phi( Psi )` needs the shapes of the *interior* flux
surfaces. MEQ computes them as part of solving; DESC has to be told, so the
conversion needs an equilibrium before it can pose the problem for DESC. The
default is the `freegs4e` 257² solve in `tools/freegs4e-benchmark/ref-n257`,
which is the same object M-147 measures MEQ against, so nothing new enters the
comparison — but it is not free, and `convert_desc.py --check` prints the bound:

```
THE OUTERMOST SURFACE: the reference's contour against Gamma itself
  area   contour 3.21595886e-01  MXH 3.21615233e-01  rel 6.016e-05
  int r  contour 3.23875825e-01  MXH 3.23895278e-01  rel 6.006e-05
  int 1/r contour 3.27745423e-01  MXH 3.27765861e-01  rel 6.236e-05
```

i.e. **6.0e-05 relative**, which is the MXH fit residual (2.7e-05 m on a
0.319 m minor radius) seen as an area. That is a floor on the DESC arm that the
MEQ arm does not pay, because `Gamma` for MEQ **is** the MXH curve.

**2. The absolute toroidal field.** Grad-Shafranov sees `g` only through
`g dg/dpsi`, so MEQ's input fixes `g^2` up to an additive constant and **MEQ's
answer for `psi` does not depend on it**. DESC's `Psi` is the toroidal flux,
which is proportional to `g`, so DESC must be given the constant. It comes from
the reference's own `fpol` at the boundary surface. The reconstruction is
checked rather than assumed: integrating MEQ's `gg'` table inward from that one
edge value gives `g` on the axis as **1.054355** against the reference's own
`fpol[0] = 1.0543556` — six figures, from a completely different route.

**3. The shape of the profiles in `rho`.** DESC wants power series. The refit is
a truncation and `convert_desc.py` prints its relative residual: at twelve even
powers, **5.8e-04 for `p` and 2.6e-06 for `I`** on `fixed-h-circular`. Even
powers only, because `Phi ~ psi_n` near the axis so `rho ~ sqrt( psi_n )` and
every flux function is smooth in `rho^2`; an odd term would put a `|rho|` cusp
on the axis that no equilibrium has.

### Why the flux integrals are contour integrals

Both integrands are of the form `a( Psi ) r + b( Psi )/r`, so

```
int_{A(u)} h( psi_n ) dA  =  h( u ) V( u ) - int_0^u h'( v ) V( v ) dv
```

with `V` the enclosed moment of `r` or `1/r`, and each moment follows from the
bounding contour alone by Green's theorem:

```
int r dA = oint ( r^2 / 2 ) dz          int ( 1/r ) dA = oint ln r dz
```

`h'` comes from the profile's own spline, so nothing differentiates `V`
numerically. **A hard mask over grid cells is `O( h )` in the boundary** and
reads the enclosed current 3.9e-04 wrong on the shipped 257² reference; the
contour route is `O( h^2 )` in the contour position. Both are computed and
printed by `--check`, and the cell sum reproduces `freegs4e`'s own `Ip` to
machine precision — which is how you know the 3.9e-04 is the *mask's* error and
not the contour's, since `freegs4e` computes its current with exactly that mask.

### The contour chords the surface, and that was the whole error

**THE SINGLE BIGGEST ITEM IN THIS CONVERSION, AND IT IS THE INSTRUMENT RATHER
THAN THE PHYSICS.** `find_contours` walks a marching-squares polygon with
linearly interpolated crossings, so a convex flux surface comes out as an
**inscribed** polygon and every enclosed moment reads LOW. The error is
`O( h^2 )` per segment, but the number of segments falls with the surface, so
it is worst on the innermost surfaces — exactly the ones the `rho` label near
the axis is built from. Measured on `fixed-h-circular`, with the field bicubic-
upsampled by `factor` before any contour is taken:

| factor | 1 | 2 | 4 |
|---|---|---|---|
| `int r dA` inside `psi_n = 0.1` | 1.760669e-02 | 1.763618e-02 | 1.764379e-02 |
| `I` at `Gamma` | 2.9974710e+05 | 2.9983552e+05 | 2.9985773e+05 |

Second order, converging from below, and the last column agrees to **4.6e-06**
with the same quantity taken off a MEQ `.nc` whose grid is 3.4× finer over the
plasma. **Two different fields agreeing is what says this is the instrument and
not the equilibrium**, and it is why the fix is `--refine` rather than "use a
better reference".

**WHAT IT IS WORTH IS THE POINT.** `psi_ax` follows the enclosed current almost
one for one, so on `fixed-h-circular` at `M = 12` the DESC arm moves

```
    --refine 1     psi_ax 6.351764849e-02     4.200e-04 from the reference
    --refine 4     psi_ax 6.354379868e-02     8.456e-06
```

**a factor of fifty**, for four hundred-odd spline evaluations. Every other
item in this conversion is one to two orders of magnitude smaller, and every
one of them was measured *first*, on a run where this one dominated — see
*hypotheses measured and killed*.

### The poloidal angle runs the other way

DESC's `theta` is **clockwise** in the `( R, z )` plane, so that `( rho, theta,
zeta )` is right-handed with `zeta = phi`. Its own default torus is
`R = 10 + cos theta`, `Z = -sin theta`, and `desc/examples/DSHAPE` carries
`Z1 = -1.47` at `m = -1` beside `R1 = +1.0` at `m = +1`. MXH's `t` is
counterclockwise. The conversion samples `Gamma` at `t = -theta`; getting it
backwards produces a left-handed surface that DESC **silently flips for you**,
after which its `theta` is not the one you fitted.

### The current's sign is not determined, and that is correct

Grad-Shafranov sees only `g dg/dpsi`, so the mirrored equilibrium solves MEQ's
problem equally well and MEQ's input cannot say which helicity the machine has.
Measured: `--current-sign +1` and `-1` give the same `nit`, the same `cost` and
the same `|chi( 1 )|` to every printed digit, with `iota` on the axis
`+2.6467` and `-2.6467`. So `psi` is the same function either way and the race
is unaffected. The reference's actual helicity is recoverable from its `Ip` sign
if anything ever needs it.

---

## What is NOT exact, in one list

Sizes on `fixed-h-circular`, with the conversion as it now stands.

| | size | who pays it |
|---|---|---|
| MXH fit of the reference's `psi_n = 0.95` contour | 2.7e-05 m, **6.0e-05** in enclosed area | the flux-label map, so DESC only |
| contour chording, at `--refine 4` | **~5e-06** in the enclosed current, from the factor-of-4 to factor-of-2 gap | DESC only |
| profile refit onto even powers of `rho` | `p` value **2.9e-04**, `p` DERIVATIVE **7e-03**; `I` 2.5e-06 and 1.0e-04 | DESC only |
| the boundary truncated to DESC's own `M` | **1.3e-04 m at M = 6, 7.1e-05 at M = 12, 2.8e-06 at M = 20** | DESC only |
| `chi`'s radial quadrature on the default grid | 1.6e-06 — removed, the fine grid is free | DESC only |
| MEQ's `Gamma` | exact — it IS the MXH curve | — |
| the reference as a yardstick | M-147 puts MEQ 5.0e-05 to 1.1e-04 from it | both |

**FIVE OF THE SEVEN ROWS ARE DESC'S ALONE, AND THAT IS STRUCTURAL RATHER THAN
UNFAIR.** MEQ is handed the problem in the coordinates it solves in; DESC is
handed the same physics in coordinates that have to be derived. Any comparison
of a Grad-Shafranov code against an inverse-equilibrium code on a case *posed*
as Grad-Shafranov pays this, and the only cure is `--from`, below.

**The boundary row is the one to watch and it is not a conversion error at all.**
`Equilibrium.__init__` calls `self._surface.change_resolution( L, M, N )`, so
**DESC cannot be given a boundary finer than its own poloidal resolution.** The
MXH curves here carry twenty harmonics, so `M < 20` is solving inside a
different curve. MEQ takes `Gamma` as an analytic curve at every degree and
pays nothing for it.

**The profile row's second column is the live one.** A least-squares fit to `p`
VALUES does not control `dp/drho`, which is what DESC's force balance uses, and
the derivative error **stops improving at twelve terms and then gets worse** —
1.5e-02, 6.7e-03, 7.4e-03, 1.1e-02, 1.7e-02 at 8, 12, 16, 20, 28 terms. That is
the monomial basis losing conditioning, not the profile, so "fit it harder" is
not available. `--profiles hermite` is the alternative and matches value and
slope at 65 knots with the axis slope enforced to zero.

## Hypotheses measured and killed

Recorded because each one is the obvious next suspect and re-deriving them is
the expensive part. All on `fixed-h-circular`, all against `psi_ax`, and **all
taken while the contour-chording error above was still in** — so each reads
"4.2e-04, unmoved", which is the point.

| hypothesis | test | verdict |
|---|---|---|
| DESC's poloidal resolution `M` | 6, 8, 10, 12, 16, 20 | **converged by M = 12.** 1.091e-03, 6.293e-04, 4.833e-04, 4.198e-04, 4.191e-04, 4.224e-04 — flat from 12 on, while the boundary residual kept falling 7.1e-05 → 2.8e-06 m |
| DESC's radial Zernike order `L`, independently of `M` | `M = 12`, `L` = 12, 18, 24 | **no.** 4.200e-04, 4.255e-04, 4.257e-04 |
| DESC's stopping tolerances | `ftol` 1e-2, 1e-6, 1e-10 | **no.** 9, 19 and 40 iterations for a 3.2e-07 relative move |
| the even-power truncation of the profiles | 8 and 12 terms; and `--profiles hermite` | **no.** 4.133e-04, 4.200e-04, 4.174e-04 |
| `chi`'s radial quadrature | default 25 nodes vs 1025 | **1.6e-06.** Real, tiny, and now removed |
| the pressure being inconsistent with `p'` in the reference | `p` against the integral of `pprime` | **no.** Agrees to 8e-16; and MEQ's own table reproduces the reference's `pprime` to 3e-13 |
| the current's SIGN | `--current-sign` ±1 | **irrelevant, as it must be.** Same `nit`, same `cost`, same `\|chi(1)\|`, `iota` on the axis +2.6467 and −2.6467. Grad-Shafranov sees only `g dg/dpsi`, so MEQ's input cannot say which helicity the machine has and `psi` is the same function either way |
| the flux-label map's SOURCE | the reference against MEQ's own converged `k = 2, r = 1` answer | **yes, 4.200e-04 → 2.201e-05** — and then `--refine 4` on the reference reached 8.456e-06 **without** MEQ's answer, which is what identified the cause as the contour and not the equilibrium |

Two independent checks that the conversion is right at all, rather than merely
self-consistent:

* **The toroidal field on the axis, rebuilt from MEQ's `gg'` table alone.**
  Integrating `g^2( u ) = g^2( Gamma ) + 2 int ( g dg/dpsi ) dpsi` inward from
  the reference's single edge value reproduces the reference's own `fpol[ 0 ]`
  on all six machines to **1.3e-08 to 2.8e-07**. That exercises the `Psi` to
  `psi_n` relabelling, the `d/dPsi` to `d/dpsi` chain rule and the sign, and
  nothing in it is fitted.
* **`q` on the axis.** The converted `Phi( u )` and the flux label give
  `q = ( 1/2pi ) ( dPhi/du )/( dpsi/du )` = 0.381 just off axis, against DESC's
  own `1/iota( 0 ) = 0.378`, and 2.02 at `Gamma` against a cylindrical estimate
  of 1.76.

### Removing the reference from the posing

`--from <meq.nc>` takes the flux-label map from MEQ's own converged field
instead of the reference's, and **is implemented** — `convert_desc.field_from_meq`,
reached as `descrun.py --from` and `race_desc.py --map-from`. Measured on
`fixed-h-circular` at `M = 12`, before `--refine` existed: 4.200e-04 with the
reference's map against **2.201e-05** with MEQ's. With `--refine 4` in, the
reference's own map reaches 8.456e-06 without MEQ's answer at all, and the two
maps agree to 4.6e-06 in the enclosed current — which is what identified the
cause as the contour and not the equilibrium.

ITERATING it on DESC's own field is the fixed point that removes both codes'
answers from the posing, and it **is implemented**:
`descrun.py --self-consistent N`. Pose from anything, solve, rebuild the map
from that field, re-pose, re-solve. What the limit depends on is the MXH curve,
the two profile tables and the two scalars Grad–Shafranov does not fix — which
is exactly the problem MEQ is handed, and nothing else. It is **untimed**, for
the reason the header gives: posing is not solving.

**IT DOES NOT CONVERGE UNDAMPED, AND THAT IS MEASURED RATHER THAN GUARDED
AGAINST.** Replacing the map outright — `--self-consistent-mix 1` — on
`fixed-h-circular` at `M = 10`:

| sweep | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| `psi_ax` | 6.353780e-02 | 6.352345e-02 | 6.353595e-02 | 6.352159e-02 | 6.353659e-02 |
| field change | — | 1.251e-04 | 8.574e-05 | 1.058e-04 | 1.076e-04 |

a **period-2 limit cycle** of amplitude 2.3e-04 relative, with the change
plateauing at 1e-04 and refusing to fall. `Psi_total` converges perfectly
throughout — 3.3051345e-01 to 3.3051059e-01 — so it is the axis flux that
oscillates, which is the ordinary marginal instability of a Picard iteration on
a self-consistent quantity. At the default `--self-consistent-mix 0.5` it is
monotone and stops in eight sweeps:

    1.251e-04 -> 4.321e-05 -> 5.763e-06 -> 2.151e-06
             -> 1.523e-06 -> 1.006e-06 -> 6.659e-07

**AND THE SELF-CONSISTENT ANSWER IS FURTHER FROM THE REFERENCE, WHICH IS THE
POINT.** Same case, same `M = 10`, `psi_ax` against the reference's:

| posing | relative |
|---|---|
| from the reference | **7.402e-05** |
| self-consistent | **2.129e-04** |

2.9× worse — because the first number was partly circular. freegs4e's own
interior surfaces were in DESC's input, so DESC was being scored against a code
that had helped write its problem. **Quote the self-consistent number for "two
independent codes on the same problem" and the other one for nothing.**

---

## What it actually reports — ROUGH numbers, on a machine that was not quiet

**EVERY SECOND BELOW IS ROUGH AND NONE OF THEM IS PUBLISHABLE.** They were taken
with another agent's solves running: the one-minute load average was **12.6**
when the race below ran, against a physical core count of 8. `CLAUDE.md` records
a 490 s / 540 s spread on identical code on an idle machine, and this is far
worse than that. **Re-run with `--require-quiet` before quoting any of it.** The
ERROR columns do not care about load and are reproducible; the seconds are here
only to show the harness works and to fix the order of magnitude.

### The six cases, DESC at `M = 12`, against the reference's own `psi_ax`

| case | DESC `psi_ax` | reference | relative | MEQ at its shipped `k = 2`, from M-147 | DESC `nit` | boundary residual |
|---|---|---|---|---|---|---|
| `fixed-f-diiid` | 2.897785899e-01 | 2.897794745e-01 | **3.05e-06** | 6.99e-06 | 44 | 4.31e-04 m |
| `fixed-h-circular` | 6.354367862e-02 | 6.354433603e-02 | **1.04e-05** | 1.14e-04 | 9 | 7.14e-05 m |
| `fixed-e-diamagnetic` | 4.354609459e-02 | 4.354531339e-02 | **1.79e-05** | 1.26e-05 | 56 | 4.51e-05 m |
| `fixed-a-testtokamak` | 4.778661435e-02 | 4.778846259e-02 | **3.87e-05** | 1.12e-05 | 12 | 1.81e-04 m |
| `fixed-g-mastu` | 8.042850571e-02 | 8.043176968e-02 | **4.06e-05** | 1.79e-05 | 43 | 4.80e-04 m |
| `fixed-d-tcv` | 1.751865449e-02 | 1.755677481e-02 | **2.17e-03** | 9.41e-06 | 72 | 3.13e-04 m |

**FIVE OF THE SIX SIT WHERE MEQ SITS**, between 3e-06 and 4e-05 of a reference
neither of them can see. `fixed-d-tcv` is fifty times worse and it is the most
ELONGATED of the six at `kappa = 1.75`; `M = 12` is not enough poloidal
resolution for that shape, which its boundary residual of 3.1e-04 m on a 0.24 m
minor radius already says before any solve happens.

### `fixed-d-tcv` does not converge in `M`, and it is the open problem here

| `M` | 12 | 16 | 20 | 24 | 24, `ftol 1e-8`, 200 iterations |
|---|---|---|---|---|---|
| relative `psi_ax` | 2.17e-03 | 9.61e-04 | **8.52e-05** | **2.91e-02** | 6.32e-04 |
| `nit` | 72 | 51 | 100 | 49 | 200 (capped) |
| final cost | 3.04e-08 | 1.00e-08 | 5.25e-10 | 8.26e-08 | 2.78e-09 |
| `R_axis` | 0.9156 | 0.9148 | 0.9150 | 0.8971 | 0.9095 |

**NOT MONOTONE, AND THE `M = 24` ROW IS NOT A RESOLUTION EFFECT.** Its boundary
residual is 1.47e-06 m — the best of the five — and its cost is two orders
WORSE than `M = 20`'s while reporting `success True` and a nested solution with
the magnetic axis 2 per cent off. That is the optimiser stopping, not the
discretisation. Tightening to `ftol 1e-8` and two hundred iterations recovers
most of it and still does not reach `M = 20`.

**SO THE TOLERANCE RESULT ABOVE IS CASE-SPECIFIC AND MUST NOT BE GENERALISED.**
On `fixed-h-circular` the default `ftol = 1e-2` and `ftol = 1e-10` agree to
3.2e-07 and the loose one takes a quarter of the iterations; on `fixed-d-tcv` at
`M = 24` the same default lands 2.9e-02 away. What a stopping rule is worth is a
property of the problem, and the honest reading is that **DESC's default
tolerances are adequate on five of these six cases and are not adequate on the
most elongated one.** A race run on TCV should not be published until this is
understood; the other five stand.

### The race, `fixed-h-circular`, error as relative L2 of `psi` over the 257² grid

| code | resolution | seconds | warm | vs `freegs4e` | vs the other code's best |
|---|---|---|---|---|---|
| MEQ | `k=1 r=0` | 7.25 | — | 1.0011e-03 | 1.0050e-03 |
| MEQ | `k=2 r=0` | 8.69 | — | 9.7266e-05 | 1.1797e-04 |
| MEQ | `k=3 r=0` | 12.08 | — | 1.1808e-04 | 1.3647e-04 |
| DESC | `M=8` | 74.10 | 1.37 | 1.3377e-03 | 1.3649e-03 |
| DESC | `M=12` | 73.34 | 1.57 | 9.3225e-05 | 1.3934e-04 |
| DESC | `M=16` | 82.41 | 3.16 | 5.8923e-05 | 1.1797e-04 |

**MEQ's `k = 2, r = 0` row reads 9.7266e-05 against M-147's published
9.712e-05**, which is the harness checking itself: the MEQ arm is the shipped
driver on the shipped file through `compare.py`'s own norm, so it has to
reproduce the published table and it does.

**THE TWO CODES AGREE TO 1.18e-04 OVER TWELVE THOUSAND NODES, AND THAT IS THE
FLOOR RATHER THAN AN ACCURACY.** It is the same size as each code's own distance
from the reference, so nothing in this comparison resolves below it and the
harness refuses targets that try. Pushing it down means a better reference, not
a better solve on either side — the MXH fit residual is 6.0e-05 in enclosed
area and M-147 already identifies it as what floors the `freegs4e` comparison.

**MEQ's `k = 3` row is WORSE than its `k = 2` row and that is the instrument.**
The reference is a 257² finite-difference free-boundary solve; once a MEQ run is
nearer to the true answer than the reference is, refining it moves it AWAY from
the reference. M-147's own order table avoids this by measuring against MEQ's
own `r = 3` answer instead. The same caution applies to every "vs freegs4e"
number here below about 1e-04.

**THE VERDICT DEPENDS ENTIRELY ON WHICH CLOCK, WHICH IS THE POINT OF HAVING
TWO.** At comparable accuracy — MEQ `k=2 r=0` against DESC `M=12` — MEQ is
**8.4× faster cold** and DESC is **5.5× faster warm**. Neither number is wrong
and neither is the answer on its own: the cold one is what solving one
equilibrium once costs, the warm one is what the second and subsequent solves of
the same shape cost in one process, and MEQ has no equivalent split because it
compiled at build time.

## Running it

```sh
cd /home/ian/projects/meq/tools/desc-benchmark

./venv/bin/python convert_desc.py fixed-h-circular --check    # the posing, audited
./venv/bin/python descrun.py fixed-h-circular -M 12      # one DESC solve
./venv/bin/python compare_desc.py fixed-h-circular \
      --desc runs/fixed-h-circular-desc-M12.npz --meq /tmp/x/fixed-h-circular.nc

./venv/bin/python race_desc.py fixed-h-circular --mode accuracy --require-quiet
./venv/bin/python race_desc.py fixed-h-circular --mode wall --match dofs
```

`race_desc.py` runs **both** races. Its header carries the argument for each and
the definition of "matched", which is a judgement and is written down as one.

---

## Two clocks for DESC, and quoting one of them is how this goes wrong

DESC is JAX: the first solve of a given *shape* compiles the residual, its
Jacobian and the optimiser loop, and the second does not. Measured on
`fixed-h-circular` at `M = 10`, one process:

* **build 11.6 s + solve 15.7 s = 27.3 s cold**
* **0.75 s warm**, the same shape solved again in the same process

Neither is "the" answer. `cold` is what a user pays to solve one equilibrium
once and is the number comparable with MEQ's whole driver run; `warm` is what a
parameter scan or an optimiser inner loop pays. MEQ has no equivalent split — it
compiles at build time — so the race times DESC as a **fresh process**, exactly
as it times MEQ's driver, and reports `warm` beside it rather than instead of
it.

**AND THE TWO CLOCKS GIVE OPPOSITE VERDICTS, WHICH IS WHY NEITHER MAY BE QUOTED
ALONE.** At comparable accuracy on `fixed-h-circular` — MEQ `k = 2, r = 0` at
9.73e-05 against DESC `M = 12` at 9.32e-05 — MEQ is **8.4x faster cold** and
DESC is **5.5x faster warm**. Both are rough, both were taken under load, and
both will survive a quiet re-run as ratios of roughly that size because the
mechanism is not a contention effect: a JAX compilation either happened in the
clock or it did not.

**DESC'S CPU THREADING IGNORES BOTH THREAD VARIABLES, WHICH IS ITS OWN TRAP.**
A DESC solve with `OMP_NUM_THREADS = MKL_NUM_THREADS = 8` was measured at
**1257% CPU** — twelve and a half cores. XLA's Eigen pool reads neither
variable. `race_desc.py` therefore pins BOTH arms with `taskset` to one logical
cpu per physical core, because pinning only one of them would be a measurement
of who ignored an environment variable.

---

## What is not done

**Two scalars that should be in the case files and are not.**
Grad–Shafranov sees `p` only through `dp/dpsi` and `g` only through
`g dg/dpsi`, so MEQ's input determines neither `p` nor `g` on Gamma — and
DESC needs both, `g` because its `Psi` is the toroidal flux and `p` because it
is a profile rather than a gradient. They are taken from the reference's
`pressure` and `fpol` arrays today, and `--p-boundary` / `--g-boundary`
override them, so the reference is not *required* for either. **What would
remove it properly is `p_at_gamma` and `g_at_gamma` in each
`<stem>-meta.json`**, beside the `pprime_at_gamma` and `ggprime_at_gamma`
already there — two numbers in
`tools/freegs4e-benchmark/make_case.py`. Until then a `fixed-*.toml` does not
describe a full MHD equilibrium, only a Grad–Shafranov one, and that is worth
knowing for any consumer and not only for DESC.

**The race's `--mode wall`.** Written, argued in `race_desc.py`'s header, and
**never run**. Neither of its two definitions of "matched" is privileged and
neither has been exercised.

**`--self-consistent` on the other five cases.** It is written and measured on
`fixed-h-circular` only. The damping constant is the default 0.5 and nothing
has been tuned; on a case whose map moves more it may want a smaller one, and
`fixed-d-tcv` should not be attempted until its `M` convergence is understood.

**`fixed-d-tcv`.** See its own section. Five cases stand, one does not.

**A quiet-machine run of anything.** Every second in this file is contended.

**`--mode wall`'s two matchings are implemented and unexercised.** `race_desc.py
--mode wall --match dofs` and `--match scale` run, but no matched-resolution
table has been taken, deliberately: a matched-resolution number taken under load
is worth less than no number, and the definition is a judgement that should be
argued over before it is measured. `--mode accuracy` is the one to believe and
this file's header says why.

**No plot.** `tools/plot_equilibrium.py` reads MEQ's two NetCDF formats; the
DESC arm writes the reference's `.npz` layout and nothing draws it.

## Being a good citizen of this machine

Both arms are pinned to the physical cores with `taskset`, which bounds them but
does not deprioritise them. Prefix the whole race with `nice -n 10` when
anything else is running — and then do not read the seconds, because the point
of `--require-quiet` is that nothing else should be.

---

## The free-boundary race, which is a different comparison and only one case reaches it

Everything above is the SIX FIXED-BOUNDARY machines, where MEQ is handed `Gamma`
and DESC is handed a surface fitted to it. `descfreeb.py` is the other
comparison: both codes solve for the plasma boundary, from a coil set.

```sh
./venv/bin/python convert_desc.py limited-tokamak-filament --check
./venv/bin/python descfreeb.py limited-tokamak-filament -M 10 --boundary reference
./venv/bin/python compare_desc.py limited-tokamak-filament \
      --desc runs/limited-tokamak-filament-descfreeb-M10-reference.npz \
      --meq ../../limited-tokamak-filament.nc
```

**ONE CASE, AND WHICH ONE IS NOT A CHOICE.** `BoundaryError` constrains the
LCFS, and the LCFS is a `FourierRZToroidalSurface` — a truncated Fourier series
in the poloidal angle. Every free-boundary machine in `examples/` but one is
DIVERTED, so its LCFS is a separatrix through an X-point and has a corner no
truncation turns; it is the same fact that forces the fixed ladder to
`psi_n = 0.95` rather than the separatrix. `examples/limited-tokamak-filament
.toml` is limited, so its edge is smooth, and it is the only free-boundary case
here.

**AND IT DOES NOT CONVERGE** → **[M-162](../../MEASUREMENTS.md#m-162)**.

This section used to publish a three-code agreement at 1.3e-04. It is
withdrawn. `descfreeb.py` never read `result["success"]`, so a run that stopped
at its starting point was reported as an answer — and a boundary seeded from the
reference's own LCFS starts at a stationary point, cannot be improved on by the
trust region, and therefore reports the reference back with a tiny error
**because it did not move**. The 1.3329e-04 and 1.2873e-04 rows both came from
that run.

Re-taken with a convergence check, **3 of 13 runs converge**, and the three land
at −2.7e-02, −3.9e-02 and −4.6e-02 — **none within 1e-02 of the reference**.
There is no resolution-converged free-boundary answer here to quote.

**USE THE FIXED-BOUNDARY LADDER.** It is monotone in `M`, agrees with MEQ at
9.35e-05, and is what the race numbers rest on — M-157 and M-163 section 3.
`tools/benchmarks/HOW_TO_DRIVE_DESC.md` is the short version of which arm to
trust and why.

### Three things this comparison needs that the fixed one does not

**1. The toroidal field, and it is the one most likely to be rediscovered.**
`BoundaryError` enforces `B_out . n = 0` **and**
`B_out^2 - B_in^2 - 2 mu0 p = 0`, and `B_in` carries `B_phi = g/R` while a set
of poloidal field coils produces none at all. A tokamak's toroidal field comes
from a TF coil that **no Grad-Shafranov input names** — GS sees `g` only
through `g dg/dpsi`, so MEQ's answer does not depend on the additive constant
and `[[coils]]` is the whole machine for MEQ and is not the whole machine for
DESC. With the PF coils alone the pressure residual starts at **0.90
normalised**, which is `B_phi^2` entire; the optimiser trades the good
normal-field residual against the impossible one and walks the boundary
**0.81 m** on a plasma of minor radius 0.34, for a `psi_ax` 92% wrong. It does
not fail. `g_at_gamma` supplies it, which is the same scalar §2 of *the three
things the conversion needs* already uses for the total toroidal flux.

**2. Profile amplitudes that are this reference's, not a sibling's.** Both this
machine's profiles are exactly `amplitude * Psi^2` at every grid, and MEQ
carries a profile scale that `[source] PlasmaCurrent` closes — so one common
factor on both tables is absorbed exactly and MEQ converges to the same
equilibrium from either pair. **The amplitudes are 6.03% apart**, not the 0.22%
the span would explain, because `freegs4e`'s profile amplitudes are an output of
its control system exactly as its coil currents are. DESC is handed `p( rho )`
and `I( rho )` as ABSOLUTE profiles, so the sibling's tables pose it a 6%
stronger plasma and it reads `enclosed current 3.192583e+05 A` against a target
of 3.0e+05 — with MEQ's own run saying nothing.
`tools/freegs4e-benchmark/make_freeb_profiles.py` writes the case its own pair,
and `--check`'s enclosed-current line is what catches this.

**3. A branch, because the objective has more than one minimum.** Started from
a circle of the machine's design size rather than from the reference's LCFS,
DESC converges in 5 iterations to a boundary **0.19 m** from the reference's
with `psi_ax` **15%** out — at a residual sum of squares of 4.760e-06 against
the reference boundary's own 1.369e-05. It is not under-converged: it found a
better minimum, at a different equilibrium. `--boundary reference` is the
matched default for that reason and `--boundary circle` is the cold number
beside it. MEQ needs the same thing and gets it from `[boundary.xpoint]`,
`[source] PsiAxis` and the initial guess; see `MEASUREMENTS.md` M-26.

### What is asymmetric here, and it is not the same list as the fixed ladder's

**The conductor model is EXACT on all three arms, which no fixed case can say.**
`freegs4e`'s `H_limited_circular` is `freegs4e.machine.Coil`, a point filament;
MEQ's `[conductors] Model = "filament"` is a point filament; DESC's
`FourierPlanarCoil` with a single `r_n` is a circular loop, checked against
`mu0 I R^2 / 2( R^2 + z^2 )^{3/2}` on the loop's own axis to nine figures. So
the row that floors the meshed sibling — a `( w/d )^2` term of about 1.6e-02 on
the near field — is absent.

**The MXH fit is still DESC's alone but it is ten times smaller**, 8.3304e-06 m
at 20 harmonics, **2.44e-05** of the minor radius, and it does not improve with
more harmonics: the floor is the 512-point LCFS polygon the 513² reference
writes, not the fit.

**AND THE TWO CODES ARE NOT GIVEN THE SAME STATEMENT OF THE PROBLEM.** MEQ is
given the coil currents, `p'( Psi )`, `g g'( Psi )`, a target `I_p` and **a
limiter contact**; DESC is given the coil currents, `p( rho )`, `I( rho )` and
**the total toroidal flux `Psi`** — there is no wall and no limiter in its
formulation at all. The plasma's size is pinned by the contact in one and by
`Psi` in the other. Both are complete statements of one equilibrium and neither
is the other's, so a disagreement can live in that difference as well as in
either code's discretisation. This is the row with no analogue above.

### Two harness defects this case found, both in the gauge

**`compare_desc.py` was applying the FIXED ladder's gauge shift to a
free-boundary MEQ run.** `compare.py` already knows the difference and takes
`free_boundary=` for it — a free-boundary MEQ run shares `freegs4e`'s gauge
exactly, both solving `psi -> 0` at infinity — and `meq_pair()` was not passing
it. Measured: **1.36e-01 with the shift against 1.4e-04 without**, on a case
whose scalars agree to 4.6e-05. `meta["free_boundary"]` is what selects it now.

**And a trap that is NOT a defect, recorded because it cost time here.** A MEQ
`.nc` written under `[conductors] Model` carries
`content = "remainder (psi - psi_c)"`, and the `psi( Z, R )` **grid variable is
the TOTAL**: `apps/meq.cpp` says so at the attribute's own site — it describes
*which field the coefficients are*. Under a split the lossy interchange format
carries the physically exact field and the exact format, the `.gf`, carries the
remainder. Adding `psi_c` back to the grid double-counts it, and on this machine
that reads 4.9e-01.
