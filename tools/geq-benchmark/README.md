# The `geq` benchmark: MEQ's rotating source against an independent implementation

`../geq` drives `../freegs4e`'s `ProfilesCentrifugalMirror` for **rotating
magnetic mirrors** — Python, finite differences on a uniform `(R, Z)` grid,
Picard with adaptive blending. It implements Abel et al. 2013,
Rep. Prog. Phys. 76 116201, **eq (136)** closed by its **(96)** and **(97)**:
the same paper, from the same equations, as `meq::RotatingSource`. Two
independent implementations of one equation is a rarer thing than it sounds,
and this directory is what spends it.

**→ [M-167](../../MEASUREMENTS.md#m-167)** is the measurement.
`tests/convergence/GeqSourceComparison.cpp` is the acceptance and it is a
registered ctest.

## What is compared, and what is not

**The SOURCE, and nothing else.** No solver runs on either side and no mesh
exists. The export writes a *prescribed* state — a vacuum-field mirror with
species, rotation and the quasineutrality potential solved — and the test
evaluates MEQ's `F` at geq's own state against geq's own `mu0 R Jtor`. A
disagreement between two equilibrium codes is three things at once: the source
each was given, the discretisation it solves with, and the iteration that got
there. This separates the first from the other two exactly, which is the same
reason `tools/freegs4e-benchmark/source_check.py` exists for the MHD source.

**The two forms agree analytically**, so a disagreement would mean something.
geq writes

```
mu0 R Jtor = mu0 R^2 Sum_s n_s [ T_s dlnN_s/dpsi + ( chi_s + 1 ) T_s dlnT_s/dpsi ]
           + mu0 R^2 Sum_s m_s n_s R^2 omega domega/dpsi
```

and MEQ writes `F = mu0 R^2 dp/dpsi|_r + g g'` with `p = Sum_s n_s T_s`.
Expanding MEQ's derivative reproduces geq's bracket term for term once the
`Sum_s Z_s n_s d( e phi_0 )/dpsi` piece is dropped — and it drops by
quasineutrality, the cancellation `docs/rotation.rst` records and the reason
MEQ's residual needs `phi_0` but never its derivative.

**Not compared**: `g g'`, which is identically zero on every geq path, so this
checks the pressure term alone; anything about the discretisation, the
boundary, or the iteration; and geq's **anisotropic and kinetic closures**
(`kinetic_closure.py`, `fast_species.py`), which MEQ has no counterpart for.
This is pinned to the isotropic `ProfilesCentrifugalMirror` path and must stay
there.

## The gauge is the conversion, and it is not a detail

| | where `phi_0` is pinned | what the density table means |
|---|---|---|
| MEQ | `phi_0( R_ref, psi ) = 0`, **one condition per flux surface** | `n_s0( psi )`, the physical density on the curve `R = R_ref` |
| geq | **one global point**, `psi_n = 0.5` on the midplane | `N_s( psi )`, in that gauge |

`N_s` absorbs the difference through `exp( Z_s e delta( psi )/T_s )`, which is a
**different factor for every species** — so the two codes' tabulated densities
mean different things and no single rescaling relates them. Compare `n_s(R,z)`
and `F`; **never** `N_s( psi )`.

**The transfer is geq's own quasineutrality, solved at `R = R_ref`.** For each
`psi`,

```
Sum_s Z_s N_s( psi ) exp( [ m_s omega^2 R_ref^2/2 - Z_s e Phi ]/T_s ) = 0
```

has a unique root — the sum is strictly decreasing in `Phi` — and the densities
it produces **are** MEQ's `n_s0`. This needs no interpolation of geq's `phi_0`
field, works at any `R_ref`, and satisfies `Sum_s Z_s n_s0 = 0` to root-find
precision, which is what `meq::RotatingSource`'s constructor checks.

**It is not circular.** The transfer fixes the densities on ONE curve; what is
then compared is MEQ's evaluation everywhere else — its own `phi_0( R, psi )` at
every other radius, the `psi`-derivative chain through `N`, `T` and `omega`, and
the `C( psi )` closure. A shared sign error would cancel at `R_ref` and survive
off it, which is where the comparison lives.

**Its `psi`-derivative is by implicit differentiation, not by differencing.**
`F` is `mu0 R^2 dp/dpsi`, so the derivative column of the density table is not a
convenience — it *is* the thing being compared.

**And the raw tables are exported too**, as `Nraw_<name>.dat`, so
`theGaugeTransferIsLoadBearing` can assert that handing them to MEQ is wrong.
It fails as a **refusal** rather than as a bad number: `Sum_s Z_s N_s` is
9.999994e-01 of the density scale at `R_ref` and the constructor throws. Without
that assertion, a transfer that happened to be near-identity would look like a
passing comparison.

## Files

| | |
|---|---|
| `export_mirror.py` | the geq side. Builds the state, does the gauge transfer, writes the reference |
| `reference/two/` | electron + deuterium at **M = 6.0**, the mirror operating point. MEQ takes `Closure::ClosedForm` |
| `reference/three/` | + carbon at `Z = 6`, at **M = 0.76**. MEQ takes `Closure::RootFind`, and the charge weighting stops being a plain sum |

**THE TWO MACH NUMBERS ARE FORCED, NOT CHOSEN**, with
`M² = m_i ω² R( z=0 )²/T_e` on the `ψ_n = 0.5` surface and `meta.txt` recording
it per case. The gauge transfer carries
`exp( m_s ω²( R_ref² − R_mid² )/2T_s )`, an exponent **linear in species mass**:
at `M = 6` that is +31.2 to −51.8 for deuterium and **+187.2 to −310.5 for
carbon**, a dynamic range of 1e216 against a double's `exp( 709.8 )` ceiling. The
three-species transfer loses its root-find bracket and **cannot be posed at
`M = 6`**. That is MEQ's own interface rather than anything about the comparison
— `MIRROR-PLAN.md` item 6 — and `R_ref` is already at its range-minimising value,
so it is not a tuning failure.

Per case: `meta.txt` (species constants, `R_ref`, `mu0`, the flux range),
`omega.dat`, `T_<name>.dat`, `n_<name>.dat` in MEQ's own three-column
`psi f f'` profile format, `Nraw_<name>.dat`, and two point sets —

| | |
|---|---|
| `points.dat` | `R z psi F` at geq's own `phi_0`, exactly as it ships |
| `points_tight.dat` | the same with `phi_0` solved pointwise to root-find precision |

**THE SECOND ONE EXISTS BECAUSE THE FIRST HAS A FLOOR THAT IS NOT MEQ'S.**
geq converges its `phi_0` *field* to about 5.5e-05 rms in `Sum_s Z_s n_s`, which
is ample for a Picard equilibrium solve and is not ample as a reference for a
source that solves the same condition to 1e-14 pointwise. Measured, that
residual alone moves **geq's own** `Jtor` by 1.324204e-03 at `M = 6`, against
MEQ's disagreement of 1.324229e-03 — four figures, with nothing of MEQ's in the
first number, and the same to six figures at `M = 0.76`. Against the tightened
reference the agreement is **4.45e-07** at `M = 6` and **6.17e-09** at `M = 0.76`,
and both fall with the table resolution, so neither implementation has an error
floor of its own.

**geq's own `φ₀` residual is 600x larger at `M = 6`** — a stiffer quasineutrality
problem is one its field solve converges less far on — which is why the loose
bound is per case in the test.

Nothing of MEQ's enters `points_tight.dat`: the species model, the profiles and
the `Jtor` expression are all geq's, and only the scalar condition both codes
state identically is solved harder.

## Regenerating the reference

Needs `../geq`, `../freegs4e` and the benchmark venv. The committed reference
means the **test** needs none of them.

```sh
cd /home/ian/projects/meq
PYTHONPATH=/home/ian/projects/geq/src \
tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_mirror.py \
    --case two --knots 2049 --omega-axis 4.758e6 \
    --out tools/geq-benchmark/reference/two

PYTHONPATH=/home/ian/projects/geq/src \
tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_mirror.py \
    --case three --knots 1025 --out tools/geq-benchmark/reference/three
```

`geq`'s own `CLAUDE.md` names a conda env and `/Users/...` paths; that is
another machine's recipe. Here both packages import from the venv at
`tools/freegs4e-benchmark/venv` with `geq/src` on `PYTHONPATH`, and geq's own
`tests/test_jtor_abel.py` passes. `pint` and `xarray` are absent, so `geq.utils`
does not import; nothing in the (136) path needs it.

**`--knots` is a genuine knob and the sweep is in M-167.** At `M = 0.76` the
error falls about 4x a doubling; at `M = 6` a clean **8x**, third order, because
the density table spans 22 orders and the interpolant has to work for it. The
committed 2049 gives 4.45e-07 at `M = 6` and 4097 gives 5.59e-08 at twice the
file size. **High Mach does not break the agreement, it makes the tables work
harder.** Do not read a change in that number as a change in the physics without
re-running the sweep.

`MEQ_GEQ_REFERENCE` points the test at another directory, which is how the sweep
is taken without touching what is committed.
