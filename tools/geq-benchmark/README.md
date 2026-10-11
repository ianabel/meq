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
| `export_vacuum.py`, `reference/vacuum/` | the **vacuum** machine, no plasma, for MR-3's solve. See *The vacuum mirror* below |
| `export_rotating.py`, `reference/rotating/` | the **solved** free-boundary rotating mirror, for MR-4. See *The rotating mirror* below |

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
`psi f f'` profile format, `Nraw_<name>.dat`, and two point sets, which are
geq's OUTPUT and are regenerated rather than tracked —

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

Needs `../geq`, `../freegs4e` and the benchmark venv. **The point sets are
geq's answer and are not tracked**, so `GeqSourceComparison` needs this run
before it can read them; the tables it is given are tracked. Every reference in
this directory was generated at geq `73b9a72` and freegs4e `9a1ee0c`, both
clean; `export_rotating.py` records that in its `meta.txt`, and a regeneration
at any other commit is a different reference.

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

## The vacuum mirror: `export_vacuum.py` and `reference/vacuum/`

`MIRROR-PLAN.md`'s MR-3 compares MEQ's **solve** with geq's field, which is the
other half of what this directory does: no plasma, geq's two packs, MEQ's
half-disc with the exterior coupling live. `tests/convergence/VacuumMirror.cpp`
is the acceptance and is a registered ctest.

**WITH NO PLASMA geq'S FIELD IS A CLOSED-FORM SUM, NOT A SOLVE**, and the script
asserts it to round-off before writing anything: freegs4e's `Equilibrium.psi()`
is then the machine's Greens sum over the coils' own filament lists. So what MEQ
is compared with is a **coil model**, read off freegs4e rather than assumed:

| | |
|---|---|
| a pack | a `MultiCoil` with `R` and `Z` given as lists, so `turns = len( R )` and the `turns` geq passes is ignored on that branch. Each filament carries `current`, and the pack's total is `current * nfil` |
| the area | never enters the field; it is a current-density limit for the self-force. With `Nr = 1` the pack is **ten filaments at the radial centre**, `z = zc ± 0.05 … 0.45` — exactly `meq::filamentStack( coil, 1, 10 )`, the midpoint rule for the rectangle |
| `psi` | weber per radian on both sides: freegs4e's `Greens` carries `mu0/2pi` and says "Webers/2pi", `meq::filamentPsi()` is `R A_phi` |
| the sign | both packs carry the same current — a mirror, not a cusp |

`reference/vacuum/` holds `meta.txt`, the machine as geq computed it — pack radii,
length, centre, current, filament count, which the test **reads** rather than
re-deriving geq's tape-and-turn arithmetic — and `points.dat`, geq's output and
not tracked, `R z psi` on every other node of geq's 65 × 65 grid (`R ∈ [0.05, 0.63]`, `|z| ≤ 1.5`, the corners
and the midplane kept). Regenerate with

```sh
cd /home/ian/projects/meq
PYTHONPATH=/home/ian/projects/geq/src \
tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_vacuum.py \
    --out tools/geq-benchmark/reference/vacuum
```

and `MEQ_GEQ_VACUUM_REFERENCE` points the test elsewhere.

**THE AGREEMENT IS LIMITED BY geq'S COIL MODEL, AND THE TEST SAYS SO BY
REMOVING IT.** MEQ meshes the pack as the rectangle it is, uniform current
density over the winding; geq's ten filaments are that rectangle's midpoint rule.
The test measures geq's filaments in MEQ's own kernel (round-off: every
convention agrees), geq against `meq::coilPsi()` with no solve at all (the model
floor), MEQ against both, and finally MEQ handed geq's own twenty filaments by
subtraction — at which point the floor is gone. → **[M-260](../../MEASUREMENTS.md#m-260)**.

## The rotating mirror: `export_rotating.py` and `reference/rotating/`

`MIRROR-PLAN.md`'s MR-4 compares MEQ's free-boundary rotating mirror with
geq's, **in the mode geq solves in**: Picard, the twenty filaments FIXED, both
labels read every step at `( 0.10, 0 )` and `( 0.35, 0 )`, no current
constraint. Unlike `reference/two/` this is a solve, not a prescribed state.
MR-4b and MR-4d's ctests read it.

**IT IS NOT SHIPPED geq**, for the three facts `MIRROR-PLAN.md` records, and
the script changes nothing in `../geq` or `../freegs4e` — every change is a
monkeypatch from the script:

1. **The direct `freeBoundary`, not von Hagenow.** Shipped Hagenow's fixed
   `eps = 1e-2` m biases `ψ` by about 1% of the plasma flux, and its own
   self-convergence cannot see that. The shipped route is solved once, at
   257², as the CONTROL: a comparison that cannot read its bias cannot resolve
   an error that size.
2. **New profiles, flat at both labels to a stated order.** geq's Gaussians
   step there, which MEQ's support cannot take. Every profile — `T_e`, `T_D`,
   `N_e = N_D`, `ω` — is
   `f_edge + ( f_peak − f_edge )·C ψ_n³( 1 − ψ_n )⁴`, `C = 1/( (3/7)³ (4/7)⁴ )`,
   peaked at `ψ_n = 3/7`, and is evaluated **analytically** inside geq rather
   than through its cubic splines, whose not-a-knot ends would leave an
   `O( h³ )` slope at the labels.
3. **`φ₀` tightened inside every Picard step**, by a pointwise Newton on
   `Σ_s Z_s n_s = 0` after geq's global solve.

**THE SHAPE IS SET BY THE EDGE ORDER, AND IT IS ASYMMETRIC ON PURPOSE.** What
an edge costs MEQ is its order `j`, `p′ ~ d^j` at distance `d` inside a label:
`CLAUDE_FB.md`, *The plasma edge caps the order*, has `ψ*` keep `k + 2` only for
`k ≤ j`. The plan's first choice, `16 ψ_n²( 1 − ψ_n )²`, is `j = 1` and would cap
MR-4c's rates from `k = 2`. This shape is `j = 2` at the inner label and `j = 3`
at the outer, and `J_φ` on the converged state falls as `d^2.000` and `d^3.000`
there. The asymmetry is the second reason for it: every table differs from its
own reflection by 35–54% of its peak, so a reader that took them in `ψ_n`
rather than `Ψ` gets different numbers, and the sense of the abscissa —
`MIRROR-PLAN.md`'s item 4 — is something this reference can catch.
`--shape quartic` and `--shape sextic` (`64 ψ_n³( 1 − ψ_n )³`, `j = 2` both sides)
remain as options; both are symmetric and neither is the reference.

**AND `N_D` IS HELD.** geq rebuilds the deuteron's flux function every step from
the electrons' midplane density. Here `AllSpecies.set_flux_density` is a no-op
and the deuteron is flagged `_density_eval`, so both `N_s` are the prescribed
profiles and `φ₀` alone carries quasineutrality. With both `N_s` prescribed no
gauge freedom is left, so geq's one-point pin of `φ₀` (`ψ_n = 0.5` on the
midplane) becomes a wrong value at ONE node, and that node is why fact 3 is
needed every step and not just at the end: tightening only the converged state
leaves `ψ` at the untightened fixed point, which is 7e-05 away at 129² — the
size of the Hagenow bias, and far above geq's own resolution floor. Off that
node geq's global solve is already neutral to round-off. The script asserts
after every solve that geq has replaced none of the prescribed profiles.

**THE DENSITY IS SET IN geq'S GAUGE, SO ITS SCALE IS NOT conftest's.**
conftest prescribes the midplane density; here `N_s` is the density at zero
potential, and the midplane density is `N_s exp( m_D ω² R²/2( T_e + T_D ) )`, an
exponent of 2.2 at the peak. So `N_peak = 5e18` puts the physical midplane peak
near conftest's 5e19. And `ω_edge` is a fifth of `ω_peak`, because the same
exponent grows as `R²` against a cold edge: at half the peak, the outer edge's
density is hundreds of times the peak's and the first Picard step diverges.

| | peak, `ψ_n = 3/7` | edge, both labels |
|---|---|---|
| `T_e` | 1.28 keV | 128 eV |
| `T_D` | 1.53 keV | 153 eV |
| `N_e = N_D` | 5e18 m⁻³ | 5e17 m⁻³ |
| `ω` | geq-Mach 2 at `R = 0.240` | `ω_peak/5` |

`R = 0.240` m is the converged midplane radius of `ψ_n = 3/7`, rounded, so geq's
Mach there is 2.00 and M-167's is 3.09.

**THE TABLES ARE IN MEQ'S CONVENTIONS**: `Ψ = 1 − ψ_n` ascending, slopes
`d/dΨ`, `T_s` in joules, `ω` in rad/s; `n_<name>.dat` is the physical density on
`R = R_ref`, carried there by `export_mirror.py`'s own transfer, and
`Nraw_<name>.dat` is geq's prescribed `N_s`. `R_ref = 0.213` m is the radius that
minimises the largest `| ln( n_s/n_s0 ) |` over the plasma's nodes, at ±1.16;
the midplane radius of the peak, 0.240, gives −1.60 to +0.98. Evaluated through
MEQ's own gauge from the tables, `ln n_D` matches geq's at random `( R, ψ_n )` to
1e-11. The 2049 knots are the fewest power-of-two count whose Hermite
interpolant is within 1e-10 of every table's values, pointwise and relative;
the slopes, which `F` is built from, are within 2.6e-09 of their largest. `meta.txt`
records the measured error of each.

**THE PROBLEM IS TRACKED AND geq'S ANSWER IS NOT.** Tracked, because they are
what MEQ is given: `meta.txt` — provenance, every setting, `R_ref`, the two
points, the machine under `export_vacuum.py`'s keys and the species under
`export_mirror.py`'s — the tables, and `filaments.dat` (`R z current`, the
twenty fixed filaments). With both `N_s` prescribed, the transfer to `R_ref` is
pointwise quasineutrality on the prescribed profiles, so the tables do not
depend on geq's solve; `R_ref` was chosen on the solved state and is fixed here
as a setting. Written by a run and ignored by git, because they are geq's
output: `solution.txt` — the labels by `psiRZ` at the two points for each
resolution, the optimal `R_ref`, the Mach numbers on the converged state and
the fitted edge orders — and `psi_direct_{129,257,513}.dat` and
`psi_hagenow_257.dat`, total `ψ` as `R z psi` on the 65 × 65 nodes every
resolution shares. A test that compares against them needs the script run
first.

**geq's four grid corners carry exactly half the plasma's own flux** — a share
of 0.500000 at `( 0.05, ±1.5 )` and 0.500002 at `( 0.63, ±1.5 )`, at every
resolution, direct and Hagenow alike, while neighbouring nodes agree with MEQ to
3e-06. The cause is not located. Drop the corners from any comparison of `ψ`;
`RotatingMirrorResolution` asserts the share, so it goes red if geq changes, and
excludes them (MEASUREMENTS.md M-263).

```sh
cd /home/ian/projects/meq
PYTHONPATH=/home/ian/projects/geq/src \
tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_rotating.py \
    --out tools/geq-benchmark/reference/rotating
```

About eleven minutes, nine and a half of them the 513² direct solve (5 s at
129², 41 s at 257², the Hagenow control under 10 s), and 2 GB of memory.
Picard runs with `atol = 0` and `rtol = 2e-09`, because geq's default
`atol = 1e-06` is what actually stops it. `--cache DIR` keeps each solve and
reuses it.
