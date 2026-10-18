# A centrifugal mirror, free boundary, in MEQ

**The target**: MEQ solves a rotating magnetic mirror free boundary, and the
answer is checked against `../geq`. `CLAUDE_FLOW.md` is the campaign record for
the rotating source, `CLAUDE_FB.md` for free boundary, and this file is the plan
that joins them. `MEASUREMENTS.md` holds the tables under stable `M-nn` anchors.

**The source half is already done and measured.**
[M-167](MEASUREMENTS.md#m-167): `meq::RotatingSource` and geq's
`ProfilesCentrifugalMirror` are two independent implementations of RoPP (136),
and they agree to the accuracy the interchange tables carry, with no error floor
on either side. **So nothing below is about whether MEQ can evaluate a rotating
source.** It is about posing a mirror as a free-boundary problem in a solver
whose every free-boundary assumption was written for a tokamak.

## The standing assumption, and it is a modelling choice rather than a fact

> `Γ` is large enough, and the confinement good enough, that the plasma's
> toroidal current is negligible on and beyond `Γ`.

**It is NOT the assumption that there is no plasma outside `Γ`** — only that
whatever is out there does not move the magnetic field. That is the right
assumption for a well-confined centrifugal mirror and it is what the whole
exterior treatment rests on: MEQ truncates at a **semicircle** centred on the
axis and represents the exterior by Gegenbauer modes decaying as `ρ^(1−n)`,
which is exact where the exterior is current-free.

**IT IS SATISFIED BY PLACING `Γ` OUTSIDE THE MACHINE, AT EVERY MACH NUMBER** —
[M-169](MEASUREMENTS.md#m-169). The plasma ends inside `ρ = 1.63`, so a
semicircle beyond it has **exactly zero** plasma current outside it at `M = 2` as
at `M = 6`. There is no Mach threshold for admissibility. **The copper is a
different matter and decides `ρ_Γ`**: geq's coil packs span `ρ = 1.27` to
`2.18`, so `Γ` has to clear 2.18 and the coils are meshed inside it
([M-260](MEASUREMENTS.md#m-260)).

**AN EARLIER VERSION OF THIS SECTION SAID `M ≳ 6` AND THAT WAS WRONG.** It rested
on measuring the current beyond `ρ = 1.2` — a radius carried out of an
exploratory run into a conclusion, on a machine half a metre larger than it.
Two tables were built on it. The lesson is general enough to keep: **a threshold
measured at one arbitrary value of a free parameter is a statement about that
value.**

**WHAT LOW MACH ACTUALLY COSTS IS MODES, AND NOT MANY.** Relative L2 on `Γ` at
`ρ_Γ = 1.7`:

| | N=8 | N=16 | N=32 | N=64 |
|---|---|---|---|---|
| `M = 2` | 6.47e-03 | 1.12e-03 | **2.65e-05** | 1.99e-07 |
| `M = 6` | 6.00e-05 | 1.08e-05 | 4.22e-07 | 1.33e-08 |

`M = 2` at 32 modes is already far below the discretisation error of any mesh MEQ
would run, and at `ρ_Γ = 2.2` it reaches 1.78e-05 with **16**. Each mode is one
border column, and [M-98](MEASUREMENTS.md#m-98) made a column cheap by blocking
the trace solves.

**THE TRADE IS MESH AREA AGAINST MODE COUNT AND BOTH ARE CHEAP.** A half-disc to
`ρ_Γ = 1.7` is 4.54 m² against the machine's 1.89 m², about **2.4x**, and the
extra is vacuum where the field is smooth and the elements can be coarse — which
is what the adaptive loop does unprompted.

**SO THERE IS NOTHING TO EXTEND.** `ρ_Γ` and the mode count are both already
parameters of `meq::ExteriorDtN`. The one thing still worth having is the
DIAGNOSTIC: `M` is a property of the configuration and the current on `Γ` is a
property of the answer, so a run should report the latter rather than assume it.

**AND IF `Γ` EVER HAS TO SIT CLOSE** — a mesh budget, or a machine whose plasma
genuinely leaves the domain — the extension that fits MEQ is already sketched by
FB-7. `ψ = ψ_ext + ψ̃` holds for ANY known exterior current, not just coils; what
makes a coil special is that its current is fixed. Generalising it to the
iterate's own exterior plasma current means evaluating that current from the
exterior model, taking its field on `Γ` by the same Green's function this study
used, and adding it to BOTH halves of the datum — which FB-7 already warns must
happen together or not at all. It needs a quadrature over an unmeshed region, so
build it only if a measurement says `ρ_Γ` cannot move. At `ρ_Γ = 1.2` and
`M = 2` the irreducible error is **2.09e-02**, which is what that machinery would
be buying back.

## What transfers unchanged, and must not be rebuilt

| | |
|---|---|
| `meq::RotatingSource` and its two closures | [M-167](MEASUREMENTS.md#m-167). Nothing owed |
| the exterior DtN | `meq::ExteriorDtN` is a **semicircle centred on the axis**, and every mode vanishes on the axis identically, so the flat side needs no treatment. It knows nothing about what is inside it |
| conductors **inside** `Γ` | geq's packs span `ρ = 1.27` to `2.18` while the plasma reaches `ρ ≈ 1.51` at the throat, so no semicircle puts them outside `Γ`, and `ρ_Γ = 1.7` straddles both (`ExteriorCoilSet::clearance()` and `ConductorField::containment()` both refuse it). The packs are **meshed** rectangles — `[[coils]]`, `CoilAugmentedSource` — inside `Γ = 3.0`. FB-7's `ExteriorCoilSet` does not apply to this machine. [M-260](MEASUREMENTS.md#m-260) |
| the bordered Newton | the `( N + 2 )` elimination already carries `ψ_ax`, `ψ_bnd` and the exterior coefficients as one system |
| `ψ_bnd` as an unknown pinned at a POINT | `setBoundaryFluxPoint( R, z )`. Under NPC the row is exactly `−e_j` and the corner exactly 1, neither differenced. **This is the mechanism the mirror needs twice** |
| meshing | `tools/mesh/halfdisc.py`, a semicircle reaching the axis with coils meshed to, driven from `[mesh.generate]` |
| `η₆` | the exterior transmission indicator, for when the coefficients stall while the interior refines |
| the normalisation column | `NormalisedRotatingSource::normalisationDerivatives()` is analytic — `meq::NormalisedMHDSource`'s algebra under `g → H`, with the closure's derivative already paid inside `inner.dFdPsi()` — so the bordered Newton never differences across the support's edge. `theRotatingNormalisationDerivativesAreAnalytic` is the test |

## What is structurally wrong for a mirror

**A mirror has no magnetic axis.** `ψ` increases monotonically outward from
`R = 0`; there is no interior extremum. The plasma is an **annulus** in the
poloidal plane bounded by two flux surfaces, and geq labels it by two
**prescribed points** — `ψ_axis` at `( R_inner, 0 )` and `ψ_bndry` at
`( R_outer, 0 )`. Five consequences, in the order they bite.

### 1. Both `AxisConstraint` modes assume an interior extremum

`NodalMaximum` constrains `ψ_ax = max_j ψ_h`, and `LocatedAxis` constrains it at
a zero of `q_h`. A mirror has neither. What it needs is a third mode,
`PrescribedPoint`, which is **structurally `setBoundaryFluxPoint` reused** — the
same `−e_j` row, fixed at setup, nothing differenced and no root find. The
cheapest item on this list and the one everything else waits on.

`CriticalPointFinder::checkAxis()` then finds no interior extremum, which
`CLAUDE.md` already records as a **warning** rather than a refusal; it must stay
one. `checkAxisSource()` — is `μ₀ j_φ` bounded on `R = 0` — is untouched and
still has teeth, and for a mirror it is satisfied for a second reason: the
plasma does not reach the axis at all.

### 2. The support is a box in `( Ψ, z )`, and `insidePlasma()` is one-sided

Two bounds are missing, one radial and one axial.

**RADIALLY, `insidePlasma()` ADMITS THE WHOLE AXIS REGION.**

```cpp
return ( psi - edge )*span > 0.0;     // Source.hpp
```

With `ψ_ax` the inner label and `ψ_bnd` the outer, `span < 0` and the test reads
`ψ < ψ_outer` — correctly excluding the region radially outside the plasma and
**admitting everything inside the inner surface**, down to `R = 0`. In a tokamak
`Ψ ≤ 1` holds automatically because `ψ_ax` is the maximum; in a mirror it does
not. **The mirror needs `Ψ ∈ [ 0, 1 ]` where the tokamak needs `Ψ ≥ 0`** — one
extra bound. Without it the profiles clamp, so the source sits at its inner-edge
value all the way to the axis rather than vanishing: the same defect shape
`CLAUDE_FLOW.md` records for `ConfineToPlasma` reaching the rotating source, and
a current density filling a region with no plasma in it.

**AXIALLY, THERE IS NO BOUND AT ALL, AND A MIRROR NEEDS ONE.** Beyond the mirror
throat the flux tube flares into the expander, the centrifugal potential stops
confining, and the density is radically reduced — **best approximated as zero
within this model**. A purely flux-based mask cannot express that, because flux
is not a function of `z` along a tube: the same `Ψ` recurs in the expander and
gets enrolled as plasma. **That is not hypothetical** — it is exactly what
[M-168](MEASUREMENTS.md#m-168)'s `z_extent` table measured on geq, whose mask is
`0 ≤ ψ_n ≤ 1` with no axial cutoff, and where the end-plane current rose to
**five orders above the midplane's** as the domain grew.

So the mirror's support is

```
    0 <= Psi <= 1        and        |z| <= z_edge
```

with `z_edge` at or near the throat, **prescribed** rather than solved.

**IT IS A MODELLING APPROXIMATION AND MEQ SHOULD SAY SO**, the way it says so
about everything else: report the current density reached at the cutoff, so a
configuration where the approximation is doing violence is visible rather than
silent. A cutoff through a region still carrying current is a discontinuity in
`F`, and the profiles' own clamping will not soften it.

**AND IT MAKES [M-169](MEASUREMENTS.md#m-169)'s MODE COUNTS CONSERVATIVE.** That
study used geq's flux-only mask, which reaches the domain end at `|z| = 1.5`.
With the current cut at the throat the source is more contained, so `ρ_Γ` may
come in and `N` may come down from the 32 quoted there.

### 3. XP-1's element flood fill — not wanted on this path, dissolved by item 2

The axial cutoff makes the support a bounded box that needs no seed and no
connected component, so the fill is off here:
→ **[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md)**, *XP-1's element flood fill on the mirror path*

**Keep the fill's own guard in view even so.** `CLAUDE_FB.md` records XP-1's fill
being silently off as a real defect on the diverted case, found only because a
`2 × 2` cross separated it from the frozen edge. Turning it off for mirrors must
be a **stated** branch that a test asserts, not a flag that happens to be false —
which is the same failure this tree has already met twice.

### 4. The sense of the normalised flux is inverted

geq's `ψ_n = ( ψ − ψ_axis )/( ψ_bndry − ψ_axis )` is 0 at the inner surface;
MEQ's `Ψ = ( ψ − ψ_bnd )/( ψ_ax − ψ_bnd )` is 1 at the axis label. So
**`Ψ = 1 − ψ_n`** and every profile table has to be re-expressed, abscissa
reversed and derivative column negated. This is the trap
`tools/freegs4e-benchmark/source_check.py` documents for the MHD case — getting
it backwards converges to a different equilibrium and never complains — and it
is worse here because a mirror's two labels are not ordered the way a tokamak's
are.

### 5. The single-`R_ref` gauge does not survive a high-Mach multi-species mirror

**FOUND BY REGENERATING THE SOURCE BENCHMARK AT `M = 6`. IT IS NOT ON MR-4's
PATH AT TWO SPECIES; SEE THE CORRECTION BELOW BEFORE SIZING IT.** `meq::Species::density` is `n_s0( ψ )` — the physical density on ONE
curve `R = R_ref`, for every flux surface. A surface whose own midplane radius is
`R_mid` therefore carries `exp( m_s ω²( R_ref² − R_mid² )/2T_s )`, an exponent
**linear in the species mass**. At `M = 6` that is +31.2 to −51.8 for deuterium
and **+187.2 to −310.5 for carbon**, so carbon's `n_s0` spans about **1e216** —
against a double's `exp( 709.8 )` ceiling — and the three-species case cannot be
posed at all. → [M-167](MEASUREMENTS.md#m-167) §7.

**`R_ref` is already at its optimum**, `0.2549` against the range-minimising
`0.257`, so this is not a tuning failure. **A MIRROR IS A HIGH-MACH
MULTI-SPECIES DEVICE BY DEFINITION**, so this is not a corner either: `M ≥ 4` is
the operating point and impurities are the reason to run three species at all.

### The fix: a fixed-`z` gauge for mirrors, and the fixed-`R` one stays for tokamaks

**BUILT AND MEASURED** → [M-278](MEASUREMENTS.md#m-278). `[source] Gauge =
"midplane"`, `meq::RotatingSource::Gauge::Midplane`, refused off a boxed
normalised source and beside `ReferenceRadius`. **geq's mirror solved in the
midplane gauge agrees with geq as MR-4d asks**: falling at 1.69 then 1.83 against
geq's 1.63, and at 513² 0.39 of geq's own 257-against-513 difference, where the
fixed-radius gauge sits at 0.44. **The same plasma in the two gauges is one
equilibrium** to 3e-10 at level 2, 25× under MEQ's own error there and falling
faster than it. **Three species at `M = 6.5` solve** and self-converge, with
every `n_s0` the physical midplane density to 5e-14 in the log, where the same
plasma's fixed-radius carbon table would span 1e117. `RotatingMirrorMidplane` is
the acceptance. Three findings the paragraphs below did not anticipate are
recorded beside the paragraph each corrects.

**THE DECISION IS TAKEN.** A tokamak keeps
`φ₀( R_ref, ψ ) = 0` — a constant radius, exactly as today. A mirror uses

```
φ₀( R, z = 0 ) = 0
```

— pinned, on every flux surface, at that surface's own **midplane crossing**
`R₀( ψ )`. `n_s0( ψ )` is then the physical density where the surface crosses
`z = 0`.

**WHY THE SPLIT IS FORCED RATHER THAN A PREFERENCE.** `z = 0` crosses a tokamak
surface **twice** — inboard and outboard of the magnetic axis — so the condition
does not define a gauge there at all without a further choice. On a mirror `ψ` is
monotone in `R` along the midplane, so `R₀( ψ )` is single-valued and the
condition is exactly one per surface, which is what a gauge needs. **That
monotonicity is a precondition MEQ can check** rather than assume, and it should
be checked where `checkAxisSource()` is checked.

**THE EXPONENT ABOVE IS AN UPPER BOUND, ABOUT TWICE THE REAL ONE.**
`m_s ω²( R_ref² − R_mid² )/2T_s` is the `φ₀ = 0` limit. With the electrons in,
quasineutrality makes the deuteron's exponent the SHARED one,
`m_D ω²ΔR²/( 2( T_e + T_i ) )`, and since carbon and deuterium have the same
`Z/m`, carbon's is exactly six times deuterium's whatever `φ₀` is. The shared
form is the one M-167 measured: at `M = 6` it gives a range of 51.6 in the
exponent, `e^51.6 = 2.6e22`, against the **2.538e22** spread
[M-167](MEASUREMENTS.md#m-167) reports for `n_s0`; §7's +31.2/−51.8 are the
upper bound. On geq's machine at the optimal `R_ref`, `n_s0` overflows a double
at about `M = 32` for deuterium and **`M = 13` for deuterium and carbon**, with
geq's temperature profiles (`M` as M-167 defines it, on `T_e`). The practical
limit before that is the spline: a Hermite table holds `n_s0` to `( E′h )^4/384`
relative, which is about 500 knots for 1e-08 on deuterium at `M = 4` and 6500
for carbon at `M = 6`. **So the fixed-`R` gauge carries MR-4's two-species
comparison, and the midplane gauge below is the route for three species at high
Mach rather than a precondition for any mirror.**

**AND IT REMOVES THE DEFECT AT THE ROOT, MEASURED.** The gauge-transfer exponent
is `m_s ω²( R_ref² − R₀( ψ )² )/2T_s`, and the fixed-`z` gauge makes it
**identically zero on every surface by construction**. On the `M = 6` state:

| | fixed-`R`, max `\|A\|` → range of `n_s0` | fixed-`z` |
|---|---|---|
| deuterium | 64.9 → 2.13e+56 | **0 → 1.0** |
| carbon | 389.1 → **inf** | **0 → 1.0** |

So `n_s0`'s dynamic range becomes the *physical* density range, and the
three-species case at `M = 6` — which cannot be posed at all today — becomes
ordinary.

**WHAT IT COSTS, AND IT IS THE PROPERTY THE HEADER ADVERTISES.** `R_ref` is a
constant; `R₀( ψ )` is **a functional of the solution**. So `F` stops being a
pointwise function of `( R, z, ψ )`, which is what
`RotatingSource.hpp` gives as the reason a non-flux-function density needs no
interface change anywhere in the solver. That reason does not survive this.

**IT IS STILL MUCH CHEAPER THAN THE AVERAGE MEQ REJECTED.** RoPP's own gauge is
`⟨φ₀⟩_ψ = 0`, a flux-surface average of the unknown on every surface. `R₀( ψ )`
is a **1D root find along one line**, `ψ( R, 0 ) = ψ`, on a monotone function —
not an integral over a surface, and cacheable as a single tabulated `R₀( ψ )` per
Newton iteration.

**THE JACOBIAN TERM IS ANALYTIC AND THE MIXED METHOD PAYS FOR IT AGAIN.**
`∂A_s/∂ψ` gains `−m_s ω² R₀ R₀′( ψ )/T_s`, and by the inverse function theorem
`R₀′( ψ ) = 1/( ∂ψ/∂R )|_{z=0}`, which is `1/( R₀ q_R( R₀, 0 ) )` since
`∇̄ψ = R q`. **`q` is a SOLVED variable at the same order as `ψ`**, so this is
read off the answer rather than differenced — the same argument that lets the
`.nc` band continue at `q`'s order.

**AND `ψ` ALONG THE MIDPLANE IS `q` INTEGRATED, NOT `ψ_h` READ — MEASURED.** `ψ_h`
is discontinuous between elements, and on the mirror's level-0 mesh at `k = 3` it
jumps by 1.3e-03 across a face — the plasma's own flux is 8.7e-04 — with a slope
inside an element 30% off `R q_R`. A crossing table read off it has knot values
its own slopes disagree with by 100%, and a source gauged on it moves by that
much at every refresh: Newton stalled near 3e-04. With `ψ( R, 0 )` the integral
of `R q_R` from the midplane's first point, panel by panel over the elements, the
mismatch is 2e-04 and the same solve converges.

**FREEZE `R₀( ψ )` WITHIN A NEWTON ITERATION — BUILT, MEASURED LINEAR, AND
RETIRED.** The tree's own precedent is `setPlasmaSupportFrozen` and XP-2's
frozen plasma edge: a reference that moves inside a Jacobian makes a difference
straddle a kink instead of measuring a derivative, so the table was refreshed
once per iteration, beside the support freeze. **It made Newton linear**: the
field's pull on `R₀` is the one term that Jacobian does not carry, so geq's
problem took 6 steps where the fixed radius takes 3, contracting about 50× a
step. **That the left-out term is what capped the pressure — three species at
`M = 6` stalling near 1.5e-05 at a midplane peak of 2e19 m⁻³ and 8e-05 at
5e19, where they converge at 5e18 — was believed and is FALSIFIED**: the cap is
the mesh's. On the vacuum's fixed curve below, with an exact Jacobian, the same
three species stall on level 0 between 1.25e19 and 1.5e19 even by continuation
in pressure, and converge at 2e19 on levels 1 and 2; differencing the border
columns, or taking `F` at the labels down to 1e-13 of its peak, leaves every
stall where it was → [M-281](MEASUREMENTS.md#m-281). The held table
with nothing moving stalls there too → [M-285](MEASUREMENTS.md#m-285).
**And the mesh's cap is the STABILISATION's**: at the default constant `τ = 1`
one element's local HDG block goes singular when the reaction `(1/R) ∂F/∂ψ`
reaches its lowest eigenvalue `(k + 1) τ |∂K|/|K|`, so the discrete problem
folds inside that element — at 1.367e19 on level 0, the step between
neighbours living in one element and growing as `(n_c − n)^(−1/2)`. The
ceiling scales with `τ`, and **a box support therefore defaults to the
mirror's stabilisation**, `τ = 3/max(R, 0.1 m)` (`meq::MirrorStabilisation`;
→ [M-293](MEASUREMENTS.md#m-293) chose the two numbers): 2e19
converges in three steps on level 0 in both midplane gauges. **What is left
differs by gauge.** On the vacuum's curve it is the posed problem's own fold,
smooth, global and mesh converged near **4.16e19**, so 5e19 there has no
equilibrium on any mesh; under the midplane gauge the only fold is `τ`'s,
linear in `τ₀` and near 2.15e20 on level 0 at `τ₀ = 3`, and 5e19 converges on
level 0 in four steps, which `RotatingMirrorVacuumGauge` asserts
→ [M-291](MEASUREMENTS.md#m-291),
[M-294](MEASUREMENTS.md#m-294).

**THE CROSSINGS ARE BORDER UNKNOWNS — BUILT, AND WHAT SHIPS**
→ [M-283](MEASUREMENTS.md#m-283). `GradShafranovSolver` makes `R₀`
at `K` knots in `Ψ` unknowns beside the labels and the exterior modes — 32 by
default, `setMidplaneGaugeKnots()` in the library and no key — the table the
`C²` spline of `R₀²/2` through them, one row per knot,
`ψ( R₀ₖ, 0 ) = ψ_bnd + Ψₖ( ψ_ax − ψ_bnd )` with `ψ` along the midplane `q`
integrated. The knots start on the starting iterate's crossings: the vacuum's
from cold, the warm start's on an adaptive cycle. geq's problem converges in 3
steps where the frozen refresh took 6, three species at 5e18 in 3 where it took
7, and at 2e19 on level 1 in 3 where it took 14, landing on the frozen refresh's
answer to 1.3e-08 in `ψ*` at level 0 — the interpolant's error, a thousandth
of that level's discretisation error. **It does not lift the stall**: at 2e19
and 5e19 on level 0 it stalls at the frozen refresh's merit from every start
tried, and by continuation in the peak both fail together between 1.63e19 and
1.67e19 with the border's condition flat. **M-278 §2's six steps with the
vacuum's table held are not this tree's**: held, it converges in the fixed
radius's 3, so a fixed table's Jacobian is the whole one. **Nothing only the
frozen refresh could do was found** — it converged nowhere the knots do not, in
2 to 5 times the steps — **so it is retired**, and with it the per-iteration
`NormalisedSource::setMidplaneField()` hook. A prescribed plasma current beside
the knots is refused (its corner against them is not assembled); the driver
refuses a current beside the box, so no run loses it.

**THE START THE KNOTS WERE GUESSED TO NEED — MEASURED, AND NOT TAKEN.** The
guess: solve the same densities first on the vacuum's crossing held — the
vacuum-midplane gauge, quadratic and cheap — so that the bordered Newton then
needs only a very few of its dearer steps. **It saves one bordered step of three
at most, none at 2e19, and costs more than it saves**: the whole held-then-
bordered route takes 1.35× to 1.66× the wall clock of the bordered Newton from
the vacuum, and still 1.06× to 1.36× with the held solve stopped at the derived
hand-over, the merit the bordered problem reads at the held answer. The
bordered Newton from the vacuum is already quadratic in three steps; its first
step does what the whole held solve does, and the held answer's distance from
the true one grows with the pressure — 1.9%, 4.3% and 22% of the plasma's flux
for geq's two species, three at 5e18 and three at 2e19 — exactly where a better
start would be wanted. Where every start stalls, the held solve stalls first.
→ [M-284](MEASUREMENTS.md#m-284)

**AND THE OUTER LOOP — Picard on the curve, whole held-table solves with the
table each time the last answer's crossings, no knots — IS THE FROZEN REFRESH
WITH ITS INNER SOLVES CONVERGED, AND NO BETTER.** It is linear at the frozen
refresh's own contraction, 0.016 for geq's two species, 0.041 for three at
5e18 and 0.18 at 2e19 on level 1, so the left-out term grows about linearly in
the pressure; it costs 3 to 6 times the bordered gauge, and its first inner
solve at 2e19 on level 0, a fixed residual, stalls where the bordered and the
frozen ones do. **So the hypothesis that the 2e19 stall was the line search
fighting a moving residual is FALSIFIED**: nothing moves in that solve. It is a
test-side harness → [M-285](MEASUREMENTS.md#m-285).

**THE INTERFACE SHAPE** is a gauge policy on the source rather than a second
class — `enum class Gauge { FixedRadius, Midplane, FixedCurve }`, `FixedRadius` carrying
today's `referenceRadius` and remaining the default so every tokamak caller is
untouched, `Midplane` carrying the interpolant through the solver's knots,
replaced at every residual, `FixedCurve` the curve it installs once.

**AND IT MAKES THE BENCHMARK EASIER, WHICH IS A CHECK ON THE DESIGN RATHER THAN
A BONUS — FALSIFIED, AND THE REASON IS MR-4a's HOLD.** MR-4a holds geq's `N_s`
as prescribed flux functions, which are geq's ZERO-POTENTIAL densities, not its
midplane ones; the midplane tables are geq's quasineutrality solved at each
surface's crossing on geq's 513² state, `max |ln( n_mid/N_s )|` = 2.37 against
1.17 for the fixed-radius transfer, and since they need geq's answer they are an
output, `nmid_*.dat`, beside it. What follows was written expecting otherwise: geq's `get_flux_density_qn` already references each surface to its own
midplane radius, so under the fixed-`z` gauge MEQ's `n_s0` and geq's `N_s` are
the same quantity up to geq's single global offset at `ψ_n = 0.5`. The transfer
in `tools/geq-benchmark/export_mirror.py` becomes near-identity, and the
reference should be regenerated under the new gauge once this lands —
`theGaugeTransferIsLoadBearing` will then be asserting something much weaker and
should be re-pointed at the fixed-`R` gauge, which is where the trap still lives.
**Re-pointed by name, and the refusal asserted under the midplane gauge too**,
since geq's raw `N_s` are no closer to the midplane densities than to `R_ref`'s.

**AND "FLAT AT A LABEL" MEANS SOMETHING DIFFERENT UNDER THIS GAUGE.** The box's
precondition is `F` vanishing at both labels. A midplane-gauged `n_s0` flat at a
label still carries `−n_s m_s ω² R₀ R₀′/T_s` in `∂n_s/∂ψ|_R` wherever `ω` does
not vanish there, and `R₀′` is the solution's. So a profile family written in
this gauge has to let the densities or the rotation fall away at the labels —
the three-species case takes its densities to 1e-6 of their peak there — or be
transferred from a fixed-radius statement of the same plasma, as geq's is.
**On a fixed curve it is checked before the solve and refused** (below); under
this gauge it cannot be, since `R₀′` at a label is the answer's and the
vacuum's is not an estimate to refuse on — measured, geq's midplane tables read
8.3e-08 of their peak at a label on their own crossings and 2.8e-05 on the
vacuum's — so the driver reports it on the answer instead.

**THE `M = 6` THREE-SPECIES CASE, MEASURED RATHER THAN BOUNDED.** On the solved
plasma — `T` to a tenth at the labels, `ω` drifting 1.3×, carbon at 1% — the best
single `R_ref` is 0.31 m, where carbon's fixed-radius `n_s0` spans **1e117**, the
electrons' 1e68 and the deuterons' 1e23. That is inside a double; it is the
Hermite table at `( E′h )⁴/384` that cannot carry it, at thousands of knots.
Under the midplane gauge each
species spans the 1e6 it was given. `examples/mirror-three-species.toml` runs
it from a file.

### A fixed curve is a gauge: the vacuum's crossing — BUILT AND MEASURED

→ [M-281](MEASUREMENTS.md#m-281), [M-282](MEASUREMENTS.md#m-282).
**A GAUGE PIN DOES NOT HAVE TO LIE ON THE SURFACE**: quasineutrality fixes
`φ₀( R, ψ )` pointwise for any `( R, ψ )`, which is why `R_ref` is a gauge on a
mirror. So any fixed, known curve `R_g( Ψ )` is one, with `F` pointwise in
`( R, z, ψ )` and an exact Jacobian. `RotatingSource::Gauge::FixedCurve` is
that gauge, and `[source] Gauge = "vacuum-midplane"` ships it on
`R_g = R₀vac( Ψ )`, the vacuum's midplane crossing: the solver tabulates it
once, before the first residual, from `ψ_c` and a zero remainder through the
same `R q_R` integration as the midplane table, normalised by `ψ_c` at the two
label points, and never refreshes it — not per iteration, not per adaptive
cycle. The exponent from a `Density` to the solved crossing is then
`m_s ω²( R₀vac² − R₀² )/2T_s`, the plasma's own displacement of its surfaces.

**WHAT IT NEEDS.** The vacuum is `ψ_c`, so it is refused without
`[conductors]` subtracting and beside data on `Γ`; a curve from anything else —
the guess, a previous answer — would make what each `Density` means depend on
it. **No key for a user's own curve**, for that reason →
[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md), *A user-supplied gauge curve*.

**THE BOX'S PRECONDITION IS CHECKED AT SETUP**, `RotatingSource::labelCheck()`:
`|F|/R` at each label over every radius that label's surface reaches, against
its peak along the curve, refused above `meq::labelTolerance` = 1e-6 as a
configuration error naming the label. That is MR-4a's amplitude, and what it
guards is the edge's order: measured, `F` at a label at 1e-3 of its peak still
leaves Newton three quadratic steps.

**WHAT IT MEASURED.** geq's plasma on the vacuum's curve is the fixed-radius
plasma — three steps, the same history to every printed digit, `ψ*` together
to 2e-13 — and the differenced border columns give the same history. **M-278
§2's six steps on a frozen vacuum table are not reproduced**: the midplane gauge
keeping its first table, and the fixed curve, both took three quadratic steps
on the same tables, so nothing is missing from the Jacobian with the table
fixed, and the refresh was the whole of the old midplane gauge's linear rate. Three
species at `M = 6`: three quadratic steps at 5e18 on every level and at 2e19 on
levels 1 and 2; carbon's exponent to the solved crossing is 0.74 at 5e18 and 3.3
at 2e19 against the best fixed radius's 134. The pressure cap on level 0 is
constant `τ`'s element-local fold, and 5e19 is past the plasma's own
(above, → [M-291](MEASUREMENTS.md#m-291)).

**OUTPUT.** Every rotating mirror's `.nc`, in any gauge, carries `R₀( Ψ )` and
`n_s( R₀( Ψ ), Ψ )` with slopes on the solved crossing, as `DensityFile`
tables; under the midplane gauge they are its own tables back at its knots and,
between them, the tables moved by the knots' interpolant's gap from the crossing
— 7e-05 relative for carbon, 1.5e-07 m of gap — and off a vacuum-gauge answer
they pose the same plasma in the midplane gauge to 1.2e-06 of `ψ*`, the same
gap's (→ [M-284](MEASUREMENTS.md#m-284) §4).

**TRUE-MIDPLANE DENSITIES WITH A QUADRATIC NEWTON — DONE, AND NOT BY THE
INTENDED ROUTE.** `Gauge = "midplane"` is quadratic itself with `R₀` at its
knots as border unknowns, from the vacuum. The route intended here — solve on
the vacuum's curve to a start near the answer, convert, then a few bordered
steps — was built and measured, and loses: see *The start the knots were
guessed to need*, above → [M-284](MEASUREMENTS.md#m-284).

### What does NOT bite, which is worth saying

**The plasma-current border is not needed**, because both labels are pinned and
the support is not free:
→ **[CLAUDE_DECLINED.md](CLAUDE_DECLINED.md)**, *A plasma-current border for the mirror*

**`g = 0` is fine.** `F = μ₀ R² ∂p/∂ψ ≠ 0` inside the annulus, so the trivial
branch that `examples/rotating-rectangle.toml` warns about is not reached. The
vacuum coil field is the initial guess.

**BUT `F` MUST VANISH AT THE INNER LABEL, NOT ONLY AT THE OUTER ONE.** Under the
box, `Ψ = 1` is an edge of the support, and it moves with `ψ_ax`.
`setPlasmaSupport()`'s precondition — `F → 0` at the edge, or the Newton
Jacobian is missing the edge term — therefore applies at both labels, so a
mirror's profiles need `p′` to vanish at `Ψ = 1` as well as at `Ψ = 0`. That
decides how MR-4 writes them, and how geq's profiles are transferred.

## The stages

Each ends at a measured acceptance, not at "it runs".

### MR-0 — pick `ρ_Γ` and `N`, and keep the diagnostic

**THIS STAGE WAS THE CAMPAIGN'S FALSIFYING EXPERIMENT AND
[M-169](MEASUREMENTS.md#m-169) HAS RUN IT.** The exterior is admissible at every
Mach number with `Γ` outside the machine, so what is left is small and is
engineering rather than physics:

* **A choice**: `ρ_Γ` must clear the **copper** (`ρ = 2.18`), not the plasma.
  `ρ_Γ = 3.0` with `N = 32` is measured and on the plateau
  ([M-260](MEASUREMENTS.md#m-260)). [M-169](MEASUREMENTS.md#m-169)'s counts hold
  for the plasma's own current; with the coils inside `Γ`, the copper's
  `( 2.18/ρ_Γ )^n` tail sets `N`. Both are existing parameters.
* **A MEQ-side diagnostic**: the fraction of `∫|F|/R` lying outside `Γ`,
  reported on the converged answer. It is the mirror's analogue of
  `checkAxisSource()` — one-sided, refusing or warning on a positive detection —
  and it exists because `M` is chosen while the current on `Γ` is solved for.

**Acceptance**: the diagnostic reads below `1e-3` on a shipped mirror
configuration, and a mode-count sweep on MEQ's own solve reproduces
[M-169](MEASUREMENTS.md#m-169)'s trend. **The risk this stage existed to retire
is retired**, so it no longer gates the others; MR-1 and MR-2 can start now.

### MR-1 — the plasma support becomes a box in `( Ψ, z )` — DONE

**DONE**, → [M-258](MEASUREMENTS.md#m-258). `meq::PlasmaSupportShape::Box`, set
by `NormalisedSource::setPlasmaSupportShape( Box, z_edge )`, makes the support
`0 < Ψ ≤ 1` and `|z| ≤ z_edge`; the upper bound is written in `Ψ`, so a mirror's
negative span needs nothing. Asking for a box switches the support on; an edge
that is not finite and positive is refused. `OneSided` stays the default and is
the tokamak's `Ψ > 0` unchanged — it must not gain `Ψ ≤ 1`, because `ψ_h` can
overshoot the located `ψ_ax`. `insidePlasma( R, z, ψ )` takes `f()`'s
arguments and has no flux-only form, because a flux-only overload would be a
gate that silently skips the axial bound. The box reaches `f`, `dF/dψ`, both
normalisation columns, the coil wrapper by forwarding, and the device
`SourceView`. XP-1's fill is off under the box by a switch in
`plasmaComponentWanted()`, asserted in `MirrorSupport` against a one-sided
control, and a mask left from an earlier fill is cleared rather than kept.
`GradShafranovSolver::checkAxialCutoff()` reports the peak `μ₀j_φ`, its value on
`|z| = z_edge` and their ratio; it is a report, since no threshold is measured.


`insidePlasma()` gains a radial upper bound (`Ψ ≤ 1`) and an axial one
(`|z| ≤ z_edge`), and the element-level fill is switched off on this path as a
stated branch rather than an unset flag. **`z_edge` needs `meq::Source::f`'s `z`,
which it already takes** — the same fact that let a non-flux-function density in
at all.

**Acceptance**: a unit test that a mirror-shaped source vanishes on all three
sides — inside the inner surface, outside the outer, and beyond the throat —
mutation-checked the way `aConfinedRotatingSourceVanishesOutsideThePlasma` was;
plus a diagnostic reporting the current density reached at the cutoff, so the
approximation is visible. Cheap, and it needs no mirror solve.

### MR-2 — `AxisConstraint::PrescribedPoint` — DONE

**DONE**, → [M-259](MEASUREMENTS.md#m-259). `setAxisFluxPoint( R, z )` selects
`AxisConstraint::PrescribedPoint`. The row is the containing element's potential
shape functions, located once per solve, corner 1, live from iteration 0.
Pinned at the located axis it reproduces `LocatedAxis` to 11 digits in `ψ_ax`
and `ψ_bnd` on the stiff high-β case and FB-3's two-border case, and is
quadratic where `LocatedAxis` is linear. The setter refuses a non-finite point,
`R ≤ 0`, a point outside the mesh, and choosing the mode without a point.
`axisWasLocated()` is false under it and `axisR()`/`axisZ()` report the pin.


`ψ_ax` pinned at a prescribed point, reusing `ψ_bnd`'s border row.
**Acceptance**: on an existing *tokamak* fixture, prescribing the point at the
known axis reproduces `LocatedAxis`'s answer to round-off. **It needs no mirror
at all**, which is what makes it a safe stage to land early.

### MR-3 — the vacuum mirror, through MEQ — DONE

**DONE**, `VacuumMirror`, → [M-260](MEASUREMENTS.md#m-260). The coils are meshed
inside `Γ = 3.0`, not outside it (see the table above).
* `ψ` converges on `coilPsi()` at `k+1` for `k = 1…3`, over the machine region
  and over `D_h`.
* `N` stops mattering over `D_h` at 24, about where the copper's trace tail,
  `0.725^N`, drops below the discretisation error. The machine region is on its
  plateau by `N = 8`.
* Against geq the agreement is 7.29e-04, and that is geq's coil model: ten
  filaments, the midpoint rule for the winding, 7.28e-04 from `coilPsi()` with
  no solve. MEQ's own error at those points is 1.6e-06, and handed geq's
  filaments MEQ agrees with geq to 3.9e-10.

No plasma. Mirror coils, half-disc mesh, DtN live. **Acceptance**: MEQ's `ψ` against the coil field computed directly from
`meq::coilPsi()`, and against geq's own vacuum field. This validates the
geometry, the mesh, the coil model and the exterior **together and without the
plasma**, which is the rung that does not exist today between the analytic
fixtures and a full machine solve. It also sizes the mode count `N` that an
elongated machine needs on `Γ` — the number MR-0 bounds and this one measures.

### MR-4 — the free-boundary rotating mirror, against geq — DONE

**DONE: MEQ AGREES WITH geq TO geq'S OWN RESOLUTION, AND THE DIFFERENCE FALLS AS
geq REFINES, AT geq'S RATE** → [M-261](MEASUREMENTS.md#m-261),
[M-262](MEASUREMENTS.md#m-262), [M-263](MEASUREMENTS.md#m-263).
`RotatingMirror` (MR-4b, seconds) and `RotatingMirrorResolution` (MR-4c and
MR-4d, about four minutes) are the acceptance. From the vacuum guess the bordered
Newton converges in three or four steps with a negative span, quadratically; `ψ*`
converges at its design rate; and against geq's direct solves at 129², 257² and
513² the difference falls at 1.64 to 1.83 against geq's own 1.60 to 1.63,
ending at 5.4e-10 rms, 0.44 of geq's own last refinement. Shipped geq reads its
8.06e-06 bias through the same comparison. Departures from what follows: the
MEQ-native profiles are gone, the problem being geq's; the packs are not meshed
on the plasma mesh, since under subtraction they are no source and a level-2
vertex landed on a filament; `ψ_h` at `k ≤ 2` is printed and not asserted,
being not yet asymptotic; and geq's four grid corners are excluded, for the
reason M-263 records.

Everything at once: rotating source, both labels prescribed at points, the box
support, the coils, the DtN, the bordered Newton — compared with geq's converged
rotating mirror **in the mode geq solves in**.

**NOTHING IN THE SOLVER IS MISSING.** Every piece is library API, so MR-4 is a
ctest in the shape of `VacuumMirror` and not a driver change. What it needs is a
reference that does not exist yet, profiles that meet the support's
precondition, and three small corrections to diagnostics.

#### What geq solves, and therefore what MEQ must solve

| | geq | MEQ, for MR-4 |
|---|---|---|
| coils | 20 filaments, currents FIXED (`control=False`) | the same 20 filaments **by subtraction**, which [M-260](MEASUREMENTS.md#m-260) puts at 3.9e-10 of geq with no plasma — so the coil model leaves the comparison entirely |
| labels | `ψ_axis` read at `( 0.10, 0 )` and `ψ_bndry` at `( 0.35, 0 )` EVERY iteration, so both float | `setAxisFluxPoint( 0.10, 0 )` and `setBoundaryFluxPoint( 0.35, 0 )` |
| amplitude | the profiles' own SI units; no current constraint, no scale | the same: `NormalisedRotatingSource` in SI, **no `setPlasmaCurrent()`** |
| support | `0 ≤ ψ_n ≤ 1`, and the box end `\|z\| = 1.5` as the de facto cutoff (boundary-row current dropped) | `PlasmaSupportShape::Box`, `z_edge = 1.5` |
| profiles | in `ψ_n`, 0 at the inner label | in `Ψ = 1 − ψ_n`: abscissa reversed, slopes negated |
| route | Picard, `solve_equilibrium( ..., bndry )`; the Newton–Krylov and hybrid paths FREEZE the labels and `ω` as values, which is a different mode and is not used | bordered Newton |

**geq re-derives two things every Picard step, and the reference must freeze
them**: `ω` from the Mach number (so pass `bndry["omega_profile"]`, a rigid
`ω( ψ_n )`), and the deuteron's flux-function density `N_D`, which it rebuilds
from the electrons' midplane density by quasineutrality — so geq prescribes the
MIDPLANE density, not `N_s`. Frozen at convergence, geq's state is still a fixed
point of the frozen problem; only its iteration changes.

#### Three facts about geq that set the reference

1. **SHIPPED geq IS BIASED BY ABOUT 1% OF THE PLASMA FLUX, AND ITS OWN
   SELF-CONVERGENCE CANNOT SEE IT.** `freegs4e/boundary.py:213`'s von Hagenow
   offset is a fixed `eps = 1e-2` m, the defect `tools/freegs4e-benchmark`'s
   WORKAROUND 5 exists for. At 257², shipped minus direct `freeBoundary`:
   `Δψ_axis` 1.15e-05 and 3.43e-05 against a plasma flux of 1.24e-03 and
   2.95e-03 at geq-Mach 2 and 4. Shipped converges cleanly — 129 against 257
   differs by 3e-07 to 3e-06 — to the wrong limit; Hagenow with `eps = 0.2h`
   converges at first order to the direct answer. **The reference is built
   with `boundary=freeBoundary`** (6–8 s at 129², 42–61 s at 257², about 500 s
   at 513²).
2. **geq'S SOURCE STEPS AT BOTH LABELS ON EVERY SHIPPED PROFILE SET.** Its
   profiles are Gaussians peaked at `ψ_n = 0.5`, with `d ln T/dψ_n ≈ ±10.6` and
   `d ln N_e/dψ_n ≈ ±12.5` at the labels; `J_φ` just inside them reads
   −1.3e-04 to −2.6e-02 of peak at the inner label and up to **−1.76e-01** at the
   outer. That violates the precondition under *What does NOT bite*, which MEQ
   cannot relax: `CLAUDE_FB.md`, *At `j = 0` the question does not arise*, records
   a source with `p′ ≠ 0` at a support read off the solution failing to converge
   at any `k`, any mesh, from the exact solution and under `PicardThenNewton`. **The reference uses new profiles.**
3. **geq'S `φ₀` RESIDUAL IS NOT ITERATION ERROR.** Its global quasineutrality
   solve uses a pointwise `N_D` with linearly interpolated `R₀`; its `J_φ` uses
   a cubic spline of `N_D` through the grid nodes. The two disagree by up to
   3.3e-02 in `ln N_D`, concentrated at the inner edge, and `Σ Z n` misses by
   2.6e-04 to 3.9e-03 rms. A pointwise Newton on the spline densities after the
   global solve closes it to 1.5e-15 at no Picard cost, and moves the solved
   field by 2.5e-07 to 3.7e-06. **The reference is tightened this way.**

**geq'S FLOOR, ON THE DIRECT PATH**: 257² against 513² differs by 5.9e-08
(geq-Mach 2) and 1.7e-07 (geq-Mach 4) in max `ψ`, `ψ_axis` by 4e-09 and 4e-08 —
about 1e-07 Wb/rad, or 2e-06 of the label span. Convergence is irregular, since
the mask is a staircase on geq's grid, and its plasma current converges at first
order. **That is the resolution any MR-4 agreement can be asserted to.**

**geq's Mach is not M-167's.** geq defines `M = ωR/√( 2T_i/m_i )`, which is
M-167's `M` (on `T_e`) divided by about 1.55: geq-Mach 2 is M-167's 3.1, and 4 is
6.2. The plan quotes geq's.

#### The profiles: `F` must vanish at both labels, to an order

`F = μ₀R² ∂p/∂ψ|_R` vanishes at a label when `N_s′`, `T_s′` and `ω′` all vanish
there (with `gg′ = 0`): then quasineutrality has no `ψ`-dependence at that label,
`φ₀′ = 0`, and `∂p/∂ψ = 0`. That is a statement about the PHYSICAL state, so it
holds in geq's gauge and in MEQ's alike, and the gauge transfer preserves it.

**WHAT AN EDGE COSTS IS ITS ORDER `j`**, `p′ ~ d^j` at distance `d` inside it:
[M-09](MEASUREMENTS.md#m-09) caps `ψ*` at `min( k+2, j+2.5 )`. So the profiles
are

```
    f( ψ_n ) = f_edge + ( f_peak − f_edge )·C ψ_n³ ( 1 − ψ_n )⁴
```

for `T_e`, `T_D`, `N_e`, `N_D` and `ω`, `C` normalising the peak at `ψ_n = 3/7`
to 1: `j = 2` at the inner label and `j = 3` at the outer, and polynomial, so
MEQ's Hermite tables carry them with explicit derivative columns. A quartic
`16 ψ_n²( 1 − ψ_n )²` is `j = 1`, which caps the rates from `k = 2`, and is
symmetric about 0.5, so a reader confusing `Ψ` with `ψ_n` reads the same numbers;
**the asymmetry is what lets the comparison see item 4's inversion**, and it
does — read in `ψ_n`, the label powers swap and the difference sits at 8.6e-04.

**geq holds `N_D` as a flux function by a monkeypatch**: `AllSpecies.set_flux_density`
is kept from rebuilding deuterium's spline, so both `N_s` are the prescribed
profiles and `φ₀` alone carries quasineutrality.

#### Stages

**MR-4a — the reference.** `tools/geq-benchmark/export_rotating.py`: two
species (e, D), the profiles above, rigid `ω`, `N_D` held, direct
`freeBoundary`, `φ₀` tightened, Picard with `atol = 0` and `rtol = 2e-09` (geq's
default `atol = 1e-06` is what actually stops it, and its returned `psi_bndry`
then lags `psiRZ` at the point by 2.3e-07 — **export the labels by `psiRZ` at
the points**). At geq-Mach 2, and at 4 if 2 is clean. Write the tables in MEQ's
`Ψ` and SI with the gauge transferred pointwise as `export_mirror.py` does, `ψ`
on geq's grid at 129², 257² and 513², and the labels — the profile tables
tracked, geq's answer regenerated and not. **Acceptance**: `J_φ` falling at
`d^2` and `d^3` inside the two labels on the converged state, where the shipped
profiles step; geq's direct self-convergence table re-taken; the
shipped-Hagenow bias re-taken as the control that shows the comparison can see
a 1e-05 error. → [M-261](MEASUREMENTS.md#m-261).

**MR-4b — MEQ converges on it.** `tests/convergence/RotatingMirror.cpp`, with
`VacuumMirror`'s machine, mesh and domain extracted into a shared header
(`MirrorMachine.hpp`, as `DivertedMachine.hpp` is for the diverted case) rather
than copied. The sequence the setters impose:

1. A vacuum solve on its own solver (one solver holds one source), as
   `VacuumMirror` does, for the guess.
2. `NormalisedRotatingSource` from the exported tables, `gg′ = 0`, wrapped in
   `CoilAugmentedNormalisedSource` with the filaments subtracted;
   `setPlasmaSupportShape( Box, 1.5 )`.
3. `setAxisFluxPoint( 0.10, 0 )`, `setBoundaryFluxPoint( 0.35, 0 )`, and
   **`setBoundaryFluxInitial( ψ_vac( 0.35, 0 ) )`**: its default is 0 on a
   fresh solver, which makes the first span `ψ_ax − 0 > 0` and puts the box
   INSIDE the inner label, around the axis. The row is linear, so the first
   step pins `ψ_bnd` and the solve recovers to the same field — in 5 steps
   rather than 3. A cost, not a failure, and the test asserts both label rows
   read zero at iteration 0 so that it is paid nowhere.
4. `setSource( wrapper, ψ_vac( 0.10, 0 ) )`, the boundary data, extension and
   DtN at `N = 32`, `setInitialGuess( vacuum GridFunction )` — the GridFunction
   form, which seeds `q` consistently.

**The edge is live first**, not frozen: with `F` vanishing at both labels the
normalisation columns are analytic and complete, and the freeze exists for a
moving fill the box does not have. If a sweep is needed, it is judged settled by
the frozen `( ψ_ax, ψ_bnd )` repeating, **not by the component's element count,
which is identically zero under the box**. **Acceptance**: convergence from the
vacuum guess, quadratic on the best triple; `F` at each label below the same
1e-06; `checkAxialCutoff()` reported; the exterior coefficients' trace tail at
`N = 32` reported for the plasma's own current, which M-260 sized for the copper
only.

**MR-4c — MEQ against itself.** No closed form exists, so `ψ` at the machine
region's points converges in `h` and `k` at the rate the method promises, under
the two-tier rule. This is the stage that separates MEQ's discretisation error
from the cross-code difference MR-4d measures.

**MR-4d — MEQ against geq.** `ψ` at geq's grid nodes in the machine region,
nodes outside `D_h` dropped, and both labels, against the direct reference at
129², 257² and 513². **Acceptance, as a MATRIX** ([M-158](MEASUREMENTS.md#m-158)
is the precedent): the difference must FALL as geq refines, at geq's own observed
rate, and at 513² sit within a small multiple of geq's own 257-against-513
difference — only a difference that falls rules out a conversion error. And the
**control**: the same MEQ answer against the shipped-Hagenow reference reads the
1e-05 bias, which says the comparison resolves an error of that size.

#### Small corrections MR-4 lands, because its run would otherwise read them wrong — DONE

* `borderDiagnosis()` no longer reads a negative span as an empty support: a
  zero or non-finite span is flagged as `Ψ` undefined, a negative one is named
  as a mirror's sense, and emptiness is COUNTED — `plasmaSupportNodes()` asks
  the source's own gate at every potential node.
* `checkAxisSource()` decides `axisInsidePlasma` by the shape's own test,
  `NormalisedSource::supportAdmits()`, so a box whose `Ψ( R = 0 ) > 1` reads
  the axis outside it; the one-sided answer is unchanged.
  `MirrorSupport`'s `aMirrorsAxisIsOutsideItsBoxAndTheDiagnosticAgrees`.
* `ExteriorDtN::coefficients()` takes `max( 30, 2 × lastMode )` Gauss points,
  bit-identical up to 14 modes, and recovers a 48-mode trace to 8.9e-15 where
  the fixed rule read 4.2e-02.

#### What the run answered

1. **The bordered Newton converges with a negative span from a vacuum guess**,
   in three or four full steps at every `k` and level, no step damped.
2. **The `N_D` hold works inside geq**: 32 rebuilds refused in the 513² solve.
3. **The stopping target is not loose**: three orders tighter moves `ψ` by at
   most 1.8e-10 of the plasma's own flux, and on geq's tables by nothing.
4. **Zero initial exterior coefficients cost nothing**: the exterior merit is
   1e-16 after the first step, and a `setCurrentScale( 0 )` first solve gives
   an identical second history.

### MR-5 — from a file — DONE

`build/meq-run examples/mirror-rotating.toml` solves geq's rotating mirror.
`RotatingMirrorDriver` runs the shipped file and agrees with the library
solving geq's problem to 0.000e+00 over every dof, in the same four Newton
steps → **[M-264](MEASUREMENTS.md#m-264)**. The example names the reference
tables in `tools/geq-benchmark/reference/rotating/` rather than copying them.

**The schema.** Every new key states the problem; none is a solver knob:

* `[boundary.electrode] R, Z` is the inner label `ψ_inner`, the flux at the
  limiting electrode, through `setAxisFluxPoint()`. It has no starting value
  to give: the driver starts it at the starting state's flux there through
  `setAxisFluxInitial()`, and `[source] PsiAxis` and `[source.axis]` are
  refused beside the box.
* `[source] SupportShape = "box"` and `AxialEdge` give the support.
* `[mesh.generate] AlignPlasma` makes the box's edges element faces.
* The outer label is `[boundary.limiter] R, Z`, which already existed.

Beside the box, the parse refuses everything that would be accepted and
ignored, or that has not been measured on a mirror. Each refusal is asserted on
both `mhd` and `rotating`. `ConfigTests`' `mirror_keys` is the list.

**The four driver assumptions, resolved:**

* **`checkAxis()`** is not asked under `PrescribedPoint`.
  `checkSolvedAxis()` answers `std::nullopt`, for both the report and the
  annulus retry.
* **The starting labels** are the starting state's total flux at each label's
  point, on the box path only. It is a fixed rule, not a key. Spreading the
  `ψ_bnd` seed to a tokamak costs `limited-tokamak` 12 steps → 129, and
  `DriverAcceptance` now bounds that count.
* **The support loop** is refused beside the box. Under the box the fill is
  off, so the loop would read settled whatever the labels did. Judging it by
  the labels would need a tolerance, which would be a knob. The live edge
  converges.
* **The Picard pre-stage** is refused beside any inner label at a point, as
  measured.
  With the labels read at their points, every sweep's Picard-then-Newton
  diverges on the mirror, and a pinned tokamak exits 2.

**The `.nc` support mask** is `node_material`, `ELEMENT-REGIONS-PLAN.md`
ER-3 ([M-270](MEASUREMENTS.md#m-270)). **Adaptivity** runs on the mirror through
the driver, and its estimator tracks the flux error there.
`RotatingMirrorAdaptive` measures `η / ‖q − q_h‖` against a `k = 4` reference
at every cycle and on the uniform L1 mesh. It holds the index in
`AdaptiveRefinement`'s `( 5, 40 )`, reading 7.7, 13.6 and 13.1 over three
cycles, with `η` and the error both falling at every cycle
→ [M-279](MEASUREMENTS.md#m-279). Judging it found two defects of the
split, both repaired. `η₅` compared against the physical modal sum rather than
the remainder's datum. `η₁` read `F` at `ψ_p` rather than at `ψ_p + ψ_c`.

**`--device cuda` reproduces the host on the mirror.** M-276 raced it at 2.0e-15
([M-276](MEASUREMENTS.md#m-276)). `RotatingMirrorDevice` runs the shipped file
under `--device cpu` and `--device cuda`. It asserts `ψ_h` and `q_h` within
1e-12 of their largest values and the same Newton count, and skips loudly
without a card → [M-280](MEASUREMENTS.md#m-280).

**The bordered retry rung is skipped beside the box**, and a failed first
attempt exits 2 as it is → **[M-277](MEASUREMENTS.md#m-277)**. Run on the
mirror, it agrees with the first attempt to 7.3e-13 where both converge. But
from two of the four free `ψ_inner` starts that reach it, starts a file cannot
pose since `ψ_inner` begins at the electrode's flux, it converges to a different
equilibrium: a second discrete root with a spike in one element of the box, and
an answer accepted against a cold residual of 1e+12. The annulus trigger cannot
fire on a box run, since `[boundary.xpoint]` is refused beside it.
`RotatingMirrorDriver`'s `aMirrorWhoseFirstAttemptFailsIsNotRetried` asserts
the branch.

**Open:**
* **A second discrete root exists, and only a far `[initialguess]` can reach it
  now.** Under the mirror's stabilisation the plain bordered Newton started at
  `ψ_inner = 0.15` converges in 27 steps to `ψ_inner` = 4.5747e-03, `ψ_bnd` =
  4.980e-02, where the shipped answer is 5.5473e-03 and 6.530e-02 (M-294).
  `ψ_inner` now starts at the electrode's flux, 5.64e-03 on a cold start, so a
  file cannot pose that start; whether the root is M-277's 0.2 root, and
  whether it survives refinement, is untried.
* **The bordered target is relative to the cold residual**, so a bad enough
  start is accepted loose: from `ψ_inner = 0.5` a first residual of 1.1e+12
  accepts 9.1e-05. That bites any first attempt, not only the rung's, and is
  not fixed; on a mirror only a far `[initialguess]` now poses such a start.
* **The device on the mirror is measured in double only**: device against host
  in single, and `DeviceParity`'s `MEQ_AXIS_DEVICE=force` arm, are not run on
  it. Neither `--device` nor the adaptive loop has been run under
  `Gauge = "midplane"` or `"vacuum-midplane"`, and the knots' rows have no
  device twin, so under a Device their dot products take the host route.
* **The three-species pressure ceiling** (risk 5) was constant `τ = 1`'s
  element-local fold of the discretisation in every stall measured but one,
  and **a box support now defaults to `τ = 3/max(R, 0.1 m)`**, under which
  the coarsest mesh converges 2e19 in both midplane gauges and 5e19 under the
  midplane gauge, asserted. The one stall that is the plasma's is the
  vacuum-midplane gauge at 5e19: that posed problem folds, mesh converged,
  near 4.16e19, and the default reaches it to about 1% on level 0
  (→ [M-291](MEASUREMENTS.md#m-291), [M-293](MEASUREMENTS.md#m-293),
  [M-294](MEASUREMENTS.md#m-294)). **Open**: in single
  the format's floor rises with `τ`'s size — constant or scaled alike — so a
  solve with a meshed coil sits above the solver's `10 ε` absolute target
  under the default; the vacuum mirror's fixture asks single for 64 `ε`, while
  a meshed-pack mirror through the driver in a single build would meet the
  default target, which is unmeasured (the shipped files' filaments are
  subtracted, and they converge in two steps). The midplane gauge's posed-problem fold, if it has one, lies above
  6.4e20; and `q_h` on the rotating mirror is reference-limited on its finest
  rungs.
* **32 knots are fixed for every mesh.** At level 2 their error is 170× under
  the discretisation error, but between knots they move the
  `.nc`'s true-midplane densities off the input tables by up to 7e-05
  (carbon), and re-posing them agrees to 1.2e-06 of `ψ*`; finer meshes than
  level 2 are not checked against them ([M-283](MEASUREMENTS.md#m-283) §5,
  [M-284](MEASUREMENTS.md#m-284) §4).
* **A prescribed plasma current beside the knots is refused**: its corner
  against them is not assembled.
* **The vacuum-midplane gauge's label refusal has no end-to-end test**: the
  driver's exit 1 on a profile set above `meq::labelTolerance` is asserted only
  at library level, and that tolerance's cost in order at the label is argued,
  not measured (Newton is unaffected up to 1e-3).
* **Self-convergence on the vacuum's curve at 5e18 is slower than M-278's
  midplane run** (`ψ*` max 7.4e-07 then 5.2e-08, against 5.7e-07 then
  1.1e-08); the cause is not looked at ([M-281](MEASUREMENTS.md#m-281) §4).
* **The estimator's uniform effectivity settling is not asserted**: 7.7, 11.7,
  14.3 over levels 0 to 2 at `k = 2`, which is pre-asymptotic there.
* **A curved boundary with subtracted conductors and no `[boundary.exterior]`**
  takes η₅ against `transferredDatum()` with `g = 0`. That matches what the
  solve imposes; whether that combination poses the intended problem is
  unexamined.

## Risks, in the order they are likely to hurt

1. ~~**MR-0 fails.**~~ **RETIRED** by [M-169](MEASUREMENTS.md#m-169): zero
   unrepresentable current at `ρ_Γ = 1.7`, at every Mach number.
2. ~~**The mode count `N` is large.**~~ **MEASURED and moderate**: the copper
   sets it, not the Mach number — 24 to 32 at `ρ_Γ = 3.0`
   ([M-260](MEASUREMENTS.md#m-260)).
3. ~~**The support's element-level half.**~~ **DISSOLVED** by the axial cutoff:
   the support is a bounded box, so no connected component and no seed is needed
   and the fill is simply off here. What remains is that turning it off must be a
   **stated branch a test asserts**, not a flag that happens to be false — the
   failure this tree has already met twice.
3b. **The cutoff is a modelling approximation.** `z_edge` is prescribed, and a
   cutoff through a region still carrying current is a discontinuity in `F` that
   the profiles' clamping will not soften. The diagnostic in MR-1 is what keeps
   it honest.
4. **geq has moved on.** `kinetic_closure.py`, `fast_species.py` and an
   anisotropic effective-current hook are in geq now and MEQ has no counterpart.
   The benchmark is pinned to the isotropic `ProfilesCentrifugalMirror` path and
   must stay there; a geq default that changes under us is a stale-reference
   failure. geq's answers are regenerated rather than committed, so what pins
   them is the commit they are regenerated at: geq `73b9a72` and freegs4e
   `9a1ee0c` for every reference so far (`tools/geq-benchmark/README.md`).
5. **The gauge**, item 5. It is the only item on this list that is a change to a
   published interface rather than to solver internals, and the only one found by
   running the benchmark rather than by reading the code. **Not on MR-4's path**:
   its exponent is half what §5 first quoted, and two species at geq-Mach 2 to 4
   are comfortably representable. **BUILT** → [M-278](MEASUREMENTS.md#m-278):
   geq's mirror in the midplane gauge meets MR-4d's matrix, and three species at
   `M = 6` solve. **The freeze's price is paid and the freeze retired**: `R₀` at
   32 knots are border unknowns, so Newton is quadratic under the midplane gauge,
   in 3 steps from the vacuum at every pressure that converges
   ([M-283](MEASUREMENTS.md#m-283), [M-284](MEASUREMENTS.md#m-284)).
   **What remained was a pressure ceiling that is not the freeze's**: the same
   three species converged at 5e18 everywhere and at 2e19 from level 1, and
   stalled at 2e19 on level 0 and at 5e19 on levels 0 and 1, in every gauge and
   from every start measured, with the Jacobian exact and nothing moving
   ([M-285](MEASUREMENTS.md#m-285)). **It was the HDG stabilisation's**: one
   element's local block going singular at constant `τ = 1`, which rises with
   the mesh. **A box support defaults to `τ = 3/max(R, 0.1 m)`**, which lifts
   it — 5e19 converges on level 0 under the midplane gauge, asserted — and
   only the vacuum-midplane gauge's 5e19 is the posed problem's own fold, mesh
   converged near 4.16e19 ([M-291](MEASUREMENTS.md#m-291),
   [M-293](MEASUREMENTS.md#m-293),
   [M-294](MEASUREMENTS.md#m-294)).
6. **geq is a moving reference in a second sense**: shipped, it is biased by its
   own von Hagenow offset, and its source steps at both labels. MR-4a builds a
   reference that is neither, so it is a geq the shipped tools do not produce,
   and the export script is what makes it reproducible.
