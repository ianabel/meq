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

**FREEZE `R₀( ψ )` WITHIN A NEWTON ITERATION.** The tree's own precedent is
`setPlasmaSupportFrozen` and XP-2's frozen plasma edge: a reference that moves
inside a Jacobian makes a difference straddle a kink instead of measuring a
derivative. Refresh it once per iteration, beside the support freeze.

**THE INTERFACE SHAPE** is a gauge policy on the source rather than a second
class — `enum class Gauge { FixedRadius, Midplane }`, `FixedRadius` carrying
today's `referenceRadius` and remaining the default so every tokamak caller is
untouched, `Midplane` carrying the tabulated `R₀( ψ )` the solver refreshes.

**AND IT MAKES THE BENCHMARK EASIER, WHICH IS A CHECK ON THE DESIGN RATHER THAN
A BONUS.** geq's `get_flux_density_qn` already references each surface to its own
midplane radius, so under the fixed-`z` gauge MEQ's `n_s0` and geq's `N_s` are
the same quantity up to geq's single global offset at `ψ_n = 0.5`. The transfer
in `tools/geq-benchmark/export_mirror.py` becomes near-identity, and the
reference should be regenerated under the new gauge once this lands —
`theGaugeTransferIsLoadBearing` will then be asserting something much weaker and
should be re-pointed at the fixed-`R` gauge, which is where the trap still lives.

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

* `[source.axis] R, Z` is the inner label, through `setAxisFluxPoint()`.
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
* **The starting `ψ_bnd`** is the starting state's total flux at the outer
  label, on the box path only. It is a fixed rule, not a key. Spreading it to a
  tokamak costs `limited-tokamak` 12 steps → 129, and `DriverAcceptance` now
  bounds that count.
* **The support loop** is refused beside the box. Under the box the fill is
  off, so the loop would read settled whatever the labels did. Judging it by
  the labels would need a tolerance, which would be a knob. The live edge
  converges.
* **The Picard pre-stage** is refused beside any `[source.axis]`, as measured.
  With the labels read at their points, every sweep's Picard-then-Newton
  diverges on the mirror, and a pinned tokamak exits 2.

**Not done:**

* the `.nc` `plasma` support mask (the docs warn that `n_s` and `e_phi_0` are
  evaluated outside the support);
* adaptivity, `--device` and the bordered retry rung, which are accepted beside
  the box but have not been run on a mirror.

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
   are comfortably representable.
6. **geq is a moving reference in a second sense**: shipped, it is biased by its
   own von Hagenow offset, and its source steps at both labels. MR-4a builds a
   reference that is neither, so it is a geq the shipped tools do not produce,
   and the export script is what makes it reproducible.
