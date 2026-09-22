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
[M-169](MEASUREMENTS.md#m-169). The machine ends at `ρ = 1.63` and the plasma
with it, so a semicircle at `ρ_Γ = 1.7` has **exactly zero** current beyond it at
`M = 2` as at `M = 6`. There is no Mach threshold for admissibility.

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
| conductors outside `Γ` | FB-7's `meq::ExteriorCoilSet` / `setExteriorConductors()`. `ψ = ψ_coil + ψ̃` with `Δ*ψ_coil = 0` inside, so the coils enter **only through `Γ`**, additively and known, with no new unknowns. Mirror coils sit at `ρ ≈ 1.7` where the plasma is inside `ρ ≈ 1.5` |
| the bordered Newton | the `( N + 2 )` elimination already carries `ψ_ax`, `ψ_bnd` and the exterior coefficients as one system |
| `ψ_bnd` as an unknown pinned at a POINT | `setBoundaryFluxPoint( R, z )`. Under NPC the row is exactly `−e_j` and the corner exactly 1, neither differenced. **This is the mechanism the mirror needs twice** |
| meshing | `tools/mesh/halfdisc.py`, a semicircle reaching the axis with coils meshed to, driven from `[mesh.generate]` |
| `η₆` | the exterior transmission indicator, for when the coefficients stall while the interior refines |

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

### 3. ~~XP-1's element flood fill has nothing to seed from~~ — dissolved by item 2

The element-level half of the support takes the connected component containing
the magnetic axis. A mirror has no axis, and a flux-labelled annulus is connected
*through the tube* to the ends of the machine — so the fill neither starts nor
terminates where it should.

**THE AXIAL CUTOFF REMOVES THE QUESTION.** With `|z| ≤ z_edge` the support is a
bounded annular box, defined pointwise and needing no seed and no connected
component. So the fill is simply not wanted on this path — which is a decision
the plan can now take rather than a `2 × 2` it has to run.

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

### 5. `normalisationDerivatives()` is not overridden for the rotating closure

Recorded in `CLAUDE_FLOW.md` and still open: the base returns false, so the
bordered Newton **differences** that column, and a difference perturbs the
normalisation, which moves the edge, so it straddles a kink instead of measuring
a derivative. With both flux labels pinned at prescribed points the support is
far less free than in a tokamak free-boundary run, so this degrades rather than
fails — but it is on the path and it is the one piece of
`meq::NormalisedMHDSource`'s treatment the rotating source still lacks.

### 6. The single-`R_ref` gauge does not survive a high-Mach multi-species mirror

**FOUND BY REGENERATING THE SOURCE BENCHMARK AT `M = 6`, AND IT IS ON THE PATH
TO MR-4.** `meq::Species::density` is `n_s0( ψ )` — the physical density on ONE
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

**The plasma-current border is not needed.** `CLAUDE_FB.md` records that with a
fixed profile amplitude a confined equilibrium is a non-linear *eigenvalue*
problem and the coupled solve fails from everywhere without `setPlasmaCurrent`.
That argument turns on the support being free. Here **both** labels are pinned at
prescribed points, so the support is not free and the eigenvalue does not arise —
geq prescribes physical densities and temperatures and never a current. Expect
to need no current border, and treat needing one as a finding.

**`g = 0` is fine.** `F = μ₀ R² ∂p/∂ψ` with `∂p/∂ψ ≠ 0` at the inner label, so
the trivial branch that `examples/rotating-rectangle.toml` warns about is not
reached. The vacuum coil field is the initial guess.

## The stages

Each ends at a measured acceptance, not at "it runs".

### MR-0 — pick `ρ_Γ` and `N`, and keep the diagnostic

**THIS STAGE WAS THE CAMPAIGN'S FALSIFYING EXPERIMENT AND
[M-169](MEASUREMENTS.md#m-169) HAS RUN IT.** The exterior is admissible at every
Mach number with `Γ` outside the machine, so what is left is small and is
engineering rather than physics:

* **A choice**: `ρ_Γ ≳ 1.7` for this machine, and `N = 32` for `M = 2` or 16 for
  `M ≥ 4`. Both are existing parameters.
* **A MEQ-side diagnostic**: the fraction of `∫|F|/R` lying outside `Γ`,
  reported on the converged answer. It is the mirror's analogue of
  `checkAxisSource()` — one-sided, refusing or warning on a positive detection —
  and it exists because `M` is chosen while the current on `Γ` is solved for.

**Acceptance**: the diagnostic reads below `1e-3` on a shipped mirror
configuration, and a mode-count sweep on MEQ's own solve reproduces
[M-169](MEASUREMENTS.md#m-169)'s trend. **The risk this stage existed to retire
is retired**, so it no longer gates the others; MR-1 and MR-2 can start now.

### MR-1 — the plasma support becomes a box in `( Ψ, z )`

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

### MR-2 — `AxisConstraint::PrescribedPoint`

`ψ_ax` pinned at a prescribed point, reusing `ψ_bnd`'s border row.
**Acceptance**: on an existing *tokamak* fixture, prescribing the point at the
known axis reproduces `LocatedAxis`'s answer to round-off. **It needs no mirror
at all**, which is what makes it a safe stage to land early.

### MR-3 — the vacuum mirror, through MEQ

No plasma. Mirror coils as `meq::ExteriorCoilSet` outside `Γ`, half-disc mesh,
DtN live. **Acceptance**: MEQ's `ψ` against the coil field computed directly from
`meq::coilPsi()`, and against geq's own vacuum field. This validates the
geometry, the mesh, the coil model and the exterior **together and without the
plasma**, which is the rung that does not exist today between the analytic
fixtures and a full machine solve. It also sizes the mode count `N` that an
elongated machine needs on `Γ` — the number MR-0 bounds and this one measures.

### MR-4 — the free-boundary rotating mirror

Everything at once: rotating source, both labels prescribed at points, two-sided
support, exterior coils, DtN, bordered Newton. **Acceptance**: agreement with
geq's converged rotating mirror, in the same mode MEQ solves in — `CLAUDE_FB.md`'s
standing rule. Expect the comparison to be limited by geq's `φ₀` convergence
before it is limited by MEQ; [M-167](MEASUREMENTS.md#m-167) measured that at
5.5e-05 rms in `Σ Z_s n_s`, and the tightened-reference trick is available here
too.

### MR-5 — from a file

`[source] Type = "rotating"` already parses. What is new is the axis constraint's
prescribed point, the two-sided support, and a `[mesh.generate]` mirror geometry.
**Acceptance**: an `examples/mirror-*.toml` that the shipped binary solves, with
a driver acceptance in the shape of `MachineFixedBoundary` — the driver on the
shipped file, because what rots about an example is the file.

## Risks, in the order they are likely to hurt

1. ~~**MR-0 fails.**~~ **RETIRED** by [M-169](MEASUREMENTS.md#m-169): zero
   unrepresentable current at `ρ_Γ = 1.7`, at every Mach number.
2. ~~**The mode count `N` is large.**~~ **MEASURED and small**: 32 modes at
   `M = 2`, 16 at `M ≥ 4`. MR-3 confirms it on MEQ's own solve.
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
   failure, which is why the reference is committed.
5. **The gauge**, item 6. It is the only item on this list that is a change to a
   published interface rather than to solver internals, and the only one found by
   running the benchmark rather than by reading the code.
6. **`normalisationDerivatives()`**, item 5. Degrades the Newton tail rather than
   the answer, and is visible as a linear rate — `CLAUDE_FB.md`'s *a linear rate
   on a Newton method is a Jacobian statement*.
