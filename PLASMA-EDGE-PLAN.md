# The plasma edge as an interior interface, coupled at a distance

**A design, not an implementation, and deliberately not to be built yet.**
Written 2026-09-05, out of FB-4's measurements. `FREE-BOUNDARY-PLAN.md` §5.3 and
§7.10 are the measurements; `CLAUDE_FB.md`'s *The plasma edge caps the order* is
the record. This file is the one route out of that cap that fits the machinery
MEQ already has, worked out far enough to be costed and argued with.

**THE PRECONDITION WAS EXPLICIT — get `j ≥ 1` right first — AND IT IS NOW MET.**
Nothing here was to be started until the ordinary path was green and measured at
`j = 1`. All three parts of that have landed: `PlasmaEdgeConvergence` green as of
2026-09-05, FB-5's bordered solve carrying the **moving support**
(`setPlasmaCurrent()` and `[source] ConfineToPlasma`, with the support a
connected component rather than a pointwise test), and a **machine case against
`../freegs4e`** — `examples/limited-tokamak.toml`, whose profile tables are
`Ψ²`, so `j = 2`.

**SO WHAT HOLDS THIS PLAN IS NO LONGER A PREREQUISITE. IT IS THE COST-BENEFIT
ITS OWN NUMBERS MAKE**, and that argument is unchanged: **`j ≥ 1` already gives
`k+2` at `k ≤ j`, with nothing built.** There is still a real risk of building
this instead of the things `ROADMAP.md` puts ahead of it — §10's diverted
plasmas among them, which this plan explicitly cannot cover, both transfer path
families giving out at a corner. **Its own central premise is now measured at all three rungs and it is TRUE** —
PE-0, §6 — and the same measurement says the plan's own machinery does not cash
it. **PE-0's `ψ*` acceptance is `k+2` where `λ` sits on `Γ_{p,h}` and `k+1`
through the transfer**, the latter being ExtensionConvergence's per-pair floor —
`ψ*` is never worse than `ψ_h`'s own design rate — applied across the sequence.

**AND THE TRANSFER DOES NOT COST AN ORDER, WHICH IS WHAT THIS PARAGRAPH USED TO
SAY.** → **[M-83](MEASUREMENTS.md#m-83)**. Read pair by pair at the fine end the
transferred `ψ*` gives **2.84 / 3.87 / 4.79** against `k+2` of 3 / 4 / 5: the
rate is intact. What the transfer costs is a **constant** — 30×, 93× and 555×
the fitted column at `k = 1, 2, 3` — because the lifting evaluates the element's
own degree-`k` polynomial *outside* its element along a path of length
`d ≈ 1.3 h`, and extrapolating a degree-`k` polynomial that far amplifies its
error by a factor exploding in `k`. A rate read endpoint to endpoint over a
finite range of `h` therefore reads `k+1`. The gate stays at `k+1` because that
is the only floor holding over the range this sequence spans; it is not a
statement that an order is gone. **The constant is set by `dist/h`**, so
`meq::AdaptiveDomain` is the lever here too.

**`k+1.5` IS THE OBVIOUS NUMBER AND IT IS MEASURED WRONG.** Mirroring
`ExtensionConvergence`'s *sequence* floor on the outer boundary would fail seven
of PE-0's nine rows: 2.23 against 2.35, 3.16 against 3.35 and 4.00 against 4.35
at `j = 0` alone. **[M-81](MEASUREMENTS.md#m-81)** also shows that no measuring
window on this mesh family rescues either column — refining floors the fitted
one at 3.6e-14 and stalls the transferred one three orders *above* the floor,
where `ψ*` and `ψ_h` land within 20% of each other.

## 1. What this fixes, and what it costs to leave unfixed

Free boundary makes the source stop on `Γ_p = ∂Ω_p`, a level set of the solution
that cuts through elements. The exact `ψ` then carries `|d|^{j+2}` across `Γ_p`,
where `j` is the order to which the profiles vanish at the edge — `p' ~ Ψ^j`,
a **modelling** choice. Measured (FB-4):

| | best approximation, any unfitted `P_k` | MEQ today |
|---|---|---|
| `ψ_h` | `min( k+1, j+2.5 )` | `min( k+1, j+1.5 )` |
| `q_h` | `min( k+1, j+1.5 )` | `min( k+1, j+1.5 )` — at its bound |
| `ψ*` | `min( k+2, j+2.5 )` | `min( k+2, j+2 )`, `j+2.5` with the rule raised |

So `ψ*` keeps `k+2` **exactly when `k ≤ j`**, and the cap is a theorem about
polynomials on a mesh that does not follow `Γ_p`, not about the equation.

**AND THE `ψ*` ROW OF WHAT THIS PLAN DELIVERS IS NOT `k+2`, WHICH PE-0
MEASURED.** Removing `j` from `ψ_h` works and is measured — 2.00 / 3.04 / 4.01 at
`j = 0` through the transfer, against the cut mesh's cap of 1.5 whatever `k` is.
`ψ*` reaches `k+2` only where `λ` sits on `Γ_{p,h}`; through the transfer it
reads 2.23 / 3.16 / 4.00, about `k+1` to `k+1.5`, which is the shortfall
`ExtensionConvergence` already asserts on the outer boundary. Read the rows below
as what an uncut *domain* admits, not as what this plan hands you.

What
this plan buys is the removal of the `j` in those expressions: with no element
straddling `Γ_p`, every subdomain solution is smooth up to its own boundary and
the rates go back to `k+1` and `k+2` **at every `j`, including `j = 0`**.

**What it costs to leave unfixed, stated honestly.** At `j ≥ 1` it costs
nothing a user of MEQ would notice below `k = j+1`. At `j = 0` it costs the
problem entirely: the assembled residual is *discontinuous* in the unknowns —
a fixed quadrature point crossing `Γ_p` makes `F` jump there — and Newton does
not converge from anywhere, the exact solution included. **That, and not the
order, is the case for doing this at all.**

## 2. Why this route and not the other two

Three routes beat an unfitted polynomial cap, and all three work by putting the
interface into the discretisation.

| | what it is | why not |
|---|---|---|
| **Fit the mesh** | remesh so element faces follow `Γ_p`, curved | `Γ_p` moves every iterate, so it is a remesh per iterate; the trace space changes shape; and a *straight-sided* fit caps at second order, which `refs/CutElementQuadratureSurvey.pdf` makes its central finding |
| **Enrich the space** | add `d_+^{j+2}` to `P_k` on cut elements (XFEM/GFEM) | cheapest in principle, and the exponent is known exactly. But HDG's `k+2` post-processing is a theorem about polynomial local problems, and it would have to be re-derived on an enriched space. A large piece of analysis MEQ cannot check by measurement alone |
| **Two subdomains, coupled at a distance** | **this plan** | it is `refs/CouplingAtADistance.pdf` applied to an interior interface rather than an exterior operator — the same paper, the same primitives, and MEQ has already built one instance of it |

**The third is the one that fits because MEQ has already built its harder
half.** Stage 5 solves on a subdomain of whole elements with a boundary
condition imposed on a *curved boundary the mesh does not follow*, transferred
along paths, and it delivers `k+1` in `ψ` and `k+2` in `ψ*`. FB-1 then couples
that subdomain to a **second operator** across a gap, with the coupling
unknowns carried in a border. This plan is those two things composed: replace
FB-1's exterior DtN with a second HDG problem, and the exterior boundary with
the plasma edge.

## 3. The formulation

**The transmission conditions are homogeneous, which is the piece of luck that
makes this tractable.** Across `Γ_p`, with `F` bounded, `ψ` is `C¹` and
`q = (1/r)∇̄ψ` is continuous — both components, not merely `q·n`. So there is no
interface data to model: the conditions are

```
[[ ψ ]] = 0        and        [[ q ]] = 0        on Γ_p
```

and the only reason the interface exists at all is that `∇̄·q` jumps.

**Three domains, two of them discretised.**

```
Ω_{p,h}   the union of elements lying ENTIRELY inside Ω_p.   Γ_{p,h} = ∂Ω_{p,h}
Ω_{v,h}   the union of elements lying ENTIRELY outside Ω_p.  Γ_{v,h} is its inner boundary
B_h       the band of elements Γ_p passes through.  DISCRETISED BY NEITHER
```

`B_h` is `O(h)` wide and is covered by the two subdomains' own polynomial
extensions, exactly as `Ω \ D_h` is in stage 5. **This is the "at a distance" of
the title**: the two solvers never share a face, and everything they say to each
other is said on `Γ_p`, which lies between them.

**The coupling unknown is modal, as FB-1's is.** `Γ_p` is a smooth closed curve
and `ψ|_{Γ_p}` is smooth on it, so a Fourier series in a curve parameter
converges spectrally and `N ~ log(1/h)` modes suffice. Let

```
λ( s ) = Σ_{m<N} λ_m C_m( s )        the common trace of ψ on Γ_p
```

Then the coupled system is

```
(a)   plasma HDG on Ω_{p,h},  with ψ = λ on Γ_p transferred from Γ_{p,h}
(b)   vacuum HDG on Ω_{v,h},  with ψ = λ on Γ_p transferred from Γ_{v,h}
(c)   N flux-balance rows:    ∫_{Γ_p} ( E_h(q_p) - E_h(q_v) )·ν C_m dΓ = 0
```

which is a Dirichlet-to-Neumann coupling: (a) and (b) each turn `λ` into a
Neumann trace, and (c) says the two agree. **It is FB-1's structure with the
exterior DtN replaced by a second interior solve**, and the count is `N`
unknowns against `N` equations.

**No mesh on `Γ_p` is needed and that is deliberate.** The integral in (c) is
swept by `mfem::ExtensionBoundaryQuadrature`, which is the routine MEQ wrote for
FB-1's transmission row and upstream merged; its acceptance is that the boundary
weights sum to `|Γ_p|`, which is also the statement that the two sides' path
families cover the interface once each. A mortar space on `Γ_p` is the
alternative and is **not** preferred: it needs the crossing points of `Γ_p` with
every element face, which is the geometry a cut rule needs, and it costs
`O(1/h)` dofs where the modal route costs `O(log 1/h)`.

**And the `1/r` is not there, for the same reason it is not there in
`exteriorTransmissionRows()`**: MEQ's `q` *is* `(1/r)∇̄ψ`, so contracting it
against the plain measure already carries the radius. Writing `dΓ/r` divides by
it twice. Expect to get this wrong once anyway.

## 4. What MEQ already has

Almost all of it, which is the argument for this route.

| piece | where | what changes |
|---|---|---|
| a subdomain of whole elements | `meq::BoundaryShape`, `mfem::SubMesh` | needs to mark **two** subdomains from one level set, and the vacuum one has a hole |
| transfer paths | `mfem::TransferPath`, `VertexConePath` | the vacuum side's paths point **inward**; the cone's orientation is the detail to get right |
| the datum on a curved boundary | `setExteriorDatum()`, `mfem::PathTraceCoefficient` | used twice, once per side, with `λ` as the `PositionFunction` |
| the extension of the flux past an element | `mfem::ElementExtension` | unchanged |
| a quadrature sweeping the true boundary | `mfem::ExtensionBoundaryQuadrature` | used twice; the tiling check becomes a **two-sided** one |
| the transmission row | `exteriorTransmissionRows()` | becomes a difference of two extensions rather than one |
| the projection `λ → trace` | FB-1's `P` | unchanged in shape |
| a bordered solve for functional unknowns | `solveWithNormalisation()` | `2×2` today, `(2+N)×(2+N)` here |
| the estimator on a transferred datum | `transferredDatum()`, `η₅` | must see the interface faces on both sides |

**Two things are genuinely new**, and they are the ones to cost:

* **N1. A vacuum subdomain with a hole in it.** `SubMesh::CreateFromDomain` on a
  marked attribute gives it, but every consumer that assumes a simply connected
  `D_h` has to be checked, `VertexConePath`'s cone construction included. Note
  `CLAUDE.md`'s standing trap: **a `SubMesh` keeps a pointer to its parent**, so
  the background mesh must outlive both.
* **N2. Paths pointing inward.** Stage 5's family lifts *outward* from an
  inscribed `Γ_h` to a circumscribed `Γ`. The vacuum side lifts *inward*. The
  analysis is symmetric; the implementation's signs and the star-shapedness
  margin are not, and `gammaHIsStarShapedAboutTheCentre…` is the case to
  mirror.

## 5. The moving interface, which is the whole difficulty

**`Γ_p` is a level set of the solution, so every piece of geometry above is
rebuilt as the iterate moves.** Three consequences, in increasing order of how
awkward they are.

**5.1 The geometry's derivative enters the Jacobian.** This is
`FREE-BOUNDARY-PLAN.md` §6.4's "one real gap", and it is not avoided by this
route — it is the *price* of every route to high order at a moving interface.
MEQ avoids it today only because a fixed quadrature rule has no geometry to
differentiate.

**5.2 The discrete problem is piecewise constant in the geometry, which is
worse than merely non-smooth.** As `Γ_p` moves, an element flips between `B_h`
and `Ω_{p,h}`, and the *spaces themselves* change. No amount of care makes that
a differentiable function of the unknowns: the coupled residual has jumps of
size `O(h²·F)` scattered through coefficient space, of the same species as the
quadrature jumps measured at `j = 0` and a factor of about `h⁻¹` fewer.

**So the outer iteration on the geometry must be a FIXED POINT, not a Newton**,
and that is a design decision rather than a fallback:

```
repeat
    freeze  Ω_{p,h}, the paths and the sweeps at the current psi_h
    solve   the coupled system by Newton, with the geometry CONSTANT
until   the support stops changing
```

The inner Newton is then differentiating a smooth residual and keeps its
quadratic convergence; the outer loop converges linearly, and terminates
**exactly** rather than to a tolerance, because the support is a finite set of
elements. That last point is worth stating: the natural stopping rule is *no
element changed subdomain*, which is checkable and cannot oscillate if the
support is monotone. **Whether it oscillates in practice is unmeasured and is
PE-4's acceptance criterion.**

**5.3 It costs the `ψ_ax` border its exactness.** Under NPC the axis border row
is exactly `−e_j` (`CLAUDE_HDGGS.md`, *The NPC port*). With the geometry frozen inside
the Newton that stays true. With the geometry live it would not.

## 6. Staging

Each stage ends at a **measured rate**, as every stage in this tree does.

| | | acceptance |
|---|---|---|
| **PE-0** | **A fixed interface, one side only.** Solve the plasma problem alone on `Ω_{p,h}` with `λ` GIVEN, against `PlasmaEdge`'s exact solution. No coupling, no vacuum side. | **RUN.** `tests/convergence/PlasmaEdgeConvergence.cpp`, three cases, `j = 0, 1, 2` × `k = 1, 2, 3`; the control is the same solve with `λ` deleted and it is flat at −0.06. **`ψ` at `k+1`: MET.** At `j = 0` it reads **2.00 / 3.04 / 4.01** where the cut mesh is capped at `j + 1.5 = 1.5` whatever `k` is, and the cap is gone at every `j`. **`ψ*` at `k+2`: NOT MET with `λ` transferred** — 2.23 / 3.16 / 4.00 at `j = 0` — and **MET with the same `λ` imposed on `Γ_{p,h}` instead**, 2.88 / 3.94 / 4.93. The two columns differ only in the transfer, so **the premise is true and the TRANSFER is what costs `ψ*` its extra order**, by about a full order, at every `j`. `j` has left the rates: the spread of `ψ*` over `j = 0, 1, 2` is 0.07 / 0.18 / 0.29 against the cut mesh's 1.11 / 2.04 / 2.18. §7's `j = 2, k = 3` worry does NOT reproduce here — see the note under this table |
| **PE-1** | **The vacuum side, with its inward paths.** Same, on `Ω_{v,h}`, `λ` still given. | the same rates from the other side, plus the mirrored star-shapedness margin `≥ 0` |
| **PE-2** | **The two-sided sweep and its tiling check.** Both path families over the same `Γ_p`; weights from each must sum to `\|Γ_p\|`. | agreement to the cone-off floor, `~5e-10`, which is the number `ExtensionBoundaryQuadrature`'s own check reached. **Refine the RULE at fixed `h`** to separate coverage from smoothness — the mistake recorded in §7.7 of the free-boundary plan |
| **PE-3** | **The coupling closed on a FIXED interface**, `λ` solved for by the `N` flux-balance rows in a border. Still no moving support: prescribe `Γ_p`. | `k+1` / `k+2` against the exact solution at `j = 0, 1, 2`, and `λ`'s truncation shown to be spectral in `N` and therefore not the binding term |
| **PE-4** | **The support read off `ψ_h`**, with the outer fixed point of §5.2. | the same rates as PE-3 on `MovingPlasmaEdge`; the outer loop terminating on "no element changed"; **and `j = 0` converging at all**, which is the finding this plan exists to overturn |
| **PE-5** | **Into the free-boundary solve**: `ψ_ax`, `ψ_bnd`, the exterior coupling and the plasma edge in one bordered system. | `HighBetaConvergence` and the FB-1 cases unchanged to every digit, which is what says the generalisation reduces |
| **PE-6** | **A machine case at `j = 0`** against `../freegs4e`. | agreement at the level FB-6 reaches at `j = 1` |

**§7'S `j = 2, k = 3` ANOMALY IS THE SUBDOMAIN'S STAIRCASE AND NOT THE PLASMA
EDGE, AND PE-0 SEPARATES THE TWO WITH A THIRD CASE.** The same equilibrium on a
**fitted rectangle** strictly inside the disc — no subdomain, no staircase, no
transfer, and no plasma edge in the domain at all — reads, at worst over the
three rungs, **2.99 in `ψ` and 3.99 in `ψ*`** at `k = 2` and **3.98 and 4.94** at
`k = 3`. So the
plasma-side solution admits `k+1` and `k+2` at every `j`, which is the premise
as a statement about the SOLUTION, and it is not in doubt. On `Ω_{p,h}` with the
same exact trace on its own boundary it reads 3.59 / 4.64 at `j = 2, k = 3` —
short, and the shortfall is what the inscribed staircase boundary costs. Raising
the source rule from `2k + 4` to `2k + 16` moves that column by not one digit,
so it is not the quadrature; the fitted rectangle rules out the fixture and the
solver. The obvious
reading, and it is a reading rather than a measurement: HDG's `k+1` in `ψ_h` and
`k+2` in `ψ*` come from a duality argument that wants the adjoint problem
`H²`-regular, and `Γ_{p,h}` has re-entrant corners. **Measured, not assumed,
and the phrase "a staircase of 270° corners" that stood here was quad language
on a simplicial mesh**: the background is `Element::TRIANGLE`, right triangles of
45/45/90 with six meeting at an interior vertex, so a boundary vertex's interior
angle is a sum drawn from { 45°, 90° }. At `n = 64` the interface carries 22 at
90°, 16 at 135°, 54 at 180°, **12 at 225° and 20 at 270°**. → **[M-84](MEASUREMENTS.md#m-84)**.

**AND THE TRANSFERRED ROUTE IS INTERMITTENT, ON A SMOOTH CIRCLE.** §9.4 records
a mesh-dependent fragility and attributes it to a corner. PE-0 meets it with no
corner anywhere: at `k = 3`, sweeping `n = 16, 24, 32, 48, 64, 96, 128`, the
transferred `L2(ψ)` reads 3.39e–07 at `n = 16` and **5.89e–07 at `n = 24`**, and
`L2(q)` improves by 0.67 of an order between `n = 64` and `n = 96` — with
`VertexConePath::NumWidened()` **zero** at every mesh and
`dist(Γ_{p,h}, Γ_p)/h` in **[1.02, 1.33]**, so P.1 holds throughout. Raising the
path quadrature to order 24 changes the answer in the sixth figure. **So the
corner sharpens that fragility rather than causing it**, and a per-pair rate is
not assertable on the transferred route — which is why PE-0's pairwise tier
asserts that the error FALLS rather than asserting a rate.

**AND A DIVERTED CASE NEEDS A STAGE BEFORE PE-0, WHICH §9.4 HAS ALREADY RUN
ONCE.** PE-0's interface is smooth, so it is silent about the corner; the
lens/union sweep of §9.4 is the corner's own premise test, it needs no coupling
and no border, and it says the corner is not free. Extending it to say *where*
the error lives is what decides whether an `O( h )` exclusion around the node
is affordable — see §9.5.

**PE-0 IS THE STAGE TO PROTECT.** It needs no coupling, no border and no moving
geometry, and it either shows `k+2` at `j = 0` or shows that the whole premise
is wrong. It is a day's work against a fixture that already exists. **Do PE-0
before believing any of the rest of this file.**

## 7. Risks

**The X-point, and it is the one that could stop this being useful. SECTION 9 IS
NOW THE ACCOUNT OF IT AND THIS PARAGRAPH'S DIAGNOSIS WAS HALF WRONG.** A
diverted plasma's edge is a separatrix with a **corner**, and this paragraph
used to say that is because `∇ψ = 0` and both transfer-path families give out.
Two claims, and only the second is measured: what the families meet is a
**degenerate and double-valued defining function**, not a corner, and
`ExtensionConvergence.cpp`'s own header says so precisely — a ray leaves the
core, crosses the vacuum and lands in the private flux region *where `ψ` is
negative again*. That is fixable, and §9.2 says how.

**What is NOT fixable by a better path family is a corner between two
transferred pieces**, and §9.4 measures what one costs: on a smooth exact
solution, with the corner as the only variable, `ψ*` falls back to `ψ` and `q`
to first order or worse — so a corner takes away exactly the two quantities this
plan exists to recover. **So this design still covers LIMITER plasmas and not
diverted ones**, which is a large scope limit and is not a detail to discover at
PE-5; but the reason is now a number rather than a hand-wave, and §9.5 lists
three geometries that would remove it. `FREE-BOUNDARY-PLAN.md` §10 is the
diverted-plasma pathway and **it does not rescue this one**: it works by never
meshing the separatrix, which is available to a free-boundary solve and is
exactly what an interior interface cannot do.

**The band can exceed one element.** Where `Γ_p` runs nearly parallel to a mesh
line, `dist(Γ_{p,h}, Γ_p)/h_loc` grows and assumption P.1 fails. Stage 5's
answer is `meq::AdaptiveDomain`, the companion mesh of GS-2 §3.3, and the same
answer should work here — but it was written for a *fixed* `Γ` and would have to
be rebuilt per outer iterate.

**The two sides can disagree about which point of `Γ_p` they are talking
about.** Both sweeps must parametrise the interface consistently, or the
flux-balance row contracts one side's `s` against the other's. The tiling check
of PE-2 does **not** catch this: two families can each cover `Γ_p` once and
still disagree pointwise. A separate check is needed — the obvious one is to
transfer a known non-constant function from both sides and difference it, which
is the two-sided analogue of `theTransferredDatumReproducesTheImposedCondition`.

**`N` is a new tuning parameter and this project distrusts those.** Spectral
convergence means it should be small and insensitive, but "should be" is not a
measurement. PE-3's acceptance is written to force it.

**The outer fixed point may oscillate.** Two elements swapping in and out on
alternate iterations is the obvious failure, and it has the same shape as the
period-4 limit cycle `CLAUDE.md` records for NPC on §4.5. There is no
globalisation for a combinatorial iteration; the answer if it happens is
probably to require the support to grow monotonically within an outer sweep.

**THE PREMISE HAD ONE PIECE OF EVIDENCE AGAINST IT AND PE-0 IS THE TEST THAT
DISPOSED OF IT.** This plan rests on *no element straddling `Γ_p` implies full
order*, which is why the uncut elements' own convergence matters.
`theCutCapsTheOrderBeforeAnyMethodIsChosen` measures the best approximation, and
at `j = 0` and `j = 1` it behaves — 2.88 / 3.88 / 4.98 and 2.96 / 3.72 / 4.76
against `k+2` — while **at `j = 2, k = 3` its pair rate FALLS**, 4.67 then 4.50,
onto the cap `j + 2.5` rather than climbing to `k+2 = 5`.

**PE-0 SOLVES THE SAME EQUILIBRIUM AND THE FALL DOES NOT REPRODUCE.** On a
FITTED rectangle strictly inside the plasma — no cut, no staircase, no transfer
— `j = 2, k = 3` reads 3.98 / 3.98 / **4.94**, and on `Ω_{p,h}` with the exact
trace the per-pair `ψ*` reads 4.70, 4.50, 4.74: scatter, and rising rather than
falling. Raising the source rule from `2k+4` to `2k+16` moves that column in **no
digit**, so it is not the quadrature either. **The premise is true at all three
rungs**; what is left of the shortfall belongs to `Γ_{p,h}`'s re-entrant 225° and 270°
re-entrant corners, where HDG's duality argument wants an `H²`-regular adjoint
and does not get one. §6's PE-0 row has the numbers.

**And this is a research project, not a port.** Unlike stage 5, there is no
paper that has done exactly this, no published table to reproduce, and no
reference implementation to check against. Every number will be MEQ's own. That
is an argument for PE-0's exact solution being the acceptance and for not
skipping it.

## 8. What this deliberately does not do

**It does not remove `j` from the physics.** A profile with `p'(0) ≠ 0` still
puts a jump in the current density at the plasma edge; this plan represents that
jump instead of smearing it. Whether a user wants such an equilibrium is a
modelling question and this changes nothing about it.

**It does not help the fixed-boundary path.** There `Γ` is given and stage 5
already does this.

**It is not a cut quadrature and does not become one.** If what is wanted is
*solvability* at `j = 0` rather than high order, a cut volume rule plus the
matching cut surface rule for the Jacobian's `∮ F φ/|∇ψ|` term is a much smaller
piece of work with a much smaller prize — about second order — and it should be
costed separately rather than folded in here.

**And it should not be started while `j ≥ 1` is unfinished.** See the top of this
file.

## 9. The diverted case, and the one measurement that says what it costs

**§7 records the X-point as this plan's scope limit and gives one sentence of
diagnosis: *"a separatrix with a corner, where `∇ψ = 0` and both transfer-path
families give out"*. That is two claims welded together, only one of them is
right, and neither is the thing that actually decides the case.** This section
is the measurement.

### 9.1 The separatrix is four analytic arcs and a node, not a singular curve

**`ψ` is analytic at an X-point** — `FREE-BOUNDARY-PLAN.md` §10.1's organising
fact — because `Δ*ψ = −F` with `F` bounded and the X-point is an ordinary
non-degenerate saddle. The Morse lemma then gives coordinates in which

```
ψ − ψ_X = a( x ) · b( x )        with ∇a, ∇b independent at the node
```

so each of the four branches is an **analytic arc terminating at the node**.
Nothing about the separatrix is singular; the *description* `ψ − ψ_X` is.

**Measured on MEQ's own diverted equilibrium**, `Soloviev::nstx()` evaluated in
closed form:

```
X-point   r = 0.699700   z = -1.716000   psi_X = -2.8e-17   |grad psi| = 1.6e-11
axis      r = 1.318168   z = +0.011089   psi   = -2.663e-01
Hessian eigenvalues  -0.4049  +0.6291
```

`ψ_X = 0` to round-off, which is a free confirmation that the **corrected**
`c₇`/`c₁₀` put the separatrix where a Solov'ev normalisation says it belongs —
`nstxAsPublished` puts it at −8.7e-3 and its zero set is not closed. The
branches leave the node at 36.11°, 113.59°, −66.41° and −143.89°, and

| | interior angle at the node | `π/ω` |
|---|---|---|
| **core** (the plasma; here `{ψ < ψ_X}`) | **77.48°** — convex | 2.323 |
| **vacuum**, taken as one region | **282.52°** — re-entrant | 0.637 |
| SOL lobe, if the node is cut | 102.5° | 1.756 |

and the two boundary branches separate at **1.25 × arc length** from the node.

### 9.2 What defeats the path families is the DEFINING FUNCTION, not the corner

`mfem::TransferPath` takes a `PositionFunction φ` whose zero set is `Γ`, and the
instinct is `φ = ψ − ψ_X`. That is what breaks, in two ways at once: `|∇φ| → 0`
linearly at the node, so a root find along a ray is degenerate; and there are
**two** roots within `O(h)`, so the nearest one is not the right one.

**And this is already measured, in this tree, in a comment nobody had connected
to §7.** `tests/convergence/ExtensionConvergence.cpp`'s header:

> *"LevelSetPath aborts because the outward normal below the plasma tip never
> meets the level set — it runs straight through the X-point into the
> private-flux region, where psi is negative again — and VertexConePath aborts
> for the same reason, its whole admissible fan missing Gamma."*

**A ray leaving the core, crossing the vacuum and landing in the private flux
region is a double-root problem and a connectivity problem. It is not a corner
problem.** `LevelSetPath` even carries `search_steps` for "where `Γ` can be
crossed twice within `search_length`", which is the same defect in its milder
form.

**The fix follows from §9.1**: give each branch its own defining function,
`φ₁ = a`, `φ₂ = b`, each with `|∇φ| ≠ 0` at the node — in practice a signed
distance to each traced arc. The node is located by
`meq::CriticalPointFinder` as a root of `q_h`, which is a **solved** field at
the potential's own order (IN-A measured the ITER saddle to **4.5e-6**), and the
launch directions are the null directions of the Hessian there.

### 9.3 The connectivity fix is a REGULARITY prerequisite, not a physics nicety

`{ψ < ψ_X}` near the node is **two opposite sectors meeting at a point** — the
core and the private flux region — and that set is not locally Lipschitz, so no
extension analysis applies to it at all. Take the core alone and both domains
become Lipschitz: 77.48° convex, 282.52° re-entrant.

So `FREE-BOUNDARY-PLAN.md` §10.3's flood fill is not a later refinement for
diverted plasmas. **It is what makes `Ω_p` a domain this plan can be posed on**,
and the plan's dependency on it should be stated at the top rather than
discovered at PE-5.

### 9.4 THE MEASUREMENT: a corner is not free, and it is not the singularity

**The design.** `Γ` is the intersection (a **lens**) or union of two equal discs,
and the exact solution is the Solov'ev NSTX equilibrium — smooth everywhere and
knowing nothing about the shape. So the corner is the only variable: no level
set of `ψ`, no vanishing gradient, no private flux region, no connectivity. The
datum `g = ψ_exact` is carried in by `setExteriorDatum()`. Two circles of radius
`a` with centres `2e` apart meet at

```
lens    omega = pi - 2 asin( e/a )        union   omega = pi + 2 asin( e/a )
```

and `e/a = 0.7800` gives **77.48°** and **282.52°** — the two angles of §9.1,
from one parameter.

**What the theory predicts, and why it is the wrong prediction.** For a Dirichlet
corner of interior angle `ω` the singular exponents are `mπ/ω`, so a general
solution carries `ψ ~ r^{π/ω}` and `q ~ r^{π/ω − 1}`, capping the `L2` rates at
`min( k+1, 1 + π/ω )` and `min( k+1, π/ω )`. **The exact solution here is
smooth, so that mode's coefficient is zero and no cap applies to it** — exactly
as it is zero at an X-point, where `ψ` is smooth across its own separatrix. Any
order loss measured here is therefore the METHOD exciting a mode the solution
does not have, and that is the only mechanism that carries over.

**The control is clean at every degree.** One disc, no corner:

| `k` | `ψ` | `q` | `ψ*` | target |
|---|---|---|---|---|
| 1 | 1.98, 1.99 | 2.00, 1.98 | 3.30, 2.97 | 2, 2, 3 |
| 2 | 2.97, 2.99 | 2.90, 2.75 | 4.34, 3.80 | 3, 3, 4 |
| 3 | 4.53, 4.02 | 4.32, 3.25 | 5.65, 4.39 | 4, 4, 5 |

**AND THE SHARPEST CONTROL IN THE STUDY IS THE SAME LENS WITH THE CORNER
ROUNDED.** `smax( u, v ) = ( u + v + √( ( u − v )² + ρ² ) )/2` at `ρ = 0.05` —
nine cells at `n = 32` and thirty-six at `n = 128`, so every mesh resolves it.
Same two discs, same thin geometry, same box, same datum, same path family; the
corner is the only thing removed:

| `k` | `ψ` | `q` | `ψ*` | target |
|---|---|---|---|---|
| 1 | 1.99, 1.99 | 2.54, 2.06 | 2.75, 3.08 | 2, 2, 3 |
| 2 | 3.17, 3.09 | 3.80, 2.80 | 3.83, 4.76 | 3, 3, 4 |
| 3 | 4.80, 6.37 | 4.39, 4.40 | 4.80, 6.36 | 4, 4, 5 |

**Full order at every degree, and five orders of magnitude smaller**: `L2( ψ )`
at `k = 3, n = 128` reads **2.1707e-12** rounded against **2.5437e-07** sharp.
So it is the corner. Nothing else in the configuration is the variable, and no
argument is needed.

**The corner is not.** At `ω = 77.48°`, `n = 32, 64, 128`:

| `k` | `ψ` rate | `q` rate | `ψ*` rate |
|---|---|---|---|
| 1 | 1.88, 2.11 | **0.53, 1.60** | 1.44, 2.63 |
| 2 | 2.88, 1.09 | **1.33, 0.70** | 2.85, 1.13 |
| 3 | 4.81, **−0.94** | **2.58, −0.72** | 4.80, **−0.94** |

`L2( q )` sits at **1e-4** where the control reaches 1e-10, and at `k = 3` the
error **rises** between the two finest meshes. **`ψ*` equals `ψ` to three
figures throughout** — 1.7969e-06 against 1.8019e-06 at `k = 2, n = 32` — so
the local post-processing, which is the whole of this plan's `k+2` claim, buys
**nothing at all** on a domain with a corner.

**AND IT IS NOT AN ORDER CAP. IT IS A MESH-DEPENDENT CATASTROPHE, WHICH IS
BOTH BETTER AND WORSE NEWS.** Widen the corner and the failures do not become
milder, they become **intermittent**:

| shape | `k` | `L2( q )` at `n = 32, 64, 128` |
|---|---|---|
| lens, 150° | 2 | 1.2261e-06 → **4.2838e-06** → 6.5829e-09 |
| lens, 150° | 3 | 1.6201e-08 → **1.0559e-04** |
| lens, 120° | 1 | 3.2358e-05 → **1.8812e-04** → 1.7624e-05 |
| union, 210° | 3 | 8.3199e-08 → 1.1503e-08 → **4.5115e-06** |
| union, 210° | 1, 2 | clean: rates 1.99/1.98 and 3.03/2.30 |

`6.58e-09` at 150°, `k = 2`, `n = 128` is **better than the control's
1.70e-08** — so the method is fully capable of `k+1` at that corner and simply
fails on particular meshes. That is a **geometric fragility**, not a regularity
bound, and at 77.48° it has stopped being intermittent and become the rule.

**Three things that are ruled out by the same runs**: `VertexConePath::NumWidened()`
is **0** at every mesh, every angle and every degree, so the admissible fan
never had to be widened; `dist(Γ_h, Γ)/h` stays in **[1.17, 1.33]**, so
assumption P.1 holds; and the control shares the mesher, the datum route, the
path family and the exact solution, so none of those is the variable.

**And this tree already carries the discriminating control from the other
side.** `theSolverReachesTheExteriorDatumOnTheHalfDisc` solves on a semicircle
whose `Γ` meets the axis at **two right-angle corners**, and it reads
**1.99 / 2.99 / 3.99** in `ψ`. The difference is that those corners are between
a **transferred** piece and a **fitted** one. So what is fragile is specifically
a corner between two TRANSFERRED pieces, which is exactly what a separatrix
node is.

### 9.5 What this means for the plan

**It is not a scope note any more, it is a numbered risk with a number on it.**
This plan exists to recover `k+2` in `ψ*` and `k+1` in `q` that a cut element
costs; a corner in `Γ` takes `ψ*` down to `ψ` and `q` to first order or worse.
**Trading the cut-element cap for a corner cap is not obviously a trade worth
making**, and PE-0's premise test is silent on it because PE-0's interface is
smooth.

**So the staging gains a stage before PE-0 for the diverted case, and it is
cheap**: the lens/union sweep above, extended to say *where* the error lives.
If the loss is confined to the `O(h)` neighbourhood of the node, the cost of
excluding it is estimable — `q·ν ~ distance` on a branch, so an excluded disc of
radius `Ch` contributes `∫₀^{O(h)} O(s) ds = O(h²)` to a flux-balance row, which
is affordable at `k = 1` and **not** at `k ≥ 2`. If it is not confined, the
diverted case needs a different geometry rather than a better path family.

**Three geometries, in increasing order of cost.**

1. **Interface strictly inside, at `ψ = ψ_X − ε`.** Smooth, closed, no node, and
   everything above becomes irrelevant. It gives up representing the current
   jump exactly at the edge, which is what this plan is for — so it is a
   fallback, not a design, and it is named here only so nobody rediscovers it.
2. **A mesh node at the X-point.** This is what the SOL community does with
   block-structured flux-aligned grids, and it converts the node into two
   ordinary corners of the mesh with nothing to transfer across. It was
   unavailable when this plan was written and is available now that FB-6 meshes
   with gmsh — and §5.2's outer loop already freezes the geometry per iterate,
   so a re-mesh per outer iteration is not absurd.
3. **Four subdomains meeting at the node.** The four sectors are core 77.5°,
   private flux 77.5°, and two SOL lobes of 102.5° — **every one convex**, every
   `π/ω > 1.75`. Treating the vacuum as one region is what manufactures the
   282.52°. This is the structural reason block-structured divertor grids exist,
   and it costs two more solves and a four-way flux balance at the node.

### 9.6 The modal border must be per-arc, and the plan inherited its smoothness from the wrong parent

**§3 expands `λ( s )` in one curve parameter and claims spectral convergence with
`N ~ log( 1/h )`. That claim comes from `CouplingAtADistance.pdf`, whose `Γ` is
an ARTIFICIAL interface chosen smooth** — a circle, in its §4.2 — and which says
so: *"the BEM is defined on a suitably chosen, **smooth** artificial interface"*,
*"let `x : R → Γ` be a **smooth regular**"* parametrisation. They pick a smooth
`Γ` precisely to avoid this question. It is the **extension** work,
`HDG-CurvedExtensions.pdf` and the Joukowsky aerofoil of its §3.4, that is built
for a sharp `Γ` — and `VertexConePath` exists because both closed-form families
give out there.

**So the transfer half of this plan has a parent that meets corners and the
modal half does not**, and a global Fourier series on a closed curve with a
corner converges at `O( 1/m² )`, not spectrally: `N ~ h⁻¹`, and §3's whole cost
argument against a mortar space collapses.

**The parametrisation that works** follows from §9.1 and needs nothing exotic,
because each arc is analytic in arc length right up to the node:

* **one arc per branch, normalised arc length from the node**, with the node's
  position from `CriticalPointFinder` and the launch directions from the
  Hessian's null directions;
* **a per-arc polynomial basis** — Chebyshev or Legendre in `t ∈ [0,1]` — with
  continuity of `λ` at the node as one linear constraint. Spectral per arc,
  algebraic globally, which is why the split is not optional;
* **testing weighted by the endpoint behaviour**: `q·ν = |q| ~ distance` on a
  branch, so the flux-balance integrand vanishes linearly at the node and
  `t·C_m`, or Jacobi `P^{(0,1)}`, stops the near-node modes being
  ill-determined. Measure it rather than assuming it;
* **branch assignment for every `Γ_h` face before any foot is searched for**, so
  a ray cannot cross to the other branch — which is the failure
  `ExtensionConvergence.cpp` recorded.

**One cap that is not removable.** The analytic Morse lemma wants `ψ` analytic.
In the **vacuum** it is (`Δ*ψ = 0`). In the **core** `ψ` is only as smooth as `F`
allows: with `F ~ Ψ^j` extended by zero, `F ∈ C^{j−1}` and `ψ ∈ C^{j+1}` — `C³`
at `j = 2`. So the per-arc basis converges algebraically at a rate set by `j`,
not spectrally, and `N` is capped by the same regularity that caps everything
else here. Consistent with §5 and with FB-4, and better stated now than
discovered at PE-3.

### 9.7 And there is a better answer than excluding the node

Inside the `O(h)` disc where the two branches are closer together than the band
is wide, you have a **closed-form local solution**: `ψ ≈ ψ_X + ½ ξᵀHξ`, three
numbers plus the node position, all available from `CriticalPointFinder` at the
flux's own order. The interface data there can be **computed rather than
transferred**. **The X-point is the one place on the separatrix where more is
known, not less**, and treating it as the hard part inverts that.
