# The plasma edge as an interior interface, coupled at a distance

**A design, not an implementation, and deliberately not to be built yet.**
Written 2026-09-05, out of FB-4's measurements. `FREE-BOUNDARY-PLAN.md` §5.3 and
§7.10 are the measurements; `CLAUDE.md`'s *The plasma edge caps the order* is
the record. This file is the one route out of that cap that fits the machinery
MEQ already has, worked out far enough to be costed and argued with.

**THE PRECONDITION IS EXPLICIT: get `j ≥ 1` right first.** Nothing here should
be started until the ordinary path is green and measured at `j = 1`.
`PlasmaEdgeConvergence` is green as of 2026-09-05 and is the first half of
that; what remains is FB-5's bordered solve carrying the moving support, and a
machine case against `../freegs4e` at `j ≥ 1`. There
is a real risk of building this instead of finishing that, and the numbers below
are exactly the argument for not doing so: **`j ≥ 1` already gives `k+2` at
`k ≤ j`, with nothing built.**

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
polynomials on a mesh that does not follow `Γ_p`, not about the equation. What
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
is exactly `−e_j` (`CLAUDE.md`, *The NPC port*). With the geometry frozen inside
the Newton that stays true. With the geometry live it would not.

## 6. Staging

Each stage ends at a **measured rate**, as every stage in this tree does.

| | | acceptance |
|---|---|---|
| **PE-0** | **A fixed interface, one side only.** Solve the plasma problem alone on `Ω_{p,h}` with `λ` GIVEN, against `PlasmaEdge`'s exact solution. No coupling, no vacuum side. | `ψ` at `k+1` and `ψ*` at `k+2` **at `j = 0`**, where the unfitted path reads 1.9 and 1.5. This is the whole claim of the plan, and PE-0 is where it is either true or not |
| **PE-1** | **The vacuum side, with its inward paths.** Same, on `Ω_{v,h}`, `λ` still given. | the same rates from the other side, plus the mirrored star-shapedness margin `≥ 0` |
| **PE-2** | **The two-sided sweep and its tiling check.** Both path families over the same `Γ_p`; weights from each must sum to `\|Γ_p\|`. | agreement to the cone-off floor, `~5e-10`, which is the number `ExtensionBoundaryQuadrature`'s own check reached. **Refine the RULE at fixed `h`** to separate coverage from smoothness — the mistake recorded in §7.7 of the free-boundary plan |
| **PE-3** | **The coupling closed on a FIXED interface**, `λ` solved for by the `N` flux-balance rows in a border. Still no moving support: prescribe `Γ_p`. | `k+1` / `k+2` against the exact solution at `j = 0, 1, 2`, and `λ`'s truncation shown to be spectral in `N` and therefore not the binding term |
| **PE-4** | **The support read off `ψ_h`**, with the outer fixed point of §5.2. | the same rates as PE-3 on `MovingPlasmaEdge`; the outer loop terminating on "no element changed"; **and `j = 0` converging at all**, which is the finding this plan exists to overturn |
| **PE-5** | **Into the free-boundary solve**: `ψ_ax`, `ψ_bnd`, the exterior coupling and the plasma edge in one bordered system. | `HighBetaConvergence` and the FB-1 cases unchanged to every digit, which is what says the generalisation reduces |
| **PE-6** | **A machine case at `j = 0`** against `../freegs4e`. | agreement at the level FB-6 reaches at `j = 1` |

**PE-0 IS THE STAGE TO PROTECT.** It needs no coupling, no border and no moving
geometry, and it either shows `k+2` at `j = 0` or shows that the whole premise
is wrong. It is a day's work against a fixture that already exists. **Do PE-0
before believing any of the rest of this file.**

## 7. Risks

**The X-point, and it is the one that could stop this being useful.** A diverted
plasma's edge is a separatrix with a **corner**, where `∇ψ = 0` and both
transfer-path families give out — `CLAUDE.md` records `ExtensionConvergence`
taking `Γ` to be `ψ = −0.03` rather than the separatrix for exactly this reason,
and `INVERSION-PLAN.md`'s IN-5 stops at the same wall. **So this design covers
LIMITER plasmas and not diverted ones**, which is a large scope limit and is not
a detail to discover at PE-5. `FREE-BOUNDARY-PLAN.md` §10 is the diverted-plasma
pathway and **it does not rescue this one**: it works by never meshing the
separatrix, which is available to a free-boundary solve and is exactly what an
interior interface cannot do. Whether an X-point can be handled by excluding a
disc around it, and what that costs the order, is unexamined.

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

**THE PREMISE ITSELF HAS ONE PIECE OF EVIDENCE AGAINST IT, AND PE-0 IS THE
TEST OF EXACTLY THAT.** This plan rests on *no element straddling `Γ_p` implies
full order*, which is why the uncut elements' own convergence matters.
`theCutCapsTheOrderBeforeAnyMethodIsChosen` measures it, and at `j = 0` and
`j = 1` it behaves — 2.88 / 3.88 / 4.98 and 2.96 / 3.72 / 4.76 against `k+2`.
**At `j = 2, k = 3` it does not**: the pair rate FALLS, 4.67 then 4.50, onto the
cap `j + 2.5 = 4.5` rather than climbing to `k+2 = 5`. A falling rate is not
pre-asymptotics, and the uncut elements carry nothing but smooth functions, so
either the classification of "uncut" is wrong or something couples the interior
to the band. **Until that is understood the plan's central premise is measured
at two rungs out of three.** PE-0 settles it directly and cheaply, which is the
reason it is first and the reason not to skip it.

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
