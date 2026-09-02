# Solution inversion: `ψ(R, z)` to `R(Ψ, l)`, `z(Ψ, l)`

Written 2026-09-02, after a literature survey whose references are indexed in
`refs/Refs.md` under *Solution inversion*. `CLAUDE.md` is the operational record
and is authoritative on anything technical; `ROADMAP.md` is the priority order.
This file is the design for one item.

**The problem.** meq solves for `ψ` on a mesh. Almost everything downstream
wants the inverse: the flux surfaces themselves, as curves parametrised by a
flux label and an in-surface coordinate, and the integrals taken over them.

**What it unblocks**, and this is why it is worth doing before either of the
things that need it:

* **`ROADMAP.md` item 10, the fixed-`q(ψ)` solver.** RoPP (142) is
  `q(ψ) = V′(ψ) I(ψ) ⟨r^{-2}⟩_ψ / 4π²` — a flux-surface average, and the
  machinery `FLOW-PLAN.md` §3.3 deliberately avoided needing.
* **`MANTA-COUPLING.md` §5 and §6**, whose own checklist calls
  `dGeometry_dpsi` **with the moving-surface term** and *making `Geometry` fast
  enough for a pointwise call pattern* "the two genuinely hard ones. Everything
  else is bookkeeping." Both are this.
* **`DRIVER-PLAN.md` §3's flux-surface `(Ψ, θ)` output grid**, deferred there as
  "worth doing; not worth blocking the driver on" precisely because it needs
  contour tracing.

**The prior art is `v0-legacy:FluxSurfaces.cpp` and it has never compiled.**
`main()` declares `unsigned int N_surfaces =` with no initialiser and uses an
undeclared `psi_values` four times. The `Trace` function is sound and is
textbook Allgower–Georg; there is an algorithm to reuse and no working tool to
extend. Read it before starting — `git show v0-legacy:FluxSurfaces.cpp`.

---

## 1. What changes, and what does not

**Nothing about the solve.** This is post-processing. The operator, the
discretisation, `τ`, the hybridization, the estimator, the adaptive loop and the
curved boundary are all untouched, in the same way `FLOW-PLAN.md` was a change
to `F` alone. What is new is a consumer of `ψ_h` and `q_h`.

**Four new pieces**, in increasing order of how much thought they need:

| | | |
|---|---|---|
| a tracer | predictor–corrector on `ψ_h = c`, with the tangent and the corrector both from `q` | **ordinary** — the prior art is most of it |
| critical points | the magnetic axis and, when there is one, the X-point, as roots of `q = 0` | **small, and it is a prerequisite for the tracer rather than a refinement of it** |
| a representation | `R(Ψ, l)`, `z(Ψ, l)` as a differentiable object, not a point cloud | **the design decision**, §4 |
| flux-surface averages | the integrals the whole thing exists for | **the accuracy question**, §3.3 |

---

## 2. The three-way error split, which is the frame for everything below

Any extraction has three independent errors and conflating them is how this
problem gets done badly. Written for a contour of `ψ_h` at level `c`:

| | what it is | size |
|---|---|---|
| **(a) field** | `{ψ_h = c}` against `{ψ = c}` | `O(h^{k+1} / \|∇ψ\|)`, by the implicit function theorem |
| **(b) point location** | how far a computed point is from `{ψ_h = c}` | corrector: machine precision. Marching squares: `O(h_cell²)` |
| **(c) representation** | how far the interpolant *between* points is from the curve | governed by the point spacing `Δs` and the interpolation order |

**(a) is the discretisation and is not ours to improve here.** Note its
`1/|∇ψ|`: the geometric error of a contour is worst exactly where the gradient
is small, which is the magnetic axis and the X-point. This is not a defect of
any method; it is the implicit function theorem, and it is why §5 treats those
two points as objects rather than as hard cases.

**(b) is what the corrector buys**, and it is the reason to prefer this family
outright — see §3.1.

**(c) IS DECOUPLED FROM `h` AND FROM `k` ENTIRELY.** `Δs` is a free parameter
and each extra point costs a few Newton steps. This is the single most
important realisation in the survey, because the obvious objection to a traced
contour — "exact points joined by chords still give a second-order line
integral" — is *true and not binding*. To make the chord error subdominant to
the field error needs `Δs ≲ h^{(k+1)/2}`; at `k = 3`, `h = 0.05` that is about
four points per element edge, affordable even naively.

**And `q` improves it for free.** See §3.2.

---

## 3. The method, and why not the alternatives

### 3.1 Predictor–corrector continuation

The contour is the solution set of one scalar equation in two unknowns,
`ψ_h(x) − c = 0`, a 1-manifold wherever `∇ψ ≠ 0`. Predict along the tangent,
correct back with the minimum-norm Newton step:

```
t     = ( −∂_z ψ, ∂_R ψ ) / |∇ψ|          the tangent -- and r q IS grad-bar psi
x*    = x_i + Δs t                        predictor (explicit Euler)
x    <- x + ∇ψ ( c − ψ(x) ) / |∇ψ|²       corrector, iterate to tolerance
```

That is exactly what `v0-legacy:FluxSurfaces.cpp` does. Allgower & Georg is the
canonical treatment (`10.1137/1.9780898719154`); ch. 3 is the corrector, ch. 5
the convergence, ch. 6 steplength, ch. 8 bifurcation. **The whole book is on
disk at `~/AllgowerBook/`, chapter by chapter.**

**AND CH. 15 IS THE ONE TO READ BEFORE §4, FOR A REASON THAT IS NOT OBVIOUS
FROM ITS TITLE.** "Approximating Implicitly Defined Manifolds" is about
`K`-dimensional manifolds with `K ≥ 2`, where the earlier chapters are the
`K = 1` case. That looks like it is about surfaces in higher dimension and
therefore not about us. It is about us, via one device:

```
H( R, z, c ) := ψ_h( R, z ) − c,        H : R³ → R¹,      N = 1,  K = 2
```

**The whole family of flux surfaces is a single implicitly defined 2-manifold**
`M = H⁻¹(0)` in `(R, z, c)` space, and ch. 15's predictor–corrector machinery
applies to it directly. That is exactly the object §4 wants a representation of
— `R(Ψ, l)`, `z(Ψ, l)` is a chart on `M` — rather than a pile of separately
traced curves.

**And the device removes the critical points.** `∂H/∂c = −1` everywhere, so
`∇H ≠ 0` *always*: `H` has **no** critical points, and `M` is smooth even over
the magnetic axis and over an X-point, where the individual level sets are not
1-manifolds at all. What degenerates there is not the manifold but the
**projection** of it onto `(R, z)`, and the `(Ψ, l)` chart with it — `l` is
undefined at the axis.

That is worth stating precisely because it is the cleanest justification for
§4.1: the object is smooth, and what needs care is the *chart*, which is a
chart on a disc with a coordinate singularity at its centre. **That is exactly
what the Zernike basis is for**, and it is a better argument for it than "DESC
does it".

**THE PREDICTOR'S ORDER BUYS STEP LENGTH, NOT ACCURACY.** The corrector is a
root-find, so every accepted point sits on the discrete level set to solver
tolerance whatever the predictor did. There is nothing for a higher-order
predictor to buy, which is why Euler-plus-Newton is the norm and why the
first-order predictor in the prior art **is not a defect**. An earlier reading
of that file listed it as one; that was wrong.

**The corrector is what matters, and the reason is the error model, not the
constant.** Without it, error off the contour *accumulates* with path length —
you end up on a different contour and the curve does not close. With it, each
point is independently on the curve and there is no accumulation at all. A
first-order predictor with a corrector beats a fifth-order integrator without
one. Allgower & Georg ch. 6 states the same preference from the other side:
*"the general opinion is that it is preferable to exploit the contractive
properties of the zero set `H⁻¹(0)` relative to such iterative methods as those
of Newton type."*

**A corollary worth having: ODE tracing with projection and predictor–corrector
continuation are the same algorithm**, reached from opposite directions. There
is no third option hiding there.

**And it makes the DG jump a near-non-issue.** `ψ_h` disagrees with itself by
`O(h^{k+1})` across a face. A tracer that steps across one lands in whichever
element it lands in and the corrector re-anchors it there. Compare an ODE
integrator, for which a face is a discontinuity violating the smoothness every
Runge–Kutta order derives from.

> **Steplength: the standard strategies optimise the wrong objective for us.**
> Allgower & Georg ch. 6 details two, both keyed on **corrector effort** —
> Georg (1983), "asymptotic estimates in the mentality of initial value
> solvers", and Den Heijer & Rheinboldt (1981), "an error model for the
> corrector iteration". Continuation software wants to *traverse a path
> cheaply*, so equidistributing corrector work is right for it. **We want a
> well-shaped curve, so we want to equidistribute interpolation error**, which
> means curvature control, `Δs ∝ 1/κ`. Take the geometric variant and say why in
> the code.

### 3.2 The two nearly-free uses of `q`, and the trap between them

This is where meq has something no other code in the survey has: **the gradient
is a solved unknown at the same order as the potential**, not a differentiated
one.

**(i) Cubic Hermite between traced points.** `q` gives the unit tangent at
every point at `O(h^{k+1})`. Cubic Hermite through two points with two tangents
is `O(Δs⁴)` in geometry and in line integrals, is purely local — no global
spline solve — and costs nothing, since the tracer already evaluates `q` for its
predictor and its corrector. The spacing requirement drops from
`Δs ≲ h^{(k+1)/2}` to `Δs ≲ h^{(k+1)/4}`, about one point per element at
`k = 3`. **This is the highest-value single use of `q` in the whole item.**

**(ii) For a CLOSED surface, spectral quadrature.** A closed flux surface is
analytic away from the separatrix. For a smooth **periodic** integrand the
equispaced trapezoidal rule converges **geometrically, not algebraically**
(Trefethen & Weideman, `10.1137/130932132`). So the chord objection does not
merely become manageable for closed surfaces; it disappears.

> **THE TRAP, AND IT IS THE ONE MOST LIKELY TO BE WALKED INTO.** It is easy to
> build a spectrally accurate rule and then feed it a **second-order Jacobian**
> obtained by differencing neighbouring node positions — at which point the
> whole scheme is second order and nothing in the output says so. The
> arc-length element must come from `q`, pointwise, at full order.
>
> The survey found this stated nowhere. It is the same species as the band
> continuation of `B` recorded in `CLAUDE.md`: a quantity that is available at
> full order being obtained by differencing something instead, with no failing
> test to show for it. **IN-1 exists to measure exactly this**, with the
> differenced version kept as the control.

### 3.3 If the deliverable is an integral, consider not extracting at all

Implicit quadrature builds a rule on the level set directly. Saye
(`10.1137/140966290`) reports convergence tests **to twentieth order with
strictly positive weights**, needing only 1-D root finding and 1-D Gauss
quadrature — the root finding being a corrector in disguise.

**The caveat is decisive for meq: Saye 2015 is for hyperrectangles and meq is
on unstructured triangles.** The simplex-friendly line is Fries & Omerović
(`10.1016/j.cma.2016.10.019`) and Fries & Schöllhammer part II
(`10.1016/j.cma.2017.07.037`), the latter specifically on integration *on* the
extracted manifold; moment-fitting (Müller, Kummer & Oberlack,
`10.1002/nme.4569`) is cell-shape agnostic and is a third route.

**But it returns a rule, not a curve.** `MANTA-COUPLING.md` needs
`dGeometry_dpsi`, and `DRIVER-PLAN.md` §3 wants a `(Ψ, θ)` grid — both want the
curve. So implicit quadrature is a **cross-check on IN-2's integrals**, not the
primary route. Worth building as a control precisely because it is independent.

### 3.4 What is rejected, and why

| | |
|---|---|
| **Marching squares / marching triangles** | Caps at `O(h_cell²)` regardless of `k` — Etiene et al. measured vertex position 1.94, normals 0.93, curvature *divergent* (`10.1109/tvcg.2009.194`). On a simplicial mesh it is at least unambiguous, the affine interpolant's level set being convex hence connected; but the order is the point and it throws away `k−1` of them. |
| **Subdivide, then march** — what ParaView, VisIt, VTK and MFEM's own exporters do | **Second order in the SUB-CELL size with no `k`-dependence at all**; Remacle et al.'s error falls as `2^{2r}` and the degree never enters. Fine for a picture, which is what it is for. Not fine for `q(ψ)`. |
| **Exact rational contours** | Wiley et al. give the contour of a bivariate quadratic over a triangle as an exact conic. **It caps at `k = 2` in principle**: a level set of a degree-`d` polynomial is an algebraic curve of degree `d`, and rational parametrisation needs genus 0 — a generic smooth plane cubic has genus 1. meq runs to `k = 4`. (Their *quartic* comes from curved element geometry; meq's solve mesh is straight-sided, so it would be a conic here.) |
| **Plain ODE tracing without projection** | Error accumulates monotonically; the curve drifts onto a neighbouring contour and does not close. With projection it *is* §3.1. |
| **Ray bisection from the axis** | Point accuracy is as good as a corrector's and it parallelises trivially — worth keeping as an **independent cross-check** and as the initialiser for a periodic parametrisation, since it produces points already indexed by angle. **Not the primary route**: it assumes star-shapedness, which fails on indented cross-sections, at and beyond the separatrix, and near the axis where the bracket degenerates. |
| **The Morse–Smale complex** | **It gives the wrong curves.** Its "separatrices" are integral curves of `∇ψ`, *orthogonal* to the level sets; the plasma separatrix is the level **set** through the X-point. Recorded because it is the natural thing to reach for by name. |

---

## 4. The representation — the design decision

Tracing produces points. `MANTA-COUPLING.md` needs a **derivative with respect
to a `ψ` that moves the surface**, which a point cloud cannot supply and which
re-tracing at a perturbed `ψ` and differencing supplies badly. So the output of
this item is a **differentiable parametrised object**, and choosing it is the
substance of the design.

### 4.1 Closed surfaces: Fourier–Zernike

Nested closed surfaces filling in to a single maximum are a **disc**,
topologically, and the map `(Ψ, l) → (R, z)` is a map from a disc. The natural
basis is the one DESC uses for the same reason
(<https://desc-docs.readthedocs.io/en/stable/theory_general.html>):

```
R( ρ, θ ) = Σ_{l,m} R_lm Z_l^m( ρ, θ ),      Z_l^m = R_l^m( ρ ) { cos, sin }( m θ )
```

In 2D axisymmetric the toroidal Fourier factor of DESC's basis drops out and
what remains is the Zernike disc basis, whose **own** structure is the coupling
that matters: `R_l^m` contains only `ρ^|m|, ρ^{|m|+2}, … ρ^l`. That
parity-and-minimum-power constraint is what makes the magnetic axis a regular
point instead of a coordinate singularity to special-case. A naive tensor
product of "polynomial in `ρ`" times "Fourier in `θ`" admits modes that are not
analytic at `ρ = 0`; Zernike excludes them by construction.

> **THE RADIAL COORDINATE MUST BE `ρ = √Ψ_N`, NOT `Ψ_N`.** DESC's theory page
> says `ρ` is "the square root of the normalized toroidal flux, which is
> proportional to the minor radius", and the reason bites us identically:
> **`ψ` has a quadratic maximum at the axis**, so `ψ − ψ_ax ∝ (distance)²`.
> Parametrise by `Ψ` directly and `R(Ψ, l)` acquires a square-root branch point
> at the axis; every basis then converges algebraically, worst near the axis,
> and no convergence table says why. The Zernike degree counts powers of `√Ψ_N`.

**There is low-order precedent for the whole approach.** MXH (Arbon, Candy &
Belli, `10.1088/1361-6587/abc63b`, in `refs/`) is a Fourier-in-angle
parametrisation of a closed surface, and its §5 measures **exactly IN-2's
quantity** — the flux-surface average `⟨B²⟩` against the number of shape
coefficients — reporting rapid convergence at about half the coefficients a
plain Fourier expansion in `R` and `z` needs. Miller is the ancestor.

### 4.2 Open surfaces: Chebyshev in `l`

Deferred with free boundary (§7, IN-5), but the choice is settled. Open
surfaces terminate on the domain boundary, so they have genuine endpoints and a
periodic basis is simply wrong. Chebyshev is the natural basis on an interval,
and there is a second reason beyond non-periodicity: an open surface approaching
the separatrix becomes **stiff near its ends**, because arc length diverges
logarithmically as it nears the X-point. Chebyshev clusters resolution exactly
there.

### 4.3 The `ψ`-varying element: THE OPEN DECISION

**This is unresolved and this file does not pretend otherwise.** Nothing in the
survey addresses it, because nobody else is trying to satisfy all three of these
at once:

1. **Axis regularity.** Wants a global Zernike in `ρ`, which gives it free.
2. **A element boundary AT the separatrix.** A global expansion in `ρ` assumes
   smoothness all the way out, and that dies at the separatrix: the surface
   develops a corner at the X-point and flux-surface quantities diverge
   logarithmically. **No polynomial basis in `Ψ` converges against a
   logarithm.**
3. **C¹ in `Ψ` at minimum**, because `MANTA-COUPLING.md` needs
   `dGeometry_dpsi`. This is a floor, not a nicety, and it eliminates some
   element choices outright.

(1) and (2) pull opposite ways. The natural resolution is a **disc element in
Zernike containing the axis, surrounded by annular elements in something else**
— a known construction, but a real decision rather than a default.

**Resolve it by measurement, not by argument** — that is what IN-4 is for, and
it is deliberately placed *after* IN-2 so that there is a measured quantity
(`q(ψ)` against a closed form) to judge the candidates by. Candidates:

| | |
|---|---|
| global Zernike | simplest; axis free; expected to fail approaching the separatrix |
| Zernike disc + Chebyshev annuli, C¹ matched | the expected answer |
| global Zernike with the last element boundary short of the separatrix | accept a gap, measure how much it costs |

---

## 5. Critical points are objects, not hard cases

The magnetic axis and any X-point are where `∇ψ = 0`, which is where the
tracer's tangent is undefined and its corrector divides by zero. **Locate them
first, explicitly, and start the tracer from knowledge of them.**

**Note this is a statement about tracing ONE contour, not about the family.**
Under §3.1's lift the family is a smooth manifold through both points; it is the
`(Ψ, l)` chart that degenerates. So there are two defensible architectures — trace
contour by contour and treat the critical points as exceptions, or work on `M`
directly and treat them as ordinary interior points of a smooth surface whose
chart is awkward. **This plan takes the first**, because it reuses the prior art
and because IN-2's integrals are per-surface anyway; the second is the thing to
reach for if IN-4 finds the chart is the binding constraint.

**`q` gives the residual directly.** `∇̄ψ = r q`, so `q = 0` is a 2×2 system in
which the residual is a *solved field* and only the Hessian needs differencing.
Compare CEDRES++, which records as an open problem that in P1 continuous
Galerkin the axis and X-point are confined to mesh vertices; meq's high-order
representation resolves them sub-element, which is what makes polynomial root
finding worth doing at all. TokaMaker (`10.1016/j.cpc.2024.109111`) notes the
same for Lagrange order ≥ 2: saddles "can exist anywhere within the mesh".

**Exhaustiveness comes from subdivision with the Bernstein convex-hull test**,
and the triangle-native version is the one to use: Reuter et al.
(`10.1007/s00371-007-0184-x`, in `refs/`) work in the **barycentric** Bernstein
basis, so the subdivision primitive is a triangle and meq's mesh is native to
it. If all Bernstein coefficients of a component share a sign on a sub-cell, no
zero lies there — discard; recurse otherwise. Only provably-empty regions are
ever discarded, so nothing can be missed.

> **TOPOLOGICAL DEGREE IS A CERTIFICATION TEST AND NEVER AN EXCLUSION TEST.**
> The degree on a cell boundary gives the **sum of the indices** of the zeros
> inside, not a count. Two roots of opposite index cancel, so **degree zero does
> not imply no root.** Use the hull test for exhaustiveness and degree only to
> certify and count what survives. This was got wrong once during the survey
> and is written down so it is not got wrong again.

**A cheap global audit, with a stated blind spot.** Σ indices = χ(Ω) costs only
boundary evaluations and catches a missed or mis-signed critical point. It is
**blind to spurious critical-point pairs** — a spurious maximum (+1) with a
spurious saddle (−1) sums to zero — and numerical noise creates them strictly in
pairs. Pair it with a persistence threshold, which need not be tuned: the
stability theorem (Cohen-Steiner, Edelsbrunner & Harer,
`10.1007/s00454-006-1276-5`) bounds every spurious feature's persistence by
`2‖ψ_h − ψ‖_∞`, a quantity meq already measures to convergence-rate precision.
**The index check and persistence are complementary, not redundant.**

---

## 6. Scoping: the global-structure work is deferred, and the reason is sharp

The original framing of this item asked for open-versus-closed classification
and saddle crossing. **For the fixed-boundary problem meq solves today, there
are no interior saddles at all**, so none of that apparatus is exercised.

The argument, and note that **the load-bearing step is analytic, not
topological**:

* Poincaré–Hopf on a disc gives `#max + #min − #saddle = χ = 1`. That is an
  *identity*, not an absence — two maxima and one saddle also sums to 1.
* What excludes an interior minimum is the **maximum principle**. With
  single-signed `F`, `−∇̄·((1/r)∇̄ψ) = F/r ≥ 0` makes `ψ` a supersolution of a
  uniformly elliptic divergence-form operator (`1/r` is bounded away from zero
  and infinity since `r > 0`), so the minimum sits on `Γ`.
* Hence `#min = 0`, `#saddle = #max − 1`, and one axis gives **zero interior
  saddles**.

The conclusion is immune to all three of meq's sign conventions: `ψ → −ψ` swaps
max and min and both have index +1 in 2D; `q = ∇̄ψ / r` has the same zeros and
indices since `r > 0`; and `DarcyForm` holding `−q` changes nothing because in
even dimension `index(−v) = index(v)`.

**FREE BOUNDARY IS WHERE THIS CHANGES, AND IT CHANGES BY FORCE.** The vacuum
region between plasma and vessel makes the domain an **annulus**, `χ = 0`, so
`#max + #min − #saddle = 0` and with one axis a saddle **must** exist. X-points
are not incidental to free boundary; they are forced by the topology of the
domain. That is a much better reason to defer this work than "we have not seen
one yet".

**Two carve-outs to keep in view.**

* **GS-2 §4.4, the current hole.** Reversed core current destroys the
  supersolution property, so the "no interior minimum" step fails and an
  interior saddle becomes possible. It is already meq's unsolved benchmark for
  unrelated reasons (`max|∂F/∂ψ|/λ₁ ≈ 26`, multi-valued).
* **Rotating runs where the object of interest is `p` rather than `ψ`.** In a
  rotating equilibrium `p` is **not** a flux function — that is the entire
  physical content of the model — so it has its own critical points in its own
  places, and contours of `p` and of `ψ` are different curve families. Nothing
  above transfers. `RotatingSource` exists, so this is live rather than
  hypothetical.

**One thing that lands in the current tree regardless.** `refs/Refs.md` records
that with HDG-GS-2's **published** coefficients the X-point sits at
`ψ = −8.7e-3`, so the zero level set is **not a closed curve**. "Is this contour
closed?" already has a wrong answer in one of meq's own fixtures.

**When the topology work does come**, at free boundary, the shape is settled:
the **join tree only** (a superlevel-set question — the split tree and the merge
are not needed, which halves the work), output-sensitive (Chiang et al.,
`10.1016/j.comgeo.2004.05.002`, sorts only the critical vertices, so with
`t ≈ 3` the sort disappears), on a gluing rule derived from the HDG trace `ψ̂`
rather than a lossy filter, with critical points from §5. Carr, Snoeyink & Axen
(`10.1016/s0925-7721(02)00093-7`) is the construction; Nucha et al.
(`10.1111/cgf.13165`) is the only paper on contour trees for 2D piecewise
*polynomial* functions and has no follow-up.

---

## 7. The staged plan

Every stage ends at a **measured rate**, not at "it runs". That is this
project's standing rule and it is what caught the Solov'ev coefficients, the `τ`
sign and the inert normalisation.

### IN-A — critical points, and nothing else

**Measurable today, with no contouring at all**, which is why it is first. Find
the magnetic axis as a root of `q = 0` by barycentric Bernstein subdivision plus
a certifying degree test; audit with Σ indices = χ.

**Acceptance.** The axis position converges against the analytic Solov'ev axis
at the rate `q`'s own order predicts. The audit returns exactly 1 on every
fixture in `tests/analytic/`. And — the control that makes it mean something —
it returns 1 on a *deliberately under-resolved* mesh where spurious pairs are
likely, or the persistence threshold explains why not.

**Note the existing `ψ_ax` is a different quantity** and stays: it is "the
largest **nodal** value", `O(h^{k+1})` from the true polynomial maximum, chosen
because it is what makes the bordered Newton's constraint differentiable. IN-A's
axis is the critical point. Do not conflate them, and do not "fix" one to match
the other.

### IN-0 — the tracer

Predictor–corrector per §3.1, curvature-controlled step, cubic Hermite from `q`
per §3.2(i). Start from IN-A's axis.

**Acceptance.** On an exact Solov'ev contour, the geometric error of the traced
curve converges at `O(Δs⁴)` with Hermite against `O(Δs²)` with chords, both
measured, **with the chord version kept as the control**. Closure error over a
full circuit at machine precision, and — the property that says the corrector is
doing its job — **independent of path length**. At fixed `Δs`, the error is
independent of `h`, which is what says (b) and (c) have been separated.

### IN-1 — arc length and the metric

The trap of §3.2 made into a measurement.

**Acceptance.** Arc length of a closed contour computed two ways: from `q`
pointwise, and by differencing node positions. The first converges spectrally in
the number of points, the second at `O(Δs²)`. **If the differenced version ever
converges spectrally the comparison is empty and the test is worthless** — say
so in the failure message, as `ExtensionConvergence` does for its pinned-zero
control.

### IN-2 — flux-surface averages

`⟨r^{-2}⟩_ψ`, `V′(ψ)`, and `q(ψ)`, by periodic parametrisation and the
trapezoidal rule per §3.2(ii).

**Acceptance.** Against the closed-form Solov'ev values: spectral convergence in
the number of angular points at fixed mesh, and `k+1` in `h` at fixed angular
resolution. **Two independent controls**, because this is the number everything
downstream rests on: ray bisection from the axis (§3.4) as a second extraction,
and an implicit-quadrature rule (§3.3) as a route that never extracts a curve at
all. All three agreeing is worth more than any one being plausible.

**This is where `ROADMAP.md` item 10 becomes reachable.**

### IN-3 — the representation

Fit `R(ρ, θ)`, `z(ρ, θ)` in Fourier–Zernike per §4.1, on `ρ = √Ψ_N`.

**Acceptance.** The fit converges spectrally in mode number on an analytic
equilibrium. And the one that matters: **`dGeometry/dΨ` obtained by
differentiating the fit agrees with a finite difference of independently traced
surfaces**, to the fit's own accuracy. That is the property `MANTA-COUPLING.md`
needs and the reason a representation exists at all.

A second, cheap check with real teeth: the fit's own `|dx/dθ|` against `q`
evaluated pointwise. Two independent routes to the same metric.

### IN-4 — the `ψ`-element decision

Measure the three candidates of §4.3 against IN-2's converged `q(ψ)`, including
how each behaves as the outermost surface approaches the separatrix.

**Acceptance.** A decision, recorded with the measurement that made it. This
stage is allowed to conclude "global Zernike is sufficient for the flux range we
care about" — that is a result, not a failure.

### IN-5 — open surfaces

Chebyshev per §4.2. **Deferred with free boundary**, per §6.

### IN-6 — the output

`DRIVER-PLAN.md` §3's flux-surface `(Ψ, θ)` NetCDF grid, and whatever
`MANTA-COUPLING.md` settles on.

---

## 8. Risks, in likely-to-bite order

1. **The metric trap of §3.2.** Highest probability, lowest visibility: a
   spectral rule fed a differenced Jacobian is silently second order. IN-1 is
   the mitigation and it is a stage rather than an assertion for that reason.
2. **The separatrix.** Everything degrades approaching it — `1/|∇ψ|` in the
   integrand, the corner in the surface, the logarithm no polynomial basis
   catches. Production codes simply stop: LIUQE cuts at `Ψ_N = 0.95`, FreeGS
   extrapolates outside `[0.01, 0.99]`. **Decide meq's cut deliberately and
   record it**, rather than discovering it as a convergence failure.
3. **`Δs` chosen once and left.** The whole §2 argument rests on `Δs` being
   *tuned to `h` and `k`*. A fixed `Δs` reintroduces the cap the method exists
   to avoid — and it will look like the method failing rather than the parameter
   being wrong.
4. **The axis.** `1/|∇ψ|` diverges there, surfaces shrink to a point, and the
   corrector's `∇ψ(c−ψ)/|∇ψ|²` is 0/0. IN-A first is the mitigation; the
   innermost surfaces may still want the elliptical Hessian approximation LIUQE
   uses.
5. **A rotating equilibrium's `p`.** §6's second carve-out. If anyone asks for
   pressure contours, none of the scoping in §6 applies and the topology work
   is back.
6. **Conflating IN-A's axis with `ψ_ax`.** Two different quantities, both
   correct, `O(h^{k+1})` apart.

---

## 9. Out of scope

| | |
|---|---|
| **Topological classification** — contour trees, open/closed, saddle crossing | §6. It belongs with free boundary, where `χ(annulus) = 0` forces X-points to exist. |
| **3D / stellarator geometry** | meq is axisymmetric. DESC's basis is borrowed; its toroidal Fourier factor is not. |
| **Straight-field-line, PEST, Boozer or Hamada angles** | Stability and gyrokinetic codes want these; flux-surface *averages* do not. Build them on top of IN-3 if a consumer appears. Note the practical warning that a straight-field-line angle resolves the low-field-side midplane poorly, which is where ballooning modes live. |
| **Contouring `p` in a rotating equilibrium** | A different curve family (§6). Not hard, but not this. |
| **The inverse-equilibrium formulation** | VMEC, DESC and CHEASE make the surfaces unknowns of the *solve*. That is a different code, not a post-processing step, and it hard-codes star-shapedness into the discretisation. |

---

## 10. Two holes in the literature that meq is positioned to fill

Recorded because they are unusual, and because both surveys reached them
independently.

* **No paper treats iso-contouring a field whose gradient is an independently
  solved unknown at the same order.** The nearest is the gradient-augmented
  level set line (Nave, Rosales & Seibold, `10.1016/j.jcp.2010.01.029`), where
  the gradient is *advected* rather than co-solved. §3.2's Hermite-with-`q` is
  well-founded by analogy and is **not** something meq can cite as done.
* **Nobody extracts flux surfaces from a `P_k` finite-element `ψ` at
  `O(h^{k+1})` and demonstrates that `q(ψ)` and the flux-surface averages
  converge at that rate.** CHEASE bisects its element cubic and preserves FE
  order in the quadrature (1996); ECOM does the fully spectral version but
  cannot handle an X-point (2015). Between them is exactly this item.

meq is the only code in either survey with `q` in hand.
