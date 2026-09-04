# Solution inversion: `ψ(R, z)` to `R(Ψ, l)`, `z(Ψ, l)`

Written 2026-09-02, after a literature survey whose references are indexed in
`refs/Refs.md` under *Solution inversion*. `CLAUDE.md` is the operational record
and is authoritative on anything technical; `ROADMAP.md` is the priority order.
This file is the design for one item.

**IN-A, IN-0, IN-1, IN-2, IN-3 and IN-4 are done and green; IN-5 is deferred
with free boundary, IN-6 is open, and IN-P is under way.** So §§2–6 are the
design the code was built from and the arguments it was built on — several of
which the measurement then falsified, which is why they are still here — and §7
is what each stage found. **Section numbers are load bearing**: `src/meq` and
`tests/` cite §2, §3.2, §3.3, §3.4, §4.1, §4.3, §4.4, §5, §6 and §11 by number,
§4.4 cites §4.3, §4.3 cites §8.2, and `CLAUDE.md` cites §5 and §11. Do not
renumber without a `grep -rn`.

**The problem.** meq solves for `ψ` on a mesh. Almost everything downstream
wants the inverse: the flux surfaces themselves, as curves parametrised by a
flux label and an in-surface coordinate, and the integrals taken over them.

**What it unblocks**, and this is why it is worth doing before either of the
things that need it:

* **`ROADMAP.md` item 10, the fixed-`q(ψ)` solver.** RoPP (142) is
  `q(ψ) = V′(ψ) I(ψ) ⟨r^{-2}⟩_ψ / 4π²` — a flux-surface average, and the
  machinery `FLOW-PLAN.md` §3.3 deliberately avoided needing. **Here and in
  §3.4, §4.3 and §4.4 that `q` is the safety factor and NOT meq's flux**, which
  is the same letter for a solved unknown of the discretisation; IN-2 settles
  what the code is allowed to call each.
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
> means curvature control. Take the geometric variant and say why in the code.
>
> **BUT `Δs ∝ 1/κ` IS NOT WHAT EQUIDISTRIBUTES INTERPOLATION ERROR, AND AN
> EARLIER DRAFT OF THIS BOX SAID IT WAS.** A circular arc of turning angle
> `θ = κ Δs` departs from its cubic Hermite by `θ⁴/128` per unit radius, so a
> segment's deviation is `κ³ Δs⁴/128`. Equidistributing **that** wants
> `Δs ∝ κ^(−3/4)`. What `Δs ∝ 1/κ` equidistributes is the **turning per step**,
> which makes the deviation proportional to `Δs` and the deviation *relative to
> the segment's own length* constant. That is the robust geometric variant and
> is what IN-0 implements — but they are different rules and the difference
> should be stated rather than blurred.

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

> **AND IT ONLY HOLDS FOR ONE OF THE TWO PAIRINGS MEQ HAS. MEASURED IN IN-0,
> AND IT IS THE FINDING THAT STAGE DID NOT EXPECT.** The interpolant is built on
> tangents from the **flux** and is measured against the level set of the
> **potential**. Those are the same curve only to the extent that the two fields
> agree — and `∇ψ_h/r` agrees with `q_h` only to `O(h^k)`, because
> differentiating an L2 potential of degree `k` loses an order while `q_h` keeps
> `k+1`. So paired with `ψ_h`, the Hermite is fourth order in `Δs` **until the
> tangent tilt takes over and second order afterwards**: rates 3.809 → 1.400 →
> 1.569 down a `Δs` sweep, against a tilt of 3.5e-5 to 7.3e-5.
>
> **The control that identifies the tilt as the cause** is a third column built
> on `∇ψ_h` instead of `q_h` — the *exact* tangent of the curve being measured
> against, which therefore stays fourth order in both pairings (3.957, 3.993,
> 3.995). It is not the column to prefer: `∇ψ_h` is the right tangent for error
> (c) and the *wrong* one for error (a), where `q` is a full order closer to the
> true curve. It exists to separate the two.
>
> Paired with `ψ*` and `q*` the tilt is a full order smaller — 3.1e-7 to 5.9e-7
> — and fourth order survives the whole sweep, 3.960 / 3.996 / 3.812. **So the
> post-processed pairing is not merely more accurate; it is what makes this
> item's claim about `q` true at all on this discretisation.** IN-0's default is
> `Potential::PostProcessed` for that reason and not for the order of `ψ*`
> alone.
>
> **The usual reason given for the pairing is wrong, and was checked rather than
> repeated.** It is not that the local post-processing is built so that
> `∇ψ*` matches `r q*`: MFEM's reconstruction drives Stenberg's local problem
> with the reconstructed **total** flux `q̂_h` in `RT_k` — the normally
> continuous field the constraint equation projects onto — so `∇ψ*` and `r q*`
> are different objects. What is true, and measured, is that they agree an order
> better than the raw pair does. Right answer, wrong mechanism, and the
> distinction matters because the wrong mechanism would predict exactness.

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

### 4.1 Closed surfaces: the Zernike disc basis

**Called "Fourier–Zernike" in the first draft of this file, and that name is
dropped.** It is DESC's, where the Fourier factor is the *toroidal* one; in 2D
axisymmetric that factor is gone, and what is left has no separate Fourier part
at all — the `θ`-dependence **is** Zernike's own angular factor. A reader given
the compound name goes looking for a second index that does not exist.

> **THE ANGLE IS NOT FREE, AND THE OBVIOUS CHOICE IS INADMISSIBLE. MEASURED IN
> IN-3, AND IT IS THE MOST CONSEQUENTIAL FINDING OF THE ITEM.**
>
> Labelling surfaces by the **geometric poloidal angle about the axis** — the
> obvious choice, and what `AngleParametrisation` produces — **makes the disc map
> non-smooth at the axis for any non-circular surface**, so no basis that is
> smooth there can converge against it.
>
> The argument is exact. A function smooth at the origin has precisely **one**
> angular harmonic multiplying `ρ¹`. Take nested ellipses of semi-axes `a`, `b`
> labelled by geometric angle: with `u = ρ cos θ`, `v = ρ sin θ`, the map is
>
> ```
> x( u, v ) = a b · u · √(u² + v²) / √(b² u² + a² v²)
> ```
>
> whose trailing factor is **homogeneous of degree zero** — a function of
> direction alone, with no limit at the origin unless `a = b`. Expanded, the `ρ¹`
> coefficient carries `cos θ, cos 3θ, cos 5θ, …`, and **every harmonic past the
> first is one the Zernike index constraint excludes**, precisely because it is
> not smooth there. The parametrisation puts content exactly where the basis
> refuses to look.
>
> Measured on nested ellipses, where the answer is known exactly, the geometric
> column decays like `L^{-1.2}` and **never converges**: 1.53e-01 at `L = 2` and
> still 1.07e-02 at `L = 20`, against 1.8e-14 for the same points relabelled.
>
> **The repair needs no field.** Near the axis every equilibrium's surfaces are
> ellipses, so a three-parameter fit of a quadratic form to the *innermost traced
> surface* recovers the tilt and flattening, and the relabelling is then an exact
> reparametrisation of the circle. On `nstx()` it recovers short/long = 0.4943
> from samples against 0.4845 from the Hessian, and is worth **45× to 660×** in
> worst fit error. `meq::relabelByAxisShape` is it.

**And `l` is overloaded in this document.** The title and §2 use `l` for the
in-surface coordinate; the literature uses `l` for the Zernike radial degree.
Below, and in `src/meq/Zernike.hpp`, **`l` is the degree and `θ` is the angle**,
which is the universal convention. `R(Ψ, l)` in the title means `R(Ψ, θ)`.

Nested closed surfaces filling in to a single maximum are a **disc**,
topologically, and the map `(Ψ, θ) → (R, z)` is a map from a disc. The natural
basis is the one DESC uses for the same reason
(<https://desc-docs.readthedocs.io/en/stable/theory_general.html>):

```
R( ρ, θ ) = Σ_{l,m} R_lm Z_l^m( ρ, θ )

Z_l^m = R_l^m( ρ ) cos( m θ )    for m >= 0
        R_l^m( ρ ) sin( |m| θ )  for m <  0        -- the sign of m selects which,
                                                      and without saying so a
                                                      coefficient vector is not
                                                      well defined
```

with `l ≥ |m|` and `l − |m|` **even**. `R_l^m` then contains only
`ρ^|m|, ρ^{|m|+2}, … ρ^l`, and that parity-and-minimum-power constraint is what
makes the magnetic axis a regular point instead of a coordinate singularity to
special-case.

**Implementation: `R_l^m` is a Jacobi polynomial under a coordinate change**,

```
R_l^m( ρ ) = (−1)^{(l−|m|)/2} ρ^|m| P_{(l−|m|)/2}^{(|m|, 0)}( 1 − 2 ρ² )
```

so `boost::math::jacobi` and `jacobi_prime` give the value and an **exact**
derivative. **The argument is `1 − 2ρ²` and not `2ρ² − 1`**, and the leading
sign is not optional: measured over 136 modes to degree 30, the flipped
argument fails `R_l^m(1) = 1` on **113** of them and the dropped sign on **64** —
but *both* reproduce every pure power `R_l^l` exactly, which is why `R_3^1 =
3ρ³ − 2ρ` is the cheapest case that separates them. The explicit factorial sum
that every reference prints is **not** an acceptable alternative: it disagrees
with the Jacobi route by 4.5e-08 at `l = 30` and 1.3e-04 at `l = 40`.

> **A NAIVE TENSOR PRODUCT ADMITS MODES THAT ARE C¹ AT THE ORIGIN AND FAIL AT
> THE SECOND DERIVATIVE — WHICH IS WHY ADMITTING ONE WOULD BE SILENT.** The
> first draft of this file said such modes are "not analytic at `ρ = 0`". True,
> but weaker than it sounds and misleading about how it would be caught:
> `ρ² cos θ = x√(x²+y²)` has a *continuous gradient* at the origin, so a C¹
> check passes it. The sharp statement is that **every admissible Zernike mode
> is a polynomial in `(x, y)`** and the inadmissible ones are not. Measured, by
> a one-sided second difference from each side of the origin: admissible modes
> give 0 or shrink at `O(h)`; `ρ² cos θ` sits at `4|cos θ|` at every step size.
> (A *centred* second difference is useless for this — the defect is odd and
> cancels itself out of the measurement.)

> **THE RADIAL COORDINATE MUST BE `ρ = √Ψ_N`, NOT `Ψ_N`.** DESC's theory page
> says `ρ` is "the square root of the normalized toroidal flux, which is
> proportional to the minor radius", and the reason bites us identically:
> **`ψ` has a quadratic maximum at the axis**, so `ψ − ψ_ax ∝ (distance)²`.
> Parametrise by `Ψ` directly and `R(Ψ, l)` acquires a square-root branch point
> at the axis; every basis then converges algebraically, worst near the axis,
> and no convergence table says why. The Zernike degree counts powers of `√Ψ_N`.

> **AND WHAT BUYS GEOMETRIC CONVERGENCE IS SMOOTHNESS IN `(x, y)`, NOT
> SMOOTHNESS IN `ρ`.** Measured through the same quadrature: a function with a
> pole just outside the disc decays geometrically, `0.42` per degree, reaching
> `6.2e-09` truncation error at degree 24. The control is `f = √(x²+y²)`, which
> **is** `ρ` and is therefore perfectly smooth as a function of `ρ` — and is
> only C⁰ at the Cartesian centre. It decays **algebraically**, at `l^{-1.92}`,
> and is seven orders worse at degree 24. So "the surfaces are smooth in the
> flux label" is *not* the hypothesis this basis needs; smoothness of the map
> into the plane is.

**Ordering matters and is not cosmetic.** `src/meq/Zernike.hpp` uses ANSI/OSA
order because it is **prefix stable**: truncating to a lower degree is
truncating the coefficient vector, with no re-indexing. A convergence study in
mode number — which is IN-3's acceptance — needs exactly that, and it is not
free from an arbitrary ordering.

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

### 4.3 The band between `Γ_h` and `Γ`, and why it reaches this item

**On the curved path the outermost flux surfaces are not in the mesh.** `Ω_h` is
the union of background elements lying *inside* `Γ`, so `Γ_h` is inscribed and
there is a band `O(h)` wide that is inside the plasma and outside the mesh.
`CLAUDE.md` records the measurement: on `examples/miller-curved.toml` about one
grid node in ten is in it.

**This is not a corner case for this item, it is the interesting part.** The
outermost closed surface is `ψ = 0`, which *is* `Γ`. So the surfaces the band
affects are exactly those with `Ψ_N → 1` — which is where `q(ψ)` and the
flux-surface averages are most wanted and least forgiving, and where every
production code already struggles for unrelated reasons (§8.2).

**The tracer must therefore accept a field defined outside the mesh, and must
know when it is using one.** That is a requirement on IN-0, not a redesign: the
corrector calls something for `ψ` and something for `∇ψ`, and in the band those
calls are answered by an extension. What it must not do is what the prior art
does — `v0-legacy:FluxSurfaces.cpp` prints *"Terminating because curve left
domain"* and returns a partial curve, which for an outer surface silently
returns an arc labelled as a closed contour.

**~~Which extension is a real question and it should be measured, not argued.~~
MEASURED, 2026-09-02: `BandExtension::TransferLift`, and it is not close.**

The decisive test traces **`Γ` itself**, which lies entirely in the band at
every mesh — `ψ_h` is strictly negative inside `Ω_h`, so every point of
`{ψ = 0}` is answered by the extension, and the exact answer is the curve `D_h`
was cut from. Worst distance from the exact `Γ`, over `n = 16/32/64`:

| `k` | flux Taylor, rate | transfer lift, rate | lift closer at `n = 64` by |
|---|---|---|---|
| 1 | 1.797 | **2.298** | **40×** |
| 2 | 2.138 | **3.995** | **1,610×** |
| 3 | 2.138 | **4.848** | **84,695×** |

**The Taylor step is second order at every `k` and the reason is structural**:
its remainder is `O(h²)` over a band of width `O(h)` however good `q` is, so it
cannot improve with the polynomial degree and it does not — 2.138 at `k = 2` and
at `k = 3` alike. The lift is the error of `q` *integrated along a path of
length `O(h)`*, which is `k+2` at `k = 2, 3`. At `k = 3` the Taylor step is five
orders of magnitude worse.

**The Taylor step is kept as the control and not as a fallback**, and the tests
say in their failure messages that if it ever matches the lift the comparison is
empty.

> **AND §4.3'S THIRD ARGUMENT FOR `ψ`-ELEMENTS IS WEAKENED BY THIS, WHICH
> MATTERS AT IN-4.** The claim below is that a global fit's accuracy is set by
> its worst region, so surfaces known only to the extension's order degrade the
> representation everywhere. **With the lift, the band is not the worst
> region** — it converges at `k+2` where the interior population of the same
> contour reads `k+1` to `k+2`. That argument was written when the extension was
> assumed to cost an order. It does not. **Two of the three arguments for
> elements survive; this one should not be leaned on.**

Three candidates were considered, and meq already owned machinery for two:

| | | |
|---|---|---|
| **the flux Taylor step** | what `GridSampler::samplePotentialWithFlux()` already does: `ψ(p) ≈ ψ(x₀) + r₀ q(x₀)·(p − x₀)` from the foot `x₀` on `Γ_h` | **nothing is ever evaluated outside an element**, which is the property that made it right for the `.nc`. But `∇ψ` in the band is then **frozen at the foot**, so the extended field is affine and **contours in the band are straight lines** — the curvature is gone, whatever `k` is |
| **the transfer-path lift** — **CHOSEN** | the extension technique's own construction. **NOT through `PathLiftCoefficient`, which this file named and which cannot do it**: that class `dynamic_cast`s its `ElementTransformation` to `FaceElementTransformations` and lifts from *that face's own* integration point, answering "what is `φ_h` on `Γ_h`" — η₅'s question, not this one. The usable primitive is one level down and public, **`mfem::PathIntegral( Cu, x, xbar, line_ir )`, which takes arbitrary endpoints**, with `mfem::ElementExtension` supplying `E_h(q_h)`. MFEM's own comment on it — it must return `p(x) − p(a(x))` "whatever the path" — is the licence | **method-native**: it is how the solver itself relates `ψ` on `Γ` to `ψ` on `Γ_h`, by integrating `r q` along the path. So the answer to "is it usable at an arbitrary band point" is **yes, but not through the class named here first** |
| **the known boundary condition** | `ψ = 0` exactly on `Γ`, and `Γ` is known analytically from `BoundaryShape` | the outermost contour is therefore known **exactly**, for free, with no tracing at all — a fact worth exploiting whatever else is chosen |

> **A WARNING FROM THE TREE, BECAUSE THE OBVIOUS COMBINATION HAS BEEN TRIED AND
> FAILED.** Blending an extrapolated value toward the known `ψ = 0` on `Γ` does
> **not** rescue a bad extension: `CLAUDE.md` records that `(1−t)·v` scales a
> positive value down and never changes its sign, so all seventeen offending
> nodes survived it. *"The error was in **where the field was evaluated**, not in
> how it was weighted."* The purpose here is different — parametrising a curve
> rather than sampling a grid — but the lesson transfers: fix the evaluation,
> and use the exact boundary as an **anchor** rather than as a correction.

**AND IT SHARPENS §4.4, WHICH IS THE PART THAT ACTUALLY BEARS ON THE DESIGN.**
A single global Zernike expansion in `ρ` is a **global** fit: its accuracy is set
by its worst region, so surfaces known only to the extension's order degrade the
representation *everywhere*, not just near `ρ = 1`. Elements in `ψ` localise
that damage — and they give the outermost element a natural place to begin, at
the flux label where surfaces start crossing into the band. **That is a third
argument for elements over one global expansion, independent of the two in
§4.4, and it is the one that arrived last.**

**Report it, the way the `.nc` already does.** That file carries a `byte
extrapolated( Z, R )` beside `inside` precisely because a band node holds real
data and is therefore indistinguishable from a solved one without being told.
The same is true of a surface: **every flux-surface quantity computed on a
surface that crosses the band must be flagged as such**, so that a consumer can
drop it before computing an error norm or differencing two runs. A count is not
enough — `CLAUDE.md` records that being the other half of the same defect.

> **THREE TRAPS MET IN THE MEASURING, AND THE FIRST NEARLY PRODUCED A WRONG
> HEADLINE.**
>
> **A contour at fixed `Ψ_N` cannot measure the extension's order.** As `h`
> falls, `Γ_h` climbs towards `Γ`, so a contour sitting a fixed distance inside
> `Γ` has its band excursion shrink *faster* than `h` — measured, `deep/h` goes
> 0.92 → 0.66 → 0.26 and the band population goes 136 points → 9. **Both**
> columns then converge faster than the extension they are built on: the Taylor
> step reads 3.75 at `k = 3` against its true 2.1. Tracing `Γ` itself is what
> fixes it, because `Γ` is in the band at every mesh by construction. The
> fixed-`Ψ_N` test is kept, asserting on the *margin* between the two extensions
> rather than on a rate.
>
> **The nearest face of `Γ_h` is not the face you are outside of.** A staircase
> `D_h` cut from a diagonally split Cartesian mesh **pinches** — two triangles
> meeting at a single vertex — so both lobes' faces are equidistant from a point
> just outside the pinch and the tie goes to loop order. Half the time that
> picks the lobe whose outward normal points the other way, the outwardness test
> fails, and a genuine band point is refused: measured at `k = 2, n = 64` the
> trace stopped after 85 points of 320. Apply the outward test **first, as a
> filter**, and take the nearest among the survivors.
>
> **`ψ_h` is not strictly negative inside `Ω_h` to machine precision**, so a
> handful of points of `{ψ = 0}` land back inside — one of 322 at `k = 1,
> n = 64`. That is the discretisation creeping above its own imposed datum near
> `Γ_h`, not a tracer fault, and the acceptance asserts ≥ 98% band rather than
> 100%.
>
> **And the order in question is in `h`, not in `Δs`.** An earlier draft asked
> for the band rate against `Δs` and expected it "worse than the `O(Δs⁴)` the
> fitted path reaches". The band error is §2's error **(a)** — the field — and
> (a) and (c) are decoupled by §2's own argument. `Δs` does not enter.

**None of this arises on the fitted path**, where `Γ` is the mesh boundary and
there is no band. So the fitted path is the right place to establish every rate
in §7 before the curved path is turned on, and IN-0's acceptance says so.
### 4.4 The `ψ`-varying element — **ANSWERED BY IN-4, AND THE QUESTION WAS WRONG**

**No `ψ`-element was built and none is needed.** This section asked how to buy a
parametrisation that varies with the surface; IN-4 measured that **solving for
the angle directly is cheaper than every way of buying it**. `meq::gaugeFreeFit()`
requires each disc node only to *land on the right surface* rather than to sit at
a prescribed angle. The section is kept because the argument that made the
question look necessary is sound, and because two of its four constraints are
still live.

**Why a fixed relabelling cannot work, which is the part that was right.** A
relabelling of the poloidal angle that does **not** depend on the flux label has
one function's worth of freedom, so it can straighten exactly **one** order:
§4.1's axis-ellipse relabelling fixes the `ρ¹` harmonics and `ρ²`, `ρ³`, … keep
whatever harmonics the shaping puts there. Measured: with the innermost surface
at `Ψ_N = 0.10` the envelope decays geometrically, and pulled in to `Ψ_N = 0.02`
**an algebraic tail appears that no degree removes.** A parametrisation smooth to
all orders must therefore vary with the surface — and that is a structural
difference between meq and the code this basis was borrowed from. **DESC *solves*
for its `R`, `z` coefficients**, so its parametrisation is part of the unknown and
it gets this for free. **meq fits after the fact**, so it has to be bought.

**The three ways to buy it were all the same basis, and the DESC papers say so.**
Read 2026-09-03, `refs/DESC-Dudt.pdf` and `refs/DESC-Panici.pdf`:

* A global expansion with one relabelling, accepting the tail and bounding the
  flux range away from the axis; **elements in `ψ`**, each with its own
  relabelling; or **a per-surface shape expanded in `ρ`**, which is MXH's
  factorisation.
* **The regularity condition has a name and a proof.** What IN-3 derived
  independently — `a_m(ρ) = ρ^m (a_{m,0} + a_{m,2}ρ² + …)` — is **Lewis &
  Bellan**, *J. Math. Phys.* **31**, 2592–2596 (1990), which Dudt §III.B cites
  for exactly this purpose. So the third option does not escape that condition,
  it **relocates** it onto the shape coefficients' `ρ`-dependence, where it must
  be imposed by hand. **With regularity imposed the three are not three bases;
  they are the same basis, and the only real variable is the angle.**
* **DESC has the same axis problem and says so.** Dudt: DESC "is restricted to
  operating in straight field-line coordinates, which **may appear to be a
  disadvantage compared to the optimal poloidal angle of the VMEC
  formulation**", and, measured, "**more modes are required in the core,
  resulting in increased error near the magnetic axis for a given resolution**".
  That is IN-3's finding, in the code this basis was borrowed from.
* **And the later paper shows the mechanism that fixes it.** Panici expands in a
  *general* computational poloidal angle with `λ` carrying the map to
  straight-field-line, and reports, with a spectral-width figure: "**DESC, while
  not explicitly enforcing any poloidal angle constraints, ends up finding an
  optimal representation through the course of the optimization procedure.**"
  **That is the whole difference between fitting and solving**, measured by them
  rather than argued by us — and it is what pointed at IN-4's answer: give the
  angle its freedom back and let a Gauss–Newton find it.

**Two of the four constraints this section listed are still live, and one of
them is IN-5's whole problem.**

1. ~~**Axis regularity.**~~ A global Zernike in `ρ` gives it free, and the
   gauge-free fit holds to `Ψ_N = 0.005` with no inner limit tried costing more
   than a factor of 1.5.
2. **An element boundary AT the separatrix.** A global expansion in `ρ` assumes
   smoothness all the way out, and that dies at the separatrix: the surface
   develops a corner at the X-point and flux-surface quantities diverge
   logarithmically. **No polynomial basis in `Ψ` converges against a logarithm.**
   Untouched by IN-4, and it is why IN-5 is deferred with free boundary rather
   than merely unstarted.
3. ~~**C¹ in `Ψ` at minimum**, because `MANTA-COUPLING.md` needs
   `dGeometry_dpsi`.~~ A global expansion is `C^∞`.
4. ~~**Localising the band's damage**, §4.3.~~ Removed by the transfer lift,
   which put the band at `k+2` rather than at the flux Taylor step's flat second
   order.

> **`∂R/∂Ψ` IS UNBOUNDED AT THE AXIS, AND THIS FILE ORIGINALLY READ AS IF IT
> WERE NOT.** Found while implementing `Zernike`, 2026-09-02, and it is a fact
> about the problem rather than about the basis.
>
> The chain rule from `ρ = √Ψ_N` carries a factor `1/(2ρ)`, so `∂/∂Ψ_N` at fixed
> `θ` diverges at `ρ = 0` for any mode whose lowest radial power is `ρ¹` — that
> is, `|m| = 1`. And `|m| = 1` is not incidental: **for a surface
> `R = R₀ + a cos θ`, the `m = 1` coefficient of `R` IS the minor radius `a`**,
> which grows like `√Ψ`. So `∂R/∂Ψ ∝ 1/(2√Ψ_N) → ∞` at the axis in **every**
> equilibrium, universally, not as a property of any configuration.
>
> (The implementing agent attributed this to the Shafranov shift. It is not
> that — a rigid shift is the `m = 0` coefficient and is well behaved. The
> divergence is the minor-radius growth itself, which makes it more fundamental
> rather than less.)
>
> **The distinction that saves the downstream consumers.** Flux-surface
> *averages* are unaffected: `V ∝ Ψ` near the axis, so `V′(Ψ)` tends to a
> constant and IN-2's quantities are fine. What diverges is the derivative of
> the surface **position**.
>
> **AND THE COUPLING'S SIDE OF IT IS ANSWERED, 2026-09-02, BY THE MANTA
> AUTHOR.** Two mechanisms, and the first is the real one:
>
> * **MaNTA enforces that the fluxes vanish on axis**, which they must — there
>   is no area for a flux to cross at the axis, so it is a physical boundary
>   condition rather than a numerical guard. A metric derivative that diverges
>   like `Ψ^{-1/2}` multiplied by a flux that vanishes gives a finite product.
>   The singularity is removable **in the coupled system**, which is the level
>   at which it has to be removable.
> * **And in the worst case MaNTA can guarantee that none of its nodes lies at
>   `Ψ = 0` exactly**, so nothing is ever evaluated at the singular point.
>
> **What that leaves is a conditioning question, not a correctness one, and it
> is worth measuring rather than assuming.** `1/(2√Ψ_N)` at the innermost node
> is finite but grows as that node approaches the axis, so the magnitude of the
> geometry block of MaNTA's Jacobian scales like `Ψ₁^{-1/2}`. Whether that is
> comfortable depends on where MaNTA puts its first node, which is a
> configuration choice rather than a property of either code. **Measure it once
> the coupling exists; do not design around it now.**
>
> **This is not a defect being worked around.** The minor radius really does
> grow like `√Ψ`, so any code working in `Ψ` near the axis meets it, and the
> standard remedy across the transport community is exactly the pair above:
> a radial coordinate that is not the flux, plus the axis boundary condition.
> Recording it because a reader meeting `1/(2ρ)` in `Zernike.hpp` will
> reasonably wonder whether it is a mistake.
>
> **Not written into `MANTA-COUPLING.md`**, which is a *received* document —
> `TODO` records the rule, the same one that makes `../mfem-hdg-dev`
> receive-only. It belongs here.
>
> `src/meq/Zernike.hpp` supplies both: `fluxDerivativeFromRadial()` carries the
> `1/(2ρ)` factor with the divergence documented, and `radialDerivative()` stays
> finite. Prefer the second wherever a caller can work in `ρ`.

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
* What excludes an interior extremum of ONE sign is the **maximum principle**,
  and *which* sign depends on the sign of `F`. `1/r` is bounded away from zero
  and infinity since `r > 0`, so the operator is uniformly elliptic in
  divergence form either way:
  * `F ≥ 0` makes `ψ` a **supersolution**, its minimum sits on `Γ`, and
    `#min = 0`.
  * `F ≤ 0` makes `ψ` a **subsolution**, its maximum sits on `Γ`, and
    `#max = 0`.
* Either way one interior extremum gives `#saddle = 0`: **zero interior
  saddles**.

> **AND MEQ'S SOLOV'EV FIXTURES ARE THE SECOND CASE, WHICH AN EARLIER DRAFT OF
> THIS SECTION DID NOT SAY.** It wrote only the `F ≥ 0` half. `F = −((1−A)r² + A)`
> is single-signed **negative**, so `ψ` is a subsolution, its maximum is on `Γ`,
> and **the magnetic axis is an interior MINIMUM** — Hessian determinant +0.693
> and trace +2.121 on `nstx`, measured. The conclusion above is untouched,
> because in 2D a maximum and a minimum both have index +1 and the sign never
> reaches the count. **But the practical instruction derived from it was wrong**:
> "seed the axis search from the largest nodal value" finds a **corner of the
> benchmark rectangle** on these fixtures. `CriticalPointFinder::findAxis()`
> seeds from both nodal extremes and refuses rather than guesses if both are
> interior extrema; `AxisSense` is for a caller who knows which they want.

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

Every stage ended at a **measured rate**, not at "it runs". **IN-A, IN-0, IN-1,
IN-2, IN-3 and IN-4 are done and green**; IN-5 is deferred, IN-6 is open, IN-P is
under way. `CLAUDE.md`'s *Solution inversion* carries the measurements; what is
kept below per stage is where the code lives, what the stage **found** that this
plan did not predict, and the few numbers that are recorded nowhere else.

### IN-A — critical points — **DONE, 2026-09-02**

The magnetic axis as a root of `q = 0` by barycentric Bernstein subdivision plus
a certifying degree test, audited with Σ indices = χ.
`src/meq/CriticalPoints.{hpp,cpp}`, `tests/convergence/CriticalPointConvergence.cpp`.

**The sharpest assertion was not in the brief**: the ratio of the position error
to `|q_h − q|` at the *exact* axis is **0.77 to 3.37** against **0.77 to 3.27**
predicted by linearising `q` about its root — a direct statement that the root
finder adds nothing to the error of the field it is rooting.

> **A DEGREE IS A SUM OF INDICES AND THE SUITE DEMONSTRATES IT RATHER THAN
> ASSERTING IT.** A box round `iterExample2`'s axis *and* its X-point reads
> **winding 0 with two critical points inside**, the saddle located to 4.5e-6 of
> the published X-point. And the Poincaré–Hopf hypothesis is **transversality of
> `q·n`**, not that the boundary is a level set: on the standard rectangle
> `min |q·n|/|q| = 0.15` with one sign, so `winding == χ == 1` *is* a theorem
> there; on a box reaching past the X-point it reads 0.00 with a sign change and
> is not.
>
> **AND A DISCONTINUOUS `q_h` CAN CARRY BOUNDARY DEGREE 1 WITH NO ZERO IN ANY
> ELEMENT.** At `h = 0.4, k = 1` the audit reads 1 and an element-by-element
> search finds nothing: each element's polynomial puts its zero just inside a
> neighbour's territory, and with a face jump of `O(h^{k+1})` against an element
> of size `h` there is a window where the zero belongs to neither.
> **Poincaré–Hopf is a theorem about continuous fields.** The window closes with
> refinement and the test pins it to that one coarse mesh, so a *finer* mesh
> losing the axis fails rather than passing quietly. Where the audit and the
> search disagree, believe the audit.

**Note the existing `ψ_ax` is a different quantity** and stays: it is "the
largest **nodal** value", `O(h^{k+1})` from the true polynomial maximum, chosen
because it is what makes the bordered Newton's constraint differentiable. IN-A's
axis is the critical point. Do not conflate them, and do not "fix" one to match
the other.

### IN-0 — the tracer — **DONE, both halves, 2026-09-02**

Predictor–corrector per §3.1, curvature-controlled step, cubic Hermite from `q`
per §3.2(i), started from IN-A's axis. `src/meq/FluxSurfaces.{hpp,cpp}`,
`tests/convergence/FluxSurfaceConvergence.cpp`. The band is
`BandExtension::TransferLift`, chosen by measuring it against the flux Taylor
step — §4.3.

**THE DEFECT THIS STAGE MET AND THE PLAN DID NOT PREDICT.** `{ψ_h = c}` is a
union of per-element arcs offset by the jump, so a point landing within
`jump/|∇ψ|` of a face is on **neither**: the Newton step computed in element A
pushes it into B, B's pushes it back, and the residual alternates without ever
meeting a tolerance tighter than the jump. It needs a point to land in a band
about 2e-5 wide on a contour of length 1.7, it gets **commoner as `Δs` falls**
simply because more points are placed, and left alone **it ends the trace** — it
ended several before it was diagnosed. The corrector now keeps its best iterate,
accepts after four non-improving steps, refuses to travel more than one predictor
step, and reports `stalledCorrections` against `correctorTarget`. A segment with
such an endpoint is excluded from the representation-error measurement, because
one endpoint accepted at the jump level poisons a segment's interpolation error
by a factor of a hundred on a coarse mesh.

**And the plan's "closure error at machine precision" was wrong twice over.** It
cannot be at the corrector tolerance even normal to the curve: both the returning
point and the start are on the level set to 1e-13, but they are separated *along*
the curve by the final step's tangential offset `g_t`, and an arc departs from its
own tangent line by `κ g_t²/2` — measured 3.573e-10 against a bound of 7.2e-10,
which is the discriminating statement, since a drifting tracer's normal error
would grow with path length and be unrelated to `g_t`. Nor is "1, 5 and 10
circuits, flat" right for the closure: the 1-circuit value is 36× smaller because
the step controller starts cold and that trace's shortened final step lands
better. **Comparing 1 against 10 compares two final steps, not two path lengths.**
The residual is what is flat.

**The band mask is asserted point by point against `mfem::Mesh::FindPoints`** —
which knows nothing of the band machinery — in both directions, because
`CLAUDE.md` records that a *count* rather than a mask was the other half of a real
defect in the `.nc`.

### IN-1 — arc length and the metric — **DONE, 2026-09-02**

§3.2's trap made into a measurement, and it measures at five orders. The metric
length against the Hermite arc length agrees to **2.5e-08**; transversality
`min|u × t|` reads 0.844 / 0.823 / 0.805 and a fit about a point *outside* the
surface is **refused**.

> **"SPECTRAL IN `N`" IS NOT ATTAINABLE ON A DISCRETE CONTOUR, AND THAT IS NOT A
> DEFECT.** `ψ_h` jumps across faces, so `ρ(θ)` is piecewise analytic with jumps
> and no quadrature is geometric on it. The column plunges — 7.24 then 14.8 — and
> then **floors at about 1.2e-9**, which is where the DG jump of `ψ*` converts to
> a distance (6.80e-10). The control that says the floor is the **field** and not
> the **rule** is the identical rule run on the *analytic* contour, which reaches
> **3.775e-15** — and independently, on the same fixture,
> `tests/analytic/FluxSurfaceReference.hpp` reaches 3.11e-15. Two
> implementations, one conclusion.

> **STILL OPEN, SMALL, AND FOUND BY IN-2 USING IT: `fitByAngle()` THROWS WHERE
> THE CORRECTOR WOULD ACCEPT.** The ray Newton demands `|ψ_h − c| ≤ tol × scale`,
> and on a discontinuous field that is sometimes **unattainable** — a ray crossing
> a face where `c` falls inside the jump has no point on it with `ψ_h = c` at all.
> This is the same phenomenon the tracer's own corrector already handles by
> keeping its best iterate; the ray Newton does not, and throws. Measured, at
> `k = 1` on the raw pairing it fails at 1e-12 on **every** mesh from `n = 12` to
> 32, and the failure probability rises with the angle count — the same contour
> that fits at 256 angles can fail at 4096. **Give the ray Newton the corrector's
> best-iterate acceptance, with a `stalledRays` count beside `worstResidual`**;
> the tolerance ladder IN-2's test currently carries then goes away. Two smaller
> ones: `AngleParametrisation` evaluates `q` at every node to build `ρ′` and
> **discards it**, forcing a consumer to re-sample; and `sampleAt()` should report
> `extended`, without which a caller flagging a quadrature point has to mark whole
> segments conservatively.

### IN-2 — flux-surface averages — **DONE on the fitted path, 2026-09-03**

`src/meq/SurfaceAverage.{hpp,cpp}`, `tests/convergence/SurfaceAverageConvergence.cpp`.
One primitive over a callable integrand in `(R, z, ψ, q)`, with every named
quantity a one-line wrapper — the shape argued for below, because
`MANTA-COUPLING.md` says the slot list "is negotiated with the transport physics
case, not fixed by MaNTA" and its illustrative set is known to be under revision.
**Do not build against the enumeration.** Conventions:
`V′ = ∮ 2πR dl/|∇ψ|`, `⟨X⟩ = (1/V′) ∮ 2πR X dl/|∇ψ|`.

> **A NAME COLLISION, SETTLED.** In this project `q` is the **flux**,
> `q = (1/r)∇̄ψ`. The safety factor is also universally written `q`. **In code the
> safety factor is `safetyFactor` and never `q`**, and in this file it is named in
> words. A reader who meets `q` in `src/meq` is entitled to assume the flux, and
> one silent conflation would be very hard to see afterwards.

> **AND THE CONSUMER'S JACOBIAN HAS THIS PROJECT'S OWN FAILURE MODE.**
> `MANTA-COUPLING.md` §6: MaNTA never assembles its Jacobian, so a missing
> surface-motion term "does not produce a wrong answer, only slow Newton
> convergence", and several defects survived a passing suite for months on it.
> That is *A wrong Jacobian is invisible to a convergence table* arriving from the
> consumer's side, and the remedy is the one meq already uses on `dFdPsi`:
> **finite-difference `dGeometry_dpsi` against `Geometry` directly**, as a unit
> test, and never take convergence of the coupled solve as evidence.

**MEQ'S SOLOV'EV FIXTURES HAVE NO CLOSED-FORM FLUX-SURFACE AVERAGES, AND AN
EARLIER DRAFT ASKED FOR THEM.** `ψ` is closed form; `⟨R^{-2}⟩`, `V′` and the
safety factor are integrals over a *contour* of it, and the contours of the
Cerfon–Freidberg twelve-coefficient family are not curves whose arc length is
elementary. **Be precise about the scope of that**: it is a statement about
*these* fixtures, not about Solov'ev equilibria in general — a deliberately
constructed one whose surfaces have elementary arc length would give a genuine
closed form, at the cost of a new fixture with its own transcription to check.

What replaces it is a **converged reference on the exact field**,
`tests/analytic/FluxSurfaceReference.hpp`: a safeguarded root solve along rays
against analytic `psi` and `gradPsi`, with the periodic trapezoidal rule. It never
touches `ψ_h`, so it separates error (a) from (b) and (c) exactly as §2 requires.
**That is a reference *value*, not a closed form**, and the distinction should be
stated wherever the number is printed. Its values are asserted in
`SurfaceAverageConvergence.cpp` and tabulated in the fixture's header.

**AND THERE IS AN EXACT IDENTITY THAT NEEDS NO REFERENCE VALUE AT ALL.** The
flux-surface average of the Grad–Shafranov equation itself. With
`⟨∇·Γ⟩ = (1/V′) ∂_ψ( V′ ⟨Γ·∇ψ⟩ )` and `Γ = ∇ψ/R²`, using
`∇·( ∇ψ/R² ) = Δ*ψ/R²` and `Δ*ψ = −F`:

```
(1/V′) d/dψ ( V′ ⟨ |∇ψ|² / R² ⟩ )  =  ⟨ Δ*ψ / R² ⟩  =  −⟨ F / R² ⟩
```

**Write the right-hand side as `−⟨F/R²⟩` and not as `−μ₀p′ − g g′⟨R^{-2}⟩`.** The
two are equal, but the first uses **the `F` the solver is actually fed** and so
applies to every fixture in `tests/analytic/` rather than to Solov'ev alone. That
is the discipline `deltaStarFD()` follows in checking the twelve-term
transcription, and for the same reason: **an independent quantity is the only
thing that catches a misread formula, and a right-hand side re-derived by hand is
not independent of the hand that derived it.** The convention is part of the
identity — it is stated for `V′ = ∮ 2πR dl/|∇ψ|`, and a per-unit-length `V′`
changes it. Verified on `ψ = (R−2)² + z²`, whose level sets are circles, to
**2.7e-13**.

> **QUOTE A DIFFERENCED RESIDUAL WITH ITS STEP OR NOT AT ALL.** An earlier draft
> gave "1.3e-13, 6.3e-13, 4.2e-12" as though those were properties of the
> routine. They are properties of the **step**: swept on the exact field the
> residual runs 9.6e-08 at 5% of `|ψ_ax|` down to 2.3e-11 at 0.6%, bottoming near
> 1e-12 around 0.15%.
>
> **And the step is NOT monotone on a discrete field.** At `k = 2, n = 96` a step
> of 2% reads 1.5e-08 where 1% reads 3.4e-07 — **the smaller step is 20× worse**,
> because the difference divides the surfaces' own DG-jump noise by the step.
> There is an optimum and it must be found rather than assumed.

**Two expectations died here**, and both are in `CLAUDE.md`: `ψ*` does **not**
buy `k+2` in an average, because the weight carries `|∇ψ| = r|q|` and there is no
`k+2` flux to divide by; and an average does **not** escape the metric trap, a
ratio cancelling a constant (about 40×) and **nothing in the order**.

**A free fourth check**: re-parametrising the *same* contour about a displaced ray
origin moves `V′` by **3.6e-10** at an offset of 0.05 and 3.9e-11 at 0.15, with
transversality falling 0.658 → 0.438. The average is a property of the surface and
not of the chart, and now that is measured.

**Two of the three cross-checks were built and the third was not.** Ray bisection
from the axis (§3.4) is kept as the second extraction; **§3.3's implicit
quadrature — a rule on the level set with no curve extracted at all — is the
missing third leg** and is deliberately not built.

**The fixture needed its own box.** `standardBox()` cannot hold these surfaces —
`Ψ_N = 0.25` on `nstx()` already spans `r ∈ [0.99, 1.57]` against a box ending at
1.4 — so the study runs on `[0.60, 1.90] × [-1.10, 1.10]`. And `Ψ_N = 0.75` is
**not measurable on a fitted rectangle at all**: enclosing it leaves under one
cell of margin at the coarsest mesh of a dyadic sweep, at which point one is
measuring the contour's distance to the mesh boundary. Its *reference* value is
asserted, which is free; the `h`-study stops at `Ψ_N = 0.50`. That is the
fixture's elongation, not a limitation of the method — and it is one more reason
the curved path is where this item actually lives.

**This is where `ROADMAP.md` item 10 becomes reachable**, and it is where cost
stops being hypothetical: a flux-surface average is the first thing here a
consumer calls in a loop. §11 is the analysis and IN-P is the harness.

### IN-3 — the representation — **DONE, 2026-09-03**

`R(ρ, θ)` and `z(ρ, θ)` in the Zernike disc basis on `ρ = √Ψ_N`.
`src/meq/SurfaceFit.{hpp,cpp}`, `tests/convergence/SurfaceFitConvergence.cpp`.
**MFEM-free**, like `Zernike`, `Profiles` and `Source`, so CI can build and test
it; `CMakeLists.txt` records that reason beside the source list.

**`ρ = √Ψ_N` is now measured rather than argued**, which was the point of keeping
`Ψ_N` as a control. Coefficient envelope on `nstx()` at `L = 20`: **4.89e-05
against 8.62e-03**, a ratio of 176; worst fit error 3.44e-05 against 7.01e-03, a
factor of 204; over `l = 10 → 20` the `ρ` column falls by 32.6 and the control by
4.5. **Conditioning is untouched by the choice** — 1.87e3 against 1.44e1 — so it
is not a conditioning artefact.

**The parity control fits the sample cloud EIGHT TIMES BETTER and is useless**:
condition number 9.15e+16, axis error 1.23e-04 against Zernike's 3.56e-07, and an
axis that moves by 4.8e-03 depending on which `θ` you approach along. **A better
residual on the data you fitted, and nothing anywhere else.**

**The axis comes out `θ`-independent for FREE and exactly** — spread `0.000e+00`
at every degree, on analytic and discrete data alike, because every mode with
`m ≠ 0` carries `ρ^|m|` and above. Asserted as an exact zero rather than a
tolerance. What is *not* free is whether that point is the magnetic axis; that is
an extrapolation into the hole, and it reads 3.56e-07 at `L = 20`.

**The derivative against independently traced surfaces** falls **276×** over
`L = 4 → 20`, and the fit's own `|dx/dθ|` against the field's `√(ρ′² + ρ²)`
reaches **3.27e-04** — two independent routes to the same metric, which is the
cheap check with real teeth.

**And the Richardson finding needed a step sweep to appear at all.** At the
natural step the extrapolated difference is worth only **1.6×**, because here the
*fit's own* derivative error is the binding constraint rather than the instrument.
Swept, it separates properly — 70.5× at a step of 0.16 — and the diagnostic is
that **the plain column falls by 50.8 across an eightfold refinement while the
extrapolated one moves 11%**: the converging column is the instrument and the flat
one is the answer. This is now the third place in the tree where a plain central
difference floors a check at the instrument rather than the quantity.

> **THE LEVER ON CONDITIONING IS THE HOLE, NOT THE LAYOUT, AND THIS FILE SAID
> OTHERWISE.** A sample set has a hole in the middle — no surface is traced at
> `Ψ_N = 0` — and the orthogonality argument needs nodes spanning the whole disc.
> Condition number against the inner limit, on the rescaled edge: 7.78 at
> `Ψ_min = 0.02`, 3.19e+02 at 0.10, **7.30e+04 at 0.25**. Four orders. The three
> sample layouts — equispaced in `Ψ_N`, equispaced in `ρ`, Gauss in `Ψ_N` — agree
> to within **25%** at every hole size, and Gauss is sometimes the worst. The test
> asserts the wrong story dead, at `layout spread < 2×`, so it cannot quietly
> return.
>
> **And rescaling the disc edge is a change of BASIS, not of model.** A Zernike
> expansion of degree `L` spans the polynomials of degree `L` in `(x, y)`, and that
> space is closed under scaling — measured, identical worst errors to better than
> 1e-6 relative at every inner limit. So the choice is *purely* conditioning and is
> worth up to **11,500×**. The default stays `1.0` only so that a change of
> coordinate is never silent.

**The number a coupling reads.** `∂(geometry)/∂Ψ_N` grows as the inner limit
falls, and the product with `1/(2√Ψ_N)` **settles at 1.148** over
`Ψ_N = 0.50 → 0.005`: an innermost node at `Ψ₁` carries a geometry derivative of
about **0.574 / √Ψ₁**. The assertion is that the product is *bounded and
settling*, which is the statement that the growth belongs to the coordinate and
not to the fit. **State an acceptance in `ρ`, or name a flux range bounded away
from the axis** — per §4.4, `∂R/∂Ψ` diverges like `1/(2√Ψ_N)`, so a criterion
written against `dGeometry/dΨ` over the whole domain would be comparing infinities
and would fail for a reason that has nothing to do with the fit.

### IN-4 — the `ψ`-element decision — **ANSWERED, 2026-09-03: it was the ANGLE**

**The question was the wrong one**, and §4.4 records why. `meq::gaugeFreeFit()` —
a geometric Gauss–Newton on `Ψ_N(x(ρ,θ)) − Ψ`, warm-started by IN-3's linear fit,
with `∇Ψ_N` from the solved flux — requires each disc node only to **land on the
right surface**. **No `ψ`-element is needed and none is implemented.**
`CLAUDE.md`'s *IN-4* has the measurements; four findings belong here.

> **THE EXPLICIT SPECTRAL-WIDTH PENALTY LOSES ON ITS OWN METRIC, AND THE REASON
> IS STRUCTURAL.** `M(p,q) = Σ m^{p+q}(R²+Z²) / Σ m^p(R²+Z²)` is a **ratio** of
> two weighted sums of the same coefficients, so a quadratic penalty is not a
> surrogate for it: twelve decades of `λ` move `M` by 1.6% **in the wrong
> direction**, and cost **44×** in surface error. Hirshman & Breslau minimise `M`
> itself, which is not a quadratic problem. Kept as the losing column.

> **AND THE GAUGE IS A SOFT TAIL WITH NO GAP, NOT A NULL SUBSPACE — THIS FILE
> SAID OTHERWISE.** Measured: the ellipse family has **exactly 3** null directions
> at every degree from 2 to 16, of 6 to 306 columns, and **`nstx` has none at
> all**; what both have is a smooth tail running down to 8e-08 of the largest
> singular value, 58 of 306 directions below `1e-4 σ_max` on the ellipses and 97
> on `nstx`. So there is nothing to project out and **the floor is a threshold
> that must be chosen rather than read off a gap**. The no-gauge control still
> fires on both fields — first step `8.6e+10 ×` the coefficient norm on the
> ellipses, Jacobian `−2.2e+13`, i.e. **folded** — but the mechanism and the
> magnitude differ by six orders between fields, so **a control measured on one
> field alone would have reported whichever it happened to meet**.

> **AND A TRUST REGION IS A THIRD TREATMENT THAT NEITHER CANDIDATE COVERED.** The
> undamped pseudo-inverse reaches round-off at `L ≤ 12` and **fails at `L = 16`
> and 20**. What makes it robust is an adaptive Levenberg–Marquardt damping —
> which is **itself a Tikhonov gauge**, so `SurfaceGauge::None` has to disable the
> damping as well as the floor. **A "no gauge" control that kept the trust region
> would pass while testing nothing.**

> **ONCE THE ANGLE IS FREE, ANY ACCEPTANCE WRITTEN AGAINST A PRESCRIBED ANGLE IS
> MEASURING THE GAUGE RATHER THAN THE FIT.** IN-3's derivative and metric checks
> compare the fit's position *at a given `θ`* against a surface traced at that `θ`
> — exactly the freedom being granted — so they are not re-asserted for the
> gauge-free fit and are replaced by **gauge-invariant** properties: distance to
> the surface, fitted perimeter against exact (3.5e-07 relative), and the sign of
> the map Jacobian. **The map must be checked for folding**: minimum Jacobian
> +5.2e-02 to +6.3e-02 with the gauge on, against negative on **every** ungauged
> run — a surface residual alone admits a beautiful number over a folded map.
>
> **This costs the consumer nothing, which is the point.** `MANTA-COUPLING.md`
> reads flux-surface *averages*, and an average over a surface does not know how
> the surface was parametrised. **The deliverable was gauge-invariant all along.**

### IN-5 — open surfaces

Chebyshev per §4.2. **Deferred with free boundary**, per §6.

### IN-6 — the output

`DRIVER-PLAN.md` §3's flux-surface `(Ψ, θ)` NetCDF grid, and whatever
`MANTA-COUPLING.md` settles on — **including the per-`ψ` cache**, which §11.1
argues is a requirement of the pointwise call pattern rather than an
optimisation, and whose invalidation contract is `MANTA-COUPLING.md` §8.

### IN-P — the performance harness

**Not a stage in the ladder and not gated on one**: it can start as soon as IN-2
exists, and it should, because every design choice after that point trades
accuracy against cost and there is currently **no cost number at all** to trade
against. §11 is the analysis; this is the measurement.

`tests/performance/`, **not a ctest**, `TraceSolverScaling` + `scan.sh` as the
pattern. Non-zero exit only for the correctness properties — chiefly that a
threaded extraction reproduces the serial one **bit for bit**, which independent
surfaces and independent rays make available and which a tolerance would
undermine.

**Do the three call-site changes of §11.3 first**, before any parallel region is
written: `Mesh::GetElementTransformation( int )` hands out shared scratch and is
used at two sites in `FluxSurfaces.cpp` and one in `CriticalPoints.cpp`. They
are correct while everything is serial and are a silent wrong answer the moment
it is not.

---

## 8. Risks, in likely-to-bite order

1. ~~**The metric trap of §3.2.**~~ **Met, measured and now a live control.** A
   spectral rule fed a differenced Jacobian is silently second order — IN-1 reads
   7.03 against 1.97 on the same trapezoid over the same points, and IN-2's
   averages carry the differenced column as a four-column control rather than as
   a remark. It was the highest-probability, lowest-visibility risk and it is the
   one the ladder was ordered around.
2. **The separatrix.** Everything degrades approaching it — `1/|∇ψ|` in the
   integrand, the corner in the surface, the logarithm no polynomial basis
   catches. Production codes simply stop: LIUQE cuts at `Ψ_N = 0.95`, FreeGS
   extrapolates outside `[0.01, 0.99]`. **Decide meq's cut deliberately and
   record it**, rather than discovering it as a convergence failure.
3. ~~**The band silently truncating an outer surface.**~~ **Closed**, and the
   prior art is what it was closed against: `v0-legacy:FluxSurfaces.cpp` returns a
   partial curve with a message on stderr, and a flux-surface average over an arc
   labelled as a closed contour is wrong by an amount nothing reports. The band is
   traced through (§4.3) and **every mask is asserted point by point against
   `mfem::Mesh::FindPoints`, in both directions** — a count is not a mask.
4. **`Δs` chosen once and left.** The whole §2 argument rests on `Δs` being
   *tuned to `h` and `k`*. A fixed `Δs` reintroduces the cap the method exists
   to avoid — and it will look like the method failing rather than the parameter
   being wrong.
5. **The axis.** `1/|∇ψ|` diverges there, surfaces shrink to a point, and the
   corrector's `∇ψ(c−ψ)/|∇ψ|²` is 0/0. IN-A first is the mitigation; the
   innermost surfaces may still want the elliptical Hessian approximation LIUQE
   uses.
6. **A rotating equilibrium's `p`.** §6's second carve-out. If anyone asks for
   pressure contours, none of the scoping in §6 applies and the topology work
   is back.
7. **Conflating IN-A's axis with `ψ_ax`.** Two different quantities, both
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

---

## 11. Performance, and where the parallelism is

**NOTHING IN THIS ITEM HAS BEEN TIMED, AND THAT IS THE FIRST THING TO SAY.**
Every number in §7 is an accuracy measurement. The only cost datum in existence
is that `FluxSurfaceConvergence` runs its whole sweep in 16.6 s, which is a
*test* and not a workload — it traces a handful of contours at several `Δs` and
throws them away. Nobody has measured the cost of the thing a consumer will
actually ask for.

### 11.1 What the cost is made of

| | |
|---|---|
| per accepted point | a few corrector iterations, each an element-local evaluation plus a walk step |
| per surface | `N` points — **and the tracer is sequential along the curve by construction**, since each predictor starts from the last accepted point |
| per family | `n_ψ` surfaces, currently traced independently and from scratch |
| per *call* | see below, and this is the one that matters |

**THE CONSUMER'S CALL PATTERN IS THE REAL PROBLEM AND `MANTA-COUPLING.md` §5
ALREADY SAYS SO.** `Geometry` is **pointwise**: called once per physics node,
per residual evaluation, handed the whole `ψ` vector each time
(`SystemSolver::evaluateGeometry` is a plain loop over points). For a
flux-surface average that means locating the surface through `x` and integrating
on it, *per node, per residual*. That document's own conclusion — "not
affordable done naively, so MEQ will need to cache per `ψ`" — is a requirement
on this item and not an aside. **The cache, and its invalidation contract, is
part of IN-6 and should be designed with the harness in front of it rather than
after.**

**AND `dGeometry_dpsi` IS WORSE BY A FACTOR OF `nFieldDOF`.** It is shaped
`(nGeometry, nFieldDOF)` — a derivative with respect to **every** `ψ` degree of
freedom. Obtained by differencing, that is `nFieldDOF` complete re-extractions
of the whole surface family. On any mesh worth solving on this dominates
everything else in the coupling by orders of magnitude, and it is the number to
design against.

### 11.2 Where the parallelism is, in decreasing order of value

1. **Surfaces are independent — embarrassingly parallel over the flux label.**
   Tracing the surface at `Ψ₁` and at `Ψ₂` share nothing but a read-only mesh
   and read-only fields. This is the obvious win and the one to take first.
   **But see §11.4: it is in direct tension with continuation.**
2. **`dGeometry_dpsi` is embarrassingly parallel over `ψ` DOFs**, and it is the
   dominant cost, so this is where a factor of the core count is actually worth
   something.
3. **The rays of the angle parametrisation are independent; the tracer's steps
   are not.** Each `θ` is a separate 1-D root solve along its own ray, sharing
   nothing. **That asymmetry is a design fact worth exploiting rather than
   regretting**: the tracer is the robust route for *establishing* a surface —
   it assumes no star-shapedness and finds the curve — and the rays are the
   parallel route for *populating* it once `AngleParametrisation` has confirmed
   `min|u × t|` is healthy. Do the sequential part once and the parallel part
   wide.
4. **Quadrature over a populated surface** is trivially parallel and is *not*
   where the time is; the cost is locating the points, not summing over them.
5. **The Zernike least squares** is one dense solve on a small system. Leave it
   to BLAS — and read §11.3 before doing even that.

### 11.3 The threading constraints are meq's own, and two of them are traps

**`MKL_NUM_THREADS=1` IS NOT NEGOTIABLE AND THIS ITEM DOES NOT GET TO RELAX
IT.** `CLAUDE.md` records `ComputeH()`'s element-local dense LU degrading by a
factor of **forty** at `k = 3` on threaded MKL, and the variable is
process-wide. So parallelism here must be **OpenMP over independent work, never
threaded BLAS**, and any harness that appears to gain from raising MKL threads
is measuring the solver getting slower somewhere else.

> **`Mesh::GetElementTransformation( int )` RETURNS A POINTER TO SHARED
> SCRATCH.** MFEM's own comment: *"The returned object is owned by the class and
> is shared, i.e., calling this function resets pointers obtained from previous
> calls."* It is one `IsoparametricTransformation` member of the `Mesh`. Two
> threads evaluating in different elements will silently corrupt each other's
> transformation — **no crash, no error, a wrong point**. The thread-safe route
> is the `( i, IsoparametricTransformation * )` overload into a thread-local.
>
> **This is not hypothetical: `src/meq/FluxSurfaces.cpp` uses the shared
> overload at two sites and `src/meq/CriticalPoints.cpp` at one.** All three are
> correct today because everything is serial. They are the concrete first item
> of any threading work, and they should be changed *before* a parallel region
> is written rather than during the debugging of one.

Three more, none of them exotic:

* **`MFEM_THREAD_SAFE` is ON in meq's build** and is what makes `FiniteElement`
  evaluation reentrant at all. Without it none of this is safe, and it is not
  meq's flag to assume elsewhere.
* **Build every cached table before the parallel region.**
  `Mesh::ElementToElementTable()` builds and caches on first call; building it
  inside a parallel region is a data race on the cache.
* **The `FindPoints` fallback must stay outside**, or be serialised. It is
  already `O(elements × points)` and the tracer reports zero fallbacks, so this
  costs nothing to honour.

**THE HARNESS MUST ASSERT BIT-EXACT REPRODUCTION OF THE SERIAL ANSWER**, at
`0.000e+00` and not at a tolerance — the precedent is
`threadedAssemblyReproducesSerialAssemblyExactly`, and the reason applies here
in a stronger form: independent surfaces and independent rays **reassociate
nothing**, so exactness is available. A tolerance would be an admission that
something is shared, which is exactly the defect being guarded against.

> **AND DO NOT MAKE IT AUTOMATIC.** `CLAUDE.md` records a gate on
> `omp_get_max_threads() > 1` being tried for assembly and **removed**: MFEM
> forks a team per call, so a caller that assembles hundreds of times inside a
> bordered Newton pays the fork every time, and `HighBetaConvergence` went 21.5 s
> → 39 s under it. The same shape is here — a parallel region per surface is
> right for a `(Ψ, θ)` grid built once and **wrong** for a single surface asked
> for repeatedly, which is precisely MaNTA's pointwise call pattern. Informed
> opt-in, like `setAssemblyMode()`.

### 11.4 Two ALGORITHMIC levers, both of which beat threading

**1. Continuation in the flux label, which falls straight out of §3.1.**
Allgower & Georg ch. 15 applied to `H(R, z, c) = ψ_h(R, z) − c` says the whole
family of surfaces is a single 2-manifold with **no critical points at all**,
since `∂H/∂c = −1`. So the surface at `Ψ_j` is an excellent predictor for the
surface at `Ψ_{j+1}`, and a corrector started there should converge in one or
two iterations where a fresh trace searches from nothing. Tracing every surface
from scratch throws that structure away.

> **AND IT IS IN DIRECT TENSION WITH §11.2's FIRST ITEM.** Continuation makes
> the surfaces **sequential**; independence makes them parallel. They cannot
> both be taken. Which wins depends on the core count, on `n_ψ`, and on how much
> the continuation actually saves — all three measurable, none of them yet
> measured. **Do not assume; this is exactly the kind of trade this project has
> repeatedly got wrong by reasoning from structure.** A hybrid is available and
> may be the answer: continue along a few chains in parallel, one chain per
> thread.

**2. The shape derivative instead of differencing, for `dGeometry_dpsi`.** A
flux-surface average is an integral over a surface whose location is itself
determined by `ψ`, so its derivative has two parts — the integrand moving and
the surface moving. Both have closed forms. Differencing costs `nFieldDOF` full
re-extractions; the closed form costs one extraction plus a quadrature per DOF.
**This is the single largest lever in the whole coupling**, larger than any
amount of threading.

> **But the differenced version is the TEST and must not be deleted.**
> `MANTA-COUPLING.md` §6: MaNTA never assembles its Jacobian, so a missing
> surface-motion term "does not produce a wrong answer, only slow Newton
> convergence", and several defects survived a passing suite for months on
> exactly that. The differenced derivative is what the closed form is checked
> against, as `dFdPsi` already is here.

### 11.5 What to measure, and where it lives

**In `tests/performance/`, and NOT as a ctest**, per the standing rule: every
number in it is a timing, and a threaded timing on this machine is a
measurement about the machine. `TraceSolverScaling` plus `scan.sh` is the
pattern to copy — one process per point, because MKL fixes its threading at
first use and an in-process sweep measures the first setting repeatedly; and a
non-zero exit **only** for the correctness properties that make the timings mean
anything.

What the harness should report:

* cost per traced point against `k`, and against the corrector tolerance
* cost per surface against `N`, separating the tracer from the ray population
* cost per family against `n_ψ`, with and without continuation
* the distribution of corrector iterations per point, and `stalledCorrections`
* `fallbackLocations`, which should stay at zero and is a cliff if it does not
* thread scaling over surfaces and over `ψ` DOFs, against the bit-exactness
  assertion
* the cache hit rate under a simulated pointwise call pattern — the MaNTA one,
  not a batch

**And no device timing without a synchronise.** If any of this ever reaches a
GPU, `CLAUDE.md` records a 148,224-unknown factorisation reading 2.0e-04 s
because the work was still queued — four orders out, and *plausible enough to
publish*.
