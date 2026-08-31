# Free boundary: HDG inside, an exact exterior operator outside

A plan, not an implementation. Nothing here has been built. Written 2026-08-29.

`CLAUDE.md` is the operational record and is authoritative on anything already
measured; `ROADMAP.md` is the order of work; `DRIVER-PLAN.md` is the stage-7
design this one follows. Free boundary is the first item after the driver that is
new physics rather than plumbing, and the ROADMAP has said for weeks that the
only thing known about it was that its prerequisite existed. This is that
sentence expanded.

## 1. What free boundary is, and what changes

Fixed boundary, which meq solves today: `Γ` is given, `ψ = 0` on it, and the
domain is what `Γ` encloses. Free boundary: **the plasma boundary is an unknown**,
the field extends to infinity through a vacuum region containing the coils, and
what is given instead is the coil currents and the machine geometry.

```
−∇̄·( (1/r) ∇̄ψ ) = F( r, z, ψ ) / r          in the plasma
−∇̄·( (1/r) ∇̄ψ ) = μ₀ j_coil                 in the coils
−∇̄·( (1/r) ∇̄ψ ) = 0                          everywhere else, out to infinity
ψ = 0                                          on the axis r = 0
ψ → 0                                          at infinity
```

with the plasma occupying `Ω_p(ψ)`, the region bounded by the last closed flux
surface, and

```
F( r, z, ψ ) = [ μ₀ r² p'(Ψ) + (g g')(Ψ) ] · χ_{Ω_p(ψ)},
Ψ = ( ψ − ψ_bnd ) / ( ψ_ax − ψ_bnd ),
ψ_ax = max ψ,   ψ_bnd = the flux at the limiter contact or the X-point.
```

**Four things are new, and they are of very different sizes.**

| | |
|---|---|
| **The unbounded domain.** | §2–§4. Structural, and the part this plan is mostly about. |
| **`ψ_bnd` as a second unknown.** | §5.2. Small — `ψ_ax` already works, and this is the same machinery once more. |
| **The plasma support `χ_{Ω_p}`.** | §5.3. The genuinely hard part, and the one CEDRES++ names as its own obstacle to higher order. |
| **The coils.** | §5.4. Ordinary. |

The nonlinearity meq already has — `ψ_ax` inside the residual, closed by a
bordered Newton — is the *pattern* for two of these four, which is the single
biggest reason this is now approachable.

## 2. The method, and why this one

**HDG on a polygonal subdomain, an exact exterior operator on a smooth
artificial boundary, coupled at a distance by the extension technique meq
already has.**

`refs/CouplingAtADistance.pdf` — Cockburn, Sayas & Solano, *Coupling at a
Distance HDG and BEM*, SIAM J. Sci. Comput. 34 (2012) A28–A47 — is the method,
and the choice of it is not arbitrary. **Its reference [5] is Cockburn & Solano,
which is meq's stage 5**: the family of paths `Σ_h`, the extension `E_h(q_h)`,
and the lifting

```
L_h(g)|_{K_ext}(x) = g( a(x) ) + ∫_σ C E_h(q_h) · m ds
```

are, term for term, `mfem::TransferPath`, `mfem::ElementExtension` and
`mfem::TransferredDatumCoefficient` — the machinery `setExtension()` drives and
`ExtensionConvergence` measures. **So the free-boundary coupling is stage 5 with
`g` unknown instead of zero.** That is the whole structural insight, and it is
why this method rather than another.

What the paper adds is: put a smooth `Γ` outside the mesh, solve the exterior
problem on it spectrally, and couple the two through the transmission condition

```
E_h(q_h)·ν + λ = 0     on Γ,        λ := the exterior normal derivative.
```

`Γ_h` need not fit `Γ` and need not even be close to it — the paper measures
optimal orders with `dist(Γ_h, Γ) = O(h)`, which is assumption P.1, which
`meq::AdaptiveDomain` already maintains through refinement.

**Why not the alternatives**, all of which meq has already looked at and
`refs/Refs.md` records at length:

* **Lackner / von Hagenow** — the classical route, and what `attic/free-boundary/`
  implemented. Rejected for one reason, in Lackner's own words: the
  precompute-and-reuse structure "will probably not be competitive with iteration
  methods if the geometry of R is changed after each calculation". meq's stage 6
  changes it every cycle.
* **A directly coupled BEM with a kernel** — CEDRES++'s `c(·,·)`, eq (3.5), a
  double surface integral over `Γ` with a complete-elliptic-integral kernel. It
  is the right operator; §3 says why meq should not assemble it that way.
* **An outer fixed point** between interior and exterior — the paper's own §4.3,
  a Richardson iteration. Rejected for the reason meq rejects every outer fixed
  point: it is outside the residual, so the Jacobian cannot see it. The
  measurement that settled that argument is in `HighBetaConvergence` and it was
  not close.

## 3. The semicircle collapses the BEM to a diagonal operator

**This is the simplification the whole plan is built on, and it was derived and
checked in the course of writing this. It is not in any of the three papers in
this form.**

Take `Γ` to be a **semicircle of radius `ρ_Γ` centred on the axis** — CEDRES++'s
choice, and Gatica & Hsiao's circular interface — so the computational domain is
a half-disc whose flat side is the axis. In three dimensions that semicircle is a
**sphere**, and the exterior operator on a sphere is diagonal.

### 3.1 The separation

In spherical coordinates `r = ρ sinθ`, `z = ρ cosθ`, `μ = cosθ`, the
Grad–Shafranov operator is

```
Δ* = ∂_ρρ + (1/ρ²)( ∂_θθ − cotθ ∂_θ )
```

— the `(2/r)∂_r` that separates `Δ*` from the axisymmetric Laplacian cancels the
`(2/ρ)∂_ρ` exactly. Separating `ψ = ρ^α f(μ)` gives

```
(1 − μ²) f'' + α(α−1) f = 0,
```

which is the **Gegenbauer equation of order −1/2**. So `α = n` or `α = 1−n`, and
the angular functions are `C_n^{−1/2}(μ) = ( P_{n−2}(μ) − P_n(μ) ) / (2n−1)`.

Three consequences, and each of them removes something:

* **`C_n^{−1/2}(±1) = 0` for every `n ≥ 2`**, so the basis satisfies `ψ = 0` on
  the axis identically. The flat side of the half-disc needs no separate
  treatment in the exterior at all.
* **The exterior modes `ρ^{1−n}` all decay**, and there is no constant mode
  compatible with `ψ(axis) = 0`. So the paper's undetermined constant `u_∞` and
  its compatibility condition `∫_Γ λ = 0` — its two fiddliest pieces, §4.3's
  whole `Ξ` apparatus — **both disappear**. The axisymmetric problem is cleaner
  than the Laplace one it is modelled on.
* **The Dirichlet-to-Neumann map is diagonal.** If `ψ|_Γ = Σ_{n≥2} a_n C_n(μ)`
  then `ψ_ext = Σ a_n (ρ/ρ_Γ)^{1−n} C_n(μ)` and

```
∂ψ/∂ρ |_Γ  =  Σ_n a_n · (1 − n)/ρ_Γ · C_n(μ).
```

### 3.2 And the weight is meq's own

The Gegenbauer functions of order `−1/2` are orthogonal in the weight
`(1 − μ²)^{−1}`. On a semicircle centred on the axis,

```
dΓ / r  =  ρ_Γ dθ / ( ρ_Γ sinθ )  =  dθ / sinθ  =  dμ / (1 − μ²).
```

**The orthogonality weight is exactly `dΓ/r`**, which is the weight the
Grad–Shafranov weak form carries anyway, and `ρ_Γ` cancels out of it. So

```
∫_Γ C_m C_n dΓ/r = δ_mn h_n,        h_n = 2 / ( n(n−1)(2n−1) ),
```

and the exterior contributes to the trace system a **diagonal** block with entry
`(n−1) h_n / ρ_Γ` — no layer potentials, no singular quadrature, no elliptic
integrals, no `O(N²)` kernel evaluations. Gatica & Hsiao's uncoupling taken all
the way.

### 3.3 Measured, because a derivation is not a result

Three checks, all run 2026-08-29 and reproducible from
`FREE-BOUNDARY-PLAN.md`'s description alone:

1. **`Δ*( ρ^α C_n^{−1/2}(cosθ) ) = 0`** for `α = n` and `α = 1−n`, by central
   differences in `(r,z)`, at `n = 2…5` and several points: residuals of
   `1e−8` to `8e−6` against `|ψ|` of `5e−3` to `0.85`, which is the difference
   floor.
2. **Orthogonality in `dΓ/r`**: off-diagonal entries at `1e−17`, and
   `h_n` agreeing with `2/(n(n−1)(2n−1))` to six digits at `n = 2…6`.
3. **The DtN against an independent exact field.** Take a circular current loop
   — `ψ` from complete elliptic integrals, which is not the formula any of the
   above came from — put a semicircle around it, expand its trace in `C_n`,
   apply the symbol `(1−n)/ρ_Γ`, and compare against its exact normal
   derivative:

   | `ρ_Γ` | modes | trace error | **DtN error** | `|a_n|` at `n = 2, 8, 16, 20` |
   |---|---|---|---|---|
   | 2.5 | 24 | 1.7e−11 | **2.3e−09** | 2.0e−1, 8.7e−4, 1.3e−6, 5.8e−8 |
   | 4.0 | 24 | 4.4e−14 | **6.8e−09** | 1.2e−1, 3.2e−5, 1.2e−9, 7.7e−12 |

   The DtN figure is limited by the finite-difference reference, not by the
   method. **This is the check that matters**, and it is deliberately the kind
   this project has learned to insist on: checking the expansion against the
   formula it was derived from would catch nothing, and the Solov'ev
   coefficients are the standing reminder.

The spectrum is the other half of that table. `|a_n|` falls geometrically, at a
rate set by `ρ_plasma / ρ_Γ`, so **`N` is small and `ρ_Γ` trades mesh against
modes**: a bigger `ρ_Γ` needs fewer modes and more elements. The paper's §5 runs
that same trade-off for its own case and it should be re-run for meq's.

### 3.4 What this claim is not

It is not a claim that CEDRES++ is doing anything wrong. Their `c(·,·)` is the
same operator in position space; they use P1 nodal elements on `Γ`, in which
basis it is dense and needs the kernel. meq is free to use a spectral trace on
`Γ` — which is what the coupling paper does, with trigonometric polynomials —
and in *that* basis the operator diagonalises. **The two must agree, and that
agreement is a test to write**: assemble CEDRES++ eq (3.5) against the
Gegenbauer basis and check it comes out diagonal with the symbol above. If it
does not, this section is wrong and the rest of the plan needs the kernel.

## 4. The coupled system, and its algebra

Unknowns, after the element-local elimination:

| | | |
|---|---|---|
| `λ` | the hybridized trace | `n_trace`, sparse |
| `a` | the Gegenbauer coefficients of `ψ` on `Γ` | `N ≈ 20–40` |
| `ψ_ax`, `ψ_bnd` | the normalisation | 2 |

Equations:

```
R( λ, a, ψ_ax, ψ_bnd ) = 0     the HDG trace residual, with the transferred
                               datum φ_h = g(a(x)) + ∫_σ r E_h(q_h)·m on Γ_h
T_m( λ, a )           = 0     ∫_Γ E_h(q_h)·ν C_m dΓ − a_m (1−m) h_m / ρ_Γ,   m = 2…N+1
G_ax                  = 0     ψ_ax − max ψ_h
G_bnd                 = 0     ψ_bnd − (limiter / X-point functional)
```

and the Jacobian is **bordered, with a sparse border**:

```
[ A    B    c_ax  c_bnd ]        A  = ∂R/∂λ, the existing hybridized Jacobian
[ Bᵀ*  D    ·     ·     ]        B  = ∂R/∂a,  supported on Γ_h faces ONLY
[ b_ax ·    d     ·     ]        D  = the DtN, DIAGONAL
[ b_bnd ·   ·     d'    ]        b_ax = the border meq already builds
```

**`B` is sparse and that is the point.** `g` enters the datum on a face of `Γ_h`
only through `g(a(x))` for `x` on that face, so `∂φ_h/∂a_n` is local to that
face, and `∂R/∂a_n` reaches only the trace dofs of the elements owning `Γ_h`
faces. Likewise `T_m` reads `E_h(q_h)` on the far end of those same paths. So the
whole augmented system is a sparse matrix of size `n_trace + N + 2` and a direct
solve handles it, provided the blocks can be assembled — which is §6.

**The fallback if they cannot** is the block elimination meq already runs for
`ψ_ax`: `N + 2` extra backsolves against one factorisation, using
`SetReuseSymbolic()`. That works, it is written, and at `N = 40` it is roughly
`40` backsolves per Newton step against one factorisation — tolerable but not
free. It is the reason §6's request is worth making rather than working around
forever.

**And this is the same shape as the bordered Newton already in the tree.** That
solver — `GradShafranovSolver::setSource( NormalisedSource &, double )` — is one
border column and one row. Free boundary is `N + 2`, and the generalisation is
mechanical.

### 4.1 Which nonlinear ordering

**`B` and `T` are assemblable, and `∂R/∂a` is not the same thing.** The raw
blocks of §6 are integrals over `Γ_h` faces and over `Γ`; they involve the
geometry and the basis and nothing else. But `R` is the residual *after* the
element-local elimination, so `∂R/∂a` is the raw block pushed through that
elimination, and that is inside MFEM either way.

**REWRITTEN 2026-08-30, because the version here was wrong.** It said that
differencing `R` in `a` works under `CondenseThenLinearise` and returns *exactly
zero* under `LineariseThenCondense`, the auxiliary unknown being invisible in a
cached `r_lin`. That came from `darcyhybridization.hpp`'s summary of the ordering
and not from the code: `Relinearise()` deliberately does **not** retain the local
residual, and `MultInvLin()`'s correction evaluates `LocalResidual()` — and so
the source — at the current fields, every time. Measured: meq's `ψ_ax` bordered
Newton, border still differenced, reaches the same `ψ_ax` under both orderings to
every digit printed, in one or two more iterations. The guard that refused the
ordering has been lifted.

So the position is simpler than it looked:

| | `CondenseThenLinearise` | `LineariseThenCondense` |
|---|---|---|
| local elimination | a nonlinear solve per element | a linear solve |
| `∂R/∂a` **differenced** | works; fragility is the local solves' convergence | works; nothing local to fail |
| `∂R/∂a` **assembled**, with §6 | the converged local Jacobian | `M`, already factored — cleaner |

**Neither ordering blocks the coupling.** What differs is cost and fragility, and
`LineariseThenCondense` is the better of the two for this: an assembled border is
a straight application of factors it already holds. `meq/NORMALISED-LINEARISE-FIRST.md`
is the design for that.

**What does still choose the ordering is the source, not the coupling.**
`LineariseThenCondense` fails at 60 on the GS-2 pedestal at every resolution
tried, including ones condense-first solves in five iterations, because its
reduced residual **is not a function of the trace**: the retained local fields are
hidden state, and at the published pedestal width `R` at one trace moves by 149%
depending on where the linearisation was last taken. Condense-first's is exactly
history-independent — and its *gradient* is the half that fails there, by 100×,
its local solves being at their cap. `CLAUDE.md` carries both tables. Neither is
a coupling problem, and a free-boundary source with a discontinuous support will
make both worse before it makes either better.

So FB-0 to FB-3 run whichever ordering the source allows, and the cut-element
source of FB-4 is where that question is settled. **The coupling is indifferent;
the plasma is not.**

## 5. The four new pieces, in meq

### 5.1 The exterior operator — small, exact, and testable alone

A class carrying the semicircle, the mode count, the symbol and the mass:

```cpp
namespace meq
{
    /// The exterior of a semicircle, as an operator on its trace.
    class ExteriorDtN
    {
        public:
            ExteriorDtN( double zCentre, double rhoGamma, int modes );

            /// C_n^{-1/2}( cos theta ) at a point of Gamma, n = 2 .. modes+1.
            double basis( int n, double r, double z ) const;
            /// ( 1 - n )/rho_Gamma, the Dirichlet-to-Neumann symbol.
            double symbol( int n ) const;
            /// h_n = 2/( n( n-1 )( 2n-1 ) ), the 1/r-weighted mass.
            double mass( int n ) const;
            /// psi at any exterior point, from the coefficients. For output,
            /// and for the field at the coils.
            double exterior( double r, double z,
                             std::vector<double> const &a ) const;
    };
}
```

Acceptance: the current-loop test of §3.3, as a unit test, plus the agreement
with CEDRES++ eq (3.5) of §3.4. **Both are exact-answer tests**, which is rare
in this subject and worth spending.

### 5.2 `ψ_bnd`, and a second border

`meq::NormalisedSource` currently takes `setNormalisation( ψ_ax )` with `ψ_bnd`
fixed at zero, and the header already says free boundary is where the second
argument goes. It becomes

```cpp
virtual void setNormalisation( double psiAxis, double psiBoundary ) = 0;
```

with `Ψ = ( ψ − ψ_bnd )/( ψ_ax − ψ_bnd )` throughout. The bordered Newton gains a
second column and row of the same shape as the first.

**`ψ_bnd` is a max of two things and that is not smooth.** It is the larger of
the flux at the limiter contact and the flux at the saddle, and the saddle's
location is itself an implicit function of `ψ`. The `ψ_ax` border is already a
semismooth Newton — the argmax dof can change between iterations — and this is
the same thing twice over. Expect it to work and expect it to need the line
search that `ψ_ax` needed; do not expect the argument to be clean.

### 5.3 The plasma support — the hard one

`χ_{Ω_p(ψ)}` makes the source discontinuous across a curve that moves with the
iterate. Three things follow.

**It stays a pointwise function of `ψ`, which is better than it looks.** With
`Ω_p` approximated as `{ Ψ > 0 }` intersected with a search region, `F` is still
`F( r, z, ψ )` given `ψ_ax` and `ψ_bnd`, so `meq::Source`'s interface survives
untouched. What is lost is that `{ Ψ > 0 }` can pick up private-flux regions and
near-coil regions that are not the plasma; CEDRES++ handles that with a
connectivity test and so must meq.

**The Jacobian acquires a surface term unless the profiles vanish at the
boundary.** `∂F/∂ψ` picks up `F·δ(Ψ)` at the plasma edge. If `p'(0) = 0` and
`(gg')(0) = 0` — which is the usual convention and which every profile in
`tests/analytic/` except `HighBetaPoloidal` satisfies — the source is continuous
and the term vanishes. **Decide this deliberately and write it down**, because a
missing surface term is exactly the kind of Jacobian error that converges to the
right answer at the wrong rate, which `CLAUDE.md` documents at length.

**Cut-element quadrature is real work and meq has none.** An element straddling
`∂Ω_p` needs a rule over `K ∩ Ω_p`, and its *derivative* with respect to the
iterate. CEDRES++'s §5 names exactly this — "quadrature over polygonal domains
with curved boundaries, plus the derivatives of that quadrature" — as what stops
them going above first order. **meq is a `k+1` code and this is where that is at
risk.** Budget accordingly, and treat "what order survives the cut" as a
measurement to make early rather than a hope.

### 5.4 The coils — ordinary, and useful early

Coil currents are data: `F_coil = μ₀ r I_k / |Ω_ck|` on each coil subdomain, or a
filament. A `[coils]` table in the TOML, a `meq::CoilSet`, and a source that adds
the coil term. Nothing structural — and it is what makes stage FB-1 below
possible, which is the acceptance test for everything in §3 and §4.

## 6. The split: what belongs in MFEM

The full request is written up as
`../mfem-hdg-dev/doc/HDG-BEM-COUPLING-FROM-MEQ.md`. In summary:

**MFEM-side, and small.** Two rectangular integrators beside the ones that exist,
both element-local in exactly the way `HDGExtensionIntegrator` is:

* `⟨ g_n ∘ a, v·n ⟩_e` on a face of `Γ_h` — the block `B`. `PathTraceCoefficient`
  already computes `g∘a` for a *given* `g`; what is missing is the rectangular
  form against a basis of them.
* `⟨ E_h(v_i)·ν, C_n ⟩_Γ` on the far end of the paths — the block `T`. The
  extension quadrature exists over `K^ext`; what is missing is its boundary
  piece, the image of a `Γ_h` face on `Γ` with the induced measure.

**MFEM-side, and structural.** `DarcyHybridization` has no notion of a global
unknown that is not a trace dof, so a rectangular block cannot survive the static
condensation. Carrying `M` auxiliary globally-coupled unknowns through the
elimination — so that `GetGradient()` returns the bordered system — is what
turns this from `M` extra backsolves into one solve. **It is also what meq's
`ψ_ax` border wants and does not have**, and it is what any global constraint on
a hybridized system would want. That makes it worth asking for on its own merits
rather than as a free-boundary special case.

**meq-side, and everything else.** The exterior operator (§5.1), the semicircular
domain and its level set, `ψ_bnd` (§5.2), the plasma support and cut quadrature
(§5.3), coils (§5.4), the augmented Newton, the configuration, and the whole
verification ladder. **The physics is meq's and the discretisation machinery is
MFEM's**, which is the same line the tree already draws.

## 7. The staged plan

Each stage ends at a measured number. That is this project's standing rule and it
is what makes free boundary tractable at all: the received wisdom, which
CEDRES++ states outright, is that no analytic free-boundary solution exists, so
validation is against a fine mesh. **That is true of the finished problem and
false of every stage below FB-4**, and the exact answers available early are what
should be spent first.

| | | acceptance |
|---|---|---|
| **FB-0** | `meq::ExteriorDtN`: the basis, the symbol, the mass. No solver. | The current-loop test of §3.3 to 1e−9, and agreement with CEDRES++ eq (3.5) |
| **FB-1** | **Vacuum only.** Coils, no plasma. The whole coupling, on a linear problem with an exact answer — the sum of the coils' loop fields. | `ψ_h` against that closed form at `k+1`, and the coupling sign pinned by it |
| **FB-2** | A **prescribed** plasma current, still linear: put a known `j_φ` inside and check the exterior. | Same rate; and `ComputeOutwardFlux` against the total current, which is the sharpest whole-assembly test available |
| **FB-3** | `ψ_bnd` as an unknown, plasma support still fixed. | Self-consistency of `ψ_bnd` to round-off, Newton order 2, exactly as `HighBetaConvergence` asserts for `ψ_ax` |
| **FB-4** | The moving plasma support and cut quadrature. | The order that survives the cut, measured, against `k+1` |
| **FB-5** | The augmented Newton as one bordered solve, and adaptivity through it. | `η` monotone through refinement with `Γ` fixed; assumption P.1 preserved |
| **FB-6** | A machine case, against a fine-mesh reference. | Convergence to the reference; agreement with CEDRES++ where a published case exists |

**FB-1 is the stage to protect.** It exercises `ExteriorDtN`, the transferred
datum with a non-zero `g`, the transmission condition, the augmented solve and
the coil sources, and it does so on a problem whose answer is known to machine
precision. If anything in §3 or §4 is wrong, FB-1 is where it shows, and it shows
as a wrong number rather than as a plausible picture.

**And the sign convention will be wrong at least once.** `DarcyForm` holds `−q`,
the papers' `τ` carries the opposite sign to what is stable, and the trace matrix
is negative definite in meq's convention. The exterior DtN enters that system
with a sign nobody will get right by reasoning, and the thing that settles it is
FB-1, not an argument. Budget a day and write down what wins.

## 8. Risks, in the order they are likely to bite

**The axis.** The half-disc includes `r = 0`, where the flux mass form `(r q, v)`
degenerates and `BoundaryShape` currently refuses to go — its constructor rejects
a surface reaching the axis, "where the operator's 1/r is not integrable". Three
things say this is survivable and none of them is a measurement: `q = (1/r)∇̄ψ` is
*bounded* at the axis because `ψ ~ r²`; the source is identically zero there in
free boundary, so `(F/r, w)` never arises; and the mass matrix is degenerate but
still positive definite on any element of positive measure. **What is unknown is
the conditioning as `h → 0`, and whether the local solves survive it.** This is
the first thing to measure and it can be measured today, with no free boundary at
all: solve a vacuum problem on a half-disc mesh touching the axis and watch the
element-local iteration counts. Do that before FB-1.

**The corner where `Γ` meets the axis.** Two right-angle junctions, and
`CLAUDE.md` records that corners are where the transfer-path analysis gives out.
Here the corner is between the arc and a fitted straight boundary rather than a
corner of `Γ` itself, and the lifting's weight `C = r` vanishes there, so the
transferred datum degenerates to `g(a(x)) → 0` — probably benign, definitely not
established.

**Cut quadrature and the order.** §5.3. The one place where a published code says
it hit a wall.

**Vertical instability.** CEDRES++ names vertically unstable plasmas as the case
where fixed-point iteration fails outright, and `refs/LacknerFreeBoundary.pdf`
describes the axis-pinning feedback that a fixed-point scheme needs to survive
it. meq's answer is that it is not a fixed-point scheme — but a Newton on an
indefinite problem is not automatically safe either, and the line search that the
bordered Newton needed is the shape of the answer.

**`N`, `ρ_Γ` and the mesh.** §3.3 shows the trade-off exists and is geometric.
It has not been run for a tokamak geometry, where the plasma is elongated and the
paper's own caveat applies: "the introduction of a circular interface may require
a large computational domain in situations where the support of source terms is
very elongated". Measure the spectrum on a real coil set before choosing `ρ_Γ`.

## 9. Deliberately out of scope

* **Iron.** Ferromagnetic structures make `K − I` nonlinear inside `Γ`. The
  method accommodates it — that is exactly the paper's hypothesis — but it is a
  second nonlinearity and belongs after FB-6.
* **The inverse problem.** Finding coil currents to achieve a shape. CEDRES++
  does it as SQP with Tikhonov regularisation reusing the same derivatives, which
  is the right design and is only available once the derivatives exist.
* **The evolution problem.** Quasi-static resistive evolution, CEDRES++'s
  Problems 5 and 6.
* **Anything from `attic/free-boundary/`.** The von Hagenow implementation there
  is the good algorithm badly amortised, and §2 says why the amortisation is what
  adaptivity destroys. It is worth reading and not worth restoring.
