# Free boundary: HDG inside, an exact exterior operator outside

A plan, not an implementation. Nothing here has been built. Written 2026-08-29,
**rewritten 2026-09-01 against the NPC API**, which changed the answer to the
two questions this plan is mostly about: what the coupled Jacobian looks like,
and what has to come from MFEM before any of it can be tried.

`CLAUDE.md` is the operational record and is authoritative on anything already
measured; `ROADMAP.md` is the order of work. `DRIVER-PLAN.md` was the stage-7
design this one follows and is now that stage's record.

**Two things have arrived since this was written and both bear on §5.3.**
`src/meq/CriticalPoints.{hpp,cpp}` locates the axis and any X-point as roots of
`q = 0`, sub-element, with a Poincaré–Hopf audit — which is machinery the
connectivity test below wants, and which did not exist when §5.3 was written.
And `INVERSION-PLAN.md`'s **IN-5, open surfaces, is deferred to arrive with this
item**: a disc chart has no meaning through a separatrix and an angle about the
axis has none on an open field line, so the in-surface coordinate there has to be
poloidal arc length normalised to `2π` from a fixed-`z` reference ray. Whoever
starts FB-4 should read `INVERSION-PLAN.md` §4.2 and §6 first.

**What the rewrite changed, in one paragraph.** The August version was written
against `NLOrdering::LineariseThenCondense`, chose it over the condensation for
the coupling, and concluded that the coupled system could not reach the reduced
solve without a new capability in `DarcyHybridization` — §6's "structural ask".
That ordering **has since been deleted from MFEM** as a condensation in
disguise, and MEQ's default is now `NonlinearOrdering::NPC`: Newton on the full
`(q, ψ, ψ̂)` system, with the Jacobian solved by hybridized elimination. Under
NPC **the residual is unreduced**, so the derivative of the residual with
respect to an auxiliary unknown is a *raw* block and there is no elimination for
it to survive. The structural ask is therefore an optimisation and not a
prerequisite, and **FB-0 to FB-3 need nothing from MFEM that is not already
there**. §6 is the section that changed most, and §6.4 is the per-stage note.

## 1. What free boundary is, and what changes

Fixed boundary, which MEQ solves today: `Γ` is given, `ψ = 0` on it, and the
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
| **`ψ_bnd` as a second unknown.** | §5.2. Small — `ψ_ax` already works, and under NPC it works more cheaply than it used to. |
| **The plasma support `χ_{Ω_p}`.** | §5.3. The genuinely hard part, and the one CEDRES++ names as its own obstacle to higher order. |
| **The coils.** | §5.4. Ordinary. |

The nonlinearity MEQ already has — `ψ_ax` inside the residual, closed by a
bordered Newton — is the *pattern* for two of these four, which is the single
biggest reason this is now approachable.

## 2. The method, and why this one

**HDG on a polygonal subdomain, an exact exterior operator on a smooth
artificial boundary, coupled at a distance by the extension technique MEQ
already has.**

`refs/CouplingAtADistance.pdf` — Cockburn, Sayas & Solano, *Coupling at a
Distance HDG and BEM*, SIAM J. Sci. Comput. 34 (2012) A28–A47 — is the method,
and the choice of it is not arbitrary. **Its reference [5] is Cockburn & Solano,
which is MEQ's stage 5**: the family of paths `Σ_h`, the extension `E_h(q_h)`,
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

**Why not the alternatives**, all of which MEQ has already looked at and
`refs/Refs.md` records at length:

* **Lackner / von Hagenow** — the classical route, and what `attic/free-boundary/`
  implemented. Rejected for one reason, in Lackner's own words: the
  precompute-and-reuse structure "will probably not be competitive with iteration
  methods if the geometry of R is changed after each calculation". MEQ's stage 6
  changes it every cycle.
* **A directly coupled BEM with a kernel** — CEDRES++'s `c(·,·)`, eq (3.5), a
  double surface integral over `Γ` with a complete-elliptic-integral kernel. It
  is the right operator; §3 says why MEQ should not assemble it that way.
* **An outer fixed point** between interior and exterior — the paper's own §4.3,
  a Richardson iteration. Rejected for the reason MEQ rejects every outer fixed
  point: it is outside the residual, so the Jacobian cannot see it. The
  measurement that settled that argument is in `HighBetaConvergence` and it was
  not close — coupled 8.3e-02 → 4.4e-15 in four steps against decoupled
  8.3e-02 → 8.2e-02 in fifteen. **It does not converge slowly; it does not move.**

## 3. The semicircle collapses the BEM to a diagonal operator

**This is the simplification the whole plan is built on, and it was derived and
checked in the course of writing this. It is not in any of the three papers in
this form.** Nothing in this section is affected by the NPC rewrite: it is
geometry and special functions, and it was measured.

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

### 3.2 And the weight is MEQ's own

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

and the exterior contributes to the coupled system a **diagonal** block with
entry `(n−1) h_n / ρ_Γ` — no layer potentials, no singular quadrature, no
elliptic integrals, no `O(N²)` kernel evaluations. Gatica & Hsiao's uncoupling
taken all the way.

### 3.3 Measured, because a derivation is not a result

Three checks, all run 2026-08-29 and reproducible from this section's
description alone:

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
that same trade-off for its own case and it should be re-run for MEQ's.

### 3.4 What this claim is not

It is not a claim that CEDRES++ is doing anything wrong. Their `c(·,·)` is the
same operator in position space; they use P1 nodal elements on `Γ`, in which
basis it is dense and needs the kernel. MEQ is free to use a spectral trace on
`Γ` — which is what the coupling paper does, with trigonometric polynomials —
and in *that* basis the operator diagonalises. **The two must agree, and that
agreement is a test to write**: assemble CEDRES++ eq (3.5) against the
Gegenbauer basis and check it comes out diagonal with the symbol above. If it
does not, this section is wrong and the rest of the plan needs the kernel.

## 4. The coupled system under NPC

**REWRITTEN 2026-09-01. The August version of this section described a border on
the condensed trace residual, and that is no longer the system MEQ solves.**

### 4.1 The ordering question is settled, and it settles itself

The August version weighed `CondenseThenLinearise` against
`LineariseThenCondense` and preferred the latter for the coupling. **MFEM
deleted `LineariseThenCondense` on 2026-08-31** as *"a condensation in
disguise"*, and `SetNonlinearOrdering()` went with it. What remains is
`CondenseThenLinearise` and **`NPC`** — `mfem::DarcyNPCOperator` with
`mfem::DarcyNPCSolver`, Nguyen, Peraire & Cockburn eqs (14)–(18) — which is
MEQ's default.

The design that section pointed at — `NORMALISED-LINEARISE-FIRST.md`, a
mechanism for carrying `ψ_ax` under the deleted ordering — **is deleted with
it**. It was a design for a mode that no longer exists, and its one measurement
was that differencing the reduced residual in an auxiliary unknown returns a
derivative rather than zero under *that* ordering: a fact about a code path that
is gone. Under NPC the question does not arise, because the residual is never
reduced before it is differenced. Recoverable from git if the argument is ever
wanted.

**Under NPC the coupling question is not a question.** The reason is one line of
the method: the unknown is the whole `(q, ψ, ψ̂)` vector and the residual is the
**unreduced** `F(q, ψ, ψ̂)`. An auxiliary unknown's border column is therefore a
derivative of an *evaluated* residual, not of one reconstructed from a
linearisation, and there is no static condensation for a rectangular block to
fail to survive. **The whole of §6's "structural ask" was a consequence of
condensing before differentiating.**

### 4.2 The unknowns and the equations

| | | |
|---|---|---|
| `x = (q, ψ, ψ̂)` | the full NPC state | `n_flux + n_pot + n_trace`, and `ψ` is an *independent* unknown |
| `a` | the Gegenbauer coefficients of `ψ` on `Γ` | `N ≈ 20–40` |
| `ψ_ax`, `ψ_bnd` | the normalisation | 2 |

```
F( x, a, ψ_ax, ψ_bnd ) = 0     the full HDG residual, with the transferred datum
                               φ_h = g(a(x)) + ∫_σ r E_h(q_h)·m imposed on Γ_h
T_m( x, a )           = 0      ∫_Γ E_h(q_h)·ν C_m dΓ − a_m (1−m) h_m / ρ_Γ
G_ax                  = 0      ψ_ax − max ψ_h
G_bnd                 = 0      ψ_bnd − (limiter / X-point functional)
```

### 4.3 Where `a` enters, and why its column is constant

**`g` enters as an essential trace value on `Γ_h`, and nowhere else.** That is
how MEQ already imposes the fixed-boundary datum: `SetEssentialBC` covers every
boundary attribute, and `prepare()` projects the datum onto the trace. Today
`Γ_h`'s trace dofs are pinned to **zero** — `ProjectBdrCoefficient` is called
with `fittedMarker` rather than `dirichletMarker`, because with `g ≡ 0` the
whole datum is the lift term that `HDGExtensionIntegrator` deposits in the flux
mass block. Free boundary un-pins them:

```
ψ̂|_{Γ_h}  =  P a,        P_{in} = the trace projection of ( C_n ∘ a )|_{Γ_h}
```

`P` is `n_{Γ_h trace} × N`, and each of its columns is one call to
`ProjectBdrCoefficient` against `mfem::PathTraceCoefficient( path, C_n )` —
which is the class that already exists, taking an arbitrary `PositionFunction`.
**The paths do not move during a solve, so `P` is assembled once.**

**And the border column is constant in the iterate, which is worth more than it
looks.** `ψ̂` enters the flux row as `⟨ψ̂, v·n⟩`, the potential row as `⟨τψ̂, w⟩`
and the trace row as `⟨τψ̂, μ⟩` — **linearly in all three**, every one of them.
All the nonlinearity is `F(r, z, ψ)` in the potential row, which depends on `ψ`
and not on `ψ̂`. So `∂F/∂a = (∂F/∂ψ̂)·P` does not depend on where the iterate is,
and it can be built once per mesh rather than once per Newton step.

**That is an argument and not a measurement, and it is a cheap one to make**:
build the column at two well-separated iterates and difference them. Do it in
FB-1, where the whole problem is linear anyway and any discrepancy is a defect
rather than a nonlinearity. It is exactly the shape of claim this file's own
history says to check — see the block-structure claim that upstream corrected in
`CLAUDE.md`'s *Why it fails*.

### 4.4 The bordered solve, which MEQ already runs at `N = 1`

```
[ J    R  ] [ dx ]     J = the NPC Jacobian, factored by NPCGradient()
[ T    D  ] [ da ]     R = −P in the Γ_h essential trace rows, zero elsewhere
                       T = ∂(transmission)/∂x, on the flux dofs of Γ_h elements
                       D = the DtN, DIAGONAL
```

with the `ψ_ax` and `ψ_bnd` rows and columns beside them. Block elimination
against one factorisation:

```
DarcyNPCSolver lin( traceSolver );
Operator &S = npc.GetGradient( x );   // ONE local factorisation, ONE Schur
lin.SetOperator( S );                 // ONE symbolic + numeric trace factor
lin.Mult( F, y );                     // backsolve 1
for ( n : columns ) lin.Mult( R_n, Z_n );   // backsolves 2 .. N+3
// then a dense (N+2) solve for da, dψ_ax, dψ_bnd, and dx = y − Z·(those)
```

**`DarcyNPCSolver::SetOperator` re-points the trace solver at `S` and
`Mult()` is reduce → backsolve → recover**, so the cost of the border is
`N + 2` extra *backsolves* and no extra factorisation. With
`SetReuseSymbolic()` the symbolic analysis is done once for the whole Newton
run — which is `A quarter of every Newton step was thrown away` in `CLAUDE.md`,
paying off in a place it was not written for.

**MEQ already does exactly this at `N = 1`.** `solveWithNormalisation()` builds
a `DarcyNPCOperator`, solves `J y = R` and `J z = c` through a
`DarcyNPCSolver`, and closes the border with `δs = (b·y − G)/(d − b·z)`. Free
boundary is the same function with the scalar corner replaced by a dense
`(N+2)×(N+2)` one. **The generalisation is mechanical and it is MEQ's own
code.**

### 4.5 Two of the three borders are exact under NPC, and that is new

The August plan treated every border entry as a finite difference. Under NPC
they are not:

| | condensation | **NPC** |
|---|---|---|
| `b_ax = −∂(max ψ_h)/∂x` | `3(k+1)` central differences over one element's trace dofs | **exactly `−e_j`** — `max ψ_h` is one entry of the unknown |
| `d = ∂G_ax/∂ψ_ax` | `1 −` a central difference | **exactly `1`** |
| `∂F/∂a` | the raw block pushed through the local elimination | **`(∂F/∂ψ̂)·P`, raw, constant** — §4.3 |
| `∂F/∂ψ_ax` | one central difference in a scalar | the same, and it stays differenced |

**And the frozen-seed hazard goes away entirely.** Under the condensation, every
differenced border is a difference of a residual whose element-local Newtons are
seeded from a vector captured at `FormLinearSystem()` time and never refreshed;
on a hard problem those hit their cap and return something that is not a
function of anything, and MEQ measured a `9e−6` perturbation moving the
recovered peak from `0.896` to `3.84`. `formSystem()` exists to work round it.
**NPC has no element-local nonlinear solve, so there is no seed, nothing to go
stale and nothing to re-form** — which removes the single sharpest edge the
August plan had to plan around.

## 5. The four new pieces, in MEQ

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
in this subject and worth spending. The loop field needs complete elliptic
integrals; C++17's `std::comp_ellint_1` / `_2` supply them and no dependency is
needed.

**It belongs in `meq_core` and must stay MFEM-free**, like `Profiles` and
`Source`: it is special functions and geometry, it is under the `naming` check,
and it is unit-testable without the library.

### 5.2 `ψ_bnd`, and a second border

`meq::NormalisedSource` currently takes `setNormalisation( ψ_ax )` with `ψ_bnd`
fixed at zero, and the header already says free boundary is where the second
argument goes. It becomes

```cpp
virtual void setNormalisation( double psiAxis, double psiBoundary ) = 0;
```

with `Ψ = ( ψ − ψ_bnd )/( ψ_ax − ψ_bnd )` throughout. The bordered Newton gains a
second column and row of the same shape as the first.

**Under NPC the limiter case is exact and the X-point case is not.** If `ψ_bnd`
is the flux at a limiter contact point it is a nodal value of `ψ`, which is an
unknown, so its border row is `−e_j` exactly, as `ψ_ax`'s is. If it is the flux
at a saddle it is `ψ` evaluated where `∇̄ψ = 0`, and the saddle's location is an
implicit function of the state — a differentiable functional in principle,
through the implicit function theorem, and one whose derivative MEQ has no
assembled route to. Expect that one to stay differenced.

**`ψ_bnd` is a max of two things and that is not smooth.** It is the larger of
the two, and the argmax can change between iterations. The `ψ_ax` border is
already a semismooth Newton for the same reason — its argmax *dof* can change —
and this is the same thing twice over. Expect it to work and expect it to need
the line search that `ψ_ax` needed; do not expect the argument to be clean.

### 5.3 The plasma support — the hard one, and MFEM has more of it than we thought

`χ_{Ω_p(ψ)}` makes the source discontinuous across a curve that moves with the
iterate. Four things follow.

**It stays a pointwise function of `ψ`, which is better than it looks.** With
`Ω_p` approximated as `{ Ψ > 0 }` intersected with a search region, `F` is still
`F( r, z, ψ )` given `ψ_ax` and `ψ_bnd`, so `meq::Source`'s interface survives
untouched. What is lost is that `{ Ψ > 0 }` can pick up private-flux regions and
near-coil regions that are not the plasma; CEDRES++ handles that with a
connectivity test and so must MEQ. **`meq::CriticalPointFinder` is the piece that
makes that test cheap** — the X-point is what separates the private flux from the
plasma, and it is now locatable sub-element as a root of `q` rather than confined
to a mesh vertex, which is what CEDRES++ records as an open problem on its own P1
discretisation.

**The Jacobian acquires a surface term unless the profiles vanish at the
boundary.** `∂F/∂ψ` picks up `F·δ(Ψ)` at the plasma edge. If `p'(0) = 0` and
`(gg')(0) = 0` — which is the usual convention and which every profile in
`tests/analytic/` except `HighBetaPoloidal` satisfies — the source is continuous
and the term vanishes. **Decide this deliberately and write it down**, because a
missing surface term is exactly the kind of Jacobian error that converges to the
right answer at the wrong rate. `CLAUDE.md`'s *A wrong Jacobian is invisible to a
convergence table* is the measurement: perturbing `∂F/∂ψ` by 5% leaves every
error and every rate unchanged to six figures and only drops Newton's observed
order to 1.000.

**MFEM HAS CUT-ELEMENT QUADRATURE AND THIS PLAN SAID IT DID NOT.** The August
version's "MEQ has none" was a statement about MEQ and was allowed to stand as a
statement about the stack. `fem/intrules_cut.hpp` carries
`mfem::CutIntegrationRules` and `mfem::MomentFittingIntRules`, giving both the
cut-volume and the cut-surface rule for the zero level set of a `Coefficient`, in
2D and 3D. It is gated on `MFEM_USE_LAPACK`, which MEQ's build has, and the
header **is installed** in `../mfem/install`. Algoim is a second backend and is
not needed. So the rule itself is available today and should be tried before
anything is written.

**What is still missing is the DERIVATIVE of the rule**, which is the half
CEDRES++ actually names — "quadrature over polygonal domains with curved
boundaries, **plus the derivatives of that quadrature**" — as what stops them
going above first order. Moving `ψ` moves the cut, which moves the points and the
weights, and `MomentFittingIntRules` has no interface for that sensitivity. Two
honest options, and the choice is a measurement rather than an argument:

* **Ignore it**, and accept an inconsistent Jacobian on cut elements only. That
  is a Jacobian error of the kind above: it costs Newton's order and not the
  answer, and the answer is what the rate table measures. Cheap, and it may be
  enough.
* **Difference it**, per cut element, which is `O(cut elements)` extra rule
  constructions per Jacobian and is affordable because the cut set is `O(h⁻¹)`.

**MEQ is a `k+1` code and this is where that is at risk.** Treat "what order
survives the cut" as a measurement to make early rather than a hope — it is
FB-4's acceptance criterion for that reason.

### 5.4 The coils — ordinary, and useful early

Coil currents are data: `F_coil = μ₀ r I_k / |Ω_ck|` on each coil subdomain, or a
filament. A `[coils]` table in the TOML, a `meq::CoilSet`, and a source that adds
the coil term. Nothing structural — and it is what makes FB-1 possible, which is
the acceptance test for everything in §3 and §4.

## 6. The split: what belongs in MFEM

**REWRITTEN 2026-09-01, AND THE ANSWER IS MUCH SMALLER THAN IT WAS.** The full
request is `../mfem-hdg-dev/doc/HDG-BEM-COUPLING-FROM-MEQ.md`, filed 2026-08-29.
That document asks for two things: §2, two rectangular integrators, and §3, a
structural capability for auxiliary globally-coupled unknowns. **Under NPC MEQ
needs neither of them to start**, and the reasons are §4.3 and §4.4. It is worth
being precise about why, because the request as filed overstates what is
blocking.

### 6.1 Why §2.1 is not needed — the block `B`

The request asks for `⟨φ_n ∘ a, v·n⟩_e`, "the datum's data half as a rectangular
form against a basis". **MEQ does not need that block, because the datum's data
half is not a form — it is an essential trace value.** §4.3: `ψ̂|_{Γ_h} = P a`,
and `P`'s columns are `ProjectBdrCoefficient` against `PathTraceCoefficient`,
which exists and takes an arbitrary `PositionFunction`. The rectangular
integrator would be the right thing if the datum entered weakly; it does not.

### 6.2 Why §2.2 is reachable, though it is the one worth asking for — the block `T`

The transmission rows need `E_h(q_h)·ν` on `Γ` — the far end of the paths —
against the basis and the induced measure there. `ExtensionRegionQuadrature`
sweeps the *region* `K^ext_e` and there is no boundary variant, so this is the
one piece with nothing directly behind it. But every primitive is public:

* `TransferPath::Endpoint( FTr, ip, xbar )` gives `a(x)` at a face quadrature
  point, virtual and public;
* the surface Jacobian of `ξ ↦ a(x(ξ))` is a central difference along the face,
  which is what `ExtensionRegionQuadrature` does for its own `t`-face and what
  its `fd_step` parameter is for;
* `ElementExtension::SetElement` + `TransformBack` evaluate the owning element's
  polynomial at a point outside it — the pattern `meq::Sampler::extendOutward`
  already uses for the `.nc` band.

So `T` is about forty lines of MEQ code. **Write it in MEQ first and ask for it
upstream afterwards, with the tiling check attached** — summing the boundary
weights over the faces must give `|Γ|`, exactly as summing the volume weights
must give `|Ω| − |D_h|`, and that check is what says the path family covers `Γ`
once. A version that has been used is a better request than a version that has
not.

### 6.3 Why §3 is an optimisation and not a prerequisite

The request's §3 asks `DarcyHybridization` to carry `M` auxiliary
globally-coupled unknowns through the static condensation, so that
`GetGradient()` returns the bordered matrix. **That is a requirement of a
condensation and not of NPC.** Under NPC the residual is unreduced, the border
column is raw, and the elimination happens once inside a block solve the caller
drives — `N + 2` backsolves against one factorisation, which is what MEQ's
`ψ_ax` border already costs and pays.

**It is still worth having**, and the request's own §3.1 makes the general case
better than free boundary does: an auxiliary unknown coupled to a hybridized
system is the shape of every global constraint, every mean-value condition and
every rank-`M` nonlocal term. What it buys here is one solve instead of `N + 2`
backsolves. At `N = 40` and a trace factorisation that dominates, that is a
modest fraction of a Newton step, and it is **not** worth blocking on.

**The request should be revised to say so.** As filed it calls §3 "structural",
which was true of the ordering it was written against and is not true of the one
MEQ runs. Leaving it standing would have another team build the larger thing
first.

### 6.4 So: what is needed from MFEM, per stage

**Nothing blocks FB-0 through FB-3. FB-4 has one real gap and it is not the one
the plan predicted.**

| stage | needs from MFEM | status |
|---|---|---|
| **FB-0** `ExteriorDtN` | **nothing** — MEQ-side, MFEM-free by design | clear |
| **FB-1** vacuum + coils + the whole coupling | **nothing.** `P` from `PathTraceCoefficient` (§6.1); `T` from `Endpoint` + `TransformBack` (§6.2); the bordered solve from `DarcyNPCOperator` / `DarcyNPCSolver`, which MEQ already drives at `N = 1` | clear |
| **FB-2** prescribed plasma current | **nothing** new beyond FB-1 | clear |
| **FB-3** `ψ_bnd` unknown | **nothing** — the `ψ_ax` border is the pattern and it is MEQ's own code | clear |
| **FB-4** moving support + cut quadrature | the **sensitivity of a cut rule** to the level set. `MomentFittingIntRules` gives the rule (§5.3) and no derivative. Not blocking — difference it per cut element, or accept an inconsistent Jacobian and measure the cost in Newton's order | **the one real gap** |
| **FB-5** one bordered solve | §3 of the request: auxiliary unknowns carried through the elimination. An optimisation over `N + 2` backsolves | wanted, not blocking |

**Two things to ask for anyway, on their own merits and not as blockers**:
§2.2's boundary quadrature on `Γ` with its tiling check, once MEQ has written
one and used it; and §3's auxiliary unknowns, which MEQ's `ψ_ax` border wants
today and which every global constraint on a hybridized system wants.

**And one thing to keep watching rather than ask for.** `DarcyNPCOperator`'s
`Jacobian` handle is **solve-only** — its `Mult()` aborts, because after
`ComputeH()` the local arrays hold factored blocks — so the coupled Jacobian
cannot be *applied*, only inverted. That is fine for the block elimination above
and it is fatal to any scheme that wants a matrix-vector product with the full
coupled operator, a Jacobian-free Krylov method over the border included. If
free boundary ever wants one, that is the constraint to design around.

**MEQ-side, and everything else.** The exterior operator (§5.1), the semicircular
domain and its level set, `ψ_bnd` (§5.2), the plasma support and cut quadrature
(§5.3), coils (§5.4), the augmented Newton, the configuration, and the whole
verification ladder. **The physics is MEQ's and the discretisation machinery is
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
| **FB-A** | **The axis.** A vacuum solve on a half-disc mesh touching `r = 0`. No free boundary, no coupling. | Element-local iteration counts and the trace condition number bounded under refinement. See §8 — this is the first thing to measure and it can be measured today |
| **FB-0** | `meq::ExteriorDtN`: the basis, the symbol, the mass. No solver. | The current-loop test of §3.3 to 1e−9, and agreement with CEDRES++ eq (3.5) |
| **FB-1** | **Vacuum only.** Coils, no plasma. The whole coupling, on a linear problem with an exact answer — the sum of the coils' loop fields. | `ψ_h` against that closed form at `k+1`; the coupling sign pinned by it; **and `∂F/∂a` constant between two well-separated iterates**, which is §4.3's unmeasured claim |
| **FB-2** | A **prescribed** plasma current, still linear: put a known `j_φ` inside and check the exterior. | Same rate; and `ComputeOutwardFlux` against the total current, which is the sharpest whole-assembly test available |
| **FB-3** | `ψ_bnd` as an unknown, plasma support still fixed. | Self-consistency of `ψ_bnd` to round-off, Newton order 2, exactly as `HighBetaConvergence` asserts for `ψ_ax` |
| **FB-4** | The moving plasma support and cut quadrature. | The order that survives the cut, measured, against `k+1` — **and the cost of an inconsistent cut Jacobian measured in Newton's observed order**, which decides §5.3's two options |
| **FB-5** | The augmented Newton as one bordered solve, and adaptivity through it. | `η` monotone through refinement with `Γ` fixed; assumption P.1 preserved |
| **FB-6** | A machine case, against a fine-mesh reference. | Convergence to the reference; agreement with CEDRES++ where a published case exists |

**FB-1 is the stage to protect.** It exercises `ExteriorDtN`, the transferred
datum with a non-zero `g`, the transmission condition, the augmented solve and
the coil sources, and it does so on a problem whose answer is known to machine
precision. If anything in §3 or §4 is wrong, FB-1 is where it shows, and it shows
as a wrong number rather than as a plausible picture.

**And the sign convention will be wrong at least once.** `DarcyForm` holds `−q`,
the papers' `τ` carries the opposite sign to what is stable, the trace matrix is
negative definite in MEQ's convention, and `HDGExtensionIntegrator`'s `+1` was
itself settled by measurement rather than by argument. The exterior DtN enters
that system with a sign nobody will get right by reasoning, and the thing that
settles it is FB-1. Budget a day and write down what wins.

**One caution on the globalisation, which is new since August.** The reactive
ladder — Newton, and on observed failure `PicardThenNewton` — is what the driver
runs, and `CLAUDE.md`'s *Should `PicardThenNewton` simply be the default?*
records that it must not be made predictive, because on an under-resolved mesh
the three routes reach discrete solutions differing by up to **9.4%**. Free
boundary starts every adaptive run on exactly such a mesh. **Do not let a
free-boundary failure be answered by quietly changing the globalisation**; it
changes which equilibrium is reported.

## 8. Risks, in the order they are likely to bite

**The axis, and it is FB-A because it can be measured now.** The half-disc
includes `r = 0`, where the flux mass form `(r q, v)` degenerates and
`BoundaryShape` currently refuses to go — its constructor rejects a surface
reaching the axis, "where the operator's 1/r is not integrable". Three things say
this is survivable and none of them is a measurement: `q = (1/r)∇̄ψ` is *bounded*
at the axis because `ψ ~ r²`; the source is identically zero there in free
boundary, so `(F/r, w)` never arises; and the mass matrix is degenerate but still
positive definite on any element of positive measure. **What is unknown is the
conditioning as `h → 0`.** Under NPC there is no element-local nonlinear solve to
watch, so the diagnostic changes: watch the trace solve and the local
factorisation rather than a local iteration count. Do it before FB-1.

**The corner where `Γ` meets the axis.** Two right-angle junctions, and
`CLAUDE.md` records that corners are where the transfer-path analysis gives out
— which is why `ExtensionConvergence` takes `Γ` to be `ψ = −0.03` rather than the
separatrix through the X-point. Here the corner is between the arc and a fitted
straight boundary rather than a corner of `Γ` itself, and the lifting's weight
`C = r` vanishes there, so the transferred datum degenerates to `g(a(x)) → 0` —
probably benign, definitely not established.

**Cut quadrature and the order.** §5.3. The one place where a published code says
it hit a wall, and the one row of §6.4 with a real gap in it.

**Vertical instability.** CEDRES++ names vertically unstable plasmas as the case
where fixed-point iteration fails outright, and `refs/LacknerFreeBoundary.pdf`
describes the axis-pinning feedback that a fixed-point scheme needs to survive
it. MEQ's answer is that it is not a fixed-point scheme — but a Newton on an
indefinite problem is not automatically safe either, and the line search that the
bordered Newton needed is the shape of the answer. Note also that a line search
on the full NPC residual has been measured making **every** MEQ case worse; see
`CLAUDE.md`'s *Why it fails*. Whatever globalisation this needs, it is not that
one.

**`N`, `ρ_Γ` and the mesh.** §3.3 shows the trade-off exists and is geometric.
It has not been run for a tokamak geometry, where the plasma is elongated and the
paper's own caveat applies: "the introduction of a circular interface may require
a large computational domain in situations where the support of source terms is
very elongated". Measure the spectrum on a real coil set before choosing `ρ_Γ`.

**The border cost, if `∂F/∂a` turns out not to be constant.** §4.3 argues it is,
from the weak form. If FB-1 says otherwise, the column is rebuilt every Newton
step at `N` residual evaluations — cheap under NPC, where a residual evaluation
carries no local nonlinear solve, and ruinous under the condensation. That
asymmetry is one more reason the coupling belongs on NPC.

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
