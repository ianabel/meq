# Free boundary: HDG inside, an exact exterior operator outside

A plan, not an implementation. Nothing here has been built. Written 2026-08-29,
**rewritten 2026-09-01 against the NPC API**, which changed the answer to the
two questions this plan is mostly about: what the coupled Jacobian looks like,
and what has to come from MFEM before any of it can be tried.

`CLAUDE.md` is the operational record and is authoritative on anything already
measured; `ROADMAP.md` is the order of work. `docs/` is the manual and carries the stage-7
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

### 3.5 CEDRES++'s form read off the rendered page, 2026-09-04

The falsifying test of §3.4 needs their operator stated exactly. It has now been
transcribed from `refs/CEDRES.pdf` at 900 dpi — **not from `pdftotext`, which
mangles this particular page in two silently fatal ways; see `CLAUDE.md`'s
tooling warning, which this extends.** What follows is what the page says.

**The equation number and the notation in §3.4 above are RIGHT.** It is (3.5) on
printed page 13, it is called `c(·,·)`, and it is introduced as *"a bilinear
form `c : V × V → ℝ`, accounting for the boundary conditions at infinity"*.

```
c( psi, xi ) := (1/mu0) ∫_Γ psi(P1) N(P1) xi(P1) dS1
              + (1/(2 mu0)) ∫_Γ ∫_Γ ( psi(P1) − psi(P2) ) M(P1,P2) ( xi(P1) − xi(P2) ) dS1 dS2

M(P1,P2) = k / ( 2π (r1 r2)^{3/2} ) · ( (2 − k²)/(2 − 2k²) E(k) − K(k) )
N(P1)    = (1/r1) ( 1/δ₊ + 1/δ₋ − 1/ρ_Γ )
δ_±      = √( r1² + ( ρ_Γ ± z1 )² )
k        = √( 4 r_j r_k / ( (r_j + r_k)² + (z_j − z_k)² ) )
```

`k` is the **modulus**, not the parameter — the paper writes the radical
explicitly, which is exactly the character `pdftotext` deletes.

**THE `dΓ/r` WEIGHT IS THERE, AND IT IS WHY THE EXPONENT IS `3/2`.** They never
write `c` in weighted form, but splitting `(r₁r₂)^{3/2} = r₁ r₂ · (r₁r₂)^{1/2}`
gives

```
M dS1 dS2 = [ k / ( 2π √(r1 r2) ) ( (2−k²)/(2−2k²) E − K ) ] · ( dS1/r1 )( dS2/r2 )
```

and `N` carries a leading `1/r₁` of its own. **So `c` is naturally an integral
against `dΓ/r` in each slot** — precisely the measure §3.2 needs, with the `3/2`
being one power per slot for the measure plus a half shared by the kernel. That
is real support for §3.2 from an independent source, and it is also the reason
the exponent matters: read as `2` or as `1/2` the weight becomes `dΓ/r²` or
`dΓ`, and the §3.4 test would fail for a transcription reason.

**THREE THINGS TO CARRY INTO THE §3.4 TEST.**

* **`M` is HYPERSINGULAR, `M ~ 1/(π r d²)`** as the two points approach —
  measured over `d = 1e−1 … 1e−5`, agreeing to seven digits. The double-difference
  form *is* the regularisation: each factor is `O(d)` and the product cancels the
  `1/d²` exactly. So the test **cannot** assemble `∫∫ C_m M C_n` directly; it must
  keep the difference structure or take a Hadamard finite part. If diagonality
  comes out wrong, suspect this before the algebra.
* **`N` carries a `−1/ρ_Γ`** with no counterpart in a pure DtN mode sum. It is
  the only place `ρ_Γ` appears besides `δ_±`, so it is presumably what makes the
  truncated form exact at finite radius. Do not drop it when comparing against
  `(n−1) h_n / ρ_Γ`.
* **They do not print Lackner's Green's function**, only the resulting `M` and
  `N`, and delegate the derivation to Grandgirard 1999 Ch. 2.4, which is not in
  `refs/`.

**TWO NEGATIVE FINDINGS, both clean rather than a failure to look.**

* **On choosing `ρ_Γ`: nothing.** The symbol occurs three times in the whole
  paper — once defining the semicircle, twice inside `N`. No guidance, no
  sensitivity study, no numerical value for any of their ITER or WEST cases, and
  no statement that the answer is independent of it. §8's *"measure the spectrum
  on a real coil set before choosing `ρ_Γ`"* has no prior art to lean on.
* **On the axis: nothing either, and they go there.** Their domain includes
  `r = 0` — `∂Ω = Γ ∪ Γ_{r=0}` — with `ψ = 0` imposed both in (2.5) and in the
  space (3.1), whose norms are weighted `r` and `r^{−1}`. Their triangulation
  reaches the axis. **They report no conditioning difficulty, no loss of order
  and no special treatment**, and the axis is absent from §5's own list of known
  accuracy limitations. So CEDRES++ neither corroborates nor contradicts §7.2's
  measured `O(1/h)` — they are P1 throughout and never took the measurement that
  would show it.

**And one correction to §3.4's wording.** It says their operator is dense
*"because they use P1 nodal elements on `Γ`"*. The paper never discusses density
at all, and the causation is slightly off: density follows from a **non-local
kernel against a basis with local support**, and P1 is merely one such basis —
P4 would be dense too. What a spectral trace changes is not the density but
whether the basis **diagonalises** the operator. The conclusion §3.4 draws is
unaffected; the reason should be stated the other way.

### 3.6 §3 verified independently, 2026-09-04. All six claims stand.

Re-derived from scratch in sympy and mpmath, deliberately without reading
`src/meq/ExteriorDtN.cpp`, so the reference values are independent of the
implementation they check. **Five of the six are exact** — symbolic or exact
rational, not to a tolerance.

| | claim | |
|---|---|---|
| 3.1 | `Δ*` separates with no `∂_ρ` term | **exact**, residual identically 0 |
| 3.1 | the Gegenbauer equation and `α ∈ {n, 1−n}` | **exact** |
| — | `C_n^{−1/2}` from Legendre = library `gegenbauer(n,−1/2,μ)` | **exact, normalisation factor 1** |
| 3.1 | `Δ*(ρ^α C_n) = 0`, both branches, **in (r,z)** | **exact**, n = 2…8, 16 cases |
| 3.1 | `C_n(±1) = 0` | **exact**, n = 2…12 |
| 3.2 | orthogonality and `h_n` | **exact rational**; 72 off-diagonals identically 0 |
| 3.2 | `dΓ/r = dμ/(1−μ²)`, `ρ_Γ` cancels | confirmed |
| 3.3 | the DtN symbol, including its sign | confirmed |

Stronger than §3.3's own run in two places: the off-diagonals are **identically
zero** rather than 1e−17, and `h_n` agrees to **all 60 digits** rather than six.

**FOUR ADDITIONS WORTH ACTING ON, NONE OF THEM A CORRECTION TO THE MATHEMATICS.**

**1. Two exact identities §3 does not record, and the implementation wants both.**

```
C_n( μ )   = ( 1 − μ² ) P'_{n−1}( μ ) / ( n( n − 1 ) )     exact, n = 2..12
dC_n/dμ    = − P_{n−1}( μ )                                exact, n = 2..12
```

The first is the accuracy fix recorded in §7.3. The second means the derivative
needs one Legendre evaluation and no difference at all — which is what
`ExteriorDtN.cpp` now does.

**2. The zero at the axis is SIMPLE, and the weight is over-cancelled.**
`C_n/(1−μ²)` is exactly `1/2` at `μ = +1` and exactly `(−1)^n/2` at `μ = −1`, for
every `n`. So `C_n ~ (1−μ²)/2` there, and the integrand `C_m C_n/(1−μ²)`
*vanishes* linearly at the endpoints rather than merely staying finite. The
singular weight is not just cancelled, it is beaten — which is why a plain
Gauss-Legendre rule is exact rather than adequate.

**3. THE `n = 0, 1` ARGUMENT IS STRONGER THAN §3 CLAIMS, AND THE COUNT IS WRONG.**
The degenerate sector is **four** functions, not two — `1`, `z/ρ`, `ρ` and `z`,
all `Δ*`-harmonic — and **not one of them vanishes on the axis**. Two of them
(`1` and `z/ρ`) also decay, so both would otherwise be admissible exterior modes
and both are killed by the axis condition. Better still: they have **infinite
norm** in `L²(dΓ/r)`, so there is nothing for a compatibility condition to be
imposed *on* — a cleaner statement than "no constant mode is compatible". And
nothing admissible is lost by starting at 2: `{C_n}_{n≥2}` is **complete** in
`L²(dμ/(1−μ²))`.

Also: **the printed formula is valid for `n ≥ 2` only, and fails silently below
it** — `(P_{n−2} − P_n)/(2n−1)` returns `1 − μ` for both `n = 0` and `n = 1`,
where the true functions are `1` and `−μ`. It does not merely divide by zero; it
returns the wrong function. Recorded in the header so nobody extends the loop
downward to see.

**4. `(n−1) h_n / ρ_Γ = 2/( n( 2n − 1 ) ρ_Γ )`**, exactly, and one factor
shorter.

**AND A METHOD WARNING FOR ANYONE RE-RUNNING CLAIM 1.** Substituting
`θ = atan2(r,z)` leaves sympy unable to reduce `atan(tan θ)` and it reports a
**nonzero** residual — an artefact, not a refutation. The control settles it: the
axisymmetric Laplacian, whose spherical form is known, shows the identical
artefact. The clean route is the forward chain rule, which never forms an
`atan`.

**ONE DOCUMENTATION GAP, AND IT IS THIS PLAN'S.** §3.3's `|a_n|` spectrum table
cannot be reproduced. The decay *rate* is confirmed and is internally consistent
— both rows imply a source radius of 1.149 and 1.132, i.e. geometric decay at
exactly `ρ_source/ρ_Γ` — but the magnitudes need the loop radius and current
normalisation, which the table does not record. Against a unit loop they differ
by an `n`-dependent factor rising to about 9, identically in both rows. **Record
the source parameters beside that table** or it cannot be re-derived; the rate,
which is the part the cost argument uses, is sound.

### 3.4 What this claim is not

It is not a claim that CEDRES++ is doing anything wrong. Their `c(·,·)` is the
same operator in position space; they use a basis with LOCAL SUPPORT on `Γ` — P1
nodal elements — and against a non-local kernel any such basis gives a dense
block and needs the kernel. (**Corrected 2026-09-04**: this used to say "dense
because P1", which puts the causation on the degree rather than on the support.
P4 would be dense too. See §3.5.) MEQ is free to use a spectral trace on
`Γ` — which is what the coupling paper does, with trigonometric polynomials —
and in *that* basis the operator diagonalises.

**THE TEST IS WRITTEN AND SECTION 3 SURVIVES IT, 2026-09-06.**
`cedres_boundary_form_is_diagonal_in_this_basis` in
`tests/unit/ExteriorDtNTests.cpp` assembles eq (3.5) — the single-layer `N`
term and the hypersingular double-layer `M` term, kept in its double-difference
form — against the Gegenbauer basis at `ρ_Γ = 1.3`:

| | |
|---|---|
| worst off-diagonal | **1.15e-10** against a diagonal scale of 2.56e-01 |
| worst diagonal, relative to `blockEntry( n )` | **3.20e-09** |

so `μ₀ c( C_m, C_n ) = δ_mn (n−1) h_n / ρ_Γ` to about ten significant figures.
**This is the most independent check in the tree**: a boundary integral with
elliptic-integral kernels against a separation of variables, sharing the
equation and nothing else — no basis, no measure written the same way, no
quadrature. It is MFEM-free and runs in 30 ms, so CI gates it.

**Two things made it computable and both are recorded in the test**: the two
distances close on the semicircle, `δ_± = 2ρ_Γ cos(t/2)` and `2ρ_Γ sin(t/2)`,
so `N` is elementary and its `1/t²` at the poles is cancelled by every `C_n`
vanishing like `t²/2`; and the inner integral must be **split at the diagonal**,
where what survives the double-difference regularisation is `Δ² log Δ` — a
tensor rule straddling it does not converge.

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

`P` is `n_{Γ_h trace} × N`, and each of its columns is one projection of
`mfem::PathTraceCoefficient( path, C_n )` — which is the class that already
exists, taking an arbitrary `PositionFunction`.
**The paths do not move during a solve, so `P` is assembled once.**

**CORRECTION, 2026-09-04: it is NOT `ProjectBdrCoefficient`, and that sentence
used to say it was.** `GridFunction::ProjectBdrCoefficient` evaluates through
the **element** transformation, and a path coefficient needs the **face** one —
the path family may want the outward normal of `Γ_h`. MFEM aborts rather than
coping, which is the good case:

```
PathTraceCoefficient must be evaluated on a face: the path family may need
the outward normal of Gamma_h
```

`mfem::TransferredDatumCoefficient`'s own header says the same about itself, so
the requirement was documented one class over. The projection is written out in
`GradShafranovSolver::projectPathTraceOntoGammaH()` — which is
`projectOntoTrace()` with the coefficient evaluated on `*ftr` rather than on
`*ftr->Elem1`, and the loop restricted to boundary elements carrying a `Γ_h`
attribute.

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

> **MEASURED 2026-09-05, AND MOST OF THIS SECTION IS SUPERSEDED BY IT.**
> `tests/convergence/PlasmaEdgeConvergence.cpp`. The two options above were
> "ignore the cut-rule derivative" and "difference it", and **neither arises**,
> because no cut rule is worth adopting.
>
> * **The cap is the profile's, not the quadrature's.** With `p' ~ Ψ^j` at the
>   edge the exact `ψ` carries `|d|^{j+2}`, so the L2 BEST APPROXIMATION — no
>   solver, no quadrature question — caps `ψ*` at `min(k+2, j+2.5)`. `k+2`
>   therefore needs `k ≤ j`, for an exact cut rule and a blind Gauss rule alike.
>   Measured across `j = 0…3` and `k = 1…4`, MEQ's plain rule crosses at exactly
>   that threshold, with `j = 3, k = 3` reading **4.989** against 5.
> * **For the ORDER, a cut rule would buy `ψ_h` alone**, `j+1.5 → j+2.5`. It
>   would buy `q_h` nothing — `q` is at its own regularity bound already — and
>   `ψ*` nothing that raising the plain Gauss order does not already buy,
>   measured 2.87 → 3.57 against a bound of 3.5.
> * **BUT IT WOULD MAKE `j = 0` SOLVABLE, WHICH IS THE REAL CASE FOR ONE.**
>   With a fixed rule the assembled residual is *discontinuous* in the unknowns
>   at `j = 0` — a quadrature point crossing the edge makes `F` jump there, and
>   the measured step does not shrink as the sampling interval is quartered
>   (2.871e-04, 2.909e-04, 2.928e-04, against an integral of 8.6e-03). A rule
>   that follows the level set makes it continuous and differentiable, with the
>   derivative carrying §5.3's own `∮ F φ/|∇ψ|`. The prize is a second-order
>   solve where there is currently none, not `k+2`.
> * **There is no missing derivative *while there is no cut rule*.** With a
>   fixed rule the quadrature points do not move, so the assembled Jacobian is
>   the exact derivative of the assembled residual — of a discontinuous function
>   at `j = 0`, which is the catch. §6.4's "one real gap" is created by adopting
>   a cut rule rather than closed by it, and the same is true of every route to
>   high order at the edge: fitting the mesh, enriching the space and
>   transferring across an interface all put the geometry into the
>   discretisation, and all three then owe the Jacobian its derivative.
> * **The surface term is real at `j = 0` and it is fatal rather than
>   expensive.** With the support read off `ψ_h`, a source with `p'(0) ≠ 0` does
>   not converge at any degree or mesh, **not from the exact solution**, and not
>   under `PicardThenNewton`. `j ≥ 1` is a precondition of MEQ's free-boundary
>   path, not a convention.
> * **MFEM's cut backends are quadrilateral-only** and MEQ's meshes are
>   triangles; the Algoim path aborts with *"supports only quads and hexes"*.
>   See `refs/CutElementQuadratureSurvey.pdf`, whose own conclusion is that
>   boundary-tessellating methods are capped at second order and only
>   higher-order boundary representations escape it.
>
> **AND THE ROUTE OUT OF THE CAP, IF ONE IS EVER WANTED, IS
> `PLASMA-EDGE-PLAN.md`** — the plasma edge as an interior interface coupled at
> a distance, which is `refs/CouplingAtADistance.pdf` applied to `Γ_p` instead
> of to the exterior operator, and which reuses stage 5 and FB-1 almost
> entirely. It is a design and explicitly not to be built until `j ≥ 1` is
> finished and measured.

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

**§6.1 IS NOW CONFIRMED BY BUILDING IT RATHER THAN ARGUED**, 2026-09-04: `P`
exists, its columns are essential trace values, and §7.4 measures them. The
rectangular integrator the request asks for would be right if the datum entered
weakly, and it does not.

### 6.1 Why §2.1 is not needed — the block `B`

The request asks for `⟨φ_n ∘ a, v·n⟩_e`, "the datum's data half as a rectangular
form against a basis". **MEQ does not need that block, because the datum's data
half is not a form — it is an essential trace value.** §4.3: `ψ̂|_{Γ_h} = P a`,
and `P`'s columns are `ProjectBdrCoefficient` against `PathTraceCoefficient`,
which exists and takes an arbitrary `PositionFunction`. The rectangular
integrator would be the right thing if the datum entered weakly; it does not.

### 6.2 Why §2.2 is reachable, though it is the one worth asking for — the block `T`

**AND IT WAS ASKED FOR, AND GRANTED — 2026-09-05. `mfem::ExtensionBoundaryQuadrature`
is in the library**, merged into `gf-hdg-subdomains-dev` from MEQ. What follows
is the argument as it stood before that, kept because it is why the request was
small enough to be worth making and because the primitives it names are still
what the row is built on. §7.5 records what the filing bought.

The transmission rows need `E_h(q_h)·ν` on `Γ` — the far end of the paths —
against the basis and the induced measure there. `ExtensionRegionQuadrature`
sweeps the *region* `K^ext_e` and there was no boundary variant, so this was the
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

**That is exactly what happened, and the tiling check is what earned the
merge.** It found that an unsigned weight counts a folded sweep twice, which
upstream then fixed in the *region* sweep as well — a defect on their aerofoil,
not on MEQ's disc, that a request without a used implementation behind it would
never have surfaced.

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
| **FB-1** vacuum + coils + the whole coupling | **nothing.** `P` from `PathTraceCoefficient` (§6.1); `T` from `ExtensionBoundaryQuadrature`, which MEQ wrote and upstream merged 2026-09-05 (§6.2, §7.5); the bordered solve from `DarcyNPCOperator` / `DarcyNPCSolver`, which MEQ already drives at `N = 1` | clear |
| **FB-2** prescribed plasma current | **nothing** new beyond FB-1 | clear |
| **FB-3** `ψ_bnd` unknown | **nothing** — the `ψ_ax` border is the pattern and it is MEQ's own code | clear |
| **FB-4** moving support + cut quadrature | **NOTHING, and the gap this row named is closed by not opening it.** MEQ adopts no cut rule, so there is no rule-sensitivity to supply: with fixed quadrature points the assembled Jacobian is already the exact derivative of the assembled residual. See §5.3's measurement | **clear** |
| **FB-5** one bordered solve | §3 of the request: auxiliary unknowns carried through the elimination. **Still an optimisation and still not taken** — MEQ's border costs `N + 2` backsolves against one factorisation, which is affordable, plus ONE RE-ASSEMBLY per accepted step because the datum is a load term and `prepare()` is where a load is built | wanted, not blocking |

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

**`../freegs4e` IS THE BENCHMARK FOR FB-6, AND IT RUNS ON THIS MACHINE.**
Recorded 2026-09-05. It self-describes as a *"free boundary tokamak plasma
equilibrium Grad–Shafranov solver"*, and it is FreeGS-derived: `TestTokamak`,
`DIIID`, `MAST`, `MASTU`, `TCV`, coils, Green's-function free boundary by von
Hagenow, 2nd/4th-order finite differences on a uniform `(R,Z)` grid, Picard with
adaptive blending and an optional Newton–Krylov polish. **That is a genuinely
different algorithm from HDG-plus-Newton**, which is the whole point of it.

**So FB-6's acceptance should be against `freegs4e`, not against a fine mesh.** A
fine-mesh self-comparison shares every convention, every sign and every misread
paper with the code being checked, and this file records what that costs three
times over. `freegs4e` shares none of them, and it solves **the same equation**:
`Δ*ψ = −μ₀RJ_φ` with `μ₀RJ_φ = μ₀R²p′ + gg′`, `ψ` in Wb/rad — MEQ's `F`, MEQ's
sign, MEQ's units.

**`GeneralPprimeFFprime` is the class to drive**, because it takes tabulated
`p′` and `ff′` directly rather than solving for a normalisation to hit a target
`p_axis` and `I_p`. That makes `F` a known function of the flux on both sides,
which is what turns a comparison into a measurement.

**BEWARE `../geq`, WHICH IS A DIFFERENT THING.** It is a thin wrapper on
`freegs4e` for rotating **magnetic mirrors**: it sets `ffprime = fpol = fvac = 0`,
so `g ≡ 0`, there is no toroidal field, and its safety factor is identically
zero. A mirror cannot stand in for FB-6's machine case. Its value is to
the rotating source instead — see *Toroidal flow* in `CLAUDE.md`, where its
independent implementation of Abel (136) is the first outside check of the
`C′(ψ)` term.

**THE INTERIM COMPARISON, AVAILABLE BEFORE FB-1 LANDS.** MEQ cannot solve a free
boundary yet, but the benchmark does not have to wait for it: run `freegs4e`
free-boundary, take its converged LCFS, `ψ_ax`, `ψ_bnd` and profiles, fit the
LCFS to MXH — MEQ has the shape and no fitter, and the fit is a dozen lines of
numpy — and solve the **fixed**-boundary problem on it with the same `F`. That
checks MEQ's elliptic solve against a free-boundary answer without needing MEQ's
free boundary to work, and it is the natural rehearsal for FB-6.

**Three things to reconcile, all measured rather than guessed.**

* **`freegs4e`'s profiles are `dp/dψ_n` and `fdf/dψ_n`** — derivatives with
  respect to **normalised** flux. MEQ's `MHDSource` takes `dp/dψ`. The
  conversion divides the derivative column by `(ψ_bnd − ψ_ax)` while the
  abscissa maps `ψ_n → ψ`. That is exactly the trap `examples/rotating-density.dat`
  and its normalised twin exist to document, met from the other side: hand over
  the wrong one and it parses, solves and converges to a plasma whose gradient
  is out by a constant factor.
* **MEQ's ψ is zero on `Γ`**, so `ψ_MEQ = ψ_fgs − ψ_bnd`, and `ψ_n` is affine in
  it. Fixing `ψ_ax` to the converged value makes MEQ's tabulated `p′(ψ)` a fixed
  function of `ψ`, which is a well-posed semi-linear problem and the right thing
  to compare — the normalisation closure is a separate question.
* **The grids are transposed**: MEQ writes `psi(Z,R)` at 129², `freegs4e` works
  in `(R,Z)`. One transpose and one interpolation, and drop MEQ's
  `extrapolated` nodes first.

**THE ENVIRONMENT IS THE HARD-WON PART, SO IT IS WRITTEN DOWN.** `freegs4e`
does not run out of the box here and `geq` does not run at all. What works:

```sh
python3 -m venv venv && venv/bin/pip install numpy scipy matplotlib h5py shapely numba
PYTHONPATH=/home/ian/projects/freegs4e venv/bin/python your_driver.py
```

Deliberately **not** `pip install -e /home/ian/projects/freegs4e`, which enforces
`requirements.txt` — `numpy<2.0`, `numba~=0.60`, `Shapely~=2.0.6` — and **none of
those has a wheel for this machine's Python 3.14**. Shapely then tries to build
from source and fails on a missing `geos_c.h`, which is not installed. Relaxed,
`numpy 2.5.2`, `shapely 2.1.2` and **`numba 0.67`** all resolve and `freegs4e`
imports and builds a `TestTokamak`.

**`numba` is not optional in practice, though it looks it.** `critical.py` wraps
its import in `try/except ImportError` and the fallback calls `warnings.warn`
**without importing `warnings`**, so the no-numba path raises `NameError` and has
evidently never been exercised. Install `numba` or patch a copy; do not expect
the fallback to work.

| | | acceptance |
|---|---|---|
| **FB-A** | **The axis.** A vacuum solve on a mesh touching `r = 0`. No free boundary, no coupling. | **DONE, 2026-09-04 — see §7.2.** `ψ` at `k+1` on a mesh reaching the axis; `q` short by half an order; the conditioning penalty `O(1/h)` and not `O(1/h²)`. `tests/convergence/AxisConvergence.cpp` |
| **FB-0** | `meq::ExteriorDtN`: the basis, the symbol, the mass. No solver. | **DONE, 2026-09-04, AND COMPLETED 2026-09-06 — see §7.3 and §3.4.** The current-loop test reads **1.4e−15** in the trace and **6.9e−14** in the DtN. §3.4's CEDRES++ agreement — the falsifying test of the whole of §3 — is now written and green: their boundary form comes out diagonal to **1.15e-10** against a scale of 2.56e-01, with the diagonal matching `blockEntry( n )` to **3.20e-09** relative |
| **FB-1** | **Vacuum only.** The whole coupling, on a linear problem with an exact answer. | **DONE 2026-09-05 — see §7.8.** `ψ` at **1.99 / 2.99 / 3.99** on the half-disc with the datum given (FB-1a), and the transmission condition recovers the exterior coefficients to **1.9e-04, converging at 3.30** (FB-1b). `∂F/∂a` needs no measurement under NPC — §7.4 |
| **FB-2** | A **prescribed** plasma current, still linear. | **DONE 2026-09-05 — see §7.9.** `ψ` at 1.99 / 2.88 / 3.01, and Ampère's law through the solve: **round-off over `Γ_h`** and `k+1`-convergent on the half-disc. Two meshing findings came out of it, both about aligning the mesh to geometry that is known in advance |
| **FB-3** | `ψ_bnd` as an unknown, plasma support still fixed. | **DONE 2026-09-05.** `setNormalisation( ψ_ax, ψ_bnd )` exists, the profiles take `(ψ − ψ_bnd)/(ψ_ax − ψ_bnd)`, the validation is on the SPAN, and the second BORDER is closed: `solveWithNormalisation()` does 2×2, with `ψ_ax`'s residual at **1.88e-16** and **0.00e+00** on two meshes. It is the cheaper of the two borders — `ψ_bnd`'s dof is fixed at setup where `ψ_ax`'s needs an argmax — and `HighBetaConvergence` is bit-identical, which is what says the generalisation reduces |
| **FB-4** | The moving plasma support and cut quadrature. | **ANSWERED 2026-09-05, AND THE ANSWER MOVED THE WORK RATHER THAN DOING IT — see §7.10.** The order is capped by the PROFILE and not by the quadrature: with `p' ~ Ψ^j` at the edge, `ψ*` keeps `k+2` exactly when **`k ≤ j`**, and that threshold is the same for an exact cut rule as for MEQ's plain one. `j = 3, k = 3` reads **4.989** against a target of 5. **No cut quadrature was built**, and there is no inconsistent-cut-Jacobian cost to measure because there is no cut rule — the quadrature points do not move, so the assembled Jacobian is exact. What IS measured is that `j = 0` does not converge at all |
| **FB-5** | The augmented Newton as one bordered solve, and adaptivity through it. | **THE BORDERED SOLVE IS BUILT AND MEASURED, 2026-09-06; adaptivity through it is not.** `setExteriorCoupling()` carries the N Gegenbauer coefficients as unknowns of the same Newton as `psi_ax` and `psi_bnd`, and `solveWithNormalisation()` now does a general `( N + 2 )` elimination against ONE factorisation. On FB-1b's half-disc it agrees with the superposition route to **2.3e-15**, takes **one** Newton step because the residual is affine in `( x, a )`, and converges to the exact coefficients at **3.32** against FB-1b's 3.30. `HighBetaConvergence` is bit-identical. What remains is the adaptive loop through the coupling: eta monotone with `Gamma` fixed, and P.1 preserved |
| **FB-6** | A machine case, against **`../freegs4e`**. **The driver pieces landed 2026-09-06**: `[[coils]]` reaches the solve through `meq::makeCoilSet` and the `CoilAugmentedSource` adapters, `[source] ConfineToPlasma` switches the source off outside `{ Ψ > 0 }`, and `tools/mesh/halfdisc.py` generates the half-disc-with-conductors mesh MFEM reads natively. What is still missing for a machine case is `ψ_bnd` as an unknown (FB-3) and the exterior coupling from a config file. | Agreement with an independent free-boundary tokamak code by a **different algorithm** — von Hagenow Green's functions, finite differences, Picard — on the same coils and the same tabulated `p′`, `ff′` through its `GeneralPprimeFFprime`. A fine-mesh self-comparison is the fallback, not the target: it shares every convention with the code it checks. See §7 |

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

### 7.1 FB-A can be measured against a closed form, which this plan did not know

**Two corrections to FB-A's acceptance as written above, made 2026-09-04.**

**The stale half.** It asked for "element-local iteration counts ... bounded
under refinement". Under `NonlinearOrdering::NPC`, which is MEQ's default, there
is no element-local non-linear solve at all — `GetNumLocalNLIterations()` is
identically zero and `SolverContract` asserts it. §8 already notes this and says
to watch the trace solve and the local factorisation instead; the table did not.
It does now. (A vacuum solve is linear anyway, so the count would have been zero
for a second, unrelated reason.)

**The weak half, and this is the useful part.** "Bounded under refinement" is a
much poorer statement than this tree accepts anywhere else, and it is not
necessary. A vacuum solve needs `Δ*ψ = 0`, and there are polynomial
`Δ*`-harmonic functions that vanish **identically on the axis**:

```
ψ = r²            ψ = r² z            ψ = r⁴ − 4 r² z²
```

Verified symbolically rather than asserted, with three controls — one of them a
candidate that **failed**, which is why the check was run:

| candidate | `Δ*` | `ψ(r=0)` | |
|---|---|---|---|
| `r²` | 0 | 0 | |
| `r² z` | 0 | 0 | |
| `r⁴ − 4 r² z²` | 0 | 0 | |
| `r² (r² − 4z²) z` | **−16 r² z** | 0 | a guess, and it was wrong |
| `r² ln r − z²` | 0 | **−z²** | harmonic, and NOT zero on the axis |
| `r⁴/8` | **r²** | 0 | `CLAUDE.md`'s own check, so the operator is MEQ's |

The last row is the control that matters: `Δ*(r⁴/8) = r²` is the value
`CLAUDE.md` prints under *The two papers disagree about the sign of the Solov'ev
source*, so the operator differentiated here is the one MEQ solves.

**So FB-A gets a rate against an exact answer, like every other stage in this
tree, instead of a boundedness claim.** `k+1` in `ψ` and `q` on a mesh reaching
`r = 0` is a far sharper statement than "the condition number did not blow up",
and it costs nothing to make.

**`r² ln r − z²` is the interesting control.** It is `Δ*`-harmonic, bounded at
the axis, and **not zero there** — one of Cerfon & Freidberg's twelve terms. A
solve that handles `r²` and fails on this one has found something about the axis
rather than about the mesh, which is exactly the discrimination §8's three
non-measurements are missing.

**And FB-A needs no new code to pose.** `MeshConfig` already permits `RMin = 0`
and says so explicitly (`src/meq/Config.hpp`: *"a box reaching r = 0 contains the
coordinate singularity. **That is allowed**, but is rarely what is wanted"*). It
is `meq::BoundaryShape` that refuses the axis, in its constructor, and FB-A does
not use one: it is the fitted path with a box whose inner edge is the axis. So
"it can be measured today" is literally true — the configuration is expressible
now.

**The one thing FB-A must also run is the control at `RMin > 0`**, on the same
sequence and degrees. Without it a condition-number column says nothing: some
growth with refinement is what any `h`-refinement does, and the question is
what the *axis* adds. This is the same requirement `theTransferredDatumRestoresEtaFive`
imposes on itself by keeping the pinned-zero column.

### 7.2 FB-A IS DONE, AND IT ANSWERED THE UNKNOWN. Measured 2026-09-04.

`tests/convergence/AxisConvergence.cpp` and `tests/analytic/VacuumHarmonic.hpp`,
registered as a ctest. Four dyadic meshes, `k = 1, 2, 3`, on a box whose inner
edge is exactly `r = 0`, against the identical study standing 0.25 clear of it.

**THE POTENTIAL IS UNHARMED AND THE FLUX IS NOT QUITE.**

| `k` | `ψ`, axis / control | `q`, axis / control |
|---|---|---|
| 1 | **2.000 / 2.000** | 1.79 / 2.00 |
| 2 | **3.000 / 3.000** | 2.51 / 3.00 |
| 3 | **4.000 / 4.000** | 3.999 / 3.00 |

`ψ` converges at `k+1` on a mesh reaching the axis, at every degree, to three
decimal places and identically with the control. `q` is short by about half an
order at `k = 2` and by a fifth and worsening at `k = 1`. `k = 3` is clean, and
that is a property of this fixture rather than of the method: its `q` is a
quadratic, so `P_3` has room to spare. **A fixture with a higher-degree flux
should be expected to show the deficit at `k = 3` too, and measuring that is the
obvious next step.**

**AND THE CONDITIONING QUESTION — §8's stated unknown — HAS AN ANSWER: `O(1/h)`.**

| `k` | on axis, over the sweep | the control | the ratio grows by |
|---|---|---|---|
| 1 | 5.48e1 → 4.61e2 | 7.28 → 8.54 | **7.169** |
| 2 | 1.24e2 → 1.05e3 | 1.00e1 → 1.13e1 | **7.536** |
| 3 | 1.83e2 → 1.53e3 | 1.07e1 → 1.17e1 | **7.631** |

The on-axis column **doubles with every halving of `h`**; the control **settles**.
So the axis costs exactly one power of `h`, at every degree, and the far field
costs nothing. Over an eightfold refinement that is 7.2 to 7.6 against the 8 a
clean `1/h` would give.

**Both findings have one mechanism, and it is the weight rather than the
singularity.** The flux mass form is `(r q, v)`, so the element touching the
axis carries a weight of order `h`: its diagonal is the smallest in the system,
which is the `1/h`; and it is the element the method controls least while an
unweighted `L2` error norm counts it in full, which is the half order. §8's
three non-measurements were about `q` being *bounded* and the mass matrix being
*positive definite* — both true, and both measured here — but neither was the
thing that gives way.

**WHAT THIS MEANS FOR THE REST OF THE PLAN, and the answer is that it does not
block FB-1.** `1/h` is survivable: MEQ solves the trace system with a **direct**
solver, whose cost and accuracy are nearly insensitive to conditioning at 2D
serial sizes, and `ψ` — which is what the coupling in §4 transmits, through
`ψ̂|_{Γ_h} = P a` — keeps full order. What would have stopped FB-1 is `1/h²`, and
it is not that. **It would matter to an iterative trace solve**, which MEQ does
not use and which `CLAUDE.md`'s *Do not reach for AMG* argues against at these
sizes anyway.

**Three smaller things FB-A settled on the way.**

* **`q` is bounded at the axis**, §8's first non-measurement, now a number:
  `q_r → 2.120000` as `r → 10⁻¹, 10⁻³, 10⁻⁶, 10⁻⁹, 0`, exactly, and finite *at*
  `r = 0` rather than a NaN. `VacuumHarmonic` carries the `1/r` cancellation
  already done for that reason — writing `gradPsi()/r` as every other fixture
  does would be `0/0` there.
* **A vacuum solve is affine and step one is exact**, `‖r₁‖/‖r₀‖` between 1e-13
  and 5e-12. Asserted as the **drop** and not as an iteration count, because
  written the obvious way it failed at the finest mesh with two iterations —
  `CLAUDE.md`'s *One more test moved from the stopping rule to the property*,
  met again from scratch.
* **`MeshConfig` already allows `RMin = 0`** and `meq::BoundaryShape` refuses
  the axis, so FB-A runs on the fitted path with no shape and needed no new
  library code at all. Only a fixture and a test.

### 7.3 FB-0 IS DONE. Measured 2026-09-04.

`src/meq/ExteriorDtN.{hpp,cpp}`, in `meq_core`, **MFEM-free** as §5.1 requires
and under the `naming` check; `tests/unit/ExteriorDtNTests.cpp` is the
acceptance and is a **unit** test, so CI — which cannot obtain the MFEM branch —
can gate this stage. `tests/analytic/CurrentLoop.hpp` is the independent
reference field.

**THE ACCEPTANCE, WHICH IS §3.3's CHECK 3 AND THE ONLY ONE THAT COULD CATCH A
SELF-CONSISTENT MISREADING.** The loop's flux is built from complete elliptic
integrals and shares no line of code with the Gegenbauer basis:

| `ρ_Γ` | modes | trace error | DtN, relative |
|---|---|---|---|
| 2.5 | 12 | 5.60e−07 | 1.95e−05 |
| 2.5 | 24 | 5.14e−12 | 3.41e−10 |
| 4.0 | 12 | 1.15e−09 | 7.12e−08 |
| **4.0** | **24** | **1.44e−15** | **6.89e−14** |

against §3.3's own 1.7e−11 / 2.3e−09 and 4.4e−14 / **6.8e−09**.

**THE DtN COLUMN IS FIVE ORDERS BETTER AND §3.3 PREDICTED WHY.** It says of its
own figure that it is *"limited by the finite-difference reference, not by the
method"*. It was right. This test compares against the loop's **analytic**
gradient, and the floor moves from 6.8e−09 to 6.9e−14. The old number was
measuring its instrument; this one measures the expansion. `CurrentLoop`'s
derivatives are themselves pinned against a Richardson-extrapolated difference
at 4e−11, which is what entitles them to be the reference.

**The other four checks, all green**: every mode vanishes on the axis at
**exactly 0.000e+00**; `h_n` matches the closed form to round-off at `n = 2..9`;
the projection is orthogonal to 1.2e−15; each exterior mode is `Δ*`-harmonic to
4.2e−08…6.5e−07 relative, **recomputed by central differences in `(r,z)`** so
that an error in the spherical separation could not hide; the symbol agrees with
`∂ψ/∂ρ` of its own field to 2.6e−08; and the modes decay as `2^{1−n}` exactly.

**THE SPECTRUM CONFIRMS THE COST ARGUMENT**, which is the other half of §3.3 and
is what makes `N` small. `|a_12|/|a_2|` for a unit loop:

| `ρ_Γ` | 2.0 | 2.5 | 4.0 | 8.0 |
|---|---|---|---|---|
| `\|a_12\|/\|a_2\|` | 2.6e−03 | 2.8e−04 | 2.6e−06 | **2.5e−09** |

Monotone, and geometric. §8's `ρ_Γ`-against-mesh trade is real and this is the
curve it is traded along — still unrun for a tokamak geometry, and §3.5 records
that CEDRES++ offers no prior art for choosing `ρ_Γ` at all.

**ONE IMPLEMENTATION DECISION THAT DEPARTS FROM THIS PLAN, DELIBERATELY.** §3.1
writes `C_n = ( P_{n−2} − P_n )/(2n−1)`, which is correct and is a bad way to
compute it: both Legendre values tend to 1 as `μ → ±1`, so it is a cancelling
difference exactly at the axis. Measured at `n = 8`, the ratio of the larger
operand to the result is 2.4 at `μ = 0.9`, 66 at 0.999 and **6666 at 0.99999**,
growing without bound. `ExteriorDtN.cpp` uses the identically equal

```
C_n( μ ) = ( 1 − μ² ) P'_{n−1}( μ ) / ( n( n − 1 ) )
```

— agreeing with the printed form to **1.4e−42** in 40-digit arithmetic — which
carries `(1 − μ²)` as an explicit factor. Two consequences, and both are load
bearing rather than cosmetic:

* **The axis is exactly zero rather than round-off.** Given §7.2 has just
  measured the axis to cost a power of `h` in the trace conditioning, a basis
  that quietly loses four digits there is not what to build on.
* **The singular weight never appears anywhere.** `C_n/(1 − μ²)` is a
  **polynomial**, so `dΓ/r = dμ/(1 − μ²)` cancels in closed form and a plain
  Gauss-Legendre rule integrates the mass integrals **exactly** — which the
  1.2e−15 orthogonality is the evidence for. The printed form would need `0/0`
  at the endpoints and a near-cancellation beside them.

`h_n = 2/(n(n−1)(2n−1))` then falls out of the standard `∫(1−μ²)[P'_m]²` integral
rather than being transcribed, so **the closed form agreeing with this plan is a
check rather than a restatement**.

**WHAT IS NOT DONE: §3.4's CEDRES++ agreement.** §3.5 now states their operator
exactly, read off the rendered page, and names three things that test must
respect — chiefly that `M` is **hypersingular, `~1/(π r d²)`**, so the
double-difference form *is* the regularisation and `∫∫ C_m M C_n` cannot be
assembled directly. That is the remaining piece of FB-0 and it is specified
rather than started.

### 7.4 FB-1's coupling matrix is built and measured. 2026-09-04.

`GradShafranovSolver::exteriorTraceColumns( ExteriorDtN const & )` returns the
columns of `P`; `tests/convergence/FreeBoundaryCoupling.cpp` is the acceptance.
This is **part** of FB-1 — the Dirichlet half of the coupling — and §7.5 below
says what remains.

**THE FOOT IS THE POINT, AND IT IS WORTH 18% ON A COARSE MESH.** `Γ_h` is the
inscribed polygon and `Γ` is the true boundary; a column must carry `C_n`
evaluated at the **foot** `a(x)` on `Γ`, not at `x` on `Γ_h`. Both choices give
a `P` that is supported in the right place, has independent columns, and couples
an exterior expansion to the solve — the wrong one is simply a different
function by `O(dist(Γ_h, Γ))`, which is `O(h)`, which is exactly the error the
transfer technique exists to remove. Measured, worst relative difference between
the two:

| `h` | 0.2125 | 0.1062 | 0.0531 |
|---|---|---|---|
| foot against `Γ_h` | **0.179** | **0.098** | **0.057** |

Falling at about `O(h)` — ratios 1.83 and 1.72 — exactly as `dist(Γ_h, Γ)` does.
`theColumnsAreTheModeAtTheFootAndNotOnGammaH` asserts both that they differ at
all (identical columns would mean the path is not being used) and that the
difference shrinks (a fixed difference would be a defect rather than the
geometry). **It is the test to read first if anything in FB-1 fails.**

`P` is also measured to live only on `Γ_h` — 135 of 1044 trace dofs at
`k = 2, n = 12` — and to have independent columns, Gram determinant 2.24e-01
with diagonals falling 1.51e+01 → 1.08e-01 across six modes, which is the modes'
own mass ordering. Dependent columns would give a singular corner block and a
Newton that cannot take a step.

**§4.3's REQUESTED MEASUREMENT IS NOT NEEDED UNDER NPC, AND THAT IS A RESULT
RATHER THAN A DODGE.** It asks for the column to be built at two well-separated
iterates and differenced, to test that `∂F/∂a` is constant. Under NPC it is
constant *by construction*, twice over:

* **`P` is a function of the geometry alone.** `exteriorTraceColumns()` takes no
  iterate and cannot — it reads the mesh, the path and the mode, and nothing
  else exists for it to read. The signature is the assertion; if it ever needs
  an iterate, the claim has failed.
* **`Γ_h`'s trace dofs are essential**, so the reduced operator masks its
  residual to zero there and puts a unit row in the Jacobian — which
  `theEssentialTraceConditionImposesTheDatum` already pins. The border column is
  therefore **exactly `−P`**, not a difference of one.

That is §4.5's *"two of the three borders are exact under NPC"* applying to this
border for the same reason it applies to `ψ_ax`'s. The differencing §4.3 asks
for is the right measurement for a **condensed** formulation, and MEQ is not one
any more.

### 7.4a FB-1's exact answer, and why the plan's proposal could not be one

`tests/analytic/ExteriorMatched.hpp`. FB-1's acceptance needs a problem whose
answer is known, and §7's *"the sum of the coils' loop fields"* **cannot be
one for an order study.** An exact loop field comes from a **filament**, which
in the `(r, z)` half-plane is a point source; its `ψ` has a logarithmic
singularity, so `ψ ∉ H¹` there, the finite element solution converges at a
reduced rate, and `k+1` is unreachable. A finite cross-section restores the
regularity but its exact field is a 2-D integral of loop fields, semi-analytic
and weakly singular inside the coil.

**So FB-1 gets a manufactured solution instead**: an exact exterior mode outside
a radius `ρ_0`, matched `C^1` to something regular inside, so the source is
compactly supported strictly inside `Γ` and the exterior expansion on any
`Γ` with `ρ_Γ > ρ_0` is **exactly the modes put in, with known coefficients**.

**AND THE OBVIOUS CONSTRUCTION HAS THE SAME DEFECT AS THE COIL, ARRIVED AT FROM
THE OTHER SIDE.** Matching value and slope with `A ρ^n + B ρ^{n+2}` works, and
leaves the SOURCE discontinuous at `ρ_0` — so `ψ ∈ H^{5/2−ε}` and the rate caps
near **2.5 whatever `k` is**. A `k = 3` study built on it would be measuring the
matching radius rather than the coupling. The fixture therefore carries a
**contact order `p`**, the order to which the source vanishes at `ρ_0`:

```
Δ*ψ = κ ρ^n ( ρ_0² − ρ² )^p C_n(μ)   inside,   0 outside
```

with `κ`, `A` and the `c_i` in closed form and `D = ∫_0^1 x^{2n}(1−x²)^p dx > 0`
so `κ` always exists. **Default `p = 4`** — `ψ ∈ C^5`, enough for `k = 1…4`.

**`p = 0` SHIPS AS THE CONTROL AND IT MEASURES THE DEFECT.** In the band around
`ρ_0` the `Δ*` residual at `p = 0` **sits at 4.81 at every `h` and does not
converge at all**, while the interior and exterior columns are untouched. That
is the discontinuity, isolated.

**Measured**, and the coefficient check is the one FB-1 will lean on:

| | |
|---|---|
| independent sympy rebuild, `Δ*ψ + f` | **exactly zero, symbolically** |
| `Δ*ψ = −f` by differences, `p = 4` | `O(h²)` to 1.7e−06 at `h = 1e−4`, the floor |
| source outside `ρ_0`, 861 points | **0 nonzero** |
| analytic gradients vs Richardson | 2.8e−11 |
| **exterior coefficients, single mode** | **0.000e+00 — bit exact**; 11 zero entries at 3.3e−15 |
| exterior coefficients, three modes | 2.0e−16 and 6.8e−15 at two radii |
| **the DtN symbol against the analytic `∂ψ/∂ρ`** | **6.3e−16** |
| `ψ(0, z)` on the axis | **literally 0.0**, and `∇ψ = (0,0)` |

Two design details worth keeping. `multiMode()` uses degrees **2, 3 and 5** —
skipping 4 deliberately, so a zero coefficient sits *between* two live ones and
an off-by-one in the mode indexing cannot hide. And both radii `1.5 ρ_0` and
`3 ρ_0` are checked, because the scaling is `a_n = α_n ρ_Γ^{1−n}` and **one
radius cannot see it**.

The axis result is a free confirmation of §7.3's evaluation choice: `ψ` and
`∇ψ` are *literally* zero on the axis because `C_n` carries `(1−μ)(1+μ)` as an
explicit factor, and that propagates through the fixture without being asked
for.

### 7.6 FB-2's current and its acceptance identity are built. 2026-09-04.

`src/meq/Coils.{hpp,cpp}` — `meq::Coil` and `meq::CoilSet`, **MFEM-free** so CI
can test them; `tests/unit/CoilsTests.cpp` is the acceptance. What is *not* here
is the solve to apply them to, which waits on FB-1's transmission row.

**THE SOURCE FACTOR IS DERIVED, NOT TRANSCRIBED.** From
`(curl B)_φ = −(1/r)Δ*ψ` in MEQ's own field convention `B = (−q_z, +q_r)`:

```
μ₀ j_φ = (curl B)_φ   ⟹   Δ*ψ = −μ₀ r j_φ   ⟹   F = μ₀ r j_φ
```

which is §5.4's `μ₀ r I/|Ω_c|` exactly, and agrees with `MHDSource`'s
`μ₀ r j_φ = μ₀r²p′ + gg′` — so the coil and plasma terms are in the same units
and simply add.

**AND IT IS CHECKED AGAINST ITS OWN FIELD, WITH THE NEIGHBOURING CONVENTIONS
REJECTED.** `f()` and `psi()` are independent statements of the same physics, so
recomputing `Δ*` of the field and comparing against `−f` is the one check that
can see a missing or extra `r` — which converges at full rate to the wrong
function. Measured at the coil centre: `Δ*_FD ψ = −2.094605e+01` against
`−F = −2.094395e+01`, **1.0e−04** at the difference floor, while the plausible
wrong answers `F/r = 1.047e+01` and `F·r = 4.189e+01` sit a **factor of two**
away. The test asserts the rejection as well as the agreement, which is what
turns "these two agree" into "these two agree and their neighbours do not".

**FB-2'S OWN ACCEPTANCE IDENTITY, ON THE EXACT FIELD:**

```
∮ (1/r) ∂ψ/∂n dl  =  −μ₀ · totalCurrent()
```

Measured **3.3e−11** relative, Richardson-extrapolated from a midpoint rule.
There is no discretisation in that identity, so establishing it here means that
when FB-2 checks it on a **solve**, any discrepancy is the solve.

**THE SIGN IS NEGATIVE AND IS ASSERTED SEPARATELY.** §7 predicts the coupling
sign will be got wrong at least once, and this is the cheapest place to pin it.
The same identity written as a counterclockwise circulation of `B` in `(r, z)`
comes out **positive**, because `φ̂ = ẑ × r̂` — so both signs are defensible
sentences about different quantities, which is exactly how a sign error survives
review.

**Two more measurements.** A shrinking coil approaches its filament at clean
second order — ratios 4.025, 4.006, 4.002, 4.001 over five halvings — which says
the finite cross-section is a *smoothing* of the filament rather than a
different object. And `f` is proportional to `r` inside the coil, not constant:
it is the current *density* that is uniform, and stating that correctly is the
difference between `μ₀ r j` and `μ₀ j`.

**ONE IMPLEMENTATION FINDING WORTH KEEPING.** The cross-section integral is
logarithmically singular when the field point is inside the coil, and the
obvious route fails there: `k = 2√(ar)/d` rounds to exactly 1.0 within ~1e−8 of
a filament, where `std::comp_ellint_1(1.0)` is **NaN** — and the quadrature must
evaluate there. Carlson's forms take `k′² = ((a−r)² + dz²)/d²`, which is the
squared source distance over `d²` and is formed with **no cancellation**, and
still work at `k′² = 1e−300`. With a panel split at the field point and cubic
grading the interior converges at about `4p` in the rule order, giving
**machine precision outside and ~3e−12 inside** at the default order 32. So the
model is *not* exterior-only, which it would have had to be otherwise.

**And it is sharper near the conductor than `tests/analytic/CurrentLoop.hpp`**,
which is worth recording because that fixture is FB-1's reference: the two agree
to **1.4e−12** over 4813 points, and walking in to the conductor `CurrentLoop`
returns NaN from `ε = 1e−9` while `Coils` tracks the line-current logarithm
exactly. Neither is wrong — `CurrentLoop` documents its own domain — but a
caller evaluating near a coil should use `Coils`.

**Deliberately not built**: no `meq::Source` adapter (`Source::f` takes `ψ` and
a coil current does not; `dFdPsi` would be identically zero, and the adapter
belongs with the assembly), and no `[coils]` TOML parsing. Both are driver work
and belong with FB-1's solve.

### 7.5 What FB-1 still needs, now that the quadrature over `Γ` has landed

**Done**: `meq::ExteriorDtN` (§7.3), the coupling matrix `P` (§7.4), — from
the two fixtures beside them — a verified current-loop field and a
manufactured exterior-matched solution to measure against, and, since
2026-09-05, `mfem::ExtensionBoundaryQuadrature`, which was the one piece that
was genuinely new code and is now upstream's.

**Not done, in the order they block each other:**

1. **The transmission row `T_m = ∫_Γ E_h(q_h)·ν C_m dΓ − a_m (1−n) h_m/ρ_Γ`.**
   This is the Neumann half of the coupling and without it `a` is undetermined.
   It is an integral over the **true `Γ`** of the element extension of the flux.

   **ITS MFEM HALF IS DONE, IS UPSTREAM'S, AND THE FILING PAID FOR ITSELF —
   2026-09-05.** This section used to be an instruction; it is now a record.
   `mfem::ExtensionBoundaryQuadrature` was written in MEQ against
   `gf-hdg-subdomains-dev`, used, given the tiling check below, and filed — and
   upstream merged it as *"Merge meq's ExtensionBoundaryQuadrature: quadrature
   over Gamma"*. It is a free function taking a visitor, the shape §6.2 argued
   for rather than an `Integrator` subclass, and it is in MEQ's install:

   ```
   void ExtensionBoundaryQuadrature(
      FaceElementTransformations &FTr, const TransferPath &path,
      const IntegrationRule &face_ir,
      const std::function<void( const ExtensionBoundaryPoint & )> &visit,
      real_t fd_step = 1e-6 );
   ```

   **AND THE TILING CHECK IS WHAT PAID.** The acceptance was that the boundary
   weights sum to `|Γ|`, as the volume weights sum to `|Ω| − |D_h|`, and writing
   it found that an **unsigned** weight integrates a swept image *with
   multiplicity*: where a boundary feature is thinner than a mesh width the foot
   map folds, and the fold adds where it should cancel. Upstream then found the
   sibling `ExtensionRegionQuadrature` had the same defect — `|det J|` — and
   signed it, taking the worst relative error on their aerofoil from `+1.13e-02`
   to `−2.29e-04`, fifty-fold. **MEQ's own geometry is a disc and does not fold,
   so MEQ's numbers are untouched**; the return on filing was to somebody else's
   test case, which is the argument for filing a thing that has been used.

   **What is still MEQ's** is the row itself: sweep `Γ` with that routine,
   evaluate `E_h(q_h)` at each point on it, contract against `C_m`, and subtract
   the diagonal exterior term. `ElementExtension::SetElement` + `TransformBack`
   evaluate the owning element's polynomial outside it — the pattern
   `meq::Sampler::extendOutward` already uses — and on a **semicircular** `Γ`
   the normal is `ρ̂` analytically, which removes the fiddliest part. It is
   assembly rather than research.

   **THE PREREQUISITE ARGUMENT STILL STANDS AND IS WHY THE TARGET WAS RIGHT.**
   The patch went to `gf-hdg-subdomains-dev` because that is where the extension
   machinery lives, and — verified rather than assumed, 2026-09-05 — that branch
   carries **no `DarcyNPCOperator` at all**. So the routine provably does not
   depend on the ordering, which is what made subdomains both the right target
   and a *sufficient* one.

   **Never commit to `meq-integration`.** `CLAUDE.md` is explicit that it exists
   only to be built against, is re-created whenever any of its four parents
   moves, and that anything committed to it is lost. The work belongs on a
   branch off subdomains — `hdg-extension-boundary-quadrature` is the one this
   went on — and reaches MEQ by fetch, never by editing the install's sources.
   **THE ROW ITSELF IS NOW WRITTEN TOO, 2026-09-05.**
   `GradShafranovSolver::exteriorTransmissionRows()` is the contraction, and
   `theTransmissionRowIsTheBoundaryIntegralItClaims` measures it against a
   closed form at **7.8e-16 to 1.4e-13**. The measure question this section left
   open — whether `T_m` wants `dΓ` or `dΓ/r` — is settled and it is `dΓ`: the
   exterior block is diagonal in `dΓ/r`, and `q` *is* `(1/r)∇̄ψ`, so testing
   `q·ν` in the plain measure already carries the radius. `dΓ/r` would divide by
   it twice.

   **What is NOT settled is the tiling.** The boundary sweep's own acceptance —
   weights summing to `|Γ|` — was exact and mesh-independent when the routine
   was written and now converges at `O(h²)` against the cone-carrying
   `VertexConePath`. That test is left red; §7.7 has the numbers.
2. **The bordered solve at `N + 2`.** §4.4 is right that this is mechanical:
   `solveWithNormalisation()` already does it at `N = 1` through a
   `DarcyNPCSolver`, and the generalisation replaces a scalar corner with a
   dense `(N+2)×(N+2)` one. `N + 2` extra backsolves, no extra factorisation.
3. **The acceptance itself** — `ψ_h` against the manufactured solution at
   `k+1`, and the coupling sign, which §7 says will be got wrong at least once.

**A DOMAIN CONSTRAINT THAT ONLY BECAME OBVIOUS ON BUILDING THIS.** The exterior
expansion is only valid on a **semicircle centred on the axis**, so FB-1's `Γ`
must be one — which means the domain reaches `r = 0` and everything §7.2
measured about the axis applies to it. FB-A was not merely a warm-up for FB-1;
it is its prerequisite, and the `O(1/h)` conditioning it found is the number
FB-1 inherits. `FreeBoundaryCoupling.cpp` deliberately uses
`ExtensionConvergence`'s Solov'ev surface instead, because `P` is a projection
and needs no semicircle — but nothing that SOLVES can take that shortcut.

### 7.7 The tiling was red, MEQ's diagnosis was wrong, and the mistake is the lesson

**Filed and answered the same day, 2026-09-05.** At a 12th-order face rule the
boundary weights summed to `|Γ|` only to 1.01e-04, 2.24e-05, 4.59e-06 at
`n = 12, 24, 48`, converging at about `O(h²)` where the expectation is a
mesh-independent floor near 1e-10. A one-variable control — the same `D_h` with
`VertexConePath`'s cone off — put it back at the floor, so the cone was
unambiguously the cause. MEQ filed that as **lost coverage**, arguing that *"a
quadrature residual that converges is measuring a geometry rather than an
instrument"*.

**THE ARGUMENT WAS WRONG AND THE HOLE IN IT IS WORTH REMEMBERING: TWO THINGS
CONVERGE.** The second is a curve the *rule* under-resolves, which straightens as
`h` falls. Refining the **rule** at fixed `h` separates them, because no
quadrature recovers coverage that is not there. Measured on MEQ's own circle at
`n = 12`, cone on:

| q8 | q12 | q20 | q40 | q80 |
|---|---|---|---|---|
| 2.61e-04 | 1.01e-04 | 1.34e-05 | 6.68e-08 | **5.40e-10** |

`5.40e-10` is the cone-off floor. **Coverage is exact either way.** What the cone
costs is the smoothness of `ξ ↦ a(x(ξ))` — it drives the two interpolated vertex
directions apart, the foot map roughens, and a fixed-order rule under-resolves
it.

**So the repair was the rule, not the gate.** The gate was never too tight. That
distinction is the whole finding: relaxing it would have made the suite green and
hidden a genuine under-resolution in MEQ's own transmission row, whose
`transmissionQuadratureOrder` was the same 12 and is now 40.

**And the acceptance beside it could not have caught it**, which is worth stating
because it reads 7.8e-16.
`theTransmissionRowIsTheBoundaryIntegralItClaims` builds its reference by
sweeping `Γ` with the **same rule the row uses**, so the quadrature error is
common to both sides and cancels exactly. That is deliberate — it is what
isolates the contraction, the vdof ordering, the sign and the measure — but it
makes the case blind to whether the rule resolves the foot map at all.

**What upstream did with it**, all of it useful: reproduced it on their own disc,
made the cone a `use_cone` constructor flag and turned it **off by default**,
added the boundary-sweep tiling case MEQ asked for, and corrected their own
commit message's *"changes nothing"*. They also measured that the cone does not
do what it was added for — the aerofoil flux order it targeted is
**pre-asymptotic** and recovers on its own — and connected it to MEQ's earlier
signed-weight finding: **the cone is what makes the foot map backtrack**, so the
`O(h)` unsigned overcount and this `O(h²)` residual are two readings of one
thing.

### 7.8 FB-1 is done, and the datum on `Γ` was the missing half

**FB-1a: the half-disc, datum given.** The first solve in this tree on a
semicircle centred on the axis — §7.5's domain constraint, which every other
extension study dodges by taking a Solov'ev surface away from the axis because a
projection needs no semicircle. Rates in `ψ` **1.971, 1.992 / 2.967, 2.991 /
3.961, 3.989** at `k = 1, 2, 3`; in `q` 1.94 / 2.67 / 3.93. So §8's corner is
benign, and `q`'s shortfall at `k = 2` is FB-A's half-order axis loss rather than
the corner.

**AND IT FOUND THAT MEQ COULD NOT IMPOSE A NON-ZERO DATUM ON `Γ` AT ALL**, which
is the one thing FB-1 could not do without. `setBoundaryData()` is projected
against `fittedMarker`; `HDGExtensionIntegrator` supplies only the
solution-dependent half of the transferred datum, which is the whole of it
exactly when `g` is homogeneous — and every extension study here wants `ψ = 0`.

**The first repair was inert, and why is the transferable part.** Setting
`Γ_h`'s trace dofs and letting `FormLinearSystem` eliminate them — which *is* how
the fitted datum is imposed — changed not one digit. The flux divergence form's
boundary face integrator carries `fittedMarker`, and `EnableHybridization`
registers a boundary flux constraint on exactly the attributes marked there;
`Γ_h` is not among them, so its trace dofs are essential in name with nothing
coupled to them.

**`miniapps/hdg/extension.cpp` is the worked example**, and it imposes the datum
as a **load on the flux equation** — `VectorBoundaryFluxLFIntegrator` with a
`PathTraceCoefficient`, which is `⟨ψ̂, v·n⟩` of (8a). It **negates**
(`pNatural = -pExact`). `setExteriorDatum()` does that now. **There was no MFEM
defect**, and a report was nearly filed: the capability exists and is exercised
upstream on all three of the miniapp's problems.

**FB-1b: `a` becomes an unknown**, recovered from the transmission condition
alone at `k = 2`:

| degree | solved | exact | relative |
|---|---|---|---|
| 2 | 6.667965e-01 | 6.666667e-01 | 1.95e-04 |
| 3 | −2.666651e-01 | −2.666667e-01 | 2.40e-06 |
| **4** | **−2.31e-05** | **0** | 3.46e-05 |
| 5 | 6.902000e-02 | 6.913580e-02 | 1.74e-04 |

**converging at 3.30** over `n = 12 → 24`. Degree 4 is the mode `multiMode()`
deliberately skips so a zero sits *between* two live ones; it comes back zero, so
the indexing is right and the off-by-TWO `ExteriorDtN` warns about is not there.

**THE RATE IS WHAT MAKES IT A RESULT.** A coupling wrong by a *constant* — a
sign, a stray `r`, a misindexed mode — sits at a fixed distance and looks like a
plausible discretisation error at any single mesh. Only refinement separates
them, which is why one mesh was not enough.

**What is deliberately NOT here.** The solve is closed by **superposition** in
the test rather than by a bordered solve in `GradShafranovSolver`. For a vacuum
problem that is exact rather than approximate — everything is affine in `a`, so
`x(a) = x_0 + Σ a_n x_n` and block elimination *is* one factorisation and `N+1`
back-substitutions. What the test throws away is the *reuse* of the
factorisation, which is a cost and not an answer. A plasma makes the interior
non-linear and then the border has to live inside Newton, which is **FB-5**, and
§4.4 is right that `solveWithNormalisation()` already does it at `N = 1`.

### 7.9 FB-2, and two things the mesh should be aligned to

**`outwardFlux()` is Ampère's law on the solve**, `∮_Γ q·ν dΓ = −μ₀ I`, with no
discretisation anywhere in the identity. `CoilsTests` already pinned it on the
exact field at 3.3e-11, so a discrepancy is the solve.

**MEQ CONSERVES CURRENT EXACTLY.** Integrated over the **mesh** boundary with no
extension, it is **round-off — 1e-14 to 6e-13** — at every degree and mesh on a
contour clear of the axis, and `k+1`-convergent to 5.6e-12 on the half-disc. So
the assembled operator, the source, the boundary condition and the trace solve
reproduce the enclosed current exactly.

**EVERYTHING LOST IS THE BAND.** Over the *true* `Γ`, through `E_h(q_h)`, the
residual is 1e-3 on the half-disc and floors flat in `h` and `k`. Away from the
axis it converges instead. The extension is an extrapolation and satisfies
`div q = 0` only to its own order.

**AND THE FLOOR IS THE GAP BETWEEN `Γ_h` AND `Γ` AT THE AXIS. CLOSING IT IS
WORTH 17×**, changing nothing but `ρ_Γ` so `D_h`'s topmost axis row is included
rather than excluded:

| `ρ_Γ` | axis gap | band | mesh boundary | `ψ` L2 |
|---|---|---|---|---|
| 1.5000 | 0.0125 | 1.07e-03 | 2.01e-08 | 1.6355e-06 |
| 1.5416 | 0.0010 | **6.41e-05** | 1.74e-08 | 1.6361e-06 |

The mesh-boundary residual and `ψ` barely move, which is what says the solve is
untouched. `D_h` is the union of elements *entirely* inside `Γ`, so the staircase
stops at the last mesh line whose outer corner still fits — putting `ρ_Γ` just
**above** a mesh line rather than just below includes that row. Making the gap
fall as `O(h²)` wants an offset of about `h²/(2ρ_Γ)` chosen per mesh.

**THIS BOUNDS FB-1'S TRANSMISSION ROW**, which is the reason to care: that row is
`∫_Γ E_h(q)·ν C_m dΓ` over this same contour, so a coupling needing better than
1e-3 near the axis wants the gap closed rather than the mesh refined.

**THE CONDUCTOR SHOULD BE MESH-ALIGNED TOO, AND THAT IS WORTH A WHOLE ORDER.**
`F = μ₀ r j` is discontinuous at a uniform-density conductor's edge, and where
that edge cuts a cell the element quadrature integrates a discontinuous
integrand with a rule assuming smoothness:

| `k` | cut cells | aligned |
|---|---|---|
| 1 | 1.330 | **1.991** |
| 2 | 1.265 | **2.876** |
| 3 | 1.086 | **3.013** |

with the `k = 3` L2 falling from 1.08e-04 to **1.97e-08**. A rate that *falls*
with `k` is the signature — a genuine regularity limit is flat in `k`, and only
a quadrature error worsens as the rest of the scheme improves.

**AND THAT IS A SCOPE REDUCTION FOR FB-4.** The two look like one problem and are
not: a conductor's geometry is **prescribed input** and can always be meshed to,
while `χ_{Ω_p}` is bounded by a boundary that **moves with the solution** and
cannot be. Cut quadrature is needed for the plasma support and **not** for the
conductors, so §6.4's one real gap is narrower than it looked.

The aligned rate caps at 3 because a rectangular conductor has **corners**, and a
corner in the forcing gives the same `r² log r` behaviour a corner in the domain
does. Alignment cannot fix that; rounding the conductor would.

### 7.10 FB-4: the plasma edge caps the order, and the profile sets the cap

**Measured 2026-09-05. `tests/analytic/PlasmaEdge.hpp`,
`tests/convergence/PlasmaEdgeConvergence.cpp`.** The full account is in
`CLAUDE.md`; this is what changes about the plan.

**The fixture pair.** `PlasmaEdge` is a `Δ*`-harmonic vacuum field plus
`c(φ_+)^m` on a prescribed circle, so `F` is supported in the plasma and nowhere
else and the cut is fixed. `MovingPlasmaEdge` makes the plasma exactly `{ψ > 0}`
of its own solution, at the price of a smooth background source outside it —
unavoidable, because `Δ*` obeys a maximum principle and a compactly contained
`{w > w₀}` cannot exist for harmonic `w`. `m = j + 2` where `j` is the order to
which the profiles vanish at the edge.

**Three measurements, and the first needs no solver.**

1. **The L2 best approximation** by `P_k` and `P_{k+1}` caps `ψ*` at
   `min(k+2, j+2.5)` for *any* method, with the composite rule on cut elements
   so that approximation and quadrature are separated. The uncut elements read
   `k+2` throughout — 2.96, 3.95, 4.99 — so the loss is a set of measure `O(h)`.
2. **The solve**, over `j = 0…3` and `k = 1…4`, crosses `k+2` exactly at
   `k ≤ j`. `j = 3, k = 3` reads **4.989**.
3. **The rule swept at fixed geometry**: `ψ_h`'s rate is pinned to three
   figures across `extra = 4…20`, which is what says the loss is blindness to a
   kink rather than too few points.

**What this does to §5.3, §6.4 and §9.** §5.3's two options — ignore the
cut-rule derivative or difference it — do not arise. §6.4's "one real gap" is
created by adopting a cut rule rather than closed by it. And **FB-4 is not the
obstacle to `k+2` that §8 lists it as**; what is, is `j`.

**The one thing that does bite, and it is worse than an order.** At `j = 0` the
moving support does not converge at all — not from the exact solution, not under
`PicardThenNewton`. `∂F/∂ψ` carries `F·δ(Ψ)` there and `meq::Source` cannot.
**`j ≥ 1` is a precondition of the free-boundary path**, and FB-6's freegs4e
cases must be checked against it: FreeGS's `(1 − Ψ_n^α)^β` gives `j = β`, and
`β = 1` is its default, which is enough for `k = 1` and not for more.

**A consequence for FB-5 and for the driver.** A configuration whose profiles do
not vanish at the edge should be **refused at parse**, the way
`[boundary] Type = "exact"` is, rather than run to a non-convergence the user
has to diagnose. That is driver work and belongs with FB-5.

### 7.11 FB-6's test problem: there is no reproducible ITER case, and there is a better answer

**Searched 2026-09-06.** The question was whether `refs/` carries an ITER
free-boundary benchmark, ideally CEDRES++'s. It does not — but the search
turned up something more useful, and it turns on a symbol clash between two
papers.

**CEDRES++ §4.1 IS AN ITER CASE AND IS NOT REPRODUCIBLE.** It specifies
`I_p = 15.10 MA`, `r₀ = 6.2 m` and the profile parameters, and then:

* **no coil currents.** Its Problem 2 takes them as *given* and §4.1 never gives
  them. Without them there is no free-boundary problem to pose.
* **no wall, limiter or vessel geometry.** "ITER geometry" is all that is said.
* **no reference field.** Their `ψ_ref` is *their own* fine-mesh answer on
  1,153,174 triangles, shown as a figure (their Fig. 3). §4.1 is a
  self-convergence study, which is precisely the fallback FB-6 rejects: it
  shares every convention with the code it checks.

`refs/MFEM-GS-Newton.pdf` (Serino, Tang, Tang, Kolev & Lipnikov) has a **15 MA
ITER baseline** and is MFEM-based, so it looked promising — but it solves the
INVERSE problem, seeded from a proprietary ITER discharge carrying the reference
number **ABT4ZL**. The coil currents are an output of a run nobody outside can
reproduce. Nothing else on disk carries a free-boundary case with data.

**WHAT THE SEARCH DID FIND: CEDRES++'s PROFILE MODEL AND `freegs4e`'s ARE THE
SAME FAMILY, TERM FOR TERM.** CEDRES++ (2.11), off the rendered page:

```
S_p'( Ψ ) = ( β/r₀ )( 1 − Ψ^α )^γ ,     S_ff'( Ψ ) = ( 1 − β ) μ₀ r₀ ( 1 − Ψ^α )^γ
```

and `freegs4e`'s `ConstrainBetapIp`:

```
J_φ = L [ β₀ R/R_axis + ( 1 − β₀ ) R_axis/R ] ( 1 − Ψ^{alpha_m} )^{alpha_n}
```

Since `μ₀ R J_φ = μ₀ R² p' + ff'`, these are the same two-parameter family with

| CEDRES++ | `freegs4e` | ITER value |
|---|---|---|
| `α`, the peakage exponent | `alpha_m` | **2** |
| `γ` | `alpha_n` | **1.395** |
| `β` | `beta0` | **0.5978** |
| `r₀` | `Raxis` | **6.2 m** |
| `λ`, fixed by `I_p` | `L`, fixed by `Ip` | `I_p = 15.10 MA` |

**So the published ITER PROFILE is exactly runnable in `freegs4e`**, and FB-6 can
use it on a machine `freegs4e` defines — which gives a published, ITER-relevant
current profile *and* an independent reference, where CEDRES++ offers the first
without the second.

**AND `α = 2` RATHER THAN THE `0.5978` CEDRES++ PRINTS.** The paper's §4.1 says
`α = 0.5978, β = 0.5978, γ = 1.395` — the same number twice. Serino et al., using
*"the same coefficients used in [11] for the ITER configuration"*, write their
(2.5) with a different symbol assignment — their `α` is the scaling constant
CEDRES++ calls `λ`, and their `δ` is CEDRES++'s peakage `α` — and set **`δ = 2`**,
`β = 0.5978`, `γ = 1.395`, `r₀ = 6.2`. **CEDRES++'s printed `α` is a repetition of
`β`.** Both readings were taken off rendered pages, not from `pdftotext`, per the
standing rule for this pair of papers. `α = 2` also gives the standard
`( 1 − Ψ² )^γ` shape and is `freegs4e`'s own convention.

**AND FB-4's PRECONDITION BITES ON THE PUBLISHED ITER PROFILE.** The vanishing
order at the plasma edge is `j = γ` exactly, since `( 1 − Ψ^α )^γ ~ ( α( 1 − Ψ ) )^γ`
there. With ITER's `γ = 1.395`:

* `ψ*` keeps `k+2` only where `k + 2 ≤ j + 2.5 = 3.895`, i.e. **`k = 1` and no
  higher**;
* and `γ < 2`, so `∂F/∂ψ ~ ( 1 − Ψ )^{γ−2}` is **unbounded at the edge** — worse
  than §7.10's `j = 1` fixture, where it merely jumps.

**So the ITER profile as published is a `k = 1` case with a singular Jacobian at
the plasma boundary.** That makes it a good late stress test and a bad first
one. **FB-6 should open at `γ = 2` or 3**, which is one number in the same
family, and bring `γ = 1.395` in afterwards as the published case — with the
expectation, from §7.10, that `k ≥ 2` will not hold `k+2` on it and that this is
the profile's property rather than MEQ's.

**The machine.** `freegs4e` defines `TestTokamak`, `DIIID`, `MAST`, `MAST_sym`,
`TCV`, `MASTU_simple` and `MASTU` — **and no ITER**. Any of them serves, and
`tools/freegs4e-benchmark/` already drives seven configurations across them for
the fixed-boundary rehearsal. Building an ITER coil set from public design data
is possible and would not be *the* CEDRES++ case anyway, its currents being
unpublished; it is worth doing only if an ITER-scale aspect ratio is wanted for
its own sake.

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

**~~The corner where `Γ` meets the axis.~~ — SETTLED 2026-09-05, AND IT IS
BENIGN.** `theSolverReachesTheExteriorDatumOnTheHalfDisc` solves on the
half-disc with the exterior datum given and reads **1.99 / 2.99 / 3.99** in `ψ`
at `k = 1, 2, 3` — full `k+1` across the corner. `q` reads 1.94 / 2.67 / 3.93,
and the `k = 2` shortfall is FB-A's half-order axis loss rather than the corner.
The reasoning that made it a risk is kept below because it was sound and only the
conclusion was unknown.

Two right-angle junctions, and
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

**`N`, `ρ_Γ` and the mesh — AND THE CAVEAT IS ANSWERED, 2026-09-06.** §3.3 shows
the trade-off exists and is geometric, and the coupling paper's own worry is that
*"the introduction of a circular interface may require a large computational
domain in situations where the support of source terms is very elongated"*.

**That is a worry about DEGREES OF FREEDOM and not about geometry, and MEQ
already has the answer to it.** A large `ρ_Γ` costs nothing if the mesh out
there is coarse: take a **coarse background covering the semicircle and every
coil**, and let the adaptive loop refine aggressively where the plasma is. The
vacuum between the plasma and `Γ` is source-free and `Δ*`-harmonic — the
smoothest thing in the problem — so it is exactly where a coarse mesh is
cheapest, and `meq::AdaptiveDomain` and the residual estimator are what put the
elements where they earn their place.

So the choice of `ρ_Γ` is decoupled from the cost of an elongated plasma, and
what is left to measure is the **spectrum**: how `N` must grow with `ρ_Γ` and
with the coil set, which §3.3 measures on a disc and nobody has measured on a
machine. That is FB-6's, and it is a sweep rather than a design question.

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
