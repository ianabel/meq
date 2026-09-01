# Toroidal flow: the generalised Grad–Shafranov equation of RoPP (136)

Written 2026-09-01. `CLAUDE.md` is the operational record and is authoritative on
anything technical; `ROADMAP.md` is the order of work; this file is the design for
**sonic toroidal rotation**, from the two-species case to `n` species.

The equation is `refs/RotatingGK.pdf` eq (136) — Abel, Plunk, Wang, Barnes,
Cowley, Dorland & Schekochihin, *Multiscale gyrokinetics for rotating tokamak
plasmas*, Rep. Prog. Phys. **76** (2013) 116201 — closed by its (96) and (97),
which determine the poloidal density variation and the electrostatic potential
`φ₀` that holds quasineutrality against it.

**Inputs are `T_s(ψ)` per species, `ω(ψ)`, `n − 1` density flux functions, and
`I(ψ)`.** `I(ψ)` is an input here, exactly as `g g′` is today. A fixed-`q(ψ)`
solver — where `q(ψ)` is given and `I(ψ)` is found from it — is a separate
extension and is §9.

This plan supersedes `TODO`'s *Sonic toroidal rotation* entry, which was right
about the gauge and about the root find and is folded in below.

---

## 1. What changes, and what does not

**What does not change is most of meq.** The operator is still `Δ*`, the
discretisation is still LDG-H on `DarcyForm`, `τ` is still constant, the
hybridization is untouched, the estimator and the adaptive loop are untouched,
and the curved-boundary extension is untouched. Rotation is a change to `F`
alone.

**And `Source` already has the right signature.** `f( r, z, ψ )` and
`dFdPsi( r, z, ψ )` carry `r`, so a source with genuine `r`-dependence at fixed
`ψ` needs no interface change. `MHDSource` simply happens to use `r` only through
the `μ₀ r²` prefactor. That is the single largest piece of luck in this whole
plan and it is worth saying out loud, because the obvious fear — that a
non-flux-function density forces an interface change through the solver — is
unfounded.

**What changes is four things:**

| | |
|---|---|
| the density is not a flux function | `n_s = n_s(r, ψ)`, by (96) |
| a new unknown appears in the source | `φ₀(r, ψ)`, closed pointwise by (97) |
| the profile set grows | per-species `T_s`, `n_s0`, plus `ω`; and the source can no longer store pre-multiplied products |
| `∂F/∂ψ` acquires a chain rule | through `ω`, `T_s`, `n_s0` **and** `φ₀` |

The fourth is the dangerous one, for the reason `CLAUDE.md` records under *A
wrong Jacobian is invisible to a convergence table*: every error and every rate
survives a wrong `∂F/∂ψ` unchanged, and only Newton's observed *order* moves.

---

## 2. The equation, in meq's convention

### 2.1 The three equations

meq solves `−∇̄·( (1/r) ∇̄ψ ) = F/r`, i.e. `−Δ*ψ = F`. RoPP is Gaussian; in SI,
with `g ≡ r B_φ` for its `I` and `μ₀` for its `4π`:

```
F(r,z,ψ) = μ₀ r² Σ_s n_s { T_s (ln N_s)′ + [ Z_s e φ₀ − ½ m_s ω² r² + T_s ] (ln T_s)′ }
         + μ₀ r⁴ ω ω′ Σ_s m_s n_s
         + g g′
```

closed by

```
(96)   n_s(r,ψ) = N_s(ψ) exp[ m_s ω²(ψ) r² / 2T_s(ψ) − Z_s e φ₀ / T_s(ψ) ]
(97)   Σ_s Z_s n_s(r,ψ) = 0                      ← determines φ₀(r, ψ)
```

`N_s`, `T_s`, `ω` and `g` are flux functions; `n_s` and `φ₀` are not.

### 2.2 The brace collapses, and that is the form to implement

**(136) is `∂p/∂ψ` at fixed `r`, written out.** Differentiating
`p = Σ_s n_s T_s` at fixed `r`, the `∂φ₀/∂ψ` terms collect into
`−e (∂φ₀/∂ψ) Σ_s Z_s n_s`, which vanishes **identically by (97)**. What is left
is exactly the brace of (136) plus its `ω ω′` term. So

```
F(r,z,ψ) = μ₀ r² ∂p/∂ψ|_r + g g′,        p(r,ψ) = Σ_s n_s(r,ψ) T_s(ψ)
```

which is `MHDSource`'s shape with an `r`-dependent `p`. Three consequences:

* **The residual needs `φ₀` but never `∂φ₀/∂ψ`.** The Jacobian does; see §4.3.
* **The implementation is `p`, then differentiate**, not the brace term by term.
  Fewer places to drop a factor, and the brace becomes a *check* on the
  derivative rather than the thing being coded.
* It is the same cancellation that makes the gauge free (§3.1) — both are
  `Σ_s Z_s n_s = 0` doing the work.

**Two independent confirmations, and this is a derivation rather than something
the paper writes out, so it wants re-checking before anything rests on it.** It
is the paper's own force balance (128), `(1/c) j×B = ∇p − ρ (ω²/2) ∇(r²)`,
projected on `∇ψ`. And at `ω → 0` it gives `F → μ₀ r² Σ_s p_s′ + g g′`, which is
the paper's low-Mach result (243) and **is meq's current equation**.

### 2.3 Conventions, pinned before any code

`TODO` flagged three and named none as settled. Two are settled here; the third
is a genuine hazard.

| | |
|---|---|
| **units** | RoPP is Gaussian (`4π` in (136), `c` in (135)). meq is SI. `4π → μ₀`, `I → g` |
| **the sign of `Δ*`** | RoPP (135) is `j·∇φ = −(c/4πr²) Δ*ψ`, so (136) reads `Δ*ψ = −4πr²{…}`, so `F = −Δ*ψ = +μ₀r²{…}`. **Positive**, and matching `MHDSource`'s `F = μ₀ r² p′ + g g′` |
| **the `2π` in `ψ`** | **checked and consistent.** RoPP (31) is `B = I∇φ + ∇ψ×∇φ`, so `B_φ = I/r` and `I = r B_φ = g`; `ψ` is poloidal flux **per radian**, which is what meq's `g g′` in `T²m²` per Wb/rad already assumes, and what EQDSK tabulates as `FF′` |
| **`Δ*` itself** | RoPP (125), Li & Zhu (7) and meq all agree: `∂²/∂r² − (1/r)∂/∂r + ∂²/∂z²` |

**Write the conversion down in the header of the new source and assert it in a
unit test.** Two sign errors have already been found in the papers meq *does*
follow, and neither was visible in a convergence rate.

---

## 3. The gauge is free, and meq takes the local one

### 3.1 The freedom is exact

RoPP fixes `φ₀` by `⟨φ₀⟩_ψ = 0` and says at (59) why: *"There is some
arbitrariness in the definition of φ₀ … as we can add any function of ψ to it.
We resolve this by requiring…"*. It is a convention, not a closure.

The transformation is, for any flux function `δ(ψ)`,

```
φ₀ → φ₀ + δ(ψ),        N_s → N_s exp( Z_s e δ / T_s )
```

which leaves `n_s` of (96) unchanged, hence leaves (97) satisfied, hence leaves
`p` and `F` unchanged. (`TODO` verified the `F` invariance symbolically with
sympy for two species; §2.2 is the same statement, since both reduce to
`Σ_s Z_s n_s = 0`.)

### 3.2 meq takes `φ₀(r_ref, ψ) = 0`

**`r_ref` is a constant** — the geometric axis `R₀`, given in the configuration —
**not the magnetic axis and not a flux-surface average.** The consequences:

* `F` is a pointwise function of `(r, z, ψ)`. No flux-surface averaging, no
  contour tracing on level sets of `ψ_h`, no second non-local closure.
* `N_s(ψ)` becomes the **physical density of species `s` on the curve
  `r = r_ref`**, which is a quantity a user can state and another code can be
  compared against.
* Li & Zhu do exactly this — their exponent is `R²/R₀² − 1`, referenced to the
  axis — which is independent support that it is the workable choice.

### 3.3 What the paper's gauge would have cost, and why it is not taken

`⟨φ₀⟩_ψ = 0` makes `φ₀ ∝ r² − ⟨r²⟩_ψ` (RoPP footnote 33 says so directly), and
`⟨r²⟩_ψ` is a flux-surface average of the unknown. That is a non-local closure on
**every** surface, not one scalar like `ψ_ax`, and meq has no flux-surface-average
machinery at all — `FluxSurfaces` was a `v0-legacy` driver and was not ported.
It would be a larger piece of work than the rest of this plan put together, and
it buys nothing physical.

**The one place it can still bite is documentation.** Two sets of `N_s` differing
by the gauge factor describe the same plasma. Anyone comparing profiles against
GS2, against a transport code, or against RoPP's own notation needs to know which
convention each is in. That belongs beside the profile documentation and in the
example TOML, not in a comment in the solver.

### 3.4 The counting: `n − 1` density flux functions

The gauge is one function's worth of freedom per `ψ`. Fixing it removes exactly
one from the set `{ N_s }`, so for `n` species there are **`n − 1` independent
density flux functions**. For two species that is one, which is the familiar
statement.

The gauge-free way to say the same thing, and the one to expose in the
configuration: **give the physical density of each species on `r = r_ref`,
subject to charge neutrality there**,

```
Σ_s Z_s n_s0(ψ) = 0
```

which is `n` profiles minus one constraint. In the local gauge `n_s0 ≡ N_s`, so
this is not a reparametrisation — it is the same numbers with a physical name.
The configuration layer should **solve** the constraint for one nominated species
rather than ask the user to satisfy it, and should reject a set that cannot
(all-positive charges, say).

---

## 4. `φ₀`: the pointwise closure

### 4.1 Two species: closed form, no root find at all

For electrons and one ion species, (97) is linear in `φ₀` after taking
logarithms. In the local gauge, and writing `T_eff ≡ T_i + Z_i T_e`:

```
e φ₀(r,ψ) = ( ω²/2 )( m_i/T_i − m_e/T_e )( r² − r_ref² ) / ( Z_i/T_i + 1/T_e )

n_i(r,ψ)  = n_i0(ψ) exp[ m_i ω²(ψ) ( r² − r_ref² ) / 2 T_eff(ψ) ]      (m_e → 0)
n_e       = Z_i n_i
p(r,ψ)    = P₀(ψ) exp[ m_i ω²(ψ) ( r² − r_ref² ) / 2 T_eff(ψ) ],   P₀ = n_i0 T_eff
```

Both species carry the **same** exponential, which is the classical result, and
`p` is Li & Zhu's (8) with `Z_i = 1`, `T = T_i + T_e`, `r_ref = R₀`.

`∂p/∂ψ` and `∂²p/∂ψ²` are then explicit — one exponential and its chain rule
through `P₀`, `ω` and `T_eff`. **The two-species production path needs no root
find, no implicit differentiation, and no inner tolerance.** Keep `m_e/T_e`
rather than dropping it: it costs one term and removes a question.

### 4.2 `n` species: a safeguarded scalar Newton with a unique bracketed root

For three or more species (97) is transcendental and is solved per evaluation
point. It is as well behaved as such a thing gets:

```
∂/∂φ₀ Σ_s Z_s N_s exp[…] = −e Σ_s ( Z_s²/T_s ) N_s exp[…]   <  0   strictly
```

— every term carries `Z_s²`, so the left-hand side is **strictly decreasing** for
any positive `N_s`, `T_s`. With at least one positive and one negative charge it
runs from `+∞` to `−∞`. So the root exists, is unique, and can be bracketed;
a safeguarded Newton or Brent cannot fail.

Three practical points:

* **Bracket from the two-species formula.** §4.1 applied to the dominant ion and
  the electrons gives a starting guess that is exact when `n = 2` and close when
  the impurity fraction is small.
* **Work in `e φ₀ / T_ref`**, not in volts. The exponents are `O(M²)` and the
  arithmetic should be too.
* **Do not cache across Newton steps.** The previous iterate at the same
  quadrature point is a good guess, but `CLAUDE.md`'s standing warning applies:
  anything that makes an evaluation depend on history stops it being a function,
  and the assembled Jacobian is then differencing something that is not one.

### 4.3 `∂φ₀/∂ψ` by implicit differentiation — never by differencing

The residual does not need it (§2.2); **the Jacobian does**, through `∂n_s/∂ψ`
and through the explicit `Z_s e φ₀` in the bracket. Write
`G(φ₀; r, ψ) ≡ Σ_s Z_s n_s = 0`; then

```
∂φ₀/∂ψ = − ( ∂G/∂ψ |_φ₀ ) / ( ∂G/∂φ₀ )
```

with `∂G/∂φ₀` the strictly negative sum above — already computed by the root
find — and `∂G/∂ψ|_φ₀ = Σ_s Z_s n_s [ (ln N_s)′ + ∂A_s/∂ψ|_φ₀ ]` where
`A_s` is (96)'s exponent. Both closed form.

**Differencing it would be the classic failure mode this project has already
paid for twice**: an inner solve differenced from outside gives a derivative
whose accuracy is the inner tolerance, and Newton degrades from quadratic to
linear with no wrong answer and no failing test. See `CLAUDE.md`, *The trap that
cost the most*.

---

## 5. The new pieces in meq

### 5.1 `meq::Species` — a value type, no MFEM

```
struct Species
{
    double mass;                                  // kg
    double charge;                                // Z_s e, or Z_s; pick one and say so
    std::shared_ptr<Profile const> temperature;   // T_s(psi), energy units
    std::shared_ptr<Profile const> density;       // n_s0(psi) on r = r_ref
};
```

`src/meq` deliberately keeps MFEM out of `Profiles` and `Source`; this stays on
the right side of that line and is unit-testable without the library.

### 5.2 `meq::RotatingSource`, and why it is not a generalised `MHDSource`

A new class, with `MHDSource` kept as the `ω = 0` path. The reason is in
`Source.hpp`'s own documentation: `MHDSource` stores the **products** `p′` and
`g g′` precisely so that `dFdPsi` is one `prime()` call per profile and no chain
rule exists to get wrong. **That trick is unavailable with rotation** — `p`
depends on `ω`, `T_s` and `n_s0` separately and on `r`, so the chain rule is
unavoidable. Two classes with different invariants, not one class with a flag.

Holds: `std::vector<Species>`, `ω(ψ)`, `g g′(ψ)`, `r_ref`, and a mode for the
two-species closed form versus the general root find. Evaluates `p`, then
`∂p/∂ψ`, then `F`.

### 5.3 `dFdPsi` is where the risk is

`∂F/∂ψ = μ₀ r² ∂²p/∂ψ² + (g g′)′`, and `∂²p/∂ψ²` runs the chain rule twice
through everything in §5.1 plus `φ₀`. **`SourceTests`' finite-difference check on
`dFdPsi` stops being box-ticking and becomes the load-bearing test of this
plan.** It must be swept over Mach number, not evaluated at one point: the terms
it is checking are the ones that vanish at `ω = 0`.

### 5.4 Normalised flux, and the bordered Newton

If `T_s`, `n_s0` and `ω` are given against normalised flux — which is how
`refs/GourdainContour.pdf` §V poses profiles and how every equilibrium code does
— then this must be a `meq::NormalisedSource` and `ψ_ax` is an unknown of the
system. **That machinery already exists and is green** (`HighBetaConvergence`),
and under NPC two of the three border quantities are exact rather than
differenced. So `meq::NormalisedRotatingSource` is a subclass away, not a
project.

`ψ_bnd` stays zero: this is still the fixed-boundary problem.

### 5.5 `Config` and `SourceFactory`

A new `SourceType::Rotating` and `RotatingParameters`, a `[[source.species]]`
array of tables, `[source] Omega`, `[source] ReferenceRadius`, and the
charge-neutrality solve of §3.4.

**`meq::NormalisedMHDSource` is still not reachable from a TOML file**, and that
gap is a prerequisite rather than a parallel task: the normalised plumbing gets
written once and both sources use it. Do that first; it is `ROADMAP.md` item 1.

Read `Config.cpp`'s `asFloat()` before adding keys — toml11's `find_or<double>`
silently returns the default for an integer node, which is how `Omega = 0` would
become `Omega = <default>` with no error.

### 5.6 Output: `n_s` is a field now

`n_s(r,z)` and `φ₀(r,z)` are genuine two-dimensional fields, not profiles, and a
rotating equilibrium is not interpretable without them. They belong in the
NetCDF alongside `ψ` and `B`, and they are cheap — both are algebraic in
`(r, ψ_h)`, so `Sampler` already has everything it needs.

---

## 6. Tests: what is a real benchmark, and what only looks like one

**The literature on rotating equilibria solves at least three different
equations, and they are easy to mistake for each other.** All three produce a
`Δ*ψ = −μ₀ r² ∂p/∂ψ|_r − g g′` with `p` varying on the surface, and they differ
only in the closure that fixes that variation:

| closure | `p(r,ψ)` | is it (136)? |
|---|---|---|
| **isothermal on a flux surface**, `T_s = T_s(ψ)` | `P₀(ψ) exp[ m_i ω² r² / 2T_eff ]` | **yes** — this is what (96)+(97) give |
| **adiabatic / polytropic**, `S(ψ) = p/ρ^γ` | algebraic in `r²`, carries `γ` | **no** |
| **constant density on a surface**, `ρ = ρ(ψ)` | linear in `r²` | **no** |

So:

### 6.1 The `ω → 0` limit — free, and first

`F` must collapse to `μ₀ r² Σ_s (n_s0 T_s)′ + g g′`, which is `MHDSource`. Run
the whole existing Solov'ev study through `RotatingSource` at `ω ≡ 0` and require
agreement with `MHDSource` to round-off and the published rates unchanged to
every digit printed. **This exercises the entire new profile, `φ₀` and chain-rule
machinery against an answer that is already known**, which is what `TODO` called
the cheap first milestone and it is right.

Not bit-for-bit: `Σ_s (n_s0 T_s)′` and a single tabulated `p′` are different
floating-point expressions. Round-off is the honest ask.

### 6.2 Li & Zhu's rotating Solov'ev — the same equation, and linear

`refs/SpectralElementGSRotation.pdf` §3.1 (Li & Zhu, *CPC* **260** (2021)
107264) gives a closed-form solution, their (14)–(15), of
`Δ*ψ = −p₁ r² exp[ M₀²( r²/R₀² − 1 ) ] − F₀` with `μ₀P₀′ = p₁`, `FF′ = F₀`,
`T = T₀` and `Ω = Ω₀` all constant. Their (8) **is** the isothermal closure with
`T = T_i + T_e`, so this coincides with (136) for two species sharing a
temperature flux function — *verify that reduction before using it*, in
particular that the per-species `(ln T_s)′` terms cancel as they must when
`T_i/T_e` is constant.

Their (16) is its `M₀ → 0` limit and is the static Solov'ev particular solution,
so §6.1 and §6.2 are the same test at two Mach numbers.

**But its source is constant in `ψ`, so `∂F/∂ψ = 0`.** It tests the
discretisation, the `r`-dependence and the `φ₀` machinery. It says nothing
whatever about the Jacobian — the same blind spot meq's existing `Soloviev.hpp`
rung has, and for the same reason.

**Check it with `deltaStarFD()` before trusting it.** That is the house pattern
and it is exactly the guard against the hazard in this section: recompute `Δ*ψ`
from their published closed form by central differences and assert it against the
`F` the solver is actually fed. If the closure differs anywhere, that goes red
immediately instead of converging beautifully to somebody else's equilibrium.

### 6.3 What is *not* a test of (136)

* **Maschke–Perrin** (Li & Zhu §3.2, their (18)–(21)). Their (20) carries `γ`,
  the ratio of specific heats: an **adiabatic** closure, not (136)'s isothermal
  one. Do not use it without re-deriving what equation it solves.
* **FLOW** (`refs/GuazzottoFLOW.pdf`, Guazzotto *et al*, *Phys. Plasmas* **11**
  (2004) 604). Arbitrary flow, poloidal **and** toroidal. With poloidal flow the
  equation changes type across the poloidal sonic surface; it is a different and
  harder problem. At best a cross-check in a pure-toroidal isothermal limit, if
  it has one.

### 6.4 A manufactured nonlinear rotating case — the one that tests the Jacobian

Neither published solution has `∂F/∂ψ ≠ 0`, so one is needed, in the shape of
`ManufacturedNonlinear.hpp`: choose `ψ`, choose profiles with real `ψ`-dependence
in all of `n_s0`, `T_s` and `ω`, and build the datum to fit. Acceptance is what
`NewtonConvergence` already asserts — observed Newton order 2, and **the
assembled Jacobian against a central difference of the assembled residual**,
which is a stronger check than the one on `dFdPsi` alone.

### 6.5 The `n`-species checks

* At `n = 2` the root find must reproduce the closed form of §4.1 to round-off.
  Same source, two code paths, one answer.
* `Σ_s Z_s n_s = 0` to round-off at every quadrature point, swept over Mach
  number and over impurity fraction.
* `∂φ₀/∂ψ` against a central difference of `φ₀` — the *inner* analogue of §5.3,
  and the thing that catches an implicit-differentiation slip.
* A three-species case (D, C⁶⁺, e) with the heavy impurity centrifugally
  enriched on the outboard side, which is the physics the whole exercise is for.

---

## 7. The staged plan

Each stage ends at a measured number.

| | | acceptance |
|---|---|---|
| **FL-0** | `meq::Species`, the profile container, the charge-neutrality solve. No solver, no MFEM. | Every profile's `prime()` against a finite difference; a neutrality set that cannot be solved is rejected with a named error |
| **FL-1** | `φ₀` and `p(r,ψ)` for two species, closed form. Still no solver. | `Σ_s Z_s n_s = 0` to round-off over a Mach sweep; `φ₀(r_ref) = 0` exactly; `p` against Li & Zhu (8) |
| **FL-2** | `meq::RotatingSource`, and the `ω → 0` collapse. | §6.1: agreement with `MHDSource` to round-off, Solov'ev rates unchanged to every digit |
| **FL-3** | `dFdPsi`, two species. | §5.3: central difference at the `O(step²)` floor, **swept over Mach number** |
| **FL-4** | Solve the rotating Solov'ev. | §6.2: `deltaStarFD()` on their closed form **first**; then `k+1` in `ψ` and `q` over four dyadic meshes, `k = 1,2,3` |
| **FL-5** | The manufactured nonlinear rotating case. | §6.4: Newton observed order 2; assembled Jacobian against a difference of the assembled residual |
| **FL-6** | `n` species: the safeguarded root find and `∂φ₀/∂ψ`. | §6.5, all four |
| **FL-7** | Normalised flux: `meq::NormalisedRotatingSource` through the bordered Newton. | `ψ_ax − max ψ_h` at machine zero and Newton order 2, exactly as `HighBetaConvergence` asserts today |
| **FL-8** | Through the driver: TOML, and `n_s`, `φ₀` written. | Driver against library on the same configuration, as `DriverAcceptance` does; a worked `examples/rotating-*.toml` |

**FL-2 is the stage to protect.** It exercises the species container, the profile
chain rule, `φ₀`, the units conversion and the sign of `Δ*` — everything
structural — against an answer meq already has to fifteen digits. If a `4π`, a
`2π` or a sign is wrong, FL-2 is where it shows, and it shows as a wrong number
rather than as a plausible equilibrium.

**FL-4 before FL-5 is deliberate.** FL-4 is linear, so a failure there is the
discretisation or the source and cannot be the Jacobian. FL-5 is the first stage
where a Jacobian error is possible, and by then everything else is pinned.

---

## 8. Risks, in the order they are likely to bite

**The source is exponential in `M²`, and meq has no experience of that.**
`p ∝ exp[ m_i ω² r² / 2T_eff ]` varies by `e^{M²/2}` across a surface — a factor
of 7 at `M = 2`. `CLAUDE.md`'s reaction-ratio diagnostic `max|∂F/∂ψ|/λ₁` is the
thing to measure, over the solution's *actual* range and not over a nominal one,
which is a mistake this project has already made once. Measure it at FL-3, before
any solve depends on it.

**Rotation is reported to make the iteration harder.** Li & Zhu note that the
optimal range of their Picard relaxation parameter *narrows* with rotation. meq
uses Newton, which should help, but the honest position is that this is unmeasured
and the reactive ladder should be exercised on the rotating cases rather than
assumed adequate.

**High-Mach rotating equilibria bifurcate.** That meets, from a third direction,
what `CLAUDE.md` already records: coarse discretisations of these sources carry
more than one solution, and three solve routes gave three answers 9.4% apart. **A
rotating failure must not be answered by quietly changing the globalisation** —
it changes which equilibrium is reported.

**The trivial branch may or may not apply.** `F(r, z, 0) ≠ 0` in general here, as
with the high-β source, so the GS-2 trivial-branch trap may not bite — but the
multiple-root problem that fixture *did* have may. `setInitialGuess()` exists and
this path will probably need it; treat that as part of the problem statement, not
as an optimisation.

**Somebody will compare profiles against another code and be off by the gauge.**
§3.3. It is the one part of the gauge freedom that is not a simplification, and
it is a documentation problem rather than a numerical one — which is precisely
the kind this project has found hardest to notice.

---

## 9. Deliberately out of scope

* **The fixed-`q(ψ)` solver.** `I(ψ)` is an input here. Taking `q(ψ)` instead and
  finding `I(ψ)` from it is a genuine extension — RoPP (142) gives
  `q = V′ I ⟨r^{-2}⟩_ψ / 4π²`, which is a **flux-surface average**, so it needs
  exactly the machinery §3.3 avoids, plus an inner iteration for `I`. It is the
  natural next thing and it is a `ROADMAP.md` item, not a section here.
* **Poloidal flow.** A different equation, of variable type. FLOW's problem, not
  this one.
* **Anisotropic pressure.** Independent of rotation and with no reference pinned;
  `TODO` carries the seeds.
* **The transport-timescale evolution.** RoPP §8 and (252)–(254) evolve `n_s`,
  `T_s` and `ω`. meq computes an instantaneous equilibrium; the evolution is
  MaNTA's, and `MANTA-COUPLING.md` is where that conversation lives.
* **`F₁s`, neoclassical and the gyrokinetic equation.** RoPP's §7.1 and §7.4 are
  the rest of the paper and are nothing to do with meq.
