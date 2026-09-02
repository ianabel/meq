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

The general two-species form, which is what is implemented and which
specialises to the ion/electron one below, is

```
e φ₀ = ω² ( r² − r_ref² )( m₁T₂ − m₂T₁ ) / 2( Z₁T₂ − Z₂T₁ )
C    = ω² ( Z₁m₂ − Z₂m₁ ) / ( Z₁T₂ − Z₂T₁ )
```

— the numerator of `φ₀` being a mass-weighted temperature difference, so that
two species with equal `m/T` leave nothing for the field to separate and `φ₀`
vanishes identically. In the ion/electron case:

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
    double charge;                                // Z_s, signed and dimensionless
    std::shared_ptr<Profile const> temperature;   // T_s(psi), JOULES
    std::shared_ptr<Profile const> density;       // n_s0(psi) on r = r_ref, m^-3
};
```

`src/meq` deliberately keeps MFEM out of `Profiles` and `Source`; this stays on
the right side of that line and is unit-testable without the library.

**Built 2026-09-01** in `src/meq/RotatingSource.{hpp,cpp}`, with
`chargeNeutralityResidual()` and `neutralisingDensity()` beside it. The charge is
`Z_s` rather than `Z_s e`, so the elementary charge appears in exactly one place;
`potential()` returns `e φ₀` in Joules for the same reason, that being the
combination every exponent actually contains.

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

### 5.3 `dFdPsi` is where the risk is, and it cost `Profile` a new method

`∂F/∂ψ = μ₀ r² ∂²p/∂ψ² + (g g′)′`, and `∂²p/∂ψ²` runs the chain rule twice
through everything in §5.1 plus `φ₀`. **The finite-difference check on `dFdPsi`
stops being box-ticking and becomes the load-bearing test of this plan.** It
must be swept over Mach number, not evaluated at one point: the terms it is
checking are the ones that vanish at `ω = 0`.

**IMPLEMENTED 2026-09-01, AND IT NEEDED AN INTERFACE CHANGE THIS PLAN DID NOT
FORESEE.** `meq::Profile` supplied `operator()` and `prime()` and nothing more —
exactly two levels, which is what `meq::MHDSource` needs, because it stores the
*products* `p′` and `g g′`, so `F` is one evaluation and `∂F/∂ψ` is one
`prime()`. A rotating source needs **three**: `p = P₀(ψ) exp(C(ψ)Δ/2)` is built
from flux functions, `F` is already `∂p/∂ψ` and spends one derivative of each,
and the Jacobian spends a second. No reparametrisation avoids it — `p` is not a
flux function, so there is no product to pre-store, and taking `P₀′` as the input
instead loses `P₀`.

So `Profile` grew **`doublePrime()`**, as a *pure* virtual rather than one with a
default: there are three subclasses in the tree, and a compile error naming each
is better than a runtime surprise in one. `ConstantProfile` returns zero and
`HermiteCubicSpline` the exact second derivative of its own cubic. **That last
one carries a caveat worth keeping**: a Hermite cubic is `C¹` and no more, so its
second derivative is piecewise linear and **jumps at every interior knot**.
Measure zero, so it cannot move a converged answer, but a Newton step landing on
a knot sees a one-sided Jacobian. Nothing has measured what that costs, and it
does not arise for the analytic and constant profiles FL-0 to FL-5 use.

### 5.4 Normalised flux, and the bordered Newton

If `T_s`, `n_s0` and `ω` are given against normalised flux — which is how
`refs/GourdainContour.pdf` §V poses profiles and how every equilibrium code does
— then this must be a `meq::NormalisedSource` and `ψ_ax` is an unknown of the
system. **That machinery already exists and is green** (`HighBetaConvergence`),
and under NPC two of the three border quantities are exact rather than
differenced. So `meq::NormalisedRotatingSource` is a subclass away, not a
project.

`ψ_bnd` stays zero: this is still the fixed-boundary problem.

### 5.5 `Config` and `SourceFactory` — built 2026-09-01

`SourceType::Rotating`, `RotatingParameters`, `SpeciesParameters`, a
`[[source.species]]` array of tables, and the charge-neutrality solve of §3.4.
The normalised plumbing was the prerequisite and is done in the same pass, so
`meq::NormalisedMHDSource` is reachable from a TOML file too — that was
`ROADMAP.md` item 1.

**`[[source.species]]` IS THE FIRST ARRAY OF TABLES IN MEQ'S SCHEMA**, and the
existing nested-table reader would have refused one: an array of tables *is* an
array, so its `is_table()` check rejects it, even though TOML nests `[[a.b]]`
under `[a]` exactly as `[a.b]` does. `Table::getTableArrayOr()` is the new
primitive, and it names each element `source.species[i]` so that
`rejectUnknownKeys()` and `fail()` qualify a diagnostic all the way to
`source.species[0].Denisty` for free — measured, that is what it prints.

**Every profile is a constant OR a table, in two keys, with a scale.**
`Temperature` / `TemperatureFile` / `TemperatureScale`, and the same triple for
`Density`, `Omega`, `GGPrime`, `PPrime`. Two keys rather than one that dispatches
on node type, because a type-dispatched key inherits exactly the trap
`Config.cpp`'s header records: TOML distinguishes `1` from `1.0`, and toml11's
`find_or<double>` returns the *default* for an integer node rather than
failing. **The scale is per profile, not global** — "that temperature ×2,
leave the density alone" is the case it exists for — and it multiplies all three
derivative levels through `meq::ScaledProfile`, so a scaled profile is exactly
the profile of the scaled quantity and `prime()` stays exact.

**Two factories, and the second returns a NON-CONST pointer.** The solver calls
`setNormalisation()` on a normalised source before every residual evaluation, so
it cannot come back from `makeSource()` as `Source const`. `makeNormalisedSource`
is the second entry point, and **`makeSource` throws on a normalised
configuration rather than returning one** — which is the point: a
`NormalisedRotatingSource` *is-a* `Source`, so returning it would compile and
solve, with `ψ_ax` frozen at the guess for ever because nothing would move it.
The answer would be a converged solution to a problem the file did not describe.

**`Neutralising = true` marks the one species whose density is derived.** That
puts §3.4's `n − 1` counting in the file and makes the parser enforce it, rather
than asking for `n` profiles that happen to sum to zero. Exactly one species
must set it; zero or two is refused, with the counting in the message.

**The NetCDF keys are reserved and refused by name.** `[source] ProfileFile` and
the per-profile `<Key>Variable` / `<Key>Fit` are recognised and rejected with
"not implemented yet", not as unknown keys — so adding a multi-profile NetCDF
reader later is purely additive. `TODO`'s NetCDF entry carries the decisions
already taken, including that a variable with values but no derivative will be
refused unless the configuration opts in per profile, and that MaNTA's `(u, q)`
pair *is* `meq::Knot`, which makes the rotating source the easy case to build it
against rather than the awkward one.

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
temperature flux function. **That reduction is now verified rather than
assumed**: `RotatingSourceConvergence`'s `theSourceAndTheFixtureAgreeUnderRotation`
builds a two-species `meq::RotatingSource` with `T_i = T_e = ½` and matches the
fixture pointwise at `M² = 0, 1, 4`, which is the statement that our
`C = ω²(Z₁m₂ − Z₂m₁)/(Z₁T₂ − Z₂T₁)` reproduces their single `T = T_i + T_e`.

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

### 6.3 Maschke–Perrin IS a test of (136), and this section said otherwise

**CHECKED DIRECTLY 2026-09-01 AND THE EARLIER READING HERE WAS WRONG.** This
section used to reject Maschke–Perrin on the grounds that its `γ` made it an
adiabatic closure. That is the right instinct applied to the wrong section of
the paper, and it is worth leaving recorded because it is this section's own
hazard biting the section that warns about it.

`refs/MaschkePerrin.pdf` is **Maschke & Perrin, *Plasma Physics* 22 (1980)
579** — not the Phys. Lett. A 102 (1984) 106 that `TODO` and Li & Zhu's [48]
cite — and it carries **two** solutions in **two** closures:

| | closure | `γ` | ours? |
|---|---|---|---|
| **§3** | entropy a surface quantity, polytrope `p = A(S)ρ^γ` | `η = γ/(γ−1)`, **singular at `γ = 1`** | **no** — a power law in `R²` |
| **§4** | **temperature** a surface quantity, `B·∇T = 0` | appears only inside `γΩ²` | **yes** |

**Li & Zhu's §3.2 is §4**, the isothermal one — their (18)–(21) reproduce M&P
(4.9), (4.12), (4.17), (4.18) exactly. And in §4 **`γ` is vestigial**: it enters
once, at (4.7), and only as the product `γΩ²`, whose whole job is to make `Ω` the
*adiabatic* Mach number. Substituting `C ≡ γΩ²/R₀²` removes `γ` and `Ω` from the
solution entirely. So it is not that `γ = 1` is a usable special case — **every
`γ` is usable**, because `γ` is a label. (§3 → §4 as `γ → 1` as one would hope,
confirmed symbolically, but §3's formula is singular there, so the limit is not
where a usable expression lives.)

**Verified by substitution rather than by re-derivation**, which is what the
question actually was: 50-digit finite differences give
`max|Δ*ψ + F| / max|F| = 8e-26` over a 25×25 grid with `g g′ ≠ 0` exercised, and
sympy gives **exactly zero**. The only constant needing care is `μ₀`: M&P work
in `j = ∇×B` units, so `p_SI = p_M&P/μ₀`. No factor of two, no sign flip.
Their single `T` maps onto our `T_i + T_e`, and our two-species
`C = ω²(Z₁m₂ − Z₂m₁)/(Z₁T₂ − Z₂T₁)` reproduces their `ω²ρ/p` with difference
**0.000e+00** — exact, with no `m_e → 0` limit needed.

**TWO CAVEATS BEFORE IT GOES IN `tests/analytic/`, AND THE FIRST IS THE
IMPORTANT ONE.** M&P's (4.7) *forces* `C` to be a constant, so the
`P₀(ψ) C′(ψ)(r² − r_ref²)/2` term of `∂p/∂ψ` is identically zero in this
fixture. **Neither published rotating benchmark can see that term** — Li & Zhu's
Solov'ev case has `T` and `Ω` constant for the same reason — and it is precisely
the term whose sign they got wrong. Only the `dFdPsi` sweep of §5.3, over
profiles with genuine `ψ`-dependence in `ω` and `T`, touches it. And second: `p_T`
is linear in `ψ` and `g g′` is constant, so `∂F/∂ψ = 0` and this sits **beside
`Soloviev.hpp` on meq's ladder, not beside `McCarthy.hpp`**. Within M&P's ansatz
it cannot be made into a Jacobian test.

* **M&P §3, the actual polytrope**, is a different closure and is **not** usable:
  `p = (P/R₀⁴)(ψ−F₁)[1 + Ω²R²/2R₀²]^η` is a power law in `R²`, not an
  exponential, and no choice of our `C` represents it. Measured best-fit
  mismatch over the benchmark box: 4.1e-3 at `γ = 5/3`, falling to 3.9e-5 at
  `γ = 1.01` and to zero only in the singular limit.
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
| **FL-0** ✔ | `meq::Species`, the profile container, the charge-neutrality solve. No solver, no MFEM. | **DONE.** `neutralisingDensity()` closes `Σ_s Z_s n_s0 = 0` to 1e-14 and is exact at all three derivative levels; every refusal — three species, same-sign charges, a violated neutrality set, a null profile, a bad radius — throws a named `invalid_argument` |
| **FL-1** ✔ | `φ₀` and `p(r,ψ)` for two species, closed form. Still no solver. | **DONE.** `φ₀(r_ref) = 0` **exactly**; `Σ_s Z_s n_s = 0` to 1e-13 over a Mach sweep; `p` against Li & Zhu (8) to 1e-14. And the closed form checked against a **brentq root of (97)** in an independent Python implementation: **1.9e-14** relative, with the two species' exponents agreeing to 3.0e-14 |
| **FL-2** ✔ | `meq::RotatingSource`, and the `ω → 0` collapse. | **DONE, and better than asked.** Pointwise against `MHDSource` at 1e-13 in both `f` and `dFdPsi`, by both routes to no rotation. Through the solver, the Solov'ev study driven from a `RotatingSource` reproduces `SolovievConvergence`'s errors **to every printed digit** — 3.596959e-05, 1.916351e-07, 5.510484e-10 in `ψ` — not merely the rates |
| **FL-3** ✔ | `dFdPsi`, two species. | **DONE.** Central difference at the `O(h²)` floor over five Mach scales and two steps. **This is the only test in the stage that touches the `C′(ψ)` term** — see §6.3; neither published benchmark can |
| **FL-4** ✔ | Solve the rotating Solov'ev. | **DONE.** `deltaStarFD` first: 1.6e-08 at `M² = 0`, 2.8e-07 at `M² = 1`. Then 1.995 / 2.995 / 3.995 in `ψ` and 1.980 / 2.985 / 3.984 in `q`, Newton = 1 throughout as an affine system must give. And the source and the fixture — two independent implementations — agree pointwise at `M² = 0, 1, 4` to 1e-12, with the solver driven from the source giving the same errors as driven from the fixture |
| **FL-5** ✔ | The manufactured nonlinear rotating case. | **DONE.** `F = rot.f + h(r,z)` with `h` carrying no `ψ`, so the whole Jacobian is the rotating source's. Assembled Jacobian against a difference of the assembled residual at **3.1e-11**; Newton order **1.980**; `k+1` at 2.007 / 2.999 / 4.002 in `ψ`. And **mutation-tested**: `1.05×dFdPsi` leaves every error and every rate unchanged to all seven digits, and moves the order to 1.055 and the Jacobian check to 2.1e-04 |
| **FL-6** ✔ | `n` species: the safeguarded root find and `∂φ₀/∂ψ`. | **DONE.** At two species the root find reproduces the closed form — `φ₀` to 1e-12, `∂φ₀/∂ψ` to 1e-11, `F` to 1e-11, `∂F/∂ψ` to 1e-9. Three species (D, C⁶⁺, e) hold `Σ_s Z_s n_s = 0` to 1e-12 at every radius, `∂φ₀/∂ψ` matches a central difference, and the carbon is centrifugally enriched outboard |
| **FL-7** ✔ | Normalised flux: `meq::NormalisedRotatingSource` through the bordered Newton. | **DONE.** `ψ_ax − max ψ_h` at **0.000e+00** on three meshes, `ψ_ax` converging 1.146743e-01 → 1.146034e-01, rotation moving it by 6.0%, and Newton's **tail** order 2.000 |
| **FL-8** ✔ | Through the driver: TOML, and `n_s`, `φ₀` written. | **DONE.** Driver against library at **1.310e-16** rotating and **1.013e-16** normalised, with `ψ_ax` = 1.039163167e-01 and its constraint at −1.388e-17 in 7 bordered steps. And rotation **reaches `ψ`**, which is the point: against the same configuration at `Omega = 0`, `‖ψ(ω) − ψ(0)‖/‖ψ(0)‖ = 1.2683e-01` over 4608 dofs and `max ψ` moves 9.755876e-02 → 1.079860e-01, gated at `> 5e-2` so it cannot go stale. `examples/rotating-rectangle.toml` and `rotating-normalised.toml` are the worked pair. **It also measured item 11**: `n_s` is algebraic in `(r, ψ)`, so the same closed form evaluated at a band node's own radius reads 5.132e-07 against **8.571e-02** at the foot on `Γ_h` — a factor of 1.7e5, which is what `B` silently takes today |

**FL-2 is the stage to protect.** It exercises the species container, the profile
chain rule, `φ₀`, the units conversion and the sign of `Δ*` — everything
structural — against an answer meq already has to fifteen digits. If a `4π`, a
`2π` or a sign is wrong, FL-2 is where it shows, and it shows as a wrong number
rather than as a plausible equilibrium.

**FL-4 before FL-5 is deliberate.** FL-4 is linear, so a failure there is the
discretisation or the source and cannot be the Jacobian. FL-5 is the first stage
where a Jacobian error is possible, and by then everything else is pinned.

**AND FL-5 EARNED ITS PLACE, MEASURED RATHER THAN ARGUED.** Perturbing
`dFdPsi` by +5% and re-running the whole study leaves **every L2 error and every
convergence rate unchanged to all seven digits printed**, at every `k` — exactly
reproducing what `CLAUDE.md` records for the static case. What moves is Newton:
3 iterations to 6, observed order 1.980 to 1.055, and the assembled-Jacobian
check from 3.5e-11 to 2.1e-04. So FL-4's rate tables, which are the instrument
the earlier stages rest on, are **blind** to the defect FL-5 exists to catch.

**AND THREE MORE FROM FL-5 TO FL-7.**

**A control on the reference curve would have been blind to the whole rotation
chain rule.** FL-5's check that `∂F/∂ψ` genuinely varies with `ψ` was first
placed at `r = r_ref` and read 0.333 against a wanted 0.5. Not a defect in the
profiles: the gauge pins `φ₀(r_ref) = 0`, so the exponent `(r² − r_ref²)/2` and
*both* its `ψ`-derivatives vanish there and `∂²p/∂ψ²` collapses to `P₀″(ψ)` —
the answer a non-rotating source gives. Moved to `r = r_max` it reads 7.300.
**The one radius at which the gauge is exact is the one radius at which a
rotating source is indistinguishable from a static one.**

**The best observed Newton order is not the order.** FL-7's bordered history
opens 6.24e-02 → 4.01e-02 → 7.44e-03, whose "order" is **3.81** — the iterate
walking into the basin, not converging in it. Taking the best triple would
report that and pass, and would go on passing with a Jacobian degraded enough to
destroy the tail. The assertion is on the last triple above the round-off floor,
**bounded on both sides**, because 1 is a broken Jacobian and 3.8 is an artefact.

**`meq::Profile`'s clamping meets a real range.** A manufactured `ψ` runs over
roughly `[−0.11, +1.00]` on the standard box, not `[0, 1]`, and a Newton iterate
overshoots further. Every temperature profile therefore has to stay strictly
positive well outside the nominal range or the closure's denominator goes
singular mid-quadrature — FL-5 asserts the range and the positivity rather than
hoping.

**THREE THINGS FL-0 TO FL-4 TURNED UP THAT THIS PLAN DID NOT PREDICT.**

`meq::Profile` needed a third derivative level — §5.3. Foreseeable in hindsight
and not foreseen here.

**The `M² → 0` cancellation is ALGEBRAIC, not asymptotic**, which this plan got
wrong when it called for a series expansion in small `M²`. Writing
`v = r²/R₀² − 1` and `u = M²v`, the particular solution is
`−(p₁R₀⁴/4) v² G₂(u)` with `G₂(u) = (e^u − u − 1)/u²`: `M²` cancels exactly,
nothing is ever divided by it, and `M² = 0` needs no branch at all. What does
need a series is **small `|u|`, which happens near `r = R₀` at every Mach
number** — a different and more common condition than the one predicted. The
fixture carries a control that measures the naive form failing where the
cancellation bites: 380% wrong at `M² = 1e-8`, exactly zero at `1e-10`.

**A benchmark can be too good for its own guard.** `deltaStarFD`'s central
difference cannot resolve `Δ*ψ` at `M² = 4`: the fourth derivative carries
`(2M²r/R₀²)⁴e^u`, so the truncation floor at `h = 1e-4` — already the optimum
against round-off — is 4.8e-05, above the 1e-5 gate. The `fastRotating()`
factory is kept and used for solves, and is excluded from the scan, because
loosening the one guard standing between a mistyped term and a beautiful wrong
answer was not worth an extra Mach number.

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
