# Toroidal flow: the generalised Grad-Shafranov equation of RoPP (136)

**MEQ's rotating-source half**, split out of `CLAUDE_HDGGS.md` because it is a
campaign rather than part of the core method — FL-0 to FL-8, **done**.
`CLAUDE.md` is the index; `CLAUDE_HDGGS.md` has the equation, the
discretisation, Newton and the bordered Newton that this source is driven
through, and `CLAUDE_FB.md` and `CLAUDE_INVERSION.md` are the other two
campaigns.

**`docs/rotation.rst` is the derivation** — the plan file `FLOW-PLAN.md` was
converted to it once everything it staged was built, and `RotatingSource.hpp`
defers to it. What is here is what a maintainer needs: the measurements, the
three errors found in Li & Zhu, and the traps.

**Nothing here changes the discretisation.** A rotating equilibrium is the same
equation with a different `F`, and `meq::RotatingSource` is a `meq::Source` —
which is why the `k+1` / `k+2` rates, the sign conventions and the Newton
machinery are all `CLAUDE_HDGGS.md`'s and are inherited rather than restated.

**Measurement tables live in `MEASUREMENTS.md`** under stable `M-nn` anchors.
**Do not renumber the anchors.**

## Toroidal flow: FL-0 to FL-8 are done and green

**`meq::RotatingSource` solves the generalised Grad-Shafranov equation of
`refs/RotatingGK.pdf` (136)**, closed by its (96) and (97), for two species in
the local gauge `φ₀(r_ref) = 0`. `docs/rotation.rst` is the derivation; this
section is only what a maintainer needs. **Three or more
species is `Closure::RootFind`** — a safeguarded scalar Newton on (97) with
`φ₀`'s two `ψ`-derivatives by implicit differentiation — and
`meq::NormalisedRotatingSource` puts the profiles in normalised flux, where
`ψ_ax` is an unknown and the existing bordered Newton closes it.

**FL-8 IS DONE: ROTATION IS REACHABLE FROM A TOML FILE.** `[source] Type =
"rotating"` takes an array of `[[source.species]]` tables — mass, charge and a
temperature profile each — plus `Omega`, one density profile and `GGPrime`.
That array is **MEQ's first array-of-tables in the schema**, read by
`Table::getTableArrayOr()`, which names its elements `source.species[i]` so a
fault in the third one says so. `examples/rotating-rectangle.toml` is the
worked example and `examples/rotating-normalised.toml` the same thing in
normalised flux, where `ψ_ax` is an unknown and the bordered Newton closes it —
`makeSource` **throws** on a normalised config and `makeNormalisedSource` is
the other door, because handing a normalised source to the plain path would
converge to a different equilibrium rather than fail.

**AND ROTATION REACHES `ψ` ITSELF THERE, WHICH IS A MEASUREMENT RATHER THAN A
CLAIM.** Against the same configuration with `Omega = 0` and nothing else
changed, `‖ψ(ω) − ψ(0)‖/‖ψ(0)‖ = 1.2683e-01` in L2 over 4608 dofs, and
`max ψ` moves 9.755876e-02 → 1.079860e-01, up 10.7%.
`theDriverSolvesARotatingEquilibrium` re-measures both every run and gates the
shift at `> 5e-2`, so an edit that shrinks the pressure term fails loudly
instead of quietly re-vacuating the example. **It is not free**: it needs
`GGPrime` down from 2.0 to 0.8 so the pressure term is a real part of `F`
rather than a correction on `g g′` — 15 / 29 / **47** / 66 / 81 % of `F` from
`RMin` to `RMax` at `Omega = 4.0e5`, which is `M = 1.33` and squarely the sonic
regime RoPP is about.

**TWO PROFILE TABLES SHIP, AND THE REASON IS A TRAP.** Re-expressing `f(ψ)`
against `Ψ = ψ/ψ_ax` divides the abscissa by the axis flux and **multiplies the
derivative column by it**. `examples/rotating-density.dat` and
`rotating-density-normalised.dat` are the same profile written both ways; hand
one to the wrong configuration and it parses, solves and converges to a plasma
whose density gradient is out by a factor of ten. Both headers say so, and they
are also the only place in `examples/` that documents the `SplineProfile` file
format at all.

**(136) collapses to `F = μ₀ r² ∂p/∂ψ|_r + g g′`** with `p = Σ_s n_s T_s`,
because the `∂φ₀/∂ψ` terms cancel identically against quasineutrality. So the
residual needs `φ₀` and never its derivative; only the Jacobian does. **Two
species need no root find at all** — (97) is linear in `φ₀` after logs, giving
`C = ω²(Z₁m₂ − Z₂m₁)/(Z₁T₂ − Z₂T₁)` as the exponent both species share, exact
with the electron mass kept. **Three or more species are handled by the root
find**, `Closure::Automatic` picks between the two, and the ceiling is
`meq::maxSpecies = 8`, enforced at parse rather than at construction.

**IT COST `meq::Profile` A THIRD DERIVATIVE LEVEL, AND THAT IS THE ONE
STRUCTURAL CHANGE.** `MHDSource` stores the *products* `p′` and `g g′`, so `F`
is one evaluation and `∂F/∂ψ` is one `prime()` — two levels, which is all
`Profile` had. A rotating `F` is *already* `∂p/∂ψ` of something built from flux
functions, so the Jacobian spends a second derivative of every input. No
reparametrisation avoids it. `doublePrime()` is a **pure** virtual, so every
subclass had to answer rather than one silently returning zero. **The caveat is
real**: a Hermite cubic is `C¹`, so its second derivative jumps at every
interior knot — `the_second_derivative_jumps_at_an_interior_knot` asserts that
rather than pretending otherwise, and nothing has measured what it costs Newton.

## What is measured

→ **[M-41](MEASUREMENTS.md#m-41)** — `φ₀` against a **brentq root of (97)**, independent Python

**THREE PAPER ERRORS WERE FOUND ON THE WAY AND ALL THREE ARE THE KIND THAT
CONVERGE BEAUTIFULLY.** Li & Zhu's (12)–(16) write `M₀²` where their prose
defines `M₀` as the group without a square root — the group is an energy ratio,
so it is a Mach number *squared*, and `RotatingSoloviev.hpp` names its member
`machSquared` and cites the exponent rather than their symbol. Their (9) carries
**two reversed signs**, on the `dΩ/dψ` and `dT/dψ` corrections, found
independently three times — by transcription, by an unrelated numerical check,
and by MEQ's own derivation agreeing with the corrected form. **Neither of their
own benchmarks can see it**, because both have `dC/dψ = 0`. And their (6) omits
the `μ₀` their (9) carries.

**AND THERE IS AN INDEPENDENT IMPLEMENTATION OF (136) ON THIS MACHINE, WHICH
WOULD BE THE FIRST OUTSIDE CHECK OF THAT GAP.** `../freegs4e`'s
`ProfilesCentrifugalMirror` — driven by `../geq`, a rotating-mirror wrapper —
implements Abel (136) closed by its (96) and (97), the same paper from the same
equations, in Python, by finite differences and Picard. So it and
`meq::RotatingSource` are two independent implementations of one equation, which
is a rarer thing than it sounds and is exactly what the `C′(ψ)` gap wants.

**The smallest comparison is source against source and needs no solve**: geq's
`Jtor × μ₀R` against MEQ's `F`, pointwise on a prescribed `ψ` and prescribed
species profiles. `geq/tests/test_jtor_abel.py` already does this against a
hand-written reference, so MEQ's would substitute for that reference.

**RECONCILE THE `φ₀` GAUGE FIRST — IT IS NOT A DETAIL.** MEQ pins
`φ₀(ReferenceRadius, ψ) = 0`, one condition per flux surface. geq pins **one
point globally**, at `ψ_n = 0.5` on the midplane. `N_s(ψ)` absorbs the difference
through `exp(Z_s e δ/T_s)`, so **the two codes' tabulated densities mean
different things** until it is fixed. Compare `n_s(r,z)` and `F`; never `N_s(ψ)`.
Smaller traps beside it: geq wants `T` in **eV** and MEQ in **Joules**; geq
derives `ω` from a Mach number or a voltage and re-derives it every Picard step,
so its `omega_profile` callable is the only way to hand both codes the same
`ω(ψ)`; and geq sets `gg′ ≡ 0` unconditionally, which on MEQ's side means
`GGPrime = 0` and the trivial-branch warning in
`examples/rotating-rectangle.toml`.

**Nothing on that side runs here yet** — `freegs4e` is not importable and geq's
paths point at another machine.

**THE `C′(ψ)` GAP IS CLOSED, BY A MANUFACTURED FIXTURE, AND NO PUBLISHED
BENCHMARK COULD HAVE CLOSED IT.** `C` constant is what collapses Maschke &
Perrin's (4.6) to its (4.8) and makes the equation solvable at all, so a varying
`C` is exactly what has no closed form — the route is manufactured by necessity
and there is no point looking for another paper.
`tests/analytic/VaryingCentrifugal.hpp` and
`tests/convergence/VaryingCentrifugalConvergence.cpp` are it.

**THE CONSTRUCTION RUNS BACKWARDS FROM `C`, AND THAT IS THE DESIGN DECISION.**
Choosing `ω` and the temperatures and reading `C` off leaves it a ratio whose two
derivatives are awkward to state; prescribing `C(ψ)` as a quadratic and deriving
`ω = √( C D / K )` from (97) makes `C`, `C′` and `C″` two lines each and leaves
the temperatures free. `C` then drifts **2.50×** over the sweep with `C″ = 2.00`
exactly, the densities vary **33.8×** across the box, and the two
implementations agree at **7.13e-16 in `F`** and **5.26e-16 in `dF/dψ`**:

→ **[M-68](MEASUREMENTS.md#m-68)** — both closures · `F` · `dF/dψ` · `p` · `φ₀`

**AND TWO TRANSCRIPTIONS AGREEING ONLY PROVES THE COPYING WAS CONSISTENT**, so
the fixture carries a second route to its own pressure: `potentialByBisection()`
solves (97) by bisection with **`C` appearing nowhere**, and differencing it
reaches the closed form's derivatives. A slip shared between MEQ's `C`-chain and
the fixture's would pass the pointwise comparison and fail this one.

→ **[M-69](MEASUREMENTS.md#m-69)** — (97)'s residual · `p` · `dp/dψ` · `d²p/dψ²`

→ **[M-70](MEASUREMENTS.md#m-70)** — the mutations, and `p` not moving under them

**THE STANDING CLAIM THAT ONLY THE `dFdPsi` SWEEP TOUCHED `C′` IS WRONG, AND THE
REAL GAP WAS NARROWER AND SHARPER.** `RotatingSourceTests`'
`thePressureMatchesTheIsothermalClosedForm` writes the exponent out by hand and
pins `pressure()` and `densityExponent()` against it to **1e-14**, at profiles
whose `C` drifts 1.21× — so **`C` itself was already independently pinned**, and
`RotatingNewtonConvergence` runs at `C′ ≠ 0` as well. What had no analytic
reference is `C′` and `C″`, which enter nothing but `dp/dψ` and `d²p/dψ²` — that
is, nothing but `f()` and `dFdPsi()` — and were checked only by central
differences at 1e-6 to 1e-7. An `O(1)` error would have been caught; a term
smaller than the differencing floor would not.

**AND IT FOUND A HOLE NOBODY HAD NAMED, IN THE CHARGE WEIGHTING.**
`closedFormState()` is built on `Z₁T₂ − Z₂T₁`, `Z₁m₂ − Z₂m₁` and `m₁T₂ − m₂T₁`,
and **every two-species rotating configuration in the tree runs at `Z = ±1`**,
where those collapse to `T₁ + T₂` and `m₁ + m₂`. `Z = 6` appears only in
three-species sets, which take `RootFind` and never form the closed-form
combinations at all. **A closed form that wrote the plain sums would have passed
everything in this tree.** The fixture is at `Z₁ = +2`, where `Z₁T₂ − Z₂T₁` is
40.0% away from `T₁ + T₂` — nine orders above the tolerance — so the agreement
above says the *weighting* is right and not only the chain rule.
`theChargeWeightedCombinationsAreNotPlainSums` records it.

**THE MANUFACTURED SOLVE IS REACHABLE AND ITS VALUE IS NARROW.** `meq::Source`
is `f( r, z, ψ )`, so a test-local source adds a `ψ`-independent remainder built
from **the fixture** rather than from MEQ's own `f()` — which is what
`RotatingNewtonConvergence` does, and is why there `ψ_e` stays exact whatever
`f()` computes. Rates are 1.999 / 2.997 / 3.999 in `ψ` at `k = 1, 2, 3`, Newton
three steps throughout. The pointwise sweep is nine orders sharper, so the solve
is not what closes the gap; what it reaches is the **assembly** — `SourceIntegrator`
evaluating this source and its Jacobian at quadrature points, at `Z₁ = +2` and a
quadratic `C`.

**THE +5% MUTATION REPRODUCES ON A THIRD SOURCE.** Perturbing `dFdPsi` by 5%
leaves both `L2` errors unchanged to seven figures and moves the observed Newton
order 1.792 → 1.173 and the count 3 → 4. **1.792 is not a shortfall**: the first
step already takes 500× off the residual, so the best triple straddles the
pre-asymptotic step — `r₂/r₁² = 1.2e-02` says the iteration is quadratic and
there is simply no clean tail. The threshold was calibrated between the two
readings rather than copied.

**`ψ`-DEPENDENT `T` AND `ω` ARE NO LONGER PART OF THAT GAP, AND THE TWO ARE
WORTH KEEPING APART.** Maschke & Perrin's §4 leaves `T(ψ)` free and `ω(ψ)` with
it, so `MaschkePerrinConvergence.cpp` drives `meq::RotatingSource` on five
varying profiles against a closed form. **`C′ ≠ 0` is what remains**, and no
exact solution can reach it — see *Maschke & Perrin is the second exact rotating
benchmark*.

**THE +5% MUTATION TEST WAS RE-RUN ON THE ROTATING CASE AND REPRODUCES THE
STATIC RESULT EXACTLY.** Perturbing `RotatingSource::dFdPsi` by 5% leaves every
L2 error and every convergence rate unchanged **to all seven digits printed**, at
`k = 1, 2, 3`; Newton goes 3 iterations to 6, observed order 1.980 to 1.055, and
the assembled-Jacobian check 3.5e-11 to 2.1e-04. So the rate tables FL-4 rests on
are blind to the defect, and the Jacobian check and the order are what see it.

**THREE TRAPS FROM FL-5 TO FL-7, ALL OF THE SAME SPECIES.**

**A control on the reference curve is blind to the entire rotation chain rule.**
The gauge pins `φ₀(r_ref) = 0`, so at `r = r_ref` the exponent and *both* its
`ψ`-derivatives vanish and `∂²p/∂ψ²` collapses to `P₀″(ψ)` — the answer a
non-rotating source gives. A check that `∂F/∂ψ` varies with `ψ`, placed there,
read 0.333 against 7.300 at the outboard edge. **The one radius where the gauge
is exact is the one radius where a rotating source is indistinguishable from a
static one.**

**The best observed Newton order is not the order.** The bordered history opens
6.24e-02 → 4.01e-02 → 7.44e-03, whose "order" reads **3.81** — the iterate
walking into the basin. Asserting on the best triple would pass on that and keep
passing with a Jacobian degraded enough to destroy the tail. The assertion is on
the last triple above the round-off floor and is bounded **both** sides: 1 is a
broken Jacobian, 3.8 is an artefact.

**A comparison against an exact zero has no relative tolerance.** `φ₀` vanishes
identically on `r = r_ref`, and identically everywhere at `ω = 0`, so the
root-find-against-closed-form check compared 1e-33 with 0.0 and failed at every
such point. The fix is a floor at the problem's own energy scale, not a
case-dependent one read off the configuration — which was the first attempt and
failed again at `ω = 0`, where the scale is itself zero.

## Maschke & Perrin is the second exact rotating benchmark, and what it buys is the SOURCE and not the discretisation

`tests/analytic/MaschkePerrin.hpp` and `tests/convergence/MaschkePerrinConvergence.cpp`.

**THE FLOW PLAN REJECTED THIS PAPER AS AN ADIABATIC CLOSURE ON THE STRENGTH OF A
`γ` IN THE EQUATIONS, AND THAT READING IS WRONG.** Wrong section:
`refs/MaschkePerrin.pdf` — *Plasma Physics* **22** (1980) 579, not the Phys.
Lett. A 102 (1984) everyone cites — carries two solutions, and its **§4** takes
the temperature as a surface quantity and is (136)'s isothermal closure. `γ`
appears there once, only inside `γΩ²`, whose job is to convert between the
adiabatic sound speed of (4.11)'s Mach number and the isothermal one the
equation actually knows, so it cancels out of the solution: **every `γ` works,
not just `γ = 1`**, and `γΩ²` is the isothermal Mach number squared at `R₀`. Its
**§3**, the actual polytrope, carries a **power law** where §4 carries an
exponential and is **not** ours. The mistake is recorded because it is this
project's standing hazard — three closures that look alike — biting the plan
that warns about it. `μ₀` is restored by `p_SI = p_M&P/μ₀`, the paper working in
`j = ∇×B`.

**AND THE PDE IS LI & ZHU'S RENAMED, WHICH IS THE FINDING RATHER THAN A
DISAPPOINTMENT.** Read (4.10) against Li & Zhu's (12):

```
Δ*ψ = −p₁ r² exp[ M₀²( r²/R₀² − 1 ) ] − F₀           Li & Zhu (12)
Δ*ψ = −( P/R₀⁴ ) r² exp[ m r²/2R₀² ]   − M/R₀²       M&P (4.10),  m := γΩ²
```

— the same equation, differing by a constant absorbed into the amplitude; and
(4.17)'s harmonic terms `1`, `R²`, `X²R² − R⁴/4` are three of Li & Zhu's four,
with the same exponential particular solution. **So a convergence study on it is
a restatement of `RotatingSolovievConvergence.cpp`**, and anything claiming this
fixture is a second independent test of the *discretisation* is wrong.

**WHAT IT BUYS IS THE ONLY EXACT SOLUTION MEQ HAS WITH `T′ ≠ 0` AND `ω′ ≠ 0`.**
(4.7) constrains only the **ratio** `ω²/(R̄T)`, leaving `T(F)` an arbitrary
surface function and `ω(F)` following it through (4.13). Every other closed-form
check of `meq::RotatingSource` in this tree runs at **constant** `T` and
**constant** `ω` — `RotatingSourceConvergence.cpp`'s species carry
`ConstantMassProfile` for both, because Li & Zhu's Solov'ev case holds `T₀` and
`Ω₀` constant and no closed form survives otherwise. What checks the varying
case is a central difference of MEQ's own `f()`, which **cannot see a term
missing from both `f()` and `dFdPsi()`**.

In MEQ's language §4 is exactly two conditions: **`C(ψ)` constant** (that is
(4.7)) and **`p₀(ψ)` linear** (that is (4.9)). Satisfying both while leaving
every profile varying takes one shape function `θ(ψ) = 1 + σψ`, with
`T_s = τ_sθ`, `ω = ω₀√θ` — which **is** (4.13) — and `n_s0 = c_s N(ψ)/θ` with
`N` linear and `Σ Z_s c_s = 0`. **Five profiles varying, and their variations
cancel to a source independent of `ψ`.** Measured over the benchmark box:

| | |
|---|---|
| `meq::RotatingSource` against (4.10) | **4.3e-16** relative |
| `dF/dψ`, which (4.9)'s linear `p_T` makes zero | **2.6e-15** |
| the shared exponent's drift with `ψ`, i.e. (4.7) | **2.2e-16** |
| the terms cancelling inside `μ₀r² d²p/dψ²` | **1.957**, so the zero above is not vacuous |
| the same with (4.7) **broken** | **2.9e-01** and `dF/dψ = 5.8e+00` |

and the production source drives the solve to **1.996 / 2.997 / 3.998** in `ψ`
and 1.987 / 2.989 / 3.989 in `q`, Newton taking **1** step everywhere.

**IT DOES NOT CLOSE THE `C′(ψ)` GAP AND NO EXACT SOLUTION CAN.** `C` constant is
precisely what collapses (4.6) to (4.8); a varying `C` is what makes the
equation unsolvable in closed form. `tests/analytic/VaryingCentrifugal.hpp` is
what closes it, and it is manufactured for exactly this reason — see *Toroidal
flow* above.

**THE GEOMETRIC CONSTANTS NEED THEIR OWN CHECK, AND SUBSTITUTION INTO THE PDE
CANNOT SUPPLY IT.** `C` of (4.18) and `ε_a` multiply `Δ*`-**harmonic** terms, so
a wrong one leaves `F`, `Δ*ψ` and every rate exact — `Soloviev.hpp`'s lesson on a
different paper. (4.19)'s axis ellipticity is the independent statement:
separately published, recoverable from `ψ`'s own Hessian, reading **2.297442550
against 2.297442541**. Mutating `C` is caught by that case and by nothing else in
the file — the `Δ*` scan and all four rate tables stay green, which is the
property being relied on rather than an accident. **(4.19) holds only at `M = 0`**,
the poloidal-current term moving `∂²ψ/∂X²` and not `∂²ψ/∂R²`: at `M = 0.3` the
true elongation is 1.767 against (4.19)'s 2.297, and that 23% is asserted as a
control so the scope is not a footnote.
