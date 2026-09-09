# Solution inversion: `ψ(R, z)` to `R(Ψ, l)`, `z(Ψ, l)`

**`CLAUDE.md`'s flux-surface half.** `CLAUDE.md` is still the maintainer's index;
`CLAUDE_HDGGS.md` has the discretisation and the solve this work post-processes,
and `CLAUDE_FB.md` has free boundary, which IN-5 waited on. `CLAUDE_FLOW.md` is
the third campaign.

**The plan is `INVERSION-PLAN.md`** and the stage names (IN-A, IN-0 … IN-P) and
section numbers cited below are its. `MANTA-COUPLING.md` is the consumer this
exists for, written from MaNTA's side. `ROADMAP.md` item 10 is where it sits.

**Measurement tables live in `MEASUREMENTS.md`** under stable `M-nn` anchors.
**Do not renumber the anchors.**

**Nothing here changes the solve.** This is post-processing, in the same way
toroidal flow — `CLAUDE_FLOW.md` — was a change to `F` alone: a new consumer of
`ψ_h` and `q_h`.

## Solution inversion: every stage is done

**`INVERSION-PLAN.md` is the design and the staged plan** — `ψ(R, z)` to
`R(Ψ, l)`, `z(Ψ, l)`, which is what `MANTA-COUPLING.md` needs, what the driver's
`(Ψ, θ)` grid is, and what `ROADMAP.md` item 10 is. **Every stage is done**,
IN-5 included — it was deferred with free boundary, and free boundary now
solves. This section is only what a reader of the code needs.

**Nothing about the solve changes.** This is post-processing, in the same way
toroidal flow was a change to `F` alone: a new consumer of `ψ_h` and `q_h`.

**AND `q` IS THE ASSET AGAIN, FOR THE THIRD TIME.** The band continuation of `ψ`
uses it, the band continuation of `B` uses it, and now the inversion does: a
critical point is a root of `q_h = 0`, which is a **solved** field converging at
the potential's own order, not a derivative of one converging an order down.
Compare CEDRES++, which records as an open problem that in P1 continuous
Galerkin the axis and the X-point are confined to mesh vertices, and TokaMaker,
which notes for Lagrange order ≥ 2 that saddles "can exist anywhere within the
mesh". MEQ resolves both sub-element by root finding.

The two pieces that exist are `src/meq/CriticalPoints.{hpp,cpp}` (stage IN-A)
and `src/meq/Zernike.{hpp,cpp}` (the basis IN-3 will fit in). Both headers carry
the long-form reasoning; what is below is what a maintainer needs to not
misread them.

### IN-A: the axis as a root of `q`, and a degree is never a count

**Measured against the analytic Solov'ev axis**, `CriticalPointConvergence.cpp`,
`k = 1, 2, 3` over `n = 4, 8, 16, 32`:

→ **[M-42](MEASUREMENTS.md#m-42)** — `k = 1` · `k = 2` · `k = 3`

against design orders of 2, 3 and 4 — clearing `k+1` itself rather than `k+1`
less slack. **The per-pair rate is not a rate here** and the test says so: the
error is *pointwise*, so it oscillates as the axis moves within its element —
4.17 / 1.54 / 1.31 at `k = 1` — and the two-tier assertion pattern
`ExtensionConvergence` needs for the same reason is what this uses.

**The sharpest assertion in the stage was not in the brief**: the ratio of the
position error to `|q_h − q|` evaluated at the *exact* axis is **0.77 to 3.37**,
against **0.77 to 3.27** predicted by linearising `q` about its root. That is a
direct statement that the root finder adds nothing to the error of the field it
is rooting — the same shape as `CLAUDE_HDGGS.md`'s *A wrong Jacobian is
invisible to a convergence table*, from the other side.

**A DEGREE IS A SUM OF INDICES AND NEVER A COUNT, and the suite demonstrates it
rather than asserting it.** `audit()` walks the mesh boundary and accumulates
the turning of `q`, which by the Poincaré index theorem is the sum of the
indices of the interior zeros — `+1` for either extremum, `−1` for a saddle. A
box drawn round `iterExample2`'s axis *and* its X-point therefore reads
**winding 0 with two critical points inside**, the saddle located to **4.5e-6**
of the published X-point. Anything that reads a zero degree as "there is nothing
here" is wrong, and `theWindingNumberIsASumOfIndicesAndNotACount` is the live
demonstration.

**The Poincaré–Hopf hypothesis is transversality of `q·n`, NOT that the boundary
is a level set**, and the two boxes in the suite are the two cases. On the
standard benchmark rectangle — which is a level set of nothing — `q·n` keeps one
sign the whole way round at `min |q·n|/|q| = 0.15`, so `winding == χ == 1` **is**
a theorem there. On a box reaching past the X-point it reads **0.00** with a
sign change, and the degree is 0 against `χ = 1` with no contradiction whatever.
`IndexAudit::transverse` records which situation the caller is in, so the
comparison is not read as a theorem where it is a coincidence.

**AND A DISCONTINUOUS `q_h` CAN CARRY BOUNDARY DEGREE 1 WITH NO ZERO IN ANY
ELEMENT.** At `h = 0.4, k = 1` the audit reads 1 and an element-by-element
search finds nothing: each element's polynomial puts its zero just inside a
neighbour's territory, and with a face jump of `O(h^{k+1})` against an element of
size `h` there is a window where the zero belongs to neither. **Poincaré–Hopf is
a theorem about continuous fields**, and this is the DG jump meeting it head on.
The window closes with refinement, so it is a property of `h` rather than a
defect — and the test pins it to that one coarse mesh, so a *finer* mesh losing
the axis fails rather than passing quietly. The practical rule: where the audit
and the search disagree, **believe the audit**.

**`CriticalPoint::overshoot` is the same phenomenon at the level of one root**,
and refusing it outright does not work: at `k = 1, n = 4`, where the axis sits
beside the mesh line `z = 0`, the two candidates are 6.6e-4 and 8.9e-2 outside
their own elements and **with no allowance at all the axis is not found**.
`setContainment()` is in *reference*-element units on purpose, so the allowance
shrinks with the mesh exactly as the ambiguity it covers does. **It is not a
tuning parameter and that was checked**: swept over 0.001, 0.01, 0.05, 0.10 and
0.20 across the whole `k × n` benchmark the located axis is identical to every
digit printed.

### `CriticalPointFinder`'s axis and `GradShafranovSolver::psiAxis()` are THE SAME POINT

`AxisConstraint::LocatedAxis` makes them one quantity: `ψ_ax` is constrained at
a zero of `q_h`, which is what this class finds. A solve and a search may still
land on the two sides of a face and differ by the `O(h^{k+1})` jump in `q_h`
across it, so expect agreement at the field's own order rather than the last bit.

**Under `AxisConstraint::NodalMaximum` — the control — they are different and
must not be reconciled.** `ψ_ax` is then *the largest nodal value*, which makes
the border row exactly `−e_j`; IN-A's axis is *the point where `q_h` vanishes*.
They differ by `O(h)` in position and `O(h²)` in value, **both independent of
`k`**, so on a refined high-order mesh the two readings *separate* rather than
converge: measured on the finest Solov'ev mesh the gap is **202×** `ψ_h`'s own L2
error at `k = 2` and **4204×** at `k = 3`. That is the reading every bordered
`ψ_ax` in this tree was published under before that change.

**AND THE NODAL DEFINITION WAS THE WEAK POINT OF THE FREE-BOUNDARY PATH**, which
is where the consequences are: nothing in *the largest nodal value of `ψ_h`* says
that value is a magnetic axis, and three configurations were found reporting one
that was not. `CLAUDE_FB.md`'s *`ψ_ax` was the weak point of the whole
free-boundary path* is that story, and the guard below is what came out of it.

**So the check exists.** `meq::CriticalPointFinder::checkAxis()` sweeps for
zeros of `q_h`, picks the O-point of the sense **the sign of the span dictates**
— which is `NormalisedSource::insidePlasma()`'s own test, and is how it avoids
`AxisSense::Either`'s refusal without guessing — and reports the **normalised
flux `Ψ` at the located axis**, which must be 1 by definition. The driver runs it
after every converged normalised solve, prints one line, and **warns** on
disagreement; `axis_normalised_flux`, `axis_r` and `axis_z` go into the `.nc`
beside `psi_axis`, because that is the interchange format and a consumer reading
`psi_axis` has nothing else to judge it by.

**Three design decisions in it are worth not undoing.** It uses `sweep()` and not
`tryFindAxis()`, because the latter seeds from the extreme nodal value — the
quantity under suspicion — and refuses outright where more than one extremum is
reachable, which the axis ridge guarantees. The test is **one-sided**: a
polynomial's peak over a closed element is ≥ its largest nodal value, so a
healthy field reads 1 **from above** and only the low side can indicate a defect.
And **the largest `Ψ` among the extrema wins**, so a spurious O-point costs a
missed detection and never a false alarm — which is the safe direction for a
warning and is the honest answer to the axis ridge, where the flux mass
`(r q, v)` degenerates and `sweep()` returns dozens of degenerate criticals.

→ **[M-46](MEASUREMENTS.md#m-46)** — `Ψ` at the located axis

Through the driver, `rotating-normalised.toml` reads **1.0003** at 768 elements
and **1.0000** at 12,288 — the `O(h²)`-in-value gap of the paragraph above,
measured. It costs 0.041 s and 0.214 s against solves of 1.0 s and 36 s.

**`findAxis()` seeds from BOTH nodal extremes, and the reason is a sign error
this file's plan carried.** MEQ's `ψ` is not sign-normalised across sources:
the Solov'ev fixtures have `F` single-signed **negative**, so `ψ` is a
*subsolution*, its maximum is on `Γ`, and **the magnetic axis is an interior
MINIMUM** (Hessian determinant +0.693, trace +2.121 on `nstx`). With `F` positive
— the high-beta source — it is a maximum. So "seed from the largest nodal value"
is right for one sign of `F` and finds **a corner of the benchmark rectangle**
for the other. `AxisSense` is there for a caller who knows which they want, and
`findAxis()` refuses rather than guesses when both are present.

**Two rings of face neighbours, not one, and that is a measurement.** Where the
critical point sits on a mesh line the extreme nodal value is the shared vertex,
and which element is credited with it is decided by an L2 jump of 1e-8: on
`iterExample2` at `k = 2, n = 24` the minimum nodal value is in element 694 and
the root is in element **696**, which is not a face neighbour of it. The seeded
path is a fast path only — where it does not produce exactly one interior
extremum the search spends a full `sweep()`, so **the answer is not allowed to
depend on the seed**.

**Handing it the raw flux block instead of `flux()` is the failure to watch
for.** The raw block holds `−q`; in even dimension `index(−v) = index(v)`, so
every winding number is unchanged and **every Maximum silently becomes a
Minimum**. The audit still passes. `sweep()` is seeded Newton and is **not
exhaustive** — the certified subdivision of `INVERSION-PLAN.md` §5 is
deliberately not built, because IN-A's acceptance needs the axis and the audit
and neither needs exhaustiveness.

### Driving by `q( ψ )`: the inversion, and the outer Newton that closes it

**`ROADMAP.md` item 10.** Every source in this tree takes the toroidal field as
input and reports `q` as an output; a transport code hands an equilibrium code
the other way round. `src/meq/SafetyFactor.{hpp,cpp}` is the inversion —
**MFEM-free, so CI gates it** — and `src/meq/SafetyFactorSolve.{hpp,cpp}` is the
loop, on KINSOL.

**THE ALGEBRA IS A DIVISION AND THE DIFFICULTY IS ALL GEOMETRY.**
`q = V′ g ⟨R⁻²⟩/4π²`, so `g = 4π² q/( V′ ⟨R⁻²⟩ )` — one division per surface at
fixed geometry. What makes it a solver is that `V′` and `⟨R⁻²⟩` are functionals
of the solution. **The round trip closes**: from a `g` 40% too large everywhere,
**6 KINSOL iterations and 20 inner solves** recover the closed form to
**2.5e-06** and `ψ_ax` to **2.35e-06**.

**THE TWO NORMALISED FLUXES RUN IN OPPOSITE DIRECTIONS**, and this is the trap
the whole file is arranged around. `FluxSurfaceFamily`'s `Ψ_N` is **zero** on the
axis; `NormalisedSource`'s `Ψ` is **one** there. So `Ψ = 1 − Ψ_N` and
`d/dΨ = −d/dΨ_N`, and a `gg′` differentiated against the family's label and
handed to the source without the sign **does not fail** — it converges, at full
order, to an equilibrium with its **shear reversed**, which is a configuration a
real machine can have. Pinned on exactly linear data, where it is an equality.

**A DAMPED PICARD CANNOT DO THIS AND THAT IS A THEOREM RATHER THAN A
MEASUREMENT.** The relaxed iteration has derivative `1 + ω( G′ − 1 )` at the
fixed point, which for `G′ > 1` exceeds one for **every** `ω > 0`:
under-relaxation stabilises a map that oscillates and can do nothing for one
that runs away. Measured, the damped loop walks the error 0.362 → **0.019** and
straight back up to 0.50, **the step never shrinking — including where the error
passes through zero**. It does not stall at the fixed point, it crosses it.

**AND THE NEWTON IS AFFORDABLE BECAUSE THE PROFILE IS FITTED.** `g²` is carried
as a handful of least-squares coefficients rather than a table, so a differenced
Jacobian costs a few map evaluations a step against `INVERSION-PLAN.md` §11.1's
**5.7 hours** for `dGeometry_dpsi`. The fit was put there for conditioning — the
innermost surface's inverted `g` moves **13% for a 1.3%** change in the profile,
where the extraction is least reliable — and it pays for the Newton as well.

**KINSOL RATHER THAN A HAND-ROLLED ITERATION, FOR THE LINE SEARCH.** An undamped
outer Newton takes a first step the **inner** bordered Newton cannot solve at,
and that inner solve has no globalisation of its own —
`GradShafranovSolver` refuses every `Globalisation` but `None` while `ψ_ax` is a
border unknown. So the outer step length is the only control there is. Writing
the search here would repeat `HDG-NPC-GLOBALISATION-FROM-MEQ.md`'s recorded
mistake: a monotone test with no sufficient-decrease constant accepts any small
enough step and creeps instead of failing honestly. `KIN_LINESEARCH` applies
Armijo.

**AND THE MAP MUST BE TOTAL, BECAUSE KINSOL IS C.** An exception raised inside
the residual unwinds through its frames and denies the line search a finite
value at the trial point. With the map throwing, the first full step lands on
`g² < 0` across the whole profile and the run dies **without backtracking
once**. An inadmissible state returns a residual pointing back to the last
admissible one instead; it can introduce no spurious root, `g² < 0` not being an
equilibrium.

**AND THE LINE SEARCH HAS A SECOND PRECONDITION, WHICH IS THAT THE JACOBIAN IS
NOT SINGULAR — AND FAILING IT DOES NOT FAIL, IT HANGS.** `KIN_LINESEARCH`
interpolates its step length on a quotient whose numerator and denominator both
carry the directional derivative `⟨F, J p⟩`. A singular `J` makes that zero
**whatever the step is**, so the quotient is `0/0`, the iterate goes to NaN —
and **KINSOL then never returns**, because every one of its convergence and
failure tests is a comparison against NaN and every comparison against NaN is
false. Measured on a two-variable map with no root at all: still calling the map,
at NaN, after **two million evaluations**, which on the real loop is two million
equilibrium solves. It was found as a test case that ran for **36 minutes**
without finishing.

**NEITHER OBVIOUS REPAIR REACHES IT, AND THAT IS THE PART WORTH KEEPING.**
Bounding the linear solver's iteration count does not: `mfem::GMRESSolver`'s own
back-substitution is `y( i ) /= h( i, i )`, unguarded, and on a singular system
`h( 0, 0 )` is zero, so it comes back **infinite on the first iteration** and a
cap changes nothing. Replacing it with a rank-revealing dense solve does not
either — **even though that is correct and returns a step of exactly zero** —
because it is `⟨F, J p⟩` and not the step that has gone to zero. Both repairs
were built and measured before that was understood. **The line search cannot be
rescued from outside it**, so the degeneracy is detected *before* KINSOL is
entered, and the Jacobian is cached on its point so KINSOL's own first
evaluation reuses it and the check costs **no map evaluations at all**.

**A ZERO JACOBIAN IS ALSO WHAT AN ALREADY-SOLVED PROBLEM HAS, AND THEY ARE
OPPOSITE ANSWERS.** `G = identity` makes every point a fixed point, so returning
at once is right there where refusing would be wrong — and the two are
**indistinguishable by their Jacobians**, both exactly zero. The *residual* is
the discriminator, and a refusal keyed on rank alone refuses both. Both cases are
asserted, side by side, for that reason.

**THE ASSERTION IS A COUNT AND NOT A TIMING**, on the same principle as the
symbolic-factorisation reuse: a wall clock is a measurement about the machine
where the map-evaluation count is a measurement about the code. A map with no
root must be refused within one Jacobian and its residual — **5 evaluations, and
the case now takes 393 µs**. The dense truncated-SVD solve is kept regardless,
being the right solver for a small dense system and the thing that bounds a
Jacobian going singular partway through, where the pre-flight check cannot see
it.

**THE ONE THAT LOOKED LIKE PHYSICS AND WAS ARITHMETIC.** The fit's knots must
span the whole of `Ψ`, not the family's own range. `FluxSurfaceFamily` refuses to
extrapolate and is right to — outside the cut a traced surface is made of
something else — but **a fit is a model, not a measurement**, and evaluating it
at `Ψ = 0` or `1` is what a model is for. Laying knots only over `[ 0.05, 0.95 ]`
and letting `SplineProfile` clamp gives a profile differing from the closed form
exactly where it clamps, so **the outer residual AT THE ANSWER was 1.3e-01
rather than zero** — and both outer methods then converged, correctly, to a
fixed point that was not the answer, from 40% away and from 5% away alike,
landing on the same wrong one to six digits.
`theOuterResidualsConditioningAtTheAnswer` found it and now guards it:
**4.9e-13**, on a Jacobian whose degeneracy measure is 0.12.

**A CONVERGED OUTER SOLVE IS NOT EVIDENCE THE ROOT IS THE ANSWER**, which is the
transferable part. Both failing configurations converged — to 3.5e-07 — and
reported success. What separated them was measuring the residual at a point
known independently to be right.

**AND IT REACHES THE SOLVE FROM A FILE**, as `[source] SafetyFactorFile`, with
`SafetyFactorDegree` and `ToroidalFieldGuess` beside it.
`examples/q-driven.toml` is the worked example and
`DriverAcceptance::theDriverSolvesForTheToroidalField` is the acceptance:
**42.5 s, 12 outer iterations, 43 equilibria**, recovering a closed-form
`g = 2.20 + 0.55Ψ` to **5.0e-06** in the worst of the three coefficients of
`g²`.

**THE TARGET IS MEASURED RATHER THAN INVENTED, WHICH IS THE WHOLE OF WHY THAT
IS AN ACCEPTANCE.** `examples/q-driven-q.dat` is the `q` of an equilibrium built
from that closed form on that box **and that mesh**, taken through
`q = V′g⟨R⁻²⟩/4π²` at 24 surfaces — so the answer is known independently of the
loop, and the discretisation error is common to the measurement and to the
inversion, which leaves what the case reads the LOOP's accuracy rather than the
mesh's. A target invented instead would only have shown that the loop reaches
**some** fixed point, which is exactly the failure the knot-span defect
produced.

**ONE SOLVER SERVES EVERY MAP EVALUATION**, through
`NormalisedMHDSource::setGGPrime()`. Only `gg′` changes between steps, so the
mesh, the spaces, the forms and the trace solver's symbolic factorisation all
survive, and each solve is warm-started from the previous one's field —
rebuilding the solver per evaluation, which is what the library test does, would
throw all of that away 43 times. **The normalisation is deliberately not reset
by that setter**: `ψ_ax` and `ψ_bnd` are unknowns of the bordered Newton and
belong to the solve rather than to the profile.

**AND THE SOLVER IS NOT HOLDING THE ANSWER WHEN THE LOOP RETURNS**, which is
the trap in wiring this to a driver. KINSOL's last call to the map is whatever
it needed last, and on a converged run that is usually a **differencing
column** — an equilibrium one step off the answer in one coefficient. Every
output reads the solver, so the driver re-solves at the converged coefficients
before writing; without that it would publish a perturbed equilibrium beside
converged coefficients, agreeing with them to the differencing step and to
nothing else.

**`meq::ToroidalFieldMap` IS THE PICARD STEP AS AN OBJECT, AND IT IS MFEM-FREE.**
The map is `c → gg′ → [ solve, extract ] → invert → fit → c`, and only the
bracketed half needs a mesh — so that half is a **callable the caller supplies**
and the rest is CI-gated arithmetic beside the inversion it wraps. It owns the
two obligations that are easy to miss and were both hand-rolled in the test
first: that the map is **total**, and that the target is asked in the SOURCE's
`Ψ` while the family is labelled in `Ψ_N`. A null return is how the caller says
"did not converge" without an exception.

**THE `.nc` CARRIES THE ANSWER**, as `toroidal_field_driven`,
`g_squared_coefficients` and `safety_factor_target`. On this route `g` is what
the run **found**, so there is no `GGPrimeFile` beside the output for a consumer
to look it up in — the same argument that put `axis_normalised_flux` there.

**WHAT IS LEFT IS THE FIXTURE AND NOT THE MACHINERY**: the example is a
rectangle with one closed plasma rather than a machine, and nothing has yet
driven `q` on the limited tokamak, where the support moves and `ψ_bnd` is an
unknown.

### IN-5: an open surface is not a loop, and a periodic basis cannot fit one

**`ContourTracer::traceOpen()` and `meq::fitOpenSurface()`.** Every surface
IN-0 to IN-6 handles is a loop — the tracer closes it, the disc chart labels it
by an angle about the axis, the averages integrate round it. **A level outside
`ψ_bnd` is none of those.** It runs to the wall and stops, so it has two genuine
endpoints, no period, and no enclosed area.

**IT WAS DEFERRED WITH FREE BOUNDARY AND IS NOT ANY MORE**, for the reason the
deferral gave: a fixed-boundary problem has no level outside its own boundary to
trace, and free boundary now solves.

**THE TRACE IS BOTH DIRECTIONS FROM THE SEED, JOINED.** `trace()` follows one
and comes back with HALF the curve wearing the label `LeftMesh` — which is what
`v0-legacy:FluxSurfaces.cpp` returned as a contour, printing *"Terminating
because curve left domain"*. The second half is the same march with the tangent
negated, which is **exact** rather than approximate: the tangent is
`( −q_z, q_r )/|q|`, so one sign reverses the whole march and leaves the curve
identical. Not a second tracer over `−q`, which does the same thing at the cost
of a copy of the field and a second set of settings to keep in step.
`ContourStatus::Open` is a **success**; a level that loops comes back `TooLong`,
the closure gate being suppressed, and the caller is told to use `trace()`.

**THREE THINGS IN THE JOIN THAT ARE SILENT IF WRONG.** The reversed half's
tangents must be negated too, or the Hermite interpolation folds back at the
seam while every point sits exactly on the level set. The arc length is
re-accumulated from one endpoint rather than patched, since each half measured
its own from the shared seed and the seam would otherwise carry two origins. And
`turning` is left at **zero** rather than summed: it is the winding of a closed
curve, and on an open arc the two halves measure it from opposite orientations,
so a sum is a number with no meaning that a reader would take for one.

**THE BASIS IS CHEBYSHEV IN NORMALISED ARC LENGTH**, `t = 2s/L − 1`, and the
measurement is on an analytic fixture the tracer is pointed at directly —
`ψ = z − a( r − r₀ )²` in an H1 space of degree 2, which represents that
quadratic **exactly**, so nothing below is the discretisation:

| modes | 4 | 12 | 20 | 28 | 32 |
|---|---|---|---|---|---|
| **Chebyshev**, off-sample | 3.31e-02 | 3.01e-04 | 6.97e-06 | 1.28e-07 | **6.24e-08** |
| **periodic control** | 2.41e-01 | 3.02e-01 | 3.07e-01 | 3.11e-01 | **3.12e-01** |

**5.3e+05 against 0.77**, and the control gets *worse*. Measured at 401 points
the fit never saw, against the closed form — a residual at its own samples is
what a least squares minimises and would flatter both columns.

**THE RATIO IS WHAT SAYS GEOMETRIC, NOT THE TOTAL.** The steps are equal and
**additive**, so a geometric `C q^m` holds its ratio while an algebraic `C m^-p`
gives `(16/12)^p`, `(20/16)^p`, `(24/20)^p` — falling, because the multiplier
falls. Measured, flat at **6.37 to 6.43**. It floors at 6e-08, which is the
traced points' own 1.6e-10 through a least squares of 32 columns, and the
assertion stops before it.

**AND THE CONTROL FAILING IS THE POINT RATHER THAN THE CHEBYSHEV SUCCEEDING.** A
periodic basis forces `R( −1 ) = R( +1 )` on a curve whose ends are 0.8 m apart,
so it is **inadmissible** rather than merely worse; given more modes it fits the
samples better while the periodicity pushes the curve further from the truth
between them. Same mode count, same decomposition, same distance measure — the
two differ in their basis and in nothing else.

**WHAT IS DELIBERATELY ABSENT.** There are no flux-surface AVERAGES over an open
surface and there should not be: `V′` and `⟨R⁻²⟩` are integrals round a closed
loop and the volume an open curve encloses is not defined. IN-5 delivers the
**geometry**, which is what a scrape-off-layer consumer asks for. And the tracer
is still not X-point aware — a level AT a separatrix stalls at the saddle, where
the level set is not a 1-manifold, and `Stalled` is the honest answer.

### The disc basis, and why `ρ = √Ψ_N` rather than `Ψ_N`

`src/meq/Zernike.{hpp,cpp}` is the basis `IN-3` fits in, landed early because it
is **MFEM-free** — plain doubles, like `Profiles` and `Source`, so CI can build
and test it without the MFEM branch it cannot obtain.

**The index constraint is the entire point.** `l − |m|` even and `l ≥ |m|` is
exactly what excludes `ρ² cos θ` and friends, which are not smooth at the
origin; what survives is that **every admissible mode is a bivariate polynomial
in `(x, y)`**, so the centre of the disc — the magnetic axis — is an ordinary
interior point and needs no special case. That is a better argument for the
basis than "DESC does it".

**GETTING THE RADIAL COORDINATE WRONG IS SILENT.** `ψ` has a quadratic maximum
at the axis, so `Ψ_N` behaves like (distance)² there and the geometry is smooth
in the *distance*: parametrise by `Ψ_N` directly and every basis converges
algebraically against a square-root branch point, worst near the axis, **with
nothing in a convergence table to say why**. Same species as `CLAUDE_HDGGS.md`'s
*A wrong Jacobian is invisible to a convergence table*. `radiusFromNormalisedFlux()` and
`fluxDerivativeFromRadial()` exist so the `1/(2ρ)` chain factor is not written
by hand at each call site, because that factor is the thing that gets dropped.

**The radial polynomial is a Jacobi polynomial under a change of variable**, so
`Zernike.cpp` calls `boost::math::jacobi()` rather than carrying a recurrence,
and **the explicit factorial sum every reference prints is not used at all**:
its terms are binomial-sized while its answer is `O(1)` — the largest is about
1e10 at `l = 30` — so it costs about `log₁₀(largest term)` digits to
cancellation. Measured at `ρ = 0.83` it disagrees with the Jacobi route by
**4.5e-8 at `l = 30` and 1.3e-4 at `l = 40`**, and it is kept in the tests as a
*control* rather than as an implementation. A flux-surface fit wanting twenty or
thirty modes would be reading noise. `find_package(Boost CONFIG REQUIRED)` and
`Boost::headers` are new dependencies of `meq_core`; the include is confined to
the `.cpp`, so a consumer of the basis takes on nothing.

**A DERIVATIVE CHECKED AGAINST A CENTRAL DIFFERENCE IS FLOORED BY THE
INSTRUMENT, NOT BY THE DERIVATIVE.** An earlier draft of the plan implied the
tolerance would tighten once the derivative was exact. It does not: a central
difference carries its own `O(h²)` truncation, so the comparison sits at
**1.3e-07** however exact the derivative is. Richardson extrapolation,
`(4D(h/2) − D(h))/3`, reaches **1.4e-11**. **That applies wherever this suite
checks a derivative against a difference**, which is several places.

### IN-0: the tracer, and the pairing that decides whether `q`'s tangent is worth anything

`src/meq/FluxSurfaces.{hpp,cpp}` is the predictor–corrector tracer and the
poloidal-angle parametrisation; `FluxSurfaceConvergence.cpp` is the acceptance.
**Fitted path only** — the band between `Γ_h` and `Γ` is deliberately not in it,
and `sampleField()` is the single seam it will be added at.

**THE DEFAULT IS `Potential::PostProcessed`, AND THE REASON IS NOT THE ONE THAT
WAS EXPECTED.** MEQ has two candidate pairs: `ψ_h` with `q_h`, both `k+1`; and
`ψ*` with `q*` from `DarcyForm::Reconstruct()`, where `ψ*` is `k+2`. Rooting
`ψ*` puts the traced curve `k+2` from the true one instead of `k+1` — measured
**60×, 54× and 83×** closer at `k = 1, 2, 3` on the same mesh, at *fewer*
corrector iterations per point.

**But the finding that matters is about the Hermite, and it was not predicted.**
The interpolant is built on tangents from the **flux** and measured against the
level set of the **potential**, and those are the same curve only so far as the
two fields agree. `∇ψ_h/r` agrees with `q_h` only to `O(h^k)` — differentiating
an L2 potential of degree `k` loses an order while `q_h` keeps `k+1`. So down a
`Δs` sweep the `ψ_h`/`q_h` pair is fourth order **until the tangent tilt takes
over and second order afterwards**: 3.809 → 1.400 → 1.569, tilt 3.5e-5 to
7.3e-5. The `ψ*`/`q*` pair has a tilt a full order smaller, 3.1e-7 to 5.9e-7,
and holds 3.960 / 3.996 / 3.812 across the whole sweep.

**So the post-processed pairing is what makes the cubic-Hermite-from-`q` claim
true at all on this discretisation**, rather than merely making it more
accurate. The control that identifies the tilt as the cause is a third column
built on `∇ψ_h`, the *exact* tangent of the curve being measured against, which
stays fourth order in both pairings — the right tangent for the representation
error and the wrong one for the field error, which is why all three are printed.

**AND THE USUAL REASON FOR PREFERRING `ψ*` IS WRONG.** It is not that the local
post-processing is built so `∇ψ*` matches `r q*`. MFEM's own documentation is
explicit: the constraint equation projects the **total** flux onto the face
restriction of an `RT_k` space, and *that* field is the source term of
Stenberg's local problem — so `∇ψ*` and `r q*` are different objects. What is
true, and measured, is that they agree an order better than the raw pair does.
**Right answer, wrong mechanism**, and the difference matters because the wrong
mechanism predicts exactness and would have made the measured tilt look like a
defect.

**THE FAILURE MODE THAT ENDS A TRACE, AND IT GETS COMMONER AS `Δs` FALLS.**
`{ψ_h = c}` is a union of per-element arcs offset by the face jump, so a point
landing within `jump/|∇ψ|` of a face is on **neither**: the Newton step computed
in element A pushes it into B, B's pushes it back, and the residual alternates
without ever meeting a tolerance tighter than the jump. The band is about 2e-5
wide on a contour of length 1.7 — rare, and *more* likely the finer the spacing,
simply because more points are placed. Left alone it ends the trace, and it
ended several before it was diagnosed. The corrector keeps its best iterate,
accepts after four non-improving steps, refuses to travel more than one
predictor step, and reports `stalledCorrections` against `correctorTarget`.
**One endpoint accepted at the jump level poisons that segment's interpolation
error by a factor of a hundred on a coarse mesh**, so the (c) measurement
excludes those segments as well as the face-crossing ones.

**Closure is not "machine precision", and the correct statement is geometric.**
Both the returning point and the start sit on the level set to 1e-13, but they
are separated *along* the curve by the final step's tangential offset `g_t`, and
an arc departs from its own tangent line by `κ g_t²/2`. Measured **3.573e-10
against a bound of 7.2e-10** — which is the discriminating assertion, because a
drifting tracer's normal error would grow with path length and have nothing to
do with `g_t`. The **residual** is what is flat with path length: 1.559e-13,
1.632e-13, 1.632e-13 over 1, 5 and 10 circuits. Do not compare the closure error
at 1 circuit against 10 — the step controller starts cold, so those two numbers
differ by 36× through their final steps rather than through their lengths.

**The rest, briefly.** `(c)` at fixed `Δs` is flat over a 16× change in dofs —
5.436e-09, 5.441e-09, 5.332e-09 — which is what says (b) and (c) are separated.
`(a)` against the exact contour converges at `k+1` on `Potential::Raw`. The
face-crossing jump converges at `k+2`. And the element walk —
`TransformBack`, then a breadth-first widening over face neighbours — reads
**zero** `FindPoints` fallbacks on every trace, which matters because
`CLAUDE.md` records `FindPoints` as `O(elements × points)` and a corrector calls
for a location once per iteration.

### IN-1: a spectral rule fed a second-order Jacobian is a second-order scheme

The trap `INVERSION-PLAN.md` §3.2 warns about, made into a measurement. Arc
length of one contour, three ways, rates in the number of angles:

→ **[M-47](MEASUREMENTS.md#m-47)** — **from `q`, pointwise**

**The rule is identical in the first two rows.** What differs is that `ρ′` comes
from `q` pointwise in one and from a central difference of neighbouring radii in
the other — and they converge to the **same limit**, so nothing in the second's
own output says it is orders worse. `ρ′ = ρ (u·t)/(u′·t)`, derived from
`ψ(a + ρ(θ)u(θ)) = c` and checked before use.

**Star-shapedness is a hypothesis, is measured, and is refused when it fails** —
`min|u × t|`, the same denominator, reading 0.844 / 0.823 / 0.805 on the
benchmarks. This is `IndexAudit::transversality` again, one stage later.

> **"SPECTRAL IN `N`" IS NOT ATTAINABLE ON A DISCRETE CONTOUR AND THAT IS NOT A
> DEFECT.** `ψ_h` jumps across faces, so `ρ(θ)` is piecewise analytic with jumps
> and no quadrature is geometric on it. The column plunges — 7.24, then 14.8 —
> and floors at about 1.2e-9, which is where the DG jump of `ψ*` converts to a
> distance (6.80e-10). **The control that says the floor is the field and not
> the rule** is the identical rule on the *analytic* contour, reaching
> 3.775e-15.

**`tests/analytic/FluxSurfaceReference.hpp` is the same statement from the other
side**, and it exists because IN-2 needs something to compare against.
Flux-surface averages on an *exact* equilibrium by rays plus the periodic
trapezoid: it reaches 3.11e-15 on the analytic contour, agreeing with the
tracer's own control, and it reproduces the trap at 1.37e-01 → 6.84e-04 for the
differenced metric while the pointwise one is at round-off by 256 angles.

**There are no closed-form Solov'ev flux-surface averages** — `ψ` is elementary,
but `V′`, `⟨R^{-2}⟩` and the safety factor are integrals over a *contour* of it,
and the Cerfon–Freidberg contours have no elementary arc length. What that file
supplies is a converged **reference value**, and the distinction should be said
wherever the number is printed.

**And it carries an identity that needs no reference value at all**: the
flux-surface average of the Grad–Shafranov equation,
`(1/V′) d/dψ ( V′ ⟨|∇ψ|²/R²⟩ ) = −⟨F/R²⟩`. **Quote it with its step or not at
all**: the residual is a property of the differencing step as much as of the
averages, running 9.6e-08 at 5% of `|ψ_ax|` down to 2.3e-11 at 0.6% on the exact
field, and on a *discrete* field it is not even monotone in the step — 2% reads
1.5e-08 where 1% reads 3.4e-07, because the difference divides the surfaces' own
DG-jump noise by the step.
Write the right-hand side with the `F` the solver is fed and **not** as
`−μ₀p′ − g g′⟨R^{-2}⟩`: the second is Solov'ev-specific and is re-derived by
hand, so it is not independent of the hand that derived it. **The `d/dψ` must be
Richardson-extrapolated** — a plain central difference floors the agreement at
8.1e-07 where `(4D(h/2) − D(h))/3` reaches 2.7e-13, which is `Zernike`'s
derivative finding for the third time in this file.

### The band, and the extension that was chosen by measuring both

`ContourTracer::setBandExtension()`. **`BandExtension::None` is the default**, so
the fitted path is bit-unchanged; the curved path is an explicit opt-in.

**THE DECIDING TEST TRACES `Γ` ITSELF**, which lies entirely in the band at every
mesh — `ψ_h` is strictly negative inside `Ω_h`, so every point of `{ψ = 0}` is
answered by the extension, and the exact answer is the curve `D_h` was cut from.
Worst distance from the exact `Γ` over `n = 16/32/64`:

→ **[M-48](MEASUREMENTS.md#m-48)** — `k` · flux Taylor · transfer lift · lift closer at `n = 64`

**The flux Taylor step is second order at every `k` and that is structural**: its
remainder is `O(h²)` over a band of width `O(h)` however good `q` is, so it
cannot improve with the polynomial degree, and it does not — 2.138 at `k = 2`
and at `k = 3` alike. The transfer lift is the error of `q` **integrated along a
path of length `O(h)`**, which is `k+2`. `BandExtension::FluxTaylor` is kept as
the **control**, not as a fallback, and the tests say so in their failure
messages.

**AND THE PRIMITIVE THE PLAN NAMED CANNOT DO IT.** `mfem::PathLiftCoefficient`
`dynamic_cast`s its `ElementTransformation` to `FaceElementTransformations` and
lifts from *that face's own* integration point — it answers "what is `φ_h` on
`Γ_h`", which is `η₅`'s question and not this one. The usable primitive is one
level down and public: **`mfem::PathIntegral( Cu, x, xbar, line_ir )` takes
arbitrary endpoints**, with `mfem::ElementExtension` supplying `E_h(q_h)`.
MFEM's own comment is the licence — fed the exact flux it must return
`p(x) − p(a(x))` *"whatever the path"*.

**THREE TRAPS IN THE MEASURING, AND THE FIRST NEARLY PRODUCED A WRONG HEADLINE.**

* **A contour at fixed `Ψ_N` cannot measure the extension's order.** As `h`
  falls, `Γ_h` climbs toward `Γ`, so a contour a fixed distance inside has its
  band excursion shrink *faster* than `h` — `deep/h` goes 0.92 → 0.66 → 0.26 and
  the band population 136 points → 9. Both columns then converge faster than the
  extension beneath them, the Taylor step reading 3.75 at `k = 3` against its
  true 2.1. Tracing `Γ` is what fixes it.
* **The nearest face of `Γ_h` is not the face you are outside of.** A staircase
  `D_h` cut from a diagonally split Cartesian mesh **pinches** — two triangles
  meeting at one vertex — so both lobes' faces are equidistant from a point just
  outside and the tie goes to loop order. Half the time that picks the lobe
  whose outward normal points the other way and a genuine band point is refused:
  a trace stopped after 85 points of 320. **The outward test must be a filter
  applied first**, with the nearest taken among the survivors.
* **`ψ_h` is not strictly negative inside `Ω_h` to machine precision**, so a
  handful of `{ψ = 0}` points land back inside — one of 322 at `k = 1, n = 64`.
  That is the discretisation creeping above its own imposed datum near `Γ_h`,
  and the acceptance asserts ≥ 98% band rather than 100%.

**The flag is per point and every mask is asserted against `FindPoints`**, in
both directions, because `CLAUDE.md` records that a *count* rather than a mask
was the other half of a real defect in the `.nc`. `ContourPoint::extended` and
`bandDepth`, `Contour::extendedPoints` / `deepestBandPoint` / `bandExtension`,
and `AngleParametrisation::extended` per node, so a consumer can report
band-crossing surfaces separately.

### IN-2: the averages, and `ψ*` does NOT buy an order here

`src/meq/SurfaceAverage.{hpp,cpp}`. **One primitive with two builders** — an
integrand in `(R, z, ψ, q)` over either the angle parametrisation or the traced
contour with Gauss points — and every named quantity a one-line wrapper. That
shape is deliberate: `MANTA-COUPLING.md` says the slot list "is negotiated with
the transport physics case", so a list that is still moving must not be a list
of functions. Conventions, which are part of the definitions:
`V′ = ∮2πR dl/|∇ψ|` and `⟨X⟩ = (1/V′)∮2πR X dl/|∇ψ|`. **The safety factor is
`safetyFactor(g)` and never `q`**, `q` being the flux.

**`ψ*` DOES NOT BUY `k+2` IN AN AVERAGE AND THE REASON IS STRUCTURAL.** Both
pairings converge at **`k+1`** — `V′` at 2.230 / 3.185 / 4.296 raw against
1.957 / 2.827 / 4.265 post-processed. The weight is `2πR dl/|∇ψ|` and
**`|∇ψ| = r|q|`**: the reconstruction buys its extra order in the *potential*,
and there is no `k+2` flux to divide by. The level set improves, the weight does
not, and the average inherits the worse. What `ψ*` buys is a **constant** —
×1.29, ×1.46, ×1.73 in `V′`.

**That is the same shape as the band continuation of `B`**: a quantity limited
by the one factor with no solved variable behind it. Third time in this item
that the answer turned on *which field the error divides by* rather than on
which field was rooted.

**AN AVERAGE DOES NOT ESCAPE THE METRIC TRAP, AND THE PLAUSIBLE ARGUMENT THAT IT
DOES IS WRONG.** A ratio looks as though it should cancel a bad metric, the same
weight appearing above and below. It cancels a **constant** — about 40× — and
**nothing in the order**: from `q` pointwise the sequence rates are 6.77 for `V′`
and 7.38 for `⟨R^{-2}⟩`; with the metric differenced, 1.92 and 1.80, a separation
of 1.9e+06 by 128 angles.

**THE IDENTITY'S CONTROL IS FLAT, WHICH IS THE SHARPEST FORM OF THE RICHARDSON
FINDING YET.** At `k = 2` over a *sixteenfold* refinement the plain central
difference reads 5.8e-05, 7.9e-06, 8.7e-06, 8.9e-06 — **it stops moving at three
figures** — while Richardson goes 8.9e-05 → 1.5e-08 at rate 4.185. A column that
does not converge under mesh refinement is measuring the instrument.

**And the step is not monotone on a discrete field.** At `k = 2, n = 96` a step
of 2% of `|ψ_ax|` reads 1.5e-08 where **1% reads 3.4e-07** — the smaller step 20×
worse, because the difference divides the surfaces' own DG-jump noise by the
step. There is an optimum; find it rather than assuming smaller is better.

**Two extractions agree to about their own error and CANNOT do better.** The
angle fit and the Hermite contour disagree by an amount converging at the
field's own order — 2.96/3.09 at `k = 2`, 4.14/4.48 at `k = 3` — because
`{ψ_h = c}` is a union of per-element arcs offset by the jump and two routes
placing nodes differently sample different arcs. **So the assertion is on the
rate**, not on the gap being small. A missing `2πR`, a metric about the wrong
point, or a gradient on the wrong side of the division would each break the
rate, and no single-route table could see any of them.

**The fixture needed its own box, and one surface is not measurable at all.**
`standardBox()` cannot hold these surfaces — `Ψ_N = 0.25` on `nstx()` already
spans `r ∈ [0.99, 1.57]` against a box ending at 1.4 — so the study runs on
`[0.60, 1.90] × [-1.10, 1.10]`. `Ψ_N = 0.75` leaves under one cell of margin at
the coarsest mesh of a dyadic sweep, so an `h`-study there measures the
contour's distance to the mesh boundary; its reference value is asserted
instead. **That is the fixture's elongation, and it is one more reason the
curved path is where this item actually lives.**

### IN-3: the fit, and the angle is not free

`src/meq/SurfaceFit.{hpp,cpp}` fits `R` and `z` as truncated Zernike expansions
to a point cloud of traced surfaces. **MFEM-free**, like `Zernike`, `Profiles`
and `Source` — plain doubles in, coefficients out — so CI can build and test it;
a caller writes a two-line loop to turn an `AngleParametrisation` into samples,
and `CMakeLists.txt` records that reason beside the source list.

**THE GEOMETRIC POLOIDAL ANGLE IS INADMISSIBLE, AND EVERYTHING ELSE IN THE STAGE
IS DOWNSTREAM OF IT.** Labelling surfaces by the geometric angle about the axis
— the obvious choice, and what `AngleParametrisation` produces — **makes the
disc map non-smooth at the axis for any non-circular surface**, so no basis that
is smooth there can converge against it.

The argument is exact, and it takes one line to check. A function smooth at the
origin has precisely **one** angular harmonic multiplying `ρ¹`. Take nested
ellipses of semi-axes `a`, `b` labelled by geometric angle; with `u = ρ cos θ`
and `v = ρ sin θ`,

```
x( u, v ) = a b · u · √(u² + v²) / √(b² u² + a² v²)
```

whose trailing factor is **homogeneous of degree zero** — a function of
direction alone, with no limit at the origin unless `a = b`. The `ρ¹`
coefficient therefore carries `cos θ, cos 3θ, cos 5θ, …`, and every harmonic
past the first is one the Zernike index constraint **excludes**, precisely
because it is not smooth there. The parametrisation puts content exactly where
the basis refuses to look. Measured on nested ellipses, where the answer is
known exactly, that column decays like `L^{-1.2}` and **never converges**:
1.53e-01 at `L = 2`, still 1.07e-02 at `L = 20`, against 1.8e-14 relabelled.

**The repair needs no field.** Near the axis every equilibrium's surfaces are
ellipses, so a three-parameter fit of a quadratic form to the *innermost traced
surface* recovers tilt and flattening, and the relabelling is an exact
reparametrisation of the circle. `meq::relabelByAxisShape`, worth **45× to
660×**. On `nstx()` it recovers short/long = 0.4943 from samples against 0.4845
from the Hessian — two independent routes to the same shape.

**`ρ = √Ψ_N` IS NOW MEASURED RATHER THAN ARGUED**, which is what the control was
for. Worst fit error **3.44e-05 against 7.01e-03** for the same points
parametrised by `Ψ_N`, a factor of 204; the envelope over `l = 10 → 20` falls by
32.6 against the control's 4.5. **Conditioning is untouched by the choice**
(1.87e3 against 1.44e1), so it is not a conditioning artefact — it is the
square-root branch point at the axis, exactly as `Zernike.hpp` claims.

**The parity control fits the sample cloud eight times better and is useless**:
condition number **9.15e+16**, axis error 1.23e-04 against 3.56e-07, and an axis
that moves by 4.8e-03 depending which `θ` you approach along. A better residual
on the data you fitted and nothing anywhere else.

**The axis comes out `θ`-independent for free and exactly** — spread
`0.000e+00` at every degree, on analytic and discrete data alike, because every
mode with `m ≠ 0` carries `ρ^|m|` and above. Asserted as an exact zero, not a
tolerance. **A caveat for anyone rebuilding that control**: a tensor product
that keeps `l ≥ |m|` *also* gives an exactly `θ`-independent axis, because those
modes vanish at `ρ = 0` too. Only admitting modes with **no radial factor**
breaks it, and a control built the obvious way demonstrates nothing.

**THE LEVER ON CONDITIONING IS THE HOLE, NOT THE SAMPLE LAYOUT.** A sample set
has a hole in the middle — no surface is traced at `Ψ_N = 0` — and the
orthogonality argument needs nodes spanning the whole disc. Condition number
against the inner limit: 7.78 at `Ψ_min = 0.02`, 3.19e+02 at 0.10, **7.30e+04 at
0.25**. Four orders. The three layouts — equispaced in `Ψ_N`, equispaced in `ρ`,
Gauss in `Ψ_N` — agree to within **25%** at every hole size, and Gauss is
sometimes the worst. The test asserts that wrong story dead, at
`layout spread < 2×`.

**And rescaling the disc edge is a change of BASIS, not of model.** A Zernike
expansion of degree `L` spans the polynomials of degree `L` in `(x, y)` and that
space is closed under scaling, so the two extents fit the *same function* —
measured, identical worst errors to better than 1e-6 relative at every inner
limit. The choice is purely conditioning and is worth up to **11,500×**. The
default stays `1.0` so that a change of coordinate is never silent, and
`majorRadiusExpansion()` refuses unless basis, coordinate and edge are all the
plain ones.

**The Richardson finding needed a step sweep to appear at all, and that is
itself the lesson.** At the natural step it is worth **1.6×**, not the 10²–10⁴
seen elsewhere in this tree, because here the *fit's own* derivative error is the
binding constraint rather than the instrument. Swept, it separates properly —
70.5× at a step of 0.16 — and the diagnostic is that **the plain column falls by
50.8 across an eightfold refinement while the extrapolated one moves 11%**: the
converging column is the instrument, the flat one is the answer.

**The number a coupling reads.** `∂(geometry)/∂Ψ_N` grows exactly as the
coordinate demands: the product with `2√Ψ_N` settles at **1.148** as the inner
limit falls to 0.005, so an innermost node at `Ψ₁` carries a geometry derivative
of about **0.574 / √Ψ₁**. The assertion is that the product is bounded and
settling, which is the statement that the growth belongs to the coordinate and
not to the fit. `MANTA-COUPLING.md`'s consumer can keep its nodes off `Ψ = 0`,
so this is a conditioning number for a coupling to read rather than a defect.

**`fitByAngle()` now accepts its best iterate** where the tolerance is
unattainable, with `AngleParametrisation::stalledRays` beside `worstResidual`,
and keeps the throw for a ray that never brackets. The acceptance asserts against
the **measured face jump** rather than a chosen tolerance: 11 of 1536 rays
accepted at their best, worst residual 6.97e-04 against a DG jump of 1.63e-03 —
a ratio of 0.43, i.e. as close as the field allows. `AngleParametrisation` also
keeps the `q` it was computing and discarding, and `sampleAt()` has an overload
reporting `extended`.

### IN-4: the answer is that the question was the wrong one

**The `ψ`-element decision is closed and no element was built.** `§4.4`'s three
candidates were three ways to buy a `ρ`-dependent poloidal angle; **solving for
the angle directly is cheaper than all of them**, and it is what DESC was
measured doing.

`meq::gaugeFreeFit()` requires each disc node only to **land on the right
surface** rather than to sit at a prescribed angle — a geometric Gauss–Newton on
`Ψ_N( x(ρ,θ) ) − Ψ`, warm-started from IN-3's linear fit, with `∇Ψ_N` from the
**solved flux**, `∇ψ = r q`. No force balance and no second solver: MEQ already
has `ψ`. Rows are scaled by `1/|∇Ψ_N|`, so the residual is a **distance in
metres** and the error measure is itself gauge invariant. `SurfaceFit` stays
**MFEM-free** — the field arrives as one callable returning value, gradient and
a refusal.

**On nested ellipses, where the answer is known exactly**, started from the bad
fit: **1.718e-01 → 8.31e-16** at `L = 2`, and 1.623e-02 → 1.86e-09 at `L = 16`.
**That is the theoretically right answer and not a lucky one**: for nested
similar ellipses `x = a u`, `z = b v` in the disc's own Cartesian coordinates, so
the family **is a degree-1 map** under the correct angle, and the solve finds it
in four iterations at round-off. The prescribed-angle fit decays at `L^{-1.2}`
and never converges.

**IN-3's algebraic tail is gone and it holds to `Ψ_N = 0.005`** — linear
1.67e-02 → 3.78e-04 against gauge-free 1.78e-03 → **4.16e-09** over
`L = 4 → 16`, the difference concentrated in the last leg (1.69× against 52×).
**No inner limit tried costs it more than a factor of 1.5.** That is the number
`MANTA-COUPLING.md` needed, and it removes the constraint that motivated
elements in the first place.

**Panici's Figure 5 reproduced on a fit, with no penalty at all**: spectral width
`M(2,2)` falls 1.607 → 1.366 under minimum-norm damping alone, and 1.481 → 1.265
on the discrete field. A solve slides the angle to what its basis represents
best without being asked.

**THE EXPLICIT SPECTRAL-WIDTH PENALTY LOSES ON ITS OWN METRIC.** `M(p,q)` is a
**ratio** of two weighted sums of the same coefficients, so a quadratic penalty
is not a surrogate for minimising it: twelve decades of `λ` move `M` by 1.6% *in
the wrong direction* and cost **44×** in surface error. Hirshman & Breslau
minimise `M` itself, which is not a quadratic problem. Kept as the losing column.

**THE GAUGE IS A SOFT TAIL WITH NO GAP, NOT A NULL SUBSPACE.** Measured, the
ellipse family has **exactly 3** null directions at every degree from 2 to 16 —
of 6 to 306 columns — and **`nstx` has none at all**. What both have is a smooth
tail running to 8e-08 of the largest singular value, with 58 of 306 directions
below `1e-4 σ_max` on the ellipses and 97 on `nstx`. So there is nothing to
project out, and **the floor is a threshold that has to be chosen rather than
read off a gap**. The no-gauge control fires on both fields — first step
`8.6e+10 ×` the coefficient norm, Jacobian `−2.2e+13`, i.e. **folded** — but the
mechanism and the magnitude differ by six orders between them, so a control
measured on one field would have reported whichever it happened to meet.

**AND THE TRUST REGION IS ITSELF A GAUGE, WHICH THE CONTROL HAS TO KNOW.** The
undamped pseudo-inverse reaches round-off at `L ≤ 12` and **fails at `L = 16` and
20**; what makes it robust is adaptive Levenberg–Marquardt damping, which is a
Tikhonov term in disguise. `SurfaceGauge::None` therefore disables the damping as
well as the floor — **a "no gauge" control that kept the trust region would pass
while testing nothing.**

**THE MAP IS CHECKED FOR FOLDING.** Minimum Jacobian **+5.2e-02 to +6.3e-02**
with the gauge on, negative on **every** ungauged run. A surface residual alone
admits a beautiful number over a folded map, which is exactly the class of quiet
wrong answer this file exists to catalogue.

**ONCE THE ANGLE IS FREE, ANY ACCEPTANCE WRITTEN AGAINST A PRESCRIBED ANGLE
MEASURES THE GAUGE RATHER THAN THE FIT.** IN-3's derivative and metric checks
compare the fit's position *at a given `θ`* against a surface traced at that `θ`
— precisely the freedom being granted — so they are **not** re-asserted for the
gauge-free fit. They are replaced by gauge-invariant properties: distance to the
surface, fitted perimeter against exact (3.5e-07 relative), and the sign of the
map Jacobian. **This costs the consumer nothing**, and that is the point:
`MANTA-COUPLING.md` reads flux-surface *averages*, and an average does not know
how its surface was parametrised. The deliverable was gauge-invariant all along.

The conversion product is unmoved at 1.156 against 1.148, the axis spread stays
an exact `0.000e+00` at every degree on both fields, and on the discrete field
the fit sits **5.71e-07** from the exact surfaces — the post-processed pairing's
own `O(h^{k+2})`, so what remains is the discretisation and not the
representation.

### IN-P: what the inversion actually costs, and the cost model was wrong

`tests/performance/InversionScaling.cpp` and `inversion-scan.sh`. **Not a
ctest**, per the standing rule that every number in `tests/performance/` is a
timing and a threaded timing on this machine is a measurement about the machine;
it exits non-zero only for the correctness properties that make its timings mean
anything. Note that **`naming` does not cover `tests/performance`** — it runs
over `MEQ_CORE_SOURCES_PRESENT`, which is `src/meq` alone.

**85% OF THE CHAIN IS `gaugeFreeFit`, AND `INVERSION-PLAN.md` §11 DOES NOT
MENTION IT** — that section was written before IN-4 existed and puts the weight
on the tracer's per-point corrector instead. Measured at `k = 2`, `n = 48`, 12
surfaces × 48 angles, `L = 10`, serial: `trace` 11.2%, `fitByAngle` 1.8%,
`surfaceAverages` 1.8%, the linear `SurfaceFit` 0.6%, `findAxis` 0.01%, and
`gaugeFreeFit` **84.6%**. For scale the solve itself is 0.37 s and
`postProcess()` 0.62 s.

**And inside `gaugeFreeFit`, 82.5% is field evaluation** — 21,888 `sampleAt`
calls — so **about 70% of the whole chain is `ContourTracer::sampleAt`** and the
linear algebra is nowhere near the cost. §11 predicted the quadrature would not
be where the time is, and that much is right: it is 1.8%.

**The corrector is not the problem either.** 97.3% of accepted points take
**exactly two** iterations, mean 1.994, worst 3, with `stalledCorrections` and
`fallbackLocations` both zero on every trace.

**Extraction cost is independent of `k`** — 0.140, 0.141, 0.140 s at `k = 1, 2,
3` — and grows like `1/h`, because the step ceiling is a fraction of the element
size. **In the angle count it is not monotone: more points can be cheaper**,
0.048 s at 24 angles against 0.023 s at 192, because closer spacing makes the
element walk hit instead of missing.

### The largest lever is one integer, and it is free

**`sampleAt` decomposes as 26% walk, 0.9% evaluation, and 73%
`Mesh::FindPoints`** — the last-resort fallback, taken on 184 of 576 points.
Four rings of *face* neighbours reach about four triangles in a straight line and
fewer diagonally, while consecutive ray nodes are one to two cells apart. Swept
on the real code path, with the answers **bit-identical at every depth**:

→ **[M-49](MEASUREMENTS.md#m-49)** — `setWalkDepth()` · seconds · fallbacks

Worth **2.0× on `fitByAngle`, 2.1× on `surfaceAverages`, 1.8× on
`gaugeFreeFit`**, and about **1.57× on the whole chain**, for no change to any
answer. **The default is 12**, and what made that safe to change on a library is
the bit-identity: the walk decides how a point is *found*, not where it is, so
the harness asserts `0.000e+00` across the sweep and the change is provably
free.

**And it is not only a speed question.** With no fallbacks there is no call into
`FindPoints`, which is what makes threading impossible today. Fixing the depth
and giving `traceFromAxis()` its axis element as a hint — `findAxis()` already
knows it — would remove the last unconditional `FindPoints` call and unblock the
shared-tracer construction below.

### Threading: available, and blocked by one non-reentrant function

Given a mesh, a solve and a tracer per thread, the parallelism §11.2 predicted is
there and is exact:

→ **[M-50](MEASUREMENTS.md#m-50)** — threads · over surfaces · over rays

**Every count reproduces serial at `0.000e+00`** — contour points, ray radii,
quadrature weights and `V′` — which is available exactly because independent
surfaces and independent rays reassociate nothing. Surfaces cap at about 3.3×
because there are only twelve of them and the outer ones are longer; rays scale
properly, which is §11.2's asymmetry confirmed: **the tracer's steps are
sequential and the rays are not.**

**What is not available is sharing one tracer**, for the `FindPoints` reason
under `CLAUDE.md`, *Traps*. That is the single blocker, and the walk-depth item above is most
of its cure.

### Two levers from §11.4, one of which does not exist

**Continuation in the flux label is worth 1.004× — nothing — and corrector
iterations went UP.** §11.4 framed this as a genuine trade against parallelism
over surfaces and warned against reasoning from structure. **There is no trade to
resolve**: the predictor for every point after the first already comes from the
previous point of the *same* surface, so continuation can only save
`traceFromAxis`'s bracket and one corrector per surface. Take the parallelism.

**The per-`ψ` cache MaNTA needs is worth `nodes/surfaces` and nothing eats it**:
60 physics nodes over 3 residual evaluations cost 6.55 s naively and 1.29 s
served from one family per `ψ`, a factor of **5.1**.

**AND HERE IS THE NUMBER THAT DECIDES THE COUPLING'S DESIGN.**
`dGeometry_dpsi` by differencing costs one extraction per `ψ` degree of
freedom: **0.447 s × 46,080 dofs ≈ 5.7 hours for a single Jacobian, serial.**
That is what §11.4's shape derivative has to beat, and it is the only part of the
chain where a core count is worth a factor rather than a few percent.

### Eigen, and what the swap did and did not move

`SurfaceFit.cpp`'s hand-rolled Householder QR and one-sided Jacobi SVD are gone,
replaced by `Eigen::JacobiSVD` at both call sites. `find_package(Eigen3 3.3
REQUIRED NO_MODULE)`, **REQUIRED rather than optional-with-fallback**, because a
second numerical path is exactly the maintenance burden the swap exists to
remove; `PRIVATE` to `meq_core` with the include confined to the `.cpp`, so
`SurfaceFit` stays MFEM-free *and* Eigen-free to its consumers.

**`JacobiSVD` and not the faster `CompleteOrthogonalDecomposition`**, for two
reasons: `gaugeFreeFit` applies `σ/(σ² + μ)` and so needs the singular values and
`V` explicitly, and the diagnostics report the design matrix's spectrum — a
rank-revealing QR gives a rank and a solve and no spectrum. It is also the same
algorithm MEQ had, with the same high *relative* accuracy in the small singular
values, which IN-4's soft tail needs.

**`EIGEN_DONT_PARALLELIZE` is set, and for a better reason than the one it was
asked about.** Read from the source, `Parallelizer.h` bails to the sequential
path when `omp_get_num_threads() > 1`, so Eigen does **not** nest inside an
OpenMP region and `setNbThreads(1)` is unnecessary for that. It is pinned because
*outside* a parallel region Eigen would take `omp_get_max_threads()`, and a
threaded GEMM blocks differently at different thread counts — **MEQ's answers
would depend on `OMP_NUM_THREADS`.** `EIGEN_USE_BLAS` and `EIGEN_USE_LAPACKE`
are not set, per the MKL rule.

**IN-4's numbers survive and two moved, reported rather than re-baselined.** The
nstx tail reads 4.163719e-09 against 4.163718e-09, the axis spread is still an
exact `0.000e+00`, and the minimum Jacobian is unmoved to five figures. The
ellipse `L = 2` figure moved 8.31e-16 → 6.80e-16 — both round-off, the headline
ratio reading 2.53e+14 instead of 2.07e+14. **The no-gauge control moved by a
factor of ten** and that is expected: it inverts singular values down to exactly
zero, so its step is dominated by the smallest resolved `σ` and two SVDs report
that differently. Its verdict is unchanged — still folded, still diverging.

### Three things in the kernel survey worth not re-doing

→ **[M-51](MEASUREMENTS.md#m-51)** — measured

**And the reason the third one fails is worth keeping**, because the brief
asserted the opposite: `GridFunction::GetValues( i, ir, vals )` is **a loop
calling `CalcShape` per point, not a GEMM**. It amortises the dof gather and the
temporaries and nothing else. The access pattern is 1.15 points per element in
any case, so there is nothing to batch.

### IN-6: the `(Ψ, θ)` file, the per-`ψ` cache, and a cut that cannot be discovered

**DONE.** `src/meq/FluxFamily.{hpp,cpp}` — **MFEM-free** — holds the
family, the flux label, the interpolation and `meq::GeometryCache`;
`src/meq/FluxExtraction.{hpp,cpp}` is the trace-fit-average loop that needs a
mesh; `meq::FluxGridWriter` in `Output.{hpp,cpp}` writes the file;
`[output] FluxSurfaces` reaches it from TOML. `tests/unit/FluxFamilyTests.cpp`,
`tests/convergence/FluxGridConvergence.cpp` and
`DriverAcceptance::theDriverWritesTheFluxSurfaceGrid` are the acceptance.

**MEQ WRITES ITS EQUILIBRIUM FOUR TIMES NOW, AND THE FOURTH IS NOT THE
EQUILIBRIUM.** The other three are one object at three resolutions — a field on
a domain. `<stem>_surfaces.nc` is a **reduction**: the surfaces themselves and
the integrals over them, against a flux label, which is what a 1-D transport
code reads and what no rasterised `ψ` can be turned into without redoing the
whole of this item at the far end, by a code that does not have `q`.
`tools/README.md` and `docs/output.rst` carry the layout.

**IT IS OFF BY DEFAULT AND THE REASON IS NOT CAUTION.** It costs a contour trace
and an angle fit per surface, and it is the **one output that can be impossible
on a run that solved perfectly well** — a level whose surface is not closed or
not star-shaped has no flux-surface average. A failure there is a warning on
stderr and does not move the exit code, because the equilibrium is already
written and is unaffected.

**THE LABEL IS `ρ = √Ψ_N` AND THE FILE SAYS SO.** IN-3 measured that choice at
**204×** in the worst fit error with the conditioning untouched, so it is the
square-root branch point at the axis and not the algebra. `Ψ_N` is carried
beside it; the interpolation is in `ρ`.

### §8's separatrix cut: nothing gives out, and that is the finding

`INVERSION-PLAN.md` §8's second risk said to decide the cut deliberately
*"rather than discovering it as a convergence failure"*. **The premise is wrong
— there is no failure to discover.** Swept on a curved Miller boundary at seven
levels on two meshes, every trace closed, every fit converged, **zero** rays
stalled, transversality moved 0.730 → 0.684, and `|ψ_h − c|` sat at 2e-13 at
every level up to `Ψ_N = 0.995`. What changes is what the surface is **made
of**:

→ **[M-52](MEASUREMENTS.md#m-52)** — `Ψ_N` · 0.50 · 0.80 · 0.90 · **0.95** · 0.98 · 0.99 · 0.995

with the deepest band excursion **halving exactly with `h`** — 3.07e-02 →
9.84e-03 at 0.95. **The extension answers as confidently as an element does and
the residual cannot tell them apart** (2.04e-13 at `Ψ_N = 0.50` against 2.23e-13
at 0.995), so the per-node mask is the only signal there is and the cut is a
decision about **data provenance** rather than about convergence.

MEQ ships `Ψ_N ∈ [0.05, 0.95]`, configurable, and **refuses rather than
extrapolating** outside it — the opposite of FreeGS, and it is §8's *fail by
throwing* applied to a boundary condition. **The band is `O(h)`, so the right
cut moves with the mesh**: the operational rule is to read `extrapolated`, not
to trust the default.

**And the two ends are cut for different reasons, which is why they are two
numbers.** The outer end is the band. The inner end is that `traceFromAxis()`
brackets a level by walking a ray and near the axis the bracket is a fraction of
an element wide — `dρ/dψ`'s `1/(2√Ψ_N)` divergence belongs to the **coordinate**
(IN-3's product settling at 1.148) and is not what fails first.

### The coarea formula is a third cross-check and it was free

`SurfaceAverage.hpp` records §3.3's implicit quadrature as the missing third leg
of IN-2's cross-check. There is a cheaper one for `V'` alone. Green's theorem
gives the enclosed volume as `V = ∮ πR² dz`, a contour integral with **no
gradient in it at all**, against `V' = ∮ 2πR dl/|∇ψ|`, which divides by the flux
at every node; the coarea formula says `V' = |dV/dψ|`. Nothing in common but the
node positions, and no reference value anywhere. At `k = 3, n = 48`:

→ **[M-53](MEASUREMENTS.md#m-53)** — best over four steps

**The differenced column is flat to three figures across a fourfold change of
step** — 4.71, 4.70, 4.70, 4.70e-05 — which is the metric trap's own signature:
a column that stops moving is measuring the instrument. The Richardson finding
again, and the **third** reading of the metric trap after IN-1's arc length and
IN-2's average, on an integrand neither of them touched.

**And `V'` IS `|dV/dψ|`, WHICH `SurfaceAverage.hpp` COULD BE READ EITHER WAY
ON.** It says both *"`V' = ∮2πR dl/|∇ψ|`"*, positive by construction, and
*"dV/dψ with V the volume enclosed"*, which carries a sign. They agree on
`nstx()`, where the axis is an interior **minimum** of `ψ`; on a fixture whose
axis is a maximum they differ by a sign. Nothing is wrong and nothing moved —
recorded so the next reader does not re-derive it.

### The cache is keyed bitwise, and a hash was rejected on §8's own wording

`MANTA-COUPLING.md` §5: `Geometry` is **pointwise** — once per physics node, per
residual evaluation, handed the whole `ψ` vector each time. §8 then asks that a
stale cache be **impossible** rather than unlikely, and that a served answer
equal a cold one **bit for bit**, zero tolerance, because the last defect of
that kind left a second run completing, plausible, and wrong in the eleventh
digit.

**So the key is the `ψ` vector compared entry by entry with `memcmp`.** A hash
would be cheaper and makes staleness *unlikely*, which is the word §8 rules out.
`memcmp` and not `==`, and the difference is in the safe direction: two vectors
differing only in the sign of a zero compare **equal** under `==` and unequal
under `memcmp`, so `memcmp` recomputes where `==` would serve; and a NaN
compares **equal to itself** under `memcmp`, which is correct and is what stops
the cache defeating itself on a state the consumer is entitled to hand it. It
costs nothing: a served query over 400 dofs is **0.47 µs** against an extraction
of **0.27 s** on the fixture beside it — five orders, and the gap widens with
the mesh, since the comparison is `O( nDOF )` and the extraction is not.

**`resetForRun()` IS NOT REDUNDANT WITH THE KEY, AND THAT IS THE HALF EASIEST TO
MISS.** The key catches a changed `ψ`. It cannot catch a changed **extractor** —
a second run whose `ψ` happens to start where the first ended, against a
different mesh, matches the key and is served the first run's geometry. That is
§8's *"keyed on the object rather than the run"*, and the cache cannot detect it
for itself. The unit suite asserts that a reset forces the rebuild, and — the
part that makes the reset necessary rather than decorative — that **without one
the cache correctly does not**.

**A FAILED EXTRACTION COMMITS NOTHING.** Leaving the previous family under the
previous key serves a different state's geometry the moment the caller retries
at the old `ψ`; leaving it under the new key serves it immediately. Both are the
fudge §8 forbids, so the cache is invalidated *before* the extractor is entered
and the key is committed only on success. And a `ψ` carrying a non-finite entry
is refused **before** the extractor, since §1 guarantees MEQ is called far from
equilibrium routinely and a NaN reaching the tracer is an unbounded search
rather than a message.

**§11.1's `nodes/surfaces` IS CONFIRMED, AND THE BASELINE HAS TO BE NAMED OR THE
NUMBER MEANS NOTHING.** Against a naive that locates the **one** surface through
each node — which is what a consumer without a cache would write, and what IN-P
measured — the saving is `nodes/surfaces`: **3.3× at 40 nodes over 12 surfaces**,
against a predicted 3.33. Against a naive that rebuilds the whole family per
node it is **42×**. Quote the first.

**AND A FAMILY WITH A HOLE IS REFUSED RATHER THAN RETURNED.**
`extractFluxSurfaces()` throws on the first level it cannot reach, naming the
level and its `Ψ_N`. A missing surface would be bridged **silently** by the
interpolation, which is a plausible number for a place MEQ could not look at —
the same failure mode as `v0-legacy:FluxSurfaces.cpp` returning a partial curve
labelled as a closed contour, which §8's risk 3 was closed against.

### What the geometry reads, and what the file does not carry

Against the converged reference on the **exact** field, `k = 2`, `n = 48`, 12
surfaces at 256 angles: `V'` to **1.4e-07**, `⟨R^{-2}⟩` to **8.9e-08**, and every
node within **5.2e-08 m** of the exact contour. Not a rate study —
`SurfaceAverageConvergence` measures the averages at `k+1` and nothing here adds
to it; what this asserts is that assembling a whole family delivers what the
one-surface route already does.

**`safety_factor` IS ABSENT RATHER THAN ZERO WHEN IT IS UNKNOWN.**
`V' g ⟨R^{-2}⟩/4π²` needs `g(ψ) = R B_φ`, and a `meq::Source` carries `g g'` and
not `g`. So the driver writes no such column, and a column of zeroes — which is
what a caller would get from a "sensible default" — is indistinguishable from a
machine with no toroidal field. The file says which by not having the variable.

**AND THE EXTRACTION HAPPENS BEFORE THE VTK STEP, WHICH IS LOAD BEARING.**
`apps/meq.cpp` bends the mesh boundary out onto `Γ` for the `.vtu`, which
changes the map from reference to physical space — and the tracer reads geometry
at every corrector step. Extracting afterwards would trace contours of a field
on a mesh that is no longer the mesh it was solved on. The same reason the
driver already gives for the sampler, one consumer further along.

### What is next

**EVERY STAGE IS DONE** — IN-A, IN-0 (both halves), IN-1, IN-2, IN-3, IN-4,
IN-5, IN-6 and IN-P. What is left is a short list of things the stages left
behind.

**IN-5 IS DONE** — see its own section above. Its in-surface coordinate is
**arc length**, normalised to `[ −1, 1 ]` for Chebyshev rather than to `2π`: a
disc chart has no meaning through a separatrix and an angle about the axis has
none on an open line. Note that arc length does **not** fix axis regularity —
for similar surfaces it is a `ρ`-independent relabelling and buys the same one
order — which is why the two concerns are handled by separate machinery.

**IN-6 is done** — see its own section above. What it did NOT take up is the
other number IN-P measured: evaluating a fit at many points by
Vandermonde-plus-GEMM is worth **34.7×**, and the `(Ψ, θ)` file is exactly the
many-point case. It is not taken because the file is written from the traced
nodes rather than from a `SurfaceFit`, so there is no Vandermonde in the path
at all; it becomes worth taking the moment a consumer asks for the geometry on
its own grid rather than on MEQ's.

**And the number that decides the coupling's design**: `dGeometry_dpsi` by
differencing is **about 5.7 hours for one Jacobian, serial**. The shape
derivative has to beat that, and it is the only part of the chain where a core
count buys a factor rather than a few percent.

**Open, small, and each found by the stage after the one that caused it:**

* **`setWalkDepth()` defaults to 12, and four was not enough** — settled, and
  recorded because the reason is not the timing. Four rings suits `trace()`,
  which steps a fraction of an element, and not the *rays* of a parametrisation,
  which are placed by angle and land one to two cells apart: depth 4 takes the
  `FindPoints` fallback on 183 of 576 rays and depth 12 on none. Worth about
  **1.57× on a whole extraction** with the answers bit-identical at every depth.
  **The threading blocker beside it is closed, and the counter that should have
  caught it could not.** `traceFromAxis()` takes `CriticalPoint::element` as its
  seed hint — one line — and the library's extraction path is then free of
  `Mesh::FindPoints` end to end: `SurfaceAverage.cpp` and `FluxExtraction.cpp`
  reach the tracer only through `traceFromAxis()`, and the two remaining seeds,
  in `trace()` and `traceOpen()`, are called from tests alone and can carry no
  hint — their signatures take a bare point.

  → **[M-71](MEASUREMENTS.md#m-71)** — `FindPoints` calls before · after

  **NO GUARD IS WRITTEN AND THAT IS DELIBERATE.** `locate()` already checks the
  hint twice — `tryElement()` refuses an index outside `[ 0, GetNE() )`, and an
  element is accepted only after inverting **its own** map and finding the point
  `Inside` — so an axis located against a different mesh, or a
  default-constructed one carrying `-1`, costs a failed walk and falls through to
  the `FindPoints` that would have run anyway. The degradation is in the cost and
  never in the answer, and a hint that is merely *near* is the ordinary case:
  `CriticalPoint::overshoot` records that a root beside a face sits slightly
  outside the element that rooted it, and the first ring absorbs it.

  **`Contour::fallbackLocations` CANNOT SEE A SEED, WHICH IS WHY A PER-SURFACE
  FULL-MESH SCAN SAT UNDER A TEST ASSERTING THAT COUNT IS ZERO.** `sampleAt()`
  declares a local `int fallbacks = 0`, hands it to `sampleField()` and throws it
  away, on both the six- and seven-argument overloads — so the public seam reports
  no fallback however it located the point, and
  `theTracerClosesAndTheElementWalkDoesNotFallBack` read 0 before the fix and 0
  after while the real count went 6 → 0. Measuring it needed a breakpoint on
  `mfem::Mesh::FindPoints`. **An instrument that cannot see the thing it is named
  for is worse than no instrument**, and it is the same species as the fixed line
  range that produced a false claim about `ComputeHDGFaceEnergy()`.

  **IT DOES NOT BY ITSELF MAKE `ContourTracer` SHAREABLE**, and saying so is the
  point of measuring rather than asserting. `fitByAngle()` and `faceJump()` still
  reach the counted last resort on real fixtures — 1 in `SurfaceAverageConvergence`,
  17 in `FluxGridConvergence`, 10 in `OpenSurfaces` — and those are data-dependent
  walk failures, not seeds. A shared tracer needs them provably unreachable or the
  fallback made reentrant.
* **`fitByAngle()` does NOT throw where the corrector would accept.** It keeps
  its best iterate, accepts it, and counts it in `stalledRays`; the only throw
  left on that path is a ray on which *every* evaluation left the field, which
  is a real failure and not a tolerance problem. **What that left behind was the
  workaround**: `SurfaceAverageConvergence` carried a tolerance ladder whose
  first rung always succeeded, and which caught `std::runtime_error` — so every
  *other* reason `fitByAngle()` throws was being swallowed, retried at eight
  loosening tolerances, and reported as the ladder's own guess. Removed
 the settled fit stalls on **0 of 2048 rays**.
* **The band mask is per NODE, and making it so was a correctness fix rather
  than a tidying.** The contour builder marks **each Gauss node for itself** through the
  seven-argument `sampleAt()`. Marking a whole segment from its endpoints was
  conservative in one direction and **wrong in the other**: the band does not
  respect the segment a node sits in, so a segment can have both endpoints
  inside `Ω_h` and still cross `Γ_h` in between, which **under-reported** — a
  band quantity presented as a solved one. And the fit builder now samples
  **nothing at all**: `AngleParametrisation` keeps the potential, the flux and
  the band flag its own ray Newton found, so the averages read them instead of
  re-deriving them. `potential` was added beside `fluxR`/`fluxZ` for this, and
  it is deliberately *not* the surface's level — on a stalled ray the node sits
  as close to the level as the field's jump allows and no closer, and this
  records where it actually is.
* **§3.3's implicit quadrature is the missing third leg** of IN-2's
  cross-check. Its acceptance said "all three agreeing is worth more than any one
  being plausible" and two were delivered.
* **Maschke & Perrin is in `tests/analytic/`**, as `MaschkePerrin.hpp` and
  `MaschkePerrinConvergence.cpp`. See `CLAUDE_FLOW.md`, *Maschke & Perrin is
  the second exact rotating benchmark*, and note what it does **not** buy: its PDE is Li & Zhu's
  renamed, so the discretisation study is a restatement. What it buys is the
  **source**.

**What is deliberately absent from the tracer**, so nobody reads more into it
than is there: it follows **one connected component** and neither finds nor
reports a disjoint island at the same level — the same disclaimer
`CriticalPoints.hpp` makes about seeded Newton, and for the same reason. It is
**not X-point aware**: the level set through a saddle is not a 1-manifold and
the tangent is undefined there, so a trace at that level will stall or turn a
corner arbitrarily. `pointAtArcLength()` parametrises segments linearly in
*polyline* length rather than true arc length, and says so.
