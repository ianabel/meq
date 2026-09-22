# Free boundary: HDG inside, an exact exterior operator outside

**`CLAUDE.md`'s free-boundary half.** It was about a third of that file.
`CLAUDE.md` is still the maintainer's index and is authoritative on anything not
about free boundary; `CLAUDE_HDGGS.md` has the discretisation, the sign
conventions, the orderings and the Newton machinery this campaign inherits
rather than owns, and `CLAUDE_INVERSION.md` and `CLAUDE_FLOW.md` are the other
two campaigns.

**The plan is `FREE-BOUNDARY-PLAN.md`** and every section number cited below
(§3.4, §7.16, §7.20, §11.3 …) is its. `ROADMAP.md` item 1 is where free boundary
sits in the priority order.

**Measurement tables live in `MEASUREMENTS.md`** under stable `M-nn` anchors,
reached by the `→ **[M-nn](MEASUREMENTS.md#m-nn)**` lines below. **Do not
renumber the anchors.**

**The standing rules are `CLAUDE.md`'s.** Every stage ends at a **measured**
number rather than at "it runs". A hypothesis that was measured and falsified
stays, because recording it is what stops it being re-derived. The history of the
file itself does not: present tense, no "used to", no dated "now fixed".

## Free boundary: MEQ SOLVES A MACHINE CASE, AND IT AGREES WITH freegs4e

**`ROADMAP.md` item 1** and the plan is `FREE-BOUNDARY-PLAN.md`. **MEQ solves a
free-boundary tokamak equilibrium and reproduces an independent code's answer.**

**THE CASE IS `freegs4e`'s LIMITED CIRCULAR TOKAMAK** — vertical-field coils
only, so no X-point exists in range and the boundary is set by a limiter, which
is the one configuration MEQ's pointwise plasma-support test can represent. Same
four coils, same profile shape, same prescribed `I_p = 300 kA`, same limiter
point. MEQ solves it on a **gmsh half-disc with the conductors meshed to**, with
the exterior DtN on a semicircle at `ρ_Γ = 2.4`, and with `ψ_ax`, `ψ_bnd`, the
profile scale and ten Gegenbauer coefficients all unknowns of **one** bordered
Newton.

**AGAINST THE CONVERGED 513² REFERENCE, WITH ITS OWN COIL CURRENTS AND ITS OWN
LIMITER CONTACT** — `k = 3`, **26375 elements**, ten modes, 24 Newton steps:

→ **[M-02](MEASUREMENTS.md#m-02)** — freegs4e, 513² · MEQ · apart

both flux values **inside** the 7.92e-04 the reference itself sits from its own
Richardson limit, and MEQ converged in its own mesh to **3.2e-05**. The scale is
a third quantity predicted from the reference's own metadata rather than fitted.
**And on the SHIPPED 129²-consistent fixture** — `k = 3`, **1601 elements**, ten
modes, 12 Newton steps — `ψ_ax = 9.484400e-02` against 9.483141e-02, **1.33e-04**,
falling to **7.9e-05** at MEQ's own limit under refinement at order **2.93**.

Free boundary by a Dirichlet-to-Neumann map against free boundary by von Hagenow
Green's functions; HDG Newton against finite-difference Picard; C++ against
Python. **They share the equation and essentially no code.**

**AND THE TWO CODES DO NOT MODEL THE CONDUCTORS ALIKE, WHICH WAS NOT NOTICED
UNTIL RECENTLY.** `freegs4e`'s default `Coil` is an **exact filament** —
`controlPsi` is `Greens( self.R, self.Z, R, Z )*turns`, a point source, and its
`area` attribute only imposes a current-density limit and never enters the field.
MEQ's four are 0.1 × 0.1 m **rectangles**. So the 1.3e-04 above was reached
*despite* a modelling difference rather than because the models agree, and with
the coils about 0.4 m from the plasma edge the leading finite-size correction
goes as `( w/d )² ≈ 1.6e-02` on the near field — three orders above the quoted
agreement before any cancellation. **Worth measuring rather than assuming
small**, and `meq::CurrentFilament` (FB-7, plan §7.19) is what would measure it:
match MEQ's conductor model to the reference's, and the difference between the
two MEQ runs is the finite-size effect on one code with one mesh.

**AND IT WAS MEASURED THE OTHER WAY ROUND INSTEAD, WHICH IS THE BETTER
DIRECTION.** Rather than degrading MEQ's rectangles to filaments, the REFERENCE
is rebuilt with finite conductors: `tools/freegs4e-benchmark/shaped.py` replaces
every filament and solenoid with a `freegs4e.shaped_coil.ShapedCoil` on
`conductors.py`'s own rectangle — the same one MEQ meshes — and `fgsref.py
--shaped` is the switch. All seven diverted machines converge with it. On the
diverted DIII-D comparison it is worth **5.0× in relative `L2` and 7.2× in
`L∞`**, and **about 2× on the nodes outside the conductors as well**, which is
the part this section did not predict: the `( w/d )²` correction does not stop
at the conductor's edge. → **[M-96](MEASUREMENTS.md#m-96)**.

**AND THE FIRST DIRECTION IS NOW MEASURED TOO, SO §0a's 2 × 2 IS COMPLETE.**
`COIL-SUBTRACTION-PLAN.md`'s split makes MEQ able to carry the reference's OWN
filaments, at the reference's own positions, and the driver now meshes without
them — so `examples/machine-f-diiid.toml` under `[conductors] Model =
"filament"` is the same machine `F_diiid_conventional` is, rather than a
rectangle approximation to it. It reproduces that reference's `psi_axis` to
**2.8e-05**, its magnetic axis to **0.3 mm** and its X-point to **0.64 mm**, in
two Newton steps on **1510 elements** against the meshed route's 4848.

**THE FIELD-WIDE NUMBER IS WHAT SETTLES [M-111](MEASUREMENTS.md#m-111).** That
entry records DIII-D's relative error flat at 5.785e-03 over a 16× range in dofs
and concludes the conductor model is what it is. Matching the model on MEQ's
side reads **2.031e-03 against the meshed route's 4.593e-03** — 2.26× lower on
3.2× fewer elements, and below a floor M-111 says no refinement buys down.
→ **[M-154](MEASUREMENTS.md#m-154)**.

**AND THE FULL COLD RACE ON THAT MATCHED MACHINE IS RUN, AND MEQ LOSES IT.**
`race.py F --filament` on a quiet machine: `freegs4e` takes **4.9 s at 129² and
12.9 s at 257²** where MEQ's converging rungs take **80 to 256 s** — 6.2× the
257² reference at best and 19.8× at the rung that gives the best answer.
`race.py`'s own docstring predicted *"a race where one arm is both faster and
more accurate"*, and the accuracy half is the only half that survived: the
prediction was inferred from M-154's **dofs**, where nothing compares work per
unknown. → **[M-155](MEASUREMENTS.md#m-155)**.

**AND THE RACE IS NOT LOST ON THE DISCRETISATION. IT IS LOST ON A SOLVE WHOSE
ANSWER IS THROWN AWAY.** Profiled, the `k2r1` rung spends **92.5% of its wall in
plasma-support sweep 0** — 114 Newton steps ending at a `psi_ax` 9.8% from the
answer and an X-point 2.7e-02 m from it, all of it superseded the moment the
support is re-decided. The eleven Newton steps that produce the answer cost
**7.69 s** against `freegs4e`'s **3.14 s**, which is 2.4× rather than 19.8×.
→ **[M-159](MEASUREMENTS.md#m-159)**.

**THE COST IS NOT WHERE A FINITE-ELEMENT CODE IS SUPPOSED TO BE EXPENSIVE.**
`gradient + trace factorisation + backsolve + NPC` is **15.9%** of the solve.
**`re-assembly` is 51.4%, at 2.11 cores** — 8.14 line-search trials per Newton
step against [M-138](MEASUREMENTS.md#m-138)'s warm 1.33 — and the axis is
relocated once per residual at 19.1 ms against the 12.9 ms residual it decorates.
**A share is a property of the problem**: M-138 sizes `re-assembly` at 5.6% and
`postProcess` at 18.5%, and cold and bordered those read 51.4% and 0.6%.

**WHY SWEEP 0 CRAWLS, AND IT IS THE STEP DIRECTION RATHER THAN THE TOLERANCE.**
The full Newton step raises the merit by **1.23 to 1.37× at every iteration for a
hundred of them**, while the field residual across the same trial ladder falls
33×; the line search halves up to twelve times and the best trial beats the
incumbent by 0.2%. Then iterations 112 to 114 read
`1.554e-05 → 9.986e-11 → 9.124e-16` — three steps of ordinary quadratic Newton
once it arrives. `NewtonRelativeTolerance` from 1e-10 to 1e-4 takes sweep 0 from
114 steps to **103**, which is the control that says the last digits are not what
is being paid for.

**AND THE SUSPECT IS NAMED: THE FIRST SUPPORT FREEZE IS POSED WITH `psi_bnd`
WRONG BY 2.26×.** `apps/meq.cpp`'s `edgeFluxOf()` reads `psi` at the bounding
point straight off the grid function, which under a subtracting conductor model
is `psi_p`; every consumer pairs it with a **total** `psi_ax`. On DIII-D `psi_c`
at the X-point seed is −8.945e-02, so the first sweep freezes the support at a
contour well inside the plasma — 517 elements against the 1159 the answer
carries. **The loop repairs itself by sweep 1**, since `axis` and `boundary` are
re-read from `psiAxis()` and `psiBoundary()` thereafter and those are totals, so
it never produced a wrong answer on this rung and presented purely as cost.
Repaired, the whole filament ladder converges in single-figure Newton steps per
sweep, `k3r0` converges where it used to fail under both routes, and `k1r0`
stops reporting a different equilibrium. → **[M-160](MEASUREMENTS.md#m-160)**.

**THE FROZEN SUPPORT'S OWN COLUMN DOES NOT SAY WHAT IT LOOKS LIKE IT SAYS.** `517/517` is `plasmaComponentElements()` over
`plasmaCandidateElements()` — the connected fill keeping every element with
`Psi > 0` — so it reports that the fill removed **nothing**, not that the mesh is
full: 517 of **6134** elements carry plasma at the cold guess. Read as a series
it says the plasma more than doubles on the way to the answer, 517 → 1082 → 1163
→ 1159, and sweep 0 is frozen on the smallest of those for all 114 iterations.
`[solver] PlasmaSupportSweeps` is the outer loop and on this case it **does not
settle** in four.

**THE ERROR STILL DOES NOT REFINE, AND NOW THE CONDUCTOR MODEL CANNOT BE THE
REASON.** 9.533e-04 → 9.473e-04 across a 3.3× range in dofs and two degrees.
M-111 met this floor at 5.785e-03 with the models differing and M-139 at 7.7e-04
with them matched as rectangles; **the filament diagonal makes the model
difference exactly zero and a floor is still there.** More than half of it is
inside the coil boxes — the off-coil column reads 4.2e-04 — and **whose error
that is has not been established**: `freegs4e` carries a filament as a source on
a finite-difference grid, so its `psi` near one is a discretisation of a
logarithm where MEQ's is the logarithm. The experiment that separates them is
the same run read against the **129²** reference as well, and it has not been
run.

**IT COST THREE MORE BORDERS READING THE REMAINDER**, none of them reachable by
any fixture in the tree, and each named by an AGREEMENT rather than a
disagreement — the reference saves `psi`, `coil_psi` and `plasma_psi`
separately, so the decomposition MEQ's split is a decomposition INTO was there
to compare against term by term. M-154a has the three and the order they have to
be fixed in, which matters: two of the three intermediate states look like
progress.

**Three things the conversion has to get right, each of which silently
converges if got wrong.** The CIRCUIT TOPOLOGY, because the references come from
an INVERSE solve and flattening a circuit into independent coils hands it more
freedom than the machine has. The `control` flag per conductor, which says
whether that current may move at all. And the rectangles must be
`shrink_to_fit`'s, not fresh ones, or the two codes are still describing
different machines — TCV narrows two conductors and MAST-U six.

**The inverse solve then moves the currents**, up to 38% relative on DIII-D's
F3A, so the MEQ side is regenerated from the shaped `.npz` rather than reused.

**AND PUTTING `Γ` BETWEEN THE PLASMA AND THE COILS, WHICH WOULD HAVE MADE THE
CONDUCTOR MODEL MOOT BY MAKING EVERY COIL EXTERIOR, IS IMPOSSIBLE ON ALL SEVEN**
— measured, not argued: the margin is negative at every centre on every machine,
−0.26 m at best. → **[M-95](MEASUREMENTS.md#m-95)**.

**AND IT IS A REGRESSION, NOT AN ANECDOTE.**
`DriverAcceptance::theDriverSolvesALimitedTokamak` drives the whole thing from
`examples/limited-tokamak.toml` — `[[coils]]`, `[boundary.limiter]`,
`[boundary.exterior]`, `[source] PlasmaCurrent`, `ConfineToPlasma`,
`[mesh] File` and a `gridfunction` guess, all live at once — and gates `ψ_ax`,
`ψ_bnd`, the profile scale and the delivered current, plus the reported
`psi_axis` against the peak of the field the run actually **wrote**. A DRIVER
test rather than a library one for a reason: all three defects found building
this case were driver-side, in keys describing a box a file mesh never builds,
and a library test would have caught none of them. **11.0 s wall, 17.0 s CPU**,
and the fixture — the gmsh half-disc with the conductors meshed to, the two
`dp/dΨ` tables and the Green's-function guess — is 216 kB in `examples/`.

**Cost against the reference, both serial at `MKL_NUM_THREADS=1`**: freegs4e
**24.77 s CPU** for 42 Picard steps over 16,641 grid points; MEQ **19.78 s** for
**7 Newton steps** over **137,100** unknowns (76,040 flux, 38,020 potential,
23,040 trace) with a **23,040**-dof trace system factorised seven times. Eight
times the unknowns and a larger factorisation, slightly cheaper in total, because
`O(n³)` Green's-function boundary conditions are what freegs4e pays per step. A
wall-clock ratio between the two is not a statement about either code; accuracy
per unknown is the only column with meaning.

**AND THERE IS A SECOND EQUILIBRIUM, REACHED BY AN INPUT ERROR NOTHING ELSE
CATCHES.** The same configuration can converge — every border at machine zero,
the current delivered to seven figures — to `ψ_ax = 2.734289e+00` against
9.48e-02, because **`ψ_ax` WAS the largest NODAL value of `ψ_h` and nothing said
the largest nodal value is a magnetic axis** — the defect option 3 closed, and
the reason it was worth closing: it latches onto a single dof
spiking at the plasma edge, next to nodal values of 0.98, and the current border
then raises the profile scale by 980 to keep `∫F/r` at `μ₀I_p`. A spurious
`ψ_ax` inflates the span, `Ψ` collapses, and the scale compensates — the three
unknowns conspire, and every constraint is satisfied by the artefact.

**WHAT REACHES IT IS A PROFILE TABLE WRONG BY A FACTOR OF THE SPAN, AND THIS
DEGREE IS NOT THE VARIABLE.** A normalised table holds
`dp/dΨ`; another code's `dp/dψ` arrays are a factor of `ψ_ax − ψ_bnd` out.
Measured with one variable changed and everything else held — same mesh, same
degree, same guess — the mis-scaled table gives **44 Newton steps,
`ψ_ax = 2.734289e+00`, scale 9.807e+02** and the corrected one **8 steps,
9.676040e-02, scale 1.0222**. So `p`-refinement is not the cure and the degree
was never the variable. **The tell is the profile scale**: `O(1)` when the table
is right, and nothing else in the output moves. `FREE-BOUNDARY-PLAN.md` §7.16 is
the record; the second, independent diagnostic — `ψ_ax` against the field's own
maximum, or against `meq::CriticalPointFinder`'s O-point, which IN-A already
built — is being wired.

**What is still open** is everything that decides WHICH equilibrium a
free-boundary solve reports: a guess that is not already the answer does not
reach it (a cold bump wanders for 200 iterations), the toy fixture of §7.14 finds
a wall-hugging annulus at every coil current tried, and §7.18 records half of
§7.14's and §7.15's recorded verdicts as false once the `ψ_bnd` defect was
repaired.

The pieces underneath, all measured:

| | |
|---|---|
| **FB-A** | the axis, below. `tests/convergence/AxisConvergence.cpp` |
| **FB-0** | **DONE, and §3.4 closed it.** `src/meq/ExteriorDtN.{hpp,cpp}` — the Gegenbauer basis, the DtN symbol and the mass, **MFEM-free** so CI gates it. Checked against a current loop built from elliptic integrals sharing no code with it: **1.4e−15** in the trace, **6.9e−14** in the DtN. And now against **CEDRES++'s own boundary form** — eq (3.5), a hypersingular double-layer kept in its double-difference form plus a single layer, elliptic-integral kernels throughout — which comes out **diagonal to 1.15e-10** against a scale of 2.56e-01, its diagonal matching `blockEntry( n )` to **3.20e-09**. A boundary integral against a separation of variables, sharing the equation and nothing else. §3.4 named this as the test that would falsify all of §3; it does not |
| **FB-1** | **DONE.** `P`, the transmission row, `setExteriorDatum()`, and both halves measured: `ψ` at 1.99/2.99/3.99 on the half-disc with the datum given, and the exterior coefficients recovered from the transmission condition to 1.9e-04, converging at 3.30. `tests/analytic/ExteriorMatched.hpp` is the exact answer — **the plan's proposed one, filament loop fields, cannot support an order study at all**, `ψ ∉ H¹` at a point source |
| **FB-2** | `src/meq/Coils.{hpp,cpp}`, MFEM-free, and the acceptance identity `∮(1/r)∂ψ/∂n dl = −μ₀I` at **3.3e−11** on the exact field, so a discrepancy on a solve is the solve. **And the coils are reachable from a TOML file** — see below |
| **FB-7** | **CONDUCTORS OUTSIDE `Γ`, THROUGH THE COUPLING RATHER THAN THE MESH.** `meq::ExteriorCoilSet` and `setExteriorConductors()`. The exterior problem is linear, so `ψ = ψ_coil + ψ̃` with `Δ*ψ_coil = 0` inside `Ω` — **the interior equation is untouched** — and the conductor enters through `Γ` alone, additively and KNOWN, on **both** halves: `ψ_coil` in the datum and `q_coil·ν` in the transmission row. No new unknowns. Rates **1.980 / 3.409 / 4.069** at `k = 1,2,3` with the datum given, against a datum-removed control **148,165×** larger; **one Newton step** on the coupled path, `\|a\|` at 4.076 and `ψ` at 4.486. And the **interior** route as a control on it, on a second geometry: **1.951 / 3.727 / 3.879** between the two routes, control 6306×. See below |

## One argument where two were meant: psi_bnd was zeroed inside the Jacobian window

**FOUND BY FOUR AGENTS DIAGNOSING THE SAME FAILURE INDEPENDENTLY, AND
IT IS ONE LINE.** `GradShafranov.cpp`, inside the Newton loop under the comment
*"LAST, and after the source is put back to s"*:

```
normalisedSource->setNormalisation( s );        // ONE argument
```

`setNormalisation( double )` forwards to `setNormalisation( psiAxis, 0.0 )`. So
from that line until `setNormalisation( s, sB )` was restored ~140 lines later,
**`ψ_bnd` was zero** — and what is assembled in that window is `GetGradient()`
and every plasma-current assembly. All of them were built against
`Ψ = ψ/ψ_ax` on a support of `{ψ > 0}` instead of the iterate's own.

**IT IS INERT WHEREVER `ψ_bnd = 0`, WHICH IS EVERY OTHER BORDERED CASE IN THE
TREE.** That is why it survived FB-3. `HighBetaConvergence`'s FB-3 case *does*
carry `ψ_bnd = 0.509`, but its profiles are **constant**, so `∂F/∂ψ ≡ 0` and the
source term never reaches the Jacobian at all — the published table is
**bit-identical** across the repair, which is the proof that nothing else moved.

**WHAT IT COST, MEASURED FOUR WAYS.**

→ **[M-03](MEASUREMENTS.md#m-03)** — the Newton direction against its own linearised system

So the right-hand side told the current row to **reduce** the current when it was
10% short, and no damping was a descent direction. The line search then sat at
its minimum trial for ever, Armijo failing every step, which is the "floors and
then creeps upward" signature.

**THE ASSUMED ZEROS REALLY ARE ZERO**, checked at the same time: all 66
off-diagonal corner entries the code assumes vanish difference to an exact
`0.000000e+00`. The exterior constraints do not depend on `ψ_ax`, `ψ_bnd` or the
scale; neither normalisation constraint depends on the scale or on `a`.

**AND THREE HYPOTHESES WERE KILLED ON THE WAY, WHICH IS WHY THEY ARE WORTH
RECORDING.**

* **Not conditioning.** The dense border matrix has `cond ≈ 1.0e3` — three digits
  of sixteen — and solves its own system to **8.4e-16**. Equilibration is worth
  36×, not eleven orders. The SI dynamic range never reaches `M`, because every
  border row is already a ratio.
* **Not units.** Non-dimensionalising the whole problem — `μ₀ = 1`, `ψ` scaled to
  `O(1)`, currents and amplitudes rescaled — reproduces the SI run **to ten
  digits**, floor and all. It neither fixes nor perturbs it.
* **Not combinatorial.** In the stalled phase the `ψ_ax` argmax is **frozen** for
  21 consecutive steps, the support creeps one-directionally (39 points in, 0
  out), and the line search returns the identical verdict every step. The control
  is sharper still: the configuration that *converges* is the combinatorially
  noisier one, with the argmax jumping three times and hundreds of support points
  flipping **both ways**.

**A SECOND, SEPARATE DEFECT WAS FOUND IN THE SAME PLACE AND IS ALSO FIXED.**
`augmentedNorm` — the merit the line search and the stopping rule both use — was
never given `constraintL`. The comment three lines above its own call site warns
about exactly this for `ψ_bnd`: *"BOTH constraints, or the line search is blind to
the one it is not told about"*. The plasma-current constraint would have been
**98.5%** of the merit had it been included. It is **not** what caused the stall —
measured, the current-aware merit rises at every damping too — but it is why a
solve delivering `∫F/r` **15.9% wrong** could report a converged-looking floor.

**AND `plasmaCurrent()` WAS PUBLISHING THE `ψ_bnd = 0` INTEGRAL**, so
`thePlasmaCurrentClosesAsABorderUnknown`'s assertion that the current is delivered
to 3e-08 was **checking the solve against the formula it used** — this tree's own
recorded trap, met from the inside. The third site is repaired too, and
`recoverPeak()` now takes `ψ_bnd` as an argument rather than being unable to
express it.

**THAT TEST IS GONE, AND WHY IT HAD TO GO IS THE LAST PART OF THE STORY.** With
the Jacobian repaired its configuration — no coils, a free profile scale, a
prescribed current — **does not converge**: the scale runs away to **27.8** and
`ψ_ax` to 3.8e-01, inflating the plasma until it carries the required current in
a larger, flatter channel. It had only ever "converged" because the wrong
Jacobian happened to pin it near scale 0.084, and it had only ever "delivered"
its current because the quantity it asserted on was the one the defect computed.
**A test that passes because of a defect is worse than no test**, so it is
removed rather than re-based; §7.17 records the state.

**WHAT A CURRENT-CONSTRAINED SOLVE STILL NEEDS is the thing that stops the
plasma inflating**, and freegs4e shows what: its `ConstrainPaxisIp` fixes **two**
quantities — the axis pressure *and* `I_p` — against two free parameters. MEQ's
border fixes one. With one constraint and one scale the counting is right but the
solution is not unique, and Newton has no reason to prefer the physical branch.

## The plasma current is a border unknown, and it is what makes a moving support solvable

`setPlasmaCurrent( μ₀ I_p )` makes the profile **scale** an
unknown of the same bordered Newton and prescribes the current instead — the
profiles give the current's shape, the border gives its size. It is what
CEDRES++ and FreeGS both do, and `FREE-BOUNDARY-PLAN.md` §7.14 is the record.

**WHY IT IS NOT A CONVENIENCE.** With the amplitude fixed, a confined
equilibrium is a non-linear **eigenvalue** problem: `Λ = A/span²` must be an
eigenvalue of the linearised operator **on the plasma region**, and the region is
itself unknown. Scaling `A` changes nothing, because `span` moves with `√A` and
the reaction ratio is amplitude-independent — the same statement recorded for the
high-β source. §7.13 measured the consequence: the coupled solve with a moving
support failed from **everywhere** — three limiter radii, two profile exponents,
warm starts, 200 iterations.

**WITH THE CONSTRAINT IT CLOSES IN 63 NEWTON STEPS**, with `ψ_ax`, `ψ_bnd`, four
exterior coefficients and the scale all unknowns and the support moving, and
delivers `∫F/r = 3.499999970e-01` against the **0.35** asked for. The scale came
out **8.40e-02**, so the border did real work.

**ALL FOUR JACOBIAN PIECES ARE ANALYTIC** — `∂R/∂λ` is the source term over `λ`,
`∂G/∂x` is `∫(∂F/∂ψ)/r φ_j`, `∂G/∂λ` is `(∫F/r)/λ` — so the constraint costs one
backsolve and three element loops, not a second factorisation.

**AND THE FOURTH WAS MISSED, WHICH COST THE RATE AND NOT THE ANSWER.** `∫F/r`
depends on `ψ_ax` and `ψ_bnd` **explicitly**, through the `Ψ` the profiles are
evaluated at, so the current row of the corner block is **not diagonal**. Without
those two entries it converged **linearly**, at a clean geometric contraction of
about **0.8 a step** — the signature `CLAUDE_HDGGS.md`'s *A wrong Jacobian is
invisible to a convergence table* describes. Adding them took the residual at iteration 60 from
**7.06e-08 to 1.69e-10**. **A linear rate on a Newton method is a Jacobian
statement**, and this is the third time that has been the diagnosis in this tree.

**THE ARGUMENT IS `μ₀ I_p`, NOT `I_p`, DELIBERATELY.** Everything here already
speaks in it: Ampère's law reads `∮q·ν = −μ₀I_p`, `outwardFlux()` returns the
left side, and the constraint is assembled as `∫F/r`, which **is** `μ₀I_p`.
Taking amperes would mean knowing `μ₀` here, and one that disagreed with the
source's own would scale two terms of one equation differently and converge, at
full order, to a machine nobody described — the trap recorded for why a coil block
carries no `Mu0` key.

**WHAT IT DOES NOT GIVE ON THAT FIXTURE IS A CORE.** The equilibrium it finds is
a wall-hugging **annulus**: `ψ` rises monotonically from 6.1e-04 at `r = 0.1` to
1.23e-01 at `r = 1.4`, so `{ψ > ψ_bnd}` is the outer shell and `ψ_ax` sits near
`Γ`. Every constraint is satisfied by it — they constrain the current and the
normalisations, and none says the plasma is a core.
`meq::CriticalPointFinder::sweep()` now certifies it: exactly **one** O-point in
the whole domain, a maximum at `( 1.3768, +0.0012 )` with `|q| = 4.1e-18`.

**AND A VERTICAL FIELD DOES NOT SUPPRESS THAT BRANCH**, which is the obvious
guess and is false. `μ₀I = 0 / −0.05 / −0.10 / −0.20` per coil **all four
converge** at
limiter 0.80, to 1e-11 or better, and at limiter 0.90 the pattern is
FAILED/FAILED/CONVERGED/CONVERGED — the same verdicts at `n = 32`, so these are
discrete solutions rather than mesh artefacts. And the vertical field **does not
make a core**: the O-point moves only `r = 1.377 → 1.298` over `0 → −0.20`, and
at `−0.35` the equilibrium flips past any core onto an axis-hugging branch with
`ψ < 0` across the whole midplane. The scale numbers themselves are confirmed
(`−1.48809e-02` at `r = 0.4` against `−1.14329e-01` at `r = 1.4`), and the
currents are **per coil**. `FREE-BOUNDARY-PLAN.md` §7.18 is the full
re-measurement, and §7.16 is what a machine case with CONSISTENT inputs does
instead.

**AND A VACUUM START IS NOT THE ANSWER, THOUGH NOT FOR THE REASON THIS
PARAGRAPH GAVE.** At `F_plasma = 0` there is no plasma, so `ψ_ax` is a maximum
of the *coil* field rather than an axis and `ψ_bnd` is the edge of nothing —
**confirmed**: the only critical points are the two coil centres and the grid
maximum of `ψ_h` is `−6.06e-06`. What is **FALSE** is that `setNormalisation()`
refuses it: measured, the vacuum solve **converges in 2 Newton steps to
7.03e-17** with a span of 3.83e-02, and reports three unknowns that mean nothing.
A solve that succeeds and describes nothing is worse than one that refuses.
**The three unknowns the
continuation would exist to carry are exactly the three that stop meaning
anything at `λ = 0`**, and ramping the profile amplitude is the same thing in
slow motion. Ramping **`I_p`** is the exception and is fine: the plasma is
present throughout, and it is prescribed input rather than something recovered
from a black-box `F`.

**THE DEEPER POINT IS THAT `../freegs4e` NEVER MEETS THIS, AND WHY IS
DECISIVE.** Its `constrain` object is a **control system that solves for the
coil currents inside every Picard iteration** — `control.py` assembles a
least-squares system over the currents against shape constraints (`Br = 0` and
`Bz = 0` at a prescribed X-point, isoflux pairs). **The coils adapt to the
plasma.** MEQ's attempt prescribed currents chosen by eye and asked a cold Newton
to find a core plasma consistent with them, which is the opposite. The cheap fix
is not continuation at all but **choosing consistent inputs** — Shafranov's
`B_v = μ₀I_p/(4πR)·[ ln(8R/a) + β_p + l_i/2 − 3/2 ]` gives the vertical field a
given `I_p` needs, so the conductors and the guess can be made to agree before
the solve rather than found to disagree during it. §7.15 is the write-up.

## FB-5: the exterior coupling is an unknown of the same Newton

`setExteriorCoupling( ExteriorDtN const & )`,
and `solveWithNormalisation()` generalised from a `2x2` border to `(N + 2)`.

FB-1b already solves for the exterior coefficients, by **superposition**: one
full solve per mode, one more for the source, and a dense `N x N` assembled out
of the answers. That is exact and it is available only because the problem is
linear. **A plasma source is not**, and superposition stops meaning anything the
moment `F` depends on `psi`. The border is the same system solved as one Newton:

```
T_m( x, a ) = ( transmission integral of x )_m + blockEntry( m ) a_m = 0
```

with `psi_ax` and `psi_bnd` beside it. Measured on FB-1b's half-disc at `k = 2`:

→ **[M-04](MEASUREMENTS.md#m-04)** — `n` · `\|a − exact\|` · `\|a − superposition\|` · Newton

converging at **3.32** against FB-1b's 3.30. **Two entirely different routes to
the same coefficients** — `N + 1` factorisations and a dense assembly against one
factorisation and `N + 2` backsolves — agreeing at round-off. And **one Newton
step**, which is what says the columns really are constant: the residual is
affine in `( x, a )` and an exact Jacobian must finish it in one.

**AND THE ADAPTIVE LOOP RUNS THROUGH THE COUPLING, WHICH COMPLETES
FB-5 — AND IT FOUND THAT `η` CANNOT SEE THE COUPLING AT ALL.**
`theCoupledSolveSurvivesTheAdaptiveLoop` drives solve → post-process → estimate
→ mark → refine with the exterior coefficients live, `Γ` fixed at `ρ_Γ = 1.5` so
the DtN is the same object every cycle and only `Γ_h` climbs toward it. `η` falls
**2.5109e-01 → 3.2079e-02** over four cycles, `L2(ψ)` 3.58e-03 → 9.33e-04, **0
widened fans** so assumption P.1 holds on a graded `Γ_h`, and **one** Newton step
per cycle because the residual stays affine.

**The exterior coefficients do not move — 1.3194e-03 at every cycle, to five
digits — and `Γ_h` keeps its 34 faces while the element count goes 314 → 640.**
Under *uniform* refinement the same quantity converges at 3.32, so it is not a
mode-truncation floor.

**IT IS NOT A MARKING ACCIDENT.** The 34 elements touching `Γ_h` are 11% of the
mesh and carry **0.00% of `η²`**, rising to 0.07% by cycle 3 only because the
interior improves around them. No threshold would mark them.

**AND `η` IS RIGHT, WHICH IS WHY THIS IS A SCOPE STATEMENT RATHER THAN A DEFECT.**
`η` estimates the **interior** discretisation error, and `η₅` on `Γ_h` compares
`ψ*` against the datum actually imposed — the repair recorded under *A separate
`η₅` problem on the extension path* — so it correctly reports the boundary as
well resolved *for the interior problem*. The coefficients are a **boundary
functional**: a transmission integral over `Γ` reached by extension from `Γ_h`.
Nothing in `η` measures it. Both are true at once.

**THE STALL IS QUANTIFIED, AND `η₆` IS THE CURE.** At cycle 3 the interior error was 9.3e-04 against a frozen 1.3e-03
in the coefficients, a factor of 1.4 from the loop reporting a falling `η` while
`ψ` stopped improving.

`GradShafranovSolver::exteriorTransmissionResidual()` is the indicator and
`ResidualEstimator::setExteriorCoupling()` puts it in `η` as a **sixth term**,
zero unless asked for. What it measures is the pointwise mismatch

```
d( x ) = q_h·ν( x ) − ( 1/r ) Σ_n a_n symbol( n ) C_n( μ )
```

integrated as `h_e ∫|d|² dΓ` per element — `η₃`'s scaling for a flux jump. The
border imposes only its **projection onto the retained modes**, so at convergence
`d` is orthogonal to `C_2..C_{N+1}` and is **not zero**: what survives is the
modes past the truncation and the discretisation error, which is exactly the
quantity that froze.

**AND ADDING IT TO `η` IS NOT ENOUGH, WHICH IS THE FINDING RATHER THAN THE
PLAN.** §7.12 said "added to the marking" and that is what was built first. It
does nothing: `η₆` is **8.58e-04** against an `η` of **2.51e-01**, so its share of
`η²` is about **1e-5**, and a Dörfler competition at `γ = 0.6` never reaches it —
`Γ_h` kept all **34** of its faces for four cycles with the term summed in.

**That is not a threshold to lower.** The two are different quantities in
different units — an interior discretisation error and a boundary functional — so
one sum over both is a comparison with no meaning however it is weighted. The
boundary term marks on **its own distribution** and the two sets are unioned.
`η₆` stays *in* `η` regardless, because the **stopping** rule does have to see
it: a loop halting on the interior error alone would report success with the
boundary unresolved.

**MEASURED, THE SAME LOOP RUN BOTH WAYS**, the only difference being one call:

→ **[M-05](MEASUREMENTS.md#m-05)** — cycle · `Γ_h` faces, off · `\|a − exact\|`, off · `Γ_h` faces, **on** · `\|a − exact\|`, **on** · `η₆`

**8.6× on the boundary functional against a control that does not move at all**,
and `η₆` falls with it.
`theBoundaryIndicatorRefinesGammaHAndMovesTheCoefficients` is that table, and the
**off column is the control**: a boundary term that did nothing would leave both
columns identical and still satisfy every assertion about the on column alone.

**ADDING THE TERM OVERRAN A BUFFER IN A TEST, WHICH IS THE TRAP TO EXPECT AGAIN.**
`EstimatorConvergence` declared `double component[ 5 ]` beside a loop running to
`ResidualEstimator::termCount`, so `termCount` going from 5 to 6 wrote off the
end of a struct — a **memory access violation**, not a wrong number, and only in
the one file that iterates the terms. Both arrays are `[ termCount ]` now: a
count that lives in the class has to be read from the class. That file also
asserts `η₆ == 0.0` there rather than a rate, which is the stronger statement —
a boundary functional has nothing to converge to on a fitted rectangle, and what
this file must guard is that the sixth term is **opt-in and inert**.

**AND `η` ROSE ON THE DRIVER'S FREE-BOUNDARY EXAMPLE — A ONE-LINE DEFECT IN THE
DRIVER.** The tempting suspect is the tabulated source's `j = 1` kink in
`∂F/∂ψ` being chased by `η₁`, and **there is no kink**:
`examples/fb-pprime.dat` is `f(Ψ) = 0.6Ψ` with `f′` **constant everywhere** over
the tabulated range, so `F` is linear in `ψ` across the whole domain. Guessing a
mechanism from the shape of the problem is what that is.

**THE CAUSE IS THE DATUM `η₅` IS COMPARED AGAINST.**
`GradShafranovSolver::transferredDatum()`'s default `g` is the **zero function**
— correct for every fixed-boundary case, and wrong the moment
`[boundary.exterior]` is present, where `Γ` carries the Gegenbauer trace
`Σ a_n C_n` that `setExteriorDatum()` deposits as a load. The driver never passed
it. The **library** test gets this right and says so in a comment; the driver was
not updated when the TOML wiring landed, and the only driver test of the coupling
was a **non-adaptive single solve**, so nothing could see it.

**IT DID NOT BIAS `η₅`, IT MADE IT DIVERGE.** `η₅²` carries an `h_e⁻¹` weight, so
an `O(1)` per-face mismatch contributes one copy of its square **per face** and
the term grows as `√(faces)`. Measured under near-uniform refinement:

→ **[M-06](MEASUREMENTS.md#m-06)** — `Γ_h` faces · 71 · 142 · 282 · 570

and the arc-length RMS of the omitted datum over `Γ` is **3.464e-02** against
`η₅/√faces` settling at **3.681e-02** — agreeing to **6%**. The accounting closes:
`η₅` *was* the omitted exterior datum, one copy per face.

**One variable changed, the datum:** `η₅` 2.864e-01 → **1.554e-04**, a factor of
**1844**, and `η` from **rising** 3.100e-01 → 5.519e-01 to **falling** 1.186e-01
→ 7.530e-02 over the same three cycles. Every other term fell throughout, which
is what localised it.

**AND IT WAS CORRUPTING THE MARKING, NOT ONLY THE NUMBER** — 33 and 68 elements
marked against 13 and 29 once repaired, the budget spent crowding `Γ_h` against
an indicator measuring nothing. That is the *refines the wrong elements* failure
`Estimator.hpp` warns about, in production. `TargetError` could never be met
either, since a diverging `η` never falls below anything.

`theDriverRefinesOverAnExteriorCoupling` is the regression, and it asserts
**monotonicity** — a real assertion rather than a formality, since the defect it
replaces produced a strictly increasing sequence. There had been **no driver test
of adaptivity over a coupling at all**, which is the gap that let a defect ship.

**`meq::AdaptiveDomain` HAD TO BE RELAXED AND THE OLD GUARD WAS A REAL
RESTRICTION.** It required `Ω` **strictly inside** the background box — exactly
one boundary attribute on the SubMesh — which every fixed-boundary case in this
tree satisfies and **the one geometry free boundary needs does not**: the
half-disc's flat side *is* the box's `r = 0` edge, because FB-A requires the
domain to reach the axis exactly. Inherited boundary is fitted and wants no
transfer, so the guard now checks that some boundary was **generated** — that
there is a `Γ_h` at all — and leaves inherited attributes out of
`gammaHMarker()`. Strictly more permissive, so no existing caller moves.

**`HighBetaConvergence` IS BIT-IDENTICAL**, every digit of the table under `CLAUDE_HDGGS.md`'s
*The measurement* including the `0.00e+00` and the `−5.55e-17`. The generalisation
keeps the **scalar** division when there is one border rather than routing it
through the dense solve, precisely so that it can be: the dense route computes
the same quotient by a different sequence of roundings, and a refactor that moves
the last bit of a published number is a refactor that has to be argued about.

**THREE DEFECTS ON THE WAY, AND THE THIRD IS THE ONE WORTH KEEPING.**

* **`prepare()` was not re-entrant with an exterior datum, and never had been.**
  Each call re-added a `VectorBoundaryFluxLFIntegrator` to `fluxRhs` while
  destroying the `PathTraceCoefficient` the previous one referenced, so the
  second `Assemble()` read through a dangling pointer. Latent for as long as
  nobody prepared twice. `fluxRhs` is a `unique_ptr` now and is rebuilt whole.
* **`DarcyNPCOperator` caches what it finds in the hybridization**, so a
  `prepare()` inside the solve leaves it stale and the process dies inside
  `solveWithNormalisation()` with no MEQ frame in the trace. It is rebuilt after
  every re-preparation.
* **THE DATUM IS A LOAD, AND A LOAD IS ASSEMBLED IN `prepare()`.** A Newton step
  that moved `a` and did not re-assemble evaluated its next residual against the
  PREVIOUS step's datum. It does not diverge, and that is what makes it worth
  recording: the transmission constraints sat at **1e-17** — the border was
  perfect — while `‖R‖` fell by a factor of **2/3 per iteration** on a problem
  that is affine. `a` and `x` chasing each other looks exactly like a
  well-conditioned solve converging slowly. **A linear rate on an affine problem
  is a Jacobian statement**, and `CLAUDE_HDGGS.md` already says so under
  *A wrong Jacobian is invisible to a convergence table*; the new part is that a
  perfectly satisfied constraint is not evidence the coupling is right.

**WHAT IT COSTS**, and it is more than §4.4 predicted: `N + 2` backsolves against
one factorisation, **plus one full re-assembly per accepted step**, because the
datum reaches the system through the right-hand side.

**THE `N + 2` ARE ONE CALL AND NOT `N + 2` CALLS.** They share a Jacobian by
construction, which is exactly `mfem::DarcyNPCSolver::ArrayMult()`'s case: one
pass over the mesh in `NPCReduce()`, one blocked trace solve, one pass back in
`NPCRecover()`. What is saved is the **traversal** — `GetElementFaces()`,
`GetFaceElements()`, `GetCtFaceMatrix()`, `GetFaceVDofs()` and the gathers run
once instead of once per column — plus, since PARDISO overrides `ArrayMult` and
is the default trace solver, one walk of the factors instead of `N + 2`. The
arithmetic on any one column is unchanged, so this is a traversal saving and not
a different method, and the agreement it owes is to round-off rather than
bitwise. **1.34× on the DIII-D solve leg at 14 columns**, and the decorator trap
that would have thrown half of it away silently:
→ **[M-98](MEASUREMENTS.md#m-98)**. That is the same price the
condensation path already pays for its own reasons.

**AND WHAT IS LEFT OF THAT TRAVERSAL IS THE BORDERED STEP'S SERIAL SPINE.**
Blocking it made it one pass instead of `N + 2`; it did not make that pass
parallel, and `NPCReduce()` / `NPCRecover()` are the only element loops in
`DarcyHybridization` without an `omp parallel`. What is left is
`O( elements × columns )` against integrator loops that are `O( elements )` and
already threaded, so on DIII-D at 10 modes it reads **0.212 s at one thread and
0.212 s at eight — 30.6% of a threaded step and the largest leg in it**. It has
a `StepProfile` leg of its own, and `other` fell from 45% to 14% when it got
one. `CLAUDE_HDGGS.md`, *Threading, measured*, and
→ **[M-126](MEASUREMENTS.md#m-126)**. Assembling the load directly
would remove it and would have to be checked against the differenced column
first, since `FormLinearSystem` transforms the right-hand side the residual is
measured against.

**AND THE ROW SIGN WENT WRONG ONCE, AS PREDICTED.** `exteriorTransmissionRows()`
is built to be contracted against `flux()`, which undoes `DarcyForm`'s `−q`;
the unknown carries the raw block. So contracting the row against the unknown
directly is right and negating it again is not. The wrong sign does not
diverge — it fails to converge, which is the same disguise as the stale load.

**THE COILS ARE WIRED, AND `[[coils]]` REACHES THE SOLVE.**
`meq::makeCoilSet` is the one door from the schema to a `CoilSet`, and
`meq::CoilAugmentedSource` / `CoilAugmentedNormalisedSource` are the adapters
that put `F_coil( r, z )` beside `F_plasma( r, z, ψ )` — both MFEM-free, both in
`Coils.hpp`, because it is the coils that need adapting: `Source::f()` takes ψ
and a coil current does not depend on it. `examples/coils-rectangle.toml` is the
worked example.

**THE CONTROL IS BIT IDENTITY AND IT IS WHAT MAKES THE WRAPPER SAFE.** Carrying
coils means every run with a `[[coils]]` block gets a *different object* handed
to `setSource()`, so the question is not only whether the coil term arrives but
whether wrapping perturbs anything else. `theDriverAddsTheCoilsToF` runs the
shipped example three ways and reads **`0.000e+00`** between zero-current coils
and no coils at all — not "agrees to round-off" — against **2.232e-01** relative
in L2 for 300 kA. An empty contribution is exactly empty.

**TWO THINGS THE WIRING CHANGED THAT ARE NOT ABOUT COILS.**
`NormalisedSource::setPlasmaSupport` is **virtual** now. It is consulted by
whichever object evaluates the profiles, and under a wrapper that is the *held*
source — so non-virtual, a caller holding a `NormalisedSource &` (which is what
the driver holds) would set the wrapper's flag, the wrapper's `f()` would
delegate to a plasma source still unconfined, and the moving support would
silently do nothing.
`the_plasma_support_reaches_the_source_that_evaluates_the_profiles` makes that
call **through a base reference** for exactly that reason. And
`SourceConfig::permeability()` exists so the coils can share the source's `μ₀`:
there is deliberately **no `Mu0` key on a coil block**, because a normalised-unit
run setting `[source] Mu0 = 1` while the coils kept the SI value would sum two
terms scaled a million-fold apart and converge, at full order, to a machine
nobody described.

**THE PLASMA SUPPORT IS A CONNECTED SET, AND THIS SECTION USED
TO CALL THE DEFECT LATENT.** `NormalisedSource::insidePlasma()` is a pointwise
test on the value, `(ψ − ψ_bnd)·span > 0`, with no connectivity. Across an
X-point the level `ψ = ψ_X` cuts a neighbourhood into four sectors and **two
opposite ones carry `Ψ > 0`** — the plasma, and the **private flux region** under
the divertor — so a diverted run describes a machine with a second current
channel nobody asked for. This file said it "cannot fire today: no shipped
example sets `ConfineToPlasma` and MEQ has no diverted case".

**IT IS LIVE ON A LIMITER CASE.** `theTwoBordersConvergeTogether`'s configuration
at limiter `R = 1.20`, reproduced to every published digit, has `{ψ > ψ_bnd}` in
**705 of 1333 elements and more than one piece** — midplane
`[0.35, 0.92] ∪ [1.25, 1.45]`; other rows of the §7.18 sweep read four and five
components. **No X-point is needed to make a level set disconnected**, and that
sentence had been generalising from the dramatic case to the only case.

`meq::PlasmaComponent` is the fill — MFEM-free, a CSR graph of ints in and a mask
out, so CI gates it — driven by `GradShafranovSolver::refreshPlasmaComponent()`,
which takes its adjacency from `Mesh::ElementToElementTable()` and refreshes
**before every residual and every Jacobian**, so both are taken at the same
support. `[source] PlasmaConnectivity` selects it, `"component"` is the default,
`"pointwise"` is the control, and it is **refused** under `CondenseThenLinearise`
rather than silently downgraded.

**AND §10.3's FACE-NEIGHBOUR PREDICTION IS HALF FALSE, WHICH IS THE FINDING.**
The half that holds is the one worth having: vertex-touching lobes really are
separated, so the element graph carries what `freegs4e`'s grid destroys and no
explicit X-point blocking is needed for *that*. The half that fails is the leak
estimate. **The lobes are not joined at the saddle — they are joined through the
BAND of elements straddling the separatrix**, every one of which carries `Ψ > 0`
at some vertex, so an inclusive candidate rule connects them all along the
divertor legs rather than at a point. Measured on `iterExample2` at 36,864
elements, a one-rule fill leaves **2,275 elements below the X-point, exactly what
the pointwise test leaves**. And §10.3's own cure, blocking the saddle's element,
is **resolution-dependent**: 161 of 2304 still leak at the coarsest of three
meshes.

**What ships needs no X-point finder and no parameter**: the fill traverses only
strictly interior elements and shares the straddling band out by a **watershed**,
reaching the blocked fill's mask to every printed digit. Fixed rings were tried
first and measured out.

**AND IT IS NOT A JUMP, WHICH DECIDES WHERE IT LIVES.** A lobe leaving whole is
`O(1)` however smooth the profile is, so this looked like FB-4's `j = 0`
discontinuity. Measured FB-4's way — slide `ψ_bnd`, refine the sampling, watch
the largest step — both the pointwise and the connected support fall like `h²`
(3.99, 4.00 against **3.94, 3.96**), so **the fill stays inside the Newton loop**.
The watershed never hands over a lobe, only a *straddling* element where `Ψ ≈ 0`
and at `j ≥ 1` the profile with it. **The ring-depth alternative DID jump** on the
same experiment — 1.68 then 1.20, floored — which is what says the property
belongs to the rule. It also helps: **24 Newton steps against 47** on a confined
rectangle, because at intermediate iterates the level really does fragment.

**AND `setPlasmaSupport()` ONLY REACHES THE FILL IF IT IS CALLED ON THE WRAPPER,
WHICH COST XP-2 THREE SESSIONS.** `plasmaComponentWanted()` asks
`normalisedSource->plasmaSupport()`, and `normalisedSource` is whatever the
solver was handed — on a free-boundary machine that is
`meq::CoilAugmentedNormalisedSource`, whose constructor does **not** adopt the
flag from the source it wraps. So `plasma->setPlasmaSupport( true )` followed by
wrapping leaves the wrapper saying **false**, `plasmaComponentWanted()` is false,
and **the fill never runs** — while `F` stays confined, because the inner
source's pointwise test is live. Nothing fails loudly and nothing is logged; the
only tell is `plasmaComponentElements()` reading `0 of 0 over 0 components`,
which is what *no fill is live* prints and not what *an empty plasma* prints.

`apps/meq.cpp` has always called it on the wrapper and is correct;
`XPointOuter`'s fixture did not, so **the one diverted case in the tree ran
without XP-1** — 333 elements of private flux region, 44% of the candidates,
carrying a current channel nobody asked for. → **[M-82](MEASUREMENTS.md#m-82)**.
The class documents the trap in the OTHER direction, a caller setting only the
wrapper's flag and leaving the wrapped source unconfined; this is its mirror
image and the documentation did not cover it. **The transferable part**: a
virtual forwarder runs one way, and the object a *third* party interrogates may
be neither the one you set nor the one that acts.

**AND THE MOVING SUPPORT IS REACHABLE TOO, AS `[source] ConfineToPlasma`** —
`F = 0` wherever `Ψ ≤ 0`, refused unless `Normalised = true` because the test is
on `Ψ`. **The coil term is OUTSIDE that support and the order the sum is taken
in is what puts it there**: a coil sits in the vacuum region by construction, so
confining it to the plasma would switch off every coil in the machine. **No
shipped example sets it**, and the reason is the precondition rather than the
plumbing: `p′( 0 ) = 0` is required or the assembled residual is discontinuous
in the unknowns, and `examples/mhd-ggprime.dat` is 1.35 at `ψ = 0`. With
`ψ_bnd` still fixed at zero (FB-3) the edge is pinned at `ψ = 0` rather than
found, so this is half of a moving boundary and the half that is built.

**A COIL THE MESH DOES NOT REACH CONTRIBUTES NOTHING, AND THE DRIVER SAYS SO.
FB-7 IS THE ROUTE THAT WOULD FIX IT, AND IT IS WRITTEN UP NOW** —
`FREE-BOUNDARY-PLAN.md` §7.19. The exterior stays linear, so
`ψ = ψ_coil + ψ̃` with `Δ*ψ_coil = 0` inside `Ω`: the interior equation is
untouched and the conductor enters as a **known** additive term on both halves of
the transmission condition, with no new unknowns. The one gap is a **gradient**
on `meq::CoilSet`, which has `psi()` and not `∇ψ`.

**IT BUYS DOMAIN REDUCTION AND NOT ACCURACY, AND THE OPPOSITE WAS NEARLY
RECORDED HERE.** A closed-form conductor field is not better than one assembled
on the mesh: a real coil has finite extent and finite current density, and near
the plasma that matters to the solution, so exactness of *evaluation* is not
adequacy of the *model*. And MEQ's coils were never filaments in the first place
— `meq::Coil` is a rectangular cross-section carrying a uniform current density,
and `coilPsi()` integrates the Green's function over that footprint at machine
precision outside it. The exterior route inherits that fidelity; what it removes
is the need to mesh out to the conductor.

`F` is assembled by quadrature over the elements, so a coil outside the `[mesh]`
box is never sampled — the run converges, writes its files, and describes a
machine with that conductor switched off, while `coil_current` in the `.nc`
reports the current whether or not it did any work. A **warning** rather than a
refusal, because `FREE-BOUNDARY-PLAN.md` §5.4 offers exactly that configuration
— coils outside `Ω`, entering through the exterior coupling — and that route is
now wired, as FB-7. **The warning stays and is still right**: it fires on a coil
in the `[[coils]]` block, which is the SOURCE, and FB-7's conductors do not go
there. The two routes are different objects on purpose — see below.

## FB-7: a conductor outside `Γ`, and the type that carries it has no `f()`

`meq::ExteriorCoilSet` in
`src/meq/Coils.{hpp,cpp}` — MFEM-free, so CI gates it — and
`GradShafranovSolver::setExteriorConductors()`.

**THE DECOMPOSITION IS THE WHOLE OF IT.** Write `ψ = ψ_coil + ψ̃`. The conductor
is outside `Γ ⊇ ∂Ω`, so `Δ*ψ_coil = 0` **inside `Ω`** and the interior equation
is untouched — there is no coil term in the source at all; and `ψ̃` is
`Δ*`-harmonic in the whole exterior and decays, its singular support having been
subtracted, so the Gegenbauer expansion is valid for it. The conductor therefore
enters **only through `Γ`**, additively and known, on both halves:
`ψ|_Γ = Σ a_n C_n + ψ_coil|_Γ` and `q·ν|_Γ = (modal) + q_coil·ν|_Γ`. **No new
unknowns and no change to the border's size.**

**BOTH HALVES LAND TOGETHER OR NOT AT ALL**, which is why it is one call. Adding
`ψ_coil` to the datum while leaving the transmission row alone is an inconsistent
pair that **converges to something** — the disguise the stale-load defect wore,
where the border sat at 1e-17 while the answer was wrong.

**THE TYPE HAS NO `f()`, AND THAT ABSENCE IS THE DESIGN.** An exterior conductor
contributes nothing to the interior equation, so there is no interior current
density to return — not merely one nobody wrote down. And **a filament has none
at all**: infinite current density on a set of measure zero is not a function, so
`CoilSet::f()` could not be honest about one however the set were arranged. That
is what had blocked a filament from reaching the coupled path, and removing the
one method neither member kind can answer is what lets a rectangle and a filament
share a set. **Handing the solver a `meq::CoilSet` is now a compile error**,
where before it compiled and converged having silently dropped that conductor's
current from the source.

**THERE IS DELIBERATELY NO CONVERSION FROM `CoilSet`.** The one error the type
exists to prevent is a conductor appearing in **both** — once in the source and
once through the coupling, double-counted while every residual and every border
converges. Making that conversion one line would make the mistake one line.

**WHAT IS MEASURED**, `tests/convergence/FreeBoundaryCoupling.cpp`:

→ **[M-07](MEASUREMENTS.md#m-07)** — the datum given, `ψ` against `ψ_coil`

**AND THE LAST ROW IS THE INTERIOR ROUTE USED AS A CONTROL ON THE EXTERIOR ONE,
WHICH NOTHING HAD DONE.** A conductor is at a fixed position, so it is inside `Ω`
or outside `Γ` and never both — the comparison needs a **second geometry**, and
it is a **fitted rectangle** `[ 0, 2.6 ] × [ −1.8, 1.8 ]` containing both the
conductor and the whole half-disc. There the conductor is a **source term**
assembled by quadrature (`meq::CoilAugmentedSource`, the `[[coils]]` route) with
`ψ_coil` as the datum; on the half-disc it is outside `Γ` **and outside the
background box**, reaching the solve through the live coupling. They are compared
at 240 fixed points inside `ρ = 1.2`. **A rectangle on purpose**: no SubMesh, no
level set, no transfer path, no extension, no exterior expansion, no border — so
the two runs share the equation, the element loops and essentially nothing else,
and run B carries no modal truncation to floor the comparison with. The control
is run B with its coil term removed and the datum kept, which converges perfectly
well to a different function.

**THE MESH IS ALIGNED TO THE CONDUCTOR AND THE CASE ASSERTS IT PER MESH**, §7.9
having measured a cut conductor costing more than a degree. **The prediction
before running it was `min( k+1, 3 )` from the `r² log r` at the corners of a
top-hat source, and it is half right**: run B's L2 over its **whole box** reads
**2.741 / 3.052** at `k = 2,3` — §7.9's aligned 2.88 / 3.01 reproduced on a
different mesh family — while the **overlap** column reads 3.879 at `k = 3`,
because the corner term is local and the comparison region is half a metre away.
Both are right and both are printed. **And the interior route is the BETTER of
the two there**, three to five times closer to `ψ_coil` at `k ≥ 2`, so the
difference tracks the *exterior* route and the extension path's own `O( h )`
geometry cost on `Γ_h`.

**THE GEOMETRY IS WRITTEN OUT TWICE, WHICH IS THE PRICE OF THE TYPE SPLIT**, and
the case asserts `CoilSet::psi` against `ExteriorCoilSet::psi` at twenty points
before comparing solves: a typo in one copy is how this would become a comparison
of two different machines while still converging. **34 s wall.**

**THAT FILAMENT ROW IS THE CROSS-CHECK AND IT IS WHAT THE FILAMENT BOUGHT.** The
datum-given pair differs only in Dirichlet data; the coupled pair differs in that
**and** in `q_coil·ν` through the transmission row, and solves for `a` besides.
Same physical quantity two ways, and **a Neumann half inconsistent with its own
Dirichlet twin would separate them while every border still converged.**

**WHY A FILAMENT AT ALL, AND IT IS NOT THAT IT IS BETTER.** `../freegs4e`'s
default `Coil` **is** an exact filament — `controlPsi` returns
`Greens( R, Z, · )*turns`, its `area` only imposing a current-density limit and
never entering the field — and so is FreeGS's. `meq::Coil` is a rectangle with
uniform current density. §7.16's 1.3e-04 agreement was reached **across** that
mismatch rather than with it removed, and `meq::CurrentFilament` is what lets the
two be matched. A real conductor has finite extent and near the plasma that
matters; the filament is for measuring the difference, not for being right.

**AND THE PRECONDITION IS ENFORCED NOW, IN BOTH ORDERS.**
`ExteriorCoilSet::clearance( centreZ, ρ_Γ )` measures the distance from `Γ`'s
centre to the **nearest point** of each member — the closest edge of a rectangle,
the ring itself for a filament — so a coil whose *centre* clears `Γ` while its
inboard edge does not is caught, which a `hypot()` at the call site is not. Both
`setExteriorConductors()` and `setExteriorCoupling()` refuse on it, whichever
arrives second, which closes the *"IT MUST BE OUTSIDE `Γ`, AND THAT IS NOT
CHECKED HERE"* the solver's own header carried. A conductor inside `Γ` is not a
small error: `ψ_coil` is then not `Δ*`-harmonic where the expansion assumes it
is, and the run converges to a machine nobody described.

**AND IT PAID FOR ITSELF THE DAY AFTER IT LANDED** — §7.12b's fixture had
*"nowhere to put a coil genuinely outside the plasma"*, and a conductor outside
`Γ` costs no mesh. See *§7.12b's fixture was the defect*.

**FB-1 IS COMPLETE** — see the entries below for the two halves
and for what building it found. What remains of free boundary is the **plasma**:
FB-2's prescribed current on a solve, FB-3's `ψ_bnd`, FB-4's moving support and
cut quadrature, and FB-5's bordered Newton, which is where the superposition
FB-1b uses stops being exact.

## FB-6 is beaten, and what was in the way was the limiter, pinned to a dof

`LimiterConstraint::ExactPoint` is the default and
`NearestDof` — what MEQ did — is the control. `FREE-BOUNDARY-PLAN.md` §7.20 is
the record; the outcome is in *MEQ against freegs4e* below.

**FB-3 PINNED `ψ_bnd` AT THE NEAREST POTENTIAL DOF TO THE LIMITER CONTACT, AND
THE HEADER CALLED THAT "A DEFINITION RATHER THAN AN APPROXIMATION … THE SAME
CHOICE `ψ_ax` MAKES".** The analogy is false and it is what hid the size of it.
A nodal maximum is wrong by `O( h^{k+1} )`, because a polynomial's peak over a
closed element is within one interpolation error of its largest nodal value. A
limiter contact is a **prescribed point**, the nearest dof is up to half a dof
spacing from it, and `ψ` there is out by `dist × |∇ψ|` — **`O( h )` at every
degree**, so `p`-refinement makes it worse per dof rather than better.

**IT IS A STAIRCASE IN THE POINT ASKED FOR**, which is the unambiguous form of
it. Sweeping `[boundary.limiter] R` on `examples/limited-tokamak.toml`, one key
changed, the **whole solve is bit-identical over `R ∈ [ 1.3250, 1.3500 ]`** — a
plateau **0.025 m** wide, 7% of the minor radius, same 11 Newton steps, same
`ψ_ax = 9.466087e-02` to the last bit — and steps by **5%** in `ψ_ax` and
**11%** in `ψ_bnd` at each end. Repaired, the same sweep is smooth:
`ψ_bnd` 2.685e-02 → 2.847e-02 monotonically over the same range.

**AND IT IS WHY THE SHIPPED FIXTURE APPEARED TO CONVERGE, WHICH IS THE PART
WORTH KEEPING.** Uniform refinement of the shipped mesh, `ψ_ax` at
1601 / 6527 / 26375 elements:

→ **[M-08](MEASUREMENTS.md#m-08)** — order · limit · vs the 129² reference

**The control converges beautifully, at order 2.89, to the wrong number**,
because this fixture's limiter happens to sit within 1e-4 of a dof, so snapping
is nearly exact *at that point on that mesh* and refinement never disturbs it.
Move the contact 0.6 m round the same limiter circle — to where the 513²
reference puts it — and the same ladder **scatters by 1.7%** instead of
converging, while `ExactPoint` converges to 3.2e-05. **A clean convergence table
is not evidence that the boundary condition is right**, which is this tree's
oldest lesson — `CLAUDE_HDGGS.md`'s *A wrong Jacobian is invisible to a
convergence table* — arriving on a boundary condition rather than on a
Jacobian.

**THE ROW IS SIMPLER THAN THE AXIS ROW IT COPIES.** `ψ_bnd = ψ_h( r, z )` inside
the containing element, the row that element's potential shape functions there,
the corner still exactly 1. `AxisConstraint::LocatedAxis` needs the **envelope
theorem** to be exact because its point moves with the solution; a prescribed
limiter contact does not move, so there is no position term to argue away.
`theLimiterConstraintIsEvaluatedWhereItIsAsked` sweeps both choices and gates on
the control **repeating a value over four consecutive samples** while the repair
repeats none, with the two **6.3%** apart in `ψ_bnd`.

**AND `TransformBack` IS SAFE HERE WHERE THE AXIS FORBIDS IT.** A zero of a
discontinuous `q_h` can lie outside its own element and the inverse map does not
converge there — measured, a residual of 5.1e-02 on an element of 5e-02. A
limiter contact is a point of `Ω` somebody asked for, so it is inside an element
by construction, and on straight-sided triangles the inverse is affine and
`Inside` is exact rather than probable.

**THE LIMITER IS A CURVE, IN THE RESTRICTED FORM WHERE IT IS A POLYGON THE MESH
IS FITTED TO.**
`setLimiterSurface( attribute )` and `LimiterConstraint::LocatedContact`: the
limiter is the boundary of an element-attribute region, so it **is** a union of
mesh faces, `ψ_bnd = max ψ_h` over it, and the contact is an output of the solve
rather than an input to it. `tests/convergence/LimiterCurve.cpp` is the
acceptance at **20.2 s** under `ctest -j4` on a quiet machine — and
31.8 and 29.5 s standalone while two agents were building, three readings of one
binary, which is `CLAUDE.md`'s own standing caution about reading a suite time
as a measurement; `FREE-BOUNDARY-PLAN.md` §7.20 has the full table.

**THE ROW IS AGAIN THE SHAPE FUNCTIONS, AND NEWTON'S ORDER IS THE ONLY THING
THAT CAN CHECK IT.** The envelope argument now has to cover two cases and does:
on the interior of an edge the tangential derivative vanishes at a maximum of
the restriction while the contact moves *tangentially*, and at a vertex the
contact does not move at all — so no sensitivity of the search enters the
Jacobian either way. Nothing in an error norm can see whether that is right —
`CLAUDE_HDGGS.md`, *A wrong Jacobian is invisible to a convergence table* — so the assertion is the
**observed order, and it reads 2.000 at `n = 24` and `n = 48`**.

**THE FIXTURE KNOWS ITS OWN ANSWER, WHICH IS WHAT MAKES THIS SHARP RATHER THAN
SELF-CONSISTENT.** Domain, source and boundary condition are all symmetric in
`z`, so the exact contact on the outboard edge is `z = 0` exactly, while
`MakeCartesian2D`'s diagonal split is **not** symmetric — so the discrete
contact converges to the midplane rather than sitting on it, at **2.217 and
2.127**. That is `O( h² )`, which is what a root of a *differentiated* quantity
gives; the axis gets `k+1` from rooting the **solved** `q` instead, and the same
trade is available here and not taken.

**AND THE CONTROL IS A WRONG POINT ON THE SAME POLYGON, WHICH MUST NOT
CONVERGE**: 7.072e-02 → 7.051e-02 → **7.049e-02**, flat, against the located
contact's 2.149e-04 → **2.382e-06** — **29,600× apart** at `n = 48`. A contact
prescribed in the wrong place is wrong by `dist × |∇ψ|` however fine the mesh,
which is this section's own defect in miniature.

**PAINTING A REGION IS NOT FITTING A MESH TO A CURVE, AND THE DIFFERENCE IS AN
ORDER.** Painting by centroid gives a closed polygon whose perimeter converges
to the **wrong number** — 2.828, 2.368, 2.611 against a circle's 1.885 — the
staircase, so `max ψ` over it converges at `O( h )`. `halfdisc.py --limiter`
fragments the circle into the geometry instead, which puts the polygon's
**vertices on the true circle to 3.3e-16** and inscribes at `O( h² )`. Verified
through the solver on a gmsh mesh: 19 faces recovered from the attribute
interface, 7 Newton steps, contact **0.345370** from the limiter centre against
a radius of 0.350 — the chord sag of a 19-segment inscribed polygon,
`a( 1 − cos( π/19 ) ) = 0.0048`, to the digit.

**WHAT IS STILL RESTRICTED** is that the limiter is a polygon: a smooth limiter
is reached only through the mesh, at whatever inscription order the mesher
gives. No element is cut and none needs to be — the same trade §7.9 made for the
conductors, and for the same reason, a limiter being prescribed input that does
not move.

**AND IT REACHES THE SOLVE FROM A FILE**, as `[boundary.limiter]
SurfaceAttribute`, which is an **alternative** to `R`/`Z` and is refused beside
them — both constraints converge, to equilibria differing by the `O( h )` the
point costs, so a precedence rule would decide which equilibrium is reported on
the strength of key order. Driven end to end on a `halfdisc.py --limiter` mesh
the driver reproduces the library to every digit: `ψ_ax = 7.459870e-01`,
`ψ_bnd = 6.698128e-01`, contact `( 1.038429, 0.343225 )`, 7 Newton steps, tail
order 1.999. The `.nc` carries `limiter_r`, `limiter_z` and
`limiter_contact_located` — **the contact that produced `ψ_bnd`, not the
configuration echoed back**, which on this route are different things and only
one of them is any use to a reader differencing two runs.

**AND THE DEFECT THAT FOUND WAS DRIVER-SIDE, WHICH A LIBRARY TEST STRUCTURALLY
COULD NOT SEE.** `buildSubdomain()` selects `Ω_h` by overwriting **every element
attribute** with 1 or 2 and cutting a `SubMesh`, so a `.msh`'s own regions —
`halfdisc.py` writes conductors as `10 + i` and the limiter's interior as `20` —
were gone before `SubMesh` copied anything, and `SurfaceAttribute` refused a
mesh that plainly carried attribute 20. `LimiterCurve` builds a plain mesh and
paints its own regions, so no cut stands between the attribute and the solver.
**Between a file mesh and the solver there is always a cut, and only a driver
test has one.** The cut's marking is scratch — `meq::AdaptiveDomain` uses 1/2/3
for the same purpose and says so — and the material attributes are not; both
meshes are restored across it. `theDriverFindsTheLimiterContact` is the
regression, on `examples/limiter-halfdisc.toml`, and it asserts the contact
lands **on the polygon**: 0.345370 from the limiter centre in the 4.8e-03 band
between the polygon's inradius and its circle. **Its control is a fixed point**
— prescribing the contact the search found must reproduce the solve, since the
maximum is attained there — and the two agree to **5.5e-12**.

**ON THE MACHINE CASE, FINDING THE CONTACT IS 11.7× CLOSER TO THE CONVERGED
REFERENCE, AND IT COSTS BRANCH SENSITIVITY.** On a rebuild of
`limited-tokamak`'s geometry with the limiter meshed in — 1607 elements against
the shipped 1601 — the prescribed 129² contact gives `ψ_ax = 9.490127e-02`,
1.95e-02 from the 513² reference, and the found one 9.293175e-02, **1.67e-03**.
The found contact lands on the **inboard shoulder** at `( 0.7739, 0.2571 )`,
which is where that file's own comment says the true circle's maximum is on both
grids. **The shipped example keeps its prescribed point and should**: it
reproduces `freegs4e`'s own 129² artefact, which is what makes the comparison
against that grid well posed. What is new is the measurement of what the
artefact costs. **From the shipped guess the found-contact run reaches a
different equilibrium** — 38 Newton steps, scale 13.6, `ψ_ax` at **1.48×** the
peak of the field it wrote, which the existing gate catches — and from a
converged nearby equilibrium it takes 8. **Refinement does not cure it**: the
same cold run at 5225 elements does not converge at all, 200 iterations with the
residual floored at 6.84e-04 and the driver exiting 2, where the prescribed-point
run on that mesh takes 7. The branch is chosen by the guess, not by the mesh.
Prescribing a point at the found contact reaches a *third*, whose own maximum
over the polygon agrees with what it was given to 0.2%, so these are genuinely
distinct fixed points and not a broken search. `FREE-BOUNDARY-PLAN.md` §7.20 has both tables.

**A WRONG SURFACE ATTRIBUTE IS REFUSED AT THE SETTER, NOT AT THE SOLVE**, which
matters because the failure is otherwise silent: `max` over an empty polygon is
minus infinity and its border row is all zeroes, so the bordered solve does not
diverge — it converges, to a `ψ_bnd` pinned by nothing. The commonest way to
reach it is a mesh built without `--limiter`. The polygon is collected again per
solve regardless, because an adaptive cycle refines the mesh under a solver that
keeps its settings and stale face indices do not throw, they name other faces.

**AND THE SMOKE TEST THAT FAILED FIRST WAS THE FIXTURE, WHICH IS THIS FILE'S
OWN RECURRING LESSON.** The first driver run creeped for 60 iterations without
converging, on a half-disc with a tabulated `p′ = 0.6Ψ`, no confinement, no
prescribed current and no coils. **The control settled it in one run**: the same
file with a *prescribed* contact fails at least as badly, and worse — so it is
§7.14's non-linear eigenvalue problem and §11.7's three missing ingredients, not
the constraint. One variable changed.

**AND IT CLOSED THE SPIKE BRANCH ON THE MACHINE CASE, NOT MERELY DETECTED IT.**
`examples/limited-tokamak.toml` carried a comment saying the degree was load
bearing — at `PolynomialDegree = 2` the run reported `ψ_ax = 2.73e+00` against a
reference 9.48e-02, a single spiking nodal value, so *"`p`-refinement is what
reaches the physical branch"*. Re-measured at `k = 2` on 1601 elements with
everything else held: **8 Newton steps to `ψ_ax = 9.349647e-02`**, the physical
branch, 1.4% from the reference, which is about what a degree costs. **A single
spiking dof is not a zero of `q_h`.** So the located-axis constraint makes that
branch unreachable rather than merely reportable, which is more than the
detection §11 set out to build.

**AND `AxisConstraint::LocatedAxis` CAN LOCK ONTO A CONDUCTOR'S OWN O-POINT**,
found while doing this. On the coarse mesh with the 513² currents the solve
**converges** — 199 Newton steps, every border satisfied, the current delivered —
to `ψ_ax = −1.232718e-01` at **( 1.7516, 0.9000 )**, the P1U coil centre, with a
profile scale of **−1.03**. The coil current is negative so its own field has a
minimum there, the span came out negative, and the constraint went looking for a
minimum and found the coil's. **It reports `normalised flux 1.0000`**, because
`Ψ` at the located axis is 1 by construction under that constraint — the
tautology recorded under *IN-A* meeting a case it cannot see. Refinement cures
it; the tells are the axis POSITION and the negative scale.

**AND THE DRIVER WARNS ON IT NOW.** `apps/meq.cpp` tests the located axis
against every `[[coils]]` conductor's own rectangle and says so, because that is
the one thing here that is **not circular**: every test phrased in terms of the
plasma is satisfied at the coil, `Ψ` being 1 there by construction. A warning
rather than a refusal, on the same precedent as the several-O-points one — which
branch is the core is the guess's choice, and refinement cures this one. Verified
firing on the case above: *"the located magnetic axis ( 1.7516, 0.9000 ) is
INSIDE conductor 0, which spans r [ 1.7000, 1.8000 ] z [ 0.8500, 0.9500 ]"*.

## The plasma edge caps the order, and it is the PROFILE that sets the cap

**FB-4's question, answered, and the answer moved the work rather
than doing it.** `tests/convergence/PlasmaEdgeConvergence.cpp` and
`tests/analytic/PlasmaEdge.hpp`.

The question was *can `ψ*` keep `k+2` across a plasma edge that cuts through
elements and moves while Newton runs*. `FREE-BOUNDARY-PLAN.md` §5.3 framed it as
a quadrature problem, named `MomentFittingIntRules` as the tool and the
**derivative** of a cut rule as FB-4's one real gap, and quoted CEDRES++ saying
this is where they stopped going above first order.

**IT IS NOT A QUADRATURE PROBLEM.** Let `j` be the order to which the profiles
vanish at the edge — `p' ~ Ψ^j`, which is a **modelling** choice a user makes
and not a numerical one; FreeGS's `(1 − Ψ_n^α)^β` gives `j = β` and defaults to
`β = 1`. Then the exact `ψ` carries `|d|^{j+2}` across the edge, and:

→ **[M-09](MEASUREMENTS.md#m-09)** — best approximation, ANY method · MEQ, plain Gauss rule

**So `ψ*` keeps `k+2` exactly when `k ≤ j`**, and that threshold is the same for
an exact cut rule as for a blind one, because `k+2 ≤ j+2.5` and `k+2 ≤ j+2`
differ only off the integers. Measured on four rungs, with MEQ's ordinary
quadrature and no cut rule anywhere:

→ **[M-10](MEASUREMENTS.md#m-10)** — `j` · `k = 1` · `k = 2` · `k = 3` · `k = 4`

The bold entries are `k+2` reached; every one of them has `k ≤ j` and no entry
with `k > j` reaches it. **`j = 3, k = 3` reads 4.989 against a target of 5**,
which is the sharpest single statement here: `k+2` survives a plasma edge
outright, at third order, with nothing built.

**THE UPPER BOUND IS MEASURED WITHOUT A SOLVER, WHICH IS WHY IT IS DECISIVE.**
`theCutCapsTheOrderBeforeAnyMethodIsChosen` takes the L2 **best approximation**
of the exact solution by `P_k` and `P_{k+1}` — the two spaces `ψ_h` and `ψ*`
live in — with a composite rule on the cut elements so that what is measured is
the approximation and not the quadrature. A discrete solution cannot beat its
own space, so `min(k+2, j+2.5)` bounds **any** method however the cut is
integrated. The mechanism is that the best polynomial approximation of `|d|^m`
on an element of size `h` is `O(h^m)` at every degree — the *constant* falls
with `k`, the order does not — over `O(1/h)` elements of area `O(h²)`.

**AND AWAY FROM THE BAND `k+2` IS ALWAYS THERE**: the uncut elements read 2.96,
3.95, 4.99 at `j = 0`, where the total is 1.9. The loss is a set of measure
`O(h)` and nothing else, which is worth knowing for a consumer who cares about
the core.

**THE DISCRIMINATOR THAT SAYS THE QUADRATURE IS BLIND AND NOT COARSE**: sweep
the rule at fixed geometry. Over `extra = 4, 8, 12, 16, 20` at `j = 1, k = 2`,
`ψ_h`'s rate reads **2.582, 2.577, 2.578, 2.577, 2.576** and `q_h`'s 2.400,
2.404, 2.405, 2.404, 2.404 — pinned to three figures. A Gauss rule cannot see a
kink between its points however many it has. `ψ*` **does** move, 2.87 → 3.57
against its bound of 3.5, so the half order the shipped default costs is
recoverable with the rule MEQ already has. `setSourceQuadratureOrder()` is that
knob and it exists for this measurement.

**SO NO CUT QUADRATURE WAS BUILT *FOR THE ORDER*, AND THAT IS A DECISION
RATHER THAN AN OMISSION.** An exact cut rule would move `ψ_h` from `j+1.5` to
`j+2.5` and `ψ*` by nothing that a higher Gauss order does not already buy; it
would buy `q_h` **nothing at all**, `q` being at its own regularity bound
already; and it would not move the `k ≤ j` threshold. **What it WOULD buy is
`j = 0` at all — see the next section, so "a cut rule buys nothing" is too
strong.** `refs/CutElementQuadratureSurvey.pdf` (Loibl et
al., arXiv:2602.18130) is the survey of the ten open-source implementations and
it adds two practical reasons: every one of its benchmarks is Cartesian, and
**MFEM's two cut backends are quadrilateral and hexahedral only** — the Algoim
path aborts with *"supports only quads and hexes"* and `MomentFittingIntRules`
builds four-vertex local meshes — where MEQ's meshes are triangles. Turning
`MFEM_USE_ALGOIM` on would therefore cost a mesh-type change as well, for a
threshold that does not move.

**`MomentFittingIntRules` WAS TRIED AND IT WORKS, WITH A CAVEAT WORTH KEEPING.**
On quadrilaterals it is exact on a straight cut (9e-16) and reached 3.3e-06 on a
disc — but on one mesh in three it produced weights of **−30 and −63 on elements
of area 1e-4**, and the area came out 1.3e-02 wrong by cancellation. That is the
known conditioning fragility of moment fitting on a nearly degenerate cut, and
it is why the survey's other family — dimension reduction, which is Algoim — is
the one that reaches its design order.

**AND THE FIRST ATTEMPT AT IT READ ZERO FROM EVERY CUT ELEMENT**, because it
held `mesh.GetElementTransformation( int )` across the call and MFEM's own cut
code resets that shared scratch underneath it. Third time this trap has been
recorded in this tree — `CLAUDE.md`'s *Traps* is where it lives; the fix is a
local `IsoparametricTransformation`, which is what `ex38.cpp` does.

**MFEM'S TRIANGLE RULES GO BAD ABOVE ORDER 25 AND NOTHING WARNS YOU.** Measured
while sweeping `setSourceQuadratureOrder()`: exact and positive-weighted to
order 25, and from **26** a construction with a least weight of **−3.6e+01**,
reaching **−1.9e+07** by order 64, with the integral of a monomial degrading
from 1e-16 to 2.7e-05. A solve at `2k + 30` returns errors of order **1e+3**.
So that knob is only meaningful while `2k + extra ≤ 25`, which at `k = 3` is
`extra ≤ 19`. `theLossIsTheRulesBlindnessAndNotItsResolution` pins the boundary
so the next caller meets it as an assertion.

## An X-point on that edge costs the approximation NOTHING — measured

**§10.1 PREDICTED IT, NOTHING HAD MEASURED IT FOR THE WHOLE CAMPAIGN, AND IT
HOLDS** → **[M-133](MEASUREMENTS.md#m-133)**. The prediction cuts against
expectation and the section says so: at a null `Ψ` vanishes **quadratically**,
so a source `F ~ Ψ^j` vanishes to order `2j` there and the crossing is
*smoother* than the branches. The worry was never that — it was that the
plasma's **support** acquires a corner at the null, and a corner is measured to
cost MEQ's extension an order.

**IT DOES NOT COST THE APPROXIMATION ONE.** Best approximation in `P_k` over
`n = 8, 16, 32, 64`, a crossed edge (`φ = dz² − dr²`, two lines through a null)
against a smooth one (`φ = a² − dr² − dz²`), with the cut as the only thing that
differs: the **worst rate drop over nine `( j, k )` pairs is 0.055 of an
order**, against the 0.25 a moving cut is allowed, and the crossed arm is
*faster* in seven of the nine.

**SO THE CAP ABOVE IS THE DIVERTED CAP TOO.** `min( k+1, j+1.5 )` for `ψ_h` and
the `k ≤ j` threshold for `ψ*` are properties of the profile's exponent, and a
null on the edge does not add a term to either. That is the half of §10.1 that
had been resting on an argument; the solved half is XP-3's and is green.

**TWO LIMITS, BOTH IN THE MEASUREMENT RATHER THAN AROUND IT.** Only the *rates*
compare — the two arms are not normalised to a common magnitude and their
constants stand up to 2000× apart, which is why the assertion is a difference of
rates. And it is best approximation, so it bounds what any method on these
spaces can do and says nothing about whether MEQ's solve attains it on a
diverted machine, where the fill, the moving support and the X-point border all
still apply.

## At `j = 0` the question does not arise: Newton cannot solve it at all

**THE SHARPEST FB-4 FINDING, AND IT IS ABOUT THE JACOBIAN RATHER THAN THE
ORDER.** With the support read off the solution — `Ω_p = { ψ_h > 0 }`, nothing
telling the solver where the edge is — a source with `p'(0) ≠ 0` **does not
converge**: not at `k = 1, 2, 3`, not at `n = 8, 16, 32, 64`, **not started from
the exact solution**, and **not under `PicardThenNewton`**, which cures every
other hard case in this tree.

Starting from the exact solution and failing is what rules out the comfortable
explanation. This is not a basin.

**AND THE MECHANISM IS MEASURED, NOT ARGUED, BECAUSE THE OBVIOUS STORY IS TOO
WEAK.** The obvious story is that the residual is Lipschitz but not
differentiable. **It is worse than that: with a FIXED rule the residual is
DISCONTINUOUS in the unknowns.** A quadrature point sits at a fixed reference
position; as the coefficients move, `ψ_h` at that point crosses `ψ_bnd` and `F`
there jumps from zero to its edge value, so the element integral jumps by
`w_q·|jump|`. Measured on one cut element, sliding the level set across it and
sampling the integral at ever finer intervals — the largest step between
neighbouring samples, as the interval is quartered:

→ **[M-11](MEASUREMENTS.md#m-11)** — 401 samples · 801 · 1601

**The plain rule's step does not shrink**, which is what distinguishes a jump
from a steep slope; the composite's does. At 3.4% of the integral's own scale
that is not a small perturbation, and Newton is chasing a root of a
discontinuous function — which is why an exact starting point and
`PicardThenNewton` are equally useless. At `j = 1` and `j = 2` the two rules
agree to the last figure, because `F → 0` continuously as a point crosses.

**SO A CUT RULE *WOULD* FIX `j = 0`, AND NOT BY IMPROVING THE ORDER.** With the
integration domain following the level set, the integral is a continuous — and
differentiable — function of the coefficients, and its derivative carries
exactly the surface term `∮ F φ/|∇ψ|` that §5.3 predicted, which the cut
*surface* rule supplies alongside the volume one. **That is the one thing worth
building a cut rule for**, and the prize is solvability at about second order
rather than `k+2`, since the approximation cap of `j + 2.5 = 2.5` stands
regardless. Set against `j ≥ 1` being one line in a profile, it is not worth
it — but the reason is the cost-benefit and not that a cut rule is useless.

**AND THE ROUTE TO `k+2` AT `j = 0` IS DESIGNED AND NOT BUILT.**
`PLASMA-EDGE-PLAN.md` is the plasma edge as an interior interface coupled at a
distance — two HDG subdomains meeting across an `O(h)` band, with the interface
trace as a modal border, which is FB-1's structure with the exterior DtN
replaced by a second interior solve. It reuses stage 5, `TransferPath`,
`ElementExtension`, `ExtensionBoundaryQuadrature` and the bordered Newton, and
its PE-0 is a day's work against a fixture that already exists. **It is not to
be started until `j ≥ 1` is finished**, and it covers limiter plasmas only —
the transfer families give out at an X-point.

`∂F/∂ψ` also acquires the surface term `F·δ(Ψ)` that `meq::Source` structurally
cannot carry. §5.3 predicted the term and said to *"decide this deliberately and
write it down"*; the decision is that `j ≥ 1` is a **precondition** of MEQ's
free-boundary path rather than a convention, and it is now measured rather than
assumed.

At `j ≥ 1` the residual is `C¹` and Newton is ordinary — 3 to 4 iterations at
every mesh and degree, `2.71e-01 → 9.22e-04 → 8.88e-07 → 4.20e-12` at `j = 1`,
and `2.71e-01 → 3.20e-06 → 5.20e-13` at `j = 3`. The observed order is
**superlinear rather than quadratic** at `j = 1`, about 1.7–1.8, which is what a
`C¹`-but-not-`C²` residual gives: `∂F/∂ψ` itself jumps at the edge, so the
Jacobian is not Lipschitz.

**AND MEQ HAS NO CUT-RULE DERIVATIVE GAP, BECAUSE IT USES NO CUT RULE.** §5.3
and §6.4 both name the sensitivity of a cut rule to the level set as FB-4's one
real gap. With a fixed Gauss rule the quadrature points do not move, so the
assembled Jacobian is the **exact** derivative of the assembled residual
whatever the edge is doing. Adopting a cut rule is what would create the gap it
was meant to close.

**THE MOVING EDGE COSTS THE RATE NOTHING.** `MovingPlasmaEdge` against
`PlasmaEdge` at the same `j` and `k`, on the same meshes: 2.000/1.995/3.012
against 2.000/1.998/3.010 at `j = 1, k = 1`, and 2.435/2.461 against
2.431/2.455 at `k = 3`. The order is a property of where the edge **converges
to**, not of its having moved to get there. So the two fixtures are a matched
pair and the fixed one is where a rate study belongs.

**THE FIXTURE PAIR, AND WHY THERE ARE TWO.** `PlasmaEdge` puts a
`Δ*`-harmonic vacuum field outside a prescribed circular edge, so `F` is
supported in the plasma **and nowhere else** — the free-boundary structure — at
the price of a cut that does not move. `MovingPlasmaEdge` makes the plasma
exactly `{ψ > 0}` of its own solution, at the price of a smooth background
source outside it. **The price is not laziness and cannot be avoided**: `Δ*`
has no zeroth-order term, so it obeys a maximum principle and `{w > w₀}` can
never be compactly contained for a `Δ*`-harmonic `w` — a real equilibrium
confines its plasma with **coils**, and a manufactured fixed-boundary problem
has to put something there instead.

**The transmission row** `∫_Γ E_h(q_h)·ν C_m dΓ` is the Neumann half of the
coupling.
**Its MFEM half was written here and is now upstream's** —
`mfem::ExtensionBoundaryQuadrature`, merged into `gf-hdg-subdomains-dev`
 Writing its tiling check (the boundary weights must sum to `|Γ|`)
found that an **unsigned** weight integrates a folded sweep with multiplicity;
upstream then found the same defect in the sibling `ExtensionRegionQuadrature`
and signed it, worth fifty-fold on their aerofoil. **MEQ's own disc does not
fold and MEQ's numbers are untouched** — the return on filing went to somebody
else's test case, which is the argument for filing a thing that has been *used*
rather than merely asked for.

**THE TRANSMISSION ROW IS WRITTEN AND IS MEASURED AGAINST A CLOSED FORM.**
`exteriorTransmissionRows()` sweeps `Γ` with that routine, evaluates `E_h(q_h)`
at each foot through `mfem::ElementExtension`, and contracts against each mode.
**The integral carries no `1/r` and that is not an omission**: the exterior block
is diagonal in the weight `dΓ/r`, and MEQ's `q` *is* `(1/r)∇̄ψ`, so `q·ν` tested
in the plain measure already carries the radius the exterior side carries in its
weight. Writing `dΓ/r` would divide by it twice. The flux is the asset again, for
the fourth time in this tree.

`theTransmissionRowIsTheBoundaryIntegralItClaims` feeds it fields the space
represents exactly — constants `(1,0)` and `(0,1)`, then the linear `(z,r)` —
so the integral collapses to a quadrature of pure geometry and the check is
against a closed form rather than a re-implementation. **7.8e-16 to 1.4e-13**
over three modes. The constants pin the vdof ordering, the sign against
`DarcyForm`'s `−q` and the measure; **the linear rung is the one that exercises
the extension**, because a constant is what an element's polynomial gives at any
point whether or not `TransformBack` found the right reference coordinates, so a
clamped inverse map — the failure `ElementExtension` exists to prevent — is
invisible to it.

**AND THE TILING CHECK WAS RED, ON A DIAGNOSIS THAT WAS WRONG. IT IS GREEN, AND
THE MISTAKE IS THE PART WORTH KEEPING.** `ExtensionBoundaryQuadrature`'s
acceptance is that the boundary weights sum to `|Γ|`. At a 12th-order face rule
MEQ read 1.01e-04, 2.24e-05, 4.59e-06 at `n = 12, 24, 48` — converging at about
`O(h²)` where the expectation is a mesh-independent floor near 1e-10 — and filed
it upstream as lost coverage, on the argument that **"a quadrature residual that
converges is measuring a geometry rather than an instrument"**.

**THAT ARGUMENT HAS A HOLE: TWO THINGS CONVERGE.** The other is a curve the
*rule* under-resolves, which straightens as `h` falls. Refining the **rule** at
fixed `h` separates them, because no quadrature recovers coverage that is not
there. Upstream made that point and MEQ reproduced it on its own circle — at
`n = 12`, cone on:

→ **[M-12](MEASUREMENTS.md#m-12)** — q8 · q12 · q20 · q40 · q80

and **5.40e-10 is the cone-off floor**. So **coverage is exact with the cone and
without it**. What the cone costs is the *smoothness* of `ξ ↦ a(x(ξ))`: it drives
the two interpolated vertex directions apart, the foot map roughens, and a
fixed-order Gauss rule under-resolves it. The rule is now 80 in that case, the
gate stands at its original value, and it is green.

**THE REPAIR THAT WOULD HAVE BEEN WRONG IS RELAXING THE GATE.** It was never too
tight; the rule was too coarse for the path family. Same shape as the Richardson
findings `CLAUDE_INVERSION.md` records half a dozen times over — a column that
keeps moving is measuring the instrument, and the way to tell is to refine the
instrument at fixed geometry.

**AND IT MOVED A DEFAULT IN MEQ'S OWN NEW CODE.**
`transmissionQuadratureOrder` was 12, which is the same rule over the same foot
map, so the transmission row was short by `O(h²)` against a coned family. It is
**40** now. `theTransmissionRowIsTheBoundaryIntegralItClaims` could not have
caught that and now says so: it builds its reference with the *same* rule the row
uses, so the quadrature error is common to both sides and cancels **exactly** —
which is what isolates the ordering, the sign and the measure, and is why it
reads 1e-16 while being blind to resolution. Checking a solve against the formula
it used, one level up.

**AND THE BLIND SPOT DEMONSTRATED ITSELF ON THE WAY, WHICH IS BETTER THAN THE
COMMENT SAYING SO.** Raising the default from 12 to 40 while that case's
reference stayed at 12 un-matched the two rules, and it failed at **5.9e-07** —
the gap between an order-40 and an order-12 sweep of the same integral, which is
exactly the under-resolution being fixed, surfacing as an apparent defect in the
row. Both are pinned to 80 now, and the case records that **same** and **high**
are separate requirements: same is what isolates the contraction, high is what
stops a coarse pin hiding an inadequate shipped default.

**AND THAT DEFAULT REACHED MEQ AND TURNED TWO CASES RED.**
`VertexConePath`'s four-argument constructor now builds **no cone at all**, so
`theConeIsWhatCostsTheTiling` compared a cone-off column against a cone-off
column — its own message said so, *"the two columns are the same experiment"* —
and `theBoundarySweepTilesGammaAndTheRegionSweepTilesTheGap` found the signed
sweep no better than the unsigned one, because signing matters only where the
foot map backtracks and it is the cone that makes it backtrack. **Both cases
were right and the calls were stale**: they now pass `use_cone = true`
explicitly, and the table below reproduces to every digit, which is what says
the repair is the right one rather than merely a green one. A default that
changes under a test whose subject IS that default is the shape to expect again.

**Upstream also turned the cone off by default** (a `use_cone` constructor flag),
having measured that it does not do what it was added for — the aerofoil's flux
order it was meant to fix turns out to be **pre-asymptotic** and recovers on its
own. And they connected it to MEQ's earlier signed-weight finding: **the cone is
what makes the foot map backtrack**, so MEQ's `O(h)` unsigned overcount and this
`O(h²)` residual are two readings of one thing.

**One assertion beside it was wrong and is corrected rather than relaxed.**
`gammaHIsStarShapedAboutTheCentreAndStaysSoUnderRefinement` asserted
`margin > 0` and reads **exactly `0.000000000e+00`** at every mesh. An exact zero
repeated across three meshes is not geometry. `Γ_h` is a union of background
element faces, so it carries axis-aligned ones, and a **horizontal** face with an
endpoint on `z = centreZ` has `(x − c)` purely radial against a purely vertical
`n` — exactly zero, and `MakeCartesian2D`'s diagonal split is not symmetric about
that line, so such a face survives at every mesh. A ray through that corner is
**tangent**, not a ray meeting `Γ_h` twice, and the margin is never negative. The
property wanted is `margin ≥ 0`; the strict inequality was a guess at how to say
it and is false for every staircase. The control — a centre outside the disc,
reading `−1` — is what keeps the `≥` from being vacuous.

**FB-1a IS DONE AND IT SETTLES §8's CORNER RISK.**
`theSolverReachesTheExteriorDatumOnTheHalfDisc` is the first thing in the tree
that **solves** on a semicircle centred on the axis, with the exterior datum
*given* rather than solved — FB-1 minus the border. Rates in `ψ`: **1.99, 2.99,
3.99** at `k = 1, 2, 3`; in `q`: 1.94, 2.67, 3.93. So the corner where `Γ` meets
the axis, where the lifting weight `C = r` vanishes and the plan called the
behaviour *"probably benign, definitely not established"*, is **benign and now
established**. `q` at `k = 2` sitting short is FB-A's half-order axis loss
appearing in a real curved-boundary solve.

**AND IT FOUND THAT MEQ COULD NOT IMPOSE A NON-ZERO DATUM ON `Γ` AT ALL.**
`setBoundaryData()` is projected against `fittedMarker`; `HDGExtensionIntegrator`
supplies only the *solution-dependent* half of the transferred datum, which is
the whole of it exactly when `g` is homogeneous. Every extension study in this
tree happens to want `ψ = 0` on `Γ`, so nothing had ever asked for more.

**THE FIRST FIX WAS WRONG AND WAS INERT, AND WHY IT WAS INERT IS THE USEFUL
PART.** Setting `Γ_h`'s trace dofs and letting `FormLinearSystem` eliminate them
— which *is* how the fitted datum is imposed — changed **not one digit**. The
flux divergence form's boundary face integrator carries `fittedMarker`, and
`EnableHybridization` registers a boundary flux constraint on exactly the
attributes it finds marked there. **`Γ_h` is not among them**, so its trace dofs
are essential in name while nothing couples to them.

**`miniapps/hdg/extension.cpp` IS THE WORKED EXAMPLE AND MEQ HAD NOT READ IT
CLOSELY ENOUGH.** It imposes the datum as a **load term on the flux equation** —
`fform->AddBdrFaceIntegrator( new VectorBoundaryFluxLFIntegrator( datum ),
bdr_gamma_h )` with `datum` a `PathTraceCoefficient` — which is `⟨ψ̂, v·n⟩` of
(8a) where a non-homogeneous `g` belongs. `setExteriorDatum()` does that now, and
**negates**, as the miniapp does (`pNatural = -pExact`, *"the datum as the flux
equation takes it"*): MEQ's flux block holds `−q` for the same reason Darcy's
does, so the conventions coincide.

**SO THERE WAS NO MFEM DEFECT TO FILE**, which is worth recording because a
report was nearly written. The capability exists and is demonstrated upstream;
MEQ had simply never built the half it did not need. The miniapp *does* exercise
a non-zero datum — checked, `pExact` is non-trivial on all three of its problems.

**And the coupled domain must reach `r = 0`**, because the exterior expansion is
only valid on a semicircle centred on the axis. So FB-A is FB-1's prerequisite
rather than its warm-up, and the `O(1/h)` conditioning below is what FB-1
inherits.

FB-A is the one stage that needs no free boundary at all: a vacuum solve on a
mesh whose inner edge is `r = 0`, where the flux mass `(r q, v)` degenerates and
the operator's `1/r` is not integrable.

**IT CAN BE MEASURED AGAINST A CLOSED FORM, WHICH THE PLAN DID NOT KNOW.** Its
acceptance was written as "the trace condition number bounded under refinement",
which is weaker than anything else in this tree and half stale besides — the
other half asked for element-local iteration counts, which are identically zero
under NPC. But there are polynomial `Δ*`-harmonic functions vanishing on the
axis: `r²`, `r² z` and `r⁴ − 4r²z²`, in `tests/analytic/VacuumHarmonic.hpp`.
**A fourth was guessed and was wrong** — `r²(r² − 4z²)z` has `Δ* = −16r²z` — which
is why the fixture recomputes `Δ*` by central differences rather than trusting
the algebra, and why `r⁴/8 → r²` is kept as the control that says the operator
being differentiated is MEQ's.

→ **[M-13](MEASUREMENTS.md#m-13)** — `k` · `ψ`, axis / control · `q`, axis / control · conditioning ratio grows by

**`ψ` is unharmed by the axis and `q` is short by about half an order**, and the
conditioning penalty is **`O(1/h)`** — the on-axis column doubles with every
halving of `h` while the control settles. Over an eightfold refinement that is
7.2 to 7.6 against the 8 a clean `1/h` gives and the 64 a `1/h²` would.

**One mechanism, and it is the weight rather than the singularity.** `(r q, v)`
gives the element touching the axis a weight of order `h`: its diagonal is the
smallest in the system, which is the `1/h`; and it is the element the method
controls least while an unweighted `L2` counts it in full, which is the half
order. The plan's three non-measurements were that `q` is *bounded* and the mass
matrix *positive definite* — both true, both now measured — and neither is what
gives way.

**It does not block FB-1.** MEQ's trace solve is direct, and a direct solver is
nearly insensitive to conditioning at 2D serial sizes; `ψ`, which is what §4's
coupling transmits, keeps full order. `1/h²` would have stopped it.

**`k = 3` escaping is a property of the fixture, not of the method** — its `q`
is a quadratic, so `P_3` has room to spare. A higher-degree flux should show the
deficit there too, and measuring that is the obvious next step.

**Two traps met while writing it, both already recorded in `CLAUDE_HDGGS.md`.**
The conditioning
had to be measured on the **linear** path, because under NPC `reducedOperator()`
is a `DarcyNPCOperator` with no entries — `theTwoPathsAgreeOnTheVacuumField`
pins the two paths at 1e-16 so the conditioning is a statement about the system
the rates were measured on. And the affine-source check had to assert the
**residual drop** rather than `newtonIterations() <= 1`, which failed at the
finest mesh: `CLAUDE_HDGGS.md`, *One more test moved from the stopping rule to the property*, met
again from scratch.

## XP-0: MEQ locates an X-point at the order its flux converges at, unseeded

**`FREE-BOUNDARY-PLAN.md` §10.6's cheapest stage, and the one everything above
XP-1 assumes.** `Soloviev::nstx()`'s X-point is fixed by the twelve
Cerfon–Freidberg constraints at `( 0.699700, −1.716000 )`, so it can be measured
against a closed form with no free boundary, no normalisation and no coupling —
the same shape as FB-A, which is the argument that made FB-A the best-value stage
in that plan.

→ **[M-72](MEASUREMENTS.md#m-72)** — the position sweep · **[M-73](MEASUREMENTS.md#m-73)** — the reference check and the audit

**NOTHING IS SEEDED, AND THAT IS THE PART WORTH KNOWING.** `findAxis()` and
`tryFindAxisFrom()` cannot reach a saddle by construction — both filter on
`AxisSense` — so what reaches the X-point is `sweep()`, which roots every element
from its own centre and its own quietest flux node. On a box holding exactly one
saddle it returns exactly one, at every `k` and every `n` tried, so the located
point carries no prior at all. The coarsest row is the proof: at `k = 1, n = 4` it
lands 3.5e-3 away. `overshoot` is exactly zero at every point, so
`setContainment()` never arbitrates here.

**THE ROOT FINDER ADDS NOTHING, AND HERE THAT CAN BE SAID TWO-SIDEDLY.**
`q_h( x_h ) = 0` and `q( x* ) = 0` give `x_h − x* = −J⁻¹( q_h − q )( x* )` with
`J = dq/dx = Hess( ψ )/r`, which is **symmetric** — a Hessian being so — with
eigenvalues `+0.899099` and `−0.578731` at this saddle. So the position error is
trapped between **1.112** and **1.728** times the pointwise flux error, and over
all fifteen points the measured ratio is **1.11 to 1.73**, touching both ends. The
axis study can only state the upper half: there the window is 0.77 to 3.27 and the
measurement never reaches either end. The test computes the window from the
fixture's own Hessian rather than quoting it.

**THE REFERENCE HAS TO BE CHECKED FIRST AND ONE FIXTURE FAILS THAT CHECK.**
`nstxAsPublished()`'s prescribed X-point carries `|∇ψ| = 2.97e-02` and its saddle
is **9.10e-02** away, at `( 0.695811, −1.806937 )` — exactly what `Soloviev.hpp`
records. A rate study against that point would converge at a clean rate to a place
the finder is right not to be, which is this tree's standing hazard. It ships as
the control in `theClosedFormXPointIsASaddleOfTheClosedForm`, so that "`nstx()` is
the right fixture" is a measurement rather than a preference.

**A FOUR-LEVEL SEQUENCE IS NOT ENOUGH AT `k = 1`, AND THAT IS THE SEQUENCE AND NOT
THE FINDER.** On the axis study's `{ 4, 8, 16, 32 }` the position rate is
**1.685**, under `k+1` less 0.25 of slack; five levels give **1.814** and `n = 128`
gives 1.916 from `n = 4` and 2.017 from `n = 8`. The pointwise error of `q` at the
same point does the same thing in the same places — 1.696 then 1.862 — while
`L2( q )` sits at 1.971 throughout. A pointwise error carries a constant that is
wherever in its element the point falls, so a short sequence measures the ratio of
two of those constants as much as it measures an order. **And there is no sequence
comfortable at all three orders at once**: `n = 128` lifts `k = 1` to 1.916 and
drops `k = 2` to 2.884, because at `n = 64` the X-point happens to fall where `q_h`
is unusually good. So the shipped case asserts the rate **and** asserts that it
agrees with the pointwise flux error's own rate to 0.3 — measured 0.048, 0.022,
0.057 — which is the assertion that separates *the finder is wrong* from *the field
is coarse here*.

**AND `audit.consistent()` CANNOT HOLD OVER A BOX ENCLOSING BOTH, WHICH IS
ARITHMETIC RATHER THAN A DISCRETISATION FAILURE.** Poincaré–Hopf says an
outward-transverse `q` gives sum-of-indices `= χ = 1`; the interior indices are
`+1` at the axis and `−1` at the X-point and sum to 0; so `q` is not
outward-transverse on any such box, at any mesh, and the degree is entitled to
disagree with χ. Measured at three orders: degree 0, χ 1, `transversality 0.00`.
The half of the theorem that needs no hypothesis — degree equals the sum of the
indices inside — holds at every order, and it is the half a diverted
free-boundary solve would use: *the sweep has found everything the boundary says
is in there*.

**WHAT XP-0 EXPOSED FOR XP-2 AND XP-3, AND IT IS CLOSED.** `sweep()` was the only
entry point that reaches a saddle, and it costs one Newton per element.
`tryFindCriticalPointFrom( r, z, AxisSense::Saddle, found )` is the seeded route
— `tryFindAxisFrom()`'s own seed-and-rings contract, one piece of code, with the
axis-named forward refusing `Saddle` so that *axis* keeps meaning axis. It is
measured in `theSeededSaddleSearchCostsLessThanASweep`: the same X-point to
round-off for a **thirtieth to a two-hundredth** of the Newton solves, the ratio
widening with refinement because rings are an absolute element count and a sweep
is the whole mesh.

**AND XP-3 TURNED OUT NOT TO NEED IT, WHICH IS WORTH SAYING BECAUSE THE PLAN SAID
IT WOULD.** The prediction here was that *"XP-3's border relocates the X-point
once per Jacobian"*. It does not relocate it at all: the two rows `q_r = q_z = 0`
ARE the root find, so Newton moves the point and all the border needs per
iterate is a **point location** — which element holds `( r_X, z_X )` — and not a
search. What the seeded search is actually for is XP-2, which uses it every
sweep, and XP-3's own acceptance, which checks the border's answer against an
independent root find on the same field.

## XP-3: the X-point is two more unknowns of the same Newton

**`FREE-BOUNDARY-PLAN.md` §10.4 and §10.6's last rung before a diverted machine
case.** `meq::GradShafranovSolver::setXPointBoundary` is the whole interface, and
what it adds to the bordered system of FB-5 is

```
q_r( r_X, z_X )             = 0
q_z( r_X, z_X )             = 0
psi_bnd - psi_h( r_X, z_X ) = 0
```

with `( r_X, z_X )` unknowns beside `ψ_ax`, `ψ_bnd`, the profile scale and the
exterior coefficients. `tests/convergence/XPointBorder.cpp` is the acceptance and
it shares XP-2's fixture through `tests/convergence/DivertedMachine.hpp` — the
same machine and not a copy of it, because the headline assertion is agreement
between the two.

→ **[M-86](MEASUREMENTS.md#m-86)** — the sweeps · where it closed · against XP-2 · the cost · the `ψ_ax` first step

**THE TWO NEW UNKNOWNS COST NO BACKSOLVE, AND THAT IS NOT OBVIOUS.** The bordered
elimination pays one backsolve per border COLUMN, `c_j = ∂R/∂p_j`, and the field
residual does not contain `( r_X, z_X )`: they reach it only through `ψ_bnd`,
which has a column already. So both columns are **exactly zero**, `z_j = J⁻¹0 = 0`,
and the system grows from `( N + 2 )` to `( N + 4 )` in the dense corner alone.
An X-point costs two rows of a 4×4 and nothing else.

**AND THE CORNER BLOCK IS EXACT, WHERE §10.4 EXPECTED IT NOT TO BE. THE ERROR IN
THE PREDICTION IS WORTH KEEPING BECAUSE IT IS A CLASS OF ERROR.** The plan reads
`∂( q_r, q_z )/∂( r_X, z_X )` as `∇q` — the potential's Hessian — and notes there
is no solved variable for it, since differentiating an L2 field of degree `k`
leaves `k−1`; it calls this *"the same wall recorded for the band continuation of
`B`"* and plans a differenced fallback. The analogy does not carry. The band
continuation needs `∇q` as an approximation of the **continuous** Hessian, and
there it really is an order down. Newton needs the derivative of the **discrete**
residual with respect to the discrete unknowns — and `q_h` is a polynomial on its
element, so `∂q_h/∂x` there is exact arithmetic, not an approximation of
anything. Nothing on this path is differenced, the fallback was not needed, and
the same goes for `∂ψ_h/∂( r_X, z_X )` in the `ψ_bnd` row. *Approximating a
continuous object and differentiating a discrete one are different questions, and
the second is the one a Jacobian asks.*

**THE `ψ_bnd` ROW'S TWO NEW ENTRIES ARE NOT `r q`.** `ψ_h` and `q_h` are separate
solved fields whose identity `∇̄ψ = r q` is only weak, so the row takes the
potential's own derivative from the same element. Both entries vanish at a
converged X-point either way, which is exactly why writing them costs nothing and
leaving them out would be invisible until the order was measured.

**WHAT IS GENUINELY NOT SMOOTH IS THE ELEMENT CHANGING.** `q_h` is discontinuous
across a face, so carrying the point over a mesh line changes the polynomial the
rows evaluate, by `O( h^{k+1} )`. The element is therefore frozen within a
Jacobian and re-decided at each accepted step — the discipline the axis row and
the limiter contact already keep — and a point a hair outside its own element is
evaluated **in** it rather than refused, so a step across a face is a jump in the
residual rather than a failure to have one. `locateFieldPoint()` is that, and it
uses MFEM's **unprojected** inverse map: the default `NewtonElementProject`
clamps, and a clamped point reports an overshoot of zero however far outside it
really is, which is the clamped-inverse-map failure `CLAUDE.md` already records.

**IT IS NOT CHEAPER IN NEWTON STEPS THAN XP-2 AND THE WIN IS STRUCTURAL** — 21
steps against 18 on the same fixture in the same process, M-86. What it buys is
that the X-point stops being an **outer state**: §10.5 names two on XP-2's fixed
point, which critical point bounds the plasma and which elements carry current,
and XP-3 removes the first by making the position differentiable, which it is. It
does not touch the second, which is not, so this case still freezes the support
and still re-decides it between solves.

**THE OBSERVED ORDER IS THE FIXTURE'S CEILING AND NOT THE BORDER'S, AND ONLY THE
CONTROL COULD SAY SO.** XP-3's bootstrap tail reads 1.664 and XP-2's — the same
bordered system less exactly these two rows — reads 1.667. Neither is quadratic,
and what caps them is not known: the candidates are the axis argmax hopping
elements, an `O( h^{k+1} )` kink, and the element-granular plasma support. The
acceptance is therefore a **comparison** rather than an absolute bound, which is
the only honest reading. A clean quadratic run on a bordered free-boundary solve
is `LimiterCurve`'s, on a fixture with none of that.

**AND THE FIRST NEWTON STEP OF EVERY BORDERED SOLVE IS TAKEN WITH THE `ψ_ax`
BORDER DECOUPLED.** Found while looking for the order ceiling and measured with a
trace at the border build: `constraintLocated` was last written by the `peakAt()`
that computes the convergence target, which runs on the **cold** state where no
axis exists — so iteration 0 builds the empty row and every iteration after it
builds the real one. **Repairing it is three lines and it changes which
equilibrium the diverted machine reports**: XP-3 goes 14 steps to 10 with the
same answer in every digit, and XP-2 goes 7 steps to 3 and lands on a different
branch, `ψ_ax` 8.052272e-02 against 8.266004e-02 with the X-point 1.26e-02 m away.
The decoupled step is `ψ_ax ← max ψ_h`, which is conservative, and on this machine
that is what keeps the iteration in the physical branch's basin. So it is a
branch-selection question rather than a Jacobian one and it is **left alone** —
recorded in M-86 and beside the code, so that the next person to find the stale
flag does not repair it without measuring the second row.

## XP-4: the diverted machine from a file, and against `freegs4e`

**TWO KEYS, AND THE SECOND IS THE ONE THAT IS NOT OBVIOUS.**
`[boundary.xpoint] R, Z` seeds `setXPointBoundary()` — and it is a SEED, which is
the whole difference from `[boundary.limiter]`: the driver carries it forward
across adaptive cycles as it carries `psiAxisGuess`, echoes nothing of it back
into the `.nc`, and reports how far the solve moved it.
`[solver] PlasmaSupportSweeps` is the other half, and it is the outer loop XP-3
deliberately left standing: freeze the support, solve, refreeze at the answer,
solve again, stop when a sweep changes nothing.

**THE LOOP'S FIRST SWEEP NEEDS A `ψ_bnd` THE SOLVER CANNOT YET SUPPLY.**
`psiBoundary()` is the CONVERGED value and there is none before the first solve,
so the driver reads `ψ` of the initial iterate at whichever bounding point the
file named — the X-point seed, the limiter contact, or the maximum over the
meshed limiter surface. Getting the iterate is `prepare()`: `setInitialGuess()`
only records the guess, and what puts it into `potential()` is the projection
`prepare()` does, which `solve()` would do a moment later anyway. That costs one
assembly per cycle and is the only route that works for the ramp and the bump,
which arrive as `Coefficient`s and have no field to read.

**AND THE FREEZE IS LEFT IN FORCE AFTERWARDS, DELIBERATELY.** The answer is the
solution of the problem the LAST sweep posed — that support, that threshold — and
everything downstream reads the source: the post-processing, the estimator, the
flux surfaces, the axis check. Thawing would evaluate `F` on a support the
residual was never driven to zero on, which on a settled loop is the same support
to round-off and on an unsettled one is a different problem wearing the answer's
numbers. So the fixed-point question is ANSWERED — `plasma_support_settled` in
the `.nc` — rather than smoothed over by a last refresh nobody solved with.

**THE DRIVER REPRODUCES `XPointBorder`'s IN-PROCESS ANSWER TO EVERY PRINTED
DIGIT**, per-sweep iteration counts included, which is the check that says the
TOML route and the API route are one computation.

**AGAINST `freegs4e`, AND THIS COMPARISON IS SHARPER THAN THE LIMITED ONE.**
`theDriverSolvesALimitedTokamak` hands BOTH codes the same limiter point, because
freegs4e's own contact is a maximum over grid CELLS that moves 0.6 m between
129² and 513²; that case deliberately takes the contact-finding out of the
comparison. Here nothing is handed over — MEQ solves for the null as three rows
of its Newton, from a seed 7.07e-02 m away, and freegs4e finds it by a
critical-point search — so the agreement in its POSITION is a result.

→ **[M-87](MEASUREMENTS.md#m-87)** — X-point 4.4e-04 m · `ψ_ax` 6.9e-04 · `ψ_bnd` 7.7e-04 · the field by region

**THE FIELD COMPARISON FOUND THE CONDUCTOR MODEL AND NOTHING ELSE.** Over every
comparable node the relative `L∞` is **2.973e-01**, which read alone would be a
failure; **all of it is inside two coils**. `P2L` and `P2U` are freegs4e
FILAMENTS — a point source with a logarithmic singularity — and MEQ meshes them
as 0.1 × 0.1 m rectangles carrying a uniform current density, so a node inside
one compares two conductor models rather than two solvers, and no refinement of
either closes it. Excluded with a 5 cm collar, 5478 nodes agree at **5.3e-04**
relative `L2` and **1.8e-03** `L∞`, worst on the outboard midplane near the
plasma edge rather than anywhere near a conductor. `P1L` and `P1U` are
`ShapedCoil`s and agree by construction — and contribute no nodes at all, sitting
outside the reference's own box.

**`compare.py --exclude-box` REPORTS THE EXCLUDED REGION AS A ROW OF ITS OWN
rather than dropping it**, which is the point: a comparison that quietly threw
away the nodes where it disagrees would be the instrument choosing the answer.
`--free-boundary` is the other new flag and says there is no gauge shift and no
`-meta.json` — both codes solve `ψ → 0` at infinity, MEQ through
`meq::ExteriorDtN` on `Γ` and freegs4e through von Hagenow, so `ψ_ax` and
`ψ_bnd` are directly comparable numbers rather than differences, and the MXH fit
that floors the fixed-boundary rehearsal at 2–4e-04 m is absent entirely.

**~7e-04 IS THE FLOOR THIS COMPARISON CAN MEAN AND THE BOUNDS SAY SO.**
`fgsref.py` fits a `UnivariateSpline` to its own analytic profile shape before
solving and the fit MOVES it — **1.707e-02** of the amplitude in `ff'` — while
the TOML's tables are the analytic shape. The two codes are therefore solving
sources that differ at the per-cent level, and an agreement much tighter than
what is measured would be evidence of a shared mistake rather than of two right
answers.

## PE-0: the premise is true, and the plan's own machinery does not cash it

**`PLASMA-EDGE-PLAN.md` §6's first stage is run, and it splits the plasma edge's
order loss into two halves that the plan had welded together.**
`tests/convergence/PlasmaEdgeConvergence.cpp` solves `PlasmaEdge` on `Ω_{p,h}` —
the elements lying entirely inside the disc — with `λ` given, which is stage 5's
curved-boundary machinery pointed at an interior circle. **The premise is true.**
On a fitted rectangle strictly inside the plasma the same equilibrium reads `k+1`
and `k+2` to two decimal places at `j = 0`, `1` and `2` alike, and on `Ω_{p,h}`
itself `ψ_h` climbs from the cut mesh's cap of `min( k+1, j+1.5 )` — **1.5 at
`j = 0` whatever `k` is** — to **2.00 / 3.04 / 4.01**. `j` leaves the rates:
`ψ*`'s spread over the three rungs falls from 1.11 / 2.04 / 2.18 on the cut mesh
to **0.07 / 0.18 / 0.29** on the uncut one.

→ **[M-74](MEASUREMENTS.md#m-74)** — the rate against `j`, both routes, against the cut mesh · **[M-75](MEASUREMENTS.md#m-75)** — the fitted-rectangle control

**WHAT IS NOT TRUE IS THAT THE PLAN'S OWN MACHINERY CASHES IT, AND THE EXPERIMENT
THAT SAYS SO IS FREE.** Each row is run twice — with `λ` on `Γ_{p,h}`, where it is
the exact trace and the problem is consistent, and with `λ` on `Γ_p`, transferred,
which is what PE-3 must do because there `λ` is the interface unknown and lives
nowhere else. `ψ*` reaches **`k+2` on the first and never on the second**: 2.88 /
3.94 / 4.93 against 2.23 / 3.16 / 4.00 at `j = 0`. The two differ only in the
transfer, so **the missing order is the transfer's and not the plasma edge's** —
the same shortfall `ExtensionConvergence` records on the OUTER boundary, where it
asserts `k+1.5` rather than `k+2` for exactly this reason. **PE-0's acceptance as
`PLASMA-EDGE-PLAN.md` §6 writes it is therefore not reachable, and it is the
acceptance that has to move.** The case is red on that assertion and stays red,
per *Testing stance*.

**AND THE MESH-DEPENDENT FRAGILITY OF §9.4 IS NOT THE CORNER.** That section
measures a transferred solve failing on particular meshes and reads the X-point's
corner as the cause. PE-0 meets the same thing on a **smooth circle**: at `k = 3`
the transferred `L2(ψ)` is 3.39e-07 at `n = 16` and **5.89e-07 at `n = 24`**, with
`VertexConePath::NumWidened()` zero at every mesh, `dist( Γ_{p,h}, Γ_p )/h` in
[1.02, 1.33] so assumption P.1 holds, and a path quadrature raised to order 24
changing the sixth figure. **The corner sharpens the fragility rather than causing
it.** A per-pair rate is not assertable on the transferred route, which is why
PE-0's pairwise tier asserts that the error FALLS.

→ **[M-76](MEASUREMENTS.md#m-76)** — the non-monotone transferred sequence

**AND THE `j = 2, k = 3` ANOMALY IS THE STAIRCASE.**
`theCutCapsTheOrderBeforeAnyMethodIsChosen`'s uncut column falls, 4.67 then 4.50,
which `PLASMA-EDGE-PLAN.md` §7 records as the one piece of evidence against the
premise. It does not reproduce on the solve: `Ω_{p,h}`'s per-pair `ψ*` reads 4.70,
**4.50**, 4.74 — scatter, and rising rather than falling — while the fitted
rectangle reads 4.94 on the same equilibrium at the same `j`. Raising the source
rule from `2k+4` to `2k+16` moves the column in **no digit**, so it is not the
quadrature. What is left is `Γ_{p,h}`'s re-entrant corners — 225° and 270° on a
simplicial background, **[M-84](MEASUREMENTS.md#m-84)** —
where HDG's duality argument for `k+1` and `k+2` wants an `H²`-regular adjoint
and does not get one. The lever this tree already has is `meq::AdaptiveDomain`,
GS-2 §3.3's companion mesh.

## `ψ_bnd`: settable, an unknown, and reachable from a file

**This is FB-3.** It lives beside `ψ_ax`'s border, which it generalises, and
that border is described in `CLAUDE_HDGGS.md` under *Newton, and the obligation
it creates* — worth reading first, because `ψ_ax` is an unknown for a
*fixed*-boundary run with normalised profiles and everything below is what it
took to make `ψ_bnd` one too.

What follows is how `ψ_bnd` became an unknown, kept because the traps in it
recur.

`ψ_bnd` **is settable AND is an unknown.** `setNormalisation( ψ_ax, ψ_bnd )` is
the interface and the profiles take
`Ψ = (ψ − ψ_bnd)/(ψ_ax − ψ_bnd)` rather than `ψ/ψ_ax` — MEQ's
fixed-boundary problem has `ψ = 0` on `Γ`, so `ψ_bnd` vanished and the
assumption was baked into the class.

**The validation moved from `ψ_ax` to the SPAN, and that is the whole
generalisation.** With a boundary flux in hand `ψ_ax = 0` is perfectly
admissible; `ψ_ax = ψ_bnd` is not. The one-argument form is kept as a non-virtual
convenience passing zero, so every fixed-boundary caller is untouched.

**Asserted bit for bit rather than to a tolerance.** `Ψ` depends on `ψ` and
`ψ_bnd` only through their difference, so `F( ψ ; ψ_ax, ψ_bnd )` and
`F( ψ − ψ_bnd ; ψ_ax − ψ_bnd, 0 )` are the *same arithmetic* —
`the_boundary_flux_enters_only_through_the_span` requires `0.0`, which is what
catches a `ψ_ax` left where the span belongs. A tolerance would let that through
wherever `ψ_bnd` happened to be small. It carries a control, since a class
ignoring `ψ_bnd` outright would pass the translation check, and `dFdPsi` is
checked separately because it carries **two** factors of the span and a missed
one costs only the convergence, not the answer.

**And `HighBetaConvergence` is unchanged to every digit** — `ψ_ax` at
3.058984e-01 and 2.834510e-01 **as they then were, under
`AxisConstraint::NodalMaximum`**, `ψ_ax − max ψ_h` at 0.00e+00 and −5.55e-17, 4 and
6 Newton iterations — which is what says the generalisation reduces exactly.

**The analytic fixtures REFUSE a non-zero `ψ_bnd` rather than ignoring it.**
Every equilibrium behind `ConvergenceHarness` is written for `Ψ = ψ/ψ_ax`, and
silently dropping a boundary flux would converge beautifully to the wrong
equilibrium.

**AND THE SECOND BORDER IS CLOSED — FB-3 IS DONE.** `setBoundaryFluxPoint( r, z )`
makes `ψ_bnd` an unknown pinned by `ψ_h` at the nearest potential dof to the
limiter contact, and `solveWithNormalisation()` does a general `( N + 2 )`
elimination of which this is the `2×2` case. It is the cheaper of the two borders,
as predicted: `ψ_bnd`'s dof is fixed at setup where `ψ_ax`'s needs an argmax, so
under NPC the row is exactly `−e_j` and the corner exactly 1, neither differenced.
`HighBetaConvergence` drives it and reads `ψ_ax`'s constraint at **1.88e-16** and
**0.00e+00** on two meshes. **NPC only, and refused otherwise** — under the
condensation `ψ` is a function of the trace through every element's source, so
both the row and the corner would have to be differenced.

**THE DRIVER HALF IS DONE, AND MEQ SOLVES FREE BOUNDARY FROM A
CONFIGURATION FILE.** `[boundary.limiter]` gives `setBoundaryFluxPoint()` its
point and `[boundary.exterior]` gives `setExteriorCoupling()` its DtN, so a
machine case is no longer a library caller. `examples/free-boundary-halfdisc.toml`
converges in **7 Newton steps**, `ψ_ax = 9.758655e-02`, with four Gegenbauer
coefficients solved for beside it, and `theDriverReachesTheExteriorCoupling`
pins the driver against the library at **7.0e-17 over 7998 dofs**.

**AND ITS CONTROL IS WHAT MAKES THAT MEAN ANYTHING**: the same problem with the
coupling removed — a zero datum on `Γ_h`, which converges perfectly well — moves
`ψ` by **76.6%** in L2. A coupling that had silently done nothing would have
passed a driver-against-library check, because both sides would have done
nothing.

**`[boundary.exterior]` DEFINES `Γ` ITSELF, WHICH `[boundary.shape]` CANNOT.**
`meq::BoundaryShape` refuses a surface reaching `r ≤ 0` — rightly, a closed
plasma surface through the axis carries a non-integrable `1/r` — and this `Γ` is
a **semicircle whose flat side IS the axis**. So the two blocks are alternatives
and naming both is refused. The driver refuses `[mesh] RMin ≠ 0` at startup, with
no tolerance, for the reason `halfdisc.py` does: a domain stopping at `r = 0.05`
is not a slightly worse semicircle, the Gegenbauer modes do not span its
exterior at all.

**AND `[boundary.exterior]` BESIDE `[mesh] File` DID NOT RUN AT ALL UNTIL
 — WHICH IS FB-6's OWN CONFIGURATION.** A gmsh half-disc with the
conductors meshed to, plus the exterior coupling, is what `tools/mesh/halfdisc.py`
and §7.9 exist for, and it was the one path on which **neither** of the two
preconditions `[boundary.exterior]` enforces was actually enforced:

* the axis test read `[mesh] RMin`, which a file mesh leaves at its **default of
  zero**, so it passed vacuously on a mesh nobody had looked at. A `.msh` whose
  inner edge sat at `r = 0.05` would have gone straight through the one check
  written to stop it, and converged at full order to a machine nobody described;
* the radius test compared `Γ` against `[mesh] RMax`, also zero, so it **refused
  every file outright** — with a message about a box the run does not have;
* and `backgroundCellSize()` computed `( 0 − 0 )/( 0 · 1 )`, so the transfer
  path's search length was **zero** and `mfem::VertexConePath` aborted on the
  first vertex of `Γ_h`, several frames deep in MFEM with nothing naming a
  configuration key.

All three read the mesh now — `mfem::Mesh::GetBoundingBox` for the two
preconditions, the largest element diameter for the search length, which is what
the adaptive path already uses. The built-box branch computes exactly what it
computed before, bit for bit, so no existing configuration moves. **The
transferable part**: `[mesh] File` had been wired for the FITTED path only —
`theDriverTakesItsGridFromAMeshItDidNotBuild` covers the output grid — and every
CURVED-path quantity was still reading keys that describe a box the run never
built. Expect more of them wherever a driver derives a length from
`[mesh] RMin .. ZMax`.

**AND `buildSubdomain()` HAD TO BE RELAXED EXACTLY AS `AdaptiveDomain` WAS**, for
the same geometry and by the same rule — some boundary must be **generated**
rather than none inherited. The half-disc's flat side is the box's `r = 0` edge,
which is fitted boundary wanting no transfer.

**A `bump` INITIAL GUESS WAS NEEDED AND IS NOT A CONVENIENCE.** A ramp is
antisymmetric in `z` and describes no plasma; it exists to keep `ψ = 0` off the
trivial branch. Once the boundary is free the guess **chooses which equilibrium
is reported**, so it is part of the problem statement — the same multiplicity the
freegs4e rehearsal measures at 9.4% across three solve routes. Measured here: a
ramp does not converge on this problem at all, and the bump does in 7 steps.

**A NORMALISED PROFILE TABLE HOLDS `dp/dΨ`, NOT `dp/dψ`, AND WITH A PRESCRIBED
CURRENT THE DIFFERENCE IS INVISIBLE.** `meq::NormalisedMHDSource::f` evaluates
`F = scale·( μ₀r²p′( Ψ ) + gg′( Ψ ) )/span`, so the file is differentiated with
respect to `Ψ` and converting another code's `dp/dψ` arrays means MULTIPLYING by
`ψ_ax − ψ_bnd`. That is `examples/rotating-density.dat`'s trap, and
`[source] PlasmaCurrent` **hides it completely**: both profiles carry the same
wrong factor, the scale is an unknown, and the border absorbs it. Measured on the
§7.16 machine case with the tables a factor of `span` too small, the run
converged to the right equilibrium in every digit and the only tell was a scale
of **6.713e-02 against a span of 6.689e-02** where it should be `O(1)`. **So
read the reported profile scale**: on a correctly converted table it is 0.9989.
With the amplitude fixed instead there is nothing to absorb it and the run is a
plasma fifteen times too weak — which for a while looked like a failure of the
amplitude-fixed formulation and was not.

**AND A TABULATED NORMALISED PROFILE MUST COVER THE RANGE THE ITERATE VISITS.**
`meq::SplineProfile` **clamps** outside its knots — deliberately, since a linear
extrapolation of a steep edge profile turns an overshoot into a NaN — but a table
that STOPS at `Ψ = 0` does not describe the same source as one that continues:
it switches `F` off outside the plasma while leaving `∂F/∂ψ` discontinuous there,
which is a moving support imposed by the *table* rather than by the solver. Three
spellings of one configuration, all of which parse:

→ **[M-32](MEASUREMENTS.md#m-32)** — table over `[0, 1]`

The third is a **different branch** — four times the axis flux on the same
profiles — not a better answer.

**AND THE TWO BORDERS TOGETHER WORK, WHICH THE RECORD SAID THEY DID NOT.**
§7.13 recorded `ψ_bnd` converging alone, the exterior coefficients converging
alone, and the combination failing — "the one combination still open". **Every
one of those attempts predates the `ψ_bnd` repair**, and a border on a quantity
the Jacobian could not see is exactly the border that would fail. **Nobody
re-ran it after the fix.** Re-measured, four limiter positions:

→ **[M-33](MEASUREMENTS.md#m-33)** — limiter `R` · Newton · final residual · `ψ_ax` · `ψ_bnd`

Machine zero at every one, `ψ_ax`'s own constraint at 1e-17, and `ψ_ax > ψ_bnd`
throughout — so it is a basin rather than a lucky point.
`theTwoBordersConvergeTogether` guards it and its message names the repair as
what to look at if it ever fails again.

**THE TRANSFERABLE PART IS THAT A FAILURE MEASURED UNDER A DEFECT IS NOT A
PROPERTY OF THE METHOD.** This file's standing rule — re-measure a claim after
the thing it was measured against changes — has always been applied to
*successes*. It applies to failures identically, and this is the first time that
has cost anything: the one item standing between MEQ and a machine case had been
fixed for hours and was still being recorded as open.

The shipped example still leaves `ψ_bnd` at zero, because an example should
demonstrate one thing; `[boundary.limiter]` beside it is a two-line addition and
the case above is where the numbers live.

`Globalisation` other than `None` is refused on this path, loudly: the KINSOL
paths drive a residual of their own and the Picard ones build no Jacobian to
border. **No ordering is refused, and that has been settled twice.** A guard once
refused `NonlinearOrdering::LineariseThenCondense`, written from
`darcyhybridization.hpp`'s summary rather than from the code; it was removed when
both orderings were measured to the same `ψ_ax`, and MFEM has since deleted that
mode outright. Under `NonlinearOrdering::NPC` the question does not arise at all
— `ψ` is an unknown, so **the border row is exactly `−e_j` and the corner is
exactly `1`**, neither of them differenced. See `CLAUDE_HDGGS.md`, *The NPC port*.

## `ψ_ax` was the weak point of the whole free-boundary path

**Two defects wearing one symptom.** The other half of the story is in
`CLAUDE_INVERSION.md` under *`CriticalPointFinder`'s axis and
`GradShafranovSolver::psiAxis()` are THE SAME POINT* — the two definitions of the
axis, what separates them, and the search that finds one of them.
`FREE-BOUNDARY-PLAN.md` §11 is the plan-side record.

**AND THAT DEFINITION IS THE WEAK POINT OF THE WHOLE FREE-BOUNDARY PATH, WHICH
TURNED UP THREE INDEPENDENT WAYS IN ONE NIGHT.** Nothing in *the largest nodal
value of `ψ_h`* says that value is a magnetic axis, and
`G = ψ_ax − max ψ_h = 0` is satisfied at machine zero by a spurious nodal spike
exactly as it is by an axis:

* **A machine case converged to a spike.** Every border at machine zero, the
  prescribed current delivered to seven figures, and `ψ_ax` **twenty-nine times
  too large** from a single dof of one element — with `ψ*` on the output grid
  peaking thirty times lower. The runaway is self-consistent: a spurious `ψ_ax`
  inflates the span, `Ψ` collapses over the real plasma, and the current border
  raises the profile scale by 980 to hold `∫F/r` at `μ₀I_p`. What reached it was
  a profile table wrong by a factor of the span; `FREE-BOUNDARY-PLAN.md` §7.16.
* **`ψ_ax` ON A PUBLISHED TABLE IS THE `r = 0` LAYER.**
  `theTwoBordersConvergeTogether`'s converged answer at limiter `R = 1.20` attains
  its `ψ_ax` in an element **touching `r = 0`**, reading **1.0916e-01** there
  against **4.4472e-02** as the largest `ψ_h` anywhere off the axis, a factor of
  **2.5** — and the bordered Newton was constraining that value. **This entry
  said "a corner artefact" and named FB-1a's corner; §11.1 measured that it is
  not a corner at all** — the largest nodal values are a flat layer along the
  **whole** axis, 168 dofs at `r = 0` whose largest ten agree to **3.5e-05** of
  1.09e-01, and the corner wins the argmax by 2e-05 and lands on the *other*
  corner at `k = 3`. It does not fall with `h` either: 1.0916e-01, 1.0982e-01,
  1.0953e-01 at `n = 24, 32, 48`, against a datum of zero imposed on that very
  boundary.

  **`ψ_ax` NO LONGER REPORTS THE LAYER, AND THE LAYER IS STILL THERE.** Under
  option 3 the same case reads `ψ_ax = 1.2595e-01` at the located axis while the
  `r = 0` layer sits at **3.08e-01** — so the definition is repaired and the
  FIELD is not, which is exactly the split §11.3 and §11.5 describe and why the
  two guards are separate. `theTwoBorderSolveReportsATrueMagneticAxis` asserts
  the field's half.
* **The sweep saw the same corner on the toy fixture**, reading `ψ_ax` at 8.12e-02
  against a field maximum of 2.50e-02 — ratio 0.31 — and 0.36 and ≈ 0 on two
  others. §7.18 item 3.

**AND OPTION 3 IS BUILT: `ψ_ax` IS THE FLUX AT THE LOCATED MAGNETIC
AXIS.** `setAxisConstraint()`, default `AxisConstraint::LocatedAxis`, constrains
`ψ_ax` at a zero of `q_h` rather than at the largest nodal value.
`AxisConstraint::NodalMaximum` is kept as the **control**, because every number
published here before that date was measured with it.

**IT IS FREE IN THE JACOBIAN AND THAT IS THE ENVELOPE THEOREM** — `∇ψ_h( x* ) = 0`
at a zero of `q_h`, so the position term vanishes and the row is the potential
shape functions of `x*`'s element, exact and undifferenced. **The warm start is
what makes the SEARCH free**: `CriticalPointFinder::tryFindAxisFrom()` seeds from
the previous iterate's axis, roots **one element** when the seed is still on the
answer and 19 of 2048 when it is an element away, so the cost is bounded by a
ring count rather than by the mesh and the ratio against a sweep *improves* with
refinement. `rotating-normalised` went 5.51 s → **2.25 s** when the sweep was
replaced by it. **NPC only, refused under the condensation.**

**TWO TRAPS, BOTH CAUGHT BY A TEST.** Re-inverting the map with `TransformBack`
is wrong *and* unnecessary — a zero of a discontinuous `q_h` can lie outside its
own element, where the inverse does not converge (measured, a residual of
**5.1e-02 on an element of 5e-02**), and the finder has the reference coordinates
exactly. And **stopping at the first root found puts the row on the WRONG
element's dofs** — the same physical point across a face is a different element,
harmless for a position and fatal for a row of shape functions.

**WHAT IT MOVED**: every bordered `ψ_ax`, by the `O( h² )` gap between the
definitions — `limited-tokamak` 9.455354e-02 → 9.466087e-02,
`rotating-normalised` 1.039163e-01 → 1.039325e-01 — and all of them now read a
normalised flux of **1.0000** at the located axis, which they did not before.

**AND IT WEAKENS `checkAxis()` DELIBERATELY.** `Ψ` there is now nearly 1 by
construction, so that guard is close to checking a solve against the formula it
used; it still catches a run with **no** O-point, and one where the largest-`Ψ`
O-point is not the one the constraint followed. **The guard that now carries the
weight is `checkAxisSource()`** — is the toroidal current density bounded on the
axis — which is untouched by any of this. The split is the right one: §11.3's
defect is about the FIELD, §11.5's about the DEFINITION, and neither guard can be
satisfied by repairing the other.

**`FREE-BOUNDARY-PLAN.md` §11 IS THE LIST AND IT IS NOW WORKED THROUGH.** Its
first item was that **the guard had never been run on the §7.12b sighting** — if
that spike were itself an O-point of `q_h` then `Ψ` would read ≈ 1 there and the
check would AGREE with it, the check being one-sided and largest-`Ψ`-wins by
design, so it misses rather than false-alarms.

**§11.1 IS DONE, AND THE GUARD CATCHES IT.** `Ψ` at the located
O-point reads **0.56 / 0.51 / 0.35 / 0.52** at §7.12b's four limiter radii where
1 is required — REFUSES at every one, and at `n = 24, 32, 48` and `k = 2, 3`
alike, so it is not one mesh's accident. **The reason is structural rather than
luck**: `ψ_ax` is attained at `r = 0` **exactly**, on the domain boundary, where
an interior extremum cannot be, and every maximum `sweep()` finds sits at
`r ≥ 0.76`. The located O-point at `R = 1.20` reads 4.420624e-02 at
`( 1.411, −0.001 )` against the 4.4472e-02 the nodal argmax found off the axis —
a root of `q_h` and an argmax reaching one feature by two routes.

**AND §11.3 IS DONE TOO: THE CAUSE IS A `1/r` POLE IN THE LOAD, AND `F/r` IS
`μ₀ j_φ`.** `meq::SourceIntegrator` assembles `−( F/r, w )`, and
`j_φ = r p′( Ψ ) + gg′( Ψ )/( μ₀ r )` — so a finite current density on the
symmetry axis **requires `F( 0, z ) = 0`**, and `F = μ₀r²p′ + gg′` leaves only
`gg′` there. `p′` is protected by its own `r²`; `gg′` is not.

**WHICH `Ψ` THE AXIS SITS AT IS THE WHOLE OF IT.** `ψ( 0, z ) = 0` exactly, the
flux through a circle of vanishing area, so `Ψ_axis = −ψ_bnd/span`. **On a fixed
boundary `ψ_bnd = 0` and the axis sits at `Ψ = 0`, where every profile in this
tree vanishes** — which is why nothing had ever met this. **FB-3's limiter border
makes `ψ_bnd` an unknown, it comes out positive, and the axis is then at negative
`Ψ`: in the VACUUM**, where the physics is `g = const` so `gg′ = 0`, and where an
unconfined profile extrapolates and returns `0.05·Ψ_axis` instead.

A 2×2 factorial on that fixture, one variable at a time:

`k = 2`, 1333 elements, limiter `R = 1.15`:

→ **[M-43](MEASUREMENTS.md#m-43)** — limiter · `gg′` · `ψ_bnd` · `Ψ_axis` · `F( 0, z )` · verdict

**THE TWO LIMITER ROWS MOVE AND THE RADIUS WITH THEM**, from
9.21e-03 and 2.13e-02 at `R = 1.20`: `ψ_bnd` stopped being snapped to the
nearest potential dof, and an `O( h )` change in `ψ_bnd` is enough on this
fixture to move it onto the other branch and flip its **sign**. At 1.20 the
`gg′ = 0` row then put the axis INSIDE the plasma, which is the one thing the
cell may not do — so the radius moved to 1.15, where both limiter rows keep
their axis in the vacuum. The mechanism the table demonstrates is unchanged and
the sweep that settled the radius is beside the cells.

**The fourth row rules out the limiter itself** — `Ψ_axis = −0.22`, deep in the
vacuum, and clean, because `F( 0, z )` is machine zero. It is neither the limiter
alone nor `gg′` alone; it is `F( 0, z ) ≠ 0`.

**AND THE DISCRETE HALF IS WHY IT IS A DEFECT RATHER THAN AN UGLINESS.** The
*continuous* problem is well posed: the energy `∫( 1/r )|∇̄ψ|²` forces its members
to vanish faster than `r` — the physical `ψ ~ r²` — and against such test
functions `∫( gg′/r ) w` converges. **The discrete space is `L2` polynomials,
free to be nonzero at `r = 0`, and against those the load functional is
UNBOUNDED**; the quadrature is the only thing making it finite. Sweeping
`setSourceQuadratureOrder()` at fixed `h` is what proves it: with `gg′ = 0` the
answer is **bit-identical at `extra` = 4, 8, 16, 20** — 1.034602461e-05 on the
axis, 1.175887576e-01 off it, ten digits, `F/r` being a polynomial then — and
with `gg′ = 0.05` nothing settles (1.09e-01, 1.15e-01, 9.18e-02, 8.76e-02), nor
does it fall with `h`.

**THREE CONTROLS SAY IT IS NOT THE GEOMETRY, THE CORNER, OR THE AXIS.** FB-A's
fixtures are `Δ*`-harmonic, `F ≡ 0`, no load at all. **FB-1a runs on the same
mesh, extension and corner** and converges at 1.99/2.99/3.99, because
`ExteriorMatched`'s `F` is built on `ExteriorDtN::basis`, which carries
`( 1 − μ )( 1 + μ )` explicitly *"so the axis is exactly zero"* — `F ~ r²`. And
`examples/free-boundary-halfdisc.toml`, with the **same `gg′` table** and the
**same `RMin = 0`**, is healthy at `Ψ = 1.0036`, because it has no limiter and so
keeps `ψ_bnd = 0`.

**SO THE REPAIR IS PHYSICAL AND ALREADY EXISTS: `[source] ConfineToPlasma`**,
which sets `F = 0` wherever `Ψ ≤ 0` — the statement that the vacuum carries no
current. `examples/limited-tokamak.toml` has every ingredient (a limiter, a free
`ψ_bnd`, a gmsh mesh reaching `r = 0`), sets it, and reads `Ψ = 1.0016`. **For a
domain reaching the axis with `ψ_bnd` free it is a PRECONDITION, not an option**,
in the same sense as `j ≥ 1` at the plasma edge.

**SO `ψ_ax` IS TWO DEFECTS WEARING ONE SYMPTOM**, which is the thing to carry
forward: on §7.12b's case the **field** is wrong on a boundary layer and the
argmax merely reports it, while on §7.16's the field is sound and only the argmax
is not. §11.5's three redefinitions address the second and would **hide** the
first.

## §7.12b's fixture was the defect, and what a red case there means

**A RED `theTwoBorderSolveReportsATrueMagneticAxis` IS THE INTENDED SIGNAL WHILE
THE PROPERTY IS ABSENT** — it asserts what is WANTED and fails until it is there,
which is `CLAUDE.md`'s *Testing stance*. **What it turns on is the fixture,
rather than `ψ_ax`'s
definition**: a fixture can ask for an equilibrium that does not exist, and the
repair is to give it the physics it is missing rather than to relax anything.

**THREE THINGS WERE MISSING AND ALL THREE ARE NECESSARY**, measured
one at a time on the same fixture:

→ **[M-44](MEASUREMENTS.md#m-44)** — removed · what happens

**AND THE COIL CURRENT IS DERIVED RATHER THAN TUNED**, which is §7.15's finding
applied: Shafranov's `B_v = μ₀I_p/(4πR)·[ln(8R/a) + β_p + l_i/2 − 3/2]` says what
field a given `I_p` needs, and the current delivering it is measured **from the
coils themselves** through `ExteriorCoilSet::gradPsi` — `B_z = (1/r)∂_rψ` at
`( R₀, 0 )`, which is *not* on the symmetry axis, so the textbook on-axis loop
formula does not apply and was not used.

**AND THE FIELD IS DERIVED PER RADIUS, WHICH IS WHAT MADE THE
SWEEP CONTIGUOUS.** Take one `a = 0.30` for the whole sweep and exactly one row
— `R = R₀ + a = 1.05` — is on design while the rest are the same coils holding a
plasma of a different size. Measured, that is not a small inconsistency: half the swept radii landed on a spurious branch with `ψ_bnd`
NEGATIVE and the "axis" at `( 0.071, 1.488 )`, which is on `Γ`. Shafranov's
formula takes `a`, so it is given `a = R_limiter − R₀`. `k = 2`, 1333 elements,
four exterior modes:

→ **[M-45](MEASUREMENTS.md#m-45)** — limiter · coil `μ₀I` · coils · Newton · `ψ_ax` · `ψ_bnd` · `\|F\|` on `r = 0` · `ψ_ax` attained at · `Ψ` at the O-point

**FIVE HEALTHY RADII WHERE THERE WERE THREE**, `ψ_ax` attained at
`r = 0.78`–`0.85` where the old fixture attained it at `r = 0.00000`, the
prescribed current delivered to every digit, and both `ψ_ax` and `ψ_bnd`
**monotone and smooth in the limiter radius** — which is itself a reading of the
repaired limiter constraint, since a dof-snapped one is a staircase in exactly
that variable. The coil-free control converges at three of the five and does not
converge at 1.15 or 1.18.

**THE COIL-FREE ROW IS THE CONTROL AND IT IS SHARPER THAN A FAILURE WOULD BE.**
It converges, `|F|` on the axis is exactly zero, and `Ψ` at its O-point reads
**1.0000** — *every health check in the case passes on it*. It is simply a
different equilibrium, and the only thing separating them is **where the axis
is**, so the control asserts on the position.

**WHERE IT GIVES OUT IS THE GEOMETRY AND IS RECORDED RATHER THAN HIDDEN.** From
about `R = 1.20` — 0.8 of `ρ_Γ` — this fixture does not reach a tokamak: `ψ_bnd`
comes out **negative**, the O-point lands on `Γ` at `( 0.071, 1.488 )` carrying
`Ψ ≈ 1.1–1.2`, and `|F|` on `r = 0` is back. It converges, in 8 to 12 steps.
**A different branch, not a worse answer** — the honest fix is a larger `Γ`, so
the sweep stops at 1.18. And `R = 1.05` no longer converges at all, which is
the one thing the per-radius field cost: at `a = 0.30` the plasma fills the
whole guess bump.

**THE CONDUCTORS COST THE FIXTURE NO MESH**, which is FB-7 being used by
something other than its own acceptance the day after it landed;
`meq::ExteriorCoilSet` is what can hold them, and meshing them in would have
meant the gmsh half-disc and a different discretisation from the one that case is
about.

**`theTwoBordersConvergeTogether` KEEPS ITS TABLE AND LOSES A CLAIM.** It closed
§7.13's *"one combination still open"* — both borders reaching a common root —
and it still does; closing on an unphysical equilibrium is still closing. What
its `ψ_ax` and `ψ_bnd` columns are **not** is a machine. **Anything quoting
§7.12b's numbers as an equilibrium is quoting the `r = 0` layer.**

**THE TRANSFERABLE PART**: a red test names a defect, and *where* the defect is
was got wrong for a day. The assertion was right, the diagnosis was right, and
the thing that needed repairing was the fixture rather than the solver. §11.7 of
`FREE-BOUNDARY-PLAN.md` is the full record.

**AND THE FILL DOES REACH THE AXIS — ON INTERMEDIATE ITERATES OF THAT VERY
FIXTURE, AND ON NO CONVERGED ANSWER THE GUARD LETS THROUGH.**
`refreshPlasmaComponent()` refuses to *seed* an axis-touching element and does
not refuse to *reach* one, so an iterate lifting `ψ_h` above `ψ_bnd` there puts
§11.3's `1/r` pole back with `ConfineToPlasma` on. Measured on the machine above
at limiter 1.15, `k = 2`, 1333 elements: `ψ_bnd` is **negative** for the first
three Newton steps, all **84** axis-touching elements are in the plasma's
component, and the assembly puts `|F|` of 3.7e-03 / 1.6e-03 / 5.0e-03 on `r = 0`.
It clears by step 2 to 4 and the converged answer has **0** of the 84 as even a
candidate, with `|F|` on the axis exactly `0.0`. Same shape at 314 and 573
elements.

**THE FILL IS NOT THE DEFECT AND CHANGING IT WOULD BE WRONG TWICE.** At those
iterates its answer is *correct* — `{Ψ > 0}` really is connected and really does
contain `r = 0` when `ψ_bnd < 0` — so a rule refusing to say so would be a fill
lying about the state it was handed. And blocking axis-touching elements outright
switches `p′` off in the one device class §11.6 protects: a levitated dipole and
a magnetic mirror both carry plasma to `r = 0`, both have `gg′ ≡ 0` there, and
the pole never existed.

**WHAT MAKES IT SAFE IS A COUPLING BETWEEN TWO READINGS, AND IT IS ASSERTED
RATHER THAN OBSERVED.** `ψ( 0, z ) = 0` exactly, so the axis is inside the plasma
precisely when `Ψ_axis = −ψ_bnd/span > 0`, which is
`checkAxisSource().axisInsidePlasma` and which the driver refuses on. At a
converged answer `ψ_h( 0, z ) → 0`, so the fill can reach the axis only when that
reading is positive too — measured on both branches, none of 84 at limiter 1.15
with the guard clean, all 84 at 1.20 where `ψ_bnd` goes negative and the guard
refuses on both counts. **So the guard does catch the converged case**, and
`theFillReachesTheAxisOnlyWhereTheAxisGuardRefuses` asserts the *implication*
rather than either column.

**THE GAP THE GUARD CANNOT SEE IS THE DISCRETE LAYER, AND THE MARGIN WIDENS
UNDER REFINEMENT.** `checkAxisSource()` asks at `ψ = 0` deliberately, so it
cannot see a `ψ_h( 0, z )` exceeding a *positive* `ψ_bnd`. That state exists —
an unconfined solve grows a layer reaching `Ψ = 2.7` at `r = 0` while the guard
reads `Ψ_axis = −2.6e-01` — and it is **not reachable from a confined solve**,
because the layer is what the pole builds and the pole is what confinement
removes. Where the layer does exist the **fill** is what separates it, as its own
component: 84 candidates, **0 reached**, `|F|` on the axis exactly `0.0`, which
`theLimiterCaseAlreadyHasMoreThanOneLobe` asserts — and with `elementInPlasma()`
ignored the same state assembles **`|F| = 1.376`** there, so the fill is the
whole of the difference between a bounded load and `μ₀ j_φ` diverging like `1/r`
along the entire axis, not a trim at the edge. Sliding `ψ_bnd` by hand is the
only way found to merge the layer with the plasma, and it takes
`ψ_bnd = 1e-08` against the 7e-04 a real solve carries. `max ψ_h` on `r = 0` runs
5.34e-06 → 2.25e-06 → 6.68e-07 over 314/573/1333 elements, i.e. 7.8e-03 → 9.5e-04
of `ψ_bnd`, **falling at 2.92 against `k+1 = 3`** while `ψ_bnd` does not move. The
two readings separate faster as the mesh is refined, so the coupling is a limit
and not a number.

**AND MACHINE C'S REFUSAL IS A FALSE ONE AFTER ALL — AND THIS SECTION SAID THE
OPPOSITE FOR AN HOUR, ON A LEAKY INSTRUMENT IT HAD ALREADY DOCUMENTED.** The
claim made here was that `C_mast_spherical`'s own reference reaches `R = 0.100`
at `Ψ_N = +0.2863` with the magnetic axis's component, so the guard names a true
property of the target and what saves `freegs4e` is its vessel mask. **The
measurement behind it was a `scipy.ndimage` flood fill on `{ Ψ_N > 0 }`, which
leaks through the X-point** — the identical failure recorded three paragraphs
above for MEQ's own fill, and the one that made machine A's core read 1.13 m²
against its true 0.784.

**ERODE THE THRESHOLD BY A TENTH OF A PER MILLE AND THE LEAK CLOSES.** A level
set touches its own saddle, so two regions meeting at a point are one component
at `ε = 0` and two at any `ε > 0`; regions that genuinely overlap survive the
erosion. On `C_mast_spherical`:

| `ε` | components | axis component reaches | area |
|---|---|---|---|
| 0 | 1 | **R = 0.1000** | 3.196 m² |
| 1e-04 | 2 | **R = 0.3078** | 1.632 m² |
| 1e-02 | 2 | R = 0.3227 | 1.546 m² |

and the reference's own `core_mask` reaches **R = 0.308**. `G_mastu_simple`
behaves identically — R = 0.1000 at `ε = 0`, **0.2781** at 1e-04, against a
`core_mask` reaching 0.278.

**THE ONE-DIMENSIONAL CHECK NEEDS NO FILL AT ALL AND IS THE ONE TO REACH FOR.**
Walk the midplane from the magnetic axis inboard and ask whether `ψ` dips below
`ψ_bnd` on the way. On MAST it falls to **−3.0876e-02 at R = 0.1594**, below
`ψ_bnd` = −2.3596e-02, before rising to −1.3990e-02 at `R = 0.100`. So the
near-axis lobe is a SEPARATE component, connected to the core only through a
saddle, and **connectivity alone separates it — no vessel mask is required.**
`setting psi_bnd IS the confinement`, and the level set plus the component is the
whole of it.

**SO THE GUARD IS CONSULTING THE WRONG COLUMN, WHICH IS WHAT THE DRIVER'S OWN
DIAGNOSTIC SAID ALL ALONG.** It refuses on `axisInsidePlasma` — the pure level
set `−ψ_bnd/span > 0` at `ψ = 0` — while printing *"plasma component kept it:
NO"* on the same line. The assembled source is what has or has not got a pole,
and the component is what decides that, so the refusal belongs on
`axisInPlasmaComponent`. **That strengthens rather than weakens the coupling
asserted above**: `theFillReachesTheAxisOnlyWhereTheAxisGuardRefuses` asserts
*fill reaches ⇒ guard refuses*, one way, and machine C is the converse case the
implication never covered.

**BOTH GUARDS NOW ASK THE SUPPORT, AND MACHINE C CONVERGES.**
`AxisSourceCheck` carries `supportReachesAxis`, `worstOnAxisInSupport` and
`boundedInSupport` beside the level-set readings, and the driver refuses on
those. **The second guard mattered as much as the first**: `bounded` reads
`worstOnAxis`, which is `f()` evaluated pointwise and knows nothing about
elements, so correcting only the topology guard would have moved the same false
refusal one branch down. MAST's shipped configuration goes **exit 1 → exit 0**
with `psi_ax`, `psi_bnd` and the axis unchanged to every digit — the same solve,
a different verdict — and lands **1.737e-02** from its reference. The level-set
reading is kept and REPORTED as a warning, because a run where the two disagree
is one where the fill is the only thing between the load and a `1/r` pole.

→ **[M-109](MEASUREMENTS.md#m-109)**

**AND freegs4e IS NOT USING A WALL ON THAT MACHINE AT ALL**, which is the last
nail in the version of this paragraph that blamed one: `C_mast_spherical.npz`
carries **zero** wall points, so `mask_inside_limiter` is `None` and nothing
masks `Jtor` geometrically there. Connectivity is the whole of what keeps the
reference off the axis.

**A THREADING NOTE, BECAUSE IT COST A WRONG CONCLUSION.** That case sits on the
knife edge CLAUDE.md documents: at `MKL_NUM_THREADS=1` it fails in the bordered
Newton with no finite damped step, and at `OMP = MKL = 16` it converges in 3.
M-103 was taken at 16, so a guard fix verified at 1 reads as no change at all.

**THE EXCEPTION THE FILL CANNOT SEE EITHER IS STILL THERE**: several O-points separated by saddles. Then the
axis-containing component is the right answer only if the OTHER O-points'
components are genuinely disconnected from it, and at the discrete level they
touch at the saddle exactly as the near-axis lobe does. The erosion above is a
diagnostic and not a rule -- it works because these lobes meet at a point, and
`PlasmaConnectivity` would need its own saddle-aware cut to be sure. `freegs4e`
carries `geom_inside_mask`, a half-plane through the active X-point, for exactly
this.

**AND THE ATTRIBUTE IS BUILT, WHICH IS WHAT COVERS THAT EXCEPTION.**
`[source] ExcludeAttributes` names mesh ELEMENT attributes that can never be
plasma; `meq::PlasmaComponent::holds()` tests them ahead of its own
constant-true shortcut, so the source is off out there with or without
`ConfineToPlasma`, and `refreshPlasmaComponent()` treats such an element as
neither a candidate nor a seed nor traversable — which is what stops it BRIDGING
two lobes that are otherwise separate, since the straddling band is shared out
by a watershed over the candidates.

**AN ATTRIBUTE AND NOT A POLYGON, AND THAT IS THE ECONOMY**: the support is
re-decided on every residual evaluation and a vessel does not move, so a
point-in-polygon per element per evaluation would be paying repeatedly for an
answer that cannot change. `tools/mesh/halfdisc.py --vessel` fragments a closed
polygon into the geometry and tags everything outside it — and not a conductor —
as attribute **30**, by the same centroid-AND-AREA test `--limiter` already
needs, since the outer region's centroid can land inside a convex vessel.
`[mesh.generate] Vessel` is the key, and `meq --mesh-command` emits it at full
precision beside the conductors.

→ **[M-110](MEASUREMENTS.md#m-110)** — the whole path · what it excludes · what
it leaves alone

**INERT WHERE IT SHOULD BE, AND THAT IS THE PROPERTY THAT WAS MEASURED.** On
machine A with its real vessel, 1373 of 2105 elements are excluded and `psi_ax`
comes out **identical to every printed digit** to the same run with the
exclusion removed. It buys the fill work rather than answers, until the
configuration is one connectivity cannot settle. An attribute the mesh does not
carry is **refused**, exit 1, listing what the mesh does have: a
silently-empty exclusion is the worst outcome available, since the run would
converge, report nothing unusual, and describe a machine with a current channel
behind its own wall.

**AND FRAGMENTING A VESSEL PERTURBS THE MESH.** On machine A the
vessel-fragmented mesh fails in the bordered Newton from the shipped `PsiAxis`
— **and so does the same mesh with the exclusion removed**, which is what says
it is the mesh and not the feature; regenerating without the vessel converges.
From `PsiAxis` = 5.0e-02, the same basin by M-107, the vessel mesh converges at
1.917e-03 against the committed mesh's 2.762e-03. That machine sits near a fold
with two roots, and a three per cent change in the mesh is what such a case is
fragile to.

**Most of the candidate restriction already existed** and this widens it:
`[mesh.generate] PlasmaRMin/RMax/ZMin/ZMax` refines a box and the support
reports candidates against it — 815 of **2148** on machine A. MAST's
configuration sets `PlasmaRMin = 0.0000`, so the axis is a candidate at all,
against machine A's 0.5925; a vessel would have bounded it independently of the
guard.

### `findAxis()` returns a COIL's O-point on a machine with meshed conductors

**A coil's O-point is a genuine critical point of `ψ`, so a general
critical-point finder is right to return it.** `meq::CriticalPointFinder` has no
conductor exclusion and cannot have one — it is handed a field, not a machine.
`apps/meq.cpp` carries a three-tier exclusion on the plasma fill's seed for
exactly this reason, which is §10.7's defect 2, and **nothing carries one on
`findAxis()`**.

Measured on `examples/diverted-tokamak.toml`'s mesh through the library fixture:
`findAxis()` with no argument **refuses**, reporting *"no unique interior
extremum … 3 maxima, 2 minima and 2 saddles"*, which is the honest answer.
Forcing it with `AxisSense::Maximum` returns `( 1.006240, -1.099327 )` — inside
coil P1L, `r ∈ [ 0.95, 1.05 ]`, `z ∈ [ -1.15, -1.05 ]` — on every solve of two
different equilibria at two resolutions. It is stable, reproducible and is not
the magnetic axis.

**So an axis position obtained that way on a machine with `[[coils]]` in the
mesh is about a conductor until something says otherwise**, and a displacement
quoted in metres from such a call needs checking against the coil list before it
is believed. → **[M-132](MEASUREMENTS.md#m-132)**. The refusal is the feature;
the forcing argument is what removes it.

### The cold failures are three different things, and the merit's balance is machine-dependent

**NOT ONE OF THE SIX COLD `exit 2` CASES HAS A NON-FINITE DIRECTION**, which is
the condition `BORDERED-GLOBALISATION-PLAN.md` §0.1 names as the one that would
make the whole globalisation campaign moot. They are three distinct failures:
two machines lose a SINGLE step of twenty-seven; one is chronic, exhausting the
halving ladder on 55 of 288; and three are not step-length failures at all but
**X-point EXCURSIONS** — TCV's null travels 1.6 m in Z from its seed and E's
flips sign in one step.

→ **[M-117](MEASUREMENTS.md#m-117)** — the classification · the traces ·
→ **[M-116](MEASUREMENTS.md#m-116)** — `PicardSweeps`, which fixes none of them

**AND THE SAME TRACES EXPLAIN WHY NO MERIT WEIGHT EVER TRANSFERRED BETWEEN
MACHINES.** The border constraints' share of `augmentedNorm` differs by three
orders of magnitude across two cases of one benchmark: on MAST `g·axis` alone is
**0.75** of a merit of 0.843, and on TCV the whole border is **~1e-03** against a
field residual of 2.953. So MAST is the end where the border already dominates —
raising it cannot do anything, which is exactly what `BorderMeritWeight`'s
upward saturation measured — and TCV is the end where nothing restrains the
X-point at all. M-113's inversion, M-114's local optimum and that saturation are
one fact seen from opposite ends of a range, not three findings.

**`PicardSweeps` IS THE CHEAPEST UNTRIED THING AND IT IS A TRAP.** It had never
been set by anything in the tree. It fixes none of the six and turns **three of
them from `exit 2` into `exit 0` with a wrong answer and no warning** — C shaped
out by a factor of 3.9, both E rows by 23%. A failure is actionable; a green run
reporting a different machine is not.

**AND E's DEFECT IS ROOT SELECTION, NOT GLOBALISATION.** The Picard pre-stage
and M-106's exact-seed run reach 1.0308e-01 and 1.0252e-01 — 0.5% apart, both
23% from the reference. Two routes sharing nothing find the same wrong branch,
so no step-length management is the repair there.

### The bordered path has a globalisation, and it is Picard in the field alone

`Globalisation::BorderedPicardThenNewton` is the second rung of the driver's
reactive ladder on a free-boundary run, and until it existed there was no rung
at all: every other globalisation MEQ has either drives a residual of its own
(the KINSOL ones) or puts the potential block on the **linear** form and builds
no Jacobian for a border to be eliminated against (`AndersonPicard`,
`PicardOnly`, and `PicardThenNewton` which is those two in sequence).

**IT CHANGES ONE BLOCK AND NOTHING ELSE.** The border, its elimination, the
Armijo loop and the augmented norm are untouched; `meq::FieldLinearisation`
chooses only what the semi-linear term contributes to the **derivative**:

```
Newton   field block = A_lin - ( 1/r )( dF/dpsi ) M
Picard   field block = A_lin
```

**And that is why there is a Jacobian in a method that is not supposed to have
one.** Undamped, from the current iterate, `A_lin d = -( A_lin u^k - b( u^k ) )`
gives `u^{k+1} = A_lin^-1 b( u^k )` — the classical Picard map written as a
defect correction. The increment form is what lets the bordered machinery be
reused verbatim, and under hybridization applying `A_lin^-1` is a `ComputeH`
plus a trace solve either way, so it costs nothing extra.

**THE RESIDUAL IS NEVER FROZEN, WHICH IS A DEPARTURE FROM
`BORDERED-GLOBALISATION-PLAN.md` §3.2 AND IS THE WHOLE OF THE DESIGN.** At the
iterate `F( r, z, psi_h )` and `F( r, z, psi^k )` are the same numbers, so
freezing the residual changes no step. It changes only what is evaluated AWAY
from the iterate, and both places it reaches are places it does harm: the Armijo
trials, where a frozen `F` makes the field part of the merit affine in the field
so the line search stops seeing the excursion it exists to catch; and the
differenced border columns, which would then be derivatives of a source the
residual is not using. §3.3's "one-line consistency requirement" does not arise
because the inconsistency does not arise. `AssembleElementVector` is untouched
and the change is an early `return` in `AssembleElementGrad`.

**IT CANNOT CHANGE WHICH EQUILIBRIUM IS REPORTED, AND THAT IS ASSERTED RATHER
THAN ARGUED.** The fixed points of `Phi( x ) = x - alpha M( x )^-1 G( x )` are
the zeros of `G` whatever non-singular `M` is, and only the field block of `M`
moves — but `by construction` is exactly the class of claim this project does
not accept on its own, and M-26 has three solve routes reaching discrete
solutions 9.4% apart. `theBorderedPicardReachesTheSameEquilibriumAsTheBordered
Newton` measures it on the diverted machine: **`psi_ax` 8.26600363e-02 and
`psi_bnd` 3.23793176e-02 from both routes, every printed digit**, 14 Newton
iterations against 15 of which 8 ran under Picard. The case also asserts that
`borderedPicardIterations()` is non-zero, because an implementation that wired
the enum through and quietly assembled the Newton Jacobian anyway would pass
every agreement assertion perfectly.

What changes is the iteration matrix, `I - M^-1 J = M^-1 diag( K, 0 )` with `K`
the reaction, so convergence is linear at a rate set by the reaction against the
elliptic operator. M-36 already records which side of that trade this is: **a
robustness route, not a faster one.**

**WHAT IT WAS EXPECTED TO BE WORTH WAS WRITTEN DOWN BEFORE THE MEASUREMENT, AND
THE MEASUREMENT FALSIFIED IT.** The prediction: the fixed point does not move, so
it helps only by having a different basin, and `A_lin` is unconditionally
invertible where `A_lin - K` is what goes near-singular — which against M-117's
classification made **B and B-shaped** (one fatal step) and **C-shaped**
(chronic) the addressable ones and D/E/E-shaped the X-point excursions it could
not reach. **Three of six, recorded as the maximum rather than the hope.**

**IT CLOSES NONE OF THE SIX, AND THE ONE CASE IT HELPS IS NOT ONE OF THE TWO THE
ARGUMENT LEANED ON.** → **[M-119](MEASUREMENTS.md#m-119)**. B and B-shaped are
made *worse*, and both arms there plateau at about half the initial residual and
wander — which is neither a fatal step nor anything a linearisation reaches.
Four of six reach a worse minimum than plain Newton, E's residual growing 6.7×
where Newton's merely stalls, because a Picard direction is longer and less well
aimed and so walks the null further on exactly the excursion cases.

**AND THE ONE ROW THAT PAYS, PAYS AS AN INSTRUMENT RATHER THAN AS A REPAIR.**
`c-mast-shaped` goes from 163 iterations diverging to 3.0e+06 down to 3.19e-01
in 35, and then throws *"the bordered Jacobian is singular in
( psi_ax, psi_bnd, a )"* — a guard standing since FB-3 that the Newton arm never
reaches. With the field block replaced by something that cannot be singular,
what is left singular is **the border's Schur complement**. So the rung's real
use is as a way of ASKING whether a failure is the field's, and on five of six
the answer is no. That is `BORDERED-GLOBALISATION-PLAN.md` §0.1's first row
reached by measurement, and it points at §6.1's degenerate axis row and at
M-115 rather than at anything else in that plan.

**The handoff tolerance is also measured useless here.** It is 1e-3 relative and
no case gets below 2.5e-01, so on all six the phase change is the budget or the
line-search rescue and never the tolerance. A tolerance is the right control for
a case that is merely slow; none of these is.

Two smaller departures from the plan, both recorded so the plan is not read as
the built thing. BG-3 proposes re-purposing `Globalisation::PicardOnly`; a new
value is used instead, because `PicardOnly` makes `usesNonlinearForms()` false
and so moves the potential block onto the linear form — the exact route §3.2
spends four reasons ruling out. And BG-4 proposes two calls to `solve()`; both
phases are iterations of **one** loop inside `solveWithNormalisation()`, because
a second call re-enters `prepare()`, which under NPC seeds the potential and the
trace and leaves the **flux** block at zero. One loop also keeps `gamma` frozen
at the first iterate across both phases, so the printed history is one merit end
to end.

### The axis row's envelope term is written, correct, and off by default

The row drops `grad psi_h( x* ) . d( x* )/du` by the envelope theorem, which is
exact only if `grad psi_h( x* ) = 0`. That holds for the CONTINUOUS fields and
not for the discrete ones — `r q_h - grad_bar psi_h` is the local lifting of the
trace jump — and `cornerEntry()`'s XP-3 arm said so about the same identity all
along. → **[M-120](MEASUREMENTS.md#m-120)** measures the dropped term at 0.03 to
0.16 per cent of the row.

`AxisRow::WithEnvelope` assembles it, and → **[M-121](MEASUREMENTS.md#m-121)**
is why `PositionDropped` is the default. It is correctly signed — XP-3's answer
does not move to any printed digit while its sweeps go 4 to 3 and its final
residual 1.702e-13 to 5.551e-16 — and it costs `FreeBoundaryCoupling`'s physical
fixture one of its five limiter radii. The endgame improves and the approach
gets worse, which is the shape of a sharper derivative of a function that is not
differentiable: M-120's other finding is that `C_Ax` JUMPS when the axis crosses
a face. **Make the constraint continuous first**, and then the term is either
unnecessary (define `x*` as a critical point of `psi*` and the envelope theorem
is exact by construction) or safe.

## The geometry: meshing a half-disc that reaches the axis

Free boundary is the one campaign here whose geometry
`mfem::Mesh::MakeCartesian2D` cannot build, so it is the one that needs a mesher.
Both sections below are `tools/mesh/halfdisc.py`: what it produces and the
defects found producing it, then why gmsh and why it is not linked into
`meq_core`.

## Meshing the half-disc, and the boundary tag nothing would have reported

**`tools/mesh/halfdisc.py` is the generator free boundary needs**, and MFEM's
built-in mesher cannot make what it makes: a **semicircle centred on the axis,
reaching `r = 0` exactly**, with the coil rectangles fragmented in as their own
subdomains, coarse over the vacuum and refined where the plasma is. gmsh 4.14 is
Debian-packaged here and MFEM reads its format natively, so the coupling is a
**file** and nothing else.

**`r = 0` IS THE REQUIREMENT AND IT IS EXACT.** Measured through `mfem::Mesh`
rather than through gmsh — `r ∈ [0, 1.5]`, **21 vertices at exactly `r == 0.0`**,
20 boundary faces on the axis at `|r| = 0.000e+00`, and Γ's vertices on the true
circle to **4.441e-16**. At geometric order 2 the mid-edge nodes land on the arc
to the same 4.441e-16 and the axis nodes stay at exactly `0.0`, which is why the
geometry is a **disc minus a half-plane** rather than an arc plus two lines: OCC
owns the curvature, so gmsh has the real circle to place them on. `--check`
asserts `r == 0.0` **without a tolerance**, because a domain stopping at
`r = 0.05` is not a slightly worse semicircle — the Gegenbauer basis of
`FREE-BOUNDARY-PLAN.md` §3 does not span its exterior, and nothing downstream
would say so.

**AND THE FIRST VERSION TAGGED TWELVE INTERIOR EDGES AS Γ.** Γ and the axis were
selected by asking each 1-D entity where its centre of mass sat; after
`occ.fragment` the model also holds **each conductor's own outline**, every one
of which is "not on the axis", so they joined Γ. A 1-D physical group is written
to the `.msh`, and **MFEM's reader turns 1-D elements into BOUNDARY elements
whether or not they are topologically on the boundary**. Measured: 64 boundary
faces, of which 32 on the arc, 20 on the axis and **12 on the coil rectangles**.
Attribute 1 carries the transferred exterior datum, so the run would have
imposed Γ's Dirichlet condition on twelve edges buried inside the conductors —
converging at full order to the wrong equilibrium.
`getBoundary( surfaces, combined = True )` is the fix and is the right question:
combined over every surface, shared edges cancel and the outer boundary
survives. 32 + 20 + **0 stray** now.

**THE INSTRUMENT WAS WRONG TOO, WHICH MAKES IT FIVE TIMES IN THIS TREE.**
`gmsh.model.mesh.getNodes()` returns coordinates ordered to match its **tag
list**, not sorted by tag — the tags come out grouped by owning entity and are
neither sorted nor contiguous. Indexing `coords[3*(tag-1)]` therefore reads a
*different node*, and the first `--check` did: it reported the axis group at
`r = 1.4928` and a coil spanning most of the disc, which looks like catastrophic
meshing and was the checker. Same shape as the cone tiling, the MXH distance,
the transmission row's shared rule and the FB-1b single-mesh check.

**COARSE DISC, REFINED PLASMA, AND THE SAVING GROWS WITH `ρ`.** Γ must sit far
enough out for the exterior expansion to converge and the plasma occupies little
of what that encloses, so a uniform mesh at the plasma's resolution spends most
of its elements on sourceless vacuum. Against a uniform mesh at `h = 0.035`:
`ρ = 1.5` reads 2,965 against 7,672 (**2.6×**) and `ρ = 3.0` reads 3,368 against
30,392 (**9.0×**) — moving Γ out costs the graded mesh 14% more elements and the
uniform one four times as many. That is the direction the coupling pushes.

**`[mesh] File` USED TO LOSE THE GRIDDED OUTPUT, AFTER SOLVING.** `RMin`…`ZMax`
describe a box that was not built, so the extent came out empty and
`meq::GridSampler` refused — the answer was lost at the output stage. **Not a
corner case**: the half-disc cannot come from `MakeCartesian2D`, so a
free-boundary run *always* reads its mesh from a file. The driver takes the
mesh's own bounding box now, which is also the right answer.
`theDriverTakesItsGridFromAMeshItDidNotBuild` pins it.

**WHY PYTHON, AND WHY GMSH IS NOT LINKED INTO `meq_core`.** gmsh ships a C++ API
and a Python module, and **the Python module is `ctypes` over the same
`libgmsh.so`**, generated from the same `api/gen.py`: identical function set,
identical meshes, the same code doing the work. Linking would buy no capability
— MFEM reads the `.msh` natively, so the file is the whole interface — and would
cost `libgmsh-dev` as a build dependency of a tree whose CI cannot build the
solver at all. The one thing that would change the answer is re-meshing *during*
a solve; MEQ's adaptive loop **refines** through `meq::AdaptiveDomain` instead,
so that is not on the plan, and the script transliterates if it arrives.

## Meshing beyond `MakeCartesian2D`: gmsh, and it is already here

**FB-6 needs a geometry `mfem::Mesh::MakeCartesian2D` cannot make** — a
semicircle centred on the axis, because §3's exterior expansion is valid only
there; reaching `r = 0` exactly, because FB-A says so; with the coil rectangles
as their own subdomains. Everything MEQ has meshed until now was a rectangle of
triangles with curved geometry reached by the extension on top.

**The answer is gmsh and it needed no procurement**: MFEM reads its format
natively (`Mesh::ReadGmshMesh`, no converter), and **gmsh 4.14 is already
installed on this machine** with its Python API. `tools/mesh/halfdisc.py` is the
generator, and the round trip is verified rather than assumed:

```
  elements 781  vertices 424        domain attributes: 1 10 11
  attr 10:   12 elements   r [0.9000, 1.1500]  z [-0.8600, -0.7400]
  attr 11:   12 elements   r [0.9000, 1.1500]  z [ 0.7400,  0.8600]
  boundary attributes: 1 (arc) 2 (axis)
  minimum r over the mesh: 0.000e+00
```

**Its physical groups arrive as MFEM attributes**, which is the whole reason it
is the right tool: one attribute per coil so a source can be restricted to it,
and *separate* boundary attributes for the arc and the axis — which are not
interchangeable, `Γ` carrying the transferred exterior datum and the axis a
plain `ψ = 0`. They are told apart by geometry rather than by tag order, which
OCC does not promise.

**AND THE COILS ARE MESHED TO RATHER THAN CUT**, which is §7.9's finding applied
rather than a convenience: a conductor's geometry is **prescribed input** and
does not move with the solution, so it can always be aligned to — and aligning
it took the measured rates from 1.33/1.27/1.09 to **1.99/2.88/3.01**. Only the
plasma support genuinely moves, and FB-4 measures what that costs.
`occ.fragment` is what makes the coil edges mesh edges.

**What was NOT chosen, and why the survey was short.** Triangle is 2D-only and
has no subdomain tagging worth the name; CGAL's `Mesh_2` is a large dependency
for a 2D job; MMG is a re-mesher rather than a generator and would sit *after*
this rather than instead of it — worth revisiting if FB-5's adaptivity ever
wants to move the geometry between cycles. Gmsh reads a `.geo` script, drives
from Python or C++, and is what MFEM's own examples assume.

## FB-6's test problem: no reproducible ITER case, and the profile that is

**The record is `FREE-BOUNDARY-PLAN.md` §7.11.**

**CEDRES++'s ITER case cannot be reproduced**: no coil currents, no wall, and a
reference that is their own fine mesh printed as a figure. Serino et al.'s 15 MA
ITER baseline is an *inverse* solve seeded from a proprietary discharge
(ABT4ZL). Nothing else in `refs/` carries a free-boundary case with data.

**But CEDRES++'s profile model (2.11) and `freegs4e`'s `ConstrainBetapIp` are
the same family term for term** — `α`↔`alpha_m`, `γ`↔`alpha_n`, `β`↔`beta0`,
`r₀`↔`Raxis`, `λ`↔`L` — so the published ITER profile is runnable against an
independent code. ITER: `α = 2`, `γ = 1.395`, `β = 0.5978`, `r₀ = 6.2 m`,
`I_p = 15.10 MA`.

**`α = 2` AND NOT THE `0.5978` CEDRES++ PRINTS.** Their §4.1 prints `β`'s value
twice; Serino et al., citing the same ITER coefficients, give the peakage
exponent as 2 under the name `δ` — their `α` being CEDRES++'s scaling `λ`, a
symbol clash rather than a disagreement. Both readings taken off **rendered
pages**, which is the standing rule for this pair.

**AND THE PUBLISHED ITER PROFILE IS A `k = 1` CASE.** The vanishing order at the
plasma edge is `j = γ` exactly, so `γ = 1.395` caps `ψ*` at `j + 2.5 = 3.895` and
`k+2` survives only at `k = 1`. Worse, `γ < 2` makes `∂F/∂ψ ~ (1−Ψ)^{γ−2}`
**unbounded** at the edge, where §7.10's fixture merely jumps. **Open FB-6 at
`γ = 2` or 3** — one number in the same family — and bring `γ = 1.395` in
afterwards as the published case, expecting `k ≥ 2` not to hold `k+2` on it.

**`freegs4e` has no ITER**: `TestTokamak`, `DIIID`, `MAST`, `MAST_sym`, `TCV`,
`MASTU_simple`, `MASTU`. Any serves, and the fixed-boundary rehearsal already
drives seven configurations across them.

## MEQ against freegs4e, and root selection is the whole difficulty

### STANDING RULE: COMPARE IN THE SAME MODE MEQ SOLVES IN

**Every comparison against `freegs4e` is to be run with `freegs4e` in the SAME
mode as MEQ — forward against forward — and the benchmark as it stands does
not do this.**

`fgsref.py` drives every one of the seven references through
`SyncConstrain( xpoints = ..., isoflux = ... )`, which is `freegs4e`'s
**INVERSE** solve: the coil currents are re-solved on every Picard pass to hit
X-point and isoflux targets. **The currents are an OUTPUT of the reference and
an INPUT to MEQ**, which solves FORWARD from them. So the two codes are not
being asked the same question, and the difference between them is not a
measurement of either solver.

**WHAT THE MISMATCH COSTS, AND IT IS THE SIZE OF EVERY DISAGREEMENT IN THE
BENCHMARK.** An inverse solve is pinned to its targets, so it reports very
nearly the same equilibrium however its conductors are modelled; a forward
solve is not, so every difference in conductor model, Green's function or
discretisation moves the answer. Measured on machine A, both arms at `k = 2`:

| | `psi_axis` |
|---|---|
| reference, filaments | 8.2717514448e-02 |
| reference, `ShapedCoil` | 8.2716849781e-02 |
| **the two references against each other** | **8e-06** |
| MEQ forward, filament currents | 8.068175e-02 |
| MEQ forward, shaped currents | 8.073704e-02 |
| **MEQ against either** | **2.4e-02** |

**The two references differ in conductor model and agree to 8e-06 BECAUSE the
inverse solve re-tunes the currents to the same targets.** MEQ is 2.4e-02 from
both. A comparison built this way measures the forward sensitivity of the
equilibrium to a current set, which is a real and interesting quantity and is
not what the benchmark claims to report.

**AND IT IS WHY THE PER-MACHINE AGREEMENT SCATTERS BY FOUR ORDERS OF
MAGNITUDE**, which nothing else explains — the machines with few conductors
close to a small plasma are the ones whose forward solve is most sensitive to
the current set, and DIII-D, with eighteen conductors far from a large one, is
the least:

| machine | MEQ `psi_ax` against its reference |
|---|---|
| DIII-D, shaped | **3.4e-06** |
| DIII-D, filament | 3.5e-04 |
| A testtokamak | 2.5e-02 |
| TCV, shaped | 3.5e-02 |
| MAST-U | 4.9e-02 |

**SEVEN OTHER CAUSES WERE ELIMINATED FIRST, ON MACHINE A, AND THE RULE IS WHAT
SURVIVED.** MEQ's own `h` and `p` refinement (converged, 8.056e-02 across a 4x
element refinement and a degree change); `freegs4e`'s grid (converged —
`psi_axis` moves 1.1e-05 relative from 129^2 to 513^2, so M-88's first-order
reading does not apply to this quantity on this case); the conductor model
(both arms, above); the exterior mode truncation (converged by `Modes = 20`);
the `ff'` amplitude (a +-1.7e-02 change moves `psi_ax` by only -+2.9e-03, so
closing 2.6e-02 would need 15 per cent); **the profile conversion entirely** —
rebuilding MEQ's tables from the reference's OWN solved `pprime` / `ffprime`
arrays rather than the analytic shape moves `psi_ax` by 1e-05 relative; and the
plasma current, which agrees to **9e-11**.

**What is left is the axis POSITION**: MEQ puts it 2.5 cm inboard and 1.5 cm
high of the reference's, with `I_p` identical — the signature of the right
total current distributed differently, which is exactly what a forward solve
does with an inverse solve's currents.

**HOW TO FIX THE COMPARISON.** Run `freegs4e` forward — no `constrain` — from
the reference's own converged currents, and difference THAT against MEQ.
`tools/freegs4e-benchmark/forward.py` is that run, and it reproduces the
inverse reference on machine A to **6.3e-08**.

**AND THE MODE MISMATCH IS NOT WHAT THE PER-MACHINE DISAGREEMENT IS. MEASURED,
AND THIS FILE SAID OTHERWISE FOR HALF A DAY.** The sentence that stood here
told a reader to treat M-87, M-96 and M-103's accuracy columns as "an upper
bound carrying an unmeasured mode mismatch". It is measured now, on A:

| | `psi_axis` |
|---|---|
| `freegs4e` INVERSE, 129² | 8.2717514448e-02 |
| `freegs4e` FORWARD, 129² | 8.2717509068e-02 |
| `freegs4e` FORWARD, 257² | 8.2726591774e-02 |
| MEQ, `k = 3`, one refinement | **8.0567700000e-02** |

**The two modes agree with each other to 6.5e-08 at the same grid, and refining
the forward solve moves it by 1.1e-04** — three orders short of the 2.6e-02 MEQ
sits away from both. So the mode was worth less than 1e-04 here and the
disagreement is something else.

**THE RULE STILL STANDS AND ITS REASON IS UNCHANGED**: an inverse solve's
currents are an output, so a comparison built on them measures the forward
sensitivity of an equilibrium to a current set. That is a methodological defect
whatever its size, and the tell that provoked it — two references differing in
conductor model and agreeing to 8e-06 — is still exactly what an inverse solve
re-tuning to fixed targets looks like. What is falsified is only that the
mismatch EXPLAINS the disagreement.

**NINE CAUSES WERE ELIMINATED ON MACHINE A AND THE TENTH IS THE ANSWER: THE
PRESCRIBED-CURRENT ROW, AT A TANGENCY.** The eliminations stand and are worth
keeping, because each one is a thing nobody now has to re-test — MEQ's `h` and
`p` refinement (`psi_ax` 8.056e-02 across a 4× element refinement AND a degree
change); `freegs4e`'s grid in BOTH modes; the conductor model, now measured
DIRECTLY rather than by the two references agreeing — spreading the reference's
two point filaments into the 10 cm squares MEQ meshes moves `psi` by **1.5e-06
at the axis and 4.7e-06 at the active X-point**, and the X-point itself by
5e-05 m, the 22 per cent maximum difference being entirely the filament's own
log singularity sitting on top of the coil; the exterior mode truncation
(converged by `Modes = 20`, 4e-04); the exterior RADIUS, which is the test the
mode count cannot do — a map with a wrong coefficient converges spectrally to
the wrong answer — `psi_ax` reading 8.059e-02, 8.101e-02, 8.065e-02, 8.060e-02
at `rho_Gamma` = 2.05, 2.2, 2.4, 2.55, every one of them 2.5e-02 from the
reference; the `ff'` amplitude; **the profile conversion, exactly** — rebuilt on
`freegs4e`'s OWN 256 knots so that `ff'( Ψ )` agrees between the codes to
**5e-16**, since a C² cubic spline IS the cubic Hermite through its own values
and derivatives and `meq::SplineProfile` is a Hermite cubic; the plasma support,
which carries the core between the X-points and NOT the private flux region; and
the source arithmetic itself, since MEQ's tables, span, scale and `mu0`
reassembled from the written field reproduce MEQ's own reported `I_p` to
**2.3e-06**.

**AND THE 2.46e-02 ITSELF IS CONDITIONAL — READ M-107 BEFORE ANY OF WHAT
FOLLOWS.** That number was taken on a configuration whose `[source] PsiAxis` had
been set to the reference's own `psi_axis` while chasing something else, and
`PsiAxis` is a STARTING VALUE that selects which of two roots is reported. From
the SHIPPED value — the cold guess's peak — the same machine, the byte-identical
seed file and the same everything else gives **2.76e-03**. The bad root's basin
is the window ( 7.2e-02, 8.6e-02 ) and the true answer 8.271751e-02 is INSIDE
it, so starting at the truth is a factor of nine worse than starting 25 per cent
away. Everything below is a correct account of the branch reached from there; it
is not an account of what the shipped configurations do.

**WHAT REMAINED WAS ONE ROW, AND REMOVING IT CLOSES THE GAP BY A FACTOR OF
SEVENTEEN.** `[source] PlasmaCurrent` makes the profile amplitude an unknown and
prescribes the current instead. Drop it — the amplitude is then fixed at one,
which is exactly the amplitude the tables were built at — and MEQ lands
**1.49e-03** from the reference in `psi_ax` at `k = 3`, 3.9e-04 in `psi_bnd` and
**3.6e-04 m** in the X-point, improving from `k = 2` to `k = 3` where the
constrained arm does not improve at all.

→ **[M-105](MEASUREMENTS.md#m-105)** — the two arms at two resolutions ·
`I_p( λ )` and its maximum · what asking for less current does

**THE MECHANISM IS A FOLD AND IT IS PHYSICS.** `I_p( λ )` on this machine PEAKS
at 2.0013e+05 A near λ = 0.967, because raising the profile amplitude raises the
current density and SHRINKS the plasma — core area falls 0.784 → 0.655 m² from
λ = 0.967 to 1.10 — and past the peak the second effect wins. That is the
equilibrium current limit. The prescribed 2.0000000e+05 sits **0.065 per cent
below that maximum**, so the row `I_p( λ ) = I_target` is a near TANGENCY: two
roots nearly on top of each other and `∂I_p/∂λ ≈ 0` between them. The one
amplitude in the sweep that fails to converge outright is 0.97, at the top of
the curve.

**AND THE REFERENCE'S OWN OPERATING POINT IS NOT A ROOT.** `freegs4e` reports
`Ip_logic` L = 1.000000 — λ = 1 carrying exactly 2.0e5 — where MEQ at λ = 1
carries 1.997272e+05 at `k = 2` and 1.998794e+05 at `k = 3`. **That deficit is
MEQ's own discretisation of the current integral and it halves under
refinement**, but while it stands the target is unreachable on the branch the
reference is on and the constrained Newton leaves for the tangency. Ask for
1.99e+05 instead and MEQ comes back to **2.42e-03**; ask for 1.98e+05 and
**1.17e-03**, with the scale flipping from 0.963 to 1.013 as the solver crosses
the peak.

**SO THE PROFILE SCALE IS A DIAGNOSTIC MEQ ALREADY PRINTS AND NOBODY READS.** On
a case whose profiles came from a converged reference, λ far from one says the
current row is working hard; here it says the row is at a tangency. The control
is DIII-D, seeded identically from its own reference: **λ = 0.9993353**, and its
constrained answer is 3.48e-04 from its reference. *That* is what M-103's
"DIII-D is the only machine that closes" is measuring. **The obvious next move
is for MEQ to report `∂I_p/∂λ`, which the bordered Newton already assembles as
part of its own Jacobian, and to say so when it is small** — the number is free
and the failure it names is otherwise invisible.

**AND THE SCREEN WORKS ACROSS THE SET, WHICH IS THE PART THAT GENERALISES.**
With every machine seeded from its own reference at `n = 128` and run on both
arms, **six of fourteen close at 3.5e-02 or better against two in M-103** — and
the amplitude separates the answers with nothing in between: every row with
∣λ − 1∣ < 0.007 lands between **3.4e-06 and 3.5e-02**, and both B rows at 1.2302
and both E rows at 0.4142 are **2.2e-01 to 6.9e-01** out. The axis guard
independently refuses both B rows, so two screens that share no reasoning agree.

→ **[M-106](MEASUREMENTS.md#m-106)** — the set, both arms, one seed each ·
→ **[M-107](MEASUREMENTS.md#m-107)** — which root `PsiAxis` selects

**THE COLD START IS THE LARGER HALF OF M-103.** Seven of its rows exit 2 with no
finite damped step; most of those converge once the seed is the reference's own
reconstruction, with nothing in the solver changed. **And the seed's RESOLUTION
is a variable**: `mkexactguess.py` defaults to `n = 32`, measured on
`examples/limited-tokamak.toml`, and at `n = 32` machine A fails on BOTH arms
where the identical construction at `n = 128` converges on both. A sweep built on
the default measures the seed.

**DROPPING THE CURRENT ROW IS NOT A GENERAL CURE, AND MACHINE A'S FOLD IS A'S
ALONE.** It is the whole answer on **C shaped** — prescribed has no finite damped
step, fixed reaches 8.335e-04 — and the opposite on A, F and G, where the fixed
arm is worse or stops converging: F shaped goes from 3.4e-06 to 2.5e-03.

**THE VERTICAL-INSTABILITY READING IS CONFIRMED AND IT IS TEXTBOOK.** Seeded at
its own converged equilibrium with the currents frozen and no constraint, the
DIII-D forward solve holds for forty passes and then grows by a factor of ten
every twenty: `Zaxis` −1.55e-04 → −0.555 while `Raxis` moves 1.4e-05 m, with
`rel` growing in LOCKSTEP at 1.122 per iteration. One unstable eigenvalue, and
it is a rigid vertical displacement. Machine A, run identically, reports a
`Zaxis` span of **0.000000**.

→ **[M-104](MEASUREMENTS.md#m-104)** — the trace · the gain scan · the whole
set pinned

**AND ONE SCALAR OF FEEDBACK COSTING LESS THAN A MILLIAMP FIXES IT**, which is
what gives the standing rule a forward reference to compare against at all. One
antisymmetric current combination driven by the axis's own height — no flux
targeted anywhere, nothing re-optimised — converges DIII-D in 12 passes to
**9.8e-08** of the inverse reference at a current increment of **8.1e-04 A**
against coil currents of 1.75e+05 A. Eight of twelve cases close that way. Read
`max |dI|` before the agreement: it is what says the pinned equilibrium IS the
unpinned one. The four that do not close are over-driven rather than unpinnable,
so the harness calibrates the gain from the machine's own response by default.

**WHAT THIS DOES NOT SETTLE**: whether MEQ's own mode sits between forward and
inverse. MEQ pins `psi_bnd` at a LOCATED X-point, and locating a saddle is a
determination rather than a constraint — `freegs4e` does the same thing through
`find_critical` — so the two look equivalent, and nothing here has tested it.

**Two of `forward.py`'s rows are its own defect and must not be read**: TCV and
MAST-U rebuild 24 and 26 conductors against 20 and 14 circuits, because
solenoids expand into windings in the `.npz` and a per-row rebuild is not
`freegs4e`'s `Solenoid`. Their flux is wrong before any iteration — 50× the
field range on TCV — and M-96 already names those as the solenoid-carrying
machines.


**`tools/freegs4e-benchmark/`.** The first check of MEQ against a
code that shares the equation and essentially no code — `../freegs4e`, a FreeGS
derivative: free boundary by von Hagenow Green's functions, 2nd/4th-order finite
differences on a uniform `(R,Z)` grid, Picard with adaptive blending, in Python.
`docs/validation.rst` is the user-facing account; this is the record.

**WHY IT IS WORTH MORE THAN A FINER MESH.** Everything else in this file checks
MEQ against closed forms, manufactured solutions or its own refinement, and all
of those share MEQ's conventions. This file records three separate occasions
where a convention was misread and the fixture checking it was misread the same
way — the Solov'ev coefficients, `τ` in eq (8e), `DarcyForm`'s `−q`. A
self-consistent tree cannot find that class of error at all.

**Seven configurations, each seeded from its own reference:**

→ **[M-60](MEASUREMENTS.md#m-60)** — case · machine · rel `L2` · rel `L∞`

over ~11,000 nodes per case, aspect ratios 1.4–3.4, elongations 1.3–2.1, one
case with `gg' < 0` throughout.

**THE COMPARISON IS ONE-DIRECTIONAL AND HAS TO BE.** MEQ solves fixed boundary
and freegs4e solves free, so MEQ is handed freegs4e's own **interior flux
surface** as a boundary and its own **profiles** as a source, and has to
reproduce the field inside.

**ROOT SELECTION IS THE FINDING, NOT THE AGREEMENT.** The fixed-boundary problem
with `p'` and `gg'` given as functions of `ψ` has **more than one solution**. On
case A, three:

→ **[M-61](MEASUREMENTS.md#m-61)** — start · `max ψ` · vs reference

The physical one lies **between** the two a ramp reaches, and no amplitude finds
it — the sweep steps from the lower root to the upper between 0.1 and 0.2. That
is the unstable middle branch of an S-curve, which Newton slides off from either
side. Seeded, MEQ takes **two Newton iterations** — 4.5e-03 → 2.7e-09 →
7.1e-15 — and then agrees best of all seven.

**AND WHAT DISTINGUISHED A SECOND SOLUTION FROM AN ERROR WAS ONE SWEEP:
REFINE.** The upper root read **+15.39% → +15.41% across refinement levels 2→4
and degrees 2→3**, converged to six digits. A discretisation error falls; a
different solution does not. Before that measurement this looked like a 15%
disagreement about the equation. It is the mirror of the cone-tiling mistake
recorded above — there a *converging* column was wrongly read as geometry; here
a *flat* one correctly identified a root.

This is the same multiplicity recorded under `CLAUDE_HDGGS.md`'s
*Should `PicardThenNewton` simply be the default?*, where three solve routes reach discrete solutions 9.4% apart,
met from outside the codebase.

**AND IT IS NOW CURED RATHER THAN NAVIGATED, WHICH IS WHAT LET THE REHEARSAL
BECOME SIX SHIPPED EXAMPLES.** → **[M-147](MEASUREMENTS.md#m-147)**. The upper
root is an artefact of the COORDINATE: against `ψ` in Wb/rad the profile table
necessarily stops at the axis, `freegs4e`'s profiles living on `ψ_n ∈ [0, 1]`,
and `meq::SplineProfile` extends a table by a **constant** past its end knots —
so above the axis flux the source is a plateau, and the plateau carries a
solution. `[source] Normalised = true` confines `Ψ` to `[0, 1]` by the
constraint `ψ_ax = max ψ`, the profile is never evaluated off its own table, and
the region that root lived in does not exist. `examples/fixed-*.toml` are six
machines on exactly this footing — `ψ_n = 0.95`, a `257²` reference, twenty
harmonics — agreeing with `freegs4e` to **1.1e-04 worst case** against M-60's
6.8e-03, and each converging from its own ramp with no seed at all.

**ONE SENTENCE ABOVE IS QUALIFIED BY THAT SWEEP AND THE QUALIFICATION MATTERS
MORE THAN THE CURE.** *"no amplitude finds it"* is true of the amplitudes
swept; at `ψ_n = 0.95` a ramp at **12×** the reference's axis height reaches the
physical root, having passed through the upper one at 4× to 8×. So the middle
branch is reachable and the selection is **not monotone in the amplitude** —
which is worse than unreachable for anything shipped, because it looks like a
tunable parameter and is not one.

**FOUR CONVERSIONS, EACH OF WHICH CONVERGES TO A WRONG ANSWER.**

* **freegs4e's `pprime` is `dp/dψ`, and its own docstring says `dp/dψ_n`.** The
  docstring is **wrong** — established four ways, decisively by rebuilding
  `Jtor` from the saved arrays and matching the solver's own to **4e-16**.
  Dividing by `ψ_bnd − ψ_ax` on its word costs a factor that for a tokamak is
  `O(1)` and does not look wrong. This tree records the identical trap from the
  other side in `examples/rotating-density.dat`.
* **Not the separatrix.** Every case is diverted, so the LCFS has an X-point
  **corner** and MXH is a truncated Fourier series that cannot turn one: 2.6e-03
  to 1.3e-02 m where the fitter does **1.9e-05** on a smooth shape.
  `ExtensionConvergence` avoids the separatrix for exactly this reason. An
  interior surface fits to 2–4e-04 m and **stops improving past 10 harmonics**,
  which is what identifies that floor as the extracted contour's rather than the
  fitter's.
* **`ψ = 0` must land on the surface actually handed over.** Mapping with
  `ψ = 0` at `ψ_n = 1` while giving MEQ the `ψ_n = 0.9` surface puts the
  boundary where the source has died — `p'` of 2.7 instead of 3.1e5 — and the
  solve converges in three Newton steps to a field **23× too small**.
* **The guess picks the branch**, above.

**THE SOURCE CONVERSION IS RULED OUT DIRECTLY RATHER THAN ARGUED**: MEQ's tables
evaluated on the reference's own `ψ` reproduce freegs4e's `μ₀ R J_φ` to
**2.3e-05**. Both codes solve the same equation.

**WHAT LIMITS THE REMAINING 0.06–0.7%** is most likely geometric: the MXH fit is
2–4e-04 m against minor radii of 0.24–0.61 m, and the contour it fits comes from
freegs4e's own 129² grid. Refining that grid is the test and is **not done**.

**TWO SELF-TESTS IN THE HARNESS, AND BOTH FOUND REAL BUGS IN THEMSELVES.**
`mxh.py` recovers a known shape to 1.9e-05 m; writing it found that the `tR`
branch switches at the **R extrema** rather than at the midplane — true only for
an up-down symmetric surface — worth 25×, and that its distance measure was
reading its own 2048-point sampling, **pinned at 1.6e-03 while the fit improved
25-fold**. `convert.py` checks its derivative column against a closed form at
1.7e-07 and carries a control that drops the chain rule and reads 0.52.

**THAT MAKES FOUR TIMES IN ONE SESSION AN INSTRUMENT WAS MISTAKEN FOR A
RESULT** — the cone tiling, the MXH distance, the transmission row's shared
quadrature rule, and the FB-1b single-mesh check. The tell is always the same: a
number that does not move when something that should move it changes.

**AND IT IS `ψ*` THAT IS COMPARED, NOT `ψ_h`.** `apps/meq.cpp` samples
`postProcessedPotential()` into the `.nc` — its own comment says *"THE POTENTIAL
SAMPLED HERE IS psi*, NOT psi_h"* — and `compare.py` reads the `.nc`, so every
benchmark number is the post-processed field. That is the right thing to compare,
since the `.nc` is what a consumer receives, but `ψ*` converges at **k+2**, so
the accuracy-per-dof result below must not be read as a statement about `k+1`.
`FreeBoundaryCoupling`'s rates use `solver.potential()` and ARE `ψ_h`'s; the two
are different measurements on purpose. Comparing `ψ_h` would need the `.gf`
route, which is a second reason to want pyMFEM.

**AND THE FLOORS CAME OFF ONE AT A TIME, ENDING AT A FACTOR OF 17.** Carrying
the refinement to 513² and then chasing what was left is a sequence in which
**each step's diagnosis was the previous step's mistake**:

→ **[M-62](MEASUREMENTS.md#m-62)** — changed · MXH fit · MEQ rel `L2` · limited by

The fit stopping at 513² identified **the fitter** rather than the contour — the
opposite of what this file and `tools/README.md` said, and they were right *at
129²*, where ten harmonics happen to reach the contour's own limit. Raising the
harmonics then bought a 4.7× better boundary and **`L2` did not move at all**,
which is what identified **MEQ's own discretisation** as the binding constraint
— the first time in this benchmark's history that has been true. The last two
rows are controls: `k=3 r=3` moves it 0.7%, and a sixteenfold denser comparison
grid moves it 0.03%.

**`L∞` behaves differently and confirms it**: it DID follow the boundary fit,
2.900e-04 → 1.347e-04 for the harmonics, while `L2` sat still. A shape error
shows up worst near the boundary, which a maximum norm sees and an `L2` over the
core does not.

**SO THE ANSWER TO "DO THEY AGREE TO ROUND-OFF" IS NO, AND CANNOT BE YES HERE.**
They agree to **8.8e-06**, and what is left is the boundary fit at 2.2e-05 m
against a minor radius of 0.42 m. Round-off is unreachable *through this
comparison* by construction, because the comparison contains a fit. **FB-6
removes it.**

**THE SATURATION WAS THE REFERENCE'S AND IT IS THE BOUNDARY FIT, MEASURED
 This file has said since it was written that MEQ saturates the
benchmark and that the next work is to refine the reference. Now measured, on
case A, by refining the reference grid and rerunning MEQ unchanged:

→ **[M-63](MEASUREMENTS.md#m-63)** — reference · MXH shape fit · MEQ rel `L2` · rel `L∞`

**MEQ's error tracks the shape fit.** The LCFS handed to MEQ is a fit to a
contour extracted from the reference's grid; the contour improves like `h`, and
the answer improves with it. Extrapolated, 2049² would reach about **8e-06** —
an order of magnitude, **not round-off**. Grid refinement alone cannot get
there, because the fit is linear in `h` and is in the budget at all *only
because this is the fixed-boundary rehearsal*. **FB-6 removes it entirely**: a
free-boundary MEQ takes the same coils and profiles and never sees an LCFS.

**THE CEILING IS 1025², AND MEMORY IS WHAT SETS IT.** The reference at 1025²
needs about **6.0 GB**, which this machine has; the cached boundary matrix there
would need **35 GB**, which it has not, so that rung falls back to the shipped
loop and is slow rather than impossible. 2049² is out on both counts. **513² is
the rung that is both exact and fast** — 72.75 s with the matrix cached — and is
the one stored.

**FOUR THINGS ABOUT `freegs4e` CAME OUT OF THE ATTEMPT.**

* **IT DOES NOT PARALLELISE, AND THAT IS NOT WHY IT IS SLOW.** Checked
  2026-09-06: three files carry `numba.njit`, **none** with `prange` or
  `parallel=True`, and there is **no `multiprocessing`, `joblib`, `concurrent`
  or `threading` anywhere in the package**. It is serial by construction apart
  from whatever BLAS numpy pulls in — and the hot loop calls no BLAS.
* **THE REAL COST IS THAT IT RECOMPUTES A GEOMETRY-ONLY MATRIX EVERY PICARD
  STEP, AND CACHING IT IS WORTH 163x.** `boundary.freeBoundary` loops over the
  `4n` boundary points and calls `Greens( R, Z, R[x,y], Z[x,y] )` over the whole
  grid for each — elliptic integrals, `O(n²)` of them, `4n` times. **`Greens`
  depends only on the geometry**, and `romb` is linear, so the whole boundary
  condition is one fixed matrix `M` applied to `Jtor`: `psi_bndry = M Jtor`.
  `M` is the same at every iteration and is rebuilt at every iteration.

  Measured at `n = 129`, against freegs4e's own loop:

  | | |
  |---|---|
  | the loop, per Picard step | **518.93 ms** |
  | building `M` once | 496.82 ms, **66 MiB** |
  | `M @ Jtor`, per Picard step | **3.19 ms** |
  | agreement | **4.6e-16** relative — the same arithmetic, reassociated |

  **163x per step**, and over a whole run — 23 to 103 Picard steps are what this
  benchmark sees — **20.9x to 64.8x**. The matrix is 0.51 GB at `n = 257` and
  **4.02 GB at `n = 513`, measured**, so caching is comfortable at both; 1025²
  would need 35 GB and is out, where `CACHE_BOUNDARY_MAX_GB` falls back to the
  loop. **A refinement not
  measured**: `Jtor` is zero outside the plasma, so only its support's columns
  are needed, which should cut the matrix by about an order of magnitude — at
  the cost of tracking a support that moves as Picard runs.

  **`../freegs4e` IS NOT EDITED — IT IS CACHED FROM `fgsref.py`, AS
  WORKAROUND 6.** `M` is built once and kept on the Equilibrium; the Romberg
  weights come from `romb( I )` rather than a reimplementation, so the rule
  cannot drift from the shipped loop's. Measured: the boundary values agree to
  **5.5e-16**, the per-step cost goes 651.5 ms → **8.07 ms** at 129² and
  4271 ms → **24.0 ms** at 257² (81× and 178×), and the whole run goes
  **29.2 s → 4.9 s** at 129² and **175.1 s → 13.3 s** at 257². Memory is the
  constraint — 0.54 GB at 257², **4.3 GB at 513²**, 35 GB at 1025² — so
  `CACHE_BOUNDARY_MAX_GB` falls back to the loop rather than risking an OOM.
  **This is what makes an exact 513² reference affordable**, which is what von
  Hagenow was wanted for and is better than it: von Hagenow is a different
  method converging at `O(h)`, this is the same arithmetic reassociated.

  **AND THE CONVERGED ANSWERS DIFFER BY 4.5e-08, WHICH IS `freegs4e`'s OWN
  FLOOR AND NOT THE REASSOCIATION.** Two hypotheses were wrong before the
  control settled it. The boundary condition on the real converged `Jtor` agrees
  to **6.8e-16**, so the arithmetic is exact. It is **not the stopping rule**:
  tightening the Picard tolerance 1e-9 → 1e-11 → 1e-13, which takes 41 → 60 →
  401 iterations, leaves the gap at 4.511e-08 **to four digits**. It is **not
  amplification**: perturbing `M` by 1e-12, 1e-10 and 1e-8 moves the field by
  3.39e-08, 6.50e-08 and 7.71e-08 — **four orders in, a factor of 2.3 out**, a
  floor rather than a gain. And each variant is **bit-reproducible against
  itself**, `0.000e+00`, which is the control that says the difference is real.
  So this reference is determinate to about **5e-08 whatever the arithmetic
  does** — five orders below its own grid error, so irrelevant to use and
  material to what may be claimed.

  **WE USED TO DRIVE IT COLD, AND `--seed-from` IS THE FIX — WORTH 1.21x, NOT
  THE ORDER OF MAGNITUDE THIS ENTRY PREDICTED.** `fgsref.py` keeps freegs4e's
  `2^n + 1` grid convention **specifically** so that 129, 257, 513 and 1025 nest
  point for point, and `--seed-from=auto` now lifts the converged coarse
  `plasma_psi` onto the fine grid with a cubic `RectBivariateSpline` and starts
  Picard there. Measured on `H_limited_circular`, 129² → 257²: **39 Picard steps
  and 152.7 s cold against 31 and 126.1 s seeded**.

  **THE PREDICTION WAS WRONG AND THE MECHANISM IS WHY.** This entry said a cold
  fine run is "hours where a seeded one should be a handful of steps". That
  assumes the cost of a Picard run is set by where it starts. It is set by how
  far it has to go and by a contraction the seed does not change: here about
  **0.60 a step**, and a *converged* coarse answer is still **1.6% wrong** on the
  fine grid — that being the reference's own discretisation error, the quantity
  the refinement exists to measure. So the seed starts a factor of ~62 closer and
  buys `log(62)/log(1/0.60) ≈ 8` steps **at the top of the run and nothing
  after**. A fixed number of steps, so proportionally less the longer the run.

  **The answers are NOT bit-identical — 3.1e-08 relative in `ψ_ax`** — which is
  four orders below the 1.6% being measured, so it does not move the benchmark.
  Picard converges to the fine grid's own solution whatever it starts from; what
  differs at the eighth digit is where `rtol = 1e-9` happened to bite. The
  nesting is **asserted**: a coarse grid that does not divide, a seed finer than
  the run, and a coarse file with different extents are all refused, while a
  *missing* seed warns and runs cold — seeding is an optimisation, so its absence
  costs time and never correctness.

  **AND ONE OF THAT FILE'S OWN COMMENTS WAS WRONG ABOUT WHY IT IS SLOW.** It
  attributed the 85 s → 656 s to *"the `n³` of a 2D sparse LU"*. It is not:
  `multigrid.MGDirect.__init__` calls `factorized( A )` **once** at construction
  and `__call__` only backsolves, so the LU is not a per-step cost at all — and
  the `HAGENOW` note a few lines below said so correctly the whole time. The two
  comments contradicted each other and the measurement settles it in favour of
  the second. **The consequence is that `VCYCLE_LEVELS` is a lever on the wrong
  term**: turning the multigrid on buys almost nothing, so leaving it off is
  right for a reason the file did not give. 
* **Its cost is the boundary condition, not the solve.**
  `Equilibrium.__init__` takes `boundary=freeBoundary`, whose own docstring
  calls it *"an integral over the area of the domain for each point"* — it loops
  the `4n` boundary points and evaluates `Greens` over the whole `n²` grid for
  each, so it is **`O(n³)` per Picard step**. Measured 85 s at 129² against
  **656 s** at 257², which is `2³` and not `2²`.
  `boundary.freeBoundaryHagenow` is in the same file, is `O(n²)`, runs 29 / 71 /
  230 s at 129² / 257² / 513², and is not the default.
* **The multigrid is built and switched off.** `createVcycle( …, nlevels=1, … )`
  — at one level a V-cycle is a direct sparse solve on the full grid.
* **`setSolverVcycle()` hard-codes the SECOND-order generator** and ignores
  `Equilibrium.order`, so calling it on a 4th-order equilibrium silently solves a
  different problem. A V-cycle built correctly on the 4th-order generator does
  not converge at all — `ValueError: No opoints found!`.
* **ITS TWO BOUNDARY CONDITIONS DISAGREE BY 3.2e-03 WITH A FLAT GAP, AND THE
  FAST ONE IS INCONSISTENT — ONE CONSTANT.**
  `boundary.freeBoundaryHagenow` displaces its observation point off the
  boundary by a hard-coded **`eps = 1e-2` METRES**, *"to avoid the singularity in
  `G(R,R')` when `R'=R`"* — a fixed **physical** distance, so the
  `O( eps·∂ψ/∂n )` error does not shrink with the grid and the method converges
  to a different answer. **A gap that does not shrink between two
  discretisations of one problem means one of them is not a discretisation of
  it**, and reading it as an open question about which to believe was the
  mistake.

  Measured on one fixed `Jtor` with no solve in the way — the direct Green's
  integral is ground truth, so what is printed is von Hagenow's own error —
  the observed order in `h` at `n = 65, 129, 257` is **0.63, 0.21, 0.07** at the
  shipped `eps` and **0.90, 1.00, 0.99** at `eps = 0.2 h`. On the real solve,
  `H_limited_circular` against the direct integral, the relative `L2` gap goes
  1.008e-02 → 1.032e-02 shipped (**rate −0.04**, slightly worse) and
  2.139e-03 → **1.093e-03** scaled (**rate 0.97**).

  **AND THE OPTIMUM IS NOT "AS SMALL AS POSSIBLE"**: swept at a fixed 129² grid
  the gap reads 7.6e-03, **3.3e-03**, 8.9e-03, 1.6e-02, 2.3e-02 at `eps` = 1e-2,
  3e-3, 1e-3, 3e-4, 1e-4 — below about `0.2 h` the log spike is narrower than the
  cell and the Romberg rule misses it.

  **FIRST ORDER IS THE CEILING AND THE RESIDUAL IS THE QUADRATURE, NOT THE
  DISPLACEMENT.** `∮G σ dl` is a single-layer potential, continuous across the
  boundary with its normal derivative jumping, so averaging `+eps` and `−eps`
  ought to cancel the `O(eps)` term. Measured, it does not — same rate, slightly
  worse — which is what locates the error in the near-singular quadrature.
  Getting past first order needs the log subtracted analytically, i.e. a
  boundary-element method, and is deliberately not attempted.

  **WHAT IT BUYS**: extrapolated, the boundary error at 513² is about 2.7e-04 in
  `ψ_ax` against that grid's own 7.9e-04 on this case, so it does not dominate —
  and 513² costs **230 s** instead of hours. That is what makes the reference
  refinable, which is the binding constraint on this benchmark.
  **`freegs4e` is not edited**: one constant, tied to the cell, reinstalled from
  `fgsref.py` as WORKAROUND 5. `--hagenow` is the usable version and
  `--hagenow-eps=shipped` reproduces `freegs4e` exactly.

**COST, AND THE ONLY COLUMN THAT MEANS ANYTHING IS ACCURACY PER UNKNOWN.** A
wall-clock ratio is not a statement about either code: freegs4e converges a
FREE-boundary equilibrium — coil Green's functions, X-point finding, a control
system, 23–103 Picard steps — while MEQ solves the FIXED-boundary problem inside
a surface it is handed, in C++ against Python. Case A, `MKL=1`, `OMP=1`:

→ **[M-64](MEASUREMENTS.md#m-64)** — `k` · refine · elements · dofs · dofs/ref pt · wall · rel L2

against freegs4e's own 84.7 s at 129² = 16,641 grid points.

**High order is worth about 5x in unknowns**: `k = 3` on 397 elements reaches
1.7e-04 with 11,910 dofs in 3.0 s, where `k = 1` needs 7,185 elements, 64,665
dofs and 13.6 s for the same thing.

**AND MEQ SATURATES THE BENCHMARK, WHICH IS THE MOST USEFUL LINE IN THE TABLE.**
Errors stop at 1.44e-04 on case A, and at 5e-03 / 2.4e-03 on MAST and DIII-D,
because that is the REFERENCE's accuracy — the MXH fit and the contour off its
129² grid. Refining MEQ past the second row buys nothing. **Future work on this
benchmark should refine the reference, not MEQ.**

Timings are the whole driver; the `.nc` sampling is ≤0.5 s of it, measured by
rerunning at a 17² output grid. Read them as an order of magnitude.

**`k = 2, refine = 1` DOES NOT CONVERGE on MAST or DIII-D while `k = 3` does on
the same mesh** — Newton fails, `PicardThenNewton` fails, and the driver's advice
to raise the degree is what works. That is `p`-refinement reaching a case neither
`h` nor a globalisation does, which is the pedestal finding from a new direction.

**Deliberately not established**: nothing about free boundary, since MEQ does
not solve it; nothing about the flux-surface machinery, since only `ψ` on a grid
is differenced and freegs4e computes no averages; and nothing about a real
reconstruction, since both codes are given analytic profile shapes.

**AND THERE IS A LIMITED CASE NOW, WHICH IS THE ONLY ONE MEQ COULD EVER COMPARE
AGAINST.** `fgsref.py`'s `H_limited_circular`, : vertical-field
coils only, so no X-point exists in range and the boundary is the flux surface
through the limiter. Every other case in that table is **diverted**, and MEQ's
moving support is a pointwise test on `Ψ` — see `FREE-BOUNDARY-PLAN.md` §10.3.
Its coil currents were **solved for** by freegs4e's own control system rather
than guessed, which is what makes them consistent input for both codes.

**IT IS NOT CONVERGED AT 129², AND THE DIVERTED CASES ARE.**

→ **[M-65](MEASUREMENTS.md#m-65)** — grid · `ψ_ax` · `ψ_bnd`

Successive differences fall by **4.98** each time — order ≈ 2.3 — and Richardson
gives `ψ_ax ≈ 9.3014e-02`, `ψ_bnd ≈ 2.6158e-02`, against which **129² is 2.0%
and 6.3% out**. A **limited** boundary is a maximum over the limiter ring, a
pointwise operation on a discrete set, where a diverted one is a saddle located
by interpolation — so it converges more slowly in the grid. Use the 513² run.

**AND THAT 513² RUN NOW EXISTS AS A STORED REFERENCE**, produced
with the exact boundary condition cached — `tools/freegs4e-benchmark/ref-n513/`,
with `baseline.json` the machine-readable record and `baseline.py` the
regenerator. The whole nested ladder, each rung seeded from the one below:

→ **[M-66](MEASUREMENTS.md#m-66)** — `n` · points · `ψ_ax` · Picard · wall · peak RSS

**72.75 s, of which 37.4 s is building the 4.02 GB matrix once** — against
*hours* for the shipped loop and 230 s for von Hagenow, so the cache is **3×
faster than the approximate method and exact**. That is why FB-6 does not need
von Hagenow after all.

**THE NUMBER FB-6 HAS TO BEAT IS 7.92e-04**, which is how far 513² sits from its
own Richardson limit. A MEQ result closer to the 513² values than that is
measuring the REFERENCE's grid error rather than MEQ. Beneath it sits freegs4e's
~5e-08 indeterminacy, four orders lower and so never binding.

**THIS MACHINE HAS 23 GB, 21 of them available**, and the WSL2 allocation can be
raised, so a memory ceiling recorded here is a claim about the machine on the
day. **Re-check `free` rather than trusting a number in this file**: the ceiling
is a property of the day, and it has already moved once.

**THE RICHARDSON SIGN WAS GOT WRONG ON THE FIRST PASS AND THIS FILE CAUGHT IT.**
The sequence decreases toward its limit, so the extrapolate is BELOW the finest
rung; adding the correction rather than subtracting it gave 0.093161 against the
9.3014e-02 recorded above. **An independently recorded number refuted a fresh
calculation**, which is the argument for writing one down the first time it is
measured.

**AND IT IS BEATEN. MEQ REPRODUCES THE CONVERGED REFERENCE TO
1.5e-04**, and what was in the way was the **limiter border**.
`FREE-BOUNDARY-PLAN.md` §7.20 is the record:

→ **[M-67](MEASUREMENTS.md#m-67)** — freegs4e, 513² · MEQ, `k = 3`, 26375 el · apart

both **inside** the 7.92e-04 above, so what is being measured is the reference's
grid error rather than MEQ; MEQ's own remaining discretisation error is
**3.2e-05** from its last two rungs. **Two inputs had to be taken from the same
run being compared against**, and neither is obvious: freegs4e's coil currents
are an OUTPUT of its control system and move **14%** on P2 between 129² and 513²,
and its limiter contact is the maximum over a ring of grid CELLS, which at 129²
sits on the outboard midplane and at 513² on the **inboard shoulder, 0.6 m
away** — the true circle's maximum is on the inboard shoulder at both, so the
129² contact is on the wrong side of the machine.

**pyMFEM would improve this and is not used.** `mkguess.py` hand-writes MFEM's
ASCII mesh and GridFunction format to seed the restart, and `compare.py` reads
the **lossy** `.nc` rather than `_psi.gf`. Reading the exact `P_k` coefficients
would take the grid sampling, the band mask and the interpolation out of the
error budget — which matters now that the residual is at 1e-04.

## MEQ solves MAST-U from a configuration file alone

**THE GUESS IS MEQ'S WORK AND NOT THE USER'S**, and until recently every machine
example in this tree said the opposite — each was handed a `.gf` that
`mkexactguess.py` had reconstructed by summing Green's functions over the
source, which is the solver asking the user to do its convergence work.
`[initialguess] Type = "conductors"` is the replacement and it reads only what
is already in the file.

**IT IS IN TWO PARTS AND THE SECOND IS WHAT SELECTS THE BRANCH.**

* **The conductors' own vacuum field**, summed over `[[coils]]` by
  `meq::coilPsi`. It is Δ*-harmonic off the conductors at measured rate 2.00,
  so it is an exact solution of the vacuum problem rather than an approximation
  of one, and it carries the machine's scale and topology.
* **`I_p` over an ELLIPTICAL COLUMN** about a guessed magnetic axis, by
  `meq::ellipsePsi`. `[source] PlasmaCurrent` is already in the file, so this
  costs the user ONE number — where the plasma is.

**THE VACUUM FIELD ALONE CONVERGES TO A BRANCH THAT IS NOT ONE**, which is the
measurement that makes the second part necessary rather than an improvement:
MAST-U converges in 9 iterations to `psi_ax` 3.5× the reference, with its axis
at ( 2.65, −0.91 ) — outside the machine.

**AND THE COLUMN IS AN ELLIPSE BECAUSE OF WHAT THE OTHER TWO SHAPES DO.** A
FILAMENT is singular at exactly the point `locateAxisPoint()` has to find —
measured, 6.667e-01 at 3 mm from the guessed axis against MAST-U's reference
`psi_axis` of 9.187e-02, and still climbing — and the solve fails. A RECTANGLE
is bounded, `meq::Coil` already is one, and a uniform current density over it
carries a logarithm in its second derivatives at each of four corners; those are
artefacts of the shape rather than anything the equilibrium puts there.

→ **[M-124](MEASUREMENTS.md#m-124)** — the three guesses against the reference ·
the approach to the axis · `ellipsePsi`'s two quadrature rules and what one rule
alone does

**AND IT IS THE MOST EXPENSIVE THING IN THE RUN IF ITS QUADRATURE IS LEFT AT THE
DEFAULT, WHICH IS A REFERENCE ORDER.** The guess evaluates every conductor at
every nodal point of the potential space AND the trace space; at
`defaultCoilQuadratureOrder = 32` that is 1392 us a point for MAST-U's 23
conductors, and it read **93.6% of that machine's whole run** — against a solve
it was helping by 26 s. `meq::guessCoilQuadratureOrder = 6` is the measured
replacement, wrong by a thousandth of the `psi_ax` it is guessing, and
`prepare()` now caches the seed so the driver's prepare and `solve()`'s do the
projection once between them. **275 s to 21 s, with every printed digit of the
answer unchanged.** The guess is not the problem statement, which is what makes
its accuracy a free choice.

→ **[M-128](MEASUREMENTS.md#m-128)** — the profile that found it · the order
against accuracy and cost · what did not move

**AND IT TAKES `BorderRegularisation` AS WELL.** M-123's repair and this one are
necessary together and neither is sufficient: without the damping MAST-U fails
whatever the guess, and without the guess it converges to the wrong branch. The
committed example was measured failing with only one of them in it.

**THE SHAPE OF THE COLUMN DOES NOT MATTER AND THAT IS WHAT MAKES IT USABLE.** A
column tuned to MAST-U's own plasma and a round one at the default `0.5*CentreR`
reach the same `psi_ax`, `psi_bnd` and X-point to every printed digit in the
same 12 Newton iterations.

### `UpDownSymmetry` is sound, and it needs a mesh made for it

**A DOUBLE NULL SHOULD NOT NEED TO BE TOLD WHICH SADDLE TO FOLLOW**, and
`[boundary.xpoint]` makes it: XP-3's border follows ONE saddle, chosen by its
seed. `[solver] UpDownSymmetry` is the alternative — project the iterate onto
the subspace of fields even in `z`, and the two saddles are exactly degenerate
when the search looks for them.

**IT IS A DOF-FOR-DOF AVERAGE, SO IT NEEDS A MIRROR-SYMMETRIC MESH, AND AN
ORDINARY ONE IS NOT.** `setUpDownSymmetry()` REFUSES such a mesh by name rather
than projecting onto something that is not a reflection. On MAST-U's committed
mesh 3624 of 4735 vertices have no mirror partner — and **the machine is not the
problem**: its 23 conductors are mirror-paired to the last digit, one pair
excepted whose currents are equal and opposite at 8e-11 of the total. **gmsh's
triangulation of a symmetric geometry is not symmetric**; it picks a diagonal
and picks freely.

**SO THE MESHER MAKES THE MESH INSTEAD OF GMSH CHOOSING IT.**
`[mesh.generate] Symmetric = true` reaches `halfdisc.py --symmetric`, which
meshes `z >= 0` and reflects it — 0 unpaired nodes and 0 involution failures on
MAST-U's own geometry, at a worst pairing discrepancy of **0.000e+00 m** rather
than a small one, because the image coordinate is the written decimal with its
sign flipped. It REFUSES a geometry that is not itself mirror-symmetric rather
than reflecting a machine into a different machine.

**AND `[solver] UpDownSymmetry` ON A GENERATED MESH WITHOUT IT IS A PARSE
ERROR**, because `buildMirrorMaps()` runs inside the bordered Newton driver and
not in `prepare()`: the refusal otherwise arrives after the mesh, the spaces,
the assembly and the guess — measured, past 150 s on MAST-U at `k = 2`. Where
the file also says how the mesh is made, MEQ knows the answer at parse time.

→ **[M-125](MEASUREMENTS.md#m-125)** — the projection as the identity, in every
block · the mesh · the two defects the acceptance case found —
and **[M-127](MEASUREMENTS.md#m-127)** — the symmetric mesher, its exactness
and its refusals

**AND `MakeCartesian2D`'s TRIANGLES ARE NOT A SYMMETRIC MESH EITHER**: it splits
every cell along one diagonal, so a box symmetric in `z` has a triangulation
that is not. `LimiterCurve.cpp` records the same fact from the other side. The
QUADRILATERAL variant has no diagonal to choose and is what
`tests/convergence/UpDownSymmetry.cpp` solves on.

## The limiter border read the remainder, and the case that found it is the free-boundary DESC race

**THE FOURTH CONSUMER OF `psi = psi_p + psi_c` TO READ THE REMAINDER, AND THE
FIRST THAT COULD ONLY BE REACHED BY THREE FEATURES AT ONCE.**
`GradShafranovSolver::solve()` splits the limiter contact into `limiterValue()`,
the border row as a linear functional, and `limiterTotal()`, the physical flux —
M-148's lesson, in the code, correctly. `limiterTotal()` then enumerated **three**
routes to the contact and there are **four**, and the one left out is the
default: a PRESCRIBED point, `[boundary.limiter] R` and `Z`, under
`LimiterConstraint::ExactPoint` with no X-point unknown. Nothing sets
`limiterContactLocatedValue` there and `xPointIsUnknown` is false, so every
branch falls through.

→ **[M-161](MEASUREMENTS.md#m-161)** — the four routes · `psi_c` at the contact
against the span · the defective run's own numbers · the cold start that was a
symptom

**IT DOES NOT PERTURB THE ANSWER, IT CHOOSES A DIFFERENT ONE.** `psi_c` at the
contact is **0.50** of the span on the 513² reference and **1.18** on the 129²
one, so the run converges — 27 Newton steps, every border at machine zero, `I_p`
to seven figures — onto a branch with `psi_bnd > psi_ax`, a negative profile
scale, no O-point and `r = 0` inside the plasma. What caught it is the
source-on-the-axis refusal, firing on the consequence.

**THE `profile scale` COLUMN IS THE DIAGNOSTIC AND IT IS FREE.** It is the
unknown `[source] PlasmaCurrent` closes, so it lands on **1** exactly when the
profile tables carry the amplitudes the reference converged to. It reads
0.9999423 on the repaired run and −2.1e-04 on the defective one. It is the one
printed number that separates *converged* from *converged to the right thing*,
and it is worth reading on every coil-subtracted run.

**AND IT ALSO SAYS THE COLD START WAS NEVER THIS CASE'S PROPERTY.**
`examples/limited-tokamak.toml`'s header records that a cold start does not
converge here, and the cold DESIGN guess behaved exactly as predicted while the
border was wrong. With it fixed the same guess converges in **18** Newton
iterations over 4 sweeps, against **127** for the exact guess built from the
reference's own `Jtor` — the same answer to every digit, and the warm start is
seven times the work.

## DESC on a free boundary: only a limited case reaches it, and the solve then does not converge

**DESC's `BoundaryError` CONSTRAINS THE LCFS, AND THE LCFS IS A
`FourierRZToroidalSurface` — A TRUNCATED FOURIER SERIES.** Every free-boundary
machine in `examples/` but one is DIVERTED, so its LCFS is a separatrix through
an X-point and has a **corner**, which no truncation turns. That is the same
fact `CLAUDE.md` records from the other side as the reason the
`examples/fixed-*.toml` ladder is posed at `psi_n = 0.95`. A LIMITED plasma's
edge is the smooth flux surface tangent to the limiter and is representable, so
`examples/limited-tokamak-filament.toml` is the case the free-boundary
comparison exists on and `tools/desc-benchmark/descfreeb.py` is its runner.

**AND THE DESC ARM OF IT DOES NOT CONVERGE.** `descfreeb.py` never read
`result["success"]`, so a run that stopped at its starting point was reported as
an answer — and a boundary seeded from the reference's own LCFS starts at a
stationary point, so it reports the reference back with a tiny error **because
it did not move**. Re-taken with the check, **3 of 13 runs converge** and none
of the three is within 1e-02. The three-pairing agreement at 1.3e-04 this
section used to carry is withdrawn. What stands is the conversion's own audits,
which are checked against the reference by routes that do not use the optimiser.
**The DESC comparison that is built on is the FIXED ladder.**

→ **[M-162](MEASUREMENTS.md#m-162)** — the retraction · the verified tally ·
the four harness defects · the conversion's self-audits · the toroidal field

**AND THE COMPARISON THAT TABLE SITS IN IS NOT ONE PROBLEM, WHICH ONLY TWO
FURTHER CODES COULD SHOW.** `freegs4e` takes `psi_bndry` to be the maximum of
`psi` over the innermost ring of grid NODES inside its wall, so the reference's
boundary flux is the value at one node `0.86 h` inside the limiter and its
plasma never touches its own limiter anywhere.
`examples/limited-tokamak-filament.toml` prescribes that node, which is what
makes MEQ's agreement on this case well posed and also means **MEQ is handed
the number it is then scored on**. NICE finds its contact on the true circle
and TSC is given a different point, so the accuracy column compares solutions
of at least two problems and is reported per posture rather than as one column.
**MEQ TAKES NICE's POSTURE TOO, AND `examples/limited-tokamak-filament-curve
.toml` IS IT** — the same machine with the contact FOUND on a meshed limiter
circle, 3 support sweeps and 3.817 s to a `psi_ax` of 9.273894e-02 with the
contact at `( 0.81699, 0.297389 )`. That is **0.38% from the point posture's
`psi_ax`**, which is the 1.31% in `psi_bnd` above propagating, so the two files
are one machine under two boundary conditions rather than one problem solved
twice. **The column stays two columns anyway**: three codes, three postures,
and each is scoreable only against the reference its own posture was taken
from.

→ **[M-164](MEASUREMENTS.md#m-164)** — NICE and TSC timed · which codes can
take eight threads and what it buys them · the staircase ring reproduced to
every digit · §4 REPAIRED, and its retraction points at M-165

## The driver's first support freeze read a different functional under the same words

**AND IT IS WHY THE CURVE POSTURE LOOKED IMPOSSIBLE.** `apps/meq.cpp`'s
`edgeFluxOf()` supplies `psi_bnd` of the iterate the first sweep starts from —
the one number the freeze needs and the one `psiBoundary()` cannot give before
a solve. Under `[boundary.limiter] SurfaceAttribute` it read *"the maximum over
the meshed limiter surface"* as the maximum over the region the limiter
**encloses**, where `setLimiterSurface()` constrains `psi_bnd` to the maximum
over the faces **bounding** it. **The enclosed region contains the magnetic
axis**, so the estimate was `psi_ax`: 7.545003e-02 against a `psiAxisGuess` of
7.543374e-02, a span of **−1.63e-05 where the true one is +4.76e-02** — wrong
by 2.9e+03 and of the wrong sign, so the support inverted and the plasma became
everything the plasma is not.

**THE FIX IS ONE IMPLEMENTATION, NOT A CORRECTED SECOND ONE.**
`meq::collectLimiterPolygon()` and `meq::limiterPolygonMaximum()` are free
functions holding what `locateLimiterContact()` had inside it; the border calls
them and so does the driver, on its own field at `offset = 0`. A free function
because the driver needs the constraint **before any solver exists** — the
Picard pre-stage's home solver is built unbordered and never hears of the
limiter at all — so no method on `GradShafranovSolver` could serve all three
call sites without an ordering rule.

**THE BOUND THAT LET IT THROUGH IS THE TRANSFERABLE PART.** That helper's own
comment says *"a bad estimate costs sweeps rather than correctness"*, and it is
right: the sweep loop is a fixed point over the support, so any estimate in the
neighbourhood is recovered. **An estimate of a DIFFERENT QUANTITY is not in the
neighbourhood**, and the loop converged — from two starting points a whole
equilibrium apart, to the same annulus in every printed digit. *A bound on how
wrong an estimate may be is a bound on the estimate of that quantity.* The same
comment recorded the escape in the next breath — *"THE SurfaceAttribute BRANCH
IS EXERCISED BY NO SHIPPED FIXTURE"* — written as a note on the cost of
building one.

→ **[M-165](MEASUREMENTS.md#m-165)** — both readings of one guess · the
repaired run against the defective one · the regression's own control, where
the region maximum reproduces `psi_ax` on a fixture that knows nothing about
the machine

**THE CONDUCTOR MODEL IS EXACT ON ALL THREE ARMS, WHICH NO OTHER CASE CAN SAY.**
`freegs4e`'s `H_limited_circular` is point filaments, MEQ's `Model = "filament"`
is point filaments, and DESC's `FourierPlanarCoil` with one `r_n` is a circular
loop — checked against the closed form on the loop's own axis to **nine
figures**. `examples/limited-tokamak.toml` carries 0.1 × 0.1 m rectangles
against the same filament reference, a `( w/d )^2` term of about 1.6e-02 on the
near field, so its agreement is reached despite a modelling difference rather
than because the models agree.

**THE EXTERNAL FIELD IS NOT THE `[[coils]]` BLOCKS, AND THIS IS THE ITEM MOST
LIKELY TO BE REDISCOVERED.** `BoundaryError`'s second residual is
`B_out² − B_in² − 2 mu0 p = 0`, and `B_in` carries `B_phi = g/R` while poloidal
field coils produce none. A tokamak's toroidal field comes from a TF coil that
**no Grad–Shafranov input names**, because GS sees `g` only through `g dg/dpsi`
and MEQ's answer does not depend on the constant. With the PF coils alone the
pressure residual starts at **0.90 normalised** — `B_phi²` entire — the
optimiser trades the good normal-field residual against the impossible one, and
the boundary walks **0.81 m on a plasma of minor radius 0.34**. It does not
fail; it converges to a different machine. `g_at_gamma` supplies it, and that is
the same scalar the conversion already needs for the total toroidal flux.

**AND THE TWO CODES ARE NOT GIVEN THE SAME STATEMENT.** MEQ is given the
currents, `p'`, `g g'`, a target `I_p` and **a limiter contact**; DESC is given
the currents, `p( rho )`, `I( rho )` and **the total toroidal flux**, with no
wall and no limiter in the formulation at all. The plasma's size is pinned by
the contact in one and by `Psi` in the other. Both are complete statements of
one equilibrium and a disagreement can live in that difference as well as in
either discretisation.

**DESC'S FREE-BOUNDARY PROBLEM HAS A SECOND BRANCH ITS OWN OBJECTIVE PREFERS.**
From a circle of the machine's design size rather than the reference's LCFS it
converges to a boundary **0.19 m** away with `psi_ax` 15% out — at a residual
**lower** than the reference boundary's. Not under-converged: a better minimum,
at a different equilibrium. M-26 in DESC's coordinates, and the reason the
matched posing hands both codes the branch rather than letting either search.
