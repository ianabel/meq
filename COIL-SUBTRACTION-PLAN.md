# The coils out of the mesh: `psi = psi_coils_analytic + psi_solved`

Ian's idea, 2026-09-18, in his words: *"possibly we could have inside-coil
approximation that uses the fact that the near field of a rectangular coil is
solvable? approximate psi_coil as psi_analytic_coil + polynomial"*, against the
requirement that *"we should be faster at 1e-3 and 1e-4 not just at very high
accuracy, high order should allow us to not use very many elements"*.

`CLAUDE.md` is the index and `CLAUDE_FB.md` owns the free boundary; this defers
to both. It is a DESIGN, not a measurement — nothing here has been built.

---

## 0. What would falsify it, and the cheapest experiment

Put first, because this plan has one load-bearing assumption and it is testable
in an afternoon without writing any of it.

**MEASURED, 2026-09-18, and the answer is split.**
`tools/freegs4e-benchmark/conductor_model.py` builds both fields from the same
18 DIII-D currents on the reference's own grid — freegs4e's point filaments
against MEQ's uniform rectangles, 24² Gauss per rectangle — with no solver and
no mesh anywhere:

| region | model difference | M-111's irreducible error | |
|---|---|---|---|
| everywhere | **6.081e-03** | **5.785e-03** | **1.05×** |
| off the conductors | 2.315e-04 | 6.3e-04 to 7.6e-04 | 0.33× |
| inside them | 4.143e-02 | — | upper bound, log kernel |

**THE GLOBAL IRREDUCIBLE ERROR IS THE CONDUCTOR MODEL.** A figure that refused
to refine across a 16× range in dofs is reproduced to **5 per cent** by two
Green's-function sums. M-111 inferred that from the error refusing to refine;
this measures it.

**OFF the conductors it explains only a third**, so about two thirds of that
floor is something else — regularity, the plasma edge, or the reference's own
accuracy. *The first version of the verdict in that script read only the
off-conductor row and concluded the opposite; both rows have to be read.*

## 0a. AND THE 2x2, WHICH MAKES THE BENCHMARK LIKE-FOR-LIKE TODAY

Ian: *"we should also allow for freegs4e's filaments so we can compare against
an equivalent model ( MEQ filament vs freegs filament + MEQ shaped vs freegs
shaped )"*.

**freegs4e ALREADY HAS BOTH SIDES.** `freegs4e/coil.py`'s `Coil` is the point
filament, `multi_coil.py`'s `MultiCoil` is explicit filaments — what freegsnke's
MAST-U uses — and **`shaped_coil.py`'s `ShapedCoil` is a polygon carrying a
UNIFORM CURRENT DENSITY**, integrated by `quadrature.polygon_quad` at 1, 3 or 6
points per triangle. That is `meq::Coil`'s model.

Measured on the same 18 DIII-D currents:

| | everywhere | off the conductors |
|---|---|---|
| filament against MEQ's rectangle — **what the benchmark does today** | 6.081e-03 | 2.315e-04 |
| **`ShapedCoil` against MEQ's rectangle — the same model both sides** | **6.580e-04** | **6.731e-05** |

**Matching the model is worth 9.2x globally and 3.4x off the conductors**, and
takes the conductor floor from ABOVE M-111's flat 5.785e-03 to nine times below
it — which is the difference between a benchmark that cannot see MEQ's
discretisation and one that can.

**AND THE SHAPED ARM NEEDS NO MEQ CHANGE.** It is a reference regeneration with
`ShapedCoil` in place of `Coil`, on freegs4e's side, and it can be done before
any of this plan is built. What is left, 6.580e-04, is the QUADRATURE alone —
6 points per triangle against MEQ's 24^2 tensor Gauss — and `ShapedCoil` caps
at 6, so that is the floor of a shaped-against-shaped race rather than something
to tune away.

**The FILAMENT arm is the one that needs MEQ**, and it falls out of CS-1 for
free: a filament is the degenerate case of `psi_c`'s quadrature, a single point,
so `meq::ConductorField` gives MEQ filament conductors as a by-product of being
built at all. That is an argument for CS-1's ordering rather than an extra
stage.

---

## 0b. THE OBJECTIVE IS ACCELERATION, NOT A CHANGE OF ANSWER

Ian, on reading the above: *"the objective of the coil subtraction is to
accelerate, not to give up solving for the exact psi including the plasma
effect"*.

That rules out the reading this plan first leaned on — "adopt freegs4e's
filament model and the floor goes away". **`psi_c` must be MEQ's own conductor
model evaluated exactly**, not somebody else's; the split is a change of
REPRESENTATION and must reproduce the same equilibrium. The benchmark agreement
that follows is a side effect of both codes then resolving their own conductors
exactly, not of MEQ adopting freegs4e's.

## 0c. AND THE VACUUM IS LINEAR, WHICH IS THE BIGGER HALF

Ian again: *"as we're outside the plasma boundary, that's linear, so we have a
way to cheat"*.

Outside the plasma `Delta* psi = 0`. The vacuum region carries no source, so
nothing there needs a non-linear solve and nothing there needs the resolution
the plasma does — it needs only to carry a harmonic extension. MEQ already
exploits this OUTSIDE `Gamma`, which is FB-5's exterior DtN coupling.

**AND THE "SHRINK GAMMA" VERSION OF THAT ARGUMENT IS WRONG. THIS SECTION SAID
2.6x AND THE GEOMETRY SAYS OTHERWISE.**

What it claimed: `Gamma` sits at 3.4 in `examples/mastu-nke.toml` only because
the conductors are meshed, so take them out and it moves to about 2.1 — 2.6x
less area. **Both halves are false**, measured off the config's own 23 blocks:

* `Gamma` must CONTAIN THE PLASMA, and MAST-U's limiter reaches `R = 1.80`,
  `|Z| = 1.51`, i.e. **radius 2.35**. 2.1 was never available at any conductor
  model.
* and MEQ's EXISTING exterior machinery — `ExteriorCoilSet`,
  `setExteriorConductors()`, which already takes a `meq::CurrentFilament` and
  already computes its field by Carlson's elliptic integrals — requires a
  conductor to CLEAR `Gamma`. At `Gamma = 2.4`, **2 of 23 conductors do.** The
  Solenoid is at 0.189 and PX at 1.05: a tokamak's coils hug the vessel and the
  plasma fills the vessel, so they interleave in radius rather than nesting.

| `Gamma` | conductors outside it |
|---|---|
| 2.0 (too small for the plasma) | 6 of 23 |
| **2.4 (the smallest that contains the plasma)** | **2 of 23** |
| 2.8 | 0 of 23 |

**So the free route does not exist and the honest win is smaller.** `Gamma` goes
3.4 -> about 2.5, which is **1.85x in area**, not 2.6x — and getting even that
requires the INTERIOR subtraction, because at any `Gamma` containing the plasma
almost every conductor is inside it. `setExteriorConductors()` handles two of
MAST-U's twenty-three and is not the answer.

**What survives of §0c**: the vacuum is still linear and the conductors still
need not be resolved by the mesh. What does not survive is the idea that MEQ's
existing exterior expansion already does it.

**What still has to be meshed**: the plasma AND the vacuum between its moving
edge and `Gamma`, because the free boundary moves during the solve and a DtN map
on a moving surface would have to be rebuilt every iteration. `Gamma` is fixed
and always contains the plasma, which is what makes its DtN a precompute.

---

---

## 1. The split

Write the total flux as

```
psi = psi_c + psi_p
```

where `psi_c` is the field of the CONDUCTORS alone, computed analytically, and
`psi_p` is what MEQ solves for. Since `Delta*` is linear and the conductor
currents are prescribed inputs of a forward solve,

```
Delta* psi_p = -mu0 r J_plasma( psi_c + psi_p )
```

with the coil term gone from the right-hand side entirely. **The coils need not
be in the mesh at all.**

`psi_c` is a Green's-function integral over each conductor. It needs no closed
form — a closed form for a rectangle is not required and probably not worth
chasing — because the integrand is smooth except for a log at the source point
and a tensor-product Gauss rule reaches machine precision on it cheaply. MEQ
already owns that quadrature in `meq::Coil`.

**AND IT IS COMPUTED ONCE PER MESH, NOT PER NEWTON STEP.** The conductor
currents do not move during a forward solve; only the profile scale does. So
`psi_c` at every quadrature point and every boundary node is a precompute, and
its cost is amortised over every iteration of every support sweep.

## 2. What it buys, in the order the evidence supports

**(a) THE DOMAIN SHRINKS, but by 1.85x and not 2.6x.** See §0c for the
correction and the geometry behind it: `Gamma` goes 3.4 -> about 2.5, capped by
the plasma's own 2.35 rather than by the conductors. On top of that, **1024 of 9361 triangles, 10.9%,
carry a coil attribute** on `examples/mastu-nke.msh`, and `CoilSize = 0.03`
against `Size = 0.30` forces a 10× graded refinement around each of them. The
two together are the case for the plan.

**(b) The conductors are resolved EXACTLY rather than to the mesh's order**, at
no mesh cost — which is where M-111's global 6.081e-03 goes. Not by adopting
anybody else's model (§0b) but by integrating MEQ's own rectangle to machine
precision instead of representing its field in a degree-k space.

**(c) Regularity, IF §0's experiment says so.** A jump in the source across a
rectangle's edge limits `psi`'s smoothness there, and a rectangle CORNER limits
it further; neither is something degree `k` recovers. Remove the coils from the
solved field and `psi_p`'s only irregularity is the plasma edge, which is
`CLAUDE_FB.md`'s `k <= j` result and is a separate known cap.

## 3. What it costs, and the two that are not obvious

* **The Dirichlet datum and the exterior coupling both move.** On `Gamma` the
  condition becomes `psi_p = psi_datum - psi_c`, and FB-5's DtN coupling then
  sees the plasma alone. That is arguably SIMPLER — the exterior problem's
  source is one connected region instead of a plasma plus twenty-three
  conductors — but it is a real change to `meq::ExteriorDtN`'s setup and to
  `setBoundaryData()`, and the two must move together or the datum is
  double-counted.
* **Every place that reads `psi` must read `psi_c + psi_p`.** The critical-point
  finder, the plasma-support fill, the limiter and X-point borders, the output
  writers, the flux-surface tracer. **This is the real cost of the plan**: the
  solved field stops being the physical field, and every consumer that forgets
  is a silent wrong answer rather than a failure. `CLAUDE.md`'s own record has
  three separate defects of exactly that shape — `psi*` against `psi_h`, the
  band continuation reading the foot, `scaledF` against `f`.
* **FB-7 already does a piece of this** — a conductor outside `Gamma` that the
  mesh does not carry — so the machinery for "a coil MEQ does not mesh" exists
  and this generalises it inward rather than inventing it.

## 4. What it does NOT fix

The plasma edge. `ConstrainPaxisIp` at `alpha_n = 1.2` makes the edge a
fractional power and caps the order at `k <= 1.2` whatever the coils do, so
**MAST-U cannot demonstrate this plan's benefit** and a smooth-edge case is
needed to test it. Nor does it touch anything in M-115, M-117 or M-119 — the
constraint discontinuities and the cold failures are about the border, not the
conductors.

## 4b. FILAMENTS FIRST, AND THE REASON IS NOT ONLY SPEED

Ian: *"we should definitely implement a filament model in MEQ as it may be
faster ( then solve for psi - psi_filaments, as that's nonsingular ), and an
error of 1e-2 may be acceptable ( to some people )"*.

**A FILAMENT CANNOT BE A DOMAIN SOURCE AT ALL, so filament support and the
subtraction are ONE piece of work rather than two.** A point source has no
finite-element representation as a current density; `psi` near it is
logarithmic. The only way MEQ can carry a filament is to subtract its known
field and solve for a remainder — which is precisely this plan. That is an
argument for building it, not a caveat.

**And the remainder is nonsingular, which is the whole trick.** `psi_c` carries
the log; `psi - psi_c` satisfies `Delta*( psi - psi_c ) = -mu0 r J_plasma`, whose
right-hand side is bounded and supported on the plasma alone.

**Cheaper than the rectangle, too.** `psi_c` for a filament is ONE Carlson
evaluation per point per conductor — `meq::Coil`'s own machinery, which
`Coils.hpp` records agreeing to 7.0e-13 — against a tensor quadrature for a
rectangle. So the filament case is both the easiest to build and the fastest to
evaluate, and the rectangle is the same code with a quadrature rule around it.

**THE ACCURACY IT COSTS IS MEASURED AND IS A PRODUCT DECISION.**
`conductor_model.py` on DIII-D: filament against rectangle is **6.081e-03**
globally and **4.143e-02 inside the conductors**, against **6.731e-05** off them
for the matched model. So a filament MEQ is a percent-level answer near the
coils and far better away from them — which, as Ian puts it, *may be acceptable
to some people*, and is a MODE to offer rather than a defect to hide. It must
say which model produced an answer, for the reason every other rung in this tree
must.

## 5. Staging

| | |
|---|---|
| **CS-0** | §0's two-Green's-function difference. **DONE** — the global irreducible error IS the conductor model, 1.05x |
| **CS-0b** | regenerate the references with freegs4e's `ShapedCoil`. **DONE** → **[M-139](MEASUREMENTS.md#m-139)**: DIII-D's benchmark floor goes **5.785e-03 → 7.7e-04, 7.5×**, and lands on §0a's predicted quadrature residual of 6.580e-04. **It still does not refine**, so the floor is lower and is still not MEQ's discretisation. §6 has what blocked it for a week |
| **CS-1** | `meq::ConductorField` for FILAMENTS FIRST — one Carlson evaluation per point, no quadrature, and the only conductor MEQ structurally cannot mesh. §4b. **BUILT**: `src/meq/ConductorField.{hpp,cpp}` and `tests/unit/ConductorFieldTests.cpp`, 7 cases, MFEM-free so CI runs it. §7 has the two decisions behind it and §7.3 what it does NOT yet do |
| **CS-1b** | the same with a quadrature rule around it, which is the rectangle |
| **CS-2** | the split on a FIXED-boundary case with coils, where nothing else moves |
| **CS-3** | the Dirichlet datum and the DtN coupling |
| **CS-4** | every consumer of `psi`, with a test per consumer that the total is read |
| **CS-5** | re-take M-111 with the conductor models matched. **[M-139](MEASUREMENTS.md#m-139) partly kills this as written** — the benchmark cannot resolve MEQ below about 7e-04 whatever either code does, so it cannot be the acceptance for a change whose whole claim is that the conductors are resolved EXACTLY. **The replacement is MEQ's own**: an expensive quadrature-based reference, §7 |

---

## 6. CS-0b, and what stood between the harness and its first run

**`race.py --shaped` WAS BUILT AND HAD NEVER BEEN RUN, AND IT COULD NOT HAVE
BEEN.** The flag, the `shaped_ref()`/`shaped_stem()` pair, the seven
`*_shaped.npz` references and the seven `machine-*-shaped.toml` configurations
all existed. What did not exist was the **initial guess** each of those
configurations names, and the failure is a long way from its cause.

**`.gitignore` ATE IT.** Lines 60–61 ignore `*.mesh` and `*.gf` across the whole
tree — correct, since both are ordinarily run output. `make_all.sh` writes each
machine's `<stem>-guess.mesh` and `<stem>-guess.gf` **beside** the TOML that
names them, from `mkcoldguess.py`, so the configuration and its guess are one
artefact; but only the TOML survives a commit. The **filament** set is tracked
because somebody `git add -f`'d it. The **shaped** set, generated later, was
written by the generator, swallowed by the ignore rule, and never committed.

**So the symptom appeared eight rungs later in a different tool**: every MEQ
rung of `race.py --shaped` reporting `FAILED`, with race.py's own failure filter
printing nothing because it greps for `MEQ:` lines containing *not converge* or
*error* and this was neither. Run by hand, the message is

    MFEM abort: Mesh file not found: examples/machine-f-diiid-shaped-guess.mesh

naming a file whose generator reports writing it. **A generated input that
`.gitignore` eats does not fail at generation**, and nothing between the two
points at the ignore rule.

**THE MEASUREMENT IT UNBLOCKED IS [M-139](MEASUREMENTS.md#m-139), AND §0a's
PREDICTION HELD.** Matching the model is worth **7.5×** on DIII-D's global
floor — 5.785e-03 to 7.7e-04 — and 2.0× off the coils. **And the error still
does not refine**: 8.484e-04 to 7.750e-04 over a 6.7× range in dofs and from
`k = 2` to `k = 3`, nine per cent and not monotone. §0a said what would be left
and put a number on it *before* the solve existed — the quadrature alone, MEQ's
24² tensor Gauss against `ShapedCoil`'s 6 points per triangle, at **6.580e-04**
— and the measured floor is **7.7e-04**. The reference is not the limiter: its
own 129²-against-257² self-difference is **9.864e-05**, an order below.

**SO CS-0b DOES NOT MAKE MEQ'S DISCRETISATION VISIBLE, AND THAT IS THE FINDING.**
It was expected to: §0a's closing sentence calls the shaped arm *"the difference
between a benchmark that cannot see MEQ's discretisation and one that can"*.
It lowers the floor by 7.5× and the floor is still a conductor artefact — a
**quadrature** one now rather than a **model** one, and one that cannot be tuned
away from either side, since `ShapedCoil` caps at 6 points. **A benchmark
against this reference can resolve MEQ no better than about 7e-04 whatever
either code does**, which is the honest statement of what this harness is for
and is worth knowing before CS-5 is planned around it.

**Fixed** by `!examples/machine-*-guess.mesh` and `!examples/machine-*-guess.gf`,
which covers both conductor models and both the existing and any future machine,
beside the exception `examples/limited-tokamak-guess.*` already had for the same
reason. The seven shaped guesses are now tracked.

**AND THE TRACKED SHAPED CONFIGURATIONS WERE STALE, THOUGH ONLY IN PROSE.**
Regenerating them changes nothing but a comment block — every number, including
`PsiAxis = 2.888899229e-01`, is reproduced to the last digit, which is the check
that the generator is deterministic. What the comment says, though, is worth
having: the tracked copy still described the guess as *"a workaround, not a
design"* and credited it to `mkexactguess.py`, and the generator now says
plainly that `mkcoldguess.py` builds it, that `mkexactguess.py` is *"the OTHER
one, which sums Green's functions over the reference's own converged Jtor"*, and
that **a run seeded from `mkexactguess` has the answer in the starting position
and is a different benchmark**. The filament set still carries the old comment
and wants the same regeneration, which is deferred only because
[M-137](MEASUREMENTS.md#m-137) was just measured on one of those files.

---

## 7. Two decisions taken, and what they settle

### 7.1 A node or vertex ON a filament is refused, and PROXIMITY is not

**The mesh is rejected if any node or vertex coincides with a filament.** That
is the whole of the policy, and the reason it can be that small is measured
rather than assumed: `filamentPsi()` is verified in `Coils.hpp` against an
independent transcription to **1.4e-12** relative, and walking in to the
conductor it keeps working where the textbook form does not —

| eps | `CurrentLoop` | `filamentPsi()` | `−½ ln eps` |
|---|---|---|---|
| 1e-07 | 8.0828343 | 8.0987690 | 8.0590478 |
| 1e-09 | **NaN** | 1.0401354e+01 | 1.0361633e+01 |
| 1e-13 | **NaN** | 1.5006924e+01 | 1.4966803e+01 |

sitting a constant 0.0397 above `−½ ln eps` at every one of the last rows, to
`k'² = 1e-300`. **So a point NEAR a filament is not a problem, and only a point
ON one is** — and a large value at one evaluation point is not a large error,
because nothing here approximates `psi_c` by a polynomial. `psi_c` is
*evaluated*; what is approximated is `psi_p`, which is smooth **because** `psi_c`
carries the whole logarithm.

**So the tolerance is about coincidence and not about conditioning**, and it is
set to catch "the mesher put a node here" — a coordinate that came from the same
double, possibly through a text round trip — rather than to enforce a clearance.
A clearance rule would be a different and much stronger claim, and nothing
measured supports needing one.

**Why refuse at all rather than perturb or take the finite part.** A perturbed
node changes the geometry the solve reports without saying so, and a finite part
is a different field from the one `psi = psi_c + psi_p` names. Both make a run
that looks like it worked. `CLAUDE.md`'s standing preference is a refusal that
names the cause, and a coincident node is a thing the user chose — `[[coils]]`
drives the mesher, so a filament at a meshed coil's centre is reachable by
writing the obvious file.

### 7.2 CS-5's acceptance is MEQ's own quadrature, not `freegs4e`

**CS-5 as staged cannot work and [M-139](MEASUREMENTS.md#m-139) is why.** The
benchmark's floor after matching the conductor models is `7.7e-04` and it does
not refine, because `freegs4e.shaped_coil.ShapedCoil` caps at 6 points per
triangle. **A change whose entire claim is that the conductors are resolved
exactly cannot be accepted against a reference that resolves them to 6 points.**

**The replacement is an expensive quadrature-based check of MEQ's own.** MEQ
already owns the kernel — `coilPsi()` at whatever order is asked for, and
`filamentPsi()` exactly — so a reference field can be built at a quadrature order
far beyond what a solve would use, on the same conductors, and the subtraction
measured against it. That is a self-comparison and it is the right one here: the
claim is about representation, not about physics, and the two arms must reach the
**same equilibrium** by construction (§0b) rather than merely agree.

**The second acceptance is the one with teeth, and it needs no reference at
all**: the same machine solved with a conductor MESHED and with it SUBTRACTED
must agree to the discretisation. That is what says the split changed the
representation and not the answer, and it is the property §0b exists to protect.

### 7.3 What CS-1 is, and the three things it is not

**`meq::ConductorField` is `psi_c` and nothing else**: a set of
`meq::CurrentFilament`, `psi`/`gradPsi`/`flux` summed over them, `totalCurrent()`
for the boundary-integral check, and `coincides()`/`indexAt()` for §7.1's
refusal. Seven unit cases, three of them asserting at **zero tolerance** rather
than a small one — superposition, `psi_c ≡ 0` on the axis, and the flux
identity — because all three are identities rather than approximations.

**IT IS MFEM-FREE AND THAT DECIDED THE REFUSAL'S SHAPE.** It sits beside
`meq::Coils` in the half of `src/meq` CI can build, so it cannot take an
`mfem::Mesh`; the check is a per-**point** predicate the caller loops over the
mesh's vertices and nodes. That is not a limitation worked around — it keeps the
physics unit-testable where a convergence study cannot run, and puts the mesh
walking where meshes already live.

**`flux()` is the gradient of the sum divided once, not a sum of quotients.**
The same number in exact arithmetic; one division rather than `N`, and NaN once
on the axis rather than per filament. Commented at the site, because it is
exactly the kind of thing a later rewrite tidies back.

Three things it deliberately does **not** do, each belonging to a later stage:

| | |
|---|---|
| **no caching** | §1 wants `psi_c` precomputed once per mesh rather than per Newton step, and that is right — but there is no caller yet, and a cache built before its access pattern is known is a cache built against a guess. **CS-2** |
| **no rectangles** | `CS-1b`, and it is the same code with a quadrature rule around it: `meq::CoilSet` already answers `psi`/`gradPsi`/`flux` for MEQ's rectangles, so CS-1b is a second `add()` and a second loop rather than new physics |
| **nothing reads it** | no solver, datum, output writer or estimator consults it yet, which is **CS-2 to CS-4** and is where the plan's real cost lives — §3's *"every place that reads `psi` must read `psi_c + psi_p`"*, with a test per consumer |

**The refusal is written but not yet CALLED**, which is the one gap a reader
should not mistake for an oversight: `coincides()` exists and is tested, and the
site that loops it over a mesh arrives with the first caller in CS-2. Until then
a coincident node is still caught — by `meq::filamentPsi()`'s own throw, at the
first evaluation instead of at setup. **The early check is a better message, not
a missing guard.**
