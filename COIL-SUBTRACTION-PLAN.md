# The coils out of the mesh: `psi = psi_coils_analytic + psi_solved`

Ian's idea, 2026-09-18, in his words: *"possibly we could have inside-coil
approximation that uses the fact that the near field of a rectangular coil is
solvable? approximate psi_coil as psi_analytic_coil + polynomial"*, against the
requirement that *"we should be faster at 1e-3 and 1e-4 not just at very high
accuracy, high order should allow us to not use very many elements"*.

`CLAUDE.md` is the index and `CLAUDE_FB.md` owns the free boundary; this defers
to both. It is a DESIGN, not a measurement — nothing here has been built.

---

## 0a-pre. A STANDING REQUIREMENT: THE MESHED FINITE COIL NEVER STOPS WORKING

**MEQ must always be able to switch back to solving finite-sized coils
accurately, and the split is an OPTION rather than a replacement.** This is a
constraint on every stage below and it outranks any of them.

**Why it has to be said rather than assumed.** A plan whose whole content is
"take the conductors out of the mesh" drifts naturally towards making the other
route second class — the meshed source stops being exercised, then stops being
measured, then stops working, and nobody notices until a case needs it. And
cases will: a finite coil carried as a domain source is the route that needs no
`psi_c` evaluation per quadrature point, that puts the conductor's own current
into the interior equation where a force calculation can see it, and that has
been MEQ's answer since FB-2.

**What it forbids, concretely.**

* `meq::CoilSet::f()` and `meq::CoilAugmentedSource` stay, stay tested, and stay
  the default for a rectangle. The split must be asked for.
* **No stage may make the meshed route pay for the split.** Today it does not:
  `setConductorField()` is off by default and every shift in the tree is guarded
  on a null pointer, so a run that does not ask for it is **bit-identical** —
  which CS-2's control case asserts rather than assumes.
* CS-6's output format must be able to say a run used **either**, and a `.nc`
  that regenerates its input regenerates whichever was asked for.

**AND IT IS WHAT MAKES §7.2's SECOND ACCEPTANCE THE IMPORTANT ONE.** The same
machine solved with a conductor **meshed** and with it **subtracted** must agree
to the discretisation. That comparison is only available while both routes work,
and it is the only check that says the split is a change of REPRESENTATION and
not of answer — §0b's requirement. **So the two routes are not rivals; the
meshed one is the subtraction's instrument.**

**A rectangle is the conductor that can go either way**, which is what makes the
comparison possible at all — CS-1b put rectangles in `meq::ConductorField` for
that reason as much as for completeness. A filament can only be subtracted, so
it can never be the cross-check.

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
| **CS-1b** | the same with a quadrature rule around it, which is the rectangle. **BUILT** — `add( Coil const & )` beside `add( CurrentFilament const & )`, summing through `meq::CoilSet`. **`coincides()` stays FILAMENT-ONLY**, which is a contract and not an omission: `meq::coilPsi()` is *"Valid EVERYWHERE, including inside the coil"*, so a rectangle has no line singularity for a mesh point to land on and refusing one would reject the ordinary configuration. The quadrature order forwards, which is what §7.2's CS-5 replacement needs |
| **CS-2** | the split on a FIXED-boundary case with coils, where nothing else moves. **BUILT AND GREEN**, and the acceptance is an IDENTITY rather than a rate — §7.4 |
| **CS-3** | the Dirichlet datum and the DtN coupling. **BUILT AND GREEN** — §10. Two halves that had to land together: `prepare()` transfers `g − psi_c` inward, and `transmissionConstraint()` ADDS the subtracted conductors' moment on Gamma where FB-7's is subtracted. The refusal in `setConductorField()` is lifted and replaced by a GEOMETRIC one — a subtracted conductor must lie strictly INSIDE Gamma, the mirror image of `setExteriorConductors()`' own precondition |
| **CS-4** | every consumer of `psi`, with a test per consumer that the total is read. **DONE** — the critical-point finder, the element fill, `peakAt`, the source and its Jacobian, the limiter search and value. **AND THE `.nc` GRID**, which §8.3 says is one of the two formats that CAN keep `psi_c` exact: the driver adds it back at every located node, and `B` with it through `meq::ConductorField::poloidalField()` — the one entry point that takes the axis limit, where `q = ( 1/r ) grad_bar psi` is `0/0` and a half-disc machine's whole first grid column sits. **AND TRANCHE TWO CLOSES IT** — `meq::ContourTracer` takes the solver's conductors in its constructor and shifts at its own seam, so `_surfaces.nc`, the flux-surface averages and the `(Psi, theta)` family are level sets of the physical flux. §12 |
| **CS-T** | **`[conductors] Model`, the key that makes any of this reachable from a file.** **BUILT** — §11. `"meshed"` (the default, unchanged), `"subtracted"`, `"filament"`, plus `QuadratureOrder` for the rectangles. The driver drops its `meq::CoilSet` when the split is taken, because the double count is the failure with no symptom |
| **CS-6** | **a restart format that self-describes.** §8.3: a `.gf` cannot say whether it holds `psi` or `psi − psi_c`, and a flag would not be enough because a remainder is only meaningful with the conductors it is a remainder from. NetCDF is the vehicle MEQ already has — §9 |
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

### 7.4 CS-2, and why its acceptance is an identity rather than a rate

**`setConductorField()` IS THE OPPOSITE CASE TO `setExteriorConductors()`, AND
THE ALGEBRA IS THE WHOLE DIFFERENCE.** FB-7's call needs its conductors
**outside** `Γ`, because there `psi_coil` is `Δ*`-harmonic in `Ω` and the
conductor can enter through the boundary alone — its own documentation says a
conductor inside `Ω` "belongs in the SOURCE... put it here and the interior
equation silently loses it". CS-2 is exactly the inside case, and what makes it
work is that `psi_c` is **not** harmonic there:

    Δ* psi_p = Δ* psi − Δ* psi_c
             = −mu0 r ( J_plasma + J_coil ) − ( −mu0 r J_coil )
             = −mu0 r J_plasma

**the conductor's delta cancels exactly.** So the conductor is in neither the
mesh nor the source, and the remainder sees a bounded right-hand side supported
on the plasma alone.

**WHAT MOVED IN THE SOLVER IS ONE COEFFICIENT AND ONE REFUSAL.**
`setBoundaryData()` keeps its meaning — it is the datum for the **physical**
`psi` — and what is projected onto `Γ` becomes `psi|_Γ − psi_c|_Γ`, through a
`RemainderDatumCoefficient` at the one existing projection site. Plus the §7.1
refusal, walked over the mesh's **vertices and nodes**, which are different sets:
a curved or high-order mesh carries nodes the vertex array does not, and it is
the nodes the datum and the output grid are evaluated at.

**AND `setConductorField()` REFUSES A CONDUCTOR THAT IS NOT INSIDE `Γ`** when a
DtN is installed — §10, which is where the coupling's own two halves are. On a
**fixed**-boundary problem there is no `Γ` and no series, so a subtracted
conductor may sit anywhere at all and no geometric refusal applies; that is
this section's own domain.

**THE ACCEPTANCE IS `0.000e+00` AND NOT A RATE.** Put a filament inside `Ω`,
give the solve no plasma, and hand it the filament's own field as the physical
datum. Then `Δ* psi_p = 0` with `psi_p|_Γ = 0`, so the remainder is identically
zero and the reported total is the filament's field **exactly** — at every
point, including arbitrarily close to the conductor, which no meshed solve of
any degree can do. That is §2(b)'s claim reduced to something a zero-tolerance
assertion can check:

| `k` | elements | `max |psi_p|` | `max |psi_c|` on `Γ` |
|---|---|---|---|
| 1 | 128 | **0.000e+00** | 1.757e-01 |
| 2 | 128 | **0.000e+00** | 1.757e-01 |
| 3 | 128 | **0.000e+00** | 1.757e-01 |

**AND IT DISCRIMINATES, WHICH IS THE HALF THAT MAKES IT WORTH HAVING.** The same
configuration with `setConductorField()` simply not called gives
`max |psi_h| = 7.854e-02` — 45% of the datum's own scale. It converges, it is
smooth, and it is the harmonic extension of `psi_c|_Γ` rather than zero: a
perfectly plausible solve of the wrong problem. Without that control the
identity above would also be satisfied by a solver that ignored the datum
entirely.

`tests/convergence/ConductorSubtraction.cpp` is the acceptance, four cases.

### 7.5 Writing the doxygen found two defects, which is the argument for writing it

Annotating the accessors whose MEANING changes under the split — the task being
"so we know these meanings, this could bite later" — required saying what each
one returns, and two of those sentences turned out to be false of the code.

**THE SOURCE WAS EVALUATED AT THE REMAINDER.**
`SourceIntegrator::sourceValue()` calls `source->f( r, z, psi )` with the solved
potential, which under the split is `psi_p` — while `J_plasma` is a function of
the **physical** flux. `f( r, z, psi_p )` where `f( r, z, psi_p + psi_c )` is
meant converges at the full rate to a different equilibrium. Fixed by
`SourceIntegrator::setConductorField()`, on the existing `setPlasmaComponent()`
pattern; null shifts by exactly zero, so nothing that does not use the split
moves by a bit.

**AND THE JACOBIAN TOOK THE WRONG ARGUMENT, WITH A PLAUSIBLE REASON FOR IT.**
The first version of that doxygen said `dFdPsi` needs no shift because
`d( psi_c + psi_p )/d( psi_p )` is exactly 1. **True of the factor and
irrelevant to the argument**: `dFdPsi` is *evaluated at* a flux, and
`dF/d(psi_p)` at `psi_p` is `dF/dpsi` at `psi_c + psi_p`. The assembly was
passing the remainder. This is `CLAUDE_HDGGS.md`'s *A wrong Jacobian is
invisible to a convergence table* in its exact form — Newton converges either
way, to a different problem or merely more slowly.

**NEITHER COULD HAVE BEEN CAUGHT BY §7.4's ACCEPTANCE**, and that is the
transferable part. Its source is a vacuum, `f ≡ 0`, so it returns the same
number whichever flux it is handed: **a test that cannot fail on the thing being
changed**, which is the same shape as `Contour::fallbackLocations` reading 0
before and after, and as the dead `of which transmission` counter. The case
added for it compares the split against a control that does the shift **in the
caller** — a source wrapper plus a caller-shifted datum, with no conductor field
set — so the control runs none of the new code:

    |split - control| = 8.979e-15   against a scale of 5.590e-02

with the control's own scale asserted non-trivial, since a comparison between
two ways of computing nothing would also agree.

**AND ONE LIMIT IS REFUSED RATHER THAN DOCUMENTED AWAY.**
`setConductorField()` rejects a `meq::NormalisedSource`: that machinery divides
by `psi_ax` and `psi_bnd`, which are functionals of the physical flux, while
this solver's critical-point search runs on `psi_p`. The normalisation would be
taken against the remainder's axis — not a failure, a different problem solved
perfectly. **CS-4 is the stage that makes every consumer read the total**, and
until then the refusal is what stops the split being used where it is wrong.

---

## 8. CS-4 scoped: the inventory, the design, and the one place exactness cannot survive

**CS-4 IS THE REAL COST OF THIS PLAN AND §3 SAID SO.** The solved field stops
being the physical field, and every consumer that forgets is a silent wrong
answer rather than a failure. This section is what makes it safe to start.

### 8.1 The inventory, counted rather than remembered

| where | reads | what |
|---|---|---|
| `apps/meq.cpp` | **19** | `potential()` ×8, `flux()` ×7, `postProcessedPotential()` ×4 — the four output formats |
| `src/meq/GradShafranov.cpp` | **40** | `potentialGf`, `darcyFlux`, `traceGf` internally: the bordered rows, the support fill, the limiter and X-point |
| `src/meq/FluxSurfaces.cpp` | 4 | the tracer and the surface averages |
| `src/meq/Estimator.cpp` | 4 | the residual estimator, all four fields |
| `src/meq/CriticalPoints.cpp` | 1 | the axis and X-point search |

**About 68 sites.** The count matters because it decides the design: a change
that has to be made correctly 68 times will not be.

### 8.2 The design: hand each consumer `psi_c`, do not pre-add it

**The rejected option is to form a total `GridFunction` once and let consumers
read that.** It is the smallest diff and it throws away the plan's central
claim: `psi_c` projected onto the discrete space is `psi_c` to the mesh's order,
which is exactly the error §2(b) exists to remove. It would also reintroduce the
singularity into a finite-element function, which for a filament is the thing
that cannot be represented at all.

**So each consumer takes a `ConductorField const *` and adds `psi_c` at its own
evaluation points**, which is the pattern `SourceIntegrator::setConductorField()`
already follows — null shifts by exactly zero, so no existing path moves by a
bit. Consumers that evaluate pointwise keep the exactness; each is one setter
and one addition at the point where it already has `( r, z )` in hand.

### 8.3 AND THE `.gf` OUTPUT CANNOT KEEP IT, WHICH IS A FINDING RATHER THAN A GAP

A `.gf` **is** a finite-element function — every `P_k` coefficient of a field in
a space — so "write the physical flux to a `.gf`" and "represent `psi_c` in the
discrete space" are the same request. For a filament that is not merely
inaccurate, it is impossible: `psi_c` is logarithmic at the conductor and no
polynomial space contains it.

**So the split's exactness reaches the outputs unevenly, and the difference is
structural:**

| format | keeps `psi_c` exact? | |
|---|---|---|
| `.nc` grid, `_surfaces.nc` | **yes** | sampled pointwise, so `psi_c` is evaluated rather than represented |
| flux surfaces, critical points, the tracer | **yes** | same reason |
| `.gf`, `.vtu` | **no** | both are the field's own coefficients; `psi_c` can only be projected |

**A `.gf` written under the split should therefore be the REMAINDER and say so
— and it CANNOT say so, which is the finding.** `mfem::GridFunction::Save()` is
`fes->Save( os )`, a newline, and the raw coefficients; the header carries
`FiniteElementCollection`, `VDim` and `Ordering` and nothing else, and the
loader reads the space and then the vector, so an extra line is consumed as
data. **There is no slot for "this is a difference" in the format**, and a file
that cannot distinguish `psi` from `psi − psi_c` is a restart waiting to read
the wrong field.

**AND A FLAG WOULD NOT BE ENOUGH EITHER, WHICH IS THE PART THAT DECIDES THE
DESIGN.** A remainder is only meaningful WITH the conductors it is a remainder
from: `psi_c` is not recoverable from the mesh, the space or the coefficients,
and no amount of provenance short of the conductor set reconstructs the physical
field. **A file holding a difference has to carry what it is a difference
from.** So the requirement is not a boolean attribute bolted to a `.gf`, it is a
container that holds data *and* metadata — which is CS-6.

### 8.4 What unlocks what

The refusals `setConductorField()` carried are the measure of CS-4's progress,
and both are now lifted:

* **a `NormalisedSource`** needs `psi_ax` and `psi_bnd` off the total, which is
  `CriticalPoints.cpp` plus the limiter and X-point rows — and `insidePlasma()`,
  which decides the support fill. That is the whole of what makes a *physical*
  case reachable, so it is CS-4's first tranche and everything else can wait.
* **an exterior coupling** was CS-3 and is a different question: the DtN's datum
  and transmission rows, not a consumer of `psi`. §10 is what it took, and the
  refusal it left behind is geometric rather than about progress.

### 8.5 The evaluation abstraction, which replaces §8.2's sixty-eight setters

**§8.2 PROPOSED HANDING EACH CONSUMER A `ConductorField const *` AND WAS THE
WRONG SHAPE.** It is correct and it is sixty-eight opportunities to forget, on a
defect whose failure mode is a plausible answer. `src/meq/FieldViews.hpp` is the
replacement: `meq::PotentialView` and `meq::FluxVectorView`, header-only.

**WHAT MAKES FORGETTING HARD IS THE NAMED CONSTRUCTOR, NOT THE DOCUMENTATION.**
There is no default constructor and no implicit conversion from a
`GridFunction`, so a call site writes one of

    PotentialView::physical( solved, conductors )   psi_c + psi_p
    PotentialView::remainder( solved )              psi_p alone, deliberately

and the choice is made in a word where it is made, rather than by knowing what
the class does. `remainder()` is spelled out because there **are** consumers
that want it — the residual estimator measures the discretisation of the problem
actually solved — and because a reviewer can grep for it.

**IT IS FREE WHEN THE SPLIT IS OFF, AND THAT IS ASSERTED RATHER THAN CLAIMED.**
Every method is inline and the conductor pointer is null on every existing path.
Measured over 128 elements at four quadrature points each:

| | error |
|---|---|
| `remainder()` and `physical( …, nullptr )` against the bare `GridFunction` | **0.000e+00** |
| `physical( …, &conductors )` against `GetValue` + `psi_c` at that point | **0.000e+00** |
| the two `value()` overloads against each other | **0.000e+00** |

against a field scale of 2.739e-01. **Bit equality is the point**: it is what
lets sixty-eight call sites migrate without any of them moving an existing
answer, so the migration can be done a file at a time and reviewed by its diff.

**AND IT IS WHERE §1's CACHE BELONGS WHEN IT IS BUILT.** `psi_c` is a Carlson
elliptic integral per filament and a cross-section quadrature per rectangle, at
every point a consumer asks about; §1 wants it once per mesh. That is now **one
type to change rather than sixty-eight call sites**, which is the second reason
for the abstraction and the one that was not the motivation.

**A TRAP IT ABSORBS ONCE INSTEAD OF SIXTY-EIGHT TIMES.** Adding `psi_c` needs
the point's `( r, z )`, which needs the element's transformation — and
`mfem::Mesh::GetElementTransformation( int )` hands out **shared scratch**,
which `CLAUDE.md` records as a silent wrong answer under threading and which
cost this project six call sites once already. The views take a caller-supplied
transformation where one is in hand and use a function-local `thread_local`
where it is not.

### 8.6 The first migration has a complication, and it is not the views' fault

`meq::CriticalPointFinder` is CS-4's first tranche — `psi_ax` is what lifts the
`NormalisedSource` refusal — and it does **not** read the field only through
`GetValue()`. Its screening passes read **raw dof coefficients**:
`potentialField( dofs[ i ] )` at `:564` and `:1131`, and
`fluxField( space.DofToVDof( i, 0 ) )` at `:268`.

**A dof coefficient of `psi_p` plus `psi_c` at that dof's node is not the
coefficient of the total**, and for a filament there is no coefficient of the
total at all. So those passes cannot be shifted the way an evaluation can, and
they need deciding one at a time:

* where the read is a **screen** — "which elements could hold a maximum" — it
  may be sound to screen on the remainder and evaluate the candidates on the
  total, but that is an argument about where extrema can be and it needs making
  rather than assuming;
* where it is a **value**, it has to become an evaluation.

**This is exactly the kind of place a hurried migration produces the silent
defect the abstraction exists to prevent**, so it is written down here rather
than attempted at the end of a long session.

---

## 9. CS-6: the exact field needs a container that carries metadata

**THE REQUIREMENT, IN ONE LINE: A FILE HOLDING A DIFFERENCE HAS TO CARRY WHAT IT
IS A DIFFERENCE FROM.** §8.3 establishes that `.gf` cannot — `GridFunction::Save`
writes the FE space header and the raw coefficients, the loader reads the space
and then the vector, and an extra line is consumed as data.

**AND THE METADATA IS NOT A FLAG.** `content = "remainder"` tells a reader it is
holding the wrong field; it does not let them fix it. `psi_c` is recoverable
from neither the mesh, nor the space, nor the coefficients — only from the
conductors. So the file must carry the **conductor set itself**: each filament's
`( r, z, I )`, each rectangle's `( centre, half-extents, I )`, and `mu0`, which
together with `meq::ConductorField` reconstruct `psi_c` exactly at any point.
That is a few dozen numbers beside a field of tens of thousands, so the cost is
nil and the alternative is a file nobody can safely read.

**NetCDF IS THE VEHICLE AND MEQ ALREADY HAS IT.** `src/meq/Output.cpp` writes
attributes through `putAtt` already — `axis_normalised_flux`, `psi_axis`,
`extrapolated_nodes` — so what is missing is not a dependency but a second kind
of NetCDF file: one holding **every `P_k` coefficient** rather than a sampled
grid. The existing `.nc` is deliberately lossy and is the interchange format;
this one is the *restart* format and must be exact.

| | today | under CS-6 |
|---|---|---|
| `.gf` | exact, and MFEM's own format, which GLVis reads | **kept, unchanged, and documented as the REMAINDER under a split** — it is MFEM's format and MEQ should not redefine it |
| the exact restart | the `.gf` pair | a NetCDF file: coefficients, the space description, `content`, the conductor set, `mu0` |
| `.nc` grid | sampled, lossy | unchanged, and it can write the **total** exactly, since it samples pointwise — one attribute says which it did |

**WHAT MAKES THIS MORE THAN BOOKKEEPING**: the `.nc` grid and the flux surfaces
keep the split's exactness because they sample, so under CS-6 the *lossy*
interchange format carries the physically exact field while the *exact* format
carries a remainder. That inversion is worth stating in `tools/README.md`, which
is the guide to which format goes with which reader, because it is the opposite
of what the names suggest.

**SETTLED: ONE NETCDF FILE CARRIES BOTH REPRESENTATIONS, AND THE `.gf` BECOMES A
DERIVED ARTEFACT FOR GLVis.**

| | |
|---|---|
| **the `.nc` carries two formats** | the **rasterization** it carries today — `psi` and `B` on a uniform `( R, Z )` grid, sampled, lossy, the interchange — **and the MFEM gridfunctions**: every `P_k` coefficient, on a `dof` dimension, with the space description in attributes. One file, two representations, and adding variables is backward compatible so existing readers are untouched |
| **the `.gf` is emitted by a small script** | from the `.nc`, carrying the **full `psi`** for GLVis. It stops being MEQ's restart format and becomes a *viewing* artefact — which is the right place for the lossy step, since for a filament `psi_c` has no representation in the space and a picture is what is wanted |
| **the input parameters go in the output** | so that the output file regenerates the input file |

**THE THIRD ONE ALSO SATISFIES THE SECOND HALF OF §8.3 BY ITSELF, WHICH IS THE
PART WORTH NOTICING.** §8.3 argues that a file holding a difference must carry
what it is a difference from — the conductor set, not a flag. **The conductors
ARE input parameters**, `[[coils]]` and whatever CS-6 adds for filaments, so an
output that regenerates its input carries them by construction. One mechanism,
not two, and the provenance requirement stops being special pleading for the
split and becomes an ordinary consequence of reproducibility.

### 9.1 What "the parameters that were read in" should mean, and it is two things

`meq::Configuration` **has no serialiser and does not keep the parsed table** —
it reads TOML into typed fields and the `toml::value` is gone. So there are two
different artefacts here and they answer different questions:

| | what it records | catches |
|---|---|---|
| **the input verbatim**, as a string attribute | the file as given | nothing, but it is trivially correct, costs a string, and makes the output self-contained **today** |
| **what MEQ actually used**, serialised from `Configuration`'s fields | the resolved configuration, defaults included | **a key accepted and ignored** — which this tree has met twice, `ProfileFile` and its kind on one source type, and `ConfineToPlasma` on the rotating source |

**Both, and they are not redundant.** The verbatim copy is provenance; the
serialised one is what the solve was. **And the pair supports a test with real
teeth**: parse a configuration, serialise it, parse the result, and require the
two `Configuration` objects to agree — which fails exactly when a key is read
into a field that nothing writes back, the accepted-and-ignored shape.

### 9.2 What is open, and it is the geometry rather than the fields

* **the mesh.** A `GridFunction` is meaningless without one, and GLVis needs
  `.mesh` beside `.gf` — so either the `.nc` embeds the mesh (vertices,
  connectivity, attributes, boundary) and the emitter script writes both, or the
  `.mesh` stays a separate file named in an attribute. **Embedding is what makes
  the file self-contained**, which is the whole design; it is also MEQ
  re-implementing a mesh serialiser, and `[mesh.generate]` already makes the
  mesh a build product of the TOML for the generated cases but not for the
  hand-made ones.
* **the profile tables.** `<stem>-pprime.dat` and `-ggprime.dat` are named by
  the configuration and read by `meq::SplineProfile`. An output that regenerates
  its input has to carry them too, or it regenerates a file that cannot be run.
* **whether the `.gf` pair survives at all** as MEQ's warm start.
  `[initialguess] File`/`MeshFile` reads `.gf` today, and a configuration naming
  one written under a split is the failure this stage exists to prevent.
  **REFUSING THE PAIR WAS PROPOSED HERE AND WAS DECIDED AGAINST.** Ian:
  *"assume the user knows what they're doing and so the warm start .gf file
  changes meaning (from psi to psi_p); warn loudly on stdout that you're doing
  it though."* So `apps/meq.cpp` prints the change of meaning whenever both are
  set, and the key is not refused. What CS-6 owes is the format that makes the
  warning unnecessary.


---

## 10. CS-3: the exterior coupling, and the two halves that had to land together

**WHAT MOVED IS TWO EXPRESSIONS, AND §3 NAMED THEM BOTH IN ADVANCE**: *"The
Dirichlet datum and the exterior coupling both move... the two must move
together or the datum is double-counted."* That is exactly what happened and
exactly what the acceptance falsifies.

| half | where | what it became |
|---|---|---|
| **Dirichlet** | `prepare()`, at the one site that builds the `PathTraceCoefficient` | what is transferred inward is `g( x ) − psi_c( x )` rather than `g( x )` |
| **Neumann** | `transmissionConstraint()` | the subtracted conductors' own moment on `Γ` is **ADDED**, where FB-7's exterior one is **SUBTRACTED** |

**THE SIGNS ARE OPPOSITE AND THAT IS A DERIVATION RATHER THAN AN ASYMMETRY TO
TIDY AWAY.** The condition is that the **physical** interior normal flux match
the exterior one. `row . state` is `∫_Γ ( q_h · ν ) C_m` over the **solved**
flux, which under the split is `q_p` while the physical flux is `q_p + q_c` — so
the interior side is short by the subtracted conductors' moment and it is put
back with a plus. FB-7's conductor is on the other side of the equation: its
field is part of the **exterior**, appears on both sides, and is removed from
the interior one with a minus. Same integral, opposite sign, and which it is
depends on **which side of `Γ` the conductor is on**.

**SO THE GEOMETRY IS NO LONGER OPTIONAL, AND THE OLD REFUSAL IS REPLACED BY A
BETTER ONE.** A subtracted conductor must lie strictly **inside** `Γ`; an
exterior one strictly **outside** it. `meq::ConductorField::containment()` is
the mirror image of `meq::ExteriorCoilSet::clearance()`, measured to the
**farthest** point of each conductor rather than to its centre, so one
straddling `Γ` is refused by both and belongs to neither. The reason is dual to
FB-7's: the exterior is a Gegenbauer series in `rho^( 1 − n )` standing for the
sources **within** `Γ`, and a subtracted conductor beyond it puts a singularity
in the region that series describes. **Both setters check the pair**, because
either may be called first.

**AND THE ESTIMATOR NEEDED THE SAME SHIFT.** `exteriorTransmissionResidual()`
compares an interior normal flux against an exterior one, and under the split
the interior one it reads off the solved field is `q_p`. Without `q_c` added,
`η₆` reports a perfectly matched boundary as mismatched by exactly the
subtracted conductors' own flux — and the adaptive loop **marks on it**,
refining `Γ_h` to chase an error that is not there.

### 10.1 The acceptance is FB-7's case turned inside out

`aConductorOutsideGammaReachesTheCoupledSolve` puts a conductor **outside** `Γ`
and the continuous answer is `a = 0`. CS-3 puts one **inside** and every part of
that inverts: `psi_p` is identically zero and `a` is **not** zero — it is
`psi_c`'s own Gegenbauer trace on `Γ`, the whole datum rather than an error.

`theSplitReachesTheExteriorCoupling`, no plasma, twelve modes, conductors at
`r = 0.50`, `z = ±0.20` inside `Γ = 1.5`:

| `n` | Newton | `max |psi_p|` | modal sum `−` `psi_c` on `Γ` |
|---|---|---|---|
| 12 | 1 | 8.7964e-09 | 1.2045e-08 |
| 24 | 1 | 1.4261e-08 | 1.2099e-08 |
| 48 | 1 | 1.4687e-08 | 1.1955e-08 |

`psi_c` on `Γ` is **6.3838e-02**, so the remainder is **2.3e-07 of the field it
is a remainder from** and the modal sum reproduces `psi_c` to **1.9e-07** of it.

**IT IS FLAT IN `h`, AND THAT IS THE RIGHT SHAPE RATHER THAN A DISAPPOINTMENT.**
`psi_p` is identically zero in the continuum and the discrete problem has
nothing to converge — the remainder equation's right-hand side is zero and its
boundary data is zero to the series truncation. What is left is the twelve-mode
residual plus round-off, neither of which depends on `h`. **A rate assertion
here would be wrong in both directions**: it would fail on correct code, and it
would pass on code merely converging towards the right answer from somewhere
large. §7.4 makes the same choice for the same reason.

**The conductor's distance from `Γ` is the experiment's precision.** The
truncated series cannot represent what falls like `( reach/rhoGamma )^n`, so the
fixture puts the conductors at 0.40 of `Γ` and twelve modes leave about `1e-08`
— below the round-off rather than beside it. At `r = 0.9` the ratio is 0.73, the
residual is 15%, and `psi_p` would be carrying the **series** rather than the
solve.

### 10.2 Both halves are falsified, not argued

→ **[M-143](MEASUREMENTS.md#m-143)** — the three arms at `n = 48` · what each
broken half leaves behind · why the gate is `1e-5` and not something loose

**NEITHER BROKEN ARM DIVERGES OR FAILS TO CONVERGE.** Both take one Newton step
and report a smooth, plausible field — the same disguise FB-7 records for the
transmission row's own sign. A residual check could not tell them from the
correct code; only a comparison against the closed form can.

---

## 11. CS-T: `[conductors] Model`, and the double count

**Nothing in this plan was reachable from a file**, so the only thing that could
run it was a convergence test — which is why *"at what point can we start
testing if this pathway is useful?"* had no answer. One key decides it:

| `Model` | | |
|---|---|---|
| `"meshed"` | **the default** | a rectangle carrying a uniform current density as a domain source, which the mesh must resolve. Unchanged, and §0a-pre is why the default is load bearing rather than a convenience |
| `"subtracted"` | | the same rectangle, integrated analytically and taken out of the mesh |
| `"filament"` | | a point at each rectangle's centre carrying the same total current — the only conductor MEQ structurally cannot mesh (§4b), and a **different machine** rather than the same one computed differently |

`QuadratureOrder` is the rectangles' and only theirs; it is refused on the other
two models rather than accepted and ignored, as is `Model` naming anything but
`"meshed"` with no `[[coils]]` to act on.

**THE DOUBLE COUNT IS THE FAILURE WITH NO SYMPTOM, AND IT IS WHAT THE DRIVER
ACCEPTANCE IS BUILT AROUND.** Leaving `meq::CoilAugmentedSource` in place beside
`psi_c` puts the same amperes into the equation **twice**: the run converges,
closes every border, and reports the file's own `coil_current` while carrying
double it. `coils.reset()` in the driver is what stops it, and the control is
`psi_p + 2 psi_c` rather than a tolerance.

→ **[M-144](MEASUREMENTS.md#m-144)** — the three figures · why the middle one is
the MESHED arm's error · the refinement study behind the gate

**THREE DIAGNOSTICS ARE ABOUT WHERE THE COPPER IS RATHER THAN HOW IT ENTERS
`F`**, and reading the driver's `meq::CoilSet` for them would have switched all
three off under the new key: the axis screen, `checkAxis()`'s exclusion, and the
axis-inside-copper warning. *"A plasma has no magnetic axis inside a conductor"*
is a statement about the machine, and `psi_c` has its own O-point in the middle
of each rectangle exactly as the meshed source does.

**AND THE `.nc` GRID CARRIES THE PHYSICAL `psi` UNDER BOTH ROUTES**, which is
§8.3's division made real: the driver adds `psi_c` back at every located node,
and `B` with it. That needed the **axis**, where `q = ( 1/r ) grad_bar psi` is
`0/0` and a half-disc machine's whole first grid column sits —
`meq::filamentAxisFlux()`, `meq::coilAxisFlux()` and
`meq::ConductorField::poloidalField()` are the closed-form limit, and
`flux()` keeps its NaN deliberately.

### 11.1 And the filament model's error profile falls out, agreeing with §4b

§4b measured filament-against-rectangle on DIII-D's eighteen conductors by two
Green's-function sums with no solver anywhere: **6.081e-03** globally and
**6.731e-05** off them. The key makes the same comparison available through a
whole MEQ solve, on a different machine at a different scale, for the cost of
one edit — and it lands in the same place.

→ **[M-145](MEASUREMENTS.md#m-145)** — the profile against distance from a
conductor · why the global figure is 1.23 and why that is the model working

**THE ONE THING WORTH SAYING OUT LOUD**: `examples/coils-rectangle.toml` is
fixed boundary on a mesh that IS the plasma, so its conductors sit *inside* the
sampled domain and the output grid has nodes millimetres from each filament.
That is the **worst case** for a filament model, not a typical one. A real
machine has its coils in a vacuum region outside the plasma.

### 11.2 And on `examples/coils-rectangle.toml` the split buys exactly 1.00×

Worth saying before somebody times it. That file's mesh is a plain Cartesian box
that was never graded to its coils — its own header notes about three cells per
coil are **cut** — so there is no coil refinement for the split to remove and
the element count is identical either way. M-142's **1.96×** lives on a mesh
built **around** its conductors, and this tree has no such case on a fixed
boundary. `TODO` carries that as its own entry: half a dozen serious
machine-relevant fixed-boundary cases, meshed to their conductors, is what turns
this key into a measurement.

---

## 12. CS-4 tranche two: the tracer, and the curve that closed cleanly

**THE LAST PLACE THE SPLIT COULD PRODUCE A WRONG ANSWER RATHER THAN A REFUSAL.**
`meq::ContourTracer` roots the field it is handed, so under `[conductors]` it
rooted `psi_p` — and a level set of the remainder is a smooth closed curve that
is not a flux surface. Everything downstream inherited it: `V'`, the safety
factor, the metric, `_surfaces.nc`'s every column, the `(Ψ, θ)` family, all
computed correctly over the wrong curve.

**IT HAD NO SYMPTOM, AND THE MEASUREMENT IS THAT RATHER THAN THE FIX.**
→ **[M-146](MEASUREMENTS.md#m-146)**. Asked for *"the surface through this
point"*, the conductor-blind tracer returns **89 points, closed**, at the same
corrector tolerance as the real one, with an ordinary turning number and no
stalled corrections — and the physical flux varies by **4.19e-02** along it
against a level of 8.39e-02. Half the level, in a file whose only per-node mask
is about the band.

### 12.1 One seam, because the class already had one

`sampleField()` is documented in `FluxSurfaces.hpp` as *"THE SEAM. The only
place psi and q are read at a physical point"* — written for IN-0's band, and it
paid for itself a second time here. The shift is one addition there and the
four public entry points, `meq::surfaceAverages`, `meq::extractFluxSurfaces` and
`AngleParametrisation` all inherit it without an edit each. **That is §8.5's
argument arriving at its second customer**: the alternative was one shift per
consumer on a defect whose failure mode is a plausible curve.

**THE CONSTRUCTOR TAKES THE CONDUCTORS RATHER THAN ASKING FOR THEM.** A caller
holding a `GradShafranovSolver` is holding something that knows whether its
field is the whole flux or a remainder, so `ContourTracer( solver )` reads
`solver.conductorField()` itself and `apps/meq.cpp` needed **no edit at all**.
The bare-field constructor cannot know, which is the only reason
`setConductorField()` is public.

### 12.2 Three details that are not the obvious ones

| | |
|---|---|
| **`q_c` comes through `poloidalField()`, not `flux()`** | `q = ( 1/r ) grad_bar psi` is `0/0` on `r = 0`, where `flux()` returns NaN deliberately. `poloidalField()` is the entry point that takes the closed-form limit, and `B_R = −q_z`, `B_Z = +q_r` inverts to what the sample wants. Off the axis the two are the same numbers through two sign flips, so nothing there moves by a bit — but a half-disc machine's contour can reach `r = 0` and a NaN would end the trace with `LeftMesh` |
| **`faceJump()` keeps reading the REMAINDER, deliberately** | it measures the DG discontinuity of `psi_h` across a face, and `psi_c` is analytic: it takes the same value on both sides and cancels exactly. Adding it would cost two Carlson evaluations per face crossing to subtract two equal numbers, in front of a jump this project measures down to 6.8e-10 |
| **the band datum is PHYSICAL and the lift is of the remainder** | `setBandExtension()`'s `g` keeps the meaning its name has, so its default of zero stays right for a fixed-boundary problem either way. `extendField()` subtracts `psi_c` at the foot and the seam adds it at the point, and `psi_c( x ) − psi_c( xbar )` **is** the line integral of `r q_c` along that path — exactly, no quadrature, and the conductor's logarithm never enters the lift |

### 12.3 And the corrector's scale had to move with it

`potentialScale()` multiplies `tolerance` to give an **absolute** residual
target, and the residual the corrector now stops on is one of `psi_p + psi_c`.
A vacuum remainder can sit six orders below the flux that is physically there —
§10 measures `max | psi_p |` at 1.4e-08 against a datum of 6.4e-02 on Γ — so
scaling by the remainder alone asks for a relative accuracy below anything a
discontinuous `psi_h` can offer anywhere.

**It would not give a wrong answer; it would stall every point and say so.**
That is the milder failure, and it is still a failure of the SCALE rather than
of the field, which is the kind that gets diagnosed as the tracer being broken.
`conductorScale()` is `max | psi_c |` over the mesh's **vertices**, cached per
tracer, and it is added: `| psi_p | + | psi_c |` bounds `| psi_p + psi_c |`, and
a bound is what a scale wants to be. Exactly zero without conductors, so every
existing path keeps the value it had.
