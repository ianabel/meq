# Where MEQ is, and what to do next

Written 2026-08-26, substantially revised 2026-08-27, 2026-08-29, 2026-09-01 and
2026-09-03. **The operational record is authoritative on anything technical and
is now four files**: `CLAUDE.md` is the index and carries the build, the
commands, the traps and the layout, with `CLAUDE_HDGGS.md` (the equation, the
discretisation, the solve, toroidal flow, the linear solves), `CLAUDE_FB.md`
(free boundary and the `freegs4e` benchmark) and `CLAUDE_INVERSION.md` (the
flux-surface work) beside it. `TODO` holds work that is understood but not
scheduled.

**This file is only about order** — what to do first, what waits on what, and
what is deliberately not being done yet. **Item numbers are cited from `TODO`,
from the plan files and from the CLAUDE files, so they do not get renumbered**; a
closed item becomes a marker rather than being removed.

The plan files. **Three are gone** under the standing rule that a plan with
nothing left in it is not a plan. `DRIVER-PLAN.md` and `FLOW-PLAN.md` were
converted to `docs/` on 2026-09-06, because everything they staged is built and
what a reader needs is the manual. `INTERPOLATORY-HDG-DECISION.md` went on
2026-09-19 for a different reason — it asked *where the code should live*,
the answer was acted on, and upstream then built the half it assigned upstream —
so there was nothing to convert and its two live questions moved into the plan
that succeeded it. Git has all three.

**A plan that is finished and a plan that is superseded both get closed, and
they are not the same thing.** The first becomes documentation; the second
becomes a pointer. What neither becomes is an entry that still reads like work.

| | |
|---|---|
| ~~`DRIVER-PLAN.md`~~ | stage 7 — **done**, and now `docs/running.rst`, `docs/output.rst` and `docs/configuration.rst` |
| ~~`FLOW-PLAN.md`~~ | item 9, FL-0 to FL-8 — **done**, and now `docs/rotation.rst`, which carries the derivation `RotatingSource.hpp` defers to |
| ~~`INTERPOLATORY-HDG-DECISION.md`~~ | where CCSZ-I's code lives — **closed 2026-09-19**, superseded by `INTERPOLATORY-HDG-PLAN.md`. Its answer was *split, with the element-local algebra upstream*; both asks were filed and **upstream built them**, as `gf-interp-hdg-dev`. Its two surviving open questions — Remark 2.2's `ℓ`, and whether `HDGPotentialPostprocessor` is reentrant under a threaded element loop — are carried into that plan's §7.4. Git has it |
| `INVERSION-PLAN.md` | item 10's machinery — **CLOSED 2026-09-19**, IN-A to IN-P, 9 of 9, with the ten acceptance tests verified passing together. **It stays rather than being converted** for a reason the other two did not have: the user-facing half already went to `docs/flux_surfaces.rst`, `docs/surface_geometry.rst` and `docs/normalised_flux.rst`, and what is left is the derivation, which `src/meq`, `tests/` and `apps/` cite by section and stage number at **125 sites across 32 files**. It is a design record, not a plan. **The one gap closing it found is now half fixed**: the driver set `options.toroidalField` nowhere, so `q` reached no output at all; the `q`-driven route now reports it to 4.5e-07 of target (M-131), and giving a prescribed-field run a `q` is a schema decision left open |
| `INTERPOLATORY-HDG-PLAN.md` | item 12.1, and it **supersedes §12.1 below**. MFEM has already built the method and MEQ links it; this is a wiring job plus two small asks. **Read its §0 before its stages** — the falsifying arithmetic deflates the headline to 1.6–2.0× fewer `F` evaluations against an Amdahl ceiling of 23–28% |
| `THREADING-PLAN.md` | **CLOSED.** Items A, B, C and D built (`0108f49`), **E and F declined standing** — if either grows it is a new plan, and that file's *When E or F would be worth revisiting* says what "grows" has to mean, since a share inflates whenever anything else is threaded — G moved to `BORDERED-GLOBALISATION-PLAN.md` §12, §0's two sub-slices read and the payoff measured at **1.16×** → **[M-137](MEASUREMENTS.md#m-137)**. The budget it was written against is superseded by **[M-138](MEASUREMENTS.md#m-138)**: about 59% of the run now threads, and the largest single-threaded leg is MFEM's `Reconstruct()` at 18.5% |
| `DEVICE-PLAN.md` | an audit of what a device-resident solve would need. **Not a commitment**, and §3.3's refusals stand. §1.2 and §4's stages 1 and 4 have since been walked → M-130 |
| `BORDERED-GLOBALISATION-PLAN.md` | **the ladder is built and closes none of the six cold failures** — M-116, M-117, M-119. Read it for the design and the four families ruled out; do not start its §8 |
| `BORDER-SCALING-PLAN.md` | **BS-0 is run** — M-117. H2 excluded, H1 and H3 split the six, and BS-5 is the stage its own gate selects |
| `COIL-SUBTRACTION-PLAN.md` | the conductors out of the mesh. **CS-0 and CS-0b measured** → **[M-139](MEASUREMENTS.md#m-139)**: matching the conductor model takes DIII-D's benchmark floor 5.785e-03 → 7.7e-04, **7.5×**, onto §0a's predicted quadrature residual — and **it still does not refine**, so the benchmark cannot resolve MEQ below about 7e-04 whatever either code does. CS-1 onward open, and CS-1 is where filaments arrive |
| `INVERSE-AND-RECONSTRUCTION-PLAN.md` | Part One an inverse solve, Part Two exploratory research on reconstruction. **Nothing run**; its §0 gate is a `freegs4e`-side experiment that costs no MEQ code |
| `FREE-BOUNDARY-PLAN.md` | item 8 — **the staged ladder is CLOSED**: FB-A to FB-7 done or answered, FB-6 met, and §10's XP-0 to XP-4 all met (M-82, M-86, M-87). The conductor model went to `COIL-SUBTRACTION-PLAN.md`. Its **§12 collected the three claims that were still marked unmeasured** as FB-R, FB-S and FB-T, and **all three are now measured**. **FB-R** → M-132: the initial `ψ_bnd` selects among discrete equilibria and **refinement does not merge them** — a factor of 6.3 in `ψ_bnd` at 43,026 elements. **FB-S** → M-133: an X-point on the plasma edge costs the approximation **nothing**, worst rate drop 0.055 of an order over nine `( j, k )` pairs, so the support-corner counter-argument does not bite. **FB-T** → M-134: a prescribed current moves confinement off zero (0/4 → 2/4) and does **not** cure it, where the same lever cures the clamped arm 4/4 — so `theTwoBordersConvergeTogether`'s §11.2 retraction stands and is not going to be lifted this way. Suite 7/7 |
| `PLASMA-EDGE-PLAN.md` | a design out of FB-4, still **not scheduled**. Its own precondition — `j ≥ 1` green, and a machine case at `j ≥ 1` — is now met, so what holds it is the cost-benefit its own numbers make: `j ≥ 1` already gives `k+2` at `k ≤ j` with nothing built |
| `MANTA-COUPLING.md` | the socket MaNTA presents, written from MaNTA's side. No field model is registered there yet |

## So what is next

Nothing is red and stages 0 to 7 are done, so the order is:

1. **Free boundary** — item 8, `FREE-BOUNDARY-PLAN.md`. **Every stage is now
   done or answered**: FB-A, FB-0, FB-1, FB-2, FB-3, FB-5 and FB-7 built and
   measured, FB-4 answered, FB-6 met. §8 below has the per-stage table. **What
   the item is down to** is the conductor model against the reference's
   filaments. **§10's diverted plasmas are done too** — XP-0 to XP-4 all met,
   M-82, M-86 and M-87 — so that is no longer a remainder either, and the
   conductor model has since become its own plan,
   `COIL-SUBTRACTION-PLAN.md`, whose CS-0 measures what this item could only
   infer. The bullets below are the record of how it was reached and are kept
   for the findings in them:

   * ~~**Wire the two borders to TOML.**~~ **DONE 2026-09-06.**
     `[boundary.limiter]` and `[boundary.exterior]` reach
     `setBoundaryFluxPoint()` and `setExteriorCoupling()`, so a machine case is
     no longer a library caller. `examples/free-boundary-halfdisc.toml` is the
     first free-boundary equilibrium MEQ has produced from a file — 7 Newton
     steps, four Gegenbauer coefficients solved for — and
     `theDriverReachesTheExteriorCoupling` pins the driver against the library
     at **7.0e-17**, with a control at **76.6%** saying the coupling is not
     inert. A `bump` initial guess came with it: once the boundary is free the
     guess chooses which equilibrium is reported, so it is part of the problem
     statement.
   * ~~**FB-5's adaptive loop.**~~ **DONE 2026-09-06.** `η` monotone with `Γ`
     fixed (2.51e-01 → 3.21e-02 over four cycles), P.1 preserved on a graded
     `Γ_h`, one Newton step per cycle. **And it found that `η` cannot see the
     coupling**: the exterior coefficients are frozen at 1.3194e-03 across the
     whole loop because the elements touching `Γ_h` carry 0.00% of `η²`. `η` is
     right — it estimates the interior error and the coefficients are a boundary
     functional — but the loop will stall once the interior error passes the
     frozen one, which at cycle 3 is a factor of 1.4 away. ~~**A boundary
     indicator is the cure and is not built.**~~ **BUILT 2026-09-06**, as `η₆`,
     and the plan's own prescription was wrong: *added to the marking* it does
     nothing, because `η₆` is 8.58e-04 against an `η` of 2.51e-01 and never wins
     a Dörfler competition. Two different quantities in different units cannot
     share one threshold however it is weighted, so the boundary term marks on
     its own distribution and the sets are unioned — `Γ_h` then refines
     34 → 46 → 57 → 64 and `|a − exact|` falls **1.32e-03 → 1.53e-04** against a
     control that does not move. Plan §7.12 and §7.12a.
   * **A first coupled free-boundary solve.** Attempted twice on 2026-09-06;
     §7.12 and §7.13 are the record. **The analytic border column is done and
     wired** — `dR/ds` is assembled from the source's own `dF/ds` rather than
     differenced, agreeing with the difference to ten digits and reaching
     **5.70e-16 against 1.69e-13** — and the `( N + 1 )` bordered system now
     closes on a genuinely non-linear source in **4 Newton steps**. **Coils were
     the wrong instinct** and the write-up says why: the flux already crosses
     zero at `r ≈ 1.03`, so an edge exists without a conductor, and both failed
     attempts had put the limiter out in the tail. **What is left is one idea**:
     a plasma-current constraint, which is how CEDRES++ and FreeGS pose it, and
     which turns an ill-conditioned non-linear eigenvalue balance into an
     ordinary unknown. **DONE 2026-09-06**: `setPlasmaCurrent()`, all four
     Jacobian pieces analytic, and it converges the moving-support solve that
     previously failed from everywhere — **63 Newton steps**, current delivered
     to 3e-08. **And what looked like branch selection was the INPUTS**, which
     §7.18 established by re-measuring §7.14 and §7.15 against the repaired
     code and finding half their verdicts false: all four coil currents
     converge where the record said three failed, a **vertical field does not
     make a core** (the O-point moves 1.377 → 1.298 over `0 → −0.20`), a vacuum
     start is not refused but converges in 2 steps reporting unknowns that mean
     nothing, and continuation in the coil current *works and is not needed*,
     reaching the failing currents on a degenerate branch. What settled it is
     §7.15's own conclusion reached by better evidence: **choose consistent
     inputs** — Shafranov's `B_v` says what vertical field a given `I_p` needs,
     so the conductors and the guess can be made to agree before the solve
     rather than found to disagree during it. §7.12b's fixture was repaired that
     way and is green at five limiter radii.
   * ~~**A machine case.**~~ **DONE 2026-09-06.**
     `examples/limited-tokamak.toml` — coils, a limiter, an exterior coupling, a
     prescribed current, a confined source and a gmsh mesh reaching `r = 0`, all
     live at once — solves in 11 Newton steps and reproduces `freegs4e`'s limited
     circular tokamak: `ψ_ax` **1.3e-04** apart, `ψ_bnd` 5.8e-05, `ψ` at
     **5.3e-03** relative `L2` over the reference's whole box. Free boundary by a
     DtN map against free boundary by von Hagenow Green's functions, HDG Newton
     against finite-difference Picard, C++ against Python. §7.16.
   * **`ψ_ax` was the weak point of all of it, and §11 is now worked through.**
     `ψ_ax` was *the largest nodal value of `ψ_h`*, and nothing in that says it
     is a magnetic axis — three configurations were found reporting one that was
     not. It is now constrained at the **located axis**, a zero of `q_h`
     (§11.5 option 3), which the envelope theorem makes free in the Jacobian; and
     two guards refuse an unphysical answer, one on the axis current density and
     one on the plasma containing the symmetry axis. §11.
   * **FB-7, conductors outside `Γ`** — plan §7.19, **DONE 2026-09-07**.
     `meq::ExteriorCoilSet` and `setExteriorConductors()`: the interior equation
     is untouched and the conductor is a known additive term on **both** halves
     of the coupling, the datum and the transmission row, with no new unknowns.
     Rates 1.980/3.409/4.069 with the datum given against a control 148,165×
     larger, and **one Newton step** on the coupled path. It buys domain
     reduction rather than accuracy — MEQ's coils already carry finite extent
     and current density, and near the plasma that is what matters. **The type
     has no `f()`** and that is the design: an exterior conductor has no
     interior current density, and neither has a filament, which is what lets
     `meq::CurrentFilament` share a set with a rectangle and lets MEQ model what
     `../freegs4e` models. **Acceptance 4 landed later the same day** — the
     interior and exterior routes agreeing where both are legal, which needed a
     second geometry, a conductor being either inside `Ω` or outside `Γ` and
     never both. It is a **fitted rectangle** containing the whole half-disc,
     where the conductor is a source term assembled by quadrature: the two
     routes agree at **1.951 / 3.727 / 3.879** against a datum-removed control
     6306× larger, and the interior route is the *better* of the two there, so
     the difference tracks the extension path's own `O( h )` geometry cost.
   * ~~**FB-6, the machine case**, against `../freegs4e`.~~ **MET 2026-09-07**,
     and plan §7.20 is the record. MEQ solves `H_limited_circular` as a free
     boundary and reproduces the **converged** 513² reference at `k = 3` on
     26375 elements: `ψ_ax` **1.5e-04**, `ψ_bnd` **4.5e-05**, and the profile
     scale **3.8e-04** against a value predicted from the reference's own
     current-scaling factor — all three inside the **7.92e-04** that reference
     sits from its own Richardson limit, with MEQ converged in its own mesh to
     3.2e-05 at order **2.93**.

     **What was in the way was not the solver.** `ψ_bnd` was pinned at the
     potential dof NEAREST the limiter contact, which its own header called *"a
     definition rather than an approximation"* — and is an `O( h )` error at
     every degree, not the `O( h^{k+1} )` nodal choice that phrase implies. It
     made whole solves bit-identical over a **0.025 m** plateau in the requested
     point. `LimiterConstraint::ExactPoint` is the repair and was worth a factor
     of twenty on every column.

     **ONE OF THE TWO THINGS LEFT OF FB-6 LANDED THE SAME DAY.** The limiter
     as a **curve** is built, in the restricted form where it is a polygon the
     mesh is fitted to: `setLimiterSurface()` and
     `LimiterConstraint::LocatedContact`, so `ψ_bnd = max ψ_h` over a union of
     mesh faces and the contact is an **output** of the solve rather than an
     input to it. The row is again just the shape functions — the envelope
     theorem covers both an edge interior and a vertex — and the only thing
     that can check that is the observed Newton **order**, which reads
     **2.000** at `n = 24` and `n = 48`.
     `tests/convergence/LimiterCurve.cpp` is the acceptance and
     `[boundary.limiter] SurfaceAttribute` with
     `theDriverFindsTheLimiterContact` is the route from a file. On the machine
     case finding the contact is **11.7× closer** to the converged reference and
     costs branch sensitivity; the shipped example keeps its prescribed point,
     because reproducing `freegs4e`'s own 129² artefact is what makes that
     comparison well posed.

     ~~**The conductor model is what is left**~~ — **DONE 2026-09-15**, and
     both halves of that sentence were tested rather than one assumed.
     → **[M-95](MEASUREMENTS.md#m-95)**, **[M-96](MEASUREMENTS.md#m-96)**.

     **No `Γ` works, and now that is measured on all seven diverted machines
     rather than argued on one.** The margin — nearest conductor minus furthest
     LCFS point, over the whole axis — is negative everywhere, −0.26 m at best
     and −1.50 m on DIII-D, and sliding the centre away makes it worse. The
     mixed case, some conductors inside and some out, is legal geometrically
     and dies on the DtN's mode convergence instead: at a usable clearance the
     conductors CASCADE, nothing can be externalised on any of the seven, and
     the `Γ` radii already configured turn out to be at the optimum.

     **So the reference-side job is the one that was done.**
     `tools/freegs4e-benchmark/shaped.py` rebuilds any of the seven machines
     with every filament and solenoid a `freegs4e.shaped_coil.ShapedCoil` on
     `conductors.py`'s own rectangle, preserving circuit topology and the
     `control` flags the inverse solve needs; `fgsref.py --shaped` is the
     switch and all seven converge. On DIII-D it is worth **5.0× in relative
     `L2` and 7.2× in `L∞`** — and about 2× on the nodes OUTSIDE the
     conductors, which M-87 did not predict, the finite-size correction not
     stopping at the conductor's edge.

     **The reference side is finished for all seven; the COMPARISON is finished
     for one.** Every machine was run both ways and five of the seven behave
     identically — B and E fail at residuals agreeing to three figures either
     way, A trips the axis guard on both — so those are M-89's cold-start
     finding restated and not the conversion. The two that differ do so in
     opposite directions, TCV gaining a solve and MAST-U losing one, and both
     carry a `Solenoid`, which is where the conversion stops being a
     perturbation and the inverse solve moves root. **What is left here is
     MEQ's cold start**, not the conductor model.

     §7.11's ITER profile family stays the next *harder* case, and its advice
     stands: **open it at `γ = 2` or 3, not the published `γ = 1.395`** — the
     profile's vanishing order at the edge *is* `γ`, so 1.395 caps `ψ*` at 3.895
     and makes `∂F/∂ψ` unbounded there.

   **The fixed-boundary rehearsal is already running and is at 8.8e-06**, which
   is now limited by MEQ's own discretisation rather than by the boundary fit —
   see *MEQ against freegs4e* in `CLAUDE_FB.md`. FB-6 removes the fit from the
   comparison entirely, which is the only way past that floor.

   **AND DIVERTED PLASMAS ARE NOW PLANNED, AS §10 OF THAT FILE.** Every real
   machine MEQ would be pointed at is diverted and all seven `freegs4e`
   benchmark configurations are, so it is on the path rather than beyond it. Its
   gate — a limiter free-boundary solve green — **has opened**, and it is still
   **not scheduled**.

   **XP-1 IS ALREADY DONE, AHEAD OF XP-0, BECAUSE THE DEFECT TURNED OUT TO BE
   LIVE ON A LIMITER CASE.** §10.3's pointwise support test was recorded as
   *latent, because nothing diverted exists to run it*, and that was
   generalising from the dramatic case to the only case: `{ψ > ψ_bnd}` is
   disconnected in **705 of 1333 elements** on `theTwoBordersConvergeTogether`'s
   own configuration, no X-point needed. `meq::PlasmaComponent` and
   `[source] PlasmaConnectivity` are the fix, and §10.3's own prediction was
   half false — the lobes join through the **band** straddling the separatrix
   rather than at the saddle, so blocking the X-point element is both
   unnecessary and resolution-dependent, and what ships is a watershed needing
   no X-point finder and no parameter.

   ~~**XP-0 is the one still open**~~ — **DONE 2026-09-08, and it needed no
   free boundary and no seed.** `Soloviev::nstx()`'s X-point is fixed in closed
   form, and `CriticalPointFinder::sweep()` locates it at **1.814 / 3.566 /
   4.223** over five dyadic meshes, against `k+1` less the axis study's own 0.25
   of slack. `findAxis()` and `tryFindAxisFrom()` cannot reach a saddle at all —
   both filter on `AxisSense` — so **XP-3 will want a seeded saddle entry point
   before its three-row border is built**, which is the one API gap XP-0
   exposed. `CLAUDE_FB.md`, *XP-0*.

   ~~**XP-2 is now the next rung**~~ — **DONE 2026-09-12, and `XPointOuter` is
   green.** `psi_bnd` taken from the located X-point as an outer fixed point
   contracts quadratically — 1.909e-03, 3.732e-06, 9.554e-10, 2.255e-13 — finds
   both nulls of the double-null machine, and sits **4.378e-04 m** from
   freegs4e's X-point with `psi_bnd` agreeing to 0.08%. → **[M-82](MEASUREMENTS.md#m-82)**.
   **It needed two things and neither was sufficient alone**: XP-1's fill
   actually running, which a fixture defect had silently disabled on the one
   diverted case in the tree, and the plasma edge held fixed within each Newton.

   ~~**XP-3, the three-row border, is now the next rung**~~ — **DONE 2026-09-13,
   and `XPointBorder` is green.** `setXPointBoundary()` makes `( r_X, z_X )` two
   more unknowns of the same Newton, so `q_r = q_z = 0` are two more rows and
   the root find IS part of the solve. It reaches XP-2's answer, in one process
   on one fixture, to **1.934e-14 m** in the X-point and 1.3e-13 of the span in
   `ψ_ax`, and its own point sits **6.4e-15 m** from an independent root find on
   the field it solved. → **[M-86](MEASUREMENTS.md#m-86)**.

   **XP-0's API gap was its stated prerequisite and was already closed** —
   `tryFindCriticalPointFrom( …, AxisSense::Saddle, … )` landed on 2026-09-08 —
   **and XP-3 turned out not to need it**: the two rows are the root find, so
   the border needs a point LOCATION per iterate and not a search. The seeded
   search is what XP-2 runs on and what XP-3's acceptance checks against.

   **Two things §10.4 predicted came out the other way.** The corner block
   `∂( q_r, q_z )/∂( r_X, z_X )` is **exact**, not an order down: `∇q` is an
   order down as an approximation of the CONTINUOUS Hessian, and what Newton
   needs is the derivative of the DISCRETE residual, where `q_h` is a polynomial
   on its element. And the two unknowns cost **no backsolve at all**, the field
   residual not containing them. What XP-3 does NOT buy is Newton steps — 21
   against XP-2's 18 on the same fixture — and what it does buy is one outer
   discrete state where XP-2 has two.

   **AND IT FOUND THAT THE FIRST STEP OF EVERY BORDERED SOLVE IS TAKEN WITH THE
   `ψ_ax` BORDER DECOUPLED**, the flag having last been written on the cold
   state that sets the convergence target. Repairing it is three lines, it takes
   XP-3 from 14 Newton steps to 10 with the same answer in every digit, **and it
   moves XP-2 onto a different branch** — `ψ_ax` 8.052e-02 against 8.266e-02.
   So it is branch selection rather than a Jacobian fix, and it is left alone
   and recorded rather than taken.

   **XP-4 IS DONE TOO, 2026-09-13**, and it is the rung that made the campaign
   usable rather than only correct: `[boundary.xpoint]` seeds the X-point border
   from a file and `[solver] PlasmaSupportSweeps` runs the support's outer loop
   in the driver, which is what neither XP-2 nor XP-3 built.
   `examples/diverted-tokamak-xpoint.toml` is the worked case and the driver
   reproduces `XPointBorder`'s in-process answer to every printed digit.

   **Against freegs4e's own diverted equilibrium**, with the null handed over in
   neither direction — MEQ solves for it from a seed 7.07e-02 m away, freegs4e
   finds it by a critical-point search — the X-point is **4.4e-04 m** apart,
   `ψ_ax` 6.9e-04 and `ψ_bnd` 7.7e-04. The FIELD agrees at **5.3e-04** relative
   `L2` over 5478 nodes, and the whole of the 2.97e-01 `L∞` is inside the two
   conductors freegs4e models as FILAMENTS and MEQ meshes as rectangles — the
   conductor model, not either solver. → **[M-87](MEASUREMENTS.md#m-87)**.

   **AND THE MESH IS NOW PART OF THE FILE, 2026-09-13.** `[mesh.generate]`
   describes the half-disc `tools/mesh/halfdisc.py` builds, and `meq-run` — a
   script `file(GENERATE)`d beside the binary — makes it and then execs `meq`.
   The conductors' rectangles are **derived from `[[coils]]`**, so a machine
   whose coils were written twice in two conventions (once in a header comment,
   once as blocks) is now written once, and the mesh aligns to the rectangle
   `meq::Coil`'s quadrature uses to the ulp. `meq` refuses such a configuration
   without `--mesh-ready`, so an edited geometry cannot be answered from the
   previous geometry's mesh; `examples/diverted-tokamak-generated.toml` reaches
   the committed mesh's answer to **1.6e-06 m in the X-point**.

   **The obvious next rung is `tools/mesh/` growing a second generator**, since
   every geometry in the tree is a half-disc with rectangles and MAST-U is not:
   a 47-point limiter polygon and a 91-point wall are what `[mesh.generate]`
   would have to name, and `Tool` exists so that arrives as a second value
   rather than as a second table.

   **What is left in §10 is a machine nobody here designed.** `tools/MAST-U*` is
   a freegsnke configuration and is the next case up in complexity; the ladder
   above is what it is to be run against.

   **`PLASMA-EDGE-PLAN.md` is NOT on this list**, deliberately. It is the route
   out of FB-4's `k ≤ j` cap. **Its preamble's precondition is now met** —
   `PlasmaEdgeConvergence` is green, FB-5's bordered solve carries the moving
   support, and the machine case runs at `j = 2` — so what holds it is no longer
   a prerequisite but the cost-benefit its own numbers make: `j ≥ 1` already
   gives `k+2` at `k ≤ j` with nothing built.

   **PE-0 IS RUN, 2026-09-08, AND IT CHANGED WHAT THE REST OF THAT PLAN IS
   WORTH.** Its central premise is now measured at **all three rungs and is
   true** — an uncut plasma subdomain admits `k+1` and `k+2` at `j = 0, 1, 2`
   alike, and `ψ_h` climbs from the cut mesh's cap of 1.5 to **2.00 / 3.04 /
   4.01**. But `ψ*` reaches `k+2` only with `λ` on `Γ_{p,h}`; **through the
   transfer, which is what PE-3 must do, it reads 2.23 / 3.16 / 4.00** — about a
   full order short, and the same shortfall `ExtensionConvergence` records on the
   outer boundary. So PE-0's acceptance as the plan writes it is not reachable
   with the plan's own machinery, the case is **red on that assertion** and stays
   red, and the acceptance is what has to move. The `ψ_h` result alone is a
   strong case at `j = 0`; the `ψ*` case is weaker than §1's table implied.
   Everything after PE-0 still waits. `CLAUDE_FB.md`, *PE-0*.

2. ~~**Finish the inversion**~~ — item 10, and **every stage is done**, with
   the threading blocker closed since: `traceFromAxis()` takes the axis's own
   element as its seed hint, so the library's extraction path makes no
   `Mesh::FindPoints` call at all. It does **not** make `ContourTracer`
   shareable on its own — `fitByAngle()` and `faceJump()` still reach the
   counted last resort on real fixtures — and `Contour::fallbackLocations` could
   never have seen the seed, which is why the scan sat under a test asserting
   that count is zero. `CLAUDE_INVERSION.md` has it. IN-5,
   open surfaces, was the last: `ContourTracer::traceOpen()` traces both ways
   from the seed and joins, and `meq::fitOpenSurface()` is Chebyshev in
   normalised arc length — **5.3e+05 over 4 → 32 modes against a periodic
   control that gets worse**, since forcing `R( −1 ) = R( +1 )` on a curve
   whose ends are 0.8 m apart is inadmissible rather than merely poor.
   `[output]
   FluxSurfaces` writes the `(Ψ, θ)` grid and `meq::GeometryCache` is the
   per-`ψ` cache `MANTA-COUPLING.md` §5's call pattern requires — confirmed at
   `nodes/surfaces`, 3.3× at 40 nodes over 12 surfaces against the naive that
   locates one surface per node. The 34.7× Vandermonde-plus-GEMM figure was NOT
   taken up and is still available: the file is written from the traced nodes
   rather than from a `SurfaceFit`, so there is no Vandermonde in that path.
3. **The fixed-`q(ψ)` solver** — also item 10, and **the round trip closes**:
   `meq::SafetyFactor` inverts `q = V′ g ⟨r^{-2}⟩/4π²` (MFEM-free, CI-gated) and
   `meq::SafetyFactorSolve` closes the loop on KINSOL. From a `g` 40% too large
   everywhere it recovers the closed form to **2.5e-06** in 6 outer iterations
   and 20 inner solves. **A damped Picard provably cannot do it** — the relaxed
   derivative is `1 + ω( G′ − 1 )`, above one for every `ω > 0` when `G′ > 1` —
   and the Newton is affordable only because the profile is *fitted*, so the
   Jacobian is a handful of coefficients rather than `nFieldDOF`.

   **And it reaches the solve from a file**: `[source] SafetyFactorFile` with
   `SafetyFactorDegree` and `ToroidalFieldGuess`, `examples/q-driven.toml` as
   the worked example, and `theDriverSolvesForTheToroidalField` recovering a
   closed-form `g` to **5.0e-06** in 12 outer iterations and 43 equilibria.
   **One solver serves all of them** — only `gg′` changes, so the mesh, the
   spaces, the forms and the symbolic factorisation survive, and each solve is
   warm-started from the last.

   **What is left is the fixture**: it is a rectangle with one closed plasma
   rather than a machine, and nothing has driven `q` on the limited tokamak,
   where the support moves and `ψ_bnd` is an unknown.

**Item 4 is done and is upstream's**, including the two legs MEQ asked for —
it wants a `meq-integration` rebuild to reach MEQ, not more work. Item 6 is
performance and is not urgent; item 5 is a defect in MFEM's local solves that
MEQ works around and has **not filed**; item 7 is a closed investigation kept
for its answer; **item 12.1 has moved into `INTERPOLATORY-HDG-PLAN.md`**. Items
1, 2, 3, 9, 11 and 12.2 are done and are markers.

**And what is not on this list at all is the standing work the plan table
carries**: MEQ's own threading (A–D built, the payoff unmeasured), the device
audit, the two border plans whose campaign came back negative, and the coil
subtraction. Those are ordered inside their own files.

## The state in one paragraph

The solver works and every claim about it is a measured convergence rate. Stages
0 to 6 are done, **stage 7 is finished** — MEQ is a program that solves on a
curved boundary, refines its own mesh, and restarts from a previous answer in one
Newton step — **toroidal flow is finished**, FL-0 to FL-8, and **solution
inversion is finished**, every stage of it including IN-6's output and IN-5's
open surfaces.

**The driver refuses SIX things, in three kinds.** About the **file**:
`[boundary] Type = "exact"`. About the **build**, checked at startup before a
mesh exists: an `AssemblyMode` the file asked for that the build cannot honour,
a `TraceSolver` naming a package this build lacks, and `TraceSolver = "cudss"`
**even where the build has it** — withheld rather than unimplemented, until
MFEM's integrator offload lands. About the **answer**, and these are the first
two that cost a solve, because `ψ_bnd` is an unknown of the bordered Newton:
**a `ψ_ax` that is not the flux at a magnetic axis**, and **a source that does
not vanish on the symmetry axis**, the second ordered first because a bad `ψ_ax`
is its consequence.

## The nonlinear question is settled

It ran through four wrong answers, which is why it is written down rather than
summarised away. MEQ's Newton was thought to fail on the stiff GS-2 sources.

* **Not Newton, and not the boundary condition.** `refs/MFEM-GS-Newton.pdf`
  reports Newton robust where "conventional Picard-based solvers fail to
  converge", and calls its own *fixed*-boundary solver "significantly easier"
  than the free one.
* **Not the hybridization ordering.** This one was believed and acted on.
  `NLOrdering::LineariseThenCondense` removes every element-local nonlinear
  solve (`GetNumLocalNLIterations()` reads 0) and the stiff cases get *worse*.
  **Confirmed twice**: that mode is now deleted, and `NPC` — the canonical
  version of the same idea, and MEQ's default — reproduces the verdict.
* **Mostly it was under-resolution.** §4.2, §4.3 and §4.5 are not stiff at all:
  raw Newton solves all three in 7 to 17 iterations once resolved, and refining
  `h` or raising `k` cures them independently. §4.2 at `k = 3, n = 48` takes 7,
  and `k = 2` converges on the same `n = 16` mesh where `k = 1` takes 42. The
  whole "stiff sources" narrative came from benchmarking at one under-resolved
  point.
* **What is really left is the coarse mesh, and §4.4.** Newton must start inside
  its basin and on a coarse mesh the Dirichlet datum is not — hence
  `Globalisation::PicardThenNewton`, which matters because an adaptive run has to
  solve *before* it has an estimator to refine with. §4.4 is a different animal:
  `∂F/∂ψ` spans `[−580, +566]` against a first Dirichlet eigenvalue of 22.3, so
  the linearised operator has swept past ~26 eigenvalues and the problem is
  multi-valued. Refinement cannot touch that.

Numbers in `CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles* and *Picard, then Newton*.
**Nothing on this path blocks anything else.**

## Division of labour

**The MFEM tree is worked by another agent** and is receive-only: fetch from it,
write requests into its `doc/`, never branch, commit, check out or build there.
MEQ owns `src/`, `tests/` and `apps/`.

**UPSTREAM NO LONGER TRACKS INTER-PROJECT CORRESPONDENCE UNDER `doc/`, SO MEQ
CANNOT LEARN WHAT IT HAS FILED BY LOOKING FOR ITS OWN DOCUMENTS THERE.**
`46ac3d3bc7` deletes every `*-FROM-MEQ.md` in one commit, on the stated
principle that working notes exchanged with a consumer are not the library's
documentation and that what they establish belongs in doxygen on the thing it is
about, or beside the code it constrains. **That is a policy and not a verdict**:
a document's absence says nothing whatever about whether what it asked for was
built. Reading an empty `doc/` as "everything closed" is the obvious inference
and it is the wrong one. What survives is upstream's own plan documents —
`doc/HDG-BEM-COUPLING.md` is where MEQ's coupling request went, absorbed rather
than answered — and the doxygen on each thing that landed.

**SO THE ONLY TEST OF "LANDED" MEQ CAN APPLY IS THE CODE MEQ BUILDS AGAINST**,
which it always was, and it is now the only one available: look for the symbol,
not for the document.

**AND DO NOT CITE A SECTION NUMBER OF ANY OF THOSE FILES.** Upstream records
that section numbers do not agree across the branch family, and this file has
cited one capability as both `§2` and `§3` of the same request. Name the
capability.

| capability | state |
|---|---|
| **Auxiliary globally-coupled unknowns** — `SetNumAuxiliaryUnknowns()` and the two per-element assemble hooks | **the one unbuilt piece, and it is no longer worth asking for.** Upstream carries it as the single open item of `doc/HDG-BEM-COUPLING.md` and confirms nothing of it exists on any branch. **What retires it is the multi-RHS work**: the entry points MEQ named are already public as `NPCReduce()` and `NPCRecover()`, and `DarcyNPCSolver::ArrayMult` now applies them to several right-hand sides in one pass, so a differenced border is `K` applications of a routine that blocks them. MEQ gets the saving without the new API, **and now takes it**: the bordered step queues every column and flushes once, 1.34× on the DIII-D solve leg at 14 columns — → **[M-98](MEASUREMENTS.md#m-98)** |

**The BEM coupling itself is not wanted and upstream says so in its own file.**
It is marked FULLY OPTIONAL, on Gatica & Hsiao (1995): a circular or spherical
artificial boundary makes the exterior operator exact and diagonal and deletes
the layer potentials outright. `src/meq/ExteriorDtN.hpp` is that route already
built — one number per Gegenbauer mode, no singular quadrature — so MEQ is
listed there as the consumer who took it. **Build nothing against the coupling
plan**, which is upstream's own instruction in it.

**What is NOT filed is the local-solve seed** — item 5. It was found from MEQ's
side and nothing has been written into that tree about it, beyond one paragraph
inside the coupling request using it as the argument for why a caller should not
have to difference a condensed residual. That paragraph survives in upstream's
plan, so the finding is on the record there even though the document is not.

---

## 0. FOUR MFEM branches, one integration branch — standing maintenance

**It is four now, not three, and this is the entry to read before rebuilding
`../mfem/install`.** MEQ needs work from

| branch | what MEQ needs from it |
|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — the curved `Γ`, and `TransferredDatumCoefficient` |
| `direct-solver-symbolic-reuse` | `SetReuseSymbolic()` on the direct solvers |
| `gf-hdg-linearise-first` | **`DarcyNPCOperator` / `DarcyNPCSolver`** — MEQ's ordering. `NLOrdering::LineariseThenCondense` was on this branch and is **deleted** |
| **`gf-hdg-dev`** | **the reconstruction fix — *"The postprocessing closes on the element average, always"*** |

**`meq-integration`** is their merge in `../mfem/mfem-src` and is what
`../mfem/install` is built from. Local only, never pushed, re-created whenever any
of the four moves — so nothing may be committed directly to it.

**HOW MANY MERGES IT TAKES HAS CHANGED THREE TIMES, SO DO NOT TRUST THIS
PARAGRAPH — RUN THE LOOP.** It has read three merges, then two (2026-08-30, when
`gf-hdg-dev` was an **ancestor** of both `gf-hdg-subdomains-dev` and
`gf-hdg-linearise-first`, so it arrived with the base), and as of **2026-09-05
`gf-hdg-dev` is ahead again and needs a merge of its own**, by a single commit
that conflicts in `tests/unit/CMakeLists.txt` alone. Each of those was true when
written. **The topology is the other tree's to change**, so the containment loop
is the authority and this table is a hint. `CLAUDE.md` has the recipe and the
conflict resolutions.

What has not changed is *why* `gf-hdg-dev` is on the list. Dropping it silently
loses the reconstruction fix, `ψ*` goes back to being a different function
wherever `∂F/∂ψ` vanishes, and the adaptive loop goes back a full order — with
no error and a green suite until
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` catches it. It is
free today rather than unnecessary.

**Run the containment loop; it has already fired once.** `CLAUDE.md` carries it
with the recipe. Against `meq-integration` at `fa65a2f932` on 2026-09-01 three of
the four report contained and `gf-hdg-linearise-first` does not, because that
branch has advanced 17 commits since the merge. **Read that as "the branch has
moved", not "the merge failed"** — what MEQ needs from it is in the tree. It
says a re-merge is available.

A standing cost of the arrangement, and it falls on MEQ's side.
## 1. Finish the driver — **done**

`docs/running.rst`, `docs/output.rst` and `docs/configuration.rst` are the
user-facing account; `CLAUDE.md` has the rates and the findings. `meq config.toml` parses, solves and writes the equilibrium three
times over, with exit codes 0/1/2/3, the curved boundary, the adaptive loop on
both paths, and the reactive non-linear ladder. The two bullets that stood open
here are closed: the normalised source is reachable from a TOML file (by FL-8,
written once for both as this item predicted), and the interpolating warm start
is wired — adaptive cycles after the first now start from the previous one
interpolated onto the refined mesh, worth **4, 2, 2 Newton iterations against
4, 4, 4 cold** on a nonlinear source, with the final L2 unmoved to 2.5e-12.

**The nonlinear path the driver ships is reactive and must stay so.** Nothing may
be inferred from `F` about which solver to run, because nothing can be: the ratio
`max|∂F/∂ψ|/λ₁` is computable black-box but the pedestal converges at 7 where the
hole fails at 26, which is two points and not a threshold. A second candidate
detector — the first Newton step making the residual worse — has since been
measured and is *anti*-correlated. **Continuation must not go in** either, for
the stronger reason that it has no black-box form at all; see item 7.

**The driver refuses SIX things**, and they fall into three kinds. About the
**file**: `[boundary] Type = "exact"`, which needs a closed form `meq::Source`
does not carry. About the **build**: an `AssemblyMode` or a `TraceSolver` the
linked library cannot honour, checked at startup before a mesh exists — plus
`TraceSolver = "cudss"` even where the build has it, withheld until MFEM's
integrator offload lands. And about the **answer**, which are the first two that
cost a solve, because `ψ_bnd` is an unknown of the bordered Newton: **a `ψ_ax`
that is not the flux at a magnetic axis**, and **a source that does not vanish on
the symmetry axis**. The last is ordered first, a bad `ψ_ax` being its
consequence. All exit 1 with an explanation rather than approximating.

**The flux block is seeded now**, which closes the one small thing that stood
here. `prepare()` used to leave it at zero when a guess was set, so under NPC a
warm start was inconsistent in exactly the row that couples `q` to `ψ` and
`‖r₀‖` went *up*. `seedFluxFromGuess()` solves the flux row of (8a) itself,
element by element with the weight `r` — which removes the `1/r` rather than
guarding it — so the seeded state **satisfies** the row. **It buys an honest
`‖r₀‖` and not one iteration**, and that is a finding rather than a
disappointment: the whole residual is affine in `q`, so an undamped Newton step
lands on the same point whatever `q` it started from — measured, bit-identical
residual histories from iteration 1 onward. See `CLAUDE_HDGGS.md`, *The flux is
seeded*.

## 2. ~~Symbolic factorisation reuse~~ — **done**

`SetReuseSymbolic()` is on for the Newton and Picard paths and off deliberately
on the linear one, which factorises once.
`theSymbolicAnalysisIsReusedAcrossNewtonSteps` asserts one analysis against one
factorisation per iteration — a count and not a timing, and the only thing that
could notice the reuse lapsing, since a lapse costs speed and nothing else.
## 3. Hygiene — **done**

All six bullets are closed and what they found is in `CLAUDE.md` rather than
here: the pedestal tripwire asserts the two
refinement cures rather than a knife edge that threaded-MKL rounding decides;
`everyNonlinearPathReachesTheSameExactSolution` runs four paths to the same L2;
`postProcess()`'s refusal is retired; `η₅` is rebuilt on
`TransferredDatumCoefficient`, going 4.07e-1 → **9.6e-5** on the coarsest
extension mesh and converging at 2.78 against 0.40; `README.md` is true; and
`docs/` is a Sphinx manual on Read the Docs.

**Two of them are worth one line each, because both changed what the driver
does.** The estimator needs `ψ*` in four of eq. (20)'s five terms, and MFEM's
reconstruction was silently wrong **per element** wherever `∂F/∂ψ` vanishes — a
pure Neumann problem whose mean-value regularisation was being skipped, so a
singular matrix was factored. **A whole-domain norm could not see it**: at an
eighth of the domain dead, individual elements were 20× wrong while the norm read
1.87. The test flipped from asserting the defect to asserting the behaviour and
went green on its own, which is the testing stance working as intended.

**And the README rewrite left behind the transferable part.** It kept a paragraph
insisting that "the environment variable is not optional" **without naming a
variable**, because the `MKL_THREADING_LAYER=GNU` it was about had been deleted
from the ctests and the driver as inert. Removing a variable from the code and
from every test left the sentence *about* it behind, pointing at nothing — the
same species of residue as a stale measurement, and **harder to notice because it
reads as a warning rather than as a claim.**

## 4. Element-local parallelism — MFEM, **DONE upstream, not yet in MEQ's install**

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`, on
`gf-hdg-linearise-first`. Twelve sequential element loops in
`darcyhybridization.cpp` over work independent by construction.

**EVERY ONE OF THEM IS NOW THREADED**, that file says so in its own words, and
the last two to go were the ones MEQ asked for: `ComputeH`, `InvertA`,
`InvertD`, `MultNL`, `ComputeSolution`, `EliminateVDofsInRHS`,
`EliminateTrueDofsInRHS`, `ReduceRHS`, **`NPCReduce` and `NPCRecover`** —
**4.35× at eight threads on the leg, bit for bit**.

**AND THE "UNDER 6%" VERDICT THIS SECTION QUOTED IS WITHDRAWN BY THE PROJECT
THAT MADE IT.** `d7ea90a538` measured the threadable loops at 5.4–6.2% of an NPC
step and concluded Amdahl caps any gain — **on four fixed-boundary,
single-right-hand-side cases**. `NPCReduce` and `NPCRecover` are
`O( elements × columns )` where the integrator-bound legs are `O( elements )`,
so a **bordered** Newton inverts the ratio, and MEQ measured the same two legs
at **30.6% of its step and 1.00× across eight threads** →
**[M-126](MEASUREMENTS.md#m-126)**. *A share measured on a problem with one
right-hand side is not a share.*

**IT IS NOT IN THE LIBRARY MEQ LINKS.** `GetNPCTraversalTime` is in
`../mfem-hdg-dev/fem/darcy/darcyhybridization.hpp` and **not** in
`../mfem/install/include/`, so this buys MEQ nothing until `meq-integration` is
re-created and reinstalled. Per the standing rule, the only test of "landed" is
the code MEQ builds against.

**A second correction from the same file, and it is a live trap.** `InvertA` and
`InvertD` are **not** OpenMP loops; what is threaded beside them is the
*batched* route, which needs `LocalFactorMode::Batched` **and** a device
backend. On the host `mfem::forall` is a serial loop, so asking for the batched
local factorisation takes `ComputeH` from 0.418 s to **1.075 s and it stops
scaling**. The default is `Serial`, so nothing is hit by default — and M-101
already measured the trade flat from MEQ's side, for a different reason.

**And §1 of it has LANDED and MEQ is not using it**: `LocalFactorMode`,
`SetLocalFactorMode` and `CanBatchLocalFactor` are in the *installed*
`darcyhybridization.hpp`, and `LocalFactorMode` appears nowhere in `src/`,
`tests/` or `apps/`. See item 6 for what it does and does not buy.

**AND IT IS UNBLOCKED, MORE CLEANLY THAN WHEN IT WAS FILED.** MEQ's default is
now `NonlinearOrdering::NPC` — `mfem::DarcyNPCOperator` over the full
`(q, ψ, ψ̂)` system — and `GetNumLocalNLIterations()` reads exactly **0** under
it against 3644 for the condensation, asserted by
`SolverContract::theOrderingsAgreeAndOnlyOneIteratesLocally`. So the batched work
really is fixed-size **linear** solves, which is what `linalg/batched/` provides
and what a device wants.

**The `ψ_ax` borders are not an obstacle either.** They were believed to require
condense-first and do not: measured 2026-08-30, both orderings reach the same
`ψ_ax` to every digit printed, and under NPC the border row is exactly `−e_j` and
the corner exactly `1`, neither differenced.

**A LONG HISTORY STOOD HERE AND IS DELETED.** It was the investigation of
`NLOrdering::LineariseThenCondense` — a mode MFEM has since retired as *"a
condensation in disguise"* — including its frozen-correction diagnosis, its 149%
cross-linearisation measurement and the case-by-case parity table. It describes
nothing MEQ runs and is recoverable from git. Three things in it are still live
and are kept in `CLAUDE.md` rather than here: the **falsified prediction** that
an ordering would fix a stiff source, which NPC has now falsified a second time;
the **cross-linearisation technique** for detecting hidden state in a reduced
operator; and the **live defect in the ordering MEQ still ships as the backup**,
where the assembled gradient disagrees with a central difference by 100× once the
element-local solves hit their iteration cap. See *The NPC port* and *Why MEQ's
Newton struggles*.

**Condense-first is retained only so the two can be measured against each
other**, which is what turns "this fails" into "this fails and the other does
not" — and on stiff under-resolved meshes it is still the one that works.

## 5. The element-local solve's frozen seed — MFEM, **found and not filed**

`DarcyHybridization::EliminateVDofsInRHS` copies the flux and potential blocks
into `darcy_u` and `darcy_p` at `FormLinearSystem()` time, and every
element-local Newton afterwards — in every residual evaluation, every
`GetGradient()` and every `ComputeSolution()` — starts from that same vector.
`ComputeSolution()` does not take the output vector as a guess, so there is no
way to supply a better one from outside.

**Two costs, and the second is the one that matters.** On an ordinary Newton path
it is iterations: measured on a high-β profile at `n = 16, k = 2`, 40,000 to
60,000 element-local iterations per outer step, most of them hitting the cap of
100. On anything that obtains a derivative by **differencing the residual** it is
correctness, because a local solve that ran out of iterations returns whatever it
had reached, which is not a function of anything: differencing `ψ_ax` by 9e-6
moved `max ψ_h` from 0.896 to 2.04 and then to 3.84.

MEQ works around it in `GradShafranovSolver::formSystem()`, re-forming the
reduced system from the recovered state once per Newton step, which puts every
local solve within an iteration or two of its answer. **That workaround is
cheap and MEQ is not blocked**, which is why this is item 5 and not item 0.

What a request would ask for is small and specific: a way to refresh the local
initial guess without re-forming the system — the obvious shape is
`ComputeSolution()` reading the vector it is handed. **It has not been written**,
and by this project's own rule it should be written only once somebody is sure
the behaviour is settled rather than mid-change in that tree.

## 6. Threaded MKL costs 140x — and threaded assembly put PARDISO's scaling back in reach — MEQ

**The link-line straddle this item was about is FIXED** (2026-09-01: MEQ builds
its own SuiteSparse, so exactly one MKL is loaded), **and fixing it exposed a
140x regression in the suite.** The stray sequential layer had been making all
of MKL sequential, which is why every timing this project recorded before
2026-08-30 was fast.

What threading MKL actually costs, with the columns separated:

* **Not UMFPACK.** `UMFPackSolver::SetOperator` degrades about 40% across the
  whole thread range, never more.
* **`ComputeH()`'s element-local dense work**, through LAPACK, on blocks of
  order 10–30. `k = 2` untouched; **`k = 3` forty-fold worse** at two threads.
  Not the dense **LU**, which is the detail this item first got wrong —
  `dgetrf` does not move at these block sizes at all — but the
  back-substitutions and the Schur-complement `dgemm`.
* So `MKL_NUM_THREADS=1` is on every ctest.

**~~And that pins PARDISO's own scaling out of reach.~~ — IT DID, AND
`AssemblyMode::Threaded` UNPINNED IT, WITH NO NEW CODE.** The problem was that
`MKL_NUM_THREADS` is process-wide, so buying PARDISO's threading meant paying 40x
on assembly, and the two routes out were `mkl_set_num_threads_local()` around the
trace solve or the element-local factorisation off threaded MKL. **Neither was
needed.** MKL suppresses its own threading inside an active OpenMP region, so
once the element loop *is* such a region the local dense work is **nested** and
`MKL_NUM_THREADS` costs it nothing, while the trace solve — on the master thread,
outside any parallel region — still takes all of them. Measured on a whole
nonlinear solve at `k = 3, n = 16`, `OMP=8`: **1285x between the two assembly
modes at `MKL=8`**.

**So the recipe is `Threaded` + PARDISO with `OMP_NUM_THREADS` and
`MKL_NUM_THREADS` at the SAME value above one**, and threaded assembly is MEQ's
default now. Two configurations to avoid, both met while measuring: **UMFPACK can
never take MKL threads**, its BLAS calls being outside any parallel region and so
un-nested; and **threaded assembly at `OMP_NUM_THREADS=1` with MKL threads on**
is catastrophic — 12.4 s against 0.069 s serial at `k = 3` — because a team of
one thread gets no nesting suppression. `apps/meq.cpp` warns about the second at
startup. `CLAUDE_HDGGS.md`'s *What to do, in order of value* no longer carries this as
its item 0. **`LocalFactorMode::Batched` is explicitly not the answer**: it
batches `InvertA()`/`InvertD()`, which run once from `Finalize()`, while the
per-linearisation factorisation is inside `ComputeElementH()`.

**Two things settled and no longer worth an item.** PARDISO's `n ≈ 3000`
ceiling was Debian's `intel-mkl` 2020.4.304 and nothing about the method; and
the trace solver is a run-time choice, `setTraceSolver()` picking among
`UMFPack`, `Pardiso` and `CuDSS`, with UMFPack still the default because
oneMKL's licence is not everybody's to accept. `TODO`'s PARDISO entry carries
the one question that is left, which is reproducibility under threading.

`CLAUDE_HDGGS.md`'s *Threading, measured* has the tables.
## 7. §4.4, the current hole — a known answer that must stay out of the driver

Solvable: adaptive continuation in the added term's amplitude walks `c₃` from 0
to −18 at `k = 2, n = 32` in **9 solves with 2 retreats**, final solve 8 Newton
iterations, no limit point on the branch. Ten uniform steps stall at −10.8, so
the step control is the whole of it.

**It ramps a constructor argument of the _test fixture_.** `Source` exposes
`f( r, z, ψ )` and `dFdPsi( r, z, ψ )` and nothing else, so a black-box `F` has
no amplitude to continue in and this route is unavailable to a driver on
principle. The expressible analogue `F_λ = λF` is a **different** homotopy —
degenerate at `λ = 0`, where this path starts from the converged pedestal — and
is untested.

So this stays a **test-only investigation**, and its value is that §4.4 is now
known to have a solution MEQ can reach. Making it general needs a
self-parameterising `Source`: an interface change to argue on its own merits, not
to smuggle in.

---

## 8. Free boundary — MEQ, and it solves a machine case

`FREE-BOUNDARY-PLAN.md` is the design and is the one real plan left in the tree.
Every stage is built or answered, each against a measured number:

| | |
|---|---|
| **FB-A** | the axis. `ψ` at `k+1` on a mesh reaching `r = 0`, `q` short by half an order, conditioning `O(1/h)` |
| **FB-0** | `meq::ExteriorDtN`. Against a current loop **1.4e−15**; against CEDRES++'s own boundary form **diagonal to 1.15e-10**, which is the test §3.4 named as the one that would falsify all of §3 |
| **FB-1** | the whole coupling on a vacuum problem. `ψ` at **1.99 / 2.99 / 3.99**, exterior coefficients recovered at **3.30** |
| **FB-2** | a prescribed plasma current, and Ampère's law through the solve at round-off over `Γ_h` |
| **FB-3** | `ψ_bnd` as a second border. `setBoundaryFluxPoint()`, constraint at **1.88e-16** |
| **FB-4** | **answered rather than done**: the order is capped by the PROFILE, not the quadrature — `ψ*` keeps `k+2` exactly when `k ≤ j`. No cut rule was built and that is a decision. `PLASMA-EDGE-PLAN.md` is the route out of the cap and is not to be started yet |
| **FB-5** | the `( N + 2 )` bordered Newton, agreeing with superposition to **2.3e-15** in **one** step, **and adaptivity through it**: `η` 2.51e-01 → 3.21e-02 over four cycles with `Γ` fixed, 0 widened fans. It found that `η` **cannot see the coupling** — the exterior coefficients froze at 1.3194e-03 while `Γ_h` kept all 34 faces — which `η₆` fixes, and only by marking on its own distribution: summed into `η` its share of `η²` is 1e-5 and a Dörfler competition never reaches it |
| **FB-6** | **MET 2026-09-07**, plan §7.20. `k = 3`, 26375 elements: `ψ_ax` 1.5e-04, `ψ_bnd` 4.5e-05 and the profile scale 3.8e-04 against the converged 513² reference, all inside its own 7.92e-04 grid error, with MEQ converged in its own mesh at order 2.93. `examples/limited-tokamak.toml` and `theDriverSolvesALimitedTokamak` |
| **FB-7** | conductors **outside `Γ`**, through the coupling rather than the mesh. `meq::ExteriorCoilSet` and `setExteriorConductors()`: the interior equation is untouched, the conductor is a known additive term on both halves of the transmission condition, and there are no new unknowns. **1.980 / 3.409 / 4.069** with the datum given against a control 148,165× larger, **one** Newton step coupled, and the interior route as a control on a second geometry at 1.951 / 3.727 / 3.879. The type has **no `f()`**, which is what lets a rectangle and a filament share a set |

**THE DISTANCE THAT WAS LEFT WAS DRIVER WORK AND IT IS DONE.** FB-3's and FB-5's
borders were library capability with no route from a TOML file;
`[boundary.limiter]` and `[boundary.exterior]` reach `setBoundaryFluxPoint()` and
`setExteriorCoupling()` now, `[[coils]]` reaches `makeCoilSet`, and
`examples/limited-tokamak.toml` drives coils, a limiter, an exterior coupling, a
prescribed current, a confined source and a gmsh mesh reaching `r = 0` all at
once. **Three of the defects found building that case were driver-side**, in keys
describing a box a file mesh never builds, which is why the regression is a
driver test rather than a library one.

The design, which is unchanged:

* **The exterior is exact, not a BEM.** With `Γ` a semicircle centred on the
  axis, the exterior Dirichlet-to-Neumann map for `Δ*` is **diagonal** in the
  Gegenbauer basis `C_n^{−1/2}(cosθ)`, with symbol `(1−n)/ρ_Γ` and mass
  `2/(n(n−1)(2n−1))` in the `dΓ/r` weight the weak form already carries. So
  CEDRES++'s double surface integral over an elliptic-integral kernel collapses
  to `N ≈ 20–40` numbers. Derived and measured — against an independent exact
  current-loop field, to 2.3e−9 — in that plan's §3, and **not taken from any of
  the papers**, which is why its §3.4 names the test that would falsify it.
* **The coupling is stage 5 with the datum unknown.**
  `refs/CouplingAtADistance.pdf`'s own reference [5] is Cockburn & Solano, which
  is what `setExtension()` already drives.
* **The augmented Newton is the bordered one, `N + 2` times over.** `ψ_ax` is one
  border column today; free boundary adds `ψ_bnd` and the `N` coefficients.
  **Under NPC two of the three border quantities are exact rather than
  differenced** — `b = −e_j` and `d = 1`, because `ψ` is an unknown — and the
  coefficients' column is a raw block, since the NPC residual is unreduced and
  has no condensation for a rectangular block to survive.
* **The order of work was FB-A, then FB-0 to FB-6**, and FB-1 — vacuum with
  coils, whose answer is a sum of loop fields known to machine precision — was
  the stage that tested everything structural against an exact answer. The
  received wisdom that free boundary has no analytic solution is true of FB-4
  upwards and false below it, and spending the early exact answers first is what
  made the plan tractable. It is the transferable half of this item.

**The two things the plan named as most likely to hurt have both been measured,
and neither hurt in the way predicted.** The axis costs `q` half an order and
`O(1/h)` in conditioning, which a direct trace solve does not care about at these
sizes. And cut-element quadrature turned out not to be the binding constraint at
all — the profile's vanishing order `j` is. See FB-A and FB-4 above.

**WHAT THIS NEEDS FROM THE OTHER TREE IS NOW NOTHING, FOR FB-0 THROUGH FB-3, AND
THIS FILE SAID OTHERWISE.** It claimed the two rectangular integrators MEQ asked
upstream for were blocking. They are not, and the reason is the NPC
port: the datum's data half is an **essential trace value**, not a weak form, so
its block is `ProjectBdrCoefficient` against the `PathTraceCoefficient` that
already exists; and the transmission block is reachable from `TransferPath::
Endpoint` and `ElementExtension::TransformBack`, both public, in about forty
lines of MEQ. The auxiliary globally-coupled unknowns are an optimisation over
`N + 2` backsolves against one factorisation — the cost MEQ's `ψ_ax` border
already pays, and one `DarcyNPCSolver::ArrayMult` now blocks. Plan §6.4 is the
per-stage table. **Neither request exists as a document any more** and that is a
policy of upstream's rather than a verdict; see *Division of labour*.

**AND FB-4's "one real gap" CLOSED BY NOT OPENING IT.** MFEM *does* have
cut-element quadrature — `mfem::MomentFittingIntRules` in `fem/intrules_cut.hpp`,
`MFEM_USE_LAPACK`-gated and installed — and what it does not have is the
**sensitivity of a cut rule to the level set**, which is the half CEDRES++
actually names. **MEQ adopts no cut rule, so there is nothing to supply**: with a
fixed Gauss rule the quadrature points do not move, so the assembled Jacobian is
the *exact* derivative of the assembled residual whatever the edge is doing.
Adopting a cut rule is what would create the gap it was meant to close.

Both were exercised rather than reasoned about. `MomentFittingIntRules` works —
exact on a straight cut, 3.3e-06 on a disc — but on one mesh in three it produced
weights of **−30 and −63** on elements of area 1e-4 with the area 1.3e-02 wrong
by cancellation, which is the known conditioning fragility of moment fitting on a
nearly degenerate cut. And **both of MFEM's cut backends are quadrilateral and
hexahedral only** where MEQ's meshes are triangles, so turning `MFEM_USE_ALGOIM`
on would cost a mesh-type change as well, for a threshold that does not move.

`NORMALISED-LINEARISE-FIRST.md` was the design for MEQ's half of this under
`LineariseThenCondense`. **Both the mode and the design are deleted** — the
mechanism it describes has no code path left to run on, and under NPC the
problem it solved does not arise. Plan §4.1; git has it if it is ever wanted.

---
## 9. Toroidal flow — MEQ, DONE, FL-0 to FL-8

**`src/meq/RotatingSource.{hpp,cpp}` solves `refs/RotatingGK.pdf` (136), closed
by its (96) and (97), and it is reachable from a TOML file.** Two species in
closed form, `n` species by a safeguarded root find, normalised flux through the
existing bordered Newton, and `[source] Type = "rotating"` with
`examples/rotating-rectangle.toml` and `rotating-normalised.toml` as the worked
examples. `docs/rotation.rst` is the derivation; `CLAUDE_FLOW.md`'s
*Toroidal flow* has every measurement, the three errors found in Li & Zhu, and
the Maschke–Perrin reading **this file previously got wrong** — its §4 is (136)'s
isothermal closure at every `γ`, and the paragraph that stood here called it an
adiabatic one.

**Maschke–Perrin is in `tests/analytic/`** as `MaschkePerrin.hpp`, with
`MaschkePerrinConvergence.cpp` driving it. Be precise about what it is worth:
its **§4 PDE is Li & Zhu's (12) renamed**, so nothing about the discretisation
is new. What is new is that its (4.7) constrains only `ω²/(R̄T)`, so it is the
**one exact solution MEQ has with `T′ ≠ 0` and `ω′ ≠ 0`** — every other
closed-form check of `meq::RotatingSource` runs at constant `T` and constant
`ω`, and the varying case was covered only by a central difference of MEQ's own
`f()`. Measured: the production source reproduces (4.10) to 4.3e-16 with a
`dF/dψ` of 2.6e-15 against `O(1)` cancelling terms, and drives the solve at
1.996 / 2.997 / 3.998.

**The one thing it turned up is now closed**, and by the manufactured route it
predicted: `C` constant is what makes the equation solvable in closed form at
all, so no paper could have supplied it.
`tests/analytic/VaryingCentrifugal.hpp` prescribes `C(ψ)` as a quadratic and
derives `ω` from (97), so `ω²/T` varies across surfaces and `C` drifts 2.50×;
two independent implementations of (96)/(97)/(136) then agree at **7.1e-16 in
`F`** and **5.3e-16 in `dF/dψ`** on both closures, with a bisection route in
which `C` appears nowhere as the third leg. Dropping `C′` moves `F` by 7.9e-01
and flipping its sign — Li & Zhu's (9) — by 1.6e+00.

**It also closed a hole nobody had named**: every two-species rotating
configuration in the tree ran at `Z = ±1`, where the closure's charge-weighted
combinations `Z₁T₂ − Z₂T₁` and `Z₁m₂ − Z₂m₁` collapse to `T₁ + T₂` and
`m₁ + m₂`, so **a closed form written with the plain sums would have passed
everything**. The fixture is at `Z₁ = +2`. `CLAUDE_FLOW.md` has both.

## 10. The fixed-`q(ψ)` solver — MEQ, and the round trip closes

**`INVERSION-PLAN.md` is the design, and every stage is done and green.** This
item became reachable at **IN-2**, where the flux-surface averages
`⟨r^{-2}⟩_ψ` and `V′(ψ)` are measured against a converged reference on the exact
field; **IN-6** writes them to a file against a flux label; and **IN-5** covers
the surfaces that do not close, which needed free boundary to exist before there
was one to trace.

**Take `q(ψ)` as input and find `I(ψ)` from it**, rather than taking `I(ψ)`
directly as items 1 and 9 both do. It is how a transport code hands an
equilibrium code its target, and it is what a coupling to MaNTA will want. RoPP
(142) is the relation:

```
q(ψ) = V′(ψ) I(ψ) ⟨r^{-2}⟩_ψ / 4π²
```

**AND THE SOLVER IS BUILT.** `meq::SafetyFactor` is the inversion —
`g = 4π²q/(V′⟨r^{-2}⟩)`, one division per surface, MFEM-free so CI gates it —
and `meq::SafetyFactorSolve` is the outer Newton that closes it, on KINSOL.
From a `g` 40% too large everywhere it recovers a closed form to **2.5e-06** in
**6 outer iterations and 20 inner solves**.

**IT IS AN OUTER LOOP AND NOT A NON-LOCAL JACOBIAN**, which is the shape this
item did not anticipate. The paragraph below asked for `∂F/∂ψ` to carry the
whole `V′⟨r^{-2}⟩` dependence — a continuum of border rows — and that is not
what was needed: the profile is **fitted**, so the outer unknown is a handful of
coefficients and a differenced Jacobian costs a few solves per step rather than
`nFieldDOF` re-extractions. The fit was put there for conditioning and it paid
for the Newton as well. The inner solve is MEQ's existing bordered Newton,
untouched.

**A damped Picard provably cannot replace it.** The relaxed iteration
`c ← c + ω( G(c) − c )` has derivative `1 + ω( G′ − 1 )`, which is above one
for **every** `ω > 0` when `G′ > 1`: under-relaxation stabilises a map that
oscillates and does nothing at all for one that runs away. Measured, this map
runs away.

**AND IT REACHES THE SOLVE FROM A FILE.** `[source] SafetyFactorFile` makes
`gg′` an output, `examples/q-driven.toml` is the worked example, and the `.nc`
carries the recovered `g^2` coefficients — on this route there is no
`GGPrimeFile` beside the output for a consumer to look the field up in.
**One solver serves every map evaluation**, through
`NormalisedMHDSource::setGGPrime()`.

**WHAT IS LEFT IS TWO THINGS, AND THE SECOND WAS FOUND BY CLOSING
`INVERSION-PLAN.md` RATHER THAN BY USING THE CODE.**

**THE FIXTURE.** It is a rectangle with one closed plasma rather than a machine,
and nothing has yet driven `q` on the limited tokamak, where the support moves
and `ψ_bnd` is an unknown. Nothing refuses the combination —
`readToroidalFieldTarget()` requires only `Normalised = true` — so this is an
experiment and not a port. Note the interaction it would test: `PlasmaCurrent`
scales the whole source through the border unknown `λ`, and `q` is a statement
about `g`, so the outer loop and the current border both act on the toroidal
field and it is not obvious they are consistent.

**AND `q` WAS NEVER REPORTED, SO THE ROUND TRIP COULD NOT BE READ OFF THE
OUTPUT. FIXED ON THE DRIVEN ROUTE** → **[M-131](MEASUREMENTS.md#m-131)**.
`FluxSurface::safetyFactor` was computed and `Output.cpp` had the column, gated
on `FluxFamilyOptions::toroidalField`, which `apps/meq.cpp` set at neither
extraction site. It now sets it from the coefficients the outer Newton solved
for, `examples/q-driven.toml` carries `[output] FluxSurfaces`, and **the
reported `q` is the target to 4.5e-07 over 24 surfaces** — an order better than
the same run's `g`, because the target was measured on this mesh.
`DriverAcceptance::theQDrivenRunReportsTheSafetyFactorItReached`.

**The refusal stands where it was right.** A prescribed-field run still reports
no `q`: a `meq::Source` carries `gg'` and not `g`, the constant of integration
has no key, and a column of zeroes is indistinguishable from a machine with no
toroidal field. **Giving it one is a schema decision and not a wiring job** — a
new `[source]` key carrying the vacuum `R₀B₀`, which if wrong yields a plausible
`q` that nothing contradicts. That is the remaining half.

Three things written down before it was started, kept because two of them held:

* **The non-local Jacobian has a precedent in the tree.** `ψ_ax` is one border
  row and column today and `HighBetaConvergence` is the acceptance criterion for
  it; a `q`-driven `I(ψ)` is the same shape with a continuum of rows rather than
  one, which is the part that is genuinely new.
* **It is independent of rotation.** A fixed-`q` static solver is useful on its
  own and is the smaller problem; doing it first and then composing is likely
  cheaper than doing it inside item 9.
* **Two things from `INVERSION-PLAN.md` changed this item's shape.**
  `v0-legacy`'s `FluxSurfaces` **has never compiled**, so there was an algorithm
  to reuse and no working tool to extend; and the global-structure work this item
  was assumed to need is **deferred**, because a fixed-boundary problem with one
  axis has no interior saddles at all.

## 11. `B` in the band — MEQ, DONE 2026-09-02

**FOUND WHILE WIRING THE ROTATING OUTPUT AND FIXED THE SAME WEEK.** Only
`samplePotentialWithFlux()` applied the Taylor step into the band between `Γ_h`
and `Γ`, so `B` — which went through `sampleComponent()` — was read at the foot
on `Γ_h`: about one node in ten of the interchange file piecewise constant,
`O(h)`, behind a mask saying `inside = 1`.
`GridSampler::sampleComponentWithGradient()` is the fix, at **rate 2.20 against
0.92 for the foot**, and the `.nc` now carries `byte extrapolated( Z, R )` beside
`inside` so a reader can drop continued nodes rather than being told only how
many there were.

**What is deliberately NOT done**, recorded so nobody re-derives it: this is
`O(h²)` at every `k` and does not reach `ψ`'s order, because `ψ` is continued
with a *solved* variable and `∇q` is a *differentiated* one. `div q = −F/r` and
`∂_r q_z − ∂_z q_r = −q_z/r` pin two of `∇q`'s four entries exactly, but leave
the symmetric traceless part differentiated — structure rather than an order, and
it would need the source plumbed into `GridSampler`. `CLAUDE.md`'s *Status*
section has the detail.

## 12. Interpolatory HDG, and two experiments beside it — MEQ, open

Came out of one question — *is there mileage in evaluating the GS right hand
side at `ψ*` rather than `ψ_h`, and could an interpolatory method save the
repeated assembly?* — and the answer turned out to be that **those are not two
questions.** The `ψ*` evaluation is what makes the interpolatory method keep its
superconvergence, and neither half is worth much without the other.

### 12.1 Interpolatory HDG — assemble the source ONCE

**`INTERPOLATORY-HDG-PLAN.md` SUPERSEDES THIS SUBSECTION AND IS WHERE THE WORK
IS.** Two things in it change the shape of the item and both postdate what
follows. **MFEM has already built the method and MEQ links it today** —
`HDGInterpolatoryReactionIntegrator`, `NodalReactionFunction`,
`HDGPostprocessBlocks` and `DarcyHybridization::Bg_data`, all on
`gf-interp-hdg-dev` and all inside `meq-integration` — so this is an MEQ wiring
job against a built, documented, unit-tested facility rather than a feature
request with a consumer attached. **And the headline below is deflated by that
plan's own arithmetic**: at `extraOrder = 0` interpolation evaluates `F` *more*
often than quadrature does, the honest figure is 1.6–2.0× fewer evaluations,
`ComputeH()` is untouched, and the Amdahl ceiling off M-80 is **23–28%**. Read
its §0 before its stages; IH-0 is a falsifier whose result may be *stop*.

What is kept here is the argument for the method, which is unchanged and is why
the plan exists.

**THE METHOD.** Interpolate the nonlinearity into a finite element space instead
of integrating it. The load becomes `A F⃗` with `A` assembled once and `F⃗` the
vector of pointwise `Source::f()` values; the Jacobian becomes a diagonal
scaling of the same fixed matrices. No quadrature of `F` anywhere and nothing
about the source reassembled between Newton steps. It is a **different
discretisation** — a variational crime — and not an implementation of the one
MEQ has, so the acceptance is a rate and never a bit comparison.

**AND THE OBVIOUS VERSION LOSES `k+2`.** Interpolating `F( ψ_h )` into `W_h` is
Cockburn, Singler & Zhang (2019), which proves optimal rates for `ψ_h` and `q_h`
and reports no superconvergence of the postprocessed `ψ*` at all. For MEQ that
would be silent and expensive: `OutputConvergence`, `Estimator` and the whole
adaptive loop rest on `ψ*` being a full order better. **The fix is to evaluate
the nonlinearity at `ψ*`** — `I_h F( ψ*_h )`, interpolating into the
postprocessing space `P^(k+1)` — which is Chen, Cockburn, Singler & Zhang I,
Remark 2.1, and it restores `k+2`: their Table 1 reads 2.95, 3.02, 3.02, 3.01 at
`k = 1`. `refs/Refs.md`, *Interpolatory HDG*, carries both with the dois checked.

Their `HDG_k` is MEQ's method — `V_h`, `W_h`, `M_h` all at degree `k`, `τ`
constant and `O(1)`, simplices — so paper I transfers and the HHO-flavoured
sequel does not without changing the spaces and the stabilisation.

**WHAT IT COSTS MEQ, IN THE ORDER THE WORK WOULD BE MET.**

* **`ψ*` MUST COME FROM THE CLASSIC POSTPROCESSING AND MEQ DOES NOT USE IT.**
  The method needs Nguyen, Peraire & Cockburn eq (25) — a pure Neumann local
  problem in `q_h` closed by the element average of `ψ_h` — because that is
  LINEAR in the unknowns, which is what makes the `B11`, `B12` of §2.2 constant
  and assemblable once. MFEM supplies exactly it as
  `mfem::HDGPotentialPostprocessor`. MEQ instead builds `ψ*` from
  `DarcyForm::Reconstruct()`, the richer mixed reconstruction, which **lifts the
  non-linear potential integrators into its local problem** and takes their
  gradient at the computed potential — so its local matrix carries `∂F/∂ψ` and
  feeding its output back into the source really would make each element's
  reconstruction an implicit local fixed point. The two answer different
  questions and `postprocess_hdg.hpp` says so; the interpolatory method wants
  the smaller one.
* **THE SOURCE ACQUIRES A JACOBIAN BLOCK AGAINST THE FLUX.** `ψ*` is built from
  `q_h` as well as `ψ_h`, so the paper's Jacobian is
  `A9 diag( F'( γ ) ) B11` against the flux coefficients and
  `A9 diag( F'( γ ) ) B12` against the potential's. MEQ's reaction term is a
  potential-potential block today. Everything stays element-local, so
  hybridization survives, but this is not a new integrator dropped into the
  existing nonlinear potential mass — `LocalNLOperator` and the parity gap both
  need re-measuring.
* **`ConfineToPlasma` IS THE REAL RISK.** `F` is discontinuous across
  `{ Ψ = 0 }`, and an interpolant of a discontinuous function on a cut element
  is an oscillating polynomial taking nodal values from the wrong side. The
  answer is a hybrid — interpolate on uncut elements, keep the cut quadrature
  where the edge is — which keeps the exact treatment exactly where
  `PLASMA-EDGE-PLAN.md`'s rate cap lives and takes the fast path in the bulk,
  the cut elements being `O( h^-1 )` against `O( h^-2 )`. The papers assume a
  smooth `F` satisfying a local Lipschitz condition and say nothing about this.
* **The axisymmetric weights and the explicit `( r, z )` dependence** are
  outside the papers' `F( u )`, but only notationally: `A` becomes
  `∫_K ( 1/r ) φ_j w_i` and `F` is evaluated at `( x_j, ψ*( x_j ) )`.

**WHAT IT WOULD BUY, AND THE HONEST SHARE.**
`../mfem-hdg-dev/doc/HDG-DEVICE-OFFLOAD.md` puts the integrators at **46–53%**
of an NPC step, but that is ALL integrators; the `ψ`-dependent part is not
apportioned and measuring it is the first thing to do. There is a second prize
that is not a share of anything: **it should retire
`setSourceQuadratureOrder()`.** The only reason that knob is swept into
Grundmann–Möller territory — negative weights reaching `−1.9e+07`,
`QUADRATURE-HIGH-ORDER-TRIANGLES-FROM-MEQ.md` — is that `F( ψ_h )` is a
nonpolynomial integrand. Under interpolation the only quadrature left is
`∫ ( 1/r ) φ_j w_i`, which the tabulated rules cover.

**ACCEPTANCE.** The ladder, not a tolerance: `k+1` in `ψ` and `q` and `k+2` in
`ψ*` on `ManufacturedNonlinear` and `SimilarityExponential` at `k = 1…3`, plus
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes`. `k = 0` is outside
paper I's theory and is where the sequel's method (B) would be needed.

### 12.2 Locate the plasma edge with `ψ*` — MEASURED AND NOT WORTH BUILDING

**ANSWERED 2026-09-14, and it cost no new code because the oracle was already in
the tree.** → **[M-94](MEASUREMENTS.md#m-94)**. `tests/analytic/PlasmaEdge.hpp`
carries two fixtures over one exact solution: `MovingPlasmaEdge` reads its
support off `ψ_h`, which is what MEQ does, and `PlasmaEdge` cuts at the **exact**
edge. Since `ψ*` is a better approximation to `ψ` than `ψ_h` is, and the fixed
cut *is* `ψ`, the second column is an **upper bound on everything this item could
buy**. `locatingTheEdgeExactlyBuysNoOrder` is the case.

**NO ORDER IS ON THE TABLE.** If the level set limited the rate, the oracle's
advantage would grow as `h` fell; over a fourfold refinement it grows by at most
**1.018** in `ψ_h` and **1.121** in `ψ*`, and in the one row where the oracle is
clearly ahead it *falls*, 2.98 → 2.03.

**AND IN THE REGIME MEQ RUNS IN IT BUYS NOTHING AT ALL.** At `k ≤ j`, where `ψ*`
keeps `k+2`, cutting at the exact edge is worth between **0.95× and 1.28×** — the
moving edge is the better of the two in half those rows. The one corner with a
real factor is `j = 2, k = 3`, about **2×**, which is `k > j`: the regime where
the cut already caps the order, which is what `PLASMA-EDGE-PLAN.md` addresses by
a different route, and where `ψ*` would capture only the part of the gap lying
between `ψ_h` and `ψ`.

**THE PREMISE BELOW IS THE PART THAT WAS WRONG**, and it is kept because it is
the reasoning to not repeat: *"the place `ψ`'s accuracy DOES bind is the level
set"*. It does not bind there either. What limits the cut is the **rule's
blindness** — a quadrature rule cannot see a kink between its points, whatever
field decided where the kink is — which
`theLossIsTheRulesBlindnessAndNotItsResolution` measures separately and which
no improvement to the level set can touch.


Substituting `ψ*` into the source **for its own sake**, with the quadrature
left in place, is not worth it: the consistency error of `F( ψ_h )` is
`O( h^(k+1) )`, the method's own rate, so `ψ*`'s `O( h^(k+2) )` removes a term
that limits nothing and buys a constant. Its value in 12.1 is structural — it
is what the interpolation needs to stay superconvergent — and not that the
source was inaccurate.

**The place `ψ`'s accuracy DOES bind is the level set.** `ConfineToPlasma`
locates `{ Ψ > 0 }` as a level set of `ψ_h`, and `PLASMA-EDGE-PLAN.md`'s rate
cap is about the geometry of that set rather than about the load. Locating the
edge with `ψ*` while still assembling `F` at `ψ_h` is independent of 12.1, needs
no new analysis and no Jacobian change, and has `PlasmaEdgeConvergence` waiting
as its acceptance criterion. It is an afternoon.

**It was an afternoon, and it was spent on the measurement instead of on the
build — which is the cheaper order whenever an oracle exists.** Bounding what a
change can buy before writing it needed two fixtures that were already there and
one test case; writing it would have needed a second support field through
`meq::SourceIntegrator` and a source-side way to evaluate `F` at one `ψ` while
testing the support at another, which is an API change across all three
`NormalisedSource` subclasses.

### 12.3 Reuse the Jacobian across Newton steps — the fallback

A chord or Shamanskii iteration. `m = 3` would skip two of every three
assemblies *and* two of every three `ComputeH()` factorisations, which is the
46–53% plus the 7–10%. KINSOL implements it (`msbset`) and MEQ already links
KINSOL for `AndersonPicard` and `PicardThenNewton`, so the machinery is in the
build. A sweep over `m = 1, 2, 3, 5` against `ManufacturedNonlinear` and
`examples/limited-tokamak.toml`, reporting wall clock **and** iteration count,
answers whether the extra iterations eat the saving. Expect the moving support
to be the trap: it makes the residual mildly discontinuous, so a stale Jacobian
may stall where a fresh one does not.

**IT IS THE FALLBACK BECAUSE 12.1 DOMINATES IT.** Interpolatory HDG buys the
same assembly saving while keeping the Jacobian EXACT, so Newton stays
quadratic; this buys it by giving the convergence rate up. Reach for it if
12.1's rates do not hold.

## Deliberately not yet

* **GPU and cuDSS.** Correctness-testable here; this card cannot say whether it
  is worth it. Consumer FP64 is 1/32 of FP32 where datacentre parts are 1/2, so a
  local timing can invert the production conclusion.

  **AND THERE IS A SECOND REASON, WHICH IS THE BINDING ONE AND IS NOT ABOUT THIS
  MACHINE.** A device trace solve is only worth having if the data STAYS on the
  device. `../mfem-hdg-dev/doc/HDG-DEVICE-OFFLOAD.md` — under construction —
  measures the shares of an NPC step and gates the whole device path on the
  **integrators** (46–53%), because the trace solve (26–31%) and the local dense
  algebra (7–10%) can be moved almost for free and moving only those is *"worse
  than doing nothing … plausibly slower than staying on the host throughout"*.
  MEQ's integrators and scatter have no device kernels.

  So `[solver] TraceSolver = "cudss"` is **refused by the driver** — withheld
  rather than unimplemented — and the library keeps cuDSS only so that
  correctness on the device path is checkable. **The key opens when MFEM's group
  2 lands**, and that is a request on `../mfem-hdg-dev` rather than MEQ work.
  Note the entry below: partial assembly for the HDG integrators is exactly what
  group 2 is, so these two items are one item seen from two sides.
* **Partial assembly.** Not implemented for the HDG integrators at all, and its
  case is sum factorisation, which wants tensor-product elements where MEQ is
  triangles. A discretisation decision, not a flag.
* **A fully device-resident solve.** `DEVICE-PLAN.md` is the audit and it is
  explicitly not a commitment. What has changed since the paragraph above is
  that *"group 2 needs a partial-assembly rewrite and is not built"* is now
  only half true: there is still no partial-assembly route and MEQ is
  simplices, so **that** economic case does not apply here at all — but what
  upstream built instead is batched element and face kernels behind `CanBatch*`
  predicates, MEQ passes most of them, and upstream's own chain measurement
  reads **1.33× at 4096 elements** where its headline used to say *worse than
  doing nothing*. **The refusal of `TraceSolver = "cudss"` from a file stands
  regardless**, on the trade argument alone and on MEQ's own leg profile: M-126
  puts the trace solve at 20.1% of a threaded step, so group-4-alone buys a
  share of a fifth while paying a round trip on the other four fifths. Three
  MEQ-specific facts keep upstream's crossover from being MEQ's — their fixture
  is quads and MEQ is triangles everywhere, five of MEQ's eight shipped meshes
  are *below* 4096 elements, and MEQ's Newton is bordered with the border being
  host arithmetic that would not go to a device.
* **The rest of the physics** — anisotropic pressure, NetCDF profiles, MaNTA
  coupling. `TODO` carries each with what has been established, and anisotropic
  pressure still has **no reference pinned**, which is the first thing it needs.
  **Sonic rotation has left this list**: it is item 9 and it is built, which is
  why `TODO`'s *Sonic toroidal rotation* entry is a stub pointing at
  `docs/rotation.rst`. `TODO`'s PARDISO entry is down to its one open question,
  reproducibility under threading, and its performance entry to what the
  threading campaign did not already answer.

## The standing rule

Every stage ends at a **measured** number, not at "it runs". That is what caught
the Solov'ev coefficients, the `τ` sign, the `ψ*` regression, the hybridization
ordering, the under-resolution and the inert normalisation — six times the code
or the story looked right and was not. Four corollaries this project has paid
for:

* **A property measured on the easy configuration is not a property of the
  code.** Symmetry held to 2e-16 on a fitted rectangle and failed at 5.4e-1 on
  the geometry MEQ exists for.
* **A difficulty measured at one resolution is not a property of the problem.**
  Three sources were called stiff for months on the strength of `k = 1,
  h = 0.05`, and a whole MFEM work item was requested to fix them. They were
  under-resolved.
* **A measurement taken on this machine may not be a measurement about the
  code.** Threaded MKL decides a marginal Newton; a consumer GPU decides nothing
  about a production one. And a *stale* one is worse than none: this file said
  PARDISO stopped working above `n ≈ 3000` for months after the MKL it was a
  verdict on had been replaced.
* **A derivative obtained by differencing is only as good as the thing being
  differenced is a function.** The bordered Newton's border read 1.6e5 where it
  should read 1, and looked exactly like a singular Jacobian, because an
  element-local solve underneath it was hitting an iteration cap and returning
  whatever it had reached. Nothing in the arithmetic was wrong. Before believing
  a finite difference of a hybridized residual, check that the local solves
  under it converged.
