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

The plan files. **`DRIVER-PLAN.md` and `FLOW-PLAN.md` are gone**, converted to
`docs/` on 2026-09-06 under the standing rule that a plan with nothing left in it
is documentation: everything they staged is built, so what a reader needs is the
manual and what a maintainer needs is `CLAUDE.md`. Git has them.

| | |
|---|---|
| ~~`DRIVER-PLAN.md`~~ | stage 7 — **done**, and now `docs/running.rst`, `docs/output.rst` and `docs/configuration.rst` |
| ~~`FLOW-PLAN.md`~~ | item 9, FL-0 to FL-8 — **done**, and now `docs/rotation.rst`, which carries the derivation `RotatingSource.hpp` defers to |
| `INVERSION-PLAN.md` | item 10's machinery — **every stage done**, IN-A to IN-P |
| `FREE-BOUNDARY-PLAN.md` | item 8 — FB-A, FB-0, FB-1, FB-2, FB-3, FB-5, FB-7 **done**, FB-4 **answered**, §10's XP-0 to XP-4 **done**, and **FB-6 met 2026-09-07**: MEQ reproduces `freegs4e`'s converged limited tokamak to 1.5e-04 in `ψ_ax`, with the limiter as a curve landing the same day. What is left of it is the conductor model — MEQ's rectangles against the reference's filaments, a reference-side job — and §10's diverted plasmas, whose gate has opened |
| `PLASMA-EDGE-PLAN.md` | a design out of FB-4, still **not scheduled**. Its own precondition — `j ≥ 1` green, and a machine case at `j ≥ 1` — is now met, so what holds it is the cost-benefit its own numbers make: `j ≥ 1` already gives `k+2` at `k ≤ j` with nothing built |
| `MANTA-COUPLING.md` | the socket MaNTA presents, written from MaNTA's side. No field model is registered there yet |

## So what is next

Nothing is red and stages 0 to 7 are done, so the order is:

1. **Free boundary** — item 8, `FREE-BOUNDARY-PLAN.md`. **Every stage is now
   done or answered**: FB-A, FB-0, FB-1, FB-2, FB-3, FB-5 and FB-7 built and
   measured, FB-4 answered, FB-6 met. §8 below has the per-stage table. **What
   the item is down to** is the conductor model against the reference's
   filaments, which is a reference-side job, and §10's diverted plasmas, which
   is a new campaign rather than a remainder. The bullets below are the record
   of how it was reached and are kept for the findings in them:

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

     **The conductor model is what is left**, MEQ's rectangles against
     `freegs4e`'s filaments, which FB-7's `CurrentFilament` cannot fix on this
     geometry because no `Γ` both encloses the plasma and excludes the inner
     coil pair. Rebuilding the reference on `freegs4e.shaped_coil.ShapedCoil`
     is the way, and it is a reference-side job.

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

Items 4 and 6 are performance and neither is urgent; item 5 is a defect in
MFEM's local solves that MEQ works around and has **not filed**; item 7 is a
closed investigation kept for its answer. Items 1, 2, 3, 9 and 11 are done and
are markers.

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

**Every request MEQ has filed is landed in the code MEQ builds against except
the three still open**, which is the only test of "landed" MEQ can apply. **Ask
`git` which documents exist, not `ls`**: the older open ones live on
`gf-hdg-linearise-first` alone, so a listing taken while that tree sits on
`gf-hdg-dev` shows an almost empty `doc/` and means nothing — and one is
**untracked**, so `git status` there is how to see it rather than `git cat-file`.
`CLAUDE.md` has the branch-by-branch table and the mistake that produced it.

| filed | state |
|---|---|
| **`HDG-BEM-COUPLING-FROM-MEQ.md` §2** | **the one thing still open.** Auxiliary globally-coupled unknowns — a request for a capability, not a finding: nothing in that tree has been measured wrong by it, and it says so |

Every other report MEQ has filed is closed, and they are not listed: a closed one
has left a fix with a test on it or a number under an `M-nn` anchor, and neither
needs a ledger entry. `CLAUDE.md` has the two operational facts about that tree
that do outlast a report.

**What is NOT filed is the local-solve seed** — item 5. It was found from MEQ's
side and nothing has been written into that tree about it, beyond one paragraph
inside the coupling request using it as the argument for why a caller should not
have to difference a condensed residual.

**And the coupling request is a request for a capability, not a finding.**
Nothing in that tree has been measured to be wrong by it; §2 and §3 of it are
things that do not exist rather than things that misbehave, and it says so.

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

## 4. Element-local parallelism — MFEM, outstanding

`../mfem-hdg-dev/doc/HDG-ELEMENT-LOCAL-PARALLELISM.md`, on
`gf-hdg-linearise-first`. Twelve sequential element loops in
`darcyhybridization.cpp` over work independent by construction.

**IT IS NOT "PROBABLY THE LARGEST PERFORMANCE WIN AVAILABLE", WHICH IS WHAT
THIS FILE SAID, AND UPSTREAM HAS MEASURED IT.** `d7ea90a538`, *"Element-local
parallelism in NPC: the structural win is worth 6%"*, breaks an NPC step down
per phase over six steps at four `(n, k)`: integrator-bound loops **59–63%**,
integrator-free loops **5.4–6.2%**, trace solve **31–35%**. Their conclusion is
that the loops threadable *today* are under 6% of the step and flat in both mesh
size and order, so Amdahl caps any gain there until the integrators are made
thread-safe. MEQ should stop quoting this as the big win.

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
THIS FILE SAID OTHERWISE.** It claimed §2 of `HDG-BEM-COUPLING-FROM-MEQ.md` — two
rectangular integrators — was blocking. It is not, and the reason is the NPC
port: the datum's data half is an **essential trace value**, not a weak form, so
its block is `ProjectBdrCoefficient` against the `PathTraceCoefficient` that
already exists; and the transmission block is reachable from `TransferPath::
Endpoint` and `ElementExtension::TransformBack`, both public, in about forty
lines of MEQ. §3, auxiliary globally-coupled unknowns, is an optimisation over
`N + 2` backsolves against one factorisation — the cost MEQ's `ψ_ax` border
already pays. Plan §6.4 is the per-stage table.

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

**WHAT IS LEFT IS THE FIXTURE.** It is a rectangle with one closed plasma rather
than a machine, and nothing has yet driven `q` on the limited tokamak, where the
support moves and `ψ_bnd` is an unknown.

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

### 12.2 Locate the plasma edge with `ψ*` — and this one is separate

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
