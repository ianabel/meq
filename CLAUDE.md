# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Do not append session logs, temporary notes, or todo lists to this file

MEQ (the Maryland Equilibrium Solver) computes axisymmetric plasma equilibria by
solving the **Grad–Shafranov** equation with a hybridizable discontinuous
Galerkin (HDG) discretisation built on MFEM.

`refs/Refs.md` indexes the papers behind the numerics, with the doi to fetch each
from. The two marked ✔ there are not background reading — they *are* the method,
and `src/meq` is an implementation of them.

**THIS FILE IS AN INDEX AND THREE COMPANIONS CARRY THE CAMPAIGNS.** It was 7,600
lines; what stayed is what a maintainer reads every session. Read the companion
for the campaign you are working on and not before:

| | |
|---|---|
| **`CLAUDE_HDGGS.md`** | the equation, the discretisation and its sign conventions, the post-processing, the NPC port, Newton and the bordered Newton, the analytic fixture ladder, and the linear solves |
| **`CLAUDE_FB.md`** | free boundary, FB-A to FB-7, and the `freegs4e` benchmark |
| **`CLAUDE_INVERSION.md`** | `ψ(R, z)` to `R(Ψ, l)` — critical points, the tracer, the averages, the disc basis, the `(Ψ, θ)` output, the `q(ψ)`-driven solve |
| **`CLAUDE_FLOW.md`** | toroidal flow, FL-0 to FL-8 — RoPP (136), the species closure, and Maschke & Perrin |

*The campaigns are their own files* below has the longer form, and the three
rules that survive the split.

**`MEASUREMENTS.md` HOLDS THE PUBLISHED TABLES AND THIS FILE POINTS AT THEM.**
Every measurement table lives there under a stable anchor,
and every place one stood now carries a `→ **[M-nn](MEASUREMENTS.md#m-nn)**`
line naming what it holds. The *argument* each table supports stayed here, which
is the half a maintainer reads every time; the digits moved, which is the half
they read occasionally. **Do not renumber the anchors.** What did NOT move is
the reference tables — the output formats, the stages, the MFEM options, the
branch list, the merge resolutions, the solver table, the fixture ladder — since
those are things you look up rather than results you check.

## Status: the solver works, and MEQ is a program

**MEQ solves the semi-linear Grad–Shafranov equation by Newton** `k+1` in both `ψ` 
and `q` for `k = 1,2,3` over four dyadic meshes, against an exact Solov'ev equilibrium on the linear path and
HDG-GS-1's Example 5 manufactured solution on the non-linear one, with Newton
converging quadratically. `tests/convergence/SolovievConvergence.cpp` and
`NewtonConvergence.cpp` are the acceptance criteria and print the tables.

* `src/meq/GradShafranov.{hpp,cpp}` is the solver — `DarcyForm` with
  `EnableHybridization`. `Config`, `Profiles`, `Source`, `BoundaryShape`,
  `Estimator`, `Field`, `Sampler`, `WarmStart` and `Output` are all in
  `meq_core` and all under the `naming` check.
* `src/meq/Solution.hpp` is **gone**, not ported: it wrapped hand-rolled block
  offsets and a second set of FE spaces around what `DarcyForm` now owns.
  `v0-legacy` has the original.
* `Estimator.hpp` called `GridFunction::GetValueFacet`, which 4.9.1 does not
  have; it is replaced by `traceFes->GetFaceElement(f)` + `GetFaceVDofs()` +
  `CalcShape()`, exploiting `DG_Interface_FECollection`'s `VALUE` map type — the
  same pattern `estimators_hdg.cpp` uses.
* **`ψ_ax` is an unknown of the non-linear system**, not an input.
  `meq::NormalisedSource` is the interface and
  `setSource( NormalisedSource &, double )` closes it by a bordered Newton, with
  the two non-local terms in the border. `HighBetaConvergence` is the acceptance
  criterion and it is green. `meq::NormalisedMHDSource` is the production source
  built on two `Profile`s and is reachable from a TOML file as `[source] Type =
  "mhd"` with `Normalised = true`, through `meq::makeNormalisedSource`.

**MEQ is runnable.** `apps/meq.cpp` is the driver, `MEQ_BUILD_APP` defaults
`ON`, and `meq config.toml` parses, builds the mesh and source, solves — with
the adaptive loop if asked — and writes the same equilibrium **three times**,
with exit codes 0/1/2/3 as `docs/running.rst` specifies — plus, on request, a
fourth file that is not the equilibrium at all but a **reduction** of it:

| | |
|---|---|
| `.mesh`, `_psi.gf`, `_grad_psi.gf` | exact — every P_k coefficient. GLVis, and MEQ's own restart |
| `<stem>/<stem>.pvd` + `Cycle000000/` | VTK, **at the degree of the field it draws**, which is `k+1` since `ψ*` became that field — `apps/meq.cpp` passes `polynomialDegree + 1` at four sites. ParaView, VisIt |
| `.nc` | ψ and **B** on a uniform `(R, Z)` grid. Lossy, and the interchange format |
| `_surfaces.nc` | the flux surfaces and the flux-surface averages, against `ρ = √Ψ_N`. **Only with `[output] FluxSurfaces`**, because it costs a trace and a fit per surface and is the one output that can be impossible on a run that solved perfectly well. IN-6, and `CLAUDE_INVERSION.md` |

**AND `meq-run` MAKES THE MESH FIRST, FROM THE SAME FILE.** `[mesh.generate]`
names a generator and its geometry, and `build/meq-run config.toml` runs it and
then execs `meq`. It is a **script**, generated by `file(GENERATE)` into the
build directory beside the binary it runs, and that is a decision rather than a
stage: `meq` links MFEM and not gmsh, so something outside the solver has to run
the generator.

**THE WRAPPER DOES NOT PARSE THE CONFIGURATION, AND THAT IS THE WHOLE DESIGN.**
There is one reader of MEQ's schema and it is `meq::Configuration`; a wrapper
with its own opinion about what `[mesh.generate]` means is how a format grows
two readers that disagree and never say so. So `meq --mesh-command` prints the
generator's argv — shell-quoted, consumed by `shlex.split` — or **nothing at
all with exit 0** when the file names a mesh somebody else made, which is what
lets the wrapper tell "nothing to do" from "I could not read that". Every number
in that command came out of the same parser the solve will use.

**THE COILS ARE DERIVED FROM `[[coils]]` AND THAT IS THE POINT OF THE BLOCK.**
`examples/diverted-tokamak.toml` writes each conductor **twice** — as
`--coil 0.95 -1.15 0.10 0.10` in a header comment and as
`CentreR = 1.00, HalfWidth = 0.05` in a `[[coils]]` block — in two different
conventions, agreeing because somebody kept them in step. `halfdisc.py`'s
`10 + i` attributes are in command-line order *precisely* so the two lists can be
one, and the driver converts centre+half-extents to corner+extents when it
prints. It emits the **exact** doubles, so the mesh aligns to the rectangle
`meq::Coil`'s quadrature uses to the ulp — `-1.1500000000000001`, not the
comment's `-1.15`.

**AND `meq` REFUSES `[mesh.generate]` WITHOUT `--mesh-ready`.** MEQ cannot make
that mesh, so it cannot check it either, and the failure it would otherwise walk
into is the quiet kind: edit a coil, forget to re-mesh, and the run converges at
full order to the machine the PREVIOUS mesh described with every printed number
looking as it should. `meq-run` writes the generator command beside the mesh as
`<mesh>.meq-mesh` and skips the meshing when it has not changed — so editing a
coil re-meshes and editing `PolynomialDegree` does not, which is why the stamp
holds the COMMAND rather than the config file's mtime.

**PUTTING THE MESH'S GEOMETRY IN THE FILE MAKES THREE SOLVE-TIME FAILURES INTO
PARSE ERRORS**, each of which otherwise costs a mesh generation to find: a
`[boundary.limiter] SurfaceAttribute` naming a region the mesh will not carry
(at the solve that reads *"psi_bnd = max psi_h over the empty set"*), a limiter
circle through the axis, and — the one that is easy to get backwards — a
`[boundary.exterior] Radius` that does not fit strictly inside the disc.
**`[mesh.generate] Radius` IS NOT GAMMA**: `D_h` is cut FROM the generated mesh
as the elements inside `[boundary.exterior] Radius`, so the arc gmsh draws is
the background's outer edge and Gamma is the smaller semicircle within it.
`examples/diverted-tokamak-generated.toml` is the worked example and
`theDriverMeshesTheMachineItSolves` is the acceptance; it reaches the committed
mesh's answer to **1.6e-06 m in the X-point and 3.1e-06 in `psi_ax`**, on a mesh
gmsh rebuilt from the TOML at 2854 elements against the committed 2870.

**The VTK is written high-order deliberately**, and it is the one thing in that
row that can be silently wrong: VTK's native cells are linear, so the default
path draws a `k = 3` solution as though it were `k = 1` — a picture that looks
like a coarse mesh rather than like a bug.
`OutputConvergence::theVtkFilesCarryTheHighOrderSolution` asserts the point
count against the vertex count for exactly that reason, and reads **320 points
against 25 vertices** at `k = 3`.

**THE CURVED PATH LEAVES A BAND BETWEEN `Γ_h` AND `Γ`, AND BOTH OUTPUT FORMATS
NOW DEAL WITH IT — DIFFERENTLY, BECAUSE THEY HAVE TO.** `Ω_h` is the union of
background elements lying *inside* `Γ`, so `Γ_h` is inscribed and there is a
band `O(h)` wide that is inside the plasma and outside the mesh.

* **The `.nc` grid continues into it USING THE FLUX**, which is the mixed
  method paying off. `q` is computed at the *same*
  order as ψ and `∇̄ψ = R q`, so a node `p` outside the mesh is reached from its
  foot `x₀` on `Γ_h` as `ψ(x₀) + R₀ q(x₀)·(p − x₀)` — **nothing is ever
  evaluated outside an element**.

  **The obvious alternative was implemented first and does NOT work -- polynomial fits outside their domain are unbounded,
  even if you blend towards the known value at Gamma**

  **`B` GETS THE BAND TOO** Only `samplePotentialWithFlux()` applied
  the Taylor step; `sample()`, `sampleComponent()` and `sampleCoefficient()`
  never read `offsetR`/`offsetZ` at all, so for a band node they return the
  value **at the foot on `Γ_h`** — and `apps/meq.cpp` sent `ψ` through the first
  and `B` through `sampleComponent`. Measured on `miller-curved.nc`:
  `extrapolated_nodes = 1667` of 129×129 = 16,641, so about **one node in ten
  carried a piecewise-constant `B`**, `O(h)`, in a field `q` itself resolves at
  `k+1` and that `ψ` receives at `O(h²)` in the same band. `B = (−q_z, +q_r)` is
  a pure relabelling, so nothing softened it.

  **`sampleComponentWithGradient()` is the fix**, and it takes the same Taylor
  step from the same foot, using the field's own `∇u` — which is read *inside*
  the element, so the guarantee that nothing is evaluated outside one survives.
  `theBandVectorContinuesAtItsGradientsOrder` measures both sides on a quadratic
  its own space represents exactly, so the band error is the truncation alone:
  **rate 2.20 with the gradient against 0.92 reading the foot, and 124× smaller
  at `n = 16`**, interior nodes exact to 2.2e-15.

  **IT DOES NOT REACH `ψ`'s ORDER AND THE GAP IS STRUCTURAL, NOT A SHORTCUT.**
  `ψ` is continued with `q`, a **solved** variable carrying the potential's own
  order — that is the mixed method paying off. There is no solved variable for
  `∇q`: differentiating an L2 field of degree `k` leaves `k−1`, so this is
  `O(h²)` at every `k`. A full order better than what it replaced, and not the
  same thing. **The route to better is known and was not taken**: `div q = −F/R`
  is the equation being solved and `∂_r q_z − ∂_z q_r = −q_z/R` follows from
  `R q = ∇̄ψ`, which pin two of `∇q`'s four entries exactly — but they leave the
  symmetric traceless part still differentiated, so they buy *structure* rather
  than an order, at the cost of plumbing the source into `GridSampler`.

  **AND THE COST IS MEASURED, NOT ESTIMATED.** Wiring the rotating output gave
  a free controlled experiment: `n_s` is algebraic in `(R, ψ)`, so the same
  closed form can be evaluated both ways over the same 356 band nodes. Against
  the node's own `R` it is **5.13e-07** wrong; against the foot's `R`, which is
  what `sampleCoefficient` hands you, it is **8.57e-02** — a factor of
  **1.7e5**. `n_s` gets the radius wrong twice over, since the exponent of (96)
  carries `r²`; `B` gets it wrong once, being a relabelling of `q`. The
  experiment is `theRotatingFieldsUseTheNodesOwnRadius` in
  `OutputConvergence.cpp`, and it exists because the trap was met and dodged for
  the new fields while `B` was left in it two lines above.

  **AND THE FILE NOW SAYS WHICH NODES THOSE ARE**, which was the other half of
  the defect and is not fixed by making `B` better. `extendOutward()` counts band
  nodes as found, so `located()` is true and `inside` says 1 there — correctly,
  since the node is in the plasma and carries real data — and
  `theMaskAgreesWithTheData` asserts mask-matches-non-NaN, which is precisely the
  property that made the band look trustworthy. `extrapolated_nodes` was a
  **count, not a mask**, so nothing downstream could tell *which*. The `.nc`
  carries a second `byte extrapolated( Z, R )` beside `inside`, from
  `GridSampler::wasExtended()`; on `miller-curved` it sums to 1667, agreeing with
  the attribute, and is a strict subset of `inside`. Drop those nodes before
  computing an error norm or differencing two runs. **This is the interchange
  format, and `B` is what most readers open it for**, which is why it was an item
  rather than a footnote.

  **AND WRITING IT COSTS ALMOST NOTHING; LOCATING THE GRID IN THE MESH DID.**
  The NetCDF write is one `putVar` per variable — a bulk call, 5 ms — and what
  stood beside it was `meq::GridSampler`'s constructor, inverting the element map
  through `ElementTransformation::TransformBack()` at **1.45 us a call** for a
  straight-sided triangle whose map is affine and whose inverse is a 2x2 solve.
  `Mesh::GetNodes() == nullptr` is exactly that condition and a curved mesh keeps
  the Newton route. **69x at 129² and 81x at 513²**, with every node landing in
  the same element and the sampled fields moving by 1e-14 relative, which is
  round-off in a different order of operations. The sampling passes then group
  the located nodes **by element** rather than walking them in grid order, so the
  dof list, the gathered coefficients and the shape vector — three heap
  allocations per node through `GridFunction::GetValue()` — are fetched once per
  element instead; that is 1.6x, and what is left is `CalcShape`, which is
  genuinely per point. **The grouping is also what would make the passes
  parallel**, each element writing a disjoint set of node indices, which a
  node-major loop over a scattered element map cannot offer.

  → **[M-92](MEASUREMENTS.md#m-92)** — the output stage per writer, at two grid
  sizes · what the two fixes are worth · what does not move

  **THE BIGGEST ITEM IN THAT PHASE IS NOT GRID-SHAPED AT ALL**, and anyone timing
  the output will meet it first: `postProcess()` is **0.141 s at 6.31 cores**
  → **[M-172](MEASUREMENTS.md#m-172)**, and it does not move with the grid
  because it has nothing to do with it. It is `DarcyForm::Reconstruct()`
  re-assembling four integrators at the enriched order per element, which is the
  price of reporting `ψ*` rather than `ψ_h`.

  **IT IS THREADED UPSTREAM, AND THAT WAS THE LARGEST SINGLE-THREADED ITEM IN A
  MEQ RUN.** At 0.632 s and 1.00 cores it was 18.5% of the DIII-D run
  ([M-138](MEASUREMENTS.md#m-138)) — a leg that had read 0.61 s when it was 8%,
  its share growing while the leg itself sat still, because the NPC traversal and
  MEQ's own serial legs were threaded around it. **A fixed-cost serial leg gets
  more expensive every time somebody threads something else**, which is the
  general form and is why it was worth asking upstream for rather than a footnote
  to the output stage. It is **4.60× on the leg and `1.177×` on the whole run**,
  against the `1.19×` an Amdahl sum off M-138's share predicted.

  **AND IT IS HOST THREADING, SO A DEVICE RUN GIVES EVERY BIT OF IT BACK — WHICH
  IS NOW THE WHOLE OF WHAT A DEVICE RUN COSTS** → **[M-223](MEASUREMENTS.md#m-223)**.
  A configured `mfem::Device` and a host OpenMP region that touches MFEM memory
  are mutually exclusive (M-188), so `ThreadHostLoops()` returns false and this
  leg runs serially there. Measured by a thread count rather than argued: under
  `--device cuda` the output stage reads **0.90 and 0.97 s at
  `OMP_NUM_THREADS=8` and the same at 1** — flat — against a host leg that goes
  from 0.89 s at one thread to 0.37–0.51 at eight. With the solve leg now at
  parity between the two, this one leg **is** the 1.07–1.13× a device run loses
  on the wall. *Nobody should size the device work against a host-threading
  gain, and this is the leg that shows why.*
* **The `.vtu` bends the mesh onto `Γ`** — a curvature is installed and each
  boundary face is moved out. Since the VTK is already Lagrange cells this
  needed **nothing further from the format**; the two features composed.
  `theBoundaryBendsOntoTheTrueGamma` asserts the nodes land on the target to
  **1.1e-16**.

**Two things about that bending were got wrong on the way and are worth not
repeating.** Smoothing the displacement into the interior on the *vdof* graph
couples R to Z — they share one index range — so a radial displacement gets
averaged against a vertical one; measured, that made tangling **worse**, 25%
surviving where moving the boundary alone managed 50%. And backing the
displacement off **globally** costs every face the worst face's limit: per-node
backoff takes `miller-curved` from 50% to **96%** of the boundary reaching `Γ`.
The gap can exceed an element's own size, so some faces genuinely cannot reach,
and the driver reports the fraction that did.

**AN ADAPTIVE RUN ALSO WRITES `<stem>_cycles/`, ONE VTK FRAME PER CYCLE**, which
ParaView scrubs through as a time series — 97 → 254 → 342 → 449 elements over
`miller-adaptive`'s four cycles. It is a separate collection from the answer
because `<stem>` gets its boundary bent onto `Γ`, and doing that mid-loop would
hand the next refinement a geometry the estimator never saw.

**One frame per cycle is easy to get ALMOST right, and the near miss is
invisible.** Rebuilding the collection per frame puts every `Cycle` directory on
disk with the correct refined mesh, and leaves the `.pvd` index listing only the
last of them — so ParaView opens the file and shows a single frame, with no
error and all the data present. `ParaViewDataCollection` appends to its `.pvd`
and does not scan the directory, so the collection has to survive between frames
and be rebound with `SetMesh()`. `theAdaptiveSeriesIndexesEveryFrame` asserts on
the **index**, not the pieces, because the pieces were never what broke.

`tools/README.md` is the guide to which format goes with which reader.
`DriverAcceptance.cpp` asserts the driver reproduces the *library* on the same
configuration — 1.189e-16 over 15,360 dofs — rather than comparing against a
closed form, for the reason recorded beside `examples/soloviev-nstx.toml`.

**The curved boundary works through the driver.** `[boundary.shape]` builds the
shape, marks `D_h`, locates `Γ_h` and hands a `VertexConePath` to
`setExtension()` — GS-2's technique end to end from a config file;
`examples/miller-curved.toml` is the worked example. **Two checks guard it,
because the failure mode is quiet**: `Γ_h` carrying no transferred datum
silently imposes zero and still converges, so `theDriverSolvesOnACurvedBoundary`
pins the driver against the library (1.6e-16) *and* against that zero-datum
solve (2.7e-1, i.e. the transfer is doing real work).

**The adaptive loop works through the driver, on both paths.** solve →
post-process → estimate → mark → refine, stopping at `TargetError` or
`MaxIterations` and saying which. On the curved path it uses
`meq::AdaptiveDomain` — the companion mesh of GS-2 §3.3, without which
`dist(Γ_h, Γ)/h_loc` doubles every cycle and the transfer silently leaves the
regime it is analysed in — and calls `setTransferredBoundary()` automatically.
`examples/miller-adaptive.toml`: η monotone 4.73e-4 → 6.87e-5 over four cycles,
97 → 449 elements, 0 transfer paths widened. `theDriverRunsTheAdaptiveLoop` pins
it against the same loop driven through the library at **1.666e-16 over 2694
dofs**, and separately asserts that it refined, that it refined *adaptively*,
that η came down, and that assumption P.1 survived the graded `Γ_h`.

**The driver's non-linear path is a reactive ladder**: Newton, and on *observed*
failure `PicardThenNewton`, rebuilding the solver because a caught
`ErrorException` leaves one unusable. **Never predictive** — nothing may be
inferred from `F` about which solver to run, and two candidate detectors have
now been measured and killed. See `CLAUDE_HDGGS.md`, *Why MEQ's Newton struggles* and *There is no
cheap discriminator*.

**THE BORDERED PATH HAS ITS OWN SECOND RUNG, AND IT IS WHERE MEQ'S ONE BRANCH
CONTROL LIVES** → **[M-247](MEASUREMENTS.md#m-247)**, **[M-248](MEASUREMENTS.md#m-248)**.
`BorderedPicardThenNewton`, rebuilt, with the X-point held within
`0.05 × [boundary.exterior] Radius` of its seed — `setXPointRelocationCap()`,
BS-5. **It is reached two ways and both are observed**: the bordered Newton
failing, or converging to a flux with no interior extremum, an annulus. In the
second the first answer stands if the retry finds no axis either, the shared
source's state put back exactly. **The cap is a branch control and is legal
ONLY here**: from the first solve it loses MAST's and MAST-U's cold solves; on
the retry rung it closes A, A shaped, E and E shaped and moves nothing, taking
MEQ's cold table from five of fourteen with a located axis to nine. **And the
rung starts from the edge the first attempt started from** —
`freezePlasmaEdge()` on the shared source otherwise carried the failed attempt's
freeze into the retry, so one configuration closed or not depending on how the
attempt before it had failed. *A second attempt is only a second attempt if it
starts where the first did.*

**What the driver refuses rather than approximates**: `[boundary] Type =
"exact"` needs a closed form `meq::Source` does not carry, and exits 1 with an
explanation. There are also two more, and they are refusals about
the **build** rather than about the file — `[solver] AssemblyMode = "threaded"`
on an MFEM without OpenMP or thread safety, and `[solver] TraceSolver` naming a
package this build lacks. Both are checked once at startup, before a mesh
exists, so a run that cannot be honoured costs milliseconds.

**AND THERE ARE TWO MORE, AND THEY ARE THE FIRST THAT COST A
SOLVE.** Both are about the **answer** rather than the file or the build, and
neither can be recognised before the equilibrium exists — `ψ_bnd` is an unknown
of the bordered Newton, so what the profiles are evaluated at on the axis is not
known until it closes. **A `ψ_ax` that is not the flux at a magnetic axis**
refuses on a *positive* detection only, `checkAxis()` being one-sided and
missing rather than false-alarming; no interior extremum found at all stays a
warning, since a wall-hugging annulus is a real thing to look at — after a
diverted run has first been retried on the bordered rung, above. **A source that
does not vanish on the symmetry axis** refuses outright — `F/R` is `μ₀ j_φ`, so
that is an infinite current density on `R = 0`. The source one is **ordered
first**, because a bad `ψ_ax` is its consequence and the `ψ_ax` message's advice
(look at the guess, look at the mesh) is wrong when the pole is the cause. Both
exit 1 and write nothing; `docs/running.rst`'s gloss on code 1 is amended, since
it said *"nothing was attempted"*.

**AND cuDSS IS NONDETERMINISTIC RUN TO RUN, WITH NO USABLE SWITCH AT MEQ'S
COLUMN COUNT — SO IT CANNOT SIT ON EITHER SIDE OF A BIT-EXACTNESS CLAIM.**
cuDSS 0.8.0 leaves `CUDSS_CONFIG_DETERMINISTIC_MODE` at **0**: same matrix, same
right-hand sides, fresh handle, five repeats, **4 of 4 differ**, at 6.7e-16 /
1.8e-15 / 3.6e-15 for `K = 1 / 4 / 16`. **`SetDeterministic( true )` is not the
answer and must not be reached for.** It is supported for **one right-hand side
only**; the solve phase refuses 4 and 16 columns with
`CUDSS_STATUS_NOT_SUPPORTED`, and MEQ's bordered step is `N + 4` columns through
`ArrayMult` — *the one configuration MEQ would want determinism in is the one
where it cannot have it*. With `CUDSS_MTYPE_SYMMETRIC` the switch is
additionally **silently wrong**, returning an exactly zero solution that every
phase reports as success. Buying determinism back with column-by-column `Mult()`
gives up M-98's 1.34× and validates a configuration that is not the one shipped.
0.8.0 is the newest NVIDIA ships, so this is the library rather than something
about to be fixed. **The consequence is for the DEVICE plan and it is a large
one**: *"`--device cuda` reproduces `--device cpu` to every printed digit"* — the
instrument that root-caused M-79, M-130 and M-187 — is available only while the
trace solve is PARDISO, so `DEVICE-PLAN.md` §3.5 now orders cuDSS strictly last
among the residency items, beside single precision, as the second
instrument-destroying step. **MEQ is clear of the two wrong-answer traps and it
was checked rather than assumed**: it builds `NONSYMMETRIC` / `FULL` and calls
`SetMatrixSymType()` before `SetMatrixViewType()` at both call sites
(`GradShafranov.cpp:451`, `TraceSolverScaling.cpp:470`), which is the order
upstream's own order-dependence defect requires.

**AND `TraceSolver = "gmres"` IS A FOURTH CHOICE, EXPOSED TO BE MEASURED AND
MEASURED AS UNVIABLE** → **[M-179](MEASUREMENTS.md#m-179)**. GMRES with a
Gauss-Seidel smoother needs **247 Krylov iterations** on the easiest case in the
tree and **4,179** on `machine-f-diiid`, where the run does not finish in 1,800 s
against PARDISO's 3.0 s. It also gives up the two things the direct solvers are
chosen for — one factorisation amortised over the bordered step's fourteen
columns (M-98's 1.34×) and the symbolic reuse kept across Newton steps. **So an
iterative trace solve is closed rather than pending**, and Ginkgo — which is
otherwise well fitted, being the one MFEM solver package that is not MPI-gated
and that takes an `mfem::Device` directly — would have to close that gap with a
preconditioner rather than with hardware. **And the driver silently ran UMFPack
for the first measurement**, because `apps/meq.cpp:1634` initialises
`TS traceSolver = TS::UMFPack` before a `switch` with no `default`; see the
*Traps* entry.

**THE TWO PERFORMANCE KEYS ARE THE ONLY SOLVER KNOBS EXPOSED TO TOML, AND THE
LINE IS DRAWN WHERE IT IS FOR A REASON.** `AssemblyMode` and `TraceSolver`
cannot change the answer **where both converge** — the assembly modes are
asserted bit for bit at `MKL_NUM_THREADS=1` and the trace solvers agree to
1e-14 — so a file may choose them freely.

**THAT QUALIFICATION IS MEASURED AND IT IS NEW.** Making PARDISO the default
turned `PedestalConvergence`'s `andersonPicardReachesTheSameSolutionAsNewton`
**red**, and not by disagreeing: at `k = 1`, `n = 16` — the coarsest and
stiffest case in the suite, unglobalised, 500 iterations allowed — Newton
*reaches* a solution under UMFPack and does not under PARDISO, so the case's
precondition failed before any comparison happened. Isolated by rebuilding with
that one line changed and nothing else. Several hundred compounded 1e-15
differences land on different sides of a convergence boundary; the same binary
prints `Newton alone : FAILED in 60 iterations` for a neighbouring case under
**both** solvers, so unglobalised non-convergence there is the problem's
property and not the package's. See **[M-78](MEASUREMENTS.md#m-78)**.

**AND THE RED WAS THE TEST'S, NOT PARDISO'S.** `k = 1, n = 16` is the knife edge
`pedestalConvergenceIsAResolutionThreshold` spends forty lines establishing must
never be asserted on — whether that one element-local Newton converges is decided
by a threaded BLAS-3's summation order at O(1e-16), so the same source converges
in 42 iterations against one MFEM build and fails at the cap against another
differing by a single flag. The comparison case was taking it as a hard
precondition, which is the one thing its sibling forbids, and in the more brittle
direction: a `BOOST_TEST_REQUIRE` abandons the whole sweep on the first mesh. It
now skips a mesh either path fails, sweeps `{ 16, 24, 32, 48 }`, and requires
three of the four to produce a comparison — so the knife edge may fall either way
and two dropouts is still a finding. **Green, and the agreement is better than
the assertion needs**: best 9.068e-14, worst 3.894e-12. **Changing the default
did not change which problems MEQ can solve; it changed which mesh a test's
precondition landed on.** This is
still much weaker than the `Globalisation` hazard below, which reports one of
two discrete solutions 9.4% apart.
`Globalisation` and `NonlinearOrdering` are **not** exposed and must not be:
`CLAUDE_HDGGS.md`, *Should `PicardThenNewton` simply be the default?* measures three solve routes
reaching discrete solutions **9.4% apart** on an under-resolved mesh, which is
exactly the mesh an adaptive run starts from. A key that silently changes which
equilibrium is reported is not a performance key.

Two asymmetries in how the pair is handled, both deliberate. A wrong
**spelling** fails at parse (`Config` is MFEM-free and can only check the
string); an unavailable **choice** fails at startup (only the linked library
knows). And an `AssemblyMode` the file *asked for* is refused, while one merely
**inherited** from the default is downgraded to `Serial` — without that, the new
`Threaded` default would make every example in the tree fail on a stock MFEM,
which is most of them and is what CI builds.

**AND THE DRIVER REFUSES `TraceSolver = "cudss"` EVEN WHERE THE BUILD HAS IT.
IT IS WITHHELD, NOT UNIMPLEMENTED, AND THE REASON IS NOT THE MISSING
`mfem::Device`.**

A device solver is only worth having if the data **stays** on the device.
`../mfem-hdg-dev`'s `doc/HDG-DEVICE-OFFLOAD.md` — on the branch MEQ builds from,
and explicitly under construction — divides the element-local work into four
groups and measures their shares of an NPC step: local dense linear algebra
7–10%, **the integrators 46–53%**, the scatter 12–17%, and **the trace solve
26–31%**, which is what cuDSS is. It then gates the whole thing on the
integrators, in its own words: *"Two of the four groups are nearly free and
doing only those is worse than doing nothing. Groups 1 and 4 leave the
integrators on the host, so every iteration would copy the local blocks
host↔device around host-side integrator work — plausibly slower than staying on
the host throughout."*

**"Group 2 needs a partial-assembly rewrite and is not built" IS HALF TRUE AND
THE OTHER HALF IS STALE.** There is still no partial-assembly route, and MEQ is
triangles everywhere, so **that** economic case — sum factorisation — does not
apply to MEQ at any point. What upstream built *instead* is batched element and
face kernels reached through `CanBatch*` predicates, MEQ already passes most of
them, and upstream's own end-to-end measurement now reads **1.33× at 4096
elements** where the headline quoted above said *worse than doing nothing*.
`DEVICE-PLAN.md` §0 and §2 are the audit.

**The refusal stands anyway, and on MEQ's own numbers rather than on that
quotation.** M-126 puts the trace solve at **20.1%** of a threaded step, so a
config-file cuDSS is still exactly the group-4-alone case: a Device would be
configured, every Vector in the process would allocate through it, four fifths
of the step would still run on the host, and each Newton iteration would pay a
round trip for the one part that moved. **And upstream's crossover is not
MEQ's** — their fixture is quads, five of MEQ's eight shipped meshes are *below*
4096 elements, and MEQ's Newton is bordered with the border being host
arithmetic that would not go to a device at all.

**There is an immediate failure as well, and it is what made the refusal urgent
rather than only principled.** Found by exercising the key rather than reasoning
about it: without a Device, `CuDSSSolver` does not fall back — it reads host
pointers through the device-aware accessors and aborts with a raw
`CUDA error … cudaMemcpyDeviceToDevice … invalid argument`, a message with
nothing in it about the key that caused it. The **library** still offers cuDSS
and must: `TraceSolverScaling` constructs the Device first and the
`cuDSSTraceSolver` ctest is how its agreement with UMFPACK is checked at all.
Only the config-file route is closed, and
`theDriverRefusesASolverItCannotHonour` keeps it closed.

**AND THE SOLVE UNDER AN `mfem::Device` WAS BROKEN BY AN ALIAS LEAVING ITS
BASE'S VALIDITY FLAGS STALE, AT FOUR SITES, TWO OF THEM MEQ'S OWN.** The
symptom was a whole NPC solve mis-solving at `OMP_NUM_THREADS=1` — `example5`
reporting **0/0 Newton iterations on a case that needs four**, the potential
coming back identically zero, and all three trace solvers agreeing to four
figures, which is what said one common fault rather than three.

→ **[M-79](MEASUREMENTS.md#m-79)** — the chain · the faulting call · cpu, cuda and debug agreeing to every digit

**THE INSTRUMENT IS `mfem::Device( "debug" )` AND IT IS THE WHOLE STORY.** It
has device memory semantics with host arithmetic and `mprotect`s the host page,
so **a raw host read of a device-valid buffer is a named fault with a backtrace
rather than a wrong number** — which is exactly the shape of a symptom that
reads "identically zero". `NpcThreadScaling --device debug` costs nothing, the
harness already takes the string, and it faults on the first Newton residual.
Each fix then moves the fault **forward** to the next site, which is the only
thing that distinguishes a fix from a coincidence. Reach for it before a
debugger on anything device-shaped.

`GetAliasDevicePtr()` ends in `AliasProtect()` on the **base's** host range and
never touches the base's flags — `Memory::SyncAlias`'s own comment says so and
says what is owed. The base still reads `VALID_HOST`, the next `ReadWrite_`
computes `copy = !( flags & VALID_DEVICE )` = true, and `GetDevicePtr()` memcpys
h2d out of the page it just protected. Under CUDA nothing is protected, so it
silently copies the **stale host copy** over the device buffer the solve had
correctly written.

**MEQ'S HALF IS THE TRANSFERABLE PART**, because it was latent from the day it
was written and could not fail until a Device existed: **the constructor**
makes three long-lived `MakeRef` aliases — `darcyFlux`, `potentialGf`,
`traceGf` — onto the blocks of one `BlockVector solution`, and every one of them
is a separate alias registration that never hears about a device write. **They
are in the constructor and not in `buildForms()`**, which is worth saying
because `buildForms()` is the one that gets re-entered and is therefore where a
reader looks. `solve()`
therefore calls `solution.SyncToBlocks()` and a `SyncMemory()` per grid function
before the first host read, and they are inert with no Device configured. This
is the **third** time a library update or a new configuration has turned a
dormant MEQ contract into a failure in somebody else's code; see the `SubMesh`
parent and `SourceIntegrator`'s shared scratch under *Traps*.

**cuDSS is a bystander and always was.** UMFPack and PARDISO never touch the
device and failed identically. With the four sites fixed, `--device cuda` and
`--device debug` both reproduce the CPU answer to **every printed digit** on
the case that chain was walked with.

**AND THE CHAIN HAD TWO MORE LINKS, WHICH ONLY A BORDERED FREE-BOUNDARY CASE
COULD REACH.** `prepare()` seeds the iterate by writing THROUGH `potentialGf`
and `traceGf`; a device write through an alias marks the ALIAS device-valid,
`AliasProtect`s the base's host range, and leaves the base's flags unchanged —
so `formSystem()`'s `darcySolution` and `traceX`, built from those flags,
inherit the lie and `DarcyHybridization` reads the blocks on the host.
`solveWithNormalisation()` then reads `GetEssentialTrueDofs()` through
`Array<int>::operator[]`, which is as raw an accessor as `Vector`'s. Under
`debug` the first is a named fault; under CUDA nothing is protected, the border
columns are differenced against a state that is not the iterate, and the run
reports *"the bordered Jacobian is singular"* — **one fault, two faces, and the
CUDA face names a border rather than a memory.** `machine-f-diiid` now solves
under `--device cuda` in 2 Newton steps to the host's every digit.
→ **[M-130](MEASUREMENTS.md#m-130)**.

**AND THE CHAIN HAS A FIFTH LINK THAT COSTS NO ACCURACY AT ALL, WHICH IS WHY IT
SURVIVED THE OTHER FOUR.** Under any `mfem::Device` — `cuda` and `debug`
alike, and with every batched key set to `"serial"` so that **nothing runs on
the GPU** — `machine-f-diiid` took **31 / 16 / 4** Newton iterations per support
sweep against the host's **7 / 3 / 2**, and reported the same `psi_ax`, the same
`psi_bnd`, the same X-point and the same support counts at every sweep. The
answer was right and the work was 4×. → **[M-187](MEASUREMENTS.md#m-187)**.

**AND IT IS CURED — BY UPSTREAM'S `b_tr`, NOT BY MEQ'S `xLoc` SYNC, AND THAT
CORRECTION IS THE HALF WORTH KEEPING** → **[M-207](MEASUREMENTS.md#m-207)**,
**[M-210](MEASUREMENTS.md#m-210)**. `machine-f-diiid` under `--device cuda`
reads **7 / 3 / 2**, the host's to every digit at every sweep. M-207 credited
MEQ's own `xLoc.SyncFromBlocks()`; upstream reproduced the defect, **added that
exact sync, and the count did not move**, and the 2 × 2 MEQ then ran says
neither of MEQ's two changes holds it separately or together. The defect is
`DarcyNPCSolver`'s reused `b_tr`, written through an alias in `NPCReduce()`'s
host route and never carried back — MEQ's own `DeviceRows::upload()` shape, a
buffer refreshed in place being device-valid from its own previous use, which
is how upstream now describes it. *A cure measured after changing four things
attributes nothing.* The `xLoc` sync stays because the routine writes those
blocks and the flags should be honest; it is not load-bearing and its comment
says so. **The walk that named the mechanism is the one M-187 prescribed and
did not keep**:
`meq::ProbedSolver` sits where `mfem::NewtonSolver` keeps `c = J⁻¹r` to itself,
and on `mhd-rectangle` it reads the **same right-hand side to twelve digits
returning an increment 36× too small** — a wrong Jacobian, not a wrong
residual, which is what points at `GetGradient()`. **`mhd-rectangle` takes
`mfem::DarcyNPCOperator` instead of MEQ's own subclass, and that half is
upstream's and is CLOSED** — `ce744b84de`, in the library MEQ links, and it
reads **5 Newton steps under `--device debug` against 5 on the host**, every
printed line agreeing but the final residual at 5.3e-19.
→ **[M-223](MEASUREMENTS.md#m-223)**. **MEQ's report was right about the
symptom and wrong about the member**: it named `x_loc`, by analogy with its own
defect, and the fault was `b_tr`. *A reproduction is not a diagnosis*, and what
separated them was upstream running the same Newton two ways — through the
operator, and written out on the legs — in one process.

**THE SYMPTOM NAMES THE JACOBIAN AND NOT THE STATE**: the host is quadratic and
done in two, the device contracts by a constant factor of about ten per step,
and linear convergence at a fixed rate is what an inexact Jacobian looks like.
Nine of the ten **exterior columns** were identically zero on the device, so the
bordered Jacobian's exterior block was nine tenths missing — while the exterior
coupling stayed exact in the RESIDUAL, which is what drove the iterate to the
right answer anyway.

**`assembleExteriorColumns()` reuses ONE `mfem::Vector load` across the modes**,
zeroes it with `load = 0.0` (device-aware: `Write()` then a `forall`), hands it
to an `mfem::LinearForm` through `Update()` — which ALIASES it, with the alias's
own flags — assembles on the host into the alias, and then reads the BASE with
`load.HostRead()`, which is still holding the zeros. **Mode 0 survives**,
because a freshly constructed `Vector` has `UseDevice()` false and its first
`= 0.0` is a host loop; every mode after it reads zero. The repair is
`form.SyncAliasMemory( load )` before the read — MFEM's own call, and the same
one `formSystem()` makes for the solution's three aliases.

**THE LESSON IS ABOUT THE FUNNELS THEMSELVES.** Every sync this project added in
its raw-read sweep is a `HostRead()` **on the object about to be read**, which
is right when that object is the one that was written and useless when the write
went through an alias. Twelve sites were asked *"is this read funnelled"* and
twelve gave the right answer to the wrong question at this one. *A sync on the
base cannot see data that is still on the alias.*

**AND `mfem::Device( "debug" )` CANNOT FIND THIS ONE**, which is the first time
that instrument has come up short. It `mprotect`s a host page that is
device-valid, so it catches a raw access to the **wrong location**; here every
flag was honest about a location and wrong only about which copy was **newer**.
A 1170-second `debug` run reproduced the defect exactly and faulted zero times.
**What found it was a stage-by-stage host-side checksum** of the three solution
blocks and of the border's columns, printed under `--device cpu` and
`--device cuda` and diffed: the first stage that disagrees names the defect, and
here it was four stages past where anyone was looking.

**AND A SEVENTH IS NOT A DEFECT AT ALL BUT PURE WASTE, FOUND BY COUNTING
TRANSFERS RATHER THAN BY A CLOCK.** `mfem::Vector::operator=( 0.0 )` is
device-aware — `Write()` then a `forall` — so a host-only scratch buffer zeroed
that way acquires a device copy, and the next host read drags it straight back.
Two sites had it: `assembleExteriorColumns()`'s `load`, where
`mfem::LinearForm::Assemble()`'s own comment says *"The above operation is
executed on device because of UseDevice(). The first use of AddElementVector()
below will move it back to host"*; and the X-point's `xFluxField` /
`xPotentialField`, read only through host `GridFunction` methods. **Together
they were 12× MEQ's transfer count and 17× its bytes** — 36 transfers and
14.9 MB down to 3 and 0.888 MB on `machine-f-diiid` — for four
`UseDevice( false )` calls and no change to any answer.
→ **[M-200](MEASUREMENTS.md#m-200)**, **[M-202](MEASUREMENTS.md#m-202)**.

**NEITHER SHOWED IN ANY LEG PROFILE AND NEITHER WAS IN THE PLAN**, because both
are once per support sweep and a per-step timing cannot see them. **The
instrument that finds this class is a count of transfers** —
`meq::countedHostRead()`, which asks `Memory::DeviceIsValid()` before the read
and so separates a funnel call from a data movement. A count of *calls* points
at the wrong sites: MEQ's busiest funnel is the bordered `( N + 4 )²` row dot at
240 reads a step, and it transfers **zero**.

**AND A SIXTH LINK IS NOT AN ALIAS AT ALL: A DEVICE AND A HOST OpenMP REGION
THAT TOUCHES MFEM MEMORY ARE MUTUALLY EXCLUSIVE.** Under any configured
`mfem::Device`, MFEM registers buffers lazily — `Memory<T>::MakeAlias()`,
`SetDeviceMemoryType()` and the device arms of `Read`/`Write`/`ReadWrite` — into
`general/mem_manager.cpp`'s `maps->memories` and `maps->aliases`, which are bare
`std::unordered_map`s carrying **no mutex, no critical section and no atomic
anywhere in the file**. Any thread team that reaches one corrupts it.
→ **[M-188](MEASUREMENTS.md#m-188)**.

**THE INSTRUMENT IS A THREAD COUNT AND IT SEPARATES TWO CONTRIBUTORS.** Same
binary, same file, only `OMP_NUM_THREADS` moving: `machine-f-diiid` under
`--device cuda` exits 0 at 1 and aborts at 8, and under `--device debug` — no GPU
arithmetic in the process — it does the same. **With every batched key
`"serial"`, so that nothing of MFEM's is threaded at all, `debug` still
SEGFAULTS at 8 threads**, in `MemoryManager::Insert()` beneath
`meq::CriticalPointFinder::sweep()`'s own `GOMP_parallel` frame: each thread's
`GridFunction::GetVectorValue()` allocates an `Array<int>` of vdofs through
`Ordering::DofsToVDofs`, and under a Device that allocation registers. **`cuda`
survives that same loop by luck** — `GetHostMemoryType()` is `HOST` there and
leaves owners unregistered where `debug`'s `HOST_DEBUG` registers them.

**MEQ'S HALF IS ONE LINE IN ONE PLACE, AND IT COPIES MFEM'S OWN ANSWER.**
`DarcyHybridization::ThreadHostLoops()` opens
`if ( Device::Allows( Backend::DEVICE_MASK ) ) { return false; }`; MEQ gates all
of its OpenMP through `MEQ_OMP()` in `src/meq/Threading.hpp`, so the same
condition is now stated once there as `meq::hostLoopsAreSafe()`, carried into the
eight element loops by OpenMP's `if` clause. **`--device cpu` is not a device by
that test**, so every host measurement in `MEASUREMENTS.md` stands, and the host
path is bit-identical.

**THE OTHER CONTRIBUTOR IS MFEM'S AND MEQ CAN ONLY REFUSE IT.** With MEQ's loops
gated, `AssemblyMode = "threaded"` under `--device cuda` at 8 threads **still
aborts**, so the driver refuses that pair above one thread — exit 1 naming both
the key and the device, where it previously aborted three frames inside MFEM
naming neither, which is the same shape as the cuDSS refusal.
`theDriverRefusesThreadedAssemblyOnADevice` asserts the **pairing**, not the
refusal: refused at two threads, honoured at one, host path untouched.

**THE PART WORTH GENERALISING IS HOW THE SECOND PAIR SURVIVED M-79.** M-129
measured them, reproduced them twice, and wrote them down as a property of the
device path — while the sentence it closed with, *a device is an instrument for
finding unsynced reads*, described exactly what it was holding. A result that is
stable and reproducible is still a bug until something explains it, and the
difference was one `gdb` session breaking on `mprotect` with the faulting page
as the condition.

**AND AN EIGHTH IS A BUFFER REFRESHED IN PLACE, WHICH IS DEVICE-VALID FROM ITS
OWN PREVIOUS USE** → **[M-211](MEASUREMENTS.md#m-211)**. `meq::DeviceRows`'
mirror is re-uploaded once per Newton step and the `Read()` that ends the
upload leaves it on the device; `SetSize()` at an unchanged size reallocates
nothing, so the next fill is a raw host write into a protected page. **The
FIRST call is the only one that survives**, a fresh `mfem::Vector` having
`UseDevice()` false — M-187's *"mode 0 survives, every mode after it reads
zero"* in MEQ's own code. `HostWrite()` before each fill, and
`theMirrorSurvivesBeingRefreshed` pins it by uploading five times.

**AND IT CAME BACK, IN NEW CODE, CREATED BY THE OPTIMISATION THAT REMOVED THE
CONSTRUCTION** → **[M-226](MEASUREMENTS.md#m-226)**. `DEVICE-PLAN.md` §3.3's
phase B is *one upload*, and MEQ was re-uploading its five per-mesh device
constants once per rebuild of `meq::DeviceNPCOperator` — **sixteen times** on
the shipped machine case. Hoisting them into a solver-owned
`meq::DeviceConstants` is 80 uploads and 1.22 MB down to **6 and 0.113 MB** —
and it shipped a **neighbouring answer**, `psi_ax = 3.759845e-01` against
`3.759851e-01`, because the constructor built those arrays **fresh** every time
and a fresh `mfem::Array` has `UseDevice()` false, so its raw `operator[]` fill
was safe. A persistent object fills in place, `SetSize()` at an unchanged size
reallocates nothing, and the `Read()` after the fill copies nothing because MFEM
still believes the device copy is current. ***A fill that is safe in a
constructor is not safe in a refresh***, and the whole of phase B is turning the
first into the second — so expect this defect wherever work is hoisted out of a
loop onto a device. Every other printed line was healthy and the bordered
constraint was still satisfied at 1e-12, which is M-206's *a bordered row
certifies consistency and not correctness* a second time.

**THE INSTRUMENT THAT FOUND IT DID NOT EXIST THAT MORNING, AND THAT IS THE OTHER
HALF.** `meq::countedHostRead()` counts **downloads**; phase B is an **upload**,
so the acceptance `DEVICE-PLAN.md` §4 stage 4b rests on could not see the
traffic its own first phase is defined by. *An instrument that can see one
direction of the bus cannot measure a plan whose first phase is an upload.*
`meq::countedDeviceRead()` is the mirror — `HostIsValid() && !DeviceIsValid()`
before the `Read()` — and its per-site table read **`2 device reads, 1 upload`**
on the flag buffer, a fill that happened beside a transfer that did not. No
clock can see a copy that did not occur, and `Device( "debug" )` is blind to it
for M-206's recorded reason: every flag was honest about a location and wrong
only about which copy was newer.

**AND SPLITTING THE COUNT BY PHASE SAID THE PUBLISHED HEADLINE WAS NOT A
MEASUREMENT OF THE SOLVE.** `meq::PhaseScope` attributes each crossing to A/B/C/
D/E, and the *0.2 transfers a step* stage 4b reports turns out to be **three
transfers in phase A**, divided by a Newton step count — a rate for a phase they
do not occur in — while **phase C reads zero over 10,710 funnels**. *That zero is
not residency*: the counter reads the same `0` for a perfectly resident solve and
for one where nothing reaches the device at all, because 10,710 host reads that
move nothing means the state is host-valid throughout. Only the two directions
together separate them — uploads in B, zero in C, one download in D is residency;
zero everywhere is an empty device. `GradShafranovSolver::downloadSolution()` is
phase D and **moves nothing today**, for precisely the reason residency has not
arrived.

**AND THE ELEMENT BLOCKS ARE THE NINTH, AND THAT ONE IS UPSTREAM'S AND IS
FIXED** → **[M-214](MEASUREMENTS.md#m-214)**. `LocalFactorMode::Batched` — MEQ's
own default under a Device — left `Ct_data` device-only and
`DarcyHybridization::ReconstructTotalFlux()` reads it in a host `GOMP_parallel`
loop. **One key isolated it**: `Serial` gave exit 0 and the host's answer to
every digit, `Batched` aborted in post-processing with the solve already
correct. Upstream's `d720f02d62` gives the two batched NPC traversal kernels the
`HostRead()` every other kernel reading those blocks is followed by, and
`DeviceParity` is green with the key left on.

**AND MEQ'S REPORT CLAIMED CUDA WOULD SHIP A WRONG `ψ*`, WHICH IS FALSE AND WAS
CHEAP TO CHECK.** `_psistar.gf` on `cpu` against `cuda`, at the pre-fix install,
agrees to **6.94e-15** — round-off. A device `Read()` only ever ADDS
`VALID_DEVICE` and never clears `VALID_HOST`, `Ct_data` is host-built and
written by no kernel, and only the debug backend protects the page, so the
fault was the only symptom. *An argument that a defect COULD corrupt something
is not a measurement that it did*, and the distance between them was two runs
and a `cmp`.

**AND A SEVENTH LINK IS A COPY CONSTRUCTION, WHICH HAS NO ACCESSOR IN IT AT
ALL.** `mfem::Vector v2( v1 )` is a raw `memcpy` of `v1`'s host pointer, and
this project's raw-read sweep asked *is this read funnelled* of every
`operator()`, `GetData()` and `HostRead()` in `src/meq` — a copy construction is
a **declaration** and contains none of them. Three sites in
`solveWithNormalisation()` copy the bordered `unknown`; all three now funnel
first, and each fix moved `Device( "debug" )`'s fault forward, which is the only
thing that separates a fix from a coincidence.
`grep -nE 'mfem::(Vector|Array<[^>]*>) (const )?[a-zA-Z]+\( *[a-zA-Z]'` finds the
rest. The lesson is M-187's one level out: **the sweep's question was right and
its SURFACE was short.** → **[M-205](MEASUREMENTS.md#m-205)**.

**AND MEQ REPORTED TWO THINGS ABOUT `DarcyNPCOperator` AND WAS RIGHT ABOUT BOTH
SYMPTOMS AND WRONG ABOUT BOTH MECHANISMS. BOTH ARE FIXED, AT `ce744b84de`.**
MEQ read the step count as `Mult()` and `GetGradient()` writing `x_loc`'s blocks
in the idiom `bf27f5a928` fixed for `yb`. **Nothing in `DarcyHybridization`
reads `x_loc` as a whole** — both entry points take `x.GetBlock( 0 )` and
`x.GetBlock( 1 )` — so that sync is inert, which is checkable by reading two
functions and was not checked. The defect is `DarcyNPCSolver`'s reused `b_tr`,
written through an alias in `NPCReduce()`'s host route and never carried back
while the previous call's `b_tr -= r_tr` had left it device-valid.
**That is MEQ's own `DeviceRows::upload()` shape** — a buffer refreshed in place
is device-valid from its own previous use — and it is how upstream now describes
the class.

**AND THE SECOND WAS `Df_data += Df_lin_data`**, in
`DarcyHybridization::NonlinearGradientFinish()`, which MEQ named from a
`RelWithDebInfo` build because without debug info gdb collapses the inlined
template frames onto `GetGradient()` and names the wrong function. **It is a
view's flags rather than an array's**: the caller-integrand residual reads
`Df_lin_data` on the device through `DenseTensor::NewMemoryAndSize()`, and a
view is a **copy** of the `Memory` object rather than an alias, so the device
read is recorded in the copy's flags alone. *A sync on the object cannot see a
read taken through a copy of its bookkeeping* — the alias lesson one level over.
It needs a residual BETWEEN two gradients to reproduce, which is a Newton loop
and not what a single-gradient case does.

**AND `examples/q-driven.toml` WAS WRONG UNDER ANY CONFIGURED DEVICE BY A FACTOR
OF TWO, WHICH IS A SEPARATE DEFECT AND IS NOW FOUND** →
**[M-206](MEASUREMENTS.md#m-206)**. `psi_ax = 1.723374e-01`
against the host's `3.457380e-01`, the outer loop reporting *"converged in 0
KINSOL iterations … final norm 0"*, **and the inner bordered solve converging in
the same 4 Newton steps and satisfying its own constraint to 1e-16 on both** —
because a bordered row is satisfied by the state it was differenced at, so it
**certifies consistency and not correctness**. Nothing printed is red.

**IT IS NOT THE FAULT ABOVE, AND THAT IS THE CONTROL RATHER THAN THE READING.**
Forcing `makeDeviceNPCOperator()` to return null removes the fault outright and
leaves the answer wrong to every digit, so the two are two defects reachable by
one command. **`Device( "debug" )` cannot see this one** — it runs the whole
43-solve outer loop to completion without faulting — which is M-187's recorded
limit: it catches a raw access to the **wrong location** and is blind to honest
flags over a **stale copy**.

**THE STAGE-BY-STAGE CHECKSUM IS WHAT FOUND IT, AND IT IS KEPT THIS TIME.**
`meq::stageChecksum()` in `src/meq/StageChecksum.{hpp,cpp}` reports each
buffer's **validity flags before it reads** — the half a plain checksum lacks,
since a stage whose host copy is not live is the defect whether or not its
value is wrong — and prints under `MEQ_STAGE_CHECKSUM`. It named
`outerResidual.y` at `H- D+` with norm **0** where the host had 4.28.

**MEQ HANDS SUNDIALS HOST-FILLED VECTORS AND SUNDIALS READS THE DEVICE ARRAY.**
`mfem::SundialsNVector` binds BOTH of SUNDIALS' pointers into one `Memory`
(`linalg/sundials.cpp:465`), and `N_VSetDeviceArrayPointer` is what KINSOL's
CUDA kernels use. **Five seams, each hidden behind the one in front of it**:
`FixedPointResidual::Mult` filling `y` rawly (`‖F‖ = 0`, converged at iteration
0), `::GetGradient`'s `plus( j ) += step` after a device-aware `plus = x` (**a
zero finite-difference Jacobian, so `-nan` in three evaluations**),
`DenseJacobianSolver::Mult`'s `x = 0.0` plus a raw fill, `::sameAs`'s raw
comparison of the two iterates, and the starting guess and answer in
`solveForToroidalField`. *Fixing only the first replaced a wrong answer with
730,196 spinning NaN iterations — a partial fix here does not look like
progress, it looks like a new defect.*

**AND THE FOURTH HAS NO SYMPTOM OF ITS OWN, WHICH MAKES IT THE HARDEST.**
`sameAs()` decides whether the cached Jacobian is reused, so a stale comparison
answers *same* at a point the gradient does not belong to and the iteration
contracts **linearly with every printed norm looking healthy**. With the other
four in and this one out, `psi_ax` was already right to 6.6e-06 and the loop
took 30 iterations against the host's 12, stopping at its cap. *A raw read on a
control decision costs convergence rate rather than correctness.*

**`HostWrite()` IS NOT ENOUGH ON ITS OWN AND THAT IS THE DURABLE PART.** It
makes the flags honest and moves nothing; what is also needed is a `Read()`
*after* the fill, whose only job is the h2d copy into the array the consumer
will actually look at. *A sync fixes a lie about where the data is; it does not
fix the data being in the wrong place.*

**AND `--device cuda` NOW REPRODUCES `--device cpu` TO EVERY PRINTED DIGIT** —
12 KINSOL iterations, 43 map evaluations, final norm `1.94945e-06`, conditioning
`1.386e-02`, `psi_ax = 3.457380e-01`, with the inner constraint at round-off on
both sides. From an answer wrong by a factor of two to parity on every number
the driver prints.

**And `q-driven` is the only example in the tree that drives `solve()` from an
outer Newton**, so it is the only one that re-enters the solver's long-lived
state dozens of times with a Device configured; *a defect in re-entry needs a
case that re-enters.*

**AND cuDSS IS NOW THE TRACE SOLVER UNDER CUDA — ITEM 6, TAKEN DELIBERATELY —
AND BOTH OF THE PLAN'S PREDICTIONS ABOUT IT CAME OUT THE OTHER WAY**
→ **[M-229](MEASUREMENTS.md#m-229)**. The driver defaults to cuDSS when
`Backend::CUDA_MASK` is allowed and the file named no trace solver. **`CUDA_MASK`
and not `DEVICE_MASK`**, because `DEVICE_MASK` includes `DEBUG_DEVICE` and
`--device debug` is a Device with no GPU execution for a CUDA library to run on;
the two batched keys beside it are right to use `DEVICE_MASK` because they are
about memory and routes rather than about a vendor library executing.

**THE INSTRUMENT SURVIVED, WHICH §3.5 SAID IT WOULD NOT.** cuDSS is live and
nondeterministic — the bordered constraint reads `-6.465e-12`, `-6.466e-12`,
`-6.466e-12` over three runs where PARDISO is stable — and over **eight** repeats
every field `DeviceParity` compares is **identical**: 2 Newton steps,
`psi_ax = 3.759851e-01`, `psi_bnd = 7.081394e-02`, X-point
`( 1.200929, -0.999491 )`. The driver prints seven significant figures of a
quantity perturbed at the fifteenth. *The exact assertion stays exact and what
it rests on has changed* — determinism before, printed precision now — so it is
documented rather than pre-emptively widened, and the tolerance to reach for if
it ever flickers is the measured spread and not a convenient one.

**AND IT BOUGHT NO RESIDENCY, BECAUSE THE BLOCKER WAS NEVER THE TRACE SOLVER.**
[M-227](MEASUREMENTS.md#m-227) found the bordered columns `H+ D-` at the flush
and named the host solve; with cuDSS in they were still **168 of 168 `H+ D-`**,
because they are declared `mfem::Vector v( n )` and **a freshly constructed
`mfem::Vector` has `UseDevice()` FALSE**. *Item 4's gate was a predicate on data
that could not satisfy it under any trace solver, and it reported the same NO for
both reasons.* `UseDevice( true )` on them works — 126 of 168 arrive device-valid
— and costs **17× the bus traffic**, 0.888 MB to 15.1 MB, because the crossing
does not vanish but MOVES: `unknown.Add( -damping, y )` makes the ITERATE
device-valid and the next step's axis constraint, a host mesh search, drags the
whole state back. **Residency is a property of a CYCLE and not of a leg** — every
leg must move together or the one host leg converts every other leg's kernel into
an upload — and the last leg is the mesh search `DEVICE-PLAN.md` §3.4 already
calls the hardest part of the page. The arm is reverted and kept as a
measurement.

**AND THAT MESH SEARCH IS NOW ON THE DEVICE, AND THE PLAN WAS RIGHT ABOUT WHICH
ITEM WAS HARDEST AND WRONG ABOUT WHY** → **[M-230](MEASUREMENTS.md#m-230)**,
**[M-231](MEASUREMENTS.md#m-231)**, **[M-232](MEASUREMENTS.md#m-232)**.
`DEVICE-PLAN.md` §3.4 read the obstacle as host control flow -- a scalar
reduction feeding a branch -- and said *"none of it is add a kernel"*. It was
one level down and it WAS add a kernel, once the thing a kernel needs existed:
**`FiniteElement::CalcShape` is a host virtual taking an
`mfem::IntegrationPoint`**, so no device code could evaluate `q_h` at an
arbitrary reference point, and an element-local Newton is nothing but such
evaluations. `src/meq/PointLocatorDevice.cpp` had named that same blocker from
the other side for as long as it existed -- *"IT LOCATES, IT DOES NOT
INTERPOLATE"* -- and neither file had connected the two.

**`meq::DeviceBasis` IS THE MISSING PIECE AND IT BORROWS NOTHING FROM MFEM'S
INTERNALS.** The host asks the element for its own shape functions at its own
nodes and solves for the constant matrix re-expressing them in barycentric
monomials; that matrix depends on the geometry, the degree and the basis type
and on **nothing about the solution**, so it is §3.3's phase B in its purest
form. **A fixed quadrature rule could not have served**: `GetDofToQuad()` gives
`phi_i` at ONE `mfem::IntegrationRule`'s points, which is everything an integral
needs and nothing a search needs -- *a rule is a set of points; a search needs a
function*. **The ceiling is degree 4 and it is a MEASURED conditioning bound
rather than an array size**: outside the element, where this evaluator's caller
deliberately spends part of its Newton, the error grows twentyfold per degree
and reaches 4.8e-07 at degree 6, which a stopping target of `1e-13 max|q|`
cannot survive.

**THE ACCEPTANCE IS THE CANDIDATE LIST AND NOT THE ANSWER.**
`meq::axisCandidatesOnDevice()` roots every element -- a sweep, because on a
card M-100's 40x economy for a seeded search inverts -- and agrees with
`CriticalPointFinder::sweep()` to **2.2e-16 in both directions**, with `psi` at
a matched point to 7.8e-16 and the stopping scale exactly. *A test on the final
answer passes a port that found a different set and picked the same one*:
`machine-f-diiid` reports the same `psi_ax` whether or not the device list is
merged, and the merge is the half that was nearly left out.

**AND THE SAME 12.5 MB MOVED THREE TIMES, WHICH IS THE FINDING.** Forced open,
the axis constraint's own download disappears and the identical 12.5 MB
reappears at `refreshXPoint()`; port that too and it reappears at
`limiterValue()`. The iterate was crossing **once a step** and every host leg
after the first rode on it for nothing. ***A leg ported in isolation does not
remove a transfer, it RENAMES one***, and the count falls only when the LAST
host reader goes. So both routes are gated on
`DeviceIsValid() && !HostIsValid()` -- a predicate about the DATA, the only form
right both before and after residency arrives -- and are **inert today**, which
is the arrangement `meq::DeviceRows` has been in since M-208. **A per-site table
attributes a transfer to whoever asks FIRST, not to whoever needs it.**

**AND THE SEARCH IS OVER A NEIGHBOURHOOD, WHOSE INDEX WAS ALREADY IN THE TREE**
-> **[M-233](MEASUREMENTS.md#m-233)**, **[M-234](MEASUREMENTS.md#m-234)**.
Sweeping every element is 36x the host's warm ring search per call, so the
device roots a BLOCK around the previous axis instead -- and *finding* that
block is free, because `meq::PointLocator` already builds a uniform bin grid
over element bounding boxes in CSR form, at about one element per bin, for the
warm start's point location. `PointLocator::grid()` exposes it so there is one
grid rather than two uploads of one table. A warm attempt that finds nothing
widens to the whole mesh once, which is the ring walk's own structure and is
what makes a bounded search safe: *a neighbourhood that silently missed the
axis would hand the selection a competition its winner was not entered in.*

**AND THE WHOLE OF IT IS MEASURED NOW, IN BOTH FORMATS, AND THE PORT BEING
CORRECTLY SWITCHED OFF IS WHAT A DEVICE RUN LOSES ON**
-> **[M-235](MEASUREMENTS.md#m-235)**. Four cases, six arms, gated on a quiet
machine. The device wins the production bordered case in both formats --
**0.81x and 0.85x on `mf`'s solve**, and the fastest of all sixteen arms is
single + CUDA at 2.445 s against the double host's 3.727 -- and loses the big
fixed ones. **The reason is not the FP64 rate**: both CUDA arms print
`0 taken, 15 declined`, because the gate is a predicate about the DATA and is
inert until residency arrives, so the axis leg runs the HOST sweep with
`meq::hostLoopsAreSafe()` false (M-188) -- **1.00 cores against 1.97 and 2.61,
and 18 to 20 per cent of a device run's wall**. The 11.9 ms kernel sits beside
it, built and idle. *The leg the port was written for is the one the device is
losing on, and it is losing on it because the port is correctly off.*

**AND THE TEN HOST READERS OF THE ITERATE ARE PORTED, TO THE BIT, AND THE BUS
SAYS THEY WERE UNDER ONE PER CENT OF IT** -> **[M-237](MEASUREMENTS.md#m-237)**.
`src/meq/ResidentIterate.hpp` gives the five plasma element loops, the
transmission gather, the `psi_bnd` row and the three copy constructions device
routes behind the axis route's gate; forced, 126 of 126 evaluations equal their
host loops and the iterate is read on the host in phase D alone. **Three things
from it generalise.** *A second expression of the physics is bit-identical only
where it is the same expression*: `hermite10Prime` factored where `Profiles.cpp`
expands, an ulp in every derivative route, invisible to `SourceViewTests`
because its derivative envelope is four ulps of `psi` wide. *Removing the readers
does not open the gate*: unforced it reads `0 taken` until the bordered columns
are device-capable, and then every route opens by itself. And ***an instrument
that counts one project's funnels measures that project's readers, not the
bus***: the crossing counter read 0.888 MB on a run `nsys` puts near 950 MB, and
the largest single item -- **every backsolved column downloaded inside upstream's
`DarcyNPCSolver::ArrayMult()`, fourteen columns a Newton step and 140 MB over
`machine-f-diiid`'s twelve, on the default path**, because
`SyncFromBlocks()` matches each block to a base MEQ hands it host-valid -- was
found with a `gdb` breakpoint on `mfem::CuMemcpyDtoH`, one call further in than
the counter looks.

**AND THE BORDERED COLUMNS NOW STAY ON THE DEVICE, AND THE ROUTES OPEN WITHOUT
BEING FORCED** -> **[M-238](MEASUREMENTS.md#m-238)**. Device-to-host over a
`machine-f-diiid` run falls from 259 MB to 101 MB and host-to-device from 692 MB
to 543 MB, with 27 of 29 shipped examples at the host's every printed digit
under CUDA. Of the other two, `diverted-tokamak` was already wrong under CUDA at
`HEAD` -- the device integrand read a stale plasma support, fixed in
[M-239](MEASUREMENTS.md#m-239) -- and `rotating-normalised` locates a different
axis on the device, the axis routes' one stated semantic difference, which is
owed. **The transferable half is one sentence:
`BlockVector::SyncFromBlocks()` makes each block match its BASE**, so a kernel
writing the blocks of a vector whose base says host has its output downloaded
on the spot -- upstream's `ArrayMult()` did it to every backsolved column, and
MEQ's own residual operator did it to its output while its comment said the
call *carries the blocks up*. `Write()` on the base before the blocks are
written is the repair, plus `UseDevice( true )` on each block by name where
`Update()` kept blocks made before the base was device-capable --
`Device( "debug" )` faulted on exactly that. And the dense fill's device
contraction is bit-identical to `rowDot()` only because `psi_bnd`'s
`-limiterValue()` is a flag on the finished sum and not a negated weight:
*`-( +0 )` is `-0`, and a sum of negated terms from `+0` is `+0`.*

**AND AN UPLOAD OF HOST STATE IS A SNAPSHOT, SO IT GOES STALE WHEN THE HOST
STATE MOVES BETWEEN REFRESHES** -> **[M-239](MEASUREMENTS.md#m-239)**. The
device integrand's element flags were uploaded at the NPC operator's rebuild;
under a MOVING plasma support `refreshPlasmaComponent()` changes the mask at
every residual, the host integrand followed it and the device one did not --
so `diverted-tokamak` under any Device, `HEAD` included, converged to a
different equilibrium. Seven host perturbations agreeing exactly and four
device arms agreeing exactly with each other but not with the host is the
signature of a semantic difference rather than arithmetic, and it is worth
reading before any code. **`DeviceParity` could not see it because both its
cases freeze the support**, which is the coverage lesson: *a device test's
cases have to exercise every way host state can move, not only every leg.* The
same entry records that this example reported M-132's OTHER branch, 2.6% from
the freegs4e reference its header quotes, until `PlasmaSupportSweeps = 4` put it
on the benchmarked one; `DeviceParity` runs a copy with that key stripped to
keep the moving support covered.

**AND THE FIRST EVALUATION OF A SOLVE IS ON THE DEVICE, AND WHAT KEPT IT OFF
WAS A QUADRATIC MERGE THE HOST PAID TOO** -> **[M-240](MEASUREMENTS.md#m-240)**.
`ff3`'s axis leg reads **0.40 s against 22.1**, its solve 21.7 s against 43.3,
with the Newton history unmoved. Uploading the iterate before its first reader
was necessary and on its own made phase C WORSE: a `Type = "ramp"` guess leaves
the flux block zero, so every element reports a root at its seed, the device
route refused the list as over its ceiling, and the host sweep that replaced it
downloaded the iterate it had just been sent. `perf` then put the 22 seconds in
`CriticalPointFinder::mergeCandidates()`, which walked the kept list for every
candidate -- the rooting did not rank. **It is grid-indexed now and returns the
walk's list to the bit**, the walk itself being the test's oracle; the ceiling
is gone, the list being bounded by the mesh; and the nodal fallback is two exact
reductions on the device, because `mfem::ArgMaxReducer` gives a tie to whichever
joined LAST. *A cost attributed to a route is the route's only if the other
route does not pay it.* Four more of MEQ's own crossings went with it -- phase D
downloading the answer twice through its aliases, the output stage restoring
mesh nodes with a device-aware `=`, `fluxGf` negated on the device, and `rLoc`
arriving host-valid as upstream's reply said -- taking `machine-f-diiid` to
508 MB up and 78.6 down from 543 and 101. **And the residual's last uploads
ARE MEQ's load, although uploading `darcyRhs` added 11 MB and moved the count
not at all** -> **[M-243](MEASUREMENTS.md#m-243)**. The control uploaded the
BASE, and what `ResidualRowsBatched()` reads is the base's block VIEWS, which
`darcyRhs.Update()` re-makes host-only at every `prepare()`. Printing the
views' flags at the call settled it: host-only at exactly the re-assembly
leg's 16 calls of 28. *A fix that does not move its number measures the
attribution only if it reached the object the reader holds.* A view's validity
flags belong to the view.

**AND SINGLE PRECISION READ AS A 1.1x TO 2.4x LOSS WHILE EVERY ARITHMETIC LEG
IN IT WAS FASTER, WHICH IS THE SEVENTH CONSTANT OF M-184'S SPECIES AND THE
FIRST THAT COSTS ONLY TIME.** Per call the threaded solve legs are 10 to 27 per
cent FASTER in single, as FP32 should be; the loss is two SINGLE-THREADED host
legs, `setup` at 2.15x and `re-assembly` at 2.8x, on identical work -- same
mesh, same elements, same 0 paths widened, same 98 per cent reaching `Gamma`.
**The suspect was `extensionFdStep()` and it was wrong**; `perf` named
`meq::BoundaryShape::shiftedAngle`, which takes a `double` and is in the
MFEM-free half, so its arithmetic is identical in the two builds and its 2.99x
self time IS a call count. **`mfem::VertexConePath`'s default `tol_` is
`1e-13`** and `FirstCrossing()`'s bisection floors at `eps*length`, six orders
ABOVE it in single -- so the loop runs to `max_iter = 100` on every ray of every
vertex, where double stops at 39 halvings. *A fixed-count consumer and a
tolerance-driven consumer of one function, with the arithmetic held equal, is
what separated it*: `MarkLevelSetSubdomain` is flat at 1.25x and
`FirstCrossing` triples. **It misreports NOTHING** -- the crossing is located
to `eps`, better than the target asks -- so no assertion in the tree can see it
and only a clock can, which is how a suite at M-222's 62-of-63 parity carried
it. **AND THE REPAIR IS NOT A TOLERANCE -- IT IS THE BISECTION'S TERMINATION,
AND IT IS UPSTREAM'S.** MEQ built `meq::pathTolerance()` first and swept it:
every value loose enough to help is loose enough to MOVE AN ANSWER, because
this constant decides where the transfer path ends and the exterior datum is
assembled there. `64 eps` takes `ConductorSubtraction` red, `4 eps` takes
`DriverAcceptance` red at *"620 elements where the same loop by hand ends at
612"* -- a deterministic element count, not a coin toss -- and
`MachineFixedBoundary` with it. **MEQ's arm is reverted whole.** The fix that
costs nothing is one line of `FirstCrossing()`: a bisection in a finite format
should stop when the midpoint stops moving, `if ( tm <= ta || tm >= tb )`.
Simulated, float reaches the SAME answer in **20 iterations against the cap's
100**, and double's `1e-13` still fires first at 39, so **the double path is
bit-identical and no caller sees a change**. *The question was never how tight
to ask for; it was that the loop had no idea when it had finished.* And the
pair is carried by THREE path classes and not four -- `LevelSetPath`,
`VertexConePath` and the `TiltedPath` added by the 2026-10-01 path merge, which
all reach the one `FirstCrossing()`. `ClosestPointPath` is a closed form with no
bisection in it; MEQ's report said four and upstream corrected it. **THE FIX IS
UPSTREAM'S NOW, TAKEN AS WRITTEN** -- `238f719f50`, contained in
`meq-integration` at `5c2324e921`. **AND THEIR PIN IS THE SHARPER HALF, BECAUSE
IT RUNS IN DOUBLE**: a single-precision defect cannot be pinned by a
single-precision test in a tree whose single tests skip, so they ask a DOUBLE
for a tolerance no double bracket can reach -- `1e-30`, a cap of 10,000, and a
COUNT of level-set calls rather than a clock -- reading **10,002 before and
under 100 after**. *A format-specific defect is pinned by putting the other
format in the same position*, which is a cheaper instrument than a second build.
Reported as `HDG-PATH-BISECTION-TERMINATION-FROM-MEQ.md`.
**AND ON THE LIBRARY CARRYING IT THE LOSS IS GONE** ->
**[M-242](MEASUREMENTS.md#m-242)**: single is 0.81x to 0.90x of double per
Newton step on the three `fixed-f-diiid` meshes, host and CUDA alike, and its
SETUP is 0.69x to 0.81x of double's -- the direction the no-progress rule
predicts, a float bisection stopping at about 20 halvings where double's
`1e-13` stops at 39.

**AND THE KERNEL WAS 65x OFF ITS OWN ARITHMETIC, TWO THIRDS OF IT IN THE
REDUCTIONS AND A LOOP BOUND.** Both reductions were `mfem::forall( 1, ... )` --
the flux scale walking **1.23 million dofs on one CUDA core** -- and the loop
bounds were runtime, so nvcc could not unroll, every coefficient read stayed a
dynamic index and the per-thread arrays went to LOCAL MEMORY: **608 bytes of
stack frame with ZERO register spills**, which is what that pairing means.
*A stack frame with no spills is not register pressure; it is an array whose
bounds the compiler cannot see.* Halving `triangleBasisMaxDegree` as a
diagnostic accounted for the frame to the byte before any refactor. Templating
on the degree takes it to 0 bytes and 96 registers, and the sweep from **297 ms
to 11.9 ms at 102,840 elements** -- with the candidate lists bit-identical and
the stopping scale still 0.0 eps from `fluxScale()`'s. **It is NOT
`mfem::reduce`**, which takes a host reference and would put the one scalar
these kernels keep on the device back on the bus; `meq::DeviceRows` refuses the
same primitive for the opposite reason, its sum having an order to reproduce.

**AND REMOVING A FUNNEL EXPOSED A DORMANT DEFECT, WHICH IS NEW.**
`mfem::Device( "debug" )` aborted on the first run: `coldState( unknown )` is a
copy construction whose SOURCE is funnelled -- M-205's repair -- while the
DESTINATION is not. `Vector::Vector( Vector const & )` propagates `UseDevice()`
and copies device-to-device, so the copy is born device-valid with its host page
protected and the raw loops after it fault. ***Funnelling the source of a copy
does not funnel the copy.*** It was unreachable while the axis funnel ran first
on every path -- M-130's *ordering has protected it so far* for the third time,
and the first time **deleting** a funnel is what exposed the next one.

**AND THE REFUSAL OF `TraceSolver = "cudss"` WITHOUT A DEVICE STANDS ON A
NARROWER MEASUREMENT THAN IT USED TO, AND INVERTS FOR A RESIDENT PATH.** The
trade argument is measured rather than argued, and **the measurement has moved
twice**. [M-173](MEASUREMENTS.md#m-173) read `--device cuda` on
`machine-f-diiid` at **2.14× slower** with exactly one leg gaining — the trace
backsolve, 1.14×, which is cuDSS's — and four host-threaded legs collapsing to
1.00 cores. That arm ran `AssemblyMode::Batched`, whose host cost alone is
**1.74×** ([M-185](MEASUREMENTS.md#m-185)); with the driver's device default now
`Serial`, the NPC traversal on kernels and MEQ's integrand on the device, **the
solve leg is at parity and the residual 1.07–1.13× on the wall is the OUTPUT
stage** → **[M-223](MEASUREMENTS.md#m-223)**. So the group-4-alone case is no
longer what a device run of MEQ is; three of upstream's four groups reach the
device today. **It is wrong about the target**: on a fully device-resident solve the
integrators are on the device, nothing is copied per iteration, and cuDSS is
the only trace solver that does not force a round trip. `DEVICE-PLAN.md` §3.5
is the reversal and it is a change of DEFAULT, not of the refusal — the gate is
already *refuse without a device* rather than *refuse always*.

**THE MFEM HALF IS UPSTREAM'S AND IS LANDED** — `bf27f5a928` on
`gf-hdg-linearise-first`, found independently the same afternoon and merged into
`meq-integration`. `DarcyNPCOperator::Mult` and `DarcyNPCSolver::Mult` were the
only two places in `darcyhybridization.cpp` writing a `BlockVector`'s blocks
without syncing, an idiom the rest of that file uses at six other sites. **In
`DarcyNPCSolver::Mult` the ORDER is a third constraint**: `Neg()` acts on `xb`'s
own `Memory`, so the blocks are carried up to `xb` *before* it and `xb` up to `x`
*after* it.

**AND UPSTREAM'S FIX ALONE IS NOT SUFFICIENT, WHICH IS THE HALF ONLY MEQ CAN
SEE.** With it in and MEQ's own two links out, the iteration count is repaired —
0/0 becomes 2/2 — and `max|ψ_h|` is **still identically zero**. The `0/0` and the
zero potential are *different links of one chain*, so a fix for either leaves a
symptom standing, and upstream's own new regression asserts on the iteration
count. A library-side test cannot see a caller's aliases.

### `AssemblyMode::Threaded` is now a PROMISE MEQ makes, and it is load bearing

**MFEM REFUSES THREADED ASSEMBLY ON ANY PROBLEM WHOSE THREADED ELEMENT LOOP
WOULD EVALUATE AN INTEGRATOR, AND MEQ IS SUCH A PROBLEM.** `MultNL()` aborts,
naming `SetAssemblyMode()`, because the library can see exactly *which*
integrator objects it will call and nothing whatever about whether they are
reentrant — MFEM's convention is `#ifndef MFEM_THREAD_SAFE` around the scratch,
which is a **build switch** rather than a property of the class, and 67 of the
99 scratch-carrying classes do not use it at all.
`DarcyHybridization::SetIntegratorsThreadSafe()` is where a caller who has read
their integrators says so. Off by default, because the default decides what
happens to a caller who has not read any of this.

**MEQ CALLS IT, IN `buildForms()` BESIDE `SetAssemblyMode()`, AND THE AUDIT IS
THE JUSTIFICATION.** Two sets are covered, because
`DarcyForm::CanThreadAssembly()` gates the assembly element loops on the same
flag and those loops exist to evaluate the **linear** domain integrators:
`meq::SourceIntegrator` — MEQ's own, whose per-point scratch is already under
`#ifndef MFEM_THREAD_SAFE` and therefore does not exist in this build — and
`VectorMassIntegrator` and `VectorDivergenceIntegrator`, which upstream guarded
in the commit that put eight classes behind that switch. **What is not covered
is not reached**, and that was checked rather than assumed: see the
`HDGExtensionIntegrator` row of the upstream table above.

**IT IS NOT DECORATIVE, AND THAT IS FALSIFIED RATHER THAN ARGUED.** With the one
line commented out and nothing else changed, `machine-f-diiid` exits **2** with
*"AssemblyMode::Threaded is unsafe for this problem: this element loop evaluates
a nonlinear integrator per element or per face"* — on both the Newton and the
Picard-then-Newton path. So this file's standing warning that an MFEM update
carrying that abort would stop every threaded MEQ run was correct, and one line
is what answers it.

### Which integrators need device offload, and why kernels alone are not enough

**THE INVENTORY, TAKEN FROM THE INSTALLED HEADERS RATHER THAN REMEMBERED.** Every
integrator MEQ installs, and whether a device entry point exists for it today.
**The sites name `buildForms()` and not a line number, deliberately**: this
table carried line numbers for months, they were wrong by eight hundred lines
when anybody checked, and the replacements offered for them were wrong again
within three commits. `grep` for the integrator.

| integrator | where MEQ installs it | device path |
|---|---|---|
| `mfem::VectorMassIntegrator( radius )` | flux mass, domain — `buildForms()` | **`AssemblePA`, `AssembleMF`, `AssembleDiagonalPA`, `AddMultPA`** — and separately `HDGElementMassBatched()` under `AssemblyMode::Batched` |
| `mfem::VectorDivergenceIntegrator` | flux divergence, domain — `buildForms()` | **`AssemblePA`** |
| `mfem::TransposeIntegrator( mfem::DGNormalTraceIntegrator )` | flux divergence, interior + boundary faces — `buildForms()` | the wrapper has PA and the wrapped integrator has none — **and neither matters**, because under hybridization these are never assembled. They are MARKERS; see the measured trap under *A trap in that table* |
| `mfem::NormalTraceJumpIntegrator` | the hybridization constraint, handed to `EnableHybridization()` | `AssembleEAInteriorFaces` **only**: element assembly, interior faces, no `AddMultPA` and nothing for boundary faces |
| **`mfem::HDGDiffusionIntegrator`** | potential mass, interior **and** boundary faces — `buildForms()` | **none** |
| **`meq::SourceIntegrator`** | potential mass NON-LINEAR form, domain — `buildForms()` | **none, and it is the one MEQ owns** |
| **`mfem::HDGExtensionIntegrator`** | flux mass, boundary faces, curved and free-boundary only — `buildForms()` | **no kernel of its own, and it needs none.** `AssembleFluxMassBdrMatricesBatched()` batches the SCATTER and never dispatches on the integrator: it asks for an element block, checks its size, and everything after that is a kernel. This integrator satisfies that check. The pair exists only in `meq-integration` and had never been executed anywhere — **it agrees to 8.65e-15** → **[M-185](MEASUREMENTS.md#m-185)** |
| `mfem::DomainLFIntegrator`, `mfem::VectorBoundaryFluxLFIntegrator` | the two right-hand sides, in the two `assemble*Rhs()` helpers | none, and it does not matter: once per mesh, not once per step |

**TWO STRUCTURAL FACTS MATTER MORE THAN THE TABLE, AND EITHER ONE ALONE MAKES
"ADD KERNELS" THE WRONG PLAN.**

**The hybridized path never reaches a PA entry point, and that is fine because
PA is not the route.** `darcyform.hpp` says the hybridized assembly goes through
`BilinearForm::ComputeElementMatrix()` — the dense, per-element, host route —
so the three integrators that *do* have PA kernels are not using them in MEQ,
and kernels added to the other three would not be reached either. **The
conclusion this file drew from that was wrong**: it read *"what is missing is a
route"* and named upstream's group 2 partial-assembly rewrite. The route exists
and it is not PA.

**`AssemblyMode = "batched"` IS THE ROUTE, AND IT IS A TOML KEY MEQ ALREADY
EXPOSES.** `DarcyForm::Assemble()` calls `AssembleFluxMassMatricesBatched`,
`AssembleDivMatricesBatched`, `AssemblePotMassMatricesBatched` and
`AssembleFluxMassBdrMatricesBatched` **in front of** each element loop; each
returns false having done nothing unless `Batched` was asked for, and only then
does `ComputeElementMatrix()` happen. The comment quoted above describes the
**default**, not a wall. `HDGElementMassCanBatch()` admits `VectorMassIntegrator`
and `HDGElementDivCanBatch()` admits `VectorDivergenceIntegrator`, which is
exactly what `buildForms()` installs. **Measured on the production free-boundary
case: all four predicates yes, `psi_ax` to every digit of the host default**
→ **[M-185](MEASUREMENTS.md#m-185)**. Partial assembly stays out of the plan
for the separate reason that its economic case is sum factorisation and MEQ is
triangles.

**THE MISREADING IS WORTH KEEPING BECAUSE OF WHERE IT CAME FROM**: a doxygen
comment on a **member declaration**, believed without reading the function it
describes or the call site in front of it, and then carried as the single
largest dependency in `DEVICE-PLAN.md` for as long as that file existed. *A
comment is documentation of an intent; the call site is the fact.* Same species
as the `ComputeHDGFaceEnergy()` line-range error under *Traps*, and this one ran
longer.

**AND MEQ IS IN THE 15 THAT THE BATCHED FACE KERNEL TAKES, NOT THE 52 IT
REFUSES — THIS FILE SAID THE OPPOSITE AND THE REASONING WAS WRONG.** The wrong
step was inferring the constraint's ROUTE from which form the integrators are
added to. MEQ does put both HDG face integrators on the potential-mass
**non-linear** form, by the load-bearing decision at `GradShafranov.cpp:2975` —
but that is not what decides `c_bfi_p` against `c_nlfi_p`. Read out of
`darcyform.cpp` rather than reasoned about:

* `GetPotentialMassForm()` is called at `GradShafranov.cpp:3156` and **only**
  there, inside the `else` of `usesNonlinearForms()` — the Picard/Anderson path.
  On the Newton path `M_p` is therefore **null**.
* `EnableHybridization()` tests `if (M_p)` first (`darcyform.cpp:367`), falls
  past it, and reaches `else if (Mnl_p && FaceIntegratorsAreLinear(...))` at
  `:416`.
* `FaceIntegratorsAreLinear()` (`:289`) asks two things: is the problem
  non-linear for some other reason — `Mnl_p->GetDNFI()->Size() > 0`, which
  `meq::SourceIntegrator` satisfies — and is **every** face integrator, interior
  and boundary, a plain `BilinearFormIntegrator`.
  `class HDGDiffusionIntegrator : public BilinearFormIntegrator`
  (`bilininteg_hdg.hpp:843`), and both of MEQ's are.

So MEQ's face constraint is installed as **`c_bfi_p`, the linear route, assembled
once per solve rather than once per step** — which is exactly
`CanBatchPotFaceAssembly()`'s domain. The one remaining gate inside
`HDGFaceScatterCanBatch()` refuses a state-dependent stabilization hook
(`bilininteg_hdg.cpp:2658`); MEQ's is `meq::ConstantStabilization`, so it passes
that too.

**The transferable part**: "which form an integrator is added to" and "which
route its constraint takes" are different questions, and the second is decided by
a predicate in somebody else's translation unit. This file inferred one from the
other and got a requirement list wrong in MEQ's own favour.

**MEQ REACHES ALL THREE BATCHED PATHS AND TWO OF THEM ARE TOML KEYS.** This file
said the opposite — that `AssemblyMode` carried only `Serial` and `Threaded` and
that reaching MFEM's batched work needed a new enum value — and every part of
that is stale. `AssemblyMode::Batched`, `LocalFactorMode` and
`TraceAssemblyMode` all exist in `GradShafranovSolver` and all three are
settable from a file. What the measurement says is which of them is worth
having: **only `TraceAssemblyMode`**, at 7 to 9 per cent of the DIII-D solve, and
it is now the default in a file and in the class alike.
**`AssemblyMode::Batched` is a 15 to 27 per cent LOSS** on the host — **and
that is a measurement of the HOST, not a verdict on the route**: the batched
entry points *are* the device route, carrying the 25 `forall` launch sites in
`fem/darcy/` and `BatchedLinAlg`, whose `GPUBlasBatchedLinAlg` backend is linked
here. Batching costs a gather and buys nothing where the host was already
threaded; on a device it is the only way the work executes. *A cost is a
property of a target, not of a route.* `DEVICE-PLAN.md` §3.2. And
**`LocalFactorMode::Batched` IS NO LONGER A TRADE** — it bought a batched local
factorisation and paid the condensation cache, since `CanCacheCondensation()`
refused outright under it, and upstream has since handed
`FactorElementsBatched()` the cache's own buffer. Re-measured on the new
library: `condensationCacheTaken()` reads **yes** on that row, and the key is now
a small consistent gain — `ComputeH` lower in all three interleaved pairs for
about 4% — rather than a wash. **MEQ's HOST default stays `Serial`**, because 1% of
the wall is inside the serial arm's own spread and is not a reason to move a
default that has never been asserted bit-for-bit against the other route.

**ON A DEVICE IT IS THE DEFAULT, AND THAT IS A ROUTE DECISION RATHER THAN A
PERFORMANCE ONE.** `CanBatchLocalSolve()` opens on
`lfac_mode == LocalFactorMode::Batched`, so with the key off the element-local
dense algebra — upstream's group 1, and the one group `GPUBlasBatchedLinAlg`
already owns outright — **stays on the host under `--device`**. MEQ's device
runs were leaving it there by a driver line whose stated justification, the
condensation cache, had expired. `apps/meq.cpp` now turns it on under a Device
exactly as it turns on `AssemblyMode::Batched`, and only where the file said
nothing.

**AND THAT PAIRING IS SETTLED: `AssemblyMode` IS `Serial` UNDER A DEVICE** →
**[M-215](MEASUREMENTS.md#m-215)**, **[M-209](MEASUREMENTS.md#m-209)**. Upstream swept the three batched keys one at
a time, bitwise: `LocalFactorMode::Batched` and `TraceAssemblyMode::Batched`
**reproduce exactly** and `AssemblyMode::Batched` does **not**, at 3.5e-16 to
4.7e-16, because its face kernels accumulate `D` through `AtomicAdd` from three
or four faces per element. So the driver's own device default is what destroys
*"`--device cuda` reproduces `--device cpu`"* — the instrument `DEVICE-PLAN.md`
§3.5 orders cuDSS last to protect, and the reason residency is built in double
at all. **It is `Serial` there now**, which keeps the instrument and still moves the
local factorisation, the NPC traversal and the trace assembly — three of the
four groups, all bit-reproducible. What it gives up is the face kernels alone,
whose host measurement is a 15 to 27 per cent loss anyway. Two consecutive
`--device cuda` runs of `machine-f-diiid` now agree exactly where they read
`-6.465e-12` then `-6.466e-12`, and against the host with the keys held equal
every printed line agrees but one residual, at 5.3e-17 absolute.
**`Serial` and not simply dropping the override**: the host default is
`Threaded`, which the driver REFUSES under a Device above one thread (M-188),
so inheriting it would turn an unstated mode into a refusal. *An inherited
default has to be honourable in the configuration that inherits it.* Worth **1.54× on `ComputeH` under CUDA** and 1.97× on the host, with
`psi_ax` unmoved → **[M-185](MEASUREMENTS.md#m-185)**.
→ **[M-136](MEASUREMENTS.md#m-136)**.
→ **[M-99](MEASUREMENTS.md#m-99)** and **[M-101](MEASUREMENTS.md#m-101)**.

**THE ORDERED LIST, WEIGHTED BY THE LEG PROFILE RATHER THAN BY COUNT.** From
**[M-80](MEASUREMENTS.md#m-80)**:

1. **`meq::SourceIntegrator`** — the `F(R, z, ψ)` domain term, evaluated per
   quadrature point per element per residual **and** per Jacobian. **MEQ'S HALF
   OF THE DEVICE INTERFACE IS BUILT** → **[M-177](MEASUREMENTS.md#m-177)**:
   `meq::SourceView` is the POD, `sourceF()`/`sourceDFdPsi()` the
   `MFEM_HOST_DEVICE` evaluators, `sourceResidualAtPoint()` the element-local
   signature, and `SourceViewDevice.cpp` compiles the lot as CUDA. **It is the
   single `CanBatch*` predicate MEQ fails** — upstream puts MEQ at YES on the
   other seven — and the refusal is correct rather than an oversight, since
   this integrator genuinely is nonlinear; what was missing is an interface,
   and the other side of it is upstream's. **AND THE POD NOW COVERS THE
   PRODUCTION CONFIGURATION, WHICH IT DID NOT** → **[M-191](MEASUREMENTS.md#m-191)**:
   `makeSourceView()` refuses a coil-augmented source, which is every machine
   case in the tree, and that refusal is right for `f()` — the SUM — and
   inverted for `scaledF()`, the plasma term alone that the prescribed-current
   constraint integrates. `meq::makePlasmaSourceView()` is a second entry point
   rather than a widened cast, pinned to the bit through a real wrapper. **The POD is a SECOND expression of
   the physics and `SourceViewTests` pins it to the bit**; that test is
   load-bearing, for the reason M-177 gives. It is the
   only one on this list MEQ can write itself, and it is the leg that dominates
   MEQ's own bordered configuration: the residual leg reads **48.2% serial and
   23.1% threaded** on the high-beta case, because the border spends four
   residual evaluations per step. **`INTERPOLATORY-HDG-PLAN.md` is the plan for
   replacing it** with CCSZ interpolatory HDG — interpolate `F` into `P^{k+1}`
   rather than integrating it — and the library half is already built and
   linked, on `gf-interp-hdg-dev`. Read its falsifying experiment before its
   stages: the arithmetic there deflates the headline to 1.6–2.0× fewer `F`
   evaluations, `ComputeH()` is untouched, and the Amdahl ceiling off M-80 is
   23–28%.
2. ~~**`AssemblyMode::Batched` in MEQ's own enum**~~ — **built, and measured a
   15 to 27 per cent LOSS on the host.** It was listed here as the whole of what
   stood between MEQ and the face kernel that already exists, and it is; the
   kernel is now reachable and does not pay on this machine. It stays a TOML key
   for the device, where the trade may invert, and the item is closed rather
   than pending. → **[M-99](MEASUREMENTS.md#m-99)**.
3. ~~**A device route through `DarcyHybridization`**~~ — **IT EXISTS AND IT IS
   `AssemblyMode = "batched"`.** The four `Assemble*Batched()` routines sit in
   front of the element loops in `DarcyForm::Assemble()` and
   `ComputeElementMatrix()` is what happens when the key is off. Run on
   `machine-f-diiid` with all three batched keys: **face assembly yes, local
   factorisation yes, local solve yes, trace assembly yes**, and `psi_ax` to
   every digit of the host default → **[M-185](MEASUREMENTS.md#m-185)**. The
   one refusal left is `CanBatchLinearResidual()`, which is item 1 being
   genuinely nonlinear — an interface, not a route.
4. **`mfem::HDGExtensionIntegrator`** — needed only for the curved and
   free-boundary paths, which is to say for every problem MEQ actually exists to
   solve, and for none of the fixtures in M-80. **It does NOT disqualify the
   condensation cache, and this file said it did.** `Bnl_data` is written at one
   site, `darcyhybridization.cpp:6487` in `ConstructGrad()`, and only from
   `grad_Aup` — the (0,1) block of a **block non-linear** integrator's element
   gradient, i.e. a flux law depending on the potential. `HDGExtensionIntegrator`
   is a `BilinearFormIntegrator` (`extension_hdg.hpp:509`) on the flux mass
   **bilinear** form, so it lands in `A` and is assembled once. It makes `A`
   asymmetric — MEQ's own 5.4e-01 — but not solution dependent, and `A⁻¹Bᵀ` is
   constant whether or not `A` is symmetric. Same error as the one above:
   asymmetry was read as state dependence.
5. **Nothing for the two linear-form integrators.** They are assembled once per
   mesh; offloading them would buy a share of `prepare()`, which is outside the
   step entirely.

**And none of it is worth starting before the Device works at all** — M-79 — or
worth *timing* here afterwards in DOUBLE, since this machine's consumer FP64
runs at 1/32 to 1/64 of its FP32 rate where a datacentre part runs it at about
1/2. **That rate is why `MFEM_PRECISION = single` is a section of `DEVICE-PLAN.md`
rather than a footnote**: under a partial offload the loss is host
serialisation and the format hardly matters
([M-173](MEASUREMENTS.md#m-173)), but a fully device-resident solve on this card
is competing against eight host cores in the one format the part is crippled
for. MFEM supports it (`MFEM_PRECISION`, `CMakeLists.txt:213`) and cuDSS already
tracks `real_t`; what it costs is MEQ's own **3,913 `double`s against 13
`real_t`s** in `src/meq`, and a tolerance sweep of every assertion below FP32's
epsilon — the tightest in the tree sit at 1e-16. **The `override` surface is
SMALLER than a count of `double`s suggests, and that is the MFEM-free design
paying off**: the 33 `double`-returning overrides in `src/meq` are all on MEQ's
*own* `Source` / `Profile` / `Coil` hierarchies, which *Layout* keeps free of
MFEM so they stay `double` on both sides. What couples to `real_t` is the twelve
classes deriving from an MFEM type.

**AND IT WORKS** → **[M-184](MEASUREMENTS.md#m-184)**. MEQ solves
`machine-f-diiid` in single precision, on the host and under `--device cuda`, at
`psi_ax = 3.759862e-01` against the double build's `3.759851e-01` — **2.9e-6
relative, which is FP32's own floor** — with the magnetic axis located at the
same point to four decimals and `checkAxis()` accepting it.

**IT TOOK SIX CONSTANTS, AND THE PATTERN IS THE FINDING.**
`PointLocator::insideTolerance` (1e-8), MFEM's `MFEM_VERIFY( err < 1e-10 )`,
`normalisationDifferenceStep` (1e-5), MFEM's `ExtensionBoundaryQuadrature`
`fd_step` (1e-6), and `CriticalPointFinder`'s `tolerance` (1e-13) and
`jacobianStep` (1e-4). **Every one was right for double, wrong for single, and
silent** — the symptoms were a wrong number, an empty quadrature that reported
success, and a correct answer refused. **A tolerance that does not scale with
`real_t` is a build-configuration bug wearing a numerical costume.** Each is now
`std::max`'d or `#ifdef`'d against its historical value, so **the double path is
bit-for-bit unchanged** — checked at `MKL = OMP = 1`, where `soloviev-nstx`
reads `9.674113e-14` exactly as before.

**THE WORST WAS MFEM'S AND IS REPORTED**: `ExtensionBoundaryQuadrature` yields
**21 points per face in double and 0 in single**, because its rejection floor
follows `eps` while its `fd_step` does not. It returns normally having done
nothing, so a caller's rows stay zero and its "did we reach it" flag stays true.
In MEQ that dropped the exterior DtN coupling out of the bordered Jacobian while
leaving it in the residual → **[M-183](MEASUREMENTS.md#m-183)**.

**AND IT WAS RUN FIRST ON THE HOST** → **[M-174](MEASUREMENTS.md#m-174)**. `../mfem/install-single`
exists and MEQ solves against it. Three things it settled: `‖r₀‖` agrees with
the double build to every printed digit and the residual floor rises to
**2.2e-06** relative; **nothing converges out of the box because
`newtonRelativeTolerance` is `1.0e-12`**, seven orders below that floor, so even
a one-step Solov'ev case fails until it is rescaled; and the port is **45 sites
in 11 files**, byte-identical in double at `MKL = OMP = 1`. **The blocker is one
nobody had listed**: `fem/gslib.hpp` is raw `double`, so `MFEM_USE_GSLIB` cannot
coexist with `MFEM_PRECISION = single` — and every bordered free-boundary case
in the tree fails cold **in double**, so the interpolating warm start is
load-bearing physics rather than a convenience.

**AND THE BORDERED QUESTION IS ANSWERED: THE BORDER CLOSES TO `eps`, AND WHAT
DOES NOT SURVIVE IS ONE CONSTANT PER FILE** → **[M-204](MEASUREMENTS.md#m-204)**.
The gslib blocker went away when `meq::PointLocator` replaced it, so a single
install has the warm start. Against the **double build stopped at the same
tolerance** — which is the control, since otherwise the format's error and two
Newton steps' progress are added together — `limited-tokamak` closes its border
at `-1.240e-07` and `machine-f-diiid` at `-1.863e-07`, both at FP32's
`eps = 1.19e-07`, and the residual history follows double's over all ten
steps — the same order jump, the same non-monotone stretch — with the
disagreement widening from 2.0e-04 to 5.5e-02 as `‖r‖` itself falls by six
orders. **`NewtonRelativeTolerance = 1.0e-10` is what fails**: the
run stalls at `3.3e-07`, the reactive ladder retries, and it burns **400
iterations** before exiting 2. **AND THE ACCURACY IS THE CASE'S, NOT THE
FORMAT'S, OVER A RANGE OF 560×** — the four-border diverted machine is
**2.2 × eps** in `psi_ax` and the one-border limited case is **1,260 × eps**, so
a single-precision equilibrium is qualified per problem and not per build.
`DEVICE-PLAN.md` §7.3, §7.4 and §7.5.

**AND THE TEST TREE NOW BUILDS AND RUNS THERE, WHICH INVERTED WHAT THE PLAN
PREDICTED** → **[M-217](MEASUREMENTS.md#m-217)**. The port is **38 files and 405
lines** in four classes, and the double build is **62/63 — bit-for-bit
unchanged** on the same working tree, which is what the `std::max` construction
is for. Single reads **44 red of 63**, and the classification is the finding:
`DEVICE-PLAN.md` §7.2 called the tolerance sweep *"the single biggest cost"* and
named the flattened rate tail as the class to expect, and those are **25 and 8**
of 44 — while **23 of the original 45 never reached an assertion at all**,
because `newtonRelativeTolerance` is `1.0e-12`, seven orders below single's
residual floor, so `solve()` REFUSES. **A tolerance tuned against a solve that
did not converge is tuned against nothing**, so no assertion sweep could have
run first, and none would have found this — it is a solver constant, not an
assertion.

**THREE THINGS FROM IT GENERALISE.** *A default is only a default on the paths
that do not override it*: scaling the constructor's value alone moved 23 to 21
and took nothing green, because `setNewtonControl()` is called at **39 sites**
carrying their own literals — now one shared `meq::newtonTolerance()`, the
`limiterPolygonMaximum()` and `extensionFdStep()` precedent again. *A suite
count is not a measurement of a suite*: the repair moved the count by **one** and
the inventory by **eighteen**. And **`install-single` is not `install` minus
precision — it is minus six things**, two of them structural (SuiteSparse and
gslib have no single build) and three incidental, so a red test there cannot say
which it met; the attribution discipline of [M-210](MEASUREMENTS.md#m-210)
applies to a failure exactly as it does to a cure.

**AND A SINGLE BUILD TYPE-CHECKS A SEAM THE DOUBLE BUILD COLLAPSES.** `double`
and `real_t` are the same type there, so a function returning both compiles;
under `MFEM_PRECISION = single` it does not, and that is how MEQ's own item-4
wiring defect was found. **`std::max`'s refusal to deduce a common type is the
same instrument one level down** — it turned all 13 of the tree's
`std::max( 1.0e-300, ‖x‖ )` division guards into compile errors, which matters
because `1e-300` is below float's smallest SUBNORMAL, so as a float comparison
the guard does not loosen, it disappears.

**AND `install-single` IS NOW AT PARITY, WHICH ANSWERED THE RATE QUESTION IN
THE NEGATIVE** → **[M-218](MEASUREMENTS.md#m-218)**. CUDA, cuDSS and SUNDIALS
are restored; SuiteSparse and gslib stay structurally absent. It cost a
single-precision SUNDIALS — MFEM has **no guard** that `SUNDIALS_PRECISION`
matches `real_t` — and six `double`/`real_t` signatures in
`linalg/sundials.{hpp,cpp}` without which **MFEM does not compile in that
pairing at all**. And it turned up a **seventh** install difference that no list
of `MFEM_USE_*` flags contains: PARDISO resolved against `libmkl_sequential`
where the reference uses `libmkl_gnu_thread`, on the default trace solver.
*A configuration stated in flags is not a configuration.*

**`k+1` ABOVE `k = 1` IS A DOUBLE-BUILD ACCEPTANCE AND A SINGLE BUILD CANNOT
CARRY IT.** Against the double build on the same tree, `Solov'ev NSTX` in the
flux is 1.09× contaminated at `k = 1` on the finest mesh, **203×** at `k = 2`,
and **1.5× on the COARSEST mesh** at `k = 3` — so coarsening leaves the
asymptotic regime and refining goes deeper into round-off, and there is no
window. **Three per-pair criteria were tried and each provably fails**: no
constant separates a 9%-contaminated pair reporting a good 1.873 from a
36%-contaminated one; *"the error must at least halve"* is necessary but lets
through a pair that halves and still reads **3.03** against 3.85. The tail
**bends before it flattens**. So the skip is by DEGREE, loudly, and the
**absolute ceiling is still asserted** — widening `rateSlack` until the tables
passed would weaken every mesh in every build to accommodate meshes that are
not measuring anything, which is what *Testing stance* forbids.

**AND THE COUNT WAS THE WRONG INSTRUMENT FOR THE THIRD TIME.** Restoring the
install took the red count UP, 44 to 46, while every *"built without
`MFEM_USE_SUNDIALS`"* and *"no Device is configured"* disappeared — because a
test that aborted at its first case now runs all of them, `SourceViewTests`
going from 15 error lines to 436. Read the distinct failing CASES. The suite is
now **40 red of 63** at **157** cases, and the double control is **62/63 at
every step**.

**AND ELEVEN OF THOSE TESTS ARE NOW GREEN, WITH THE SAME LESSON EACH TIME** →
**[M-219](MEASUREMENTS.md#m-219)**. Single reads **29 red of 63** at 125 cases,
against M-218's 40 and 157. **Every one of them was asserting something that is
true in double and an ARTEFACT of double there** — an iteration count, a
bit-equality, a threshold with no scale in it, an order estimate needing fifteen
decades — and in every case the claim underneath survives the format. Asserting
it directly is often a STRONGER test: *"Newton takes one step"* becomes
*"one step reaches the round-off floor"*, which a merely-close Jacobian could
pass under a loose target and cannot reach in one step. The reduction reads
~500 `eps` in BOTH builds while the count reads 1 and 2.

**THREE DISTINCTIONS FROM IT GENERALISE.** *A difference is not a derivative*:
where two routes are both exact `eps` is the right floor, and where one is a
finite difference the band is `eps^(2/3)` to `sqrt( eps )` — measured at
1.8e-05 to 1.8e-04 on `HighBetaConvergence`'s corner route. *The robust half of
a rate claim survives where the `k+1` half does not*: `theXPointConverges...`
compares the finder against the FIELD it roots, which degrades on both sides and
holds at any precision, and only the `k+1` half is skipped. And *a controlled
experiment whose signal is below the floor is a double-build acceptance, not a
tolerance to widen*: `aPerturbedJacobianCosts...` measures an order gap of 0.6
that single reads as 0 on both arms, and loosening it would leave a case that
cannot tell a 5% Jacobian error from an exact one, looking green.

**AND TWO FAILURES IN `SourceViewTests` WERE NOT ROUND-OFF AT ALL.**
`meq::KnotView::psi` is `real_t` where `meq::Profile`'s abscissae are `double`,
so an interval BOUNDARY moves by an ulp and a `dF/dpsi` probed four ulps either
side reads 0 and 2.04e-03 — the two routes evaluate a rapidly varying function a
few ulps apart. The same thing at a coil edge switches **the whole coil term**,
-9.2153 against 0. Neither is a tolerance; both are the MFEM-free seam reaching
a boundary predicate.

**AND THE FLOORS ARE NOT ALL `eps`. THE FOUR THAT ARE NOT ARE THE NEXT
FINDING** → **[M-220](MEASUREMENTS.md#m-220)**. *A tolerance a solver was GIVEN
is the floor on what it can achieve*, and seven sites asked
`newtonResiduals().back() < 1e-8` while handing `solve()` a target
`meq::newtonTolerance()` floors at ten `eps` — M-217's lesson one level out,
with the override in the caller. *A fixed-point iteration cannot be asked for a
Newton's tolerance*: Newton's last step squares the residual and crosses a
target near the floor, Picard contracts by a fixed factor and never does —
Anderson–Picard runs its full 500 iterations and stalls at **5.80e-05** against
a target of 1.19e-06, so `meq::picardTolerance()` is `4096 eps` in single and
the caller's own value in double. *A stationary point's POSITION costs half the
digits and its VALUE costs none*: `sqrt( eps )` = 3.45e-04, measured on the
limiter contact flooring at 1.5e-04 and the located axis reading
`R = 0.99999583` — and the envelope argument that hides the position is the same
one that keeps the value at `eps`. And *an estimator differentiates, so its
floor RISES as the mesh falls*: `eta` reads 8.93e-05, 9.73e-05, 1.97e-04,
3.93e-04 over four refinements in single where double's finest is 1.927e-08, and
one derivative further a flux-surface average floors at `eps/h²`.

**AND THAT CROSSING IS WHY AN ADAPTIVE LOOP IN A NARROW FORMAT STOPS BEING A
DESCENT.** Past it, refining genuinely RAISES the error — the conditioning
growing like `h^-2` while the truncation falls like `h^(k+1)` — so the marker is
steering on round-off and neither `eta` nor the true error is monotone. It is
the estimator reporting what the arithmetic gives it, not the estimator failing.

**AND THE COMPENSATED-SUMMATION ROUTE INTO KINSOL DOES NOT EXIST, WHICH WAS
CHECKED RATHER THAN ASSUMED.** SUNDIALS 7.5.0's compensated sums are an
**ARKODE SPRK** runtime option — `ARKodeSetUseCompensatedSums()` — for a
symplectic stepper's time accumulation. There is no build flag, and
`src/kinsol/` contains the string nowhere. `picardTolerance()` is the answer
instead.

**AND SEVENTEEN TESTS WENT GREEN: SINGLE READS 12 RED OF 63 AT 61 CASES**,
against M-219's 29 and 125, with the double control at **62/63** — the one red
there being the standing deliberate `PlasmaEdgeConvergence`, so it is one of the
twelve rather than a finding. **What is left is a different kind of thing**: of
the eleven others, the remaining cases are dominated by solves that do not
CONVERGE rather than by assertions that cannot be read — `PlasmaConnectivity`,
`XPointOuter`, `DeviceParity`'s `mhd-rectangle`, `MachineFixedBoundary` (which
is additionally M-186's coin toss), and most of `FreeBoundaryCoupling`'s 13 and
`DriverAcceptance`'s 14. The tolerance sweep `DEVICE-PLAN.md` §7.2 called the
single biggest cost is substantially done; what stands between single precision
and a green suite is now the solver.

**AND THE SOLVER IS IN BETTER SHAPE THAN THE COUNT SUGGESTED** →
**[M-221](MEASUREMENTS.md#m-221)**. Measuring the EXAMPLES rather than the
suite: of the 22 that run without a mesh generator, **17 give the same verdict
in both builds** and one fails in both. Every shipped fixed-boundary machine
converges, and `psi_ax` agrees over a range of 1,600 — **6.9e-07** on
`fixed-f-diiid`, the largest problem in the tree, to 1.1e-03 on the
`q(psi)`-driven outer loop. That is M-204's *"qualified per problem and not per
build"* widened from two cases to eleven.

**THE SHARPEST ROW IS `diverted-tokamak-xpoint` SOLVING WHILE
`diverted-tokamak` DOES NOT.** Same machine, same mesh: the first puts the
X-point in the border and the second prescribes a limiter at it. So what fails
is not the diverted geometry, the free-boundary coupling, the moving support or
the coil subtraction — the working row exercises all four — but the BRANCH the
prescribed-limiter formulation lands on. At `MKL_NUM_THREADS` above 1 it
converges to `psi_bnd = -1.20e-01` against `+3.23e-02` and the driver's
axis-source guard refuses it: *the safety net works; the branch does not.*

**AND ONLY ONE OF THE THREE CAUSES WAS NUMERICAL.** `TraceSolver = "umfpack"`
in two shipped examples is a STRUCTURAL absence — UMFPACK is double only — and
it was six failing cases; no named trace solver exists that every supported
build has, which is why `meq::defaultTraceSolver()` is build-conditional and
why an example must name none. `locatePotentialPoint()` refused a point its own
TOML prescribes, because `TransformBack()`'s third answer is `Unknown` and
MFEM's `phys_rtol` is `4e-15` in double against `1e-6` in single; it now keeps
the best near miss and accepts it within `meq::samplerInsideTolerance()`.

**AND WHAT WAS LEFT AFTER THAT WAS NOT A TOLERANCE AND MOSTLY NOT MEQ'S: THE
TWO BUILDS ARE NOW AT PARITY, 62 OF 63 IN EACH, AND IT IS THE SAME TEST RED IN
BOTH** -> **[M-222](MEASUREMENTS.md#m-222)** -- `PlasmaEdgeConvergence`, the
standing deliberate one. Against M-221's 11 red and M-217's 44. The double
control held at 62/63 at every step, which is the whole point of the `std::max`
and `#ifdef` construction: a guard applied in the wrong place cannot weaken the
double build.

**THE LARGEST SINGLE ITEM WAS ONE LINE OF SOMEBODY ELSE'S HEADER.**
`ElementTransformation::TransformBack`'s default argument is `tol_0`, which is
`1e-15` in double and **`1e-7` in single** -- **0.84 of one `float` epsilon**,
so the physical residual it asks for is under an ulp. It is also TIGHTER than
the `InverseElementTransformation` it constructs, whose own `phys_rtol` is
`1e-6`. For a straight-sided triangle the Newton then reaches neither exit --
the reference one needs `eps|pt|/h < 4e-7`, which fails for every `h` below 0.31
-- runs to `max_iter`, and returns **`Unknown`**, the third answer, which a
caller testing `== Inside` reads as *not this element* for the element the point
is in. `meq::inverseTransformTolerance()` is the repair at MEQ's five
verdict-consuming sites, `1e-15` in double so every double call is unchanged.
Reported as `HDG-TRANSFORMBACK-SINGLE-TOLERANCE-FROM-MEQ.md`.

**IT WAS A WRONG ANSWER AND NOT ONLY A SLOW ONE, AND ONLY A COUNT COULD SEE
IT.** `ContourTracer`'s element walk fell back on the `O(elements)`
`Mesh::FindPoints` **108, 240 and 95** times over three surfaces where the
double build falls back none -- the trace still closing, because `FindPoints`
finds the point. And two band measurements were wrong rather than slow: the
deepest band excursion one fixture reports moved from **2.027 h to 0.628 h** on
the same mesh with only this tolerance changed, so the band's geometry was being
read off failed inversions. It also explains, and corrects, a result this file's
`PointLocator.cpp` already recorded: the warm start losing **458 of 25616**
nodes was cured by CONSTRUCTING the inverse rather than by the better initial
guess it was credited to. *A fix that works is not an explanation of why.*

**AND A SECOND MFEM DEFECT IS A DISPLACEMENT BELOW ONE ULP.**
`mfem::VertexConePath` builds a transfer direction by DIFFERENCING the mean of
its ray endpoints from the vertex, so a vertex lying ON `Gamma` -- which a
Cartesian grid puts at a circle's own extremum -- gives the zero vector and the
constructor aborts. `meq::markingOffset()` is MEQ's half: every
`MarkLevelSetSubdomain()` call in the tree passed `offset = 0`, so the set the
marking decides on is exactly the set round-off cannot decide. **Exactly zero in
double**, so `D_h` is the set it has always been. It closed a second symptom
nobody had connected to the first -- `DriverAcceptance` marking one shape twice
and reading **451 elements against 450**.

**THE REST IS M-186'S COIN TOSS ON FOUR MORE CASES, MEASURED WITH M-186'S OWN
INSTRUMENT.** `MKL_NUM_THREADS` is a numerically neutral perturbation applied to
one binary, and `tests/CMakeLists.txt` pins it at 1 -- one draw.
`pedestalSelfConverges` fails at 1 and converges at 2 and 4;
`fixed-h-circular` exits 2 at 1 and solves at 2, 4 and 8;
`PlasmaConnectivity`'s confined half fails at 1 and 2 and converges in 11 at 4,
to double's `psi_ax` to seven digits. **Every one is repaired by relaxing a
PRECONDITION and none by relaxing a tolerance**, and two of them are skipped
**by name** so that any other case failing still fails.

**AND THE MERIT-WEIGHT ROW TAUGHT SOMETHING NEW ABOUT THE PRECONDITION
ITSELF.** `XPointBorder` had skipped a weight whose solve *did not finish*. At
`MKL = 1` in single, weight 4 finishes -- in **195 bootstrap iterations against
the other three weights' 11** -- and comes to rest in a NEIGHBOURING BASIN,
reporting `converged`, at `psi_ax = 8.05e-02` against `8.27e-02`. At `MKL = 2`
the same weight runs past the cap and is skipped correctly. So *a run that
converges is not a run that converged to your fixed point*, and the
discriminator is the WORK: eighteen times the reference's iterations is not the
same descent.

**AND A `#ifdef MFEM_USE_SINGLE` IN A TRANSLATION UNIT THAT INCLUDES NO MFEM
HEADER IS SILENTLY FALSE.** `MachineFixedBoundary.cpp` drives the shipped binary
and needs nothing of the library, so it included none -- and the guard compiled,
ran, and took the other branch. `mfem/config/config.hpp` is what defines the
precision macros and is three lines of preprocessor. Same species as every other
instrument-not-answer finding here, and it cost one full suite run to notice.

**AND IT IS THE LAST ITEM OF THAT PLAN, NOT THE FIRST, WHICH IS THE ORDERING
DECISION WORTH KEEPING.** Everything on the device is validated by
`--device cuda` reproducing `--device cpu` **to every printed digit** — that is
how M-79 and M-130 were root-caused. A single-precision build destroys that
instrument: a disagreement is then either a wrong kernel or a narrower format
and nothing distinguishes them. **So residency is built in double, where it is
checkable, and converted afterwards.** The one piece that can run early is the
*experiment* rather than the sweep — a single-precision **host** build needs no
device and answers whether a bordered Newton closing `psi_ax` against
`-6.466e-12` survives a 1e-7 floor at all. `DEVICE-PLAN.md` §3.3 and §7.

**The transferable part**: exposing a knob makes reachable, by somebody who has
not read the library's documentation, every precondition that documentation
records. The Device requirement had been written down in `GradShafranov.hpp`
all along.

**THE INTERPOLATING WARM START IS WIRED, AND IT NO LONGER GOES THROUGH gslib.**
`meq::FieldTransfer` is **`meq::PointLocator`** in `WarmStart.{hpp,cpp}` —
MEQ's own point locator, about 620 lines, reproducing gslib to every printed
digit on both cases that exercise it and slightly faster
→ **[M-176](MEASUREMENTS.md#m-176)**. **`MFEM_USE_GSLIB` is now optional**,
which is the point: gslib has no floating-point abstraction at all, so it
cannot be built beside `MFEM_PRECISION = single`; its simplex path splits every
triangle into **three** quads; and MFEM's device findpts is gated on
`tensor_product_only`, so MEQ's triangles can never reach it
→ **[M-175](MEASUREMENTS.md#m-175)**. **The device half is MEQ's own and is one
of the SEVEN files in `src/meq` a device compiler owns** —
`PointLocatorDevice.cpp`, `SourceViewDevice.cpp`, `DeviceNPCDevice.cpp`,
`DeviceQuadratureDevice.cpp`, `DeviceRowsDevice.cpp`, `AxisSearchDevice.cpp` and
`ResidentIterateDevice.cpp` — because `mfem::forall`
`static_assert`s under `g++`; CMake marks that set and reads the architecture
out of MFEM's `config.mk` rather than guessing. Each compiles to a **host stub**
without CUDA, so the build is host-only and correct rather than short of a
route, and that is what CI builds.

`MFEM_USE_GSLIB` may still be `YES`,
`meq::FieldTransfer` was `FindPointsGSLIB` in `WarmStart.{hpp,cpp}`,
`WarmStartConvergence` is a registered ctest, and the driver reaches it in two
places:

* **A stored guess on a DIFFERENT mesh** transfers, so
  restarting from a run at another resolution — the ordinary way to use a stored
  answer — works. The exact restart is still taken when the meshes match, since
  it is every coefficient rather than an interpolation.
* **Adaptive cycles after the first start from the previous cycle**
  interpolated onto the refined mesh, rather than cold — which would throw away
  a converged answer at every refinement. The Dirichlet datum is the fallback
  at nodes the coarse mesh does not cover — and on the curved path there really
  are such nodes, since the computational domain grows as it refines: measured,
  48 and 360 of them on `miller-adaptive`.

**WHAT IT BUYS, AND THE MEASUREMENT NEEDED A NONLINEAR SOURCE TO EXIST AT ALL.**
MEQ's adaptive examples run a Solov'ev source, whose `∂F/∂ψ` is zero, so Newton
takes one step whatever it starts from — `miller-adaptive` shows `1` in every
cycle warm or cold, and the warm start is unmeasurable there by construction. On
`ManufacturedNonlinear` over three cycles at `k = 2`,
`carryingTheAnswerAcrossCyclesCutsTheWork` reads:

→ **[M-01](MEASUREMENTS.md#m-01)** — per cycle · total · final L2

**A third of the Newton work, and the same answer to 2.5e-12.** Both halves are
asserted: strictly fewer iterations, and an L2 that does not move — a warm start
must change the work and not the answer.

**`theDriverRunsTheAdaptiveLoop` now reads 4.4e-14 rather than 1.7e-16** for
exactly this reason, and that is a stronger assertion rather than a weaker one:
the driver warm-starts each cycle and the library loop it is pinned against
deliberately does not, so the number is now a statement that the warm start left
the equilibrium alone.

**The non-linear ordering is `NonlinearOrdering::NPC`**, since MFEM deleted the
mode. `CondenseThenLinearise` is kept as the backup and
is still the one that converges on stiff under-resolved meshes, at three to four
times the wall clock everywhere else. See `CLAUDE_HDGGS.md`, *The NPC port*.

**`docs/` is the user-facing account** and `README.md` points at
it; this file stays the maintainer's one. See *Layout* for why the numbers live
here and not there.

**`ROADMAP.md` is the priority order**, and says which items belong to MEQ and
which are requests on `../mfem-hdg-dev`, where another agent is working. This
file stays the technical record; that one is only about what to do first.

### The stages, and where we are

| | | |
|---|---|---|
| 0 | Git reconciliation, tag `v0-legacy` | **done** |
| 1 | Tree, CMake, `Config` / `Profiles` / `Source` | **done** |
| 2 | Linear `Δ*` on `DarcyForm`, fitted polygonal domain | **done** |
| 3 | Local post-processing `ψ*_h` | **dropped as MEQ code** — `DarcyForm::Reconstruct()` supplies it, at `k+2`. `CLAUDE_HDGGS.md`, *Post-processing is back* |
| 4 | Newton on the semi-linear source | **done** |
| 5 | Curved `Γ` by extension from subdomains | **done** |
| 6 | Adaptivity: the residual estimator and mesh update | **done** |

Beyond the port, four campaigns have their own plans, their own staging and
their own file:

| | | |
|---|---|---|
| toroidal flow | FL-0 to FL-8, **done** | `docs/rotation.rst` is the derivation, and `CLAUDE_FLOW.md` the record |
| solution inversion | IN-A to IN-P, **every stage done** | `INVERSION-PLAN.md`, and `CLAUDE_INVERSION.md` |
| free boundary | FB-A to FB-7 **done or answered**, FB-6 **met**, and §10's XP-0 to XP-4 **done** — the X-point is two unknowns of the bordered Newton, driveable from a file, and MEQ reproduces `freegs4e`'s DIVERTED equilibrium | `FREE-BOUNDARY-PLAN.md`, and `CLAUDE_FB.md` |
| inverse solve | IV-0 to IV-3 **done**: coil currents are border unknowns of the same Newton and isoflux targets its rows, in the library, M-250. IV-4, the X-point rows, is next; the driver refuses `Free` until IV-5 | `INVERSE-AND-RECONSTRUCTION-PLAN.md`, and `CLAUDE_FB.md` |

The two plan files that had nothing left in them are gone, converted to `docs/`:
`DRIVER-PLAN.md` was stage 7 and `FLOW-PLAN.md` was the flow campaign. Git has
them. What is left open across all four is short and `ROADMAP.md` is the
priority order for it.

### The campaigns are their own files

**`CLAUDE.md` GOT TO 7,600 LINES AND IS NOW AN INDEX WITH FOUR COMPANIONS.**
What stayed here is what a maintainer reads on every session — the status, the
build, the commands, the traps, the testing stance and the layout. What moved is
the core method and the three campaigns, each read when working on that one:

| | |
|---|---|
| **`CLAUDE_HDGGS.md`** | the equation, the two papers' sign errors, the discretisation and its conventions (`DarcyForm` holds `−q`, `τ` carries the opposite sign to eq (8e)), the post-processing, **the NPC port and the parity gap**, Newton and the obligation it creates, the bordered Newton on `ψ_ax`, the analytic fixture ladder, why MEQ's Newton struggles, `PicardThenNewton`, SUNDIALS, and the linear solves with the threading measurements |
| **`CLAUDE_FLOW.md`** | toroidal flow — `meq::RotatingSource` on RoPP (136), the two-species closure in closed form and `n` species by root find, the third derivative level it cost `meq::Profile`, the three errors in Li & Zhu, and Maschke & Perrin |
| **`CLAUDE_FB.md`** | free boundary end to end — FB-A to FB-7, the exterior DtN, the two borders, the plasma current, the moving support and the plasma edge's `k ≤ j` cap, the limiter as a point and as a curve, `ψ_ax`'s own defect, the half-disc mesher, and **the `freegs4e` benchmark** |
| **`CLAUDE_INVERSION.md`** | `ψ(R, z)` to `R(Ψ, l)` — critical points as roots of `q`, the contour tracer, the flux-surface averages, the Zernike disc basis, the gauge-free fit, the `(Ψ, θ)` output, and the `q(ψ)`-driven solve |

**The split is by campaign and not by size**, so a claim that touches two of them
is written in the file that owns the *subject* and cross-referenced from the
other. Three rules survive the split unchanged and apply in all five files:
**`MEASUREMENTS.md` holds the published tables under stable `M-nn` anchors and
they are not to be renumbered**; a hypothesis that was **measured and falsified**
stays, because recording it is what stops it being re-derived; and the history of
the *files* does not — present tense, no "used to", no dated "now fixed".

Each stage ends at a **measured convergence rate**, not at "it runs". See
*Testing stance* below for why that is the acceptance criterion.

## Commands

```sh
git submodule update --init --recursive     # extern/toml11
cmake -B build
cmake --build build -j6
cd build && OMP_NUM_THREADS=4 ctest -j4      # 52/53, about 480 s
```

**RUN IT `-j4` WITH `OMP_NUM_THREADS=4`, WHICH IS 3.2x FASTER AND MEASURED.**
Since `AssemblyMode::Threaded` became the default, each solver test already takes
every core through OpenMP, so a naive `ctest -j` oversubscribes. Holding the
product at the core count is what pays:

→ **[M-14](MEASUREMENTS.md#m-14)** — wall · CPU

37/37 in every configuration when that table was taken, and **62 of 63 today**,
the one red being `PlasmaEdgeConvergence` while the staircase `Γ_{p,h}` stands.
**`DeviceParity` NOW RUNS THREE ARMS, cpu AND debug AND cuda** →
**[M-228](MEASUREMENTS.md#m-228)**. The CUDA gap was left open on the grounds
that a third arm would be a second slow case, and that is backwards: `debug`
`mprotect`s the host page and is the expensive one, so the card arm costs
**35 s on a 427 s test**. It is not a repeat — `debug` catches a raw read of the
WRONG LOCATION and is blind to honest flags over a STALE COPY, which is what
cost `q-driven` a factor of two past a clean `debug` run and an unchanged Newton
count (M-206). Skipped loudly, on MFEM's own `Device configuration:` line, where
no device is configured. **Added before item 6 deliberately**: cuDSS as the
device default ends the digit-for-digit instrument, so this was the last moment
an exact CUDA arm could be pinned at all.

**`DeviceParity` is the newest case and the slowest after
`FreeBoundaryCoupling`** — 361 s, of which 329 s is `machine-f-diiid` under
`--device debug` — and it went green when upstream fixed
[M-214](MEASUREMENTS.md#m-214); it was red for exactly as long as that defect
stood, which is the stance working —
the count moves as cases are added, so read the table's ratios rather than its
absolute seconds. Nothing in the suite depends on a thread count, which is the
correctness half.

**WHERE A CASE LIKE `theTwoBorderSolveReportsATrueMagneticAxis` GOES RED, SUSPECT
THE FIXTURE BEFORE THE SOLVER.** `examples/limited-tokamak.toml` has **every** ingredient the case
has — a domain reaching `R = 0`, `[boundary.limiter]`, `[boundary.exterior]` and
`ConfineToPlasma` — and converges in **11 Newton steps** to a `Ψ` of 1.0000,
reproducing freegs4e to 1.3e-04. What it had that the fixture did not is **coils
and a prescribed current**. The fixture now has both, plus the confinement, and
it is **green at five limiter radii in a row**, 1.08 to 1.18, with the coil-free
control converging to the wrong topology at three of them; see `CLAUDE_FB.md`,
*§7.12b's fixture was the defect*, and `FREE-BOUNDARY-PLAN.md` §11.3 for the diagnosis and §11.7 for
the repair.

**AND "THE FIXTURE CANNOT BE REPAIRED IN PLACE" IS THE WRONG CONCLUSION TO DRAW
FROM IT BY ONE INGREDIENT.** `ConfineToPlasma` alone fails 4/4 and clamped
profiles alone fail 4/4; both measurements stand, and neither is a measurement
about the pair — the same run records "clamped WITH a prescribed current
converges in 22 steps". **A measurement that rules two things out is not a
measurement that rules everything out.**

**THE FOURTH CELL IS NOW RUN AND IT SAYS THE PAIR IS NOT ENOUGH EITHER** →
**[M-134](MEASUREMENTS.md#m-134)**. Confinement **with** a prescribed current,
swept over 250x in the target, reaches **2 of 4 at best** — 0 of 4 below
`μ₀I_p = 0.20`, 2 of 4 at 0.50 and again at 2.00, and back to 0 of 4 at 5.00, so
the maximum is interior and the sweep is not truncated. The clamped arm is cured
outright by that same lever. So a prescribed current moves confinement off zero
and does not repair it, and `theTwoBordersConvergeTogether` keeps
`FREE-BOUNDARY-PLAN.md` §11.2's retraction of its `ψ_ax` and `ψ_bnd` as physics.
**The lesson survives its own example**: it took the full 2 x 2 to say this, and
each of the four one-key runs on its own supported a different wrong reading.

**IT MUST NOT ASSERT `checkAxis().agrees`, WHICH IS A TAUTOLOGY UNDER THE
LOCATED-AXIS CONSTRAINT**: `Ψ` there reads 1 by construction — measured, that
very configuration reports `Ψ = 1.0000 AGREES` on a solve whose axis source is a
pole. A test that cannot fail is worse than no test. It asserts four things, and
`checkAxisSource().bounded` is the one with teeth.

**AND IT HAS NOW HAPPENED TWICE, ON A CASE THAT HAD BEEN RED FOR SESSIONS.**
`XPointOuter` was failing at its bootstrap solve with a long, careful list of
ruled-out causes and two surviving suspects, both of them about the solver. It
was neither: the fixture called `setPlasmaSupport()` on the plasma source
**before** constructing the coil wrapper, so the wrapper's own flag stayed false,
`plasmaComponentWanted()` read false, and **XP-1's flood fill never ran on the
only diverted case in the tree** — 333 elements of private flux region, 44% of
the candidates, carrying a current channel nobody asked for. `apps/meq.cpp` has
always done it the other way round and is correct. See `CLAUDE_FB.md`,
*`setPlasmaSupport()` only reaches the fill if it is called on the wrapper*, and
→ **[M-82](MEASUREMENTS.md#m-82)**.

**THE PART WORTH GENERALISING IS HOW IT HID.** Every experiment on that case
changed one key and read the outcome, which is the right instinct and is exactly
what could not find this: the fill was off in **both** arms of every one of them,
so a real measurement — "turning `ConfineToPlasma` off makes it converge" —
supported a conclusion that was wrong. It took a **2 × 2 cross** to see that the
fill and the frozen plasma edge are each necessary and neither sufficient. *A
one-key experiment separates two hypotheses only if everything else is where you
think it is.*

Per *Testing stance*, a red suite is the intended signal while a defect stands.
**Do not "fix" it by relaxing an assertion** — the repair in both cases above was
to give a fixture something it was missing, and the assertions never moved.

**AND A THIRD CASE WENT RED WHERE HEAD WAS GREEN, AND IT WAS NEITHER THE FIXTURE
NOR THE SOLVER — IT WAS A COIN TOSS, AND `HEAD` LOSES IT TOO.**
`XPointBorder`'s `theXPointMeritWeightChangesTheWorkAndNotTheAnswer` failed at
`XPointMeritWeight = 4` with *"weight 4 did not solve"* while a `HEAD` worktree
passed it in 81 Newton iterations against a cap of 200. **A worktree comparison
says *this is not how the tree was*; it does not say *this is a regression***,
and the step that separates them is a **numerically neutral perturbation applied
to BOTH binaries**. `MKL_NUM_THREADS` is that perturbation — it reassociates a
blocked BLAS-3 and changes nothing else. Bootstrap iterations at weight 4,
against the file's cap of 200:

| | MKL=1 | MKL=2 | MKL=4 | MKL=8 | MKL=16 |
|---|---|---|---|---|---|
| `HEAD` | 81 | 65 | 115 | **over the cap** | **over the cap** |
| working tree | **over the cap** | 84 | 75 | 103 | 103 |

**`HEAD` FAILS THE SAME CASE IN TWO OF FIVE CONFIGURATIONS**, and every arm that
converges reports `psi_ax` to eleven digits and the X-point to six decimals —
the invariant the case exists for, which never moved. `MKL_NUM_THREADS=1` — what
`tests/CMakeLists.txt` pins, and the only value ctest ever runs — is one draw,
not a property of a tree. → **[M-186](MEASUREMENTS.md#m-186)**.

**The repair is the same shape as `PedestalConvergence`'s and the assertions did
not move.** The case took `BOOST_TEST_REQUIRE( got.sweeps > 0 )` per weight,
which is a **precondition** on a chaotic plateau — the residual sits near
`1.8e-04` for sixty iterations, RISES to `4.3e-04`, then falls to `1e-12` in
five — and in the more brittle direction, since a `REQUIRE` abandons the weights
after it as well. It now sweeps four weights either side of one, prints and
skips a weight that does not solve inside the cap, and requires three of the four
to compare. **What was relaxed is a precondition, not a tolerance**, and that
distinction is the whole of what makes this a repair rather than the thing this
section forbids. **That the majority is a requirement and not a second coin toss
had to be measured before it could be asserted**: the other three weights read
`13, 13, 14` iterations at every `MKL_NUM_THREADS` tried, none of them failing,
so weight 4 is the only one near the plateau and at most one dropout is
reachable.

**AND THE ROUND-OFF HAD NO SEMANTIC CAUSE TO FIND, WHICH IS WHY BISECTING WOULD
HAVE BURNED THE SESSION.** Diffing the two binaries' full transcripts, 192 lines
of 377 differ and **every one is a residual at `1e-12` or below** differing in
its fifth figure — `1.5e-17` absolute against a problem scale of `7e-2`. No
iteration count, no `psi_ax`, no support count differs anywhere else in the
binary. The tolerances that had been rewritten as
`std::max( <the double value>, <something following real_t's epsilon> )` each
**are** their historical value in double, checked one at a time. What moves the
last bits is the compilation: a literal moved behind a function call or a
`std::max` changes which FMA contractions `g++` emits at `-O3`. *Identity of a
constant is not identity of the emitted code, and on a plateau that is enough.*

**`OMP_NUM_THREADS` is deliberately NOT pinned in
`tests/CMakeLists.txt`** — a plain `ctest` must still exercise threaded assembly
as it ships, and pinning it would quietly change what the bit-exactness cases
test.

**AND THE FLOOR IS `naming`, NOT A SOLVER.** It is the slowest test in the suite
at **199.7 s** — one `clang-tidy` invocation over every file in `src/meq`,
single-threaded — against `PedestalConvergence`'s 178 s. So 223 s is very nearly
"the suite costs one lint run", and going below it means parallelising
clang-tidy (`run-clang-tidy`) rather than anything about the solver.

**AND `naming` IS NO LONGER AT THE TOP AT ALL.** At **51** tests the two longest
are **`FreeBoundaryCoupling`** and **`naming`**, and back-to-back runs of the
same tree read them at 281.8 s / 262.6 s and 311.7 s / 294.0 s — within a factor
of 1.1 of each other in both runs, so `-j4` cannot overlap them once everything
else has finished and the run is very nearly the cost of the longest chain
rather than of the total work. `DriverAcceptance` read 77.3 s, 94.4 s and 80.3 s;
`XPointBorder`, XP-3's acceptance, 35.4 s, 49.7 s and 34.0 s for five solves of
the diverted machine. **Three readings of one tree, 477 s, 545 s and 443 s**,
which is the paragraph below making its own point rather than anything having
changed. **The third of those carries XP-4's driver run and is the FASTEST**:
`DriverAcceptance` grew a whole diverted free-boundary solve, 11.5 s standalone,
and read 14 s LOWER than its previous figure. Read the ratios.

**AND THE `FreeBoundaryCoupling` FIGURE HAS BEEN 276 s, 433.8 s AND 381.7 s
WITHOUT THE CODE BEING THE VARIABLE ANY OF THOSE TIMES.** 276 s is stale — the
case grew FB-7's fourth acceptance and the limiter sweep since. Of the other two,
433.8 s and 381.7 s are the SAME content: 440.5 s standalone under load, 433.8 s
in a contended run, 381.7 s under `-j4` on a quiet machine. `LimiterCurve` does
the same thing at a third of the size — **20.2 s** in the quiet run against 29.5
and 31.8 s standalone while two agents were building. **Three readings of one
binary, spread by 1.6x**, which is the whole of why the paragraph below says to
read ratios.

**THOSE SECONDS ARE APPROXIMATE ON PURPOSE AND THE RUNS THAT PRODUCED
THEM WERE BOTH UNDER CONTENTION** — one shared the machine with a stray
busy-wait loop of the session's own making, the next with another agent's build.
**A suite time is only a measurement on an idle machine**, and this file already
records a 490 s / 540 s spread on identical code. Read the ratios; re-time on an
idle machine before reading anything into a change of tens of percent.

**AND IT HAS JUST GROWN BY 14 MINUTES, WHICH IS WORTH KNOWING BEFORE SOMEBODY
TIMES THE SUITE AND CONCLUDES SOMETHING.** FB-T's
`confinementWithAPrescribedCurrentIsTheFourthCell` is **forty solves** — a 2 x 2
re-run whole with the target current swept across eight values — and it is kept
at its full sweep deliberately, because a trimmed one would no longer reproduce
[M-134](MEASUREMENTS.md#m-134)'s published table. Measured on two `-j4` runs of the
same tree: `FreeBoundaryCoupling` **1149.60 s and 991.99 s**, against 288.83 s
before FB-T and against `naming` at 364.44 s and 346.69 s. **So the suite is now
that one binary**, not the longest chain of several — the whole run reads
exactly `FreeBoundaryCoupling` and nothing else, and its own 1.16x spread
between two runs of identical code is this machine rather than anything having
changed, and a suite time taken after this is not comparable with one taken
before it.

**`FreeBoundaryCoupling` overtakes the lint and it is buying
something.** It carries FB-7's four acceptances, a coupled filament run, and
since later that day a sweep of the limiter constraint against its own control;
its two-border case now solves a **physical** equilibrium eight times at `k = 2`
on 1333 elements — five limiter radii with the vertical field, three without —
where the unphysical one it replaced solved a cheaper problem badly. It did not
get slower on the day it grew, which is the interesting part: the trim that keeps
the coil-free control off the two radii where it cannot converge saves about 300
Newton steps and pays for the new sweep. `PlasmaConnectivity`'s expensive halves are the
diverted-fixture fill at three resolutions and the jump study's 6401-sample
sweep, both of which are the measurements the case exists for.

**Read those seconds as one machine's one run.** `CLAUDE_HDGGS.md`'s *PARDISO
and the MKL link line* records a 490 s / 540 s spread on identical code — so a change of a few
tens of percent in a suite time is this machine, not a regression.

**ctest needs no environment set by hand.** `tests/CMakeLists.txt` puts
`MKL_NUM_THREADS=1` on every registered test, without which the suite takes well
over an hour instead of nine. It is the only variable MEQ needs; for the one
`MKL_THREADING_LAYER` — see *Traps* for why it is not one. Running a test
binary **directly** needs the same one thing:

```sh
MKL_NUM_THREADS=1 ./tests/SolovievConvergence
```

The performance harness is separate and is not a ctest:

```sh
tests/performance/scan.sh build 3            # both thread axes, ~40 min
tests/performance/npc-scan.sh build 3        # the NONLINEAR one
```

**`npc-scan.sh` sweeps the two thread axes TOGETHER, which is the opposite of
what `scan.sh` does, and the inversion is the point.** `scan.sh` holds one axis
at 1 while sweeping the other because sweeping both was measured to be
misleading — and that rule is correct for *serial* assembly. Under threaded
assembly the element-local work is nested inside an OpenMP region, where MKL
suppresses its own threading, so `OMP = MKL = N` is the configuration that
matters. See `CLAUDE_HDGGS.md`, *Threading, measured*.

`MFEM_DIR` defaults to `../mfem/install` and also reads the environment, so the
`-D` is usually unnecessary. Never a bare `make -j` or `cmake --build -j`; see
*Traps*.

### Coverage, and what CI can and cannot check

```sh
cmake -B build-cov -DMEQ_ENABLE_COVERAGE=ON && cmake --build build-cov -j6
cd build-cov && ctest
gcovr --root .. --filter 'src/meq/' --print-summary        # or --html-details
```

`MEQ_ENABLE_COVERAGE` is an option, not a build type, on purpose: a coverage
build wants `-O0` and a release build wants `-O3`, and that difference should be
chosen rather than inherited from whatever `CMAKE_BUILD_TYPE` happens to hold.

**AND THE TEMPTING REASON IS WRONG IN A WAY THAT COSTS SOMEBODY AN AFTERNOON.**
Routing `-O0` through `CMAKE_BUILD_TYPE` does **not** drop `NDEBUG` in a way that
changes which `MFEM_ASSERT`s are live. Audited:

* `MFEM_ASSERT` is gated on **`MFEM_DEBUG`**, never on `NDEBUG`
  (`general/error.hpp`), and nothing anywhere in the installed tree derives one
  from the other.
* The installed MFEM has **`MFEM_DEBUG = NO`**, so **every `MFEM_ASSERT` is
  already dead in every MEQ build**, release or debug.
* MEQ's own sources contain **no `assert()` and no `<cassert>`**, so `-DNDEBUG`
  disables nothing of MEQ's either. It is inert here.

So `cmake -DCMAKE_BUILD_TYPE=Debug` on MEQ buys the debugger and nothing else,
and the standing advice under *Traps* to take a debug build "for this reason
alone" **did not work as written** and has been corrected. Recovering those
assertions needs a *second MFEM install* configured `MFEM_DEBUG=YES` — a
different and much larger undertaking, and one nobody has done.

**CI CANNOT BUILD THE SOLVER, AND THIS IS STRUCTURAL.** The MFEM MEQ needs is
`meq-integration`, a local merge published on no remote, so a hosted runner
cannot obtain it and caching does not help — there is nothing to fetch.
`find_package(MFEM)` is therefore **not `REQUIRED`**: without it `CMakeLists.txt`
builds the MFEM-free half — `Config`, `Profiles`, `Source`, `SourceFactory` — and
`tests/CMakeLists.txt`'s existing "unit not present" guard skips every
convergence test without needing to know why. So `ci.yml` runs four unit suites,
the `naming` check and a coverage gate at 90% lines. **Treat a green CI badge as
evidence about the configuration and profile layers and about nothing else** —
the `full` job is written but `if: false` until the branch is published or a
self-hosted runner exists.

**NO COVERAGE PERCENTAGE IS QUOTED HERE, DELIBERATELY** — a stale one is worse
than none, and this number moves with every solve path added. Two things about
the *measurement* are stable, and they are the transferable part.

**A line-coverage percentage on this codebase is partly a measurement of comment
density.** On the same run, `GradShafranov.cpp` read **92.7%** to `gcov`
(executable lines only, 586 of them) and **64%** to `gcovr` (its own line set,
838). This project comments heavily, and `gcovr` counts roughly 250 lines that
`gcov` does not consider executable at all — then reports most of them as
uncovered. The CI gate uses `gcovr`, so that is the operative number for the
gate; `gcov`'s is the operative number for "is this code tested".

**A claim this file carried, and which was wrong, is worth leaving recorded**:
that the solver sat at 60% against 92–100% for the configuration layer, so "the
best-covered code in MEQ parses TOML and the least-covered is the part every
physical claim rests on". The 60% was `gcovr`'s number compared against
`gcov`-flavoured intuitions. On one denominator the solver and `Config.cpp` were
the same. **A measurement about the tool, not about the code** — the same trap
as the threaded-MKL one.

What is true regardless of the denominator is that **22 of
`GradShafranovSolver`'s 51 methods were called by nothing**, and
`SolverContract.cpp` was written for that rather than for a percentage — it
found `setSource()` silently accepting a second source of the same kind, and
value faults throwing `logic_error` where the constructor throws
`invalid_argument`. **Branch coverage, last seen at 48%, is the figure with real
headroom in it and the one nobody has looked at.**

**THE DOCS BUILD, AND `.venv-docs` IS HOW. USE IT.**

```sh
.venv-docs/bin/sphinx-build -W -b html docs docs/_build/html
```

**Verified: `build succeeded`, 31 pages, zero warnings under `-W`**,
which is what `.readthedocs.yaml`'s `fail_on_warning` requires. The venv is
gitignored (`.gitignore:110`) and carries exactly what `docs/requirements.txt`
pins — **sphinx 7.4.7** and **sphinx-material 0.0.36** — so a local pass and an
RTD pass mean the same thing, which is the whole point of the pin.

**DO NOT REACH FOR `sphinx-build` ON `PATH`, AND DO NOT CREATE A SECOND VENV.**
The system one is **8.2.3**, and `docs/requirements.txt` caps sphinx below 8 on
purpose: sphinx-material has been unmaintained since 2023 and predates
Sphinx 8. `python3 -c "import sphinx_material"` also fails against the system
interpreter, because the theme is only in the venv — and `docs/conf.py:75`
imports it unconditionally, so the failure is an unhandled `ModuleNotFoundError`
at config-eval time rather than anything naming a theme.

**AND THIS FILE SAID "THE DOCS STILL DO NOT BUILD HERE" FOR ONE SESSION, WHICH
IS WHY THE VENV IS NAMED HERE RATHER THAN LEFT TO BE FOUND.** The check that
produced that claim ran `command -v sphinx-build` and `python3 -c "import ..."`
— both of which answer for `PATH` and the system interpreter and neither of
which can see a venv nobody activated. A commit went out carrying an
UNVERIFIED docs edit on the strength of it, and the recipe offered was to build
a venv that already existed. **An absent tool and an unactivated venv are the
same two commands and different facts**; `ls -d .venv*` separates them and costs
nothing.

For the record, since it is genuinely useful when a theme is what is missing:
`sphinx-build -W -b dummy` parses every page and enforces `-W` **without needing
a theme or writing output**. It is the right tool for checking prose and
cross-references in an environment that cannot render. It is not needed here.

**`gcovr` is not installed on this machine**, so the recipe above fails with
`command not found`; a venv is the way round it, and it needs
`--merge-mode-functions=separate` or it aborts on `Config.hpp`'s inline
accessors appearing at two line numbers across translation units. Run it from
the repo root: `--root .` after a `cd` into the build directory silently matches
nothing and reports 0%.

### Why toml11 is a submodule and MFEM is not

**toml11 is `extern/toml11`**, pinned to a commit, as MaNTA does it. It is header
only, small, and its version is something MEQ should control rather than inherit
from whatever is installed — the `find_or<double>` behaviour recorded under
*Traps* is version-dependent, and a silent config bug is exactly what a floating
dependency buys you. Currently `b32a2ff`, v4.4.0-31, the same commit MaNTA pins.

**MFEM stays out of tree**, hand-built, found through `MFEM_DIR`. Its history is
enormous, so a submodule would make every clone of MEQ expensive; it needs an
out-of-tree configure-and-build of its own anyway; and `../mfem-hdg-dev` is a
working tree with its own active development on `gf-hdg-subdomains-dev`, which a
submodule pin would fight rather than help.

If configure fails with toml11 not found, the cause is almost always a clone
without `--recursive`; the error message says so.

### Prefer a maintained library to a hand-rolled algorithm

**Standing preference, and it is about maintenance rather
than speed.** Where a well-known algorithm has a well-maintained
implementation — a Householder QR, a rank-revealing least squares, a special
function, a spline — **take the library**, and take it even when a measurement
says the hand-rolled version is no slower. *"The chance that you found the
perfect Householder implementation and are able to maintain it indefinitely
into the future is low."*

A performance tie is therefore **not** an argument for keeping bespoke code. The
bar for bringing in a header-only, well-known, well-maintained library is that
it does not break something else — not that it wins a benchmark.

**AND FOR DEVICE CODE THE PREFERENCE IS ORDERED AND THE MIDDLE RUNG IS THE ONE
WORTH SPELLING OUT: `mfem` ABSTRACTION > CUB > A MANUAL KERNEL.** MFEM's own
facilities come first — `mfem::forall` for a launch, `mfem::reduce` with
`SumReducer` and friends (`general/reducers.hpp`) for a reduction,
`mfem::GeometricFactors` and `mfem::DofToQuad` for element data,
`mfem::BatchedLinAlg` for dense blocks. **CUB comes before a hand-written kernel
and after MFEM** for one reason that is not about quality: MEQ builds for CUDA
today and `MFEM_USE_HIP` is a real target, so a CUB-only route is a route that
has to be written twice, while MFEM's abstractions dispatch to both. A manual
kernel is the last resort and needs a reason in a comment.

**The one thing that outranks the order is a contract MEQ has already made.**
`meq::CompressedRows` documents bit-identity with a dense loop's every partial
sum, and `assemblePlasmaCurrent()` refuses `reduction( + : )` so that the answer
does not depend on the thread count — so `mfem::reduce` is the right primitive
for a single scalar and the *wrong* one where a reassociation would cost a
digit somebody has published. `meq::DeviceRows` is the worked example:
`mfem::forall` over rows with a sequential ascending sum inside, because
`mfem::reduce` is single-output and would reassociate
→ **[M-190](MEASUREMENTS.md#m-190)**.

What this does **not** license: a dependency that is heavy, obscure, unmaintained
or hard to obtain, and it does not override the reasons recorded above for why
MFEM is out of tree and toml11 is pinned. And it is not a reason to rewrite
working code on sight — it decides which way to go when the question is live.

Already taken on these grounds: **Boost.Math** for the Zernike radial
polynomials, where the Jacobi route is both more accurate than the textbook
factorial sum and somebody else's to maintain — see `CLAUDE_INVERSION.md`, *The disc basis*. **Eigen**
(3.4.0, `/usr/include/eigen3`, header only) is the natural home for the dense
linear algebra in `src/meq/SurfaceFit.cpp`, which currently hand-rolls a
Householder QR and a one-sided Jacobi SVD in about 150 lines.

**The one thing to hold fixed when swapping in a library**: `SurfaceFit`'s
truncation threshold and its Levenberg–Marquardt damping together *are* the
gauge that makes the gauge-free fit well posed, and a different rank-revealing
decomposition can truncate different directions. The ellipse result, the `nstx`
tail, the exactly-zero axis spread and the positive minimum Jacobian are the
numbers that say the swap was clean. If one moves, that is a finding, not a new
baseline.

## The equation, the discretisation and the solve → `CLAUDE_HDGGS.md`

The Grad-Shafranov equation as MEQ poses it, the two papers' disagreement about
the sign of the Solov'ev source, eq (20)'s `η₅` erratum, the LDG-H spaces, and
the two conventions that cost time if rediscovered — **`DarcyForm` holds `−q`**
and **the papers' `τ` carries the wrong sign** — are all in `CLAUDE_HDGGS.md`,
together with the post-processing, the NPC port, the Newton chapter and the
linear solves.

### Which MFEM, and why not master

**MEQ builds against `../mfem/install`** — MFEM **4.10.1** on branch
**`meq-integration`**, CMake-built in `../mfem/build` from sources in
`../mfem/mfem-src`. The options actually set, read out of
`share/mfem/config.mk` and re-checked 2026-09-01 rather than remembered:

| | |
|---|---|
| `MFEM_USE_SUNDIALS` | KINSol, so `KINSolver(KIN_LINESEARCH)` is reachable. **Now `../sundials/cuda-install`, not `../sundials/install`** — see the CUDA row |
| `MFEM_USE_GSLIB` | `FindPointsGSLIB`; gslib v1.0.9 built alongside at `../mfem/gslib` |
| `MFEM_USE_SUITESPARSE` | UMFPACK, the direct solver MEQ's own solver runs on. **`../suitesparse/install`, not Debian's** — v7.12.2, the same version Debian ships, built against oneAPI MKL. Debian's carried a `NEEDED` on `libblas.so.3` -> `libmkl_rt.so`; see `CLAUDE_HDGGS.md`, *PARDISO and the MKL link line* |
| `MFEM_USE_MKL_PARDISO` | `PardisoSolver`, oneAPI MKL 2026.1, **on the threaded `mkl_gnu_thread` layer**. Resolving it to Debian's `mkl_sequential` instead is a bigger deal than it sounds — see *Traps* |
| `MFEM_USE_EXCEPTIONS` | `MFEM_ERROR_THROW` by default, which is what makes the driver's exit code 2 reachable |
| `MFEM_USE_OPENMP` + `MFEM_THREAD_SAFE` | **both, and both load bearing.** `DarcyHybridization::SetAssemblyMode( Threaded )` needs both and **aborts** rather than falling back without them. MEQ reaches it through `setAssemblyMode()` |
| `MFEM_USE_LAPACK`, `MFEM_USE_ZLIB`, `MFEM_USE_DOUBLE`, `MFEM_USE_MEMALLOC` | |
| `MFEM_USE_MPI = NO` | MEQ is serial throughout; `../mfem-hdg-dev-par` exists for parallel work |
| **`MFEM_USE_CUDA = YES`**, `CUDA_ARCH = sm_75` | CUDA 13.3. **`sm_75` because the card is compute capability 7.5**; MFEM's cache defaulted to `sm_60`, which is wrong here and nothing warns you |
| **`MFEM_USE_CUDSS = YES`** | cuDSS 0.8.0 from Debian (`/usr/include/cudss.h`), found with `-DCUDSS_DIR=/usr`. `mfem::CuDSSSolver` agrees with UMFPACK to 3.5e-14; its **timings on this machine are not reproducible** and MEQ's solver does not use it. See `CLAUDE_HDGGS.md`, *Threading, measured* |
| `RAJA`, `OCCA`, `CEED`, `SIMD`, `AMGX` all `NO` | |

**Turning CUDA on is not one flag, and two of the three obstacles are silent.**

* **SUNDIALS must itself be a CUDA build.** `linalg/sundials.hpp` hard-errors
  with `MFEM_USE_CUDA=TRUE requires SUNDIALS to be built with CUDA support`,
  testing `SUNDIALS_NVECTOR_CUDA`. MEQ cannot drop SUNDIALS — `KINSolver` backs
  `AndersonPicard` and `PicardThenNewton` — so the answer is a second SUNDIALS.
  `../sundials/cuda-install` is 7.5.0 configured `-DENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=75`, otherwise matching the old one.
* **`MFEM_ENABLE_TESTING` is ON, and it needed an upstream fix to be.** A
  serial CUDA build would not generate at all: the bug was that a device test loop merged the serial and parallel example lists and
  never re-split them, registering a test against `sundials_ex9p` in a build
  that never created it. Both loops — `examples/CMakeLists.txt` and
  `examples/sundials/CMakeLists.txt` — now guard the parallel branch with
  `elseif (MFEM_USE_MPI)`, so a serial build skips it. **Verified**:
  reconfigure, full rebuild and install all returned 0, and ctest registers
  **131 tests**. `MFEM_ENABLE_TESTING` does not appear in the installed
  `config.mk`, so it changes nothing MEQ's finder reads.

  **This file claimed the bug was "filed as
  `../mfem-hdg-dev/doc/CMAKE-SERIAL-CUDA-SUNDIALS.md`". IT WAS NOT.** That path
  has never existed in that tree's history and exists nowhere on disk — the
  report was never written. Treat a "filed as" claim in this file as needing a
  check, not as evidence.
* **Stale derived cache variables survive a `-DSUNDIALS_DIR=` change.**
  `SUNDIALS_INCLUDE_DIR` and friends are separate cache entries, so repointing
  the directory leaves the old include path on `TPL_INCLUDE_DIRS` *in front of
  or behind* the new one, and which wins is an ordering accident. Delete every
  `SUNDIALS_*` entry but `SUNDIALS_DIR` and reconfigure.
* **AND THE SAME CACHE SURVIVES A NEW SUNDIALS *COMPONENT*, WHICH IS WORSE
  BECAUSE THE SYMPTOM NAMES NOTHING.** `meq-integration` was
  re-created with a `sundials-ida-integration` branch merged in, which adds
  `mfem::IDASolver` and adds `IDAS` to `SUNDIALS_COMPONENTS` in MFEM's
  `CMakeLists.txt`. `SUNDIALS_LIBRARIES` is a cached STRING, so
  `find_package(SUNDIALS)` short-circuited on the old list and
  `-lsundials_idas` never reached the link line — while `libmfem.a` now carries
  `sundials.cpp.o` referencing `IDASVtolerances`, `IDAPrintAllStats` and eight
  more. **Every MEQ test binary failed to link**, with an error naming IDA and
  nothing naming a cache. MFEM's source was correct throughout. Same cure:
  delete every `SUNDIALS_*` entry but `SUNDIALS_DIR`, reconfigure, reinstall.
  **Not filed upstream**, per the standing rule about findings against work that
  has just landed — and there was nothing to file, since the defect was entirely
  in MEQ's own build directory.

**And a CUDA build broke MEQ's own `FindMFEM.cmake`, which is worth knowing
because every serial build hides it.** MFEM puts the toolkit headers on the line
as `-isystem <dir>` inside `MFEM_CXXFLAGS`, not as `-I<dir>` in
`MFEM_TPLFLAGS`, while `general/backends.hpp` includes `<cuda_runtime.h>`
unconditionally once `MFEM_USE_CUDA` is defined. MEQ's finder took only
`TPLFLAGS`, so the first `#include <mfem.hpp>` failed on a missing CUDA header —
a confusing place to meet a build-configuration problem. The finder now also
parses `-isystem` out of `MFEM_CXXFLAGS`.

**WHAT CUDA BUYS IS THREE THINGS AND `fem/darcy/`'s MISSING PARTIAL ASSEMBLY IS
NOT ONE OF THEM.** Partial assembly is still absent and still irrelevant — its
economic case is sum factorisation and MEQ is triangles. What the flag enables
is `BatchedLinAlg`'s `gpu_blas` backend and `CuDSSSolver`; **MFEM's batched
element and face kernels**, which are reached through `AssemblyMode`,
`LocalFactorMode` and `TraceAssemblyMode` and **are** the device route; and
**MEQ's own device translation units** — `DeviceNPCDevice.cpp`,
`DeviceQuadratureDevice.cpp`, `DeviceRowsDevice.cpp`, `PointLocatorDevice.cpp`,
`SourceViewDevice.cpp`, `AxisSearchDevice.cpp` and `ResidentIterateDevice.cpp`,
which CMake hands to `nvcc` at the architecture it reads out of `config.mk`, and
which hold every `mfem::forall` launch in `src/meq`.
Without CUDA those seven compile to host stubs and the build is host-only and
correct, which is what CI builds. → **[M-223](MEASUREMENTS.md#m-223)**.

**AND IT IS NOT FREE — IT COSTS 7% OF THE WALL CLOCK ON THE CPU PATH, MEASURED
AGAINST AN OTHERWISE IDENTICAL INSTALL.** `../mfem/install-nocuda` differs from
`../mfem/install` in `MFEM_USE_CUDA` and the `MFEM_USE_CUDSS` that must follow
it, and in nothing else; the two produce **byte-identical output** on the machine
case, so this is cost and not answers.

→ **[M-77](MEASUREMENTS.md#m-77)** — wall · user · allocation count

The mechanism is visible and only partly the obvious one: `mfem::forall` on the
host path builds an `__nv_hdl_wrapper_t`, nvcc's host-lambda wrapper, whose
`manager::do_call` shows in `perf` at 1.9–2.9% — and interposing `operator new`
at its call site counts **226.7 million allocations in a 30-second run**, of
which **55.4 million (24%) go away with CUDA off**. The other 171 million are
MFEM's ordinary per-element temporaries and are not CUDA's doing. So the flag is
a prerequisite that currently charges rent, and a CPU-only production build has
a measured reason to be a separate install.

**`../mfem-hdg-dev` is now the development tree and not what MEQ links.** That
split exists because the alternative was demonstrated: mid-session that tree was
switched to `gf-hdg-dev`, which has no `fem/darcy/extension_hdg.{hpp,cpp}`, and
MEQ's `naming` check failed on a tree MEQ had not touched. A library MEQ depends
on should not move under it while somebody is working on it. Point `MFEM_DIR` at
the dev tree deliberately when testing a fix; otherwise leave it alone.

### `meq-integration`: the branch MEQ builds from, and why it is local only

MEQ needs work from **fourteen branches**, thirteen of them `hdgdev`'s and one off
MFEM master, and for as long as this section existed it named four. See the warning under the containment loop for how the
other five were found and why no check built from this table could find them:

| branch | what MEQ needs from it | |
|---|---|---|
| `gf-hdg-subdomains-dev` | `fem/darcy/extension_hdg.*` — stage 5's curved `Γ`, and `TransferredDatumCoefficient`. **AND IT IS THE TARGET FOR ANYTHING MEQ SENDS BACK ABOUT THE EXTENSION MACHINERY**, free boundary's transmission quadrature included: it carries `TransferPath`, `ElementExtension`, `PathTraceCoefficient` and `ExtensionRegionQuadrature`, and carries **no NPC at all** — verified 2026-09-05 — so a patch against it provably does not depend on the ordering | merge base |
| `direct-solver-symbolic-reuse` | `UMFPackSolver` and `PardisoSolver` keeping their symbolic factorisation across Newton steps | merge 1 |
| `gf-hdg-linearise-first` | **`DarcyNPCOperator` / `DarcyNPCSolver`** — the NPC method, which is what MEQ's default ordering now is — plus `SetAssemblyMode`, `SetGradientMode`, `SetLocalFactorMode`. **`SetNonlinearOrdering` was on this list and is DELETED**; see `CLAUDE_HDGGS.md`, *The NPC port* | merge 2 |
| `gf-hdg-dev` | ***"The postprocessing closes on the element average, always"*** — the reconstruction fix, without which `ψ*` is a different function wherever `∂F/∂ψ` vanishes | **merge 3, again** |
| **`pardiso-multi-rhs`** | **`PardisoSolver::ArrayMult`** — the multi-RHS trace solve. MEQ's bordered step depends on it for **1.34×**, → **[M-98](MEASUREMENTS.md#m-98)**, and **losing it does not break the build**: `Operator::ArrayMult`'s base loops `Mult()`, so a recipe that drops this branch compiles, answers identically and is silently slower | merge 4 |
| **`pardiso-device-host-sync`** | ***"PardisoSolver: synchronise the vectors, not just the matrix"*** — the other half of M-79's device chain, on the default trace solver | merge 5 |
| **`sundials-ida-integration`** | `mfem::IDASolver`, and the `IDAS` component that broke every MEQ test binary's link when the cache went stale — see the CUDA section's fourth bullet | merge 6 |
| **`gf-hdg-p-adaptivity`** | upstream's own §2 work; MEQ calls nothing from it, and it is in this table because it is in the merge, which is the point of the table | merge 7 |
| **`memory-manager-registration-refcount`** | **NOT an `hdgdev` branch — it is off MFEM master**, which is why a recipe built from the `hdgdev` list cannot see it. `internal::Memory` gains the reference count `internal::Alias` always had, so registering one address twice no longer lets the first holder's `Delete` evict the second's entry. Four routines register lazily and reach it: `MakeAlias` on an unregistered base, the device-class arms of `Read`/`Write`/`ReadWrite`, and `SetDeviceMemoryType`. MEQ proposed this shape after both projects' first mechanism proved wrong; see the `Reconstruct()` row above | merge 9 |
| **`cmake-find-package-component-cache`** | the fix to the defect MEQ reported when `SUNDIALS_LIBRARIES` short-circuited on a stale cached list and `-lsundials_idas` never reached the link line — see the CUDA section's fourth bullet. **It is in the merge because MEQ asked for it**, and its absence is a link failure rather than a silent one | merge 10 |
| **`cudss-solver-configuration`** | `CuDSSSolver`'s setters, including the one that refuses `CUDSS_MTYPE_SYMMETRIC` by name rather than returning an exactly zero solution every phase reports as success. MEQ does not call cuDSS from a config file and still wants this: `TraceSolverScaling` constructs one | merge 11 |
| **`intrules-triangle-high-order`** | `c6693d5f09`, the simplex-rule ceiling. Triangle rules are tabulated to **order 25** and above that `GrundmannMollerSimplexRule()` returns partly NEGATIVE weights — correct in exact arithmetic, destructive in double, and returned with a one-time warning rather than refused. **MEQ is clear and it was checked rather than assumed**: `SourceIntegrator::rule()` asks for `2k + OrderW + 4`, which is 12 at the suite's `k = 4`, and `PlasmaEdgeConvergence`'s triangle rules are `2k + 14` over `k <= 3`, so 20. `setSourceQuadratureOrder()`'s doxygen already stated the bound | merge 12 |
| **`gf-interp-hdg-dev`** | **CCSZ INTERPOLATORY HDG** — `fem/darcy/reaction_hdg.{hpp,cpp}` with `HDGInterpolatoryReactionIntegrator`, `NodalReactionFunction` and a quadrature control, plus `HDGPostprocessBlocks` in `postprocess_hdg.*` and `DarcyHybridization::Bg_data`. **15,024 insertions across `fem/darcy/`**, verified in `libmfem.a` by `nm`. MEQ does not call any of it yet; `INTERPOLATORY-HDG-PLAN.md` is the plan | merge 8 |
| **`sundials-real-t`** | **`linalg/sundials.{hpp,cpp}` in `real_t` rather than `double`**, plus a `static_assert` that `real_t` and `sunrealtype` agree and a makefile that stops refusing SUNDIALS in a single-precision build. **MEQ reported it and MEQ needs it**: without it MFEM does not compile at all in the `MFEM_PRECISION = single` + `MFEM_USE_SUNDIALS` pairing, which is `install-single`, and KINSOL is what backs `AndersonPicard` and `PicardThenNewton`. MEQ's own six sites were the ones that stop the compile; upstream found **73**. It is a branch off `origin/master` rather than one of the HDG line, so it is the second one the `hdgdev/`-anchored pattern below cannot see | merge 13 |

**THE TOPOLOGY HAS NOW CHANGED THREE TIMES AND THIS ROW HAS BEEN WRONG TWICE.**
It has read "pairwise independent", then "an ANCESTOR of both, so it needs no
merge of its own" — verified 2026-08-30 and true then — and as of **2026-09-05
`gf-hdg-dev` is once more ahead of the merge and needs a merge of its own**, by
a single commit (*"Four Darcy test files the trunk's CMake build never ran"*),
which conflicts in `tests/unit/CMakeLists.txt` alone. **The standing instruction
is therefore to run the containment loop and believe it over this table**, which
is a claim about someone else's branch and goes stale without MEQ doing
anything. The reconstruction fix itself has been contained throughout; a "not
contained" here has always meant *the branch has moved*, never *the fix is
missing*.

What has not changed is *why* `gf-hdg-dev` is in the list. Dropping it gives a
silently degraded `ψ*` on any element where the Jacobian's reaction term
vanishes, an adaptive loop a full order worse, and a suite that stays green
until `thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` runs. It is
free today rather than unnecessary.

`meq-integration` is their merge, and it exists **only** to be built against. It
is not a development branch and must not be pushed anywhere: it will be
re-created whenever any of them moves, so anything committed directly to it is
lost. Do the work on the branch it belongs to and re-merge.

**It lives in `../mfem/mfem-src`, which is MEQ's own clone — never in
`../mfem-hdg-dev`.** See the rule below.

**Expect the fetch to be a forced update.** Two of the four were rebased between
2026-08-29 and 08-30, so a merge *into* the existing `meq-integration` is not
the operation you want; re-create it. Tag the old one first — the tree carries
`meq-integration-before-2026-08-30` for exactly that reason — because a bad
re-merge is otherwise unrecoverable without redoing the conflict resolution.

Refreshing it, from `../mfem/mfem-src`:

```sh
git fetch hdgdev gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
          gf-hdg-linearise-first gf-hdg-dev
git tag -f meq-integration-before-$(date +%F) meq-integration
git checkout -B meq-integration hdgdev/gf-hdg-subdomains-dev
git merge hdgdev/direct-solver-symbolic-reuse     # 3 conflicts; see below
git merge hdgdev/gf-hdg-linearise-first           # 3 conflicts; see below
git merge hdgdev/gf-hdg-dev                       # 1 conflict as of 2026-09-05
git merge hdgdev/pardiso-multi-rhs                # M-98 depends on this
git merge hdgdev/pardiso-device-host-sync
git merge hdgdev/sundials-ida-integration
git merge hdgdev/gf-hdg-p-adaptivity
git merge hdgdev/gf-interp-hdg-dev                # the interpolatory half
git merge hdgdev/sundials-real-t                 # single precision needs it
git merge memory-manager-registration-refcount   # NOT hdgdev's -- off master
# ALWAYS finish with this, and believe it over the table above:
for b in gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
         gf-hdg-linearise-first gf-hdg-dev pardiso-multi-rhs \
         pardiso-device-host-sync sundials-ida-integration \
         gf-hdg-p-adaptivity gf-interp-hdg-dev sundials-real-t \
         cmake-find-package-component-cache cudss-solver-configuration \
         intrules-triangle-high-order; do
  git merge-base --is-ancestor hdgdev/$b HEAD && echo "$b contained" \
                                              || echo "$b NOT CONTAINED"
done
git merge-base --is-ancestor memory-manager-registration-refcount HEAD \
  && echo "memory-manager contained" || echo "memory-manager NOT CONTAINED"
```

**Run that last loop.** It is nine cheap commands and it is the only thing that
catches a branch having moved out from under the recipe — which is precisely
what happened to the three-branch version of this section.

**AND IT CANNOT CATCH A BRANCH THE RECIPE NEVER NAMED, WHICH IS THE WORSE
FAILURE AND IS THE ONE THAT HAPPENED.** The loop asks *is each branch I listed
contained*, so it is only ever as complete as the list — and the list was
missing **five of nine**. `gf-interp-hdg-dev` alone is 15,024 lines of
`fem/darcy/`, the whole interpolatory-HDG half; `pardiso-multi-rhs` is what
M-98's 1.34× rests on. Re-creating `meq-integration` from the recipe as it stood
would have deleted all five, and **four of them without a compile error**:
`Operator::ArrayMult` has a working base implementation, and MEQ calls nothing
in `reaction_hdg.hpp` yet. The build stays green, the answers do not move, and
the solver gets slower.

**Ask the tree the other question instead.** It needs no list, which is the
whole point, and it is what found this — **and it has now found a tenth branch
the same way**, `memory-manager-registration-refcount`, which is off MFEM master
rather than `hdgdev` and which the `hdgdev/`-anchored pattern below would itself
have missed. Match both:

```sh
git -C ../mfem/mfem-src log --merges --format=%s meq-integration \
  | grep -oE "hdgdev/[a-z0-9-]+|memory-manager-[a-z-]+" | sort -u
```

**The pattern is as much a list as the loop was**, which is the lesson repeating
one level down: an `hdgdev/`-anchored grep encodes the assumption that every
merge comes from that remote, and the first branch that does not is invisible to
it. Widen the pattern or read the merge subjects unfiltered.

**BUT THE LOOP CHECKS THE BRANCH, NOT THE INSTALL, AND THAT GAP HAS ALREADY
BITTEN.** Every command in it runs in `../mfem/mfem-src` and says nothing
whatever about `../mfem/install`, which is what MEQ actually links. On
2026-09-01 the loop reported **all four contained** while the installed library
was `fa65a2f932` and `meq-integration` had moved to `d244329d8c` — **514 changed
lines in `fem/darcy/` across seven files**, including *"The HDG face quadrature
never saw the trace element"* and *"DarcyOperator dereferenced a null
prolongation on a hanging-node-free NC mesh"*, the second in exactly the
non-conforming meshes the adaptive loop makes. So the one check this project
trusts passed on a tree whose measurements were all taken against a library a
day out of date.

**A REFRESH IS A CLEAN REBUILD OF BOTH TREES, AND THE REASON IS THE CLASS
LAYOUT.** Anything that moves a member behind `#ifndef MFEM_THREAD_SAFE` — which
MEQ's build defines, and MEQ constructs such objects itself — changes an object's
size without changing a signature, so an incremental link is an ODR violation
with no diagnostic. **MEQ must be fully rebuilt too**, `cmake --build build -j6
--clean-first`, and the absolute timing tables are historical afterwards:
[M-138](MEASUREMENTS.md#m-138)'s budget is the one that wants re-taking, since
it is absolute where [M-137](MEASUREMENTS.md#m-137)'s payoff is a ratio between
two MEQ binaries on one library and survives.

**BOTH INSTALLS ARE CURRENT.** Both checks below pass on `../mfem/install` and
on `../mfem/install-single` against `meq-integration` at `4929772e95`, and the
DIII-D machine reads `psi_ax = 3.759850e-01` and a constraint of `-6.466e-12`
under `--device cuda` on the double one, as it did on the install before --
[M-243](MEASUREMENTS.md#m-243), where the double suite is 65 of 66. **The
single suite read 60 of 66 at `969569440e` and no answer moved**, and it has
not been re-run since -> [M-242](MEASUREMENTS.md#m-242): two device-route tests
assert bit-identity that only a double build can give, and three solves sit on
M-222's knife edge at `MKL_NUM_THREADS = 1`. The single install is a
**second** tree with the same obligation and its own reason to move first:
`sundials-real-t` is in the merge because MFEM does not compile **at all** in
the `MFEM_PRECISION = single` + `MFEM_USE_SUNDIALS` pairing without it.

**AND REFRESHING IT IS NOT FREE: THE LAST REFRESH TOOK A CASE RED BY SEVEN
ORDERS, AND THE ANSWER IT MEASURES HAD NOT MOVED AT ALL**
→ **[M-189](MEASUREMENTS.md#m-189)**. `../mfem/install` rebuilt from
`meq-integration` at `6d05189467` takes `FreeBoundaryCoupling` red at
`theTransmissionRowIsTheBoundaryIntegralItClaims`, which has read 7.8e-16 to
1.4e-13 since it was written and now reads **1.27e-06**. Bisected to one default
argument: `e38906d853` — upstream's fix to MEQ's own M-183 — changes
`ExtensionBoundaryQuadrature`'s `fd_step` default from the literal `1e-6` to
`ExtensionFDStep()`, which is `cbrt( eps )` and therefore **6.06e-06 in
double**. **MEQ's four production sweeps pass their own step and were never on
that default**; what took it was the reference the TEST builds, so the case
contracted a row at `1e-6` against a reference at `6.06e-6`. Give both halves
the same step and it returns to 1.4e-16.

**THE LESSON IS ABOUT THE TEST AND NOT ABOUT THE LIBRARY.** A central
difference's error is a deterministic function of its step, so it is common to
two sweeps that share one and **cancels exactly** — which is the only way that
case could ever have read 7.8e-16, orders below what a difference at `1e-6` can
resolve. **Where a caller computes one quantity two ways and compares them, both
sides must be given the same step**, and `meq::extensionFdStep()` is now
`inline` in `GradShafranov.hpp` so they can share one expression rather than two
constants kept in step by hand — the `meq::limiterPolygonMaximum()` precedent,
M-165. And MEQ's first report of this blamed the production Jacobian: it was one
`grep` from being checked and was not checked.

**Check the install too, and it is one command:**

```sh
diff <( git -C ../mfem/mfem-src show meq-integration:fem/darcy/darcyhybridization.hpp ) \
     ../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp \
  && echo "install is current"
```

**AND THAT CHECK IS BLIND TO A `.cpp`-ONLY FIX, WHICH IS MOST OF THEM.** It
diffs an installed HEADER, and only headers are installed — so a library that
is a bug fix behind reports *current*. Measured: `meq-integration` moved from
`93952626e2` to `dfb5525a35` carrying `d720f02d62`, the fix to MEQ's own
M-214, and every installed header in `fem/darcy/` was byte-identical across it.
The check passed on a library missing the fix MEQ had asked for.

**Ask the TIMESTAMPS instead, which cannot miss a translation unit:**

```sh
test ../mfem/install/lib/libmfem.a -nt \
     <( git -C ../mfem/mfem-src log -1 --format=%cd --date=iso meq-integration >/dev/null; \
        echo ) 2>/dev/null
# or, plainly:
[ "$( date -r ../mfem/install/lib/libmfem.a +%s )" -ge \
  "$( git -C ../mfem/mfem-src log -1 --format=%ct meq-integration )" ] \
  && echo "install is newer than the branch tip" || echo "INSTALL IS STALE"
```

The header diff still earns its place — it catches a *layout* change, which is
the case that silently produces an ODR violation rather than a missing fix — so
run both. They fail in different directions.

A branch that another agent re-creates is a branch that moves without MEQ doing
anything, so this belongs beside the containment loop rather than in a rebuild
checklist nobody reads.

**And it has already fired again. Run 2026-09-01 against `meq-integration` at
`fa65a2f932`: three of the four report contained and `gf-hdg-linearise-first`
does not**, because that branch has advanced **17 commits** since the merge.
**Read a "not contained" as *the branch has moved*, not as *the merge failed*** —
what MEQ needs from it is in the tree (`DarcyNPCOperator` is in the installed
`fem/darcy/darcyhybridization.hpp`, and NPC is what MEQ runs on). It says a
re-merge is available, not that the build is wrong. Re-create rather than merge
into the existing branch, per the recipe above.

**The second merge conflicts, and WHICH THREE FILES IT IS HAS ALREADY CHANGED
ONCE.** This file recorded `gf-hdg-linearise-first` as merging clean, then as
bringing `HDG-DEFECTS-FROM-MEQ.md`, `HDG-ROADMAP.md` and the makefile. As of
**2026-09-05 the defects report no longer conflicts** — it is gone from both
sides — and `tests/unit/CMakeLists.txt` has taken its place:

| file | resolution |
|---|---|
| `doc/HDG-ROADMAP.md` | **Not a disjoint pair.** Four regions, two of them large *divergent rewrites of the same prose*, so keeping both sides duplicates whole sections and produces a worse file than either. Take `--ours`, the subdomains side: it is the merge base, it is the newer of the two, and each dev branch keeps its own copy of that roadmap anyway. MEQ neither builds nor reads it |
| `miniapps/hdg/makefile` | a union: `extension` from one side, `navierstokes` from the other. Merge the `SEQ_MINIAPPS` and `PAR_MINIAPPS` lines into one, keep the `ifeq (MFEM_USE_SUNDIALS,NO)` filter, and keep both test targets |
| `tests/unit/CMakeLists.txt` | a union of a sorted source list — `test_darcy_extension.cpp` from subdomains, `test_darcy_npc.cpp` from linearise-first. **`test_darcy_system.cpp` appears INSIDE one hunk and again in the common region below it**, so a mechanical keep-both duplicates it and CMake then compiles it twice. Take the union and assert no entry repeats |

Only the makefile and the CMakeLists affect a build. **The `HDG-ROADMAP.md`
marker count is not scriptable**; count them each time, for the `CHANGELOG`
reason below.

**AND THE SAME UNION TRAP APPEARS AGAIN IN MERGE 3.** `gf-hdg-dev` conflicts in
`tests/unit/CMakeLists.txt` alone and by then the merged file already carries
the union, so `--ours` is right there — but verify it rather than assuming, by
checking every darcy entry from all four parents is present exactly once.

**The first merge conflicts and the resolution is always the same.** The
symbolic-reuse branch is based on a much later upstream (767 commits past where
the subdomains branch left), so `CHANGELOG`, `doc/CodeDocumentation.dox` and
`fem/nonlininteg.cpp` collide where both sides appended. In every case the two
sides are **disjoint** — `nonlininteg.cpp` gets `SumNLFIntegrator` and
`SumBlockNLFIntegrator` from one side and the Navier–Stokes PA kernels
`ConvectiveVectorConvectionNLFIntegrator` and
`SkewSymmetricVectorConvectionNLFIntegrator` from the other — so **keep both
sides** of all three. Check the function names are still disjoint before assuming
that; it is true today and is a property of the branches rather than a rule.

**And "keep both sides" is not quite enough for `nonlininteg.cpp`.** Git leaves
the file's trailing

```
}

}
```

as *common* content after the `>>>>>>>` marker — one brace closing
`SumBlockNLFIntegrator::~SumBlockNLFIntegrator`, one closing `namespace mfem`.
That is correct for whichever single side you keep, and one brace short when you
keep both: the first side's destructor is then left open and the second side's
functions sit at namespace depth 2. **Add the destructor's closing brace at the
seam.** It fails to compile rather than failing silently, but the error appears
forty lines later on the *other* side's first function, which is a misleading
place to start looking. `git show <branch>:fem/nonlininteg.cpp | tail` on both
parents shows what the ending should be, and counting braces across the three
versions finds it in seconds — both parents balance and a bad merge does not.

*A trap when resolving `CHANGELOG` programmatically*: it contains long `======`
heading underlines, so a script that strips lines *beginning* with `=======`
destroys the file. Match the conflict marker exactly, alone on its line.

Then rebuild:

```sh
cd /home/ian/projects/mfem && cmake --build build -j4 && cmake --install build
```

### `../mfem-hdg-dev` is receive-only

**Write MFEM development *requests* into `../mfem-hdg-dev/doc/` and nothing
else.** No branches, no commits, no checkouts, no builds. Another agent works in
that tree and it is theirs; MEQ consumes it through `git fetch` from
`../mfem/mfem-src` and in no other way.

This is written down because it was got wrong: `meq-integration` was first
created *in* the dev tree, given two merge commits and a fix commit, and that
tree was left checked out on it — so the next person to touch it would have been
on someone else's branch. Fetching from the dev tree is reading; everything else
is not.

**WHAT IS OPEN WITH UPSTREAM, AND NOTHING ELSE ABOUT IT.**

| | |
|---|---|
| **Auxiliary globally-coupled unknowns** | `SetNumAuxiliaryUnknowns()` and the two per-element assemble hooks. **The one unbuilt piece, and no longer worth asking for**: the entry points MEQ named are public already as `NPCReduce()` and `NPCRecover()`, and `DarcyNPCSolver::ArrayMult` applies them to several right-hand sides in one pass, so a differenced border is `K` applications of a routine that blocks them. **MEQ CALLS IT** — the bordered step queues every column and flushes once, worth **1.34×** on the DIII-D solve leg at 14 columns, → **[M-98](MEASUREMENTS.md#m-98)** |
| **`SetIntegratorsThreadSafe()`'s SCOPE, which is the durable half and is now upstream's doxygen** | **CLOSED.** MEQ reported that `mfem::NormalTraceJumpIntegrator` carries `Vector face_shape, normal, shape1_n, shape2_n` and two `DenseMatrix` unguarded (`fem/bilininteg.hpp:4872`), that MEQ installs exactly one as the object handed to `EnableHybridization()`, and that `DarcyForm::ReconstructFluxAndPot()` takes it per face inside its element loop. **Guarded upstream, on the trunk**, so all five branches take it — and a clean MFEM rebuild is owed here, since the class layout moves in a `MFEM_THREAD_SAFE` build. **MEQ's report overstated one thing and the correction is the tense**: it said the class is reached by "the loop you are threading", and there is no such loop — `darcyform.cpp`'s three `#pragma omp parallel` regions are all assembly and both reconstruction entry points are serial. MEQ inferred that from being told the work was in progress and wrote the inference as an observation. What upstream took as the real defect was **theirs and not the class**: `SetIntegratorsThreadSafe()` made a promise whose scope was written down nowhere, so threading one more loop silently widens what a caller who set it months earlier has undertaken. Its doxygen now says the promise is about **every integrator installed on the hybridization** rather than about whichever loops happen to be threaded |
| ~~**`mfem::HDGExtensionIntegrator` carries unguarded scratch**~~ | **CLOSED, and MEQ's own files went on asserting it afterwards.** Upstream's `26fccc9f3d` put `shape, shape_ext, nor, x, xbar, m, y, CTm`, `Cmat`, `L` and `ext` behind `#ifndef MFEM_THREAD_SAFE` (`extension_hdg.hpp:24`), and `git merge-base --is-ancestor` confirms it is in the `meq-integration` MEQ links. The row that stood here, and the matching comment at MEQ's `SetIntegratorsThreadSafe()` call, both described the window between the eight-class commit and that one — **a true conclusion resting on a premise that had expired**. Neither was re-read against the tree. *A claim about somebody else's code is only as fresh as the last time it was checked* |
| ~~**`Reconstruct()` aborts under `Device("debug")`**~~ — **CLOSED, and neither project's mechanism was right** | Fixed at `86b0b0fd13`, **contained and installed**, and `DeviceParity` is green on it. The abort is `CheckHostMemoryType_` at `mem_manager.cpp:1778`, from `Vector b_ze( b_z.GetData() + e*nd_ut, nd_ut )` with `e == 0` aliasing `b_z`'s base pointer exactly. **MEQ's mechanism was wrong and so was upstream's first one, in the same way.** Both said `Wrap( ptr, n, own = false )` takes `GetHostMemoryType()` so `std_delete` is false and `Delete_` evicts — but `Delete_` opens `if (!mm.exists \|\| !registered) { return; }`, so an unregistered wrap is **inert** and `h_mt` decides nothing. What actually registers the view is **`Memory<T>::MakeAlias()`**, whose "Register 'base'" branch fires whenever `IsDeviceMemory( GetDeviceMemoryType() )` — true under `debug` and `cuda` — setting `Registered \| OWNS_INTERNAL` on the **view**; `MemoryManager::Insert()`'s `emplace` then keeps the existing entry **silently**, and the view's destructor takes the `Known` branch and erases the **owner's**. In `ReconstructTotalFlux` the trigger is `b_zi.MakeRef( b_ze, ... )` and needs no device flag. **Both of MEQ's candidate fixes fail**, measured: keying on `flags & Registered` is a no-op because `registered` IS true by then, and making `Wrap` take `HOST` when `!own` aborts on a different assertion. What landed is MEQ's **third** option, the one it did not recommend — build the sub-views with `MakeRef` into the owner rather than wrapping a raw `GetData()` — which also syncs, where a raw `GetData()` is a host read of a buffer whose live copy may be on the device. **Scope: 2 of 5 branches**, `gf-hdg-linearise-first` and `gf-interp-hdg-dev`, both of which `meq-integration` merges. **"Why `cuda` never met this" has now been wrong three times and this is the version with evidence**: not that the destructor is inert (`h_mt` decides nothing), and not that the eviction "erases a key that is not there" (`Erase` does `if (mem_map_iter == maps->memories.end()) { mfem_error( "Unknown pointer!" ); }`, so a missing key **aborts**). Under `cuda`, `GetHostMemoryType() == HOST` leaves the **owner** unregistered, so the view's own lazy registration is what **creates** the entry and erasing it on destruction is **balanced**. The immunity is a property of the host memory type, and it is the same observation that makes a reference count the right shape: one registration one erase is already correct, and it is two registrations of one address where today's first erase takes the second's entry. **A second question is live and is upstream's user's**: whether to close the trap in `Memory` rather than per caller. MEQ said yes, in `../mfem-hdg-dev/doc/HDG-MEMORY-OWNERSHIP-FROM-MEQ.md`, and its addendum argues the change is **smaller than either project first thought** — `maps->aliases` is already reference-counted against exactly this hazard (`InsertAlias` does `alias.counter++` on a collision, `EraseAlias` does `if (--alias.counter) return;`) and `maps->memories` is not (`Insert` discards the `emplace` result, and the `auto res =` is itself inside `#ifdef MFEM_DEBUG`; `Erase` erases unconditionally). The lazy registrations the counter must be correct over are **four routes across five call sites** — `MakeAlias`'s "Register 'base'", `SetDeviceMemoryType()`, and `Read`/`Write`/`ReadWrite`; the two `Wrap` sites need a caller ownership error first. So a counter on `internal::Memory` closes it **without touching ownership semantics or any call site**, where the proposed `own_internal` flag cannot: the call site cannot know whether the address was already registered, so no constant value of it is right. Upstream verified the shape on their trunk, confirmed `mem_manager.cpp` is unmodified from MFEM master in that family, and **withdrew the flag**. **The counter's own cost is bounded rather than hypothetical**: an unbalanced register/erase already leaks a `memories` entry today — `Memory<T>::Reset()` sets `flags = 0` without going near `Delete_`, and `Vector::StealData` is `*p = data; data.Reset(); size = 0;` on a `Memory` that may be registered, with `linalg/hypre.cpp` and `linalg/petsc.cpp` doing the same — so a counter **extends** that class rather than creating it, and turns an unbalanced pair from *an entry erased out from under its owner* into *an entry that stays*. Strictly milder and still silent |

**AND THE BATCHED LOCAL RESIDUAL THAT TAKES A CALLER'S INTEGRAND IS BUILT ON
BOTH SIDES AND CLOSED** → **[M-197](MEASUREMENTS.md#m-197)**,
**[M-223](MEASUREMENTS.md#m-223)**. `machine-f-diiid` under `--device cuda`
reports `batched nonlinear residual yes, gradient yes, device integrand yes`,
so **all eight `CanBatch*` predicates** read yes on MEQ's production bordered
configuration, the whole NPC traversal included. **A default device run prints
`face assembly NO` and that is the KEY and not the predicate**: the driver
leaves `AssemblyMode` Serial under a Device to keep bit-reproducibility
([M-215](MEASUREMENTS.md#m-215)), and asking for `AssemblyMode = "batched"` on
the same file reads `face assembly yes` — with the constraint moving from
`-6.465e-12` to `-6.466e-12`, which is exactly the last-bit difference that
default is protecting. *A predicate printed under a mode nobody asked for
answers about the mode.* What it was asked for
is below, because the SHAPE of the seam is the transferable part.
`HDG-DEVICE-RESIDENCY-FROM-MEQ.md`, request 3.
`CanBatchLinearResidual()` refuses MEQ on two clauses — `m_nlfi_p` not being
linear, and `lop_type == LocalOpType::PotNL` — and **both refusals are
correct**, so MEQ asks for a *sibling* rather than a loosening: a batched
residual for the `PotNL` case in which the element-local nonlinearity is
supplied by the caller. The shape proposed is a **header template over a
callable**, instantiated in MEQ's own device translation unit, so that upstream
owns the element loop and the block layout while MEQ owns only the integrand
and neither header names the other's types — the same seam
`meq::PointLocator` already uses one level down. MEQ's half is built and tested
→ **[M-177](MEASUREMENTS.md#m-177)**. **It is the last route MEQ's
configuration needs**: upstream's own audit puts MEQ at YES on the other seven
`CanBatch*` predicates.

**AND IT IS NO LONGER BLOCKED ON ANYTHING OF UPSTREAM'S EXCEPT THE BUILD ITSELF.**
Upstream's reply is explicit that this does not wait on their request 1 — *"it is
not blocked on request 1; request 1 is done. It is blocked on an interface, and
the interface is ours"* — and ends by asking MEQ for the one thing only MEQ can
supply: *"Send us the POD and the signature you want to write `F` and `F'` in and
we will build the other side."* **Sent**, as
`../mfem-hdg-dev/doc/HDG-BATCHED-RESIDUAL-SEAM-FROM-MEQ.md`: `meq::SourceView`
and its two POD companions, the `MFEM_HOST_DEVICE` evaluators, the
element-local `sourceResidualAtPoint()` / `sourceGradientAtPoint()` pair, and the
four constraints that are not decoration — `radius` unpre-divided because `F/R`
is `μ₀ j_φ` and MEQ's meshes reach `R = 0` exactly, the `psi_c` shift, the
per-ELEMENT connectivity flag, and the view being a **per-call** argument because
`psi_ax`, `psi_bnd` and the current scale move between residual evaluations
inside one Newton step.

**AND ONE HAS CLOSED, WHICH IS THE ROW THAT USED TO SIT HERE.** Threading
`NPCReduce()` and `NPCRecover()` — MEQ's largest measured leg, 0.212 s at one
thread and 0.212 s at eight, **30.6% of a threaded step** →
**[M-126](MEASUREMENTS.md#m-126)** — is **done upstream**, at 4.35× on eight
threads, bit for bit. Two things about how it closed are worth keeping.

**The reason it had been declined was a share measured on the wrong problem.**
Upstream's doxygen declined both loops at *"under 6% of the step"*, taken on
four **fixed-boundary, single-right-hand-side** cases. A bordered Newton applies
`J⁻¹` to `N + 4` columns, so the traversal is `O( elements × columns )` where
the integrator loops are `O( elements )` — the ratio inverts, and the same two
legs that are 6% there are 30.6% here. *A share is a property of a problem, not
of a loop.*

**AND IT IS NOW IN THE LIBRARY MEQ LINKS, AND MEASURED** →
**[M-135](MEASUREMENTS.md#m-135)**. `../mfem/install` is rebuilt from
`meq-integration` at `3424f33bd5`, the install is verified against the branch,
and an interleaved pair on a quiet machine — a pre-upgrade binary kept aside,
`libmfem.a` being static, so the library is the only variable — reads the leg at
**1.315 s at 1.00 cores against 0.503 s at 3.32 cores, 2.61×**, with every other
leg flat within 2% and the **whole run 1.25×**. Upstream predicted 1.24× from
MEQ's own arithmetic. `psi_ax` is `3.759851e-01` in all twelve runs.

**THE PREDICTION LANDED FOR A REASON WORTH KEEPING: THE SHARE WAS RIGHT EVEN
THOUGH THE FACTOR WAS NOT.** The leg factor here is 2.61× where upstream
measured 4.35× on their own fixture — MEQ's traversal is twelve calls of
fourteen columns and theirs is not — and the *run* number still came out where
the arithmetic said, because what that arithmetic rested on was M-126's
**share**. A leg speedup quoted from somebody else's problem does not transfer;
a share measured on yours does.

**AND IT IS A 5% LOSS AT ONE THREAD**, which a speedup table hides: 1.375 →
1.446 s on the leg and 6.983 → 7.121 s on the run at `OMP = MKL = 1`, the
OpenMP region and the colouring being paid for and unusable. Upstream's own
table has the same sign at the same size. A serial deployment of MEQ is
very slightly worse off for this change.

**MEQ'S OWN HALF OF THE THREADING IS ALSO MEASURED NOW, AND `THREADING-PLAN.md`
IS CLOSED.** Items A, B and C — the sparse border rows and the five plasma
element loops and the critical-point sweep — are worth **`1.16×` on the whole
DIII-D run at eight threads**, against the `1.14×` that plan predicted, and
`1.03×` at one thread where only the sparsity item has anything to give
→ **[M-137](MEASUREMENTS.md#m-137)**. The `before` arm is the pre-threading
commit carrying **one** cherry-picked line, `SetIntegratorsThreadSafe()`,
without which it exits 2 against today's library — both arms on the shipping
`libmfem.a`, so MEQ's own source is the only variable.

**AND THE BUDGET EVERYTHING WAS SIZED AGAINST HAS INVERTED** →
**[M-138](MEASUREMENTS.md#m-138)**. About **59% of a MEQ run now threads**
against the ~30% M-126 measured, and the largest single-threaded item was
`DarcyForm::Reconstruct()` at 18.5% — MFEM's, unmoved by MEQ's own threading,
and ahead of the NPC traversal it used to sit behind. **Upstream has since
threaded it**, worth `1.177×` on the whole run
→ **[M-172](MEASUREMENTS.md#m-172)**, so the largest single-threaded item MEQ
pays for is now `re-assembly` at about 6%, which MEQ owns: the line search
re-assembling a right-hand side that is linear in the exterior coefficients, and
it lives in `BORDERED-GLOBALISATION-PLAN.md` §12. **Size anything on this path
against M-172, then M-138, and not against M-126.**

Everything else MEQ has sent is closed. **A CLOSED REPORT NEEDS NO ENTRY HERE**:
what it changed is in the code with a test on it, or it is a measurement under an
`M-nn` anchor, and either way this file is the wrong place to re-tell it. Three
operational facts are worth the space and the rest is not:

* **UPSTREAM DOES NOT TRACK INTER-PROJECT CORRESPONDENCE UNDER `doc/` AT ALL, SO
  AN ABSENT DOCUMENT MEANS NOTHING.** `46ac3d3bc7` deletes every `*-FROM-MEQ.md`
  in one commit, on the principle that working notes exchanged with a consumer
  are not the library's documentation and that what they establish belongs in
  doxygen on the thing it is about. **That is a policy and not a verdict** —
  reading an empty `doc/` as "every request closed" is the obvious inference and
  is wrong. What a request established survives in upstream's own plan documents
  and in the doxygen; the correspondence does not.
* **SO THE ONLY TEST OF "LANDED" IS THE CODE MEQ BUILDS AGAINST.** Look for the
  symbol in `../mfem/install/include`, not for the document. That was always the
  only test MEQ could apply and it is now the only one that exists.
* **Which documents exist is a question for `git`, not for `ls`** — a listing
  reports whatever branch that tree is checked out on, which is not MEQ's to
  control and has changed under this file more than once. **And a request MEQ
  writes will never appear in `git` at all**: that tree's `.gitignore` carries
  `doc/*-FROM-*.md`, so a report lands as an untracked working-tree file by
  design and is read there. `git status` not showing it is the policy working,
  not the write having failed — check with `ls`, which is the one question `ls`
  answers better.

**And do not file findings against unfinished work.** A branch that exists is not
a branch that is done. Measure it if it is useful to know, keep the numbers in
MEQ, and wait to be told it is ready before writing anything into that tree about
how it behaves.

Rebuilding the install, after fetching whatever is wanted into `../mfem/mfem-src`:

```sh
cd /home/ian/projects/mfem
cmake --build build -j4 && cmake --install build
```

The finder reads `share/mfem/config.mk` from an install and `config/config.mk`
from an in-source make build, so either layout works and neither needs a flag.

`mfem/master` will not do. MEQ needs `DarcyForm` and the HDG integrators in
`fem/darcy/`, and for stage 5 it needs the transfer-path machinery in
`fem/darcy/extension_hdg.{hpp,cpp}` — which is an implementation of exactly the
curved-boundary technique the GS papers use.

**The old MFEM's HDG API is gone, and that is why this is a port and not a
recompile.** `HDGBilinearForm`, `HDGDomainIntegratorGS`, `HDGFaceIntegratorGS`,
`AssembleSC` and `Reconstruct` have no counterparts; the replacement is
`DarcyForm` with `EnableHybridization`. The mapping, term for term:

| Weak term | Old | New |
|---|---|---|
| `(R q_h, v)` | `HDGDomainIntegratorGS` | flux mass: `VectorMassIntegrator` with an `R` coefficient |
| `(ψ_h, ∇̄·v)`, `(q_h, ∇̄w)` | `HDGDomainIntegratorGS` | flux divergence: `VectorDivergenceIntegrator` |
| `⟨ψ̂_h, v·n⟩` | `HDGFaceIntegratorGS` | the transpose of `NormalTraceJumpIntegrator` — **not** the face integrators on `B`, see below |
| `⟨τ(ψ_h − ψ̂_h), w⟩` | `HDGFaceIntegratorGS` | potential mass: `HDGDiffusionIntegrator` |
| `⟨q̂_h·n, μ⟩ = 0` | `AssembleSC` | `EnableHybridization(trace_space, new NormalTraceJumpIntegrator(), ess_list)` |
| condense / reconstruct | `AssembleSC` + `Reconstruct` | `FormLinearSystem` / `RecoverFEMSolution` |
| `φ_h` on `Γ_h` (stage 5) | — | `HDGExtensionIntegrator` on the **flux mass form**, `C = R`, sign `+1` |

`miniapps/hdg/convdiff.cpp` in that tree is the worked example to copy from;
lines 445–469 build exactly the three spaces above.

**A trap in that table, measured.** The `⟨ψ̂_h, v·n⟩` coupling does **not** come
from `TransposeIntegrator(DGNormalTraceIntegrator)` added to `B` on interior and
boundary faces, which is the natural reading. Under hybridization
`DarcyForm::Assemble()` builds `B` through `ComputeElementMatrix()`, which sums
**domain integrators only**, so those face integrators are never assembled —
changing the boundary one's coefficient from `−2` to `−0.7`, and deleting the
interior one outright, does not move a single digit of the answer. The coupling
comes from the transpose of `NormalTraceJumpIntegrator`, supplied to
`EnableHybridization`.

The boundary-face integrator on `B` still has to be there, but as a **marker**:
`EnableHybridization()` reads `B`'s boundary-face markers to decide where to
register the flux constraint. Remove it and the Dirichlet faces get no
constraint at all — the error goes flat at 1.5e-1. So it is load-bearing for a
reason unrelated to the integral it appears to compute, which is exactly the kind
of thing to leave a note about.

**And `DarcyForm::GetOffsets()` returns three entries, not four** — it never
learns about the trace space. The miniapps get the 4-entry version from
`DarcyOperator::ConstructOffsets()`, which is not in `libmfem.a`, so
`GradShafranovSolver` builds its own.

## Post-processing, NPC, Newton and the linear solves → `CLAUDE_HDGGS.md`

`DarcyForm::Reconstruct()` at `k+2` and the per-element defect a global norm
cannot see; the NPC port and the parity gap that came back the other way; Newton
and the obligation it creates, including the bordered Newton that makes `ψ_ax` an
unknown; the Solov'ev coefficients; *A wrong Jacobian is invisible to a
convergence table*; why MEQ's Newton struggles and what `PicardThenNewton` is
for; and the three linear solvers with the threading measurements. Toroidal
flow moved on again, to `CLAUDE_FLOW.md`.

## Solution inversion → `CLAUDE_INVERSION.md`

`ψ(R, z)` to `R(Ψ, l)`, `z(Ψ, l)`: the critical points as roots of the **solved**
`q`, the predictor-corrector tracer, the flux-surface averages and the metric
trap, the Zernike disc basis and why `ρ = √Ψ_N`, the gauge-free fit, the
`(Ψ, θ)` output file with its per-`ψ` cache, and the `q(ψ)`-driven outer Newton.
`INVERSION-PLAN.md` is the plan and every stage is done.

## Traps

**Every run needs `MKL_NUM_THREADS=1`**, which every registered ctest sets. It
is a no-op wherever MKL is absent, so it costs a non-MKL machine nothing.

**`MKL_THREADING_LAYER=GNU` was the other half of this trap and is now gone**,
from the ctests, the performance harness and the driver alike. It guarded the
`libmkl_rt` **dispatcher** — what `/usr/lib/x86_64-linux-gnu/libblas.so.3`
resolves to on this machine — whose default threading layer silently corrupted
UMFPACK's BLAS-3: you got numbers, and they were wrong. MEQ no longer loads it,
building against its own SuiteSparse; see `CLAUDE_HDGGS.md`, *PARDISO and the MKL link line* for
why the variable was dropped rather than kept. **The trap is still live for
anyone using Debian's SuiteSparse**, and `../mfem-hdg-dev/CLAUDE.md` records it.

`MKL_NUM_THREADS=1` is the **factor of 140** one. Measured, one
process, nothing else running: `SolovievConvergence` **1.24 s** against
**177.83 s** at `=16`, and `FieldConvergence` 1.11 s against 198.35 s — at 273%
of *sixteen* cores, which is the shape of the answer: the threads are spinning
on barriers, not working. `SamplerConvergence`, which never reaches a direct
solver, does not move.

**THE CULPRIT IS NOT UMFPACK, AND THE FIRST VERSION OF THIS SECTION SAID IT
WAS.** The plausible story — UMFPACK calling threaded BLAS on the small dense
frontal matrices of a sparse LU, thousands of times, each call costing a fork
and a barrier — fits every whole-test fact and is **wrong**. Separated out,
`UMFPackSolver::SetOperator` degrades by about **40%** across the whole thread
range, never more.

**It is MFEM's element-local dense LU.** `ComputeH()` factors `A`, forms the
Schur complement and factors that, on blocks of order 10–30, once per element,
through `LUFactors` → LAPACK → MKL:

→ **[M-59](MEASUREMENTS.md#m-59)** — assembly + reduction · `MKL=1` · `MKL=2`

**`k = 2` does not move and `k = 3` degrades by a factor of forty**, putting
MKL's threading threshold between those two block sizes. The suite runs `k` up
to 4 and every element of every mesh pays it, which is why the whole-test
factors dwarf anything in the solver columns. **A production run at `k ≥ 3` with
MKL threading left on would be unusable and would look like a solver problem.**

**It is MKL and not OpenMP**, which the obvious cross settles: `FieldConvergence`
at `OMP=16 MKL=1` is **1.19 s at 99% CPU**, and at `OMP=1 MKL=16` is 183.54 s at
266%.

**Keep the wrong answer in view, because of HOW it was wrong.** It fitted every
fact then available — the collapse tracked exactly the tests that reach a direct
solver — and **no whole-test timing could have refuted it**, because in those
tests UMFPACK and `ComputeH()` always appear together. Only separating the
columns did. Same lesson as the Solov'ev coefficients: checking a story against
the evidence that suggested it proves nothing.

**A `libmkl_sequential` behind SuiteSparse hides all of this**, since MKL then
resolves sequential whatever it is asked for — see `CLAUDE_HDGGS.md`,
*PARDISO and the MKL link line*.

**AND THE TENSION THAT USED TO SIT HERE — `MKL_NUM_THREADS=1` BEING RIGHT FOR
MEQ'S SOLVER AND WRONG FOR PARDISO — IS RESOLVED BY THE ASSEMBLY MODE, AND THE
TWO DEFAULTS NOW MOVE TOGETHER.** Everything above is a measurement of the
*serial* element loop. Under `AssemblyMode::Threaded`, which is the default, the
element-local dense work is nested inside an active OpenMP region and MKL
suppresses its own threading there, so `MKL_NUM_THREADS` costs it nothing —
while the trace solve runs on the master thread outside every region and takes
every thread. **So PARDISO is the default trace solver wherever
`MFEM_USE_MKL_PARDISO` is set**, chosen by `defaultTraceSolver()`, which is
build-conditional for exactly the reason `defaultAssemblyMode()` is: the
constructor bypasses `setTraceSolver()`, and an unhonourable default would reach
`makeTraceSolver()` and throw out of a solve.

→ **[M-78](MEASUREMENTS.md#m-78)** — driver wall by solver × MKL · per-stage solve · what is free and what is not

The short form: **the solver change is worth 8–10% and is free**, user time
falling with it; **the MKL threads are worth a further 7–9% and cost 75% more
CPU**; and UMFPack moves the *other* way under them, 6% slower for 30% more CPU,
which is why the pair is chosen together and not separately. The 1.9× that
`setTraceSolver()`'s documentation quotes for PARDISO's threads is an isolated
trace-solve figure and does not survive to the wall clock.

**The registered tests still pin `MKL_NUM_THREADS=1` and must**, because the
bit-exactness assertions between the two assembly modes hold at that value and
not above it: at `MKL=8` the modes hand MKL different thread counts and a
blocked BLAS-3 reassociates, measured at 1.4e-15 in ψ and 1.8e-13 in the flux.
That is arithmetic, not a race, and `NpcThreadScaling` now says so in its own
diagnostic rather than sending the reader to hunt for shared scratch in MEQ's
integrators — which is what its message used to advise.

**This machine's GPU is for development, not for performance conclusions.**
There is an RTX 2070 SUPER with CUDA 13.3, which is enough to write and debug a
device path and to check that it gives the same answers. It is **not** a
representative part for MEQ: MFEM is built `MFEM_USE_DOUBLE`, and consumer
NVIDIA cards run FP64 at 1/32 to 1/64 of their FP32 rate where datacentre parts
run it at about 1/2. A GPU timing taken here can invert the conclusion a
production part would give, so treat local device runs as correctness evidence
and take timings elsewhere. `TODO` carries the detail. This is the same species
of warning as the threaded-MKL note under `CLAUDE_HDGGS.md`, *On SUNDIALS* — a measurement on this
machine that is not a measurement about the code.

**Never run a bare `make -j`, anywhere.** With no argument it is unbounded, and
this machine is WSL2 — the job count goes to the host's core count with a fraction
of the host's memory behind it, and the whole VM falls over rather than the build
merely failing. Always give a number: **`-j6` for MEQ's own build**, and for
the MFEM tree specifically **`make -j4`, never more**, since its translation
units are large enough that even 8 exhausts memory here. The two limits are
different because the trees are: MEQ's translation units are small enough that
six fit, and `cmake --build ... -j6` is the ordinary build command.

**And `make clean` after editing any MFEM header.** MFEM's
makefiles have no `.d` files and no header dependency tracking. That trap has
produced heap corruption in unrelated functions and "unimplemented" aborts for
methods that had just been added. MEQ's own CMake build tracks headers properly;
this applies only to the MFEM tree.

**AND `cmake --install` PRESERVES THE SOURCE FILE'S MTIME, SO A REFRESHED MFEM
INSTALL LEAVES MEQ SILENTLY UP TO DATE ACROSS A CLASS-LAYOUT CHANGE.** The
installed header carries the mtime of the file in `../mfem/mfem-src`, which is
when that tree was last checked out — routinely EARLIER than MEQ's existing
objects. Measured after one ordinary reinstall: **all 317 installed headers read
older than MEQ's objects**, `cmake --build build` reported nothing to do, and
`DarcyHybridization` had just gained a member. That links MEQ objects built
against the old layout to a `libmfem.a` built against the new one — an ODR
violation with no error anywhere, which is the failure MEQ's own header tracking
is supposed to make impossible. **Touch the installed headers after every
install** — `find ../mfem/install/include -name '*.hpp' -exec touch {} +` — and
check the dry run reports work: `make -n` in `build` going from 0 to 57 is the
confirmation. This is why "MEQ must then be fully rebuilt too" appears wherever
the install is refreshed; the mechanism is the mtime, not the flags.

**AND A FRESH CONFIGURE RE-SEARCHES FOR WHAT A LONG-LIVED BUILD DIRECTORY WAS
GIVEN BY HAND, WHICH IS HOW A BASELINE STOPS BEING AN ARM.** Configuring a
second MFEM build from the `MFEM_USE_*` flags and TPL directories alone — the
list a reader would write down — resolved **PARDISO to Debian's
`libmkl_sequential`** and **BLAS/LAPACK to Debian's MKL with `libiomp5`**, where
the reference build has oneAPI and `libgomp` throughout. The first fails to
link, which is the good case; the second links two OpenMP runtimes into one
process and runs. `BLAS_LIBRARIES`, `LAPACK_LIBRARIES`, `MKL_PARDISO_LIBRARIES`
and `MKL_PARDISO_INCLUDE_DIRS` are cached **results**, not inputs, and the same
is true one level up: **MEQ's own build directory caches the link line it read
out of `config.mk`**, so repairing the library and rebuilding MEQ in place left
Debian's MKL on MEQ's `link.txt`. Copy the resolved values out of the reference
cache, wipe rather than repoint, and **diff `MFEM_EXT_LIBS` token by token
before believing any comparison between two installs**.
→ **[M-172](MEASUREMENTS.md#m-172)**.

**AND ANYTHING A CMakeLists READS WITH `execute_process()` IS READ ONCE, WHEN
THE BUILD DIRECTORY IS CONFIGURED — WHICH IS HOW A BINARY STOPS KNOWING WHICH
TREE IT CAME FROM** → **[M-224](MEASUREMENTS.md#m-224)**. The commit was
captured that way and `configure_file`d into `meq_version.hpp`, so a long-lived
build directory stamped **every** binary it produced, and every `.nc` those
binaries wrote, with whatever `HEAD` was the day somebody first ran
`cmake -B build`. Measured at five commits stale, and the string is *plausible*
— a real commit of this repository, merely not this one.

**IT IS A PROVENANCE DEFECT AND THREE THINGS IN THIS FILE REST ON IT**: the
`git worktree` recipe for separating *this was already failing* from *I broke
it*; the device programme, which is a comparison between two builds throughout;
and `docs/running.rst`'s promise that the answer to *which equilibrium did this
plot come from* is a file that can be re-run. `cmake/GitCommit.cmake` runs as a
**build step** and writes through `copy_if_different`, so the relink happens
when the commit moves rather than once per build, and a `-dirty` suffix comes
from `git status --porcelain --untracked-files=no` — tracked files only, because
this repository's root holds hundreds of untracked run products and a flag that
is always set means nothing.

**AND THE FILE SAID NOTHING ABOUT THE BUILD AT ALL, WHICH IS THE LARGER HALF.**
`meq --version` now reads `MEQ 0.1.0 (6a078e2) (MFEM 4.10.1, double, cuda)` and
the `.nc` carries `meq_precision` and `meq_device_backends` beside
`meq_version`, because a double build and a single build are otherwise two
binaries with one name answering to floors three orders apart
([M-204](MEASUREMENTS.md#m-204)). `meq::realPrecision()` and
`meq::deviceBackends()` are in `src/meq/BuildInfo.cpp` **and not in its header**:
they are answers about MFEM's own preprocessor, and a precision macro tested in
a translation unit that has not included MFEM is not false but UNDEFINED — which
is `MachineFixedBoundary.cpp`'s trap, one section down.

**A convergence target scaled to `‖r₀‖` makes a good initial guess FAIL, and
the better the guess the more certain the failure.** MFEM's `NewtonSolver` stops
at `‖r‖ ≤ max(rel_tol·‖r₀‖, abs_tol)` with `‖r₀‖` measured at the iterate it was
handed. Warm-start from a converged answer and `‖r₀‖` is already small, so the
target shrinks with it — and past a point it falls below the round-off floor,
where nothing can meet it. Measured on Example 5 at `k = 3, n = 8`, restarting
from the exact answer: the residual reaches the floor in **two** iterations and
is then reported as a failure at the thirtieth.

```
   it 0   1.180694e-03
   it 1   8.782203e-13
   it 2   4.151551e-14      <- converged; the floor
   ...    ~4e-14 for 28 more iterations, then FAIL
```

The target was `max(1e-12 × 1.18e-3, 1e-14) = 1e-14`, under the 3.7e-14 this
problem can reach. Cold, `‖r₀‖ = 11.2` gives a target of 1.1e-11 and it converges
in four. **This is the failure mode that matters most for how MEQ will be used** —
moving to an adjacent equilibrium should cost one or two Newton steps, and it
instead threw away a converged answer.

`solve()` therefore takes the reference from the **cold** iterate — the Dirichlet
datum alone, where this solve would have started with no guess — and sets a pure
absolute target from it. `rel_tol` keeps exactly the meaning it always had, a
cold solve is bit-identical because there the reference *is* `‖r₀‖`, and only a
warm one changes: from failing to converging in one step. It costs one extra
residual evaluation, and only when a guess was set.

**`NewtonSolver` needs `iterative_mode = true`, and the failure is silent.** The
Dirichlet values ride in the iterate, so with `iterative_mode` left false
`NewtonSolver` zeroes `x` on entry and throws the boundary data away — silently,
because the residual is masked on exactly those rows, so nothing complains and
the answer is merely wrong.

**A `BlockVector` size mismatch that Release builds do not catch.** Stage 2
passed a three-block vector (flux, potential, trace) where `DarcyHybridization`
expects two: `GetOffsets()` stops at the potential, and `ReduceRHS` does
`darcy_rhs = b_t`, whose `BlockVector::operator=` calls
`mfem_error("Number of Blocks don't match")`. Found in stage 4 and
fixed by passing views over `darcy->GetOffsets()`, as
`DarcyOperator::ImplicitSolve` does — **a latent defect in the linear path, not
something Newton introduced.**

**THE EXPLANATION THIS ENTRY USED TO CARRY WAS WRONG.** It said the linear path
survived because its checks are `MFEM_ASSERT`, "compiled out with `NDEBUG`", and
concluded that a debug build was "worth doing now and then for this reason
alone". Neither half holds: `MFEM_ASSERT` is gated on `MFEM_DEBUG`, not on
`NDEBUG`, and the installed MFEM sets `MFEM_DEBUG = NO` — so those checks are
dead in **every** MEQ build and **a MEQ debug build revives none of them**. What
actually caught the semi-linear path was `mfem_error`, which is unconditional.
See *Coverage* above for the audit. The failure mode is still real; only the
advice was.

**A nonlinear potential mass and a linear one do not mix, and how it fails
depends on where the integrators sit.** With *face* integrators on the linear
form (they become `c_bfi_p`) MFEM aborts loudly: "Non-linear mass cannot work
with a linear constraint". With *domain* integrators it is silent —
`LocalNLOperator::AddMultDE` and `ConstructGrad` simply drop them. So leaving the
HDG stabilisation on `M_p` beside a nonlinear source aborts, which is the good
case; the silent one is the reason `buildForms()` documents both.

**Every GS-2 §4.2–4.5 source vanishes at `ψ = 0`, so the paper's own problem has
a trivial branch — and MEQ falls into it.** `F(R, 0) = 0` for eqs (24), (25),
(26) and (27); for (25) because `b = 2` makes the bracket `O(ψ²)`. With
homogeneous Dirichlet data `ψ ≡ 0` therefore *solves* the problem, and Newton —
which starts from the Dirichlet data — lands on it and stops in **zero
iterations** with an identically zero residual. Not a hypothesis: GS-1's
Algorithm 2 literally opens `ψ⁰ ; // Non-trivial initial guess`.

**`setInitialGuess()` exists** — in `Coefficient` and `GridFunction` overloads —
and it is what those four cases need; without it `prepare()` resets the iterate
and `solve()` calls `prepare()`, so no caller could work round the trivial
branch. They are posed with a non-homogeneous ramp that puts `ψ = 0` in the
interior. **Under NPC the guess seeds the potential and the trace and leaves the
flux block at zero**, which is a real difference; see `CLAUDE_HDGGS.md`,
*The flux is seeded*.

**`DarcyHybridization` freezes the element-local Newton's initial guess at
`FormLinearSystem()` time, and `ComputeSolution()` will not take another**, so
seeding the vector you hand `RecoverFEMSolution()` is inert. On an ordinary
Newton path that costs iterations; on anything that **differentiates the
residual by differencing it** it costs correctness — measured, it moved
`max ψ_h` from 0.896 to 3.84 for a perturbation of 9e-6. `formSystem()` is the
fix. Full account under `CLAUDE_HDGGS.md`, *The trap that cost the most*.

**`DarcyForm::Reconstruct()` RETURNS A DIFFERENT FUNCTION WHERE `∂F/∂ψ` VANISHES
UNLESS THE FIX IS IN.** It fails silently and per element, so a profile with a
flat segment corrupts part of the domain and leaves the rest exact, and a
whole-domain check misses it. Full account and the transferable lesson under
`CLAUDE_HDGGS.md`, *Post-processing is back*; the part that matters here is that **the fix is on `gf-hdg-dev` and on no
other branch**, so a `meq-integration` rebuilt without it silently loses this.

**`mfem::Mesh::GetElementTransformation( int )` HANDS OUT SHARED SCRATCH, AND
THREADING IT IS A SILENT WRONG ANSWER.** MFEM's own comment: *"The returned
object is owned by the class and is shared, i.e., calling this function resets
pointers obtained from previous calls."* It is one `IsoparametricTransformation`
member of the `Mesh`. Two threads evaluating in different elements overwrite
each other's transformation — no crash, no error, a point transformed by the
wrong element. The reentrant route is the `( i, IsoparametricTransformation * )`
overload into a thread-local.

**SIX CALL SITES IN MEQ TOOK THE SHARED OVERLOAD, WHERE A COUNT MADE BY EYE
FOUND THREE — THE UNDERCOUNT IS THE POINT.** Five in `src/meq/FluxSurfaces.cpp` — `setBandExtension`'s face
loop, `elementSize()`, `locate()`, and **both** branches of `extendField()` —
and one in `src/meq/CriticalPoints.cpp`. All six are **fixed**, each into a
function-local `thread_local` so a transformation held live across a call cannot
be reset underneath it, and every printed number in the four affected
convergence tests is **byte-identical** afterwards.

**AND THERE IS A FOURTH KIND OF SHARED SCRATCH THAT NOTHING HAD NAMED.**
`Mesh::GetBdrFaceTransformations( int )` returns the mesh's own
`FaceElementTransformations`, with the same hazard. The caller-allocated variant
signals failure by `GetGeometryType() == Geometry::INVALID` where the pointer
version returns `nullptr`.

**A `SubMesh` KEEPS A POINTER TO ITS PARENT, AND THREE FIXTURES LET THE PARENT
DIE. LATENT FOR MONTHS, THEN A SEGFAULT IN SOMEBODY ELSE'S CONSTRUCTOR.**
`mfem::SubMesh::CreateFromDomain( background, ... )` stores `&background`.
`ExtensionConvergence`, `FluxSurfaceConvergence` and `FreeBoundaryCoupling` each
built the background as a **local** and returned the `SubMesh` out of the
function, so the parent was dangling the moment they returned. That is undefined
behaviour from the first day and cost nothing for as long as nothing asked a
`SubMesh` where it came from.

`mfem::VertexConePath` now does. Its cone `C(x)` reads the **parent's** edges at
each vertex of `Γ_h` — `HasCone()` is documented as *"whether the mesh handed to
the constructor was a SubMesh with a parent to read edges from"* — so the
constructor walks freed memory and dies in `Mesh::GetVertexToVertexTable`. **It
presents as an MFEM regression and is MEQ's own bug**, which is the reason to
write it down: the backtrace names three MFEM frames and no MEQ frame at all.

The production code was never exposed — `Estimator` holds its background as a
member and `apps/meq.cpp` takes it by reference — so this is a fixture defect,
and the fix is a `static` pool that lives as long as the process. **Anything
that returns a `SubMesh` must return, or otherwise outlive, its parent.**

**AND THE TRIGGER IS NOW GATED, WHICH CHANGES WHAT IS LATENT AND NOT WHAT IS
OWED.** The constructor opens
`const SubMesh *sub = use_cone_ ? dynamic_cast<const SubMesh *>( mesh ) : NULL;`
— so it reaches the parent **only when the cone is asked for**, and `use_cone_`
is `false` by default with MEQ passing it at none of its construction sites. So
MEQ has never used `C(x)`, and this backtrace is not reachable on MEQ's own
paths today. **The cone being off is upstream's measured decision rather than an
oversight, and half the measurement is MEQ's**: it restricts every vertex at
every refinement, `π/2` becoming `π/8`, and changes neither the tiling residual
nor the flux rates — while it costs the far face's quadrature, which is the half
MEQ found. The cost is not a coverage failure but the SMOOTHNESS of
`ξ ↦ a( x( ξ ) )` along a face, and refining the RULE rather than the mesh is
what separates those two: the same cone-on sums fall to the instrument's own
floor at a 40th-order rule. *A caller who wants the cone must raise the face
rule with it.* The lifetime rule above stands on its own — a dangling parent is
undefined behaviour whether or not anything currently dereferences it, and what
the gate bought is time, not safety.

**AND IT IS THE SECOND TIME AN MFEM UPDATE HAS TURNED A LATENT MEQ DEFECT INTO A
CRASH IN MFEM'S OWN CODE.** The first was `SourceIntegrator`'s shared scratch
meeting the threaded element loop. Both were dormant contracts — *the parent
outlives the child*, *an integrator on a threaded loop is reentrant* — that MEQ
had never had to honour. Expect the next library update to find a third.

**AND A MESH ACQUIRES NODES WITHOUT BECOMING CURVED, WHICH COSTS `GridSampler`
TWO ORDERS OF MAGNITUDE AND CHANGES NO ANSWER.** The affine fast path keyed on
`Mesh::GetNodes() == nullptr`; `meq::FieldTransfer`'s constructor calls
`Mesh::EnsureNodes()`, because `FindPointsGSLIB` refuses a mesh without one. So
constructing a `FieldTransfer` dropped every later `GridSampler` on that mesh
onto the Newton route for the life of the process — **159× at 129² and 83× at
513²**, with every located node and every sampled value identical, so only a
clock can see it. The fix is to ask the geometry's **degree**: an order-1 nodal
field on a triangle is the same affine map the vertices are, checked against the
vertex array rather than assumed equal to it, since MFEM lets a caller deform a
mesh through its nodes while leaving the vertices behind.

**It was latent rather than live, and nothing was keeping it that way.** Both
`FieldTransfer`s in `apps/meq.cpp` are built on the stored guess's mesh or on
`previousMesh`, which is a **copy** — so the mesh the driver samples was never
the mesh gslib was set up on. That is an accident of how the adaptive loop
happens to be written, not a rule anybody stated, and it is the second entry in
this section where one component's requirement silently disabled another's
optimisation. → **[M-93](MEASUREMENTS.md#m-93)**

**AND THE THIRD IS A DECORATOR, WHICH IS PASS-THROUGH ONLY FOR THE METHODS IT
NAMES.** `mfem::Operator::ArrayMult()` has a base implementation that loops
`Mult()`, so it is never missing and never errors — it is merely slow. A
decorator that overrides `Mult()` and not `ArrayMult()` therefore un-blocks
whatever it wraps: `PardisoSolver`, `CuDSSSolver`, `MUMPSSolver`,
`SuperLUSolver` and `STRUMPACKSolver` all override `ArrayMult` to walk their
factors once for every column, and the base loop sends them `K` separate walks
instead, **with identical answers and identical call counts**. MEQ wraps its
trace solver in `TimedSolver` on every bordered path, so the blocked trace solve
would have been bought and thrown away at the wrapper. The rule generalises past
this one class: **when a library adds a batched entry point with a working
default, every decorator in the chain has to be revisited, and nothing will tell
you.** → **[M-98](MEASUREMENTS.md#m-98)**

**A RESIDUAL EVALUATION MAY NOT SIT BETWEEN `NPCGradient()` AND THE SOLVES IT
ENTITLES.** `NPCResidual()` is non-`const` on `DarcyHybridization` where
`NPCReduce()` and `NPCRecover()` are `const`, and under a moving support MEQ's
own `fieldResidual()` calls `refreshPlasmaComponent()` as well — so neither the
factored local blocks nor the *problem they were factored for* survive it by
contract. This is why the bordered step's differenced `psi_bnd` column flushes
its queue before differencing rather than solving across it. The failure it
guards has no symptom: a column solved against a Jacobian belonging to a
different plasma support is a plausible number, and Newton converges to
something.


**THE ONE THAT ACTUALLY BLOCKS THREADING IS `Mesh::FindPoints`, AND IT CANNOT BE
FIXED LOCALLY.** It loops over every element through that same shared
transformation *and* builds a vertex-to-element table on the way, so it is not
reentrant — an attempt to share one `ContourTracer` across threads aborts with
*"the axis is not in the mesh"*. Every entry point can reach it. The
**unconditional** per-surface call is gone: `traceFromAxis()` seeds the walk with
`CriticalPoint::element`, taking `theTracerClosesAndTheElementWalkDoesNotFallBack`
from 6 calls to 0 over its six surfaces and `SurfaceAverageConvergence` from 196
to 1 over the whole binary. What is left is the data-dependent
last resort inside `fitByAngle()` and `faceJump()`, so a shared tracer still
needs those unreachable or serialised.

**AND "THE TRACER REPORTS ZERO FALLBACKS" COULD NOT HAVE SEEN ANY OF THAT.**
`Contour::fallbackLocations` never counts a **seed** call —  `sampleAt()` hands
`sampleField()` a local counter and discards it on both overloads — so a
per-surface full-mesh scan sat under a test asserting that count is zero,
reading 0 before and 0 after while the real count moved. The numbers above are
from a breakpoint on `mfem::Mesh::FindPoints`, which is the only instrument that
can see them. `INVERSION-PLAN.md` §11.3 has the detail.

**`mfem::Mesh::FindPoints` is `O(elements × points)`** — a brute-force scan over
element centres. It caps sample-cloud sizes in any off-grid error measure.

**MEQ will be an early user of a thinly-tested MFEM combination.** Fixed-boundary
Grad–Shafranov is a Dirichlet problem, so the trace carries an essential BC; and
Newton makes it nonlinear. That pairing —
`DarcyHybridization::SetEssentialBC` together with a nonlinear reduced operator —
was broken on the MFEM branch until recently: `EliminateTraceTrueDofsInRHS`
returned early for nonlinear problems and the essential BC was silently ignored,
so the solver was solving a different problem than asked. **It is fixed**
(`fem/darcy/darcyhybridization.cpp`, and the reasoning is commented there), but
`../mfem-hdg-dev/CLAUDE.md` records that **no regression covers the combination**.
If stage 4 produces a converged-but-wrong answer near the boundary, look here
first rather than at MEQ's assembly.

**toml11's `find_or<double>` silently returns the default when the node is an
integer.** Paid for once. The original bug was that everything was read through
`.as_floating()`, so `RMin = 0` threw instead of converting — and the obvious fix
is worse than the bug. In toml11 4.4.0, `toml::find<double>` on an integer node
throws just the same, and `toml::find_or<double>(v, "RMin", 1.0)` **returns
1.0**, because a failed conversion is indistinguishable from a missing key. That
turns a loud failure into a silent wrong answer in the configuration, which is
the last place you want one. `Config.cpp` therefore reads numbers through its own
`asFloat()`, accepting `is_floating()` or `is_integer()` explicitly.

**A RESERVED KEY IS ONLY RESERVED ON THE PATHS THAT CALL THE REFUSAL, AND
LISTING IT AS ACCEPTED IS WHAT MAKES THE GAP SILENT.** `ProfileFile` and every
`<Profile>Variable` / `<Profile>Fit` name a facility that is not written —
several profiles read out of one NetCDF, which MaNTA produces. They are listed
in `rejectUnknownKeys` **on purpose**: `rejectUnknownKeys` runs first, so a key
has to be accepted before it can reach a message of its own, and "not
implemented yet" is a better answer than "is not a key of `[source]`".

That listing is also what turns a missing refusal into an accepted-and-ignored
key. `readEitherProfile()` carried the refusal inline; `Type = "mhd"` requires
`PPrimeFile` and `GGPrimeFile` and offers **no constant form**, so it reads them
directly and never reached it. Measured against the built driver, the same key
name gave two answers:

| | `PPrimeVariable = "Var0"` |
|---|---|
| `Type = "rotating"` | refused, naming the unwritten facility |
| `Type = "mhd"` | **parsed, ignored, ran on to open the text table** |

**Nothing tested the reserved keys on either path**, which is how it survived.
The refusal is now `refuseReservedVariableKeys()`, a free function both paths
call, taking a flag for whether a constant form exists so the advice at the end
of the message is true either way; and
`the_reserved_profile_keys_are_refused_on_every_source` asserts **the pairing**
rather than the refusal, since a single-path test would have passed throughout.
**Recorded because the same shape can recur** wherever one source type reads a
key directly and another goes through a shared helper. **And it recurred**, in
the same pair of source types: `ConfineToPlasma` was accepted for
`Type = "rotating"`, set its flag, and was consulted by none of
`meq::NormalisedRotatingSource`'s gates, so the key parsed and the solve was
unconfined — with XP-1's element-level fill running anyway, since that reads
`plasmaSupport()`. `CLAUDE_FLOW.md`, *`ConfineToPlasma` reaches the rotating
source*.

**And `ConfigError::what()` carried a `MEQ: ` banner while `apps/meq.cpp`
prefixes every line with the same thing**, so every configuration error read
`MEQ: MEQ: configuration error in …`. It is not new — it was `meq: meq:` for
months and nobody read it at that size; the rename to `MEQ` is what made it
visible. The banner is the caller's now. Every other exception MEQ throws
identifies itself by a namespace-qualified function name —
`meq::SplineProfile::fromFile` and its kind — so `ConfigError` was the odd one
out rather than the model.

**A BOUND ON HOW WRONG AN ESTIMATE MAY BE IS A BOUND ON AN ESTIMATE OF THAT
QUANTITY, AND SAYS NOTHING ABOUT AN ESTIMATE OF A DIFFERENT ONE.** The driver's
plasma-support loop is a fixed point over the support, so `apps/meq.cpp`'s own
comment is right that *"a bad estimate costs sweeps rather than correctness"* —
the next sweep re-decides from the answer just reached. That licence is what
made a second, informal reading of a constraint look safe to write, and the
second reading computed something else: `max psi` over the region a limiter
ENCLOSES rather than over the faces BOUNDING it. The region holds the magnetic
axis, so the estimate was `psi_ax`, the first freeze was posed with a span of
**−1.63e-05 where the true one is +4.76e-02** — negative, so the support
inverted — and the loop **converged**, from two starting points a whole
equilibrium apart, to one wall-hugging annulus in every printed digit.
`meq::limiterPolygonMaximum()` is the single implementation both now call.
→ **[M-165](MEASUREMENTS.md#m-165)**, and `CLAUDE_FB.md`, *The driver's first
support freeze read a different functional under the same words*.

**AND THE ESCAPE WAS WRITTEN DOWN BESIDE IT, AS A COST NOTE.** The same comment
said *"THE SurfaceAttribute BRANCH IS EXERCISED BY NO SHIPPED FIXTURE"* and then
explained why building one was not small. **A sentence naming what nothing
tests is a finding, not a caveat** — it is the same shape as the reserved-key
pairing above, where a refusal that only one of two source types reached was
listed as accepted on both.

**A `switch` OVER AN ENUM WITH NO `default`, AFTER A VARIABLE IS ALREADY
INITIALISED, SILENTLY KEEPS THE INITIALISER — AND `-Wswitch` DID NOT CATCH IT.**
`apps/meq.cpp:1634` is `TS traceSolver = TS::UMFPack;` followed by a `switch`
over `config->getSolver().traceSolver` with one case per choice and no
`default`. Add a `TraceSolverType` and forget the case, and every run naming it
**runs UMFPack instead**, converging normally to the right answer with no
diagnostic anywhere. Measured: a run with `TraceSolver = "gmres"` before the
case existed clocked 5.621 s against UMFPack's 5.547 s and PARDISO's 3.006 s,
and **the wall time was the only tell**.

`CMakeLists.txt:556` sets `-Wall -Wextra`, the subject is a scoped enum
(`Config.hpp:1136`), and the switch has no `default` — so `-Wswitch` should
warn and, tested directly with the case removed and a forced recompile, **it
does not**. Why is not established. Until it is, **a new enum value needs its
switches found by hand**, and the same shape exists wherever a driver maps a
`Config` enum onto a library one. → **[M-179](MEASUREMENTS.md#m-179)**.

**AND A BRACELESS `for` ADOPTS THE NEXT BRACED BLOCK AS ITS BODY, WITH NOTHING
IN `-Wall -Wextra` SAYING SO — NOT EVEN `-Wmisleading-indentation`.** Left in
`GradShafranov.cpp` by a session that was experimenting with the bordered rows:

```cpp
        flush();

        for ( int mode = 0; mode < nModes && mode < 3; ++mode )


        // rowDot() AND THE INVERSE TOGETHER, ...
        {
                LegTimer const denseTimer( ... );
                ...                       // the WHOLE bordered dense solve
        }
```

The loop's body is that block, so the bordered step's row contractions and its
dense solve ran `min( nModes, 3 )` times — and **zero times whenever
`nModes == 0`**, which is every bordered case without exterior DtN coupling.
`rotating-normalised` then stalls with `||r||` CREEPING UP at observed order
1.000, and the q-driven outer loop and the meshed limiter go with it.

**Tested rather than assumed: `g++ -Wall -Wextra`, `g++ -Wmisleading-indentation`
and `clang++ -Wall -Wextra` all compile the shape silently.** A braced compound
statement as a loop body is ordinary C++ and the indentation heuristic does not
fire on it, so **there is no compiler flag that would have caught this**.

**HOW IT SURVIVED IS THE HALF WORTH KEEPING.** The session that left it verified
on three example solves — `machine-f-diiid`, `limited-tokamak`,
`soloviev-nstx` — and five unit binaries, and every one of those is immune:
the first two carry exterior modes so the loop runs at least once, and the third
takes a single Newton step. **Three cases chosen for being quick to compare were
chosen against exactly the configurations the defect reaches.** Same shape as the
six-of-fifty-eight sample under the rename entry above.

**AND IT IS GREPPABLE, WHICH IS THE CHEAP HALF.** The shape is a braceless
`for`/`if`/`while` head whose next *code* line is a bare `{` with blank or
comment lines between them — Allman style puts the brace on the very next line,
so any gap at all is the tell, whatever the indentation. Thirty lines of python
over `src/meq`, `apps` and `tests` reads **zero** today. Run it after any
session that edited a control-flow head; it costs a second and no compiler
offers it.

**The instrument that settled it is the one this file already prescribes**:
`git worktree add <dir> HEAD --detach`, configure with `-DMFEM_DIR=` and
`-DTOML11_DIR=`, build `meq` alone, and run the failing example against both
binaries. HEAD reached `psi_ax = 1.039325e-01` in 7 Newton steps and the working
tree did not converge in 40 — which is the difference between *this was already
red* and *I broke it*, and no amount of reading the diff establishes it as
cheaply.

**Four tooling traps that cost one session six idle polling loops**, and the
reason they are worth writing down is that a polling loop which can never
terminate is *invisible* when the thing it polls for is also reported some other
way.

* **`cd` persists between Bash tool calls**, so a `cd x && ...` that assumes the
  repo root silently runs somewhere else. Use absolute paths.
* **Shell VARIABLES do not persist**, only the working directory does — the same
  trap read the other way round, and the more dangerous half. `L=...` in one call
  is gone by the next, so `until grep -q '^EXIT=' "$LOG"` with an unset `LOG`
  becomes `grep -q '^EXIT=' ""`, and **`grep -q` with no file argument reads
  stdin**. In a detached background process stdin never reaches EOF, so the loop
  blocks silently for the life of the session.
* **A completion marker appended to a log is not necessarily at the start of a
  line.** `echo "EXIT=$?" >> log` after a Boost test binary lands on the end of
  its final ANSI reset sequence, because Boost's coloured output does not
  terminate with a newline — so `grep -q '^EXIT='` never matches. ctest's output
  *does* end with a newline, so the same waiter works on `ctest` and hangs on a
  bare test binary, which gets diagnosed as "the run is slow". Drop the `^`.
* **AND THE WAITER MUST GREP FOR THE MARKER THE WRITER ACTUALLY WROTE**, which
  is the same trap one level up and bites hardest when two parties agree on a
  log file but not on its marker. Measured here: writers emitting
  `echo "CTEST EXIT=$?"` against waiters greping `^EXIT=` — the line begins
  `CTEST `, so the pattern never matches and **fifteen shells polled a finished
  run indefinitely**, at a `sleep` each, for as long as the session lasted. The
  run had succeeded; nothing was wrong but the pattern. **A shared log file is a
  shared protocol**: write the marker and the pattern in one place, or have the
  waiter match the loosest thing that cannot appear early — `EXIT=` unanchored.
* **`grep -c` exits 1 on zero matches**, so a background check reports failure
  spuriously.
* **`pgrep -f` AND `pkill -f` MATCH THE INVOKING SHELL'S OWN COMMAND LINE.**
  This is the one that has now fired three times in a day, twice on the same
  afternoon, and it fails in two different directions. A waiter built on
  `pgrep -f "cmake --install build"` matches *itself*, so it reports the install
  running for as long as the session lasts and the thing it waits for is never
  seen — the install had in fact finished minutes earlier. And `pkill -f ctest`
  in a shell whose own command line contains `ctest` **kills that shell**, which
  surfaces as an unexplained exit code 144 (128 + SIGTERM) rather than as
  anything mentioning the pattern. Use `pgrep -x`/`pkill -x`, which match the
  process *name*, or check for the artefact — a file's timestamp, an exit
  marker — rather than for a process at all.
* **AND `-x` SILENTLY MATCHES NOTHING WHEN THE BINARY'S NAME IS LONGER THAN 15
  CHARACTERS**, which is what the advice above sends you into. `-x` compares
  against the kernel's `comm` field, which is capped at 15 bytes, so
  `pgrep -x PlasmaEdgeConvergence` — 23 characters — never matches and never
  errors: it just answers *not running*. A waiter built on it falls through
  immediately, and what follows reads a **stale binary** and concludes something
  false about it. Measured here: a `--run_test` filter reported "no test cases
  matching" for a case that was in the source, because the waiter had returned
  before the link finished.

  **So both forms fail and in opposite directions** — `-f` never stops, `-x`
  never starts — and `MEQ`'s test binaries are nearly all over 15 characters, so
  the shorter form is the one that looks right here and is not. **Wait on the
  ARTEFACT**: `until grep -q '<marker>' out.txt; do sleep 5; done`, or let the
  harness's own completion notification do it. Five waiters spun for up to four
  hours in one session on `until ! pgrep -f 'build/tests/PlasmaEdge...'`, every
  one of them matching its own command line.

* **AND THE WRITER MAY NEVER WRITE AT ALL, WHICH NO PATTERN CAN FIX.** The
  bullet above is about a waiter and a writer disagreeing on the marker; this
  is the case where the writer is *gone* — killed mid-flight, timed out, or
  restarted under a different filename because the first attempt was being
  redone. The file then sits at its last byte forever and an unbounded waiter
  sleeps on it for the life of the session. Measured here: **three waiters on
  one such file, alive at 19.5 hours**, after the run producing it was killed
  and restarted writing somewhere else. They cost nothing but a `sleep` each
  and were invisible because every answer they were waiting for had already
  arrived by another route.

  **So give a waiter a BOUND**: `for i in $( seq 1 120 ); do <check> && break;
  sleep 30; done` surfaces a dead producer as *gave up* rather than as
  silence. **And when a long run is killed and restarted, kill its waiters
  too** — they are not addressed to the run, they are addressed to a filename.

**If a waiter is used, check it actually fired.**

**AN ELEMENT-LOCAL NEWTON EVALUATES ITS POLYNOMIAL OUTSIDE ITS ELEMENT ON
PURPOSE, AND ANYTHING ADDED TO THAT POLYNOMIAL HAS TO SURVIVE IT.**
`CriticalPointFinder::rootInElement()` lets the iterate leave the reference
element by up to 2 — that is how a root near a face is found — and a polynomial
is defined everywhere, so this cost nothing for as long as the field being
rooted was one. Under `COIL-SUBTRACTION-PLAN.md`'s split the field is
`q_p + q_c`, and `q_c` is a function on the OPEN half-plane: NaN on `R = 0`, and
`meq::coilGradPsi` REFUSES `R < 0`. On a half-disc machine, whose elements reach
`R = 0` exactly, the iterate leaves the half-plane and the run dies with
*"the field point radius must not be negative"* — three frames down, naming a
radius and nothing else. **A point off the half-plane is not a point of the
machine**, so the seam returns a bool and every caller abandons rather than
carrying on with `q_p` alone, which would root a different function. The
transferable form: **a term added to an extrapolating evaluation must be defined
wherever the extrapolation goes, or the evaluation has to be able to refuse.**

**AND THE SAME THROW CAME BACK AT A SECOND SEAM, WHERE THE FIRST BRANCH OF THAT
RULE IS THE RIGHT ONE — THE TERM *IS* DEFINED THERE, AND THE PARITY IS WHAT
DECIDES WHICH SEAM GETS WHICH ANSWER.** The exterior datum is assembled at
transfer-path TARGETS on `Γ`, and on a half-disc machine `Γ` is a semicircle
whose two endpoints lie exactly ON `R = 0` — so a target near an endpoint lands
either side of the axis by an amount that belongs to the path map and not to
whether the mesh is valid. Measured, `( R, z ) = ( −1.0984e-03, 3.4 )`, about 1%
of `h`. **`psi = R A_φ` is EVEN in `R`**, both factors changing sign together, so
`meq::ConductorField::psi` reflects and the point is answered exactly rather than
refused or clamped. `gradPsi`, `flux` and `poloidalField` still refuse, because
`∂_r psi` is **odd** where `psi` is even and the continuation therefore differs
in sign **between the two entries of one vector** — which is why
`CriticalPointFinder::totalFlux` above keeps its bool and `totalPotential`, the
scalar sibling, needs no guard at all. **The two seams differ by the parity of
what they evaluate, not by how carefully each was written.**
→ **[M-156](MEASUREMENTS.md#m-156)**.

**AND NEVER COMMIT A FILE A SUBAGENT OWNS ON THE STRENGTH OF ITS COMPILING.**
an agent was searching for a fixture in
`FreeBoundaryCoupling.cpp`, the file was staged after checking it built, and the
commit captured 99 lines of the agent's temporary `zzExperiment` scaffolding.
**Scaffolding compiles.** The check that was wanted is a `git diff` against what
the file is supposed to contain, or a `grep` for the marker names an agent uses
for its own workings — not a build. Staging early does lock the content against
a later write, which was the right instinct; it locks in whatever is there at
that moment, which was the wrong assumption. **Ask the agent to confirm the file
is clean, and read the diff.**

**AND NEVER `cmake --build` WHILE `ctest` IS RUNNING.** 
for a comment-only header change, which relinked ten test binaries *including
the one ctest was starting*. Nothing errored — replacing a running executable's
file leaves the running process on its old inode — so the run continued and
would have reported a pass over a mix of old and new binaries. A suite result is
only an acceptance measurement if every binary in it is the one being committed;
that run was killed and re-run rather than believed.

**A MECHANICAL IDENTIFIER RENAME IS NOT VERIFIED BY A BUILD AND A LINT, AND
`-Wshadow` IS THE INSTRUMENT THAT SAYS SO.** Renaming the cylindrical radius
from `r` to `radius` across the tree touched about 1,300 sites; the tree built
clean, the `naming` ctest passed, `sphinx-build -W` passed, and **six sites had
changed an answer**. C++ is what makes it silent: an inner declaration may
shadow an outer one of the same type, and a declaration may initialise itself.
So `double const r = radius*sin( theta )` inside a loop over `radius` became
`double const radius = radius*sin( theta )` — legal, and `-nan` at runtime — and
a lambda parameter renamed onto a captured circle radius quietly redirected
`phi = radius*radius - dr*dr - dz*dz` from the circle to the field point, which
is a different level set that still converges.

Configure a second tree with `-DCMAKE_CXX_FLAGS="-Wshadow -Winit-self"` and grep
its build log for the renamed identifier. **Build it CLEAN**: an incremental
rebuild only re-warns about the files just edited, so the second look at a log
is systematically narrower than the first and says so nowhere.

**THE WORST ONE WAS IN THE LOG FROM THE FIRST RUN AND WAS FILTERED OUT, WHICH IS
THE PART WORTH KEEPING.** The warnings split into `shadows a previous local` and
`shadows a member`. The member ones were checked in `GradShafranovSolver`, whose
`radius` is an `mfem::FunctionCoefficient` — no conversion to or from `double`
exists, so a stolen binding cannot compile and the shadowing is provably inert —
and that finding was then applied to the CATEGORY. `meq::analytic::PlasmaEdge`
has a member `radius` too and it is a **`double`**: the disc's radius `a`. Its
fifteen sites were in the lines already filtered away, so
`levelSet()` returned `R² − dr² − dz²` instead of `a² − dr² − dz²` and **the
analytic fixture stopped being the equilibrium it claims to be** — `f()` off
from `−Δ*ψ` by 3.2e+01 against a 1.0e-08 floor. *A reason that makes one
instance safe is not a reason about the category, and the grep that encodes it
is applied to the category.*

**AND THE FIXTURE'S OWN TEST COULD NOT REPORT IT, BECAUSE THAT BINARY WAS
ALREADY THE DELIBERATE RED.** `PlasmaEdgeConvergence` is red by design while the
staircase `Γ_{p,h}` stands, so "PlasmaEdgeConvergence fails" carried no
information. What settled it was **building `HEAD` in a `git worktree` and
diffing the failure inventories**: 4 errors in one case against 318 across seven.
`git worktree add <dir> HEAD --detach` plus `-DMFEM_DIR=` and `-DTOML11_DIR=`
pointed at the main checkout is the whole recipe, and it is the only check that
distinguishes *this was already failing* from *I broke it*. **A red test cannot
also serve as a regression on something else**, and a tree-wide rename touches
every test including the red ones.

**WHAT NO COMPILER CAN SEE IS `r` THAT NEVER MEANT A RADIUS**, and the tree had
nine: loop counters for a **repeat**, a **refinement level** and a table **row**,
a range-for over `Reading`s, two `r = rate( ... )`, and the `||r||` residual norm
— which the driver had been printing inconsistently as `||R||` in the bordered
step table while the four test harnesses printed `||r||`. Ask `git` for the
originals instead: `git grep -nE "\b(double|int|auto)\s+(const\s+)?r\s*="
HEAD` and read the initialisers.

**THE PYTHON HALF FAILS DIFFERENTLY AND `py_compile` IS FREE.** `!r` in an
f-string is the **repr conversion**, not an identifier, so `{float(v)!r}` became
`{float(v)!R}` and **two tools stopped parsing** — `SyntaxError: f-string:
invalid conversion character`. Nothing in `ctest` builds the tools, so nothing
would have said so. The same pass also rewrote the `k3r0` case-stem scheme that
`race.py` and `MEASUREMENTS.md`'s published tables share, and renamed two result
dicts to `R`. `for f in $(git ls-files '*.py'); do python3 -m py_compile "$f";
done` catches the first class outright.

**AND A FENCED BLOCK HOLDING CAPTURED OUTPUT IS DATA, NOT PROSE.** Three
transcripts were rewritten into records of something no program ever printed —
a Newton table's `||r||` header, a mesh check's `minimum r over the mesh`, an
X-point locator's `r = 0.699700`. Extract every fenced block from `HEAD` and
from the working tree and diff them pairwise; the maths blocks SHOULD move under
a convention change and the transcripts must not, so the two are told apart by
asking which program prints the line.

**A VERIFICATION THAT COMPARES PRINTED VALUES MUST NAME WHICH BINARIES IT RAN.**
This rename was reported verified on "six value-printing tests byte-identical to
baseline", which was true and did not cover any of the cases that had moved. Six
of fifty-eight is a sample, and a sample chosen for being easy to diff is chosen
against exactly the tests whose output is geometry.

**The four defects MEQ reported to MFEM are closed**, though
`HDG-DEFECTS-FROM-MEQ.md` itself is **not** gone — it is alive on `gf-hdg-dev`
and `gf-hdg-subdomains-dev` and deleted only on the symbolic-reuse line, which
is the modify/delete conflict in the merge recipe — a working-tree listing taken
while that tree sits on the branch that deletes it will say it is gone. Checked
one at a time, because "closed" arrived in four different ways: `Reconstruct()`
on a singular local matrix was **fixed** (MEQ's test flipped red to green on
it); `φ_h` unreachable after a solve was **fixed**, as
`mfem::TransferredDatumCoefficient` in `extension_hdg.hpp` — MEQ still calls
`setTransferredBoundary()` — and rebuilding `η₅` on it, which this sentence used
to say was "MEQ's work and is not done", **is done**: see *A separate `η₅`
problem on the extension path*, where it converges at 2.78 against the pinned
version's 0.40; `ReconstructFluxAndPot()` lifting only domain integrators was
**withdrawn as not a defect**, and MEQ had measured it harmless, which was the
right answer for the wrong reason; and `ComputeHDGFaceEnergy()` ignoring an
installed `HDGStabilization` is **FIXED — `StabValue()` is called in the
function body**, with a comment saying what it returns when no hook is
installed.

**THAT LAST ROW SAID THE OPPOSITE UNTIL RECENTLY AND THE REASON IS WORTH
KEEPING.** It read *"a code read today still shows it computing the `{h⁻¹Q}`
form with no call into an installed hook"*, and that code read used a **fixed
line range** which stopped short of the call — the function is 194 lines and the
call is at 153. Re-read by extracting the function body,
`awk '/ComputeHDGFaceEnergy/,/^}/'`, it is there. A line range is an instrument
and the function is the answer, which is the same species as every other
instrument-not-answer finding in this file — in the one place where it produced
a false claim about somebody else's code. **All four defects are now closed and
the report is deletable** — and upstream has since deleted it from every working
branch, along with the three other closed ones.

## The `freegs4e` benchmark → `CLAUDE_FB.md`

The comparison against an independent free-boundary code — a different algorithm,
von Hagenow Green's functions and finite-difference Picard in Python — plus the
fixed-boundary rehearsal that preceded it, the four conversions each of which
converges to a wrong answer, and the four times an instrument was mistaken for a
result. `tools/freegs4e-benchmark/` is the harness and `docs/validation.rst` is
the user-facing account.

## Testing stance

**A test asserts the behaviour that is wanted, and fails until it is there.**
Never the reverse. A test that asserts a known defect passes while the defect
stands, and that makes a green suite compatible with a broken solver — which
destroys the only property the suite is for: **100% green must mean there is no
lurking defect.** So a defect gets a *failing* test naming it, not a passing test
recording it.

**This was got wrong, corrected, and then vindicated.**
`thePostProcessedPotentialIsCorrectWhereTheJacobianVanishes` was first written
to assert the *corruption*, and passed. Rewritten to assert that `ψ*` is a
post-processing of `ψ_h` on every element, it went **red** — and stayed red
until MFEM fixed the reconstruction, at which point it flipped green on its own
and became the regression. **That flip is what said to put the driver back on
`ψ*`.** A test asserting the defect would have said nothing.

Two consequences. The suite is expected to be **red while a known defect
stands**, and that is the intended signal rather than a broken build — the
failing test's message is the record. And a failing test must say what fixes it,
because a red suite nobody can action is one people learn to ignore.

**A convergence table finds what unit tests do not.** This is the lesson from the
sibling MFEM branch, where reproducing a published table turned up three library
defects a passing unit suite had walked past for years — two of them silent.

**AND A STAGED PLAN'S ACCEPTANCES CAN EACH BE SHARP AND JOINTLY MISS THE CELL
WHERE TWO FEATURES MEET.** The coil-subtraction campaign has an acceptance per
stage, every one of them an identity or a rate rather than a tolerance. CS-3's
is a conductor inside `Γ` with **no plasma**; CS-4's driver acceptance is a
**fixed-boundary box**. Neither has an X-point border and neither has a plasma
current, so when the first machine was posed free boundary with its conductors
subtracted, **both of those borders turned out to read the remainder rather than
the total** — and the second failed in a way that named the wrong thing entirely,
`int F/R` coming out exactly zero and the run reporting a singular bordered
Jacobian. → **[M-148](MEASUREMENTS.md#m-148)**. The acceptance that was missing
now exists and is green at a rate of 3.78, → **[M-151](MEASUREMENTS.md#m-151)**. The same shape as the `2 × 2`
cross under *Commands* above, one level up: there a one-key experiment separated
two hypotheses only if everything else was where you thought it was; here a
per-stage acceptance covers a stage only if no other stage's feature is
interacting with it. **When a plan's stages compose, the acceptance that matters
is the one on the composition**, and enumerating consumers — CS-4's staging entry
names nine — is exactly the kind of list that looks complete and is not.

The specific hazard here is worth stating plainly, because it defeats the obvious
test design: **a wrong sign convention converges, at the right rate, to the wrong
function.** An order-of-accuracy study cannot catch it. Only a comparison against
a closed form can. That is why stage 3 is pinned to
`refs/HDG-GradShafranov-Adaptive.pdf` §4.1 — an exact NSTX Solov'ev solution with
`A = −0.52` whose twelve coefficients are published to 15 digits — and not to a
self-convergence study.

So the test ladder is, in order:

1. **Unit**, in `tests/unit/` — Boost.Test. Config parsing, spline value *and*
   derivative against closed forms, `dFdPsi` against a finite difference.
   And, in `NewtonConvergence.cpp`, the *assembled* Jacobian against a central
   difference of the assembled residual — which is a different and stronger
   check than the one on `dFdPsi` alone.
2. **An exact solution**, as an absolute-error regression.
3. **Convergence rates**, in `tests/convergence/` — asserting a *rate*, not a
   tolerance: `k+1` for `ψ` and `∇ψ`, `k+2` for `ψ*` at `k ≥ 1`, over `k = 1…4`
   and several dyadic refinements, against Tables 1–5 of the first paper and
   Table 1 of the second.

Expect rates to stop improving beyond `k ≈ 5` or `h/8`: the papers report
round-off dominating there, so a table that flattens at that point is behaving
correctly and is not a bug to chase.

**Unfitted convergence needs a two-tier rate assertion.** On a fitted mesh a
single-pair rate is tight enough to assert at 0.15 of slack. On the extension
path it is not: `Ω^h` is the union of background elements lying inside `Ω`, and
*which* elements those are is not a smooth function of `h`. At `k = 3` the error
is not even monotone in the mesh count — 2.02e-5, 3.03e-5, 3.88e-5 at
`n = 6, 7, 8` — and `{12,24,48,96}` is a worse sequence than `{8,16,32,64}`.
That is geometry, not the transfer: all four path families give the same
numbers. So `ExtensionConvergence` allows 0.30 per pair and asserts 0.15 on the
rate across the whole sequence.

**Newton is not monotone on a stiff source, so assert on the best triple, not on
every one.** Example 5 gives a clean four-step quadratic run; the GS-2 sources do
not. A typical history wanders for four to eight steps with observed orders of
0.02, 0.33, −0.37, −5.9 and *then* enters the quadratic regime and finishes
1e-9 → 1e-13 in two. §4.5 at `k = 2` even rises mid-run, 1.7e-2 → 5.3e-1 at
iteration 3, before recovering. Print the whole history; assert on the best
triple above the round-off floor.

**A self-convergence study on a rectangle cannot demonstrate `k+1`, and the cap
has nothing to do with the physics.** Measured on a *Solov'ev* source — constant
in `ψ`, `dF/dψ ≡ 0`, one Newton step — with homogeneous data on the benchmark
box, the self-difference rate is 2.00/2.86/2.96 in `ψ` and 1.88/2.25/2.18 in `q`:
flat at about 3 and 2.2 from `k = 2` on, with nothing nonlinear anywhere. That is
the `r² log R` corner term of a right angle — the solution sits in `H^{3−ε}` and
its gradient in `H^{2−ε}`, and no polynomial degree recovers it. The
exact-solution studies are immune because there the datum *is* the trace of a
smooth solution.

It is worse on a polygon approximating a curved boundary. GS-1 Example 6 on a
fixed 40-gon gives 2.12 in `ψ` at `k = 3`, against 1.995/3.000/4.000 for a
Solov'ev control **on exactly the same meshes** — the singular exponent at an
interior angle of `π − 2π/40` is `40/38 = 1.0526`, and interior pollution goes at
about twice that. **This is precisely the difficulty GS-2's curved-boundary
technique exists to remove**, and it is worth knowing before anyone designs
another fitted-polygon study.

**Mutation-test a suite you are relying on.** `ProfilesTests` and `SourceTests`
were checked by deliberately introducing fifteen defects — a dropped `r²`, `p'`
replaced by `p`, `μ₀` on the wrong term, the Solov'ev sign flipped, an
off-by-one in the interval lookup at a knot, `write()` truncated to six
significant figures — and confirming each was caught. That is a cheap way to
find out whether a green suite means anything, and worth repeating for the
convergence tests, where the risk of a test that passes regardless is higher.

Analytic solutions live in `tests/analytic/`, and they form a deliberate ladder
in how the source depends on `ψ` — which is to say, in what they demand of a
Newton Jacobian:

| fixture | source | `∂F/∂ψ` | what it catches |
|---|---|---|---|
| `Soloviev.hpp` | constant in `ψ` | `0` | the discretisation alone; Newton must converge in one step |
| `McCarthy.hpp` | linear in `ψ` | `T`, a nonzero constant | a Jacobian missing its mass term |
| `ManufacturedNonlinear.hpp` | `ψ²` and `e^(−ψ)` | varies | a Jacobian that is present but wrong |
| `SimilarityExponential.hpp` | `f₀e^{nψ}(1+εr²)` | `nF` | the same, against an *exact* rather than manufactured solution |

The last two are both nonlinear but are not redundant, and the difference is
worth keeping straight. **`ManufacturedNonlinear` is manufactured**: a convenient
`ψ` was chosen and `F` built to fit it, so its `F` is not of a form any physical
profile produces. **`SimilarityExponential` is exact**: the free function is
chosen — `f(u) = f₀e^{nu}` — and the solution follows from a similarity
reduction (Kaltsas & Throumoulopoulos; `refs/Refs.md`). So it tests the solver
against a nonlinear equation somebody might pose, rather than one reverse
engineered from an answer. It is also a cleaner shape: the whole `ψ`-dependence
is one exponential, so `∂F/∂ψ = nF` exactly, and at `n = 3` that exponential
varies by a factor of 29 across the benchmark rectangle.

The middle rung earns its place: with `∂F/∂ψ = 0` a Solov'ev run cannot tell a
correct Jacobian from an absent one, and in the nonlinear case the algebra is
messy enough to hide a factor. McCarthy's `∂F/∂ψ` is a single constant, so the
mass term is either there or it is not.

Each fixture checks its own transcription rather than trusting it: `deltaStarFD()`
recomputes `Δ*ψ` by central differences and the suite asserts it against `−f()`.
For McCarthy's eighteen-term Bessel and Neumann expansion that agrees to 3.6e-6,
and its analytic gradients match finite differences to 5e-9.

### And a rung above them all with no closed form at all: six machines, fixed boundary

**THE LADDER ABOVE IS ANALYTIC AND THE TREE HAD NOTHING BETWEEN IT AND THE WHOLE
FREE-BOUNDARY APPARATUS.** Every machine-geometry case here is free boundary —
the half-disc reaching the axis, the DtN coupling, the two borders, the moving
plasma edge — and the fixed-boundary examples are the fixtures above or one
768-element box. `examples/fixed-*.toml` are that missing rung: **six real
machines posed fixed boundary, 1,086 to 102,840 elements**, each with a
tabulated `p'` and `gg'` and an independent code's answer to check against.

**Ian's recipe, and it is what makes them a REFERENCE rather than only a source
of geometry.** Solve the machine free boundary in `freegs4e` at `257²`; fit MXH
to its `psi_n = 0.95` surface; hand MEQ that curve with `psi = 0` on it. `psi` is
defined only up to an additive constant, so that is the **same equilibrium in a
different gauge**, and the free solution shifted by the surface's flux is the
answer MEQ must reproduce. `tools/freegs4e-benchmark/make_case.py` is the
generator, `make_fixed.sh` the driver.

**AND `psi_n = 0.95` IS WHAT LIFTS THE ORDER CAP, WHICH IS THE WHOLE REASON
THESE EXIST.** On every free-boundary machine case `ConstrainPaxisIp` at
`alpha_n = 1.2` makes `p'` a fractional power of `1 − Psi` at the plasma edge and
the rate is capped at about **1.2** whatever the degree. At 0.95 that
singularity is outside the domain and the profiles are analytic across Γ — so
`psi*` reaches **`k + 2` on machine geometry**, measured at 3.15, 4.05 and 5.50
for `k = 1, 2, 3`. The separatrix is also excluded for a second reason worth
keeping: it passes through an X-point and has a **corner**, and MXH is a
truncated Fourier series that cannot turn one.

→ **[M-147](MEASUREMENTS.md#m-147)** — the six, their cost, their agreement
with `freegs4e` · the rate table · and the cure for
[M-61](MEASUREMENTS.md#m-61)'s three roots

**THE COMPARISON IS AGAINST [M-60](MEASUREMENTS.md#m-60) AND NOT AGAINST THE
FREE-BOUNDARY ONE**, which is worth saying because the machines are the same
and the two are easy to confuse. M-60 IS this comparison — `CLAUDE_FB.md`'s
fixed-boundary rehearsal — at `psi_n = 0.90`, a `129²` reference and ten
harmonics, reading 1.5e-04 to 6.8e-03. Moving to 0.95, `257²` and twenty takes
it to **5.0e-05 to 1.1e-04**, and [M-62](MEASUREMENTS.md#m-62) is the study that
said where that order would come from: it swept the grid and the harmonic count
on one machine and named the limiter at each step, the contour, then the fitter,
then MEQ.

**AND THE ROOT SELECTION IS M-61'S AND NOT NEW.** What is new is that
`Normalised = true` removes the upper root rather than navigating around it, so
these need no seed; and one sentence of M-61 is qualified — a ramp at 12× the
axis height *does* reach the physical root, so the branch is reachable and the
selection is **not monotone in the amplitude**. `CLAUDE_FB.md`, *Root selection
is the finding, not the agreement*, carries both.

**They carry no conductors and that is structural.** Γ is the plasma edge, so
every coil is outside the computational domain and its whole influence is a
Dirichlet datum that is identically zero — `[conductors] Model` has nothing to
act on. A fixed-boundary case that **can** measure coil subtraction needs a
domain larger than the plasma, which is what `examples/coils-rectangle.toml` is
and why its split reads exactly 1.00×: that box was never graded to its coils.
`TODO` carries what is still owed.

`tests/convergence/MachineFixedBoundary.cpp` is the acceptance, and it drives
the **shipped binary on the shipped file** on purpose: what rots about an
example is the FILE, and a library-level fixture would re-implement the curved
path and notice none of it.

## Layout

```
src/meq/     the library. Config, Profiles, Source, SourceFactory,
             GradShafranov, BoundaryShape, Estimator, Field, Sampler,
             WarmStart, Output -- all ported and all under the naming check.
             ConductorStore is psi_c ON DISK and is MFEM-free on purpose: it
             takes plain vectors, so the format is unit-testable without the
             library and netcdf stays inside two .cpp files. It is a CACHE and
             not a field -- values at the solver's own points rather than an
             mgrid-style grid -- and COIL-SUBTRACTION-PLAN.md CS-6 says why
             that trade goes this way and what the refusal has to guard.
             CriticalPoints, FluxSurfaces, SurfaceAverage, SurfaceFit and
             Zernike are the solution-inversion work and were written here
             rather than ported; see that section.
             DeviceBasis and AxisSearch are the keystone of the device
             programme and its first consumer. DeviceBasis evaluates an L2
             grid function at an ARBITRARY reference point on a straight
             triangle, which MFEM cannot do on a device at all --
             FiniteElement::CalcShape is a host virtual -- and which every
             SEARCH needs where an integral is served by GetDofToQuad()'s
             fixed rule. It is a change of basis built once per space in
             phase B. AxisSearch is what it was written for: the magnetic
             axis found by rooting every element on the device, with the host
             keeping the selection, the conductor filter and the merge --
             the half carrying the physics. M-230, M-231, M-232.
             ResidentIterate is the rest of the iterate's readers: the five
             plasma element loops, the psi_bnd row's ordered dot and the cold
             state's mask, each reproducing its host loop TO THE BIT on the
             host's own radii and weights. Gated like AxisSearch, so inert
             until the iterate stays on the device. M-237.
             Crossings and StageChecksum are the two DEVICE instruments and
             neither is part of the solve. Crossings counts host reads of a
             device-written state and, separately, how many of them MOVED
             data -- read the transfer column, not the read column. It counts
             the UPLOAD direction too, which it did not until phase B needed
             it: countedDeviceRead() asks HostIsValid() && !DeviceIsValid()
             before the Read(), and phase B is defined as "one upload", so an
             instrument seeing downloads alone could not measure it. And every
             crossing is attributed to a DEVICE-PLAN.md 3.3 residency phase
             through meq::PhaseScope, because a transfer in phase B is the
             design and the same transfer in phase C is the defect, and one
             total is an average over populations with opposite meanings.
             M-226. Stage-
             Checksum is the other half: it hashes a buffer at a named stage
             and reports its VALIDITY FLAGS BEFORE the read, so `H- D+` names
             a stage where anything not going through a funnel gets stale
             data. It prints only under MEQ_STAGE_CHECKSUM, an environment
             variable rather than a key because a walk has to be reachable on
             a shipped examples/*.toml without editing it. M-206, M-207.
apps/        drivers. meq.cpp, and MEQ_BUILD_APP defaults ON -- plus
             meq-run.in, the mesh-then-solve wrapper, which CMake
             configures and file(GENERATE)s into the build directory as
             `meq-run`. It is not runnable from apps/ and is not meant
             to be: it is a build product exactly as meq is, with both
             the binary's path and tools/'s baked in so it searches for
             nothing. See "MEQ is runnable" above for why it is a script.
tests/       unit/ (Boost.Test), convergence/ (rate assertions),
             analytic/ (closed-form solutions used by both),
             performance/ (TraceSolverScaling + scan.sh for the LINEAR path,
             NpcThreadScaling + npc-scan.sh for the NONLINEAR one,
             NewtonStepProfile + newton-step-profile.sh for WHERE a step's
             time goes, and InversionScaling -- all built, NONE a ctest,
             because every number in them is a timing. The NPC one exists
             separately because a linear solve never calls MultNL(), so it
             cannot see the loop AssemblyMode::Threaded now spends most of
             its time in. The PROFILE one exists separately from THAT because
             the two want different statistics: a best-of is right for a
             minimum wall time and a median for a leg SHARE, and it splits
             the step into legs by call site -- residual, gradient, trace
             factorisation, backsolve, plasma fill -- with `other` left as
             an explicit remainder so nothing the legs miss inflates one of
             them. GradShafranovSolver::StepProfile is the accessor and is
             always on)
tools/       plotting and visualisation. plot_equilibrium.py reads BOTH
             NetCDF files -- the (R, Z) grid and the (Psi, theta) flux
             surfaces -- and tells them apart by their own dimensions rather
             than by the name, so a --what belonging to the other one is
             refused rather than ignored. On the surfaces file it draws the
             band PER NODE and shades the two cut regions, which is the whole
             of why it exists: a traced curve cannot be read for which of its
             arcs is solved data, and a family over Psi_N in [ 0.05, 0.95 ]
             drawn on its own range looks like an answer rather than like a
             cut. tools/README.md says which of the four output formats
             goes with which reader, and why they are not interchangeable.
             freegs4e-benchmark/ is the independent-code comparison.
             desc-benchmark/ is the SECOND one, against DESC -- a
             fixed-boundary inverse spectral code, so it races the
             examples/fixed-*.toml machines rather than the free-boundary
             ones. Its profiles are labelled by TOROIDAL flux, which is not
             known until the equilibrium is, so race_desc.py poses each
             resolution self-consistently UNTIMED and clocks only the one
             forward cold-start solve that follows. M-157.
             chease-benchmark/ is the THIRD, against CHEASE -- which solves
             the same equation from the same two free functions, so its
             conversion needs nothing from anybody's converged answer and
             the one scalar MEQ hides, psi_ax, comes out of a fixed point
             that is exact in one step. converge.py is the acceptance and
             it is a MATRIX rather than a pair: a conversion error is a
             fixed offset, so only a difference that FALLS as each code is
             refined separately rules one out. M-158, and it is what says
             M-147's column is the reference's error rather than MEQ's.
             mesh/halfdisc.py is FB-6's geometry: a semicircle reaching the
             axis with rectangular coils MESHED TO, written against gmsh's
             python API. See `CLAUDE_FB.md`,
             *Meshing beyond MakeCartesian2D*. --symmetric meshes z >= 0 and
             REFLECTS it, which is what [solver] UpDownSymmetry needs and is
             reached from a file as [mesh.generate] Symmetric; --symmetry-check
             asks the question of a mesh this tool did not make. M-127
examples/    TOML run configurations
MEASUREMENTS.md  the published tables, under stable `M-nn` anchors that
             all five CLAUDE files point at. Measurements only; the reference
             tables stay inline where they are read.
CLAUDE_HDGGS.md      the core solver -- the equation, the discretisation and
             its conventions, Newton and NPC, the linear solves -- and then the
CLAUDE_FB.md         three campaigns: free boundary, the flux-surface
CLAUDE_INVERSION.md  inversion, and toroidal flow. CLAUDE.md is the index and
CLAUDE_FLOW.md       keeps the build, the commands, the traps, the testing
             stance and this layout.

             THE SPLIT IS BY CAMPAIGN AND NOT BY SIZE, so a claim that touches
             two of them is written in the file that owns the SUBJECT and
             cross-referenced from the other -- which is why a reference in
             any of them names the file when it leaves it. Three rules apply
             in all five: MEASUREMENTS.md's anchors are not renumbered; a
             hypothesis that was MEASURED AND FALSIFIED stays; and the history
             of the FILES does not.
refs/        Refs.md is tracked; the PDFs are gitignored, fetch by doi
             ( attic/ is GONE. It held the original
             von Hagenow / Lackner free-boundary code, unported and unbuilt,
             and it is superseded rather than pending: MEQ's free boundary is
             the exterior DtN coupling of FB-0..FB-5, which is a different
             method. Recover it with
             `git checkout 635aa3d -- attic/` if the Green's-function route is
             ever wanted; refs/Refs.md keeps the analysis of what it did. )
docs/        the Sphinx manual, published to Read the Docs. Built by
             `make -C docs html` or the `docs` CMake target, both under
             -W to match .readthedocs.yaml's fail_on_warning. Citations
             are sphinxcontrib-bibtex over docs/references.bib, which is
             written by hand from refs/Refs.md and is meant to agree with
             it. docs/manual/ is the pre-Sphinx LaTeX manual, moved intact
             and keeping its own Makefile.

             THE SPLIT IS DELIBERATE AND IT IS ALSO A STANDING RULE
             ABOUT WHAT MAY GO IN docs/:

               **Documentation reflects the state of the CODE, not a
               historical record. On release, docs/ must not carry the
               path not travelled or the mistaken ideas had along the
               way.**

             So a docs/ page says what MEQ does and, where a choice is
             exposed as an option, which way it came out and that it was
             measured.

             UNTIL A RELEASE IS TAGGED, THIS FILE FOLLOWS THE SAME RULE
             FOR THE SAME REASON -- everything is allowed to change
             instantly, so nothing here records that the CODE or this
             FILE once said something else. No "used to", no "an earlier
             version", no dated "now fixed". Write the present tense.

             WHAT IS NOT HISTORY, AND MUST STAY: a hypothesis that was
             MEASURED AND FALSIFIED. "An ordering does not fix a stiff
             source", "a line search makes every case worse", "blending
             the extrapolation does not change a sign" are findings about
             the PHYSICS AND THE CODE, and recording them is what stops
             them being re-derived -- most of the value in this file is
             exactly that. So is a live trap: "freegs4e's docstring says
             normalised flux and is wrong" belongs here and in docs/,
             because a reader meets it today. The line is between "this
             was once different" and "this does not work, here is the
             measurement".

             These pages are what a USER needs;
             this file is what a maintainer needs. Almost every choice in
             meq was settled by measurement, so the docs say which way it
             came out and that it was measured, and then tell the reader
             to measure it themselves wherever the choice is exposed as an
             option -- the trace solver, the assembly mode, the ordering,
             the globalisation, MKL_NUM_THREADS. The NUMBERS stay here,
             beside the tests that produce them, because a measurement in
             a manual goes stale silently.
```

Code style follows the sibling project MaNTA: **tab indentation**, Allman braces,
C++17. Naming, as decided for `../gffp` and carried here:

| | |
|---|---|
| Types — class, struct, enum, alias | `UpperCamelCase` |
| Enum values | `UpperCamelCase` |
| Functions and methods | `lowerCamelCase` |
| Variables, parameters, members | `lowerCamelCase` |
| **TOML configuration keys** | **`UpperCamelCase`** |

That last row is a deliberate mismatch with the C++ rule, not an oversight — it
comes from the same MaNTA convention set and applies to the key names in
`examples/*.toml` and to the string literals the parser looks them up by. Do not
reconcile the two.

**THE PROJECT IS `MEQ` IN PROSE AND `meq` IN CODE, AND THE SPLIT IS THE WHOLE
RULE.** It is an acronym and it is set the way its neighbours are — TOQ, EFIT,
CHEASE, ECOM — so any sentence about the code says **MEQ**. What stays lower
case is everything the compiler or the shell sees: the namespace `meq::`, the
directory `src/meq/`, the target `meq_core`, the umbrella header `meq/meq.hpp`,
and the command `meq config.toml`. `\meq` in `docs/manual/` expands to plain
`MEQ`.

**External names keep their author's capitalisation, whatever the table
says.** `TraceSolver::cuDSS` is spelled the way NVIDIA spells it and carries a
`// NOLINT(readability-identifier-naming)` saying so; `UMFPack` and `Pardiso`
are spelled as MFEM's wrappers spell them and happen to need no suppression.
The house rule governs MEQ's own identifiers and does not extend to renaming
other people's products. This is the same exemption `.clang-tidy` already
records for MFEM-imposed overrides like `Eval` and `Mult`, and it is written in
both places

**This is enforced**, by `.clang-tidy`'s `readability-identifier-naming` and a
ctest named `naming`. It runs over `MEQ_CORE_SOURCES_PRESENT`, which is now
**every** file in `src/meq` **except the seven a device compiler owns**.

**AND THAT EXCEPTION IS AN INSTRUMENT FAILING, NOT A RULE BEING RELAXED.**
clang-tidy reads each file's flags out of the compile database, and a CUDA entry
carries nvcc's: `--expt-extended-lambda`, `-forward-unknown-to-host-compiler`
and `--generate-code=arch=compute_75,...` are all `unknown argument` to clang,
after which `__clang_cuda_runtime_wrapper.h` cannot find
`texture_fetch_functions.h`. **Sixteen `clang-diagnostic-error`s and a check
that has run on nothing** — reported as a failure, which is the right report of
the wrong thing, and it took the `naming` ctest red for a session.
`CMakeLists.txt` fills `MEQ_CORE_SOURCES_DEVICE_LANGUAGE` only where a device
language was actually enabled, so **on a host build — which is what CI builds —
all seven are checked like everything else**. The convention still binds them; what is missing is a tool that can read
them here.

`src/meq` deliberately keeps MFEM out of `Profiles` and `Source` — plain `double`
arguments, no `mfem::Vector` — so both are unit-testable without the library, and
the `mfem::Coefficient` adapters live with the assembly that needs them.

## Git

Branch **`main`**. The authoritative remote is **`origin` =
`github:ianabel/meq.git`** (private; `github` is a `Host` alias in
`~/.ssh/config`).

**Tag `v0-legacy` is the pre-modernisation tree**, and it is where everything the
restructure deleted still lives. Reach for it before concluding something was lost.

**History is not rewritten**, deliberately: the deletions above are only
recoverable because it is not.

Commit messages end with the `Co-Authored-By` trailer.
