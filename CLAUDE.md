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
  order as ψ and `∇̄ψ = r q`, so a node `p` outside the mesh is reached from its
  foot `x₀` on `Γ_h` as `ψ(x₀) + r₀ q(x₀)·(p − x₀)` — **nothing is ever
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
  same thing. **The route to better is known and was not taken**: `div q = −F/r`
  is the equation being solved and `∂_r q_z − ∂_z q_r = −q_z/r` follows from
  `r q = ∇̄ψ`, which pin two of `∇q`'s four entries exactly — but they leave the
  symmetric traceless part still differentiated, so they buy *structure* rather
  than an order, at the cost of plumbing the source into `GridSampler`.

  **AND THE COST IS MEASURED, NOT ESTIMATED.** Wiring the rotating output gave
  a free controlled experiment: `n_s` is algebraic in `(r, ψ)`, so the same
  closed form can be evaluated both ways over the same 356 band nodes. Against
  the node's own `r` it is **5.13e-07** wrong; against the foot's `r`, which is
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
  the output will meet it first: `postProcess()` is **0.62 s, 67% of the output at
  the default grid and 8% of the whole DIII-D run**, and it does not move with the
  grid because it has nothing to do with it. It is `DarcyForm::Reconstruct()`
  re-assembling four integrators at the enriched order per element, which is the
  price of reporting `ψ*` rather than `ψ_h`.
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
warning, since a wall-hugging annulus is a real thing to look at. **A source that
does not vanish on the symmetry axis** refuses outright — `F/r` is `μ₀ j_φ`, so
that is an infinite current density on `r = 0`. The source one is **ordered
first**, because a bad `ψ_ax` is its consequence and the `ψ_ax` message's advice
(look at the guess, look at the mesh) is wrong when the pole is the cause. Both
exit 1 and write nothing; `docs/running.rst`'s gloss on code 1 is amended, since
it said *"nothing was attempted"*.

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

Group 2 needs a partial-assembly rewrite and is not built. **So a config-file
cuDSS today is exactly the group-4-alone case that plan says not to do**: a
Device would be configured, every Vector in the process would allocate through
it, 58–70% of an NPC step would still run on the host, and each Newton
iteration would pay a round trip for the one part that moved. The key opens when
the offload work lands, not before.

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
was written and could not fail until a Device existed: `buildForms()` makes
three long-lived `MakeRef` aliases — `darcyFlux`, `potentialGf`, `traceGf` —
onto the blocks of one `BlockVector solution`, and every one of them is a
separate alias registration that never hears about a device write. `solve()`
therefore calls `solution.SyncToBlocks()` and a `SyncMemory()` per grid function
before the first host read, and they are inert with no Device configured. This
is the **third** time a library update or a new configuration has turned a
dormant MEQ contract into a failure in somebody else's code; see the `SubMesh`
parent and `SourceIntegrator`'s shared scratch under *Traps*.

**cuDSS is a bystander and always was.** UMFPack and PARDISO never touch the
device and failed identically. With the four sites fixed, `--device cuda` and
`--device debug` both reproduce the CPU answer to **every printed digit**.

**AND THE REFUSAL OF `TraceSolver = "cudss"` STANDS, ON THE TRADE ARGUMENT
ALONE.** The second and stronger reason that used to sit here — that the solve
did not survive a Device — is gone, and what is left is the group-4-alone
argument above, which is unchanged: the integrators are 46–53% of a step and
have no kernels. The flag stays in the harness so that the day they get some,
this costs one command.

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

### Which integrators need device offload, and why kernels alone are not enough

**THE INVENTORY, TAKEN FROM THE INSTALLED HEADERS RATHER THAN REMEMBERED.** Every
integrator MEQ installs, and whether a device entry point exists for it today:

| integrator | where MEQ installs it | device path |
|---|---|---|
| `mfem::VectorMassIntegrator( radius )` | flux mass, domain — `GradShafranov.cpp:3077` | **`AssemblePA`, `AssembleMF`, `AssembleDiagonalPA`, `AddMultPA`** — and separately `HDGElementMassBatched()` under `AssemblyMode::Batched` |
| `mfem::VectorDivergenceIntegrator` | flux divergence, domain — `:3176` | **`AssemblePA`** |
| `mfem::TransposeIntegrator( mfem::DGNormalTraceIntegrator )` | flux divergence, interior + boundary faces — `:3177`, `:3179` | the wrapper has PA and the wrapped integrator has none — **and neither matters**, because under hybridization these are never assembled. They are MARKERS; see the measured trap under *A trap in that table* |
| `mfem::NormalTraceJumpIntegrator` | the hybridization constraint — `:3186` | `AssembleEAInteriorFaces` **only**: element assembly, interior faces, no `AddMultPA` and nothing for boundary faces |
| **`mfem::HDGDiffusionIntegrator`** | potential mass, interior **and** boundary faces — `:3109`, `:3151`, `:3152` | **none** |
| **`meq::SourceIntegrator`** | potential mass NON-LINEAR form, domain — `:3140`, `:3150` | **none, and it is the one MEQ owns** |
| **`mfem::HDGExtensionIntegrator`** | flux mass, boundary faces, curved and free-boundary only — `:3098` | **none** |
| `mfem::DomainLFIntegrator`, `mfem::VectorBoundaryFluxLFIntegrator` | the two right-hand sides — `:3465`, `:3414` | none, and it does not matter: once per mesh, not once per step |

**TWO STRUCTURAL FACTS MATTER MORE THAN THE TABLE, AND EITHER ONE ALONE MAKES
"ADD KERNELS" THE WRONG PLAN.**

**The hybridized path never reaches a PA entry point.** `darcyform.hpp:169` says
the hybridized assembly goes through `BilinearForm::ComputeElementMatrix()` —
the dense, per-element, host route. So the three integrators that *do* have PA
kernels are not using them in MEQ, and kernels added to the other three would
not be reached either. What is missing is a **route**, which is upstream's group
2 partial-assembly rewrite; kernels are necessary and nowhere near sufficient.

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
**`AssemblyMode::Batched` is a 15 to 27 per cent LOSS** on the host, and
**`LocalFactorMode::Batched` is a trade that comes out flat** — it buys a batched
local factorisation and pays the condensation cache, since
`CanCacheCondensation()` refuses outright under it.
→ **[M-99](MEASUREMENTS.md#m-99)** and **[M-101](MEASUREMENTS.md#m-101)**.

**THE ORDERED LIST, WEIGHTED BY THE LEG PROFILE RATHER THAN BY COUNT.** From
**[M-80](MEASUREMENTS.md#m-80)**:

1. **`meq::SourceIntegrator`** — the `F(r, z, ψ)` domain term, evaluated per
   quadrature point per element per residual **and** per Jacobian. It is the
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
3. **A device route through `DarcyHybridization`** that reaches those kernels
   instead of `ComputeElementMatrix()`. Without this, 1 and 2 are unreachable.
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
worth *timing* here afterwards, since this machine's consumer FP64 runs at 1/32
to 1/64 of its FP32 rate where a datacentre part runs it at about 1/2.

**The transferable part**: exposing a knob makes reachable, by somebody who has
not read the library's documentation, every precondition that documentation
records. The Device requirement had been written down in `GradShafranov.hpp`
all along.

**THE INTERPOLATING WARM START IS WIRED.** `MFEM_USE_GSLIB` is `YES`,
`meq::FieldTransfer` is `FindPointsGSLIB` in `WarmStart.{hpp,cpp}`,
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

Beyond the port, three campaigns have their own plans, their own staging and
their own file:

| | | |
|---|---|---|
| toroidal flow | FL-0 to FL-8, **done** | `docs/rotation.rst` is the derivation, and `CLAUDE_FLOW.md` the record |
| solution inversion | IN-A to IN-P, **every stage done** | `INVERSION-PLAN.md`, and `CLAUDE_INVERSION.md` |
| free boundary | FB-A to FB-7 **done or answered**, FB-6 **met**, and §10's XP-0 to XP-4 **done** — the X-point is two unknowns of the bordered Newton, driveable from a file, and MEQ reproduces `freegs4e`'s DIVERTED equilibrium | `FREE-BOUNDARY-PLAN.md`, and `CLAUDE_FB.md` |

The two plan files that had nothing left in them are gone, converted to `docs/`:
`DRIVER-PLAN.md` was stage 7 and `FLOW-PLAN.md` was the flow campaign. Git has
them. What is left open across all three is short and `ROADMAP.md` is the
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

37/37 in every configuration when that table was taken, and **52 of 53 today**,
the one red being `PlasmaEdgeConvergence` and deliberate —
the count moves as cases are added, so read the table's ratios rather than its
absolute seconds. Nothing in the suite depends on a thread count, which is the
correctness half.

**WHERE A CASE LIKE `theTwoBorderSolveReportsATrueMagneticAxis` GOES RED, SUSPECT
THE FIXTURE BEFORE THE SOLVER.** `examples/limited-tokamak.toml` has **every** ingredient the case
has — a domain reaching `r = 0`, `[boundary.limiter]`, `[boundary.exterior]` and
`ConfineToPlasma` — and converges in **11 Newton steps** to a `Ψ` of 1.0000,
reproducing freegs4e to 1.3e-04. What it had that the fixture did not is **coils
and a prescribed current**. The fixture now has both, plus the confinement, and
it is **green at five limiter radii in a row**, 1.08 to 1.18, with the coil-free
control converging to the wrong topology at three of them; see `CLAUDE_FB.md`,
*§7.12b's fixture was the defect*, and `FREE-BOUNDARY-PLAN.md` §11.3 for the diagnosis and §11.7 for
the repair.

**AND "THE FIXTURE CANNOT BE REPAIRED IN PLACE" IS THE WRONG CONCLUSION TO DRAW
FROM IT BY ONE INGREDIENT.** `ConfineToPlasma` alone fails 4/4 and clamped
profiles alone fail 4/4; both measurements stand, and what they show is that
confinement is **one third** of the repair rather than that there is none — the
same run records "clamped WITH a prescribed current converges in 22 steps".
**A measurement that rules two things out is not a measurement that rules
everything out.**

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

**What CUDA does NOT buy MEQ today**, so nobody reads the flag as a capability:
`fem/darcy/` still has no partial-assembly kernels, MEQ uses no `mfem::forall`,
and MEQ's own sources still compile with `g++` rather than `nvcc`. What it
enables is `BatchedLinAlg`'s `gpu_blas`/`magma` backends and `CuDSSSolver`. It
is a prerequisite, not a speedup.

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

MEQ needs work from **nine MFEM branches**, and for as long as this section
existed it named four. See the warning under the containment loop for how the
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
| **`gf-interp-hdg-dev`** | **CCSZ INTERPOLATORY HDG** — `fem/darcy/reaction_hdg.{hpp,cpp}` with `HDGInterpolatoryReactionIntegrator`, `NodalReactionFunction` and a quadrature control, plus `HDGPostprocessBlocks` in `postprocess_hdg.*` and `DarcyHybridization::Bg_data`. **15,024 insertions across `fem/darcy/`**, verified in `libmfem.a` by `nm`. MEQ does not call any of it yet; `INTERPOLATORY-HDG-PLAN.md` is the plan | merge 8 |

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
# ALWAYS finish with this, and believe it over the table above:
for b in gf-hdg-subdomains-dev direct-solver-symbolic-reuse \
         gf-hdg-linearise-first gf-hdg-dev pardiso-multi-rhs \
         pardiso-device-host-sync sundials-ida-integration \
         gf-hdg-p-adaptivity gf-interp-hdg-dev; do
  git merge-base --is-ancestor hdgdev/$b HEAD && echo "$b contained" \
                                              || echo "$b NOT CONTAINED"
done
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
whole point, and it is what found this:

```sh
git -C ../mfem/mfem-src log --merges --format=%s meq-integration \
  | grep -oE "hdgdev/[a-z0-9-]+" | sort -u
```

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

**Check the install too, and it is one command:**

```sh
diff <( git -C ../mfem/mfem-src show meq-integration:fem/darcy/darcyhybridization.hpp ) \
     ../mfem/install/include/mfem/fem/darcy/darcyhybridization.hpp \
  && echo "install is current"
```

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
| **`AssemblyMode::Threaded` must not be refused on `CopyLinearGradBlocks()` declining** | Upstream is adding an abort for a real hazard — 67 of 99 MFEM integrators with scratch are unguarded — keyed on a predicate that catches MEQ for an unrelated reason. MEQ is `LocalOpType::PotNL`, so the cache **always** declines, and on that branch `ConstructGrad()` skips the flux block and evaluates only `meq::SourceIntegrator`, which is MEQ's own and reentrant. **An MFEM update carrying that abort stops every threaded MEQ run.** Filed as `../mfem-hdg-dev/doc/HDG-THREADED-REFUSAL-FROM-MEQ.md`; `CLAUDE_HDGGS.md`, *Threading, measured* |
| **Threading `NPCReduce()` and `NPCRecover()`** | The two element loops either side of the trace solve are the only ones in `DarcyHybridization` with no `omp parallel`, and upstream's own doxygen declines them at *"under 6% of the step"* — measured fixed boundary, one right-hand side. MEQ's bordered step applies `J^-1` to `N + 4` columns, so the traversal is `O( elements × columns )` against integrator loops that are `O( elements )` and already threaded: **0.212 s at one thread and 0.212 s at eight, 30.6% of a threaded step and the largest leg in it.** `NPCRecover` writes only the calling element's own L2 dofs and needs neither colouring nor atomics. Filed as `../mfem-hdg-dev/doc/HDG-NPC-TRAVERSAL-FROM-MEQ.md`, → **[M-126](MEASUREMENTS.md#m-126)** |

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
| `(r q_h, v)` | `HDGDomainIntegratorGS` | flux mass: `VectorMassIntegrator` with an `r` coefficient |
| `(ψ_h, ∇̄·v)`, `(q_h, ∇̄w)` | `HDGDomainIntegratorGS` | flux divergence: `VectorDivergenceIntegrator` |
| `⟨ψ̂_h, v·n⟩` | `HDGFaceIntegratorGS` | the transpose of `NormalTraceJumpIntegrator` — **not** the face integrators on `B`, see below |
| `⟨τ(ψ_h − ψ̂_h), w⟩` | `HDGFaceIntegratorGS` | potential mass: `HDGDiffusionIntegrator` |
| `⟨q̂_h·n, μ⟩ = 0` | `AssembleSC` | `EnableHybridization(trace_space, new NormalTraceJumpIntegrator(), ess_list)` |
| condense / reconstruct | `AssembleSC` + `Reconstruct` | `FormLinearSystem` / `RecoverFEMSolution` |
| `φ_h` on `Γ_h` (stage 5) | — | `HDGExtensionIntegrator` on the **flux mass form**, `C = r`, sign `+1` |

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
a trivial branch — and MEQ falls into it.** `F(r, 0) = 0` for eqs (24), (25),
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

**If a waiter is used, check it actually fired.**

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
the `r² log r` corner term of a right angle — the solution sits in `H^{3−ε}` and
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

## Layout

```
src/meq/     the library. Config, Profiles, Source, SourceFactory,
             GradShafranov, BoundaryShape, Estimator, Field, Sampler,
             WarmStart, Output -- all ported and all under the naming check.
             CriticalPoints, FluxSurfaces, SurfaceAverage, SurfaceFit and
             Zernike are the solution-inversion work and were written here
             rather than ported; see that section.
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
**every** file in `src/meq`

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
