# Stage 7: the driver, warm starts, curved boundaries, and NetCDF output

**Stage 7 is done.** `meq config.toml` parses, builds the mesh and the source,
solves — with the adaptive loop and the curved boundary if asked — and writes
the equilibrium three times over, with exit codes 0/1/2/3. `CLAUDE.md` is the
operational record and carries every rate; `apps/meq.cpp` is the driver;
`examples/` holds the worked configurations.

**What is left here is what the code does not say for itself**: the findings the
staging produced, three acceptance criteria that were *wrong* and had to be
struck, and the two decisions that are still open. Section numbers are stable —
`src/meq/Output.hpp`, `src/meq/Sampler.hpp` and `src/meq/Config.hpp` cite §3 and
§4 by number, and `INVERSION-PLAN.md` and `CLAUDE.md` cite §3 and §5.

**One thing this plan was written around is worth keeping in view.** It recorded
the MFEM defect fixes as done on the strength of the other tree's measurements,
and said meq-side verification was the first job. **Insisting on that paid**: the
reconstruction fix reported here as complete was found from meq's side to be
silently wrong on **every element where `∂F/∂ψ` vanishes** — a per-element defect
no whole-domain norm can see — and needed a second upstream change to close. See
`CLAUDE.md`, *Post-processing is back*.

---

## 1. `setInitialGuess()`

**Built**: `setInitialGuess()`, `clearInitialGuess()` and `projectOntoTrace()`,
with three tests in `PedestalConvergence.cpp`.

### The obvious call is the wrong one

**Newton's unknown is the trace, not the potential**, so a guess expressed as
`ψ(r,z)` has to reach `M_h`. `GridFunction::ProjectCoefficient` loops over
`fes->GetNE()` — *volume* elements — so on a trace space it never reaches a face
dof, and `ProjectBdrCoefficient` iterates boundary faces only and so cannot seed
the interior trace, which is where nearly all the dofs are. Checked in
`fem/gridfunc.cpp:2529`. The pattern that works is the per-face loop
`Estimator.cpp` already uses on the same space —
`GetFaceElement` / `GetFaceVDofs` / `Project` — which is why `projectOntoTrace()`
is a helper and not three lines. `Estimator.cpp` records that `GetFaceVDofs` and
`GetFaceElement` agree on orientation, which is the detail that would otherwise
cost an afternoon.

The essential dofs are overwritten by the Dirichlet datum afterwards, **in that
order**, so the boundary condition always wins over the guess.

### A guess alone does not leave the trivial branch — the acceptance criterion was wrong

The brief asked for all four GS-2 §4.2–4.5 sources posed with *homogeneous* data
and a non-trivial guess converging to `‖ψ‖ > 0`. **Not achieved, and not
achievable with a guess alone.** Measured on §4.3: a bump of amplitude 0.05 or
0.20 converges in 4–5 iterations straight back to `ψ ~ 1e-16`, and 0.40 and above
fail in the element-local solves. **Newton goes to the root nearest its iterate
and `ψ ≡ 0` is a root**; Picard is carried *away* from zero because it evaluates
`F(ψ⁰) ≠ 0`, which is why the papers need only a guess and meq needs a guess
**and** globalisation. The non-homogeneous ramp stays the workaround.

§4's acceptance block wrote the same wrong criterion a second time, and it is
struck there for the same reason. **It should not have been written twice.**

The two criteria that were right and are asserted: the existing convergence
tables unchanged to six significant figures with a guess supplied — a starting
point must move the path and not the answer — and a guess that is already the
exact solution finishing in **one** Newton step.

### Still open: the flux block is not seeded

This plan hypothesised that seeding the **potential** block as well as the trace
would fix the pedestal's element-local failures. That hypothesis is moot under
`NonlinearOrdering::NPC`, which has no element-local non-linear solve at all —
but a sharper version of it survives. `prepare()` projects the guess onto the
potential and the trace and leaves the **flux block at zero**, because the guess
arrives as a bare `mfem::Coefficient` and cannot be differentiated. Under NPC
`q` is an unknown, so a warm start is inconsistent in exactly the row that
couples them and `‖r₀‖` goes *up*. The fix is to seed
`darcyFlux = −(1/r)∇̄ψ_guess` through `GradientGridFunctionCoefficient`, for the
`setInitialGuess( GridFunction const & )` overload that has something
differentiable to work with. **Not done, and it wants its own measurement.** See
`CLAUDE.md`, *A warm start no longer shows up in `‖r₀‖`*.

---

## 2. Configuring a curved boundary

**Built**: `src/meq/BoundaryShape.{hpp,cpp}` and the `[boundary.shape]` /
`[boundary.extension]` schema, with ten tests in `BoundaryShapeTests.cpp` and six
in `ConfigTests.cpp` — all MFEM-free, so **CI runs them**, which is a first for
anything geometric here. `examples/miller-curved.toml` is the worked example.

### Miller and MXH are one parametrisation, not two

`refs/MXH.pdf` eq (4) states the reduction explicitly: Turnbull–Miller is
recovered from MXH by keeping only `s_1` and `s_2`, with `s_1 = arcsin δ` and
`s_2 = −ζ` (`ζ` the squareness). **So Miller is MXH truncated at `N = 1`**, there
is one evaluator, and Miller is sugar over it — two implementations of the same
curve is how the `p`-versus-`F` inconsistency in Example 6 happened. Physical
readings, from MXH §2: `s_1` triangularity, `−s_2` squareness, `c_0` tilt, `c_1`
ovality; the `n`th harmonic makes an `(n−2)`-sided polygon.

**A transcription note worth keeping.** HDG-GS-1 Example 6 prints
`arcsin( δ sin t )` where Cerfon–Freidberg print `arcsin( δ ) sin t`. Miller
eq (34) is the third source and agrees with Cerfon–Freidberg, so **Example 6 is
the outlier**. `MillerDShape.hpp` keeps Example 6's form because it reproduces
that paper's table; the library does not. Measured, the MXH evaluator reproduces
`boundaryPointCerfonFreidberg()` **exactly** — worst difference `0` over 720
samples — and differs from Example 6's form by **6.25e-4**, matching the ≤6e-4
that fixture's header quotes. Both bounds are asserted, the second from *below*
as well, so a future edit that silently makes the two conventions one fails.

### The two coefficient arrays have different origins

`CosCoefficients` starts at `c_0` and `SinCoefficients` at `s_1`, because
`sin(0·θ)` is identically zero and an `s_0` slot would be a place to put a number
that silently does nothing. The parser says so and validates the lengths.

### Bisection needs star-shapedness, and this is the likeliest way to a wrong domain

The level set is the **radial gap** — polar angle about `(R0, Z0)`, curve
parameter by bisection, distance difference — not a signed distance. Larger than
the distance, and that is fine: `MarkLevelSetSubdomain` needs a sign at vertices
and `VertexConePath` needs a root along a ray, and neither needs a metric.

**But bisection is legitimate only if the curve is star-shaped about
`(R0, Z0)`.** For Miller it is, at moderate parameters. For MXH with a large tilt
`c_0`, a large `c_1`, or high harmonics it need not be, and the bisection would
then return one of several roots without complaining. So it is validated at load
over a fine `θ` sample and throws a `ConfigError` naming the `θ` where
monotonicity fails. **This is the single most likely way a plausible-looking
config produces a wrong domain**, and if a genuinely non-star-shaped `Γ` is ever
wanted the honest answer is a signed-distance field built from a sampled polygon
— a larger piece of work than anything in this plan.

The other three validations the config layer owes, all of them live: the
background box must contain `Γ` with a margin of a few `h` for the paths to reach
into; `RMin > 0` whenever a shape is given, since the `1/r` is not integrable
through the axis and a D-shape that crosses it is a mis-typed `R0`; and
`SearchLength` defaults **per path family**, `ExtensionConvergence.cpp` having
measured `VertexConePath` needing `6h` and `LevelSetPath` `12h` on the Solov'ev
level set.

---

## 3. Output

**Built**: three artefacts, two audiences — `.mesh` plus `_psi.gf` and
`_grad_psi.gf` for GLVis and meq's own exact restart; a ParaView collection at
the polynomial degree of the solve; and the `.nc` interchange file.
`tools/README.md` is the guide to which format goes with which reader, and
`src/meq/Output.{hpp,cpp}` and `src/meq/Sampler.{hpp,cpp}` are the code.

### `B_poloidal` is a relabelling of `q`, not a differentiation

From `refs/Miller.pdf` eq (1), `B = ∇φ × ∇ψ + f(ψ)∇φ` with `∇φ = ê_φ / R`:

```
B_R = -( 1/R ) dpsi/dZ          B_Z = +( 1/R ) dpsi/dR
```

and meq solves for `q = (1/r) ∇̄ψ` directly, so

```
B_R = -q_z          B_Z = +q_r
```

pointwise, with no derivative taken anywhere. **This is the entire argument for
the mixed method made concrete**: `q` converges at `k+1`, so `B_pol` does too.

**Do not trust those two signs on the strength of the derivation.** This
project's record on sign conventions in derivations is poor — the Solov'ev
source, the `τ` in eq (8e) and `DarcyForm`'s `−q` were each settled by
measurement after an argument said otherwise — and `flux()` already undoes
`DarcyForm`'s negation once, so there are two conventions in play and a sign
error here is invisible to every convergence rate. The acceptance was therefore a
measurement: `B_R`, `B_Z` from `q` against central differences of the *exact*
Solov'ev `ψ` on a sample cloud, at the finite-difference floor.

### The grid, and two decisions with reasons

A uniform tensor-product `(R, Z)` grid over the bounding box of `Γ`. It is the
EFIT-shaped thing every downstream tool expects and — decisively for §4 — it is a
format a warm start can interpolate from in `O(1)` per point with no mesh search
at all. The variables and attributes are in `Output.cpp` and `tools/README.md`.

* **`(Z, R)` with `R` fastest** is C row-major and matches MaNTA's `(t, x)`.
* **Points outside `Γ` get `_FillValue = NaN` *and* a zero in `inside`**: some
  tools honour the fill attribute and some do not, and the mask is what a reader
  can rely on. `byte extrapolated( Z, R )` was added later beside `inside`, for
  the reason `CLAUDE.md` records under *`B` gets the band too* — a *count* of
  continued nodes is not a mask, and nothing downstream could tell which nodes
  they were.

**Deferred, and still open: a flux-surface `(ψ_N, θ)` grid.** More useful
physically, and it needs contour tracing. That machinery now exists —
`INVERSION-PLAN.md` stages IN-0 to IN-4 — so the grid is stage **IN-6** rather
than an unbounded piece of work, and it is where this deferral now lives.

### Sampling without hitting the `FindPoints` trap

`mfem::Mesh::FindPoints` is `O(elements × points)` — a brute-force scan over
element centres. A 129×129 grid against a 20k-element mesh is 3.3e8 element tests
for a job that should be instant.

**Invert the loop.** For each element take its bounding box, convert that
directly to a grid index range — `O(1)`, because the grid is uniform — and test
only the handful of grid points inside it with
`ElementTransformation::TransformBack`. Total cost is linear in both. **This is
the one sampling job that does not want GSLIB** (§6): the grid is uniform, so
locating a point is index arithmetic and the search is already free; what remains
is the inverse element map, which MFEM supplies.

A grid point landing exactly on an inter-element face is ambiguous, `ψ_h` being
discontinuous; whichever element claims it first wins, and the jump is
`O(h^{k+1})` and converges away. **GSLIB has a principled answer where this has a
shrug** — `FindPointsGSLIB::SetL2AvgType`, arithmetic or harmonic averaging over
the elements meeting at the point, exactly for L2 fields that are multi-valued on
element boundaries. Both `ψ_h` and `q_h` are L2. Now that GSLIB is enabled, that
is the thing to reach for here, and it has not been.

---

## 4. Warm starts

### Three routes, and they are not the same problem

**Exact restart** — same mesh, same degree. Read the stored grid function and
hand it to `setInitialGuess()`. Bitwise resumable.

**Full-order interpolating restart** — a different mesh or degree, same code.
`meq::FieldTransfer` over `FindPointsGSLIB`, in `WarmStart.{hpp,cpp}`. The guess
carries the full `k+1` accuracy of the solve it came from. This is what needed
§6.

**Interchange restart** — another code entirely, or a mesh that no longer exists.
Read `psi(R,Z)` from the NetCDF and interpolate bilinearly, as an
`mfem::Coefficient`. `Config.hpp` names this one and says it is not the
interpolating route.

**An earlier draft had only the first and third, and that was a compromise
dressed up as a design.** The NetCDF route was chosen *because* `FindPoints` is
`O(elements × points)` and there was no other way to interpolate from a foreign
mesh cheaply. So the third route now keeps its job and stops pretending to be the
general one: it is the **interchange** format, and it is genuinely the best thing
for that — an outside code has to produce `ψ` on a rectangle and nothing else, no
MFEM, no mesh format, no agreement about element types.

### The comparison measurement, and the claim it corrected

The plan stated the two interpolating routes' ratio as `h^{k+1}` to `h^2`, which
**quietly compares two different `h`'s**: the mesh's and the *grid's*. Measured,
sampling a `k = 3, n = 8` solve (own error 6.00e-6) onto grids of 65, 129 and 257
nodes and reading it back bilinearly gives 1.18e-4, 3.36e-5, 1.01e-5 — **second
order in the GRID spacing**, rates 1.813 and 1.735, drifting below two only
because by 257 nodes the grid error is within a factor of two of the solve's own.

So at 257² the grid route is only **1.7×** worse than full order, not the
order-of-magnitude the plan implied. **The real statement is about scaling, not a
ratio**: the grid route's accuracy is a property of the **file**, the full-order
route's is a property of the **solve**. Refine the mesh or raise the degree and
the grid has to be refined quadratically to keep pace; the full-order route never
has to be. That is what justifies the dependency, and
`fullOrderCarriesMoreThanAStructuredGrid` asserts the scaling rather than the
ratio.

The same-mesh and refined-mesh restarts measured: `‖r₀‖` falls from 1.569e+01 to
1.803e-03, a factor of **8704**, Newton finishes in **1** step against 4 cold,
and the converged `L2` is 3.663930e-07 against 3.663932e-07 — seven figures.

### Points the old grid does not cover

A new domain can extend beyond the stored grid, and the stored grid is masked
outside `Γ`. Both give an interpolation query with no data. The guess is only a
guess, so the fallback need not be clever, but it must be *stated*: fall back to
the Dirichlet datum, count the misses, and report the count. **A restart that
silently found no data for 90% of the domain is a restart that quietly became a
cold start**, and the iteration count will not obviously say so. On the curved
path there really are such nodes, the computational domain growing as it refines
— measured, 48 and 360 of them on `miller-adaptive`.

### The criterion that was struck

~~The GS-2 sources warm-started from a coarse solve, converging where the
homogeneous cold start lands on `ψ ≡ 0`.~~ **Struck**, for §1's reason: a warm
start is a guess, and a coarse solve of the same homogeneous problem *is* the
trivial root, so warming from it cannot change which root is found. This bullet
contradicted §1's acceptance block on the same page.

---

## 5. The driver

Single binary, `meq config.toml`. No subcommands, and specifically not MFEM's
`OptionsParser`, which wants to own argument parsing for the whole program.

### What the driver owes the user

* **Exit codes that mean something.** 0 solved; 1 configuration error; 2 Newton
  did not converge; 3 output could not be written. A shell script driving
  parameter scans is a first-class caller. `MFEM_USE_EXCEPTIONS` is what makes 2
  reachable rather than a SIGABRT.
* **A residual history on stdout** — iteration, `‖r‖`, `‖r‖/‖r₀‖`, observed
  order. It is the diagnostic that separates a wrong Jacobian from a hard
  problem, and it costs nothing to print.
* ~~**`MKL_THREADING_LAYER=GNU`.**~~ **DROPPED, 2026-09-01.** This asked the
  driver to detect MKL and set the variable or refuse to start, on the grounds
  that the failure is silent wrong numbers. It was implemented as a warning and
  then deleted: meq builds against its own SuiteSparse and never loads the
  `libmkl_rt` dispatcher the variable configured, so it is inert — measured
  bit-identical across three value-asserting suites — and **most builds link no
  MKL at all**, so the warning told the majority of users to set a variable
  naming a library they do not have. Choosing a threading layer is CMake's job,
  since it can see what the build links.

### Three things the plan did not anticipate

* **The curved path needed `meq::AdaptiveDomain`, not the driver's one-shot
  `buildSubdomain()`.** Refining `D_h` alone leaves `Γ_h` where it is while
  `h_loc` halves, so `dist/h_loc` doubles every cycle and the transfer quietly
  leaves the regime it is analysed in. Both constructions are in the driver and
  the comment at the branch says why, rather than tidying one away.
* **`ψ*` was a hard prerequisite, and retiring `postProcess()`'s refusal was
  listed in `ROADMAP.md` as an independent item.** It is not: four of eq. (20)'s
  five terms are built on `ψ*` and the driver's path is *always* semi-linear.
* **And the refusal could not be retired outright.** Measuring it found the MFEM
  fix worked where `∂F/∂ψ ≠ 0` and returned a different function on **any**
  element where it vanishes. The hole is closed upstream as *"The postprocessing
  closes on the element average, always"*; `Potential::Raw` survives as the
  *measurement* of the one order it costs, not as anything the driver runs.

### The transferred datum is rebuilt every cycle, and that is not fussiness

It lifts the *solved* flux along the transfer paths, so it goes stale the moment
the solver is solved again — which in an adaptive loop is every cycle, on a new
mesh with a new path family. Hoisting it out would compare `ψ*` against the
previous cycle's boundary condition.

### What the driver still refuses

`[boundary] Type = "exact"`, which needs a closed form `meq::Source` does not
carry. It exits 1 with an explanation rather than approximating. **That is the
only thing it refuses.**

---

## 6. GSLIB

**Taken.** `MFEM_USE_GSLIB` is `YES` in the install, gslib v1.0.9 is built
alongside at `../mfem/gslib`, and `meq::FieldTransfer` is the consumer.

**The argument for taking it, which generalises.** The dependency lands on the
**MFEM tree, not on meq**: meq's `CMakeLists.txt` is untouched,
`find_package(MFEM)` is unchanged, and the MFEM-free half of the library is
unaffected. CI already cannot build MFEM at all, since the branch is unpublished,
so this adds nothing to what CI cannot do. That is a different calculation from
adding a dependency to meq itself, and it is the same reasoning `CLAUDE.md`
records for why toml11 is a submodule and MFEM is not.

It is packaged by no distribution — checked: nothing in Debian, Ubuntu, Fedora
43–45, EPEL 8–10.4 or Rawhide — but it is pure C with no dependencies of its own
and builds as a sibling of the MFEM tree in three lines.

**The caveat is half-measured, and the missing half is the one open item here.**
meq meshes are **triangles** and gslib's `findpts` is tensor-product; MFEM
handles that by splitting simplices into quads internally (`mesh_split`,
`ir_split`, `split_element_map` in `fem/gslib.hpp`). It works, and it is a code
path a quad-based user never exercises. `WarmStartConvergence` does run it on
triangles and does measure the result — but **against the exact solution, not
against §3's sampler**. The two are never compared to each other, on the same
points, on the same mesh. That check is what the inverted loop was kept for, and
**two independent samplers agreeing is worth more than either one being
plausible**.

---

## Risks that are still live

**Star-shapedness.** §2. The validation turns a wrong domain into an error
message, which is the most that can be done without a different level-set
construction.

**`FindPoints`.** §3's inverted loop is a complexity claim, and complexity claims
rot. Time it rather than asserting it in a comment.

**GSLIB's simplex path.** §6. The fallback if it is slow or fragile is §4's third
route, and nothing is lost but the order of the warm start.
