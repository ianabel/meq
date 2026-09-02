# Stage 7: the driver, warm starts, curved boundaries, and NetCDF output

A plan, not an implementation. Nothing here has been built.

meq closes stage 6 with every claim about rates measured and the solver
unusable by anyone who is not running `ctest`. `apps/meq.cpp` is still the
*legacy* driver — it solves the vacuum coil field and calls `WriteOutputMFEM()`
on a class the port deleted. Nothing in `src/` writes a file. This plan is the
work that turns a measured library into a program.

**Two restrictions this plan was written around have since been removed.** All
four defects of `../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md` are fixed on
`gf-hdg-subdomains-dev`, and one of them was worse than reported from outside:
the local problem behind `Reconstruct()` had no face constraint either, so the
matrix was singular and was factored and solved anyway. What follows for this
plan:

* **`ψ*` works through Newton**, so §3 may write it unconditionally and §5's
  adaptive loop is no longer linear-only. Both sections are marked below.
* **`φ_h` is reachable** — `TransferredDatumCoefficient` is the datum itself,
  with `PathLiftCoefficient` for the solution-dependent half. So meq's
  `setTransferredBoundary()` can stop excluding those faces and `η₅` can be
  built as eq. (20) writes it. That is a repair of an omission, not a new
  feature, and it is stage-6 work rather than stage-7 — but it belongs on the
  list.

~~**None of it is verified from meq's side yet**, which is the first thing to do
once meq builds against the new library. The claims above are the MFEM tree's
measurements, not this project's.~~ **BOTH ARE VERIFIED FROM MEQ'S SIDE NOW**,
the first in stage 6 and the second on **2026-09-02**.

* `postProcess()`'s refusal is **retired**, and the comment where it stood
  records what it used to guard. `ψ*` is measured through Newton at `k+2` —
  3.05, 4.05, 5.03 — by `NewtonConvergence.cpp`'s
  `thePostProcessedPotentialSurvivesNewton`.
* `η₅` is **rebuilt on the datum**. `ResidualEstimator::setTransferredBoundary()`
  takes `φ_h` as a second argument, `GradShafranovSolver::transferredDatum()`
  builds it over `mfem::TransferredDatumCoefficient`, and the two measurements
  are `ExtensionConvergence.cpp`'s `theTransferredDatumRestoresEtaFive` and
  `theTransferredDatumReproducesTheImposedCondition`.

**And insisting on meq-side verification paid, which is why the sentence above
is struck rather than deleted.** The reconstruction fix reported here as done
was measured from meq's side and found to be silently wrong on **every element
where `∂F/∂ψ` vanishes** — a per-element defect no whole-domain norm can see. It
took a second upstream change, *"The postprocessing closes on the element
average, always"*, to close. See `CLAUDE.md`, *Post-processing is back*.

## The pieces, and why this order

| | | ends at |
|---|---|---|
| 7a | `setInitialGuess()` | a source that vanishes at `ψ = 0` converging to something that is not zero |
| 7b | Boundary shapes: Miller and MXH | **done** — `BoundaryShape`, exact against the fixture's Cerfon–Freidberg form |
| 7c | Output: mesh, grid functions, NetCDF | `B` from `q` matching a finite difference of the exact `ψ` |
| 7d | Warm start | a restart from 7c's file cutting Newton to one or two steps |
| 7e | The driver | `meq examples/*.toml` writing files, end to end, as a ctest — **done**, and the adaptive loop and the curved boundary with it |
| — | GSLIB, in the SUNDIALS rebuild (§6) | a full-order warm start measurably beating the interpolated one |

7a before 7e because a driver without it is *actively misleading*: handed a
physical profile it converges in zero iterations and writes a file full of
zeros, which looks exactly like success. 7c before 7d because the interchange
warm-start input format **is** the output format — see §4. 7b is independent and could be
done at any point; it is placed third because 7e needs it and 7c does not.

Each stage ends at a measurement, per the standing rule in `CLAUDE.md`.

---

## 1. `setInitialGuess()`

### The problem it solves

`F(r, 0) = 0` for GS-2 eqs (24), (25), (26) and (27) — every source in §§4.2–4.5.
With homogeneous Dirichlet data `ψ ≡ 0` therefore *solves* the discrete problem,
and Newton starts from the Dirichlet data and stops there in zero iterations with
an identically zero residual. GS-1's Algorithm 2 opens
`ψ⁰ ; // Non-trivial initial guess` for exactly this reason.

Today no caller can supply one: `prepare()` does `solution = 0.0` at
`GradShafranov.cpp:459` and `solve()` calls `prepare()`.

### The design point that is easy to get wrong

**Newton's unknown is the trace, not the potential.** `solve()` runs
`newton.Mult( traceB, traceX )` on the condensed system; the volume unknowns are
reconstructed afterwards by `RecoverFEMSolution`. So an initial guess expressed
as `ψ(r,z)` cannot simply be written into the potential block — it has to reach
`M_h`.

The route is an `L²(e)` projection of the guess onto each face:

```
psihat_0 |_e  =  P_{M_h} ( psi_guess |_e )
```

**and there is no library call that does it.** `GridFunction::ProjectCoefficient`
loops over `fes->GetNE()` — *volume* elements — so on a trace space it never
reaches a face dof. `ProjectBdrCoefficient`, which `prepare()` already uses,
iterates boundary faces only and so cannot seed the interior trace, which is
where nearly all the dofs are. Checked in `fem/gridfunc.cpp:2529`; this is the
one place in this plan where the obvious call is the wrong one.

The pattern that works is the one `Estimator.cpp:325` already uses for the same
space:

```
for each face f:
    traceFe    = traceFes.GetFaceElement( f )     // DG_Interface, VALUE map type
    traceFes.GetFaceVDofs( f, vdofs )
    traceFe->Project( guessOnThatFace, faceTransformation, vals )
    traceGf.SetSubVector( vdofs, vals )
```

with the guess restricted to the face through
`Mesh::GetFaceElementTransformations`. `Estimator.cpp:317` records that
`GetFaceVDofs` and `GetFaceElement` agree on orientation, which is the detail
that would otherwise cost an afternoon. So this is a **small helper, not three
lines** — call it `projectOntoTrace()`, put it beside the estimator's copy of the
same loop, and consider whether the two should share.

The essential dofs are then overwritten by the Dirichlet datum, **in that
order**, so the boundary condition always wins over the guess.

### Proposed interface

```cpp
/// Start Newton from psiGuess rather than from zero. Ignored on the linear
/// path, where the solve is direct. The guess is projected onto the trace
/// space; the Dirichlet datum is applied after it and therefore overrides it
/// on essential dofs. Borrowed, and must outlive the next solve().
void setInitialGuess( mfem::Coefficient &psiGuess );

/// Start from a potential computed elsewhere -- a previous solve on this or
/// another mesh. See Warm starts.
void setInitialGuess( mfem::GridFunction const &psiGuess );

/// Forget the guess; the next solve() starts from the Dirichlet data alone.
void clearInitialGuess();
```

`prepare()` changes from `solution = 0.0` to: zero, project the guess if there is
one, then project the Dirichlet datum on `fittedMarker` — three lines in
`prepare()`, over the helper above.

### The second thing to seed, and why it is a hypothesis rather than a claim

`CLAUDE.md` records that the pressure pedestal fails at `k = 1` for `h ≥ 0.05`
with MFEM's *element-local* nonlinear solves giving up —
`el: N not convered in 100 iters`. Those local solves have their own Newton
iteration and their own starting point, which comes from the volume blocks of the
block vector that `prepare()` currently zeroes.

So seeding the **potential block** as well as the trace may fix a failure that
looks like a globalisation problem and is not one. That is worth measuring before
reaching for SUNDIALS, because it is a three-line change against a rebuild of
another tree. **It is a hypothesis.** It might do nothing. Measure it on the
pedestal at `k = 1, h = 0.05`, which is the case that fails today.

### Acceptance

**Done, and one criterion below was wrong.** `setInitialGuess()`,
`clearInitialGuess()` and `projectOntoTrace()` are in, with three tests in
`PedestalConvergence.cpp`.

* ~~All four GS-2 §4.2–4.5 sources posed with *homogeneous* data and a
  non-trivial guess, converging to `‖ψ‖ > 0`.~~ **Not achieved, and not
  achievable with a guess alone.** Measured on §4.3: a bump of amplitude 0.05
  or 0.20 converges in 4–5 iterations straight back to `ψ ~ 1e-16`, and 0.40
  and above fail in the element-local solves. Newton goes to the root nearest
  its iterate and `ψ ≡ 0` is a root; Picard is carried *away* from zero because
  it evaluates `F(ψ⁰) ≠ 0`, which is why the papers need only a guess and meq
  needs a guess **and** globalisation. The ramp stays the workaround.
* The existing convergence tables unchanged to six significant figures with a
  guess supplied. A starting point must not move the answer — only the path. This
  is the same invariance the `+5%` Jacobian experiment relied on, and it is the
  check that says the guess is entering as a guess and not as data.
* A guess that is already the exact solution finishing in **one** Newton step.

---

## 2. Configuring a curved boundary

### What exists

The extension machinery is stage-5 work and is done. `ExtensionConvergence.cpp`
drives it:

```
levelSet(x)                                  a function, negative inside
  -> mfem::MarkLevelSetSubdomain( background, levelSet, 0.0, ... )
  -> mfem::SubMesh                           D_h, the computational domain
  -> mfem::VertexConePath( *sub, gammaH, levelSet, 6*h )
  -> solver.setExtension( path, gammaHMarker )
```

So the whole configuration problem reduces to one question: **how does a TOML
file produce `levelSet`?** Everything downstream is already wired.

### Miller and MXH are one parametrisation, not two

`refs/Miller.pdf` eq (34):

```
R_s = R0 + r cos[ theta + (arcsin delta) sin theta ]
Z_s = kappa r sin theta
```

`refs/MXH.pdf` eqs (1)–(3):

```
R( theta ) = R0 + r cos( theta_R )
Z( theta ) = Z0 + kappa r sin( theta )
theta_R    = theta + c_0 + sum_{n=1}^{N} [ c_n cos( n theta ) + s_n sin( n theta ) ]
```

and MXH eq (4) states the reduction explicitly: Turnbull–Miller is recovered by
keeping only `s_1` and `s_2`, with `s_1 = arcsin delta` and `s_2 = -zeta`
(`zeta` the squareness). **So Miller is MXH truncated at `N = 1`.** There should
be one evaluator and Miller should be sugar over it — two implementations of the
same curve is how the `p`-versus-`F` inconsistency in Example 6 happened.

Physical readings, from MXH §2: `s_1` triangularity, `-s_2` squareness, `c_0`
tilt, `c_1` ovality; the `n`th harmonic makes an `(n-2)`-sided polygon.

**A transcription note worth keeping.** `MillerDShape.hpp` records that HDG-GS-1
Example 6 prints `arcsin( delta sin t )` where Cerfon–Freidberg print
`arcsin( delta ) sin t`, and calls the difference harmless (≤ 6e-4 in `r`).
Miller eq (34) is now the third source and agrees with Cerfon–Freidberg, so
**Example 6 is the outlier and the `arcsin(delta) * sin(theta)` form is the
right one to implement.** `MillerDShape.hpp` keeps Example 6's form because it
reproduces that paper's table; the library should not.

### Proposed schema

```toml
[boundary]
Type = "zero"              # unchanged: zero | exact

[boundary.shape]
Type = "mxh"               # none | miller | mxh | levelset | polygon
R0 = 1.0                   # geometric centre, metres
Z0 = 0.0
MinorRadius = 0.32         # r
Elongation = 1.7           # kappa

# Type = "miller": physical parameters, converted to MXH on load
Triangularity = 0.33       # delta   -> s_1 = arcsin( delta )
Squareness = 0.0           # zeta    -> s_2 = -zeta          (optional)

# Type = "mxh": the harmonics directly
CosCoefficients = [ 0.0, -0.03, 0.01 ]    # c_0, c_1, ... c_N   (c_0 is the tilt)
SinCoefficients = [ 0.3576, 0.02 ]        # s_1, s_2, ... s_N   (there is no s_0)

[boundary.extension]
PathFamily = "vertexcone"  # vertexcone | levelset
SearchLength = 6.0         # in units of h
ConeRays = 64              # vertexcone only
```

**The two arrays have different origins and that is a trap**, so the parser must
say so: `CosCoefficients` starts at `c_0`, `SinCoefficients` at `s_1`, because
`sin(0*theta)` is identically zero and an `s_0` slot would be a place to put a
number that silently does nothing. Validate that `CosCoefficients.size()` is
`SinCoefficients.size() + 1`, or accept either and document the padding.

`Type = "none"` keeps today's behaviour: solve on the mesh boundary directly,
fitted, no extension. That is what every current convergence test does and it
must stay reachable.

### From parametrisation to level set

`MillerDShape::levelSet()` already proves the pattern and the library should lift
it: for a point `(r,z)`, take the polar angle about `(R0, Z0)`, find the curve
parameter `t` with that polar angle by bisection, and return

```
levelSet( r, z ) = |(r,z) - (R0,Z0)| - |boundaryPoint(t) - (R0,Z0)|
```

This is the **radial gap**, not a signed distance — larger than the distance, and
that is fine: `MarkLevelSetSubdomain` needs a sign at vertices and
`VertexConePath` needs a root along a ray. Neither needs a metric.

**Bisection is legitimate only if the curve is star-shaped about `(R0, Z0)`** —
the polar angle must be strictly increasing in `theta`. For Miller it is, for
moderate parameters. For MXH with a large tilt `c_0`, a large `c_1`, or high
harmonics it need not be, and the bisection would then return one of several
roots without complaining.

**So validate it at load, over a fine `theta` sample, and throw a `ConfigError`
naming the `theta` where monotonicity fails.** `MillerDShape` checks this in a
test rather than at construction; a configured shape comes from a user and must
check itself. This is the single most likely way a plausible-looking config
produces a wrong domain.

### The other validation the config layer owes

* **The background box must contain `Γ` with room to spare.** The extension needs
  `Ω^h ⊂ Ω` *and* a layer of background elements outside `Γ_h` for the paths to
  reach into. Check the shape's bounding box against `[RMin,RMax] × [ZMin,ZMax]`
  with a margin of a few `h`, and fail with both boxes printed.
* **`RMin > 0` when a shape is given.** The `1/r` is not integrable through the
  axis and a D-shape that crosses it is a mis-typed `R0`.
* **`SearchLength`.** `ExtensionConvergence.cpp` measured `VertexConePath` needing
  `6h` and `LevelSetPath` needing `12h` on the Solov'ev level set. Default per
  family rather than to one number, and say in the error message which family is
  in use when a path fails to find its endpoint.

### Acceptance

**Done.** `src/meq/BoundaryShape.{hpp,cpp}` and `[boundary.shape]`, with ten
tests in `tests/unit/BoundaryShapeTests.cpp` and six in `ConfigTests.cpp` — all
MFEM-free, so **CI runs them**, which is a first for anything geometric here.

* The MXH evaluator reproduces `MillerDShape::boundaryPointCerfonFreidberg()`
  **exactly** — worst difference `0` over 720 samples — and differs from Example
  6's `arcsin( δ sin t )` form by **6.25e-4**, matching the ≤6e-4 that fixture's
  header quotes. Both bounds are asserted, the second from *below* as well, so a
  future edit that silently makes the two conventions one fails.
* A configured Miller boundary and the hand-built one in `MillerConvergence.cpp`
  giving the same `levelSet` to round-off on a sample cloud.
* A deliberately non-star-shaped MXH set rejected at load, with the `theta` named.
* An end-to-end Solov'ev run on a configured curved `Γ` reproducing
  `ExtensionConvergence.cpp`'s rates.

---

## 3. Output

### What has to be written, and in what

Three artefacts, two audiences:

| file | for | when |
|---|---|---|
| `<prefix>.mesh` + `<prefix>_psi.gf` + `<prefix>_grad_psi.gf` | MFEM, GLVis, exact restart | always |
| `<prefix>.nc` | analysis, other codes, interpolating restart | always |
| `<prefix>_psistar.gf` | the estimator, and accuracy studies | linear path only — see below |

`Config::OutputConfig` already names the first three (`Config.cpp:341`) and
nothing acts on it. The NetCDF file is new and needs a `getNetCDFFile()` beside
them.

**`ψ*` used not to survive the nonlinear path**, returning ~1e15 silently, and
`postProcess()` throws rather than pass that on. **That is fixed** — §1 of
`../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md` — so `ψ*` may be written on both
paths and the driver need not hedge.

`postProcess()`'s refusal should nevertheless stay until meq has *measured* the
fix on its own benchmarks, because the failure it guards against is silent and
the guard costs nothing. Retire it in the same change that adds the measurement,
not before.

### `B_poloidal` is a relabelling of `q`, not a differentiation

From `refs/Miller.pdf` eq (1), `B = ∇φ × ∇ψ + f(ψ)∇φ` with `∇φ = ê_φ / R`:

```
B_R = -( 1/R ) dpsi/dZ
B_Z = +( 1/R ) dpsi/dR
```

and meq solves for `q = ( 1/r ) grad_bar( psi ) = ( q_r, q_z )` directly. So

```
B_R = -q_z          B_Z = +q_r
```

pointwise, with no derivative taken anywhere. **This is the entire argument for
the mixed method made concrete** — `README.md` says the physically interesting
output is the field and that a mixed method resolves it at the same order as `ψ`
rather than one lower, and this is where that is cashed in. `q` converges at
`k+1`, so `B_pol` does too.

**Do not trust the two signs above.** They are a derivation, and this project's
record on sign conventions in derivations is poor — the Solov'ev source, the `τ`
in eq (8e), and `DarcyForm`'s `-q` were each settled by measurement after an
argument said otherwise. `flux()` already undoes `DarcyForm`'s negation once, so
there are two conventions in play and a sign error here is invisible to every
convergence rate. **Acceptance is a measurement**: `B_R`, `B_Z` from `q` against
central differences of the *exact* Solov'ev `ψ` on a sample cloud, at the
finite-difference floor. That test is cheap and it is the only thing standing
between a released file and a field that points the wrong way.

### The grid

A uniform tensor-product `(R, Z)` grid over the bounding box of `Γ`, with points
outside `Γ` masked. This is the EFIT-shaped thing every downstream tool expects,
and — decisively for §4's interchange route — it is a format a warm start can
interpolate from in `O(1)` per point with no mesh search at all.

```
dimensions:
    R = NR,  Z = NZ,  boundary = NB

variables:
    double R( R )                "Major radius",                 m
    double Z( Z )                "Height",                       m
    double psi( Z, R )           "Poloidal flux per radian",     Wb/rad
    double B_R( Z, R )           "Radial magnetic field",        T
    double B_Z( Z, R )           "Vertical magnetic field",      T
    byte   inside( Z, R )        1 inside Gamma, 0 outside
    double boundary_R( boundary )  the prescribed Gamma, sampled
    double boundary_Z( boundary )

global attributes:
    title, meq_version, git_commit, config_file, source_type,
    polynomial_degree, refinement_levels, newton_iterations,
    final_residual, solved_at
```

`(Z, R)` with `R` fastest is C row-major and matches MaNTA's `(t, x)`. Points
outside `Γ` get `_FillValue = NaN` **and** a zero in `inside`: some tools honour
the fill attribute and some do not, and the mask is what a reader can rely on.

`[output] GridNR`, `GridNZ` configure the resolution, defaulting to something
like 129×129. The boundary polyline is written so that a plot can draw `Γ`
without re-deriving it from the config.

**Not in version one:** a flux-surface `(ψ_N, θ)` grid. It is more useful
physically, and it needs contour tracing — the `FluxSurfaces` driver that went to
`v0-legacy`. Worth doing; not worth blocking the driver on.

### Sampling without hitting the `FindPoints` trap

`mfem::Mesh::FindPoints` is `O(elements × points)` — a brute-force scan over
element centres, recorded in `CLAUDE.md`. A 129×129 grid against a 20k-element
mesh is 3.3e8 element tests for a job that should be instant.

**Invert the loop.** For each element, take its bounding box, convert that
directly to a grid index range — `O(1)`, because the grid is uniform — and test
only the handful of grid points inside it with
`ElementTransformation::TransformBack`. Total cost is `O(elements × points per
element)`, which is linear in both. This is the standard scatter-rather-than-
gather fix and it should be written once, in the output layer.

**This is the one sampling job that does not want GSLIB** (§6). The grid is
uniform, so locating a point is index arithmetic and the search — the part
`FindPointsGSLIB` exists to make fast — is already free. What remains is the
inverse element map, and MFEM supplies that. Roughly thirty lines, no new
dependency, and it stays the right answer even after GSLIB is enabled.

A grid point landing exactly on an inter-element face is ambiguous, `ψ_h` being
discontinuous. Take whichever element claims it first; the jump is `O(h^{k+1})`
and converges away. **GSLIB has a principled answer where this has a shrug** —
`FindPointsGSLIB::SetL2AvgType`, arithmetic or harmonic averaging over the
elements meeting at the point, exactly for L2 fields that are multi-valued on
element boundaries. Both `ψ_h` and `q_h` are L2. If GSLIB is enabled, prefer it
here and delete the shrug.

### Acceptance

* `B_R`, `B_Z` against central differences of the exact Solov'ev `ψ`, at the
  finite-difference floor, on a cloud inside `Γ`.
* The sampled `ψ` against the exact solution at grid points, at the discretisation
  error, on the same benchmark used for the rate tables.
* Sampling cost linear in the grid size and in the element count — timed, since
  the whole point of the inverted loop is a complexity claim.
* A written file read back by `ncdump` and by the warm-start reader of §4.

---

## 4. Warm starts

### Three routes, and they are not the same problem

**Exact restart** — same mesh, same degree. Read the stored trace or potential
grid function and hand it to `setInitialGuess()`. Bitwise resumable, and the
right thing for continuing an adaptive run or re-solving after a profile tweak.

**Full-order interpolating restart** — a different mesh or degree, same code.
Read the stored `ψ_h` and its mesh, and evaluate it at the new space's points
with `FindPointsGSLIB`. The guess then carries the full `k+1` accuracy of the
solve it came from. This needs §6.

**Interchange restart** — another code entirely, or a mesh that no longer exists.
Read `psi(R,Z)` from a NetCDF file and evaluate it by bilinear interpolation on
the structured grid, as an `mfem::Coefficient`.

**An earlier draft of this plan had only the first and third, and that was a
compromise dressed up as a design.** The NetCDF route was chosen *because*
`Mesh::FindPoints` is `O(elements × points)` and there was no other way to
interpolate from a foreign mesh cheaply. But bilinear interpolation on a 129×129
grid is **second order**, so restarting a `k = 3` solve through it discards
almost everything the previous solve computed. It still converges — a guess is
only a guess — but "warm" is doing less work than it looks. With GSLIB the
same-code case never leaves the finite element representation and the compromise
disappears.

So the third route keeps its job and stops pretending to be the general one: it
is the **interchange** format, and it is genuinely the best thing for that. An
outside code has to produce `ψ` on a rectangle and nothing else — no MFEM, no
mesh format, no agreement about element types. That is why §3's file is designed
the way it is, and it is a virtue for foreign input specifically.

```toml
[initialguess]
Type = "gridfunction"      # none | gridfunction | netcdf | ramp
File = "previous_psi.gf"   # gridfunction: full order, needs GSLIB unless the
MeshFile = "previous.mesh" #   mesh matches exactly, in which case it is direct
# Type = "netcdf":
# File = "previous.nc"     # or a file from any other code
# Variable = "psi"         # so a foreign file can name its own field
```

`Type = "ramp"` keeps the interim device the GS-2 benchmarks use today — a
non-homogeneous ramp putting `ψ = 0` in the interior — as a named, documented
option rather than a workaround buried in a test.

### The awkward part: points the old grid does not cover

A new domain can extend beyond the stored grid, and the stored grid is masked
outside `Γ`. Both give an interpolation query with no data. The guess is only a
guess, so the fallback need not be clever, but it must be *stated*: fall back to
the Dirichlet datum, count the misses, and report the count. A restart that
silently found no data for 90% of the domain is a restart that quietly became a
cold start, and the iteration count will not obviously say so.

### Acceptance

* ~~Solve, write, restart on the same mesh: Newton finishes in one step.~~
  **Done, and it did not work when first measured** — see `CLAUDE.md`'s trap on a
  convergence target scaled to `‖r₀‖`. Now **one** Newton iteration.
* ~~Solve, write, restart on a mesh refined once~~ **Done**: `‖r₀‖` falls from
  1.569e+01 to 1.803e-03, a factor of **8704**, Newton finishes in **1** step
  against 4 cold, and the converged `L2` is 3.663930e-07 against 3.663932e-07 —
  seven figures, not six.
* ~~**The two interpolating routes compared on the same restart**, at `k = 3`~~ —
  **done, and the claim needed correcting.** The ratio was stated as `h^{k+1}` to
  `h^2`, which quietly compares two different `h`'s: the mesh's and the *grid's*.
  Measured, sampling a `k = 3, n = 8` solve (own error 6.00e-6) onto grids of 65,
  129 and 257 nodes and reading it back bilinearly gives 1.18e-4, 3.36e-5,
  1.01e-5 — **second order in the GRID spacing**, rates 1.813 and 1.735, drifting
  below two only because by 257 nodes the grid error is within a factor of two of
  the solve's own and the sequence is leaving the asymptotic regime.

  So at 257² the grid route is only **1.7×** worse than full order, not the
  order-of-magnitude the plan implied. The real statement is about *scaling*, not
  a ratio: the grid route's accuracy is a property of the **file**, the
  full-order route's is a property of the **solve**. Refine the mesh or raise the
  degree and the grid has to be refined quadratically to keep pace; the
  full-order route never has to be. That is what justifies the dependency, and
  `fullOrderCarriesMoreThanAStructuredGrid` asserts the scaling rather than the
  ratio.
* A restart whose grid covers nothing reports its miss count and still converges
  from the fallback.
* ~~The GS-2 sources warm-started from a coarse solve, converging where the
  homogeneous cold start lands on `ψ ≡ 0`.~~ **Struck for the same reason §1's
  own version of it was struck, and it should never have been written twice.**
  A warm start is a guess, and a guess alone does not leave the trivial branch:
  §1 measured a bump of 0.05 or 0.20 on §4.3 converging in 4–5 iterations
  straight back to `ψ ~ 1e-16`, because **Newton goes to the root nearest its
  iterate and `ψ ≡ 0` is a root**. A coarse solve of the same homogeneous
  problem is that root, so warming from it cannot change which root is found.
  This bullet contradicted §1's acceptance block on the same page. **The ramp
  stays the workaround** — see §1, and `Type = "ramp"` above.

---

## 5. The driver

Single binary, `meq config.toml`, as agreed. No subcommands. `--help` and a
`--version` that prints the git commit; nothing else, and specifically not
MFEM's `OptionsParser`, which the legacy file used and which wants to own
argument parsing for the whole program.

```
parse config                        -> Configuration
build background mesh               -> box from [mesh], or [mesh].File; refine
build boundary shape                -> levelSet, or none for the fitted path
  if a shape is configured:
    MarkLevelSetSubdomain -> SubMesh D_h -> identify Gamma_h -> TransferPath
build source                        -> makeSource( config.getSource(), ... )
construct GradShafranovSolver
  setSource, setBoundaryData, [ setExtension ], [ setInitialGuess ]
solve                               -> Newton, or a direct solve if linear
  [ adaptive loop: estimate, mark, refine, re-solve ]
postProcess                         -> only on the linear path; log why if not
write                               -> mesh, gf, NetCDF
report                              -> iterations, residual history, timings
```

**DONE.** The adaptive loop is stage-6 work that already existed and has been
exposed rather than rebuilt. What follows is the design as planned; the
differences that measurement forced are recorded after it.

```toml
[adaptivity]
Enabled = false
MaxIterations = 10
Strategy = "doerfler"      # doerfler | maximum
Theta = 0.6
TargetError = 1.0e-6
```

**It inherits an omission that is now repairable.** On the extension path `η₅`
compares `ψ*` against a trace pinned to zero rather than the `φ_h` actually
imposed, so `setTransferredBoundary()` excludes those faces. That was §4 of the
MFEM defects document, and **`φ_h` is now reachable**:
`TransferredDatumCoefficient` is the datum, `PathLiftCoefficient` the
solution-dependent half.

~~So the honest sequence is: the driver calls `setTransferredBoundary()`
automatically whenever an extension is configured and says so in the log, and
that stays true until `η₅` is rebuilt on `TransferredDatumCoefficient` — which
is stage-6 work, not stage-7, and wants its own convergence measurement. Until
then a user should not have to know about the exclusion to get a correct
refinement pattern, but should be able to find out it is in play.~~ **THE
REBUILD HAPPENED, 2026-09-02**, so the exclusion is history and there is nothing
left for the log to warn about. `setTransferredBoundary()` takes the datum as a
second argument — `GradShafranovSolver::transferredDatum()`, which is
`mfem::TransferredDatumCoefficient` — and those faces stay **in**, compared
against the `φ_h` actually imposed. `η₅` goes from 4.07e-1 to **9.6e-5** on the
coarsest extension mesh and converges at **2.78** against 0.40, so the term is
two percent of `η₁` rather than four orders larger than it.
`theTransferredDatumRestoresEtaFive` is the convergence measurement this
paragraph asked for, and it keeps the pinned column as its control.

**The driver rebuilds the datum every cycle, and that is not fussiness.** It
lifts the *solved* flux along the transfer paths, so it goes stale the moment
the solver is solved again — which in an adaptive loop is every cycle, on a new
mesh with a new path family. Hoisting it out would compare `ψ*` against the
previous cycle's boundary condition.

**And `η` itself barely moves**, which is the point rather than a
disappointment: 97 → 254 → 342 → 449 and `η` 4.7352e-04 → 6.8668e-05, unchanged.
What changed is that boundary elements now *have* an `η₅` contribution instead of
none, so the marking can see them.

**That is what was built.** Three things this plan did not anticipate:

* **The curved path needed `meq::AdaptiveDomain`, not the driver's existing
  one-shot `buildSubdomain()`.** Refining `D_h` alone leaves `Γ_h` where it is
  while `h_loc` halves, so `dist/h_loc` doubles every cycle and the transfer
  quietly leaves the regime it is analysed in. Both constructions are therefore
  in the driver, and the comment at the branch says why rather than tidying one
  away.
* **`ψ*` was a hard prerequisite, and retiring `postProcess()`'s refusal was
  listed in `ROADMAP.md` as an independent item.** It is not: four of eq. (20)'s
  five terms are built on `ψ*`, and the driver's path is *always* semi-linear.
* **And the refusal could not be retired outright.** Measuring it found the MFEM
  fix works where `∂F/∂ψ ≠ 0` and returns a different function on **any element**
  where it vanishes — a singular local matrix, factored anyway, because the
  mean-value regularisation is skipped by a flag set on the wrong condition.
  ~~So the loop builds `η` on `Potential::Raw`: one order down, correct, and a
  standing decision rather than a runtime check.~~ **IT WAS NOT A STANDING
  DECISION AND THIS BULLET IS THE STALE HALF OF ITS OWN FINDING.** The hole is
  closed upstream, unconditionally, as *"The postprocessing closes on the
  element average, always"* — the local problem is a pure Neumann one by
  construction, so there was never a condition to test. `apps/meq.cpp` calls
  `postProcess()` and takes the estimator's default `Potential::PostProcessed`;
  the case that read 20.3, 64.1 and 61.6 on dead elements now reads **1.0069**
  against 1.0048 where `∂F/∂ψ` does not vanish at all, and the loop reaches
  `η = 6.87e-5` on **449** elements where the degraded estimator needed 1069 to
  reach 4.77e-4. `Potential::Raw` survives as the *measurement* of the one order
  it costs, not as anything the driver runs. See `CLAUDE.md`, *Post-processing is
  back*; the request itself,
  `../mfem-hdg-dev/doc/HDG-RECONSTRUCT-DEGENERATE-POTENTIAL-MASS.md`, is landed
  and retired and is now on no branch, so that pointer resolves to nothing.

### What the driver owes the user

* **Exit codes that mean something.** 0 solved; 1 configuration error; 2 Newton
  did not converge; 3 output could not be written. A shell script driving
  parameter scans is a first-class caller.
* **A residual history on stdout**, in the shape `CLAUDE.md` records — iteration,
  `‖r‖`, `‖r‖/‖r_0‖`, observed order. It is the diagnostic that separates a wrong
  Jacobian from a hard problem, and it costs nothing to print.
* ~~**`MKL_THREADING_LAYER=GNU`.**~~ **DROPPED, 2026-09-01.** This asked the
  driver to detect MKL and set the variable or refuse to start, on the grounds
  that the failure is silent wrong numbers. It was implemented as a warning and
  has been deleted. meq builds against its own SuiteSparse now and never loads
  the `libmkl_rt` dispatcher the variable configured, so it is inert — measured
  bit-identical across three value-asserting suites — and **most builds link no
  MKL at all**, so the warning told the majority of users to set a variable
  naming a library they do not have. Choosing a threading layer is CMake's job,
  since it can see what the build links. See `CLAUDE.md`, *PARDISO and the MKL
  link line*.

### Acceptance

* `meq examples/solovev-nstx.toml` writing mesh, grid functions and NetCDF, with
  the sampled `ψ` matching the exact solution to the discretisation error — as a
  ctest, so it stays true.
* An example TOML per source type, each one runnable.
* A configuration error in each of the new sections producing a message naming
  the key, checked in `ConfigTests`.
* `README.md` rewritten. `CLAUDE.md` currently warns that its claims about
  testing are aspirational; once there is a driver, most of that file's
  description of meq becomes true for the first time and the warning should go.

---

## 6. GSLIB, and whether to take the dependency

**Recommendation: yes, enabled in the same rebuild as SUNDIALS.** Not as its own
errand, and not blocking anything in §§1–5.

### What it is, and what it costs

`MFEM_USE_GSLIB=YES` brings in `FindPointsGSLIB`: robust evaluation of a grid
function at an arbitrary cloud of points. Per point it returns the element, the
MPI rank, the reference coordinates, a code distinguishing *inside* from *on an
element border* from *not found*, and a distance to the border so that points
just outside the domain can be rejected honestly rather than clamped. It also
brings `SetL2AvgType` (§3), surface and edge variants, a device path, and the
gather-scatter machinery meq has no use for today.

The algorithm is a hash grid over oriented element bounding boxes followed by a
Newton solve in reference space — `bb_t`, `newt_tol` and `npt_max` are exposed
directly in the `FindPoints` signature. Expected `O(1)` per point after an
`O(elements)` setup, against `Mesh::FindPoints`'s `O(elements)` *per point*. It
is the standard tool for this job in both Nek5000 and MFEM's own `field-interp`.

**It is packaged by no distribution.** Checked: no `gslib` in Debian or Ubuntu,
and none in Fedora 43–45, EPEL 8–10.4 or Rawhide. It is an HPC-ecosystem library
obtained from source or through Spack. MFEM vendors device kernels under
`fem/gslib/`, but they are all inside `#ifdef MFEM_USE_GSLIB` and include
`gslib.h`, so there is no partial path — the external library or nothing.

The build is nevertheless small. `../mfem-hdg-dev/INSTALL:831`:

```sh
# as a sibling of the MFEM tree
tar xf gslib-1.0.9.tar.gz && ln -s gslib-1.0.9 gslib
cd gslib && make CC=gcc MPI=0
# then rebuild MFEM with MFEM_USE_GSLIB=YES
```

Pure C, no dependencies of its own, the same sibling-build pattern MFEM already
uses for every optional package.

### Why the cost is lower than "not packaged" makes it sound

**The dependency lands on the MFEM tree, not on meq.** meq's `CMakeLists.txt` is
untouched, `find_package(MFEM)` is unchanged, and the MFEM-free half of the
library — `Config`, `Profiles`, `Source`, `SourceFactory` — is unaffected. CI
already cannot build MFEM at all, since the branch is unpublished, so this adds
nothing to what CI cannot do.

That is a different calculation from adding a dependency to meq itself, and it is
the same reasoning `CLAUDE.md` records for why toml11 is a submodule and MFEM is
not: things meq controls are vendored and pinned; things meq consumes out of tree
are built by hand.

### Where meq would use it

| use | needs GSLIB |
|---|---|
| NetCDF sampling on the uniform `(R,Z)` grid (§3) | **No** — the grid is uniform, so the search is index arithmetic |
| Full-order warm start from a grid function (§4) | **Yes.** This is the case that justifies it |
| Multi-valued `ψ_h`, `q_h` at element boundaries (§3) | No, but `SetL2AvgType` replaces a shrug with an answer |
| Off-grid error measures in the convergence suite | Yes — lifts a cap `CLAUDE.md` already records against `FindPoints` |
| A flux-surface `(ψ_N, θ)` grid; MaNTA coupling | Yes, when either happens |

### The caveat to measure first

**meq meshes are triangles and gslib's findpts is tensor-product.** MFEM handles
this by splitting simplices into quads internally — `mesh_split`, `ir_split`,
`split_element_map` in `fem/gslib.hpp`. It works, and it is a code path that a
quad-based user never exercises. Measure it on a Solov'ev solve before relying on
it, and compare against the inverted loop of §3 on the same points: two
independent samplers agreeing is worth more than either one being plausible.

### Sequencing

`../mfem-hdg-dev` needs `make clean && make -j4` anyway when it returns to
`gf-hdg-subdomains-dev`, and a SUNDIALS rebuild is already on the table for the
globalisation failures. Doing all three at once makes it one rebuild instead of
three. **Never a bare `make -j` in that tree** — see `CLAUDE.md`.

Nothing in §§1–5 waits on this. The driver lands either way; the warm start gets
better when GSLIB arrives, and §4's comparison measurement is how anyone will
know it did.

---

## Risks, and the two that are not mine to close

**The MFEM branch.** Nothing in §2, §3's `ψ*`, or §5's adaptivity can be built
while `../mfem-hdg-dev` is on `gf-hdg-dev`. This is the immediate blocker and it
resolves by a `git checkout` and a `make clean && make -j4` in that tree.

**`Reconstruct()` on the nonlinear path — fixed, and unverified here.** §1 of
the defects document. The adaptive loop is no longer linear-only in principle.
In practice meq has measured none of it against the new library, and a fix
verified only by the tree that made it is exactly the situation this project's
testing stance exists for. Re-run `EstimatorConvergence` on the Newton path
before believing it.

**Sign conventions, again.** `B_R = -q_z`, `B_Z = +q_r` is a derivation from
Miller eq (1), and three of this project's four convention questions were settled
against the derivation. The finite-difference test in §3 is not optional.

**Star-shapedness.** The level-set-by-bisection route is what makes a configured
boundary cheap, and it fails for exactly the strongly-shaped surfaces MXH exists
to describe. The validation in §2 turns that from a wrong answer into an error
message, which is the most that can be done without a different level-set
construction. If a genuinely non-star-shaped `Γ` is ever wanted, the honest
answer is a signed-distance field built from a sampled polygon, and that is a
larger piece of work than anything in this plan.

**`FindPoints`.** The inverted loop in §3 is a complexity claim and complexity
claims rot. Time it in the test rather than asserting it in a comment.

**GSLIB's simplex path.** §6 takes a dependency partly on the strength of a code
path — findpts on split triangles — that a tensor-product user never runs. If it
turns out to be slow or fragile on meq's meshes, the fallback is the third route
of §4 and nothing is lost but the order of the warm start. Measure before
relying, and keep the inverted loop of §3 as the independent check.
