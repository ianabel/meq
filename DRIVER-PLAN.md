# Stage 7: the driver, warm starts, curved boundaries, and NetCDF output

A plan, not an implementation. Nothing here has been built.

meq closes stage 6 with every claim about rates measured and the solver
unusable by anyone who is not running `ctest`. `apps/meq.cpp` is still the
*legacy* driver — it solves the vacuum coil field and calls `WriteOutputMFEM()`
on a class the port deleted. Nothing in `src/` writes a file. This plan is the
work that turns a measured library into a program.

**Written against a tree that cannot currently build it.**
`../mfem-hdg-dev` is on `gf-hdg-dev`, which has no
`fem/darcy/extension_hdg.{hpp,cpp}`, so `GradShafranov.cpp` does not compile and
§4 of this plan has nothing to stand on. Everything below assumes the tree is
back on `gf-hdg-subdomains-dev`. That is a scheduling fact, not a design
constraint.

## The five pieces, and why this order

| | | ends at |
|---|---|---|
| 7a | `setInitialGuess()` | a source that vanishes at `ψ = 0` converging to something that is not zero |
| 7b | Boundary shapes: Miller and MXH | a configured `Γ` reproducing `MillerDShape.hpp` to round-off |
| 7c | Output: mesh, grid functions, NetCDF | `B` from `q` matching a finite difference of the exact `ψ` |
| 7d | Warm start | a restart from 7c's file cutting Newton to one or two steps |
| 7e | The driver | `meq examples/*.toml` writing files, end to end, as a ctest |

7a before 7e because a driver without it is *actively misleading*: handed a
physical profile it converges in zero iterations and writes a file full of
zeros, which looks exactly like success. 7c before 7d because the warm-start
input format **is** the output format — see §5. 7b is independent and could be
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

* All four GS-2 §4.2–4.5 sources posed with *homogeneous* data — as their own
  paper poses them — and a non-trivial guess, converging to a solution with
  `‖ψ‖ > 0`. Today they need a non-homogeneous ramp to avoid the trivial branch;
  that workaround should become unnecessary and the tests should say so.
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

* The MXH evaluator reproducing `MillerDShape::boundaryPointCerfonFreidberg()` to
  round-off at `N = 1`, `s_1 = arcsin( delta )` — a direct check of the reduction
  MXH eq (4) asserts.
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

**`ψ*` cannot be written on the nonlinear path**, and the reason is MFEM's, not
meq's: `ReconstructFluxAndPot()` reads only the linear `M_p`, so through Newton
it returns ~1e15 and `postProcess()` throws rather than pass it on. This is §1 of
`../mfem-hdg-dev/doc/HDG-DEFECTS-FROM-MEQ.md`. The driver must therefore *not*
promise `ψ*` unconditionally — write it when the solve was linear, and say
plainly in the log why it is absent otherwise. Writing a 1e15 field would be the
worst available outcome.

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
and — decisively for §5 — it is the format a warm start can interpolate from in
`O(1)` per point.

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
mesh is 3.3e8 element tests for a job that should be instant, and MFEM here is
built with `MFEM_USE_GSLIB = NO`, so `FindPointsGSLIB` is not available as an
escape.

**Invert the loop.** For each element, take its bounding box, convert that
directly to a grid index range — `O(1)`, because the grid is uniform — and test
only the handful of grid points inside it with
`ElementTransformation::TransformBack`. Total cost is `O(elements × points per
element)`, which is linear in both. This is the standard scatter-rather-than-
gather fix and it should be written once, in the output layer, and reused by the
warm-start reader.

A grid point landing exactly on an inter-element face is ambiguous, `ψ_h` being
discontinuous. Take whichever element claims it first; the jump is `O(h^{k+1})`
and converges away. Say so in a comment rather than pretending the choice is
principled.

### Acceptance

* `B_R`, `B_Z` against central differences of the exact Solov'ev `ψ`, at the
  finite-difference floor, on a cloud inside `Γ`.
* The sampled `ψ` against the exact solution at grid points, at the discretisation
  error, on the same benchmark used for the rate tables.
* Sampling cost linear in the grid size and in the element count — timed, since
  the whole point of the inverted loop is a complexity claim.
* A written file read back by `ncdump` and by the warm-start reader of §5.

---

## 4. Warm starts

### Two routes, and they are not the same problem

**Exact restart** — same mesh, same degree. Read the stored trace or potential
grid function and hand it to `setInitialGuess()`. Bitwise resumable, and the
right thing for continuing an adaptive run or re-solving after a profile tweak.

**Interpolating restart** — different mesh, different degree, different geometry,
or another code entirely. Read `psi(R,Z)` from a NetCDF file and evaluate it by
bilinear interpolation on the structured grid, as an `mfem::Coefficient`.

The second is why §3's file is designed the way it is. A structured grid makes
interpolation `O(1)` per point with no mesh search at all, which is what makes
"warm start from elsewhere" cheap enough to be routine. It also means the *only*
thing an outside code has to produce is `ψ` on a rectangular grid — no MFEM, no
mesh format, no agreement about element types.

```toml
[initialguess]
Type = "netcdf"            # none | netcdf | gridfunction | ramp
File = "previous.nc"       # netcdf: reads psi(R,Z) from a meq output file
Variable = "psi"           # so a file from elsewhere can name its own field
MeshFile = "previous.mesh" # gridfunction: exact restart, with .gf beside it
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

* Solve, write, restart on the same mesh: Newton finishes in one step.
* Solve, write, restart on a mesh refined once: Newton finishes in fewer steps
  than from cold, with the converged answer unchanged to six figures.
* A restart whose grid covers nothing reports its miss count and still converges
  from the fallback.
* The GS-2 sources warm-started from a coarse solve, converging where the
  homogeneous cold start lands on `ψ ≡ 0`.

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

The adaptive loop is stage-6 work that already exists and should be exposed
rather than rebuilt:

```toml
[adaptivity]
Enabled = false
MaxIterations = 10
Strategy = "doerfler"      # doerfler | maximum
Theta = 0.6
TargetError = 1.0e-6
```

**It inherits a known omission.** On the extension path `η₅` compares `ψ*`
against a trace pinned to zero rather than the `φ_h` actually imposed, so
`setTransferredBoundary()` excludes those faces — §4 of the MFEM defects
document. The driver should call it automatically whenever an extension is
configured, and the log should say that it did. A user should not have to know
about an MFEM limitation to get a correct refinement pattern, but they should be
able to find out that one is in play.

### What the driver owes the user

* **Exit codes that mean something.** 0 solved; 1 configuration error; 2 Newton
  did not converge; 3 output could not be written. A shell script driving
  parameter scans is a first-class caller.
* **A residual history on stdout**, in the shape `CLAUDE.md` records — iteration,
  `‖r‖`, `‖r‖/‖r_0‖`, observed order. It is the diagnostic that separates a wrong
  Jacobian from a hard problem, and it costs nothing to print.
* **`MKL_THREADING_LAYER=GNU`.** CMake sets it on registered ctests; a user
  running the binary by hand gets no such help, and the failure is silent wrong
  numbers. The driver should detect MKL and either set it or refuse to start.
  This is the one trap in `CLAUDE.md` that a released binary can actually
  protect its user from.

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

## Risks, and the two that are not mine to close

**The MFEM branch.** Nothing in §2, §3's `ψ*`, or §5's adaptivity can be built
while `../mfem-hdg-dev` is on `gf-hdg-dev`. This is the immediate blocker and it
resolves by a `git checkout` and a `make clean && make -j4` in that tree.

**`Reconstruct()` on the nonlinear path.** §1 of the defects document. Until it
is fixed, no nonlinear run can write `ψ*` and no nonlinear run can drive the
estimator, which means **the adaptive loop is linear-only**. That is a real
restriction on §5 and it should be stated in the driver's log rather than
discovered.

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
