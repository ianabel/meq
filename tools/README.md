# Looking at what MEQ wrote

A run writes the same equilibrium three times, in three formats, for three
different readers, and on request a fourth file that is not the equilibrium at
all but a **reduction** of it. **They are not interchangeable, and picking the
wrong one is the usual way to waste an afternoon.**

| file | what it is | read it with |
|---|---|---|
| `<stem>.mesh`, `<stem>_psi.gf`, `<stem>_grad_psi.gf` | the discrete solution **exactly** — same spaces, same degree, every coefficient | GLVis; MEQ itself, for an exact restart |
| `<stem>/<stem>.pvd` | VTK, at the solve's own polynomial degree | ParaView, VisIt |
| `<stem>.nc` | ψ and **B** on a uniform `(R, Z)` grid, plus the boundary MEQ was given | `plot_equilibrium.py`; any downstream tool |
| `<stem>_surfaces.nc` | the flux surfaces and the flux-surface averages, against a flux label. **Only when asked for** | a 1-D transport code; `plot_equilibrium.py` |

## The three, and when each is the right one

**The MFEM files are the only lossless ones.** `_psi.gf` carries the P_k
coefficients, so a restart from them is bit-identical to the answer that was
written. Nothing outside MFEM reads them. Note `_grad_psi.gf` is the HDG flux
`q`, in `DarcyForm`'s sign convention — **not** the poloidal field. The two
differ by a relabelling (`B_R = −q_z`, `B_Z = +q_r`); see `src/meq/Field.hpp`.

```sh
glvis -m run.mesh -g run_psi.gf
```

**The VTK files are for looking at.** They carry ψ and **B** — the physical
field, already relabelled — written as VTK Lagrange cells at the polynomial
degree of the solve. That last part matters more here than in most codes: VTK's
native cells are linear, so the default path would draw a `k = 3` solution as
though it were `k = 1`, and the result looks like a coarse mesh rather than
like a bug. Open the `.pvd`, which is the index:

```sh
paraview run/run.pvd
```

**An adaptive run also writes `<stem>_cycles/`, one frame per cycle**, which
ParaView reads as a time series and will scrub through — mesh and all, so the
refinement can be watched rather than inferred from the table of element counts.
Measured on `examples/miller-adaptive.toml`: 97 → 254 → 342 → 449 elements over
four frames.

```sh
paraview miller-adaptive_cycles/miller-adaptive_cycles.pvd
```

It is a **separate collection from the answer**, and deliberately: `<stem>` is
the converged state and has its boundary bent onto Γ, which mutates the mesh —
doing that mid-loop would hand the next refinement a geometry the estimator
never saw. So the frames are **uncurved**, Γ_h as actually solved on, which is
the honest thing to animate.

**The `.pvd` is inside the collection directory, not beside it.** That is
`ParaViewDataCollection`'s layout, not a choice MEQ made; the `Cycle000000/`
directory underneath it is an implementation detail and is not meant to be
opened piece by piece.

**On the curved path the VTK mesh is bent out onto Γ.** The solve happens on
Ω_h, whose boundary is a polygon inscribed in Γ, so the drawn domain would
otherwise have a faceted edge that is not the boundary anybody asked for. MEQ
installs a curvature on the mesh and moves each boundary face onto Γ — and
because the VTK is already written as Lagrange cells, **this needs nothing
further from the format**; the two features compose.

Moving a boundary by `O(h)` can turn an element inside out, and the gap can
exceed an element's own size, so the displacement is backed off **per node**
wherever an element would fold. The driver reports the fraction that reached Γ:
96% on `miller-curved`, 88% on `miller-adaptive`. The rest sits between Γ_h and
Γ. Comparing the VTK edge against the `.nc` `boundary_R/Z` will show that
difference, and it is the mesh being coarse there rather than an error.

**The NetCDF file is the interchange format**, and the only lossy one — a `k+1`
field sampled onto a rectangle. Its grid is `[output] GridNR × GridNZ`, which
has nothing to do with `[mesh] NR/NZ`. It is the format every downstream tool
actually wants, and the one for warm-starting MEQ from a foreign code: a
structured grid interpolates back in `O(1)` per point with no mesh search. What
that costs is second order in the *grid* spacing at every `k`, so its accuracy is
a property of the file rather than of the solve — `docs/output.rst` and
`docs/running.rst` have the statement and the reason it decides which route a
restart takes.

**Nodes outside the domain carry both a NaN and a zero in `inside`.** Both,
deliberately: some readers honour the fill attribute and some do not, and
`inside` is the one that can always be relied on.

**`boundary_R` / `boundary_Z` are the boundary MEQ was actually given**, so a
plot can show the answer against what was asked for rather than against the
mesh's own edge. `boundary_source` says which it is: `shape (smooth Gamma)` on
the curved path, sampled from the shape itself at 512 points; `mesh boundary` on
the fitted path, where Γ *is* the mesh boundary and is walked out of it as an
ordered loop. The loop is **not** closed — a reader appends the first point,
which `plot_equilibrium.py` does.

### The band between Γ_h and Γ, and how far to trust it

On the curved path Ω_h is the union of background elements lying **inside** Γ,
so Γ_h is inscribed and there is a band `O(h)` wide that is inside the plasma
and outside the mesh. **The solve makes no claim there** — the discrete problem
is posed on Ω_h — so anything drawn in that band is a continuation, and
`extrapolated_nodes` says how many nodes it covers.

**The continuation is carried by the FLUX, and that is the mixed method paying
off somewhere unexpected.** `q` is computed at the *same* order as ψ, and
`∇̄ψ = r q`, so for a node `p` outside the mesh MEQ takes its foot `x₀` on Γ_h —
a point on the owning element's own boundary — and steps out:

```
psi(p) = psi(x0) + r0 q(x0) . (p - x0)
```

**Nothing is ever evaluated outside an element.** The obvious alternative,
continuing ψ_h's own polynomial past its element, is bounded by nothing:
measured on `examples/miller-curved.toml`, where `[boundary] Type = "zero"`
makes ψ exactly zero on Γ and strictly negative inside, it put **17 nodes at
positive ψ**, worst +1.06e-02 against a peak of 2.5e-01. The flux version puts
**none**, with a maximum of −5.9e-05. Its band error converges at the flux's own
rate — measured 3.75 — where an extrapolation does not converge in the band at
all.

## The fourth file: `<stem>_surfaces.nc`, and why the other three cannot do its job

**Off by default.** `[output] FluxSurfaces = true` turns it on, and
`FluxSurfaceCount`, `FluxAngleCount`, `FluxInnerCut` and `FluxOuterCut` size it.

The other three files are all the same object in different resolutions: a field
on a domain. A 1-D transport code does not want a field on a domain. It wants
**scalar functions of a flux label** — `V′(ρ)`, `⟨R^{-2}⟩(ρ)`, the shape of each
surface — and getting those out of a rasterised ψ means contour tracing,
flux-surface quadrature and a critical-point search at the far end, by a code
that does not have `q`. MEQ does have `q`, so it does the reduction once and
writes the answer.

The layout is `flux × theta`, θ fastest, and the variables are named in
`src/meq/Output.hpp`. Three things about it are worth knowing before reading
one.

**The label is `ρ = √Ψ_N`, and the file says so** in the `flux_label`
attribute. `Ψ_N` is carried beside it so a consumer with its own normalised-flux
grid need not square anything, but ρ is what the geometry is smooth in: ψ has a
quadratic maximum at the axis, so a surface's minor radius grows like `√Ψ_N` and
parametrising by `Ψ_N` puts a square-root branch point on the axis.

**Both ends are cut, for different reasons, and the defaults are
`Ψ_N ∈ [0.05, 0.95]`.** At the inner end a surface shrinks to a point and
`dρ/dψ` is unbounded; at the outer end the surface stops being made of solved
data. Nothing *fails* at either end — which is exactly why the cut has to be a
decision rather than something discovered. See below.

**`extrapolated(flux, theta)` is a mask and `extrapolated_nodes` is a count.**
On the curved path Ω_h is inscribed in Γ, so the outer surfaces cross the band
and the extension answers for those nodes. A surface can be inside Ω_h at one θ
and outside it at the next, so `band(flux)` — the per-surface summary — under-
reports in the middle of a band excursion, which is where a `q(ψ)` profile is
being read. Drop nodes with `extrapolated = 1` before computing an error norm or
differencing two runs.

### Where the outer cut has to go, and why nothing tells you

Measured on a curved Miller boundary at two resolutions, tracing at seven
levels: the share of each surface's nodes that are **band data** rather than
solved data.

| `Ψ_N` | 0.50 | 0.80 | 0.90 | **0.95** | 0.98 | 0.99 | 0.995 |
|---|---|---|---|---|---|---|---|
| `h = 0.0425` | 0% | 1% | 15% | **41%** | 80% | 94% | 98% |
| `h = 0.0212` | 0% | 0% | 2% | **28%** | 46% | 78% | 91% |

**Every one of those traces closed, every fit converged, no ray stalled, and
`|ψ_h − c|` sat at 2e-13 throughout.** The residual does not know the difference
between an element and an extension, because the extension answers as
confidently as an element does. So the cut cannot be found by pushing outward
until something breaks; the mask is the only signal there is, and 0.95 is where
MEQ stops by default — close to LIUQE's 0.95 and inside FreeGS's 0.99, and
chosen here because it is where the outermost surface is still mostly solved
data at production resolutions.

A query outside the cut is **refused**, not extrapolated. FreeGS extrapolates;
MEQ does not, because a plausible `V′` past the boundary is worse than no answer.

## `plot_equilibrium.py`

**It reads both NetCDF files**, and tells them apart by their own dimensions
rather than by the name — the surfaces file is the one with a `flux` dimension —
so a renamed copy still plots and a `--what` meant for the other file is refused
rather than quietly drawing something else.

```sh
tools/plot_equilibrium.py run.nc                            # both, to a window
tools/plot_equilibrium.py run.nc --what surfaces -o psi.png
tools/plot_equilibrium.py run.nc --what field --levels 40 -o b.png
```

* `--what surfaces` — ψ contours, which **are** the flux surfaces, with a
  colourbar. `ψ = 0` is drawn heavier: on the fixed-boundary problem that level
  set is the plasma boundary by definition.
* `--what field` — `|B_pol|` as a colour map.
* `--what both` — the surfaces in white over the field. The default.

Needs `numpy`, `matplotlib` and `netCDF4`. It picks the `Agg` backend
automatically when writing to a file, so it works headless without `--show`.

The figure's subtitle is provenance read out of the file's own attributes —
degree, element count, Newton iterations, whether the boundary was fitted or
curved, the final residual. A directory of scan output is unreadable otherwise.

### The same script on `<stem>_surfaces.nc`

```sh
tools/plot_equilibrium.py run_surfaces.nc -o family.png
tools/plot_equilibrium.py run_surfaces.nc --what geometry -o poloidal.png
tools/plot_equilibrium.py run_surfaces.nc --what profiles -o averages.png
tools/plot_equilibrium.py run_surfaces.nc --profiles V_prime,inverse_R_squared
```

* `--what geometry` — the traced surfaces in the poloidal plane, coloured by
  `ρ`, with the magnetic axis from the file's own attributes.
* `--what profiles` — one panel per flux-surface average, against `ρ`.
* `--what both` — the plane beside the panels. The default.

**The panel list comes out of the file, not out of the script.** Every variable
on the `flux` dimension that is not a coordinate gets a panel, in the order the
file carries them, so `safety_factor` — which is written only when a caller
supplied `g(ψ)`, and is **absent rather than zero** otherwise — appears without
the script being told it exists. `--profiles` takes a comma-separated list to
pick a few, or `all` to add the per-surface diagnostics `worst_residual` and
`transversality`. Naming a variable the file does not carry lists the ones it
does, which is the quickest way to find out whether a run has a safety factor.

**The band is drawn per node, in red, over the surface it belongs to.** A
surface can be inside Ω_h at one θ and outside it at the next, so which arc of a
curve is solved data is not something the curve can be read for — the extension
answers as confidently as an element does. Nodes as well as a line, because the
innermost surface to reach the band does so at a handful of isolated nodes.
Beside the averages there is a **band fraction** panel, which is the per-node
mask summed per surface: `band(flux)` says a surface is affected and this says
by how much, and a surface one node into the extension and one nine tenths of
the way through it are not the same statement about `V′(ρ)`. On the profile
panels the affected surfaces are ringed — there `band(flux)` is the right flag,
since an average is one number over a whole surface and a single continued node
contaminates all of it.

**The cut is drawn as a cut.** The profile axes span the whole of `ρ ∈ [0, 1]`
with the two ends shaded, so a family over `Ψ_N ∈ [0.05, 0.95]` looks like what
it is rather than like an answer that happens to start somewhere. Nothing is
drawn in the shaded region: MEQ refuses a query there and so does the plot.

**The label on the axes is `ρ = √Ψ_N`**, with `Ψ_N` on a secondary axis along
the top row so neither has to be squared by eye. Which one a curve is against
changes its shape near the axis and nothing in the numbers says so, which is why
both are on the picture.

## What is not here

**Nothing plots the `.gf` files**, and nothing should: GLVis already does, and
a Python reader for MFEM's format would be a second implementation of a file
format MEQ does not own.

**Nothing plots a scan.** A directory of runs at different resolutions or
profile amplitudes is a common thing to want a figure of, and each `.nc` carries
the provenance to build one from (`polynomial_degree`, `elements`,
`adaptive_eta`, ...). Nothing reads across files yet.

## Making a mesh

`mesh/halfdisc.py` generates the half-disc-with-conductors that free boundary
needs and MFEM's built-in mesher cannot make: a semicircle **reaching r = 0
exactly**, with the coil rectangles meshed to as their own subdomains, coarse
over the vacuum and refined where the plasma is. `mesh/README.md` says why each
of those is a requirement rather than a preference — the exterior expansion is
a statement about that exact geometry — and how to read the attributes it
writes.

```sh
python3 tools/mesh/halfdisc.py --rho 1.5 --size 0.15 \
    --coil 0.90 -0.86 0.25 0.12 --plasma 0.30 -0.55 0.70 1.10 \
    --plasma-size 0.035 --check -o machine.msh
```

Pass `--check` whenever the geometry changes: it re-reads the written file and
asserts the properties MEQ depends on. It exists because the first version of
the tool tagged twelve interior conductor edges as Γ, which nothing downstream
would have reported.

## Clearing up afterwards

**A MEQ run leaves seven files and a directory in whatever directory it was
started from**, which is usually the repo root, and a few dozen runs turn that
into several hundred paths. `clean.sh` removes them.

```sh
tools/clean.sh --dry-run      # say what would go
tools/clean.sh                # run artefacts and editor scratch
tools/clean.sh --builds       # also the build-*/ experiments, keeping ./build
```

**WHAT MAKES IT SAFE IS `git check-ignore`, NOT THE PATTERN LIST.** Every
candidate is confirmed ignored before it is removed, so a file git tracks — or
one that is merely untracked, which is what work in progress and an incoming
report from another project look like — is skipped whatever it is called. That
is why `examples/limited-tokamak-guess.{mesh,gf}` survive a rule that deletes
`*.mesh` and `*.gf`: `.gitignore` negates them, so `check-ignore` says they are
not ignored, and the script never has to know they are special. A directory goes
only when git tracks nothing inside it *and* nothing inside it is
untracked-but-not-ignored.

`refs/*.pdf`, `.venv-docs/`, the benchmark's `venv/` and its `ref-n*/`
reference solutions are kept even though all four are ignored — they are
fetched or built one at a time rather than regenerated by a command, which is
a different thing from being large. `./build` is kept unless you ask twice,
with `--all-builds`, because every recipe in `CLAUDE.md` names it.
