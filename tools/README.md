# Looking at what meq wrote

A run writes the same equilibrium three times, in three formats, for three
different readers. **They are not interchangeable, and picking the wrong one is
the usual way to waste an afternoon.**

| file | what it is | read it with |
|---|---|---|
| `<stem>.mesh`, `<stem>_psi.gf`, `<stem>_grad_psi.gf` | the discrete solution **exactly** — same spaces, same degree, every coefficient | GLVis; meq itself, for an exact restart |
| `<stem>/<stem>.pvd` | VTK, at the solve's own polynomial degree | ParaView, VisIt |
| `<stem>.nc` | ψ and **B** on a uniform `(R, Z)` grid, plus the boundary meq was given | `plot_equilibrium.py`; any downstream tool |

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
`ParaViewDataCollection`'s layout, not a choice meq made; the `Cycle000000/`
directory underneath it is an implementation detail and is not meant to be
opened piece by piece.

**On the curved path the VTK mesh is bent out onto Γ.** The solve happens on
Ω_h, whose boundary is a polygon inscribed in Γ, so the drawn domain would
otherwise have a faceted edge that is not the boundary anybody asked for. meq
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
actually wants, and the one `DRIVER-PLAN.md` §4 specifies for warm-starting meq
from a foreign code: a structured grid interpolates back in `O(1)` per point
with no mesh search.

**Nodes outside the domain carry both a NaN and a zero in `inside`.** Both,
deliberately: some readers honour the fill attribute and some do not, and
`inside` is the one that can always be relied on.

**`boundary_R` / `boundary_Z` are the boundary meq was actually given**, so a
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
`∇̄ψ = r q`, so for a node `p` outside the mesh meq takes its foot `x₀` on Γ_h —
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

## `plot_equilibrium.py`

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

## What is not here

**Nothing plots the `.gf` files**, and nothing should: GLVis already does, and
a Python reader for MFEM's format would be a second implementation of a file
format meq does not own.

**Nothing plots a scan.** A directory of runs at different resolutions or
profile amplitudes is a common thing to want a figure of, and each `.nc` carries
the provenance to build one from (`polynomial_degree`, `elements`,
`adaptive_eta`, ...). Nothing reads across files yet.
