# Meshing

MEQ's own meshes come from `mfem::Mesh::MakeCartesian2D` — a rectangle of
triangles — with curved and non-rectangular geometry reached by the extension
technique on top of it. That is enough for every fixed-boundary configuration in
`examples/`, and it is not enough for free boundary.

`halfdisc.py` is the generator for what free boundary needs.

```sh
python3 tools/mesh/halfdisc.py --rho 1.5 --size 0.15 \
    --coil 0.90 -0.86 0.25 0.12 --coil 0.90 0.74 0.25 0.12 \
    --plasma 0.30 -0.55 0.70 1.10 --plasma-size 0.035 \
    --check -o machine.msh
```

`File = "machine.msh"` under `[mesh]` is how a run reads it — MFEM
reads gmsh's format natively, with no converter.

## Or drive it from the configuration itself

`[mesh.generate]` puts the geometry above in the same TOML file as the
equilibrium, and `meq-run` makes the mesh and then solves:

```sh
./build/meq-run examples/diverted-tokamak-generated.toml
```

**The conductors are then written once.** The `--coil` rectangles are derived
from the file's `[[coils]]` blocks, in file order — which is exactly the
agreement *What the attributes mean* below records, and which is what makes the
derivation possible at all. The command that results is printed by `meq
--mesh-command <config>`, and written beside the mesh as `<mesh>.meq-mesh` so a
re-run that would produce the same command skips the meshing.

Nothing about this file changes: it is still a plain gmsh script with the CLI
below, `meq-run` invokes it through that CLI, and `[mesh.generate]`'s schema is
MEQ's own — four bounds where `--plasma` takes a corner and two extents, a
centre and half-extents where `--coil` takes a corner and two extents. The
driver converts.

## Why r = 0 is the requirement, not a preference

`FREE-BOUNDARY-PLAN.md` §3 expands the exterior field in Gegenbauer functions
of order −1/2 on a **semicircle centred on the axis**. The expansion is a
statement about that geometry and about no other, and `meq::ExteriorDtN`'s
diagonal symbol is the DtN map of that shape.

So a domain stopping at `r = 0.05` is not a slightly worse version of it — the
basis does not span the exterior of that shape at all, and nothing downstream
would report the difference: the solve converges, at full order, to the wrong
exterior problem. `--check` therefore asserts on `r == 0.0` **exactly** rather
than within a tolerance. The axis is a straight geometric edge at `x = 0`, so
every node OCC places on it lands there bit exactly, at geometric order 1 and
at order 2 alike.

Read `--rho` as *where Γ is*, and note that Γ must **enclose every conductor**:
outside it the field has to be Δ\*-harmonic, and a current there is not. The
tool refuses a coil that reaches or crosses Γ, and refuses one that reaches the
axis for the same reason `meq::Coil` and `[[coils]]` do.

**With the exterior coupling on, `--rho` is one step further out than that.**
`[boundary.exterior] Radius` is Γ, and Ω_h is cut *from* this mesh as the
elements lying inside it — so the arc drawn here is the **background's** outer
edge and Γ is a smaller semicircle within it, with the band between the
resulting staircase and Γ bridged by the same Cockburn–Solano transfer the
curved path uses. The driver refuses a Γ that does not fit strictly inside the
mesh, and with `[mesh.generate]` it refuses it before gmsh runs.

## What the attributes mean

| | |
|---|---|
| element attribute `1` | everything that is not a conductor |
| element attribute `10 + i` | the `i`-th `--coil`, **in command-line order** |
| boundary attribute `1` | Γ, the arc — carries the transferred exterior datum |
| boundary attribute `2` | the axis — carries a plain essential condition |

The coil numbering is by command-line position rather than by whatever tag OCC
assigned, so a caller writing `[[coils]]` blocks in the same order has the two
lists lined up. That is the only way the attributes are usable.

**Γ and the axis must not share an attribute**, since they carry different
conditions, and they are told apart by geometry rather than by tag order —
which OCC does not promise and which `fragment` reshuffles.

## The coils are meshed to, not cut

A conductor's geometry is **prescribed input**: it does not move with the
solution, so it can always be aligned to. FB-2 measured what that is worth —
rates of 1.33 / 1.27 / 1.09 with the conductor cut by elements, against
1.99 / 2.88 / 3.01 aligned. `occ.fragment` is what does it.

Only the plasma support genuinely moves, and FB-4 measures what *that* costs.

## Coarse disc, refined plasma

Γ has to sit far enough out for the exterior expansion to converge, and the
plasma occupies a small part of what that encloses. Meshing the whole disc at
the plasma's resolution spends most of the elements on vacuum that carries no
source and needs no accuracy. `--plasma`/`--plasma-size` refine a box, with a
graded transition out of it; `--coil-size` does the same inside the conductors.

Measured, two coils, plasma box `r ∈ [0.30, 1.00]`, `z ∈ [-0.55, 0.55]` at
`h = 0.035`:

| `--rho` | graded | uniform at `h = 0.035` | saving |
|---|---|---|---|
| 1.5 | 2,965 | 7,672 | 2.6× |
| 3.0 | 3,368 | 30,392 | **9.0×** |

The saving grows with `rho`, which is the direction the exterior coupling
pushes: moving Γ from 1.5 to 3.0 costs the graded mesh 14% more elements and
the uniform one four times as many.

## Always pass `--check` when the geometry changes

It re-reads the written file — not the model still in memory, which would share
every assumption with the code that built it — and asserts that `r` reaches 0
exactly, that Γ and the axis are the outer boundary and nothing else, and that
each coil attribute covers exactly its rectangle.

That last group exists because the first version of this tool got it wrong.
Γ and the axis were selected by asking each 1-D entity where its centre of mass
sat, and after `fragment` the model also contains **each conductor's own
outline** — every one of which is "not on the axis", so they joined Γ. A 1-D
physical group is written to the `.msh`, and MFEM's reader turns 1-D elements
into *boundary* elements whether or not they are topologically on the boundary.
Measured: 64 boundary faces, of which 32 on the arc, 20 on the axis and **12 on
the coil rectangles**, carrying Γ's attribute in the middle of the mesh.

That is not a cosmetic mis-tag: attribute 1 carries the transferred exterior
datum, so the run would have imposed Γ's Dirichlet condition on twelve edges
buried inside the conductors. `getBoundary( surfaces, combined = True )` is the
fix, and the check is what says it stayed fixed.

## Why Python, and why gmsh is not linked into `meq_core`

gmsh ships a C++ API and a Python module, and **the Python module is `ctypes`
over the same `libgmsh.so`**, generated from the same `api/gen.py`: identical
function set, identical meshes, the same code doing the work. There is no
capability or speed argument between them.

What linking would buy is in-process access to gmsh's *in-memory* model — and
nothing MEQ does needs it, because MFEM reads the `.msh` natively and the file
is the whole interface. What it would cost is `libgmsh-dev` as a build
dependency of a tree whose CI cannot build the solver at all, and whose
meshing would then become untestable there too.

The one thing that would change the answer is a run that re-meshes *during* the
solve — moving conductors, or a shape optimiser. MEQ's adaptive loop *refines*
an existing mesh through `meq::AdaptiveDomain` rather than re-meshing, so that
is not on the plan; and if it arrives, this script transliterates, because it
is the same API.
