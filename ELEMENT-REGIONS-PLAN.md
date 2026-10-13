# Element regions in the output

**What this adds.** For every element of the written mesh, the fraction of it
occupied by each material: vacuum, plasma, each coil, and, reserved for later,
passive and ferromagnetic structures. Each coil's circuit is named too. From
these the file derives a region per element for plotting, with the cells the
plasma boundary cuts marked. It all goes in the `.nc` and the VTK. The grid
also gets a **pointwise** `node_material` mask, which is MR-5's deferred `.nc`
support mask.

**Status: DONE.** ER-0 → **[M-266](MEASUREMENTS.md#m-266)**,
**[M-267](MEASUREMENTS.md#m-267)**; ER-1 → **[M-268](MEASUREMENTS.md#m-268)**;
ER-2 → **[M-269](MEASUREMENTS.md#m-269)**; ER-3 →
**[M-270](MEASUREMENTS.md#m-270)**.

## 1. What a reader cannot get from the rest of the output

* **Plasma.** The `.nc`'s `inside` says which grid nodes are in `Ω_h`, but not
  which carry plasma current. A mirror's `n_s` and `e_phi_0` are evaluated
  outside the support, and need a mask to be read by. A tokamak reader
  recovers the plasma by thresholding `psi` at `psi_bnd`. That threshold is
  wrong wherever the connectivity fill or `ExcludeAttributes` removed a lobe,
  and wrong on a mirror, where the support is a box in `( Ψ, z )`.
* **Conductors.** The element attributes `10 + i` survive in the `.mesh` and in
  the VTK's `attribute` cell data, but only for a mesh `halfdisc.py` made, and
  only for meshed coils. Under `Model = "subtracted"` or `"filament"` the copper
  is not separately tagged. Circuits are resolved away at parse time, so
  nothing downstream knows them.
* **The plasma boundary.** Nothing says which elements it cuts. Those are the
  elements whose source integrand switches inside them: they cap the order
  (M-09, M-10) and are where refinement pays.

## 2. The design

### 2.1 Materials and fractions

**Every element stores the fraction of it occupied by each material.** A
material is one row of a table, and the table is built to hold kinds MEQ does
not produce yet:

| `material_kind` | name | rows | today |
|---|---|---|---|
| 0 | `vacuum` | one | produced |
| 1 | `plasma` | one | produced |
| 2 | `coil` | one per `[[coils]]` block | produced |
| 3 | `passive` | one per passive structure | **reserved**: no such input yet |
| 4 | `ferromagnetic` | one per iron region | **reserved**: no such input yet |

Each row has these columns. A column that does not apply to a kind holds its
fill value, so adding a kind adds rows and no new variables:

* `material_name` (string);
* `material_circuit` (int, −1 for none). Passive structures have none;
* `material_current` (A; the resolved current of a coil);
* `material_resistivity` (Ω m; passive);
* `material_relative_permeability` (1 except linear iron);
* `material_permeability_curve` (int, −1 for none), an index into the
  permeability curve table, for non-linear iron.

Circuits get their own table: `circuit_name`, `circuit_current` and
`circuit_supplied`.

**Permeability curves are stored in the file**, so that a `.nc` describes its
own equilibrium without a library of steels beside it. A curve is iron's
`B(H)`, which saturates, so that `μ` depends on the solution and the operator is
non-linear inside the iron. A ferromagnetic row holds either a constant
`material_relative_permeability` or a curve, never both; the other is its fill
value. The curve table is on a dimension `permeability_curve`:

* `permeability_curve_material` (string), **the name of the material it
  belongs to**. It repeats `material_name` on purpose, so a reader can check
  that the index and the name agree. A curve belongs to one material, and two
  regions of the same steel hold two copies;
* `permeability_curve_point_count`, with
  `sample_dimension = "permeability_curve_point"`;
* `permeability_curve_h` (A/m) and `permeability_curve_b` (T), on
  `permeability_curve_point`, increasing in `H` from `( 0, 0 )`. An
  `interpolation` attribute names the rule the solver uses between the points.

**Every coil outside a `[[circuits]]` block belongs to one circuit named
`default`.** It is a grouping, not a supply: its members carry independent
currents, `circuit_supplied = 0`, `circuit_current` is the fill value, and
each coil's own current is its `material_current`. A `[[circuits]]` block is
one current with fixed multipliers, so `circuit_supplied = 1`. The table has a
`default` row only when some coil is unwired, so a file whose coils are all
wired has no empty group.

**`default` becomes a reserved name.** `[[circuits]] Name = "default"` is
refused at parse, since it would otherwise merge a supplied circuit with the
unwired coils. This is the one configuration change the plan makes.
`CircuitParameters`' "a coil in no circuit is its own one-member circuit" is
unchanged for the solve: the `default` group exists only in the output.

**Fractions sum to one** over an element's materials, at the solve's quadrature
points. Each fraction is the quadrature-weighted share of the element's area
(`ip.weight × |J|`) where that material holds. The derived per-element `region`
is a convenience for plotting, not the record:

| `region` | when |
|---|---|
| `vacuum`, `plasma`, `coil`, `passive`, `ferromagnetic` | that material's fraction is 1 |
| `plasma_boundary` | plasma and vacuum share the element |
| `conductor_boundary` | a conductor and vacuum share the element (§2.1b) |
| `mixed` | anything else, which §2.1b says is a fault or a mesh to fix |

### 2.1b Which cells may be mixed

Each pairing has a different cause, and the classifier treats them differently:

* **Plasma and vacuum: unavoidable.** The plasma boundary is a level set of the
  solution, which the mesh cannot be fitted to in advance. These are the
  `plasma_boundary` cells, and they are what caps the order (M-09, M-10).
* **Conductor and vacuum: only under a subtracting conductor model.** A meshed
  coil is fragmented into the mesh by `halfdisc.py`, and its `Check` asserts
  each coil attribute covers its rectangle, so its cells are pure. Under
  `Model = "subtracted"` or `"filament"` the mesh command emits no `--coil` at
  all (`apps/meq.cpp`, `meshCommand`). That is deliberate:
  `COIL-SUBTRACTION-PLAN.md` exists so that the mesh need not resolve the
  copper, since `ψ_c` is evaluated analytically. In the solve those cells are
  vacuum for the remainder equation, so the split costs nothing numerically.
  The fraction is a statement about the drawing, not the current. geq's packs
  (0.07 m wide, in 0.4 m elements) are the case.
* **Plasma and a conductor at one point: a fault, never a mesh question.** It
  means plasma current inside copper. The classifier counts such elements
  (`region_faults`), the driver warns, and the file keeps the fractions.
* **All three in one cell, without sharing a point: a cell to split.** A
  plasma boundary passes through the vacuum part of a cell that holds copper.
  The remedy is a geometry-only `[mesh.generate] AlignConductors`, sketched
  below. A fourth state of the classifier is not the remedy.

**`AlignConductors`, the same idea as `AlignPlasma`.** Fragment the coil
rectangles into the mesh as **geometry only**: no coil attribute, no
`CoilSize` grading, and the remainder equation sees no source there. That
makes the conductor cells pure under subtraction, at the price of conforming
to the rectangle's edges. It is off by default, to keep CS's purpose, and is
measured in ER-0. A run that reports any `mixed` cell names this key.

### 2.2 One implementation, in the solver

**The classification is the solver's own support decision, asked at the
solver's own points, and is not re-derived in the driver.** Two readers of "is
this plasma" is how a written mask comes to disagree with the current that was
integrated. CS-4 and `totalFluxAt()` are the precedents.

* **The new method** is `GradShafranovSolver::elementMaterials()`. It
  returns, per element, the non-zero `( material, fraction )` pairs.
* **The predicate** is exactly `SourceIntegrator`'s:
  `elementCarriesPlasma( e )` (the component mask and the exclusions) and the
  source's `insidePlasma( R, z, ψ_h + ψ_c )` (the one-sided or box support,
  `supportAdmits()`).
* **The points** are the source rule's (`2k + OrderW + extra`), weighted by
  `ip.weight × |J|`.
* **Without a support** (an unconfined source, or no normalisation), plasma is
  every element where the plasma term is evaluated. That is all of `Ω_h`, so
  `plasma_fraction = 1`, which is what such a source integrates.
* **Conductor membership is geometric**: `meq::Coil::contains()` at the same
  points, whatever the conductor model. The cross-section is the copper whether
  its current is meshed, subtracted or carried as filaments. A future passive
  or iron region is a shape with a `contains()`, asked the same way.
* **Each point is assigned to one material.** A conductor wins over vacuum. A
  point that is both plasma and conductor is a fault: it is counted, and goes
  to the conductor.
* **It is computed on the host after the solve**, reading the potential through
  `HostRead()`. The element loop is threaded with caller-allocated
  transformations into `thread_local` (the `GetElementTransformation` trap).

### 2.3 Where it is written

* **`.nc`, element-indexed (exact).**
  * **The element index** is the element's index in `<stem>.mesh`, which the
    file already names (`mesh_file`, `mesh_digest`) and indexes `psi_dof`
    against.
  * **The fractions are a CF contiguous ragged array**, because a machine with
    hundreds of passive filaments makes a dense `element × material` matrix
    hundreds of megabytes, while almost every element has one entry:
    * `element_material_count( element )`, with
      `sample_dimension = "element_entry"`;
    * `element_material( element_entry )`, an index into the material table;
    * `element_fraction( element_entry )`.
  * **Dense, one per element:**
    * `element_region` (byte, CF `flag_values` / `flag_meanings`);
    * `element_plasma_fraction`, the one fraction every reader wants, kept
      out of the ragged array;
    * `element_excluded`;
    * `element_attribute`.
  * **The material table**, on a dimension `material`, carries the columns of
    §2.1. The circuit table is on `circuit` and the permeability curves are on
    `permeability_curve`, as a second ragged array. The tables use NetCDF-4
    strings, which the writer gains.
  * **Global attributes:** `region_elements_<name>`, `material_area_<name>`
    and `region_faults`.
  * **`material` is not CS-6's `conductor` dimension.** That one holds the
    `ConductorField`'s members, expanded to filaments under
    `Model = "filament"` (geq's two packs are twenty there). A coil row here is
    a `[[coils]]` block.
* **`.nc`, on the grid (pointwise).** `node_material( Z, R )` is an int, the
  material table's index at that node, and −1 where `inside` is 0. Plasma is
  judged by the same predicate, `GradShafranovSolver::materialAt()`, on the
  field the grid carries (`ψ* + ψ_c`) and in the element holding the node (its
  foot's, in the band); conductors by `contains()`. No fraction and no boundary
  exist at a point. It is not named `material`, which is the table's
  dimension and would read as its coordinate variable. This is MR-5's F8, and
  `docs/output.rst` points a reader of `n_s` to it.
* **VTK.** The cell data on `<stem>.pvd`, and on each `<stem>_cycles` frame so
  the plasma boundary can be watched under adaptivity, are `region`,
  `plasma_fraction`, `conductor_material` (the dominant conductor's row, −1
  for none) and `conductor_fraction`: `meq::ElementCellData`, L2 order-0
  `GridFunction`s registered beside `psi` and `B`. `ParaViewDataCollection`
  writes them as point data at each cell's own Lagrange points, which no two
  cells share, so each is one constant per cell. `attribute` is the
  collection's own cell data. The full ragged record stays in the `.nc`, since
  VTK has no ragged cell data.
* **Not the `.gf` trio.** Those three files are the exact restart and stay
  byte for byte what they are. The element arrays in the `.nc` are already
  exact against the `.mesh`.

## 3. Stages

Each stage ends at a measured number.

**ER-0: the classification in the library. DONE** →
**[M-266](MEASUREMENTS.md#m-266)**, **[M-267](MEASUREMENTS.md#m-267)**.
`elementMaterials()` takes a `meq::MaterialTable` (`src/meq/Materials.hpp`,
MFEM-free: the kinds, the rows with every §2.1 column, the region rule) and
returns a `meq::ElementMaterials`. `tests/convergence/ElementRegions.cpp`
asserts each bullet below. Three are met in a nearby form, for these reasons:

* The fill and the exclusions use a Solov'ev private-flux lobe, ITER Example 2
  below its X-point, in place of `PlasmaConnectivity`'s solved half-disc lobe.
* The mirror is classified on the vacuum field at MR-4's labels and not on the
  rotating solve. The cutoff is the same line either way.
* "The shipped examples count 0" is measured through the driver's `regions:`
  line, which ER-0 adds with §2.1b's warnings. It is 0 on every example but
  `coils-rectangle`, whose unnormalised source fills a box that holds its coils
  and faults in 32 elements (M-267). That is the problem it poses, so its
  header says so and `theDriverAddsTheCoilsToF` pins the count.

The bullets, each asserted against something the classifier cannot fake:

* **The partition.** In every element the fractions sum to 1 to round-off.
  The plasma area, `Σ plasma_fraction × |e|`, converges to the area of
  `{ ψ > ψ_bnd }` on the exact Solov'ev, computed by fine quadrature of the
  closed form. The `plasma_boundary` element count grows like `1/h`.
* **The conductors.** On a meshed-coil `halfdisc` mesh no cell is
  `conductor_boundary`, and each coil's area equals its rectangle's to
  round-off. **Control:** the same machine under `Model = "subtracted"`,
  meshed without `--coil`, has `conductor_boundary` cells on every rectangle
  edge.
* **`AlignConductors`.** Under subtraction it makes those cells pure. Its
  element count is measured against the coil-free mesh, on geq's mirror and on
  `limited-tokamak-filament`: that count is the cost of the key and decides
  whether it is worth offering.
* **The mirror.** On the aligned mesh no element straddles `|z| = AxialEdge`.
  **Control:** M-262's unaligned mesh has `plasma_boundary` elements on the
  cutoff.
* **The fill and the exclusions.** `PlasmaConnectivity`'s disconnected lobe
  reads `plasma_fraction = 0`; turning the fill off makes it 1. An excluded
  attribute reads 0.
* **The fault.** A plasma that reaches a coil, posed deliberately, is counted
  in `region_faults`. The shipped examples count 0.
* **The predicate is the integrated one.** Summing `F/R` over the points the
  classifier calls plasma reproduces the assembled plasma load. Mutation:
  classify on `ψ_h` without `ψ_c` under a split, and this must go red.

**ER-1: the `.nc`. DONE** → **[M-268](MEASUREMENTS.md#m-268)**.
`NetCDFWriter::materials()` writes the record the driver's `regions:` line
reports, classified once on the converged solve; `meq::buildMaterialTable()`
makes the table from `[[coils]]` and `[[circuits]]`, with the solved currents.
The `DriverAcceptance` bullets are met on the three runs named. The diverted
tokamak is wired in the test, not in `examples/`: `P2L` and `P2U` on one
supply and `P1L` and `P1U` unwired. The mirror's coil areas are not asserted,
since its mesh is not aligned to them. It writes the ragged array, the
material and circuit tables, and the attributes of §2.3. The writer gains
string vectors, and the parser
reserves `default`. `DriverAcceptance` reads them back on `limited-tokamak`,
`diverted-tokamak` (with circuits and unwired coils, so both kinds of circuit
row appear) and `mirror-rotating`:

* the element count equals the `.mesh`'s, and `Σ count` equals the length of
  `element_entry`;
* the material areas sum to the domain's;
* the coil areas match the rectangles where the mesh is aligned;
* every coil's circuit is either its `[[circuits]]` block or `default`;
* the mirror has no fault and a non-empty plasma.

**The reserved kinds must be readable before they are produced**, so a unit
test writes a material table holding a `passive` row and two `ferromagnetic`
rows, one linear and one with a permeability curve, through the writer and
reads them back. It asserts that the curve's `permeability_curve_material`
equals the `material_name` of the row that points to it, and that the curve's
points round-trip exactly. `docs/output.rst` tells a reader to make the same
check. Also `docs/output.rst`, and
`tools/plot_equilibrium.py` drawing the regions.

**ER-2: VTK. DONE** → **[M-269](MEASUREMENTS.md#m-269)**. The cell data on the
answer and on every cycle frame, asserted by reading the `.vtu` back by hand
(`tests/convergence/VtuReader.hpp`) in `OutputTests` and `DriverAcceptance`.
**The high-order output does not smear an L2-0 field**: every element's
Lagrange points are its own, so the point count is the cell count times the
points of a cell, and each array is one value per cell, on curved cells too.
The adaptive frames are the rotating mirror with `[adaptivity]` added, which
holds both boundary kinds. Its estimator refines away from the plasma
boundary, so the frames show the record per cycle without showing that
boundary move.

**ER-3: the grid mask. DONE** → **[M-270](MEASUREMENTS.md#m-270)**.
`node_material( Z, R )`. On the mirror the plasma entries agree with
`0 < Ψ ≤ 1, |z| ≤ AxialEdge`, evaluated independently from the written `psi`
and the labels, and the copper entries with the rectangles, at every node more
than one grid spacing from an edge. No node disagrees anywhere on the shipped
run. The `docs/output.rst` note about `n_s` outside the support points to the
mask.

## 4. Decided

1. **Fractions are stored**, as a ragged array per element. The categories are
   derived from them for plotting.
2. **A cell with plasma and copper at one point is a fault**, counted and
   warned about. A cell holding all three without such a point is a cell to
   split, which is `AlignConductors`' job.
3. **Unwired coils form one `default` circuit**, a grouping and not a supply,
   and `default` is reserved in `[[circuits]]`.
4. **Passive and ferromagnetic kinds are reserved** in the material table and
   its columns, and are readable before MEQ produces them. **Permeability
   curves are stored in the file**, each carrying its material's name. The vessel's
   attribute 30 is `excluded`, not a material.
5. **No new `.gf` files**; the `.nc` and the VTK carry the regions.
