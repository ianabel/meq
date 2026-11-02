# MEQ from Python: a third front end, for transport workflows

**The goal.** A Python process calls MEQ in-process, hands it profiles as
arrays, gets an equilibrium and its geometry back as arrays, and does it again
with the next profiles. **The interface is general, not MaNTA's**: any Python
caller that has species profiles can drive it. MaNTA is its first user, and
for MaNTA it is the LOOSE coupling of `MANTA-COUPLING.md` (its table, the left
column): an outer fixed-point loop in Python, with no derivatives and no
`FieldModel`. It is the file-based `TODO` entry's in-memory twin, and it is
what the tight coupling will later be validated against (`MANTA-COUPLING.md`,
*Which coupling is this?*).

**What it is not.** No `dGeometry_dpsi`, no `A_2`, no `FieldModel` plugin, no
adjoint. Nothing here closes an item of `MANTA-COUPLING.md` §13 except the
geometry served cheaply, which `meq::GeometryCache` already does.

**Where it runs.** MEQ and MaNTA are compiled on the same machine, and both
extensions are built for **the system Python**, 3.14.4 here (headers in
`/usr/include/python3.14`), so that one interpreter imports both.

## The design rule: a front end, not a second driver

**`meq`, `meq-convert` and the Python module are three front ends on ONE run.**
The run is `meq::Run` (`src/meq/Run.cpp`, about 8,100 lines), in four phases
that share state: configuration, setup, solve and refine, and write. Every
refusal, the reactive ladder, the retry rung, the label checks, the warm
start, the exterior-source diagnostic and the outputs live there and nowhere
else. `apps/meq.cpp` is the command line over it.

**A binding over the library pieces (`Configuration`, `makeSource`,
`GradShafranovSolver`) would be a second driver**, and two drivers can each
converge, to different equilibria, which is exactly what M-264 asserts the
driver and the library do not do. It would also lose every refusal the driver
makes on the user's behalf. So the first stage, PY-0, moved the run out of
`main()`, and the module binds that.

**Python must reproduce `meq` to the bit** on the same configuration. That is
PY-1's acceptance, and it is why the rule exists: a race or a coupling run
through Python is then a run through MEQ's front end
(`races-go-through-user-front-ends`).

## The Python surface

```python
import meq

eq = meq.Equilibrium("examples/mirror-rotating.toml")        # the file, as `meq` reads it
eq = meq.Equilibrium.from_dict(cfg, base_dir="examples/")     # or a dict, MaNTA's PyToml route

# a species profile: values and their derivatives in the abscissa named
eq.set_profile("species.deuterium.Temperature", x, T, dT_dx, abscissa="R_midplane")
eq.set_profile("species.electron.Density", x, n, dn_dx, abscissa="sqrt_psi_n")
eq.set_profile("Omega", x, omega, domega_dx, abscissa="psi")

r = eq.solve()            # raises meq.ConfigError (exit 1) or meq.SolveError (exit 2)
r.psi_inner, r.psi_boundary, r.newton_steps, r.longest_crawl, r.exterior_fraction

eq.grid()                 # the .nc's (R, Z) fields as numpy, its names and masks
eq.psi(R, Z), eq.B(R, Z)  # point evaluation at arrays of points
eq.surfaces(rho)          # the _surfaces.nc quantities at these labels (closed surfaces)
eq.write("stem")          # what `meq` would have written, [output] as configured
```

Seven decisions in it, each with its reason:

1. **The configuration is MEQ's TOML, whichever way it arrives.** A dict
   becomes a `toml::value` by MaNTA's `tomlFromPy()` (`../../manta/src/PyToml.hpp`:
   `bool` before `int`, an integral float stays a float, a numpy scalar is
   unwrapped through `item()`), and goes to the existing
   `Configuration( toml::value const &, std::string const & )`. One schema, one
   parser, one set of refusals. `base_dir` does for a dict what a file's own
   directory does in `resolveInputPaths()`.
2. **Profiles are species profiles, on every machine.** Arrays go only into a
   `[source] Type = "rotating"` source: `species.<Name>.Density` and
   `species.<Name>.Temperature` for each species, `Omega`, and `GGPrime`. With
   `Omega` absent the source is the static multi-species equation
   (`Config.cpp`: both keys absent means `ω = 0`), so a tokamak takes the same
   route as a mirror, and `p = Σ n_s T_s` is closed in MEQ by quasineutrality,
   as it is from files. **A slot is named by the key that names its file, less
   `File`.** Refused: a slot the configuration did not give (a `Neutralising`
   species has no density), and arrays for any other source type, which keep
   working from files. Once set, the slot holds those knots in place of the
   file, and the file is never read.
3. **The abscissa is named every time a profile is passed, and MEQ does the
   mapping.** Four choices:

   | `abscissa` | is | moves with the solve |
   |---|---|---|
   | `"psi"` | poloidal flux, Wb/rad | yes, through the labels |
   | `"psi_n"` | normalised poloidal flux, 0 at the axis or the inner label, 1 at `ψ_bnd`: geq's `ψ_n`, MEQ's `1 − Ψ` | no |
   | `"sqrt_psi_n"` | `√ψ_n` | no |
   | `"R_midplane"` | the radius, m, where the surface crosses `z = 0` | yes, through the crossing |

   **A name, not a default**, because the `ψ_n`/`Ψ` swap is a known silent
   trap: tables read in the wrong one converge, at full order, to an answer
   8.56e-04 away (M-263's mutation). **The table stays a Hermite spline in its
   own variable and is composed with `x( Ψ )`**, with `df/dΨ = f′( x ) x′( Ψ )`
   evaluated at the quadrature point. It is never re-knotted in `Ψ`: `√ψ_n`
   has an infinite slope at the axis, and a cubic in `Ψ` cannot carry it.
   **The derivative is required, in `x`.** A missing `f′` is refused rather
   than differenced, because the Hermite table's accuracy is the derivative's
   (the `TODO` entry's question 1: MaNTA's `q` is that derivative, for free).
4. **`solve()` after the first warm-starts from the last converged state, on
   the same mesh.** It takes the whole iterate `( q, ψ, ψ̂ )`, labels and
   exterior coefficients, exactly as a `.nc` restart on its own mesh does.
   **The mesh is frozen after the first solve unless `solve( adapt = True )`**:
   an outer fixed-point loop converges to a fixed point of a fixed discrete
   map, and a mesh that moves every outer iteration makes the loop's residual
   the mesh's.
5. **A failed solve raises and leaves the object as it was**, so that an outer
   loop can damp and retry. The exception carries the failed attempt's report.
   "Throw, never fudge" is MaNTA's own contract (`MANTA-COUPLING.md` §8), and a
   Picard driver relies on it more than IDA does. Exit 4, an inverse best fit,
   is a status on the result, not an exception, as it is a status on `meq`.
6. **Results carry the `.nc`'s names, never new ones.** `grid()` is the `.nc`'s
   uniform grid with `inside` and `extrapolated`. `surfaces()` takes the
   `_surfaces.nc` variables. A Python user who later reads a file finds the
   same words.
7. **What a density means is the gauge's, and the result says which gauge.**
   Under the fixed-radius gauge a `Density` is the density at
   `ReferenceRadius`. Under the midplane gauge it is the density where the
   surface crosses `z = 0`, so `R_midplane` with the midplane gauge is exactly
   a midplane profile. That is the M-167 trap (*compare `n_s( R, z )`, never
   `N_s( ψ )`*), and it is why the gauge goes back to the caller with the
   answer.

**The GIL is held for the whole solve, and the module does not declare
`mod_gil_not_used`.** This is MaNTA's reasoning (`Python.cpp`, the comment
above `PYBIND11_MODULE`) and it applies here unchanged: MEQ has statics nobody
has audited for a free-threaded caller, and OpenMP threads are not Python
threads, so holding the GIL costs a sequential outer loop nothing.
**Threads come from the environment exactly as for `meq`**: `MKL_NUM_THREADS`,
`OMP_NUM_THREADS`, and the performance keys in the TOML. The module adds no
knob, for the reason `CLAUDE.md` gives about `Globalisation`.

## Stages

Each stage ends at a measured, asserted result, as every stage here does.

| stage | what | acceptance |
|---|---|---|
| PY-0 | `meq::Run`, the run out of `main()` | `meq` unchanged, to the byte — **done** |
| PY-1 | the `_meq` extension and the `meq` package; `Equilibrium( file )`, `from_dict`, `solve`, `write` | Python's `.nc` equals `meq`'s, every coefficient — **done** |
| PY-2 | species profiles from arrays, on the fixed abscissae `psi_n` and `sqrt_psi_n` | a slot set from a file's knots equals the file, to the bit |
| PY-2b | the moving abscissae, `psi` and `R_midplane` | Newton stays quadratic; a moving table equals the fixed table it maps to at the answer |
| PY-3 | repeated solves: warm, frozen mesh, transactional failure | the three properties below, asserted |
| PY-4 | results in memory: scalars, `grid()`, `psi`/`B` at points, `surfaces()` | equal to the `.nc` and `_surfaces.nc` the same run writes |
| PY-5 | a mirror's transport geometry | closed-form and finite-difference checks, below |
| PY-6 | MEQ and MaNTA in one process, then the first outer loop | each code's answer in the joint process equals its answer alone |

**PY-0: `meq::Run`.** The four phases are methods of a private
`meq::Run::Impl` in `src/meq/Run.cpp` — `configure()`, `setUp()`, `solve()`
and `write()` — and the locals they share are its members, declared in
`main()`'s order so that they are destroyed in it. `main()` keeps `argv`,
`--help`, `--version`, `--mesh-command` and the exit code, and hands the rest
to `meq::RunOptions`. The bodies moved as text. Every edit to them is one of:
a declaration hoisted to a member; the two clocks and three type aliases
moved to namespace scope; the two lambdas the solve calls (`makeSolver`,
`totalFluxAt`) made members; one by-name capture made an init-capture of the
same value; the one `std::exit` made a throw that `execute()` turns into the
same code; and `MEQ_VERSION` handed in by the front end, whose link stamps the
commit. **The report still goes to stdout and stderr**, as `meq`'s always has;
capturing it for Python is PY-1's.

**This was the largest stage and it is behaviour-free by construction**, held
to that by four checks. The hoisted names: `-Wshadow -Winit-self` gives the
same five warnings, on the same names, as the original, and no member name is
used before the line `main()` declared it on, so none now resolves to
something else. The whole suite must pass as before. And the `.nc`, exit code
and masked report of every shipped example must equal the ones the installed
binary made before the move.

**PY-1: the extension.** Laid out as MaNTA's, and built:

* `extern/pybind11` is a submodule **pinned to MaNTA's commit, `b0fa039`**
  (v3.0.4, `stable`), for the reason toml11 is pinned to MaNTA's: two modules
  loaded into one interpreter should share a pybind11 internals ABI.
* `python/meq/` holds the package, whose `__init__.py` re-exports `_meq`;
  `python/_meq.cpp` is the extension, built by
  `pybind11_add_module( _meq NO_EXTRAS ... )` -- **NO_EXTRAS for MaNTA's
  reasons**: no second LTO, and no strip, because a backtrace from inside MKL
  needs symbols. It is staged into `builds/<name>/python/meq/` by symlink and
  installed to `install/<name>/lib/python3.14/site-packages/meq`.
* **The interpreter is the system `python3`**, found once in
  `python/CMakeLists.txt`, with the headers, the ABI suffix and the test's
  interpreter derived from it. No venv is preferred over it.
* `MEQ_PYTHON` defaults to ON beside `MEQ_BUILD_APP` and needs
  `extern/pybind11` populated; CI, which builds no MFEM, never reaches it.
* **`meq.Equilibrium( path )` is `meq path` a phase at a time** --
  construction is `meq::Run::prepare()`, then `solve()` and `write()`, which
  return `True` or `False` where an inverse solve exits 0 or 4 -- and
  `Equilibrium.from_dict( tables, base_dir=, name= )` takes the TOML's tables:
  `tomlFromPy()`, MaNTA's, makes them TOML text, which `meq::RunOptions`
  carries as `configText` beside `baseDirectory`, and which the `.nc` records
  as `input_toml`. `mesh_ready`, `grid_functions`, `vtk` and `profile` are
  the command line's flags. **The device is not an argument**: `mfem::Device`
  is process-wide and configured once, which one equilibrium among several in
  a process cannot own.
* Failures raise `meq.ConfigError`, `meq.SolveError` and `meq.OutputError`,
  all `meq.Error`, where `meq` exits 1, 2 and 3; a phase out of order, again,
  or after a failed one raises `RuntimeError` (`meq::Run`'s `logic_error`).
* **The report goes to the process's file descriptors 1 and 2**, as `meq`'s
  does, with Python's buffered stdout flushed before each phase so the two
  interleave in order. A notebook does not capture it; capturing it is still
  owed, and it is a sink in `meq::Run` beside its 220 `printf`s.
* **No stub yet**: `pybind11-stubgen` is not in the system Python, so the
  `stubs-check` target MaNTA has is not built. Owed with the interface's first
  change after this one.
* `meq.build_info()` reports version, precision and backends, as `meq
  --version` and the `.nc` do; `meq/_libraries.py`, written after every link
  by `python/record_libraries.py` from `ldd`, records the files the module
  resolves to.

**THE MODULE NEEDS POSITION-INDEPENDENT MFEM AND meq_core, AND HAS BOTH.**
`libmfem.a` and `meq_core` are static and the module is a shared object. The
compiler's default `-fPIE` is position independent for an executable only:
MFEM's objects reach libstdc++'s data PC-relatively, and `meq_core`'s
`thread_local`s use a TLS model a `dlopen`ed object cannot. So
`../../mfem/build` is configured `-DCMAKE_POSITION_INDEPENDENT_CODE=ON`
(`CLAUDE_MFEM.md`'s table), `python/CMakeLists.txt` gives `meq_core` the same
with `-fno-semantic-interposition`, and the module links with
`--exclude-libs,ALL`, which makes MFEM's and `meq_core`'s symbols local to it
and out of the process's global scope.

**AND IT IS THE RUN `meq` RUNS, TO THE BIT, IN A PROCESS THAT LOADS NOTHING
ELSE OF ITS LIBRARIES.** `PythonDriverParity` runs `soloviev-nstx`,
`mhd-rectangle`, `limited-tokamak` (free boundary, coils, a `.gf` guess) and
`mirror-rotating` (a generated mesh, the box, species) through the build's
`meq` and through the module in a child process that imports `meq` alone,
and requires every variable and attribute of the two `.nc` files equal but
`meq_version`; then `from_dict()` against the file it was read from, a
relative table without its `base_dir` refused, each exception against its
exit code, and the phase order.

**DEBIAN'S numpy BRINGS A SECOND MKL, AND THAT IS PY-6's HAZARD MET IN
PY-1.** The system numpy reaches `/usr/lib/x86_64-linux-gnu/libmkl_rt.so`,
the 2020-era system MKL, through the BLAS alternative, and `libmkl_rt` loads
MKL's layers into the global scope. Measured, with `LD_DEBUG=bindings`:

* **numpy imported first**, `_meq`, UMFPACK and CHOLMOD bind `dgemm_`,
  `dgetrf_`, `dpotrf_` and PARDISO to the system MKL, which MEQ was not built
  against. `import meq` refuses it, naming the two files.
* **`meq` imported first**, its own calls stay bound to oneAPI 2026.2, but
  MKL's lazily loaded kernels mix with the second MKL: `soloviev-nstx`'s
  `psi` moves 7.5e-14 of its 0.44 peak and `limited-tokamak`'s 3.0e-14, with
  the same Newton counts. Neither a warm-up of oneAPI's kernels before numpy
  arrives nor MKL's `MKL_CBWR=COMPATIBLE` brings the bits back (both tried),
  and at one thread it is still there, so it is not OpenMP's reduction order.
  meq warns, once, at the next phase.

One MKL in the process is the fix, and it is PY-6's: the bits are not wrong
by more than round-off, but they are not `meq`'s.

**PY-2: species profiles from arrays, on the fixed maps.** One C++ change:
`makeSource()` and `makeNormalisedSource()` take an optional
`meq::ProfileTables`, a map from slot to a table and its abscissa, consulted
before the file. A file's table is the abscissa `Ψ`, which is the identity
map, so **the file path and the array path build the same `SplineProfile`
from the same knots**. That is the first acceptance: a slot set from
`SplineProfile::fromFile()`'s own knots equals the file run to the bit. The
second: a table in `psi_n` and the same function in `Ψ` agree to the
tabulation. The third: on a tokamak, the species route at `ω = 0` reproduces
the `mhd` source with `p′` from `Σ n_s T_s`, to the tabulation. Refused by
name, each with a test: an unknown slot, a slot the configuration did not
give, a source type with no species, a non-increasing or non-finite abscissa,
a missing derivative, and a length mismatch. Setting a slot rebuilds the
source and hands it to the solver again on the next `solve()`. Nothing else is
rebuilt.

**PY-2b: the moving maps.** `psi` makes `x` depend on the labels, which are
already border unknowns, so it adds terms to the label columns the bordered
Newton carries. `R_midplane` makes `x` depend on the surface's crossing, which
the midplane gauge already carries as border unknowns, `R₀` at 32 knots in
`Ψ` (`GradShafranovSolver::setMidplaneGaugeKnots()`, M-283). The abscissa
reuses those knots rather than building a second crossing table.
**`R_midplane` is refused beside a located axis**: a tokamak surface crosses
`z = 0` twice, and which crossing is meant is a choice no default should make.
The acceptance is that Newton stays quadratic, as M-283's was, and that a
table in a moving abscissa, re-expressed in `psi_n` through the answer's own
map, gives the same answer to the tabulation. **If quadratic Newton is not
reached, the fallback holds the map at the warm start's state and says so**:
a held map is a fixed point of the outer loop and not of the solve, which is
the smaller claim.

**PY-3: repeated solves.** Three properties, each asserted:

* **Unchanged profiles re-solve in at most one Newton step, to the same
  answer.**
* **Perturbed profiles re-solve, warm, to the cold solve on the same mesh**,
  within the solver's tolerance, and in fewer steps. Both counts are reported.
  A warm solve that crawls is a defect, not a speed-up that is smaller than
  hoped (`crawling-is-a-bug`).
* **A solve that fails leaves the object bit-identical to before it**: the
  iterate, the labels, the source and the mesh. Asserted by solving a profile
  set the mesh cannot carry (GF1 from the vacuum on level 0 is one, `TODO`)
  and then the original again.

**PY-4: results in memory.** The scalars are the `.nc`'s attributes, with the
gauge among them. `grid()` is assembled by the code that samples the `.nc`, not
by a second sampler. `psi( R, Z )` and `B( R, Z )` go through
`meq::PointLocator`, with points outside `D_h` refused or masked, matching the
`.nc`'s `extrapolated` rule. `surfaces( rho )` is `GeometryCache` behind the
`_surfaces.nc` names, and refused beside a box support, where the field lines
are open. The acceptance is equality with the files the same run writes.

**PY-5: a mirror's transport geometry.** This is new physics output, not
binding work, and a mirror transport case will ask for it first. The quantities
MaNTA's mirror physics read before field models existed (`MANTA-COUPLING.md`
§5) were `ψ(V)`, `B(V, z)`, `R(ψ, s)`, `dR/dV`, `V′`, the mirror ratio, `L_∥`,
`R_max` and `R_min`. Each is served against any of the four abscissae.
Open surfaces are traced and fitted already (IN-5:
`ContourTracer::traceOpen()`, `meq::fitOpenSurface()`). **Built in C++,
written to the `.nc` and served to Python by one implementation**, so the
binary and the module cannot disagree. Checked against the vacuum mirror's
closed-form field where one exists, and against a finite difference of `V( ψ )`
for `V′`.

**PY-6: two codes in one process, then the loop.** Two hazards, measured on
2026-11-02 by `ldd`:

* **SUNDIALS by SONAME.** MaNTA's `_manta` links `libsundials_core.so.7`,
  `kinsol.so.7` and `ida.so.7` from `../../sundials/install`; MEQ links the
  same names from `../../sundials/cuda-install`. The loader binds a SONAME
  once per process, so whichever module is imported first serves both. If the
  builds differ in `sunrealtype` or index size, the second module runs on a
  library it was not built for, silently. A single build of MEQ makes that
  likely.
* **Two MKLs.** `_manta` links the system `libmkl_rt.so`; MEQ links oneAPI's
  `latest` `lp64`/`gnu_thread`/`core`. Two MKL versions in one process are not
  supported by Intel, and `mkl-2026-2-moves-single-precision` says this
  machine's answers can move with the MKL version. **Measured in PY-1, and
  numpy alone is enough**: Debian's numpy is on the same system `libmkl_rt`,
  and with it in the process MEQ's answers move about 1e-13 relative (PY-1,
  above). So the one MKL has to serve numpy too.

**The fix is in the build systems, not in the module: one SUNDIALS and one MKL
serve both codes.** The candidate is MEQ's pair, `cuda-install` and oneAPI,
and whether MaNTA builds against them is being checked on MaNTA's side. What
MEQ builds, so that a mismatch can never be silent:

* **At configure**, `-DMEQ_MANTA_PREFIX=` reads that MaNTA's `manta.pc` and
  fails if its SUNDIALS or MKL is not the one MEQ resolved, naming both.
* **At import**, `meq` reads `/proc/self/maps`, and a library already mapped
  under the name of one the module links, from another file, refuses the
  import, naming both files; one that arrives after the import is warned
  about at the next phase. **Built in PY-1** (`meq/__init__.py`, against
  `meq/_libraries.py`). A refusal at import costs nothing, where the wrong
  library costs a wrong answer.

The acceptance is that each code's answer in the joint process, imported in
either order, equals its answer alone, to the bit.

**Then the loop.** A small Python fixed-point iteration between a MaNTA case
and MEQ on a mirror: MaNTA's `n_s`, `T_s` and `ω` go into MEQ's species slots,
on whichever abscissa MaNTA's case uses, and PY-5's geometry comes back. It
reports the outer contraction per iteration and the damping it needed. That
number is what says whether the loose coupling is a usable workflow, or only
the validation tool for the tight one.

## Deliberately not in it

* **Derivatives of any kind.** That is the point of the loose coupling.
* **Releasing the GIL**, until someone does the audit MaNTA's comment asks for.
* **Interrupting a solve with Ctrl-C.** It would need a per-Newton-step hook
  in the solver calling `PyErr_CheckSignals()`. Worth adding if a solve long
  enough to want interrupting turns up; GF1 on level 2 takes seconds.
* **A Python API for meshing.** `[mesh.generate]` runs outside `meq`
  (`meq-run`), because `meq` links MFEM and not gmsh. The module takes a mesh
  file as `meq` does, and a workflow runs `meq-run`'s generator first, or
  `meq --mesh-command`.
* **Arrays for the `mhd`, `soloviev` and `manufactured` sources.** They keep
  working from files; species are the one array route.
* **Writing a source or a field model in Python.** MaNTA cannot do that either
  (`MANTA-COUPLING.md` §11).

## Open

1. **Whether MaNTA builds against oneAPI MKL and `sundials/cuda-install`**
   (PY-6). Being checked on MaNTA's side. If it cannot, PY-6 needs another
   common pair.
2. **Which abscissa MaNTA's mirror case uses.** It no longer changes the
   design, since all four are served. It sets which of PY-2 and PY-2b the
   first loop needs.
