"""PythonResults: the solved state in memory, under the .nc's names (PY-4).

PYTHON-PLAN.md, PY-4: `eq.attributes()`, `eq.grid()` and `eq.dataset()` are
the .nc's root group, recorded by the calls that write the file
(meq::NetCDFWriter::dataset()) rather than sampled a second time. Asserted:
what an Equilibrium holds in memory before any write() is, bit for bit, every
dimension, variable and attribute of the .nc its write() then writes; grid()
is the ( Z, R ) lattice of it; nothing is held before a solve; and after a
re-solve, or a failed one, what is held is the state the Equilibrium holds.
In this process, which imports meq first (meq/__init__.py).
"""

import os
import shutil
import tomllib

import meq

import netCDF4
import numpy
import pytest

from harness import EXAMPLES, STAMPS, fresh, run_binary

CASE = "rotating-normalised"
TOML = str(EXAMPLES / f"{CASE}.toml")


def same(mine, theirs):
    """Bit for bit and of the same kind: NaN equal to NaN, an int attribute
    not a double one, a byte variable not an int one."""
    if isinstance(mine, str) or isinstance(theirs, str):
        return mine == theirs
    if isinstance(mine, (int, float)):
        # A Python scalar for a scalar attribute: its kind is the file's.
        kind = numpy.asarray(theirs).dtype.kind
        if kind != ("f" if isinstance(mine, float) else "i") or numpy.ndim(theirs) != 0:
            return False
        return mine == theirs or (mine != mine and theirs != theirs)
    if isinstance(mine, list) or (isinstance(theirs, numpy.ndarray) and theirs.dtype == object):
        return list(mine) == list(numpy.asarray(theirs).ravel())
    a, b = numpy.asarray(mine), numpy.asarray(theirs)
    if a.shape != b.shape or a.dtype.kind != b.dtype.kind or a.dtype.itemsize != b.dtype.itemsize:
        return False
    return numpy.array_equal(a, b, equal_nan=a.dtype.kind == "f")


def against_file(memory, path):
    """Every way the dict `memory` (eq.dataset()) differs from the root group
    of the file at `path`."""
    out = []
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_maskandscale(False)
        dims = {name: len(d) for name, d in nc.dimensions.items()}
        if dims != memory["dimensions"]:
            out.append(f"dimensions {memory['dimensions']} against {dims}")
        attributes = memory["attributes"]
        for name in sorted(set(nc.ncattrs()) | set(attributes)):
            if name in STAMPS:
                continue
            if name not in nc.ncattrs() or name not in attributes:
                out.append(f":{name} in one only")
            elif not same(attributes[name], nc.getncattr(name)):
                out.append(f":{name} {attributes[name]!r} against {nc.getncattr(name)!r}")
        variables = memory["variables"]
        for name in sorted(set(nc.variables) | set(variables)):
            if name not in nc.variables or name not in variables:
                out.append(f"{name} in one only")
                continue
            v, mine = nc.variables[name], variables[name]
            if tuple(v.dimensions) != mine["dimensions"]:
                out.append(f"{name} on {mine['dimensions']} against {v.dimensions}")
            elif not same(mine["values"], v[...] if v.dimensions else v.getValue()):
                out.append(f"{name} differs")
            for attr in sorted(set(v.ncattrs()) | set(mine["attributes"])):
                if attr not in v.ncattrs() or attr not in mine["attributes"] \
                        or not same(mine["attributes"][attr], v.getncattr(attr)):
                    out.append(f"{name}:{attr} differs")
    return out


@pytest.fixture
def where(monkeypatch):
    """Run in a fresh directory of the test's own, where write() writes."""
    def into(name):
        path = fresh(name)
        monkeypatch.chdir(path)
        return path
    return into


def test_the_answer_in_memory_is_the_file(where):
    """rotating-normalised: dataset() before any write() -- held from the state
    file a re-solve would restart from -- is the .nc write() then writes,
    every dimension, variable and attribute, to the bit."""
    here = where("file")
    eq = meq.Equilibrium(TOML)
    eq.solve()
    memory = eq.dataset()
    assert eq.write()
    found = against_file(memory, here / f"{CASE}.nc")
    assert not found, "; ".join(found[:12])
    assert len(memory["variables"]) > 10


@pytest.fixture(scope="module")
def mirror_mesh():
    """mirror-rotating's mesh, made once by meq-run."""
    made = fresh("mirror/mesh")
    assert run_binary(EXAMPLES / "mirror-rotating.toml", made, wrapper=True) == 0
    return made


def with_mesh(made, here):
    for one in made.iterdir():
        if one.suffix in (".msh", ".meq-mesh"):
            shutil.copy2(one, here / one.name)


def test_a_mirror_in_memory_is_its_file(where, mirror_mesh):
    """mirror-rotating, whose .nc carries the species densities and the gauge
    besides: the same, after write() this time."""
    here = where("mirror/module")
    with_mesh(mirror_mesh, here)
    eq = meq.Equilibrium(str(EXAMPLES / "mirror-rotating.toml"), mesh_ready=True)
    eq.solve()
    assert eq.write()
    memory = eq.dataset()
    found = against_file(memory, here / "mirror-rotating.nc")
    assert not found, "; ".join(found[:12])
    assert {"n_deuterium", "n_electron"} <= set(eq.grid())


def test_the_grid_is_the_lattice(where):
    """grid() is every variable on ( Z, R ), and R and Z: shaped ( Z, R ), with
    NaN wherever `inside` is 0, as the file has it."""
    where("grid")
    eq = meq.Equilibrium(TOML)
    eq.solve()
    grid, memory = eq.grid(), eq.dataset()
    lattice = {name for name, v in memory["variables"].items()
               if v["dimensions"] in (("Z", "R"), ("R",), ("Z",))}
    assert set(grid) == lattice
    assert {"R", "Z", "inside", "extrapolated", "psi", "B_R", "B_Z"} <= lattice
    shape = (len(grid["Z"]), len(grid["R"]))
    assert grid["psi"].shape == shape and grid["inside"].dtype == numpy.int8
    outside = grid["inside"] == 0
    assert outside.any() or grid["inside"].all()
    assert numpy.isnan(grid["psi"][outside]).all()
    assert not numpy.isnan(grid["psi"][~outside]).any()
    attributes = eq.attributes()
    assert attributes == memory["attributes"]
    assert isinstance(attributes["newton_iterations"], int)
    assert isinstance(attributes["psi_axis"], float)


def test_nothing_is_held_before_a_solve(where):
    where("before")
    eq = meq.Equilibrium(TOML)
    for read in (eq.attributes, eq.grid, eq.dataset):
        with pytest.raises(meq.Error, match="solve\\(\\) first"):
            read()
    assert not list(os.listdir("."))


def omega(value):
    return meq.Table([-0.5, 1.5], [value, value], [0.0, 0.0], abscissa="Psi")


def test_what_is_held_is_the_state(where):
    """After a re-solve the results are the new state's, and the file its
    write() writes; after a re-solve that fails they are the last good
    state's, unchanged."""
    here = where("state")
    eq = meq.Equilibrium(TOML)
    eq.solve()
    first = eq.dataset()

    eq.set_profile("Omega", omega(4.4e5))
    eq.solve()
    second = eq.dataset()
    assert not numpy.array_equal(first["variables"]["psi"]["values"],
                                 second["variables"]["psi"]["values"], equal_nan=True)
    assert second["attributes"]["warm_start"] != first["attributes"].get("warm_start")
    assert eq.write()
    found = against_file(second, here / f"{CASE}.nc")
    assert not found, "; ".join(found[:12])

    eq.set_profile("Omega", omega(3.0e6))
    with pytest.raises(meq.SolveError):
        eq.solve()
    after = eq.dataset()
    assert numpy.array_equal(after["variables"]["psi"]["values"],
                             second["variables"]["psi"]["values"], equal_nan=True)
    assert after["attributes"] == second["attributes"]


# ---- fields at points --------------------------------------------------------

@pytest.mark.parametrize("case", [CASE, "miller-curved", "limited-tokamak-filament"])
def test_at_the_lattice_nodes_the_points_are_the_grid(case, where):
    """at(R, Z) at every node of grid()'s lattice is grid() to the bit: psi,
    B_R, B_Z and both masks. miller-curved has a band between Gamma_h and
    Gamma, which at() continues as the .nc does; limited-tokamak-filament
    adds psi_c back under [conductors]."""
    toml = EXAMPLES / f"{case}.toml"
    if case == "limited-tokamak-filament":
        made = fresh(f"{case}/mesh")
        assert run_binary(toml, made, wrapper=True) == 0
        here = where(f"{case}/module")
        for one in made.iterdir():
            if one.suffix in (".msh", ".meq-mesh"):
                shutil.copy2(one, here / one.name)
        eq = meq.Equilibrium(str(toml), mesh_ready=True)
    else:
        where(case)
        eq = meq.Equilibrium(str(toml))
    eq.solve()
    grid = eq.grid()
    R, Z = numpy.meshgrid(grid["R"], grid["Z"])
    found = eq.at(R, Z)
    for name in ("psi", "B_R", "B_Z", "inside", "extrapolated"):
        assert found[name].shape == grid[name].shape
        assert found[name].dtype == grid[name].dtype
        assert numpy.array_equal(found[name], grid[name], equal_nan=True), name
    if case == "miller-curved":
        assert grid["extrapolated"].any(), "no band node: the case does not test the band"


def test_points_off_the_lattice_are_the_field(where):
    """Between the nodes, at() is the solved field: outside the domain NaN with
    inside 0, and inside it B_R = -(1/R) dpsi/dZ and B_Z = (1/R) dpsi/dR to
    the discretisation's error -- B is the solved flux q relabelled and psi
    is psi*, two fields that agree to the solve's accuracy and not to
    rounding: 1.1e-05 and 1.8e-04 of |B| here."""
    where("off")
    eq = meq.Equilibrium(TOML)
    eq.solve()
    # BEFORE anything is written, and so before the write phase has made psi*:
    # at() makes it itself, and reads what it reads after.
    mesh = tomllib.loads(open(TOML).read())["mesh"]
    early_R = numpy.linspace(mesh["RMin"], mesh["RMax"], 7)[1:-1]
    early_Z = numpy.full_like(early_R, 0.5 * (mesh["ZMin"] + mesh["ZMax"]))
    early = eq.at(early_R, early_Z)
    assert numpy.isfinite(early["psi"]).all() and (early["psi"] != 0).all()
    grid = eq.grid()
    again = eq.at(early_R, early_Z)
    for name in early:
        assert numpy.array_equal(early[name], again[name], equal_nan=True), name
    rows, cols = numpy.nonzero((grid["inside"] == 1) & (grid["extrapolated"] == 0))
    middle = len(rows) // 2
    R0 = 0.5 * (grid["R"][cols[middle]] + grid["R"][cols[middle] + 1])
    Z0 = 0.5 * (grid["Z"][rows[middle]] + grid["Z"][rows[middle] + 1])
    h = 1e-5
    R = numpy.array([R0, R0 + h, R0 - h, R0, R0, 50.0])
    Z = numpy.array([Z0, Z0, Z0, Z0 + h, Z0 - h, 0.0])
    found = eq.at(R, Z)
    assert list(found["inside"]) == [1, 1, 1, 1, 1, 0]
    assert numpy.isnan(found["psi"][-1]) and numpy.isnan(found["B_R"][-1])
    psi = found["psi"]
    dpsi_dR = (psi[1] - psi[2]) / (2 * h)
    dpsi_dZ = (psi[3] - psi[4]) / (2 * h)
    B = numpy.hypot(found["B_R"][0], found["B_Z"][0])
    print(f"\n  B from psi's differences: {abs(-dpsi_dZ / R0 - found['B_R'][0]) / B:.1e}, "
          f"{abs(dpsi_dR / R0 - found['B_Z'][0]) / B:.1e} of |B|")
    assert abs(-dpsi_dZ / R0 - found["B_R"][0]) < 1e-3 * B
    assert abs(dpsi_dR / R0 - found["B_Z"][0]) < 1e-3 * B
    with pytest.raises(meq.Error, match="one shape"):
        eq.at(numpy.zeros(3), numpy.zeros(2))


# ---- flux surfaces -------------------------------------------------------------

SURFACES = "flux-surfaces"


def test_the_surfaces_in_memory_are_the_file(where):
    """flux-surfaces, whose [output] writes a _surfaces.nc: surfaces() before
    any write() -- extracted then, into a private file -- is the _surfaces.nc
    write() then writes, every dimension, variable and attribute."""
    here = where("surfaces/file")
    eq = meq.Equilibrium(str(EXAMPLES / f"{SURFACES}.toml"))
    eq.solve()
    memory = eq.surfaces()
    assert eq.write()
    written = list(here.rglob("*_surfaces.nc"))
    assert len(written) == 1, written
    found = against_file(memory, written[0])
    assert not found, "; ".join(found[:12])


def test_the_surfaces_at_labels_are_the_familys(where):
    """surfaces(rho) at the family's own labels is its columns, to rounding;
    between them its monotone cubic; outside the cut refused."""
    where("surfaces/labels")
    eq = meq.Equilibrium(str(EXAMPLES / f"{SURFACES}.toml"))
    eq.solve()
    family = eq.surfaces()["variables"]
    rho = family["rho"]["values"]
    at = eq.surfaces(rho)
    for name in ("V_prime", "volume", "cross_section_area", "arc_length", "surface_area",
                 "inverse_R_squared", "grad_psi_squared_over_R_squared", "abs_grad_psi",
                 "grad_psi_squared", "psi", "normalised_flux"):
        column = family[name]["values"]
        gap = numpy.max(numpy.abs(at[name] - column)) / numpy.max(numpy.abs(column))
        assert gap < 1e-14, f"{name}: {gap:.2e}"
    assert set(at) >= {"rho", "safety_factor"} or "safety_factor" not in family
    between = 0.5 * (rho[1:] + rho[:-1])
    volume = eq.surfaces(between)["volume"]
    assert numpy.all(volume > family["volume"]["values"][:-1])
    assert numpy.all(volume < family["volume"]["values"][1:])
    assert eq.surfaces(numpy.array([[rho[2], rho[3]]]))["volume"].shape == (1, 2)
    with pytest.raises(meq.Error, match="outside the family's cut"):
        eq.surfaces([0.5 * rho[0]])


def test_a_box_has_no_surfaces(where, mirror_mesh):
    here = where("surfaces/mirror")
    with_mesh(mirror_mesh, here)
    eq = meq.Equilibrium(str(EXAMPLES / "mirror-rotating.toml"), mesh_ready=True)
    eq.solve()
    with pytest.raises(meq.ConfigError, match="open"):
        eq.surfaces()
    with pytest.raises(meq.ConfigError, match="open"):
        eq.surfaces([0.5])
