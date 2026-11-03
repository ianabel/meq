"""PythonProfiles: species profiles handed in as arrays (PY-2).

PYTHON-PLAN.md, PY-2: a slot set from a file's own knots is the file, to the
bit; the same profile in psi_n is re-knotted exactly; in sqrt( psi_n ) it is
composed, and agrees to its tabulation, at the Hermite cubic's order; and what
a table cannot do is refused by name. Every solve runs in a child process that
imports meq alone, for the reason harness.run_module() gives.
"""

import json
import math
import re
import shutil
import subprocess
import sys
import tomllib

# meq FIRST (meq/__init__.py).
import meq

import netCDF4
import numpy
import pytest

from harness import EXAMPLES, compare, fresh, run_binary, run_module

CASE = "rotating-normalised"
SLOT = "species.D.Density"


def read_table(path):
    """A profile file's knots, as SplineProfile::fromStream reads them:
    `psi f f'` a line, '#' comments, a blank line ends it. float() and the
    C++ reader both round the decimal correctly, so these are its doubles."""
    x, f, d = [], [], []
    for line in open(path):
        line = line.strip()
        if line.startswith("#"):
            continue
        if not line:
            if x:
                break
            continue
        a, b, c = (float(t) for t in line.split()[:3])
        x.append(a), f.append(b), d.append(c)
    return x, f, d


@pytest.fixture(scope="module")
def unscaled():
    """rotating-normalised's tables with the density's DensityScale folded into
    a file of its own, written to round-trip exactly: a table takes no scale,
    and the file and the table must be the same knots for the comparison to be
    about the route and not the rounding of a multiplication."""
    tables = tomllib.loads((EXAMPLES / f"{CASE}.toml").read_text())
    species = tables["source"]["species"][0]
    assert species["Name"] == "D" and "DensityScale" in species
    scale = species.pop("DensityScale")
    x, f, d = read_table(EXAMPLES / species["DensityFile"])
    f = [scale*v for v in f]
    d = [scale*v for v in d]
    where = fresh("unscaled")
    density = where / "density.dat"
    density.write_text("".join(f"{a!r} {b!r} {c!r}\n" for a, b, c in zip(x, f, d)) + "\n")
    species["DensityFile"] = str(density)
    config = where / "config.json"
    config.write_text(json.dumps(tables))
    return config, (x, f, d)


def solve(name, config, profiles=None):
    """One arm in a child process: the dict from `config`, `profiles` as
    {slot: (x, value, derivative, abscissa)}."""
    where = fresh(name)
    preamble = f"import json\ntables = json.load(open({str(config)!r}))\n"
    preamble += "profiles = {}\n"
    for slot, (x, f, d, abscissa) in (profiles or {}).items():
        preamble += (f"profiles[{slot!r}] = meq.Table({list(x)!r}, {list(f)!r}, "
                     f"{list(d)!r}, abscissa={abscissa!r})\n")
    run_module(where, "meq.Equilibrium.from_dict(tables, name='rotating', profiles=profiles)",
               preamble)
    return where


def psi_of(where, stem):
    with netCDF4.Dataset(where / f"{stem}.nc") as nc:
        nc.set_auto_mask(False)
        inside = nc["inside"][...] == 1
        return numpy.asarray(nc["psi"][...])[inside]


def psi(where):
    return psi_of(where, CASE)


def test_a_table_of_a_files_own_knots_is_the_file(unscaled):
    config, (x, f, d) = unscaled
    from_file = solve("own/file", config)
    from_table = solve("own/table", config, {SLOT: (x, f, d, "Psi")})
    compare(from_file, from_table, CASE, skip={"profile_tables"})
    with netCDF4.Dataset(from_table / f"{CASE}.nc") as nc:
        assert nc.getncattr("profile_tables") == f"{SLOT} (Psi, {len(x)} knots)"
    with netCDF4.Dataset(from_file / f"{CASE}.nc") as nc:
        assert "profile_tables" not in nc.ncattrs()


def test_a_table_in_psi_n_is_the_same_profile(unscaled):
    """Knots reversed at 1 - Psi with their slopes negated: re-knotted exactly,
    so only the rounding of 1 - x separates it from the file."""
    config, (x, f, d) = unscaled
    from_file = psi(solve("psi_n/file", config))
    in_psi_n = psi(solve("psi_n/table", config,
                         {SLOT: ([1.0 - a for a in reversed(x)], list(reversed(f)),
                                 [-c for c in reversed(d)], "psi_n")}))
    gap = numpy.max(numpy.abs(in_psi_n - from_file)) / numpy.max(numpy.abs(from_file))
    print(f"\n  psi_n against the file: {gap:.2e} of max |psi|")
    assert gap < 1e-12


def density(psi):
    """n( Psi ) = 1e20 ( 0.5 + 2 Psi + 0.6 sin 2 Psi ) and its slope: smooth,
    and no polynomial, so a resampled table has a tabulation error to fall."""
    return (1e20*(0.5 + 2.0*psi + 0.6*math.sin(2.0*psi)),
            1e20*(2.0 + 1.2*math.cos(2.0*psi)))


def test_a_table_in_sqrt_psi_n_agrees_to_its_tabulation(unscaled):
    """density() as a fine Psi table against it in rho = sqrt( psi_n ) on 8,
    16 and 32 intervals, each with its exact slope in rho,
    dn/drho = -2 rho n'( Psi ), and re-knotted in Psi by MEQ: the answers
    close on the fine one at THIRD order, because F carries the profile's
    slope and a Hermite cubic tabulates a slope at third order. Measured
    from 8 to 256 intervals: 2.45e-6 to 8.4e-11, 2.97 a halving on average
    and 3.5x to 16x a step, so the order is asserted over two halvings.

    THE TWO KINDS OF TABLE MUST POSE THE SAME PROBLEM OUTSIDE THE PLASMA, so
    the rho tables run to rho = sqrt( 1.2 ), Psi = -0.2, as the fine one does:
    this source is not confined to the plasma and its solution dips below
    psi = 0 near the wall, and clamping the profile at Psi = 0 instead moves
    the answer 5.7e-4. PAST THE AXIS a rho table cannot reach, and MEQ
    continues it linearly, so the fine reference is continued the same way;
    and that the answer barely reads past the axis is asserted on its own,
    against the same fine table carrying the profile itself there."""
    config, _ = unscaled
    inside = [-0.2 + 1.2*i/480 for i in range(481)]
    axis_value, axis_slope = density(1.0)
    continued = {SLOT: (inside + [1.5], [density(p)[0] for p in inside]
                        + [axis_value + 0.5*axis_slope],
                        [density(p)[1] for p in inside] + [axis_slope], "Psi")}
    analytic_past = [-0.2 + 1.7*i/680 for i in range(681)]
    itself = {SLOT: (analytic_past, [density(p)[0] for p in analytic_past],
                     [density(p)[1] for p in analytic_past], "Psi")}
    reference = psi(solve("sqrt/continued", config, continued))
    past = psi(solve("sqrt/itself", config, itself))
    beyond = numpy.max(numpy.abs(past - reference)) / numpy.max(numpy.abs(reference))

    gaps = []
    for intervals in (8, 16, 32):
        rho = [i / intervals for i in range(math.ceil(intervals*math.sqrt(1.2)) + 1)]
        values = [density(1.0 - r*r)[0] for r in rho]
        slopes = [-2.0*r*density(1.0 - r*r)[1] for r in rho]
        answer = psi(solve(f"sqrt/{intervals}", config,
                           {SLOT: (rho, values, slopes, "sqrt_psi_n")}))
        gaps.append(numpy.max(numpy.abs(answer - reference)) / numpy.max(numpy.abs(reference)))
    print("\n  sqrt_psi_n against the fine Psi table, 8/16/32 intervals: "
          + ", ".join(f"{g:.2e}" for g in gaps)
          + f"; the profile past the axis moves the answer {beyond:.2e}")
    assert beyond < 1e-7, "the equilibrium reads the profile past the axis"
    assert gaps[2] < 1e-7
    assert gaps[0] / gaps[2] > 2.0**5, \
        "the tabulation error does not fall at third order over two halvings"


# ---- what a table cannot do ------------------------------------------------------

def test_a_table_says_its_abscissa_and_its_lengths_agree():
    with pytest.raises(TypeError):
        meq.Table([0.0, 1.0], [1.0, 1.0], [0.0, 0.0])
    with pytest.raises(ValueError, match="no abscissa"):
        meq.Table([0.0, 1.0], [1.0, 1.0], [0.0, 0.0], abscissa="rho")
    with pytest.raises(ValueError, match="lengths"):
        meq.Table([0.0, 1.0], [1.0, 1.0], [0.0], abscissa="Psi")
    table = meq.Table([0.0, 1.0], [1.0, 2.0], [1.0, 1.0], abscissa="psi_n")
    assert table.abscissa == "psi_n" and table.x == [0.0, 1.0] and table.value == [1.0, 2.0]


@pytest.mark.parametrize("slot", ["species.D.Densty", "species.e.Density",
                                  "species.D.Temperature", "PPrime"])
def test_a_slot_the_source_cannot_take_is_refused(slot, unscaled):
    """A misspelt slot; the Neutralising species' density; a slot the file
    scales (TemperatureScale); a slot of another source type."""
    config, _ = unscaled
    tables = json.loads(config.read_text())
    flat = meq.Table([0.0, 1.0], [1.0, 1.0], [0.0, 0.0], abscissa="Psi")
    with pytest.raises(meq.ConfigError):
        meq.Equilibrium.from_dict(tables, profiles={slot: flat})


def test_a_source_without_species_takes_no_table():
    tables = tomllib.loads((EXAMPLES / "mhd-rectangle.toml").read_text())
    flat = meq.Table([0.0, 1.0], [1.0, 1.0], [0.0, 0.0], abscissa="Psi")
    with pytest.raises(meq.ConfigError):
        meq.Equilibrium.from_dict(tables, base_dir=str(EXAMPLES), profiles={"PPrime": flat})


# ---- tables in psi itself (PY-2b) ----------------------------------------------

def newton_steps(report):
    found = re.search(r"converged in (\d+) Newton iterations", report.read_text())
    assert found, f"no convergence line in {report}"
    return int(found.group(1))


MIRROR = EXAMPLES / "mirror-rotating.toml"


@pytest.fixture(scope="module")
def mirror_from_file():
    """mirror-rotating from its files, through meq-run: its mesh and its answer."""
    where = fresh("mirror/file")
    assert run_binary(MIRROR, where, wrapper=True) == 0
    return where


def mirror_slots():
    """mirror-rotating's four profile slots, and the file each reads."""
    source = tomllib.loads(MIRROR.read_text())["source"]
    species = {one["Name"]: one for one in source["species"]}
    return {"Omega": source["OmegaFile"],
            "species.electron.Temperature": species["electron"]["TemperatureFile"],
            "species.electron.Density": species["electron"]["DensityFile"],
            "species.deuterium.Temperature": species["deuterium"]["TemperatureFile"]}


def with_mesh_of(made, name):
    """A fresh directory holding the mesh `made` generated."""
    where = fresh(name)
    for one in made.iterdir():
        if one.suffix in (".msh", ".meq-mesh"):
            shutil.copy2(one, where / one.name)
    return where


def test_a_mirror_in_psi_is_the_mirror(mirror_from_file):
    """mirror-rotating from its files, through meq-run, and then with every
    profile -- both temperatures, the electron density and omega -- as a table
    in psi, mapped from the file's at the labels that run converged to: the
    same equilibrium, in a solve whose labels float, with Newton taking as many
    steps. At fixed psi the tables do not move with the labels, and the
    source reports their derivatives as zero (SourceFactoryTests)."""
    toml = MIRROR
    from_file = mirror_from_file
    with netCDF4.Dataset(from_file / "mirror-rotating.nc") as nc:
        axis, boundary = nc.getncattr("psi_axis"), nc.getncattr("psi_boundary")
    span = axis - boundary

    def in_psi(name):
        x, f, d = read_table(EXAMPLES / name)
        knots = sorted(zip((boundary + span*a for a in x), f, (c/span for c in d)))
        return [k[0] for k in knots], [k[1] for k in knots], [k[2] for k in knots]

    profiles = {slot: (*in_psi(name), "psi") for slot, name in mirror_slots().items()}

    in_flux = with_mesh_of(from_file, "mirror/psi")
    (in_flux / "tables.json").write_text(json.dumps(profiles))
    preamble = ("import json\n"
                "profiles = {slot: meq.Table(x, f, d, abscissa=a) for slot, (x, f, d, a)"
                " in json.load(open('tables.json')).items()}\n")
    run_module(in_flux, f"meq.Equilibrium({str(toml)!r}, mesh_ready=True, profiles=profiles)",
               preamble)

    with netCDF4.Dataset(in_flux / "mirror-rotating.nc") as nc:
        assert abs(nc.getncattr("psi_axis") - axis) < 1e-12*abs(axis)
        assert abs(nc.getncattr("psi_boundary") - boundary) < 1e-12*abs(boundary)
    a, b = psi_of(from_file, "mirror-rotating"), psi_of(in_flux, "mirror-rotating")
    gap = numpy.max(numpy.abs(a - b)) / numpy.max(numpy.abs(a))
    steps = newton_steps(from_file / "report.txt"), newton_steps(in_flux / "report.txt")
    print(f"\n  the mirror in psi against its files: {gap:.2e} of max |psi|, "
          f"Newton {steps[1]} steps against {steps[0]}")
    assert gap < 1e-11
    assert abs(steps[1] - steps[0]) <= 1, "Newton takes more steps on tables in psi"


def test_a_table_in_psi_refuses_a_midplane_gauge():
    """mirror-three-species is in the midplane gauge, which tabulates its
    crossings in Psi; a table in psi there is refused, by name."""
    flat = "meq.Table([0.0, 1.0], [1.0, 1.0], [0.0, 0.0], abscissa='psi')"
    where = fresh("midplane")
    done = subprocess.run(
        [sys.executable, "-c",
         "import meq\n"
         f"meq.Equilibrium({str(EXAMPLES / 'mirror-three-species.toml')!r}, mesh_ready=True,"
         f" profiles={{'Omega': {flat}}})\n"],
        cwd=where, capture_output=True, text=True)
    assert done.returncode != 0
    assert "ConfigError" in done.stderr and "midplane gauge" in done.stderr, done.stderr[-800:]


# ---- tables in R at the midplane (PY-2b) -----------------------------------------

def hermite(x, f, d, at):
    """The Hermite cubic through knots (x, f, d) at `at`, and its slope: the
    crossing midplane_crossing() returns, as the solver interpolates it."""
    i = max(0, min(len(x) - 2, int(numpy.searchsorted(x, at)) - 1))
    h = x[i + 1] - x[i]
    t = (at - x[i]) / h
    value = ((2*t**3 - 3*t**2 + 1)*f[i] + (t**3 - 2*t**2 + t)*h*d[i]
             + (-2*t**3 + 3*t**2)*f[i + 1] + (t**3 - t**2)*h*d[i + 1])
    slope = ((6*t**2 - 6*t)*f[i] + (3*t**2 - 4*t + 1)*h*d[i]
             + (-6*t**2 + 6*t)*f[i + 1] + (3*t**2 - 2*t)*h*d[i + 1]) / h
    return value, slope


def in_r(crossing, name, scale=1.0):
    """A profile file's knots in Psi re-expressed in R through `crossing`, its
    own knots mapped and its slopes divided by dR/dPsi there, `scale` times."""
    x, f, d = read_table(EXAMPLES / name)
    knots = []
    for a, b, c in zip(x, f, d):
        if 0.0 <= a <= 1.0:
            r, dr = hermite(crossing["Psi"], crossing["R"], crossing["dR_dPsi"], a)
            knots.append((r, scale*b, scale*c/dr))
    knots.sort()
    return meq.Table([k[0] for k in knots], [k[1] for k in knots], [k[2] for k in knots],
                     abscissa="R_midplane")


def test_a_mirror_in_R_midplane_is_held_through_its_crossing(mirror_from_file, monkeypatch,
                                                             capfd):
    """A table in R_midplane is mapped into Psi through the crossing of the
    last answer, held for the solve. Before any solve there is none, and it is
    refused. After the solve from the files, the four profiles re-expressed in
    R through its crossing re-solve to that answer in at most one Newton step,
    and again in none. With Omega 10% up in R, each re-solve moves the crossing
    and the next solve follows it: the steps between successive answers shrink
    by more than ten each time, to rounding -- the held map is a fixed-point
    iteration, and on this mirror it contracts."""
    where = with_mesh_of(mirror_from_file, "mirror/R_midplane")
    monkeypatch.chdir(where)
    slots = mirror_slots()

    early = meq.Equilibrium(str(MIRROR), mesh_ready=True)
    flat = meq.Table([0.1, 0.4], [1.0, 1.0], [0.0, 0.0], abscissa="R_midplane")
    early.set_profile("Omega", flat)
    capfd.readouterr()
    with pytest.raises(meq.ConfigError):
        early.solve()
    assert "no answer before the first solve" in capfd.readouterr().err

    reference = psi_of(mirror_from_file, "mirror-rotating")
    scale = numpy.max(numpy.abs(reference))

    def solved(eq):
        eq.solve()
        eq.write()
        with netCDF4.Dataset(where / "mirror-rotating.nc") as nc:
            steps = nc.getncattr("newton_iterations")
        return psi_of(where, "mirror-rotating"), steps

    eq = meq.Equilibrium(str(MIRROR), mesh_ready=True)
    first, _ = solved(eq)
    assert numpy.max(numpy.abs(first - reference)) / scale < 1e-12
    crossing = eq.midplane_crossing()
    assert crossing["R"][0] > crossing["R"][-1], "the crossing is not the mirror's, outward at Psi 0"

    for slot, name in slots.items():
        eq.set_profile(slot, in_r(crossing, name))
    for allowed in (1, 0):
        now, steps = solved(eq)
        gap = numpy.max(numpy.abs(now - reference)) / scale
        print(f"\n  in R through its own crossing: {gap:.2e} of max |psi|, Newton {steps}")
        assert gap < 1e-11
        assert steps <= allowed

    eq.set_profile("Omega", in_r(crossing, slots["Omega"], scale=1.1))
    last, _ = solved(eq)
    moves = []
    for _ in range(8):
        now, _ = solved(eq)
        moves.append(numpy.max(numpy.abs(now - last)) / scale)
        last = now
        if moves[-1] < 1e-13:
            break
    print("  Omega 10% up in R, each re-solve's move: " + ", ".join(f"{m:.1e}" for m in moves))
    assert moves[-1] < 1e-13, "the held map did not settle in eight re-solves"
    for before, after in zip(moves, moves[1:]):
        assert after < before / 10 or after < 1e-13, "the held map contracts by less than ten"


def test_a_table_in_R_midplane_takes_a_box(monkeypatch, capfd):
    """R_midplane names one crossing of each surface, a box mirror's; a
    tokamak's surfaces cross z = 0 twice, and a re-solve with one is refused
    by its key. midplane_crossing() refuses too."""
    monkeypatch.chdir(fresh("R_midplane/tokamak"))
    eq = meq.Equilibrium(str(EXAMPLES / f"{CASE}.toml"))
    eq.solve()
    with pytest.raises(meq.ConfigError, match="not a box"):
        eq.midplane_crossing()
    eq.set_profile("Omega", meq.Table([0.5, 1.5], [1.0, 1.0], [0.0, 0.0], abscissa="R_midplane"))
    capfd.readouterr()
    with pytest.raises(meq.ConfigError):
        eq.solve()
    assert "key 'source.SupportShape'" in capfd.readouterr().err
