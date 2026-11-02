"""PythonDriverParity: the meq Python module is the run `meq` runs.

PY-1's acceptance (PYTHON-PLAN.md). Each case runs one configuration twice,
once through the installed-alike `meq` binary of this build and once through
`meq.Equilibrium` in this process, each in a directory of its own, and requires
the two .nc files to agree in every variable and every attribute but the
version stamp, bit for bit. Then what only the module has: from_dict(), the
exceptions in place of exit codes, and the phase order.

Run by ctest with PYTHONPATH at this build's package and MEQ_DRIVER_PATH,
MEQ_RUN_PATH and MEQ_SOURCE_DIR in the environment.
"""

import os
import pathlib
import shutil
import subprocess
import sys
import tomllib

# meq FIRST, before numpy and netCDF4, as its docstring asks: numpy brings the
# system MKL into the global scope, and meq imported after it is refused --
# test_meq_refuses_a_process_holding_another_mkl says so.
import meq

import netCDF4
import numpy
import pytest

SOURCE = pathlib.Path(os.environ["MEQ_SOURCE_DIR"])
EXAMPLES = SOURCE / "examples"
DRIVER = os.environ["MEQ_DRIVER_PATH"]
WRAPPER = os.environ["MEQ_RUN_PATH"]
WORK = pathlib.Path.cwd()

# THE ATTRIBUTES THAT MAY DIFFER, and why each may. meq_version: the binary and
# the module are linked separately and each stamps its own commit, which is
# the same commit with or without "-dirty" and nothing else.
STAMPS = {"meq_version"}


def fresh(name):
    """An empty directory for one arm of one case, in the test's own."""
    path = WORK / name
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def run_binary(toml, where, wrapper=False):
    """`meq toml` (or `meq-run toml`) in `where`, and its exit code."""
    command = [WRAPPER if wrapper else DRIVER, str(toml)]
    done = subprocess.run(command, cwd=where, capture_output=True, text=True)
    (where / "report.txt").write_text(done.stdout + done.stderr)
    return done.returncode


def run_module(where, construct):
    """meq.Equilibrium in `where` -- prepare, solve, write -- IN A CHILD PROCESS
    THAT IMPORTS meq AND NOTHING ELSE, which is the process PY-1's claim is
    about. This one has numpy and netCDF4 loaded to read the files, and
    Debian's numpy brings a second MKL, with which meq's mixes in the last
    digits (meq/__init__.py). `construct` is the expression that makes the
    Equilibrium. True when write() reports the targets met."""
    script = ("import meq\n"
              f"eq = {construct}\n"
              "eq.solve()\n"
              "print('MET' if eq.write() else 'UNMET')\n")
    done = subprocess.run([sys.executable, "-c", script], cwd=where,
                          capture_output=True, text=True)
    (where / "report.txt").write_text(done.stdout + done.stderr)
    assert done.returncode == 0, \
        f"the module did not run; its report is {where / 'report.txt'}:\n{done.stderr[-1500:]}"
    return "MET" in done.stdout.split()


def same_value(a, b):
    """Bit for bit, with NaN equal to NaN and masked equal to masked."""
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    a = numpy.ma.asarray(a)
    b = numpy.ma.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if not numpy.array_equal(numpy.ma.getmaskarray(a), numpy.ma.getmaskarray(b)):
        return False
    da = numpy.ma.getdata(a)
    db = numpy.ma.getdata(b)
    if da.dtype.kind == "f":
        return numpy.array_equal(da, db, equal_nan=True)
    return numpy.array_equal(da, db)


def differences(a, b, skip=frozenset(), path="/"):
    """Every attribute and variable of two netCDF groups that differs."""
    out = []
    for name in sorted(set(a.ncattrs()) | set(b.ncattrs())):
        if name in STAMPS or name in skip:
            continue
        if name not in a.ncattrs() or name not in b.ncattrs():
            out.append(f"{path}:{name} present in one file only")
        elif not same_value(a.getncattr(name), b.getncattr(name)):
            out.append(f"{path}:{name} differs")
    for name in sorted(set(a.variables) | set(b.variables)):
        if name not in a.variables or name not in b.variables:
            out.append(f"{path}{name} present in one file only")
            continue
        va, vb = a.variables[name], b.variables[name]
        va.set_auto_maskandscale(False)
        vb.set_auto_maskandscale(False)
        if va.dimensions != vb.dimensions:
            out.append(f"{path}{name} has dimensions {va.dimensions} against {vb.dimensions}")
        elif not same_value(va[...], vb[...]):
            out.append(f"{path}{name} differs")
        for attr in sorted(set(va.ncattrs()) | set(vb.ncattrs())):
            if attr not in va.ncattrs() or attr not in vb.ncattrs() \
                    or not same_value(va.getncattr(attr), vb.getncattr(attr)):
                out.append(f"{path}{name}:{attr} differs")
    for name in sorted(set(a.groups) | set(b.groups)):
        if name not in a.groups or name not in b.groups:
            out.append(f"{path}{name}/ present in one file only")
        else:
            out += differences(a.groups[name], b.groups[name], skip, f"{path}{name}/")
    return out


def compare(binary_dir, module_dir, stem, skip=frozenset()):
    a = binary_dir / f"{stem}.nc"
    b = module_dir / f"{stem}.nc"
    assert a.exists(), f"meq wrote no {a.name}; its report is {binary_dir / 'report.txt'}"
    assert b.exists(), f"the module wrote no {b.name}"
    with netCDF4.Dataset(a) as fa, netCDF4.Dataset(b) as fb:
        found = differences(fa, fb, skip)
    assert not found, (f"{stem}: the module's .nc is not meq's -- " + "; ".join(found[:12])
                       + (f"; and {len(found) - 12} more" if len(found) > 12 else ""))


# ---- the binary and the module, file for file ----------------------------------

FILE_CASES = [
    # name, generated mesh
    ("soloviev-nstx", False),        # fixed boundary, linear
    ("mhd-rectangle", False),        # profile tables, Newton
    ("limited-tokamak", False),      # free boundary: coils, limiter, exterior, a .gf guess
    ("mirror-rotating", True),       # a [mesh.generate] mirror, species, the box
]


@pytest.mark.parametrize("name,generated", FILE_CASES, ids=[c[0] for c in FILE_CASES])
def test_the_module_writes_what_meq_writes(name, generated):
    toml = EXAMPLES / f"{name}.toml"
    binary_dir = fresh(f"{name}/meq")
    module_dir = fresh(f"{name}/module")

    assert run_binary(toml, binary_dir, wrapper=generated) == 0, \
        f"meq did not solve {name}; its report is {binary_dir / 'report.txt'}"

    if generated:
        # THE MESH IS A BUILD PRODUCT meq-run MADE, and the module, like meq,
        # cannot make one: it is handed the same file and told it is ready.
        for made in binary_dir.iterdir():
            if made.suffix in (".msh", ".meq-mesh"):
                shutil.copy2(made, module_dir / made.name)

    met = run_module(module_dir, f"meq.Equilibrium({str(toml)!r}, mesh_ready={generated})")
    assert met is True
    compare(binary_dir, module_dir, name)


# ---- what only the module has --------------------------------------------------

def test_a_dict_is_the_file_it_was_read_from():
    """from_dict() on the file's own tables, with base_dir its directory, is
    the file -- except input_toml, which is the dict written out as TOML
    rather than the file's text, and the path the file was given by."""
    toml = EXAMPLES / "mhd-rectangle.toml"
    binary_dir = fresh("from_dict/meq")
    module_dir = fresh("from_dict/module")
    assert run_binary(toml, binary_dir) == 0

    tables = tomllib.loads(toml.read_text())
    run_module(module_dir,
               f"meq.Equilibrium.from_dict(__import__('tomllib').loads(open({str(toml)!r}).read()), "
               f"base_dir={str(EXAMPLES)!r}, name={str(toml)!r})")
    compare(binary_dir, module_dir, "mhd-rectangle", skip={"input_toml"})

    with netCDF4.Dataset(module_dir / "mhd-rectangle.nc") as written:
        assert tomllib.loads(written.getncattr("input_toml")) == tables, \
            "the .nc's input_toml does not read back as the dict that was run"


def test_a_relative_table_without_its_directory_is_refused():
    """With no base_dir a dict's relative paths are the working directory's,
    where mhd-rectangle's profile tables are not: refused as the binary
    refuses a file whose tables are missing."""
    tables = tomllib.loads((EXAMPLES / "mhd-rectangle.toml").read_text())
    where = fresh("relative")
    previous = os.getcwd()
    os.chdir(where)
    try:
        with pytest.raises(meq.ConfigError):
            meq.Equilibrium.from_dict(tables)
    finally:
        os.chdir(previous)


def test_a_refused_configuration_raises_config_error():
    tables = tomllib.loads((EXAMPLES / "soloviev-nstx.toml").read_text())
    tables["solver"]["NoSuchKey"] = 1
    toml = fresh("refused/meq") / "refused.toml"
    toml.write_text((EXAMPLES / "soloviev-nstx.toml").read_text()
                    .replace("[solver]\n", "[solver]\nNoSuchKey = 1\n"))
    assert run_binary(toml, toml.parent) == 1
    with pytest.raises(meq.ConfigError):
        meq.Equilibrium.from_dict(tables)
    assert issubclass(meq.ConfigError, meq.Error)


def test_a_solve_that_does_not_converge_raises_solve_error():
    """mhd-rectangle capped at two Newton steps and held to a tolerance no
    solve reaches, so the retry ladder cannot rescue it either: meq exits 2,
    and the module raises meq.SolveError, refuses to write, and writes
    nothing. A cap of two, not one: at one KINSOL's Anderson depth clamps to
    zero and the Picard fallback dereferences workspace it never allocated
    (TODO)."""
    text = (EXAMPLES / "mhd-rectangle.toml").read_text()
    for old, new in (("NewtonMaxIterations = 20\n", "NewtonMaxIterations = 2\n"),
                     ("NewtonRelativeTolerance = 1.0e-10\n", "NewtonRelativeTolerance = 1.0e-30\n"),
                     ("NewtonAbsoluteTolerance = 1.0e-14\n", "NewtonAbsoluteTolerance = 1.0e-30\n")):
        assert old in text, f"mhd-rectangle.toml no longer says {old.strip()}"
        text = text.replace(old, new)

    # The capped file beside copies of its two tables, so both arms read
    # them from their own directory.
    binary_dir = fresh("unconverged/meq")
    for table in ("mhd-pprime.dat", "mhd-ggprime.dat"):
        shutil.copy2(EXAMPLES / table, binary_dir / table)
    capped = binary_dir / "capped.toml"
    capped.write_text(text)
    assert run_binary(capped, binary_dir) == 2
    assert not list(binary_dir.glob("*.nc")), "meq wrote output from a failed solve"

    module_dir = fresh("unconverged/module")
    previous = os.getcwd()
    os.chdir(module_dir)
    try:
        eq = meq.Equilibrium.from_dict(tomllib.loads(text), base_dir=str(binary_dir))
        with pytest.raises(meq.SolveError):
            eq.solve()
        with pytest.raises(RuntimeError, match="earlier phase"):
            eq.write()
    finally:
        os.chdir(previous)
    assert not list(module_dir.glob("*.nc")), "the module wrote output from a failed solve"


def test_meq_refuses_a_process_holding_another_mkl():
    """numpy imported first brings the system MKL into the process's global
    scope; meq would bind BLAS to it, so the import is refused and says why.
    In a child process, since this one imported meq first on purpose. Skipped
    where numpy's BLAS is not an MKL meq does not link, since then there is
    nothing to refuse."""
    probe = subprocess.run(
        [sys.executable, "-c",
         "import numpy, os\n"
         "maps = open('/proc/self/maps').read().split()\n"
         "print(sorted({os.path.realpath(p) for p in maps if 'libmkl_core' in p}))"],
        capture_output=True, text=True, check=True)
    numpy_mkl = eval(probe.stdout)
    from meq._libraries import LINKED
    ours = {path for path in LINKED.values() if "libmkl_core" in path}
    if not numpy_mkl or set(numpy_mkl) <= ours:
        pytest.skip("numpy here brings no other MKL")

    done = subprocess.run([sys.executable, "-c", "import numpy; import meq"],
                          capture_output=True, text=True)
    assert done.returncode != 0, "meq imported after numpy's MKL and did not refuse"
    assert "Import meq before numpy" in done.stderr, done.stderr[-600:]


def test_the_phases_run_once_each_in_order():
    where = fresh("order")
    previous = os.getcwd()
    os.chdir(where)
    try:
        eq = meq.Equilibrium(str(EXAMPLES / "soloviev-nstx.toml"))
        with pytest.raises(RuntimeError, match="out of order"):
            eq.write()
        eq.solve()
        with pytest.raises(RuntimeError, match="out of order"):
            eq.solve()
        assert eq.write() is True
        with pytest.raises(RuntimeError, match="out of order"):
            eq.write()
    finally:
        os.chdir(previous)


def test_the_build_is_described():
    info = meq.build_info()
    assert info["version"] == meq.__version__
    assert info["precision"] in ("double", "single")
    assert info["description"]
