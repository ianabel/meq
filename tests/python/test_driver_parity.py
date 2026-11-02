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
import pytest

from harness import (EXAMPLES, compare, fresh, run_binary, run_module)

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
