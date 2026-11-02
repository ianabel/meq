"""The machinery the tests/python suites share: directories of their own, the
binary and the module run in child processes, and .nc files compared bit for
bit. Imported by each test file after `import meq`, which must come first
(meq/__init__.py).
"""

import os
import pathlib
import shutil
import subprocess
import sys

import netCDF4
import numpy

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


def run_module(where, construct, preamble=""):
    """meq.Equilibrium in `where` -- prepare, solve, write -- IN A CHILD PROCESS
    THAT IMPORTS meq AND NOTHING ELSE, which is the process PY-1's claim is
    about. This one has numpy and netCDF4 loaded to read the files, and
    Debian's numpy brings a second MKL, with which meq's mixes in the last
    digits (meq/__init__.py). `construct` is the expression that makes the
    Equilibrium. True when write() reports the targets met."""
    script = ("import meq\n"
              + preamble +
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
