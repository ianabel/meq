"""MEQ from Python: the run `meq` runs, in-process.

    import meq
    eq = meq.Equilibrium("examples/soloviev-nstx.toml")
    eq.solve()
    eq.write()

is `meq examples/soloviev-nstx.toml`, a phase at a time, and writes the same
files. `meq.Equilibrium.from_dict(tables, base_dir=...)` takes the TOML's
tables as a dict instead of a file.

Profiles can be handed in as arrays in place of the files a rotating
source's configuration names: `profiles={"species.deuterium.Temperature":
meq.Table(x, T, dT_dx, abscissa="psi_n"), ...}`. The abscissa is always
named -- "Psi" (MEQ's own normalised flux, 1 at the axis), "psi_n" (1 - Psi),
"sqrt_psi_n", "psi", the poloidal flux itself, which a normalised source
takes with every profile in psi or none, or "R_midplane", a box mirror's
radius at z = 0, mapped through the last answer's crossing
(`eq.midplane_crossing()`) and so only after a first solve -- and the
derivative is always given. The .nc lists them as `profile_tables`.

An Equilibrium solves again after `eq.set_profile(slot, table)`: warm,
restarted exactly from the last state on the same mesh, and a solve that
fails leaves it holding the last good state. A configuration with
[adaptivity] solves once.

The solved state reads back under the .nc's names, without a file:
`eq.attributes()`, `eq.grid()` (the ( Z, R ) lattice), `eq.dataset()` (the
whole .nc), `eq.at(R, Z)` (psi and B at points, the grid's own sampling) and
`eq.surfaces()` or `eq.surfaces(rho)` (the _surfaces.nc, closed surfaces only).

A configuration with [mesh.generate] needs its mesh made first, by `meq-run`
or by the command `meq --mesh-command` prints, and then `mesh_ready=True`, as
`meq --mesh-ready` does.

Failures raise meq.ConfigError, meq.SolveError or meq.OutputError, all
meq.Error, where `meq` would exit 1, 2 or 3; the run's own diagnostics are on
stderr. The report goes to the process's stdout, as `meq`'s does.

IMPORT meq BEFORE numpy, OR ANYTHING ELSE THAT BRINGS AN OpenMP RUNTIME.
meq carries its own MKL and SuiteSparse, static and hidden, so nothing loaded
before or after it can reach them (python/CMakeLists.txt). OpenMP it cannot
carry: its threaded loops call libgomp, and the calls are bound when meq is
imported, to whichever OpenMP the process holds first. Debian's numpy brings
Intel's, through the system MKL, and imported before meq it would run meq's
loops on that runtime -- so that import is refused, naming the library.
Imported first, meq binds to libgomp, and numpy's runtime arriving later
beside it changes nothing of meq's. pip's numpy and JAX on the CPU bring no
MKL and no OpenMP -- their BLAS is OpenBLAS on pthreads, under prefixed names
-- so with those either order is fine; measured, in a venv of this
interpreter, meq's files are the binary's beside a jitted JAX computation.

PYTHON-PLAN.md in the source tree is the design and what comes next.
"""

import os as _os

#: Every OpenMP runtime exports the GOMP_* entry points meq's loops call, so
#: for the guard they are one library whatever each is named.
_OPENMP = ("libgomp", "libiomp5", "libomp", "libiomp")


def _family(path):
    name = _os.path.basename(path).split(".so")[0]
    return "OpenMP" if name in _OPENMP else name


def _foreign_libraries():
    """Each mapped library that is another file of one _meq links.

    By name less its version -- libsundials_core.so and .so.7 are one library,
    and every OpenMP runtime is one -- against the files recorded at build
    time. Linux only, where /proc/self/maps lists what is mapped; elsewhere
    there is nothing to read.
    """
    from ._libraries import LINKED
    try:
        with open("/proc/self/maps") as maps:
            mapped = {line.split()[-1] for line in maps if ".so" in line}
    except OSError:
        return []
    ours = {}
    for path in LINKED.values():
        ours.setdefault(_family(path), set()).add(_os.path.realpath(path))
    clashes = []
    for path in sorted(mapped):
        real = _os.path.realpath(path)
        known = ours.get(_family(path))
        if known and real not in known:
            clashes.append(f"{real} is loaded, and meq was built against "
                           f"{' or '.join(sorted(known))}")
    return clashes


_clashes = _foreign_libraries()
if _clashes:
    raise ImportError(
        "meq: this process already holds another copy of a library meq links, "
        "and meq would bind to it rather than to its own: " + "; ".join(_clashes)
        + ". Import meq before numpy, scipy or anything else that loads an "
          "OpenMP runtime or an MKL.")


from ._meq import (
    ConfigError,
    Equilibrium,
    Error,
    OutputError,
    SolveError,
    Table,
    __version__,
    build_info,
)

__all__ = [
    "ConfigError",
    "Equilibrium",
    "Error",
    "OutputError",
    "SolveError",
    "Table",
    "__version__",
    "build_info",
]
