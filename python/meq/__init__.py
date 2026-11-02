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
named -- "Psi" (MEQ's own normalised flux, 1 at the axis), "psi_n" (1 - Psi)
or "sqrt_psi_n" -- and the derivative is always given. The .nc lists them as
`profile_tables`.

A configuration with [mesh.generate] needs its mesh made first, by `meq-run`
or by the command `meq --mesh-command` prints, and then `mesh_ready=True`, as
`meq --mesh-ready` does.

Failures raise meq.ConfigError, meq.SolveError or meq.OutputError, all
meq.Error, where `meq` would exit 1, 2 or 3; the run's own diagnostics are on
stderr. The report goes to the process's stdout, as `meq`'s does.

IMPORT meq BEFORE numpy, OR ANYTHING ELSE THAT BRINGS AN MKL. Debian's numpy
reaches the system MKL through libmkl_rt, which loads MKL's layers into the
process's global scope, and anything loaded after that binds dgemm_, PARDISO
and the rest of MKL to it: MEQ, UMFPACK and CHOLMOD would run on an MKL they
were not built against. That import is refused, naming the library.

Imported first, meq's own calls stay bound to its own MKL; but an MKL that
arrives later still mixes with it inside MKL, and the answers then move in
their last digits (about 1e-13 relative, measured), so the run is no longer
`meq`'s to the bit. meq warns when that has happened. One MKL for every
module in the process is the fix: PYTHON-PLAN.md, PY-6.

PYTHON-PLAN.md in the source tree is the design and what comes next.
"""

import os as _os
import warnings as _warnings


def _foreign_libraries():
    """Each mapped library that is another file of one _meq links.

    By name less its version -- libmkl_core.so and libmkl_core.so.3 are one
    library -- against the files recorded at build time. Linux only, where
    /proc/self/maps lists what is mapped; elsewhere there is nothing to read.
    """
    from ._libraries import LINKED
    try:
        with open("/proc/self/maps") as maps:
            mapped = {line.split()[-1] for line in maps if ".so" in line}
    except OSError:
        return []
    family = lambda path: _os.path.basename(path).split(".so")[0]
    ours = {}
    for path in LINKED.values():
        ours.setdefault(family(path), set()).add(_os.path.realpath(path))
    clashes = []
    for path in sorted(mapped):
        real = _os.path.realpath(path)
        known = ours.get(family(path))
        if known and real not in known:
            clashes.append(f"{real} is loaded, and meq was built against "
                           f"{' or '.join(sorted(known))}")
    return clashes


_clashes = _foreign_libraries()
if _clashes:
    raise ImportError(
        "meq: this process already holds another copy of a library meq links, "
        "and meq would bind to it rather than to its own: " + "; ".join(_clashes)
        + ". Import meq before numpy, scipy or anything else that loads an MKL, "
          "or build them against the same library.")

_warned = False


def _warn_if_foreign_libraries():
    """Called by each phase: warn, once, if another copy has arrived since."""
    global _warned
    if _warned:
        return
    clashes = _foreign_libraries()
    if clashes:
        _warned = True
        _warnings.warn(
            "meq: another copy of a library meq links has been loaded since meq "
            "was imported: " + "; ".join(clashes) + ". meq's own calls stay bound "
            "to its own, but the two mix inside it, and this run's answers are "
            "no longer meq's to the last digits. One MKL for every module in the "
            "process is the fix (PYTHON-PLAN.md, PY-6).", RuntimeWarning, stacklevel=3)


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
