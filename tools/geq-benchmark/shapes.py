#!/usr/bin/env python3
"""THE PROFILE-SHAPE CASES: one mirror, four profile shapes, each written as a
geq user's script and a MEQ user's configuration.

Each case directory, tools/geq-benchmark/cases/<shape>/, holds

  geq_<shape>.py   a geq script in the shape of geq's own tools/run.py: the
                   machine, the profiles and the solve stated in code, and only
                   geq's and freegs4e's interface called -- nothing is patched.
                   geq has no input file; this is how a geq user poses a case.
  meq.toml         the same problem as MEQ's configuration, read by `meq-run`.
  *.dat            MEQ's profile tables, in Psi = 1 - psi_n with their slopes.

race_geq.py runs the two through their own front ends, changing only the lines
a user would change to refine a run. This script is the ONE place the cases are
defined: it writes all four directories, and `--check` asserts the tracked files
are what it would write, which race_geq.py does before it times anything.

THE PROBLEM, every case. geq's mirror machine -- two packs at z = -1.5 and
+1.5, field 4 T, inner radius 0.7855 m, length 1 m, ten filaments each -- with
electrons and deuterons, the labels at the electrode ( 0.10, 0 ) and the
limiter ( 0.35, 0 ), and g = 0. THE ELECTRON DENSITY IS GIVEN WHERE EACH FLUX
SURFACE CROSSES z = 0 and the deuterons follow by quasineutrality: geq's own
statement (SingleSpecies density_Z0, AllSpecies.set_flux_density), and MEQ's
Gauge = "midplane" with the deuteron Neutralising.

Every profile is f_edge + ( f_peak - f_edge ) g( psi_n ), with ONE shape g per
case, 1 at its highest point and 0 at both labels, for T_e, T_D, n_e( z = 0 ) and
omega. THE SOURCE MUST VANISH AT BOTH LABELS (docs/rotation.rst, the box's
precondition). g'( 0 ) = g'( 1 ) = 0 makes every profile flat there, but under
the midplane statement a flat density still carries -n m omega^2 R_0 R_0'/T in
dp/dpsi wherever omega does not vanish, so THE DENSITY FALLS TO 1e-6 OF ITS PEAK
at both labels -- geq's own habit, tools/run.py's edge being 1e14 against 5e19 --
which leaves the source there about 1e-7 of its peak.

OMEGA DOES NOT VANISH AT THE LABELS, AND CANNOT HERE: geq's first phi_0 guess
(AllSpecies.solve_quasineutrality_global_no_phi0) is a shape times
m_e omega^2 R^2/2e, which vanishes with omega. Measured on the skew at 129^2,
omega's edge at 0.2 and 0.01 of its peak solves in 13 and 14 Picard steps, and at
1e-4 raises "Can't find a phi0 guess" before the first. So omega's edge is 0.1
of its peak.

  flattop  ( 1 - u^4 )^2, u = 2 psi_n - 1. Broad and flat across the middle,
           symmetric; j = 1 at both labels. A parabola 1 - u^2 has slope at the
           labels, so it cannot be a profile here: this keeps its breadth.
  quartic  16 psi_n^2 ( 1 - psi_n )^2, singly peaked and symmetric; j = 1.
  skew     C psi_n^3 ( 1 - psi_n )^4, peak at 3/7; j = 2 inner, 3 outer. MR-4's
           shape (export_rotating.py).
  double   C psi_n^3 ( 1 - psi_n )^3 [ ( psi_n - 0.55 )^2 + 0.01 ]: two interior
           peaks, 1 at psi_n = 0.280 and 0.485 at 0.750, with a dip to 0.222 at
           0.557; j = 2 at both. A mirror's profiles are not single humps.

j is the EDGE ORDER, p' ~ d^j at distance d inside a label, and it caps MEQ's
convergence order (CLAUDE_FB.md, "The plasma edge caps the order").

USAGE, from source/:

    python3 tools/geq-benchmark/shapes.py            # write the four cases
    python3 tools/geq-benchmark/shapes.py --check    # assert they are current
"""

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(HERE, "cases")

EV = 1.602176634e-19
ELECTRON_MASS = 9.1093837139e-31        # scipy's CODATA 2022, as geq takes it
DEUTERON_MASS = 3.3435837768e-27
MU0 = 4.0e-7 * math.pi                  # freegs4e.gradshafranov's, in J_phi

TE_PEAK, TE_EDGE = 1.28e3, 1.28e2       # eV: geq's tools/run.py
TD_PEAK, TD_EDGE = 1.53e3, 1.53e2       # eV: run.py's td_scale = 1.53/1.28
NE_PEAK, NE_EDGE = 2.0e19, 2.0e13       # m^-3 at z = 0: mach_equilibria.py's peak
OMEGA_EDGE_FRACTION = 0.1               # geq cannot start from omega = 0 at a label
#: geq's Mach omega R/sqrt( 2 T_D/m_D ) = 2 at R = 0.25 m and the peak T_D.
OMEGA_PEAK = 2.0 * math.sqrt(2.0 * TD_PEAK * EV / DEUTERON_MASS) / 0.25
KNOTS = 1025

SKEW_C = 1.0 / ((3.0 / 7.0) ** 3 * (4.0 / 7.0) ** 4)
DOUBLE_DIP, DOUBLE_WIDTH2 = 0.55, 0.01


def _double_raw(x):
    return (x * (1.0 - x)) ** 3 * ((x - DOUBLE_DIP) ** 2 + DOUBLE_WIDTH2)


def _double_peak():
    """The higher peak, by bisection on g' = x^2 ( 1 - x )^2 r( x ) inside the
    first hump: r( 0.1 ) > 0 > r( 0.4 )."""
    a, b = DOUBLE_DIP, DOUBLE_WIDTH2

    def r(x):
        return 3.0 * (1.0 - 2.0 * x) * ((x - a) ** 2 + b) + 2.0 * x * (1.0 - x) * (x - a)
    lo, hi = 0.1, 0.4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if r(mid) > 0.0 else (lo, mid)
    return 0.5 * (lo + hi)


DOUBLE_C = 1.0 / _double_raw(_double_peak())

#: name: ( one-line description, the shape as Python source in x, its slope ).
#: The source strings ARE the definition: they are evaluated here for MEQ's
#: tables and written verbatim into the geq script.
SHAPES = {
    "flattop": ("( 1 - u^4 )^2, u = 2 psi_n - 1: flat-topped, symmetric, j = 1",
                "(1.0 - (2.0 * x - 1.0) ** 4) ** 2",
                "-16.0 * (2.0 * x - 1.0) ** 3 * (1.0 - (2.0 * x - 1.0) ** 4)"),
    "quartic": ("16 psi_n^2 ( 1 - psi_n )^2: singly peaked, symmetric, j = 1",
                "16.0 * (x * (1.0 - x)) ** 2",
                "32.0 * x * (1.0 - x) * (1.0 - 2.0 * x)"),
    "skew": ("C psi_n^3 ( 1 - psi_n )^4, peak at 3/7: j = 2 inner, 3 outer",
             f"{SKEW_C!r} * x ** 3 * (1.0 - x) ** 4",
             f"{SKEW_C!r} * x ** 2 * (1.0 - x) ** 3 * (3.0 - 7.0 * x)"),
    "double": ("C psi_n^3 ( 1 - psi_n )^3 [ ( psi_n - 0.55 )^2 + 0.01 ]: two"
               " interior peaks, j = 2",
               f"{DOUBLE_C!r} * (x * (1.0 - x)) ** 3"
               f" * ((x - {DOUBLE_DIP!r}) ** 2 + {DOUBLE_WIDTH2!r})",
               f"{DOUBLE_C!r} * (x * (1.0 - x)) ** 2 * (3.0 * (1.0 - 2.0 * x)"
               f" * ((x - {DOUBLE_DIP!r}) ** 2 + {DOUBLE_WIDTH2!r})"
               f" + 2.0 * x * (1.0 - x) * (x - {DOUBLE_DIP!r}))"),
}

#: ( MEQ table stem, peak, edge, header ). Temperatures in eV, read through
#: TemperatureScale; the density is where each surface crosses z = 0.
PROFILES = (
    ("omega", OMEGA_PEAK, OMEGA_EDGE_FRACTION * OMEGA_PEAK, "omega( Psi ) in rad/s"),
    ("T_electron", TE_PEAK, TE_EDGE, "T_e( Psi ) in eV"),
    ("T_deuterium", TD_PEAK, TD_EDGE, "T_D( Psi ) in eV"),
    ("n_electron", NE_PEAK, NE_EDGE,
     "n_e( Psi ) in m^-3 WHERE EACH FLUX SURFACE CROSSES z = 0"),
)


def evaluate(source, x):
    return eval(source, {"__builtins__": {}}, {"x": x})


def table(name, peak, edge):
    """MEQ's three columns, Psi f df/dPsi, ascending in Psi = 1 - psi_n."""
    _, g, dg = SHAPES[name]
    rows = []
    for i in range(KNOTS):
        psi = i / (KNOTS - 1.0)
        x = 1.0 - psi
        rows.append((psi, edge + (peak - edge) * evaluate(g, x),
                     -(peak - edge) * evaluate(dg, x)))
    return rows


def table_text(name, stem, peak, edge, what):
    lines = [f"# {what}. Case {name}: f_edge + ( f_peak - f_edge ) g( psi_n ),",
             f"# g = {SHAPES[name][0]}.",
             f"# f_peak {peak:.17e}, f_edge {edge:.17e}. Psi = 1 - psi_n.",
             "# WRITTEN BY tools/geq-benchmark/shapes.py. Columns: Psi f df/dPsi."]
    lines += [f"{p:.17e} {f:.17e} {d:.17e}" for p, f, d in table(name, peak, edge)]
    return "\n".join(lines) + "\n"


GEQ_SCRIPT = '''#!/usr/bin/env python
"""geq ON THE MIRROR, PROFILE SHAPE "{name}": {description}.

WRITTEN BY tools/geq-benchmark/shapes.py, which writes meq.toml beside it from
the same profiles. Edit shapes.py and re-run it; do not edit this file.

A geq user's script, in the shape of geq's own tools/run.py: the machine, the
profiles and the solve are stated in code, and only geq's and freegs4e's own
interface is called. Run it from a run directory,

    PYTHONPATH=<geq>/src python <this file>

and it writes OUTPUT there: psi on geq's grid, the labels, and the solve time.

THE PROBLEM: tools/geq-benchmark/shapes.py's docstring. The electron density is
given at z = 0 on each flux surface and the deuterons follow by quasineutrality;
the density falls to 1e-6 of its peak at both labels, the inner read at the
electrode ( 0.10, 0 ) and the outer at the limiter ( 0.35, 0 ), every Picard
step.
The boundary is freegs4e's default, the direct Green's-function integral.
"""
import time

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from freegs4e.boundary import freeBoundary
from freegs4e.equilibrium import Equilibrium
from freegs4e.jtor import ProfilesCentrifugalMirror
from freegs4e.species import AllSpecies, SingleSpecies

from geq.machine import make_mirror_machine_from_field
from geq.solver import solve_equilibrium

# THE RUN: the grid, the Picard tolerances and the output. race_geq.py changes
# these lines and nothing else. ATOL = 0 leaves RTOL alone to decide.
N = 257
RTOL = 1.0e-10
ATOL = 0.0
MAXITS = 400
OUTPUT = "geq-{name}.npz"

# The machine: two packs, centres 3 m apart, as tools/run.py builds them.
FIELD = 4.0
COIL_INNER_RADIUS = 0.7855
COIL_LENGTH = 1.0
COIL_SEPARATION = 3.0
R_INNER, R_OUTER = 0.10, 0.35          # the electrode and the limiter, at z = 0

# The profiles, f_edge + ( f_peak - f_edge ) g( psi_n ).
ELECTRON_MASS = {electron_mass!r}
DEUTERON_MASS = {deuteron_mass!r}
TE_PEAK, TE_EDGE = {te_peak!r}, {te_edge!r}          # eV
TD_PEAK, TD_EDGE = {td_peak!r}, {td_edge!r}          # eV
NE_PEAK, NE_EDGE = {ne_peak!r}, {ne_edge!r}          # m^-3 at z = 0
OMEGA_PEAK, OMEGA_EDGE = {omega_peak!r}, {omega_edge!r}   # rad/s
KNOTS = {knots}


def shape(x):
    """{description}"""
    return {shape}


def main():
    psi_n = np.linspace(0.0, 1.0, KNOTS)
    g = shape(psi_n)
    te = TE_EDGE + (TE_PEAK - TE_EDGE) * g
    td = TD_EDGE + (TD_PEAK - TD_EDGE) * g
    ne = NE_EDGE + (NE_PEAK - NE_EDGE) * g
    omega = OMEGA_EDGE + (OMEGA_PEAK - OMEGA_EDGE) * g

    zwall = [-1.5 * COIL_SEPARATION, -1.5 * COIL_SEPARATION,
             1.5 * COIL_SEPARATION, 1.5 * COIL_SEPARATION]
    rwall = [0.0, 0.5, 0.5, 0.0]
    machine, _ = make_mirror_machine_from_field(
        field=FIELD, L=COIL_LENGTH, ri=COIL_INNER_RADIUS, dz=COIL_SEPARATION,
        Nr=1, Nz=10, rwall=rwall, zwall=zwall)

    plasma_psi_guess = np.zeros((N, N))
    eq = Equilibrium(tokamak=machine, Rmin=0.05, Rmax=0.63,
                     Zmin=-0.5 * COIL_SEPARATION, Zmax=0.5 * COIL_SEPARATION,
                     nx=N, ny=N, boundary=freeBoundary, psi=plasma_psi_guess,
                     order=4)
    eq._find_critical = False
    machine_psi = eq.tokamak.calcPsiFromGreens(eq._pgreen)
    eq._updatePlasmaPsi(plasma_psi_guess)

    electron = SingleSpecies(name="electron", z=-1.0, mass=ELECTRON_MASS,
                             psi_n=psi_n, temp=te, density_Z0=ne,
                             flux_density=ne)
    deuterium = SingleSpecies(name="deuterium", z=1.0, mass=DEUTERON_MASS,
                              psi_n=psi_n, temp=td)
    all_species = AllSpecies(species_list=[electron, deuterium])
    all_species.zeff_dict = {{}}

    omega_profile = InterpolatedUnivariateSpline(psi_n, omega, k=3, ext=3)
    psi_axis = eq.psiRZ(R_INNER, 0.0)
    psi_bndry = eq.psiRZ(R_OUTER, 0.0)
    R_guess = np.linspace(R_INNER, R_OUTER, 101)
    psi_n_guess = (machine.psi(R_guess, 0.0) - psi_axis) / (psi_bndry - psi_axis)
    all_species.set_omega(psi_n=psi_n_guess, R=R_guess,
                          psi_range=psi_bndry - psi_axis,
                          omega_profile=omega_profile)
    bndry = dict(r_inner=R_INNER, z_inner=0.0, r_outer=R_OUTER, z_outer=0.0,
                 mach=None, voltage=None, omega_profile=omega_profile)

    profiles = ProfilesCentrifugalMirror(all_species=all_species)
    started = time.perf_counter()
    solve_equilibrium(eq, profiles, machine_psi, solve_method="picard",
                      maxits=MAXITS, rtol=RTOL, atol=ATOL, psi_axis=psi_axis,
                      psi_bndry=psi_bndry, bndry=bndry)
    seconds = time.perf_counter() - started

    np.savez(OUTPUT, psi=np.asarray(eq.psi()), R=np.asarray(eq.R),
             Z=np.asarray(eq.Z), psi_axis=eq.psi_axis, psi_bndry=eq.psi_bndry,
             solve_seconds=seconds)
    print(f"geq: {{N}} x {{N}}, solved in {{seconds:.2f}} s; psi_axis"
          f" {{eq.psi_axis:.12e}}, psi_bndry {{eq.psi_bndry:.12e}}; wrote {{OUTPUT}}")


if __name__ == "__main__":
    main()
'''


MEQ_CONFIG = '''# MEQ ON THE MIRROR, PROFILE SHAPE "{name}": {description}.
#
# WRITTEN BY tools/geq-benchmark/shapes.py, which writes geq_{name}.py beside it
# from the same profiles. Edit shapes.py and re-run it; do not edit this file.
# From a run directory:
#
#     meq-run <this file>
#
# THE PROBLEM: tools/geq-benchmark/shapes.py's docstring, the one geq_{name}.py
# states in code. The electron density is given where each flux surface
# crosses z = 0 -- Gauge = "midplane" -- and the deuterons are Neutralising, as
# geq derives them. psi_inner is read at the electrode, psi_bnd at the limiter,
# and the density falls to 1e-6 of its peak at both, which is what leaves the
# source zero there under this gauge: a flat density alone would not, where
# omega does not vanish.

[mesh]
File = "mirror-{name}.msh"

[mesh.generate]
Tool = "halfdisc"
Radius = 3.4
Size = 0.4
PlasmaRMin = 0.0
PlasmaRMax = 0.75
PlasmaZMin = -1.5
PlasmaZMax = 1.5
PlasmaSize = 0.1
AlignPlasma = true

[discretisation]
PolynomialDegree = 2

[solver]
NewtonRelativeTolerance = 1.0e-12
NewtonAbsoluteTolerance = 1.0e-14
NewtonMaxIterations = 30

[source]
Type = "rotating"
Normalised = true
ConfineToPlasma = true
SupportShape = "box"
AxialEdge = 1.5
Gauge = "midplane"
Mu0 = {mu0!r}
GGPrime = 0.0
OmegaFile = "omega.dat"

[[source.species]]
Name = "electron"
Mass = {electron_mass!r}
Charge = -1.0
TemperatureFile = "T_electron.dat"
TemperatureScale = {ev!r}
DensityFile = "n_electron.dat"

[[source.species]]
Name = "deuterium"
Mass = {deuteron_mass!r}
Charge = 1.0
TemperatureFile = "T_deuterium.dat"
TemperatureScale = {ev!r}
Neutralising = true

[boundary]
Type = "zero"

[boundary.electrode]
R = 0.1
Z = 0.0

[boundary.limiter]
R = 0.35
Z = 0.0

[boundary.exterior]
Radius = 3.0
Modes = 32

# geq's packs as geq carries them, ten filaments each at the radial centre,
# subtracted rather than meshed: examples/mirror-rotating.toml's machine.
[conductors]
Model = "filament"

[[coils]]
Name = "lower"
CentreR = 0.82078441790868473
CentreZ = -1.5
HalfWidth = 0.035284417908684695
HalfHeight = 0.5
Current = 5927782.2086590333
FilamentsR = 1
FilamentsZ = 10

[[coils]]
Name = "upper"
CentreR = 0.82078441790868473
CentreZ = 1.5
HalfWidth = 0.035284417908684695
HalfHeight = 0.5
Current = 5927782.2086590333
FilamentsR = 1
FilamentsZ = 10

[output]
Directory = "."
Prefix = "mirror-{name}"
GridNR = 129
GridNZ = 129
'''


def files(name):
    """{ file name: text } for one case."""
    description, g, _ = SHAPES[name]
    values = dict(name=name, description=description, shape=g, knots=KNOTS,
                  electron_mass=ELECTRON_MASS, deuteron_mass=DEUTERON_MASS,
                  te_peak=TE_PEAK, te_edge=TE_EDGE, td_peak=TD_PEAK,
                  td_edge=TD_EDGE, ne_peak=NE_PEAK, ne_edge=NE_EDGE,
                  omega_peak=OMEGA_PEAK,
                  omega_edge=OMEGA_EDGE_FRACTION * OMEGA_PEAK, mu0=MU0, ev=EV)
    out = {f"geq_{name}.py": GEQ_SCRIPT.format(**values),
           "meq.toml": MEQ_CONFIG.format(**values)}
    for stem, peak, edge, what in PROFILES:
        out[f"{stem}.dat"] = table_text(name, stem, peak, edge, what)
    return out


def check(names):
    """The names of every tracked case file that differs from what this script
    writes, or is missing."""
    stale = []
    for name in names:
        for fname, text in files(name).items():
            path = os.path.join(CASES, name, fname)
            try:
                with open(path) as fh:
                    if fh.read() != text:
                        stale.append(path)
            except FileNotFoundError:
                stale.append(path)
    return stale


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any tracked case is not what this writes")
    ap.add_argument("shapes", nargs="*", default=list(SHAPES))
    args = ap.parse_args(argv)
    for name in args.shapes:
        if name not in SHAPES:
            raise SystemExit(f"unknown shape {name}; one of {', '.join(SHAPES)}")
    if args.check:
        stale = check(args.shapes)
        for path in stale:
            print(f"STALE: {os.path.relpath(path)}")
        return 1 if stale else 0
    for name in args.shapes:
        os.makedirs(os.path.join(CASES, name), exist_ok=True)
        for fname, text in files(name).items():
            with open(os.path.join(CASES, name, fname), "w") as fh:
                fh.write(text)
        # The two edges, which the source's vanishing there rests on.
        _, g, dg = SHAPES[name]
        print(f"{name:8s} g( 0 ) {evaluate(g, 0.0):.1e} g( 1 ) {evaluate(g, 1.0):.1e}"
              f"  g'( 0 ) {evaluate(dg, 0.0):.1e} g'( 1 ) {evaluate(dg, 1.0):.1e}"
              f"  -> {os.path.relpath(os.path.join(CASES, name))}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
