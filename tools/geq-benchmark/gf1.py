#!/usr/bin/env python3
"""THE GF1 CASE: geq's GF1 working point with its profiles tapered at the
labels, written as a geq user's script and a MEQ user's configuration.

cases/gf1/ holds what shapes.py's cases hold -- geq_gf1.py, meq.toml and MEQ's
tables -- and race_geq.py runs it the same way. This script is the ONE place
the case is defined, and `--check` asserts the tracked files are what it
writes.

THE MACHINE is shapes.py's, which is GF1's: geq's two-pack at 4 T, inner radius
0.7855 m, length 1 m, centres 3 m apart, ten filaments each, the labels at the
electrode ( 0.10, 0 ) and the limiter ( 0.35, 0 ), and g = 0.

THE PLASMA is GF1's 0.4 MV point at n_e = 1e20 (MEASUREMENTS.md M-341), held
and not power balanced: electrons, deuterons and fully stripped carbon at
Z_eff = 1.55, carbon at 6 m_D as geq's power_balance.generate_all_species
builds it.

  temperatures  geq's GF1 profiles, geq.profiles.generate_temperature_profile(
                psi_n, peak, 10 eV, alpha = 1 ): a Gaussian in psi_n about 0.5
                of width 0.2. The peaks are M-341's balanced ones, held here
                as inputs: 4911, 7817 and 8020 eV.
  n_e( z = 0 )  geq's GF1 profile, generate_density_z0_profile( psi_n, 1e20,
                1e14, alpha = 1 ), TAPERED: n_e = 1e14 + ( G - 1e14 ) t( psi_n ).
  n_C( z = 0 )  ( Z_eff - 1 )/( Z ( Z - 1 ) ) n_e = 0.55/30 n_e, as freegs4e's
                AllSpecies.set_flux_density derives it; the deuterons take the
                rest by quasineutrality.
  omega         M-341's solved omega( psi_n ), held: geq's voltage drive moves
                omega with the equilibrium every Picard step, which MEQ does
                not pose. A degree-8 polynomial fitted to it, within 1.0e-3 of
                its peak; the coefficients below ARE the profile.

THE TAPER. geq's profiles do not vanish at the labels -- the densities sit at
4.4% of their peak there and geq zeroes them outside -- and MEQ's box needs the
source to vanish at both (docs/rotation.rst, the box's precondition). So the
density, and only the density, is multiplied by

    t( psi_n ) = S( psi_n/w ) S( ( 1 - psi_n )/w ),  w = 0.1,
    S( u ) = u^3 ( 10 - 15 u + 6 u^2 ) for u < 1, and 1 beyond,

which is 1 on 0.1 <= psi_n <= 0.9, joins it with two continuous derivatives,
and vanishes at each label to third order: the edge order j = 3. There the
densities fall to geq's own floor, 1e14, a millionth of the peak. The
temperatures and omega are not tapered: under the midplane statement the source
vanishes with the density, and geq's first phi_0 guess needs omega nonzero
at the labels (shapes.py).

THE STABILISATION IS THE CASE'S OWN: tau = 30/max( R, 0.1 m ), not the box
default tau0 = 3. GF1 is ten times M-294's pressure, and at tau0 = 3 the coarse
rungs fold inside one element -- the reaction (1/R) dF/dpsi reaches the local
block's lowest eigenvalue (k + 1) tau |dK|/|K| (MEASUREMENTS.md M-345) -- so
they do not solve from the vacuum. Of 3, 10 and 30, 30 is the smallest that
solves every race rung from the vacuum, and it moves the L2 k = 4 reference by
2e-7 of a max |psi| of 0.88 (M-347). This is a statement about the benchmark,
not a default: a general tau that follows the reaction is TODO's.

USAGE, from source/:

    python3 tools/geq-benchmark/gf1.py            # write cases/gf1/
    python3 tools/geq-benchmark/gf1.py --check    # assert it is current
"""

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import shapes                                                   # noqa: E402

NAME = "gf1"
DIRECTORY = os.path.join(shapes.CASES, NAME)

CARBON_MASS = 6.0 * shapes.DEUTERON_MASS    # generate_all_species' carbon_mass
CARBON_CHARGE = 6.0
ZEFF_CARBON = 0.55                          # gf1_scan_points.ZEFF: Z_eff 1.55
CARBON_FRACTION = ZEFF_CARBON / (CARBON_CHARGE ** 2 - CARBON_CHARGE)

TE_PEAK, TD_PEAK, TC_PEAK = 4911.0, 7817.0, 8020.0     # eV, M-341
T_EDGE = 10.0                               # eV, generate_all_species' t_edge
NE_PEAK, NE_EDGE = 1.0e20, 1.0e14           # m^-3 at z = 0: ne_axis, ne_edge
CENTRE, WIDTH = 0.5, 0.2                    # generate_gaussian_profile's
TAPER_WIDTH = 0.1
#: omega( psi_n ) = sum c_i psi_n^i in rad/s: M-341's, geq c86ced9.
OMEGA = (6781989.6621, -2125179.7763, 166873913.29, -781369920.96,
         2222605962.6, -4221170741.0, 4674349401.8, -2681912084.4,
         617912791.57)
KNOTS = 1025


def gaussian(x):
    """generate_gaussian_profile's shape at alpha = 1, and its slope."""
    e = math.exp(-((x - CENTRE) ** 2) / (2.0 * WIDTH ** 2))
    return e, -e * (x - CENTRE) / WIDTH ** 2


def step(u):
    """S( u ) and S'( u )."""
    if u >= 1.0:
        return 1.0, 0.0
    return u ** 3 * (10.0 - 15.0 * u + 6.0 * u ** 2), 30.0 * u ** 2 * (1.0 - u) ** 2


def taper(x):
    a, da = step(x / TAPER_WIDTH)
    b, db = step((1.0 - x) / TAPER_WIDTH)
    return a * b, (da * b - a * db) / TAPER_WIDTH


def temperature(peak, x):
    g, dg = gaussian(x)
    return T_EDGE + (peak - T_EDGE) * g, (peak - T_EDGE) * dg


def electron_density(x):
    g, dg = gaussian(x)
    t, dt = taper(x)
    amplitude = NE_PEAK - NE_EDGE
    return NE_EDGE + amplitude * g * t, amplitude * (dg * t + g * dt)


def omega(x):
    value = slope = 0.0
    for c in reversed(OMEGA):
        slope = slope * x + value
        value = value * x + c
    return value, slope


def carbon_density(x):
    n, dn = electron_density(x)
    return CARBON_FRACTION * n, CARBON_FRACTION * dn


#: ( MEQ table stem, f( psi_n ) -> ( value, d/dpsi_n ), header ).
PROFILES = (
    ("omega", omega, "omega( Psi ) in rad/s, M-341's, held"),
    ("T_electron", lambda x: temperature(TE_PEAK, x), "T_e( Psi ) in eV"),
    ("T_deuterium", lambda x: temperature(TD_PEAK, x), "T_D( Psi ) in eV"),
    ("T_carbon", lambda x: temperature(TC_PEAK, x), "T_C( Psi ) in eV"),
    ("n_electron", electron_density,
     "n_e( Psi ) in m^-3 WHERE EACH FLUX SURFACE CROSSES z = 0, tapered"),
    ("n_carbon", carbon_density,
     "n_C( Psi ) = 0.55/30 n_e( Psi ) in m^-3 at z = 0, tapered"),
)


def table_text(profile, what):
    lines = [f"# {what}. Case gf1: tools/geq-benchmark/gf1.py's docstring.",
             "# Psi = 1 - psi_n.",
             "# WRITTEN BY tools/geq-benchmark/gf1.py. Columns: Psi f df/dPsi."]
    for i in range(KNOTS):
        psi = i / (KNOTS - 1.0)
        value, slope = profile(1.0 - psi)
        lines.append(f"{psi:.17e} {value:.17e} {-slope:.17e}")
    return "\n".join(lines) + "\n"


GEQ_SCRIPT = '''#!/usr/bin/env python
"""geq ON GF1 AT n_e = 1e20, ITS PROFILES TAPERED AT THE LABELS.

WRITTEN BY tools/geq-benchmark/gf1.py, which writes meq.toml beside it from the
same profiles. Edit gf1.py and re-run it; do not edit this file.

A geq user's script, in the shape of geq's own tools/run.py: the machine, the
profiles and the solve are stated in code, and only geq's and freegs4e's own
interface is called. No power balance: the temperatures, the densities and
omega are held. Run it from a run directory,

    PYTHONPATH=<geq>/src python <this file>

and it writes OUTPUT there: psi on geq's grid, the labels, and the solve time.

THE PROBLEM: tools/geq-benchmark/gf1.py's docstring. The electron density is
given at z = 0 on each flux surface, carbon at 0.55/30 of it through Z_eff,
and the deuterons follow by quasineutrality.
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
from geq.profiles import generate_density_z0_profile, generate_temperature_profile
from geq.solver import solve_equilibrium

# THE RUN: the grid, the Picard tolerances and the output. race_geq.py changes
# these lines and nothing else. ATOL = 0 leaves RTOL alone to decide.
N = 257
RTOL = 1.0e-10
ATOL = 0.0
MAXITS = 400
OUTPUT = "geq-gf1.npz"

# The machine: two packs, centres 3 m apart, as tools/run.py builds them.
FIELD = 4.0
COIL_INNER_RADIUS = 0.7855
COIL_LENGTH = 1.0
COIL_SEPARATION = 3.0
R_INNER, R_OUTER = 0.10, 0.35          # the electrode and the limiter, at z = 0

# The plasma: GF1 at n_e = 1e20, held (MEASUREMENTS.md M-341).
ELECTRON_MASS = {electron_mass!r}
DEUTERON_MASS = {deuteron_mass!r}
CARBON_MASS = {carbon_mass!r}
ZEFF = {{"carbon": {zeff_carbon!r}}}
TE_PEAK, TD_PEAK, TC_PEAK, T_EDGE = {te_peak!r}, {td_peak!r}, {tc_peak!r}, {t_edge!r}
NE_PEAK, NE_EDGE = {ne_peak!r}, {ne_edge!r}          # m^-3 at z = 0
TAPER_WIDTH = {taper_width!r}
OMEGA = {omega!r}                                    # rad/s, in powers of psi_n
KNOTS = {knots}


def taper(x):
    """gf1.py's t( psi_n ): 1 on [ w, 1 - w ], zero to third order at each label."""
    def step(u):
        u = np.minimum(u, 1.0)
        return u ** 3 * (10.0 - 15.0 * u + 6.0 * u ** 2)
    return step(x / TAPER_WIDTH) * step((1.0 - x) / TAPER_WIDTH)


def main():
    psi_n = np.linspace(0.0, 1.0, KNOTS)
    te = generate_temperature_profile(psi_n, TE_PEAK, T_EDGE, 1.0)
    td = generate_temperature_profile(psi_n, TD_PEAK, T_EDGE, 1.0)
    tc = generate_temperature_profile(psi_n, TC_PEAK, T_EDGE, 1.0)
    gaussian = generate_density_z0_profile(psi_n, NE_PEAK, NE_EDGE, 1.0)
    ne = NE_EDGE + (gaussian - NE_EDGE) * taper(psi_n)
    omega = np.polynomial.polynomial.polyval(psi_n, OMEGA)

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
    carbon = SingleSpecies(name="carbon", z=6.0, mass=CARBON_MASS,
                           psi_n=psi_n, temp=tc)
    all_species = AllSpecies(species_list=[electron, deuterium, carbon])
    all_species.zeff_dict = dict(ZEFF)

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


#: GF1's tau0, in m^-1, under TauScaling = "diffusion" (the docstring).
TAU0 = 30.0


def meq_config():
    """shapes.py's configuration for this case: the same machine, labels, mesh
    and output, with carbon added and the profile comments replaced."""
    text = shapes.MEQ_CONFIG.format(
        name=NAME, description="GF1 at n_e = 1e20, tapered (gf1.py)",
        mu0=shapes.MU0, electron_mass=shapes.ELECTRON_MASS,
        deuteron_mass=shapes.DEUTERON_MASS, ev=shapes.EV)
    for old, new in (
            ("# WRITTEN BY tools/geq-benchmark/shapes.py, which writes geq_gf1.py beside it\n"
             "# from the same profiles. Edit shapes.py and re-run it; do not edit this file.",
             "# WRITTEN BY tools/geq-benchmark/gf1.py, which writes geq_gf1.py beside it\n"
             "# from the same profiles. Edit gf1.py and re-run it; do not edit this file."),
            ("# THE PROBLEM: tools/geq-benchmark/shapes.py's docstring, the one geq_gf1.py",
             "# THE PROBLEM: tools/geq-benchmark/gf1.py's docstring, the one geq_gf1.py"),
            ("# geq derives them. psi_inner is read at the electrode, psi_bnd at the limiter,\n"
             "# and the density falls to 1e-6 of its peak at both,",
             "# geq derives them; carbon is 0.55/30 of the electrons, Z_eff = 1.55.\n"
             "# psi_inner is read at the electrode, psi_bnd at the limiter, and the\n"
             "# taper takes the densities to 1e-6 of their peak at both,"),
            ('TemperatureFile = "T_deuterium.dat"\n'
             f"TemperatureScale = {shapes.EV!r}\n"
             "Neutralising = true\n",
             'TemperatureFile = "T_deuterium.dat"\n'
             f"TemperatureScale = {shapes.EV!r}\n"
             "Neutralising = true\n\n"
             "[[source.species]]\n"
             'Name = "carbon"\n'
             f"Mass = {CARBON_MASS!r}\n"
             f"Charge = {CARBON_CHARGE!r}\n"
             'TemperatureFile = "T_carbon.dat"\n'
             f"TemperatureScale = {shapes.EV!r}\n"
             'DensityFile = "n_carbon.dat"\n'),
            ("[discretisation]\nPolynomialDegree = 2\n",
             "[discretisation]\nPolynomialDegree = 2\n"
             "# GF1's own stabilisation, tau = 30/max( R, 0.1 m ): at the box default\n"
             "# tau0 = 3 the coarse rungs fold inside one element and do not solve from\n"
             "# the vacuum (gf1.py's docstring; MEASUREMENTS.md M-345, M-347).\n"
             'TauScaling = "diffusion"\n'
             f"Tau = {TAU0!r}\n")):
        if text.count(old) != 1:
            raise SystemExit(f"gf1.py: shapes.MEQ_CONFIG no longer holds '{old[:50]}'")
        text = text.replace(old, new)
    return text


def files():
    """{ file name: text } for the case."""
    values = dict(electron_mass=shapes.ELECTRON_MASS,
                  deuteron_mass=shapes.DEUTERON_MASS, carbon_mass=CARBON_MASS,
                  zeff_carbon=ZEFF_CARBON, te_peak=TE_PEAK, td_peak=TD_PEAK,
                  tc_peak=TC_PEAK, t_edge=T_EDGE, ne_peak=NE_PEAK,
                  ne_edge=NE_EDGE, taper_width=TAPER_WIDTH, omega=OMEGA,
                  knots=KNOTS)
    out = {f"geq_{NAME}.py": GEQ_SCRIPT.format(**values), "meq.toml": meq_config()}
    for stem, profile, what in PROFILES:
        out[f"{stem}.dat"] = table_text(profile, what)
    return out


def check():
    """The path of every tracked case file that differs from what this writes."""
    stale = []
    for fname, text in files().items():
        path = os.path.join(DIRECTORY, fname)
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
                    help="exit 1 if the tracked case is not what this writes")
    args = ap.parse_args(argv)
    if args.check:
        stale = check()
        for path in stale:
            print(f"STALE: {os.path.relpath(path)}")
        return 1 if stale else 0
    os.makedirs(DIRECTORY, exist_ok=True)
    for fname, text in files().items():
        with open(os.path.join(DIRECTORY, fname), "w") as fh:
            fh.write(text)
    n0, dn0 = electron_density(0.0)
    n1, dn1 = electron_density(1.0)
    print(f"gf1  n_e( 0 ) {n0:.3e} n_e( 1 ) {n1:.3e}  n_e'( 0 ) {dn0:.1e}"
          f" n_e'( 1 ) {dn1:.1e}  omega( 0 ) {omega(0.0)[0]:.4e}"
          f" omega( 1 ) {omega(1.0)[0]:.4e}  -> {os.path.relpath(DIRECTORY)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
