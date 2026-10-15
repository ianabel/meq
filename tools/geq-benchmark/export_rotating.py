#!/usr/bin/env python
"""EXPORT geq'S FREE-BOUNDARY ROTATING MIRROR, FOR MIRROR-PLAN.md MR-4.

A SOLVED reference, unlike export_mirror.py's prescribed state: geq's Picard
iteration on its two-pack mirror machine (coils fixed, 2 x 10 filaments, labels
read at ( 0.10, 0 ) and ( 0.35, 0 ) every step), with electrons and deuterons,
rigid rotation omega( psi_n ) and the quasineutrality potential phi_0. MEQ's
MR-4b and MR-4d ctests read what this writes.

IT IS NOT SHIPPED geq, FOR THREE REASONS MIRROR-PLAN.md, MR-4 MEASURES:

  1. SHIPPED geq IS BIASED BY ABOUT 1% OF THE PLASMA FLUX, and its own
     self-convergence cannot see it: freegs4e/boundary.py's von Hagenow offset
     is a fixed eps = 1e-2 m. So the reference uses the DIRECT `freeBoundary`
     Green's-function integral, and the shipped route is solved once, at 257^2,
     as the CONTROL that shows a comparison can resolve that bias.
  2. geq'S SHIPPED PROFILES STEP AT BOTH LABELS -- Gaussians peaked at
     psi_n = 0.5 -- and MEQ's support cannot take a source that is non-zero
     where the support is read off the solution. So every profile here is
         f( psi_n ) = f_edge + ( f_peak - f_edge ) C psi_n^3 ( 1 - psi_n )^4,
     flat to order 2 at the inner label and 3 at the outer, and deliberately
     asymmetric (see SHAPES below); it is evaluated ANALYTICALLY inside
     geq (below), not through geq's cubic splines, whose not-a-knot end
     conditions would leave a slope of order h^3 at the labels.
  3. geq'S phi_0 RESIDUAL IS NOT ITERATION ERROR, so phi_0 is TIGHTENED -- a
     pointwise Newton on Sum_s Z_s n_s = 0 after geq's global solve -- inside
     EVERY Picard step, so that the fixed point is a fixed point of the
     consistent problem.

AND geq NEEDS ONE MORE HOLD. It rebuilds the deuteron's flux-function density
N_D every step from the electrons' MIDPLANE density (AllSpecies.set_flux_density,
called from both quasineutrality solves). Here both N_e and N_D are PRESCRIBED
flux functions and phi_0 alone carries quasineutrality, so set_flux_density is
made a no-op and the deuteron is flagged `_density_eval`, which makes geq's
global solve read its prescribed N_D rather than re-deriving one. With both
N_s prescribed there is no gauge freedom left: phi_0 is fixed pointwise by
quasineutrality, and geq's one-point pin ( phi_0 at psi_n = 0.5 on the
midplane ) becomes a wrong value at one node -- which the tightening corrects.

geq and freegs4e are NOT modified; everything is a monkeypatch from here.

THE EXPORT IS IN MEQ'S CONVENTIONS: Psi = 1 - psi_n (abscissa reversed,
derivative columns negated), T_s in JOULES, omega in rad/s, masses in kg,
charge numbers, and densities n_s0( Psi ) ON THE CURVE R = R_ref -- MEQ's
fixed-R gauge phi_0( R_ref, psi ) = 0 -- transferred pointwise exactly as
export_mirror.py does. geq's own-gauge N_s are written too, as Nraw_*.dat.

USAGE, from source/

    PYTHONPATH=/home/ian/projects/geq/src \\
    tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_rotating.py

About eleven minutes on this machine, most of it the 513^2 direct solve.
`--cache DIR` keeps each solve as an .npz and reuses it on a re-run.

AND THE MIDPLANE GAUGE'S TABLES, nmid_<name>.dat: the same densities where each
surface crosses z = 0, MEQ's Gauge::Midplane (MIRROR-PLAN.md section 5). geq's
N_s are its zero-potential densities, so the transfer needs the crossings of
geq's solved state -- the finest direct one -- and these tables are an OUTPUT,
written beside geq's answer rather than with the tracked problem.

geq's ANSWER, psi_direct_<n>.dat, psi_hagenow_<n>.dat and solution.txt, goes to
--out, default runs/geq-benchmark/reference/rotating/, the path the C++ tests
read. THE PROBLEM, meta.txt, filaments.dat and the tables, goes to --tables,
default the same directory, and is reported against the tracked
tools/geq-benchmark/reference/rotating/, which
`--tables tools/geq-benchmark/reference/rotating` refreshes.
"""

import argparse
import contextlib
import io
import os
import subprocess
import sys
import time

import numpy as np
from scipy import constants
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

import matplotlib

matplotlib.use("Agg")

from freegs4e.boundary import freeBoundary, freeBoundaryHagenow
from freegs4e.equilibrium import Equilibrium
from freegs4e.gradshafranov import Greens
from freegs4e.jtor import ProfilesCentrifugalMirror
from freegs4e.species import AllSpecies, SingleSpecies

from geq.machine import make_mirror_machine_from_field
from geq.solver import solve_equilibrium

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from export_mirror import (against_tracked, answers_dir,  # noqa: E402
                           transfer_to_reference_radius)

qe = constants.elementary_charge
MU0 = 4.0e-7 * np.pi
ELECTRON_MASS = constants.electron_mass
DEUTERON_MASS = constants.physical_constants["deuteron mass"][0]

# The machine: geq's tests/conftest.py, export_vacuum.py, VacuumMirror.
FIELD = 4.0
COIL_LENGTH = 1.0
COIL_INNER_RADIUS = 0.7855
COIL_SEPARATION = 3.0
R_MIN = 0.05
R_INNER = 0.10
R_OUTER = 0.35
R_MAX = R_OUTER * 1.8
Z_EXTENT = 0.5 * COIL_SEPARATION

# The profiles. Temperature peaks are conftest's, edges one order down.
#
# N_s IS NOT conftest's 5e19, AND THE REASON IS THE GAUGE. conftest prescribes
# the MIDPLANE density; here N_s is the flux function itself, the density at
# zero potential, and the midplane density is N_s exp( Theta ) with
#     Theta = m_D omega^2 R^2/( 2( T_e + T_D ) )
# ( two species, m_e negligible ), which is about 2.2 at the peak at geq-Mach 2.
# So N_peak = 5e18 puts the PHYSICAL midplane peak near 5e19, conftest's.
#
# omega_edge = omega_peak/5, because Theta grows as R^2 and the edge is cold:
# at the outer label ( R = 0.35, T_e + T_D = 281 eV ) Theta is 39.5 f^2 for
# omega_edge = f omega_peak. With the quartic, f = 0.5 made the outer edge's
# density 200x the peak's and the first Picard step diverged ( measured );
# f = 0.2 converges for every shape here.
TE_PEAK, TE_EDGE = 1.28e3, 1.28e2        # eV
TD_PEAK, TD_EDGE = 1.53e3, 1.53e2        # eV
N_PEAK, N_EDGE = 5.0e18, 5.0e17          # m^-3, N_D = N_e
OMEGA_EDGE_FRACTION = 0.2
R_MACH = {"value": None}                 # set from SHAPES or --r-mach

NAMES = ("electron", "deuterium")
SPECIES = {
    "electron": dict(z=-1.0, mass=ELECTRON_MASS, peak=TE_PEAK, edge=TE_EDGE),
    "deuterium": dict(z=1.0, mass=DEUTERON_MASS, peak=TD_PEAK, edge=TD_EDGE),
}

COARSE = 65                              # the grid every resolution shares


def omega_peak_for(mach):
    """omega at the profiles' peak such that geq's Mach, omega R/sqrt( 2 T_D/m_D ),
    is `mach` there, at the nominal midplane radius R_MACH of the peak."""
    return mach * np.sqrt(2.0 * TD_PEAK * qe / DEUTERON_MASS) / R_MACH["value"]


# ---------------------------------------------------------------------------
# The analytic profiles, in the interface geq's splines present
# ---------------------------------------------------------------------------


# The shape g( x ) of every profile, 1 at its peak x_peak. What matters at a
# label is the EDGE ORDER j, p' ~ d^j at distance d from it: CLAUDE_FB.md, "The
# plasma edge caps the order", has psi* keep k+2 only for k <= j, and the
# best-approximation bound min( k+2, j+2.5 ).
#
#   "skew"     C x^3 ( 1 - x )^4, C = 1/( (3/7)^3 (4/7)^4 ): j = 2 at the inner
#              label, j = 3 at the outer, peak at x = 3/7. THE REFERENCE.
#              Deliberately asymmetric, so that tables read in psi_n instead of
#              Psi give different numbers -- the sense of the abscissa is tested.
#   "quartic"  16 x^2 ( 1 - x )^2, MIRROR-PLAN.md's first choice: j = 1, which
#              would cap MR-4c's rates at k >= 2, and symmetric about 1/2.
#   "sextic"   64 x^3 ( 1 - x )^3: j = 2 at both, symmetric.
SKEW_C = 1.0 / ((3.0 / 7.0) ** 3 * (4.0 / 7.0) ** 4)
SHAPES = {
    # name: ( peak location, default R_MACH -- the midplane radius of the
    #         peak on the converged 513^2 state, rounded to 1 mm )
    "skew": (3.0 / 7.0, 0.240),
    "quartic": (0.5, 0.26),
    "sextic": (0.5, 0.26),
}
SHAPE = {"name": "skew"}


def shape(x):
    if SHAPE["name"] == "skew":
        return SKEW_C * x ** 3 * (1.0 - x) ** 4
    if SHAPE["name"] == "sextic":
        return 64.0 * (x * (1.0 - x)) ** 3
    return 16.0 * (x * (1.0 - x)) ** 2


def shape_slope(x):
    if SHAPE["name"] == "skew":
        return SKEW_C * x ** 2 * (1.0 - x) ** 3 * (3.0 - 7.0 * x)
    if SHAPE["name"] == "sextic":
        return 192.0 * (x * (1.0 - x)) ** 2 * (1.0 - 2.0 * x)
    return 32.0 * x * (1.0 - x) * (1.0 - 2.0 * x)


class Bump:
    """f( x ) = edge + ( peak - edge ) g( x ) on [ 0, 1 ], clamped
    outside as geq's ext=3 splines are. `log=True` presents ln f, which is what
    geq's _flux_density spline holds. Callable, with .derivative(), so it stands
    wherever geq calls an InterpolatedUnivariateSpline."""

    def __init__(self, peak, edge, log=False):
        self.peak, self.edge, self.log = float(peak), float(edge), log

    def value(self, x):
        x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
        return self.edge + (self.peak - self.edge) * shape(x)

    def slope(self, x):
        x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
        return (self.peak - self.edge) * shape_slope(x)

    def __call__(self, x):
        return np.log(self.value(x)) if self.log else self.value(x)

    def derivative(self):
        if self.log:
            return lambda x: self.slope(x) / self.value(x)
        return self.slope


# ---------------------------------------------------------------------------
# The monkeypatches: N_s held, omega analytic, phi_0 tightened
# ---------------------------------------------------------------------------

_ORIG_SET_FLUX_DENSITY = AllSpecies.set_flux_density
_ORIG_SET_OMEGA = AllSpecies.set_omega
_ORIG_QN_GLOBAL = AllSpecies.solve_quasineutrality_global

STATE = {"tighten": True, "rebuilds_refused": 0, "last_qn": None}


def held_set_flux_density(self, psi_n=None, R=None, Z=None, phi0=None):
    """N_s IS PRESCRIBED for every species, so geq must not rebuild one."""
    if all(isinstance(s._flux_density, Bump) for s in self.species_list):
        STATE["rebuilds_refused"] += 1
        return None
    return _ORIG_SET_FLUX_DENSITY(self, psi_n=psi_n, R=R, Z=Z, phi0=phi0)


def analytic_set_omega(self, psi_n, R, mach=None, mach_species="deuterium",
                       voltage=None, psi_range=None, omega_profile=None):
    """geq's omega_profile hook re-splines omega on 101 midplane points every
    Picard step; an analytic profile is installed as itself instead."""
    if isinstance(omega_profile, Bump):
        for species in self.species_list:
            species.omega_spline = omega_profile
        self.omega_mach = True
        return None
    return _ORIG_SET_OMEGA(self, psi_n, R, mach=mach, mach_species=mach_species,
                           voltage=voltage, psi_range=psi_range,
                           omega_profile=omega_profile)


def tighten(all_species, psi_n, R, phi0):
    """phi_0 by a pointwise Newton on ln( Sum_ions Z n ) - ln n_e = 0 at every
    node with 0 <= psi_n <= 1, from geq's own N_s, T_s and omega. Strictly
    monotone in phi_0, and exactly linear for two species."""
    mask = (psi_n >= 0.0) & (psi_n <= 1.0)
    pn = np.clip(psi_n[mask], 0.0, 1.0)
    r2w2 = all_species.get_omega(pn) ** 2 * np.asarray(R)[mask] ** 2
    phi = np.array(phi0, dtype=float, copy=True)
    p = phi[mask]
    sp = all_species.species_list
    logN = [np.log(s.get_flux_density(pn)) for s in sp]
    T = [s.get_temp(pn) * qe for s in sp]
    for _ in range(100):
        ions, dions, elec, delec = [], [], None, None
        for k, s in enumerate(sp):
            ln = logN[k] - s.charge * p / T[k] + s.mass * r2w2 / (2.0 * T[k])
            if s.name == "electron":
                elec, delec = ln, -s.charge / T[k]
            else:
                ions.append(np.log(s.z) + ln)
                dions.append(-s.charge / T[k])
        stack = np.stack(ions)
        lse = logsumexp(stack, axis=0)
        w = np.exp(stack - lse)
        g = lse - elec
        dg = np.sum(w * np.stack(dions), axis=0) - delec
        step = -g / dg
        p = p + step
        if np.max(np.abs(step)) <= 1e-14 * max(1.0, np.max(np.abs(p))):
            break
    phi[mask] = p
    return phi


def qn_residual(all_species, psi_n, R, phi0):
    """Sum_s Z_s n_s/Sum_s |Z_s| n_s on the plasma nodes, from geq's densities."""
    mask = (psi_n >= 0.0) & (psi_n <= 1.0)
    pn = np.clip(psi_n[mask], 0.0, 1.0)
    r = np.asarray(R)[mask]
    ph = np.asarray(phi0)[mask]
    num, den = 0.0, 0.0
    for s in all_species.species_list:
        n = s.get_flux_density(pn) * np.exp(-s.get_potential(pn, r, phi0=ph))
        num = num + s.z * n
        den = den + abs(s.z) * n
    return num / den


def tightened_qn_global(self, psi_n, R=None, Z=None, maxiter=100,
                        verbose=False, multiwell=False):
    phi = _ORIG_QN_GLOBAL(self, psi_n, R, Z, maxiter=maxiter, verbose=verbose,
                          multiwell=multiwell)
    before = qn_residual(self, psi_n, R, phi)
    # geq's pin: phi_0 at the midplane node nearest psi_n = 0.5.
    Z0 = np.argmin(abs(Z[0, :]))
    pin = (np.argmin(np.abs(psi_n[:, Z0] - 0.5)), Z0)
    inside = (psi_n >= 0.0) & (psi_n <= 1.0)
    nopin = np.array(inside, copy=True)
    nopin[pin] = False
    before_nopin = qn_residual(self, np.where(nopin, psi_n, -1.0), R, phi)
    if STATE["tighten"]:
        phi = tighten(self, psi_n, R, phi)
    after = qn_residual(self, psi_n, R, phi)
    STATE["last_qn"] = dict(
        before_rms=float(np.sqrt(np.mean(before ** 2))),
        before_max=float(np.max(np.abs(before))),
        before_max_nopin=float(np.max(np.abs(before_nopin))),
        after_rms=float(np.sqrt(np.mean(after ** 2))),
        after_max=float(np.max(np.abs(after))),
    )
    return phi


AllSpecies.set_flux_density = held_set_flux_density
AllSpecies.set_omega = analytic_set_omega
AllSpecies.solve_quasineutrality_global = tightened_qn_global


# ---------------------------------------------------------------------------
# The state
# ---------------------------------------------------------------------------


def build_species(omega_peak):
    knots = np.linspace(0.0, 1.0, 129)
    density = Bump(N_PEAK, N_EDGE)
    species = []
    for name in NAMES:
        c = SPECIES[name]
        temp = Bump(c["peak"], c["edge"])
        s = SingleSpecies(name=name, z=c["z"], mass=c["mass"], psi_n=knots,
                          temp=temp(knots), flux_density=density(knots),
                          density_Z0=density(knots))
        # The analytic profiles REPLACE the splines geq just built.
        s._temp_spline = temp
        s._flux_density = Bump(N_PEAK, N_EDGE, log=True)
        s._density_Z0 = density
        s._density_eval = True
        species.append(s)
    all_species = AllSpecies(species_list=species)
    all_species.zeff_dict = {}
    omega = Bump(omega_peak, OMEGA_EDGE_FRACTION * omega_peak)
    return all_species, omega


def make_machine():
    return make_mirror_machine_from_field(
        field=FIELD, L=COIL_LENGTH, ri=COIL_INNER_RADIUS, dz=COIL_SEPARATION,
        Nr=1, Nz=10, rwall=[0.0, 0.5, 0.5, 0.0],
        zwall=[-COIL_SEPARATION * 1.5, -COIL_SEPARATION * 1.5,
               COIL_SEPARATION * 1.5, COIL_SEPARATION * 1.5])


def build(n, boundary, omega_peak):
    machine, info = make_machine()
    eq = Equilibrium(tokamak=machine, Rmin=R_MIN, Rmax=R_MAX, Zmin=-Z_EXTENT,
                     Zmax=Z_EXTENT, nx=n, ny=n,
                     boundary=freeBoundary if boundary == "direct"
                     else freeBoundaryHagenow,
                     psi=np.zeros((n, n)), order=4)
    eq._find_critical = False
    machine_psi = eq.tokamak.calcPsiFromGreens(eq._pgreen)
    eq._updatePlasmaPsi(np.zeros((n, n)))
    all_species, omega = build_species(omega_peak)
    psi_axis = eq.psiRZ(R_INNER, 0.0)
    psi_bndry = eq.psiRZ(R_OUTER, 0.0)
    bndry = dict(r_inner=R_INNER, z_inner=0.0, r_outer=R_OUTER, z_outer=0.0,
                 mach=None, voltage=None, omega_profile=omega)
    all_species.set_omega(psi_n=None, R=None, omega_profile=omega)
    profiles = ProfilesCentrifugalMirror(all_species=all_species)
    return eq, machine_psi, profiles, psi_axis, psi_bndry, bndry, info


def check_hold(profiles):
    """The hold, asserted rather than trusted: after a solve every species
    still carries the prescribed analytic profiles."""
    for s in profiles.all_species.species_list:
        if not (isinstance(s._flux_density, Bump) and isinstance(s._temp_spline, Bump)
                and isinstance(s.omega_spline, Bump) and s._density_eval):
            raise RuntimeError(f"geq replaced a prescribed profile of {s.name}")


def solve(n, boundary, omega_peak, rtol, maxits, cache=None, tight=True):
    """One converged state. Returns the Equilibrium, the profiles and a record."""
    tag = f"{SHAPE['name']}_{boundary}_{n}_{'tight' if tight else 'loose'}_w{omega_peak:.6e}_r{rtol:.1e}"
    eq, machine_psi, profiles, pa, pb, bndry, info = build(n, boundary, omega_peak)
    path = os.path.join(cache, tag + ".npz") if cache else None
    if path and os.path.exists(path):
        data = np.load(path, allow_pickle=True)
        eq._updatePlasmaPsi(np.array(data["plasma_psi"]))
        record = data["record"].item()
        record["cached"] = True
    else:
        STATE["tighten"] = tight
        STATE["rebuilds_refused"] = 0
        buf = io.StringIO()
        t0 = time.time()
        with contextlib.redirect_stdout(buf):
            solve_equilibrium(eq, profiles, machine_psi, solve_method="picard",
                              maxits=maxits, rtol=rtol, atol=0.0,
                              psi_axis=pa, psi_bndry=pb, bndry=bndry)
        wall = time.time() - t0
        STATE["tighten"] = True
        lines = [l for l in buf.getvalue().splitlines() if l.startswith("iter")]
        relchange = [float(l.split("psi_relchange=")[1].split(",")[0]) for l in lines]
        record = dict(n=n, boundary=boundary, tight=tight, wall=wall,
                      iterations=len(lines), relchange=relchange[-1],
                      qn=dict(STATE["last_qn"]),
                      rebuilds_refused=STATE["rebuilds_refused"],
                      coil_currents=[float(c.current) for _, c in eq.tokamak.coils],
                      vacuum_current=float(info["current_per_filament_A"]),
                      cached=False)
        check_hold(profiles)
        if path:
            os.makedirs(cache, exist_ok=True)
            np.savez(path, plasma_psi=eq.plasma_psi, record=np.array(record, dtype=object))
    record["psi_axis"] = float(eq.psiRZ(R_INNER, 0.0))
    record["psi_bndry"] = float(eq.psiRZ(R_OUTER, 0.0))
    stride = (n - 1) // (COARSE - 1)
    record["psi65"] = np.asarray(eq.psi())[::stride, ::stride].copy()
    record["R65"] = np.asarray(eq.R)[::stride, ::stride].copy()
    record["Z65"] = np.asarray(eq.Z)[::stride, ::stride].copy()
    record["plasma65"] = np.asarray(eq.plasma_psi)[::stride, ::stride].copy()
    return eq, profiles, record


# ---------------------------------------------------------------------------
# Diagnostics on a converged state
# ---------------------------------------------------------------------------


def midplane_current(eq, profiles, pa, pb, R):
    """geq's own Jtor at midplane radii R, phi_0 solved pointwise there."""
    R = np.atleast_1d(np.asarray(R, dtype=float))
    Z = np.zeros_like(R)
    psi = eq.psiRZ(R, Z)
    pn = (psi - pa) / (pb - pa)
    phi = tighten(profiles.all_species, pn, R, np.zeros_like(R))
    return profiles.Jtor(R, Z, psi, psi_bndry=pb, psi_axis=pa, phi0=phi), pn


def edge_currents(eq, profiles, pa, pb, offsets=(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)):
    """J_phi/peak on the midplane at psi_n = d and 1 - d."""
    Rs = np.linspace(R_INNER, R_OUTER, 8001)
    J, _ = midplane_current(eq, profiles, pa, pb, Rs)
    peak = float(np.max(np.abs(J)))

    def radius_of(target):
        f = lambda r: (eq.psiRZ(r, 0.0) - pa) / (pb - pa) - target
        return brentq(f, R_INNER - 0.02, R_OUTER + 0.02, xtol=1e-15, rtol=1e-15)

    out = {}
    for d in offsets:
        for label, target in (("inner", d), ("outer", 1.0 - d)):
            Jd, _ = midplane_current(eq, profiles, pa, pb, [radius_of(target)])
            out[(label, d)] = float(Jd[0]) / peak
    return peak, out


def machine_record(eq, info):
    filaments = []
    for name, coil in eq.tokamak.coils:
        for r_fil, z_fil in zip(coil.Rfil, coil.Zfil):
            filaments.append((float(r_fil), float(z_fil), float(coil.current)))
    return filaments


# ---------------------------------------------------------------------------
# The gauge transfer, R_ref and the tables
# ---------------------------------------------------------------------------


def species_arrays(psi_n, omega_peak):
    """N_s, T_s (J), dlnN_s/dpsi_n, dT_s/dpsi_n, from the analytic profiles."""
    density = Bump(N_PEAK, N_EDGE)
    out = {}
    for name in NAMES:
        c = SPECIES[name]
        temp = Bump(c["peak"], c["edge"])
        out[name] = dict(mass=c["mass"], z=c["z"],
                         flux_density=density.value(psi_n),
                         temp=temp.value(psi_n) * qe,
                         dlnN_dpsi=density.slope(psi_n) / density.value(psi_n),
                         dT_dpsi=temp.slope(psi_n) * qe)
    omega = Bump(omega_peak, OMEGA_EDGE_FRACTION * omega_peak)
    return out, omega.value(psi_n), omega.slope(psi_n)


def log_density_closed_form(R, psi_n, omega_peak):
    """ln n_s( R, psi_n ) for the two-species set, phi_0 in closed form:
    ln N_D - e phi/T_D + m_D w^2 R^2/2T_D = ln N_e + e phi/T_e + m_e w^2 R^2/2T_e."""
    arrays, omega, _ = species_arrays(psi_n, omega_peak)
    e, d = arrays["electron"], arrays["deuterium"]
    c = omega ** 2 * R ** 2 / 2.0
    phi = (np.log(d["flux_density"] / e["flux_density"])
           + c * (d["mass"] / d["temp"] - e["mass"] / e["temp"])) \
        / (qe * (1.0 / d["temp"] + 1.0 / e["temp"]))
    out = {}
    for name, a in arrays.items():
        out[name] = (np.log(a["flux_density"]) - a["z"] * qe * phi / a["temp"]
                     + a["mass"] * c / a["temp"])
    return out


def optimal_reference_radius(R, psi_n, omega_peak):
    """R_ref minimising max_s max_nodes | ln( n_s( R, psi )/n_s0( psi ) ) | --
    the largest factor between a density and its table value -- over the
    plasma nodes of the converged state."""
    here = log_density_closed_form(R, psi_n, omega_peak)

    def worst(r_ref):
        there = log_density_closed_form(np.full_like(R, r_ref), psi_n, omega_peak)
        return max(np.max(np.abs(here[s] - there[s])) for s in NAMES)

    scan = np.linspace(0.15, 0.40, 251)
    best = scan[np.argmin([worst(r) for r in scan])]
    res = minimize_scalar(worst, bounds=(best - 0.002, best + 0.002),
                          method="bounded", options=dict(xatol=1e-7))
    exponents = {}
    for r_ref in (res.x, round(res.x, 3)):
        there = log_density_closed_form(np.full_like(R, r_ref), psi_n, omega_peak)
        exponents[r_ref] = {s: (float(np.min(here[s] - there[s])),
                                float(np.max(here[s] - there[s]))) for s in NAMES}
    return float(res.x), exponents


def tables(knots, omega_peak, r_ref):
    """Every MEQ table on `knots` points of Psi = 1 - psi_n, ascending in Psi.
    Returns {file stem: ( Psi, f, df/dPsi )}."""
    Psi = np.linspace(0.0, 1.0, knots)
    return tables_at(Psi, omega_peak, r_ref)


def tables_at(Psi, omega_peak, r_ref):
    psi_n = 1.0 - Psi
    arrays, omega, domega = species_arrays(psi_n, omega_peak)
    n0, dn0, _ = transfer_to_reference_radius(arrays, NAMES, omega, domega, r_ref)
    out = {"omega": (Psi, omega, -domega)}
    for name in NAMES:
        a = arrays[name]
        out[f"T_{name}"] = (Psi, a["temp"], -a["dT_dpsi"])
        out[f"n_{name}"] = (Psi, n0[name], -dn0[name])
        out[f"Nraw_{name}"] = (Psi, a["flux_density"],
                               -a["dlnN_dpsi"] * a["flux_density"])
    return out


def midplane_crossings(eq, pa, pb, psi_n):
    """R_0( psi_n ), where each surface of geq's solved state crosses z = 0, by
    a root find on psiRZ along the midplane, and dR_0/dpsi_n = 1/( dpsi_n/dR )
    there, with dpsi/dR = R B_z read off geq's own field rather than
    differenced. MIRROR-PLAN.md section 5's crossing, on geq's state."""
    span = pb - pa
    lo, hi = R_INNER - 0.02, R_OUTER + 0.02
    r0 = np.empty_like(psi_n)
    for i, target in enumerate(psi_n):
        f = lambda r, t=target: (eq.psiRZ(r, 0.0) - pa) / span - t
        r0[i] = brentq(f, lo, hi, xtol=1e-15, rtol=1e-15)
    slope = np.array([r * float(eq.Bz(r, 0.0)) for r in r0]) / span
    return r0, 1.0 / slope


def transfer_to_midplane(psi_n, omega_peak, r0, dr0):
    """The densities MEQ's MIDPLANE gauge wants, n_s( R_0( psi ), psi ), and
    their slope d/dpsi_n along the crossing.

    geq's N_s are its ZERO-POTENTIAL densities, held as prescribed flux
    functions, so this is NOT the near-identity MIRROR-PLAN.md section 5
    expected: it is geq's quasineutrality solved at each surface's own crossing,
    transfer_to_reference_radius() point by point with R_ref = R_0, plus the
    crossing's own motion,

        dn_s/dpsi_n = dn_s/dpsi_n |_R + dn_s/dR |_psi R_0',
        dn_s/dR |_psi = n_s ( m_s omega^2 R - Z_s e dphi_0/dR )/T_s,

    with dphi_0/dR by implicit differentiation of Sum_s Z_s n_s = 0 at fixed
    psi. It depends on R_0, which is geq's ANSWER, so these tables are an
    output and are not tracked."""
    arrays, omega, domega = species_arrays(psi_n, omega_peak)
    n0 = {name: np.zeros_like(psi_n) for name in NAMES}
    dn0 = {name: np.zeros_like(psi_n) for name in NAMES}
    for i in range(len(psi_n)):
        at = {name: {k: (v[i:i + 1] if isinstance(v, np.ndarray) else v)
                     for k, v in arrays[name].items()} for name in NAMES}
        n, dn, phi = transfer_to_reference_radius(at, NAMES, omega[i:i + 1],
                                                  domega[i:i + 1], float(r0[i]))
        wsq = omega[i] ** 2
        num, den = 0.0, 0.0
        for name in NAMES:
            a = at[name]
            num += a["z"] * n[name][0] * a["mass"] * wsq * r0[i] / a["temp"][0]
            den += a["z"] ** 2 * qe * n[name][0] / a["temp"][0]
        dphi_dr = num / den
        for name in NAMES:
            a = at[name]
            dn_dr = n[name][0] * (a["mass"] * wsq * r0[i]
                                  - a["z"] * qe * dphi_dr) / a["temp"][0]
            n0[name][i] = n[name][0]
            dn0[name][i] = dn[name][0] + dn_dr * dr0[i]
    return n0, dn0


def midplane_tables(eq, pa, pb, omega_peak, knots):
    """The midplane gauge's density tables on `knots` points of Psi, ascending,
    and the Hermite error of each against the transfer at the midpoints."""
    Psi = np.linspace(0.0, 1.0, knots)
    mid = 0.5 * (Psi[1:] + Psi[:-1])
    out, exact = {}, {}
    for where, x in (("knots", Psi), ("mid", mid)):
        psi_n = 1.0 - x
        r0, dr0 = midplane_crossings(eq, pa, pb, psi_n)
        n0, dn0 = transfer_to_midplane(psi_n, omega_peak, r0, dr0)
        for name in NAMES:
            (out if where == "knots" else exact)[f"nmid_{name}"] = (x, n0[name], -dn0[name])
    errors = {}
    for stem, (xk, fk, dk) in out.items():
        v, s = hermite(exact[stem][0], xk, fk, dk)
        ev, es = exact[stem][1], exact[stem][2]
        errors[stem] = (float(np.max(np.abs(v - ev) / np.abs(ev))),
                        float(np.max(np.abs(s - es)) / np.max(np.abs(es))))
    return out, errors



def hermite(x, xk, fk, dk):
    """meq::SplineProfile's interpolant: the cubic Hermite through ( f, f' )."""
    i = np.clip(np.searchsorted(xk, x, side="right") - 1, 0, len(xk) - 2)
    h = xk[i + 1] - xk[i]
    t = (x - xk[i]) / h
    h00 = 2 * t ** 3 - 3 * t ** 2 + 1
    h10 = t ** 3 - 2 * t ** 2 + t
    h01 = -2 * t ** 3 + 3 * t ** 2
    h11 = t ** 3 - t ** 2
    value = h00 * fk[i] + h10 * h * dk[i] + h01 * fk[i + 1] + h11 * h * dk[i + 1]
    d00 = (6 * t ** 2 - 6 * t) / h
    d10 = 3 * t ** 2 - 4 * t + 1
    d01 = (-6 * t ** 2 + 6 * t) / h
    d11 = 3 * t ** 2 - 2 * t
    slope = d00 * fk[i] + d10 * dk[i] + d01 * fk[i + 1] + d11 * dk[i + 1]
    return value, slope


def interpolation_errors(knots, omega_peak, r_ref):
    """Max pointwise relative error of each table's Hermite interpolant, and of
    its slope relative to the table's largest slope, at three interior points
    of every interval."""
    tk = tables(knots, omega_peak, r_ref)
    xk = tk["omega"][0]
    x = np.sort(np.concatenate([xk[:-1] + f * (xk[1] - xk[0])
                                for f in (0.25, 0.5, 0.75)]))
    exact = tables_at(x, omega_peak, r_ref)
    errors = {}
    for stem, (_, fk, dk) in tk.items():
        v, s = hermite(x, xk, fk, dk)
        ev, es = exact[stem][1], exact[stem][2]
        errors[stem] = (float(np.max(np.abs(v - ev) / np.abs(ev))),
                        float(np.max(np.abs(s - es)) / np.max(np.abs(es))))
    return tk, errors


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_table(path, x, value, slope, header):
    with open(path, "w") as fh:
        for line in header.strip().splitlines():
            fh.write(f"# {line}".rstrip() + "\n")
        fh.write("#\n#    Psi            f( Psi )        df/dPsi\n")
        for a, b, c in zip(x, value, slope):
            fh.write(f"{a:.17e} {b:.17e} {c:.17e}\n")


def write_psi(path, record, what):
    with open(path, "w") as fh:
        fh.write(f"# geq's total psi, Wb/rad, {what}, on the {COARSE} x {COARSE}\n")
        fh.write("# nodes every resolution shares ( R in [ 0.05, 0.63 ], |z| <= 1.5 ).\n")
        fh.write(f"# psi_axis  = psiRZ( 0.10, 0 ) = {record['psi_axis']:.17e}\n")
        fh.write(f"# psi_bndry = psiRZ( 0.35, 0 ) = {record['psi_bndry']:.17e}\n")
        fh.write("# See export_rotating.py and README.md.\n")
        fh.write("#\n#    R              z               psi\n")
        R, Z, P = record["R65"], record["Z65"], record["psi65"]
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                fh.write(f"{R[i, j]:.17e} {Z[i, j]:.17e} {P[i, j]:.17e}\n")


def git_commit(path):
    try:
        return subprocess.check_output(["git", "-C", path, "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:  # noqa: BLE001 -- provenance is best effort
        return "unknown"


def git_dirty(path):
    try:
        out = subprocess.check_output(["git", "-C", path, "status", "--porcelain",
                                       "--untracked-files=no"], text=True)
        return "dirty" if out.strip() else "clean"
    except Exception:  # noqa: BLE001
        return "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=None,
                    help="geq's answer, psi_*.dat and solution.txt; default "
                         "runs/geq-benchmark/reference/rotating")
    ap.add_argument("--tables", default=None,
                    help="the problem, meta.txt, filaments.dat and the tables; "
                         "default --out. tools/geq-benchmark/reference/rotating "
                         "refreshes the tracked copy")
    ap.add_argument("--mach", type=float, default=2.0,
                    help="geq's Mach, omega R/sqrt( 2 T_D/m_D ), at the "
                         "profiles' peak and R = R_MACH")
    ap.add_argument("--r-mach", type=float, default=None,
                    help="the R that Mach is set at; default per shape, the "
                         "converged midplane radius of the peak")
    ap.add_argument("--resolutions", default="129,257,513")
    ap.add_argument("--control", type=int, default=257,
                    help="resolution of the shipped-Hagenow control; 0 skips it")
    ap.add_argument("--rtol", type=float, default=2e-9)
    ap.add_argument("--maxits", type=int, default=400)
    ap.add_argument("--knots", type=int, default=None,
                    help="table knots; default the smallest power of two plus one "
                         "whose Hermite error is below 1e-10 in every table")
    ap.add_argument("--R-ref", type=float, default=None,
                    help="default: the exponent-range optimum, rounded to 1 mm")
    ap.add_argument("--loose-study", type=int, default=129,
                    help="also solve WITHOUT the phi_0 tightening at this "
                         "resolution, to measure what it moves; 0 skips")
    ap.add_argument("--cache", default=None,
                    help="directory for per-solve .npz files, reused on re-runs")
    ap.add_argument("--shape", choices=tuple(SHAPES), default="skew",
                    help="profile shape; the reference is the skew")
    args = ap.parse_args(argv)
    args.out = args.out or answers_dir("rotating")
    args.tables = args.tables or args.out
    SHAPE["name"] = args.shape
    x_peak = SHAPES[args.shape][0]
    R_MACH["value"] = args.r_mach if args.r_mach is not None else SHAPES[args.shape][1]

    omega_peak = omega_peak_for(args.mach)
    ns = sorted(int(x) for x in args.resolutions.split(","))
    for n in ns + ([args.control] if args.control else []):
        if (n - 1) % (COARSE - 1):
            raise SystemExit(f"resolution {n} does not contain the {COARSE}-grid")

    # --- the solves ---------------------------------------------------------
    runs, states = {}, {}
    for n in ns:
        eq, prof, rec = solve(n, "direct", omega_peak, args.rtol, args.maxits, args.cache)
        runs[("direct", n)] = rec
        states[("direct", n)] = (eq, prof)
        print(f"direct   {n:4d}^2  {rec['wall']:8.1f} s  {rec['iterations']:3d} its"
              f"  relchange {rec['relchange']:.2e}  psi_axis {rec['psi_axis']:.12e}"
              f"  psi_bndry {rec['psi_bndry']:.12e}"
              f"{'  (cached)' if rec['cached'] else ''}", flush=True)
    if args.control:
        eq, prof, rec = solve(args.control, "hagenow", omega_peak, args.rtol,
                              args.maxits, args.cache)
        runs[("hagenow", args.control)] = rec
        print(f"hagenow  {args.control:4d}^2  {rec['wall']:8.1f} s  {rec['iterations']:3d} its"
              f"  relchange {rec['relchange']:.2e}  psi_axis {rec['psi_axis']:.12e}"
              f"  psi_bndry {rec['psi_bndry']:.12e}", flush=True)
    if args.loose_study:
        _, _, rec = solve(args.loose_study, "direct", omega_peak, args.rtol,
                          args.maxits, args.cache, tight=False)
        runs[("loose", args.loose_study)] = rec
        print(f"loose    {args.loose_study:4d}^2  {rec['wall']:8.1f} s"
              f"  {rec['iterations']:3d} its", flush=True)

    finest = ns[-1]
    eqf, proff = states[("direct", finest)]
    recf = runs[("direct", finest)]
    pa, pb = recf["psi_axis"], recf["psi_bndry"]
    span = pb - pa

    # --- (a) J_phi at the labels -------------------------------------------
    print("\n(a) J_phi on the midplane just inside each label, relative to peak")
    edges = {}
    for n in ns:
        eq, prof = states[("direct", n)]
        r = runs[("direct", n)]
        peak, e = edge_currents(eq, prof, r["psi_axis"], r["psi_bndry"])
        edges[n] = (peak, e)
        print(f"  {n:4d}^2  peak |J| {peak:.4e} A/m^2")
        for d in sorted({k[1] for k in e}):
            print(f"     d = {d:.0e}   inner {e[('inner', d)]:+.3e}"
                  f"   outer {e[('outer', d)]:+.3e}")

    # The edge order j, J ~ d^j: a least-squares slope in log-log over every
    # offset, and over the three smallest, on the finest state.
    _, e = edges[finest]
    ds = np.array(sorted({k[1] for k in e}))
    edge_orders = {}
    for label in ("inner", "outer"):
        y = np.log(np.abs([e[(label, d)] for d in ds]))
        p_all = float(np.polyfit(np.log(ds), y, 1)[0])
        p_small = float(np.polyfit(np.log(ds[:3]), y[:3], 1)[0])
        edge_orders[label] = (p_all, p_small)
        print(f"  fitted exponent of |J| against d, {label}: {p_all:.4f} over all,"
              f" {p_small:.4f} over the three smallest")

    # --- (b) self-convergence ------------------------------------------------
    print("\n(b) direct self-convergence on the 65-grid nodes, against", finest)
    ref = runs[("direct", finest)]
    plasma_scale = float(np.max(np.abs(ref["plasma65"])))
    conv = {}
    for n in ns[:-1]:
        r = runs[("direct", n)]
        d = r["psi65"] - ref["psi65"]
        conv[n] = (float(np.max(np.abs(d))), r["psi_axis"] - ref["psi_axis"],
                   r["psi_bndry"] - ref["psi_bndry"])
        print(f"  {n:4d} vs {finest}: max|dpsi| {conv[n][0]:.3e}"
              f"  dpsi_axis {conv[n][1]:+.3e}  dpsi_bndry {conv[n][2]:+.3e}")
    if len(ns) >= 3:
        a, b, c = (runs[("direct", n)] for n in ns[-3:])
        e1 = np.max(np.abs(a["psi65"] - b["psi65"]))
        e2 = np.max(np.abs(b["psi65"] - c["psi65"]))
        print(f"  successive max|dpsi| {e1:.3e}, {e2:.3e}: order {np.log2(e1 / e2):.2f}")
        for k in ("psi_axis", "psi_bndry"):
            e1, e2 = a[k] - b[k], b[k] - c[k]
            print(f"  successive {k} {e1:+.3e}, {e2:+.3e}: ratio {e1 / e2:.2f}")
    print(f"  max |plasma psi| {plasma_scale:.4e}, label span {span:.4e}")

    # --- (c) the shipped-Hagenow bias ----------------------------------------
    if args.control:
        h = runs[("hagenow", args.control)]
        dref = runs[("direct", args.control)] if ("direct", args.control) in runs else ref
        d = h["psi65"] - dref["psi65"]
        print(f"\n(c) shipped Hagenow minus direct at {args.control}^2:"
              f" max|dpsi| {np.max(np.abs(d)):.3e}  dpsi_axis"
              f" {h['psi_axis'] - dref['psi_axis']:+.3e}  dpsi_bndry"
              f" {h['psi_bndry'] - dref['psi_bndry']:+.3e}")

    # --- (d) phi_0 residual --------------------------------------------------
    print("\n(d) Sum Z n/Sum |Z| n on the plasma nodes, final global solve")
    for key, r in runs.items():
        q = r["qn"]
        print(f"  {key[0]:8s} {key[1]:4d}^2  geq: rms {q['before_rms']:.3e} max"
              f" {q['before_max']:.3e} (max off the pin node {q['before_max_nopin']:.3e})"
              f"   tightened: rms {q['after_rms']:.3e} max {q['after_max']:.3e}")
    if args.loose_study:
        lo = runs[("loose", args.loose_study)]
        ti = runs[("direct", args.loose_study)]
        d = lo["psi65"] - ti["psi65"]
        print(f"  untightened fixed point minus tightened at {args.loose_study}^2:"
              f" max|dpsi| {np.max(np.abs(d)):.3e}  dpsi_axis"
              f" {lo['psi_axis'] - ti['psi_axis']:+.3e}  dpsi_bndry"
              f" {lo['psi_bndry'] - ti['psi_bndry']:+.3e}")

    # --- R_ref ---------------------------------------------------------------
    R, Z = np.asarray(eqf.R), np.asarray(eqf.Z)
    psi_n_2d = (np.asarray(eqf.psi()) - pa) / span
    inside = (psi_n_2d >= 0.0) & (psi_n_2d <= 1.0)
    r_opt, exps = optimal_reference_radius(R[inside], psi_n_2d[inside], omega_peak)
    r_ref = args.R_ref if args.R_ref is not None else round(r_opt, 3)
    if r_ref not in exps:
        there = log_density_closed_form(np.full_like(R[inside], r_ref),
                                        psi_n_2d[inside], omega_peak)
        here = log_density_closed_form(R[inside], psi_n_2d[inside], omega_peak)
        exps[r_ref] = {s: (float(np.min(here[s] - there[s])),
                           float(np.max(here[s] - there[s]))) for s in NAMES}
    print(f"\nR_ref: optimum {r_opt:.6f} m, used {r_ref:.6f} m")
    for rr, e in exps.items():
        print(f"  R_ref {rr:.6f}: " + "  ".join(
            f"{s} [{lo:+.3f}, {hi:+.3f}]" for s, (lo, hi) in e.items()))

    # --- the tables and their interpolation error ---------------------------
    print("\nHermite interpolation error, pointwise relative (value) and relative"
          " to max|f'| (slope)")
    candidates = [args.knots] if args.knots else [257, 513, 1025, 2049, 4097]
    chosen = None
    for k in candidates:
        tk, err = interpolation_errors(k, omega_peak, r_ref)
        worst = max(v for v, _ in err.values())
        print(f"  {k:5d} knots: " + "  ".join(f"{s} {v:.1e}/{d:.1e}"
                                             for s, (v, d) in err.items()))
        chosen = (k, tk, err)
        if worst < 1e-10 and not args.knots:
            break
    knots, tk, err = chosen

    # Every slope at both labels, which is the precondition F = 0 there rests on.
    end_slopes = {s: (float(v[2][0]), float(v[2][-1])) for s, v in tk.items()}

    # --- the midplane gauge's tables, on geq's finest crossings -------------
    # MIRROR-PLAN.md section 5: n_s0 is the density where each surface crosses
    # z = 0. On geq's held N_s that needs geq's solved R_0( psi ), so these are
    # made from geq's answer and go beside it, untracked.
    mk, merr = midplane_tables(eqf, pa, pb, omega_peak, knots)
    print("\nthe midplane gauge's tables, on geq's " f"{finest}^2 crossings:")
    midplane_from = {}
    for name in NAMES:
        nmid = mk[f"nmid_{name}"][1]
        midplane_from[name] = (
            float(np.max(np.abs(np.log(nmid / tk[f"Nraw_{name}"][1])))),
            float(np.max(np.abs(np.log(nmid / tk[f"n_{name}"][1])))))
        v, d = merr[f"nmid_{name}"]
        print(f"  {name:10s} Hermite {v:.1e}/{d:.1e};  max | ln( n_mid/N_geq ) |"
              f" {midplane_from[name][0]:.4f},  max | ln( n_mid/n_R_ref ) |"
              f" {midplane_from[name][1]:.4f}")
    nmid_e, nmid_d = mk["nmid_electron"][1], mk["nmid_deuterium"][1]
    midplane_neutrality = float(np.max(np.abs(nmid_d - nmid_e) / (nmid_d + nmid_e)))
    print(f"  neutrality at the crossings {midplane_neutrality:.3e}")

    # Neutrality at R_ref, what meq::RotatingSource's constructor checks.
    ne, nd = tk["n_electron"][1], tk["n_deuterium"][1]
    neutrality = float(np.max(np.abs(nd - ne) / (nd + ne)))

    # Mach numbers on the converged state.
    Z0 = np.argmin(abs(Z[0, :]))
    r_half = float(np.interp(x_peak, psi_n_2d[:, Z0], R[:, Z0]))
    omega_half = omega_peak
    mach_geq = omega_half * r_half / np.sqrt(2.0 * TD_PEAK * qe / DEUTERON_MASS)
    mach_167 = np.sqrt(DEUTERON_MASS * omega_half ** 2 * r_half ** 2 / (TE_PEAK * qe))

    # --- write ---------------------------------------------------------------
    out, tables = args.out, args.tables
    os.makedirs(out, exist_ok=True)
    os.makedirs(tables, exist_ok=True)

    write_table(os.path.join(tables, "omega.dat"), *tk["omega"],
                "omega( Psi ) in rad/s, Psi = 1 - psi_n ( 1 at the inner label\n"
                "( 0.10, 0 ), 0 at the outer ( 0.35, 0 ) ). Rigid rotation, the\n"
                "analytic profile geq was run with; the slope is d/dPsi.")
    for name in NAMES:
        write_table(os.path.join(tables, f"T_{name}.dat"), *tk[f"T_{name}"],
                    f"T_{name}( Psi ) in JOULES, Psi = 1 - psi_n. The analytic\n"
                    "profile geq was run with ( in eV there ); slope d/dPsi.")
        write_table(os.path.join(tables, f"n_{name}.dat"), *tk[f"n_{name}"],
                    f"n_{name}( Psi ) in m^-3 ON THE CURVE R = R_ref = {r_ref:.17e},\n"
                    "MEQ's gauge phi_0( R_ref, psi ) = 0: geq's prescribed N_s carried\n"
                    "to R_ref by quasineutrality there, slope by implicit\n"
                    "differentiation. Psi = 1 - psi_n; slope d/dPsi.")
        write_table(os.path.join(tables, f"Nraw_{name}.dat"), *tk[f"Nraw_{name}"],
                    f"N_{name}( Psi ) IN geq'S OWN GAUGE: the prescribed flux function,\n"
                    "n_s = N_s exp( -( Z_s e phi_0 - m_s omega^2 R^2/2 )/T_s ).\n"
                    "NOT what meq::Species wants. Psi = 1 - psi_n; slope d/dPsi.")

    for n in ns:
        write_psi(os.path.join(out, f"psi_direct_{n}.dat"), runs[("direct", n)],
                  f"direct freeBoundary, {n} x {n}")
    for name in NAMES:
        write_table(os.path.join(out, f"nmid_{name}.dat"), *mk[f"nmid_{name}"],
                    f"n_{name}( Psi ) in m^-3 WHERE EACH FLUX SURFACE CROSSES z = 0, MEQ's\n"
                    "MIDPLANE gauge phi_0( R_0( psi ), psi ) = 0: geq's prescribed N_s\n"
                    f"carried by quasineutrality to the crossings R_0 of geq's {finest}^2\n"
                    "direct state, slope along the crossing. MADE FROM geq's ANSWER, so\n"
                    "an output beside psi_*.dat and not tracked. Psi = 1 - psi_n; d/dPsi.")
    if args.control:
        write_psi(os.path.join(out, f"psi_hagenow_{args.control}.dat"),
                  runs[("hagenow", args.control)],
                  f"SHIPPED von Hagenow ( eps = 1e-2 m ), {args.control} x "
                  f"{args.control} -- THE CONTROL, biased")

    machine, info = make_machine()
    filaments = machine_record(eqf, info)
    with open(os.path.join(tables, "filaments.dat"), "w") as fh:
        fh.write("# geq's 2 x 10 point filaments, FIXED ( control = False ), as\n")
        fh.write("# given: R z current-per-filament ( A ). The same machine as\n")
        fh.write("# reference/vacuum, i.e. meq::filamentStack( pack, 1, 10 ).\n")
        fh.write("#\n#    R              z               current\n")
        for r_fil, z_fil, cur in filaments:
            fh.write(f"{r_fil:.17e} {z_fil:.17e} {cur:.17e}\n")

    geq_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(sys.modules["geq.solver"].__file__))))
    fgs_dir = os.path.dirname(os.path.dirname(
        os.path.abspath(sys.modules["freegs4e.species"].__file__)))
    with open(os.path.join(tables, "meta.txt"), "w") as fh:
        fh.write("case rotating\n")
        fh.write(f"geq_commit {git_commit(geq_dir)} {git_dirty(geq_dir)}\n")
        fh.write(f"freegs4e_commit {git_commit(fgs_dir)} {git_dirty(fgs_dir)}\n")
        fh.write("# Units: psi in Wb/rad; R, z in m; T in J; n in m^-3; omega in\n")
        fh.write("# rad/s; mass in kg; charge a charge number. Tables are in\n")
        fh.write("# Psi = 1 - psi_n = ( psi - psi_bndry )/( psi_axis - psi_bndry ).\n")
        fh.write(f"mu0 {MU0:.17e}\n")
        fh.write(f"R_ref {r_ref:.17e}\n")
        fh.write(f"axis_point {R_INNER:.17e} 0.0\n")
        fh.write(f"boundary_point {R_OUTER:.17e} 0.0\n")
        fh.write(f"z_edge {Z_EXTENT:.17e}\n")
        fh.write("ggprime 0.0\n")
        fh.write(f"omega_peak {omega_peak:.17e}\n")
        fh.write(f"omega_edge {OMEGA_EDGE_FRACTION * omega_peak:.17e}\n")
        fh.write(f"Te_peak_eV {TE_PEAK:.17e}\nTe_edge_eV {TE_EDGE:.17e}\n")
        fh.write(f"TD_peak_eV {TD_PEAK:.17e}\nTD_edge_eV {TD_EDGE:.17e}\n")
        fh.write(f"N_peak {N_PEAK:.17e}\nN_edge {N_EDGE:.17e}\n")
        fh.write({"skew": "profile f_edge+(f_peak-f_edge)*C*psi_n^3*(1-psi_n)^4"
                          f" C {SKEW_C:.17e}\n",
                  "quartic": "profile f_edge+(f_peak-f_edge)*16*psi_n^2*(1-psi_n)^2\n",
                  "sextic": "profile f_edge+(f_peak-f_edge)*64*psi_n^3*(1-psi_n)^3\n",
                  }[args.shape])
        fh.write(f"profile_peak_psi_n {x_peak:.17e}\n")
        fh.write(f"mach_set_at_R {R_MACH['value']:.17e}\n")
        fh.write(f"knots {knots}\n")
        for s, (v, d) in err.items():
            fh.write(f"hermite_error_{s} {v:.3e} {d:.3e}\n")
        fh.write(f"neutrality_at_R_ref {neutrality:.3e}\n")
        fh.write("picard atol 0 rtol {:.3e} maxits {}\n".format(args.rtol, args.maxits))
        fh.write("boundary freeBoundary ( direct ); control freeBoundaryHagenow\n")
        fh.write("phi0 tightened_every_picard_step\n")
        fh.write("N_s held_as_prescribed_flux_functions\n")
        fh.write("grid R 0.05 0.63 z -1.5 1.5 common 65\n")
        # The machine, keyed as export_vacuum.py writes it.
        fh.write(f"field {FIELD:.17e}\n")
        fh.write(f"coil_inner_radius {info['ri_m']:.17e}\n")
        fh.write(f"coil_outer_radius {info['ro_m']:.17e}\n")
        fh.write(f"coil_length {info['L_m']:.17e}\n")
        fh.write(f"coil_centre_z {0.5 * info['dz_m']:.17e}\n")
        fh.write(f"pack_current {info['total_current_per_pack_A']:.17e}\n")
        fh.write(f"filaments_per_pack {info['nfil_per_pack']}\n")
        fh.write(f"species {len(NAMES)}\n")
        for name in NAMES:
            c = SPECIES[name]
            fh.write(f"name {name} mass {c['mass']:.17e} charge {c['z']:.17e}"
                     f" temperature T_{name}.dat density n_{name}.dat\n")

    # WHAT geq SOLVED, apart from what it was given. meta.txt and the tables are
    # the problem; this file and psi_*.dat are geq's answer to it, regenerated
    # by a run into runs/ and never committed.
    with open(os.path.join(out, "solution.txt"), "w") as fh:
        fh.write("# geq's converged state for meta.txt's problem. An OUTPUT of\n")
        fh.write("# export_rotating.py, beside psi_*.dat; not tracked.\n")
        fh.write(f"psi_axis {pa:.17e}\n")
        fh.write(f"psi_boundary {pb:.17e}\n")
        for (kind, n), r in sorted(runs.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            if kind == "loose":
                continue
            fh.write(f"labels_{kind}_{n} {r['psi_axis']:.17e} {r['psi_bndry']:.17e}\n")
        fh.write(f"R_ref_optimum {r_opt:.17e}\n")
        fh.write(f"mach_geq {mach_geq:.17e}\n")
        fh.write(f"mach_m167 {mach_167:.17e}\n")
        fh.write(f"mach_evaluated_at_R {r_half:.17e}\n")
        for label, (p_all, p_small) in edge_orders.items():
            fh.write(f"edge_order_{label} {p_all:.4f} {p_small:.4f}\n")
        fh.write(f"midplane_tables_on {finest}\n")
        fh.write(f"midplane_neutrality {midplane_neutrality:.3e}\n")
        for name in NAMES:
            v, d = merr[f"nmid_{name}"]
            fh.write(f"midplane_hermite_error_{name} {v:.3e} {d:.3e}\n")
            fh.write(f"midplane_max_abs_log_from_geq_{name} {midplane_from[name][0]:.17e}\n")
            fh.write(f"midplane_max_abs_log_from_R_ref_{name} {midplane_from[name][1]:.17e}\n")

    print(f"\nknots {knots}; slopes at Psi = 0, 1: " + "  ".join(
        f"{s} {a:+.1e}/{b:+.1e}" for s, (a, b) in end_slopes.items()))
    print(f"neutrality at R_ref {neutrality:.3e}")
    print(f"geq Mach at psi_n = {x_peak:.6f}: {mach_geq:.4f} ( R = {r_half:.5f} );"
          f" M-167's {mach_167:.4f}")
    print(f"coil currents unchanged: "
          f"{all(abs(c - recf['vacuum_current']) == 0.0 for c in recf['coil_currents'])}")
    print(f"N_s rebuilds refused in the last solve: {recf['rebuilds_refused']}")
    print(f"answer written to {out}")
    print(f"problem written to {tables}")
    against_tracked(tables, "rotating")
    return 0


if __name__ == "__main__":
    sys.exit(main())
