#!/usr/bin/env python
"""EXPORT A ROTATING-MIRROR STATE FROM geq, IN THE GAUGE meq::RotatingSource USES.

NO SOLVER RUNS ON EITHER SIDE. This writes a PRESCRIBED state -- a vacuum-field
mirror with species, rotation and the quasineutrality potential solved -- so that
MEQ's F and geq's mu0 R Jtor can be compared pointwise. The same shape of
question tools/freegs4e-benchmark/source_check.py asks of the MHD source, and it
separates the SOURCE from the discretisation and the iteration exactly.

WHAT IS BEING COMPARED. Both codes implement Abel et al. 2013,
Rep. Prog. Phys. 76 116201, eq (136), closed by its (96) and (97). They write it
differently and the two forms agree analytically:

    geq   mu0 R Jtor = mu0 R^2 Sum_s n_s [ T_s dlnN_s/dpsi
                                           + ( chi_s + 1 ) T_s dlnT_s/dpsi ]
                     + mu0 R^2 Sum_s m_s n_s R^2 omega domega/dpsi
    MEQ   F          = mu0 R^2 dp/dpsi|_r + g g',   p = Sum_s n_s T_s

Expanding MEQ's dp/dpsi|_r gives geq's bracket term for term ONCE the
Sum_s Z_s n_s d(e phi_0)/dpsi term is dropped -- and it drops by
quasineutrality, which is the cancellation docs/rotation.rst records. So a
disagreement here is an implementation fault in one of the two, not a
difference of formulation.

THE GAUGE IS THE WHOLE OF THE CONVERSION AND IT IS NOT A DETAIL.

    geq  pins phi_0 at ONE GLOBAL POINT, psi_n = 0.5 on the midplane, and its
         N_s( psi ) is defined by  n_s = N_s exp( [ m_s omega^2 R^2/2
                                                    - Z_s e phi_0 ]/T_s ).
    MEQ  pins phi_0( R_ref, psi ) = 0, ONE CONDITION PER FLUX SURFACE, and its
         n_s0( psi ) is the PHYSICAL density on the curve R = R_ref.

N_s absorbs the difference through exp( Z_s e delta( psi )/T_s ), which is a
DIFFERENT factor for every species -- so the two codes' tabulated densities mean
different things and no single rescaling relates them. Hand geq's N_s to MEQ
unchanged and it parses, solves, and converges to a different plasma.

THE TRANSFER IS geq'S OWN QUASINEUTRALITY, SOLVED AT R = R_ref. For each psi,

    Sum_s Z_s N_s( psi ) exp( [ m_s omega^2 R_ref^2/2 - Z_s e Phi ]/T_s ) = 0

has a unique root Phi = phi_0^geq( R_ref, psi ) -- the sum is strictly
decreasing in Phi -- and the densities it produces ARE MEQ's n_s0. This needs no
interpolation of geq's phi_0 field, works at any R_ref, and satisfies
Sum_s Z_s n_s0 = 0 to root-find precision, which is what
meq::RotatingSource's constructor checks.

IT IS NOT CIRCULAR, and the reason is worth stating. The transfer fixes the
densities on ONE curve; what is then compared is MEQ's evaluation everywhere
else -- its own phi_0( R, psi ) solve at every other radius, the psi-derivative
chain through N, T and omega, and the C( psi ) closure. A shared sign error in
the exponent would cancel at R_ref and survive off it, which is precisely where
the comparison lives.

ITS psi-DERIVATIVE IS BY IMPLICIT DIFFERENTIATION, NOT BY DIFFERENCING. F is
mu0 R^2 dp/dpsi, so the derivative column of the density table is not a
convenience -- it IS the thing being compared. Differentiating the quasineutrality
condition gives Phi' in closed form, so the tables carry exact slopes and
meq::Profile's Hermite cubic reproduces geq's own splines rather than a
difference of them.

USAGE

    PYTHONPATH=/home/ian/projects/geq/src \
    tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_mirror.py \
        --case two --out tools/geq-benchmark/reference/two

Written per case: meta.txt, omega.dat, T_<name>.dat, n_<name>.dat, points.dat.
See README.md for the file formats and for what the numbers came out at.
"""

import argparse
import os
import sys

import numpy as np
from scipy import constants
from scipy.optimize import brentq

import matplotlib

matplotlib.use("Agg")

from freegs4e.boundary import freeBoundaryHagenow
from freegs4e.equilibrium import Equilibrium
from freegs4e.jtor import ProfilesCentrifugalMirror
from freegs4e.species import AllSpecies, SingleSpecies

from geq.machine import make_mirror_machine_from_field
from geq.profiles import generate_density_z0_profile, generate_temperature_profile

qe = constants.elementary_charge
MU0 = 4.0e-7 * np.pi

ELECTRON_MASS = constants.electron_mass
DEUTERON_MASS = constants.physical_constants["deuteron mass"][0]
TRITON_MASS = constants.physical_constants["triton mass"][0]

CATALOGUE = {
    "electron": dict(z=-1.0, mass=ELECTRON_MASS),
    "deuterium": dict(z=1.0, mass=DEUTERON_MASS),
    "tritium": dict(z=1.0, mass=TRITON_MASS),
    "carbon": dict(z=6.0, mass=6.0 * DEUTERON_MASS),
}

# The mirror machine, as geq's own run.py and tests/conftest.py pose it.
R_INNER = 0.1
R_OUTER = 0.35
N_PSI_SPECIES = 129


# ---------------------------------------------------------------------------
# The state
# ---------------------------------------------------------------------------


def build_state(case, nx, ny, omega_axis, te_axis, ti_axis, t_edge, ne_axis,
                ne_edge, z_extent=None):
    """A vacuum-field rotating mirror with species, omega and phi_0 solved.

    THE ROTATION IS HANDED IN AS omega( psi_n ) RATHER THAN AS A MACH NUMBER OR
    A VOLTAGE, through geq's own omega_profile hook. Both of the other routes
    re-derive omega from the temperatures at every Picard step, so the rotation
    would be a DERIVED quantity on one side and an input on the other; the hook
    is the only way to give the two codes the same omega( psi ) by construction.
    """
    names = {
        "two": ("electron", "deuterium"),
        "three": ("electron", "deuterium", "carbon"),
    }[case]
    zeff = {"three": {"carbon": 1.0}}.get(case, {})

    L, dz = 1.0, 3.0
    if z_extent is None:
        z_extent = 0.5 * dz
    machine, _ = make_mirror_machine_from_field(
        field=4.0, L=L, ri=0.7855, dz=dz, Nr=1, Nz=10,
        rwall=[0.0, 0.5, 0.5, 0.0], zwall=[-dz * 1.5, -dz * 1.5, dz * 1.5, dz * 1.5],
    )

    eq = Equilibrium(
        tokamak=machine, Rmin=0.05, Rmax=R_OUTER * 1.8,
        Zmin=-z_extent, Zmax=z_extent, nx=nx, ny=ny,
        boundary=freeBoundaryHagenow, psi=np.zeros((nx, ny)), order=4,
    )
    eq._find_critical = False
    eq.tokamak.calcPsiFromGreens(eq._pgreen)
    eq._updatePlasmaPsi(np.zeros((nx, ny)))

    psi_n_knots = np.linspace(0.0, 1.0, N_PSI_SPECIES)
    electron_temp = generate_temperature_profile(psi_n_knots, te_axis, t_edge, 1.0)
    ion_temp = generate_temperature_profile(psi_n_knots, ti_axis, t_edge, 1.0)
    electron_dens = generate_density_z0_profile(
        psi_n_knots, ne_axis, ne_edge, 1.0, width=0.2
    )

    species_list = []
    for name in names:
        kwargs = dict(
            name=name, z=CATALOGUE[name]["z"], mass=CATALOGUE[name]["mass"],
            psi_n=psi_n_knots,
            temp=electron_temp if name == "electron" else ion_temp,
        )
        if name == "electron":
            kwargs["density_Z0"] = electron_dens
            kwargs["flux_density"] = electron_dens
        species_list.append(SingleSpecies(**kwargs))

    all_species = AllSpecies(species_list=species_list)
    all_species.zeff_dict = dict(zeff)

    psi_axis = eq.psiRZ(R_INNER, 0.0)
    psi_bndry = eq.psiRZ(R_OUTER, 0.0)

    # omega( psi_n ), rad/s: a smooth hollow-to-edge ramp. Varying, so that
    # C( psi ) = omega^2 ( Z1 m2 - Z2 m1 )/( Z1 T2 - Z2 T1 ) genuinely drifts --
    # which is the term no published rotating benchmark reaches.
    def omega_profile(psi_n):
        p = np.clip(np.asarray(psi_n, dtype=float), 0.0, 1.0)
        return omega_axis * (1.0 + 0.6 * p - 0.3 * p * p)

    all_species.set_omega(psi_n=psi_n_knots, R=None, omega_profile=omega_profile)

    eq._profiles = ProfilesCentrifugalMirror(all_species=all_species)
    eq.psi_axis = psi_axis
    eq.psi_bndry = psi_bndry

    psi_n_2D = (eq.psi() - psi_axis) / (psi_bndry - psi_axis)
    phi0 = all_species.solve_quasineutrality_global(psi_n=psi_n_2D, R=eq.R, Z=eq.Z)
    all_species.update_phi0(psi_n_2D, eq.R, eq.Z, phi0, interpolator=True)

    return eq, all_species, psi_axis, psi_bndry, names


# ---------------------------------------------------------------------------
# The gauge transfer
# ---------------------------------------------------------------------------


def solve_phi0_pointwise(eq, all_species, names, psi_n_2D, inside):
    """phi_0 on geq's own grid, solved to root-find precision at every node.

    WHY A SECOND REFERENCE EXISTS AT ALL. geq's solve_quasineutrality_global
    converges its phi_0 FIELD to about 5e-05 in Sum_s Z_s n_s, which is ample for
    a Picard equilibrium solve and is NOT ample as a reference for a source that
    solves the same condition to 1e-14 pointwise. Measured, that residual alone
    moves geq's own Jtor by 2.2e-06 (two species) and 6.3e-06 (three) -- so a
    comparison against the as-shipped field bottoms out there and cannot say
    anything sharper, whatever MEQ does.

    NOTHING OF MEQ'S ENTERS HERE. The species model, the profiles and the Jtor
    expression are all geq's; only the scalar condition Sum_s Z_s n_s = 0, which
    both codes state identically, is solved harder. A sign error in this root
    find would show as a LARGER disagreement, not a cancelling one, because what
    consumes phi_0 downstream is geq's own formula.
    """
    species = [all_species.species_dict[n] for n in names]
    omega = all_species.get_omega(np.clip(psi_n_2D, 0.0, 1.0))
    phi0 = np.array(all_species.get_phi0_RZ(eq.R, eq.Z), dtype=float, copy=True)

    for (i, j) in np.argwhere(inside):
        pn = float(np.clip(psi_n_2D[i, j], 0.0, 1.0))
        R = float(eq.R[i, j])
        wsq = float(omega[i, j]) ** 2
        N = [sp.get_flux_density(pn) for sp in species]
        T = [sp.get_temp(pn) * qe for sp in species]

        def residual(phi):
            expo = [(sp.mass * wsq * R * R / 2.0 - sp.z * qe * phi) / t
                    for sp, t in zip(species, T)]
            shift = max(expo)
            return sum(sp.z * n * np.exp(e - shift)
                       for sp, n, e in zip(species, N, expo))

        lo, hi = -1.0, 1.0
        for _ in range(100):
            if residual(lo) * residual(hi) < 0.0:
                break
            lo, hi = lo * 2.0, hi * 2.0
        else:
            raise RuntimeError(f"no bracket for phi_0 at node {(i, j)}")
        phi0[i, j] = brentq(residual, lo, hi, rtol=8.9e-16, maxiter=200)

    return phi0


def species_arrays(all_species, names, psi_n, span):
    """N_s, T_s (J), dlnN_s/dpsi and dT_s/dpsi on a psi_n grid, from geq's splines."""
    out = {}
    for name in names:
        sp = all_species.species_dict[name]
        temp = sp.get_temp(psi_n) * qe
        out[name] = dict(
            mass=sp.mass,
            z=sp.z,
            flux_density=sp.get_flux_density(psi_n),
            temp=temp,
            dlnN_dpsi=sp.get_log_flux_density_derivative(psi_n) / span,
            dT_dpsi=sp._temp_spline.derivative()(psi_n) * qe / span,
        )
    return out


def neutrality_residual(phi, arrays, names, wsq, R_ref):
    """Sum_s Z_s N_s exp( a_s - b_s phi ), with the largest exponent factored out.

    Strictly decreasing in phi -- the derivative is -e Sum_s Z_s^2 N_s e^(..)/T_s
    -- so the root is unique and bracketing cannot pick a spurious one.
    """
    expo = []
    for name in names:
        a = arrays[name]
        expo.append(
            (a["mass"] * wsq * R_ref * R_ref / 2.0 - a["z"] * qe * phi) / a["temp"]
        )
    shift = max(expo)
    total = 0.0
    for name, e in zip(names, expo):
        a = arrays[name]
        total += a["z"] * a["flux_density"] * np.exp(e - shift)
    return total


def transfer_to_reference_radius(arrays, names, omega, domega_dpsi, R_ref):
    """phi_0^geq( R_ref, psi ) and the densities MEQ's gauge wants, with slopes.

    Returns ( n_s0, dn_s0/dpsi ) per species. The derivative is the implicit
    derivative of the quasineutrality condition, not a difference of the values.
    """
    npts = len(omega)
    n0 = {name: np.zeros(npts) for name in names}
    dn0 = {name: np.zeros(npts) for name in names}
    phi = np.zeros(npts)

    for i in range(npts):
        at = {name: {k: (v[i] if isinstance(v, np.ndarray) else v)
                     for k, v in arrays[name].items()} for name in names}
        wsq = omega[i] ** 2

        def f(p):
            return neutrality_residual(p, at, names, wsq, R_ref)

        # Expand a bracket outward from zero. The scale is set by the largest
        # centrifugal energy in the set, which is what phi_0 has to balance.
        scale = max(
            abs(at[n]["mass"] * wsq * R_ref * R_ref / (2.0 * at[n]["z"] * qe))
            for n in names
        )
        lo, hi = -max(scale, 1.0), max(scale, 1.0)
        for _ in range(80):
            if f(lo) * f(hi) < 0.0:
                break
            lo, hi = lo * 2.0, hi * 2.0
        else:
            raise RuntimeError(f"no bracket for phi_0 at index {i}")
        phi[i] = brentq(f, lo, hi, xtol=1.0e-300, rtol=8.9e-16, maxiter=200)

        # Densities at R_ref, then Phi' by implicit differentiation of
        # Sum_s Z_s n_s = 0:  Phi' = Sum Z n ( dlnN + a' - b' Phi ) / Sum Z n b.
        num = 0.0
        den = 0.0
        pieces = {}
        for name in names:
            a = at[name]
            expo = (a["mass"] * wsq * R_ref * R_ref / 2.0
                    - a["z"] * qe * phi[i]) / a["temp"]
            n_s = a["flux_density"] * np.exp(expo)
            n0[name][i] = n_s
            # a_s = m W R_ref^2/( 2 T ),  b_s = Z e/T
            dwsq = 2.0 * omega[i] * domega_dpsi[i]
            a_prime = a["mass"] * R_ref * R_ref / 2.0 * (
                dwsq / a["temp"] - wsq * a["dT_dpsi"] / a["temp"] ** 2
            )
            b = a["z"] * qe / a["temp"]
            b_prime = -a["z"] * qe * a["dT_dpsi"] / a["temp"] ** 2
            pieces[name] = (a_prime, b, b_prime)
            num += a["z"] * n_s * (a["dlnN_dpsi"] + a_prime - b_prime * phi[i])
            den += a["z"] * n_s * b
        phi_prime = num / den

        for name in names:
            a = at[name]
            a_prime, b, b_prime = pieces[name]
            dn0[name][i] = n0[name][i] * (
                a["dlnN_dpsi"] + a_prime - b_prime * phi[i] - b * phi_prime
            )

    return n0, dn0, phi


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_table(path, x, value, slope, header):
    with open(path, "w") as fh:
        for line in header.strip().splitlines():
            fh.write(f"# {line}\n".rstrip() + "\n")
        fh.write("#\n#    psi            f( psi )        f'( psi )\n")
        for a, b, c in zip(x, value, slope):
            fh.write(f"{a:.17e} {b:.17e} {c:.17e}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", choices=("two", "three"), default="two")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nx", type=int, default=65)
    ap.add_argument("--ny", type=int, default=65)
    ap.add_argument("--knots", type=int, default=513,
                    help="knots in the exported profile tables")
    ap.add_argument("--R-ref", type=float, default=None,
                    help="default: the midplane radius of the psi_n = 0.5 surface")
    ap.add_argument("--omega-axis", type=float, default=6.0e5)
    ap.add_argument("--te-axis", type=float, default=1.28e3)
    ap.add_argument("--ti-axis", type=float, default=1.53e3)
    ap.add_argument("--t-edge", type=float, default=3.0e2)
    ap.add_argument("--ne-axis", type=float, default=5.0e19)
    ap.add_argument("--ne-edge", type=float, default=5.0e18)
    ap.add_argument("--z-extent", type=float, default=None,
                    help="half-height of the domain; the default ends it at the "
                         "coil planes, larger keeps the expander inside")
    args = ap.parse_args(argv)

    eq, all_species, psi_ax, psi_bnd, names = build_state(
        args.case, args.nx, args.ny, args.omega_axis,
        args.te_axis, args.ti_axis, args.t_edge, args.ne_axis, args.ne_edge,
        args.z_extent,
    )
    span = psi_bnd - psi_ax

    if args.R_ref is None:
        Z0 = np.argmin(abs(eq.Z[0, :]))
        psi_n_mid = (eq.psi()[:, Z0] - psi_ax) / span
        args.R_ref = float(np.interp(0.5, psi_n_mid, eq.R[:, Z0]))

    # The profile tables, on a psi grid spanning the plasma exactly.
    psi_n = np.linspace(0.0, 1.0, args.knots)
    psi = psi_ax + psi_n * span
    arrays = species_arrays(all_species, names, psi_n, span)
    omega = all_species.get_omega(psi_n)
    domega = all_species.get_domega_dpsin(psi_n) / span

    n0, dn0, phi_ref = transfer_to_reference_radius(
        arrays, names, omega, domega, args.R_ref
    )

    # M^2 = m_i omega^2 R( z=0 )^2/T_e on the psi_n = 0.5 surface -- the number a
    # mirror is specified by, and the one that decides whether the exterior
    # assumption of MIRROR-PLAN.md holds. Reported and recorded, because the
    # gauge transfer's conditioning is a steep function of it.
    Z0mid = np.argmin(abs(eq.Z[0, :]))
    psin_mid = (eq.psi()[:, Z0mid] - psi_ax)/span
    r_half = float(np.interp(0.5, psin_mid, eq.R[:, Z0mid]))
    ion = next(n for n in names if n != "electron")
    mach_half = float(np.sqrt(
        all_species.species_dict[ion].mass*all_species.get_omega(0.5)**2*r_half**2
        /(all_species.species_dict["electron"].get_temp(0.5)*qe)))

    os.makedirs(args.out, exist_ok=True)

    # geq's own current density at its own converged state -- the reference.
    psi2d = eq.psi()
    psi_n_2D = (psi2d - psi_ax) / span
    phi0 = all_species.get_phi0_RZ(eq.R, eq.Z)
    jtor = eq._profiles.Jtor(
        eq.R, eq.Z, psi2d, psi_bndry=psi_bnd, psi_axis=psi_ax, phi0=phi0
    )
    inside = (psi_n_2D >= 0.0) & (psi_n_2D <= 1.0)

    f_ref = MU0 * eq.R * jtor
    sel = inside & np.isfinite(f_ref)

    with open(os.path.join(args.out, "points.dat"), "w") as fh:
        fh.write("# geq's own mu0 R Jtor at its own state -- the reference MEQ's\n")
        fh.write("# F is compared against. Only 0 <= psi_n <= 1 nodes appear:\n")
        fh.write("# geq masks Jtor to zero outside and an unconfined\n")
        fh.write("# meq::RotatingSource does not, so the vacuum region is a\n")
        fh.write("# difference of conventions rather than of physics.\n")
        fh.write("#\n#    R              z               psi             F_ref\n")
        for R, z, p, f in zip(eq.R[sel], eq.Z[sel], psi2d[sel], f_ref[sel]):
            fh.write(f"{R:.17e} {z:.17e} {p:.17e} {f:.17e}\n")

    phi0_tight = solve_phi0_pointwise(eq, all_species, names, psi_n_2D, inside)
    f_tight = MU0 * eq.R * eq._profiles.Jtor(
        eq.R, eq.Z, psi2d, psi_bndry=psi_bnd, psi_axis=psi_ax, phi0=phi0_tight
    )
    with open(os.path.join(args.out, "points_tight.dat"), "w") as fh:
        fh.write("# THE SAME REFERENCE WITH geq's phi_0 SOLVED POINTWISE TO\n")
        fh.write("# ROOT-FIND PRECISION. geq's own field is neutral to about\n")
        fh.write("# 5e-05, which moves its own Jtor by 2.2e-06 and is therefore\n")
        fh.write("# the floor of any comparison against points.dat. Expression,\n")
        fh.write("# profiles and species model are geq's throughout; only the\n")
        fh.write("# scalar condition both codes share is solved harder.\n")
        fh.write("#\n#    R              z               psi             F_ref\n")
        for R, z, p, f in zip(eq.R[sel], eq.Z[sel], psi2d[sel], f_tight[sel]):
            fh.write(f"{R:.17e} {z:.17e} {p:.17e} {f:.17e}\n")

    write_table(
        os.path.join(args.out, "omega.dat"), psi, omega, domega,
        "omega( psi ) in rad/s, sampled from geq's own spline so that both\n"
        "codes read one definition of the rotation.",
    )
    for name in names:
        a = arrays[name]
        write_table(
            os.path.join(args.out, f"T_{name}.dat"), psi, a["temp"], a["dT_dpsi"],
            f"T_{name}( psi ) in JOULES -- geq tabulates eV and meq::Species\n"
            "documents joules, so the conversion is here and not in the TOML.",
        )
        write_table(
            os.path.join(args.out, f"n_{name}.dat"), psi, n0[name], dn0[name],
            f"n_{name}( psi ) in m^-3 ON THE CURVE R = {args.R_ref:.17e}, which\n"
            "is what meq's gauge phi_0( R_ref, psi ) = 0 makes N_s of (96) mean.\n"
            "NOT geq's N_s: see the module docstring for why no single rescaling\n"
            "relates the two.",
        )

        # geq's OWN N_s, in geq's OWN gauge, written so that the comparison can
        # show the transfer is load bearing rather than decorative. Handing THIS
        # to meq is the mistake the module docstring warns about, and a test that
        # cannot tell the two apart is not testing the gauge.
        write_table(
            os.path.join(args.out, f"Nraw_{name}.dat"), psi,
            a["flux_density"], a["dlnN_dpsi"] * a["flux_density"],
            f"N_{name}( psi ) IN geq'S OWN GAUGE -- phi_0 pinned at one global\n"
            "point, psi_n = 0.5 on the midplane. NOT what meq::Species wants.\n"
            "Present so the comparison can assert that using it is WRONG.",
        )

    with open(os.path.join(args.out, "meta.txt"), "w") as fh:
        fh.write(f"case {args.case}\n")
        fh.write(f"mach_at_half {mach_half:.17e}\n")
        fh.write(f"R_ref {args.R_ref:.17e}\n")
        fh.write(f"mu0 {MU0:.17e}\n")
        fh.write(f"psi_axis {psi_ax:.17e}\n")
        fh.write(f"psi_boundary {psi_bnd:.17e}\n")
        fh.write("ggprime 0.0\n")
        fh.write(f"species {len(names)}\n")
        for name in names:
            a = arrays[name]
            fh.write(f"name {name} mass {a['mass']:.17e} charge {a['z']:.17e}"
                     f" temperature T_{name}.dat density n_{name}.dat\n")

    # Diagnostics: what the state is, and whether the transfer is conditioned.
    expo_max = 0.0
    for name in names:
        a = arrays[name]
        expo_max = max(expo_max, float(np.max(np.abs(
            a["mass"] * omega ** 2 * (eq.R.max() ** 2 - args.R_ref ** 2)
            / (2.0 * a["temp"])
        ))))
    resid = np.array([
        sum(arrays[n]["z"] * n0[n][i] for n in names) for i in range(len(psi))
    ])
    scale = np.array([
        sum(abs(arrays[n]["z"]) * n0[n][i] for n in names) for i in range(len(psi))
    ])
    print(f"case                {args.case}   species {names}")
    print(f"M at psi_n = 0.5    {mach_half:.3f}"
          f"   ( m_i omega^2 R(z=0)^2/T_e )")
    print(f"grid                {args.nx} x {args.ny}, {int(sel.sum())} points inside")
    print(f"psi                 [{psi_ax:.6e}, {psi_bnd:.6e}]  span {span:.6e}")
    print(f"R_ref               {args.R_ref:.6f} m")
    print(f"omega               {omega.min():.4e} .. {omega.max():.4e} rad/s"
          f"  drift {omega.max()/omega.min():.3f}x")
    for name in names:
        a = arrays[name]
        print(f"  {name:<10} T {a['temp'].min()/qe:9.2f} .. {a['temp'].max()/qe:9.2f} eV"
              f"   n_s0 {n0[name].min():.4e} .. {n0[name].max():.4e} m^-3"
              f"  ({n0[name].max()/max(n0[name].min(),1e-300):.3e}x)")
    print(f"max |exponent|      {expo_max:.3f}   (over the grid's full radial reach)")
    print(f"neutrality at R_ref {np.max(np.abs(resid)/np.maximum(scale,1e-300)):.3e}"
          "  relative, which is what meq's constructor checks")
    shift = np.linalg.norm(f_tight[sel] - f_ref[sel])/np.linalg.norm(f_ref[sel])
    print(f"|F_ref|             max {np.max(np.abs(f_ref[sel])):.6e}")
    print(f"phi_0 tightening    {shift:.6e} relative in geq's OWN F -- the floor"
          " of any comparison against points.dat")
    print(f"written to          {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
