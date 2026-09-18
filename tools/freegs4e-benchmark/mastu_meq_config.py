#!/usr/bin/env python
"""THE MEQ CONFIGURATION FOR freegsnke'S MAST-U FORWARD CASE.

Emits `[[coils]]`, the two profile tables and a TOML for the case
`mastu_reference.py` solves, so `race_nke.py` has a MEQ arm.
`examples/machine-g-mastu.toml` is a DIFFERENT MAST-U -- it comes from the
freegs4e benchmark's G case, with different conductors and a different limiter
-- and must not be substituted.

THREE CONVENTIONS, EACH SETTLED BY MEASUREMENT RATHER THAN BY READING, because
this tree records four separate occasions when a sign or a factor here was
settled by argument and was wrong.

  AMPERE-TURNS.  Every filament of a freegsnke MultiCoil carries the FULL coil
  current, not a share of it.  Checked against freegs4e's own Greens(): summing
  the Green's function over D1's 70 filaments at 1 A each reproduces
  tokamak.psi() to a ratio of 1.000000, where current/nfil is out by exactly 70.
  Confirmed independently against the worked G precedent -- G's Solenoid is
  -3445.589 A in the npz and -1.1163709182e+06 in machine-g-mastu.toml, a ratio
  of exactly 324.0, its filament count.

  THE SIGN IS UNCHANGED.  Same precedent, same direction: negative stays
  negative.

  THE PROFILE ABSCISSA FLIPS AND THE VALUES SCALE BY ( psi_ax - psi_bnd ).
  MEQ's normalised flux is 1 on the axis where freegs4e's is 0, so
  `Psi = 1 - psi_n`; and meq::NormalisedMHDSource evaluates the tables at Psi_N
  while freegs4e's pprime is d/dpsi against UNNORMALISED flux, so the chain rule
  puts ( psi_ax - psi_bnd ) on the value column.  Measured against
  machine-g-mastu-pprime.dat at three points: the ratio to the npz's own pprime
  is 0.083339, 0.083338, 0.083317 against that case's psi_ax - psi_bnd of
  8.3338e-02.

A CAVEAT THAT IS THE PROFILE'S AND NOT THE CONVERSION'S.  `alpha_n = 1.2` makes
the edge a FRACTIONAL power, `( 1 - Psi_N^1.8 )^1.2`.  convert.py's own
self-test measures what that costs a tabulated derivative -- 5.6e-02 relative at
the last node against 2.6e-04 in the interior -- and CLAUDE_FB.md's plasma-edge
result caps the achievable order at `k <= j` with `j = 1.2` here, i.e. `k = 1`.
So this case cannot be expected to show MEQ's `k = 2` or `k = 3` rates, and a
disappointing order on it is the PROFILE rather than the solver.
"""
import argparse, os, pickle, sys
import numpy as np

M = "/home/ian/projects/freegsnke/machine_configs/MAST-U"
CURRENTS = "/home/ian/projects/freegsnke/examples/data/simple_diverted_currents_PaxisIp.pk"
MU0 = 1.25663706212e-6


def subgroups(name, entry):
    """A coil's filament groups: ( label, R, Z, dR, dZ, multiplier )."""
    if "R" in entry:
        return [(name, np.asarray(entry["R"], float), np.asarray(entry["Z"], float),
                 float(entry["dR"]), float(entry["dZ"]),
                 float(entry.get("polarity", 1)) * float(entry.get("multiplier", 1)))]
    out = []
    for key, sub in entry.items():
        out.append((f"{name}{key}", np.asarray(sub["R"], float),
                    np.asarray(sub["Z"], float), float(sub["dR"]), float(sub["dZ"]),
                    float(sub.get("polarity", 1)) * float(sub.get("multiplier", 1))))
    return out


def solve_reference(nx, ny):
    """The converged forward solve, for L, Beta0, Raxis and the normalisation."""
    from freegsnke import build_machine, equilibrium_update, GSstaticsolver
    from freegsnke.jtor_update import ConstrainPaxisIp
    tok = build_machine.tokamak(
        active_coils_path=f"{M}/MAST-U_like_active_coils.pickle",
        passive_coils_path=f"{M}/MAST-U_like_passive_coils.pickle",
        limiter_path=f"{M}/MAST-U_like_limiter.pickle",
        wall_path=f"{M}/MAST-U_like_wall.pickle")
    eq = equilibrium_update.Equilibrium(tokamak=tok, Rmin=0.1, Rmax=2.0,
                                        Zmin=-2.2, Zmax=2.2, nx=nx, ny=ny)
    profiles = ConstrainPaxisIp(eq=eq, paxis=8e3, Ip=6e5, fvac=0.5,
                                alpha_m=1.8, alpha_n=1.2)
    solver = GSstaticsolver.NKGSsolver(eq, gs_operator_order=4)
    with open(CURRENTS, "rb") as handle:
        for label, value in pickle.load(handle).items():
            eq.tokamak.set_coil_current(coil_label=label, current_value=value)
    solver.solve(eq=eq, profiles=profiles, constrain=None,
                 target_relative_tolerance=1e-9, verbose=False)
    return eq, profiles


def write_table(path, psi, value, derivative, title):
    with open(path, "w") as f:
        f.write(f"# {title}\n#\n")
        f.write("# Psi = 1 - psi_n: freegs4e's normalised flux is 0 on the AXIS\n"
                "# and 1 on the boundary, and MEQ's is the other way up.\n#\n")
        f.write("# The value column is d/dPsi and carries ( psi_ax - psi_bnd )\n"
                "# from the chain rule; see this generator's header for the\n"
                "# measurement that settled it.\n#\n")
        f.write(f"# {'Psi':>10} {'d/dPsi':>26} {'d2/dPsi2':>26}\n")
        for a, b, c in zip(psi, value, derivative):
            f.write(f"  {a:+10.6f}   {b:+.15e}   {c:+.15e}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/home/ian/projects/meq/examples")
    ap.add_argument("--stem", default="mastu-nke")
    ap.add_argument("--nx", type=int, default=65)
    ap.add_argument("--ny", type=int, default=129)
    ap.add_argument("--degree", type=int, default=2)
    args = ap.parse_args()

    eq, profiles = solve_reference(args.nx, args.ny)
    psi_ax, psi_bnd = float(eq.psi_axis), float(eq.psi_bndry)
    span = psi_ax - psi_bnd
    print(f"reference: psi_axis {psi_ax:.9e}  psi_bndry {psi_bnd:.9e}  "
          f"span {span:.9e}  Ip {eq.plasmaCurrent():.6e}")
    print(f"           L {profiles.L:.9e}  Beta0 {profiles.Beta0:.9e}  "
          f"Raxis {profiles.Raxis:.6f}")

    # ---- the profile tables ----
    # Psi from -0.3 to 1.0, matching the G case's range so the two are read the
    # same way.  Below Psi = 0 the shape clips to zero, which is the plasma edge
    # and is why ConfineToPlasma can be trusted here.
    Psi = np.round(np.arange(-0.30, 1.0 + 1e-9, 0.01), 6)
    pn = 1.0 - Psi
    for label, values in (("pprime", profiles.pprime(pn) * span),
                          ("ggprime", profiles.ffprime(pn) * span)):
        # d/dPsi by central differences of the tabulated column, which is what
        # MEQ's spline will differentiate anyway.
        derivative = np.gradient(values, Psi)
        path = os.path.join(args.out, f"{args.stem}-{label}.dat")
        write_table(path, Psi, values, derivative,
                    f"d{'p' if label=='pprime' else 'gg'}/dPsi( Psi ) for "
                    f"freegsnke MAST-U forward (example02)")
        print(f"wrote {path}")

    # ---- the coils ----
    active = pickle.load(open(f"{M}/MAST-U_like_active_coils.pickle", "rb"))
    currents = pickle.load(open(CURRENTS, "rb"))
    blocks, extent_r, extent_z = [], [], []
    for name, entry in active.items():
        base = float(currents.get(name, 0.0))
        for label, R, Z, dR, dZ, mult in subgroups(name, entry):
            lo_r, hi_r = R.min() - dR / 2, R.max() + dR / 2
            lo_z, hi_z = Z.min() - dZ / 2, Z.max() + dZ / 2
            extent_r += [lo_r, hi_r]; extent_z += [lo_z, hi_z]
            total = base * len(R) * mult          # AMPERE-TURNS, sign unchanged
            blocks.append(f"""[[coils]]
Name = "{label}"
CentreR = {0.5*(lo_r+hi_r):.9f}
CentreZ = {0.5*(lo_z+hi_z):.9f}
HalfWidth = {0.5*(hi_r-lo_r):.9f}
HalfHeight = {0.5*(hi_z-lo_z):.9f}
Current = {total:+.10e}
""")
    radius = float(np.hypot(max(extent_r), max(abs(min(extent_z)), max(extent_z))))
    print(f"conductors span R [{min(extent_r):.4f}, {max(extent_r):.4f}] "
          f"Z [{min(extent_z):.4f}, {max(extent_z):.4f}]; "
          f"they need a disc of radius {radius:.3f}")

    # ---- the null the psi_bnd border follows ----
    # A GENUINE DOUBLE NULL: the reference's two primary saddles sit at
    # ( 0.5985, -1.0972 ) and ( 0.5985, +1.0972 ) at IDENTICAL flux, and that
    # flux IS psi_bndry. MEQ's [boundary.xpoint] follows ONE saddle chosen by
    # its seed, so the lower is seeded arbitrarily -- the two are equivalent in
    # psi_bnd, which is what makes an arbitrary choice safe here and would not
    # in an up-down asymmetric machine.
    saddles = np.asarray(eq.xpt, float)
    lower = saddles[np.argmin(saddles[:, 1])] if len(saddles) else None
    primary = min((row for row in saddles),
                  key=lambda row: abs(row[2] - psi_bnd))
    print(f"nulls: {len(saddles)}; seeding ( {primary[0]:.6f}, {primary[1]:.6f} ) "
          f"at psi {primary[2]:.9e} against psi_bndry {psi_bnd:.9e}")

    limiter = pickle.load(open(f"{M}/MAST-U_like_limiter.pickle", "rb"))
    lim = np.asarray( [ [ p["R"], p["Z"] ] for p in limiter ], float )
    print(f"limiter: {len(lim)} points, R [{lim[:,0].min():.4f}, {lim[:,0].max():.4f}] "
          f"Z [{lim[:,1].min():.4f}, {lim[:,1].max():.4f}]")
    np.savetxt(os.path.join(args.out, f"{args.stem}-limiter.dat"), lim,
               header="R Z -- freegsnke MAST-U limiter, 47 points")

    disc = float(np.ceil((radius + 0.6) * 10) / 10)
    toml = f"""# freegsnke's MAST-U FORWARD case (examples/example02), as MEQ solves it.
#
# GENERATED by tools/freegs4e-benchmark/mastu_meq_config.py -- edit that and
# regenerate rather than editing this, because the coil currents are
# AMPERE-TURNS derived from a filament count and the profile tables carry a
# ( psi_ax - psi_bnd ) from the chain rule. Both are easy to "correct" wrongly.
#
# The reference is mastu-fwd-<nx>x<ny>.npz from mastu_reference.py. Its own
# accuracy ceiling is about 1e-05 in psi_axis (MEASUREMENTS.md M-122), so
# agreement below that is the reference's noise.
#
# alpha_n = 1.2 makes the plasma edge a FRACTIONAL power, which caps the
# achievable order at k <= 1.2 -- see the generator's header.

[mesh]
File = "examples/{args.stem}.msh"

[mesh.generate]
Tool = "halfdisc"
Radius = {disc:.1f}
Size = 0.30
PlasmaRMin = {max(0.05, lim[:,0].min()):.4f}
PlasmaRMax = {lim[:,0].max():.4f}
PlasmaZMin = {lim[:,1].min():.4f}
PlasmaZMax = {lim[:,1].max():.4f}
PlasmaSize = 0.06
CoilSize = 0.03

[boundary]
Type = "zero"

# THE TARGET NULL as an INITIAL VALUE of two unknowns, not a prescription:
# XP-3's border closes q_r = q_z = 0 and psi_bnd = psi_h( r_X, z_X ) on the
# same Newton, and the seed selects WHICH saddle is followed. This case is a
# DOUBLE NULL at identical flux, so either seed gives the same psi_bnd.
[boundary.xpoint]
R = {primary[0]:.12f}
Z = {primary[1]:.12f}

# Strictly inside [mesh.generate] Radius, which is what Gamma being the smaller
# semicircle within the generated disc means.
[boundary.exterior]
Radius = {disc - 0.4:.1f}
CentreZ = 0.0
Modes = 10

[discretisation]
PolynomialDegree = {args.degree}

[solver]
NewtonMaxIterations = 200
NewtonRelativeTolerance = 1.0e-10
PlasmaSupportSweeps = 4

[source]
Type = "mhd"
Normalised = true
ConfineToPlasma = true
PsiAxis = {psi_ax:.9e}
PlasmaCurrent = {float(eq.plasmaCurrent()):.10e}
Mu0 = {MU0}
PPrimeFile = "examples/{args.stem}-pprime.dat"
GGPrimeFile = "examples/{args.stem}-ggprime.dat"

{''.join(blocks)}"""
    path = os.path.join(args.out, f"{args.stem}.toml")
    with open(path, "w") as handle:
        handle.write(toml)
    print(f"wrote {path} with {len(blocks)} coil blocks, disc radius {disc:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
