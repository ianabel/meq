#!/usr/bin/env python
"""EXPORT geq'S VACUUM MIRROR FIELD, FOR MIRROR-PLAN.md MR-3.

No plasma. The field of geq's two-pack mirror machine on geq's own grid, exactly
as geq poses it -- `make_mirror_machine_from_field( field=4.0, L=1.0,
ri=0.7855, dz=3.0, Nr=1, Nz=10 )` on a freegs4e `Equilibrium` with
R in [ 0.05, 0.63 ], Z in [ -1.5, 1.5 ], the recipe of geq's own
tests/conftest.py build_mirror_eq() and of export_mirror.py. It is what
tests/convergence/VacuumMirror.cpp compares MEQ's vacuum solve against.

WHAT geq'S COIL IS, READ FROM freegs4e RATHER THAN ASSUMED, because three
conventions decide whether two codes describe one machine:

  * A PACK IS A STACK OF POINT FILAMENTS. `MultiCoil` with R and Z given as
    lists sets `turns = len( R )` -- the turns argument geq passes is ignored on
    that branch -- and `controlPsi()` sums `Greens( R_fil, Z_fil, R, Z )` over
    the filaments. So the field is `current * sum_fil Greens`, each filament
    carrying `current`, and the pack's total is `current * nfil`, which is
    geq's `total_current_per_pack_A` by its own construction
    ( dI = total/nfil ).
  * THE AREA NEVER ENTERS THE FIELD. It is a current-density limit for the
    self-force and nothing else, as FB-7 records of freegs4e's default `Coil`.
    So geq's pack is ten filaments at R = ri + dr/2 and z = zc +- 0.05, 0.15,
    ..., 0.45 -- with Nr = 1 there is ONE radial filament, at the pack's radial
    centre -- and that is the MIDPOINT RULE meq::filamentStack( coil, 1, 10 )
    builds from the rectangle [ ri, ro ] x [ zc - L/2, zc + L/2 ]. The two
    describe one machine filament for filament.
  * psi IS WEBER PER RADIAN on both sides. freegs4e's Greens carries
    mu0/( 2 pi ) and documents "Webers/2pi"; meq::filamentPsi() is R A_phi.
    Both packs carry the SAME sign of current -- a mirror, not a cusp.

So the reference is a closed-form sum and not a solve: with no plasma geq's
`Equilibrium.psi()` is the machine's Greens sum on the grid, and this script
asserts that rather than trusting it. That is what lets the test attribute
whatever MEQ's disagreement with it turns out to be.

USAGE, from source/

    PYTHONPATH=/home/ian/projects/geq/src \\
    tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_vacuum.py

Written: points.dat (R z psi on every --stride'th node of geq's grid), geq's
ANSWER, to --out, default runs/geq-benchmark/reference/vacuum/, the path the
C++ tests read; and meta.txt (the machine, as geq computed it), THE PROBLEM, to
--tables, default the same directory, reported against the tracked
tools/geq-benchmark/reference/vacuum/meta.txt.
"""

import argparse
import os
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")

from freegs4e.boundary import freeBoundaryHagenow
from freegs4e.equilibrium import Equilibrium
from freegs4e.gradshafranov import Greens

from geq.machine import make_mirror_machine_from_field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from export_mirror import against_tracked, answers_dir  # noqa: E402

MU0 = 4.0e-7 * np.pi

# geq's machine, as tests/conftest.py and export_mirror.py pose it.
FIELD = 4.0
COIL_LENGTH = 1.0
COIL_INNER_RADIUS = 0.7855
COIL_SEPARATION = 3.0
R_OUTER = 0.35


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--nx", type=int, default=65,
                   help="geq's grid, radial nodes (conftest's default)")
    p.add_argument("--ny", type=int, default=65,
                   help="geq's grid, axial nodes (conftest's default)")
    p.add_argument("--stride", type=int, default=2,
                   help="write every stride'th node in each direction; 2 keeps "
                        "the corners and the midplane of a 65 x 65 grid")
    p.add_argument("--out", default=None,
                   help="geq's answer, points.dat; default "
                        "runs/geq-benchmark/reference/vacuum")
    p.add_argument("--tables", default=None,
                   help="the problem, meta.txt; default --out. "
                        "tools/geq-benchmark/reference/vacuum refreshes the "
                        "tracked copy")
    args = p.parse_args()
    args.out = args.out or answers_dir("vacuum")
    args.tables = args.tables or args.out

    machine, info = make_mirror_machine_from_field(
        field=FIELD, L=COIL_LENGTH, ri=COIL_INNER_RADIUS, dz=COIL_SEPARATION,
        Nr=1, Nz=10,
        rwall=[0.0, 0.5, 0.5, 0.0],
        zwall=[-COIL_SEPARATION * 1.5, -COIL_SEPARATION * 1.5,
               COIL_SEPARATION * 1.5, COIL_SEPARATION * 1.5],
    )
    z_extent = 0.5 * COIL_SEPARATION

    eq = Equilibrium(
        tokamak=machine, Rmin=0.05, Rmax=R_OUTER * 1.8,
        Zmin=-z_extent, Zmax=z_extent, nx=args.nx, ny=args.ny,
        boundary=freeBoundaryHagenow, psi=np.zeros((args.nx, args.ny)), order=4,
    )
    eq._find_critical = False
    eq.tokamak.calcPsiFromGreens(eq._pgreen)
    eq._updatePlasmaPsi(np.zeros((args.nx, args.ny)))
    psi = np.asarray(eq.psi(), dtype=float)

    # THE CLAIM THE DOCSTRING MAKES, CHECKED: with no plasma, geq's field is
    # the filament sum and nothing else. Rebuilt here from the coils' own
    # filament lists, so a coil model this script misread would show.
    direct = np.zeros_like(psi)
    filaments = []
    for name, coil in machine.coils:
        for r_fil, z_fil in zip(coil.Rfil, coil.Zfil):
            direct += coil.current * Greens(r_fil, z_fil, eq.R, eq.Z)
            filaments.append((name, float(r_fil), float(z_fil),
                              float(coil.current)))
    agreement = np.max(np.abs(psi - direct)) / np.max(np.abs(direct))
    if agreement > 1e-12:
        print(f"geq's vacuum psi is not its own filament sum: {agreement:.3e}")
        return 1

    total = sum(f[3] for f in filaments) / len(machine.coils)
    if abs(total - info["total_current_per_pack_A"]) > 1e-9 * total:
        print(f"the filaments carry {total:.6e} A per pack, not geq's "
              f"{info['total_current_per_pack_A']:.6e}")
        return 1

    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.tables, exist_ok=True)
    with open(os.path.join(args.tables, "meta.txt"), "w") as fh:
        fh.write("case vacuum\n")
        fh.write(f"mu0 {MU0:.17e}\n")
        fh.write(f"field {FIELD:.17e}\n")
        fh.write(f"coil_inner_radius {info['ri_m']:.17e}\n")
        fh.write(f"coil_outer_radius {info['ro_m']:.17e}\n")
        fh.write(f"coil_length {info['L_m']:.17e}\n")
        fh.write(f"coil_centre_z {0.5 * info['dz_m']:.17e}\n")
        fh.write(f"pack_current {info['total_current_per_pack_A']:.17e}\n")
        fh.write(f"filaments_per_pack {info['nfil_per_pack']}\n")
        fh.write(f"grid {args.nx} {args.ny} stride {args.stride}\n")
        fh.write(f"filament_sum_agreement {agreement:.3e}\n")

    sel_r = np.arange(0, args.nx, args.stride)
    sel_z = np.arange(0, args.ny, args.stride)
    with open(os.path.join(args.out, "points.dat"), "w") as fh:
        fh.write("# geq's vacuum mirror: psi in Wb/rad on every "
                 f"{args.stride}'th node of its {args.nx} x {args.ny} grid.\n")
        fh.write("# No plasma, so this is freegs4e's Greens sum over geq's\n")
        fh.write("# 2 x 10 point filaments -- see export_vacuum.py.\n")
        fh.write("#\n#    R              z               psi\n")
        for i in sel_r:
            for j in sel_z:
                fh.write(f"{eq.R[i, j]:.17e} {eq.Z[i, j]:.17e} "
                         f"{psi[i, j]:.17e}\n")

    print(f"pack current        {info['total_current_per_pack_A']:.6e} A"
          f" in {info['nfil_per_pack']} filaments at R = "
          f"{filaments[0][1]:.6f}")
    print(f"pack                R [{info['ri_m']:.6f}, {info['ro_m']:.6f}]"
          f"  z +-[{0.5 * info['dz_m'] - 0.5 * info['L_m']:.3f}, "
          f"{0.5 * info['dz_m'] + 0.5 * info['L_m']:.3f}]")
    print(f"rho of the copper   {np.hypot(info['ri_m'], 1.0):.4f} .. "
          f"{np.hypot(info['ro_m'], 2.0):.4f}")
    print(f"psi on the grid     {psi.min():.6e} .. {psi.max():.6e}")
    print(f"filament sum        {agreement:.3e} relative to geq's own psi()")
    print(f"written             {len(sel_r) * len(sel_z)} points to {args.out}")
    print(f"problem written to  {args.tables}")
    against_tracked(args.tables, "vacuum")
    return 0


if __name__ == "__main__":
    sys.exit(main())
