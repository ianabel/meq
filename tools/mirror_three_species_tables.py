#!/usr/bin/env python3
"""THE PROFILE TABLES OF examples/mirror-three-species.toml, written once.

Electrons, deuterium and fully stripped carbon on geq's mirror at M = 6, in
MEQ's MIDPLANE gauge (MIRROR-PLAN.md section 5): each density is the density
where its flux surface crosses z = 0. tests/convergence/RotatingMirror.cpp's
`threeSpecies()` poses the same plasma in closed form; this writes it as the
three-column tables a configuration reads.

Every profile is export_rotating.py's skew shape in Psi = 1 - psi_n,

    f( Psi ) = f_edge + ( f_peak - f_edge ) C x^3 ( 1 - x )^4,   x = 1 - Psi,

flat to order 2 at the inner label and 3 at the outer, with its slope in
closed form. The densities fall to 1e-6 of their peak at the labels: under the
midplane gauge a flat density still carries -n m omega^2 R_0 R_0'/T in its
psi-derivative where omega does not vanish, so it is the density that must.

Temperatures are written in eV and converted by the configuration's
TemperatureScale.

USAGE, from source/:    python3 tools/mirror_three_species_tables.py
"""

import math
import os

KNOTS = 257
SKEW_C = 1.0 / ((3.0 / 7.0) ** 3 * (4.0 / 7.0) ** 4)
EV = 1.602176634e-19
DEUTERON_MASS = 3.3435837768e-27
MACH, MACH_RADIUS, OMEGA_DRIFT = 6.0, 0.24, 1.30
TE_PEAK, TI_PEAK = 1.28e3, 1.53e3
DENSITY_PEAK, CARBON_FRACTION, FLOOR = 5.0e18, 0.01, 1.0e-6

OMEGA_PEAK = MACH * math.sqrt(TE_PEAK * EV / DEUTERON_MASS) / MACH_RADIUS

TABLES = {
    "omega": (OMEGA_PEAK, OMEGA_PEAK / OMEGA_DRIFT,
              "omega( Psi ) in rad/s: M = 6 on T_e at R = 0.24 m, drifting 1.30x"),
    "T-electron": (TE_PEAK, 0.1 * TE_PEAK, "T_e( Psi ) in eV"),
    "T-ion": (TI_PEAK, 0.1 * TI_PEAK, "T_D( Psi ) = T_C( Psi ) in eV"),
    "n-electron": (DENSITY_PEAK, FLOOR * DENSITY_PEAK,
                   "n_e( Psi ) in m^-3 WHERE EACH SURFACE CROSSES z = 0"),
    "n-carbon": (CARBON_FRACTION * DENSITY_PEAK, FLOOR * CARBON_FRACTION * DENSITY_PEAK,
                 "n_C( Psi ) in m^-3 WHERE EACH SURFACE CROSSES z = 0"),
}


def write(path, peak, edge, what):
    with open(path, "w") as fh:
        fh.write(f"# {what}. Written by tools/mirror_three_species_tables.py:\n")
        fh.write("# f_edge + ( f_peak - f_edge ) C x^3 ( 1 - x )^4, x = 1 - Psi.\n")
        fh.write(f"# f_peak {peak:.17e}, f_edge {edge:.17e}. Columns: Psi f df/dPsi.\n")
        for i in range(KNOTS):
            psi = i / (KNOTS - 1.0)
            x = 1.0 - psi
            value = edge + (peak - edge) * SKEW_C * x ** 3 * (1.0 - x) ** 4
            slope = -(peak - edge) * SKEW_C * x ** 2 * (1.0 - x) ** 3 * (3.0 - 7.0 * x)
            fh.write(f"{psi:.17e} {value:.17e} {slope:.17e}\n")


def main():
    examples = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "examples")
    for stem, (peak, edge, what) in TABLES.items():
        write(os.path.join(examples, f"mirror-three-species-{stem}.dat"), peak, edge, what)
    print(f"omega_peak {OMEGA_PEAK:.6e} rad/s; tables written to {examples}")


if __name__ == "__main__":
    main()
