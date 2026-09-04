#!/usr/bin/env python3
"""Plot a MEQ equilibrium from its NetCDF output.

Reads the ( R, Z ) grid file MEQ writes -- psi, B_R, B_Z and the `inside`
mask -- and draws flux surfaces, field magnitude, or both.

    tools/plot_equilibrium.py run.nc
    tools/plot_equilibrium.py run.nc --what surfaces --levels 30 -o psi.png
    tools/plot_equilibrium.py run.nc --what both --show

WHY THE MASK IS USED RATHER THAN THE FILL VALUE.  MEQ writes NaN outside the
domain *and* a zero in `inside`, deliberately, because some readers honour the
fill attribute and some do not.  This script trusts `inside`: a NaN that
arrived some other way -- a solve that diverged, say -- is then still visible
as a gap rather than being silently indistinguishable from "outside".

THE FILE ALSO CARRIES AN `extrapolated` MASK, which this script does not use
and which matters for anything quantitative.  A node in the band between
Gamma_h and the true boundary is inside the plasma and carries a real value --
so `inside` is 1 there, correctly -- but that value was continued outward from
the mesh boundary rather than solved on, and is an order less accurate than its
neighbours.  For a picture that is invisible.  For an error norm, a fitted flux
surface, or a difference between two resolutions, drop `extrapolated != 0`.
"""

import argparse
import sys

import numpy as np


def read(path):
    """Pull the fields out of a MEQ NetCDF file, masked to the domain."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        sys.exit("plot_equilibrium: needs netCDF4 (pip install netCDF4)")

    with Dataset(path) as ds:
        missing = [n for n in ("R", "Z", "psi", "B_R", "B_Z")
                   if n not in ds.variables]
        if missing:
            sys.exit(f"plot_equilibrium: {path} has no {', '.join(missing)} -- "
                     "is it a MEQ output file?")

        data = {
            "R": np.asarray(ds.variables["R"][:], dtype=float),
            "Z": np.asarray(ds.variables["Z"][:], dtype=float),
            "psi": np.asarray(ds.variables["psi"][:], dtype=float),
            "B_R": np.asarray(ds.variables["B_R"][:], dtype=float),
            "B_Z": np.asarray(ds.variables["B_Z"][:], dtype=float),
        }
        # `inside` is optional only so that a file written by something other
        # than MEQ still plots; MEQ always writes it.
        if "inside" in ds.variables:
            data["inside"] = np.asarray(ds.variables["inside"][:]) != 0
        else:
            data["inside"] = np.isfinite(data["psi"])

        if "boundary_R" in ds.variables and "boundary_Z" in ds.variables:
            data["boundary"] = (
                np.asarray(ds.variables["boundary_R"][:], dtype=float),
                np.asarray(ds.variables["boundary_Z"][:], dtype=float),
            )
        else:
            data["boundary"] = None

        data["attributes"] = {k: ds.getncattr(k) for k in ds.ncattrs()}

    outside = ~data["inside"]
    for key in ("psi", "B_R", "B_Z"):
        data[key] = np.where(outside, np.nan, data[key])
    return data


def subtitle(attributes):
    """One line of provenance, from whichever attributes the file carries."""
    bits = []
    for key, label in (("polynomial_degree", "k"),
                       ("elements", "elements"),
                       ("newton_iterations", "Newton"),
                       ("adaptive_cycles", "cycles")):
        if key in attributes:
            bits.append(f"{label} {attributes[key]}")
    if "boundary" in attributes:
        bits.append(str(attributes["boundary"]).split(" (")[0])
    # Which potential the psi variable holds. Since 2026-09-02 it is psi*, the
    # post-processed field, one degree richer than the solve; before that it was
    # psi_h and there is no attribute. So the ABSENCE of this is what identifies
    # an older file, and two runs differenced across that change measure the
    # post-processing rather than the physics -- worth having on the picture.
    if "potential" in attributes:
        bits.append(str(attributes["potential"]).split(" (")[0])
    if "final_residual" in attributes:
        bits.append(f"|r| {float(attributes['final_residual']):.1e}")
    return "   ".join(bits)


def draw(ax, data, what, levels):
    R, Z, psi = data["R"], data["Z"], data["psi"]
    magnitude = np.hypot(data["B_R"], data["B_Z"])

    if what in ("field", "both"):
        mesh = ax.pcolormesh(R, Z, magnitude, shading="auto")
        bar = ax.figure.colorbar(mesh, ax=ax)
        bar.set_label(r"$|B_{pol}|$  [T]")

    if what in ("surfaces", "both"):
        # Flux surfaces ARE the contours of psi, so a contour plot is the
        # physics here rather than a way of colouring it in.
        #
        # linestyles is forced solid because matplotlib dashes negative
        # contours by default. That convention encodes the sign of psi, which
        # sounds useful and is not: with a datum that changes sign most of the
        # surfaces come out dashed and the plot just looks noisy. The sign is
        # on the colourbar instead, where it can be read.
        if what == "both":
            ax.contour(R, Z, psi, levels=levels, colors="white",
                       linewidths=0.6, linestyles="solid")
        else:
            lines = ax.contour(R, Z, psi, levels=levels, linewidths=0.8,
                               linestyles="solid")
            bar = ax.figure.colorbar(lines, ax=ax)
            bar.set_label(r"$\psi$  [Wb/rad]")

        # psi = 0 heavier: on the fixed-boundary problem that level set is the
        # plasma boundary by definition, so it is the one contour that means
        # something different from its neighbours.
        if np.nanmin(psi) < 0.0 < np.nanmax(psi):
            ax.contour(R, Z, psi, levels=[0.0],
                       colors="white" if what == "both" else "k",
                       linewidths=1.6, linestyles="solid")

    if data["boundary"] is not None:
        bR, bZ = data["boundary"]
        ax.plot(np.append(bR, bR[:1]), np.append(bZ, bZ[:1]),
                color="C3", linewidth=1.2, label=r"$\Gamma$")
        ax.legend(loc="upper right", fontsize="small")

    ax.set_xlabel("R  [m]")
    ax.set_ylabel("Z  [m]")
    # Equal aspect is not cosmetic: a poloidal cross-section drawn with
    # unequal axes misrepresents elongation and triangularity, which are
    # usually the first things anybody looks at.
    ax.set_aspect("equal")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot a MEQ equilibrium from its NetCDF output.")
    parser.add_argument("file", help="the .nc file MEQ wrote")
    parser.add_argument("--what", default="both",
                        choices=("surfaces", "field", "both"),
                        help="flux surfaces, |B_pol|, or both (default: both)")
    parser.add_argument("--levels", type=int, default=25,
                        help="number of psi contours (default: 25)")
    parser.add_argument("-o", "--output", default=None,
                        help="write a figure here instead of showing it")
    parser.add_argument("--show", action="store_true",
                        help="open a window even when -o is given")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args(argv)

    data = read(args.file)

    import matplotlib
    # Choose the backend BEFORE pyplot is imported, and only when there is
    # nothing to show: on a headless machine the default backend fails at
    # import rather than at draw time, which reads as a broken script.
    if args.output and not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, ax = plt.subplots(figsize=(5.5, 7.0))
    draw(ax, data, args.what, args.levels)

    title = data["attributes"].get("title", "MEQ equilibrium")
    line = subtitle(data["attributes"])
    ax.set_title(f"{title}\n{line}" if line else title, fontsize="medium")
    figure.tight_layout()

    if args.output:
        figure.savefig(args.output, dpi=args.dpi)
        print(f"wrote {args.output}")
    if args.show or not args.output:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
