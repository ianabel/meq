#!/usr/bin/env python3
"""Plot a MEQ equilibrium from its NetCDF output.

Reads either of the two NetCDF files a run writes, and tells them apart by
their own dimensions rather than by the name:

  * the ( R, Z ) grid file -- psi, B_R, B_Z and the `inside` mask -- drawn as
    flux surfaces, field magnitude, or both;
  * the ( Psi, theta ) flux-surface file `<stem>_surfaces.nc`, written when
    `[output] FluxSurfaces` is set -- the traced surfaces in the poloidal
    plane, and the flux-surface averages against the label rho = sqrt( Psi_N ).

    tools/plot_equilibrium.py run.nc
    tools/plot_equilibrium.py run.nc --what surfaces --levels 30 -o psi.png
    tools/plot_equilibrium.py run_surfaces.nc -o family.png
    tools/plot_equilibrium.py run_surfaces.nc --what profiles -o averages.png
    tools/plot_equilibrium.py run.nc --what both --show

WHY THE MASK IS USED RATHER THAN THE FILL VALUE.  MEQ writes NaN outside the
domain *and* a zero in `inside`, deliberately, because some readers honour the
fill attribute and some do not.  This script trusts `inside`: a NaN that
arrived some other way -- a solve that diverged, say -- is then still visible
as a gap rather than being silently indistinguishable from "outside".

THE GRID FILE ALSO CARRIES AN `extrapolated` MASK, which this script does not
use *there* and which matters for anything quantitative.  A node in the band
between Gamma_h and the true boundary is inside the plasma and carries a real
value -- so `inside` is 1 there, correctly -- but that value was continued
outward from the mesh boundary rather than solved on, and is an order less
accurate than its neighbours.  For a picture that is invisible.  For an error
norm, a fitted flux surface, or a difference between two resolutions, drop
`extrapolated != 0`.

ON THE SURFACE FILE THE BAND MASK IS DRAWN, AND HAS TO BE.  A surface can be
inside Omega_h at one theta and outside it at the next, so which arc of a
closed curve is solved data cannot be inferred from the curve: the extension
answers as confidently as an element does.  So band nodes are drawn per node,
in red, over the surface they belong to -- markers as well as a line, because
the first surface to cross the band does so at a handful of isolated nodes and
those make no line segment at all -- and the per-surface band fraction gets a
panel of its own beside the averages.  The profile panels ring the surfaces
`band( flux )` marks, which is the right flag *there*: a flux-surface average
over a surface with one band node is contaminated as a whole.

THE FAMILY IS CUT AT BOTH ENDS AND THE PLOT SAYS SO.  The profile axes span
the whole of rho in [ 0, 1 ] with the two cut regions shaded, so the absence of
an answer at the magnetic axis and at the plasma boundary is visible rather
than being papered over by an axis that starts where the data does.  MEQ
refuses a query outside the cut; it does not extrapolate, and neither does
this.
"""

import argparse
import re
import sys

import numpy as np


# Variables on the flux dimension that are a coordinate or a per-surface
# diagnostic rather than a quantity to plot against rho.  Everything else the
# file carries on that dimension gets a panel, so a family carrying
# `safety_factor` -- which needs a g( psi ) the driver has none of, a
# meq::Source carrying g g' and not g -- plots one without this script having
# to be told the variable exists.
FLUX_COORDINATES = ("rho", "normalised_flux", "psi")
FLUX_DIAGNOSTICS = ("worst_residual", "transversality", "band")


def read(path):
    """Pull a MEQ NetCDF file apart, whichever of the two kinds it is."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        sys.exit("plot_equilibrium: needs netCDF4 (pip install netCDF4)")

    with Dataset(path) as ds:
        # Told apart by the file's own dimensions and not by its name, so a
        # renamed copy, or one a pipeline moved, still reads.
        if "flux" in ds.dimensions and "theta" in ds.dimensions:
            return read_surfaces(path, ds)
        return read_grid(path, ds)


def read_grid(path, ds):
    """The ( R, Z ) grid file, masked to the domain."""
    missing = [n for n in ("R", "Z", "psi", "B_R", "B_Z")
               if n not in ds.variables]
    if missing:
        sys.exit(f"plot_equilibrium: {path} has no {', '.join(missing)} -- "
                 "is it a MEQ output file?")

    data = {
        "kind": "grid",
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


def read_surfaces(path, ds):
    """The ( Psi, theta ) flux-surface file: the surfaces and the averages."""
    missing = [n for n in ("rho", "normalised_flux", "theta", "R", "Z")
               if n not in ds.variables]
    if missing:
        sys.exit(f"plot_equilibrium: {path} has a flux dimension but no "
                 f"{', '.join(missing)} -- is it a MEQ flux-surface file?")

    def column(name):
        return np.asarray(ds.variables[name][:], dtype=float)

    def described(name):
        var = ds.variables[name]
        attributes = var.ncattrs()
        return {
            "values": column(name),
            "long_name": (var.getncattr("long_name")
                          if "long_name" in attributes else name),
            "units": var.getncattr("units") if "units" in attributes else "",
        }

    data = {
        "kind": "surfaces",
        "rho": column("rho"),
        "normalised_flux": column("normalised_flux"),
        "theta": column("theta"),
        "R": column("R"),
        "Z": column("Z"),
        "psi": column("psi") if "psi" in ds.variables else None,
    }

    # THE BAND, PER NODE.  A count is not a mask: the `extrapolated_nodes`
    # attribute says how many nodes are continued and nothing whatever about
    # which, and the per-surface `band` flag is a summary that under-reports a
    # partial excursion.  Both are read, because each is right somewhere, and
    # neither substitutes for the per-node mask.
    if "extrapolated" in ds.variables:
        data["extrapolated"] = np.asarray(ds.variables["extrapolated"][:]) != 0
    else:
        data["extrapolated"] = np.zeros(data["R"].shape, dtype=bool)

    if "band" in ds.variables:
        data["band"] = np.asarray(ds.variables["band"][:]) != 0
    else:
        data["band"] = data["extrapolated"].any(axis=1)

    # Whatever the file carries on the flux dimension, in the file's own
    # order, which is the writer's.  Deciding here which averages exist would
    # be a second copy of that list, and a stale one the moment the family
    # gains a column.
    data["profiles"] = {
        name: described(name) for name, var in ds.variables.items()
        if var.dimensions == ("flux",)
        and name not in FLUX_COORDINATES and name not in FLUX_DIAGNOSTICS
    }
    data["diagnostics"] = {
        name: described(name) for name in FLUX_DIAGNOSTICS
        if name in ds.variables and ds.variables[name].dimensions == ("flux",)
        and name != "band"
    }

    data["attributes"] = {k: ds.getncattr(k) for k in ds.ncattrs()}
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


def wrapped(bits, budget):
    """Join provenance fragments, breaking between them to fit a width.

    The estimate is crude on purpose: a mathtext command renders as one glyph
    and is several characters of source, so the length is taken after the
    markup is stripped.  Being a little conservative costs a line break and
    being wrong the other way runs the caption off the edge of the figure,
    which is what this exists to stop.
    """
    def shown(text):
        return len(re.sub(r"\\[a-zA-Z]+|[${}]", "", text))

    lines, line = [], ""
    for bit in bits:
        if line and shown(line) + 3 + shown(bit) > budget:
            lines.append(line)
            line = bit
        else:
            line = f"{line}   {bit}" if line else bit
    if line:
        lines.append(line)
    return "\n".join(lines)


def surface_subtitle(data):
    """Provenance for the flux-surface file, as fragments to be wrapped."""
    attributes = data["attributes"]
    surfaces, angles = data["R"].shape
    inner = float(np.min(data["normalised_flux"]))
    outer = float(np.max(data["normalised_flux"]))
    band = int(np.count_nonzero(data["extrapolated"]))

    bits = [f"{surfaces} surfaces $\\times$ {angles} nodes",
            f"$\\Psi_N \\in [{inner:.3f}, {outer:.3f}]$"]
    # The levels themselves, in MEQ's own units, so that the label on the axes
    # and the flux the solver reports can be lined up without arithmetic.
    if data["psi"] is not None:
        bits.append(f"$\\psi \\in [{float(np.min(data['psi'])):.4f}, "
                    f"{float(np.max(data['psi'])):.4f}]$ Wb/rad")
    bits.append(f"band {band}/{surfaces*angles} nodes")
    for key, label in (("polynomial_degree", "k"),
                       ("source_type", "source")):
        if key in attributes:
            bits.append(f"{label} {attributes[key]}")
    if "potential" in attributes:
        bits.append(str(attributes["potential"]).split(" (")[0])
    # The worst | psi_h - level | over the whole family.  FluxFamily.hpp calls
    # it "the thing to read before believing any of the rest", so it belongs
    # on the picture rather than only in the file.
    if "worst_residual" in attributes:
        bits.append(r"worst $|\psi_h - c|$ "
                    f"{float(attributes['worst_residual']):.1e}")
    return bits


def profile_panels(data, requested):
    """What to draw against rho: the file's own columns, plus the band.

    The default is every quantity the file carries on the flux dimension, in
    the order it carries them, which is how `safety_factor` appears when a
    caller supplied a g( psi ) and stays absent -- rather than zero -- when
    nobody did.  The per-surface diagnostics are left out unless asked for:
    `worst_residual` is already on the figure as a single number and
    `transversality` is a property of the fit rather than of the equilibrium.
    """
    available = dict(data["profiles"])
    available.update(data["diagnostics"])

    if requested is None:
        names = list(data["profiles"])
    elif requested == "all":
        names = list(available)
    else:
        names = [n.strip() for n in requested.split(",") if n.strip()]
        unknown = [n for n in names if n not in available]
        if unknown:
            sys.exit(f"plot_equilibrium: no such profile: {', '.join(unknown)}"
                     f" -- the file carries {', '.join(available)}")

    panels = [{"name": name,
               "values": available[name]["values"],
               "units": available[name]["units"]} for name in names]

    # THE BAND FRACTION IS DERIVED AND IS THE POINT OF THE PANEL.  Every other
    # panel is a number the file states; this one is the per-node mask summed
    # per surface, and it is what says how much of the curve beside it is
    # solved data.  A surface flagged by `band` may be one node deep into the
    # extension or nine tenths of the way through it, and the two are not the
    # same statement about V'( rho ).
    angles = data["extrapolated"].shape[1]
    panels.append({
        "name": "band fraction",
        "values": data["extrapolated"].sum(axis=1)/float(angles),
        "units": "1",
        "derived": True,
    })
    return panels


def draw_surface_geometry(ax, data):
    """The traced surfaces in the poloidal plane, with the band drawn on."""
    import matplotlib as mpl

    rho = data["rho"]
    R, Z = data["R"], data["Z"]
    band = data["extrapolated"]

    norm = mpl.colors.Normalize(vmin=float(rho.min()), vmax=float(rho.max()))
    cmap = mpl.colormaps["viridis"]

    for i in range(len(rho)):
        # theta runs over 2 pi j / N, so the loop is open in the file and a
        # reader closes it -- the same convention `boundary_R/Z` uses on the
        # grid file.
        r = np.append(R[i], R[i, 0])
        z = np.append(Z[i], Z[i, 0])
        ax.plot(r, z, color=cmap(norm(rho[i])), linewidth=0.9)

        mask = np.append(band[i], band[i, 0])
        if mask.any():
            # Markers as well as a line: the innermost surface to cross the
            # band does so at a few isolated nodes, and an isolated node makes
            # no segment for a line to draw.
            ax.plot(np.where(mask, r, np.nan), np.where(mask, z, np.nan),
                    color="C3", linewidth=1.5, linestyle="solid",
                    marker=".", markersize=2.0)

    if band.any():
        ax.plot([], [], color="C3", linewidth=1.5, marker=".", markersize=2.0,
                label=r"band data, continued $\Gamma_h \to \Gamma$")

    attributes = data["attributes"]
    if "axis_r" in attributes and "axis_z" in attributes:
        ax.plot(float(attributes["axis_r"]), float(attributes["axis_z"]),
                marker="+", color="k", markersize=9, markeredgewidth=1.4,
                linestyle="none", label="magnetic axis")

    bar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                             ax=ax, fraction=0.055, pad=0.03)
    bar.set_label(r"$\rho = \sqrt{\Psi_N}$")

    inner = float(np.min(data["normalised_flux"]))
    outer = float(np.max(data["normalised_flux"]))
    # The cut on the title rather than in a footnote: the innermost curve is
    # not the axis and the outermost is not the boundary, and a picture of
    # nested closed curves does not say so by itself.
    ax.set_title(f"traced surfaces, $\\Psi_N \\in [{inner:.3f}, {outer:.3f}]$\n"
                 "neither the axis nor the boundary is in the family",
                 fontsize="small")
    ax.set_xlabel("R  [m]")
    ax.set_ylabel("Z  [m]")
    ax.set_aspect("equal")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best", fontsize="x-small")


def draw_surface_profiles(axes, data, panels, columns):
    """The flux-surface averages against rho, with the cut shown as a cut."""
    from matplotlib.patches import Patch

    rho = data["rho"]
    band = data["band"]
    inner, outer = float(rho.min()), float(rho.max())

    for index, (ax, panel) in enumerate(zip(axes, panels)):
        values = panel["values"]
        ax.plot(rho, values, color="C0", linewidth=1.0, marker=".",
                markersize=3.5, zorder=3)
        if band.any():
            # `band( flux )` and not the node mask, and that is the right flag
            # HERE: an average is one number over a whole surface, so a single
            # continued node contaminates all of it.  How much is continued is
            # the band-fraction panel's job.
            ax.plot(rho[band], values[band], linestyle="none", marker="o",
                    markersize=5.0, markerfacecolor="none",
                    markeredgecolor="C3", markeredgewidth=1.0, zorder=4)

        # THE CUT, DRAWN AS A CUT.  The axes span the whole of rho so that the
        # two ends MEQ does not answer at are visible as gaps rather than
        # being cropped away.  A query there is refused, not extrapolated.
        ax.axvspan(0.0, inner, color="0.88", linewidth=0.0, zorder=0)
        ax.axvspan(outer, 1.0, color="0.88", linewidth=0.0, zorder=0)
        ax.set_xlim(0.0, 1.0)

        ax.set_title(panel["name"], fontsize="small")
        if panel["units"] and panel["units"] != "1":
            ax.set_ylabel(f"[{panel['units']}]", fontsize="x-small")
        elif "derived" in panel:
            ax.set_ylabel("fraction of nodes", fontsize="x-small")
        ax.tick_params(labelsize="x-small")

        if index < columns:
            # Psi_N on the top of the top row, because the difference between
            # the two labels is the whole of Zernike.hpp's argument and a
            # reader who wants Psi_N should not be squaring by eye.
            secondary = ax.secondary_xaxis(
                "top", functions=(lambda r: r**2,
                                  lambda p: np.sqrt(np.clip(p, 0.0, None))))
            secondary.set_xlabel(r"$\Psi_N$", fontsize="x-small")
            secondary.tick_params(labelsize="xx-small")
        if index >= len(panels) - columns:
            ax.set_xlabel(r"$\rho = \sqrt{\Psi_N}$", fontsize="x-small")

    if panels:
        handles = [Patch(facecolor="0.88", label="outside the cut")]
        if band.any():
            handles.append(band_legend_entry())
        axes[0].legend(handles=handles, loc="best", fontsize="xx-small")


def band_legend_entry():
    """A legend entry for the band ring, built without touching an axes."""
    from matplotlib.lines import Line2D

    return Line2D([], [], linestyle="none", marker="o", markersize=5.0,
                  markerfacecolor="none", markeredgecolor="C3",
                  markeredgewidth=1.0, label="crosses the band")


def draw_surface_figure(plt, data, what, requested):
    """Lay the flux-surface file out: the plane, the averages, or both."""
    panels = profile_panels(data, requested) if what != "geometry" else []

    if what == "geometry":
        figure, ax = plt.subplots(figsize=(5.5, 7.0))
        draw_surface_geometry(ax, data)
    elif what == "profiles":
        columns = 3
        rows = -(-len(panels)//columns)
        figure, grid = plt.subplots(rows, columns, squeeze=False,
                                    figsize=(3.4*columns, 2.7*rows))
        axes = list(grid.ravel())
        draw_surface_profiles(axes[:len(panels)], data, panels, columns)
        for spare in axes[len(panels):]:
            spare.set_axis_off()
    else:
        # The poloidal plane is drawn at equal aspect, so it is limited by
        # whichever of its two allowances is smaller.  Three columns of panels
        # keeps the grid short enough that the plane is not squeezed into a
        # sliver beside a very tall stack.
        columns = 3
        rows = -(-len(panels)//columns)
        figure = plt.figure(figsize=(4.6 + 3.3*columns, 2.9*rows))
        grid = figure.add_gridspec(rows, columns + 1,
                                   width_ratios=[2.0] + [1.0]*columns)
        draw_surface_geometry(figure.add_subplot(grid[:, 0]), data)
        axes = [figure.add_subplot(grid[i//columns, 1 + i % columns])
                for i in range(len(panels))]
        draw_surface_profiles(axes, data, panels, columns)

    title = data["attributes"].get("title", "MEQ flux surfaces")
    # About eleven characters to the inch at this font size, measured off a
    # rendered caption rather than assumed.
    line = wrapped(surface_subtitle(data), int(11*figure.get_figwidth()))
    figure.suptitle(f"{title}\n{line}" if line else title, fontsize="medium")
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    return figure


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot a MEQ equilibrium from its NetCDF output.")
    parser.add_argument("file", help="the .nc file MEQ wrote, of either kind")
    parser.add_argument("--what", default="both",
                        choices=("surfaces", "field", "geometry", "profiles",
                                 "both"),
                        help="grid file: flux surfaces, |B_pol| or both; "
                             "flux-surface file: geometry, profiles or both "
                             "(default: both, which reads on either)")
    parser.add_argument("--levels", type=int, default=25,
                        help="number of psi contours, grid file only "
                             "(default: 25)")
    parser.add_argument("--profiles", default=None,
                        help="flux-surface file only: a comma-separated list "
                             "of variables to plot against rho, or `all` to "
                             "include the per-surface diagnostics (default: "
                             "every average the file carries)")
    parser.add_argument("-o", "--output", default=None,
                        help="write a figure here instead of showing it")
    parser.add_argument("--show", action="store_true",
                        help="open a window even when -o is given")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args(argv)

    data = read(args.file)

    # The two files answer different questions, so a --what meant for one of
    # them is refused rather than silently ignored: the file kind is not
    # something the caller chose, and a picture that quietly became a
    # different picture is the failure this whole script is written against.
    for kind, wrong, right in (("grid", ("geometry", "profiles"), "surfaces"),
                               ("surfaces", ("surfaces", "field"), "geometry")):
        if data["kind"] == kind and args.what in wrong:
            sys.exit(f"plot_equilibrium: --what {args.what} is for the "
                     f"{'flux-surface' if kind == 'grid' else '(R, Z) grid'} "
                     f"file; {args.file} is the "
                     f"{'(R, Z) grid' if kind == 'grid' else 'flux-surface'} "
                     f"file. Try --what {right}, or --what both.")

    import matplotlib
    # Choose the backend BEFORE pyplot is imported, and only when there is
    # nothing to show: on a headless machine the default backend fails at
    # import rather than at draw time, which reads as a broken script.
    if args.output and not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if data["kind"] == "surfaces":
        figure = draw_surface_figure(plt, data, args.what, args.profiles)
    else:
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
