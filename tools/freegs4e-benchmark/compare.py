"""Compare MEQ's .nc against the freegs4e reference it was built from."""
import argparse
import json
import os
import sys
import numpy as np
from netCDF4 import Dataset


# The same one-sided test meq::CriticalPointFinder::checkAxis applies, with the
# same default -- src/meq/CriticalPoints.hpp:531 declares `toleranceIn = 0.10`
# and CriticalPoints.cpp:868 tests `normalisedFlux >= 1.0 - toleranceIn`.
# Copied rather than derived because this script cannot link the library; if the
# library's default moves, this is the line that has to move with it.
AXIS_TOLERANCE = 0.10


def axis_status(ds):
    """Is psi_axis the flux at a magnetic axis, or the largest number in the
    potential vector?

    Psi at the magnetic axis is 1 BY DEFINITION when psi_ax is the axis flux, so
    the driver writes `axis_normalised_flux` -- Psi at the located O-point of
    q_h -- and it must read 1. The test is ONE SIDED and only the low side is a
    defect: psi_ax is the largest NODAL value, and the peak of a polynomial over
    a closed element is at least its largest nodal value, so a healthy field
    approaches 1 from ABOVE.

    Returns ( verdict, flux, note ) with verdict in { "ok", "bad", "absent" }.
    """
    attributes = ds.ncattrs()
    if "axis_normalised_flux" in attributes:
        flux = float(ds.getncattr("axis_normalised_flux"))
        where = ""
        if "axis_r" in attributes and "axis_z" in attributes:
            where = (f" at ( {float(ds.getncattr('axis_r')):.4f}, "
                     f"{float(ds.getncattr('axis_z')):.4f} )")
        if flux >= 1.0 - AXIS_TOLERANCE:
            return "ok", flux, where
        return "bad", flux, where

    # ABSENT, AND THE REASON MATTERS -- the driver writes all three attributes
    # only when the source is normalised AND an O-point was located, so there
    # are three ways to get here and they are not equally bad:
    #
    #   * no `psi_axis` either: the run was not normalised, psi_ax is not an
    #     answer, and there is nothing to judge;
    #   * `psi_axis` present and the axis absent: apps/meq.cpp is explicit that
    #     "their ABSENCE is informative too" -- no interior extremum was found
    #     anywhere, which is the wall-hugging annulus branch;
    #   * or the file simply predates the attribute.
    #
    # The last two cannot be told apart from the file alone, so this WARNS and
    # lets the comparison run rather than refusing: an old .nc must not become
    # unreadable, and refusing would make the two cases indistinguishable in the
    # other direction.
    if "psi_axis" in attributes:
        return "absent", None, ("normalised run, no axis written: either the "
                                "annulus branch (no O-point anywhere) or a "
                                ".nc predating the attribute")
    return "absent", None, "not a normalised run, so psi_axis is not an answer"


def compare(npz_path, nc_path, meta_path):
    ref = np.load(npz_path, allow_pickle=True)
    meta = json.load(open(meta_path))
    # The surface MEQ was actually given, not the separatrix.
    psi_surface = float(meta["psi_surface"])

    with Dataset(nc_path) as ds:
        Rm = np.array(ds["R"][:], float)
        Zm = np.array(ds["Z"][:], float)
        psi_m = np.array(ds["psi"][:], float)          # ( Z, R )
        inside = np.array(ds["inside"][:]).astype(bool)
        extrap = (np.array(ds["extrapolated"][:]).astype(bool)
                  if "extrapolated" in ds.variables
                  else np.zeros_like(inside))
        # Read here rather than in a second open: the caller needs it BEFORE it
        # quotes any of the norms below.
        verdict, axis_flux, axis_note = axis_status(ds)

    Rf = np.array(ref["R"], float)
    Zf = np.array(ref["Z"], float)
    psi_f = np.array(ref["psi"], float)
    # freegs4e stores psi[ iR, iZ ] -- R IS AXIS 0 -- and MEQ stores ( Z, R ).
    # UNCONDITIONAL, because a shape test cannot tell them apart on a square
    # grid and every reference case here is 129 x 129. Asserted rather than
    # inferred: a silent transpose of a nearly up-down symmetric equilibrium
    # gives a plausible wrong answer rather than an obvious one.
    assert psi_f.shape == (len(Rf), len(Zf)), (
        f"reference psi has shape {psi_f.shape}, expected "
        f"( nR, nZ ) = ( {len(Rf)}, {len(Zf)} )")
    psi_f = psi_f.T
    psi_f = psi_f - psi_surface          # into MEQ's gauge: zero on Gamma

    # Bilinear interpolation of the reference onto MEQ's grid.
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((Zf, Rf), psi_f, bounds_error=False,
                                     fill_value=np.nan)
    ZZ, RR = np.meshgrid(Zm, Rm, indexing="ij")
    ref_on_meq = interp(np.stack([ZZ.ravel(), RR.ravel()], axis=-1)).reshape(ZZ.shape)

    # DROP THE BAND. `extrapolated` nodes are continued past Gamma_h by a Taylor
    # step rather than solved, so they are not the same quantity the reference
    # computes there. CLAUDE.md is explicit that a count is not a mask and that
    # these must come out before any error norm.
    use = inside & ~extrap & np.isfinite(ref_on_meq) & np.isfinite(psi_m)
    if use.sum() == 0:
        raise RuntimeError("no comparable nodes")

    diff = psi_m[use] - ref_on_meq[use]
    scale = np.abs(ref_on_meq[use]).max()
    return dict(nodes=int(use.sum()), inside=int(inside.sum()),
                dropped_band=int((inside & extrap).sum()),
                linf=float(np.abs(diff).max()),
                l2=float(np.sqrt(np.mean(diff ** 2))),
                scale=float(scale),
                rel_linf=float(np.abs(diff).max() / scale),
                rel_l2=float(np.sqrt(np.mean(diff ** 2)) / scale),
                axis_verdict=verdict, axis_flux=axis_flux, axis_note=axis_note)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare MEQ's .nc against the freegs4e reference .npz it "
                    "was built from.")
    parser.add_argument("npz", nargs="+",
                        help="reference .npz files; the .nc and -meta.json are "
                             "found beside each")
    parser.add_argument(
        "--allow-bad-axis", action="store_true",
        help="quote the norms even where psi_ax is NOT the flux at a magnetic "
             "axis. You probably do not want this: psi_ax is what the profiles "
             "are normalised by, so a run that fails the check did not solve "
             "the equilibrium [source] describes, and differencing it against "
             "the reference measures two different problems. It exists for "
             "deliberately inspecting a known-bad run.")
    args = parser.parse_args()

    rows = []
    refused = []
    print(f"\n  MEQ AGAINST freegs4e\n")
    print(f"    {'case':22s} {'nodes':>7s} {'band':>6s} {'rel L2':>11s} "
          f"{'rel Linf':>11s} {'scale':>11s}")
    for npz in args.npz:
        stem = os.path.splitext(os.path.basename(npz))[0]
        d = os.path.dirname(npz) or "."
        nc = os.path.join(d, f"{stem}.nc")
        meta = os.path.join(d, f"{stem}-meta.json")
        if not os.path.exists(nc):
            print(f"    {stem:22s} {'-':>7s} MEQ produced no .nc")
            continue
        r = compare(npz, nc, meta)

        # THE AXIS GATE, BEFORE ANY NUMBER IS QUOTED. A run whose psi_ax is not
        # the flux at a magnetic axis solved a different equilibrium from the
        # one [source] describes -- psi_ax is what the profiles are normalised
        # by -- so its norms are a comparison between two different problems and
        # printing them is worse than printing nothing.
        if r["axis_verdict"] == "bad" and not args.allow_bad_axis:
            refused.append((stem, r))
            print(f"    {stem:22s} {'-':>7s} {'-':>6s} {'REFUSED':>11s} "
                  f"{'-':>11s} {'-':>11s}")
            continue
        if r["axis_verdict"] == "absent":
            # stdout carries the table and stderr the warnings, so without the
            # flush the two arrive interleaved wherever the pair is piped to one
            # file -- and this warning is about the row printed next to it.
            # apps/meq.cpp does the same thing at its own axis warning and says
            # so in the same words.
            sys.stdout.flush()
            print(f"    {stem:22s} -- warning: no axis check in the .nc "
                  f"( {r['axis_note']} )", file=sys.stderr)

        rows.append((stem, r))
        flag = "  BAD AXIS" if r["axis_verdict"] == "bad" else ""
        print(f"    {stem:22s} {r['nodes']:7d} {r['dropped_band']:6d} "
              f"{r['rel_l2']:11.3e} {r['rel_linf']:11.3e} {r['scale']:11.3e}"
              f"{flag}")
    if rows:
        worst = max(r['rel_l2'] for _, r in rows)
        print(f"\n    worst relative L2 across {len(rows)} cases: {worst:.3e}\n")

    sys.stdout.flush()              # see the note at the "absent" warning above
    for stem, r in refused:
        print(f"\n  REFUSED {stem}: psi_ax is NOT the flux at a magnetic axis.\n"
              f"    The .nc reports a normalised flux of {r['axis_flux']:.4f} at "
              f"the located O-point{r['axis_note']},\n"
              f"    where a magnetic axis must read 1 and this check accepts "
              f"anything at or above\n"
              f"    {1.0 - AXIS_TOLERANCE:.2f}. psi_ax is the largest NODAL value "
              f"of psi_h and nothing makes that\n"
              f"    an axis, so the profiles were evaluated over a normalised "
              f"flux the plasma never\n"
              f"    reaches: this is not the equilibrium [source] describes and "
              f"differencing it\n"
              f"    against the reference compares two different problems. MEQ "
              f"warns about this on\n"
              f"    stderr when it writes the file. Pass --allow-bad-axis to "
              f"see the norms anyway.",
              file=sys.stderr)
    if refused:
        sys.exit(1)
