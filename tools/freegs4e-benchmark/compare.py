"""Compare MEQ's .nc against the freegs4e reference it was built from."""
import json
import os
import numpy as np
from netCDF4 import Dataset


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
                rel_l2=float(np.sqrt(np.mean(diff ** 2)) / scale))


if __name__ == "__main__":
    import sys
    rows = []
    print(f"\n  MEQ AGAINST freegs4e\n")
    print(f"    {'case':22s} {'nodes':>7s} {'band':>6s} {'rel L2':>11s} "
          f"{'rel Linf':>11s} {'scale':>11s}")
    for npz in sys.argv[1:]:
        stem = os.path.splitext(os.path.basename(npz))[0]
        d = os.path.dirname(npz) or "."
        nc = os.path.join(d, f"{stem}.nc")
        meta = os.path.join(d, f"{stem}-meta.json")
        if not os.path.exists(nc):
            print(f"    {stem:22s} {'-':>7s} MEQ produced no .nc")
            continue
        r = compare(npz, nc, meta)
        rows.append((stem, r))
        print(f"    {stem:22s} {r['nodes']:7d} {r['dropped_band']:6d} "
              f"{r['rel_l2']:11.3e} {r['rel_linf']:11.3e} {r['scale']:11.3e}")
    if rows:
        worst = max(r['rel_l2'] for _, r in rows)
        print(f"\n    worst relative L2 across {len(rows)} cases: {worst:.3e}\n")
