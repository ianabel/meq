#!/usr/bin/env python
"""ONE ERROR NORM, THREE CODES, AND NO CODE SCORED AGAINST ITSELF.

Every number the race reports is a relative L2 of the poloidal flux over the
( R, z ) nodes both sides can answer for, in MEQ's gauge -- psi = 0 on Gamma.
Three pairings exist and they answer different questions:

  MEQ  vs freegs4e   what MEASUREMENTS.md M-147 publishes, 5.0e-05 to 1.1e-04
  DESC vs freegs4e   the same question of the other code
  MEQ  vs DESC       the one that needs no third party, and the one that says
                     whether the CONVERSION in convert_desc.py is right at all

**THE THIRD IS THE FLOOR AND IT IS NOT AN ACCURACY.** Two codes agreeing to X
does not make either of them accurate to X, and neither one is the truth: MEQ
is a fixed-boundary HDG solve on the exact MXH curve, DESC is an inverse
spectral solve of a problem whose PROFILES were re-expressed through an
approximate flux-label map, and freegs4e is a 257^2 finite-difference free
boundary solve whose psi_n = 0.95 contour is only the MXH curve to 6.0e-05.
What the three pairings do is bound the disagreement from three directions, so
that a race target BELOW the mutual floor is refused rather than reported.

THE BAND IS DROPPED WHEREVER MEQ IS AN ARM.  `compare.py` already does it, from
the .nc's own `extrapolated` mask, and this calls that function rather than
writing a second norm that could drift from it.  For a DESC-against-freegs4e
pair there is no band -- neither is an unfitted mesh -- and the mask is
"where DESC answered", which is the interior of Gamma.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "freegs4e-benchmark"))
sys.path.insert(0, HERE)
import compare as meq_compare                                  # noqa: E402
import convert_desc as convert                                 # noqa: E402


def grid_pair(npz_a, npz_b, meta, core=None):
	"""Relative L2 and Linf between two ( R, z ) grid files, in MEQ's gauge.

	Both files are in the reference's own layout -- psi[ iR, iZ ] and the
	freegs4e gauge -- so the shift to MEQ's gauge is the same subtraction for
	both and cancels out of the difference; it is applied anyway because the
	SCALE the norm is relative to is psi in MEQ's gauge, which is the quantity
	the race's targets are quoted against.
	"""
	a = np.load(npz_a, allow_pickle=True)
	b = np.load(npz_b, allow_pickle=True)
	Ra, Za = np.array(a["R"], float), np.array(a["Z"], float)
	Rb, Zb = np.array(b["R"], float), np.array(b["Z"], float)
	shift = float(meta["psi_surface"])
	pa = np.array(a["psi"], float) - shift
	pb = np.array(b["psi"], float) - shift
	if not (np.array_equal(Ra, Rb) and np.array_equal(Za, Zb)):
		from scipy.interpolate import RegularGridInterpolator
		interp = RegularGridInterpolator((Rb, Zb), pb, bounds_error=False,
		                                 fill_value=np.nan)
		RR, ZZ = np.meshgrid(Ra, Za, indexing="ij")
		pb = interp(np.stack([RR.ravel(), ZZ.ravel()], -1)).reshape(RR.shape)
	use = np.isfinite(pa) & np.isfinite(pb)
	if core is not None:
		psi_ax = float(meta["psi_axis"]) - shift
		use &= (pa >= (1.0 - core) * psi_ax)
	if use.sum() == 0:
		raise RuntimeError(f"{npz_a} and {npz_b} share no node")
	d = pa[use] - pb[use]
	scale = float(np.abs(pa[use]).max())
	return dict(nodes=int(use.sum()),
	            rel_l2=float(np.sqrt(np.mean(d ** 2)) / scale),
	            rel_linf=float(np.abs(d).max() / scale), scale=scale)


def meq_pair(npz, nc, meta_path, core=None, free_boundary=False):
	"""MEQ's .nc against any grid file, through compare.py's own norm.

	`free_boundary` PICKS THE GAUGE AND IT IS NOT COSMETIC.  On the fixed
	ladder MEQ is handed the reference's psi_n = 0.95 surface with psi = 0 on
	it, so its psi is the reference's shifted by that surface's flux and
	`compare.py` subtracts `psi_surface`.  A free-boundary MEQ run shares the
	reference's gauge exactly -- both solve the same exterior problem, psi -> 0
	at infinity -- so there is nothing to subtract, and subtracting anyway
	offsets one arm of the difference by the whole boundary flux.  Measured on
	the limited machine: 1.36e-01 with the shift applied against 1.4e-04
	without, on a case whose scalars agree to 4.6e-05.
	"""
	return meq_compare.compare(npz, nc, None if free_boundary else meta_path,
	                           free_boundary=free_boundary, core=core)


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(convert.CASES))
	ap.add_argument("--desc", required=True, help="a descrun.py .npz")
	ap.add_argument("--meq", default="", help="a MEQ .nc")
	ap.add_argument("--core", type=float, default=None)
	ap.add_argument("--json", default="")
	args = ap.parse_args()

	meta = json.load(open(os.path.join(convert.MEQ, "examples",
	                                   f"{args.stem}-meta.json")))
	reference = os.path.join(convert.REFDIRS.get(args.stem, convert.REFDIR), f"{convert.CASES[args.stem]}.npz")
	out = {}
	out["desc_vs_freegs4e"] = grid_pair(args.desc, reference, meta, args.core)
	free = bool(meta.get("free_boundary"))
	meta_path = os.path.join(convert.MEQ, "examples",
	                         f"{args.stem}-meta.json")
	if args.meq:
		out["meq_vs_freegs4e"] = meq_pair(reference, args.meq, meta_path,
		                                  core=args.core, free_boundary=free)
		out["meq_vs_desc"] = meq_pair(args.desc, args.meq, meta_path,
		                              core=args.core, free_boundary=free)
	print(f"  {args.stem}")
	for key, row in out.items():
		print(f"    {key:20s} nodes {row['nodes']:7d}  rel L2 {row['rel_l2']:.4e}"
		      f"  rel Linf {row['rel_linf']:.4e}")
	if args.json:
		json.dump(out, open(args.json, "w"), indent=2)
	return 0


if __name__ == "__main__":
	sys.exit(main())
