#!/usr/bin/env python
"""MEQ against CHEASE on one case, through the norm that already exists.

**NO THIRD ERROR NORM.**  `tools/freegs4e-benchmark/compare.py` defines the
relative L2 of psi over the nodes both sides can answer for, in MEQ's gauge;
`tools/desc-benchmark/compare_desc.py` wraps it for a second code; this wraps
that.  The band is dropped wherever MEQ is an arm, from the `.nc`'s own
`extrapolated` mask, because `compare.py` does it.

Three pairings, and they answer different questions:

    MEQ    vs freegs4e   what MEASUREMENTS.md M-147 publishes
    CHEASE vs freegs4e   the same question of the other code
    MEQ    vs CHEASE     the one that needs no third party

**THE THIRD IS A FLOOR AND NOT AN ACCURACY**, for the reason `compare_desc.py`
argues at length: two codes agreeing to X does not make either accurate to X.
What makes THIS pairing sharper than the DESC one is that the two codes are
solving the same equation from the same two tables on the same analytic curve,
so a disagreement is discretisation on one side or the other and not a
conversion -- which is exactly why the difference is expected to FALL under
refinement, and why `converge.py` exists to check that it does.

CHEASE's psi outside Gamma is an extrapolation with a floor clamp
(`src-f90/psibox.f90:388-468`) and is written as NaN by `convert_chease.py`, so
it never reaches a norm.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "desc-benchmark"))
sys.path.insert(0, os.path.join(HERE, "..", "freegs4e-benchmark"))
sys.path.insert(0, HERE)

import compare_desc                                            # noqa: E402
import convert_desc as convert                                 # noqa: E402
from convert_desc import CASES, MEQ, REFDIR                    # noqa: E402


def meta_path(stem):
	return os.path.join(MEQ, "examples", f"{stem}-meta.json")


def reference(stem):
	return os.path.join(REFDIR, f"{CASES[stem]}.npz")


def pairings(stem, chease_npz, meq_nc=None, core=None):
	"""The three numbers, or the one that does not need a MEQ run."""
	meta = json.load(open(meta_path(stem)))
	out = {}
	out["chease_vs_freegs4e"] = compare_desc.grid_pair(
		chease_npz, reference(stem), meta, core)
	if meq_nc:
		out["meq_vs_freegs4e"] = compare_desc.meq_pair(
			reference(stem), meq_nc, meta_path(stem), core=core)
		out["meq_vs_chease"] = compare_desc.meq_pair(
			chease_npz, meq_nc, meta_path(stem), core=core)
	return out


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(CASES))
	ap.add_argument("--chease", required=True, help="a cheaserun.py .npz")
	ap.add_argument("--meq", default="", help="a MEQ .nc")
	ap.add_argument("--core", type=float, default=None)
	ap.add_argument("--json", default="")
	args = ap.parse_args()

	out = pairings(args.stem, args.chease, args.meq, args.core)
	print(f"  {args.stem}")
	for key, row in out.items():
		print(f"    {key:20s} nodes {row['nodes']:7d}  rel L2 {row['rel_l2']:.4e}"
		      f"  rel Linf {row['rel_linf']:.4e}")
	if args.json:
		json.dump(out, open(args.json, "w"), indent=2)
	return 0


if __name__ == "__main__":
	sys.exit(main())
