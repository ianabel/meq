#!/usr/bin/env python
"""Does MEQ's answer and CHEASE's approach each other under refinement?

    ../freegs4e-benchmark/venv/bin/python converge.py fixed-h-circular

**THIS IS THE ACCEPTANCE AND NOT THE DEMONSTRATION.**  `CLAUDE.md`'s *Testing
stance* records the hazard that makes a single-resolution agreement worth
nothing: *a wrong sign convention converges, at the right rate, to the wrong
function*.  Two codes agreeing at one resolution can agree for that reason.  Two
codes whose difference FALLS as each is refined cannot -- a conversion error is
a fixed offset and does not care about either mesh.

So this runs a matrix, not a pair: MEQ's ladder down one axis, CHEASE's along
the other, and the relative L2 between them in every cell.  Three readings come
out of one matrix and they are different claims:

  * **a ROW** -- MEQ fixed, CHEASE refined -- falls until it reaches MEQ's own
    discretisation error and then stops.  Where it stops IS MEQ's error at that
    row, measured without a reference.
  * **a COLUMN** is the same statement with the codes swapped.
  * **the DIAGONAL** is the only one that can fall all the way, and where IT
    stops is the floor of the comparison: the MXH curve, the two profile tables
    and the two scalars, i.e. the problem statement rather than either solver.

**THE KNOBS.**  MEQ's are `[discretisation] PolynomialDegree` and `[mesh]
RefinementLevels`, taken through `race_desc.meq_config` so this ladder and the
DESC one are the same ladder.  CHEASE's is `NS x NT`, the bicubic Hermite finite
element mesh of the Grad-Shafranov solve itself.  `NPSI`, `NCHI` and `NISO` are
held FIXED across the sweep and checked separately with `--mapping-control`,
because they size the flux-surface MAPPING rather than the solve and sweeping
them together would leave two knobs moving at once.

**EVERY SECOND PRINTED IS ROUGH.**  The machine was not quiet.  The error
columns do not depend on load and are reproducible; the seconds are there to
say what a cell costs, and a real timing needs `CLAUDE.md`'s idle machine.
"""
import argparse
import json
import os
import shutil
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "desc-benchmark"))
sys.path.insert(0, os.path.join(HERE, "..", "freegs4e-benchmark"))
sys.path.insert(0, HERE)

import race_desc                                               # noqa: E402
import compare_chease                                          # noqa: E402
import cheaserun                                               # noqa: E402
import convert_chease as conv                                  # noqa: E402
from convert_desc import CASES, MEQ, load_case                 # noqa: E402


def meq_rungs(text):
	"""`"1:0,2:0,3:0"` -> [ ( degree, refinement ), ... ]."""
	out = []
	for item in text.split(","):
		k, _, R = item.partition(":")
		out.append((int(k), int(R or 0)))
	return out


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(CASES))
	ap.add_argument("--meq", default="1:0,2:0,3:0,2:1,3:1",
	                help="degree:refinement rungs")
	ap.add_argument("--chease", default="20,30,40,60,80",
	                help="NS = NT rungs")
	ap.add_argument("--npsi", type=int, default=200,
	                help="the MAPPING mesh, held fixed across the sweep")
	ap.add_argument("--nrbox", type=int, default=999,
	                help="CHEASE's EQDSK box.  999 AND NOT MORE: the header "
	                     "format is 3I4 and a four-digit box size runs the "
	                     "three integers together.  It is the LIMITING "
	                     "instrument in the deepest cells -- see the README")
	ap.add_argument("--grid", type=int, default=257,
	                help="MEQ's output grid, which is the comparison grid")
	ap.add_argument("--mapping-control", action="store_true",
	                help="re-run the finest CHEASE rung at twice NPSI and at "
	                     "twice NRBOX, to show neither is what is being measured")
	ap.add_argument("--out",
	                default=os.path.join(tempfile.gettempdir(),
	                                     "chease-converge"))
	args = ap.parse_args()

	case = load_case(args.stem)
	os.makedirs(args.out, exist_ok=True)
	report = dict(stem=args.stem, meq=args.meq, chease=args.chease,
	              npsi=args.npsi, nrbox=args.nrbox, grid=args.grid)

	# ---------------------------------------------------------- the CHEASE arm
	chease = []
	trial = 1.0
	for ns in [int(x) for x in args.chease.split(",")]:
		tag = f"{args.stem}-chease-NS{ns}"
		rundir = os.path.join(args.out, "scratch", tag)
		shutil.rmtree(rundir, ignore_errors=True)
		run = cheaserun.solve(case, rundir, trial=trial, ns=ns, nt=ns,
		                      npsi=args.npsi, nchi=args.npsi, niso=args.npsi,
		                      nrbox=args.nrbox, nzbox=args.nrbox)
		# SEEDING THE NEXT RUNG COSTS NOTHING AND CHANGES NOTHING.  The fixed
		# point is exact in one step from ANY trial -- demonstrated at
		# trial = 1.0, sixteen times away -- so a warm trial only saves a sweep.
		trial = run["psi_ax"]
		grid = conv.eqdsk_to_grid(run["eqdsk_path"], case["meta"])
		npz = os.path.join(args.out, f"{tag}.npz")
		np.savez(npz, R=grid["R"], Z=grid["Z"], psi=grid["psi"])
		chease.append(dict(ns=ns, npz=npz, psi_ax=run["psi_ax"], ip=run["ip"],
		                   q_axis=run["q_axis"], r_axis=run["r_axis"],
		                   sweeps=run["sweeps"], wall=run["wall_total"]))
		shutil.rmtree(rundir, ignore_errors=True)
		print(f"  CHEASE NS=NT={ns:3d}  psi_ax {run['psi_ax']:.10e}  "
		      f"Ip {run['ip']:.6e}  {run['sweeps']} sweeps  "
		      f"{run['wall_total']:.1f} s ROUGH")
	report["chease_runs"] = chease

	# ------------------------------------------------------------- the MEQ arm
	meq = []
	for degree, refine in meq_rungs(args.meq):
		config, nc = race_desc.meq_config(args.stem, args.out, degree, refine,
		                                  grid=args.grid)
		run = race_desc.run_meq(config)
		meq.append(dict(degree=degree, refine=refine, nc=nc, ok=run["ok"],
		                elements=run["elements"], newton=run["newton"],
		                wall=run["seconds"]))
		print(f"  MEQ    k={degree} r={refine}    {run['elements']:6d} elements  "
		      f"{run['newton']} Newton  {run['seconds']:.1f} s ROUGH"
		      f"{'' if run['ok'] else '   FAILED: ' + run['tail']}")
	report["meq_runs"] = meq

	# ------------------------------------------------------------- the matrix
	matrix = {}
	for m in meq:
		if not m["ok"]:
			continue
		row = {}
		for c in chease:
			key = f"NS{c['ns']}"
			row[key] = compare_chease.pairings(args.stem, c["npz"], m["nc"])
		matrix[f"k{m['degree']}r{m['refine']}"] = row
	report["matrix"] = matrix

	print()
	print("  MEQ vs CHEASE, relative L2 of psi over MEQ's own nodes")
	head = "".join(f"{'NS' + str(c['ns']):>12s}" for c in chease)
	print(f"    {'':10s}{head}")
	for name, row in matrix.items():
		cells = "".join(f"{row['NS' + str(c['ns'])]['meq_vs_chease']['rel_l2']:12.3e}"
		                for c in chease)
		print(f"    {name:10s}{cells}")

	print()
	print("  each arm against the freegs4e reference neither of them saw")
	for m in meq:
		if not m["ok"]:
			continue
		key = f"k{m['degree']}r{m['refine']}"
		v = matrix[key][f"NS{chease[0]['ns']}"]["meq_vs_freegs4e"]["rel_l2"]
		print(f"    MEQ    {key:10s} {v:.4e}")
	for c in chease:
		first = next(iter(matrix))
		v = matrix[first][f"NS{c['ns']}"]["chease_vs_freegs4e"]["rel_l2"]
		print(f"    CHEASE NS{c['ns']:<9d} {v:.4e}")

	# ---------------------------------------------------- the two instrument
	# controls, so that the matrix above is a statement about the SOLVES
	if args.mapping_control:
		ns = int(args.chease.split(",")[-1])
		# **THE SECOND CONTROL HALVES THE BOX RATHER THAN DOUBLING IT**,
		# because `NRBOX >= 1000` writes an unparseable EQDSK header (see
		# `convert_chease.namelist_text`).  Halving is the better experiment
		# anyway: the error should GROW by four, and a control that grows by
		# the predicted factor says more than one that shrinks by it.
		for label, kw in (("NPSI x2", dict(npsi=2 * args.npsi,
		                                   nchi=2 * args.npsi,
		                                   niso=2 * args.npsi,
		                                   nrbox=args.nrbox)),
		                  ("NRBOX /2", dict(npsi=args.npsi, nchi=args.npsi,
		                                    niso=args.npsi,
		                                    nrbox=(args.nrbox + 1) // 2))):
			rundir = os.path.join(args.out, "scratch", "control")
			shutil.rmtree(rundir, ignore_errors=True)
			run = cheaserun.solve(case, rundir, trial=trial, ns=ns, nt=ns,
			                      nzbox=kw.get("nrbox"), **kw)
			grid = conv.eqdsk_to_grid(run["eqdsk_path"], case["meta"])
			npz = os.path.join(args.out, f"control-{label.replace(' ', '')}.npz")
			np.savez(npz, R=grid["R"], Z=grid["Z"], psi=grid["psi"])
			best = [m for m in meq if m["ok"]][-1]
			p = compare_chease.pairings(args.stem, npz, best["nc"])
			base = matrix[f"k{best['degree']}r{best['refine']}"][f"NS{ns}"]
			print(f"  control {label:9s} NS={ns}: MEQ vs CHEASE "
			      f"{p['meq_vs_chease']['rel_l2']:.4e}   against "
			      f"{base['meq_vs_chease']['rel_l2']:.4e} at the sweep's value")
			report.setdefault("controls", {})[label] = dict(
				rel_l2=p["meq_vs_chease"]["rel_l2"],
				baseline=base["meq_vs_chease"]["rel_l2"],
				psi_ax=run["psi_ax"])
			shutil.rmtree(rundir, ignore_errors=True)

	path = os.path.join(args.out, f"{args.stem}-converge.json")
	json.dump(report, open(path, "w"), indent=2)
	print(f"\n  wrote {path}")
	return 0


if __name__ == "__main__":
	sys.exit(main())
