#!/usr/bin/env python
"""One CHEASE solve of one MEQ case, with the amplitude closed self-consistently.

    ../freegs4e-benchmark/venv/bin/python cheaserun.py fixed-h-circular --ns 40 --nt 40

writes `runs/<stem>-chease-NS<ns>NT<nt>.npz` -- 1-D `R`, `Z` and `psi[ iR, iZ ]`
in the freegs4e gauge, which is the layout `tools/freegs4e-benchmark/compare.py`
and `tools/desc-benchmark/compare_desc.py` already read -- plus a `.json` of
everything the run reported.

**THE AMPLITUDE IS SOLVED FOR, NOT SUPPLIED.**  `convert_chease.py`'s module
docstring has the derivation; the mechanics are that CHEASE is run at a trial
`A0`, returns `A_C`, and the answer is `A = sqrt( A0 * A_C )` exactly, for any
`A0`.  This runs it at least twice and stops when the returned `A_C` reproduces
its own `A0`.  **The convergence of that loop is itself a check on the
conversion**: the reciprocal law `A_C( A0 ) = K/A0` holds only if the profiles
really were scaled by `1/A0` and by nothing else, so a loop that wanders says
something upstream is wrong.  `--trial` starts it somewhere else if you want to
see that it does not matter; the default is 1.0, which is 16x away from the
answer on `fixed-h-circular` and still converges in two.

**CHEASE WRITES INTO ITS WORKING DIRECTORY AND READS `EXPEQ` AND
`chease_namelist` FROM IT**, both by hard-coded name (`src-f90/iodisk.f90:361`,
`src-f90/chease_prog.f90:116`).  There is no argv.  So every run gets its own
directory and `cwd` is how it is told anything.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import convert_chease as conv                                   # noqa: E402
from convert_chease import CASES, MEQ, CHEASE, load_case        # noqa: E402


def run_once(case, psi_ax, rundir, ns=40, nt=40, npsi=100, nchi=100, niso=100,
             nrbox=513, nzbox=513, box=None, nice=10, timeout=1800):
	"""One CHEASE invocation at one trial axis flux.  Returns its EQDSK.

	`niso` is forced up to `npsi` because CHEASE aborts with "NISO IN NIN TOO
	LARGE" the other way round (`iodisk.f90:109`), and the two are separate
	namelist keys that a caller sweeping one will forget about.
	"""
	os.makedirs(rundir, exist_ok=True)
	niso = max(niso, npsi)
	open(os.path.join(rundir, "EXPEQ"), "w").write(
		conv.expeq_text(case, psi_ax))
	open(os.path.join(rundir, "chease_namelist"), "w").write(
		conv.namelist_text(case, psi_ax, ns=ns, nt=nt, npsi=npsi, nchi=nchi,
		                   niso=niso, nrbox=nrbox, nzbox=nzbox, box=box))
	started = time.time()
	with open(os.path.join(rundir, "o.chease"), "w") as out, \
	     open(os.path.join(rundir, "err.chease"), "w") as err:
		code = subprocess.call(["nice", "-n", str(nice), CHEASE],
		                       cwd=rundir, stdout=out, stderr=err,
		                       timeout=timeout)
	wall = time.time() - started
	eqdsk = os.path.join(rundir, "EQDSK_COCOS_02.OUT")
	if code != 0 or not os.path.exists(eqdsk):
		tail = open(os.path.join(rundir, "o.chease")).read()[-2000:]
		raise RuntimeError(f"CHEASE failed in {rundir} ( exit {code} )\n{tail}")
	g = conv.read_eqdsk(eqdsk)
	return dict(eqdsk_path=eqdsk, g=g, wall=wall, rundir=rundir,
	            psi_ax_trial=psi_ax,
	            # CHEASE's psi is zero on Gamma and negative on the axis, so the
	            # amplitude is |SIMAG| and `SIBRY` is a literal zero.
	            psi_ax_returned=abs(float(g["simagx"])),
	            psi_bndry=float(g["sibdry"]),
	            ip=float(g["cpasma"]),
	            r_axis=float(g["rmagx"]), z_axis=float(g["zmagx"]),
	            q_axis=float(g["qpsi"][0]), q_edge=float(g["qpsi"][-1]),
	            g_axis=float(g["fpol"][0]), g_edge=float(g["fpol"][-1]),
	            p_axis=float(g["pres"][0]), p_edge=float(g["pres"][-1]))


def solve(case, rundir, trial=1.0, tol=1.0e-9, max_sweeps=6, **kw):
	"""The amplitude fixed point.  Returns the converged run plus its history.

	The update is `A <- sqrt( A0 * A_C( A0 ) )`, which is the EXACT solution of
	the reciprocal law rather than a relaxation of it, so a correct conversion
	lands in one step and the second sweep only confirms it.  Compare
	`tools/desc-benchmark`'s `--self-consistent`, which iterates a map that is
	NOT exact and needs damping to avoid a period-2 limit cycle; this one needs
	none, and that difference is a property of the two conversions rather than
	of the two codes.
	"""
	history = []
	A0 = trial
	run = None
	for sweep in range(max_sweeps):
		run = run_once(case, A0, os.path.join(rundir, f"sweep{sweep}"), **kw)
		AC = run["psi_ax_returned"]
		history.append(dict(sweep=sweep, trial=A0, returned=AC,
		                    residual=abs(AC - A0) / AC, wall=run["wall"]))
		if abs(AC - A0) / AC < tol:
			break
		A0 = float(np.sqrt(A0 * AC))
	run["history"] = history
	run["psi_ax"] = run["psi_ax_returned"]
	run["sweeps"] = len(history)
	run["wall_total"] = sum(h["wall"] for h in history)
	return run


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(CASES))
	ap.add_argument("--ns", type=int, default=40, help="radial FE mesh")
	ap.add_argument("--nt", type=int, default=40, help="poloidal FE mesh")
	ap.add_argument("--npsi", type=int, default=100)
	ap.add_argument("--nchi", type=int, default=100)
	ap.add_argument("--niso", type=int, default=100)
	ap.add_argument("--nrbox", type=int, default=513)
	ap.add_argument("--nzbox", type=int, default=513)
	ap.add_argument("--match-box", action="store_true",
	                help="put the EQDSK box on the case's own mesh rectangle, "
	                     "so the grid comparison needs no interpolation")
	ap.add_argument("--trial", type=float, default=1.0)
	ap.add_argument("--out", default="", help="the .npz; default runs/<stem>-...")
	ap.add_argument("--rundir", default="", help="scratch; default <tempdir>/chease-bench/<tag>")
	ap.add_argument("--keep", action="store_true", help="keep the scratch dir")
	args = ap.parse_args()

	case = load_case(args.stem)
	meta = case["meta"]
	box = None
	if args.match_box:
		b = meta["box"]
		box = (b["rmin"], b["rmax"], b["zmin"], b["zmax"])

	tag = f"{args.stem}-chease-NS{args.ns}NT{args.nt}NPSI{args.npsi}"
	# TEMPDIR AND NOT A LITERAL /tmp PATH: CHEASE writes into its working
	# directory by hard-coded name, so this is real scratch and belongs
	# wherever the machine puts scratch.  --rundir overrides it.
	rundir = args.rundir or os.path.join(tempfile.gettempdir(),
	                                     "chease-bench", tag)
	if os.path.isdir(rundir):
		shutil.rmtree(rundir)
	run = solve(case, rundir, trial=args.trial, ns=args.ns, nt=args.nt,
	            npsi=args.npsi, nchi=args.nchi, niso=args.niso,
	            nrbox=args.nrbox, nzbox=args.nzbox, box=box)

	grid = conv.eqdsk_to_grid(run["eqdsk_path"], meta)
	out = args.out or os.path.join(HERE, "runs", f"{tag}.npz")
	os.makedirs(os.path.dirname(out), exist_ok=True)
	np.savez(out, R=grid["R"], Z=grid["Z"], psi=grid["psi"])

	record = {k: v for k, v in run.items() if k not in ("g", "eqdsk_path")}
	record.update(stem=args.stem, ns=args.ns, nt=args.nt, npsi=args.npsi,
	              nchi=args.nchi, niso=args.niso, nrbox=args.nrbox,
	              nzbox=args.nzbox, npz=out,
	              # NOT an input to anything -- printed so the closure can be
	              # scored against a number it never saw.
	              reference_psi_ax=float(meta["psi_axis"])
	              - float(meta["psi_surface"]))
	json.dump(record, open(out.replace(".npz", ".json"), "w"), indent=2)

	print(f"  {args.stem}  NS {args.ns}  NT {args.nt}  NPSI {args.npsi}")
	for h in run["history"]:
		print(f"    sweep {h['sweep']}  trial {h['trial']:.10e}  "
		      f"returned {h['returned']:.10e}  residual {h['residual']:.3e}"
		      f"   ( {h['wall']:.1f} s, ROUGH )")
	print(f"    psi_ax  {run['psi_ax']:.10e}   reference "
	      f"{record['reference_psi_ax']:.10e}   relative "
	      f"{abs(run['psi_ax'] - record['reference_psi_ax']) / record['reference_psi_ax']:.4e}")
	print(f"    Ip {run['ip']:.6e} A   axis ( {run['r_axis']:.6f}, "
	      f"{run['z_axis']:.6f} )   q0 {run['q_axis']:.4f}  qa {run['q_edge']:.4f}")
	print(f"    g on axis {run['g_axis']:.8f}   on Gamma {run['g_edge']:.8f}"
	      f"   ( g_at_gamma {meta['g_at_gamma']:.8f} )")
	print(f"    wrote {out}")
	if not args.keep:
		shutil.rmtree(rundir, ignore_errors=True)
	return 0


if __name__ == "__main__":
	sys.exit(main())
