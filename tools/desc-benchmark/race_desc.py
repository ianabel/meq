#!/usr/bin/env python
"""MEQ AGAINST DESC ON THE SIX SHIPPED FIXED-BOUNDARY MACHINE CASES.

Two races, because the user asked for both and they answer different questions.

  `--mode wall`      straight wall clock at a MATCHED resolution.  The simpler
                     number, and the one whose definition is a judgement.
  `--mode accuracy`  time to accuracy: the cheapest run of each code reaching
                     each error target, a Pareto front rather than a pair of
                     seconds.  `tools/freegs4e-benchmark/time_to_accuracy.py`
                     is the same instrument against freegsnke and its header is
                     the argument for why this is the one that can DECIDE
                     anything.

WHAT "MATCHED RESOLUTION" MEANS FOR TWO CODES WITH DIFFERENT DISCRETISATIONS,
SAID OUT LOUD BECAUSE IT IS A CHOICE AND NOT A FACT.

MEQ is an HDG method on a triangulation of the ( R, z ) PLANE: the unknowns are
a flux, a potential and a trace, the trace is what is globally coupled, and the
resolution knobs are the element count and the polynomial degree.  DESC is a
spectral method in FLUX COORDINATES: the unknowns are the Fourier-Zernike
coefficients of R( rho, theta ), Z( rho, theta ) and lambda, the surfaces are
nested by construction, and the resolution knobs are the mode numbers L and M.
There is no quantity that means the same thing on both sides.  Three candidate
matchings, all implemented, none privileged:

  `dofs`      MEQ's globally-coupled trace unknowns against DESC's `dim_x`.
              MEQ's count is estimated as faces * ( k + 1 ) with faces taken
              from Euler for a triangulation, ~1.5 * elements; it is an
              ESTIMATE and is printed as one.  The objection: MEQ's
              element-local unknowns are condensed out and cost real time, so
              its trace count understates its work, while DESC's `dim_x` is
              reduced again by the fixed-boundary constraints.  Comparable in
              spirit, not in kind.
  `scale`     the smallest resolved length.  MEQ: h / k with h from the element
              count and the plasma area.  DESC: the minor radius over M, which
              is its finest poloidal feature.  The objection: it charges DESC
              nothing for knowing the surfaces are nested, which is most of
              what an inverse solver buys.
  `accuracy`  the pairing that needs no analogy at all -- and is `--mode
              accuracy` rather than a matching, which is the point.

THE HONEST ANSWER IS THAT ONLY THE THIRD DECIDES ANYTHING, and the first two
are reported because "which is faster at the same resolution" is a question
people ask and refusing to answer it is not better than answering it with the
definition attached.

THE PROTOCOL is `tools/freegs4e-benchmark/race_nke.py`'s and each clause of it
is a measurement this tree has already got wrong once: INTERLEAVED arms, a
MEDIAN of repeats rather than a best-of, the load average checked just BEFORE
each level and not after, and provenance on every row.  Two clauses are new
here:

  DESC IS TIMED AS A FRESH PROCESS.  JAX compiles the residual, the Jacobian
  and the optimiser loop on the first solve of a given shape, and MEQ has no
  equivalent -- it compiled at build time.  Timing DESC's second solve against
  MEQ's only solve would hand DESC a leg MEQ pays too and is not credited for.
  So the wall clock is `subprocess` around `descrun.py`, exactly as MEQ's is
  `subprocess` around the driver, and the WARM figure is reported beside it
  from inside that process rather than instead of it.

  NEITHER CODE IS THE TRUTH.  See compare_desc.py.
"""
import argparse
import json
import os
import re
import statistics
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "freegs4e-benchmark"))
import convert_desc as convert                                 # noqa: E402
import compare_desc                                            # noqa: E402
from cores import physical_cores                               # noqa: E402

THREADS = str(physical_cores())
VENV = os.path.join(HERE, "venv", "bin", "python")
DRIVER = os.path.join(convert.MEQ, "build", "meq")
TARGETS = [1e-3, 1e-4, 1e-5, 1e-6]


def physical_cpu_list(n=None):
	"""One logical CPU per physical core, as a taskset list.

	SMT SIBLINGS SHARE AN FPU, so handing eight "threads" to four cores and
	their siblings is not eight cores' worth of FP64 -- cores.py's docstring
	makes the same point for the thread COUNT and this is the same point for
	WHICH cpus.  On this machine `lscpu -e=CPU,CORE` pairs cpu 2k and 2k+1 on
	core k, so the physical set is the even ones; the pairing is read out of
	/proc rather than assumed.
	"""
	seen, order = {}, []
	core = None
	with open("/proc/cpuinfo") as handle:
		for line in handle:
			if line.startswith("processor"):
				cpu = int(line.split(":")[1])
			elif line.startswith("core id"):
				core = int(line.split(":")[1])
			elif not line.strip() and core is not None:
				if core not in seen:
					seen[core] = cpu; order.append(cpu)
				core = None
	order = order or list(range(int(THREADS)))
	return order[:(n or int(THREADS))]


def environment():
	"""Both thread axes at the physical core count, for both arms.

	cores.py's docstring is the argument: MKL_NUM_THREADS = 1 is the ctest
	setting, pinned there because the bit-exactness assertions between assembly
	modes hold at that value, and it is not the production one.

	AND THE TWO ARMS ARE NOT LIMITED BY THE SAME MECHANISM, WHICH IS MEASURED
	RATHER THAN ASSUMED.  MEQ takes its threads from OMP_NUM_THREADS and
	MKL_NUM_THREADS and honours both.  DESC's CPU work goes through XLA's own
	Eigen thread pool, which reads NEITHER: a DESC solve left alone was
	measured at **1257% CPU**, i.e. twelve and a half cores, with
	OMP_NUM_THREADS set to 8.  So the DESC arm is pinned with `taskset` to one
	logical cpu per physical core, which is a cap the process cannot talk its
	way out of, and the MEQ arm is given the same set for the same reason.
	Pinning only one of them would be a measurement of who ignored an
	environment variable.
	"""
	return {**os.environ, "OMP_NUM_THREADS": THREADS,
	        "MKL_NUM_THREADS": THREADS, "MPLBACKEND": "Agg",
	        "JAX_PLATFORMS": "cpu"}


def pinned(argv):
	"""`argv` under taskset on the physical cores, when taskset exists."""
	cpus = ",".join(str(c) for c in physical_cpu_list())
	if os.path.exists("/usr/bin/taskset"):
		return ["/usr/bin/taskset", "-c", cpus] + list(argv)
	return list(argv)


def load_average():
	with open("/proc/loadavg") as handle:
		return float(handle.read().split()[0])


#: Where DESC's source lives, for the provenance SHA only -- the venv holds the
#: installed copy and this is never imported from.  `DESC_DIR` overrides it.
DESC_DIR = os.environ.get("DESC_DIR",
                          os.path.join(os.path.dirname(convert.MEQ), "DESC"))

#: The MFEM MEQ links.  The same default and the same environment variable
#: MEQ's own CMake finder uses -- CLAUDE.md, *Which MFEM, and why not master* --
#: so a tree configured with a different MFEM_DIR is recorded correctly here
#: without being told twice.
MFEM_DIR = os.environ.get("MFEM_DIR",
                          os.path.join(os.path.dirname(convert.MEQ), "mfem",
                                       "install"))


def provenance():
	"""What was raced, so a number can be traced back to two trees.

	**ABSOLUTE PATHS WERE BAKED IN HERE AND ARE NOT ANY MORE.** A race result
	is only worth keeping if it says which MEQ and which DESC produced it, and
	a hardcoded `/home/ian/...` makes the harness silently report `?` on any
	other machine -- which looks like a missing checkout rather than a wrong
	path.  `MFEM_DIR` is the variable MEQ's own finder reads, and `DESC_DIR`
	is offered for the same reason; both default beside the MEQ tree, which is
	where this repository's sibling checkouts live.
	"""
	def sha(path):
		if not os.path.isdir(path):
			return "?"
		out = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
		                     capture_output=True, text=True)
		return out.stdout.strip() or "?"
	lib = os.path.join(MFEM_DIR, "lib", "libmfem.a")
	return dict(meq=sha(convert.MEQ), desc=sha(DESC_DIR),
	            desc_dir=DESC_DIR, mfem_dir=MFEM_DIR,
	            threads=THREADS,
	            libmfem_mtime=(time.strftime("%Y-%m-%d %H:%M",
	                                         time.localtime(os.path.getmtime(lib)))
	                           if os.path.exists(lib) else "?"))


# ---------------------------------------------------------------- the MEQ arm
def meq_config(stem, out, degree, refine, grid=257):
	"""The shipped TOML with three keys edited and the paths absolutised.

	`RefinementLevels` IS THE h KNOB and NR/NZ are not, for the reason the
	shipped file's own comment gives: bisection of triangles is exact, so a
	level is NESTED in the last and Gamma does not move, while changing NR
	changes WHICH background elements fall inside Gamma and therefore poses a
	different problem rather than a finer one.
	"""
	text = open(os.path.join(convert.MEQ, "examples", f"{stem}.toml")).read()
	text = re.sub(r"^PolynomialDegree = .*$", f"PolynomialDegree = {degree}",
	              text, flags=re.M)
	text = re.sub(r"^RefinementLevels = .*$", f"RefinementLevels = {refine}",
	              text, flags=re.M)
	text = text.replace('"examples/', f'"{convert.MEQ}/examples/')
	text = re.sub(r"^Directory = .*$", f'Directory = "{out}"', text, flags=re.M)
	stemmed = f"{stem}-k{degree}-r{refine}"
	text = re.sub(r"^Prefix = .*$", f'Prefix = "{stemmed}"', text, flags=re.M)
	text = re.sub(r"^GridNR = .*$", f"GridNR = {grid}", text, flags=re.M)
	text = re.sub(r"^GridNZ = .*$", f"GridNZ = {grid}", text, flags=re.M)
	path = os.path.join(out, f"{stemmed}.toml")
	open(path, "w").write(text)
	return path, os.path.join(out, f"{stemmed}.nc")


def run_meq(config, timeout=7200):
	started = time.perf_counter()
	done = subprocess.run(pinned([DRIVER, config]), capture_output=True,
	                      text=True, env=environment(), timeout=timeout,
	                      cwd=convert.MEQ)
	seconds = time.perf_counter() - started
	elements = newton = 0
	found = re.search(r"converged in (\d+) Newton iterations on (\d+) elements",
	                  done.stdout)
	if found:
		newton, elements = int(found.group(1)), int(found.group(2))
	legs = re.search(r"wall ([0-9.]+) s = setup ([0-9.]+) \+ solve ([0-9.]+)"
	                 r" \+ output ([0-9.]+)", done.stdout)
	return dict(seconds=seconds, ok=done.returncode == 0, exit=done.returncode,
	            elements=elements, newton=newton,
	            wall=float(legs.group(1)) if legs else None,
	            setup=float(legs.group(2)) if legs else None,
	            solve=float(legs.group(3)) if legs else None,
	            output=float(legs.group(4)) if legs else None,
	            tail=(done.stdout or done.stderr).strip().split("\n")[-1][:160])


def meq_trace_dofs(elements, degree):
	"""AN ESTIMATE of MEQ's globally-coupled unknowns, and it says so.

	The hybridised system is solved for the trace alone: one scalar space of
	degree k on every face.  For a triangulation with E elements, 3E/2 faces is
	Euler's count in the interior and an over-count by the boundary's own
	faces, which is O( sqrt E ).  So this reads high by a per cent or two on
	the coarse meshes and less on the fine ones, and it is never the number a
	conclusion rests on -- `--mode accuracy` is.
	"""
	return int(round(1.5 * elements * (degree + 1)))


# --------------------------------------------------------------- the DESC arm
def pose_self_consistent(stem, out, M, sweeps, mix, terms=12, timeout=7200):
	"""Iterate the flux-label map to a fixed point, UNTIMED, and return the
	converged field as a file `run_desc` can be pointed at with `--from`.

	**THE SWEEPS ARE POSING AND POSING IS NOT SOLVING**, which is the line
	descrun.py's header draws and the reason MEQ's meshing sits outside
	race.py's clock too.  What is being raced is a forward solve given the
	pressure and the field as functions of the toroidal flux; arriving at those
	functions is what this does, once, and a user who had them from anywhere
	else would not pay for it.

	It is also the only way the race is between two codes rather than three.
	Posed from the reference, DESC's input carries freegs4e's own interior
	surfaces and a MEQ-against-DESC number is partly MEQ against freegs4e --
	measured, and it flatters DESC: on fixed-h-circular at M = 10 the
	reference-posed `psi_ax` is 7.402e-05 from the reference and the
	self-consistent one is 2.129e-04, because the first was partly circular.

	Returns the .npz path, or "" if the fixed point was not reached -- in
	which case the caller must NOT time a run on it and pretend otherwise.
	"""
	row_path = os.path.join(out, f"{stem}-desc-M{M}-sc.json")
	argv = [VENV, os.path.join(HERE, "descrun.py"), stem, "-M", str(M),
	        "--terms", str(terms), "--warm", "0", "--out", out,
	        "--self-consistent", str(sweeps),
	        "--self-consistent-mix", str(mix), "--json", row_path]
	done = subprocess.run(pinned(argv), capture_output=True, text=True,
	                      env=environment(), timeout=timeout)
	if done.returncode != 0 or not os.path.exists(row_path):
		return "", dict(ok=False, tail=(done.stdout or done.stderr)
		                .strip().split("\n")[-1][:200])
	row = json.load(open(row_path))
	field = row.get("npz", "")
	if not row.get("self_consistent_reached"):
		return "", dict(ok=False, history=row.get("self_consistent_history"),
		                tail="the posing did not reach its tolerance")
	return field, dict(ok=True, sweeps=len(row.get("self_consistent_history", [])),
	                   history=row.get("self_consistent_history"))


def run_desc(stem, out, M, terms=12, warm=1, grid_from="", map_from="",
             timeout=7200):
	"""A FRESH PROCESS, so the JAX compilation is inside the clock."""
	row_path = os.path.join(out, f"{stem}-desc-M{M}.json")
	argv = [VENV, os.path.join(HERE, "descrun.py"), stem, "-M", str(M),
	        "--terms", str(terms), "--warm", str(warm), "--out", out,
	        "--json", row_path]
	if grid_from:
		argv += ["--grid-from", grid_from]
	if map_from:
		argv += ["--from", map_from]
	started = time.perf_counter()
	done = subprocess.run(pinned(argv), capture_output=True, text=True,
	                      env=environment(), timeout=timeout)
	seconds = time.perf_counter() - started
	row = json.load(open(row_path)) if os.path.exists(row_path) else {}
	return dict(seconds=seconds, ok=done.returncode == 0 and bool(row),
	            exit=done.returncode, process_seconds=seconds, **row,
	            tail=(done.stdout or done.stderr).strip().split("\n")[-1][:160])


def desc_scale(stem, M):
	"""The minor radius over M: DESC's finest poloidal feature."""
	meta = json.load(open(os.path.join(convert.MEQ, "examples",
	                                   f"{stem}-meta.json")))
	return float(meta["minor_radius"]) / M


def meq_scale(stem, elements, degree):
	"""h / k, with h from the element count and the plasma's own area."""
	meta = json.load(open(os.path.join(convert.MEQ, "examples",
	                                   f"{stem}-meta.json")))
	area = np.pi * meta["minor_radius"] ** 2 * meta["kappa"]
	return float(np.sqrt(2.0 * area / max(elements, 1)) / degree)


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stems", nargs="*", default=["fixed-h-circular"])
	ap.add_argument("--mode", choices=("wall", "accuracy"), default="accuracy")
	ap.add_argument("--out", default="/tmp/race-desc")
	ap.add_argument("--degrees", default="1,2,3")
	ap.add_argument("--refinements", default="0,1")
	ap.add_argument("--modes", default="6,8,10,12,16,20",
	                help="DESC's poloidal resolution M")
	ap.add_argument("--repeats", type=int, default=1)
	ap.add_argument("--grid", type=int, default=257,
	                help="the common ( R, z ) grid both fields are written on")
	ap.add_argument("--core", type=float, default=None)
	ap.add_argument("--terms", type=int, default=12)
	ap.add_argument("--map-from", default="",
	                help="a MEQ .nc supplying the flux-label map instead of "
	                     "the freegs4e reference; convert.py's field_from_meq "
	                     "says why that is a different and better posing")
	ap.add_argument("--match", choices=("dofs", "scale"), default="dofs",
	                help="only for --mode wall; see this file's header")
	ap.add_argument("--self-consistent", type=int, default=0, metavar="N",
	                help="pose each case self-consistently first, UNTIMED, at "
	                     "most N sweeps, and then time one cold forward solve "
	                     "on the converged profiles. 0 poses from --map-from "
	                     "or the reference and leaves that code in the posing")
	ap.add_argument("--self-consistent-mix", type=float, default=0.5,
	                help="see descrun.py: 1 does not converge")
	ap.add_argument("--require-quiet", action="store_true")
	ap.add_argument("--quiet-threshold", type=float, default=0.5)
	ap.add_argument("--settle-seconds", type=float, default=20.0,
	                help="with --require-quiet, how long to wait between load "
	                     "checks before each case")
	ap.add_argument("--settle-tries", type=int, default=30)
	ap.add_argument("--json", default="")
	args = ap.parse_args()

	os.makedirs(args.out, exist_ok=True)
	before = load_average()
	if args.require_quiet and before >= args.quiet_threshold:
		print(f"load is {before:.2f}, threshold {args.quiet_threshold}: "
		      "these seconds would be a measurement of the machine, not of the "
		      "codes", file=sys.stderr)
		return 2
	# **AND THE LOAD IS RECORDED PER CASE, BECAUSE ONE CHECK AT THE TOP IS NOT
	# THE GUARANTEE THE FLAG LOOKS LIKE.**
	#
	# The check above fires once.  Measured on the first full run of this
	# harness: the five cases started at load 0.39, 3.36, 7.83, 5.29 and 7.72,
	# so only the FIRST was timed under the condition --require-quiet was asked
	# for -- the rest started while the previous case's arms were still
	# draining.  Those numbers are the race's OWN work and no other workload was
	# present, so they are not wrong; what was wrong is that nothing said so and
	# the flag implied otherwise.
	#
	# A SETTLE RATHER THAN A REFUSAL, because refusing mid-matrix would throw
	# away the cases already run, and because the load a race leaves behind is
	# its own and will drain on its own.  It waits, and it records what it
	# waited for, so a reader can see which rows were taken clean.
	def settle(tag):
		if not args.require_quiet:
			return load_average()
		for _ in range(args.settle_tries):
			now = load_average()
			if now < args.quiet_threshold:
				return now
			time.sleep(args.settle_seconds)
		now = load_average()
		print(f"    {tag}: load {now:.2f} did not fall below "
		      f"{args.quiet_threshold} in "
		      f"{args.settle_tries * args.settle_seconds} s -- THESE SECONDS "
		      f"ARE NOT A QUIET-MACHINE MEASUREMENT", file=sys.stderr)
		return now

	degrees = [int(v) for v in args.degrees.split(",") if v.strip()]
	refines = [int(v) for v in args.refinements.split(",") if v.strip()]
	modes = [int(v) for v in args.modes.split(",") if v.strip()]

	everything = dict(provenance=provenance(), load_before=before,
	                  mode=args.mode, cases={})
	for stem in args.stems:
		reference = os.path.join(convert.REFDIR, f"{convert.CASES[stem]}.npz")
		meta_path = os.path.join(convert.MEQ, "examples", f"{stem}-meta.json")
		meta = json.load(open(meta_path))
		print(f"\n  {stem}   load {settle(stem):.2f}")
		print(f"    {'code':5s} {'resolution':>12s} {'dofs':>9s} {'scale m':>9s} "
		      f"{'seconds':>9s} {'warm':>8s} {'vs ref':>10s} {'vs other':>10s}")

		meq_rows, desc_rows = [], []
		# INTERLEAVED: one MEQ rung, then one DESC rung, alternating, so that a
		# machine drifting through the run drifts through both arms.  M-102.
		meq_plan = [(k, r) for k in degrees for r in refines]
		plan = []
		for i in range(max(len(meq_plan), len(modes))):
			if i < len(meq_plan):
				plan.append(("meq",) + meq_plan[i])
			if i < len(modes):
				plan.append(("desc", modes[i]))

		for item in plan:
			if item[0] == "meq":
				_, degree, refine = item
				config, nc = meq_config(stem, args.out, degree, refine, args.grid)
				trials = [run_meq(config) for _ in range(args.repeats)]
				result = trials[0]
				result["seconds"] = statistics.median(t["seconds"] for t in trials)
				if not result["ok"]:
					print(f"    MEQ   k={degree} r={refine}  FAILED  {result['tail']}")
					continue
				dofs = meq_trace_dofs(result["elements"], degree)
				scale = meq_scale(stem, result["elements"], degree)
				row = dict(result)
				row.update(code="meq", degree=degree, refine=refine,
				           dofs=dofs, scale=scale, nc=nc)
				meq_rows.append(row)
				print(f"    MEQ   k={degree} r={refine}{'':>4s} {dofs:9d} "
				      f"{scale:9.4f} {result['seconds']:9.2f} {'-':>8s} "
				      f"{'':>10s} {'':>10s}")
			else:
				M = item[1]
				# ---- PHASE A, UNTIMED: the posing ------------------------
				# See pose_self_consistent().  A race target given the
				# profiles as functions of the toroidal flux should not be
				# charged for arriving at those functions, and the sweeps are
				# how they are arrived at.  Skipped entirely at
				# --self-consistent 0, which is the old behaviour.
				posed = args.map_from
				if args.self_consistent:
					posed, sc = pose_self_consistent(
						stem, args.out, M, args.self_consistent,
						args.self_consistent_mix, terms=args.terms)
					if not sc["ok"]:
						print(f"    DESC  M={M}  POSING FAILED  {sc['tail']}"
						      "  -- not timed")
						continue
					print(f"    DESC  M={M:<3d}  posed self-consistently in "
					      f"{sc['sweeps']} sweeps ( untimed )")
				# ---- PHASE B, TIMED: one cold forward solve --------------
				trials = [run_desc(stem, args.out, M, warm=1,
				                   grid_from=reference,
				                   map_from=posed,
				                   terms=args.terms)
				          for _ in range(args.repeats)]
				result = trials[0]
				result["seconds"] = statistics.median(t["seconds"] for t in trials)
				if not result["ok"]:
					print(f"    DESC  M={M}  FAILED  {result['tail']}")
					continue
				row = dict(result)
				row.update(code="desc", M=M, dofs=result.get("dim_x", 0),
				           scale=desc_scale(stem, M))
				desc_rows.append(row)
				print(f"    DESC  M={M:<3d}{'':>6s} {row['dofs']:9d} "
				      f"{row['scale']:9.4f} {result['seconds']:9.2f} "
				      f"{result.get('warm_seconds') or float('nan'):8.2f} "
				      f"{'':>10s} {'':>10s}")
			sys.stdout.flush()

		# ---- errors, after both arms are done, so nothing is timed while a
		# comparison runs.
		for row in desc_rows:
			row["vs_reference"] = compare_desc.grid_pair(
				row["npz"], reference, meta, args.core)["rel_l2"]
		for row in meq_rows:
			try:
				row["vs_reference"] = compare_desc.meq_pair(
					reference, row["nc"], meta_path, core=args.core)["rel_l2"]
			except Exception as why:                       # noqa: BLE001
				row["vs_reference"] = float("nan"); row["note"] = str(why)[:120]
		# the cross pairing, every MEQ rung against the FINEST DESC run and
		# every DESC rung against the finest MEQ run: neither code is the
		# truth, so each is scored against the other's best.
		if desc_rows and meq_rows:
			best_desc = min(desc_rows, key=lambda r: r["vs_reference"])
			best_meq = min(meq_rows, key=lambda r: r["vs_reference"])
			for row in meq_rows:
				row["vs_other"] = compare_desc.meq_pair(
					best_desc["npz"], row["nc"], meta_path,
					core=args.core)["rel_l2"]
			for row in desc_rows:
				row["vs_other"] = compare_desc.meq_pair(
					row["npz"], best_meq["nc"], meta_path,
					core=args.core)["rel_l2"]
			floor = compare_desc.meq_pair(best_desc["npz"], best_meq["nc"],
			                              meta_path, core=args.core)["rel_l2"]
		else:
			floor = float("nan")

		print(f"\n    {'code':5s} {'resolution':>12s} {'seconds':>9s} "
		      f"{'vs freegs4e':>12s} {'vs the other':>13s}")
		for row in sorted(meq_rows + desc_rows, key=lambda r: r["seconds"]):
			label = (f"k={row['degree']} r={row['refine']}" if row["code"] == "meq"
			         else f"M={row['M']}")
			print(f"    {row['code'].upper():5s} {label:>12s} "
			      f"{row['seconds']:9.2f} {row['vs_reference']:12.4e} "
			      f"{row.get('vs_other', float('nan')):13.4e}")
		print(f"\n    THE FLOOR: the two codes' own best answers differ by "
		      f"{floor:.4e}.")
		print("    A target below that is not resolvable by this comparison and "
		      "is refused.")

		if args.mode == "accuracy":
			print(f"\n    {'target':>9s} {'MEQ s':>9s} {'at':>12s} "
			      f"{'DESC s':>9s} {'at':>8s} {'DESC warm':>10s} {'ratio':>7s}")
			for target in TARGETS:
				if target < floor:
					print(f"    {target:9.0e}   below the floor -- refused")
					continue
				m = [r for r in meq_rows if r["vs_reference"] <= target]
				d = [r for r in desc_rows if r["vs_reference"] <= target]
				bm = min(m, key=lambda r: r["seconds"]) if m else None
				bd = min(d, key=lambda r: r["seconds"]) if d else None
				print(f"    {target:9.0e} "
				      + (f"{bm['seconds']:9.2f} {'k=%d r=%d' % (bm['degree'], bm['refine']):>12s} "
				         if bm else f"{'never':>9s} {'-':>12s} ")
				      + (f"{bd['seconds']:9.2f} {'M=%d' % bd['M']:>8s} "
				         f"{bd.get('warm_seconds') or float('nan'):10.2f} "
				         if bd else f"{'never':>9s} {'-':>8s} {'-':>10s} ")
				      + (f"{bd['seconds']/bm['seconds']:7.2f}" if bm and bd else ""))
		else:
			key = "dofs" if args.match == "dofs" else "scale"
			print(f"\n    MATCHED ON {key.upper()} -- see this file's header for "
			      "why this is a choice and not a fact")
			print(f"    {'MEQ':>16s} {'DESC':>16s} {'MEQ s':>8s} {'DESC s':>8s} "
			      f"{'warm':>8s} {'ratio':>7s}")
			for row in meq_rows:
				if not desc_rows:
					break
				near = min(desc_rows, key=lambda r: abs(np.log(max(r[key], 1e-30))
				                                        - np.log(max(row[key], 1e-30))))
				print(f"    k={row['degree']} r={row['refine']} {row[key]:8.4g} "
				      f"   M={near['M']:<3d} {near[key]:8.4g} "
				      f"{row['seconds']:8.2f} {near['seconds']:8.2f} "
				      f"{near.get('warm_seconds') or float('nan'):8.2f} "
				      f"{near['seconds']/row['seconds']:7.2f}")

		everything["cases"][stem] = dict(meq=meq_rows, desc=desc_rows,
		                                 floor=floor)

	everything["load_after"] = load_average()
	path = args.json or os.path.join(args.out, "race-desc.json")
	json.dump(everything, open(path, "w"), indent=2, default=str)
	print(f"\n  wrote {path}")
	print("  THESE SECONDS ARE ONLY A MEASUREMENT ON AN IDLE MACHINE.  CLAUDE.md "
	      "records a\n  490 s / 540 s spread on identical code; re-run before "
	      "reading anything into a\n  change of tens of per cent.")
	return 0


if __name__ == "__main__":
	sys.exit(main())
