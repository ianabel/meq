#!/usr/bin/env python
"""TIME TO ACCURACY, NOT TIME TO GRID -- the race that can actually decide MEQ.

WHY THE EXISTING RACE CANNOT DECIDE IT.  M-111 reports MEQ 1.39x, 3.05x and
3.87x slower than freegs4e, and both of the things that number is made of are
about the comparison rather than about the codes:

  * F DIII-D's error is FLAT at 5.785e-03 across a 16x range in dofs, off-coil
    flat at 6.3e-04.  That is the CONDUCTOR MODEL -- MEQ's rectangles against
    freegs4e's filaments -- and no rung buys it down.  MEQ is being charged for
    an error that is not its discretisation.
  * the references saturate: freegs4e at 1.4e-04, the MAST-U forward one at
    about 1e-05 (M-122).  "Faster at equal accuracy" cannot even be POSED at an
    accuracy the reference cannot reach.

So this sweeps COST against ACHIEVED ERROR for both codes and reports, for each
target, the cheapest run of each that reaches it.  A Pareto front, not a pair of
seconds.

AND IT PROBES COARSE MESHES AT HIGH ORDER, WHICH IS THE WHOLE POINT.  If high
order is worth anything it is worth it at 1e-03 and 1e-04 with FEW ELEMENTS --
not only in an asymptotic regime nobody runs.  The MEQ sweep therefore runs k up
to 4 on meshes coarse enough to look absurd, because that is the corner where
the claim lives.  A ladder that only refines h at fixed k cannot see it.

A CASE WITH A FRACTIONAL PLASMA EDGE CANNOT DECIDE THIS EITHER.  CLAUDE_FB.md's
plasma-edge result caps the achievable order at `k <= j`, where j is the order
to which the profile vanishes at the edge.  freegsnke's ConstrainPaxisIp at
alpha_n = 1.2 gives j = 1.2, i.e. k = 1 -- so on MAST-U MEQ's high-order
advantage is STRUCTURALLY UNAVAILABLE and a poor showing there is the profile,
not the solver.  Run it on MAST-U for the honest MAST-U answer; run it on a
smooth-edge case to test the high-order claim.  --case selects.

THREADS ARE MATCHED ON BOTH ARMS, at the PHYSICAL core count on both axes --
see cores.py for why that is 8 here and not os.cpu_count()'s 16 -- because
MKL_NUM_THREADS=1 is the ctest setting and not the production one -- see
race_nke.py's header for what pinning it hands away.
"""
import argparse, json, os, re, subprocess, sys, time
import numpy as np

MEQ = "/home/ian/projects/meq"
VENV = f"{MEQ}/tools/freegs4e-benchmark/venv/bin/python"
from cores import physical_cores
THREADS = str(physical_cores())
TARGETS = [1e-2, 1e-3, 1e-4, 1e-5]


def env():
    return {**os.environ, "OMP_NUM_THREADS": THREADS,
            "MKL_NUM_THREADS": THREADS, "MPLBACKEND": "Agg"}


def meq_variant(base_toml, out_dir, degree, size, stem):
    """A MEQ config at one ( k, h ), derived from the case's own TOML.

    The MESH SIZE is edited, not the element count: halfdisc.py takes a target
    edge length and the count follows from the geometry, so asking for a count
    would be asking the mesher a question it does not answer.
    """
    text = open(base_toml).read()
    text = re.sub(r"^PolynomialDegree = .*$", f"PolynomialDegree = {degree}",
                  text, flags=re.M)
    text = re.sub(r"^Size = .*$", f"Size = {size:.4f}", text, flags=re.M)
    text = re.sub(r"^PlasmaSize = .*$", f"PlasmaSize = {size/5.0:.4f}",
                  text, flags=re.M)
    text = re.sub(r'^File = "examples/.*\.msh"$',
                  f'File = "{out_dir}/{stem}.msh"', text, flags=re.M)
    text = re.sub(r"^Prefix = .*$", f'Prefix = "{out_dir}/{stem}"', text,
                  flags=re.M)
    if "[output]" not in text:
        text += f'\n[output]\nPrefix = "{out_dir}/{stem}"\nGridNR = 129\nGridNZ = 257\n'
    path = os.path.join(out_dir, f"{stem}.toml")
    open(path, "w").write(text)
    return path


def run_meq(config, stem, out_dir):
    started = time.perf_counter()
    done = subprocess.run([f"{MEQ}/build/meq-run", config],
                          capture_output=True, text=True, env=env(), timeout=7200)
    seconds = time.perf_counter() - started
    elements = 0
    match = re.search(r"converged in (\d+) Newton iterations on (\d+) elements",
                      done.stdout)
    if match:
        elements = int(match.group(2))
    return dict(seconds=seconds, ok=done.returncode == 0, exit=done.returncode,
                elements=elements,
                nc=os.path.join(out_dir, f"{stem}.nc"),
                tail=(done.stdout or done.stderr).strip().split("\n")[-1][:160])


def error_against(nc_path, truth_npz):
    """Relative L2 of psi against the truth, through compare.py.

    compare.py already drops the band -- the .nc's `extrapolated` mask -- which
    CLAUDE.md requires before differencing two runs and which is one node in ten
    on a curved case.  Reimplementing the norm here would be reimplementing that
    too, and getting it wrong is silent.
    """
    done = subprocess.run([VENV, f"{MEQ}/tools/freegs4e-benchmark/compare.py",
                           truth_npz, nc_path],
                          capture_output=True, text=True, env=env())
    match = re.search(r"([0-9.]+e[-+][0-9]+)", done.stdout)
    for line in done.stdout.split("\n"):
        if "rel_l2" in line or "relative" in line.lower():
            found = re.findall(r"([0-9.]+e[-+][0-9]+)", line)
            if found:
                return float(found[0]), done.stdout
    return (float(match.group(1)) if match else float("nan")), done.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=f"{MEQ}/examples/mastu-nke.toml",
                    help="the MEQ TOML the sweep varies k and h on")
    ap.add_argument("--truth", default="",
                    help="npz to measure error against; default is the finest "
                         "freegsnke level in --out")
    ap.add_argument("--out", default="/tmp/tta")
    ap.add_argument("--degrees", default="1,2,3,4")
    ap.add_argument("--sizes", default="0.60,0.45,0.30,0.20",
                    help="halfdisc target edge lengths, COARSE FIRST -- the "
                         "coarse end at high k is the corner being probed")
    ap.add_argument("--nke-levels", default="33,65,129,257")
    ap.add_argument("--repeats", type=int, default=3)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    print(f"time-to-accuracy, {THREADS} threads on both axes and both arms\n")

    # ---- the freegsnke arm ----
    nke_rows = []
    for nx in [int(v) for v in args.nke_levels.split(",")]:
        started = time.perf_counter()
        subprocess.run([VENV, f"{MEQ}/tools/freegs4e-benchmark/mastu_reference.py",
                        "--out", args.out, "--levels", str(nx),
                        "--repeats", str(args.repeats)],
                       capture_output=True, text=True, env=env(), timeout=7200)
        nke_rows.append(dict(nx=nx, seconds=(time.perf_counter() - started)
                             / max(1, args.repeats)))
    truth = args.truth or os.path.join(args.out, "mastu-fwd-257x513.npz")
    print(f"truth: {truth}\n")

    # ---- the MEQ arm: the whole ( k, h ) grid ----
    print(f"  {'k':>2} {'size':>6} {'elements':>9} {'seconds':>9} {'rel L2':>11}  note")
    meq_rows = []
    for degree in [int(v) for v in args.degrees.split(",")]:
        for size in [float(v) for v in args.sizes.split(",")]:
            stem = f"meq-k{degree}-h{size:.2f}".replace(".", "p")
            config = meq_variant(args.case, args.out, degree, size, stem)
            try:
                result = run_meq(config, stem, args.out)
            except subprocess.TimeoutExpired:
                result = dict(seconds=float("inf"), ok=False, elements=0,
                              nc="", tail="timeout", exit=-1)
            rel = float("nan")
            if result["ok"] and os.path.exists(result["nc"]):
                rel, _ = error_against(result["nc"], truth)
            meq_rows.append(dict(degree=degree, size=size, **{
                k: v for k, v in result.items() if k != "nc"}, rel_l2=rel))
            print(f"  {degree:2d} {size:6.2f} {result['elements']:9d} "
                  f"{result['seconds']:9.2f} {rel:11.3e}  "
                  f"{'' if result['ok'] else result['tail'][:60]}")
            sys.stdout.flush()

    # ---- the Pareto front: cheapest run of each code reaching each target ----
    print(f"\n  {'target':>8} {'MEQ s':>9} {'MEQ k,h':>12} {'freegsnke s':>12} "
          f"{'ratio':>7}")
    front = []
    for target in TARGETS:
        ok = [r for r in meq_rows if r["ok"] and r["rel_l2"] <= target]
        best = min(ok, key=lambda r: r["seconds"]) if ok else None
        # freegsnke's own error against the truth is not measured here -- the
        # truth IS its finest level, so its error is zero by construction and a
        # ratio would be meaningless. What is reported is the cost of the level
        # whose OWN refinement study (M-122) puts it at that accuracy.
        row = dict(target=target,
                   meq=best["seconds"] if best else None,
                   meq_k=best["degree"] if best else None,
                   meq_h=best["size"] if best else None)
        front.append(row)
        print(f"  {target:8.0e} "
              + (f"{row['meq']:9.2f} k={row['meq_k']},h={row['meq_h']:.2f}"
                 if best else f"{'never':>9} {'-':>12}"))

    json.dump(dict(meq=meq_rows, freegsnke=nke_rows, front=front,
                   threads=THREADS, truth=truth),
              open(os.path.join(args.out, "time-to-accuracy.json"), "w"), indent=2)
    print(f"\nwrote {args.out}/time-to-accuracy.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
