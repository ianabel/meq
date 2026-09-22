"""Cost of a given accuracy: MEQ against the freegs4e reference it reproduces.

THE TWO CODES DO NOT SOLVE THE SAME PROBLEM, so a wall-clock ratio is not a
statement about either.  freegs4e converges a FREE-boundary equilibrium -- coil
Green's functions, X-point finding, a control system and 23 to 103 Picard steps
-- while MEQ solves the FIXED-boundary problem inside a surface it is handed.
MEQ is doing strictly less work, and it is in C++ against Python.  A raw ratio
flatters MEQ for reasons that have nothing to do with the discretisation.

SO THE COLUMN THAT MEANS SOMETHING IS ACCURACY PER UNKNOWN, which is
language-neutral and is the real question: does a high-order hybridized DG reach
a given accuracy with fewer degrees of freedom than a 2nd/4th-order finite
difference grid?  Wall clock is reported beside it and should be read only as an
order of magnitude, per this project's standing rule that a timing here is a
measurement about this machine.
"""
import json
import os
import re
import subprocess
import time

import numpy as np

from compare import compare

MEQ = "/home/ian/projects/meq/build/meq"


def run_one(stem, outdir, degree, refine, grid=129, env=None):
    """One MEQ solve at a given degree and refinement; returns its cost and error."""
    base = os.path.join(outdir, f"{stem}.toml")
    text = open(base).read()
    text = re.sub(r"^PolynomialDegree = .*$", f"PolynomialDegree = {degree}",
                  text, flags=re.M)
    text = re.sub(r"^RefinementLevels = .*$", f"RefinementLevels = {refine}",
                  text, flags=re.M)
    text = re.sub(r"^GridNR = .*$", f"GridNR = {grid}", text, flags=re.M)
    text = re.sub(r"^GridNZ = .*$", f"GridNZ = {grid}", text, flags=re.M)
    text = text.replace(f'Prefix = "{stem}"', 'Prefix = "perf"')
    path = os.path.join(outdir, "perf.toml")
    open(path, "w").write(text)

    e = dict(os.environ)
    e["MKL_NUM_THREADS"] = "1"
    e["OMP_NUM_THREADS"] = "1"
    if env:
        e.update(env)

    t0 = time.perf_counter()
    proc = subprocess.run([MEQ, path], capture_output=True, text=True,
                          cwd=outdir, env=e, timeout=1800)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        return dict(ok=False, wall=wall, log=proc.stdout + proc.stderr)

    m = re.search(r"converged in (\d+) Newton iterations on (\d+) elements",
                  proc.stdout)
    iters = int(m.group(1)) if m else -1
    elements = int(m.group(2)) if m else -1

    # P_k on a triangle has ( k+1 )( k+2 )/2 dofs; the potential is one such
    # field and the flux is two. The trace is what is actually solved after
    # hybridization, but the total is the honest measure of the discrete
    # problem's size against a grid point count.
    per = (degree + 1) * (degree + 2) // 2
    dofs = elements * per * 3

    result = dict(ok=True, wall=wall, iters=iters, elements=elements, dofs=dofs)
    if grid > 17:
        comparison = compare(os.path.join(outdir, f"{stem}.npz"),
                    os.path.join(outdir, "perf.nc"),
                    os.path.join(outdir, f"{stem}-meta.json"))
        result["rel_l2"] = comparison["rel_l2"]
    return result


if __name__ == "__main__":
    import sys
    outdir = sys.argv[1]
    cases = sys.argv[2:]

    print("\n  COST OF AN ACCURACY: MEQ AGAINST ITS freegs4e REFERENCE\n")
    print("  MKL_NUM_THREADS=1, OMP_NUM_THREADS=1. Timings include mesh")
    print("  generation, the solve, post-processing and all three output")
    print("  formats -- i.e. the whole driver, not the solve alone.\n")

    for stem in cases:
        ref = json.load(open(os.path.join(outdir, f"{stem}.json")))
        refwall = ref["wall_time_s"]
        print(f"  {stem}")
        print(f"    freegs4e: {refwall:.1f} s, 129^2 = 16641 grid points, "
              f"{ref.get('stage2_iterations', 0)} Picard steps, free boundary\n")
        print(f"    {'k':>2} {'ref':>4} {'elements':>9} {'dofs':>9} "
              f"{'newton':>7} {'wall/s':>8} {'rel L2':>11} {'dofs/ref pt':>12}")
        for degree, refine in [(1, 2), (1, 3), (2, 1), (2, 2), (2, 3),
                               (3, 1), (3, 2)]:
            run = run_one(stem, outdir, degree, refine)
            if not run["ok"]:
                print(f"    {degree:>2} {refine:>4}   FAILED")
                continue
            print(f"    {degree:>2} {refine:>4} {run['elements']:>9} "
                  f"{run['dofs']:>9} {run['iters']:>7} {run['wall']:>8.2f} "
                  f"{run['rel_l2']:>11.3e} {run['dofs']/16641:>12.2f}")
        print()
