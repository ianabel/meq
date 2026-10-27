#!/usr/bin/env python3
"""The device timing kit's results as MEASUREMENTS.md M-242 prints them.

    python3 summarize.py results/            tables on stdout
    python3 summarize.py results/ --csv      every run, one row each

For each case and arm, the best round by wall (the M-135 convention: best of
the rounds, interleaved), with setup, solve and output, the Newton count,
solve/N and psi_ax; then the solve's main legs for that round, and the GPU's
solve against each host arm. psi_ax is checked across the arms of a case, and
a disagreement is flagged rather than averaged.

Standard library only, so it runs on any Python 3 the cluster has.
"""

import csv
import re
import sys
from pathlib import Path

WALL = re.compile(r"^MEQ: wall ([0-9.]+) s = setup ([0-9.]+) \+ solve ([0-9.]+) \+ output ([0-9.]+)")
NEWTON = re.compile(r"converged in (\d+) Newton iterations on (\d+) elements")
PSI_AX = re.compile(r"psi_ax = ([0-9.eE+-]+)")
LEG = re.compile(r"^MEQ:   (\S[^\n]*?)\s{2,}([0-9.]+)\s+([0-9.]+)%\s+([0-9.]+|-)\s+(\d+)\s*$")
LEGS = ["residual", "gradient", "trace factorisation", "trace backsolve",
        "NPC reduce+recover", "constraint location", "border assembly",
        "re-assembly", "outside solve()"]


def parse(log):
    run = {"legs": {}}
    for line in log.read_text(errors="replace").splitlines():
        if m := WALL.match(line):
            run.update(wall=float(m[1]), setup=float(m[2]), solve=float(m[3]), output=float(m[4]))
        elif (m := NEWTON.search(line)) and "newton" not in run:
            run.update(newton=int(m[1]), elements=int(m[2]))
        elif (m := PSI_AX.search(line)) and "psi_ax" not in run:
            run["psi_ax"] = m[1]
        elif (m := LEG.match(line)) and not m[1].startswith("of which"):
            cores = None if m[4] == "-" else float(m[4])
            run["legs"][m[1].strip()] = (float(m[2]), cores, int(m[5]))
    return run


def collect(root):
    runs = []
    for log in sorted(root.glob("*/*/r*/run.log")):
        case, arm, rnd = log.parts[-4], log.parts[-3], log.parts[-2]
        code = (log.parent / "exit").read_text().strip() if (log.parent / "exit").exists() else "?"
        runs.append({"case": case, "arm": arm, "round": rnd, "exit": code, **parse(log)})
    return runs


def arm_key(arm):
    return (arm == "cuda", int(arm[3:]) if arm.startswith("cpu") and arm[3:].isdigit() else 0)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    root = Path(sys.argv[1])
    runs = collect(root)
    if not runs:
        sys.exit(f"no runs under {root}")

    if "--csv" in sys.argv:
        out = csv.writer(sys.stdout)
        out.writerow(["case", "arm", "round", "exit", "elements", "newton", "wall", "setup",
                      "solve", "output", "psi_ax"] + LEGS)
        for r in runs:
            out.writerow([r["case"], r["arm"], r["round"], r["exit"], r.get("elements"),
                          r.get("newton"), r.get("wall"), r.get("setup"), r.get("solve"),
                          r.get("output"), r.get("psi_ax")]
                         + [r["legs"].get(leg, (None,))[0] for leg in LEGS])
        return

    if (env := root / "env.txt").exists():
        print("\n".join(env.read_text().splitlines()[:5]) + "\n")
    for case in sorted({r["case"] for r in runs}, key=lambda c: (c != "mf", c)):
        mine = [r for r in runs if r["case"] == case]
        arms = sorted({r["arm"] for r in mine}, key=arm_key)
        best = {}
        failed = []
        for arm in arms:
            ok = [r for r in mine if r["arm"] == arm and r["exit"] == "0" and "wall" in r]
            failed += [f'{arm} {r["round"]} exit {r["exit"]}' for r in mine
                       if r["arm"] == arm and r not in ok]
            if ok:
                best[arm] = min(ok, key=lambda r: r["wall"]) | {"rounds": len(ok)}
        elements = next((r.get("elements") for r in best.values()), "?")
        print(f"## {case}, {elements} elements")
        print("| arm | rounds | wall | setup | solve | output | N | solve/N | psi_ax |")
        print("|---|---|---|---|---|---|---|---|---|")
        for arm, r in best.items():
            n = r.get("newton") or 0
            per = f'{r["solve"] / n:.3f}' if n else "-"
            print(f'| {arm} | {r["rounds"]} | {r["wall"]:.2f} | {r["setup"]:.2f} | {r["solve"]:.2f} '
                  f'| {r["output"]:.2f} | {n} | {per} | {r.get("psi_ax", "-")} |')
        labels = {r.get("psi_ax") for r in best.values()}
        if len(labels) > 1:
            print(f"\n**psi_ax DISAGREES between arms: {sorted(labels)}** -- not one problem; "
                  "read check.txt before any timing in this case.")
        if failed:
            print(f"\nfailed or incomplete: {', '.join(failed)}")
        if "cuda" in best:
            gpu = best["cuda"]
            ratios = [f'{gpu["solve"] / r["solve"]:.2f}x of {arm}' for arm, r in best.items()
                      if arm != "cuda" and r["solve"] > 0]
            if ratios:
                print(f'\nGPU solve: {", ".join(ratios)}')
        print("\n| leg, s (cores) | " + " | ".join(best) + " |")
        print("|---|" + "---|" * len(best))
        for leg in LEGS:
            cells = []
            for r in best.values():
                sec, cores, _ = r["legs"].get(leg, (None, None, 0))
                cells.append("-" if sec is None else f"{sec:.3f}" + (f" ({cores:.2f})" if cores else ""))
            print(f"| {leg} | " + " | ".join(cells) + " |")
        print()


if __name__ == "__main__":
    main()
