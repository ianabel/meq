"""THE INVERSE RACE'S PROBLEM FILE: one posed problem, read by every arm.

`RACE-ARMS.md` is the contract and `INVERSE-RACE-PLAN.md` the design. This
writes the list of configurations every arm reads, from a set of `freegs4e`
references -- runs/freegs4e-benchmark/ref-interp/ by default, `fgsref.py
--shaped --interpolating-profiles`, so every code solves the analytic profile
MEQ's tables carry (MEASUREMENTS.md M-254, M-255).

WHAT EACH CONFIGURATION FIXES, AND WHY.

  * THE FREE SUPPLIES ARE THE DETERMINED SUBSET. freegs4e frees every supply
    and its targets do not determine them, so all-free currents differ between
    codes by a selection and not an error. A column-pivoted QR of the
    row-normalised vacuum target response picks as many as its rank
    (inverse_compare.free_subset); every other supply is held at the
    reference's value.
  * THE COLD START'S FREE CURRENTS ARE THE VACUUM STAGE: the minimum-norm
    currents meeting the targets against the design blob AND the held supplies'
    own field -- freegs4e's first control step, restricted to the free ones.
    Every code can be handed them; nothing in them was solved for.
  * THE WARM START'S TARGETS ARE THE DESIGN ONES MOVED 0.01 m, each point along
    the unit vector from the reference's magnetic axis to it: a shape step, the
    workload of a scan or a controller.

USAGE

    venv/bin/python inverse_problem.py OUT.json [--reference DIR] [C D ...]

OUT.json is conventionally runs/freegs4e-benchmark/ref-interp/problem.json,
which race_inverse.py reads by default.

The default cases are the five INVERSE-RACE-PLAN.md section 2 races.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

import benchpaths                                              # noqa: E402
from fgsref import CASES                                       # noqa: E402
from circuits import MACHINES, wiring                          # noqa: E402
from inverse import target_response                            # noqa: E402
from inverse_compare import free_subset, supply_currents       # noqa: E402

RACE = ["C", "D", "E", "F", "G"]
SHIFT = 0.01


def shifted(point, axis):
    """`point` moved SHIFT metres along the unit vector from `axis` to it."""
    d = np.asarray(point, float) - np.asarray(axis, float)
    return (np.asarray(point, float) + SHIFT*d/np.linalg.norm(d)).tolist()


def up_down_symmetric(xpoints, isoflux):
    """The target set maps onto itself under Z -> -Z."""
    def key(p):
        return (round(p[0], 9), round(abs(p[1]), 9), np.sign(round(p[1], 9)))
    pts = {(round(R, 9), round(Z, 9)) for R, Z in xpoints}
    if any((R, round(-Z, 9)) not in pts for R, Z in pts):
        return False
    pairs = {frozenset([(round(a, 9), round(b, 9)), (round(c, 9), round(d, 9))])
             for a, b, c, d in isoflux}
    for pair in pairs:
        mirrored = frozenset((R, round(-Z, 9)) for R, Z in pair)
        if mirrored not in pairs:
            return False
    return True


def configuration(case, stem, refdir, wires):
    npz_name = case["name"] + "_shaped"
    npz_path = os.path.join(refdir, npz_name + ".npz")
    json_path = os.path.join(refdir, npz_name + ".json")
    d = np.load(npz_path, allow_pickle=True)
    j = json.load(open(json_path))
    A, supplies, blob = target_response(npz_path, case, wires)
    current = supply_currents(json_path, supplies)
    pick, rank, _ = free_subset(A)
    names = [n for n, _ in supplies]
    free = [names[i] for i in pick]
    held = [i for i in range(len(names)) if i not in pick]

    # The vacuum stage over the free supplies, against the blob and the held
    # supplies' own field.
    rhs = -(blob + A[:, held] @ np.array([current[names[i]] for i in held]))
    x = np.linalg.lstsq(A[:, pick], rhs, rcond=None)[0]
    vacuum = {names[i]: float(v) for i, v in zip(pick, x)}

    xpts = [list(p) for p in case["xpoints"]]
    iso = [list(p) for p in case["isoflux"]]
    symmetric = up_down_symmetric(xpts, iso)
    # Shifted from the axis ON THE MIDPLANE when the targets are symmetric:
    # the reference's axis sits ~1e-4 m off it, and shifting from there would
    # break a symmetry three of these configurations' row sets need to stay
    # consistent ( rank below the row count: C, E, G ).
    axis = (float(j["Raxis"]), 0.0 if symmetric else float(j["Zaxis"]))
    rows = 2*len(xpts) + len(iso)
    return dict(
        config="machine-%s-shaped" % stem[len("machine-"):],
        fgsref_case=case["name"], shaped=True,
        reference_npz=npz_path, reference_json=json_path,
        meq_example=os.path.join(benchpaths.EXAMPLES,
                                 "machine-%s-shaped.toml" % stem[len("machine-"):]),
        supplies=[dict(name=n, members=[dict(conductor=l, multiplier=m)
                                        for l, m in mem],
                       reference_current=current[n], free=(n in free))
                  for n, mem in supplies],
        conductors=j["conductors"],
        rank=rank, rows=rows,
        xpoint_targets=xpts, isoflux_targets=iso,
        shifted_xpoint_targets=[shifted(p, axis) for p in xpts],
        shifted_isoflux_targets=[shifted(p[:2], axis) + shifted(p[2:], axis)
                                 for p in iso],
        symmetric=symmetric,
        vacuum_seed=vacuum,
        Ip=float(j["Ip"]), profile_shape=j["profile_shape"],
        pprime_amplitude=float(np.asarray(d["pprime_tabulated"])[0]),
        ffprime_amplitude=float(np.asarray(d["ffprime_tabulated"])[0]),
        fvac=float(j["fvac"]), grid=j["grid"], design=j["design"],
        reference=dict(psi_axis=float(j["psi_axis"]),
                       psi_bndry=float(j["psi_bndry"]),
                       Raxis=axis[0], Zaxis=axis[1], xpoints=j["xpoints"]))


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out = sys.argv[1]
    args = sys.argv[2:]
    refdir = os.path.join(benchpaths.OUT, "ref-interp")
    if "--reference" in args:
        at = args.index("--reference")
        refdir = os.path.abspath(args[at + 1])
        del args[at:at + 2]
    want = args or RACE
    cases = {c["name"]: c for c in CASES}
    problem = []
    for ref, stem in MACHINES:
        if ref[0] not in want:
            continue
        case = cases[ref]
        problem.append(configuration(case, stem, refdir,
                                     wiring(case["make"]())))
        p = problem[-1]
        print("%-28s free %s  symmetric %s  vacuum seed %s" % (
            p["config"], [s["name"] for s in p["supplies"] if s["free"]],
            p["symmetric"], {k: "%.3g" % v for k, v in p["vacuum_seed"].items()}),
            flush=True)
    json.dump(problem, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
