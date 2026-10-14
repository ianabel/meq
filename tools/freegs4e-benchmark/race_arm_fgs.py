"""THE INVERSE RACE'S freegs4e ARM. RACE-ARMS.md is the contract it implements.

    race_arm_fgs.py --problem P.json --config NAME --rung R
                    --start cold|seeded|warm [--state DIR] --out DIR
    race_arm_fgs.py --problem P.json --config NAME --list-rungs

THE SAME SOLVE fgsref.py RUNS FOR ITS REFERENCES, WITHOUT ITS ROBUSTNESS
LADDER: the shaped machine through shaped.py ( turns folded, every conductor on
a tensor Gauss rule, M-255 ), the cached boundary matrix ( fgsref.py's
WORKAROUND 6 ), GeneralPprimeFFprime with INTERPOLATING splines so the profile
is the analytic one every arm solves, Ip_logic closing I_p with one scale, and
SyncConstrain at gamma = 1e-12 on the problem's targets. What is left out is
fgsref.py's retry with blending and its ConstrainPaxisIp reseed: those are a
robustness device no other arm is given, and a run that needs them is a DNF.

ONLY THE FREE SUPPLIES ARE UNDER CONTROL. freegs4e's `control` flag is per
top-level item, which is exactly the problem's supply; every other supply sits
at the reference's current.

THE STARTS.
  cold    freegs4e's ordinary plasma guess ( Equilibrium's own initial
          plasma_psi -- what fgsref.py starts from ), free supplies at the
          vacuum seed, the design targets
  seeded  the same, free supplies at the reference's currents
  warm    this arm's own converged plasma_psi and currents at the same grid,
          and the SHIFTED targets

ROWS ARE eq.Br / eq.Bz / eq.psiRZ ON THE CONVERGED FIELD, which sum over the
coil objects analytically -- not the Picard loop's step size.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

import numpy as np                                               # noqa: E402

import benchpaths                                                # noqa: E402
import fgsref                                                    # noqa: E402
import shaped                                                    # noqa: E402
import conductors                                                # noqa: E402

LADDER = [65, 129, 257, 513]
fgsref.INTERPOLATING_PROFILES = True


def rungs():
    return [dict(rung=i, label="%d^2" % n, n=n, finest=(i == len(LADDER) - 1))
            for i, n in enumerate(LADDER)]


def exponents(shape):
    """( alpha, beta ) out of fgsref.py's 'P0*(1-psi_n**a)**b'."""
    m = re.search(r"psi_n\*\*([0-9.eE+-]+)\)\*\*([0-9.eE+-]+)", shape)
    return float(m.group(1)), float(m.group(2))


def provenance():
    def sha(path):
        r = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True)
        return r.stdout.strip() or None
    import freegs4e
    return dict(code_version="freegs4e %s" % getattr(freegs4e, "__version__", "?"),
                tree=sha(os.path.dirname(os.path.dirname(freegs4e.__file__))),
                meq_tree=sha(benchpaths.SOURCE),
                notes="fgsref.py's solve; shaped.py turns fold and "
                      "rectangle_rule; interpolating profile splines; cached "
                      "boundary matrix; no blend or reseed retries")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--rung", type=int)
    ap.add_argument("--start", choices=("cold", "seeded", "warm"))
    ap.add_argument("--state")
    ap.add_argument("--out")
    ap.add_argument("--list-rungs", action="store_true")
    a = ap.parse_args()
    if a.list_rungs:
        print(json.dumps([dict(rung=r["rung"], label=r["label"],
                               finest=r["finest"]) for r in rungs()]))
        return
    t0 = time.time()
    problem = {p["config"]: p for p in json.load(open(a.problem))}[a.config]
    rung = rungs()[a.rung]
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    case = dict({c["name"]: c for c in fgsref.CASES}[problem["fgsref_case"]])
    case["grid"] = dict(case["grid"], nx=rung["n"], ny=rung["n"])
    free = [s["name"] for s in problem["supplies"] if s["free"]]
    currents = {s["name"]: s["reference_current"] for s in problem["supplies"]}
    state = None
    if a.start == "warm":
        state = np.load(os.path.join(a.state, "state.npz"), allow_pickle=True)
        currents.update(json.load(open(os.path.join(a.state, "state.json")))
                        ["currents"])
        xpts, iso = (problem["shifted_xpoint_targets"],
                     problem["shifted_isoflux_targets"])
    else:
        if a.start == "cold":
            currents.update(problem["vacuum_seed"])
        xpts, iso = problem["xpoint_targets"], problem["isoflux_targets"]

    tok, _ = shaped.shaped_machine(case["make"](), half=conductors.DEFAULT_HALF)
    for label, item in tok.coils:
        item.current = float(currents[label])
        item.control = label in free
    tok.getCurrentsVec()

    eq = fgsref.build_eq(case, tok)
    eq._applyBoundary = fgsref._free_boundary_cached
    if state is not None:
        eq.plasma_psi = np.asarray(state["plasma_psi"], float)
        eq._updatePlasmaPsi(eq.plasma_psi)

    psi_n = np.linspace(0.0, 1.0, fgsref.NPSI)
    pa, pb = exponents(problem["profile_shape"]["pprime"])
    fa, fb = exponents(problem["profile_shape"]["ffprime"])
    prof = fgsref.GeneralProfile(
        Ip=problem["Ip"], fvac=problem["fvac"], psi_n=psi_n,
        pprime_data=problem["pprime_amplitude"]*fgsref.pshape(psi_n, pa, pb),
        ffprime_data=problem["ffprime_amplitude"]*fgsref.pshape(psi_n, fa, fb),
        Raxis=1.0, Ip_logic=True)
    fgsref.attach(prof, case, eq)
    ctrl = fgsref.SyncConstrain(xpoints=[tuple(p) for p in xpts],
                                isoflux=[tuple(p) for p in iso], gamma=1e-12)

    ts = time.time()
    status = "exception"
    try:
        st = fgsref.picard_loop(eq, prof, ctrl, rtol=fgsref.RTOL, atol=1e-14,
                                maxits=fgsref.MAXITS)
        status = st["status"]
        iterations = int(st["iterations"])
    except Exception as exc:                                     # noqa: BLE001
        status = "exception: %r" % (exc,)
        iterations = None
    solve = time.time() - ts

    record = dict(code="freegs4e", config=a.config, rung=a.rung,
                  rung_label=rung["label"], start=a.start,
                  unknowns=rung["n"]**2, finished=False, status=status,
                  wall_s=None, solve_s=solve, currents={}, psi_axis=None,
                  psi_bndry=None, axis=None, xpoints=[], rows=[],
                  iterations=iterations, state=None, provenance=provenance())
    if status in ("rtol", "atol"):
        from freegs4e import critical
        psi = eq.psi()
        opt, xpt = critical.find_critical(eq.R, eq.Z, psi)
        rows = []
        for R, Z in xpts:
            rows += [float(np.asarray(eq.Br(R, Z))), float(np.asarray(eq.Bz(R, Z)))]
        for R1, Z1, R2, Z2 in iso:
            rows.append(float(np.asarray(eq.psiRZ(R1, Z1))
                              - np.asarray(eq.psiRZ(R2, Z2))))
        record.update(
            finished=bool(len(opt)), status="converged" if len(opt) else "no axis",
            currents={label: float(item.current) for label, item in tok.coils
                      if label in free},
            psi_axis=float(opt[0][2]) if len(opt) else None,
            psi_bndry=float(getattr(prof, "psi_bndry", np.nan)),
            axis=[float(opt[0][0]), float(opt[0][1])] if len(opt) else None,
            xpoints=[[float(x[0]), float(x[1])] for x in xpt[:4]],
            rows=rows, state=out)
        np.savez(os.path.join(out, "state.npz"), plasma_psi=eq.plasma_psi)
        json.dump(dict(currents={label: float(item.current)
                                 for label, item in tok.coils}),
                  open(os.path.join(out, "state.json"), "w"), indent=1)
    record["wall_s"] = time.time() - t0
    json.dump(record, open(os.path.join(out, "record.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
