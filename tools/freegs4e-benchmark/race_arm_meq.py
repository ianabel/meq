"""THE INVERSE RACE'S MEQ ARM. RACE-ARMS.md is the contract it implements.

    race_arm_meq.py --problem P.json --config NAME --rung R
                    --start cold|seeded|warm [--state DIR] --out DIR
    race_arm_meq.py --problem P.json --config NAME --list-rungs

The configuration is the shipped `examples/machine-*-shaped.toml` the problem
names, changed in exactly four places: the coil currents ( held supplies at the
reference's, free supplies at the start's ) with the free ones `Free`; an
`[inverse]` table of the problem's targets; the rung's degree and mesh sizes;
and the initial guess. Everything else -- profiles, `I_p`, the exterior, the
support sweeps, the retry rung -- is the shipped file's.

THE LADDER is degree against mesh size, every [mesh.generate] size scaled
together, because MEQ's claim (M-157) is high order on coarse meshes and a
ladder that only refines h at one degree cannot see it. The last rung is the
accuracy datum.

THE STARTS.
  cold    the design blob carrying I_p plus the conductors at the start's
          currents (mkcoldguess.py --currents), [source] PsiAxis at that
          guess's peak, free supplies at the problem's vacuum seed
  seeded  the same with the free supplies at the reference's currents
  warm    this arm's own converged run at the same rung: its psi, its mesh
          ( so nothing is re-meshed on the clock ), its psi_ax and X-point as
          the border's starting values, its free currents -- and the SHIFTED
          targets

THE ROWS come from the run's `.nc`, `inverse_target_residuals`, MEQ's own
evaluation of each target row on its converged field; MEQ poses the isoflux
rows first, so they are reordered to the contract's X-point-first order.
`finished` means the driver exited 0 or 4 ( 4: a best fit, which on a
consistent row set is a symmetry the mesh broke -- `symmetric` in the problem
turns on [mesh.generate] Symmetric where the rank is below the row count ).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
os.environ.setdefault("FGSREF_OUT", "/tmp")

from fgsref import CASES                                        # noqa: E402
from circuits import wiring                                     # noqa: E402
from inverse import rewire, set_key                             # noqa: E402
from inverse_compare import put_key                             # noqa: E402
from make_diverted_case import disc_radius                      # noqa: E402

#: ( degree, scale on every [mesh.generate] size ), coarse to fine.
#:
#: NO k = 1 RUNG, AND THE REASON IS A FINDING. At k = 1 the flux is linear in
#: each element, so the X-point border's Newton step rests on a grad q that is
#: constant per element and jumps between them, and on DIII-D shaped the first
#: bordered Newton throws the X-point unknown out of the mesh EVEN RESTARTED
#: FROM ITS OWN CONVERGED ANSWER WITH THE SAME TARGETS -- 41 iterations on the
#: retry rung against 2 at k = 2. A k = 1 rung would time that defect, not the
#: cost of an accuracy.
LADDER = [(2, 1.0), (2, 0.5), (3, 1.0), (3, 0.5), (4, 0.5), (4, 0.35)]
MESH_KEYS = ("Size", "PlasmaSize", "CoilSize")


def rungs():
    return [dict(rung=i, label="k%d h%.2f" % (k, s), degree=k, scale=s,
                 finest=(i == len(LADDER) - 1))
            for i, (k, s) in enumerate(LADDER)]


def table_value(text, table, key):
    m = re.search(r"(?ms)^\[%s\]\s*$(.*?)(?=^\[|\Z)" % re.escape(table), text)
    v = re.search(r"(?m)^%s = (\S+)" % re.escape(key), m.group(1)) if m else None
    return v.group(1) if v else None


def provenance():
    sha = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain",
                            "--untracked-files=no"],
                           capture_output=True, text=True).stdout.strip()
    version = subprocess.run([os.path.join(ROOT, "build", "meq"), "--version"],
                             capture_output=True, text=True).stdout.strip()
    lib = os.path.join(ROOT, "..", "mfem", "install", "lib", "libmfem.a")
    return dict(code_version=version, tree=sha + ("-dirty" if dirty else ""),
                libmfem_mtime=os.path.getmtime(lib) if os.path.exists(lib) else None,
                notes="shipped shaped example; [inverse] Track = false rows")


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
    rel = os.path.relpath(out, ROOT)
    stem = "meq"

    case = {c["name"]: c for c in CASES}[problem["fgsref_case"]]
    wires = wiring(case["make"]())
    free = [s["name"] for s in problem["supplies"] if s["free"]]
    currents = {s["name"]: s["reference_current"] for s in problem["supplies"]}
    text = open(problem["meq_example"]).read()

    if a.start == "warm":
        state = json.load(open(os.path.join(a.state, "state.json")))
        currents.update(state["currents"])
        xpts, iso = (problem["shifted_xpoint_targets"],
                     problem["shifted_isoflux_targets"])
    else:
        if a.start == "cold":
            currents.update(problem["vacuum_seed"])
        xpts, iso = problem["xpoint_targets"], problem["isoflux_targets"]

    text = rewire(text, wires, currents, free=set(free))
    text = set_key(text, "discretisation", "PolynomialDegree",
                   str(rung["degree"]))
    for key in MESH_KEYS:
        value = table_value(text, "mesh.generate", key)
        if value is not None:
            text = set_key(text, "mesh.generate", key,
                           "%.6g" % (float(value)*rung["scale"]))
    if problem["symmetric"] and problem["rank"] < problem["rows"]:
        text = put_key(text, "mesh.generate", "Symmetric", "true")

    if a.start == "warm":
        text = set_key(text, "mesh", "File", '"%s"' % state["mesh"])
        text = set_key(text, "initialguess", "File", '"%s"' % state["psi"])
        text = set_key(text, "initialguess", "MeshFile",
                       '"%s"' % state["guess_mesh"])
        text = set_key(text, "source", "PsiAxis", "%.12e" % state["psi_axis"])
        text = set_key(text, "boundary.xpoint", "R", "%.10f" % state["xpoint"][0])
        text = set_key(text, "boundary.xpoint", "Z", "%.10f" % state["xpoint"][1])
    else:
        text = set_key(text, "mesh", "File",
                       '"%s"' % os.path.join(rel, stem + ".msh"))
        # The cold guess: the design blob plus every conductor at the start's
        # currents, which mkcoldguess.py takes as per-conductor totals.
        totals = {}
        for s in problem["supplies"]:
            for m in s["members"]:
                totals[m["conductor"]] = m["multiplier"]*currents[s["name"]]
        jpath = os.path.join(out, "currents.json")
        json.dump(totals, open(jpath, "w"), indent=1)
        gmesh, ggf = os.path.join(out, "guess.mesh"), os.path.join(out, "guess.gf")
        import numpy as np
        npz = np.load(problem["reference_npz"], allow_pickle=True)
        made = subprocess.run(
            [sys.executable, os.path.join(HERE, "mkcoldguess.py"),
             problem["reference_npz"], gmesh, ggf, "%g" % disc_radius(npz),
             "32", "--currents", jpath], capture_output=True, text=True)
        if made.returncode != 0:
            raise SystemExit(made.stderr)
        peak = [l for l in made.stdout.splitlines() if l.startswith("PsiAxis ")]
        text = set_key(text, "initialguess", "File", '"%s"' % ggf)
        text = set_key(text, "initialguess", "MeshFile", '"%s"' % gmesh)
        text = set_key(text, "source", "PsiAxis",
                       "%.10e" % float(peak[0].split()[1]))
    text = set_key(text, "output", "Prefix", '"%s"' % os.path.join(rel, stem))

    targets = ["", "[inverse]", ""]
    for R, Z in xpts:
        targets += ["[[inverse.xpoint]]", "R = %r" % R, "Z = %r" % Z,
                    "Track = false", ""]
    for R1, Z1, R2, Z2 in iso:
        targets += ["[[inverse.isoflux]]", "R1 = %r" % R1, "Z1 = %r" % Z1,
                    "R2 = %r" % R2, "Z2 = %r" % Z2, ""]
    text = ("# race_arm_meq.py: %s, rung %s, %s start. Generated.\n\n"
            % (a.config, rung["label"], a.start)
            + text.rstrip("\n") + "\n" + "\n".join(targets))
    toml = os.path.join(out, stem + ".toml")
    open(toml, "w").write(text)

    log = os.path.join(out, stem + ".log")
    with open(log, "w") as fh:
        code = subprocess.run([os.path.join(ROOT, "build", "meq-run"), toml],
                              cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT
                              ).returncode
    wall = time.time() - t0
    body = open(log).read()

    record = dict(code="meq", config=a.config, rung=a.rung,
                  rung_label=rung["label"], start=a.start, unknowns=None,
                  finished=False, status="exit %d" % code, wall_s=wall,
                  solve_s=None, currents={}, psi_axis=None, psi_bndry=None,
                  axis=None, xpoints=[], rows=[], iterations=None, state=None,
                  provenance=provenance())
    m = re.search(r"converged in (\d+) Newton iterations on (\d+) elements, "
                  r"degree (\d+)", body)
    if m:
        k, n = int(m.group(3)), int(m.group(2))
        record["iterations"] = int(m.group(1))
        # Element-interior potential dofs, (k+1)(k+2)/2 per triangle: the
        # size the Pareto axis needs, comparable across MEQ's own rungs.
        record["unknowns"] = n*(k + 1)*(k + 2)//2
    nc = os.path.join(out, stem + ".nc")
    if code in (0, 4) and os.path.exists(nc):
        from netCDF4 import Dataset
        d = Dataset(nc)
        names = d.getncattr("inverse_circuits").split()
        values = [float(v) for v in str(d.getncattr("inverse_currents")).split()]
        record["currents"] = dict(zip(names, values))
        res = [float(v) for v in
               str(d.getncattr("inverse_target_residuals")).split()]
        niso = len(iso)
        record["rows"] = res[niso:] + res[:niso]
        record["psi_axis"] = float(d.getncattr("psi_axis"))
        record["psi_bndry"] = float(d.getncattr("psi_boundary"))
        record["axis"] = [float(d.getncattr("axis_r")), float(d.getncattr("axis_z"))]
        record["xpoints"] = [[float(d.getncattr("xpoint_r")),
                              float(d.getncattr("xpoint_z"))]]
        # A run with no located magnetic axis writes no axis_normalised_flux.
        located = "axis_normalised_flux" in d.ncattrs()
        # A SUPPORT LOOP THAT DID NOT SETTLE IS NOT AN ANSWER, and neither is
        # one that never ran: the file asks for sweeps, so a missing attribute
        # means the solve bypassed them. That is how the retry rung reported
        # its cold guess's support on five of F shaped's six cold rungs, every
        # row met and 8e-03 from the settled answer.
        swept = re.search(r"^\s*PlasmaSupportSweeps\s*=\s*([1-9]\d*)", text,
                          re.MULTILINE) is not None
        settled = (not swept) or (
            "plasma_support_settled" in d.ncattrs()
            and int(d.getncattr("plasma_support_settled")) == 1)
        record["support_sweeps"] = (int(d.getncattr("plasma_support_sweeps"))
                                    if "plasma_support_sweeps" in d.ncattrs()
                                    else None)
        record["finished"] = bool(located) and settled and \
            set(free) <= set(names)
        record["status"] = ("converged" if code == 0 else "best fit") \
            if record["finished"] else ("no located axis" if not located
                                        else "support not settled"
                                        if not settled else "circuits missing")
        # The state a warm start restarts from: this run's own field and mesh.
        record["state"] = out
        json.dump(dict(currents=record["currents"], psi_axis=record["psi_axis"],
                       xpoint=record["xpoints"][0],
                       psi=os.path.join(out, stem + "_psi.gf"),
                       guess_mesh=os.path.join(out, stem + ".mesh"),
                       mesh=(state["mesh"] if a.start == "warm"
                             else os.path.join(out, stem + ".msh"))),
                  open(os.path.join(out, "state.json"), "w"), indent=1)
    json.dump(record, open(os.path.join(out, "record.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
