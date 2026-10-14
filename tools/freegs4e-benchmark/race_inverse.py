"""THE INVERSE RACE: four codes, one posed problem, time to accuracy.

INVERSE-RACE-PLAN.md is the design and RACE-ARMS.md the contract every arm
implements; inverse_problem.py writes the problem file they all read. This is
the orchestrator: it decides nothing about any code's answer except by the
acceptance test below, and it records everything.

    venv/bin/python race_inverse.py run STAGE --config machine-f-diiid-shaped
        [--problem PROBLEM] [--codes meq,nice,freegs4e,freegsnke]
        [--starts cold,warm] [--repeats 3] [--max-rung N] [--row-tol 1e-6]
        [--load-max 1.0] [--quiet-wait 900] [--dry-run]
    venv/bin/python race_inverse.py summarise STAGE

A bare STAGE name is runs/freegs4e-benchmark/race-runs/STAGE/; anything with a
slash in it is a directory. `summarise` also finds a stage recorded in this
directory's tracked race-runs/, and writes its summary under runs/ rather than
over the record. PROBLEM defaults to runs/freegs4e-benchmark/ref-interp/
problem.json, inverse_problem.py's output.

THE PROTOCOL, each clause of which this tree learned once ( race.py,
race_nke.py, MEASUREMENTS.md M-102, M-111, M-155, M-157 ).

  * ONE CONFIGURATION PER INVOCATION. A multi-case run times each case on a
    machine the last one left at full load -- 9.49 s against 5.41 s settled,
    M-111. `--config` takes one name.
  * INTERLEAVED PER RUNG. For each rung index every arm that has one runs
    before any arm goes to the next, and the arm order rotates by rung, so a
    machine that drifts drifts through all of them.
  * QUIET BEFORE EVERY RUN. The one-minute load average is sampled just before
    each run and the run waits, at most `--quiet-wait` seconds, for it to fall
    under `--load-max`. A run that could not wait is flagged `noisy` in its
    record rather than dropped.
  * MEDIAN OF `--repeats`. Each timed run is repeated in a fresh directory;
    the median wall is the figure, the first repeat's answer the answer.
  * THE FINEST RUNG IS THE DATUM, run once and not timed as a contestant.
  * THREADS MATCHED, at the physical core count ( cores.py ), OMP and MKL alike.
  * ACCEPTANCE IS THE ROWS. A run counts only if the arm reports `finished`
    AND every target row it evaluated from its own field is within
    `--row-tol`. freegsnke reports SUCCESS with rows at 0.1 T (M-255).

THE TWO STARTS ( RACE-ARMS.md ). `cold` is timed for every rung. `warm` needs
each arm's own converged design-target state at the same rung: the first
accepted cold repeat if there is one, else an UNTIMED `seeded` run; an arm with
neither is a DNF for the warm start too.

THE SUMMARY. For each start and code: per rung the median wall, the size, and
the error of the free currents against the SAME code's finest accepted run
( max |dI| over the free supplies over the largest free current ) -- and for
each accuracy target the cheapest rung that reaches it. Separately, the
cross-code floor: the same measure between every pair of codes' finest runs.
Below that floor no claim about which code is more ACCURATE is supported, and
the summary says so.
"""

import argparse
import itertools
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import benchpaths                                               # noqa: E402
from cores import physical_cores                                # noqa: E402

VENV = benchpaths.VENV
ARMS = {
    "meq": [VENV, os.path.join(HERE, "race_arm_meq.py")],
    "nice": [VENV, os.path.join(benchpaths.NICE, "meq-race", "inverse",
                                "race_arm.py")],
    "freegs4e": [VENV, os.path.join(HERE, "race_arm_fgs.py")],
    "freegsnke": [VENV, os.path.join(HERE, "race_arm_nke.py")],
}
TARGETS = [1e-2, 1e-3, 1e-4, 1e-5]


def env():
    threads = str(physical_cores())
    return dict(os.environ, OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads,
                PYTHONPATH=benchpaths.FREEGS4E, MPLBACKEND="Agg",
                FGSREF_OUT=os.environ.get("FGSREF_OUT", "/tmp"))


def ladder(code, problem, config):
    out = subprocess.run(ARMS[code] + ["--problem", problem, "--config", config,
                                       "--list-rungs"],
                         capture_output=True, text=True, env=env())
    if out.returncode != 0:
        raise SystemExit("%s --list-rungs failed:\n%s" % (code, out.stderr))
    return json.loads(out.stdout.strip().splitlines()[-1])


def quiet(load_max, wait):
    """Wait, bounded, for a quiet machine. Returns ( load, waited, noisy )."""
    t0 = time.time()
    while True:
        load = os.getloadavg()[0]
        if load <= load_max:
            return load, time.time() - t0, False
        if time.time() - t0 > wait:
            return load, time.time() - t0, True
        time.sleep(15)


def accepted(record, row_tol):
    return bool(record.get("finished")) and bool(record.get("rows")) and \
        max(abs(r) for r in record["rows"]) <= row_tol


def run_arm(code, problem, config, rung, start, outdir, state, args, log):
    os.makedirs(outdir, exist_ok=True)
    cmd = ARMS[code] + ["--problem", problem, "--config", config,
                        "--rung", str(rung), "--start", start, "--out", outdir]
    if state:
        cmd += ["--state", state]
    load, waited, noisy = quiet(args.load_max, args.quiet_wait)
    t0 = time.time()
    with open(os.path.join(outdir, "arm.log"), "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env(),
                            cwd=outdir).returncode
    outer = time.time() - t0
    path = os.path.join(outdir, "record.json")
    record = json.load(open(path)) if os.path.exists(path) else dict(
        code=code, config=config, rung=rung, start=start, finished=False,
        status="arm exited %d with no record" % rc, rows=[], currents={})
    record.update(load_before=load, waited_s=waited, noisy=noisy,
                  outer_wall_s=outer, arm_exit=rc, outdir=outdir,
                  accepted=accepted(record, args.row_tol))
    log.write(json.dumps(record) + "\n")
    log.flush()
    print("  %-9s %-6s rung %d %-10s %s  wall %7.1f s  %s%s" % (
        code, start, rung, record.get("rung_label", ""),
        "ACCEPTED" if record["accepted"] else "rejected",
        record.get("wall_s") or outer, record.get("status"),
        "  NOISY" if noisy else ""), flush=True)
    return record


def plan(args):
    rungs = {c: ladder(c, args.problem, args.config) for c in args.codes}
    top = max(len(r) for r in rungs.values())
    if args.max_rung is not None:
        top = min(top, args.max_rung + 1)
    return rungs, top


def stage(name, recorded=False):
    """A stage's directory: a bare name is runs/freegs4e-benchmark/race-runs/
    <name>/, anything with a separator a path. `recorded` takes the tracked
    race-runs/<name>/ beside this script instead when only it holds the
    records -- a stage summarised again, whose run products are under runs/."""
    if os.sep in name or (os.altsep and os.altsep in name):
        return os.path.abspath(name)
    here = os.path.join(benchpaths.OUT, "race-runs", name)
    tracked = os.path.join(HERE, "race-runs", name)
    records = "records.jsonl"
    if (recorded and not os.path.exists(os.path.join(here, records))
            and os.path.exists(os.path.join(tracked, records))):
        return tracked
    return here


def cmd_run(args):
    out = stage(args.outdir)
    os.makedirs(out, exist_ok=True)
    rungs, top = plan(args)
    print("config %s, codes %s, starts %s, repeats %d" % (
        args.config, args.codes, args.starts, args.repeats))
    for c in args.codes:
        print("  %-9s %s" % (c, ", ".join("%d:%s%s" % (r["rung"], r["label"],
                                                       " (datum)" if r["finest"]
                                                       else "")
                                          for r in rungs[c])))
    if args.dry_run:
        return
    json.dump(dict(vars(args), rungs=rungs, host=dict(zip(("sysname", "nodename", "release", "version", "machine"),
                                   os.uname())),
                   physical_cores=physical_cores(), started=time.time()),
              open(os.path.join(out, "plan.json"), "w"), indent=1, default=str)
    log = open(os.path.join(out, "records.jsonl"), "a")
    states = {}
    for start in args.starts:
        print("== %s start" % start, flush=True)
        for r in range(top):
            order = args.codes[r % len(args.codes):] + args.codes[:r % len(args.codes)]
            for code in order:
                if r >= len(rungs[code]):
                    continue
                datum = rungs[code][r]["finest"]
                state = None
                if start == "warm":
                    state = states.get((code, r))
                    if state is None:
                        seed_dir = os.path.join(out, args.config, "seeded", code,
                                                "r%02d" % r)
                        rec = run_arm(code, args.problem, args.config, r,
                                      "seeded", seed_dir, None, args, log)
                        state = rec.get("state") if rec["accepted"] else None
                        states[(code, r)] = state
                    if state is None:
                        print("  %-9s warm   rung %d: DNF -- no converged "
                              "design-target state" % (code, r), flush=True)
                        continue
                for k in range(1 if datum else args.repeats):
                    d = os.path.join(out, args.config, start, code,
                                     "r%02d" % r, "rep%d" % k)
                    rec = run_arm(code, args.problem, args.config, r, start, d,
                                  state, args, log)
                    if start == "cold" and rec["accepted"] and (code, r) not in states:
                        states[(code, r)] = rec.get("state")
    log.close()


def load_records(out):
    with open(os.path.join(out, "records.jsonl")) as fh:
        return [json.loads(l) for l in fh if l.strip()]


def current_error(a, b, free):
    scale = max(abs(b["currents"][f]) for f in free)
    return max(abs(a["currents"][f] - b["currents"][f]) for f in free)/scale


def cmd_summarise(args):
    out = stage(args.outdir, recorded=True)
    recs = [r for r in load_records(out) if r.get("start") != "seeded"]
    plan_ = json.load(open(os.path.join(out, "plan.json")))
    # A recorded plan names the problem file where it was written, which for
    # the tracked stages is this directory's old ref-interp/.
    path = plan_["problem"]
    if not os.path.exists(path):
        path = os.path.join(benchpaths.OUT, "ref-interp", os.path.basename(path))
    problem = {p["config"]: p for p in json.load(open(path))}
    # NOTHING IS WRITTEN INTO THE SOURCE TREE: a tracked stage's summary goes
    # to its namesake under runs/.
    written = out
    if os.path.commonpath([out, benchpaths.SOURCE]) == benchpaths.SOURCE:
        written = benchpaths.out("race-runs", os.path.basename(out))
    lines, summary = [], {}
    for config in sorted({r["config"] for r in recs}):
        free = [s["name"] for s in problem[config]["supplies"] if s["free"]]
        lines.append("## %s\n" % config)
        for start in ("cold", "warm"):
            group = [r for r in recs if r["config"] == config and r["start"] == start]
            if not group:
                continue
            lines.append("### %s start\n" % start)
            finest = {}
            for code in sorted({r["code"] for r in group}):
                mine = [r for r in group if r["code"] == code]
                # The datum is the arm's OWN finest rung, from its ladder --
                # not the finest that happened to be run.
                top = [x["rung"] for x in plan_["rungs"][code] if x["finest"]][0]
                datum = [r for r in mine if r["rung"] == top and r["accepted"]]
                finest[code] = datum[0] if datum else None
                lines.append("**%s**%s\n" % (code, "" if datum else
                                             " -- the finest rung has no accepted "
                                             "run: no self-error can be reported"))
                lines.append("| rung | size | median wall s | repeats | accepted | "
                             "dI vs own finest | psi_ax vs own finest |")
                lines.append("|---|---|---|---|---|---|---|")
                rows = []
                for rung in sorted({r["rung"] for r in mine}):
                    rr = [r for r in mine if r["rung"] == rung]
                    ok = [r for r in rr if r["accepted"]]
                    walls = [r["wall_s"] for r in ok if r.get("wall_s")]
                    wall = statistics.median(walls) if walls else None
                    err = perr = None
                    if ok and datum and rung != top:
                        err = current_error(ok[0], datum[0], free)
                        perr = abs(ok[0]["psi_axis"] - datum[0]["psi_axis"]) / \
                            abs(datum[0]["psi_axis"])
                    rows.append(dict(rung=rung, label=rr[0].get("rung_label"),
                                     size=rr[0].get("unknowns"), wall=wall,
                                     err=err, accepted=len(ok), runs=len(rr)))
                    lines.append("| %s%s | %s | %s | %d | %d/%d | %s | %s |" % (
                        rr[0].get("rung_label"), " (datum)" if rung == top else "",
                        rr[0].get("unknowns"),
                        "%.2f" % wall if wall else "--", len(rr), len(ok), len(rr),
                        "%.2e" % err if err is not None else "--",
                        "%.2e" % perr if perr is not None else "--"))
                lines.append("")
                cheapest = {}
                for t in TARGETS:
                    hit = [x for x in rows if x["err"] is not None and x["err"] <= t
                           and x["wall"]]
                    cheapest[t] = min(hit, key=lambda x: x["wall"]) if hit else None
                summary.setdefault(config, {}).setdefault(start, {})[code] = dict(
                    rungs=rows, cheapest={str(t): v for t, v in cheapest.items()})
            lines.append("**Cheapest run reaching each self-consistent accuracy "
                         "(currents), median wall s:**\n")
            codes = sorted(summary[config][start])
            lines.append("| target | " + " | ".join(codes) + " |")
            lines.append("|---|" + "---|"*len(codes))
            for t in TARGETS:
                cells = []
                for code in codes:
                    c = summary[config][start][code]["cheapest"][str(t)]
                    cells.append("%.2f (%s)" % (c["wall"], c["label"]) if c else "--")
                lines.append("| %.0e | %s |" % (t, " | ".join(cells)))
            lines.append("")
            lines.append("**Cross-code floor** (finest accepted runs, dI over "
                         "the free supplies):\n")
            floor = {}
            for a, b in itertools.combinations(sorted(finest), 2):
                if finest[a] and finest[b]:
                    floor["%s-%s" % (a, b)] = current_error(finest[a], finest[b], free)
                    lines.append("* %s - %s: %.2e" % (a, b, floor["%s-%s" % (a, b)]))
            summary[config][start]["floor"] = floor
            if not floor:
                lines.append("* none: no two codes have an accepted run on their "
                             "finest rung")
            if floor:
                lines.append("\nNo accuracy target below %.1e distinguishes the "
                             "codes' ANSWERS; the table above compares each "
                             "code with itself.\n" % max(floor.values()))
    text = "\n".join(lines)
    open(os.path.join(written, "summary.md"), "w").write(text)
    json.dump(summary, open(os.path.join(written, "summary.json"), "w"), indent=1,
              default=str)
    print(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("outdir")
    r.add_argument("--config", required=True)
    r.add_argument("--problem",
                   default=os.path.join(benchpaths.OUT, "ref-interp",
                                        "problem.json"))
    r.add_argument("--codes", default="meq,nice,freegs4e,freegsnke")
    r.add_argument("--starts", default="cold,warm")
    r.add_argument("--repeats", type=int, default=3)
    r.add_argument("--max-rung", type=int)
    r.add_argument("--row-tol", type=float, default=1e-6)
    r.add_argument("--load-max", type=float, default=1.0)
    r.add_argument("--quiet-wait", type=float, default=900.0)
    r.add_argument("--dry-run", action="store_true")
    s = sub.add_parser("summarise")
    s.add_argument("outdir")
    a = ap.parse_args()
    if a.cmd == "run":
        a.problem = os.path.abspath(a.problem)
        a.codes = [c for c in a.codes.split(",") if c]
        a.starts = [s for s in a.starts.split(",") if s]
        unknown = [c for c in a.codes if c not in ARMS]
        if unknown:
            raise SystemExit("unknown codes %s" % unknown)
        cmd_run(a)
    else:
        cmd_summarise(a)


if __name__ == "__main__":
    main()
