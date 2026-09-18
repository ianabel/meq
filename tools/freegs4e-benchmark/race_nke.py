#!/usr/bin/env python
"""RACE MEQ AGAINST freegsnke ON THE MAST-U FORWARD CASE.

`race.py` races MEQ against freegs4e and its docstring carries the warning that
its freegs4e arm is an INVERSE solve while MEQ's is forward, so a multi-case run
there measures the machine as much as the codes.  This one does not have that
defect by construction: `mastu_reference.py`'s case is `constrain=None`, so both
arms are handed the same twelve coil currents and both solve forward.

THE PROTOCOL, and every clause of it is a measurement this tree has already got
wrong once.

  INTERLEAVED, NEVER ARM-BY-ARM.  MEASUREMENTS.md M-102: a 2.2x "regression"
  was reported from timings taken after a race against a control read before
  it.  The arms alternate per level, so a machine that drifts drifts through
  both.

  QUIET FIRST, AND CHECKED AGAIN BEFORE EVERY LEVEL -- NOT AFTERWARDS.  The
  first version of this compared the load average before AND after and would
  have refused every successful race, because after a multi-minute run the load
  average is dominated by the race itself.  Measured: the reference ladder left
  the load at 1.36 having started at 0.08, entirely its own doing.  So the
  sample that means anything is taken just BEFORE each level, where the only
  thing running is this script, and the threshold allows for one busy arm.

  MEDIAN OF REPEATS, NOT BEST-OF.  tests/performance/ draws the distinction:
  best-of is right for a minimum wall time, a median for a leg share.  This is
  comparing two codes' whole solves, so the median is what does not reward
  whichever arm happened to get a quiet slice.

  PROVENANCE ON EVERY ROW.  Both trees' SHAs and libmfem's mtime, because this
  session measured an XP-3 iteration count move by a whole step on a library
  change with no MEQ edit at all.

WHAT IS NOT BUILT YET, said plainly rather than left to be discovered: the MEQ
arm needs a TOML for THIS machine and there is not one.  `examples/machine-g-mastu.toml`
is a DIFFERENT MAST-U -- it comes from the freegs4e benchmark's G case and has
neither these conductors nor this limiter.  What a new one needs:

  * a half-disc mesh reaching r = 0 with all twelve active coils inside it
    (`tools/mesh/halfdisc.py`), since FB-5's exterior coupling wants the
    semicircle.  freegsnke's grid starts at Rmin = 0.1 and does not reach the
    axis, so the domains are not the same shape and MEQ's must be the larger.
  * the 47-point limiter polygon out of `MAST-U_like_limiter.pickle`, as a
    meshed surface for `[boundary.limiter] SurfaceAttribute`.
  * the twelve currents below, which this script prints with `--dump-currents`
    in the exact form `[[coils]]` wants.
  * `[source] Type = "mhd"` with `PlasmaCurrent = 6e5` and p'/ff' tabulated
    from freegsnke's `ConstrainPaxisIp` at alpha_m 1.8, alpha_n 1.2.

The passives need nothing: all 138 carry exactly zero current in this case,
checked rather than assumed, which is what makes it portable to a code with no
passive-conductor model.
"""
import argparse, json, os, pickle, statistics, subprocess, sys, time

# BOTH AXES, MATCHED, FOR BOTH ARMS -- which is what race.py has always done
# and what an earlier draft of this file got wrong by pinning the MEQ arm to
# MKL_NUM_THREADS=1 and leaving the freegsnke arm alone.
#
# MKL_NUM_THREADS=1 is the CTEST setting. It is pinned there because the
# bit-exactness assertions between assembly modes hold at that value and not
# above it; it is NOT the production setting. CLAUDE.md: under
# AssemblyMode::Threaded -- the default -- the element-local dense work is
# nested inside an active OpenMP region where MKL suppresses its own threading,
# so MKL_NUM_THREADS costs the assembly nothing, while the TRACE SOLVE runs on
# the master thread outside every region and takes every thread it is given.
# M-78 sizes it: the MKL threads are worth a further 7 to 9 per cent. Pinning
# them to 1 for a timed run hands that away and flatters whatever it is raced
# against.
from cores import physical_cores
THREADS = str(physical_cores())
MEQ = "/home/ian/projects/meq"
VENV = f"{MEQ}/tools/freegs4e-benchmark/venv/bin/python"
CURRENTS = "/home/ian/projects/freegsnke/examples/data/simple_diverted_currents_PaxisIp.pk"
LADDER = [(33, 65), (65, 129), (129, 257), (257, 513)]


def load_average():
    with open("/proc/loadavg") as handle:
        return float(handle.read().split()[0])


def wait_for_quiet(threshold, settle, timeout):
    """Block until the load has been under `threshold` for `settle` seconds.

    Not a process check.  CLAUDE.md records five waiters in one session spinning
    on `pgrep -f`, which matches the waiting shell's own command line, and `-x`
    failing silently on any binary over fifteen characters.  The load average is
    an artefact, which is what this tree's own rule says to wait on.
    """
    deadline = time.time() + timeout
    quiet_since = None
    while time.time() < deadline:
        if load_average() < threshold:
            quiet_since = quiet_since or time.time()
            if time.time() - quiet_since >= settle:
                return True
        else:
            quiet_since = None
        time.sleep(10)
    return False


def provenance():
    def sha(path):
        out = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True)
        return out.stdout.strip() or "?"
    lib = "/home/ian/projects/mfem/install/lib/libmfem.a"
    return {"meq": sha(MEQ), "freegsnke": sha("/home/ian/projects/freegsnke"),
            "freegs4e": sha("/home/ian/projects/freegs4e"),
            "libmfem_mtime": time.strftime(
                "%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(lib)))
            if os.path.exists(lib) else "?"}


def run_freegsnke(nx, ny, tolerance, out):
    """One freegsnke forward solve, timed, through mastu_reference.py."""
    started = time.perf_counter()
    done = subprocess.run(
        [VENV, f"{MEQ}/tools/freegs4e-benchmark/mastu_reference.py",
         "--out", out, "--levels", str(nx), "--repeats", "1",
         "--tolerance", str(tolerance)],
        capture_output=True, text=True,
        env={**os.environ, "MPLBACKEND": "Agg"})
    seconds = time.perf_counter() - started
    return dict(seconds=seconds, ok=done.returncode == 0,
                tail=done.stdout.strip().split("\n")[-1] if done.stdout else done.stderr[-200:])


def run_meq(config, out):
    """One MEQ forward solve, timed, through the driver.

    The WHOLE driver, not the solve leg alone: the freegsnke arm is timed the
    same way and a comparison of a solve against a solve-plus-output would
    flatter whichever side was measured more narrowly.  `mastu_reference.py`
    starts its clock after the machine is built for the same reason -- MEQ pays
    that cost as meshing, outside this timer.
    """
    started = time.perf_counter()
    done = subprocess.run([f"{MEQ}/build/meq-run", config],
                          capture_output=True, text=True,
                          env={**os.environ, "OMP_NUM_THREADS": THREADS,
                               "MKL_NUM_THREADS": THREADS})
    seconds = time.perf_counter() - started
    return dict(seconds=seconds, ok=done.returncode == 0, exit=done.returncode,
                tail=done.stdout.strip().split("\n")[-1] if done.stdout else done.stderr[-200:])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meq-config", default="",
                    help="MEQ TOML for THIS machine. Absent, only the freegsnke "
                         "arm runs and the race says so rather than comparing "
                         "against a different machine's config")
    ap.add_argument("--out", default="/tmp/race-nke")
    ap.add_argument("--tolerance", type=float, default=1e-9)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--levels", default="33,65,129")
    ap.add_argument("--require-quiet", action="store_true",
                    help="wait for a quiet machine first and refuse to publish "
                         "if it is loaded at the end")
    ap.add_argument("--quiet-threshold", type=float, default=0.5)
    ap.add_argument("--settle", type=float, default=300.0,
                    help="seconds the load must stay under the threshold")
    ap.add_argument("--dump-currents", action="store_true",
                    help="print the twelve active currents and exit")
    args = ap.parse_args()

    if args.dump_currents:
        with open(CURRENTS, "rb") as handle:
            currents = pickle.load(handle)
        for label, value in currents.items():
            if abs(float(value)) > 0.0:
                print(f"{label:12s} {float(value): .12e}")
        return 0

    if args.require_quiet:
        print(f"waiting for load < {args.quiet_threshold} sustained for "
              f"{args.settle:.0f}s ...", flush=True)
        if not wait_for_quiet(args.quiet_threshold, args.settle, 3600.0):
            print("the machine never went quiet; refusing to time anything",
                  file=sys.stderr)
            return 2
    before = load_average()
    # ONE BUSY ARM'S WORTH OF SLACK, because between levels this script's own
    # previous level is still decaying out of the average and that is not
    # contention -- it is the measurement.
    per_level_limit = args.quiet_threshold + 1.5
    intruded = []

    os.makedirs(args.out, exist_ok=True)
    wanted = [int(v) for v in args.levels.split(",") if v.strip()]
    ladder = [g for g in LADDER if g[0] in wanted]

    rows = []
    for nx, ny in ladder:
        if args.require_quiet:
            level_load = load_average()
            if level_load >= per_level_limit:
                intruded.append((nx, level_load))
        nke, meq = [], []
        for _ in range(max(1, args.repeats)):
            # INTERLEAVED WITHIN THE REPEAT, not arm after arm.  M-102.
            nke.append(run_freegsnke(nx, ny, args.tolerance, args.out))
            if args.meq_config:
                meq.append(run_meq(args.meq_config, args.out))
        row = dict(nx=nx, ny=ny,
                   freegsnke=statistics.median(r["seconds"] for r in nke),
                   freegsnke_ok=all(r["ok"] for r in nke))
        if meq:
            row.update(meq=statistics.median(r["seconds"] for r in meq),
                       meq_ok=all(r["ok"] for r in meq),
                       ratio=statistics.median(r["seconds"] for r in meq)
                             / max(statistics.median(r["seconds"] for r in nke), 1e-12))
        rows.append(row)
        print(f"  {nx:4d}x{ny:<5d} freegsnke {row['freegsnke']:8.2f}s"
              + (f"   MEQ {row['meq']:8.2f}s   ratio {row['ratio']:5.2f}x"
                 if meq else "   MEQ  (no config)"), flush=True)

    after = load_average()
    verdict = "clean"
    if args.require_quiet and intruded:
        verdict = ("REFUSED: something else was running. Load at the start of "
                   + ", ".join(f"{nx}x: {l:.2f}" for nx, l in intruded)
                   + f" against a per-level limit of {per_level_limit:.2f}. "
                     "These seconds are a measurement of the machine.")
        print(verdict, file=sys.stderr)

    with open(os.path.join(args.out, "race-nke.json"), "w") as handle:
        json.dump(dict(provenance=provenance(), rows=rows, verdict=verdict,
                       load_before=before, load_after=after,
                       repeats=args.repeats, tolerance=args.tolerance),
                  handle, indent=2)
    if not args.meq_config:
        print("\nOnly the freegsnke arm ran. See this file's docstring for what "
              "a MEQ config for THIS machine still needs; examples/machine-g-mastu.toml "
              "is a different MAST-U and must not be substituted.")
    return 0 if verdict == "clean" else 1


if __name__ == "__main__":
    sys.exit(main())
