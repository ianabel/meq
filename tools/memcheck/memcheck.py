#!/usr/bin/env python3
"""Run MEQ's shipped examples under valgrind's memcheck, and say what popped.

    tools/memcheck/memcheck.py                  every example within the budget
    tools/memcheck/memcheck.py --cases a,b      just these (stems of examples/*.toml)
    tools/memcheck/memcheck.py --native-only    the control pass alone: exit and time

WORTH RUNNING PERIODICALLY, AND THE SUITE IS NOT A SUBSTITUTE. The suite runs
test harnesses; this runs the driver on the inputs a user runs, through the same
front end -- meq-run makes the mesh, meq solves -- and RotatingMirrorAnisotropic's
crash (an out-of-bounds read decided by allocation order, c5a1ecd) is the kind of
defect that passes a suite for months and that memcheck reports every time.

EACH CASE RUNS TWICE, IN runs/memcheck/<stamp>/<case>/:

  native     the same binary without valgrind -- the control: its exit code and
             its wall time, which decides whether the case is run under valgrind
             at all (--budget, native seconds; memcheck costs 20-50x).
  memcheck   --track-origins, definite leaks counted as errors, every error
             a failure (--error-exitcode).

A case POPS when memcheck reports an error, or when its exit code under valgrind
differs from the native one -- a solve that takes another branch under valgrind's
serialised threads and exact x87-free arithmetic is worth reading too. A case
whose native exit is not 0 is reported, not judged: some examples are meant to
exit 2 or 4, and a native failure is the driver's report, not memcheck's.

THE BINARY IS AN INSTALL, as every tool's is (tools/meqpaths.py), and it should
carry debug information or the reports name functions without lines:

    cmake -S source -B builds/valgrind -DCMAKE_BUILD_TYPE=RelWithDebInfo \\
          -DMEQ_BUILD_TESTS=OFF -DMEQ_PYTHON=OFF \\
          -DCMAKE_INSTALL_PREFIX=$PWD/install/valgrind  [+ the usual -DMFEM_DIR=]
    cmake --build builds/valgrind -j6 && cmake --install builds/valgrind
    MEQ_BUILD=valgrind source/tools/memcheck/memcheck.py

memcheck.supp beside this file holds the suppressions for code MEQ does not own;
each says why. A suppression for MEQ's own code is a defect hidden, not noise.
"""

import argparse
import concurrent.futures
import datetime
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import meqpaths  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = meqpaths.source("examples")
SUPPRESSIONS = os.path.join(HERE, "memcheck.supp")


def environment(threads):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    # As ctest sets it: PARDISO's own threads are not what is being checked.
    env["MKL_NUM_THREADS"] = "1"
    return env


def run(command, directory, env, log, timeout):
    start = time.monotonic()
    with open(log, "w") as out:
        try:
            code = subprocess.run(command, cwd=directory, env=env, stdout=out,
                                  stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            code = "timeout"
    return code, time.monotonic() - start


def prepare(toml, directory, env, timeout):
    """Make @a directory a place the run can start: the mesh, through the
    user's front end, when the configuration asks for one. Returns the flags
    meq then needs, or None when the mesh could not be made."""
    os.makedirs(directory, exist_ok=True)
    asked = subprocess.run([meqpaths.meq(), "--mesh-command", toml], cwd=directory,
                           capture_output=True, text=True)
    if asked.returncode != 0 or not asked.stdout.strip():
        return []
    code, _ = run([meqpaths.meq_run(), "--mesh-only", toml], directory, env,
                  os.path.join(directory, "mesh.log"), timeout)
    return ["--mesh-ready"] if code == 0 else None


def one_case(stem, args, root):
    toml = os.path.join(EXAMPLES, stem + ".toml")
    result = {"case": stem}

    native = os.path.join(root, stem, "native")
    env = environment(args.threads)
    flags = prepare(toml, native, env, args.timeout)
    if flags is None:
        result["verdict"] = "MESH FAILED"
        return result
    code, seconds = run([meqpaths.meq()] + flags + [toml], native, env,
                        os.path.join(native, "meq.log"), args.timeout)
    result.update(native=code, seconds=seconds)
    if args.native_only:
        result["verdict"] = "native"
        return result
    if code == "timeout" or seconds > args.budget:
        result["verdict"] = "skipped: %.0f s native" % seconds
        return result

    checked = os.path.join(root, stem, "memcheck")
    env = environment(args.valgrind_threads)
    flags = prepare(toml, checked, env, args.timeout)
    report = os.path.join(checked, "valgrind.log")
    command = ["valgrind", "--tool=memcheck", "--error-exitcode=99",
               "--track-origins=yes", "--num-callers=40",
               "--leak-check=full", "--show-leak-kinds=definite",
               "--errors-for-leak-kinds=definite", "--log-file=" + report]
    if os.path.exists(SUPPRESSIONS):
        command.append("--suppressions=" + SUPPRESSIONS)
    code, seconds = run(command + [meqpaths.meq()] + flags + [toml], checked, env,
                        os.path.join(checked, "meq.log"), args.timeout*50)
    result.update(memcheck=code, vseconds=seconds)
    errors = None
    if os.path.exists(report):
        with open(report) as f:
            found = re.findall(r"ERROR SUMMARY: (\d+) errors from (\d+) contexts", f.read())
        if found:
            errors = int(found[-1][0])
            result["summary"] = "%s errors / %s contexts" % found[-1]
    # Clean only when memcheck read its summary as zero AND the driver exited
    # as it did natively; anything else -- an error, a crash, a timeout, a
    # different branch -- is for reading.
    if errors is None:
        # No summary read: valgrind did not run, or died before it wrote one.
        # Not a verdict on MEQ -- the harness is what to read (meq.log).
        result["verdict"] = "NO REPORT"
        return result
    clean = errors == 0 and code == result["native"]
    result["verdict"] = "clean" if clean else "POPPED"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cases", help="comma-separated example stems; default every one")
    parser.add_argument("--budget", type=float, default=60.0,
                        help="largest native wall time, s, run under valgrind (60)")
    parser.add_argument("--jobs", type=int, default=6, help="cases at once (6)")
    parser.add_argument("--threads", type=int, default=4,
                        help="OMP_NUM_THREADS for the native run (4)")
    parser.add_argument("--valgrind-threads", type=int, default=2,
                        help="OMP_NUM_THREADS under valgrind, which serialises them; "
                             "more than one keeps the threaded paths in play (2)")
    parser.add_argument("--timeout", type=float, default=1800.0,
                        help="native timeout, s; memcheck gets 50x it (1800)")
    parser.add_argument("--native-only", action="store_true",
                        help="the control pass alone")
    parser.add_argument("--out", help="output directory; default runs/memcheck/<stamp>")
    args = parser.parse_args()

    if not os.access(meqpaths.meq(), os.X_OK):
        sys.exit("no meq at %s: set MEQ_BUILD or MEQ_INSTALL" % meqpaths.meq())
    if not args.native_only and shutil.which("valgrind") is None:
        sys.exit("valgrind is not on PATH")
    stems = sorted(f[:-5] for f in os.listdir(EXAMPLES) if f.endswith(".toml"))
    if args.cases:
        wanted = args.cases.split(",")
        unknown = [c for c in wanted if c not in stems]
        if unknown:
            sys.exit("not an example: " + ", ".join(unknown))
        stems = wanted
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    # Absolute: every path below is handed to a process running in a case's
    # own directory, valgrind's --log-file among them.
    root = os.path.abspath(args.out or meqpaths.runs("memcheck", stamp))
    os.makedirs(root, exist_ok=True)
    version = subprocess.run([meqpaths.meq(), "--version"],
                             capture_output=True, text=True).stdout.strip()
    print("%s -> %s" % (version, root), flush=True)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(one_case, s, args, root): s for s in stems}
        for future in concurrent.futures.as_completed(futures):
            r = future.result()
            results.append(r)
            print("  %-36s native %-7s %6.1f s  %s  %s" % (
                r["case"], r.get("native"), r.get("seconds", 0.0),
                r.get("summary", ""), r["verdict"]), flush=True)

    lines = ["%s\n" % version]
    for r in sorted(results, key=lambda r: r["case"]):
        lines.append("%-36s native %-7s %7.1f s   memcheck %-7s %8.1f s   %-26s %s\n" % (
            r["case"], r.get("native"), r.get("seconds", 0.0), r.get("memcheck", "-"),
            r.get("vseconds", 0.0), r.get("summary", "-"), r["verdict"]))
    with open(os.path.join(root, "summary.txt"), "w") as f:
        f.writelines(lines)
    popped = [r["case"] for r in results if r["verdict"] in ("POPPED", "MESH FAILED", "NO REPORT")]
    print("\n%d cases, %d popped%s" % (len(results), len(popped),
                                       (": " + ", ".join(sorted(popped))) if popped else ""))
    return 1 if popped else 0


if __name__ == "__main__":
    sys.exit(main())
