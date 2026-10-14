#!/usr/bin/env python
"""MEQ AGAINST geq ON THE ROTATING MIRROR: A STRAIGHT RACE ON ONE INPUT.

The input is MR-4's problem, `reference/rotating/`: geq's two packs as twenty
fixed filaments, the skew profiles, the two labels at ( 0.10, 0 ) and
( 0.35, 0 ), R_ref and Gamma.  Each code is driven the way a user drives it:

  geq   `export_rotating.solve( n, "direct", ... )` in a fresh process: geq's
        direct free-boundary Picard on an n x n grid, with export_rotating.py's
        three changes (the direct boundary, analytic profiles, phi_0 tightened
        every step).  See README.md, *The rotating mirror*.
  MEQ   `examples/mirror-rotating.toml` through `build/meq`, with only the
        resolution ( [mesh] RefinementLevels, [discretisation]
        PolynomialDegree ), the threads and the output grid changed.  The mesh
        is generated once by `meq-run --mesh-only` and is not timed.

TWO RACES, as race_desc.py has them, and only the second decides anything.

  matched   each code at the resolution whose unknown count is nearest the
            other's: geq's n^2 grid values against MEQ's globally coupled trace
            unknowns, faces x ( k + 1 ).  A choice, printed with its definition.
  accuracy  the cheapest run of each code reaching each error target: a Pareto
            front rather than a pair of seconds.

THE ERROR IS MEASURED AT geq'S OWN NODES, the 65 x 65 every geq resolution
shares, R in [ 0.05, 0.63 ] and |z| <= 1.5, less geq's four grid corners
(README.md: there geq carries half the plasma's flux, M-263).  MEQ samples
them exactly: [output] GridRMin..GridZMax put its grid on those nodes, read from
psi*, with no interpolation.  The reference is MEQ at RefinementLevels = 2,
k = 4, MR-4c's reference, bounded at 2.4e-08 max (M-263), and it is NOT a race
entry.  So targets below about 1e-07 are not resolvable here and none is set.
NEITHER CODE IS THE TRUTH; this one is the finer of the two by two orders, and
geq's 513 is checked against it before anything is timed.

THE PROTOCOL is race_nke.py's: fresh processes, the arms interleaved, the
median of repeats, quiet waited for before anything runs, provenance on every
row.  Before each run the gate is the cores busy NOW, from /proc/stat over two
seconds, and not the load average, which after a multi-threaded run is mostly
that run (race_nke.py's own finding) and would add a minute's wait to each.  Both arms run at one thread and at the physical core count, with
OMP, MKL and OpenBLAS all set to it.

    cd /home/ian/projects/meq
    tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/race_geq.py

Results go to tools/geq-benchmark/runs/<stamp>/, which git ignores.
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
MEQ = os.path.abspath(os.path.join(HERE, "..", ".."))
GEQ_SRC = os.path.abspath(os.path.join(MEQ, "..", "geq", "src"))
VENV = os.path.join(MEQ, "tools", "freegs4e-benchmark", "venv", "bin", "python")
EXAMPLE = os.path.join(MEQ, "examples", "mirror-rotating.toml")
sys.path.insert(0, os.path.join(MEQ, "tools", "freegs4e-benchmark"))
from cores import physical_cores                                # noqa: E402

COARSE = 65
GRID = dict(rmin=0.05, rmax=0.63, zmin=-1.5, zmax=1.5)
TARGETS = (1e-4, 1e-5, 1e-6, 1e-7)


# ---------------------------------------------------------------------------
# The geq arm, in its own process
# ---------------------------------------------------------------------------

def geq_child(n, rtol, out):
    """One geq solve: export_rotating.py's solve(), exactly as the reference
    is generated, timed inside the process as well as around it."""
    sys.path.insert(0, HERE)
    import export_rotating as er
    er.SHAPE["name"] = "skew"
    er.R_MACH["value"] = er.SHAPES["skew"][1]
    omega_peak = er.omega_peak_for(2.0)
    _, _, rec = er.solve(n, "direct", omega_peak, rtol, 400, cache=None)
    np.savez(out, psi65=rec["psi65"], R65=rec["R65"], Z65=rec["Z65"],
             wall=rec["wall"], iterations=rec["iterations"],
             relchange=rec["relchange"], psi_axis=rec["psi_axis"],
             psi_bndry=rec["psi_bndry"])


def thread_env(threads):
    t = str(threads)
    return {**os.environ, "OMP_NUM_THREADS": t, "MKL_NUM_THREADS": t,
            "OPENBLAS_NUM_THREADS": t, "MPLBACKEND": "Agg"}


def run_geq(n, rtol, threads, out):
    env = thread_env(threads)
    env["PYTHONPATH"] = GEQ_SRC + os.pathsep + env.get("PYTHONPATH", "")
    started = time.perf_counter()
    done = subprocess.run([VENV, os.path.abspath(__file__), "--geq-child",
                           str(n), str(rtol), out],
                          capture_output=True, text=True, env=env, cwd=MEQ)
    seconds = time.perf_counter() - started
    if done.returncode != 0:
        return dict(ok=False, seconds=seconds, why=done.stderr[-400:])
    data = np.load(out)
    return dict(ok=True, seconds=seconds, solve=float(data["wall"]),
                iterations=int(data["iterations"]),
                relchange=float(data["relchange"]),
                psi=np.asarray(data["psi65"]), R=np.asarray(data["R65"]),
                Z=np.asarray(data["Z65"]), unknowns=n*n)


# ---------------------------------------------------------------------------
# The MEQ arm, through the driver
# ---------------------------------------------------------------------------

def meq_config(directory, mesh, level, order, prefix, plasma_size=0.1):
    """The shipped example with the resolution, the mesh's path and the output
    changed and nothing else.  Every substitution must match exactly once."""
    with open(EXAMPLE) as handle:
        text = handle.read()
    subs = [
        ('File = "examples/mirror-rotating.msh"',
         f'File = "{mesh}"\nRefinementLevels = {level}'),
        ("PolynomialDegree = 2", f"PolynomialDegree = {order}"),
        ("PlasmaSize = 0.1", f"PlasmaSize = {plasma_size}"),
        ('Directory = "."', f'Directory = "{directory}"'),
        ('Prefix = "mirror-rotating"', f'Prefix = "{prefix}"'),
        ("GridNR = 129", f"GridNR = {COARSE}\n"
                         f"GridRMin = {GRID['rmin']}\nGridRMax = {GRID['rmax']}"),
        ("GridNZ = 129", f"GridNZ = {COARSE}\n"
                         f"GridZMin = {GRID['zmin']}\nGridZMax = {GRID['zmax']}"),
    ]
    for old, new in subs:
        if text.count(old) != 1:
            raise SystemExit(f"{EXAMPLE}: '{old}' matched {text.count(old)} times")
        text = text.replace(old, new)
    path = os.path.join(directory, prefix + ".toml")
    with open(path, "w") as handle:
        handle.write(text)
    return path


def trace_unknowns(mesh_path, order):
    """faces x ( k + 1 ) from the written .mesh: Euler for a triangulation,
    ( 3 NE + boundary edges ) / 2, counted rather than estimated."""
    with open(mesh_path) as handle:
        lines = handle.read().split("\n")
    elements = int(lines[lines.index("elements") + 1])
    boundary = int(lines[lines.index("boundary") + 1])
    return elements, ((3*elements + boundary)//2)*(order + 1)


def run_meq(directory, mesh, level, order, threads, tag, plasma_size=0.1):
    os.makedirs(directory, exist_ok=True)
    config = meq_config(directory, mesh, level, order, tag, plasma_size)
    started = time.perf_counter()
    # --mesh-ready: the mesh was made once by meq-run --mesh-only from the same
    # [mesh.generate] block, which no variant changes.
    done = subprocess.run([os.path.join(MEQ, "build", "meq"), "--mesh-ready", config],
                          capture_output=True, text=True,
                          env=thread_env(threads), cwd=MEQ)
    seconds = time.perf_counter() - started
    if done.returncode != 0:
        return dict(ok=False, seconds=seconds, why=done.stdout[-400:] + done.stderr[-400:])
    m = re.search(r"MEQ: wall ([0-9.]+) s = setup ([0-9.]+) \+ solve ([0-9.]+) "
                  r"\+ output ([0-9.]+)", done.stdout)
    newton = re.search(r"converged in (\d+) Newton", done.stdout)
    import netCDF4
    with netCDF4.Dataset(os.path.join(directory, tag + ".nc")) as nc:
        psi = np.asarray(nc.variables["psi"][:], dtype=float).T   # ( R, Z ), as geq
        R = np.asarray(nc.variables["R"][:], dtype=float)
        Z = np.asarray(nc.variables["Z"][:], dtype=float)
        inside = np.asarray(nc.variables["inside"][:]).T
    elements, unknowns = trace_unknowns(os.path.join(directory, tag + ".mesh"), order)
    return dict(ok=True, seconds=seconds,
                solve=float(m.group(3)) if m else float("nan"),
                setup=float(m.group(2)) if m else float("nan"),
                output=float(m.group(4)) if m else float("nan"),
                iterations=int(newton.group(1)) if newton else -1,
                psi=psi, R=R, Z=Z, inside=inside, elements=elements,
                unknowns=unknowns)


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

def corner_mask():
    keep = np.ones((COARSE, COARSE), dtype=bool)
    for i in (0, COARSE - 1):
        for j in (0, COARSE - 1):
            keep[i, j] = False
    return keep


def error(psi, reference):
    keep = corner_mask()
    d = (psi - reference)[keep]
    return float(np.max(np.abs(d))), float(np.sqrt(np.mean(d*d)))


def load_average():
    with open("/proc/loadavg") as handle:
        return float(handle.read().split()[0])


def wait_for_quiet(threshold, settle, timeout):
    deadline = time.time() + timeout
    quiet_since = None
    while time.time() < deadline:
        if load_average() < threshold:
            quiet_since = quiet_since or time.time()
            if time.time() - quiet_since >= settle:
                return True
        else:
            quiet_since = None
        time.sleep(5)
    return False


def busy_cores(window=2.0):
    """Cores busy over the last @a window seconds, from /proc/stat: what is
    running NOW, where the one-minute load average is mostly the race's own
    previous run (race_nke.py's finding)."""
    def sample():
        with open("/proc/stat") as handle:
            fields = [float(x) for x in handle.readline().split()[1:]]
        idle = fields[3] + fields[4]
        return sum(fields), idle
    total0, idle0 = sample()
    time.sleep(window)
    total1, idle1 = sample()
    cores = os.cpu_count() or 1
    if total1 <= total0:
        return 0.0
    return cores*(1.0 - (idle1 - idle0)/(total1 - total0))


def wait_for_idle(threshold, timeout):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if busy_cores() < threshold:
            return True
    return False


def provenance():
    def sha(path):
        out = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True)
        dirty = subprocess.run(["git", "-C", path, "status", "--porcelain",
                                "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return (out.stdout.strip() or "?") + ("-dirty" if dirty else "")
    lib = "/home/ian/projects/mfem/install/lib/libmfem.a"
    version = subprocess.run([os.path.join(MEQ, "build", "meq"), "--version"],
                             capture_output=True, text=True).stdout.strip()
    return {"meq": sha(MEQ), "meq_version": version,
            "geq": sha(os.path.join(MEQ, "..", "geq")),
            "freegs4e": sha(os.path.join(MEQ, "..", "freegs4e")),
            "libmfem_mtime": time.strftime("%Y-%m-%d %H:%M",
                                           time.localtime(os.path.getmtime(lib)))
            if os.path.exists(lib) else "?",
            "physical_cores": physical_cores()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--geq-child", nargs=3, metavar=("N", "RTOL", "OUT"),
                    help=argparse.SUPPRESS)
    ap.add_argument("--out", default=None)
    ap.add_argument("--geq", default="65,129,257,513",
                    help="geq grids; each must contain the 65 grid")
    ap.add_argument("--geq-rtol", default="1e-6,2e-9",
                    help="geq's Picard rtol: 1e-6 is geq's own atol scale, "
                         "2e-9 the reference's")
    ap.add_argument("--levels", default="0,1,2",
                    help="MEQ's uniform refinements; L2 at k = 4 is the reference "
                         "and is never a race entry")
    ap.add_argument("--orders", default="1,2,3,4")
    ap.add_argument("--plasma-size", type=float, default=0.1,
                    help="[mesh.generate] PlasmaSize for MEQ's ladder, the element "
                         "size inside the plasma box; the shipped 0.1 by default. "
                         "The reference is always the shipped mesh at L2, k = 4")
    ap.add_argument("--threads", default=None,
                    help="default: 1 and the physical core count")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--repeats-slow", type=int, default=1,
                    help="repeats for any run whose first take exceeds 120 s")
    ap.add_argument("--quiet", type=float, default=0.5,
                    help="load average to wait under ONCE, before anything runs")
    ap.add_argument("--idle", type=float, default=0.5,
                    help="busy cores, from /proc/stat over 2 s, to wait under "
                         "before each run")
    args = ap.parse_args(argv)

    if args.geq_child:
        geq_child(int(args.geq_child[0]), float(args.geq_child[1]), args.geq_child[2])
        return

    print(f"waiting for the load average to stay under {args.quiet} for 60 s ...",
          flush=True)
    if not wait_for_quiet(args.quiet, 60, 3600):
        print("  (it never did; running anyway, and saying so)", flush=True)
    out = args.out or os.path.join(HERE, "runs", time.strftime("%Y%m%d-%H%M%S"))
    os.makedirs(out, exist_ok=True)
    threads = ([int(t) for t in args.threads.split(",")] if args.threads
               else sorted({1, physical_cores()}))
    record = dict(provenance=provenance(), threads=threads, runs=[])
    print(json.dumps(record["provenance"], indent=1), flush=True)

    # THE MESH, ONCE, UNTIMED: meq-run's generator on the shipped block.
    def make_mesh(path, size):
        probe = meq_config(out, path, 0, 2, f"mesh-probe-{size:g}", size)
        done = subprocess.run([os.path.join(MEQ, "build", "meq-run"), "--mesh-only", probe],
                              capture_output=True, text=True, cwd=MEQ)
        if done.returncode != 0:
            raise SystemExit("meq-run --mesh-only failed:\n" + done.stdout + done.stderr)
    mesh = os.path.join(out, "mirror.msh")
    make_mesh(mesh, 0.1)
    ladder = mesh
    if args.plasma_size != 0.1:
        ladder = os.path.join(out, f"mirror-ps{args.plasma_size:g}.msh")
        make_mesh(ladder, args.plasma_size)

    # THE REFERENCE, UNTIMED: MEQ at level 2, k = 4.
    print("\nreference: MEQ, RefinementLevels = 2, k = 4", flush=True)
    ref = run_meq(os.path.join(out, "reference"), mesh, 2, 4, max(threads), "reference")
    if not ref["ok"]:
        raise SystemExit("the reference did not solve:\n" + ref["why"])
    if not np.all(ref["inside"] == 1):
        raise SystemExit("a geq node is outside MEQ's domain")
    reference = ref["psi"]

    # geq's NODES ARE MEQ's NODES, checked rather than assumed, on one cheap run.
    first = run_geq(65, 1e-6, 1, os.path.join(out, "geq-nodes.npz"))
    if not first["ok"]:
        raise SystemExit("geq did not run:\n" + first["why"])
    node_gap = max(float(np.max(np.abs(first["R"] - ref["R"][:, None]))),
                   float(np.max(np.abs(first["Z"] - ref["Z"][None, :]))))
    print(f"geq's 65 x 65 nodes against MEQ's grid: {node_gap:.2e} m", flush=True)
    if node_gap > 1e-12:
        raise SystemExit("MEQ's grid is not on geq's nodes")

    entries = []
    for n in (int(x) for x in args.geq.split(",") if x):
        if (n - 1) % (COARSE - 1):
            raise SystemExit(f"geq grid {n} does not contain the {COARSE} grid")
        for rtol in (float(x) for x in args.geq_rtol.split(",")):
            for t in threads:
                entries.append(("geq", dict(n=n, rtol=rtol, threads=t)))
    for level in (int(x) for x in args.levels.split(",")):
        for order in (int(x) for x in args.orders.split(",")):
            if ( level, order ) == ( 2, 4 ) and args.plasma_size == 0.1:
                continue
            for t in threads:
                entries.append(("meq", dict(level=level, order=order, threads=t,
                                            plasma_size=args.plasma_size)))
    # INTERLEAVED: alternate the arms, cheap to expensive within each.
    geq_list = [e for e in entries if e[0] == "geq"]
    meq_list = [e for e in entries if e[0] == "meq"]
    order = []
    while geq_list or meq_list:
        if meq_list:
            order.append(meq_list.pop(0))
        if geq_list:
            order.append(geq_list.pop(0))

    slow = set()
    for repeat in range(args.repeats):
        for arm, cfg in order:
            key = json.dumps([arm, cfg], sort_keys=True)
            if repeat >= args.repeats_slow and key in slow:
                continue
            if not wait_for_idle(args.idle, 600):
                print("  (the machine never went idle; running anyway, and saying so)")
            load = busy_cores()
            if arm == "geq":
                r = run_geq(cfg["n"], cfg["rtol"], cfg["threads"],
                            os.path.join(out, f"geq-{cfg['n']}-{cfg['rtol']:g}-t{cfg['threads']}.npz"))
            else:
                tag = (f"meq-ps{cfg['plasma_size']:g}-L{cfg['level']}-k{cfg['order']}"
                       f"-t{cfg['threads']}-r{repeat}")
                r = run_meq(os.path.join(out, tag), ladder, cfg["level"], cfg["order"],
                            cfg["threads"], tag, cfg["plasma_size"])
            row = dict(arm=arm, repeat=repeat, load=load, **cfg)
            if r["ok"]:
                row.update(seconds=r["seconds"], solve=r["solve"],
                           iterations=r["iterations"], unknowns=r["unknowns"])
                row["max"], row["rms"] = error(r["psi"], reference)
                if arm == "meq":
                    row.update(elements=r["elements"], setup=r["setup"], output=r["output"])
                else:
                    row["relchange"] = r["relchange"]
                if r["seconds"] > 120.0:
                    slow.add(key)
            else:
                row.update(failed=r["why"])
            record["runs"].append(row)
            with open(os.path.join(out, "race.json"), "w") as handle:
                json.dump(record, handle, indent=1)
            label = (f"geq {cfg['n']:4d}^2 rtol {cfg['rtol']:.0e}" if arm == "geq"
                     else f"MEQ ps{cfg['plasma_size']:g} L{cfg['level']} k={cfg['order']}")
            print(f"  [{repeat}] {label:24s} t={cfg['threads']:2d}  "
                  + (f"{row['seconds']:8.2f} s (solve {row['solve']:8.2f})  "
                     f"max {row['max']:.3e}  rms {row['rms']:.3e}  busy {load:.2f}"
                     if r["ok"] else "FAILED " + r["why"][-120:]), flush=True)

    summarise(record, out)


def summarise(record, out):
    """Medians per configuration, then the two races."""
    groups = {}
    for row in record["runs"]:
        if "seconds" not in row:
            continue
        key = ((row["arm"], row["n"], row["rtol"], row["threads"]) if row["arm"] == "geq"
               else (row["arm"], row["level"], row["order"], row["threads"]))
        groups.setdefault(key, []).append(row)
    table = []
    for key, rows in sorted(groups.items(), key=lambda kv: str(kv[0])):
        first = rows[0]
        table.append(dict(key=key, arm=key[0], threads=key[3], takes=len(rows),
                          seconds=statistics.median(r["seconds"] for r in rows),
                          solve=statistics.median(r["solve"] for r in rows),
                          max=first["max"], rms=first["rms"],
                          unknowns=first["unknowns"],
                          ps=first.get("plasma_size", float("nan"))))
    lines = ["# MEQ against geq on the rotating mirror", "",
             "```", json.dumps(record["provenance"], indent=1), "```", "",
             "| arm | resolution | threads | takes | wall s | solve s | unknowns | max err | rms err |",
             "|---|---|---|---|---|---|---|---|---|"]
    for t in table:
        k = t["key"]
        res = (f"{k[1]}^2, rtol {k[2]:.0e}" if t["arm"] == "geq"
               else f"ps {t['ps']:g}, L{k[1]}, k = {k[2]}")
        lines.append(f"| {t['arm']} | {res} | {t['threads']} | {t['takes']} | "
                     f"{t['seconds']:.2f} | {t['solve']:.2f} | {t['unknowns']} | "
                     f"{t['max']:.2e} | {t['rms']:.2e} |")
    for threads in sorted({t["threads"] for t in table}):
        lines += ["", f"## Matched unknowns, {threads} thread(s): each geq grid against the "
                  "MEQ run whose faces x ( k + 1 ) is nearest its n^2", "",
                  "| geq | unknowns | s | max err | MEQ | unknowns | s | max err |",
                  "|---|---|---|---|---|---|---|---|"]
        meqs = [t for t in table if t["arm"] == "meq" and t["threads"] == threads]
        for g in (t for t in table if t["arm"] == "geq" and t["threads"] == threads):
            if not meqs:
                break
            m = min(meqs, key=lambda t: abs(np.log(t["unknowns"]/g["unknowns"])))
            lines.append(f"| {g['key'][1]}^2 rtol {g['key'][2]:.0e} | {g['unknowns']} | "
                         f"{g['seconds']:.2f} | {g['max']:.2e} | ps {m['ps']:g} L{m['key'][1]} k={m['key'][2]} | "
                         f"{m['unknowns']} | {m['seconds']:.2f} | {m['max']:.2e} |")
    for threads in sorted({t["threads"] for t in table}):
        lines += ["", f"## Time to accuracy, {threads} thread(s): the cheapest run "
                  "reaching each max-error target", "",
                  "| target | geq | s | MEQ | s | ratio geq/MEQ |", "|---|---|---|---|---|---|"]
        for target in TARGETS:
            best = {}
            for arm in ("geq", "meq"):
                hits = [t for t in table if t["arm"] == arm and t["threads"] == threads
                        and t["max"] <= target]
                best[arm] = min(hits, key=lambda t: t["seconds"]) if hits else None
            def cell(t):
                if t is None:
                    return "not reached", "-"
                k = t["key"]
                res = (f"{k[1]}^2 rtol {k[2]:.0e}" if t["arm"] == "geq"
                       else f"ps {t['ps']:g} L{k[1]} k={k[2]}")
                return res, f"{t['seconds']:.2f}"
            g, m = cell(best["geq"]), cell(best["meq"])
            ratio = (f"{best['geq']['seconds']/best['meq']['seconds']:.1f}"
                     if best["geq"] and best["meq"] else "-")
            lines.append(f"| {target:.0e} | {g[0]} | {g[1]} | {m[0]} | {m[1]} | {ratio} |")
    with open(os.path.join(out, "race.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
