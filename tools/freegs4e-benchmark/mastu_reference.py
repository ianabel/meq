#!/usr/bin/env python
"""THE MAST-U FORWARD REFERENCE, OVER A RESOLUTION LADDER.

freegsnke's `example02 - static_forward_solve_MASTU.ipynb` with the notebook
removed and a grid sweep put round it.  It is the benchmark's first SAME-MODE
case: both arms solve a FORWARD problem with the coil currents prescribed,
where `tools/MAST-U-freegsnke.py` is an INVERSE solve and so cannot be compared
against a MEQ run at all (memory: *compare in the same mode MEQ solves in*).

WHY THIS CASE AND NOT THE INVERSE ONE
  * `constrain=None` -- no control loop, so the currents are inputs on both
    sides rather than something one code solves for and the other is handed.
  * the 138 passive structures carry EXACTLY ZERO current here.  Checked, not
    assumed: the currents pickle has 150 entries and only the 12 active coils
    are non-zero.  So MEQ needs no passive-conductor model to run this, which
    is what made the inverse case unportable.
  * the profile is `ConstrainPaxisIp`, which is a prescribed p' and ff' family
    plus a prescribed I_p -- exactly MEQ's `[source] Type = "mhd"` with
    `PlasmaCurrent`.

THE LADDER is 2^n + 1 on both axes so that each grid is a point-for-point
subset of the next, which is what a reference refinement study needs and what
`ladder.sh` already requires of the freegs4e one.  The aspect follows the
domain: R spans 1.9 m and Z spans 4.4 m, so ny = 2 nx - 1 keeps the cells
roughly square.

WHAT IT WRITES, per level:
  * `mastu-fwd-<nx>x<ny>.npz`, in the shape `compare.py` already reads, so a
    MEQ `.nc` can be differenced against it with the existing tool.
  * a row in `mastu-fwd-timings.json` carrying wall clock, peak RSS, the
    Picard/Newton-Krylov split and the provenance of both trees.

READ THE TIMINGS AS ONE MACHINE'S ONE RUN.  CLAUDE.md records a 490/540 s
spread on identical code, and this session added that an ITERATION COUNT is a
measurement against one library too.  `--repeats` defaults to 3 and the median
is reported, because a median is right for a leg share where a best-of is right
for a minimum wall time.
"""
import argparse, contextlib, io as _io, json, os, pickle, resource, subprocess, sys, time
import numpy as np

MACHINE = "/home/ian/projects/freegsnke/machine_configs/MAST-U"
CURRENTS = "/home/ian/projects/freegsnke/examples/data/simple_diverted_currents_PaxisIp.pk"

# example02's profile, verbatim.  Changing any of these makes a different
# reference and the npz files must be regenerated together.
PROFILE = dict(paxis=8e3, Ip=6e5, fvac=0.5, alpha_m=1.8, alpha_n=1.2)
DOMAIN = dict(Rmin=0.1, Rmax=2.0, Zmin=-2.2, Zmax=2.2)
LADDER = [(33, 65), (65, 129), (129, 257), (257, 513)]


def provenance():
    """Both trees' SHAs and the MFEM install's mtime.

    Recorded because a reference generated against one freegs4e is not the same
    reference as one generated against another, and the pin warning on this
    venv -- freegsnke asks for freegs4e >= 0.14 and gets 0.12 -- makes that a
    live possibility rather than a formality.
    """
    def sha(path):
        try:
            return subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                                  capture_output=True, text=True).stdout.strip() or "?"
        except Exception:
            return "?"
    out = {"freegsnke": sha("/home/ian/projects/freegsnke"),
           "freegs4e": sha("/home/ian/projects/freegs4e"),
           "meq": sha("/home/ian/projects/meq")}
    try:
        out["libmfem_mtime"] = time.strftime(
            "%Y-%m-%d %H:%M",
            time.localtime(os.path.getmtime("/home/ian/projects/mfem/install/lib/libmfem.a")))
    except OSError:
        out["libmfem_mtime"] = "?"
    return out


def solve_at(nx, ny, tolerance, verbose=False):
    """One forward solve.  Returns (eq, profiles, seconds, log)."""
    from freegsnke import build_machine, equilibrium_update, GSstaticsolver
    from freegsnke.jtor_update import ConstrainPaxisIp

    tokamak = build_machine.tokamak(
        active_coils_path=f"{MACHINE}/MAST-U_like_active_coils.pickle",
        passive_coils_path=f"{MACHINE}/MAST-U_like_passive_coils.pickle",
        limiter_path=f"{MACHINE}/MAST-U_like_limiter.pickle",
        wall_path=f"{MACHINE}/MAST-U_like_wall.pickle")
    eq = equilibrium_update.Equilibrium(tokamak=tokamak, nx=nx, ny=ny, **DOMAIN)
    profiles = ConstrainPaxisIp(eq=eq, **PROFILE)
    solver = GSstaticsolver.NKGSsolver(eq, gs_operator_order=4)

    with open(CURRENTS, "rb") as handle:
        currents = pickle.load(handle)
    for label, value in currents.items():
        eq.tokamak.set_coil_current(coil_label=label, current_value=value)

    # THE CLOCK STARTS AFTER THE MACHINE IS BUILT, deliberately.  Building the
    # Green's functions is a per-machine cost that a MEQ run pays as meshing,
    # outside its own solve timer -- so including it here would compare a
    # solve against a solve-plus-setup.
    # THE LOG IS CAPTURED RATHER THAN READ OFF THE SOLVER, because there is
    # nothing to read: forward_solve() builds `log` as a LOCAL, prints it and
    # clears it every iteration.  An attribute that looked like it held the
    # history would have reported zero Picard and zero Newton-Krylov steps
    # forever, which is the shape of instrument error this tree records most.
    captured = _io.StringIO()
    started = time.perf_counter()
    if verbose:
        with contextlib.redirect_stdout(captured):
            solver.solve(eq=eq, profiles=profiles, constrain=None,
                         target_relative_tolerance=tolerance, verbose=True)
    else:
        solver.solve(eq=eq, profiles=profiles, constrain=None,
                     target_relative_tolerance=tolerance, verbose=False)
    seconds = time.perf_counter() - started
    return eq, profiles, seconds, solver, captured.getvalue()


def write_npz(path, eq, profiles, nx, ny):
    """The interchange shape `compare.py` reads, plus what this case adds."""
    psi = np.asarray(eq.psi())
    np.savez_compressed(
        path,
        R=np.asarray(eq.R)[:, 0].astype(float),
        Z=np.asarray(eq.Z)[0, :].astype(float),
        psi=psi.astype(float),
        psi_axis=np.float64(eq.psi_axis),
        psi_bndry=np.float64(eq.psi_bndry),
        Ip=np.float64(eq.plasmaCurrent()),
        Raxis=np.float64(eq.Rmagnetic()),
        Zaxis=np.float64(eq.Zmagnetic()),
        nx=np.int64(nx), ny=np.int64(ny),
        machine=np.array("MAST-U (freegsnke, forward)"),
        conductor_model=np.array("filament, actives only; passives at zero"),
        fvac=np.float64(PROFILE["fvac"]),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=".", help="where to write the npz and the timings")
    ap.add_argument("--tolerance", type=float, default=1e-9)
    ap.add_argument("--repeats", type=int, default=3,
                    help="timed repeats per level; the MEDIAN is reported")
    ap.add_argument("--levels", default="",
                    help="comma-separated nx values to restrict the ladder to")
    args = ap.parse_args()

    wanted = {int(v) for v in args.levels.split(",") if v.strip()} if args.levels else None
    ladder = [g for g in LADDER if wanted is None or g[0] in wanted]

    os.makedirs(args.out, exist_ok=True)
    rows = []
    print(f"{'grid':>12} {'median s':>10} {'peak MB':>9} {'picard':>7} "
          f"{'NK':>4} {'rel':>10} {'psi_axis':>14} {'Ip':>12}")
    for nx, ny in ladder:
        times = []
        log = ""
        eq = profiles = solver = None
        for repeat in range(max(1, args.repeats)):
            # VERBOSE ON THE FIRST REPEAT ONLY -- the Picard/NK split is a
            # property of the solve and printing it three times would only add
            # I/O to the two timed repeats that follow.
            eq, profiles, seconds, solver, text = solve_at(
                nx, ny, args.tolerance, verbose=(repeat == 0))
            if repeat == 0:
                log = text
            times.append(seconds)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0
        median = float(np.median(times))

        # The split, counted from the captured verbose output.
        picard = log.count("Picard iteration")
        krylov = log.count("Newton-Krylov iteration")

        stem = os.path.join(args.out, f"mastu-fwd-{nx}x{ny}")
        write_npz(stem + ".npz", eq, profiles, nx, ny)
        rows.append(dict(nx=nx, ny=ny, seconds=median, all_seconds=times,
                         peak_mb=peak, picard=picard, krylov=krylov,
                         relative=float(getattr(solver, "relative_change", float("nan"))),
                         psi_axis=float(eq.psi_axis), psi_bndry=float(eq.psi_bndry),
                         Ip=float(eq.plasmaCurrent())))
        print(f"{nx:5d}x{ny:<6d} {median:10.2f} {peak:9.0f} {picard:7d} "
              f"{krylov:4d} {rows[-1]['relative']:10.2e} "
              f"{eq.psi_axis:14.8e} {eq.plasmaCurrent():12.5e}")
        sys.stdout.flush()

    with open(os.path.join(args.out, "mastu-fwd-timings.json"), "w") as handle:
        json.dump(dict(provenance=provenance(), profile=PROFILE, domain=DOMAIN,
                       tolerance=args.tolerance, repeats=args.repeats,
                       rows=rows), handle, indent=2)
    print(f"\nwrote {len(rows)} level(s) and mastu-fwd-timings.json to {args.out}")


if __name__ == "__main__":
    main()
