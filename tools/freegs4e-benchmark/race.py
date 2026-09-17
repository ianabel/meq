"""COST OF AN ACCURACY: MEQ against freegs4e, both cold, on the same machine.

WHAT A RACE BETWEEN TWO EQUILIBRIUM CODES CAN AND CANNOT MEAN.

perf.py is this file's predecessor and it compares MEQ solving a FIXED boundary
against freegs4e solving a free one -- strictly less work, and it says so at
length.  That comparison is gone: MEQ now solves the same free-boundary problem
these references do, from the same machine description, so a wall clock is
finally a comparable number.

WHAT IS HELD EQUAL, and each of these took real work to arrange:

  * BOTH CODES START COLD.  MEQ's [initialguess] is mkcoldguess.py -- the
    conductors at their given currents plus a parabolic current blob inscribed
    in the DESIGN plasma footprint -- and its [boundary.xpoint] seed is the
    X-point the control system was ASKED for.  Nothing in either comes from a
    converged equilibrium.  ( mkexactguess.py, which sums Green's functions
    over the reference's own converged Jtor, would put the answer in the
    starting position; it is kept for what it was written for and is not used
    here. )
  * BOTH RUN ON EVERY CORE.  freegs4e's inner loops are numpy and scipy over
    the whole grid, so they take the BLAS's threads; MEQ is given
    AssemblyMode::Threaded and the PARDISO trace solver, which are its defaults
    and which MEASUREMENTS.md M-78 measures as the fast pair.  OMP and MKL are
    both set to the core count, which is the configuration M-78 says to use
    once assembly is threaded.
  * BOTH ARE TIMED END TO END.  MEQ's number includes generating its mesh with
    gmsh, building the source, the solve, and writing four output formats;
    freegs4e's includes building its boundary matrix, the Picard loop and its
    own diagnostics.  Neither is a solve time and both are what a user waits.

WHAT IS NOT EQUAL AND CANNOT BE.  MEQ is C++ and freegs4e is Python; MEQ writes
four files and freegs4e writes one; freegs4e's conductors are FILAMENTS and
MEQ's are rectangles of uniform current density.  So read the accuracy-per-cost
SHAPE rather than any single ratio, which is the same rule perf.py states.

THE ACCURACY EACH CODE IS MEASURED AGAINST IS ITS OWN FINEST RUN, not the other
code's.  Measuring MEQ against freegs4e's grid charges MEQ for the conductor
model -- MEASUREMENTS.md M-87 puts that at 2.97e-01 relative L-infinity inside
two filament coils against 5.3e-04 outside them -- which is a modelling
difference and not a discretisation error.  So each code's own convergence is
reported, and their agreement at the finest is reported separately, as the floor
it is.

    python3 race.py <scratch-dir> [case-letters]

ONE CASE PER INVOCATION IF THE WALL CLOCKS ARE THE POINT, AND LET THE MACHINE
SETTLE BETWEEN THEM.  This script runs its cases back to back and its rungs back
to back, so on a multi-case run every case is timed on a machine the previous
ones have been holding at sixteen threads.  Measured on DIII-D's k2r0 -- one
binary, one configuration, psi_ax 3.759851e-01 either way -- 9.49 s inside a
three-case `race.py C F G` against 5.41 s settled and interleaved.  A factor of
two of pure machine state.  M-97 and M-102 never met it because they are
single-case, and MEASUREMENTS.md M-111 is the three-case table taken properly.

The ACCURACY columns are unaffected: those are answers, not clocks.

AND THE freegs4e ARM IS THE INVERSE SOLVE.  run_freegs() calls fgsref.py, which
builds through SyncConstrain, so this races MEQ's FORWARD problem against
freegs4e's INVERSE one -- not a like-for-like comparison of solvers, and the
reason freegs4e converges here from cold on all seven while its forward arm
fails all fourteen.  M-112.  It is not a defect of this script so much as the
only comparison available: there is no case where both codes solve the same
problem from cold.
"""

import json
import os
import re
import subprocess
import sys
import time

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conductors

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MEQ_RUN = os.path.join(ROOT, "build", "meq-run")
EXAMPLES = os.path.join(ROOT, "examples")
VENV = os.path.join(HERE, "venv", "bin", "python")

# reference stem -> MEQ example stem
CASES = [
    ("A", "A_testtokamak_classic", "machine-a-testtokamak"),
    ("B", "B_testtokamak_peaked_ffdom", "machine-b-ffprime"),
    ("C", "C_mast_spherical", "machine-c-mast"),
    ("D", "D_tcv_conventional", "machine-d-tcv"),
    ("E", "E_testtokamak_diamagnetic", "machine-e-diamagnetic"),
    ("F", "F_diiid_conventional", "machine-f-diiid"),
    ("G", "G_mastu_simple", "machine-g-mastu"),
]

# ( polynomial degree, uniform refinement levels ).  RefinementLevels rather
# than a coarser [mesh.generate] Size, so every rung of one case shares one
# gmsh mesh and the comparison is a pure h-refinement of it.
# ( degree, uniform refinement levels [, adaptive cycles ] ).
# RefinementLevels rather than a coarser [mesh.generate] Size, so every rung of
# one case shares one gmsh mesh and the h-sweep is a pure refinement of it; the
# adaptive rows use the driver's own loop, which warm starts each cycle from
# the last and refines where the residual estimator says to.
MEQ_RUNGS = [(1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1),
             (2, 0, 3), (3, 0, 3)]

# freegs4e grids.  Its multigrid wants 2^n + 1 and the convention is kept even
# for the direct solver, because a refinement study wants the coarse grid to be
# a SUBSET of the fine one -- which is what lets the two be differenced without
# interpolating either.
# TIMED COLD, and only as far as 257^2.  A cold 513^2 is an afternoon -- the
# boundary condition alone scales about 6.9x per doubling -- and M-88 measures
# what it would buy: freegs4e converges at about FIRST order here, so 513^2 is
# a factor of two on 257^2 for five to seven times the cost.  The 513^2 the
# errors below are measured against is the SEEDED one from that ladder, which
# is the same equilibrium however it was reached.
FGS_GRIDS = [129, 257]

# Where the seeded 513^2 ladder of M-88 landed. Accuracy does not care how a
# reference was reached; only the TIMED rows have to be cold.
TRUTH = os.environ.get("MEQ_RACE_TRUTH", "")

THREADS = str(os.cpu_count() or 1)


def environment():
    e = dict(os.environ)
    e["OMP_NUM_THREADS"] = THREADS
    e["MKL_NUM_THREADS"] = THREADS
    return e


# ---------------------------------------------------------------------------
# MEQ
# ---------------------------------------------------------------------------
def meq_toml(stem, degree, refine, scratch, sample=513, adaptive=0):
    """One rung's configuration, derived from the case's own file."""
    text = open(os.path.join(EXAMPLES, "%s.toml" % stem)).read()
    text = re.sub(r"^PolynomialDegree = .*$", "PolynomialDegree = %d" % degree,
                  text, flags=re.M)
    # RefinementLevels is not in the generated file, so it is inserted into
    # [mesh] rather than substituted.
    text = re.sub(r'^(File = "examples/%s\.msh")$' % re.escape(stem),
                  r"\1\nRefinementLevels = %d" % refine, text, flags=re.M)
    # THE SAMPLING GRID IS THE COMPARISON'S AND NOT THE RUN'S. MEQ's .nc
    # covers the whole half-disc, so 129^2 over a radius of 3.7 m is 0.029 m in
    # the plasma against the 513^2 reference's 0.0035 m -- the sampling would
    # be the error rather than the solve. It is raised to match, which inflates
    # the `output` column and nothing else; that column is reported separately
    # for exactly this reason.
    text = re.sub(r"^GridNR = .*$", "GridNR = %d" % sample, text, flags=re.M)
    text = re.sub(r"^GridNZ = .*$", "GridNZ = %d" % sample, text, flags=re.M)

    label = "%s-k%dr%d" % (stem, degree, refine)
    if adaptive:
        label += "a%d" % adaptive
        text += ("\n[adaptivity]\nEnabled = true\nMaxIterations = %d\n"
                 "Theta = 0.6\nTargetError = 1.0e-9\n" % adaptive)
    text = re.sub(r'^Prefix = ".*"$', 'Prefix = "%s"' % label, text, flags=re.M)
    text = re.sub(r'^(\[output\])$', r'\1\nDirectory = "%s"' % scratch,
                  text, flags=re.M)
    path = os.path.join(scratch, "%s.toml" % label)
    open(path, "w").write(text)
    return path, label


def run_meq(stem, degree, refine, scratch, sample=513, adaptive=0):
    path, label = meq_toml(stem, degree, refine, scratch, sample, adaptive)
    started = time.perf_counter()
    done = subprocess.run([MEQ_RUN, path], capture_output=True, text=True,
                          cwd=ROOT, env=environment(), timeout=7200)
    wall = time.perf_counter() - started
    row = dict(label=label, degree=degree, refine=refine, wall=wall,
               ok=done.returncode == 0, log=done.stdout + done.stderr)
    if not row["ok"]:
        return row

    phases = re.search(r"MEQ: wall ([\d.]+) s = setup ([\d.]+) \+ solve "
                       r"([\d.]+) \+ output ([\d.]+)", done.stdout)
    if phases:
        row.update(driver=float(phases.group(1)), setup=float(phases.group(2)),
                   solve=float(phases.group(3)), output=float(phases.group(4)))
    counts = re.search(r"converged in (\d+) Newton iterations on (\d+) "
                       r"elements, degree (\d+)", done.stdout)
    if counts:
        row.update(iterations=int(counts.group(1)),
                   elements=int(counts.group(2)))
        # P_k on a triangle carries ( k+1 )( k+2 )/2 coefficients; the potential
        # is one such field and the flux is two. The trace is what is actually
        # solved after hybridization, but the total is the honest measure of the
        # discrete problem's size against a grid's point count.
        per = (degree + 1)*(degree + 2)//2
        row["dofs"] = row["elements"]*per*3
    row["nc"] = os.path.join(scratch, "%s.nc" % label)
    return row


def meq_field(path):
    """psi on the run's own ( R, Z ) grid, with the nodes it does not speak
    for masked. Every rung of one case shares a grid, since the disc and
    GridNR/GridNZ are the same file's."""
    with Dataset(path) as ds:
        psi = np.array(ds["psi"][:], float)
        inside = np.array(ds["inside"][:]).astype(bool)
        extrapolated = (np.array(ds["extrapolated"][:]).astype(bool)
                        if "extrapolated" in ds.variables
                        else np.zeros_like(inside))
        return (np.array(ds["R"][:], float), np.array(ds["Z"][:], float),
                psi, inside & ~extrapolated,
                {k: float(getattr(ds, k)) for k in
                 ("psi_axis", "psi_boundary", "plasma_current",
                  "xpoint_r", "xpoint_z")
                 if hasattr(ds, k)})


# ---------------------------------------------------------------------------
# freegs4e
# ---------------------------------------------------------------------------
def run_freegs(ref, nx, scratch):
    """One reference resolution, COLD. Returns its record and where it went."""
    out = os.path.join(scratch, "fgs")
    env = environment()
    env["FGSREF_OUT"] = out
    env["PYTHONPATH"] = os.path.abspath(os.path.join(ROOT, "..", "freegs4e"))
    where = os.path.join(out, "n%d" % nx)
    npz = os.path.join(where, "%s.npz" % ref)
    if not os.path.exists(npz):
        started = time.perf_counter()
        done = subprocess.run([VENV, os.path.join(HERE, "fgsref.py"),
                               "--nx=%d" % nx, ref[0]],
                              capture_output=True, text=True, cwd=HERE,
                              env=env, timeout=14400)
        wall = time.perf_counter() - started
        if not os.path.exists(npz):
            return dict(nx=nx, ok=False, wall=wall,
                        log=done.stdout[-4000:] + done.stderr[-4000:])
        open(os.path.join(where, "%s.race" % ref), "w").write("%.6f\n" % wall)
    wall = float(open(os.path.join(where, "%s.race" % ref)).read())
    rec = json.load(open(os.path.join(where, "%s.json" % ref)))
    return dict(nx=nx, ok=True, wall=wall, npz=npz, rec=rec,
                inner=float(rec["wall_time_s"]))


def freegs_field(npz):
    d = np.load(npz, allow_pickle=True)
    return (np.array(d["R"], float), np.array(d["Z"], float),
            np.array(d["psi"], float))       # psi[ iR, iZ ]


# ---------------------------------------------------------------------------
# errors
# ---------------------------------------------------------------------------
def relative(a, b, mask):
    """|a - b| over |b|, in L2 and L-infinity, over `mask`."""
    da = (a - b)[mask]
    scale = np.sqrt(np.mean(b[mask]**2))
    return (np.sqrt(np.mean(da**2))/scale,
            np.max(np.abs(da))/np.max(np.abs(b[mask])))


def conductor_mask(R, Z, flat, collar):
    """True where a node is well away from every conductor.

    M-87: the whole of a 2.97e-01 relative L-infinity between the two codes on
    case A lives inside two coils freegs4e models as filaments and MEQ as
    rectangles.  That is a modelling difference and it does not refine away in
    either code, so it has no place in a convergence study of either.
    """
    keep = np.ones((Z.size, R.size), bool)
    RR, ZZ = np.meshgrid(R, Z)
    for c in flat:
        keep &= ~((RR > c["R"] - c["half_width"] - collar)
                  & (RR < c["R"] + c["half_width"] + collar)
                  & (ZZ > c["Z"] - c["half_height"] - collar)
                  & (ZZ < c["Z"] + c["half_height"] + collar))
    return keep


# ---------------------------------------------------------------------------
# the table
# ---------------------------------------------------------------------------
def sweep(tag, ref, stem, scratch, rungs, grids, sample=513, collar=0.05):
    fine = None
    print("\n" + "="*78)
    print("  CASE %s -- %s" % (tag, ref))
    print("="*78)

    # ---- freegs4e ------------------------------------------------------
    rows = []
    for nx in grids:
        row = run_freegs(ref, nx, scratch)
        rows.append(row)
        if not row["ok"]:
            print("  freegs4e %d^2 FAILED" % nx)
            continue
        fine = row["npz"]
    print("\n  freegs4e, cold at each grid")
    print("    %6s %10s %9s %15s %15s" %
          ("grid", "wall/s", "picard", "psi_ax", "psi_bnd"))
    # THE REFERENCE'S OWN ACTIVE NULL, so MEQ's column has something to be read
    # against rather than only against the previous rung. Taken from the
    # FINEST grid that ran, since that is what the errors are measured to.
    # THE ACTIVE ONE IS THE ONE AT psi_bndry, not the first in the list.
    # fgsref records every saddle it found -- 23 of them on MAST-U, most of
    # them on the centre column -- so picking by index would name a null the
    # plasma has nothing to do with.
    for row in reversed(rows):
        if not row["ok"]:
            continue
        xs = row["rec"].get("xpoints") or []
        if not xs:
            continue
        b = row["rec"]["psi_bndry"]
        r, z, _ = min(xs, key=lambda t: abs(t[2] - b))
        print("      the reference's ACTIVE X-point at %d^2: ( %.4f, %.4f )"
              "   [ %d saddles found ]" % (row["nx"], r, z, len(xs)))
        break
    for row in rows:
        if not row["ok"]:
            continue
        r = row["rec"]
        print("    %6d %10.1f %9d %15.9e %15.9e" %
              (row["nx"], row["wall"], r.get("stage2_iterations", 0),
               r["psi_axis"], r["psi_bndry"]))

    if TRUTH:
        candidate = os.path.join(TRUTH, "%s.npz" % ref)
        if os.path.exists(candidate):
            fine = candidate
    if fine is None:
        print("  no reference to measure against")
        return
    print("\n  measured against %s" % fine)

    truth = json.load(open(fine.replace(".npz", ".json")))
    flat = conductors.from_npz(np.load(os.path.join(HERE, "%s.npz" % ref),
                                       allow_pickle=True))
    boxes = [(c["R"], c["Z"], c["half_width"] + collar,
              c["half_height"] + collar) for c in flat]

    # ---- MEQ -----------------------------------------------------------
    from compare import compare
    print("\n  MEQ, cold")
    # THE X-POINT IS REPORTED BECAUSE A RUNG CAN SELECT THE WRONG NULL AND
    # LOOK FINE. Measured on G MAST-U: k2r0 converges cleanly, reports a
    # plausible psi_ax and has picked the X-point at ( 0.425, -1.932 ) where
    # k3r0 picks ( 0.621, -1.130 ) -- a different equilibrium on a super-X
    # geometry with several candidate nulls. Nothing in the other columns says
    # so; psi_ax alone reads 4.79e-02 against 1.29e-04 and could be excused as
    # coarseness. MEASUREMENTS.md M-111.
    print("    %-10s %9s %9s %7s %8s %8s %8s %11s %11s %10s %10s %19s" %
          ("rung", "elements", "dofs", "newton", "setup", "solve", "wall",
           "rel L2", "no coils", "psi_ax", "psi_bnd", "X-point"))
    for rung in rungs:
        row = run_meq(stem, rung[0], rung[1], scratch, sample=sample,
                      adaptive=rung[2] if len(rung) > 2 else 0)
        if not row["ok"] or "nc" not in row or not os.path.exists(row["nc"]):
            print("    %-10s FAILED" % row["label"].split("-")[-1])
            tail = [l for l in row["log"].splitlines()
                    if "MEQ:" in l and ("not converge" in l or "error" in l)]
            for line in tail[:2]:
                print("      " + line.strip())
            continue
        try:
            got = compare(fine, row["nc"], None, free_boundary=True,
                          boxes=boxes)
        except Exception as exc:
            print("    %-10s compare failed: %r" % (row["label"], exc))
            continue
        with Dataset(row["nc"]) as ds:
            psi_ax = float(getattr(ds, "psi_axis"))
            psi_bnd = float(getattr(ds, "psi_boundary"))
            xr = float(getattr(ds, "xpoint_r", float("nan")))
            xz = float(getattr(ds, "xpoint_z", float("nan")))
        ax = abs(psi_ax - truth["psi_axis"])/abs(truth["psi_axis"])
        bn = abs(psi_bnd - truth["psi_bndry"])/max(abs(truth["psi_bndry"]),
                                                   1e-300)
        keep = got["kept"]["rel_l2"] if got["kept"] else float("nan")
        print("    %-10s %9d %9d %7d %8.2f %8.2f %8.2f %11.3e %11.3e %10.2e "
              "%10.2e  ( %7.4f, %7.4f )" %
              (row["label"].split(stem + "-")[-1], row.get("elements", -1),
               row.get("dofs", -1), row.get("iterations", -1),
               row.get("setup", float("nan")), row.get("solve", float("nan")),
               row["wall"], got["rel_l2"], keep, ax, bn, xr, xz))


if __name__ == "__main__":
    scratch = sys.argv[1]
    want = [a.upper() for a in sys.argv[2:]] or [c[0] for c in CASES]
    os.makedirs(scratch, exist_ok=True)
    print("\n  A RACE, BOTH CODES COLD, %s threads each" % THREADS)
    print("  MEQ: PARDISO trace solver and threaded assembly, which are its")
    print("  defaults; its wall includes gmsh, the solve and four output")
    print("  formats. freegs4e's includes its boundary matrix, the Picard")
    print("  loop and its own diagnostics. Neither is a solve time.")
    for tag, ref, stem in CASES:
        if tag in want:
            sweep(tag, ref, stem, scratch, MEQ_RUNGS, FGS_GRIDS)
