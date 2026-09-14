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
MEQ_RUNGS = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1)]

# freegs4e grids.  Its multigrid wants 2^n + 1 and the convention is kept even
# for the direct solver, because a refinement study wants the coarse grid to be
# a SUBSET of the fine one -- which is what lets the two be differenced without
# interpolating either.
FGS_GRIDS = [129, 257, 513]

THREADS = str(os.cpu_count() or 1)


def environment():
    e = dict(os.environ)
    e["OMP_NUM_THREADS"] = THREADS
    e["MKL_NUM_THREADS"] = THREADS
    return e


# ---------------------------------------------------------------------------
# MEQ
# ---------------------------------------------------------------------------
def meq_toml(stem, degree, refine, scratch):
    """One rung's configuration, derived from the case's own file."""
    text = open(os.path.join(EXAMPLES, "%s.toml" % stem)).read()
    text = re.sub(r"^PolynomialDegree = .*$", "PolynomialDegree = %d" % degree,
                  text, flags=re.M)
    # RefinementLevels is not in the generated file, so it is inserted into
    # [mesh] rather than substituted.
    text = re.sub(r'^(File = "examples/%s\.msh")$' % re.escape(stem),
                  r"\1\nRefinementLevels = %d" % refine, text, flags=re.M)
    label = "%s-k%dr%d" % (stem, degree, refine)
    text = re.sub(r'^Prefix = ".*"$', 'Prefix = "%s"' % label, text, flags=re.M)
    text = re.sub(r'^(\[output\])$', r'\1\nDirectory = "%s"' % scratch,
                  text, flags=re.M)
    path = os.path.join(scratch, "%s.toml" % label)
    open(path, "w").write(text)
    return path, label


def run_meq(stem, degree, refine, scratch):
    path, label = meq_toml(stem, degree, refine, scratch)
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
                 ("psi_axis", "psi_boundary", "plasma_current")
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
