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
    python3 race.py <scratch-dir> F --filament      # the matched race

THE RACE TO RUN ON A QUIET MACHINE, and it is the `--filament` one:

    cd tools/freegs4e-benchmark
    nice -n 10 ./venv/bin/python race.py /tmp/race-fil F --filament

**COLD-STARTED FORWARD FILAMENTS, BOTH SIDES.**  freegs4e's DIII-D reference
carries its 18 conductors as point filaments and MEQ now carries the same 18 at
the same positions -- `[conductors] Model = "filament"`, the conductors out of
the mesh and evaluated analytically -- so for the first time the two arms are
the SAME MACHINE rather than two approximations to one.  Both start cold: MEQ
from `mkcoldguess.py --remainder`, which is the design footprint carrying the
target current and no converged field anywhere in it.

What makes it worth the machine time is that the matched model breaks below a
floor this project has measured twice and called irreducible.
MEASUREMENTS.md M-111 reads DIII-D flat at 5.785e-03 over a 16x range in dofs;
M-154 reads the filament arm at **2.031e-03 against the meshed route's
4.593e-03, on 1510 elements against 4848**.

**AND THIS DOCSTRING PREDICTED "A RACE WHERE ONE ARM IS BOTH FASTER AND MORE
ACCURATE", WHICH IS HALF RIGHT AND THE WRONG HALF.**  Run -- M-155 -- MEQ is
**6.2x to 19.8x SLOWER**: freegs4e takes 4.9 s at 129^2 and 12.9 s at 257^2
where MEQ's converging rungs take 80 to 256 s.  The accuracy half stands, and
the right shape for the expectation was always the one M-154 actually
measured -- fewer elements for a better answer -- which is a statement about
DOFS and says nothing whatever about seconds.  The two arms do different
amounts of work per unknown and nothing here had compared that.

The race is still the one to run, for the accuracy and because a wall clock is
the only thing that can say how much the accuracy costs.  It is not a race MEQ
is expected to win.

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
# --shaped: RACE THE SAME CONDUCTOR MODEL ON BOTH SIDES.
#
# Without it the freegs4e arm is a FILAMENT reference and MEQ's conductors are
# meshed rectangles carrying a uniform current density, so the two arms are not
# solving the same problem and part of every accuracy figure is that rather than
# either solver. tools/freegs4e-benchmark/conductor_model.py measures it on
# DIII-D with no solver anywhere: 6.081e-03 globally against M-111's own
# irreducible 5.785e-03, a 1.05x match, and freegs4e's ShapedCoil against MEQ's
# exact rectangle is 6.580e-04 -- NINE TIMES smaller, and below the error the
# race is trying to see.
#
# BOTH HALVES HAVE TO MOVE TOGETHER. The reference is regenerated with
# `fgsref.py --shaped`, AND the MEQ config becomes the `-shaped` one, because
# the guess and the profile tables are derived FROM the reference and a MEQ run
# built from the filament conversion compared against a shaped reference would
# be measuring the conversion. MEASUREMENTS.md says so in terms and re-runs
# make_diverted_case.py from the shaped npz for exactly this reason.
#
# IT IS NOT A FREE UPGRADE AND THE TREE ALREADY MEASURED THAT. Of the seven
# machines run both ways, five behave identically, D TCV goes from FAILED to
# CONVERGED and G MAST-U goes the other way. That table predates M-109's axis
# guard and M-111 races C successfully where it reads FAILED, so it wants
# re-taking -- which is the point of this flag.
SHAPED = "--shaped" in sys.argv
if SHAPED:
    sys.argv.remove("--shaped")

# --filament: THE OTHER WAY TO MATCH THE CONDUCTOR MODEL, AND THE CHEAPER ONE.
#
# `--shaped` above moves the REFERENCE to MEQ's rectangles.  This moves MEQ to
# the reference's FILAMENTS, and it is the direction `COIL-SUBTRACTION-PLAN.md`
# exists for: `[conductors] Model = "filament"` takes the conductors out of the
# mesh and evaluates them analytically at exactly the points freegs4e's own
# `Coil` sits at, so the two arms are the same machine rather than two
# approximations to one.
#
# IT IS ALSO A SMALLER MESH, which is the point of the plan rather than a side
# effect: the conductors are not meshed, not graded to, and not fragmented.
# Measured on DIII-D, MEASUREMENTS.md M-154 -- 5054 triangles to 1702, and the
# field-wide error against the filament reference 4.593e-03 meshed against
# 2.031e-03 subtracted, which is BELOW M-111's 5.785e-03 floor.
#
# FOUR THINGS MOVE TOGETHER AND ALL FOUR ARE REQUIRED.  `[mesh.generate]
# CoilSize` is refused beside a subtracting model and has to go; the mesh needs
# its own filename so the two arms do not fight over `<mesh>.meq-mesh`; the
# guess has to be built as psi_p rather than converted from a total, because
# psi_c diverges at a filament and any stored psi cannot; and the table needs
# `[conductors] Model = "filament"` appended.  Leaving out the third is the
# quiet one -- it converges, to a different equilibrium.
FILAMENT = "--filament" in sys.argv
if FILAMENT:
    sys.argv.remove("--filament")
if FILAMENT and SHAPED:
    raise SystemExit("--shaped and --filament are the two DIRECTIONS of the "
                     "same matching and asking for both is asking for neither: "
                     "one moves the reference to MEQ's rectangles and the other "
                     "moves MEQ to the reference's filaments")


def shaped_ref(ref):
    return ref + "_shaped" if SHAPED else ref


def shaped_stem(stem):
    return stem + "-shaped" if SHAPED else stem


CASES = [
    ("A", "A_testtokamak_classic", "machine-a-testtokamak"),
    ("B", "B_testtokamak_peaked_ffdom", "machine-b-ffprime"),
    ("C", "C_mast_spherical", "machine-c-mast"),
    ("D", "D_tcv_conventional", "machine-d-tcv"),
    ("E", "E_testtokamak_diamagnetic", "machine-e-diamagnetic"),
    ("F", "F_diiid_conventional", "machine-f-diiid"),
    ("G", "G_mastu_simple", "machine-g-mastu"),
    # THE ONE LIMITED MACHINE, AND IT IS HERE BECAUSE IT IS THE ONLY CASE
    # EVERY COMPETITOR CAN TAKE.  A to G are diverted, which rules out DESC's
    # free-boundary objective -- its LCFS is a truncated Fourier series and a
    # separatrix has an X-point corner -- and TSC's own deck needs a rectangle
    # holding the plasma while excluding the coils, which DIII-D's inboard F
    # coils at R = 0.8608 very nearly forbid.  This machine has neither
    # problem, so it is the common yardstick.
    ("H", "H_limited_circular", "limited-tokamak"),
]

#: stem -> a SHIPPED filament configuration, for cases where DERIVING one from
#: the meshed file does not work.
#:
#: meq_toml()'s four edits below assume the case carries `[mesh.generate]`:
#: edit (2) points the mesh at a scratch path that the generator then has to
#: make, and filament_guess() reads `Radius` out of the same block to size the
#: guess.  examples/limited-tokamak.toml names a COMMITTED mesh and has no
#: such block -- tools/freegs4e-benchmark/README.md records that the shipped
#: .msh is the fixture and that the generator command does NOT reproduce it,
#: 1775 triangles against 1807 with MEQ's own scatter at that resolution
#: 1.5-2% in psi_ax.  So the filament arm is a file rather than a derivation.
FILAMENT_FILE = {
	"limited-tokamak": "limited-tokamak-filament",
}

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

# --grids 129,257,513: THE REFERENCE RESOLUTIONS, AND THE FINEST IS WHAT MEQ IS
# SCORED AGAINST.  It is a flag because a case's MEQ configuration is posed
# against ONE reference grid and has to be scored against that one: freegs4e's
# coil currents and profile amplitudes are OUTPUTS of its control system, so
# they converge in the grid like everything else -- P2 moves 14% between 129^2
# and 513^2 on the limited machine, and its psi_bndry contact moves 0.6 m to
# the other side of the plasma.  examples/limited-tokamak-filament.toml carries
# ref-n513/'s, so racing it against a 257^2 reference is scoring one machine
# against another.  MEASUREMENTS.md M-161's own case file says the same thing.
for _i, _a in enumerate(sys.argv):
    if _a == "--grids" and _i + 1 < len(sys.argv):
        FGS_GRIDS = [int(g) for g in sys.argv[_i + 1].split(",")]
        del sys.argv[_i:_i + 2]
        break

# Where the seeded 513^2 ladder of M-88 landed. Accuracy does not care how a
# reference was reached; only the TIMED rows have to be cold.
TRUTH = os.environ.get("MEQ_RACE_TRUTH", "")

# PHYSICAL cores, not os.cpu_count()'s LOGICAL 16. This machine is an 8-core
# Ryzen with 2 threads per core, so asking for 16 gave half the "threads" to SMT
# siblings sharing an FPU with a thread already saturating it.
#
# M-111 WAS TAKEN AT 16 AND ITS ABSOLUTE SECONDS ARE THEREFORE NOT MEQ'S BEST.
# Both arms got 16, so the RATIOS are probably close to right -- but MEQ's
# threaded assembly plus PARDISO and freegs4e's numpy are not equally SMT
# sensitive, so "probably" is doing real work in that sentence and the race
# wants re-taking before 3.87x is quoted again.
from cores import physical_cores
THREADS = str(physical_cores())


def environment():
    e = dict(os.environ)
    e["OMP_NUM_THREADS"] = THREADS
    e["MKL_NUM_THREADS"] = THREADS
    return e


# ---------------------------------------------------------------------------
# MEQ
# ---------------------------------------------------------------------------
def filament_guess(stem, ref, scratch):
    """A COLD psi_p guess, for the filament arm, built once per case.

    `mkcoldguess.py --remainder` sums the design plasma blob and NOT the
    conductors, so this is a guess at the remainder the split solves for.  It
    is COLD in the sense race.py's header requires -- the design footprint
    carrying the target current, and no converged field anywhere in it.

    CONVERTING THE SHIPPED TOTAL GUESS INSTEAD DOES NOT WORK, and it is worth
    saying because `[initialguess] Content = "total"` makes it look available:
    psi_c diverges logarithmically at each filament while any stored psi
    cannot, so `psi - psi_c` at a node near one is a large negative number
    where the true psi_p is smooth.  Measured on DIII-D, the largest such shift
    is 6.790e-01 Wb/rad against a reference psi_axis of 3.759e-01.
    """
    mesh = os.path.join(scratch, "%s-filament-guess.mesh" % stem)
    gf = os.path.join(scratch, "%s-filament-guess.gf" % stem)
    if os.path.exists(gf) and os.path.exists(mesh):
        return mesh, gf
    # The disc the guess has to cover is [mesh.generate] Radius, out of the
    # case's own file rather than a table here.
    text = open(os.path.join(EXAMPLES, "%s.toml" % stem)).read()
    rho = float(re.search(r"^\[mesh\.generate\][^\[]*?^Radius = ([-\d.eE+]+)",
                          text, re.M | re.S).group(1))
    argv = [sys.executable, os.path.join(HERE, "mkcoldguess.py"), "--remainder",
            os.path.join(HERE, "%s.npz" % ref), mesh, gf, "%g" % rho, "32"]
    done = subprocess.run(argv, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit("the filament guess could not be built for %s:\n%s"
                         % (stem, done.stdout + done.stderr))
    return mesh, gf


def meq_toml(stem, degree, refine, scratch, sample=513, adaptive=0, ref=""):
    """One rung's configuration, derived from the case's own file."""
    # WHICH FILE, AND IT IS NOT ALWAYS THE CASE'S.  See FILAMENT_FILE.
    shipped = FILAMENT_FILE.get(stem) if FILAMENT else None
    source = shipped or stem
    text = open(os.path.join(EXAMPLES, "%s.toml" % source)).read()
    text = re.sub(r"^PolynomialDegree = .*$", "PolynomialDegree = %d" % degree,
                  text, flags=re.M)
    # RefinementLevels is not in the generated file, so it is inserted into
    # [mesh] rather than substituted.
    text = re.sub(r'^(File = "examples/%s\.msh")$' % re.escape(source),
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
    if shipped:
        # ALREADY ALL FOUR, and the file says so in its own header.  Applying
        # the edits below to it would strip a CoilSize that is not there,
        # repoint a mesh whose generator block IS the point, and overwrite a
        # cold guess built from this case's own reference with one built from
        # a different grid's.
        label += "-fil"
    elif FILAMENT:
        # (1) CoilSize is refused beside a subtracting model -- it grades the
        # mesh around conductors the generator is no longer told about.
        text = re.sub(r"^CoilSize = .*\n", "", text, flags=re.M)
        # (2) its own mesh, so the two arms do not share a .meq-mesh stamp.
        text = re.sub(r'^File = "examples/%s\.msh"$' % re.escape(stem),
                      'File = "%s/%s-filament.msh"' % (scratch, stem),
                      text, flags=re.M)
        # (3) the guess is psi_p and is COLD; see filament_guess().
        mesh_path, gf_path = filament_guess(stem, ref, scratch)
        text = re.sub(r'^File = "examples/[^"]*-guess\.gf"$',
                      'File = "%s"' % gf_path, text, flags=re.M)
        text = re.sub(r'^MeshFile = "examples/[^"]*-guess\.mesh"$',
                      'MeshFile = "%s"' % mesh_path, text, flags=re.M)
        # (4) and the key itself.
        text += '\n[conductors]\nModel = "filament"\n'
        label += "-fil"
    text = re.sub(r'^Prefix = ".*"$', 'Prefix = "%s"' % label, text, flags=re.M)
    text = re.sub(r'^(\[output\])$', r'\1\nDirectory = "%s"' % scratch,
                  text, flags=re.M)
    path = os.path.join(scratch, "%s.toml" % label)
    open(path, "w").write(text)
    return path, label


def rung_label(label):
    """"k2r0" out of "<stem>-k2r0" or "<stem>-k2r0-fil"."""
    parts = [p for p in label.split("-") if p and p[0] == "k"]
    return parts[-1] if parts else label


def run_meq(stem, degree, refine, scratch, sample=513, adaptive=0, ref=""):
    path, label = meq_toml(stem, degree, refine, scratch, sample, adaptive, ref)
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
                               "--nx=%d" % nx]
                              + (["--shaped"] if SHAPED else []) + [ref[0]],
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
    # THROUGH with_design(), BECAUSE THE OLDER REFERENCES CARRY NO COIL
    # POSITIONS.  fgsref.py only started saving `coil_R`/`coil_Z` and the
    # `design_*` block when make_diverted_case.py needed them, so
    # H_limited_circular.npz -- and every file under ref-n257/ and ref-n513/ --
    # has the labels and the currents and not the geometry.  with_design()
    # recovers both from fgsref.py's own CASES table, which is what wrote them
    # in the first place, and asserts the label order against the machine
    # rather than assuming it.
    from mkcoldguess import with_design
    reference = os.path.join(HERE, "%s.npz" % ref)
    flat = conductors.from_npz(with_design(np.load(reference,
                                                   allow_pickle=True),
                                           reference))
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
                      adaptive=rung[2] if len(rung) > 2 else 0, ref=ref)
        if not row["ok"] or "nc" not in row or not os.path.exists(row["nc"]):
            # THE RUNG, NOT THE SUFFIX.  `label` is "<stem>-k2r0[-fil]", so
            # taking the last dash-separated field printed "fil" for every
            # failure under --filament and three rungs were indistinguishable
            # in the first run of that flag.
            print("    %-10s FAILED" % rung_label(row["label"]))
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
    # THE BANNER HAS TO NAME THE CONDUCTOR MODEL CORRECTLY, because it is the
    # one line a reader of a pasted table sees and there are now THREE cases.
    # It said "NOT the same conductor model, see --shaped" under --filament
    # once, which is the exact opposite of what --filament does.
    if SHAPED:
        models = ", SHAPED conductors on both sides -- the models MATCH"
    elif FILAMENT:
        models = (", FILAMENT conductors on both sides -- the models MATCH, "
                  "MEQ's subtracted at the reference's own positions")
    else:
        models = (", filament reference against MEQ's rectangles -- NOT the "
                  "same conductor model, see --shaped and --filament")
    print("\n  A RACE, BOTH CODES COLD, %s threads each%s" % (THREADS, models))
    print("  MEQ: PARDISO trace solver and threaded assembly, which are its")
    print("  defaults; its wall includes gmsh, the solve and four output")
    print("  formats. freegs4e's includes its boundary matrix, the Picard")
    print("  loop and its own diagnostics. Neither is a solve time.")
    for tag, ref, stem in CASES:
        if tag in want:
            # BOTH HALVES TOGETHER -- see SHAPED's comment for why a shaped
            # reference against a filament-derived MEQ config would be measuring
            # the conversion rather than either solver.
            sweep(tag, shaped_ref(ref), shaped_stem(stem), scratch,
                  MEQ_RUNGS, FGS_GRIDS)
