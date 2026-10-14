#!/usr/bin/env python
"""THE freegsnke ARM OF THE INVERSE RACE (INVERSE-RACE-PLAN.md, RACE-ARMS.md).

    race_arm_nke.py --problem P.json --config NAME --rung R --start START
                    [--state DIR] --out DIR
    race_arm_nke.py --problem P.json --config NAME --list-rungs

One freegsnke static inverse solve per invocation, writing DIR/record.json in
RACE-ARMS.md's shape and, for a later warm start, DIR/state/.  The machine,
the rows and the profile are `inverse_nke.py`'s (M-255): every supply one
circuit with only the free ones controllable, every shaped conductor a 6 x 6
tensor Gauss-Legendre rule over its rectangle, freegsnke's own null-point and
isoflux rows at weight 1, and GeneralPprimeFFprime with the analytic shape and
one Ip scale.

THE LADDER is the grid, square, 2^n + 1: 65, 129, 257.  257 is the accuracy
datum.  513 is not on it: freegsnke was killed at 16 GB on F there (M-255).

THE STARTS
  cold    fixed supplies at reference_current, free at `vacuum_seed`; the
          plasma is the DESIGN BLOB the vacuum seed was computed against
          (mkcoldguess.plasma_filaments' ellipse and profile-shaped current,
          blob()), represented as freegsnke represents a plasma -- that
          current's free-boundary flux on its grid (blob_plasma_psi()) --
          and then freegsnke's inverse directly.  Against this blob the seed
          meets the rows to 3e-4 at the first step, and the inverse walks to
          the reference branch.  Two alternatives were measured on F at nx65
          and both land on M-255's second solution (F4A ~ 2 MA, axis at
          z = +0.74): freegs4e's PsiGuessParabolicJtor (its own ellipse, psi
          = 0 on the grid edge: rows 0.13-0.3 T off at the start), and a
          forward solve at the seed currents before the inverse.
  seeded  free supplies at reference_current; ONE freegsnke forward solve at
          those currents from the same blob, then the inverse.  The forward
          solve first is what made this start robust in M-255 (every square
          configuration converged).
  warm    --state DIR: this arm's own converged design-target run at the same
          rung -- its plasma_psi and every supply current -- and the SHIFTED
          targets; the inverse only, handing over to freegsnke's full
          (plasma-aware) Jacobian from the first step (START_KW).

THE ROWS DECIDE, NOT THE STEP RATCHET.  freegsnke's inverse loop stops when
the GS residual and the size of its last current step are small, and that
step is capped by a ratchet and damping**k, so it reports SUCCESS with rows at
0.1 T (M-255).  Here the rows are evaluated from the converged FIELD
(eq.Br, eq.Bz, eq.psiRZ at the targets, not freegsnke's constraint vector)
and the run is `finished` only if every row is within ROW_TOL (1e-8 T and
Wb/rad, INVERSE-RACE-PLAN.md section 6) and the GS residual within 10*GS_TOL.
The timed rungs ask freegsnke for GS_TOL; the datum asks for 10*GS_TOL, the
bar itself, because at nx257 the residual can floor between the two and every
forward solve inside the inverse then runs to its cap (DATUM_GS_TOL).
If freegsnke stops short, the inverse is called again from where it stopped
-- which resets the ratchet -- up to MAX_RESTARTS times; the count is in
the record.

THE LIMITER is a harness device: freegsnke requires one and restricts its
O-point search to it, and conductors of the plasma's sign inside the grid are
local maxima of psi.  Here it is a star polygon about the design ellipse's
centre, cut back from the grid edge and from every conductor rectangle, built
from the problem's geometry and design targets only.  It is inert iff
freegsnke never reports the plasma limited; a limited run is not finished.

No patch to freegsnke or freegs4e.  Threads come from the environment.
"""
import time

T0 = time.perf_counter()

import argparse                                                 # noqa: E402
import contextlib                                               # noqa: E402
import io                                                       # noqa: E402
import json                                                     # noqa: E402
import os                                                       # noqa: E402
import subprocess                                               # noqa: E402
import sys                                                      # noqa: E402

import numpy as np                                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import benchpaths                                               # noqa: E402

LADDER = [65, 129, 257]
ROW_TOL = 1e-8
GS_TOL = 1e-10
# THE DATUM ASKS FOR THE ACCEPTANCE BAR, NOT TEN TIMES BELOW IT.  freegsnke
# hands its inverse's tolerance to every forward solve inside it -- one before
# each full Jacobian, one per controlled coil for its columns, one after the
# step -- and on D shaped at nx257 the relative GS residual floors at 4e-10:
# the rows were met to 1e-16 by iteration 29, and the inverse and every inner
# solve then ran to their caps, 100 iterations at 700 s each, 70,461 s.  A run
# is accepted at 10*GS_TOL whatever it asks for, so the untimed datum asks for
# exactly that; the timed rungs keep GS_TOL, so their walls stay comparable
# with every record already taken.
DATUM_GS_TOL = 10*GS_TOL
PSIT_TOL = 1e-9
MAX_ITS = 100
MAX_RESTARTS = 4
NQUAD = 6
# freegsnke inverse-solve keywords beyond the tolerances; tuned in the
# docstring's "THE CONTROLS" paragraph
SOLVE_KW = {}
# per start.  WARM begins at a converged state and the rows are a small
# perturbation away, which is where a full-Jacobian (Gauss-Newton) step is the
# right one; freegsnke's default handover ([1e-5, 7e-3]: GS residual, flux
# step) keeps it on the frozen-plasma Green's Jacobian while each forward
# re-solve leaves the GS residual at 1e-5 to 3e-5, and on F at nx65 that
# iteration walks AWAY from the shifted targets (rows 7e-3 -> 0.5 T).
START_KW = {"warm": {"full_jacobian_handover": [1e-3, 5e-2]}}
# run one forward solve at the start currents before the inverse (cold/seeded)
FORWARD_FIRST = {"cold": False, "seeded": True}


def list_rungs():
    return [dict(rung=k, label="nx%d" % n, finest=(k == len(LADDER) - 1))
            for k, n in enumerate(LADDER)]


def sha(path):
    try:
        out = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", path, "status", "--porcelain",
                                "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return out + ("-dirty" if dirty else "")
    except Exception:                                       # noqa: BLE001
        return "?"


# ---------------------------------------------------------------------------
def design_targets(conf, start):
    if start == "warm":
        return (np.asarray(conf["shifted_xpoint_targets"], float).reshape(-1, 2),
                [list(p) for p in conf["shifted_isoflux_targets"]])
    return (np.asarray(conf["xpoint_targets"], float).reshape(-1, 2),
            [list(p) for p in conf["isoflux_targets"]])


def design_footprint(conf):
    """mkcoldguess.design_footprint(pad=0) from the problem's DESIGN block:
    the box the design's own boundary points name, an unnamed radial edge
    half the X-point half-height beyond the outermost point (these shaped
    references carry no wall)."""
    d = conf["design"]
    pts = [tuple(p) for p in np.asarray(d["xpoints"], float).reshape(-1, 2)]
    for R1, Z1, R2, Z2 in d["isoflux"]:
        pts += [(R1, Z1), (R2, Z2)]
    pts = np.asarray(pts, float)
    r0 = float(d["R0"])
    zmin, zmax = pts[:, 1].min(), pts[:, 1].max()
    zc = 0.5*(zmin + zmax)
    zhalf = max(0.5*(zmax - zmin), 1e-6)
    mid = [(R, z) for R, z in pts if abs(z - zc) < 0.15*zhalf]
    named_in = [R for R, _ in mid if R < r0]
    named_out = [R for R, _ in mid if R > r0]
    Rmin = min(named_in) if named_in else max(0.02, pts[:, 0].min() - 0.5*zhalf)
    Rmax = max(named_out) if named_out else pts[:, 0].max() + 0.5*zhalf
    return Rmin, Rmax, zmin, zmax, [tuple(p) for p in pts]


def blob(conf):
    """THE DESIGN BLOB, as the problem's `vacuum_seed` was computed against
    it (mkcoldguess.plasma_filaments): Ip spread over the ellipse inscribed
    in the design footprint with the design's own profile shapes and split,

        j(R, x) = R (1 - x^pa)^pb + f0 (1 - x^fa)^fb / (mu0 R),
        f0 = mu0 R0^2 (1 - frac_p)/frac_p,   x = rho^2 elliptical,

    clipped at zero, normalised to the design Ip.  freegsnke represents it
    as its plasma flux: that current's free-boundary field on its own grid
    (blob_plasma_psi).  Returns (jtor(R, Z) function, ellipse, design
    points)."""
    d = conf["design"]
    Rmin, Rmax, zmin, zmax, pts = design_footprint(conf)
    rc, zc = 0.5*(Rmin + Rmax), 0.5*(zmin + zmax)
    ra, za = 0.5*(Rmax - Rmin), 0.5*(zmax - zmin)
    mu0 = 4e-7*np.pi
    frac_p = float(d["frac_p"])
    f0 = mu0*float(d["R0"])**2*(1.0 - frac_p)/frac_p
    pa, pb, fa, fb = (float(d[k]) for k in ("pa", "pb", "fa", "fb"))

    def jtor(R, Z):
        rho2 = ((R - rc)/ra)**2 + ((Z - zc)/za)**2
        x = np.clip(rho2, 0.0, 1.0)
        w = (R*np.power(1.0 - np.power(x, pa), pb)
             + f0*np.power(1.0 - np.power(x, fa), fb)/(mu0*R))
        return np.where((rho2 < 1.0) & (w > 0.0), w, 0.0)

    return jtor, (rc, zc, ra, za), pts


def star_limiter(conf, centre, keep, inset_cells=3, nrays=1440):
    """Star polygon about `centre`: the grid edge less a few cells, cut back
    to every conductor rectangle.  `keep` are points that must be inside;
    returns (R, Z, worst clearance of `keep`, negative if one is outside)."""
    from matplotlib.path import Path
    g = conf["grid"]
    R0, Z0 = centre
    cell = np.array([(g["Rmax"] - g["Rmin"])/128.0, (g["Zmax"] - g["Zmin"])/128.0])
    lo = np.array([g["Rmin"], g["Zmin"]]) + inset_cells*cell
    hi = np.array([g["Rmax"], g["Zmax"]]) - inset_cells*cell
    rects = [(c["R"] - c["half_width"], c["R"] + c["half_width"],
              c["Z"] - c["half_height"], c["Z"] + c["half_height"])
             for c in conf["conductors"]]
    pts = []
    for th in np.linspace(0.0, 2*np.pi, nrays, endpoint=False):
        d = np.array([np.cos(th), np.sin(th)])
        t = np.inf
        for k in range(2):
            if d[k] > 1e-14:
                t = min(t, (hi[k] - (R0, Z0)[k])/d[k])
            elif d[k] < -1e-14:
                t = min(t, (lo[k] - (R0, Z0)[k])/d[k])
        for r0, r1, z0, z1 in rects:
            tmin, tmax, ok = 0.0, np.inf, True
            for o, dd, a, b in ((R0, d[0], r0, r1), (Z0, d[1], z0, z1)):
                if abs(dd) < 1e-14:
                    if o < a or o > b:
                        ok = False
                        break
                else:
                    t1, t2 = (a - o)/dd, (b - o)/dd
                    tmin, tmax = max(tmin, min(t1, t2)), min(tmax, max(t1, t2))
            if ok and tmin <= tmax:
                t = min(t, tmin)
        pts.append((R0 + t*d[0], Z0 + t*d[1]))
    pts = np.array(pts)
    keep = np.asarray(keep, float).reshape(-1, 2)
    inside = Path(pts).contains_points(keep)
    a, b = pts, np.roll(pts, -1, axis=0)
    ab = b - a
    worst = np.inf
    for p, ins in zip(keep, inside):
        tt = np.clip(np.sum((p - a)*ab, 1)/np.maximum(np.sum(ab*ab, 1), 1e-30), 0, 1)
        dist = float(np.min(np.hypot(*(a + tt[:, None]*ab - p).T)))
        worst = min(worst, dist if ins else -dist)
    return pts[:, 0], pts[:, 1], worst


def mirror_symmetric(conf):
    """True if the grid, the conductors and the held currents are up-down
    symmetric, so forcing a symmetric plasma is consistent."""
    g = conf["grid"]
    if abs(g["Zmin"] + g["Zmax"]) > 1e-12:
        return False
    cur = {}
    for s in conf["supplies"]:
        for m in s["members"]:
            cur[m["conductor"]] = (s["name"], m["multiplier"], s["free"],
                                   s["reference_current"])
    conds = conf["conductors"]
    for c in conds:
        twin = [d for d in conds if abs(d["R"] - c["R"]) < 1e-6
                and abs(d["Z"] + c["Z"]) < 1e-6
                and abs(d["half_width"] - c["half_width"]) < 1e-6
                and abs(d["half_height"] - c["half_height"]) < 1e-6]
        if not twin:
            return False
        a, b = cur[c["label"]], cur[twin[0]["label"]]
        if a[0] != b[0] and not (a[2] or b[2]):
            if abs(a[1]*a[3] - b[1]*b[3]) > 1e-9*max(1.0, abs(a[1]*a[3])):
                return False
        if a[0] == b[0] and a[1] != b[1] and not a[2]:
            return False             # a HELD anti-series pair
        # a FREE anti-series pair (MAST's P6, MAST-U's P61/P62) is allowed:
        # with the plasma forced symmetric the only consistent current in it
        # is zero, which is what the references carry (1e-3 A and below)
    return True


def blob_plasma_psi(solver, eq, jfun, Ip):
    """The blob's FREE-BOUNDARY plasma flux on freegsnke's grid: its current,
    its edge flux from freegsnke's own Green's boundary matrix, one linear GS
    solve -- what NKGSsolver.freeboundary does each step.  (freegs4e's
    PsiGuessParabolicJtor puts psi = 0 on the grid edge instead, which is not
    a plasma flux in freegsnke's sense, and its ellipse is not the design's:
    measured on F at nx65 that left the vacuum seed's rows 0.13-0.3 T off.)"""
    jtor = jfun(eq.R, eq.Z)
    jtor = jtor*solver.plasma_source_mask
    jtor *= Ip/(np.sum(jtor)*eq.dR*eq.dZ)
    bnd = np.zeros_like(eq.R)
    flux = solver._boundary_flux_from_jtor(jtor)
    nx, ny = eq.R.shape
    bnd[:, 0] = flux[:nx]
    bnd[:, -1] = flux[nx:2*nx]
    bnd[0, 1:ny - 1] = flux[2*nx:2*nx + ny - 2]
    bnd[-1, 1:ny - 1] = flux[2*nx + ny - 2:]
    rhs = solver.rhs_before_jtor*jtor
    rhs[0, :], rhs[:, 0] = bnd[0, :], bnd[:, 0]
    rhs[-1, :], rhs[:, -1] = bnd[-1, :], bnd[:, -1]
    return np.asarray(solver.linear_GS_solver(bnd, rhs)).reshape(nx, ny)


def gs_residual(solver, eq, profiles):
    """freegsnke's own relative GS residual of the current state."""
    gs = solver.F_function(tokamak_psi=eq.tokamak.getPsitokamak(
        vgreen=eq._vgreen).reshape(-1), plasma_psi=eq.plasma_psi.reshape(-1),
        profiles=profiles)
    return float(solver.relative_del_residual(gs, eq.plasma_psi)[0])


def field_rows(eq, xt, iso):
    rows = []
    for R, Z in xt:
        rows.append(float(eq.Br(R, Z)))
        rows.append(float(eq.Bz(R, Z)))
    for R1, Z1, R2, Z2 in iso:
        rows.append(float(eq.psiRZ(R1, Z1) - eq.psiRZ(R2, Z2)))
    return rows


# ---------------------------------------------------------------------------
def run(args):
    from inverse_nke import active_coils, make_profile_class, locate_null
    from freegsnke import build_machine, equilibrium_update, GSstaticsolver
    from freegsnke.inverse import Inverse_optimizer

    P = {c["config"]: c for c in json.load(open(args.problem))}
    conf = P[args.config]
    nx = LADDER[args.rung]
    gs_tol = DATUM_GS_TOL if args.rung == len(LADDER) - 1 else GS_TOL
    start = args.start
    xt, iso = design_targets(conf, start)
    notes = []
    rec = dict(code="freegsnke", config=args.config, rung=args.rung,
               rung_label="nx%d" % nx, unknowns=nx*nx, start=start,
               finished=False, status="", wall_s=None, solve_s=None,
               currents={}, psi_axis=None, psi_bndry=None, axis=None,
               xpoints=[], rows=[], iterations=0, state=None)

    free = [s["name"] for s in conf["supplies"] if s["free"]]
    state = None
    if start == "warm":
        if not args.state:
            raise SystemExit("--start warm needs --state DIR")
        sdir = args.state
        if os.path.isfile(os.path.join(sdir, "record.json")):
            prior = json.load(open(os.path.join(sdir, "record.json")))
            sdir = prior.get("state") or os.path.join(sdir, "state")
        state = np.load(os.path.join(sdir, "state.npz"), allow_pickle=True)
        meta = json.load(open(os.path.join(sdir, "state.json")))
        if meta["config"] != args.config or int(meta["nx"]) != nx:
            raise SystemExit("--state is %s at nx %s, not %s at nx %d" % (
                meta["config"], meta["nx"], args.config, nx))
        start_currents = meta["supply_currents"]
    else:
        start_currents = {}
        for s in conf["supplies"]:
            if s["free"] and start == "cold":
                start_currents[s["name"]] = float(conf["vacuum_seed"][s["name"]])
            else:
                start_currents[s["name"]] = float(s["reference_current"])

    # the machine
    jfun, ell, sep = blob(conf)
    centre = (ell[0], ell[1])
    keep = list(sep) + [tuple(p) for p in np.asarray(
        conf["shifted_xpoint_targets"], float).reshape(-1, 2)]
    for R1, Z1, R2, Z2 in conf["shifted_isoflux_targets"]:
        keep += [(R1, Z1), (R2, Z2)]
    LR, LZ, clearance = star_limiter(conf, centre, keep)
    rec["limiter_clearance_of_targets_m"] = clearance
    quiet = io.StringIO()
    with contextlib.redirect_stdout(quiet):
        tok = build_machine.tokamak(
            active_coils_data=active_coils(conf, NQUAD),
            limiter_data=[dict(R=float(r), Z=float(z)) for r, z in zip(LR, LZ)])
    for s in conf["supplies"]:
        tok.set_coil_current(s["name"], float(start_currents[s["name"]]))
        tok[s["name"]].control = bool(s["free"])
    tok.getCurrentsVec()

    g = conf["grid"]
    with contextlib.redirect_stdout(quiet):
        eq = equilibrium_update.Equilibrium(
            tokamak=tok, Rmin=g["Rmin"], Rmax=g["Rmax"], Zmin=g["Zmin"],
            Zmax=g["Zmax"], nx=nx, ny=nx,
            psi=(np.asarray(state["plasma_psi"]) if state is not None else None))
    Prof = make_profile_class(conf)
    profiles = Prof(eq=eq, Ip=float(conf["Ip"]), fvac=float(conf["fvac"]),
                    psi_n=np.linspace(0.0, 1.0, 256),
                    pprime_data=np.ones(256), ffprime_data=np.ones(256),
                    Raxis=1.0, Ip_logic=True)
    solver = GSstaticsolver.NKGSsolver(eq)
    if state is None:
        eq._updatePlasmaPsi(blob_plasma_psi(solver, eq, jfun,
                                            float(conf["design"]["Ip"])))
    symmetric = bool(conf.get("symmetric")) and int(conf.get("rank", 0)) < int(
        conf.get("rows", 0)) and mirror_symmetric(conf)
    rec["force_up_down_symmetric"] = symmetric

    constrain = Inverse_optimizer(
        null_points=[list(xt[:, 0]), list(xt[:, 1])] if len(xt) else None,
        isoflux_set=[[[R1, R2], [Z1, Z2]] for R1, Z1, R2, Z2 in iso] or None)
    constrain.prepare_for_solve(eq)
    constrain.build_plasma_vals(np.zeros_like(eq.plasma_psi))
    constrain.build_lsq(eq.tokamak.getCurrentsVec())
    sv = np.linalg.svd(np.array(constrain.A), compute_uv=False)
    reg = float(1e-8*sv[0]**2)

    log = io.StringIO()
    solve_s = 0.0
    try:
        if FORWARD_FIRST.get(start, False):
            t = time.perf_counter()
            with contextlib.redirect_stdout(log):
                solver.solve(eq=eq, profiles=profiles, constrain=None,
                             target_relative_tolerance=gs_tol,
                             max_solving_iterations=200,
                             force_up_down_symmetric=symmetric, verbose=True)
            solve_s += time.perf_counter() - t
            rec["forward_relative_residual"] = float(solver.relative_change)
        restarts = 0
        attempts = []
        best = None
        while True:
            t = time.perf_counter()
            with contextlib.redirect_stdout(log):
                solver.solve(eq=eq, profiles=profiles, constrain=constrain,
                             target_relative_tolerance=gs_tol,
                             target_relative_psit_update=PSIT_TOL,
                             max_solving_iterations=MAX_ITS,
                             l2_reg=reg, l2_reg_fj=reg,
                             force_up_down_symmetric=symmetric, verbose=True,
                             **{**START_KW.get(start, {}), **SOLVE_KW})
            solve_s += time.perf_counter() - t
            rows = field_rows(eq, xt, iso)
            worst = max(abs(r) for r in rows)
            gs_rel = gs_residual(solver, eq, profiles)
            attempts.append(dict(max_row=worst, gs_relative_residual=gs_rel))
            # a re-entry may make things WORSE (a Green's-function step on a
            # plasma it then fails to re-converge), so keep the best state
            # whose plasma is converged and restore it at the end
            if gs_rel <= 10*GS_TOL and (best is None or worst < best[0]):
                best = (worst, np.copy(eq.plasma_psi),
                        np.copy(eq.tokamak.getCurrentsVec()))
            if (worst <= ROW_TOL and gs_rel <= 10*GS_TOL) or restarts >= MAX_RESTARTS:
                break
            if gs_rel > 1e-3:
                break                # the re-entry has left the basin
            restarts += 1
        if best is not None and not (attempts[-1]["gs_relative_residual"] <= 10*GS_TOL
                                     and attempts[-1]["max_row"] <= best[0]):
            eq.tokamak.set_all_coil_currents(best[2])
            eq.tokamak.getCurrentsVec()
            eq._updatePlasmaPsi(best[1])
            rec["restored_best_attempt"] = True
        rec["attempts"] = attempts
        rec["restarts"] = restarts
        status = None
    except Exception as e:                                  # noqa: BLE001
        import traceback
        status = "exception: %s: %s" % (type(e).__name__, str(e)[:300])
        rec["traceback"] = traceback.format_exc()[-1500:]
    rec["solve_s"] = solve_s
    rec["gs_target"] = gs_tol
    text = log.getvalue()
    open(os.path.join(args.out, "freegsnke.log"), "w").write(text)
    rec["iterations"] = text.count("\nIteration: ") + text.startswith("Iteration: ")
    rec["forward_iterations"] = text.count("Picard iteration") + text.count(
        "Newton-Krylov iteration")
    if status is not None:
        rec["status"] = status
        return rec

    # the state reported: re-evaluate the profile on the final psi
    psi = eq.psi()
    profiles.Jtor(eq.R, eq.Z, psi)
    gs_rel = gs_residual(solver, eq, profiles)
    profiles.Jtor(eq.R, eq.Z, psi)
    opt = np.atleast_2d(profiles.opt)
    cur = eq.tokamak.getCurrents()
    rows = field_rows(eq, xt, iso)
    rec.update(
        currents={f: float(cur[f]) for f in free},
        psi_axis=float(opt[0][2]), psi_bndry=float(profiles.psi_bndry),
        axis=[float(opt[0][0]), float(opt[0][1])],
        xpoints=[x[:2] if x else None for x in
                 (locate_null(eq, R, Z) for R, Z in xt)],
        rows=rows, gs_relative_residual=gs_rel,
        flag_limiter=bool(profiles.flag_limiter),
        Ip_scale_L=float(profiles.L),
        Ip=float(np.sum(profiles.jtor)*eq.dR*eq.dZ))
    worst = max(abs(r) for r in rows)
    why = []
    if worst > ROW_TOL:
        why.append("rows not met: max |row| %.2e > %.0e" % (worst, ROW_TOL))
    if gs_rel > 10*GS_TOL:
        why.append("GS residual %.2e" % gs_rel)
    if rec["flag_limiter"]:
        why.append("plasma limited by the harness limiter")
    if any(x is None for x in rec["xpoints"]):
        why.append("no field null within 0.2 m of a target")
    rec["finished"] = not why
    rec["status"] = "converged" if not why else "; ".join(why)

    # the state a warm start restarts from
    sdir = os.path.join(args.out, "state")
    os.makedirs(sdir, exist_ok=True)
    np.savez_compressed(os.path.join(sdir, "state.npz"),
                        plasma_psi=eq.plasma_psi, R=eq.R[:, 0], Z=eq.Z[0, :])
    json.dump(dict(config=args.config, nx=nx, start=start,
                   finished=rec["finished"],
                   supply_currents={k: float(v) for k, v in cur.items()}),
              open(os.path.join(sdir, "state.json"), "w"), indent=1)
    rec["state"] = sdir
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--rung", type=int)
    ap.add_argument("--start", choices=("cold", "seeded", "warm"))
    ap.add_argument("--state")
    ap.add_argument("--out")
    ap.add_argument("--list-rungs", action="store_true")
    args = ap.parse_args()
    if args.list_rungs:
        print(json.dumps(list_rungs()))
        return 0
    if args.rung is None or args.start is None or not args.out:
        ap.error("--rung, --start and --out are required")
    if not 0 <= args.rung < len(LADDER):
        ap.error("--rung is 0..%d" % (len(LADDER) - 1))
    os.makedirs(args.out, exist_ok=True)
    try:
        rec = run(args)
    except SystemExit:
        raise
    except Exception as e:                                  # noqa: BLE001
        rec = dict(code="freegsnke", config=args.config, rung=args.rung,
                   rung_label="nx%d" % LADDER[args.rung],
                   unknowns=LADDER[args.rung]**2, start=args.start,
                   finished=False, status="exception: %s: %s" % (
                       type(e).__name__, str(e)[:300]),
                   solve_s=None, currents={}, rows=[], iterations=0, state=None)
    import freegsnke
    import freegs4e
    rec["provenance"] = dict(
        code_version="freegsnke %s on freegs4e %s" % (
            getattr(freegsnke, "__version__", "?"),
            getattr(freegs4e, "__version__", "?")),
        tree=dict(freegsnke=sha(benchpaths.FREEGSNKE),
                  freegs4e=sha(benchpaths.FREEGS4E),
                  meq=sha(benchpaths.SOURCE)),
        notes=("no patches to freegsnke or freegs4e; freegsnke 3.1.0 pins "
               "freegs4e>=0.14 and runs on 0.12. Harness: analytic-profile "
               "subclass of GeneralPprimeFFprime; shaped conductors a %dx%d "
               "Gauss-Legendre rule as single-filament circuit members; a "
               "star limiter cut at the grid edge and conductors; l2_reg "
               "1e-8 sigma_max^2 of the vacuum response; rows from eq.Br/Bz/"
               "psiRZ at the targets; cold starts from the design blob's "
               "free-boundary flux; forward solve before the inverse: %s; "
               "per-start solve keywords %s; GS tolerance %.0e asked "
               "(%.0e on the datum rung), %.0e accepted; inverse "
               "re-entered up to %d times while rows > %.0e." % (
                   NQUAD, NQUAD, FORWARD_FIRST, START_KW, GS_TOL,
                   DATUM_GS_TOL, 10*GS_TOL, MAX_RESTARTS, ROW_TOL)))
    rec["wall_s"] = time.perf_counter() - T0
    json.dump(rec, open(os.path.join(args.out, "record.json"), "w"), indent=1)
    print(json.dumps(dict(config=rec["config"], rung=rec["rung"],
                          start=rec["start"], finished=rec["finished"],
                          status=rec["status"], wall_s=round(rec["wall_s"], 2))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
