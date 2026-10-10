#!/usr/bin/env python
"""THE SAME INVERSE PROBLEM IN freegsnke: a third code for M-254's consensus.

`inverse_compare.py` poses fourteen inverse problems to MEQ and compares the
currents against freegs4e's.  This poses the IDENTICAL problem -- read from the
brief's problem.json, never re-derived -- to freegsnke's static inverse solve
(`GSstaticsolver.NKGSsolver.inverse_solve` driven by `inverse.Inverse_optimizer`)
and writes one record per (configuration, resolution, arm).

WHAT IS MAPPED, AND HOW
  * every SUPPLY is one freegsnke circuit, and only the `free` ones have
    `control = True`; every other supply is held at its reference current.
  * every CONDUCTOR is a member of its supply's circuit with the problem's
    multiplier, modelled as its `kind` says:
      filament  one filament at the centre
      shaped    an n x n tensor Gauss-Legendre rule over the rectangle, one
                single-filament member per point carrying multiplier x weight,
                the weights summing to one: a uniform current density to the
                rule's accuracy (n = 6 by default)
      solenoid  freegs4e's own Solenoid geometry -- Ns filaments equally
                spaced from Z - hh to Z + hh at R, each carrying I/Ns.  Ns is
                NOT in problem.json and is taken from freegs4e.machine
                (MAST P1 100, TCV OH1 100, MAST-U Solenoid 324, Pc 142)
  * the rows are freegsnke's own: null_points [[R...],[Z...]] gives B_R and
    B_Z at each design X-point (T) and each isoflux pair is its OWN isoflux
    set of two points, one row psi(1) - psi(2) (Wb/rad).  All weights 1.
  * the profile is freegsnke's GeneralPprimeFFprime with Ip_logic = True and
    Raxis = 1 -- Jtor = L ( R p'(psi_n) + ff'(psi_n)/(mu0 R) ), one scale L
    fixed by Ip -- with the SPLINE REPLACED BY THE ANALYTIC SHAPE: the stock
    class fits a smoothing UnivariateSpline to a table, which is exactly the
    1.3e-02 profile difference fgsref.py --interpolating-profiles removed from
    the reference.  L should therefore come out at 1 to the reference's own
    Ip error, and is reported.

TWO DEVICES THAT ARE NOT THE PROBLEM, AND ARE CHECKED TO BE INERT
  * freegsnke REQUIRES a limiter; the reference has none.  It also restricts
    its O-point search to the limiter, and several machines put a conductor
    carrying current of the plasma's sign inside the grid (MAST-U's D1), which
    is a local maximum of psi.  So the limiter is a star polygon about the
    reference axis, cut back from the grid edge and from every conductor
    rectangle grown by a gap.  It is inert iff freegsnke never reports the
    plasma limited (`flag_limiter`), which every record carries, and the
    polygon's clearance from the reference LCFS is recorded.
  * freegsnke's regularisation is an ABSOLUTE Tikhonov weight on the current
    STEP (A^2 units); its defaults (1e-9, 1e-8) are sized for multi-turn
    machines and are four orders above A^T A for these single-turn amp
    currents, which turns the update into a crawl.  It penalises the step, not
    the current, so it does not move the fixed point; it is set here RELATIVE
    to the vacuum response, `--reg` x sigma_max^2.

USAGE
    venv/bin/python inverse_nke.py OUTDIR [--configs a,b] [--nx 129,257]
                                   [--arms seeded,cold] [--jobs 4] [--tol 1e-9]
                                   [--psit 1e-9] [--reg 1e-8] [--quad 6]
    venv/bin/python inverse_nke.py OUTDIR --table
    venv/bin/python inverse_nke.py OUTDIR --one CONFIG --nx N --arm ARM
"""
import argparse
import contextlib
import io
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROBLEM = os.environ.get(
    "INVCMP_PROBLEM",
    "/tmp/claude-1000/-home-ian-projects-meq/91553931-83b4-4f7b-8d5d-9a1e25e12e16"
    "/scratchpad/invcmp/problem.json")
MEQ_RUNS = os.path.join(os.path.dirname(PROBLEM), "interp")

# freegs4e.machine's solenoid turn counts; problem.json gives the extent only.
SOLENOID_NS = {"machine-c-mast": {"P1": 100},
               "machine-d-tcv": {"OH.OH1": 100},
               "machine-g-mastu": {"Solenoid": 324, "Pc": 142}}

# set by --fgs-rule: shaped conductors on freegs4e's 12-point ShapedCoil rule
FGS_RULE = False

# Configurations whose problem.json entry is known to be wrong.  The
# coordinator found machine-g-mastu-shaped's conductor table reproduces about
# 22x freegs4e's own stored coil_psi, so it is not the machine freegs4e
# solved; it is reported, not run.
INVALID = {}

THREAD_ENV = {k: "1" for k in ("MKL_NUM_THREADS", "OMP_NUM_THREADS",
                               "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS")}


def load_problem():
    return {c["config"]: c for c in json.load(open(PROBLEM))}


# ---------------------------------------------------------------------------
# the machine
# ---------------------------------------------------------------------------
def conductor_filaments(conf, cond, nquad, as_rect=False):
    """[(R, Z, weight)] with weights summing to one, by the conductor's kind.

    `as_rect` models EVERY conductor as a uniform rectangle of its own extent,
    which is MEQ's conductor model (it refuses filaments in an inverse solve):
    the arm that compares freegsnke with MEQ like for like on the plain
    configurations, and sizes the filament-versus-rectangle difference inside
    one code."""
    R, Z = cond["R"], cond["Z"]
    hw, hh = cond["half_width"], cond["half_height"]
    kind = "shaped" if as_rect else cond["kind"]
    if kind == "filament":
        return [(R, Z, 1.0)]
    if kind == "solenoid":
        ns = SOLENOID_NS[conf["config"]][cond["label"]]
        return [(R, z, 1.0/ns) for z in np.linspace(Z - hh, Z + hh, ns)]
    if kind == "shaped" and FGS_RULE:
        # freegs4e's OWN ShapedCoil rule: the rectangle as two triangles with
        # a 6-point rule on each.  Only for sizing how much of a difference
        # is the reference's quadrature, never the default.
        from freegs4e import quadrature
        sq = [(R - hw, Z - hh), (R + hw, Z - hh), (R + hw, Z + hh), (R - hw, Z + hh)]
        return [(float(r), float(z), float(w))
                for r, z, w in quadrature.polygon_quad(sq, n=6)]
    if kind == "shaped":
        x, w = np.polynomial.legendre.leggauss(nquad)
        # a tall winding (a solenoid as a rectangle) gets proportionally more
        # points along its height
        nz = max(nquad, int(np.ceil(nquad*hh/max(hw, 0.05))))
        y, v = np.polynomial.legendre.leggauss(nz)
        out = []
        for xi, wi in zip(x, w):
            for yj, wj in zip(y, v):
                out.append((R + hw*xi, Z + hh*yj, 0.25*wi*wj))
        return out
    raise ValueError("unknown conductor kind %r" % kind)


def active_coils(conf, nquad, as_rect=False):
    """freegsnke's active-coil dictionary: one circuit per supply."""
    cond = {c["label"]: c for c in conf["conductors"]}
    out = {}
    for s in conf["supplies"]:
        members = {}
        for m in s["members"]:
            c = cond[m["conductor"]]
            fil = conductor_filaments(conf, c, nquad, as_rect)
            if c["kind"] == "solenoid" and not as_rect:
                # one MultiCoil of Ns filaments, each carrying I x mult / Ns
                members["#" + c["label"]] = dict(
                    R=[f[0] for f in fil], Z=[f[1] for f in fil],
                    dR=2*c["half_width"], dZ=2*c["half_height"]/len(fil),
                    resistivity=1e-7, polarity=1.0,
                    multiplier=float(m["multiplier"])*fil[0][2])
            else:
                for k, (r, z, w) in enumerate(fil):
                    members["#%s#%d" % (c["label"], k)] = dict(
                        R=[r], Z=[z], dR=2*c["half_width"]/max(1, nquad),
                        dZ=2*c["half_height"]/max(1, nquad),
                        resistivity=1e-7, polarity=1.0,
                        multiplier=float(m["multiplier"])*w)
        out[s["name"]] = members
    return out


def star_limiter(conf, gap=0.03, inset_cells=3, nrays=1440):
    """A star polygon about the reference axis, cut back from the grid edge and
    from every conductor rectangle grown by `gap`.  Returns (R, Z, clearance of
    the reference LCFS and design X-points from the polygon)."""
    g = conf["grid"]
    ref = np.load(conf["reference_npz"], allow_pickle=True)
    R0, Z0 = float(ref["Raxis"]), float(ref["Zaxis"])
    h = max((g["Rmax"] - g["Rmin"])/128.0, (g["Zmax"] - g["Zmin"])/128.0)
    cell = np.array([(g["Rmax"] - g["Rmin"])/128.0, (g["Zmax"] - g["Zmin"])/128.0])
    lo = np.array([g["Rmin"], g["Zmin"]]) + inset_cells*cell
    hi = np.array([g["Rmax"], g["Zmax"]]) - inset_cells*cell
    lth = np.arctan2(ref["lcfs_Z"] - Z0, ref["lcfs_R"] - R0)
    lr = np.hypot(ref["lcfs_R"] - R0, ref["lcfs_Z"] - Z0)
    order = np.argsort(lth)
    lth, lr = lth[order], lr[order]
    margin = 2.0*float(cell.min())

    def lcfs_radius(th, window=np.radians(3.0)):
        # the LARGEST LCFS radius within a few degrees, not an interpolant:
        # an X-point corner is not single valued in angle
        dth = np.abs((lth - th + np.pi) % (2*np.pi) - np.pi)
        near = lr[dth < window]
        return float(near.max()) if near.size else float(np.interp(
            (th + np.pi) % (2*np.pi) - np.pi, lth, lr, period=2*np.pi))

    rects = [(c["R"] - c["half_width"] - gap, c["R"] + c["half_width"] + gap,
              c["Z"] - c["half_height"] - gap, c["Z"] + c["half_height"] + gap)
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
        tbox = t
        for r0, r1, z0, z1 in rects:
            # slab test of the ray against the grown rectangle
            tmin, tmax = 0.0, np.inf
            ok = True
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
        # NEVER INSIDE THE REFERENCE PLASMA: MAST-U's inboard LCFS sits 0.024
        # m from its solenoid's invented rectangle, so a conductor cut alone
        # crosses the plasma there and freegsnke reports it LIMITED.  Keep the
        # polygon `margin` outside the reference LCFS along every ray, which
        # may enter a grown conductor rectangle but never the grid edge.
        t = max(t, lcfs_radius(th) + margin)
        t = min(t, tbox)
        pts.append((R0 + t*d[0], Z0 + t*d[1]))
    pts = np.array(pts)
    from matplotlib.path import Path
    poly = Path(pts)
    probe = np.r_[np.c_[ref["lcfs_R"], ref["lcfs_Z"]],
                  np.asarray(conf["xpoint_targets"], float).reshape(-1, 2)]
    inside = poly.contains_points(probe)
    # distance from each probe point to the polygon's edges
    a, b = pts, np.roll(pts, -1, axis=0)
    ab = b - a
    best = np.inf
    for p in probe:
        tt = np.clip(np.sum((p - a)*ab, 1)/np.maximum(np.sum(ab*ab, 1), 1e-30),
                     0, 1)
        best = min(best, float(np.min(np.hypot(*(a + tt[:, None]*ab - p).T))))
    clearance = best if inside.all() else -best
    return pts[:, 0], pts[:, 1], clearance


def shape_function(expr, amp_name):
    """'P0*(1-psi_n**1.5)**2' with the amplitude set to one, vectorised."""
    code = compile(expr, "<profile>", "eval")

    def f(pn):
        pn = np.clip(np.asarray(pn, float), 0.0, 1.0)
        return eval(code, {"__builtins__": {}},
                    {amp_name: 1.0, "psi_n": pn}) + 0.0*pn
    return f


def make_profile_class(conf):
    """GeneralPprimeFFprime with its splines replaced by the analytic shapes.

    A CLASS, not an instance patch, because freegsnke copies profiles
    (`eq._profiles = profiles.copy()`) and the copy calls initialize_profile()
    again, which would rebuild the smoothing spline on an instance patch."""
    from freegsnke.jtor_update import GeneralPprimeFFprime
    fp = shape_function(conf["profile_shape"]["pprime"], "P0")
    ff = shape_function(conf["profile_shape"]["ffprime"], "F0")
    P, F = float(conf["pprime_amplitude"]), float(conf["ffprime_amplitude"])

    class AnalyticProfile(GeneralPprimeFFprime):
        def initialize_profile(self):
            self.pprime_func = lambda pn: P*fp(pn)
            self.ffprime_func = lambda pn: F*ff(pn)
            self.p_func = None
            self.f_func = None

    return AnalyticProfile


def vacuum_response(eq, constrain):
    """The rows' response to the free supplies, plasma frozen: freegsnke's own
    Green's Jacobian (`constrain.A` after one build)."""
    constrain.prepare_for_solve(eq)
    constrain.build_plasma_vals(np.zeros_like(eq.plasma_psi))
    constrain.build_lsq(eq.tokamak.getCurrentsVec())
    return np.array(constrain.A), -np.array(constrain.b)


def row_labels(conf):
    labels = []
    for k, (R1, Z1, R2, Z2) in enumerate(conf["isoflux_targets"]):
        labels.append("iso%d psi(%g,%g)-psi(%g,%g)" % (k, R1, Z1, R2, Z2))
    for R, Z in conf["xpoint_targets"]:
        labels.append("B_R(%g,%g)" % (R, Z))
    for R, Z in conf["xpoint_targets"]:
        labels.append("B_Z(%g,%g)" % (R, Z))
    return labels


# ---------------------------------------------------------------------------
# one solve
# ---------------------------------------------------------------------------
def run_one(conf, nx, arm, args):
    if conf["config"] in INVALID:
        return dict(config=conf["config"], arm=arm, nx=nx, ny=nx,
                    status=INVALID[conf["config"]], converged=False, wall=0.0)
    from freegsnke import build_machine, equilibrium_update, GSstaticsolver
    from freegsnke.inverse import Inverse_optimizer

    g = conf["grid"]
    LR, LZ, clearance = star_limiter(conf)
    limiter = [dict(R=float(r), Z=float(z)) for r, z in zip(LR, LZ)]
    quiet = io.StringIO()
    with contextlib.redirect_stdout(quiet):
        tok = build_machine.tokamak(active_coils_data=active_coils(conf, args.quad,
                                                                     args.as_rect),
                                    limiter_data=limiter)
    free = [s["name"] for s in conf["supplies"] if s["free"]]
    for s in conf["supplies"]:
        cur = 0.0 if (arm == "cold" and s["free"]) else s["reference_current"]
        tok.set_coil_current(s["name"], float(cur))
        tok[s["name"]].control = bool(s["free"])
    tok.getCurrentsVec()

    eq = equilibrium_update.Equilibrium(
        tokamak=tok, Rmin=g["Rmin"], Rmax=g["Rmax"], Zmin=g["Zmin"],
        Zmax=g["Zmax"], nx=nx, ny=nx)
    Prof = make_profile_class(conf)
    profiles = Prof(eq=eq, Ip=float(conf["Ip"]), fvac=float(conf["fvac"]),
                    psi_n=np.linspace(0.0, 1.0, 256),
                    pprime_data=np.ones(256), ffprime_data=np.ones(256),
                    Raxis=1.0, Ip_logic=True)
    solver = GSstaticsolver.NKGSsolver(eq)

    xt = np.asarray(conf["xpoint_targets"], float).reshape(-1, 2)
    iso = [[[R1, R2], [Z1, Z2]] for R1, Z1, R2, Z2 in conf["isoflux_targets"]]
    constrain = Inverse_optimizer(null_points=[list(xt[:, 0]), list(xt[:, 1])]
                                  if len(xt) else None,
                                  isoflux_set=iso if iso else None)
    A, b_vac = vacuum_response(eq, constrain)
    s = np.linalg.svd(A, compute_uv=False)
    reg = float(args.reg*s[0]**2)
    reg_fj = float(args.reg_fj*s[0]**2)

    log = io.StringIO()
    t0 = time.perf_counter()
    status = "ok"
    fwd = None
    if arm == "fwd":
        # freegsnke's ordinary FORWARD solve at the seeded currents, from the
        # same Gaussian, and only then the inverse: the plasma the inverse
        # starts from is consistent with the seed currents.  Uses nothing the
        # reference solved for beyond the currents the seed already is.
        flog = io.StringIO()
        try:
            with contextlib.redirect_stdout(flog):
                solver.solve(eq=eq, profiles=profiles, constrain=None,
                             target_relative_tolerance=args.tol,
                             max_solving_iterations=200, verbose=True)
            fwd = dict(status="ok", rel=float(solver.relative_change),
                       wall=time.perf_counter() - t0,
                       psi_axis=float(profiles.psi_axis)
                       if hasattr(profiles, "psi_axis") else None)
        except Exception as e:                              # noqa: BLE001
            fwd = dict(status="exception: %s: %s" % (type(e).__name__,
                                                       str(e)[:300]))
        open(os.path.join(args.out, "%s-%s-n%d%s-forward.log" % (
            conf["config"], arm, nx, ("-" + args.tag) if args.tag else "")),
             "w").write(flog.getvalue())
    try:
        with contextlib.redirect_stdout(log):
            solver.solve(eq=eq, profiles=profiles, constrain=constrain,
                         target_relative_tolerance=args.tol,
                         target_relative_psit_update=args.psit,
                         max_solving_iterations=args.maxits,
                         l2_reg=reg, l2_reg_fj=reg_fj,
                         use_full_Jacobian=not args.no_fj,
                         damping_factor=args.damping,
                         verbose=True)
    except Exception as e:                                  # noqa: BLE001
        status = "exception: %s: %s" % (type(e).__name__, str(e)[:300])
    wall = time.perf_counter() - t0
    text = log.getvalue()

    rec = dict(config=conf["config"], arm=arm, nx=nx, ny=nx, status=status,
               wall=wall, free=free, l2_reg=reg, l2_reg_fj=reg_fj,
               full_jacobian_allowed=not args.no_fj,
               tol=args.tol, psit=args.psit, quad=args.quad,
               as_rect=args.as_rect, damping=args.damping,
               fgs_rule=args.fgs_rule,
               limiter_clearance=clearance, forward_first=fwd,
               vacuum_singular_values=[float(v) for v in s])
    m = re.search(r"Inverse static solve (SUCCESS|DID NOT CONVERGE)\. Tolerance "
                  r"(\S+) \(vs\. requested (\S+)\) reached in (\d+)/(\d+)", text)
    if m:
        rec.update(freegsnke_verdict=m.group(1), gs_rel_residual=float(m.group(2)),
                   iterations=int(m.group(4)), max_iterations=int(m.group(5)))
    else:
        rec.update(freegsnke_verdict=None, iterations=None)
    rec["full_jacobian_steps"] = text.count("Using full Jacobian")
    rec["green_jacobian_steps"] = text.count("Using simplified Green")
    rec["resizes"] = text.count("Resizing of the control current update")
    losses = re.findall(r"Constraint losses = (\S+)", text)
    rec["loss_history"] = [float(v) for v in losses]
    psit = re.findall(r"Relative update of tokamak psi \(in plasma core\): (\S+)", text)
    rec["psit_history"] = [float(v) for v in psit]
    relerr = re.findall(r"Relative error =\s+(\S+)", text)
    rec["gs_history"] = [float(v) for v in relerr]
    open(os.path.join(args.out, "%s-%s-n%d%s.log" % (
        conf["config"], arm, nx, ("-" + args.tag) if args.tag else "")),
         "w").write(text)

    if status != "ok":
        rec["converged"] = False
        return rec

    cur = eq.tokamak.getCurrents()
    rec["currents"] = {k: float(v) for k, v in cur.items()}
    rec["free_currents"] = {f: float(cur[f]) for f in free}

    # the rows, by freegsnke's own constraint machinery, at the final state
    constrain.build_plasma_vals(eq.plasma_psi)
    constrain.build_lsq(eq.tokamak.getCurrentsVec())
    rows = -np.array(constrain.b)
    rec["rows"] = dict(zip(row_labels(conf), [float(v) for v in rows]))
    # and independently, through the equilibrium's own field evaluation
    indep = []
    for R1, Z1, R2, Z2 in conf["isoflux_targets"]:
        indep.append(float(eq.psiRZ(R1, Z1) - eq.psiRZ(R2, Z2)))
    for R, Z in xt:
        indep.append(float(eq.Br(R, Z)))
    for R, Z in xt:
        indep.append(float(eq.Bz(R, Z)))
    rec["rows_independent"] = dict(zip(row_labels(conf), indep))
    # IS THE END A STATIONARY POINT OF THE LEAST SQUARES?  freegsnke stops on
    # the SIZE of its last current step, and that step is capped by a ratchet
    # (the mean of the last six) and a damping factor**iteration, so a small
    # step need not mean a small gradient.  The cosine between the residual
    # and the range of the last full Jacobian freegsnke built answers it: 0 at
    # a stationary point of min ||rows||, 1 if the rows could still all fall.
    if getattr(solver, "dbdI", None) is not None:
        J = np.array(solver.dbdI)
        g = J.T @ rows
        rec["stationarity_cosine"] = float(
            np.linalg.norm(np.linalg.lstsq(J, rows, rcond=None)[0] @ J.T)
            / max(np.linalg.norm(rows), 1e-300))
        rec["gauss_newton_step_rel"] = float(
            np.max(np.abs(np.linalg.lstsq(J, -rows, rcond=None)[0]))
            / max(abs(v) for v in rec["free_currents"].values()))
        rec["gradient_norm"] = float(np.linalg.norm(g))

    # the equilibrium: re-evaluate the profile on the final psi so opt/xpt and
    # the boundary belong to the state reported
    psi = eq.psi()
    profiles.Jtor(eq.R, eq.Z, psi)
    opt, xpt = np.atleast_2d(profiles.opt), np.atleast_2d(profiles.xpt)
    rec["psi_axis"] = float(opt[0][2])
    rec["psi_bndry"] = float(profiles.psi_bndry)
    rec["axis"] = [float(opt[0][0]), float(opt[0][1])]
    rec["xpoints"] = [[float(v) for v in row] for row in xpt] if xpt.size else []
    rec["flag_limiter"] = bool(profiles.flag_limiter)
    rec["Ip_scale_L"] = float(profiles.L)
    rec["Ip"] = float(np.sum(profiles.jtor)*eq.dR*eq.dZ)
    # refine every X-point candidate to a zero of the total field, starting
    # from freegsnke's own location, so "located X-point" is a field zero
    rec["xpoint_located"] = [locate_null(eq, *t) for t in xt]
    gs_ok = (rec.get("gs_rel_residual") is not None
             and rec["gs_rel_residual"] <= args.tol)
    it_ok = (rec.get("iterations") is not None
             and rec["iterations"] < rec["max_iterations"])
    rec["converged"] = bool(gs_ok and it_ok)
    np.savez_compressed(os.path.join(args.out, "%s-%s-n%d%s.npz" % (
        conf["config"], arm, nx, ("-" + args.tag) if args.tag else "")), R=eq.R[:, 0], Z=eq.Z[0, :], psi=psi,
        plasma_psi=eq.plasma_psi, jtor=profiles.jtor)
    return rec


def locate_null(eq, R, Z, its=40):
    """Newton on (B_R, B_Z) = 0 from (R, Z), derivatives by central difference
    of the equilibrium's own field; returns [R, Z, psi] or None."""
    x = np.array([R, Z], float)
    h = 1e-5
    for _ in range(its):
        f = np.array([eq.Br(*x), eq.Bz(*x)], float).ravel()
        J = np.zeros((2, 2))
        for k in range(2):
            d = np.zeros(2)
            d[k] = h
            J[:, k] = (np.array([eq.Br(*(x + d)), eq.Bz(*(x + d))]).ravel()
                       - np.array([eq.Br(*(x - d)), eq.Bz(*(x - d))]).ravel())/(2*h)
        try:
            dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            return None
        x += dx
        if np.hypot(*dx) < 1e-11:
            break
    if np.hypot(x[0] - R, x[1] - Z) > 0.2:
        return None
    return [float(x[0]), float(x[1]), float(eq.psiRZ(*x))]


# ---------------------------------------------------------------------------
# the driver and the table
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--configs", default="")
    ap.add_argument("--nx", default="129,257")
    ap.add_argument("--arms", default="seeded")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--psit", type=float, default=1e-9)
    ap.add_argument("--maxits", type=int, default=200)
    ap.add_argument("--reg", type=float, default=1e-8)
    ap.add_argument("--reg-fj", dest="reg_fj", type=float, default=1e-8)
    ap.add_argument("--quad", type=int, default=6)
    ap.add_argument("--no-fj", dest="no_fj", action="store_true")
    ap.add_argument("--as-rect", dest="as_rect", action="store_true")
    ap.add_argument("--fgs-rule", dest="fgs_rule", action="store_true")
    ap.add_argument("--damping", type=float, default=0.98,
                    help="freegsnke's damping_factor; 0.98 is solve()'s default")
    ap.add_argument("--tag", default="")
    ap.add_argument("--timeout", type=float, default=5400.0)
    ap.add_argument("--one", default="")
    ap.add_argument("--arm", default="seeded")
    ap.add_argument("--table", action="store_true")
    args = ap.parse_args()
    args.out = os.path.abspath(args.out)
    global FGS_RULE
    FGS_RULE = args.fgs_rule
    os.makedirs(args.out, exist_ok=True)
    P = load_problem()

    if args.one:
        rec = run_one(P[args.one], int(args.nx), args.arm, args)
        rec["tag"] = args.tag
        name = "%s-%s-n%s%s.json" % (args.one, args.arm, args.nx,
                                     ("-" + args.tag) if args.tag else "")
        json.dump(rec, open(os.path.join(args.out, name), "w"), indent=1)
        print("DONE %s %s %s converged=%s status=%s wall=%.1f" % (
            args.one, args.arm, args.nx, rec.get("converged"), rec["status"],
            rec["wall"]), flush=True)
        return 0

    if args.table:
        return table(args, P)

    want = [c for c in P if not args.configs
            or any(c.endswith(w) or c == w for w in args.configs.split(","))]
    jobs = [(c, int(n), a) for c in want for n in args.nx.split(",")
            for a in args.arms.split(",")]
    env = {**os.environ, **THREAD_ENV, "MPLBACKEND": "Agg"}
    passthrough = ["--tol", str(args.tol), "--psit", str(args.psit),
                   "--maxits", str(args.maxits), "--reg", str(args.reg),
                   "--reg-fj", str(args.reg_fj), "--quad", str(args.quad),
                   "--damping", str(args.damping)]
    if args.no_fj:
        passthrough.append("--no-fj")
    if args.as_rect:
        passthrough.append("--as-rect")
    if args.fgs_rule:
        passthrough.append("--fgs-rule")
    if args.tag:
        passthrough += ["--tag", args.tag]
    running = []
    while jobs or running:
        while jobs and len(running) < args.jobs:
            c, n, a = jobs.pop(0)
            cmd = [sys.executable, os.path.abspath(__file__), args.out,
                   "--one", c, "--nx", str(n), "--arm", a] + passthrough
            errf = open(os.path.join(args.out, "%s-%s-n%d%s.stderr" % (
                c, a, n, ("-" + args.tag) if args.tag else "")), "w")
            running.append((subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                             stderr=errf, text=True),
                            time.time(), c, n, a, errf))
        time.sleep(2)
        for item in list(running):
            p, t0, c, n, a, errf = item
            if p.poll() is not None:
                out = p.stdout.read().strip().splitlines()
                print(out[-1] if out else "%s %s %d: exit %d, no output"
                      % (c, a, n, p.returncode), flush=True)
                errf.close()
                running.remove(item)
            elif time.time() - t0 > args.timeout:
                p.kill()
                print("%s %s %d: KILLED after %.0f s" % (c, a, n, args.timeout),
                      flush=True)
                errf.close()
                running.remove(item)
                json.dump(dict(config=c, arm=a, nx=n, status="timeout after %.0f s"
                               % args.timeout, converged=False, tag=args.tag),
                          open(os.path.join(args.out, "%s-%s-n%d%s.json" % (
                              c, a, n, ("-" + args.tag) if args.tag else "")),
                               "w"), indent=1)
    print("NKE DONE", flush=True)
    return 0


NUM = r"([-+]?\d+\.\d+e[-+]\d+)"


def meq_state(config):
    """MEQ's psi_ax, psi_bnd, axis and X-point from its inverse log."""
    path = os.path.join(MEQ_RUNS, config + "-inverse.log")
    if not os.path.exists(path):
        return None
    text = open(path).read()
    out = {}
    m = re.search(r"psi_ax = " + NUM, text)
    out["psi_axis"] = float(m.group(1)) if m else None
    m = re.search(r"the axis, as a zero of q_h: psi = " + NUM
                  + r" at \(\s*([-\d.]+),\s*([-\d.]+)\s*\)", text)
    out["axis"] = [float(m.group(2)), float(m.group(3))] if m else None
    m = re.search(r"psi_bnd = " + NUM + r".*?\(\s*([-\d.]+),\s*([-\d.]+)\s*\)",
                  text)
    out["psi_bndry"] = float(m.group(1)) if m else None
    out["xpoint"] = [float(m.group(2)), float(m.group(3))] if m else None
    return out


def reference_nulls(conf):
    """The reference's X-points as zeros of ITS OWN field, by Newton on a
    quintic spline of its psi from each design target: freegs4e's stored
    X-points come from its grid critical-point finder and sit about 1e-4 m
    from the null its own rows put at the target."""
    from scipy.interpolate import RectBivariateSpline
    d = np.load(conf["reference_npz"], allow_pickle=True)
    sp = RectBivariateSpline(d["R"], d["Z"], d["psi"], kx=5, ky=5)
    out = []
    for R, Z in conf["xpoint_targets"]:
        x = np.array([R, Z], float)
        for _ in range(50):
            g = np.array([sp(*x, dx=1, grid=False), sp(*x, dy=1, grid=False)])
            H = np.array([[sp(*x, dx=2, grid=False), sp(*x, dx=1, dy=1, grid=False)],
                          [sp(*x, dx=1, dy=1, grid=False), sp(*x, dy=2, grid=False)]])
            dx = np.linalg.solve(H, -g)
            x += dx
            if np.hypot(*dx) < 1e-12:
                break
        out.append([float(x[0]), float(x[1])])
    return out


def table(args, P):
    if not args.tag:
        # results.json: EVERY run record, every arm and every variant tag, in
        # one file beside problem.json's directory as the brief asks
        allrec = []
        for f in sorted(os.listdir(args.out)):
            if f.endswith(".json") and not f.startswith(("results", "table")):
                r = json.load(open(os.path.join(args.out, f)))
                r.setdefault("tag", "")
                r["resolution"] = "%dx%d" % (r["nx"], r.get("ny", r["nx"]))
                r["exit"] = r.get("status")
                allrec.append(r)
        json.dump(allrec, open(os.path.join(args.out, "results.json"), "w"),
                  indent=1)
    recs = []
    for f in sorted(os.listdir(args.out)):
        if f.endswith(".json") and not f.startswith(("results", "table")):
            r = json.load(open(os.path.join(args.out, f)))
            if (r.get("tag") or "") == args.tag:
                recs.append(r)
    by = {}
    for r in recs:
        by[(r["config"], r["arm"], r["nx"])] = r
    rows = []
    for name, conf in P.items():
        ref = conf["reference"]
        free = [s["name"] for s in conf["supplies"] if s["free"]]
        Iref = {s["name"]: s["reference_current"] for s in conf["supplies"]}
        scale = max(abs(Iref[f]) for f in free)
        Imeq = conf["meq"].get("solved") or {}
        # MEQ's absolute state only where problem.json still carries its run
        meq = (meq_state(name) or {}) if Imeq else {}
        for arm in ("seeded", "fwd", "cold"):
            ns = sorted(n for (c, a, n) in by if c == name and a == arm)
            for n in ns:
                r = by[(name, arm, n)]
                row = dict(config=name, arm=arm, nx=n, converged=r.get("converged"),
                           status=r.get("status"), iterations=r.get("iterations"),
                           wall=r.get("wall"), flag_limiter=r.get("flag_limiter"))
                if r.get("free_currents"):
                    I = r["free_currents"]
                    row["dI_fgs"] = max(abs(I[f] - Iref[f]) for f in free)/scale
                    row["dI_meq"] = (max(abs(I[f] - Imeq[f]) for f in free)/scale
                                     if Imeq else None)
                    row["dI_meq_fgs"] = (max(abs(Imeq[f] - Iref[f]) for f in free)
                                         / scale if Imeq else None)
                    other = [m for m in ns if m != n]
                    if other and by[(name, arm, other[-1])].get("free_currents"):
                        J = by[(name, arm, max(ns))]["free_currents"]
                        K = by[(name, arm, min(ns))]["free_currents"]
                        row["dI_res"] = max(abs(J[f] - K[f]) for f in free)/scale
                    row["psi_ax_fgs"] = abs(r["psi_axis"] - ref["psi_axis"])/abs(ref["psi_axis"])
                    row["psi_ax_meq"] = (abs(r["psi_axis"] - meq["psi_axis"])
                                         / abs(ref["psi_axis"])
                                         if meq.get("psi_axis") is not None else None)
                    row["axis_fgs"] = float(np.hypot(r["axis"][0] - ref["Raxis"],
                                                     r["axis"][1] - ref["Zaxis"]))
                    row["axis_meq"] = (float(np.hypot(r["axis"][0] - meq["axis"][0],
                                                      r["axis"][1] - meq["axis"][1]))
                                       if meq.get("axis") else None)
                    xl = [x for x in r.get("xpoint_located") or [] if x]
                    refx = np.asarray(reference_nulls(conf), float)
                    if xl and refx.size:
                        # each located null against the nearest reference null
                        row["X_fgs"] = max(float(np.min(np.hypot(refx[:, 0] - x[0],
                                                                 refx[:, 1] - x[1])))
                                           for x in xl)
                    tx = np.asarray(conf["xpoint_targets"], float)
                    if xl and tx.size:
                        row["X_target"] = max(float(np.min(np.hypot(tx[:, 0] - x[0],
                                                                    tx[:, 1] - x[1])))
                                              for x in xl)
                    if xl and meq.get("xpoint"):
                        row["X_meq"] = min(float(np.hypot(x[0] - meq["xpoint"][0],
                                                          x[1] - meq["xpoint"][1]))
                                           for x in xl)
                    row["max_row"] = max(abs(v) for v in r["rows"].values())
                    row["meq_dI_fgs_published"] = conf["meq"].get("dI")
                    row["meq_exit"] = conf["meq"].get("exit")
                rows.append(row)
    json.dump(rows, open(os.path.join(args.out, "table%s.json" % (
        ("-" + args.tag) if args.tag else "")), "w"), indent=1)
    fmt = lambda v, f="%9.2e": (f % v) if isinstance(v, (int, float)) and v is not None else "%9s" % "-"
    print("%-28s %-6s %4s %4s %4s | %9s %9s %9s %9s | %9s %9s %9s %9s %9s %9s %9s" % (
        "config", "arm", "nx", "conv", "its", "nke-fgs", "nke-meq", "meq-fgs",
        "nke res", "psax fgs", "psax meq", "ax fgs", "ax meq", "X fgs", "X meq",
        "max row"))
    for r in rows:
        print("%-28s %-6s %4d %4s %4s | %s %s %s %s | %s %s %s %s %s %s %s%s" % (
            r["config"], r["arm"], r["nx"], "y" if r["converged"] else "N",
            r.get("iterations") if r.get("iterations") is not None else "-",
            fmt(r.get("dI_fgs")), fmt(r.get("dI_meq")), fmt(r.get("dI_meq_fgs")),
            fmt(r.get("dI_res")), fmt(r.get("psi_ax_fgs")), fmt(r.get("psi_ax_meq")),
            fmt(r.get("axis_fgs")), fmt(r.get("axis_meq")), fmt(r.get("X_fgs")),
            fmt(r.get("X_meq")), fmt(r.get("max_row")),
            "" if not r.get("flag_limiter") else "  LIMITED"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
