#!/usr/bin/env python
"""EXPORT geq'S ANISOTROPIC ROTATING MIRROR, FOR ANISOTROPY-PLAN.md AN-4.

MR-4's reference, export_rotating.py's -- the same machine, profiles, labels,
direct freeBoundary, phi_0 tightened every Picard step and N_D held -- with a
PRESCRIBED pressure anisotropy handed to freegs4e's ProfilesCentrifugalMirror
through geq's PrescribedAnisotropy:

    Delta( psi_n, z ) = delta0 g( psi_n ) h( z ),
    g = ( 4 psi_n ( 1 - psi_n ) )^3   on 0 < psi_n < 1, zero outside,
    h = ( 1 - ( z/1.2 )^2 )^3         on |z| < 1.2,     zero outside,

dp_par_dpsin = fast_fraction = 0, which isolates the operator: geq's Jtor is
then the isotropic bracket with the tension folded in, and MEQ's source is
unchanged. g is symmetric under psi_n -> 1 - psi_n, so the same callable is
MEQ's Delta( Psi ), Psi = 1 - psi_n, and it vanishes at both labels to
second order; h vanishes before geq's domain edge and MEQ's box end at
|z| = 1.5. That is the plan's fact 4: a Delta non-zero at either would be a
coefficient jump in MEQ and invisible to geq. tests/convergence/
AnisotropicMirror.cpp's MirrorProfile is the same function.

BOTH OF freegs4e'S REARRANGEMENTS, each lagged one Picard step:

  historical   Jtor_eff = ( Jtor - grad psi . grad Delta/( mu0 R ) )/( 1 - Delta )
  divergence   Jtor_eff = Jtor - div( ( Delta/R ) grad psi )/mu0

with np.gradient for grad Delta and the divergence. They agree only at
convergence and only as the grid resolves Delta's gradient, so their
difference is a control on geq's own discretisation of it.

WHAT IS WRITTEN. geq's answer, to --out (default
runs/geq-benchmark/reference/anisotropic/, untracked):
psi_<form>_<delta0>_<n>.dat on the 65 x 65 common nodes, export_rotating's
layout, and solution.txt. The PROBLEM is reference/rotating/'s tables and
filaments plus meta.txt here, written to --tables (default --out) and
reported against the tracked tools/geq-benchmark/reference/anisotropic/.

AND THE THREE CHECKS AN-4 IS ACCEPTED ON, printed and in solution.txt:

  (a) delta0 = 0 through the hook is M-261's solve: bit for bit against the
      hook absent, in either form, and against M-261's psi_direct_<n>.dat if
      runs/geq-benchmark/reference/rotating/ holds it;
  (b) geq's self-convergence at each delta0 and form, against the finest;
  (c) the two rearrangements' difference at each resolution, which must fall.

AGAINST THE LATEST geq AND THE freegs4e IT RUNS ON, the tip of
feature/rotation_equilibrium (README.md, Regenerating the reference).
meta.txt records what was imported.

    PYTHONPATH=<freegs4e at its tip>:/home/ian/projects/geq/src \\
    tools/freegs4e-benchmark/venv/bin/python tools/geq-benchmark/export_anisotropic.py

About an hour, most of it the four 513^2 solves. `--cache DIR` keeps each
solve and reuses it.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import export_rotating as rot  # noqa: E402  -- installs MR-4's monkeypatches
from export_mirror import answers_dir, TRACKED  # noqa: E402


class _NoUnits:
    """geq.fast_species imports `units` from geq.utils, which imports xarray
    and pint, and the benchmark venv has neither (README.md): installing them
    could move the numpy every reference here was taken with. Only the
    slowing-down classes use `units`, and PrescribedAnisotropy is not one of
    them, so geq.utils is stood in for by this, which refuses any use."""

    def __getattr__(self, name):
        raise RuntimeError(f"geq.utils.units.{name}: pint is not in the benchmark "
                           "venv, and export_anisotropic.py uses no unit")


if "geq.utils" not in sys.modules:
    import types  # noqa: E402

    _utils = types.ModuleType("geq.utils")
    _utils.units = _NoUnits()
    sys.modules["geq.utils"] = _utils

from geq.fast_species import PrescribedAnisotropy  # noqa: E402

Z_CUT = 1.2
FORMS = ("historical", "divergence")


def g(psi_n):
    x = np.asarray(psi_n, dtype=float)
    inside = (x > 0.0) & (x < 1.0)
    u = np.where(inside, 4.0 * x * (1.0 - x), 0.0)
    return u ** 3


def h(z):
    v = np.asarray(z, dtype=float) / Z_CUT
    inside = np.abs(v) < 1.0
    w = np.where(inside, 1.0 - v * v, 0.0)
    return w ** 3


class Anisotropy:
    """One posed Delta: an amplitude and a rearrangement. export_rotating.build()
    asks it for freegs4e's anisotropy object; solve() for a cache tag."""

    def __init__(self, delta0, form):
        if form not in FORMS:
            raise ValueError(f"form {form!r} is not one of {FORMS}")
        self.delta0 = float(delta0)
        self.form = form
        self.divergence = form == "divergence"
        self.last = None

    def delta(self, psi_n, R, Z, B):  # noqa: ARG002 -- geq's signature
        return self.delta0 * g(psi_n) * h(Z)

    def prescribed(self):
        # Paired zeros: dp_par_dpsin and fast_fraction together or neither, and
        # fast_fraction = 0 leaves the Maxwellian bracket unscaled.
        self.last = PrescribedAnisotropy(delta=self.delta, dp_par_dpsin=0.0,
                                         fast_fraction=0.0)
        return self.last

    def tag(self):
        return f"aniso_{self.form}_d{self.delta0:+.3f}"

    def name(self):
        return f"{self.form}_{self.delta0:+.1f}"


def max_abs(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


def read_psi(path):
    data = np.loadtxt(path)
    return data[:, 2].reshape(rot.COARSE, rot.COARSE)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=None,
                    help="geq's answer; default runs/geq-benchmark/reference/anisotropic")
    ap.add_argument("--tables", default=None,
                    help="the problem, meta.txt; default --out. "
                         "tools/geq-benchmark/reference/anisotropic refreshes the tracked copy")
    ap.add_argument("--deltas", default="-0.1,-0.3")
    ap.add_argument("--resolutions", default="129,257,513")
    ap.add_argument("--zero-at", default="129,257",
                    help="resolutions at which delta0 = 0 is checked against the hook absent")
    ap.add_argument("--mach", type=float, default=2.0)
    ap.add_argument("--rtol", type=float, default=2e-9)
    ap.add_argument("--maxits", type=int, default=400)
    ap.add_argument("--cache", default=None)
    args = ap.parse_args(argv)
    args.out = args.out or answers_dir("anisotropic")
    args.tables = args.tables or args.out

    # MR-4's settings, exactly: the skew shape and its Mach radius.
    rot.SHAPE["name"] = "skew"
    rot.R_MACH["value"] = rot.SHAPES["skew"][1]
    omega_peak = rot.omega_peak_for(args.mach)
    deltas = [float(x) for x in args.deltas.split(",")]
    ns = sorted(int(x) for x in args.resolutions.split(","))
    zero_ns = sorted(int(x) for x in args.zero_at.split(",")) if args.zero_at else []
    for n in ns + zero_ns:
        if (n - 1) % (rot.COARSE - 1):
            raise SystemExit(f"resolution {n} does not contain the {rot.COARSE}-grid")

    def run(n, anisotropy):
        eq, prof, rec = rot.solve(n, "direct", omega_peak, args.rtol, args.maxits,
                                  args.cache, anisotropy=anisotropy)
        calls = anisotropy.last.n_calls if anisotropy is not None and anisotropy.last else 0
        label = "isotropic" if anisotropy is None else anisotropy.name()
        print(f"{label:18s} {n:4d}^2  {rec['wall']:8.1f} s  {rec['iterations']:3d} its"
              f"  relchange {rec['relchange']:.2e}  psi_axis {rec['psi_axis']:.12e}"
              f"  psi_bndry {rec['psi_bndry']:.12e}  Delta evaluations {calls}"
              f"{'  (cached)' if rec['cached'] else ''}", flush=True)
        return rec

    # --- (a) delta0 = 0 through the hook is M-261's solve --------------------
    print("(a) delta0 = 0 through the hook against the hook absent")
    zero = {}
    m261 = answers_dir("rotating")
    for n in zero_ns:
        plain = run(n, None)
        row = {}
        for form in FORMS:
            r = run(n, Anisotropy(0.0, form))
            row[form] = (max_abs(r["psi65"], plain["psi65"]),
                         r["psi_axis"] - plain["psi_axis"],
                         r["psi_bndry"] - plain["psi_bndry"])
        path = os.path.join(m261, f"psi_direct_{n}.dat")
        row["M-261"] = (max_abs(plain["psi65"], read_psi(path)), None, None) \
            if os.path.exists(path) else None
        zero[n] = row
        for key, v in row.items():
            if v is None:
                print(f"  {n:4d}^2  against M-261's file: absent, {path}")
                continue
            print(f"  {n:4d}^2  {key:10s} max|dpsi| {v[0]:.3e}"
                  + ("" if v[1] is None else f"  dpsi_axis {v[1]:+.3e}  dpsi_bndry {v[2]:+.3e}"))

    # --- the anisotropic solves ----------------------------------------------
    runs = {}
    for d in deltas:
        for form in FORMS:
            for n in ns:
                runs[(d, form, n)] = run(n, Anisotropy(d, form))
    finest = ns[-1]
    iso = run(finest, None)
    path = os.path.join(m261, f"psi_direct_{finest}.dat")
    if os.path.exists(path):
        zero.setdefault(finest, {})["M-261"] = (max_abs(iso["psi65"], read_psi(path)), None, None)
        print(f"  {finest:4d}^2  isotropic against M-261's file: max|dpsi|"
              f" {zero[finest]['M-261'][0]:.3e}")

    # --- (b) self-convergence --------------------------------------------------
    print(f"\n(b) self-convergence on the 65-grid nodes, against {finest}^2")
    conv = {}
    for d in deltas:
        for form in FORMS:
            ref = runs[(d, form, finest)]
            for n in ns[:-1]:
                r = runs[(d, form, n)]
                conv[(d, form, n)] = (max_abs(r["psi65"], ref["psi65"]),
                                      r["psi_axis"] - ref["psi_axis"],
                                      r["psi_bndry"] - ref["psi_bndry"])
                c = conv[(d, form, n)]
                print(f"  delta0 {d:+.1f} {form:10s} {n:4d}: max|dpsi| {c[0]:.3e}"
                      f"  dpsi_axis {c[1]:+.3e}  dpsi_bndry {c[2]:+.3e}")
            if len(ns) >= 3:
                a, b, c = (runs[(d, form, n)] for n in ns[-3:])
                e1, e2 = max_abs(a["psi65"], b["psi65"]), max_abs(b["psi65"], c["psi65"])
                print(f"  delta0 {d:+.1f} {form:10s} successive {e1:.3e}, {e2:.3e}:"
                      f" order {np.log2(e1 / e2):.2f}")

    # --- (c) the rearrangements' difference -----------------------------------
    print("\n(c) historical minus divergence, at each resolution")
    between = {}
    for d in deltas:
        for n in ns:
            a, b = runs[(d, "historical", n)], runs[(d, "divergence", n)]
            between[(d, n)] = (max_abs(a["psi65"], b["psi65"]),
                               a["psi_axis"] - b["psi_axis"],
                               a["psi_bndry"] - b["psi_bndry"])
            v = between[(d, n)]
            print(f"  delta0 {d:+.1f} {n:4d}^2: max|dpsi| {v[0]:.3e}"
                  f"  dpsi_axis {v[1]:+.3e}  dpsi_bndry {v[2]:+.3e}")

    # --- what the anisotropy moves ----------------------------------------------
    print(f"\nwhat Delta moves, at {finest}^2, against the isotropic solve")
    signal = {}
    for d in deltas:
        for form in FORMS:
            r = runs[(d, form, finest)]
            signal[(d, form)] = (max_abs(r["psi65"], iso["psi65"]),
                                 r["psi_axis"] - iso["psi_axis"],
                                 r["psi_bndry"] - iso["psi_bndry"])
            v = signal[(d, form)]
            print(f"  delta0 {d:+.1f} {form:10s} max|dpsi| {v[0]:.3e}"
                  f"  dpsi_axis {v[1]:+.3e}  dpsi_bndry {v[2]:+.3e}")
    plasma_scale = float(np.max(np.abs(iso["plasma65"])))
    print(f"  max |plasma psi|, isotropic, {plasma_scale:.4e}")

    # --- write -------------------------------------------------------------------
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.tables, exist_ok=True)
    for (d, form, n), r in runs.items():
        rot.write_psi(os.path.join(args.out, f"psi_{form}_{d:+.1f}_{n}.dat"), r,
                      f"direct freeBoundary, {n} x {n}, delta0 = {d:+.1f}, {form} form")

    geq_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(sys.modules["geq.solver"].__file__))))
    fgs_dir = os.path.dirname(os.path.dirname(
        os.path.abspath(sys.modules["freegs4e.species"].__file__)))
    with open(os.path.join(args.tables, "meta.txt"), "w") as fh:
        fh.write("case anisotropic\n")
        fh.write(f"geq_commit {rot.git_commit(geq_dir)} {rot.git_dirty(geq_dir)}\n")
        fh.write(f"freegs4e_commit {rot.git_commit(fgs_dir)} {rot.git_dirty(fgs_dir)}\n")
        fh.write("# The problem is reference/rotating/'s -- meta.txt, the tables and\n")
        fh.write("# filaments.dat -- with the anisotropy below. Delta is in geq's\n")
        fh.write("# psi_n = 1 - Psi; g is symmetric, so it reads the same in MEQ's Psi.\n")
        fh.write("profiles reference/rotating\n")
        fh.write("delta delta0*(4*psi_n*(1-psi_n))^3*(1-(z/z_cut)^2)^3,"
                 " zero outside 0<psi_n<1 and |z|<z_cut\n")
        fh.write(f"z_cut {Z_CUT:.17e}\n")
        fh.write("delta0 " + " ".join(f"{d:+.17e}" for d in deltas) + "\n")
        fh.write("forms " + " ".join(FORMS) + "\n")
        fh.write("dp_par_dpsin 0.0\nfast_fraction 0.0\n")
        fh.write(f"picard atol 0 rtol {args.rtol:.3e} maxits {args.maxits}\n")
        fh.write("boundary freeBoundary ( direct )\n")
        fh.write("phi0 tightened_every_picard_step\n")
        fh.write("N_s held_as_prescribed_flux_functions\n")
        fh.write(f"resolutions {' '.join(str(n) for n in ns)}\n")

    with open(os.path.join(args.out, "solution.txt"), "w") as fh:
        fh.write("# geq's converged anisotropic states for meta.txt's problem. An\n")
        fh.write("# OUTPUT of export_anisotropic.py, beside psi_*.dat; not tracked.\n")
        for (d, form, n), r in sorted(runs.items()):
            fh.write(f"labels_{form}_{d:+.1f}_{n} {r['psi_axis']:.17e} {r['psi_bndry']:.17e}"
                     f" iterations {r['iterations']} relchange {r['relchange']:.3e}\n")
        fh.write(f"labels_isotropic_{finest} {iso['psi_axis']:.17e} {iso['psi_bndry']:.17e}\n")
        for n, row in zero.items():
            for key, v in row.items():
                if v is not None:
                    fh.write(f"zero_{key}_{n} {v[0]:.3e}\n")
        for (d, form, n), v in sorted(conv.items()):
            fh.write(f"selfconv_{form}_{d:+.1f}_{n} {v[0]:.3e} {v[1]:+.3e} {v[2]:+.3e}\n")
        for (d, n), v in sorted(between.items()):
            fh.write(f"forms_{d:+.1f}_{n} {v[0]:.3e} {v[1]:+.3e} {v[2]:+.3e}\n")
        for (d, form), v in sorted(signal.items()):
            fh.write(f"signal_{form}_{d:+.1f} {v[0]:.3e} {v[1]:+.3e} {v[2]:+.3e}\n")

    tracked = os.path.join(TRACKED, "anisotropic", "meta.txt")
    written = os.path.join(args.tables, "meta.txt")
    if os.path.realpath(written) != os.path.realpath(tracked):
        if not os.path.exists(tracked):
            print(f"no tracked meta.txt at {tracked}; --tables {os.path.dirname(tracked)} writes it")
        elif open(tracked).read() != open(written).read():
            print(f"THE PROBLEM MOVED: meta.txt differs from the tracked {tracked}")
        else:
            print("problem meta.txt identical to the tracked copy")
    print(f"answer written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
