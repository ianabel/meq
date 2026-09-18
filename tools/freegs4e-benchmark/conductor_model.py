#!/usr/bin/env python
"""HOW MUCH OF THE IRREDUCIBLE BENCHMARK ERROR IS THE CONDUCTOR MODEL?

COIL-SUBTRACTION-PLAN.md's CS-0, and it decides whether that plan is worth
building.  No solver runs here and no mesh exists: both fields are
Green's-function sums over the SAME currents, differing only in how the current
is distributed.

  filament   freegs4e's default Coil -- `controlPsi` is `Greens( R, Z )*turns`,
             a POINT source at the coil's centre with a log singularity.
  rectangle  MEQ's meq::Coil -- a uniform current density over the cross-section
             the mesh is aligned to, integrated by quadrature.

M-111 measures F DIII-D's relative error FLAT at 5.785e-03 across a 16x range in
dofs, off-coil flat at 6.3e-04 to 7.6e-04, and concludes *"that is the conductor
model ... it is not a discretisation error and no rung buys it down"*.  THAT
CONCLUSION HAS NEVER BEEN MEASURED DIRECTLY -- it is an inference from the error
refusing to refine.  This measures it.

  * about 5.785e-03  -> the model is the whole story. Taking the coils out of
    the mesh and computing their field analytically lets MEQ adopt freegs4e's
    own filament positions, and the floor becomes a modelling CHOICE.
  * much smaller     -> the flat error is something else, most likely
    REGULARITY, and the subtraction still helps but for a different reason.
  * much larger      -> neither story holds and the plan is premature.

THE NORM IS compare.py'S, so the number lands on the same scale as M-111's: the
relative L2 is against `max | psi_ref |` over the compared nodes, and the
conductor interiors are reported SEPARATELY rather than dropped -- a comparison
that quietly threw away the nodes where it disagrees would be the instrument
choosing the answer.

INSIDE a conductor the rectangle's kernel is log-singular and the quadrature is
the limiting error rather than the model, so that row is an upper bound and is
labelled as one.  The row that answers the question is the off-conductor one,
which is also the one M-111 reports.
"""
import argparse, json, os, re, sys
import numpy as np


def read_coils(toml_path):
    """The `[[coils]]` blocks: ( name, R, Z, halfWidth, halfHeight, current )."""
    text = open(toml_path).read()
    out = []
    for block in text.split("[[coils]]")[1:]:
        def field(key):
            m = re.search(rf"^{key}\s*=\s*([-+0-9.eE]+)", block, re.M)
            return float(m.group(1)) if m else None
        name = re.search(r'^Name\s*=\s*"([^"]*)"', block, re.M)
        values = [field(k) for k in
                  ("CentreR", "CentreZ", "HalfWidth", "HalfHeight", "Current")]
        if all(v is not None for v in values):
            out.append((name.group(1) if name else "?", *values))
    return out


def psi_filament(Rc, Zc, current, RR, ZZ):
    from freegs4e.gradshafranov import Greens
    return current * Greens(Rc, Zc, RR, ZZ)


def psi_rectangle(Rc, Zc, hw, hh, current, RR, ZZ, order):
    """Uniform current density over the rectangle, Gauss-Legendre tensor rule."""
    from freegs4e.gradshafranov import Greens
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = np.zeros_like(RR)
    wsum = 0.0
    for a, wa in zip(nodes, weights):
        for b, wb in zip(nodes, weights):
            total += wa * wb * Greens(Rc + hw * a, Zc + hh * b, RR, ZZ)
            wsum += wa * wb
    return current * total / wsum


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("toml", help="the MEQ config whose [[coils]] are the rectangles")
    ap.add_argument("npz", help="the reference, for the grid and the error scale")
    ap.add_argument("--order", type=int, default=24,
                    help="Gauss points per direction over each rectangle")
    args = ap.parse_args()

    ref = np.load(args.npz, allow_pickle=True)
    R = np.asarray(ref["R"], float)
    Z = np.asarray(ref["Z"], float)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    scale = float(np.abs(np.asarray(ref["psi"], float)).max())

    coils = read_coils(args.toml)
    if not coils:
        print(f"no [[coils]] in {args.toml}", file=sys.stderr)
        return 2
    print(f"{len(coils)} conductors, grid {RR.shape}, "
          f"scale max|psi_ref| = {scale:.6e}, Gauss {args.order}^2 per rectangle\n")

    fil = np.zeros_like(RR)
    rec = np.zeros_like(RR)
    inside = np.zeros(RR.shape, bool)
    for name, Rc, Zc, hw, hh, current in coils:
        fil += psi_filament(Rc, Zc, current, RR, ZZ)
        rec += psi_rectangle(Rc, Zc, hw, hh, current, RR, ZZ, args.order)
        inside |= (np.abs(RR - Rc) <= hw) & (np.abs(ZZ - Zc) <= hh)

    diff = rec - fil
    def norms(mask, label, note=""):
        if mask.sum() == 0:
            print(f"  {label:26s} (no nodes)"); return None
        d = diff[mask]
        rel2 = float(np.sqrt(np.mean(d ** 2)) / scale)
        reli = float(np.abs(d).max() / scale)
        print(f"  {label:26s} {int(mask.sum()):8d} {rel2:12.3e} {reli:12.3e}  {note}")
        return rel2

    print(f"  {'region':26s} {'nodes':>8} {'rel L2':>12} {'rel Linf':>12}")
    finite = np.isfinite(diff)
    whole = norms(finite, "everywhere")
    off = norms(finite & ~inside, "off the conductors",
                "<- the row M-111 reports")
    norms(finite & inside, "inside them",
          "upper bound: log kernel, quadrature limited")

    print(f"\n  M-111's F DIII-D figures, for comparison:"
          f"\n    rel L2 5.785e-03 flat across a 16x range in dofs"
          f"\n    off the coils 6.3e-04 to 7.6e-04, also flat")
    if off is not None and whole is not None:
        print(f"\n  VERDICT on COIL-SUBTRACTION-PLAN.md CS-0, ON BOTH ROWS --")
        print(f"  reading only the off-conductor one is how the first version of")
        print(f"  this script got the answer backwards.")
        g = whole / 5.785e-03
        o = off / 7.0e-04
        print(f"\n    GLOBAL   {whole:.3e} against M-111's 5.785e-03  ->  "
              f"{g:5.2f}x")
        print(f"    OFF-COIL {off:.3e} against M-111's ~7.0e-04  ->  {o:5.2f}x")
        if 0.5 <= g <= 2.0:
            print(f"\n    The GLOBAL irreducible error IS the conductor model. A figure")
            print(f"    that refused to refine across a 16x range in dofs is matched to")
            print(f"    {abs(1-g)*100:.0f} per cent by two Green's-function sums with no solver")
            print(f"    anywhere, which is as direct as this gets.")
        if o < 0.6:
            print(f"\n    OFF the conductors it explains only about {o*100:.0f} per cent, so")
            print(f"    roughly {100-o*100:.0f} per cent of the off-coil floor is something")
            print(f"    else -- regularity, the plasma edge, or the reference's own")
            print(f"    accuracy. Subtracting the coils removes the dominant global")
            print(f"    term and a minority of the off-coil one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
