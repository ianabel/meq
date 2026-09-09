"""One freegs4e DIVERTED reference -> a MEQ FREE-BOUNDARY run: the two profile
tables and the TOML.

WHY THIS IS NOT make_case.py.  That script writes the FIXED-boundary comparison:
freegs4e's converged LCFS is fitted to MXH and handed to MEQ as a boundary, so
psi is zero on Gamma and the profile tables are tabulated against
psi_MEQ = psi_fgs - psi_bndry.  This one writes the FREE-boundary problem --
Gamma is an artificial semicircle in the vacuum, the conductors are meshed in,
and psi_bnd is an unknown -- so the tables are against NORMALISED flux and the
coils and the plasma current are inputs.  examples/limited-tokamak.toml is the
same shape for a LIMITED machine; this is its diverted sibling.

THE TWO CONVENTIONS THAT MUST BE GOT RIGHT, AND BOTH ARE INVERTIBLE MISTAKES
THAT CONVERGE.

1.  THE NORMALISED FLUX RUNS THE OTHER WAY.  freegs4e's
    psi_n = ( psi - psi_axis )/( psi_bndry - psi_axis ) is 0 on the axis and 1
    on the boundary.  MEQ's Psi = ( psi - psi_bnd )/( psi_ax - psi_bnd ) is 1 on
    the axis and 0 on the boundary.  So Psi = 1 - psi_n, and a table written
    without the flip is a profile turned inside out -- hollow where it should be
    peaked -- which still solves.

2.  THE ABSCISSA MOVES SO THE DERIVATIVE SCALES.  freegs4e's saved pprime and
    ffprime are d/dpsi against the unnormalised flux; its own docstring says
    d/dpsi_n and is wrong, which tools/freegs4e-benchmark/README.md records at
    length.  meq::NormalisedMHDSource wants d/dPsi, so the values multiply by
    the span psi_ax - psi_bnd.  Both profiles carry the same factor, and
    [source] PlasmaCurrent makes the amplitude an unknown, so getting this
    wrong is INVISIBLE -- the scale absorbs it.  It is done correctly anyway,
    because the scale MEQ reports is then 1 by construction and is a check.

THE AMPLITUDES ARE READ RATHER THAN RE-DERIVED.  fgsref.py builds them with
split_amplitudes( frac_p, R0 ), and R0 there is the case's INPUT 1.35 rather
than the achieved 1.3156 -- a 5% difference in the p'/ff' ratio if re-derived
from the wrong one.  pshape( 0 ) is exactly 1, so P0 and F0 are the first
entries of the saved tabulations, and the assertion below says so.

USAGE

    python3 make_diverted_case.py <case.npz> <case.json> <outdir> <stem>
"""

import json
import os
import sys

import numpy as np


def pshape(x, alpha, beta):
    """freegs4e's own profile shape, against ITS psi_n."""
    return (1.0 - np.clip(x, 0.0, 1.0) ** alpha) ** beta


def shape_in_meq_flux(psi, alpha, beta):
    """The same shape and its d/dPsi, against MEQ's Psi = 1 - psi_n.

    UNCLIPPED BELOW Psi = 0, deliberately.  The table has to carry values in the
    vacuum, which the solve evaluates before ConfineToPlasma switches the source
    off, and freegs4e's clip would make it CONSTANT there -- a kink at the
    plasma edge, which is exactly where FB-4 says the rate is decided.  The
    analytic continuation is smooth and has the edge behaviour MEQ needs:
    ( 1 - u^a )^b with u = 1 - Psi vanishes like ( a Psi )^b as Psi -> 0, so
    j = b = 2 and ConfineToPlasma leaves the residual C^1.

    The grid stops at Psi = 1 -- the magnetic axis -- because u^a is not real
    for u < 0 at fractional a.  meq::SplineProfile clamps to the endpoint value
    outside its table, which is the right continuation there anyway.
    """
    u = 1.0 - np.asarray(psi, dtype=float)
    if np.any(u < 0.0):
        raise SystemExit("Psi > 1 is inside the magnetic axis and u^alpha is "
                         "not real there; stop the grid at 1")
    inner = 1.0 - np.power(u, alpha)
    value = np.power(inner, beta)
    # d/dPsi = ( d/du ) * ( du/dPsi ), and du/dPsi = -1
    derivative = beta * np.power(inner, beta - 1.0) * alpha * np.power(u,
                                                                      alpha - 1.0)
    return value, derivative


def write_profile(path, psi, value, derivative, what, note):
    with open(path, "w") as f:
        f.write("# %s\n#\n" % what)
        for line in note.splitlines():
            f.write("# %s\n" % line)
        f.write("#\n#   Psi            value                    d/dPsi\n")
        for p, v, d in zip(psi, value, derivative):
            f.write("  %+.6f   %+.15e   %+.15e\n" % (p, v, d))


def main():
    npz_path, json_path, outdir, stem = sys.argv[1:5]
    d = np.load(npz_path)
    j = json.load(open(json_path))

    psi_axis = float(d["psi_axis"])
    psi_bndry = float(d["psi_bndry"])
    span = psi_axis - psi_bndry          # MEQ's, and POSITIVE for these cases

    psi_n = np.asarray(d["psi_n"], dtype=float)
    if abs(psi_n[0]) > 1e-12:
        raise SystemExit("psi_n does not start at 0; P0 and F0 cannot be read "
                         "off the tabulations")
    # THE TABULATED ARRAYS AND NOT THE USED ONES, AND THE DIFFERENCE IS
    # REPORTED RATHER THAN SILENT.  fgsref.py fits a UnivariateSpline to its own
    # analytic shape before solving, and the fit MOVES it: measured on this
    # case, 2.514e-05 of the amplitude in p' and 1.707e-02 in ff'.  The
    # reference equilibrium is therefore the SMOOTHED profile's, not the
    # analytic one's.
    #
    # MEQ takes the analytic shape anyway, for two reasons that matter more here
    # than agreeing to the last digit: it is smooth and extends below Psi = 0
    # into the vacuum, which a tabulation on psi_n in [ 0, 1 ] cannot; and its
    # edge behaviour is exactly j = 2, which is what FB-4 says decides the rate
    # and what ConfineToPlasma needs.  So the reference's X-point is a
    # CROSS-CHECK at the per-cent level rather than a tight comparison, and the
    # caller is told the number.
    P0 = float(np.asarray(d["pprime_tabulated"])[0])
    F0 = float(np.asarray(d["ffprime_tabulated"])[0])

    shape = j["profile_shape"]
    pa, pb = _exponents(shape["pprime"])
    fa, fb = _exponents(shape["ffprime"])

    # THE TABULATION IS CHECKED AGAINST THE REFERENCE'S OWN before it is used,
    # on psi_n where both are defined.  This is the whole conversion in one
    # assertion: shape, amplitude and abscissa together.
    for name, amp, a, b in (("pprime", P0, pa, pb), ("ffprime", F0, fa, fb)):
        check = amp * pshape(psi_n, a, b)
        tab = np.asarray(d[name + "_tabulated"])
        rel = np.max(np.abs(check - tab))/max(abs(amp), 1e-300)
        if rel > 1e-10:
            raise SystemExit("the analytic %s shape does not reproduce the "
                             "saved tabulation: %.3e relative" % (name, rel))
        moved = np.max(np.abs(np.asarray(d[name]) - tab))/max(abs(amp), 1e-300)
        print("%-8s amplitude %+.15e, spline moved it %.3e of that"
              % (name, amp, moved))

    grid = np.linspace(-0.30, 1.00, 131)
    for name, amplitude, (a, b) in (("pprime", P0, (pa, pb)),
                                    ("ggprime", F0, (fa, fb))):
        value, derivative = shape_in_meq_flux(grid, a, b)
        path = os.path.join(outdir, "%s-%s.dat" % (stem, name))
        write_profile(
            path, grid, amplitude * span * value, amplitude * span * derivative,
            "d%s/dPsi( Psi ) for %s, from %s"
            % ("p" if name == "pprime" else "gg", j["name"],
               os.path.basename(npz_path)),
            "Psi = 1 - psi_n: freegs4e's normalised flux is 0 on the AXIS and\n"
            "1 on the boundary, and MEQ's is the other way up.\n"
            "The values are d/dpsi in the reference and d/dPsi here, so they\n"
            "carry a factor of the span %.15e.\n"
            "Amplitude %.15e, shape ( 1 - psi_n^%g )^%g."
            % (span, amplitude, a, b))
        print("wrote %s" % path)

    print("psi_axis  %.15e" % psi_axis)
    print("psi_bndry %.15e" % psi_bndry)
    print("span      %.15e" % span)
    print("Ip        %.15e" % float(d["Ip"]))
    for lab, cr, cz, cur in zip(d["coil_labels"], d["coil_R"], d["coil_Z"],
                                d["coil_currents"]):
        print("coil %-5s R %.6f Z %+.6f  I %+.15e" % (lab, cr, cz, cur))
    for x in np.asarray(d["xpoints"]):
        print("xpoint    R %.9f Z %+.9f  psi %.15e" % (x[0], x[1], x[2]))


def _exponents(expr):
    """'P0*(1-psi_n**1.5)**2' -> ( 1.5, 2.0 )."""
    head, tail = expr.split("**", 1)
    alpha = float(tail.split(")")[0])
    beta = float(expr.rsplit("**", 1)[1])
    return alpha, beta


if __name__ == "__main__":
    main()
