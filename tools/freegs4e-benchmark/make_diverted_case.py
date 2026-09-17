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
import math
import os
import subprocess
import sys

import numpy as np

import conductors
from mkcoldguess import design_footprint


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

    # A FRACTIONAL beta IS NOT REAL BELOW Psi = 0, AND FIVE OF THESE SEVEN
    # CASES HAVE ONE.  For Psi < 0 the vacuum is u > 1, so `inner` is NEGATIVE,
    # and `inner**1.5` is a NaN -- which arrives as a profile table MEQ refuses
    # to parse rather than as anything naming a profile exponent.
    #
    # THE CONTINUATION IS ODD IN `inner`: sign( x ) |x|^beta.  It is real
    # everywhere, it is C^1 for every beta > 1 (the derivative carries
    # |x|^(beta-1), which vanishes at x = 0), and for INTEGER beta below the
    # axis it is identical to the power -- so nothing about the cases that
    # already work moves.  Above Psi = 0, where `inner` is positive, the two
    # agree by definition and this is the analytic shape exactly.
    #
    # What it does NOT do is claim to be the physical profile out there. It is
    # a smooth, monotone, small extension into a region where ConfineToPlasma
    # switches the source off; what matters is that the table is finite and has
    # no kink at the plasma edge, which is where FB-4 says the rate is decided.
    magnitude = np.abs(inner)
    sign = np.sign(inner)
    even = abs(beta - round(beta)) < 1e-12 and int(round(beta)) % 2 == 0
    value = np.power(magnitude, beta)*(1.0 if even else sign)
    # d/dPsi = ( d/du ) * ( du/dPsi ), and du/dPsi = -1
    derivative = (beta*np.power(magnitude, beta - 1.0)*alpha
                  *np.power(u, alpha - 1.0))
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
    for c in conductors.from_npz(d):
        print("coil %-10s R %.6f Z %+.6f  %.4f x %.4f  I %+.15e  %s"
              % (c["label"], c["R"], c["Z"], c["half_width"], c["half_height"],
                 c["current"], c["kind"]))
    for x in np.asarray(d["xpoints"]):
        print("xpoint    R %.9f Z %+.9f  psi %.15e" % (x[0], x[1], x[2]))

    # THE GUESS FIRST, because the TOML needs its peak for [source] PsiAxis --
    # and the guess needs the disc radius, which is geometry and does not
    # depend on it.  Built by mkcoldguess.py and NOT by mkexactguess.py: the
    # latter sums Green's functions over the reference's CONVERGED Jtor, which
    # puts the answer in the starting position.
    disc = disc_radius(d)
    here = os.path.dirname(os.path.abspath(__file__))
    guess = [sys.executable, os.path.join(here, "mkcoldguess.py"), npz_path,
             os.path.join(outdir, "%s-guess.mesh" % stem),
             os.path.join(outdir, "%s-guess.gf" % stem), "%g" % disc, "32"]
    print("  " + " ".join(guess))
    done = subprocess.run(guess, capture_output=True, text=True)
    sys.stdout.write(done.stdout)
    sys.stderr.write(done.stderr)
    if done.returncode != 0:
        raise SystemExit("the guess could not be built")
    peak = [line for line in done.stdout.splitlines()
            if line.startswith("PsiAxis ")]
    if not peak:
        raise SystemExit("the guess did not report its peak")

    write_toml(os.path.join(outdir, "%s.toml" % stem), d, j, stem,
               float(peak[-1].split()[1]))


# ---------------------------------------------------------------------------
# THE TOML
# ---------------------------------------------------------------------------
#
# EVERY NUMBER BELOW COMES OUT OF THE REFERENCE, and that is the point of
# writing it rather than hand-copying it.  examples/diverted-tokamak.toml is
# the hand-written version of exactly this file for case A, and its header
# carries the generator command in a COMMENT -- four conductor rectangles
# written a second time, in a second convention, agreeing with the [[coils]]
# blocks because somebody kept them in step.  [mesh.generate] removed the need
# for the comment; this removes the need for the hand.
#
# WHAT IS CHOSEN HERE RATHER THAN READ, with the reason, because these are the
# only four judgements in the file:
#
#   Gamma          1.25 x the smallest circle containing every conductor. The
#                  exterior must be Delta*-harmonic, so a current outside Gamma
#                  is not admissible; the margin is what keeps the truncated
#                  Gegenbauer expansion converging on a trace that is not
#                  hugging a source.  Case A's hand-written 2.40 against a
#                  conductor radius of 1.914 is 1.254, so this reproduces it.
#   the disc       1.08 x Gamma.  D_h is CUT from the background mesh at Gamma,
#                  so the disc has to be strictly larger; 8% is about half an
#                  element at these sizes.  A's hand-written pair is 2.6 / 2.40.
#   the plasma box the reference's own LCFS bounding box UNION its X-points,
#                  grown by 0.3 of the minor radius.  The X-points matter: a
#                  box that stops at the separatrix leaves the private flux
#                  region coarse, and that is where ConfineToPlasma has to
#                  decide which elements carry current.
#   element sizes  7% of the disc in the background, 15% of the minor radius in
#                  the plasma box, and inside the conductors 80% of the
#                  SMALLEST half-extent in the machine -- which on TCV is
#                  T2/T3's capped 0.0166 m rather than the usual 0.05.
#
# THE X-POINT SEED IS THE REFERENCE'S ACTIVE NULL, unoffset.  examples/
# diverted-tokamak-xpoint.toml deliberately seeds 5 cm away to demonstrate that
# the border MOVES the point, and that demonstration belongs in one file rather
# than in seven: here the seed is the best available starting value and the
# question being asked is what the converged equilibrium is.
def _wrap(text, width):
    """Greedy wrap, so a machine's notes do not become a 200-column comment."""
    out, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = word if not line else line + " " + word
    if line:
        out.append(line)
    return out


def _round_up(x, step):
    """x rounded up to a multiple of `step`, so the files carry readable
    numbers rather than 2.3925000000000005."""
    return math.ceil(x/step - 1e-12)*step


def gamma_radius(d):
    return _round_up(1.25*conductors.bounding_radius(conductors.from_npz(d)),
                     0.05)


def disc_radius(d):
    return _round_up(1.08*gamma_radius(d), 0.05)


def write_toml(path, d, j, stem, psi_axis_guess):
    flat = conductors.from_npz(d)
    rho_conductors = conductors.bounding_radius(flat)
    gamma = gamma_radius(d)
    disc = disc_radius(d)

    box = design_footprint(d)
    minor = 0.5*min(box[1] - box[0], box[3] - box[2])
    xpts = np.atleast_2d(np.asarray(d["xpoints"], dtype=float))

    smallest = min(min(c["half_width"], c["half_height"]) for c in flat)

    # THE ACTIVE NULL IS THE ONE AT psi_bndry, and on a diverted plasma that is
    # what makes it active.  Several of these machines are up-down SYMMETRIC
    # double nulls, where the two carry the same flux to rounding and either is
    # a correct answer -- so the seed selects which saddle MEQ's border
    # follows, exactly as [source] PsiAxis selects one O-point.
    psi_bndry = float(d["psi_bndry"])
    # BY INDEX, not by identity: iterating a 2-D array yields a fresh view each
    # time, so `row is active` is false for the very row min() returned.
    active_at = int(np.argmin(np.abs(xpts[:, 2] - psi_bndry)))
    active = xpts[active_at]

    # THE SEED IS THE TARGET, NOT THE ACHIEVED NULL.  fgsref.py's case table
    # names the X-point locations the control system was told to hit, and those
    # are an INPUT -- what the operator asked the shape to be.  Where the null
    # ended up is an OUTPUT, and seeding MEQ with it would make this a
    # measurement of how well MEQ polishes an answer it was handed.
    #
    # WHICH target, when a machine names two: not "the one nearest the achieved
    # null", which is reading the answer, but the LOWER one -- the convention a
    # single-null scenario means, and the one every asymmetric case here ends
    # up on anyway.
    targets = np.asarray(d["design_xpoints"], dtype=float).reshape(-1, 2)
    if not targets.size:
        raise SystemExit("the reference names no target X-points, so there is "
                         "no legitimate seed for [boundary.xpoint]")
    seed = targets[int(np.argmin(targets[:, 1]))]

    lines = []
    w = lines.append
    w("# %s -- A FREE-BOUNDARY MEQ RUN OF" % stem.upper())
    w("# freegs4e's %s." % j["name"])
    w("#")
    w("# WRITTEN BY tools/freegs4e-benchmark/make_diverted_case.py FROM THE")
    w("# REFERENCE ITSELF -- every number here, the conductor rectangles and")
    w("# the mesh included, is read out of the .npz and .json that script was")
    w("# given.  Do not hand-edit it: re-run tools/freegs4e-benchmark/")
    w("# make_all.sh.")
    w("#")
    for line in _wrap(" ".join(j["notes"].split()), 72):
        w("# %s" % line)
    w("#")
    w("# THE REFERENCE'S OWN ANSWER, for anything comparing against it:")
    w("#")
    w("#   psi_axis   %.15e   axis at ( %.6f, %+.6f )"
      % (float(d["psi_axis"]), float(j["Raxis"]), float(j["Zaxis"])))
    w("#   psi_bndry  %.15e" % psi_bndry)
    w("#   I_p        %.15e A" % float(d["Ip"]))
    w("#   R0 %.4f  a %.4f  A %.3f  kappa %.4f  delta %+.4f"
      % (float(j["R0"]), minor, float(j["aspect_ratio"]), float(j["kappa"]),
         float(j["delta"])))
    # A MACHINE CAN HAVE MANY SADDLES AND MOST OF THEM ARE NOT THE PLASMA'S.
    # MAST-U's critical-point search finds fourteen, out at the coil nulls; the
    # ones worth writing down are those near the plasma edge, and the ACTIVE
    # one is the plasma edge by definition.
    order = np.argsort(np.abs(xpts[:, 2] - psi_bndry))
    for at in order[:4]:
        x = xpts[at]
        w("#   X-point    ( %.9f, %+.9f ) psi %.15e%s"
          % (x[0], x[1], x[2], "  ACTIVE" if at == active_at else ""))
    if len(order) > 4:
        w("#   ( and %d more saddles further out, at the coil nulls )"
          % (len(order) - 4))
    w("#")
    w("# ONE DELIBERATE DEVIATION FROM THE REFERENCE, and it is the same one")
    w("# examples/diverted-tokamak.toml records: fgsref.py fits a")
    w("# UnivariateSpline to its own analytic profile shape before solving and")
    w("# the fit MOVES it, by %.3e of the amplitude in p' and %.3e in ff'."
      % (float(j["spline_smoothing_pprime"]),
         float(j["spline_smoothing_ffprime"])))
    w("# The tables here are the ANALYTIC shape, because it extends below")
    w("# Psi = 0 into the vacuum -- which a tabulation on psi_n in [ 0, 1 ]")
    w("# cannot -- and because its edge behaviour is exactly j = 2, which FB-4")
    w("# says decides the rate and which ConfineToPlasma needs.  So this is a")
    w("# comparison at the per-cent level and not to the last digit.")
    w("#")
    w("# RUN IT WITH meq-run, WHICH MAKES THE MESH:")
    w("#")
    w("#     ./build/meq-run examples/%s.toml" % stem)
    w("")
    w("[mesh]")
    w('File = "examples/%s.msh"' % stem)
    w("")
    w("# THE DISC IS NOT GAMMA.  D_h is cut from this mesh at")
    w("# [boundary.exterior] Radius = %g, which has to fit strictly inside." % gamma)
    w("# Every conductor lies within %.4f m of the origin." % rho_conductors)
    w("[mesh.generate]")
    w('Tool = "halfdisc"')
    w("Radius = %g" % disc)
    w("Size = %g" % round(0.07*disc, 4))
    w("PlasmaRMin = %.4f" % box[0])
    w("PlasmaRMax = %.4f" % box[1])
    w("PlasmaZMin = %.4f" % box[2])
    w("PlasmaZMax = %.4f" % box[3])
    w("PlasmaSize = %.4f" % round(0.15*minor, 4))
    # A FIXED SIZE AND NOT THE SMALLEST CONDUCTOR'S.  Deriving it from the
    # thinnest rectangle lets one conductor set the whole mesh's cost: MAST-U's
    # solenoid is 0.02 m wide, so four fifths of that is 0.016 m over 3.16 m of
    # winding, and the graded transition out of it reaches most of the machine.
    # A thin conductor gets thin elements from its own geometry anyway; this is
    # only the EXTRA refinement asked for inside one.
    w("# The smallest conductor here is %.4f m in its shorter half-extent."
      % smallest)
    w("CoilSize = 0.04")
    w("")
    w("[discretisation]")
    w("PolynomialDegree = 2")
    w("")
    w("[solver]")
    w("NewtonMaxIterations = 200")
    w("NewtonRelativeTolerance = 1.0e-10")
    w("# The plasma support is a SET, so its derivative is a surface term the")
    w("# Jacobian does not have: it is frozen within each solve and re-decided")
    w("# between them.  MEASUREMENTS.md M-82 -- without this the first solve")
    w("# stalls at the iteration cap rather than running slowly.")
    w("PlasmaSupportSweeps = 4")
    w("")
    w("[source]")
    w('Type = "mhd"')
    w("Normalised = true")
    w("ConfineToPlasma = true")
    w("# THE STARTING VALUE OF AN UNKNOWN, taken from the peak of the cold")
    w("# guess below.  psi_ax is a functional of the solution -- a bordered")
    w("# Newton row and not data -- so what a file supplies is where to START,")
    w("# and the best estimate available before solving is the best available")
    w("# guess's own peak.  The reference converged to %.6e, which is %.0f%%"
      % (float(d["psi_axis"]),
         100.0*abs(psi_axis_guess - float(d["psi_axis"]))
         /max(abs(float(d["psi_axis"])), 1e-300)))
    w("# away: this is a basin and not an answer.")
    w("PsiAxis = %.9e" % psi_axis_guess)
    w("PlasmaCurrent = %.10e" % float(d["Ip"]))
    w("Mu0 = 1.25663706212e-6")
    w('PPrimeFile = "examples/%s-pprime.dat"' % stem)
    w('GGPrimeFile = "examples/%s-ggprime.dat"' % stem)
    w("")
    w("# THE CONDUCTORS, WRITTEN ONCE.  [mesh.generate] derives the mesher's")
    w("# --coil rectangles from these blocks, in this order, so the mesh is")
    w("# aligned to the conductors the solve integrates over.  Currents are the")
    w("# TOTAL through each cross-section: freegs4e's turns, circuit")
    w("# multipliers and solenoid windings are all resolved into them by")
    w("# tools/freegs4e-benchmark/conductors.py.")
    kinds = sorted(set(c["kind"] for c in flat))
    w("# This machine's are %s." % ", ".join(kinds))
    for c in flat:
        w("")
        w("[[coils]]")
        w('Name = "%s"' % c["label"])
        w("CentreR = %.6f" % c["R"])
        w("CentreZ = %.6f" % c["Z"])
        w("HalfWidth = %.6f" % c["half_width"])
        w("HalfHeight = %.6f" % c["half_height"])
        w("Current = %.10e" % c["current"])
    w("")
    w("[boundary]")
    w('Type = "zero"')
    w("")
    w("# THE TARGET NULL, as an INITIAL VALUE of two unknowns rather than a")
    w("# prescription: XP-3's border closes q_r = q_z = 0 and")
    w("# psi_bnd = psi_h( r_X, z_X ) on the same Newton, and the seed selects")
    w("# WHICH saddle is followed.")
    w("#")
    w("# THIS IS THE SHAPE THE OPERATOR ASKED FOR and not the null the")
    w("# reference achieved, which is %.4f m away.  That distinction is the"
      % float(np.hypot(seed[0] - active[0], seed[1] - active[1])))
    w("# whole of what makes this a cold start.")
    w("[boundary.xpoint]")
    w("R = %.12f" % seed[0])
    w("Z = %.12f" % seed[1])
    w("")
    w("[boundary.exterior]")
    w("Radius = %g" % gamma)
    w("CentreZ = 0.0")
    w("Modes = 10")
    w("")
    w("# NOT A NICETY.  A cold start on a free-boundary machine case wanders")
    w("# and does not converge; the guess is part of the problem statement.")
    w("#")
    w("# IT IS COLD.  tools/freegs4e-benchmark/mkcoldguess.py builds it from")
    w("# the MACHINE -- the conductors and their currents, the target plasma")
    w("# current, and the shape the operator asked for -- and from nothing")
    w("# that was solved for.  mkexactguess.py is the OTHER one, which sums")
    w("# Green's functions over the reference's own converged Jtor; this")
    w("# comment named it for a while and was wrong, which matters: a run")
    w("# seeded from mkexactguess has the answer in the starting position,")
    w("# and is a different benchmark.")
    w("[initialguess]")
    w('Type = "gridfunction"')
    w('File = "examples/%s-guess.gf"' % stem)
    w('MeshFile = "examples/%s-guess.mesh"' % stem)
    w("")
    w("[output]")
    w('Prefix = "%s"' % stem)
    w("GridNR = 129")
    w("GridNZ = 129")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote %s" % path)
    print("  disc %g, Gamma %g, conductors within %.4f, %d coils"
          % (disc, gamma, rho_conductors, len(flat)))
    return disc


def _exponents(expr):
    """'P0*(1-psi_n**1.5)**2' -> ( 1.5, 2.0 )."""
    head, tail = expr.split("**", 1)
    alpha = float(tail.split(")")[0])
    beta = float(expr.rsplit("**", 1)[1])
    return alpha, beta


if __name__ == "__main__":
    main()
