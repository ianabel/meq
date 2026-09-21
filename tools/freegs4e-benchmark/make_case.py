"""One freegs4e reference -> a MEQ FIXED-BOUNDARY run: TOML plus two profile tables.

THE RECIPE, WHICH IS IAN'S AND IS WHAT MAKES THIS A VALID PROBLEM RATHER THAN AN
APPROXIMATION OF ONE.  Solve a machine FREE boundary in freegs4e at whatever
resolution can be afforded; take its psi_n = 0.95 surface; fit MXH to it; and
hand MEQ that curve as Gamma with psi = 0 on it.  psi is defined only up to an
additive constant, so `psi_MEQ = psi_freegs4e - psi( the 0.95 surface )` is the
SAME equilibrium in a different gauge -- and inside that surface it satisfies
the fixed-boundary Grad-Shafranov equation exactly, with the machine's own
p' and ff' as the source.  Nothing is thrown away but the region between the
0.95 surface and the wall, which is not part of a fixed-boundary problem.

WHY 0.95 AND NOT THE SEPARATRIX.  Every reference here is DIVERTED, so its last
closed surface passes through an X-point and has a CORNER.  MXH is a truncated
Fourier series in tR and cannot represent one -- fitting the separatrix reads
2.6e-03 to 1.3e-02 m against 1.9e-05 on a smooth shape -- and the source there
is at the plasma edge, where `ConstrainPaxisIp`'s alpha_n = 1.2 makes a
fractional power and caps any convergence rate at about 1.2 whatever the
discretisation does.  At 0.95 the surface is smooth and closed, and the profiles
are analytic across it with p' and ff' both comfortably nonzero.  **That is the
whole reason these cases exist**: they are the only machine-shaped problems in
this tree on which a high-order claim can be made at all.  surface.py has the
same argument from the fitting side, and MEQ's own ExtensionConvergence makes
the same choice for the same reason.

WHAT A CASE WRITTEN BY THIS SCRIPT IS AND IS NOT.

  * IT IS self-contained.  Gamma is the MXH coefficients in the TOML -- an
    analytic, smooth, exactly-known curve -- so a convergence study on it is
    valid at any order however the coefficients were arrived at.  The fit
    residual bounds the comparison against freegs4e and bounds NOTHING about
    MEQ's own rates.
  * IT IS NOT a place to measure coil subtraction.  Gamma is the plasma edge,
    so every conductor is OUTSIDE the computational domain and enters only
    through a Dirichlet datum that is identically zero.  `[conductors] Model`
    has nothing to act on here, and that is structural rather than an omission
    -- see the TODO entry these cases were written for.
"""
import json
import os
import numpy as np

from scipy.interpolate import CubicSpline

from mxh import fit_mxh, worst_distance
from surface import interior_surface

from convert import (to_meq_normalised_table, write_meq_profile,
                     NORMALISED_ABSCISSA)

# THE BOX IS CUT, SO ITS ELEMENTS ARE NOT THE PROBLEM'S ELEMENTS.
#
# Omega_h is the union of background elements lying inside Gamma, and the box
# below pads the fitted surface by 35% of a and of kappa*a in each direction --
# so the box is ( 2.7 a ) x ( 2.7 kappa a ) and the plasma is pi kappa a^2.  The
# ratio is pi/7.29 = 0.431 for an ELLIPSE, and a shaped surface is not one --
# MEQ reported 3427 of 9176 on machine A, which is 0.373.  It is a starting
# point for picking NR and NZ and not a claim: every shipped header quotes
# the count MEQ actually reported.
KEPT_FRACTION = 0.373

# The padding itself, as a fraction of a and of kappa*a.  It has to be enough
# that Gamma is strictly inside the box with room for the transfer paths, and
# every element outside Gamma is thrown away, so more is waste.
PAD = 0.35


def choose_mesh(target_elements, kappa, base_low=600.0, base_high=2600.0):
    """( NR, NZ, RefinementLevels ) for about `target_elements` in the PLASMA.

    Cells are kept near square -- the box is 2.7a by 2.7 kappa a, so NZ/NR is
    kappa -- because the extension technique's transfer paths are measured
    against the local element size and a strongly anisotropic background makes
    `dist( Gamma_h, Gamma )/h_loc` a different number on the two faces of one
    element.

    AND THE BACKGROUND IS HELD MODEST SO THAT RefinementLevels IS THE SIZE KNOB.
    A target is reached by refining a few-thousand-element background rather
    than by a huge NR x NZ, and the reason is that the two are NOT equivalent:
    `mfem::Mesh::UniformRefinement` on triangles is exact bisection, so every
    refinement of one background is NESTED in the last, while changing NR x NZ
    changes WHICH elements fall inside Gamma and therefore changes the geometry.
    Sweeping RefinementLevels on a shipped case is a refinement study on nested
    meshes; sweeping NR is not.
    """
    best = None
    for refine in range(0, 6):
        base = target_elements / 4.0 ** refine
        if base < base_low or base > base_high:
            continue
        cells = base / (KEPT_FRACTION * 2.0)
        nr = max(2, int(round(np.sqrt(cells / kappa))))
        nz = max(2, int(round(nr * kappa)))
        got = 2 * nr * nz * 4 ** refine * KEPT_FRACTION
        score = abs(np.log(got / target_elements))
        if best is None or score < best[0]:
            best = (score, nr, nz, refine, got)

    if best is None:
        raise ValueError(
            f"no ( background, RefinementLevels ) reaches {target_elements} "
            f"elements with a background in [ {base_low:.0f}, {base_high:.0f} ]")

    _, nr, nz, refine, got = best
    return nr, nz, refine, got


def build(npz_path, outdir, stem=None, degree=2, target_elements=4000,
          harmonics=20, grid=129, level=0.95, title=None, shape_note=None,
          flux_surfaces=False):
    d = np.load(npz_path, allow_pickle=True)
    reference = os.path.splitext(os.path.basename(npz_path))[0]
    stem = stem or reference
    os.makedirs(outdir, exist_ok=True)

    psi_axis = float(d["psi_axis"])
    psi_bndry = float(d["psi_bndry"])
    D = psi_bndry - psi_axis
    grid_n = len(np.asarray(d["R"], float))

    sR, sZ, psi_surface = interior_surface(
        d["R"], d["Z"], d["psi"], psi_axis, psi_bndry,
        float(d["Raxis"]), float(d["Zaxis"]), level)

    fit = fit_mxh(sR, sZ, n_harmonics=harmonics)
    shape_err = worst_distance(np.asarray(sR, float), np.asarray(sZ, float), fit)

    note = (f"Converted from {os.path.basename(npz_path)}, a {grid_n} x {grid_n}\n"
            f"freegs4e free-boundary solve.\n"
            f"psi_axis = {psi_axis:.9e}, psi_surface = {psi_surface:.9e},\n"
            f"the psi_n = {level} surface, chosen to avoid the X-point corner\n"
            f"and the fractional power at the plasma edge.\n"
            f"BOTH abscissae are normalised flux and they run OPPOSITE WAYS:\n"
            f"freegs4e's psi_n is 0 on the axis and 1 on the separatrix, and\n"
            f"MEQ's Psi is 1 on the axis and 0 on Gamma, so\n"
            f"Psi = 1 - psi_n/{level}.\n"
            f"freegs4e's pprime/ffprime are d/dpsi -- its docstring says\n"
            f"d/dpsi_n and is wrong -- so the value column carries ONE factor\n"
            f"of dpsi/dPsi = psi_ax and nothing else. Do not pre-divide by\n"
            f"psi_ax: meq::NormalisedMHDSource does that itself.")

    # psi_ax IN MEQ'S GAUGE, WHICH IS THE CHAIN FACTOR AND THE [source] SEED.
    # MEQ's psi is freegs4e's shifted so that Gamma is zero, so its axis flux is
    # what the reference's axis sits above that surface.
    psi_ax = psi_axis - psi_surface

    pp = to_meq_normalised_table(d["psi_n"], d["pprime"], level)
    gg = to_meq_normalised_table(d["psi_n"], d["ffprime"], level)
    # d/dPsi rather than d/dpsi: one factor of dpsi/dPsi = psi_ax on BOTH
    # columns of each table, which is the whole of the conversion's arithmetic.
    pp = (pp[0], pp[1] * psi_ax, pp[2] * psi_ax)
    gg = (gg[0], gg[1] * psi_ax, gg[2] * psi_ax)

    pp_name = f"{stem}-pprime.dat"
    gg_name = f"{stem}-ggprime.dat"
    write_meq_profile(os.path.join(outdir, pp_name), *pp,
                      "dp/dPsi [Pa per unit normalised flux]", note,
                      NORMALISED_ABSCISSA)
    write_meq_profile(os.path.join(outdir, gg_name), *gg,
                      "g dg/dPsi [T^2 m^2 per unit normalised flux]", note,
                      NORMALISED_ABSCISSA)

    # The two profiles AT GAMMA, which is what says this is not the plasma edge.
    # Quoted as d/dpsi -- the physical quantity -- rather than as the table's
    # own d/dPsi, so that the number is comparable with any other profile.
    ppEdge = float(np.interp(0.0, pp[0], pp[1])) / psi_ax
    ggEdge = float(np.interp(0.0, gg[0], gg[1])) / psi_ax
    ppAxis = float(np.interp(1.0, pp[0], pp[1])) / psi_ax

    # ---- p AND g ON GAMMA, WHICH GRAD-SHAFRANOV DOES NOT DETERMINE ----------
    #
    # **THIS IS THE ONE THING A `fixed-*.toml` DOES NOT SAY ABOUT ITS OWN
    # EQUILIBRIUM, AND IT IS NOT AN OVERSIGHT -- IT IS WHAT THE EQUATION IS.**
    # Grad-Shafranov contains `dp/dpsi` and `g dg/dpsi` and nothing else, so a
    # case file fixes `p` and `g^2` only up to an additive constant each, and
    # MEQ's answer for `psi` is independent of both.  A FULL MHD equilibrium is
    # not: `B_phi = g/r` and the toroidal flux are proportional to `g`, and the
    # pressure is a profile rather than a gradient.
    #
    # So any code that solves more than Grad-Shafranov needs these two numbers
    # and MEQ's input cannot supply them.  Measured consequence, on DESC:
    # `tools/desc-benchmark` had to read them out of the freegs4e reference's
    # own `pressure` and `fpol` arrays, which put a third code into the posing
    # of a two-code comparison -- so a benchmark that was meant to be MEQ
    # against DESC was partly MEQ against freegs4e.  Two numbers here remove
    # that.
    #
    # THE INTERPOLATION IS convert.py'S OWN, at the same `level` and by the
    # same cubic spline that puts the profile tables' endpoint exactly on
    # Gamma, so these are consistent with the tables beside them rather than
    # merely close to them.
    pEdge = float(CubicSpline(d["psi_n"], d["pressure"])(level))
    gEdge = float(CubicSpline(d["psi_n"], d["fpol"])(level))

    R0, Z0, a, kappa = fit["R0"], fit["Z0"], fit["a"], fit["kappa"]
    nr, nz, refine, predicted = choose_mesh(target_elements, kappa)

    pad_r, pad_z = PAD * a, PAD * kappa * a
    box = dict(rmin=max(1e-3, R0 - a - pad_r), rmax=R0 + a + pad_r,
               zmin=Z0 - kappa * a - pad_z, zmax=Z0 + kappa * a + pad_z)

    def wrap(values, per_line=4):
        """The coefficient list across several lines, which TOML allows and a
        400-character line does not deserve."""
        items = [f"{v: .10e}" for v in values]
        rows = [", ".join(items[i:i + per_line])
                for i in range(0, len(items), per_line)]
        return "[\n\t" + ",\n\t".join(rows) + "\n]"

    cos_list = wrap(fit["cos"])
    sin_list = wrap(fit["sin"])

    # THE UP-DOWN ASYMMETRY IS THE COS COEFFICIENTS AND THE SHAPING IS THE SIN
    # ONES, WHICH IS THE OPPOSITE WAY ROUND FROM THE OBVIOUS GUESS.
    #
    # Z( t ) is odd in t whatever the harmonics do, so up-down symmetry is
    # R( -t ) = R( t ), which needs cos( tR( -t ) ) = cos( tR( t ) ). Writing
    # tR out, that holds exactly when c_0 and every c_n vanish -- the s_n drop
    # out, being odd. So s_1 is the TRIANGULARITY of an up-down symmetric
    # surface (arcsin delta) and c_0 is the TILT. mxh.py names c_0 the tilt for
    # the same reason.
    asymmetry = float(abs(np.asarray(fit["cos"], float)).max())
    triangularity = float(np.sin(fit["sin"][0])) if len(fit["sin"]) else 0.0
    shaping = float(abs(np.asarray(fit["sin"][1:], float)).max()) \
        if len(fit["sin"]) > 1 else 0.0

    header = title or f"A MACHINE-SHAPED FIXED-BOUNDARY EQUILIBRIUM: {stem}"

    toml = f"""# {header}
#
#     meq examples/{stem}.toml
#
# GENERATED -- edit tools/freegs4e-benchmark/make_case.py, not this file, and
# re-run tools/freegs4e-benchmark/make_fixed.sh.
#
# WHAT THIS IS. A real machine's equilibrium, posed as a FIXED-boundary problem
# on a smooth closed flux surface. The reference
#
#     {reference}
#
# was solved FREE boundary by freegs4e on a {grid_n} x {grid_n} grid -- a
# different code, a different algorithm, its own conductors and its own inverse
# solve for their currents -- and its psi_n = {level} surface was fitted to MXH.
# psi is defined only up to an additive constant, so subtracting that surface's
# flux gives the SAME equilibrium with psi = 0 on Gamma. Inside it, this is the
# fixed-boundary Grad-Shafranov equation with the machine's own p' and gg',
# exactly.
#
# THE REFERENCE'S OWN ANSWER, for anything comparing against it:
#
#     psi_axis   {psi_axis:.12e}   at ( {float(d['Raxis']):.6f}, {float(d['Zaxis']):.6f} )
#     psi_bndry  {psi_bndry:.12e}   the separatrix
#     psi at the psi_n = {level} surface   {psi_surface:.12e}
#     I_p        {float(d['Ip']):.6e} A
#
# so in MEQ's gauge -- psi = 0 on Gamma -- the axis flux is
#
#     psi_ax     {psi_ax:.12e} Wb/rad
#
# and [source] PsiAxis below is that number, as the starting value of an
# unknown rather than as an input. MEASUREMENTS.md M-147 is what MEQ reports
# on this file and how far it sits from the reference.
#
# WHY psi_n = {level} AND NOT THE SEPARATRIX, which is the one decision in this
# file that is not mechanical. The reference is diverted, so its last closed
# surface passes through an X-point and has a CORNER; MXH is a truncated
# Fourier series and cannot turn one. Worse, the source AT the separatrix is
# the plasma
# edge, where the reference's alpha_n = 1.2 makes p' a fractional power of
# 1 - Psi and caps any convergence rate at about 1.2 whatever the polynomial
# degree. At {level} the boundary is smooth and closed and the profiles are
# analytic across it:
#
#     p'  at Gamma  {ppEdge: .6e} Pa per Wb/rad
#                   which is {ppEdge / ppAxis if ppAxis else float('nan'):.4f} of its value on the axis
#     gg' at Gamma  {ggEdge: .6e} T^2 m^2 per Wb/rad
#
# Both comfortably nonzero, so this is NOT the trivial branch of CLAUDE.md's
# Traps section -- and see [initialguess] for the OTHER branch, which is.
#
# THE SHAPE, read out of the fit rather than asserted:
#
#     R0 {R0:.4f} m   a {a:.4f} m   aspect {R0 / a:.3f}   elongation {kappa:.4f}
#     triangularity sin( s_1 )                    {triangularity:+.4f}
#     largest further shaping | s_n |, n >= 2     {shaping:.4f}
#     up-down asymmetry, max | c_n | over n >= 0  {asymmetry:.4f}
#
# The last one is zero to machine precision for an exactly up-down symmetric
# machine, and it is the COS coefficients rather than the sin ones because
# Z( t ) is odd whatever they do -- see make_case.py where these are computed.
#
# MXH fit residual {shape_err:.3e} m on a minor radius of {a:.4f} m, i.e.
# {shape_err / a:.2e} relative, at {harmonics} harmonics. THAT NUMBER BOUNDS THE
# COMPARISON AGAINST freegs4e AND BOUNDS NOTHING ABOUT MEQ. Gamma here IS the
# MXH curve below -- analytic, smooth and exactly known -- so a convergence
# study on this file is valid at any order however the coefficients were
# arrived at.
#
# WHAT THIS FILE DOES NOT CARRY, said plainly because their absence is a
# property of the problem rather than an omission:
#
#   * NO CONDUCTORS. Gamma is the plasma edge, so every coil of the real
#     machine is OUTSIDE the computational domain and its whole influence is
#     in the Dirichlet datum -- which is zero. `[conductors] Model` would have
#     nothing to act on. A fixed-boundary case that CAN measure coil
#     subtraction needs a domain larger than the plasma;
#     examples/coils-rectangle.toml is the small one that does.
#   * NO EXTERIOR COUPLING, NO X-POINT, NO LIMITER AND NO PLASMA-CURRENT
#     BORDER. This is one Dirichlet problem, solved once, which is exactly why
#     it is cheap enough to be worth having beside the free-boundary set.

[mesh]
# The background the curved boundary is cut FROM. Omega_h is the union of these
# elements lying inside Gamma, so about {KEPT_FRACTION:.0%} of them survive and
# {predicted:.0f} elements are predicted here; MEASUREMENTS.md M-147 has what MEQ
# actually reported. Cells are near square, NZ/NR following the elongation.
#
# RefinementLevels IS THE SIZE KNOB AND THE REFINEMENT STUDY. Bisection of
# triangles is exact, so raising it gives a mesh NESTED in this one and the
# geometry does not move; raising NR or NZ instead changes which elements fall
# inside Gamma, which is a different problem rather than a finer one.
RMin = {box['rmin']:.6f}
RMax = {box['rmax']:.6f}
ZMin = {box['zmin']:.6f}
ZMax = {box['zmax']:.6f}
NR = {nr}
NZ = {nz}
RefinementLevels = {refine}

[discretisation]
PolynomialDegree = {degree}

[boundary]
# psi = 0 on Gamma, which is the MXH curve below and NOT the mesh's rectangle.
Type = "zero"

[boundary.shape]
# GS-2 section 3's extension technique: Omega_h is inscribed in Gamma, the
# datum is carried outward along transfer paths, and the band between them is
# continued with the SOLVED flux. This is what MEQ is for, and this curve is
# the only thing in this file that is not a number out of the reference.
Type = "mxh"
R0 = {R0:.10f}
Z0 = {Z0:.10f}
MinorRadius = {a:.10f}
Elongation = {kappa:.10f}
CosCoefficients = {cos_list}
SinCoefficients = {sin_list}

# THE RAMP IS PART OF THE PROBLEM STATEMENT, NOT A TUNING PARAMETER.
#
# These sources do not vanish at psi = 0 -- see the two numbers at the top of
# this file -- so this is not the trivial branch of CLAUDE.md's Traps section.
# It is the SECOND SOLUTION beside it, which that file records for the
# high-beta source: "this equation has a small positive solution and a large
# one ... Newton from the Dirichlet datum walks straight onto it". The
# amplitude is the axis height the reference itself has, so the iterate starts
# at the right SCALE and Psi starts inside its own table.
[initialguess]
Type = "ramp"
Amplitude = {psi_ax:.10e}

[source]
Type = "mhd"
# dp/dPsi and g dg/dPsi against MEQ's NORMALISED flux -- one on the axis, zero
# on Gamma -- and both already the derivative quantities. See the tables' own
# headers, and tools/freegs4e-benchmark/convert.py, which records both the
# docstring error in freegs4e that this conversion survives and the reason the
# normalised form is the one that can be shipped.
PPrimeFile = "examples/{pp_name}"
GGPrimeFile = "examples/{gg_name}"

# psi_ax IS AN UNKNOWN OF THE SOLVE AND THIS IS ITS STARTING VALUE.
#
# Normalised = true closes the system with psi_ax - max psi = 0 in the border,
# so Psi is confined to [ 0, 1 ] BY CONSTRUCTION and the profiles are never
# evaluated off their own tables. That is not a refinement -- it is what makes
# this file shippable. Against psi in Wb/rad the table has to stop at the axis,
# meq::SplineProfile extends it by a constant above that, and the plateau
# carries a solution of its own: sweeping only the guess amplitude on machine A
# reached 0.0105, 1.2604 and 0.9996 of the reference's axis height, in that
# order, as the ramp grew. convert.py has the table.
Normalised = true
PsiAxis = {psi_ax:.10e}

[solver]
NewtonMaxIterations = 40
NewtonRelativeTolerance = 1.0e-10
NewtonAbsoluteTolerance = 1.0e-14

[output]
Directory = "."
Prefix = "{stem}"
GridNR = {grid}
GridNZ = {grid}
"""
    if flux_surfaces:
        toml += """FluxSurfaces = true
FluxSurfaceCount = 33
FluxAngleCount = 128
"""

    toml_path = os.path.join(outdir, f"{stem}.toml")
    with open(toml_path, "w") as f:
        f.write(toml)

    meta = dict(stem=stem, reference=reference, reference_grid=grid_n,
                psi_axis=psi_axis, psi_bndry=psi_bndry, D=D,
                psi_surface=psi_surface, level=level,
                shape_error_m=shape_err, minor_radius=a, R0=R0, Z0=Z0,
                kappa=kappa, aspect=R0 / a, shaping=shaping,
                asymmetry=asymmetry, triangularity=triangularity, box=box,
                nr=nr, nz=nz, refine=refine, degree=degree,
                predicted_elements=predicted,
                pprime_at_gamma=ppEdge, ggprime_at_gamma=ggEdge,
                # The two Grad-Shafranov does not fix; see where they are
                # computed. Pa and T m.
                p_at_gamma=pEdge, g_at_gamma=gEdge,
                # EVERY QUANTITY'S UNITS, BECAUSE TWO OF THEM ARE A TRAP AND
                # JSON HAS NOWHERE ELSE TO SAY SO.
                #
                # `pprime_at_gamma` is d/dpsi -- the PHYSICAL derivative, per
                # Wb/rad -- and the value column of `<stem>-pprime.dat` beside
                # it is d/dPsi, per unit NORMALISED flux.  They differ by
                # psi_ax, which is 0.0635 on fixed-h-circular, so the two
                # numbers a reader meets next to each other are a factor of
                # 15.7 apart and both are right.  Each file documents its own
                # units correctly and that was not enough: an agent converting
                # between them read a CORRECT conversion as 15.7x wrong before
                # spotting it, which is exactly the failure a units field
                # prevents and a units convention does not.
                units=dict(
                    psi_axis="Wb/rad", psi_bndry="Wb/rad",
                    psi_surface="Wb/rad", D="Wb/rad",
                    pprime_at_gamma="Pa per Wb/rad -- d/dpsi, NOT the "
                                    "-pprime.dat value column's d/dPsi; they "
                                    "differ by psi_ax",
                    ggprime_at_gamma="T^2 m^2 per Wb/rad -- d/dpsi, and the "
                                     "same caution as pprime_at_gamma",
                    p_at_gamma="Pa", g_at_gamma="T m",
                    minor_radius="m", R0="m", Z0="m", shape_error_m="m",
                    kappa="dimensionless", level="dimensionless"),
                cos=list(map(float, fit["cos"])),
                sin=list(map(float, fit["sin"])))
    with open(os.path.join(outdir, f"{stem}-meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    return toml_path, meta


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("npz", nargs="+")
    parser.add_argument("-o", "--outdir", default=None)
    parser.add_argument("--stem", default=None)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--elements", type=int, default=4000)
    parser.add_argument("--level", type=float, default=0.95)
    parser.add_argument("--harmonics", type=int, default=20)
    parser.add_argument("--flux-surfaces", action="store_true")
    args = parser.parse_args()

    for path in args.npz:
        out = args.outdir or os.path.dirname(path) or "."
        t, m = build(path, out, stem=args.stem, degree=args.degree,
                     target_elements=args.elements, level=args.level,
                     harmonics=args.harmonics,
                     flux_surfaces=args.flux_surfaces)
        print(f"{m['stem']:26s} fit {m['shape_error_m']:.3e} m "
              f"(a = {m['minor_radius']:.4f}, kappa = {m['kappa']:.3f})  "
              f"{m['nr']}x{m['nz']} r{m['refine']} ~ {m['predicted_elements']:.0f} el"
              f"  -> {t}")
