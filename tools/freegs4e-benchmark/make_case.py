"""Turn one freegs4e reference .npz into a MEQ run: TOML plus two profile tables."""
import json
import os
import numpy as np

from mxh import fit_mxh, worst_distance
from surface import interior_surface
from convert import to_meq_table, write_meq_profile


def build(npz_path, outdir, degree=2, nr=12, nz=12, refine=2, harmonics=10,
          grid=129, level=0.90):
    d = np.load(npz_path, allow_pickle=True)
    stem = os.path.splitext(os.path.basename(npz_path))[0]
    os.makedirs(outdir, exist_ok=True)

    psi_axis = float(d["psi_axis"])
    psi_bndry = float(d["psi_bndry"])
    D = psi_bndry - psi_axis

    # AN INTERIOR SURFACE, NOT THE SEPARATRIX. Every reference case is diverted,
    # so its LCFS has an X-point CORNER and MXH -- a truncated Fourier series in
    # tR -- cannot represent one: fitting it reads 2.6e-03 to 1.3e-02 m against
    # a fitter measured at 1.9e-05 on a smooth shape. MEQ's own tree makes the
    # same choice for the same reason; see surface.py.
    #
    # 10 harmonics because the fit PLATEAUS there at 2-4e-04 m, which is the
    # floor of the contour extracted from a 129^2 grid rather than the fitter's.
    # 16 and 24 buy nothing, and that flatness is what says where the floor is.
    sR, sZ, psi_surface = interior_surface(
        d["R"], d["Z"], d["psi"], psi_axis, psi_bndry,
        float(d["Raxis"]), float(d["Zaxis"]), level)

    fit = fit_mxh(sR, sZ, n_harmonics=harmonics)
    shape_err = worst_distance(np.asarray(sR, float), np.asarray(sZ, float), fit)

    note = (f"Converted from {os.path.basename(npz_path)}.\n"
            f"freegs4e tabulates against NORMALISED flux; this is against psi.\n"
            f"psi_axis = {psi_axis:.9e}, psi_surface = {psi_surface:.9e},\n"
            f"the psi_n = {level} surface, chosen to avoid the X-point corner.\n"
            f"freegs4e's pprime/ffprime are ALREADY d/dpsi -- its docstring\n"
            f"says d/dpsi_n and is wrong -- so the value column is unchanged\n"
            f"and only the derivative column carries a chain-rule factor.")

    # psi_MEQ = psi_fgs - psi_surface, so MEQ's psi vanishes on the
    # boundary it is actually given rather than on the separatrix.
    pp = to_meq_table(d["psi_n"], d["pprime"], psi_axis, psi_bndry, psi_surface)
    gg = to_meq_table(d["psi_n"], d["ffprime"], psi_axis, psi_bndry, psi_surface)
    pp_path = os.path.join(outdir, f"{stem}-pprime.dat")
    gg_path = os.path.join(outdir, f"{stem}-ggprime.dat")
    write_meq_profile(pp_path, *pp, "dp/dpsi [Pa per Wb/rad]", note)
    write_meq_profile(gg_path, *gg, "g dg/dpsi [T^2 m^2 per Wb/rad]", note)

    # A box comfortably around the fitted surface.
    R0, Z0, a, kappa = fit["R0"], fit["Z0"], fit["a"], fit["kappa"]
    pad_r, pad_z = 0.35 * a, 0.35 * kappa * a
    box = dict(rmin=max(1e-3, R0 - a - pad_r), rmax=R0 + a + pad_r,
               zmin=Z0 - kappa * a - pad_z, zmax=Z0 + kappa * a + pad_z)

    cos_list = ", ".join(f"{v: .10e}" for v in fit["cos"])
    sin_list = ", ".join(f"{v: .10e}" for v in fit["sin"])

    toml = f"""# MEQ against freegs4e: {stem}
#
# GENERATED -- edit make_case.py, not this file.
#
# The boundary is freegs4e's converged LCFS fitted to MXH; the fit is good to
# {shape_err:.3e} m, against a minor radius of {a:.4f} m. The source is
# freegs4e's own tabulated p' and ff', converted from normalised flux to psi.
# psi is MEQ's here: ZERO ON GAMMA, so psi_MEQ = psi_freegs4e - psi_bndry.

[mesh]
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
Type = "zero"

[boundary.shape]
Type = "mxh"
R0 = {R0:.10f}
Z0 = {Z0:.10f}
MinorRadius = {a:.10f}
Elongation = {kappa:.10f}
CosCoefficients = [ {cos_list} ]
SinCoefficients = [ {sin_list} ]

# THE RAMP IS NOT A NICETY HERE, IT PICKS THE BRANCH.
#
# These sources do not vanish at psi = 0 -- p' at the boundary surface is 2.7
# and gg' is 4.7e-3 -- so this is not the trivial branch CLAUDE.md's Traps
# section is about. It is the SECOND SOLUTION beside it, which that file records
# for the high-beta source: "this equation has a small positive solution and a
# large one ... Newton from the Dirichlet datum walks straight onto it".
#
# Measured here without a ramp: case A converged in THREE Newton steps, cleanly,
# to a field 23x too small; E to one 3.7x too small; G did not converge at all.
# Case B converged to the right one unaided. So the amplitude is part of the
# problem statement rather than a tuning parameter, and it is set to the axis
# height the reference actually has.
[initialguess]
Type = "ramp"
Amplitude = {6.0 * (psi_axis - psi_surface):.10e}

[source]
Type = "mhd"
PPrimeFile = "{os.path.abspath(pp_path)}"
GGPrimeFile = "{os.path.abspath(gg_path)}"

[output]
Directory = "{os.path.abspath(outdir)}"
Prefix = "{stem}"
GridNR = {grid}
GridNZ = {grid}
"""
    toml_path = os.path.join(outdir, f"{stem}.toml")
    with open(toml_path, "w") as f:
        f.write(toml)

    meta = dict(stem=stem, psi_axis=psi_axis, psi_bndry=psi_bndry, D=D,
                psi_surface=psi_surface, level=level,
                shape_error_m=shape_err, minor_radius=a, R0=R0, Z0=Z0,
                kappa=kappa, box=box,
                cos=list(map(float, fit["cos"])),
                sin=list(map(float, fit["sin"])))
    with open(os.path.join(outdir, f"{stem}-meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    return toml_path, meta


if __name__ == "__main__":
    import sys
    for path in sys.argv[1:]:
        t, m = build(path, os.path.dirname(path) or ".")
        print(f"{m['stem']:24s} shape fit {m['shape_error_m']:.3e} m "
              f"(a = {m['minor_radius']:.4f})  D = {m['D']:.4e}  -> {t}")
