"""Miller eXtended Harmonic fitting, and its self-test.

MEQ's [boundary.shape] Type = "mxh" takes R0, Z0, MinorRadius, Elongation plus
CosCoefficients c_0..c_N and SinCoefficients s_1..s_N.  The parametrisation,
Arbon-Candy-Belli:

    R( t ) = R0 + a cos( tR( t ) )        Z( t ) = Z0 + kappa a sin( t )
    tR( t ) = t + c_0 + sum_n [ c_n cos( n t ) + s_n sin( n t ) ]

so t is set by Z and tR by R.  Note c STARTS AT c_0 (the tilt) and s STARTS AT
s_1 -- that asymmetry is real and MEQ checks it on load.
"""
import numpy as np


def mxh_boundary(R0, Z0, a, kappa, cos_coeffs, sin_coeffs, nt=512):
    """Points on an MXH surface, counterclockwise from the outboard midplane."""
    t = np.linspace(0.0, 2.0 * np.pi, nt, endpoint=False)
    tR = t + cos_coeffs[0]
    for n in range(1, len(cos_coeffs)):
        tR = tR + cos_coeffs[n] * np.cos(n * t)
    for n in range(1, len(sin_coeffs) + 1):
        tR = tR + sin_coeffs[n - 1] * np.sin(n * t)
    return R0 + a * np.cos(tR), Z0 + kappa * a * np.sin(t)


def fit_mxh(R, Z, n_harmonics=6):
    """Fit an MXH surface to an ordered closed contour.

    THE BRANCH IS THE WHOLE DIFFICULTY.  sin t = ( Z - Z0 )/( kappa a ) fixes t
    only up to t <-> pi - t, and cos tR = ( R - R0 )/a fixes tR only up to sign.
    Both are resolved by UNWRAPPING along the contour rather than by quadrant
    rules, which fail wherever the surface is non-convex enough that R is not
    monotone between the extrema -- exactly the shaped cases this is for.
    """
    R = np.asarray(R, float)
    Z = np.asarray(Z, float)

    R0 = 0.5 * (R.max() + R.min())
    a = 0.5 * (R.max() - R.min())
    Z0 = 0.5 * (Z.max() + Z.min())
    kappa = 0.5 * (Z.max() - Z.min()) / a

    # Order counterclockwise starting at the outboard midplane, so that both
    # angles increase together and the unwrap below has a direction to follow.
    geom = np.arctan2(Z - Z0, R - R0)
    order = np.argsort(np.mod(geom, 2.0 * np.pi))
    R, Z = R[order], Z[order]

    y = np.clip((Z - Z0) / (kappa * a), -1.0, 1.0)
    x = np.clip((R - R0) / a, -1.0, 1.0)

    # THE BRANCHES SWITCH AT THE EXTREMA, NOT AT THE MIDPLANE.
    #
    # cos tR = x, and tR runs 0 -> 2pi once, so sin tR changes sign exactly
    # where x is stationary -- at the R extrema. sin t = y likewise, so cos t
    # changes sign at the Z extrema. Using sign( Z - Z0 ) for the tR branch is
    # right only for an up-down symmetric surface, and it is what put spurious
    # c_3..c_6 at 5e-4 into a fit whose truth was exactly zero: Z0 offset and
    # s_1 both move the R extrema off the midplane.
    npt = len(R)
    i_rmax, i_rmin = int(np.argmax(R)), int(np.argmin(R))
    i_zmax, i_zmin = int(np.argmax(Z)), int(np.argmin(Z))

    def on_arc(i, lo, hi):
        """True if index i lies on the ccw arc from lo to hi."""
        return ((i - lo) % npt) < ((hi - lo) % npt)

    upper = np.array([on_arc(i, i_rmax, i_rmin) for i in range(npt)])
    tR = np.where(upper, np.arccos(x), 2.0 * np.pi - np.arccos(x))

    rising = np.array([on_arc(i, i_zmin, i_zmax) for i in range(npt)])
    t = np.where(rising, np.arcsin(y), np.pi - np.arcsin(y))

    tR = np.unwrap(tR)
    t = np.unwrap(t)
    # Both start near zero at the outboard midplane; remove any 2pi offset
    # introduced by where the sort happened to begin.
    tR = tR - 2.0 * np.pi * np.round((tR[0] - t[0]) / (2.0 * np.pi))
    t = t - 2.0 * np.pi * np.round(t[0] / (2.0 * np.pi))

    delta = tR - t

    # ONTO A UNIFORM t FIRST, AND THAT IS NOT A REFINEMENT BUT THE WHOLE
    # ACCURACY.  t comes from the geometry, so it is NOT uniformly spaced, and a
    # trapezoid over non-uniform abscissae is second order -- which shows up as
    # spurious high harmonics at 1e-4 where the truth is exactly zero.  The
    # periodic trapezoid on a UNIFORM grid is spectral, so interpolate first.
    t0 = t[0]
    tw = np.mod(t - t0, 2.0 * np.pi)
    idx = np.argsort(tw)
    tw, dw = tw[idx], delta[idx]
    # Close the period for the interpolation.
    tw = np.concatenate([tw, [tw[0] + 2.0 * np.pi]])
    dw = np.concatenate([dw, [dw[0]]])

    m = 2048
    tu = np.linspace(0.0, 2.0 * np.pi, m, endpoint=False)
    du = np.interp(tu, tw, dw)
    tabs = tu + t0

    c = [du.mean()]
    s = []
    for n in range(1, n_harmonics + 1):
        c.append(2.0 * np.mean(du * np.cos(n * tabs)))
        s.append(2.0 * np.mean(du * np.sin(n * tabs)))
    return dict(R0=R0, Z0=Z0, a=a, kappa=kappa,
                cos=np.array(c), sin=np.array(s))


def worst_distance(R, Z, fit, nt=200000):
    """Worst distance from each given point to the fitted surface.

    THE SAMPLING IS THE INSTRUMENT AND IT HAS TO BE FINER THAN THE ANSWER.
    A min-distance to a curve sampled at nt points cannot read below about
    2 pi a / nt -- at nt = 2048 on a = 0.6 that is 1.8e-3 m, and a fit good to
    1e-5 reads 1.6e-3 no matter how good it is.  Measured: improving the
    coefficients 25-fold left this number unchanged, which is the signature.
    Chunked so the pairwise array stays bounded.
    """
    fR, fZ = mxh_boundary(fit["R0"], fit["Z0"], fit["a"], fit["kappa"],
                          fit["cos"], fit["sin"], nt)
    worst = 0.0
    for i in range(0, len(R), 64):
        d = np.hypot(R[i:i + 64, None] - fR[None, :],
                     Z[i:i + 64, None] - fZ[None, :])
        worst = max(worst, float(d.min(axis=1).max()))
    return worst


def pointwise_error(fit, truth, nt=400):
    """Compare two MXH surfaces at the SAME parameter, which needs no search."""
    aR, aZ = mxh_boundary(fit["R0"], fit["Z0"], fit["a"], fit["kappa"],
                          fit["cos"], fit["sin"], nt)
    bR, bZ = mxh_boundary(truth["R0"], truth["Z0"], truth["a"], truth["kappa"],
                          truth["cos"], truth["sin"], nt)
    return float(np.hypot(aR - bR, aZ - bZ).max())


if __name__ == "__main__":
    # SELF-TEST: fit a surface whose coefficients are known, and demand them
    # back.  Without this the fitter is untested code sitting between two
    # solvers, and a wrong shape converges perfectly well to a wrong answer.
    truth = dict(R0=1.7, Z0=0.05, a=0.6, kappa=1.8,
                 cos=np.array([0.02, -0.01, 0.004, 0.0, 0.0, 0.0, 0.0]),
                 sin=np.array([0.31, -0.06, 0.012, 0.0, 0.0, 0.0]))
    R, Z = mxh_boundary(truth["R0"], truth["Z0"], truth["a"], truth["kappa"],
                        truth["cos"], truth["sin"], 400)
    got = fit_mxh(R, Z, n_harmonics=6)

    print("  MXH FITTER SELF-TEST\n")
    print(f"    {'quantity':14s} {'fitted':>14s} {'truth':>14s} {'error':>12s}")
    # THE TOLERANCE IS SET BY WHAT LIMITS THE FIT, WHICH IS THE EXTREMA.
    # R0, a, Z0 and kappa are read off the DISCRETE samples, so each is wrong by
    # about ( d theta )^2 * curvature -- 4e-6 to 1.2e-5 here at 400 points --
    # and that normalisation error propagates into every coefficient. It is not
    # the branch logic, which was worth 25x when it was fixed, and it is not the
    # quadrature, which is spectral on a uniform grid.
    #
    # 1e-4 is therefore honest rather than lenient, and it is far below what the
    # input carries: an LCFS extracted from a 129-point grid is good to nothing
    # like this. Tightening it means refining the extrema sub-sample, which is
    # worth doing only if the input ever justifies it.
    ok = True
    for key in ("R0", "Z0", "a", "kappa"):
        err = abs(got[key] - truth[key])
        ok &= err < 1e-4
        print(f"    {key:14s} {got[key]:14.8f} {truth[key]:14.8f} {err:12.2e}")
    for n in range(len(truth["cos"])):
        err = abs(got["cos"][n] - truth["cos"][n])
        ok &= err < 1e-4
        print(f"    c_{n:<12d} {got['cos'][n]:14.8f} {truth['cos'][n]:14.8f} {err:12.2e}")
    for n in range(len(truth["sin"])):
        err = abs(got["sin"][n] - truth["sin"][n])
        ok &= err < 1e-4
        print(f"    s_{n+1:<12d} {got['sin'][n]:14.8f} {truth['sin'][n]:14.8f} {err:12.2e}")
    d = worst_distance(R, Z, got)
    e = pointwise_error(got, truth)
    print(f"\n    worst distance, point to fitted surface: {d:.3e} m")
    print(f"    worst separation of the two surfaces:    {e:.3e} m")
    print("\n    SELF-TEST", "PASSED" if ok and e < 1e-4 else "FAILED")
