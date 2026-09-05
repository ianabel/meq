"""Extract a closed interior flux surface, avoiding the X-point corner.

WHY NOT THE SEPARATRIX. Every reference case is DIVERTED, so its last closed
surface passes through an X-point and has a CORNER there. MXH is a smooth
parametrisation -- tR( t ) is a truncated Fourier series -- so it cannot
represent one, and fitting the separatrix reads 2.6e-03 to 1.3e-02 m against a
fitter that does 1.9e-05 on a smooth shape. The comparison would then be
measuring MXH's inability to turn a corner.

MEQ's own tree makes the same choice for the same reason: ExtensionConvergence
takes Gamma to be psi = -0.03 rather than the separatrix because the separatrix
"passes through an X-point -- a CORNER of Gamma, where both transfer-path
families give out and the Cockburn-Solano analysis does not reach".

So take an interior surface, psi_n = 0.95 by default. It is smooth, closed, and
strictly inside the core, so the source there is the profile with no masking.
"""
import numpy as np
from contourpy import contour_generator


def interior_surface(R, Z, psi, psi_axis, psi_bndry, Raxis, Zaxis, level=0.95):
    """The closed psi_n = `level` contour enclosing the magnetic axis."""
    psi_target = psi_axis + level * (psi_bndry - psi_axis)

    # psi is stored ( R, Z ); contour_generator wants ( y, x ) = ( Z, R ).
    gen = contour_generator(x=R, y=Z, z=np.asarray(psi, float).T)
    lines = gen.lines(psi_target)

    best = None
    for seg in lines:
        seg = np.asarray(seg, float)
        if len(seg) < 8:
            continue
        closed = np.hypot(*(seg[0] - seg[-1])) < 1e-9 * max(1.0, np.ptp(seg[:, 0]))
        if not closed:
            continue
        # Must enclose the axis: a diverted case also produces private-flux
        # loops at the same level, and picking one of those silently gives a
        # tiny plausible-looking boundary.
        if not _encloses(seg, Raxis, Zaxis):
            continue
        span = np.ptp(seg[:, 0]) * np.ptp(seg[:, 1])
        if best is None or span > best[0]:
            best = (span, seg)

    if best is None:
        raise RuntimeError(f"no closed psi_n = {level} contour encloses the axis")

    seg = best[1][:-1]              # drop the repeated endpoint
    return seg[:, 0], seg[:, 1], psi_target


def _encloses(seg, r, z):
    """Winding of the polygon about ( r, z ), by the turning of the ray."""
    ang = np.unwrap(np.arctan2(seg[:, 1] - z, seg[:, 0] - r))
    return abs(ang[-1] - ang[0]) > np.pi
