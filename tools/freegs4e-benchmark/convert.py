"""freegs4e profiles -> MEQ profile tables.

THE UNITS ARE NOT WHAT THE DOCSTRING SAYS, AND THIS IS THE WHOLE FILE.

freegs4e's GeneralPprimeFFprime documents its inputs as "dp/dpsi_n" and
"fdf/dpsi_n" -- against NORMALISED flux.  THAT DOCSTRING IS WRONG.  They are
d/dpsi against the unnormalised psi in Wb/rad.  Four independent checks, the
last decisive:

  * the class docstring's own Jtor = R p'( psi ) + FF'( psi )/( R mu0 );
  * equilibrium.py's rhs = -mu0 R Jtor with Delta* as the operator, which is
    Grad-Shafranov only if p' = dp/dpsi;
  * base Profile.pressure() integrates pprime over psi_n and then multiplies by
    ( psi_axis - psi_bndry ), commented "convert from integral in normalised psi
    to integral in psi" -- that factor exists precisely because pprime() is
    already d/dpsi;
  * MEASURED: rebuilding Jtor from the saved arrays reproduces the solver's own
    Jtor to 4e-16, while reading them as d/dpsi_n is wrong by 4.3x to 55x.

So the VALUE columns carry over UNCHANGED.  Only the abscissa moves, and only
the derivative column divides by D:

    psi_MEQ = psi_fgs - psi_bndry = ( psi_n - 1 ) D,   D = psi_bndry - psi_axis
    dp/dpsi        -> unchanged
    d2p/dpsi2      = ( d/dpsi_n of the column ) / D

An earlier version of this file divided the value column by D as well, on the
strength of the docstring.  It would have parsed, solved and converged to an
equilibrium wrong by a factor of D -- which for a tokamak is order 1 and would
not have looked wrong.  CLAUDE.md records the same trap from the other side in
examples/rotating-density.dat.

INTERPOLATE WITH A CUBIC SPLINE.  freegs4e fits its own UnivariateSpline with
scipy's default ABSOLUTE smoothing budget, so what the solver evaluated is the
spline's values rather than the analytic table fed in -- for ff', whose data is
O(1), the budget never binds and the fit collapses to a single global cubic 1.3
to 1.7 percent from the tabulation.  The npz's `pprime`/`ffprime` are the
spline's own values and are what to use; `*_tabulated` is the analytic input and
is NOT what was solved.
"""
import numpy as np
from scipy.interpolate import CubicSpline


def to_meq_table(psi_n, values, psi_axis, psi_bndry, psi_zero=None):
    """(psi_n, dX/dpsi) -> (psi_MEQ, dX/dpsi, d2X/dpsi2), ascending in psi_MEQ.

    `psi_zero` is the flux MEQ's psi = 0 corresponds to -- the surface actually
    handed to MEQ as Gamma, which is NOT the separatrix here.  Defaults to
    psi_bndry.

    GETTING THIS WRONG DOES NOT LOOK WRONG.  An earlier version hardcoded
    psi_MEQ = ( psi_n - 1 ) D, i.e. psi = 0 at psi_n = 1, while the boundary
    handed to MEQ was psi_n = 0.9.  The profile then sits on an abscissa
    stretched by 1/0.9 and, far worse, MEQ's boundary lands where the source
    has nearly died: p' at the mapped psi = 0 was 2.7 instead of 3.1e5.  The
    solve converged in three Newton steps to a field 220x too small, entirely
    self-consistently -- psi stays near zero, so the source stays near zero.
    A trivial branch reached by a coordinate error rather than by the physics.
    """
    psi_n = np.asarray(psi_n, float)
    values = np.asarray(values, float)
    D = psi_bndry - psi_axis
    if D == 0.0:
        raise ValueError("psi_bndry equals psi_axis: the normalisation is degenerate")
    if psi_zero is None:
        psi_zero = psi_bndry
    level = (psi_zero - psi_axis) / D

    spline = CubicSpline(psi_n, values)
    dvalues = spline(psi_n, 1)                 # d/dpsi_n of the column

    psi_meq = (psi_n - level) * D
    column = values                            # ALREADY d/dpsi -- see the header
    dcolumn = dvalues / D                      # chain rule, once

    order = np.argsort(psi_meq)
    return psi_meq[order], column[order], dcolumn[order]


def to_meq_normalised_table(psi_n, values, level):
    """(psi_n, dX/dpsi) -> (Psi, dX/dPsi, d2X/dPsi2) for `Normalised = true`.

    WHY THE NORMALISED FORM IS THE RIGHT ONE HERE, AND IT IS NOT A PREFERENCE.

    `to_meq_table` above writes the profile against psi in Wb/rad, which is
    `Type = "mhd"` with `Normalised = false`.  That table necessarily STOPS at
    the axis: freegs4e's profiles are defined on psi_n in [ 0, 1 ] and there is
    no psi_n < 0.  meq::SplineProfile's documented out-of-range policy is to
    extend by a CONSTANT, so above the axis flux the source becomes a constant
    -- and that plateau manufactures a solution of its own.

    MEASURED on machine A, sweeping only the [initialguess] ramp amplitude and
    changing nothing else, as a ratio to the reference's own axis height:

        1x, 2x, 3x -> 0.0105        a small solution
        4x to 8x   -> 1.2604        the plateau's solution
        12x        -> 0.9996        the equilibrium

    Three branches, chosen non-monotonically by a guess.  That is not something
    to ship in an example.

    THE NORMALISED FORM HAS NO ABOVE-THE-AXIS REGION AT ALL.  Psi = psi/psi_ax
    with psi_ax an UNKNOWN of the bordered Newton, constrained to be max psi --
    so Psi <= 1 by construction and the profile is never evaluated off its own
    table.  It is also the form freegs4e itself solves in, and the form an EQDSK
    carries.  The clamp plateau is gone because the region it lived in does not
    exist.

    THE ALGEBRA, which has one factor in it that is easy to drop.  MEQ's
    boundary is the psi_n = `level` surface, so

        Psi = 1 - psi_n/level,      psi_MEQ = Psi * psi_ax_MEQ

    and meq::NormalisedMHDSource evaluates

        F = [ mu0 r^2 ( dp/dPsi )( Psi ) + ( g dg/dPsi )( Psi ) ] / psi_ax

    -- so the VALUE column is d/dPsi, not d/dpsi, and the two differ by
    psi_ax_MEQ = -level * D.  The division by psi_ax inside the source then
    takes it straight back to dp/dpsi, which is why this is the same physical
    source written in the coordinate the solve actually uses.

    @param level  the psi_n of the surface handed to MEQ as Gamma.
    @return       ( Psi, dX/dPsi / |scale|, d2X/dPsi2 ) with Psi ASCENDING, over
                  Psi in [ 0, 1 ] only -- the table is CUT at Gamma rather than
                  carried past it, because psi_n > level is outside the
                  computational domain and tabulating it would invite exactly
                  the extrapolation this function exists to remove.
    """
    psi_n = np.asarray(psi_n, float)
    values = np.asarray(values, float)
    if len(psi_n) < 4:
        raise ValueError("fewer than four profile points")

    # THE BOUNDARY IS INTERPOLATED IN RATHER THAN ROUNDED TO, WHICH IS WHAT
    # MAKES "never evaluated off its own table" EXACTLY TRUE.
    #
    # `level` lands between two of freegs4e's 256 psi_n samples -- 0.95 x 255
    # is 242.25 -- so simply keeping psi_n <= level leaves the table's lowest
    # Psi at 1.03e-03 instead of zero, and MEQ then CLAMPS across that sliver
    # at Gamma. It is a tiny window and the constant it clamps to is within a
    # per cent of the truth, so nothing would look wrong; but the whole reason
    # for the normalised form is that the profile is never extrapolated, and
    # "never, except in a sliver at the boundary" is a different claim.
    endpoint = CubicSpline(psi_n, values)(level)
    keep = psi_n < level - 1e-12
    psi_n = np.append(psi_n[keep], level)
    values = np.append(values[keep], endpoint)

    Psi = 1.0 - psi_n / level
    # And nail the endpoint against round-off, since the caller's abscissa has
    # to reach 0 exactly for the clamp never to fire.
    Psi[-1] = 0.0

    # d/dPsi of the VALUE, which is currently d/dpsi. The chain factor is
    # dpsi/dPsi = psi_ax_MEQ, and it is applied by the CALLER through `scale`
    # so that this function stays a pure coordinate change; see build().
    spline = CubicSpline(psi_n, values)
    dvalues_dpsin = spline(psi_n, 1)
    # dPsi/dpsi_n = -1/level, so d/dPsi = -level d/dpsi_n.
    dvalues = -level * dvalues_dpsin

    order = np.argsort(Psi)
    return Psi[order], values[order], dvalues[order]


DEFAULT_ABSCISSA = (
    "Three columns: psi [Wb/rad], value, d(value)/dpsi.\n"
    "psi is MEQ's, i.e. ZERO ON GAMMA and negative inside for a\n"
    "normal-sign tokamak, NOT normalised flux.")

NORMALISED_ABSCISSA = (
    "Three columns: Psi, value, d(value)/dPsi, ascending in Psi and\n"
    "reaching 0 and 1 EXACTLY -- so meq::SplineProfile interpolates\n"
    "everywhere the solve can reach and its clamp never fires.")


def write_meq_profile(path, psi, value, derivative, what, source_note,
                      abscissa=DEFAULT_ABSCISSA):
    with open(path, "w") as f:
        f.write(f"# {what}\n#\n")
        for line in source_note.strip().split("\n"):
            f.write(f"# {line}\n")
        f.write("#\n")
        for line in abscissa.strip().split("\n"):
            f.write(f"# {line}\n")
        f.write("#\n")
        for p, v, d in zip(psi, value, derivative):
            f.write(f"{p: .12e}  {v: .12e}  {d: .12e}\n")


if __name__ == "__main__":
    # SELF-TEST: the derivative column must be the derivative of the value
    # column IN psi_MEQ.  That is exactly what a dropped chain rule breaks, and
    # it is checkable without MEQ by differencing the table against itself.
    print("  PROFILE CONVERSION SELF-TEST\n")
    psi_axis, psi_bndry = -0.37, 0.11          # deliberately not 0 and 1
    psi_n = np.linspace(0.0, 1.0, 256)

    # An analytic dp/dpsi_n whose psi_n-derivative is known exactly.
    # A SMOOTH profile: beta integer, so every derivative is bounded on [0,1].
    alpha, beta, P0 = 2.0, 2.0, 3.7e5
    pprime_n = P0 * (1.0 - psi_n ** alpha) ** beta
    dpprime_n = (P0 * beta * (1.0 - psi_n ** alpha) ** (beta - 1.0)
                 * (-alpha * psi_n ** (alpha - 1.0)))

    psi, col, dcol = to_meq_table(psi_n, pprime_n, psi_axis, psi_bndry)
    D = psi_bndry - psi_axis

    want_col = np.interp(psi, (psi_n - 1.0) * D, pprime_n)
    e_val = np.abs(col - want_col).max() / np.abs(col).max()

    want_d = np.interp(psi, (psi_n - 1.0) * D, dpprime_n / D)
    e_der = np.abs(dcol - want_d).max() / np.abs(want_d).max()

    # The check that catches a dropped chain rule: difference the value column
    # in psi_MEQ and compare against the derivative column.
    mid = 0.5 * (psi[1:] + psi[:-1])
    fd = np.diff(col) / np.diff(psi)
    interp_d = np.interp(mid, psi, dcol)
    e_fd = np.abs(fd - interp_d).max() / np.abs(interp_d).max()

    print(f"    D = psi_bndry - psi_axis = {D:.6f}")
    print(f"    value column unchanged, vs input      {e_val:.3e}")
    print(f"    deriv column vs closed form / D       {e_der:.3e}")
    print(f"    deriv column vs difference of value   {e_fd:.3e}")

    # CONTROL 1: the chain rule. Drop the division by D and the difference
    # check must fail, or it is not testing anything.
    bad = dcol * D
    e_bad = np.abs(fd - np.interp(mid, psi, bad)).max() / np.abs(interp_d).max()
    print(f"    control, D dropped from the deriv     {e_bad:.3e}  (must be large)")

    # CONTROL 2: the units. Had the docstring been right and the value column
    # needed dividing by D, the two would differ by exactly D -- which is the
    # error this file's header exists to prevent.
    print(f"    control, value/D would differ by      {abs(1.0/D - 1.0):.3e}"
          f"  (the factor the docstring would have cost)")

    # AND THE PROFILE-SHAPE CONTROL: a fractional power at the edge cannot have
    # its derivative recovered from a table -- 5.6e-02 at the last node against
    # 2.6e-04 in the interior. Choose smooth tabulations.
    rough = P0 * (1.0 - psi_n ** 2.0) ** 1.5
    drough = P0 * 1.5 * (1.0 - psi_n ** 2.0) ** 0.5 * (-2.0 * psi_n)
    _, _, dr = to_meq_table(psi_n, rough, psi_axis, psi_bndry)
    wr = np.interp(psi, (psi_n - 1.0) * D, drough / D)
    rel = np.abs(dr - wr) / np.abs(wr).max()
    print(f"    fractional-power edge: last node      {rel[-1]:.3e}"
          f"   interior {rel[3:-3].max():.3e}")

    ok = e_val < 1e-12 and e_der < 1e-5 and e_fd < 1e-2 and e_bad > 0.1
    print("\n    SELF-TEST", "PASSED" if ok else "FAILED")
