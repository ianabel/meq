"""
Free-boundary Grad-Shafranov reference equilibria from freegs4e.

Produces converged free-boundary equilibria with a TABULATED (p', ff') source
(jtor.GeneralPprimeFFprime), so the current source is an explicitly known
function of normalised flux, and dumps each to a portable .npz + JSON sidecar.

Run as:
    PYTHONPATH=/home/ian/projects/freegs4e <venv-python> fgsref.py [case ...]

Three things about freegs4e had to be worked around; all are documented at the
point of use below:

  1. Machine.controlAdjust() updates coil.current but NOT Machine.current_vec,
     while Equilibrium.psi() reads the flux through current_vec.  Left alone,
     every equilibrium is solved with ZERO coil flux on the grid.
  2. Profile classes in this fork have no .Jtor() entry point (only the
     FreeGSNKE-style Jtor_part1/Jtor_part2 pair), which Equilibrium.solve()
     requires.
  3. Equilibrium.__init__ never calls _updatePlasmaPsi, so eq.psi_func does
     not exist and eq.Br()/eq.Bz() raise until it is called by hand.
"""

import json
import os
import sys
import time
import warnings

import numpy as np
from matplotlib.path import Path as MplPath
from scipy import interpolate as sinterp

from freegs4e import machine, jtor, equilibrium, control, critical

import conductors
from freegs4e.gradshafranov import GSElliptic, GSsparse4thOrder, mu0

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Multigrid levels for the linear solve, or None for freegs4e's own default.
#
# ITS DEFAULT IS nlevels = 1, WHICH IS NOT MULTIGRID AT ALL.
# Equilibrium.__init__ calls multigrid.createVcycle( ..., nlevels=1, ncycle=1,
# niter=2, direct=True ), and at one level a V-cycle is a direct sparse solve on
# the full grid.  The whole hierarchy is built and switched off.
#
# BUT TURNING IT ON BUYS ALMOST NOTHING, AND THIS COMMENT USED TO SAY OTHERWISE.
# It attributed case A's 85 s at 129^2 against 656 s at 257^2 to "the n^3 of a
# 2D sparse LU".  That is wrong, and the HAGENOW note below had it right all
# along: multigrid.MGDirect.__init__ calls scipy factorized( A ) ONCE at
# construction and __call__ only backsolves, so the LU is not a per-step cost at
# all.  Measured 2026-09-06, the per-Picard-step boundary condition alone is
# 0.512 s at 129^2 and 3.540 s at 257^2 -- against those 85 s and 656 s totals,
# it is most of the run.  The V-cycle is a lever on the wrong term.
#
# AND setSolverVcycle() IS NOT THE WAY TO TURN IT ON, which is a trap rather
# than an inconvenience: it hard-codes GSsparse -- the SECOND-order operator --
# and ignores Equilibrium.order, so calling it on a 4th-order equilibrium
# silently solves a different problem.  This builds the V-cycle directly on the
# generator the case actually asked for.
VCYCLE_LEVELS = None

# Use freegs4e's von Hagenow free-boundary condition instead of its default.
#
# ITS DEFAULT IS THE NAIVE ONE AND THAT IS WHAT COSTS THE REFINEMENT.
# Equilibrium.__init__ takes boundary=freeBoundary, whose own docstring calls it
# "an integral over the area of the domain for each point": it loops over the
# 4n boundary points and, for each, evaluates Greens over the whole n^2 grid.
# That is O( n^3 ) PER PICARD STEP, in numpy, and it is what makes the reference
# cost 7.7x for a fourfold rise in unknowns -- not the linear solve, which is a
# sparse LU factorised once.
#
# boundary.freeBoundaryHagenow is in the same file, is the method the paper this
# benchmark cites is built on, and is O( n^2 ) -- boundary work only.  It is not
# the default and nothing points at it.
HAGENOW = False

# The observation-point offset von Hagenow's boundary integral uses, as a
# MULTIPLE OF THE CELL, or None to leave freegs4e's shipped constant alone.
# See WORKAROUND 5.
HAGENOW_EPS = 0.2

# Cache the DEFAULT boundary condition as the fixed matrix it is.  See
# WORKAROUND 6.  This changes no answer -- it is the same arithmetic
# reassociated -- so it is on by default and CACHE_BOUNDARY = False is the
# control.
CACHE_BOUNDARY = True

# Above this the matrix is not built and the shipped loop runs instead: it is
# 2( nx + ny ) x nx ny doubles, which is 69 MB at 129^2, 0.54 GB at 257^2,
# 4.3 GB at 513^2 and 35 GB at 1025^2.  The cap is a guard against turning a
# slow run into an OOM on a shared machine, not a tuning parameter.
CACHE_BOUNDARY_MAX_GB = 6.0

# --seed-from=DIR: the directory holding a COARSER run's .npz files, whose
# converged answer seeds this one.  None means every case starts cold, which
# is what it did before 2026-09-06.  See seed_from_coarse().
SEED_FROM = None


def install_vcycle(eq, order, levels):
    """A V-cycle on the generator matching @a order, which is what
    Equilibrium.setSolverVcycle does not do."""
    from freegs4e import multigrid
    from freegs4e.gradshafranov import GSsparse, GSsparse4thOrder
    gen = (GSsparse4thOrder if order == 4 else GSsparse)(
        eq.Rmin, eq.Rmax, eq.Zmin, eq.Zmax)
    eq._solver = multigrid.createVcycle(
        eq.nx, eq.ny, gen, nlevels=levels, ncycle=2, niter=20, direct=True)
    return gen
NPSI = 256
RTOL = 1.0e-9      # Picard: relative change in psi
MAXITS = 400


# ----------------------------------------------------------------------------
# WORKAROUND 0
# critical.inside_mask() tests `if len(xpoint > 1)` where `if len(xpoint) > 1`
# was meant.  `xpoint > 1` is an elementwise comparison, so len() of it is the
# NUMBER of X-points -- truthy for exactly one -- and the body then indexes
# xpoint[1, 2] and raises IndexError.  Any single-null equilibrium therefore
# crashes, including inside Equilibrium._updatePlasmaPsi, which is on the
# Picard path.  Reinstalled here with the one character fixed and nothing else
# changed.
# ----------------------------------------------------------------------------
_inside_mask_orig = critical.inside_mask


def _inside_mask_fixed(R, Z, psi, opoint, xpoint=[], mask_outside_limiter=None,
                       psi_bndry=None, use_geom=True):
    mask = critical.inside_mask_(
        R, Z, psi, opoint, xpoint, mask_outside_limiter, psi_bndry
    )
    if use_geom and len(xpoint) > 0:
        mask = mask * critical.geom_inside_mask(R, Z, opoint, xpoint)
        if len(xpoint) > 1:                      # <-- was len(xpoint > 1)
            if (
                np.abs((xpoint[0, 2] - xpoint[1, 2])
                       / (opoint[0, 2] - xpoint[0, 2])) < 0.1
            ):
                mask = mask * critical.geom_inside_mask(R, Z, opoint,
                                                        xpoint[1:])
    return mask


critical.inside_mask = _inside_mask_fixed


# ----------------------------------------------------------------------------
# WORKAROUND 1
# ----------------------------------------------------------------------------
class SyncConstrain(control.constrain):
    """control.constrain that keeps Machine.current_vec in step with the coils.

    machine.Machine.controlAdjust() writes coil.current on each coil object.
    Machine.getPsitokamak(), which is what Equilibrium.psi() uses to add the
    coil flux to the grid, reads Machine.current_vec instead -- a separate
    array only refreshed by getCurrentsVec()/set_*_coil_current*().  So the
    grid psi keeps whatever coil currents current_vec was last set to (zeros,
    from Machine.__init__), while eq.Br()/eq.Bz()/eq.psiRZ(), which sum over
    the coil objects analytically, see the new ones.

    Symptom if this is not done: the constraint solve reports Br = Bz = 0 at
    the requested X-points to 1e-12, and critical.find_critical() finds no
    X-points at all, because the psi it is handed has no coil field in it.
    The Picard iteration then converges happily to a small plasma pinned
    against the outboard wall.
    """

    def __call__(self, eq):
        super().__call__(eq)
        eq.tokamak.getCurrentsVec()


# ----------------------------------------------------------------------------
# WORKAROUND 2
# ----------------------------------------------------------------------------
class FreeGSProfileMixin:
    """Supplies the Profile.Jtor() entry point this fork dropped.

    freegs4e splits the current calculation into Jtor_part1 (critical points
    and core mask) and Jtor_part2 (the current itself) and expects the driver
    -- FreeGSNKE -- to join them.  equilibrium.Equilibrium.solve(), which
    freegs4e's own picard.solve calls, still expects profiles.Jtor(R, Z, psi,
    psi_bndry=...).  No profile class in the fork provides it, so nothing in
    freegs4e's standalone solve path runs at all.

    The default behaviour here is exactly original FreeGS: the plasma boundary
    is the flux surface through the primary X-point, and if there is no
    X-point the current is left unmasked.  attach_limiter() optionally adds a
    FreeGSNKE-style wall limiter, in which case the boundary is whichever of
    the diverted and wall-limited surfaces encloses less plasma.
    """

    mask_inside_limiter = None
    mask_outside_limiter = None
    limiter_ring = None
    edge_mask = None
    diverted_core_mask = None
    flag_limiter = False

    def attach_limiter(self, R, Z, wall_R, wall_Z):
        """Build limiter masks from a wall polygon.

        Note the convention critical.inside_mask_() requires: the array handed
        in is copied into the flood-fill working array, blocks the fill
        wherever it is >= 0.5, and the answer is (work == 1).  The blocking
        value must therefore be 2, not 1, or every point outside the limiter
        comes back marked as core.
        """
        pts = np.column_stack([R.ravel(), Z.ravel()])
        poly = MplPath(np.column_stack([wall_R, wall_Z]))
        inside = poly.contains_points(pts).reshape(R.shape)
        self.mask_inside_limiter = inside
        self.mask_outside_limiter = 2.0 * (~inside)
        ring = inside & ~(
            np.roll(inside, 1, 0) & np.roll(inside, -1, 0)
            & np.roll(inside, 1, 1) & np.roll(inside, -1, 1)
        )
        self.limiter_ring = np.where(ring)
        self.edge_mask = ring

    def boundary_and_mask(self, R, Z, psi):
        sgn = np.sign(self.Ip)
        opt, xpt = critical.find_critical(
            R, Z, psi, self.mask_inside_limiter, self.Ip
        )
        if len(opt) == 0:
            raise RuntimeError("no O-point found")
        cands = []
        if len(xpt) > 0:
            cands.append(("diverted", float(xpt[0][2])))
        if self.limiter_ring is not None:
            ring_psi = psi[self.limiter_ring]
            cands.append(("limited",
                          float(ring_psi.max() if sgn > 0 else ring_psi.min())))
        if not cands:
            # original FreeGS cold-start fallback: no boundary, no mask
            return opt, xpt, float(psi[0, 0]), None, "unbounded"

        # psi_axis is an extremum of sign sgn, so the surface enclosing less
        # plasma is the one whose psi is further from the axis value
        if sgn > 0:
            kind, psi_b = max(cands, key=lambda t: t[1])
        else:
            kind, psi_b = min(cands, key=lambda t: t[1])

        if kind == "diverted":
            mask = critical.inside_mask(
                R, Z, psi, opt, xpt, self.mask_outside_limiter, psi_b,
                use_geom=True,
            )
        else:
            # the X-point, if any, is outside the plasma, so geom_inside_mask's
            # half-plane cut through it would amputate the core
            mask = critical.inside_mask_(
                R, Z, psi, opt, np.zeros((0, 3)),
                self.mask_outside_limiter, psi_b,
            )
        if self.mask_inside_limiter is not None:
            mask = mask * self.mask_inside_limiter
        return opt, xpt, psi_b, mask, kind

    def Jtor(self, R, Z, psi, psi_bndry=None):
        opt, xpt, psi_b, mask, kind = self.boundary_and_mask(R, Z, psi)
        if psi_bndry is not None:
            psi_b = psi_bndry
        self.opt, self.xpt = opt, xpt
        self.psi_axis = float(opt[0][2])
        self.diverted_core_mask = mask
        self.limiter_core_mask = mask
        self.boundary_kind = kind
        self.flag_limiter = kind == "limited"
        return self.Jtor_part2(R, Z, psi, self.psi_axis, psi_b, mask)


class GeneralProfile(FreeGSProfileMixin, jtor.GeneralPprimeFFprime):
    pass


class PaxisProfile(FreeGSProfileMixin, jtor.ConstrainPaxisIp):
    pass


# ----------------------------------------------------------------------------
# Picard loop.  Semantics identical to freegs4e.picard.solve (same convergence
# test, same ordering of constrain/solve/blend) but it returns a status instead
# of raising, so a non-converged case can still be reported rather than lost.
# ----------------------------------------------------------------------------
def good_enough(st):
    """Accept a Picard result.  The core mask is a set of grid cells, so one
    cell flipping in and out puts a floor of order (cell current)/(flux span)
    under the step size; below 1e-5 relative that floor, not the equilibrium,
    is what is being measured.  The real test of convergence is the GS
    residual computed afterwards from the final psi, which is reported
    separately and gated at 1e-6."""
    return st["status"] in ("rtol", "atol") or st["rel_history"][-1] < 1e-5


def picard_loop(eq, profiles, constrain, rtol=1e-9, atol=1e-12, blend=0.0,
                maxits=200):
    if constrain is not None:
        constrain(eq)
    psi = eq.psi()
    rel_hist, abs_hist = [], []
    it = 0
    status = "maxits"
    while True:
        psi_last = psi.copy()
        eq.solve(profiles, psi=psi)
        psi = eq.psi()
        d = psi_last - psi
        amax = float(np.max(np.abs(d)))
        arel = amax / (np.max(psi) - np.min(psi))
        abs_hist.append(amax)
        rel_hist.append(arel)
        if amax < atol:
            status = "atol"
            break
        if arel < rtol:
            status = "rtol"
            break
        if constrain is not None:
            constrain(eq)
        psi = (1.0 - blend) * eq.psi() + blend * psi_last
        it += 1
        if maxits and it > maxits:
            break
        eq._profiles = profiles
    return dict(status=status, iterations=len(rel_hist),
                rel_history=np.array(rel_hist), abs_history=np.array(abs_hist))


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# WORKAROUND 6
# boundary.freeBoundary REBUILDS A GEOMETRY-ONLY MATRIX ON EVERY PICARD STEP,
# and that -- not the linear solve -- is what makes a refined reference
# expensive.
#
# Its own docstring calls it "an integral over the area of the domain for each
# point on the boundary": for each of the 2( nx + ny ) boundary points it
# evaluates Greens over the whole nx x ny grid -- elliptic integrals, O( n^2 )
# of them, O( n ) times -- and contracts against Jtor with a Romberg rule.  So
# it is O( n^3 ) PER STEP.
#
# BUT Greens DEPENDS ONLY ON THE GEOMETRY AND romb IS LINEAR, so the whole
# boundary condition is ONE FIXED MATRIX applied to Jtor:
#
#     psi_bndry = M Jtor,     M[ b, ij ] = dR dZ w_i w_j G_b[ i, j ]
#
# with w the Romberg weights and G_b[ b ] = 0, exactly as the loop zeroes the
# self term.  M is the same at every iteration and is rebuilt at every
# iteration.  Building it once costs about one pass of the loop and every
# subsequent step is a matrix-vector product.
#
# THE WEIGHTS ARE TAKEN FROM scipy RATHER THAN REIMPLEMENTED.  romb is linear,
# so romb( I ) IS the weight vector: one call on an identity, and the rule can
# never drift from the one the shipped loop uses.  Reimplementing Romberg here
# would be a second definition of the quadrature and the whole point is that
# there is only one.
#
# WHAT IT COSTS IS MEMORY, WHICH IS WHY THERE IS A CAP.  M is
# 2( nx + ny ) x nx ny doubles: 69 MB at 129^2, 0.54 GB at 257^2, 4.3 GB at
# 513^2, 35 GB at 1025^2.  Above CACHE_BOUNDARY_MAX_GB the shipped loop runs
# instead, so the only consequence of a grid too large is that it is slow again.
#
# A REFINEMENT NOT TAKEN: Jtor is zero outside the plasma, so only its support's
# columns are needed and the matrix would shrink by about an order of magnitude.
# It costs tracking a support that MOVES as Picard runs, and the whole appeal of
# this is that M is a constant.
# ----------------------------------------------------------------------------
def _romb_weights(n):
    """The Romberg weight vector scipy applies, at unit spacing.

    romb is linear, so integrating the identity gives its weights directly.
    """
    from scipy.integrate import romb
    return romb(np.eye(n), axis=-1)


def _boundary_indices(nx, ny):
    """The boundary points, in freegs4e's own order.

    The corners appear TWICE, once from a horizontal edge and once from a
    vertical one. That is the shipped behaviour and it is harmless -- both
    write the same value to the same cell -- so it is reproduced rather than
    tidied: a different index set would be a different matrix.
    """
    return np.concatenate([
        [(x, 0) for x in range(nx)],
        [(x, ny - 1) for x in range(nx)],
        [(0, y) for y in range(ny)],
        [(nx - 1, y) for y in range(ny)],
    ]).astype(int)


def _boundary_matrix(eq):
    """M, built once and kept on the Equilibrium.

    Cached on the object rather than in a module dict so that it cannot outlive
    the grid it was built for: a stale M is a wrong answer, and a case-keyed
    global would be one edit away from producing one.
    """
    R = eq.R
    Z = eq.Z
    nx, ny = R.shape
    have = getattr(eq, "_meq_boundary_matrix", None)
    if have is not None and have[0].shape[1] == nx*ny:
        return have

    gb = 2.0*(nx + ny)*nx*ny*8.0/(1024.0**3)
    if gb > CACHE_BOUNDARY_MAX_GB:
        print("   boundary matrix would be %.1f GB > %.1f GB cap: "
              "using the shipped loop" % (gb, CACHE_BOUNDARY_MAX_GB),
              flush=True)
        eq._meq_boundary_matrix = None
        return None

    from freegs4e.gradshafranov import Greens

    dR = R[1, 0] - R[0, 0]
    dZ = Z[0, 1] - Z[0, 0]
    weight = np.outer(_romb_weights(nx), _romb_weights(ny))*(dR*dZ)

    idx = _boundary_indices(nx, ny)
    started = time.time()
    M = np.empty((idx.shape[0], nx*ny), dtype=float)
    for row, (x, y) in enumerate(idx):
        g = Greens(R, Z, R[x, y], Z[x, y])
        g[x, y] = 0.0                      # the self term, as the loop zeroes it
        M[row, :] = (g*weight).ravel()

    # A DELIBERATE PERTURBATION, for measuring how far a round-off change in
    # the boundary condition moves the converged state. Off unless asked.
    probe = os.environ.get("FGSREF_M_PERTURB")
    if probe:
        M *= (1.0 + float(probe))
        print("   boundary matrix perturbed by %s relative" % probe, flush=True)

    print("   boundary matrix cached: %d x %d, %.2f GB, built in %.1f s"
          % (M.shape[0], M.shape[1], gb, time.time() - started), flush=True)
    eq._meq_boundary_matrix = (M, idx)
    return eq._meq_boundary_matrix


def _free_boundary_cached(eq, Jtor, psi):
    """boundary.freeBoundary, reassociated. Same arithmetic, one matvec."""
    have = _boundary_matrix(eq)
    if have is None:
        from freegs4e import boundary as _bnd
        return _bnd.freeBoundary(eq, Jtor, psi)
    M, idx = have
    values = M @ np.asarray(Jtor, dtype=float).ravel()
    psi[idx[:, 0], idx[:, 1]] = values


# ----------------------------------------------------------------------------
# WORKAROUND 5
# boundary.freeBoundaryHagenow displaces its observation point off the boundary
# by a HARD-CODED eps = 1e-2 METRES, "to avoid the singularity in G(R,R') when
# R'=R".  That is a fixed PHYSICAL distance, so it does not shrink with the
# grid, and the von Hagenow boundary condition is therefore INCONSISTENT: it
# converges, but to a different answer from the direct Green's integral.
#
# THAT IS WHY THE TWO BOUNDARY CONDITIONS DISAGREE BY A GAP THAT DOES NOT
# SHRINK.  README.md records 6.460e-04 at 129^2 and 6.404e-04 at 257^2 in
# psi_ax and reads it as an open question about which one to believe.  It is
# not: measured directly, on ONE fixed Jtor with no solve in the way -- so the
# direct integral is ground truth and what is printed is von Hagenow's own
# error -- the observed order in h is
#
#     eps = 1e-2 fixed :  0.63, 0.21, 0.07     <- stalling, i.e. inconsistent
#     eps = 0.2 h      :  0.90, 1.00, 0.99     <- clean first order
#     eps = 0.5 h      :  0.97, 0.99, 0.99     <- the same, 1.8x larger
#
# at n = 65, 129, 257 on a 0.1..2.0 by -1..1 box.
#
# AND THE OPTIMUM IS NOT "AS SMALL AS POSSIBLE", which is why the constant is
# 0.2 rather than something tiny.  Swept at a FIXED grid of 129^2 the gap reads
# 7.6e-03, 3.3e-03, 8.9e-03, 1.6e-02, 2.3e-02 at eps = 1e-2, 3e-3, 1e-3, 3e-4,
# 1e-4: below about 0.2 h the log spike is NARROWER THAN THE CELL and the
# Romberg rule misses it, above it the displacement error dominates.  Tying eps
# to h is what keeps both terms on the same footing.
#
# WHAT THIS DOES NOT DO IS MAKE IT HIGH ORDER.  First order is the ceiling for
# a displaced observation point, and the residual is the QUADRATURE of the
# near-singular kernel rather than the displacement itself -- established by
# trying the obvious cure and measuring that it does not work.  oint G sigma dl
# is a SINGLE-LAYER potential, continuous across the boundary with its normal
# derivative jumping, so averaging the evaluations at +eps and -eps ought to
# cancel the O( eps ) term exactly.  Measured, it does not: 1.5e-02, 8.1e-03,
# 4.0e-03, 2.0e-03 against the one-sided 1.2e-02, 6.5e-03, 3.3e-03, 1.6e-03 --
# the same rate of 1.00 and slightly WORSE.  Getting past first order needs the
# log singularity subtracted and integrated in closed form, which is a real
# boundary-element method and is not this.
#
# Transcribed from ../freegs4e/freegs4e/boundary.py with the one constant
# lifted out and nothing else changed.  freegs4e is not MEQ's to edit.
# ----------------------------------------------------------------------------
def _hagenow_scaled_eps(eq, Jtor, psi):
    from freegs4e.gradshafranov import Greens
    from scipy.integrate import romb

    R = eq.R
    Z = eq.Z
    nx, ny = psi.shape
    dR = R[1, 0] - R[0, 0]
    dZ = Z[0, 1] - Z[0, 0]

    rhs = eq.R*Jtor
    rhs[0, :] = 0.0
    rhs[:, 0] = 0.0
    rhs[-1, :] = 0.0
    rhs[:, -1] = 0.0
    psi_fixed = eq.callSolver(psi, rhs)

    coeffs = [(0, 25.0/12), (1, -4.0), (2, 3.0), (3, -16.0/12), (4, 1.0/4)]
    dUdn_L = sum([w*psi_fixed[i, :] for i, w in coeffs])/dR
    dUdn_R = sum([w*psi_fixed[-(1 + i), :] for i, w in coeffs])/dR
    dUdn_D = sum([w*psi_fixed[:, i] for i, w in coeffs])/dZ
    dUdn_U = sum([w*psi_fixed[:, -(1 + i)] for i, w in coeffs])/dZ

    dd = np.sqrt(dR**2 + dZ**2)
    dUdn_L[0] = dUdn_D[0] = sum(
        [w*psi_fixed[i, i] for i, w in coeffs])/dd
    dUdn_L[-1] = dUdn_U[0] = sum(
        [w*psi_fixed[i, -(1 + i)] for i, w in coeffs])/dd
    dUdn_R[0] = dUdn_D[-1] = sum(
        [w*psi_fixed[-(1 + i), i] for i, w in coeffs])/dd
    dUdn_R[-1] = dUdn_U[-1] = sum(
        [w*psi_fixed[-(1 + i), -(1 + i)] for i, w in coeffs])/dd

    # THE ONE CHANGED LINE. eps was 1e-2, a length in metres.
    eps = HAGENOW_EPS*min(abs(dR), abs(dZ))

    idx = np.concatenate([
        [(x, 0, 0.0, -eps) for x in range(nx)],
        [(x, ny - 1, 0.0, eps) for x in range(nx)],
        [(0, y, -eps, 0.0) for y in range(ny)],
        [(nx - 1, y, eps, 0.0) for y in range(ny)],
    ])

    for x, y, Reps, Zeps in idx:
        x = int(round(x))
        y = int(round(y))
        Rpos = R[x, y] + Reps
        Zpos = Z[x, y] + Zeps

        result = romb(Greens(R[0, :], Z[0, :], Rpos, Zpos)
                      * dUdn_L/R[0, :])*dZ
        result += romb(Greens(R[-1, :], Z[-1, :], Rpos, Zpos)
                       * dUdn_R/R[-1, :])*dZ
        result += romb(Greens(R[:, 0], Z[:, 0], Rpos, Zpos)
                       * dUdn_D/R[:, 0])*dR
        result += romb(Greens(R[:, -1], Z[:, -1], Rpos, Zpos)
                       * dUdn_U/R[:, -1])*dR
        psi[x, y] = result


# WORKAROUND 4
# machine.DIIID() builds its coil list as a list of DICTS, while
# Machine.__init__ -> getCurrentsVec() iterates it as (label, coil) pairs, so
# freegs4e.machine.DIIID() raises at construction:
#     ValueError: too many values to unpack (expected 2, got 4)
# The coil table itself is fine; only its container is wrong.  Rebuilt here
# from the same numbers, as (label, Coil) tuples.
# ----------------------------------------------------------------------------
_DIIID_COILS = [
    ("F1A", 0.8608, 0.16830), ("F2A", 0.8614, 0.50810),
    ("F3A", 0.8628, 0.84910), ("F4A", 0.8611, 1.1899),
    ("F5A", 1.0041, 1.5169), ("F6A", 2.6124, 0.4376),
    ("F7A", 2.3733, 1.1171), ("F8A", 1.2518, 1.6019),
    ("F9A", 1.6890, 1.5874), ("F1B", 0.8608, -0.1737),
    ("F2B", 0.8607, -0.5135), ("F3B", 0.8611, -0.8543),
    ("F4B", 0.8630, -1.1957), ("F5B", 1.0025, -1.5169),
    ("F6B", 2.6124, -0.44376), ("F7B", 2.3834, -1.1171),
    ("F8B", 1.2524, -1.6027), ("F9B", 1.6889, -1.578),
]


def DIIID_asCoils():
    return machine.Machine([(lab, machine.Coil(R, Z))
                            for lab, R, Z in _DIIID_COILS])


# ----------------------------------------------------------------------------
# configurations
# ----------------------------------------------------------------------------
def pshape(x, alpha, beta):
    return (1.0 - np.clip(x, 0.0, 1.0) ** alpha) ** beta


def split_amplitudes(frac_p, R0, P0=1.0):
    """Choose (P0, F0) so that at psi_n = 0 and R = R0 a fraction frac_p of
    Jtor = R p' + ff'/(mu0 R) comes from the p' term.

    frac_p > 1 gives ff' < 0 (a diamagnetic plasma).  The overall size of the
    pair is irrelevant: Ip_logic rescales both by a common factor, so only the
    ratio is a physical choice.
    """
    F0 = mu0 * R0 * R0 * P0 * (1.0 - frac_p) / frac_p
    return P0, F0


# ----------------------------------------------------------------------------
# A LIMITED CIRCULAR TOKAMAK, for the comparison MEQ can actually make.
#
# Every other machine in this file is DIVERTED, and MEQ's moving plasma support
# is a pointwise test on the normalised flux, which picks up the private flux
# region beyond an X-point -- FREE-BOUNDARY-PLAN.md section 10.3.  So a
# like-for-like free-boundary comparison needs a limited plasma, and no machine
# in freegs4e is one.
#
# VERTICAL FIELD ONLY.  Two coil pairs, placed to give a field that holds the
# plasma at R0 without producing a null anywhere in the box: with no divertor
# coils there is no X-point for the boundary logic to find, and the plasma can
# only be bounded by the limiter.
R0_LIM, A_LIM = 1.00, 0.35
_th = np.linspace(0.0, 2.0*np.pi, 257)
LIMITER = (R0_LIM + A_LIM*np.cos(_th), A_LIM*np.sin(_th))


LIMITED_COILS = [("P1U", 1.75, 0.90), ("P1L", 1.75, -0.90),
                 ("P2U", 0.55, 1.10), ("P2L", 0.55, -1.10)]

# MEQ's conductors are RECTANGLES of half-width 0.05, so this is their extent.
# examples/limited-tokamak.toml carries the same four numbers.
LIMITED_COIL_HALFWIDTH = 0.05
LIMITED_COIL_HALFHEIGHT = 0.05


def make_limited():
    from freegs4e.machine import Coil, Machine, Wall
    coils = [(lab, Coil(R, Z)) for lab, R, Z in LIMITED_COILS]
    return Machine(coils, Wall(LIMITER[0], LIMITER[1]))


# ----------------------------------------------------------------------------
# THE SAME MACHINE WITH THE CONDUCTORS GIVEN THEIR EXTENT, WHICH IS THE ONE
# MODELLING DIFFERENCE LEFT BETWEEN THIS REFERENCE AND MEQ.
#
# freegs4e's default Coil is an exact FILAMENT: controlPsi is
# Greens( self.R, self.Z, R, Z )*turns, a point source, and its `area`
# attribute only imposes a current-density limit and never enters the field.
# MEQ's four conductors are 0.1 x 0.1 m rectangles carrying a uniform current
# density -- meq::Coil::currentDensity is I/area and the source is assembled by
# quadrature over the meshed rectangle.  So FB-6's agreement was reached DESPITE
# a modelling difference rather than because the models agree, and the leading
# finite-size correction goes as ( w/d )^2, about 1.6e-02 on the near field with
# the coils 0.4 m from the plasma edge -- three orders above the 1.3e-04 quoted.
#
# freegs4e.shaped_coil.ShapedCoil is freegs4e's OWN answer to this and is what
# closes it: a polygon, triangulated, with a Gauss rule per triangle whose
# weights sum to one, so controlPsi is the AVERAGE Greens over the cross-section
# and the total current is spread uniformly across it.  That is MEQ's model term
# for term, in the reference's own code rather than in a MEQ-side correction.
#
# THE RULE IS 12 POINTS AND THAT IS ENOUGH HERE, MEASURED RATHER THAN ASSUMED.
# polygon_quad splits the square into two triangles and puts freegs4e's 6-point
# degree-4 rule on each.  Against a 64 x 64 tensor Gauss-Legendre rule over the
# same square, the worst relative error in controlPsi over this case's grid --
# its whole boundary ring plus the grid point nearest each coil -- is 3.0e-05,
# attained at ( 1.7500, 0.8000 ), the boundary point two half-heights below P1U.
# It is the same at 129^2 and 513^2, being a property of the rule and the
# geometry rather than of the grid.  Everywhere the plasma is it is far smaller.
def make_limited_shaped(hw=LIMITED_COIL_HALFWIDTH, hh=LIMITED_COIL_HALFHEIGHT):
    from freegs4e.machine import Machine, Wall
    from freegs4e.shaped_coil import ShapedCoil

    def square(Rc, Zc):
        return [(Rc - hw, Zc - hh), (Rc + hw, Zc - hh),
                (Rc + hw, Zc + hh), (Rc - hw, Zc + hh)]

    coils = [(lab, ShapedCoil(square(R, Z))) for lab, R, Z in LIMITED_COILS]
    return Machine(coils, Wall(LIMITER[0], LIMITER[1]))


CASES = [
    dict(
        name="A_testtokamak_classic",
        machine="TestTokamak", make=machine.TestTokamak,
        grid=dict(Rmin=0.1, Rmax=2.0, Zmin=-1.0, Zmax=1.0, nx=129, ny=129),
        order=4, Ip=2.0e5, fvac=2.0,
        R0=1.35, frac_p=0.50, pa=1.5, pb=2.0, fa=1.2, fb=2.0,
        xpoints=[(1.1, -0.6), (1.1, 0.8)],
        isoflux=[(1.1, -0.6, 1.1, 0.8)],
        limiter=None, seed_paxis=1.0e3,
        notes=("The FreeGS worked example's geometry: up-down asymmetric "
               "double null.  Half the on-axis current density comes from p' "
               "and half from ff'."),
    ),
    dict(
        name="B_testtokamak_peaked_ffdom",
        machine="TestTokamak", make=machine.TestTokamak,
        grid=dict(Rmin=0.1, Rmax=2.0, Zmin=-1.0, Zmax=1.0, nx=129, ny=129),
        order=4, Ip=2.5e5, fvac=2.5,
        R0=1.25, frac_p=0.35, pa=2.0, pb=3.0, fa=1.0, fb=1.0,
        xpoints=[(1.1, -0.6), (1.1, 0.8)],
        isoflux=[(1.1, -0.6, 1.1, 0.8)],
        limiter=None, seed_paxis=1.0e3,
        notes=("Same machine and shape control as A, so A/B/E form a "
               "controlled scan over the source alone.  Here p' is strongly "
               "peaked (exponents 2, 3) and supplies only 35% of the on-axis "
               "current, and ff' is LINEAR in psi_n: ff'-dominated, and the "
               "two profile shapes are as different from each other as the "
               "tabulation allows.  A symmetric double null was tried first "
               "and limit-cycles: see the report."),
    ),
    dict(
        name="C_mast_spherical",
        machine="MAST_sym", make=machine.MAST_sym,
        grid=dict(Rmin=0.1, Rmax=2.0, Zmin=-2.0, Zmax=2.0, nx=129, ny=129),
        order=4, Ip=5.0e5, fvac=0.5,
        R0=0.9, frac_p=0.55, pa=1.5, pb=2.0, fa=1.5, fb=2.0,
        xpoints=[(0.7, -1.1), (0.7, 1.1)],
        isoflux=[(0.7, -1.1, 1.45, 0.0), (0.7, 1.1, 1.45, 0.0),
                 (0.7, -1.1, 0.7, 1.1)],
        limiter=None, seed_paxis=1.0e3,
        notes=("Spherical tokamak (MAST, up-down symmetric circuits): aspect "
               "ratio about 1.5, elongation above 2, low fvac."),
    ),
    dict(
        name="D_tcv_conventional",
        machine="TCV", make=machine.TCV,
        grid=dict(Rmin=0.4, Rmax=1.4, Zmin=-0.8, Zmax=0.8, nx=129, ny=129),
        order=4, Ip=1.2e5, fvac=1.4,
        R0=0.89, frac_p=0.40, pa=1.5, pb=2.0, fa=1.5, fb=2.0,
        xpoints=[(0.75, -0.55), (0.75, 0.55)],
        isoflux=[(0.75, -0.55, 0.75, 0.55), (0.68, 0.0, 1.12, 0.0)],
        limiter=None, seed_paxis=1.0e3,
        notes=("Conventional aspect ratio, strongly elongated, twenty small "
               "shaping coils."),
    ),
    dict(
        name="E_testtokamak_diamagnetic",
        machine="TestTokamak", make=machine.TestTokamak,
        grid=dict(Rmin=0.1, Rmax=2.0, Zmin=-1.0, Zmax=1.0, nx=129, ny=129),
        order=4, Ip=2.0e5, fvac=0.6,
        R0=1.35, frac_p=1.15, pa=1.2, pb=1.5, fa=1.5, fb=2.0,
        xpoints=[(1.1, -0.6), (1.1, 0.6)],
        isoflux=[(1.1, -0.6, 1.1, 0.6)],
        limiter=None, seed_paxis=1.0e3,
        notes=("Same machine and Ip as A, symmetric double null, but "
               "frac_p = 1.15 so ff' is "
               "NEGATIVE -- a diamagnetic, pressure-driven plasma in which "
               "the two terms of Jtor oppose each other, and 115% of the "
               "on-axis current comes from p'.  fvac is 0.6 rather than 2.0; "
               "note that in this formulation fvac does not enter the GS "
               "solve at all, only f(psi_n) (see the report)."),
    ),
    dict(
        name="F_diiid_conventional",
        machine="DIIID (coil table rebuilt, see WORKAROUND 4)",
        make=DIIID_asCoils,
        grid=dict(Rmin=0.8, Rmax=2.6, Zmin=-1.4, Zmax=1.4, nx=129, ny=129),
        order=4, Ip=1.0e6, fvac=3.5,
        R0=1.7, frac_p=0.45, pa=1.5, pb=2.0, fa=1.5, fb=2.0,
        xpoints=[(1.2, -1.0), (1.2, 1.0)],
        isoflux=[(1.2, -1.0, 1.2, 1.0), (1.05, 0.0, 2.35, 0.0)],
        limiter=None, seed_paxis=1.0e4,
        notes=("Conventional large tokamak, MA-class plasma current, high "
               "fvac."),
    ),
    dict(
        name="G_mastu_simple",
        machine="MASTU_simple", make=machine.MASTU_simple,
        grid=dict(Rmin=0.1, Rmax=2.0, Zmin=-2.2, Zmax=2.2, nx=129, ny=129),
        order=4, Ip=6.0e5, fvac=0.6,
        R0=0.85, frac_p=0.60, pa=2.0, pb=2.0, fa=1.2, fb=1.5,
        xpoints=[(0.62, -1.13), (0.62, 1.13)],
        isoflux=[(0.62, -1.13, 0.62, 1.13), (0.36, 0.0, 1.38, 0.0)],
        limiter=None, seed_paxis=1.0e3,
        notes=("MAST-U (simplified single-strand coil set): spherical "
               "tokamak, different coil topology from C."),
    ),
    dict(
        name="H_limited_circular",
        machine="LimitedCircular", make=lambda: make_limited(),
        grid=dict(Rmin=0.30, Rmax=1.90, Zmin=-0.80, Zmax=0.80, nx=129, ny=129),
        order=4, Ip=3.0e5, fvac=1.0,
        R0=1.00, frac_p=0.50, pa=1.0, pb=2.0, fa=1.0, fb=2.0,
        xpoints=[],
        isoflux=[(0.65, 0.0, 1.35, 0.0), (1.00, 0.35, 1.35, 0.0),
                 (1.00, -0.35, 1.35, 0.0)],
        limiter=LIMITER, seed_paxis=1.0e3,
        notes=("A LIMITED plasma, and the only one here: vertical-field coils "
               "only, so no X-point exists in range and the boundary is the "
               "flux surface through the limiter contact rather than a "
               "separatrix.  It is the case MEQ can actually compare against "
               "-- FREE-BOUNDARY-PLAN.md section 10 records that MEQ's plasma "
               "support test is pointwise and so is limiter-only, and that "
               "every other case in this table is diverted."),
    ),
    dict(
        name="I_limited_shaped",
        machine="LimitedCircular (ShapedCoil, 0.1 x 0.1 m)",
        make=lambda: make_limited_shaped(),
        grid=dict(Rmin=0.30, Rmax=1.90, Zmin=-0.80, Zmax=0.80, nx=129, ny=129),
        order=4, Ip=3.0e5, fvac=1.0,
        R0=1.00, frac_p=0.50, pa=1.0, pb=2.0, fa=1.0, fb=2.0,
        xpoints=[],
        isoflux=[(0.65, 0.0, 1.35, 0.0), (1.00, 0.35, 1.35, 0.0),
                 (1.00, -0.35, 1.35, 0.0)],
        limiter=LIMITER, seed_paxis=1.0e3,
        notes=("H WITH THE CONDUCTORS GIVEN THEIR EXTENT, and identical to it "
               "in every other field -- same grid, same order, same Ip, same "
               "profile shape, same isoflux constraints, same wall.  It exists "
               "so that the ONE remaining modelling difference between this "
               "reference and MEQ can be measured rather than estimated: H's "
               "coils are exact filaments and MEQ's are 0.1 x 0.1 m rectangles "
               "carrying a uniform current density, which is what ShapedCoil "
               "makes freegs4e do too.  Compare I against H at the same grid "
               "and the difference is the conductor model alone."),
    ),
]


# ----------------------------------------------------------------------------
def build_eq(case, tok):
    eq = equilibrium.Equilibrium(tokamak=tok, order=case["order"],
                                 **case["grid"])
    # WORKAROUND 3: __init__ has the _updatePlasmaPsi call commented out, so
    # psi_func / mask / psi_axis do not exist and eq.Br() raises.
    eq.mask_outside_limiter = None
    eq.mask_inside_limiter = None
    eq._updatePlasmaPsi(eq.plasma_psi)
    return eq


def seed_from_coarse(eq, case, say):
    """Start Picard from a converged COARSER grid, using the nesting.

    THE GRIDS ARE 2^n + 1 SO THAT THEY NEST, AND UNTIL 2026-09-06 NOTHING USED
    IT.  Every --nx run began from whatever Equilibrium.__init__ left in psi --
    a Gaussian bump -- so a fine grid paid the full cold Picard count, 23 to 103
    steps, at its own O( n^3 ) per-step boundary cost.  The boundary condition
    alone scales about 6.9x per grid doubling, so that is the difference between
    a handful of steps and an afternoon.

    THE SEED CANNOT MOVE THE ANSWER, WHICH IS THE WHOLE REASON IT IS SAFE.
    Picard converges to the FINE grid's own solution whatever it starts from, so
    an interpolated coarse answer changes the number of iterations and nothing
    else.  That is a claim to check rather than assume: run one case both ways
    and difference psi_axis.  The same control is already recorded for the
    ConstrainPaxisIp seed in run_case ("cold and seeded give bit-identical
    answers on A").

    A CUBIC SPLINE, NOT A BILINEAR LIFT, and not because of accuracy -- a seed
    does not need accuracy.  It is freegs4e's own representation of psi
    (Equilibrium._updatePlasmaPsi builds exactly this), so the seed is smooth in
    the same sense the solver's own iterate is, and a kinked seed would spend
    its first iterations being smoothed rather than being converged.

    THE NESTING IS ASSERTED AND NOT ASSUMED.  A coarse file from a different
    Rmin/Rmax, or at a resolution that does not divide, is a silent wrong answer
    of exactly the kind this benchmark exists to catch: it would interpolate,
    converge, and describe the machine it was given rather than the one asked
    for.  Both are refused.
    """
    path = os.path.join(SEED_FROM, case["name"] + ".npz")
    if not os.path.exists(path):
        say("no coarse seed at %s -- starting cold" % path)
        return None

    with np.load(path, allow_pickle=True) as d:
        Rc = np.asarray(d["R"], float)
        Zc = np.asarray(d["Z"], float)
        psi_c = np.asarray(d["plasma_psi"], float)

    Rf, Zf = eq.R[:, 0], eq.Z[0, :]
    nc, nf = len(Rc), len(Rf)
    mc, mf = len(Zc), len(Zf)
    if nc > nf or mc > mf:
        raise ValueError("seed grid %dx%d is FINER than the run's %dx%d"
                         % (nc, mc, nf, mf))
    if (nf - 1) % (nc - 1) or (mf - 1) % (mc - 1):
        raise ValueError("seed grid %dx%d does not nest inside %dx%d: the "
                         "2^n + 1 convention is what makes the coarse points a "
                         "subset of the fine ones" % (nc, mc, nf, mf))
    sR, sZ = (nf - 1)//(nc - 1), (mf - 1)//(mc - 1)

    # The extents have to agree too.  Same count, different Rmin, still divides.
    off = max(np.abs(Rf[::sR] - Rc).max(), np.abs(Zf[::sZ] - Zc).max())
    scale = max(Rf[-1] - Rf[0], Zf[-1] - Zf[0])
    if off > 1e-12*scale:
        raise ValueError("seed grid points do not lie on this run's grid: "
                         "worst offset %.3e against an extent of %.3e -- the "
                         "coarse run used different Rmin/Rmax/Zmin/Zmax"
                         % (off, scale))

    spline = sinterp.RectBivariateSpline(Rc, Zc, psi_c, kx=3, ky=3, s=0)
    seed = spline(Rf, Zf)

    # A cubic spline reproduces its own data, so this is a check on the GRIDS
    # rather than on the interpolant: it fails if the strides above are wrong.
    worst = float(np.abs(seed[::sR, ::sZ] - psi_c).max())
    span = float(np.abs(psi_c).max()) or 1.0
    if worst > 1e-9*span:
        raise ValueError("the coarse points did not come back unchanged: "
                         "%.3e against a psi scale of %.3e" % (worst, span))

    eq._updatePlasmaPsi(seed)
    say("seeded from %s (%dx%d, stride %dx%d); coincident points reproduce "
        "to %.2e of |psi|max" % (path, nc, mc, sR, sZ, worst/span))
    return dict(path=path, nx=nc, ny=mc, stride=[sR, sZ],
                coincident_worst_rel=worst/span)


def get_wall(case, eq):
    """The machine wall, for reference only -- it is NOT used as a limiter
    unless case["limiter"] says so.  Returns empty arrays if there is none."""
    w = eq.tokamak.wall
    if w is None:
        return np.zeros(0), np.zeros(0)
    return np.asarray(w.R, float), np.asarray(w.Z, float)


def attach(prof, case, eq):
    """Optionally give the profile a limiter.  case["limiter"] is None for
    every case here: with Machine.current_vec kept in sync (WORKAROUND 1) the
    plain FreeGS boundary logic converges from a cold start on all of them,
    and a limiter would only add a wall-conforming boundary nobody wants."""
    lim = case.get("limiter")
    if lim is not None:
        prof.attach_limiter(eq.R, eq.Z, np.asarray(lim[0], float),
                            np.asarray(lim[1], float))
    return prof


def lcfs_from_rays(eq, opt, psi, psi_axis, psi_bndry, ntheta=512,
                   offset=0.0):
    """LCFS by bisection along rays from the magnetic axis.

    Assumes the plasma boundary is star-shaped about the axis, which is true
    for every case here (checked: the returned curve is single-valued in the
    geometric angle and closes on itself).
    """
    psin = (psi - psi_axis) / (psi_bndry - psi_axis)
    f = sinterp.RectBivariateSpline(eq.R[:, 0], eq.Z[0, :], psin)
    r0, z0 = float(opt[0][0]), float(opt[0][1])
    Rlo, Rhi = eq.R[0, 0], eq.R[-1, 0]
    Zlo, Zhi = eq.Z[0, 0], eq.Z[0, -1]
    ds = 0.25 * min(eq.dR, eq.dZ)
    out, missed = [], 0
    # A ray aimed exactly down a separatrix leg can have psi_n stay just below
    # 1 all the way out of the grid; offsetting the whole theta grid by half a
    # step moves it off the singular direction.  freegs4e's own
    # critical.find_separatrix does the same thing for the same reason.
    for t in np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False) + offset:
        dr, dz = np.cos(t), np.sin(t)
        s, last_in, hi, ok = 0.0, 0.0, None, False
        while True:
            s += ds
            r, z = r0 + s * dr, z0 + s * dz
            if not (Rlo < r < Rhi and Zlo < z < Zhi):
                break
            if float(f(r, z, grid=False)) >= 1.0:
                hi, ok = s, True
                break
            last_in = s
        if not ok:
            missed += 1
            continue
        a, b = last_in, hi
        for _ in range(90):
            m = 0.5 * (a + b)
            if float(f(r0 + m * dr, z0 + m * dz, grid=False)) < 1.0:
                a = m
            else:
                b = m
        m = 0.5 * (a + b)
        out.append((r0 + m * dr, z0 + m * dz))
    return np.array(out), missed


def gs_residual(eq, psi, Jtor, mask):
    res = {}
    rhs = -mu0 * eq.R * Jtor
    core = (mask > 0.5) if mask is not None else np.ones_like(psi, bool)
    inner = np.zeros_like(psi, dtype=bool)
    inner[2:-2, 2:-2] = True
    m = core & inner
    scale = float(np.max(np.abs(rhs[m])))
    # (a) freegs4e's 2nd-order Delta*, applied to the TOTAL psi.  Inside the
    # core there are no coils, so Delta* psi_coil = 0 there and the total psi
    # must satisfy the same equation as the plasma part.  This is an
    # independent, discretisation-level check.
    ds = GSElliptic(eq.Rmin)(psi, eq.dR, eq.dZ)
    res["gs_resid_core_linf_2nd"] = float(np.max(np.abs((ds - rhs)[m])) / scale)
    res["gs_resid_core_l2_2nd"] = float(
        np.sqrt(np.mean(((ds - rhs)[m]) ** 2)) / scale)
    # (b) the 4th-order operator freegs4e actually inverted, applied to the
    # plasma part.  This is ~machine precision iff the converged Jtor is the
    # same one the last linear solve used, so it is a convergence check.
    A = GSsparse4thOrder(eq.Rmin, eq.Rmax, eq.Zmin, eq.Zmax)(eq.nx, eq.ny)
    b = A.dot(eq.plasma_psi.ravel()).reshape(eq.nx, eq.ny)
    res["gs_resid_core_linf_4th"] = float(np.max(np.abs((b - rhs)[m])) / scale)
    res["gs_resid_scale"] = scale
    return res


def shape_params(lR, lZ):
    Rmin, Rmax = float(lR.min()), float(lR.max())
    Zmin, Zmax = float(lZ.min()), float(lZ.max())
    R0 = 0.5 * (Rmin + Rmax)
    a = 0.5 * (Rmax - Rmin)
    return dict(
        lcfs_Rmin=Rmin, lcfs_Rmax=Rmax, lcfs_Zmin=Zmin, lcfs_Zmax=Zmax,
        R0=float(R0), a=float(a), aspect_ratio=float(R0 / a),
        kappa=float((Zmax - Zmin) / (2.0 * a)),
        delta_upper=float((R0 - lR[np.argmax(lZ)]) / a),
        delta_lower=float((R0 - lR[np.argmin(lZ)]) / a),
        delta=float(0.5 * ((R0 - lR[np.argmax(lZ)]) / a
                           + (R0 - lR[np.argmin(lZ)]) / a)),
    )


# ----------------------------------------------------------------------------
def run_case(case):
    t0 = time.time()
    rec = dict(name=case["name"], machine=case["machine"], notes=case["notes"])
    log = []

    def say(*a):
        s = " ".join(str(x) for x in a)
        print("   " + s, flush=True)
        log.append(s)

    tok = case["make"]()
    eq = build_eq(case, tok)
    if HAGENOW:
        from freegs4e import boundary as _bnd
        eq._applyBoundary = (_bnd.freeBoundaryHagenow if HAGENOW_EPS is None
                             else _hagenow_scaled_eps)
        print("   von Hagenow boundary condition installed%s"
              % ("" if HAGENOW_EPS is None
                 else "  ( eps = %.3g h, WORKAROUND 5 )" % HAGENOW_EPS),
              flush=True)
    elif CACHE_BOUNDARY:
        # THE DEFAULT BOUNDARY CONDITION, REASSOCIATED. Not an approximation
        # and not an alternative method: the same Greens, the same Romberg
        # rule, the same self term zeroed, contracted in a different order.
        # So it is on unless asked otherwise, and von Hagenow -- which IS a
        # different method and a different answer -- takes precedence when it
        # is asked for.
        eq._applyBoundary = _free_boundary_cached
        print("   default boundary condition cached as a matrix "
              "( WORKAROUND 6 )", flush=True)
    if VCYCLE_LEVELS:
        install_vcycle(eq, int(case["order"]), VCYCLE_LEVELS)
        print("   V-cycle installed: %d levels on the order-%d generator"
              % (VCYCLE_LEVELS, int(case["order"])), flush=True)
    ctrl = SyncConstrain(xpoints=case["xpoints"], isoflux=case["isoflux"],
                         gamma=1e-12)
    rec["seed_from_coarse"] = (seed_from_coarse(eq, case, say)
                               if SEED_FROM else None)

    psi_n = np.linspace(0.0, 1.0, NPSI)
    P0, F0 = split_amplitudes(case["frac_p"], case["R0"])
    pprime_tab = P0 * pshape(psi_n, case["pa"], case["pb"])
    ffprime_tab = F0 * pshape(psi_n, case["fa"], case["fb"])
    rec["frac_p_on_axis"] = float(case["frac_p"])
    rec["profile_shape"] = dict(pprime="P0*(1-psi_n**%g)**%g"
                                % (case["pa"], case["pb"]),
                                ffprime="F0*(1-psi_n**%g)**%g"
                                % (case["fa"], case["fb"]))

    # ---- stage 1+2: the tabulated profile from a cold start, Ip_logic on,
    # to discover the scale factor Ip_logic applies.  If that fails, retry
    # after seeding with ConstrainPaxisIp, whose Ip/paxis normalisation is
    # more forgiving from the default Gaussian guess (measured: not needed by
    # any case here -- cold and seeded give bit-identical answers on A).
    def make_prof():
        p = GeneralProfile(Ip=case["Ip"], fvac=case["fvac"], psi_n=psi_n,
                           pprime_data=pprime_tab, ffprime_data=ffprime_tab,
                           Raxis=1.0, Ip_logic=True)
        return attach(p, case, eq)

    rec["seed_used"] = False
    rec["seed_status"] = None
    rec["blend_used"] = 0.0
    prof = make_prof()
    s1 = picard_loop(eq, prof, ctrl, rtol=RTOL, atol=1e-14, maxits=MAXITS)
    if not good_enough(s1):
        say("cold start hit maxits at rel %.2e; retrying with blend=0.4"
            % s1["rel_history"][-1])
        rec["blend_used"] = 0.4
        s1 = picard_loop(eq, prof, ctrl, rtol=RTOL, atol=1e-14,
                         maxits=MAXITS, blend=0.4)
    if not good_enough(s1):
        say("still maxits at rel %.2e; restarting from a ConstrainPaxisIp "
            "seed" % s1["rel_history"][-1])
        eq = build_eq(case, tok)
        ctrl = SyncConstrain(xpoints=case["xpoints"],
                             isoflux=case["isoflux"], gamma=1e-12)
        # build_eq() returns a COLD equilibrium, so re-seed or this branch
        # silently throws the coarse answer away.
        if SEED_FROM:
            rec["seed_from_coarse"] = seed_from_coarse(eq, case, say)
        seed = PaxisProfile(case["seed_paxis"], case["Ip"], case["fvac"],
                            alpha_m=1.0, alpha_n=2.0, Raxis=1.0)
        attach(seed, case, eq)
        s0 = picard_loop(eq, seed, ctrl, rtol=1e-6, maxits=200)
        say("seed (ConstrainPaxisIp): %s after %d its"
            % (s0["status"], s0["iterations"]))
        rec["seed_used"] = True
        rec["seed_status"] = s0["status"]
        prof = make_prof()
        s1 = picard_loop(eq, prof, ctrl, rtol=RTOL, atol=1e-14, maxits=MAXITS,
                         blend=rec["blend_used"])
        say("GeneralPprimeFFprime (Ip_logic=True): %s after %d its, rel %.3e"
        % (s1["status"], s1["iterations"], s1["rel_history"][-1]))
    L1 = float(prof.L)
    say("Ip_logic scale factor L = %.12g" % L1)
    rec["stage2_status"] = s1["status"]
    rec["stage2_iterations"] = s1["iterations"]
    rec["stage2_rel"] = float(s1["rel_history"][-1])
    rec["Ip_logic_L"] = L1

    # ---- stage 3: fold L into the table and re-solve with Ip_logic=False, so
    # the saved arrays ARE the source with no hidden rescaling.
    psi_before = eq.psi().copy()
    Lk, passes = L1, []
    for _ in range(5):
        prof2 = GeneralProfile(Ip=case["Ip"], fvac=case["fvac"], psi_n=psi_n,
                               pprime_data=Lk * pprime_tab,
                               ffprime_data=Lk * ffprime_tab,
                               Raxis=1.0, Ip_logic=False)
        attach(prof2, case, eq)
        s2 = picard_loop(eq, prof2, ctrl, rtol=RTOL, atol=1e-14, maxits=MAXITS)
        Ip_now = float(np.sum(prof2.jtor) * eq.dR * eq.dZ)
        err = abs(Ip_now / case["Ip"] - 1.0)
        passes.append((float(Lk), Ip_now, float(err), s2["status"],
                       s2["iterations"]))
        say("  frozen-L pass: L=%.10g -> Ip=%.8g (%.2e off), %s in %d its"
            % (Lk, Ip_now, err, s2["status"], s2["iterations"]))
        if err < 1e-8:
            break
        # Ip_logic would have applied exactly this factor; do it explicitly so
        # the saved table needs no further scaling.
        Lk = Lk * case["Ip"] / Ip_now
    L1 = float(Lk)
    rec["frozen_L_passes"] = passes
    drift = float(np.max(np.abs(eq.psi() - psi_before))
                  / (np.max(psi_before) - np.min(psi_before)))
    say("frozen-L re-solve (Ip_logic=False): %s after %d its, rel %.3e; "
        "psi moved %.3e from the Ip_logic solution"
        % (s2["status"], s2["iterations"], s2["rel_history"][-1], drift))
    rec["stage3_status"] = s2["status"]
    rec["stage3_iterations"] = s2["iterations"]
    rec["stage3_rel"] = float(s2["rel_history"][-1])
    rec["psi_drift_on_freezing_L"] = drift
    rec["picard_status"] = s2["status"]
    rec["converged"] = bool(s2["status"] in ("rtol", "atol")
                            or s2["rel_history"][-1] < 1e-5)
    rec["Ip_logic_for_saved_profile"] = False
    prof = prof2

    # ---- final consistent evaluation on the converged psi
    psi = eq.psi()
    opt, xpt, psi_bndry, mask, kind = prof.boundary_and_mask(eq.R, eq.Z, psi)
    psi_axis = float(opt[0][2])
    Jtor = prof.Jtor_part2(eq.R, eq.Z, psi, psi_axis, psi_bndry, mask)
    Ip_ach = float(np.sum(Jtor) * eq.dR * eq.dZ)
    say("boundary: %s   psi_axis %.10g   psi_bndry %.10g   core cells %d"
        % (kind, psi_axis, psi_bndry, int(mask.sum())))
    say("Ip achieved %.8g A (target %.8g), L on this pass = %.12g"
        % (Ip_ach, case["Ip"], prof.L))
    rec["boundary_kind"] = kind
    rec["L_on_final_pass"] = float(prof.L)

    eq.opt, eq.xpt = opt, xpt
    eq.psi_axis, eq.psi_bndry, eq.mask = psi_axis, psi_bndry, mask

    # ---- the profiles as actually used
    pprime_used = np.asarray(prof.pprime(psi_n), float)
    ffprime_used = np.asarray(prof.ffprime(psi_n), float)

    # Decisive check on the UNITS of the saved arrays: rebuild Jtor through
    # Jtor = R p' + ff'/(mu0 R), the form the GS equation demands, with p'
    # read as dp/dpsi (NOT dp/dpsi_n).  Agreement to interpolation error is
    # what says the arrays are dp/dpsi.
    psin_grid = np.clip((psi - psi_axis) / (psi_bndry - psi_axis), 0.0, 1.0)
    mk = mask if mask is not None else 1.0
    scaleJ = np.max(np.abs(Jtor))
    for tag, ip in (("linear", lambda y: np.interp(psin_grid, psi_n, y)),
                    ("cubic", lambda y: sinterp.CubicSpline(psi_n, y)(
                        psin_grid))):
        pp, ff = ip(pprime_used), ip(ffprime_used)
        J = (eq.R * pp + ff / (mu0 * eq.R)) * mk
        e = float(np.max(np.abs(J - Jtor)) / scaleJ)
        say("Jtor rebuilt as R*p' + ff'/(mu0 R) from the SAVED arrays, "
            "%-6s interpolation: rel err %.3e" % (tag, e))
        rec["jtor_rebuild_relerr_" + tag] = e
        if tag == "cubic":
            # the control: if the arrays were dp/dpsi_n they would have to be
            # divided by (psi_bndry - psi_axis) first.  Report how wrong that
            # would be.
            rec["jtor_rebuild_relerr_if_treated_as_dpsi_n"] = float(
                np.max(np.abs(J / (psi_bndry - psi_axis) - Jtor)) / scaleJ)
    rec["jtor_rebuild_relerr"] = rec["jtor_rebuild_relerr_cubic"]
    say("  (the same arrays read as dp/dpsi_n instead would be wrong by "
        "%.3g)" % rec["jtor_rebuild_relerr_if_treated_as_dpsi_n"])

    # How much did GeneralPprimeFFprime's UnivariateSpline SMOOTH the table it
    # was handed?  It is built with scipy's default s = len(psi_n) = 256, an
    # ABSOLUTE residual budget, so the answer depends on the amplitude of the
    # data: a profile of order 1e5 is fitted to about 4e-5, one of order 1 is
    # smoothed down to a single global cubic and is 1.3e-2 away from the
    # tabulation.  The arrays saved here are the SPLINE's values, i.e. what
    # the solver used, not the analytic table.
    rec["spline_smoothing_pprime"] = float(
        np.max(np.abs(pprime_used - L1 * pprime_tab))
        / max(np.max(np.abs(L1 * pprime_tab)), 1e-300))
    rec["spline_smoothing_ffprime"] = float(
        np.max(np.abs(ffprime_used - L1 * ffprime_tab))
        / max(np.max(np.abs(L1 * ffprime_tab)), 1e-300))
    say("UnivariateSpline(default s) moved the tabulation by %.3e (p') and "
        "%.3e (ff') relative" % (rec["spline_smoothing_pprime"],
                                 rec["spline_smoothing_ffprime"]))

    # ---- derived 1D profiles, from the saved p' and ff' exactly as the
    # base Profile class defines them (integrate in psi_n, then multiply by
    # dpsi/dpsi_n = psi_bndry - psi_axis; see the units discussion in the
    # report).  p(1) = 0 and f(1) = fvac are the boundary conditions.
    dpsi = psi_bndry - psi_axis
    cum_p = np.concatenate([[0.0], np.cumsum(
        0.5 * (pprime_used[1:] + pprime_used[:-1]) * np.diff(psi_n))])
    pressure = (cum_p[-1] - cum_p) * (-dpsi)          # = int_{psi_n}^{1} p' ds * (psi_axis-psi_bndry)
    cum_f = np.concatenate([[0.0], np.cumsum(
        0.5 * (ffprime_used[1:] + ffprime_used[:-1]) * np.diff(psi_n))])
    f2 = case["fvac"] ** 2 + 2.0 * (cum_f[-1] - cum_f) * (-dpsi)
    fpol = np.sign(case["fvac"]) * np.sqrt(np.maximum(f2, 0.0))
    say("p(0) = %.6g Pa,  f(0) = %.6g T m (fvac = %.6g)"
        % (pressure[0], fpol[0], case["fvac"]))
    rec["p_axis"] = float(pressure[0])
    rec["f_axis"] = float(fpol[0])
    rec["f2_went_negative"] = bool(np.any(f2 < 0.0))
    B0 = case["fvac"] / float(opt[0][0])
    rec["beta_tor_axis"] = float(2.0 * mu0 * pressure[0] / (case["fvac"] /
                                                            opt[0][0]) ** 2)
    rec["B_tor_at_axis_vacuum"] = float(B0)

    # ---- LCFS
    lcfs, missed = lcfs_from_rays(eq, opt, psi, psi_axis, psi_bndry, 512)
    if missed:
        lcfs2, missed2 = lcfs_from_rays(eq, opt, psi, psi_axis, psi_bndry, 512,
                                        offset=np.pi / 512.0)
        say("%d of 512 rays found no psi_n = 1 crossing; half-step offset "
            "leaves %d" % (missed, missed2))
        if missed2 < missed:
            lcfs, missed = lcfs2, missed2
    fpsi = sinterp.RectBivariateSpline(eq.R[:, 0], eq.Z[0, :], psi)
    lpsi = fpsi(lcfs[:, 0], lcfs[:, 1], grid=False)
    dpsi = float(np.max(np.abs(lpsi - psi_bndry)))
    span = abs(psi_axis - psi_bndry)
    gap = float(np.hypot(lcfs[0, 0] - lcfs[-1, 0], lcfs[0, 1] - lcfs[-1, 1]))
    step = float(np.max(np.hypot(np.diff(lcfs[:, 0], append=lcfs[0, 0]),
                                 np.diff(lcfs[:, 1], append=lcfs[0, 1]))))
    inside_grid = bool(lcfs[:, 0].min() > eq.Rmin
                       and lcfs[:, 0].max() < eq.Rmax
                       and lcfs[:, 1].min() > eq.Zmin
                       and lcfs[:, 1].max() < eq.Zmax)
    say("LCFS: %d/%d rays, max |psi-psi_bndry| = %.3e (%.2e of the flux span);"
        " closure gap %.2e m, max point spacing %.2e m; inside grid: %s"
        % (len(lcfs), 512, dpsi, dpsi / span, gap, step, inside_grid))
    rec.update(lcfs_npoints=int(len(lcfs)), lcfs_rays_missed=int(missed),
               lcfs_max_dpsi=dpsi, lcfs_max_dpsi_rel=float(dpsi / span),
               lcfs_closure_gap_m=gap, lcfs_max_point_spacing_m=step,
               lcfs_inside_grid=inside_grid)
    rec.update(shape_params(lcfs[:, 0], lcfs[:, 1]))

    rec.update(gs_residual(eq, psi, Jtor, mask))
    say("GS residual over the core, |Delta* psi + mu0 R Jtor| / max|mu0 R "
        "Jtor|:  2nd-order op Linf %.3e L2 %.3e ; 4th-order op Linf %.3e"
        % (rec["gs_resid_core_linf_2nd"], rec["gs_resid_core_l2_2nd"],
           rec["gs_resid_core_linf_4th"]))

    # THE CONDUCTORS, FLATTENED, at the currents this solve ended on. Taken
    # AFTER the solve for that reason: the control system moves them, and a
    # table read before it would describe the machine nobody converged.
    flat = conductors.flatten(tok)
    shrunk = conductors.shrink_to_fit(flat)
    for label, hw0, hh0, hw, hh in shrunk:
        say("conductor %s capped from %.4f x %.4f to %.4f x %.4f m half-extents "
            "so it does not overlap a neighbour or reach the axis"
            % (label, hw0, hh0, hw, hh))
    say("%d conductors, total current %+.6e A, all inside rho = %.4f"
        % (len(flat), sum(c["current"] for c in flat),
           conductors.bounding_radius(flat)))

    rec.update(
        psi_axis=psi_axis, psi_bndry=float(psi_bndry), Ip=Ip_ach,
        Ip_target=float(case["Ip"]), fvac=float(case["fvac"]),
        Raxis=float(opt[0][0]), Zaxis=float(opt[0][1]),
        order=int(case["order"]), nx=int(eq.nx), ny=int(eq.ny),
        ngrid=int(eq.nx * eq.ny),
        grid={k: float(v) for k, v in case["grid"].items()},
        dR=float(eq.dR), dZ=float(eq.dZ),
        n_xpoints=int(len(xpt)),
        xpoints=[[float(a) for a in row] for row in np.atleast_2d(xpt)]
        if len(xpt) else [],
        core_cells=int(mask.sum()),
        # THE CIRCUIT'S OWN CURRENT, which is what the control system solved
        # for and what a reader of this file wants to see. The per-CONDUCTOR
        # totals -- turns and multipliers resolved, solenoids summed -- are in
        # the .npz, because that is what a mesh and a [[coils]] block need and
        # they are a different question. conductors.py is the one authority on
        # the second.
        coil_currents={l: float(c.current) for l, c in tok.coils},
        conductors=[dict(c) for c in flat],
        # THE CASE'S OWN INPUTS, KEPT SEPARATE FROM ITS ANSWER.
        #
        # Everything else in this record is what the solve PRODUCED. These are
        # what it was ASKED FOR: the X-point locations the control system was
        # told to hit, the isoflux pairs it was told to equalise, the design
        # major radius and the target current. They are available before any
        # equilibrium exists -- they are a machine's operating point and a
        # scenario, which is what an experiment has.
        #
        # THE DISTINCTION IS THE WHOLE POINT. A comparison in which MEQ is
        # seeded from the reference's converged axis, its converged null and a
        # Green's sum over its converged Jtor is not a comparison of two
        # solvers; it is a measurement of how well MEQ polishes an answer it was
        # given. Anything MEQ is allowed to start from has to come from here.
        design=dict(
            xpoints=[[float(v) for v in p] for p in case["xpoints"]],
            isoflux=[[float(v) for v in p] for p in case["isoflux"]],
            R0=float(case["R0"]), Ip=float(case["Ip"]),
            fvac=float(case["fvac"]), frac_p=float(case["frac_p"]),
            pa=float(case["pa"]), pb=float(case["pb"]),
            fa=float(case["fa"]), fb=float(case["fb"]),
            seed_paxis=float(case["seed_paxis"]),
            limiter=(None if case["limiter"] is None
                     else [[float(v) for v in p] for p in case["limiter"]]),
        ),
        wall_time_s=float(time.time() - t0),
    )
    rec["sanity_psi_axis_ne_bndry"] = bool(
        abs(psi_axis - psi_bndry) > 1e-6 * max(1.0, abs(psi_axis)))
    rec["sanity_Ip_sign"] = bool(np.sign(Ip_ach) == np.sign(case["Ip"]))
    rec["sanity_Ip_relerr"] = float(abs(Ip_ach - case["Ip"]) / abs(case["Ip"]))
    rec["sanity_ok"] = bool(rec["converged"]
                            and rec["gs_resid_core_linf_4th"] < 1e-6
                            and rec["lcfs_inside_grid"]
                            and rec["sanity_psi_axis_ne_bndry"]
                            and rec["sanity_Ip_sign"]
                            and rec["lcfs_max_dpsi_rel"] < 1e-8
                            and rec["jtor_rebuild_relerr"] < 1e-4
                            and rec["boundary_kind"] == "diverted")

    rec["conventions"] = {
        "psi": "poloidal flux per radian, Wb/rad; psi[iR, iZ], R is axis 0",
        "grad_shafranov": "Delta* psi = -mu0 * R * Jtor, "
                          "Delta* = d2/dR2 + d2/dZ2 - (1/R) d/dR",
        "Jtor": "R*pprime(psi_n) + ffprime(psi_n)/(mu0*R) inside core_mask, "
                "zero outside",
        "psi_n": "(psi - psi_axis)/(psi_bndry - psi_axis)",
        "pprime": "dp/dpsi  (NOT dp/dpsi_n).  Do not divide by "
                  "(psi_bndry - psi_axis).",
        "ffprime": "f df/dpsi  (NOT f df/dpsi_n).  Same warning.",
        "pressure": "p(psi_n), Pa, from integrating pprime; p(1) = 0",
        "fpol": "f(psi_n) = R*B_tor, T m; f(1) = fvac",
        "fvac": "R*B_tor in vacuum.  It does NOT enter the GS solve in this "
                "formulation -- only f(psi_n) and hence B_tor.",
    }

    npz = os.path.join(OUTDIR, case["name"] + ".npz")
    np.savez_compressed(
        npz,
        R=eq.R[:, 0].astype(float),
        Z=eq.Z[0, :].astype(float),
        psi=psi.astype(float),
        psi_axis=np.float64(psi_axis),
        psi_bndry=np.float64(psi_bndry),
        lcfs_R=lcfs[:, 0].astype(float),
        lcfs_Z=lcfs[:, 1].astype(float),
        psi_n=psi_n.astype(float),
        pprime=pprime_used,
        ffprime=ffprime_used,
        pressure=pressure.astype(float),
        fpol=fpol.astype(float),
        fvac=np.float64(case["fvac"]),
        Ip=np.float64(Ip_ach),
        Raxis=np.float64(opt[0][0]),
        Zaxis=np.float64(opt[0][1]),
        machine=np.array(case["machine"]),
        order=np.int64(case["order"]),
        # extras
        Jtor=Jtor.astype(float),
        core_mask=mask.astype(np.int8),
        plasma_psi=eq.plasma_psi.astype(float),
        coil_psi=(psi - eq.plasma_psi).astype(float),
        Ip_logic_L=np.float64(L1),
        pprime_tabulated=(L1 * pprime_tab).astype(float),
        ffprime_tabulated=(L1 * ffprime_tab).astype(float),
        profile_Raxis=np.float64(1.0),
        wall_R=get_wall(case, eq)[0],
        wall_Z=get_wall(case, eq)[1],
        xpoints=(np.atleast_2d(xpt).astype(float) if len(xpt)
                 else np.zeros((0, 3))),
        # THE CONDUCTORS, SO THAT NOTHING DOWNSTREAM HAS TO KNOW THE MACHINE.
        # mkexactguess.py used to carry its own table of the limited case's four
        # coils and zip it against coil_currents POSITIONALLY, which is wrong
        # the moment a machine's labels come out in a different order -- and
        # TestTokamak's do: ['P1L','P1U','P2L','P2U'] against the limited
        # machine's ['P1U','P1L','P2U','P2L'].
        #
        # AND `float( c.R )` OFF EVERY ENTRY IS WRONG TOO, which is what this
        # replaces. Machine.coils holds four different classes and only two of
        # them have an R: MAST, TCV and MAST-U are built from Circuits and
        # Solenoids, so four of the seven diverted references SOLVED and could
        # not be saved. conductors.py flattens the tree into the rectangles MEQ
        # takes, resolving turns, circuit multipliers and solenoid windings
        # into one TOTAL current each -- see its header for why that number is
        # not `c.current`.
        **conductors.to_arrays(flat),
        # The design inputs again, as arrays, so a consumer reading only the
        # .npz can build a cold start without the .json. Empty arrays where a
        # case names none.
        design_xpoints=np.array(case["xpoints"], dtype=float).reshape(-1, 2),
        design_isoflux=np.array(case["isoflux"], dtype=float).reshape(-1, 4),
        design_R0=np.array(float(case["R0"])),
        design_Ip=np.array(float(case["Ip"])),
        boundary_kind=np.array(kind),
        psi_index_order=np.array("psi[iR, iZ]; R is axis 0, Z is axis 1"),
        pprime_units=np.array("dp/dpsi, Pa per (Wb/rad)"),
        ffprime_units=np.array("f df/dpsi, T^2 m^2 per (Wb/rad)"),
        interpolation_note=np.array(
            "Use a CUBIC SPLINE through (psi_n, pprime) and "
            "(psi_n, ffprime): that reproduces the solver's own Jtor to "
            "~4e-16.  Linear interpolation of the same points is ~3e-5."),
        source_form=np.array(
            "Jtor = R*pprime(psi_n) + ffprime(psi_n)/(mu0*R), "
            "psi_n = (psi-psi_axis)/(psi_bndry-psi_axis), zero outside "
            "core_mask; Delta* psi = -mu0*R*Jtor"),
    )
    rec["npz"] = npz
    rec["json"] = os.path.join(OUTDIR, case["name"] + ".json")
    with open(rec["json"], "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    say("wrote " + npz)
    rec["log"] = log
    return rec


def main():
    argv = list(sys.argv[1:])

    # RESOLUTION OVERRIDE, so the reference can be refined without editing the
    # case table.  --nx=N sets every case's grid to N x N.
    #
    # WHY THIS EXISTS. tools/README.md and docs/validation.rst both record that
    # MEQ SATURATES this benchmark -- its error stops falling at about 1.4e-04
    # because that is the REFERENCE's accuracy, set by the boundary fit and by
    # the contour extracted from a 129^2 grid.  "Refine the reference, not MEQ"
    # has been the standing next step since; this is the knob for it.
    #
    # freegs4e's grids are conventionally 2^n + 1, which is what its multigrid
    # solver wants; the direct 4th-order sparse solver does not care, but the
    # convention is kept so that a run at "2048" means 2049 and nests with 1025,
    # 513, 257 and 129 point for point on the coarse grid.  That nesting is the
    # whole reason to keep it: a reference refinement study wants the coarse
    # grid to be a subset of the fine one.
    # AND --seed-from IS WHAT SPENDS THE NESTING.  Without it every --nx run
    # starts COLD -- picard_loop() begins from whatever Equilibrium.__init__ left
    # in psi, at every resolution -- so a fine grid pays the full cold Picard
    # count, 23 to 103 steps, at its own O( n^3 ) per-step boundary cost.  The
    # boundary condition alone scales about 6.9x per grid doubling, so a cold
    # 1025^2 run is an afternoon where a seeded one is a handful of steps.
    #
    #     --seed-from=DIR    take the coarse answer from DIR/<case>.npz
    #     --seed-from=auto   take it from this run's own naming rule one
    #                        doubling down: 513 seeds from 257, 257 from 129
    #
    # A MISSING SEED IS A WARNING AND A WRONG ONE IS A REFUSAL, which is the
    # right way round: seeding is an optimisation, so its absence should cost
    # time and never correctness, while a coarse file from a different geometry
    # would interpolate, converge, and describe the wrong machine.
    # seed_from_coarse() checks the nesting and the extents and raises on both.
    nx = None
    vcycle = None
    seed_from = None
    rest = []
    for a in argv:
        if a.startswith("--nx="):
            nx = int(a.split("=", 1)[1])
        elif a.startswith("--vcycle="):
            vcycle = int(a.split("=", 1)[1])
        elif a.startswith("--seed-from="):
            seed_from = a.split("=", 1)[1]
        elif a == "--hagenow":
            global HAGENOW
            HAGENOW = True
        elif a.startswith("--rtol="):
            global RTOL
            RTOL = float(a.split("=", 1)[1])
        elif a == "--no-cache-boundary":
            global CACHE_BOUNDARY
            CACHE_BOUNDARY = False
        elif a.startswith("--hagenow-eps="):
            # "shipped" is freegs4e's own 1e-2 metres, kept so the two can be
            # measured against each other on one problem. Anything else is a
            # multiple of the cell. See WORKAROUND 5.
            global HAGENOW_EPS
            HAGENOW = True
            v = a.split("=", 1)[1]
            HAGENOW_EPS = None if v == "shipped" else float(v)
        else:
            rest.append(a)
    global VCYCLE_LEVELS
    VCYCLE_LEVELS = vcycle
    if nx is not None:
        n = nx if nx % 2 == 1 else nx + 1
        for case in CASES:
            case["grid"]["nx"] = n
            case["grid"]["ny"] = n
        global OUTDIR
        # A refined reference is large -- a 2049^2 case carries several
        # hundred megabytes of psi and B -- so it goes wherever FGSREF_OUT
        # says and not into the repository beside the 129^2 ones.
        base = os.environ.get("FGSREF_OUT")
        OUTDIR = (os.path.join(base, "n%d" % n) if base
                  else os.path.join(os.path.dirname(OUTDIR) or ".",
                                    os.path.basename(OUTDIR) + "-n%d" % n))
        os.makedirs(OUTDIR, exist_ok=True)
        print("resolution override: %d x %d, writing to %s" % (n, n, OUTDIR),
              flush=True)

    if seed_from is not None:
        global SEED_FROM
        if seed_from == "auto":
            # One doubling down, by the same naming rule OUTDIR uses above --
            # 513 -> 257 -> 129.  Deriving it rather than taking a path is what
            # makes a scan a loop over --nx with nothing else to keep in step.
            if nx is None:
                sys.exit("--seed-from=auto needs --nx: there is no coarser "
                         "grid to derive from the case table's own 129^2")
            n = nx if nx % 2 == 1 else nx + 1
            coarse = (n - 1)//2 + 1
            base = os.environ.get("FGSREF_OUT")
            here = os.path.dirname(os.path.abspath(__file__))
            SEED_FROM = (os.path.join(base, "n%d" % coarse) if base
                         else os.path.join(os.path.dirname(here) or ".",
                                           os.path.basename(here)
                                           + "-n%d" % coarse))
            if coarse == 129:
                # 129^2 is the case table's own resolution, so its output is
                # the undecorated directory rather than a -n129 one.
                # ANY case's own 129^2 output, not one case's by name. This
                # read "H_limited_circular.npz" while H was the only case that
                # had ever been laddered, which silently seeds the wrong case
                # -- or refuses to seed at all -- the moment a second one is.
                cand = [os.path.join(here, c["name"] + ".npz")
                        for c in CASES]
                if not os.path.isdir(SEED_FROM) and any(map(os.path.exists,
                                                            cand)):
                    SEED_FROM = here
        else:
            SEED_FROM = seed_from
        if not os.path.isdir(SEED_FROM):
            sys.exit("--seed-from: %s is not a directory" % SEED_FROM)
        print("seeding from %s" % SEED_FROM, flush=True)

    only = rest
    results = []
    for case in CASES:
        if only and not any(case["name"].startswith(o) for o in only):
            continue
        print("=" * 78, flush=True)
        print("CASE " + case["name"], flush=True)
        try:
            results.append(run_case(case))
        except Exception as exc:
            import traceback
            traceback.print_exc()
            results.append(dict(name=case["name"], machine=case["machine"],
                                failed=repr(exc)))
    print("=" * 78, flush=True)
    with open(os.path.join(OUTDIR, "summary.json"), "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True, default=str)
    for r in results:
        if "failed" in r:
            print("FAILED %-30s %s" % (r["name"], r["failed"]), flush=True)
        else:
            print("%-6s %-30s Ip=%9.4g psi_ax=%+.5g psi_b=%+.5g R0=%.3f "
                  "a=%.3f A=%.2f kappa=%.3f delta=%+.3f %s"
                  % ("OK" if r["sanity_ok"] else "CHECK", r["name"], r["Ip"],
                     r["psi_axis"], r["psi_bndry"], r["R0"], r["a"],
                     r["aspect_ratio"], r["kappa"], r["delta"],
                     r["boundary_kind"]), flush=True)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
