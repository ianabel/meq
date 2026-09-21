"""MEQ's fixed-boundary MXH cases, POSED FOR DESC -- and what the posing costs.

NAMED `convert_desc` AND NOT `convert`, WHICH COST AN AFTERNOON'S CONFUSION
ONCE.  `tools/freegs4e-benchmark/convert.py` already exists and every script
here puts BOTH directories on sys.path -- for `compare.py`, for `mxh.py` and for
`cores.py` -- so a module called `convert` in this directory is shadowed or
shadowing depending on the order of two `sys.path.insert` calls, and the symptom
is an AttributeError naming a constant rather than anything naming a path.

THE TWO CODES ARE NOT HANDED THE SAME OBJECT AND CANNOT BE.  MEQ takes the
Grad-Shafranov free functions against its own normalised POLOIDAL flux,

    F( r, z, psi ) = [ mu0 r^2 ( dp/dPsi )( Psi ) + ( g dg/dPsi )( Psi ) ] / psi_ax,
    Psi = psi / psi_ax,   Psi = 1 on the magnetic axis and 0 on Gamma

-- src/meq/Source.hpp, meq::NormalisedMHDSource -- which is the whole problem
statement and needs no equilibrium to write down.  DESC takes a boundary
surface, a total TOROIDAL flux, and two profiles against

    rho = sqrt( Phi / Phi_edge ),   Phi the enclosed toroidal flux,

namely the pressure p( rho ) and the enclosed toroidal current I( rho ).  The
*physics* is the same equilibrium; the *parameterisation* is not, and the map
between the two labels -- psi_n to rho -- is a functional of the solution.  So
the conversion cannot be done without an equilibrium, and this file is honest
about which one it uses.

THREE THINGS THE CONVERSION NEEDS THAT MEQ'S INPUT DOES NOT CARRY, each of
which is a real asymmetry and not a detail:

  1. THE FLUX-LABEL MAP.  Phi( Psi ) requires the shapes of the interior flux
     surfaces.  MEQ computes them as part of solving; DESC has to be told.
  2. THE ABSOLUTE TOROIDAL FIELD.  Grad-Shafranov depends on g only through
     g dg/dpsi, so MEQ's input fixes g^2 up to an additive constant and MEQ's
     answer for psi is independent of it.  DESC's `Psi` is the toroidal flux,
     which is proportional to g, so DESC must be given the constant.  It comes
     from the reference's `fpol`.  A different g0 gives DESC a different
     `Psi` and a different rho label -- and, at convergence, the SAME psi( R, Z ).
     That invariance is asserted in `--check`.
  3. THE SHAPE OF THE PROFILES IN rho.  DESC wants p and I as even power
     series in rho (or splines).  MEQ's tables are cubic splines in Psi.  The
     re-fit is a truncation and its residual is measured and printed.

WHICH EQUILIBRIUM SUPPLIES THE MAP, and the two answers this file implements:

  * `--from reference` (the default).  The freegs4e 257^2 solve the case was
    generated from -- `tools/freegs4e-benchmark/ref-n257` -- supplies psi( R, Z )
    and therefore Phi( Psi ) and I( Psi ).  This is the SAME object M-147
    measures MEQ against, so nothing new enters the comparison.  It is not free:
    DESC's input then inherits the reference's own geometry error, which the
    outermost surface bounds at 6.0e-05 relative in the enclosed area for
    `fixed-h-circular` -- the psi_n = 0.95 contour against the MXH curve that
    was fitted to it, printed by `--check`.
  * `--from <equilibrium.npz>`.  Any psi( R, Z ) on a grid, in MEQ's gauge or
    the reference's.  Feeding DESC's own converged answer back in is the
    fixed-point iteration that removes the reference from the posing
    altogether; `--self-consistent` in descrun.py drives it.

WHY THE MAP IS BUILT BY GREEN'S THEOREM AND NOT BY SUMMING CELLS.  Both
integrands are of the form a( Psi ) r + b( Psi )/r, so

    int_{A(u)} h( psi_n ) dA  =  h( u ) V( u ) - int_0^u h'( v ) V( v ) dv

with V( u ) the enclosed moment of r or 1/r, and each moment follows from the
bounding CONTOUR alone:

    int r dA = oint ( r^2 / 2 ) dz,        int ( 1/r ) dA = oint ln r dz.

A hard mask over grid cells is O( h ) in the boundary and reads the enclosed
current 3.9e-04 wrong on the shipped reference; the contour route is O( h^2 )
in the contour position and needs no differentiation of V.  Measured both ways
in `--check`.
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from skimage.measure import find_contours

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "freegs4e-benchmark"))
from mxh import mxh_boundary                                  # noqa: E402

MU0 = 4.0e-7 * np.pi
HERE = os.path.dirname(os.path.abspath(__file__))
MEQ = os.path.abspath(os.path.join(HERE, "..", ".."))
REFDIR = os.path.join(MEQ, "tools", "freegs4e-benchmark", "ref-n257")

#: stem -> the freegs4e reference it was generated from, out of make_fixed.sh.
CASES = {
	"fixed-h-circular":    "H_limited_circular_shaped",
	"fixed-a-testtokamak": "A_testtokamak_classic_shaped",
	"fixed-e-diamagnetic": "E_testtokamak_diamagnetic_shaped",
	"fixed-d-tcv":         "D_tcv_conventional_shaped",
	"fixed-g-mastu":       "G_mastu_simple_shaped",
	"fixed-f-diiid":       "F_diiid_conventional_shaped",
}


def read_profile_table(path):
	"""( Psi, value ) out of one of MEQ's three-column .dat tables.

	The third column is d( value )/dPsi, which meq::SplineProfile uses and this
	does not: a cubic spline through the pairs reproduces it, and reading a
	derivative that was itself produced by a spline would double-count the
	smoothing rather than improve on it.
	"""
	data = np.loadtxt(path)
	return data[:, 0], data[:, 1]


def load_case(stem, examples=None, refdir=REFDIR):
	"""Everything one shipped case carries: meta, both tables, the reference."""
	examples = examples or os.path.join(MEQ, "examples")
	meta = json.load(open(os.path.join(examples, f"{stem}-meta.json")))
	psi_p, dp = read_profile_table(os.path.join(examples, f"{stem}-pprime.dat"))
	psi_g, gg = read_profile_table(os.path.join(examples, f"{stem}-ggprime.dat"))
	ref = np.load(os.path.join(refdir, f"{CASES[stem]}.npz"), allow_pickle=True)
	return dict(stem=stem, meta=meta, pprime=(psi_p, dp), ggprime=(psi_g, gg),
	            ref=ref)


def mxh_curve(meta, n=4096):
	"""Gamma itself, counterclockwise from the outboard midplane.

	This is the curve MEQ solves inside -- `[boundary.shape] Type = "mxh"` --
	so it is exact data rather than a fit to anything, whatever the residual
	against the contour it was fitted FROM.
	"""
	return mxh_boundary(meta["R0"], meta["Z0"], meta["minor_radius"],
	                    meta["kappa"], np.array(meta["cos"]),
	                    np.array(meta["sin"]), n)


def moments(rr, zz):
	"""( area, int r dA, int 1/r dA ) for the region a closed curve encloses.

	Green's theorem on the closed polygon, orientation-corrected, so the caller
	need not know which way round the contour came out.
	"""
	dz = np.roll(zz, -1) - zz
	rmid = 0.5 * (np.roll(rr, -1) + rr)
	area = float(np.sum(rmid * dz))
	sign = 1.0 if area > 0.0 else -1.0
	v_r = sign * float(np.sum(0.25 * (np.roll(rr, -1) ** 2 + rr ** 2) * dz))
	v_i = sign * float(np.sum(np.log(rmid) * dz))
	return abs(area), v_r, v_i


class FluxGeometry:
	"""The enclosed moments V_r( u ), V_1/r( u ) of one equilibrium's surfaces.

	`u` is the reference's own normalised poloidal flux psi_n -- 0 on the axis,
	1 on the separatrix -- because that is the label the grid carries.  MEQ's
	Psi is 1 - u/level and the caller converts.

	THE OUTERMOST LEVEL IS TAKEN FROM THE MXH CURVE AND NOT FROM THE CONTOUR,
	when one is supplied.  Gamma is analytic and exactly known; the contour at
	that level is a 257^2 extraction of the curve it was fitted from, and using
	it would put the fit residual into the total toroidal flux -- the one number
	that is not a shape and cannot be corrected later.  The gap between the two
	is reported rather than smoothed over.
	"""

	def __init__(self, R, Z, psin, axis, level, n_levels=120,
	             outer_curve=None, sample_max=0.99):
		self.level = level
		self.axis = axis
		self._R, self._Z, self._psin = R, Z, psin
		# STOP SHORT OF Gamma AND TAKE THE LAST SURFACE FROM THE CURVE.
		# A field supplied by MEQ has NO DATA outside Gamma -- `inside` is
		# exactly that region -- so a contour asked for AT Gamma runs along the
		# edge of the data and a contour asked for just inside it is clean.
		# Gamma's own moments are then the MXH curve's, which is exact.  The
		# reference's field does extend past Gamma and does not need this; it
		# is done for both so that the two paths differ in the FIELD alone.
		levels = np.linspace(level / n_levels, sample_max * level, n_levels)
		us, v_r, v_i, ar = [0.0], [0.0], [0.0], [0.0]
		for value in levels:
			curve = self._contour(value)
			if curve is None:
				continue
			a, vr, vi = moments(*curve)
			us.append(float(value)); v_r.append(vr); v_i.append(vi); ar.append(a)
		at_gamma = self._contour(level)
		self.contour_outer = moments(*at_gamma) if at_gamma else (np.nan,) * 3
		if outer_curve is not None:
			a, vr, vi = moments(*outer_curve)
			self.curve_outer = (a, vr, vi)
			us.append(level); v_r.append(vr); v_i.append(vi); ar.append(a)
		else:
			self.curve_outer = self.contour_outer
			us.append(level); v_r.append(self.contour_outer[1])
			v_i.append(self.contour_outer[2]); ar.append(self.contour_outer[0])
		self.u = np.array(us)
		self.V_r = CubicSpline(self.u, np.array(v_r))
		self.V_i = CubicSpline(self.u, np.array(v_i))
		self.area = CubicSpline(self.u, np.array(ar))
		self.n_levels_used = len(us) - 1

	def _contour(self, value):
		"""The closed psi_n = value contour that encircles the magnetic axis.

		A level set can have several components -- the private flux beyond an
		X-point is the one that matters here, and it carries the same psi_n by
		VALUE while being no part of the plasma.  Selecting on the winding
		number about the axis is what separates them; a largest-component rule
		does not, because on a near-double-null case the private region can be
		the larger.
		"""
		Rgrid, Zgrid, psin = self._R, self._Z, self._psin
		best, best_len = None, 0
		for c in find_contours(psin, value):
			rr = Rgrid[0] + c[:, 0] * (Rgrid[1] - Rgrid[0])
			zz = Zgrid[0] + c[:, 1] * (Zgrid[1] - Zgrid[0])
			if np.hypot(rr[0] - rr[-1], zz[0] - zz[-1]) > 1e-9:
				continue                        # open: it leaves the box
			ang = np.unwrap(np.arctan2(zz - self.axis[1], rr - self.axis[0]))
			if abs(ang[-1] - ang[0]) < 6.0:     # does not enclose the axis
				continue
			if len(rr) > best_len:
				best, best_len = (rr, zz), len(rr)
		return best

	def enclosed(self, profile, which, u):
		"""int_{A(u)} profile( psi_n ) * ( r or 1/r ) dA, by parts.

		`profile` is a CubicSpline in psi_n, so its derivative is exact and no
		numerical differentiation of V is needed anywhere.
		"""
		V = self.V_r if which == "r" else self.V_i
		grid = np.linspace(0.0, u, 2001)
		return float(profile(u) * V(u)
		             - np.trapezoid(profile.derivative()(grid) * V(grid), grid))


def refine_field(field, factor):
	"""Bicubic upsampling of psi_n before any contour is taken from it.

	THE CONTOUR CHORDS THE SURFACE AND THE BIAS IS ONE-SIDED, WHICH IS WHY THIS
	IS NOT A REFINEMENT BUT A CORRECTION.  `find_contours` walks a marching-
	squares polygon with linearly interpolated crossings, so a convex surface
	comes out as an inscribed polygon and every enclosed moment reads LOW.  The
	error is O( h^2 ) per segment but the number of segments falls with the
	surface, so it is worst on the innermost surfaces -- exactly the ones the
	rho label near the axis is built from.  Measured on `fixed-h-circular`,
	`int r dA` inside psi_n = 0.1 and the enclosed current at Gamma:

	    factor      1          2          4
	    V_r( 0.1 )  1.760669e-02  1.763618e-02  1.764379e-02
	    I( rho=1 )  2.9974710e+05 2.9983552e+05 2.9985773e+05

	which is second order and converging from below, and the last column agrees
	with the SAME quantity taken off a MEQ .nc whose grid is 3.4x finer over
	the plasma ( 2.9985911e+05 ) to 4.6e-06.  Two different fields agreeing is
	what says this is the instrument rather than the equilibrium.

	**IT COSTS FOUR HUNDRED-ODD SPLINE EVALUATIONS AND IT IS WORTH 4.2e-04 IN
	DESC'S ANSWER**, which is an order of magnitude more than any other item in
	the conversion.  psi_ax follows the enclosed current almost one for one.

	NOT APPLIED TO A FIELD WITH NO DATA OUTSIDE Gamma.  A MEQ .nc is filled
	with a sentinel beyond Gamma so that find_contours has somewhere to stop,
	and a bicubic spline rings across that jump for several cells INWARD.  The
	default is therefore 4 for the reference and 1 for a .nc, and a .nc does
	not need it: its grid is already finer over the plasma than the
	reference's.
	"""
	if factor <= 1:
		return field
	from scipy.interpolate import RectBivariateSpline
	R, Z, psin = field["R"], field["Z"], field["psin"]
	spline = RectBivariateSpline(R, Z, psin, kx=3, ky=3)
	Rn = np.linspace(R[0], R[-1], (len(R) - 1) * factor + 1)
	Zn = np.linspace(Z[0], Z[-1], (len(Z) - 1) * factor + 1)
	out = dict(field)
	out.update(R=Rn, Z=Zn, psin=spline(Rn, Zn), refine=factor)
	return out


def field_from_reference(ref, level, refine=4):
	"""( R, Z, psi_n[ iR, iZ ], axis, psi_ax ) out of a freegs4e reference."""
	psi_axis, psi_bndry = float(ref["psi_axis"]), float(ref["psi_bndry"])
	psi = np.array(ref["psi"], float)
	return refine_field(
		dict(R=np.array(ref["R"], float), Z=np.array(ref["Z"], float),
		     psin=(psi_axis - psi) / (psi_axis - psi_bndry),
		     axis=(float(ref["Raxis"]), float(ref["Zaxis"])),
		     psi_ax=level * (psi_axis - psi_bndry), source="reference",
		     refine=1), refine)


def field_from_grid(R, Z, psi_meq, level, axis=None, source="grid"):
	"""A flux-label field from psi ON A GRID IN MEQ's GAUGE -- 0 on Gamma, the
	axis flux at the axis, NaN or non-finite outside Gamma.

	The in-memory half of field_from_meq() and field_from_desc(), which differ
	only in where the grid came from.  Having one implementation is the point:
	a self-consistent iteration compares a map built from DESC's field against
	one built from MEQ's, and two transcriptions of `psin` would put a
	difference between them that is neither code's.

	OUTSIDE Gamma THE LABEL IS PUT BEYOND IT rather than left at NaN, for
	find_contours' sake: it walks an array and cannot step over a hole, and
	every contour this is used for is interior.
	"""
	psi_meq = np.asarray(psi_meq, float)
	good = np.isfinite(psi_meq)
	if not good.any():
		raise ValueError("%s: the field is entirely non-finite" % source)
	psi_ax = float(np.nanmax(psi_meq))
	if not psi_ax > 0.0:
		raise ValueError("%s: psi_ax came out %.3e, so this is not a MEQ-gauge "
		                 "field with psi = 0 on Gamma" % (source, psi_ax))
	psin = np.where(good, level * (1.0 - psi_meq / psi_ax), 2.0 * level)
	if axis is None:
		flat = int(np.nanargmax(np.where(good, psi_meq, -np.inf)))
		axis = (float(R[flat // len(Z)]), float(Z[flat % len(Z)]))
	return dict(R=np.asarray(R, float), Z=np.asarray(Z, float), psin=psin,
	            axis=axis, psi_ax=psi_ax, source=source)


def field_from_desc(path, level):
	"""The same, out of descrun.py's own .npz.

	**THIS IS WHAT MAKES THE POSING SELF-CONSISTENT AND IT IS THE ONLY ROUTE
	THAT REMOVES BOTH CODES' ANSWERS FROM IT.**  `--from reference` inherits
	freegs4e's interior surfaces and `--from <meq.nc>` inherits MEQ's; feeding
	DESC its own converged field and re-solving is a fixed-point iteration
	whose limit depends on neither.  What survives is the MXH curve, the two
	profile tables, and the two scalars of section "what MEQ's input does not
	fix" -- which is exactly the problem MEQ is handed.

	descrun.py writes `psi_meq` beside the reference-gauge `psi` for this, so
	nothing here has to undo a gauge shift it did not apply.
	"""
	d = np.load(path, allow_pickle=True)
	if "psi_meq" not in d.files:
		raise SystemExit(
			"%s predates the self-consistent loop: it carries no `psi_meq`. "
			"Re-run descrun.py to write one." % path)
	axis = ((float(d["axis_r"]), float(d["axis_z"]))
	        if "axis_r" in d.files else None)
	return field_from_grid(d["R"], d["Z"], d["psi_meq"], level, axis=axis,
	                       source=os.path.basename(path))


def field_from_meq(path, level):
	"""The same, out of a MEQ .nc -- so the map comes from the EXACT problem.

	WHY THIS IS A DIFFERENT AND BETTER POSING, not merely another option.  The
	reference solves the FREE-boundary problem by finite differences on a
	257^2 grid, and M-147 measures MEQ against it at 5.0e-05 to 1.1e-04; its
	interior surfaces are therefore its own, not those of the problem MEQ and
	DESC are being raced on.  A .nc is the answer to the problem ACTUALLY
	POSED -- the MXH curve, those two tables -- so the flux-label map built
	from it is the map that problem has.  What it costs is that DESC's input
	now depends on MEQ's output, which is why `descrun.py --self-consistent`
	exists: iterating the map on DESC's OWN field removes both codes'
	answers from the posing and leaves only the curve and the tables.

	THE BAND IS NOT DROPPED HERE AND THAT IS DELIBERATE.  It is a Taylor
	continuation of the solved flux using the SOLVED q, so it is a smooth
	extension rather than junk, and dropping it would leave the outermost
	contours running along the edge of the data instead of through it.  The
	outermost surface is taken from the MXH curve regardless, so nothing the
	band decides reaches Gamma's own moments.
	"""
	from netCDF4 import Dataset
	with Dataset(path) as ds:
		R = np.array(ds["R"][:], float)
		Z = np.array(ds["Z"][:], float)
		psi = np.array(ds["psi"][:], float).T              # ( Z, R ) -> [ iR, iZ ]
		inside = np.array(ds["inside"][:]).astype(bool).T
		psi_ax = float(ds.getncattr("psi_axis"))
		axis = (float(ds.getncattr("axis_r")), float(ds.getncattr("axis_z")))
	# OUTSIDE Gamma THERE IS NO SOLUTION, so put those nodes BEYOND Gamma in
	# the label rather than at NaN, which find_contours cannot walk past.
	psin = np.where(inside & np.isfinite(psi),
	                level * (1.0 - psi / psi_ax), 2.0 * level)
	return dict(R=R, Z=Z, psin=psin, axis=axis, psi_ax=psi_ax,
	            source=os.path.basename(path))


def build(case, level=None, n_levels=120, n_fit=10, sym_tol=1e-6,
          field=None, g_boundary=None, p_boundary=None):
	"""The DESC posing of one MEQ case, plus everything needed to audit it."""
	meta, ref = case["meta"], case["ref"]
	level = float(meta["level"]) if level is None else level
	field = field or field_from_reference(ref, level)
	psi_axis, psi_bndry = float(ref["psi_axis"]), float(ref["psi_bndry"])
	R, Z, psin, axis = field["R"], field["Z"], field["psin"], field["axis"]

	geom = FluxGeometry(R, Z, psin, axis, level, n_levels,
	                    outer_curve=mxh_curve(meta))

	# ---- MEQ's own two tables, re-expressed against psi_n --------------------
	# Psi = 1 - psi_n/level, and MEQ's tables are d/dPsi, so one chain-rule
	# factor of psi_ax converts to d/dpsi.  psi_ax is the MEQ-gauge axis flux
	# OF THE EQUILIBRIUM SUPPLYING THE MAP -- the reference's is
	# level*( psi_axis - psi_bndry ), a .nc carries its own.
	psi_ax = float(field["psi_ax"])
	Psi_p, dpdPsi = case["pprime"]
	Psi_g, ggdPsi = case["ggprime"]
	u_p, u_g = level * (1.0 - Psi_p), level * (1.0 - Psi_g)
	order_p, order_g = np.argsort(u_p), np.argsort(u_g)
	# d p / d psi and ( g dg/d psi ), both per Wb/rad, as functions of psi_n.
	dpdpsi = CubicSpline(u_p[order_p], dpdPsi[order_p] / psi_ax)
	ggdpsi = CubicSpline(u_g[order_g], ggdPsi[order_g] / psi_ax)

	# ---- the toroidal field, which MEQ's input does not fix ------------------
	# g^2( u ) = g^2( level ) + 2 int_level^u ( g dg/dpsi ) dpsi, and
	# dpsi = -( psi_axis - psi_bndry ) d( psi_n ).  The constant comes from the
	# reference's own fpol at the boundary surface: MEQ's answer does not
	# depend on it, DESC's `Psi` does, and `--check` asserts the invariance.
	# THE CASE FILE FIRST, THE REFERENCE ONLY AS A FALLBACK.  `g_at_gamma` is
	# in every <stem>-meta.json since make_case.py started writing it, for
	# exactly this: taking it from the reference put a THIRD code into the
	# posing of a two-code comparison.  An explicit --g-boundary still wins
	# over both.
	if g_boundary is not None:
		g_edge, g_source = float(g_boundary), "--g-boundary"
	elif "g_at_gamma" in meta:
		g_edge, g_source = float(meta["g_at_gamma"]), "meta"
	else:
		g_edge = float(CubicSpline(np.array(ref["psi_n"], float),
		                           np.array(ref["fpol"], float))(level))
		g_source = "reference fpol"
	grid = np.linspace(0.0, level, 4001)
	dpsi_du = -psi_ax / level
	# int_level^u ( . ) dpsi = -int_u^level ( . ) ( dpsi/du ) du
	g2 = np.array([g_edge ** 2 - 2.0 * np.trapezoid(
		ggdpsi(np.linspace(u, level, 512)) * dpsi_du,
		np.linspace(u, level, 512)) for u in grid])
	g_of_u = CubicSpline(grid, np.sqrt(np.maximum(g2, 1e-30)))

	# ---- the two DESC profiles, against psi_n first --------------------------
	sample = np.linspace(0.0, level, 241)
	phi = np.array([geom.enclosed(g_of_u, "inv", u) for u in sample])
	cur = np.array([geom.enclosed(dpdpsi, "r", u)
	                + geom.enclosed(ggdpsi, "inv", u) / MU0 for u in sample])
	if p_boundary is None and "p_at_gamma" in meta:
		p_boundary, p_source = float(meta["p_at_gamma"]), "meta"
	elif p_boundary is not None:
		p_source = "--p-boundary"
	else:
		p_source = "reference pressure"
	if p_boundary is None:
		p_of_u = CubicSpline(np.array(ref["psi_n"], float),
		                     np.array(ref["pressure"], float))
		pres = p_of_u(sample)
	else:
		# p FROM MEQ'S OWN TABLE, WHICH DETERMINES IT UP TO THE VALUE ON Gamma.
		#
		#     p( u ) = p( level ) + int_level^u ( dp/dpsi ) dpsi,
		#     dpsi = ( dpsi/du ) du
		#
		# BY THE SPLINE'S OWN ANTIDERIVATIVE AND NOT BY QUADRATURE.  `dpdpsi`
		# IS a cubic spline, so its integral is exact arithmetic and a
		# trapezoid over it is strictly worse for more work.
		#
		# **AND IT DOES NOT REPRODUCE THE REFERENCE'S OWN `pressure` ARRAY,
		# WHICH IS THE REFERENCE'S FAULT AND NOT THIS ROUTE'S.**  Measured on
		# fixed-h-circular, the axis pressure three ways:
		#
		#   MEQ's table integrated              44430.091294 Pa
		#   the reference's own pprime integrated   44430.091294 Pa
		#   the reference's own pressure array      44430.415850 Pa
		#
		# The first two agree to **4.9e-15** -- so MEQ's re-tabulation of
		# `pprime` is exact and loses nothing -- and both sit **7.3e-06** from
		# the stored array, which is freegs4e's `pressure` not being the
		# integral of its own `pprime`.  So this is a BETTER-defined number
		# than the array, being consistent with the `p'` MEQ is actually
		# given, and the switch to it is not a loss of accuracy.  It matters
		# because this is now the DEFAULT: `p_at_gamma` is in every meta file,
		# so the reference's `pressure` array is not read at all.
		anti = dpdpsi.antiderivative()
		pres = (float(p_boundary)
		        + dpsi_du * (anti(sample) - anti(level)))

	phi_edge = float(phi[-1])
	rho = np.sqrt(np.clip(phi / phi_edge, 0.0, None))

	# ---- the same two profiles as knot data in rho ---------------------------
	# WHY A SPLINE IS OFFERED BESIDE THE POWER SERIES, AND IT IS NOT A
	# PREFERENCE.  DESC's force balance uses dp/drho, and a least-squares fit
	# to p VALUES does not control it.  Measured on fixed-h-circular, worst
	# relative error in dp/drho over rho > 0.02 against the true profile:
	#
	#     terms      8       12      16      20      28
	#     value   7.7e-04  5.2e-04 3.8e-04 3.0e-04 2.9e-04
	#     d/drho  1.5e-02  6.7e-03 7.4e-03 1.1e-02 1.7e-02
	#
	# The VALUE improves with order and the DERIVATIVE stops improving at
	# twelve terms and then gets WORSE, which is the monomial basis losing
	# conditioning rather than anything about the profile.  So "fit it harder"
	# is not available and a different representation is.
	#
	# THE AXIS SLOPE IS ENFORCED RATHER THAN FITTED.  Every flux function here
	# is smooth in rho^2, so df/drho vanishes at rho = 0 exactly; a spline
	# through values alone puts a small nonzero slope there and an O( 1 )
	# relative error in dp/drho on the innermost surfaces, where r^2 p' is the
	# whole source.
	knots = np.linspace(0.0, 1.0, 65)
	rho_sorted = np.argsort(rho)
	p_of_rho = CubicSpline(rho[rho_sorted], pres[rho_sorted])
	c_of_rho = CubicSpline(rho[rho_sorted], cur[rho_sorted])
	p_knots, c_knots = p_of_rho(knots), c_of_rho(knots)
	dp_knots, dc_knots = p_of_rho.derivative()(knots), c_of_rho.derivative()(knots)
	dp_knots[0] = 0.0
	dc_knots[0] = 0.0
	c_knots[0] = 0.0                       # no current is enclosed by nothing

	# ---- refit onto even powers of rho ---------------------------------------
	# EVEN POWERS BECAUSE THE PHYSICS IS EVEN.  Phi ~ psi_n near the axis, so
	# rho ~ sqrt( psi_n ) and every flux function is a smooth function of rho^2.
	# An odd term would put a |rho| cusp on the axis that no equilibrium has.
	p_coef, p_res = fit_even(rho, pres, n_fit, start=0)
	c_coef, c_res = fit_even(rho, cur, n_fit, start=1)

	return dict(
		stem=case["stem"], meta=meta, level=level, geom=geom,
		psi_ax_map=psi_ax, map_source=field["source"],
		g_source=g_source, p_source=p_source, g_boundary=g_edge,
		map_refine=int(field.get("refine", 1)),
		psi_ax_reference=level * (psi_axis - psi_bndry),
		psi_axis=psi_axis, psi_bndry=psi_bndry,
		Psi_total=phi_edge, g_edge=g_edge, g_axis=float(g_of_u(0.0)),
		u=sample, rho=rho, pressure=pres, current=cur, toroidal_flux=phi,
		p_coefficients=p_coef, current_coefficients=c_coef,
		p_residual=p_res, current_residual=c_res,
		knots=knots, p_knots=p_knots, current_knots=c_knots,
		dp_knots=dp_knots, dcurrent_knots=dc_knots,
		dpdpsi=dpdpsi, ggdpsi=ggdpsi, g_of_u=g_of_u,
		up_down_symmetric=bool(abs(meta["asymmetry"]) < sym_tol),
	)


def fit_even(rho, values, n_terms, start=0):
	"""Least-squares fit to rho^(2*start), ..., rho^(2*(start+n_terms-1)).

	Returns ( coefficients indexed by the POWER, relative residual ).  The
	residual is reported because it is one of the three inexactnesses this
	conversion has and the only one the caller can trade against cost.
	"""
	powers = 2 * np.arange(start, start + n_terms)
	A = rho[:, None] ** powers[None, :]
	coef, *_ = np.linalg.lstsq(A, values, rcond=None)
	fit = A @ coef
	scale = np.abs(values).max() or 1.0
	out = np.zeros(int(powers[-1]) + 1)
	out[powers] = coef
	return out, float(np.abs(fit - values).max() / scale)


def desc_boundary_samples(meta, n=512):
	"""Gamma sampled at DESC's poloidal angle, which runs the OTHER WAY.

	DESC's own default torus is R = 10 + cos( theta ), Z = -sin( theta ) --
	`FourierRZToroidalSurface.__init__`'s defaults, and `desc/examples/DSHAPE`
	has Z1 = -1.47 at m = -1 beside R1 = +1.0 at m = +1 -- so theta descends
	from the outboard midplane and is CLOCKWISE in the ( R, z ) plane.  MXH's
	t is counterclockwise.  Sampling at t = -theta is the whole of the
	conversion, and getting it backwards gives a left-handed surface that DESC
	silently flips for you, after which its theta is not the one you think.
	"""
	theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
	t = -theta
	tR = t + meta["cos"][0]
	for n_ in range(1, len(meta["cos"])):
		tR = tR + meta["cos"][n_] * np.cos(n_ * t)
	for n_ in range(1, len(meta["sin"]) + 1):
		tR = tR + meta["sin"][n_ - 1] * np.sin(n_ * t)
	a, kappa = meta["minor_radius"], meta["kappa"]
	return theta, meta["R0"] + a * np.cos(tR), meta["Z0"] + kappa * a * np.sin(t)


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", nargs="?", default="fixed-h-circular",
	                choices=sorted(CASES) + [None])
	ap.add_argument("--levels", type=int, default=120)
	ap.add_argument("--terms", type=int, default=10)
	ap.add_argument("--refine", type=int, default=4,
	                help="bicubic upsampling of the field before contouring; "
	                     "see refine_field for why it is a correction rather "
	                     "than a refinement, and why a .nc gets 1")
	ap.add_argument("--from", dest="source", default="",
	                help="a MEQ .nc whose psi supplies the flux-label map, "
	                     "instead of the freegs4e reference")
	ap.add_argument("--check", action="store_true",
	                help="the cross-checks: contour against the MXH curve, "
	                     "Green's against a cell sum, and the fit residuals")
	args = ap.parse_args()

	case = load_case(args.stem)
	level = float(case["meta"]["level"])
	field = (field_from_meq(args.source, level) if args.source
	         else field_from_reference(case["ref"], level, args.refine))
	out = build(case, n_levels=args.levels, n_fit=args.terms, field=field)
	print(f"  {args.stem}   flux-label map from {out['map_source']}"
	      f" x{out['map_refine']}")
	print(f"    Gamma at psi_n            {out['level']}")
	print(f"    psi_ax the map was built with  {out['psi_ax_map']:.9e} Wb/rad")
	print(f"    psi_ax the reference reports  {out['psi_ax_reference']:.9e} Wb/rad")
	print(f"    total toroidal flux Psi   {out['Psi_total']:.9e} Wb")
	# WHERE THE TWO GRAD-SHAFRANOV-BLIND SCALARS CAME FROM, printed because
	# "reference" here means a third code is in the posing of a two-code
	# comparison and nothing else on this page would say so.
	print(f"    p on Gamma from {out['p_source']:18s}  "
	      f"g on Gamma from {out['g_source']}")
	print(f"    g on Gamma / on the axis  {out['g_edge']:.6f} / {out['g_axis']:.6f} T m")
	print(f"    enclosed current at rho=1 {out['current'][-1]:.6e} A")
	print(f"    p on Gamma / on the axis  {out['pressure'][-1]:.4f} / {out['pressure'][0]:.4f} Pa")
	print(f"    up-down symmetric         {out['up_down_symmetric']}")
	print(f"    even-power fit residual   p {out['p_residual']:.3e}   I {out['current_residual']:.3e}")
	print(f"    surfaces resolved         {out['geom'].n_levels_used} of {args.levels}")

	if args.check:
		a_c, vr_c, vi_c = out["geom"].contour_outer
		a_m, vr_m, vi_m = out["geom"].curve_outer
		print("\n    THE OUTERMOST SURFACE: the reference's contour against Gamma itself")
		print(f"      area   contour {a_c:.8e}  MXH {a_m:.8e}  rel {abs(a_c-a_m)/a_m:.3e}")
		print(f"      int r  contour {vr_c:.8e}  MXH {vr_m:.8e}  rel {abs(vr_c-vr_m)/vr_m:.3e}")
		print(f"      int 1/r contour {vi_c:.8e}  MXH {vi_m:.8e}  rel {abs(vi_c-vi_m)/vi_m:.3e}")

		ref = case["ref"]
		R, Z = np.array(ref["R"]), np.array(ref["Z"])
		RR, _ = np.meshgrid(R, Z, indexing="ij")
		psin = ((float(ref["psi_axis"]) - np.array(ref["psi"]))
		        / (float(ref["psi_axis"]) - float(ref["psi_bndry"])))
		mask = np.array(ref["core_mask"]).astype(bool) & (psin <= out["level"])
		clipped = np.clip(psin, 0.0, 1.0)
		jphi = (RR * out["dpdpsi"](clipped)
		        + out["ggdpsi"](clipped) / (MU0 * RR))
		cells = float((jphi * mask).sum() * (R[1] - R[0]) * (Z[1] - Z[0]))
		green = out["current"][-1]
		print("\n    THE ENCLOSED CURRENT: Green's theorem against a cell sum")
		print(f"      Green   {green:.8e} A")
		print(f"      cells   {cells:.8e} A     rel {abs(green-cells)/abs(green):.3e}")
		print("      the cell sum is the O( h ) one; this is its error, not Green's")
	return 0


if __name__ == "__main__":
	sys.exit(main())
