#!/usr/bin/env python
"""MEQ's six shipped fixed-boundary cases, posed for CHEASE.

CHEASE is the EPFL/SPC fixed-boundary Grad-Shafranov code -- Fortran, bicubic
Hermite finite elements on a flux-aligned mesh, Luetjens/Bondeson/Sauter CPC 97
(1996) 219.  It solves the SAME equation MEQ does, from the SAME two free
functions, which makes it a closer comparison than either `freegs4e` (different
algorithm, same equation) or DESC (different unknown entirely).

**THE CONVERSION NEEDS NOTHING FROM MEQ'S ANSWER AND NOTHING FROM THE
REFERENCE.**  That is the one structural advantage this arm has over
`tools/desc-benchmark`, whose flux-label map has to come from somebody's
converged equilibrium.  CHEASE takes `dp/dpsi` and `T dT/dpsi` against the
normalised poloidal flux, which is exactly what `examples/fixed-*-pprime.dat`
and `-ggprime.dat` carry, and the one scalar MEQ's formulation hides -- `psi_ax`
-- comes out of a two-run fixed point that is EXACT in one step.  See
*The amplitude closure* below.

Every convention below was checked against CHEASE's source AND against a round
trip through the built binary; `README.md` records which evidence is which.


THE FIVE THINGS THE CONVERSION HAS TO GET RIGHT
-----------------------------------------------

**1. The flux label.**  CHEASE's profile abscissa in an `EXPEQ` file is

    s = sqrt( psi_N ),   psi_N = 1 - psi/psi_axis,   s = 0 on axis, 1 on Gamma

(`src-f90/ppspln.f90:76-79` converts the internal psi to exactly this before
interpolating, and `src-f90/premap.f90:62` builds the iso-mesh from its square).
MEQ's `Psi` is 1 on the axis and 0 on Gamma, so

    s = sqrt( 1 - Psi )

and the table has to be REVERSED as well as transformed, since MEQ's tables
ascend in `Psi` and CHEASE wants ascending `s`.

**2. The sign of psi.**  CHEASE's internal psi is zero on the boundary and
NEGATIVE on the axis, increasing outward (COCOS 2; `src-f90/norept.f90:245-258`
for the gauge, and the shipped EQDSK writes `SIMAG < 0`, `SIBRY = 0`).  MEQ's is
zero on Gamma and POSITIVE on the axis.  So

    psi_CHEASE = - psi_MEQ

and every `d/dpsi` flips with it: the `p'` and `TT'` handed to CHEASE are
NEGATIVE where MEQ's tables are positive.  This is the failure `CLAUDE.md`
warns about -- a wrong sign convention converges at the right rate to the wrong
function -- and it is why `--check` prints the reconstructed axis values beside
the tables they came from.

**3. The normalisation.**  A plain `EXPEQ` carries CHEASE units and the reader
applies NO conversion factors of its own (`src-f90/iodisk.f90:369-518` is
straight `READ`s).  The factors, read out of the EQDSK branch where CHEASE does
state them (`iodisk.f90:1103-1109`) and confirmed against the run's own output
to eight figures:

    R, Z    / R0EXP              psi   / ( R0EXP^2 B0EXP )
    p       * mu0 / B0EXP^2      p'    * mu0 R0EXP^2 / B0EXP
    T = g   / ( R0EXP B0EXP )    TT'   / B0EXP
    I       * mu0 / ( R0EXP B0EXP )

**4. The two scalars Grad-Shafranov does not fix.**  `g_at_gamma` and
`p_at_gamma`, now in every `<stem>-meta.json`.  `g_at_gamma` enters as the
UNITS choice: with `NTMF0 = 0` CHEASE anchors `T = 1` at the edge
(`src-f90/isofun.f90:327-331` sets `TMF(KN) = 0.5` and takes `sqrt( 2 x )`), so
`T_SI( Gamma ) = R0EXP * B0EXP` identically, and choosing

    B0EXP = g_at_gamma / R0EXP

is what makes CHEASE's `F` and `q` physical.  It does NOT touch psi: the
Grad-Shafranov operator sees `g` only through `TT'`.  `p_at_gamma` goes in as
`PREDGE`, record 3, and likewise does not touch psi.

**5. What CHEASE would otherwise do to the data.**  `TENSPROF` and `TENSBND`
both default to **-0.3**, not zero (`src-f90/preset.f90:232-233`), and at that
value CHEASE runs an `interpos` smoothing spline over BOTH input profiles
(`iodisk.f90:490-495, 530-535`) and a periodic smoothing fit over the boundary
that also DOUBLES the point count (`iodisk.f90:409-414`).  Nothing says so below
`NVERBOSE = 3`.  Both are set to zero here.  `NCSCAL` defaults to 2, which
rescales both profiles to hit `CURRT`; `NCSCAL = 4` is the no-rescale value
(`norept.f90`'s chain ends at `:224` with no branch for it) and is what this
conversion needs, since the amplitude is the answer rather than an input.


THE AMPLITUDE CLOSURE, WHICH IS THE ONE PIECE OF REAL ARITHMETIC
----------------------------------------------------------------

MEQ and CHEASE are handed DIFFERENT functions, and the difference is exactly
`psi_ax`.

    MEQ      Delta* psi = -[ mu0 r^2 P( Psi ) + G( Psi ) ] / psi_ax,
             Psi = psi / psi_ax,  psi_ax = max psi  (an unknown, bordered)

    CHEASE   Delta* psi = -mu0 r^2 p'( Flux ) - TT'( Flux ),
             Flux = ( psi_edge - psi ) / ( psi_edge - psi_axis )

`Flux` is MEQ's `Psi` exactly.  Matching the two needs `p' = -P/A` and
`TT' = -G/A` with `A` the axis flux of the answer -- so a trial value has to be
supplied and then made self-consistent.

**It is a one-step fixed point, because the map is exactly reciprocal.**  Scaling
both profiles by `lambda` scales `psi` by `lambda` and leaves `Flux` invariant,
so with a trial `A0` CHEASE returns `A_C( A0 ) = K / A0` for a constant `K`.
Self-consistency is `A0 = A_C( A0 )`, i.e.

    A = sqrt( A0 * A_C( A0 ) )

for ANY `A0`.  So two runs suffice: one to measure, one to answer -- and the
second run's own `A_C` reproducing its `A0` is a free check that the reciprocal
law holds, which `cheaserun.py` asserts rather than assumes.

**The second run is not optional even though psi could be rescaled by hand**,
because `T` does not obey the same law: `T^2 = T^2( edge ) + 2 int TT' dpsi`
picks up `lambda^2` where `T^2( edge )` is fixed, so `q`, `F` and `beta` from
the first run are wrong while psi is merely scaled.  Cheaper to run it twice.

Nothing in that closure reads MEQ's answer or the reference.  `A` is therefore
CHEASE's OWN independent measurement of `psi_ax`, and comparing it with the
reference's is a result rather than an input.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "desc-benchmark"))
sys.path.insert(0, os.path.join(HERE, "..", "freegs4e-benchmark"))

import convert_desc as desc                                    # noqa: E402
from convert_desc import CASES, MEQ, REFDIR, MU0, load_case     # noqa: E402
from mxh import mxh_boundary                                    # noqa: E402

#: Where the built binary lives, defaulted BESIDE the MEQ checkout in the same
#: way `race_desc.py` defaults `DESC_DIR` and `MFEM_DIR` -- so a tree cloned
#: anywhere finds it without an absolute path being baked in, and `CHEASE_EXE`
#: overrides it.  `README.md` has the build commands.
CHEASE = os.environ.get(
	"CHEASE_EXE",
	os.path.join(os.path.dirname(MEQ), "chease", "src-f90", "chease"))


def gamma(meta, n=512):
	"""Gamma itself, counterclockwise from the outboard midplane, in metres.

	The SAME curve MEQ solves inside -- `[boundary.shape] Type = "mxh"` -- so
	this is exact data and not a fit to anything.  CHEASE splines it and then
	represents it on its own `NT` poloidal mesh, which is the one place the two
	codes are handed geometry at different resolutions; `--check` measures it.
	"""
	return mxh_boundary(meta["R0"], meta["Z0"], meta["minor_radius"],
	                    meta["kappa"], np.array(meta["cos"]),
	                    np.array(meta["sin"]), n)


def units(meta):
	"""( R0EXP, B0EXP ) -- CHEASE's two normalising constants, in SI.

	`R0EXP` is the GEOMETRIC centre of the boundary, ( Rmax + Rmin )/2, and not
	MXH's `R0` parameter.  **The two are nearly equal and the reason to take
	the geometric one is not their size**: measured by `--check`, they differ
	by 1.1e-15 m on the circular case and by 3.0e-08 m on DIII-D, which is the
	sampling of the extremum rather than a real gap.  It has to be the
	geometric one because `NRSCAL = 0` normalises the boundary so
	that its geometric major radius is 1 in CHEASE units
	(`doc/Chease_Help:176`, and the run reports `RGEOM = 1.0000000`), so a
	different `R0EXP` would put the plasma somewhere else in metres.

	`B0EXP` follows from `g_at_gamma`, per point 4 of the module docstring.
	"""
	R, _ = gamma(meta, 4096)
	R0 = 0.5 * (R.max() + R.min())
	return R0, meta["g_at_gamma"] / R0


def profiles(case, psi_ax, n=257):
	"""( s, p', TT' ) in CHEASE units on an even mesh in s = sqrt( 1 - Psi ).

	Even in `s` rather than in `Psi` because that is the mesh CHEASE builds for
	itself when it has to invent one (`iodisk.f90:1042-1046` writes
	`FCSM( i ) = sqrt( ( i - 1 )/( n - 1 ) )`, i.e. even in `s`), and because
	`s` is the variable its spline interpolates in -- so an even-in-`s` table is
	the one whose interpolation error is uniform.

	The tables are read as ( Psi, value ) pairs and resampled by the same cubic
	spline `convert_desc` uses, for the reason recorded there: the third column
	is a derivative the spline reproduces, and reading it back would
	double-count the smoothing rather than improve on it.
	"""
	from scipy.interpolate import CubicSpline
	R0, B0 = units(case["meta"])
	s = np.linspace(0.0, 1.0, n)
	Psi = 1.0 - s ** 2                       # 1 on axis (s = 0), 0 on Gamma

	psi_p, dp = case["pprime"]
	psi_g, gg = case["ggprime"]
	P = CubicSpline(psi_p, dp)(Psi)          # dp/dPsi   [ Pa per unit Psi ]
	G = CubicSpline(psi_g, gg)(Psi)          # g dg/dPsi [ T^2 m^2 per unit ]

	# MEQ's tables are d/dPsi; CHEASE wants d/dpsi in its OWN gauge, which is
	# the negative of MEQ's.  One factor of 1/psi_ax and one sign, and then the
	# CHEASE normalisation.
	pprime_si = -P / psi_ax                  # dp/dpsi_CHEASE   [ Pa / (Wb/rad) ]
	ttprime_si = -G / psi_ax                 # T dT/dpsi_CHEASE [ T^2 m^2 / (Wb/rad) ]
	return s, pprime_si * MU0 * R0 ** 2 / B0, ttprime_si / B0


def expeq_text(case, psi_ax, n_boundary=2048, n_profile=257):
	"""The `EXPEQ` file: boundary and both profiles, in CHEASE units.

	The layout is `src-f90/iodisk.f90:369-518`, and every read there is
	list-directed, so the whitespace below is free.  Record by record:

	    ASPCT        inverse aspect ratio ( Rmax - Rmin )/( Rmax + Rmin )
	    RZ0C         Z of the geometric centre, normalised
	    PREDGE       p on Gamma, CHEASE units
	    NBPS         number of boundary points
	    R, Z         NBPS pairs, normalised by R0EXP
	    NPPF1 NPPFUN  points in the profile tables, and 4 = "p' is given"
	    NSTTP         1 = "the second table is T dT/dpsi"
	    s             NPPF1 abscissae
	    p'            NPPF1 values
	    TT'           NPPF1 values

	ASPCT and RZ0C are RECOMPUTED by CHEASE from the boundary and silently
	replaced if they disagree by more than 20% / 0.05 * ASPCT respectively
	(`iodisk.f90:394-406`), so they are diagnostics rather than inputs -- but
	they are written correctly here so that the override never fires and the
	file means what it says.
	"""
	meta = case["meta"]
	R0, B0 = units(meta)
	R, Z = gamma(meta, n_boundary)
	aspct = (R.max() - R.min()) / (R.max() + R.min())
	rz0c = 0.5 * (Z.min() + Z.max()) / R0
	predge = meta["p_at_gamma"] * MU0 / B0 ** 2

	s, pp, tt = profiles(case, psi_ax, n_profile)

	out = [f"{aspct:24.16e}", f"{rz0c:24.16e}", f"{predge:24.16e}",
	       f"{len(R):6d}"]
	out += [f"{r / R0:24.16e}{z / R0:24.16e}" for r, z in zip(R, Z)]
	out.append(f"{len(s):6d}{4:6d}")          # NPPF1, NPPFUN = 4 -> p' given
	out.append(f"{1:6d}")                     # NSTTP = 1 -> TT' given
	for block in (s, pp, tt):
		out += [f"{v:24.16e}" for v in block]
	out.append("")
	return "\n".join(out)


def namelist_text(case, psi_ax, ns=40, nt=40, npsi=100, nchi=100, niso=100,
                  nrbox=513, nzbox=513, box=None, nideal=6):
	"""The `chease_namelist`, and every key in it is there for a reason.

	**THE FIRST FOUR LINES ARE NOT COMMENTS TO CHEASE.**
	`src-f90/chease_prog.f90:133-141` reads them into `LABEL1..LABEL4` and only
	then rewinds and reads the namelist, so a file that opens with `&EQDATA`
	loses its first four namelist lines.  They are titles, and they are used.

	The settings that make this the MEQ problem rather than a CHEASE one:

	    NSURF  = 6, NEQDSK = 0   boundary point by point from EXPEQ
	    NPPFUN = 4               the first table is p', not p
	    NFUNC  = 4, NSTTP = 1    the second is tabulated, and it is TT'
	    NCSCAL = 4               DO NOT RESCALE.  1, 2 and 3 all renormalise
	                             both profiles to hit a target q or current
	    NRSCAL = 0               do not rescale lengths to put the MAGNETIC
	                             axis at 1; the geometric one is already there
	    NTMF0  = 0               anchor T at the edge, which is what makes
	                             B0EXP mean g_at_gamma / R0EXP
	    NBLOPT = 0, NBSOPT = 0   do not rewrite p' for ballooning or bootstrap
	    CPRESS = 1, CFNRESS = 1  no one-shot multipliers on either table
	    TENSPROF = 0, TENSBND = 0  do not smooth the inputs.  Both default to
	                             -0.3 and the smoothing is silent

	`NCSCAL = 4` must be paired with `NSTTP = 1` and not with 3 or 4:
	`src-f90/guess.f90:424-445` has no `NCSCAL = 4` branch, so `SCALE` keeps its
	module initialiser of 0 and the current profile is multiplied by ZERO.
	Latent, silent, and avoided here rather than worked around.
	"""
	# **NRBOX AND NZBOX MUST STAY BELOW 1000, AND CHEASE DOES NOT SAY SO.**
	# The EQDSK header is written `(A28,I2.2,A10,A8,3I4)`
	# (`src-f90/iodisk.f90:2092`, and `9380 FORMAT(A40,A8,3I4)` at `:2667` for
	# the `_POS` file), so a four-digit `NRBOX` fills its field and the three
	# integers run together: `... UNITS20260920   310251025`.  CHEASE exits 0
	# and writes a G-EQDSK no standard reader can parse -- `freegs4e`'s dies on
	# `int( 'UNITS20260920' )`.  Refused here rather than met downstream.
	if max(nrbox, nzbox) > 999:
		raise ValueError(
			f"NRBOX = {nrbox}, NZBOX = {nzbox}: CHEASE's EQDSK header format "
			f"is 3I4 and cannot separate four-digit box sizes")
	R0, B0 = units(case["meta"])
	lines = [f"*** CHEASE for MEQ case {case['stem']}",
	         f"*** written by tools/chease-benchmark/convert_chease.py",
	         f"*** trial psi_ax = {psi_ax:.16e} Wb/rad, see the amplitude closure",
	         "***",
	         "&EQDATA"]
	kv = [("NSURF", 6), ("NEQDSK", 0), ("NPPFUN", 4), ("NFUNC", 4),
	      ("NSTTP", 1), ("NFUNRHO", 0),
	      ("NCSCAL", 4), ("NRSCAL", 0), ("NTMF0", 0),
	      ("NBLOPT", 0), ("NBSOPT", 0), ("NBSEXPQ", 0), ("NOPT", 0),
	      ("NIDEAL", nideal), ("NPROPT", 1), ("NPLOT", 0), ("NSMOOTH", 1),
	      ("NS", ns), ("NT", nt), ("NPSI", npsi), ("NCHI", nchi),
	      ("NISO", niso), ("NRBOX", nrbox), ("NZBOX", nzbox),
	      ("NEQDXTPO", 1), ("NVERBOSE", 1)]
	lines += [f" {k} = {v}," for k, v in kv]
	fkv = [("R0EXP", R0), ("B0EXP", B0), ("CPRESS", 1.0), ("CFNRESS", 1.0),
	       ("PSISCL", 1.0), ("TENSPROF", 0.0), ("TENSBND", 0.0),
	       ("SIGNIPXP", 1.0), ("SIGNB0XP", 1.0),
	       ("EPSLON", 1.0e-10), ("RELAX", 0.0), ("GAMMA", 1.6666666666666667)]
	if box is not None:
		rmin, rmax, zmin, zmax = box
		fkv += [("RBOXLFT", rmin), ("RBOXLEN", rmax - rmin),
		        ("ZBOXMID", 0.5 * (zmin + zmax)), ("ZBOXLEN", zmax - zmin)]
	# `float()` before `repr()` because numpy 2 prints a float64 as
	# `np.float64( 1.011... )`, which Fortran reads as a namelist object
	# name and dies with "Cannot match namelist object name np.float64".
	lines += [f" {k} = {float(v)!r}," for k, v in fkv]
	lines += ["/", ""]
	return "\n".join(lines)


def read_eqdsk(path):
	"""CHEASE's G-EQDSK, through freegs4e's own reader.

	**The file is a G-EQDSK with a non-standard tail**: `iodisk.f90:2220` calls
	`OUTMKSA` to append the same human-readable summary table that goes on
	`EXPEQ.OUT`.  A reader that stops after the limiter block -- which
	`freegs4e._geqdsk.read` does -- is unaffected, and one that keeps reading is
	not.  Recorded because it is the obvious way for a second consumer to break.
	"""
	from freegs4e import _geqdsk
	with open(path) as handle:
		return _geqdsk.read(handle)


def inside_gamma(R, Z, meta, n=4096):
	"""Mask of ( R, Z ) grid nodes strictly inside Gamma, even-odd crossings.

	CHEASE's EQDSK carries psi everywhere in its box, and OUTSIDE the plasma
	that psi is an EXTRAPOLATION (`psibox.f90:388-462`, linear by default) with
	a floor clamp at 1e-5 (`psibox.f90:468`).  It is not a solution of anything
	and must not reach an error norm, so the exterior is written as NaN and the
	two comparison routines drop it on `isfinite` exactly as they drop MEQ's
	own band.
	"""
	gr, gz = gamma(meta, n)
	RR, ZZ = np.meshgrid(np.asarray(R, float), np.asarray(Z, float),
	                     indexing="ij")
	inside = np.zeros(RR.shape, bool)
	r1, z1 = gr, gz
	r2, z2 = np.roll(gr, -1), np.roll(gz, -1)
	for a, b, c, d in zip(r1, z1, r2, z2):
		crosses = ((b > ZZ) != (d > ZZ))
		with np.errstate(divide="ignore", invalid="ignore"):
			at = a + (ZZ - b) * (c - a) / (d - b)
		inside ^= crosses & (RR < at)
	return inside


def eqdsk_to_grid(path, meta, mask=True):
	"""CHEASE's answer in the layout `compare.py` and `compare_desc.py` read.

	The reference's own layout: 1-D `R` and `Z`, `psi[ iR, iZ ]`, in FREEGS4E's
	gauge.  Two transformations and both are recorded in the module docstring:
	CHEASE's psi is the negative of MEQ's, and MEQ's gauge is the reference's
	shifted by `psi_surface`.

	    psi_ref = -psi_CHEASE + psi_surface

	`psi_CHEASE` is already zero on Gamma exactly (`iodisk.f90:2110` writes a
	literal zero for `SIBRY`), so the subtraction has nothing to do and the
	whole gauge is one negation plus one shift.
	"""
	g = read_eqdsk(path)
	R = g["rleft"] + g["rdim"] * np.linspace(0.0, 1.0, g["nx"])
	Z = g["zmid"] + g["zdim"] * (np.linspace(0.0, 1.0, g["ny"]) - 0.5)
	psi = -np.array(g["psi"], float) + float(meta["psi_surface"])
	if mask:
		psi = np.where(inside_gamma(R, Z, meta), psi, np.nan)
	return dict(R=R, Z=Z, psi=psi, eqdsk=g)


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(CASES))
	ap.add_argument("--psi-ax", type=float, default=1.0,
	                help="the TRIAL axis flux; any value is admissible, see "
	                     "the amplitude closure")
	ap.add_argument("--check", action="store_true",
	                help="audit the conversion and print what it is worth")
	ap.add_argument("--out", default="", help="write the EXPEQ here")
	args = ap.parse_args()

	case = load_case(args.stem)
	meta = case["meta"]
	R0, B0 = units(meta)
	text = expeq_text(case, args.psi_ax)
	if args.out:
		open(args.out, "w").write(text)
		print(f"  wrote {args.out}  ({len(text.splitlines())} records)")

	if not args.check:
		return 0

	print(f"  {args.stem}")
	print(f"    R0EXP  {R0:.10f} m      the boundary's GEOMETRIC centre")
	print(f"    MXH R0 {meta['R0']:.10f} m      the shape parameter, which is "
	      f"{abs(R0 - meta['R0']):.2e} away")
	print(f"    B0EXP  {B0:.10f} T      = g_at_gamma / R0EXP")
	print(f"    T( Gamma ) = R0EXP B0EXP = {R0 * B0:.10f}  against "
	      f"g_at_gamma {meta['g_at_gamma']:.10f}")
	print(f"    PREDGE {meta['p_at_gamma'] * MU0 / B0 ** 2:.6e}   "
	      f"= p_at_gamma {meta['p_at_gamma']:.6f} Pa")

	# THE PROFILES, AT BOTH ENDS, AGAINST THE META FILE'S OWN TWO NUMBERS.
	# This is the check that the reversal, the sign and the normalisation are
	# all right at once: `pprime_at_gamma` and `ggprime_at_gamma` are in the
	# meta file as d/dPsi, and s = 1 is Gamma.
	#
	# **THE META SCALAR AND THE TABLE ARE IN DIFFERENT UNITS AND DIFFER BY
	# `psi_ax`.**  `pprime_at_gamma` is `dp/dpsi` in Pa per Wb/rad -- the toml
	# says so in its units string -- and the `.dat` table's value column is
	# `dp/dPsi`, which the table's own header describes as carrying "ONE factor
	# of dpsi/dPsi = psi_ax".  Comparing them without the factor reads a clean
	# conversion as 15.7x wrong on this case, which is what happened here
	# first.  The reference's `psi_ax` is used ONLY in this print.
	s, pp, tt = profiles(case, 1.0, 257)
	ref_ax = float(meta["psi_axis"]) - float(meta["psi_surface"])
	print(f"    at s = 1 ( Gamma ):  p' {pp[-1] / (MU0 * R0 ** 2 / B0):+.10e}"
	      f"   against -psi_ax * pprime_at_gamma "
	      f"{-ref_ax * meta['pprime_at_gamma']:+.10e}")
	print(f"                        TT' {tt[-1] * B0:+.10e}"
	      f"   against -psi_ax * ggprime_at_gamma "
	      f"{-ref_ax * meta['ggprime_at_gamma']:+.10e}")
	print(f"    at s = 0 ( axis   ):  p' {pp[0] / (MU0 * R0 ** 2 / B0):+.6e}"
	      f"   TT' {tt[0] * B0:+.6e}   ( both must be NEGATIVE )")

	# THE BOUNDARY, AGAINST THE ONE THING THE TWO CODES SHARE EXACTLY.
	R, Z = gamma(meta, 2048)
	area = 0.5 * np.abs(np.sum(R * np.roll(Z, -1) - np.roll(R, -1) * Z))
	Rf, Zf = gamma(meta, 8192)
	area_f = 0.5 * np.abs(np.sum(Rf * np.roll(Zf, -1) - np.roll(Rf, -1) * Zf))
	print(f"    Gamma at {len(R)} points encloses {area:.10e} m^2, "
	      f"{abs(area - area_f) / area_f:.3e} below the exact curve")
	return 0


if __name__ == "__main__":
	sys.exit(main())
