"""A COLD initial guess for a MEQ free-boundary run: the machine, and nothing
that was solved for.

WHY THIS EXISTS AND mkexactguess.py DOES NOT DO.  That script builds the guess
by summing Green's functions over the REFERENCE'S CONVERGED Jtor -- every core
cell of the answer, at the current density the answer ended on.  It is the right
tool for what it was written for: showing that MEQ, handed a reconstruction of
the reference, lands on the same equilibrium.  It is the wrong tool for a race,
because the thing being raced to is already in the starting position.

WHAT IS ADMISSIBLE HERE, AND WHY EACH ITEM IS.

    the conductors and their currents    a machine's operating point.  MEQ is
                                         solving the FORWARD problem: given
                                         these currents, what equilibrium?  The
                                         reference found them with a control
                                         system, which is a different question
                                         and one MEQ is not being asked
    the target plasma current            a scenario parameter, and the same
                                         number [source] PlasmaCurrent
                                         constrains the solve to
    the target X-points and isoflux      what the operator ASKED the shape to
                                         be, in fgsref.py's own case table.
                                         Not where the null ended up
    the profile shapes                   inputs, tabulated from the analytic
                                         forms rather than from the answer

WHAT IS NOT, and each of these is in the reference and is not read:

    psi_axis, psi_bndry                  outputs
    the converged Jtor or core mask      outputs
    the achieved LCFS and its X-points   outputs

SO THE GUESS IS THE VACUUM FIELD PLUS A BLOB.  psi_coils is exact and free --
Green's functions of known currents.  The plasma is a filled ellipse inscribed
in the design shape, carrying Ip with a parabolic current density.  That is
about as much as anybody knows before solving, and it is roughly what
freegs4e's own cold start amounts to.

IT IS NOT MEANT TO BE ACCURATE.  A guess has to put Newton in the right basin
and nothing more; MEASUREMENTS.md and mkexactguess.py's own header both record
that the answer does not move over a sixteenfold change in the guess's
resolution.  What it MUST do is have an O-point roughly where the plasma is,
because on a free boundary the guess chooses which equilibrium is reported.

USAGE

    python3 mkcoldguess.py <case.npz> <mesh-out> <gf-out> [rho] [n]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conductors
from mkexactguess import greens
from mkguess import write_guess


def design_footprint(d, pad=0.15):
    """The plasma's refinement box, from what the run was ASKED for.

    ONE REGION, TWO CONSUMERS.  The cold guess inscribes its current blob in
    this, and make_diverted_case.py refines the mesh inside it; they are the
    same question -- where is the plasma expected to be -- so they are one
    function, and `pad` is the only thing that differs.

    THE DESIGN CONSTRAINS THE BOUNDARY UNEVENLY, and coping with that is the
    whole of this function.  Every isoflux pair is two points the control
    system was told to put on one flux surface, so each is a point of the
    INTENDED boundary, and the target X-points are on it by construction.  For
    TCV, DIII-D and MAST-U that includes an inboard AND an outboard midplane
    point, which pins the radial extent exactly.  For TestTokamak it does not:
    A, B and E name only R = 1.1, the X-point radius, which is nowhere near the
    inboard edge at 0.85.

    So: trust the design where it speaks, and be generous only where it is
    silent.  The wall, where a machine has one, is the other thing an
    experiment knows before it fires; TestTokamak's is r [0.75, 1.80]
    z [+-0.85], which is the envelope those three cases are missing.
    """
    pts = []
    xpt = np.asarray(d["design_xpoints"], dtype=float).reshape(-1, 2)
    pts.extend([(float(r), float(z)) for r, z in xpt])
    iso = np.asarray(d["design_isoflux"], dtype=float).reshape(-1, 4)
    for r1, z1, r2, z2 in iso:
        pts.append((float(r1), float(z1)))
        pts.append((float(r2), float(z2)))
    r0 = float(np.asarray(d["design_R0"]))
    pts = np.asarray(pts, dtype=float)

    zmin, zmax = pts[:, 1].min(), pts[:, 1].max()
    zc = 0.5*(zmin + zmax)
    zhalf = max(0.5*(zmax - zmin), 1e-6)

    # IS EACH RADIAL EDGE NAMED?  A design point on the midplane, inboard or
    # outboard of R0, names one; an X-point at the same radius does not.
    midplane = [(r, z) for r, z in pts if abs(z - zc) < 0.15*zhalf]
    named_in = [r for r, _ in midplane if r < r0]
    named_out = [r for r, _ in midplane if r > r0]

    # THE WALL IS A FALLBACK AND NEVER AN EXPANSION, which is the distinction
    # MAST-U forces.  TestTokamak's wall is a hexagon drawn round its plasma,
    # r [0.75, 1.80], and is exactly the envelope A, B and E's design leaves
    # unsaid.  MAST-U's is the VESSEL, r [0.244, 2.00] z [+-2.169] -- four
    # times the plasma's height, and unioning it would refine the whole
    # divertor.  So the wall answers only where the design is silent.
    wall_r = np.asarray(d["wall_R"], dtype=float)
    wall_z = np.asarray(d["wall_Z"], dtype=float)
    have_wall = wall_r.size > 0

    if named_in:
        rmin = min(named_in)
    elif have_wall:
        rmin = float(wall_r.min())
    else:
        rmin = max(0.02, pts[:, 0].min() - 0.5*zhalf)

    if named_out:
        rmax = max(named_out)
    elif have_wall:
        rmax = float(wall_r.max())
    else:
        rmax = pts[:, 0].max() + 0.5*zhalf

    # The X-point targets bound the height on every machine here, so the design
    # always names z; the wall only clips it.
    if have_wall:
        zmin = max(zmin, float(wall_z.min()))
        zmax = min(zmax, float(wall_z.max()))

    padr = pad*(rmax - rmin)
    padz = pad*(zmax - zmin)
    return (max(0.0, rmin - padr), rmax + padr, zmin - padz, zmax + padz)


MU0 = 4.0e-7*np.pi


def _design_record(npz_path):
	"""The case's profile design, out of the .json beside the .npz.

	THE PROFILE SHAPE DOES NOT FIT IN AN .npz sensibly -- it is four exponents
	and a split, not an array -- so fgsref.py writes it into the record it
	saves as JSON, and this reads it from the file of the same stem.  The two
	are written by one call and are never apart.
	"""
	import json
	side = os.path.splitext(npz_path)[0] + ".json"
	if not os.path.exists(side):
		raise SystemExit("%s has no .json beside it; the profile design lives "
		                 "there" % npz_path)
	return json.load(open(side))["design"]


def plasma_filaments(d, npz_path, cells=61):
	"""The blob: Ip over the design footprint, shaped by the design PROFILES.

	A PARABOLIC BLOB IS NOT PEAKED ENOUGH, and that is not an aesthetic
	complaint.  MEQ's axis border constrains psi_ax at `max psi_h` whenever no
	magnetic axis has been located, and MEQ's domain is a half-disc that
	CONTAINS the conductors where freegs4e's box stops short of them -- so if
	the guess's largest value is inside a coil, the first frozen support and
	the first psi_ax are that coil's.  Measured on case A with j ~ 1 - rho^2:
	the guess peaks at 7.97e-02 inside P1L against 5.56e-02 at the plasma, and
	the run converges in four sweeps to a wall-hugging annulus with no interior
	extremum at all.

	The fix is not a bigger blob but a more PEAKED one, and the peaking is an
	input: fgsref.py's case table names the profile shapes and the split
	between them, so

	    j( R, x ) = R P0 ( 1 - x^pa )^pb + F0 ( 1 - x^fa )^fb / ( mu0 R )

	with x = rho^2 the elliptical flux label and ( P0, F0 ) in the design's own
	ratio, is the current distribution the case ASKED for.  Nothing here is
	read from an equilibrium; the overall size is fixed by Ip, so only the
	shape and the split matter, which is exactly what split_amplitudes() says
	is the physical choice.
	"""
	rmin, rmax, zmin, zmax = design_footprint(d, pad=0.0)
	rc, zc = 0.5*(rmin + rmax), 0.5*(zmin + zmax)
	ra, za = 0.5*(rmax - rmin), 0.5*(zmax - zmin)
	Ip = float(np.asarray(d["design_Ip"]))

	design = _design_record(npz_path)
	r0 = float(design["R0"])
	frac_p = float(design["frac_p"])
	pa, pb = float(design["pa"]), float(design["pb"])
	fa, fb = float(design["fa"]), float(design["fb"])
	p0 = 1.0
	f0 = MU0*r0*r0*p0*(1.0 - frac_p)/frac_p

	r = np.linspace(rmin, rmax, cells)
	z = np.linspace(zmin, zmax, cells)
	RR, ZZ = np.meshgrid(r, z, indexing="ij")
	rho2 = ((RR - rc)/ra)**2 + ((ZZ - zc)/za)**2
	x = np.clip(rho2, 0.0, 1.0)
	shape = (RR*p0*np.power(1.0 - np.power(x, pa), pb)
	         + f0*np.power(1.0 - np.power(x, fa), fb)/(MU0*RR))
	weight = np.where(rho2 < 1.0, shape, 0.0)

	# A DIAMAGNETIC CASE CAN MAKE THIS CHANGE SIGN.  E has frac_p = 1.15, so
	# ff' is negative and the two terms oppose; the blob is still a sensible
	# guess where the sum is positive and is nonsense where it is not.
	weight = np.where(weight > 0.0, weight, 0.0)
	total = weight.sum()
	if total <= 0.0:
		raise SystemExit("the design profile split gives no positive current "
		                 "anywhere in the footprint")

	cell = max((r[1] - r[0]), (z[1] - z[0]))
	inside = weight > 0.0
	return (RR[inside], ZZ[inside], Ip*weight[inside]/total,
	        0.44705*cell, (rc, zc, ra, za))


def main():
	npz = sys.argv[1]
	mesh_path, gf_path = sys.argv[2], sys.argv[3]
	rho = float(sys.argv[4]) if len(sys.argv) > 4 else 2.60
	n = int(sys.argv[5]) if len(sys.argv) > 5 else 32

	d = np.load(npz)

	rs = np.linspace(0.0, rho, n + 1)
	zs = np.linspace(-rho, rho, n + 1)
	RR, ZZ = np.meshgrid(rs, zs, indexing="ij")
	psi = np.zeros_like(RR)
	# r = 0 exactly makes Greens singular in sqrt( R Rc ) ( ... ); the value
	# there is zero, which is what the array already holds.
	inner = RR > 0.0

	# ---- the machine ---------------------------------------------------
	flat = conductors.from_npz(d)
	filaments = []
	for c in flat:
		# A tall conductor is not a point: MAST-U's solenoid is 3.16 m of
		# winding carrying more current than the rest of the machine together,
		# and its midpoint is a poor stand-in for it.
		nz = max(1, int(round(c["half_height"]/max(c["half_width"], 1e-9))))
		edges = np.linspace(c["Z"] - c["half_height"], c["Z"] + c["half_height"],
		                    nz + 1)
		for zmid in 0.5*(edges[:-1] + edges[1:]):
			filaments.append((c["R"], float(zmid), c["current"]/nz,
			                  max(c["half_width"], c["half_height"]/nz)))
	for rc, zc, current, soft in filaments:
		psi[inner] += current*greens(rc, zc, RR[inner], ZZ[inner], soft=soft)
	coil_only = psi.copy()
	print("%d conductors as %d filaments, total %+.6e A"
	      % (len(flat), len(filaments), sum(c["current"] for c in flat)))

	# ---- the blob ------------------------------------------------------
	src_R, src_Z, src_I, soft, ellipse = plasma_filaments(d, npz)
	print("design footprint r [%.3f, %.3f] z [%.3f, %.3f], blob at "
	      "( %.3f, %+.3f ) semi-axes ( %.3f, %.3f ), %d filaments carrying "
	      "%+.6e A"
	      % (ellipse[0] - ellipse[2], ellipse[0] + ellipse[2],
	         ellipse[1] - ellipse[3], ellipse[1] + ellipse[3],
	         ellipse[0], ellipse[1], ellipse[2], ellipse[3],
	         src_I.size, src_I.sum()))

	chunk = 200
	for a in range(0, src_I.size, chunk):
		b = min(a + chunk, src_I.size)
		g = greens(src_R[a:b][None, :], src_Z[a:b][None, :],
		           RR[inner][:, None], ZZ[inner][:, None], soft=soft)
		psi[inner] += g @ src_I[a:b]

	bad = ~np.isfinite(psi)
	if bad.any():
		raise SystemExit("%d guess nodes came out non-finite; the softening "
		                 "floor is meant to make that impossible" % bad.sum())

	box = dict(rmin=0.0, rmax=rho, zmin=-rho, zmax=rho)
	lo, hi = write_guess(mesh_path, gf_path, rs, zs, psi.T, box, n=n)
	print("guess psi in [%.6e, %.6e] on a %d x %d grid over [0, %g] x [%g, %g]"
	      % (lo, hi, n + 1, n + 1, rho, -rho, rho))
	print("the coils alone would give [%.6e, %.6e]; the blob is what puts an "
	      "O-point in the plasma" % (coil_only.min(), coil_only.max()))

	# THE ONE NUMBER A CALLER NEEDS BACK, and it is deliberately the guess's
	# own and not the reference's: [source] PsiAxis is the STARTING VALUE of an
	# unknown, so the best available estimate before solving is the peak of the
	# best available guess.
	#
	# INSIDE THE DESIGN FOOTPRINT, AND THE GLOBAL PEAK IS A COIL'S.  Measured
	# on case A: the guess's largest value is 8.2398e-02 at ( 1.036, -1.097 ),
	# which is the middle of P1L -- a 137 kA conductor whose own O-point beats
	# the plasma's.  Seeding [source] PsiAxis with it makes the frozen support
	# of the first sweep include the coil, and MEQ then converges in seven
	# steps to psi_ax = -8.9e-02 and an X-point at ( 0.62, +0.42 ): a real
	# solution of the problem that was posed, and not the machine.
	#
	# This is the same trap apps/meq.cpp warns about from the other end -- "the
	# axis is inside a coil, so psi_ax is a coil's own O-point" -- met on the
	# INPUT side, where nothing was looking for it.  The footprint is a
	# legitimate input and the plasma is inside it by construction, so the
	# restriction costs nothing.
	rmin, rmax, zmin, zmax = design_footprint(d, pad=0.0)
	region = ((RR >= rmin) & (RR <= rmax) & (ZZ >= zmin) & (ZZ <= zmax))
	if not region.any():
		raise SystemExit("the design footprint contains no guess node; raise n")
	at = np.argmax(np.where(region, psi, -np.inf))
	peak = float(psi.flat[at])
	print("global peak %.6e at ( %.4f, %+.4f )"
	      % (psi.max(), RR.flat[int(np.argmax(psi))],
	         ZZ.flat[int(np.argmax(psi))]))
	print("peak inside the design footprint %.6e at ( %.4f, %+.4f )"
	      % (peak, RR.flat[at], ZZ.flat[at]))
	print("PsiAxis %.9e" % peak)

	# AND HOW FAR OFF IT IS, printed for the reader and used by nothing. A cold
	# guess that happened to be excellent would be worth knowing about; this
	# one is not, which is the point.
	if "psi_axis" in d.files:
		ref = float(d["psi_axis"])
		print("( the reference's psi_axis is %.6e, so this guess's peak is "
		      "%.1f%% off -- it is a basin, not an answer )"
		      % (ref, 100.0*abs(peak - ref)/max(abs(ref), 1e-300)))


if __name__ == "__main__":
	main()
