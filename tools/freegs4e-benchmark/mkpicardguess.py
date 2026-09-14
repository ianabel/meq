"""A cold initial guess that has the RIGHT TOPOLOGY: free-boundary Picard on a
coarse grid, from the machine alone.

WHY THIS EXISTS, AND IT IS A FINDING RATHER THAN A CONVENIENCE.

mkcoldguess.py builds the best guess anybody can write down without iterating:
the conductors at their given currents, plus a current blob shaped by the
design's own profiles and inscribed in the design's own plasma footprint.
Measured on case A that guess peaks at 6.02e-02 against the reference's
8.27e-02 -- 27% low, and in the right place.

MEQ DOES NOT CONVERGE FROM IT.  It converges, in four support sweeps, to a
DIFFERENT equilibrium: psi_bnd 3.13e-02 against the reference's 3.24e-02, which
is right, and psi_ax 1.41e-01 against 8.27e-02, which is 70% high, with a
profile scale of 1.885 making up the difference so that I_p still comes out at
the 2.0e+05 A asked for.  Both are solutions of the same posed problem.  The
free-boundary Grad-Shafranov equation has several, tools/freegs4e-benchmark/
README.md's ramp sweep finds three on a fixed boundary with the physical one
BETWEEN the other two, and which one a Newton reports is decided by where it
starts.  Seeding psi_ax with the reference's converged value changes nothing --
measured, the same 1.406280e-01 to every digit -- so this is the basin and not
the normalisation.

freegs4e CONVERGES ON ALL SEVEN FROM LESS.  Its cold start is
PsiGuessGaussian: a Gaussian bump on the grid, with no coil field in it at all.
What it has that MEQ has not is BLENDED PICARD -- a fixed-point iteration that
recomputes the current from the flux and under-relaxes, which is globally
convergent on this problem where a bordered Newton is locally convergent.

SO THIS IS THAT ITERATION, ON A COARSE GRID, AS AN INITIALISER.  It is
freegs4e's algorithm and it is deliberately not freegs4e: about eighty lines,
Green's functions rather than von Hagenow, a grid that covers the plasma rather
than the machine, and no control system -- the coil currents are given.  Its
cost belongs to MEQ in any comparison and the harness counts it.

WHAT IT READS, AND EVERY ITEM IS AN INPUT:

    the conductors and their currents      the machine's operating point
    the target plasma current              a scenario parameter
    the design X-point and isoflux targets what the shape was ASKED to be
    the design profile shapes and split    inputs, the same ones the TOML's
                                           tables are tabulated from

It reads no psi_axis, no psi_bndry, no converged Jtor and no achieved LCFS.

WHAT IT IS NOT.  It is not accurate and is not meant to be: a coarse grid, a
filament model of its own source cells, and psi_bnd pinned at the DESIGN
X-point rather than at a null it searches for.  What it has to get right is
which equilibrium, which is a topological question and survives a coarse grid.

USAGE

    python3 mkpicardguess.py <case.npz> <mesh-out> <gf-out> [rho] [n] [cells]
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conductors
from mkcoldguess import design_footprint, _design_record, MU0
from mkexactguess import greens
from mkguess import write_guess


def conductor_filaments(flat):
	"""Every conductor as filaments, subdivided until roughly square."""
	out = []
	for c in flat:
		nz = max(1, int(round(c["half_height"]/max(c["half_width"], 1e-9))))
		edges = np.linspace(c["Z"] - c["half_height"], c["Z"] + c["half_height"],
		                    nz + 1)
		for zmid in 0.5*(edges[:-1] + edges[1:]):
			out.append((c["R"], float(zmid), c["current"]/nz,
			            max(c["half_width"], c["half_height"]/nz)))
	return out


def field_of(filaments, R, Z):
	"""psi at ( R, Z ) of a list of ( r, z, I, soft ) filaments."""
	psi = np.zeros_like(R)
	live = R > 0.0
	for rc, zc, current, soft in filaments:
		psi[live] += current*greens(rc, zc, R[live], Z[live], soft=soft)
	return psi


def shapes(design):
	"""p' and gg' against MEQ's Psi, in the design's own ratio.

	fgsref.py's split_amplitudes: P0 = 1 and F0 = mu0 R0^2 ( 1 - f )/f put a
	fraction f of the on-axis current density into the p' term at R = R0. Only
	the RATIO is physical -- the overall size is fixed by I_p below -- which is
	what that function's own docstring says.
	"""
	r0 = float(design["R0"])
	frac_p = float(design["frac_p"])
	p0 = 1.0
	f0 = MU0*r0*r0*p0*(1.0 - frac_p)/frac_p
	pa, pb = float(design["pa"]), float(design["pb"])
	fa, fb = float(design["fa"]), float(design["fb"])

	def current(R, psi_n):
		x = np.clip(psi_n, 0.0, 1.0)
		return (R*p0*np.power(1.0 - np.power(x, pa), pb)
		        + f0*np.power(1.0 - np.power(x, fa), fb)/(MU0*R))

	return current


def flood(candidate, axis):
	"""The connected component of `candidate` containing `axis`.

	A BOX IS NOT A PLASMA, AND THE DIFFERENCE IS A BRANCH.  { Psi > 0 } inside
	the design footprint also picks up whatever lies beyond the null and
	whatever a nearby conductor raises above psi_bnd, and a Picard fed that
	drags its own axis toward the conductor: measured on case A, the axis walks
	from ( 1.30, +0.16 ) to ( 1.20, -0.52 ) over 200 sweeps and lands on
	psi_ax = 1.24e-01 -- which is the SAME non-physical branch MEQ finds from
	the un-iterated guess, at 1.41e-01, rather than the reference's 8.27e-02.

	This is XP-1's flood fill, on a grid instead of a mesh and for the same
	reason: the plasma is the component the magnetic axis is in.
	"""
	from scipy import ndimage
	labels, count = ndimage.label(candidate)
	if count == 0:
		return np.zeros_like(candidate)
	mine = labels[axis]
	if mine == 0:
		return np.zeros_like(candidate)
	return labels == mine


def saddle_near(r, z, psi, target, reach):
	"""The discrete saddle of `psi` nearest `target`, or None.

	WHY A SADDLE SEARCH IS NOT OPTIONAL HERE.  Pinning psi_bnd at the DESIGN
	X-point instead -- which is an input, and is what examples/diverted-
	tokamak.toml does -- makes the iteration collapse, and the collapse is
	instructive: as the current concentrates toward the axis, psi at a point
	NEAR the plasma rises with it, so the core { psi > psi_bnd } shrinks, which
	concentrates the current further.  Measured on case A: the core goes
	310 -> 224 -> 153 -> ... -> 1 cell over 150 sweeps while psi_ax and psi_bnd
	climb together.  A prescribed contact is only stable when the contact is
	genuinely ON the boundary, and a target is not.
	"""
	i, j = np.arange(1, psi.shape[0] - 1), np.arange(1, psi.shape[1] - 1)
	here = psi[1:-1, 1:-1]
	up, down = psi[2:, 1:-1], psi[:-2, 1:-1]
	left, right = psi[1:-1, :-2], psi[1:-1, 2:]
	# A saddle is an extremum of the OPPOSITE sense along the two grid axes.
	sad = (((here > up) & (here > down) & (here < left) & (here < right))
	       | ((here < up) & (here < down) & (here > left) & (here > right)))
	if not sad.any():
		return None
	ii, jj = np.nonzero(sad)
	rr, zz = r[i[ii]], z[j[jj]]
	gap = np.hypot(rr - target[0], zz - target[1])
	within = gap <= reach
	if not within.any():
		return None
	at = int(np.argmin(np.where(within, gap, np.inf)))
	return float(rr[at]), float(zz[at]), float(here[ii[at], jj[at]])


def picard(d, npz_path, cells=97, blend=0.5, sweeps=200, tol=1.0e-7,
           say=print):
	"""The iteration. Returns ( R, Z, psi, psi_ax, psi_bnd, taken )."""
	design = _design_record(npz_path)
	flat = conductors.from_npz(d)
	filaments = conductor_filaments(flat)
	Ip = float(np.asarray(d["design_Ip"]))
	current_of = shapes(design)

	# THE ITERATION'S OWN GRID IS THE PLASMA'S, not the machine's.  What the
	# loop has to resolve is where psi = psi_bnd runs and where the axis is,
	# both inside the design footprint; the field far away is carried exactly by
	# the Green's sum and needs no grid at all.
	rmin, rmax, zmin, zmax = design_footprint(d, pad=0.35)
	rmin = max(rmin, 1.0e-3)
	r = np.linspace(rmin, rmax, cells)
	z = np.linspace(zmin, zmax, cells)
	RR, ZZ = np.meshgrid(r, z, indexing="ij")
	area = (r[1] - r[0])*(z[1] - z[0])
	soft = 0.44705*np.sqrt(area)

	coil_psi = field_of(filaments, RR, ZZ)

	# THE BOUNDARY IS PINNED AT THE DESIGN X-POINT, which is an input and is
	# what makes this robust.  A saddle search on an iterate that is not yet an
	# equilibrium is the thing MEQ's own X-point border struggles with here, and
	# an initialiser has no business doing anything delicate.
	targets = np.asarray(d["design_xpoints"], dtype=float).reshape(-1, 2)
	seed = targets[int(np.argmin(targets[:, 1]))] if targets.size else None

	# The footprint the core is allowed to occupy: outside it lie the private
	# flux region beyond the null and the conductors' own O-points, neither of
	# which is plasma.
	fr0, fr1, fz0, fz1 = design_footprint(d, pad=0.0)
	inside_design = ((RR >= fr0) & (RR <= fr1) & (ZZ >= fz0) & (ZZ <= fz1))
	# TIGHT, AND THE TIGHTNESS IS THE DESIGN SPEAKING.  A machine asked for a
	# null at a place; a saddle half a metre away is a different scenario, and
	# on case A adopting one at z = -0.867 rather than the target's -0.6 lands
	# the iteration on the branch with 1.5x the axis flux.
	reach = 0.15*max(fr1 - fr0, fz1 - fz0)

	# Start from the same blob mkcoldguess.py writes, so the two initialisers
	# differ in the iteration and in nothing else.
	rc, zc = 0.5*(fr0 + fr1), 0.5*(fz0 + fz1)
	ra, za = 0.5*(fr1 - fr0), 0.5*(fz1 - fz0)
	rho2 = ((RR - rc)/ra)**2 + ((ZZ - zc)/za)**2
	seedj = np.where(rho2 < 1.0, current_of(RR, rho2), 0.0)
	seedj = np.where(seedj > 0.0, seedj, 0.0)
	psi = coil_psi + plasma_field(RR, ZZ, RR, ZZ, seedj*area, Ip, soft)

	def bilinear(field, point):
		i = np.clip(np.searchsorted(r, point[0]) - 1, 0, cells - 2)
		j = np.clip(np.searchsorted(z, point[1]) - 1, 0, cells - 2)
		u = (point[0] - r[i])/(r[i + 1] - r[i])
		v = (point[1] - z[j])/(z[j + 1] - z[j])
		return ((1 - u)*(1 - v)*field[i, j] + u*(1 - v)*field[i + 1, j]
		        + (1 - u)*v*field[i, j + 1] + u*v*field[i + 1, j + 1])

	psi_ax = psi_bnd = 0.0
	taken = 0
	for taken in range(1, sweeps + 1):
		masked = np.where(inside_design, psi, -np.inf)
		at = np.unravel_index(int(np.argmax(masked)), psi.shape)
		psi_ax = float(psi[at])
		# THE SADDLE IF THERE IS ONE NEAR THE TARGET, else the target itself.
		# `reach` is a third of the footprint, which is generous enough to
		# follow a null that moves and tight enough not to adopt a coil's.
		found = (saddle_near(r, z, psi, seed, reach) if seed is not None
		         else None)
		psi_bnd = (found[2] if found is not None
		           else (float(bilinear(psi, seed)) if seed is not None
		                 else float(np.min(psi[inside_design]))))
		span = psi_ax - psi_bnd
		if not np.isfinite(span) or span <= 0.0:
			raise SystemExit("the iterate has no positive flux span; the design "
			                 "footprint or the X-point target is wrong")

		normalised = (psi - psi_bnd)/span            # MEQ's Psi: 1 on axis
		core = flood(inside_design & (normalised > 0.0), at)
		if not core.any():
			raise SystemExit("the iterate has no core")
		j = np.zeros_like(psi)
		j[core] = current_of(RR[core], 1.0 - normalised[core])
		j[j < 0.0] = 0.0
		total = (j*area).sum()
		if total <= 0.0:
			raise SystemExit("the iterate carries no current")
		j *= Ip/total

		fresh = coil_psi + plasma_field(RR, ZZ, RR[core], ZZ[core],
		                                (j*area)[core], Ip, soft)
		move = np.max(np.abs(fresh - psi))/max(abs(span), 1e-300)
		psi = (1.0 - blend)*psi + blend*fresh
		if taken <= 6 or taken % 25 == 0:
			say("    %3d  psi_ax %+.6e  psi_bnd %+.6e  core %5d  move %.3e%s"
			    % (taken, psi_ax, psi_bnd, int(core.sum()), move,
			       "" if found is None
			       else "  saddle ( %.3f, %+.3f )" % (found[0], found[1])))
		if move < tol:
			break

	say("  picard: %d sweeps, psi_ax %.6e psi_bnd %.6e at ( %.4f, %+.4f ), "
	    "core %d/%d cells" % (taken, psi_ax, psi_bnd, RR[at], ZZ[at],
	                          int(core.sum()), core.size))
	return r, z, psi, psi_ax, psi_bnd, filaments, RR[core], ZZ[core], (j*area)[core], soft


def plasma_field(R, Z, srcR, srcZ, srcI, Ip, soft, chunk=400):
	psi = np.zeros_like(R)
	live = R > 0.0
	flat = np.flatnonzero(live.ravel())
	Rl, Zl = R.ravel()[flat], Z.ravel()[flat]
	out = np.zeros(flat.size)
	sR, sZ, sI = np.ravel(srcR), np.ravel(srcZ), np.ravel(srcI)
	for a in range(0, sI.size, chunk):
		b = min(a + chunk, sI.size)
		g = greens(sR[a:b][None, :], sZ[a:b][None, :], Rl[:, None], Zl[:, None],
		           soft=soft)
		out += g @ sI[a:b]
	psi.ravel()[flat] = out
	return psi


def main():
	npz = sys.argv[1]
	mesh_path, gf_path = sys.argv[2], sys.argv[3]
	rho = float(sys.argv[4]) if len(sys.argv) > 4 else 2.60
	n = int(sys.argv[5]) if len(sys.argv) > 5 else 32
	cells = int(sys.argv[6]) if len(sys.argv) > 6 else 97

	started = time.perf_counter()
	d = np.load(npz)
	(_, _, _, psi_ax, psi_bnd, filaments,
	 coreR, coreZ, coreI, soft) = picard(d, npz, cells=cells)

	# THE GUESS MEQ READS IS ON THE HALF-DISC, and it is rebuilt from the
	# converged SOURCE rather than interpolated from the iteration's grid --
	# which covers the plasma and not the machine, so most of MEQ's mesh is
	# outside it.  The Green's sum is exact everywhere, which is the same
	# argument mkexactguess.py's header makes.
	rs = np.linspace(0.0, rho, n + 1)
	zs = np.linspace(-rho, rho, n + 1)
	RR, ZZ = np.meshgrid(rs, zs, indexing="ij")
	psi = field_of(filaments, RR, ZZ)
	psi += plasma_field(RR, ZZ, coreR, coreZ, coreI, 0.0, soft)
	if not np.isfinite(psi).all():
		raise SystemExit("the guess came out non-finite")

	box = dict(rmin=0.0, rmax=rho, zmin=-rho, zmax=rho)
	lo, hi = write_guess(mesh_path, gf_path, rs, zs, psi.T, box, n=n)
	print("guess psi in [%.6e, %.6e] on a %d x %d grid over [0, %g] x [%g, %g]"
	      % (lo, hi, n + 1, n + 1, rho, -rho, rho))
	print("initialiser wall %.3f s" % (time.perf_counter() - started))
	print("PsiAxis %.9e" % psi_ax)
	print("PsiBoundary %.9e" % psi_bnd)
	if "psi_axis" in d.files:
		ref, refb = float(d["psi_axis"]), float(d["psi_bndry"])
		print("( the reference is psi_ax %.6e psi_bnd %.6e, so this is %.1f%% "
		      "and %.1f%% off )"
		      % (ref, refb, 100.0*abs(psi_ax - ref)/max(abs(ref), 1e-300),
		         100.0*abs(psi_bnd - refb)/max(abs(refb), 1e-300)))


if __name__ == "__main__":
	main()
