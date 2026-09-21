"""The reference equilibrium reconstructed by Green's functions over MEQ's
half-disc, as an [initialguess] Type = "gridfunction".

WHY THE GUESS IS BUILT FROM THE SOURCE RATHER THAN INTERPOLATED FROM THE ANSWER.
freegs4e's psi exists only on its own 1.6 x 1.6 box, while MEQ's domain for the
free-boundary case is a half-disc of radius 2.6 reaching the axis, so most of
MEQ's mesh has no reference value to interpolate.  What IS available everywhere
is the SOURCE -- Jtor on the core cells and the four coil currents -- and psi is
the Green's-function sum of it, in the same psi -> 0 at infinity gauge both
codes use.  So the guess is as good far from the plasma as it is inside it.

IT IS ALSO A FREE CHECK ON THE WHOLE CONVERSION, because the sum and the PDE
solve share nothing but the source: summed onto a 129^2 grid the reconstruction
reproduces the reference's own psi_axis to 2.5e-04.  The script prints that
comparison every run.

USAGE

    python3 mkexactguess.py <case.npz> <mesh-out> <gf-out> [rho] [n]

rho defaults to 2.6, the radius of the background half-disc
tools/mesh/halfdisc.py cuts for this case, and n to 32, which is the coarsest
grid measured to land the solve on the physical branch -- see below.

WHY n = 32 AND NOT 128.  The guess only has to put Newton in the right basin,
and a Green's-function reconstruction is smooth, so it does not need the
reference's own resolution.  Measured on examples/limited-tokamak.toml, one key
changed, everything else fixed:

    n     mesh + gf   Newton   psi_ax
    128   1.6 MB      7        9.455354e-02
    64    395 kB      8        9.455354e-02
    48    223 kB      7        9.455354e-02
    32     99 kB      7        9.455354e-02

so the answer does not move in ANY printed digit over a sixteenfold change in
the guess's own resolution, while the work moves by one iteration.  That is the
property a guess should have -- an optimisation that moves the answer is a
different problem -- and on a case with more than one solution it is worth
measuring rather than assuming.  n = 32 is what
examples/limited-tokamak-guess.* carries.

WHAT DOES NOT WORK is a cold start: the analytic bump of
examples/free-boundary-halfdisc.toml's kind wanders for 200 iterations around
||r|| = 1.3 and never converges on this problem.  So on a free-boundary machine
case the guess is part of the problem statement, which is the same finding the
three-root sweep in tools/freegs4e-benchmark/README.md reports from the
fixed-boundary side.
"""
import os
import sys

import numpy as np
from scipy.special import ellipe, ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mkguess import write_guess
import conductors

MU0 = 4.0e-7 * np.pi


def greens(Rc, Zc, R, Z, soft=0.0):
	"""psi at (R, Z) of a unit current filament at (Rc, Zc), Wb/rad.

	`soft` REPLACES THE FILAMENT BY THE CELL IT STANDS FOR, and it is not a
	numerical fudge.  Writing the denominator out,

	    ( R + Rc )^2 + ( Z - Zc )^2 = ( R - Rc )^2 + ( Z - Zc )^2 + 4 R Rc

	so 1 - k^2 carries the SEPARATION squared, and this FLOORS that separation
	at `soft` -- it does not add to it.  The distinction is not pedantic: an
	additive soft^2 is felt by every term, and measured on the 129^2 reference
	it took the reconstructed axis flux from 4.6e-04 of the reference to
	2.5e-03, five times WORSE, because the axis sits about half a cell from its
	nearest filament and the whole near field was being damped.  A floor is
	exact everywhere beyond `soft` and finite inside it.

	The floor is the geometric mean distance of a square of side d from itself,
	0.44705 d, which is the textbook self-inductance result: the source really
	is a cell of finite area carrying a uniform current density, and its own
	flux is finite.

	WHY IT IS NEEDED RATHER THAN MERELY TIDY.  Both grids are uniform, so a
	guess node landing on a source filament is systematic and not rare -- and at
	the 513^2 reference it is TOTAL: the guess step is exactly 26 source cells
	in R and 52 in Z, so 20 of 33 R values and 9 of 33 Z values coincide
	exactly, every node inside the core is infinite, and the neighbour fill this
	script used to do instead flattened the peak to 2.6e-02 against a reference
	axis of 9.3e-02 -- 72% WRONG, in the one quantity a guess has to get right.
	At 129^2 the same step is 6.5 cells, nothing aligns, and it never showed.
	"""
	gap2 = (R - Rc) ** 2 + (Z - Zc) ** 2
	if soft > 0.0:
		gap2 = np.maximum(gap2, soft * soft)
	k2 = 4.0 * R * Rc / (gap2 + 4.0 * R * Rc)
	k = np.sqrt(k2)
	return (MU0 / (2.0 * np.pi)) * np.sqrt(R * Rc) * (
		(2.0 / k - k) * ellipk(k2) - (2.0 / k) * ellipe(k2))


def main():
	# --remainder BEFORE the positional reads, so every existing caller is
	# untouched.  See the block comment on `remainder` below for what it is for.
	argv = [a for a in sys.argv[1:] if a != "--remainder"]
	remainder = len(argv) != len(sys.argv) - 1

	npz = argv[0]
	mesh_path, gf_path = argv[1], argv[2]
	rho = float(argv[3]) if len(argv) > 3 else 2.60
	n = int(argv[4]) if len(argv) > 4 else 32

	d = np.load(npz)
	R1, Z1 = d["R"], d["Z"]
	dA = (R1[1] - R1[0]) * (Z1[1] - Z1[0])
	Rg, Zg = np.meshgrid(R1, Z1, indexing="ij")
	core = d["core_mask"] > 0
	src_R, src_Z = Rg[core], Zg[core]
	src_I = d["Jtor"][core] * dA
	print("core filaments %d, total plasma current %.6e A"
	      % (src_I.size, src_I.sum()))

	rs = np.linspace(0.0, rho, n + 1)
	zs = np.linspace(-rho, rho, n + 1)
	RR, ZZ = np.meshgrid(rs, zs, indexing="ij")
	psi = np.zeros_like(RR)

	# r = 0 exactly would make Greens singular in sqrt(R*Rc)*(...); the value
	# there is zero, which is what the array is initialised to.
	inner = RR > 0.0

	# THE COIL POSITIONS COME FROM THE REFERENCE, AND USED NOT TO.
	#
	# This carried its own table of the limited machine's four coils and zipped
	# it against coil_currents POSITIONALLY. That is right for exactly one
	# machine and silently wrong for any other: the currents are saved in
	# tok.coils order, and TestTokamak's labels come out ['P1L','P1U','P2L',
	# 'P2U'] where the limited machine's are ['P1U','P1L','P2U','P2L'] -- so the
	# table would have paired every coil with its opposite number's current and
	# produced a plausible, wrong guess.
	#
	# fgsref.py now saves coil_R and coil_Z beside the labels. The fallback is
	# kept for the two reference files that predate them, and asserts the label
	# order it assumes rather than trusting it.
	names = [str(x) for x in d["coil_labels"]]
	if "coil_R" in d.files:
		flat = conductors.from_npz(d)
	else:
		legacy = ["P1U", "P1L", "P2U", "P2L"]
		if names != legacy:
			raise SystemExit(
				"%s predates coil_R/coil_Z and its coil labels are %s, not the "
				"limited machine's %s. Regenerate it with fgsref.py."
				% (npz, names, legacy))
		legacy_pos = [(1.75, 0.90), (1.75, -0.90), (0.55, 1.10), (0.55, -1.10)]
		flat = [dict(label=n, R=r, Z=z, half_width=0.05, half_height=0.05,
		             current=float(I), kind="filament")
		        for n, (r, z), I in zip(names, legacy_pos, d["coil_currents"])]

	# A CONDUCTOR IS A RECTANGLE AND A TALL ONE IS NOT A POINT.  MEQ's coils
	# carry a uniform current density over a rectangle, and for most of these
	# machines that rectangle is 0.1 x 0.1 m -- small enough that its own centre
	# is a fine stand-in in a guess.  A SOLENOID IS NOT: MAST-U's is 0.04 m wide
	# and 3.16 m tall, carrying more current than everything else in the machine
	# put together, and collapsing it to its midpoint puts all of that at
	# Z = 0.  So a conductor is subdivided until its pieces are roughly square.
	filaments = []
	for c in flat:
		nz = max(1, int(round(c["half_height"]/max(c["half_width"], 1e-9))))
		edges = np.linspace(c["Z"] - c["half_height"], c["Z"] + c["half_height"],
		                    nz + 1)
		centres = 0.5*(edges[:-1] + edges[1:])
		for z in centres:
			filaments.append((c["R"], float(z), c["current"]/nz,
			                  max(c["half_width"], c["half_height"]/nz)))

	'''--remainder: THE PLASMA'S OWN FLUX, WHICH IS WHAT A SUBTRACTED MEQ RUN
	SOLVES FOR.

	Under `[conductors] Model = "subtracted"` or `"filament"` MEQ's unknown is
	psi_p = psi - psi_c, so a guess has to be a guess at psi_p.  Converting a
	total is the obvious route and it is the WRONG one for a FILAMENT machine:
	psi_c diverges logarithmically at each conductor while any stored psi is a
	grid or a finite-element function that cannot, so psi - psi_c at a node
	near a filament is a large negative number where the true psi_p is smooth
	and modest.  Measured on F_diiid_conventional through
	`[initialguess] Content = "total"`: the largest shift is 6.790e-01 Wb/rad
	against a reference psi_axis of 3.759e-01 -- the guess is nearly twice the
	axis flux wrong, at one node, and the bordered Newton walks off.

	Here there is nothing to cancel.  This script builds psi by summing
	Green's functions over the SOURCE, and the source comes in two parts
	already: the conductors, and Jtor on the core cells.  psi_p is the second
	sum alone -- smooth everywhere, including AT a conductor, because no
	conductor is in it.

	The conductors are still read and still printed, because the run they are
	a guess for still has them: what changes is only whether their field is
	added.
	'''
	if not remainder:
		for rc, zc, I, soft in filaments:
			psi[inner] += I * greens(rc, zc, RR[inner], ZZ[inner], soft=soft)
	for c in flat:
		print("  coil %-10s %+.6e A at (%.4f, %+.4f) %.4f x %.4f  %s"
		      % (c["label"], c["current"], c["R"], c["Z"], c["half_width"],
		         c["half_height"], c["kind"]))
	print("  %d conductors as %d filaments%s"
	      % (len(flat), len(filaments),
	         ", NOT SUMMED -- this is psi_p" if remainder else ""))

	# The geometric mean distance of a source CELL from itself.
	cell = 0.44705 * np.sqrt(dA)

	CH = 200
	for a in range(0, src_I.size, CH):
		b = min(a + CH, src_I.size)
		g = greens(src_R[a:b][None, :], src_Z[a:b][None, :],
		           RR[inner][:, None], ZZ[inner][:, None], soft=cell)
		psi[inner] += g @ src_I[a:b]
		print("  %d / %d" % (b, src_I.size), end="\r")
	print()

	# A GUESS NODE LANDING EXACTLY ON A SOURCE FILAMENT GIVES k = 1 AND AN
	# INFINITE GREEN'S FUNCTION.  Both grids are uniform, so coincidences are
	# not rare -- they are systematic.  The filament model is wrong at that
	# range anyway (the real source is a cell of finite area), and this is a
	# guess, so the handful of affected nodes are filled from their finite
	# neighbours rather than modelled.
	bad = ~np.isfinite(psi)
	if bad.any():
		print("filling %d singular node( s )" % bad.sum())
		for _ in range(4):
			if not bad.any():
				break
			fill = np.zeros_like(psi)
			cnt = np.zeros_like(psi)
			good = np.where(bad, 0.0, psi)
			ok = (~bad).astype(float)
			for ax in (0, 1):
				for sh in (-1, 1):
					fill += np.roll(good, sh, ax)
					cnt += np.roll(ok, sh, ax)
			psi = np.where(bad & (cnt > 0), fill / np.maximum(cnt, 1.0), psi)
			bad = ~np.isfinite(psi)

	box = dict(rmin=0.0, rmax=rho, zmin=-rho, zmax=rho)
	lo, hi = write_guess(mesh_path, gf_path, rs, zs, psi.T, box, n=n)
	print("guess psi in [%.6e, %.6e] on a %d x %d grid over [0, %g] x [%g, %g]"
	      % (lo, hi, n + 1, n + 1, rho, -rho, rho))
	print("reference psi_axis %.6e  psi_bndry %.6e"
	      % (d["psi_axis"], d["psi_bndry"]))

	if remainder:
		# AGAINST THE REFERENCE'S OWN DECOMPOSITION, which it saves: `psi` is
		# the total and `plasma_psi` is the plasma's half, so this sum has
		# something exact to be checked against rather than only a peak to
		# compare.  A reconstruction of psi_p cannot be checked against
		# psi_axis at all -- the axis flux is a property of the total.
		if "plasma_psi" in d.files:
			ref = np.asarray(d["plasma_psi"], float)
			scale = max(abs(ref.max()), abs(ref.min()))
			print("reference plasma_psi in [%.6e, %.6e]"
			      % (ref.min(), ref.max()))
			print("this reconstruction of psi_p in [%.6e, %.6e]" % (lo, hi))
			print("peak against the reference's own plasma_psi: %.2e relative"
			      % (abs(max(abs(hi), abs(lo)) - scale)/scale))
		else:
			print("this reconstruction of psi_p in [%.6e, %.6e]" % (lo, hi))
			print("%s carries no plasma_psi, so there is nothing exact to "
			      "check the sum against" % npz)
		print("USE IT WITH [initialguess] Content = \"remainder\", which is "
		      "the default -- this file IS psi_p and must not be converted")
	else:
		print("reconstruction peak against the reference axis: %.2e relative"
		      % (abs(hi - d["psi_axis"]) / abs(d["psi_axis"])))


if __name__ == "__main__":
	main()
