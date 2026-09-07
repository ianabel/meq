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

MU0 = 4.0e-7 * np.pi


def greens(Rc, Zc, R, Z):
	"""psi at (R, Z) of a unit current filament at (Rc, Zc), Wb/rad."""
	k2 = 4.0 * R * Rc / ((R + Rc) ** 2 + (Z - Zc) ** 2)
	k = np.sqrt(k2)
	return (MU0 / (2.0 * np.pi)) * np.sqrt(R * Rc) * (
		(2.0 / k - k) * ellipk(k2) - (2.0 / k) * ellipe(k2))


def main():
	npz = sys.argv[1]
	mesh_path, gf_path = sys.argv[2], sys.argv[3]
	rho = float(sys.argv[4]) if len(sys.argv) > 4 else 2.60
	n = int(sys.argv[5]) if len(sys.argv) > 5 else 32

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

	# The coil positions are the case's own, and are asserted against the
	# reference's saved coil field by the caller rather than trusted here.
	pos = [(1.75, 0.90), (1.75, -0.90), (0.55, 1.10), (0.55, -1.10)]
	names = ["P1U", "P1L", "P2U", "P2L"]
	currents = d["coil_currents"]
	for (rc, zc), I in zip(pos, currents):
		psi[inner] += I * greens(rc, zc, RR[inner], ZZ[inner])
	for name, I in zip(names, currents):
		print("  coil %-4s %+.6e A at %s" % (name, I, pos[names.index(name)]))

	CH = 200
	for a in range(0, src_I.size, CH):
		b = min(a + CH, src_I.size)
		g = greens(src_R[a:b][None, :], src_Z[a:b][None, :],
		           RR[inner][:, None], ZZ[inner][:, None])
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
	print("reconstruction peak against the reference axis: %.2e relative"
	      % (abs(hi - d["psi_axis"]) / abs(d["psi_axis"])))


if __name__ == "__main__":
	main()
