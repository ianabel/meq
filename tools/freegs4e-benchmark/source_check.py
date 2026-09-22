#!/usr/bin/env python
"""DOES MEQ'S SOURCE, AT THE REFERENCE'S OWN STATE, REPRODUCE THE REFERENCE'S
OWN Jtor?

NO SOLVER RUNS HERE AND NO MESH EXISTS -- the same shape of question CS-0 asks
of the conductor model, and it is worth having for the same reason.  A
disagreement between MEQ and freegs4e is three things at once: the SOURCE MEQ
was given, the DISCRETISATION it solves with, and the ITERATION that got there.
This separates the first from the other two, exactly, by evaluating MEQ's own
expression at freegs4e's own converged psi and comparing against freegs4e's own
current density.

WHAT IT CHECKS, WHICH IS EVERY STEP OF THE CONVERSION AT ONCE:

    * the SENSE of the normalised flux.  MEQ's Psi is 1 on the axis and 0 on
      the boundary and freegs4e's psi_n is the other way up, so the tables are
      written against `Psi = 1 - psi_n`.  Getting that backwards converges to a
      different equilibrium and never complains.
    * the SPAN.  meq::NormalisedMHDSource::f is
      `( mu0 R^2 p'( Psi ) + gg'( Psi ) )/span` with `span = psi_ax - psi_bnd`,
      and the tables are d/dPsi where the reference's arrays are d/dpsi -- so a
      factor of the span appears twice, once in each direction, and dropping
      either is a clean factor this catches.
    * the GEOMETRY factor.  `mu0 R^2` on p' and nothing on gg'.
    * the TABULATION itself, spline knots and all.

WHAT IT CANNOT CHECK: the discretisation, the borders, the boundary condition,
or anything about how the answer is reached.  That is the point -- a red here
is a conversion fault and a green here moves the search somewhere else.

    python3 source_check.py <case.npz> <pprime.dat> <ggprime.dat>

MEASURED on F_diiid_conventional against examples/machine-f-diiid-*.dat:
relative L2 1.49e-03 and a core-current ratio of 0.99965, which is the
tabulation and not a fault.  MEASUREMENTS.md M-152.
"""
import sys

import numpy as np

mu0 = 4.0e-7 * np.pi


def table(path):
	"""Psi and the value, from one of MEQ's two-column profile files."""
	a = np.loadtxt(path)
	return a[:, 0], a[:, 1]


def main():
	if len(sys.argv) != 4:
		raise SystemExit(__doc__)

	npz = np.load(sys.argv[1])
	R1, Z1 = npz["R"], npz["Z"]
	psi = np.asarray(npz["psi"], float)
	psi_ax = float(npz["psi_axis"])
	psi_bn = float(npz["psi_bndry"])
	Jtor = np.asarray(npz["Jtor"], float)
	core = npz["core_mask"] > 0

	Pp, pp = table(sys.argv[2])
	Pg, gg = table(sys.argv[3])

	Rg, _ = np.meshgrid(R1, Z1, indexing="ij")
	span = psi_ax - psi_bn
	Psi = (psi - psi_bn) / span

	# meq::NormalisedMHDSource::f, transcribed.  Delta* psi = -F and
	# F = mu0 R j, so the current density is F/( mu0 R ).
	F = (mu0 * Rg**2 * np.interp(Psi, Pp, pp) + np.interp(Psi, Pg, gg)) / span
	j = F / (mu0 * Rg)

	ref = Jtor[core]
	got = j[core]
	scale = np.abs(ref).max()
	dA = (R1[1] - R1[0]) * (Z1[1] - Z1[0])

	print("grid %d x %d, core cells %d" % (len(R1), len(Z1), core.sum()))
	print("Psi over the core %.4f .. %.4f   (1 on the axis, 0 on the boundary)"
	      % (Psi[core].min(), Psi[core].max()))
	print("table Psi range   p' %.3f..%.3f   gg' %.3f..%.3f"
	      % (Pp.min(), Pp.max(), Pg.min(), Pg.max()))
	# THE TABLES MUST REACH BELOW Psi = 0 AND THEY DO.  A free-boundary run
	# evaluates the source in the vacuum too, where Psi is negative, and
	# meq::SplineProfile extends a table by a CONSTANT past its end knots --
	# so a table stopping at zero puts a plateau in the vacuum.
	if min(Pp.min(), Pg.min()) > 0.0:
		print("  NOTE: neither table reaches Psi < 0, so a free-boundary run "
		      "would read a clamped value in the vacuum")
	print()
	print("reference Jtor over the core %+.6e .. %+.6e" % (ref.min(), ref.max()))
	print("MEQ's source over the core   %+.6e .. %+.6e" % (got.min(), got.max()))
	print()
	print("relative L2 against max | reference |   %.4e"
	      % (np.sqrt(np.mean((got - ref) ** 2)) / scale))
	print("ratio of the two core integrals         %.6f"
	      % (got.sum() / ref.sum()))
	print("core current, reference %.6e A   MEQ's source %.6e A"
	      % (ref.sum() * dA, got.sum() * dA))


if __name__ == "__main__":
	main()
