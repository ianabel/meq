#!/usr/bin/env python
"""THE TWO MEQ PROFILE TABLES FOR A FREE-BOUNDARY CASE, FROM ONE REFERENCE.

WHY THIS IS NOT "REUSE THE SIBLING'S TABLES", WHICH IS WHAT WAS TRIED FIRST AND
IS WRONG FOR EVERY CONSUMER BUT MEQ.  `examples/limited-tokamak-*.dat` were
converted from the 129^2 reference.  Both profiles there are exactly
`amplitude * Psi^2` and so are the 513^2 reference's, so the SHAPES are the
same function -- and MEQ carries a profile scale as an unknown that
`[source] PlasmaCurrent` closes, so one common factor on both tables is
absorbed exactly and MEQ reaches the same equilibrium from either pair.  It is
very tempting to stop there.

THE AMPLITUDES DIFFER BY 6.03%, NOT BY THE 0.22% THE SPAN WOULD EXPLAIN, and
that is the point.  freegs4e's coil currents and its profile amplitudes are
BOTH outputs of its control system, so refining its grid moves the equilibrium
it converges to: `p` on the axis is 46791.95 Pa at 129^2 and 43969.09 Pa at
513^2.  A code that is handed `p( rho )` and `I( rho )` as ABSOLUTE profiles --
DESC, whose toroidal flux depends on `g`; NICE's parametric route; anything
reading `p_at_gamma` -- is then posed a 6% stronger plasma than the reference
it is being compared against, and nothing in MEQ's own run says so, because
MEQ's answer is right either way.

So a shipped free-boundary case gets its own tables, built from its own
reference, and the scale invariance stops being load bearing.

USAGE

    python3 make_freeb_profiles.py <reference.npz> <stem> [psi_min] [psi_max] [n]

which writes `../../examples/<stem>-{pprime,ggprime}.dat`.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.abspath(os.path.join(HERE, "..", "..", "examples"))

#: The abscissa `examples/limited-tokamak-*.dat` carry, kept so the two pairs
#: are differenced row for row.  It runs BELOW zero and ABOVE one on purpose:
#: Psi is the solved normalisation and a Newton iterate leaves [0, 1] freely,
#: so meq::SplineProfile must have real data where it is asked rather than an
#: extrapolation.
PSI_MIN, PSI_MAX, N_ROWS = -0.300, 1.325, 66


def main():
	if len(sys.argv) < 3:
		raise SystemExit("usage: make_freeb_profiles.py <reference.npz> "
		                 "<stem> [psi_min] [psi_max] [n]")
	npz_path, stem = sys.argv[1], sys.argv[2]
	lo = float(sys.argv[3]) if len(sys.argv) > 3 else PSI_MIN
	hi = float(sys.argv[4]) if len(sys.argv) > 4 else PSI_MAX
	rows = int(sys.argv[5]) if len(sys.argv) > 5 else N_ROWS

	d = np.load(npz_path, allow_pickle=True)
	span = float(d["psi_axis"]) - float(d["psi_bndry"])
	psi_n = np.asarray(d["psi_n"], float)

	# THE SHAPE IS ASSERTED, NOT ASSUMED.  This writer only knows how to write
	# `amplitude * Psi^2`, which is what fgsref.py's ( pa, pb ) = ( 1, 2 ) makes
	# -- `( 1 - psi_n^1 )^2` -- and Psi = 1 - psi_n.  A case with any other
	# exponents would be written wrong and silently, so it is refused here.
	out = {}
	for key, name, units, symbol in (
			("pprime", "pprime", "Pa", "dp/dPsi"),
			("ffprime", "ggprime", "T^2 m^2", "g dg/dPsi")):
		raw = np.asarray(d[key], float)
		model = raw[0] * (1.0 - psi_n) ** 2
		residual = float(np.max(np.abs(raw - model)) / abs(raw[0]))
		if residual > 1e-12:
			raise SystemExit(
				"%s's %s is not amplitude*( 1 - psi_n )^2 to 1e-12 -- the "
				"worst residual is %.3e.  This writer emits that one shape "
				"and would produce a wrong table for any other; fgsref.py's "
				"( pa, pb ) / ( fa, fb ) for this case are not ( 1, 2 )."
				% (os.path.basename(npz_path), key, residual))
		out[name] = (float(raw[0]) * span, residual, units, symbol)

	psi = np.linspace(lo, hi, rows)
	for name, (amplitude, residual, units, symbol) in out.items():
		path = os.path.join(EXAMPLES, "%s-%s.dat" % (stem, name))
		value = amplitude * psi ** 2
		derivative = 2.0 * amplitude * psi
		with open(path, "w") as f:
			f.write(
"""# %s( Psi ) = %.15e * Psi^2, for the freegs4e limited tokamak.
#
# BUILT FROM %s, WHICH IS NOT THE GRID examples/limited-tokamak-%s.dat
# CAME FROM, and the two are 6.03%% apart in amplitude rather than the 0.22%%
# the span would explain.  freegs4e's profile amplitudes are an OUTPUT of its
# control system, like its coil currents, so they converge in the grid: p on
# the axis is 46791.95 Pa at 129^2 and 43969.09 Pa at 513^2.  MEQ does not
# care -- the profile scale is an unknown that [source] PlasmaCurrent closes,
# and one common factor on both tables is absorbed exactly -- but DESC, NICE
# and anything else handed ABSOLUTE profiles does.
# tools/freegs4e-benchmark/make_freeb_profiles.py writes this pair.
#
# THE ABSCISSA IS Psi AND SO IS THE DERIVATIVE, WHICH IS THE HALF THAT COSTS
# A FACTOR OF THE SPAN.  meq::NormalisedMHDSource::f evaluates
#     F = scale * ( mu0 R^2 pprime( Psi ) + ggprime( Psi ) ) / span
# so what a table holds is d/dPsi, not d/dpsi.  freegs4e's saved arrays are
# d/dpsi -- its own docstring says otherwise and is wrong, which is the trap
# tools/freegs4e-benchmark/README.md records -- so converting them means
# MULTIPLYING by the span, here %.6e = psi_ax - psi_bnd of the reference.
#
# freegs4e's psi_n runs the other way, psi_n = 1 - Psi, so its shape
# ( 1 - psi_n^1 )^2 is Psi^2 here.  j = 2, which meets FB-4's precondition.
# The shape is asserted against the reference's own array to %.1e before this
# file is written, so it is the reference's profile and not a model of it.
#
# THE RANGE REACHES OUTSIDE [ 0, 1 ] ON PURPOSE.  Psi is the solved
# normalisation, so a Newton iterate leaves the physical interval freely and
# meq::SplineProfile must find data where it is asked rather than extrapolate.
#
#   Psi            f( Psi )   [%s]        f'( Psi )
""" % (symbol, amplitude, os.path.basename(npz_path), name, span,
	   max(residual, 1e-16), units))
			for p, v, dv in zip(psi, value, derivative):
				f.write("%12.6f   %.15e   %.15e\n" % (p, v, dv))
		print("wrote %s" % path)
		print("    %s( Psi ) = %.15e * Psi^2   ( shape residual %.3e )"
		      % (symbol, amplitude, residual))


if __name__ == "__main__":
	main()
