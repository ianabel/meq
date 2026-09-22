#!/usr/bin/env python
"""THE `-meta.json` FOR A **FREE-BOUNDARY** CASE, WHICH make_case.py DOES NOT WRITE.

`make_case.py` builds the six `examples/fixed-*.toml` machines: it fits MXH to
the reference's `psi_n = 0.95` surface and hands MEQ that curve as `Gamma` with
`psi = 0` on it.  The meta it writes beside each is what
`tools/desc-benchmark/convert_desc.py` reads to pose the case for DESC.

A FREE-BOUNDARY CASE HAS NO SUCH CURVE AND STILL NEEDS THE META, which is the
asymmetry this script exists for.  MEQ solves a free-boundary case inside an
artificial `Gamma` in the vacuum and finds the plasma edge itself, so nothing in
the TOML describes the LCFS.  DESC is an inverse code and its unknown IS the
LCFS, so its posing needs one whatever the problem -- as the initial surface for
a free-boundary solve, and, through `FluxGeometry`, as the outermost level of
the flux-label map `Phi( psi_n )`.

SO THE META IS TAKEN AT `level = 1.0` RATHER THAN AT 0.95, AND THAT IS ONLY
LEGITIMATE BECAUSE THE CASE IS LIMITED.  The 0.95 in the fixed cases buys two
things -- profiles analytic across `Gamma`, and a surface with no corner in it.
The second is why `make_case.py` cannot be pointed at a diverted case's
separatrix: it passes through an X-point, and MXH is a truncated Fourier series
that cannot turn a corner.  A LIMITED plasma's LCFS is the flux surface tangent
to the limiter, which is smooth, so the corner objection does not apply and the
LCFS itself is representable.  The first is not bought and is a real cost on
the DESC arm, recorded by `--check`: the profiles are evaluated right up to the
edge.

WHAT IT IS FOR.  `tools/desc-benchmark/` races MEQ against DESC, and until this
existed it could only do so on the fixed-boundary ladder -- because every
free-boundary machine in `examples/` is diverted and DESC's `BoundaryError`
constrains an LCFS it cannot represent.  `examples/limited-tokamak-filament
.toml` is the one that is not.

IT READS THE CASE FILE, AND IT CHECKS ITSELF AGAINST THE REFERENCE FOR DOING
SO.  `CLAUDE.md`'s standing rule is that there is one reader of MEQ's schema and
it is `meq::Configuration`, because a format with two readers that disagree is a
format that grows a silent bug.  This script is a second reader of a narrow
slice -- the two profile file names and the `[[coils]]` blocks -- and what makes
that admissible is that every number it takes out of the TOML is ASSERTED
against the same number in the reference the meta is built from.  A coil moved
in the TOML and not in the reference, or the other way round, is a hard error
here rather than a quiet disagreement downstream.

USAGE

    python3 make_freeb_meta.py <reference.npz> <case.toml> [n_harmonics]

which writes `../../examples/<stem>-meta.json`, the stem being the TOML's.
"""
import json
import os
import sys
import tomllib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mxh import fit_mxh, worst_distance                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.abspath(os.path.join(HERE, "..", "..", "examples"))


def main():
	if len(sys.argv) < 3:
		raise SystemExit("usage: make_freeb_meta.py <reference.npz> "
		                 "<case.toml> [n_harmonics]")
	npz_path, toml_path = sys.argv[1], sys.argv[2]
	harmonics = int(sys.argv[3]) if len(sys.argv) > 3 else 20
	stem = os.path.splitext(os.path.basename(toml_path))[0]

	with open(toml_path, "rb") as handle:
		config = tomllib.load(handle)

	d = np.load(npz_path, allow_pickle=True)
	kind = str(d["boundary_kind"])
	if kind != "limited":
		raise SystemExit(
			"%s is a %s equilibrium.  Its LCFS is a separatrix through an "
			"X-point and has a CORNER, which a truncated Fourier series "
			"cannot turn -- so an MXH fit of it is not a description of it at "
			"any harmonic count.  This script refuses rather than writing a "
			"meta that would be believed." % (npz_path, kind))

	lr = np.asarray(d["lcfs_R"], float)
	lz = np.asarray(d["lcfs_Z"], float)
	fit = fit_mxh(lr, lz, n_harmonics=harmonics)
	shape_error = float(worst_distance(lr, lz, fit))

	psi_axis, psi_bndry = float(d["psi_axis"]), float(d["psi_bndry"])

	# g AND p ON THE SURFACE, AND AT level = 1 BOTH ARE EXACT RATHER THAN
	# INTERPOLATED.  fgsref.py's own conventions record f( 1 ) = fvac and
	# p( 1 ) = 0, so neither needs the spline the fixed cases evaluate at 0.95.
	# Grad-Shafranov fixes neither -- it sees g only through g dg/dpsi and p
	# only through dp/dpsi -- so MEQ's answer does not depend on them and
	# DESC's toroidal flux does.
	g_gamma = float(d["fvac"])
	p_gamma = float(np.asarray(d["pressure"], float)[-1])

	# ---- the machine, out of the case file and checked against the reference
	coils = config.get("coils", [])
	if not coils:
		raise SystemExit("%s has no [[coils]] blocks; a free-boundary case "
		                 "with no conductors has no machine" % toml_path)
	ref_labels = [str(x) for x in d["coil_labels"]]
	ref_current = {lab: float(cur)
	               for lab, cur in zip(ref_labels, d["coil_currents"])}
	if sorted(c["Name"] for c in coils) != sorted(ref_labels):
		raise SystemExit("%s names coils %s and %s has %s"
		                 % (toml_path, sorted(c["Name"] for c in coils),
		                    os.path.basename(npz_path), sorted(ref_labels)))
	for coil in coils:
		want = ref_current[coil["Name"]]
		got = float(coil["Current"])
		# THE TOLERANCE IS THE FILE'S OWN ROUNDING AND NOTHING MORE.  The
		# currents are written to ten figures; anything looser would let a
		# case built against one grid's control-system solution be posed
		# against another's, which is the exact mistake this reference choice
		# exists to avoid.
		if abs(got - want) > 1e-9 * max(abs(want), 1.0):
			raise SystemExit(
				"%s gives %s a current of %.10e and %s solved for %.10e.  "
				"freegs4e's currents are an OUTPUT of its control system, so "
				"they differ between grids -- check the case file is posed "
				"against THIS reference."
				% (toml_path, coil["Name"], got, os.path.basename(npz_path),
				   want))

	source = config.get("source", {})
	conductors = config.get("conductors", {})
	psi_n = np.asarray(d["psi_n"], float)
	meta = dict(
		stem=stem,
		reference=os.path.splitext(os.path.basename(npz_path))[0],
		reference_grid=int(np.asarray(d["R"]).size),
		reference_path=os.path.relpath(os.path.abspath(npz_path), HERE),
		boundary_kind=kind,
		free_boundary=True,
		psi_axis=psi_axis,
		psi_bndry=psi_bndry,
		D=psi_bndry - psi_axis,
		psi_surface=psi_bndry,
		level=1.0,
		shape_error_m=shape_error,
		shape_harmonics=harmonics,
		minor_radius=float(fit["a"]),
		R0=float(fit["R0"]),
		Z0=float(fit["Z0"]),
		kappa=float(fit["kappa"]),
		aspect=float(fit["R0"] / fit["a"]),
		# THE SAME THREE SHAPE SCALARS make_case.py WRITES, AND `asymmetry`
		# IS THE ONE THAT MUST NOT BE GUESSED AT: convert_desc.py decides
		# DESC's `sym` from it, and in MXH it is the COSINE coefficients that
		# break up-down symmetry -- R = R0 + a cos( t + sum c_n cos nt +
		# sum s_n sin nt ) against Z = Z0 + kappa a sin t, where the sines are
		# triangularity and squareness.  Taking max | s_n | instead reads this
		# machine as asymmetric at 5.7e-02, costs DESC a factor of two in
		# unknowns, and is wrong.
		asymmetry=float(np.abs(np.asarray(fit["cos"], float)).max()),
		triangularity=(float(np.sin(fit["sin"][0])) if len(fit["sin"]) else 0.0),
		shaping=(float(np.abs(np.asarray(fit["sin"][1:], float)).max())
		         if len(fit["sin"]) > 1 else 0.0),
		Ip=float(d["Ip"]),
		p_at_gamma=p_gamma,
		g_at_gamma=g_gamma,
		# THE MACHINE, IN THE CASE FILE'S OWN ORDER, so that anything posing
		# this case for another code takes the same conductors MEQ solves
		# with rather than re-deriving them.
		conductor_model=conductors.get("Model", "meshed"),
		coils=[dict(name=c["Name"], R=float(c["CentreR"]),
		            Z=float(c["CentreZ"]),
		            half_width=float(c["HalfWidth"]),
		            half_height=float(c["HalfHeight"]),
		            current=float(c["Current"])) for c in coils],
		# AND THE TWO PROFILE TABLES, BECAUSE A FREE-BOUNDARY CASE NEED NOT
		# OWN THEM.  examples/limited-tokamak-filament.toml reads
		# examples/limited-tokamak-*.dat: the shapes are analytic and the
		# amplitudes are absorbed by [source] PlasmaCurrent, so there is one
		# pair of tables for two case files and the stem does not name them.
		pprime_file=source.get("PPrimeFile"),
		ggprime_file=source.get("GGPrimeFile"),
		plasma_current=float(source["PlasmaCurrent"]),
		units=dict(
			psi_axis="Wb/rad", psi_bndry="Wb/rad", psi_surface="Wb/rad",
			D="Wb/rad", p_at_gamma="Pa", g_at_gamma="T m", Ip="A",
			coil_currents="A, total through each cross-section"),
		note=("level = 1.0: this case is LIMITED, so its LCFS is a smooth "
		      "closed flux surface and MXH describes it.  A diverted case's "
		      "separatrix has an X-point corner and does not, which is why "
		      "make_case.py poses the fixed-boundary ladder at 0.95."),
		cos=[float(x) for x in fit["cos"]],
		sin=[float(x) for x in fit["sin"]],
	)

	out = os.path.join(EXAMPLES, "%s-meta.json" % stem)
	with open(out, "w") as handle:
		json.dump(meta, handle, indent=2)
	print("wrote %s" % out)
	print("  MXH at %d harmonics: worst %.4e m on a minor radius of %.6f "
	      "( %.2e relative )"
	      % (harmonics, shape_error, meta["minor_radius"],
	         shape_error / meta["minor_radius"]))
	print("  R0 %.9f  a %.9f  kappa %.9f  asymmetry %.3e"
	      % (meta["R0"], meta["minor_radius"], meta["kappa"],
	         meta["asymmetry"]))
	print("  psi_axis %.15e  psi_bndry %.15e" % (psi_axis, psi_bndry))
	print("  g( Gamma ) %.9f from fvac, p( Gamma ) %.6e from the profile's "
	      "own last point ( psi_n = %.6f )" % (g_gamma, p_gamma, psi_n[-1]))
	print("  %d conductors, model %r, every current agreeing with %s to 1e-09"
	      % (len(coils), meta["conductor_model"], os.path.basename(npz_path)))
	print("  profiles %s and %s" % (meta["pprime_file"], meta["ggprime_file"]))


if __name__ == "__main__":
	main()
