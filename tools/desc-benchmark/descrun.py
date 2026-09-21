#!/usr/bin/env python
"""ONE DESC SOLVE OF ONE MEQ CASE, TIMED, AND ITS psi( R, z ) ON A COMMON GRID.

The DESC arm of the race.  `convert_desc.py` poses the case; this runs it, times it,
and writes the one thing the two codes can be differenced on -- the poloidal
flux on a fixed ( R, z ) grid, in MEQ's gauge, so that
`tools/freegs4e-benchmark/compare.py` reads it unchanged.

TWO CLOCKS, AND QUOTING ONLY ONE OF THEM IS THE EASIEST WAY TO GET THIS RACE
WRONG.  DESC is JAX, so the first solve of a given SHAPE compiles the residual,
its Jacobian and the whole optimiser loop, and the second solve of that same
shape does not.  Measured on `fixed-h-circular` at M = 10: **14.1 s for the
first and 0.4 s for the second in one process.**  Neither is "the" answer:

  * `cold` -- a fresh process, one solve.  What a user pays to solve one
    equilibrium once, and the number to compare against MEQ's whole driver run.
  * `warm` -- the same shape solved again in the same process.  What a
    parameter scan or an optimiser inner loop pays, and the number to compare
    against MEQ's own warm-started solve leg.

Both are reported on every row.  MEQ has no equivalent split: it compiles at
build time, and its `setup` leg is mesh and assembly rather than compilation.

WHAT IS TIMED.  Everything from "the inputs exist" to "the equilibrium is
solved": surface fit, `Equilibrium` construction and `solve()`.  NOT the
conversion in `convert_desc.py`, which is posing the problem rather than solving it
-- the same reason MEQ's meshing sits outside `race_nke.py`'s clock -- and NOT
the field evaluation below, which is output.  `--breakdown` prints the legs.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import convert_desc as convert                                # noqa: E402


def build_equilibrium(spec, M, L=None, M_grid=None, L_grid=None,
                      boundary_M=None, current_sign=+1.0, profiles="power"):
	"""The DESC object for one converted case at one resolution.

	`boundary_M` defaults to M: the surface is re-fitted at every resolution
	rather than fitted once and truncated, because a surface carried at a
	higher M than the equilibrium is silently projected down by
	`Equilibrium.__init__` and the boundary you measured is not the one solved.
	"""
	from desc.equilibrium import Equilibrium
	from desc.geometry import FourierRZToroidalSurface
	from desc.profiles import PowerSeriesProfile

	L = M if L is None else L
	M_grid = 2 * M if M_grid is None else M_grid
	L_grid = 2 * L if L_grid is None else L_grid
	boundary_M = M if boundary_M is None else boundary_M

	theta, Rb, Zb = convert.desc_boundary_samples(spec["meta"],
	                                              max(512, 8 * boundary_M))
	surface = FourierRZToroidalSurface.from_values(
		np.column_stack([Rb, np.zeros_like(Rb), Zb]), theta,
		M=boundary_M, N=0, NFP=1, sym=spec["up_down_symmetric"])

	if profiles == "power":
		# sym = "auto" rather than False: DESC warns on a profile it cannot see
		# as even, and the zero odd coefficients are exactly what lets it.
		pressure = PowerSeriesProfile(
			params=spec["p_coefficients"],
			modes=np.arange(len(spec["p_coefficients"])))
		current = PowerSeriesProfile(
			params=current_sign * spec["current_coefficients"],
			modes=np.arange(len(spec["current_coefficients"])))
	else:
		from desc.profiles import HermiteSplineProfile
		pressure = HermiteSplineProfile(spec["p_knots"], spec["dp_knots"],
		                                knots=spec["knots"])
		current = HermiteSplineProfile(
			current_sign * spec["current_knots"],
			current_sign * spec["dcurrent_knots"], knots=spec["knots"])

	eq = Equilibrium(surface=surface, Psi=spec["Psi_total"], NFP=1,
	                 L=L, M=M, N=0, L_grid=L_grid, M_grid=M_grid, N_grid=0,
	                 sym=spec["up_down_symmetric"],
	                 pressure=pressure, current=current)
	return eq, (theta, Rb, Zb, surface)


def boundary_residual(surface, theta, Rb, Zb):
	"""Worst distance between the fitted surface and Gamma, at the same theta.

	AT THE SAME PARAMETER, which needs no search -- `mxh.pointwise_error` makes
	the same choice for the same reason.  This is one of the three inexactnesses
	of the posing and it is the one that moves with resolution, so it belongs on
	every row rather than in the conversion's own report.
	"""
	from desc.grid import LinearGrid
	grid = LinearGrid(rho=np.array([1.0]), theta=theta, zeta=np.array([0.0]))
	data = surface.compute(["R", "Z"], grid=grid)
	return float(np.hypot(np.asarray(data["R"]) - Rb,
	                      np.asarray(data["Z"]) - Zb).max())


def solve(eq, ftol=1e-8, xtol=1e-8, gtol=1e-10, maxiter=100, verbose=0):
	started = time.perf_counter()
	eq, result = eq.solve(objective="force", optimizer="lsq-exact", ftol=ftol,
	                      xtol=xtol, gtol=gtol, maxiter=maxiter,
	                      verbose=verbose, copy=False)
	return eq, result, time.perf_counter() - started


def poloidal_flux(eq, n_rho=None):
	"""( rho, psi in MEQ's gauge ) -- zero on Gamma, psi_ax on the axis.

	DESC's `chi` is the poloidal flux normalised by 2 pi, which is MEQ's psi in
	MEQ's units, and it is zero on the AXIS and grows outward.  MEQ's is zero on
	Gamma and grows inward, so the two differ by a reflection about the
	half-span:  psi_MEQ( rho ) = | chi( 1 ) | - | chi( rho ) |.

	The absolute values are not sloppiness.  `chi_r = psi_r * iota` carries the
	sign of iota, which carries the sign the current profile was given, and the
	current's sign is NOT determined by MEQ's input: Grad-Shafranov sees only
	g dg/dpsi, so the mirrored equilibrium solves it equally well.  psi is the
	same function either way, and `--current-sign` asserts that rather than
	assuming it.
	"""
	# A FINE RADIAL GRID, BECAUSE `chi` IS A QUADRATURE AND NOT A FIELD.
	# DESC integrates `chi_r` with `cumulative_simpson` over whatever radial
	# nodes the grid has, and `LinearGrid( L = eq.L_grid )` gives 25 of them at
	# M = 12.  Measured: 6.351762845e-02 on the default grid against
	# 6.351752988e-02 converged, i.e. 1.6e-06 of quadrature sitting in the one
	# number this file reports.  257 points costs milliseconds and removes it.
	from desc.grid import LinearGrid
	grid = LinearGrid(rho=np.linspace(1e-12, 1.0, n_rho or 257),
	                  M=eq.M_grid, N=0, NFP=1)
	data = eq.compute(["rho", "chi"], grid=grid)
	rho = np.asarray(grid.compress(np.asarray(data["rho"])))
	chi = np.abs(np.asarray(grid.compress(np.asarray(data["chi"]))))
	return rho, chi[-1] - chi


def flux_on_grid(eq, R, Z, n_rho=192, n_theta=384):
	"""psi_MEQ on a fixed ( R, z ) grid, from the forward map alone.

	NO COORDINATE INVERSION.  `Equilibrium.map_coordinates` would answer this
	exactly, and it is a JIT'd vmap of a Newton solve per point that returns NaN
	where it fails -- which is every node outside the plasma, i.e. most of a
	rectangle containing it.  The forward map ( rho, theta ) -> ( R, z ) is a
	transform evaluation, it cannot fail, and psi is a function of rho alone, so
	interpolating the point cloud is enough.  `--check-inverse` measures the two
	against each other on a subsample.

	C1 CUBIC RATHER THAN LINEAR, because this interpolation is a FLOOR on every
	error the race reports and a linear one sits at O( h_cloud^2 ) -- about
	8e-07 relative at the default cloud, which is the same size as the numbers
	being resolved at the fine end.  `--cloud` moves it and `--check-cloud`
	measures the movement.
	"""
	from scipy.interpolate import CloughTocher2DInterpolator
	from desc.grid import LinearGrid

	rho_edges = np.linspace(0.0, 1.0, n_rho)
	theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
	grid = LinearGrid(rho=rho_edges, theta=theta, zeta=np.array([0.0]), NFP=1)
	data = eq.compute(["R", "Z", "rho"], grid=grid)
	Rc = np.asarray(data["R"]); Zc = np.asarray(data["Z"])
	rho_c = np.asarray(data["rho"])

	rho_line, psi_line = poloidal_flux(eq)
	psi_c = np.interp(rho_c, rho_line, psi_line)

	interp = CloughTocher2DInterpolator(np.column_stack([Rc, Zc]), psi_c,
	                                    fill_value=np.nan)
	RR, ZZ = np.meshgrid(R, Z, indexing="ij")
	return interp(RR, ZZ), (Rc, Zc, psi_c)


def write_npz(path, R, Z, psi_meq, spec, eq, result, timing, extra):
	"""The DESC answer in the REFERENCE's own file format and gauge.

	`tools/freegs4e-benchmark/compare.py` reads `psi[ iR, iZ ]`, `psi_axis` and
	`psi_bndry` and shifts by the meta file's `psi_surface` to reach MEQ's
	gauge.  Writing this file the same way means the MEQ-against-DESC
	comparison is that script, unchanged, with its band mask and its axis check
	-- rather than a second norm written here that could differ from it
	silently.
	"""
	shift = float(spec["meta"]["psi_surface"])
	psi_axis = float(np.nanmax(psi_meq)) + shift
	# `psi_meq` BESIDE `psi`, AND IT IS NOT REDUNDANT.  `psi` is in the
	# reference's gauge so that compare.py reads this file unchanged;
	# `psi_meq` is the same field in MEQ's, which is what
	# convert_desc.field_from_desc() needs to build a flux label from.  Storing
	# both means the self-consistent loop never has to undo a shift this
	# function applied, which is one transcription of the gauge rather than
	# two.
	np.savez(path, R=R, Z=Z, psi=psi_meq + shift, psi_meq=psi_meq,
	         psi_surface=shift,
	         psi_axis=psi_axis, psi_bndry=float(spec["psi_bndry"]),
	         psi_index_order="psi[iR, iZ]; R is axis 0, Z is axis 1",
	         **extra)


def grid_for(case, grid_from):
	"""The ( R, Z ) the answer is reported on: the reference's, or another
	file's."""
	if grid_from:
		ref = np.load(grid_from, allow_pickle=True)
		return np.array(ref["R"], float), np.array(ref["Z"], float)
	return (np.array(case["ref"]["R"], float),
	        np.array(case["ref"]["Z"], float))


def iterate_map(case, args, field, R, Z):
	"""Re-pose on DESC's OWN field until the posing stops moving.

	**WHY THIS EXISTS: WITHOUT IT THE TWO CODES ARE NOT HANDED THE SAME
	PROBLEM.**  DESC's profiles live against `rho = sqrt( Phi/Phi_edge )`, a
	TOROIDAL flux label, and MEQ's live against normalised POLOIDAL flux.  The
	map between them is a functional of the solution, so posing the case for
	DESC needs an equilibrium -- and whichever one is used leaks into DESC's
	input.  `--from reference` leaks freegs4e's interior surfaces; `--from
	<meq.nc>` leaks MEQ's.  Either way a MEQ-against-DESC number is bounded
	below by the posing rather than by the solvers, and a race is then partly
	a measurement of the third code.

	The fixed point is the way out.  Pose from anything, solve, rebuild the
	map from THAT field, re-pose, re-solve.  What the limit depends on is the
	MXH curve, the two profile tables and the two scalars Grad-Shafranov does
	not fix -- which is exactly the problem MEQ is handed, and nothing else.
	The starting map becomes a starting guess, which is what it always should
	have been.

	**UNTIMED, AND DELIBERATELY.**  This file's header draws the line at "the
	inputs exist": posing is not solving, the same reason MEQ's meshing sits
	outside `race_nke.py`'s clock.  The caller times one final solve on the
	converged posing.  The sweeps cost one DESC solve each and there are two
	or three of them.

	Returns ( field, history ); `history` is one row per sweep and is written
	into the result so a run can be audited for whether it actually converged.
	"""
	history = []
	previous = None
	blended = None
	for sweep in range(1, args.self_consistent + 1):
		spec = convert.build(case, n_levels=args.levels, n_fit=args.terms,
		                     field=field, g_boundary=args.g_boundary,
		                     p_boundary=args.p_boundary)
		eq, _ = build_equilibrium(spec, args.M, L=args.L,
		                          current_sign=args.current_sign,
		                          profiles=args.profiles)
		eq, result, _ = solve(eq, args.ftol, args.xtol, args.gtol,
		                      args.maxiter, args.verbose)
		psi_meq, _ = flux_on_grid(eq, R, Z, n_rho=args.cloud,
		                          n_theta=2 * args.cloud)

		# THE CHANGE IS IN THE FIELD AND NOT IN A SCALAR, because a scalar can
		# sit still while the surfaces the map is built from move.  Over the
		# nodes BOTH sweeps answered for, which is the plasma.
		change = float("inf")
		if previous is not None:
			both = np.isfinite(psi_meq) & np.isfinite(previous)
			if both.any():
				scale = np.max(np.abs(previous[both]))
				change = float(np.sqrt(np.mean(
					(psi_meq[both] - previous[both]) ** 2)) / max(scale, 1e-300))

		row = dict(sweep=sweep, source=spec["map_source"],
		           psi_ax=float(np.nanmax(psi_meq)),
		           Psi_total=float(spec["Psi_total"]),
		           nit=int(result["nit"]), field_change=change)
		history.append(row)
		print(f"    map sweep {sweep}: from {row['source']:28s} "
		      f"Psi {row['Psi_total']:.9e}  psi_ax {row['psi_ax']:.9e}  "
		      f"change {change:.3e}")

		# ---- DAMPED, AND THE UNDAMPED ITERATION DOES NOT CONVERGE --------
		#
		# MEASURED on fixed-h-circular at M = 10, replacing the map outright
		# ( mix = 1 ) over five sweeps:
		#
		#   Psi_total  3.305134474e-01 -> 3.305105922e-01, settling at 1e-07
		#   psi_ax     6.353780e-02, 6.352345e-02, 6.353595e-02, 6.352159e-02,
		#              6.353659e-02
		#
		# -- a PERIOD-2 LIMIT CYCLE of amplitude 2.3e-04 relative, with the
		# field change plateauing at 1e-04 and refusing to fall.  The map
		# itself converges; what oscillates is the axis flux, which is the
		# ordinary marginal instability of a Picard iteration on a
		# self-consistent quantity and is exactly what a mixing parameter
		# exists for.  Half-and-half collapses it.
		#
		# BLENDED ON THE FIELD AND NOT ON THE LABEL, because `psin` is
		# nonlinear in psi -- it carries a 1/psi_ax that moves between sweeps
		# -- so averaging two labels is not the label of an averaged field.
		mixed = psi_meq
		if blended is not None:
			both = np.isfinite(psi_meq) & np.isfinite(blended)
			mixed = np.where(both,
			                 (1.0 - args.self_consistent_mix) * blended
			                 + args.self_consistent_mix * psi_meq,
			                 psi_meq)
		blended = mixed

		previous = psi_meq
		field = convert.field_from_grid(R, Z, mixed, float(case["meta"]["level"]),
		                                source=f"desc self-consistent sweep {sweep}")
		if change < args.self_consistent_tol:
			print(f"    the posing has stopped moving at {change:.3e} "
			      f"< {args.self_consistent_tol:.0e}")
			break
	else:
		if args.self_consistent:
			print(f"    WARNING: {args.self_consistent} sweeps did not reach "
			      f"{args.self_consistent_tol:.0e}; the posing is still moving "
			      f"and the number below is NOT self-consistent")
	return field, history


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", choices=sorted(convert.CASES))
	ap.add_argument("-M", type=int, default=10, help="poloidal/radial resolution")
	ap.add_argument("-L", type=int, default=None)
	ap.add_argument("--refine", type=int, default=4,
	                help="bicubic upsampling of the reference field before "
	                     "contouring; convert_desc.refine_field has why")
	ap.add_argument("--from", dest="source", default="",
	                help="a MEQ .nc whose psi supplies the flux-label map, "
	                     "instead of the freegs4e reference -- see convert.py's "
	                     "field_from_meq for why that is a different posing")
	ap.add_argument("--terms", type=int, default=12,
	                help="even powers of rho in each refitted profile")
	ap.add_argument("--levels", type=int, default=120,
	                help="flux surfaces the conversion integrates over")
	# DESC'S OWN DEFAULTS, AND TIGHTENING THEM WAS MEASURED TO BUY NOTHING.
	# `ftol` is a RELATIVE COST REDUCTION, so it is not the same kind of number
	# as MEQ's NewtonRelativeTolerance and the two cannot be set to "the same"
	# value; what can be done is to put each code past the point where its own
	# answer stops moving, and to show that it did.  On fixed-h-circular at
	# M = 12, `--from` the reference:
	#
	#     ftol, xtol, gtol      nit   psi_ax          vs the reference
	#     1e-2, 1e-6, 1e-6        9   6.351762845e-02   4.203e-04
	#     1e-6, 1e-8, 1e-8       19   6.351764849e-02   4.200e-04
	#     1e-10,1e-10,1e-12      40   6.351764927e-02   4.200e-04
	#
	# Four times the iterations for a 3.2e-07 relative move, while the thing
	# actually limiting that row is 4.2e-04 and lives in the conversion.  So
	# the defaults are DESC's, and a race run on tighter ones would be charging
	# DESC for iterations its own recommended settings do not take.
	ap.add_argument("--ftol", type=float, default=1e-2)
	ap.add_argument("--xtol", type=float, default=1e-6)
	ap.add_argument("--gtol", type=float, default=1e-6)
	ap.add_argument("--maxiter", type=int, default=100)
	ap.add_argument("--current-sign", type=float, default=1.0)
	ap.add_argument("--profiles", choices=("power", "hermite"), default="power",
	                help="how p( rho ) and I( rho ) are represented to DESC; "
	                     "convert_desc.py's knot block has the measurement that "
	                     "makes this a real choice rather than a preference")
	ap.add_argument("--out", default="")
	ap.add_argument("--grid-from", default="",
	                help="an .npz or .nc whose ( R, z ) grid the field is "
	                     "written on; default is the reference's 257^2")
	ap.add_argument("--cloud", type=int, default=192,
	                help="radial points in the forward-map cloud")
	ap.add_argument("--warm", type=int, default=1,
	                help="extra solves of the same shape, for the warm clock")
	ap.add_argument("--no-field", action="store_true")
	ap.add_argument("--verbose", type=int, default=0)
	ap.add_argument("--self-consistent", type=int, default=0, metavar="N",
	                help="iterate the flux-label map on DESC's own field, at "
	                     "most N sweeps; 0 (the default) poses once from "
	                     "--from and leaves that code in the posing")
	ap.add_argument("--self-consistent-tol", type=float, default=1e-6,
	                help="stop when the field moves less than this, relative")
	ap.add_argument("--self-consistent-mix", type=float, default=0.5,
	                metavar="W",
	                help="how much of each new map to take, 0 < W <= 1. W = 1 "
	                     "replaces it outright and does NOT converge -- see "
	                     "iterate_map() for the period-2 cycle that measures")
	ap.add_argument("--g-boundary", type=float, default=None,
	                help="g = R B_phi on Gamma, in T m. Grad-Shafranov sees g "
	                     "only through g dg/dpsi, so MEQ's input fixes it only "
	                     "up to this constant and DESC's toroidal flux needs "
	                     "it. Default: the reference's own fpol")
	ap.add_argument("--p-boundary", type=float, default=None,
	                help="the pressure on Gamma, in Pa. The same asymmetry: "
	                     "MEQ's table is dp/dPsi. Default: the reference's own "
	                     "pressure array")
	ap.add_argument("--json", default="", help="write the row here")
	args = ap.parse_args()

	case = convert.load_case(args.stem)
	level = float(case["meta"]["level"])
	field = (convert.field_from_meq(args.source, level) if args.source
	         else convert.field_from_reference(case["ref"], level, args.refine))
	R, Z = grid_for(case, args.grid_from)

	# THE MAP FIRST, IF IT IS TO BE ITERATED AT ALL.  Untimed: see iterate_map.
	sc_history = []
	if args.self_consistent:
		field, sc_history = iterate_map(case, args, field, R, Z)

	t_convert = time.perf_counter()
	spec = convert.build(case, n_levels=args.levels, n_fit=args.terms,
	                     field=field, g_boundary=args.g_boundary,
	                     p_boundary=args.p_boundary)
	t_convert = time.perf_counter() - t_convert

	# BUILD AND SOLVE TOGETHER ON BOTH CLOCKS.  `Equilibrium.__init__` and the
	# surface fit are themselves JIT'd, so the first one costs seconds and the
	# second costs milliseconds; charging the cold clock for the solve alone
	# and the warm clock for the same thing would hide a leg that is 43% of a
	# first run and 3% of a later one.  Measured on fixed-h-circular at M = 10:
	# build 11.6 s cold against 0.06 s warm.
	t0 = time.perf_counter()
	eq, (theta, Rb, Zb, surface) = build_equilibrium(
		spec, args.M, L=args.L, current_sign=args.current_sign,
		profiles=args.profiles)
	t_build = time.perf_counter() - t0
	eq, result, t_solve = solve(eq, args.ftol, args.xtol, args.gtol,
	                            args.maxiter, args.verbose)
	t_cold = t_build + t_solve

	warm, warm_build = [], []
	for _ in range(max(0, args.warm)):
		t1 = time.perf_counter()
		again, _ = build_equilibrium(spec, args.M, L=args.L,
		                             current_sign=args.current_sign,
		                             profiles=args.profiles)
		t1 = time.perf_counter() - t1
		_, _, seconds = solve(again, args.ftol, args.xtol, args.gtol,
		                      args.maxiter, args.verbose)
		warm.append(t1 + seconds); warm_build.append(t1)

	rho, psi_line = poloidal_flux(eq)
	psi_ax = float(psi_line[0])
	R_axis = float(np.asarray(eq.get_axis().R_n)[0])
	residual = boundary_residual(surface, theta, Rb, Zb)

	row = dict(stem=args.stem, M=args.M, L=int(eq.L), M_grid=int(eq.M_grid),
	           L_grid=int(eq.L_grid), dim_x=int(eq.dim_x),
	           nit=int(result["nit"]), cost=float(result["cost"]),
	           success=bool(result["success"]),
	           cold_seconds=t_cold, build_seconds=t_build,
	           solve_seconds=t_solve,
	           warm_build_seconds=(float(np.median(warm_build)) if warm else None),
	           warm_seconds=(float(np.median(warm)) if warm else None),
	           convert_seconds=t_convert,
	           psi_ax=psi_ax, R_axis=R_axis,
	           psi_ax_reference=spec["psi_ax_reference"],
	           boundary_residual_m=residual,
	           p_fit_residual=spec["p_residual"],
	           current_fit_residual=spec["current_residual"],
	           Psi_total=spec["Psi_total"],
	           map_source=spec["map_source"], map_refine=spec["map_refine"],
	           profiles=args.profiles,
	           ftol=args.ftol, xtol=args.xtol, gtol=args.gtol,
	           nested=bool(eq.is_nested()),
	           self_consistent=args.self_consistent,
	           self_consistent_mix=args.self_consistent_mix,
	           self_consistent_tol=args.self_consistent_tol,
	           self_consistent_history=sc_history,
	           self_consistent_reached=(bool(sc_history)
	                                    and sc_history[-1]["field_change"]
	                                        < args.self_consistent_tol),
	           g_boundary=args.g_boundary, p_boundary=args.p_boundary)

	out = args.out or os.path.join(HERE, "runs")
	os.makedirs(out, exist_ok=True)
	if not args.no_field:
		field, _ = flux_on_grid(eq, R, Z, n_rho=args.cloud,
		                        n_theta=2 * args.cloud)
		path = os.path.join(out, f"{args.stem}-desc-M{args.M}.npz")
		write_npz(path, R, Z, field, spec, eq, result, row,
		          dict(desc_M=args.M, desc_nit=row["nit"],
		               axis_r=R_axis, axis_z=0.0))
		row["npz"] = path

	print(f"  {args.stem}  M = {args.M}  dim_x = {row['dim_x']}  "
	      f"map from {spec['map_source']}")
	print(f"    nit {row['nit']:3d}  cost {row['cost']:.3e}  nested {row['nested']}"
	      f"  success {row['success']}")
	print(f"    cold {t_cold:8.2f} s ( build {t_build:.2f} + solve {t_solve:.2f} )"
	      f"   warm {(row['warm_seconds'] if warm else float('nan')):8.2f} s"
	      f"   convert {t_convert:5.2f} s")
	print(f"    psi_ax {psi_ax:.9e}   reference {spec['psi_ax_reference']:.9e}"
	      f"   rel {abs(psi_ax-spec['psi_ax_reference'])/spec['psi_ax_reference']:.3e}")
	print(f"    R_axis {R_axis:.6f}   boundary residual {residual:.3e} m")
	if args.json:
		json.dump(row, open(args.json, "w"), indent=2)
	return 0


if __name__ == "__main__":
	sys.exit(main())
