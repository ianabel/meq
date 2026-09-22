#!/usr/bin/env python
"""ONE **FREE-BOUNDARY** DESC SOLVE OF A MEQ FREE-BOUNDARY CASE, TIMED.

`descrun.py` solves inside a prescribed `Gamma`, which is what the six
`examples/fixed-*.toml` machines are.  This solves for the boundary as well,
from a coil set, which is what MEQ's `examples/*-filament.toml` machines are --
and it is a different objective, a different optimiser and a different set of
things the two codes are and are not given.

WHY ONLY A LIMITED CASE REACHES THIS FILE.  DESC's `BoundaryError` enforces
`B_out . n = 0` and pressure balance ON THE LCFS, and the LCFS is
`equilibrium.surface`, a `FourierRZToroidalSurface` -- a truncated Fourier
series in the poloidal angle.  Every free-boundary machine in `examples/` but
one is DIVERTED, so its LCFS is a separatrix through an X-point and has a
CORNER, which no truncation turns; `CLAUDE.md` records the same fact from the
other side, as the reason the fixed ladder is posed at `psi_n = 0.95` rather
than at the separatrix.  A LIMITED plasma's edge is the smooth flux surface
tangent to the limiter and is representable, so
`examples/limited-tokamak-filament.toml` is the case this exists for.

WHAT EACH CODE IS GIVEN, WHICH IS NOT THE SAME THING AND CANNOT BE MADE SO.

    MEQ          the coil currents, `p'( Psi )` and `g g'( Psi )`, a target
                 `I_p`, and A LIMITER CONTACT.  The plasma's size is pinned by
                 the contact: the edge is the flux surface through that point.
    DESC         the coil currents, `p( rho )` and `I( rho )`, and the TOTAL
                 TOROIDAL FLUX `Psi`.  There is no wall and no limiter in the
                 formulation at all; the plasma's size is pinned by `Psi`.

Both are complete statements of one equilibrium and neither is the other's.
The conversion carries `Psi` over from the reference's own flux map, so what
DESC is told about the size is what MEQ's limiter contact implies -- but it is
told it as a different scalar, and a disagreement in the answer can live there
as well as in either code's discretisation.  That is an asymmetry of the
COMPARISON and is listed with the others in README.md.

THE POSING IS UNTIMED AND THE SOLVE IS TIMED, exactly as `descrun.py` does it
and for the same reason: `Phi( Psi )` needs the shapes of the interior flux
surfaces, which is an equilibrium, so DESC cannot be handed this problem
without one.  `race_desc.py --self-consistent` is the machinery that removes
the reference from that loop and it applies here unchanged.

WHERE THE BOUNDARY STARTS IS A CHOICE AND IT IS THE ONE THAT DECIDES WHAT THE
CLOCK MEANS.  `--boundary reference` starts DESC at the reference's own LCFS,
which is a warm start; `--boundary circle` starts it at a circle through the
limiter contact, which knows only the machine.  MEQ's own arm on this case is
warm-started -- its cold guess does not converge, and
`examples/limited-tokamak-filament.toml` says so in its own header -- so
`reference` is the matched default and `circle` is the honest cold number
beside it.

USAGE

    ./venv/bin/python descfreeb.py limited-tokamak-filament -M 12
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
import descrun                                                # noqa: E402

MU0 = 4.0e-7 * np.pi


def coil_field(meta):
	"""The machine, as DESC sees it: the PF filaments AND the toroidal field.

	**THE TOROIDAL FIELD IS NOT OPTIONAL HERE AND LEAVING IT OUT IS NOT A SMALL
	ERROR.**  `BoundaryError` enforces two residuals on the LCFS, `B_out . n =
	0` and `B_out^2 - B_in^2 - 2 mu0 p = 0`, and `B_in` carries the toroidal
	field `g/R` while a set of poloidal field coils produces none at all.  A
	tokamak's toroidal field comes from a TF coil that no Grad-Shafranov input
	names, because GS sees `g` only through `g dg/dpsi` and MEQ's answer does
	not depend on the additive constant -- so `[[coils]]` is the whole machine
	for MEQ and is not the whole machine for DESC.

	MEASURED, with the PF coils alone: the solve starts with a normal-field
	error of 1.7e-03 normalised, which is small because the reference's LCFS
	really is a free-boundary surface, and a magnetic-pressure error of
	**9.0e-01 normalised**, which is `B_phi^2` entire.  The optimiser then
	trades the good residual against the impossible one and walks the boundary
	**0.81 m** -- on a plasma of minor radius 0.34 -- for a `psi_ax` 92% wrong.
	A free-boundary run that is missing the toroidal field does not fail, it
	converges to a different machine.

	`g_at_gamma` is the scalar that supplies it, and it is the same one the
	conversion already uses to fix the total toroidal flux -- README.md's
	second row of *the three things the conversion needs*.  Here it is used a
	second time and for a different reason.

	ONE CIRCULAR FILAMENT PER CONDUCTOR, for the poloidal half.

	MEQ's `[conductors] Model = "filament"` collapses each conductor to a point
	filament at its centre -- `Config.cpp` refuses a `QuadratureOrder` beside it
	saying so -- and an axisymmetric point filament IS a circular current loop.
	`FourierPlanarCoil` with a single `r_n` is that loop exactly, so this half
	of the posing costs nothing: checked against the closed form for a loop's
	field on its own axis, `mu0 I R^2 / 2( R^2 + z^2 )^{3/2}`, it agrees to
	nine figures.

	THAT IS NOT TRUE OF THE MESHED SIBLING and the difference is why this races
	the filament file.  `examples/limited-tokamak.toml` gives the same four
	conductors a 0.1 x 0.1 m extent carrying a uniform current density, which
	neither freegs4e's reference nor a single `FourierPlanarCoil` is; racing it
	would put a ( w/d )^2 modelling difference of about 1.6e-02 on the near
	field into a comparison resolving 1e-04.
	"""
	from desc.coils import CoilSet, FourierPlanarCoil
	from desc.magnetic_fields import SumMagneticField, ToroidalMagneticField

	if meta.get("conductor_model") != "filament":
		raise SystemExit(
			"%s has [conductors] Model = %r.  A DESC coil is a curve, so only "
			"a filament model converts without a quadrature this script does "
			"not do -- race the -filament sibling."
			% (meta["stem"], meta.get("conductor_model")))
	coils = [FourierPlanarCoil(current=float(c["current"]),
	                           center=[0.0, 0.0, float(c["Z"])],
	                           normal=[0.0, 0.0, 1.0],
	                           r_n=float(c["R"]), name=c["name"])
	         for c in meta["coils"]]
	# B0*R0/R with B0*R0 = g( Gamma ), i.e. B_phi = g/R exactly.  R0 = 1 makes
	# B0 the number `g_at_gamma` holds, in T m, rather than a field strength at
	# some other radius that would have to be converted.
	toroidal = ToroidalMagneticField(B0=float(meta["g_at_gamma"]), R0=1.0)
	return SumMagneticField(CoilSet(*coils), toroidal), toroidal


def initial_surface(spec, meta, M, kind, blend=None):
	"""The LCFS the free-boundary solve starts from.

	`reference` is the MXH curve the conversion already fitted, i.e. the
	reference's own answer.  `circle` is a circle of the machine's design minor
	radius about its design centre, which knows the limiter and nothing else --
	the analogue of MEQ's cold guess, and the thing to use if the clock is
	meant to be a cold one.

	`blend` IS THE INSTRUMENT FOR THE QUESTION M-162 LEAVES OPEN, and it is the
	only reason this takes a third form.  That sweep found `psi_ax` non-monotone
	in `M` and unmoved by a four-order tightening of the optimiser tolerance,
	which says the runs are converged and converged SOMEWHERE ELSE -- multiple
	minima rather than resolution.  If that is right then the SEED decides which
	one is reached, and `M` is only a proxy for it, because the starting surface
	is re-fitted at every `M` and a different fit lands in a different basin.

	So: hold `M` fixed and move the seed continuously.  `blend = t` samples
	`( 1 - t )` of the circle plus `t` of the reference's curve at the same
	poloidal parameter and fits that, so `t = 0` and `t = 1` reproduce the two
	forms above and everything between is a seed no equilibrium corresponds to.
	A monotone answer in `t` says the objective has one minimum and the seed
	only decides how far the optimiser has to walk; a scattered one says basins.
	"""
	from desc.geometry import FourierRZToroidalSurface

	if kind == "reference" and blend is None:
		theta, Rb, Zb = convert.desc_boundary_samples(meta, max(512, 8 * M))
		return FourierRZToroidalSurface.from_values(
			np.column_stack([Rb, np.zeros_like(Rb), Zb]), theta,
			M=M, N=0, NFP=1, sym=spec["up_down_symmetric"]), (theta, Rb, Zb)

	R0, a = float(meta["R0"]), float(meta["minor_radius"])
	if blend is not None:
		# THE SAME PARAMETER ON BOTH CURVES, which needs no search and is the
		# choice mxh.pointwise_error makes for the same reason.  Both are
		# sampled at DESC's own theta, so the circle is written with the sign
		# convention `desc_boundary_samples` already uses rather than with a
		# `+sin` that would blend a curve against its own reflection.
		theta, Rb, Zb = convert.desc_boundary_samples(meta, max(512, 8 * M))
		Rc = R0 + a * np.cos(theta)
		Zc = -a * np.sin(theta)
		t = float(blend)
		Rb = (1.0 - t) * Rc + t * Rb
		Zb = (1.0 - t) * Zc + t * Zb
		return FourierRZToroidalSurface.from_values(
			np.column_stack([Rb, np.zeros_like(Rb), Zb]), theta,
			M=M, N=0, NFP=1, sym=spec["up_down_symmetric"]), (theta, Rb, Zb)

	surface = FourierRZToroidalSurface(
		R_lmn=[R0, a], modes_R=[[0, 0], [1, 0]],
		Z_lmn=[-a], modes_Z=[[-1, 0]], NFP=1, sym=spec["up_down_symmetric"])
	theta = np.linspace(0.0, 2.0 * np.pi, 512, endpoint=False)
	# -sin, MATCHING THE SURFACE JUST BUILT.  `modes_Z = [[-1, 0]]` with
	# coefficient -a is Z = -a sin( theta ), and this tuple is what
	# `boundary_residual` measures the solved boundary against -- so a +sin
	# here compares the answer with the seed's own REFLECTION and reports a
	# displacement wrong by up to 2a on an up-down symmetric machine, where
	# the two curves are the same SET of points traversed oppositely and
	# nothing looks amiss.  The blend branch above already had this right.
	return surface, (theta, R0 + a * np.cos(theta), -a * np.sin(theta))


def solve_free(eq, field, M_boundary, ftol=1e-6, xtol=1e-8, gtol=1e-10,
               maxiter=50, verbose=0, sheet=False):
	"""The free-boundary step: DESC's own tutorial recipe, on a tokamak.

	`field_fixed=True` because the coils are given and only the plasma moves;
	the alternative is single-stage coil optimisation, which is a different
	question.  No sheet current unless asked: `BoundaryError`'s third residual
	is only present when the surface is a `FourierCurrentPotentialField`, and
	its own docstring says it is negligible when the edge pressure is zero --
	which it is here, `p( Gamma ) = 0` exactly.
	"""
	from desc.objectives import (BoundaryError, FixBoundaryR, FixBoundaryZ,
	                             FixCurrent, FixPressure, FixPsi, ForceBalance,
	                             ObjectiveFunction)

	if sheet:
		from desc.magnetic_fields import FourierCurrentPotentialField
		eq.surface = FourierCurrentPotentialField.from_surface(eq.surface,
		                                                       M_Phi=4)
	objective = ObjectiveFunction(
		BoundaryError(eq=eq, field=field, field_fixed=True))
	constraints = (ForceBalance(eq=eq), FixPressure(eq=eq), FixCurrent(eq=eq),
	               FixPsi(eq=eq))
	# THE HIGH POLOIDAL MODES ARE HELD, which the tutorial does too.  The modes
	# the objective can move are the ones the boundary genuinely has; leaving
	# every mode free lets a nearly circular surface grow structure the data
	# does not support, at the cost of conditioning.
	modes_R = eq.surface.R_basis.modes[
		np.max(np.abs(eq.surface.R_basis.modes), 1) > M_boundary, :]
	modes_Z = eq.surface.Z_basis.modes[
		np.max(np.abs(eq.surface.Z_basis.modes), 1) > M_boundary, :]
	constraints = constraints + (FixBoundaryR(eq=eq, modes=modes_R),
	                             FixBoundaryZ(eq=eq, modes=modes_Z))
	started = time.perf_counter()
	# `ftol=None` HANDS DESC ITS OWN DEFAULT, which is not the same as passing
	# a tight number: a trust-region method whose tolerances are set below what
	# its model can resolve stops on "failure to predict improvement" rather
	# than on any of them, and then reports a point that is not a minimum.
	# Distinguishing that from the problem being hard needs the library's own
	# settings as a control.
	kwargs = {k: v for k, v in
	          dict(ftol=ftol, xtol=xtol, gtol=gtol).items() if v is not None}
	eq, result = eq.optimize(objective, constraints,
	                         optimizer="proximal-lsq-exact",
	                         maxiter=maxiter, verbose=verbose, copy=False,
	                         **kwargs)
	return eq, result, time.perf_counter() - started


def main():
	ap = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument("stem", nargs="?", default="limited-tokamak-filament")
	ap.add_argument("-M", type=int, default=12)
	ap.add_argument("--terms", type=int, default=12)
	ap.add_argument("--refine", type=int, default=4)
	ap.add_argument("--boundary", choices=("reference", "circle"),
	                default="reference")
	ap.add_argument("--blend", type=float, default=None,
	                help="seed the boundary at ( 1 - t ) circle + t reference, "
	                     "at fixed M; see initial_surface() for the question "
	                     "this answers")
	ap.add_argument("--free-modes", type=int, default=2,
	                help="poloidal modes the boundary may move in; higher "
	                     "ones are held")
	ap.add_argument("--maxiter", type=int, default=50)
	ap.add_argument("--ftol", type=float, default=1e-6)
	ap.add_argument("--xtol", type=float, default=1e-8)
	ap.add_argument("--gtol", type=float, default=1e-10)
	ap.add_argument("--desc-defaults", action="store_true",
	                help="pass none of the three tolerances, so DESC uses its "
	                     "own -- the control for a trust-region collapse")
	ap.add_argument("--current-sign", type=float, default=+1.0)
	ap.add_argument("--warm", type=int, default=0,
	                help="solve the free-boundary step this many more times "
	                     "in the same process, to separate JAX compilation "
	                     "from arithmetic")
	ap.add_argument("--out", default="")
	ap.add_argument("--verbose", type=int, default=3)
	args = ap.parse_args()

	case = convert.load_case(args.stem)
	meta = case["meta"]
	if not meta.get("free_boundary"):
		raise SystemExit("%s-meta.json does not say free_boundary; descrun.py "
		                 "is the entry point for a fixed case" % args.stem)

	field = convert.field_from_reference(case["ref"], float(meta["level"]),
	                                     args.refine)
	spec = convert.build(case, n_fit=args.terms, field=field)
	coils, toroidal = coil_field(meta)

	print("  %s   free boundary, M = %d, boundary from %s"
	      % (args.stem, args.M,
	         args.boundary if args.blend is None
	         else "blend t = %.3f" % args.blend))
	print("    %d filament coils, %+.6e A in total, plus a toroidal field "
	      "B_phi = %.6f / R T"
	      % (len(meta["coils"]), sum(c["current"] for c in meta["coils"]),
	         float(meta["g_at_gamma"])))
	print("    total toroidal flux Psi   %.9e Wb" % spec["Psi_total"])
	print("    enclosed current at rho=1 %.6e A against the reference's "
	      "%.6e" % (spec["current"][-1], float(meta["Ip"])))

	started = time.perf_counter()
	surface, truth = initial_surface(spec, meta, args.M, args.boundary,
	                                 blend=args.blend)
	eq, _ = descrun.build_equilibrium(spec, args.M,
	                                  current_sign=args.current_sign)
	eq.surface = surface
	# A CONSISTENT INTERIOR FIRST.  `BoundaryError` moves the boundary of an
	# equilibrium; starting it from an unsolved one asks the optimiser to do
	# both jobs at once, and the tutorial does not.
	eq, fixed_result = eq.solve(objective="force", optimizer="lsq-exact",
	                            ftol=1e-8, xtol=1e-8, gtol=1e-10, maxiter=100,
	                            verbose=args.verbose, copy=False)
	fixed_seconds = time.perf_counter() - started

	eq, result, free_seconds = solve_free(eq, coils, args.free_modes,
	                                      ftol=None if args.desc_defaults else args.ftol,
	                                      xtol=None if args.desc_defaults else args.xtol,
	                                      gtol=None if args.desc_defaults else args.gtol,
	                                      maxiter=args.maxiter,
	                                      verbose=args.verbose)
	cold = fixed_seconds + free_seconds

	warm = []
	for _ in range(args.warm):
		again = eq.copy()
		_, _, seconds = solve_free(again, coils, args.free_modes,
		                           ftol=args.ftol, maxiter=args.maxiter,
		                           verbose=0)
		warm.append(seconds)

	rho, psi_line = descrun.poloidal_flux(eq)
	psi_ax_desc = float(psi_line[0])
	# WHY THE OPTIMISER STOPPED, AND IT IS NOT THE ITERATION COUNT.  A sweep
	# of --ftol that leaves both the answer and the count unmoved does NOT
	# establish that the run was converged rather than halted: if some OTHER
	# termination test fired -- xtol, gtol, a step-size floor, a cap -- then
	# ftol was never what stopped it and tightening it could not have changed
	# anything.  The exit reason names which, and nothing else does.
	# A FAILED OPTIMISATION IS NOT AN ANSWER, AND IT PRINTS ONE.  Every
	# quantity below is computed from `eq` whether or not the optimiser
	# reached a minimum, so without this the run reports a plausible psi_ax,
	# a plausible boundary and a wall clock, with nothing saying that the
	# trust region collapsed four iterations in.  MEASUREMENTS.md M-162 is
	# what that cost: a whole campaign of resolution and seed sweeps built on
	# runs that had all failed, and three published revisions chasing
	# structure in the scatter.
	succeeded = bool(result.get("success", False))
	print("\n    optimiser exit        %r"
	      % str(result.get("message", "?")).replace("\n", " "))
	if not succeeded:
		print("    *** THE OPTIMISER DID NOT SUCCEED.  Nothing below is an "
		      "answer: it is the iterate the solve stopped at.  The usual "
		      "cause here is starting AT a stationary point -- a reference "
		      "seed -- where the trust region cannot predict improvement.")
	print("    fixed-boundary seed   %8.2f s" % fixed_seconds)
	print("    free-boundary step    %8.2f s   ( %d iterations )"
	      % (free_seconds, result.get("nit", -1)))
	print("    cold total            %8.2f s" % cold)
	if warm:
		print("    warm free-boundary    %8.2f s   ( best of %d )"
		      % (min(warm), len(warm)))
	print("    psi_ax  DESC %.9e   reference %.9e   rel %.3e"
	      % (psi_ax_desc, float(meta["psi_axis"]) - float(meta["psi_surface"]),
	         abs(psi_ax_desc - (float(meta["psi_axis"])
	                            - float(meta["psi_surface"])))
	         / abs(float(meta["psi_axis"]) - float(meta["psi_surface"]))))

	# HOW FAR THE BOUNDARY MOVED, which is the one diagnostic a fixed-boundary
	# run has no analogue of and the thing this solve is for.
	moved = descrun.boundary_residual(eq.surface, *truth)
	print("    boundary moved from its start by %.4e m" % moved)
	theta, Rb, Zb = convert.desc_boundary_samples(meta, 512)
	against_reference = descrun.boundary_residual(eq.surface, theta, Rb, Zb)
	print("    boundary against the reference's own LCFS   %.4e m"
	      % against_reference)

	R, Z = np.array(case["ref"]["R"]), np.array(case["ref"]["Z"])
	psi_meq, _ = descrun.flux_on_grid(eq, R, Z)
	out = args.out or os.path.join(
		HERE, "runs", "%s-descfreeb-M%d-%s.npz"
		% (args.stem, args.M,
		   args.boundary if args.blend is None
		   else "blend%.3f" % args.blend))
	os.makedirs(os.path.dirname(out), exist_ok=True)
	descrun.write_npz(out, R, Z, psi_meq, spec, eq, result,
	                  dict(cold=cold, free=free_seconds, fixed=fixed_seconds,
	                       warm=warm),
	                  dict(free_boundary=True,
	                       # THE SEED THAT ACTUALLY RAN, not the flag.  --blend
	                       # OVERRIDES --boundary, so recording args.boundary
	                       # here labelled every blended run with whatever
	                       # --boundary happened to default to: `blend = 0.000`
	                       # is the CIRCLE and was stored as "reference".  The
	                       # filename was right throughout and this field was
	                       # not, which is the same defect as the one the
	                       # optimiser_success guard below exists to close --
	                       # a harness reporting something other than what ran.
	                       boundary_start=(
	                           args.boundary if args.blend is None
	                           else "blend t = %.3f" % args.blend),
	                       optimiser_success=succeeded,
	                       optimiser_message=str(result.get("message", "")),
	                       boundary_moved=moved,
	                       boundary_vs_reference=against_reference,
	                       seconds_cold=cold, seconds_free=free_seconds,
	                       seconds_fixed=fixed_seconds,
	                       seconds_warm=np.array(warm, dtype=float),
	                       # WHAT THE RUN WAS GIVEN.  Without these an archived
	                       # file from the tolerance study cannot be told from
	                       # one from the seed study, and neither can be
	                       # compared with a run taken later -- which is what
	                       # made M-162's first sweep unrescuable by re-reading
	                       # it.
	                       ftol=(np.nan if args.desc_defaults else args.ftol),
	                       xtol=(np.nan if args.desc_defaults else args.xtol),
	                       gtol=(np.nan if args.desc_defaults else args.gtol),
	                       maxiter=args.maxiter,
	                       desc_defaults=bool(args.desc_defaults),
	                       M=args.M, free_modes=args.free_modes))
	print("    wrote %s" % out)


if __name__ == "__main__":
	main()
