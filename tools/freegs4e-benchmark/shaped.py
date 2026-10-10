"""Every conductor of a freegs4e machine as a ShapedCoil on MEQ's own rectangle.

WHY THIS EXISTS.  MEQ's `meq::Coil` is a rectangle carrying a uniform current
density.  Six of the seven diverted references model most of their conductors as
FILAMENTS -- `freegs4e.machine.Coil` is a point source with a log singularity in
psi -- and one as a stack of them.  MEASUREMENTS.md M-87 measures what that
costs on the diverted DIII-D comparison: the whole of a 2.97e-01 relative
L-infinity, all of it inside the conductors, against 5.3e-04 over the 5478 nodes
outside them.  That is a MODELLING difference between the two codes and not a
discretisation error in either, so it cannot be refined away and it cannot be
fixed on MEQ's side.

`freegs4e.shaped_coil.ShapedCoil` is freegs4e's OWN answer: a polygon,
triangulated, with a Gauss rule per triangle whose weights sum to one, so
`controlPsi` is the AVERAGE Greens over the cross-section and the total current
is spread uniformly across it.  That is MEQ's model term for term, in the
reference's own code rather than as a MEQ-side correction.  `fgsref.py`'s
`make_limited_shaped` already does this for the LIMITED case by hand; this
module does it for an arbitrary machine, which is what the diverted seven need
because four of them are built from Circuits and Solenoids rather than a flat
list.

THE RECTANGLES ARE conductors.py's, NOT NEW ONES, and that is the whole point.
A filament has no size, so the size it is GIVEN has to come from somewhere, and
it has to be the same size MEQ meshes or the two codes are still solving
different machines.  `conductors.flatten` is what MEQ's `[[coils]]` blocks are
written from, `shrink_to_fit` is what keeps TCV's T2/T3 pair and MAST-U's
solenoid from overlapping, and this module rebuilds the machine from THAT list.
Pass the same `half` both places or the comparison means nothing.

WHAT IS PRESERVED, AND WHY EACH MATTERS.

  * THE CIRCUIT TOPOLOGY.  A Circuit is N conductors on one supply, and the
    references are produced by an INVERSE solve -- freegs4e adjusts coil
    currents to hit the xpoint and isoflux constraints in `fgsref.py`'s CASES.
    Flattening a circuit into independent coils would hand that solve more
    degrees of freedom than the machine has and converge to a different
    equilibrium, which would look like a conductor-model result and would not
    be one.  So circuits are rebuilt as circuits, with their multipliers.
  * THE `control` FLAG, per conductor, for the same reason: it is what says
    whether the inverse solve may move that current at all.
  * THE TOTAL CURRENT, with `turns` FOLDED IN.  ShapedCoil's docstring says
    its total is `current*turns`, but its field ignores `turns` -- see
    `substitute` below -- so every ShapedCoil here has turns = 1 and the turns
    are carried by the current (a lone coil) or the circuit multiplier (a
    member). A Solenoid's Ns turns are folded the same way, which keeps its
    total `Ns*I` as freegs4e's stack of Ns filaments has it.

WHAT CHANGES, DELIBERATELY: a filament acquires the extent it never had, and a
solenoid becomes one uniform winding instead of Ns discrete filaments.  Both are
the point.  A ShapedCoil already in the machine is LEFT ALONE -- its polygon is
what the machine says, and rewriting it as a bounding rectangle would make the
comparison about a machine nobody described.
"""

import conductors


def _square(rect):
	"""The rectangle as a polygon, anticlockwise from its lower inboard corner."""
	r0, r1 = rect.Rmin, rect.Rmax
	z0, z1 = rect.zmin, rect.zmax
	return [(r0, z0), (r1, z0), (r1, z1), (r0, z1)]


#: Gauss-Legendre points across a rectangle's SHORTER side; the longer side
#: gets proportionally more. 6 puts every conductor of the seven machines
#: within 1e-07 of a 48 x 48 rule at the plasma.
RULE_POINTS = 6


def rectangle_rule(rect, n=RULE_POINTS):
	"""A tensor Gauss-Legendre rule over the rectangle, as ShapedCoil's
	( R, Z, weight ) list, the weights summing to one.

	WHY NOT ShapedCoil's OWN RULE. It triangulates the polygon -- two
	triangles for a rectangle -- and puts a 6-point rule on each: 12 points
	whatever the conductor's shape. That is accurate to 3e-05 on a square
	(fgsref.py measured it on the limited case) and NOT on a tall conductor:
	measured against a 48 x 48 rule at plasma points, a conductor's own flux is
	off by up to 2.3e-01 on MAST, 1.8e-01 on TCV and 3.0e-01 on MAST-U, whose
	solenoid is 79 times taller than it is wide -- refining the polygon does
	not help, the fan triangulation keeps every triangle full-height. So the
	reference was not solving a uniform rectangle. This rule is, to 1e-07:
	MEASUREMENTS.md M-255. Assigned to `_points`, which is all controlPsi,
	controlBr and controlBz read.
	"""
	import math
	import numpy as np
	hw, hh = float(rect["half_width"]), float(rect["half_height"])
	nr, nz = n, n
	if hh > hw:
		nz = max(n, int(math.ceil(n*hh/hw)))
	else:
		nr = max(n, int(math.ceil(n*hw/hh)))
	xr, wr = np.polynomial.legendre.leggauss(nr)
	xz, wz = np.polynomial.legendre.leggauss(nz)
	return [(float(rect["R"]) + hw*a, float(rect["Z"]) + hh*b, 0.25*p*q)
	        for a, p in zip(xr, wr) for b, q in zip(xz, wz)]


def shaped_machine(mach, half=conductors.DEFAULT_HALF, shrink=True,
                   npoints=6):
	"""`mach` rebuilt with every filament and solenoid a ShapedCoil.

	@param half     the half-extent given to a conductor with no geometry of its
	                own.  Must match what the MEQ side was written with.
	@param shrink   apply conductors.shrink_to_fit, which is what keeps the
	                invented rectangles from overlapping on TCV and MAST-U.
	@param npoints  ShapedCoil's quadrature points per triangle: 1, 3 or 6.
	@return         ( machine, moved ) -- the new Machine, and
	                shrink_to_fit's own report of what it reduced,
	                ( label, hw_before, hh_before, hw, hh ) per entry, so a
	                caller can print it rather than discover it.
	"""
	from freegs4e.machine import Circuit, Machine
	from freegs4e.shaped_coil import ShapedCoil

	# shrink_to_fit MUTATES the list it is given and RETURNS the report of what
	# it moved, which is the opposite of what its name suggests and is why this
	# keeps `rects` rather than rebinding it.
	rects = conductors.flatten(mach, half=half)
	moved = conductors.shrink_to_fit(rects) if shrink else []
	byLabel = {c["label"]: c for c in rects}

	# SHAPEDCOIL IGNORES ITS OWN `turns` IN THE FIELD, so the turns are folded
	# into the current instead. Its docstring says the block carries
	# current*turns, and controlPsi/Br/Bz sum Greens over the quadrature points
	# with weights summing to ONE and never multiply by turns -- so a ShapedCoil
	# with turns = N produces the field of `current`, not of N*current, where
	# freegs4e's filament Coil does multiply. conductors.flatten exports
	# turns*current, as the docstring promises, so every multi-turn conductor
	# of a shaped reference was written N times stronger than the field
	# freegs4e solved with: MAST-U's coils (2 to 324 turns) read about 22x
	# freegs4e's own coil_psi, and MAST's P1 and TCV's OH1 (100 turns each)
	# put those references on other equilibria. MEASUREMENTS.md M-255.
	#
	# So every ShapedCoil here has turns = 1. A lone coil's control variable
	# becomes its TOTAL current; a circuit member's turns go into its circuit
	# multiplier, so the circuit current keeps the filament machine's meaning
	# and the shaped machine is physically the filament one with extent.
	def substitute(label, item):
		"""( replacement, turns folded into it ) for one item."""
		cls = type(item).__name__

		if cls == "Circuit":
			inner = []
			for name, coil, mult in item.coils:
				new, turns = substitute("%s.%s" % (label, name), coil)
				inner.append((name, new, float(mult)*turns))
			return Circuit(inner, current=float(item.current),
			               control=bool(item.control)), 1

		# A polygon the machine itself described: left exactly as it is.
		if cls == "ShapedCoil":
			return item, 1

		rect = byLabel[label]
		turns = int(item.Ns) if cls == "Solenoid" else \
			int(getattr(item, "turns", 1))
		control = bool(getattr(item, "control", True))
		# A Solenoid becomes one uniform winding rather than Ns stacked
		# filaments; a filament acquires its rectangle.
		coil = ShapedCoil(_square(rect), current=float(item.current),
		                  turns=1, control=control, npoints=npoints)
		coil._points = rectangle_rule(rect)
		return coil, turns

	coils = []
	for label, item in mach.coils:
		new, turns = substitute(label, item)
		if turns != 1:
			# A lone multi-turn coil: its control variable is the total.
			new.current = float(new.current)*turns
		coils.append((label, new))
	return Machine(coils, getattr(mach, "wall", None)), moved
