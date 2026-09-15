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
  * `current` AND `turns` SEPARATELY.  ShapedCoil's total is `current*turns`,
    the same convention Coil uses, so a substitution that carries both across
    is exact -- and a Solenoid's Ns turns become the ShapedCoil's turns, which
    is what makes its total `Ns*I` as freegs4e's stack of Ns filaments is.

WHAT CHANGES, DELIBERATELY: a filament acquires the extent it never had, and a
solenoid becomes one uniform winding instead of Ns discrete filaments.  Both are
the point.  A ShapedCoil already in the machine is LEFT ALONE -- its polygon is
what the machine says, and rewriting it as a bounding rectangle would make the
comparison about a machine nobody described.
"""

import conductors


def _square(rect):
	"""The rectangle as a polygon, anticlockwise from its lower inboard corner."""
	r0, r1 = rect.rmin, rect.rmax
	z0, z1 = rect.zmin, rect.zmax
	return [(r0, z0), (r1, z0), (r1, z1), (r0, z1)]


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

	def substitute(label, item):
		cls = type(item).__name__

		if cls == "Circuit":
			inner = [(name, substitute("%s.%s" % (label, name), coil), mult)
			         for name, coil, mult in item.coils]
			return Circuit(inner, current=float(item.current),
			               control=bool(item.control))

		# A polygon the machine itself described: left exactly as it is.
		if cls == "ShapedCoil":
			return item

		rect = byLabel[label]
		if cls == "Solenoid":
			# Ns turns of a uniform winding rather than Ns stacked filaments.
			return ShapedCoil(_square(rect), current=float(item.current),
			                  turns=int(item.Ns), control=bool(item.control),
			                  npoints=npoints)

		# A filament, which is every remaining case.
		return ShapedCoil(_square(rect), current=float(item.current),
		                  turns=int(getattr(item, "turns", 1)),
		                  control=bool(getattr(item, "control", True)),
		                  npoints=npoints)

	coils = [(label, substitute(label, item)) for label, item in mach.coils]
	return Machine(coils, getattr(mach, "wall", None)), moved
