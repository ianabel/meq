"""A freegs4e machine's conductors, flattened into the rectangles MEQ takes.

WHY THIS EXISTS.  `Machine.coils` is a list of whatever the machine definition
wrote down, and that is four different classes with four different geometries:

    Coil( R, Z, turns )            a FILAMENT -- a point source, no extent
    ShapedCoil( polygon, turns )   a polygon, meshed by freegs4e itself
    Circuit( [ ( name, coil, multiplier ) ... ] )   N conductors on one supply
    Solenoid( Rs, Zmin, Zmax, Ns ) Ns filaments in a vertical line at Rs

Only the first two have `.R`, and fgsref.py's saver read `float( c.R )` off
every entry -- so four of the seven diverted references could not be SAVED,
though all seven solve.  Nothing downstream could ask a machine where its
conductors are.

MEQ's `[[coils]]` is one thing: a RECTANGLE carrying a uniform current density,
which is what meq::Coil integrates and what tools/mesh/halfdisc.py fragments
into the mesh.  So the question this module answers is "what rectangles is this
machine", and it answers it once for every consumer -- the reference saver, the
TOML writer and the Green's-function guess.

THE CURRENT IS THE TOTAL THROUGH THE CROSS-SECTION, WHICH IS NOT `c.current`.
freegs4e's Greens functions carry the turns:

    Coil.controlPsi     = Greens( R, Z, ... ) * turns
    Circuit.psi         sets each sub-coil's current to multiplier * I_circuit
    Solenoid.controlPsi = sum over Ns filaments, each at the SAME current

so a conductor's total is `turns * multiplier * I`, and a solenoid's is
`Ns * I`.  MEQ's `[[coils]] Current` is documented as the total through the
cross-section, so these line up -- and getting it wrong is a silent factor of
up to 324 (MAST-U's solenoid), which would still converge to something.

TWO MODELLING DIFFERENCES THIS CANNOT REMOVE, AND BOTH ARE RECORDED RATHER THAN
HIDDEN.

  * A FILAMENT HAS NO SIZE.  freegs4e's `Coil` is a point source with a log
    singularity in psi; MEQ's is a rectangle of uniform current density, which
    is finite everywhere.  They agree away from the conductor and cannot agree
    at it -- MEASUREMENTS.md M-87 measures that as the whole of a 2.97e-01
    relative L-infinity, all of it inside two filament coils, against 5.3e-04
    over the 5478 nodes outside them.  `half` below is the size given to a
    filament and 0.05 m is the default, because that is the half-extent of the
    ShapedCoils in the same machine.
  * A SOLENOID IS A STACK OF FILAMENTS in freegs4e and one tall rectangle here.
    At Ns = 100 to 324 turns the stack is very nearly the uniform winding, so
    this is the direction where MEQ's model is the BETTER one; the difference
    is again confined to the conductor.

CONDUCTORS MUST NOT OVERLAP, and on two of these machines they would.  TCV's T2
and T3 sit 0.037 m apart at the same height, so the default 0.05 m half-width
makes two rectangles that intersect -- and a region that is inside two coils is
not something an element attribute can express, nor something `occ.fragment`
will produce.  MAST-U is worse and is what makes this a real algorithm rather
than a nearest-neighbour cap: its solenoid is 0.04 m wide and 3.16 m TALL, so it
intersects the Px coils 1.23 m away in Z, which no test on centre distances
sees.

`shrink_to_fit` therefore tests rectangles against rectangles, and only shrinks
extents this module INVENTED -- a filament's size, a solenoid's width.  A
ShapedCoil's polygon and a solenoid's winding length are what the machine says,
and editing those would make the comparison about a machine nobody described.
It says what it moved.
"""

import numpy as np


# The default half-extent given to a conductor with no geometry of its own,
# metres.  It is the ShapedCoils' own half-extent in TestTokamak, which is the
# machine both codes have been compared on most.
DEFAULT_HALF = 0.05

# A solenoid is thin: freegs4e gives it no radial extent at all, so anything
# here is invention.  Kept small, and kept off the axis.
SOLENOID_HALF_WIDTH = 0.02

# The fraction of the gap to a neighbour a conductor may occupy.  0.45 leaves a
# tenth of the gap between two rectangles that have both been shrunk, which is
# enough for the mesher to put an element in.
PACKING = 0.45


class Conductor(dict):
	"""One rectangle: label, centre, half-extents, total current."""

	@property
	def Rmin(self):
		return self["R"] - self["half_width"]

	@property
	def Rmax(self):
		return self["R"] + self["half_width"]

	@property
	def zmin(self):
		return self["Z"] - self["half_height"]

	@property
	def zmax(self):
		return self["Z"] + self["half_height"]


def _rect(label, R, z, half_width, half_height, current, kind):
	return Conductor(label=str(label), R=float(R), Z=float(z),
	                 half_width=float(half_width),
	                 half_height=float(half_height),
	                 current=float(current), kind=kind)


def _shaped_extent(coil):
	"""The bounding half-extents of a ShapedCoil's polygon, and its centre."""
	shape = np.asarray(coil.shape, dtype=float)
	r0, r1 = shape[:, 0].min(), shape[:, 0].max()
	z0, z1 = shape[:, 1].min(), shape[:, 1].max()
	return 0.5*(r0 + r1), 0.5*(z0 + z1), 0.5*(r1 - r0), 0.5*(z1 - z0)


def flatten(machine, half=DEFAULT_HALF):
	"""Every conductor of a freegs4e Machine, as rectangles with total currents.

	The order is the machine's own, depth first through circuits, so a caller
	writing [[coils]] blocks in this order has them lined up with the `10 + i`
	element attributes tools/mesh/halfdisc.py assigns in command-line order.
	"""
	out = []

	def walk(label, item, scale):
		cls = type(item).__name__

		if cls == "Circuit":
			# The sub-coil's current is multiplier * the circuit's, and each
			# sub-coil then has turns of its own.
			for name, coil, multiplier in item.coils:
				outer = 1.0 if scale is None else scale
				walk("%s.%s" % (label, name), coil,
				     outer*float(multiplier)*float(item.current))
			return

		if cls == "Solenoid":
			# ONE RECTANGLE FOR THE WHOLE WINDING, carrying Ns turns' worth.
			# freegs4e stacks Ns filaments between Zsmin and Zsmax; at the turn
			# counts these machines use that is a uniform winding to any
			# accuracy that matters outside it.
			rs = float(item.Rs)
			width = min(SOLENOID_HALF_WIDTH, 0.25*rs)
			out.append(_rect(label, rs, 0.5*(item.Zsmin + item.Zsmax),
			                 width, 0.5*abs(item.Zsmax - item.Zsmin),
			                 (1.0 if scale is None else scale)
			                 *float(item.Ns)*float(item.current),
			                 "solenoid"))
			return

		# `scale` carries the circuit's contribution, and None -- not zero --
		# is what says there is no circuit. A circuit sitting at exactly zero
		# current is a real state of a converged machine, and testing the
		# number's truthiness would silently replace it with the sub-coil's own
		# stale current.
		turns = float(getattr(item, "turns", 1.0))
		current = float(item.current) if scale is None else scale
		if cls == "ShapedCoil":
			R, z, hw, hh = _shaped_extent(item)
			out.append(_rect(label, R, z, hw, hh, turns*current, "shaped"))
			return

		out.append(_rect(label, item.R, item.Z, half, half, turns*current,
		                 "filament"))

	for label, item in machine.coils:
		walk(label, item, None)
	return out


def _overlaps(a, b):
	return (a.Rmin < b.Rmax and b.Rmin < a.Rmax
	        and a.zmin < b.zmax and b.zmin < a.zmax)


# Which half-extents are REAL geometry and which are this module's invention.
# Only an invention may be shrunk: a ShapedCoil's polygon and a solenoid's
# winding length are what the machine says, and quietly editing them would make
# the comparison about a machine nobody described.
_INVENTED = {
	"filament": ("half_width", "half_height"),
	"solenoid": ("half_width",),
	"shaped": (),
}


def shrink_to_fit(conductors, packing=PACKING, axis_margin=0.9, rounds=12):
	"""Cap INVENTED half-extents so no two rectangles intersect and none
	reaches R = 0.

	Returns ( label, hw_before, hh_before, hw, hh ) for whatever moved, so a
	caller can report it rather than a geometry silently becoming a different
	machine.

	A CENTRE-DISTANCE TEST IS NOT AN OVERLAP TEST, and MAST-U is where that
	shows. Its solenoid is 0.04 m wide and 3.16 m tall and its Px coils sit at
	Z = +-1.2285, so the two are 1.23 m apart in Z and their rectangles
	INTERSECT anyway -- the solenoid spans both. Separating them has to happen
	in R, and it has to leave the solenoid's height alone.

	THE PASS IS COLLECTIVE RATHER THAN PAIR-BY-PAIR, and on an up-down
	symmetric machine that is visible: resolving ( solenoid, PxU ) and then
	( solenoid, PxL ) shrinks the solenoid in between, so the two Px coils come
	out at DIFFERENT widths -- 0.0294 and 0.0333 -- on a machine whose own
	description is symmetric. Every pair is measured against the same geometry
	and the caps are applied together.
	"""
	before = [(c["half_width"], c["half_height"]) for c in conductors]

	# THE AXIS FIRST. meq::Coil, [[coils]] and halfdisc.py all refuse a
	# conductor reaching R = 0, because the operator's 1/R is not integrable
	# through it.
	for c in conductors:
		if "half_width" in _INVENTED.get(c["kind"], ()):
			c["half_width"] = min(c["half_width"], axis_margin*c["R"])

	for _ in range(rounds):
		caps = {}

		def cap(index, key, value):
			was = caps.get((index, key))
			caps[(index, key)] = value if was is None else min(was, value)

		clashing = False
		for i in range(len(conductors)):
			for jj in range(i + 1, len(conductors)):
				a, b = conductors[i], conductors[jj]
				if not _overlaps(a, b):
					continue
				clashing = True
				key, share_a, share_b = _separation(a, b, packing)
				if share_a is not None:
					cap(i, key, share_a)
				if share_b is not None:
					cap(jj, key, share_b)
		if not clashing:
			break
		for (index, key), value in caps.items():
			conductors[index][key] = min(conductors[index][key], value)

	moved = []
	for c, (hw0, hh0) in zip(conductors, before):
		if (c["half_width"] < hw0 - 1e-15) or (c["half_height"] < hh0 - 1e-15):
			moved.append((c["label"], hw0, hh0, c["half_width"],
			              c["half_height"]))
	return moved


def _separation(a, b, packing):
	"""Where two overlapping rectangles should be pulled apart, and to what.

	Returns ( key, cap_a, cap_b ), either cap being None for an extent that is
	real geometry and must not move.  The axis is the one that can be separated
	with the least loss, and the room is shared by WATER-FILLING rather than
	proportionally: everybody gets an equal share, anybody already smaller than
	their share keeps what they have, and the slack goes to the rest.  That is
	what stops a thin conductor being made thinner still -- MAST-U's solenoid
	is 0.02 m wide against Px's 0.05, and a proportional split would take the
	solenoid to 0.008 and leave a 3.16 m tall sliver for the mesher to fill.
	"""
	options = []
	for key, gap in (("half_width", abs(a["R"] - b["R"])),
	                 ("half_height", abs(a["Z"] - b["Z"]))):
		free_a = key in _INVENTED.get(a["kind"], ())
		free_b = key in _INVENTED.get(b["kind"], ())
		if not (free_a or free_b) or gap <= 0.0:
			continue
		fixed = (0.0 if free_a else a[key]) + (0.0 if free_b else b[key])
		room = 2.0*packing*gap - fixed
		if room <= 0.0:
			continue

		free = [(x[key], which) for x, which in ((a, "a"), (b, "b"))
		        if (free_a if which == "a" else free_b)]
		shares = _water_fill([v for v, _ in free], room)
		if sum(shares) <= 0.0:
			continue
		got = dict(zip([w for _, w in free], shares))
		loss = 1.0 - sum(shares)/max(sum(v for v, _ in free), 1e-300)
		options.append((loss, key, got.get("a"), got.get("b")))

	if not options:
		raise SystemExit(
			"conductors %s and %s overlap and neither has an extent this "
			"module is allowed to shrink -- both are real geometry, so it is "
			"the machine's own description that has them intersecting"
			% (a["label"], b["label"]))

	_, key, cap_a, cap_b = min(options)
	return key, cap_a, cap_b


def _water_fill(wants, room):
	"""Share `room` among `wants`: an equal share each, anybody wanting less
	than their share takes what they want, and the slack is shared again."""
	out = [0.0]*len(wants)
	live = list(range(len(wants)))
	while live:
		share = room/len(live)
		small = [i for i in live if wants[i] <= share]
		if not small:
			for i in live:
				out[i] = share
			break
		for i in small:
			out[i] = wants[i]
			room -= wants[i]
		live = [i for i in live if i not in small]
	return out


def _overlaps(a, b):
	return (a.Rmin < b.Rmax and b.Rmin < a.Rmax
	        and a.zmin < b.zmax and b.zmin < a.zmax)


# Which half-extents are REAL geometry and which are this module's invention.
# Only an invention may be shrunk: a ShapedCoil's polygon and a solenoid's
# winding length are what the machine says, and quietly editing them would make
# the comparison about a machine nobody described.
_INVENTED = {
	"filament": ("half_width", "half_height"),
	"solenoid": ("half_width",),
	"shaped": (),
}


def to_arrays(conductors):
	"""The flat table fgsref.py saves into the .npz, as parallel arrays."""
	return dict(
		coil_labels=np.array([c["label"] for c in conductors]),
		coil_R=np.array([c["R"] for c in conductors], dtype=float),
		coil_Z=np.array([c["Z"] for c in conductors], dtype=float),
		coil_half_width=np.array([c["half_width"] for c in conductors],
		                         dtype=float),
		coil_half_height=np.array([c["half_height"] for c in conductors],
		                          dtype=float),
		coil_currents=np.array([c["current"] for c in conductors], dtype=float),
		coil_kind=np.array([c["kind"] for c in conductors]),
	)


def from_npz(d, half=DEFAULT_HALF):
	"""Read the table back. Older references carry no half-extents; they get
	the default, which is what the hand-written examples/diverted-tokamak.toml
	used and so keeps that file reproducible from its own reference."""
	n = len(d["coil_labels"])
	half_width = (d["coil_half_width"] if "coil_half_width" in d.files
	              else np.full(n, half))
	half_height = (d["coil_half_height"] if "coil_half_height" in d.files
	               else np.full(n, half))
	kind = (d["coil_kind"] if "coil_kind" in d.files
	        else np.array(["filament"]*n))
	return [_rect(d["coil_labels"][i], d["coil_R"][i], d["coil_Z"][i],
	              half_width[i], half_height[i], d["coil_currents"][i],
	              str(kind[i]))
	        for i in range(n)]


def bounding_radius(conductors):
	"""The smallest circle about the origin containing every conductor."""
	return max(np.hypot(max(abs(c.Rmin), abs(c.Rmax)),
	                    max(abs(c.zmin), abs(c.zmax)))
	           for c in conductors)
