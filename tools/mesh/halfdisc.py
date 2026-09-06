#!/usr/bin/env python3
"""A half-disc reaching the axis, with rectangular coils MESHED TO, for FB-6.

WHY A MESHING LIBRARY AT ALL.  Everything MEQ has meshed so far came from
`mfem::Mesh::MakeCartesian2D` -- a rectangle of triangles -- with curved and
non-rectangular geometry reached by the extension technique on top of it.  Free
boundary needs something that grid cannot be: a SEMICIRCLE centred on the axis,
reaching r = 0 EXACTLY, with the coil rectangles as their own subdomains.

    r = 0 IS NOT A CONVENIENCE, IT IS WHAT MAKES THE DECOUPLING LEGAL.
    FREE-BOUNDARY-PLAN.md section 3 expands the exterior field in Gegenbauer
    functions of order -1/2 on a semicircle CENTRED ON THE AXIS; the expansion
    is a statement about that geometry and about no other.  A domain stopping
    at r = 0.05 is not a slightly worse version of it -- the basis does not
    span the exterior of that shape, and meq::ExteriorDtN's diagonal symbol is
    not the DtN map of it.  So `--rho` sets Gamma's radius and the flat side is
    the axis itself, and `--check` asserts on `r == 0.0` EXACTLY rather than on
    a tolerance, because a mesh that merely comes close has a different
    exterior problem and nothing downstream would say so.

WHY GMSH, AND WHY PYTHON.  MFEM reads gmsh's format natively --
`Mesh::ReadGmshMesh`, no converter -- so the coupling between MEQ and the
mesher is a FILE and nothing else.  gmsh ships a C++ API and a Python module,
and the Python module is ctypes over the same `libgmsh.so`, generated from the
same `api/gen.py`: identical function set, identical meshes, the same code doing
the work.  So linking it into `meq_core` would buy no capability and would cost
`libgmsh-dev` as a build dependency, on a tree whose CI cannot build the solver
at all.  This is a tool that writes a file, which is the whole interface.

Nothing here is MEQ-specific beyond the defaults; it is a plain gmsh script.

    python3 halfdisc.py --rho 1.5 --size 0.15 \
        --coil 0.90 -0.86 0.25 0.12 --coil 0.90 0.74 0.25 0.12 \
        --plasma 0.6 0.0 0.45 0.65 --plasma-size 0.04 -o machine.msh


THE COILS ARE MESHED TO RATHER THAN CUT
---------------------------------------

FB-2's finding applied rather than a convenience.  A conductor's geometry is
PRESCRIBED INPUT -- it does not move with the solution -- so it can always be
aligned to, and aligning it took the measured rates from 1.33 / 1.27 / 1.09 to
1.99 / 2.88 / 3.01.  Only the plasma support genuinely moves, and FB-4 measures
what that costs.  `occ.fragment` is what does it: it splits the disc along the
coil rectangles so their edges are mesh edges.


THE BOUNDARY GROUPS ARE THE OUTER BOUNDARY AND NOTHING ELSE, AND THE FIRST
VERSION OF THIS FILE GOT THAT WRONG IN A WAY NOTHING WOULD HAVE REPORTED
-------------------------------------------------------------------------

Gamma and the axis were selected by asking each 1-D entity where its centre of
mass was -- on the axis, or not.  After `occ.fragment` the model also contains
each CONDUCTOR'S OWN OUTLINE as 1-D entities, and every one of them is "not on
the axis", so they joined Gamma.  A 1-D physical group is written to the `.msh`,
and MFEM's reader turns 1-D elements into BOUNDARY elements whether or not they
are topologically on the boundary.  Measured on the example above: 64 boundary
faces, of which 32 on the arc, 20 on the axis and **12 on the coil rectangles**,
carrying Gamma's attribute in the middle of the mesh.

That is not a cosmetic mis-tag.  Attribute 1 is what carries the transferred
exterior datum, so the run would have imposed Gamma's Dirichlet condition on
twelve edges buried inside the conductors -- an interior constraint, converging
at full order to the wrong equilibrium, with nothing in any output saying so.

`gmsh.model.getBoundary( surfaces, combined = True )` is the fix and is the
right question in the first place: combined over ALL the surfaces, every edge
shared by two of them cancels, and what survives is the outer boundary.  The
axis is then told from the arc by geometry, which is what `--check` re-asserts
from the written file.


WHAT IT PRODUCES, MEASURED THROUGH MFEM RATHER THAN THROUGH GMSH
----------------------------------------------------------------

`--check` re-reads the `.msh`, but MFEM is the reader that matters, so the
numbers below were taken by loading the file with `mfem::Mesh` and walking it.
The example at the top of this file, `--rho 1.5 --size 0.15` with two coils:

    vertices 264, elements 474, boundary elements 52
    r over the mesh                          [0, 1.5]
    vertices at EXACTLY r == 0.0             21
    element attr  1  r [0.0000, 1.5000]  z [-1.5000, 1.5000]
    element attr 10  r [0.9000, 1.1500]  z [-0.8600, -0.7400]
    element attr 11  r [0.9000, 1.1500]  z [ 0.7400,  0.8600]
    bdr attr 1 (Gamma)  32 faces, worst | |x| - rho | = 4.441e-16
    bdr attr 2 (axis)   20 faces, worst |r|           = 0.000e+00
    faces on neither the arc nor the axis    0     <- was 12

`r == 0.0` and `|r| = 0.000e+00` are EXACT ZEROS and not small numbers.  The
axis is a straight geometric edge at x = 0, so every node OCC places on it lands
there bit exactly; that is the property the exterior expansion needs and the one
`--check` asserts without a tolerance.

AT `--order 2` THE ARC IS THE ARC AND THE AXIS IS STILL EXACT.  The mid-edge
nodes are placed on the OCC curve rather than on the chord: 48 arc nodes at
4.441e-16 from the true circle, 30 axis nodes at exactly 0.0.  That is why the
geometry is built as a disc-minus-halfplane rather than as an arc plus two
lines -- OCC owns the curvature, and gmsh has something to place them on.


THE GRADING, AND WHAT IT IS WORTH
---------------------------------

A coarse disc aggressively refined where the plasma is.  Gamma has to sit far
enough out that the exterior expansion converges, and the plasma occupies a
small part of what that encloses, so meshing the whole disc at the plasma's
resolution spends most of the elements on vacuum carrying no source.  Measured,
same two coils, plasma box r [0.30, 1.00] z [-0.55, 0.55] at h = 0.035:

    rho    graded    uniform at h = 0.035    saving
    1.5     2,965            7,672           2.6x
    3.0     3,368           30,392           9.0x

-- and the saving GROWS with rho, which is the direction the exterior coupling
pushes.  The graded mesh costs almost nothing extra to move Gamma out from 1.5
to 3.0 (2,965 to 3,368 elements) where the uniform one quadruples.
"""
import argparse
import math
import sys

import gmsh


AXIS_TOLERANCE = 1.0e-9

VACUUM_ATTRIBUTE = 1
FIRST_COIL_ATTRIBUTE = 10
GAMMA_ATTRIBUTE = 1
AXIS_ATTRIBUTE = 2


def _coil_of(centre, coils):
    """The index of the coil whose rectangle contains `centre`, or None."""
    x, y, _ = centre
    for i, (r, z, w, h) in enumerate(coils):
        if r <= x <= r + w and z <= y <= z + h:
            return i
    return None


def build(rho, coils, size, out, order=1, plasma=None, plasma_size=None,
          coil_size=None, transition=None):
    """A half-disc of radius `rho` about the origin, r >= 0, with `coils` a list
    of (rmin, zmin, width, height) rectangles fragmented into it.

    `plasma` is an optional (rmin, zmin, width, height) box to refine inside, to
    `plasma_size`; `coil_size` refines inside the conductors.  Both default to
    the background `size`, so a caller who asks for neither gets a uniform mesh
    and the size field is not installed at all.

    Returns a dict of counts, which `main` prints and `--check` re-derives from
    the written file rather than trusting.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("halfdisc")
        occ = gmsh.model.occ

        disc = occ.addDisk(0, 0, 0, rho, rho)
        # A disc plus a cut rather than an arc plus two lines: OCC then owns the
        # curvature, so `--order > 1` puts the high-order nodes ON the true arc
        # instead of on the chord.
        box = occ.addRectangle(-rho, -rho, 0, rho, 2 * rho)
        half, _ = occ.cut([(2, disc)], [(2, box)])

        tags = [occ.addRectangle(r, z, 0, w, h) for (r, z, w, h) in coils]
        if tags:
            occ.fragment(half, [(2, t) for t in tags])
        occ.synchronize()

        surfaces = [(2, t) for _, t in gmsh.model.getEntities(2)]

        # ---- element attributes: one per coil, one for everything else ----
        #
        # Ordered by the COIL's index rather than by the surface tag, so that
        # attribute 10 + i is the i-th --coil on the command line whatever OCC
        # decided to call it.  A caller writing [[coils]] blocks in the same
        # order then has the two lists lined up, which is the only way the
        # attributes are usable at all.
        owner = {}
        for dim, tag in surfaces:
            owner[tag] = _coil_of(occ.getCenterOfMass(dim, tag), coils)

        vacuum = sorted(t for t in owner if owner[t] is None)
        gmsh.model.addPhysicalGroup(2, vacuum, VACUUM_ATTRIBUTE)
        gmsh.model.setPhysicalName(2, VACUUM_ATTRIBUTE, "vacuum")

        coil_surfaces = {}
        for tag, which in owner.items():
            if which is not None:
                coil_surfaces.setdefault(which, []).append(tag)
        for which in sorted(coil_surfaces):
            attribute = FIRST_COIL_ATTRIBUTE + which
            gmsh.model.addPhysicalGroup(2, sorted(coil_surfaces[which]), attribute)
            gmsh.model.setPhysicalName(2, attribute, "coil%d" % which)

        # ---- boundary attributes: the OUTER boundary, and only that ----
        #
        # combined = True is the whole point: every edge shared by two surfaces
        # cancels, so the conductors' own outlines -- which fragment created and
        # which are interior -- do not appear.  See the module docstring for
        # what happened when they did.
        outer = gmsh.model.getBoundary(surfaces, combined=True, oriented=False)

        # Gamma and the axis must not share an attribute: Gamma carries the
        # transferred exterior datum and the axis carries a plain essential
        # condition.  They are told apart by geometry -- the axis is the only
        # boundary at r = 0 -- rather than by tag order, which OCC does not
        # promise and which fragment reshuffles.
        arc, axis = [], []
        for dim, tag in outer:
            x, _, _ = occ.getCenterOfMass(1, abs(tag))
            (axis if abs(x) < AXIS_TOLERANCE else arc).append(abs(tag))
        gmsh.model.addPhysicalGroup(1, sorted(arc), GAMMA_ATTRIBUTE)
        gmsh.model.setPhysicalName(1, GAMMA_ATTRIBUTE, "Gamma")
        gmsh.model.addPhysicalGroup(1, sorted(axis), AXIS_ATTRIBUTE)
        gmsh.model.setPhysicalName(1, AXIS_ATTRIBUTE, "axis")

        # ---- the size field ----
        _install_size_field(rho, coils, size, plasma, plasma_size, coil_size,
                            transition)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.model.mesh.generate(2)
        if order > 1:
            gmsh.model.mesh.setOrder(order)
        gmsh.write(out)

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        elements = sum(len(block) // gmsh.model.mesh.getElementProperties(kind)[3]
                       for kind, block
                       in zip(*[gmsh.model.mesh.getElements(2)[i] for i in (0, 2)]))
        return dict(surfaces=len(surfaces), coils=len(coil_surfaces),
                    arc_curves=len(arc), axis_curves=len(axis),
                    nodes=len(node_tags), elements=elements)
    finally:
        gmsh.finalize()


def _install_size_field(rho, coils, size, plasma, plasma_size, coil_size,
                        transition):
    """A background size, refined inside the plasma box and inside the coils.

    THIS IS THE COARSE-DISC-AND-REFINED-PLASMA SHAPE, and it is what makes the
    half-disc affordable.  Gamma has to sit far enough out that the exterior
    expansion converges, and the plasma occupies a small part of what that
    encloses; meshing the whole disc at the plasma's resolution would spend most
    of the elements on vacuum that carries no source and needs no accuracy.

    Boxes rather than a distance field, because both regions ARE boxes and a
    Box field is exact on them -- `Thickness` gives the graded transition
    outward, so there is no jump in element size at the boundary of the refined
    region for the mesher to struggle with.
    """
    fields = []

    def box_field(rmin, zmin, width, height, inside):
        tag = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(tag, "XMin", rmin)
        gmsh.model.mesh.field.setNumber(tag, "XMax", rmin + width)
        gmsh.model.mesh.field.setNumber(tag, "YMin", zmin)
        gmsh.model.mesh.field.setNumber(tag, "YMax", zmin + height)
        gmsh.model.mesh.field.setNumber(tag, "VIn", inside)
        gmsh.model.mesh.field.setNumber(tag, "VOut", size)
        gmsh.model.mesh.field.setNumber(tag, "Thickness", transition)
        fields.append(tag)

    if plasma is not None and plasma_size is not None:
        box_field(*plasma, inside=plasma_size)
    if coil_size is not None:
        for (r, z, w, h) in coils:
            box_field(r, z, w, h, inside=coil_size)

    if not fields:
        # No refinement asked for: a plain uniform mesh, and the field machinery
        # is left out entirely rather than installed with every value equal.
        gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", size / 8.0)
        return

    smallest = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(smallest, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(smallest)

    # WITHOUT THESE THREE THE FIELD IS ONLY A CEILING.  gmsh otherwise also
    # takes sizes from the geometry -- from point sizes, from the curvature of
    # the arc, and by extending the boundary spacing inward -- and takes the
    # MINIMUM of all of them, so a coarse background over a curved boundary
    # quietly refines anyway and `--size` stops meaning what it says.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def check(path, rho, coils):
    """Re-read the written file and assert the properties MEQ depends on.

    THE FILE, NOT THE MODEL IN MEMORY.  A check run against the live gmsh model
    would share every assumption with the code that built it -- CLAUDE.md
    records three separate occasions in this tree where a fixture was misread
    the same way as the thing it checked.  This opens the `.msh` again from
    disk, which is the artefact MFEM will read.

    Returns a list of failure strings; empty means every property holds.
    """
    problems = []
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(path)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        # BY TAG, NOT BY POSITION.  getNodes() returns `coords` ordered to match
        # `tags`, and the tags are neither sorted nor necessarily contiguous --
        # they come out grouped by the entity that owns them.  Indexing `coords`
        # with `tag - 1` therefore reads a DIFFERENT node, and the first version
        # of this check did exactly that: it reported the axis group at
        # r = 1.4928 and a coil spanning most of the disc, which looked like a
        # catastrophic meshing failure and was the checker being wrong.
        at = {int(t): (coords[3 * i], coords[3 * i + 1])
              for i, t in enumerate(tags)}
        r = coords[0::3]
        z = coords[1::3]

        # r = 0 EXACTLY, on real nodes rather than within a tolerance.  See the
        # module docstring: a domain that merely comes close to the axis has a
        # different exterior problem.
        on_axis = sum(1 for value in r if value == 0.0)
        if on_axis == 0:
            problems.append("no node is at exactly r = 0.0; the exterior "
                            "expansion is not valid on this domain")
        if min(r) < 0.0:
            problems.append("a node is at r = %.17g < 0" % min(r))
        if max(r) > rho * (1.0 + 1.0e-12):
            problems.append("a node is at r = %.17g, outside rho = %g"
                            % (max(r), rho))

        # Every boundary element belongs to the arc or to the axis and to
        # nothing else.  This is the check the first version of this file did
        # not have and needed.
        stray = {}
        counts = {}
        for group_dim, group_tag in gmsh.model.getPhysicalGroups(1):
            for entity in gmsh.model.getEntitiesForPhysicalGroup(1, group_tag):
                types, _, node_tags = gmsh.model.mesh.getElements(1, entity)
                for kind, block in zip(types, node_tags):
                    per = gmsh.model.mesh.getElementProperties(kind)[3]
                    counts[group_tag] = counts.get(group_tag, 0) + len(block) // per
                    for node in block:
                        x, y = at[int(node)]
                        near_axis = abs(x) <= AXIS_TOLERANCE
                        near_arc = abs(math.hypot(x, y) - rho) <= 1.0e-9 * rho
                        if group_tag == AXIS_ATTRIBUTE and not near_axis:
                            stray.setdefault("axis", []).append((x, y))
                        if group_tag == GAMMA_ATTRIBUTE and not near_arc:
                            stray.setdefault("Gamma", []).append((x, y))
        for where, points in stray.items():
            problems.append(
                "%d boundary nodes carry the %s attribute without lying on it, "
                "the first at r = %.4f z = %.4f -- an interior edge has been "
                "tagged as a boundary"
                % (len(points), where, points[0][0], points[0][1]))
        for wanted, name in ((GAMMA_ATTRIBUTE, "Gamma"), (AXIS_ATTRIBUTE, "axis")):
            if counts.get(wanted, 0) == 0:
                problems.append("boundary attribute %d (%s) has no elements"
                                % (wanted, name))

        # Each coil attribute covers its own rectangle and no more.
        for i, (rmin, zmin, width, height) in enumerate(coils):
            attribute = FIRST_COIL_ATTRIBUTE + i
            entities = gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
            if len(entities) == 0:
                problems.append("coil %d has no element attribute %d"
                                % (i, attribute))
                continue
            lo = [float("inf")] * 2
            hi = [float("-inf")] * 2
            for entity in entities:
                _, _, node_tags = gmsh.model.mesh.getElements(2, entity)
                for block in node_tags:
                    for node in block:
                        for axis_index, value in enumerate(at[int(node)]):
                            lo[axis_index] = min(lo[axis_index], value)
                            hi[axis_index] = max(hi[axis_index], value)
            want = (rmin, rmin + width, zmin, zmin + height)
            got = (lo[0], hi[0], lo[1], hi[1])
            if max(abs(a - b) for a, b in zip(want, got)) > 1.0e-9:
                problems.append(
                    "coil %d (attribute %d) covers r [%.6f, %.6f] z [%.6f, %.6f], "
                    "not the requested r [%.6f, %.6f] z [%.6f, %.6f]"
                    % ((i, attribute) + got[:2] + got[2:] + want[:2] + want[2:]))
    finally:
        gmsh.finalize()
    return problems


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rho", type=float, default=1.5,
                   help="radius of the semicircle Gamma, metres. Gamma must "
                        "enclose every conductor: the exterior expansion is "
                        "only valid where the field is Delta*-harmonic")
    p.add_argument("--size", type=float, default=0.15,
                   help="background element size, metres")
    p.add_argument("--order", type=int, default=1,
                   help="geometric order; >1 curves the arc")
    p.add_argument("--coil", type=float, nargs=4, action="append",
                   metavar=("RMIN", "ZMIN", "WIDTH", "HEIGHT"), default=None,
                   help="a rectangular coil, repeatable. Meshed TO, and given "
                        "element attribute 10 + i in this order")
    p.add_argument("--coil-size", type=float, default=None,
                   help="element size inside the conductors; defaults to --size")
    p.add_argument("--plasma", type=float, nargs=4, default=None,
                   metavar=("RMIN", "ZMIN", "WIDTH", "HEIGHT"),
                   help="a box to refine inside, where the plasma is expected")
    p.add_argument("--plasma-size", type=float, default=None,
                   help="element size inside --plasma; required with it")
    p.add_argument("--transition", type=float, default=None,
                   help="width of the graded transition out of a refined "
                        "region, metres; defaults to four background sizes")
    p.add_argument("--check", action="store_true",
                   help="re-read the written file and assert MEQ's "
                        "preconditions on it; non-zero exit if any fails")
    p.add_argument("-o", "--out", default="halfdisc.msh")
    a = p.parse_args()

    coils = a.coil or []
    for (r, z, w, h) in coils:
        if r <= 0.0:
            p.error("a coil at RMIN = %g reaches or crosses the axis; the "
                    "operator's 1/r is not integrable through r = 0, which is "
                    "the same refusal meq::Coil and [[coils]] make" % r)
        if w <= 0.0 or h <= 0.0:
            p.error("a coil of WIDTH %g HEIGHT %g has no cross-section" % (w, h))
        if math.hypot(max(abs(r), abs(r + w)), max(abs(z), abs(z + h))) >= a.rho:
            p.error("a coil reaches r = %g z = %g, outside or on Gamma at "
                    "rho = %g. The exterior expansion assumes the field is "
                    "Delta*-harmonic outside Gamma, and a current there is not"
                    % (r + w, max(abs(z), abs(z + h)), a.rho))
    if (a.plasma is None) != (a.plasma_size is None):
        p.error("--plasma and --plasma-size go together: a region with no size "
                "refines nothing, and a size with no region has nowhere to act")

    transition = a.transition if a.transition is not None else 4.0 * a.size
    rep = build(a.rho, coils, a.size, a.out, a.order, a.plasma, a.plasma_size,
                a.coil_size, transition)
    print("  %s: %d elements, %d nodes, %d surfaces (%d coil), "
          "%d arc curves, %d axis curves"
          % (a.out, rep["elements"], rep["nodes"], rep["surfaces"],
             rep["coils"], rep["arc_curves"], rep["axis_curves"]))

    if a.check:
        problems = check(a.out, a.rho, coils)
        for problem in problems:
            print("  FAIL: %s" % problem, file=sys.stderr)
        if problems:
            return 1
        print("  check: r reaches 0 exactly, Gamma and the axis are the outer "
              "boundary and nothing else, every coil covers its rectangle")
    return 0


if __name__ == "__main__":
    sys.exit(main())
