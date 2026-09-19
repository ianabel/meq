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


`--limiter` MESHES THE LIMITER CIRCLE IN, AND THE SOLVER HALF IS NOT BUILT
------------------------------------------------------------------------

Same argument as the coils, one line further: a limiter is prescribed input too,
so it can be aligned to.  What that is FOR is the limiter as a CURVE --
`psi_bnd = max psi` over it, with the contact found rather than prescribed,
which is what every production free-boundary code does and what would take
`freegs4e`'s own grid artefact out of MEQ's comparison against it
(FREE-BOUNDARY-PLAN.md section 7.20).

`meq::GradShafranovSolver` does not do that yet: `setBoundaryFluxPoint()` takes
a POINT and `LimiterConstraint::ExactPoint` evaluates `psi_h` there.  So the
option is the enabling half and buys nothing on its own today -- it is here
because the mesh is the part that has to exist first, and because an interior
curve is exactly the thing this file already records getting wrong once.

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


`--symmetric` MESHES z >= 0 AND REFLECTS IT, BECAUSE A SYMMETRIC GEOMETRY DOES
NOT GIVE A SYMMETRIC TRIANGULATION
-----------------------------------------------------------------------------

`[solver] UpDownSymmetry` projects the iterate onto the fields even in z, dof
for dof, so the two saddles of a double null are exactly degenerate when the
X-point search looks for them.  Being a DOF-FOR-DOF AVERAGE it needs a mesh
that is its own mirror image, and `meq::GradShafranovSolver::buildMirrorMaps`
refuses one that is not -- by name, rather than projecting onto something that
is not a reflection.

A MESHER HANDED A MIRROR-SYMMETRIC MACHINE DOES NOT RETURN A MIRROR-SYMMETRIC
MESH.  Every triangulation picks diagonals and gmsh picks them freely; measured
on MAST-U, 3624 of 4735 vertices have no partner while all 23 conductors are
paired to the last printed digit.  So the symmetry has to be MADE rather than
hoped for: mesh the half z >= 0, reflect it, and share the nodes on the seam.

WHAT THE REFLECTION ACTS ON IS THE WRITTEN FILE, AND THAT IS NOT SQUEAMISHNESS
ABOUT THE API.  A node's image has to be the EXACT negation of its coordinate,
and the coordinate that matters is the decimal gmsh wrote rather than the
double it held.  Negating the TOKEN -- dropping a leading `-`, or adding one --
is exact whatever precision gmsh chose, where negating the double and
reformatting is exact only if the formatter is.  The half is meshed and written
by gmsh exactly as an ordinary run is, and its header and physical names pass
through untouched.

A REFLECTED COIL TAKES THE ATTRIBUTE OF THE COIL IT IS THE MIRROR OF, WHICH IS
NOT THE ATTRIBUTE OF THE HALF IT WAS COPIED FROM.  `10 + i` is the i-th --coil
in COMMAND-LINE ORDER and a machine's [[coils]] blocks are written to match, so
the image of coil i has to come out as coil i's PARTNER j.  Get that backwards
and every conductor in z < 0 carries its opposite number's current: a mesh of a
machine nobody built, converging at full order, with nothing downstream able to
say so.

AND THE REFUSALS ARE ABOUT THE GEOMETRY AND NOT ABOUT THE SIZE FIELD, which is
the one line here worth drawing carefully.  An unpaired conductor, a limiter
circle off the midplane and a vessel outline that is not its own mirror are
REFUSED -- reflecting any of them makes a mesh of a different machine.  An
asymmetric `--plasma` box is reported instead: it says where elements are
small rather than where anything IS, and the half that is meshed gets exactly
the refinement that was asked for.

THE SEAM IS SHARED AND IS NOT A BOUNDARY.  Nodes at z == 0 exactly have no
image and are used by both halves; a node merely NEAR the seam is refused,
because duplicating it gives a mesh that looks perfect and is two disconnected
halves.  The z = 0 line is also part of the meshed half's outer boundary, so it
is told from the arc by geometry and given no physical group at all -- MFEM's
reader turns 1-D elements into BOUNDARY elements wherever they sit, and this
one sits across the middle of the machine.

AND A REFLECTION REVERSES ORIENTATION, so an image element's nodes are
permuted.  The permutation is read out of gmsh's own reference-element tables
rather than tabulated here: `getElementProperties` gives the local coordinates,
the reference reversal is ( u, v ) -> ( v, u ) on a triangle and u -> -u on an
edge, and matching those against the node list gives the permutation at any
order.

`--symmetry-check FILE` IS THE CHECKER, AND IT IS A MODE RATHER THAN A PRIVATE
HELPER because the question it answers -- is this mesh its own mirror image --
is one worth asking of a mesh SOMEBODY ELSE made, and is the question MEQ asks
in `prepare()` a long way into a run.  It pairs every node, asserts the pairing
is an INVOLUTION, maps the element set onto itself, and reports the attribute
pairing and the element orientations it finds.  The involution is the check
with teeth: every way a pairing can go wrong short of leaving a node unmatched
leaves a map that is not its own inverse, which is to say not a reflection.
"""
import argparse
import math
import sys

import gmsh


AXIS_TOLERANCE = 1.0e-9

# THE TOLERANCE meq::GradShafranovSolver::buildMirrorMaps USES, and it is
# relative for the reason recorded there: a mesh generated symmetric mirrors to
# round-off, but a machine's geometry arrives as decimal text and a conductor
# at -1.0972 need not be the exact negation of one at +1.0972.  MEQ scales by
# max( |r| + |z| ) over the mesh, which on a disc of radius rho is rho times
# sqrt( 2 ), so scaling by rho alone is the TIGHTER of the two -- a mesh this
# accepts is one MEQ accepts, and not the other way about.
SYMMETRY_TOLERANCE = 1.0e-9

VACUUM_ATTRIBUTE = 1
FIRST_COIL_ATTRIBUTE = 10
LIMITER_ATTRIBUTE = 20

# EVERYTHING OUTSIDE THE VESSEL, WHICH IS AN ATTRIBUTE AND NOT A GEOMETRY.
#
# `--vessel` fragments a closed polygon into the half-disc and tags what lies
# OUTSIDE it, so that `[source] ExcludeAttributes = [ 30 ]` can say "the plasma
# cannot be here" once, at mesh time, instead of the solve deciding
# inside/outside per element on every residual evaluation.  A vessel does not
# move; the support does.
#
# IT IS NOT CONFINEMENT AND MUST NOT BE READ AS SUCH.  `psi_bnd` confines the
# plasma -- F is zero wherever Psi <= 0 -- and the connectivity fill separates
# the lobes of { Psi > 0 } that a level set leaves joined.  This is for the one
# case neither can settle: several O-points across a saddle, where connectivity
# cannot say which of them is the plasma.
OUTSIDE_VESSEL_ATTRIBUTE = 30

GAMMA_ATTRIBUTE = 1
AXIS_ATTRIBUTE = 2


def _inside_polygon(x, y, polygon):
    """Ray casting, on a closed polygon given as [ (r, z), ... ].

    A centroid test alone is not enough to name the vessel's own fragments --
    the OUTER region's centroid can land inside a convex shape, which is the
    trap `--limiter` already records and pays for with an AREA check beside it.
    """
    inside = False
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            crossing = x1 + (y - y1)*(x2 - x1)/(y2 - y1)
            if x < crossing:
                inside = not inside
    return inside


def _polygon_area(polygon):
    total = 0.0
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        total += x1*y2 - x2*y1
    return abs(total)/2.0


def _coil_of(centre, coils):
    """The index of the coil whose rectangle contains `centre`, or None.

    A None entry is a conductor that lies entirely in z < 0 under `--symmetric`
    and is therefore not built: its elements arrive as the reflection of its
    partner's.  It keeps its place in the list because the index IS the
    attribute -- `10 + i` is the i-th --coil on the command line, and a list
    that closed up over the missing ones would renumber every conductor after
    the first one below the midplane.
    """
    x, y, _ = centre
    for i, rectangle in enumerate(coils):
        if rectangle is None:
            continue
        r, z, w, h = rectangle
        if r <= x <= r + w and z <= y <= z + h:
            return i
    return None


def _mirror_partners(coils, tolerance):
    """For each coil, the index of the coil it is the mirror of about z = 0.

    Greedy and CONSUMING, so two conductors at the same place pair with each
    other rather than both claiming the first match -- which is what makes the
    result an INVOLUTION, the property meq::GradShafranovSolver's own mirror
    map asserts and the one a partner table can get wrong without leaving
    anything unmatched.

    Returns (partner, worst): `partner` a list with partner[partner[i]] == i,
    and `worst` the largest coordinate discrepancy over every pair, which the
    caller prints rather than hides -- a machine paired at 1e-10 is meshed as
    though it were paired exactly, and the number says by how much.

    Raises ValueError naming the first conductor with no partner.
    """
    def discrepancy(i, j):
        ri, zi, wi, hi = coils[i]
        rj, zj, wj, hj = coils[j]
        return max(abs(ri - rj), abs(wi - wj), abs(hi - hj),
                   abs(zj + (zi + hi)))

    partner = [None]*len(coils)
    worst = 0.0
    for i, (r, z, w, h) in enumerate(coils):
        if partner[i] is not None:
            continue
        # A conductor centred on z = 0 is its OWN mirror, and it is the case a
        # partner search cannot find: nothing else is where it is.  The
        # solenoid of every spherical tokamak is one.
        if abs(z + (z + h)) <= tolerance:
            partner[i] = i
            worst = max(worst, abs(z + (z + h)))
            continue
        found = None
        for j in range(len(coils)):
            if (j != i and partner[j] is None
                    and discrepancy(i, j) <= tolerance):
                found = j
                break
        if found is None:
            raise ValueError(
                "conductor %d, r [%.6f, %.6f] z [%.6f, %.6f], has no mirror "
                "partner about z = 0 within %g m. --symmetric meshes z >= 0 "
                "and reflects it, so a machine that is not itself symmetric "
                "would be meshed as a DIFFERENT machine -- every conductor "
                "below the midplane replaced by the image of one above it -- "
                "and nothing downstream could tell"
                % (i, r, r + w, z, z + h, tolerance))
        partner[i] = found
        partner[found] = i
        worst = max(worst, discrepancy(i, found))
    return partner, worst


def _upper_half_of(rectangle):
    """A coil rectangle clipped to z >= 0, or None if it lies wholly below.

    The comparisons are exact rather than toleranced on purpose: a conductor
    whose bottom edge is a rounding error below zero clips to a rectangle whose
    top is the same rounding error off, which is far inside every tolerance
    downstream, where a toleranced clip would leave a rectangle poking a
    rounding error INTO the half-plane that was cut away and hand OCC a sliver.
    """
    r, z, w, h = rectangle
    if z + h <= 0.0:
        return None
    if z >= 0.0:
        return (r, z, w, h)
    return (r, 0.0, w, z + h)


def _upper_half_polygon(polygon):
    """Sutherland-Hodgman against z >= 0, on a closed polygon."""
    clipped = []
    for i in range(len(polygon)):
        px, pz = polygon[i - 1]
        qx, qz = polygon[i]
        if (pz >= 0.0) != (qz >= 0.0):
            t = pz/(pz - qz)
            clipped.append((px + t*(qx - px), 0.0))
        if qz >= 0.0:
            clipped.append((qx, qz))
    # A vertex ON the seam is emitted once as itself and once as the crossing,
    # so consecutive duplicates are the normal case rather than a degeneracy.
    tidy = [point for i, point in enumerate(clipped)
            if point != clipped[i - 1]]
    return tidy


def _polygon_is_mirrored(polygon, tolerance):
    """Whether the closed polygon is its own image under z -> -z.

    THE VERTEX SET BEING SYMMETRIC IS NOT ENOUGH and that is the whole reason
    this is not three lines: a polygon can have a mirror-symmetric set of
    vertices and join them up asymmetrically, which clips to an upper half
    whose reflection is a different outline.  So the test is on the CURVE:
    reflect it, which reverses its orientation, and ask whether the result is a
    cyclic rotation of the original.
    """
    n = len(polygon)
    mirrored = [(r, -z) for (r, z) in reversed(polygon)]
    for shift in range(n):
        rotated = [mirrored[(i + shift) % n] for i in range(n)]
        if all(abs(a[0] - b[0]) <= tolerance and abs(a[1] - b[1]) <= tolerance
               for a, b in zip(polygon, rotated)):
            return True
    return False


def _reversal(kind):
    """The node permutation that reverses a gmsh element of type `kind`.

    Read out of gmsh's own reference-element table rather than tabulated here,
    so it is right at every geometric order instead of at the two somebody
    thought to write down.  `getElementProperties` gives the local coordinates
    in node order; the reference reversal is (u, v) -> (v, u) on a triangle and
    u -> -u on an edge, and the permutation is the one that carries each node
    onto the node at its reflected local coordinate.

    Returns (permutation, dimension), where the image element's j-th node is
    the image of the original's permutation[j]-th.  See the derivation in
    `_reflect_msh`.
    """
    _, dimension, _, count, local, _ = \
        gmsh.model.mesh.getElementProperties(kind)
    local = [float(value) for value in local]
    nodes = [tuple(local[dimension*i:dimension*(i + 1)]) for i in range(count)]
    if dimension == 1:
        span = min(u for (u,) in nodes) + max(u for (u,) in nodes)
        flipped = [(span - u,) for (u,) in nodes]
    elif dimension == 2:
        flipped = [(v, u) for (u, v) in nodes]
    else:
        raise ValueError(
            "gmsh element type %d has dimension %d, and --symmetric knows how "
            "to reverse an edge and a triangle" % (kind, dimension))

    permutation = []
    for want in flipped:
        matches = [k for k, have in enumerate(nodes)
                   if all(abs(a - b) <= 1.0e-9 for a, b in zip(want, have))]
        if len(matches) != 1:
            raise ValueError(
                "gmsh element type %d has no single node at the reflection of "
                "its node %d, so its reversal is not a permutation of its "
                "nodes" % (kind, len(permutation)))
        permutation.append(matches[0])
    return permutation, dimension


def _negated(token):
    """The exact negation of a decimal as WRITTEN, by its sign alone.

    -0 is not reachable here: only a node off the seam is ever negated, and a
    node on the seam is at exactly 0 and keeps its tag.
    """
    return token[1:] if token.startswith("-") else "-" + token


def _reflect_msh(path, attributes, tolerance):
    """Rewrite a z >= 0 mesh in place as the whole mesh, mirrored about z = 0.

    `attributes` maps (dimension, physical tag) to the physical tag its image
    carries -- which is the identity everywhere except the conductors, where
    the image of coil i is coil i's PARTNER and not coil i.

    THE PERMUTATION, DERIVED RATHER THAN ASSERTED.  An element is the map
    Phi(xi) = sum_k N_k(xi) x_k over its nodes.  Reflecting the nodes alone
    gives M . Phi, whose Jacobian has the opposite sign, so the image element
    is M . Phi . R with R the reference reversal.  Its own nodes are then
    y_j = M( x_p(j) ) where R carries reference node j onto reference node
    p(j), since N_k( R(n_j) ) = delta_{k, p(j)}.  `_reversal` returns p, and p
    is its own inverse because R is.

    Returns a dict of what the written file now holds.
    """
    with open(path) as source:
        lines = source.read().splitlines()

    def block(name):
        return lines.index("$" + name), lines.index("$End" + name)

    node_start, node_stop = block("Nodes")
    element_start, element_stop = block("Elements")
    node_lines = lines[node_start + 2:node_stop]
    element_lines = lines[element_start + 2:element_stop]

    # ---- the nodes, and the seam that is shared rather than copied ----
    image = {}
    extra_nodes = []
    seam = 0
    next_node = max(int(line.split()[0]) for line in node_lines)
    for line in node_lines:
        tag, r, z, third = line.split()
        value = float(z)
        if value == 0.0:
            image[int(tag)] = int(tag)
            seam += 1
            continue
        if abs(value) <= tolerance:
            raise ValueError(
                "a node sits at z = %.17g, within %g of the seam without being "
                "exactly on it. Reflecting it would leave TWO nodes where the "
                "mesh needs one, and a mesh whose halves are joined nowhere "
                "reads as a perfectly good mesh of two disconnected pieces"
                % (value, tolerance))
        next_node += 1
        image[int(tag)] = next_node
        extra_nodes.append("%d %s %s %s" % (next_node, r, _negated(z), third))

    # ---- the elements ----
    #
    # The image's elementary tag is shifted past every tag in the file, so the
    # two halves are distinct geometric entities.  MFEM reads the PHYSICAL tag
    # and skips the rest, so this is for every other reader of the file.
    shift = 1
    for line in element_lines:
        field = line.split()
        for value in field[3:3 + int(field[2])]:
            shift = max(shift, 1 + int(value))
    reversal = {}
    extra_elements = []
    surfaces = 0
    next_element = max(int(line.split()[0]) for line in element_lines)
    for line in element_lines:
        field = line.split()
        kind = int(field[1])
        count = int(field[2])
        tags = field[3:3 + count]
        nodes = [int(value) for value in field[3 + count:]]
        if kind not in reversal:
            reversal[kind] = _reversal(kind)
        permutation, dimension = reversal[kind]
        if dimension == 2:
            surfaces += 1
        mirrored = [image[nodes[k]] for k in permutation]
        if sorted(mirrored) == sorted(nodes):
            raise ValueError(
                "element %s lies entirely on the seam, so it is its own image "
                "and writing the reflection would duplicate it. The z = 0 line "
                "carries no physical group, so nothing on it should have been "
                "written at all" % field[0])
        next_element += 1
        if tags:
            tags = list(tags)
            tags[0] = str(attributes.get((dimension, int(tags[0])),
                                         int(tags[0])))
            if len(tags) > 1:
                tags[1] = str(int(tags[1]) + shift)
        extra_elements.append(
            " ".join([str(next_element), str(kind), str(count)] + tags
                     + [str(value) for value in mirrored]))

    whole = (lines[:node_start + 1]
             + [str(len(node_lines) + len(extra_nodes))]
             + node_lines + extra_nodes
             + lines[node_stop:element_start + 1]
             + [str(len(element_lines) + len(extra_elements))]
             + element_lines + extra_elements
             + lines[element_stop:])
    with open(path, "w") as target:
        target.write("\n".join(whole) + "\n")

    return dict(nodes=len(node_lines) + len(extra_nodes),
                elements=2*surfaces, seam_nodes=seam)


def build(rho, coils, size, out, order=1, plasma=None, plasma_size=None,
          coil_size=None, transition=None, limiter=None, vessel=None,
          symmetric=False, partner=None):
    """A half-disc of radius `rho` about the origin, r >= 0, with `coils` a list
    of (rmin, zmin, width, height) rectangles fragmented into it.

    `plasma` is an optional (rmin, zmin, width, height) box to refine inside, to
    `plasma_size`; `coil_size` refines inside the conductors.  Both default to
    the background `size`, so a caller who asks for neither gets a uniform mesh
    and the size field is not installed at all.

    `limiter` is an optional (R0, Z0, a) circle FRAGMENTED IN, for the same
    reason the conductors are: a limiter is PRESCRIBED INPUT and does not move
    with the solution, so it can be aligned to.  What that buys is that the
    limiter contact -- where psi_bnd is pinned -- lies on mesh entities, so the
    contact can be FOUND on the curve rather than prescribed as a point, and
    the trace space lives exactly there.

    ITS INTERIOR TAKES AN ELEMENT ATTRIBUTE AND NOT A 1-D GROUP, WHICH IS THE
    ONE DECISION IN THIS THAT COULD GO SILENTLY WRONG.  Tagging the circle as a
    1-D physical group would be the obvious thing and is exactly the defect this
    module's docstring records: MFEM's reader turns 1-D elements into BOUNDARY
    elements whether or not they are topologically on the boundary, and
    EnableHybridization registers a flux constraint on every marked boundary
    attribute -- so an interior curve tagged that way would impose a boundary
    condition through the middle of the plasma.  A 2-D attribute is safe, and
    the limiter faces are then recoverable as the interface between attribute
    20 and its neighbours, which is a question about element attributes rather
    than about boundary ones.

    `vessel` is an optional closed polygon [ (r, z), ... ] fragmented in the
    same way, whose OUTSIDE takes attribute 30 -- everything within Gamma that
    is not inside the vessel and is not a conductor.  It exists so that
    `[source] ExcludeAttributes` can name a region that can never be plasma
    ONCE, at mesh time, rather than the solve deciding inside/outside per
    element on every residual evaluation.  A vessel does not move and the
    support does, so the work belongs here.

    `symmetric` meshes the half z >= 0 and reflects it, so that the result is
    its own mirror image dof for dof and `[solver] UpDownSymmetry` will take
    it; `partner` is then the conductor pairing from `_mirror_partners`, which
    says what attribute each reflected conductor carries.  The module docstring
    has the argument.  It changes the GEOMETRY that is built -- a conductor
    below the midplane is not built at all and a straddling one is clipped --
    and leaves the SIZE FIELD at its full extent, so the half that is meshed
    gets exactly the refinement the caller asked for.

    Returns a dict of counts, which `main` prints and `--check` re-derives from
    the written file rather than trusting.
    """
    if symmetric and partner is None:
        # `main` works the pairing out first so that it can REFUSE before
        # anything is built, and passes it in; a caller who uses this as a
        # library still gets the refusal, one step later.
        partner, _ = _mirror_partners(coils,
                                      SYMMETRY_TOLERANCE*max(rho, 1.0))

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
        cutters = [(2, box)]
        if symmetric:
            # And the half below the midplane goes the same way, for the same
            # reason: the cut plane is a straight geometric edge at z = 0, so
            # every node OCC places on the seam lands at EXACTLY 0.0 -- which
            # is what lets the two halves share those nodes rather than
            # duplicate them, and is the same property the axis relies on.
            cutters.append((2, occ.addRectangle(-rho, -rho, 0, 2 * rho, rho)))
        half, _ = occ.cut([(2, disc)], cutters)

        # The conductors as BUILT, which under --symmetric is the list clipped
        # to z >= 0 with the ones wholly below it left as None.  The index is
        # the attribute, so the list keeps its length; see `_coil_of`.
        meshed = ([_upper_half_of(rectangle) for rectangle in coils]
                  if symmetric else coils)
        tags = []
        for rectangle in meshed:
            if rectangle is None:
                continue
            r, z, w, h = rectangle
            tags.append(occ.addRectangle(r, z, 0, w, h))
        if limiter is not None:
            lr, lz, la = limiter
            disk = occ.addDisk(lr, lz, 0, la, la)
            if symmetric:
                # A limiter is only ever built centred on the midplane here --
                # main refuses any other -- so its upper half is a half-disc
                # and OCC still owns the curvature of the part that is kept.
                kept, _ = occ.cut([(2, disk)],
                                  [(2, occ.addRectangle(lr - la, lz - la, 0,
                                                        2 * la, la))])
                tags.extend(tag for (_, tag) in kept)
            else:
                tags.append(disk)
        if vessel is not None:
            # A PLANE SURFACE FROM ITS OWN POINTS, not addPolygon: the OCC
            # kernel has no polygon primitive, and building the loop explicitly
            # is also what puts the vertices on the vessel EXACTLY rather than
            # wherever a helper rounds them to -- the same property --coil
            # depends on and halfdisc's own --check asserts.
            outline = _upper_half_polygon(vessel) if symmetric else vessel
            points = [occ.addPoint(r, z, 0) for (r, z) in outline]
            lines = [occ.addLine(points[i], points[(i + 1) % len(points)])
                     for i in range(len(points))]
            tags.append(occ.addPlaneSurface([occ.addCurveLoop(lines)]))
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
            owner[tag] = _coil_of(occ.getCenterOfMass(dim, tag), meshed)

        # THE LIMITER INTERIOR IS TOLD BY AREA AND NOT BY ITS CENTRE OF MASS,
        # and that is not fastidiousness: the OUTER region's centroid on this
        # geometry sits at about r = 1.11, which is INSIDE a limiter circle of
        # R0 = 1.00, a = 0.35, so a centre-of-mass test misclassifies the
        # vacuum as the limiter and every element in the mesh changes
        # attribute.  The areas differ by a factor of twenty-five.
        inside_limiter = []
        if limiter is not None:
            lr, lz, la = limiter
            for _, tag in surfaces:
                if owner[tag] is not None:
                    continue
                x, y, _ = occ.getCenterOfMass(2, tag)
                if (math.hypot(x - lr, y - lz) < la
                        and occ.getMass(2, tag) < 1.05*math.pi*la*la):
                    inside_limiter.append(tag)

        # THE VESSEL'S OUTSIDE, BY THE SAME TWO-PART TEST THE LIMITER USES AND
        # FOR THE SAME REASON: a centroid alone misclassifies the big outer
        # region, whose centre of mass can sit inside a convex vessel.  Area is
        # what separates them, and the two differ by a large factor here as
        # well.  A conductor keeps its own attribute either way -- a coil is
        # outside the vessel on every machine and saying so twice would lose
        # the coil.
        outside_vessel = []
        if vessel is not None:
            area = _polygon_area(vessel)
            for _, tag in surfaces:
                if owner[tag] is not None or tag in inside_limiter:
                    continue
                x, y, _ = occ.getCenterOfMass(2, tag)
                if not (_inside_polygon(x, y, vessel)
                        and occ.getMass(2, tag) < 1.05*area):
                    outside_vessel.append(tag)

        vacuum = sorted(t for t in owner
                        if owner[t] is None and t not in inside_limiter
                        and t not in outside_vessel)
        gmsh.model.addPhysicalGroup(2, vacuum, VACUUM_ATTRIBUTE)
        gmsh.model.setPhysicalName(2, VACUUM_ATTRIBUTE, "vacuum")

        if outside_vessel:
            gmsh.model.addPhysicalGroup(2, sorted(outside_vessel),
                                        OUTSIDE_VESSEL_ATTRIBUTE)
            gmsh.model.setPhysicalName(2, OUTSIDE_VESSEL_ATTRIBUTE,
                                       "outside_vessel")

        if inside_limiter:
            gmsh.model.addPhysicalGroup(2, sorted(inside_limiter),
                                        LIMITER_ATTRIBUTE)
            gmsh.model.setPhysicalName(2, LIMITER_ATTRIBUTE, "limiter")

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
        #
        # AND UNDER --symmetric THERE IS A THIRD KIND, WHICH MUST GET NO GROUP
        # AT ALL.  The z = 0 line bounds the half that is meshed and is
        # INTERIOR to the mesh that is written, so tagging it would impose
        # Gamma's transferred datum along the midplane of the machine -- the
        # same defect this module's docstring records for the conductor
        # outlines, one geometry further on, and the same silence.  It is told
        # from the arc by geometry, which is also how the axis is.
        arc, axis, seam = [], [], []
        for dim, tag in outer:
            x, y, _ = occ.getCenterOfMass(1, abs(tag))
            if abs(x) < AXIS_TOLERANCE:
                axis.append(abs(tag))
            elif symmetric and abs(y) < AXIS_TOLERANCE:
                seam.append(abs(tag))
            else:
                arc.append(abs(tag))
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
        report = dict(surfaces=len(surfaces), coils=len(coil_surfaces),
                      arc_curves=len(arc), axis_curves=len(axis),
                      limiter_surfaces=len(inside_limiter),
                      outside_vessel_surfaces=len(outside_vessel),
                      nodes=len(node_tags), elements=elements,
                      seam_curves=len(seam), seam_nodes=0)

        if symmetric:
            if not seam:
                raise ValueError(
                    "--symmetric found no boundary curve on z = 0, so the "
                    "half-plane cut did not happen and there is nothing to "
                    "reflect about")
            # The image of conductor i is conductor i's PARTNER, which is what
            # makes the reflected half a mesh of THIS machine rather than of
            # one with its lower conductors relabelled.  Everything else --
            # the vacuum, the limiter, outside the vessel, Gamma, the axis --
            # is its own image and needs no entry.
            attributes = {(2, FIRST_COIL_ATTRIBUTE + i):
                          FIRST_COIL_ATTRIBUTE + partner[i]
                          for i in range(len(coils))}
            report.update(_reflect_msh(out, attributes,
                                       SYMMETRY_TOLERANCE*max(rho, 1.0)))
            # The conductor count is of the FILE and not of the half-model,
            # since a conductor below the midplane exists only as the image of
            # its partner and would otherwise be reported missing.
            report["coils"] = len(set(coil_surfaces)
                                  | {partner[i] for i in coil_surfaces})
        return report
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



def check_symmetry(path, tolerance=None):
    """Is this mesh its own mirror image about z = 0?  The file, re-read.

    THE QUESTION meq::GradShafranovSolver::buildMirrorMaps ASKS, asked here
    where it costs a second rather than there where it costs a run: that one
    runs in `prepare()`, so a mesh it will refuse is refused after the spaces
    are built.  It is also worth asking of a mesh this tool did not make, which
    is why it is a mode of its own and not a private half of `--check`.

    Three findings, in ascending order of teeth:

      * every node has a partner within `tolerance` -- the count of those that
        do not is the headline, and is 3624 of 4735 on a gmsh triangulation of
        a perfectly symmetric machine;
      * the pairing is an INVOLUTION, map[ map[ i ] ] == i.  Every way a
        pairing can go wrong without leaving anything unmatched -- two nodes
        claiming one partner, a tolerance wide enough to reach a neighbour --
        leaves a map that is not its own inverse, which is to say not a
        reflection;
      * the ELEMENT set maps onto itself under that map.  Symmetric nodes do
        not make a symmetric mesh: the diagonal a triangulation picks is
        exactly what a node-only check cannot see.

    And two things a reflection has to have got right, reported alongside: the
    physical tag each attribute's image carries, which for a conductor is its
    PARTNER's and not its own, and the sign of every element, since a
    reflection reverses orientation.

    Returns (problems, report); `problems` empty means every property holds.
    """
    problems = []
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(path)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        # BY TAG, NOT BY POSITION, for the reason `check` records at length.
        at = {int(tag): (coords[3*i], coords[3*i + 1])
              for i, tag in enumerate(tags)}
        if not at:
            return ["%s carries no nodes" % path], dict(nodes=0)
        if tolerance is None:
            extent = max(abs(r) + abs(z) for (r, z) in at.values())
            tolerance = SYMMETRY_TOLERANCE*max(extent, 1.0)

        def cell(value):
            return int(round(value/tolerance))

        table = {}
        for tag, (r, z) in at.items():
            table.setdefault((cell(r), cell(z)), []).append(tag)

        mirror = {}
        unpaired = []
        seam = 0
        worst = 0.0
        for tag, (r, z) in at.items():
            if abs(z) <= tolerance:
                seam += 1
            found, offset = None, 0.0
            # The 3 x 3 sweep is not fastidiousness: a coordinate a hair either
            # side of a cell boundary hashes to a neighbour of the cell its
            # partner is in, and the lookup would miss a pair that is exact.
            for dr in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in table.get((cell(r) + dr, cell(-z) + dz), ()):
                        gap = max(abs(at[other][0] - r), abs(at[other][1] + z))
                        if gap <= tolerance:
                            found, offset = other, gap
            if found is None:
                unpaired.append(tag)
            else:
                mirror[tag] = found
                worst = max(worst, offset)

        if unpaired:
            r, z = at[unpaired[0]]
            problems.append(
                "%d of %d nodes have no mirror partner within %g m, the first "
                "at r = %.6f z = %.6f. meq::GradShafranovSolver refuses this "
                "mesh by name rather than projecting onto something that is "
                "not a reflection"
                % (len(unpaired), len(at), tolerance, r, z))

        involution = sum(1 for tag in mirror
                         if mirror.get(mirror[tag], -1) != tag)
        if involution:
            problems.append(
                "%d nodes reflect to a node that does not reflect back, so the "
                "pairing is not an involution and is therefore not a "
                "reflection" % involution)

        # ---- the elements, and the attribute each image carries ----
        elements = []
        physical_of = {}
        for dimension in (1, 2):
            for _, entity in gmsh.model.getEntities(dimension):
                groups = gmsh.model.getPhysicalGroupsForEntity(dimension,
                                                               entity)
                physical = int(groups[0]) if len(groups) else 0
                kinds, _, blocks = gmsh.model.mesh.getElements(dimension,
                                                               entity)
                for kind, block in zip(kinds, blocks):
                    per = gmsh.model.mesh.getElementProperties(kind)[3]
                    for start in range(0, len(block), per):
                        nodes = tuple(int(value)
                                      for value in block[start:start + per])
                        elements.append((dimension, kind, physical, nodes))
                        physical_of[(dimension, kind,
                                     tuple(sorted(nodes)))] = physical

        stray = 0
        pairs = {}
        for dimension, kind, physical, nodes in elements:
            if any(node not in mirror for node in nodes):
                stray += 1
                continue
            key = (dimension, kind,
                   tuple(sorted(mirror[node] for node in nodes)))
            if key not in physical_of:
                stray += 1
                continue
            pairs.setdefault((dimension, physical), set()).add(
                physical_of[key])
        if stray:
            problems.append(
                "%d of %d elements have no image in the element set, so the "
                "nodes may be paired but the triangulation is not"
                % (stray, len(elements)))
        for (dimension, physical), images in sorted(pairs.items()):
            if len(images) != 1:
                problems.append(
                    "dimension %d attribute %d reflects onto attributes %s -- "
                    "one attribute cannot be the image of a region and of "
                    "part of another"
                    % (dimension, physical, sorted(images)))

        # ---- orientation, which a reflection reverses ----
        negative = 0
        for dimension, kind, physical, nodes in elements:
            if dimension != 2:
                continue
            (r0, z0), (r1, z1), (r2, z2) = (at[node] for node in nodes[:3])
            if 0.5*((r1 - r0)*(z2 - z0) - (r2 - r0)*(z1 - z0)) <= 0.0:
                negative += 1
        if negative:
            problems.append(
                "%d of %d surface elements have a non-positive signed area. A "
                "reflection reverses orientation, so an image element's nodes "
                "have to be permuted; MFEM would integrate these with the "
                "wrong sign" % (negative, sum(1 for e in elements if e[0] == 2)))

        report = dict(nodes=len(at), unpaired=len(unpaired), seam_nodes=seam,
                      elements=len(elements),
                      surface_elements=sum(1 for e in elements if e[0] == 2),
                      involution_failures=involution, unmatched_elements=stray,
                      negative_elements=negative, worst_pairing=worst,
                      tolerance=tolerance,
                      attribute_pairs={key: sorted(value)
                                       for key, value in sorted(pairs.items())})
    finally:
        gmsh.finalize()
    return problems, report


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
    p.add_argument("--vessel", type=float, nargs="+", default=None,
                   metavar="R Z",
                   help="a closed vessel polygon, given as alternating R and Z "
                        "in metres. Everything inside Gamma, outside this and "
                        "not a conductor takes element attribute %d, which "
                        "[source] ExcludeAttributes can then name as a region "
                        "that can never be plasma. The polygon is fragmented "
                        "in, so its edges are mesh faces."
                        % OUTSIDE_VESSEL_ATTRIBUTE)
    p.add_argument("--limiter", type=float, nargs=3, default=None,
                   metavar=("R0", "Z0", "A"),
                   help="a circular limiter to mesh TO, given element "
                        "attribute 20 inside it. A limiter is prescribed input "
                        "and does not move with the solution, so aligning to "
                        "it is the same argument as aligning to a conductor -- "
                        "and it is what lets the contact be FOUND on the curve "
                        "rather than prescribed as a point")
    p.add_argument("--transition", type=float, default=None,
                   help="width of the graded transition out of a refined "
                        "region, metres; defaults to four background sizes")
    p.add_argument("--symmetric", action="store_true",
                   help="mesh the half z >= 0 and REFLECT it, so the result "
                        "is its own mirror image node for node and element "
                        "for element. [solver] UpDownSymmetry is a "
                        "dof-for-dof average and needs one; a triangulation "
                        "of a symmetric geometry is not one, because every "
                        "mesher picks diagonals. The geometry must itself be "
                        "symmetric and an unpaired conductor is refused")
    p.add_argument("--symmetry-check", metavar="FILE", default=None,
                   help="read FILE and report whether it is mirror-symmetric "
                        "about z = 0 -- the question meq's buildMirrorMaps "
                        "asks. Generates nothing and ignores every other "
                        "option; non-zero exit if the mesh is not")
    p.add_argument("--check", action="store_true",
                   help="re-read the written file and assert MEQ's "
                        "preconditions on it; non-zero exit if any fails")
    p.add_argument("-o", "--out", default="halfdisc.msh")
    a = p.parse_args()

    if a.symmetry_check is not None:
        return _report_symmetry(a.symmetry_check)

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

    # THE CONDUCTOR ATTRIBUTES RUN 10, 11, ... AND THE OTHER TWO ARE FIXED
    # NUMBERS IN THAT RANGE, so enough conductors collide with them.  The
    # eleventh coil is attribute 20 and the twenty-first is 30, and an element
    # carrying one of those is then two things at once -- which nothing
    # downstream can tell, `[source] ExcludeAttributes = [ 30 ]` least of all.
    #
    # REFUSED RATHER THAN RENUMBERED.  Moving LIMITER_ATTRIBUTE out of the way
    # would change every mesh that already has one and every config that names
    # it, to buy a case nobody has posed: the machines with twenty conductors
    # are free boundary with no meshed limiter, and the limiter cases have two
    # coils.  The refusal costs those nothing and says what to do.
    coils = a.coil or []
    if a.limiter is not None and len(coils) + 10 > LIMITER_ATTRIBUTE:
        p.error("%d conductors reach attribute %d, which is the limiter's, so "
                "the %dth --coil and the --limiter interior would be the same "
                "element attribute. Mesh the limiter as a coil, or drop "
                "--limiter and recover its faces from the geometry"
                % (len(coils), LIMITER_ATTRIBUTE,
                   LIMITER_ATTRIBUTE - 10 + 1))
    if a.vessel is not None and len(coils) + 10 > OUTSIDE_VESSEL_ATTRIBUTE:
        p.error("%d conductors reach attribute %d, which is the outside of the "
                "vessel, so the %dth --coil and everything beyond the vessel "
                "would be the same element attribute -- and [source] "
                "ExcludeAttributes would then exclude a conductor"
                % (len(coils), OUTSIDE_VESSEL_ATTRIBUTE,
                   OUTSIDE_VESSEL_ATTRIBUTE - 10 + 1))

    if a.limiter is not None:
        lr, lz, la = a.limiter
        if la <= 0.0:
            p.error("a limiter of radius %g has no interior" % la)
        if lr - la <= 0.0:
            p.error("a limiter at R0 = %g a = %g reaches the axis; a closed "
                    "plasma surface through r = 0 carries a non-integrable "
                    "1/r, which is meq::BoundaryShape's own refusal" % (lr, la))
        if math.hypot(lr + la, abs(lz) + la) >= a.rho:
            p.error("a limiter reaching r = %g z = %g is outside Gamma at "
                    "rho = %g; the plasma has to be inside the coupled domain"
                    % (lr + la, abs(lz) + la, a.rho))
        for (r, z, w, h) in coils:
            nearest_r = min(max(lr, r), r + w)
            nearest_z = min(max(lz, z), z + h)
            if math.hypot(nearest_r - lr, nearest_z - lz) < la:
                p.error("the limiter circle ( %g, %g ) a = %g cuts the "
                        "conductor at ( %g, %g ); fragmenting both would make "
                        "a surface that is inside the limiter AND inside a "
                        "coil, which the attributes cannot express"
                        % (lr, lz, la, r, z))

    vessel = None
    if a.vessel is not None:
        if len(a.vessel) < 6 or len(a.vessel) % 2 != 0:
            p.error("--vessel takes alternating R and Z and needs at least "
                    "three points, so an even count of six or more; got %d "
                    "numbers" % len(a.vessel))
        vessel = [(a.vessel[i], a.vessel[i + 1])
                  for i in range(0, len(a.vessel), 2)]
        for (r, z) in vessel:
            if r < 0.0:
                p.error("--vessel reaches r = %g, and the half-disc is r >= 0"
                        % r)
            # AND IT HAS TO FIT INSIDE THE DISC, which is the refusal --coil
            # and --limiter already make and this one did not.  `occ.fragment`
            # keeps the parts of a tool that lie OUTSIDE the shape it is
            # fragmented into, so a vessel reaching past rho adds surfaces to
            # the model that are not the domain: measured, a polygon out to
            # r = 2.5 on rho = 1.5 writes a mesh reaching r = 2.5 with 32
            # nodes carrying Gamma's attribute in the middle of it.  Gamma's
            # attribute is what carries the transferred exterior datum, so
            # that is an interior Dirichlet condition converging at full order
            # to the wrong equilibrium.  The disc is convex, so testing the
            # VERTICES tests the whole outline.
            if math.hypot(r, z) >= a.rho:
                p.error("--vessel reaches ( %g, %g ), on or outside Gamma at "
                        "rho = %g. The vessel is fragmented INTO the disc, so "
                        "the part of it that sticks out becomes mesh of its "
                        "own -- a domain bigger than the one the exterior "
                        "expansion is written for, with Gamma's attribute on "
                        "edges inside it" % (r, z, a.rho))
        if _polygon_area(vessel) <= 0.0:
            p.error("--vessel has zero area: its points are collinear or "
                    "repeated")

    partner = None
    if a.symmetric:
        # THE REFUSALS ARE ABOUT THE MACHINE AND THEY COME FIRST.  Reflecting
        # a geometry that is not its own mirror image produces a mesh of a
        # DIFFERENT machine -- conductors below the midplane replaced by the
        # images of the ones above -- which converges at full order and which
        # nothing downstream can question, since the mesh is all it gets.
        tolerance = SYMMETRY_TOLERANCE*max(a.rho, 1.0)
        try:
            partner, worst = _mirror_partners(coils, tolerance)
        except ValueError as failure:
            p.error(str(failure))
        if a.limiter is not None and a.limiter[1] != 0.0:
            p.error("--limiter names ONE circle and --symmetric would need a "
                    "mirror pair of them, so a limiter at Z0 = %g has no "
                    "partner the option can express. A limiter on the "
                    "midplane, Z0 = 0, is its own" % a.limiter[1])
        if vessel is not None and not _polygon_is_mirrored(vessel, tolerance):
            p.error("--vessel's outline is not its own image under z -> -z "
                    "within %g m, so --symmetric would mesh the upper half of "
                    "one vessel and report the reflection as the other half "
                    "of a vessel nobody described" % tolerance)
        # AND THE SIZE FIELD IS NOT REFUSED, WHICH IS THE LINE WORTH DRAWING
        # CAREFULLY.  --plasma says where elements are SMALL rather than where
        # anything is, so an asymmetric box makes a valid mesh of the right
        # machine; the half that is meshed gets exactly the refinement asked
        # for and the other half gets its reflection.  Said out loud rather
        # than left to be discovered from an element count.
        if (a.plasma is not None
                and abs(2.0*a.plasma[1] + a.plasma[3]) > tolerance):
            print("  note: --plasma z [%.4f, %.4f] is not symmetric about "
                  "z = 0. It is a size field and not geometry, so the machine "
                  "is unaffected; z >= 0 is refined as asked and z < 0 gets "
                  "the reflection of that"
                  % (a.plasma[1], a.plasma[1] + a.plasma[3]))
        if coils:
            print("  symmetric: %d conductors paired about z = 0, %d of them "
                  "their own mirror, worst discrepancy %.3e m"
                  % (len(coils), sum(1 for i, j in enumerate(partner)
                                     if i == j), worst))

    transition = a.transition if a.transition is not None else 4.0 * a.size
    rep = build(a.rho, coils, a.size, a.out, a.order, a.plasma, a.plasma_size,
                a.coil_size, transition, a.limiter, vessel, a.symmetric,
                partner)
    print("  %s: %d elements, %d nodes, %d surfaces (%d coil, %d limiter, "
          "%d outside the vessel), %d arc curves, %d axis curves"
          % (a.out, rep["elements"], rep["nodes"], rep["surfaces"],
             rep["coils"], rep["limiter_surfaces"],
             rep["outside_vessel_surfaces"], rep["arc_curves"],
             rep["axis_curves"]))
    if a.symmetric:
        print("  reflected about z = 0: %d curves on the seam carry no "
              "physical group, %d nodes are shared by both halves. The "
              "surface and curve counts above are the half gmsh meshed; the "
              "elements, the nodes and the conductors are the whole file"
              % (rep["seam_curves"], rep["seam_nodes"]))

    if a.check:
        problems = check(a.out, a.rho, coils)
        for problem in problems:
            print("  FAIL: %s" % problem, file=sys.stderr)
        if problems:
            return 1
        print("  check: r reaches 0 exactly, Gamma and the axis are the outer "
              "boundary and nothing else, every coil covers its rectangle")
        if a.symmetric and _report_symmetry(a.out) != 0:
            return 1
    return 0


def _report_symmetry(path):
    """`--symmetry-check`, and the tail of `--symmetric --check`."""
    problems, rep = check_symmetry(path)
    for problem in problems:
        print("  FAIL: %s" % problem, file=sys.stderr)
    print("  %s: %d nodes (%d unpaired, %d on the seam), %d elements (%d "
          "surface, %d without an image), involution failures %d, "
          "non-positive elements %d, worst pairing %.3e m at a tolerance of "
          "%.3e m"
          % (path, rep["nodes"], rep["unpaired"], rep["seam_nodes"],
             rep["elements"], rep["surface_elements"],
             rep["unmatched_elements"], rep["involution_failures"],
             rep["negative_elements"], rep["worst_pairing"],
             rep["tolerance"]))
    print("  attributes: %s"
          % ", ".join("%dD %d -> %s" % (dimension, physical,
                                        ",".join(str(v) for v in images))
                      for (dimension, physical), images
                      in rep["attribute_pairs"].items()))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
