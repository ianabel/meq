#!/usr/bin/env python3
"""A half-disc reaching the axis, with rectangular coils MESHED TO, for FB-6.

WHY A MESHING LIBRARY AT ALL.  Everything MEQ has meshed so far came from
`mfem::Mesh::MakeCartesian2D` -- a rectangle of triangles -- with curved and
non-rectangular geometry reached by the extension technique on top of it.  Free
boundary needs something that grid cannot be: a SEMICIRCLE centred on the axis
(the exterior expansion of FREE-BOUNDARY-PLAN.md section 3 is only valid there),
reaching r = 0 exactly (FB-A), with the coil rectangles as their own subdomains.

WHY GMSH, AND WHY THIS IS NOT A RECOMMENDATION BUT A DEMONSTRATION.
MFEM reads Gmsh's format natively -- `Mesh::ReadGmshMesh`, no converter -- and
gmsh has a Python and a C++ API, so the geometry is a program rather than a file
somebody edited.  Its PHYSICAL GROUPS arrive in MFEM as element attributes and
boundary attributes, which is exactly the two things this needs: one attribute
per coil so a source can be restricted to it, and separate boundary attributes
for the arc (Gamma) and the flat side (the axis).  Verified end to end:

    surfaces 3   coil attrs 10, 11   plasma attr 1
    boundary attrs 1 (arc), 2 (axis)
    attr 10:  r [0.9000, 1.1500]  z [-0.8600, -0.7400]
    minimum r over the mesh: 0.000e+00

THE COILS ARE MESHED TO RATHER THAN CUT, and that is FB-2's finding applied
rather than a convenience.  A conductor's geometry is PRESCRIBED INPUT -- it
does not move with the solution -- so it can always be aligned to, and aligning
it took the measured rates from 1.33/1.27/1.09 to 1.99/2.88/3.01.  Only the
plasma support genuinely moves, and FB-4 measures what that costs.  `occ.fragment`
is what does it: it splits the disc along the coil rectangles so their edges are
mesh edges.

Nothing here is MEQ-specific beyond the defaults; it is a plain gmsh script.

    python3 halfdisc.py --rho 1.5 --h 0.1 --coil 0.9 0.8 0.25 0.12 -o machine.msh
"""
import argparse
import sys

import gmsh


def build(rho, coils, size, out, order=1):
    """A half-disc of radius `rho` about the origin, r >= 0, with `coils` a list
    of (rmin, zmin, width, height) rectangles fragmented into it."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("halfdisc")
    occ = gmsh.model.occ

    disc = occ.addDisk(0, 0, 0, rho, rho)
    # Cut the r < 0 half away.  A disc plus a cut is used rather than a circle
    # arc plus two lines because OCC then owns the curvature, and gmsh will put
    # high-order nodes ON the true arc if asked for order > 1.
    box = occ.addRectangle(-rho, -rho, 0, rho, 2 * rho)
    half, _ = occ.cut([(2, disc)], [(2, box)])

    tags = [occ.addRectangle(r, z, 0, w, h) for (r, z, w, h) in coils]
    if tags:
        occ.fragment(half, [(2, t) for t in tags])
    occ.synchronize()

    surfaces = [t[1] for t in gmsh.model.getEntities(2)]
    centres = {s: occ.getCenterOfMass(2, s) for s in surfaces}

    def in_coil(s):
        x, y, _ = centres[s]
        return any(r <= x <= r + w and z <= y <= z + h for (r, z, w, h) in coils)

    coil_surfaces = sorted(s for s in surfaces if in_coil(s))
    vacuum = [s for s in surfaces if s not in coil_surfaces]

    # Attribute 1 is everything that is not a coil.  Coils take 10, 11, ...
    # which leaves room and keeps them obviously distinct in a .mesh dump.
    gmsh.model.addPhysicalGroup(2, vacuum, 1)
    gmsh.model.setPhysicalName(2, 1, "vacuum")
    for i, s in enumerate(coil_surfaces):
        gmsh.model.addPhysicalGroup(2, [s], 10 + i)
        gmsh.model.setPhysicalName(2, 10 + i, "coil%d" % i)

    # THE TWO BOUNDARIES ARE NOT INTERCHANGEABLE and must not share an
    # attribute: Gamma carries the transferred exterior datum and the axis
    # carries psi = 0, which is a plain essential condition.  They are told
    # apart by geometry -- the axis is the only boundary at r = 0 -- rather than
    # by tag order, which OCC does not promise.
    arc, axis = [], []
    for _, e in gmsh.model.getEntities(1):
        x, _, _ = occ.getCenterOfMass(1, e)
        (axis if abs(x) < 1e-9 else arc).append(e)
    gmsh.model.addPhysicalGroup(1, arc, 1)
    gmsh.model.setPhysicalName(1, 1, "Gamma")
    gmsh.model.addPhysicalGroup(1, axis, 2)
    gmsh.model.setPhysicalName(1, 2, "axis")

    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size / 8.0)
    # 2.2 rather than 4.1: MFEM's reader handles both, and 2.2 is the one whose
    # ASCII is readable when something needs looking at by eye.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(2)
    if order > 1:
        gmsh.model.mesh.setOrder(order)
    gmsh.write(out)

    report = dict(surfaces=len(surfaces), coils=len(coil_surfaces),
                  arc_edges=len(arc), axis_edges=len(axis))
    gmsh.finalize()
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rho", type=float, default=1.5,
                   help="radius of the semicircle Gamma, metres")
    p.add_argument("--size", type=float, default=0.15, help="target element size")
    p.add_argument("--order", type=int, default=1,
                   help="geometric order; >1 curves the arc")
    p.add_argument("--coil", type=float, nargs=4, action="append",
                   metavar=("RMIN", "ZMIN", "WIDTH", "HEIGHT"), default=None,
                   help="a rectangular coil, repeatable")
    p.add_argument("-o", "--out", default="halfdisc.msh")
    a = p.parse_args()

    coils = a.coil or []
    for (r, z, w, h) in coils:
        if r <= 0.0:
            p.error("a coil at RMIN = %g reaches or crosses the axis; the "
                    "operator's 1/r is not integrable through r = 0, which is "
                    "the same refusal meq::Coil and [[coils]] make" % r)
    rep = build(a.rho, coils, a.size, a.out, a.order)
    print("  %s: %d surfaces (%d coil), %d arc edges, %d axis edges"
          % (a.out, rep["surfaces"], rep["coils"], rep["arc_edges"],
             rep["axis_edges"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
