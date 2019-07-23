#!/usr/bin/python3
#

import toml
import gmsh
gmsh.initialize("-format vtk")


conf_dict = toml.load("meq.conf")

coil_list = conf_dict['coils']

points = []


geo = gmsh.model.geo
model = gmsh.model

h = 0.25
coil_h = .05

RMin = conf_dict['domain']['RMin']
RMax = conf_dict['domain']['RMax']
ZMin = conf_dict['domain']['ZMin']
ZMax = conf_dict['domain']['ZMax']

LowerLeft  = geo.addPoint(RMin,ZMin,0,h)
LowerRight = geo.addPoint(RMax,ZMin,0,h)
UpperRight = geo.addPoint(RMax,ZMax,0,h)
UpperLeft  = geo.addPoint(RMin,ZMax,0,h)

domainCurve = []
domainCurve.append( geo.addLine(LowerLeft,LowerRight) )
domainCurve.append( geo.addLine(LowerRight,UpperRight) )
domainCurve.append( geo.addLine(UpperRight,UpperLeft) )
domainCurve.append( geo.addLine(UpperLeft,LowerLeft) )

domainCurveTag = geo.addCurveLoop(domainCurve)

coilCurveTags = []
coilSurfaceTags = []
for coil in coil_list:
    R = coil['R']
    Z = coil['Z']
    w = coil['Width']
    h = coil['Height']

    RMin = R - w/2
    RMax = R + w/2
    ZMin = Z - h/2
    ZMax = Z + h/2
    delta = min(h,w)/4

    LowerLeft  = geo.addPoint(RMin,ZMin,0,coil_h)
    LowerRight = geo.addPoint(RMax,ZMin,0,coil_h)
    UpperRight = geo.addPoint(RMax,ZMax,0,coil_h)
    UpperLeft  = geo.addPoint(RMin,ZMax,0,coil_h)

    coilCurve = []
    coilCurve.append( geo.addLine(LowerLeft,LowerRight) )
    coilCurve.append( geo.addLine(LowerRight,UpperRight) )
    coilCurve.append( geo.addLine(UpperRight,UpperLeft) )
    coilCurve.append( geo.addLine(UpperLeft,LowerLeft) )
    coilCurveTag = geo.addCurveLoop(coilCurve)

    coilCurveTags.append(coilCurveTag)

    coilSurfaceTags.append(geo.addPlaneSurface([coilCurveTag]))

domainSurfaceTag = geo.addPlaneSurface([domainCurveTag] + coilCurveTags)

model.addPhysicalGroup(2,coilSurfaceTags + [domainSurfaceTag],1)
model.mesh.generate(2)

filename = conf_dict['options']['MeshFile']
gmsh.write(filename)
gmsh.finalize()
