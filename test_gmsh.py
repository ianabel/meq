#!/usr/bin/python3
#

import toml
import gmsh

foo = []

gmsh.initialize("-format vtk")

geo = gmsh.model.geo
model = gmsh.model

h = 0.2
coil_h = .01

geo.addPoint(.05,-1.5,0,h,1)
geo.addPoint(1.5,-1.5,0,h,2)
geo.addPoint(1.5, 1.5,0,h,3)
geo.addPoint(.05, 1.5,0,h,4)

geo.addLine(1,2,1)
geo.addLine(2,3,2)
geo.addLine(3,4,3)
geo.addLine(4,1,4)

geo.addPoint(0.975,0.975,0,coil_h,5);
geo.addPoint(1.025,0.975,0,coil_h,6);
geo.addPoint(1.025,1.025,0,coil_h,7);
geo.addPoint(0.975,1.025,0,coil_h,8);

geo.addLine(5,6,5)
geo.addLine(6,7,6)
geo.addLine(7,8,7)
geo.addLine(8,5,8)


geo.addPoint(0.975,-0.975,0,coil_h,9);
geo.addPoint(1.025,-0.975,0,coil_h,10);
geo.addPoint(1.025,-1.025,0,coil_h,11);
geo.addPoint(0.975,-1.025,0,coil_h,12);

geo.addLine(9,10,9)
geo.addLine(10,11,10)
geo.addLine(11,12,11)
geo.addLine(12,9,12)

geo.addCurveLoop([1,2,3,4],1)
geo.addCurveLoop([5,6,7,8],2)
geo.addCurveLoop([9,10,11,12],3)

geo.addPlaneSurface([1,2,3],1)
geo.addPlaneSurface([2],2)
geo.addPlaneSurface([3],3)

model.addPhysicalGroup(2,[1,2,3],1)
model.mesh.generate(2)

gmsh.write("test_gmsh.vtk")
gmsh.finalize()






