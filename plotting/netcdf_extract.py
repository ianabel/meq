#!/usr/bin/python3

from netCDF4 import Dataset
import matplotlib
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sys

if( len(sys.argv) == 1 ):
    nc_root = Dataset("data.nc","r",format="NETCDF4")
elif ( len(sys.argv) == 2 ):
    nc_root = Dataset(sys.argv[1],"r",format="NETCDF4")
else:
    print("Incorrect number of arguments")
    sys.exit()



R_variable = nc_root.variables["R"]
Z_variable = nc_root.variables["Z"]
psi_variable = nc_root.variables["psi"]

R = numpy.array(R_variable)
Z = numpy.array(Z_variable)
psi = numpy.array(psi_variable)

for i in range(0,3):
    for j in range(0,3):
        print("Psi(",i,",",j,") = ", psi[i][j])
