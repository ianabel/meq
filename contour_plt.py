#!/usr/bin/python3

from netCDF4 import Dataset
import matplotlib
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sys

nc_root = Dataset("contours.nc","r",format="NETCDF4")


R_variable = nc_root.variables["R"]
Z_variable = nc_root.variables["Z"]
psi_variable = nc_root.variables["Psi"]

psi_values = numpy.array(psi_variable)


for i in range(0,psi_values.size):
    psi_val = psi_values[i]
    Rvec = numpy.asarray(R_variable[i])
    Zvec = numpy.asarray(Z_variable[i])
    if( Rvec.shape[0] != Zvec.shape[0] ):
        sys.exit(-1)
    nPts = Rvec.shape[0]
    print("")
    print("# Psi = ",psi_val," / #pts = ", nPts)
    for j in range(0,nPts):
        print(Rvec[j],"\t",Zvec[j])
    print("")



  
