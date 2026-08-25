#!/usr/bin/python3

from netCDF4 import Dataset
import matplotlib
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sys

nc_root1 = Dataset(sys.argv[1],"r",format="NETCDF4")
nc_root2 = Dataset(sys.argv[2],"r",format="NETCDF4")



R1_variable = nc_root1.variables["R"]
Z1_variable = nc_root1.variables["Z"]
psi1_variable = nc_root1.variables["psi"]

R2_variable = nc_root2.variables["R"]
Z2_variable = nc_root2.variables["Z"]
psi2_variable = nc_root2.variables["psi"]

i = int(sys.argv[3])
j = int(sys.argv[4])

psi1 = psi1_variable[i][j]
psi2 = psi2_variable[i][j]

print("Psi_1(",R1_variable[i],",",Z1_variable[j],") = ", psi1)
print("Psi_2(",R2_variable[i],",",Z2_variable[j],") = ", psi2)
print("Relative Error is ", abs(psi1 - psi2)/abs(max(psi1,psi2)))
