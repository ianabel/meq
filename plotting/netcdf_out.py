#!/usr/bin/python3 

from netCDF4 import Dataset
import numpy

nc_root = Dataset("data.nc","w",format="NETCDF4")

N_R = 51
N_Z = 51
# Add dimensions R,z
R_dim = nc_root.createDimension("R", N_R )
Z_dim = nc_root.createDimension("Z", N_Z )

R_variable = nc_root.createVariable("R","f8",("R",))
Z_variable = nc_root.createVariable("Z","f8",("Z",))

R_min = 0.1
R_max = 1
R_variable[:] = numpy.linspace( R_min, R_max, N_R )

Z_min = -1
Z_max = 1
Z_variable[:] = numpy.linspace( Z_min, Z_max, N_Z )

psi_variable = nc_root.createVariable("psi","f8",("R","Z"))


for i in range(0,N_R):
    for j in range(0,N_Z):
        r = R_variable[ i ] 
        z = Z_variable[ j ] 
        psi_variable[ i, j ] = numpy.sqrt( (r - .5)**2 + z ** 2 )

nc_root.close()
