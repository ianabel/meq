#!/usr/bin/python3

from netCDF4 import Dataset
import matplotlib
import numpy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sys

if( len(sys.argv) == 1 ):
    config_file = 'meq.conf'
elif ( len(sys.argv) == 2 ):
    config_file = argv[1]
else:
    print("Incorrect number of arguments")
    sys.exit()

import toml
config = toml.load(config_file)
    
nc_root = Dataset(config['options']['NetCDFFile'],"r",format="NETCDF4")

fig, ax = plt.subplots()

from matplotlib.patches import Polygon
vacuum_vessel = config['domain']['VacuumVessel']
vv_poly = Polygon( vacuum_vessel, True, fill=False )
ax.add_patch( vv_poly )

ax.set_ylim(0,1.1)
ax.set_xlim(-1.5,1.5)


R_variable = nc_root.variables["R"]
Z_variable = nc_root.variables["Z"]
psi_variable = nc_root.variables["psi"]

R = numpy.array(R_variable)
Z = numpy.array(Z_variable)
psi = numpy.array(psi_variable)

psi_mask = numpy.zeros_like( psi, dtype=bool )

nPts = R.size * Z.size
nMasked = 0

vacuum_vessel_polygon = Polygon( vacuum_vessel, True )
for i in range(0,R.size):
    for j in range(0,Z.size):
        if not vacuum_vessel_polygon.contains_point((Z[j],R[i]),1e-6):
            psi_mask[ i ][ j ] = True
            nMasked += 1

psi_masked = numpy.ma.array( psi, mask=psi_mask )

from matplotlib.patches import Rectangle

for c in config['coils']:
    Z_c = c['Z']
    R_c = c['R']
    w = c['Width']
    h = c['Height']
    box_patch = Rectangle( (Z_c-h/2,R_c-w/2), h, w, facecolor='red', alpha=.7 )
    ax.add_patch( box_patch )
    

CS = ax.contour(Z, R, psi_masked, 25)

ax.clabel(CS, inline=1, fontsize=10)

ax.set_title('Magnetic Flux Contours')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.gca().set_aspect('equal', adjustable='box')

plt.show()
