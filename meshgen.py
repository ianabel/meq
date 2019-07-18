#!/usr/bin/python3
#

import toml
import gmsh

conf_dict = toml.load("meq.conf")

coil_list = conf_dict['coils']

points = []

for coil in coil_list:
    R = coil['R']
    Z = coil['Z']
    w = coil['Width']
    h = coil['Height']

    x_l = R - w/2
    x_u = R + w/2
    y_l = Z - h/2
    y_u = Z + h/2
    delta = min(h,w)/4



