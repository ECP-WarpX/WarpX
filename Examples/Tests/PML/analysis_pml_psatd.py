#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

############################
### INITIAL LASER ENERGY ###
############################
energy_start = 7.297301362346945e-08 # electromagnetic energy at iteration 50

##########################
### FINAL LASER ENERGY ###
##########################
ds = yt.load( filename )
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Bx = all_data_level_0['boxlib', 'Bx'].v.squeeze()
By = all_data_level_0['boxlib', 'By'].v.squeeze()
Bz = all_data_level_0['boxlib', 'Bz'].v.squeeze()
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()
rho = all_data_level_0['boxlib','rho'].v.squeeze()
divE = all_data_level_0['boxlib','divE'].v.squeeze()
energyE = np.sum(0.5 * scc.epsilon_0 * (Ex**2 + Ey**2 + Ez**2))
energyB = np.sum(0.5 / scc.mu_0 * (Bx**2 + By**2 + Bz**2))
energy_end = energyE + energyB

reflectivity     = energy_end / energy_start
reflectivity_max = 1e-06

print("reflectivity     = " + str(reflectivity))
print("reflectivity_max = " + str(reflectivity_max))

assert(reflectivity < reflectivity_max)

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
