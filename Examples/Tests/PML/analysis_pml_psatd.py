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
energy_start = 7.297324456269243e-08 # electromagnetic energy at iteration 50

# Check consistency of field energy diagnostics with initial energy above
ds = yt.load('pml_x_psatd_plt00050')
all_data_level_0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Bx = all_data_level_0['boxlib', 'Bx'].v.squeeze()
By = all_data_level_0['boxlib', 'By'].v.squeeze()
Bz = all_data_level_0['boxlib', 'Bz'].v.squeeze()
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()
energyE = np.sum(0.5 * scc.epsilon_0 * (Ex**2 + Ey**2 + Ez**2))
energyB = np.sum(0.5 / scc.mu_0 * (Bx**2 + By**2 + Bz**2))
energy_start_diags = energyE + energyB
error = abs(energy_start - energy_start_diags) / energy_start
tolerance = 1e-14
print("energy_start expected = " + str(energy_start))
print("energy_start_diags    = " + str(energy_start_diags))
print("relative error    = " + str(error))
assert (error < tolerance)

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
energyE = np.sum(0.5 * scc.epsilon_0 * (Ex**2 + Ey**2 + Ez**2))
energyB = np.sum(0.5 / scc.mu_0 * (Bx**2 + By**2 + Bz**2))
energy_end = energyE + energyB

reflectivity = energy_end / energy_start
reflectivity_max = 1e-06

print("reflectivity     = " + str(reflectivity))
print("reflectivity_max = " + str(reflectivity_max))

assert(reflectivity < reflectivity_max)

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
