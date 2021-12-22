#!/usr/bin/env python3
"""
This script is used to test the results of the Galilean PSATD method and
averaged Galilean PSATD method in WarpX.
It compares the energy of the electric field with precalculated reference energy.
  1) Galilean PSATD test: reference energy was calculated with
     standard PSATD (v_galilean = (0.,0.,0.)):
         * if 'v_galilean == 0': simulation is unstable because of the arosen NCI;
         * if 'v_galilean != 0 : NCI is suppressed => simulation is stable.
  2) Averaged Galilean PSATD with large timestep dz/dx = 3. and c*dt = dz:
     reference energy was calculated with Galilean PSATD (v_galilean = (0.,0.,0.99498743710662):
         * if standard Galilean PSATD is used (psatd.do_time_averaging == 0'):
           simulation is unstable because of the arosen NCI.
         * if averaged Galilean PSATD is used ('psatd.do_time_averaging == 1) :
           NCI is suppressed => simulation is stable.
"""
import os
import re
import sys

import numpy as np
import scipy.constants as scc

import yt ; yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# Parse test name
averaged = True if re.search( 'averaged', filename ) else False
current_correction = True if re.search( 'current_correction', filename ) else False

ds = yt.load( filename )

all_data = ds.covering_grid(level = 0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)
Ex = all_data['boxlib', 'Ex'].squeeze().v
Ey = all_data['boxlib', 'Ey'].squeeze().v
Ez = all_data['boxlib', 'Ez'].squeeze().v

if (averaged):
    # energyE_ref was calculated with Galilean PSATD method (v_galilean = (0,0,0.99498743710662))
    energyE_ref = 6.816182771544472
    tolerance_rel = 1e-4
elif (current_correction):
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.)):
    energyE_ref = 75333.81851879464
    tolerance_rel = 5e-8;
else:
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.))
    energyE_ref = 8218.678808709019
    tolerance_rel = 1e-6;

energyE = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

error_rel = energyE / energyE_ref

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check charge conservation (relative L-infinity norm of error) with current correction
if current_correction:
    rho  = all_data['boxlib', 'rho' ].squeeze().v
    divE = all_data['boxlib', 'divE'].squeeze().v
    error_rel = np.amax( np.abs( divE - rho/scc.epsilon_0 ) ) / np.amax( np.abs( rho/scc.epsilon_0 ) )
    tolerance = 1e-9
    print("Check charge conservation:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert( error_rel < tolerance )

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
