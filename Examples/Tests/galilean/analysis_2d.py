#! /usr/bin/env python
"""
This script is used to test the results of the Galilean PSATD method and
averaged Galilean PSATD method in WarpX.
It compares the energy of the electric field with precalculated reference energy.
  1) Galilean PSATD test: reference energy was calculated with
     standard PSATD (v_galilean = (0.,0.,0.)):
         * if 'v_galilean == 0': simulation is unstable because of the arosen NCI;
         * if 'v_galilean != 0 : NCI is suppressed => simulation is stable.
  2) Averaged Galilean PSATD with large timestep dz/dx = 4. and c*dt = dz:
     reference energy was calculated with Galilean PSATD (v_galilean = (0.,0.,0.99498743710662):
         * if standard Galilean PSATD is used (psatd.do_time_averaging == 0'):
           simulation is unstable because of the arosen NCI.
         * if averaged Galilean PSATD is used ('psatd.do_time_averaging == 1) :
           NCI is suppressed => simulation is stable.
"""
import sys
import re
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# Parse test name
averaged = True if re.search( 'averaged', filename ) else False
current_correction = True if re.search( 'current_correction', filename ) else False
dims_RZ  = True if re.search('rz', filename) else False

ds = yt.load( filename )

# yt 4.0+ has rounding issues with our domain data:
# RuntimeError: yt attempted to read outside the boundaries
# of a non-periodic domain along dimension 0.
if 'force_periodicity' in dir(ds): ds.force_periodicity()

all_data = ds.covering_grid(level = 0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)
Ex = all_data['boxlib', 'Ex'].squeeze().v
Ey = all_data['boxlib', 'Ey'].squeeze().v
Ez = all_data['boxlib', 'Ez'].squeeze().v

if (averaged):
    # energyE_ref was calculated with Galilean PSATD method (v_galilean = (0,0,0.99498743710662))
    energyE_ref = 32532.00882239954
    tolerance_rel = 1e-6
elif (not dims_RZ and not current_correction):
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.))
    energyE_ref = 35657.99361677053
    tolerance_rel = 1e-8
elif (not dims_RZ and current_correction):
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.)):
    energyE_ref = 35024.02751955393
    tolerance_rel = 2e-8
elif (dims_RZ and not current_correction):
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.))
    energyE_ref = 239019.10670780553
    tolerance_rel = 1e-8
elif (dims_RZ and current_correction):
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.))
    energyE_ref = 471730.0524143545
    tolerance_rel = 1e-9

energyE = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

error_rel = energyE / energyE_ref

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check charge conservation (relative L-infinity norm of error) with current correction
if current_correction:
    divE = all_data['boxlib', 'divE'].squeeze().v
    rho  = all_data['boxlib', 'rho' ].squeeze().v / scc.epsilon_0
    error_rel = np.amax(np.abs(divE - rho)) / max(np.amax(divE), np.amax(rho))
    tolerance = 1e-9
    print("Check charge conservation:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert( error_rel < tolerance )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
