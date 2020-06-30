#! /usr/bin/env python
"""
This script tests the result of the Galilen method in WarpX.
It compares the energy of the electric field calculated using Galilean method with
'v_galiean = (0.,0., 0.99498743710662)' versus standard PSATD (v_galiean = (0.,0.,0.)):
    * if 'v_galilean == 0': simulation is unstable because of the arosen NCI;
    * if 'v_galilean != 0 : NCI is suppresed => simulation is stable.
"""
import sys
import re
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
from scipy.constants import c, epsilon_0
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# Parse test name and check if current correction (psatd.do_current_correction=1) is applied
current_correction = True if re.search( 'current_correction', filename ) else False

ds = yt.load( filename )

Ex= ds.index.grids[0]['boxlib', 'Ex'].squeeze().v
Ey= ds.index.grids[0]['boxlib', 'Ey'].squeeze().v
Ez= ds.index.grids[0]['boxlib', 'Ez'].squeeze().v

#E field energy calculated with Galilean method (v_galilean = (0,0,0.99498743710662))
energyE_gal_psatd = np.sum(epsilon_0/2*(Ex**2+Ey**2+Ez**2))

#E field energy precalculated with standard PSATD (v_galilean = (0,0,0))
energyE_psatd = 30719.555920

error_rel = energyE_gal_psatd / energyE_psatd
tolerance_rel = 1e-7

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check relative L-infinity spatial norm of div(E) - rho/epsilon_0 when
# current correction (psatd.do_current_correction=1) is applied
if current_correction:
    rho  = ds.index.grids[0]['boxlib','rho' ].squeeze().v
    divE = ds.index.grids[0]['boxlib','divE'].squeeze().v
    error_rel = np.amax( np.abs( divE - rho/epsilon_0 ) ) / np.amax( np.abs( rho/epsilon_0 ) )
    tolerance = 1.e-9
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert( error_rel < tolerance )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
