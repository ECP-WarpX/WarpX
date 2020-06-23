#! /usr/bin/env python
"""
This script tests the result of the Galilen method in WarpX.
It compares the energy of the electric field calculated using Galilean method with
'v_galiean = (0.,0., 0.99498743710662)' versus standard PSATD (v_galiean = (0.,0.,0.)):
    * if 'v_galilean == 0': simulation is unstable because of the arosen NCI;
    * if 'v_galilean != 0 : NCI is suppresed => simulation is stable.
"""
import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]


ds = yt.load( filename )

Ex= ds.index.grids[0]['boxlib', 'Ex'].squeeze().v
Ey= ds.index.grids[0]['boxlib', 'Ey'].squeeze().v
Ez= ds.index.grids[0]['boxlib', 'Ez'].squeeze().v

#E field energy calculated with Galilean method (v_galilean = (0,0,0.99498743710662))
energyE_gal_psatd = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

#E field energy precalculated with standard PSATD (v_galilean = (0,0,0))
energyE_psatd = 30719.555920

error_rel = energyE_gal_psatd / energyE_psatd
tolerance_rel = 1e-7

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
