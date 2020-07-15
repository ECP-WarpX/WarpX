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


ds = yt.load( filename )

averaged = True if re.search( 'averaged', filename ) else False

Ex= ds.index.grids[0]['boxlib', 'Ex'].squeeze().v
Ey= ds.index.grids[0]['boxlib', 'Ey'].squeeze().v
Ez= ds.index.grids[0]['boxlib', 'Ez'].squeeze().v

if (averaged):
    # energyE_ref was calculated with Galilean PSATD method (v_galilean = (0,0,0.99498743710662))
    energyE_ref = 26913.546573259937
    tolerance_rel = 1e-5
else:
    # energyE_ref was calculated with standard PSATD method (v_galilean = (0.,0.,0.))
    energyE_ref = 38362.88743899688
    tolerance_rel = 1e-8;

energyE = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

error_rel = energyE / energyE_ref

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
