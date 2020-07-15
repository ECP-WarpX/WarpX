#! /usr/bin/env python
"""
This script tests the result of the averaged Galilean PSATD with large
timestep: dz/dx = 2. and c*dt = dz in WarpX.
It compares the energy of the electric field of the uniform plasma calculated
using the averaged Galilean PSATD versus standard Galilean PSATD
with v_galilean = (0.,0.,0.99498743710662):

    * if standard Galilean PSATD is used (psatd.do_time_averaging == 0'):
      simulation is unstable because of the arosen NCI.

    * if averaged Galilean PSATD is used ('psatd.do_time_averaging == 1) :
      NCI is suppresed => simulation is stable.

Accepted relative tolerance: 1e-4
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

#E field energy calculated with averaged Galilean PSATD method (v_galilean = (0,0,0.99498743710662))
energyE_averaged_psatd = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

#E field energy precalculated with standard Galilean PSATD (v_galilean = (0,0,0.99498743710662))
energyE_galilean_psatd = 460315.9845556079

error_rel = energyE_averaged_psatd / energyE_galilean_psatd
tolerance_rel = 1e-4

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
