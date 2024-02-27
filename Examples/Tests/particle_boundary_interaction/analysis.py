#!/usr/bin/env python

"""
This script tests the last coordinate after adding an electron.
The sphere is centered on O and has a radius of 0.2 (EB)
The electron is initially at: (0,0,-0.25) and moves with a velocity: 
(0.5e10,0,1.0e10) with a time step of 1e-11.
An input file PICMI_inputs_rz.py is used.
"""
import os
import sys

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
import yt

yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Open plotfile specified in command line
filename = sys.argv[1]
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename, output_format='openpmd')

ts = OpenPMDTimeSeries('./diags/diag1/')

it=ts.iterations
x,y,z=ts.get_particle(['x','y','z'], species='electrons', iteration=it[-1])

# Analytical results calculated
x_analytic=0.03532
y_analytic=0.00000
z_analytic=-0.20531

print('NUMERICAL coordinates of the point of contact:')
print('x=%5.5f, y=%5.5f, z=%5.5f' % (x[0], y[0], z[0]))
print('\n')
print('ANALYTICAL coordinates of the point of contact:')
print('x=%5.5f, y=%5.5f, z=%5.5f' % (x_analytic, y_analytic, z_analytic))

tolerance=1e-5

diff_x=np.abs((x[0]-x_analytic)/x_analytic)
diff_z=np.abs((z[0]-z_analytic)/z_analytic)

print('\n')
print("percentage error for x = %5.4f %%" %(diff_x *100))
print("percentage error for z = %5.4f %%" %(diff_z *100))

assert (diff_x < tolerance) and (y[0] < 1e-8) and (diff_z < tolerance), 'Test particle_boundary_interaction did not pass'
