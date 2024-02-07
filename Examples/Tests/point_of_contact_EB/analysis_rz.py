#!/usr/bin/env python

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

ts_scraping = OpenPMDTimeSeries('./diags/diag2/particles_at_eb/')


it=ts_scraping.iterations
r,theta,z=ts_scraping.get_particle( ['x','theta','z'], species='electron', iteration=it )

r_analytic=0.2
theta_analytic=3.0
z_analytic=0.0000

print('NUMERICAL coordinates of the point of contact:')
print('r=%5.4f, theta=%5.4f, z=%5.4f' % (r[0], theta[0], z[0]))
print('\n')
print('ANALYTICAL coordinates of the point of contact:')
print('x=%5.4f, y=%5.4f, z=%5.4f' % (r_analytic, theta_analytic, z_analytic))

tolerance=0.001
print("tolerance = "+ str(tolerance *100) + '%')

diff_x=np.abs((x[0]-x_analytic)/x_analytic)
diff_y=np.abs((y[0]-y_analytic)/y_analytic)

print("percentage error for x = %5.4f %%" %(diff_x *100))
print("percentage error for y = %5.4f %%" %(diff_y *100))

assert (diff_x < tolerance) and (diff_y < tolerance) and (np.abs(z[0]) < tolerance), 'Test point_of_contact did not pass'
