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

diff_r=np.abs((r[0]-r_analytic)/r_analytic)
diff_theta=np.abs((theta[0]-theta_analytic)/theta_analytic)

print("percentage error for r = %5.4f %%" %(diff_r *100))
print("percentage error for theta = %5.4f %%" %(diff_theta *100))

assert (diff_r < tolerance) and (diff_theta < tolerance) and (np.abs(z[0]) < tolerance), 'Test point_of_contact did not pass'
