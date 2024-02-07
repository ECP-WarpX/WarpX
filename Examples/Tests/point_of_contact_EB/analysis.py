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
timestamp,x,y,z,nx,ny,nz=ts_scraping.get_particle( ['timestamp','x','y','z','nx','ny','nz'], species='electron', iteration=it )

timestamp_analytic=3.58
x_analytic=-0.1983
y_analytic=0.02584
z_analytic=0.0000
nx_analytic=-0.99
ny_analytic=0.13
nz_analytic=0.0

timestamp=timestamp[0]*1e10
print('\n')
print('NUMERICAL coordinates of the point of contact:')
print('time_stamp=%5.4f e-10, x=%5.4f, y=%5.4f, z=%5.4f, nx=%5.4f, ny=%5.4f, nz=%5.4f' % (timestamp, x[0], y[0], z[0], nx[0], ny[0], nz[0]))
print('\n')
print('ANALYTICAL coordinates of the point of contact:')
print('time_stamp=%5.2f, x=%5.4f, y=%5.4f, z=%5.4f, nx=%5.4f, ny=%5.4f, nz=%5.4f' % (timestamp_analytic, x_analytic, y_analytic, z_analytic, nx_analytic, ny_analytic, nz_analytic))

tolerance=0.001
tolerance_t=0.003
tolerance_n=0.01
print('\n')
print("tolerance = "+ str(tolerance *100) + '%')
print("tolerance for the time_stamp = "+ str(tolerance_t *100) + '%')
print("tolerance for the normal components = "+ str(tolerance_n *100) + '%')
diff_timestamp=np.abs((timestamp-timestamp_analytic)/timestamp_analytic)
diff_x=np.abs((x[0]-x_analytic)/x_analytic)
diff_y=np.abs((y[0]-y_analytic)/y_analytic)
diff_nx=np.abs((nx[0]-nx_analytic)/nx_analytic)
diff_ny=np.abs((ny[0]-ny_analytic)/ny_analytic)

print('\n')
print("percentage error for timestamp = %5.4f %%" %(diff_timestamp *100))
print("percentage error for x = %5.4f %%" %(diff_x *100))
print("percentage error for y = %5.4f %%" %(diff_y *100))
print("percentage error for nx = %5.2f %%" %(diff_nx *100))
print("percentage error for ny = %5.2f %%" %(diff_ny *100))
print("nz = %5.2f " %(nz[0]))


assert (diff_timestamp < tolerance_t) and (diff_x < tolerance) and (diff_y < tolerance) and (diff_nx < tolerance_n) and (diff_ny < tolerance_n) and (np.abs(nz) < tolerance_n), 'Test point_of_contact did not pass'
