import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.optimize import curve_fit


ts_scraping = OpenPMDTimeSeries('../../../../warpx/Examples/Tests/point_of_contact_EB/diags/diag3/particles_at_eb/')
x,y,z=ts_scraping.get_particle( ['x','y','z'], species='electron', iteration=it_stamp )

x_analytic=-0.1983
y_analytic=0.02584
z_analytic=0.0000

print('NUMERICAL coordinates of the point of contact:')
print('x=%5.4f, y=%5.4f, z=%5.4f' % (x, y, z))
print('\n')
print('ANALYTICAL coordinates of the point of contact:')
print('x=%5.4f, y=%5.4f, z=%5.4f' % (x_analytic, y_analytic, z_analytic))

tolerance=0.001
print("tolerance = "+ str(tolerance *100) + '%')

diff_x=np.abs((x-x_analytic)/x_analytic)
diff_y=np.abs((y-y_analytic)/y_analytic)
diff_z=np.abs((z-z_analytic)/z_analytic)

print("pourcentage error for x = "+ str(diff_x *100) + '%')
print("pourcentage error for y = "+ str(diff_y*100) + '%')
print("pourcentage error for z = "+ str(diff_z*100) + '%')

assert (diff_x < tolerance) and (diff_y < tolerance) and (diff_z < tolerance), 'Test point_of_contact did not pass'


