<<<<<<< HEAD
import numpy as np
import math
=======
>>>>>>> ef1e962595decf0b2be4c24004b44886094e2cc5
import matplotlib.pyplot as plt
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.optimize import curve_fit

ts_scraping = OpenPMDTimeSeries('../../../../warpx/Examples/Physics_applications/point_of_contact_EB/diags/diag3/particles_at_eb/')
ts_full = OpenPMDTimeSeries('../../../../warpx/Examples/Physics_applications/point_of_contact_EB/diags/diag2/')
it_stamp=ts_scraping.get_particle( ['timestamp'], species='electron' )
x,y,z=ts_full.get_particle( ['x','y','z'], species='electron', iteration=it_stamp )




print('Coordinates of the point of contact:')
print('x=%5.3f, y=%5.9f, z=%5.9f' % (x, y, z))
x_round = round(x, 3)
y_cut = float('%5.4f' %y)

assert (x_round!=-0.1983) and (y_cut!=0.0258), 'Test point of contact at EB did not pass'
