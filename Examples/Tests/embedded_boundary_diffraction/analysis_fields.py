#!/usr/bin/env python3
"""
This test checks the implementation of the embedded boundary in cylindrical geometry,
by checking the diffraction of a laser by an embedded boundary here.
We then check that the first minimum of the diffracted intensity pattern
occurs along the angle given by the theoretical Airy pattern, i.e.
theta_diffraction = 1.22 * lambda / d
"""
import os
import sys

import numpy as np
from scipy.constants import c, mu_0, pi
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

ts = OpenPMDTimeSeries('./diags/diag1/')

# Extract the intensity as a function of r and z
Ex, info = ts.get_field('E', 'x', iteration=300)
I = gaussian_filter1d(Ex**2, sigma=5, axis=0)
plt.imshow( I, extent=info.imshow_extent, origin='lower')
irmax = np.argmax( I, axis=-1)

# Find the radius of the first minimum, as a function of z
def r_first_minimum(iz):
    ir = len(info.r)//2
    while I[iz, ir+1] < I[iz, ir]:
        ir += 1
    return info.r[ir]
r = np.array([ r_first_minimum(iz) for iz in range(len(info.z)) ])

# Check that this corresponds to the prediction from the Airy pattern
theta_diffraction = np.arcsin(1.22*0.1/0.4)/2
assert np.all( abs(r[50:] - theta_diffraction*info.z[50:]) < 0.03 )

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename, )
