#!/usr/bin/env python3

# Copyright 2022 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This script tests the particle pusher (HC)
# using a force-free field,
# in which position x should remain 0.
# An initial velocity Vy corresponding to
# Lorentz factor = 20 is used.
# Bz is fixed at 1 T.
# Ex = -Vy*Bz.

# Possible errors:
# Boris:     2321.3958529
# Vay:       0.00010467
# HC:        0.00011403
# tolerance: 0.001
# Possible running time: ~ 4.0 s

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries
import scipy.constants as sc

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 0.001
x0 = 1.140180338454055003e-1
y0 = -4.908434106235705363e-2
z0 = 4.379090207432996706

filename = sys.argv[1]
print(filename)
ts = OpenPMDTimeSeries(filename)
it = ts.iterations
x,y,z = ts.get_particle(var_list=['x','y','z','ux','uy','uz'], iteration=it[len(it)-1], species='proton')

error = np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
