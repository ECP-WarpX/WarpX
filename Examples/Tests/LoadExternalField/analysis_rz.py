#!/usr/bin/env python3

# Copyright 2022 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os
import sys

import numpy as np
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 0.001
r0 = 0.12402005
z0 = 4.3632492

filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()
r  = ad['proton','particle_position_x'].to_ndarray()
z  = ad['proton','particle_position_y'].to_ndarray()

error = np.min(np.sqrt((r-r0)**2+(z-z0)**2))

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
