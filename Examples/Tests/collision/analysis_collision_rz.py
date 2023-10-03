#!/usr/bin/env python3

# Copyright 2019-2021 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the collision module in RZ.
# Electrons are set with the same vr,
# and thus particles are traveling locally in the same direction.
# Physically, the collision rate should therefore be very low.
# Thus the initial electron momentum will be compared with the final momentum.

# Possible errors:
# tolerance: 1.0e-30
# Possible running time: ~ 1.0 s

from glob import glob
import os
import sys

import numpy as np
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 1.0e-15

last_plotfile = sys.argv[1]
if (last_plotfile[-1] == "/"): last_plotfile = last_plotfile[:-1]
plotfile_list = glob(last_plotfile[:-5] + "?????")
opmdfile = './diags/diag2'

for plotfile in plotfile_list:
    # get time index j
    j = int(plotfile[-5:])
    if j==0:
        # load file
        ds = yt.load(plotfile)
        ad = ds.all_data()
        px1 = ad['particle_momentum_x'].to_ndarray()
        py1 = ad['particle_momentum_y'].to_ndarray()
    if j==150:
        # load file
        ds = yt.load(plotfile)
        ad = ds.all_data()
        px2 = ad['particle_momentum_x'].to_ndarray()
        py2 = ad['particle_momentum_y'].to_ndarray()

error = np.max( abs(px1-px2)+abs(py1-py2) )

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, output_file=plotfile, output_format='plotfile', do_particles=False)
checksumAPI.evaluate_checksum(test_name, output_file=opmdfile, output_format='openpmd', do_particles=False)
