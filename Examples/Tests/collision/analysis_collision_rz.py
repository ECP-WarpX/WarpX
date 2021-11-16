#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the collision module in RZ.
# Electrons are set with the same vr,
# such that correct collisions will not lead to any velocity changes.
# Thus the initial electron momentum will be compared with the final momentum.

# Possible errors:
# tolerance: 1.0e-30
# Possible running time: ~ 1.0 s

import sys
import yt
import numpy as np
from glob import glob
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 1.0e-15

last_fn = sys.argv[1]
if (last_fn[-1] == "/"): last_fn = last_fn[:-1]
fn_list = glob(last_fn[:-5] + "?????")

for fn in fn_list:
    # get time index j
    j = int(fn[-5:])
    if j==0:
        # load file
        ds = yt.load( fn )
        ad = ds.all_data()
        px1 = ad['particle_momentum_x'].to_ndarray()
        py1 = ad['particle_momentum_y'].to_ndarray()
    if j==150:
        # load file
        ds = yt.load( fn )
        ad = ds.all_data()
        px2 = ad['particle_momentum_x'].to_ndarray()
        py2 = ad['particle_momentum_y'].to_ndarray()

error = np.max( abs(px1-px2)+abs(py1-py2) )

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)

test_name = last_fn[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, do_particles=False)
