#! /usr/bin/env python

# Copyright 2021
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
"""
This tests verifies the Silver-Mueller boundary condition.
A laser pulse is emitted and propagates towards the boundaries ; the
test check that the reflected field at the boundary is negligible.
"""

import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

ds = yt.load( filename )
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()
# The peak of the initial laser pulse is on the order of 6 V/m
# Check that the amplitude after reflection is less than 0.01 V/m
max_reflection_amplitude = 0.01
assert np.all( abs(Ex) < max_reflection_amplitude )
assert np.all( abs(Ey) < max_reflection_amplitude )
assert np.all( abs(Ez) < max_reflection_amplitude )

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
