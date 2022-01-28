#!/usr/bin/env python

# Copyright 2019-2021 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the particle scraping for the embedded boundary in RZ.
# Particles are initialized between r=0.15 and r=0.2
# having a negative radial velocity.
# A cylindrical embedded surface is placed at r=0.1.
# Upon reaching the surface, particles should be removed.
# At the end of the simulation, i.e., at time step 37,
# there should be 512 particles left.

# Possible errors: 0
# tolerance: 0
# Possible running time: < 1 s

import os
import sys

import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 0

fn = sys.argv[1]
ds = yt.load( fn )
ad = ds.all_data()
x = ad['electron', 'particle_position_x'].v

error = len(x)-512
print('error = ', error)
print('tolerance = ', tolerance)
assert(error==tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, do_particles=False)
