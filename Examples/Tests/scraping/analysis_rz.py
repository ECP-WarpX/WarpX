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
# In addition, the test checks the boundary scraping diagnostic
# by making sure that all removed particles are properly recorded.

# Possible errors: 0
# tolerance: 0
# Possible running time: < 1 s

import os
import sys

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
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

# Check that all the removed particles are properly recorded
# by making sure that, at each iteration, the sum of the number of
# remaining particles and scraped particles is equal to the
# original number of particles
ts_full = OpenPMDTimeSeries('./diags/diag2/')
ts_scraping = OpenPMDTimeSeries('./diags/diag3/')

def n_remaining_particles( iteration ):
    w, = ts_full.get_particle(['w'], iteration=iteration)
    return len(w)
def n_scraped_particles( iteration ):
    timestamp = ts_scraping.get_particle( ['timestamp'] )
    return (timestamp <= iteration).sum()
n_remaining = np.array([ n_remaining_particles(iteration) for iteration in ts_full.iterations ])
n_scraped = np.array([ n_scraped_particles(iteration) for iteration in ts_full.iterations ])
n_total = n_remaining[0]
assert np.all( n_scraped+n_remaining == n_total)

# Checksum test
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, do_particles=False)
