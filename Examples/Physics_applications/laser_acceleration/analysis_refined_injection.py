#!/usr/bin/env python3

# Copyright 2019 Andrew Myers
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the "warpx.refine_plasma=1" option by comparing
# the actual number of electrons at step 200 to the expected value
import os
import sys

import yt

yt.funcs.mylog.setLevel(50)
import numpy as np

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Read the file
ds = yt.load(fn)

# count the number of particles
ad = ds.all_data()
np = ad['electrons', 'particle_id'].size

# the number of coarse particle streams
n_coarse = 10

# the number of fine particle streams
n_fine = 64

# num particles per stream at time 0
n_0 = 15

# number of times moving window moves
n_move = 99

# ref ratio
rr = 2

np_expected = (n_coarse + n_fine*rr)*(n_0 + n_move)

assert( np == np_expected )

# Test uniformity of rho, by taking a slice of rho that
# crosses the edge of the refined injection region
# (but is ahead of the mesh refinement patch)
ds.force_periodicity()
ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
rho = ad['rho'].to_ndarray().squeeze()
rho_slice = rho[13:51, 475]
# Test uniformity up to 0.5% relative variation
assert( rho_slice.std() < 0.005*abs(rho_slice.mean()) )

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
