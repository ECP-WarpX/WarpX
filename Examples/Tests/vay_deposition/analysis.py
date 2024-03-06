#!/usr/bin/env python3

# Copyright 2019-2022
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os
import sys

import numpy as np
from scipy.constants import epsilon_0
import yt

yt.funcs.mylog.setLevel(50)

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Plotfile data set
fn = sys.argv[1]
ds = yt.load(fn)

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E)
data = ds.covering_grid(
    level=0,
    left_edge=ds.domain_left_edge,
    dims=ds.domain_dimensions)
rho  = data[('boxlib','rho')].to_ndarray()
divE = data[('boxlib','divE')].to_ndarray()
error_rel = np.amax(np.abs(divE-rho/epsilon_0))/np.amax(np.abs(rho/epsilon_0))
tolerance = 1e-3
print("Error on charge conservation:")
print("error_rel = {}".format(error_rel))
print("tolerance = {}".format(tolerance))
assert( error_rel < tolerance )

# Checksum analysis
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
