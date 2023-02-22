#!/usr/bin/env python3

# Copyright 2019-2023
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# Initialize a deuterium-deuterium thermal plasma as in section 4.1 of [1]
# and check that the simulation reactivity matches the analytical value
# plotted in figure 3b of [1].
#
# Reference:
# [1] Higginson, D.P., Link, A. and Schmidt, A., 2019.
#     A pairwise nuclear fusion algorithm for weighted particle-in-cell
#     plasma simulations. Journal of Computational Physics, 388, pp.439-453.
#     DOI: https://doi.org/10.1016/j.jcp.2019.03.020


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import yt

yt.funcs.mylog.setLevel(0)

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Name of the plotfile
fn = sys.argv[1]

# Load data from reduced diagnostics (physical time and neutron weights)
time = np.loadtxt('./reduced_diags/particle_number.txt', usecols=1)
neutron = np.loadtxt('./reduced_diags/particle_number.txt', usecols=9)

# Compute reactivity [cm^3]/[s]
dY_dt = np.abs(neutron[-1]-neutron[0])/(time[-1]-time[0])
delta_ij = 1
nD = 1e26
V = (2e-3)**3
sigma = dY_dt*(1+delta_ij)/(nD**2)/V*(1e2)**3

# Compare checksums with benchmark
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
