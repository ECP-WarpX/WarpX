#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the initial Gaussian momentum distribution.
# The setup is a uniform plasma of electrons.
# The distribution is obtained through reduced diagnostic ParticleHistogram.

import numpy as np
import scipy.constants as scc

# get results from reduced diagnostics

# load data
hx = np.genfromtxt("hx.txt")
hy = np.genfromtxt("hy.txt")
hz = np.genfromtxt("hz.txt")

# parameters of bin
bin_min = -1.2e7
bin_max = +1.2e7
bin_num = 50
bin_size = (bin_max-bin_min)/bin_num

# parameters of theory
u_rms = 0.01
gamma = np.sqrt(1+u_rms*u_rms)
v_rms = u_rms / gamma * scc.c

# compute the analytical solution from the theory
f = np.zeros( bin_num )
for i in range(bin_num):
    v = (i+0.5)*bin_size + bin_min
    f[i] = np.exp(-0.5*(v/v_rms)**2)/(v_rms*np.sqrt(2.0*scc.pi))

# normalization
f_max = np.amax(f)
for i in range(bin_num):
    f[i] = f[i] / f_max

# compute error
f_error = 0
for i in range(bin_num):
    f_error = f_error + abs(f[i] - hx[i+2])
    f_error = f_error + abs(f[i] - hy[i+2])
    f_error = f_error + abs(f[i] - hz[i+2])
f_error = f_error/bin_num

print('difference:', f_error)
print('tolerance:', 0.008)

assert(f_error < 0.008)
