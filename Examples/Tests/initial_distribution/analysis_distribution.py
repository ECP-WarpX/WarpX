#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the initial momentum distributions.
# 1 denotes gaussian distribution.
# 2 denotes maxwell-boltzmann distribution.
# The setup is a uniform plasma of electrons.
# The distribution is obtained through reduced diagnostic ParticleHistogram.

import numpy as np
import scipy.constants as scc

# load data
h1x = np.genfromtxt("h1x.txt")
h1y = np.genfromtxt("h1y.txt")
h1z = np.genfromtxt("h1z.txt")
h2x = np.genfromtxt("h2x.txt")
h2y = np.genfromtxt("h2y.txt")
h2z = np.genfromtxt("h2z.txt")

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
# note that parameters are chosen such that gaussian and
# maxwell-boltzmann distributions are identical
f1_error = 0.0
f2_error = 0.0
for i in range(bin_num):
    f1_error = f1_error + abs(f[i] - h1x[i+2])
    f1_error = f1_error + abs(f[i] - h1y[i+2])
    f1_error = f1_error + abs(f[i] - h1z[i+2])
    f2_error = f2_error + abs(f[i] - h2x[i+2])
    f2_error = f2_error + abs(f[i] - h2y[i+2])
    f2_error = f2_error + abs(f[i] - h2z[i+2])

f1_error = f1_error/bin_num
f2_error = f2_error/bin_num

print('Gaussian difference:', f1_error)
print('Maxwell-Boltzmann difference:', f2_error)
print('tolerance:', 0.008)

assert(f1_error < 0.008)
assert(f2_error < 0.008)
