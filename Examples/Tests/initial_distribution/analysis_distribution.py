#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests initial distributions.
# 1 denotes gaussian distribution.
# 2 denotes maxwell-boltzmann distribution.
# 3 denotes maxwell-juttner distribution.
# 4 denotes gaussian position distribution.
# The distribution is obtained through reduced diagnostic ParticleHistogram.

import numpy as np
import scipy.constants as scc
import scipy.special as scs

# print tolerance
tolerance = 0.009
print('tolerance:', tolerance)

#===============================
# gaussian and maxwell-boltzmann
#===============================

# load data
h1x = np.genfromtxt("h1x.txt")
h1y = np.genfromtxt("h1y.txt")
h1z = np.genfromtxt("h1z.txt")
h2x = np.genfromtxt("h2x.txt")
h2y = np.genfromtxt("h2y.txt")
h2z = np.genfromtxt("h2z.txt")

# parameters of bin
bin_min = -4.0e-2 * scc.c
bin_max = +4.0e-2 * scc.c
bin_num = 50
bin_size = (bin_max-bin_min)/bin_num

# parameters of theory
u_rms = 0.01
gamma = np.sqrt(1+u_rms*u_rms)
v_rms = u_rms / gamma * scc.c

# compute the analytical solution
f = np.zeros( bin_num )
for i in range(bin_num):
    v = (i+0.5)*bin_size + bin_min
    f[i] = np.exp(-0.5*(v/v_rms)**2)/(v_rms*np.sqrt(2.0*scc.pi))

# normalization
f = f / np.amax( f )

# compute error
# note that parameters are chosen such that gaussian and
# maxwell-boltzmann distributions are identical
f1_error = np.sum(np.abs(f-h1x[2:])+np.abs(f-h1y[2:])+np.abs(f-h1z[2:]))/bin_num
f2_error = np.sum(np.abs(f-h2x[2:])+np.abs(f-h2y[2:])+np.abs(f-h2z[2:]))/bin_num

print('Gaussian distribution difference:', f1_error)
print('Maxwell-Boltzmann distribution difference:', f2_error)

assert(f1_error < tolerance)
assert(f2_error < tolerance)

#================
# maxwell-juttner
#================

# load data
h3  = np.genfromtxt("h3.txt")

# parameters of bin
bin_min = 1.0
bin_max = 12.0
bin_num = 50
bin_size = (bin_max-bin_min)/bin_num
theta   = 1.0
K2      = scs.kn(2,1.0/theta)

# compute the analytical solution
f = np.zeros( bin_num )
for i in range(bin_num):
    gamma = (i+0.5)*bin_size + bin_min
    beta  = np.sqrt(1.0-1.0/gamma**2)
    f[i]  = gamma**2 * beta / (theta*K2) * np.exp(-gamma/theta)

# normalization
f = f / np.amax( f )

# compute error
f3_error = np.sum( np.abs(f-h3[2:]) ) / bin_num

print('Maxwell-Juttner distribution difference:', f3_error)

assert(f3_error < tolerance)

#==============
# gaussian beam
#==============

# load data
h4x = np.genfromtxt("h4x.txt")
h4y = np.genfromtxt("h4y.txt")
h4z = np.genfromtxt("h4z.txt")

# parameters of bin
bin_min = -1.0
bin_max = +1.0
bin_num = 50
bin_size = (bin_max-bin_min)/bin_num

# parameters of theory
x_rms = 0.25

# compute the analytical solution
f = np.zeros( bin_num )
for i in range(bin_num):
    x = (i+0.5)*bin_size + bin_min
    f[i] = np.exp(-0.5*(x/x_rms)**2)/(x_rms*np.sqrt(2.0*scc.pi))

# normalization
f = f / np.amax( f )

# compute error
f4_error = np.sum(np.abs(f-h4x[2:])+np.abs(f-h4y[2:])+np.abs(f-h4z[2:]))/bin_num

print('Gaussian position distribution difference:', f4_error)

assert(f4_error < tolerance)
