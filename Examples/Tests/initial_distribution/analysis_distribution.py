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
# 5 denotes maxwell-juttner distribution w/ spatially varying temperature
# The distribution is obtained through reduced diagnostic ParticleHistogram.

import numpy as np
import scipy.constants as scc
import scipy.special as scs
from read_raw_data import read_reduced_diags_histogram, read_reduced_diags
import sys
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# print tolerance
tolerance = 0.02
print('Tolerance:', tolerance)

#===============================
# gaussian and maxwell-boltzmann
#===============================

# load data
bin_value, h1x = read_reduced_diags_histogram("h1x.txt")[2:]
h1y = read_reduced_diags_histogram("h1y.txt")[3]
h1z = read_reduced_diags_histogram("h1z.txt")[3]
h2x = read_reduced_diags_histogram("h2x.txt")[3]
h2y = read_reduced_diags_histogram("h2y.txt")[3]
h2z = read_reduced_diags_histogram("h2z.txt")[3]

# parameters of theory
u_rms = 0.01
gamma = np.sqrt(1.0+u_rms*u_rms)
v_rms = u_rms / gamma * scc.c
n     = 1.0e21
V     = 8.0
db    = 0.0016

# compute the analytical solution
f = n*V*scc.c*db*np.exp(-0.5*(bin_value*scc.c/v_rms)**2)/(v_rms*np.sqrt(2.0*scc.pi))
f_peak = np.amax(f)

# compute error
# note that parameters are chosen such that gaussian and
# maxwell-boltzmann distributions are identical
f1_error = np.sum(np.abs(f-h1x)+np.abs(f-h1y)+np.abs(f-h1z))/bin_value.size / f_peak
f2_error = np.sum(np.abs(f-h2x)+np.abs(f-h2y)+np.abs(f-h2z))/bin_value.size / f_peak

print('Gaussian distribution difference:', f1_error)
print('Maxwell-Boltzmann distribution difference:', f2_error)

assert(f1_error < tolerance)
assert(f2_error < tolerance)

#================
# maxwell-juttner
#================

# load data
bin_value, bin_data = read_reduced_diags_histogram("h3.txt")[2:]
bin_data_filtered = read_reduced_diags_histogram("h3_filtered.txt")[3]

# parameters of theory
theta = 1.0
K2    = scs.kn(2,1.0/theta)
n     = 1.0e21
V     = 8.0
db    = 0.22

# compute the analytical solution
f = n*V*db * bin_value**2 * np.sqrt(1.0-1.0/bin_value**2) / \
    (theta*K2) * np.exp(-bin_value/theta)
f_peak = np.amax(f)

# analytical solution for the filtered histogram: we just filter out gamma values < 5.5
f_filtered = f*(bin_value > 5.5)

# compute error
f3_error = np.sum( np.abs(f-bin_data) + np.abs(f_filtered-bin_data_filtered) ) \
           / bin_value.size / f_peak

print('Maxwell-Juttner distribution difference:', f3_error)

assert(f3_error < tolerance)

#==============
# gaussian beam
#==============

# load data
bin_value, h4x = read_reduced_diags_histogram("h4x.txt")[2:]
h4y = read_reduced_diags_histogram("h4y.txt")[3]
h4z = read_reduced_diags_histogram("h4z.txt")[3]
_, bmmntr = read_reduced_diags("bmmntr.txt")
charge = bmmntr['charge'][0]

# parameters of theory
x_rms = 0.25
z_cut = 2.
q_tot = -1.0e-20
q_e   = -1.602176634e-19
npart = q_tot/q_e
db    = bin_value[1]-bin_value[0]

# compute the analytical solution
f_xy = npart*db * np.exp(-0.5*(bin_value/x_rms)**2)/(x_rms*np.sqrt(2.0*scc.pi)) * scs.erf(z_cut/np.sqrt(2.))
f_z = npart*db * np.exp(-0.5*(bin_value/x_rms)**2)/(x_rms*np.sqrt(2.0*scc.pi))
f_z[ np.absolute(bin_value) > z_cut * x_rms ] = 0.
f_peak = np.amax(f_z)
q_tot_cut = q_tot * scs.erf(z_cut/np.sqrt(2.))

# compute error
f4_error = np.sum(np.abs(f_xy-h4x)+np.abs(f_xy-h4y)+np.abs(f_z-h4z))/bin_value.size / f_peak
charge_error = np.abs((q_tot_cut - charge) / q_tot)

do_plot = False
if do_plot:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.subplot(121)
    plt.plot(bin_value, f_xy, '+-', label='ref')
    plt.plot(bin_value, h4x, '+--', label='sim')
    plt.legend()
    plt.subplot(122)
    plt.plot(bin_value, f_z, '+-', label='ref')
    plt.plot(bin_value, h4z, '+--', label='sim')
    plt.legend()
    plt.savefig('toto.pdf', bbox_inches='tight')

print('Gaussian position distribution difference:', f4_error)
assert(f4_error < tolerance)

print('Relative beam charge difference:', charge_error)
assert(charge_error < tolerance)

#=============================================
# maxwell-juttner with temperature from parser
#=============================================

# load data
bin_value, bin_data_neg = read_reduced_diags_histogram("h5_neg.txt")[2:]
bin_data_pos = read_reduced_diags_histogram("h5_pos.txt")[3]

# parameters of theory
# _neg denotes where x<0, _pos where x>0
theta_neg = 1.0
theta_pos = 2.0
K2_neg    = scs.kn(2,1.0/theta_neg)
K2_pos    = scs.kn(2,1.0/theta_pos)
n     = 1.0e21
V     = 8.0 / 2 # because each of these are for half the domain
db    = 0.22

# compute the analytical solution for each half of the domain
f_neg = n*V*db * bin_value**2 * np.sqrt(1.0-1.0/bin_value**2) / \
    (theta_neg*K2_neg) * np.exp(-bin_value/theta_neg)
f_neg_peak = np.amax(f_neg)
f_pos = n*V*db * bin_value**2 * np.sqrt(1.0-1.0/bin_value**2) / \
    (theta_pos*K2_pos) * np.exp(-bin_value/theta_pos)
f_pos_peak = np.amax(f_pos)
f_peak = max(f_neg_peak, f_pos_peak)

# compute error
f5_error = np.sum( np.abs(f_neg-bin_data_neg) + np.abs(f_pos-bin_data_pos) ) \
           / bin_value.size / f_peak

print('Maxwell-Juttner parser temperature difference:', f5_error)

assert(f5_error < tolerance)

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
