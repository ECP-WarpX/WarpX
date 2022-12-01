#!/usr/bin/env python3
#
# Copyright 2021 Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script tests the Gaussian-flux injection (and in particular
the rejection method in WarpX that we use to generate the right
velocity distribution).

Two population of particles are injected, with a slightly different
ratio of u_m/u_th. (This is in order to test the two different
rejection methods implemented in WarpX, which depend on the u_m/u_th ratio.)

After the particles are emitted with flux injection, this script produces
histograms of the velocity distribution and compares it with the expected
velocity distibution (Gaussian or Gaussian-flux depending on the direction
of space)
"""
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, m_e, m_p
from scipy.special import erf
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

yt.funcs.mylog.setLevel(0)

# Open plotfile specified in command line
fn = sys.argv[1]
ds = yt.load( fn )
ad = ds.all_data()
t_max = ds.current_time.item()  # time of simulation

# Total number of electrons expected:
# Simulation parameters determine the total number of particles emitted (Ntot)
flux = 1. # in m^-2.s^-1, from the input script
emission_surface = 8*8 # in m^2
Ntot = flux * emission_surface * t_max

# Parameters of the histogram
hist_bins = 50
hist_range = [-0.5, 0.5]

# Define function that histogram and check the data

def gaussian_dist(u, u_th):
    return 1./((2*np.pi)**.5*u_th) * np.exp(-u**2/(2*u_th**2) )

def gaussian_flux_dist(u, u_th, u_m):
    normalization_factor = u_th**2 * np.exp(-u_m**2/(2*u_th**2)) + (np.pi/2)**.5*u_m*u_th * (1 + erf(u_m/(2**.5*u_th)))
    return 1./normalization_factor * np.where( u>0, u * np.exp(-(u-u_m)**2/(2*u_th**2)), 0 )

def compare_gaussian(u, w, u_th, label=''):
    du = (hist_range[1]-hist_range[0])/hist_bins
    w_hist, u_hist = np.histogram(u, bins=hist_bins, weights=w/du, range=hist_range)
    u_hist = 0.5*(u_hist[1:]+u_hist[:-1])
    w_th = Ntot*gaussian_dist(u_hist, u_th)
    plt.plot( u_hist, w_hist, label=label+': simulation' )
    plt.plot( u_hist, w_th, '--', label=label+': theory' )
    assert np.allclose( w_hist, w_th, atol=0.07*w_th.max() )

def compare_gaussian_flux(u, w, u_th, u_m, label=''):
    du = (hist_range[1]-hist_range[0])/hist_bins
    w_hist, u_hist = np.histogram(u, bins=hist_bins, weights=w/du, range=hist_range)
    u_hist = 0.5*(u_hist[1:]+u_hist[:-1])
    w_th = Ntot*gaussian_flux_dist(u_hist, u_th, u_m)
    plt.plot( u_hist, w_hist, label=label+': simulation' )
    plt.plot( u_hist, w_th, '--', label=label+': theory' )
    assert np.allclose( w_hist, w_th, atol=0.05*w_th.max() )

# Load data and perform check

plt.figure(figsize=(5,7))

plt.subplot(211)
plt.title('Electrons')

ux = ad['electron','particle_momentum_x'].to_ndarray()/(m_e*c)
uy = ad['electron','particle_momentum_y'].to_ndarray()/(m_e*c)
uz = ad['electron','particle_momentum_z'].to_ndarray()/(m_e*c)
w = ad['electron', 'particle_weight'].to_ndarray()

compare_gaussian(ux, w, u_th=0.1, label='u_x')
compare_gaussian_flux(uy, w, u_th=0.1, u_m=0.07, label='u_y')
compare_gaussian(uz, w, u_th=0.1, label='u_z')
plt.legend(loc=0)

plt.subplot(212)
plt.title('Protons')

ux = ad['proton','particle_momentum_x'].to_ndarray()/(m_p*c)
uy = ad['proton','particle_momentum_y'].to_ndarray()/(m_p*c)
uz = ad['proton','particle_momentum_z'].to_ndarray()/(m_p*c)
w = ad['proton', 'particle_weight'].to_ndarray()

compare_gaussian_flux(ux, w, u_th=0.1, u_m=0.05, label='u_x')
compare_gaussian(uy, w, u_th=0.1, label='u_y')
compare_gaussian(uz, w, u_th=0.1, label='u_z')
plt.legend(loc=0)

plt.savefig('Distribution.png')

# Verify checksum
test_name = os.path.split(os.getcwd())[1]
if re.search( 'single_precision', fn ):
    checksumAPI.evaluate_checksum(test_name, fn, rtol=1.e-3)
else:
    checksumAPI.evaluate_checksum(test_name, fn)
