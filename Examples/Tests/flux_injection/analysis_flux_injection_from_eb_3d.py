#!/usr/bin/env python3
#
# Copyright 2024 Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script tests the emission of particles from the embedded boundary.
(In this case, the embedded boundary is a sphere.) We check that the
embedded boundary emits the correct number of particles, and that the
the particle distributions are consistent with the expected distributions.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy.constants import c, m_e
from scipy.special import erf

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")

yt.funcs.mylog.setLevel(0)

# Open plotfile specified in command line
fn = sys.argv[1]
ds = yt.load(fn)
ad = ds.all_data()
t_max = ds.current_time.item()  # time of simulation

# Total number of electrons expected:
# Simulation parameters determine the total number of particles emitted (Ntot)
flux = 1.0  # in m^-2.s^-1, from the input script
R = 2.0  # in m, radius of the sphere
emission_surface = 4 * np.pi * R**2  # in m^2
Ntot = flux * emission_surface * t_max

# Parameters of the histogram
hist_bins = 50
hist_range = [-0.5, 0.5]


# Define function that histogram and check the data
def gaussian_dist(u, u_th):
    return 1.0 / ((2 * np.pi) ** 0.5 * u_th) * np.exp(-(u**2) / (2 * u_th**2))


def gaussian_flux_dist(u, u_th, u_m):
    normalization_factor = u_th**2 * np.exp(-(u_m**2) / (2 * u_th**2)) + (
        np.pi / 2
    ) ** 0.5 * u_m * u_th * (1 + erf(u_m / (2**0.5 * u_th)))
    result = (
        1.0
        / normalization_factor
        * np.where(u > 0, u * np.exp(-((u - u_m) ** 2) / (2 * u_th**2)), 0)
    )
    return result


def compare_gaussian(u, w, u_th, label=""):
    du = (hist_range[1] - hist_range[0]) / hist_bins
    w_hist, u_hist = np.histogram(u, bins=hist_bins, weights=w / du, range=hist_range)
    u_hist = 0.5 * (u_hist[1:] + u_hist[:-1])
    w_th = Ntot * gaussian_dist(u_hist, u_th)
    plt.plot(u_hist, w_hist, label=label + ": simulation")
    plt.plot(u_hist, w_th, "--", label=label + ": theory")
    assert np.allclose(w_hist, w_th, atol=0.07 * w_th.max())


def compare_gaussian_flux(u, w, u_th, u_m, label=""):
    du = (hist_range[1] - hist_range[0]) / hist_bins
    w_hist, u_hist = np.histogram(u, bins=hist_bins, weights=w / du, range=hist_range)
    u_hist = 0.5 * (u_hist[1:] + u_hist[:-1])
    w_th = Ntot * gaussian_flux_dist(u_hist, u_th, u_m)
    plt.plot(u_hist, w_hist, label=label + ": simulation")
    plt.plot(u_hist, w_th, "--", label=label + ": theory")
    assert np.allclose(w_hist, w_th, atol=0.05 * w_th.max())


# Load data and perform check

plt.figure()

x = ad["electron", "particle_position_x"].to_ndarray() / (m_e * c)
y = ad["electron", "particle_position_y"].to_ndarray() / (m_e * c)
z = ad["electron", "particle_position_z"].to_ndarray() / (m_e * c)
ux = ad["electron", "particle_momentum_x"].to_ndarray() / (m_e * c)
uy = ad["electron", "particle_momentum_y"].to_ndarray() / (m_e * c)
uz = ad["electron", "particle_momentum_z"].to_ndarray() / (m_e * c)
w = ad["electron", "particle_weight"].to_ndarray()

# Check that the total number of particles emitted is correct
Ntot_sim = np.sum(w)
print("Ntot_sim = ", Ntot_sim)
print("Ntot = ", Ntot)
assert np.isclose(Ntot_sim, Ntot, rtol=0.01)

# Check that none of the particles are inside the EB
assert np.all(x**2 + y**2 + z**2 > R**2)

# Check that the normal component of the velocity is consistent with the expected distribution
r = np.sqrt(x**2 + y**2 + z**2)
un = ux * (x / r) + uy * (y / r) + uz * (z / r)  # normal component
compare_gaussian_flux(un, w, u_th=0.1, u_m=0.07, label="u_n")

# compare_gaussian(ux, w, u_th=0.1, label="u_x")
# compare_gaussian(uz, w, u_th=0.1, label="u_z")

plt.tight_layout()
plt.savefig("Distribution.png")

# Verify checksum
test_name = os.path.split(os.getcwd())[1]
# if re.search("single_precision", fn):
#    checksumAPI.evaluate_checksum(test_name, fn, rtol=1.0e-3)
# else:
#    checksumAPI.evaluate_checksum(test_name, fn)
