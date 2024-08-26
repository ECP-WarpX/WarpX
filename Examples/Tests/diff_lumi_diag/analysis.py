#!/usr/bin/env python3
# This test checks the differential luminosity of the beam-beam interaction
# in the case of two Gaussian beams crossing rigidly at the interaction point.
# The beams have a Gaussian distribution both in energy and in transverse positions.
# In that case, the differential luminosity can be calculated analytically.

import os
import sys

import numpy as np
from read_raw_data import read_reduced_diags_histogram
from scipy.constants import eV

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

# Extract the differential luminosity from the file
_, _, E_bin, bin_data = read_reduced_diags_histogram(
    "./diags/reducedfiles/DifferentialLuminosity_beam1_beam2.txt"
)
dL_dE_sim = bin_data[-1]  # Differential luminosity at the end of the simulation

# Beam parameters
N = 1.2e10
E_beam = 125e9 * eV
sigma_x = 500e-9
sigma_y = 10e-9

# Compute the analytical differential luminosity for 2 Gaussian beams
sigma_E1 = 0.02 * E_beam
sigma_E2 = 0.03 * E_beam
sigma_E = np.sqrt(
    sigma_E1**2 + sigma_E2**2
)  # Expected width of the differential luminosity
dL_dE_th = (
    N**2
    / (2 * (2 * np.pi) ** 1.5 * sigma_x * sigma_y * sigma_E)
    * np.exp(-((E_bin - 2 * E_beam) ** 2) / (2 * sigma_E**2))
)

# Check that the simulation result and analytical result match
error = abs(dL_dE_sim - dL_dE_th).max() / abs(dL_dE_th).max()
tol = 1e-2
print('Relative error: ', error)
print('Tolerance: ', tol )
assert error < tol

# Get name of the test
fn = sys.argv[1]
test_name = os.path.split(os.getcwd())[1]

# Run checksum regression test
checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-2)
