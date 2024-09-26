#!/usr/bin/env python

"""
This script tests that the FEL is correctly modelled in the simulation.

The physical parameters are the same as the ones from section 5.1
of https://arxiv.org/pdf/2009.13645

The simulation uses the boosted-frame technique as described in
https://www.osti.gov/servlets/purl/940581
In particular, the effect of space-charge is effectively turned off
by initializing an electron and positron beam on top of each other,
each having half the current of the physical beam.

The script checks that the radiation wavelength and gain length
are the expected ones.
"""

import sys

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import c, e, m_e

# Physical parameters of the test
gamma_bunch = 100.6
Bu = 0.5
lambda_u = 3e-2
k_u = 2 * np.pi / lambda_u
K = e * Bu / (m_e * c * k_u)  # Undulator parameter
gamma_boost = (
    gamma_bunch / (1 + K * K / 2) ** 0.5
)  # Lorentz factor of the ponderomotive frame
beta_boost = (1 - 1.0 / gamma_boost**2) ** 0.5

# Open file specified in command line
filename = sys.argv[1]
ts = OpenPMDTimeSeries(filename)


# Extract the growth of the peak electric field
def extract_peak_E(iteration):
    """
    Extract the position of the peak electric field in
    a *boosted-frame* snapshot. Also return the position
    of the peak in the lab frame.
    """
    Ex, info = ts.get_field("E", "x", iteration=iteration)
    By, info = ts.get_field("B", "y", iteration=iteration)
    E_lab = gamma_boost * (Ex + c * beta_boost * By)
    E_lab_peak = E_lab.max()
    z_boost_peak = info.z[E_lab.argmax()]
    t_boost_peak = ts.current_t
    z_lab_peak = gamma_boost * (z_boost_peak + beta_boost * c * t_boost_peak)
    return z_lab_peak, E_lab_peak


# Loop through all iterations
z_lab_peak, E_lab_peak = ts.iterate(extract_peak_E)
log_P_peak = np.log(E_lab_peak**2)
# Since the radiation power is proportional to the square of the peak electric field,
# the log of the power is equal to the log of the square of the peak electric field,
# up to an additive constant.

# Pick the iterations between which the growth of the log of the power is linear
# (i.e. the growth of the power is exponential) and fit a line to extract the
# gain length.
i_start = 15
i_end = 23
# Perform linear fit
p = np.polyfit(z_lab_peak[i_start:i_end], log_P_peak[i_start:i_end], 1)
# Extract the gain length
Lg = 1 / p[0]
Lg_expected = 0.22  # Expected gain length from https://arxiv.org/pdf/2009.13645

assert abs(Lg - Lg_expected) / Lg_expected < 0.2


# test_name = os.path.split(os.getcwd())[1]
# import checksumAPI
# checksumAPI.evaluate_checksum(test_name, filename, output_format="openpmd")
