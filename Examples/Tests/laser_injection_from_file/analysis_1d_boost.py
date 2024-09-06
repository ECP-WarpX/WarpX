#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe, Ilian Kara-Mostefa
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external lasy file:
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation in 1D, for both envelope and central frequency

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy.constants import c, epsilon_0
from scipy.signal import hilbert

yt.funcs.mylog.setLevel(50)

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

# Maximum acceptable error for this test
relative_error_threshold = 0.065

# Physical parameters
um = 1.0e-6
fs = 1.0e-15

# Parameters of the gaussian beam
wavelength = 1.0 * um
w0 = 12.0 * um
tt = 10.0 * fs
t_c = 20.0 * fs

laser_energy = 1.0
E_max = np.sqrt(
    2 * (2 / np.pi) ** (3 / 2) * laser_energy / (epsilon_0 * w0**2 * c * tt)
)


# Function for the envelope
def gauss_env(T, Z):
    # Function to compute the theory for the envelope
    inv_tau2 = 1.0 / tt / tt
    exp_arg = -inv_tau2 / c / c * (Z - T * c) * (Z - T * c)
    return E_max * np.real(np.exp(exp_arg))


filename = sys.argv[1]
compname = "comp_unf.pdf"
ds = yt.load(filename)
dz = (ds.domain_right_edge[0].v - ds.domain_left_edge[0].v) / ds.domain_dimensions[0]
dt = dz / c

z = np.linspace(
    ds.domain_left_edge[0].v, ds.domain_right_edge[0].v, ds.domain_dimensions[0]
)

# Compute the theory for envelope
env_theory = gauss_env(+t_c - ds.current_time.to_value(), z) + gauss_env(
    -t_c + ds.current_time.to_value(), z
)

# Read laser field in PIC simulation, and compute envelope
all_data_level_0 = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
F_laser = all_data_level_0["boxlib", "Ey"].v.squeeze()
env = abs(hilbert(F_laser))

# Plot results
plt.figure(figsize=(8, 8))
plt.subplot(221)
plt.title("PIC field")
plt.plot(z, F_laser)
plt.subplot(222)
plt.title("PIC envelope")
plt.plot(z, env)
plt.subplot(223)
plt.title("Theory envelope")
plt.plot(z, env_theory)
plt.subplot(224)
plt.title("Difference")
plt.plot(z, env - env_theory)

plt.tight_layout()
plt.savefig(compname, bbox_inches="tight")

relative_error_env = np.sum(np.abs(env - env_theory)) / np.sum(np.abs(env))
print("Relative error envelope: ", relative_error_env)
assert relative_error_env < relative_error_threshold

fft_F_laser = np.fft.fftn(F_laser)

freq_z = np.fft.fftfreq(F_laser.shape[0], dt)

pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

freq = np.abs(freq_z[pos_max[0]])
exp_freq = c / wavelength
relative_error_freq = np.abs(freq - exp_freq) / exp_freq
print("Relative error frequency: ", relative_error_freq)
assert relative_error_freq < relative_error_threshold

# Do the checksum test
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
