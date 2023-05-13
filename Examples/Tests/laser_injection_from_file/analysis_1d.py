#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external binary file.
#
# - Generate an input binary file with a gaussian laser pulse.
# - Run the WarpX simulation for time T, when the pulse is fully injected
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation, for both envelope and central frequency

import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

import yt ; yt.funcs.mylog.setLevel(50)

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

#Maximum acceptable error for this test
relative_error_threshold = 0.065

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#Parameters of the gaussian beam
wavelength = 1.*um
w0 = 12.*um
tt = 10.*fs
t_c = 20.*fs

E_max= 16282454014843.37

# Function for the envelope
def gauss_env(T,Z):
    # Function to compute the theory for the envelope
    inv_tau2 = 1./tt/tt
    exp_arg = - inv_tau2 / c/c * (Z-T*c)*(Z-T*c)
    return E_max * np.real(np.exp(exp_arg))

def do_analysis(fname, compname, steps):
    ds = yt.load(fname)
    dt = ds.current_time.to_value()/steps

    z = np.linspace(
        ds.domain_left_edge[0].v,
        ds.domain_right_edge[0].v,
        ds.domain_dimensions[0])

    # Compute the theory for envelope
    env_theory = gauss_env(+t_c-ds.current_time.to_value(), z)+gauss_env(-t_c+ds.current_time.to_value(),z)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser))

    # Plot results
    plt.figure(figsize=(8,8))
    plt.subplot(221)
    plt.title('PIC field')
    plt.plot(z, F_laser)
    plt.subplot(222)
    plt.title('PIC envelope')
    plt.plot(z, env)
    plt.subplot(223)
    plt.title('Theory envelope')
    plt.plot(z,env_theory)
    plt.subplot(224)
    plt.title('Difference')
    plt.plot(z,env-env_theory)

    plt.tight_layout()
    plt.savefig(compname, bbox_inches='tight')

    relative_error_env = np.sum(np.abs(env-env_theory)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < relative_error_threshold)

    fft_F_laser = np.fft.fftn(F_laser)

    freq_z = np.fft.fftfreq(F_laser.shape[0],dt)

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = freq_z[pos_max[0]]
    exp_freq = c/wavelength
    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < relative_error_threshold)



def launch_analysis(executable):
    os.system("./" + executable + " inputs.1d_test diag1.file_prefix=diags/plotfiles/plt")
    do_analysis("diags/plotfiles/plt000251/", "comp_unf.pdf", 251)


def main() :

    from lasy.laser import Laser
    from lasy.profiles import GaussianProfile

    # Create a laser using lasy
    pol = (1, 0)
    laser_energy = 1.0  # J
    profile = GaussianProfile(wavelength, pol, laser_energy, w0, tt, t_peak=0)
    dim = "xyt"
    lo = (-25e-6, -25e-6, -20e-15)
    hi = (+25e-6, +25e-6, +20e-15)
    npoints = (100, 100, 100)
    laser = Laser(dim, lo, hi, npoints, profile)
    laser.normalize(laser_energy, kind="energy")
    laser.write_to_file("gaussianlaser3d")
    executables = glob.glob("*.ex")
    if len(executables) == 1 :
        launch_analysis(executables[0])
    else :
        assert(False)

    # Do the checksum test
    filename_end = "diags/plotfiles/plt000251/"
    test_name = os.path.split(os.getcwd())[1]
    checksumAPI.evaluate_checksum(test_name, filename_end)
    print('Passed')

if __name__ == "__main__":
    main()
