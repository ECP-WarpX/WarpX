#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe, Ilian Kara-Mostefa
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external lasy file.
#
# - Generate an input lasy file with a gaussian laser pulse.
# - Run the WarpX simulation for time T, when the pulse is fully injected
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation in 2D, for both envelope and central frequency

import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, epsilon_0
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

laser_energy = 1.0
E_max = np.sqrt( 2*(2/np.pi)**(3/2)*laser_energy / (epsilon_0*w0**2*c*tt) )

# Function for the envelope
def gauss_env(T, X, Y, Z):
    # Function to compute the theory for the envelope
    inv_tau2 = 1./tt/tt
    inv_w_2 = 1.0/(w0*w0)
    exp_arg = - (X*X)*inv_w_2 - (Y*Y)*inv_w_2- inv_tau2 / c/c * (Z-T*c)*(Z-T*c)
    return E_max * np.real(np.exp(exp_arg))

def do_analysis(fname, compname, steps):
    ds = yt.load(fname)
    dt = ds.current_time.to_value()/steps

    # Define 3D meshes
    x = np.linspace(
        ds.domain_left_edge[0],
        ds.domain_right_edge[0],
        ds.domain_dimensions[0]).v
    y = np.linspace(
        ds.domain_left_edge[1],
        ds.domain_right_edge[1],
        ds.domain_dimensions[1]).v
    z = np.linspace(
        ds.domain_left_edge[ds.dimensionality-1],
        ds.domain_right_edge[ds.dimensionality-1],
        ds.domain_dimensions[ds.dimensionality-1]).v
    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')

    # Compute the theory for envelope
    env_theory = gauss_env(+t_c-ds.current_time.to_value(), X,Y,Z)+gauss_env(-t_c+ds.current_time.to_value(), X,Y,Z)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser))
    extent = [ds.domain_left_edge[ds.dimensionality-1], ds.domain_right_edge[ds.dimensionality-1],
          ds.domain_left_edge[0], ds.domain_right_edge[0] ]
    env_theory_slice= env_theory[:,env_theory.shape[1]//2, :]

    # Plot results
    plt.figure(figsize=(8,6))
    plt.subplot(221)
    plt.title('PIC field')
    plt.imshow(F_laser, extent=extent)
    plt.colorbar()
    plt.subplot(222)
    plt.title('PIC envelope')
    plt.imshow(env, extent=extent)
    plt.colorbar()
    plt.subplot(223)
    plt.title('Theory envelope')
    plt.imshow(env_theory_slice, extent=extent)
    plt.colorbar()
    plt.subplot(224)
    plt.title('Difference')
    plt.imshow(env-env_theory_slice, extent=extent)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(compname, bbox_inches='tight')

    relative_error_env = np.sum(np.abs(env-env_theory_slice)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < relative_error_threshold)

    fft_F_laser = np.fft.fftn(F_laser)

    freq_x = np.fft.fftfreq(F_laser.shape[0],dt)
    freq_z = np.fft.fftfreq(F_laser.shape[1],dt)

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = np.sqrt((freq_x[pos_max[0]])**2 + (freq_z[pos_max[1]])**2)
    exp_freq = c/wavelength
    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < relative_error_threshold)



def launch_analysis(executable):
    os.system("./" + executable + " inputs.2d_test diag1.file_prefix=diags/plotfiles/plt")
    do_analysis("diags/plotfiles/plt000251/", "comp_unf.pdf", 251)


def main() :

    from lasy.laser import Laser
    from lasy.profiles import GaussianProfile

    # Create a laser using lasy
    pol = (1, 0)
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
