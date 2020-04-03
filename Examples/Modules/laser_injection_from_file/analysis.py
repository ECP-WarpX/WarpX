#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external binary file for DIM=2 and DIM=3.
#
# - Generate an input binary file with a gaussian laser pulse.
# - Run the WarpX simulation for time T, when the pulse is fully injected
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation, for both envelope and central frequency

import yt ; yt.funcs.mylog.setLevel(50)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import glob
import os
import argparse

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#################### PARAMS FOR 3D TEST#########################################
class params_3d:
    pass

#Maximum acceptable error for this test
params_3d.relative_error_threshold_envelope = 0.12
params_3d.relative_error_threshold_frequency = 0.05

#Parameters of the gaussian beam
params_3d.wavelength = 1.*um
params_3d.w0 = 2.0*um
params_3d.tt = 10.*fs
params_3d.x_c = 0.*um
params_3d.y_c = 0.*um
params_3d.t_c = 20.*fs
params_3d.foc_dist = 10*um
params_3d.E_max = 1e12

#Parameters of the txy grid
params_3d.x_l = -5.0*um
params_3d.x_r = 5.0*um
params_3d.x_points = 160
params_3d.y_l = -5.0*um
params_3d.y_r = 5.0*um
params_3d.y_points = 160
params_3d.t_l = 0.0*fs
params_3d.t_r = 40.0*fs
params_3d.t_points = 200

################################################################################

#################### PARAMS FOR 2D TEST#########################################
class params_2d:
    pass

#Maximum acceptable error for this test
params_2d.relative_error_threshold_envelope = 0.65
params_2d.relative_error_threshold_frequency = 0.65

#Parameters of the gaussian beam
params_2d.wavelength = 1.*um
params_2d.w0 = 6.0*um
params_2d.tt = 10.*fs
params_2d.x_c = 0.*um
params_2d.t_c = 20.*fs
params_2d.foc_dist = 10*um
params_2d.E_max = 1e12
params_2d.rot_angle = -np.pi/4.0

#Parameters of the txy grid
params_2d.x_l = -12.0*um
params_2d.x_r = 12.0*um
params_2d.x_points = 480
params_2d.t_l = 0.0*fs
params_2d.t_r = 40.0*fs
params_2d.t_points = 400

################################################################################

def gauss(T,X,Y, params, opt):
    """Compute the electric field for a Gaussian laser pulse.
       This is used to write the binary input file.
    """

    k0 = 2.0*np.pi/params.wavelength
    inv_tau2 = 1./params.tt/params.tt
    osc_phase = k0*c*(T-params.t_c)

    diff_factor = 1.0 + 1.0j* params.foc_dist * 2/(k0*params.w0*params.w0)
    inv_w_2 = 1.0/(params.w0*params.w0*diff_factor)

    pre_fact = np.exp(1.0j * osc_phase)

    if opt == '3d':
        pre_fact = pre_fact/diff_factor
    else:
        pre_fact = pre_fact/np.sqrt(diff_factor)

    exp_arg = - (X*X + Y*Y)*inv_w_2 - inv_tau2 * (T-params.t_c)*(T-params.t_c)

    return np.real(pre_fact * np.exp(exp_arg))


def write_file(fname, x, y, t, E):
    """ For a given filename fname, space coordinates x and y, time coordinate t
    and field E, write a WarpX-compatible input binary file containing the
    profile of the laser pulse
    """

    with open(fname, 'wb') as file:
        flag_unif = 0
        file.write(flag_unif.to_bytes(1, byteorder='little'))
        file.write((len(t)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(x)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(y)).to_bytes(4, byteorder='little', signed=False))
        file.write(t.tobytes())
        file.write(x.tobytes())
        file.write(y.tobytes())
        file.write(E.tobytes())


def write_file_unf(fname, x, y, t, E):
    """ For a given filename fname, space coordinates x and y, time coordinate t
    and field E, write a WarpX-compatible input binary file containing the
    profile of the laser pulse. This function should be used in the case
    of a uniform spatio-temporal mesh
    """

    with open(fname, 'wb') as file:
        flag_unif = 1
        file.write(flag_unif.to_bytes(1, byteorder='little'))
        file.write((len(t)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(x)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(y)).to_bytes(4, byteorder='little', signed=False))
        file.write(t[0].tobytes())
        file.write(t[-1].tobytes())
        file.write(x[0].tobytes())
        file.write(x[-1].tobytes())
        if len(y) == 1 :
            file.write(y[0].tobytes())
        else :
            file.write(y[0].tobytes())
            file.write(y[-1].tobytes())
        file.write(E.tobytes())


# Function for the envelope
def gauss_env_2d(T,XX,ZZ,params):
    '''Function to compute the theory for the envelope
    '''
    X = np.cos(params.rot_angle)*XX + np.sin(params.rot_angle)*ZZ
    Z = -np.sin(params.rot_angle)*XX + np.cos(params.rot_angle)*ZZ

    inv_tau2 = 1./params.tt/params.tt
    inv_w_2 = 1.0/(params.w0*params.w0)
    exp_arg = - (X*X)*inv_w_2 - inv_tau2 / c/c * (Z-T*c)*(Z-T*c)
    return params.E_max * np.real(np.exp(exp_arg))

# Function for the envelope
def gauss_env_3d(T,XX,YY,ZZ,params):
    '''Function to compute the theory for the envelope
    '''

    inv_tau2 = 1./params.tt/params.tt
    inv_w_2 = 1.0/(params.w0*params.w0)
    exp_arg = - (YY*YY + ZZ*ZZ)*inv_w_2 - inv_tau2 / c/c * (XX-T*c)*(XX-T*c)
    return params.E_max * np.real(np.exp(exp_arg))


def create_gaussian_2d():
   tcoords = np.linspace(params_2d.t_l, params_2d.t_r, params_2d.t_points)
   xcoords = np.linspace(params_2d.x_l, params_2d.x_r, params_2d.x_points)

   T, X, Y = np.meshgrid(tcoords, xcoords, np.array([0.0]), indexing='ij')
   E_t = gauss(T,X,Y,params_2d,'2d')
   write_file("gauss_2d.txye", xcoords, np.array([0.0]), tcoords, E_t)
   write_file_unf("gauss_2d_unf.txye", xcoords, np.array([0.0]), tcoords, E_t)


def create_gaussian_3d():
   tcoords = np.linspace(params_3d.t_l, params_3d.t_r, params_3d.t_points)
   xcoords = np.linspace(params_3d.x_l, params_3d.x_r, params_3d.x_points)
   ycoords = np.linspace(params_3d.y_l, params_3d.y_r, params_3d.y_points)

   T, X, Y = np.meshgrid(tcoords, xcoords, ycoords, indexing='ij')
   E_t = gauss(T,X,Y,params_3d,'3d')
   write_file("gauss_3d.txye", xcoords, ycoords, tcoords, E_t)
   write_file_unf("gauss_3d_unf.txye", xcoords, ycoords, tcoords, E_t)


def do_analysis_2d(fname, compname, steps):
    ds = yt.load(fname)

    dt = ds.current_time.to_value()/steps

    # Define 2D meshes
    x = np.linspace(
        ds.domain_left_edge[0],
        ds.domain_right_edge[0],
        ds.domain_dimensions[0]).v
    z = np.linspace(
        ds.domain_left_edge[ds.dimensionality-1],
        ds.domain_right_edge[ds.dimensionality-1],
        ds.domain_dimensions[ds.dimensionality-1]).v
    X, Z = np.meshgrid(x, z, sparse=False, indexing='ij')

    # Compute the theory for envelope
    env_theory = gauss_env_2d(+params_2d.t_c-ds.current_time.to_value(), X,Z,params_2d)+gauss_env_2d(-params_2d.t_c+ds.current_time.to_value(), X,Z,params_2d)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser))
    extent = [ds.domain_left_edge[ds.dimensionality-1], ds.domain_right_edge[ds.dimensionality-1],
          ds.domain_left_edge[0], ds.domain_right_edge[0] ]

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
    plt.imshow(env_theory, extent=extent)
    plt.colorbar()
    plt.subplot(224)
    plt.title('Difference')
    plt.imshow(env-env_theory, extent=extent)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(compname, bbox_inches='tight')

    relative_error_env = np.sum(np.abs(env-env_theory)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < params_2d.relative_error_threshold_envelope)

    fft_F_laser = np.fft.fft2(F_laser)

    freq_rows = np.fft.fftfreq(F_laser.shape[0],dt)
    freq_cols = np.fft.fftfreq(F_laser.shape[1],dt)

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = np.sqrt((freq_rows[pos_max[0]])**2 +  (freq_cols[pos_max[1]]**2))
    exp_freq = c/params_2d.wavelength

    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < params_2d.relative_error_threshold_frequency)


def do_analysis_3d(fname, compname, steps):
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
        ds.domain_left_edge[2],
        ds.domain_right_edge[2],
        ds.domain_dimensions[2]).v

    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')

    # Compute the theory for envelope
    env_theory = gauss_env_3d(+params_3d.t_c-ds.current_time.to_value(), X, Y, Z,params_3d) + gauss_env_3d(-params_3d.t_c+ds.current_time.to_value(), X, Y,Z, params_3d)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser, axis=0))

    relative_error_env = np.sum(np.abs(env-env_theory)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < params_3d.relative_error_threshold_envelope)

    fft_F_laser = np.fft.fftn(F_laser)

    freq = np.array([np.fft.fftfreq(s,dt) for s in F_laser.shape])

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = np.sqrt((freq[0][pos_max[0]])**2 +  (freq[1][pos_max[1]]**2) +  (freq[2][pos_max[2]]**2))
    exp_freq = c/params_3d.wavelength

    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < params_3d.relative_error_threshold_frequency)


def launch_analysis_3d(executable):
    create_gaussian_3d()
    os.system("./" + executable + " inputs.3d_test_txye")
    do_analysis_3d("diags/plotfiles/plt00150/", "comp_unf.pdf", 150)
    os.system("sed 's/gauss_3d_unf.txye/gauss_3d.txye/g' inputs.3d_test_txye > inputs.3d_test_txye_non_unf")
    os.system("./" + executable + " inputs.3d_test_txye_non_unf")
    do_analysis_3d("diags/plotfiles/plt00150/", "comp_non_unf.pdf", 150)


def launch_analysis_2d(executable):
    create_gaussian_2d()
    os.system("./" + executable + " inputs.2d_test_txye")
    do_analysis_2d("diags/plotfiles/plt00250/", "comp_unf.pdf", 250)
    os.system("sed 's/gauss_2d_unf.txye/gauss_2d.txye/g' inputs.2d_test_txye > inputs.2d_test_txye_non_unf")
    os.system("./" + executable + " inputs.2d_test_txye_non_unf")
    do_analysis_2d("diags/plotfiles/plt00250/", "comp_non_unf.pdf", 250)


def do_2d():
    executables = glob.glob("main2d*")
    if len(executables) == 1 :
        launch_analysis_2d(executables[0])
    else :
        assert(False)

def do_3d():
    executables = glob.glob("main3d*")
    if len(executables) == 1 :
        launch_analysis_3d(executables[0])
    else :
        assert(False)


def main() :

    parser = argparse.ArgumentParser(description='Performs tests for laser injection from file in 2D and 3D')
    parser.add_argument('dimensionality', type=str, help='Dimensionality of the test. Possible options: 2D, 3D.')

    args = parser.parse_args()

    if args.dimensionality == "2D":
        do_2d()
    elif args.dimensionality == "3D":
        do_3d()
    else:
        assert(False)

    print('Passed')


if __name__ == "__main__":
    main()
