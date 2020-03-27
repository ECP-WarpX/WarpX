#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external binary file for DIM=3.
#
# - Generate an input binary file with a gaussian laser pulse.
# - Run the WarpX simulation for time T, when the pulse is fully injected
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation, for both envelope and central frequency

import yt ; yt.funcs.mylog.setLevel(50)
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.signal import hilbert
import glob
import os

#Maximum acceptable error for this test
relative_error_threshold_envelope = 0.12
relative_error_threshold_frequency = 0.05

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#Parameters of the gaussian beam
wavelength = 1.*um
w0 = 2.0*um
tt = 10.*fs
x_c = 0.*um
y_c = 0.*um
t_c = 20.*fs
foc_dist = 10*um
E_max = 1e12

#Parameters of the txy grid
x_l = -5.0*um
x_r = 5.0*um
x_points = 160
y_l = -5.0*um
y_r = 5.0*um
y_points = 160
t_l = 0.0*fs
t_r = 40.0*fs
t_points = 200
tcoords = np.linspace(t_l, t_r, t_points)
xcoords = np.linspace(x_l, x_r, x_points)
ycoords = np.linspace(y_l, y_r, y_points)

def gauss(T,X,Y,opt):
    """Compute the electric field for a Gaussian laser pulse.
       This is used to write the binary input file.
    """

    k0 = 2.0*np.pi/wavelength
    inv_tau2 = 1./tt/tt
    osc_phase = k0*c*(T-t_c)

    diff_factor = 1.0 + 1.0j* foc_dist * 2/(k0*w0*w0)
    inv_w_2 = 1.0/(w0*w0*diff_factor)

    pre_fact = np.exp(1.0j * osc_phase)

    if opt == '3d':
        pre_fact = pre_fact/diff_factor
    else:
        pre_fact = pre_fact/np.sqrt(diff_factor)

    exp_arg = - (X*X + Y*Y)*inv_w_2 - inv_tau2 * (T-t_c)*(T-t_c)

    return np.real(pre_fact * np.exp(exp_arg))

# Function for the envelope
def gauss_env(T,XX,YY,ZZ):
    '''Function to compute the theory for the envelope
    '''

    inv_tau2 = 1./tt/tt
    inv_w_2 = 1.0/(w0*w0)
    exp_arg = - (YY*YY + ZZ*ZZ)*inv_w_2 - inv_tau2 / c/c * (XX-T*c)*(XX-T*c)
    return E_max * np.real(np.exp(exp_arg))

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

def create_gaussian_3d():
   T, X, Y = np.meshgrid(tcoords, xcoords, ycoords, indexing='ij')
   E_t = gauss(T,X,Y,'3d')
   write_file("gauss_3d.txye", xcoords, ycoords, tcoords, E_t)
   write_file_unf("gauss_3d_unf.txye", xcoords, ycoords, tcoords, E_t)


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
        ds.domain_left_edge[2],
        ds.domain_right_edge[2],
        ds.domain_dimensions[2]).v

    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')

    # Compute the theory for envelope
    env_theory = gauss_env(+t_c-ds.current_time.to_value(), X, Y, Z) + gauss_env(-t_c+ds.current_time.to_value(), X, Y ,Z)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser, axis=0))

    relative_error_env = np.sum(np.abs(env-env_theory)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < relative_error_threshold_envelope)



    fft_F_laser = np.fft.fftn(F_laser)

    freq = np.array([np.fft.fftfreq(s,dt) for s in F_laser.shape])

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = np.sqrt((freq[0][pos_max[0]])**2 +  (freq[1][pos_max[1]]**2) +  (freq[2][pos_max[2]]**2))
    exp_freq = c/wavelength

    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < relative_error_threshold_frequency)



def launch_analysis(executable):
    create_gaussian_3d()
    os.system("./" + executable + " inputs.3d_test_txye")
    do_analysis("diags/plotfiles/plt00150/", "comp_unf.pdf", 150)
    os.system("sed 's/gauss_3d_unf.txye/gauss_3d.txye/g' inputs.3d_test_txye > inputs.3d_test_txye_non_unf")
    os.system("./" + executable + " inputs.3d_test_txye_non_unf")
    do_analysis("diags/plotfiles/plt00150/", "comp_non_unf.pdf", 150)


def main() :
    executables = glob.glob("main3d*")
    if len(executables) == 1 :
        launch_analysis(executables[0])
    else :
        assert(False)
    print('Passed')

if __name__ == "__main__":
    main()
