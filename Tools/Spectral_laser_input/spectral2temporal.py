"""
Utility to convert spectral information to temporal profile for use in WarpX

Based on the corresponding functionality in FBPIC:
https://github.com/fbpic/fbpic/blob/ce3d530b0ffffc14b334b76850825dafc9d18942/fbpic/lpa_utils/laser/longitudinal_laser_profiles.py#L190
"""

import os

from matplotlib import pyplot as plt
import numpy as np
import openpmd_api as io
from scipy.constants import c, e, epsilon_0, m_e
from scipy.interpolate import interp1d


def get_central_frequency(wavelength, intensity):

    lambda0 = np.trapz(wavelength*intensity, wavelength) * \
                            1. / np.trapz(intensity, wavelength)
    return lambda0

def spectral_to_temporal(spectral_file_name,
                         Ninterp = 1000):
    """
    """
    data = np.loadtxt(spectral_file_name, delimiter=',',skiprows=1)
    wavelength, intensity, spectral_phase = data.T
    wavelength = 1e-9 * wavelength
    omega_axis = 2 * np.pi * c / wavelength

    lambda0 = get_central_frequency(wavelength, intensity)

    # Ninterp = 1000
    lambda_resolution = lambda0 / Ninterp  # spectral resolution
    dt = lambda_resolution / c
    time_window = lambda0 * lambda0 / c / lambda_resolution
    Nt = np.round(time_window/dt).astype(int)



    # Define the time array and its corresponding frequency array
    time_arr = -0.5*time_window + dt*np.arange(Nt)
    omega_arr = 2*np.pi * np.fft.fftfreq( Nt, dt )

    spectral_inten_fn = interp1d( omega_axis,
                                intensity*wavelength**2,
                                fill_value=0, bounds_error=False)
    spectral_phase_fn = interp1d( omega_axis, spectral_phase,
                                fill_value=0, bounds_error=False)

    spectral_Efield = np.sqrt( spectral_inten_fn(omega_arr) ) * \
                        np.exp( 1j*spectral_phase_fn(omega_arr) )

    # Ek = wavelength**2 * np.sqrt(intensity)*np.exp(1j*spectral_phase)
    Ex = np.fft.fftshift(np.fft.fft(spectral_Efield))

    return time_arr, Ex

def save_to_openpmd(time_arr,
                    Ex_arr,
                    temporal_file_name,
                    author_name = None,
                    author_email = None):
    """
    Save `E(t_i), t_i` arrays in an openPMD file

    Notes
    ---
    This function assumes a regularly spaced t_i array

    The file structure is E/Ex and the information in t array is
     stored in the grid attributes of the mesh record E
    """
    Ex_n_arr = Ex_arr / max(abs(Ex_arr))
    dt = time_arr[1] - time_arr[0]

    # write to openPMD format

    series = io.Series(temporal_file_name, io.Access.create)
    i = series.iterations[1]
    if author_name is not None or author_email is not None:
        author_str = ''
        if author_name is not None:
            author_str += author_name
        if author_email is not None:
            author_str += ' ' + author_email
        series.author = author_str


    E = i.meshes["E"]
    # [io.Mesh_Record_Component.SCALAR]
    E.geometry = io.Geometry.cartesian
    # Ez = E[io.Mesh_Record_Component.SCALAR]
    E.grid_spacing = [dt]
    print('dt=',dt)
    E.grid_global_offset = [time_arr[0]]
    E.axis_labels = ['z']
    E.data_order = "C"
    E.unit_dimension = {io.Unit_Dimension.I: 1.0,
                        io.Unit_Dimension.J: 2.0}
    E.grid_unit_SI = 1.0

    Ex = E["x"]
    Ex.unit_SI = 1.
    Ex.position = [0]

    Ex.reset_dataset(io.Dataset(Ex_arr.dtype, Ex_arr.shape))
    Ex.store_chunk(Ex_n_arr)
    series.flush()

    del series


if __name__ == '__main__':


    import argparse
    parser = argparse.ArgumentParser("Command line utility for converting laser spectral intensity to temporal profile")
    parser.add_argument('spectral_file', help='location of spectral data .csv file')
    parser.add_argument('--temporal_file', '-t', help='location of output openPMD file')
    parser.add_argument('--author', '-a', help="Author's name to save in openPMD metadata")
    parser.add_argument('--email', '-e', help="Author's e-mail to save in openPMD metadata")
    args = parser.parse_args()

    if args.temporal_file is not None:
        temporal_file = args.temporal_file
    else:
        temporal_file = os.path.dirname(args.spectral_file) + '/temporal_profile.h5'
    times, Exs = spectral_to_temporal(args.spectral_file)
    save_to_openpmd(times,
                    Exs,
                    temporal_file,
                    args.author,
                    args.email)
