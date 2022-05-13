"""
Python script to plot the evolution of the total electromagnetic field
energy over time (using data from WarpX's reduced diagnostics) and the
growth rate computed from the electric field component Ez (using data
from WarpX's standard diagnostics).

Users can set the input parameters defined in the main environment below,
providing the correct paths to the diagnostics data to be post-processed
as well as the start and end time steps between which the growth rate is
computed.
"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import yt

yt.funcs.mylog.setLevel(50)
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# Constants
c = scipy.constants.speed_of_light

def set_axes(ax, iter_min, iter_max):
    """
    Set axes properties.

    Parameters
    ----------
    ax : Axes
        Axes object.
    iter_min : int
        Start time step to print in title.
    iter_max : int
        End time step to print in title.
    """
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    title = 'Growth Rate (time step {} to {})'.format(iter_min, iter_max)
    ax.set_title(title)

def set_colorbar(cb):
    """
    Set colorbar properties.

    Parameters
    ----------
    cb : Colorbar
        Colorbar object.
    """
    cb.ax.tick_params()
    cb.ax.yaxis.get_offset_text().set()
    cb.ax.yaxis.offsetText.set_ha('center')
    cb.ax.yaxis.offsetText.set_va('bottom')
    cb.formatter.set_powerlimits((0,0)) # comment out if using log scale
    cb.update_ticks()

def load_data(ds, field):
    """
    Load plotfile data.

    Parameters
    ----------
    ds : yt.frontends.boxlib.data_structures.WarpXDataset
        Yt dataset extracted from plotfile.

    field : string tuple
        Tuple of strings to specify the grid field to be loaded
        (e.g., ('boxlib', 'Ez')).

    Returns
    -------
    data : numpy.ndarray
        Array of data.
    """
    level = 0
    grids = [grid for grid in ds.index.grids if grid.Level == level]
    LeftEdge = np.min(np.array([grid.LeftEdge.v for grid in grids]), axis=0)
    all_data = ds.covering_grid(level=level, left_edge=LeftEdge,
                                dims=ds.domain_dimensions)
    data = all_data[field].v.squeeze()
    return data

if __name__ == '__main__':

    #--------------------------------------------------------------------------
    # Users can set the input parameters defined in the main environment here,
    # providing the correct paths to the diagnostics data to be post-processed
    # as well as the start and end time steps between which the growth rate is
    # computed.
    #--------------------------------------------------------------------------

    # Path to standard diagnostics,
    # path to reduced diagnostics,
    # path to save output figures
    path_diags         = './path_to_standard_diagnostics/'
    path_diags_reduced = './path_to_reduced_diagnostics/'
    path_save          = './path_to_save_output_figures/'

    # Name of standard diagnostics
    # and reduced diagnostics (for total electromagnetic field energy)
    diag_name = 'diag'
    diag_name_reduced = 'field_energy.txt'

    # Start and end iteration to compute growth rate
    iter_min = 500
    iter_max = 1000

    #--------------------------------------------------------------------------

    # Get list of output files
    fn_l = glob.glob(os.path.join(path_diags, diag_name + '*[0-9]'))
    fn_l.sort()

    # Measure diagnostics frequency from two plotfiles and
    # add 1 if last iteration does not match diagnostics frequency
    dfreq = (int((re.search(os.path.join(path_diags, diag_name + '(.*)'), fn_l[1])).group(1))
            - int((re.search(os.path.join(path_diags, diag_name + '(.*)'), fn_l[0])).group(1)))
    itmin = iter_min // dfreq
    itmax = iter_max // dfreq

    # Get list of data at itmin and itmax
    # (type: yt.frontends.boxlib.data_structures.WarpXDataset)
    ds_itmin = yt.load(fn_l[itmin])
    ds_itmax = yt.load(fn_l[itmax])

    # Get data
    field = ('boxlib', 'Ez')
    data_itmin = load_data(ds_itmin, field)
    data_itmax = load_data(ds_itmax, field)

    # Get parameters from first data set
    # Dt is the time between the min and max iteration, not the time step dt
    Dt = (ds_itmax.current_time - ds_itmin.current_time).to_ndarray().item()
    dxz = ds_itmin.domain_width / ds_itmin.domain_dimensions
    dx = dxz[0]
    dz = dxz[1]
    kxmax = np.pi / dx
    kzmax = np.pi / dz
    Nx = data_itmin.shape[0]
    Nz = data_itmin.shape[1]

    # Compute spectrum
    spectrum_itmin = np.fft.fftshift(np.fft.fft2(data_itmin))[Nx//2:,Nz//2:]
    spectrum_itmax = np.fft.fftshift(np.fft.fft2(data_itmax))[Nx//2:,Nz//2:]

    # Compute growth rate
    growth_rate = np.log(np.abs(spectrum_itmax)) - np.log(np.abs(spectrum_itmin))
    growth_rate /= c*Dt

    # Plot energy over time (based on reduced diagnostics, if available)
    energy = np.loadtxt(os.path.join(path_diags_reduced, diag_name_reduced), usecols=2)
    fig, ax = plt.subplots(dpi=100)
    ax.plot(energy, '-')
    ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel(r'Time Steps')
    ax.set_title('Total Field Energy')
    fig.tight_layout()
    name = os.path.join(path_save, 'total_energy')
    fig.savefig(name + '.png', dpi=100)
    fig.show()

    # Plot growth rate
    fig, ax = plt.subplots(dpi=100)
    cax = make_axes_locatable(ax).append_axes('right', size='5%', pad='5%')
    # Min and max values for the colorbar of the plot
    # (possibly change these here to get better plots)
    vmax = np.max(growth_rate)
    vmin = np.min(growth_rate)
    im = ax.imshow(growth_rate, origin='lower', extent=[0,1,0,1],
                   vmin=vmin, vmax=vmax, aspect='auto', cmap='RdBu_r')
    cb = fig.colorbar(im, cax=cax)
    set_axes(ax, iter_min, iter_max)
    set_colorbar(cb)
    fig.tight_layout()
    name = os.path.join(path_save, 'growth_rate')
    fig.savefig(name + '.png', dpi=100)
    fig.show()

    print('\nOutput saved in ' + path_save + '\n')
