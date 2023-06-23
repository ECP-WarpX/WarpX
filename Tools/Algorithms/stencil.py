"""
Python script to compute the minimum number of guard cells for a given
error threshold, based on the measurement of the PSATD stencil extent
(that is, the minimum number of guard cells such that the stencil
measure is not larger than the error threshold).
Reference: https://doi.org/10.1016/j.cpc.2022.108457

How to run the script:
$ python stencil.py --path path_to_input_file
or, using IPython,
$ run stencil.py --path path_to_input_file
"""

import argparse
import os
from parser import parse_input

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c

plt.style.use('tableau-colorblind10')
plt.rcParams.update({'font.size': 14})

sp = np.finfo(np.float32).eps
dp = np.finfo(np.float64).eps

def get_Fornberg_coeffs(order, staggered):
    """
    Compute the centered or staggered Fornberg coefficients at finite order.

    Parameters
    ----------
    order : int
        Finite order of the approximation.
    staggered : bool
        Whether to compute the centered or staggered Fornberg coefficients.

    Returns
    -------
    coeffs : numpy.ndarray
        Array of centered or staggered Fornberg coefficients.
    """
    m = order//2
    coeffs = np.zeros(m+1)

    # Compute Fornberg coefficients by recurrence
    if staggered:
        prod = 1.
        for k in range(1, m+1):
            prod = prod*(m+k)/(4*k)
        coeffs[0] = 4*m*prod**2
        for n in range(1, m+1):
            coeffs[n] = -(((2*n-3)*(m+1-n))/((2*n-1)*(m-1+n))*coeffs[n-1])
    else:
        coeffs[0] = -2.
        for n in range(1, m+1):
            coeffs[n] = -(m+1-n)/(m+n)*coeffs[n-1]

    return coeffs

def modified_k(kx, dx, order, staggered):
    """
    Compute the centered or staggered modified wave vector at finite order.

    Parameters
    ----------
    kx : numpy.ndarray
        Standard wave vector.
    dx : float
        Cell size in real space.
    order : int
        Finite order of the approximation.
    staggered : bool
        Whether to compute the centered or staggered modified wave vector.

    Returns
    -------
    k_mod : numpy.ndarray
        Centered or staggered modified wave vector.
    """
    m = order//2
    coeffs = get_Fornberg_coeffs(order, staggered)

    # Array of values for n: from 1 to m
    n = np.arange(1, m+1)

    # Array of values of sin
    # (first axis corresponds to k and second axis to n)
    if staggered:
        sin_kn = (np.sin(kx[:,np.newaxis]*(n[np.newaxis,:]-0.5)*dx)/((n[np.newaxis,:]-0.5)*dx))
    else:
        sin_kn = (np.sin(kx[:,np.newaxis]*n[np.newaxis,:]*dx)/(n[np.newaxis,:]*dx))

    # Modified k
    k_mod = np.tensordot(sin_kn, coeffs[1:], axes=(-1,-1))

    return k_mod

def func_cosine(om, w_c, dt):
    """
    Compute the leading spectral coefficient of the general PSATD equations:
    theta_c**2*cos(om*dt), where theta_c = exp(i*w_c*dt/2), w_c = v_gal*[kz]_c,
    om_s = c*|[k]| (and [k] or [kz] denote the centered or staggered modified
    wave vector or vector component).

    Parameters
    ----------
    om : numpy.ndarray
        Array of centered or staggered modified frequencies.
    w_c : numpy.ndarray
        Array of values of v_gal*[kz]_c.
    dt : float
        Time step.

    Returns
    -------
    coeff : numpy.ndarray
        Leading spectral coefficient of the general PSATD equations.
    """
    theta_c = np.exp(1.j*w_c*dt*0.5)
    coeff = theta_c**2*np.cos(om*dt)
    return coeff

def compute_stencils(coeff_nodal, coeff_stagg, axis):
    """
    Compute nodal and staggered stencils along a given direction.

    Parameters
    ----------
    coeff_nodal : numpy.ndarray
        Leading spectral nodal coefficient of the general PSATD equations.
    coeff_stagg : numpy.ndarray
        Leading spectral staggered coefficient of the general PSATD equations.
    axis : int
        Axis or direction (must be 0, 1, or 2).

    Returns
    -------
    (stencil_avg_nodal, stencil_avg_stagg) : tuple
        Nodal and staggered stencils along a given direction.
    """
    # Inverse FFTs of the spectral coefficient along the chosen axis
    stencil_nodal = np.fft.ifft(coeff_nodal, axis=axis)
    stencil_stagg = np.fft.ifft(coeff_stagg, axis=axis)

    # Average results over remaining axes in spectral space
    if dims == 1:
        stencil_avg_nodal = stencil_nodal
        stencil_avg_stagg = stencil_stagg
    elif dims == 2:
        if axis == 0:
              # Averaged over kz
              stencil_avg_nodal  = np.sum(stencil_nodal, axis=-1)
              stencil_avg_nodal /= stencil_nodal.shape[-1]
              stencil_avg_stagg  = np.sum(stencil_stagg, axis=-1)
              stencil_avg_stagg /= stencil_stagg.shape[-1]
        elif axis == 1 or axis == -1:
              # Averaged over kx
              stencil_avg_nodal  = np.sum(stencil_nodal, axis=0)
              stencil_avg_nodal /= stencil_nodal.shape[0]
              stencil_avg_stagg  = np.sum(stencil_stagg, axis=0)
              stencil_avg_stagg /= stencil_stagg.shape[0]
    elif dims == 3:
        if axis == 0:
              # Averaged over ky and kz
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=-1), axis=1)
              stencil_avg_nodal /= (stencil_nodal.shape[-1]*stencil_nodal.shape[1])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=-1), axis=1)
              stencil_avg_stagg /= (stencil_stagg.shape[-1]*stencil_stagg.shape[1])
        elif axis == 1:
              # Averaged over kx and kz
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=-1), axis=0)
              stencil_avg_nodal /= (stencil_nodal.shape[-1]*stencil_nodal.shape[0])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=-1), axis=0)
              stencil_avg_stagg /= (stencil_stagg.shape[-1]*stencil_stagg.shape[0])
        elif axis == 2 or axis == -1:
              # Averaged over kx and ky
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=1), axis=0)
              stencil_avg_nodal /= (stencil_nodal.shape[1]*stencil_nodal.shape[0])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=1), axis=0)
              stencil_avg_stagg /= (stencil_stagg.shape[1]*stencil_stagg.shape[0])

    stencils = dict()
    stencils['nodal'] = abs(stencil_avg_nodal)
    stencils['stagg'] = abs(stencil_avg_stagg)

    return stencils

def compute_all(dx_boosted, dt, psatd_order, v_gal, nx=None):
    """
    Compute nodal and staggered stencils along all directions.

    Parameters
    ----------
    dx_boosted : np.ndarray (float)
        Cell size along each direction.
    dt : float
        Time step.
    psatd_order : np.ndarray (int)
        Spectral order along each direction.
    v_gal : float
        Galilean velocity.
    nx : np.ndarray (int), optional (default = 256 in each direction)
        Number of mesh points along each direction.

    Returns
    -------
    (stencil_nodal, stencil_stagg) : tuple
        Nodal and staggered stencils along all directions.
    """
    # Number of dimensions
    dims = len(dx_boosted)

    # Default value for nx
    if nx is None:
        nx = np.full(shape=dims, fill_value=256)

    # k vectors and modified k vectors
    k_arr   = []
    k_arr_c = []
    k_arr_s = []
    for i in range(dims):
        k_arr.append(2*np.pi*np.fft.fftfreq(nx[i], dx_boosted[i]))
        if psatd_order[i] != 'inf':
            k_arr_c.append(modified_k(k_arr[i], dx_boosted[i], psatd_order[i], False))
            k_arr_s.append(modified_k(k_arr[i], dx_boosted[i], psatd_order[i], True))
        else:
            k_arr_c.append(k_arr[i])
            k_arr_s.append(k_arr[i])

    # Mesh in k space
    k_c = np.meshgrid(*k_arr_c)
    k_s = np.meshgrid(*k_arr_s)
    kk_c = 0.
    kk_s = 0.
    for k in k_c:
        kk_c += k**2
    for k in k_s:
        kk_s += k**2
    kk_c = np.sqrt(kk_c)
    kk_s = np.sqrt(kk_s)

    # Frequencies
    om_c = c*kk_c
    om_s = c*kk_s
    w_c = v_gal*k_c[-1]

    # Spectral coefficient
    coeff_nodal = func_cosine(om_c, w_c, dt)
    coeff_stagg = func_cosine(om_s, w_c, dt)

    # Stencils
    stencils = dict()
    stencils['x'] = compute_stencils(coeff_nodal, coeff_stagg, axis=0)
    if dims == 3:
        stencils['y'] = compute_stencils(coeff_nodal, coeff_stagg, axis=1)
    if dims > 1:
        stencils['z'] = compute_stencils(coeff_nodal, coeff_stagg, axis=-1)

    return stencils

def compute_guard_cells(errmin, errmax, stencil):
    """
    Compute the minimum number of guard cells for a given error threshold
    (number of guard cells such that the stencil measure is not larger
    than the error threshold).

    Parameters
    ----------
    error : float
        Error threshold.
    stencil : numpy.ndarray
        Stencil array.

    Returns
    -------
    guard_cells : numpy.int64
        Number of cells.
    """
    diff = stencil - errmin
    v = next(d for d in diff if d < 0)
    gcmin = np.argwhere(diff == v)[0,0]
    diff = stencil - errmax
    try:
        v = next(d for d in diff if d < 0)
        gcmax = np.argwhere(diff == v)[0,0] - 1
    except StopIteration:
        gcmin, gcmax = compute_guard_cells(errmin, errmax*10, stencil)
    return (gcmin, gcmax)

def plot_stencil(cells, stencil_nodal, stencil_stagg, label, path, name):
    """
    Plot stencil extent for nodal and staggered/hybrid solver,
    as a function of the number of cells.

    Parameters
    ----------
    cells : numpy.ndarray
        Array of cell numbers.
    stencil_nodal : numpy.ndarray
        Stencil array for the nodal solver.
    stencil_stagg : numpy.ndarray
        Stencil array for the staggered or hybrid solver.
    label : str
        Label for plot annotations.
    name : str
        Label for figure name.
    """
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(111)
    ax.plot(cells, stencil_nodal, linestyle='-', label='nodal')
    ax.plot(cells, stencil_stagg, linestyle='-', label='staggered or hybrid')
    # Plot single and double precision machine epsilons
    ax.axhline(y=sp, c='grey', ls='dashed', label='machine epsilon (single precision)')
    ax.axhline(y=dp, c='grey', ls='dotted', label='machine epsilon (double precision)')
    # Shade regions between single and double precision machine epsilons
    xmin, xmax = compute_guard_cells(sp, dp, stencil_nodal)
    ax.fill_between(cells[xmin:xmax+1], stencil_nodal[xmin:xmax+1], alpha=0.5)
    xmin, xmax = compute_guard_cells(sp, dp, stencil_stagg)
    ax.fill_between(cells[xmin:xmax+1], stencil_stagg[xmin:xmax+1], alpha=0.5)
    #
    ax.set_yscale('log')
    ax.set_xticks(cells, minor=True)
    ax.grid(which='minor', linewidth=0.2)
    ax.grid(which='major', linewidth=0.4)
    ax.legend()
    ax.set_xlabel('number of cells')
    ax.set_ylabel('signal to be truncated')
    ax.set_title(r'Stencil extent along ${:s}$'.format(label))
    fig.tight_layout()
    fig_name = os.path.join(path, 'figure_stencil_' + label)
    if name:
        fig_name += '_' + name
    fig.savefig(fig_name + '.png', dpi=150)

def run_main(dims, dx_boosted, dt, psatd_order, gamma=1., galilean=False, path='.', name=''):
    """
    Main function.

    Parameters
    ----------
    dims : int
        Number of dimensions.
    dx_boosted : numpy.ndarray (float)
        Cell size along each direction.
    dt : float
        Time step.
    psatd_order : numpy.ndarray (int)
        Spectral order along each direction.
    gamma : float, optional (default = 1.)
        Lorentz factor.
    galilean : bool, optional (default = False)
        Galilean scheme.
    path : str, optional (default = '.')
        Path where figures are saved.
    name : str, optional (default = '')
        Common label for figure names.
    """

    # Galilean velocity (default = 0.)
    v_gal = 0.
    if galilean:
        v_gal = -np.sqrt(1.-1./gamma**2)*c

    # Display some output
    print('\nCell size:')
    print(f'- dx = {dx_boosted}')
    if dims > 1:
        print(f'- dx[1:]/dx[0] = {dx_boosted[1:]/dx_boosted[0]}')
    print('\nTime step:')
    print(f'- dt = {dt}')
    print(f'- c*dt/dx = {c*dt/dx_boosted}')
    print('\nSpectral order:')
    print(f'- order = {psatd_order}')
    print('\nLorentz boost, Galilean velocity:')
    print(f'- gamma = {gamma}')
    print(f'- v_gal = {v_gal}')

    stencils = compute_all(dx_boosted, dt, psatd_order, v_gal)

    # Maximum number of cells
    ncx = 65
    if dims == 3:
        ncy = 65
    if dims > 1:
        ncz = 65

    # Array of cell numbers
    cx = np.arange(ncx)
    if dims == 3:
        cy = np.arange(ncy)
    if dims > 1:
        cz = np.arange(ncz)

    # Arrays of stencils
    stencils['x']['nodal'] = stencils['x']['nodal'][:ncx]
    stencils['x']['stagg'] = stencils['x']['stagg'][:ncx]
    if dims == 3:
        stencils['y']['nodal'] = stencils['y']['nodal'][:ncy]
        stencils['y']['stagg'] = stencils['y']['stagg'][:ncy]
    if dims > 1:
        stencils['z']['nodal'] = stencils['z']['nodal'][:ncz]
        stencils['z']['stagg'] = stencils['z']['stagg'][:ncz]

    # Plot stencils
    plot_stencil(cx, stencils['x']['nodal'], stencils['x']['stagg'], 'x', path, name)
    if dims == 3:
        plot_stencil(cy, stencils['y']['nodal'], stencils['y']['stagg'], 'y', path, name)
    if dims > 1:
        plot_stencil(cz, stencils['z']['nodal'], stencils['z']['stagg'], 'z', path, name)

    # Compute min and max numbers of guard cells
    gcmin_x_nodal, gcmax_x_nodal = compute_guard_cells(sp, dp, stencils['x']['nodal'])
    if dims == 3:
        gcmin_y_nodal, gcmax_y_nodal = compute_guard_cells(sp, dp, stencils['y']['nodal'])
    if dims > 1:
        gcmin_z_nodal, gcmax_z_nodal = compute_guard_cells(sp, dp, stencils['z']['nodal'])
    #
    gcmin_x_stagg, gcmax_x_stagg = compute_guard_cells(sp, dp, stencils['x']['stagg'])
    if dims == 3:
        gcmin_y_stagg, gcmax_y_stagg = compute_guard_cells(sp, dp, stencils['y']['stagg'])
    if dims > 1:
        gcmin_z_stagg, gcmax_z_stagg = compute_guard_cells(sp, dp, stencils['z']['stagg'])

    fig_path = os.path.abspath(path)
    print(f'\nFigures saved in {fig_path}/.')
    print('\nThe plots show the extent of the signal to be truncated (y-axis)'
        + '\nby choosing a given number of cells (x-axis) for the ghost regions'
        + '\nof each simulation grid, along x, y, and z.')
    print('\nIt is recommended to choose a number of ghost cells that corresponds to'
        + '\na truncation of the signal between single and double machine precision.'
        + '\nThe more ghost cells, the more accurate, yet expensive, results.'
        + '\nFor each stencil the region of accuracy between single and double precision'
        + '\nis shaded to help you identify a suitable number of ghost cells.')
    print('\nFor a nodal simulation, choose:')
    print(f'- between {gcmin_x_nodal} and {gcmax_x_nodal} ghost cells along x')
    if dims == 3:
        print(f'- between {gcmin_y_nodal} and {gcmax_y_nodal} ghost cells along y')
    if dims > 1:
        print(f'- between {gcmin_z_nodal} and {gcmax_z_nodal} ghost cells along z')
    print('\nFor a staggered or hybrid simulation, choose:')
    print(f'- between {gcmin_x_stagg} and {gcmax_x_stagg} ghost cells along x')
    if dims == 3:
        print(f'- between {gcmin_y_stagg} and {gcmax_y_stagg} ghost cells along y')
    if dims > 1:
        print(f'- between {gcmin_z_stagg} and {gcmax_z_stagg} ghost cells along z')
    print()

    return

if __name__ == '__main__':

    # Parse path to input file from command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--path')
    args = parser.parse_args()
    input_file = args.path

    # Parse input file
    input_dict = parse_input(input_file)

    # TODO Handle RZ
    dims = int(input_dict['geometry.dims'][0])

    # Notation considering x as vector of coordinates (x,y,z)
    nx = np.array([int(w) for w in input_dict['amr.n_cell']])
    xmin = np.array([float(w) for w in input_dict['geometry.prob_lo']])
    xmax = np.array([float(w) for w in input_dict['geometry.prob_hi']])

    # Cell size in the lab frame and boosted frame (boost along z)
    ## lab frame
    dx = (xmax-xmin) / nx
    ## boosted frame
    gamma = 1.
    if 'warpx.gamma_boost' in input_dict:
        gamma = float(input_dict['warpx.gamma_boost'][0])
    beta = np.sqrt(1. - 1./gamma**2)
    dx_boosted = np.copy(dx)
    dx_boosted[-1] = (1. + beta) * gamma * dx[-1]

    # Time step for pseudo-spectral scheme
    cfl = 0.999
    if 'warpx.cfl' in input_dict:
        cfl = float(input_dict['warpx.cfl'][0])
    dt = cfl * np.min(dx_boosted) / c

    # Pseudo-spectral order
    psatd_order = np.full(shape=dims, fill_value=16)
    if 'psatd.nox' in input_dict:
        psatd_order[0] = int(input_dict['psatd.nox'][0])
    if 'psatd.noy' in input_dict:
        psatd_order[1] = int(input_dict['psatd.noy'][0])
    if 'psatd.noz' in input_dict:
        psatd_order[-1] = int(input_dict['psatd.noz'][0])

    # Galilean flag
    # TODO Extend for alternative syntax psatd.v_galilean
    galilean = False
    if 'psatd.use_default_v_galilean' in input_dict:
        galilean = bool(input_dict['psatd.use_default_v_galilean'][0])
    if 'psatd.v_galilean' in input_dict:
        galilean = bool(input_dict['psatd.v_galilean'][-1])

    # Run main function (some arguments are optional,
    # see definition of run_main function for help)
    run_main(dims, dx_boosted, dt, psatd_order, gamma, galilean)
