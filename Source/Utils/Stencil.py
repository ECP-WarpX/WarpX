"""
Python script to compute the minimum number of guard cells for a given error threshold,
based on the measurement of the PSATD stencil extent (that is, the minimum number of
guard cells such that the stencil measure is not larger than the error threshold).
Reference: https://arxiv.org/abs/2106.12919

Example of how to run the script for a 2D case:
    python Stencil.py --dx 1e-06 --dz 1e-06 --dt 1e-14 --nox 16 --noz 16
                      --gamma 30 --galilean --path /path/to/output --name this_test

Example of how to run the script for a 3D case:
    python Stencil.py --dx 1e-06 --dy 1e-06 --dz 1e-06 --dt 1e-14 --nox 16 --noy 16 --noz 16
                      --gamma 30 --galilean --path /path/to/output --name this_test

Some command line arguments are optional. For help, run: python Stencil.py --help
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.constants import c

plt.style.use('tableau-colorblind10')
plt.rcParams.update({'font.size': 14})
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

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
    m = order // 2
    coeffs = np.zeros(m+1)

    if (staggered):
        prod = 1.
        for k in range(1, m+1):
            prod = prod * (m+k) / (4*k)
        coeffs[0] = 4 * m * prod**2
        for n in range(1, m+1):
            coeffs[n] = - ((2*n-3) * (m+1-n)) / ((2*n-1) * (m-1+n)) * coeffs[n-1]
    else:
        coeffs[0] = -2.
        for n in range(1, m+1):
            coeffs[n] = - (m+1-n) / (m+n) * coeffs[n-1]

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
    m = order // 2
    coeffs = get_Fornberg_coeffs(order, staggered)
        
    # Array of values for n: from 1 to m
    n = np.arange(1, m+1)
    
    # Array of values of sin (first axis corresponds to k and second axis to n)
    if staggered:
        sin = np.sin(kx[:,np.newaxis] * (n[np.newaxis,:]-0.5) * dx) / ((n[np.newaxis,:]-0.5) * dx)
    else:
        sin = np.sin(kx[:,np.newaxis] * n[np.newaxis,:] * dx) / (n[np.newaxis,:] * dx)
            
    # Modified k
    k_mod = np.tensordot(sin, coeffs[1:], axes=(-1,-1))
    
    return k_mod

def func_cosine(om, w_c, dt):
    """
    Compute the leading spectral coefficient of the general PSATD equations: theta_c**2*cos(om*dt),
    where theta_c = exp(i*w_c*dt/2), w_c = v_gal*[kz]_c, om_s = c*|[k]| (and [k] or [kz]
    denote the centered or staggered modified wave vector or vector component).

    Parameters
    ----------
    om : numpy.ndarray
        Array of centered or staggered modified frequencies.
    w_c : numpy.ndarray
        Array of values of v_gal * [kz]_c.
    dt : float
        Time step.

    Returns
    -------
    coeff : numpy.ndarray
        Leading spectral coefficient of the general PSATD equations.
    """
    theta_c = np.exp(1.j * w_c * dt * 0.5)
    coeff = theta_c**2 * np.cos(om * dt)
    return coeff

def compute_stencils(coeff_nodal, coeff_stagg, dim, axis):
    """
    Compute nodal and staggered stencils along a given direction.

    Parameters
    ----------
    coeff_nodal : numpy.ndarray
        Leading spectral nodal coefficient of the general PSATD equations.
    coeff_stagg : numpy.ndarray
        Leading spectral staggered coefficient of the general PSATD equations.
    dim : int
        Number of dimensions.
    axis : int
        Axis or direction.

    Returns
    -------
    (stencil_avg_nodal, stencil_avg_stagg) : tuple
        Nodal and staggered stencils along a given direction.
    """
    # Inverse FFTs of the spectral coefficient along the chosen axis
    stencil_nodal = np.fft.ifft(coeff_nodal, axis = axis)
    stencil_stagg = np.fft.ifft(coeff_stagg, axis = axis)

    # Average results over remaining axes in spectral space
    if (dim == 2):
        if (axis == 0):
              # Averaged over kz
              stencil_avg_nodal  = np.sum(stencil_nodal, axis = 1)
              stencil_avg_nodal /= (stencil_nodal.shape[1])
              stencil_avg_stagg  = np.sum(stencil_stagg, axis = 1)
              stencil_avg_stagg /= (stencil_stagg.shape[1])
        elif (axis == 1):
              # Averaged over kx
              stencil_avg_nodal  = np.sum(stencil_nodal, axis = 0)
              stencil_avg_nodal /= (stencil_nodal.shape[0])
              stencil_avg_stagg  = np.sum(stencil_stagg, axis = 0)
              stencil_avg_stagg /= (stencil_stagg.shape[0])
    elif (dim == 3):
        if (axis == 0):
              # Averaged over ky and kz
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis = 2), axis = 1)
              stencil_avg_nodal /= (stencil_nodal.shape[2] * stencil_nodal.shape[1])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis = 2), axis = 1)
              stencil_avg_stagg /= (stencil_stagg.shape[2] * stencil_stagg.shape[1])
        elif (axis == 1):
              # Averaged over kx and kz
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis = 2), axis = 0)
              stencil_avg_nodal /= (stencil_nodal.shape[2] * stencil_nodal.shape[0])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis = 2), axis = 0)
              stencil_avg_stagg /= (stencil_stagg.shape[2] * stencil_stagg.shape[0])
        elif (axis == 2):
              # Averaged over kx and ky
              stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis = 1), axis = 0)
              stencil_avg_nodal /= (stencil_nodal.shape[1] * stencil_nodal.shape[0])
              stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis = 1), axis = 0)
              stencil_avg_stagg /= (stencil_stagg.shape[1] * stencil_stagg.shape[0])

    return (stencil_avg_nodal, stencil_avg_stagg)

def compute_all(dx, dy, dz, dt, dim, nox, noy, noz, Nx, Ny, Nz, v_gal):
    """
    Compute nodal and staggered stencils along all directions.

    Parameters
    ----------
    dx : float
        Cell size along x.
    dy : float
        Cell size along y (0 in 2D).
    dz : float
        Cell size along z.
    dt : float
        Time step.
    dim : int
        Number of dimensions.
    nox : int
        Spectral order along x.
    noy : int
        Spectral order along y (0 in 2D).
    noz : int
        Spectral order along z.
    Nx : int
        Number of mesh points along x.
    Ny : int
        Number of mesh points along y (0 in 2D).
    Nz : int
        Number of mesh points along z.
    v_gal : float
        Galilean ublic/nb/default/home?option=ngCLBeneficiaryvelocity.

    Returns
    -------
    (stencil_nodal, stencil_stagg) : tuple
        Nodal and staggered stencils along all directions.
    """

    # k vectors
    kx_arr = 2 * np.pi * np.fft.fftfreq(Nx, dx)
    kz_arr = 2 * np.pi * np.fft.fftfreq(Nz, dz)
    if (dim == 3):
        ky_arr = 2 * np.pi * np.fft.fftfreq(Ny, dy)
    
    # Centered modified k vectors
    kx_arr_c = modified_k(kx_arr, dx, nox, False) if (nox != 'inf') else kx_arr
    kz_arr_c = modified_k(kz_arr, dz, noz, False) if (noz != 'inf') else kz_arr
    if (dim == 3):
        ky_arr_c = modified_k(ky_arr, dy, noy, False) if (noy != 'inf') else ky_arr
    
    # Staggered modified k vectors
    kx_arr_s = modified_k(kx_arr, dx, nox, True) if (nox != 'inf') else kx_arr
    kz_arr_s = modified_k(kz_arr, dz, noz, True) if (noz != 'inf') else kz_arr
    if (dim == 3):
        ky_arr_s = modified_k(ky_arr, dy, noy, True) if (noy != 'inf') else ky_arr
    
    # Mesh in k space
    if (dim == 2):
        kx_c, kz_c = np.meshgrid(kx_arr_c, kz_arr_c)
        kx_s, kz_s = np.meshgrid(kx_arr_s, kz_arr_s)
    elif (dim == 3):
        kx_c, ky_c, kz_c = np.meshgrid(kx_arr_c, ky_arr_c, kz_arr_c)
        kx_s, ky_s, kz_s = np.meshgrid(kx_arr_s, ky_arr_s, kz_arr_s)
    
    # Frequencies
    if (dim == 2):
        kk_c = np.sqrt(kx_c**2 + kz_c**2)
        kk_s = np.sqrt(kx_s**2 + kz_s**2)
    elif (dim == 3):
        kk_c = np.sqrt(kx_c**2 + ky_c**2 + kz_c**2)
        kk_s = np.sqrt(kx_s**2 + ky_s**2 + kz_s**2)
    om_c = c * kk_c
    om_s = c * kk_s
    w_c = v_gal * kz_c
    
    # Spectral coefficient
    coeff_nodal = func_cosine(om_c, w_c, dt)
    coeff_stagg = func_cosine(om_s, w_c, dt)
    
    # Stencils
    stencil_x_nodal, stencil_x_stagg = compute_stencils(coeff_nodal, coeff_stagg, dim, axis = 0)
    if (dim == 2):
        stencil_z_nodal, stencil_z_stagg = compute_stencils(coeff_nodal, coeff_stagg, dim, axis = 1)
    elif (dim == 3):
        stencil_y_nodal, stencil_y_stagg = compute_stencils(coeff_nodal, coeff_stagg, dim, axis = 1)
        stencil_z_nodal, stencil_z_stagg = compute_stencils(coeff_nodal, coeff_stagg, dim, axis = 2)

    if (dim == 2):
        stencil_nodal = (stencil_x_nodal, stencil_z_nodal)
        stencil_stagg = (stencil_x_stagg, stencil_z_stagg)
    elif (dim == 3):
        stencil_nodal = (stencil_x_nodal, stencil_y_nodal, stencil_z_nodal)
        stencil_stagg = (stencil_x_stagg, stencil_y_stagg, stencil_z_stagg)

    return (stencil_nodal, stencil_stagg)

def compute_guard_cells(error, stencil):
    """
    Compute the minimum number of guard cells for a given error threshold
    (number of guard cells such that the stencil measure is not larger than the error threshold).

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
    diff = stencil - error
    minv = diff[diff > 0.].min()
    guard_cells = np.argwhere(diff == minv)[0,0]
    return guard_cells

def plot_stencil(cells, stencil_nodal, stencil_stagg, label, name):
    """
    Plot stencil extent for nodal and staggered/hybrid solver, as a function of the number of cells,
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
    fig = plt.figure(dpi = 100)
    ax = fig.add_subplot(111)
    ax.plot(cells, stencil_nodal, '-', color = cc[0], label = r'nodal')
    ax.plot(cells, stencil_stagg, '-', color = cc[1], label = r'staggered, hybrid')
    ax.set_yscale('log')
    ax.set_xticks(cells, minor = True)
    ax.grid(which = 'minor', linewidth = 0.2)
    ax.grid(which = 'major', linewidth = 0.4)
    ax.legend()
    ax.set_xlabel('number of cells')
    ax.set_ylabel(r'$\Gamma({:s})$'.format(label))
    ax.set_title(r'Stencil extent along ${:s}$'.format(label))
    fig.tight_layout()
    fig_name = args.path + './figure_stencil_' + label
    if (name):
        fig_name += '_' + name
    fig.savefig(fig_name + '.pdf', dpi = 100)
    fig.savefig(fig_name + '.png', dpi = 100)

if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--dx', type = float, required = True,
                        help = 'cell size along x')
    parser.add_argument('--dy', type = float, required = False, default = 0.,
                        help = 'cell size along y')
    parser.add_argument('--dz', type = float, required = True,
                        help = 'cell size along z')
    parser.add_argument('--dt', type = float, required = True,
                        help = 'time step')
    parser.add_argument('--nox', type = int, required = True,
                        help = 'spectral order along x')
    parser.add_argument('--noy', type = int, required = False, default = 0,
                        help = 'spectral order along y')
    parser.add_argument('--noz', type = int, required = True,
                        help = 'spectral order along z')
    parser.add_argument('--gamma', type = float, required = False, default = 1.,
                        help = 'Lorentz factor')
    parser.add_argument('--galilean', action = 'store_true', default = False,
                        help = 'add this for Galilean schemes')
    parser.add_argument('--path', type = str, required = False, default = './',
                        help = 'path where figures are saved')
    parser.add_argument('--name', type = str, required = False, default = '',
                        help = 'common label for figure names')
    args = parser.parse_args()

    # Add trailing '/' to args.path, if necessary
    if (args.path[-1] != '/'):
        args.path += '/'

    if (args.dy != 0 and args.noy == 0):
        sys.exit('You entered --dy, but not --noy. Please enter either both or none.')

    if (args.dy == 0 and args.noy != 0):
        sys.exit('You entered --noy, but not --dy. Please enter either both or none.')

    # Number of dimensions
    dim = 2 if (args.dy == 0 and args.noy == 0) else 3
    zid = dim - 1

    # Galilean velocity (default = 0.)
    v_gal = 0.
    if (args.galilean):
        v_gal = - np.sqrt(1. - 1./args.gamma**2) * c
    
    # Display some output
    print('\n >> Numerical Setup')
    print('    ---------------')
    print('\n    cell size:')
    print('    - dx = {:g}'.format(args.dx))
    if (dim == 3):
        print('    - dy = {:g}'.format(args.dy))
    print('    - dz = {:g}'.format(args.dz))
    print('    - dz/dx = {:g}'.format(args.dz / args.dx))
    if (dim == 3):
        print('    - dz/dy = {:g}'.format(args.dz / args.dy))
    print('\n    time step:')
    print('    - dt = {:g}'.format(args.dt))
    print('    - c*dt/dx = {:g}'.format(c * args.dt / args.dx))
    if (dim == 3):
        print('    - c*dt/dy = {:g}'.format(c * args.dt / args.dy))
    print('    - c*dt/dz = {:g}'.format(c * args.dt / args.dz))
    print('\n    spectral order:')
    print('    - nox = {:d}'.format(args.nox))
    if (dim == 3):
        print('    - noy = {:d}'.format(args.noy))
    print('    - noz = {:d}'.format(args.noz))
    print('\n    Lorentz boost, Galilean velocity:')
    print('    - gamma = {:g}'.format(args.gamma))
    print('    - v_gal = {:g}'.format(v_gal))

    # Number of grid points (arbitrary mesh)
    Nx = 256
    Nz = 256
    Ny = 256 if (dim == 3) else 0

    print('\n >> Compute Stencil')
    print('    ---------------')

    stencil_nodal, stencil_stagg = compute_all(args.dx, args.dy, args.dz, args.dt, dim,
                                               args.nox, args.noy, args.noz, Nx, Ny, Nz, v_gal)

    stencil_x_nodal = stencil_nodal[0]
    stencil_z_nodal = stencil_nodal[zid]
    if (dim == 3):
        stencil_y_nodal = stencil_nodal[1]
    
    stencil_x_stagg = stencil_stagg[0]
    stencil_z_stagg = stencil_stagg[zid]
    if (dim == 3):
        stencil_y_stagg = stencil_stagg[1]
 
    # Maximum number of cells
    nx = 65
    nz = 65
    if (dim == 3):
        ny = 65

    # Array of cell numbers
    cx = np.arange(nx)
    cz = np.arange(nz)
    if (dim == 3):
        cy = np.arange(ny)

    # Arrays of stencils
    sx_nodal = abs(stencil_x_nodal[:nx])
    sx_stagg = abs(stencil_x_stagg[:nx])
    sz_nodal = abs(stencil_z_nodal[:nz])
    sz_stagg = abs(stencil_z_stagg[:nz])
    if (dim == 3):
        sy_nodal = abs(stencil_y_nodal[:ny])
        sy_stagg = abs(stencil_y_stagg[:ny])

    # Compute minimum number of guard cells for given error threshold
    # (number of guard cells such that the stencil measure is not larger than the error threshold)

    # 1) Error threshold
    err_x = 1e-07
    err_z = 1e-07
    if (dim == 3):
        err_y = 1e-07

    print('\n    error threshold:')
    print('    - err_x = {:g}'.format(err_x))
    if (dim == 3):
        print('    - err_y = {:g}'.format(err_y))
    print('    - err_z = {:g}'.format(err_z))

    # 2) Nodal solver
    gc_nodal_x = compute_guard_cells(err_x, sx_nodal)
    gc_nodal_z = compute_guard_cells(err_z, sz_nodal)
    if (dim == 3):
        gc_nodal_y = compute_guard_cells(err_y, sy_nodal)

    print('\n    nodal solver:')
    print('    - {:d} guard cells along x'.format(gc_nodal_x))
    if (dim == 3):
        print('    - {:d} guard cells along y'.format(gc_nodal_y))
    print('    - {:d} guard cells along z'.format(gc_nodal_z))

    # 3) Staggered solver
    gc_stagg_x = compute_guard_cells(err_x, sx_stagg)
    gc_stagg_z = compute_guard_cells(err_z, sz_stagg)
    if (dim == 3):
        gc_stagg_y = compute_guard_cells(err_y, sy_stagg)

    print('\n    staggered solver:')
    print('    - {:d} guard cells along x'.format(gc_stagg_x))
    if (dim == 3):
        print('    - {:d} guard cells along y'.format(gc_stagg_y))
    print('    - {:d} guard cells along z'.format(gc_stagg_z))

    # Plot stencils
    plot_stencil(cx, sx_nodal, sx_stagg, 'x', args.name)
    plot_stencil(cz, sz_nodal, sz_stagg, 'z', args.name)
    if (dim == 3):
        plot_stencil(cy, sy_nodal, sy_stagg, 'y', args.name)

    print('\n >> Path to Figures: \'' + os.path.abspath(args.path) + '\'\n')
