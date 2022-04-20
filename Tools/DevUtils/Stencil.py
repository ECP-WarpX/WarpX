"""
Python script to compute the minimum number of guard cells for a given
error threshold, based on the measurement of the PSATD stencil extent
(that is, the minimum number of guard cells such that the stencil
measure is not larger than the error threshold).
Reference: https://arxiv.org/abs/2106.12919

Run the script simply with "python Stencil.py" (or with "run Stencil.py"
using IPython). The user can modify the input parameters set in the main
function at the end of the file.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c

plt.style.use('tableau-colorblind10')
plt.rcParams.update({'font.size': 14})

def get_Fornberg_coeffs(order, staggered):
    """
    Compute the centered or staggered
    Fornberg coefficients at finite order.

    Parameters
    ----------
    order : int
        Finite order of the approximation.
    staggered : bool
        Whether to compute the centered or staggered
        Fornberg coefficients.

    Returns
    -------
    coeffs : numpy.ndarray
        Array of centered or staggered Fornberg coefficients.
    """
    m = order // 2
    coeffs = np.zeros(m+1)

    # Compute Fornberg coefficients by recurrence
    if staggered:
        prod = 1.
        for k in range(1, m+1):
            prod = prod * (m+k) / (4*k)
        coeffs[0] = 4 * m * prod**2
        for n in range(1, m+1):
            coeffs[n] = - (((2*n-3) * (m+1-n))
                        / ((2*n-1) * (m-1+n)) * coeffs[n-1])
    else:
        coeffs[0] = -2.
        for n in range(1, m+1):
            coeffs[n] = - (m+1-n) / (m+n) * coeffs[n-1]

    return coeffs

def modified_k(kx, dx, order, staggered):
    """
    Compute the centered or staggered
    modified wave vector at finite order.

    Parameters
    ----------
    kx : numpy.ndarray
        Standard wave vector.
    dx : float
        Cell size in real space.
    order : int
        Finite order of the approximation.
    staggered : bool
        Whether to compute the centered or staggered
        modified wave vector.

    Returns
    -------
    k_mod : numpy.ndarray
        Centered or staggered modified wave vector.
    """
    m = order // 2
    coeffs = get_Fornberg_coeffs(order, staggered)

    # Array of values for n: from 1 to m
    n = np.arange(1, m+1)

    # Array of values of sin
    # (first axis corresponds to k and second axis to n)
    if staggered:
        sin_kn = (np.sin(kx[:,np.newaxis]*(n[np.newaxis,:]-0.5) * dx)
                 / ((n[np.newaxis,:]-0.5) * dx))
    else:
        sin_kn = (np.sin(kx[:,np.newaxis]*n[np.newaxis,:] * dx)
                 / (n[np.newaxis,:] * dx))

    # Modified k
    k_mod = np.tensordot(sin_kn, coeffs[1:], axes=(-1,-1))

    return k_mod

def func_cosine(om, w_c, dt):
    """
    Compute the leading spectral coefficient of the general
    PSATD equations: theta_c**2*cos(om*dt), where
    theta_c = exp(i*w_c*dt/2), w_c = v_gal*[kz]_c,
    om_s = c*|[k]| (and [k] or [kz] denote the centered or
    staggered modified wave vector or vector component).

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
    theta_c = np.exp(1.j*w_c*dt*0.5)
    coeff = theta_c**2 * np.cos(om*dt)
    return coeff

def compute_stencils(coeff_nodal, coeff_stagg, axis):
    """
    Compute nodal and staggered stencils along a given direction.

    Parameters
    ----------
    coeff_nodal : numpy.ndarray
        Leading spectral nodal coefficient of the general
        PSATD equations.
    coeff_stagg : numpy.ndarray
        Leading spectral staggered coefficient of the general
        PSATD equations.
    axis : int
        Axis or direction.

    Returns
    -------
    (stencil_avg_nodal, stencil_avg_stagg) : tuple
        Nodal and staggered stencils along a given direction.
    """
    # Inverse FFTs of the spectral coefficient along the chosen axis
    stencil_nodal = np.fft.ifft(coeff_nodal, axis=axis)
    stencil_stagg = np.fft.ifft(coeff_stagg, axis=axis)

    # Average results over remaining axes in spectral space
    if axis == 0:
          # Averaged over ky and kz
          stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=2), axis=1)
          stencil_avg_nodal /= (stencil_nodal.shape[2]*stencil_nodal.shape[1])
          stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=2), axis=1)
          stencil_avg_stagg /= (stencil_stagg.shape[2]*stencil_stagg.shape[1])
    elif axis == 1:
          # Averaged over kx and kz
          stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=2), axis=0)
          stencil_avg_nodal /= (stencil_nodal.shape[2]*stencil_nodal.shape[0])
          stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=2), axis=0)
          stencil_avg_stagg /= (stencil_stagg.shape[2]*stencil_stagg.shape[0])
    elif axis == 2:
          # Averaged over kx and ky
          stencil_avg_nodal  = np.sum(np.sum(stencil_nodal, axis=1), axis=0)
          stencil_avg_nodal /= (stencil_nodal.shape[1]*stencil_nodal.shape[0])
          stencil_avg_stagg  = np.sum(np.sum(stencil_stagg, axis=1), axis=0)
          stencil_avg_stagg /= (stencil_stagg.shape[1]*stencil_stagg.shape[0])

    stencils = dict()
    stencils['nodal'] = abs(stencil_avg_nodal)
    stencils['stagg'] = abs(stencil_avg_stagg)

    return stencils

def compute_all(dx, dy, dz, dt, nox, noy, noz, v_gal, nx=256, ny=256, nz=256):
    """
    Compute nodal and staggered stencils along all directions.

    Parameters
    ----------
    dx : float
        Cell size along x.
    dy : float
        Cell size along y.
    dz : float
        Cell size along z.
    dt : float
        Time step.
    nox : int
        Spectral order along x.
    noy : int
        Spectral order along y.
    noz : int
        Spectral order along z.
    v_gal : float
        Galilean velocity.
    nx : int, optional (default = 256)
        Number of mesh points along x.
    ny : int, optional (default = 256)
        Number of mesh points along y.
    nz : int, optional (default = 256)
        Number of mesh points along z.

    Returns
    -------
    (stencil_nodal, stencil_stagg) : tuple
        Nodal and staggered stencils along all directions.
    """
    # k vectors
    kx_arr = 2 * np.pi * np.fft.fftfreq(nx, dx)
    ky_arr = 2 * np.pi * np.fft.fftfreq(ny, dy)
    kz_arr = 2 * np.pi * np.fft.fftfreq(nz, dz)

    # Centered modified k vectors
    kx_arr_c = modified_k(kx_arr, dx, nox, False) if nox != 'inf' else kx_arr
    ky_arr_c = modified_k(ky_arr, dy, noy, False) if noy != 'inf' else ky_arr
    kz_arr_c = modified_k(kz_arr, dz, noz, False) if noz != 'inf' else kz_arr

    # Staggered modified k vectors
    kx_arr_s = modified_k(kx_arr, dx, nox, True) if nox != 'inf' else kx_arr
    ky_arr_s = modified_k(ky_arr, dy, noy, True) if noy != 'inf' else ky_arr
    kz_arr_s = modified_k(kz_arr, dz, noz, True) if noz != 'inf' else kz_arr

    # Mesh in k space
    kx_c, ky_c, kz_c = np.meshgrid(kx_arr_c, ky_arr_c, kz_arr_c)
    kx_s, ky_s, kz_s = np.meshgrid(kx_arr_s, ky_arr_s, kz_arr_s)

    # Frequencies
    kk_c = np.sqrt(kx_c**2 + ky_c**2 + kz_c**2)
    kk_s = np.sqrt(kx_s**2 + ky_s**2 + kz_s**2)
    om_c = c * kk_c
    om_s = c * kk_s
    w_c = v_gal * kz_c

    # Spectral coefficient
    coeff_nodal = func_cosine(om_c, w_c, dt)
    coeff_stagg = func_cosine(om_s, w_c, dt)

    # Stencils
    stencils = dict()
    stencils['x'] = compute_stencils(coeff_nodal, coeff_stagg, axis=0)
    stencils['y'] = compute_stencils(coeff_nodal, coeff_stagg, axis=1)
    stencils['z'] = compute_stencils(coeff_nodal, coeff_stagg, axis=2)

    return stencils

def compute_guard_cells(error, stencil):
    """
    Compute the minimum number of guard cells for a given
    error threshold (number of guard cells such that the
    stencil measure is not larger than the error threshold).

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
    fig = plt.figure(dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(cells, stencil_nodal, '-', label='nodal')
    ax.plot(cells, stencil_stagg, '-', label='staggered, hybrid')
    ax.set_yscale('log')
    ax.set_xticks(cells, minor=True)
    ax.grid(which='minor', linewidth=0.2)
    ax.grid(which='major', linewidth=0.4)
    ax.legend()
    ax.set_xlabel('number of cells')
    ax.set_ylabel(r'$\Gamma({:s})$'.format(label))
    ax.set_title(r'Stencil extent along ${:s}$'.format(label))
    fig.tight_layout()
    fig_name = os.path.join(path, 'figure_stencil_' + label)
    if name:
        fig_name += '_' + name
    fig.savefig(fig_name + '.pdf', dpi=100)
    fig.savefig(fig_name + '.png', dpi=100)

def run_main(dx, dy, dz, dt, nox, noy, noz, gamma=1., galilean=False,
             ex=1e-07, ey=1e-07, ez=1e-07, path='.', name=''):
    """
    Main function.

    Parameters
    ----------
    dx : float
        Cell size along x.
    dy : float
        Cell size along y.
    dz : float
        Cell size along z.
    dt : float
        Time step.
    nox : int
        Spectral order along x.
    noy : int
        Spectral order along y.
    noz : int
        Spectral order along z.
    gamma : float, optional (default = 1.)
        Lorentz factor.
    galilean : bool, optional (default = False)
        Galilean scheme.
    ex : float, optional (default = 1e-07)
        Error threshold along x.
    ey : float, optional (default = 1e-07)
        Error threshold along y.
    ez : float, optional (default = 1e-07)
        Error threshold along z.
    path : str, optional (default = '.')
        Path where figures are saved.
    name : str, optional (default = '')
        Common label for figure names.

    Returns
    -------
    stencils : dict
        Dictionary of nodal and staggered stencils along all directions.
        Its element of the dictionary is a dictionary itself,
        containing numpy.ndarray objects.
        Keys: stencils.keys() = dict_keys(['x', 'y', 'z'])
              stencils['x'].keys() = dict_keys(['nodal', 'stagg'])
              stencils['y'].keys() = dict_keys(['nodal', 'stagg'])
              stencils['z'].keys() = dict_keys(['nodal', 'stagg'])
    """
    # Galilean velocity (default = 0.)
    v_gal = 0.
    if galilean:
        v_gal = - np.sqrt(1.-1./gamma**2) * c

    # Display some output
    print('\n >> Numerical Setup')
    print('    ---------------')
    print('\n    cell size:')
    print('    - dx = {:g}'.format(dx))
    print('    - dy = {:g}'.format(dy))
    print('    - dz = {:g}'.format(dz))
    print('    - dz/dx = {:g}'.format(dz / dx))
    print('    - dz/dy = {:g}'.format(dz / dy))
    print('\n    time step:')
    print('    - dt = {:g}'.format(dt))
    print('    - c*dt/dx = {:g}'.format(c * dt / dx))
    print('    - c*dt/dy = {:g}'.format(c * dt / dy))
    print('    - c*dt/dz = {:g}'.format(c * dt / dz))
    print('\n    spectral order:')
    print('    - nox = {:d}'.format(nox))
    print('    - noy = {:d}'.format(noy))
    print('    - noz = {:d}'.format(noz))
    print('\n    Lorentz boost, Galilean velocity:')
    print('    - gamma = {:g}'.format(gamma))
    print('    - v_gal = {:g}'.format(v_gal))

    print('\n >> Compute Stencil')
    print('    ---------------')

    stencils = compute_all(dx, dy, dz, dt, nox, noy, noz, v_gal)

    # Maximum number of cells
    ncx = 65
    ncy = 65
    ncz = 65

    # Array of cell numbers
    cx = np.arange(ncx)
    cy = np.arange(ncy)
    cz = np.arange(ncz)

    # Arrays of stencils
    stencils['x']['nodal'] = stencils['x']['nodal'][:ncx]
    stencils['x']['stagg'] = stencils['x']['stagg'][:ncx]
    stencils['y']['nodal'] = stencils['y']['nodal'][:ncy]
    stencils['y']['stagg'] = stencils['y']['stagg'][:ncy]
    stencils['z']['nodal'] = stencils['z']['nodal'][:ncz]
    stencils['z']['stagg'] = stencils['z']['stagg'][:ncz]

    # Compute minimum number of guard cells for given error threshold
    # (number of guard cells such that the stencil measure
    # is not larger than the error threshold)

    print('\n    error threshold:')
    print('    - ex = {:g}'.format(ex))
    print('    - ey = {:g}'.format(ey))
    print('    - ez = {:g}'.format(ez))

    # 1) Nodal solver
    gc_nodal_x = compute_guard_cells(ex, stencils['x']['nodal'])
    gc_nodal_y = compute_guard_cells(ey, stencils['y']['nodal'])
    gc_nodal_z = compute_guard_cells(ez, stencils['z']['nodal'])

    print('\n    nodal solver:')
    print('    - {:d} guard cells along x'.format(gc_nodal_x))
    print('    - {:d} guard cells along y'.format(gc_nodal_y))
    print('    - {:d} guard cells along z'.format(gc_nodal_z))

    # 2) Staggered solver
    gc_stagg_x = compute_guard_cells(ex, stencils['x']['stagg'])
    gc_stagg_y = compute_guard_cells(ey, stencils['y']['stagg'])
    gc_stagg_z = compute_guard_cells(ez, stencils['z']['stagg'])

    print('\n    staggered solver:')
    print('    - {:d} guard cells along x'.format(gc_stagg_x))
    print('    - {:d} guard cells along y'.format(gc_stagg_y))
    print('    - {:d} guard cells along z'.format(gc_stagg_z))

    # Plot stencils
    plot_stencil(cx, stencils['x']['nodal'], stencils['x']['stagg'],
                 'x', path, name)
    plot_stencil(cy, stencils['y']['nodal'], stencils['y']['stagg'],
                 'y', path, name)
    plot_stencil(cz, stencils['z']['nodal'], stencils['z']['stagg'],
                 'z', path, name)

    print('\n >> Output Summary')
    print('    --------------')

    print('\n    path to figures: \'' + os.path.abspath(path) + '\'')

    print('\n    stencil arrays:')
    print('    - sx_nodal')
    print('    - sx_stagg')
    print('    - sy_nodal')
    print('    - sy_stagg')
    print('    - sz_nodal')
    print('    - sz_stagg\n')

    return stencils

if __name__ == '__main__':

    # --
    # User can modify these input parameters
    # --
    # Cell size
    dx = 1e-06
    dy = 1e-06
    dz = 2e-06
    # Time step
    dt = 1e-14
    # Spectral order
    nox = 8
    noy = 8
    noz = 16
    # Lorentz boost
    gamma = 30.
    # Galilean flag
    galilean = True
    # Error threshold
    # (might need to be smaller in z than (x,y) w/ Galilean algorithms)
    ex = 1e-07
    ey = 1e-07
    ez = 1e-07
    # Output path
    path = '.'
    # Output name tag (can be empty: '')
    name = 'test'
    # --

    # Run main function (some arguments are optional,
    # see definition of run_main function for help)
    stencils = run_main(dx, dy, dz, dt, nox, noy, noz, gamma,
                        galilean, ex, ey, ez, path, name)

    # Make stencil arrays available for inspection
    sx_nodal = stencils['x']['nodal']
    sx_stagg = stencils['x']['stagg']
    sy_nodal = stencils['y']['nodal']
    sy_stagg = stencils['y']['stagg']
    sz_nodal = stencils['z']['nodal']
    sz_stagg = stencils['z']['stagg']
