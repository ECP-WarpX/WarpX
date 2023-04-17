"""
Python script to compute the minimum number of guard cells for a given
error threshold, based on the measurement of the PSATD stencil extent
(that is, the minimum number of guard cells such that the stencil
measure is not larger than the error threshold).
Reference: https://doi.org/10.1016/j.cpc.2022.108457

Run the script simply with "python stencil.py" (or with "run stencil.py"
using IPython). The user can modify the input parameters set in the main
function at the end of the file.
"""

import os

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
    kx_arr = 2*np.pi*np.fft.fftfreq(nx, dx)
    ky_arr = 2*np.pi*np.fft.fftfreq(ny, dy)
    kz_arr = 2*np.pi*np.fft.fftfreq(nz, dz)

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
    kk_c = np.sqrt(kx_c**2+ky_c**2+kz_c**2)
    kk_s = np.sqrt(kx_s**2+ky_s**2+kz_s**2)
    om_c = c*kk_c
    om_s = c*kk_s
    w_c = v_gal*kz_c

    # Spectral coefficient
    coeff_nodal = func_cosine(om_c, w_c, dt)
    coeff_stagg = func_cosine(om_s, w_c, dt)

    # Stencils
    stencils = dict()
    stencils['x'] = compute_stencils(coeff_nodal, coeff_stagg, axis=0)
    stencils['y'] = compute_stencils(coeff_nodal, coeff_stagg, axis=1)
    stencils['z'] = compute_stencils(coeff_nodal, coeff_stagg, axis=2)

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

def run_main(dx, dy, dz, dt, nox, noy, noz, gamma=1., galilean=False, path='.', name=''):
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
    print('- dx = {:g}'.format(dx))
    print('- dy = {:g}'.format(dy))
    print('- dz = {:g}'.format(dz))
    print('- dz/dx = {:g}'.format(dz/dx))
    print('- dz/dy = {:g}'.format(dz/dy))
    print('\nTime step:')
    print('- dt = {:g}'.format(dt))
    print('- c*dt/dx = {:g}'.format(c*dt/dx))
    print('- c*dt/dy = {:g}'.format(c*dt/dy))
    print('- c*dt/dz = {:g}'.format(c*dt/dz))
    print('\nSpectral order:')
    print('- nox = {:d}'.format(nox))
    print('- noy = {:d}'.format(noy))
    print('- noz = {:d}'.format(noz))
    print('\nLorentz boost, Galilean velocity:')
    print('- gamma = {:g}'.format(gamma))
    print('- v_gal = {:g}'.format(v_gal))

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

    # Plot stencils
    plot_stencil(cx, stencils['x']['nodal'], stencils['x']['stagg'], 'x', path, name)
    plot_stencil(cy, stencils['y']['nodal'], stencils['y']['stagg'], 'y', path, name)
    plot_stencil(cz, stencils['z']['nodal'], stencils['z']['stagg'], 'z', path, name)

    # Compute min and max numbers of guard cells
    gcmin_x_nodal, gcmax_x_nodal = compute_guard_cells(sp, dp, stencils['x']['nodal'])
    gcmin_y_nodal, gcmax_y_nodal = compute_guard_cells(sp, dp, stencils['y']['nodal'])
    gcmin_z_nodal, gcmax_z_nodal = compute_guard_cells(sp, dp, stencils['z']['nodal'])
    #
    gcmin_x_stagg, gcmax_x_stagg = compute_guard_cells(sp, dp, stencils['x']['stagg'])
    gcmin_y_stagg, gcmax_y_stagg = compute_guard_cells(sp, dp, stencils['y']['stagg'])
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
    print(f'- between {gcmin_y_nodal} and {gcmax_y_nodal} ghost cells along y')
    print(f'- between {gcmin_z_nodal} and {gcmax_z_nodal} ghost cells along z')
    print('\nFor a staggered or hybrid simulation, choose:')
    print(f'- between {gcmin_x_stagg} and {gcmax_x_stagg} ghost cells along x')
    print(f'- between {gcmin_y_stagg} and {gcmax_y_stagg} ghost cells along y')
    print(f'- between {gcmin_z_stagg} and {gcmax_z_stagg} ghost cells along z')
    print()

    return

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
    # Output path
    path = '.'
    # Output name tag (can be empty: '')
    name = 'test'
    # --

    # Run main function (some arguments are optional,
    # see definition of run_main function for help)
    run_main(dx, dy, dz, dt, nox, noy, noz, gamma, galilean, path, name)
