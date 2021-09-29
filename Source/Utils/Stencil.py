import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
from scipy.constants import c

plt.style.use('tableau-colorblind10')
plt.rcParams.update({'font.size': 14})
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

def get_stencil_coeffs(order, staggered):
 
    m = order // 2
    
    # Calculate the stencil coefficients by recurrence
    stencil_coef = np.zeros(m+1)
    
    # Calculation for a staggered scheme
    if staggered == True:
        # Initial coefficient
        prod = 1.
        for k in range(1, m+1):
            prod = prod * (m+k) / (4*k)
        stencil_coef[0] = 4 * m * prod**2
        # Get the other coefficients by recurrence
        for n in range(1, m+1):
            stencil_coef[n] = - ((2*n-3) * (m+1-n)) / ((2*n-1) * (m-1+n)) * stencil_coef[n-1]

    # Calculation for a centered scheme
    elif staggered == False:
        # Initial coefficient
        stencil_coef[0] = -2.
        # Get the other coefficients by recurrence
        for n in range(1, m+1):
            stencil_coef[n] = - (m+1-n) / (m+n) * stencil_coef[n-1]
        
    return(stencil_coef)

def modified_k(kx, dx, order, staggered):
 
    m = order // 2
    
    # Get the stencil coefficients
    stencil_coef = get_stencil_coeffs(order, staggered)
        
    # Array of values for n: from 1 to m
    n_array = np.arange(1, m+1)
    
    # Array of values of sin (first axis corresponds to k and second axis to n)
    if staggered:
        sin_array = np.sin(kx[:,np.newaxis] * (n_array[np.newaxis,:]-0.5) * dx) \
                  / ((n_array[np.newaxis,:]-0.5) * dx)
    else:
        sin_array = np.sin(kx[:,np.newaxis] * n_array[np.newaxis,:] * dx) \
                  / (n_array[np.newaxis,:] * dx)
            
    # Modified k
    k_array = np.tensordot(sin_array, stencil_coef[1:], axes=(-1,-1))
    
    return(k_array)

def func_cosine(om_s, w_c, dt):
    theta_c = np.exp(1.j * w_c * dt * 0.5)
    f = theta_c**2 * np.cos(om_s * dt)
    return f

def compute_stencils(coeff_nodal, coeff_stagg, dim, axis):
    stencil_nodal = np.fft.ifft(coeff_nodal, axis = axis)
    stencil_stagg = np.fft.ifft(coeff_stagg, axis = axis)
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

def compute_all(dx, dy, dz, dt, nox, noy, noz, Nx, Ny, Nz):

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
    w_c = vgal * kz_c
    
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
    parser.add_argument('--path', type = str, required = False, default = './',
                        help = 'path where figures are saved')
    parser.add_argument('--name', type = str, required = False, default = '',
                        help = 'common label for figure names')
    args = parser.parse_args()

    if (args.dy != 0 and args.noy == 0):
        sys.exit('You entered --dy, but not --noy. Please enter either both or none.')

    if (args.dy == 0 and args.noy != 0):
        sys.exit('You entered --noy, but not --dy. Please enter either both or none.')

    dim = 2 if (args.dy == 0 and args.noy == 0) else 3

    # Galilean velocity (default = 0.)
    vgal = - np.sqrt(1. - 1./args.gamma**2) * c
    
    # Display some output
    print('\n >> Numerical Setup:')
    print('\n    dx = {:g}'.format(args.dx))
    if (dim == 3):
        print('    dy = {:g}'.format(args.dy))
    print('    dz = {:g}'.format(args.dz))
    print('    dt = {:g}'.format(args.dt))
    print('    nox = {:d}'.format(args.nox))
    if (dim == 3):
        print('    noy = {:d}'.format(args.noy))
    print('    noz = {:d}'.format(args.noz))
    print('    gamma = {:g}'.format(args.gamma))

    # Number of grid points (arbitrary mesh)
    Nx = 256
    Ny = 256 if (dim == 3) else 0
    Nz = 256

    print('\n >> Compute Stencil:')

    stencil_nodal, stencil_stagg = compute_all(args.dx, args.dy, args.dz, args.dt,
                                               args.nox, args.noy, args.noz, Nx, Ny, Nz)

    stencil_x_nodal = stencil_nodal[0]
    if (dim == 2):
        stencil_z_nodal = stencil_nodal[1]
    elif (dim == 3):
        stencil_y_nodal = stencil_nodal[1]
        stencil_z_nodal = stencil_nodal[2]
    
    stencil_x_stagg = stencil_stagg[0]
    if (dim == 2):
        stencil_z_stagg = stencil_stagg[1]
    elif (dim == 3):
        stencil_y_stagg = stencil_stagg[1]
        stencil_z_stagg = stencil_stagg[2]
 
    # Maximum number of cells
    nx = 65
    if (dim == 3):
        ny = 65
    nz = 65

    # Array of cell numbers
    cx = np.arange(nx)
    if (dim == 3):
        cy = np.arange(ny)
    cz = np.arange(nz)

    # Arrays of stencils
    sx_nodal = abs(stencil_x_nodal[:nx])
    sx_stagg = abs(stencil_x_stagg[:nx])
    if (dim == 3):
        sy_nodal = abs(stencil_y_nodal[:ny])
        sy_stagg = abs(stencil_y_stagg[:ny])
    sz_nodal = abs(stencil_z_nodal[:nz])
    sz_stagg = abs(stencil_z_stagg[:nz])

    # Compute number of guard cells for given error threshold
    # (number of guard cells such that stencil measure is not smaller than the error threshold)
    ex = 1e-7
    if (dim == 3):
        ey = 1e-7
    ez = 1e-10

    # Nodal solver
    print('\n    nodal solver:')
    ddx = sx_nodal - ex
    mmx = ddx[ddx > 0.].min()
    gcx = np.argwhere(ddx == mmx)
    print('    - {:d} guard cells along x'.format(gcx[0,0]))
    if (dim == 3):
        ddy = sy_nodal - ey
        mmy = ddy[ddy > 0.].min()
        gcy = np.argwhere(ddy == mmy)
        print('    - {:d} guard cells along y'.format(gcy[0,0]))
    ddz = sz_nodal - ez
    mmz = ddz[ddz > 0.].min()
    gcz = np.argwhere(ddz == mmz)
    print('    - {:d} guard cells along z'.format(gcz[0,0]))

    # Staggered solver
    print('\n    staggered solver:')
    ddx = sx_stagg - ex
    mmx = ddx[ddx > 0.].min()
    gcx = np.argwhere(ddx == mmx)
    print('    - {:d} guard cells along x'.format(gcx[0,0]))
    if (dim == 3):
        ddy = sy_stagg - ey
        mmy = ddy[ddy > 0.].min()
        gcy = np.argwhere(ddy == mmy)
        print('    - {:d} guard cells along y'.format(gcy[0,0]))
    ddz = sz_stagg - ez
    mmz = ddz[ddz > 0.].min()
    gcz = np.argwhere(ddz == mmz)
    print('    - {:d} guard cells along z'.format(gcz[0,0]))

    # Plot stencil along x
    fig = plt.figure(dpi = 100)
    ax = fig.add_subplot(111)
    ax.plot(cx, sx_nodal, '-', color = cc[0], label = r'nodal')
    ax.plot(cx, sx_stagg, '-', color = cc[1], label = r'staggered, hybrid')
    ax.set_yscale('log')
    ax.set_xticks(cx, minor = True)
    ax.grid(which = 'minor', linewidth = 0.2)
    ax.grid(which = 'major', linewidth = 0.4)
    ax.legend()
    ax.set_xlabel('number of cells')
    ax.set_ylabel(r'$\Gamma(x)$')
    ax.set_title(r'Stencil extent along x')
    fig.tight_layout()
    fig_name = args.path + './figure_stencil_x'
    if (args.name):
        fig_name += '_' + args.name
    fig.savefig(fig_name + '.pdf', dpi = 100)
    fig.savefig(fig_name + '.png', dpi = 100)
    
    # Plot stencil along y
    if (dim == 3):
        fig = plt.figure(dpi = 100)
        ax = fig.add_subplot(111)
        ax.plot(cy, sy_nodal, '-', color = cc[0], label = r'nodal')
        ax.plot(cy, sy_stagg, '-', color = cc[1], label = r'staggered, hybrid')
        ax.set_yscale('log')
        ax.set_xticks(cx, minor = True)
        ax.grid(which = 'minor', linewidth = 0.2)
        ax.grid(which = 'major', linewidth = 0.4)
        ax.legend()
        ax.set_xlabel('number of cells')
        ax.set_ylabel(r'$\Gamma(y)$')
        ax.set_title(r'Stencil extent along y')
        fig.tight_layout()
        fig_name = args.path + './figure_stencil_y'
        if (args.name):
            fig_name += '_' + args.name
        fig.savefig(fig_name + '.pdf', dpi = 100)
        fig.savefig(fig_name + '.png', dpi = 100)
    
    # Plot stencil along z
    fig = plt.figure(dpi = 100)
    ax = fig.add_subplot(111)
    ax.plot(cz, sz_nodal, label = 'nodal')
    ax.plot(cz, sz_stagg, label = 'staggered, hybrid')
    ax.set_yscale('log')
    ax.set_xticks(cx, minor = True)
    ax.grid(which = 'minor', linewidth = 0.2)
    ax.grid(which = 'major', linewidth = 0.4)
    ax.legend()
    ax.set_xlabel('number of cells')
    ax.set_ylabel(r'$\Gamma(z)$')
    ax.set_title(r'Stencil extent along z')
    fig.tight_layout()
    fig_name = args.path + './figure_stencil_z'
    if (args.name):
        fig_name += '_' + args.name
    fig.savefig(fig_name + '.pdf', dpi = 100)
    fig.savefig(fig_name + '.png', dpi = 100)
    
    print('\n >> Path to Figures: \'' + os.path.abspath(args.path) + '\'\n')
