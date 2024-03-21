#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from openpmd_viewer import OpenPMDTimeSeries

series = OpenPMDTimeSeries('./diags/diag1')
steps = series.iterations
print(steps)

ylabel = 'y [m]'
xlabel = 'z [m]'
#slice_axis = 'x'
slice_axis = 'y'

#loop through E,B,Rho
for n in steps:

    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(40, 10), dpi=100., sharex=True, sharey=True)

    #E field
    Ex, info = series.get_field(field='E', coord='x', iteration=n, plot=False, slice_across=slice_axis)
    Ey, info = series.get_field(field='E', coord='y', iteration=n, plot=False, slice_across=slice_axis)
    Ez, info = series.get_field(field='E', coord='z', iteration=n, plot=False, slice_across=slice_axis)

    #B field
    Bx, info = series.get_field(field='B', coord=slice_axis, iteration=n, plot=False, slice_across=slice_axis)
    By, info = series.get_field(field='B', coord='y', iteration=n, plot=False, slice_across=slice_axis)
    Bz, info = series.get_field(field='B', coord='z', iteration=n, plot=False, slice_across=slice_axis)

    # Rho
    rho_beam1, info = series.get_field(field='rho_beam1', iteration=n, plot=False, slice_across=slice_axis)
    rho_beam2, info = series.get_field(field='rho_beam2', iteration=n, plot=False, slice_across=slice_axis)
    rho_ele1, info = series.get_field(field='rho_ele1', iteration=n, plot=False, slice_across=slice_axis)
    rho_pos1, info = series.get_field(field='rho_pos1', iteration=n, plot=False, slice_across=slice_axis)
    rho_ele2, info = series.get_field(field='rho_ele2', iteration=n, plot=False, slice_across=slice_axis)
    rho_pos2, info = series.get_field(field='rho_pos2', iteration=n, plot=False, slice_across=slice_axis)

    xmin = info.z.min()
    xmax = info.z.max()
    if slice_axis == 'x':
        ymin = info.y.min()
        ymax = info.y.max()
    elif slice_axis == 'y':
        ymin = info.x.min()
        ymax = info.x.max()


     #E field plots
    im1 = ax[0,0].imshow(np.transpose(Ex), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[0,0].set_title(f'E$_x$', fontsize=20)
    divider1 = make_axes_locatable(ax[0,0])
    cax = divider1.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')


    im2 = ax[0,1].imshow(np.transpose(Ey), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[0,1].set_title(f'E$_y$', fontsize=20)
    divider2 = make_axes_locatable(ax[0,1])
    cax = divider2.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')

    im3 = ax[0,2].imshow(np.transpose(Ez), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[0,2].set_title(f'E$_z$', fontsize=20)
    divider3 = make_axes_locatable(ax[0,2])
    cax = divider3.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3, cax=cax, orientation='vertical')

    #B field plots
    im4 = ax[1,0].imshow(np.transpose(Bx), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[1,0].set_title(f'B$_x$', fontsize=20)
    divider4 = make_axes_locatable(ax[1,0])
    cax = divider4.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im4, cax=cax, orientation='vertical')


    im5 = ax[1,1].imshow(np.transpose(By), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[1,1].set_title(f'B$_y$', fontsize=20)
    divider5 = make_axes_locatable(ax[1,1])
    cax = divider5.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im5, cax=cax, orientation='vertical')

    im6 = ax[1,2].imshow(np.transpose(Bz), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[1,2].set_title(f'B$_z$', fontsize=20)
    divider6 = make_axes_locatable(ax[1,2])
    cax = divider6.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im6, cax=cax, orientation='vertical')

    #Rho plots
    im7 = ax[0,3].imshow(np.transpose(rho_beam1+rho_beam2), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[0,3].set_title(f'rho$_{{\tbeams{{}}}}$', fontsize=20)
    divider7 = make_axes_locatable(ax[0,3])
    cax = divider7.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im7, cax=cax, orientation='vertical')


    im8 = ax[1,3].imshow(np.transpose(rho_ele1+rho_pos1+rho_ele2+rho_pos2), cmap='seismic', extent=[xmin , xmax, ymin, ymax])
    ax[1,3].set_title(f'rho$_{{\tsecondaries{{}}}}$', fontsize=20)
    divider8 = make_axes_locatable(ax[1,3])
    cax = divider8.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im8, cax=cax, orientation='vertical')


    ax[1,0].set_ylabel(ylabel, fontsize=20)
    ax[0,0].set_ylabel(ylabel, fontsize=20)
    ax[1,1].set_xlabel(xlabel, fontsize=20)
    ax[1,2].set_xlabel(xlabel, fontsize=20)
    ax[1,3].set_xlabel(xlabel, fontsize=20)
    ax[1,0].set_xlabel(xlabel, fontsize=20)


    fig.suptitle(f'Iteration {n:0}', fontsize=20)
    plt.tight_layout()

    image_file_name = 'FULL_'+slice_axis+f'{n:03d}.png'
    plt.savefig(image_file_name, dpi=100, bbox_inches='tight')
    plt.close()
