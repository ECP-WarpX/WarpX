import os

from matplotlib import cm, rcParams, use
import matplotlib.colors
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import openpmd_api as io
from scipy.constants import (c, centi, e, eV, femto, m_e, micron, milli, pi,
                             pico)

use('AGG')
rcParams.update({'font.size': 8})
MeV=1e6*eV

# create plot directory
plot_dir = './plots'
if os.path.exists(plot_dir) is False:
    os.mkdir(plot_dir)

# warpx
warpx_dir='./diags/particles/'
series = io.Series(warpx_dir+"/openpmd_%T.bp",io.Access.read_only)

# colormaps

nodes = [0 , 0.5, 1.]
colors = ["tomato", "white", "dodgerblue"]
cmap_rho = matplotlib.colors.LinearSegmentedColormap.from_list("cmap_rho", list(zip(nodes, colors)))

nodes = [0 , 1.]
colors =  ["white", "tomato"]
cmap_ele = matplotlib.colors.LinearSegmentedColormap.from_list("cmap_ele", list(zip(nodes, colors)))

nodes = [0 , 1.]
colors = ["white", "dodgerblue"]
cmap_pos = matplotlib.colors.LinearSegmentedColormap.from_list("cmap_pos", list(zip(nodes, colors)))

nodes = [0 , 1.]
colors = ["black", "white"]
cmap_pho = matplotlib.colors.LinearSegmentedColormap.from_list("cmap_pos", list(zip(nodes, colors)))


cmapE = matplotlib.colors.LinearSegmentedColormap.from_list("", ["orange","white","lightseagreen"])

my_dpi = 300.

def plot_one_panel(ax,data,cmap,title):
    im=ax.imshow(data, cmap=cmap, aspect='equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.001)
    fig.colorbar(im, cax=cax, orientation='vertical')
    ax.set_title(title)


for n,ts in enumerate(np.asarray(series.iterations)):
    print("timestep = %d" % ts)

    # plot images
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(1000./my_dpi, 1000./my_dpi), dpi=my_dpi)
    divider = make_axes_locatable(ax)

    j = series.iterations[ts]
    time = j.time

    '''
    rho_e = j.meshes["rho_beam_e"][io.Mesh_Record_Component.SCALAR]
    rho_p = j.meshes["rho_beam_p"][io.Mesh_Record_Component.SCALAR]
    rho_e_data = rho_e.load_chunk()
    rho_p_data = rho_p.load_chunk()
    series.flush()
    charge_density = rho_e_data+rho_p_data
    print(np.min(charge_density), np.max(charge_density))
    im=ax.imshow(charge_density, cmap=cmap_rho, aspect='equal', extent=[xmin,xmax,zmin,zmax], vmin=-20, vmax=20)
    cax = divider.append_axes('right', size='2%', pad=0.1)
    fig.colorbar(im, cax=cax, orientation='vertical', label=r'$\rho$ [m$^{-3}$]')
    E = j.meshes["E"]
    Ez = E["z"]
    data = Ez.load_chunk()
    series.flush()
    emax = data.max()
    emin = data.min()
    alphas = Normalize(0,emax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 0.8)
    colors = Normalize(emin, emax)(data)
    colors = cmapE(colors) #cm.seismic(colors)
    colors[..., -1] = alphas
    im=ax.imshow(colors, cmap=cmapE, aspect='equal', extent=[xmin,xmax,zmin,zmax], vmin=emin, vmax=emax)
    cax = divider.append_axes('right', size='2%', pad=0.7)
    fig.colorbar(im, cax=cax, orientation='vertical', label=r'$E_z [V/m]$')
    del charge_density, colors
    '''

    for ps in j.particles:
        print("\t {0}".format(ps))

    ps = j.particles["photon1"]
    x = ps["position"]["x"].load_chunk()
    z = ps["position"]["z"].load_chunk()
    px = ps["momentum"]["x"].load_chunk()
    py = ps["momentum"]["y"].load_chunk()
    pz = ps["momentum"]["z"].load_chunk()
    series.flush()
    ekin = c*np.sqrt(px**2+py**2+pz**2)/MeV
    index = np.argsort(ekin)

    if ekin.size:
        print('array is not empty')
        print(np.max(ekin))

        im = ax.scatter(x[index], z[index], s=1) #, c=ekin[index], cmap=cmap_pho, vmin=0, vmax=ekin.max(), zorder=12, alpha=0.8)
        #cax = divider.append_axes('right', size='2%', pad=0.75)
        #fig.colorbar(im, cax=cax, orientation='vertical', label=r'photon energy [MeV]')


    del x, z, px, py, pz, ekin


    ax.set_xlabel(r'x [mm]')
    ax.set_ylabel(r'y [mm]')
    ax.set_aspect('equal')

    ax.set_title("time = %.03f ps" %(time/pico));
    #plt.subplots_adjust(left=0.,bottom=0., right=0.5, top=1., hspace=0., wspace=0.)
    #__________________________________________________________
    # save image
    image_file_name =plot_dir+'/frame_%04d.png' % n
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=my_dpi)
    plt.close("all")

del series
