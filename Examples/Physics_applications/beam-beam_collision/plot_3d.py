#!/usr/bin/env python3

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Arianna Formenti
# License: BSD-3-Clause-LBNL
#
# This script plots some diagnostics of a beam-beam collision.

import sys
import os 

import matplotlib.pyplot as plt
import openpmd_api as io
import numpy as np 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.constants import e as q_e
import pandas as pd 

my_dpi=300.

def plot_bbc(): 
    # this will be the name of the plot file
    #fn = sys.argv[1]

    '''
    # Read the file
    ds = yt.load(fn)

    # plot the laser field and absolute density
    fields = ["Ey", "rho"]
    normal = "y"
    sl = yt.SlicePlot(ds, normal=normal, fields=fields)
    for field in fields:
        sl.set_log(field, False)

    sl.set_figure_size((4, 8))
    fig = sl.export_to_mpl_figure(nrows_ncols=(2, 1))
    fig.tight_layout()
    plt.show()
    '''
    
    path=os.path.join(sys.argv[1],'diags/field_zx')
    series = io.Series(path+"/openpmd_%T.bp", io.Access.read_only)
    iterations = np.asarray(series.iterations)


    for n,ts in enumerate(iterations[0:3]):
        it = series.iterations[ts]
        print(f"Iteration {ts} contains {len(it.particles)} particle species")


    for n,ts in enumerate(iterations):

        fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(1000./my_dpi, 2000./my_dpi), dpi=my_dpi, sharex=True, sharey=True)
         
        it = series.iterations[ts]
        fig.suptitle(f"step = {ts}, t = {it.time} s")

        print(ts)
        print(f"Iteration {ts} contains {len(it.particles)} particle species")
        
        E = it.meshes["E"]
        Ex = E["x"].load_chunk()        
        rho1 = it.meshes["rho_beam1"][io.Mesh_Record_Component.SCALAR].load_chunk()
        rho2 = it.meshes["rho_beam2"][io.Mesh_Record_Component.SCALAR].load_chunk()

        series.flush()
        
        dz, dy, dx = E.grid_spacing
        nz, ny, nx = Ex.shape
        z0, y0, x0 = E.grid_global_offset

        extent = [z0, (z0+nz*dz), x0, (x0+nx*dx)]

        print(np.shape(Ex))

        Ex = 0.5*(Ex[:,0,:]+Ex[:,1,:])
        rho1 = 0.5*(rho1[:,0,:]+rho1[:,1,:])
        rho2 = 0.5*(rho2[:,0,:]+rho2[:,1,:])
        
        im=ax[0].imshow(np.transpose(Ex), cmap='RdGy', extent=extent) #, vmin=-5e13, vmax=5e13)
        ax[0].set_title(r'E$_x$')
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        
        
        im=ax[1].imshow(np.transpose(rho1+rho2), cmap='rainbow', extent=extent) #, vmin=-2e11, vmax=0)
        ax[1].set_title(r'$n_1+n_2$ beams')
        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        
        it.close()
        
        '''
        it2 = series2.iterations[ts]
        print(f"Iteration {ts} contains {len(it2.particles)} particle species")
        ps = it2.particles["beam_e"]
        x = ps["position"]["x"].load_chunk()
        y = ps["position"]["y"].load_chunk()
        #z = ps["position"]["z"].load_chunk()
        #px = ps["momentum"]["x"].load_chunk()
        #py = ps["momentum"]["y"].load_chunk()
        #pz = ps["momentum"]["z"].load_chunk()
        id = ps["id"][io.Mesh_Record_Component.SCALAR].load_chunk()
        #w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk() 
        #gamma = np.sqrt(1.+(px**2+py**2+pz**2)/(m_e*c)**2)
        #E = (gamma-1.)*m_e*c**2/GeV
        
        series2.flush()

        index = np.argsort(id)
        
        ax[3][0].scatter(x[index],y[index],c=id, s=0.2, zorder=11, cmap='rainbow')


        ps = it2.particles["beam_p"]
        x = ps["position"]["x"].load_chunk()
        y = ps["position"]["y"].load_chunk()
        id = ps["id"][io.Mesh_Record_Component.SCALAR].load_chunk()

        series2.flush()

        index = np.argsort(id)
        
        ax[3][1].scatter(x[index],y[index],c=id, s=0.2, zorder=11, cmap='rainbow')

        ps = it2.particles["pho"]
        x = ps["position"]["x"].load_chunk()
        y = ps["position"]["y"].load_chunk()
        id = ps["id"][io.Mesh_Record_Component.SCALAR].load_chunk()
        
        series2.flush()
        it2.close()

        index = np.argsort(id)
        ax[3][2].scatter(x[index],y[index], c=id, s=0.2, zorder=11, cmap='rainbow')
        '''
        

        for a in ax.reshape(-1):
            a.set_xlabel(r'z [$\mu$m]')
            a.set_ylabel(r'x [nm]')
            a.set_aspect('equal')
            

        image_file_name =os.path.join(sys.argv[1], 'plots/full_%04d.png' % n) 
        plt.tight_layout()
        plt.savefig(image_file_name,dpi=300, bbox_inches='tight') 
        plt.close("all")


    del series
        
def check_energy_conservation():  

    species_list = ['beam1', 'beam2', 'pho1', 'pho2', 'ele_nlbw1', 'pos_nlbw1', 'ele_nlbw2', 'pos_nlbw2', 'pho', 'ele', 'pos']

    kinetic_energy = np.zeros((57,))

    df = pd.read_csv('diags/reducedfiles/ParticleEnergy.txt', sep=" ", header=0)
    for k,s in enumerate(species_list):
        kinetic_energy += df[[col for col in df.columns if f']'+s+'(J)' in col]].to_numpy()[:,0]
             
    df = pd.read_csv('diags/reducedfiles/FieldEnergy.txt', sep=" ", header=0)
    field_energy = df[[col for col in df.columns if f']total_lev0(J)' in col]].to_numpy()[:,0]
    
    total_energy = kinetic_energy + field_energy
    

    plt.plot(total_energy)
    plt.show()
    #times = df_cr[[col for col in df_cr.columns if f']time' in col]].to_numpy()
    #steps = df_cr[[col for col in df_cr.columns if f']step' in col]].to_numpy()



#def plot_luminosity():
#    df_cr = pd.read_csv(rdir+'ColliderRelevant_beam1_beam2.txt', sep=" ", header=0)
# luminosity
    #x = df_cr[[col for col in df_cr.columns if f']dL_dt' in col]].to_numpy()
    #coll_index = np.argmax(x)
    #coll_time = times[coll_index]
    #print(coll_index)
    #L = [np.trapz(x[:k,0], times[:k,0]) for k in range(len(times))]
    #ax[0][0].plot(times, np.asarray(L)/1e4, lw=lw)
    #ax[0][0].set_title('L [cm$^{-2}$]')
    #ax[0][0].axhline(y=L_geom, lw=3, color='mediumseagreen', label='L geom')
    #ax[0][0].legend()


    
    

if __name__ == "__main__":
    plot_bbc()
    #check_energy_conservation()
