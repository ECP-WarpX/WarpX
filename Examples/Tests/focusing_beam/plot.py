import os 
import sys 
import numpy as np
import openpmd_api as io
from scipy.constants import m_e, c, eV, micro, nano, milli, e as q_e
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, Normalize, SymLogNorm, ListedColormap, LinearSegmentedColormap
from matplotlib import use, cm, rcParams

mydir='.'
path=os.path.join(mydir, 'diags/full')
series = io.Series(path+"/openpmd_%T.bp", io.Access.read_only)
iterations = np.asarray(series.iterations)

my_dpi=300.


for n,ts in enumerate(iterations):

    fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(1500./my_dpi, 4500./my_dpi), dpi=my_dpi)
    it = series.iterations[ts]
    
    ps = it.particles["beam1"]
    x = ps["position"]["x"].load_chunk()
    y = ps["position"]["y"].load_chunk()
    z = ps["position"]["z"].load_chunk()
    w= ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()
    print(len(w))
    ax[0].scatter(z,x,color='blue')
    ax[1].scatter(z,y,color='blue')
    ax[2].scatter(x,y,color='blue')

    it.close()

    for a in ax.reshape(-1):
        a.set_aspect('auto')

    ax[0].set_xlabel(r'z')
    ax[0].set_ylabel(r'x')

    ax[1].set_xlabel(r'z')
    ax[1].set_ylabel(r'y')

    ax[2].set_xlabel(r'x')
    ax[2].set_ylabel(r'y')


    image_file_name =os.path.join('TEST_%04d.png' % n) 
    plt.tight_layout()
    plt.savefig(image_file_name,dpi=my_dpi, bbox_inches='tight') 
    plt.close("all")
