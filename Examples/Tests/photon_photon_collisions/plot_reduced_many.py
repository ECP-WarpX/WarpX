import os 
import sys 
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import micron, c, pi, centi, femto, e, milli, eV
import openpmd_api as io
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, Normalize
from matplotlib import use, cm
import matplotlib.colors

my_dpi=300
E0 = 0.511 
Lx = 4*milli
Lz = 2*milli
grid_x = np.linspace(-0.5*Lx, 0.5*Lx, 4096)
grid_z = np.linspace(-0.5*Lz, 0.5*Lz, 2048)
MeV=1e6*eV

fig, ax = plt.subplots(ncols=4, nrows=3, figsize=(4000./my_dpi, 3000./my_dpi), dpi=my_dpi, sharex=True)

lw = len(sys.argv[1:])+1 
for d in sys.argv[1:]:
    # warpx
    rdir= d+'diags/reducedfiles/'
    print(rdir)

    # number of photon1
    x = np.loadtxt(rdir+'ParticleNumber.txt')
    ax[0][0].plot(x[:,1], x[:,3], lw=lw)
    ax[0][0].set_title(r'macro photon1 number')

    # number of photon2    
    ax[0][1].plot(x[:,1], x[:,4], lw=lw)
    ax[0][1].set_title(r'macro photon2 number')
   
    print(np.max(x[:,5]),np.max(x[:,6]))

    # number of electrons 
    ax[0][2].plot(x[:,1], x[:,5], lw=lw)
    ax[0][2].set_title(r'macro electrons')

    # number of positrons  
    x = np.loadtxt(rdir+'ParticleNumber.txt')
    ax[0][3].plot(x[:,1], x[:,6], lw=lw)
    ax[0][3].set_title('macro positrons')


    
    
    lw = lw - 1 
    

for a in ax.reshape(-1):
    a.set_xlabel('time [s]')

image_file_name ='many_reduced.png'
plt.tight_layout()
plt.savefig(image_file_name,dpi=300, bbox_inches='tight') 
plt.close("all")
