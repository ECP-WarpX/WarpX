#!/usr/bin/env python3

import os
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import yt
import scipy.constants as scc
from openpmd_viewer import OpenPMDTimeSeries

import yt ; yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI
filename = sys.argv[1]

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)

#Visualization of the Ey field
#ts = OpenPMDTimeSeries('./diags/diag1')
#Ey_2d, info_Ey_2d = ts.get_field( field='E', coord='y', iteration=40, slice_across=None )
#print( Ey_2d )
#plt.figure()
#matplotlib widget
#ts.slider()



#Visualization of the Ey field
yt.funcs.mylog.setLevel(50)

def get_contribution( is_cos, k ):
    du = (xmax-xmin)/Nx
    u = xmin + du*( 0.5 + np.arange(Nx) )
    if is_cos == 1:
        return( np.cos(k*u) )
    else:
        return( np.sin(k*u) )

def get_theoretical_field( field, t ):
    amplitude = epsilon * (m_e*c**2*k[field])/e * np.sin(wp*t)
    cos_flag = cos[field]
    x_contribution = get_contribution( cos_flag[0], kx )
    z_contribution = get_contribution( cos_flag[2], kz )

    E = amplitude * x_contribution[:, np.newaxis ] \
                  * z_contribution[np.newaxis, :]

    return( E )

for field in ['Ex', 'Ez']:
    E_sim = data[('mesh',field)].to_ndarray()[:,:,0]
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim-E_th).max()/abs(E_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    error_rel = max( error_rel, max_error )


# Read the file
ds = yt.load('./diags/diag1000050/')
t0 = ds.current_time.to_value()
data = ds.covering_grid(level = 0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)
edge = np.array([(ds.domain_left_edge[1]).item(), (ds.domain_right_edge[1]).item(), \
                 (ds.domain_left_edge[0]).item(), (ds.domain_right_edge[0]).item()])

vmin = E_sim.min()
vmax = E_sim.max()
cax = make_axes_locatable(ax1).append_axes('right', size = '5%', pad = '5%')
im = ax1.imshow(E_sim, origin = 'lower', extent = edge, vmin = vmin, vmax = vmax)
cb = fig.colorbar(im1, cax = cax)
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$x$')
ax.set_title(r'$E_z$ (sim)')

# Save figure
fig.tight_layout()
fig.savefig('laser_envelope_model.png', dpi = 200)