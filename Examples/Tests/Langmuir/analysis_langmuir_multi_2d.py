#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This is a script that analyzes the simulation results of the following CI tests:
# - Langmuir_multi_2d_nodal
# - Langmuir_multi_2d_psatd
# - Langmuir_multi_2d_psatd_nodal
# - Langmuir_multi_2d_psatd_momentum_conserving
# - Langmuir_multi_2d_psatd_current_correction
# - Langmuir_multi_2d_psatd_current_correction_nodal
# - Langmuir_multi_2d_psatd_Vay_deposition
# - Langmuir_multi_2d_psatd_Vay_deposition_nodal
# - Langmuir_multi_2d_MR
# The tests simulate a 3D periodic plasma wave. The theoretical electric field is given by
# E_x = epsilon \frac{me c2 k_x}{q_e} sin(kx x) cos(ky y) cos(kz z) sin(omega_p t)
# E_y = epsilon \frac{me c2 k_y}{q_e} cos(kx x) sin(ky y) cos(kz z) sin(omega_p t)
# E_z = epsilon \frac{me c2 k_z}{q_e} cos(kx x) cos(ky y) sin(kz z) sin(omega_p t)

import sys
import re
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
from scipy.constants import e, m_e, epsilon_0, c
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Name of the plotfile
fn = sys.argv[1]

# Parse test name and check if current correction (psatd.current_correction=1) is applied
current_correction = True if re.search( 'current_correction', fn ) else False

# Parse test name and check if Vay current deposition (algo.current_deposition=vay) is used
vay_deposition = True if re.search( 'Vay_deposition', fn ) else False

# Parameters (these parameters must match the parameters in the relevant input files as well as the
# additional parameters specified as runtime_params for the relevant tests in Regression/WarpX-tests.ini)
epsilon = 0.01
n = 4.e24
n_osc_x = 2
n_osc_z = 2
xmin = -20e-6; xmax = 20.e-6; Nx = 128
zmin = -20e-6; zmax = 20.e-6; Nz = 128

Lx = xmax - xmin
Lz = zmax - zmin

dx = Lx / Nx
dz = Lz / Nz

xarr = xmin + dx * (0.5 + np.arange(Nx))
zarr = zmin + dz * (0.5 + np.arange(Nz))

# Wave vector of the wave
kx = 2. * np.pi * n_osc_x / Lx
kz = 2. * np.pi * n_osc_z / Lz
# Plasma frequency
wp = np.sqrt((n * e**2) / (m_e * epsilon_0))

k = {'Ex': kx, 'Ez': kz}
cos = {'Ex': (0,1,1), 'Ez': (1,1,0)}

def get_theoretical_field(field, t):
    amplitude = (epsilon * m_e * c**2 * k[field] / e) * np.sin(wp*t)
    cos_flag = cos[field]
    x_contribution = np.cos(kx * xarr) if cos_flag[0] else np.sin(kx * xarr)
    z_contribution = np.cos(kz * zarr) if cos_flag[2] else np.sin(kz * zarr)
    E = amplitude * x_contribution[:,np.newaxis] * z_contribution[np.newaxis,:]
    return E

# Read the file and load data from main grid
ds = yt.load(fn)
t  = ds.current_time.to_value()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
edge = np.array([np.array([(ds.domain_left_edge[1]).item(), (ds.domain_right_edge[1]).item(), \
                           (ds.domain_left_edge[0]).item(), (ds.domain_right_edge[0]).item()])])

# Check the validity of the fields
error_rel = 0.
for field in ['Ex', 'Ez']:
    E_sim = data[field].to_ndarray()[:,:,0]
    E_th = get_theoretical_field(field, t)
    max_error = abs(E_sim - E_th).max() / abs(E_th).max()
    print('{:s}: max error = {}'.format(field, max_error))
    error_rel = max(error_rel, max_error)

tolerance = 0.05
print("error_rel = {}".format(error_rel))
print("tolerance = {}".format(tolerance))
assert(error_rel < tolerance)

# Plot Ez
E_sim = data['Ez'].to_ndarray()[:,:,0]

fig, (ax1,ax2) = plt.subplots(1,2, dpi=100)
# Absolute min and max field values
vmin = min(E_sim.min(), E_th.min())
vmax = max(E_sim.max(), E_th.max())
# First plot
cax1 = make_axes_locatable(ax1).append_axes('right', size='5%', pad='5%')
im1 = ax1.imshow(E_sim, origin = 'lower', extent = edge[0], vmin = vmin, vmax = vmax)
cb1 = fig.colorbar(im1, cax = cax1)
ax1.set_xlabel(r'$z$')
ax1.set_ylabel(r'$x$')
ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
ax1.set_title('Ez (sim)')
cb1.ax.tick_params()
cb1.ax.yaxis.get_offset_text().set()
cb1.ax.yaxis.offsetText.set_ha('center')
cb1.ax.yaxis.offsetText.set_va('bottom')
cb1.formatter.set_powerlimits((0,0))
cb1.update_ticks()
# Second plot
cax2 = make_axes_locatable(ax2).append_axes('right', size='5%', pad='5%')
im2 = ax2.imshow(E_th, origin = 'lower', extent = edge[0], vmin = vmin, vmax = vmax)
cb2 = fig.colorbar(im2, cax = cax2)
ax2.set_xlabel(r'$z$')
ax2.set_ylabel(r'$x$')
ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
ax2.set_title('Ez (theory)')
cb2.ax.tick_params()
cb2.ax.yaxis.get_offset_text().set()
cb2.ax.yaxis.offsetText.set_ha('center')
cb2.ax.yaxis.offsetText.set_va('bottom')
cb2.formatter.set_powerlimits((0,0))
cb2.update_ticks()
# Save figure
fig.tight_layout()
fig.savefig('analysis_langmuir_multi_2d.png', dpi=200)

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E) when
# current correction (psatd.do_current_correction=1) is applied or when
# Vay current deposition (algo.current_deposition=vay) is used
if current_correction or vay_deposition:
    rho  = data['rho' ].to_ndarray()
    divE = data['divE'].to_ndarray()
    error_rel = np.amax(np.abs(divE - rho/epsilon_0)) / np.amax(np.abs(rho/epsilon_0))
    tolerance = 1.e-9
    print("Check charge conservation:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert(error_rel < tolerance)

# Regression analysis on checksum benchmark
test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
