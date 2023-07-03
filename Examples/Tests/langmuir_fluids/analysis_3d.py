#!/usr/bin/env python3

# Copyright 2019-2022 Jean-Luc Vay, Maxence Thevenet, Remi Lehe, Axel Huebl
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = \epsilon \,\frac{m_e c^2 k_x}{q_e}\sin(k_x x)\cos(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_y = \epsilon \,\frac{m_e c^2 k_y}{q_e}\cos(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
import os
import re
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import yt

yt.funcs.mylog.setLevel(50)

import numpy as np
from scipy.constants import c, e, epsilon_0, m_e

sys.path.insert(1, '/Users/gjohnson/src/warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Parse test name and check if current correction (psatd.current_correction=1) is applied
current_correction = True if re.search( 'current_correction', fn ) else False

# Parse test name and check if Vay current deposition (algo.current_deposition=vay) is used
vay_deposition = True if re.search( 'Vay_deposition', fn ) else False

# Parse test name and check if div(E)/div(B) cleaning (warpx.do_div<e,b>_cleaning=1) is used
div_cleaning = True if re.search('div_cleaning', fn) else False

# Parameters (these parameters must match the parameters in `inputs.multi.rt`)
epsilon = 0.01
n = 4.e24
n_osc_x = 2
n_osc_y = 2
n_osc_z = 2
lo = [-20.e-6, -20.e-6, -20.e-6]
hi = [ 20.e-6,  20.e-6,  20.e-6]
Ncell = [64, 64, 64]

# Wave vector of the wave
kx = 2.*np.pi*n_osc_x/(hi[0]-lo[0])
ky = 2.*np.pi*n_osc_y/(hi[1]-lo[1])
kz = 2.*np.pi*n_osc_z/(hi[2]-lo[2])
# Plasma frequency
wp = np.sqrt((n*e**2)/(m_e*epsilon_0))

k = {'Ex':kx, 'Ey':ky, 'Ez':kz, 'Jx':kx, 'Jy':ky, 'Jz':kz}
cos = {'Ex': (0,1,1), 'Ey':(1,0,1), 'Ez':(1,1,0),'Jx': (0,1,1), 'Jy':(1,0,1), 'Jz':(1,1,0)}
cos_rho = {'rho': (1,1,1)}

def get_contribution( is_cos, k, idim ):
    du = (hi[idim]-lo[idim])/Ncell[idim]
    u = lo[idim] + du*( 0.5 + np.arange(Ncell[idim]) )
    if is_cos[idim] == 1:
        return( np.cos(k*u) )
    else:
        return( np.sin(k*u) )

def get_theoretical_field( field, t ):
    amplitude = epsilon * (m_e*c**2*k[field])/e * np.sin(wp*t)
    cos_flag = cos[field]
    x_contribution = get_contribution( cos_flag, kx, 0 )
    y_contribution = get_contribution( cos_flag, ky, 1 )
    z_contribution = get_contribution( cos_flag, kz, 2 )

    E = amplitude * x_contribution[:, np.newaxis, np.newaxis] \
                  * y_contribution[np.newaxis, :, np.newaxis] \
                  * z_contribution[np.newaxis, np.newaxis, :]

    return( E )

def get_theoretical_J_field( field, t ):
    amplitude = - epsilon_0 * wp * epsilon * (m_e*c**2*k[field])/e * np.cos(wp*t)
    cos_flag = cos[field]
    x_contribution = get_contribution( cos_flag, kx, 0 )
    y_contribution = get_contribution( cos_flag, ky, 1 )
    z_contribution = get_contribution( cos_flag, kz, 2 )

    J = amplitude * x_contribution[:, np.newaxis, np.newaxis] \
                  * y_contribution[np.newaxis, :, np.newaxis] \
                  * z_contribution[np.newaxis, np.newaxis, :]

    return( J )

def get_theoretical_rho_field( field, t ):
    amplitude = epsilon_0 * epsilon * (m_e*c**2*(kx*kx+ky*ky+kz*kz))/e * np.sin(wp*t)
    cos_flag = cos_rho[field]
    x_contribution = get_contribution( cos_flag, kx, 0 )
    y_contribution = get_contribution( cos_flag, ky, 1 )
    z_contribution = get_contribution( cos_flag, kz, 2 )

    rho = amplitude * x_contribution[:, np.newaxis, np.newaxis] \
                  * y_contribution[np.newaxis, :, np.newaxis] \
                  * z_contribution[np.newaxis, np.newaxis, :]

    return( rho )

# Read the file
ds = yt.load(fn)


t0 = ds.current_time.to_value()
data = ds.covering_grid(level = 0, left_edge = ds.domain_left_edge, dims = ds.domain_dimensions)
edge = np.array([(ds.domain_left_edge[2]).item(), (ds.domain_right_edge[2]).item(), \
                 (ds.domain_left_edge[0]).item(), (ds.domain_right_edge[0]).item()])

# Check the validity of the fields
error_rel = 0
for field in ['Ex', 'Ey', 'Ez']:
    E_sim = data[('mesh',field)].to_ndarray()
    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim-E_th).max()/abs(E_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    error_rel = max( error_rel, max_error )
    print("error_rel    : " + str(error_rel))


# Check the validity of the currents
error_rel = 0
for field in ['Jx', 'Jy', 'Jz']:
    J_sim = data[('mesh',field)].to_ndarray()
    J_th = get_theoretical_J_field(field, t0)
    max_error = abs(J_sim-J_th).max()/abs(J_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    error_rel = max( error_rel, max_error )
    print("error_rel    : " + str(error_rel))

# Check the validity of the charge
error_rel = 0
for field in ['rho']:
    rho_sim = data[('boxlib',field)].to_ndarray()
    rho_th = get_theoretical_rho_field(field, t0)
    max_error = abs(rho_sim-rho_th).max()/abs(rho_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    error_rel = max( error_rel, max_error )
    print("error_rel    : " + str(error_rel))


# Plot the last field from the loop (Ez at iteration 40)
fig, (ax1, ax2) = plt.subplots(1, 2, dpi = 100)
# First plot (slice at y=0)
E_plot = E_sim[:,Ncell[1]//2+1,:]
vmin = E_plot.min()
vmax = E_plot.max()
cax1 = make_axes_locatable(ax1).append_axes('right', size = '5%', pad = '5%')
im1 = ax1.imshow(E_plot, origin = 'lower', extent = edge, vmin = vmin, vmax = vmax)
cb1 = fig.colorbar(im1, cax = cax1)
ax1.set_xlabel(r'$z$')
ax1.set_ylabel(r'$x$')
ax1.set_title(r'$E_z$ (sim)')
# Second plot (slice at y=0)
E_plot = E_th[:,Ncell[1]//2+1,:]
vmin = E_plot.min()
vmax = E_plot.max()
cax2 = make_axes_locatable(ax2).append_axes('right', size = '5%', pad = '5%')
im2 = ax2.imshow(E_plot, origin = 'lower', extent = edge, vmin = vmin, vmax = vmax)
cb2 = fig.colorbar(im2, cax = cax2)
ax2.set_xlabel(r'$z$')
ax2.set_ylabel(r'$x$')
ax2.set_title(r'$E_z$ (theory)')
# Save figure
fig.tight_layout()
fig.savefig('Langmuir_fluid_multi_analysis.png', dpi = 200)

tolerance_rel = 5e-2

print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E)
# with current correction (and periodic single box option) or with Vay current deposition
if current_correction:
    tolerance = 1e-9
elif vay_deposition:
    tolerance = 1e-3
if current_correction or vay_deposition:
    rho  = data[('boxlib','rho')].to_ndarray()
    divE = data[('boxlib','divE')].to_ndarray()
    error_rel = np.amax( np.abs( divE - rho/epsilon_0 ) ) / np.amax( np.abs( rho/epsilon_0 ) )
    print("Check charge conservation:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert( error_rel < tolerance )

if div_cleaning:
    ds_old = yt.load('Langmuir_multi_psatd_div_cleaning_plt000038')
    ds_mid = yt.load('Langmuir_multi_psatd_div_cleaning_plt000039')
    ds_new = yt.load(fn) # this is the last plotfile

    ad_old = ds_old.covering_grid(level = 0, left_edge = ds_old.domain_left_edge, dims = ds_old.domain_dimensions)
    ad_mid = ds_mid.covering_grid(level = 0, left_edge = ds_mid.domain_left_edge, dims = ds_mid.domain_dimensions)
    ad_new = ds_new.covering_grid(level = 0, left_edge = ds_new.domain_left_edge, dims = ds_new.domain_dimensions)

    rho   = ad_mid['rho'].v.squeeze()
    divE  = ad_mid['divE'].v.squeeze()
    F_old = ad_old['F'].v.squeeze()
    F_new = ad_new['F'].v.squeeze()

    # Check max norm of error on dF/dt = div(E) - rho/epsilon_0
    # (the time interval between the old and new data is 2*dt)
    dt = 1.203645751e-15
    x = F_new - F_old
    y = (divE - rho/epsilon_0) * 2 * dt
    error_rel = np.amax(np.abs(x - y)) / np.amax(np.abs(y))
    tolerance = 1e-2
    print("Check div(E) cleaning:")
    print("error_rel = {}".format(error_rel))
    print("tolerance = {}".format(tolerance))
    assert(error_rel < tolerance)

test_name = os.path.split(os.getcwd())[1]

if re.search( 'single_precision', fn ):
    checksumAPI.evaluate_checksum(test_name, fn, rtol=1.e-3)
else:
    checksumAPI.evaluate_checksum(test_name, fn)
