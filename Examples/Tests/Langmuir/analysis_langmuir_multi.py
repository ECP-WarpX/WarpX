#! /usr/bin/env python

# Copyright 2019 Jean-Luc Vay, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This is a script that analyses the simulation results from
# the script `inputs.multi.rt`. This simulates a 3D periodic plasma wave.
# The electric field in the simulation is given (in theory) by:
# $$ E_x = \epsilon \,\frac{m_e c^2 k_x}{q_e}\sin(k_x x)\cos(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_y = \epsilon \,\frac{m_e c^2 k_y}{q_e}\cos(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
# $$ E_z = \epsilon \,\frac{m_e c^2 k_z}{q_e}\cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
yt.funcs.mylog.setLevel(50)
import numpy as np
from scipy.constants import e, m_e, epsilon_0, c

# this will be the name of the plot file
fn = sys.argv[1]

# Parse test name and check if current correction (psatd.do_current_correction=1) is applied
current_correction = True if re.search( 'current_correction', fn ) else False

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

k = {'Ex':kx, 'Ey':ky, 'Ez':kz}
cos = {'Ex': (0,1,1), 'Ey':(1,0,1), 'Ez':(1,1,0)}

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

# Read the file
ds = yt.load(fn)

# Check that the particle selective output worked:
species = 'electrons'
print('ds.field_list', ds.field_list)
for field in ['particle_weight',
              'particle_momentum_x']:
    print('assert that this is in ds.field_list', (species, field))
    assert (species, field) in ds.field_list
for field in ['particle_momentum_y',
              'particle_momentum_z']:
    print('assert that this is NOT in ds.field_list', (species, field))
    assert (species, field) not in ds.field_list
species = 'positrons'
for field in ['particle_Ey']:
    print('assert that this is in ds.field_list', (species, field))
    assert (species, field) in ds.field_list
for field in ['particle_momentum_y',
              'particle_momentum_z']:
    print('assert that this is NOT in ds.field_list', (species, field))
    assert (species, field) not in ds.field_list

t0 = ds.current_time.to_ndarray().mean()
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                                    dims=ds.domain_dimensions)

# Check the validity of the fields
error_rel = 0
for field in ['Ex', 'Ey', 'Ez']:
    E_sim = data[field].to_ndarray()

    E_th = get_theoretical_field(field, t0)
    max_error = abs(E_sim-E_th).max()/abs(E_th).max()
    print('%s: Max error: %.2e' %(field,max_error))
    error_rel = max( error_rel, max_error )

# Plot the last field from the loop (Ez at iteration 40)
plt.subplot2grid( (1,2), (0,0) )
plt.imshow( E_sim[:,:,Ncell[2]//2] )
plt.colorbar()
plt.title('Ez, last iteration\n(simulation)')
plt.subplot2grid( (1,2), (0,1) )
plt.imshow( E_th[:,:,Ncell[2]//2] )
plt.colorbar()
plt.title('Ez, last iteration\n(theory)')
plt.tight_layout()
plt.savefig('langmuir_multi_analysis.png')

tolerance_rel = 0.15

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

# Check relative L-infinity spatial norm of rho/epsilon_0 - div(E) when
# current correction (psatd.do_current_correction=1) is applied
if current_correction:
    rho  = data['rho' ].to_ndarray()
    divE = data['divE'].to_ndarray()
    Linf_norm = np.amax( np.abs( rho/epsilon_0 - divE ) ) / np.amax( np.abs( rho/epsilon_0 ) )
    print("error: " + str(Linf_norm))
    print("tolerance: 1.e-9")
    assert( Linf_norm < 1.e-9 )
