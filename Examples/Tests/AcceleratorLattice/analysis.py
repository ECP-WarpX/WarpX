#!/usr/bin/env python3

# Copyright 2022 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


"""
This script tests the quad lattice
The input file sets up a series of quadrupoles and propagates two particles through them.
One particle is in the X plane, the other the Y plane.
The final positions are compared to the analytic solutions.
The motion is slow enough that relativistic effects are ignored.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, m_e
import yt

yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()

# Fetch the final particle position
xx_sim = ad['electron', 'particle_position_x'].v[0]
zz_sim = ad['electron', 'particle_position_z'].v[0]
ux_sim = ad['electron', 'particle_momentum_x'].v[0]/m_e
uz_sim = ad['electron', 'particle_momentum_z'].v[0]/m_e

# Fetch the quadrupole lattice data
quad_starts = np.array([float(x) for x in ds.parameters.get('quad.zstarts').split()])
quad_ends   = np.array([float(x) for x in ds.parameters.get('quad.zends').split()])
quad_lengths = quad_ends - quad_starts
quad_strengths_E = np.array([eval(x) for x in ds.parameters.get('quad.dEdx').split()])

# Fetch the initial position of the particle
x0 = [float(x) for x in ds.parameters.get('electron.single_particle_pos').split()]
ux0 = [float(x)*c for x in ds.parameters.get('electron.single_particle_vel').split()]

xx = x0[0]
zz = x0[2]
ux = ux0[0]
uz = ux0[2]

gamma = np.sqrt(uz**2/c**2 + 1.)
vz = uz/gamma

def applylens(x0, vx0, vz0, gamma, lens_length, lens_strength):
    """Use analytic solution of a particle with a transverse dependent field"""
    kb0 = np.sqrt(e/(m_e*gamma*vz0**2)*abs(lens_strength))
    if lens_strength >= 0.:
        x1 = x0*np.cos(kb0*lens_length) + (vx0/vz0)/kb0*np.sin(kb0*lens_length)
        vx1 = vz0*(-kb0*x0*np.sin(kb0*lens_length) + (vx0/vz0)*np.cos(kb0*lens_length))
    else:
        x1 = x0*np.cosh(kb0*lens_length) + (vx0/vz0)/kb0*np.sinh(kb0*lens_length)
        vx1 = vz0*(+kb0*x0*np.sinh(kb0*lens_length) + (vx0/vz0)*np.cosh(kb0*lens_length))
    return x1, vx1

# Integrate the particle using the analytic solution
for i in range(len(quad_starts)):
    z_lens = quad_starts[i]
    vx = ux/gamma
    dt = (z_lens - zz)/vz
    xx = xx + dt*vx
    xx, vx = applylens(xx, vx, vz, gamma, quad_lengths[i], quad_strengths_E[i])
    dt = quad_lengths[i]/vz
    ux = gamma*vx
    zz = z_lens + quad_lengths[i]

dt = (zz_sim - zz)/vz
vx = ux/gamma
xx = xx + dt*vx

# Compare the analytic to the simulated final values
print(f'Error in x position is {abs(np.abs((xx - xx_sim)/xx))}, which should be < 0.01')
print(f'Error in x velocity is {abs(np.abs((ux - ux_sim)/ux))}, which should be < 0.002')

assert abs(np.abs((xx - xx_sim)/xx)) < 0.01, Exception('error in x particle position')
assert abs(np.abs((ux - ux_sim)/ux)) < 0.002, Exception('error in x particle velocity')

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
