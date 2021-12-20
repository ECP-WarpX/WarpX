#!/usr/bin/env python3

# Copyright 2021 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


"""
This script tests the plasma lens.
The input file sets up a series of plasma lens and propagates two particles through them.
One particle is in the X plane, the other the Y plane.
The final positions are compared to the analytic solutions.
The motion is slow enough that relativistic effects are ignored.
"""

import sys
import os
import yt
import numpy as np
from scipy.constants import e, m_e, c
yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]
ds = yt.load( filename )
ad = ds.all_data()

# Get final position of the particles.
# There are two particles, one moves in x, the other in y.
# The particles may not be in order, so determine which is which
# by looking at their max positions in the respective planes.
i0 = np.argmax(np.abs(ad['electrons', 'particle_position_x'].v))
i1 = np.argmax(np.abs(ad['electrons', 'particle_position_y'].v))

xx_sim = ad['electrons', 'particle_position_x'].v[i0]
yy_sim = ad['electrons', 'particle_position_y'].v[i1]
zz_sim0 = ad['electrons', 'particle_position_z'].v[i0]
zz_sim1 = ad['electrons', 'particle_position_z'].v[i1]

ux_sim = ad['electrons', 'particle_momentum_x'].v[i0]/m_e
uy_sim = ad['electrons', 'particle_momentum_y'].v[i1]/m_e

if 'warpx.gamma_boost' in ds.parameters:
    gamma_boost = float(ds.parameters.get('warpx.gamma_boost'))
    uz_boost = np.sqrt(gamma_boost*gamma_boost - 1.)*c
    time = ds.current_time.to_value()
    zz_sim0 = gamma_boost*zz_sim0 + uz_boost*time
    zz_sim1 = gamma_boost*zz_sim1 + uz_boost*time


def applylens(x0, vx0, vz0, gamma, lens_length, lens_strength):
    kb0 = np.sqrt(e/(m_e*gamma*vz0**2)*lens_strength)
    x1 = x0*np.cos(kb0*lens_length) + (vx0/vz0)/kb0*np.sin(kb0*lens_length)
    vx1 = vz0*(-kb0*x0*np.sin(kb0*lens_length) + (vx0/vz0)*np.cos(kb0*lens_length))
    return x1, vx1

clight = c
vel_z = eval(ds.parameters.get('my_constants.vel_z'))

plasma_lens_period = float(ds.parameters.get('particles.repeated_plasma_lens_period'))
plasma_lens_starts = [float(x) for x in ds.parameters.get('particles.repeated_plasma_lens_starts').split()]
plasma_lens_lengths = [float(x) for x in ds.parameters.get('particles.repeated_plasma_lens_lengths').split()]
plasma_lens_strengths_E = [eval(x) for x in ds.parameters.get('particles.repeated_plasma_lens_strengths_E').split()]
plasma_lens_strengths_B = [eval(x) for x in ds.parameters.get('particles.repeated_plasma_lens_strengths_B').split()]


x0 = float(ds.parameters.get('electrons.multiple_particles_pos_x').split()[0])
y0 = float(ds.parameters.get('electrons.multiple_particles_pos_y').split()[1])
z0 = float(ds.parameters.get('electrons.multiple_particles_pos_z').split()[0])
ux0 = float(ds.parameters.get('electrons.multiple_particles_vel_x').split()[0])*c
uy0 = float(ds.parameters.get('electrons.multiple_particles_vel_y').split()[1])*c
uz0 = eval(ds.parameters.get('electrons.multiple_particles_vel_z').split()[0])*c

tt = 0.
xx = x0
yy = y0
zz = z0
ux = ux0
uy = uy0
uz = uz0

gamma = np.sqrt(uz0**2/c**2 + 1.)
vz = uz/gamma

for i in range(len(plasma_lens_starts)):
    z_lens = i*plasma_lens_period + plasma_lens_starts[i]
    vx = ux/gamma
    vy = uy/gamma
    dt = (z_lens - zz)/vz
    tt = tt + dt
    xx = xx + dt*vx
    yy = yy + dt*vy
    lens_strength = plasma_lens_strengths_E[i] + plasma_lens_strengths_B[i]*vel_z
    xx, vx = applylens(xx, vx, vz, gamma, plasma_lens_lengths[i], lens_strength)
    yy, vy = applylens(yy, vy, vz, gamma, plasma_lens_lengths[i], lens_strength)
    dt = plasma_lens_lengths[i]/vz
    tt = tt + dt
    ux = gamma*vx
    uy = gamma*vy
    zz = z_lens + plasma_lens_lengths[i]

dt0 = (zz_sim0 - zz)/vz
dt1 = (zz_sim1 - zz)/vz
vx = ux/gamma
vy = uy/gamma
xx = xx + dt0*vx
yy = yy + dt1*vy

print(f'Error in x position is {abs(np.abs((xx - xx_sim)/xx))}, which should be < 0.02')
print(f'Error in y position is {abs(np.abs((yy - yy_sim)/yy))}, which should be < 0.02')
print(f'Error in x velocity is {abs(np.abs((ux - ux_sim)/ux))}, which should be < 0.002')
print(f'Error in y velocity is {abs(np.abs((uy - uy_sim)/uy))}, which should be < 0.002')

assert abs(np.abs((xx - xx_sim)/xx)) < 0.02, Exception('error in x particle position')
assert abs(np.abs((yy - yy_sim)/yy)) < 0.02, Exception('error in y particle position')
assert abs(np.abs((ux - ux_sim)/ux)) < 0.002, Exception('error in x particle velocity')
assert abs(np.abs((uy - uy_sim)/uy)) < 0.002, Exception('error in y particle velocity')

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
