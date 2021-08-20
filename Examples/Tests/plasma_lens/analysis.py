#! /usr/bin/env python

# Copyright 2021 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


"""
This script tests the plasma lens.
The input file sets up a series of plasma lens and propagates a particle through them.
The final position is compared to the analytic solution.
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

# Get final position of the particles
# The sorting by id is to get them in the correct order
# There are two particles, one moves in x, the other in y.
r_id = ad['electrons', 'particle_id'].v

xx_sim = ad['electrons', 'particle_position_x'].v[np.argsort(r_id)][0]
yy_sim = ad['electrons', 'particle_position_y'].v[np.argsort(r_id)][1]
zz_sim0 = ad['electrons', 'particle_position_z'].v[np.argsort(r_id)][0]
zz_sim1 = ad['electrons', 'particle_position_z'].v[np.argsort(r_id)][1]

ux_sim = ad['electrons', 'particle_momentum_x'].v[np.argsort(r_id)][0]/m_e
uy_sim = ad['electrons', 'particle_momentum_y'].v[np.argsort(r_id)][1]/m_e


def applylens(x0, vx0, vz0, lens_length, lens_strength):
    """Given the initial position at the start of the lens,
    calculate the end position at the exit.
    This assumes harmonic motion."""
    w = np.sqrt(e/m_e*lens_strength)
    A = np.sqrt(x0**2 + (vx0/w)**2)
    phi = np.arctan2(-vx0, w*x0)
    t = lens_length/vz0
    x1 = A*np.cos(w*t + phi)
    vx1 = -w*A*np.sin(w*t + phi)
    return x1, vx1

vel_z = eval(ds.parameters.get('my_constants.vel_z'))

plasma_lens_period = float(ds.parameters.get('particles.repeated_plasma_lens_period'))
plasma_lens_starts = [float(x) for x in ds.parameters.get('particles.repeated_plasma_lens_starts').split()]
plasma_lens_lengths = [float(x) for x in ds.parameters.get('particles.repeated_plasma_lens_lengths').split()]
plasma_lens_strengths_E = [eval(x) for x in ds.parameters.get('particles.repeated_plasma_lens_strengths_E').split()]
plasma_lens_strengths_B = [eval(x) for x in ds.parameters.get('particles.repeated_plasma_lens_strengths_B').split()]

clight = c

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

for i in range(len(plasma_lens_starts)):
    z_lens = i*plasma_lens_period + plasma_lens_starts[i]
    dt = (z_lens - zz)/uz
    tt = tt + dt
    xx = xx + dt*ux
    yy = yy + dt*uy
    lens_strength = plasma_lens_strengths_E[i] + plasma_lens_strengths_B[i]*vel_z
    xx, ux = applylens(xx, ux, uz, plasma_lens_lengths[i], lens_strength)
    yy, uy = applylens(yy, uy, uz, plasma_lens_lengths[i], lens_strength)
    dt = plasma_lens_lengths[i]/uz
    tt = tt + dt
    zz = z_lens + plasma_lens_lengths[i]

dt0 = (zz_sim0 - zz)/uz
dt1 = (zz_sim1 - zz)/uz
xx = xx + dt0*ux
yy = yy + dt1*uy

assert abs(np.abs((xx - xx_sim)/xx)) < 0.003, Exception('error in x particle position')
assert abs(np.abs((yy - yy_sim)/yy)) < 0.003, Exception('error in y particle position')
assert abs(np.abs((ux - ux_sim)/ux)) < 5.e-5, Exception('error in x particle velocity')
assert abs(np.abs((uy - uy_sim)/uy)) < 5.e-5, Exception('error in y particle velocity')

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)

