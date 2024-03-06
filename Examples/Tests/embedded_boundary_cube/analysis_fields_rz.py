#!/usr/bin/env python3

import os
import sys

import numpy as np
from scipy.constants import c, pi
from scipy.special import j0, j1
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# This is a script that analyses the simulation results from
# the script `inputs_rz`. This simulates an EM mode in a PEC cylindrical resonator.
# The magnetic field in the simulation is given (in theory) by:
# $$ B_z = J_0(k_r r)\sin(k_z z)\cos(\omega t)$$
# $$ B_r = - k_z/k_r J_1(k_r r)\sin(k_z z)\cos(\omega t)
# with
# $$ k_r = \frac{2\alpha}{L}$$
# $$ k_z = \frac{p\pi}{L}$$
# $$ \omega^2 = c^2 (k_r^2 + k_z^2)$$
# where $\alpha$ is the first zero of the Bessel function J_0

hi = [0.8, 0.8]
lo = [-0.8, -0.8]
ncells = [32, 32, 1]
dx = (hi[0] - lo[0]) / ncells[0]
dz = (hi[1] - lo[1]) / ncells[1]
m = 0
n = 1
L = 1

# Open the right plot file
filename = sys.argv[1]
ds = yt.load(filename)
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)

t = ds.current_time.to_value()

alpha = 2.40482556
omega = c * np.sqrt( np.pi**2 + (2*alpha)**2 )/ L

# Compute the analytic solution
Br_th = np.zeros(ncells)
Bz_th = np.zeros(ncells)
for i in range(ncells[0]):
    for j in range(ncells[1]):
        x = (i+0.5)*dx + lo[0]
        z = j*dz + lo[1]
        Br_th[i, j, 0] = - pi/(2*alpha) * j1( 2*alpha/L * x ) * \
                        np.cos(n * pi / L * (z - L / 2)) * \
                        (-L / 2 <= x < L / 2) * \
                        (-L / 2 <= z < L / 2) * \
                        np.cos(omega * t)
        x = i*dx + lo[0]
        z = (j+0.5)*dz + lo[1]
        Bz_th[i, j, 0] = j0( 2*alpha/L * x ) * \
                        np.cos(n * pi / L * (z - L / 2)) * \
                        (-L / 2 <= x < L / 2) * \
                        (-L / 2 <= z < L / 2) * \
                        np.cos(omega * t)

rel_tol_err = 1e-1

# Compute relative l^2 error on Br
Br_sim = data[('mesh','Bx')].to_ndarray()
rel_err_r = np.sqrt( np.sum(np.square(Br_sim - Br_th)) / np.sum(np.square(Br_th)))
assert(rel_err_r < rel_tol_err)
# Compute relative l^2 error on Bz
Bz_sim = data[('mesh','Bz')].to_ndarray()
rel_err_z = np.sqrt( np.sum(np.square(Bz_sim - Bz_th)) / np.sum(np.square(Bz_th)))
assert(rel_err_z < rel_tol_err)

test_name = os.path.split(os.getcwd())[1]

checksumAPI.evaluate_checksum(test_name, filename)
