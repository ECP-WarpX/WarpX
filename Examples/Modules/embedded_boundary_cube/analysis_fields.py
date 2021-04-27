#! /usr/bin/env python

import yt
import os, sys
from scipy.constants import mu_0, pi, c
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# This is a script that analyses the simulation results from
# the script `inputs_3d`. This simulates a TMmnp mode in a PEC cubic resonator.
# The magnetic field in the simulation is given (in theory) by:
# $$ B_x = \frac{-2\mu}{h^2}\, k_x k_z \sin(k_x x)\cos(k_y y)\cos(k_z z)\cos( \omega_p t)$$
# $$ B_y = \frac{-2\mu}{h^2}\, k_y k_z \cos(k_x x)\sin(k_y y)\cos(k_z z)\cos( \omega_p t)$$
# $$ B_z = \cos(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
# with
# $$ h^2 = k_x^2 + k_y^2 + k_z^2$$
# $$ k_x = \frac{m\pi}{L}$$
# $$ k_y = \frac{n\pi}{L}$$
# $$ k_z = \frac{p\pi}{L}$$

hi = [0.8, 0.8, 0.8]
lo = [-0.8, -0.8, -0.8]
ncells = [48, 48, 48]
dx = (hi[0] - lo[0])/ncells[0]
dy = (hi[1] - lo[1])/ncells[1]
dz = (hi[2] - lo[2])/ncells[2]
m = 0
n = 1
p = 1
Lx = 1
Ly = 1
Lz = 1
h_2 = (m * pi / Lx) ** 2 + (n * pi / Ly) ** 2 + (p * pi / Lz) ** 2
t = 1.3342563807926085e-08

# Compute the analytic solution
Bx_th = np.zeros(ncells)
By_th = np.zeros(ncells)
Bz_th = np.zeros(ncells)
for i in range(ncells[0]):
    for j in range(ncells[1]):
        for k in range(ncells[2]):
            x = i*dx + lo[0]
            y = (j+0.5)*dy + lo[1]
            z = k*dz + lo[2]

            By_th[i, j, k] = -2/h_2*mu_0*(n * pi/Ly)*(p * pi/Lz) * (np.cos(m * pi/Lx * (x - Lx/2)) *
                                                                    np.sin(n * pi/Ly * (y - Ly/2)) *
                                                                    np.cos(p * pi/Lz * (z - Lz/2)) *
                                                                    (-Lx/2 <= x < Lx/2) *
                                                                    (-Ly/2 <= y < Ly/2) *
                                                                    (-Lz/2 <= z < Lz/2) *
                                                                    np.cos(np.sqrt(2) *
                                                                    np.pi / Lx * c * t))

            x = i*dx + lo[0]
            y = j*dy + lo[1]
            z = (k+0.5)*dz + lo[2]
            Bz_th[i, j, k] = mu_0*(np.cos(m * pi/Lx * (x - Lx/2)) *
                                   np.cos(n * pi/Ly * (y - Ly/2)) *
                                   np.sin(p * pi/Lz * (z - Lz/2)) *
                                   (-Lx/2 <= x < Lx/2) *
                                   (-Ly/2 <= y < Ly/2) *
                                   (-Lz/2 <= z < Lz/2) *
                                   np.cos(np.sqrt(2) * np.pi / Lx * c * t))

# Open the right plot file
filename = sys.argv[1]
ds = yt.load(filename)
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)

rel_tol_err = 1e-1

# Compute relative l^2 error on By
By_sim = data['By'].to_ndarray()
rel_err_y = np.sqrt( np.sum(np.square(By_sim - By_th)) / np.sum(np.square(By_th)))
assert(rel_err_y < rel_tol_err)
# Compute relative l^2 error on Bz
Bz_sim = data['Bz'].to_ndarray()
rel_err_z = np.sqrt( np.sum(np.square(Bz_sim - Bz_th)) / np.sum(np.square(Bz_th)))
assert(rel_err_z < rel_tol_err)

test_name = os.path.split(os.getcwd())[1]

checksumAPI.evaluate_checksum(test_name, filename)
