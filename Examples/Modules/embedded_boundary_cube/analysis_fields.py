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
            y = j*dy + lo[1]
            z = k*dz + lo[2]

            Bx_th[i, j, k] = -2/h_2*mu_0*(m * pi/Lx)*(p * pi/Lz) * (np.sin(m * pi/Lx * (x - Lx/2)) *
                                                                    np.sin(n * pi/Ly * (y - Ly/2)) *
                                                                    np.cos(p * pi/Lz * (z - Lz/2)) *
                                                                    (-Lx/2 < x < Lx/2) *
                                                                    (-Ly/2 < y < Ly/2) *
                                                                    (-Lz/2 < z < Lz/2) *
                                                                    np.cos(np.sqrt(2) *
                                                                           np.pi / Lx * c * t))

            By_th[i, j, k] = -2/h_2*mu_0*(n * pi/Ly)*(p * pi/Lz) * (np.cos(m * pi/Lx * (x - Lx/2)) *
                                                                    np.sin(n * pi/Ly * (y - Ly/2)) *
                                                                    np.cos(p * pi/Lz * (z - Lz/2)) *
                                                                    (-Lx/2 < x < Lx/2) *
                                                                    (-Ly/2 < y < Ly/2) *
                                                                    (-Lz/2 < z < Lz/2) *
                                                                    np.cos(np.sqrt(2) *
                                                                    np.pi / Lx * c * t))

            Bz_th[i, j, k] = mu_0*(n * pi/Ly)*(p * pi/Lz) * (np.cos(m * pi/Lx * (x - Lx/2)) *
                                                             np.cos(n * pi/Ly * (y - Ly/2)) *
                                                             np.sin(p * pi/Lz * (z - Lz/2)) *
                                                             (-Lx/2 < x < Lx/2) *
                                                             (-Ly/2 < y < Ly/2) *
                                                             (-Lz/2 < z < Lz/2) *
                                                             np.cos(np.sqrt(2) *
                                                                    np.pi / Lx * c * t))

# Compute the analytic solution
path, dirs, files = next(os.walk("diags/"))
tsteps = np.zeros_like(dirs, dtype=int)
for i, d in enumerate(dirs):
    # Strip off diags
    d = d[4:]
    # Strip off leading zeros
    while d[0] == 0:
        d = d[1:]
    tsteps[i] = int(d)

# Open the right plot file
filename = 'diags/diag00208/'

ds = yt.load(filename)
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)

tol_err = 1e-5

# Compute l^2 error on Bx
Bx_sim = data['Bx'].to_ndarray()
err_x = np.sqrt(np.sum(np.square(Bx_sim - Bx_th))*dx*dy*dz)
assert(err_x < tol_err)
# Compute l^2 error on By
By_sim = data['By'].to_ndarray()
err_y = np.sqrt(np.sum(np.square(By_sim - By_th))*dx*dy*dz)
assert(err_y < tol_err)

# Compute l^2 error on Bz
Bz_sim = data['Bz'].to_ndarray()
err_z = np.sqrt(np.sum(np.square(Bz_sim - Bz_th))*dx*dy*dz)
assert(err_z < tol_err)

test_name = os.path.split(os.getcwd())[1]

checksumAPI.evaluate_checksum(test_name, filename)
