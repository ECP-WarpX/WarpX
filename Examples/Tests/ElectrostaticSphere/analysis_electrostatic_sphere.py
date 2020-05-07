#! /usr/bin/env python
#
# Copyright 2019-2020 David Bizzozero
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script tests the expansion of a uniformly charged sphere of electrons.

The input file inputs_3d/inputs is used: it defines a sphere of charge with
given initial conditions. The sphere of charge starts at rest and will expand
due to Coulomb interactions. The test will check if the sphere has expanded at
the correct speed and that the electric field is accurately modeled against a
known analytic solution. While the radius r(t) is not analytically known, its
inverse t(r) can be solved for exactly.
"""
import numpy as np
from scipy.optimize import fsolve
import sys
import yt
yt.funcs.mylog.setLevel(0)

# Constants and initial conditions
xmin = ymin = zmin = -0.5  #Box dimensions
xmax = ymax = zmax = 0.5
eps_0 = 8.8541878128e-12  #Vacuum Permittivity in C/(V*m)
m_e = 9.10938356e-31  #Electron mass in kg
q_e = -1.60217662e-19  #Electron charge in C
pi = np.pi  #Circular constant of the universe
r_0 = 0.1  #Initial radius of sphere
dt = 1e-6  #Timestep used for simulation
t_max = 3e-5  #Final time of simulation
q_tot = -1e-15  #Total charge of sphere in C

# Define functions for exact forms of v(r), t(r), Er(r) with r as the radius of
# the sphere. The sphere starts with initial radius r_0 and this radius expands
# with zero initial velocity. Note: v(t) and r(t) are not known analytically but
# v(r) and t(r) can be solved analytically.
#
# The solution r(t) solves the ODE: r''(t) = a/(r(t)**2) with initial conditions
# r(0) = r_0, r'(0) = 0, and a = q_e*q_tot/(4*pi*eps_0*m_e)
#
# Note: E might be saved on a staggered time grid shifted by dt/2 (check this!)
v_exact = lambda r: np.sqrt((q_e*q_tot)/(2*pi*m_e*eps_0)*(1/r_0-1/r))
t_exact = lambda r: np.sqrt(r_0**3*2*pi*m_e*eps_0/(q_e*q_tot)) \
    * (np.sqrt(r/r_0-1)*np.sqrt(r/r_0) \
       + np.log(np.sqrt(r/r_0-1)+np.sqrt(r/r_0)))
func = lambda rho: t_exact(rho)-t_max+dt/2  #Objective function to find r(t_max)
r_end = fsolve(func,r_0)[0]  #Numerically solve for r(t_max)
E_exact = lambda r: np.sign(r)*(q_tot/(4*pi*eps_0*r**2)*(abs(r)>=r_end) \
    + q_tot*abs(r)/(4*pi*eps_0*r_end**3)*(abs(r)<r_end))

# Open plotfile specified in command line
filename = sys.argv[1]
ds = yt.load( filename )

# Load data pertaining to fields
data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                        dims=ds.domain_dimensions)
Ex = data['Ex'].to_ndarray()
Ey = data['Ey'].to_ndarray()
Ez = data['Ez'].to_ndarray()
nx, ny, nz = Ex.shape

# Extract longitudinal fields along x-, y-, and z-axes
Ex_axis = Ex[:,int(ny/2),int(nz/2)]
Ey_axis = Ey[int(nx/2),:,int(nz/2)]
Ez_axis = Ez[int(nx/2),int(ny/2),:]

# Compute cell centers for grid
x_cell_edges = np.linspace(xmin,xmax,nx+1)
x_cell_centers = (x_cell_edges[1::]+x_cell_edges[0:-1])/2
y_cell_edges = np.linspace(ymin,ymax,ny+1)
y_cell_centers = (y_cell_edges[1::]+y_cell_edges[0:-1])/2
z_cell_edges = np.linspace(zmin,zmax,nz+1)
z_cell_centers = (z_cell_edges[1::]+z_cell_edges[0:-1])/2

# Extract subgrid away from boundary (exact solution assumes infinite/open
# domain but WarpX solution assumes perfect conducting walls)
x_sub_grid = x_cell_centers[int(nx/4):int(3*nx/4)]
y_sub_grid = y_cell_centers[int(ny/4):int(3*ny/4)]
z_sub_grid = z_cell_centers[int(nz/4):int(3*nz/4)]

# Exact solution of field along Cartesian axes
Ex_exact_grid = E_exact(x_sub_grid)
Ey_exact_grid = E_exact(y_sub_grid)
Ez_exact_grid = E_exact(z_sub_grid)

# WarpX solution of field along Cartesian axes
Ex_grid = Ex_axis[int(nx/4):int(3*nx/4)]
Ey_grid = Ey_axis[int(ny/4):int(3*nz/4)]
Ez_grid = Ez_axis[int(ny/4):int(3*nz/4)]

# Define approximate L2 norm error between exact and numerical solutions
L2_error_x = np.sqrt(sum((Ex_exact_grid-Ex_grid)**2)) \
    / np.sqrt(sum((Ex_exact_grid)**2))
L2_error_y = np.sqrt(sum((Ey_exact_grid-Ey_grid)**2)) \
    / np.sqrt(sum((Ey_exact_grid)**2))
L2_error_z = np.sqrt(sum((Ez_exact_grid-Ez_grid)**2)) \
    / np.sqrt(sum((Ez_exact_grid)**2))

print("L2 error along x-axis = %s" %L2_error_x)
print("L2 error along y-axis = %s" %L2_error_y)
print("L2 error along z-axis = %s" %L2_error_z)

assert L2_error_x < 0.05
assert L2_error_y < 0.05
assert L2_error_z < 0.05
