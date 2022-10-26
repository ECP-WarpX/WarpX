#!/usr/bin/env python3

import numpy as np
import scipy.constants as con

from pywarpx import picmi
from pywarpx.Amrex import amrex
from pywarpx.WarpX import warpx
from pywarpx.Algo import algo

GB = 1024**3

warpx.add_new_attr("numprocs", [1,1,1])
warpx.add_new_attr("do_current_centering", 1)

amrex.add_new_attr("the_arena_is_managed", 0)
amrex.add_new_attr("the_arena_init_size", 45*GB)
# amrex.add_new_attr("throw_exception", 0)
# amrex.add_new_attr("signal_handling", 0)

##########################
# physics parameters
##########################

V_domain_boundary = 0.0
V_embedded_boundary = 1.0


##########################
# numerics parameters
##########################

dt = 1e-11

# --- Nb time steps

max_steps = 1

# --- grid

nx = 64
ny = 64
nz = 64

xmin = -0.25
xmax = 0.25
ymin = -0.25
ymax = 0.25
zmin = 0.
zmax = 1.0


##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(
    number_of_cells = [nx, ny, nz],
    lower_bound = [xmin, ymin, zmin],
    upper_bound = [xmax, ymax, zmax],
    lower_boundary_conditions = ['neumann', 'neumann', 'dirichlet'],
    upper_boundary_conditions = ['neumann', 'neumann', 'neumann'],
    lower_boundary_conditions_particles = ['absorbing', 'absorbing', 'absorbing'],
    upper_boundary_conditions_particles = ['absorbing', 'absorbing', 'absorbing'],
    # warpx_potential_lo_x = V_domain_boundary,
    # warpx_potential_hi_x = V_domain_boundary,
    # warpx_potential_lo_y = V_domain_boundary,
    # warpx_potential_hi_y = V_domain_boundary,
    warpx_potential_lo_z = V_domain_boundary,
    # warpx_potential_hi_z = V_domain_boundary,
    warpx_blocking_factor=8,
    warpx_max_grid_size = 128
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-7,warpx_magnetostatic=True,warpx_self_fields_verbosity=3
)

embedded_boundary = picmi.EmbeddedBoundary(
    implicit_function="(x**2+y**2-radius**2)",
    potential=V_embedded_boundary,
    radius = 0.15
)


# Beam Current Density
current = 1000  # A
beam_r = 0.1 # m

J = current/(np.pi * beam_r**2)
beam_gamma = 10.
vz = con.c*np.sqrt(1. - 1./beam_gamma**2)
n0 = J/(con.e*vz)


beam_dist = picmi.AnalyticDistribution(
    density_expression='((x**2+y**2)<rmax**2)*n0',
    rmax=beam_r,
    n0=n0,
    directed_velocity=[0., 0., beam_gamma*vz],
)
beam = picmi.Species(particle_type='electron', name='beam', initial_distribution=beam_dist)
beam_layout = picmi.PseudoRandomLayout(n_macroparticles_per_cell=[2,2,2], grid=grid)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diag1',
    grid = grid,
    period = 1,
    data_list = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz','Ax', 'Ay', 'Az', 'Jx', 'Jy', 'Jz', 'phi', 'rho'],
    write_dir = '.',
    warpx_file_prefix = 'Python_MagnetostaticEB_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    time_step_size = dt,
    max_steps = max_steps,
    warpx_embedded_boundary=embedded_boundary,
    warpx_field_gathering_algo='momentum-conserving',
    warpx_current_deposition_algo='direct'
)

sim.add_species(beam, layout=beam_layout, initialize_self_field=True)

sim.add_diagnostic(field_diag)

##########################
# simulation run
##########################

sim.step(max_steps)
