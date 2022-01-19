#!/usr/bin/env python3
#
# --- Input file to test the saving of old particle positions

import numpy as np
from pywarpx import picmi

constants = picmi.constants

##########################
# numerics parameters
##########################

dt = 7.5e-10

# --- Nb time steps

max_steps = 10

# --- grid

nx = 64
nz = 64

xmin = 0
xmax = 0.03
zmin = 0
zmax = 0.03


##########################
# numerics components
##########################

grid = picmi.Cartesian2DGrid(
    number_of_cells = [nx, nz],
    lower_bound = [xmin, zmin],
    upper_bound = [xmax, zmax],
    lower_boundary_conditions = ['dirichlet', 'periodic'],
    upper_boundary_conditions = ['dirichlet', 'periodic'],
    lower_boundary_conditions_particles = ['absorbing', 'periodic'],
    upper_boundary_conditions_particles = ['absorbing', 'periodic'],
    moving_window_velocity = None,
    warpx_max_grid_size = 32
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-6,
    warpx_self_fields_verbosity=0
)

##########################
# physics components
##########################

uniform_plasma_elec = picmi.UniformDistribution(
    density = 1e15,
    upper_bound = [None] * 3,
    rms_velocity = [np.sqrt(constants.kb * 1e3 / constants.m_e)] * 3,
    directed_velocity = [0.] * 3
)

electrons = picmi.Species(
    particle_type='electron', name='electrons',
    initial_distribution=uniform_plasma_elec,
    warpx_save_previous_position=True
)

##########################
# diagnostics
##########################

field_diag = picmi.ParticleDiagnostic(
    species=electrons,
    name = 'diag1',
    data_list=['previous_positions'],
    period = 10,
    write_dir = '.',
    warpx_file_prefix = 'Python_prev_positions_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    time_step_size = dt,
    max_steps = max_steps,
    verbose = 1
)

sim.add_species(
    electrons,
    layout = picmi.GriddedLayout(
        n_macroparticle_per_cell=[1, 1], grid=grid
    )
)
sim.add_diagnostic(field_diag)

##########################
# simulation run
##########################

sim.step(max_steps - 1)

##########################
# check that the new PIDs
# exist
##########################

assert (sim.extension.get_particle_comp_index('electrons', 'prev_x') > 0)
assert (sim.extension.get_particle_comp_index('electrons', 'prev_z') > 0)

prev_z_vals = sim.extension.get_particle_arrays(
    'electrons', 'prev_z', 0
)
for z_vals in prev_z_vals:
    assert np.all(z_vals < zmax)

##########################
# take the final sim step
##########################

sim.step(1)
