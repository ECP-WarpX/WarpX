#!/usr/bin/env python3

from pywarpx import picmi

# Number of time steps
max_steps = 3

# Number of cells
nx = 128
ny = 32
nz = 128

# Physical domain
xmin = 0.
xmax = 4.
ymin = 0.
ymax = 4.
zmin = 0.
zmax = 4.

# Create grid
grid = picmi.Cartesian3DGrid(
    number_of_cells=[nx, ny, nz],
    warpx_max_grid_size=32,
    lower_bound=[xmin, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    lower_boundary_conditions=['periodic', 'periodic', 'periodic'],
    upper_boundary_conditions=['periodic', 'periodic', 'periodic']
)

# Electromagnetic solver
solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',
    cfl=0.99999
)

# Particles
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=picmi.UniformDistribution(
        density=1e14,
        rms_velocity=[0]*3,
        upper_bound=[xmax, ymax, 1.0]
    )
)
layout = picmi.GriddedLayout(
    n_macroparticle_per_cell=[1, 1, 1], grid=grid
)

# Reduced diagnostic
reduced_diag = picmi.ReducedDiagnostic(
    diag_type='LoadBalanceCosts',
    period=1,
    name='LBC'
)

# Diagnostic
diag = picmi.FieldDiagnostic(
    grid=grid,
    period=3,
    write_dir='.',
    warpx_file_prefix='Python_reduced_diags_loadbalancecosts_timers_plt'
)

# Set up simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    verbose=1,
    particle_shape=1,
    warpx_current_deposition_algo='esirkepov',
    warpx_field_gathering_algo='energy-conserving',
    warpx_load_balance_intervals=2
)

# Add species
sim.add_species(electrons, layout=layout)

# Add reduced diagnostics
sim.add_diagnostic(reduced_diag)

# Add diagnostics
sim.add_diagnostic(diag)

# Advance simulation until last time step
# sim.write_input_file("test_input")
sim.step(max_steps)
