# --- Input file to test the particle scraper and the Python wrappers
# --- to access the buffer of scraped particles.

import numpy as np
from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

D_CA = 0.067 # m

N_INERT = 9.64e20 # m^-3
T_INERT = 300.0 # K

M_ION = 6.67e-27 # kg

PLASMA_DENSITY = 2.56e14 # m^-3
T_ELEC = 30000.0 # K

DT = 1e-10

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 5
diagnostic_intervals = "::5"

# --- Grid
nx = 128
ny = 16

xmin = 0.0
ymin = 0.0
xmax = D_CA
ymax = D_CA / nx * ny

number_per_cell_each_dim = [32, 16]

##########################
# physics components
##########################

v_rms_elec = np.sqrt(constants.kb * T_ELEC / constants.m_e)
v_rms_ion = np.sqrt(constants.kb * T_INERT / M_ION)

uniform_plasma_elec = picmi.UniformDistribution(
    density = PLASMA_DENSITY,
    upper_bound = [None] * 3,
    rms_velocity = [v_rms_elec] * 3,
    directed_velocity = [0.] * 3
)

uniform_plasma_ion = picmi.UniformDistribution(
    density = PLASMA_DENSITY,
    upper_bound = [None] * 3,
    rms_velocity = [v_rms_ion] * 3,
    directed_velocity = [0.] * 3
)

electrons = picmi.Species(
    particle_type='electron', name='electrons',
    initial_distribution=uniform_plasma_elec,
    warpx_scrape_lo_x=1, warpx_scrape_hi_x=1
)
ions = picmi.Species(
    particle_type='He', name='he_ions',
    charge='q_e',
    initial_distribution=uniform_plasma_ion,
    warpx_scrape_lo_x=1, warpx_scrape_hi_x=1
)

##########################
# numerics components
##########################

grid = picmi.Cartesian2DGrid(
    number_of_cells = [nx, ny],
    lower_bound = [xmin, ymin],
    upper_bound = [xmax, ymax],
    bc_xmin = 'dirichlet',
    bc_xmax = 'dirichlet',
    bc_ymin = 'periodic',
    bc_ymax = 'periodic',
    warpx_potential_hi_x=5,
    lower_boundary_conditions_particles=['absorbing', 'periodic'],
    upper_boundary_conditions_particles=['absorbing', 'periodic'],
    warpx_max_grid_size = nx//4
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-6
)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diag1',
    grid = grid,
    period = diagnostic_intervals,
    data_list = ['rho_electrons', 'rho_he_ions', 'phi'],
    write_dir = '.',
    warpx_file_prefix = 'Python_particle_scraper_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    time_step_size = DT,
    max_steps = max_steps,
    verbose=True
)

sim.add_species(
    electrons,
    layout = picmi.GriddedLayout(
        n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
    )
)
sim.add_species(
    ions,
    layout = picmi.GriddedLayout(
        n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid
    )
)

sim.add_diagnostic(field_diag)

sim.step(max_steps - 1)

# check scraped particle buffer part way through the simulation
from pywarpx import _libwarpx

n = _libwarpx.get_num_particles_impacted_boundary("electrons", 'x_hi')
print("Number of electrons scraped:", n)
weights = _libwarpx.get_particle_boundary_buffer("electrons", 'x_hi', 'w', 0)
assert n == len(weights[0])+len(weights[1])

##########################
# simulation run
##########################

sim.step(1)
