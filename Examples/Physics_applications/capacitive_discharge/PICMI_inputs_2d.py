# --- Input file for MCC testing. There is already a test of the MCC
# --- functionality so this just tests the PICMI interface

import numpy as np
from pywarpx import picmi

constants = picmi.constants

##########################
# physics parameters
##########################

D_CA = 0.067 # m

N_INERT = 9.64e20 # m^-3
T_INERT = 300.0 # K

FREQ = 13.56e6 # MHz

VOLTAGE = 450.0

M_ION = 6.67e-27 # kg

PLASMA_DENSITY = 2.56e14 # m^-3
T_ELEC = 30000.0 # K

DT = 1.0 / (400 * FREQ)

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 10
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
    initial_distribution=uniform_plasma_elec
)
ions = picmi.Species(
    particle_type='He', name='he_ions',
    charge='q_e',
    initial_distribution=uniform_plasma_ion
)

# MCC collisions
cross_sec_direc = '../../../../warpx-data/MCC_cross_sections/He/'
mcc_electrons = picmi.MCCCollisions(
    name='coll_elec',
    species=electrons,
    background_density=N_INERT,
    background_temperature=T_INERT,
    background_mass=ions.mass,
    scattering_processes={
        'elastic' : {
            'cross_section' : cross_sec_direc+'electron_scattering.dat'
        },
        'excitation1' : {
            'cross_section': cross_sec_direc+'excitation_1.dat',
            'energy' : 19.82
        },
        'excitation2' : {
            'cross_section': cross_sec_direc+'excitation_2.dat',
            'energy' : 20.61
        },
        'ionization' : {
            'cross_section' : cross_sec_direc+'ionization.dat',
            'energy' : 24.55,
            'species' : ions
        },
    }
)

mcc_ions = picmi.MCCCollisions(
    name='coll_ion',
    species=ions,
    background_density=N_INERT,
    background_temperature=T_INERT,
    scattering_processes={
        'elastic' : {
            'cross_section' : cross_sec_direc+'ion_scattering.dat'
        },
        'back' : {
            'cross_section' : cross_sec_direc+'ion_back_scatter.dat'
        },
        # 'charge_exchange' : {
        #    'cross_section' : cross_sec_direc+'charge_exchange.dat'
        # }
    }
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
    warpx_potential_hi_x = "%.1f*sin(2*pi*%.5e*t)" % (VOLTAGE, FREQ),
    lower_boundary_conditions_particles=['absorbing', 'periodic'],
    upper_boundary_conditions_particles=['absorbing', 'periodic'],
    moving_window_velocity = None,
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
    warpx_file_prefix = 'Python_background_mcc_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    time_step_size = DT,
    max_steps = max_steps,
    warpx_collisions=[mcc_electrons, mcc_ions]
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

##########################
# simulation run
##########################

sim.step(max_steps)
