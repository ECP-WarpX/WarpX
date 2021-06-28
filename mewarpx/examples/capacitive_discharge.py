"""
Monte-Carlo Collision script benchmark against case 1 results from
Turner et al. (2013) - https://doi.org/10.1063/1.4775084
"""
from pywarpx import picmi
import pywarpx

import numpy as np
import time

import shutil
import yt

from minerva import util as minutil

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

##########################
# numerics parameters
##########################

# --- Grid
nx = 128
ny = 16

xmin = 0.0
ymin = 0.0
xmax = D_CA
ymax = D_CA / nx * ny

number_per_cell_each_dim = [32, 16]

DT = 1.0 / (400 * FREQ)

# Total simulation time in seconds
TOTAL_TIME = 10.0 * DT # 1280 / FREQ
# Time (in seconds) between diagnostic evaluations
DIAG_INTERVAL = 2.0 * DT # 32 / FREQ

# --- Number of time steps
max_steps = int(TOTAL_TIME / DT)
diag_steps = int(DIAG_INTERVAL / DT)
diagnostic_intervals = "::%i" % diag_steps #"%i:" % (max_steps - diag_steps + 1)

print('Setting up simulation with')
print('  dt = %.3e s' % DT)
print('  Total time = %.3e s (%i timesteps)' % (TOTAL_TIME, max_steps))
print('  Diag time = %.3e s (%i timesteps)' % (DIAG_INTERVAL, diag_steps))

##########################
# physics components
##########################

v_rms_elec = np.sqrt(minutil.kb_J * T_ELEC / minutil.m_e)
v_rms_ion = np.sqrt(minutil.kb_J * T_INERT / M_ION)

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
cross_sec_direc = '../../../warpx-data/MCC_cross_sections/He/'
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
    name = 'diags',
    grid = grid,
    period = diagnostic_intervals,
    data_list = ['rho_electrons', 'rho_he_ions', 'phi'],
    write_dir = 'diags/',
    # warpx_file_prefix = 'diags',
    #warpx_format = 'openpmd',
    #warpx_openpmd_backend='h5'
)
'''
restart_dumps = picmi.FieldDiagnostic(
    name = 'checkpoints',
    warpx_format = 'checkpoint',
    grid = grid,
    period = diagnostic_intervals//2,
    write_dir = './restarts',
    # warpx_file_prefix = 'diags',
)
'''
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
# sim.add_diagnostic(restart_dumps)

##########################
# simulation run
##########################

#sim.write_input_file(file_name = 'input2d')

# step the simulation to the point where output should start
sim.step(max_steps)
