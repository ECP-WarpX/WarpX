"""
Monte-Carlo Collision script benchmark against case 1 results from
Turner et al. (2013) - https://doi.org/10.1063/1.4775084
"""

from mewarpx import util as mwxutil
mwxutil.init_libwarpx(ndim=2, rz=False)

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store import diag_base
from mewarpx import mcc_wrapper

from pywarpx import picmi
import pywarpx
from pywarpx import callbacks

import numpy as np
import time

import shutil
import yt

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
mcc_electrons, mcc_ions = mcc_wrapper.init_mcc(electrons, ions, N_INERT, T_INERT)

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
# WarpX and mewarpx initialization
##########################

mwxrun.init_run(simulation=sim)

print('Set up simulation with')
print('  dt = %.3e s' % DT)
print('  Total time = %.3e s (%i timesteps)' % (TOTAL_TIME, max_steps))
print('  Diag time = %.3e s (%i timesteps)' % (DIAG_INTERVAL, diag_steps))

##########################
# Add ME diagnostic
##########################

diag_base.TextDiag(diag_steps=diag_steps, preset_string='perfdebug')

##########################
# simulation run
##########################

import matplotlib.pyplot as plt

#my_solver = PoissonSolverPseudo1D(nx, ny, D_CA / nx)

def direct_solve():
    rho_data = mwxrun.get_rho_grid()[0]
    rho = rho_data[:,:,0].T
    phi = my_solver.solve(rho,
        VOLTAGE * np.sin(2.0*np.pi*FREQ* (mwxrun.get_it()-1.0)*mwxrun.get_dt())
    )
    mwxrun.set_phi_grid(phi)

def plot_phi():
    phi_data = mwxrun.get_gathered_phi_grid()
    # phi_data = mwxrun.get_phi_grid()
    if mwxrun.me == 0:
        print(phi_data)
        phi_data = phi_data[0]
        print(mwxrun.me, phi_data.shape)
        plt.plot(np.mean(phi_data[1:-1], axis=1), 'o-')
    # plt.plot(phi_data[:,4], 'o-')
    # plt.plot(phi[:,4], 'o-')

# comment line below to use the multigrid solver
# callbacks.installfieldsolver(direct_solve)

callbacks.installafterstep(plot_phi)

sim.step(5)
if mwxrun.me == 0:
    plt.ylim(-2, 28)
    plt.xlabel('Cell number')
    plt.ylabel('$\phi$ (eV)')
    plt.title(
        'Electrostatic potential at different times gathered\n'
        'to the root proc from a 2 proc simulation'
    )
    plt.grid()
    plt.savefig('gathered_phi.png')
    plt.show()
