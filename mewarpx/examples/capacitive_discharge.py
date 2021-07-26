"""
Monte-Carlo Collision script based on case 1 from
Turner et al. (2013) - https://doi.org/10.1063/1.4775084
"""

from mewarpx import util as mwxutil
mwxutil.init_libwarpx(ndim=2, rz=False)

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store import diag_base
from mewarpx import mepicmi, emission, mcc_wrapper, poisson_pseudo_1d

from pywarpx import picmi
import pywarpx

import numpy as np

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

SEED_NPPC = 16 * 32

##########################
# numerics parameters
##########################

# --- Grid
nx = 8
nz = 128

xmin = 0.0
zmin = 0.0
xmax = D_CA / nz * nx
zmax = D_CA

DT = 1.0 / (400 * FREQ)

# Total simulation time in seconds
TOTAL_TIME = 500.0 * DT # 1280 / FREQ
# Time (in seconds) between diagnostic evaluations
# DIAG_INTERVAL = 100.0 * DT # 32 / FREQ

# --- Number of time steps
max_steps = int(TOTAL_TIME / DT)
diag_steps = 100
diagnostic_intervals = "400::10"

print('Setting up simulation with')
print('  dt = %.3e s' % DT)
print('  Total time = %.3e s (%i timesteps)' % (TOTAL_TIME, max_steps))

##########################
# physics components
##########################

anode_voltage = lambda t: VOLTAGE * np.sin(2.0 * np.pi * FREQ * t)
# anode_voltage = f"{VOLTAGE}*sin(2*pi*{FREQ:.5e}*t)"

##########################
# declare the simulation grid
##########################

mwxrun.init_grid(xmin, xmax, zmin, zmax, nx, nz)
mwxrun.grid.potential_zmax = anode_voltage

##########################
# declare species
##########################

electrons = mepicmi.Species(particle_type='electron', name='electrons')
ions = mepicmi.Species(particle_type='He', name='he_ions', charge='q_e')

##########################
# neutral plasma injection
##########################

vol_emitter = emission.UniformDistributionVolumeEmitter(T=T_ELEC)

plasma_injector = emission.PlasmaInjector(
    emitter=vol_emitter, species1=electrons, species2=ions,
    npart=2 * SEED_NPPC * nx * nz,
    T_2=T_INERT, plasma_density=PLASMA_DENSITY
)

##########################
# collision physics
##########################

# MCC collisions
mcc = mcc_wrapper.MCC(
    electrons, ions, T_INERT=T_INERT, N_INERT=N_INERT,
    exclude_collisions=['charge_exchange']
)

##########################
# declare solver
##########################

# solver = picmi.ElectrostaticSolver(
#    grid=mwxrun.grid, method='Multigrid', required_precision=1e-12
# )
solver = poisson_pseudo_1d.PoissonSolverPseudo1D(grid=mwxrun.grid)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diags',
    grid = mwxrun.grid,
    period = diagnostic_intervals,
    data_list = ['rho_electrons', 'rho_he_ions', 'phi'],
    write_dir = 'diags/',
)

##########################
# simulation setup
##########################

mwxrun.simulation.solver = solver
mwxrun.simulation.time_step_size = DT
mwxrun.simulation.max_steps = max_steps

mwxrun.simulation.add_diagnostic(field_diag)

##########################
# WarpX and mewarpx initialization
##########################

mwxrun.init_run()

##########################
# Add ME diagnostic
##########################

diag_base.TextDiag(diag_steps=diag_steps, preset_string='perfdebug')

##########################
# simulation run
##########################

mwxrun.simulation.step()

##########################
# collect diagnostics
##########################

if mwxrun.me == 0:
    import glob
    import yt

    data_dirs = glob.glob('diags/diags*')
    if len(data_dirs) == 0:
        raise RuntimeError("No data files found.")

    for ii, data in enumerate(data_dirs):

        datafolder = data
        print('Reading ', datafolder, '\n')
        ds = yt.load( datafolder )
        grid_data = ds.covering_grid(
            level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
        )
        if ii == 0:
            rho_data = np.mean(
                grid_data['rho_he_ions'].to_ndarray()[:,:,0], axis=0
            ) / constants.q_e
        else:
            rho_data += np.mean(
                grid_data['rho_he_ions'].to_ndarray()[:,:,0], axis=0
            ) / constants.q_e
    rho_data /= (ii + 1)
    np.save('direct_solver_avg_rho_data.npy', rho_data)
