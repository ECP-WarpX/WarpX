from mewarpx import util as mwxutil
mwxutil.init_libwarpx(ndim=2, rz=False)

from pywarpx import picmi

from mewarpx import assemblies, emission, mepicmi
from mewarpx.poisson_pseudo_1d import PoissonSolverPseudo1D

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store import diag_base
from mewarpx.mcc_wrapper import MCC

import mewarpx.mwxconstants as constants

####################################
# physics parameters
####################################

P_INERT = 2.0 # torr
T_INERT = 300.0 # K
N_INERT = (P_INERT * constants.torr_SI) / (constants.kb_J * T_INERT) # m^-3

D_CA = 1e-4 # m
V_bias = 30 # V

###################################
# numerics parameters
##################################

# --- Grid
nx = 8
ny = 128

xmin = 0.0
ymin = 0.0

xmax = D_CA / ny * nx
ymax = D_CA
number_per_cell_each_dim = [16, 16]

TOTAL_TIME = 1e-7 # s
DIAG_INTERVAL = 1.0e-9
DT = 0.5e-12 # s

max_steps = int(TOTAL_TIME / DT)
diag_steps = int(DIAG_INTERVAL / DT)
diagnostic_intervals = '::%i' % diag_steps

print('Setting up simulation with')
print(' dt = %.3e s' % DT)
print(' Total time = %.3e s (%i timesteps)' % (TOTAL_TIME, max_steps))
print(' Diag time = %.3e (%i timesteps)' % (DIAG_INTERVAL, diag_steps))

#####################################
# grid and solver
#####################################

grid = picmi.Cartesian2DGrid(
    number_of_cells=[nx, ny],
    lower_bound=[xmin, ymin],
    upper_bound=[xmax, ymax],
    bc_xmin='periodic',
    bc_xmax='periodic',
    bc_ymin='dirichlet',
    bc_ymax='dirichlet',
    bc_xmin_particles='periodic',
    bc_xmax_particles='periodic',
    bc_ymin_particles='absorbing',
    bc_ymax_particles='absorbing',
    moving_window_velocity=None,
    warpx_max_grid_size=128,
    warpx_potential_lo_z=0.0,
    warpx_potential_hi_z=V_bias,
)

solver = PoissonSolverPseudo1D(grid=grid)

#################################
# physics components
################################

electrons = mepicmi.Species(
    particle_type='electron',
    name='electrons',
    warpx_grid=grid,
    warpx_n_macroparticle_per_cell=number_per_cell_each_dim
)

ions = mepicmi.Species(
    particle_type='Ar',
    name='ar_ions',
    charge='q_e',
    warpx_grid=grid,
    warpx_n_macroparticle_per_cell=number_per_cell_each_dim
)

# MCC Collisions
mcc_wrapper = MCC(
    electrons, ions, T_INERT=T_INERT, N_INERT=N_INERT,
    exclude_collisions=['charge_exchange']
)

###################################
# diagnostics
##################################

field_diag = picmi.FieldDiagnostic(
    name='diags',
    grid=grid,
    period=diagnostic_intervals,
    data_list=['rho_electrons', 'rho_ar_ions', 'phi', 'J'],
    write_dir='.',
    warpx_file_prefix='diags'
)

#################################
# simulation setup
################################

mwxrun.simulation.solver = solver
mwxrun.simulation.time_step_size = DT
mwxrun.simulation.max_steps = max_steps

mwxrun.simulation.add_diagnostic(field_diag)
mwxrun.init_run()

######################################
# Add ME emission
#####################################
T_cathode = 1100.0 + 273.15 # K
WF_cathode = 2.0 # eV

cathode = assemblies.ZPlane(z=1e-10, zsign=-1, V=0, T=T_cathode,
                            WF=WF_cathode,
                            name='cathode')
emitter = emission.ZPlaneEmitter(conductor=cathode, T=T_cathode,
                                 use_Schottky=False)
injector = emission.ThermionicInjector(emitter=emitter, species=electrons,
                                       npart_per_cellstep=50,
                                       T=T_cathode, WF=WF_cathode,
                                       A=6e5)

####################################
# Add ME diagnostic
###################################

diag_base.TextDiag(diag_steps=diag_steps, preset_string='perfdebug')

##################################
# Simulation run
#################################
mwxrun.simulation.step(max_steps)
