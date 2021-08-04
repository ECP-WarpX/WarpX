from mewarpx import util as mwxutil
mwxutil.init_libwarpx(ndim=2, rz=False)

from mewarpx import assemblies, emission, mepicmi
from mewarpx.poisson_pseudo_1d import PoissonSolverPseudo1D

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store import diag_base
from mewarpx.mcc_wrapper import MCC
from mewarpx.diags_store.field_diagnostic import FieldDiagnostic

import mewarpx.mwxconstants as constants

import sys
import argparse

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
nz = 128

xmin = 0.0
zmin = 0.0

xmax = D_CA / nz * nx
zmax = D_CA

TOTAL_TIME = 1e-11 # s
DIAG_INTERVAL = 1.0e-9
DT = 0.5e-12 # s

max_steps = int(TOTAL_TIME / DT)
diag_steps = int(DIAG_INTERVAL / DT)

parser = argparse.ArgumentParser()
parser.add_argument('--steps', help='set the number of simulation steps manually', type=int)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

print('Setting up simulation with')
print(f'  dt = {DT:3e} s')
if args.steps:
    print(f' Total time = {DT*max_steps:3e} ({args.steps} timesteps)')
else:
    print(f' Total time = {TOTAL_TIME:3e} s ({max_steps} timesteps)')

#####################################
# grid and solver
#####################################

mwxrun.init_grid(xmin, xmax, zmin, zmax, nx, nz)
mwxrun.grid.potential_zmin = 0.0
mwxrun.grid.potential_zmax = V_bias

solver = PoissonSolverPseudo1D(grid=mwxrun.grid)

#################################
# physics components
################################

electrons = mepicmi.Species(
    particle_type='electron',
    name='electrons',
)

ions = mepicmi.Species(
    particle_type='Ar',
    name='ar_ions',
    charge='q_e',
)

# MCC Collisions
mcc_wrapper = MCC(
    electrons, ions, T_INERT=T_INERT, N_INERT=N_INERT,
    exclude_collisions=['charge_exchange']
)

###################################
# diagnostics
##################################

field_diag = FieldDiagnostic(
    diag_steps=diag_steps,
    diag_data_list=['rho_electrons', 'rho_ar_ions', 'phi', 'J'],
    grid=mwxrun.grid,
    name='diags',
    write_dir='diags/'
)

#################################
# simulation setup
################################

mwxrun.simulation.solver = solver
mwxrun.simulation.time_step_size = DT
mwxrun.simulation.max_steps = max_steps

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

if args.steps:
    mwxrun.simulation.step(args.steps)
else:
    mwxrun.simulation.step(max_steps)

print(electrons.get_array_from_pid("w"))
print(electrons.get_array_from_pid("E_total"))
