# This is a script that adds particle components at runtime,
# then performs checkpoint / restart and compares the result
# to the original simulation.

from pywarpx import picmi
import numpy as np
import sys

##########################
# physics parameters
##########################

dt = 7.5e-10

##########################
# numerics parameters
##########################

max_steps = 10

nx = 64
ny = 64

xmin = 0
xmax = 0.03
ymin = 0
ymax = 0.03

##########################
# numerics components
##########################

grid = picmi.Cartesian2DGrid(
    number_of_cells = [nx, ny],
    lower_bound = [xmin, ymin],
    upper_bound = [xmax, ymax],
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

electrons = picmi.Species(
    particle_type='electron', name='electrons'
)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diag1',
    grid = grid,
    period = 10,
    data_list = ['phi'],
    write_dir = '.',
    warpx_file_prefix = f'Python_restart_runtime_components_plt'
)

checkpoint = picmi.Checkpoint(
    name = 'chkpoint',
    period = 5,
    write_dir = '.',
    warpx_file_prefix = f'Python_restart_runtime_components_chk'
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
        n_macroparticle_per_cell=[0, 0], grid=grid
    )
)

for arg in sys.argv:
    if arg.startswith("amr.restart"):
        restart_file_name = arg.split("=")[1]
        sim.amr_restart = restart_file_name

sim.add_diagnostic(field_diag)
sim.add_diagnostic(checkpoint)
sim.initialize_inputs()
sim.initialize_warpx()

##########################
# python particle data access
##########################

from pywarpx import _libwarpx, callbacks

_libwarpx.add_real_comp('electrons', 'newPid')

def add_particles():

    nps = 10
    x = np.random.rand(nps) * 0.03
    y = np.zeros(nps)
    z = np.random.random(nps) * 0.03
    ux = np.random.normal(loc=0, scale=1e3, size=nps)
    uy = np.random.normal(loc=0, scale=1e3, size=nps)
    uz = np.random.normal(loc=0, scale=1e3, size=nps)
    w = np.ones(nps) * 2.0
    newPid = 5.0

    _libwarpx.add_particles(
        species_name='electrons', x=x, y=y, z=z, ux=ux, uy=uy, uz=uz,
        w=w, newPid=newPid
    )

callbacks.installbeforestep(add_particles)

##########################
# simulation run
##########################

step_number = _libwarpx.libwarpx.warpx_getistep(0)
sim.step(max_steps - 1 - step_number)

##########################
# check that the new PIDs are properly set
##########################

assert(_libwarpx.get_particle_count('electrons') == 90)
assert (_libwarpx.get_particle_comp_index('electrons', 'w') == 0)
assert (_libwarpx.get_particle_comp_index('electrons', 'newPid') == 4)

new_pid_vals = _libwarpx.get_particle_arrays(
    'electrons', 'newPid', 0
)
for vals in new_pid_vals:
    assert np.allclose(vals, 5)

##########################
# take the final sim step
##########################

sim.step(1)
