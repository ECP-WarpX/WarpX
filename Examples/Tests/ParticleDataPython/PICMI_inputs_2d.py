from mpi4py import MPI
from pywarpx import picmi
import numpy as np

##########################
# MPI communicator setup
##########################

# split processor 0 into separate communicator from others
comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
if rank == 0:
    color = 0
else:
    color = 1
new_comm = comm_world.Split(color)

##########################
# numerics parameters
##########################

dt = 7.5e-10

# --- Nb time steps

max_steps = 10

# --- grid

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
    warpx_file_prefix = f'Python_particle_attr_access_plt_{color}'
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
sim.add_diagnostic(field_diag)

sim.initialize_inputs()
sim.initialize_warpx(mpi_comm=new_comm)

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
        w=w, newPid=newPid, unique_particles=(not color)
    )

callbacks.installbeforestep(add_particles)

##########################
# simulation run
##########################

sim.step(max_steps - 1, mpi_comm=new_comm)

##########################
# check that the new PIDs are properly set
##########################

if color == 0:
    assert (_libwarpx.get_particle_count('electrons') == 90)
else:
    assert (_libwarpx.get_particle_count('electrons') == 90)
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
