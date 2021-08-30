# --- Input file to test the particle scraper and the Python wrappers
# --- to access the buffer of scraped particles.

from pywarpx import picmi

##########################
# numerics parameters
##########################

# --- Number of time steps
max_steps = 60
diagnostic_intervals = 20

# --- Grid
nx = 64
ny = 64
nz = 128

cfl = 0.99

xmin = -125e-6
ymin = -125e-6
zmin = -149e-6
xmax = 125e-6
ymax = 125e-6
zmax = 1e-6

##########################
# physics components
##########################

uniform_plasma_elec = picmi.UniformDistribution(
    density = 1e23, # number of electrons per m^3
    lower_bound = [-1e-5, -1e-5, -149e-6],
    upper_bound = [1e-5, 1e-5, -129e-6],
    directed_velocity = [0., 0., 2000.*picmi.constants.c] # uth the std of the (unitless) momentum
)

electrons = picmi.Species(
    particle_type='electron', name='electrons',
    initial_distribution=uniform_plasma_elec,
    warpx_save_particles_at_xhi=1, warpx_save_particles_at_eb=1
)

##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(
    number_of_cells = [nx, ny, nz],
    lower_bound = [xmin, ymin, zmin],
    upper_bound = [xmax, ymax, zmax],
    lower_boundary_conditions=['none', 'none', 'none'],
    upper_boundary_conditions=['none', 'none', 'none'],
    lower_boundary_conditions_particles=['open', 'open', 'open'],
    upper_boundary_conditions_particles=['open', 'open', 'open'],
    warpx_max_grid_size = 128
)

solver = picmi.ElectromagneticSolver(
    grid=grid, cfl=cfl
)

embedded_boundary = picmi.EmbeddedBoundary(
    implicit_function="-max(max(max(x-12.5e-6,-12.5e-6-x),max(y-12.5e-6,-12.5e-6-y)),max(z-(-6.15e-5),-8.65e-5-z))"
)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diag1',
    grid = grid,
    period = diagnostic_intervals,
    data_list = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'],
    write_dir = '.',
    warpx_file_prefix = 'Python_particle_scrape_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    max_steps = max_steps,
    warpx_embedded_boundary=embedded_boundary,
    verbose=True
)

sim.add_species(
    electrons,
    layout = picmi.GriddedLayout(
        n_macroparticle_per_cell=[1, 1, 1], grid=grid
    )
)

sim.add_diagnostic(field_diag)

##########################
# simulation run
##########################

# sim.write_input_file(file_name = 'inputs_from_PICMI')
sim.step(max_steps)

################################################
# check that the wrappers to access the particle
# buffer functions as intended
################################################

from pywarpx import _libwarpx

n = _libwarpx.get_particle_boundary_buffer_size("electrons", 'eb')
print("Number of electrons in buffer:", n)
assert n == 612

scraped_steps = _libwarpx.get_particle_boundary_buffer("electrons", 'eb', 'step_scraped', 0)
for arr in scraped_steps:
    assert all(arr > 40)

weights = _libwarpx.get_particle_boundary_buffer("electrons", 'eb', 'w', 0)
assert sum(len(arr) for arr in weights) == 612

# clear the particle buffer
_libwarpx.libwarpx.warpx_clearParticleBoundaryBuffer()
# confirm that the buffer was cleared
n = _libwarpx.get_particle_boundary_buffer_size("electrons", 'eb')
print("Number of electrons in buffer:", n)
assert n == 0
