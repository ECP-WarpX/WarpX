from pywarpx import picmi
#from warp import picmi

# Note: Very low resolution is used in this test case since it's main
# purpose is testing the operation of the code.

constants = picmi.constants

nx = 32
ny = 32
nz = 16

xmin = -100.e-6
xmax = +100.e-6
ymin = -100.e-6
ymax = +100.e-6
zmin = -100.e-6
zmax = +0.e-6

moving_window_velocity = [0., 0., constants.c]

number_per_cell_each_dim = [1, 1, 1]

grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                             lower_bound = [xmin, ymin, zmin],
                             upper_bound = [xmax, ymax, zmax],
                             lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                             upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                             moving_window_velocity = moving_window_velocity,
                             warpx_max_grid_size=32)

solver = picmi.ElectromagneticSolver(grid=grid, cfl=1)

beam_distribution = picmi.UniformDistribution(density = 1.e23,
                                              lower_bound = [-20.e-6, -20.e-6, -75.e-6],
                                              upper_bound = [+20.e-6, +20.e-6, -25.e-6],
                                              directed_velocity = [0., 0., 1.e9])

plasma_distribution = picmi.UniformDistribution(density = 1.e22,
                                                lower_bound = [None, None, 0.],
                                                upper_bound = [None, None, None],
                                                fill_in = True)

beam = picmi.Species(particle_type='electron', name='beam', initial_distribution=beam_distribution)
plasma = picmi.Species(particle_type='electron', name='plasma', initial_distribution=plasma_distribution)

sim = picmi.Simulation(solver = solver,
                       max_steps = 10,
                       verbose = 1,
                       warpx_plot_int = 10,
                       warpx_current_deposition_algo = 'esirkepov')

sim.add_species(beam, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))
sim.add_species(plasma, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
#sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()

