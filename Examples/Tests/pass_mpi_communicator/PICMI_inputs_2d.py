# --- Test if MPI_communicator is correctly passed from Python to CPP
# --- This test should be run with MPI enabled
# --- Inputs taken from langmuir_2d. Runs 1 step, and will check
# --- if the correct amount of processors are initialized in AMReX.

from mpi4py import MPI
from pywarpx import picmi
constants = picmi.constants
# from pywarpx import geometry
# # need to set geometry before importing _libwarpx
# geometry.coord_sys = 0
# geometry.prob_lo = [0]*2


##########################
# MPI communicator setup
##########################

# must initialize everything separately on separate processors
# so we need to create comm first before setting params,
# and importing packages

# split processor 0 into separate communicator from others
comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
if rank == 0:
    color = 0
else:
    color = 1
new_comm = comm_world.Split(color)

if color == 0:

    ##########################
    # physics parameters
    ##########################

    plasma_density = 1.e18
    plasma_xmin = 0.
    plasma_x_velocity = 0.1*constants.c

    ##########################
    # numerics parameters
    ##########################

    # --- Number of time steps
    max_steps = 10
    diagnostic_intervals = "::10"

    # --- Grid
    nx = 64
    ny = 64

    xmin = -20.e-6
    ymin = -20.e-6
    xmax = +20.e-6
    ymax = +20.e-6

    number_per_cell_each_dim = [2,2]

    ##########################
    # physics components
    ##########################

    uniform_plasma = picmi.UniformDistribution(density = plasma_density,
                                            upper_bound = [0., None, None],
                                            directed_velocity = [0.1*constants.c, 0., 0.])

    electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=uniform_plasma)

    ##########################
    # numerics components
    ##########################

    grid = picmi.Cartesian2DGrid(number_of_cells = [nx, ny],
                                lower_bound = [xmin, ymin],
                                upper_bound = [xmax, ymax],
                                lower_boundary_conditions = ['periodic', 'periodic'],
                                upper_boundary_conditions = ['periodic', 'periodic'],
                                moving_window_velocity = [0., 0., 0.],
                                warpx_max_grid_size = 32)

    solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

    ##########################
    # diagnostics
    ##########################

    field_diag1 = picmi.FieldDiagnostic(name = 'diag1',
                                        grid = grid,
                                        period = diagnostic_intervals,
                                        data_list = ['Ex', 'Jx'],
                                        write_dir = '.',
                                        warpx_file_prefix = 'Python_pass_mpi_comm_plt1_')

    part_diag1 = picmi.ParticleDiagnostic(name = 'diag1',
                                        period = diagnostic_intervals,
                                        species = [electrons],
                                        data_list = ['weighting', 'ux'])

    ##########################
    # simulation setup
    ##########################

    sim = picmi.Simulation(solver = solver,
                        max_steps = max_steps,
                        verbose = 1,
                        warpx_current_deposition_algo = 'direct')

    sim.add_species(electrons,
                    layout = picmi.GriddedLayout(n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid))

    sim.add_diagnostic(field_diag1)
    sim.add_diagnostic(part_diag1)

    ##########################
    # simulation run and test
    ##########################

    # NOTE: Some tests are done in this input PICMI file. These tests
    # check that amrex initialized the correct amount of procs and
    # that the procs ranks in amrex are correct.
    # If any of these tests fail, the terminal will print that the
    # program crashed.

    # pass communicators to amrex
    # Manually initialize inputs and warpx, to pass mpi communicator
    sim.step(max_steps, mpi_comm=new_comm)

    comm_world_size = comm_world.size
    new_comm_size = new_comm.size

    # only import wx after splitting up processors otherwise
    # throws a segmentation fault. Don't know why this happens.
    from pywarpx import wx

    # verify that communicator contains correct number of procs (1)
    assert wx.libwarpx.warpx_getNProcs() == comm_world_size - 2
    assert wx.libwarpx.warpx_getNProcs() == new_comm_size


else:

    ##########################
    # physics parameters
    ##########################

    plasma_density = 1.e19
    plasma_xmin = 0.
    plasma_x_velocity = 0.1*constants.c

    ##########################
    # numerics parameters
    ##########################

    # --- Number of time steps
    max_steps = 10
    diagnostic_intervals = "::10"

    # --- Grid
    nx = 64
    ny = 64

    xmin = -20.e-6
    ymin = -20.e-6
    xmax = +20.e-6
    ymax = +20.e-6

    number_per_cell_each_dim = [2,2]

    ##########################
    # physics components
    ##########################

    uniform_plasma = picmi.UniformDistribution(density = plasma_density,
                                            upper_bound = [0., None, None],
                                            directed_velocity = [0.1*constants.c, 0., 0.])

    electrons = picmi.Species(particle_type='electron', name='electrons', initial_distribution=uniform_plasma)

    ##########################
    # numerics components
    ##########################

    grid = picmi.Cartesian2DGrid(number_of_cells = [nx, ny],
                                lower_bound = [xmin, ymin],
                                upper_bound = [xmax, ymax],
                                lower_boundary_conditions = ['periodic', 'periodic'],
                                upper_boundary_conditions = ['periodic', 'periodic'],
                                moving_window_velocity = [0., 0., 0.],
                                warpx_max_grid_size = 32)

    solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.)

    ##########################
    # diagnostics
    ##########################

    field_diag2 = picmi.FieldDiagnostic(name = 'diag2',
                                        grid = grid,
                                        period = diagnostic_intervals,
                                        data_list = ['Ex', 'Jx'],
                                        write_dir = '.',
                                        warpx_file_prefix = 'Python_pass_mpi_comm_plt2_')

    part_diag2 = picmi.ParticleDiagnostic(name = 'diag2',
                                        period = diagnostic_intervals,
                                        species = [electrons],
                                        data_list = ['weighting', 'ux'])

    ##########################
    # simulation setup
    ##########################

    sim = picmi.Simulation(solver = solver,
                        max_steps = max_steps,
                        verbose = 1,
                        warpx_current_deposition_algo = 'direct')

    sim.add_species(electrons,
                    layout = picmi.GriddedLayout(n_macroparticle_per_cell=number_per_cell_each_dim, grid=grid))

    sim.add_diagnostic(field_diag2)
    sim.add_diagnostic(part_diag2)

    ##########################
    # simulation run and test
    ##########################

    # NOTE: Some tests are done in this input PICMI file to verify
    # that the communicators were passed correctly to the
    # , in this section
    # the analysis file is only there to satisfy the test script
    # requirements. The test will crash if any asserts fail.

    # pass communicators to amrex
    # Manually initialize inputs and warpx, to pass mpi communicator
    sim.step(max_steps, mpi_comm=new_comm)


    comm_world_size = comm_world.size
    new_comm_size = new_comm.size

    # only import wx after splitting up processors otherwise
    # throws a segmentation fault. Don't know why this happens.
    from pywarpx import wx

    # verify that amrex initialized with 1 fewer proc than comm world
    assert wx.libwarpx.warpx_getNProcs() == comm_world_size - 1
    assert wx.libwarpx.warpx_getNProcs() == new_comm_size

    # verify that amrex proc ranks are offset by -1 from
    # world comm proc ranks
    assert wx.libwarpx.warpx_getMyProc() == rank - 1
