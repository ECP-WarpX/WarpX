from pywarpx import picmi

grid = picmi.Cartesian2DGrid(
    number_of_cells = [64, 8],
    lower_bound = [0, 0],
    upper_bound = [0.032, 0.004],
    lower_boundary_conditions = ['pec', 'periodic'],
    upper_boundary_conditions = ['pec', 'periodic'],
    lower_boundary_conditions_particles = ['absorbing', 'periodic'],
    upper_boundary_conditions_particles = ['absorbing', 'periodic'],
    warpx_potential_lo_x = "150.0*sin(2*pi*6.78e6*t)",
    warpx_potential_hi_x = "450.0*sin(2*pi*13.56e6*t)",
    moving_window_velocity = None,
    warpx_max_grid_size = 32
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-6
)

field_diag = picmi.FieldDiagnostic(
    name='diag1',
    grid=grid,
    period=4,
    data_list=['phi'],
    write_dir='.',
    warpx_file_prefix='dirichletbc_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    time_step_size=7.5e-10,
    max_steps=100,
    verbose=0
)

sim.add_diagnostic(field_diag)
sim.step()
