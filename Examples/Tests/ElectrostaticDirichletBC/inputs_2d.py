from pywarpx import picmi

grid = picmi.Cartesian2DGrid(
    number_of_cells = [64, 8],
    lower_bound = [0, 0],
    upper_bound = [0.032, 0.004],
    bc_xmin = 'pec',
    bc_xmax = 'pec',
    bc_ymin = 'periodic',
    bc_ymax = 'periodic',
    warpx_potential_lo_x = "150.0*sin(2*pi*6.78e6*t)",
    warpx_potential_hi_x = "450.0*sin(2*pi*13.56e6*t)",
    bc_xmin_particles = 'absorbing',
    bc_xmax_particles = 'absorbing',
    bc_ymin_particles = 'periodic',
    bc_ymax_particles = 'periodic',
    moving_window_velocity = None,
    warpx_max_grid_size = 32
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-6
)

field_diag = picmi.FieldDiagnostic(
    name='dirichletbc_plt',
    grid=grid,
    period=4,
    data_list=['phi'],
    write_dir='.',
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
