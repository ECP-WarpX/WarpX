from pywarpx import picmi

##########################
# physics parameters
##########################

V_domain_boundary = 0.0
V_embedded_boundary = 1.0

##########################
# numerics parameters
##########################

dt = 1e-6

# --- Nb time steps

max_steps = 1

# --- grid

nx = 64
ny = 64
nz = 64

xmin = -0.5
xmax = 0.5
ymin = -0.5
ymax = 0.5
zmin = -0.5
zmax = 0.5


##########################
# numerics components
##########################

grid = picmi.Cartesian3DGrid(
    number_of_cells = [nx, ny, nz],
    lower_bound = [xmin, ymin, zmin],
    upper_bound = [xmax, ymax, zmax],
    lower_boundary_conditions = ['dirichlet', 'dirichlet', 'dirichlet'],
    upper_boundary_conditions = ['dirichlet', 'dirichlet', 'dirichlet'],
    lower_boundary_conditions_particles = ['absorbing', 'absorbing', 'absorbing'],
    upper_boundary_conditions_particles = ['absorbing', 'absorbing', 'absorbing'],
    warpx_potential_lo_x = V_domain_boundary,
    warpx_potential_hi_x = V_domain_boundary,
    warpx_potential_lo_y = V_domain_boundary,
    warpx_potential_hi_y = V_domain_boundary,
    warpx_potential_lo_z = V_domain_boundary,
    warpx_potential_hi_z = V_domain_boundary,
    warpx_blocking_factor=8,
    warpx_max_grid_size = 128
)

solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid', required_precision=1e-7
)

embedded_boundary = picmi.EmbeddedBoundary(
    implicit_function="-(x**2+y**2+z**2-0.3**2)",
    potential=V_embedded_boundary
)

##########################
# diagnostics
##########################

field_diag = picmi.FieldDiagnostic(
    name = 'diag1',
    grid = grid,
    period = 1,
    data_list = ['Ex', 'Ey', 'Ez', 'phi', 'rho'],
    write_dir = '.',
    warpx_file_prefix = 'Python_ElectrostaticSphereEB_plt'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver = solver,
    time_step_size = dt,
    max_steps = max_steps,
    warpx_embedded_boundary=embedded_boundary,
    warpx_field_gathering_algo='momentum-conserving'
)

sim.add_diagnostic(field_diag)

##########################
# simulation run
##########################

sim.step(max_steps)
