#!/usr/bin/env python3

from pywarpx import picmi

c = picmi.constants.c
q_e = picmi.constants.q_e
m_e = picmi.constants.m_e
m_p = picmi.constants.m_p
ep0 = picmi.constants.ep0

nx = 64
nz = 64

xmin = -100.e-6
xmax = +100.e-6
zmin = -149.e-6
zmax = +1.e-6

moving_window_velocity = [0., 0., c]
gamma_boost = 10.0

npc_x = 1
npc_z = 1

max_steps = 10

field_bc_z='open'
field_bc_x='open'
grid = picmi.Cartesian3DGrid(
    number_of_cells = [nx, nx, nz],
    lower_bound = [xmin, xmin, zmin],
    upper_bound = [xmax, xmax, zmax],
    lower_boundary_conditions = [field_bc_x, field_bc_x, field_bc_z],
    upper_boundary_conditions = [field_bc_x, field_bc_x, field_bc_z],
    lower_boundary_conditions_particles = ['absorbing']*3,
    upper_boundary_conditions_particles = ['absorbing']*3,
    moving_window_velocity = moving_window_velocity,
    warpx_blocking_factor = 32,
    warpx_max_grid_size=128
)

smoother = picmi.BinomialSmoother(n_pass=[1]*3)
solver = picmi.ElectromagneticSolver(
    method='CKC',
    grid=grid,
    cfl=0.99,
    source_smoother=smoother,
    warpx_pml_ncell=10)

drive_charge = -10e-11
driver_gaussian_distribution = picmi.GaussianBunchDistribution(
    n_physical_particles= abs(drive_charge) / q_e,
    rms_bunch_size=[2e-6, 2e-6, 4.e-6],
    rms_velocity=[2.*c, 2.*c, 2.e4*c],
    centroid_position=[0., 0., -20.e-6],
    centroid_velocity=[0., 0., 2.e5*c],
)
beam_charge = -3.e-13
beam_z0 = -90e-6
beam_gaussian_distribution = picmi.GaussianBunchDistribution(
    n_physical_particles= abs(beam_charge) / q_e,
    rms_bunch_size=[0.5e-6, 0.5e-6, 1.e-6],
    rms_velocity=[2.*c, 2.*c, 200.*c],
    centroid_position=[0., 0., beam_z0],
    centroid_velocity=[0., 0., 4000.*c],
)

dens = 1.e23
lramp = 8.e-3
plasma_distribution = picmi.AnalyticDistribution(
    density_expression=
        f'(z<lramp)*0.5*(1-cos(pi*z/lramp))*dens+(z>lramp)*dens',
    pi=3.141592653589793,
    dens=dens,
    lramp = lramp,
    lower_bound=[-70.e-6, -70.e-6, 0.],
    upper_bound=[70.e-6, 70.e-6, 0.2],
    fill_in=True)

beam = picmi.Species(
    particle_type='electron',
    name='beam',
    initial_distribution=beam_gaussian_distribution)
driver = picmi.Species(
    particle_type='electron',
    name='driver',
    initial_distribution=driver_gaussian_distribution)
plasma_e = picmi.Species(
    particle_type='electron',
    name='plasma_e',
    initial_distribution=plasma_distribution)
plasma_i = picmi.Species(
    particle_type='proton',
    name='plasma_p',
    initial_distribution=plasma_distribution)

sim = picmi.Simulation(
    solver = solver,
    max_steps = max_steps,
    verbose = 1,
    warpx_use_filter = 0,
    gamma_boost = gamma_boost,
    warpx_use_fdtd_nci_corr = True,
    particle_shape='cubic'
)

sim.add_species_through_plane(
    beam,
    layout=picmi.PseudoRandomLayout(
        grid=grid,
        n_macroparticles=1000
    ),
    injection_plane_position=0.,
    injection_plane_normal_vector=[0.,0.,1.]
)
sim.add_species_through_plane(
    driver,
    layout=picmi.PseudoRandomLayout(
        grid=grid,
        n_macroparticles=1000
    ),
    injection_plane_position=0.,
    injection_plane_normal_vector=[0.,0.,1.]
)
sim.add_species(
    plasma_e,
    layout=picmi.GriddedLayout(
        grid=grid,
        n_macroparticle_per_cell=[npc_x, npc_x, npc_z]
    )
)
sim.add_species(
    plasma_i,
    layout=picmi.GriddedLayout(
        grid=grid,
        n_macroparticle_per_cell=[npc_x, npc_x, npc_z]
    )
)

field_diag = picmi.FieldDiagnostic(
    name = 'Python_PlasmaAcceleration_plt',
    grid = grid,
    period = max_steps,
    data_list = ['Ex', 'Ey', 'Ez', 'Jx', 'Jy', 'Jz', 'part_per_cell'],
    # write_dir='Python_PlasmaAcceleration_plt',
    warpx_format='openpmd',
)

part_diag = picmi.ParticleDiagnostic(
    name = 'Python_PlasmaAcceleration_plt',
    period = max_steps,
    species = [driver, beam, plasma_e],
    # write_dir='Python_PlasmaAcceleration_plt',
    warpx_format='openpmd',
)

beamrel_red_diag = picmi.ReducedDiagnostic(
    diag_type='BeamRelevant',
    name='beamrel',
    species=beam,
    period=10)

sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)
sim.add_diagnostic(beamrel_red_diag)

# write_inputs will create an inputs file that can be used to run
# with the compiled version.
sim.write_input_file(file_name = 'inputs_from_PICMI')

# Alternatively, sim.step will run WarpX, controlling it from Python
sim.step()
