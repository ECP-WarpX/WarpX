#!/usr/bin/env python3
import numpy as np
from pywarpx import picmi

# Physical constants
c = picmi.constants.c

max_steps = 200

# Create grid
grid = picmi.Cartesian3DGrid(
    number_of_cells=[32, 32, 32],
    lower_bound=[-1., -1., -1.],
    upper_bound=[1., 1., 1.],
    lower_boundary_conditions=['periodic', 'periodic','periodic'],
    upper_boundary_conditions=['periodic', 'periodic','periodic'],
    lower_boundary_conditions_particles=['periodic', 'periodic','periodic'],
    upper_boundary_conditions_particles=['periodic', 'periodic','periodic'],
    warpx_max_grid_size=32)

electrons = picmi.Species(
        particle_type='electron',
        name=f'electrons',
        initial_distribution= picmi.UniformDistribution(
                                density=1.e14,
                                rms_velocity = [0.035*c]*3
        ))
photons = picmi.Species(
        particle_type='photon',
        name='photons',
        initial_distribution=picmi.UniformDistribution(
                                density=1.e14,
                                rms_velocity = [0.2*c]*3)
)

protons = picmi.Species(
        particle_type='proton',
        name=f'protons',
        initial_distribution=picmi.UniformDistribution(
                                density=1.e14
        ))

species_list = [electrons, protons, photons]


# reduced diagnostics
EP = picmi.ReducedDiagnostic(name='EP', type='ParticleEnergy',intervals=200)
NP = picmi.ReducedDiagnostic(name='NP', type='ParticleNumber',intervals=200)
EF = picmi.ReducedDiagnostic(name='EF', type='FieldEnergy',intervals=200)
PP = picmi.ReducedDiagnostic(name='PP', type='ParticleMomentum',intervals=200)
PF = picmi.ReducedDiagnostic(name='PF', type='FieldMomentum',intervals=200)
MF = picmi.ReducedDiagnostic(name='MF', type='FieldMaximum',intervals=200)
MR = picmi.ReducedDiagnostic(name='MR', type='RhoMaximum',intervals=200)
FP = picmi.ReducedDiagnostic(name='FP',
                                type='FieldProbe',
                                intervals=200,
                                #warpx_probe_geometry='Point',
                                warpx_x_probe = 0.53125,
                                warpx_y_probe = 0.53125,
                                warpx_z_probe = 0.53125)
FP_integrate = picmi.ReducedDiagnostic(name='FP_integrate',
                                type='FieldProbe',
                                intervals=20,
                                warpx_probe_geometry='Point',
                                warpx_x_probe = 0.53125,
                                warpx_y_probe = 0.53125,
                                warpx_z_probe = 0.53125,
                                warpx_integrate=1)
FP_line = picmi.ReducedDiagnostic(name='FP_line',
                                type='FieldProbe',
                                intervals=200,
                                warpx_probe_geometry='Line',
                                warpx_x_probe = 0.53125,
                                warpx_y_probe = 0.53125,
                                warpx_z_probe = 0.53125,
                                warpx_x1_probe = 0.70225,
                                warpx_y1_probe = 0.70225,
                                warpx_z1_probe = 0.70225,
                                warpx_resolution = 100)
FP_plane = picmi.ReducedDiagnostic(name='FP_plane',
                                type='FieldProbe',
                                intervals=200,
                                warpx_probe_geometry='Plane',
                                warpx_x_probe = 0.5,
                                warpx_y_probe = 0.5,
                                warpx_z_probe = 0.5,
                                warpx_target_normal_x = 0,
                                warpx_target_normal_y = 0,
                                warpx_target_normal_z = 1,
                                warpx_target_up_x = 0,
                                warpx_target_up_y = 1,
                                warpx_target_up_z = 0,
                                warpx_detector_radius = 0.25,
                                warpx_resolution = 10)
FR_Max = picmi.ReducedDiagnostic(name='FR_Max',
                                type='FieldReduction',
                                intervals=200,
                                warpx_reduction_type='Maximum',
                                warpx_reduced_function='sqrt(Bx**2 + By**2 + Bz**2)'
                                )
FR_Min = picmi.ReducedDiagnostic(name='FR_Min',
                                type='FieldReduction',
                                intervals=200,
                                warpx_reduction_type='Minimum',
                                warpx_reduced_function='x*Ey*Bz'
                                )
FR_Integral = picmi.ReducedDiagnostic(name='FR_Integral',
                                type='FieldReduction',
                                intervals=200,
                                warpx_reduction_type='Integral',
                                warpx_reduced_function='"if(y > 0 and z < 0, '
                                    + '0.5*((Ex**2 + Ey**2 + Ez**2)*epsilon0+(Bx**2 + By**2 + Bz**2)/mu0), 0)"'
                                )

reduced_diags = [EP, NP, EF, PP, PF, MF, MR, FP, FP_integrate, FP_line, FP_plane, FR_Max,FR_Min,FR_Integral]

smoother = picmi.BinomialSmoother()
# Electromagnetic solver
solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',
    cfl=0.9999,
    source_smoother=smoother
    )

# Diagnostics
diag_field_list = ['B', 'E', 'J', 'rho', 'rho_electrons', 'rho_protons']
field_diag = picmi.FieldDiagnostic(
    name='diag1',
    grid=grid,
    period=200,
    data_list=diag_field_list,
    write_dir='diags')

# Set up simulation
sim = picmi.Simulation(
    solver=solver,
    max_steps=max_steps,
    # max_time=max_time,
    verbose=2,
    particle_shape='linear',
    warpx_charge_deposition_algo='standard',
    warpx_current_deposition_algo='esirkepov',
    warpx_field_gathering_algo='energy-conserving',
    warpx_particle_pusher_algo='vay')

for species in species_list:
    sim.add_species(
        species,
        layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[1,1,1]))

# Add diagnostics
sim.add_diagnostic(field_diag)
for rd in reduced_diags:
    sim.add_reduced_diagnostic(rd)

sim.initialize_inputs()
sim.initialize_warpx()

sim.write_input_file('inputs_PICMI_reduced')

# Advance simulation until last time step
sim.step(max_steps)
