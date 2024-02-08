#!/usr/bin/env python3
import argparse
import sys

import numpy as np

from pywarpx import callbacks, libwarpx, particle_containers, picmi

# Create the parser and add the argument
parser = argparse.ArgumentParser()
parser.add_argument(
    '-u', '--unique', action='store_true',
    help="Whether injected particles should be treated as unique"
)

# Parse the input
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1] + left

##########################
# numerics parameters
##########################

dt = 1.0e-11

# --- Nb time steps

max_steps = 100
diagnostic_interval = 10

# --- grid

nr = 64
nz= 64

rmin = 0.0
rmax = 2
zmin = -2
zmax = 2

##########################
# numerics components
##########################

grid = picmi.CylindricalGrid(
    number_of_cells = [nr, nz],
    n_azimuthal_modes = 1,
    lower_bound = [rmin, zmin],
    upper_bound = [rmax, zmax],
    lower_boundary_conditions = ['none', 'dirichlet'],
    upper_boundary_conditions =  ['dirichlet', 'dirichlet'],
    lower_boundary_conditions_particles = ['absorbing', 'reflecting'],
    upper_boundary_conditions_particles =  ['absorbing', 'reflecting']
)


solver = picmi.ElectrostaticSolver(
    grid=grid, method='Multigrid',
    warpx_absolute_tolerance=1e-7
)

embedded_boundary = picmi.EmbeddedBoundary(
    implicit_function="-(x**2+y**2+z**2-radius**2)",
    radius = 0.3277
)

##########################
# physics components
##########################


n = 7.0e9 #plasma density #particles/m^3
Te = 85 #Electron temp in eV
qe = picmi.constants.q_e #elementary charge
m_e = picmi.constants.m_e #electron mass
v_eth = (qe * Te / m_e) ** 0.5

# nothing to change in the distribution function?
#e_dist = picmi.UniformDistribution(density = n, rms_velocity=[v_eth, v_eth, v_eth] )
e_dist=picmi.ParticleListDistribution(x=0.0, y=0.0, z=-0.4, ux=0.5e10, uy=0.0, uz=1.0e10, weight=1)

electrons = picmi.Species(
    particle_type='electron', name='electrons', initial_distribution=e_dist, warpx_save_particles_at_eb=1
)

##########################
# diagnostics
##########################

part_diag = picmi.ParticleDiagnostic(name = 'diag1',
    period = diagnostic_interval,
    species = [electrons],
    warpx_format = 'openpmd'
)

##########################
# simulation setup
##########################

sim = picmi.Simulation(
    solver=solver,
    time_step_size = dt,
    max_steps = max_steps,
    warpx_embedded_boundary=embedded_boundary,
    warpx_amrex_the_arena_is_managed=1,
)

sim.add_species(
    electrons,
    layout = picmi.GriddedLayout(
        n_macroparticle_per_cell=[10, 1, 1], grid=grid
    )
)
sim.add_diagnostic(part_diag)


sim.initialize_inputs()
sim.initialize_warpx()


##########################
# python particle data access
##########################

#angle=
#damping_coef=
#number=


def mirror_reflection():
    #global nb_steps_scraping
    buffer = particle_containers.ParticleBoundaryBufferWrapper()
    if (len(buffer.get_particle_boundary_buffer("electrons", 'eb', 'timestamp', 0))!=0): #otherwise np.concatenate doesnt work
        print('the particle hitted the EB')
        time_steps = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'timestamp', 0))

        #step 1: extract the different parameters of the scraping buffer (normal, time, position)
        x = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'x', 0))
        y = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'theta', 0))
        z = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'z', 0))
        ux = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'ux', 0))
        uy = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'uy', 0))
        uz = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'uz', 0))
        w = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'w', 0))
        nx = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'nx', 0))
        ny = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'ny', 0))
        nz = np.concatenate(buffer.get_particle_boundary_buffer("electrons", 'eb', 'nz', 0))

        #step 2: use these parameters to inject from the same position electrons in the plasma
        elect_pc = particle_containers.ParticleContainerWrapper('electrons')

        un=ux*nx+uy*ny+uz*nz
        ux_reflect=2*un*nx-ux
        uy_reflect=2*un*ny-uy
        uz_reflect=2*un*nz-uz

        elect_pc.add_particles(
            #for a "mirror reflection" u(sym)=2(u.n)n-u
            x=x, y=y, z=z, ux=ux_reflect, uy=uy_reflect, uz=uz_reflect,
            w=w,
            unique_particles=args.unique
            )

        buffer.clear_buffer()


callbacks.installafterstep(mirror_reflection)

##########################
# simulation run
##########################

sim.step(max_steps)
