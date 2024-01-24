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

print('code debute')

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

electrons = picmi.Species(
    particle_type='electron', name='electrons'
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
    #step 1: extract the different parameters of the scraping buffer (normal, time, position)
    buffer = particle_containers.ParticleBoundaryBufferWrapper()
    #this doesn't work because we are not sure that at each step particles have been scraped --> we have to ensure that just before
    time = buffer.get_particle_boundary_buffer("electrons", 'eb', 'step_scraped', 0)[-1] #real
    x = buffer.get_particle_boundary_buffer("electrons", 'eb', 'x', 0)[-1] #array of reals
    y = buffer.get_particle_boundary_buffer("electrons", 'eb', 'y', 0)[-1] #array of reals
    z = buffer.get_particle_boundary_buffer("electrons", 'eb', 'z', 0)[-1] #array of reals
    ux = buffer.get_particle_boundary_buffer("electrons", 'eb', 'ux', 0)[-1] #array of reals
    uy = buffer.get_particle_boundary_buffer("electrons", 'eb', 'uy', 0)[-1] #array of reals
    uz = buffer.get_particle_boundary_buffer("electrons", 'eb', 'uz', 0)[-1] #array of reals
    w = buffer.get_particle_boundary_buffer("electrons", 'eb', 'w', 0)[-1] #array of reals
    nx = buffer.get_particle_boundary_buffer("electrons", 'eb', 'nx', 0)[-1] #array of reals
    ny = buffer.get_particle_boundary_buffer("electrons", 'eb', 'ny', 0)[-1] #array of reals
    nz = buffer.get_particle_boundary_buffer("electrons", 'eb', 'nz', 0)[-1] #array of reals


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
    

    

callbacks.installafterstep(mirror_reflection)

##########################
# simulation run
##########################

sim.step(max_steps)
