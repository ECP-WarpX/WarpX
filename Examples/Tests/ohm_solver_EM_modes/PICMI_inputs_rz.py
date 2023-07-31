#!/usr/bin/env python3
#
# --- Test script for the kinetic-fluid hybrid model in WarpX wherein ions are
# --- treated as kinetic particles and electrons as an isothermal, inertialess
# --- background fluid. The script is set up to produce parallel normal EM modes
# --- in a metallic cylinder and is run in RZ geometry.
# --- See Section xxxx for a discussion and theoretical solution for this test.
# --- As a CI test only a small number of steps are taken.

import argparse
import os
import sys

import dill
from mpi4py import MPI as mpi
import numpy as np

from pywarpx import callbacks, fields, picmi

constants = picmi.constants

comm = mpi.COMM_WORLD

simulation = picmi.Simulation(verbose=0)
# make a shorthand for simulation.extension since we use it a lot
sim_ext = simulation.extension


class CylindricalNormalModes(object):
    '''The following runs a simulation of an uniform plasma at a set ion
    temperature (Te = 0) with an external magnetic field applied in the
    z-direction (parallel to domain).
    The analysis script (in this same directory) analyzes the output field
    data for EM modes.
    '''
    # Applied field parameters
    B0          = 0.5 # Initial magnetic field strength (T)
    beta        = 0.05 # Plasma beta, used to calculate temperature

    # Plasma species parameters
    m_ion       = 400.0 # Ion mass (electron masses)
    vA_over_c   = 1e-3 # ratio of Alfven speed and the speed of light

    # Spatial domain
    Nz          = 8 # 512 # 1024 # number of cells in z direction
    Nr          = 128 # 64 # 256 # number of cells in r direction

    # Temporal domain (if not run as a CI test)
    LT          = 100 # 500.0 # Simulation temporal length (ion cyclotron periods)

    # Numerical parameters
    NPPC        = 4*1024 # Seed number of particles per cell
    DZ          = 1.0 / 4.0 # Cell size (ion skin depths)
    DT          = 1e-6 # 5e-3 # Time step (ion cyclotron periods)

    # Plasma resistivity - used to dampen the mode excitation
    eta = 1e-7
    # Number of substeps used to update B
    substeps = 250

    def __init__(self, test, verbose):
        """Get input parameters for the specific case desired."""
        self.test = test
        self.verbose = verbose or self.test

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        self.dz = self.DZ * self.l_i
        self.Lz = self.Nz * self.dz
        self.Lr = self.Nr * self.dz

        self.dt = self.DT * self.t_ci

        if not self.test:
            self.total_steps = int(self.LT / self.DT)
        else:
            # if this is a test case run for only a small number of steps
            self.total_steps = 250
        # output diagnostics 20 times per cyclotron period
        self.diag_steps = int(1.0/20 / self.DT)

        # dump all the current attributes to a dill pickle file
        if comm.rank == 0:
            with open(f'sim_parameters.dpkl', 'wb') as f:
                dill.dump(self, f)

        # print out plasma parameters
        if comm.rank == 0:
            print(
                f"Initializing simulation with input parameters:\n"
                f"\tT = {self.T_plasma:.3f} eV\n"
                f"\tn = {self.n_plasma:.1e} m^-3\n"
                f"\tB0 = {self.B0:.2f} T\n"
                f"\tM/m = {self.m_ion:.0f}\n"
            )
            print(
                f"Plasma parameters:\n"
                f"\tl_i = {self.l_i:.1e} m\n"
                f"\tt_ci = {self.t_ci:.1e} s\n"
                f"\tv_ti = {self.v_ti:.1e} m/s\n"
                f"\tvA = {self.vA:.1e} m/s\n"
            )
            print(
                f"Numerical parameters:\n"
                f"\tdz = {self.dz:.1e} m\n"
                f"\tdt = {self.dt:.1e} s\n"
                f"\tdiag steps = {self.diag_steps:d}\n"
                f"\ttotal steps = {self.total_steps:d}\n"
            )

        self.setup_run()

    def get_plasma_quantities(self):
        """Calculate various plasma parameters based on the simulation input."""
        # Ion mass (kg)
        self.M = self.m_ion * constants.m_e

        # Cyclotron angular frequency (rad/s) and period (s)
        self.w_ci = constants.q_e * abs(self.B0) / self.M
        self.t_ci = 2.0 * np.pi / self.w_ci

        # Alfven speed (m/s): vA = B / sqrt(mu0 * n * (M + m)) = c * omega_ci / w_pi
        self.vA = self.vA_over_c * constants.c
        self.n_plasma = (
            (self.B0 / self.vA)**2 / (constants.mu0 * (self.M + constants.m_e))
        )

        # Ion plasma frequency (Hz)
        self.w_pi = np.sqrt(
            constants.q_e**2 * self.n_plasma / (self.M * constants.ep0)
        )

        # Skin depth (m)
        self.l_i = constants.c / self.w_pi

        # Ion thermal velocity (m/s) from beta = 2 * (v_ti / vA)**2
        self.v_ti = np.sqrt(self.beta / 2.0) * self.vA

        # Temperature (eV) from thermal speed: v_ti = sqrt(kT / M)
        self.T_plasma = self.v_ti**2 * self.M / constants.q_e # eV

        # Larmor radius (m)
        self.rho_i = self.v_ti / self.w_ci

    def setup_run(self):
        """Setup simulation components."""

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        self.grid = picmi.CylindricalGrid(
            number_of_cells=[self.Nr, self.Nz],
            warpx_max_grid_size=self.Nr, # self.Nz,
            lower_bound=[0, 0],
            upper_bound=[self.Lr, self.Lz],
            lower_boundary_conditions = ['none', 'periodic'],
            upper_boundary_conditions = ['dirichlet', 'periodic'],
            lower_boundary_conditions_particles = ['absorbing', 'periodic'],
            # upper_boundary_conditions_particles = ['absorbing', 'periodic'],
            upper_boundary_conditions_particles = ['reflecting', 'periodic']
        )
        simulation.time_step_size = self.dt
        simulation.max_steps = self.total_steps
        simulation.current_deposition_algo = 'direct'
        simulation.particle_shape = 1
        simulation.verbose = self.verbose

        #######################################################################
        # Field solver and external field                                     #
        #######################################################################

        self.solver = picmi.HybridPICSolver(
            grid=self.grid,
            Te=0.0, n0=self.n_plasma, plasma_resistivity=self.eta,
            substeps=2,#self.substeps,
            n_floor=self.n_plasma*0.05
            # n_floor=self.n_plasma*1.1 # check if boundary handling causes the problems
        )
        simulation.solver = self.solver

        B_ext = picmi.AnalyticInitialField(
            Bz_expression=self.B0
        )
        simulation.add_applied_field(B_ext)

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        self.ions = picmi.Species(
            name='ions', charge='q_e', mass=self.M,
            initial_distribution=picmi.UniformDistribution(
                density=self.n_plasma,
                rms_velocity=[self.v_ti]*3,
            )
        )
        simulation.add_species(
            self.ions,
            layout=picmi.PseudoRandomLayout(
                grid=self.grid, n_macroparticles_per_cell=0.0 #self.NPPC
            )
        )
        # use custom particle loader to avoid numerical Ji_r initially
        # due to macro-particle diffusion
        callbacks.installparticleloader(self._load_particles)

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        field_diag = picmi.FieldDiagnostic(
            name='field_diag',
            grid=self.grid,
            period=(self.total_steps if self.test else self.diag_steps),
            data_list=['B', 'E', 'J', 'rho'],
            # data_list=(['B', 'E'] if not self.test else ['B', 'E', 'j', 'rho']),
            # write_dir=('.' if self.test else 'diags'),
            write_dir='.',
            warpx_file_prefix=(
                'Python_ohms_law_solver_EM_modes_1d_plt' if self.test
                else 'diags'
            ),
            warpx_format=('openpmd' if not self.test else None),
            warpx_openpmd_backend=('h5' if not self.test else None),
            warpx_write_species=(True if self.test else False)
        )
        simulation.add_diagnostic(field_diag)

        callbacks.installafterstep(self.inspect_fields)

        #######################################################################
        # Initialize simulation                                               #
        #######################################################################

        simulation.initialize_inputs()
        simulation.initialize_warpx()

    def _load_particles(self):

        # determine how many particles to inject in total
        parts_to_inject = self.Nr * self.Nz * self.NPPC

        # determine how many particles to inject per process
        my_parts_to_inject = parts_to_inject // sim_ext.getNProcs()

        # determine particle weights
        weight = (
            self.n_plasma * np.pi * self.Lr**2 * self.Lz
            / (my_parts_to_inject * sim_ext.getNProcs())
        )

        # generate uniformly random positions
        rp = np.sqrt(np.random.rand(my_parts_to_inject)) * self.Lr
        tp = np.random.rand(my_parts_to_inject) * 2.0 * np.pi

        xp = rp * np.cos(tp)
        yp = rp * np.sin(tp)
        zp = np.random.rand(my_parts_to_inject) * self.Lz

        uxp = self.v_ti * np.random.randn(my_parts_to_inject)
        uyp = self.v_ti * np.random.randn(my_parts_to_inject)
        uzp = self.v_ti * np.random.randn(my_parts_to_inject)

        # add rigid rotation and radial velocity
        self.ur = 1.0
        self.ut = 1.0
        uxp += xp / rp * self.ur - yp / rp * self.ut
        uyp += xp / rp * self.ut + yp / rp * self.ur

        # add axial velocity
        self.uz = 1.0
        uzp += self.uz

        sim_ext.add_particles(
            self.ions.name, xp, yp, zp, uxp, uyp, uzp, weight
        )
        print("Particles injected.")

    def inspect_fields(self):
        # if sim_ext.getistep() % 20:
        #     return

        # nodal = np.arange(-simulation.particle_shape-1, self.Nr+2+simulation.particle_shape)
        nodal = np.arange(-simulation.particle_shape-1, self.Nr+2+simulation.particle_shape)
        cell_centered = nodal[:-1] + 0.5

        import matplotlib.pyplot as plt

        # Jir = fields.JxWrapper(include_ghosts=True)[...]
        # Jir = np.mean(
        #     Jir[:,simulation.particle_shape+1:-simulation.particle_shape-1],
        #     axis=1
        # )
        # plt.plot(cell_centered[1:-1], Jir / (self.n_plasma * constants.q_e * self.ur), 'o--')
        # plt.plot(cell_centered, Jir[:,3] / (self.n_plasma * constants.q_e * self.ur), 'o--')
        # plt.plot(cell_centered, Jir[:,5] / (self.n_plasma * constants.q_e * self.ur), 'o--')

        Jit = fields.JyWrapper(include_ghosts=True)[...]
        Jit = np.mean(
            Jit[:,simulation.particle_shape+1:-simulation.particle_shape-1],
            axis=1
        )
        plt.plot(nodal[1:-1], Jit / (self.n_plasma * constants.q_e * self.ut), 'o--')

        # Jiz = fields.JzWrapper(include_ghosts=True)[...]
        # Jiz = np.mean(
        #     Jiz[:,simulation.particle_shape+1:-simulation.particle_shape-1],
        #     axis=1
        # )
        # plt.plot(nodal[1:-1], Jiz / (self.n_plasma * constants.q_e * self.uz), 'o--')

        # rho = fields.RhoFPWrapper(include_ghosts=True)[...]
        # print(rho.shape)
        # rho = np.mean(
        #     rho[:,simulation.particle_shape+1:-simulation.particle_shape-1],
        #     axis=1
        # )
        # plt.plot(nodal, rho/constants.q_e/self.n_plasma, 'o--')


        # Jr = fields.JxFPAmpereWrapper(include_ghosts=True)[...]
        # # Jir = np.mean(
        # #     Jir[:,simulation.particle_shape+1:-simulation.particle_shape-1],
        # #     axis=1
        # # )
        # plt.plot(cell_centered[1:-1], Jr[:,3], 'o--')
        # plt.plot(cell_centered[1:-1], Jr[:,5], 'o--')
        # print(Jir.shape)

        # Bt = fields.ByWrapper(include_ghosts=True)[...]
        # # Jir = np.mean(
        # #     Jir[:,simulation.particle_shape+1:-simulation.particle_shape-1],
        # #     axis=1
        # # )
        # plt.plot(cell_centered, Bt[:,3], 'o--')
        # plt.plot(cell_centered, Bt[:,5], 'o--')

        # plt.plot(np.var(Jir, axis=1), 'o--')
        plt.show()




##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-t', '--test', help='toggle whether this script is run as a short CI test',
    action='store_true',
)
parser.add_argument(
    '-v', '--verbose', help='Verbose output', action='store_true',
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

run = CylindricalNormalModes(test=args.test, verbose=args.verbose)
simulation.step(1)
