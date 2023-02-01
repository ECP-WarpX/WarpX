#!/usr/bin/env python3
#
# --- Test script for the kinetic-fluid hybrid model in WarpX wherein ions are
# --- treated as kinetic particles and electrons as an isothermal, inertialess
# --- background fluid. The script is set up to produce either parallel or
# --- perpendicular (Bernstein) EM modes and can be run in 1d, 2d or 3d
# --- Cartesian geometries. As a CI test only a small number of steps are taken
# --- using the 1d version.

import argparse
import dill
import sys

import numpy as np
from mpi4py import MPI as mpi

from pywarpx import picmi

constants = picmi.constants

comm = mpi.COMM_WORLD

simulation = picmi.Simulation(verbose=0)
# make a shorthand for simulation.extension since we use it a lot
sim_ext = simulation.extension


class EMModes(object):
    '''The following runs a simulation of an uniform plasma at a set
    temperature (Te = Ti) with an external magnetic field applied in either the
    z-direction (parallel to domain) or x-direction (perpendicular to domain).
    The analysis script (in this same directory) analyzes the output field data
    for EM modes. This input is based on the EM modes tests as described by
    Munoz et al. (2018) and tests done by Scott Nicks at TAE Technologies.
    '''

    # Applied field parameters
    B0          = 0.25 # Initial magnetic field strength (T)
    beta        = [0.01, 0.1] # Plasma beta, used to calculate temperature

    # Plasma species parameters
    m_ion       = [100.0, 400.0] # Ion mass (electron masses)
    vA_over_c  = [1e-4, 1e-3] # ratio of Alfven speed and the speed of light

    # Spatial domain
    Nz          = [1024, 1920] # number of cells in z direction
    Nx          = 8 # number of cells in x (and y) direction for >1 dimensions

    # Temporal domain (if not run as a CI test)
    LT          = 200.0 # 600.0 # Simulation temporal length (ion cyclotron periods)

    # Numerical parameters
    NPPC        = [1024, 512, 128] # Seed number of particles per cell
    DZ          = [1.0/10.0, 1.0/30.0] # Cell size (ion skin depths)
    DT          = [5e-3, 4e-3] # Time step (ion cyclotron periods)

    # Plasma resistivity - used to dampen the parallel mode excitation
    eta = [1e-7, 0.0]
    # Number of substeps used to update B
    substeps = 200

    def __init__(self, test, dim, B_dir):
        """Get input parameters for the specific case desired."""
        self.test = test
        self.dim = int(dim)
        self.B_dir = B_dir

        # sanity check
        assert (dim > 0 and dim < 4), f"{dim}-dimensions not a valid input"

        # get simulation parameters from the defaults given the direction of
        # the initial B-field and the dimensionality
        self.get_simulation_parameters()

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        self.dz = self.DZ * self.l_i
        self.Lz = self.Nz * self.dz
        self.Lx = self.Nx * self.dz

        self.dt = self.DT * self.t_ci
        # self.dt = self.DT * self.dz / self.vA
        # self.dt = self.CFL_FACTOR * self.dz / constants.c

        if not self.test:
            self.total_steps = int(np.ceil(self.LT * self.t_ci / self.dt))
            self.diag_steps = int(1/20 * self.t_ci / self.dt)
        else:
            # if this is a test case run for only 1000 steps
            self.total_steps = 100
            self.diag_steps = 50

        # dump all the current attributes to a dill pickle file
        if comm.rank == 0:
            with open('sim_parameters.dpkl', 'wb') as f:
                dill.dump(self, f)

        # print out plasma parameters
        if comm.rank == 0:
            print(
                f"Initializing simulation with input parameters:\n"
                f"\tT = {self.T_plasma*1e3:.1f} meV\n"
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

    def get_simulation_parameters(self):
        """Pick appropriate parameters from the defaults given the direction
        of the B-field and the simulation dimensionality."""
        if self.B_dir == 'z':
            idx = 0
            self.Bx = 0.0
            self.Bz = self.B0
        else:
            idx = 1
            self.Bx = self.B0
            self.Bz = 0.0

        self.beta = self.beta[idx]
        self.m_ion = self.m_ion[idx]
        self.vA_over_c = self.vA_over_c[idx]
        self.Nz = self.Nz[idx]
        self.NPPC = self.NPPC[self.dim-1]
        self.DZ = self.DZ[idx]
        self.DT = self.DT[idx]
        self.eta = self.eta[idx]

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
        self.l_i = 2.0 * np.pi * constants.c / self.w_pi

        # Ion thermal velocity (m/s) from beta = 2 * (v_ti / vA)**2
        self.v_ti = np.sqrt(self.beta / 2.0) * self.vA

        # Temperature (eV) from thermal speed: v_ti = sqrt(kT / M)
        self.T_plasma = self.v_ti**2 * self.M / constants.q_e # eV

        # Larmor radius (m)
        self.rho_i = self.v_ti / self.w_ci

        # Coulomb_log = 15.0
        # nu = (
        #     constants.q_e**4 * self.n_plasma * Coulomb_log
        #     / (2.0*np.pi*constants.ep0)**2
        #     / np.sqrt(constants.m_e * (self.T_plasma * constants.q_e)**3)
        # )
        # print(nu)
        # print(4.18e-12*self.n_plasma*Coulomb_log/(self.T_plasma**1.5))
        # self.eta = constants.m_e * nu / (constants.q_e**2 * self.n_plasma)
        # print(self.eta)

        # print("spitzer: ",
        #     4.0 * np.sqrt(2.0*np.pi*constants.m_e) / 3.0 * constants.q_e**2 * Coulomb_log
        #     / ((4.0 * np.pi * constants.ep0)**2*(self.T_plasma*constants.q_e)**1.5)
        #     / 1.96
        # )
        # exit()

    def setup_run(self):
        """Setup simulation components."""

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        if self.dim == 1:
            grid_object = picmi.Cartesian1DGrid
        elif self.dim == 2:
            grid_object = picmi.Cartesian2DGrid
        elif self.dim == 3:
            grid_object = picmi.Cartesian3DGrid

        self.grid = grid_object(
            number_of_cells=[self.Nx, self.Nx, self.Nz][-self.dim:],
            warpx_max_grid_size=self.Nz,
            lower_bound=[-self.Lx/2.0, -self.Lx/2.0, 0][-self.dim:],
            upper_bound=[self.Lx/2.0, self.Lx/2.0, self.Lz][-self.dim:],
            lower_boundary_conditions=['periodic']*self.dim,
            upper_boundary_conditions=['periodic']*self.dim
        )
        simulation.time_step_size = self.dt
        simulation.max_steps = self.total_steps
        simulation.load_balance_intervals = self.total_steps // 200
        simulation.verbose = True #self.test

        #######################################################################
        # Field solver and external field                                     #
        #######################################################################

        self.solver = picmi.EMSolver(
            grid=self.grid, method='hybrid',
            Te=self.T_plasma, n0=self.n_plasma, plasma_resistivity=self.eta,
            substeps=self.substeps
        )
        simulation.solver = self.solver

        B_ext = picmi.AnalyticInitialField(
            Bx_expression=self.Bx,
            By_expression=0.0,
            Bz_expression=self.Bz
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
                grid=self.grid, n_macroparticles_per_cell=self.NPPC
            )
        )

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        # field_diag = picmi.FieldDiagnostic(
        #     name='field_diag',
        #     grid=self.grid,
        #     period=self.diag_steps,
        #     data_list=['B', 'E'],
        #     # write_dir=('.' if self.test else 'diags'),
        #     # warpx_file_prefix='Python_hybrid_PIC_plt',
        #     warpx_format = 'openpmd',
        #     warpx_openpmd_backend = 'h5'
        # )
        # simulation.add_diagnostic(field_diag)

        line_diag = picmi.ReducedDiagnostic(
            diag_type='FieldProbe',
            probe_geometry='Line',
            z_probe=0,
            z1_probe=self.Lz,
            resolution=self.Nz - 1,
            name='lineprobe',
            period=self.diag_steps,
            path='diags/'
        )
        simulation.add_diagnostic(line_diag)

        #######################################################################
        # Initialize simulation                                               #
        #######################################################################

        # simulation.current_deposition_algo = 'direct'

        simulation.initialize_inputs()
        simulation.initialize_warpx()


##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-t', '--test', help='toggle whether this script is run as a short CI test',
    action='store_true',
)
parser.add_argument(
    '-d', '--dim', help='Simulation dimension', required=False, type=int,
    default=1
)
parser.add_argument(
    '--bdir', help='Direction of the B-field', required=False,
    choices=['x', 'z'], default='z'
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

run = EMModes(test=args.test, dim=args.dim, B_dir=args.bdir)
simulation.step()