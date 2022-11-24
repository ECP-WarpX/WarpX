#!/usr/bin/env python3
#
# --- Test script for the kinetic-fluid hybrid model in WarpX wherein ions are
# --- treated as kinetic particles and electrons as an isothermal background
# --- fluid. The script is set up to produce ion Bernstein modes and can be run
# --- in 1d, 2d or 3d Cartesian geometries. As a CI test only a small number of
# --- steps are taken using the 1d version.

import argparse
import dill
import sys

import numpy as np
from mpi4py import MPI as mpi

from pywarpx import callbacks, fields, picmi

constants = picmi.constants

comm = mpi.COMM_WORLD

simulation = picmi.Simulation(verbose=0)
# make a shorthand for simulation.extension since we use it a lot
sim_ext = simulation.extension


class HybridPICExample(object):
    '''The following runs a simulation of an uniform plasma at a set
    temperature (Te = Ti) with an external magnetic field applied in the
    x-direction. The analysis script (in this same directory) analyzes the
    output field data for ion-Bernstein modes.
    This example is based on work done by Scott Nicks at TAE Technologies.
    '''

    # Applied field parameters
    B0          = 0.25 # Initial magnetic field strength (T)

    # Plasma species parameters
    m_ion       = 1000.0 # Ion mass (electron masses)
    n_plasma    = 1e20 # Plasma density (m^-3)
    T_plasma    = 1e3 # Electron and ion temperature (eV)

    # Spatial domain
    LZ          = 80.0 # Simulation domain length (ion skin depths)
    Nx          = 16 # number of cells in x (and y) direction for >1 dimensions

    # Temporal domain (if not run as a CI test)
    LT          = 40.0 # Simulation temporal length (ion cyclotron periods)

    # Numerical parameters
    NPPC        = [100, 50, 25] # Seed number of particles per cell
    DZ          = 0.05 # Cell size (ion skin depths)
    DT          = 0.1 # Time step (dz / Alfven speed)

    def __init__(self, test=False, dim=1):
        """Get input parameters for the specific case desired."""
        self.test = test
        self.dim = int(dim)

        # sanity check
        assert (dim > 0 and dim < 4), f"{dim}-dimensions not a valid input"

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        Nz_ideal = self.LZ / self.DZ
        self.Nz = int(2**(np.ceil(np.log2(Nz_ideal))))

        self.dz = self.DZ * self.l_i
        self.Lz = self.Nz * self.dz
        self.Lx = self.Nx * self.dz

        self.dt = self.DT * self.dz / self.vA

        if not self.test:
            self.total_steps = int(np.ceil(self.LT * self.t_ci / self.dt))
            diag_period = 1 / 32 # Output interval (ion cyclotron periods)
            self.diag_steps = int(diag_period * self.t_ci / self.dt)
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
                f"\tT = {self.T_plasma*1e-3:.1f} keV\n"
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

        # Plasma frequency (Hz)
        self.w_pi = np.sqrt(
            constants.q_e**2 * self.n_plasma / (self.M * constants.ep0)
        ) # / (2.0 * np.pi)

        # Skin depth
        self.l_i = constants.c / self.w_pi

        # Cyclotron frequency (Hz) and period (s)
        self.f_ci = constants.q_e * abs(self.B0) / (2.0 * np.pi * self.M)
        self.t_ci = 1.0 / self.f_ci

        # Thermal speed
        self.v_ti = np.sqrt(self.T_plasma * constants.q_e / self.M)

        # Larmor radius (m)
        self.rho_i = self.v_ti / self.f_ci

        # Alfven speed (m/s)
        self.vA = abs(self.B0) / np.sqrt(
            constants.mu0 * self.n_plasma * (constants.m_e + self.M)
        )

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
            warpx_max_grid_size=128,
            lower_bound=[-self.Lx/2.0, -self.Lx/2.0, 0][-self.dim:],
            upper_bound=[self.Lx/2.0, self.Lx/2.0, self.Lz][-self.dim:],
            lower_boundary_conditions=['periodic']*self.dim,
            upper_boundary_conditions=['periodic']*self.dim
        )
        simulation.time_step_size = self.dt
        simulation.max_steps = self.total_steps
        simulation.load_balance_intervals = self.total_steps // 1000
        simulation.verbose = self.test

        #######################################################################
        # Field solver and external field                                     #
        #######################################################################

        self.solver = picmi.EMSolver(
            grid=self.grid, method='hybrid',
            Te=self.T_plasma, n0=self.n_plasma, plasma_resistivity=0.0
        )
        simulation.solver = self.solver

        B_ext = picmi.AnalyticInitialField(
            Bx_expression=self.B0,
            By_expression=0.0,
            Bz_expression=0.0
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
                grid=self.grid,
                n_macroparticles_per_cell=self.NPPC[self.dim-1]
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
            x_probe=0, y_probe=0, z_probe=0,
            x1_probe=0, y1_probe=0, z1_probe=self.Lz,
            resolution=self.Nz - 1,
            name='lineprobe',
            period=self.diag_steps,
            path='diags/'
        )
        simulation.add_diagnostic(line_diag)

        #######################################################################
        # Initialize simulation                                               #
        #######################################################################

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
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

run = HybridPICExample(test=args.test, dim=args.dim)
simulation.step()