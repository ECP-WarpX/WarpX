#!/usr/bin/env python3
#
# --- Test script for the semi-implicit Poisson solver in which a plasma column
# --- expands inside a conducting cylinder.

import argparse
import sys

import numpy as np
from mpi4py import MPI as mpi

from pywarpx import picmi

constants = picmi.constants

comm = mpi.COMM_WORLD

simulation = picmi.Simulation(warpx_serialize_initial_conditions=True, verbose=0)


class PlasmaExpansionSimulation(object):
    """Simulation setup for an expanding plasma."""

    # Plasma parameters
    m_ion = 250.0  # Ion mass (electron masses)
    n_plasma = 1e18  # Plasma density (m^-3)
    T_e = 2.0  # Electron temperature (eV)
    T_i = 1.0  # Ion temperature (eV)

    scale_fac = 4.0  # scaling factor used to move between explicit and SI

    # Plasma size
    plasma_radius = 20  # Plasma column radius (Debye lengths)

    # Spatial domain
    NZ = 256 / scale_fac  # Number of cells in x and z directions

    # Temporal domain (if not run as a CI test)
    LT = 0.1  # Simulation temporal length in ion crossing times

    # Numerical parameters
    NPPC = 750  # Seed number of particles per cell
    DZ = 1.0 * scale_fac  # Cell size (Debye lengths)
    DT = 0.75  # Time step (electron CFL)

    # Solver parameter
    C_SI = 2  # Semi-implicit factor - marginal stability threshold = 1

    def __init__(self, test, verbose):
        """Get input parameters for the specific case desired."""
        self.test = test
        self.verbose = verbose or self.test

        # calculate various plasma parameters based on the simulation input
        self.get_plasma_quantities()

        self.dz = self.DZ * self.lambda_e
        self.Lz = self.NZ * self.dz

        self.dt = self.DT * self.dz / self.v_te

        self.plasma_radius *= self.lambda_e

        self.total_steps = int(self.LT * self.Lz / (self.v_ti * self.dt))
        self.diag_steps = max(50, self.total_steps // 10)

        # print out plasma parameters
        if comm.rank == 0:
            print(
                f"Initializing simulation with input parameters:\n"
                f"\tT_e = {self.T_e:.1f} eV\n"
                f"\tT_i = {self.T_i:.1f} eV\n"
                f"\tn = {self.n_plasma:.1e} m^-3\n"
                f"\tM/m = {self.m_ion:.0f}\n"
            )
            print(
                f"Plasma parameters:\n"
                f"\tlambda_e = {self.lambda_e:.1e} m\n"
                f"\tt_pe = {1.0/self.f_pe:.1e} s\n"
                f"\tv_ti = {self.v_ti:.1e} m/s\n"
            )
            print(
                f"Numerical parameters:\n"
                f"\tdz = {self.dz:.1e} m\n"
                f"\tdt*f_pe = {self.dt*self.f_pe:.2f}\n"
                f"\tdiag steps = {self.diag_steps:d}\n"
                f"\ttotal steps = {self.total_steps:d}\n"
            )
        self.setup_run()

    def get_plasma_quantities(self):
        """Calculate various plasma parameters based on the simulation input."""
        # Ion mass (kg)
        self.M = self.m_ion * constants.m_e

        # Electron plasma frequency (Hz)
        self.f_pe = np.sqrt(
            constants.q_e**2 * self.n_plasma / (constants.m_e * constants.ep0)
        ) / (2.0 * np.pi)

        # Debye length (m)
        self.lambda_e = np.sqrt(
            constants.ep0 * self.T_e / (self.n_plasma * constants.q_e)
        )

        # Thermal velocities (m/s) from v_th = np.sqrt(kT / m)
        self.v_ti = np.sqrt(self.T_i * constants.q_e / self.M)
        self.v_te = np.sqrt(self.T_e * constants.q_e / constants.m_e)

    def setup_run(self):
        """Setup simulation components."""

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        self.grid = picmi.Cartesian2DGrid(
            number_of_cells=[self.NZ, self.NZ],
            lower_bound=[-self.Lz / 2.0, -self.Lz / 2.0],
            upper_bound=[self.Lz / 2.0, self.Lz / 2.0],
            lower_boundary_conditions=["dirichlet"] * 2,
            upper_boundary_conditions=["dirichlet"] * 2,
            lower_boundary_conditions_particles=["absorbing"] * 2,
            upper_boundary_conditions_particles=["absorbing"] * 2,
            warpx_max_grid_size=self.NZ,
        )
        simulation.time_step_size = self.dt
        simulation.max_steps = self.total_steps
        simulation.current_deposition_algo = "direct"
        simulation.particle_shape = 1
        simulation.verbose = self.verbose

        #######################################################################
        # Insert probe as embedded boundary                                   #
        #######################################################################

        embedded_boundary = picmi.EmbeddedBoundary(
            implicit_function=f"(x**2+z**2-{self.Lz/2.-8*self.lambda_e}**2)",
            potential=0.0,
        )
        simulation.embedded_boundary = embedded_boundary

        #######################################################################
        # Field solver and external field                                     #
        #######################################################################

        solver = picmi.ElectrostaticSolver(
            grid=self.grid,
            method="Multigrid",
            warpx_semi_implicit=True,
            warpx_semi_implicit_factor=self.C_SI,
            warpx_self_fields_verbosity=0,
        )
        simulation.solver = solver

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        density_expr = f"if(x**2+z**2<{self.plasma_radius**2},{self.n_plasma},0)"
        self.electrons = picmi.Species(
            name="electrons",
            particle_type="electron",
            initial_distribution=picmi.AnalyticDistribution(
                density_expression=density_expr,
                warpx_momentum_spread_expressions=[self.v_te] * 3,
                warpx_density_min=self.n_plasma / 10.0,
            ),
        )
        simulation.add_species(
            self.electrons,
            layout=picmi.PseudoRandomLayout(
                grid=self.grid, n_macroparticles_per_cell=self.NPPC
            ),
        )

        self.ions = picmi.Species(
            name="ions",
            charge="q_e",
            mass=self.M,
            initial_distribution=picmi.AnalyticDistribution(
                density_expression=density_expr,
                warpx_momentum_spread_expressions=[self.v_ti] * 3,
                warpx_density_min=self.n_plasma / 10.0,
            ),
        )
        simulation.add_species(
            self.ions,
            layout=picmi.PseudoRandomLayout(
                grid=self.grid, n_macroparticles_per_cell=self.NPPC
            ),
        )

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        field_diag = picmi.FieldDiagnostic(
            name="field_diag",
            grid=self.grid,
            period=self.diag_steps,
            data_list=[
                "E",
                "J",
                "T_electrons",
                "T_ions",
                "phi",
                "rho_electrons",
                "rho_ions",
            ],
        )
        simulation.add_diagnostic(field_diag)

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
    "-t",
    "--test",
    help="toggle whether this script is run as a short CI test",
    action="store_true",
)
parser.add_argument(
    "-v",
    "--verbose",
    help="Verbose output",
    action="store_true",
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1] + left

run = PlasmaExpansionSimulation(test=args.test, verbose=args.verbose)
simulation.step()
