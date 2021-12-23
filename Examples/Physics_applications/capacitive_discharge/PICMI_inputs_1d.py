#!/usr/bin/env python3
#
# --- Copyright 2021 Modern Electron
# --- Monte-Carlo Collision script based on Turner et al. (2013)
# --- https://doi.org/10.1063/1.4775084

import argparse
from functools import partial
import sys

import numpy as np
from pywarpx import callbacks, picmi, fields, _libwarpx
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as sla

constants = picmi.constants


class PoissonSolver1D(picmi.ElectrostaticSolver):

    def __init__(self, grid, **kwargs):
        """Direct solver for the Poisson equation using superLU. This solver is
        useful for 1D cases.

        Arguments:
            grid (picmi.Cartesian1DGrid): Instance of the grid on which the
            solver will be installed.
        """
        # Sanity check that this solver is appropriate to use
        if not isinstance(grid, picmi.Cartesian1DGrid):
            raise RuntimeError('Direct solver can only be used on a 1D grid.')

        super(PoissonSolver1D, self).__init__(
            grid=grid, method=kwargs.pop('method', 'Multigrid'),
            required_precision=1, **kwargs
        )
        self.rho_wrapper = None
        self.phi_wrapper = None

    def initialize_inputs(self):
        """Grab geometrical quantities from the grid.
        """
        # grab the boundary potentials from the grid object
        self.right_voltage = self.grid.potential_zmax

        # set WarpX boundary potentials to None since we will handle it
        # ourselves in this solver
        self.grid.potential_xmin = None
        self.grid.potential_xmax = None
        self.grid.potential_ymin = None
        self.grid.potential_ymax = None
        self.grid.potential_zmin = None
        self.grid.potential_zmax = None

        super(PoissonSolver1D, self).initialize_inputs()

        self.nz = self.grid.nx
        self.dz = (self.grid.xmax - self.grid.xmin) / self.nz

        self.nxguardphi = 1
        self.nzguardphi = 1

        self.phi = np.zeros(self.nz + 1 + 2*self.nzguardphi)

        self.decompose_matrix()

        callbacks.installpoissonsolver(self._run_solve)

    def decompose_matrix(self):
        """Function to build the superLU object used to solve the linear
        system."""
        self.nsolve = self.nz + 1

        # Set up the tridiagonal computation matrix in order to solve A*phi =
        # rho for phi.
        self.A_ldiag = np.ones(self.nsolve-1) / self.dz**2
        self.A_mdiag = -2.*np.ones(self.nsolve) / self.dz**2
        self.A_udiag = np.ones(self.nsolve-1) / self.dz**2

        self.A_mdiag[0] = 1.
        self.A_udiag[0] = 0.0

        self.A_mdiag[-1] = 1.
        self.A_ldiag[-1] = 0.0

        # Set up the computation matrix in order to solve A*phi = rho
        A = np.zeros((self.nsolve, self.nsolve))
        kk = 0
        for ii in range(self.nsolve):
            for jj in range(self.nsolve):

                if ii == jj:
                    A[ii, jj] = -2.0
                elif ii == jj - 1 or ii == jj + 1:
                    A[ii, jj] = 1.0

        A[0, 1] = 0.0
        A[-1, -2] = 0.0
        A[0, 0] = 1.0
        A[-1, -1] = 1.0

        A = csc_matrix(A, dtype=np.float32)
        self.lu = sla.splu(A)

    def _run_solve(self):
        """Function run on every step to perform the required steps to solve
        Poisson's equation."""

        # get rho from WarpX
        if self.rho_wrapper is None:
            self.rho_wrapper = fields.RhoFPWrapper(0, False)
        self.rho_data = self.rho_wrapper[Ellipsis][:,0]

        # run superLU solver to get phi
        self.solve()

        # write phi to WarpX
        if self.phi_wrapper is None:
            self.phi_wrapper = fields.PhiFPWrapper(0, True)
        self.phi_wrapper[Ellipsis] = self.phi

    def solve(self):
        """The solution step. Includes getting the boundary potentials and
        calculating phi from rho."""

        left_voltage = 0.0
        right_voltage = eval(
            self.right_voltage,
            {'t':_libwarpx.libwarpx.gett_new(0), 'sin':np.sin, 'pi':np.pi}
        )

        # Construct b vector
        rho = -self.rho_data / constants.ep0
        b = np.zeros(rho.shape[0], dtype=np.float32)
        b[:] = rho * self.dz**2

        b[0] = left_voltage
        b[-1] = right_voltage

        phi = self.lu.solve(b)

        self.phi[self.nzguardphi:-self.nzguardphi] = phi

        self.phi[:self.nzguardphi] = left_voltage
        self.phi[-self.nzguardphi:] = right_voltage


class CapacitiveDischargeExample(object):

    #######################################################################
    # Begin global user parameters                                        #
    #######################################################################

    D_CA = 0.067 # m

    FREQ = 13.56e6 # Hz
    VOLTAGE = [450.0, 200.0, 150.0, 120.0] # V

    N_INERT = [9.64e20, 32.1e20, 96.4e20, 321e20] # m^-3
    T_INERT = 300.0 # K
    M_ION = 6.67e-27 # kg

    PLASMA_DENSITY = [2.56e14, 5.12e14, 5.12e14, 3.84e14] # m^-3
    T_ELEC = 30000.0 # K

    SEED_NPPC = 16 * np.array([32, 16, 8, 4])

    NZ = [128, 256, 512, 512]

    DT = 1.0 / (np.array([400, 800, 1600, 3200]) * FREQ)

    # Total simulation time in seconds
    TOTAL_TIME = np.array([1280, 5120, 5120, 15360]) / FREQ
    # Time (in seconds) between diagnostic evaluations
    DIAG_INTERVAL = 32 / FREQ

    def setup_run(self, n=0, test=False):
        """Setup run for the specific case (n) desired."""

        # Case specific input parameters
        self.VOLTAGE = f"{self.VOLTAGE[n]}*sin(2*pi*{self.FREQ:.5e}*t)"

        self.N_INERT = self.N_INERT[n]
        self.PLASMA_DENSITY = self.PLASMA_DENSITY[n]
        self.SEED_NPPC = self.SEED_NPPC[n]

        self.NZ = self.NZ[n]

        self.DT = self.DT[n]
        self.MAX_STEPS = int(self.TOTAL_TIME[n] / self.DT)
        self.DIAG_STEPS = int(self.DIAG_INTERVAL / self.DT)

        if test:
            self.MAX_STEPS = 20
            self.DIAG_STEPS = 5

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        self.grid = picmi.Cartesian1DGrid(
            number_of_cells=[self.NZ],
            warpx_max_grid_size=128,
            lower_bound=[0],
            upper_bound=[self.D_CA],
            lower_boundary_conditions=['dirichlet'],
            upper_boundary_conditions=['dirichlet'],
            lower_boundary_conditions_particles=['absorbing'],
            upper_boundary_conditions_particles=['absorbing'],
            warpx_potential_hi_z=self.VOLTAGE,
        )

        #######################################################################
        # Field solver                                                        #
        #######################################################################

        # self.solver = picmi.ElectrostaticSolver(
        #    grid=mwxrun.grid, method='Multigrid', required_precision=1e-12,
        #    warpx_self_fields_verbosity=0
        # )
        self.solver = PoissonSolver1D(grid=self.grid)

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        self.electrons = picmi.Species(
            particle_type='electron', name='electrons',
            initial_distribution=picmi.UniformDistribution(
                density=self.PLASMA_DENSITY,
                rms_velocity=[np.sqrt(constants.kb * self.T_ELEC / constants.m_e)]*3,
            )
        )
        self.ions = picmi.Species(
            particle_type='He', name='he_ions',
            charge='q_e', mass=self.M_ION,
            initial_distribution=picmi.UniformDistribution(
                density=self.PLASMA_DENSITY,
                rms_velocity=[np.sqrt(constants.kb * self.T_INERT / self.M_ION)]*3,
            )
        )

        #######################################################################
        # Collision  initialization                                           #
        #######################################################################

        cross_sec_direc = '../../../../warpx-data/MCC_cross_sections/He/'
        cross_sec_direc = '/home/roelof/software/warpx-data/MCC_cross_sections/He/'
        mcc_electrons = picmi.MCCCollisions(
            name='coll_elec',
            species=self.electrons,
            background_density=self.N_INERT,
            background_temperature=self.T_INERT,
            background_mass=self.ions.mass,
            scattering_processes={
                'elastic' : {
                    'cross_section' : cross_sec_direc+'electron_scattering.dat'
                },
                'excitation1' : {
                    'cross_section': cross_sec_direc+'excitation_1.dat',
                    'energy' : 19.82
                },
                'excitation2' : {
                    'cross_section': cross_sec_direc+'excitation_2.dat',
                    'energy' : 20.61
                },
                'ionization' : {
                    'cross_section' : cross_sec_direc+'ionization.dat',
                    'energy' : 24.55,
                    'species' : self.ions
                },
            }
        )

        mcc_ions = picmi.MCCCollisions(
            name='coll_ion',
            species=self.ions,
            background_density=self.N_INERT,
            background_temperature=self.T_INERT,
            scattering_processes={
                'elastic' : {
                    'cross_section' : cross_sec_direc+'ion_scattering.dat'
                },
                'back' : {
                    'cross_section' : cross_sec_direc+'ion_back_scatter.dat'
                },
                # 'charge_exchange' : {
                #    'cross_section' : cross_sec_direc+'charge_exchange.dat'
                # }
            }
        )

        #######################################################################
        # Initialize simulation                                               #
        #######################################################################

        self.sim = picmi.Simulation(
            solver=self.solver,
            time_step_size=self.DT,
            max_steps=self.MAX_STEPS,
            warpx_collisions=[mcc_electrons, mcc_ions],
            warpx_load_balance_intervals=self.MAX_STEPS//5000,
            verbose=test
        )

        self.sim.add_species(
            self.electrons,
            layout = picmi.GriddedLayout(
                n_macroparticle_per_cell=[self.SEED_NPPC], grid=self.grid
            )
        )
        self.sim.add_species(
            self.ions,
            layout = picmi.GriddedLayout(
                n_macroparticle_per_cell=[self.SEED_NPPC], grid=self.grid
            )
        )

        #######################################################################
        # Add diagnostics for the CI test to be happy                         #
        #######################################################################

        field_diag = picmi.FieldDiagnostic(
            name='diag1',
            grid=self.grid,
            period=0,
            data_list=['rho_electrons', 'rho_he_ions'],
            write_dir='.',
            warpx_file_prefix='Python_background_mcc_1d_plt'
        )
        self.sim.add_diagnostic(field_diag)

        #######################################################################
        # Declare array to save ion density                                       #
        #######################################################################

        self.ion_density_array = np.zeros(self.NZ + 1)

    def _get_rho_ions(self):
        # deposit the ion density in rho_fp
        self.sim.extension.depositChargeDensity('he_ions', 0)
        rho_data = self.solver.rho_wrapper[Ellipsis][:,0]
        self.ion_density_array += rho_data / constants.q_e / self.DIAG_STEPS

    def run_sim(self):

        self.sim.step(self.MAX_STEPS - self.DIAG_STEPS)
        callbacks.installafterstep(self._get_rho_ions)
        self.sim.step(self.DIAG_STEPS)

        if self.sim.extension.getMyProc() == 0:
            np.save('avg_ion_density.npy', self.ion_density_array)

##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-t', '--test', help='toggle whether this script is run as a test',
    action='store_true',
)
parser.add_argument(
    '-n', help='Test number to run (1 to 4)', required=False, type=int,
    default=1
)
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

if args.n < 1 or args.n > 4:
    raise AttributeError('Test number must be an integer from 1 to 4.')
CASE_NUM = args.n - 1

run = CapacitiveDischargeExample()
run.setup_run(n=CASE_NUM, test=args.test)
run.run_sim()
