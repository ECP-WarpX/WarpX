"""
Monte-Carlo Collision script based on Turner et al. (2013)
https://doi.org/10.1063/1.4775084
"""

import argparse
from functools import partial
import logging
import sys

import numpy as np
from pywarpx import callbacks, picmi
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as sla

from mewarpx import emission, mcc_wrapper, mespecies
from mewarpx.diags_store import diag_base
from mewarpx.diags_store.field_diagnostic import FieldDiagnostic
from mewarpx.mwxrun import mwxrun

constants = picmi.constants

logger = logging.getLogger(__name__)


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

    def initialize_inputs(self):
        """Grab geometrical quantities from the grid. The boundary potentials
        are also obtained from the grid using 'warpx_potential_zmin' for the
        left_voltage and 'warpx_potential_zmax' for the right_voltage.
        These can be given as floats or strings that can be parsed by the
        WarpX parser.
        """
        # grab the boundary potentials from the grid object
        left_voltage = self.grid.potential_zmin
        if left_voltage is None:
            left_voltage = 0.0
        self.left_voltage = partial(mwxrun.eval_expression_t, left_voltage)

        right_voltage = self.grid.potential_zmax
        if right_voltage is None:
            right_voltage = 0.0
        self.right_voltage = partial(mwxrun.eval_expression_t, right_voltage)

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

        logger.info("Using direct solver.")
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
        if not mwxrun.initialized:
            return

        # get rho from WarpX
        self.rho_data = mwxrun.get_gathered_rho_grid()[:,0]
        # run superLU solver to get phi
        self.solve()
        # write phi to WarpX
        mwxrun.set_phi_grid(self.phi)

    def solve(self):
        """The solution step. Includes getting the boundary potentials and
        calculating phi from rho."""

        left_voltage = self.left_voltage()
        right_voltage = self.right_voltage()

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

    PLASMA_DENSITY = [2.56e14, 5.12e14, 5.12e14, 3.84e14] # m^-3
    T_ELEC = 30000.0 # K

    SEED_NPPC = 16 * np.array([32, 16, 8, 4])

    NZ = [128, 256, 512, 512]

    DT = 1.0 / (np.array([400, 800, 1600, 3200]) * FREQ)

    # Total simulation time in seconds
    TOTAL_TIME = np.array([1280, 5120, 5120, 15360]) / FREQ
    # Time (in seconds) between diagnostic evaluations
    DIAG_INTERVAL = 32 / FREQ

    def setup_run(self, n=0, steps=None):
        """Setup run for the specific case (n) desired."""

        self.steps = steps

        # Case specific input parameters
        self.ANODE_VOLTAGE = f"{self.VOLTAGE[n]}*sin(2*pi*{self.FREQ:.5e}*t)"

        self.N_INERT = self.N_INERT[n]
        self.PLASMA_DENSITY = self.PLASMA_DENSITY[n]
        self.SEED_NPPC = self.SEED_NPPC[n]

        self.NZ = self.NZ[n]

        self.DT = self.DT[n]
        self.MAX_STEPS = int(self.TOTAL_TIME[n] / self.DT)
        self.DIAG_STEPS = int(self.DIAG_INTERVAL / self.DT)

        if self.steps is not None:
            self.MAX_STEPS = self.steps
            self.DIAG_STEPS = self.MAX_STEPS // 5

        #######################################################################
        # Set geometry, boundary conditions and timestep                      #
        #######################################################################

        mwxrun.init_grid(
            lower_bound=[0], upper_bound=[self.D_CA],
            number_of_cells=[self.NZ], min_tiles=16
        )
        mwxrun.grid.potential_zmax = self.ANODE_VOLTAGE
        mwxrun.init_timestep(DT=self.DT)
        mwxrun.simulation.max_steps = self.MAX_STEPS
        mwxrun.simulation.load_balance_intervals = self.MAX_STEPS // 5000

        #######################################################################
        # Field solver                                                        #
        #######################################################################

        # self.solver = picmi.ElectrostaticSolver(
        #    grid=mwxrun.grid, method='Multigrid', required_precision=1e-12,
        #    warpx_self_fields_verbosity=0
        # )
        self.solver = PoissonSolver1D(grid=mwxrun.grid)
        mwxrun.simulation.solver = self.solver

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        self.electrons = mespecies.Species(
            particle_type='electron', name='electrons'
        )
        self.ions = mespecies.Species(
            particle_type='He', name='he_ions', charge='q_e'
        )

        #######################################################################
        # Collision  initialization                                           #
        #######################################################################

        self.mcc = mcc_wrapper.MCC(
            self.electrons, self.ions, T_INERT=self.T_INERT,
            N_INERT=self.N_INERT, exclude_collisions=['charge_exchange']
        )

        #######################################################################
        # Neutral plasma injection                                            #
        #######################################################################

        self.vol_emitter = emission.UniformDistributionVolumeEmitter(
            T=self.T_ELEC
        )

        self.plasma_injector = emission.PlasmaInjector(
            emitter=self.vol_emitter, species1=self.electrons,
            species2=self.ions, npart=2 * self.SEED_NPPC * self.NZ,
            T_2=self.T_INERT, plasma_density=self.PLASMA_DENSITY
        )

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        self.text_diag = diag_base.TextDiag(
            diag_steps=self.DIAG_STEPS, preset_string='perfdebug'
        )
        '''
        self.field_diag = FieldDiagnostic(
            diag_steps=self.DIAG_STEPS, barrier_slices=[0],
            save_pdf=False, style='roelof', min_dim=2.0,
            plot=(self.steps is not None)
        )
        '''
        # array to save ion density
        self.rho_array = np.zeros(self.NZ + 1)

        #######################################################################
        # Initialize run and print diagnostic info                            #
        #######################################################################

        mwxrun.init_run()


    def _get_rho_ions(self):
        rho_data = mwxrun.get_gathered_rho_grid('he_ions', False)
        if mwxrun.me == 0:
            self.rho_array += rho_data[:,0] / constants.q_e / self.DIAG_STEPS

    def run_sim(self):

        if self.steps is not None:
            mwxrun.simulation.step(self.steps)
            self.text_diag.print_performance_summary()
        else:
            mwxrun.simulation.step(self.MAX_STEPS - self.DIAG_STEPS)
            callbacks.installafterstep(self._get_rho_ions)
            mwxrun.simulation.step(self.DIAG_STEPS)

            if mwxrun.me == 0:
                np.save('avg_rho_data.npy', self.rho_array)

##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', '--steps', help='set the number of simulation steps manually',
    type=int
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

if args.steps:
    steps = args.steps
else:
    steps = None

run = CapacitiveDischargeExample()
run.setup_run(n=CASE_NUM, steps=steps)
run.run_sim()
