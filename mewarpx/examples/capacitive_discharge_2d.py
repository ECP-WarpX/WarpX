"""
Monte-Carlo Collision script based on Turner et al. (2013)
https://doi.org/10.1063/1.4775084
"""

from mewarpx.utils_store import util as mwxutil

mwxutil.init_libwarpx(ndim=2, rz=False)

import argparse
import sys

import numpy as np
from pywarpx import callbacks, picmi

from mewarpx import emission, mcc_wrapper, mespecies, poisson_pseudo_1d
from mewarpx.diags_store import diag_base
from mewarpx.diags_store.field_diagnostic import FieldDiagnostic
from mewarpx.mwxrun import mwxrun

constants = picmi.constants


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

    NX = 8
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
            lower_bound=[-self.D_CA / self.NZ * self.NX / 2.0, 0],
            upper_bound=[self.D_CA / self.NZ * self.NX / 2.0, self.D_CA],
            number_of_cells=[self.NX, self.NZ],
            min_tiles=16
        )
        mwxrun.grid.potential_zmax = self.ANODE_VOLTAGE
        mwxrun.init_timestep(DT=self.DT)
        mwxrun.simulation.max_steps = self.MAX_STEPS
        mwxrun.simulation.load_balance_intervals = self.MAX_STEPS // 5000

        #######################################################################
        # Field solver                                                        #
        #######################################################################

        # self.solver = picmi.ElectrostaticSolver(
        #    grid=mwxrun.grid, method='Multigrid', required_precision=1e-12
        # )
        self.solver = poisson_pseudo_1d.PoissonSolverPseudo1D(grid=mwxrun.grid)
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
            species2=self.ions, npart=2 * self.SEED_NPPC * self.NX * self.NZ,
            T_2=self.T_INERT, plasma_density=self.PLASMA_DENSITY
        )

        #######################################################################
        # Initialize run and print diagnostic info                            #
        #######################################################################

        mwxrun.init_run()

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        self.text_diag = diag_base.TextDiag(
            diag_steps=self.DIAG_STEPS, preset_string='perfdebug'
        )

        self.field_diag = FieldDiagnostic(
            diag_steps=self.DIAG_STEPS, barrier_slices=[0],
            save_pdf=False, style='roelof', min_dim=2.0
        )

        # array to save ion density
        self.rho_array = np.zeros(self.NZ + 1)

    def _get_rho_ions(self):
        rho_data = mwxrun.get_gathered_rho_grid('he_ions', False)
        if mwxrun.me == 0:
            self.rho_array += (
                np.mean(rho_data[:,:,0], axis=0) / constants.q_e
                / self.DIAG_STEPS
            )

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
