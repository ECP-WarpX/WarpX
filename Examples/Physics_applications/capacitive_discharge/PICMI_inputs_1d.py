#!/usr/bin/env python3
#
# --- Copyright 2021 Modern Electron
# --- Monte-Carlo Collision script to reproduce the benchmark tests from
# --- Turner et al. (2013) - https://doi.org/10.1063/1.4775084

import argparse
import sys

import numpy as np
from pywarpx import callbacks, fields, picmi

constants = picmi.constants


class CapacitiveDischargeExample(object):
    '''The following runs a simulation of a parallel plate capacitor seeded
    with a plasma in the spacing between the plates. A time varying voltage is
    applied across the capacitor. The groups of 4 values below correspond to
    the 4 cases simulated by Turner et al. (2013) in their benchmarks of
    PIC-MCC codes.
    '''

    gap = 0.067 # m

    freq = 13.56e6 # Hz
    voltage = [450.0, 200.0, 150.0, 120.0] # V

    gas_density = [9.64e20, 32.1e20, 96.4e20, 321e20] # m^-3
    gas_temp = 300.0 # K
    m_ion = 6.67e-27 # kg

    plasma_density = [2.56e14, 5.12e14, 5.12e14, 3.84e14] # m^-3
    elec_temp = 30000.0 # K

    seed_nppc = 16 * np.array([32, 16, 8, 4])

    nz = [128, 256, 512, 512]

    dt = 1.0 / (np.array([400, 800, 1600, 3200]) * freq)

    # Total simulation time in seconds
    total_time = np.array([1280, 5120, 5120, 15360]) / freq
    # Time (in seconds) between diagnostic evaluations
    diag_interval = 32 / freq

    def __init__(self, n=0, test=False):
        """Get input parameters for the specific case (n) desired."""
        self.test = test

        # Case specific input parameters
        self.voltage = f"{self.voltage[n]}*sin(2*pi*{self.freq:.5e}*t)"

        self.gas_density = self.gas_density[n]
        self.plasma_density = self.plasma_density[n]
        self.seed_nppc = self.seed_nppc[n]

        self.nz = self.nz[n]

        self.dt = self.dt[n]
        self.max_steps = int(self.total_time[n] / self.dt)
        self.diag_steps = int(self.diag_interval / self.dt)

        if self.test:
            self.max_steps = 20
            self.diag_steps = 5

        self.rho_wrapper = None
        self.ion_density_array = np.zeros(self.nz + 1)

        self.setup_run()

    def setup_run(self):
        """Setup simulation components."""

        #######################################################################
        # Set geometry and boundary conditions                                #
        #######################################################################

        self.grid = picmi.Cartesian1DGrid(
            number_of_cells=[self.nz],
            warpx_max_grid_size=128,
            lower_bound=[0],
            upper_bound=[self.gap],
            lower_boundary_conditions=['dirichlet'],
            upper_boundary_conditions=['dirichlet'],
            lower_boundary_conditions_particles=['absorbing'],
            upper_boundary_conditions_particles=['absorbing'],
            warpx_potential_hi_z=self.voltage,
        )

        #######################################################################
        # Field solver                                                        #
        #######################################################################

        self.solver = picmi.ElectrostaticSolver(
            grid=self.grid, method='Multigrid', required_precision=1e-12,
            warpx_self_fields_verbosity=0
        )

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        self.electrons = picmi.Species(
            particle_type='electron', name='electrons',
            initial_distribution=picmi.UniformDistribution(
                density=self.plasma_density,
                rms_velocity=[np.sqrt(constants.kb * self.elec_temp / constants.m_e)]*3,
            )
        )
        self.ions = picmi.Species(
            particle_type='He', name='he_ions',
            charge='q_e', mass=self.m_ion,
            initial_distribution=picmi.UniformDistribution(
                density=self.plasma_density,
                rms_velocity=[np.sqrt(constants.kb * self.gas_temp / self.m_ion)]*3,
            )
        )

        #######################################################################
        # Collision  initialization                                           #
        #######################################################################

        cross_sec_direc = '../../../../warpx-data/MCC_cross_sections/He/'
        mcc_electrons = picmi.MCCCollisions(
            name='coll_elec',
            species=self.electrons,
            background_density=self.gas_density,
            background_temperature=self.gas_temp,
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
            background_density=self.gas_density,
            background_temperature=self.gas_temp,
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
            time_step_size=self.dt,
            max_steps=self.max_steps,
            warpx_collisions=[mcc_electrons, mcc_ions],
            warpx_load_balance_intervals=self.max_steps//5000,
            verbose=self.test
        )

        self.sim.add_species(
            self.electrons,
            layout = picmi.GriddedLayout(
                n_macroparticle_per_cell=[self.seed_nppc], grid=self.grid
            )
        )
        self.sim.add_species(
            self.ions,
            layout = picmi.GriddedLayout(
                n_macroparticle_per_cell=[self.seed_nppc], grid=self.grid
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

    def _get_rho_ions(self):
        # deposit the ion density in rho_fp
        self.sim.extension.depositChargeDensity('he_ions', 0)

        rho_data = fields.RhoFPWrapper(0, False)[...][:,0]
        self.ion_density_array += rho_data / constants.q_e / self.diag_steps

    def run_sim(self):

        self.sim.step(self.max_steps - self.diag_steps)
        callbacks.installafterstep(self._get_rho_ions)
        self.sim.step(self.diag_steps)

        if self.sim.extension.getMyProc() == 0:
            np.save('avg_ion_density.npy', self.ion_density_array)

##########################
# parse input parameters
##########################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-t', '--test', help='toggle whether this script is run as a short CI test',
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

run = CapacitiveDischargeExample(n=args.n-1, test=args.test)
run.run_sim()
