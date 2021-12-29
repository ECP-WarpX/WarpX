import argparse
import sys

from mewarpx.utils_store import util as mwxutil

mwxutil.init_libwarpx(ndim=2, rz=True)

import numpy as np
from pywarpx import picmi

from mewarpx import assemblies, diags, emission, mespecies, runinfo
from mewarpx.mwxrun import mwxrun


class CylinderVacuumTEC(object):

    #######################################################################
    # Begin global user parameters                                        #
    #######################################################################

    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_A = 6e5              # A/m^2/K^2
    CATHODE_PHI = 2.11           # eV
    USE_SCHOTTKY = False
    V_CATHODE = -CATHODE_PHI     # cathode is grounded

    ANODE_TEMP = 200             # K
    ANODE_PHI = 1.4              # eV

    D_CA = 50e-6                 # m

    NZ = 32
    NR = 512

    DT = 0.5e-12                 # s

    NPPC = 4

    #######################################################################
    # End global user parameters and user input                           #
    #######################################################################

    def __init__(self, V_ANODE_CATHODE, TOTAL_TIMESTEPS, SAVE, USE_EB=False,
                 DIAG_STEPS=None):
        self.V_ANODE_CATHODE = V_ANODE_CATHODE
        self.TOTAL_TIMESTEPS = TOTAL_TIMESTEPS
        self.SAVE = SAVE
        self.USE_EB = USE_EB
        self.DIAG_STEPS = DIAG_STEPS

    def setup_run(self):

        self.rmax = 1e-3
        self.rmin = self.rmax - self.D_CA - self.USE_EB * 10e-6
        self.zmax = (self.rmax - self.rmin) / self.NR * self.NZ / 2.0
        self.zmin = -self.zmax

        self.V_ANODE = self.V_CATHODE + self.V_ANODE_CATHODE

        if self.DIAG_STEPS is None:
            self.DIAG_STEPS = (
                (self.TOTAL_TIMESTEPS // 5) if self.TOTAL_TIMESTEPS > 10
                else self.TOTAL_TIMESTEPS
            )

        #######################################################################
        # Set geometry, boundary conditions and timestep                      #
        #######################################################################

        mwxrun.init_grid(
            lower_bound=[self.rmin, self.zmin],
            upper_bound=[self.rmax, self.zmax],
            number_of_cells=[self.NR, self.NZ], # min_tiles=1,
            bc_fields_z_min='periodic', bc_fields_z_max='periodic',
            bc_particles_z_min='periodic', bc_particles_z_max='periodic',
            bc_fields_r_min='dirichlet', bc_fields_r_max='dirichlet',
        )
        mwxrun.init_timestep(DT=self.DT)
        mwxrun.simulation.max_steps = self.TOTAL_TIMESTEPS
        mwxrun.simulation.load_balance_intervals = self.TOTAL_TIMESTEPS // 5000

        #######################################################################
        # Field solver                                                        #
        #######################################################################

        self.solver = picmi.ElectrostaticSolver(
            grid=mwxrun.grid, method='Multigrid', required_precision=1e-6,
            warpx_self_fields_verbosity=0, maximum_iterations=1000
        )
        mwxrun.simulation.solver = self.solver

        #######################################################################
        # Conductors setup and installation                                   #
        #######################################################################

        self.cathode = assemblies.CylinderZ(
            V=self.V_CATHODE, T=self.CATHODE_TEMP, WF=self.CATHODE_PHI,
            name='cathode', r_outer=self.rmax+1e-6, r_inner=self.rmax
        )
        self.anode = assemblies.CylinderZ(
            V=self.V_ANODE, T=self.ANODE_TEMP, WF=self.ANODE_PHI,
            name='anode',
            r_outer=(self.rmin if not self.USE_EB else self.rmax - self.D_CA)
        )

        #######################################################################
        # Particle types setup                                                #
        #######################################################################

        self.electrons = mespecies.Species(
            particle_type='electron', name='electrons'
        )

        #######################################################################
        # Scraper setup                                                       #
        #######################################################################

        for species in [self.electrons]:
            species.save_particles_at_xhi = 1
            if self.USE_EB:
                species.save_particles_at_eb = 1
            else:
                species.save_particles_at_xlo = 1

        #######################################################################
        # Thermionic emission                                                 #
        #######################################################################

        self.cathode_emitter = emission.ZCylinderEmitter(
            conductor=self.cathode, rdir=-1)
        self.cathode_injector = emission.ThermionicInjector(
            emitter=self.cathode_emitter, species=self.electrons,
            A=self.CATHODE_A, npart_per_cellstep=self.NPPC,
            use_Schottky=self.USE_SCHOTTKY
        )

        #######################################################################
        # Clean up and output RunInfo                                         #
        #######################################################################

        runvars = CylinderVacuumTEC.__dict__.copy()
        runvars.update(self.__dict__)

        self.runinfo = runinfo.RunInfo(
            local_vars=runvars,
            electrode_params={"CATHODE_A": self.CATHODE_A,
                              "CATHODE_TEMP": self.CATHODE_TEMP,
                              "CATHODE_PHI": self.CATHODE_PHI,
                              "ANODE_PHI": self.ANODE_PHI},
            surface_dict={'cathode': self.cathode, 'anode': self.anode},
            injector_dict={'cathode': self.cathode_injector},
        )
        # overwrite runinfo area to properly normalize currents
        self.runinfo.area = self.cathode_emitter.area

        #######################################################################
        # Add diagnostics                                                     #
        #######################################################################

        self.text_diag = diags.TextDiag(
            diag_steps=self.DIAG_STEPS, preset_string='perfdebug'
        )
        # self.field_diag = diags.FieldDiagnostic(
        #    diag_steps=self.DIAG_STEPS, style='roelof', save_pdf=False
        # )
        self.fluxdiag = diags.FluxDiagnostic(
            diag_steps=self.DIAG_STEPS,
            runinfo=self.runinfo,
            check_charge_conservation=False,
            overwrite=True
        )

        #######################################################################
        # Initialize run and print diagnostic info                            #
        #######################################################################

        mwxrun.init_run()


    def run_sim(self):

        mwxrun.simulation.step(self.TOTAL_TIMESTEPS)

        if self.SAVE and mwxrun.me == 0:
            key = ('scrape', 'anode', 'electrons')
            J_diode = self.fluxdiag.ts_dict[key].get_averagevalue_by_key('J')
            print(f'{self.V_ANODE_CATHODE}    {J_diode}')
            with open(f'results_d_{int(self.D_CA*1e6)}.dat', 'a') as f:
                f.write(f'{self.V_ANODE_CATHODE}    {J_diode}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--V', help='bias voltage in Volt', type=float, default=0)
    parser.add_argument('--steps', help='set the number of simulation steps',
                        type=int, default=3000)
    parser.add_argument('--save', help='save voltage and current pairs',
                        default=False, action='store_true')
    args, left = parser.parse_known_args()
    sys.argv = sys.argv[:1]+left

    run = CylinderVacuumTEC(args.V, args.steps, args.save)

    run.setup_run()
    run.run_sim()
