import argparse
import sys

import numpy as np

from mewarpx import assemblies, emission, runinfo
from mewarpx.mwxrun import mwxrun
from mewarpx.setups_store import diode_setup


class PatchyVacuumTEC(object):

    #######################################################################
    # Begin global user parameters                                        #
    #######################################################################

    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_A = 6e5              # A/m^2/K^2
    CATHODE_PHI = 2.11           # eV
    V_CATHODE = -CATHODE_PHI     # cathode is grounded

    ANODE_TEMP = 200             # K
    ANODE_PHI = 1.4              # eV

    D_CA = 50e-6                 # m

    NZ = 128
    NX = 64

    DT = 0.5e-12                 # s

    NPPC = 100

    FIELD_DIAG_PLOT = False

    #######################################################################
    # End global user parameters and user input                           #
    #######################################################################

    def __init__(self, V_ANODE_CATHODE, TOTAL_TIMESTEPS, SAVE, DIAG_STEPS=None,
                 DELTA_PHI=0.2, PATCH_SIZE=5e-6, USE_SCHOTTKY=False):

        self.V_ANODE_CATHODE = V_ANODE_CATHODE
        self.TOTAL_TIMESTEPS = TOTAL_TIMESTEPS
        self.SAVE = SAVE
        self.DELTA_PHI = DELTA_PHI
        self.PATCH_SIZE = PATCH_SIZE
        self.USE_SCHOTTKY = USE_SCHOTTKY

        self.DIAG_STEPS = DIAG_STEPS
        if self.DIAG_STEPS is None:
            self.DIAG_STEPS = self.TOTAL_TIMESTEPS // 5

    def setup_run(self):

        ####################################
        # Diode setup
        ####################################

        self.run = diode_setup.DiodeRun_V1(
            GEOM_STR='XZ',
            CATHODE_TEMP=self.CATHODE_TEMP,
            CATHODE_A=self.CATHODE_A,
            CATHODE_PHI=self.CATHODE_PHI,
            USE_SCHOTTKY=self.USE_SCHOTTKY,
            ANODE_TEMP=self.ANODE_TEMP,
            ANODE_PHI=self.ANODE_PHI,
            V_ANODE_CATHODE=self.V_ANODE_CATHODE,
            D_CA=self.D_CA,
            DT=self.DT,
            NX=self.NX,
            NZ=self.NZ,
            DIRECT_SOLVER=False,
            NPPC=self.NPPC,
            TOTAL_TIMESTEPS=self.TOTAL_TIMESTEPS,
            DIAG_STEPS=self.DIAG_STEPS,
            FIELD_DIAG_PLOT=self.FIELD_DIAG_PLOT
        )

        self.run.init_base()
        self.run.init_electrons()

        ####################################
        # Conductor setup
        ####################################

        self.patchy_cathode = assemblies.PatchyCathode(
            V=self.run.V_CATHODE, T=self.run.CATHODE_TEMP,
            phi_bar=self.run.CATHODE_PHI, delta_phi=self.DELTA_PHI,
            patch_size=self.PATCH_SIZE
        )
        self.run.cathode = self.patchy_cathode.cathode_list
        self.run.electrons.save_particles_at_eb = 1

        self.run.anode = assemblies.Anode(
            z=self.run.D_CA, V=self.run.V_ANODE, T=self.run.ANODE_TEMP,
            WF=self.run.ANODE_PHI
        )
        self.run.electrons.save_particles_at_zhi = 1

        self.run.surface_list = self.run.cathode + [self.run.anode]

        self.run.init_solver()

        ####################################
        # Injector setup
        ####################################

        self.run.emitter = emission.ZPlanePatchyEmitter(
            patchy_cathode=self.patchy_cathode,
            ymin=0.0, ymax=0.0,
            transverse_fac=self.run.TRANSVERSE_FAC
        )

        self.run.injector1 = emission.ThermionicInjector(
            self.run.emitter.high_wf, self.run.electrons, self.run.NPPC,
            A=self.run.CATHODE_A, use_Schottky=self.run.USE_SCHOTTKY
        )
        self.run.injector2 = emission.ThermionicInjector(
            self.run.emitter.low_wf, self.run.electrons, self.run.NPPC,
            A=self.run.CATHODE_A, use_Schottky=self.run.USE_SCHOTTKY
        )

        ####################################
        # Runinfo
        ####################################

        injector_dict = {
            'cathode': [self.run.injector1, self.run.injector2]
        }
        surface_dict = {
            'cathode': self.run.cathode[0],
            'cathode2': self.run.cathode[1],
            'anode': self.run.anode
        }

        self.run.runinfo = runinfo.RunInfo(
            injector_dict=injector_dict,
            surface_dict=surface_dict,
            local_vars=self.run.__dict__,
            run_file=__file__,
            # TODO: implement collecting run parameters
            # run_param_dict=self.setupinfo.run_param_dict,
            electrode_params={
                "CATHODE_A": self.run.CATHODE_A,
                "CATHODE_TEMP": self.run.CATHODE_TEMP,
                "CATHODE_PHI": self.run.CATHODE_PHI,
                "ANODE_PHI": self.run.ANODE_PHI
            },
        )

        self.run.init_fluxdiag()
        self.run.init_field_diag()
        self.run.init_simcontrol()
        self.run.init_simulation()
        self.run.init_warpx()

    def run_sim(self):

        #################################
        # Simulation run
        #################################

        mwxrun.simulation.step(self.TOTAL_TIMESTEPS)

        #################################
        # Save IV results
        #################################

        if self.SAVE and mwxrun.me == 0:
            key = ('scrape', 'anode', 'electrons')
            J_diode = self.run.fluxdiag.ts_dict[key].get_averagevalue_by_key('J')
            print(f'{self.V_ANODE_CATHODE}    {J_diode}')
            with open(f'results_d_{int(self.run.D_CA*1e6)}.dat', 'a') as f:
                f.write(f'{self.V_ANODE_CATHODE}    {J_diode}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--V', help='bias voltage in Volt', type=float,
                        default=0)
    parser.add_argument('--steps', help='set the number of simulation steps',
                        type=int, default=3000)
    parser.add_argument('--save', help='save voltage and current pairs',
                        default=False, action='store_true')
    parser.add_argument('--schottky', help='use Schottky enhancement',
                        default=False, action='store_true')
    args, left = parser.parse_known_args()
    sys.argv = sys.argv[:1]+left

    run = PatchyVacuumTEC(args.V, args.steps, args.save,
                          USE_SCHOTTKY=args.schottky)
    run.setup_run()
    run.run_sim()
