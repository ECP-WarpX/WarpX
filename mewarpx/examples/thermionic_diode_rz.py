import argparse
import sys

import numpy as np

from mewarpx.utils_store import util as mwxutil

mwxutil.init_libwarpx(ndim=2, rz=True)

from mewarpx import emission
from mewarpx.mwxrun import mwxrun
from mewarpx.setups_store import diode_setup


def run_simulation(V_bias, steps, save_current):

    radius_frac = 0.8

    ####################################
    # Diode setup
    ####################################

    run = diode_setup.DiodeRun_V1(
        CATHODE_TEMP = 1100.0 + 273.15, # K
        CATHODE_A = 6e5, # A/m^2/K^2
        CATHODE_PHI = 2.11, # eV
        USE_SCHOTTKY = False,
        ANODE_TEMP = 200, # K
        ANODE_PHI = 1.4, # eV
        V_ANODE_CATHODE = V_bias, # V
        D_CA = 50e-6, # m
        DT = 0.5e-12, # s
        NX = 32,
        NZ = 128,
        DIRECT_SOLVER = False,
        NPPC = 1,
        TOTAL_TIMESTEPS = steps,
        DIAG_STEPS = ((steps // 5) if steps > 10 else steps),
    )

    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_injectors=False
    )

    #################################
    # Setup cathode injector
    #################################

    run.emitter = emission.ZDiscEmitter(
        conductor=run.cathode, T=run.CATHODE_TEMP,
        outer_emission_radius=run.PERIOD*radius_frac,
        transverse_fac=run.TRANSVERSE_FAC
    )
    run.injector = emission.ThermionicInjector(
        run.emitter, run.electrons, run.NPPC, run.CATHODE_TEMP,
        run.CATHODE_PHI, run.CATHODE_A, run.USE_SCHOTTKY
    )

    #################################
    # Initialize WarpX and fluxdiags
    #################################

    run.init_warpx()
    run.init_runinfo()
    run.init_fluxdiag()

    #################################
    # Simulation run
    #################################

    mwxrun.simulation.step(steps)

    #################################
    # Save IV results
    #################################

    if save_current and mwxrun.me == 0:
        key = ('scrape', 'anode', 'electrons')
        J_diode = run.fluxdiag.ts_dict[key].get_averagevalue_by_key('J')
        # normalize appropriate
        J_diode *= 1.0 / radius_frac**2
        print(f'{V_bias}    {J_diode}')
        with open(f'results_d_{int(run.D_CA*1e6)}.dat', 'a') as f:
            f.write(f'{V_bias}    {J_diode}\n')

parser = argparse.ArgumentParser()
parser.add_argument('--V', help='bias voltage in Volt', type=float, default=0)
parser.add_argument('--steps', help='set the number of simulation steps',
                    type=int, default=3000)
parser.add_argument('--save', help='save voltage and current pairs',
                    default=False, action='store_true')
args, left = parser.parse_known_args()
sys.argv = sys.argv[:1]+left

run_simulation(args.V, args.steps, args.save)
