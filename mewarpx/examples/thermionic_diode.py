import argparse
import sys

import numpy as np

from mewarpx.utils_store import util as mwxutil

mwxutil.init_libwarpx(ndim=2, rz=False)

from mewarpx.mwxrun import mwxrun
from mewarpx.setups_store import diode_setup


def run_simulation(V_bias, steps, save_current):

    ####################################
    # Diode setup
    ####################################

    run = diode_setup.DiodeRun_V1(
        CATHODE_TEMP = 1100.0 + 273.15, # K
        CATHODE_A = 6e5, # A/m^2/K^2
        CATHODE_PHI = 2.11, # eV
        USE_SCHOTTKY = True,
        ANODE_TEMP = 200, # K
        ANODE_PHI = 1.4, # eV
        V_ANODE_CATHODE = V_bias, # V
        D_CA = 50e-6, # m
        DT = 0.5e-12, # s
        NX = 8,
        NZ = 128,
        DIRECT_SOLVER = True,
        NPPC = 100,
        TOTAL_TIMESTEPS = steps,
        DIAG_STEPS = ((steps // 5) if steps > 10 else steps),
    )

    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_runinfo=True,
        init_fluxdiag=True,
        init_simcontrol=True,
        init_warpx=True
    )

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
