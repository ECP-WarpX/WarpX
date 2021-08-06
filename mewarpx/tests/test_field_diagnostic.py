"""
Test that verifies data is plotted correctly from the
Monte-Carlo Collision script based on case 1 from
Turner et al. (2013) - https://doi.org/10.1063/1.4775084
"""

import pytest
import os.path

from mewarpx import util as mwxutil

from pywarpx import picmi

import numpy as np

constants = picmi.constants

# yt seems to be causing this error when a dataset is loaded
@pytest.mark.filterwarnings("ignore::ResourceWarning")
@pytest.mark.parametrize("plot_on_diag_steps",
[
    True,
    False
])
def test_field_diag(plot_on_diag_steps):
    dim = 2
    use_rz = False
    # We test either post processing or plotting on diag steps, not both.
    post_processing = not plot_on_diag_steps

    mwxutil.init_libwarpx(ndim=dim, rz=use_rz)

    from mewarpx.mwxrun import mwxrun
    from mewarpx.setups_store import diode_setup
    from mewarpx import testing_util

    plot_diag_str = "_with_diag_plotting" if plot_on_diag_steps else ""
    post_proc_str = "_with_post_processing" if post_processing else ""
    test_name = "field_diag_test" + plot_diag_str + post_proc_str

    testing_util.initialize_testingdir(test_name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(83197410)

    STEPS = 10
    D_CA = 0.067 # m
    FREQ = 13.56e6 # MHz
    VOLTAGE = 450.0
    NX = 8
    NZ = 128
    NUMBER_PER_CELL_EACH_DIM = [16, 32]
    DT = 1.0 / (400 * FREQ)
    DIAG_STEPS = 2
    DIAG_DATA_LIST = ['rho_electrons', 'rho_he_ions', 'phi']

    run = diode_setup.DiodeRun_V1(
        dim=dim,
        rz=use_rz,
        V_ANODE_CATHODE=VOLTAGE,
        V_ANODE_EXPRESSION="%.1f*sin(2*pi*%.5e*t)" % (VOLTAGE, FREQ),
        D_CA=D_CA,
        INERT_GAS_TYPE='He',
        N_INERT=9.64e20,  # m^-3
        T_INERT=300.0,  # K
        PLASMA_DENSITY=2.56e14,  # m^-3
        T_ELEC=30000.0,  # K
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=STEPS,
        DIAG_STEPS=DIAG_STEPS,
        NUMBER_PARTICLES_PER_CELL=NUMBER_PER_CELL_EACH_DIM,
        FIELD_DIAG_DATA_LIST=DIAG_DATA_LIST,
        FIELD_DIAG_PLOT_ON_DIAG_STEPS=plot_on_diag_steps,
        FIELD_DIAG_PLOT_DATA_LIST=['rho', 'phi'],
        FIELD_DIAG_PLOT_AFTER_RUN=post_processing
    )

    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_injectors=False,
        init_inert_gas=True,
        init_field_diag=True,
        init_simcontrol=True,
        init_warpx=True,
        init_simulation=True
    )
    while run.control.check_criteria():
        mwxrun.simulation.step()

    # verify that the plot images were created.
    if plot_on_diag_steps:
        print("Verifying that all plots were created...", flush=True)
        for i in range(1,6):
            assert os.path.isfile("phi_after_diag_step_" + str(i) + ".png")
            assert os.path.isfile("rho_after_diag_step_" + str(i) + ".png")
        print("All plots exist!", flush=True)

    # verify that the post processing image was created
    if post_processing:
        print("Verifying that all plots were created...", flush=True)
        for i in range(0, STEPS + 1, DIAG_STEPS):
            for param in DIAG_DATA_LIST:
                assert os.path.isfile("diags/fields/" + param + "_" + f"{i:05d}.png"), param + "_" + f"{i:06d}.png doesn't exist"
        print("All plots exist!", flush=True)


if __name__ == '__main__':
    test_field_diag(False)
