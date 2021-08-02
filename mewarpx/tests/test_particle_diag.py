"""
Test that verifies Particle Diagnostics is correctly outputting the correct
plots from the post processed yt data from the picmi.ParticleDiagnostics
"""

import os
import numpy as np
import pytest

from mewarpx import util as mwxutil

@pytest.mark.filterwarnings("ignore::ResourceWarning")
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_particle_diag():
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    import mewarpx.mwxconstants as constants

    test_name = "particle_diag_test_with_post_processing"

    testing_util.initialize_testingdir(test_name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239475)

    TOTAL_TIME = 1e-10 # s
    DIAG_INTERVAL = 5e-10 # s
    DT = 0.5e-12 # s

    P_INERT = 1 # torr
    T_INERT = 300 # K

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV
    NX = 8
    NZ = 128

    max_steps = 10
    diag_steps = 2

    DATA_LIST = ['position', 'momentum', 'weighting']
    PLOT_SPECIES = ['electrons']
    DIAG_PLOT_DATA_LIST = ["particle_position_x", "particle_momentum_x"]

    run = diode_setup.DiodeRun_V1(
        dim=dim,
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        P_INERT=P_INERT,
        T_INERT=T_INERT,
        NPPC=50,
        NX=NX,
        NZ=NZ,
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=max_steps,
        DIAG_STEPS=diag_steps,
        PARTICLE_DIAG_DATA_LIST=DATA_LIST,
        PARTICLE_PLOT_SPECIES=PLOT_SPECIES,
        PARTICLE_DIAG_PLOT_DATA_LIST=DIAG_PLOT_DATA_LIST,
        PARTICLE_DIAG_PLOT_AFTER_RUN=True
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_particle_diag=True,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    # verify that the plot images were created.

    print('Verifying that all plots were created...')
    for i in range(diag_steps, max_steps + 1, diag_steps):
        for specimen in PLOT_SPECIES:
            for param in DIAG_PLOT_DATA_LIST:
                file_name = os.path.join(
                    os.getcwd(),
                    'diags',
                    specimen+'_'+param+'_diag'+f'{i:05d}.jpg'
                )
                print(file_name)
                assert os.path.isfile(file_name)
        print('All plots exist!')

if __name__ == "__main__":
    test_particle_diag()
