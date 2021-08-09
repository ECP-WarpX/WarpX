"""Test pseudo 1D diode run with the superLU solver."""
import os
import numpy as np

from mewarpx.utils_store import util as mwxutil


def test_superLU_solver():
    name = "superLU_solver"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # set to False to regenerate the reference data
    DIRECT_SOLVER = True

    # Specific numbers match older run for consistency
    FREQ = 13.56e6  # MHz
    DT = 1.0 / (400 * FREQ)
    DIAG_STEPS = 50
    DIAG_INTERVAL = DIAG_STEPS*DT
    VOLTAGE = 450.0
    D_CA = 0.067  # m
    NX = 16
    NZ = 128
    run = diode_setup.DiodeRun_V1(
        dim=dim,
        DIRECT_SOLVER=DIRECT_SOLVER,
        V_ANODE_EXPRESSION=f"{VOLTAGE}*sin(2*pi*{FREQ:.5e}*t)",
        D_CA=D_CA,
        INERT_GAS_TYPE='He',
        N_INERT=9.64e20,  # m^-3
        T_INERT=300.0,  # K
        PLASMA_DENSITY=2.56e14,  # m^-3
        T_ELEC=30000.0,  # K
        SEED_NPPC=16*32,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=50,
        DIAG_STEPS=DIAG_STEPS,
        DIAG_INTERVAL=DIAG_INTERVAL
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_injectors=False,
        init_inert_gas=True,
        init_mcc=True,
        init_neutral_plasma=True,
        init_field_diag=False,
        init_simcontrol=True,
        init_warpx=True
    )

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check phi results against reference data from MLMG solver           #
    #######################################################################

    data = np.mean(mwxrun.get_gathered_phi_grid()[0], axis=0)
    # uncomment to generate reference data
    # np.save('reference_data.npy', data)

    # the MLMG solver has the last data point set to 0
    data[-1] = 0

    ref_data = np.load(os.path.join(
        testing_util.test_dir, 'direct_solver', 'reference_data.npy'
    ))

    assert np.allclose(data, ref_data, rtol=0.002)
