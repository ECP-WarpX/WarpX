"""Test pseudo 1D diode run with the superLU solver."""
import pytest
import numpy as np
import pandas

from mewarpx import util as mwxutil

def test_direct_solver():
    name = "Direct_solver"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import testing_util
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
        V_ANODE_CATHODE=VOLTAGE,
        V_ANODE_EXPRESSION=(
            lambda t: VOLTAGE * np.sin(2.0 * np.pi * FREQ * t) if DIRECT_SOLVER
            else "%.1f*sin(2*pi*%.5e*t)" % (VOLTAGE, FREQ)
        ),
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
        TOTAL_TIMESTEPS=50,
        DIAG_STEPS=DIAG_STEPS,
        DIAG_INTERVAL=DIAG_INTERVAL,
        NUMBER_PARTICLES_PER_CELL=[16, 32],
        FIELD_DIAG_DATA_LIST=['rho_electrons', 'rho_he_ions', 'phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_injectors=False,
        init_inert_gas_ionization=True,
        init_field_diag=True,
        init_simcontrol=True,
        init_warpx=True
    )

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check rho and phi results against reference data from MGML solver   #
    #######################################################################

    df = pandas.DataFrame(columns=['rho', 'phi'])
    df['rho'] = np.mean(mwxrun.get_gathered_rho_grid()[0], axis=0)[1:-1,0]
    df['phi'] = np.mean(mwxrun.get_gathered_phi_grid()[0], axis=0)
    assert testing_util.test_df_vs_ref(name, df)
