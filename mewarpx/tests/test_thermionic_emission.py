""" Test pseudo 1D diode run with thermionic emission"""
import os
import numpy as np

from mewarpx import util as mwxutil

def test_thermionic_emission():
    name = "Thermionic Emission"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    import mewarpx.mwxconstants as constants

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239475)

    TOTAL_TIME = 1e-10 # s
    DIAG_INTERVAL = 5e-10 # s
    DT = 0.5e-12 # s

    P_INERT = 1 # torr
    T_INERT = 300 # K
    N_INERT = (P_INERT * constants.torr_SI) / (constants.kb_J * T_INERT)

    USE_SCHOTTKY = False

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV
    NX = 8
    NZ = 128

    DIRECT_SOLVER = True

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    run = diode_setup.DiodeRun_V1(
        dim=dim,
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        USE_SCHOTTKY=USE_SCHOTTKY,
        N_INERT=N_INERT,
        T_INERT=T_INERT,
        T_ELEC=T_INERT,
        PLASMA_DENSITY=0,
        NPPC=50,
        NX=NX,
        NZ=NZ,
        DIRECT_SOLVER=DIRECT_SOLVER,
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=max_steps,
        DIAG_STEPS=diag_steps,
        DIAG_INTERVAL=DIAG_INTERVAL,
        NUMBER_PARTICLES_PER_CELL=[16, 16]
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    net_rho_grid = np.array(mwxrun.get_gathered_rho_grid()[0][:, :, 0])
    ref_path = os.path.join(testing_util.test_dir,
                            "thermionic_emission",
                            "thermionic_emission.npy")
    ref_rho_grid = np.load(ref_path)

    assert np.allclose(net_rho_grid, ref_rho_grid)
