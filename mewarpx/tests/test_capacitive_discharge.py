"""Test full 1D diode run with diagnostics."""
import pytest
import numpy as np

from mewarpx import util as mwxutil


@pytest.mark.parametrize(
    ("name"),
    [
        'Run2D',
        # For these two to work we'll need to allow python to choose which of
        # multiple libwarpx...so files it should load.
        # 'Run2D_RZ',
        # 'Run3D'
    ]
)
def test_capacitive_discharge_multigrid(capsys, name):
    basename = "Run"
    use_rz = 'RZ' in name
    dim = int(name.replace(basename, '')[0])

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=use_rz)
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

    # Specific numbers match older run for consistency
    FREQ = 13.56e6  # MHz
    DT = 1.0 / (400 * FREQ)
    DIAG_STEPS = 2
    DIAG_INTERVAL = DIAG_STEPS*DT
    VOLTAGE = 450.0
    D_CA = 0.067  # m
    NX = 16
    NZ = 128
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
        SEED_NPPC=16*32,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=10,
        DIAG_STEPS=DIAG_STEPS,
        DIAG_INTERVAL=DIAG_INTERVAL,
        FIELD_DIAG_DATA_LIST=['rho_electrons', 'rho_he_ions', 'phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_injectors=False,
        init_inert_gas=True,
        init_mcc=True,
        init_field_diag=True,
        init_simcontrol=True,
        init_warpx=True
    )

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Cleanup and final output                                            #
    #######################################################################

    if capsys is not None:
        out, _ = capsys.readouterr()

        print(out)
        # make sure out isn't empty
        outstr = "SimControl: Termination from criteria: eval_max_steps"
        assert outstr in out
    else:
        assert False, "This wasn't run in pytest, but it passed otherwise."


if __name__ == '__main__':
    test_capacitive_discharge_multigrid(None, "Run2D")
