"""Tests for functionality in mwxrun.py"""
import pytest
import numpy as np
import yt

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
@pytest.mark.filterwarnings("ignore::ResourceWarning")
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
        init_neutral_plasma=True,
        init_mcc=True,
        init_field_diag=True,
        init_simcontrol=True,
        init_warpx=True
    )

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Compare mwxrun data gathering to diagnostic output                  #
    #######################################################################

    # the edges are treated differently in our gather rho versus the yt
    # dump so we only look at the interior here
    rho_data = np.mean(mwxrun.get_gathered_rho_grid(
        species_name='he_ions', include_ghosts=False
    )[0][:,1:,0], axis=0)[2:-2]

    ds = yt.load( "diags/fields/fields00010" )
    grid_data = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    rho_data_yt = np.mean(
        grid_data['rho_he_ions'].to_ndarray()[:,:,0], axis=0
    )[2:-2]

    assert np.allclose(rho_data, rho_data_yt, rtol=0.01)

    #######################################################################
    # Cleanup and final output                                            #
    #######################################################################

    out, _ = capsys.readouterr()

    print(out)
    # make sure out isn't empty
    outstr = "SimControl: Termination from criteria: eval_max_steps"
    assert outstr in out
