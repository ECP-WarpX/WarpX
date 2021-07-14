"""Test the custom criteria using capacitive discharge"""
import pytest
import numpy as np

from mewarpx import util as mwxutil

def check_particle_nums():
    from mewarpx import mwxrun
    return mwxrun.mwxrun.get_npart() < 40000

def add_particles():
    from pywarpx import _libwarpx
    # add 10,000 particles
    _libwarpx.add_particles(x=np.zeros((10000)), y=np.zeros((10000)), z=np.zeros((10000)))

@pytest.mark.parametrize(
    ("name"),
    [
        'Run2D',
        # 'Run2D_RZ',
        # 'Run3D'
    ]
)
def test_capacitive_discharge_custom_criteria(capsys, name):
    basename = "Run"
    use_rz = 'RZ' in name
    dim = int(name.replace(basename, '')[0])

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=use_rz)
    from mewarpx import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun
    from pywarpx import callbacks

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
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * 16 / 128.0,
        DT=DT,
        TOTAL_TIMESTEPS=10,
        DIAG_STEPS=DIAG_STEPS,
        DIAG_INTERVAL=DIAG_INTERVAL,
        NUMBER_PARTICLES_PER_CELL=[2, 4],
        FIELD_DIAG_DATA_LIST=['rho_electrons', 'rho_he_ions', 'phi']
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
    run.control.add_checker(check_particle_nums)

    callbacks.installafterstep(add_particles)

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step(1)

    #######################################################################
    # Cleanup and final output                                            #
    #######################################################################

    if capsys is not None:
        out, _ = capsys.readouterr()

        print(out)
        # make sure out isn't empty
        outstr = "SimControl: Termination from criteria: check_particle_nums"
        assert outstr in out
    else:
        assert False, "This wasn't run in pytest, but it passed otherwise."


if __name__ == '__main__':
    test_capacitive_discharge_custom_criteria(None, "Run2D")
