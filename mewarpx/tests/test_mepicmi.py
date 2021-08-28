"""Test the functionality in mewarpx/mepicmy.py"""
import logging
import numpy as np

from mewarpx.utils_store import util as mwxutil


def test_extra_pid(caplog):
    caplog.set_level(logging.INFO)
    name = "mepicmi_extra_pid"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun
    from pywarpx import _libwarpx

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(9216001)

    # Specific numbers match older run for consistency
    DIAG_STEPS = 2
    D_CA = 0.05  # m
    NX = 16
    NZ = 128
    run = diode_setup.DiodeRun_V1(
        dim=dim,
        rz=False,
        V_ANODE_CATHODE=450.0,
        D_CA=D_CA,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=1.0e-10,
        TOTAL_TIMESTEPS=25,
        DIAG_STEPS=DIAG_STEPS,
        DIRECT_SOLVER=True
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_injectors=False,
        init_simcontrol=True,
        init_warpx=True
    )

    # add new pid for the ions
    run.electrons.add_pid('extra_pid')

    def check_particle_nums():
        return not (mwxrun.get_npart() < 950)

    nps = 1000
    w = np.random.randint(low=1, high=100, size=nps)
    _libwarpx.add_particles(
        run.electrons.name,
        x=np.random.random(nps) * D_CA / NZ * NX, y=np.zeros(nps),
        z=np.random.random(nps) * D_CA,
        ux=np.random.normal(scale=1e4, size=nps),
        uy=np.random.normal(scale=1e4, size=nps),
        uz=np.random.normal(scale=1e4, size=nps),
        w=w, extra_pid=np.copy(w)*10.0
    )

    run.control.add_checker(check_particle_nums)

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step(DIAG_STEPS)

    #######################################################################
    # Cleanup and final output                                            #
    #######################################################################

    all_log_output = ""
    records = caplog.records

    for record in records:
        all_log_output += record.msg + "\n"

    # make sure out isn't empty
    outstr = "SimControl: Termination from criteria: check_particle_nums"
    assert outstr in all_log_output

    weights = run.electrons.get_array_from_pid('w')
    extra_pid = run.electrons.get_array_from_pid('extra_pid')
    for ii in range(len(weights)):
        assert np.allclose(extra_pid[ii] / weights[ii], 10.0)
