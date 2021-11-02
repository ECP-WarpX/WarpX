import os
import pytest
import logging
import numpy as np

from mewarpx.utils_store import util as mwxutil

VOLTAGE = 25 # V
CATHODE_TEMP = 1100 + 273.15 # K
CATHODE_PHI = 2.1 # work function in eV

DIAG_STEPS = 4
MAX_STEPS = 8
CHECKPOINT_NAME = "checkpoint"
D_CA = 0.067  # m
NX = 64
NZ = 64
DT = 7.5e-10

P_INERT = 10

def get_run():
    """Utility function to get the same run setup for all tests below."""
    from mewarpx.setups_store import diode_setup

    # Initialize and import only when we know dimension
    run = diode_setup.DiodeRun_V1(
        dim=2,
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NPPC=4,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=MAX_STEPS,
        DIAG_STEPS=DIAG_STEPS,
        P_INERT=P_INERT,
        T_ELEC=1100 + 273.15,
        INERT_GAS_TYPE="Ar",
        FIELD_DIAG_DATA_LIST=['phi'],
        PLASMA_DENSITY=10.2e12,
        SEED_NPPC=10
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_injectors=False,
        init_inert_gas=True,
        init_neutral_plasma=True,
        init_mcc=False,
        init_field_diag=True,
        init_simcontrol=False,
        init_warpx=False
    )
    return run

def test_create_checkpoints():

    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.mwxrun import mwxrun
    from mewarpx.diags_store.checkpoint_diagnostic import CheckPointDiagnostic
    from mewarpx.utils_store import testing_util

    testing_util.initialize_testingdir("test_create_checkpoints")

    # use a fixed random seed
    np.random.seed(47239475)

    run = get_run()

    checkpoint = CheckPointDiagnostic(
        DIAG_STEPS, CHECKPOINT_NAME, write_dir='diags'
    )

    mwxrun.init_run(restart=False)

    # Run the main WARP loop
    mwxrun.simulation.step(MAX_STEPS)

    checkpoint_names = [
        f"{CHECKPOINT_NAME}{i:05}"
        for i in range(0, MAX_STEPS + 1, DIAG_STEPS)
    ]

    for name in checkpoint_names:
        print(f"Looking for checkpoint file 'diags/{name}'...")
        assert os.path.isdir(os.path.join("diags", name))


@pytest.mark.parametrize("force, files_exist",
[
    (True, True), # should restart from the newest checkpoint
    (False, True), # should restart from the newest checkpoint
    (True, False), # should throw an error and not start a new run
    (False, False), # should throw an error but start a fresh run
])
def test_restart_from_checkpoint(caplog, force, files_exist):

    caplog.set_level(logging.WARNING)
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.mwxrun import mwxrun
    from mewarpx.utils_store import testing_util

    testing_util.initialize_testingdir(
        f"test_restart_from_checkpoint_{force}_{files_exist}"
    )

    run = get_run()

    if force:
        restart = True
    else:
        restart = None

    if not files_exist:
        prefix = "nonexistent_prefix"
    else:
        prefix = "checkpoint"

    try:
        mwxrun.init_run(
            restart=restart,
            checkpoint_dir=os.path.join(testing_util.test_dir, "checkpoint"),
            checkpoint_prefix=prefix
        )
    except RuntimeError as e:
        if (
            "There were no checkpoint directories starting with "
            "nonexistent_prefix!" in str(e)
        ):
            # There should only be an exception if this restart was forced
            assert force

    # if the files didn't exist then we didn't restart,
    # so there's no need to verify the correct number of steps passed
    if files_exist:
        new_max_steps = mwxrun.simulation.max_steps

        start_step = mwxrun.get_it()

        if not force and not files_exist:
            log_lines = [r.msg for r in caplog.records]
            assert any(
                [f"There were no checkpoint directories starting with {prefix}!"
                in l for l in log_lines]
            )

        mwxrun.simulation.step()

        end_step = mwxrun.get_it()

        assert end_step - start_step == new_max_steps

def test_extra_steps_after_restart():

    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.mwxrun import mwxrun
    from mewarpx.utils_store import testing_util

    testing_util.initialize_testingdir("test_extra_steps_after_restart")

    # use a fixed random seed
    np.random.seed(47239475)

    run = get_run()

    additional_steps = 8

    mwxrun.init_run(
        restart=True,
        checkpoint_dir=os.path.join(testing_util.test_dir, "checkpoint"),
        additional_steps=additional_steps
    )

    start_step = mwxrun.get_it()

    mwxrun.simulation.step()

    end_step = mwxrun.get_it()

    assert start_step + additional_steps == end_step

    restart_net_charge_density = np.load(os.path.join(
        run.field_diag.write_dir, "Net_charge_density_0000000008.npy"
    ))
    original_net_charge_density = np.load(os.path.join(
        testing_util.test_dir, "checkpoint", "Net_charge_density_0000000008.npy"
    ))

    assert np.allclose(restart_net_charge_density,
                       original_net_charge_density, rtol=0.1)

def test_checkpoints_fluxdiag():

    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.mwxrun import mwxrun
    from mewarpx.diags_store.checkpoint_diagnostic import CheckPointDiagnostic
    from mewarpx.diags_store.flux_diagnostic import FluxDiagFromFile
    from mewarpx.utils_store import testing_util

    testing_util.initialize_testingdir("test_checkpoints_fluxdiag")

    # use a fixed random seed
    np.random.seed(47239475)

    run = get_run()

    run.init_injectors()

    mwxrun.init_run(
        restart=True,
        checkpoint_dir=os.path.join(testing_util.test_dir, "checkpoint"),
    )

    run.init_runinfo()
    run.init_fluxdiag()

    # Run the main WARP loop
    mwxrun.simulation.step()

    # load flux diagnostic from a completed run for comparison
    basedir = os.path.join(testing_util.test_dir, "checkpoint")
    original_flux_file = os.path.join(
        testing_util.test_dir, "checkpoint/fluxes/fluxdata_0000000008.dpkl"
    )
    original_flux = FluxDiagFromFile(basedir, original_flux_file)

    # compare injected flux from the restarted history with the original history
    flux_key = ('inject', 'cathode', 'electrons')
    for key in ['J', 'P', 'n']:
        old = original_flux.fullhist_dict[flux_key].get_timeseries_by_key(key)
        new = run.fluxdiag.fullhist_dict[flux_key].get_timeseries_by_key(key)
        assert np.allclose(old, new)
