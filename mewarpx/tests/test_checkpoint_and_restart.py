import os
import pytest
import logging

from mewarpx.utils_store import util as mwxutil


def test_create_checkpoints():

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun
    from mewarpx.diags_store.checkpoint_diagnostic import CheckPointDiagnostic

    DIAG_STEPS = 2
    CHECKPOINT_PERIOD = 5
    CHECKPOINT_NAME = "checkpoint"
    D_CA = 0.067  # m
    NX = 64
    NZ = 64
    DT = 7.5e-10

    MAX_STEPS = 10

    nx = 64
    ny = 64

    run = diode_setup.DiodeRun_V1(
        dim=2,
        rz=False,
        D_CA=D_CA,
        NX=nx,
        NZ=ny,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=MAX_STEPS,
        DIAG_STEPS=DIAG_STEPS,
        FIELD_DIAG_DATA_LIST=['phi']
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_injectors=False,
        init_inert_gas=False,
        init_neutral_plasma=False,
        init_mcc=False,
        init_field_diag=True,
        init_simcontrol=False,
        init_warpx=False
    )

    checkpoint = CheckPointDiagnostic(CHECKPOINT_PERIOD, CHECKPOINT_NAME)

    mwxrun.simulation.initialize_inputs()
    mwxrun.simulation.initialize_warpx()

    # Run the main WARP loop
    mwxrun.simulation.step(MAX_STEPS)

    checkpoint_names = [f"{CHECKPOINT_NAME}{i:05}" for i in range(0, MAX_STEPS + 1, CHECKPOINT_PERIOD)]

    for name in checkpoint_names:
        print(f"Looking for checkpoint file {name}...")
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
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    DIAG_STEPS = 2
    D_CA = 0.067  # m
    NX = 64
    NZ = 64
    DT = 7.5e-10

    MAX_STEPS = 10

    nx = 64
    ny = 64

    if not files_exist:
        prefix = "nonexistent_prefix"
    else:
        prefix = "checkpoint"

    run = diode_setup.DiodeRun_V1(
        dim=2,
        rz=False,
        D_CA=D_CA,
        NX=nx,
        NZ=ny,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=MAX_STEPS,
        DIAG_STEPS=DIAG_STEPS,
        FIELD_DIAG_DATA_LIST=['phi']
    )

    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_injectors=False,
        init_inert_gas=False,
        init_neutral_plasma=False,
        init_mcc=False,
        init_field_diag=True,
        init_simcontrol=False,
        init_warpx=False
    )
    if force:
        restart = True
    else:
        restart = None
    try:
        mwxrun.init_run(
            restart=restart,
            checkpoint_dir="tests/test_files/checkpoint",
            checkpoint_prefix=prefix
        )
    except RuntimeError as e:
        if "There were no checkpoint directories starting with nonexistent_prefix!" in str(e):
            assert force # There should only be an exception if this restart was forced

    # if the files didn't exist then we didn't restart, so there's no need to verify the correct number of steps passed
    if files_exist:
        new_max_steps = mwxrun.simulation.max_steps

        start_step = mwxrun.get_it()

        if not force and not files_exist:
            log_lines = [r.msg for r in caplog.records]
            assert any(
                [f"There were no checkpoint directories starting with {prefix}!" in l for l in log_lines]
            )

        mwxrun.simulation.step()

        end_step = mwxrun.get_it()

        assert end_step - start_step == new_max_steps

def test_extra_steps_after_restart():
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    DIAG_STEPS = 2
    D_CA = 0.067  # m
    NX = 64
    NZ = 64
    DT = 7.5e-10

    MAX_STEPS = 10

    nx = 64
    ny = 64

    run = diode_setup.DiodeRun_V1(
        dim=2,
        rz=False,
        D_CA=D_CA,
        NX=nx,
        NZ=ny,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=MAX_STEPS,
        DIAG_STEPS=DIAG_STEPS,
        FIELD_DIAG_DATA_LIST=['phi']
    )

    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_injectors=False,
        init_inert_gas=False,
        init_neutral_plasma=False,
        init_mcc=False,
        init_field_diag=True,
        init_simcontrol=False,
        init_warpx=False
    )
    additional_steps = 20

    mwxrun.init_run(
        restart=True,
        checkpoint_dir="tests/test_files/checkpoint",
        additional_steps=additional_steps
    )

    start_step = mwxrun.get_it()

    mwxrun.simulation.step()

    end_step = mwxrun.get_it()

    assert start_step + additional_steps == end_step
