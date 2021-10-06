"""
Utility functions to start a run from a checkpoint or restart.
"""
import os
import shutil
import logging

logger = logging.getLogger(__name__)


default_checkpoint_name = "checkpoint"

def clean_old_checkpoints(prefix, directory):
    if os.path.isdir(directory):
        for d in next(os.walk(directory))[1]:
            if d.startswith(prefix):
                print(f"Removing old checkpoint file {d}")
                shutil.rmtree(os.path.join(directory, d))

def run_restart(checkpoint_directory="diags", checkpoint_prefix=default_checkpoint_name,
                force=False, additional_steps=None):
    '''
    Attempts to restart a run by looking for checkpoint files
    starting with a prefix in the given directory.

    Arguments:
        checkpoint_directory (string): Look in this directory for checkpoint directories. Default is 'diags'.
        checkpoint_prefix (string): Look for a checkpoint directory starting with this prefix
          to restart from.
        force (bool): If true, a problem with restarting from a checkpoint will cause an error,
          otherwise simply print a warning.
        additional_steps (int): The number of steps to run after restarting from the checkpoint.
          If this is None then it will run to the current value of mwxrun.simulation.max_steps.
    '''
    # import must be done here to avoid a circular import
    from mewarpx.mwxrun import mwxrun

    logger.info(
        "Attempting to " + ("force a " if force else "") +
        f"restart from the most recent checkpoint in {checkpoint_directory} "
        f"starting with '{checkpoint_prefix}'"
    )

    if not os.path.isdir(checkpoint_directory):
        if force:
            raise RuntimeError(f"{checkpoint_directory} directory does not exist!")
        else:
            logger.warning(f"{checkpoint_directory} directory does not exist!")
            return False

    # using next() gives only the first layer of subdirectories
    checkpoints = sorted(
        [f for f in next(os.walk(checkpoint_directory))[1]
        if "old" not in f and f.startswith(checkpoint_prefix)]
    )

    if not checkpoints:
        if force:
            raise RuntimeError(
                "There were no checkpoint directories "
                f"starting with {checkpoint_prefix}!"
            )
        else:
            logger.warning(
                "There were no checkpoint directories "
                f"starting with {checkpoint_prefix}!"
            )
            return False

    checkpoint = checkpoints[-1]
    max_steps = mwxrun.simulation.max_steps
    checkpoint_step = int(checkpoint.replace(checkpoint_prefix, ""))

    if additional_steps is None:
        mwxrun.simulation.max_steps = max_steps - checkpoint_step
        if mwxrun.simulation.max_steps == 0:
            logger.warning(
                f"The checkpoint directory was created at step {checkpoint_step}, "
                f"but the max steps is also {max_steps}, so the simulation will "
                f"only rerun step {checkpoint_step}."
            )
        if mwxrun.simulation.max_steps < 0:
            raise RuntimeError(
                "The checkpoint directory was created at a later step "
                f"({checkpoint_step}) than the current max steps ({max_steps})!"
            )
        logger.info(f"Running until step {max_steps}")
    else:
        mwxrun.simulation.max_steps = additional_steps
        logger.info(f"Running for {additional_steps} steps after restarting")


    logger.info(f"Restarting from {checkpoint}")

    mwxrun.simulation.amr_restart = os.path.join(checkpoint_directory, checkpoint)
    return True
