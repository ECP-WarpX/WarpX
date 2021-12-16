"""
Utility functions to start a run from a checkpoint or restart.
"""
import os
import shutil
import logging

logger = logging.getLogger(__name__)

DEFAULT_CHECKPOINT_NAME = "checkpoint"


def get_sorted_checkpoints(checkpoint_directory, checkpoint_prefix):
    """Function to return valid checkpoints for restarting or removing
    checkpoints.

    Also warns if any checkpoints with .old. exist, which shouldn't in normal
    workflow.

    Arguments:
        checkpoint_directory (str): Look in this directory for checkpoint
            directories. Default is ``diags``.
        checkpoint_prefix (str): Look for a checkpoint directory starting
            with this prefix to restart from.
    Returns:
        checkpoint_dirnames (list of str): List of checkpoint directory
        strings, sorted by timestep from earliest to latest.
    """
    if ".old." in checkpoint_prefix:
        raise RuntimeError(
            f".old. in checkpoint_prefix {checkpoint_prefix} will break "
            "restarts."
        )

    # Warn on checkpoints with old in name
    # Using next() gives only the first layer of subdirectories
    checkpoints_old = [
        f for f in next(os.walk(checkpoint_directory))[1]
        if ".old." in f and f.startswith(checkpoint_prefix)
    ]
    for d in checkpoints_old:
        logger.warning(
            f"A stale checkpoint file {d} exists, which is created when "
            "the same checkpoint step is saved again. Please inspect run "
            "history manually."
        )

    # Get list of valid checkpoints
    checkpoints = [
        f for f in next(os.walk(checkpoint_directory))[1]
        if ".old." not in f and f.startswith(checkpoint_prefix)
    ]

    # Naturally sort checkpoints by extracting step number and converting to
    # int
    checkpoints = sorted(
        checkpoints, key=lambda fname: int(fname.strip(checkpoint_prefix))
    )
    return checkpoints


def clean_old_checkpoints(checkpoint_directory="diags",
                          checkpoint_prefix=DEFAULT_CHECKPOINT_NAME,
                          num_to_keep=1):
    """Utility function to remove old checkpoints.

    Arguments:
        checkpoint_directory (str): Look in this directory for checkpoint
            directories. Default is ``diags``.
        checkpoint_prefix (str): Look for a checkpoint directory starting
            with this prefix to restart from.
        num_to_keep (int): Keep this many of the newest checkpoints. Default 1.
    """
    # handle the case where num_to_keep is 0 or None
    if not num_to_keep:
        num_to_keep = None
    else:
        num_to_keep *= -1

    # Handle potentially good checkpoints
    checkpoints = get_sorted_checkpoints(
        checkpoint_directory=checkpoint_directory,
        checkpoint_prefix=checkpoint_prefix
    )

    for d in checkpoints[:num_to_keep]:
        dirpath = os.path.join(checkpoint_directory, d)
        logger.info(f"Removing old checkpoint file {dirpath}")
        shutil.rmtree(dirpath)


def run_restart(checkpoint_directory="diags",
                checkpoint_prefix=DEFAULT_CHECKPOINT_NAME,
                force=False, additional_steps=None):
    """Attempts to restart a run by looking for checkpoint files
    starting with a prefix in the given directory.

    Arguments:
        checkpoint_directory (str): Look in this directory for checkpoint
            directories. Default is ``diags``.
        checkpoint_prefix (str): Look for a checkpoint directory starting
            with this prefix to restart from.
        force (bool): If true, a problem with restarting from a checkpoint will
            cause an error, otherwise simply print a warning.
        additional_steps (int): The number of steps to run after restarting from
            the checkpoint. If this is None then it will run to the current
            value of mwxrun.simulation.max_steps.
    """
    # import must be done here to avoid a circular import
    from mewarpx.mwxrun import mwxrun

    logger.info(
        "Attempting to " + ("force a " if force else "") +
        f"restart from the most recent checkpoint in {checkpoint_directory} "
        f"starting with '{checkpoint_prefix}'"
    )

    if not os.path.isdir(checkpoint_directory):
        if force:
            raise RuntimeError(
                f"{checkpoint_directory} directory does not exist!"
            )
        else:
            logger.warning(f"{checkpoint_directory} directory does not exist!")
            return False, None, None

    checkpoints = get_sorted_checkpoints(
        checkpoint_directory=checkpoint_directory,
        checkpoint_prefix=checkpoint_prefix
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
            return False, None, None

    checkpoint = checkpoints[-1]
    max_steps = mwxrun.simulation.max_steps
    checkpoint_step = int(checkpoint.replace(checkpoint_prefix, ""))

    if checkpoint_step == 0:
        return False, None, None

    logger.info(f"Restarting from {checkpoint}")

    if additional_steps is None:
        mwxrun.simulation.max_steps = max_steps - checkpoint_step
        if mwxrun.simulation.max_steps == 0:
            logger.warning(
                f"The checkpoint directory was created at step "
                f"{checkpoint_step}, but the max steps is also {max_steps}, so "
                f"the simulation will only rerun step {checkpoint_step}."
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

    mwxrun.simulation.amr_restart = os.path.join(
        checkpoint_directory, checkpoint
    )

    return True, checkpoint_directory, checkpoint
