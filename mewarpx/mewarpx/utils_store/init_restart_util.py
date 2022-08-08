"""
Utility functions to start a run from a checkpoint or restart.
"""
import logging
import os
import shutil

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
                          num_to_keep=2):
    """Utility function to remove old checkpoints.

    Arguments:
        checkpoint_directory (str): Look in this directory for checkpoint
            directories. Default is ``diags``.
        checkpoint_prefix (str): Look for a checkpoint directory starting
            with this prefix to restart from.
        num_to_keep (int): Keep this many of the newest checkpoints. Default 2,
            so that one being judged corrupt will never ruin all checkpoints.
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


def _eval_checkpoint_validity(checkpoint_dir, checkpoint):
    """Determine if a checkpoint appears to be valid by checking if
    fluxdata.ckpt is present.

    Arguments:
        checkpoint_dir (str): Look in this directory for checkpoint
            directories.
        checkpoint (str): Checkpoint folder

    Returns:
        checkpoint_ok (bool): True if the checkpoint appears good
        (fluxdata.ckpt is present), False if fluxdata.ckpt is missing.
    """
    if not os.path.isfile(
        os.path.join(checkpoint_dir, checkpoint, "fluxdata.ckpt")
    ):
        logger.warning(
            f"Checkpoint {checkpoint} does not contain a flux diag "
            "checkpoint."
        )
        return False

    return True


def _remove_corrupt_checkpoint(checkpoint_dir, checkpoint):
    """Remove a checkpoint that has already been determined to be corrupt and
    appropriate to remove.

    Arguments:
        checkpoint_dir (str): Look in this directory for checkpoint
            directories.
        checkpoint (str): Checkpoint folder
    """
    dirpath = os.path.join(checkpoint_dir, checkpoint)
    logger.info(f"Removing corrupt checkpoint {dirpath}")
    shutil.rmtree(dirpath)


def _handle_corrupt_checkpoints(checkpoint_dir, checkpoint_list):
    """Utility function to remove the last checkpoint if it appears to be
    corrupt. Error out if two or more checkpoints are corrupt.

    Note:
        A checkpoint is judged corrupt if a flux diag checkpoint does not exist
        (which is always written after the checkpoint). If runs are ever used
        without flux diagnostics this logic should be adapted.

    Arguments:
        checkpoint_dir (str): Look in this directory for checkpoint
            directories.
        checkpoint_list (list): List of the checkpoints returned by
            get_sorted_checkpoints.

    Returns:
        checkpoint_list (list): The same checkpoint_list, but if the final
        checkpoint was corrupt it will be removed. If multiple appear to be
        corrupt an error is raised instead of returning.
    """
    last_checkpoint_ok = _eval_checkpoint_validity(
        checkpoint_dir, checkpoint_list[-1]
    )
    # Short-circuit the corrupt handling logic: We don't have an issue, carry
    # on as before.
    if last_checkpoint_ok:
        return checkpoint_list

    # If there's only one checkpoint and we know it's corrupt, remove it and
    # then execute the non-restart logic by returning an empty checkpoint
    # list.
    if len(checkpoint_list) == 1:
        _remove_corrupt_checkpoint(checkpoint_dir, checkpoint_list[-1])
        return []

    # If there are multiple checkpoints, our action depends on whether only the
    # final one is corrupt
    penultimate_checkpoint_ok = _eval_checkpoint_validity(
        checkpoint_dir, checkpoint_list[-2]
    )
    # If only the final one is corrupt, we remove it and carry on
    if penultimate_checkpoint_ok:
        _remove_corrupt_checkpoint(checkpoint_dir, checkpoint_list[-1])
        return checkpoint_list[:-1]

    # If multiple are corrupt, we raise an error and don't do anything
    raise RuntimeError(
        "Multiple checkpoints lacked fluxdata.ckpt, indicating they are "
        "corrupt. This should never occur, so the simulation is terminating "
        "now."
    )


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

    if checkpoints:
        checkpoints = _handle_corrupt_checkpoints(
            checkpoint_directory, checkpoints
        )

    # Note this can't be an else clause! checkpoints can be changed by
    # _handle_corrupt_checkpoints
    if not checkpoints:
        if force:
            raise RuntimeError(
                "There were no valid checkpoint directories "
                f"starting with {checkpoint_prefix}!"
            )
        else:
            logger.warning(
                "There were no valid checkpoint directories "
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
