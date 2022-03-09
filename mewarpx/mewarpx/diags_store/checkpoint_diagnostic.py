"""Class for installing a checkpoint diagnostic"""
import logging
import os

from pywarpx import callbacks, picmi

from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store import init_restart_util

logger = logging.getLogger(__name__)


class CheckPointDiagnostic(WarpXDiagnostic):
    def __init__(self, diag_steps,
                 name=init_restart_util.DEFAULT_CHECKPOINT_NAME,
                 clear_old_checkpoints=True, **kwargs):
        """
        This class is a wrapper for creating checkpoints from which
        to restart a simulation from. Adding flux diagnostic data to
        checkpoints are supported, but since this class has to be initialized
        before the simulation and flux diagnostics are initialized after
        the simulation, the user is responsible for adding the ``flux_diag``
        attribute to this object in a simulation input file.

        Arguments:
            diag_steps (int): Run the diagnostic with this period.
                Also plot on this period if enabled.
            name (str): The name of the diagnostic to be passed into the
                picmi checkpoint diagnostic.
            clear_old_checkpoints (bool): If True old checkpoints will be
                deleted after new ones are created.
            kwargs: For a list of valid keyword arguments see
                :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
        """
        self.checkpoint_steps = diag_steps
        self.name = name
        self.clear_old_checkpoints = clear_old_checkpoints
        self.flux_diag = None

        if kwargs.pop('diag_step_offset', 1) != 1:
            raise ValueError(
                "diag_step_offset in CheckPointDiagnostic must be 1.")
        super(CheckPointDiagnostic, self).__init__(
            diag_steps=diag_steps, diag_step_offset=1, **kwargs)

        self.write_dir = self.DIAG_DIR
        self.add_checkpoint()

        # note that WarpX diagnostics (including checkpoints) are output
        # after the "afterstep" functions, therefore we install this function
        # before steps and add 1 to the diag_steps above
        callbacks.installbeforestep(self.checkpoint_manager)

    def add_checkpoint(self):
        diagnostic = picmi.Checkpoint(
            period=self.checkpoint_steps,
            name=self.name,
            write_dir=self.write_dir
        )

        mwxrun.simulation.add_diagnostic(diagnostic)

    def checkpoint_manager(self):
        """Function executed on checkpoint steps to perform various tasks
        related to checkpoint management. These include copying the flux
        diagnostic data needed for a restart as well as deleting old
        checkpoints.
        """
        # only the root processor should execute this function but not on the
        # first step
        if not self.check_timestep() or mwxrun.me != 0 or mwxrun.get_it() == 1:
            return

        # Save a copy of flux diagnostics, if present, to load when restarting.
        if (
            self.flux_diag is not None
            and self.flux_diag.last_run_step == mwxrun.get_it() - 1
        ):
            # We use the unorthodox file extension .ckpt (for checkpoint) so
            # that we can continue to blindly move all .dpkl files from EFS
            # to S3 when running on AWS
            dst = os.path.join(
                self.write_dir, f"{self.name}{(mwxrun.get_it() - 1):06d}",
                "fluxdata.ckpt"
            )
            self.flux_diag.save(filepath=dst)
        else:
            logger.warning(
                "Flux diagnostic data will not be saved with checkpoint and "
                "therefore won't be available at restart."
            )

        if self.clear_old_checkpoints:
            init_restart_util.clean_old_checkpoints(
                checkpoint_prefix=self.name, num_to_keep=1
            )
