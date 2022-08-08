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
                 clear_old_checkpoints=True, num_to_keep=2, **kwargs):
        """
        This class is a wrapper for creating checkpoints from which a
        simulation can be restarted. Adding flux diagnostic data to
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
            num_to_keep (int): Number of checkpoints to keep. Default 1.
            kwargs: For a list of valid keyword arguments see
                :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
        """
        self.checkpoint_steps = diag_steps
        self.name = name
        self.clear_old_checkpoints = clear_old_checkpoints
        self.num_to_keep = num_to_keep
        self.flux_diag = None

        super(CheckPointDiagnostic, self).__init__(
            diag_steps=diag_steps, **kwargs)

        self.write_dir = self.DIAG_DIR
        self.add_checkpoint()

        # if checkpoints will only be created with an interrupt signal or
        # the end of the simulation, we don't need to install the callback
        if self.checkpoint_steps != mwxrun.simulation.max_steps:
            callbacks.installafterdiagnostics(self.checkpoint_manager)

    def add_checkpoint(self):
        diagnostic = picmi.Checkpoint(
            period=self.checkpoint_steps,
            name=self.name,
            write_dir=self.write_dir
        )
        mwxrun.simulation.add_diagnostic(diagnostic)

    def checkpoint_manager(self, force_run=False):
        """Function executed on checkpoint steps to perform various tasks
        related to checkpoint management. These include copying the flux
        diagnostic data needed for a restart as well as deleting old
        checkpoints.
        """
        if not force_run and not self.check_timestep():
            return

        # Save a copy of flux diagnostics, if present, to load when restarting.
        if self.flux_diag is not None:
            # If the timeseries were not updated on timestep, do so now.
            if self.flux_diag.last_run_step != mwxrun.get_it():
                self.flux_diag.update_ts_dict()
                self.flux_diag.last_run_step = mwxrun.get_it()
                if mwxrun.me == 0:
                    self.flux_diag.update_fullhist_dict()

            if mwxrun.me == 0:
                # We use the unorthodox file extension .ckpt (for checkpoint)
                # so that we can continue to blindly move all .dpkl files from
                # EFS to S3 when running on AWS
                dst = os.path.join(
                    self.write_dir,
                    f"{self.name}{self.flux_diag.last_run_step:06d}",
                    "fluxdata.ckpt"
                )
                self.flux_diag.save(filepath=dst)

        if self.clear_old_checkpoints and mwxrun.me == 0:
            init_restart_util.clean_old_checkpoints(
                checkpoint_prefix=self.name, num_to_keep=self.num_to_keep
            )
