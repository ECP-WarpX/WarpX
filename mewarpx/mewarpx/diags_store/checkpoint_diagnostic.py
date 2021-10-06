"""Class for installing a checkpoint diagnostic"""

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx.utils_store import init_restart_util

from pywarpx import picmi

import logging

logger = logging.getLogger(__name__)


class CheckPointDiagnostic(WarpXDiagnostic):
    def __init__(self, diag_steps, name=init_restart_util.default_checkpoint_name,
                 write_dir=".", **kwargs):
        """
        This class is a wrapper for creating checkpoints from which
        to restart a simulation from

        Arguments:
            diag_steps (int): Run the diagnostic with this period.
                Also plot on this period if enabled.
            name (str): The name of the diagnostic to be passed into the
                picmi checkpoint diagnostic.
            write_dir (str): The directory where diagnostic data will be written
                from the checkpoint diagnostic.
            kwargs: For a list of valid keyword arguments see diag_base.WarpXDiagnostic
        """
        self.diag_steps = diag_steps
        self.name = name
        self.write_dir = write_dir

        super(CheckPointDiagnostic, self).__init__(diag_steps, **kwargs)

        self.add_checkpoint()

    def add_checkpoint(self):
        diagnostic = picmi.Checkpoint(
            period=self.diag_steps,
            name=self.name,
            write_dir=self.write_dir
        )

        mwxrun.simulation.add_diagnostic(diagnostic)
