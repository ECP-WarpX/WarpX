""" Control code used to terminate simulation based on
a set of user defined criteria
"""
import logging
import os
import signal

from pywarpx import callbacks

from mewarpx.diags_store.checkpoint_diagnostic import CheckPointDiagnostic
from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store import util as mwxutil

logger = logging.getLogger(__name__)


class SimControl(WarpXDiagnostic):
    """
    The main simulation driving class. It evaluates whether to continue or
    terminate the simulation based on a set of user defined criteria.

    """
    def __init__(self, total_steps, diag_steps, criteria=None,
                 checkpoint=False, dump_period=None, **kwargs):
        """
        Generate and install functions to perform after a step #.

        Arguments:

            total_steps (int): Total number of steps to perform in the
                simulation, when step == total_steps check_criteria returns
                False
            diag_steps (int): Steps between diagnostic output
            criteria (list): list of user defined functions or list of user
                defined tuples(function, kwargs) that each return a True or
                False value
            checkpoint (bool): Whether or not simulation checkpoints should be
                made
            dump_period (int): Checkpoints will be created every
                diag_steps*dump_period steps. If the dump_period is not given,
                checkpoints will only be created when the simulation ends (
                including interrupts due to a TERMINATE signal)
        """
        if total_steps < 1:
            raise AssertionError("total_steps must be >= 1")
        # ensure that simulation runs through last step specified
        mwxrun.simulation.max_steps = int(total_steps)

        self.crit_list = []
        self.crit_args_list = []
        self._write_func = None

        self.diag_steps = diag_steps
        self.checkpoint = checkpoint

        self.initialize_criteria(criteria)

        self.sim_done = False

        super(SimControl, self).__init__(diag_steps=diag_steps)

        logger.info(
            f"Diagnostic time set to {self.diag_steps*mwxrun.dt:.2e} "
            f"({self.diag_steps} steps)."
        )
        logger.info(
            f"Total simulation time set to {total_steps*mwxrun.dt:.2e} "
            f"({total_steps} steps)."
        )

        # register the TERM signal to be handled in C++
        mwxrun.simulation.break_signals = "SIGTERM"

        # Install checkpointing if required
        if self.checkpoint:
            if dump_period is None:
                dump_period = total_steps
            else:
                dump_period = self.diag_steps*dump_period
            self.checkpoint_diag = CheckPointDiagnostic(dump_period, **kwargs)

        # install a callback to check whether any termination criteria is met
        callbacks.installafterstep(self.check_criteria)

    def add_checker(self, criterion):
        """Install a single function to check.

        Arguments:
            criterion (func or tuple): Either a function or a tuple of (func,
                kwargs_dict) where the kwargs_dict will be passed to the func.
        """
        if callable(criterion):
            self.crit_list.append(criterion)
            self.crit_args_list.append({})
        else:
            self.crit_list.append(criterion[0])
            self.crit_args_list.append(criterion[1])

    def initialize_criteria(self, criteria_list):
        """Install the full initial list of criteria."""
        if criteria_list:
            for crit in criteria_list:
                self.add_checker(crit)

    def check_criteria(self):
        """Sends an interrupt signal to the current process if a termination
        criteria is satisfied."""
        if self.check_timestep():
            terminate_statement = 'SimControl: Termination from criteria: '
            for i, criteria in enumerate(self.crit_list):
                continue_flag = criteria(**self.crit_args_list[i])
                self.sim_done = self.sim_done or not continue_flag
                if not continue_flag:
                    add_statement = f"{criteria.__name__} "
                    terminate_statement += add_statement

            if self.sim_done:
                logger.info(terminate_statement)
                os.kill(os.getpid(), signal.SIGTERM)

    def write_results(self):
        """Create results.txt file, and write to it if write_func is set.
        The file signifies that the simulation ran to completion."""
        results_string = ""
        if callable(self._write_func):
            results_string = self._write_func()

        if mwxrun.me == 0:
            mwxutil.mkdir_p(WarpXDiagnostic.DIAG_DIR)
            with open(
                os.path.join(WarpXDiagnostic.DIAG_DIR, "results.txt"), 'a'
            ) as results_file:
                results_file.write(results_string)

    def set_write_func(self, func):
        """Sets a function for writing to results.txt file.

        Arguments:
            func (function): Returns a string that will be written to the
            results.txt file
        """
        if not callable(func):
            raise ValueError("The write func is not callable")
        self._write_func = func

    def run(self):
        """Executes the WarpX loop."""
        mwxrun.simulation.step()

        # check if the simulation completed the total number of steps
        if mwxrun.get_it() >= mwxrun.simulation.max_steps:
            self.sim_done = True
            logger.info("SimControl: Total steps reached.")

        if self.sim_done:
            # create file to signal that simulation ran to completion
            self.write_results()
        else:
            # create fluxdiag checkpoint file if simulation was interrupted
            if self.checkpoint:
                self.checkpoint_diag.checkpoint_manager(mwxrun.get_it())
