""" Control code used to terminate simulation based on
a set of user defined criteria
"""
import logging
import os

from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store import util as mwxutil

logger = logging.getLogger(__name__)


class SimControl:
    """
    Evaluate whether to continue or terminate the simulation based on a
    set of user defined functions and criteria

    """
    def __init__(self, total_steps, criteria=None):
        """
        Generate and install functions to perform after a step #.

        Arguments:

            total_steps (int): Total number of steps to perform in the
                simulation, when step == total_steps check_criteria returns
                False
            criteria (list): list of user defined functions or list of user
                defined tuples(function, kwargs) that each return a True or
                False value

        Functions:

            check_criteria: Evaluates each criteria in the crit_list and
            returns True if all criteria are True otherwise returns False

        Example Sim run::

            steps_per_loop = 10
            while (sim_control.check_criteria()):
                sim.step(steps_per_loop)
        """

        self.crit_list = []
        self.crit_args_list = []
        self._write_func = None
        self.total_steps = total_steps
        if self.total_steps < 1:
            raise AssertionError("total_steps must be >= 1")
        self.add_checker(self.eval_total_steps)
        self.initialize_criteria(criteria)

        mwxrun.simulation.max_steps = self.total_steps
        logger.info(f"Total simulation steps set to {self.total_steps}")

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

    def eval_total_steps(self):
        """Evaluate whether the total timesteps have been reached."""
        if mwxrun.get_it() >= self.total_steps:
            logger.info("SimControl: Total steps reached")
            return False
        return True

    def check_criteria(self):
        """Return True if all installed criteria are True, else False."""
        all_good = True
        terminate_statement = 'SimControl: Termination from criteria: '
        for i, criteria in enumerate(self.crit_list):
            continue_flag = criteria(**self.crit_args_list[i])
            all_good = all_good and continue_flag
            if not continue_flag:
                add_statement = "{criteria_name} ".format(
                    criteria_name=criteria.__name__)
                terminate_statement += add_statement

        if not all_good:
            logger.info(terminate_statement)

        return all_good

    def write_results(self):
        """Create results.txt file, and write to it if write_func is set.
        The file signifies that the simulation ran to completion."""
        results_string = ""
        if callable(self._write_func):
            results_string = self._write_func()

        if mwxrun.me == 0:
            mwxutil.mkdir_p(WarpXDiagnostic.DIAG_DIR)
            with open(os.path.join(WarpXDiagnostic.DIAG_DIR, "results.txt"), 'a') as results_file:
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
