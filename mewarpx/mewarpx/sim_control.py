""" Control code used to terminate simulation based on
a set of user defined criteria
"""
import logging

from mewarpx.mwxrun import mwxrun

logger = logging.getLogger(__name__)


class SimControl:
    """
    Evaluate whether to continue or terminate the simulation based on a
    set of user defined functions and criteria

    """
    def __init__(self, max_steps, criteria=None):
        """
        Generate and install functions to perform after a step #.

        Arguments:

            max_steps (int): Maximum number of steps to perform in the
                simulation, when step == max_steps check_criteria returns False

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
        self.max_steps = max_steps
        if self.max_steps < 1:
            raise AssertionError("max_steps must be >= 1")
        self.add_checker(self.eval_max_steps)
        self.initialize_criteria(criteria)

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

    def eval_max_steps(self):
        if mwxrun.get_it() >= self.max_steps:
            logger.info("SimControl: Max steps reached!")
            return False
        return True

    def check_criteria(self):
        """Return True if all installed criteria are True, else False."""
        all_good = True
        terminate_statement = 'SimControl: Termination from criteria: '
        for i, criteria in enumerate(self.crit_list):
            all_good = all_good and criteria(**self.crit_args_list[i])
            if not criteria(**self.crit_args_list[i]):
                add_statement = "{criteria_name} ".format(
                    criteria_name=criteria.__name__)
                terminate_statement += add_statement

        if not all_good:
            logger.info(terminate_statement)

        return all_good
