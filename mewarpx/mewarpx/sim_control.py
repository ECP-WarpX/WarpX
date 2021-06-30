""" Control code used to terminate simulation based on
a set of user defined criteria
"""
from mewarpx.mwxrun import mwxrun


class SimControl:
    """
    Evaluate whether to continue or terminate the simulation based on a
    set of user defined functions and criteria

    """
    def __init__(self, max_loops, criteria=None, criteria_args=None, **kwargs):
        """
        Generate and install functions to perform after a step #.
        Arguments:

        criteria: list of user defined functions that return True or False value

        criteria_args: list of kwargs for each criteria
        """

        self.crit_list = criteria
        self.crit_args_list = criteria_args
        self.max_loops = max_loops
        if not self.max_loops >= 1:
            raise AssertionError("max_loops must be >= 1")
        self.check_increment = 0

    def add_checker(self, criteria, crit_args=None):
        self.crit_list.append(criteria)
        if crit_args:
            self.crit_args_list.append(crit_args)
        else:
            self.crit_list.append(crit_args)

    def check_criteria(self):
        """
        Evaluates each criteria in the crit_list
        returns True if all criteria are True otherwise returns False

        Example Sim run:

        steps_per_loop = 10
        while (sim_control.check_criteria()):
            sim.step(steps_per_loop)

        """

        self.check_increment += 1
        if self.check_increment >= self.max_loops:
            print('Max loops reached!')

        for i, func in enumerate(self.crit_list):
            if not self.crit_list[i](**self.crit_args_list[i]):
                print(('Termination due to criteria {}')).format(func.__name__)
                return False

        return True
