""" Control code used to terminate simulation based on
a set of user defined criteria
"""

import logging
import time

import numpy as np
import psutil

from pywarpx import _libwarpx
from mewarpx.mwxrun import mwxrun

# Get module-level logger
logger = logging.getLogger(__name__)

class SimControl:
    """
    Evaluate whether to continue or terminate the simulation based on a
    set of user defined functions and criteria

    """
    def __init__(self, criteria, criteria_args, max_loops=None, **kwargs):
        """
        Generate and install functions to perform after a step #.
        Arguments:

        criteria: list of user defined functions that return True or False value

        criteria_args: list of kwargs for each criteria
        """

        self.crit_list = criteria
        self.crit_args_list = criteria_args
        self.max_loops = max_loops
        self.check_increment = 0

    def add_checker(self, criteria=None, crit_args=None):
        self.crit_list.append(criteria)
        if crit_args:
            self.crit_args_list.append(crit_args)
        else:
            self.crit_args_list.append({})

    def check_criteria(self):
        """
        Evaluates each criteria in the crit_list
        returns True is all criteria are True otherwise returns False

        Example Sim run:

        steps_per_loop = 10
        while (sim_control.check_criteria()):
            sim.step(steps_per_loop)

        """
        if self.check_increment > 0:
            keep_going = False
            for i, func in enumerate(self.crit_list):
                keep_going = self.crit_list[i](**self.crit_args_list[i])
                if not keep_going:
                    print(("Termination due to criteria {}").format(func.__name__))
                    return keep_going
            self.check_increment += 1
            if self.max_loops:
                if self.check_increment >= self.max_loops:
                    keep_going = False
            return keep_going
        else:
            self.check_increment += 1
            return True
