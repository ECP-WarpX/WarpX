""" Control code used to terminate simulation based on 
a set of user defined criteria
"""

import logging
import time

import numpy as np
import psutil

from pywarpx._libwarpx import get_particle_id
from mewarpx.mwxrun import mwxrun

# Get module-level logger
logger = logging.getLogger(__name__)

class SimControl:
    """
    Evaluate whether to continue or terminate the simulation based on a
    set of user defined functions and criteria

    """
    def __init__(self, criteria=None,
                criteria_args=None, **kwargs):
        """
        Generate and install functions to perform after a step #.
        Arguments:
        diag_steps (int): Number of steps between each output
        """
        self.crit_list = []
        self.crit_args_list = []
        self.check_increment = 0

    def add_checker(self, criteria=None, crit_args=None):
        self.crit_list.append(criteria)
        if crit_args:
            self.crit_args_list.append(crit_args)
        else:
            self.crit_args_list.append({})

    def check_criteria(self):
        if self.check_increment > 0:
            keep_going = False
            for i, func in enumerate(self.crit_list):
                keep_going = self.crit_list[i](**self.crit_args_list[i])
                if not keep_going:
                    print(("Termination due to criteria {}").format(func.__name__))
                    return keep_going
            return keep_going
        else:
            self.check_increment += 1
            return True


def particle_checker(species_number, minimum, maximum, level=0):
    particle_ids = get_particle_id(species_number, level)
    print("HELLO: "+str(len(particle_ids)))
    if (len(particle_ids) > minimum) & (len(particle_ids) < maximum):
        return True
    else:
        return False
