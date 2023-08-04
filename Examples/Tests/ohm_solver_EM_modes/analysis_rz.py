#!/usr/bin/env python3
#
# --- Analysis script for the hybrid-PIC example producing EM modes.

import dill
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from pywarpx import picmi

constants = picmi.constants

matplotlib.rcParams.update({'font.size': 20})

# load simulation parameters
with open(f'sim_parameters.dpkl', 'rb') as f:
    sim = dill.load(f)

# TODO add mode analysis logic

if sim.test:
    import os
    import sys
    sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
    import checksumAPI

    # this will be the name of the plot file
    fn = sys.argv[1]
    test_name = os.path.split(os.getcwd())[1]
    checksumAPI.evaluate_checksum(test_name, fn)
