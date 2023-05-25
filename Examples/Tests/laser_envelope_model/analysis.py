import os
import sys

import numpy as np
import scipy.constants as scc

import yt ; yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
