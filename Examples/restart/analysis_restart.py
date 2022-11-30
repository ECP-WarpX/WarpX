#!/usr/bin/env python3

import os
import sys

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# Check restart data v. original data
sys.path.insert(0, '../../../../warpx/Examples/')
from analysis_default_restart import check_restart

check_restart(filename)

# Check-sum analysis
filename = sys.argv[1]
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
