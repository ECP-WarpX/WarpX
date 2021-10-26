#! /usr/bin/env python

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
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
