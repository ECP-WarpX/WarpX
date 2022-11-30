#!/usr/bin/env python3

import os
import sys

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Get name of the test
test_name = os.path.split(os.getcwd())[1]

# Run checksum regression test
checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-2)
