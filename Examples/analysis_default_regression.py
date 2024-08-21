#!/usr/bin/env python3

import os
import re
import sys

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Get name of the test
test_name = os.path.split(os.getcwd())[1]

restart = "restart" in test_name

if restart:
    test_name = test_name.replace("_restart", "")

# Run checksum regression test
if re.search("single_precision", fn):
    checksumAPI.evaluate_checksum(test_name, fn, rtol=2.0e-6)
else:
    if restart:
        checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-12)
    else:
        # using default relative tolerance
        checksumAPI.evaluate_checksum(test_name, fn)
