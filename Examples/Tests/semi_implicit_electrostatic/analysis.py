#!/usr/bin/env python3

import os
import sys

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

filename = sys.argv[1]

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
