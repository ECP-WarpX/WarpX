#! /usr/bin/env python

import sys
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import os
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Get name of the test
test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]

# Reset all benchmarks?
reset = ( os.getenv('CHECKSUM_RESET', 'False').lower() in
          ['true', '1', 't', 'y', 'yes', 'on'] )

# Run checksum regression test or reset
if reset:
    checksumAPI.reset_benchmark(test_name, fn)
else:
    checksumAPI.evaluate_checksum(test_name, fn)
