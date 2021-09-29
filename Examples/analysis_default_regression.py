#! /usr/bin/env python

import sys
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# Get name of the test
test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]

# Run checksum regression test
if re.search( 'single_precision', fn ):
    checksumAPI.evaluate_checksum(test_name, fn, rtol=1.e-7)
else:
    checksumAPI.evaluate_checksum(test_name, fn)

