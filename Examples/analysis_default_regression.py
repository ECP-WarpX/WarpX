#!/usr/bin/env python3

import os
import re
import sys

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the output file
plotfile = sys.argv[1]
opmdfile = './diags/diag2'

# Get name of the test
test_name = os.path.split(os.getcwd())[1]

# Run checksum regression test
if re.search( 'single_precision', test_name ):
    checksumAPI.evaluate_checksum(test_name, output_file=plotfile, output_format='plotfile', rtol=2.e-6)
    checksumAPI.evaluate_checksum(test_name, output_file=opmdfile, output_format='openpmd', rtol=2.e-6)
else:
    checksumAPI.evaluate_checksum(test_name, output_file=plotfile, output_format='plotfile')
    checksumAPI.evaluate_checksum(test_name, output_file=opmdfile, output_format='openpmd')
