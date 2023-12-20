#!/usr/bin/env python3

# Copyright 2019-2023
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs_1d`. This simulates a 1D WFA with Pondermotive Envelope:
# REF: (Equations 20-23) https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.81.1229
import os
import sys

import openpmd_api as io

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
