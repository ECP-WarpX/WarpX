#!/usr/bin/env python3

# Copyright 2023 David Grote
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# This is a script that analyses the simulation results from
# the script `inputs_1d`. This simulates a 1D periodic plasma using the implicit solver.
import os
import re
import sys

import numpy as np

sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

field_energy = np.loadtxt("diags/reducedfiles/field_energy.txt", skiprows=1)
particle_energy = np.loadtxt("diags/reducedfiles/particle_energy.txt", skiprows=1)

total_energy = field_energy[:, 2] + particle_energy[:, 2]

delta_E = (total_energy - total_energy[0]) / total_energy[0]
max_delta_E = np.abs(delta_E).max()

test_name = os.path.split(os.getcwd())[1]
if re.match("semi_implicit_picard_1d", test_name):
    tolerance_rel = 2.5e-5
elif re.match("theta_implicit_picard_1d", test_name):
    # This case should have near machine precision conservation of energy
    tolerance_rel = 1.0e-14

print(f"max change in energy: {max_delta_E}")
print(f"tolerance: {tolerance_rel}")

assert max_delta_E < tolerance_rel

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
