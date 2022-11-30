#!/usr/bin/env python3

# Copyright 2021 Modern Electron

# This script checks that the PICMI_inputs_2d.py run more-or-less matches the
# results from the non-PICMI run. The PICMI run is using an external Poisson
# solver that directly solves the Poisson equation using matrix inversion
# rather than the iterative approach from the MLMG solver.

import sys

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI

my_check = checksumAPI.evaluate_checksum(
    'background_mcc', 'Python_background_mcc_plt000050',
    do_particles=True, rtol=3.5e-3
)
