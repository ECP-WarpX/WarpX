#!/usr/bin/env python3

# Run the default regression test for the PICMI version of the EB test
# using the same reference file as for the non-PICMI test since the two
# tests are otherwise the same.

import sys

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI

# The first step is the same as in inputs_3d
my_check = checksumAPI.evaluate_checksum(
    'ElectrostaticSphereEB', 'Python_ElectrostaticSphereEB_plt000001',
    do_particles=False, atol=1e-12
)

# The second step has the EB potential modified via Python
my_check = checksumAPI.evaluate_checksum(
    'Python_ElectrostaticSphereEB', 'Python_ElectrostaticSphereEB_plt000002',
    do_particles=False, atol=1e-12
)
