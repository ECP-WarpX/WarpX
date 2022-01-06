#!/usr/bin/env python3

# Run the default regression test for the PICMI version of the EB test
# using the same reference file as for the non-PICMI test since the two
# tests are otherwise the same.

import sys

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI

my_check = checksumAPI.evaluate_checksum(
    'ElectrostaticSphereEB', 'Python_ElectrostaticSphereEB_plt00001',
    do_particles=False, atol=1e-12
)
