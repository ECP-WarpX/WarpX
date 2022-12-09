#!/usr/bin/env python3

# Run the default regression test for the PICMI version of the EB test
# using the same reference file as for the non-PICMI test since the two
# tests are otherwise the same.

import sys

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI

my_check = checksumAPI.evaluate_checksum(
    'Python_magnetostatic_eb_3d', 'Python_MagnetostaticEB_plt000000', rtol=1e-3
)
