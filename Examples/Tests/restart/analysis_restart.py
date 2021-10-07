#! /usr/bin/env python

import sys
import yt
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Load output data generated after restart
filename = sys.argv[1]
ds_restart = yt.load(filename)
ds_restart.force_periodicity()
ad_restart = ds_restart.covering_grid(level = 0,
                                      left_edge = ds_restart.domain_left_edge,
                                      dims = ds_restart.domain_dimensions)

# Load output data generated from initial run
benchmark = 'orig_' + filename
ds_benchmark = yt.load(benchmark)
ds_benchmark.force_periodicity()
ad_benchmark = ds_benchmark.covering_grid(level = 0,
                                          left_edge = ds_benchmark.domain_left_edge,
                                          dims = ds_benchmark.domain_dimensions)

# Loop over all fields (all particle species, all particle attributes, all grid fields)
# and compare output data generated from initial run with output data generated after restart
tolerance = 1e-12
print('\ntolerance = {:g}'.format(tolerance))
print()
for field in ds_benchmark.field_list:
    dr = ad_restart[field].squeeze().v
    db = ad_benchmark[field].squeeze().v
    error = np.amax(np.abs(dr - db))
    if (np.amax(np.abs(db)) != 0.):
        error /= np.amax(np.abs(db))
    print('field: {}; error = {:g}'.format(field, error))
    assert(error < tolerance)
print()

# Check-sum analysis
filename = sys.argv[1]
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
