#! /usr/bin/env python

import sys
import yt
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = sys.float_info.epsilon
print('tolerance = ', tolerance)

filename = sys.argv[1]
ds_restart = yt.load( filename )
ad_restart = ds_restart.all_data()

benchmark = 'orig_' + filename
ds_benchmark = yt.load( benchmark )
ad_benchmark = ds_benchmark.all_data()

all_species = ['beam', 'plasma_e', 'plasma_p', 'driver', 'driverback']
all_attrs = ['particle_position_x', 'particle_position_y', 'particle_position_z',
             'particle_momentum_x', 'particle_momentum_y', 'particle_momentum_z']

for species in all_species:
    for attr in all_attrs:
        dr = ad_restart[species, attr].to_ndarray()
        db = ad_benchmark[species, attr].to_ndarray()
        dr.sort()
        db.sort()
        error = np.max(np.abs(dr - db))
        assert(error < tolerance)

filename = sys.argv[1]
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
