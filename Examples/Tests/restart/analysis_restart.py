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

species = ['beam', 'plasma_e', 'plasma_p', 'driver', 'driverback']
coords  = ['x', 'y', 'z']

for s in species:
    for c in coords:
        p = 'particle_position_' + c
        print("species : ", s)
        print("position: ", p)
        dr = ad_restart[s, p].to_ndarray()
        db = ad_benchmark[s, p].to_ndarray()
        dr.sort()
        db.sort()
        assert(np.max(abs(dr-db)) < tolerance)

filename = sys.argv[1]
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
