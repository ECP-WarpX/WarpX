#! /usr/bin/env python

import sys
import yt
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]
ds_restart = yt.load( filename )

benchmark = 'orig_' + filename
ds_benchmark = yt.load( benchmark )

# Check particle output
ad_restart = ds_restart.all_data()
ad_benchmark = ds_benchmark.all_data()

all_species = ['beam', 'plasma_e', 'plasma_p', 'driver', 'driverback']
all_attrs = ['particle_position_x', 'particle_position_y', 'particle_position_z',
             'particle_momentum_x', 'particle_momentum_y', 'particle_momentum_z',
             'particle_weight']

tolerance = 1e-15
print('\ntolerance (particle output) = {:g}'.format(tolerance))
for species in all_species:
    print()
    print('species: {:s}'.format(species))
    for attr in all_attrs:
        dr = ad_restart[species, attr].to_ndarray()
        db = ad_benchmark[species, attr].to_ndarray()
        dr.sort()
        db.sort()
        error = np.max(np.abs(dr - db))
        print('attribute: {:s}; error = {:g}'.format(attr, error))
        assert(error < tolerance)
print()

# Check field output
ad_restart = ds_restart.covering_grid(level = 0,
        left_edge = ds_restart.domain_left_edge, dims = ds_restart.domain_dimensions)
ad_benchmark = ds_benchmark.covering_grid(level = 0,
        left_edge = ds_benchmark.domain_left_edge, dims = ds_benchmark.domain_dimensions)

all_fields = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'jx', 'jy', 'jz']

tolerance = 1e-5
print('tolerance (field output) = {:g}'.format(tolerance))
print()
for field in all_fields:
    dr = ad_restart['boxlib', field].squeeze().v
    db = ad_benchmark['boxlib', field].squeeze().v
    error = np.amax(np.abs(dr - db))
    if (np.amax(np.abs(db)) != 0.):
        error /= np.amax(np.abs(db))
    print('field: {:s}; relative error = {:g}'.format(field, error))
    assert(error < tolerance)
print()

# Check-sum analysis
filename = sys.argv[1]
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
