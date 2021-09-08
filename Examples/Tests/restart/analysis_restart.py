#! /usr/bin/env python

import sys
import re
import yt
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = sys.float_info.epsilon
print('tolerance = ', tolerance)

filename = sys.argv[1]
psatd = True if re.search('psatd', filename) else False

ds  = yt.load( filename )
ad  = ds.all_data()
xb  = ad['beam',     'particle_position_x'].to_ndarray()
xe  = ad['plasma_e', 'particle_position_x'].to_ndarray()
zb  = ad['beam',     'particle_position_z'].to_ndarray()
ze  = ad['plasma_e', 'particle_position_z'].to_ndarray()

filename = 'orig_restart_psatd_plt00010' if (psatd) else 'orig_restart_plt00010'
ds  = yt.load( filename )
ad  = ds.all_data()
xb0 = ad['beam',     'particle_position_x'].to_ndarray()
xe0 = ad['plasma_e', 'particle_position_x'].to_ndarray()
zb0 = ad['beam',     'particle_position_z'].to_ndarray()
ze0 = ad['plasma_e', 'particle_position_z'].to_ndarray()

xb.sort()
xb0.sort()
assert(np.max(abs(xb-xb0))<tolerance)

xe.sort()
xe0.sort()
assert(np.max(abs(xe-xe0))<tolerance)

zb.sort()
zb0.sort()
assert(np.max(abs(zb-zb0))<tolerance)

ze.sort()
ze0.sort()
assert(np.max(abs(ze-ze0))<tolerance)

filename = sys.argv[1]
test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
