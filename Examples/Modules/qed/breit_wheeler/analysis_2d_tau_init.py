#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import yt
import numpy as np
import scipy.stats as st
import sys
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# This script checks if photons initialized with Breit Wheeler process enabled
# do actually have an exponentially distributed optical depth

# Tolerance
tolerance_rel = 1e-2

def check():
    filename = sys.argv[1]
    data_set = yt.load(filename)

    all_data = data_set.all_data()
    res_tau = all_data["photons", 'particle_optical_depth_BW']

    loc, scale = st.expon.fit(res_tau)

    # loc should be very close to 0, scale should be very close to 1

    error_rel = np.abs(loc - 0)
    print("error_rel for location: " + str(error_rel))
    print("tolerance_rel: " + str(tolerance_rel))
    assert( error_rel < tolerance_rel )

    error_rel = np.abs(scale - 1)
    print("error_rel for scale: " + str(error_rel))
    print("tolerance_rel: " + str(tolerance_rel))
    assert( error_rel < tolerance_rel )

    test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
    checksumAPI.evaluate_checksum(test_name, filename)

def main():
    check()

if __name__ == "__main__":
    main()
