#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the collision module
# using electron-ion temperature relaxation in 2D.
# Initially, electrons and ions are both in equilibrium
# (gaussian) distributions, but have different temperatures.
# Relaxation occurs to bring the two temperatures to be
# a final same temperature through collisions.
# The code was tested to be valid, more detailed results
# were used to obtain an exponential fit with
# coefficients a and b.
# This automated test compares the results with the fit.
# Unrelated to the collision module, we also test the plotfile particle filter function in this
# analysis script.

# Possible errors:
# tolerance: 0.001
# Possible running time: ~ 30.0 s

import sys
import yt
import math
import numpy
from glob import glob
import post_processing_utils
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 0.001

ng = 64
ne = ng * 200
ni = ng * 200
np = ne + ni

c  = 299792458.0
me = 9.10938356e-31
mi = me * 5.0

## In the first part of the test we verify that the output data is consistent with the exponential
## fit.

# exponential fit coefficients
a =  0.041817463099883
b = -0.083851393560288

last_fn = sys.argv[1]
if (last_fn[-1] == "/"): last_fn = last_fn[:-1]
fn_list = glob(last_fn[:-5] + "?????")

error = 0.0
nt = 0
for fn in fn_list:
    # load file
    ds  = yt.load( fn )
    ad  = ds.all_data()
    px  = ad['particle_momentum_x'].to_ndarray()
    # get time index j
    j = int(fn[-5:])
    # compute error
    vxe = numpy.mean(px[ 0:ne])/me/c
    vxi = numpy.mean(px[ne:np])/mi/c
    vxd = vxe - vxi
    fit = a*math.exp(b*j)
    error = error + abs(fit-vxd)
    nt = nt + 1

error = error / nt

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)


## In the second past of the test, we verify that the diagnostic particle filter function works as
## expected. For this, we only use the last simulation timestep.

dim = "2d"
species_name = "electron"

parser_filter_fn = "diags/diag_parser_filter00150"
parser_filter_expression = "(x>200) * (z<200) * (px-3*pz>0)"
post_processing_utils.check_particle_filter(last_fn, parser_filter_fn, parser_filter_expression,
                                            dim, species_name)

uniform_filter_fn = "diags/diag_uniform_filter00150"
uniform_filter_expression = "ids%6 == 0"
post_processing_utils.check_particle_filter(last_fn, uniform_filter_fn, uniform_filter_expression,
                                            dim, species_name)

random_filter_fn = "diags/diag_random_filter00150"
random_fraction = 0.77
post_processing_utils.check_random_filter(last_fn, random_filter_fn, random_fraction,
                                          dim, species_name)

test_name = last_fn[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, do_particles=False)
