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

## This is a generic function to test a particle filter. We reproduce the filter in python and
## verify that the results are the same as with the filtered diagnostic.
def check_particle_filter(fn, filtered_fn, filter_expression):
    ds  = yt.load( fn )
    ds_filtered  = yt.load( filtered_fn )
    ad  = ds.all_data()
    ad_filtered  = ds_filtered.all_data()

    ## Load arrays from the unfiltered diagnostic
    ids = ad['electron', 'particle_id'].to_ndarray()
    cpus = ad['electron', 'particle_cpu'].to_ndarray()
    x = ad['electron', 'particle_position_x'].to_ndarray()
    px  = ad['electron', 'particle_momentum_x'].to_ndarray()
    z = ad['electron', 'particle_position_y'].to_ndarray()
    pz  = ad['electron', 'particle_momentum_z'].to_ndarray()
    py  = ad['electron', 'particle_momentum_y'].to_ndarray()
    w  = ad['electron', 'particle_weight'].to_ndarray()

    ## Load arrays from the filtered diagnostic
    ids_filtered_warpx = ad_filtered['particle_id'].to_ndarray()
    cpus_filtered_warpx = ad_filtered['particle_cpu'].to_ndarray()
    x_filtered_warpx = ad_filtered['particle_position_x'].to_ndarray()
    px_filtered_warpx  = ad_filtered['particle_momentum_x'].to_ndarray()
    z_filtered_warpx = ad_filtered['particle_position_y'].to_ndarray()
    pz_filtered_warpx  = ad_filtered['particle_momentum_z'].to_ndarray()
    py_filtered_warpx  = ad_filtered['particle_momentum_y'].to_ndarray()
    w_filtered_warpx  = ad_filtered['particle_weight'].to_ndarray()

    ## Reproduce the filter in python: this returns the indices of the filtered particles in the
    ## unfiltered arrays.
    ind_filtered_python, = numpy.where(eval(filter_expression))

    ## Sort the indices of the filtered arrays by particle id.
    sorted_ind_filtered_python = ind_filtered_python[numpy.argsort(ids[ind_filtered_python])]
    sorted_ind_filtered_warpx = numpy.argsort(ids_filtered_warpx)

    ## Check that the sorted ids are exactly the same with the warpx filter and the filter
    ## reproduced in python
    assert(numpy.array_equal(ids[sorted_ind_filtered_python],
                             ids_filtered_warpx[sorted_ind_filtered_warpx]))
    assert(numpy.array_equal(cpus[sorted_ind_filtered_python],
                             cpus_filtered_warpx[sorted_ind_filtered_warpx]))

    ## Finally, we check that the sum of the particles quantities are the same to machine precision
    tolerance_checksum = 1.e-12
    check_array_sum(x[sorted_ind_filtered_python],
                    x_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(px[sorted_ind_filtered_python],
                    px_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(z[sorted_ind_filtered_python],
                    z_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(pz[sorted_ind_filtered_python],
                    pz_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(py[sorted_ind_filtered_python],
                    py_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)
    check_array_sum(w[sorted_ind_filtered_python],
                    w_filtered_warpx[sorted_ind_filtered_warpx], tolerance_checksum)

## This function checks that the absolute sums of two arrays are the same to a required precision
def check_array_sum(array1, array2, tolerance_checksum):
    sum1 = numpy.sum(numpy.abs(array1))
    sum2 = numpy.sum(numpy.abs(array2))
    assert(abs(sum2-sum1)/sum1 < tolerance_checksum)

## This function is specifically used to test the random filter. First, we check that the number of
## dumped particles is as expected. Next, we call the generic check_particle_filter function.
def check_random_filter(fn, filtered_fn, random_fraction):
    ds  = yt.load( fn )
    ds_filtered  = yt.load( filtered_fn )
    ad  = ds.all_data()
    ad_filtered  = ds_filtered.all_data()

    ## Check that the number of particles is as expected
    numparts = ad['electron', 'particle_id'].to_ndarray().shape[0]
    numparts_filtered = ad_filtered['particle_id'].to_ndarray().shape[0]
    expected_numparts_filtered = random_fraction*numparts
    # 5 sigma test that has an intrinsic probability to fail of 1 over ~2 millions
    std_numparts_filtered = numpy.sqrt(expected_numparts_filtered)
    error = abs(numparts_filtered-expected_numparts_filtered)
    print("Random filter: difference between expected and actual number of dumped particles: " \
          + str(error))
    print("tolerance: " + str(5*std_numparts_filtered))
    assert(error<5*std_numparts_filtered)

    ## Dirty trick to find particles with the same ID + same CPU (does not work with more than 10
    ## MPI ranks)
    random_filter_expression = 'numpy.isin(ids + 0.1*cpus,' \
                                          'ids_filtered_warpx + 0.1*cpus_filtered_warpx)'
    check_particle_filter(fn, filtered_fn, random_filter_expression)

parser_filter_fn = "diags/diag_parser_filter00150"
parser_filter_expression = "(x>200) * (z<200) * (px-3*pz>0)"
check_particle_filter(last_fn, parser_filter_fn, parser_filter_expression)

uniform_filter_fn = "diags/diag_uniform_filter00150"
uniform_filter_expression = "ids%6 == 0"
check_particle_filter(last_fn, uniform_filter_fn, uniform_filter_expression)

random_filter_fn = "diags/diag_random_filter00150"
random_fraction = 0.77
check_random_filter(last_fn, random_filter_fn, random_fraction)

test_name = last_fn[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, do_particles=False)
