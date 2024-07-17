#!/usr/bin/env python

# Copyright 2024 Olga Shapoval, Edoardo Zoni
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the embedded boundary in RZ.
# A cylindrical surface (r=0.1) has a fixed potential 1 V.
# The outer surface has 0 V fixed.
# Thus the analytical solution has the form:
# phi(r) = A+B*log(r), Er(r) = -B/r.

# Possible errors:
# tolerance: 0.004
# Possible running time: < 1 s

import os
import sys

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 0.0041

fn = sys.argv[1]

def find_first_non_zero_from_bottom_left(matrix):
    nj = len(matrix[0,:])
    ni = len(matrix[:,0])
    for i in range(ni):
        for j in range(nj):
            if (matrix[i][j] != 0) and (matrix[i][j] != np.nan):
                return (i, j)
    return i, j

def find_first_non_zero_from_upper_right(matrix):
    nj = len(matrix[0,:])
    ni = len(matrix[:,0])

    for i in range(ni-1, -1, -1):
        for j in range(nj-1, -1, -1):
            if (matrix[i][j] != 0) and (matrix[i][j] != np.nan):
                return (i, j)
    return i,j

def get_fields(ts, level):
    er = []
    phi = []
    info = []
    if level == 0:
        er, info = ts.get_field('E', 'r', iteration=0)
        phi, info = ts.get_field('phi', iteration=0)
    else:
        er, info = ts.get_field(f'E_lvl{level}', 'r', iteration=0)
        phi, info = ts.get_field(f'phi_lvl{level}', iteration=0)

    return er, phi, info

def get_error_per_lev(ts,level):
    er, phi, info = get_fields(ts, level)

    nj = info.r.shape[0] // 2
    dr = info.dr

    er_patch = er[:,nj:]
    phi_patch = phi[:,nj:]
    r1 = info.r[nj:]
    patch_left_lower_i, patch_left_lower_j = find_first_non_zero_from_bottom_left(er_patch)
    patch_right_upper_i, patch_right_upper_j = find_first_non_zero_from_upper_right(er_patch)

    # phi and Er field on the MR patch
    phi_sim = phi_patch[patch_left_lower_i:patch_right_upper_i+1, patch_left_lower_j:patch_right_upper_j+1]
    er_sim =  er_patch[patch_left_lower_i:patch_right_upper_i+1, patch_left_lower_j:patch_right_upper_j+1]
    r =  r1[patch_left_lower_j:patch_right_upper_j+1]

    B = 1.0/np.log(0.1/0.5)
    A = -B*np.log(0.5)

    errmax_phi = 0.0
    errmax_Er = 0.0

    for j in range(len(r)):
        # outside EB and last cutcell
        if r[j] > 0.1 + dr:
            phi_theory = A + B*np.log(r[j])
            Er_theory = -B / float(r[j])
            err = abs( phi_theory - phi_sim [:,j] ).max() / phi_theory
            if err > errmax_phi:
                errmax_phi = err
            err = abs( Er_theory - er_sim [:,j] ).max() / Er_theory
            if err > errmax_Er and j < len(r) - 1:
                errmax_Er = err

    print('max error of phi = ', errmax_phi)
    print('max error of Er = ', errmax_Er)
    print('tolerance = ', tolerance)
    assert(errmax_phi < tolerance)
    assert(errmax_Er < tolerance)

ts = OpenPMDTimeSeries(fn)
level_fields = [field for field in ts.avail_fields if 'lvl' in field]
nlevels = 0 if level_fields == [] else int(level_fields[-1][-1])
for level in range(nlevels+1):
    print('At level =', level, ':')
    get_error_per_lev(ts,level)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn, output_format="openpmd", do_particles=False)
