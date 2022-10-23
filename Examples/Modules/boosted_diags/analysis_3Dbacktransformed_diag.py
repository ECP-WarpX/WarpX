#!/usr/bin/env python3

# Copyright 2019 Maxence Thevenet, Revathi Jambunathan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


'''
Analysis script of a WarpX simulation in a boosted frame.

The simulation runs in a boosted frame, and the analysis is done in the lab
frame, i.e., on the back-transformed diagnostics for the full 3D simulation and
an x-z slice at y=y_center. The field-data, Ez, along z, at (x_center,y_center,:) is compared
between the full back-transformed diagnostic and the reduced diagnostic (i.e., x-z slice) .
'''

import os
import sys

import numpy as np
import openpmd_api as io
import read_raw_data
import yt

yt.funcs.mylog.setLevel(0)

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

# Tolerances to check consistency between legacy BTD and new BTD
rtol = 1e-16
atol = 1e-16

# Read data from legacy back-transformed diagnostics (entire domain)
snapshot = './lab_frame_data/snapshots/snapshot00003'
header   = './lab_frame_data/snapshots/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
Ez_legacy = allrd['Ez']
print(f'Ez_legacy.shape = {Ez_legacy.shape}')
Ez_legacy_1D = np.squeeze(Ez_legacy[Ez_legacy.shape[0]//2,Ez_legacy.shape[1]//2,:])

# Read data from reduced back-transformed diagnostics (slice)
snapshot_slice = './lab_frame_data/slices/slice00003'
header_slice   = './lab_frame_data/slices/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot_slice, header_slice)
Ez_legacy_slice = allrd['Ez']
print(f'Ez_legacy_slice.shape = {Ez_legacy_slice.shape}')
Ez_legacy_slice_1D = np.squeeze(Ez_legacy_slice[Ez_legacy_slice.shape[0]//2,1,:])

# Read data from new back-transformed diagnostics (plotfile)
ds_plotfile = yt.load(filename)
data = ds_plotfile.covering_grid(
    level=0,
    left_edge=ds_plotfile.domain_left_edge,
    dims=ds_plotfile.domain_dimensions)
Ez_plotfile = data[('mesh', 'Ez')].to_ndarray()

# Read data from new back-transformed diagnostics (openPMD)
series = io.Series("./diags/diag2/openpmd_%T.h5", io.Access.read_only)
ds_openpmd = series.iterations[3]
Ez_openpmd = ds_openpmd.meshes['E']['z'].load_chunk()
Ez_openpmd = Ez_openpmd.transpose()
series.flush()

# Compare arrays to check consistency between new BTD formats (plotfile and openPMD)
assert(np.allclose(Ez_plotfile, Ez_openpmd, rtol=rtol, atol=atol))

# Check slicing
err = np.max(np.abs(Ez_legacy_slice_1D-Ez_legacy_1D)) / np.max(np.abs(Ez_legacy_1D))
tol = 1e-16
print(f'error = {err}')
print(f'tolerance = {tol}')
assert(err < tol)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
