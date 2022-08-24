#!/usr/bin/env python3

# Copyright 2019-2020 Luca Fedeli, Maxence Thevenet, Revathi Jambunathan
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


'''
Analysis script of a WarpX simulation of rigid injection in a boosted frame.

A Gaussian electron beam starts from -5 microns, propagates rigidly up to
20 microns after which it expands due to emittance only (the focal position is
20 microns). The beam width is measured after ~50 microns, and compared with
the theory (with a 5% error allowed).

The simulation runs in a boosted frame, and the analysis is done in the lab
frame, i.e., on the back-transformed diagnostics.
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

# Read data from legacy back-transformed diagnostics
snapshot = './lab_frame_data/snapshots/snapshot00001'
header   = './lab_frame_data/snapshots/Header'
allrd, info = read_raw_data.read_lab_snapshot(snapshot, header)
z_btd_legacy = read_raw_data.get_particle_field(snapshot, 'beam', 'z')
x_btd_legacy = read_raw_data.get_particle_field(snapshot, 'beam', 'x')
z_btd_legacy_mean = np.mean(z_btd_legacy)
x_btd_legacy_std = np.std(x_btd_legacy)

# Read data from new back-transformed diagnostics (plotfile)
ds_plotfile = yt.load(filename)
z_btd_plotfile = ds_plotfile.all_data()['beam', 'particle_position_y'].v
x_btd_plotfile = ds_plotfile.all_data()['beam', 'particle_position_x'].v
z_btd_plotfile_mean = np.mean(z_btd_plotfile)
x_btd_plotfile_std = np.std(x_btd_plotfile)

# Read data from new back-transformed diagnostics (openPMD)
series = io.Series("./diags/diag2/openpmd_%T.h5", io.Access.read_only)
ds_openpmd = series.iterations[1]
z_btd_openpmd = ds_openpmd.particles['beam']['position']['z'][:]
x_btd_openpmd = ds_openpmd.particles['beam']['position']['x'][:]
series.flush()
z_btd_openpmd_mean = np.mean(z_btd_openpmd)
x_btd_openpmd_std = np.std(x_btd_openpmd)

# Sort and compare arrays to check consistency between legacy BTD and new BTD
rtol = 1e-16
atol = 1e-16
assert(np.allclose(np.sort(z_btd_legacy), np.sort(z_btd_plotfile), rtol=rtol, atol=atol))
assert(np.allclose(np.sort(x_btd_legacy), np.sort(x_btd_plotfile), rtol=rtol, atol=atol))
assert(np.allclose(np.sort(z_btd_legacy), np.sort(z_btd_openpmd), rtol=rtol, atol=atol))
assert(np.allclose(np.sort(x_btd_legacy), np.sort(x_btd_openpmd), rtol=rtol, atol=atol))

# initial parameters
z0 = 20.e-6
x0 = 1.e-6
theta0 = np.arcsin(0.1)

# Theoretical beam width after propagation if rigid ON
z = z_btd_legacy_mean
x = x_btd_legacy_std
xth = np.sqrt( x0**2 + (z-z0)**2*theta0**2 )
error_rel = np.abs((x-xth)/xth)
tolerance_rel = 1e-2

# Print error and assert small error
print("Beam position: " + str(z))
print("Beam width   : " + str(x))

print("error_rel    : " + str(error_rel))
print("tolerance_rel: " + str(tolerance_rel))

assert( error_rel < tolerance_rel )

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
