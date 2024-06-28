#!/usr/bin/env python3

# Copyright 2022 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This test tests the external field loading feature.
# A magnetic mirror field is loaded, and a single particle
# in the mirror will be reflected by the magnetic mirror effect.
# At the end of the simulation, the position of the particle
# is compared with known correct results.

# Possible errors: 6.235230443866285e-9
# tolerance: 1.0e-8
# Possible running time: 0.327827743 s

import os
import sys

import numpy as np
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

tolerance = 1.0e-12

filename = sys.argv[1]

ds = yt.load( filename )
grid0 = ds.index.grids[0]


r_min = grid0.LeftEdge.v[0]
r_max = grid0.RightEdge.v[0]
nr = grid0.shape[0]
dr = grid0.dds.v[0]
rv = np.linspace(r_min, r_max, num=int(nr+1))
rbar = 0.5*(rv[1:] + rv[:-1])
r_low = rv[:-1]
r_hi = rv[1:]

z_min = grid0.LeftEdge.v[1]
z_max = grid0.RightEdge.v[1]
nz = grid0.shape[1]
dz = grid0.dds.v[1]
zv = np.linspace(z_min, z_max, num=int(nz+1))
zbar = 0.5*(zv[1:] + zv[:-1])

RL, ZL = np.meshgrid(r_low, zbar, indexing='ij')
RH, ZH = np.meshgrid(r_hi, zbar, indexing='ij')
RB, ZB = np.meshgrid(rbar, zbar, indexing='ij')

dBrdr = (RH*grid0['raw', 'Br_aux'].v[:,:,1] - RL*grid0['raw', 'Br_aux'].v[:,:,0])/(dr*RB)
dBzdz = (grid0['raw', 'Bz_aux'].v[:,:,1] - grid0['raw', 'Bz_aux'].v[:,:,0])/dz

divB = dBrdr + dBzdz

import matplotlib.pyplot as plt

plt.imshow(np.log10(np.abs(divB[:,24,:])))
plt.title('log10(|div(B)|)')
plt.colorbar()
plt.savefig('divb.png')

error = np.sqrt((divB[2:-2,2:-2]**2).sum())

print('error = ', error)
print('tolerance = ', tolerance)
assert(error < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
