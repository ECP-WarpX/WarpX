#!/usr/bin/env python3

# Copyright 2019 Maxence Thevenet, Revathi Jambunathan, Arianna Formenti
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import os
import sys

import numpy as np
import openpmd_api as io
from scipy.constants import c
from scipy.constants import e as q_e
from scipy.constants import eV, m_e, micro, nano

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

GeV=1e9*eV
energy = 125.*GeV
gamma = energy/(m_e*c**2)
sigmax = 516.0*nano
sigmay = 7.7*nano
sigmaz = 300.*micro
nz = 256
Lz = 20*sigmaz
gridz = np.linspace(-0.5*Lz, 0.5*Lz, nz)
tol = gridz[1] - gridz[0]
emitx = 5*micro
emity = 35*nano
focal_distance = 4*sigmaz


def s(z, sigma0, emit):
    return np.sqrt(sigma0**2 + emit**2 * (z - focal_distance)**2 / sigma0**2)


filename = sys.argv[1]
series = io.Series("./diags/full/openpmd_%T.bp", io.Access.read_only)


it = series.iterations[0]
ps = it.particles["beam1"]
x = ps["position"]["x"].load_chunk()
y = ps["position"]["y"].load_chunk()
z = ps["position"]["z"].load_chunk()
w = ps["weighting"][io.Mesh_Record_Component.SCALAR].load_chunk()
series.flush()

it.close()
del series

imin = np.argmin(np.sqrt((gridz+focal_distance)**2))
imax = np.argmin(np.sqrt((gridz-focal_distance)**2))


sx, sy = [], []

subgrid = gridz[imin:imax]
for d in subgrid:
    i = np.sqrt((z - d)**2) < tol
    if (np.sum(i)!=0):
        mux = np.average(x[i], weights=w[i])
        muy = np.average(y[i], weights=w[i])
        sx.append(np.sqrt(np.average((x[i]-mux)**2, weights=w[i])))
        sy.append(np.sqrt(np.average((y[i]-muy)**2, weights=w[i])))

sx_theory = s(subgrid, sigmax, emitx/gamma)
sy_theory = s(subgrid, sigmay, emity/gamma)

print('XXXXXXXXXXXXXXXXXXXXXXXX')
assert(np.allclose(sx, sx_theory, rtol=0.12, atol=0))
assert(np.allclose(sy, sy_theory, rtol=0.12, atol=0))


test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
