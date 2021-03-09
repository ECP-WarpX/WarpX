#! /usr/bin/env python

# Copyright 2021 Roelof Groenewald

# This script tests the time dependent Dirichlet boundary
# conditions in a 2D electrostatic simulation. An empty
# domain of 64 x 8 cells is simulated with periodic boundary
# conditions in the x directions and Dirichlet boundary
# conditions in the y direction with specified potentials
# of 0 V on the lo side and a sine wave on the hi side.
# One period of the sine wave is simulated and the potentials
# at the boundaries compared to expectation.

# Possible running time: ~ 19 s

import h5py
import numpy as np

import glob

dt = 7.5e-10

files = sorted(glob.glob('diags/openpmd_*.h5'))[1:]

times = np.zeros(len(files))
potentials = np.zeros(len(files))

for ii, file in enumerate(files):
    with h5py.File(file, 'r') as f:
        step = str(int(file.split('_')[-1][:-3]))
        times[ii] = (float(step) - 1.0) * dt
        potentials[ii] = np.mean(
            np.array(f['data'][step]['fields']['phi']), axis=0
        )[-1]

expected_potentials = 450.0 * np.sin(2.0 * np.pi * 13.56e6 * times)
assert np.allclose(potentials, expected_potentials, rtol=0.01)
