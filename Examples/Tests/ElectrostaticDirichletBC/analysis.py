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

import numpy as np

import yt
import glob

files = sorted(glob.glob('dirichletbc_plt*'))[1:]

times = np.zeros(len(files))
potentials = np.zeros(len(files))

for ii, file in enumerate(files):
    ds = yt.load( file )
    times[ii] = (
        ds.current_time.item() - float(ds.parameters.get('warpx.const_dt'))
    )
    data = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    potentials[ii] = np.mean(data['phi'].to_ndarray()[-1][:])

expected_potentials = 450.0 * np.sin(2.0 * np.pi * 13.56e6 * times)
assert np.allclose(potentials, expected_potentials, rtol=0.01)
