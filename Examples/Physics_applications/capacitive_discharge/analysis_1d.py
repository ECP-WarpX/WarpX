#!/usr/bin/env python3

# Copyright 2021 Modern Electron

import numpy as np

density_data = np.load( 'avg_ion_density.npy' )
assert np.isclose(np.mean(density_data), 258267505634481.94)
