#!/usr/bin/env python3

# Copyright 2021 Modern Electron

import numpy as np

density_data = np.load( 'avg_ion_density.npy' )
assert np.isclose(np.mean(density_data), 2.53530e14)
