#!/usr/bin/env python3

import os
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
import scipy.constants as scc
import yt

import yt ; yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

filename = sys.argv[1]

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)

diag_name = 'diag'
iteration = 0
plotfile = './diags/diag1000200'
field = 'Ey'
ds = yt.load( plotfile )



all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ey = all_data_level_0['boxlib', field].v.squeeze()
Dx = ds.domain_width/ds.domain_dimensions
extent = [ds.domain_left_edge[ds.dimensionality-1], ds.domain_right_edge[ds.dimensionality-1],
          ds.domain_left_edge[0], ds.domain_right_edge[0] ]

plt.figure()
plt.imshow(Ey, extent=extent, aspect='auto')
plt.savefig("Ey_visualization.png")

#Plot A potential vector:
diag_name = 'diag'
iteration = 0
plotfile = './diags/diag1000200'
field = 'A'
ds = yt.load( plotfile )

all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
A_laser_envelope = all_data_level_0['boxlib', field].v.squeeze()
Dx = ds.domain_width/ds.domain_dimensions
extent = [ds.domain_left_edge[ds.dimensionality-1], ds.domain_right_edge[ds.dimensionality-1],
          ds.domain_left_edge[0], ds.domain_right_edge[0] ]
