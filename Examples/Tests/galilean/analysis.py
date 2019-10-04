#! /usr/bin/env python

import sys
import yt ; yt.funcs.mylog.setLevel(0)
import numpy as np
import scipy.constants as scc

filename = sys.argv[1]


ds = yt.load( filename )
all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = all_data_level_0['boxlib', 'Ex'].v.squeeze()
Ey = all_data_level_0['boxlib', 'Ey'].v.squeeze()
Ez = all_data_level_0['boxlib', 'Ez'].v.squeeze()
energyE_gal_psatd = np.sum(scc.epsilon_0/2*(Ex**2+Ey**2+Ez**2))

energyE_psatd = 2.99915183857e-20 #E field energy calculated with PSATD (v_galilean = (0,0,0))
assert( energyE_gal_psatd < 1e-18 )
