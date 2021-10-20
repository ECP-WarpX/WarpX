#! /usr/bin/env python

# Copyright 2021 David Grote
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""
This script tests the absorption of fields in the PML in RZ geometry.

The input file inputs_particle_rz is used: it features an electron
moving radially that launches a pulse. This scripts runs until
most of the pulse escapes the radial boundary. If the PML fails,
the pulse will remain with in the domain.
"""
import sys
import numpy as np
import yt
yt.funcs.mylog.setLevel(0)
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Open plotfile specified in command line
filename = sys.argv[1]
ds = yt.load( filename )

# Check that the field is low enough
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex_array = ad0['boxlib', 'Ex'].to_ndarray()
Ez_array = ad0['boxlib', 'Ez'].to_ndarray()
max_Ex = np.abs(Ex_array).max()
max_Ez = np.abs(Ez_array).max()
print( f'max Ex = {max_Ex}' )
print( f'max Ez = {max_Ez}' )
max_Efield = max(max_Ex, max_Ez)

tolerance_abs = 2.
print('tolerance_abs: ' + str(tolerance_abs))
assert max_Efield < tolerance_abs

test_name = filename[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename)
