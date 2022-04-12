#!/usr/bin/env python3

import openpmd_api as io

series = io.Series("LaserAccelerationRZ_opmd_plt/openpmd_%T.h5", io.Access.read_only)

assert len(series.iterations) == 3, 'improper number of iterations stored'

ii = series.iterations[20]

assert len(ii.meshes) == 7, 'improper number of meshes'

# select j_t
jt = ii.meshes['j']['t']

# this is in C (Python) order; r is the fastest varying index
(Nm, Nz, Nr) = jt.shape

assert Nm == 3, 'Wrong number of angular modes stored or possible incorrect ordering when flushed'
assert Nr == 64, 'Wrong number of radial points stored or possible incorrect ordering when flushed'
assert Nz == 512, 'Wrong number of z points stored or possible incorrect ordering when flushed'

assert ii.meshes['part_per_grid'][io.Mesh_Record_Component.SCALAR].shape == [512,64], 'problem with part_per_grid'
assert ii.meshes['rho_electrons'][io.Mesh_Record_Component.SCALAR].shape == [3, 512, 64], 'problem with rho_electrons'
