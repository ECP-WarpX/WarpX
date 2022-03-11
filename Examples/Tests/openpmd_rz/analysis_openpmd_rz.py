#!/usr/bin/env python3

import openpmd_api as io

series = io.Series("diags/diag1/openpmd_%T.h5", io.Access.read_only)

assert len(series.iterations) == 3, 'improper number of iterations stored'

ii = series.iterations[2]

assert len(ii.meshes) == 9, 'improper number of meshes'

jt = ii.meshes['j']['t']
(Nm, Nr, Nz) = jt.shape

assert Nm == 3, 'Wrong number of angular modes stored or possible incorrect ordering when flushed'
assert Nr == 64, 'Wrong number of radial points stored or possible incorrect ordering when flushed'
assert Nz == 512, 'Wrong number of z points stored or possible incorrect ordering when flushed'

assert ii.meshes['part_per_grid'][io.Mesh_Record_Component.SCALAR].shape == [64,512], 'problem with part_per_grid'
assert ii.meshes['divE'][io.Mesh_Record_Component.SCALAR].shape == [3, 64, 512], 'problem with rho_electrons'
assert ii.meshes['divB'][io.Mesh_Record_Component.SCALAR].shape == [3, 64, 512], 'problem with rho_electrons'
assert ii.meshes['rho_electrons'][io.Mesh_Record_Component.SCALAR].shape == [3, 64, 512], 'problem with rho_electrons'
