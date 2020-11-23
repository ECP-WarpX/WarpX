#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics.
# The setup is a uniform plasma with electrons, protons and photons.
# Particle energy, field energy and maximum field values will be
# outputed using the reduced diagnostics.
# And they will be compared with the data in the plotfiles.

# Tolerance: 1.0e-8 for particle energy, 1.0e-3 for field energy,
# 1.0e-9 for the maximum electric field and 1.0e-18 for the
# maximum magnetic field.
# The difference of the field energy is relatively large,
# because fields data in plotfiles are cell-centered,
# but fields data in reduced diagnostics are staggered.

# Possible running time: ~ 2 s

import sys
import yt
import numpy as np
import scipy.constants as scc
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

fn = sys.argv[1]

ds = yt.load(fn)
ad = ds.all_data()

# PART1: get results from plotfiles

EPyt = 0.0
# electron
px = ad['electrons','particle_momentum_x'].to_ndarray()
py = ad['electrons','particle_momentum_y'].to_ndarray()
pz = ad['electrons','particle_momentum_z'].to_ndarray()
w  = ad['electrons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_e**2*scc.c**4)-scc.m_e*scc.c**2)*w )
num_electron = w.shape[0]
sum_weight_electron = np.sum(w)

# proton
px = ad['protons','particle_momentum_x'].to_ndarray()
py = ad['protons','particle_momentum_y'].to_ndarray()
pz = ad['protons','particle_momentum_z'].to_ndarray()
w  = ad['protons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_p**2*scc.c**4)-scc.m_p*scc.c**2)*w )
num_proton = w.shape[0]
sum_weight_proton = np.sum(w)

# photon
px = ad['photons','particle_momentum_x'].to_ndarray()
py = ad['photons','particle_momentum_y'].to_ndarray()
pz = ad['photons','particle_momentum_z'].to_ndarray()
w  = ad['photons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt(px**2+py**2+pz**2)*scc.c)*w )
num_photon = w.shape[0]
sum_weight_photon = np.sum(w)

num_total = num_electron + num_proton + num_photon
sum_weight_total = sum_weight_electron + sum_weight_proton + sum_weight_photon

ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex = ad['Ex'].to_ndarray()
Ey = ad['Ey'].to_ndarray()
Ez = ad['Ez'].to_ndarray()
Bx = ad['Bx'].to_ndarray()
By = ad['By'].to_ndarray()
Bz = ad['Bz'].to_ndarray()
Es = np.sum(Ex**2)+np.sum(Ey**2)+np.sum(Ez**2)
Bs = np.sum(Bx**2)+np.sum(By**2)+np.sum(Bz**2)
N  = np.array( ds.domain_width / ds.domain_dimensions )
dV = N[0]*N[1]*N[2]
EFyt = 0.5*Es*scc.epsilon_0*dV + 0.5*Bs/scc.mu_0*dV
max_Ex = np.amax(np.abs(Ex))
max_Ey = np.amax(np.abs(Ey))
max_Ez = np.amax(np.abs(Ez))
max_Bx = np.amax(np.abs(Bx))
max_By = np.amax(np.abs(By))
max_Bz = np.amax(np.abs(Bz))
max_E = np.sqrt(np.amax(Ex**2+Ey**2+Ez**2))
max_B = np.sqrt(np.amax(Bx**2+By**2+Bz**2))

# PART2: get results from reduced diagnostics

EFdata = np.genfromtxt("./diags/reducedfiles/EF.txt")
EPdata = np.genfromtxt("./diags/reducedfiles/EP.txt")
MFdata = np.genfromtxt("./diags/reducedfiles/MF.txt")
NPdata = np.genfromtxt("./diags/reducedfiles/NP.txt")
EF = EFdata[1][2]
EP = EPdata[1][2]
max_Exdata = MFdata[1][2]
max_Eydata = MFdata[1][3]
max_Ezdata = MFdata[1][4]
max_Edata  = MFdata[1][5]
max_Bxdata = MFdata[1][6]
max_Bydata = MFdata[1][7]
max_Bzdata = MFdata[1][8]
max_Bdata  = MFdata[1][9]
num_total_data = NPdata[1][2]
num_electron_data = NPdata[1][3]
num_proton_data = NPdata[1][4]
num_photon_data = NPdata[1][5]
sum_weight_total_data = NPdata[1][6]
sum_weight_electron_data = NPdata[1][7]
sum_weight_proton_data = NPdata[1][8]
sum_weight_photon_data = NPdata[1][9]

# PART3: print and assert

max_diffEmax = max(abs(max_Exdata-max_Ex),abs(max_Eydata-max_Ey),
                   abs(max_Ezdata-max_Ez),abs(max_Edata-max_E))
max_diffBmax = max(abs(max_Bxdata-max_Bx),abs(max_Bydata-max_By),
                   abs(max_Bzdata-max_Bz),abs(max_Bdata-max_B))
max_diff_number = max(abs(num_total_data-num_total),abs(num_electron_data-num_electron),
                   abs(num_proton_data-num_proton),abs(num_photon_data-num_photon))
max_diff_sum_weight = max(abs(sum_weight_total_data-sum_weight_total),
                          abs(sum_weight_electron_data-sum_weight_electron),
                          abs(sum_weight_proton_data-sum_weight_proton),
                          abs(sum_weight_photon_data-sum_weight_photon))

print('difference of field energy:', abs(EFyt-EF))
print('tolerance of field energy:', 1.0e-3)
print('difference of particle energy:', abs(EPyt-EP))
print('tolerance of particle energy:', 1.0e-8)
print('maximum difference of maximum electric field:', max_diffEmax)
print('tolerance of maximum electric field difference:', 1.0e-9)
print('maximum difference of maximum magnetic field:', max_diffBmax)
print('tolerance of maximum magnetic field difference:', 1.0e-18)
print('maximum difference of particle weight sum:', max_diff_sum_weight)
print('tolerance of particle weight sum:', 0.5)

assert(abs(EFyt-EF) < 1.0e-3)
assert(abs(EPyt-EP) < 1.0e-8)
assert(max_diffEmax < 1.0e-9)
assert(max_diffBmax < 1.0e-18)
assert(max_diff_number < 0.5)
assert(max_diff_sum_weight < 0.5)

test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, fn)
