#! /usr/bin/env python

# Copyright 2019-2020 Yinjian Zhao
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics.
# The setup is a uniform plasma with electrons, protons and photons.
# Particle energy and field energy will be outputed
# using the reduced diagnostics.
# And they will be compared with the data in the plotfiles.

# Tolerance: 1.0e-8 for particle energy, 1.0e-3 for field energy.
# The difference of the field energy is relatively large,
# because fields data in plotfiles are cell-centered,
# but fields data in reduced diagnostics are staggered.

# Possible running time: ~ 2 s

import sys
import yt
import numpy as np
import scipy.constants as scc

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
# proton
px = ad['protons','particle_momentum_x'].to_ndarray()
py = ad['protons','particle_momentum_y'].to_ndarray()
pz = ad['protons','particle_momentum_z'].to_ndarray()
w  = ad['protons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt((px**2+py**2+pz**2)*scc.c**2+scc.m_p**2*scc.c**4)-scc.m_p*scc.c**2)*w )

# photon
px = ad['photons','particle_momentum_x'].to_ndarray()
py = ad['photons','particle_momentum_y'].to_ndarray()
pz = ad['photons','particle_momentum_z'].to_ndarray()
w  = ad['photons','particle_weight'].to_ndarray()
EPyt = EPyt + np.sum( (np.sqrt(px**2+py**2+pz**2)*scc.c)*w )

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
max_E = np.sqrt(np.amax(Ex**2+Ey**2+Ez**2))
max_Bx = np.amax(np.abs(Bx))
max_By = np.amax(np.abs(By))
max_Bz = np.amax(np.abs(Bz))
max_B = np.sqrt(np.amax(Bx**2+By**2+Bz**2))

print('Ex :',Ex)

# PART2: get results from reduced diagnostics

EFdata = np.genfromtxt("./diags/reducedfiles/EF.txt")
EPdata = np.genfromtxt("./diags/reducedfiles/EP.txt")
MFdata = np.genfromtxt("./diags/reducedfiles/MF.txt")
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

# PART3: print and assert

print('difference of field energy:', abs(EFyt-EF))
print('tolerance of field energy:', 1.0e-3)
print('difference of particle energy:', abs(EPyt-EP))
print('tolerance of particle energy:', 1.0e-8)

print('max Ex:', max_Ex,max_Exdata)
print('max Ey:', max_Ey,max_Eydata)
print('max Ez:', max_Ez,max_Ezdata)
print('max E:', max_E,max_Edata)
print('max Bx:', max_Bx,max_Bxdata)
print('max By:', max_By,max_Bydata)
print('max Bz:', max_Bz,max_Bzdata)
print('max B:', max_B,max_Bdata)

assert(abs(EFyt-EF) < 1.0e-3)
assert(abs(EPyt-EP) < 1.0e-8)
