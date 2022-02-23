#!/usr/bin/env python3

# Copyright 2022 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the ion stopping model.
# It uses the same stopping power formula that
# is used in the C++ to check the resulting
# particle energies.

import os
import re
import sys

import numpy as np
from scipy.constants import c, e, m_e, epsilon_0, m_p
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Tolerance on the error in the final energy (in eV)
tolerance = 1.e-7

last_filename = sys.argv[1]

# Remove trailing '/' from file name, if necessary
last_filename.rstrip('/')
# Find last iteration in file name, such as 'test_name_plt000001' (last_it = '000001')
last_it = re.search('\d+$', last_filename).group()
# Find output prefix in file name, such as 'test_name_plt000001' (prefix = 'test_name_plt')
prefix = last_filename[:-len(last_it)]

def stopping_from_electrons(ne, Te, Zb, ion_mass):
    """Calculate the coefficient in equation 14.13 from
    "Introduction to Plasma Physics", Goldston and Rutherford.
    ne: electron density
    Te: electron temperature (eV)
    Zb: ion charge state
    ion_mass: (kg)
    """
    vthe = np.sqrt(3.*Te*e/m_e)
    wpe = np.sqrt(ne*e**2/(epsilon_0*m_e))
    lambdadb = vthe/wpe
    loglambda = np.log((12.*np.pi/Zb)*(ne*lambdadb**3))
    # Goldston's equation 14.13
    dEdt = - np.sqrt(2.)*ne*Zb**2*e**4*np.sqrt(m_e)*loglambda/(6.*np.pi**1.5*epsilon_0**2*ion_mass*(Te*e)**1.5)
    return dEdt

# Fetch background parameters and inital particle data
ds0 = yt.load(f'{prefix}{len(last_it)*"0"}')
ad0 = ds0.all_data()
ne = float(ds0.parameters['ion_stopping.background_density'])
Te = float(ds0.parameters['ion_stopping.background_temperature'])
Zb = 1.  # Ion charge state
ion_mass = m_p

vx = ad0[('ions', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad0[('ions', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad0[('ions', 'particle_momentum_z')].to_ndarray()/m_p
EE = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

ds = yt.load(last_filename)
ad = ds.all_data()
dt = ds.current_time.to_value()/int(last_it)

# Step through the same number of time steps
a_EE = EE
for it in range(int(last_it)):
    dEdt = stopping_from_electrons(ne, Te, Zb, ion_mass)
    a_EE *= np.exp(dEdt*dt)

# Fetch the final particle data
vx = ad[('ions', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad[('ions', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad[('ions', 'particle_momentum_z')].to_ndarray()/m_p
EE = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

error = np.abs(EE - a_EE)
print('error = ', error)
print('tolerance = ', tolerance)

assert np.all(error < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.reset_benchmark(test_name, last_filename)
