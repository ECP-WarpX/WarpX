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
from scipy.constants import c, e, epsilon_0, m_e, m_p
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
ne = float(ds0.parameters['ion_stopping_constant.background_density'])
Te = float(ds0.parameters['ion_stopping_constant.background_temperature'])
Zb = 1.  # Ion charge state
ion_mass = m_p

# For ions1, the background parameters are constants
vx = ad0[('ions1', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad0[('ions1', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad0[('ions1', 'particle_momentum_z')].to_ndarray()/m_p
EE1 = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

# For ions2, the background parameters are parsed
xx = ad0[('ions2', 'particle_position_x')].to_ndarray()/m_p
yy = ad0[('ions2', 'particle_position_y')].to_ndarray()/m_p
ne2 = np.where(xx > 0., 1.e20, 1.e21)
Te2 = np.where(yy > 0., 1., 2.)

vx = ad0[('ions2', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad0[('ions2', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad0[('ions2', 'particle_momentum_z')].to_ndarray()/m_p
EE2 = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

ds = yt.load(last_filename)
ad = ds.all_data()
dt = ds.current_time.to_value()/int(last_it)

# Step through the same number of time steps
a_EE1 = EE1
a_EE2 = EE2
for it in range(int(last_it)):
    dEdt1 = stopping_from_electrons(ne, Te, Zb, ion_mass)
    a_EE1 *= np.exp(dEdt1*dt)
    dEdt2 = stopping_from_electrons(ne2, Te2, Zb, ion_mass)
    a_EE2 *= np.exp(dEdt2*dt)

# Fetch the final particle data
vx = ad[('ions1', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad[('ions1', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad[('ions1', 'particle_momentum_z')].to_ndarray()/m_p
EE1 = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

vx = ad[('ions2', 'particle_momentum_x')].to_ndarray()/m_p
vy = ad[('ions2', 'particle_momentum_y')].to_ndarray()/m_p
vz = ad[('ions2', 'particle_momentum_z')].to_ndarray()/m_p
EE2 = 0.5*m_p*(vx**2 + vy**2 + vz**2)/e

error1 = np.abs(EE1 - a_EE1)
error2 = np.abs(EE2 - a_EE2)
print('error with constant = ', error1)
print('error with parsed = ', error2)
print('tolerance = ', tolerance)

assert np.all(error1 < tolerance)
assert np.all(error2 < tolerance)

test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, last_filename)
