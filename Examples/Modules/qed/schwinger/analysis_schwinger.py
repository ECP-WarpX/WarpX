#! /usr/bin/env python

# Copyright 2020 Luca Fedeli, Neil Zaim
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

## This module is imported by each of the scripts used to analyze the Schwinger tests.

## The pair production rate is calculated using the formula described in
## Bulanov, S. S., et al. Physical review letters 104.22 (2010): 220404.

import yt
import numpy as np
import sys
import re

# define some parameters

c = 299792458.
m_e = 9.10938356e-31
e = 1.6021766208e-19
hbar = 1.054571800e-34
E_S = m_e**2*c**3/e/hbar # Schwinger field

dV = (1.e-6)**3 # total simulation volume
dt = 2.407291502e-16
filename = sys.argv[1]

Ex_test = 0.
Ey_test = 0.
Ez_test = 0.
Bx_test = 0.
By_test = 0.
Bz_test = 0.

# Find which test we are doing
test_number = re.search( 'qed_schwinger([1234])', filename ).group(1)
if test_number == '1':
    Ex_test = 1.e16
    Bx_test = 16792888.570516706
    By_test = 5256650.141557486
    Bz_test = 18363530.799561853
elif test_number == '2':
    Ex_test = 1.e18
    Bx_test = 1679288857.0516706
    By_test = 525665014.1557486
    Bz_test = 1836353079.9561853
elif test_number == '3':
    Ey_test = 1.0321239524474501e+17
elif test_number == '4':
    Ez_test = 2.5e+20
    By_test = 833910154604.3563
else:
    assert(False)

def calculate_rate(Ex,Ey,Ez,Bx,By,Bz):
## Calculate theoretical pair production rate from EM field value

    E_squared = Ex**2 + Ey**2 + Ez**2
    H_squared = c**2*(Bx**2 + By**2 + Bz**2)

    F = (E_squared - H_squared)/2.
    G = c*(Ex*Bx + Ey*By + Ez*Bz)

    epsilon = np.sqrt(np.sqrt(F**2+G**2)+F)/E_S
    eta = np.sqrt(np.sqrt(F**2+G**2)-F)/E_S

    if(epsilon != 0. and eta != 0.):
        return  e**2*E_S**2/4./np.pi**2/c/hbar**2*epsilon*eta/np.tanh(np.pi*eta/epsilon)*np.exp(-np.pi/epsilon)
    elif (epsilon == 0.):
        return 0.
    else:
        return  e**2*E_S**2/4./np.pi**2/c/hbar**2*epsilon**2/np.pi*np.exp(-np.pi/epsilon)


def do_analysis(Ex,Ey,Ez,Bx,By,Bz):

    data_set = yt.load(filename)

    expected_total_physical_pairs_created = dV*dt*calculate_rate(Ex,Ey,Ez,Bx,By,Bz)
    if expected_total_physical_pairs_created < 0.01:
        np_ele = data_set.particle_type_counts["ele_schwinger"] if \
            "ele_schwinger" in data_set.particle_type_counts.keys() else 0
        np_pos = data_set.particle_type_counts["pos_schwinger"] if \
            "pos_schwinger" in data_set.particle_type_counts.keys() else 0
        assert(np_ele == 0 and np_pos == 0)
        ## Assert whether pairs are created or not.

    else:
        all_data = data_set.all_data()

        ele_data = all_data["ele_schwinger",'particle_weight']
        pos_data = all_data["pos_schwinger",'particle_weight']

        std_total_physical_pairs_created = np.sqrt(expected_total_physical_pairs_created)

        # This first assert only works if a single tile is used in the simulation
        assert(np.array_equal(ele_data,pos_data))
        # 5 sigma test that has an intrisic probability to fail of 1 over ~2 millions
        assert(np.abs(np.sum(ele_data)-expected_total_physical_pairs_created)<5*std_total_physical_pairs_created)

do_analysis(Ex_test, Ey_test, Ez_test, Bx_test, By_test, Bz_test)
