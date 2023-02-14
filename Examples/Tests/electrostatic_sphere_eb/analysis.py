#!/usr/bin/env python3

# Run the default regression test for the PICMI version of the EB test
# using the same reference file as for the non-PICMI test since the two
# tests are otherwise the same.

import sys

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI
# Check reduced diagnostics for charge on EB
import numpy as np
from scipy.constants import epsilon_0

data = np.loadtxt('diags/reducedfiles/eb_charge.txt')
data_eighth = np.loadtxt('diags/reducedfiles/eb_charge_one_eighth.txt')
q_sim = data[1,2]
q_sim_eighth = data_eighth[1,2]
# Theoretical charge on the embedded boundary, for sphere at potential phi_0
phi_0 = 1. # V
R = 0.1 # m
q_th = -4*np.pi*epsilon_0*phi_0*R
assert abs((q_sim-q_th)/q_th) < 0.03
assert abs((q_sim_eighth-q_th/8)/(q_th/8)) < 0.03

# The first step is the same as in inputs_3d
my_check = checksumAPI.evaluate_checksum(
    'ElectrostaticSphereEB', 'Python_ElectrostaticSphereEB_plt000001',
    do_particles=False, atol=1e-12
)

# The second step has the EB potential modified via Python
my_check = checksumAPI.evaluate_checksum(
    'Python_ElectrostaticSphereEB', 'Python_ElectrostaticSphereEB_plt000002',
    do_particles=False, atol=1e-12
)
