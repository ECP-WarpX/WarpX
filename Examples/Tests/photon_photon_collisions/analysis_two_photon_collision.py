#!/usr/bin/env python3

import os
import sys

import numpy as np
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
from analysis_base import *
import checksumAPI
from scipy.integrate import cumtrapz


# get input parameters from warpx_used_inputs
def get_input_parameters():
    with open('./warpx_used_inputs', 'rt') as f:
        lines = f.readlines()
        for line in lines:
            if 'photonA.single_particle_weight' in line:
                w1 = find_num_in_line(line)
            if 'photonB.single_particle_weight' in line:
                w2 = find_num_in_line(line)
    return (w1, w2)

# check that the photons have been completely transformed into pairs:
# because the fusion multiplier is 1, as soon as a linear Breit-Wheeler event occurs,
# the two photons must disappear and 2 electron-positron pairs must be generated
def check_final_macroparticles():
    (w1, w2) = get_input_parameters()
    macro_photonA_number = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,3]
    macro_photonA_weight = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,8]
    macro_photonB_number = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,4]
    macro_photonB_weight = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,9]
    assert(macro_photonA_number == macro_photonB_number == 0.)
    assert(macro_photonA_weight == macro_photonB_weight == 0.)
    macro_positron_number = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,6]
    macro_positron_weight = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,11]
    macro_electron_number = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,5]
    macro_electron_weight = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[-1,10]
    assert(macro_electron_number == macro_positron_number == 2.)
    assert(macro_electron_weight == macro_positron_weight == w1 == w2)

def main():
    check_final_macroparticles()
    check_energy_conservation()
    check_momentum_conservation()
    check_charge_conservation()

    test_name = os.path.split(os.getcwd())[1]
    filename_end = sys.argv[1]
    checksumAPI.evaluate_checksum(test_name, filename_end)

if __name__ == "__main__":
    main()
