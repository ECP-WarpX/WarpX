#! /usr/bin/env python

# Copyright 2021 Roelof Groenewald

# This script tests the MCC implementation for particles colliding
# with background neutrals. The test includes elastic scattering, excitation
# and ionization for electrons and elastic- and back-scattering for ions.
# The test simply runs a simulation seeded with a neutral plasma for a set
# number of steps and compares the density profiles of the electrons and ions
# against previously generated results. The true test of the MCC is done
# through comparison with literature results. In the same folder as this file
# there are 4 benchmark cases given comparing to results from Turner et al.
# (Phys. Plasmas 20, 013507, 2013).

# Possible running time: ~ 90 s

import numpy as np
import yt

ds = yt.load('background_mcc_plt00050')

data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)
elec_data = data['rho_electrons'].to_ndarray()
ion_data = data['rho_he_ions'].to_ndarray()

# this section saves the data for comparison
# should be commented out unless the reference data
# needs to be regenerated for some reason
# ref_data = np.zeros((np.shape(elec_data)[0], np.shape(elec_data)[1], 2))
# ref_data[:,:,0] = elec_data[:,:,0]
# ref_data[:,:,1] = ion_data[:,:,0]
# np.save('reference_results/test_ref_data.npy', ref_data)

# load reference data
ref_data = np.load('reference_results/test_ref_data.npy')

# compare results to reference
elec_rms = np.sqrt(np.mean((ref_data[:,:,0] - elec_data[:,:,0])**2))
ion_rms = np.sqrt(np.mean((ref_data[:,:,1] - ion_data[:,:,0])**2))

assert np.isclose(elec_rms, 0)
assert np.isclose(ion_rms, 0)
