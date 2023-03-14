#!/usr/bin/env python3

#####################################################
# examples of usage of the reduced diagnostic tools #
#####################################################

import sys 

sys.path.append('../../../Python/pywarpx/diagnostics/data/')  
mydir = '.'

##################
# ParticleEnergy #
##################

import ParticleEnergy as PE
fname = 'EP' # do not add .txt 
diag = PE.ParticleEnergyData(mydir, fname)

# extract some general info
times = diag.get_times()
print('available times\n {}'.format(times))

steps = diag.get_steps()
print('available steps\n {}'.format(steps))

nspecies = diag.get_nspecies()
print('number of species = {}'.format(nspecies))

species_names = diag.get_species_names()
print('species names = {}'.format(species_names))

print('\n')

# extract data 
data = diag.get_data('electrons', times='all')
print('asking for total electron energy at all available times')
print(data)
print('\n')

data = diag.get_data('protons', steps='all')
print('asking for total proton energy at all available steps')
print(data)
print('\n')

data = diag.get_data('photons', 'electrons', times=None)
print('asking for total photon and electron energy at all available times')
print(data)
print('\n')

data = diag.get_data('electrons_mean', 'electrons', steps=None)
print('asking for mean and total electron energy at all available times')
print(data)
print('\n')

data = diag.get_data('protons_mean', 'electrons')
print('asking for mean proton energy and total electron energy at all available times')
print(data)
print('\n')

data = diag.get_data('protons', 'electrons', steps=200)
print('asking for proton and electron total energy at step 200')
print(data)
print('\n')

data = diag.get_data('photons_mean', 'photons', steps=[0,200])
print('asking for mean and total photon energy at steps 0 and 200')
print(data)
print('\n')

data = diag.get_data('photons', 'protons', 'electrons', times=10e-9)
print('asking for total photon, proton and electron energy at the time closest to 10 ns')
print(data)
print('\n')

data = diag.get_data('electrons', times=[2e-9, 4e-9, 9e-9])
print('asking for total electron energy at the times closest to 2, 4 and 9 ns')
print(data)
print('\n')

data = diag.get_data('electrons', times=[0,20,40,100])
print('asking for total electron energy at steps 0, 20, 40 and 100')
print(data)
print('\n')


####################
# ParticleMomentum #
####################

import ParticleMomentum as PP
fname = 'PP' # do not add .txt 
diag = PP.ParticleMomentumData(mydir, fname)

print(diag.get_valid_args())

data = diag.get_data('electrons_x', 'electrons_mean_x', steps='all')
print('asking for total and mean electron momentum x at all steps')
print(data)
print('\n')


###############
# FieldEnergy #
###############

import FieldEnergy as FE
fname = 'EF' # do not add .txt 
diag = FE.FieldEnergyData(mydir, fname)

print(diag.get_valid_args())

data = diag.get_data('E_lev1', 'B_lev0', steps='all')
print('asking for energy stored in E at level 1 and B at level 0 at all steps')
print(data)
print('\n')


################
# FieldMaximum #
################

import FieldMaximum as FM
fname = 'MF' # do not add .txt 
diag = FM.FieldMaximumData(mydir, fname)

print(diag.get_valid_args())

data = diag.get_data('E_lev1', 'B_lev0', steps='all')
print('asking for total electron energy at all steps')
print(data)
print('\n')

##################
# ParticleNumber #
##################

import ParticleNumber as PN
fname = 'NP' # do not add .txt 
diag = PN.ParticleNumberData(mydir, fname)

print(diag.get_valid_args())

species_names = diag.get_species_names()
print('species names = {}'.format(species_names))


data = diag.get_data('protons_macroparticles', 'photons_weight', steps='all')
print('asking for total electron energy at all steps')
print(data)
print('\n')







