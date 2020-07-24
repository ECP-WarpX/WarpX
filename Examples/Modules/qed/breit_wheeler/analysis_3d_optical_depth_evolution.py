#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# -*- coding: utf-8 -*-

import yt
import numpy as np
import sys
import math as m
import scipy.special as spe
import scipy.integrate as integ
import scipy.stats as st
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# This script checks if the optical depth of photons undergoing the
# Breit Wheeler process behaves as expected. Four populations of photons
# are initialized with different momenta in a background EM field.
# The expected number of lost photons n_lost is computed for each case and is
# compared with the simulation results. The maximum tolerated error is 3*sqrt(n_lost).
# The script checks also that the optical depth distribution is still exponential.
#
# References:
# 1) R. Duclous et al 2011 Plasma Phys. Control. Fusion 53 015009
# 2) A. Gonoskov et al. 2015 Phys. Rev. E 92, 023305
# 3) M. Lobet. PhD thesis "Effets radiatifs et d'electrodynamique
#    quantique dans l'interaction laser-matiere ultra-relativiste"
#    URL: https://tel.archives-ouvertes.fr/tel-01314224


# Tolerance
tol = 1.e-2

# EM fields
E_f = np.array([-2433321316961438, 973328526784575, 1459992790176863])
B_f = np.array([2857142.85714286, 4285714.28571428, 8571428.57142857])

# Physical constants
electron_mass = 9.10938356e-31
elementary_charge = 1.6021766208e-19
speed_of_light = 299792458
reduced_plank = 1.054571800e-34
vacuum_permittivity =  8.854187817e-12
fine_structure_constant =  0.0072973525664
classical_elec_radius = (1./4./np.pi/vacuum_permittivity)*( elementary_charge**2 / (electron_mass * speed_of_light**2))
lambda_ref = 1.0e-6
field_reference = 2.0 * np.pi * electron_mass*speed_of_light * speed_of_light / (elementary_charge*lambda_ref)
schwinger_field_SI = electron_mass**2 * speed_of_light**3/(reduced_plank*elementary_charge)
schwinger_field_norm = electron_mass*speed_of_light*lambda_ref/(2.0*reduced_plank*m.pi)
#______________

# Initial parameters
spec_names = ["p1", "p2", "p3", "p4"]
mec = electron_mass*speed_of_light
p_begin = {
    "p1": np.array([2000.0,0,0])*mec,
    "p2": np.array([0.0,5000.0,0.0])*mec,
    "p3": np.array([0.0,0.0,10000.0])*mec,
    "p4": np.array([57735.02691896, 57735.02691896, 57735.02691896])*mec
}
initial_particle_number = 65536
#______________

def calc_chi_gamma(p, E, B):
    p = p / electron_mass / speed_of_light
    E = E / field_reference
    B = B * speed_of_light / field_reference
    gamma_phot = np.linalg.norm(p)
    c = p/gamma_phot
    loc_field = gamma_phot * np.linalg.norm( E - np.dot(c,E)*c + np.cross(c,B))
    return loc_field/schwinger_field_norm
    
    gamma_phot = np.linalg.norm(p, axis=1)/mec  
    u = p/(np.linalg.norm(p, axis=1)[:,None])
    udotE = np.einsum('ij,ij->i',u,E)
    Eperp = E - udotE[:,None]*u
    ccrossB = c*np.cross(u,B)
    loc_field = gamma_phot * np.linalg.norm( Eperp + ccrossB, axis=1) 
    return loc_field/E_s

#Auxiliary functions
@np.vectorize
def BW_inner(x):
    return integ.quad(lambda s: np.sqrt(s)*spe.kv(1./3., 2./3. * s**(3./2.)), x, np.inf)[0] 

def BW_X(chi_phot, chi_ele):
    res = np.zeros(np.shape(chi_phot))
    div = (chi_ele*(chi_phot-chi_ele))
    div = np.where(np.logical_and(chi_phot > chi_ele, chi_ele != 0), div, 1.0);
    res = np.where(np.logical_and(chi_phot > chi_ele, chi_ele != 0), np.power(chi_phot/div, 2./3.), np.inf)
    return res

def BW_F(chi_phot, chi_ele):
    X = BW_X(chi_phot, chi_ele)
    res = np.where(np.logical_or(chi_phot == chi_ele, chi_ele == 0), 0.0, 
         BW_inner(X) - (2.0 - chi_phot* X**(3./2.))*spe.kv(2./3., 2./3. * X**(3./2.)) )
    return res

@np.vectorize
def BW_T(chi_phot):
    coeff = 1./(np.pi * np.sqrt(3.) * (chi_phot**2))
    return coeff*integ.quad(lambda chi_ele: BW_F(chi_phot, chi_ele), 0, chi_phot)[0]
#__________________

# Breit-Wheeler total cross section
def dNBW_dt(chi_phot, energy_phot):
    energy_phot = energy_phot/electron_mass/speed_of_light/speed_of_light
    return ((electron_mass*(speed_of_light)**2)*fine_structure_constant/reduced_plank)*(chi_phot/energy_phot)*BW_T(chi_phot)
#__________________

def check():
    filename_end = sys.argv[1]
    data_set_end = yt.load(filename_end)

    sim_time = data_set_end.current_time.to_value()
    all_data_end = data_set_end.all_data()

    for name in spec_names:
        opt_depth = all_data_end[name, 'particle_optical_depth_BW']

        #check that the distribution is still exponential with scale 1 and loc 0
        opt_depth_loc, opt_depth_scale = st.expon.fit(opt_depth)
        exp_loc = 0.0
        exp_scale = 1.0
        loc_discrepancy = np.abs(opt_depth_loc-exp_loc)
        scale_discrepancy = np.abs(opt_depth_scale-exp_scale)
        print("tolerance_rel: " + str(tol))
        print("species " + name)
        print("exp distrib loc tol = " + str(tol))
        print("exp distrib loc discrepancy = " + str(loc_discrepancy))
        assert(loc_discrepancy < tol)
        print("exp distrib scale tol = " + str(tol))
        print("exp distrib scale discrepancy = " + str(scale_discrepancy/exp_scale))
        assert(scale_discrepancy/exp_scale < tol)
        ###

        #check if number of lost photons is (n0* (1 - exp(-rate*t)) )
        dNBW_dt_theo = dNBW_dt(
            calc_chi_gamma(p_begin[name], E_f, B_f),
                np.linalg.norm(p_begin[name]*speed_of_light))
        exp_lost= initial_particle_number*(1.0 - np.exp(-dNBW_dt_theo*sim_time))
        lost =  initial_particle_number-np.size(opt_depth)
        discrepancy_lost = np.abs(exp_lost-lost)
        print("lost = {:d}, exp_lost = {:f}, discrepancy = {:f}, max_discrepancy = {:f}".format(lost, exp_lost, lost-exp_lost, 3.0*np.sqrt(exp_lost)) )
        assert(np.abs(exp_lost-lost) < 3.0*np.sqrt(exp_lost))
        ###

    test_name = filename_end[:-9] # Could also be os.path.split(os.getcwd())[1]
    checksumAPI.evaluate_checksum(test_name, filename_end)

def main():
    check()

if __name__ == "__main__":
    main()
