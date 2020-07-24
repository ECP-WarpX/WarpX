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
#import checksumAPI

# This script checks if the optical depth of electron and positrons undergoing the
# Quantum Synchrotron process behaves as expected. Two populations of positrons and
# two populations of electrons are are initialized with different momenta
# in a background EM field.
# The expected number of generated photons n_phot is computed for each case and is
# compared with the simulation results. The maximum tolerated error is 3*sqrt(n_phot).
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
schwinger_field_SI = electron_mass**2 * speed_of_light**3/(reduced_plank*elementary_charge)
#______________

# Initial parameters
spec_names = ["p1", "p2", "p3", "p4"]
mec = electron_mass*speed_of_light
p_begin = {
    "p1": np.array([100.0,0,0])*mec,
    "p2": np.array([0.0,1000.0,0.0])*mec,
    "p3": np.array([0.0,0.0,10000.0])*mec,
    "p4": np.array([57735.02691896, 57735.02691896, 57735.02691896])*mec
}
initial_particle_number = 65536
csign = [-1,-1,1,1]
#______________

def calc_chi_part(p, E, B):
    p = p / electron_mass / speed_of_light
    gamma = np.sqrt( 1 + np.dot(p,p))
    b = p/gamma
    v = speed_of_light*b
    EvdotB = E + np.cross(v,B)
    bdotE = np.dot(b,E)
    loc_field = gamma * np.sqrt(np.dot(EvdotB,EvdotB) - bdotE**2)
    return loc_field/schwinger_field_SI

#Auxiliary functions
@np.vectorize
def IC_inner(y):
    return integ.quad(lambda s: spe.kv(5./3.,s), y, np.inf)[0]

@np.vectorize
def IC_inner_alternative(y):
    ff = lambda x : np.exp(-y*(1+(4*x**2)/3)*np.sqrt(1+x*x/3))*(9+36*x**2 + 16*x**4)/(3 + 4*x**2)/np.sqrt(1+(x**2)/3)
    return integ.quad(ff, 0, np.inf)[0]/np.sqrt(3)

def IC_Y(chi_ele, xi):
    res = np.zeros(np.shape(chi_ele))
    div = (chi_ele*(1-xi))
    div = np.where(np.logical_and(xi < 1, chi_ele != 0), div, 1.0)
    res = (2/3)*np.where(np.logical_and(xi < 1, chi_ele != 0), xi/div, np.inf)
    return res

def IC_S(chi_ele, xi):
    Y = IC_Y(chi_ele, xi)
    coeff = np.sqrt(3)/2.0/np.pi
    first = IC_inner_alternative(Y)
    div = np.where(xi == 1, 1.0, 1.0/(1-xi)  )
    res = np.where(np.logical_or(xi == 1, xi == 0), 0.0,
        coeff*xi*( first  + (xi**2 * spe.kv(2./3.,Y)*div )  ) )
    return res

def IC_SXI(chi_ele, xi):
    div = np.where(xi != 0, xi, 1.0)
    return np.where(xi != 0, IC_S(chi_ele, xi)/div, np.inf)

@np.vectorize
def IC_G(chi_ele):
    return integ.quad(lambda xi: IC_SXI(chi_ele, xi), 0, 1)[0]
#__________________

# Quantum-Synchrotron total cross section
def dNQS_dt(chi_ele, gamma_ele):
    coeff_IC = (2./3.) * fine_structure_constant * electron_mass*speed_of_light**2/reduced_plank
    return coeff_IC*IC_G(chi_ele)/gamma_ele
#__________________

def check():
    filename_end = sys.argv[1]
    data_set_end = yt.load(filename_end)

    sim_time = data_set_end.current_time.to_value()
    all_data_end = data_set_end.all_data()

    for i,name in enumerate(spec_names):
        opt_depth = all_data_end[name, 'particle_optical_depth_QSR']
        which_spec = "qsp_"+str(i+1)
        phot_num = data_set_end.particle_type_counts[which_spec]
        print(i,phot_num)

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

        #check if number of generated photons is n0*(1 - exp(-rate*t))
        print("*",p_begin[name], E_f, B_f)

        print(sim_time)

        def boris(pp, dt, charge_sign):
            econst = 0.5*elementary_charge*dt*charge_sign/electron_mass
            u = pp/(electron_mass)
            u += econst*E_f
            inv_gamma = 1/np.sqrt(1 + np.dot(u,u)/speed_of_light*2)
            t = econst*B_f*inv_gamma
            s = 2*t/(1 + np.dot(t,t))
            u_p = u + np.cross(u,t)
            u += np.cross(u_p, s)
            u += econst*E_f
            return u *electron_mass

        p_begin[name] = boris(p_begin[name],-sim_time*0.5,csign[i])

        p = p_begin[name] / electron_mass / speed_of_light
        print(p_begin[name], E_f, B_f)
        gamma = np.sqrt( 1 + np.dot(p,p))
        dNQS_dt_theo = dNQS_dt(
            calc_chi_part(p_begin[name], E_f, B_f), gamma)
        exp_phot= dNQS_dt_theo*sim_time*initial_particle_number
        discrepancy = np.abs(phot_num-exp_phot)
        print("photons = {:d}, exp_photons= {:f}, discrepancy = {:f}, max_discrepancy = {:f}".format(phot_num, exp_phot, discrepancy, 3.0*np.sqrt(exp_phot)) )
        #assert(np.abs(exp_lost-lost) < 3.0*np.sqrt(exp_lost))
        ###

    test_name = filename_end[:-9] # Could also be os.path.split(os.getcwd())[1]
    #checksumAPI.evaluate_checksum(test_name, filename_end)

def main():
    check()

if __name__ == "__main__":
    main()
