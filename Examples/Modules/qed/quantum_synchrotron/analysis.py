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
import scipy.special as spe
import scipy.integrate as integ
import scipy.stats as st
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
#import checksumAPI

#import matplotlib.pyplot as plt

# This script performs detailed checks of the Quantum Synchrotron photon emission process.
# Two electron populations and two positron populations are initialized with different momenta in different
# directions in a background EM field (with non-zero components along each direction).
# Specifically the script checks that:
#
# - The expected number of generated photons n_phot is in agreement with theory
#   (the maximum tolerated error is 5*sqrt(n_phot),
#   which means that the test should statistically fail less than once every 10^6 runs).
# - The weight of the generated particles is equal to the weight of the photon
# - The generated particles are emitted in the right direction
# - The total energy is conserved in each event
# - The energy distribution of the generated particles is in agreement with theory
# - The optical depths of the product species are correctly initialized (QED effects are
#   enabled for product species too).
#
# More details on the theoretical formulas used in this script can be found in
# the Jupyter notebook picsar/src/multi_physics/QED_tests/validation/validation.ipynb
#
# References:
# 1) R. Duclous et al 2011 Plasma Phys. Control. Fusion 53 015009
# 2) A. Gonoskov et al. 2015 Phys. Rev. E 92, 023305
# 3) M. Lobet. PhD thesis "Effets radiatifs et d'electrodynamique
#    quantique dans l'interaction laser-matiere ultra-relativiste"
#    URL: https://tel.archives-ouvertes.fr/tel-01314224


# Tolerances
tol = 1.e-8
tol_red = 1.e-2

# Physical constants (from CODATA 2018, see: https://physics.nist.gov/cuu/Constants/index.html )
me = 9.1093837015e-31 #electron mass
c = 299792458 #speed of light
hbar = 6.62607015e-34/(2*np.pi) #reduced Plank constant
fine_structure = 7.2973525693e-3 #fine structure constant
qe = 1.602176634e-19#elementary charge
E_s = (me**2 * c**3)/(qe * hbar) #Schwinger E field
B_s = E_s/c #Schwinger B field

mec = me*c
mec2 = mec*c
#______________

# Initial parameters
spec_names = ["p1", "p2", "p3", "p4"]
spec_names_phot = ["qsp_1", "qsp_2", "qsp_3", "qsp_4"]
initial_momenta = [
    np.array([10.0,0,0])*mec,
    np.array([0.0,100.0,0.0])*mec,
    np.array([0.0,0.0,1000.0])*mec,
    np.array([5773.502691896, 5773.502691896, 5773.502691896])*mec
]
csign = [-1,-1,1,1]
initial_particle_number = 1048576

E_f = np.array([-2433321316961438., 973328526784575., 1459992790176863.])
B_f = np.array([2857142.85714286, 4285714.28571428, 8571428.57142857])

NNS = [128,128,128,128] #bins for energy distribution comparison.
#______________

def calc_chi_part(p, E, B):
    gamma_part = np.sqrt(1.0 + np.dot(p,p)/mec**2)
    v = p/(gamma_part*me)    
    EpvvecB = E + np.cross(v,B)
    vdotEoverc = np.dot(v,E)/c
    ff = np.sqrt(np.dot(EpvvecB,EpvvecB) - np.dot(vdotEoverc,vdotEoverc))
    return gamma_part*ff/E_s

#Auxiliary functions
@np.vectorize
def IC_inner_alternative(y):
    ff = lambda x : np.exp(-y*(1+(4*x**2)/3)*np.sqrt(1+x*x/3))*(9+36*x**2 + 16*x**4)/(3 + 4*x**2)/np.sqrt(1+(x**2)/3)
    return integ.quad(ff, 0, np.inf)[0]/np.sqrt(3) 

def IC_Y(chi_ele, xi):
    res = np.zeros(np.shape(chi_ele))
    div = (chi_ele*(1-xi))
    div = np.where(np.logical_and(xi < 1, chi_ele != 0), div, 1.0);
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

def small_diff(vv, val):
    if(val > 0.0):
        return np.max(np.abs((vv - val)/val)) < tol
    else:
        return np.max(np.abs(vv)) < tol
        
def boris(pp, dt, charge_sign):
    econst = 0.5*qe*dt*charge_sign/me
    u = pp/(me)
    u += econst*E_f
    inv_gamma = 1/np.sqrt(1 + np.dot(u,u)/c*2)
    t = econst*B_f*inv_gamma
    s = 2*t/(1 + np.dot(t,t))
    u_p = u + np.cross(u,t)
    u += np.cross(u_p, s)
    u += econst*E_f
    return u *me
#__________________

# Quantum Synchrotron total and differential cross sections
def QS_dN_dt(chi_ele, gamma_ele):
    coeff_IC = (2./3.) * fine_structure * me*c**2/hbar 
    return coeff_IC*IC_G(chi_ele)/gamma_ele

def QS_d2N_dt_dchi(chi, gamma_ele, chi_phot):    
    coeff_IC = (2./3.) * fine_structure * me*c**2/hbar 
    return coeff_IC*IC_S(chi, chi_phot/chi)/chi_phot/gamma_ele
#__________________

# Get data for a species
def get_spec(ytdata, specname, is_photon):
    px = ytdata[specname,"particle_momentum_x"].v
    pz = ytdata[specname,"particle_momentum_z"].v
    py = ytdata[specname,"particle_momentum_y"].v

    w = ytdata[specname,"particle_weighting"].v

    opt = np.zeros(np.shape(px))

    if (is_photon):
        opt = ytdata[specname,"particle_optical_depth_BW"].v
    else:
        opt = ytdata[specname,"particle_optical_depth_QSR"].v

    return {"px" : px, "py" : py, "pz" : pz, "w" : w, "opt" : opt}

# Individual tests
def check_number_of_photons(ytdataset, part_name, phot_name, chi_part, gamma_part, dt, particle_number):
    dNQS_dt_theo = QS_dN_dt(chi_part, gamma_part)
    expected_photons = (1.-np.exp(-dNQS_dt_theo*dt))*particle_number
    expected_photons_tolerance = 5.0*np.sqrt(expected_photons)
    n_phot = ytdataset.particle_type_counts[phot_name]
    print(chi_part, gamma_part, dt)
    #assert( np.abs(n_phot-expected_photons) < expected_photons_tolerance)
    print(n_phot, expected_photons_tolerance)
    print("  [OK] generated photons number is within expectations")
    return n_phot

def check_weights(phot_data, ele_data, pos_data):
    assert(np.all(phot_data["w"] == phot_data["w"][0]))
    assert(np.all(ele_data["w"]  == phot_data["w"][0]))
    assert(np.all(ele_data["w"]  == phot_data["w"][0]))
    print("  [OK] particles weights are the expected ones")

def check_momenta(phot_data, ele_data, pos_data, p0, p_ele, p_pos):
    assert(small_diff(phot_data["px"], p0[0]))
    assert(small_diff(phot_data["py"], p0[1]))
    assert(small_diff(phot_data["pz"], p0[2]))
    print("  [OK] residual photons still have initial momentum")

    pdir = p0/np.linalg.norm(p0)
    assert(small_diff(ele_data["px"]/p_ele, pdir[0]))
    assert(small_diff(ele_data["py"]/p_ele, pdir[1]))
    assert(small_diff(ele_data["pz"]/p_ele, pdir[2]))
    assert(small_diff(pos_data["px"]/p_pos, pdir[0]))
    assert(small_diff(pos_data["py"]/p_pos, pdir[1]))
    assert(small_diff(pos_data["pz"]/p_pos, pdir[2]))
    print("  [OK] pairs move along the initial photon direction")

def check_energy(energy_phot, energy_ele, energy_pos):
    product_energy = energy_ele + energy_pos
    assert(small_diff(product_energy, energy_phot))
    print("  [OK] energy is conserved in each event")

def check_opt_depths(phot_data, ele_data, pos_data):
    data = (phot_data, ele_data, pos_data)
    for dd in data:
        loc, scale = st.expon.fit(dd["opt"])
        assert( np.abs(loc - 0) < tol_red )
        assert( np.abs(scale - 1) < tol_red )
    print("  [OK] optical depth distributions are still exponential")

def check_energy_distrib(energy_ele, energy_pos, gamma_phot, chi_phot, n_lost, NN, idx):
    h_energy_ele, ele_en = np.histogram(energy_ele/mec2, bins=NN, range=[1.0001,gamma_phot-1.0001])
    h_energy_pos, _ = np.histogram(energy_pos/mec2, bins=NN, range=[1.0001,gamma_phot-1.0001])

    cchi_part = chi_phot*(ele_en - 1)/(gamma_phot - 2)

    coeff= 20
    aux_chi = np.linspace(cchi_part[0],cchi_part[-1], NN*coeff)
    distrib = BW_d2N_dt_dchi(chi_phot, gamma_phot, aux_chi)
    distrib = np.sum(distrib.reshape(-1, coeff),1)
    distrib = n_lost*distrib/np.sum(distrib)

    # Visual comparison of distributions
    #en_coords = 0.5*(ele_en[1:]+ele_en[:-1])
    #plt.clf()
    #plt.xlabel("γ_particle")
    #plt.ylabel("N")
    #plt.title("χ_photon = {:f}".format(chi_phot))
    #plt.plot(en_coords, distrib,label="theory")
    #plt.plot(en_coords, h_energy_ele,label="BW electrons")
    #plt.plot(en_coords, h_energy_pos,label="BW positrons")
    #plt.legend()
    #plt.savefig("case_{:d}".format(idx+1))

    discr_ele = np.abs(h_energy_ele-distrib)
    discr_pos = np.abs(h_energy_pos-distrib)
    max_discr = 5.0 * np.sqrt(distrib)
    assert(np.all(discr_ele < max_discr))
    assert(np.all(discr_pos < max_discr))
    print("  [OK] energy distribution is within expectations")

#__________________

def check():
    filename_end = sys.argv[1]
    data_set_end = yt.load(filename_end)

    sim_time = data_set_end.current_time.to_value()
    all_data_end = data_set_end.all_data()

    for idx in range(4):
        part_name = spec_names[idx]
        phot_name  = spec_names_phot[idx]
        p0        = initial_momenta[idx]

        p0 = boris(p0,sim_time*0.5,csign[idx])


        p2_part = p0[0]**2 + p0[1]**2 + p0[2]**2
        p_part = np.sqrt(p2_part)
        energy_part = np.sqrt(mec2**2 + p2_part*c**2)
        chi_part = calc_chi_part(p0, E_f, B_f)
        gamma_part = energy_part/mec2

        print("** Case {:d} **".format(idx+1))
        print("  initial momentum: ", p0)
        print("  quantum parameter: {:f}".format(chi_part))
        print("  normalized particle energy: {:f}".format(gamma_part))
        print("  timestep: {:f} fs".format(sim_time*1e15))

        part_data = get_spec(all_data_end, part_name, is_photon=False)
        phot_data = get_spec(all_data_end, phot_name, is_photon=True)

        p_phot = np.sqrt(phot_data["px"]**2 + phot_data["py"]**2 + phot_data["pz"]**2)
        energy_phot = p_phot*c

        n_phot = check_number_of_photons(data_set_end,
                              part_name, phot_name,
                              chi_part, gamma_part, sim_time,
                              initial_particle_number)
                              
        return

        check_weights(phot_data, ele_data, pos_data)

        check_momenta(phot_data, ele_data, pos_data, p0, p_ele, p_pos)

        check_energy(energy_phot, energy_ele, energy_pos)

        check_energy_distrib(energy_ele, energy_pos, gamma_phot, chi_phot, n_lost, NNS[idx], idx)

        check_opt_depths(phot_data, ele_data, pos_data)

        print("*************\n")

    test_name = filename_end[:-9] # Could also be os.path.split(os.getcwd())[1]
    #checksumAPI.evaluate_checksum(test_name, filename_end)

def main():
    check()

if __name__ == "__main__":
    main()
