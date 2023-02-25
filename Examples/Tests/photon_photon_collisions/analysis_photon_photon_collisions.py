#!/usr/bin/env python3

import os
import re
import sys

import numpy as np
from scipy.constants import (c, m_e, physical_constants, pi)
import yt

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI
from scipy.integrate import cumtrapz

# some constants
r_e = physical_constants["classical electron radius"][0]
default_tol = 5.e-11

# check the type of test
warpx_used_inputs = open('./warpx_used_inputs', 'r').read()
if re.search('photonA.injection_style = SingleParticle', warpx_used_inputs):
    test = 'two_photon_collision_3d'
elif re.search('photonA.momentum_distribution_type = constant', warpx_used_inputs):
    test = 'many_photons_collisions_3d'

# extract numbers from a string
def find_num_in_line(line):
    items = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', line)
    fitems = [float(it) for it in items]
    if len(fitems)==1:
        return fitems[0]
    else:
        return fitems

# get input parameters from warpx_used_inputs
def get_input_parameters(test):
    if test == 'two_photon_collision_3d':
        with open('./warpx_used_inputs', 'rt') as f:
            lines = f.readlines()
            for line in lines:
                if 'photonA.single_particle_weight' in line:
                    w1 = find_num_in_line(line)
                if 'photonB.single_particle_weight' in line:
                    w2 = find_num_in_line(line)
        return w1, w2
    elif test == 'many_photons_collisions_3d':
        with open('./warpx_used_inputs', 'rt') as f:
            lines = f.readlines()
            for line in lines:
                if 'warpx.cfl' in line:
                    cfl = find_num_in_line(line)
                if 'max_step' in line:
                    num_steps = find_num_in_line(line)
                if 'geometry.prob_lo' in line:
                    xmin, ymin, zmin = find_num_in_line(line)
                if 'geometry.prob_hi' in line:
                    xmax, ymax, zmax = find_num_in_line(line)
                if 'photonA.ux' in line:
                    uAx = find_num_in_line(line)
                if 'photonA.uy' in line:
                    uAy = find_num_in_line(line)
                if 'photonA.uz' in line:
                    uAz = find_num_in_line(line)
                if 'photonA.density' in line:
                    dens = find_num_in_line(line)
                if 'photonB.density' in line:
                    dens = find_num_in_line(line)
                if 'photonB.ux' in line:
                    uBx = find_num_in_line(line)
                if 'photonB.uy' in line:
                    uBy = find_num_in_line(line)
                if 'photonB.uz' in line:
                    uBz = find_num_in_line(line)
                if 'amr.n_cell' in line:
                    nx, ny, nz = find_num_in_line(line)
                if 'LBW.event_multiplier' in line:
                    event_multiplier = find_num_in_line(line)

        pAx, pAy, pAz = uAx * m_e * c, uAy * m_e * c, uAz * m_e * c
        pBx, pBy, pBz = uBx * m_e * c, uBy * m_e * c, uBz * m_e * c
        EA, EB = np.sqrt(pAx**2+pAy**2+pAz**2) * c, np.sqrt(pBx**2+pBy**2+pBz**2) * c
        dx = (xmax-xmin)/nx
        dy = (ymax-ymin)/ny
        dz = (zmax-zmin)/nz
        dt = cfl / c / np.sqrt(1./dx**2+1./dy**2+1./dz**2) # works for Yee solver
        pA = np.sqrt(pAx**2 + pAy**2 + pAz**2)
        pB = np.sqrt(pBx**2 + pBy**2 + pBz**2)
        cos_ang = (pAx*pBx+pAy*pBy+pAz*pBz)/(pA*pB)
        theta = np.arccos(cos_ang)
        E_star = np.sqrt(0.5*c**2*pA*pB*(1.- cos_ang))
        V = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
        NA0 = dens * V
        NB0 = dens * V
        return EA, EB, theta, dt, V, num_steps, NA0, NB0

# check that the photons have been completely transformed into pairs:
# because the fusion multiplier is 1, as soon as a linear Breit-Wheeler event occurs,
# the two photons must disappear and 2 electron-positron pairs must be generated
def check_final_macroparticles(test):
    if test == 'two_photon_collision_3d':
        w1, w2 = get_input_parameters(test)
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

def check_energy_conservation():
    ekin_data = np.loadtxt('diags/reducedfiles/ParticleEnergy.txt')
    ekin_photonA = ekin_data[:,3]
    ekin_photonB = ekin_data[:,4]
    ekin_electron = ekin_data[:,5]
    ekin_positron = ekin_data[:,6]
    num_data = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')
    num_phys_electron = num_data[:,10]
    num_phys_positron = num_data[:,11]
    field_energy = np.loadtxt('diags/reducedfiles/FieldEnergy.txt')[:,2]
    total_energy = ekin_photonA + ekin_photonB + ekin_electron + ekin_positron + m_e*c**2*(num_phys_electron + num_phys_positron) + field_energy
    assert(np.all(np.isclose(total_energy, total_energy[0], rtol=default_tol, atol=0.)))

def check_momentum_conservation():
    total_momentum_x = np.loadtxt('diags/reducedfiles/ParticleMomentum.txt')[:,2]
    total_momentum_y = np.loadtxt('diags/reducedfiles/ParticleMomentum.txt')[:,3]
    total_momentum_z = np.loadtxt('diags/reducedfiles/ParticleMomentum.txt')[:,4]
    assert(np.all(np.isclose(total_momentum_x, total_momentum_x[0], rtol=default_tol, atol=0.)))
    assert(np.all(np.isclose(total_momentum_y, total_momentum_y[0], rtol=default_tol, atol=0.)))
    assert(np.all(np.isclose(total_momentum_z, total_momentum_z[0], rtol=default_tol, atol=0.)))

def check_charge_conservation():
    rho_max = np.loadtxt('diags/reducedfiles/RhoMaximum.txt')[:,2]
    assert(np.all(np.isclose(rho_max, rho_max[0], rtol=default_tol, atol=0.)))
    #series = io.Series("diags/fields_particles/openpmd_%T.bp",io.Access.read_only)
    #iterations = np.asarray(series.iterations)
    #start = series.iterations[iterations[0]]
    #rho_start = start.meshes["rho"][io.Mesh_Record_Component.SCALAR].load_chunk()
    #end = series.iterations[iterations[-1]]
    #rho_end = end.meshes["rho"][io.Mesh_Record_Component.SCALAR].load_chunk()
    #series.flush()

    filename_end = sys.argv[1]
    filename_start = filename_end[:-4] + '0000'
    ds_end = yt.load(filename_end)
    ds_start = yt.load(filename_start)
    field_data_end = ds_end.covering_grid(level=0, left_edge=ds_end.domain_left_edge,
                                          dims=ds_end.domain_dimensions)
    field_data_start = ds_start.covering_grid(level=0, left_edge=ds_start.domain_left_edge,
                                              dims=ds_start.domain_dimensions)
    rho_start = field_data_start["rho"].to_ndarray()
    rho_end = field_data_end["rho"].to_ndarray()
    assert(np.all(np.isclose(rho_start, rho_end, rtol=default_tol, atol=0.)))

def cross_section(E1_lab, E2_lab, theta):
    s = E1_lab*E2_lab/(2.*m_e**2*c**4)*(1.-np.cos(theta))
    beta = np.sqrt(1.-1./s)
    factor1 = 0.5*pi*r_e**2*(1.-beta**2)
    term1 = (3.-beta**4)*np.log((1.+beta)/(1.-beta))
    term2 = 2*beta*(beta**2-2.)
    factor2 = term1 + term2
    sigma = factor1*factor2
    return sigma

def check_pair_rate(test):
    if test == 'many_photons_collisions_3d':
        EA_lab, EB_lab, theta, dt, V, num_steps, NA0, NB0 = get_input_parameters(test)

        t = np.arange(num_steps+1)*dt
        sigma = cross_section(EA_lab, EB_lab, theta)

        # estimated number of real photons of species photonA in time
        NA_est = NA0 / (1. + 2.*sigma* c * t / V * NA0)
        # number of <<real>> photons of species photonA in time from simulation
        NA = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[:,8]

        # estimated number of real photons of species photonB in time
        NB_est = NB0 / (1. + 2.*sigma* c * t / V * NB0)
        # number of <<real>> photons of species photonA in time from simulation
        NB = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[:,9]

        Nplus_est = 2.*sigma*c/V*cumtrapz(NA_est*NB_est, x=t, dx=dt, initial=0)
        # number of <<real>> positrons in time from simulation
        Nplus = np.loadtxt('diags/reducedfiles/ParticleNumber.txt')[:,11]

        assert(np.all(np.isclose(Nplus_est, Nplus, rtol=1e-1, atol=0.)))
        assert(np.all(np.isclose(NA_est, NA, rtol=1e-1, atol=0.)))
        assert(np.all(np.isclose(NB_est, NB, rtol=1e-1, atol=0.)))

def main():
    check_final_macroparticles(test)
    check_energy_conservation()
    check_momentum_conservation()
    check_charge_conservation()
    check_pair_rate(test)

    test_name = os.path.split(os.getcwd())[1]
    filename_end = sys.argv[1]
    checksumAPI.evaluate_checksum(test_name, filename_end)

if __name__ == "__main__":
    main()
