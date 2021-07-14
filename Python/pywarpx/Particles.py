# Copyright 2017-2020 Andrew Myers, David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

particles = Bucket('particles', species_names=[], rigid_injected_species=[])
particles_list = []

electrons = Bucket('electrons')
electrons.charge = "-q_e"
electrons.mass = "m_e"
electrons.injection_style = None

positrons = Bucket('positrons')
positrons.charge = "q_e"
positrons.mass = "m_e"
positrons.injection_style = None

protons = Bucket('protons')
protons.charge = "q_e"
protons.mass = "m_p"
protons.injection_style = None

particle_dict = {'electrons':electrons,
                 'positrons':positrons,
                 'protons':protons
                 }

def newspecies(name):
    result = Bucket(name)
    particles_list.append(result)
    return result
