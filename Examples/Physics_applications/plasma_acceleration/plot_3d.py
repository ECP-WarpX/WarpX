#!/usr/bin/env python3

# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Ryan Sandberg
# License: BSD-3-Clause-LBNL
#
# This is a script plots the wakefield and energy gain of a PWFA simulation.


from matplotlib import pyplot as plt

plt.rcParams.update({'font.size':16})
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

ts = OpenPMDTimeSeries('diags/diag1')

iteration=ts.iterations[-1]
select = {'y':[-3e-6,3e-6]}

pez, pex = ts.get_particle(species='plasma_e',
                           iteration=iteration,
                           var_list=['z','x'],
                           select=select)
pdz, pdx = ts.get_particle(species='driver',
                           iteration=iteration,
                           var_list=['z','x'])
bz, bx = ts.get_particle(species='beam',
                         iteration=iteration,
                         var_list=['z','x'])

fez, info = ts.get_field('E',
                         coord='z',
                         iteration=iteration,
                         slice_across='y'
                        )

psize = 0.2
xscale = 1e6
zscale = 1e2
ext = info.imshow_extent
ext = [ext[2]*zscale,ext[3]*zscale,ext[0]*xscale,ext[1]*xscale]
plt.figure()
ezplt = plt.imshow(fez.T/1e9,
                   extent=ext,
                   origin='lower',
                   aspect='auto',
                   cmap='RdBu_r',
                   vmin=-20,
                   vmax=20,
                   interpolation='sinc',
                  )
cb = plt.colorbar(ezplt,ax=plt.gca())
cb.set_label(r'$E_z$ (GV/m)')
plt.scatter(pez*zscale,pex*xscale,s=psize,c='g')
plt.scatter(pdz*zscale,pdx*xscale,s=psize,c='tab:purple')
plt.scatter(bz*zscale,bx*xscale,s=psize,c='k')
plt.ylim(-60,60)

plt.xlabel(r'Propagation distance z (cm)')
plt.ylabel(r'Transverse direction x ($\mu$m)')
plt.tight_layout()
plt.savefig('pwfa_wakefield.png')
plt.close()


data = np.genfromtxt("diags/reducedfiles/beamrel.txt")
plt.figure()
c = 2.998e8
plt.plot(c*data[:,1]*1e3,data[:,8]*0.511)
plt.ylabel('Mean energy, boosted frame (MeV)')
plt.xlabel(r'Propagation distance z (mm)')
plt.tight_layout()
plt.savefig('pwfa_beam_energy.png')
plt.close()
