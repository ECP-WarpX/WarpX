#!/usr/bin/env python3

# Copyright 2024 Arianna Formenti
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import c, eV, m_e, micro, milli, nano

#sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
#import checksumAPI


sigmax = 2*micro
sigmay = 1*micro
sigmaz = 3.*micro
emitx = 100*milli
emity = 20*milli
gamma = 125*1e9*eV/(m_e*c**2)
x_m = 10*sigmax
y_m = 10*sigmay
z_m = 10*sigmaz
focal_distance = 4*sigmaz
theta = 0.25*np.pi

def s(z, sigma, emit, center):
    '''The theoretical size of a focusing beam (in the absence of space charge),
    at position z, given its emittance and size at focus.'''
    return np.sqrt(sigma**2 + emit**2 * (z - z_m - focal_distance)**2 / (sigma**2)) + center

filename = sys.argv[1]

series = OpenPMDTimeSeries('./diags/openpmd/')

# Rotated beam
xrot, yrot, zrot, w, uxrot, uyrot, uzrot = series.get_particle( ['x', 'y', 'z', 'w', 'ux', 'uy', 'uz'], species='beam1', iteration=0, plot=False)

# Rotate the beam back so that it propagates along z
z = z_m + np.cos(-theta)*(zrot-z_m) - np.sin(-theta)*(xrot-x_m)
x = x_m + np.sin(-theta)*(zrot-z_m) + np.cos(-theta)*(xrot-x_m)
y = yrot

# Compute the size of the beam in different z slices
gridz = np.linspace(z_m-focal_distance*0.8, z_m+focal_distance*0.8, 64)
tol = 0.5*(gridz[1]-gridz[0])
sx, sy = [], []
for z_g in gridz:
    i = np.abs(z - z_g) < tol
    if (np.sum(i)!=0):
        mux = np.average(x[i], weights=w[i])
        muy = np.average(y[i], weights=w[i])
        sx.append(x_m+np.sqrt(np.average((x[i]-mux)**2, weights=w[i])))
        sy.append(y_m+np.sqrt(np.average((y[i]-muy)**2, weights=w[i])))
    else:
        sx.append(0.)
        sy.append(0.)

# Theoretical prediction for the size of the beam in each z slice
sx_theory = s(gridz, sigmax, emitx/gamma, x_m)
sy_theory = s(gridz, sigmay, emity/gamma, y_m)

assert(np.allclose(sx, sx_theory, rtol=0.073, atol=0))
assert(np.allclose(sy, sy_theory, rtol=0.081, atol=0))

# Rotate the beam back so that it propagates along z
uz = np.cos(-theta)*uzrot - np.sin(-theta)*uxrot
ux = np.sin(-theta)*uzrot + np.cos(-theta)*uxrot
uy = uyrot

uz_m = np.average(uz, weights=w)
uz_th = np.sqrt(np.average((uz-uz_m)**2, weights=w))

ux_m = np.average(ux, weights=w)
ux_th = np.sqrt(np.average((ux-ux_m)**2, weights=w))

err1 = np.abs( (uz_m - gamma) / gamma)
err2 = np.abs( uz_th )
err3 = np.abs( ux_m )
err4 = np.abs( (ux_th - emitx/sigmax) / (emitx/sigmax))

print(ux_m, ux_th)
print(err3, err4)

print(np.isclose(uz_m, gamma, rtol=1e-8))
print(np.isclose(uz_th,0., rtol=1e-8))

print(np.isclose(ux_m, 0.))


print(np.isclose(ux_th, emitx / sigmax, rtol=1e-3))

#errx =  np.average( np.asarray(ux) )
#errz =  np.abs( (uz - gamma) / gamma )
#print(np.max(errx), np.max(errz))
#assert(np.allclose(uz, uz_m, rtol=0.081, atol=0))


#import matplotlib.pyplot as plt
#N = 100
#plt.quiver(zrot[::N], xrot[::N], uzrot[::N], uxrot[::N])
#plt.quiver(z[::N], x[::N], uz[::N], ux[::N])
#plt.show()

#test_name = os.path.split(os.getcwd())[1]
#checksumAPI.evaluate_checksum(test_name, filename)
