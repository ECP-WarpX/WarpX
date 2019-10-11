#! /usr/bin/env python
"""
This script checks the space-charge initialization routine, by
verifying that the space-charge field of a Gaussian beam corresponds to
the expected theoretical field.
"""
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scc
yt.funcs.mylog.setLevel(0)

# Parameters from the Simulation
Qtot = -1.e-20
r0 = 2.e-6

# Open data file
filename = sys.argv[1]
ds = yt.load( filename )
# Extract data
ad0 = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
Ex_array = ad0['Ex'].to_ndarray().squeeze()
Ez_array = ad0['Ez'].to_ndarray().squeeze()

if ds.dimensionality == 2:
    # Compute theoretical field
    Nx, Nz, _ =  ds.domain_dimensions
    xmin, zmin, _ = ds.domain_left_edge.v
    Lx, Lz, _ = ds.domain_width.v
    x = xmin + Lx/Nx*(0.5+np.arange(Nx))
    z = zmin + Lz/Nz*(0.5+np.arange(Nz))
    x_2d, z_2d = np.meshgrid(x, z, indexing='ij')
    factor = (Qtot/r0)/(2*np.pi*scc.epsilon_0*(x_2d**2+z_2d**2)) * \
                (1-np.exp(-(x_2d**2+z_2d**2)/(2*r0**2)))
    Ex_th = x_2d * factor
    Ez_th = z_2d * factor

# Plot theory and data
plt.figure(figsize=(10,10))
plt.subplot(221)
plt.title('Ex: Theory')
plt.imshow(Ex_th)
plt.subplot(222)
plt.title('Ex: Simulation')
plt.imshow(Ex_array)
plt.subplot(223)
plt.title('Ez: Theory')
plt.imshow(Ez_th)
plt.subplot(224)
plt.title('Ez: Simulation')
plt.imshow(Ez_array)
plt.savefig('Comparison.png')

# Automatically check the results
assert np.allclose( Ex_array, Ex_th, atol=0.1*Ex_th.max() )
assert np.allclose( Ez_array, Ez_th, atol=0.1*Ez_th.max() )
