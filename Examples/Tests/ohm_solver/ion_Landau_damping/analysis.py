#!/usr/bin/env python3
#
# --- Analysis script for the hybrid-PIC example of ion Landau damping.

import dill
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from pywarpx import picmi

constants = picmi.constants

matplotlib.rcParams.update({'font.size': 20})

# load simulation parameters
with open(f'sim_parameters.dpkl', 'rb') as f:
    sim = dill.load(f)

if sim.test:
    import os
    import sys
    sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
    import checksumAPI

    # this will be the name of the plot file
    fn = sys.argv[1]
    test_name = os.path.split(os.getcwd())[1]
    checksumAPI.evaluate_checksum(test_name, fn)

data = np.loadtxt("diags/field_data.txt", skiprows=1)
field_idx_dict = {'z': 2, 'Ez': 3}

step = data[:,0]

num_steps = len(np.unique(step))

# get the spatial resolution
resolution = len(np.where(step == 0)[0]) - 1

# reshape to separate spatial and time coordinates
sim_data = data.reshape((num_steps, resolution+1, data.shape[1]))

z_grid = sim_data[1, :, field_idx_dict['z']]
idx = np.argsort(z_grid)[1:]
dz = np.mean(np.diff(z_grid[idx]))
dt = np.mean(np.diff(sim_data[:,0,1]))

data = np.zeros((num_steps, resolution))
for i in range(num_steps):
    data[i,:] = sim_data[i,idx,field_idx_dict['Ez']]

print(f"Data file contains {num_steps} time snapshots.")
print(f"Spatial resolution is {resolution}")

field_kt = np.fft.fft(data[:, :], axis=1)

t_norm = 2.0 * np.pi * sim.m / sim.Lz * sim.v_ti

# Plot the 4th Fourier mode
fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

ax1.plot(np.arange(num_steps)*dt*t_norm, np.abs(field_kt[:, sim.m] / field_kt[0, sim.m]), 'r', label=f'$T_i/T_e$ = {sim.T_ratio:.2f}')

ax1.grid()
ax1.legend()
ax1.set_yscale('log')
ax1.set_ylabel('$|E_z|/E_0$')
ax1.set_xlabel('t $(k_mv_{th,i})$')
ax1.set_xlim(0, 18)

ax1.set_title(f"Ion Landau damping - {sim.dim}d")
plt.tight_layout()
plt.savefig(f"diags/ion_Landau_damping_T_ratio_{sim.T_ratio}.png")
