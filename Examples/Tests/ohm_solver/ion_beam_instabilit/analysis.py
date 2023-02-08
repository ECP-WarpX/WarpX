#!/usr/bin/env python3
#
# --- Analysis script for the hybrid-PIC example producing EM modes.

import dill
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from pywarpx import picmi

constants = picmi.constants

matplotlib.rcParams.update({'font.size': 20})

# load simulation parameters
with open(f'sim_parameters.dpkl', 'rb') as f:
    sim = dill.load(f)

if sim.resonant:
    resonant_str = 'resonant'
else:
    resonant_str = 'non resonant'

data = np.loadtxt("diags/field_data.txt", skiprows=1)
if sim.dim == 1:
    field_idx_dict = {'z': 4, 'Ez': 7, 'Bx': 8, 'By': 9}
else:
    field_idx_dict = {'z': 2, 'Ez': 3, 'Bx': 4, 'By': 5}

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

data = np.zeros((num_steps, resolution, 3))
for i in range(num_steps):
    data[i,:,0] = sim_data[i,idx,field_idx_dict['Bx']]
    data[i,:,1] = sim_data[i,idx,field_idx_dict['By']]
    data[i,:,2] = sim_data[i,idx,field_idx_dict['Ez']]

print(f"Data file contains {num_steps} time snapshots.")
print(f"Spatial resolution is {resolution}")

# Create the stack time plot
fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

extent = [0, sim.Lz/sim.l_i, 0, num_steps*dt*sim.w_ci] # num_steps*dt/sim.t_ci]
im = ax1.imshow(
    data[:,:,1]/sim.B0, extent=extent, origin='lower',
    cmap='seismic'#, vmin=vmin, vmax=vmax, aspect="equal",

)

# Colorbar
fig.subplots_adjust(right=0.825)
cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='$B_y/B_0$')

ax1.set_xlabel("$x/l_i$")
ax1.set_ylabel("$t \Omega_i$ (rad)")

ax1.set_title(f"Ion beam R instability - {resonant_str} case")
plt.savefig(f"diags/ion_beam_R_instability_{resonant_str}_eta_{sim.eta}_substeps_{sim.substeps}.png")
plt.close()

# Plot the 4th, 5th and 6th Fourier modes
field_kt = np.fft.fft(data[:, :, 1], axis=1)

k_norm = 1.0 / sim.l_i

k = 2*np.pi * np.fft.fftfreq(resolution, dz) / k_norm

print(field_kt.shape)

plt.plot(np.arange(num_steps)*dt*sim.w_ci, np.abs(field_kt[:, 3] / sim.B0), 'r', label=f'm = 4, $kl_i={k[3]:.2f}$')
plt.plot(np.arange(num_steps)*dt*sim.w_ci, np.abs(field_kt[:, 4] / sim.B0), 'b', label=f'm = 5, $kl_i={k[4]:.2f}$')
plt.plot(np.arange(num_steps)*dt*sim.w_ci, np.abs(field_kt[:, 5] / sim.B0), 'k', label=f'm = 6, $kl_i={k[5]:.2f}$')

plt.grid()
plt.legend()
plt.yscale('log')
plt.ylabel('$|B_y/B_0|$')
plt.xlabel('$t\Omega_i$ (rad)')
plt.tight_layout()
plt.savefig(f"diags/ion_beam_R_instability_{resonant_str}_eta_{sim.eta}_substeps_{sim.substeps}_low_modes.png")
plt.close()

with h5py.File('diags/Python_hybrid_PIC_plt/openpmd_004000.h5', 'r') as data:

    timestep = str(np.squeeze([key for key in data['data'].keys()]))

    z = np.array(data['data'][timestep]['particles']['ions']['position']['z'])
    vy = np.array(data['data'][timestep]['particles']['ions']['momentum']['y'])
    w = np.array(data['data'][timestep]['particles']['ions']['weighting'])

fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

im = ax1.hist2d(
    z/sim.l_i, vy/sim.M/sim.vA, weights=w, density=True,
    range=[[0, 250], [-10, 10]], bins=250, cmin=1e-5
)

# Colorbar
fig.subplots_adjust(bottom=0.15, right=0.815)
cbar_ax = fig.add_axes([0.83, 0.2, 0.03, 0.6])
fig.colorbar(im[3], cax=cbar_ax, orientation='vertical', format='%.0e', label='$f(z, v_y)$')

ax1.set_xlabel("$x/l_i$")
ax1.set_ylabel("$v_{y}/v_A$")

ax1.set_title(f"Ion beam R instability - {resonant_str} case")
plt.savefig(f"diags/ion_beam_R_instability_{resonant_str}_eta_{sim.eta}_substeps_{sim.substeps}_core_phase_space.png")
plt.close()

with h5py.File('diags/Python_hybrid_PIC_plt/openpmd_004000.h5', 'r') as data:

    timestep = str(np.squeeze([key for key in data['data'].keys()]))

    z = np.array(data['data'][timestep]['particles']['beam_ions']['position']['z'])
    vy = np.array(data['data'][timestep]['particles']['beam_ions']['momentum']['y'])
    w = np.array(data['data'][timestep]['particles']['beam_ions']['weighting'])

fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

im = ax1.hist2d(
    z/sim.l_i, vy/sim.M/sim.vA, weights=w, density=True,
    range=[[0, 250], [-10, 10]], bins=250, cmin=1e-5
)

# Colorbar
fig.subplots_adjust(bottom=0.15, right=0.815)
cbar_ax = fig.add_axes([0.83, 0.2, 0.03, 0.6])
fig.colorbar(im[3], cax=cbar_ax, orientation='vertical', format='%.0e', label='$f(z, v_y)$')

ax1.set_xlabel("$x/l_i$")
ax1.set_ylabel("$v_{y}/v_A$")

ax1.set_title(f"Ion beam R instability - {resonant_str} case")
plt.savefig(f"diags/ion_beam_R_instability_{resonant_str}_eta_{sim.eta}_substeps_{sim.substeps}_beam_phase_space.png")
plt.show()