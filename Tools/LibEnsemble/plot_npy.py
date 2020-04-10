#!/usr/bin/env python

"""
Plotting script for LibEnsemble run on WarpX simulations.
input: A directory which which the script will read the
latest .npy and .pickle files generated by LibEnsemble
"""

import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt

generator_type = 'aposmm'  # or 'random'

print("Reading latest .npy and .pickle files in " + sys.argv[1])

# Read latest .npy file, and store the result in results
# * means all if need specific format then *.csv
list_of_files = glob.glob(sys.argv[1] + '*npy')
latest_file = max(list_of_files, key=os.path.getctime)
results = np.load(latest_file)

# Read latest .pickle file, and store the result in pickle_data
# * means all if need specific format then *.csv
list_of_files = glob.glob(sys.argv[1] + '*pickle')
latest_file = max(list_of_files, key=os.path.getctime)
with open(latest_file, 'rb') as f:
    pickle_data = pickle.load(f)
if generator_type == 'aposmm':
    nbatches = len(pickle_data[1]['run_order'])

# Remove un-returned (i.e., not submitted) simulations
results = results[results['returned'] == 1]

# Re-organize results into a dictionary. Each key is a field in libE_output,
# each value is an array with the value of this field for all WarpX simulations
names = results.dtype.names
results_dict = {}
for i in range(len(names)):
    results_dict[names[i]] = np.array([task[i] for task in results])

# Sanity checks: exit if a simulation has null charge, as the script
# doesn't handle this (shouldn't be too hard to add though).
if np.min(results_dict['charge'] == 0.0):
    sys.exit()
if np.min(results_dict['f'] == 0.0):
    sys.exit()

# Custom: which output parameters are plotted
# as a function of which input parameters
plot_output = ['f', 'emittance', 'charge', 'energy_avg', 'energy_std']
plot_input = ['ramp_down_1', 'zlens_1', 'adjust_factor', 'ramp_down_2']

# For convenience, dictionary with pyplot.ylim
# tuples for all elements in plot_output
d_ylim = {
    'f': (3.e-7, 5.e-3),
    'emittance': (3.e-7, 5.e-3),
    'charge': (1.e-7, 1.e-5),
    'energy_avg': (5.e1, 5.e2),
    'energy_std': (1.e-2, 1.e-1)
}

# For convenience, dictionary with units for all elements in plot_output
d_yunits = {
    'f': " (a.u.)",
    'emittance': " (m rad)",
    'charge': " (pC/m)",
    'energy_avg': " (MeV)",
    'energy_std': " (dimless)"
}

# Print input and output parameters for the optimal run
ind_best = np.argmin(results_dict['f'])
for key in plot_input + plot_output:
    print(key, results_dict[key][ind_best])
print("charge_f/charge_i ",
      results[ind_best]['charge']/np.max(results_dict['charge']))

# And now plotting begins
plt.figure(figsize=(16, 10))
plt.style.use('dark_background')
# Loop over all input parameters in plot_input
for icount, iname in enumerate(plot_input):
    # Loop over all output parameters in plot_output
    for count, name in enumerate(plot_output):
        plt.subplot(len(plot_input), len(plot_output),
                    len(plot_output)*icount + count+1)
        if generator_type == 'aposmm':
            # Plot 1 point per run. Runs that are part of the same local
            # minimum search are batched together and plotted with the
            # same color
            for batch in pickle_data[1]['run_order']:
                my_sims = np.isin(results_dict['sim_id'],
                                  pickle_data[1]['run_order'][batch])
                plt.scatter(np.squeeze(results_dict[iname][my_sims]),
                            np.squeeze(results_dict[name][my_sims]),
                            s=2)
        elif generator_type == 'random':
            # Plot 1 point per run. Colorbar stands for sim ID
            # (larger value = later run)
            plt.scatter(np.squeeze(results_dict[iname]),
                        np.squeeze(results_dict[name]),
                        c=results_dict['sim_id'], s=2, cmap='viridis')
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('start time (arb. units)')
        else:
            print("generator_type should be either random or aposmm")
            sys.exit()
        plt.yscale('log')
        plt.xlim(np.min(results_dict[iname]), np.max(results_dict[iname]))
        plt.ylim(d_ylim[name])
        plt.grid()
        plt.title(name)
        plt.xlabel(iname)
        plt.ylabel(name + d_yunits[name])

plt.tight_layout()

# Save figure
plt.savefig('results.pdf', bbox_inches='tight')
plt.savefig('results.jpg', bbox_inches='tight')
