#!/usr/bin/env python

import os, sys, glob, pickle
import numpy as np
import matplotlib.pyplot as plt

color_mode = 'sim_id' # 'batch_id'

print("Reading latest npy and pickle files in " + sys.argv[1])
list_of_files = glob.glob(sys.argv[1] + '*npy') # * means all if need specific format then *.csv
latest_file = max(list_of_files, key=os.path.getctime)
results = np.load(latest_file)
list_of_files = glob.glob(sys.argv[1] + '*pickle') # * means all if need specific format then *.csv
latest_file = max(list_of_files, key=os.path.getctime)
with open(latest_file, 'rb') as f:
    data = pickle.load(f)
if color_mode == 'batch_id':
    nbatches = len(data[1]['run_order'])

# Remove un-returned (i.e., not submitted) simulations
results = results[ results['returned'] == 1 ]

names = results.dtype.names
results_dict = {}
for i in range(len(names)):
    results_dict[names[i]] = np.array([task[i] for task in results])
    
plot_output = ['f', 'emittance', 'charge', 'energy_avg', 'energy_std']
plot_input = ['ramp_down_1', 'zlens_1', 'adjust_factor', 'ramp_down_2']

d_ylim = {
    'f': (3.e-7, 5.e-3),
    'emittance': (3.e-7, 5.e-3),
    'charge': (1.e-7, 1.e-5),
    'energy_avg': (5.e1, 5.e2),
    'energy_std': (1.e-2,1.e-1),
}

d_yunits = {
    'f': " (a.u.)",
    'emittance': " (m rad)",
    'charge': " (pC/m)",
    'energy_avg': " (MeV)",
    'energy_std': " (dimless)",
}

if np.min( results_dict['charge'] == 0.0 ):
    sys.exit()
if np.min( results_dict['f'] == 0.0 ):
    sys.exit()

# Print properties for optimum
ind_best = np.argmin( results_dict['f'] )
# for key in plot_input + plot_output + ['local_pt']:
for key in plot_input + plot_output:
    print(key, results_dict[key][ind_best])
print("charge_f/charge_i ", results[ind_best]['charge']/np.max(results_dict['charge']))

plt.figure(figsize=(16,10))
plt.style.use('dark_background')
for icount, iname in enumerate(plot_input):
    for count, name in enumerate(plot_output):
        plt.subplot(4,5,5*icount + count+1)
        if color_mode == 'batch_id':
            for batch in data[1]['run_order']:
                my_sims = np.isin(results_dict['sim_id'], data[1]['run_order'][batch])
                plt.scatter(np.squeeze(results_dict[iname][my_sims]),
                            np.squeeze(results_dict[ name ][my_sims]),
                            s=2)
                # c=[1.-batch/nbatches,0,batch/nbatches], s=1)
        elif color_mode == 'sim_id':
            plt.scatter(np.squeeze(results_dict[iname]), np.squeeze(results_dict[ name ]), c=results_dict['sim_id'], s=2, cmap='viridis')
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('start time (arb. units)')
        plt.yscale('log')
        plt.xlim(np.min(results_dict[iname]), np.max(results_dict[iname]))
        plt.ylim( d_ylim[name])
        plt.grid()
        plt.title( name )
        plt.xlabel( iname )
        plt.ylabel( name + d_yunits[name] )

plt.tight_layout()

plt.savefig('results.pdf', bbox_inches='tight')
plt.savefig('results.jpg', bbox_inches='tight')
