#!/usr/bin/env python
"""
Execute locally via the following command:
   python run_libE_on_warpX.py --comms local --nworkers 3

The number of concurrent evaluations of the objective function will be 4-2=2
as one MPI rank for the manager and one MPI rank for the persistent gen_f.
"""

generator_type = 'aposmm' # 'random' # 'aposmm'
machine = 'local' # 'summit'

import numpy as np
from warpX_simf import run_warpX # Sim function from current directory

# Import libEnsemble modules
from libensemble.libE import libE
if generator_type == 'random':
    from libensemble.gen_funcs.sampling import uniform_random_sample as gen_f
    from libensemble.alloc_funcs.give_sim_work_first import give_sim_work_first as alloc_f
elif generator_type == 'aposmm':
    import libensemble.gen_funcs
    libensemble.gen_funcs.rc.aposmm_optimizer = 'nlopt'
    from libensemble.gen_funcs.persistent_aposmm import aposmm as gen_f
    from libensemble.alloc_funcs.persistent_aposmm_alloc import persistent_aposmm_alloc as alloc_f
else:
    print("you shouldn' hit that") ; sys.exit()

from libensemble.tools import parse_args, save_libE_output, add_unique_random_streams
from libensemble import libE_logger
from libensemble.executors.mpi_executor import MPIExecutor

# Import machine-specific run parameters
if machine == 'local':
    from Local import machine_specs
elif machine == 'aposmm':
    from Summit import machine_specs
else:
    print("you shouldn' hit that") ; sys.exit()

libE_logger.set_level('INFO')

nworkers, is_master, libE_specs, _ = parse_args()

# Set to full path of warp executable
sim_app = machine_specs['sim_app']

# Problem dimension. This is the number of input parameters exposed,
# that LibEnsemble will vary in order to minimize a single output parameter.
n = 4

exctr = MPIExecutor(central_mode=True)
exctr.register_calc(full_path=sim_app, calc_type='sim')

# State the objective function, its arguments, output, and necessary parameters (and their sizes)
# Here, the 'user' field is for the user's (in this case, the simulation) convenience.
# Feel free to use it to pass number of nodes, number of ranks per note, time limit per simulation etc.
sim_specs = {'sim_f': run_warpX, # Function whose output is being minimized. The parallel WarpX run is launched from run_WarpX.
             'in': ['x'], # Name of input for sim_f, that LibEnsemble is allowed to modify. May be a 1D array.
             'out': [('f', float), # f is the single float output that LibEnsemble minimizes.
                     ('energy_std', float, (1,)), # All parameters below are not used for calculation, just output for convenience. Final relative energy spread.
                     ('energy_avg', float, (1,)), # Final average energy, in MeV.
                     ('charge', float, (1,)), # Final beam charge.
                     ('emittance', float, (1,)), # Final beam emittance.
                     ('ramp_down_1', float, (1,)), # input parameter: length of first downramp.
                     ('ramp_down_2', float, (1,)), # input parameter: Length of second downramp.
                     ('zlens_1', float, (1,)), # input parameter: position of the focusing lens.
                     ('adjust_factor', float, (1,)), # Relative stength of the lens (1. is from back-of-the-envelope calculation)
             ],
             'user': {'input_filename': 'inputs', # name of input file
                      'sim_kill_minutes': 3, # Run timeouts after 3 mins
                      'machine_specs': machine_specs # machine-specific parameters
             }
}

# State the generating function, its arguments, output, and necessary parameters.
if generator_type == 'random':
    # Here, the 'user' field is for the user's (in this case, the RNG) convenience.
    gen_specs = {'gen_f': gen_f, # Generator function. Will randomly generate new sim inputs 'x'.
                 'in': [], # Generator input. This is a RNG, no need for inputs.
                 'out': [('x', float, (n,))], # parameters to input into the simulation.
                 'user': {'gen_batch_size': nworkers, # Total max number of sims running concurrently.
                          'lb': np.array([2.e-3, 2.e-3, 0.005, .1]),  # Lower bound for the n parameters.
                          'ub': np.array([2.e-2, 2.e-2, 0.028, 3.]), # Upper bound for the n parameters.
                 }
    }

    alloc_specs = {'alloc_f': alloc_f, # Allocator function, decides what a worker should do. We use a LibEnsemble allocator.
                   'out': [('allocated', bool)],
                   'user': {'batch_mode': True, # If true wait for all sims to process before generate more
                            'num_active_gens': 1 # Only one active generator at a time
                   }
    }

elif generator_type == 'aposmm':
    # Here, the 'user' field is for the user's (in this case, the optimizer) convenience.
    gen_specs = {'gen_f': gen_f, # Generator function. Will randomly generate new sim inputs 'x'.
                 'in': [], # Generator input. This is a RNG, no need for inputs.
                 'out': [('x', float, (n,)), # parameters to input into the simulation.
                         ('x_on_cube', float, (n,)), # x scaled to a unique cube.
                         ('sim_id', int), # unique ID of simulation.
                         ('local_min', bool), # Whether this point is a local minimum.
                         ('local_pt', bool) # whether the point is from a local optimization run or a random sample point.
                 ],
                 'user': {'initial_sample_size': max(2*(nworkers-2), 1), # Number of sims for initial random sampling. Optimizer starts afterwards.
                          'localopt_method': 'LN_BOBYQA', # APOSMM/NLOPT optimization method
                          'num_pts_first_pass': nworkers, # 
                          'xtol_rel': 1e-3, # Relative tolerance of inputs
                          'ftol_abs': 3e-8, # Absolute tolerance of output 'f'. Determines when local optimization stops.
                          'lb': np.array([2.e-3, 2.e-3, 0.005, .1]),  # Lower bound for the n input parameters.
                          'ub': np.array([2.e-2, 2.e-2, 0.028, 3.]), # Upper bound for the n input parameters.
                 }
    }

    alloc_specs = {'alloc_f': alloc_f, # Allocator function, decides what a worker should do. We use a LibEnsemble allocator.
                   'out': [('given_back', bool)],
                   'user': {}}

else:
    print("you shouldn' hit that") ; sys.exit()
    
libE_specs['save_every_k_sims'] = 100 # Save H to file every N simulation evaluations
libE_specs['sim_input_dir'] = 'sim' # Sim directory to be copied for each worker

sim_max = 3 # Maximum number of simulations
exit_criteria = {'sim_max': sim_max} # Exit after running sim_max simulations

# Create a different random number stream for each worker and the manager
persis_info = add_unique_random_streams({}, nworkers + 1)

# Run LibEnsemble, and store results in history array H
H, persis_info, flag = libE(sim_specs, gen_specs, exit_criteria,
                            persis_info, alloc_specs, libE_specs)

# Save results to numpy file
if is_master:
    save_libE_output(H, persis_info, __file__, nworkers)
