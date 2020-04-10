import os

machine_specs = {
    'name': 'local',  # Machine name
    'cores': 1,  # Number of cores per simulation
    'sim_app': os.environ['HOME'] + '/warpx/Bin/main2d.gnu.TPROF.MPI.OMP.ex',
    'extra_args': '',  # extra arguments passed to mpirun/mpiexec at execution
    'OMP_NUM_THREADS': '2'
}
