import os

machine_specs = \
    {
        'name': 'local',
        'cores': 1,
        'sim_app': os.environ['HOME'] + '/warpx/Bin/main2d.gnu.TPROF.MPI.OMP.ex',
        'extra_args': '',
        'OMP_NUM_THREADS': '2',
    }
