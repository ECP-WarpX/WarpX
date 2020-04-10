import os

machine_specs = \
    {
        'name': 'summit',
        'sim_app': os.environ['HOME'] + '/warpx/Bin/main2d.gnu.TPROF.MPI.CUDA.ex',
        'extra_args': '-n 1 -a 1 -g 1 -c 1 --bind=packed:1 --smpiargs="-gpu"',
        'OMP_NUM_THREADS': '1'
    }
