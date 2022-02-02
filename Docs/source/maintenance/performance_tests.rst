.. _developers-performance_tests:

Automated performance tests
===========================

WarpX has automated performance test scripts, which run weak scalings for various tests on a weekly basis. The results are stored in the `perf_logs repo <https://github.com/ECP-WarpX/perf_logs>`_ and plots of the performance history can be found on `this page <https://ecp-warpx.github.io/perf_logs/>`_.

These performance tests run automatically, so they need to do ``git`` operations etc. For this reason, they need a separate clone of the source repos, so they don't conflict with one's usual operations. This is typically in a sub-directory in the ``$HOME``, with variable ``$AUTOMATED_PERF_TESTS`` pointing to it. Similarly, a directory is needed to run the simulations and store the results. By default, it is ``$SCRATCH/performance_warpx``.

The test runs a weak scaling (1,2,8,64,256,512 nodes) for 6 different tests ``Tools/PerformanceTests/automated_test_{1,2,3,4,5,6}_*``, gathered in 1 batch job per number of nodes to avoid submitting too many jobs.

Setup on Summit @ OLCF
----------------------

Here is an example setup for Summit:

.. code-block:: sh

   # I put the next three lines in $HOME/my_bashrc.sh
   export proj=aph114  # project for job submission
   export AUTOMATED_PERF_TESTS=$HOME/AUTOMATED_PERF_TESTS/
   export SCRATCH=/gpfs/alpine/scratch/$(whoami)/$proj/

   mkdir $HOME/AUTOMATED_PERF_TESTS
   cd $AUTOMATED_PERF_TESTS
   git clone https://github.com/ECP-WarpX/WarpX.git warpx
   git clone https://github.com/ECP-WarpX/picsar.git
   git clone https://github.com/AMReX-Codes/amrex.git
   git clone https://github.com/ECP-WarpX/perf_logs.git

Then, in ``$AUTOMATED_PERF_TESTS``, create a file ``run_automated_performance_tests_512.sh`` with the following content:

.. code-block:: sh

   #!/bin/bash -l
   #BSUB -P APH114
   #BSUB -W 00:15
   #BSUB -nnodes 1
   #BSUB -J PERFTEST
   #BSUB -e err_automated_tests.txt
   #BSUB -o out_automated_tests.txt

   module load nano
   module load cmake/3.20.2
   module load gcc/9.3.0
   module load cuda/11.0.3
   module load blaspp/2021.04.01
   module load lapackpp/2021.04.00
   module load boost/1.76.0
   module load adios2/2.7.1
   module load hdf5/1.10.7

   module unload darshan-runtime

   export AMREX_CUDA_ARCH=7.0
   export CC=$(which gcc)
   export CXX=$(which g++)
   export FC=$(which gfortran)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)

   # Make sure all dependencies are installed and loaded
   cd $HOME
   module load python/3.8.10
   module load freetype/2.10.4     # matplotlib
   module load openblas/0.3.5-omp
   export BLAS=$OLCF_OPENBLAS_ROOT/lib/libopenblas.so
   export LAPACK=$OLCF_OPENBLAS_ROOT/lib/libopenblas.so
   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv
   python3 -m venv $HOME/sw/venvs/warpx-perftest
   source $HOME/sw/venvs/warpx-perftest/bin/activate
   # While setting up the performance tests for the first time,
   # execute the lines above this comment and then the commented
   # lines below this comment once, before submission.
   # The commented lines take too long for the job script.
   #python3 -m pip install --upgrade pip
   #python3 -m pip install --upgrade cython
   #python3 -m pip install --upgrade numpy
   #python3 -m pip install --upgrade markupsafe
   #python3 -m pip install --upgrade pandas
   #python3 -m pip install --upgrade matplotlib==3.2.2  # does not try to build freetype itself
   #python3 -m pip install --upgrade bokeh
   #python3 -m pip install --upgrade gitpython
   #python3 -m pip install --upgrade tables

   # Run the performance test suite
   cd $AUTOMATED_PERF_TESTS/warpx/Tools/PerformanceTests/
   python run_automated.py --n_node_list='1,2,8,64,256,512' --automated

   # submit next week's job
   cd $AUTOMATED_PERF_TESTS/
   next_date=`date -d "+7 days" '+%Y:%m:%d:%H:%M'`
   bsub -b $next_date ./run_automated_performance_tests_512.sh

Then, running

.. code-block:: sh

   bsub run_automated_performance_tests_512.sh

will submit this job once, and all the following ones. It will:

 - Create directory ``$SCRATCH/performance_warpx`` if doesn't exist.
 - Create 1 sub-directory per week per number of nodes (1,2,8,64,256,512).
 - Submit one job per number of nodes. It will run 6 different tests, each twice (to detect fluctuations).
 - Submit an analysis job, that will read the results ONLY AFTER all runs are finished. This uses the dependency feature of the batch system.
 - This job reads the Tiny Profiler output for each run, and stores the results in a pandas file at the hdf5 format.
 - Execute ``write_csv.py`` from the ``perf_logs`` repo to append a csv and a hdf5 file with the new results.
 - Commit the results (but DO NOT PUSH YET)

Then, the user periodically has to

.. code-block:: sh

   cd $AUTOMATED_PERF_TESTS/perf_logs
   git pull # to get updates from someone else, or from another supercomputer
   git push

This will update the database but not the online plots. For this, you need to periodically run something like

.. code-block:: sh

   cd $AUTOMATED_PERF_TESTS/perf_logs
   git pull
   python generate_index_html.py
   git add -u
   git commit -m "upload new html page"
   git push

Setup on Cori @ NERSC
---------------------

Still to be written!
