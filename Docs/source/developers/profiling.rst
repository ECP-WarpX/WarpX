.. _developers-profiling:

Profiling
=========

.. _developers-profiling-tiny-profiler:

Tiny Profiler
-------------

The so-called "Tiny Profiler" is a component of lightweight profiling that we perform with each WarpX run by default.
A table of the summary of most often called and most costly functions is printed to ``stdout`` at the end of every WarpX run.

`See the full documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools_Chapter.html>`__.

There is a script located `here <https://github.com/AMReX-Codes/amrex/tree/development/Tools/TinyProfileParser>`__ that parses the Tiny Profiler output and generates a JSON file that can be used with `Hatchet <https://hatchet.readthedocs.io/en/latest/>`__ in order to analyze performance.


.. _developers-profiling-nsight-systems:

Nvidia Nsight-Systems
---------------------

`Vendor homepage <https://developer.nvidia.com/nsight-systems>`__ and `product manual <https://docs.nvidia.com/nsight-systems/>`__.

Example on how to create traces on a multi-GPU system that uses the Slurm scheduler:

.. code-block:: bash

   rm -rf profiling*

   # adjust if needed
   #export GPUS_PER_SOCKET=2
   #export GPUS_PER_NODE=4
   export OMP_NUM_THREADS=1

   # record
   srun --ntasks=4 --gpus=4 --cpu-bind=cores \
       nsys profile -f true \
         -o profiling_%q{SLURM_TASK_PID} \
         -t mpi,cuda,nvtx,osrt,openmp   \
         --mpi-impl=openmpi \
       ./warpx.3d.MPI.CUDA.DP.QED inputs_3d \
         warpx.numprocs=1 1 4 amr.n_cell=512 512 2048 max_step=10

In this example, the individual lines for recording a trace profile are:

* ``srun``: execute multi-GPU runs with ``srun`` (Slurm's ``mpiexec`` wrapper), here for four GPUs
* ``-f true`` overwrite previously written trace profiles
* ``-o``: record one profile file per MPI rank (per GPU); if you run ``mpiexec``/``mpirun`` with OpenMPI directly, replace ``SLURM_TASK_PID`` with ``OMPI_COMM_WORLD_RANK``
* ``-t``: select a couple of APIs to trace
* ``--mpi--impl``: optional, hint the MPI flavor
* ``./warpx...``: select the WarpX executable and a good inputs file
* ``warpx.numprocs=...``: make the run short, reasonably small, and run only a few steps

Now open the created trace files (per rank) in the Nsight-Systems GUI.
This can be done on another system than the one that recorded the traces.
For example, if you record on a cluster and open the analysis GUI on your laptop, it is recommended to make sure that versions of Nsight-Systems match on the remote and local system.
