.. _developers-profiling:

Profiling the code
==================

Profiling allows us to find the bottle-necks of the code as it is currently implemented.
Bottle-necks are the parts of the code that may delay the simulation, making it more computationally expensive.
Once found, we can update the related code sections and improve its efficiency.
Profiling tools can also be used to check how load balanced the simulation is, i.e. if the work is well distributed across all MPI ranks used.
Load balancing can be activated in WarpX by setting input parameters, see the :ref:`parallelization input parameter section <running-cpp-parameters-parallelization>`.

.. _developers-profiling-tiny-profiler:

Profiling with AMReX's built-in profiling tools
-----------------------------------------------

By default, WarpX uses the AMReX baseline tool, the TINYPROFILER, to evaluate the time information for different parts of the code (functions) between the different MPI ranks.
The results, timers, are stored into four tables in the standard output, stdout, that are located below the simulation steps information and above the warnings regarding unused input file parameters (if there were any).

The timers are displayed in tables for which the columns correspond to:

* name of the function
* number of times it is called in total
* minimum of time spent exclusively/inclusively in it, between all ranks
* average of time, between all ranks
* maximum time, between all ranks
* maximum percentage of time spent, across all ranks

If the simulation is well load balanced the minimum, average and maximum times should be identical.

The top two tables refer to the complete simulation information.
The bottom two are related to the Evolve() section of the code (where each time step is computed).

Each set of two timers show the exclusive, top, and inclusive, bottom, information depending on whether the time spent in nested sections of the codes are included.

For more detailed information please visit the `AMReX profiling documentation <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools_Chapter.html>`__.

Please note that the `AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`__ for ``AMReX_TINY_PROFILE`` (our default: ``ON``) and full profiling traces via ``AMReX_BASE_PROFILE`` are mutually exclusive.
Further tracing options are sub-options of ``AMReX_BASE_PROFILE``.
If you want to switch from tiny profiling defaults to trace profiling, configure like this:

.. code-block:: bash

   cmake -S . -B build -DAMReX_BASE_PROFILE=YES -DAMReX_TRACE_PROFILE=YES  -DAMReX_COMM_PROFILE=YES -DAMReX_TINY_PROFILE=OFF

To turn on the tiny profiler again, remove the ``build`` directory or turn off ``AMReX_BASE_PROFILE`` again:

.. code-block:: bash

   cmake -S . -B build -DAMReX_BASE_PROFILE=OFF -DAMReX_TINY_PROFILE=ON

.. note::

   When creating performance-related issues on the WarpX GitHub repo, please include Tiny Profiler tables (besides the usual issue description, input file and submission script), or (even better) the whole standard output.

There is a script located `here <https://github.com/AMReX-Codes/amrex/tree/development/Tools/TinyProfileParser>`__ that parses the Tiny Profiler output and generates a JSON file that can be used with `Hatchet <https://hatchet.readthedocs.io/en/latest/>`__ in order to analyze performance.


.. _developers-profiling-nsight-systems:

Nvidia Nsight-Systems
---------------------

`Vendor homepage <https://developer.nvidia.com/nsight-systems>`__ and `product manual <https://docs.nvidia.com/nsight-systems/>`__.

Example on how to create traces on a multi-GPU system that uses the Slurm scheduler (e.g., NERSC's Perlmutter system):

.. code-block:: bash

   rm -rf profiling*
   export OMP_NUM_THREADS=1

   # record
   srun --ntasks=4 --gpus=4 --cpu-bind=cores \
       nsys profile -f true \
         -o profiling_%q{SLURM_TASK_PID} \
         -t mpi,cuda,nvtx,osrt,openmp   \
         --mpi-impl=openmpi \
       ./warpx.3d.MPI.CUDA.DP.QED inputs_3d \
         warpx.numprocs=1 1 4 amr.n_cell=512 512 2048 max_step=10

 Example on how to create traces on a multi-GPU system that uses the ``jsrun`` scheduler (e.g., `OLCF's Summit system <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#optimizing-and-profiling>`__):

 .. code-block:: bash

    rm -rf profiling*
    export OMP_NUM_THREADS=1

    # record
    jsrun -n 4 -a 1 -g 1 -c 7 --bind=packed:$OMP_NUM_THREADS \
        nsys profile -f true \
          -o profiling_%p \
          -t mpi,cuda,nvtx,osrt,openmp   \
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
