.. _developers-profiling:

Profiling the Code
==================

Profiling allows us to find the bottle-necks of the code as it is currently implemented.
Bottle-necks are the parts of the code that may delay the simulation, making it more computationally expensive.
Once found, we can update the related code sections and improve its efficiency.
Profiling tools can also be used to check how load balanced the simulation is, i.e. if the work is well distributed across all MPI ranks used.
Load balancing can be activated in WarpX by setting input parameters, see the :ref:`parallelization input parameter section <running-cpp-parameters-parallelization>`.

.. _developers-profiling-tiny-profiler:

AMReX's Tiny Profiler
---------------------

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

.. note::

   When creating performance-related issues on the WarpX GitHub repo, please include Tiny Profiler tables (besides the usual issue description, input file and submission script), or (even better) the whole standard output.

For more detailed information please visit the `AMReX profiling documentation <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools_Chapter.html>`__.
There is a script located `here <https://github.com/AMReX-Codes/amrex/tree/development/Tools/TinyProfileParser>`__ that parses the Tiny Profiler output and generates a JSON file that can be used with `Hatchet <https://hatchet.readthedocs.io/en/latest/>`__ in order to analyze performance.

AMReX's Full Profiler
---------------------

The Tiny Profiler provides a summary across all MPI ranks. However, when analyzing
load-balancing, it can be useful to have more detailed information about the
behavior of *each* individual MPI rank. The workflow for doing so is the following:

- Compile WarpX with full profiler support:

    .. code-block:: bash

       cmake -S . -B build -DAMReX_BASE_PROFILE=YES -DAMReX_TRACE_PROFILE=YES  -DAMReX_COMM_PROFILE=YES -DAMReX_TINY_PROFILE=OFF
       cmake --build build -j 4

    .. warning::

        Please note that the `AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`__ for ``AMReX_TINY_PROFILE`` (our default: ``ON``) and full profiling traces via ``AMReX_BASE_PROFILE`` are mutually exclusive.
        Further tracing options are sub-options of ``AMReX_BASE_PROFILE``.

        To turn on the tiny profiler again, remove the ``build`` directory or turn off ``AMReX_BASE_PROFILE`` again:

        .. code-block:: bash

            cmake -S . -B build -DAMReX_BASE_PROFILE=OFF -DAMReX_TINY_PROFILE=ON

- Run the simulation to be profiled. Note that the WarpX executable will create
and new folder `bl_prof`, which contains the profiling data.

    .. note::

        When using the full profiler, it is usually useful to profile only
        a few PIC iterations (e.g. 10-20 PIC iterations), in order to improve
        readability. If the interesting PIC iterations occur only late in a
        simulation, you can run the first part of the simulation without
        profiling, the create a checkpoint, and then restart the simulation for
        10-20 steps with the full profiler on.

.. note::

    The next steps can be done on a local computer (even if
    the simulation itself ran on an HPC cluster). In this
    case, simply copy the folder `bl_prof` to your local computer.

- In order, to visualize the profiling data, install `amrvis` using `spack`:

    .. code-block:: bash

        spack install amrvis dims=2 +profiling

- Then create timeline database from the `bl_prof` data and open it:

    .. code-block:: bash

        <amrvis-executable> -timelinepf bl_prof/
        <amrvis-executable> pltTimeline/

   In the above, `<amrvis-executable>` should be replaced by the actual of your
   `amrvis` executable, which can be found starting to type `amrvis` and then
   using Tab completion, in a Terminal.

- This will pop-up a window with the timeline. Here are few guidelines to navigate it:
    - Use the horizontal scroller to find the area where the 10-20 PIC steps occur.
    - In order to zoom on an area, you can drag and drop with the mouse, and the hit `Ctrl-S` on a keyboard.
    - You can directly click on the timeline to see which actual MPI call is being perform. (Note that the colorbar can be misleading.)

.. _developers-profiling-nsight-systems:

Nvidia Nsight-Systems
---------------------

`Vendor homepage <https://developer.nvidia.com/nsight-systems>`__ and `product manual <https://docs.nvidia.com/nsight-systems/>`__.


Perlmutter Example
""""""""""""""""""

Example on how to create traces on a multi-GPU system that uses the Slurm scheduler (e.g., NERSC's Perlmutter system):

.. code-block:: bash

   # GPU-aware MPI
   export MPICH_GPU_SUPPORT_ENABLED=1
   # 1 OpenMP thread
   export OMP_NUM_THREADS=1

   export TMPDIR="$PWD/tmp"
   rm -rf ${TMPDIR} profiling*
   mkdir -p ${TMPDIR}

   # 2021.5.1 (broken: lacks most trace info)
   #NSYS="/global/common/software/nersc/pm-2021q4/easybuild/software/Nsight-Systems/2021.5.1/bin/nsys"
   # 2021.4.1 (working)
   NSYS="/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/bin/nsys"

   # record
   srun --ntasks=4 --gpus=4 --cpu-bind=cores \
       ${NSYS} profile -f true               \
         -o profiling_%q{SLURM_TASK_PID}     \
         -t mpi,cuda,nvtx,osrt,openmp        \
         --mpi-impl=mpich                    \
       ./warpx.3d.MPI.CUDA.DP.QED            \
         inputs_3d                           \
           warpx.numprocs=1 1 4 amr.n_cell=512 512 2048 max_step=10

.. warning::

   March 23rd, 2022 (INC0182505):
   Currently, the environment pre-loads a ``Nsight-Systems/2021.5.1`` module that ships ``nsys`` version 2021.5.1.
   This version does not record all trace information.
   You need to use the one directly shipped with the NVHPC base system, version 2021.4.1, located in ``/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/bin/nsys``.


Summit Example
""""""""""""""

 Example on how to create traces on a multi-GPU system that uses the ``jsrun`` scheduler (e.g., `OLCF's Summit system <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#optimizing-and-profiling>`__):

.. code-block:: bash

    # nsys: remove old traces
    rm -rf profiling* tmp-traces
    # nsys: a location where we can write temporary nsys files to
    export TMPDIR=$PWD/tmp-traces
    mkdir -p $TMPDIR
    # WarpX: one OpenMP thread per MPI rank
    export OMP_NUM_THREADS=1

    # record
    jsrun -n 4 -a 1 -g 1 -c 7 --bind=packed:$OMP_NUM_THREADS \
        nsys profile -f true \
          -o profiling_%p \
          -t mpi,cuda,nvtx,osrt,openmp   \
          --mpi-impl=openmpi             \
        ./warpx.3d.MPI.CUDA.DP.QED inputs_3d \
          warpx.numprocs=1 1 4 amr.n_cell=512 512 2048 max_step=10

.. warning::

   Sep 10th, 2021 (OLCFHELP-3580):
   The Nsight-Compute (``nsys``) version installed on Summit does not record details of GPU kernels.
   This is reported to Nvidia and OLCF.

Details
"""""""

In these examples, the individual lines for recording a trace profile are:

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
