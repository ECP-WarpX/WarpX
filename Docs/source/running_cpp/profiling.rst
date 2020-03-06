Profiling the code
==================

Profiling allows us to find the bottle-necks of the code as it is currently implemented.
Bottle-necks are the parts of the code that may delay the simulation, making it more computationally expensive.
Once found, we can update the related code sections and improve its efficiency.
Profiling tools can also be used to check how load balanced the simulation is, i.e. if the work is well distributed accross all MPI ranks used.
Load balancing can be activated in WarpX by setting the input parameter `warpx.load_balance_int <https://warpx.readthedocs.io/en/latest/running_cpp/parameters.html#distribution-across-mpi-ranks-and-parallelization>`__.

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
* maximum percentage of time spent, accross all ranks

If the simulation is well load balanced the minimum, average and maximum times should be identical.

The top two tables refer to the complete simulation information.
The bottom two are related to the Evolve() section of the code (where each time step is computed).

Each set of two timers show the exclusive, top, and inclusive, bottom, information depending on wether the time spent in nested sections of the codes are included.

For more detailed information please visit the `AMReX profiling documentation <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools_Chapter.html>`__.

.. node:
   We encourage all users to include the TINYPROFILERS simulation information in their github issues, if they verified strong load imbalance.


Profiling the code with Intel Advisor on NERSC
----------------------------------------------

Follow these steps:

- Instrument the code during compilation

  ::

     module swap craype-haswell craype-mic-knl
     make -j 16 COMP=intel USE_VTUNE=TRUE

  (where the first line is only needed for KNL)

- In your SLURM submission script, use the following
  lines in order to run the executable. (In addition
  to setting the usual ``OMP`` environment variables.)

  ::

     module load advisor
     export ADVIXE_EXPERIMENTAL=roofline
     srun -n <n_mpi> -c <n_logical_cores_per_mpi> --cpu_bind=cores advixe-cl -collect survey -project-dir advisor -trace-mpi -- <warpx_executable> inputs
     srun -n <n_mpi> -c <n_logical_cores_per_mpi> --cpu_bind=cores advixe-cl -collect tripcounts -flop -project-dir advisor -trace-mpi -- <warpx_executable> inputs

  where ``<n_mpi>`` and ``<n_logical_cores_per_mpi>`` should be replaced by
  the proper values, and ``<warpx_executable>`` should be replaced by the
  name of the WarpX executable.

- Launch the Intel Advisor GUI

  ::

     module load advisor
     advixe-gui

  (Note: this requires to use ``ssh -XY`` when connecting to Cori.)

