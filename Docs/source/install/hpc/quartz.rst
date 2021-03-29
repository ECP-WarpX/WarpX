.. _building-quartz:

Quartz (LLNL)
=============

The `Quartz Intel CPU cluster <https://hpc.llnl.gov/hardware/platforms/quartz>`_ is located at LLNL.

If you are new to this system, please see the following resources:

* `LLNL user account <https://lc.llnl.gov/lorenz/mylc/mylc.cgi>`_
* `Quartz user guide <https://computing.llnl.gov/tutorials/linux_clusters/>`_
* Batch system: `Slurm <https://computing.llnl.gov/tutorials/moab/>`_
* `Production directories <https://hpc.llnl.gov/hardware/file-systems>`_:

  * ``/p/lustre1/$(whoami)`` and ``/p/lustre2/$(whoami)``: personal directory on the parallel filesystem
  * Note that the ``$HOME`` directory and the ``/usr/workspace/$(whoami)`` space are NFS mounted and *not* suitable for production quality data generation.


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/warpx.profile``).

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load cmake/3.18.0
   module load intel/19.1.2
   module load mvapich2/2.3

   # optional: for PSATD support
   module load fftw/3.3.8

   # optional: for QED lookup table generation support
   module load boost/1.73.0

   # optional: for openPMD support
   # TODO ADIOS2
   module load hdf5-parallel/1.10.2

   # optional: for PSATD in RZ geometry support
   # TODO: blaspp lapackpp

   # optional: for Python bindings
   module load python/3.8.2

   # optional: an alias to request an interactive node for two hours
   alias getNode="srun --time=0:30:00 --nodes=1 --ntasks-per-node=2 --cpus-per-task=18 -p pdebug --pty bash"

   # fix system defaults: do not escape $ with a \ on tab completion
   shopt -s direxpand

   # compiler environment hints
   export CC=$(which icc)
   export CXX=$(which icpc)
   export FC=$(which ifort)
   # we need a newer libstdc++:
   export CFLAGS="-gcc-name=/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc ${CFLAGS}"
   export CXXFLAGS="-gxx-name=/usr/tce/packages/gcc/gcc-8.3.1/bin/g++ ${CXXFLAGS}"


We recommend to store the above lines in a file, such as ``$HOME/quartz_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/quartz_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_OPENPMD=ON
   cmake --build build -j 6

This will build an executable in ``build/bin/``.
The other :ref:`general compile-time options <building-cmake>` apply as usual.


Running
-------

.. _running-cpp-quartz-CPUs:

Intel Xeon E5-2695 v4 CPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Quartz at LLNL.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.

.. literalinclude:: ../../../../Tools/BatchScripts/batch_quartz.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_quartz.sh`` and run

.. code-block:: bash

   sbatch batch_quartz.sh

to submit the job.
