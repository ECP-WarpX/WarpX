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

We use the following modules and environments on the system.

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load cmake/3.16.8
   module load intel/19.1.2
   module load mvapich2/2.3

   # optional: for PSATD support
   module load fftw/3.3.8

   # optional: for QED support
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

   rm -rf build/
   cmake -B build -DWarpX_OPENPMD=ON
   cmake --build build -j 6

This will build an executable in ``build/bin/``.
The other :ref:`general compile-time options <building-cmake>` apply as usual.


Running
-------

Please see :ref:`our example job scripts <running-cpp-quartz>` on how to run WarpX on Quartz.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.
