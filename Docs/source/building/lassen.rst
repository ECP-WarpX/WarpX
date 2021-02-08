.. _building-lassen:

Lassen (LLNL)
=============

The `Lassen V100 GPU cluster <https://hpc.llnl.gov/hardware/platforms/lassen>`_ is located at LLNL.

If you are new to this system, please see the following resources:

* `LLNL user account <https://lc.llnl.gov/lorenz/mylc/mylc.cgi>`_
* `Lassen user guide <https://hpc.llnl.gov/training/tutorials/using-lcs-sierra-system>`_
* Batch system: `LSF <https://hpc.llnl.gov/training/tutorials/using-lcs-sierra-system#batch-system>`_
* `Production directories <https://hpc.llnl.gov/hardware/file-systems>`_:

  * ``/p/gpfs1/$(whoami)``: personal directory on the parallel filesystem
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
   module load cmake/3.18.0
   module load gcc/8.3.1
   module load cuda/11.0.2

   # optional: for PSATD support
   module load fftw/3.3.8

   # optional: for QED lookup table generation support
   module load boost/1.70.0

   # optional: for openPMD support
   # TODO ADIOS2
   module load hdf5-parallel/1.10.4

   # optional: for PSATD in RZ geometry support
   # TODO: blaspp lapackpp

   # optional: for Python bindings
   module load python/3.8.2

   # optional: an alias to request an interactive node for two hours
   alias getNode="bsub -G $proj -W 2:00 -nnodes 1 -Is /bin/bash"

   # fix system defaults: do not escape $ with a \ on tab completion
   shopt -s direxpand

   # optimize CUDA compilation for V100
   export AMREX_CUDA_ARCH=7.0

   # compiler environment hints
   export CC=$(which gcc)
   export CXX=$(which g++)
   export FC=$(which gfortran)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)


We recommend to store the above lines in a file, such as ``$HOME/lassen_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/lassen_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   rm -rf build/
   cmake -B build -DWarpX_COMPUTE=CUDA -DWarpX_OPENPMD=ON
   cmake --build build -j 10

This will build an executable in ``build/bin/``.
The other :ref:`general compile-time options <building-cmake>` apply as usual.


Running
-------

Please see :ref:`our example job scripts <running-cpp-lassen>` on how to run WarpX on Lassen.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.
