.. _building-perlmutter:

Perlmutter (NERSC)
==================

.. warning::

   Perlmutter is still in acceptance testing.
   This page documents our internal testing workflow only.

The `Perlmutter cluster <https://www.nersc.gov/systems/perlmutter/>`_ is located at NERSC.

If you are new to this system, please see the following resources:

* `NERSC user guide <https://docs.nersc.gov/>`__
* Batch system: `Slurm <https://docs.nersc.gov/jobs/>`__
* `Jupyter service <https://docs.nersc.gov/services/jupyter/>`__
* `Production directories <https://www.nersc.gov/users/storage-and-file-systems/>`__:

  * ``$SCRATCH``: per-user production directory (20TB)
  * ``/global/cscratch1/sd/m3239``: shared production directory for users in the project ``m3239`` (50TB)
  * ``/global/cfs/cdirs/m3239/``: community file system for users in the project ``m3239`` (100TB)


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/perlmutter_warpx.profile``).

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load cmake
   module swap PrgEnv-nvidia PrgEnv-gnu
   module swap gcc gcc/9.3.0

   # newer CMake (3.22+) for HPE/Cray FindMPI.cmake improvement
   export PATH=$HOME/sw/cmake-3.22.0-dev/bin:$PATH

   # optional: just an additional text editor
   # module load nano

   # optional: FFT, I/O, Python, ...
   # TODO

   # optional: an alias to request an interactive node for two hours
   # TODO

   # optimize CUDA compilation for A100
   export AMREX_CUDA_ARCH=8.0

   # compiler environment hints
   export CC=$(which gcc)
   export CXX=$(which g++)
   export FC=$(which gfortran)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)


We recommend to store the above lines in a file, such as ``$HOME/perlmutter_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/perlmutter_warpx.profile

Furthermore, until Perlmutter provides a CMake 3.22.0+ module, we need to execute once:

.. code-block:: bash

   git clone https://gitlab.kitware.com/cmake/cmake.git $HOME/src/cmake
   cmake -S $HOME/src/cmake -B $HOME/src/cmake/build -DCMAKE_INSTALL_PREFIX=$HOME/sw/cmake-3.22.0-dev
   cmake --build $HOME/src/cmake/build --target install -j 16

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_OPENPMD=ON -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
   cmake --build build -j 10

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-perlmutter:

Running
-------

.. _running-cpp-perlmutter-A100-GPUs:

A100 GPUs
^^^^^^^^^

The batch script below can be used to run a WarpX simulation on multiple nodes on the supercomputer Perlmutter at NERSC.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that the only option so far is to run with one MPI rank per GPU.

.. literalinclude:: ../../../../Tools/BatchScripts/batch_perlmutter.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_perlmutter.sh`` and run

.. code-block:: bash

   bsub batch_perlmutter.sh

to submit the job.
