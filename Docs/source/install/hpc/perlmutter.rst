.. _building-perlmutter:

Perlmutter (NERSC)
==================

.. warning::

   Perlmutter is still in acceptance testing and environments change often.
   Please reach visit this page often for updates and reach out to us if something needs an update.

The `Perlmutter cluster <https://docs.nersc.gov/systems/perlmutter/>`_ is located at NERSC.

If you are new to this system, please see the following resources:

* `NERSC user guide <https://docs.nersc.gov/>`__
* Batch system: `Slurm <https://docs.nersc.gov/systems/perlmutter/#running-jobs>`__
* `Jupyter service <https://docs.nersc.gov/services/jupyter/>`__
* `Production directories <https://docs.nersc.gov/filesystems/perlmutter-scratch/>`__:

  * ``$PSCRATCH``: per-user production directory (<TBD>TB)
  * ``/global/cscratch1/sd/m3239``: shared production directory for users in the project ``m3239`` (50TB)
  * ``/global/cfs/cdirs/m3239/``: community file system for users in the project ``m3239`` (100TB)


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/perlmutter_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/perlmutter-nersc/perlmutter_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/perlmutter-nersc/perlmutter_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/perlmutter_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/perlmutter_warpx.profile

And since Perlmutter does not yet provide a module for it, install ADIOS2:

.. code-block:: bash

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-pm-build
   cmake -S src/c-blosc -B src/c-blosc-pm-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/c-blosc-1.21.1
   cmake --build src/c-blosc-pm-build --target install --parallel 32

   # ADIOS2
   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-pm-build
   cmake -S src/adios2 -B src/adios2-pm-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/adios2-2.7.1
   cmake --build src/adios2-pm-build --target install -j 32

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
   cmake --build build -j 32

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-perlmutter:

Running
-------

.. _running-cpp-perlmutter-A100-GPUs:

A100 GPUs (40 GB)
^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on multiple nodes (change ``-N`` accordingly) on the supercomputer Perlmutter at NERSC.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that we run one MPI rank per GPU.


.. literalinclude:: ../../../../Tools/machines/perlmutter-nersc/perlmutter.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/perlmutter-nersc/perlmutter.sbatch``.

To run a simulation, copy the lines above to a file ``perlmutter.sbatch`` and run

.. code-block:: bash

   sbatch perlmutter.sbatch

to submit the job.


.. _post-processing-perlmutter:

Post-Processing
---------------

For post-processing, most users use Python via NERSC's `Jupyter service <https://jupyter.nersc.gov>`__ (`Docs <https://docs.nersc.gov/services/jupyter/>`__).

Please follow the same process as for :ref:`NERSC Cori post-processing <post-processing-cori>`.
**Important:** The *environment + Jupyter kernel* must separate from the one you create for Cori.

The Perlmutter ``$PSCRATCH`` filesystem is currently not yet available on Jupyter.
Thus, store or copy your data to Cori's ``$SCRATCH`` or use the Community FileSystem (CFS) for now.
