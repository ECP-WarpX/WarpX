.. _building-lumi:

LUMI (CSC)
==========

The `LUMI cluster <https://www.lumi-supercomputer.eu>`_ is located at CSC (Finland).
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

Introduction
------------

If you are new to this system, **please see the following resources**:

* `Lumi user guide <https://docs.lumi-supercomputer.eu>`_
* Batch system: `Slurm <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/>`_
* `Data analytics and visualization <https://docs.lumi-supercomputer.eu/hardware/lumid/>`__
* `Production directories <https://docs.lumi-supercomputer.eu/storage/>`__:

  * ``$HOME``: single user, intended to store user configuration files and personal data (20GB default quota)
  * ``/project/$proj``: shared with all members of a project, purged at the end of a project (50 GB default quota)
  * ``/scratch/$proj``: temporary storage, main storage to be used for disk I/O needs when running simulations on LUMI, purged every 90 days (50TB default quota)


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/lumi_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/lumi-csc/lumi_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lumi-csc/lumi_warpx.profile.example``.


We recommend to store the above lines in a file, such as ``$HOME/lumi_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/lumi_warpx.profile

And since LUMI does not yet provide a module for them, install c-blosc and ADIOS2:

.. code-block:: bash

   export CMAKE_PREFIX_PATH=${HOME}/sw/lumi/gpu/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
   export CMAKE_PREFIX_PATH=${HOME}/sw/lumi/gpu/adios2-2.8.3:$CMAKE_PREFIX_PATH

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-build
   cmake -S src/c-blosc -B src/c-blosc-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/lumi/gpu/c-blosc-1.21.1
   cmake --build src/c-blosc-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-build
   cmake -S src/adios2 -B src/adios2-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/lumi/gpu//adios2-2.8.3
   cmake --build src/adios2-build --target install -j 16

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS="1;2;3" -DWarpX_COMPUTE=HIP -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON
   cmake --build build -j 16

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
WarpX executables are now in ``build/bin/`` and :ref:`can be run <running-cpp-lumi>` with matching :ref:`example inputs files <usage-examples>`.
Most people execute the binary directly or copy it out to a location in ``/scratch/<project>``.


.. _running-cpp-lumi:

Running
-------

.. _running-cpp-lumi-MI250X-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

In non-interactive runs:

.. literalinclude:: ../../../../Tools/machines/lumi-csc/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lumi-csc/submit.sh``.

.. _post-processing-lumi:

Post-Processing
---------------

.. note::

   TODO: Document any Jupyter or data services.

Known System Issues
-------------------

.. warning::

   December 12th, 2022:
   There is a caching bug in libFabric that causes WarpX simulations to occasionally hang on LUMI on more than 1 node.

   As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

   .. code-block:: bash

      #export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
      # or, less invasive:
      export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

.. warning::

   January, 2023:
   We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.

   June, 2023:
   Although a fix was planned for ROCm 5.5, we still see the same issue in this release and continue to exchange with AMD and HPE on the issue.

   Stay with the ROCm 5.2 module to avoid a 2x slowdown.

.. warning::

   May 2023:
   rocFFT in ROCm 5.1-5.3 tries to `write to a cache <https://rocfft.readthedocs.io/en/latest/#runtime-compilation>`__ in the home area by default.
   This does not scale, disable it via:

   .. code-block:: bash

      export ROCFFT_RTC_CACHE_PATH=/dev/null
