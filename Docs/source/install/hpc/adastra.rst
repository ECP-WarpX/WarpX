.. _building-adastra:

Adastra (CINES)
===============

The `Adastra cluster <https://www.cines.fr/calcul/adastra/>`_ is located at CINES (France).
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Adastra user guide <https://dci.dci-gitlab.cines.fr/webextranet/user_support/>`_
* Batch system: `Slurm <https://dci.dci-gitlab.cines.fr/webextranet/user_support/index.html?highlight=sbatch#running-jobs>`_
* `Production directories <https://dci.dci-gitlab.cines.fr/webextranet/data_and_storage/index.html#data-and-storage>`_:

  * ``$SHAREDSCRATCHDIR``: meant for short-term data storage, shared with all members of a project, purged every 30 days (17.6 TB default quota)
  * ``$SCRATCHDIR``: meant for short-term data storage, single user, purged every 30 days
  * ``$SHAREDWORKDIR``: meant for mid-term data storage, shared with all members of a project, never purged (4.76 TB default quota)
  * ``$WORKDIR``: meant for mid-term data storage, single user, never purged
  * ``$STORE`` : meant for long term storage, single user, never purged, backed up
  * ``$SHAREDHOMEDIR`` : meant for scripts and tools, shared with all members of a project, never purged, backed up
  * ``$HOME`` : meant for scripts and tools, single user, never purged, backed up


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/adastra_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/adastra-cines/adastra_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/adastra-cines/adastra_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/adastra_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/adastra_warpx.profile

And since Adastra does not yet provide a module for them, install c-blosc and ADIOS2:

.. code-block:: bash

   export CMAKE_PREFIX_PATH=${HOME}/sw/adastra/gpu/c-blosc-1.21.1:$CMAKE_PREFIX_PATH
   export CMAKE_PREFIX_PATH=${HOME}/sw/adastra/gpu/adios2-2.8.3:$CMAKE_PREFIX_PATH

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-pm-build
   cmake -S src/c-blosc -B src/c-blosc-pm-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/adastra/gpu/c-blosc-1.21.1
   cmake --build src/c-blosc-pm-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-pm-build
   cmake -S src/adios2 -B src/adios2-pm-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/adastra/gpu/adios2-2.8.3
   cmake --build src/adios2-pm-build --target install -j 16

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS="1;2;3" -DWarpX_COMPUTE=HIP -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON
   cmake --build build -j 32

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
A 3D WarpX executable is now in ``build/bin/`` and :ref:`can be run <running-cpp-adastra>` with a :ref:`3D example inputs file <usage-examples>`.
Most people execute the binary directly or copy it out to a location in ``$WORKDIR`` or ``$SCRATCHDIR``.


.. _running-cpp-adastra:

Running
-------

.. _running-cpp-adastra-MI250X-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

In non-interactive runs:

.. literalinclude:: ../../../../Tools/machines/adastra-cines/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/adastra-cines/submit.sh``.


.. _post-processing-adastra:

Post-Processing
---------------

.. note::

   TODO: Document any Jupyter or data services.

.. _known-adastra-issues:

Known System Issues
-------------------

.. warning::

   May 16th, 2022:
   There is a caching bug in Libfabric that causes WarpX simulations to occasionally hang on on more than 1 node.

   As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

   .. code-block:: bash

      #export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
      # or, less invasive:
      export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

.. warning::

   Sep 2nd, 2022:
   rocFFT in ROCm 5.1-5.3 tries to `write to a cache <https://rocfft.readthedocs.io/en/latest/#runtime-compilation>`__ in the home area by default.
   This does not scale, disable it via:

   .. code-block:: bash

      export ROCFFT_RTC_CACHE_PATH=/dev/null

.. warning::

   January, 2023:
   We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.
   Reported to AMD and fixed for the next release of ROCm.

   Stay with the ROCm 5.2 module to avoid.
