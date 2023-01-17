.. _building-crusher:

Crusher (OLCF)
==============

The `Crusher cluster <https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html>`_ is located at OLCF.
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Crusher user guide <https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html>`_
* Batch system: `Slurm <https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html#running-jobs>`_
* `Production directories <https://docs.olcf.ornl.gov/data/index.html#data-storage-and-transfers>`_:

  * ``$PROJWORK/$proj/``: shared with all members of a project, purged every 90 days (recommended)
  * ``$MEMBERWORK/$proj/``: single user, purged every 90 days (usually smaller quota)
  * ``$WORLDWORK/$proj/``: shared with all users, purged every 90 days
  * Note that the ``$HOME`` directory is mounted as read-only on compute nodes.
    That means you cannot run in your ``$HOME``.


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/crusher_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/crusher-olcf/crusher_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/crusher-olcf/crusher_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/crusher_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/crusher_warpx.profile

And since Crusher does not yet provide a module for them, install BLAS++ and LAPACK++:

.. code-block:: bash

   # BLAS++ (for PSATD+RZ)
   git clone https://github.com/icl-utk-edu/blaspp.git src/blaspp
   rm -rf src/blaspp-crusher-build
   cmake -S src/blaspp -B src/blaspp-crusher-build -Duse_openmp=OFF -Dgpu_backend=hip -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/crusher/blaspp-master
   cmake --build src/blaspp-crusher-build --target install --parallel 10

   # LAPACK++ (for PSATD+RZ)
   git clone https://github.com/icl-utk-edu/lapackpp.git src/lapackpp
   rm -rf src/lapackpp-crusher-build
   cmake -S src/lapackpp -B src/lapackpp-crusher-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$HOME/sw/crusher/lapackpp-master
   cmake --build src/lapackpp-crusher-build --target install --parallel 10

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=HIP
   cmake --build build -j 10

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-crusher:

Running
-------

.. _running-cpp-crusher-MI100-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

ECP WarpX project members, use the ``aph114`` project ID.

After requesting an interactive node with the ``getNode`` alias above, run a simulation like this, here using 8 MPI ranks and a single node:

.. code-block:: bash

   runNode ./warpx inputs

Or in non-interactive runs:

.. literalinclude:: ../../../../Tools/machines/crusher-olcf/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/crusher-olcf/submit.sh``.


.. _post-processing-crusher:

Post-Processing
---------------

For post-processing, most users use Python via OLCFs's `Jupyter service <https://jupyter.olcf.ornl.gov>`__ (`Docs <https://docs.olcf.ornl.gov/services_and_applications/jupyter/index.html>`__).

Please follow the same guidance as for :ref:`OLCF Summit post-processing <post-processing-summit>`.

.. _known-crusher-issues:

Known System Issues
-------------------

.. warning::

   May 16th, 2022 (OLCFHELP-6888):
   There is a caching bug in Libfrabric that causes WarpX simulations to occasionally hang on Crusher on more than 1 node.

   As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

   .. code-block:: bash

      #export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
      # or, less invasive:
      export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

.. warning::

   Sep 2nd, 2022 (OLCFDEV-1079):
   rocFFT in ROCm 5.1+ tries to `write to a cache <https://rocfft.readthedocs.io/en/latest/library.html#runtime-compilation>`__ in the home area by default.
   This does not scale, disable it via:

   .. code-block:: bash

      export ROCFFT_RTC_CACHE_PATH=/dev/null
