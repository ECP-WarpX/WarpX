.. _building-frontier:

Frontier (OLCF)
===============

The `Frontier cluster (see: Crusher) <https://docs.olcf.ornl.gov/systems/crusher_quick_start_guide.html>`_ is located at OLCF.
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Crusher user guide <https://docs.olcf.ornl.gov/systems/frontier_user_guide.html>`_
* Batch system: `Slurm <https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#running-jobs>`_
* `Production directories <https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#data-and-storage>`_:

  * ``$PROJWORK/$proj/``: shared with all members of a project, purged every 90 days (recommended)
  * ``$MEMBERWORK/$proj/``: single user, purged every 90 days (usually smaller quota, 50TB default quota)
  * ``$WORLDWORK/$proj/``: shared with all users, purged every 90 days (50TB default quota)
  * Note that the ``$HOME`` directory is mounted as read-only on compute nodes.
    That means you cannot run in your ``$HOME``.
    It's default quota is 50GB.

Note: the Orion lustre filesystem on Frontier and the older Alpine GPFS filesystem on Summit are not mounted on each others machines.
Use `Globus <https://www.globus.org>`__ to transfer data between them if needed.


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/frontier_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/frontier-olcf/frontier_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/frontier-olcf/frontier_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/frontier_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/frontier_warpx.profile


Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_COMPUTE=HIP
   cmake --build build -j 32

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
A 3D WarpX executable is now in ``build/bin/`` and :ref:`can be run <running-cpp-frontier>` with a :ref:`3D example inputs file <usage-examples>`.
Most people execute the binary directly or copy it out to a location in ``$PROJWORK/$proj/``.


.. _running-cpp-frontier:

Running
-------

.. _running-cpp-frontier-MI250X-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

After requesting an interactive node with the ``getNode`` alias above, run a simulation like this, here using 8 MPI ranks and a single node:

.. code-block:: bash

   runNode ./warpx inputs

Or in non-interactive runs:

.. literalinclude:: ../../../../Tools/machines/frontier-olcf/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/frontier-olcf/submit.sh``.


.. _post-processing-frontier:

Post-Processing
---------------

For post-processing, most users use Python via OLCFs's `Jupyter service <https://jupyter.olcf.ornl.gov>`__ (`Docs <https://docs.olcf.ornl.gov/services_and_applications/jupyter/index.html>`__).

Please follow the same guidance as for :ref:`OLCF Summit post-processing <post-processing-summit>`.

.. _known-frontier-issues:

Known System Issues
-------------------

.. warning::

   May 16th, 2022 (OLCFHELP-6888):
   There is a caching bug in Libfabric that causes WarpX simulations to occasionally hang on Frontier on more than 1 node.

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

.. warning::

   January, 2023 (OLCFDEV-1284, AMD Ticket: ORNLA-130):
   We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.
   Reported to AMD and fixed for the 5.5 release of ROCm.

   Upgrade ROCm or stay with the ROCm 5.2 module to avoid.
