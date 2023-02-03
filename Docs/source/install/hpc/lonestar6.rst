.. _building-lonestar6:

Lonestar6 (TACC)
================

The `Lonestar6 cluster <https://portal.tacc.utexas.edu/user-guides/lonestar6>`_ is located at `TACC <https://www.tacc.utexas.edu>`__.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `TACC user guide <https://portal.tacc.utexas.edu/user-guides/>`__
* Batch system: `Slurm <https://portal.tacc.utexas.edu/user-guides/lonestar6#job-management>`__
* `Jupyter service <https://tacc.github.io/ctls2017/docs/intro_to_python/intro_to_python_011_jupyter.html>`__
* `Filesystem directories <https://portal.tacc.utexas.edu/user-guides/lonestar6#managing-files-on-lonestar6>`__:

  * ``$HOME``: per-user home directory, backed up (10 GB)
  * ``$WORK``: per-user production directory, not backed up, not purged, Lustre (1 TB)
  * ``$SCRATCH``: per-user production directory, not backed up, purged every 10 days, Lustre (no limits, 8PByte total)


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $WORK/src/warpx

We use the following modules and environments on the system (``$HOME/lonestar6_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/lonestar6-tacc/lonestar6_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lonestar6-tacc/lonestar6_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/lonestar6_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/lonestar6_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $WORK/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA -DWarpX_PSATD=ON
   cmake --build build -j 16

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
A 3D WarpX executable is now in ``build/bin/`` and :ref:`can be run <running-cpp-lonestar6>` with a :ref:`3D example inputs file <usage-examples>`.
Most people execute the binary directly or copy it out to a location in ``$WORK`` or ``SCRATCH``.


.. _running-cpp-lonestar6:

Running
-------

.. _running-cpp-lonestar6-A100-GPUs:

A100 GPUs (40 GB)
^^^^^^^^^^^^^^^^^

`32 GPU nodes, each with 3 A100 GPUs (40 GB) <https://portal.tacc.utexas.edu/user-guides/lonestar6#system-gpu>`__.

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
