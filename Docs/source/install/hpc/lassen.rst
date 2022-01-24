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

We use the following modules and environments on the system (``$HOME/lassen_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/lassen-llnl/lassen_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lassen-llnl/lassen_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/lassen_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/lassen_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_COMPUTE=CUDA
   cmake --build build -j 10

This will build an executable in ``build/bin/``.
The other :ref:`general compile-time options <building-cmake>` apply as usual.


.. _running-cpp-lassen:

Running
-------

.. _running-cpp-lassen-V100-GPUs:

V100 GPUs (16GB)
^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Lassen at LLNL.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that the only option so far is to run with one MPI rank per GPU.

.. literalinclude:: ../../../../Tools/machines/lassen-llnl/lassen.bsub
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lassen-llnl/lassen.bsub``.

To run a simulation, copy the lines above to a file ``lassen.bsub`` and run

.. code-block:: bash

   bsub lassen.bsub

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on V100 GPUs for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=256`` and ``amr.blocking_factor=128``.

* **One MPI rank per GPU** (e.g., 4 MPI ranks for the 4 GPUs on each Lassen
  node)

* **Two `128x128x128` grids per GPU**, or **one `128x128x256` grid per GPU**.
