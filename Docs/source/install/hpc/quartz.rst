.. _building-quartz:

Quartz (LLNL)
=============

The `Quartz Intel CPU cluster <https://hpc.llnl.gov/hardware/platforms/quartz>`_ is located at LLNL.


Introduction
------------

If you are new to this system, **please see the following resources**:

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

We use the following modules and environments on the system (``$HOME/quartz_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/quartz-llnl/quartz_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/quartz-llnl/quartz_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/quartz_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/quartz_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build
   cmake --build build -j 6

This will build an executable in ``build/bin/``.
The other :ref:`general compile-time options <building-cmake>` apply as usual.


Running
-------

.. _running-cpp-quartz-CPUs:

Intel Xeon E5-2695 v4 CPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Quartz at LLNL.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.

.. literalinclude:: ../../../../Tools/machines/quartz-llnl/quartz.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/quartz-llnl/quartz.sbatch``.

To run a simulation, copy the lines above to a file ``quartz.sbatch`` and run

.. code-block:: bash

   sbatch quartz.sbatch

to submit the job.
