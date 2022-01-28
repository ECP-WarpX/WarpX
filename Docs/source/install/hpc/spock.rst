.. _building-spock:

Spock (OLCF)
=============

The `Spock cluster <https://docs.olcf.ornl.gov/systems/spock_quick_start_guide.html>`_ is located at OLCF.

If you are new to this system, please see the following resources:

* `Spock user guide <https://docs.olcf.ornl.gov/systems/spock_quick_start_guide.html>`_
* Batch system: `Slurm <https://docs.olcf.ornl.gov/systems/spock_quick_start_guide.html#running-jobs>`_
* `Production directories <https://docs.olcf.ornl.gov/data/storage_overview.html>`_:

  * ``$PROJWORK/$proj/``: shared with all members of a project (recommended)
  * ``$MEMBERWORK/$proj/``: single user (usually smaller quota)
  * ``$WORLDWORK/$proj/``: shared with all users
  * Note that the ``$HOME`` directory is mounted as read-only on compute nodes.
    That means you cannot run in your ``$HOME``.


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/spock_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/spock-olcf/spock_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/spock-olcf/spock_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/spock_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/spock_warpx.profile


Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=HIP -DAMReX_AMD_ARCH=gfx908 -DMPI_CXX_COMPILER=$(which CC) -DMPI_C_COMPILER=$(which cc) -DMPI_COMPILER_FLAGS="--cray-print-opts=all"
   cmake --build build -j 10

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-spock:

Running
-------

.. _running-cpp-spock-MI100-GPUs:

MI100 GPUs (32 GB)
^^^^^^^^^^^^^^^^^^

After requesting an interactive node with the ``getNode`` alias above, run a simulation like this, here using 4 MPI ranks:

.. code-block:: bash

   srun -n 4 -c 2 --ntasks-per-node=4 ./warpx inputs

Or in non-interactive runs started with ``sbatch``:

.. literalinclude:: ../../../../Tools/machines/spock-olcf/spock_mi100.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/spock-olcf/spock_mi100.sbatch``.

We can currently use up to ``4`` nodes with ``4`` GPUs each (maximum: ``-N 4 -n 16``).


.. _post-processing-spock:

Post-Processing
---------------

For post-processing, most users use Python via OLCFs's `Jupyter service <https://jupyter.olcf.ornl.gov>`__ (`Docs <https://docs.olcf.ornl.gov/services_and_applications/jupyter/index.html>`__).

Please follow the same guidance as for :ref:`OLCF Summit post-processing <post-processing-summit>`.
