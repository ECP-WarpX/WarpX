.. _building-juwels:

Juwels (JSC)
============

.. note::

   For the moment, WarpX doesn't run on Juwels with MPI_THREAD_MULTIPLE.
   Please compile with this compilation flag: ``MPI_THREAD_MULTIPLE=FALSE``.

The `Juwels supercomputer <https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/JUWELS_node.html>`_ is located at JSC.

See `this page <https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html>`_ for a quick introduction.

* Batch system: `Slurm <https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html#batch-system-on-system-name>`_
* `Production directories <https://apps.fz-juelich.de/jsc/hps/juwels/environment.html?highlight=scratch#available-filesystems>`_:

  * ``$SCRATCH``: Scratch filesystem for temporary data (90 day purge)
  * ``$FASTDATA/``: Storage location for large data (backuped)
  * Note that the ``$HOME`` directory is not designed for simulation runs and producing output there will impact performance.

Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   mkdir ~/src
   cd ~/src

   git clone https://github.com/ECP-WarpX/WarpX.git warpx
   git clone https://github.com/ECP-WarpX/picsar.git
   git clone https://github.com/AMReX-Codes/amrex.git

We use the following modules and environments on the system.

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load GCC
   module load OpenMPI
   module load CUDA
   module load Python

   # JUWELS' job scheduler may not map ranks to GPUs,
   # so we give a hint to AMReX about the node layout.
   # This is usually done in Make.<supercomputing center> files in AMReX
   # but there is no such file for JSC yet.
   export GPUS_PER_SOCKET=2
   export GPUS_PER_NODE=4

   # optimize CUDA compilation for V100 (CMake hint)
   export AMREX_CUDA_ARCH=7.0

Note that for now WarpX must rely on OpenMPI instead of the recommended MPI implementation on this platform MVAPICH2.

We recommend to store the above lines in a file, such as ``$HOME/warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   make -j 16 COMP=gcc USE_GPU=TRUE MPI_THREAD_MULTIPLE=FALSE

The other :ref:`general compile-time options <building-source>` apply as usual.

Running
-------

An example submission script reads

.. literalinclude:: ../../../Tools/BatchScripts/batch_juwels.sh
   :language: bash

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.

Using openPMD
-------------

In order to compile for openPMD output and HDF5 format, please use (in addition to modules above)

.. code-block:: bash

   module load HDF5

and use CMake for the compilation: in the WarpX repo

.. code-block:: bash

   rm -rf build # clean previous build
   mkdir build
   cd build
   cmake .. -DWarpX_COMPUTE=CUDA -DWarpX_MPI_THREAD_MULTIPLE=OFF -DWarpX_OPENPMD=ON
   make -j 16

This will take a while. The executable will be generated in ``<WarpX repo>/build/bin/``.

.. note::

   Currently, you need to add ``OMPI_MCA_io=romio321 `` at the beginning of your ``srun`` command to use HDF5 output.
