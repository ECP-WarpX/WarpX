.. _building-juwels:

Juwels (JSC)
============

The `Juwels supercomputer <https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUWELS/JUWELS_node.html>`_ is located at JSC.

See `this page <https://apps.fz-juelich.de/jsc/hps/juwels/quickintro.html>`__ for a quick introduction.

Installation
------------

We use the following modules and environments on the system.

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load GCC
   module load OpenMPI
   module load CUDA

Note that for now WarpX must rely on OpenMPI instead of the recommended MPI implementation on this platform MVAPICH2.

We recommend to store the above lines in a file, such as ``$HOME/warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/warpx.profile

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   mkdir ~/src
   cd ~/src

   git clone https://github.com/ECP-WarpX/WarpX.git warpx
   git clone --branch QED https://github.com/ECP-WarpX/picsar.git
   git clone --branch development https://github.com/AMReX-Codes/amrex.git

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   make -j 16 COMP=gcc USE_GPU=TRUE

The other :ref:`general compile-time options <building-source>` apply as usual.

Running
-------

An example submission script reads

.. literalinclude:: ../../../Tools/BatchScripts/batch_juwels.sh
   :language: bash

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.
