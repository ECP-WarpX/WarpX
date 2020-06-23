.. _building-summit:

Building WarpX on Summit (OLCF)
===============================

The `Summit cluster <https://www.olcf.ornl.gov/summit/>`_ is located at OLCF.

If you are new to this system, please see the following resources:

* `Summit user guide <https://docs.olcf.ornl.gov/systems/summit_user_guide.html>`_
* Batch system: `LSF <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#running-jobs>`_
* `Production directories <https://docs.olcf.ornl.gov/data/storage_overview.html>`_:

  * ``$PROJWORK/$proj/``: shared with all members of a project (recommended)
  * ``$MEMBERWORK/$proj/``: single user (usually smaller quota)
  * ``$WORLDWORK/$proj/``: shared with all users
  * Note that the ``$HOME`` directory is mounted as read-only on compute nodes.
    That means you cannot run in your ``$HOME``.


Installation
------------

We use the following modules and environments on the system.

.. code-block:: bash

   # please set your project account
   export proj=<yourProject>

   # required dependencies
   module load gcc/6.4.0
   module load cuda

   # optional: faster re-builds
   module load ccache

   # optional: for PSATD support
   module load fftw

   # optional: for QED support
   module load boost/1.66.0

   # optional: for openPMD support
   module load cmake
   module load hdf5/1.10.4
   module load adios2/2.5.0
   export PKG_CONFIG_PATH=$HOME/sw/openPMD-api-install/lib64/pkgconfig:$PKG_CONFIG_PATH
   export CMAKE_PREFIX_PATH=$HOME/sw/openPMD-api-install:$CMAKE_PREFIX_PATH

   # optional: Ascent in situ support
   #   note: you cannot yet use openPMD with parallel HDF5 and
   #         Ascent (with serial HDF5) at the same time
   export Alpine=/gpfs/alpine/world-shared/csc340/software/ascent/0.5.3-pre/summit/cuda/gnu
   export Ascent_DIR=$Alpine/ascent-install
   export Conduit_DIR=$Alpine/conduit-install
   #   work-around for ill-placed AscentConfig.cmake
   export Ascent_DIR=$Ascent_DIR/lib/cmake

   # optional: for Ascent support (work-in-progress)
   export ASCENT_DIR=/gpfs/alpine/world-shared/csc340/software/ascent/0.5.3-pre/summit/cuda/gnu/ascent-install/

   # optional: just an additional text editor
   module load nano

   # optional: an alias to request an interactive node for two hours
   alias getNode="bsub -P $proj -W 2:00 -nnodes 1 -Is /bin/bash"

   # fix system defaults: do not escape $ with a \ on tab completion
   shopt -s direxpand

   # compiler environment hints
   export CC=$(which gcc)
   export CXX=$(which g++)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)


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

Optionally, download and build openPMD-api for I/O:

.. code-block:: bash

   git clone https://github.com/openPMD/openPMD-api.git
   mkdir openPMD-api-build
   cd openPMD-api-build
   cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/openPMD-api-install/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_RPATH='$ORIGIN' -DMPIEXEC_EXECUTABLE=$(which jsrun)
   cmake --build . --target install --parallel 16

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   make -j 16 COMP=gcc USE_GPU=TRUE USE_OPENPMD=TRUE

The other :ref:`general compile-time options <building-source>` apply as usual.


Running
-------

Please see :ref:`our example job scripts <running-cpp-summit>` on how to run WarpX on Summit.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.
