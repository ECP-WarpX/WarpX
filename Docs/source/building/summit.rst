.. _building-summit:

Building WarpX on Summit (OLCF)
================================

For the `Summit cluster
<https://www.olcf.ornl.gov/summit/>`__ at OLCF,
use the following commands to download the source code, and switch to the
correct branch:

::

    mkdir warpx_directory
    cd warpx_directory

    git clone https://github.com/ECP-WarpX/WarpX.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone --branch development https://github.com/AMReX-Codes/amrex.git

Then, ``cd`` into the directory ``WarpX`` and use the following set of commands to compile:

::

    module load gcc
    module load cuda
    make -j 16 COMP=gcc USE_GPU=TRUE

See :doc:`../running_cpp/platforms` for more information on how to run WarpX on Summit.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.


.. _building-cori-openPMD:

Building WarpX with openPMD support
-----------------------------------

First, load the appropriate modules:

.. code-block:: bash

    module load gcc
    module load cuda
    module load cmake
    module load hdf5/1.10.4

Then, in the ``warpx_directory``, download and build openPMD-api:

.. code-block:: bash

   git clone https://github.com/openPMD/openPMD-api.git
   mkdir openPMD-api-build
   cd openPMD-api-build
   cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=../openPMD-install/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_RPATH='$ORIGIN'
   cmake --build . --target install --parallel 16

.. note:

   On Summit, only compute nodes provide the infiniband hardware that Summit's MPI module expects, ``jsrun`` must be used on Summit instead of ``mpiexec``, and ``$HOME`` directories are read-only when computing.
   In order to run openPMD-api unit tests, run on a compute node inside ``$PROJWORK``, e.g. via ``bsub -P <addYourProjectID> -W 2:00 -nnodes 1 -Is /bin/bash``, and add ``-DMPIEXEC_EXECUTABLE=$(which jsrun)`` to the CMake options.

Finally, compile WarpX:

.. code-block:: bash

   cd ../WarpX
   export PKG_CONFIG_PATH=$PWD/../openPMD-install/lib64/pkgconfig:$PKG_CONFIG_PATH
   export CMAKE_PREFIX_PATH=$PWD/../openPMD-install:$CMAKE_PREFIX_PATH
   make -j 16 COMP=gcc USE_GPU=TRUE USE_OPENPMD=TRUE
