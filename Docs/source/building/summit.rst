.. _building-summit:

Summit (OLCF)
=============

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
   module load cmake
   module load gcc/6.4.0
   module load cuda

   # optional: faster re-builds
   module load ccache

   # optional: for PSATD support
   module load fftw

   # optional: for QED lookup table generation support
   module load boost/1.66.0

   # optional: for openPMD support
   module load ums
   module load ums-aph114
   module load openpmd-api/0.12.0

   # optional: for PSATD in RZ geometry support
   #   note: needs the ums modules above
   module load blaspp
   module load lapackpp

   # optional: Ascent in situ support
   #   note: build WarpX with CMake
   export Alpine=/gpfs/alpine/csc340/world-shared/software/ascent/2020-11-07/summit/cuda/gnu
   export Ascent_DIR=$Alpine/ascent-install
   export Conduit_DIR=$Alpine/conduit-install

   # optional: for Python bindings or libEnsemble
   module load python/3.7.0

   # optional: for libEnsemble
   module load openblas/0.3.9-omp
   module load netlib-lapack/3.8.0
   if [ -d "$HOME/sw/venvs/warpx-libE" ]
   then
     source $HOME/sw/venvs/warpx-libE/bin/activate
   fi

   # optional: just an additional text editor
   module load nano

   # optional: an alias to request an interactive node for two hours
   alias getNode="bsub -P $proj -W 2:00 -nnodes 1 -Is /bin/bash"

   # fix system defaults: do not escape $ with a \ on tab completion
   shopt -s direxpand

   # optimize CUDA compilation for V100
   export AMREX_CUDA_ARCH=7.0

   # compiler environment hints
   export CC=$(which gcc)
   export CXX=$(which g++)
   export FC=$(which gfortran)
   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)


We recommend to store the above lines in a file, such as ``$HOME/warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/warpx.profile

Optionally, download and build openPMD-api for I/O (only needed if you did not load our module above):

.. code-block:: bash

   git clone https://github.com/openPMD/openPMD-api.git
   mkdir openPMD-api-build
   cd openPMD-api-build
   cmake ../openPMD-api -DopenPMD_USE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/openPMD-api-install/ -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_RPATH='$ORIGIN' -DMPIEXEC_EXECUTABLE=$(which jsrun)
   cmake --build . --target install --parallel 16

Optionally, download and install :ref:`libEnsemble <libensemble>` for dynamic ensemble optimizations:

.. code-block:: bash

   export BLAS=$OLCF_OPENBLAS_ROOT/lib/libopenblas.so
   export LAPACK=$OLCF_NETLIB_LAPACK_ROOT/lib64/liblapack.so
   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv
   python3 -m venv $HOME/sw/venvs/warpx-libE
   source $HOME/sw/venvs/warpx-libE/bin/activate
   python3 -m pip install --upgrade pip
   python3 -m pip install --upgrade cython
   python3 -m pip install --upgrade numpy
   python3 -m pip install --upgrade scipy
   python3 -m pip install --upgrade mpi4py --no-binary mpi4py
   python3 -m pip install --upgrade -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   mkdir -p build
   cd build
   cmake .. -DWarpX_OPENPMD=ON -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
   make -j 10

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


Running
-------

Please see :ref:`our example job scripts <running-cpp-summit>` on how to run WarpX on Summit.

See :doc:`../visualization/yt` for more information on how to visualize the simulation results.
