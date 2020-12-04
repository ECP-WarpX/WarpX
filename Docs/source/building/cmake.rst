.. _building-cmake:

CMake Transition
================

We are currently transitioning to support `CMake <https://cmake.org>`_ as our primary build system.
Until we have transitioned our documentation and functionality completely, please read the following preview section for details.

Progress status: `see on GitHub <https://github.com/ECP-WarpX/WarpX/projects/10>`_

Introduction to CMake
=====================

If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_ from the HEP Software foundation is the perfect place to get started with it.

If you just want to use CMake to build the project, jump into sections *1. Introduction*, *2. Building with CMake* and *9. Finding Packages*.

Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`_ compiler: e.g. GCC 5, Clang 3.6 or newer
- `CMake 3.14.0+ <https://cmake.org>`_
- `AMReX <https://amrex-codes.github.io>`_: we automatically download and compile a copy of AMReX

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`_: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 9.0+ <https://developer.nvidia.com/cuda-downloads>`_: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`_: for threaded CPU execution (currently not fully accelerated)
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support
- `Boost 1.66.0+ <https://www.boost.org/>`_: for QED support
- `openPMD-api 0.12.0+ <https://github.com/openPMD/openPMD-api>`_: we automatically download and compile a copy of openPMD-api for openPMD I/O support

  - see `optional I/O backends <https://github.com/openPMD/openPMD-api#dependencies>`_
- `CCache <https://ccache.dev>`_: to speed up rebuilds (needs 3.7.9+ for CUDA)

Install Dependencies
====================

macOS/Linux:

.. code-block:: bash

   spack env create warpx-dev
   spack env activate warpx-dev
   spack add adios2
   spack add ccache
   spack add cmake
   spack add fftw
   spack add hdf5
   spack add mpi
   spack add pkgconfig  # for fftw
   # optional:
   # spack add cuda
   spack install

(in new terminals, re-activate the environment with ``spack env activate warpx-dev`` again)

or macOS/Linux:

.. code-block:: bash

   brew update
   brew install adios2
   brew install ccache
   brew install cmake
   brew install fftw
   brew install hdf5-mpi
   brew install libomp
   brew install pkg-config  # for fftw
   brew install open-mpi

Now, ``cmake --version`` should be at version 3.14.0 or newer.

Configure your compiler
=======================

For example, using a GCC on macOS:

.. code-block:: bash

   export CC=$(which gcc)
   export CXX=$(which g++)

If you also want to select a CUDA compiler:

.. code-block:: bash

   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which g++)

Build & Test
============

From the base of the WarpX source directory, execute:

.. code-block:: bash

   # find dependencies & configure
   cmake -S . -B build

   # build using up to four threads
   cmake --build build -j 4

   # run tests (todo)

You can inspect and modify build options after running ``cmake -S . -B build`` with either

.. code-block:: bash

   ccmake build

or by providing arguments to the CMake call: ``cmake -S . -B build -D<OPTION_A>=<VALUE_A> -D<OPTION_B>=<VALUE_B>``

============================= ============================================ =======================================================
CMake Option                  Default & Values                             Description
============================= ============================================ =======================================================
``CMAKE_BUILD_TYPE``          **RelWithDebInfo**/Release/Debug             Type of build, symbols & optimizations
``WarpX_APP``                 **ON**/OFF                                   Build the WarpX executable application
``WarpX_ASCENT``              ON/**OFF**                                   Ascent in situ visualization
``WarpX_COMPUTE``             NOACC/**OMP**/CUDA/DPCPP/HIP                 On-node, accelerated computing backend
``WarpX_DIMS``                **3**/2/RZ                                   Simulation dimensionality
``WarpX_LIB``                 ON/**OFF**                                   Build WarpX as a shared library
``WarpX_MPI``                 **ON**/OFF                                   Multi-node support (message-passing)
``WarpX_MPI_THREAD_MULTIPLE`` **ON**/OFF                                   MPI thread-multiple support, i.e. for ``async_io``
``WarpX_OPENPMD``             ON/**OFF**                                   openPMD I/O (HDF5, ADIOS)
``WarpX_PARSER_DEPTH``        **24**                                       Maximum parser depth for input file functions
``WarpX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``WarpX_PSATD``               ON/**OFF**                                   Spectral solver
``WarpX_QED``                 ON/**OFF**                                   PICSAR QED (requires Boost and PICSAR)
``WarpX_amrex_repo``          ``https://github.com/AMReX-Codes/amrex.git`` Repository URI to pull and build AMReX from
``WarpX_amrex_branch``        ``development``                              Repository branch for ``WarpX_amrex_repo``
``WarpX_amrex_internal``      **ON**/OFF                                   Needs a pre-installed AMReX library if set to ``OFF``
``WarpX_openpmd_internal``    **ON**/OFF                                   Needs a pre-installed openPMD library if set to ``OFF``
============================= ============================================ =======================================================

For example, one can also build against a local AMReX git repo.
Assuming AMReX' source is located in ``$HOME/src/amrex`` and changes are committed into a branch such as ``my-amrex-branch`` then pass to ``cmake`` the arguments: ``-DWarpX_amrex_repo=file://$HOME/src/amrex -DWarpX_amrex_branch=my-amrex-branch``.

For developers, WarpX can be configured in further detail with options from AMReX, which are `documented in the AMReX manual <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`_.

Run
===

An executable WarpX binary with the current compile-time options encoded in its file name will be created in ``build/bin/``.

Additionally, a `symbolic link <https://en.wikipedia.org/wiki/Symbolic_link>`_ named ``warpx`` can be found in that directory, which points to the last built WarpX executable.

Python Wrapper
==============

The Python wrapper library can be built by pre-building WarpX into one or more shared libraries.
For full functionality in 2D, 3D and RZ geometry, the following workflow can be executed:

.. code-block:: bash

   # build
   for d in 2 3 RZ; do
     cmake -S . -B build -DWarpX_DIMS=$d -DWarpX_LIB=ON
     cmake --build build -j 4
   done

   # package
   PYWARPX_LIB_DIR=$PWD/build/lib python3 -m pip wheel Python/

   # install
   python3 -m pip install pywarpx-*whl
