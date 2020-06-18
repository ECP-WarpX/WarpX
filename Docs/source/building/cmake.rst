.. _building-cmake:

CMake Transition
================

We are currently transitioning to support `CMake <https://cmake.org>`_ as our primary build system.
Until we have transitioned our documentation and functionality completely, please read the following preview section for details.

Progress status: `see on GitHub <https://github.com/ECP-WarpX/WarpX/projects/10>`_

Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`_ compiler: e.g. GCC 5, Clang 3.6 or newer
- `CMake 3.14.0+ <https://cmake.org>`_
- `AMReX (development) <https://amrex-codes.github.io>`_: we automatically download and compile a copy of AMReX

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`_: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 9.0+ <https://developer.nvidia.com/cuda-downloads>`_: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`_: for threaded CPU execution (currently not fully accelerated)
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support
- `Boost 1.66.0+ <https://www.boost.org/>`_: for QED support
- `openPMD-api 0.11.1+ <https://github.com/openPMD/openPMD-api>`_: we automatically download and compile a copy of openPMD-api for openPMD I/O support

  - see `optional I/O backends <https://github.com/openPMD/openPMD-api#dependencies>`_
- `CCache <https://ccache.dev>`_: to speed up rebuilds (needs 3.7.9+ for CUDA)

Install Dependencies
====================

macOS/Linux:

.. code-block:: bash

   spack env create warpx-dev
   spack env activate warpx-dev
   spack add ccache
   spack add cmake
   spack add fftw
   spack add mpi
   spack add openpmd-api
   spack add pkgconfig  # for fftw
   # optional:
   # spack add cuda
   spack install

(in new terminals, re-activate the environment with ``spack env activate warpx-dev`` again)

or macOS/Linux:

.. code-block:: bash

   brew update
   brew install ccache
   brew install cmake
   brew install fftw
   brew install libomp
   brew install pkg-config  # for fftw
   brew install open-mpi
   brew install openpmd-api

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

   mkdir -p build
   cd build

   # find dependencies & configure
   cmake ..

   # build using up to four threads
   make -j 4

   # run tests (todo)

You can inspect and modify build options after running ``cmake ..`` with either

.. code-block:: bash

   ccmake .

or by providing arguments to the CMake call: ``cmake .. -D<OPTION_A>=<VALUE_A> -D<OPTION_B>=<VALUE_B>``

=========================== ============================================ =======================================================
CMake Option                Default & Values                             Description
=========================== ============================================ =======================================================
``CMAKE_BUILD_TYPE``        **RelWithDebInfo**/Release/Debug             Type of build, symbols & optimizations
``WarpX_ASCENT``            ON/**OFF**                                   Ascent in situ visualization
``WarpX_COMPUTE``           **NONE**/CUDA/OMP                            Parallel, on-node computing backend
``WarpX_DIMS``              **3**/2/RZ                                   Simulation dimensionality
``WarpX_OPENPMD``           ON/**OFF**                                   openPMD I/O
``WarpX_PRECISION``         **double**/single                            Floating point precision (single/double)
``WarpX_PSATD``             ON/**OFF**                                   Spectral solver
``WarpX_QED``               ON/**OFF**                                   PICSAR QED (requires Boost)
``WarpX_amrex_repo``        ``https://github.com/AMReX-Codes/amrex.git`` Repository URI to pull and build AMReX from
``WarpX_amrex_branch``      ``development``                              Repository branch for ``WarpX_amrex_repo``
``WarpX_amrex_internal``    **ON**/OFF                                   Needs a pre-installed AMReX library if set to ``OFF``
``WarpX_openpmd_internal``  **ON**/OFF                                   Needs a pre-installed openPMD library if set to ``OFF``
=========================== ============================================ =======================================================

For example, one can also build against a local AMReX git repo.
Assuming AMReX' source is located in ``$HOME/src/amrex`` and changes are committed into a branch such as ``my-amrex-branch`` then pass to ``cmake`` the arguments: ``-DWarpX_amrex_repo=file://$HOME/src/amrex -DWarpX_amrex_branch=my-amrex-branch``.

WarpX benefits from further standardized options in AMReX, which are `documented in detail in the AMReX manual <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`_.
A commonly needed option is:

=========================== ============================================ =======================================================
CMake Option                Default & Values                             Description
=========================== ============================================ =======================================================
``ENABLE_MPI``              **ON**/OFF                                   multi-node (MPI)
=========================== ============================================ =======================================================
