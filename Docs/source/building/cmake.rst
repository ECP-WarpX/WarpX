.. _building-cmake:

CMake Transition
================

We are currently transitioning to support `CMake <https://cmake.org>`_ as our primary build system.
Until we have transitioned our documentation and functionality completely, please read the following preview section for details.

Progress status: `see on GitHub <https://github.com/ECP-WarpX/WarpX/projects/10>`_


.. _building-cmake-intro:

Introduction to CMake
=====================

If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_ from the HEP Software foundation is the perfect place to get started with it.

If you just want to use CMake to build the project, jump into sections *1. Introduction*, *2. Building with CMake* and *9. Finding Packages*.


Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`_ compiler: e.g. GCC 5, Clang 3.6 or newer
- `CMake 3.15.0+ <https://cmake.org>`_
- `AMReX <https://amrex-codes.github.io>`_: we automatically download and compile a copy of AMReX
- `PICSAR <https://github.com/ECP-WarpX/picsar>`_: we automatically download and compile a copy of PICSAR

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`_: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 9.0+ <https://developer.nvidia.com/cuda-downloads>`_: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`_: for threaded CPU execution (currently not fully accelerated)
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support
- `Boost 1.66.0+ <https://www.boost.org/>`_: for QED lookup tables generation support
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

============================= ============================================ ========================================================
CMake Option                  Default & Values                             Description
============================= ============================================ ========================================================
``CMAKE_BUILD_TYPE``          **RelWithDebInfo**/Release/Debug             Type of build, symbols & optimizations
``WarpX_APP``                 **ON**/OFF                                   Build the WarpX executable application
``WarpX_ASCENT``              ON/**OFF**                                   Ascent in situ visualization
``WarpX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``WarpX_DIMS``                **3**/2/RZ                                   Simulation dimensionality
``WarpX_EB``                  ON/**OFF**                                   Embedded boundary support
``WarpX_LIB``                 ON/**OFF**                                   Build WarpX as a shared library
``WarpX_MPI``                 **ON**/OFF                                   Multi-node support (message-passing)
``WarpX_MPI_THREAD_MULTIPLE`` **ON**/OFF                                   MPI thread-multiple support, i.e. for ``async_io``
``WarpX_OPENPMD``             ON/**OFF**                                   openPMD I/O (HDF5, ADIOS)
``WarpX_PARSER_DEPTH``        **24**                                       Maximum parser depth for input file functions
``WarpX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``WarpX_PSATD``               ON/**OFF**                                   Spectral solver
``WarpX_QED``                 **ON**/OFF                                   QED support (requires PICSAR)
``WarpX_QED_TABLE_GEN``       ON/**OFF**                                   QED table generation support (requires PICSAR and Boost)
============================= ============================================ ========================================================

WarpX can be configured in further detail with options from AMReX, which are `documented in the AMReX manual <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`_.

**Developers** might be interested in additional options that control dependencies of WarpX.
By default, the most important dependencies of WarpX are automatically downloaded for convenience:

============================= ============================================== ===========================================================
CMake Option                  Default & Values                               Description
============================= ============================================== ===========================================================
``WarpX_amrex_src``           *None*                                         Path to AMReX source directory (preferred if set)
``WarpX_amrex_repo``          ``https://github.com/AMReX-Codes/amrex.git``   Repository URI to pull and build AMReX from
``WarpX_amrex_branch``        ``development``                                Repository branch for ``WarpX_amrex_repo``
``WarpX_amrex_internal``      **ON**/OFF                                     Needs a pre-installed AMReX library if set to ``OFF``
``WarpX_openpmd_src``         *None*                                         Path to openPMD-api source directory (preferred if set)
``WarpX_openpmd_repo``        ``https://github.com/openPMD/openPMD-api.git`` Repository URI to pull and build openPMD-api from
``WarpX_openpmd_branch``      ``0.13.2``                                     Repository branch for ``WarpX_openpmd_repo``
``WarpX_openpmd_internal``    **ON**/OFF                                     Needs a pre-installed openPMD-api library if set to ``OFF``
``WarpX_picsar_src``          *None*                                         Path to PICSAR source directory (preferred if set)
``WarpX_picsar_repo``         ``https://github.com/ECP-WarpX/picsar.git``    Repository URI to pull and build PICSAR from
``WarpX_picsar_branch``       ``development``                                Repository branch for ``WarpX_picsar_repo``
``WarpX_picsar_internal``     **ON**/OFF                                     Needs a pre-installed PICSAR library if set to ``OFF``
============================= ============================================== ===========================================================

For example, one can also build against a local AMReX copy.
Assuming AMReX' source is located in ``$HOME/src/amrex``, add the ``cmake`` argument ``-DWarpX_amrex_src=$HOME/src/amrex``.
Relative paths are also supported, e.g. ``-DWarpX_amrex_src=../amrex``.

Or build against an AMReX feature branch of a colleague.
Assuming your colleague pushed AMReX to ``https://github.com/WeiqunZhang/amrex/`` in a branch ``new-feature`` then pass to ``cmake`` the arguments: ``-DWarpX_amrex_repo=https://github.com/WeiqunZhang/amrex.git -DWarpX_amrex_branch=new-feature``.

You can speed up the install further if you pre-install these dependencies, e.g. with a package manager.
Set ``-DWarpX_<dependency-name>_internal=OFF`` and add installation prefix of the dependency to the environment variable `CMAKE_PREFIX_PATH <https://cmake.org/cmake/help/latest/envvar/CMAKE_PREFIX_PATH.html>`__.
Please see the :ref:`introduction to CMake <building-cmake-intro>` if this sounds new to you.


Run
===

An executable WarpX binary with the current compile-time options encoded in its file name will be created in ``build/bin/``.

Additionally, a `symbolic link <https://en.wikipedia.org/wiki/Symbolic_link>`_ named ``warpx`` can be found in that directory, which points to the last built WarpX executable.


Python Bindings
===============

Build and install ``pywarpx`` from the root of the WarpX source tree:

.. code-block:: bash

   python3 -m pip wheel -v .
   python3 -m pip install *whl

Environment variables can be used to control the build step:

============================= ============================================ ================================================================
Environment Variable          Default & Values                             Description
============================= ============================================ ================================================================
``WarpX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``WarpX_DIMS``                ``"2;3;RZ"``                                 Simulation dimensionalities (semicolon-separated list)
``WarpX_MPI``                 ON/**OFF**                                   Multi-node support (message-passing)
``WarpX_OPENPMD``             ON/**OFF**                                   openPMD I/O (HDF5, ADIOS)
``WarpX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``WarpX_PSATD``               ON/**OFF**                                   Spectral solver
``WarpX_QED``                 **ON**/OFF                                   PICSAR QED (requires PICSAR)
``WarpX_QED_TABLE_GEN``       ON/**OFF**                                   QED table generation (requires PICSAR and Boost)
``BUILD_PARALLEL``            ``2``                                        Number of threads to use for parallel builds
``BUILD_SHARED_LIBS``         ON/**OFF**                                   Build shared libraries for dependencies
``HDF5_USE_STATIC_LIBRARIES`` ON/**OFF**                                   Prefer static libraries for HDF5 dependency (openPMD)
``ADIOS_USE_STATIC_LIBS``     ON/**OFF**                                   Prefer static libraries for ADIOS1 dependency (openPMD)
``WarpX_amrex_src``           *None*                                       Absolute path to AMReX source directory (preferred if set)
``WarpX_amrex_internal``      **ON**/OFF                                   Needs a pre-installed AMReX library if set to ``OFF``
``WarpX_openpmd_src``         *None*                                       Absolute path to openPMD-api source directory (preferred if set)
``WarpX_openpmd_internal``    **ON**/OFF                                   Needs a pre-installed openPMD-api library if set to ``OFF``
``WarpX_picsar_src``          *None*                                       Absolute path to PICSAR source directory (preferred if set)
``WarpX_picsar_internal``     **ON**/OFF                                   Needs a pre-installed PICSAR library if set to ``OFF``
``PYWARPX_LIB_DIR``           *None*                                       If set, search for pre-built WarpX C++ libraries (see below)
============================= ============================================ ================================================================


Python Bindings (Developers)
============================

**Developers** might have all Python dependencies already installed and want to just overwrite an already installed ``pywarpx`` package:

.. code-block:: bash

   python3 -m pip install --force-reinstall -v .

Some Developers like to code directly against a local copy of AMReX, changing both code-bases at a time:

.. code-block:: bash

   WarpX_amrex_src=$PWD/../amrex python3 -m pip install --force-reinstall -v .

Additional environment control as common for CMake (:ref:`see above <building-cmake-intro>`) can be set as well, e.g. ``CC``, `CXX``, and ``CMAKE_PREFIX_PATH`` hints.
So another sophisticated example might be: use Clang as the compiler, build with local source copies of PICSAR and AMReX, support the PSATD solver, MPI and openPMD, hint a parallel HDF5 installation in ``$HOME/sw/hdf5-parallel-1.10.4``, and only build 3D geometry:

.. code-block:: bash

   CC=$(which clang) CXX=$(which clang++) WarpX_amrex_src=$PWD/../amrex WarpX_picsar_src=$PWD/../picsar WarpX_PSATD=ON WarpX_MPI=ON WarpX_OPENPMD=ON WarpX_DIMS=3 CMAKE_PREFIX_PATH=$HOME/sw/hdf5-parallel-1.10.4:$CMAKE_PREFIX_PATH python3 -m pip install --force-reinstall -v .

Here we wrote this all in one line, but one can also set all environment variables in a development environment and keep the pip call nice and short as in the beginning.
Note that you need to use absolute paths for external source trees, because pip builds in a temporary directory, e.g. ``export WarpX_amrex_src=$HOME/src/amrex``.

The Python library ``pywarpx`` can also be created by pre-building WarpX into one or more shared libraries externally.
For example, a package manager might split WarpX into a C++ package and a Python package.
If the C++ libraries are already pre-compiled, we can pick them up in the Python build step instead of compiling them again:

.. code-block:: bash

   # build WarpX executables and libraries
   for d in 2 3 RZ; do
     cmake -S . -B build -DWarpX_DIMS=$d -DWarpX_LIB=ON
     cmake --build build -j 4
   done

   # Python package
   PYWARPX_LIB_DIR=$PWD/build/lib python3 -m pip wheel .

   # install
   python3 -m pip install pywarpx-*whl


WarpX release managers might also want to generate a self-contained source package that can be distributed to exotic architectures:

.. code-block:: bash

   python setup.py sdist --dist-dir .
   python3 -m pip wheel -v pywarpx-*.tar.gz
   python3 -m pip install *whl

The above steps can also be executed in one go to build from source on a machine:

.. code-block:: bash

   python3 setup.py sdist --dist-dir .
   python3 -m pip install -v pywarpx-*.tar.gz

Last but not least, you can uninstall ``pywarpx`` as usual with:

   python3 -m pip uninstall pywarpx
