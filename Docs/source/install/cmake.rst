.. _install-developers:
.. _building-cmake:
.. _building-cmake-intro:

Developers
==========

`CMake <https://cmake.org>`_ is our primary build system.
If you are new to CMake, `this short tutorial <https://hsf-training.github.io/hsf-training-cmake-webpage/>`_ from the HEP Software foundation is the perfect place to get started.
If you just want to use CMake to build the project, jump into sections `1. Introduction <https://hsf-training.github.io/hsf-training-cmake-webpage/01-intro/index.html>`__, `2. Building with CMake <https://hsf-training.github.io/hsf-training-cmake-webpage/02-building/index.html>`__ and `9. Finding Packages <https://hsf-training.github.io/hsf-training-cmake-webpage/09-findingpackages/index.html>`__.


Dependencies
------------

Before you start, you will need a copy of the WarpX source code:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx
   cd $HOME/src/warpx

WarpX depends on popular third party software.

* On your development machine, :ref:`follow the instructions here <install-dependencies>`.
* If you are on an HPC machine, :ref:`follow the instructions here <install-hpc>`.

.. toctree::
   :hidden:

   dependencies


Compile
-------

From the base of the WarpX source directory, execute:

.. code-block:: bash

   # find dependencies & configure
   #   see additional options below, e.g.
   #                   -DCMAKE_INSTALL_PREFIX=$HOME/sw/warpx
   cmake -S . -B build

   # compile, here we use four threads
   cmake --build build -j 4

**That's it!**
A 3D WarpX binary is now in ``build/bin/`` and :ref:`can be run <usage_run>` with a :ref:`3D example inputs file <usage-examples>`.
Most people execute the binary directly or copy it out.

If you want to install the executables in a programmatic way, run this:

.. code-block:: bash

   # for default install paths, you will need administrator rights, e.g. with sudo:
   cmake --build build --target install

You can inspect and modify build options after running ``cmake -S . -B build`` with either

.. code-block:: bash

   ccmake build

or by adding arguments with ``-D<OPTION>=<VALUE>`` to the first CMake call.
For example, this builds WarpX all geometries and enables Nvidia GPU (CUDA) support:

.. code-block:: bash

   cmake -S . -B build -DWarpX_DIMS="1;2;RZ;3" -DWarpX_COMPUTE=CUDA


.. _building-cmake-options:

Build Options
-------------

============================= ============================================ =========================================================
CMake Option                  Default & Values                             Description
============================= ============================================ =========================================================
``CMAKE_BUILD_TYPE``          RelWithDebInfo/**Release**/Debug             `Type of build, symbols & optimizations <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`__
``CMAKE_INSTALL_PREFIX``      system-dependent path                        `Install path prefix <https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html>`__
``CMAKE_VERBOSE_MAKEFILE``    ON/**OFF**                                   `Print all compiler commands to the terminal during build <https://cmake.org/cmake/help/latest/variable/CMAKE_VERBOSE_MAKEFILE.html>`__
``PYINSTALLOPTIONS``                                                       Additional options for ``pip install``, e.g., ``-v --user``
``WarpX_APP``                 **ON**/OFF                                   Build the WarpX executable application
``WarpX_ASCENT``              ON/**OFF**                                   Ascent in situ visualization
``WarpX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``WarpX_DIMS``                **3**/2/1/RZ                                 Simulation dimensionality. Use ``"1;2;RZ;3"`` for all.
``WarpX_EB``                  ON/**OFF**                                   Embedded boundary support (not supported in RZ yet)
``WarpX_GPUCLOCK``            **ON**/OFF                                   Add GPU kernel timers (cost function, +4 registers/kernel)
``WarpX_IPO``                 ON/**OFF**                                   Compile WarpX with interprocedural optimization (aka LTO)
``WarpX_LIB``                 ON/**OFF**                                   Build WarpX as a shared library, e.g., for PICMI Python
``WarpX_MPI``                 **ON**/OFF                                   Multi-node support (message-passing)
``WarpX_MPI_THREAD_MULTIPLE`` **ON**/OFF                                   MPI thread-multiple support, i.e. for ``async_io``
``WarpX_OPENPMD``             **ON**/OFF                                   openPMD I/O (HDF5, ADIOS)
``WarpX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``WarpX_PARTICLE_PRECISION``  SINGLE/**DOUBLE**                            Particle floating point precision (single/double), defaults to WarpX_PRECISION value if not set
``WarpX_PSATD``               ON/**OFF**                                   Spectral solver
``WarpX_QED``                 **ON**/OFF                                   QED support (requires PICSAR)
``WarpX_QED_TABLE_GEN``       ON/**OFF**                                   QED table generation support (requires PICSAR and Boost)
``WarpX_SENSEI``              ON/**OFF**                                   SENSEI in situ visualization
============================= ============================================ =========================================================

WarpX can be configured in further detail with options from AMReX, which are documented in the AMReX manual:

* `general AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`__
* `GPU-specific options <https://amrex-codes.github.io/amrex/docs_html/GPU.html#building-gpu-support>`__.

**Developers** might be interested in additional options that control dependencies of WarpX.
By default, the most important dependencies of WarpX are automatically downloaded for convenience:

============================= ============================================== ===========================================================
CMake Option                  Default & Values                               Description
============================= ============================================== ===========================================================
``BUILD_SHARED_LIBS``         ON/**OFF**                                     `Build shared libraries for dependencies <https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html>`__
``CCACHE_PROGRAM``            First found ``ccache`` executable.             Set to ``-DCCACHE_PROGRAM=NO`` to disable CCache.
``AMReX_CUDA_PTX_VERBOSE``    ON/**OFF**                                     Print CUDA code generation statistics from ``ptxas``.
``WarpX_amrex_src``           *None*                                         Path to AMReX source directory (preferred if set)
``WarpX_amrex_repo``          ``https://github.com/AMReX-Codes/amrex.git``   Repository URI to pull and build AMReX from
``WarpX_amrex_branch``        *we set and maintain a compatible commit*      Repository branch for ``WarpX_amrex_repo``
``WarpX_amrex_internal``      **ON**/OFF                                     Needs a pre-installed AMReX library if set to ``OFF``
``WarpX_openpmd_src``         *None*                                         Path to openPMD-api source directory (preferred if set)
``WarpX_openpmd_repo``        ``https://github.com/openPMD/openPMD-api.git`` Repository URI to pull and build openPMD-api from
``WarpX_openpmd_branch``      ``0.15.1``                                     Repository branch for ``WarpX_openpmd_repo``
``WarpX_openpmd_internal``    **ON**/OFF                                     Needs a pre-installed openPMD-api library if set to ``OFF``
``WarpX_picsar_src``          *None*                                         Path to PICSAR source directory (preferred if set)
``WarpX_picsar_repo``         ``https://github.com/ECP-WarpX/picsar.git``    Repository URI to pull and build PICSAR from
``WarpX_picsar_branch``       *we set and maintain a compatible commit*      Repository branch for ``WarpX_picsar_repo``
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

If you re-compile often, consider installing the `Ninja <https://github.com/ninja-build/ninja/wiki/Pre-built-Ninja-packages>`__ build system.
Pass ``-G Ninja`` to the CMake configuration call to speed up parallel compiles.


.. _building-cmake-envvars:

Configure your compiler
-----------------------

If you don't want to use your default compiler, you can set the following environment variables.
For example, using a Clang/LLVM:

.. code-block:: bash

   export CC=$(which clang)
   export CXX=$(which clang++)

If you also want to select a CUDA compiler:

.. code-block:: bash

   export CUDACXX=$(which nvcc)
   export CUDAHOSTCXX=$(which clang++)

We also support adding `additional compiler flags via environment variables <https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#cmake-language-environment-variables>`__ such as `CXXFLAGS <https://cmake.org/cmake/help/latest/envvar/CXXFLAGS.html>`__/`LDFLAGS <https://cmake.org/cmake/help/latest/envvar/LDFLAGS.html>`__:

.. code-block:: bash

   # example: treat all compiler warnings as errors
   export CXXFLAGS="-Werror"

.. note::

   Please clean your build directory with ``rm -rf build/`` after changing the compiler.
   Now call ``cmake -S . -B build`` (+ further options) again to re-initialize the build configuration.


Run
---

An executable WarpX binary with the current compile-time options encoded in its file name will be created in ``build/bin/``.
Note that you need separate binaries to run 1D, 2D, 3D, and RZ geometry inputs scripts.
Additionally, a `symbolic link <https://en.wikipedia.org/wiki/Symbolic_link>`__ named ``warpx`` can be found in that directory, which points to the last built WarpX executable.

More details on running simulations are in the section :ref:`Run WarpX <usage_run>`.
Alternatively, read on and also build our PICMI Python interface.


.. _building-cmake-python:

PICMI Python Bindings
---------------------

.. note::

   Preparation: make sure you work with up-to-date Python tooling.

   .. code-block:: bash

      python3 -m pip install -U pip setuptools wheel
      python3 -m pip install -U cmake

For PICMI Python bindings, configure WarpX to produce a library and call our ``pip_install`` *CMake target*:

.. code-block:: bash

   # find dependencies & configure for all WarpX dimensionalities
   cmake -S . -B build -DWarpX_DIMS="1;2;RZ;3" -DWarpX_LIB=ON

   # build and then call "python3 -m pip install ..."
   cmake --build build --target pip_install -j 4

**That's it!**
You can now :ref:`run a first 3D PICMI script <usage-picmi>` from our :ref:`examples <usage-examples>`.

Developers could now change the WarpX source code and then call the build line again to refresh the Python installation.

.. tip::

   If you do *not* develop with :ref:`a user-level package manager <install-dependencies>`, e.g., because you rely on a HPC system's environment modules, then consider to set up a virtual environment via `Python venv <https://docs.python.org/3/library/venv.html>`__.
   Otherwise, without a virtual environment, you likely need to add the CMake option ``-DPYINSTALLOPTIONS="--user"``.


.. _building-pip-python:

Python Bindings (Package Management)
------------------------------------

This section is relevant for Python package management, mainly for maintainers or people that rather like to interact only with ``pip``.

One can build and install ``pywarpx`` from the root of the WarpX source tree:

.. code-block:: bash

   python3 -m pip wheel -v .
   python3 -m pip install pywarpx*whl

This will call the CMake logic above implicitly.
Using this workflow has the advantage that it can build and package up multiple libraries with varying ``WarpX_DIMS`` into one ``pywarpx`` package.

Environment variables can be used to control the build step:

============================= ============================================ ================================================================
Environment Variable          Default & Values                             Description
============================= ============================================ ================================================================
``WARPX_COMPUTE``             NOACC/**OMP**/CUDA/SYCL/HIP                  On-node, accelerated computing backend
``WARPX_DIMS``                ``"1;2;3;RZ"``                               Simulation dimensionalities (semicolon-separated list)
``WARPX_EB``                  ON/**OFF**                                   Embedded boundary support (not supported in RZ yet)
``WARPX_MPI``                 ON/**OFF**                                   Multi-node support (message-passing)
``WARPX_OPENPMD``             **ON**/OFF                                   openPMD I/O (HDF5, ADIOS)
``WARPX_PRECISION``           SINGLE/**DOUBLE**                            Floating point precision (single/double)
``WARPX_PARTICLE_PRECISION``  SINGLE/**DOUBLE**                            Particle floating point precision (single/double), defaults to WarpX_PRECISION value if not set
``WARPX_PSATD``               ON/**OFF**                                   Spectral solver
``WARPX_QED``                 **ON**/OFF                                   PICSAR QED (requires PICSAR)
``WARPX_QED_TABLE_GEN``       ON/**OFF**                                   QED table generation (requires PICSAR and Boost)
``BUILD_PARALLEL``            ``2``                                        Number of threads to use for parallel builds
``BUILD_SHARED_LIBS``         ON/**OFF**                                   Build shared libraries for dependencies
``HDF5_USE_STATIC_LIBRARIES`` ON/**OFF**                                   Prefer static libraries for HDF5 dependency (openPMD)
``ADIOS_USE_STATIC_LIBS``     ON/**OFF**                                   Prefer static libraries for ADIOS1 dependency (openPMD)
``WARPX_AMREX_SRC``           *None*                                       Absolute path to AMReX source directory (preferred if set)
``WARPX_AMREX_REPO``          *None (uses cmake default)*                  Repository URI to pull and build AMReX from
``WARPX_AMREX_BRANCH``        *None (uses cmake default)*                  Repository branch for ``WARPX_AMREX_REPO``
``WARPX_AMREX_INTERNAL``      **ON**/OFF                                   Needs a pre-installed AMReX library if set to ``OFF``
``WARPX_OPENPMD_SRC``         *None*                                       Absolute path to openPMD-api source directory (preferred if set)
``WARPX_OPENPMD_INTERNAL``    **ON**/OFF                                   Needs a pre-installed openPMD-api library if set to ``OFF``
``WARPX_PICSAR_SRC``          *None*                                       Absolute path to PICSAR source directory (preferred if set)
``WARPX_PICSAR_INTERNAL``     **ON**/OFF                                   Needs a pre-installed PICSAR library if set to ``OFF``
``WARPX_CCACHE_PROGRAM``      First found ``ccache`` executable.           Set to ``NO`` to disable CCache.
``PYWARPX_LIB_DIR``           *None*                                       If set, search for pre-built WarpX C++ libraries (see below)
============================= ============================================ ================================================================

Note that we currently change the ``WARPX_MPI`` default intentionally to ``OFF``, to simplify a first install from source.

Some hints and workflows follow.
Developers, that want to test a change of the source code but did not change the ``pywarpx`` version number, can force a reinstall via:

.. code-block:: bash

   python3 -m pip install --force-reinstall --no-deps -v .

Some Developers like to code directly against a local copy of AMReX, changing both code-bases at a time:

.. code-block:: bash

   WARPX_AMREX_SRC=$PWD/../amrex python3 -m pip install --force-reinstall --no-deps -v .

Additional environment control as common for CMake (:ref:`see above <building-cmake-intro>`) can be set as well, e.g. ``CC``, `CXX``, and ``CMAKE_PREFIX_PATH`` hints.
So another sophisticated example might be: use Clang as the compiler, build with local source copies of PICSAR and AMReX, support the PSATD solver, MPI and openPMD, hint a parallel HDF5 installation in ``$HOME/sw/hdf5-parallel-1.10.4``, and only build 2D and 3D geometry:

.. code-block:: bash

   CC=$(which clang) CXX=$(which clang++) WARPX_AMREX_SRC=$PWD/../amrex WARPX_PICSAR_SRC=$PWD/../picsar WARPX_PSATD=ON WARPX_MPI=ON WARPX_DIMS="2;3" CMAKE_PREFIX_PATH=$HOME/sw/hdf5-parallel-1.10.4:$CMAKE_PREFIX_PATH python3 -m pip install --force-reinstall --no-deps -v .

Here we wrote this all in one line, but one can also set all environment variables in a development environment and keep the pip call nice and short as in the beginning.
Note that you need to use absolute paths for external source trees, because pip builds in a temporary directory, e.g. ``export WARPX_AMREX_SRC=$HOME/src/amrex``.

The Python library ``pywarpx`` can also be created by pre-building WarpX into one or more shared libraries externally.
For example, a package manager might split WarpX into a C++ package and a Python package.
If the C++ libraries are already pre-compiled, we can pick them up in the Python build step instead of compiling them again:

.. code-block:: bash

   # build WarpX executables and libraries
   cmake -S . -B build -DWarpX_DIMS="1;2;RZ;3" -DWarpX_LIB=ON
   cmake --build build -j 4

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

.. code-block:: bash

   python3 -m pip uninstall pywarpx
