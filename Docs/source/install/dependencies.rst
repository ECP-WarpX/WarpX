.. _install-dependencies:

Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++17 <https://en.wikipedia.org/wiki/C%2B%2B17>`__ compiler, e.g., GCC 8, Clang 7, NVCC 11.0, MSVC 19.15 or newer
- `CMake 3.20.0+ <https://cmake.org>`__
- `Git 2.18+ <https://git-scm.com>`__
- `AMReX <https://amrex-codes.github.io>`__: we automatically download and compile a copy of AMReX
- `PICSAR <https://github.com/ECP-WarpX/picsar>`__: we automatically download and compile a copy of PICSAR

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`__: for multi-node and/or multi-GPU execution
- for on-node accelerated compute *one of either*:

  - `OpenMP 3.1+ <https://www.openmp.org>`__: for threaded CPU execution or
  - `CUDA Toolkit 11.0+ (11.3+ recommended) <https://developer.nvidia.com/cuda-downloads>`__: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_) or
  - `ROCm 5.2+ (5.5+ recommended) <https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-rocm-installation-readme/>`__: for AMD GPU support
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support when running on CPU or SYCL

  - also needs the ``pkg-config`` tool on Unix
- `BLAS++ <https://github.com/icl-utk-edu/blaspp>`_ and `LAPACK++ <https://github.com/icl-utk-edu/lapackpp>`_: for spectral solver (PSATD) support in RZ geometry
- `Boost 1.66.0+ <https://www.boost.org/>`__: for QED lookup tables generation support
- `openPMD-api 0.15.1+ <https://github.com/openPMD/openPMD-api>`__: we automatically download and compile a copy of openPMD-api for openPMD I/O support

  - see `optional I/O backends <https://github.com/openPMD/openPMD-api#dependencies>`__, i.e., ADIOS2 and/or HDF5
- `Ascent 0.8.0+ <https://ascent.readthedocs.io>`__: for in situ 3D visualization
- `SENSEI 4.0.0+ <https://sensei-insitu.org>`__: for in situ analysis and visualization
- `CCache <https://ccache.dev>`__: to speed up rebuilds (For CUDA support, needs version 3.7.9+ and 4.2+ is recommended)
- `Ninja <https://ninja-build.org>`__: for faster parallel compiles
- `Python 3.7+ <https://www.python.org>`__

  - `mpi4py <https://mpi4py.readthedocs.io>`__
  - `numpy <https://numpy.org>`__
  - `periodictable <https://periodictable.readthedocs.io>`__
  - `picmistandard <https://picmi-standard.github.io>`__
  - `lasy <https://lasydoc.readthedocs.io>`__
  - see our ``requirements.txt`` file for compatible versions


Install
-------

Pick *one* of the installation methods below to install all dependencies for WarpX development in a consistent manner.

Conda (Linux/macOS/Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. tip::

   We recommend to configure your conda to use the faster `libmamba` `dependency solver <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`__.

   .. code-block:: bash

      conda update -n base conda
      conda install -n base conda-libmamba-solver
      conda config --set solver libmamba

With MPI (only Linux/macOS):

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge blaspp boost ccache cmake compilers git lapackpp "openpmd-api=*=mpi_mpich*" python numpy pandas scipy yt "fftw=*=mpi_mpich*" pkg-config matplotlib mamba ninja mpich pip virtualenv
   source activate warpx-dev

Without MPI:

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge blaspp boost ccache cmake compilers git lapackpp openpmd-api python numpy pandas scipy yt fftw pkg-config matplotlib mamba ninja pip virtualenv
   source activate warpx-dev

   # compile WarpX with -DWarpX_MPI=OFF

For legacy ``GNUmake`` builds, after each ``source activate warpx-dev``, you also need to set:

.. code-block:: bash

   export FFTW_HOME=${CONDA_PREFIX}
   export BLASPP_HOME=${CONDA_PREFIX}
   export LAPACKPP_HOME=${CONDA_PREFIX}

.. tip::

   A general option to deactivate that conda self-activates its base environment.
   This `avoids interference with the system and other package managers <https://collegeville.github.io/CW20/WorkshopResources/WhitePapers/huebl-working-with-multiple-pkg-mgrs.pdf>`__.

   .. code-block:: bash

      conda config --set auto_activate_base false


Spack (macOS/Linux)
^^^^^^^^^^^^^^^^^^^

First, download a `WarpX Spack desktop development environment <https://github.com/ECP-WarpX/WarpX/blob/development/Tools/machines/desktop>`__ of your choice.
For most desktop developments, pick the OpenMP environment for CPUs unless you have a supported GPU.

* **Debian/Ubuntu** Linux:

  * OpenMP: ``system=ubuntu; compute=openmp`` (CPUs)
  * CUDA: ``system=ubuntu; compute=cuda`` (Nvidia GPUs)
  * ROCm: ``system=ubuntu; compute=rocm`` (AMD GPUs)
  * SYCL: *todo* (Intel GPUs)
* **macOS**: first, prepare with ``brew install gpg2; brew install gcc``

  * OpenMP: ``system=macos; compute=openmp``

If you already `installed Spack <https://spack.io>`__, we recommend to activate its `binary caches <https://spack.io/spack-binary-packages/>`__ for faster builds:

.. code-block:: bash

   spack mirror add rolling https://binaries.spack.io/develop
   spack buildcache keys --install --trust

Now install the WarpX dependencies in a new WarpX development environment:

.. code-block:: bash

   # download environment file
   curl -sLO https://raw.githubusercontent.com/ECP-WarpX/WarpX/development/Tools/machines/desktop/spack-${system}-${compute}.yaml

   # create new development environment
   spack env create warpx-${compute}-dev spack-${system}-${compute}.yaml
   spack env activate warpx-${compute}-dev

   # installation
   spack install
   python3 -m pip install jupyter matplotlib numpy openpmd-api openpmd-viewer pandas scipy virtualenv yt

In new terminal sessions, re-activate the environment with

.. code-block:: bash

   spack env activate warpx-openmp-dev

again.
Replace ``openmp`` with the equivalent you chose.

For legacy ``GNUmake`` builds, after each ``source activate warpx-openmp-dev``, you also need to set:

.. code-block:: bash

   export FFTW_HOME=${SPACK_ENV}/.spack-env/view
   export BLASPP_HOME=${SPACK_ENV}/.spack-env/view
   export LAPACKPP_HOME=${SPACK_ENV}/.spack-env/view


Brew (macOS/Linux)
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   brew update
   brew tap openpmd/openpmd
   brew install adios2      # for openPMD
   brew install ccache
   brew install cmake
   brew install fftw        # for PSATD
   brew install git
   brew install hdf5-mpi    # for openPMD
   brew install libomp
   brew unlink gcc
   brew link --force libomp
   brew install pkg-config  # for fftw
   brew install open-mpi
   brew install openblas    # for PSATD in RZ
   brew install openpmd-api # for openPMD

If you also want to compile with PSATD in RZ, you need to manually install BLAS++ and LAPACK++:

.. code-block:: bash

   sudo mkdir -p /usr/local/bin/
   sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
   sudo chmod a+x /usr/local/bin/cmake-easyinstall

   cmake-easyinstall --prefix=/usr/local git+https://github.com/icl-utk-edu/blaspp.git \
       -Duse_openmp=OFF -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON
   cmake-easyinstall --prefix=/usr/local git+https://github.com/icl-utk-edu/lapackpp.git \
       -Duse_cmake_find_lapack=ON -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON


Apt (Debian/Ubuntu)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-pandas python3-pip python3-scipy python3-venv

   # optional:
   # for CUDA, either install
   #   https://developer.nvidia.com/cuda-downloads (preferred)
   # or, if your Debian/Ubuntu is new enough, use the packages
   #   sudo apt install nvidia-cuda-dev libcub-dev
