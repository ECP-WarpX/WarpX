.. _install-dependencies:

Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`__ compiler: e.g. GCC 5, Clang 3.6 or newer
- `CMake 3.15.0+ <https://cmake.org>`__
- `Git 2.18+ <https://git-scm.com>`__
- `AMReX <https://amrex-codes.github.io>`__: we automatically download and compile a copy of AMReX
- `PICSAR <https://github.com/ECP-WarpX/picsar>`__: we automatically download and compile a copy of PICSAR

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`__: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 9.0+ <https://developer.nvidia.com/cuda-downloads>`__: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`__: for threaded CPU execution (currently not fully accelerated)
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support
- `BLAS++ <https://bitbucket.org/icl/blaspp>`_ and `LAPACK++ <https://bitbucket.org/icl/lapackpp>`_: for spectral solver (PSATD) support in RZ geometry
- `Boost 1.66.0+ <https://www.boost.org/>`__: for QED lookup tables generation support
- `openPMD-api 0.12.0+ <https://github.com/openPMD/openPMD-api>`__: we automatically download and compile a copy of openPMD-api for openPMD I/O support

  - see `optional I/O backends <https://github.com/openPMD/openPMD-api#dependencies>`__
- `CCache <https://ccache.dev>`__: to speed up rebuilds (needs 3.7.9+ for CUDA)
- `Ninja <https://ninja-build.org>`__: for faster parallel compiles


Install
-------

Pick *one* of the installation methods below to install all dependencies for WarpX development in a consistent manner.

Spack (macOS/Linux)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   spack env create warpx-dev
   spack env activate warpx-dev
   spack add adios2        # for openPMD
   spack add blaspp        # for PSATD in RZ
   spack add ccache
   spack add cmake
   spack add fftw          # for PSATD
   spack add hdf5          # for openPMD
   spack add lapackpp      # for PSATD in RZ
   spack add mpi
   spack add openpmd-api   # for openPMD
   spack add pkgconfig     # for fftw
   # optional:
   # spack add python
   # spack add py-pip
   # spack add cuda
   spack install

(in new terminals, re-activate the environment with ``spack env activate warpx-dev`` again)

If you also want to run runtime tests and added Python (``spack add python`` and ``spack add py-pip``) above, install also these additional Python packages in the active Spack environment:

.. code-block:: bash

   python -m pip install matplotlib==3.2.2 yt scipy numpy openpmd-api


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
   brew install pkg-config  # for fftw
   brew install open-mpi
   brew install openblas    # for PSATD in RZ
   brew install openpmd-api # for openPMD

If you also want to compile with PSATD in RZ, you need to manually install BLAS++ and LAPACK++:

.. code-block:: bash

   sudo mkdir -p /usr/local/bin/
   sudo curl -L -o /usr/local/bin/cmake-easyinstall https://git.io/JvLxY
   sudo chmod a+x /usr/local/bin/cmake-easyinstall

   cmake-easyinstall --prefix=/usr/local git+https://bitbucket.org/icl/blaspp.git \
       -Duse_openmp=OFF -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON
   cmake-easyinstall --prefix=/usr/local git+https://bitbucket.org/icl/lapackpp.git \
       -Duse_cmake_find_lapack=ON -Dbuild_tests=OFF -DCMAKE_VERBOSE_MAKEFILE=ON


Conda (Linux/macOS/Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Without MPI:

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge blaspp ccache cmake compilers git lapackpp openpmd-api python numpy scipy yt fftw matplotlib mamba ninja
   conda activate warpx-dev

   # compile WarpX with -DWarpX_MPI=OFF

With MPI (only Linux/macOS):

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge blaspp ccache cmake compilers git lapackpp openpmd-api=*=mpi_openmpi* python numpy scipy yt fftw=*=mpi_openmpi* matplotlib mamba ninja openmpi
   conda activate warpx-dev


Apt (Debian/Ubuntu)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-scipy

