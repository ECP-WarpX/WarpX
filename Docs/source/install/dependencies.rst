.. _install-dependencies:

Dependencies
============

WarpX depends on the following popular third party software.
Please see installation instructions below.

- a mature `C++14 <https://en.wikipedia.org/wiki/C%2B%2B14>`__ compiler: e.g. GCC 5, Clang 3.6 or newer
- `CMake 3.15.0+ <https://cmake.org>`__
- `AMReX <https://amrex-codes.github.io>`__: we automatically download and compile a copy of AMReX
- `PICSAR <https://github.com/ECP-WarpX/picsar>`__: we automatically download and compile a copy of PICSAR

Optional dependencies include:

- `MPI 3.0+ <https://www.mpi-forum.org/docs/>`__: for multi-node and/or multi-GPU execution
- `CUDA Toolkit 9.0+ <https://developer.nvidia.com/cuda-downloads>`__: for Nvidia GPU support (see `matching host-compilers <https://gist.github.com/ax3l/9489132>`_)
- `OpenMP 3.1+ <https://www.openmp.org>`__: for threaded CPU execution (currently not fully accelerated)
- `FFTW3 <http://www.fftw.org>`_: for spectral solver (PSATD) support
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
   spack add adios2
   spack add ccache
   spack add cmake
   spack add fftw
   spack add hdf5
   spack add mpi
   spack add pkgconfig  # for fftw
   # optional:
   # spack add python
   # spack add py-pip
   # spack add cuda
   spack install

(in new terminals, re-activate the environment with ``spack env activate warpx-dev`` again)

If you also want to run runtime tests and added Python above, install also these additional Python packages in the active Spack environment:

.. code-block:: bash

   python -m pip install matplotlib==3.2.2 yt scipy numpy


Brew (macOS/Linux)
^^^^^^^^^^^^^^^^^^

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


Conda (Linux/macOS/Windows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Without MPI:

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge ccache cmake compilers git openpmd-api python numpy scipy yt fftw matplotlib mamba ninja
   conda activate warpx-dev

   # compile WarpX with -DWarpX_MPI=OFF

With MPI (only Linux/macOS):

.. code-block:: bash

   conda create -n warpx-dev -c conda-forge ccache cmake compilers git openpmd-api=*=mpi_openmpi* python numpy scipy yt fftw=*=mpi_openmpi* matplotlib mamba ninja openmpi
   conda activate warpx-dev


Apt (Debian/Ubuntu)
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential ccache cmake g++ libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-scipy

