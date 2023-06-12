.. _building-fugaku:

Fugaku (Riken)
==============

The `Fugaku cluster <https://docs.nersc.gov/systems/perlmutter/>`_ is located at the Riken Center for Computational Science (Japan).


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Fugaku user guide <https://www.r-ccs.riken.jp/en/fugaku/user-guide/>`__


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx


Compiling WarpX on Fugaku is more pratical on a compute node. Use the following commands to acquire a compute node for one hour:

.. code-block:: bash

   pjsub --interact -L "elapse=02:00:00" -L "node=1" --sparam "wait-time=300" --mpi "max-proc-per-node=48" --all-mount-gfscache


Then, load ``cmake`` and ``ninja`` using ``spack``:

.. code-block:: bash

   . /vol0004/apps/oss/spack/share/spack/setup-env.sh
   spack load cmake@3.21.4%fj@4.8.1 arch=linux-rhel8-a64fx

   # optional: faster builds
   spack load ninja@1.11.1%fj@4.8.1

   # avoid harmless warning messages "[WARN] xos LPG [...]"
   export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH


At this point we need to download and compile the libraries required for OpenPMD support:

.. code-block:: bash
   export CC=$(which mpifcc)
   export CXX=$(which mpiFCC)
   export CFLAGS="-O3 -Nclang -Nlibomp -Klib -g -DNDEBUG"
   export CXXFLAGS="${CFLAGS}"

   export CMAKE_PREFIX_PATH=${HOME}/sw/a64fx-fj490/c-blosc-1.21.1-install:$CMAKE_PREFIX_PATH
   export CMAKE_PREFIX_PATH=${HOME}/sw/a64fx-fj490/adios2-2.8.3-install:$CMAKE_PREFIX_PATH

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-fg-build
   cmake -S src/c-blosc -B src/c-blosc-fg-build -DBUILD_SHARED_LIBS=OFF -DBUILD_SHARED=OFF -DBUILD_STATIC=ON -DBUILD_TESTS=OFF -DBUILD_FUZZERS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/a64fx-fj490/c-blosc-1.21.1-install
   cmake --build src/c-blosc-fg-build --target install --parallel 24

   # ADIOS2
   git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-fg-build
   cmake -S src/adios2 -B src/adios2-fg-build -DBUILD_SHARED_LIBS=OFF -DADIOS2_USE_Blosc=ON -DBUILD_TESTING=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/a64fx-fj490/adios2-2.8.3-install
   cmake --build src/adios2-fg-build --target install -j 24

Finally, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   export CC=$(which mpifcc)
   export CXX=$(which mpiFCC)
   export CFLAGS="-Nclang"
   export CXXFLAGS="-Nclang"

   cmake -S . -B build -DWarpX_COMPUTE=OMP \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_FLAGS_RELEASE="-Ofast -mllvm -polly -mllvm -polly-parallel" \
   -DAMReX_DIFFERENT_COMPILER=ON \
   -DWarpX_MPI_THREAD_MULTIPLE=OFF

   cmake --build build -j 48

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
A 3D WarpX executable is now in ``build/bin/`` and :ref:`can be run <running-cpp-fugaku>` with a :ref:`3D example inputs file <usage-examples>`.

.. _running-cpp-fugaku:

Running
-------

.. _running-cpp-fugaku-A64FX-CPUs:

A64FX CPUs
^^^^^^^^^^

In non-interactive runs, you can use `pjsub submit.sh` where `submit.sh` can be adapted from:

.. literalinclude:: ../../../../Tools/machines/fugaku-riken/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/fugaku-riken/submit.sh``.

Note: the ``Boost Eco Mode`` mode that is set in this example increases the default frequency of the A64FX
from 2 GHz to 2.2 GHz, while at the same time switching off one of the two floating-point arithmetic
pipelines. Some preliminary tests with WarpX show that this mode achieves performances similar to those of
the normal mode but with a reduction of the energy consumption of approximately 20%.
