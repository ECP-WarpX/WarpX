.. _building-perlmutter:

Perlmutter (NERSC)
==================

.. warning::

   Perlmutter is still in acceptance testing and environments change often.
   Please reach visit this page often for updates and reach out to us if something needs an update.

The `Perlmutter cluster <https://docs.nersc.gov/systems/perlmutter/>`_ is located at NERSC.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `NERSC user guide <https://docs.nersc.gov/>`__
* Batch system: `Slurm <https://docs.nersc.gov/systems/perlmutter/#running-jobs>`__
* `Jupyter service <https://docs.nersc.gov/services/jupyter/>`__
* `Production directories <https://docs.nersc.gov/filesystems/perlmutter-scratch/>`__:

  * ``$PSCRATCH``: per-user production directory, purged every 30 days (<TBD>TB)
  * ``/global/cscratch1/sd/m3239``: shared production directory for users in the project ``m3239``, purged every 30 days (50TB)
  * ``/global/cfs/cdirs/m3239/``: community file system for users in the project ``m3239`` (100TB)


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/perlmutter_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/perlmutter-nersc/perlmutter_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/perlmutter-nersc/perlmutter_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/perlmutter_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/perlmutter_warpx.profile

And since Perlmutter does not yet provide a module for them, install ADIOS2, BLAS++ and LAPACK++:

.. code-block:: bash

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-pm-build
   cmake -S src/c-blosc -B src/c-blosc-pm-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/c-blosc-1.21.1
   cmake --build src/c-blosc-pm-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-pm-build
   cmake -S src/adios2 -B src/adios2-pm-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/adios2-2.7.1
   cmake --build src/adios2-pm-build --target install -j 16

   # BLAS++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/blaspp.git src/blaspp
   rm -rf src/blaspp-pm-build
   CXX=$(which CC) cmake -S src/blaspp -B src/blaspp-pm-build -Duse_openmp=ON -Dgpu_backend=CUDA -Duse_cmake_find_blas=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/blaspp-master
   cmake --build src/blaspp-pm-build --target install --parallel 16

   # LAPACK++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/lapackpp.git src/lapackpp
   rm -rf src/lapackpp-pm-build
   CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S src/lapackpp -B src/lapackpp-pm-build -Duse_cmake_find_lapack=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DLAPACK_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$HOME/sw/perlmutter/lapackpp-master
   cmake --build src/lapackpp-pm-build --target install --parallel 16

Optionally, download and install Python packages for :ref:`PICMI <usage-picmi>` or dynamic ensemble optimizations (:ref:`libEnsemble <libensemble>`):

.. code-block:: bash

   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv
   python3 -m pip cache purge
   rm -rf $HOME/sw/perlmutter/venvs/warpx
   python3 -m venv $HOME/sw/perlmutter/venvs/warpx
   source $HOME/sw/perlmutter/venvs/warpx/bin/activate
   python3 -m pip install --upgrade pip
   python3 -m pip install --upgrade wheel
   python3 -m pip install --upgrade cython
   python3 -m pip install --upgrade numpy
   python3 -m pip install --upgrade pandas
   python3 -m pip install --upgrade scipy
   MPICC="cc -target-accel=nvidia80 -shared" python3 -m pip install --upgrade mpi4py --no-build-isolation --no-binary mpi4py
   python3 -m pip install --upgrade openpmd-api
   python3 -m pip install --upgrade matplotlib
   python3 -m pip install --upgrade yt
   # optional: for libEnsemble
   python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
   cmake --build build -j 16

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

For a *full PICMI install*, follow the :ref:`instructions for Python (PICMI) bindings <building-cmake-python>`:

.. code-block:: bash

   # PICMI build
   cd $HOME/src/warpx

   # install or update dependencies
   python3 -m pip install -r requirements.txt

   # compile parallel PICMI interfaces in 3D, 2D, 1D and RZ
   WARPX_MPI=ON WARPX_COMPUTE=CUDA WARPX_PSATD=ON BUILD_PARALLEL=16 python3 -m pip install --force-reinstall --no-deps -v .

Or, if you are *developing*, do a quick PICMI install of a *single geometry* (see: :ref:`WarpX_DIMS <building-cmake-options>`) using:

.. code-block:: bash

   # find dependencies & configure
   cmake -S . -B build -DWarpX_COMPUTE=CUDA -DWarpX_PSATD=ON -DWarpX_LIB=ON -DWarpX_DIMS=RZ

   # build and then call "python3 -m pip install ..."
   cmake --build build --target pip_install -j 16


.. _running-cpp-perlmutter:

Running
-------

.. _running-cpp-perlmutter-A100-GPUs:

A100 GPUs (40 GB)
^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on multiple nodes (change ``-N`` accordingly) on the supercomputer Perlmutter at NERSC.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that we run one MPI rank per GPU.


.. literalinclude:: ../../../../Tools/machines/perlmutter-nersc/perlmutter.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/perlmutter-nersc/perlmutter.sbatch``.

To run a simulation, copy the lines above to a file ``perlmutter.sbatch`` and run

.. code-block:: bash

   sbatch perlmutter.sbatch

to submit the job.


.. _post-processing-perlmutter:

Post-Processing
---------------

For post-processing, most users use Python via NERSC's `Jupyter service <https://jupyter.nersc.gov>`__ (`Docs <https://docs.nersc.gov/services/jupyter/>`__).

Please follow the same process as for :ref:`NERSC Cori post-processing <post-processing-cori>`.
**Important:** The *environment + Jupyter kernel* must separate from the one you create for Cori.

The Perlmutter ``$PSCRATCH`` filesystem is currently not yet available on Jupyter.
Thus, store or copy your data to Cori's ``$SCRATCH`` or use the Community FileSystem (CFS) for now.
