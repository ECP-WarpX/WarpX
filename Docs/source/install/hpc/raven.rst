.. _building-raven:

Raven (MPCDF)
=============

The `Raven cluster <https://docs.mpcdf.gov/systems/raven/>`_ is located at MPCDF.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `MPCDF user guide <https://docs.mpcdf.mpg.de/doc/computing/raven-user-guide.html>`__
* Batch system: `Slurm <https://docs.mpcdf.mpg.de/doc/computing/raven-user-guide.html#slurm-batch-system>`__
* `Jupyter service <TODO>`__
* `Production directories <https://docs.mpcdf.mpg.de/doc/computing/raven-user-guide.html#file-systems>`__:

  * HOME ...
  * AFS ...
  * GPFS ...


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/raven_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/raven-mpcdf/raven_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/raven-mpcdf/raven_warpx.profile.example``.

We recommend to store the above lines in a file, such as ``$HOME/raven_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/raven_warpx.profile

And since raven does not yet provide a module for them, install ADIOS2, BLAS++ and LAPACK++:

.. code-block:: bash

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-pm-build
   cmake -S src/c-blosc -B src/c-blosc-pm-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/raven/c-blosc-1.21.1
   cmake --build src/c-blosc-pm-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-pm-build
   cmake -S src/adios2 -B src/adios2-pm-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/raven/adios2-2.7.1
   cmake --build src/adios2-pm-build --target install -j 16

   # BLAS++ (for PSATD+RZ)
   git clone https://github.com/icl-utk-edu/blaspp.git src/blaspp
   rm -rf src/blaspp-pm-build
   CXX=$(which CC) cmake -S src/blaspp -B src/blaspp-pm-build -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/raven/blaspp-master
   cmake --build src/blaspp-pm-build --target install --parallel 16

   # LAPACK++ (for PSATD+RZ)
   git clone https://github.com/icl-utk-edu/lapackpp.git src/lapackpp
   rm -rf src/lapackpp-pm-build
   CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S src/lapackpp -B src/lapackpp-pm-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=$HOME/sw/raven/lapackpp-master
   cmake --build src/lapackpp-pm-build --target install --parallel 16

Optionally, download and install Python packages for :ref:`PICMI <usage-picmi>` or dynamic ensemble optimizations (:ref:`libEnsemble <libensemble>`):

.. code-block:: bash

   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv
   python3 -m pip cache purge
   rm -rf $HOME/sw/raven/venvs/warpx
   python3 -m venv $HOME/sw/raven/venvs/warpx
   source $HOME/sw/raven/venvs/warpx/bin/activate
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


.. _running-cpp-raven:

Running
-------

.. _running-cpp-raven-A100-GPUs:

A100 GPUs (40 GB)
^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on multiple nodes (change ``-N`` accordingly) on the supercomputer raven at MPCDF.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that we run one MPI rank per GPU.


.. literalinclude:: ../../../../Tools/machines/raven-mpcdf/raven.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/raven-mpcdf/raven.sbatch``.

To run a simulation, copy the lines above to a file ``raven.sbatch`` and run

.. code-block:: bash

   sbatch raven.sbatch

to submit the job.


.. _post-processing-raven:

Post-Processing
---------------

.. note::

   TODO - any jupyter service?
