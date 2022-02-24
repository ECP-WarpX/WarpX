.. _building-cori:

Cori (NERSC)
============

The `Cori cluster <https://docs.nersc.gov/systems/cori/>`_ is located at NERSC.

If you are new to this system, please see the following resources:

* `GPU nodes <https://docs-dev.nersc.gov/cgpu/access>`__

* `Cori user guide <https://docs.nersc.gov/>`__
* Batch system: `Slurm <https://docs.nersc.gov/jobs/>`__
* `Jupyter service <https://docs.nersc.gov/services/jupyter/>`__
* `Production directories <https://www.nersc.gov/users/storage-and-file-systems/>`__:

  * ``$SCRATCH``: per-user production directory (20TB)
  * ``/global/cscratch1/sd/m3239``: shared production directory for users in the project ``m3239`` (50TB)
  * ``/global/cfs/cdirs/m3239/``: community file system for users in the project ``m3239`` (100TB)

Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

KNL
^^^

We use the following modules and environments on the system (``$HOME/knl_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/cori-nersc/knl_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/knl_warpx.profile.example``.

And install ADIOS2, BLAS++ and LAPACK++:

.. code-block:: bash

   source $HOME/knl_warpx.profile

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-knl-build
   cmake -S src/c-blosc -B src/c-blosc-knl-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/c-blosc-1.12.1-knl-install
   cmake --build src/c-blosc-knl-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-knl-build
   cmake -S src/adios2 -B src/adios2-knl-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/adios2-2.7.1-knl-install
   cmake --build src/adios2-knl-build --target install --parallel 16

   # BLAS++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/blaspp.git src/blaspp
   rm -rf src/blaspp-knl-build
   cmake -S src/blaspp -B src/blaspp-knl-build -Duse_openmp=ON -Duse_cmake_find_blas=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/blaspp-master-knl-install
   cmake --build src/blaspp-knl-build --target install --parallel 16

   # LAPACK++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/lapackpp.git src/lapackpp
   rm -rf src/lapackpp-knl-build
   CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S src/lapackpp -B src/lapackpp-knl-build -Duse_cmake_find_lapack=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DLAPACK_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/lapackpp-master-knl-install
   cmake --build src/lapackpp-knl-build --target install --parallel 16

For PICMI and Python workflows, also install a virtual environment:

.. code-block:: bash

   # establish Python dependencies
   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv

   python3 -m venv $HOME/sw/venvs/knl_warpx
   source $HOME/sw/venvs/knl_warpx/bin/activate

   python3 -m pip install --upgrade pip
   MPICC="cc -shared" python3 -m pip install -U --no-cache-dir -v mpi4py

Haswell
^^^^^^^

We use the following modules and environments on the system (``$HOME/haswell_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/cori-nersc/haswell_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/haswell_warpx.profile.example``.

And install ADIOS2, BLAS++ and LAPACK++:

.. code-block:: bash

   source $HOME/haswell_warpx.profile

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-haswell-build
   cmake -S src/c-blosc -B src/c-blosc-haswell-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/c-blosc-1.12.1-haswell-install
   cmake --build src/c-blosc-haswell-build --target install --parallel 16

   # ADIOS2
   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-haswell-build
   cmake -S src/adios2 -B src/adios2-haswell-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/adios2-2.7.1-haswell-install
   cmake --build src/adios2-haswell-build --target install --parallel 16

   # BLAS++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/blaspp.git src/blaspp
   rm -rf src/blaspp-haswell-build
   cmake -S src/blaspp -B src/blaspp-haswell-build -Duse_openmp=ON -Duse_cmake_find_blas=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/blaspp-master-haswell-install
   cmake --build src/blaspp-haswell-build --target install --parallel 16

   # LAPACK++ (for PSATD+RZ)
   git clone https://bitbucket.org/icl/lapackpp.git src/lapackpp
   rm -rf src/lapackpp-haswell-build
   CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S src/lapackpp -B src/lapackpp-haswell-build -Duse_cmake_find_lapack=ON -DBLAS_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DLAPACK_LIBRARIES=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$HOME/sw/lapackpp-master-haswell-install
   cmake --build src/lapackpp-haswell-build --target install --parallel 16

For PICMI and Python workflows, also install a virtual environment:

.. code-block:: bash

   # establish Python dependencies
   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv

   python3 -m venv $HOME/sw/venvs/haswell_warpx
   source $HOME/sw/venvs/haswell_warpx/bin/activate

   python3 -m pip install --upgrade pip
   MPICC="cc -shared" python3 -m pip install -U --no-cache-dir -v mpi4py

GPU (V100)
^^^^^^^^^^

Cori provides a partition with `18 nodes that include V100 (16 GB) GPUs <https://docs-dev.nersc.gov/cgpu/>`__.
We use the following modules and environments on the system (``$HOME/gpu_warpx.profile``).
You can copy this file from ``Tools/machines/cori-nersc/gpu_warpx.profile.example``:

.. literalinclude:: ../../../../Tools/machines/cori-nersc/gpu_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/gpu_warpx.profile.example``.

And install ADIOS2:

.. code-block:: bash

   source $HOME/gpu_warpx.profile

   # c-blosc (I/O compression)
   git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git src/c-blosc
   rm -rf src/c-blosc-gpu-build
   cmake -S src/c-blosc -B src/c-blosc-gpu-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/c-blosc-1.12.1-gpu-install
   cmake --build src/c-blosc-gpu-build --target install --parallel 16

   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git src/adios2
   rm -rf src/adios2-gpu-build
   cmake -S src/adios2 -B src/adios2-gpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=$HOME/sw/adios2-2.7.1-gpu-install
   cmake --build src/adios2-gpu-build --target install --parallel 16

For PICMI and Python workflows, also install a virtual environment:

.. code-block:: bash

   # establish Python dependencies
   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv

   python3 -m venv $HOME/sw/venvs/gpu_warpx
   source $HOME/sw/venvs/gpu_warpx/bin/activate

   python3 -m pip install --upgrade pip
   python3 -m pip install -U --no-cache-dir -v mpi4py

Building WarpX
--------------

We recommend to store the above lines in individual ``warpx.profile`` files, as suggested above.
If you want to run on either of the three partitions of Cori, open a new terminal, log into Cori and *source* the environment you want to work with:

.. code-block:: bash

   # KNL:
   source $HOME/knl_warpx.profile

   # Haswell:
   #source $HOME/haswell_warpx.profile

   # GPU:
   #source $HOME/gpu_warpx.profile

.. warning::

   Consider that all three Cori partitions are *incompatible*.

   Do not *source* multiple ``...warpx.profile`` files in the same terminal session.
   Open a new terminal and log into Cori again, if you want to switch the targeted Cori partition.

   If you re-submit an already compiled simulation that you ran on another day or in another session, *make sure to source* the corresponding ``...warpx.profile`` again after login!

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   #                       append if you target GPUs:    -DWarpX_COMPUTE=CUDA
   cmake -S . -B build -DWarpX_DIMS=3
   cmake --build build -j 16

The general :ref:`cmake compile-time options and instructions for Python (PICMI) bindings <building-cmake-python>` apply as usual:

.. code-block:: bash

   # PICMI build
   cd $HOME/src/warpx

   # compile parallel PICMI interfaces with openPMD support and 3D, 2D and RZ
   WARPX_MPI=ON BUILD_PARALLEL=16 python3 -m pip install --force-reinstall -v .

.. _running-cpp-cori:

Running
-------

Navigate (i.e. ``cd``) into one of the production directories (e.g. ``$SCRATCH``) before executing the instructions below.

KNL
^^^

The batch script below can be used to run a WarpX simulation on 2 KNL nodes on
the supercomputer Cori at NERSC. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<job name>`` could be ``laserWakefield``.

Do not forget to first ``source $HOME/knl_warpx.profile`` if you have not done so already for this terminal session.

For PICMI Python runs, the ``<path/to/executable>`` has to read ``python3`` and the ``<input file>`` is the path to your PICMI input script.

.. literalinclude:: ../../../../Tools/machines/cori-nersc/cori_knl.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/cori_knl.sbatch``.

To run a simulation, copy the lines above to a file ``cori_knl.sbatch`` and run

.. code-block:: bash

   sbatch cori_knl.sbatch

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Cori KNL for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=64`` and ``amr.blocking_factor=64`` so that the size of
  each grid is fixed to ``64**3`` (we are not using load-balancing here).

* **8 MPI ranks per KNL node**, with ``OMP_NUM_THREADS=8`` (that is 64 threads
  per KNL node, i.e. 1 thread per physical core, and 4 cores left to the
  system).

* **2 grids per MPI**, *i.e.*, 16 grids per KNL node.

Haswell
^^^^^^^

The batch script below can be used to run a WarpX simulation on 1 `Haswell node <https://docs.nersc.gov/systems/cori/>`_ on the supercomputer Cori at NERSC.

Do not forget to first ``source $HOME/haswell_warpx.profile`` if you have not done so already for this terminal session.

.. literalinclude:: ../../../../Tools/machines/cori-nersc/cori_haswell.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/cori_haswell.sbatch``.

To run a simulation, copy the lines above to a file ``cori_haswell.sbatch`` and
run

.. code-block:: bash

   sbatch cori_haswell.sbatch

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Cori Haswell for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* **4 MPI ranks per Haswell node** (2 MPI ranks per `Intel Xeon E5-2698 v3 <https://ark.intel.com/content/www/us/en/ark/products/81060/intel-xeon-processor-e5-2698-v3-40m-cache-2-30-ghz.html>`_), with ``OMP_NUM_THREADS=16`` (which uses `2x hyperthreading <https://docs.nersc.gov/jobs/affinity/>`_)

GPU (V100)
^^^^^^^^^^

Do not forget to first ``source $HOME/gpu_warpx.profile`` if you have not done so already for this terminal session.

Due to the limited amount of GPU development nodes, just request a single node with the above defined ``getNode`` function.
For single-node runs, try to run one grid per GPU.

A multi-node batch script template can be found below:

.. literalinclude:: ../../../../Tools/machines/cori-nersc/cori_gpu.sbatch
   :language: bash
   :caption: You can copy this file from ``Tools/machines/cori-nersc/cori_gpu.sbatch``.


.. _post-processing-cori:

Post-Processing
---------------

For post-processing, most users use Python via NERSC's `Jupyter service <https://jupyter.nersc.gov>`__ (`Docs <https://docs.nersc.gov/services/jupyter/>`__).

As a one-time preparatory setup, `create your own Conda environment as described in NERSC docs <https://docs.nersc.gov/services/jupyter/#conda-environments-as-kernels>`__.
In this manual, we often use this ``conda create`` line over the officially documented one:

.. code-block:: bash

   conda create -n myenv -c conda-forge python mamba ipykernel ipympl==0.8.6 matplotlib numpy pandas yt openpmd-viewer openpmd-api h5py fast-histogram

We then follow the `Customizing Kernels with a Helper Shell Script <https://docs.nersc.gov/services/jupyter/#customizing-kernels-with-a-helper-shell-script>`__ section to finalize the setup of using this conda-environment as a custom Jupyter kernel.

When opening a Jupyter notebook, just select the name you picked for your custom kernel on the top right of the notebook.

Additional software can be installed later on, e.g., in a Jupyter cell using ``!mamba install -c conda-forge ...``.
Software that is not available via conda can be installed via ``!python -m pip install ...``.

.. warning::

   Jan 6th, 2022 (NERSC-INC0179165 and `ipympl #416 <https://github.com/matplotlib/ipympl/issues/416>`__):
   Above, we fixated the ``ipympl`` version to *not* take the latest release of `Matplotlib Jupyter Widgets <https://github.com/matplotlib/ipympl>`__.
   This is an intentional work-around; the ``ipympl`` version needs to exactly fit the version pre-installed on the Jupyter base system.
