.. _building-cori:

Cori (NERSC)
============

The `Cori cluster <http://www.nersc.gov/users/computational-systems/cori>`_ is located at NERSC.

If you are new to this system, please see the following resources:

* `GPU nodes <https://docs-dev.nersc.gov/cgpu/access>`__

* `Cori user guide <https://docs.nersc.gov/>`_
* Batch system: `Slurm <https://docs.nersc.gov/jobs/>`_
* `Production directories <https://www.nersc.gov/users/storage-and-file-systems/>`_:

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

.. code-block:: bash

   module swap craype-haswell craype-mic-knl
   module swap PrgEnv-intel PrgEnv-gnu
   module load cmake/3.18.2
   module load cray-hdf5-parallel/1.10.5.2
   module load cray-fftw
   module load cray-python/3.7.3.2

   export CMAKE_PREFIX_PATH=$PWD/adios2-2.7.1-knl-install:$CMAKE_PREFIX_PATH

   export CXXFLAGS="-march=knl"
   export CFLAGS="-march=knl"

And install ADIOS2:

.. code-block:: bash

   source $HOME/knl_warpx.profile

   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git adios2
   cmake -S adios2 -B adios2-build -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=adios2-2.7.1-knl-install
   cmake --build adios2-build --target install --parallel 16

Haswell
^^^^^^^

We use the following modules and environments on the system (``$HOME/haswell_warpx.profile``).

.. code-block:: bash

   module load cmake/3.18.2
   module load cray-hdf5-parallel/1.10.5.2
   module load cray-fftw
   module load cray-python/3.7.3.2

   export CMAKE_PREFIX_PATH=$PWD/adios2-2.7.1-haswell-install:$CMAKE_PREFIX_PATH

And install ADIOS2:

.. code-block:: bash

   source $HOME/haswell_warpx.profile

   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git adios2
   cmake -S adios2 -B adios2-build -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=adios2-2.7.1-haswell-install
   cmake --build adios2-build --target install --parallel 16

GPU
^^^

We use the following modules and environments on the system (``$HOME/gpu_warpx.profile``).

.. code-block:: bash

   export proj="m1759"

   module purge
   module load cgpu gcc cuda cmake
   module load mvapich2
   # OpenMPI-UCX instead of mvapich:
   #module load openmpi/4.0.1-ucx-1.6

   export CMAKE_PREFIX_PATH=$PWD/adios2-2.7.1-gpu-install:$CMAKE_PREFIX_PATH

   # allocate a GPU, e.g. to compile on
   #   10 logical cores (5 physical), 1 GPU
   function getNode() {
       salloc -C gpu -N 1 -t 30 -c 10 --gres=gpu:1 -A $proj
   }

And install ADIOS2:

.. code-block:: bash

   source $HOME/gpu_warpx.profile

   git clone -b v2.7.1 https://github.com/ornladios/ADIOS2.git adios2
   cmake -S adios2 -B adios2-build -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DCMAKE_INSTALL_PREFIX=adios2-2.7.1-gpu-install
   cmake --build adios2-build --target install --parallel 16


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

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   #                           append if you target GPUs:    -DWarpX_COMPUTE=CUDA
   cmake -S build -B build -DWarpX_OPENPMD=ON -DWarpX_DIMS=3
   cmake --build build -j 16

The general :ref:`cmake compile-time options and instructions for Python (PICMI) bindings <building-cmake>` apply as usual.

.. warning::

   Consider that all three Cori partitions are *incompatible*.

   Do not *source* multiple ``...warpx.profile`` files in the same terminal session.
   Open a new terminal and log into Cori again, if you want to switch the targeted Cori partition.

   If you re-submit an already compiled simulation that you ran on another day or in another session, *make sure to source* the corresponding ``...warpx.profile`` again after login!

.. _running-cpp-cori:

Running
-------

KNL
^^^

The batch script below can be used to run a WarpX simulation on 2 KNL nodes on
the supercomputer Cori at NERSC. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<job name>`` could be ``laserWakefield``.

.. literalinclude:: ../../../../Tools/BatchScripts/batch_cori.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_cori.sh`` and run

.. code-block:: bash

   sbatch batch_cori.sh

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

.. literalinclude:: ../../../../Tools/BatchScripts/batch_cori_haswell.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_cori_haswell.sh`` and
run

.. code-block:: bash

   sbatch batch_cori_haswell.sh

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Cori Haswell for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* **4 MPI ranks per Haswell node** (2 MPI ranks per `Intel Xeon E5-2698 v3 <https://ark.intel.com/content/www/us/en/ark/products/81060/intel-xeon-processor-e5-2698-v3-40m-cache-2-30-ghz.html>`_), with ``OMP_NUM_THREADS=16`` (which uses `2x hyperthreading <https://docs.nersc.gov/jobs/affinity/>`_)

GPU
^^^

Due to the limited amount of GPU development nodes, just request a single node with the above defined ``getNode`` function.
For single-node runs, try to run one grid per GPU.
