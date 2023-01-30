.. _building-perlmutter:

Perlmutter (NERSC)
==================

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

We recommend to store the above lines in a file, such as ``$HOME/perlmutter_warpx.profile``.
Change the line that reads ``export proj="<yourProject>_g"  # change me`` to your NERSC project, e.g.,

.. code-block:: bash

   export proj="m3239_g"

and load it into your shell *after each login*:

.. code-block:: bash

   source $HOME/perlmutter_warpx.profile

*Only once on first setup*, we need to install packages not provided by NERSC:

.. code-block:: bash

   # this might take around 15min
   spack install

   python3 -m pip install h5py libensemble matplotlib nlopt openpmd-api openpmd-viewer pandas pytest scipy yt

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
