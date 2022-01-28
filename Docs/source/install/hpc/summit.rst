.. _building-summit:

Summit (OLCF)
=============

The `Summit cluster <https://www.olcf.ornl.gov/summit/>`_ is located at OLCF.

If you are new to this system, please see the following resources:

* `Summit user guide <https://docs.olcf.ornl.gov/systems/summit_user_guide.html>`_
* Batch system: `LSF <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#running-jobs>`_
* `Jupyter service <https://jupyter.olcf.ornl.gov>`__
* `Production directories <https://docs.olcf.ornl.gov/data/storage_overview.html>`_:

  * ``$PROJWORK/$proj/``: shared with all members of a project (recommended)
  * ``$MEMBERWORK/$proj/``: single user (usually smaller quota)
  * ``$WORLDWORK/$proj/``: shared with all users
  * Note that the ``$HOME`` directory is mounted as read-only on compute nodes.
    That means you cannot run in your ``$HOME``.


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/summit_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/summit-olcf/summit_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/summit-olcf/summit_warpx.profile.example``.


We recommend to store the above lines in a file, such as ``$HOME/summit_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/summit_warpx.profile

Optionally, download and install Python packages for :ref:`PICMI <usage-picmi>` or dynamic ensemble optimizations (:ref:`libEnsemble <libensemble>`):

.. code-block:: bash

   python3 -m pip install --user --upgrade pip
   python3 -m pip install --user virtualenv
   python3 -m pip cache purge
   rm -rf $HOME/sw/venvs/warpx
   python3 -m venv $HOME/sw/venvs/warpx
   source $HOME/sw/venvs/warpx/bin/activate
   python3 -m pip install --upgrade pip
   python3 -m pip install --upgrade wheel
   python3 -m pip install --upgrade cython
   python3 -m pip install --upgrade numpy
   python3 -m pip install --upgrade pandas
   python3 -m pip install --upgrade scipy
   python3 -m pip install --upgrade mpi4py --no-binary mpi4py
   python3 -m pip install --upgrade openpmd-api
   python3 -m pip install --upgrade matplotlib==3.2.2  # does not try to build freetype itself
   python3 -m pip install --upgrade yt
   # WIP: issues with nlopt
   # python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=CUDA
   cmake --build build -j 6

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

For a full PICMI install, follow the :ref:`instructions for Python (PICMI) bindings <building-cmake-python>`.
We only prefix it to request a node for the compilation (``runNode``), so we can compile faster:

.. code-block:: bash

   # PICMI build
   cd $HOME/src/warpx

   # compile parallel PICMI interfaces with openPMD support and 3D, 2D and RZ
   runNode WARPX_MPI=ON WARPX_COMPUTE=CUDA WARPX_PSATD=ON BUILD_PARALLEL=32 python3 -m pip install --force-reinstall --no-deps -v .


.. _running-cpp-summit:

Running
-------

.. _running-cpp-summit-V100-GPUs:

V100 GPUs (16GB)
^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on
the supercomputer Summit at OLCF. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<input file>`` could be
``plasma_mirror_inputs``.
Note that WarpX runs with one MPI rank per GPU and there are 6 GPUs per node:

.. literalinclude:: ../../../../Tools/machines/summit-olcf/summit_v100.bsub
   :language: bash
   :caption: You can copy this file from ``Tools/machines/summit-olcf/summit_v100.bsub``.

To run a simulation, copy the lines above to a file ``summit_v100.bsub`` and
run
::

  bsub summit_v100.bsub

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Summit for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=256`` and ``amr.blocking_factor=128``.

* **One MPI rank per GPU** (e.g., 6 MPI ranks for the 6 GPUs on each Summit
  node)

* **Two `128x128x128` grids per GPU**, or **one `128x128x256` grid per GPU**.

A batch script with more options regarding profiling on Summit can be found at
:download:`Summit batch script <../../../../Tools/machines/summit-olcf/summit_profiling.bsub>`

.. _running-cpp-summit-Power9-CPUs:

Power9 CPUs
^^^^^^^^^^^

Similar to above, the batch script below can be used to run a WarpX simulation on
1 node on the supercomputer Summit at OLCF, on Power9 CPUs (i.e., the GPUs are
ignored).

.. literalinclude:: ../../../../Tools/machines/summit-olcf/summit_power9.bsub
   :language: bash
   :caption: You can copy this file from ``Tools/machines/summit-olcf/summit_power9.bsub``.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Summit for a well load-balanced problem, the following set of
parameters provided good performance:

* ``amr.max_grid_size=64`` and ``amr.blocking_factor=64``

* **Two MPI ranks per node** (i.e. 2 resource sets per node; equivalently, 1
  resource set per socket)

* **21 physical CPU cores per MPI rank**

* **21 OpenMP threads per MPI rank** (i.e. 1 OpenMP thread per physical core)

* **SMT 1 (Simultaneous Multithreading level 1)**

* **Sixteen `64x64x64` grids per MPI rank** (with default tiling in WarpX, this
  results in ~49 tiles per OpenMP thread)

.. _building-summit-io-performance:

I/O Performance Tuning
----------------------

.. _building-summit-large-blocks:

GPFS Large Block I/O
^^^^^^^^^^^^^^^^^^^^

Setting ``IBM_largeblock_io`` to ``true`` disables data shipping, saving overhead when writing/reading large contiguous I/O chunks.

.. code-block:: bash

   export IBM_largeblock_io=true

.. _building-summit-romio-hints:

ROMIO MPI-IO Hints
^^^^^^^^^^^^^^^^^^

You might notice some parallel HDF5 performance improvements on Summit by setting the appropriate ROMIO hints for MPI-IO operations.

.. code-block:: bash

   export OMPI_MCA_io=romio321
   export ROMIO_HINTS=./romio-hints

You can generate the ``romio-hints`` by issuing the following command. Remember to change the number of ``cb_nodes`` to match the number of compute nodes you are using (example here: ``64``).

.. code-block:: bash

   cat > romio-hints << EOL
   romio_cb_write enable
   romio_ds_write enable
   cb_buffer_size 16777216
   cb_nodes 64
   EOL

The ``romio-hints`` file contains pairs of key-value hints to enable and tune collective
buffering of MPI-IO operations. As Summit's Alpine file system uses a 16MB block size,
you should set the collective buffer size to 16GB and tune the number of aggregators
(``cb_nodes``) to the number of compute nodes you are using, i.e., one aggregator per node.

Further details are available at `Summit's documentation page <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#slow-performance-using-parallel-hdf5-resolved-march-12-2019>`__.

.. _building-summit-issues:

Known System Issues
-------------------

.. warning::

   Sep 16th, 2021 (OLCFHELP-3685):
   The **Jupyter** service cannot open HDF5 files without hanging, due to a filesystem mounting problem.

   `Please apply this work-around <https://github.com/openPMD/openPMD-api/pull/1106>`__ in a Jupyter cell before opening any HDF5 files for read:

   .. code-block:: python3

      import os
      os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"

.. warning::

   Aug 27th, 2021 (OLCFHELP-3442):
   Created simulation files and directories are no longer accessible by your team members, even if you create them on ``$PROJWORK``.
   Setting the proper "user mask" (``umask``) does not yet work to fix this.

   Please run those commands *after* running a simulation to fix this.
   You can also append this to the end of your job scripts *after* the ``jsrun`` line:

   .. code-block:: bash

      # cd your-simulation-directory
      find . -type d -exec chmod g+rwx {} \;
      find . -type f -exec chmod g+rw {} \;

.. warning::

   Sep 3rd, 2021 (OLCFHELP-3545):
   The implementation of barriers in IBM's MPI fork is broken and leads to crashes at scale.
   This is seen with runs using 200 nodes and above.

   Our batch script templates above `apply this work-around <https://github.com/ECP-WarpX/WarpX/pull/2283>`__ *before* the call to ``jsrun``, which avoids the broken routines from IBM and trades them for an OpenMPI implementation of collectives:

   .. code-block:: bash

      export OMPI_MCA_coll_ibm_skip_barrier=true

.. warning::

   Sep 3rd, 2021 (OLCFHELP-3319):
   If you are an active developer and compile middleware libraries (e.g., ADIOS2) yourself that use MPI and/or infiniband, be aware of ``libfabric``: IBM forks the open source version of this library and ships a patched version.

   Avoid conflicts with mainline versions of this library in MPI that lead to crashes at runtime by loading alongside the system MPI module:

   .. code-block:: bash

      module load libfabric/1.12.1-sysrdma

   For instance, if you compile large software stacks with Spack, make sure to register ``libfabric`` with that exact version as an external module.

   If you load the documented ADIOS2 module above, this problem does not affect you, since the correct ``libfabric`` version is chosen for this one.

.. warning::

   Related to the above issue, the fabric selection in ADIOS2 was designed for libfabric 1.6.
   With newer versions of libfabric, a workaround is needed to guide the selection of a functional fabric for RDMA support.
   Details are discussed in `ADIOS2 issue #2887 <https://github.com/ornladios/ADIOS2/issues/2887>`__.

   The following environment variables can be set as work-arounds, when working with ADIOS2 SST:

   .. code-block:: bash

      export FABRIC_IFACE=mlx5_0   # ADIOS SST: select interface (1 NIC on Summit)
      export FI_OFI_RXM_USE_SRX=1  # libfabric: use shared receive context from MSG provider

.. warning::

   Oct 12th, 2021 (OLCFHELP-4242):
   There is currently a problem with the pre-installed Jupyter extensions, which can lead to connection splits at long running analysis sessions.

   Work-around this issue by running in a single Jupyter cell, before starting analysis:

   .. code-block:: bash

      !jupyter serverextension enable --py --sys-prefix dask_labextension


.. _post-processing-summit:

Post-Processing
---------------

For post-processing, most users use Python via OLCFs's `Jupyter service <https://jupyter.olcf.ornl.gov>`__ (`Docs <https://docs.olcf.ornl.gov/services_and_applications/jupyter/index.html>`__).

We usually just install our software on-the-fly on Summit.
When starting up a post-processing session, run this in your first cells:

.. code-block:: bash

   # work-around for OLCFHELP-4242
   !jupyter serverextension enable --py --sys-prefix dask_labextension

   # next Jupyter cell: install a faster & better conda package manager
   !conda install -c conda-forge -y mamba

   # next cell: the software you want
   !mamba install -c conda-forge -y openpmd-api openpmd-viewer ipympl ipywidgets fast-histogram yt

   # restart notebook
