.. _running-cpp:

Running on specific platforms
=============================

.. _running-cpp-cori:

Running on Cori KNL at NERSC
----------------------------

The batch script below can be used to run a WarpX simulation on 2 KNL nodes on
the supercomputer Cori at NERSC. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<job name>`` could be ``laserWakefield``.

.. literalinclude:: ../../../Tools/BatchScripts/batch_cori.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_cori.sh`` and
run
::

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

Running on Cori Haswell at NERSC
--------------------------------

The batch script below can be used to run a WarpX simulation on 1 `Haswell node <https://docs.nersc.gov/systems/cori/>`_ on the supercomputer Cori at NERSC.

.. literalinclude:: ../../../Tools/BatchScripts/batch_cori_haswell.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_cori_haswell.sh`` and
run
::

  sbatch batch_cori_haswell.sh

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on Cori Haswell for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* **4 MPI ranks per Haswell node** (2 MPI ranks per `Intel Xeon E5-2698 v3 <https://ark.intel.com/content/www/us/en/ark/products/81060/intel-xeon-processor-e5-2698-v3-40m-cache-2-30-ghz.html>`_), with ``OMP_NUM_THREADS=16`` (which uses `2x hyperthreading <https://docs.nersc.gov/jobs/affinity/>`_)

.. _running-cpp-summit:

Running on Summit at OLCF
-------------------------

.. _running-cpp-summit-V100-GPUs:

V100 GPUs
^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on
the supercomputer Summit at OLCF. Replace descriptions between chevrons ``<>``
by relevant values, for instance ``<input file>`` could be
``plasma_mirror_inputs``. Note that the only option so far is to run with one
MPI rank per GPU.

.. literalinclude:: ../../../Tools/BatchScripts/batch_summit.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_summit.sh`` and
run
::

  bsub batch_summit.sh

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
:download:`Summit batch script<../../../Tools/BatchScripts/script_profiling_summit.sh>`

.. _running-cpp-summit-Power9-CPUs:

Power9 CPUs
^^^^^^^^^^^

Similar to above, the batch script below can be used to run a WarpX simulation on
1 node on the supercomputer Summit at OLCF, on Power9 CPUs (i.e., the GPUs are
ignored).

.. literalinclude:: ../../../Tools/BatchScripts/batch_summit_power9.sh
   :language: bash

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

.. _running-cpp-lassen:

Running on Lassen at LLNL
-------------------------

.. _running-cpp-lassen-V100-GPUs:

V100 GPUs
^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Lassen at LLNL.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.
Note that the only option so far is to run with one MPI rank per GPU.

.. literalinclude:: ../../../Tools/BatchScripts/batch_lassen.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_lassen.sh`` and run
::

  bsub batch_lassen.sh

to submit the job.

For a 3D simulation with a few (1-4) particles per cell using FDTD Maxwell
solver on V100 GPUs for a well load-balanced problem (in our case laser
wakefield acceleration simulation in a boosted frame in the quasi-linear
regime), the following set of parameters provided good performance:

* ``amr.max_grid_size=256`` and ``amr.blocking_factor=128``.

* **One MPI rank per GPU** (e.g., 4 MPI ranks for the 4 GPUs on each Lassen
  node)

* **Two `128x128x128` grids per GPU**, or **one `128x128x256` grid per GPU**.


Running on Quartz at LLNL
-------------------------

.. _running-cpp-quartz-CPUs:

Intel Xeon E5-2695 v4 CPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^

The batch script below can be used to run a WarpX simulation on 2 nodes on the supercomputer Quartz at LLNL.
Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<input file>`` could be ``plasma_mirror_inputs``.

.. literalinclude:: ../../../Tools/BatchScripts/batch_quartz.sh
   :language: bash

To run a simulation, copy the lines above to a file ``batch_quartz.sh`` and run
::

  sbatch batch_quartz.sh

to submit the job.
