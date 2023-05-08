.. _building-lumi:

LUMI (CSC)
==========

The `LUMI cluster <https://www.lumi-supercomputer.eu>`_ is located at CSC.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `Lumi user guide <https://docs.lumi-supercomputer.eu>`_
* Batch system: `Slurm <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/>`_
* `Jupyter service ? <TODO>`__
* Production directories:

  * `LUMI-P <https://docs.lumi-supercomputer.eu/hardware/storage/lumip/>`__: 4 independent [Lustre](https://docs.lumi-supercomputer.eu/hardware/storage/lumip/#lustre) file systems
  * `LUMI-F <https://docs.lumi-supercomputer.eu/hardware/storage/lumif/>`__: a fast Lustre file system
  * `LUMI-O <https://docs.lumi-supercomputer.eu/hardware/storage/lumio/>`__: object storage


Installation
------------

Use the following commands to download the WarpX source code and switch to the correct branch:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use the following modules and environments on the system (``$HOME/lumi_warpx.profile``).

.. literalinclude:: ../../../../Tools/machines/lumi-csc/lumi_warpx.profile.example
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lumi-csc/lumi_warpx.profile.example``.


We recommend to store the above lines in a file, such as ``$HOME/lumi_warpx.profile``, and load it into your shell after a login:

.. code-block:: bash

   source $HOME/lumi_warpx.profile

Then, ``cd`` into the directory ``$HOME/src/warpx`` and use the following commands to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=HIP -DWarpX_PSATD=ON
   cmake --build build -j 6

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.

**That's it!**
A 3D WarpX executable is now in ``build/bin/`` and :ref:`can be run <running-cpp-lumi>` with a :ref:`3D example inputs file <usage-examples>`.
Most people execute the binary directly or copy it out to a location in LUMI-P.


.. _running-cpp-lumi:

Running
-------

.. _running-cpp-lumi-MI250X-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

In non-interactive runs:

.. literalinclude:: ../../../../Tools/machines/lumi-csc/submit.sh
   :language: bash
   :caption: You can copy this file from ``Tools/machines/lumi-csc/submit.sh``.

.. _post-processing-lumi:

Post-Processing
---------------

.. note::

   TODO: Document any Jupyter or data services.

Known System Issues
-------------------

.. warning::

   December 12th, 2022:
   There is a caching bug in libFabric that causes WarpX simulations to occasionally hang on LUMI on more than 1 node.

   As a work-around, please export the following environment variable in your job scripts until the issue is fixed:

   .. code-block:: bash

      #export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
      # or, less invasive:
      export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

.. warning::

   January, 2023:
   We discovered a regression in AMD ROCm, leading to 2x slower current deposition (and other slowdowns) in ROCm 5.3 and 5.4.
   Reported to AMD and fixed for the 5.5 release of ROCm.

   Upgrade ROCm or stay with the ROCm 5.2 module to avoid.
