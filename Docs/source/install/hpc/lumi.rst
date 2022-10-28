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

   cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_COMPUTE=HIP
   cmake --build build -j 6

The general :ref:`cmake compile-time options <building-cmake>` apply as usual.


.. _running-cpp-lumi:

Running
-------

.. _running-cpp-lumi-MI250X-GPUs:

MI250X GPUs (2x64 GB)
^^^^^^^^^^^^^^^^^^^^^

.. note::

   TODO: Add batch script template.


.. _post-processing-lumi:

Post-Processing
---------------

.. note::

   TODO: Document any Jupyter or data services.
