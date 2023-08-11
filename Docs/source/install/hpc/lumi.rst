.. _building-lumi:

LUMI (CSC)
==========

The `LUMI cluster <https://www.lumi-supercomputer.eu>`_ is located at CSC (Finland).
Each node contains 4 AMD MI250X GPUs, each with 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node.
You can think of the 8 GCDs as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

Introduction
------------

If you are new to this system, **please see the following resources**:

* `Lumi user guide <https://docs.lumi-supercomputer.eu>`_
* Batch system: `Slurm <https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/slurm-quickstart/>`_
* `Data analytics and visualization <https://docs.lumi-supercomputer.eu/hardware/lumid/>`__
* `Production directories <https://docs.lumi-supercomputer.eu/storage/>`__:

  * ``$HOME``: single user, intended to store user configuration files and personal data (20GB default quota)
  * ``/project/$proj``: shared with all members of a project, purged at the end of a project (50 GB default quota)
  * ``/scratch/$proj``: temporary storage, main storage to be used for disk I/O needs when running simulations on LUMI, purged every 90 days (50TB default quota)


.. _building-lumi-preparation:

Preparation
-----------

Use the following commands to download the WarpX source code:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

We use system software modules, add environment hints and further dependencies via the file ``$HOME/lumi_warpx.profile``.
Create it now:

.. code-block:: bash

   cp $HOME/src/warpx/Tools/machines/lumi-csc/lumi_warpx.profile.example $HOME/lumi_warpx.profile

.. dropdown:: Script Details
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../Tools/machines/lumi-csc/lumi_warpx.profile.example
      :language: bash

Edit the 2nd line of this script, which sets the ``export proj=""`` variable using a text editor
such as ``nano``, ``emacs``, or ``vim`` (all available by default on
LUMI login nodes).

.. important::

   Now, and as the first step on future logins to LUMI, activate these environment settings:

   .. code-block:: bash

      source $HOME/lumi_warpx.profile

Finally, since LUMI does not yet provide software modules for some of our dependencies, install them once:

.. code-block:: bash

   bash $HOME/src/warpx/Tools/machines/lumi-csc/install_dependencies.sh
   source $HOME/sw/lumi/gpu/venvs/warpx-lumi/bin/activate

.. dropdown:: Script Details
   :color: light
   :icon: info
   :animate: fade-in-slide-down

   .. literalinclude:: ../../../../Tools/machines/lumi-csc/install_dependencies.sh
      :language: bash


.. _building-lumi-compilation:

Compilation
-----------

Use the following :ref:`cmake commands <building-cmake>` to compile:

.. code-block:: bash

   cd $HOME/src/warpx
   rm -rf build_lumi

   cmake -S . -B build_lumi -DWarpX_COMPUTE=HIP -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_LIB=ON -DWarpX_DIMS="1;2;RZ;3"
   cmake --build build_lumi -j 16
   cmake --build build_lumi -j 16 --target pip_install

**That's it!**
WarpX executables are now in ``build_lumi/bin/`` and :ref:`can be run <running-cpp-lumi>` with matching :ref:`example inputs files <usage-examples>`.
Most people execute the binary directly or copy it out to a location in ``/scratch/<project>``.


.. _building-lumi-update:

Update WarpX & Dependencies
---------------------------

If you already installed WarpX in the past and want to update it, start by getting the latest source code:

.. code-block:: bash

   cd $HOME/src/warpx

   # read the output of this command - does it look ok?
   git status

   # get the latest WarpX source code
   git fetch
   git pull

   # read the output of these commands - do they look ok?
   git status
   git log     # press q to exit

And, if needed,

- :ref:`update the lumi_warpx.profile file <building-lumi-preparation>`,
- log out and into the system, activate the now updated environment profile as usual,
- :ref:`execute the dependency install scripts <building-lumi-preparation>`.

As a last step, clean the build directory ``rm -rf $HOME/src/warpx/build_lumi`` and rebuild WarpX.


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


.. _known-lumi-issues:

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

   June, 2023:
   Although a fix was planned for ROCm 5.5, we still see the same issue in this release and continue to exchange with AMD and HPE on the issue.

   Stay with the ROCm 5.2 module to avoid a 2x slowdown.

.. warning::

   May 2023:
   rocFFT in ROCm 5.1-5.3 tries to `write to a cache <https://rocfft.readthedocs.io/en/latest/#runtime-compilation>`__ in the home area by default.
   This does not scale, disable it via:

   .. code-block:: bash

      export ROCFFT_RTC_CACHE_PATH=/dev/null
