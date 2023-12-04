.. _building-karolina:

Karolina (IT4I)
===============

The `Karolina cluster <https://docs.it4i.cz/karolina/introduction/>`_ is located at `IT4I, Technical University of Ostrava <https://www.it4i.cz/en>`__.


Introduction
------------

If you are new to this system, **please see the following resources**:

* `IT4I user guide <https://docs.it4i.cz>`__
* Batch system: `PBS <https://docs.it4i.cz/general/job-submission-and-execution/>`__
* Jupyter service: not provided/documented (yet)
* `Filesystems <https://docs.it4i.cz/karolina/storage/>`__:

  * ``$HOME``: per-user directory, use only for inputs, source and scripts; backed up (25GB default quota)
  * ``/scatch/``: `production directory <https://docs.it4i.cz/karolina/storage/#scratch-file-system>`__; very fast for parallel jobs (20TB default)


.. _building-karolina-preparation:

Preparation
-----------

Use the following commands to download the WarpX source code:

.. code-block:: bash

   git clone https://github.com/ECP-WarpX/WarpX.git $HOME/src/warpx

On Karolina, you can run either on GPU nodes with fast A100 GPUs (recommended) or CPU nodes.

.. tab-set::

   .. tab-item:: A100 GPUs

      We use system software modules, add environment hints and further dependencies via the file ``$HOME/karolina_gpu_warpx.profile``.
      Create it now:

      .. code-block:: bash

         cp $HOME/src/warpx/Tools/machines/karolina-it4i/karolina_gpu_warpx.profile.example $HOME/karolina_gpu_warpx.profile

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: ../../../../Tools/machines/karolina-it4i/karolina_gpu_warpx.profile.example
            :language: bash

      Edit the 2nd line of this script, which sets the ``export proj=""`` variable.
      For example, if you are member of the project ``DD-23-83``, then run ``vi $HOME/karolina_gpu_warpx.profile``.
      Enter the edit mode by typing ``i`` and edit line 2 to read:

      .. code-block:: bash

         export proj="DD-23-83"

      Exit the ``vi`` editor with ``Esc`` and then type ``:wq`` (write & quit).

      .. important::

         Now, and as the first step on future logins to Karolina, activate these environment settings:

         .. code-block:: bash

            source $HOME/karolina_gpu_warpx.profile

      Finally, since Karolina does not yet provide software modules for some of our dependencies, install them once:

      .. code-block:: bash

         bash $HOME/src/warpx/Tools/machines/karolina-it4i/install_gpu_dependencies.sh
         source $HOME/sw/karolina/gpu/venvs/warpx-gpu/bin/activate

      .. dropdown:: Script Details
         :color: light
         :icon: info
         :animate: fade-in-slide-down

         .. literalinclude:: ../../../../Tools/machines/karolina-it4i/install_gpu_dependencies.sh
            :language: bash


   .. tab-item:: CPU Nodes

      CPU usage is documentation is TODO.


.. _building-karolina-compilation:

Compilation
-----------

Use the following :ref:`cmake commands <building-cmake>` to compile the application executable:

.. tab-set::

   .. tab-item:: A100 GPUs

      .. code-block:: bash

         cd $HOME/src/warpx
         rm -rf build_gpu

         cmake -S . -B build_gpu -DWarpX_COMPUTE=CUDA -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
         cmake --build build_gpu -j 12

      The WarpX application executables are now in ``$HOME/src/warpx/build_gpu/bin/``.
      Additionally, the following commands will install WarpX as a Python module:

      .. code-block:: bash

         rm -rf build_gpu_py

         cmake -S . -B build_gpu_py -DWarpX_COMPUTE=CUDA -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
         cmake --build build_gpu_py -j 12 --target pip_install

   .. tab-item:: CPU Nodes

      .. code-block:: bash

         cd $HOME/src/warpx
         rm -rf build_cpu

         cmake -S . -B build_cpu -DWarpX_COMPUTE=OMP -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_DIMS="1;2;RZ;3"
         cmake --build build_cpu -j 12

      The WarpX application executables are now in ``$HOME/src/warpx/build_cpu/bin/``.
      Additionally, the following commands will install WarpX as a Python module:

      .. code-block:: bash

         cd $HOME/src/warpx
         rm -rf build_cpu_py

         cmake -S . -B build_cpu_py -DWarpX_COMPUTE=OMP -DWarpX_PSATD=ON -DWarpX_QED_TABLE_GEN=ON -DWarpX_APP=OFF -DWarpX_PYTHON=ON -DWarpX_DIMS="1;2;RZ;3"
         cmake --build build_cpu_py -j 12 --target pip_install

Now, you can :ref:`submit Karolina compute jobs <running-cpp-karolina>` for WarpX :ref:`Python (PICMI) scripts <usage-picmi>` (:ref:`example scripts <usage-examples>`).
Or, you can use the WarpX executables to submit Karolina jobs (:ref:`example inputs <usage-examples>`).
For executables, you can reference their location in your :ref:`job script <running-cpp-karolina>` or copy them to a location in ``/scatch/``.


.. _building-karolina-update:

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
   git log # press q to exit

And, if needed,

- :ref:`update the karolina_gpu_warpx.profile or karolina_cpu_warpx.profile files <building-karolina-preparation>`,
- log out and into the system, activate the now updated environment profile as usual,
- :ref:`execute the dependency install scripts <building-karolina-preparation>`.

As a last step, clean the build directory ``rm -rf $HOME/src/warpx/build_*`` and rebuild WarpX.


.. _running-cpp-karolina:

Running
-------

.. tab-set::

   .. tab-item:: A100 (40GB) GPUs

      The batch script below can be used to run a WarpX simulation on multiple GPU nodes (change ``#PBS -l select=`` accordingly) on the supercomputer Karolina at IT4I.
      This partition as up to `72 nodes <https://docs.it4i.cz/karolina/hardware-overview/>`__.
      Every node has 8x A100 (40GB) GPUs and 2x AMD EPYC 7763, 64-core, 2.45 GHz processors.

      Replace descriptions between chevrons ``<>`` by relevant values, for instance ``<proj>`` could be ``DD-23-83``.
      Note that we run one MPI rank per GPU.

      .. literalinclude:: ../../../../Tools/machines/karolina-it4i/karolina_gpu.qsub
         :language: bash
         :caption: You can copy this file from ``$HOME/src/warpx/Tools/machines/karolina-it4i/karolina_gpu.qsub``.

      To run a simulation, copy the lines above to a file ``karolina_gpu.qsub`` and run

      .. code-block:: bash

         qsub karolina_gpu.qsub

      to submit the job.


   .. tab-item:: CPU Nodes

      CPU usage is documentation is TODO.


.. _post-processing-karolina:

Post-Processing
---------------

.. note::

   This section was not yet written.
   Usually, we document here how to use a Jupyter service.
