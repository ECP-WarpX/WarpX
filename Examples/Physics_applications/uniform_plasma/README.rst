.. _examples-uniform-plasma:

Uniform Plasma
==============

This example evolves a uniformly distributed, hot plasma over time.


Run
---

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: 3D

      .. tab-set::

         .. tab-item:: Python: Script

            .. note::

               TODO: This input file should be created following the ``inputs_3d`` file.

            .. literalinclude:: PICMI_inputs_3d.py
               :language: python3
               :caption: You can copy this file from ``usage/examples/lwfa/PICMI_inputs_3d.py``.

         .. tab-item:: Executable: Input File

            This example can be run **either** as WarpX **executable** using an input file: ``warpx.3d inputs_3d``

             .. literalinclude:: inputs_3d
                :language: ini
                :caption: You can copy this file from ``usage/examples/lwfa/inputs_3d``.

   .. tab-item:: 2D

      .. tab-set::

         .. tab-item:: Python: Script

            .. note::

               TODO: This input file should be created following the ``inputs_2d`` file.

            .. literalinclude:: PICMI_inputs_2d.py
               :language: python3
               :caption: You can copy this file from ``usage/examples/lwfa/PICMI_inputs_2d.py``.

         .. tab-item:: Executable: Input File

            This example can be run **either** as WarpX **executable** using an input file: ``warpx.2d inputs_2d``

             .. literalinclude:: inputs_2d
                :language: ini
                :caption: You can copy this file from ``usage/examples/lwfa/inputs_2d``.

Analyze
-------

.. note::

   This section is TODO.


Visualize
---------

.. note::

   This section is TODO.
