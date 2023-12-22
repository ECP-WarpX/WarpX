.. _examples-laser-ion:

Laser-Ion Acceleration with a Planar Target
===========================================

This example shows how to model laser-ion acceleration with planar targets of solid density :cite:p:`ex-Wilks2001,ex-Bulanov2008,ex-Macchi2013`.
The acceleration mechanism in this scenario depends on target parameters

Although laser-ion acceleration requires full 3D modeling for adequate description of the acceleration dynamics, especially the acceleration field lengths and decay times, this example models a 2D example.
2D modeling can often hint at a qualitative overview of the dynamics, but mostly saves computational costs since the plasma frequency (and Debye length) of the plasma determines the resolution need in laser-solid interaction modeling.

.. note::

   TODO: The Python (PICMI) input file needs to be created.

.. note::

   The resolution of this 2D case is extremely low by default.
   You will need a computing cluster for adequate resolution of the target density, see comments in the input file.

.. warning::

   It is strongly advised to set the parameters ``<species>.zmin / zmax / xmin / ...`` when working with highly dense targets that are limited in one or multiple dimensions.
   The particle creation routine will first create particles everywhere between these limits (`defaulting to box size if unset`), setting particles to invalid only afterwards based on the density profile.
   Not setting these parameters can quickly lead to memory overflows.


Run
---

This example can be run **either** as:

* **Python** script: (*TODO*) or
* WarpX **executable** using an input file: ``warpx.2d inputs_2d``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. note::

         TODO: This input file should be created following the ``inputs_2d`` file.

   .. tab-item:: Executable: Input File

      .. literalinclude:: inputs_2d
         :language: ini
         :caption: You can copy this file from ``Examples/Physics_applications/laser_ion/inputs_2d``.

Analyze
-------

.. note::

   This section is TODO.


Visualize
---------

.. note::

   This section is TODO.
