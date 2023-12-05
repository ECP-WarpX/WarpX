.. _examples-pwfa:

Beam-Driven Wakefield Acceleration of Electrons
===============================================

This example shows how to model a beam-driven plasma-wakefield accelerator (PWFA) :cite:p:`TajimaDawson1982,Esarey1996`.

PWFA is best performed in 3D and quasi-cylindrical (RZ) geometry, which ensures that the plasma wavelength of the wakefield is modelled with the right scale lengths.
RZ modeling enables efficient modeling if effects of asymmetry shall be ignored (e.g., asymmetric beams and transverse profiles, hosing of the injected beam, etc.).

Additionally, to speed up computation, this example uses the :ref:`boosted frame method <theory-boostedframe>` to effectively model long acceleration lengths.

Alternatively, an other common approximation for PWFAs is quasi-static modeling, e.g., if effects such as self-injection can be ignored.
In the Beam, Plasma & Accelerator Simulation Toolkit (BLAST), `HiPACE++ <https://hipace.readthedocs.io>`__ provides such methods.

.. note::

   TODO: The Python (PICMI) input file should use the boosted frame method, like the ``inputs_3d_boost`` file.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 PICMI_inputs_plasma_acceleration.py`` or
* WarpX **executable** using an input file: ``warpx.3d inputs_3d_boost``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. note::

         TODO: This input file should use the boosted frame method, like the ``inputs_3d_boost`` file.

      .. literalinclude:: PICMI_inputs_3d.py
         :language: python3
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/PICMI_inputs_plasma_acceleration.py``.

   .. tab-item:: Executable: Input File

      .. literalinclude:: inputs_3d
         :language: ini
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/inputs_3d_boost``.

Analyze
-------

.. note::

   This section is TODO.

We run the following script to analyze correctness:

.. dropdown:: Script ``analysis_3d.py``

   .. literalinclude:: analysis_3d.py
      :language: python3
      :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/analysis_3d.py``.


Visualize
---------

.. note::

   This section is TODO.
