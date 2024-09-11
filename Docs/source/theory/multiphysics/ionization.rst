.. _multiphysics-ionization:

Ionization
==========

Field Ionization
----------------

Under the influence of a sufficiently strong external electric field atoms become ionized.
Particularly the dynamics of interactions between ultra-high intensity laser pulses and matter, e.g., Laser-Plasma Acceleration (LPA) with ionization injection, or Laser-Plasma Interactions with solid density targets (LPI) can depend on field ionization dynamics as well.

WarpX models field ionization based on a description of the Ammosov-Delone-Krainov model:cite:p:`mpion-Ammosov1986` following :cite:t:`mpion-ChenPRSTAB13`.

Implementation Details and Assumptions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    The current implementation makes the following assumptions

    * Energy for ionization processes is not removed from the electromagnetic fields
    * Only one single-level ionization process can occur per macroparticle and time step
    * Ionization happens at the beginning of the PIC loop before the field solve
    * Angular momentum quantum number :math:`l = 0` and magnetic quantum number :math:`m = 0`

Run
---

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: lab frame

        This example can be run **either** as:

        * **Python** script: ``python3 inputs_test_2d_ionization_picmi.py`` or
        * WarpX **executable** using an input file: ``warpx.2d inputs_test_2d_ionization_lab max_step=1600``

        .. tab-set::

            .. tab-item:: Python: Script

                .. literalinclude:: inputs_test_2d_ionization_picmi.py
                :language: python3
                :caption: You can copy this file from ``Examples/Tests/ionization/inputs_test_2d_ionization_picmi.py``.

            .. tab-item:: Executable: Input File

                .. literalinclude:: inputs_test_2d_ionization_lab
                :language: ini
                :caption: You can copy this file from ``Examples/Tests/ionization/inputs_test_2d_ionization_lab``.

   .. tab-item:: boosted frame

        This example can be run as:

        * WarpX **executable** using an input file: ``warpx.2d inputs_test_2d_ionization_boost max_step=420``

        .. literalinclude:: inputs_test_2d_ionization_boost
        :language: ini
        :caption: You can copy this file from ``Examples/Tests/ionization/inputs_test_2d_ionization_boost``.

Analyze
-------

.. note::

.. dropdown:: Script ``analysis.py``

   .. literalinclude:: analysis.py
      :language: python3
      :caption: You can copy this file from ``Examples/Tests/ionization/analysis.py diags/diag1001600/``.

.. bibliography::
    :keyprefix: mpion-