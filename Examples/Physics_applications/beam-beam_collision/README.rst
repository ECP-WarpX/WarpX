.. _examples-beam-beam_collision:

Beam-beam collision
====================

This example shows how to simulate the collision between two ultra-relativistic particle beams.
This is representative of what happens at the interaction point of a linear collider.
We consider a right-propagating electron bunch colliding against a left-propagating positron bunch.

We turn on the Quantum Synchrotron QED module for photon emission (also known as beamstrahlung in the collider community) and
the Breit-Wheeler QED module for the generation of electron-positron pairs (also known as coherent pair generation in the collider community).

This solver computes the average velocity of each species, and solves the corresponding relativistic Poisson equation (see the WarpX documentation for ``warpx.do_electrostatic = relativistic`` for more details).
This solver accurately reproduces the subtle cancellation that occurs for some components of ``E + v x B``, which is crucial in simulations of relativistic particles.

This example is based on the following paper :cite:t:`ex-Yakimenko2019`.


Run
---

The PICMI input file is not available for this example yet.

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: inputs
   :language: ini
   :caption: You can copy this file from ``Examples/Physics_applications/beam-beam_collision/inputs``.


Visualize
---------

The figure below shows the number of photons emitted per beam particle (left) and the number of secondary pairs generated per beam particle (right).
We compare different results for the reduced diagnostics with some literature:

* (red) simplified WarpX simulation as the example stored in the directory ``/Examples/Physics_applications/beam-beam_collision``;
* (blue) large-scale WarpX simulation (high resolution and ad hoc generated tables ;
* (black) literature results from :cite:t:`ex-Yakimenko2019`.

The small-scale simulation has been performed with a resolution of ``nx = 64, ny = 64, nz = 128`` grid cells, while the large-scale one with a much higher resolution of ``nx = 512, ny = 512, nz = 1024``.
Moreover, the large-scale simulation uses dedicated QED lookup tables instead of the builtin ones.
To generate the tables within WarpX, the code must be compiled with the flag ``-DWarpX_QED_TABLE_GEN=ON``.
For the large-scale simulation we have used the following options:

.. code-block:: ini

   qed_qs.lookup_table_mode = generate
   qed_bw.lookup_table_mode = generate
   qed_qs.tab_dndt_chi_min=1e-3
   qed_qs.tab_dndt_chi_max=2e3
   qed_qs.tab_dndt_how_many=512
   qed_qs.tab_em_chi_min=1e-3
   qed_qs.tab_em_chi_max=2e3
   qed_qs.tab_em_chi_how_many=512
   qed_qs.tab_em_frac_how_many=512
   qed_qs.tab_em_frac_min=1e-12
   qed_qs.save_table_in=my_qs_table.txt
   qed_bw.tab_dndt_chi_min=1e-2
   qed_bw.tab_dndt_chi_max=2e3
   qed_bw.tab_dndt_how_many=512
   qed_bw.tab_pair_chi_min=1e-2
   qed_bw.tab_pair_chi_max=2e3
   qed_bw.tab_pair_chi_how_many=512
   qed_bw.tab_pair_frac_how_many=512
   qed_bw.save_table_in=my_bw_table.txt

.. figure:: https://user-images.githubusercontent.com/17280419/291749626-aa61fff2-e6d2-45a3-80ee-84b2851ea0bf.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDMwMzQzNTEsIm5iZiI6MTcwMzAzNDA1MSwicGF0aCI6Ii8xNzI4MDQxOS8yOTE3NDk2MjYtYWE2MWZmZjItZTZkMi00NWEzLTgwZWUtODRiMjg1MWVhMGJmLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFJV05KWUFYNENTVkVINTNBJTJGMjAyMzEyMjAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjMxMjIwVDAxMDA1MVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWFiYzY2MGQyYzIyZGIzYzUxOWI3MzNjZTk5ZDM1YzgyNmY4ZDYxOGRlZjAyZTIwNTAyMTc3NTgwN2Q0YjEwNGMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0JmFjdG9yX2lkPTAma2V5X2lkPTAmcmVwb19pZD0wIn0.I96LQpjqmFXirPDVnBlFQIkCuenR6IuOSY0OIIQvtCo
   :alt: Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.
   :width: 100%

   Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.


Below are two visualizations scripts that provide examples to graph the field and reduced diagnostics.
They are available in the ``Examples/Physics_applications/beam-beam_collision/`` folder and can be run as simply as ``python3 plot_fields.py`` and ``python3 plot_reduced.py``.

.. tab-set::

   .. tab-item:: Field Diagnostics

      This script visualizes the evolution of the fields (:math:`|E|, |B|, \rho`) during the collision between the two ultra-relativistic lepton beams.
      The magnitude of E and B and the charge densities of the primary beams and of the secondary pairs are sliced along either one of the two transverse coordinates (:math:`x` and `:math:`y`).

      .. literalinclude:: plot_fields.py
         :language: python3
         :caption: You can copy this file from ``Examples/Physics_applications/beam-beam_collision/plot_fields.py``.

      .. figure:: https://gist.github.com/user-attachments/assets/04c9c0ec-b580-446f-a11a-437c1b244a41
         :alt: Slice across :math:`x` of different fields (:math:`|E|, |B|, \rho`) at timestep 45, in the middle of the collision.
         :width: 100%

         Slice across :math:`x` of different fields (:math:`|E|, |B|, \rho`) at timestep 45, in the middle of the collision.


   .. tab-item:: Reduced Diagnostics

      A similar script to the one below was used to produce the image showing the benchmark against :cite:t:`ex-Yakimenko2019`.

      .. literalinclude:: plot_reduced.py
         :language: python3
         :caption: You can copy this file from ``Examples/Physics_applications/beam-beam_collision/plot_reduced.py``.

      .. figure:: https://gist.github.com/user-attachments/assets/c280490a-f1f2-4329-ad3c-46817d245dc1
         :alt: Photon and pair production rates in time throughout the collision.
         :width: 100%

         Photon and pair production rates in time throughout the collision.
