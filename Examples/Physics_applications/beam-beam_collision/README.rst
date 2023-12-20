.. _examples-beam-beam_collision:

Beam-beam collision
====================

This example shows how to simulate the collision between two ultra-relativistic particle beams. 
This is representative of what happens at the interaction point of a linear collider. 
We consider a right-propagating electron bunch colliding against a left-propagating positron bunch. 

We turn on the Quantum Synchrotron QED module for photon emission (also known as beamstrahlung in the collider community) and 
the Breit-Wheeler QED module for the generation of electron-positron pairs (also known as coherent pair generation in the collider community). 

To solve for the electromagnetic field we use the nodal version of the electrostatic relativistic solver. 
This solver computes the average velocity of each species, moves to the rest frame of each species where it solves the relativistic Poisson equation, 
and goes back to the lab frame where it adds up all the fields of all species. 
This solver showed good field cancellation properties that are crucial in simulations of relativistic particles.


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

We compare different results:
* (red) simplified WarpX simulation as the example stored in the directory ``/Examples/Physics_applications/beam-beam_collision``;
* (blue) large-scale WarpX simulation (high resolution and ad hoc generated tables ;
* (black) literature results from :cite:t:`ex-Yakimenko2019`.

The small-scale simulation has been performed with a resolution of ``nx = 64, ny = 64, nz = 128`` grid cells, while the large-scale one has a much higher resolution of ``nx = 512, ny = 512, nz = 1024``. Moreover, the large-scale simulation uses dedicated QED lookup tables instead of the builtin tables. To generate the tables within WarpX, the code must be compiled with the flag ``-DWarpX_QED_TABLE_GEN=ON``. For the large-scale simulation we have used the following options:

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



.. figure:: https://user-images.githubusercontent.com/40245517/171573007-f7d733c7-c0de-490c-9ed6-ff4c02154358.png
   :alt: Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.
   :width: 80%

   Beam-beam collision benchmark against :cite:t:`ex-Yakimenko2019`.
